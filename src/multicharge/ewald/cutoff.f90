! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

!> @file multicharge/ewald/cutoff.f90
!> Ewald splitting parameter, cutoff search and lattice translation generation

!> Ewald splitting parameter and tolerance-based direct/reciprocal translations
!> for 1D, 2D, and 3D periodic boundary conditions
module multicharge_ewald_cutoff
   use mctc_cutoff, only : get_inverse_metric, get_lattice_points, &
      & get_periodic_axes, reduce_basis
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : crossprod, matdet_3x3
   use multicharge_ewald_pbc1d, only : get_rec_term_1d
   use multicharge_ewald_pbc2d, only : get_rec_term_2d
   use multicharge_ewald_pbc3d, only : get_rec_term_3d
   implicit none
   private

   public :: get_alpha, get_dir_cutoff, get_rec_cutoff, get_ewald_trans
   public :: compute_alpha, get_ewald_geometry, build_ewald_trans

   real(wp), parameter :: twopi = 2.0_wp * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Dimension-dependent scaling of alpha = s_d * sqrt(pi) / measure^(1/d):
   !> - 1D (0.50): shifts work toward linear-scaling direct space (erf) and away
   !>   from per-pair incomplete Bessel quadrature (fast_k012_1d), ~2x faster than 1.0.
   !> - 2D (1.00): balanced split where direct (erf) and per-pair reciprocal
   !>   (erfc_exp_2d) kernels have comparable per-vector cost and quadratic scaling.
   !> - 3D (1.25): shifts work toward reciprocal space, whose Gaussian weights are
   !>   precomputed once per cell so pairs only evaluate trigonometric factors.
   real(wp), parameter :: alpha_scale(3) = [0.5_wp, 1.0_wp, 1.25_wp]

   !> Tolerance for a single term of the real- or reciprocal-space potential sum.
   !> Potential-based cutoffs at machine epsilon already converge energies,
   !> gradients, strain derivatives, and Hessians to floating-point / linear-solve
   !> precision across 1D, 2D, and 3D without heuristic derivative prefactors.
   real(wp), parameter :: term_tol = epsilon(0.0_wp)

   abstract interface
      !> Max. value of a real or reciprocal space term for a given vector length
      pure function get_term_gen(dist, alpha, measure) result(val)
         import :: wp

         !> Length of the real or reciprocal space vector
         real(wp), intent(in) :: dist

         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha

         !> Generalized cell measure (length, area, or volume)
         real(wp), intent(in) :: measure

         !> Value of the term
         real(wp) :: val

      end function get_term_gen
   end interface


contains


!> Ewald parameter balancing the cost of the real-space and reciprocal sums
subroutine get_alpha(lattice, alpha, error, periodic)

   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Parameter of the Ewald summation
   real(wp), intent(out) :: alpha

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional periodicity mask along each lattice axis (defaults to 3D)
   logical, intent(in), optional :: periodic(:)

   logical :: mask(3)
   integer :: ndim
   real(wp) :: measure, axis_vec(3), red_lat(3, 3), rec_lat(3, 3)

   alpha = 0.0_wp
   call get_ewald_geometry(lattice, mask, ndim, measure, axis_vec, red_lat, &
      & rec_lat, error, periodic)
   if (allocated(error)) return

   alpha = compute_alpha(ndim, measure)

end subroutine get_alpha


!> Closed-form Ewald splitting parameter from dimensionality and cell measure
pure function compute_alpha(ndim, measure) result(alpha)

   !> Number of periodic dimensions (1, 2, or 3)
   integer, intent(in) :: ndim

   !> Generalized cell measure: length (1D), area (2D), or volume (3D)
   real(wp), intent(in) :: measure

   !> Ewald splitting parameter
   real(wp) :: alpha

   alpha = alpha_scale(ndim) * sqrtpi / measure**(1.0_wp / real(ndim, wp))

end function compute_alpha


!> Extract active periodic subspace geometry, reduced basis, and reciprocal basis
pure subroutine get_ewald_geometry(lattice, mask, ndim, measure, axis_vec, &
   & red_lat, rec_lat, error, periodic)

   !> Lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Normalized 3-element periodicity mask
   logical, intent(out) :: mask(3)

   !> Number of active periodic dimensions (1, 2, or 3)
   integer, intent(out) :: ndim

   !> Generalized cell measure: length (1D), area (2D), or volume (3D)
   real(wp), intent(out) :: measure

   !> Unit wire direction (1D) or unit slab normal (2D); zero in 3D
   real(wp), intent(out) :: axis_vec(3)

   !> Pairwise size-reduced active lattice vectors (zero in nonperiodic columns)
   real(wp), intent(out) :: red_lat(3, 3)

   !> Reciprocal lattice basis in the active periodic subspace (zero in inactive columns)
   real(wp), intent(out) :: rec_lat(3, 3)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional periodicity mask (1-element or 3-element; defaults to 3D)
   logical, intent(in), optional :: periodic(:)

   integer :: p_axes(3), k
   real(wp) :: lat_chk(3, 3), ginv(3, 3), n_vec(3)

   mask(:) = .false.
   ndim = 0
   measure = 0.0_wp
   axis_vec(:) = 0.0_wp
   red_lat(:, :) = 0.0_wp
   rec_lat(:, :) = 0.0_wp

   if (size(lattice, 1) < 3 .or. size(lattice, 2) < 3) then
      call fatal_error(error, "Lattice vectors must be a 3x3 matrix")
      return
   end if

   if (present(periodic)) then
      call get_periodic_axes(periodic, ndim, p_axes)
   else
      ndim = 3
      p_axes(:) = [1, 2, 3]
   end if

   if (ndim <= 0) then
      call fatal_error(error, &
         & "Periodic Ewald setup requires at least one periodic dimension")
      return
   end if

   mask(p_axes(:ndim)) = .true.
   do k = 1, ndim
      red_lat(:, p_axes(k)) = lattice(:3, p_axes(k))
   end do

   call reduce_basis(red_lat)
   call get_inverse_metric(mask, red_lat, lat_chk, ginv, error)
   if (allocated(error)) return

   select case (ndim)
   case (1)
      measure = norm2(lattice(:3, p_axes(1)))
      if (measure < eps) then
         call fatal_error(error, "Lattice vectors are singular or degenerate")
         return
      end if
      axis_vec(:) = lattice(:3, p_axes(1)) / measure
   case (2)
      n_vec(:) = crossprod(lattice(:3, p_axes(1)), lattice(:3, p_axes(2)))
      measure = norm2(n_vec)
      if (measure < eps) then
         call fatal_error(error, "Lattice vectors are singular or degenerate")
         return
      end if
      axis_vec(:) = n_vec / measure
   case (3)
      measure = abs(matdet_3x3(lattice(:3, :3)))
      if (measure < eps) then
         call fatal_error(error, "Lattice vectors are singular or degenerate")
         return
      end if
   case default
      call fatal_error(error, "Unsupported periodic dimensionality")
      return
   end select

   rec_lat(:, :) = twopi * matmul(red_lat, ginv)

end subroutine get_ewald_geometry


!> Real-space cutoff
pure function get_dir_cutoff(alpha, conv) result(cutoff)

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Tolerance for a single real-space term
   real(wp), intent(in) :: conv

   !> Real-space cutoff
   real(wp) :: cutoff

   cutoff = search_cutoff(get_dir_term, alpha, 0.0_wp, conv)

end function get_dir_cutoff


!> Reciprocal-space cutoff for a full-space single-vector potential contribution
pure function get_rec_cutoff(alpha, measure, conv, ndim) result(cutoff)

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Generalized cell measure: length (1D), area (2D), or volume (3D)
   real(wp), intent(in) :: measure

   !> Tolerance for a single full-space reciprocal-space term
   real(wp), intent(in) :: conv

   !> Optional periodic dimensionality (1, 2, or 3; defaults to 3)
   integer, intent(in), optional :: ndim

   !> Reciprocal-space cutoff
   real(wp) :: cutoff

   integer :: dim_eff

   dim_eff = 3
   if (present(ndim)) dim_eff = ndim

   select case (dim_eff)
   case (1)
      cutoff = search_cutoff(get_rec_term_1d, alpha, measure, conv)
   case (2)
      cutoff = search_cutoff(get_rec_term_2d, alpha, measure, conv)
   case (3)
      cutoff = search_cutoff(get_rec_term_3d, alpha, measure, conv)
   case default
      cutoff = 0.0_wp
   end select

end function get_rec_cutoff


!> Direct and half-space reciprocal translations within the single-term cutoffs
pure subroutine get_ewald_trans(lattice, alpha, gam, rmax, dtrans, rtrans, &
   & error, periodic)

   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Smallest Gaussian charge exponent of all pairs
   real(wp), intent(in) :: gam

   !> Largest projected length of the pair vectors
   real(wp), intent(in) :: rmax

   !> Real-space lattice translations
   real(wp), allocatable, intent(out) :: dtrans(:, :)

   !> Reciprocal lattice translations in a half space
   real(wp), allocatable, intent(out) :: rtrans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional periodicity mask along each lattice axis (defaults to 3D)
   logical, intent(in), optional :: periodic(:)

   logical :: mask(3)
   integer :: ndim
   real(wp) :: measure, axis_vec(3), red_lat(3, 3), rec_lat(3, 3)

   call get_ewald_geometry(lattice, mask, ndim, measure, axis_vec, red_lat, &
      & rec_lat, error, periodic)
   if (allocated(error)) return

   call build_ewald_trans(mask, ndim, measure, red_lat, rec_lat, alpha, gam, &
      & rmax, dtrans, rtrans, error)

end subroutine get_ewald_trans


!> Build direct and half-space reciprocal translations from prepared lattice geometry
pure subroutine build_ewald_trans(mask, ndim, measure, red_lat, rec_lat, alpha, &
   & gam, rmax, dtrans, rtrans, error, tol_scale)

   !> Normalized 3-element periodicity mask
   logical, intent(in) :: mask(3)

   !> Number of periodic dimensions (1, 2, or 3)
   integer, intent(in) :: ndim

   !> Generalized cell measure: length (1D), area (2D), or volume (3D)
   real(wp), intent(in) :: measure

   !> Reduced active real-space lattice vectors (zero in inactive columns)
   real(wp), intent(in) :: red_lat(3, 3)

   !> Reciprocal lattice basis in the active subspace (zero in inactive columns)
   real(wp), intent(in) :: rec_lat(3, 3)

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Smallest Gaussian charge exponent of all pairs
   real(wp), intent(in) :: gam

   !> Largest projected length of the pair vectors
   real(wp), intent(in) :: rmax

   !> Real-space lattice translations
   real(wp), allocatable, intent(out) :: dtrans(:, :)

   !> Reciprocal lattice translations in a half space
   real(wp), allocatable, intent(out) :: rtrans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional multiplicative scale on term_tol (defaults to 1.0)
   real(wp), intent(in), optional :: tol_scale

   integer :: rep(3), nrec, itr, ix, iy, iz, jy, jz
   real(wp) :: tol, rcut, gcut, vec(3), extent(3)
   real(wp), allocatable :: tmp(:, :)

   if (.not. (alpha > 0.0_wp .and. gam > 0.0_wp .and. rmax >= 0.0_wp)) then
      call fatal_error(error, "Invalid Ewald cutoff parameters")
      return
   end if

   tol = term_tol
   if (present(tol_scale)) tol = tol_scale * term_tol

   ! The real-space kernel decays with the smaller of alpha and gam
   rcut = get_dir_cutoff(min(alpha, gam), tol)

   ! Each half-space vector carries its +G/-G pair (2x the full-space single term)
   gcut = get_rec_cutoff(alpha, measure, 0.5_wp * tol, ndim)

   ! Pair vectors shift the real-space sum within the periodic subspace by at most rmax
   call get_lattice_points(mask, red_lat, rcut + rmax, tmp, error)
   if (allocated(error)) return

   dtrans = tmp(:, pack([(itr, itr = 1, size(tmp, 2))], &
      & norm2(tmp, dim=1) <= rcut + rmax))
   deallocate(tmp)

   ! Reciprocal vectors within the cutoff sphere, one of each +G/-G pair.
   ! Inactive columns of red_lat are zero, giving extent = 0 and rep = 0.
   extent(:) = gcut * norm2(red_lat, dim=1) / twopi
   if (too_many_points(extent)) then
      call fatal_error(error, "Too many reciprocal-space vectors for the Ewald sum, "//&
         & "lattice is likely degenerate")
      return
   end if
   rep(:) = floor(extent)
   allocate(tmp(3, (product(2 * rep + 1) - 1) / 2))
   nrec = 0
   do ix = 0, rep(1)
      do iy = 0, rep(2)
         do iz = 0, rep(3)
            if (ix == 0 .and. iy == 0 .and. iz == 0) cycle
            do jy = 1, merge(-1, 1, ix > 0 .and. iy > 0), -2
               do jz = 1, merge(-1, 1, (ix > 0 .or. iy > 0) .and. iz > 0), -2
                  vec(:) = rec_lat(:, 1) * real(ix, wp) &
                     & + rec_lat(:, 2) * real(iy * jy, wp) &
                     & + rec_lat(:, 3) * real(iz * jz, wp)
                  if (norm2(vec) > gcut) cycle
                  nrec = nrec + 1
                  tmp(:, nrec) = vec
               end do
            end do
         end do
      end do
   end do
   rtrans = tmp(:, :nrec)

end subroutine build_ewald_trans


!> Check whether a lattice box of the given extents exceeds the integer range
pure function too_many_points(extent) result(too_many)

   !> Number of repetitions along each lattice direction
   real(wp), intent(in) :: extent(3)

   !> Box with (2*extent + 1) points per direction cannot be indexed
   logical :: too_many

   ! Truncate in real arithmetic, floor would overflow for the extents rejected here
   too_many = product(2.0_wp * aint(extent) + 1.0_wp) > real(huge(1), wp)

end function too_many_points


!> Find the cutoff for an Ewald term
pure function search_cutoff(get_term, alpha, measure, conv) result(cutoff)

   !> Procedure to evaluate the Ewald term
   procedure(get_term_gen) :: get_term

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Generalized cell measure (length, area, or volume)
   real(wp), intent(in) :: measure

   !> Tolerance for a single term
   real(wp), intent(in) :: conv

   !> Cutoff for the vector length
   real(wp) :: cutoff

   integer, parameter :: niter = 30
   real(wp) :: xl, xr, x
   integer :: iter

   xr = eps
   do while (get_term(xr, alpha, measure) > conv .and. xr <= 0.5_wp * huge(1.0_wp))
      xr = 2.0_wp * xr
   end do
   xl = 0.5_wp * xr

   do iter = 1, niter
      if (xr - xl <= eps * xr) exit
      x = 0.5_wp * (xl + xr)
      if (get_term(x, alpha, measure) > conv) then
         xl = x
      else
         xr = x
      end if
   end do
   cutoff = xr

end function search_cutoff


!> Returns the max. value of a term in the real space
pure function get_dir_term(rr, alpha, measure) result(dval)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Generalized cell measure (unused)
   real(wp), intent(in) :: measure

   !> Real space term
   real(wp) :: dval

   dval = erfc(alpha * rr) / rr

end function get_dir_term


end module multicharge_ewald_cutoff
