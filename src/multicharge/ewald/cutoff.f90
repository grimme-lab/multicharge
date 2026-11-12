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
!> Ewald splitting parameter and lattice cutoff generation for 1D, 2D, and 3D PBC

!> Ewald splitting parameter search and converged direct/reciprocal lattice cutoffs
module multicharge_ewald_cutoff
   use mctc_cutoff, only : get_inverse_metric, get_lattice_points, get_periodic_axes
   use mctc_env, only : error_type, wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : crossprod, matdet_3x3
   implicit none
   private

   public :: get_alpha, get_ewald_dir_trans, get_ewald_rec_trans

   real(wp), parameter :: twopi = 2.0_wp * pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   abstract interface
      !> Returns the maximum value of a reciprocal-space Ewald term for a given vector length
      pure function get_rec_term_gen(gg, alpha, vol) result(g_term)
         import :: wp

         !> Length of the reciprocal space vector
         real(wp), intent(in) :: gg

         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha

         !> Volume (or area/length) of the real space unit cell
         real(wp), intent(in) :: vol

         !> Reciprocal term
         real(wp) :: g_term

      end function get_rec_term_gen
   end interface


contains


!> Determine the optimal Ewald splitting parameter alpha
subroutine get_alpha(lattice, alpha, error, periodic, vol, rec_lat)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Ewald splitting parameter
   real(wp), intent(out) :: alpha

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional periodicity mask (defaults to 3D periodic)
   logical, intent(in), optional :: periodic(:)

   !> Optional unit cell volume (or area/length)
   real(wp), intent(out), optional :: vol

   !> Optional reciprocal lattice vectors
   real(wp), intent(out), optional :: rec_lat(3, 3)

   integer :: ndim, p_axes(3)
   real(wp) :: volume, red_lat(3, 3), rlat(3, 3)

   alpha = 0.0_wp

   if (present(periodic)) then
      call prepare_periodic_lattice(periodic, lattice, ndim, p_axes, &
         & red_lat, rlat, volume, error)
   else
      call prepare_periodic_lattice([.true., .true., .true.], lattice, &
         & ndim, p_axes, red_lat, rlat, volume, error)
   end if
   if (allocated(error)) return

   ! Balance the decay rates of the direct- and reciprocal-space Ewald sums
   select case (ndim)
   case (0)
      alpha = 0.25_wp
   case (1)
      call search_alpha(red_lat(:, p_axes(:1)), rlat(:, p_axes(:1)), &
         & volume, eps, alpha, get_rec_term_1d)
   case (2)
      call search_alpha(red_lat(:, p_axes(:2)), rlat(:, p_axes(:2)), &
         & volume, eps, alpha, get_rec_term_2d)
   case (3)
      call search_alpha(red_lat, rlat, volume, eps, alpha, get_rec_term_3d)
      ! Prevent small alpha and excessive direct-space boxes in anisotropic 3D cells
      if (volume > eps) then
         alpha = max(alpha, sqrt(pi) / volume**(1.0_wp / 3.0_wp))
      end if
   case default
      alpha = 0.25_wp
   end select

   if (present(vol)) vol = volume
   if (present(rec_lat)) rec_lat = rlat

end subroutine get_alpha


!> Generate converged direct-space lattice translations for 1D, 2D, or 3D Ewald summation
subroutine get_ewald_dir_trans(periodic, lattice, alpha, gam_min, trans, error)

   !> Periodicity mask along each lattice axis
   logical, intent(in) :: periodic(:)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Minimum Gaussian exponent min_ij(gamma_ij)
   real(wp), intent(in) :: gam_min

   !> Direct-space lattice translations (3, ntrans)
   real(wp), allocatable, intent(out) :: trans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: tol_wsc = 0.3_wp

   integer :: ndim, p_axes(3), i
   real(wp) :: red_lat(3, 3), rec_lat(3, 3), vol
   real(wp) :: kappa, kappa_max, tol_val, u_cut, rcut, len_i, min_len, r_disp, r_dir

   call prepare_periodic_lattice(periodic, lattice, ndim, p_axes, red_lat, rec_lat, &
      & vol, error)
   if (allocated(error)) return

   if (ndim == 0) then
      allocate(trans(3, 1), source=0.0_wp)
      return
   end if

   ! Determine direct-space cutoff radius from slowest Gaussian/Ewald decay rate
   tol_val = epsilon(1.0_wp)
   kappa = min(alpha, gam_min)
   kappa_max = max(alpha, gam_min)
   u_cut = sqrt(max(42.25_wp, -log(tol_val) &
      & + log(max(1.0_wp, 4.0_wp * kappa_max**3)) + 4.0_wp))
   rcut = u_cut / kappa
   do while (get_dir_term(rcut, kappa) > tol_val)
      rcut = 1.1_wp * rcut
   end do

   ! Bound parallel Wigner-Seitz pair displacements (including competing images)
   r_disp = 0.0_wp
   min_len = huge(1.0_wp)
   do i = 1, ndim
      len_i = norm2(red_lat(:, p_axes(i)))
      r_disp = r_disp + 0.5_wp * len_i
      min_len = min(min_len, len_i)
   end do
   r_disp = sqrt(max(r_disp, min_len)**2 + tol_wsc)
   r_dir = rcut + r_disp

   call get_lattice_points(periodic, red_lat, r_dir, trans, error)
   if (allocated(error)) return

end subroutine get_ewald_dir_trans


!> Generate converged nonzero reciprocal-space lattice vectors for 1D, 2D, or 3D Ewald
subroutine get_ewald_rec_trans(periodic, lattice, alpha, trans, error, vol)

   !> Periodicity mask along each lattice axis
   logical, intent(in) :: periodic(:)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Nonzero reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: trans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional generalized cell measure (length in 1D, area in 2D, volume in 3D)
   real(wp), intent(out), optional :: vol

   integer :: ndim, p_axes(3), rep(3), i, ax
   real(wp) :: red_lat(3, 3), rec_lat(3, 3), cell_vol
   real(wp) :: tol_val, u_cut, gcut

   call prepare_periodic_lattice(periodic, lattice, ndim, p_axes, red_lat, rec_lat, &
      & cell_vol, error)
   if (present(vol)) vol = cell_vol
   if (allocated(error)) return

   if (ndim == 0) then
      allocate(trans(3, 0))
      return
   end if

   ! Determine reciprocal-space cutoff radius for each periodic dimensionality
   tol_val = epsilon(1.0_wp)
   select case (ndim)
   case (1)
      u_cut = sqrt(max(42.25_wp, -log(tol_val) &
         & + log(max(1.0_wp, 4.0_wp * max(1.0_wp, alpha**2) / cell_vol)) + 4.0_wp))
      gcut = 2.0_wp * alpha * u_cut
      do while (get_rec_term_1d(gcut, alpha, cell_vol) > tol_val)
         gcut = 1.1_wp * gcut
      end do
   case (2)
      u_cut = sqrt(max(42.25_wp, -log(tol_val) &
         & + log(max(1.0_wp, 8.0_wp * sqrt(pi) * max(alpha, alpha**3) / cell_vol)) &
         & + 4.0_wp))
      gcut = 2.0_wp * alpha * u_cut
      do while (get_rec_term_2d(gcut, alpha, cell_vol) > tol_val)
         gcut = 1.1_wp * gcut
      end do
   case (3)
      u_cut = sqrt(max(42.25_wp, -log(tol_val) &
         & + log(max(1.0_wp, 16.0_wp * pi * max(1.0_wp, alpha**2) / cell_vol)) &
         & + 4.0_wp))
      gcut = 2.0_wp * alpha * u_cut
      do while (get_rec_term_3d(gcut, alpha, cell_vol) > tol_val)
         gcut = 1.1_wp * gcut
      end do
   case default
      gcut = 0.0_wp
   end select

   ! Project reciprocal cutoff sphere onto real-space basis: |m_i| <= gcut * |u_i| / (2*pi)
   rep = 0
   do i = 1, ndim
      ax = p_axes(i)
      rep(ax) = max(2, ceiling(gcut * norm2(red_lat(:, ax)) / twopi))
   end do

   call get_lattice_points(rec_lat, rep, .false., trans, error)
   if (allocated(error)) return

end subroutine get_ewald_rec_trans


!> Prepare reduced direct and reciprocal lattice vectors for a periodic subspace
subroutine prepare_periodic_lattice(periodic, lattice, ndim, p_axes, red_lat, &
   & rec_lat, vol, error)

   !> Periodicity mask along each lattice axis
   logical, intent(in) :: periodic(:)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Number of periodic dimensions
   integer, intent(out) :: ndim

   !> Active periodic axis indices
   integer, intent(out) :: p_axes(3)

   !> Size-reduced real-space lattice vectors (3, 3)
   real(wp), intent(out) :: red_lat(3, 3)

   !> Reciprocal-space lattice vectors (3, 3)
   real(wp), intent(out) :: rec_lat(3, 3)

   !> Generalized cell measure (length in 1D, area in 2D, volume in 3D)
   real(wp), intent(out) :: vol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: basis(3, 3), ginv(3, 3)

   red_lat = 0.0_wp
   rec_lat = 0.0_wp
   vol = 0.0_wp

   call get_periodic_axes(periodic, ndim, p_axes)

   call get_inverse_metric(periodic, lattice, basis, ginv, error)
   if (allocated(error)) return

   ! Reduce the periodic basis to shorten direct and reciprocal box bounds
   call reduce_basis(ndim, p_axes, basis)

   ! Compute inverse metric G^(-1) and reciprocal basis B = 2*pi * L * G^(-1)
   call get_inverse_metric(periodic, basis, red_lat, ginv, error)
   if (allocated(error)) return

   rec_lat = twopi * matmul(red_lat, ginv)

   select case (ndim)
   case (1)
      vol = norm2(red_lat(:, p_axes(1)))
   case (2)
      vol = norm2(crossprod(red_lat(:, p_axes(1)), red_lat(:, p_axes(2))))
   case (3)
      vol = abs(matdet_3x3(red_lat))
   case default
      vol = 0.0_wp
   end select

end subroutine prepare_periodic_lattice


!> Reduce periodic lattice basis by iterative pairwise size reduction
pure subroutine reduce_basis(ndim, p_axes, lat)

   !> Number of periodic dimensions
   integer, intent(in) :: ndim

   !> Active periodic axes
   integer, intent(in) :: p_axes(3)

   !> Lattice vectors, reduced in-place on periodic columns
   real(wp), intent(inout) :: lat(3, 3)

   integer :: iter, i, j, ai, aj
   real(wp) :: ni, nj, dot_ij, mu, tmp(3)
   logical :: changed

   if (ndim < 2) return

   ! Integer shears preserve the lattice while shortening the translation search
   do iter = 1, 100
      changed = .false.
      do i = 1, ndim - 1
         do j = i + 1, ndim
            ai = p_axes(i)
            aj = p_axes(j)
            ni = dot_product(lat(:, ai), lat(:, ai))
            nj = dot_product(lat(:, aj), lat(:, aj))
            if (ni > nj * (1.0_wp + 1.0e-12_wp)) then
               tmp = lat(:, ai)
               lat(:, ai) = lat(:, aj)
               lat(:, aj) = tmp
               ni = nj
               changed = .true.
            end if
            if (ni <= tiny(1.0_wp)) cycle

            dot_ij = dot_product(lat(:, ai), lat(:, aj))
            mu = anint(dot_ij / ni)
            if (abs(mu) >= 1.0_wp &
               & .and. abs(dot_ij) > 0.5_wp * ni * (1.0_wp + 1.0e-12_wp)) then
               lat(:, aj) = lat(:, aj) - mu * lat(:, ai)
               changed = .true.
            end if
         end do
      end do
      if (.not. changed) exit
   end do

end subroutine reduce_basis


!> Get optimal alpha-parameter for the Ewald summation by finding alpha where
!> decline of real and reciprocal part of Ewald are equal.
subroutine search_alpha(lattice, rec_lat, volume, tolerance, alpha, get_rec_term)

   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Reciprocal vectors
   real(wp), intent(in) :: rec_lat(:, :)

   !> Volume (or area/length) of the unit cell
   real(wp), intent(in) :: volume

   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance

   !> Optimal alpha
   real(wp), intent(out) :: alpha

   !> Reciprocal term evaluator
   procedure(get_rec_term_gen) :: get_rec_term

   real(wp), parameter :: alpha0 = 1.0e-8_wp
   integer, parameter :: niter = 30

   real(wp) :: alpl, alpr, rlen, dlen, diff
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:, :)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:, :)**2, dim=1)))

   ! Balance decay of real- and reciprocal-space terms; fall back on failure
   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      stat = 1
   else if (alpha <= alpha0) then
      stat = 2
   end if

   if (stat == 0) then
      alpl = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0.0_wp) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
         ibs = ibs + 1
      end do
      if (ibs > niter) then
         stat = 4
      end if
   end if

   if (stat /= 0) then
      alpha = 0.25_wp
   end if

end subroutine search_alpha


!> Returns the difference in the decrease of the real and reciprocal parts of the
!> Ewald sum. In order to make the real space part shorter than the reciprocal
!> space part, the values are taken at different distances for the real and the
!> reciprocal space parts.
pure function rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume) result(diff)

   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha

   !> Procedure pointer to reciprocal routine
   procedure(get_rec_term_gen) :: get_rec_term

   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: rlen

   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: dlen

   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume

   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = (get_rec_term(4.0_wp * rlen, alpha, volume) &
      & - get_rec_term(5.0_wp * rlen, alpha, volume)) &
      & - (get_dir_term(2.0_wp * dlen, alpha) - get_dir_term(3.0_wp * dlen, alpha))

end function rec_dir_diff


!> Returns the maximum value of a term in the real space part of the Ewald summation
!> for a given vector length.
pure function get_dir_term(rr, alpha) result(dval)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Real space term
   real(wp) :: dval

   dval = erfc(alpha * rr) / rr

end function get_dir_term


!> Returns the maximum value of a term in the 3D reciprocal space part of the Ewald
!> summation for a given vector length.
pure function get_rec_term_3d(gg, alpha, vol) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp * pi * (exp(-0.25_wp * gg * gg / (alpha**2)) / (vol * gg * gg))

end function get_rec_term_3d


!> Returns the maximum value of a term in the 2D reciprocal space part of the Ewald
!> summation for a given vector length (at z = 0).
pure function get_rec_term_2d(gg, alpha, area) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Area of the 2D real space unit cell
   real(wp), intent(in) :: area

   !> Reciprocal term
   real(wp) :: rval

   rval = twopi * erfc(0.5_wp * gg / alpha) / (area * gg)

end function get_rec_term_2d


!> Returns the maximum value of a term in the 1D reciprocal space part of the Ewald
!> summation for a given vector length (at rho = 0).
pure function get_rec_term_1d(gg, alpha, len1) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Length of the 1D real space unit cell
   real(wp), intent(in) :: len1

   !> Reciprocal term
   real(wp) :: rval

   real(wp) :: u

   u = 0.25_wp * gg * gg / (alpha * alpha)
   rval = 2.0_wp * exp(-u) / (len1 * (u + 1.0_wp))

end function get_rec_term_1d


end module multicharge_ewald_cutoff
