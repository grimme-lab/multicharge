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
   use mctc_cutoff, only : get_inverse_metric, get_lattice_points, &
      & get_periodic_axes, reduce_basis
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : crossprod, matdet_3x3
   implicit none
   private

   public :: get_alpha, get_ewald_trans

   real(wp), parameter :: twopi = 2.0_wp * pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))


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
   call search_alpha(ndim, red_lat(:, p_axes(:ndim)), rlat(:, p_axes(:ndim)), &
      & volume, eps, alpha)

   ! Prevent small alpha and excessive direct-space boxes in anisotropic 3D cells
   if (ndim == 3 .and. volume > eps) then
      alpha = max(alpha, sqrt(pi) / volume**(1.0_wp / 3.0_wp))
   end if

   if (present(vol)) vol = volume
   if (present(rec_lat)) rec_lat = rlat

end subroutine get_alpha


!> Generate converged direct-space and half-space reciprocal Ewald vectors
subroutine get_ewald_trans(periodic, lattice, alpha, gam_min, dtrans, rtrans, error)

   !> Periodicity mask along each lattice axis
   logical, intent(in) :: periodic(:)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Minimum Gaussian exponent min_ij(gamma_ij)
   real(wp), intent(in) :: gam_min

   !> Direct-space lattice translations (3, ndir)
   real(wp), allocatable, intent(out) :: dtrans(:, :)

   !> Open half-space reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: rtrans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: tol_val = epsilon(1.0_wp)
   real(wp), parameter :: tol_wsc = 0.3_wp
   real(wp), parameter :: max_cut = sqrt(huge(1.0_wp))

   integer :: ndim, p_axes(3), rep(3), i, ax, nrec, itr, ix, iy, iz, jy, jz
   real(wp) :: red_lat(3, 3), rec_lat(3, 3), vol
   real(wp) :: kappa, lo, hi, rcut, len_i, min_len, r_disp, gcut, npts

   call prepare_periodic_lattice(periodic, lattice, ndim, p_axes, red_lat, rec_lat, &
      & vol, error)
   if (allocated(error)) return

   ! Bracket and bisect real-space cutoff for potential and second derivatives
   kappa = min(alpha, gam_min)
   if (.not. (alpha > 0.0_wp .and. kappa > 0.0_wp)) then
      call fatal_error(error, "Ewald exponents must be positive")
      return
   end if

   lo = 0.0_wp
   hi = 1.0_wp / kappa
   do while (dir_bound(hi, kappa) > tol_val .and. hi < max_cut)
      lo = hi
      hi = 2.0_wp * hi
   end do
   do i = 1, 20
      rcut = 0.5_wp * (lo + hi)
      if (dir_bound(rcut, kappa) > tol_val) then
         lo = rcut
      else
         hi = rcut
      end if
   end do
   rcut = hi

   ! Add bound on parallel Wigner-Seitz pair displacements (including competing images)
   r_disp = 0.0_wp
   min_len = huge(1.0_wp)
   do i = 1, ndim
      len_i = norm2(red_lat(:, p_axes(i)))
      r_disp = r_disp + 0.5_wp * len_i
      min_len = min(min_len, len_i)
   end do
   r_disp = sqrt(max(r_disp, min_len)**2 + tol_wsc)

   call get_lattice_points(periodic, red_lat, rcut + r_disp, dtrans, error)
   if (allocated(error)) return

   ! Bracket and bisect reciprocal cutoff for potential and second derivatives
   lo = 0.0_wp
   hi = 2.0_wp * alpha
   do while ((1.0_wp + hi**2) * get_rec_term(ndim, hi, alpha, vol) > tol_val &
      & .and. hi < max_cut)
      lo = hi
      hi = 2.0_wp * hi
   end do
   do i = 1, 20
      gcut = 0.5_wp * (lo + hi)
      if ((1.0_wp + gcut**2) * get_rec_term(ndim, gcut, alpha, vol) > tol_val) then
         lo = gcut
      else
         hi = gcut
      end if
   end do
   gcut = hi

   ! Project reciprocal cutoff sphere onto real-space basis: |m_i| <= gcut * |u_i| / (2*pi)
   rep = 0
   do i = 1, ndim
      ax = p_axes(i)
      len_i = gcut * norm2(red_lat(:, ax)) / twopi
      if (.not. (len_i < real(huge(nrec), wp))) then
         call fatal_error(error, "Number of lattice points exceeds the integer range")
         return
      end if
      rep(ax) = max(2, ceiling(len_i))
   end do

   npts = 0.5_wp * (product(2.0_wp * real(rep, wp) + 1.0_wp) - 1.0_wp)
   if (.not. (npts <= real(huge(nrec), wp))) then
      call fatal_error(error, "Number of lattice points exceeds the integer range")
      return
   end if

   ! Keep one representative per (+G, -G) antipodal pair (first nonzero index > 0)
   nrec = (product(2 * rep + 1) - 1) / 2
   allocate(rtrans(3, nrec))
   itr = 0
   do ix = 0, rep(1)
      do iy = 0, rep(2)
         do iz = 0, rep(3)
            if (ix == 0 .and. iy == 0 .and. iz == 0) cycle
            do jy = 1, merge(-1, 1, ix > 0 .and. iy > 0), -2
               do jz = 1, merge(-1, 1, (ix > 0 .or. iy > 0) .and. iz > 0), -2
                  itr = itr + 1
                  rtrans(:, itr) = rec_lat(:, 1) * real(ix, wp) &
                     & + rec_lat(:, 2) * real(iy * jy, wp) &
                     & + rec_lat(:, 3) * real(iz * jz, wp)
               end do
            end do
         end do
      end do
   end do

end subroutine get_ewald_trans


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

   integer :: i
   real(wp) :: basis(3, 3), ginv(3, 3)

   red_lat = 0.0_wp
   rec_lat = 0.0_wp
   vol = 0.0_wp

   call get_periodic_axes(periodic, ndim, p_axes)
   if (ndim == 0) then
      call fatal_error(error, "Ewald summation requires at least one periodic dimension")
      return
   end if

   if (size(lattice, 1) < 3 .or. size(lattice, 2) < 3) then
      call fatal_error(error, "Periodic lattice vectors must have shape (3, 3)")
      return
   end if

   basis = 0.0_wp
   do i = 1, ndim
      basis(:, p_axes(i)) = lattice(:3, p_axes(i))
   end do

   ! Reduce the periodic basis to shorten direct and reciprocal box bounds
   call reduce_basis(basis)

   ! Validate linear independence and compute reciprocal basis B = 2*pi * L * G^(-1)
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


!> Find alpha balancing the decay rates of the direct and reciprocal Ewald sums
subroutine search_alpha(ndim, lattice, rec_lat, volume, tolerance, alpha)

   !> Number of periodic dimensions
   integer, intent(in) :: ndim

   !> Reduced periodic lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Reciprocal periodic vectors
   real(wp), intent(in) :: rec_lat(:, :)

   !> Volume (or area/length) of the unit cell
   real(wp), intent(in) :: volume

   !> Tolerance for difference in real and reciprocal decay
   real(wp), intent(in) :: tolerance

   !> Optimal Ewald splitting parameter
   real(wp), intent(out) :: alpha

   real(wp), parameter :: alpha0 = 1.0e-8_wp
   integer, parameter :: niter = 30

   real(wp) :: alpl, alpr, rlen, dlen, diff
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:, :)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:, :)**2, dim=1)))

   ! Balance decay of real- and reciprocal-space terms; fall back on failure
   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(ndim, alpha, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(ndim, alpha, rlen, dlen, volume)
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
         diff = rec_dir_diff(ndim, alpha, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(ndim, alpha, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0.0_wp) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(ndim, alpha, rlen, dlen, volume)
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


!> Difference between reciprocal-space and real-space Ewald step decay
pure function rec_dir_diff(ndim, alpha, rlen, dlen, volume) result(diff)

   !> Number of periodic dimensions
   integer, intent(in) :: ndim

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Shortest reciprocal vector norm
   real(wp), intent(in) :: rlen

   !> Shortest real-space lattice vector norm
   real(wp), intent(in) :: dlen

   !> Cell length (1D), area (2D), or volume (3D)
   real(wp), intent(in) :: volume

   !> Decay-rate difference
   real(wp) :: diff

   diff = (get_rec_term(ndim, 4.0_wp * rlen, alpha, volume) &
      & - get_rec_term(ndim, 5.0_wp * rlen, alpha, volume)) &
      & - (get_dir_term(2.0_wp * dlen, alpha) - get_dir_term(3.0_wp * dlen, alpha))

end function rec_dir_diff


!> Direct-space Ewald tail bound covering potential, gradient, and Hessian
pure function dir_bound(rr, kappa) result(dval)

   !> Real-space distance
   real(wp), intent(in) :: rr

   !> Slowest Gaussian exponent
   real(wp), intent(in) :: kappa

   !> Combined real-space tail bound
   real(wp) :: dval

   dval = (1.0_wp + (kappa * rr)**2 * (1.0_wp + kappa**2)) * get_dir_term(rr, kappa)

end function dir_bound


!> Maximum real-space Ewald potential term at distance rr
pure function get_dir_term(rr, alpha) result(dval)

   !> Real-space vector norm
   real(wp), intent(in) :: rr

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Real-space term
   real(wp) :: dval

   dval = erfc(alpha * rr) / rr

end function get_dir_term


!> Maximum 1D, 2D, or 3D reciprocal-space Ewald potential term at norm gg
pure function get_rec_term(ndim, gg, alpha, vol) result(rval)

   !> Number of periodic dimensions
   integer, intent(in) :: ndim

   !> Reciprocal-space vector norm
   real(wp), intent(in) :: gg

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Cell length (1D), area (2D), or volume (3D)
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: rval

   real(wp) :: u

   select case (ndim)
   case (1)
      u = 0.25_wp * gg * gg / (alpha * alpha)
      rval = 2.0_wp * exp(-u) / (vol * (u + 1.0_wp))
   case (2)
      rval = twopi * erfc(0.5_wp * gg / alpha) / (vol * gg)
   case (3)
      rval = 4.0_wp * pi * exp(-0.25_wp * gg * gg / (alpha * alpha)) / (vol * gg * gg)
   case default
      rval = 0.0_wp
   end select

end function get_rec_term


end module multicharge_ewald_cutoff
