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
!> Ewald splitting parameter and lattice translation generation

!> Ewald splitting parameter search and direct/reciprocal lattice translations
module multicharge_ewald_cutoff
   use mctc_cutoff, only : get_lattice_points
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   implicit none
   private

   public :: get_alpha, get_ewald_trans

   real(wp), parameter :: twopi = 2.0_wp * pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))


contains


!> Determine the optimal Ewald splitting parameter alpha
subroutine get_alpha(lattice, alpha, error)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Ewald splitting parameter
   real(wp), intent(out) :: alpha

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: vol, rec_lat(3, 3)

   vol = abs(matdet_3x3(lattice))
   if (vol < eps) then
      call fatal_error(error, "Lattice vectors are singular or degenerate")
      return
   end if
   rec_lat = twopi * transpose(matinv_3x3(lattice))

   call search_alpha(lattice, rec_lat, vol, eps, alpha)

end subroutine get_alpha


!> Generate direct-space and half-space reciprocal Ewald lattice vectors
pure subroutine get_ewald_trans(lattice, dtrans, rtrans)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Direct-space lattice translations (3, ndir)
   real(wp), allocatable, intent(out) :: dtrans(:, :)

   !> Open half-space reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: rtrans(:, :)

   integer, parameter :: rep(3) = [2, 2, 2]

   integer :: nrec, itr, ix, iy, iz, jy, jz
   real(wp) :: rec_lat(3, 3)
   type(error_type), allocatable :: error

   call get_lattice_points(lattice, rep, .true., dtrans, error)

   rec_lat = twopi * transpose(matinv_3x3(lattice))

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


!> Find alpha balancing the decay rates of the direct and reciprocal Ewald sums
subroutine search_alpha(lattice, rec_lat, volume, tolerance, alpha)

   !> Real-space lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Reciprocal-space lattice vectors
   real(wp), intent(in) :: rec_lat(:, :)

   !> Volume of the unit cell
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

   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, rlen, dlen, volume)
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
         diff = rec_dir_diff(alpha, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0.0_wp) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, rlen, dlen, volume)
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
pure function rec_dir_diff(alpha, rlen, dlen, volume) result(diff)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Shortest reciprocal vector norm
   real(wp), intent(in) :: rlen

   !> Shortest real-space lattice vector norm
   real(wp), intent(in) :: dlen

   !> Unit cell volume
   real(wp), intent(in) :: volume

   !> Decay-rate difference
   real(wp) :: diff

   diff = (get_rec_term(4.0_wp * rlen, alpha, volume) &
      & - get_rec_term(5.0_wp * rlen, alpha, volume)) &
      & - (get_dir_term(2.0_wp * dlen, alpha) - get_dir_term(3.0_wp * dlen, alpha))

end function rec_dir_diff


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


!> Maximum 3D reciprocal-space Ewald potential term at norm gg
pure function get_rec_term(gg, alpha, vol) result(rval)

   !> Reciprocal-space vector norm
   real(wp), intent(in) :: gg

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Unit cell volume
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp * pi * exp(-0.25_wp * gg * gg / (alpha * alpha)) / (vol * gg * gg)

end function get_rec_term


end module multicharge_ewald_cutoff
