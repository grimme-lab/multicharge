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
module multicharge_ewald_cutoff
   use mctc_cutoff, only : get_lattice_points, reduce_basis
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use multicharge_ewald_pbc3d, only : get_rec_term_3d
   implicit none
   private

   public :: get_alpha, get_dir_cutoff, get_rec_cutoff, get_ewald_trans

   real(wp), parameter :: twopi = 2.0_wp * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Scaling of the Ewald parameter, shifting work to the reciprocal space
   real(wp), parameter :: alpha_scale = 1.25_wp

   !> Tolerance for a single term of the real- or reciprocal-space sum
   real(wp), parameter :: term_tol = epsilon(0.0_wp)

   abstract interface
      !> Max. value of a real or reciprocal space term for a given vector length
      pure function get_term_gen(dist, alpha, vol) result(val)
         import :: wp

         !> Length of the real or reciprocal space vector
         real(wp), intent(in) :: dist

         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha

         !> Volume of the real space unit cell
         real(wp), intent(in) :: vol

         !> Value of the term
         real(wp) :: val

      end function get_term_gen
   end interface


contains


!> Ewald parameter balancing the cost of the real-space and reciprocal sums
subroutine get_alpha(lattice, alpha, error)

   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Parameter of the Ewald summation
   real(wp), intent(out) :: alpha

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: vol

   vol = abs(matdet_3x3(lattice))
   if (vol < eps) then
      call fatal_error(error, "Lattice vectors are singular or degenerate")
      return
   end if

   alpha = alpha_scale * sqrtpi / vol**(1.0_wp / 3.0_wp)

end subroutine get_alpha


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


!> Reciprocal-space cutoff
pure function get_rec_cutoff(alpha, volume, conv) result(cutoff)

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume

   !> Tolerance for a single reciprocal-space term
   real(wp), intent(in) :: conv

   !> Reciprocal-space cutoff
   real(wp) :: cutoff

   cutoff = search_cutoff(get_rec_term_3d, alpha, volume, conv)

end function get_rec_cutoff


!> Direct and half-space reciprocal translations within the single-term cutoffs
pure subroutine get_ewald_trans(lattice, alpha, gam, rmax, dtrans, rtrans, error)

   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Smallest Gaussian charge exponent of all pairs
   real(wp), intent(in) :: gam

   !> Largest length of the pair vectors
   real(wp), intent(in) :: rmax

   !> Real-space lattice translations
   real(wp), allocatable, intent(out) :: dtrans(:, :)

   !> Reciprocal lattice translations in a half space
   real(wp), allocatable, intent(out) :: rtrans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: rep(3), nrec, itr, ix, iy, iz, jy, jz
   real(wp) :: lat(3, 3), rec_lat(3, 3), rcut, gcut, vec(3), extent(3)
   real(wp), allocatable :: tmp(:, :)

   ! The real-space kernel decays with the smaller of alpha and gam
   rcut = get_dir_cutoff(min(alpha, gam), term_tol)

   ! Each half-space vector carries its +G/-G pair, twice a single term
   gcut = get_rec_cutoff(alpha, abs(matdet_3x3(lattice)), 0.5_wp * term_tol)

   ! A shorter, more orthogonal basis of the same lattice brings the enumerated
   ! boxes in both spaces closer to the cutoff spheres, the selected vectors
   ! are the same for any basis
   lat(:, :) = lattice(:3, :3)
   call reduce_basis(lat)

   ! Pair vectors shift the real-space sum off the origin by at most rmax
   call get_lattice_points([.true., .true., .true.], lat, rcut + rmax, tmp, error)
   if (allocated(error)) return
   dtrans = tmp(:, pack([(itr, itr = 1, size(tmp, 2))], &
      & norm2(tmp, dim=1) <= rcut + rmax))
   deallocate(tmp)

   ! Reciprocal vectors within the cutoff sphere, one of each +G/-G pair
   rec_lat = twopi * transpose(matinv_3x3(lat))
   extent(:) = gcut * norm2(lat, dim=1) / twopi
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

end subroutine get_ewald_trans


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
pure function search_cutoff(get_term, alpha, volume, conv) result(cutoff)

   !> Procedure to evaluate the Ewald term
   procedure(get_term_gen) :: get_term

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume

   !> Tolerance for a single term
   real(wp), intent(in) :: conv

   !> Cutoff for the vector length
   real(wp) :: cutoff

   integer, parameter :: niter = 30
   real(wp) :: xl, xr, x
   integer :: iter

   xr = eps
   do while (get_term(xr, alpha, volume) > conv .and. xr <= huge(1.0_wp))
      xr = 2.0_wp * xr
   end do
   xl = 0.5_wp * xr

   do iter = 1, niter
      if (xr - xl <= eps * xr) exit
      x = 0.5_wp * (xl + xr)
      if (get_term(x, alpha, volume) > conv) then
         xl = x
      else
         xr = x
      end if
   end do
   cutoff = xr

end function search_cutoff


!> Returns the max. value of a term in the real space
pure function get_dir_term(rr, alpha, vol) result(dval)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell (unused)
   real(wp), intent(in) :: vol

   !> Real space term
   real(wp) :: dval

   dval = erfc(alpha * rr) / rr

end function get_dir_term


end module multicharge_ewald_cutoff
