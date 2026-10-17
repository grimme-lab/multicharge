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

module test_wignerseitz
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge_wignerseitz, only : wignerseitz_cell_type, new_wignerseitz_cell
   implicit none
   private

   public :: collect_wignerseitz

   real(wp), parameter :: thr = 100 * epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_wignerseitz(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("lattice-points-0d", test_latticepoints_0d), &
      & new_unittest("lattice-points-3d", test_latticepoints_3d), &
      & new_unittest("wignerseitz-cell-0d", test_wsc_0d), &
      & new_unittest("wignerseitz-cell-3d", test_wsc_3d), &
      & new_unittest("wignerseitz-self-images", test_self_images) &
      & ]

end subroutine collect_wignerseitz

subroutine test_latticepoints_0d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: trans(:, :)
   real(wp), parameter :: cutoff = 25.0_wp

   call get_structure(mol, "MB16-43", "01")

   call get_lattice_points(mol%periodic, mol%lattice, thr2, trans)

   call check(error, size(trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(trans, 2), 1)
   if (allocated(error)) return

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   call check(error, size(trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(trans, 2), 1)
   if (allocated(error)) return

end subroutine test_latticepoints_0d

subroutine test_latticepoints_3d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), allocatable :: trans(:, :)
   real(wp), parameter :: cutoff = 25.0_wp

   call get_structure(mol, "ICE10", "ii")

   call get_lattice_points(mol%periodic, mol%lattice, thr2, trans)

   call check(error, size(trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(trans, 2), 27)
   if (allocated(error)) return

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   call check(error, size(trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(trans, 2), 343)
   if (allocated(error)) return

end subroutine test_latticepoints_3d

subroutine test_wsc_0d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell_type) :: wsc

   call get_structure(mol, "MB16-43", "02")

   call new_wignerseitz_cell(wsc, mol)

   call check(error, size(wsc%trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(wsc%trans, 2), 1)
   if (allocated(error)) return

end subroutine test_wsc_0d

subroutine test_wsc_3d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(wignerseitz_cell_type) :: wsc

   call get_structure(mol, "X23", "pyrazole")

   call new_wignerseitz_cell(wsc, mol)

   call check(error, size(wsc%trans, 1), 3)
   if (allocated(error)) return

   call check(error, size(wsc%trans, 2), 27)
   if (allocated(error)) return

end subroutine test_wsc_3d

!> Excluding the zero image must preserve the original translation indices.
subroutine test_self_images(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   type(wignerseitz_cell_type) :: wsc
   real(wp) :: lattice(3, 3), vec(3), total(3)
   integer :: icase, i, expected

   do icase = 1, 2
      lattice = reshape([5.1_wp, 0.0_wp, 0.0_wp, 0.4_wp, 6.2_wp, 0.0_wp, &
         & 0.2_wp, 0.7_wp, 7.1_wp], [3, 3])
      expected = 2
      if (icase == 2) then
         lattice = 0.0_wp
         do i = 1, 3
            lattice(i, i) = 5.1_wp
         end do
         expected = 6
      end if
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), &
         & lattice=lattice, periodic=[.true., .true., .true.])
      call new_wignerseitz_cell(wsc, mol)
      call check(error, wsc%nimg(1, 1), expected)
      if (allocated(error)) return
      call check(error, wsc%nimg_max, expected)
      if (allocated(error)) return
      total = 0.0_wp
      do i = 1, wsc%nimg(1, 1)
         vec = wsc%trans(:, wsc%tridx(i, 1, 1))
         if (abs(norm2(vec)-5.1_wp) > thr) then
            call test_failed(error, "Selected self image is not a closest lattice point")
            return
         end if
         total = total+vec
      end do
      if (maxval(abs(total)) > thr) then
         call test_failed(error, "Closest self images do not occur in opposite pairs")
         return
      end if
   end do
end subroutine test_self_images

end module test_wignerseitz
