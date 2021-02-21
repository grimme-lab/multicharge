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
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_wignerseitz
   implicit none
   private

   public :: collect_wignerseitz

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
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
      & new_unittest("wignerseitz-cell-3d", test_wsc_3d) &
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


end module test_wignerseitz
