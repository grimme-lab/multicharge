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

!> @file multicharge/model/cache.f90
!> Contains the cache container for mutable data of the charge models

!> Cache for charge models
module multicharge_model_cache
   use mctc_env, only : wp
   use mctc_wignerseitz, only : wignerseitz_cell
   implicit none
   private

   public :: mchrg_cache, ensure_shape


   !> Cache for the charge model
   type :: mchrg_cache

      !> Coordination number array
      real(wp), allocatable :: cn(:)

      !> Ewald separation parameter
      real(wp) :: alpha

      !> Wigner-Seitz cell
      type(wignerseitz_cell) :: wsc

      !> Translation matrix
      real(wp), allocatable :: trans(:, :)

      !> Local charges
      real(wp), allocatable :: qloc(:)

      !> Full Maxwell capacitance matrix
      real(wp), allocatable :: cmat(:, :)

      !> Store tmp array from xvec calculation for reuse
      real(wp), allocatable :: xtmp(:)

      !> Electronegativity vector
      real(wp), allocatable :: xvec(:)

      !> Coulomb matrix
      real(wp), allocatable :: amat(:, :)

      !> Inversed amat
      real(wp), allocatable :: ainv(:, :)

      !> Solution of the ES equation
      real(wp), allocatable :: vrhs(:)

      !> Coordination number gradient w.r.t the positions
      real(wp), allocatable :: dcndr(:, :, :)

      !> Coordination number gradient w.r.t the lattice vectors
      real(wp), allocatable :: dcndL(:, :, :)

      !> Local charge derivatives w.r.t positions
      real(wp), allocatable :: dqlocdr(:, :, :)

      !> Local charge derivatives w.r.t lattice vectors
      real(wp), allocatable :: dqlocdL(:, :, :)

      !> Derivative of Maxwell capacitance matrix w.r.t positions
      real(wp), allocatable :: dcdr(:, :, :)

      !> Derivative of Maxwell capacitance matrix w.r.t positions, compressed
      real(wp), allocatable :: dcdrdiag(:, :)

      !> Derivative of Maxwell capacitance matrix w.r.t lattice vectors
      real(wp), allocatable :: dcdL(:, :, :)

      !> Coulomb matrix derivatives w.r.t positions
      real(wp), allocatable :: dadr(:, :, :)

      !> Coulomb matrix derivatives w.r.t lattice vectors
      real(wp), allocatable :: dadL(:, :, :)

      !> Electronegativity derivatives w.r.t positions
      real(wp), allocatable :: dxdr(:, :, :)

      !> Electronegativity derivatives w.r.t lattice vectors
      real(wp), allocatable :: dxdL(:, :, :)

      !> Charge derivatives w.r.t positions retained for the charge Hessian
      real(wp), allocatable :: dqdr(:, :, :)

      !> Logical flag for gradient calculation
      logical :: grad

   end type mchrg_cache

   !> Allocate a cache array with the requested shape, reallocating on mismatch
   interface ensure_shape
      module procedure :: ensure_shape_1d
      module procedure :: ensure_shape_2d
      module procedure :: ensure_shape_3d
   end interface ensure_shape


contains


!> Allocate a rank-1 array, reallocate it if its size differs
subroutine ensure_shape_1d(array, n1)
   !> Array to (re)allocate, contents are undefined after reallocation
   real(wp), allocatable, intent(inout) :: array(:)
   !> Requested size
   integer, intent(in) :: n1

   if (allocated(array)) then
      if (size(array) == n1) return
      deallocate(array)
   end if
   allocate(array(n1))
end subroutine ensure_shape_1d


!> Allocate a rank-2 array, reallocate it if its shape differs
subroutine ensure_shape_2d(array, n1, n2)
   !> Array to (re)allocate, contents are undefined after reallocation
   real(wp), allocatable, intent(inout) :: array(:, :)
   !> Requested extent of the first dimension
   integer, intent(in) :: n1
   !> Requested extent of the second dimension
   integer, intent(in) :: n2

   if (allocated(array)) then
      if (all(shape(array) == [n1, n2])) return
      deallocate(array)
   end if
   allocate(array(n1, n2))
end subroutine ensure_shape_2d


!> Allocate a rank-3 array, reallocate it if its shape differs
subroutine ensure_shape_3d(array, n1, n2, n3)
   !> Array to (re)allocate, contents are undefined after reallocation
   real(wp), allocatable, intent(inout) :: array(:, :, :)
   !> Requested extent of the first dimension
   integer, intent(in) :: n1
   !> Requested extent of the second dimension
   integer, intent(in) :: n2
   !> Requested extent of the third dimension
   integer, intent(in) :: n3

   if (allocated(array)) then
      if (all(shape(array) == [n1, n2, n3])) return
      deallocate(array)
   end if
   allocate(array(n1, n2, n3))
end subroutine ensure_shape_3d

end module multicharge_model_cache
