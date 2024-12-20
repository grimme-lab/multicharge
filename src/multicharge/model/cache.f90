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
!> Contains the cache baseclass for the charge models and a container for mutable cache data

!> Cache for charge models
module multicharge_model_cache
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use multicharge_wignerseitz, only: wignerseitz_cell_type, new_wignerseitz_cell
   use multicharge_ewald, only: get_alpha
   implicit none
   private

   type, public :: cache_container
      !> Mutable data attribute
      class(*), allocatable :: raw
   end type cache_container

   !> Cache for the charge model
   type, abstract, public :: model_cache
      !> CN array
      real(wp), allocatable :: cn(:)
      !> CN dr gradient
      real(wp), allocatable :: dcndr(:, :, :)
      !> CN dL gradient
      real(wp), allocatable :: dcndL(:, :, :)
      !> Ewald separation parameter
      real(wp) :: alpha
      type(wignerseitz_cell_type) :: wsc
   contains
      !> Create WSC
      procedure :: update
   end type model_cache

contains
   subroutine update(self, mol)
      class(model_cache), intent(inout) :: self
      type(structure_type), intent(in) :: mol

      ! Create WSC
      if (any(mol%periodic)) then
         call new_wignerseitz_cell(self%wsc, mol)
         call get_alpha(mol%lattice, self%alpha)
      end if

   end subroutine update

end module multicharge_model_cache
