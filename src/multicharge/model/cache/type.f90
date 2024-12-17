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

!> @file multicharge/model/cache/type.f90
!> Contains the cache baseclass for the charge models

!> Cache for charge models
module multicharge_model_cache
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use multicharge_wignerseitz, only: wignerseitz_cell_type
   implicit none
   private

   !> Cache for the charge model
   type, abstract, public :: mchrg_cache
      !> CN array
      real(wp), allocatable :: cn(:)
      !> Gradients
      real(wp), allocatable :: dcndr(:, :, :)
      real(wp), allocatable :: dcndL(:, :, :)
      real(wp) :: alpha
      type(wignerseitz_cell_type) :: wsc
   contains
      !>
      procedure(update), deferred :: update
   end type mchrg_cache

   abstract interface
      subroutine update(self, mol)
         import mchrg_cache, structure_type
         class(mchrg_cache), intent(inout) :: self
         type(structure_type), intent(in) :: mol
      end subroutine update

   end interface

end module multicharge_model_cache
