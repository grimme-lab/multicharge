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

!> @file multicharge/cache.f90
!> Contains the cache baseclass for the charge models and a container for mutable cache data

!> Cache for charge models
module multicharge_model_cache
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use multicharge_wignerseitz, only: wignerseitz_cell_type
   implicit none
   private

   type, public :: cache_container
      !> Mutable data attribute
      class(*), allocatable :: raw
   end type cache_container

   !> Cache for the charge model
   type, abstract, public :: model_cache
      !> Coordination number array
      real(wp), allocatable :: cn(:)
      !> Coordination number gradient w.r.t the positions
      real(wp), allocatable :: dcndr(:, :, :)
      !> Coordination number gradient w.r.t the lattice vectors
      real(wp), allocatable :: dcndL(:, :, :)
      !> Ewald separation parameter
      real(wp) :: alpha
      !> Wigner-Seitz cell
      type(wignerseitz_cell_type) :: wsc
   end type model_cache

end module multicharge_model_cache
