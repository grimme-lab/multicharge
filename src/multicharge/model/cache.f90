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
!> Contains the cache for the charge models 

!> Cache for charge models
module multicharge_model_cache
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   implicit none
   private

   !> Cache for the charge model
   type, public :: mchrg_cache
      !> Constraint matrix
      real(wp), allocatable :: cmat(:, :)
      !> Derivative of constraint matrix w.r.t positions
      real(wp), allocatable :: dcdr(:, :, :)
      !> Derivative of constraint matrix w.r.t lattice vectors
      real(wp), allocatable :: dcdL(:, :, :)
   end type mchrg_cache

end module multicharge_model_cache
