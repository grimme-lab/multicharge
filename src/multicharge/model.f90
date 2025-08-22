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

!> @dir multicharge/model
!> Contains the implementation of the charge models

!> @file multicharge/model.f90
!> Provides a reexport of the charge model implementations

!> Proxy module to reexport the charge model implementations
module multicharge_model
   use multicharge_model_type, only : mchrg_model_type
   use multicharge_model_eeq, only : eeq_model, new_eeq_model
   use multicharge_model_eeqbc, only : eeqbc_model, new_eeqbc_model
   implicit none
   private

   public :: mchrg_model_type
   public :: eeq_model, new_eeq_model
   public :: eeqbc_model, new_eeqbc_model

end module multicharge_model
