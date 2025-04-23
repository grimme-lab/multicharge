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

module multicharge
   use multicharge_charge, only : get_charges, get_eeq_charges, get_eeqbc_charges
   use multicharge_model, only : mchrg_model_type
   use multicharge_output, only : write_ascii_model, write_ascii_properties, &
      & write_ascii_results
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2024_model, mchargeModel
   use multicharge_version, only : get_multicharge_version
   implicit none
   public


end module multicharge
