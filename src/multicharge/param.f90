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

module multicharge_param
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_data, only: get_covalent_rad
   use multicharge_model, only : mchrg_model_type, new_mchrg_model
   use multicharge_param_eeq2019, only : get_eeq_chi, get_eeq_eta, &
      & get_eeq_rad, get_eeq_kcn
   implicit none
   private

   public :: new_eeq2019_model

contains

subroutine new_eeq2019_model(mol, model, error)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Electronegativity equilibration model
   type(mchrg_model_type), intent(out) :: model
   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: cutoff = 25.0_wp, cn_exp = 7.5_wp, cn_max = 8.0_wp

   real(wp), allocatable :: chi(:), eta(:), kcn(:), rad(:), rcov(:)

   chi = get_eeq_chi(mol%num)
   eta = get_eeq_eta(mol%num)
   kcn = get_eeq_kcn(mol%num)
   rad = get_eeq_rad(mol%num)
   rcov = get_covalent_rad(mol%num)

   call new_mchrg_model(model, mol, error=error, chi=chi, &
      & rad=rad, eta=eta, kcn=kcn, cutoff=cutoff, &
      & cn_exp=cn_exp, rcov=rcov, cn_max=cn_max)

end subroutine new_eeq2019_model

end module multicharge_param
