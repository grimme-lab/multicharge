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

!> @file multicharge/charge.f90
!> Contains functions to calculate the partial charges with and without
!> separate charge model setup

!> Interface to the charge models
module multicharge_charge
   use mctc_cutoff, only : get_lattice_points
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use multicharge_model, only : mchrg_model_type
   use multicharge_model_cache, only : cache_container
   use multicharge_param, only : new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private

   public :: get_charges, get_eeq_charges, get_eeqbc_charges, add_charge_hessian


contains


!> Classical electronegativity equilibration charges
subroutine get_charges(mchrg_model, mol, error, qvec, dqdr, dqdL, cache)

   !> Multicharge model
   class(mchrg_model_type), intent(in) :: mchrg_model

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   !> Retain the solution and Cartesian response for add_charge_hessian.
   !> Recompute after any change to the geometry, charge, or model.
   type(cache_container), intent(inout), optional :: cache

   logical :: grad
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: trans(:, :)

   grad = (present(dqdr) .and. present(dqdL)) .or. present(cache)

   allocate(cn(mol%nat), qloc(mol%nat))
   if (grad) then
      allocate(dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat))
      allocate (dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat))
   end if

   call get_lattice_points(mol%periodic, mol%lattice, mchrg_model%ncoord%cutoff, trans)
   call mchrg_model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call mchrg_model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call mchrg_model%solve(mol, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
      & qvec=qvec, dqdr=dqdr, dqdL=dqdL, response=cache)

end subroutine get_charges


!> Add sum_i weights(i) d2q_i/dR2, holding weights fixed in the contraction.
!> The cache must come from get_charges with the same structure and model.
subroutine add_charge_hessian(mchrg_model, mol, error, weights, hessian, cache)
   class(mchrg_model_type), intent(in) :: mchrg_model
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   real(wp), intent(in) :: weights(:)
   real(wp), intent(inout) :: hessian(:, :)
   type(cache_container), intent(inout) :: cache

   call mchrg_model%add_charge_hessian(mol, error, weights, hessian, cache)
end subroutine add_charge_hessian


!> Obtain charges from electronegativity equilibration model
subroutine get_eeq_charges(mol, error, qvec, dqdr, dqdL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   class(mchrg_model_type), allocatable :: eeq_model

   call new_eeq2019_model(mol, eeq_model, error)

   call get_charges(eeq_model, mol, error, qvec, dqdr, dqdL)

end subroutine get_eeq_charges


!> Obtain charges from bond capacity electronegativity equilibration model
subroutine get_eeqbc_charges(mol, error, qvec, dqdr, dqdL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic partial charges
   real(wp), intent(out), contiguous :: qvec(:)

   !> Derivative of the partial charges w.r.t. the Cartesian coordinates
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)

   !> Derivative of the partial charges w.r.t. strain deformations
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

   class(mchrg_model_type), allocatable :: eeqbc_model

   call new_eeqbc2025_model(mol, eeqbc_model, error)

   call get_charges(eeqbc_model, mol, error, qvec, dqdr, dqdL)

end subroutine get_eeqbc_charges


end module multicharge_charge
