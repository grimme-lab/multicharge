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

!> @dir multicharge/ewald
!> Contains the periodic Ewald summation implementation

!> @file multicharge/ewald.f90
!> Provides the periodic Ewald summation interface

!> Model-independent Ewald summation for periodic Gaussian charges.
!>
!> Provides the splitting parameter, lattice translations, and per-vector
!> direct and reciprocal-space kernels. Charge models assemble their own
!> interaction matrices from these building blocks.
module multicharge_ewald
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_math, only : matdet_3x3
   use multicharge_ewald_cutoff, only : get_alpha, get_ewald_trans
   use multicharge_ewald_direct, only : get_amat_dir, get_damat_dir, get_hess_dir
   use multicharge_ewald_pbc3d, only : get_rec_weights, get_amat_rec_3d, &
      & get_damat_rec_3d, get_hess_rec_3d
   implicit none
   private

   public :: get_alpha, get_ewald_trans, get_ewald_setup
   public :: get_amat_dir, get_damat_dir, get_hess_dir
   public :: get_rec_weights, get_amat_rec_3d, get_damat_rec_3d, get_hess_rec_3d


contains


!> Prepare lattice translations and reciprocal-space factors for one evaluation.
!> Only the requested reciprocal-space tensors are allocated.
pure subroutine get_ewald_setup(mol, alpha, dtrans, rtrans, weight, strain, hess)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Direct-space lattice translations (3, ndir)
   real(wp), allocatable, intent(out) :: dtrans(:, :)

   !> Half-space reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: rtrans(:, :)

   !> Optional reciprocal-space scalar weights
   real(wp), allocatable, intent(out), optional :: weight(:)

   !> Optional reciprocal-space strain tensors
   real(wp), allocatable, intent(out), optional :: strain(:, :, :)

   !> Optional reciprocal-space Hessian tensors
   real(wp), allocatable, intent(out), optional :: hess(:, :, :)

   real(wp) :: vol

   call get_ewald_trans(mol%lattice, dtrans, rtrans)
   vol = abs(matdet_3x3(mol%lattice))
   call get_rec_weights(vol, alpha, rtrans, weight=weight, strain=strain, hess=hess)

end subroutine get_ewald_setup


end module multicharge_ewald
