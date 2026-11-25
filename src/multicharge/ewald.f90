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

!> Ewald splitting parameter, lattice translations and reciprocal-space kernels
module multicharge_ewald
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_ewald_cutoff, only : get_alpha, get_ewald_trans
   use multicharge_ewald_pbc3d, only : get_rec_weights, get_amat_rec_3d, &
      & get_damat_rec_3d, get_hess_rec_3d
   implicit none
   private

   public :: get_alpha, get_ewald_setup
   public :: get_rec_weights, get_amat_rec_3d, get_damat_rec_3d, get_hess_rec_3d


contains


!> Prepare converged direct and reciprocal lattice translations
subroutine get_ewald_setup(mol, wsc, alpha, gam, dtrans, rtrans, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images of all pairs
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Smallest Gaussian charge exponent of all pairs
   real(wp), intent(in) :: gam

   !> Direct-space lattice translations (3, ndir)
   real(wp), allocatable, intent(out) :: dtrans(:, :)

   !> Half-space reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: rtrans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, jat, img
   real(wp) :: vec(3), rmax2

   ! Longest pair vector, which shifts the real-space sum off the origin
   rmax2 = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, wsc) private(iat, jat, img, vec) &
   !$omp reduction(max:rmax2)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) &
               & + wsc%trans(:, wsc%tridx(img, jat, iat))
            rmax2 = max(rmax2, dot_product(vec, vec))
         end do
      end do
   end do
   !$omp end parallel do

   call get_ewald_trans(mol%lattice, alpha, gam, sqrt(rmax2), dtrans, rtrans, error)

end subroutine get_ewald_setup


end module multicharge_ewald
