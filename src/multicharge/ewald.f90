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
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_ewald_cutoff, only : get_alpha, get_dir_cutoff, get_rec_cutoff, &
      & get_ewald_trans
   use multicharge_ewald_direct, only : get_amat_dir, get_damat_dir, get_hess_dir
   use multicharge_ewald_pbc3d, only : get_rec_weights, get_amat_rec_3d, &
      & get_damat_rec_3d, get_hess_rec_3d
   implicit none
   private

   public :: get_alpha, get_dir_cutoff, get_rec_cutoff, get_ewald_trans
   public :: get_ewald_setup
   public :: get_amat_dir, get_damat_dir, get_hess_dir
   public :: get_rec_weights, get_amat_rec_3d, get_damat_rec_3d, get_hess_rec_3d


contains


!> Prepare converged lattice translations and requested reciprocal-space factors
subroutine get_ewald_setup(mol, wsc, alpha, gam, dtrans, rtrans, weight, &
      & strain, hess)

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

   !> Optional reciprocal-space scalar weights
   real(wp), allocatable, intent(out), optional :: weight(:)

   !> Optional reciprocal-space strain tensors
   real(wp), allocatable, intent(out), optional :: strain(:, :, :)

   !> Optional reciprocal-space Hessian tensors
   real(wp), allocatable, intent(out), optional :: hess(:, :, :)

   integer :: iat, jat, img
   real(wp) :: vol, inv_lat(3, 3), vec(3), smax(3), rmax2

   ! Extent of the pair vectors, which shift the real-space sum off the origin
   inv_lat = matinv_3x3(mol%lattice)
   smax(:) = 0.0_wp
   rmax2 = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, wsc, inv_lat) private(iat, jat, img, vec) &
   !$omp reduction(max:smax, rmax2)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) &
               & + wsc%trans(:, wsc%tridx(img, jat, iat))
            smax(:) = max(smax, abs(matmul(inv_lat, vec)))
            rmax2 = max(rmax2, dot_product(vec, vec))
         end do
      end do
   end do
   !$omp end parallel do

   call get_ewald_trans(mol%lattice, alpha, gam, smax, sqrt(rmax2), dtrans, rtrans)
   vol = abs(matdet_3x3(mol%lattice))
   call get_rec_weights(vol, alpha, rtrans, weight=weight, strain=strain, hess=hess)

end subroutine get_ewald_setup


end module multicharge_ewald
