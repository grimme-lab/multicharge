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

!> Stateless Ewald summation routines for periodic Gaussian charges
module multicharge_ewald
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_ewald_cutoff, only : get_alpha, get_ewald_trans
   use multicharge_ewald_direct, only : get_amat_dir, get_damat_dir, get_hess_dir
   use multicharge_ewald_pbc3d, only : get_rec_weights, get_amat_rec_3d, &
      & get_damat_rec_3d, get_hess_rec_3d
   implicit none
   private

   public :: get_alpha
   public :: get_ewald_amat, get_ewald_damat, add_ewald_hessian

   real(wp), parameter :: sqrtpi = sqrt(pi)


contains


!> Evaluate the periodic Gaussian Ewald interaction matrix
subroutine get_ewald_amat(mol, wsc, alpha, rad, amat)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Atomic Gaussian charge radii indexed by species ID
   real(wp), intent(in) :: rad(:)

   !> Periodic Ewald Coulomb matrix (filled on 1:nat, 1:nat)
   real(wp), intent(inout) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: vol
   real(wp) :: rij(3), vec(3), gam, wsw, dtmp, rtmp, tmp
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_weight(:)

   call get_ewald_trans(mol%lattice, dtrans, rtrans)
   vol = abs(matdet_3x3(mol%lattice))
   call get_rec_weights(vol, alpha, rtrans, weight=rec_weight)

   ! Every matrix element belongs to one pair; no thread-local matrix is needed.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(amat, mol, wsc, alpha, dtrans, rtrans, rad, rec_weight, npair) &
   !$omp& private(iat, jat, izp, jzp, img, gam, wsw, rij, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(rad(izp)**2 + rad(jzp)**2)
      rij = mol%xyz(:, jat) - mol%xyz(:, iat)
      wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
      tmp = 0.0_wp
      do img = 1, wsc%nimg(jat, iat)
         vec = rij + wsc%trans(:, wsc%tridx(img, jat, iat))
         call get_amat_dir(vec, gam, alpha, dtrans, dtmp)
         call get_amat_rec_3d(vec, rtrans, rec_weight, rtmp)
         tmp = tmp + dtmp + rtmp
      end do
      tmp = tmp * wsw
      if (iat == jat) tmp = tmp - 2.0_wp * alpha / sqrtpi
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

end subroutine get_ewald_amat


!> Evaluate periodic Gaussian Ewald position and strain derivatives
subroutine get_ewald_damat(mol, wsc, alpha, rad, qvec, dadr, dadL, atrace)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Atomic Gaussian charge radii indexed by species ID
   real(wp), intent(in) :: rad(:)

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. atomic positions
   real(wp), intent(out) :: dadr(:, :, :)

   !> Derivative of Coulomb matrix w.r.t. lattice strain
   real(wp), intent(out) :: dadL(:, :, :)

   !> On-site trace of position derivatives
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: vol
   real(wp) :: gam, wsw, rij(3), vec(3)
   real(wp) :: dg(3), ds(3, 3), dgd(3), dsd(3, 3), dgr(3), dsr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)
   real(wp), allocatable :: rec_weight(:), rec_strain(:, :, :)

   call get_ewald_trans(mol%lattice, dtrans, rtrans)
   vol = abs(matdet_3x3(mol%lattice))
   call get_rec_weights(vol, alpha, rtrans, weight=rec_weight, strain=rec_strain)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   if (size(dadr, 3) > mol%nat) then
      dadr(:, :, mol%nat + 1:) = 0.0_wp
   end if

   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(mol, wsc, alpha, dtrans, rtrans, rad, qvec, dadr, &
   !$omp& rec_weight, rec_strain, npair) &
   !$omp& private(iat, jat, izp, jzp, img, gam, wsw, rij, vec, &
   !$omp& dg, ds, dgd, dsd, dgr, dsr) &
   !$omp& reduction(+:atrace, dadL)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(rad(izp)**2 + rad(jzp)**2)
      rij = mol%xyz(:, jat) - mol%xyz(:, iat)
      wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
      dg(:) = 0.0_wp
      ds(:, :) = 0.0_wp
      do img = 1, wsc%nimg(jat, iat)
         vec = rij + wsc%trans(:, wsc%tridx(img, jat, iat))
         call get_damat_dir(vec, gam, alpha, dtrans, dgd, dsd)
         call get_damat_rec_3d(vec, rtrans, rec_weight, rec_strain, dgr, dsr)
         dg = dg + dgd + dgr
         ds = ds + dsd + dsr
      end do
      dg = dg * wsw
      ds = ds * wsw
      if (iat == jat) then
         dadr(:, iat, iat) = 0.0_wp
         dadL(:, :, iat) = dadL(:, :, iat) + ds * qvec(iat)
      else
         atrace(:, iat) = atrace(:, iat) - dg * qvec(jat)
         atrace(:, jat) = atrace(:, jat) + dg * qvec(iat)
         dadr(:, iat, jat) = -dg * qvec(iat)
         dadr(:, jat, iat) = +dg * qvec(jat)
         dadL(:, :, jat) = dadL(:, :, jat) + ds * qvec(iat)
         dadL(:, :, iat) = dadL(:, :, iat) + ds * qvec(jat)
      end if
   end do
   !$omp end parallel do

end subroutine get_ewald_damat


!> Accumulate periodic Gaussian Ewald Cartesian Hessian blocks
subroutine add_ewald_hessian(mol, wsc, alpha, rad, qvec, adj, hessian, diagonal)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Atomic Gaussian charge radii indexed by species ID
   real(wp), intent(in) :: rad(:)

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Adjoint vector
   real(wp), intent(in) :: adj(:)

   !> Off-diagonal Cartesian Hessian blocks
   real(wp), intent(inout) :: hessian(:, :)

   !> Diagonal Cartesian Hessian blocks
   real(wp), intent(inout) :: diagonal(:, :, :)

   integer :: iat, jat, izp, jzp, ii, jj, img, ipair, npair
   real(wp) :: vol
   real(wp) :: rij(3), vec(3), box(3, 3), dbox(3, 3), rbox(3, 3), gam, weight, wsw
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_hess(:, :, :)

   call get_ewald_trans(mol%lattice, dtrans, rtrans)
   vol = abs(matdet_3x3(mol%lattice))
   call get_rec_weights(vol, alpha, rtrans, hess=rec_hess)

   ! Each pair owns its off-diagonal blocks. Only the diagonal blocks are
   ! reduced, avoiding a full (3N, 3N) Hessian copy for every OpenMP thread.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(mol, wsc, alpha, dtrans, rtrans, rad, qvec, adj, hessian, &
   !$omp& rec_hess, npair) &
   !$omp& private(iat, jat, izp, jzp, ii, jj, img, rij, vec, box, dbox, rbox, &
   !$omp& gam, weight, wsw) &
   !$omp& reduction(+:diagonal)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      if (iat == jat) cycle

      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(rad(izp)**2 + rad(jzp)**2)
      weight = -(adj(iat) * qvec(jat) + adj(jat) * qvec(iat))
      rij = mol%xyz(:, jat) - mol%xyz(:, iat)
      wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
      box = 0.0_wp
      do img = 1, wsc%nimg(jat, iat)
         vec = rij + wsc%trans(:, wsc%tridx(img, jat, iat))
         call get_hess_dir(vec, gam, alpha, dtrans, dbox)
         call get_hess_rec_3d(vec, rtrans, rec_hess, rbox)
         box = box + dbox + rbox
      end do
      box = (weight * wsw) * box
      ii = 3 * (iat - 1) + 1
      jj = 3 * (jat - 1) + 1
      diagonal(:, :, iat) = diagonal(:, :, iat) + box
      diagonal(:, :, jat) = diagonal(:, :, jat) + box
      hessian(ii:ii + 2, jj:jj + 2) = hessian(ii:ii + 2, jj:jj + 2) - box
      hessian(jj:jj + 2, ii:ii + 2) = hessian(jj:jj + 2, ii:ii + 2) - box
   end do
   !$omp end parallel do

end subroutine add_ewald_hessian


end module multicharge_ewald
