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
!> Contains the 1D, 2D, and 3D Ewald summation implementations
!>
!> The 1D and 2D Ewald decompositions follow:
!> - F. Nestler and D. Potts, *Fast Ewald summation under 2d- and 1d-periodic
!>   boundary conditions based on NFFTs*, Proc. Sampling Theory and Applications
!>   (SampTA), 329-332 (2013).
!> - F. Nestler, M. Pippig, and D. Potts, *Fast Ewald summation based on NFFT with
!>   mixed periodicity*, J. Comput. Phys. **285**, 280-315 (2015),
!>   DOI: [10.1016/j.jcp.2014.12.052](https://dx.doi.org/10.1016/j.jcp.2014.12.052).
!> This implementation evaluates their 1D and 2D Ewald formulations by direct
!> lattice summation.

!> @file multicharge/ewald.f90
!> Provides the unified 1D, 2D, and 3D Ewald summation interface

!> Stateless 1D, 2D, and 3D Ewald summation routines for periodic Gaussian charges
module multicharge_ewald
   use mctc_cutoff, only : get_periodic_axes
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : crossprod, matdet_3x3
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_ewald_cutoff, only : get_alpha, get_ewald_trans
   use multicharge_ewald_direct, only : get_amat_dir, get_damat_dir, get_hess_dir
   use multicharge_ewald_pbc1d, only : get_amat_rec_1d, get_damat_rec_1d, &
      & get_hess_rec_1d
   use multicharge_ewald_pbc2d, only : get_amat_rec_2d, get_damat_rec_2d, &
      & get_hess_rec_2d
   use multicharge_ewald_pbc3d, only : get_rec_weights, get_amat_rec_3d, &
      & get_damat_rec_3d, get_hess_rec_3d
   implicit none
   private

   public :: get_alpha, get_ewald_trans
   public :: get_ewald_amat, get_ewald_damat, add_ewald_hessian

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: zero(3) = 0.0_wp


contains


!> Evaluate the 1D, 2D, or 3D periodic Gaussian Ewald interaction matrix
subroutine get_ewald_amat(mol, wsc, alpha, dtrans, rtrans, rad, amat)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Direct-space Ewald lattice translations
   real(wp), intent(in) :: dtrans(:, :)

   !> Open half-space reciprocal Ewald vectors
   real(wp), intent(in) :: rtrans(:, :)

   !> Atomic Gaussian charge radii indexed by species ID
   real(wp), intent(in) :: rad(:)

   !> Periodic Ewald Coulomb matrix (filled on 1:nat, 1:nat)
   real(wp), intent(inout) :: amat(:, :)

   integer :: ndim, iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: vol, axis_vec(3), rself
   real(wp) :: rij(3), vec(3), gam, wsw, dtmp, rtmp, tmp
   real(wp), allocatable :: rec_weight(:)

   call get_ewald_cell(mol%periodic, mol%lattice, ndim, vol, axis_vec)

   select case (ndim)
   case (1)
      call get_amat_rec_1d(zero, alpha, vol, axis_vec, rtrans, rself)
   case (2)
      call get_amat_rec_2d(zero, alpha, vol, axis_vec, rtrans, rself)
   case (3)
      call get_rec_weights(vol, alpha, rtrans, weight=rec_weight)
      rself = sum(rec_weight)
   case default
      rself = 0.0_wp
   end select
   rself = rself - 2.0_wp * alpha / sqrtpi

   ! Every matrix element belongs to one pair; no thread-local matrix is needed.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(amat, mol, wsc, alpha, dtrans, rtrans, rad, rec_weight, &
   !$omp& ndim, vol, axis_vec, rself, npair) &
   !$omp& private(iat, jat, izp, jzp, img, gam, wsw, rij, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(rad(izp)**2 + rad(jzp)**2)
      if (iat == jat) then
         call get_amat_dir(zero, gam, alpha, dtrans, dtmp)
         tmp = dtmp + rself
      else
         rij = mol%xyz(:, jat) - mol%xyz(:, iat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         tmp = 0.0_wp
         do img = 1, wsc%nimg(jat, iat)
            vec = rij + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir(vec, gam, alpha, dtrans, dtmp)
            select case (ndim)
            case (1)
               call get_amat_rec_1d(vec, alpha, vol, axis_vec, rtrans, rtmp)
            case (2)
               call get_amat_rec_2d(vec, alpha, vol, axis_vec, rtrans, rtmp)
            case (3)
               call get_amat_rec_3d(vec, rtrans, rec_weight, rtmp)
            case default
               rtmp = 0.0_wp
            end select
            tmp = tmp + dtmp + rtmp
         end do
         tmp = tmp * wsw
      end if
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

end subroutine get_ewald_amat


!> Evaluate 1D, 2D, or 3D periodic Gaussian Ewald position and strain derivatives
subroutine get_ewald_damat(mol, wsc, alpha, dtrans, rtrans, rad, qvec, dadr, &
   & dadL, atrace)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Direct-space Ewald lattice translations
   real(wp), intent(in) :: dtrans(:, :)

   !> Open half-space reciprocal Ewald vectors
   real(wp), intent(in) :: rtrans(:, :)

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

   integer :: ndim, iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: vol, axis_vec(3), ds_self(3, 3)
   real(wp) :: gam, wsw, rij(3), vec(3)
   real(wp) :: dg(3), ds(3, 3), dgd(3), dsd(3, 3), dgr(3), dsr(3, 3)
   real(wp), allocatable :: rec_weight(:), rec_strain(:, :, :)

   call get_ewald_cell(mol%periodic, mol%lattice, ndim, vol, axis_vec)

   select case (ndim)
   case (1)
      call get_damat_rec_1d(zero, alpha, vol, axis_vec, rtrans, dgr, ds_self)
   case (2)
      call get_damat_rec_2d(zero, alpha, vol, axis_vec, rtrans, dgr, ds_self)
   case (3)
      call get_rec_weights(vol, alpha, rtrans, weight=rec_weight, strain=rec_strain)
      ds_self = sum(rec_strain, dim=3)
   case default
      ds_self = 0.0_wp
   end select

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   if (size(dadr, 3) > mol%nat) then
      dadr(:, :, mol%nat + 1:) = 0.0_wp
   end if

   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(mol, wsc, alpha, dtrans, rtrans, rad, qvec, dadr, &
   !$omp& rec_weight, rec_strain, ndim, vol, axis_vec, ds_self, npair) &
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
      if (iat == jat) then
         call get_damat_dir(zero, gam, alpha, dtrans, dgd, dsd)
         dadr(:, iat, iat) = 0.0_wp
         dadL(:, :, iat) = dadL(:, :, iat) + (dsd + ds_self) * qvec(iat)
      else
         rij = mol%xyz(:, jat) - mol%xyz(:, iat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         dg(:) = 0.0_wp
         ds(:, :) = 0.0_wp
         do img = 1, wsc%nimg(jat, iat)
            vec = rij + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_damat_dir(vec, gam, alpha, dtrans, dgd, dsd)
            select case (ndim)
            case (1)
               call get_damat_rec_1d(vec, alpha, vol, axis_vec, rtrans, dgr, dsr)
            case (2)
               call get_damat_rec_2d(vec, alpha, vol, axis_vec, rtrans, dgr, dsr)
            case (3)
               call get_damat_rec_3d(vec, rtrans, rec_weight, rec_strain, dgr, dsr)
            case default
               dgr = 0.0_wp
               dsr = 0.0_wp
            end select
            dg = dg + dgd + dgr
            ds = ds + dsd + dsr
         end do
         dg = dg * wsw
         ds = ds * wsw
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


!> Accumulate 1D, 2D, or 3D periodic Gaussian Ewald Cartesian Hessian blocks
subroutine add_ewald_hessian(mol, wsc, alpha, dtrans, rtrans, rad, qvec, adj, &
   & hessian, diagonal)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Direct-space Ewald lattice translations
   real(wp), intent(in) :: dtrans(:, :)

   !> Open half-space reciprocal Ewald vectors
   real(wp), intent(in) :: rtrans(:, :)

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

   integer :: ndim, iat, jat, izp, jzp, ii, jj, img, ipair, npair
   real(wp) :: vol, axis_vec(3)
   real(wp) :: rij(3), vec(3), box(3, 3), dbox(3, 3), rbox(3, 3), gam, weight, wsw
   real(wp), allocatable :: rec_hess(:, :, :)

   call get_ewald_cell(mol%periodic, mol%lattice, ndim, vol, axis_vec)
   if (ndim == 3) then
      call get_rec_weights(vol, alpha, rtrans, hess=rec_hess)
   end if

   ! Each pair owns its off-diagonal blocks. Only the diagonal blocks are
   ! reduced, avoiding a full (3N, 3N) Hessian copy for every OpenMP thread.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(mol, wsc, alpha, dtrans, rtrans, rad, qvec, adj, hessian, &
   !$omp& rec_hess, ndim, vol, axis_vec, npair) &
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
         select case (ndim)
         case (1)
            call get_hess_rec_1d(vec, alpha, vol, axis_vec, rtrans, rbox)
         case (2)
            call get_hess_rec_2d(vec, alpha, vol, axis_vec, rtrans, rbox)
         case (3)
            call get_hess_rec_3d(vec, rtrans, rec_hess, rbox)
         case default
            rbox = 0.0_wp
         end select
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


!> Extract periodic dimensionality, generalized cell measure, and wire/plane unit vector
pure subroutine get_ewald_cell(periodic, lattice, ndim, vol, axis_vec)

   !> Periodicity mask along each lattice axis
   logical, intent(in) :: periodic(:)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Number of periodic dimensions
   integer, intent(out) :: ndim

   !> Cell length (1D), area (2D), or volume (3D)
   real(wp), intent(out) :: vol

   !> Unit wire direction (1D) or slab normal (2D)
   real(wp), intent(out) :: axis_vec(3)

   integer :: p_axes(3)
   real(wp) :: n_vec(3)

   call get_periodic_axes(periodic, ndim, p_axes)
   vol = 0.0_wp
   axis_vec = 0.0_wp

   select case (ndim)
   case (1)
      vol = norm2(lattice(:, p_axes(1)))
      axis_vec = lattice(:, p_axes(1)) / vol
   case (2)
      n_vec = crossprod(lattice(:, p_axes(1)), lattice(:, p_axes(2)))
      vol = norm2(n_vec)
      axis_vec = n_vec / vol
   case (3)
      vol = abs(matdet_3x3(lattice))
   case default
      continue
   end select

end subroutine get_ewald_cell


end module multicharge_ewald
