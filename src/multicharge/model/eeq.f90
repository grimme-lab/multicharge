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

!> @file multicharge/model/eeq.f90
!> Provides implementation of the electronegativity equilibration model (EEQ)

!> Electronegativity equlibration charge model published in
!>
!> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth
!> and S. Grimme, *J. Chem. Phys.*, **2019**, 150, 154122.
!> DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
module multicharge_model_eeq
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matdet_3x3
   use mctc_ncoord, only: new_ncoord, cn_count
   use multicharge_blas, only: gemm
   use multicharge_ewald, only: get_alpha
   use multicharge_model_cache, only: cache_container, model_cache
   use multicharge_model_hessian, only: coulomb_kernel, radial_hessian, add_pair_hessian, add_pair_block
   use multicharge_model_type, only: mchrg_model_type, get_dir_trans, get_rec_trans
   use multicharge_wignerseitz, only: wignerseitz_cell_type, new_wignerseitz_cell
   implicit none
   private

   public :: eeq_model, new_eeq_model

   type, extends(model_cache), public :: eeq_cache
   end type eeq_cache

   type, extends(mchrg_model_type) :: eeq_model
   contains
      !> Update and allocate cache
      procedure :: update
      !> Calculate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Calculate derivatives of Coulomb matrix
      procedure :: get_coulomb_derivs
      !> Calculate right-hand side (electronegativity vector)
      procedure :: get_xvec
      !> Calculate EN vector derivatives
      procedure :: get_xvec_derivs
      !> Cartesian second-response curvature at fixed lattice
      procedure :: add_model_hessian
   end type eeq_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains

subroutine new_eeq_model(self, mol, error, chi, rad, eta, kcnchi, &
   & cutoff, cn_exp, rcov, cn_max)
   !> Electronegativity equilibration model
   type(eeq_model), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Error handling
   type(error_type), allocatable, intent(out) :: error
   !> Electronegativity
   real(wp), intent(in) :: chi(:)
   !> Exponent gaussian charge
   real(wp), intent(in) :: rad(:)
   !> Chemical hardness
   real(wp), intent(in) :: eta(:)
   !> CN scaling factor for electronegativity
   real(wp), intent(in) :: kcnchi(:)
   !> Cutoff radius for coordination number
   real(wp), intent(in), optional :: cutoff
   !> Steepness of the CN counting function
   real(wp), intent(in), optional :: cn_exp
   !> Covalent radii for CN
   real(wp), intent(in), optional :: rcov(:)
   !> Maximum CN cutoff for CN
   real(wp), intent(in), optional :: cn_max

   self%chi = chi
   self%rad = rad
   self%eta = eta
   self%kcnchi = kcnchi

   call new_ncoord(self%ncoord, mol, cn_count%erf, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, cut=cn_max)

end subroutine new_eeq_model

subroutine update(self, mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in), optional :: qloc(:)
   real(wp), intent(in), optional :: dcndr(:, :, :)
   real(wp), intent(in), optional :: dcndL(:, :, :)
   real(wp), intent(in), optional :: dqlocdr(:, :, :)
   real(wp), intent(in), optional :: dqlocdL(:, :, :)

   type(eeq_cache), pointer :: ptr

   call taint(cache, ptr)

   ! Refer CN arrays in cache
   ptr%cn = cn
   if (present(dcndr) .and. present(dcndL)) then
      ptr%dcndr = dcndr
      ptr%dcndL = dcndL
   end if

   if (any(mol%periodic)) then
      ! Create WSC
      call new_wignerseitz_cell(ptr%wsc, mol)
      call get_alpha(mol%lattice, ptr%alpha)
   end if

end subroutine update

subroutine get_xvec(self, mol, cache, xvec)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out) :: xvec(:)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, self, xvec, ptr) private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) / sqrt(ptr%cn(iat) + reg)
      xvec(iat) = -self%chi(izp) + tmp * ptr%cn(iat)
   end do
   xvec(mol%nat + 1) = mol%charge

end subroutine get_xvec

subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out), contiguous :: dxdr(:, :, :)
   real(wp), intent(out), contiguous :: dxdL(:, :, :)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   dxdr(:, :, :) = 0.0_wp
   dxdL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, self, ptr, dxdr, dxdL) &
   !$omp private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) * (0.5_wp*ptr%cn(iat) + reg) / (ptr%cn(iat) + reg)**1.5_wp
      dxdr(:, :, iat) = tmp * ptr%dcndr(:, :, iat)
      dxdL(:, :, iat) = tmp * ptr%dcndL(:, :, iat)
   end do

end subroutine get_xvec_derivs

!> u^T(x_ab-A_ab q), with u and q held fixed.
subroutine add_model_hessian(self, mol, error, cache, adj, hessian)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: adj(:)
   real(wp), intent(inout) :: hessian(:, :)

   type(eeq_cache), pointer :: ptr
   integer :: iat, jat, izp, jzp, nat, ipair, npair, ii
   real(wp), parameter :: reg = 1.0e-14_wp
   real(wp) :: c, vec(3), r, gam, val, dr, drr, weight
   real(wp), allocatable :: dEdcn(:), dcndr(:, :), work(:, :), trans(:, :)
   real(wp), allocatable :: diagonal(:, :, :)

   call view(cache, ptr)
   if (.not. associated(ptr)) then
      call fatal_error(error, "Incompatible EEQ response cache")
      return
   end if
   nat = mol%nat
   allocate(dEdcn(nat), work(3*nat, nat))
   dcndr = reshape(ptr%dcndr, [3*nat, nat])
   do iat = 1, nat
      izp = mol%id(iat)
      c = ptr%cn(iat)
      dEdcn(iat) = adj(iat)*self%kcnchi(izp)*(0.5_wp*c+reg)/(c+reg)**1.5_wp
      weight = -adj(iat)*self%kcnchi(izp)*(0.25_wp*c+reg)/(c+reg)**2.5_wp
      work(:, iat) = weight*dcndr(:, iat)
   end do
   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans)
   call self%ncoord%add_coordination_number_hessian(mol, trans, dEdcn, hessian)
   call gemm(work, dcndr, hessian, transb="t", beta=1.0_wp)

   allocate(diagonal(3, 3, nat), source=0.0_wp)
   if (any(mol%periodic)) then
      call add_periodic_hessian(self, mol, ptr%wsc, ptr%alpha, cache%qvec, adj, hessian, diagonal)
   else
      npair = nat*(nat+1)/2
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(self, mol, cache, adj, hessian, npair) &
      !$omp private(iat, jat, izp, jzp, vec, r, gam, val, dr, drr, weight) &
      !$omp reduction(+:diagonal)
      do ipair = 1, npair
         iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
         if (iat*(iat+1)/2 < ipair) iat = iat+1
         jat = ipair-iat*(iat-1)/2
         if (iat == jat) cycle
         izp = mol%id(iat)
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat)-mol%xyz(:, jat)
         r = norm2(vec)
         gam = 1.0_wp/sqrt(self%rad(izp)**2+self%rad(jzp)**2)
         call coulomb_kernel(r, gam, val, dr, drr)
         weight = -(adj(iat)*cache%qvec(jat)+adj(jat)*cache%qvec(iat))
         call add_pair_hessian(iat, jat, vec, weight*dr, weight*drr, hessian, diagonal)
      end do
      !$omp end parallel do
   end if

   !$omp parallel do default(none) schedule(static) shared(nat, hessian, diagonal) private(ii)
   do iat = 1, nat
      ii = 3*(iat-1)
      hessian(ii+1:ii+3, ii+1:ii+3) = hessian(ii+1:ii+3, ii+1:ii+3)+diagonal(:, :, iat)
   end do
   !$omp end parallel do
end subroutine add_model_hessian

!> Curvature of the same finite Ewald sums and WS average used for A.
!> The lattice, Ewald parameter, and self-image distances are fixed here.
subroutine add_periodic_hessian(self, mol, wsc, alpha, qvec, adj, hessian, diagonal)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha, qvec(:), adj(:)
   real(wp), intent(inout) :: hessian(:, :)
   real(wp), intent(inout) :: diagonal(:, :, :)
   integer :: iat, jat, izp, jzp, img, itr, ipair, npair
   real(wp) :: rij(3), vec(3), box(3, 3), gam, r, weight, wsw
   real(wp) :: val, dr, drr, ar, arr
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_weight(:)

   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)
   call get_rec_weights(abs(matdet_3x3(mol%lattice)), alpha, rtrans, rec_weight)

   ! Each pair owns its off-diagonal blocks. Only the diagonal blocks are
   ! reduced, avoiding a full (3N, 3N) Hessian copy for every OpenMP thread.
   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, wsc, alpha, qvec, adj, hessian, dtrans, rtrans, rec_weight, npair) &
   !$omp private(iat, jat, izp, jzp, img, itr, rij, vec, box, gam, r, weight, wsw, val, dr, drr, ar, arr) &
   !$omp reduction(+:diagonal)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp/sqrt(self%rad(izp)**2+self%rad(jzp)**2)
      weight = -(adj(iat)*qvec(jat)+adj(jat)*qvec(iat))
      wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
      box = 0.0_wp
      do img = 1, wsc%nimg(jat, iat)
         rij = mol%xyz(:, iat)-mol%xyz(:, jat)-wsc%trans(:, wsc%tridx(img, jat, iat))
         do itr = 1, size(dtrans, 2)
            vec = rij+dtrans(:, itr)
            r = norm2(vec)
            if (r < eps) cycle
            call coulomb_kernel(r, gam, val, dr, drr)
            call coulomb_kernel(r, alpha, val, ar, arr)
            box = box+wsw*radial_hessian(vec, dr-ar, drr-arr)
         end do
         do itr = 1, size(rtrans, 2)
            vec = rtrans(:, itr)
            box = box-wsw*rec_weight(itr)*cos(dot_product(rij, vec)) &
               & *spread(vec, 1, 3)*spread(vec, 2, 3)
         end do
      end do
      call add_pair_block(iat, jat, weight*box, hessian, diagonal)
   end do
   !$omp end parallel do
end subroutine add_periodic_hessian


subroutine get_coulomb_matrix(self, mol, cache, amat)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out) :: amat(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_amat_3d(self, mol, ptr%wsc, ptr%alpha, amat)
   else
      call get_amat_0d(self, mol, amat)
   end if
end subroutine get_coulomb_matrix

subroutine get_amat_0d(self, mol, amat)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r2, gam, tmp

   ! Every matrix element belongs to one pair; no thread-local matrix is needed.
   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, self, npair) &
   !$omp private(iat, jat, izp, jzp, gam, vec, r2, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      if (iat == jat) then
         tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      else
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2 * gam)) / sqrt(r2)
      end if
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
end subroutine get_amat_0d

subroutine get_amat_3d(self, mol, wsc, alpha, amat)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: vec(3), gam, wsw, dtmp, rtmp, tmp
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_weight(:)

   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)
   call get_rec_weights(abs(matdet_3x3(mol%lattice)), alpha, rtrans, rec_weight)

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, self, wsc, dtrans, rtrans, rec_weight, alpha, npair) &
   !$omp private(iat, jat, izp, jzp, img, gam, wsw, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
      tmp = 0.0_wp
      do img = 1, wsc%nimg(jat, iat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
         call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
         call get_amat_rec_3d(vec, rtrans, rec_weight, rtmp)
         tmp = tmp + (dtmp + rtmp) * wsw
      end do
      if (iat == jat) tmp = tmp + self%eta(izp) + sqrt2pi / self%rad(izp) - 2 * alpha / sqrtpi
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
end subroutine get_amat_3d

subroutine get_amat_dir_3d(rij, gam, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      tmp = erf(gam * r1) / r1 - erf(alp * r1) / r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d

subroutine get_amat_rec_3d(rij, trans, weight, amat)
   real(wp), intent(in) :: rij(3), trans(:, :), weight(:)
   real(wp), intent(out) :: amat
   integer :: itr

   amat = 0.0_wp
   do itr = 1, size(trans, 2)
      amat = amat + weight(itr) * cos(dot_product(rij, trans(:, itr)))
   end do
end subroutine get_amat_rec_3d

subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_damat_3d(self, mol, ptr%wsc, ptr%alpha, qvec, dadr, dadL, atrace)
   else
      call get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   end if
end subroutine get_coulomb_derivs

subroutine get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r2, gam, arg, dtmp, dG(3), dS(3, 3)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   npair = mol%nat*(mol%nat+1)/2
   ! Pair derivatives have unique owners. Only atom-wise traces and strain
   ! derivatives need reduction, with O(N) storage per thread.
   !$omp parallel default(none) shared(dadr, mol, self, qvec, npair) &
   !$omp private(iat, jat, izp, jzp, ipair, gam, r2, vec, dG, dS, dtmp, arg) &
   !$omp reduction(+:atrace, dadL)
   !$omp do schedule(static)
   do iat = 1, size(dadr, 3)
      dadr(:, :, iat) = 0.0_wp
   end do
   !$omp end do
   !$omp do schedule(runtime)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      vec = mol%xyz(:, jat) - mol%xyz(:, iat)
      r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      arg = gam * gam * r2
      dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) - erf(sqrt(arg)) / (r2 * sqrt(r2))
      dG = dtmp * vec
      dS = spread(dG, 1, 3) * spread(vec, 2, 3)
      atrace(:, iat) = atrace(:, iat) - dG * qvec(jat)
      atrace(:, jat) = atrace(:, jat) + dG * qvec(iat)
      dadr(:, iat, jat) = -dG * qvec(iat)
      dadr(:, jat, iat) = +dG * qvec(jat)
      dadL(:, :, jat) = dadL(:, :, jat) + dS * qvec(iat)
      dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(jat)
   end do
   !$omp end do
   !$omp end parallel
end subroutine get_damat_0d

subroutine get_damat_3d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha, qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: gam, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_weight(:), rec_strain(:, :, :)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)
   call get_rec_weights(abs(matdet_3x3(mol%lattice)), alpha, rtrans, rec_weight, rec_strain)

   npair = mol%nat*(mol%nat+1)/2
   !$omp parallel default(none) &
   !$omp shared(mol, self, wsc, alpha, dtrans, rtrans, rec_weight, rec_strain, qvec, dadr, npair) &
   !$omp private(iat, jat, izp, jzp, img, ipair, gam, wsw, vec, dG, dS, dGr, dSr, dGd, dSd) &
   !$omp reduction(+:atrace, dadL)
   !$omp do schedule(static)
   do iat = 1, size(dadr, 3)
      dadr(:, :, iat) = 0.0_wp
   end do
   !$omp end do
   !$omp do schedule(runtime)
   do ipair = 1, npair
      iat = int(0.5_wp*(sqrt(8.0_wp*real(ipair, wp)+1.0_wp)-1.0_wp))
      if (iat*(iat+1)/2 < ipair) iat = iat+1
      jat = ipair-iat*(iat-1)/2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      dG(:) = 0.0_wp
      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
      do img = 1, wsc%nimg(jat, iat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
         call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
         call get_damat_rec_3d(vec, rtrans, rec_weight, rec_strain, dGr, dSr)
         dG = dG + (dGd + dGr) * wsw
         dS = dS + (dSd + dSr) * wsw
      end do
      if (iat == jat) then
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(iat)
      else
         atrace(:, iat) = atrace(:, iat) - dG * qvec(jat)
         atrace(:, jat) = atrace(:, jat) + dG * qvec(iat)
         dadr(:, iat, jat) = -dG * qvec(iat)
         dadr(:, jat, iat) = +dG * qvec(jat)
         dadL(:, :, jat) = dadL(:, :, jat) + dS * qvec(iat)
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(jat)
      end if
   end do
   !$omp end do
   !$omp end parallel
end subroutine get_damat_3d

subroutine get_damat_dir_3d(rij, gam, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, atmp, gam2, alp2

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   gam2 = gam * gam
   alp2 = alp * alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1 * r1
      gtmp = +2 * gam * exp(-r2 * gam2) / (sqrtpi * r2) - erf(r1 * gam) / (r2 * r1)
      atmp = -2 * alp * exp(-r2 * alp2) / (sqrtpi * r2) + erf(r1 * alp) / (r2 * r1)
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * spread(vec, 1, 3) * spread(vec, 2, 3)
   end do

end subroutine get_damat_dir_3d

subroutine get_damat_rec_3d(rij, trans, weight, strain, dg, ds)
   real(wp), intent(in) :: rij(3), trans(:, :), weight(:), strain(:, :, :)
   real(wp), intent(out) :: dg(3), ds(3, 3)
   integer :: itr
   real(wp) :: gv

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   do itr = 1, size(trans, 2)
      gv = dot_product(rij, trans(:, itr))
      dg(:) = dg - weight(itr) * sin(gv) * trans(:, itr)
      ds(:, :) = ds + cos(gv) * strain(:, :, itr)
   end do
end subroutine get_damat_rec_3d

!> Reciprocal-space factors independent of the atom pair.
!> Prepare once for each matrix/response evaluation, including strain when needed.
subroutine get_rec_weights(vol, alpha, trans, weight, strain)
   real(wp), intent(in) :: vol, alpha, trans(:, :)
   real(wp), allocatable, intent(out) :: weight(:)
   real(wp), allocatable, intent(out), optional :: strain(:, :, :)
   integer :: itr
   real(wp) :: g2, fac, vec(3)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])

   fac = 4 * pi / vol
   allocate(weight(size(trans, 2)), source=0.0_wp)
   if (present(strain)) allocate(strain(3, 3, size(trans, 2)), source=0.0_wp)
   do itr = 1, size(trans, 2)
      vec = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      weight(itr) = fac * exp(-0.25_wp * g2 / alpha**2) / g2
      if (present(strain)) then
         strain(:, :, itr) = weight(itr) * ((2.0_wp / g2 + 0.5_wp / alpha**2) &
            & * spread(vec, 1, 3) * spread(vec, 2, 3) - unity)
      end if
   end do
end subroutine get_rec_weights

!> Inspect cache and reallocate it in case of type mismatch
subroutine taint(cache, ptr)
   !> Instance of the cache
   type(cache_container), target, intent(inout) :: cache
   !> Reference to the cache
   type(eeq_cache), pointer, intent(out) :: ptr

   if (allocated(cache%raw)) then
      call view(cache, ptr)
      if (associated(ptr)) return
      deallocate(cache%raw)
   end if

   if (.not. allocated(cache%raw)) then
      block
         type(eeq_cache), allocatable :: tmp
         allocate(tmp)
         call move_alloc(tmp, cache%raw)
      end block
   end if

   call view(cache, ptr)
end subroutine taint

!> Return reference to cache after resolving its type
subroutine view(cache, ptr)
   !> Instance of the cache
   type(cache_container), target, intent(inout) :: cache
   !> Reference to the cache
   type(eeq_cache), pointer, intent(out) :: ptr
   nullify(ptr)
   select type(target => cache%raw)
   type is(eeq_cache)
      ptr => target
   end select
end subroutine view

end module multicharge_model_eeq
