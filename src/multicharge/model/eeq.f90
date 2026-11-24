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

!> Electronegativity equilibration charge model published in
!>
!> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth
!> and S. Grimme, *J. Chem. Phys.*, **2019**, 150, 154122.
!> DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
module multicharge_model_eeq
   use mctc_cutoff, only : get_lattice_points
   use mctc_env, only : error_type, fatal_error, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_ncoord, only : new_ncoord, cn_count
   use mctc_wignerseitz, only : wignerseitz_cell, new_wignerseitz_cell
   use multicharge_blas, only : gemm
   use multicharge_ewald, only : get_alpha, get_ewald_setup, get_amat_dir, &
      & get_damat_dir, get_hess_dir, get_amat_rec_3d, get_damat_rec_3d, &
      & get_hess_rec_3d
   use multicharge_model_cache, only : cache_container, model_cache
   use multicharge_model_hessian, only : coulomb_kernel, add_pair_hessian, &
      & add_pair_block
   use multicharge_model_type, only : mchrg_model_type
   implicit none
   private

   public :: eeq_model, new_eeq_model


   !> Cache for the electronegativity equilibration model
   type, extends(model_cache), public :: eeq_cache

      !> Ewald splitting parameter
      real(wp) :: alpha = 0.0_wp

   end type eeq_cache


   !> Electronegativity equilibration model
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
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)

   !> Regularization of the square root CN dependence of the electronegativity
   real(wp), parameter :: reg = 1.0e-14_wp


contains


!> Create a new electronegativity equilibration model
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


!> Store the coordination numbers and set up the periodic environment
subroutine update(self, mol, error, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Cache for the charge model
   type(cache_container), intent(inout) :: cache

   !> Coordination number
   real(wp), intent(in) :: cn(:)

   !> Local atomic partial charges (unused by EEQ)
   real(wp), intent(in), optional :: qloc(:)

   !> Derivative of the coordination number w.r.t. atomic positions
   real(wp), intent(in), optional :: dcndr(:, :, :)

   !> Derivative of the coordination number w.r.t. lattice strain
   real(wp), intent(in), optional :: dcndL(:, :, :)

   !> Derivative of the local charges w.r.t. atomic positions (unused by EEQ)
   real(wp), intent(in), optional :: dqlocdr(:, :, :)

   !> Derivative of the local charges w.r.t. lattice strain (unused by EEQ)
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
      call new_wignerseitz_cell(ptr%wsc, mol, error)
      if (allocated(error)) return

      call get_alpha(mol%lattice, ptr%alpha, error)
      if (allocated(error)) return
   end if

end subroutine update


!> Evaluate the right-hand side (electronegativity vector)
subroutine get_xvec(self, mol, cache, xvec)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache for the charge model
   type(cache_container), intent(inout) :: cache

   !> Electronegativity vector, total charge in the last entry
   real(wp), intent(out) :: xvec(:)

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


!> Evaluate the derivatives of the electronegativity vector
subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache for the charge model
   type(cache_container), intent(inout) :: cache

   !> Derivative of the EN vector w.r.t. atomic positions
   real(wp), intent(out), contiguous :: dxdr(:, :, :)

   !> Derivative of the EN vector w.r.t. lattice strain
   real(wp), intent(out), contiguous :: dxdL(:, :, :)

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
      tmp = self%kcnchi(izp) * (0.5_wp * ptr%cn(iat) + reg) &
         & / (ptr%cn(iat) + reg)**1.5_wp
      dxdr(:, :, iat) = tmp * ptr%dcndr(:, :, iat)
      dxdL(:, :, iat) = tmp * ptr%dcndL(:, :, iat)
   end do

end subroutine get_xvec_derivs


!> Evaluate the Coulomb matrix
subroutine get_coulomb_matrix(self, mol, cache, amat)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache for the charge model
   type(cache_container), intent(inout) :: cache

   !> Coulomb matrix, bordered by the charge constraint
   real(wp), intent(out) :: amat(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_amat_3d(self, mol, ptr%wsc, ptr%alpha, amat)
   else
      call get_amat_0d(self, mol, amat)
   end if

end subroutine get_coulomb_matrix


!> Coulomb matrix of Gaussian charges for a molecular structure
subroutine get_amat_0d(self, mol, amat)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Coulomb matrix, bordered by the charge constraint
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r2, gam, tmp

   ! Every matrix element belongs to one pair; no thread-local matrix is needed.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, self, npair) &
   !$omp private(iat, jat, izp, jzp, gam, vec, r2, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
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


!> Coulomb matrix of Gaussian charges for a 3D periodic structure (Ewald sum)
subroutine get_amat_3d(self, mol, wsc, alpha, amat)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Coulomb matrix, bordered by the charge constraint
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: rij(3), vec(3), gam, wsw, dtmp, rtmp, tmp
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_weight(:)

   call get_ewald_setup(mol, alpha, dtrans, rtrans, weight=rec_weight)

   ! Every matrix element belongs to one pair; no thread-local matrix is needed.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(amat, mol, self, wsc, alpha, dtrans, rtrans, rec_weight, npair) &
   !$omp private(iat, jat, izp, jzp, img, gam, wsw, rij, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
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
      if (iat == jat) then
         tmp = tmp - 2.0_wp * alpha / sqrtpi &
            & + (self%eta(izp) + sqrt2pi / self%rad(izp))
      end if
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_3d


!> Evaluate the derivatives of the Coulomb matrix contracted with a charge vector
subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache for the charge model
   type(cache_container), intent(inout) :: cache

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Derivative of the Coulomb matrix w.r.t. atomic positions
   real(wp), intent(out) :: dadr(:, :, :)

   !> Derivative of the Coulomb matrix w.r.t. lattice strain
   real(wp), intent(out) :: dadL(:, :, :)

   !> On-site trace of the position derivatives
   real(wp), intent(out) :: atrace(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_damat_3d(self, mol, ptr%wsc, ptr%alpha, qvec, dadr, dadL, atrace)
   else
      call get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   end if

end subroutine get_coulomb_derivs


!> Coulomb matrix derivatives for a molecular structure
subroutine get_damat_0d(self, mol, qvec, dadr, dadL, atrace)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Derivative of the Coulomb matrix w.r.t. atomic positions
   real(wp), intent(out) :: dadr(:, :, :)

   !> Derivative of the Coulomb matrix w.r.t. lattice strain
   real(wp), intent(out) :: dadL(:, :, :)

   !> On-site trace of the position derivatives
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r2, gam, arg, dtmp, dG(3), dS(3, 3)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   npair = mol%nat * (mol%nat + 1) / 2
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
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      vec = mol%xyz(:, jat) - mol%xyz(:, iat)
      r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      arg = gam * gam * r2
      dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) &
         & - erf(sqrt(arg)) / (r2 * sqrt(r2))
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


!> Coulomb matrix derivatives for a 3D periodic structure (Ewald sum)
subroutine get_damat_3d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Derivative of the Coulomb matrix w.r.t. atomic positions
   real(wp), intent(out) :: dadr(:, :, :)

   !> Derivative of the Coulomb matrix w.r.t. lattice strain
   real(wp), intent(out) :: dadL(:, :, :)

   !> On-site trace of the position derivatives
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: gam, wsw, rij(3), vec(3)
   real(wp) :: dg(3), ds(3, 3), dgd(3), dsd(3, 3), dgr(3), dsr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)
   real(wp), allocatable :: rec_weight(:), rec_strain(:, :, :)

   call get_ewald_setup(mol, alpha, dtrans, rtrans, weight=rec_weight, &
      & strain=rec_strain)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   npair = mol%nat * (mol%nat + 1) / 2
   ! Pair derivatives have unique owners. Only atom-wise traces and strain
   ! derivatives need reduction, with O(N) storage per thread.
   !$omp parallel default(none) &
   !$omp shared(dadr, mol, self, wsc, alpha, dtrans, rtrans, rec_weight, &
   !$omp& rec_strain, qvec, npair) &
   !$omp private(iat, jat, izp, jzp, ipair, img, gam, wsw, rij, vec, &
   !$omp& dg, ds, dgd, dsd, dgr, dsr) &
   !$omp reduction(+:atrace, dadL)
   !$omp do schedule(static)
   do iat = 1, size(dadr, 3)
      dadr(:, :, iat) = 0.0_wp
   end do
   !$omp end do
   !$omp do schedule(runtime)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
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
         ! Self-images only strain the cell, they do not move the atom
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
   !$omp end do
   !$omp end parallel

end subroutine get_damat_3d


!> u^T(x_ab - A_ab q), with u and q held fixed
subroutine add_model_hessian(self, mol, error, cache, adj, hessian)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Cache for the charge model
   type(cache_container), intent(inout) :: cache

   !> Adjoint vector
   real(wp), intent(in) :: adj(:)

   !> Cartesian Hessian matrix
   real(wp), intent(inout) :: hessian(:, :)

   integer :: iat, izp, nat, ii
   real(wp) :: c, weight
   real(wp), allocatable :: dEdcn(:), dcndr(:, :), work(:, :), trans(:, :)
   real(wp), allocatable :: diagonal(:, :, :)
   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)
   if (.not. associated(ptr)) then
      call fatal_error(error, "Incompatible EEQ response cache")
      return
   end if

   ! Coordination number contribution through the electronegativity vector
   nat = mol%nat
   allocate(dEdcn(nat), work(3 * nat, nat))
   dcndr = reshape(ptr%dcndr, [3 * nat, nat])
   do iat = 1, nat
      izp = mol%id(iat)
      c = ptr%cn(iat)
      dEdcn(iat) = adj(iat) * self%kcnchi(izp) * (0.5_wp * c + reg) / (c + reg)**1.5_wp
      weight = -adj(iat) * self%kcnchi(izp) * (0.25_wp * c + reg) / (c + reg)**2.5_wp
      work(:, iat) = weight * dcndr(:, iat)
   end do
   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans, &
      & error)
   if (allocated(error)) return

   call self%ncoord%add_coordination_number_hessian(mol, trans, dEdcn, hessian)
   call gemm(work, dcndr, hessian, transb="t", beta=1.0_wp)

   ! Coulomb contribution, diagonal blocks are accumulated separately
   allocate(diagonal(3, 3, nat), source=0.0_wp)
   if (any(mol%periodic)) then
      call add_coulomb_hessian_3d(self, mol, ptr%wsc, ptr%alpha, cache%qvec, adj, &
         & hessian, diagonal)
   else
      call add_coulomb_hessian_0d(self, mol, cache%qvec, adj, hessian, diagonal)
   end if

   !$omp parallel do default(none) schedule(static) &
   !$omp shared(nat, hessian, diagonal) private(ii)
   do iat = 1, nat
      ii = 3 * (iat - 1)
      hessian(ii + 1:ii + 3, ii + 1:ii + 3) = hessian(ii + 1:ii + 3, ii + 1:ii + 3) &
         & + diagonal(:, :, iat)
   end do
   !$omp end parallel do

end subroutine add_model_hessian


!> Coulomb Hessian blocks for a molecular structure
subroutine add_coulomb_hessian_0d(self, mol, qvec, adj, hessian, diagonal)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Adjoint vector
   real(wp), intent(in) :: adj(:)

   !> Off-diagonal Cartesian Hessian blocks
   real(wp), intent(inout) :: hessian(:, :)

   !> Diagonal Cartesian Hessian blocks
   real(wp), intent(inout) :: diagonal(:, :, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r, gam, val, dr, drr, weight

   ! Each pair owns its off-diagonal blocks. Only the diagonal blocks are
   ! reduced, avoiding a full (3N, 3N) Hessian copy for every OpenMP thread.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, qvec, adj, hessian, npair) &
   !$omp private(iat, jat, izp, jzp, vec, r, gam, val, dr, drr, weight) &
   !$omp reduction(+:diagonal)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      vec = mol%xyz(:, iat) - mol%xyz(:, jat)
      r = norm2(vec)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      call coulomb_kernel(r, gam, val, dr, drr)
      weight = -(adj(iat) * qvec(jat) + adj(jat) * qvec(iat))
      call add_pair_hessian(iat, jat, vec, weight * dr, weight * drr, hessian, &
         & diagonal)
   end do
   !$omp end parallel do

end subroutine add_coulomb_hessian_0d


!> Coulomb Hessian blocks for a 3D periodic structure (Ewald sum)
subroutine add_coulomb_hessian_3d(self, mol, wsc, alpha, qvec, adj, hessian, diagonal)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images
   type(wignerseitz_cell), intent(in) :: wsc

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Adjoint vector
   real(wp), intent(in) :: adj(:)

   !> Off-diagonal Cartesian Hessian blocks
   real(wp), intent(inout) :: hessian(:, :)

   !> Diagonal Cartesian Hessian blocks
   real(wp), intent(inout) :: diagonal(:, :, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: rij(3), vec(3), box(3, 3), dbox(3, 3), rbox(3, 3), gam, weight, wsw
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), rec_hess(:, :, :)

   call get_ewald_setup(mol, alpha, dtrans, rtrans, hess=rec_hess)

   ! Each pair owns its off-diagonal blocks. Only the diagonal blocks are
   ! reduced, avoiding a full (3N, 3N) Hessian copy for every OpenMP thread.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(self, mol, wsc, alpha, dtrans, rtrans, rec_hess, qvec, adj, &
   !$omp& hessian, npair) &
   !$omp private(iat, jat, izp, jzp, img, rij, vec, box, dbox, rbox, &
   !$omp& gam, weight, wsw) &
   !$omp reduction(+:diagonal)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      if (iat == jat) cycle
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
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
      call add_pair_block(iat, jat, (weight * wsw) * box, hessian, diagonal)
   end do
   !$omp end parallel do

end subroutine add_coulomb_hessian_3d


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
