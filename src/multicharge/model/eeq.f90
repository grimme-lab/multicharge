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
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_ncoord, only: new_ncoord, cn_count
   use mctc_wignerseitz, only: new_wignerseitz_cell
   use multicharge_blas, only: gemm
   use multicharge_ewald, only: get_alpha, get_ewald_trans, get_ewald_amat, &
      & get_ewald_damat, add_ewald_hessian
   use multicharge_model_cache, only: cache_container, model_cache
   use multicharge_model_hessian, only: coulomb_kernel, add_pair_hessian
   use multicharge_model_type, only: mchrg_model_type
   implicit none
   private

   public :: eeq_model, new_eeq_model

   !> Reusable data container for EEQ model evaluations
   type, extends(model_cache), public :: eeq_cache

      !> Direct-space Ewald lattice translations
      real(wp), allocatable :: dtrans(:, :)

      !> Open half-space reciprocal Ewald vectors
      real(wp), allocatable :: rtrans(:, :)

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

subroutine update(self, mol, error, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in), optional :: qloc(:)
   real(wp), intent(in), optional :: dcndr(:, :, :)
   real(wp), intent(in), optional :: dcndL(:, :, :)
   real(wp), intent(in), optional :: dqlocdr(:, :, :)
   real(wp), intent(in), optional :: dqlocdL(:, :, :)

   type(eeq_cache), pointer :: ptr
   real(wp) :: gam_min

   call taint(cache, ptr)

   ! Refer CN arrays in cache
   ptr%cn = cn
   if (present(dcndr) .and. present(dcndL)) then
      ptr%dcndr = dcndr
      ptr%dcndL = dcndL
   end if

   if (any(mol%periodic)) then
      if (.not. all(mol%periodic) .and. abs(mol%charge) > eps) then
         call fatal_error(error, &
            & "Charged 1D and 2D periodic systems are not supported by EEQ")
         return
      end if

      ! Create WSC
      call new_wignerseitz_cell(ptr%wsc, mol, error)
      if (allocated(error)) return

      call get_alpha(mol%lattice, ptr%alpha, error, mol%periodic)
      if (allocated(error)) return

      gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
      call get_ewald_trans(mol%periodic, mol%lattice, ptr%alpha, gam_min, &
         & ptr%dtrans, ptr%rtrans, error)
      if (allocated(error)) return
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
   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans, error)
   if (allocated(error)) return
   call self%ncoord%add_coordination_number_hessian(mol, trans, dEdcn, hessian)
   call gemm(work, dcndr, hessian, transb="t", beta=1.0_wp)

   allocate(diagonal(3, 3, nat), source=0.0_wp)
   if (any(mol%periodic)) then
      call add_ewald_hessian(mol, ptr%wsc, ptr%alpha, ptr%dtrans, ptr%rtrans, &
         & self%rad, cache%qvec, adj, hessian, diagonal)
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

subroutine get_coulomb_matrix(self, mol, cache, amat)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(out) :: amat(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_amat_pbc(self, mol, ptr, amat)
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


!> Evaluate the 1D, 2D, or 3D Ewald Coulomb matrix
subroutine get_amat_pbc(self, mol, ptr, amat)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> EEQ model cache
   type(eeq_cache), intent(in) :: ptr

   !> Coulomb matrix
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, izp

   call get_ewald_amat(mol, ptr%wsc, ptr%alpha, ptr%dtrans, ptr%rtrans, &
      & self%rad, amat)

   do iat = 1, mol%nat
      izp = mol%id(iat)
      amat(iat, iat) = amat(iat, iat) + self%eta(izp) + sqrt2pi / self%rad(izp)
   end do

   if (size(amat, 1) == mol%nat + 1) then
      amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
      amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
      amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
   end if

end subroutine get_amat_pbc

subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
   class(eeq_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(cache_container), intent(inout) :: cache
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_ewald_damat(mol, ptr%wsc, ptr%alpha, ptr%dtrans, ptr%rtrans, &
         & self%rad, qvec, dadr, dadL, atrace)
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
