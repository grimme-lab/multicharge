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
   use mctc_io_math, only: matdet_3x3, matinv_3x3, outerprod
   use mctc_ncoord, only: new_ncoord, cn_count
   use mctc_wignerseitz, only: new_wignerseitz_cell
   use multicharge_blas, only: gemm
   use multicharge_ewald, only: get_alpha
   use multicharge_model_cache, only: cache_container, model_cache
   use multicharge_model_hessian, only: coulomb_kernel, radial_hessian, &
      & add_pair_hessian, add_pair_block
   use multicharge_model_type, only: mchrg_model_type
   implicit none
   private

   public :: eeq_model, new_eeq_model

   !> Reusable data container for EEQ model evaluations
   type, extends(model_cache), public :: eeq_cache

      !> Cell volume
      real(wp) :: vol = 0.0_wp

      !> Direct-space Ewald translations
      real(wp), allocatable :: dtrans(:, :)

      !> Nonzero reciprocal-space Ewald vectors
      real(wp), allocatable :: rtrans(:, :)

   end type eeq_cache

   !> Electronegativity equilibration charge model
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

   real(wp), parameter :: twopi = 2.0_wp * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: unity(3, 3) = reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])


contains


!> Initialize a new EEQ model instance
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


!> Update and populate the EEQ model cache
subroutine update(self, mol, error, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Cache container
   type(cache_container), intent(inout) :: cache

   !> Coordination number
   real(wp), intent(in) :: cn(:)

   !> Local atomic partial charges
   real(wp), intent(in), optional :: qloc(:)

   !> Derivative of coordination number w.r.t. atomic positions
   real(wp), intent(in), optional :: dcndr(:, :, :)

   !> Derivative of coordination number w.r.t. lattice vectors
   real(wp), intent(in), optional :: dcndL(:, :, :)

   !> Derivative of local charges w.r.t. atomic positions
   real(wp), intent(in), optional :: dqlocdr(:, :, :)

   !> Derivative of local charges w.r.t. lattice vectors
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

      call setup_ewald(mol, ptr, error)
      if (allocated(error)) return
   end if

end subroutine update


!> Set up Ewald splitting parameter, cell geometry, and lattice translations
subroutine setup_ewald(mol, cache, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> EEQ model cache
   type(eeq_cache), intent(inout) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: rep(3) = [2, 2, 2]
   real(wp) :: rec_lat(3, 3)

   call get_alpha(mol%lattice, cache%alpha, error)
   if (allocated(error)) return

   cache%vol = abs(matdet_3x3(mol%lattice))
   call get_lattice_points(mol%lattice, rep, .true., cache%dtrans, error)
   if (allocated(error)) return

   rec_lat = twopi * transpose(matinv_3x3(mol%lattice))
   call get_lattice_points(rec_lat, rep, .false., cache%rtrans, error)

end subroutine setup_ewald


!> Evaluate the electronegativity right-hand side vector
subroutine get_xvec(self, mol, cache, xvec)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache container
   type(cache_container), intent(inout) :: cache

   !> Right-hand side vector
   real(wp), intent(out) :: xvec(:)

   real(wp), parameter :: reg = 1.0e-14_wp
   integer :: iat, izp
   real(wp) :: tmp
   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(mol, self, xvec, ptr) private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) / sqrt(ptr%cn(iat) + reg)
      xvec(iat) = -self%chi(izp) + tmp * ptr%cn(iat)
   end do
   xvec(mol%nat + 1) = mol%charge

end subroutine get_xvec


!> Evaluate derivatives of the electronegativity right-hand side vector
subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache container
   type(cache_container), intent(inout) :: cache

   !> Derivative of right-hand side vector w.r.t. atomic positions
   real(wp), intent(out), contiguous :: dxdr(:, :, :)

   !> Derivative of right-hand side vector w.r.t. lattice vectors
   real(wp), intent(out), contiguous :: dxdL(:, :, :)

   real(wp), parameter :: reg = 1.0e-14_wp
   integer :: iat, izp
   real(wp) :: tmp
   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   dxdr(:, :, :) = 0.0_wp
   dxdL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(mol, self, ptr, dxdr, dxdL) private(iat, izp, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      tmp = self%kcnchi(izp) * (0.5_wp * ptr%cn(iat) + reg) &
         & / (ptr%cn(iat) + reg)**1.5_wp
      dxdr(:, :, iat) = tmp * ptr%dcndr(:, :, iat)
      dxdL(:, :, iat) = tmp * ptr%dcndL(:, :, iat)
   end do

end subroutine get_xvec_derivs


!> Contract model-specific second-response curvature u^T(x_ab - A_ab q) at fixed u and q
subroutine add_model_hessian(self, mol, error, cache, adj, hessian)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Cache container
   type(cache_container), intent(inout) :: cache

   !> Adjoint vector
   real(wp), intent(in) :: adj(:)

   !> Cartesian Hessian matrix
   real(wp), intent(inout) :: hessian(:, :)

   type(eeq_cache), pointer :: ptr
   integer :: iat, jat, izp, jzp, nat, ipair, npair, ii
   real(wp), parameter :: reg = 1.0e-14_wp
   real(wp) :: c, vec(3), r, gam, val, dr, drr, weight
   real(wp), allocatable :: dedcn(:), dcndr(:, :), work(:, :), trans(:, :)
   real(wp), allocatable :: diagonal(:, :, :)

   call view(cache, ptr)
   if (.not. associated(ptr)) then
      call fatal_error(error, "Incompatible EEQ response cache")
      return
   end if

   nat = mol%nat
   allocate(dedcn(nat), work(3 * nat, nat))
   dcndr = reshape(ptr%dcndr, [3 * nat, nat])
   do iat = 1, nat
      izp = mol%id(iat)
      c = ptr%cn(iat)
      dedcn(iat) = adj(iat) * self%kcnchi(izp) * (0.5_wp * c + reg) &
         & / (c + reg)**1.5_wp
      weight = -adj(iat) * self%kcnchi(izp) * (0.25_wp * c + reg) &
         & / (c + reg)**2.5_wp
      work(:, iat) = weight * dcndr(:, iat)
   end do
   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, &
      & trans, error)
   if (allocated(error)) return

   call self%ncoord%add_coordination_number_hessian(mol, trans, dedcn, hessian)
   call gemm(work, dcndr, hessian, transb="t", beta=1.0_wp)

   allocate(diagonal(3, 3, nat), source=0.0_wp)
   if (any(mol%periodic)) then
      call add_periodic_hessian(self, mol, ptr, cache%qvec, adj, hessian, diagonal)
   else
      npair = nat * (nat + 1) / 2
      !$omp parallel do default(none) schedule(runtime) &
      !$omp& shared(self, mol, cache, adj, hessian, npair) &
      !$omp& private(iat, jat, izp, jzp, vec, r, gam, val, dr, drr, weight) &
      !$omp& reduction(+:diagonal)
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
         weight = -(adj(iat) * cache%qvec(jat) + adj(jat) * cache%qvec(iat))
         call add_pair_hessian(iat, jat, vec, weight * dr, weight * drr, &
            & hessian, diagonal)
      end do
      !$omp end parallel do
   end if

   !$omp parallel do default(none) schedule(static) &
   !$omp& shared(nat, hessian, diagonal) private(ii)
   do iat = 1, nat
      ii = 3 * (iat - 1)
      hessian(ii + 1:ii + 3, ii + 1:ii + 3) = &
         & hessian(ii + 1:ii + 3, ii + 1:ii + 3) + diagonal(:, :, iat)
   end do
   !$omp end parallel do

end subroutine add_model_hessian


!> Curvature of the finite Ewald sums and Wigner-Seitz average used for A.
!> The lattice, Ewald parameter, and self-image distances are fixed here.
subroutine add_periodic_hessian(self, mol, ptr, qvec, adj, hessian, diagonal)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> EEQ model cache
   type(eeq_cache), intent(in) :: ptr

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Adjoint vector
   real(wp), intent(in) :: adj(:)

   !> Off-diagonal Cartesian Hessian blocks
   real(wp), intent(inout) :: hessian(:, :)

   !> Diagonal Cartesian Hessian blocks
   real(wp), intent(inout) :: diagonal(:, :, :)

   integer :: iat, jat, izp, jzp, img, itr, ipair, npair
   real(wp) :: rij(3), vec(3), box(3, 3), rbox(3, 3), gam, r, weight, wsw
   real(wp) :: val, dr, drr, ar, arr
   real(wp), allocatable :: rec_weight(:)

   call get_rec_weights(ptr%vol, ptr%alpha, ptr%rtrans, rec_weight)

   ! Each pair owns its off-diagonal blocks. Only the diagonal blocks are
   ! reduced, avoiding a full (3N, 3N) Hessian copy for every OpenMP thread.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(self, mol, ptr, qvec, adj, hessian, rec_weight, npair) &
   !$omp& private(iat, jat, izp, jzp, img, itr, rij, vec, box, rbox, gam, r, &
   !$omp& weight, wsw, val, dr, drr, ar, arr) &
   !$omp& reduction(+:diagonal)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      if (iat == jat) cycle

      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      weight = -(adj(iat) * qvec(jat) + adj(jat) * qvec(iat))
      wsw = 1.0_wp / real(ptr%wsc%nimg(jat, iat), wp)
      box = 0.0_wp
      do img = 1, ptr%wsc%nimg(jat, iat)
         rij = mol%xyz(:, iat) - mol%xyz(:, jat) &
            & - ptr%wsc%trans(:, ptr%wsc%tridx(img, jat, iat))
         do itr = 1, size(ptr%dtrans, 2)
            vec = rij + ptr%dtrans(:, itr)
            r = norm2(vec)
            call coulomb_kernel(r, gam, val, dr, drr)
            call coulomb_kernel(r, ptr%alpha, val, ar, arr)
            box = box + wsw * radial_hessian(vec, dr - ar, drr - arr)
         end do
         call get_hess_rec_3d(rij, ptr%rtrans, rec_weight, rbox)
         box = box + wsw * rbox
      end do
      call add_pair_block(iat, jat, weight * box, hessian, diagonal)
   end do
   !$omp end parallel do

end subroutine add_periodic_hessian


!> Reciprocal-space pair Hessian block for 3D Ewald
pure subroutine get_hess_rec_3d(rij, trans, weight, rbox)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space weights
   real(wp), intent(in) :: weight(:)

   !> Reciprocal-space Hessian block
   real(wp), intent(out) :: rbox(3, 3)

   integer :: itr
   real(wp) :: vec(3)

   rbox = 0.0_wp
   do itr = 1, size(trans, 2)
      vec = trans(:, itr)
      rbox = rbox - weight(itr) * cos(dot_product(rij, vec)) &
         & * outerprod(vec, vec)
   end do

end subroutine get_hess_rec_3d


!> Evaluate the EEQ Coulomb interaction matrix
subroutine get_coulomb_matrix(self, mol, cache, amat)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache container
   type(cache_container), intent(inout) :: cache

   !> Coulomb matrix
   real(wp), intent(out) :: amat(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_amat_pbc(self, mol, ptr, amat)
   else
      call get_amat_0d(self, mol, amat)
   end if

end subroutine get_coulomb_matrix


!> Evaluate the molecular (0D) Coulomb matrix
subroutine get_amat_0d(self, mol, amat)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Coulomb matrix
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r2, gam, tmp

   ! Every matrix element belongs to one pair; no thread-local matrix is needed.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(amat, mol, self, npair) &
   !$omp& private(iat, jat, izp, jzp, gam, vec, r2, tmp)
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


!> Evaluate the periodic Ewald Coulomb matrix
subroutine get_amat_pbc(self, mol, ptr, amat)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> EEQ model cache
   type(eeq_cache), intent(in) :: ptr

   !> Coulomb matrix
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: vec(3), gam, wsw, dtmp, rtmp, tmp
   real(wp), allocatable :: rec_weight(:)

   amat = 0.0_wp
   call get_rec_weights(ptr%vol, ptr%alpha, ptr%rtrans, rec_weight)

   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(amat, mol, self, ptr, rec_weight, npair) &
   !$omp& private(iat, jat, izp, jzp, img, gam, wsw, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      wsw = 1.0_wp / real(ptr%wsc%nimg(jat, iat), wp)
      tmp = 0.0_wp
      do img = 1, ptr%wsc%nimg(jat, iat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat) &
            & + ptr%wsc%trans(:, ptr%wsc%tridx(img, jat, iat))
         call get_amat_dir(vec, gam, ptr%alpha, ptr%dtrans, dtmp)
         call get_amat_rec_3d(vec, ptr%rtrans, rec_weight, rtmp)
         tmp = tmp + (dtmp + rtmp) * wsw
      end do
      if (iat == jat) then
         tmp = tmp + self%eta(izp) + sqrt2pi / self%rad(izp) &
            & - 2.0_wp * ptr%alpha / sqrtpi
      end if
      amat(jat, iat) = tmp
      amat(iat, jat) = tmp
   end do
   !$omp end parallel do

   amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
   amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
   amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

end subroutine get_amat_pbc


!> Direct-space Ewald contribution to the Coulomb matrix.
!> This short-range sum is identical for 1D, 2D, and 3D periodicity.
subroutine get_amat_dir(rij, gam, alp, trans, amat)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Effective Gaussian charge width parameter
   real(wp), intent(in) :: gam

   !> Ewald splitting parameter
   real(wp), intent(in) :: alp

   !> Direct-space lattice translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Direct-space Coulomb matrix element contribution
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

end subroutine get_amat_dir


!> Reciprocal-space contribution to the Coulomb matrix for 3D Ewald
pure subroutine get_amat_rec_3d(rij, trans, weight, amat)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space weights
   real(wp), intent(in) :: weight(:)

   !> Reciprocal-space Coulomb matrix element contribution
   real(wp), intent(out) :: amat

   integer :: itr

   amat = 0.0_wp
   do itr = 1, size(trans, 2)
      amat = amat + weight(itr) * cos(dot_product(rij, trans(:, itr)))
   end do

end subroutine get_amat_rec_3d


!> Evaluate derivatives of the EEQ Coulomb interaction matrix
subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Cache container
   type(cache_container), intent(inout) :: cache

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. atomic positions
   real(wp), intent(out) :: dadr(:, :, :)

   !> Derivative of Coulomb matrix w.r.t. lattice strain
   real(wp), intent(out) :: dadL(:, :, :)

   !> On-site trace of position derivatives
   real(wp), intent(out) :: atrace(:, :)

   type(eeq_cache), pointer :: ptr

   call view(cache, ptr)

   if (any(mol%periodic)) then
      call get_damat_pbc(self, mol, ptr, qvec, dadr, dadL, atrace)
   else
      call get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
   end if

end subroutine get_coulomb_derivs


!> Evaluate molecular (0D) Coulomb matrix derivatives
subroutine get_damat_0d(self, mol, qvec, dadr, dadL, atrace)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. atomic positions
   real(wp), intent(out) :: dadr(:, :, :)

   !> Derivative of Coulomb matrix w.r.t. lattice strain
   real(wp), intent(out) :: dadL(:, :, :)

   !> On-site trace of position derivatives
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, ipair, npair
   real(wp) :: vec(3), r2, gam, arg, dtmp, dg(3), ds(3, 3)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   npair = mol%nat * (mol%nat + 1) / 2

   ! Pair derivatives have unique owners. Only atom-wise traces and strain
   ! derivatives need reduction, with O(N) storage per thread.
   !$omp parallel default(none) shared(dadr, mol, self, qvec, npair) &
   !$omp& private(iat, jat, izp, jzp, ipair, gam, r2, vec, dg, ds, dtmp, arg) &
   !$omp& reduction(+:atrace, dadL)
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
      dg = dtmp * vec
      ds = outerprod(vec, dg)
      atrace(:, iat) = atrace(:, iat) - dg * qvec(jat)
      atrace(:, jat) = atrace(:, jat) + dg * qvec(iat)
      dadr(:, iat, jat) = -dg * qvec(iat)
      dadr(:, jat, iat) = +dg * qvec(jat)
      dadL(:, :, jat) = dadL(:, :, jat) + ds * qvec(iat)
      dadL(:, :, iat) = dadL(:, :, iat) + ds * qvec(jat)
   end do
   !$omp end do
   !$omp end parallel

end subroutine get_damat_0d


!> Evaluate periodic Ewald Coulomb matrix derivatives
subroutine get_damat_pbc(self, mol, ptr, qvec, dadr, dadL, atrace)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> EEQ model cache
   type(eeq_cache), intent(in) :: ptr

   !> Partial charge vector
   real(wp), intent(in) :: qvec(:)

   !> Derivative of Coulomb matrix w.r.t. atomic positions
   real(wp), intent(out) :: dadr(:, :, :)

   !> Derivative of Coulomb matrix w.r.t. lattice strain
   real(wp), intent(out) :: dadL(:, :, :)

   !> On-site trace of position derivatives
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: gam, wsw, vec(3), dg(3), ds(3, 3)
   real(wp) :: dgd(3), dsd(3, 3), dgr(3), dsr(3, 3)
   real(wp), allocatable :: rec_weight(:), rec_strain(:, :, :)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   call get_rec_weights(ptr%vol, ptr%alpha, ptr%rtrans, rec_weight, rec_strain)

   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel default(none) &
   !$omp& shared(mol, self, ptr, rec_weight, rec_strain, qvec, dadr, npair) &
   !$omp& private(iat, jat, izp, jzp, img, ipair, gam, wsw, vec, dg, ds, &
   !$omp& dgr, dsr, dgd, dsd) &
   !$omp& reduction(+:atrace, dadL)
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
      dg(:) = 0.0_wp
      ds(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      wsw = 1.0_wp / real(ptr%wsc%nimg(jat, iat), wp)
      do img = 1, ptr%wsc%nimg(jat, iat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat) &
            & + ptr%wsc%trans(:, ptr%wsc%tridx(img, jat, iat))
         call get_damat_dir(vec, gam, ptr%alpha, ptr%dtrans, dgd, dsd)
         call get_damat_rec_3d(vec, ptr%rtrans, rec_weight, rec_strain, dgr, dsr)
         dg = dg + (dgd + dgr) * wsw
         ds = ds + (dsd + dsr) * wsw
      end do
      if (iat == jat) then
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

end subroutine get_damat_pbc


!> Direct-space Ewald derivatives of the Coulomb matrix.
!> This short-range derivative sum is identical for 1D, 2D, and 3D periodicity.
subroutine get_damat_dir(rij, gam, alp, trans, dg, ds)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Effective Gaussian charge width parameter
   real(wp), intent(in) :: gam

   !> Ewald splitting parameter
   real(wp), intent(in) :: alp

   !> Direct-space lattice translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Position derivative contribution
   real(wp), intent(out) :: dg(3)

   !> Strain derivative contribution
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
      gtmp = +2.0_wp * gam * exp(-r2 * gam2) / (sqrtpi * r2) &
         & - erf(r1 * gam) / (r2 * r1)
      atmp = -2.0_wp * alp * exp(-r2 * alp2) / (sqrtpi * r2) &
         & + erf(r1 * alp) / (r2 * r1)
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * outerprod(vec, vec)
   end do

end subroutine get_damat_dir


!> Reciprocal-space Ewald derivatives of the Coulomb matrix for 3D periodicity
pure subroutine get_damat_rec_3d(rij, trans, weight, strain, dg, ds)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space weights
   real(wp), intent(in) :: weight(:)

   !> Reciprocal-space strain tensors
   real(wp), intent(in) :: strain(:, :, :)

   !> Position derivative contribution
   real(wp), intent(out) :: dg(3)

   !> Strain derivative contribution
   real(wp), intent(out) :: ds(3, 3)

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

   !> Unit cell volume
   real(wp), intent(in) :: vol

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space weights
   real(wp), allocatable, intent(out) :: weight(:)

   !> Optional reciprocal-space strain factors
   real(wp), allocatable, intent(out), optional :: strain(:, :, :)

   integer :: itr
   real(wp) :: g2, inv_g2, alp2, fac, vec(3)

   fac = 4.0_wp * pi / vol
   alp2 = alpha**2
   allocate(weight(size(trans, 2)), source=0.0_wp)
   if (present(strain)) then
      allocate(strain(3, 3, size(trans, 2)), source=0.0_wp)
   end if

   do itr = 1, size(trans, 2)
      vec = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle

      inv_g2 = 1.0_wp / g2
      weight(itr) = fac * exp(-0.25_wp * g2 / alp2) * inv_g2
      if (present(strain)) then
         strain(:, :, itr) = weight(itr) * ((2.0_wp * inv_g2 + 0.5_wp / alp2) &
            & * outerprod(vec, vec) - unity)
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
