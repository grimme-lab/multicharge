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

!> @file multicharge/model/eeqbc.f90
!> Provides implementation of the bond capacitor electronegativity equilibration model (EEQ_BC)

!> Bond capacitor electronegativity equilibration charge model published in
!>
!> Thomas Froitzheim, Marcel Müller, Andreas Hansen, and Stefan Grimme,
!> *J. Chem. Phys.*, **2025**, 162, 214109.
!> DOI: [10.1063/5.0268978](https://dx.doi.org/10.1063/5.0268978)
module multicharge_model_eeqbc
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: outerprod
   use mctc_ncoord, only: new_ncoord, cn_count
   use mctc_wignerseitz, only: wignerseitz_cell, new_wignerseitz_cell
   use multicharge_blas, only: gemv, gemm
   use multicharge_model_cache, only: mchrg_cache, ensure_shape
   use multicharge_model_hessian, only: coulomb_kernel, add_pair_hessian, add_chain_hessian
   use multicharge_model_type, only: mchrg_model_type, get_dir_trans
   implicit none
   private

   public :: eeqbc_model, new_eeqbc_model

   type, extends(mchrg_model_type) :: eeqbc_model
      !> Bond capacitance
      real(wp), allocatable :: cap(:)
      !> Average coordination number
      real(wp), allocatable :: avg_cn(:)
      !> Exponent of error function in bond capacitance
      real(wp) :: kbc
      !> Exponent of the distance/CN normalization
      real(wp) :: norm_exp
      !> vdW radii
      real(wp), allocatable :: rvdw(:, :)
   contains
      !> Update and allocate cache
      procedure :: update
      !> Calculate capacitance matrix
      procedure :: get_capacitance_matrix
      !> Calculate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Calculate derivatives of Coulomb matrix
      procedure :: get_coulomb_derivs
      !> Calculate right-hand side (electronegativity vector)
      procedure :: get_xvec
      !> Calculate derivatives of EN vector
      procedure :: get_xvec_derivs
      !> Add u^T (x_ab - A_ab q) at fixed lattice
      procedure :: add_model_hessian
      !> Calculate constraint matrix (molecular)
      procedure :: get_cmat_0d
      !> Calculate full constraint matrix (periodic)
      procedure :: get_cmat_3d
      !> Calculate constraint matrix derivatives (molecular)
      procedure :: get_dcmat_0d
      !> Calculate constraint matrix derivatives (periodic)
      procedure :: get_dcmat_3d
   end type eeqbc_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp / pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Default exponent of distance/CN normalization
   real(wp), parameter :: default_norm_exp = 1.0_wp

   !> Default exponent of error function in bond capacitance
   real(wp), parameter :: default_kbc = 0.65_wp

contains

subroutine new_eeqbc_model(self, mol, error, chi, rad, &
   & eta, kcnchi, kqchi, kqeta, kcnrad, cap, avg_cn, rvdw, &
   & kbc, cutoff, cn_exp, rcov, en, cn_max, norm_exp)
   !> Bond capacitor electronegativity equilibration model
   type(eeqbc_model), intent(out) :: self
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
   !> Local charge scaling factor for electronegativity
   real(wp), intent(in) :: kqchi(:)
   !> Local charge scaling factor for chemical hardness
   real(wp), intent(in) :: kqeta(:)
   !> CN scaling factor for charge width
   real(wp), intent(in) :: kcnrad
   !> Bond capacitance
   real(wp), intent(in) :: cap(:)
   !> Average coordination number
   real(wp), intent(in) :: avg_cn(:)
   !> Van-der-Waals radii
   real(wp), intent(in) :: rvdw(:, :)
   !> Exponent of error function in bond capacitance
   real(wp), intent(in), optional :: kbc
   !> Exponent of the distance normalization
   real(wp), intent(in), optional :: norm_exp
   !> Cutoff radius for coordination number
   real(wp), intent(in), optional :: cutoff
   !> Steepness of the CN counting function
   real(wp), intent(in), optional :: cn_exp
   !> Covalent radii for CN
   real(wp), intent(in), optional :: rcov(:)
   !> Maximum CN cutoff for CN
   real(wp), intent(in), optional :: cn_max
   !> Pauling electronegativities normalized to fluorine
   real(wp), intent(in), optional :: en(:)

   self%chi = chi
   self%rad = rad
   self%eta = eta
   self%kcnchi = kcnchi
   self%kqchi = kqchi
   self%kqeta = kqeta
   self%kcnrad = kcnrad
   self%cap = cap
   self%avg_cn = avg_cn
   self%rvdw = rvdw

   if (present(kbc)) then
      self%kbc = kbc
   else
      self%kbc = default_kbc
   end if

   if (present(norm_exp)) then
      self%norm_exp = norm_exp
   else
      self%norm_exp = default_norm_exp
   end if

   ! Coordination number
   call new_ncoord(self%ncoord, mol, cn_count%erf, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, cut=cn_max, &
      & norm_exp=self%norm_exp)
   if (allocated(error)) return

   ! Electronegativity weighted coordination number for local charge
   call new_ncoord(self%ncoord_en, mol, cn_count%erf_en, error, &
      & cutoff=cutoff, kcn=cn_exp, rcov=rcov, en=en, cut=cn_max, &
      & norm_exp=self%norm_exp)

end subroutine new_eeqbc_model

!> Update coordination numbers and local charges, set up the Wigner-Seitz cell
subroutine update(self, mol, error, cache, trans, grad)

   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self

   !> Structure type
   type(structure_type), intent(in) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   !> Lattice vectors
   real(wp), intent(in) :: trans(:, :)

   !> Flag to compute derivatives
   logical, intent(in) :: grad

   ! Allocate CN and local charge arrays in cache
   call ensure_shape(cache%cn, mol%nat)
   call ensure_shape(cache%qloc, mol%nat)

   cache%grad = grad

   ! Refer CN and local charge arrays in cache
   if (grad) then
      call ensure_shape(cache%dcndr, 3, mol%nat, mol%nat)
      call ensure_shape(cache%dcndL, 3, 3, mol%nat)
      call ensure_shape(cache%dqlocdr, 3, mol%nat, mol%nat)
      call ensure_shape(cache%dqlocdL, 3, 3, mol%nat)
      call self%ncoord%get_coordination_number(mol, trans, cache%cn, &
         & cache%dcndr, cache%dcndL)
      call self%local_charge(mol, trans, cache%qloc, cache%dqlocdr, cache%dqlocdL)
   else
      call self%ncoord%get_coordination_number(mol, trans, cache%cn)
      call self%local_charge(mol, trans, cache%qloc)
   end if

   if (any(mol%periodic)) then
      ! Set up the Wigner-Seitz cell
      call new_wignerseitz_cell(cache%wsc, mol, error)
      if (allocated(error)) return
   end if

end subroutine update


!> Compute the capacitance matrix (and its derivatives if needed).
subroutine get_capacitance_matrix(self, mol, ndim, cache)

   !> EEQBC model type
   class(eeqbc_model), intent(in) :: self

   !> Structure type
   type(structure_type), intent(in) :: mol

   !> System size
   integer, intent(in) :: ndim

   !> Multicharge cache
   type(mchrg_cache), intent(inout) :: cache

   call ensure_shape(cache%cmat, ndim, ndim)

   if (cache%grad) then
      call ensure_shape(cache%dcdr, 3, mol%nat, ndim)
      call ensure_shape(cache%dcdL, 3, 3, ndim)
   end if

   if (any(mol%periodic)) then
      ! Sum over all Wigner-Seitz images
      call get_cmat_3d(self, mol, cache%wsc, cache%cmat)
      if (cache%grad) then
         call get_dcmat_3d(self, mol, cache%wsc, cache%dcdr, cache%dcdL)
      end if
   else
      call get_cmat_0d(self, mol, cache%cmat)
      if (cache%grad) then
         call get_dcmat_0d(self, mol, cache%dcdr, cache%dcdL)
      end if
   end if

end subroutine get_capacitance_matrix

subroutine get_xvec(self, mol, ndim, cache, xvec)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ndim
   type(mchrg_cache), intent(inout) :: cache
   real(wp), intent(out) :: xvec(:)

   integer :: iat, izp, img
   real(wp) :: ctmp, csum, vec(3), rvdw, capi, wsw
   real(wp), allocatable :: dtrans(:, :)


   call ensure_shape(cache%xtmp, ndim)

   xvec(:) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(mol, self, cache) private(iat, izp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      cache%xtmp(iat) = -self%chi(izp) + self%kcnchi(izp) * cache%cn(iat) &
         & + self%kqchi(izp) * cache%qloc(iat)
   end do
   !$omp end parallel do
   if (ndim == mol%nat + 1) cache%xtmp(mol%nat + 1) = mol%charge

   call gemv(cache%cmat, cache%xtmp, xvec)

   if (any(mol%periodic)) then
      call get_dir_trans(mol, dtrans, self%ncoord%cutoff)
      !$omp parallel do default(none) schedule(runtime) &
      !$omp& shared(mol, self, cache, xvec, dtrans) &
      !$omp& private(iat, izp, img, wsw, capi, vec, rvdw, ctmp, csum)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         ! eliminate self-interaction (quasi off-diagonal)
         rvdw = self%rvdw(izp, izp)
         wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
         csum = 0.0_wp
         do img = 1, cache%wsc%nimg(iat, iat)
            vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))

            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
            csum = csum + ctmp
         end do
         xvec(iat) = xvec(iat) - wsw * csum * cache%xtmp(iat)
      end do
      !$omp end parallel do
   end if
end subroutine get_xvec

subroutine get_xvec_derivs(self, mol, ndim, cache, dxdr, dxdL)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ndim
   type(mchrg_cache), intent(inout) :: cache
   real(wp), intent(out), contiguous :: dxdr(:, :, :)
   real(wp), intent(out), contiguous :: dxdL(:, :, :)

   integer :: iat, izp, jat, jzp, img
   real(wp) :: capi, capj, wsw, rij(3), vec(3), ctmp, csum, rvdw, xij
   real(wp) :: dG(3), dS(3, 3), dS_sum(3, 3)
   real(wp), allocatable :: dtmpdr(:, :, :), dtmpdL(:, :, :)
   real(wp), allocatable :: dtrans(:, :)
   real(wp), allocatable :: dxdiag(:, :), dxdL_red(:, :, :)


   ! Leave uninitialized to preserve first-touch page locality;
   ! slices 1:mol%nat are fully assigned in the parallel loop.
   allocate(dtmpdr(3, mol%nat, ndim), dtmpdL(3, 3, ndim))

   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(mol, self, cache, dxdr, dxdL, dtmpdr, dtmpdL) &
   !$omp& private(iat, izp)
   do iat = 1, mol%nat
      izp = mol%id(iat)

      ! Zero intent(out) arrays for first-touch locality and intent(inout) gemm
      dxdr(:, :, iat) = 0.0_wp
      dxdL(:, :, iat) = 0.0_wp

      ! CN and effective charge derivative
      dtmpdr(:, :, iat) = self%kcnchi(izp) * cache%dcndr(:, :, iat) &
         & + self%kqchi(izp) * cache%dqlocdr(:, :, iat)
      dtmpdL(:, :, iat) = self%kcnchi(izp) * cache%dcndL(:, :, iat) &
         & + self%kqchi(izp) * cache%dqlocdL(:, :, iat)
   end do
   !$omp end parallel do

   if (ndim > mol%nat) then
      dxdr(:, :, mol%nat + 1) = 0.0_wp
      dxdL(:, :, mol%nat + 1) = 0.0_wp
      dtmpdr(:, :, mol%nat + 1) = 0.0_wp
      dtmpdL(:, :, mol%nat + 1) = 0.0_wp
   end if

   call gemm(dtmpdr, cache%cmat, dxdr(:, :, :mol%nat))
   call gemm(dtmpdL, cache%cmat, dxdL(:, :, :mol%nat))

   allocate(dxdiag(3, mol%nat), dxdL_red(3, 3, mol%nat), source=0.0_wp)

   if (any(mol%periodic)) then
      call get_dir_trans(mol, dtrans, self%ncoord%cutoff)

      ! Capacitance terms for i = j, T != 0
      !$omp parallel do default(none) schedule(runtime) &
      !$omp& shared(mol, self, cache, dtrans, dtmpdr, dtmpdL, dxdr, dxdL) &
      !$omp& private(iat, izp, img, capi, rvdw, wsw, vec, ctmp, csum)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         rvdw = self%rvdw(izp, izp)
         wsw = 1.0_wp / real(cache%wsc%nimg(iat, iat), wp)
         csum = 0.0_wp
         do img = 1, cache%wsc%nimg(iat, iat)
            vec = cache%wsc%trans(:, cache%wsc%tridx(img, iat, iat))
            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
            csum = csum + ctmp
         end do
         csum = csum * wsw
         ! EN derivative
         dxdr(:, :, iat) = dxdr(:, :, iat) - csum * dtmpdr(:, :, iat)
         dxdL(:, :, iat) = dxdL(:, :, iat) - csum * dtmpdL(:, :, iat)
      end do
      !$omp end parallel do

      ! Off-diagonal pairs (jat < iat) have unique owners; reduce only O(N) terms
      !$omp parallel do default(none) schedule(guided) &
      !$omp& shared(mol, self, cache, dtrans, dxdr) &
      !$omp& private(iat, izp, jat, jzp, img, capi, capj, rvdw, wsw, rij, vec, &
      !$omp& dG, dS, dS_sum, xij) &
      !$omp& reduction(+:dxdiag, dxdL_red)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(izp, jzp)
            xij = cache%xtmp(iat) - cache%xtmp(jat)

            ! Diagonal elements
            dxdiag(:, iat) = dxdiag(:, iat) + cache%xtmp(jat) * cache%dcdr(:, iat, jat)
            dxdiag(:, jat) = dxdiag(:, jat) + cache%xtmp(iat) * cache%dcdr(:, jat, iat)

            ! Derivative of capacitance matrix
            dxdr(:, iat, jat) = dxdr(:, iat, jat) + xij * cache%dcdr(:, iat, jat)
            dxdr(:, jat, iat) = dxdr(:, jat, iat) - xij * cache%dcdr(:, jat, iat)

            rij = mol%xyz(:, jat) - mol%xyz(:, iat)
            wsw = 1.0_wp / real(cache%wsc%nimg(jat, iat), wp)
            dS_sum(:, :) = 0.0_wp
            do img = 1, cache%wsc%nimg(jat, iat)
               vec = rij + cache%wsc%trans(:, cache%wsc%tridx(img, jat, iat))
               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, &
                  & dG, dS)
               dS_sum = dS_sum + dS
            end do
            dS_sum = dS_sum * wsw
            dxdL_red(:, :, iat) = dxdL_red(:, :, iat) + xij * dS_sum
            dxdL_red(:, :, jat) = dxdL_red(:, :, jat) - xij * dS_sum
         end do
         dxdiag(:, iat) = dxdiag(:, iat) + cache%xtmp(iat) * cache%dcdr(:, iat, iat)
      end do
      !$omp end parallel do
   else
      ! Off-diagonal pairs (jat < iat) have unique owners; reduce only O(N) terms
      !$omp parallel do default(none) schedule(guided) &
      !$omp& shared(mol, cache, dxdr) &
      !$omp& private(iat, jat, vec, xij) &
      !$omp& reduction(+:dxdiag, dxdL_red)
      do iat = 1, mol%nat
         do jat = 1, iat - 1
            xij = cache%xtmp(iat) - cache%xtmp(jat)

            ! Diagonal elements
            dxdiag(:, iat) = dxdiag(:, iat) + cache%xtmp(jat) * cache%dcdr(:, iat, jat)
            dxdiag(:, jat) = dxdiag(:, jat) + cache%xtmp(iat) * cache%dcdr(:, jat, iat)

            ! Derivative of capacitance matrix
            dxdr(:, iat, jat) = dxdr(:, iat, jat) + xij * cache%dcdr(:, iat, jat)
            dxdr(:, jat, iat) = dxdr(:, jat, iat) - xij * cache%dcdr(:, jat, iat)

            vec = mol%xyz(:, iat) - mol%xyz(:, jat)
            dxdL_red(:, :, iat) = dxdL_red(:, :, iat) &
               & + cache%xtmp(jat) * outerprod(vec, cache%dcdr(:, iat, jat))
            dxdL_red(:, :, jat) = dxdL_red(:, :, jat) &
               & + cache%xtmp(iat) * outerprod(-vec, cache%dcdr(:, jat, iat))
         end do
         dxdiag(:, iat) = dxdiag(:, iat) + cache%xtmp(iat) * cache%dcdr(:, iat, iat)
         dxdL_red(:, :, iat) = dxdL_red(:, :, iat) &
            & + cache%xtmp(iat) * cache%dcdL(:, :, iat)
      end do
      !$omp end parallel do
   end if

   !$omp parallel do default(none) schedule(static) &
   !$omp& shared(mol, dxdr, dxdL, dxdiag, dxdL_red) private(iat)
   do iat = 1, mol%nat
      dxdr(:, iat, iat) = dxdr(:, iat, iat) + dxdiag(:, iat)
      dxdL(:, :, iat) = dxdL(:, :, iat) + dxdL_red(:, :, iat)
   end do
   !$omp end parallel do

end subroutine get_xvec_derivs

!> Differentiate u^T(x-Aq), holding u and q fixed. Each bond contributes
!> c_ij [(u_i-u_j)(v_i-v_j) - u_i q_i h_i - u_j q_j h_j
!>       + (u_i q_j+u_j q_i) J_ij]. A self image uses half this weight,
!> yielding c_ii u_i q_i (J_ii-h_i). Its position is fixed but its width
!> depends on CN. Collect chain rules before Cartesian contraction.
subroutine add_model_hessian(self, mol, error, cache, adj, hessian)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(error_type), allocatable, intent(out) :: error
   type(mchrg_cache), intent(inout) :: cache
   real(wp), intent(in) :: adj(:)
   real(wp), intent(inout) :: hessian(:, :)

   integer :: nat, iat, jat, izp, jzp, img, itr, nimg
   logical :: periodic
   real(wp) :: rij(3), vec(3), unit(3), r, gam, gi, gj, gii, gjj, gij, wsw
   real(wp) :: val, dr, drr, dg, drg, dgg, cap, dcap, ddcap, k, t
   real(wp) :: udiff, piq, pjq, cross, term, ti, tj, li, lj, ri, rj, cij
   real(wp), allocatable :: width(:), dw(:), v(:), hard(:), dh(:), ddh(:)
   real(wp), allocatable :: dEdcn(:), dEdqloc(:), second(:, :), trans(:, :), dtrans(:, :)
   real(wp), allocatable :: mixed(:, :, :), mixed_loc(:, :, :)

   if (.not. allocated(cache%dcndr) .or. .not. allocated(cache%dqlocdr)) then
      call fatal_error(error, "Incompatible EEQBC response cache")
      return
   end if
   nat = mol%nat
   periodic = any(mol%periodic)
   if (periodic) then
      call get_dir_trans(mol, dtrans, self%ncoord%cutoff)
   else
      allocate(dtrans(3, 1), source=0.0_wp)
   end if
   allocate(width(nat), dw(nat), v(nat), hard(nat), dh(nat), ddh(nat))
   allocate(dEdcn(nat), dEdqloc(nat), second(nat, nat))
   allocate(mixed(3, nat, nat), mixed_loc(3, nat, nat))
   dEdcn = 0.0_wp
   dEdqloc = 0.0_wp
   second = 0.0_wp
   mixed = 0.0_wp
   mixed_loc = 0.0_wp
   do iat = 1, nat
      izp = mol%id(iat)
      dw(iat) = -self%rad(izp)*self%kcnrad/self%avg_cn(izp)**self%norm_exp
      width(iat) = self%rad(izp)+dw(iat)*cache%cn(iat)
      v(iat) = -self%chi(izp)+self%kcnchi(izp)*cache%cn(iat) &
         & + self%kqchi(izp)*cache%qloc(iat)
      hard(iat) = self%eta(izp)+self%kqeta(izp)*cache%qloc(iat)+sqrt2pi/width(iat)
      dh(iat) = -sqrt2pi*dw(iat)/width(iat)**2
      ddh(iat) = 2.0_wp*sqrt2pi*dw(iat)**2/width(iat)**3
   end do

   do iat = 1, nat
      izp = mol%id(iat)
      do jat = 1, iat
         if (iat == jat .and. .not. periodic) cycle
         jzp = mol%id(jat)
         gam = 1.0_wp/sqrt(width(iat)**2+width(jat)**2)
         gi = -width(iat)*dw(iat)*gam**3
         gj = -width(jat)*dw(jat)*gam**3
         gii = -dw(iat)**2*gam**3+3.0_wp*(width(iat)*dw(iat))**2*gam**5
         gjj = -dw(jat)**2*gam**3+3.0_wp*(width(jat)*dw(jat))**2*gam**5
         gij = 3.0_wp*width(iat)*width(jat)*dw(iat)*dw(jat)*gam**5
         udiff = adj(iat)-adj(jat)
         piq = adj(iat)*cache%vrhs(iat)
         pjq = adj(jat)*cache%vrhs(jat)
         cross = adj(iat)*cache%vrhs(jat)+adj(jat)*cache%vrhs(iat)
         li = udiff*self%kqchi(izp)-piq*self%kqeta(izp)
         lj = -udiff*self%kqchi(jzp)-pjq*self%kqeta(jzp)
         k = self%kbc/self%rvdw(izp, jzp)
         nimg = 1
         if (periodic) nimg = cache%wsc%nimg(jat, iat)
         wsw = 1.0_wp/real(nimg, wp)
         if (iat == jat) wsw = 0.5_wp*wsw
         do img = 1, nimg
            rij = mol%xyz(:, iat)-mol%xyz(:, jat)
            if (periodic) rij = rij-cache%wsc%trans(:, cache%wsc%tridx(img, jat, iat))
            do itr = 1, size(dtrans, 2)
               vec = rij+dtrans(:, itr)
               r = norm2(vec)
               if (r < eps) cycle
               call coulomb_kernel(r, gam, val, dr, drr, dg, drg, dgg)
               t = -k*(r-self%rvdw(izp, jzp))
               cap = wsw*sqrt(self%cap(izp)*self%cap(jzp))
               dcap = -cap*k/sqrtpi*exp(-t*t)
               ddcap = 2.0_wp*k*t*dcap
               cap = 0.5_wp*cap*(1.0_wp+erf(t))

               ti = udiff*self%kcnchi(izp)-piq*dh(iat)+cross*dg*gi
               tj = -udiff*self%kcnchi(jzp)-pjq*dh(jat)+cross*dg*gj
               dEdcn(iat) = dEdcn(iat)+cap*ti
               dEdcn(jat) = dEdcn(jat)+cap*tj
               dEdqloc(iat) = dEdqloc(iat)+cap*li
               dEdqloc(jat) = dEdqloc(jat)+cap*lj
               second(iat, iat) = second(iat, iat)+cap*(-piq*ddh(iat)+cross*(dgg*gi*gi+dg*gii))
               second(jat, jat) = second(jat, jat)+cap*(-pjq*ddh(jat)+cross*(dgg*gj*gj+dg*gjj))
               cij = cap*cross*(dgg*gi*gj+dg*gij)
               second(iat, jat) = second(iat, jat)+cij
               second(jat, iat) = second(jat, iat)+cij

               ! Self images contribute through CN and local charge only.
               if (iat == jat) cycle
               term = udiff*(v(iat)-v(jat))-piq*hard(iat)-pjq*hard(jat)+cross*val
               call add_pair_hessian(iat, jat, vec, dcap*term+cap*cross*dr, &
                  & ddcap*term+2.0_wp*dcap*cross*dr+cap*cross*drr, hessian)
               unit = vec/r
               ri = dcap*ti+cap*cross*drg*gi
               rj = dcap*tj+cap*cross*drg*gj
               mixed(:, iat, iat) = mixed(:, iat, iat)+unit*ri
               mixed(:, jat, iat) = mixed(:, jat, iat)-unit*ri
               mixed(:, iat, jat) = mixed(:, iat, jat)+unit*rj
               mixed(:, jat, jat) = mixed(:, jat, jat)-unit*rj
               mixed_loc(:, iat, iat) = mixed_loc(:, iat, iat)+unit*dcap*li
               mixed_loc(:, jat, iat) = mixed_loc(:, jat, iat)-unit*dcap*li
               mixed_loc(:, iat, jat) = mixed_loc(:, iat, jat)+unit*dcap*lj
               mixed_loc(:, jat, jat) = mixed_loc(:, jat, jat)-unit*dcap*lj
            end do
         end do
      end do
   end do

   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, trans, error)
   if (allocated(error)) return
   call self%ncoord%add_coordination_number_hessian(mol, trans, dEdcn, hessian)
   call self%ncoord_en%add_coordination_number_hessian(mol, trans, dEdqloc, hessian)
   call add_chain_hessian(cache%dcndr, mixed, hessian, second)
   call add_chain_hessian(cache%dqlocdr, mixed_loc, hessian)
end subroutine add_model_hessian

subroutine get_coulomb_matrix(self, mol, ndim, cache)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ndim
   type(mchrg_cache), intent(inout) :: cache

   call ensure_shape(cache%amat, ndim, ndim)

   if (any(mol%periodic)) then
      call get_amat_3d(self, mol, cache%wsc, cache%cn, cache%qloc, cache%cmat, cache%amat)
   else
      call get_amat_0d(self, mol, cache%cn, cache%qloc, cache%cmat, cache%amat)
   end if
end subroutine get_coulomb_matrix

subroutine get_amat_0d(self, mol, cn, qloc, cmat, amat)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: qloc(:)
   real(wp), intent(in) :: cmat(:, :)
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam2, tmp, norm_cn, radi, radj

   ! Each iat iteration owns the off-diagonal pairs (jat < iat) and diagonal (iat, iat)
   !$omp parallel do default(none) schedule(guided) &
   !$omp& shared(amat, mol, self, cn, qloc, cmat) &
   !$omp& private(iat, izp, jat, jzp, gam2, vec, r2, tmp, norm_cn, radi, radj)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
      radi = self%rad(izp) * (1.0_wp - self%kcnrad * cn(iat) * norm_cn)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! Effective charge width of j
         norm_cn = cn(jat) / self%avg_cn(jzp)**self%norm_exp
         radj = self%rad(jzp) * (1.0_wp - self%kcnrad * norm_cn)
         ! Coulomb interaction of Gaussian charges
         gam2 = 1.0_wp / (radi**2 + radj**2)
         tmp = erf(sqrt(r2 * gam2)) / sqrt(r2) * cmat(jat, iat)
         amat(jat, iat) = tmp
         amat(iat, jat) = tmp
      end do
      ! Effective hardness
      tmp = self%eta(izp) + self%kqeta(izp) * qloc(iat) + sqrt2pi / radi
      amat(iat, iat) = tmp * cmat(iat, iat) + 1.0_wp
   end do
   !$omp end parallel do

   if (size(amat, 1) == mol%nat + 1) then
      amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
      amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
      amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
   end if

end subroutine get_amat_0d

subroutine get_amat_3d(self, mol, wsc, cn, qloc, cmat, amat)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: cn(:), qloc(:), cmat(:, :)
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: rij(3), vec(3), gam, dtmp, tmp, capi, capj, radi, radj, norm_cn, rvdw, wsw
   real(wp), allocatable :: dtrans(:, :)

   call get_dir_trans(mol, dtrans, self%ncoord%cutoff)

   ! Each iat iteration owns the off-diagonal pairs (jat < iat) and diagonal (iat, iat)
   !$omp parallel do default(none) schedule(guided) &
   !$omp& shared(amat, cmat, mol, cn, qloc, self, wsc, dtrans) &
   !$omp& private(iat, izp, jat, jzp, img, gam, rij, vec, dtmp, tmp, norm_cn) &
   !$omp& private(radi, radj, capi, capj, rvdw, wsw)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
      radi = self%rad(izp) * (1.0_wp - self%kcnrad * cn(iat) * norm_cn)
      capi = self%cap(izp)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         ! vdw distance in Angstrom (approximate factor 2)
         rvdw = self%rvdw(izp, jzp)
         ! Effective charge width of j
         norm_cn = cn(jat) / self%avg_cn(jzp)**self%norm_exp
         radj = self%rad(jzp) * (1.0_wp - self%kcnrad * norm_cn)
         capj = self%cap(jzp)
         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         rij = mol%xyz(:, jat) - mol%xyz(:, iat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         tmp = 0.0_wp
         do img = 1, wsc%nimg(jat, iat)
            vec = rij + wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capj, dtmp)
            tmp = tmp + dtmp
         end do
         tmp = tmp * wsw
         amat(jat, iat) = tmp
         amat(iat, jat) = tmp
      end do

      ! diagonal Coulomb interaction terms
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      tmp = 0.0_wp
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capi, dtmp)
         tmp = tmp + dtmp
      end do

      ! Effective hardness
      dtmp = self%eta(izp) + self%kqeta(izp) * qloc(iat) + sqrt2pi / radi
      amat(iat, iat) = tmp * wsw + cmat(iat, iat) * dtmp + 1.0_wp
   end do
   !$omp end parallel do

   if (size(amat, 1) == mol%nat + 1) then
      amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
      amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
      amat(mol%nat + 1, mol%nat + 1) = 0.0_wp
   end if

end subroutine get_amat_3d

subroutine get_amat_dir_3d(rij, gam, trans, kbc, rvdw, capi, capj, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: kbc
   real(wp), intent(in) :: rvdw
   real(wp), intent(in) :: capi
   real(wp), intent(in) :: capj
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp, ctmp

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      call get_cpair(kbc, ctmp, r1, rvdw, capi, capj)
      tmp = -ctmp * erf(gam * r1) / r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d

subroutine get_coulomb_derivs(self, mol, ndim, cache, qvec, dadr, dadL)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   integer, intent(in) :: ndim
   type(mchrg_cache), intent(inout) :: cache
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out), contiguous :: dadr(:, :, :)
   real(wp), intent(out), contiguous :: dadL(:, :, :)

   integer :: iat
   real(wp), allocatable :: atrace(:, :)

   allocate(atrace(3, mol%nat))

   if (any(mol%periodic)) then
      call get_damat_3d(self, mol, cache%wsc, cache%cn, &
         & cache%qloc, qvec, cache%dcndr, cache%dcndL, cache%dqlocdr, &
         & cache%dqlocdL, cache%cmat, cache%dcdr, cache%dcdL, dadr, dadL, atrace)

   else
      call get_damat_0d(self, mol, cache%cn, &
         & cache%qloc, qvec, cache%dcndr, cache%dcndL, cache%dqlocdr, &
         & cache%dqlocdL, cache%cmat, cache%dcdr, cache%dcdL, dadr, dadL, atrace)
   end if
   do iat = 1, mol%nat
      dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
   end do
end subroutine get_coulomb_derivs

subroutine get_damat_0d(self, mol, cn, qloc, qvec, dcndr, dcndL, &
      & dqlocdr, dqlocdL, cmat, dcdr, dcdL, dadr, dadL, atrace)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: qloc(:)
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)
   real(wp), intent(in) :: dqlocdr(:, :, :)
   real(wp), intent(in) :: dqlocdL(:, :, :)
   real(wp), intent(in) :: cmat(:, :)
   real(wp), intent(in) :: dcdr(:, :, :)
   real(wp), intent(in) :: dcdL(:, :, :)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), dgamdL(3, 3)
   real(wp), allocatable :: dgamdr(:, :)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: atrace_local(:, :)
   real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

   allocate(dgamdr(3, mol%nat))

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(atrace, dadr, dadL, mol, self, cn, qloc, qvec) &
   !$omp shared(cmat, dcdr, dcdL, dcndr, dcndL, dqlocdr, dqlocdL) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
   !$omp private(radi, radj, dradi, dradj, dgamdr, dgamdL, dG, dS) &
   !$omp private(atrace_local, dadr_local, dadL_local)
   allocate(atrace_local, source=atrace)
   allocate(dadr_local, source=dadr)
   allocate(dadL_local, source=dadL)
   !$omp do schedule(guided)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
      radi = self%rad(izp) * (1.0_wp - self%kcnrad * cn(iat) * norm_cn)
      dradi = -self%rad(izp) * self%kcnrad * norm_cn
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
         radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cn(jat) * norm_cn)
         dradj = -self%rad(jzp) * self%kcnrad * norm_cn

         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         dgamdr(:, :) = -(radi * dradi * dcndr(:, :, iat) + radj * dradj * dcndr(:, :, jat)) &
                        & * gam**3.0_wp
         dgamdL(:, :) = -(radi * dradi * dcndL(:, :, iat) + radj * dradj * dcndL(:, :, jat)) &
                        & * gam**3.0_wp

         ! Explicit derivative
         arg = gam * gam * r2
         dtmp = 2.0_wp * gam * exp(-arg) / (sqrtpi * r2) &
            & - erf(sqrt(arg)) / (r2 * sqrt(r2))
         dG(:) = dtmp * vec
         dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
         atrace_local(:, iat) = -dG * qvec(jat) * cmat(jat, iat) + atrace_local(:, iat)
         atrace_local(:, jat) = +dG * qvec(iat) * cmat(iat, jat) + atrace_local(:, jat)
         dadr_local(:, iat, jat) = -dG * qvec(iat) * cmat(iat, jat) + dadr_local(:, iat, jat)
         dadr_local(:, jat, iat) = +dG * qvec(jat) * cmat(jat, iat) + dadr_local(:, jat, iat)
         dadL_local(:, :, iat) = +dS * qvec(jat) * cmat(jat, iat) + dadL_local(:, :, iat)
         dadL_local(:, :, jat) = +dS * qvec(iat) * cmat(iat, jat) + dadL_local(:, :, jat)

         ! Effective charge width derivative
         dtmp = 2.0_wp * exp(-arg) / (sqrtpi)
         ! The widths depend on CN, so every atom can contribute to this pair.
         dadr_local(:, :, iat) = dadr_local(:, :, iat) + dtmp * qvec(jat) * cmat(jat, iat) * dgamdr
         dadr_local(:, :, jat) = dadr_local(:, :, jat) + dtmp * qvec(iat) * cmat(iat, jat) * dgamdr
         dadL_local(:, :, iat) = +dtmp * qvec(jat) * dgamdL(:, :) * cmat(jat, iat) + dadL_local(:, :, iat)
         dadL_local(:, :, jat) = +dtmp * qvec(iat) * dgamdL(:, :) * cmat(iat, jat) + dadL_local(:, :, jat)

         ! Capacitance derivative off-diagonal
         dtmp = erf(sqrt(r2) * gam) / (sqrt(r2))
         atrace_local(:, iat) = -dtmp * qvec(jat) * dcdr(:, jat, iat) + atrace_local(:, iat)
         atrace_local(:, jat) = -dtmp * qvec(iat) * dcdr(:, iat, jat) + atrace_local(:, jat)
         dadr_local(:, iat, jat) = +dtmp * qvec(iat) * dcdr(:, iat, jat) + dadr_local(:, iat, jat)
         dadr_local(:, jat, iat) = +dtmp * qvec(jat) * dcdr(:, jat, iat) + dadr_local(:, jat, iat)
         dadL_local(:, :, iat) = -dtmp * qvec(jat) * spread(dcdr(:, iat, jat), 2, 3) * spread(vec, 1, 3) &
               & + dadL_local(:, :, iat)
         dadL_local(:, :, jat) = -dtmp * qvec(iat) * spread(dcdr(:, iat, jat), 2, 3) * spread(vec, 1, 3) &
               & + dadL_local(:, :, jat)

         ! Capacitance derivative diagonal
         dtmp = (self%eta(izp) + self%kqeta(izp) * qloc(iat) + sqrt2pi / radi) * qvec(iat)
         dadr_local(:, jat, iat) = -dtmp * dcdr(:, jat, iat) + dadr_local(:, jat, iat)

         dtmp = (self%eta(jzp) + self%kqeta(jzp) * qloc(jat) + sqrt2pi / radj) * qvec(jat)
         dadr_local(:, iat, jat) = -dtmp * dcdr(:, iat, jat) + dadr_local(:, iat, jat)
      end do

      ! Hardness derivative
      dtmp = self%kqeta(izp) * qvec(iat) * cmat(iat, iat)
      dadr_local(:, :, iat) = +dtmp * dqlocdr(:, :, iat) + dadr_local(:, :, iat)
      dadL_local(:, :, iat) = +dtmp * dqlocdL(:, :, iat) + dadL_local(:, :, iat)

      ! Effective charge width derivative
      dtmp = -sqrt2pi * dradi / (radi**2) * qvec(iat) * cmat(iat, iat)
      dadr_local(:, :, iat) = +dtmp * dcndr(:, :, iat) + dadr_local(:, :, iat)
      dadL_local(:, :, iat) = +dtmp * dcndL(:, :, iat) + dadL_local(:, :, iat)

      ! Capacitance derivative
      dtmp = (self%eta(izp) + self%kqeta(izp) * qloc(iat) + sqrt2pi / radi) * qvec(iat)
      dadr_local(:, iat, iat) = +dtmp * dcdr(:, iat, iat) + dadr_local(:, iat, iat)
      dadL_local(:, :, iat) = +dtmp * dcdL(:, :, iat) + dadL_local(:, :, iat)

   end do
   !$omp end do
   !$omp critical (get_damat_0d_)
   atrace(:, :) = atrace + atrace_local
   dadr(:, :, :) = dadr + dadr_local
   dadL(:, :, :) = dadL + dadL_local
   !$omp end critical (get_damat_0d_)
   deallocate(dadL_local, dadr_local, atrace_local)
   !$omp end parallel

end subroutine get_damat_0d

subroutine get_damat_3d(self, mol, wsc, cn, qloc, qvec, dcndr, dcndL, dqlocdr, &
   & dqlocdL, cmat, dcdr, dcdL, dadr, dadL, atrace)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: qloc(:)
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)
   real(wp), intent(in) :: dqlocdr(:, :, :)
   real(wp), intent(in) :: dqlocdL(:, :, :)
   real(wp), intent(in) :: cmat(:, :)
   real(wp), intent(in) :: dcdr(:, :, :)
   real(wp), intent(in) :: dcdL(:, :, :)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: rij(3), vec(3), r2, gam, arg, dtmp, norm_cn, rvdw, wsw, dgam
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
   real(wp) :: dgamdL(3, 3), capi, capj
   real(wp), allocatable :: dgamdr(:, :), dtrans(:, :)

   ! Thread-private arrays for reduction
   real(wp), allocatable :: atrace_local(:, :)
   real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

   call get_dir_trans(mol, dtrans, self%ncoord%cutoff)

   allocate(dgamdr(3, mol%nat))

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   !$omp parallel default(none) &
   !$omp shared(self, mol, cn, qloc, qvec, wsc, dadr, dadL, atrace) &
   !$omp shared (cmat, dcdr, dcdL, dcndr, dcndL, dqlocdr, dqlocdL, dtrans) &
   !$omp private(iat, izp, jat, jzp, img, gam, rij, vec, r2, dtmp, norm_cn, arg, rvdw) &
   !$omp private(radi, radj, dradi, dradj, capi, capj, dgamdr, dgamdL, dG, dS, wsw) &
   !$omp private(dgam, dadr_local, dadL_local, atrace_local)
   allocate(atrace_local, source=atrace)
   allocate(dadr_local, source=dadr)
   allocate(dadL_local, source=dadL)
   !$omp do schedule(guided)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
      radi = self%rad(izp) * (1.0_wp - self%kcnrad * cn(iat) * norm_cn)
      dradi = -self%rad(izp) * self%kcnrad * norm_cn
      capi = self%cap(izp)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)

         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
         radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cn(jat) * norm_cn)
         dradj = -self%rad(jzp) * self%kcnrad * norm_cn

         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         dgamdr(:, :) = -(radi * dradi * dcndr(:, :, iat) + radj * dradj * dcndr(:, :, jat)) &
                        & * gam**3.0_wp
         dgamdL(:, :) = -(radi * dradi * dcndL(:, :, iat) + radj * dradj * dcndL(:, :, jat)) &
                        & * gam**3.0_wp

         rij = mol%xyz(:, jat) - mol%xyz(:, iat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = rij + wsc%trans(:, wsc%tridx(img, jat, iat))

            call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)
            dG = dG * wsw
            dS = dS * wsw
            dgam = dgam * wsw

            ! Explicit derivative
            atrace_local(:, iat) = -dG * qvec(jat) + atrace_local(:, iat)
            atrace_local(:, jat) = +dG * qvec(iat) + atrace_local(:, jat)
            dadr_local(:, iat, jat) = -dG * qvec(iat) + dadr_local(:, iat, jat)
            dadr_local(:, jat, iat) = +dG * qvec(jat) + dadr_local(:, jat, iat)
            dadL_local(:, :, jat) = +dS * qvec(iat) + dadL_local(:, :, jat)
            dadL_local(:, :, iat) = +dS * qvec(jat) + dadL_local(:, :, iat)

            ! CN-dependent widths respond to every atom, including third atoms.
            dadr_local(:, :, iat) = dadr_local(:, :, iat) - dgam * qvec(jat) * dgamdr
            dadr_local(:, :, jat) = dadr_local(:, :, jat) - dgam * qvec(iat) * dgamdr
            dadL_local(:, :, iat) = -dgam * qvec(jat) * dgamdL(:, :) + dadL_local(:, :, iat)
            dadL_local(:, :, jat) = -dgam * qvec(iat) * dgamdL(:, :) + dadL_local(:, :, jat)

            call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
            dG = dG * wsw
            dS = dS * wsw

            ! Capacitance derivative off-diagonal
            atrace_local(:, iat) = +qvec(jat) * dG(:) + atrace_local(:, iat)
            atrace_local(:, jat) = -qvec(iat) * dG(:) + atrace_local(:, jat)
            dadr_local(:, jat, iat) = -qvec(jat) * dG(:) + dadr_local(:, jat, iat)
            dadr_local(:, iat, jat) = +qvec(iat) * dG(:) + dadr_local(:, iat, jat)
            dadL_local(:, :, jat) = -qvec(iat) * dS(:, :) + dadL_local(:, :, jat)
            dadL_local(:, :, iat) = -qvec(jat) * dS(:, :) + dadL_local(:, :, iat)

            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
            dG = dG * wsw

            ! Capacitance derivative diagonal
            dtmp = (self%eta(izp) + self%kqeta(izp) * qloc(iat) + sqrt2pi / radi) * qvec(iat)
            dadr_local(:, jat, iat) = +dtmp * dG(:) + dadr_local(:, jat, iat)
            dtmp = (self%eta(jzp) + self%kqeta(jzp) * qloc(jat) + sqrt2pi / radj) * qvec(jat)
            dadr_local(:, iat, jat) = -dtmp * dG(:) + dadr_local(:, iat, jat)
         end do
      end do

      ! diagonal explicit, charge width, and capacitance derivative terms
      gam = 1.0_wp / sqrt(2.0_wp * radi**2)
      dgamdr = -2.0_wp * radi * dradi * gam**3 * dcndr(:, :, iat)
      dgamdL = -2.0_wp * radi * dradi * gam**3 * dcndL(:, :, iat)
      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
         dgam = dgam * wsw

         ! Explicit derivative
         dadL_local(:, :, iat) = +dS * wsw * qvec(iat) + dadL_local(:, :, iat)

         ! Self-image distances are fixed under atomic displacements, but their
         ! Gaussian width still changes through the CN of this atom.
         dadr_local(:, :, iat) = dadr_local(:, :, iat) - dgam * qvec(iat) * dgamdr
         dadL_local(:, :, iat) = dadL_local(:, :, iat) - dgam * qvec(iat) * dgamdL

         ! Capacitance derivative
         call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
         dadL_local(:, :, iat) = -qvec(iat) * dS * wsw + dadL_local(:, :, iat)
      end do

      ! Hardness derivative
      dtmp = self%kqeta(izp) * qvec(iat) * cmat(iat, iat)
      dadr_local(:, :, iat) = +dtmp * dqlocdr(:, :, iat) + dadr_local(:, :, iat)
      dadL_local(:, :, iat) = +dtmp * dqlocdL(:, :, iat) + dadL_local(:, :, iat)

      ! Effective charge width derivative
      dtmp = -sqrt2pi * dradi / (radi**2) * qvec(iat) * cmat(iat, iat)
      dadr_local(:, :, iat) = +dtmp * dcndr(:, :, iat) + dadr_local(:, :, iat)
      dadL_local(:, :, iat) = +dtmp * dcndL(:, :, iat) + dadL_local(:, :, iat)

      dtmp = (self%eta(izp) + self%kqeta(izp) * qloc(iat) + sqrt2pi / radi) * qvec(iat)
      dadr_local(:, iat, iat) = +dtmp * dcdr(:, iat, iat) + dadr_local(:, iat, iat)
      dadL_local(:, :, iat) = +dtmp * dcdL(:, :, iat) + dadL_local(:, :, iat)

   end do
   !$omp end do
   !$omp critical (get_damat_3d_)
   atrace(:, :) = atrace + atrace_local
   dadr(:, :, :) = dadr + dadr_local
   dadL(:, :, :) = dadL + dadL_local
   !$omp end critical (get_damat_3d_)
   deallocate(dadL_local, dadr_local, atrace_local)
   !$omp end parallel

end subroutine get_damat_3d

subroutine get_damat_dir(rij, trans, capi, capj, rvdw, kbc, gam, dG, dS, dgam)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: capi, capj, rvdw, kbc
   real(wp), intent(out) :: dG(3)
   real(wp), intent(out) :: dS(3, 3)
   real(wp), intent(out) :: dgam

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, gam2, cmat

   dG(:) = 0.0_wp
   dS(:, :) = 0.0_wp
   dgam = 0.0_wp

   gam2 = gam * gam

   do itr = 1, size(trans, 2)
      vec(:) = rij(:) + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1 * r1
      call get_cpair(kbc, cmat, r1, rvdw, capi, capj)
      gtmp = 2.0_wp * gam * exp(-r2 * gam2) / (sqrtpi * r2) - erf(r1 * gam) / (r2 * r1)
      dG(:) = dG - cmat * gtmp * vec
      dS(:, :) = dS - cmat * gtmp * spread(vec, 1, 3) * spread(vec, 2, 3)
      dgam = dgam + cmat * 2.0_wp * exp(-gam2 * r2) / sqrtpi
   end do

end subroutine get_damat_dir

subroutine get_damat_dc_dir(rij, trans, capi, capj, rvdw, kbc, gam, dG, dS)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: capi, capj, rvdw, kbc
   real(wp), intent(out) :: dG(3)
   real(wp), intent(out) :: dS(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, gtmp(3), stmp(3, 3), tmp

   dG(:) = 0.0_wp
   dS(:, :) = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij(:) + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      call get_dcpair(kbc, vec, rvdw, capi, capj, gtmp, stmp)
      tmp = erf(gam * r1) / r1
      dG(:) = dG(:) + tmp * gtmp
      dS(:, :) = dS(:, :) + tmp * stmp
   end do

end subroutine get_damat_dc_dir

subroutine get_cmat_0d(self, mol, cmat)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: cmat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), rvdw, tmp, capi, capj, r1
   real(wp), allocatable :: cdiag(:)

   allocate(cdiag(mol%nat), source=0.0_wp)

   ! Off-diagonal pairs (jat < iat) have a unique owner; reduce only the 1D diagonal
   !$omp parallel default(none) &
   !$omp& shared(cmat, cdiag, mol, self) &
   !$omp& private(iat, izp, jat, jzp, vec, r1, rvdw, tmp, capi, capj)
   !$omp do schedule(guided) reduction(+:cdiag)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      capi = self%cap(izp)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         rvdw = self%rvdw(izp, jzp)
         capj = self%cap(jzp)

         call get_cpair(self%kbc, tmp, r1, rvdw, capi, capj)

         ! Off-diagonal elements
         cmat(jat, iat) = -tmp
         cmat(iat, jat) = -tmp
         ! Diagonal elements
         cdiag(iat) = cdiag(iat) + tmp
         cdiag(jat) = cdiag(jat) + tmp
      end do
   end do
   !$omp end do
   !$omp do schedule(static)
   do iat = 1, mol%nat
      cmat(iat, iat) = cdiag(iat)
   end do
   !$omp end do
   !$omp end parallel

   if (size(cmat, 1) == mol%nat + 1) then
      cmat(mol%nat + 1, 1:mol%nat) = 0.0_wp
      cmat(1:mol%nat, mol%nat + 1) = 0.0_wp
      cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp
   end if

end subroutine get_cmat_0d

subroutine get_cmat_3d(self, mol, wsc, cmat)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(out) :: cmat(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: rij(3), vec(3), rvdw, tmp, csum, capi, capj, wsw
   real(wp), allocatable :: dtrans(:, :), cdiag(:)

   call get_dir_trans(mol, dtrans, self%ncoord%cutoff)

   allocate(cdiag(mol%nat), source=0.0_wp)

   ! Off-diagonal pairs (jat < iat) have a unique owner; reduce only the 1D diagonal
   !$omp parallel default(none) &
   !$omp& shared(cmat, cdiag, mol, self, wsc, dtrans) &
   !$omp& private(iat, izp, jat, jzp, img, rij, vec, rvdw, tmp, csum, capi, capj, wsw)
   !$omp do schedule(guided) reduction(+:cdiag)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      capi = self%cap(izp)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         rvdw = self%rvdw(izp, jzp)
         capj = self%cap(jzp)
         rij = mol%xyz(:, iat) - mol%xyz(:, jat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         csum = 0.0_wp
         do img = 1, wsc%nimg(jat, iat)
            vec = rij - wsc%trans(:, wsc%tridx(img, jat, iat))

            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, tmp)
            csum = csum + tmp
         end do
         csum = csum * wsw

         ! Off-diagonal elements
         cmat(jat, iat) = -csum
         cmat(iat, jat) = -csum
         ! Diagonal elements
         cdiag(iat) = cdiag(iat) + csum
         cdiag(jat) = cdiag(jat) + csum
      end do

      ! diagonal capacitance (interaction with images)
      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      csum = 0.0_wp
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, tmp)
         csum = csum + tmp
      end do
      cdiag(iat) = cdiag(iat) + csum * wsw
   end do
   !$omp end do
   !$omp do schedule(static)
   do iat = 1, mol%nat
      cmat(iat, iat) = cdiag(iat)
   end do
   !$omp end do
   !$omp end parallel

   if (size(cmat, 1) == mol%nat + 1) then
      cmat(mol%nat + 1, 1:mol%nat) = 0.0_wp
      cmat(1:mol%nat, mol%nat + 1) = 0.0_wp
      cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp
   end if

end subroutine get_cmat_3d

subroutine get_cpair(kbc, cpair, r1, rvdw, capi, capj)
   real(wp), intent(in) :: kbc
   real(wp), intent(in) :: r1
   real(wp), intent(in) :: capi
   real(wp), intent(in) :: capj
   real(wp), intent(in) :: rvdw
   real(wp), intent(out) :: cpair

   real(wp) :: arg

   ! Capacitance of bond between atom i and j
   arg = -kbc * (r1 - rvdw) / rvdw
   cpair = sqrt(capi * capj) * 0.5_wp * (1.0_wp + erf(arg))
end subroutine get_cpair

subroutine get_cpair_dir(kbc, rij, trans, rvdw, capi, capj, cpair)
   real(wp), intent(in) :: kbc
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(in) :: rvdw
   real(wp), intent(in) :: capi
   real(wp), intent(in) :: capj
   real(wp), intent(out) :: cpair

   integer :: itr
   real(wp) :: vec(3), r1, tmp

   cpair = 0.0_wp
   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      call get_cpair(kbc, tmp, r1, rvdw, capi, capj)
      cpair = cpair + tmp
   end do
end subroutine get_cpair_dir

subroutine get_dcpair(kbc, vec, rvdw, capi, capj, dgpair, dspair)
   real(wp), intent(in) :: kbc
   real(wp), intent(in) :: vec(3)
   real(wp), intent(in) :: rvdw
   real(wp), intent(in) :: capi
   real(wp), intent(in) :: capj
   real(wp), intent(out) :: dgpair(3)
   real(wp), intent(out) :: dspair(3, 3)

   real(wp) :: r1, arg, dtmp

   dgpair(:) = 0.0_wp
   dspair(:, :) = 0.0_wp

   r1 = norm2(vec)
   ! Capacitance of bond between atom i and j
   arg = -(kbc * (r1 - rvdw) / rvdw)**2
   dtmp = -sqrt(capi * capj) * kbc * exp(arg) / (sqrtpi * rvdw)
   dgpair = dtmp * vec / r1
   dspair = spread(dgpair, 1, 3) * spread(vec, 2, 3)
end subroutine get_dcpair

subroutine get_dcmat_0d(self, mol, dcdr, dcdL)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: dcdr(:, :, :)
   real(wp), intent(out) :: dcdL(:, :, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), rvdw, dG(3), dS(3, 3), capi, capj
   real(wp), allocatable :: dcdiag(:, :)

   allocate(dcdiag(3, mol%nat), source=0.0_wp)
   if (size(dcdr, 3) > mol%nat) dcdr(:, :, mol%nat + 1) = 0.0_wp
   dcdL(:, :, :) = 0.0_wp

   ! Off-diagonal pairs (jat < iat) have a unique owner; reduce only O(N) terms
   !$omp parallel default(none) &
   !$omp& shared(dcdr, dcdL, dcdiag, mol, self) &
   !$omp& private(iat, izp, jat, jzp, vec, rvdw, dG, dS, capi, capj)
   !$omp do schedule(guided) reduction(+:dcdiag, dcdL)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      capi = self%cap(izp)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)

         call get_dcpair(self%kbc, vec, rvdw, capi, capj, dG, dS)

         ! Off-diagonal elements
         dcdr(:, iat, jat) = +dG
         dcdr(:, jat, iat) = -dG
         ! Diagonal elements
         dcdiag(:, iat) = dcdiag(:, iat) - dG
         dcdiag(:, jat) = dcdiag(:, jat) + dG
         dcdL(:, :, iat) = dcdL(:, :, iat) + dS
         dcdL(:, :, jat) = dcdL(:, :, jat) + dS
      end do
   end do
   !$omp end do
   !$omp do schedule(static)
   do iat = 1, mol%nat
      dcdr(:, iat, iat) = dcdiag(:, iat)
   end do
   !$omp end do
   !$omp end parallel

end subroutine get_dcmat_0d

subroutine get_dcmat_3d(self, mol, wsc, dcdr, dcdL)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell), intent(in) :: wsc
   real(wp), intent(out) :: dcdr(:, :, :)
   real(wp), intent(out) :: dcdL(:, :, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: rij(3), vec(3), rvdw, dG(3), dS(3, 3), dG_sum(3), dS_sum(3, 3)
   real(wp) :: capi, capj, wsw
   real(wp), allocatable :: dtrans(:, :), dcdiag(:, :)

   call get_dir_trans(mol, dtrans, self%ncoord%cutoff)

   allocate(dcdiag(3, mol%nat), source=0.0_wp)
   if (size(dcdr, 3) > mol%nat) dcdr(:, :, mol%nat + 1) = 0.0_wp
   dcdL(:, :, :) = 0.0_wp

   ! Off-diagonal pairs (jat < iat) have a unique owner; reduce only O(N) terms
   !$omp parallel default(none) &
   !$omp& shared(dcdr, dcdL, dcdiag, mol, self, dtrans, wsc) &
   !$omp& private(iat, izp, jat, jzp, img, rij, vec, rvdw, dG, dS, &
   !$omp& dG_sum, dS_sum, capi, capj, wsw)
   !$omp do schedule(guided) reduction(+:dcdiag, dcdL)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      capi = self%cap(izp)
      do jat = 1, iat - 1
         jzp = mol%id(jat)
         capj = self%cap(jzp)
         rvdw = self%rvdw(izp, jzp)
         rij = mol%xyz(:, jat) - mol%xyz(:, iat)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         dG_sum(:) = 0.0_wp
         dS_sum(:, :) = 0.0_wp
         do img = 1, wsc%nimg(jat, iat)
            vec = rij + wsc%trans(:, wsc%tridx(img, jat, iat))

            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
            dG_sum = dG_sum + dG
            dS_sum = dS_sum + dS
         end do
         dG_sum = dG_sum * wsw
         dS_sum = dS_sum * wsw

         ! Off-diagonal elements
         dcdr(:, iat, jat) = +dG_sum
         dcdr(:, jat, iat) = -dG_sum
         ! Diagonal elements
         dcdiag(:, iat) = dcdiag(:, iat) - dG_sum
         dcdiag(:, jat) = dcdiag(:, jat) + dG_sum
         dcdL(:, :, jat) = dcdL(:, :, jat) + dS_sum
         dcdL(:, :, iat) = dcdL(:, :, iat) + dS_sum
      end do

      rvdw = self%rvdw(izp, izp)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      dS_sum(:, :) = 0.0_wp
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))

         call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)
         dS_sum = dS_sum + dS
      end do
      ! Positive diagonal elements
      dcdL(:, :, iat) = dcdL(:, :, iat) + dS_sum * wsw
   end do
   !$omp end do
   !$omp do schedule(static)
   do iat = 1, mol%nat
      dcdr(:, iat, iat) = dcdiag(:, iat)
   end do
   !$omp end do
   !$omp end parallel

end subroutine get_dcmat_3d

subroutine get_dcpair_dir(kbc, rij, trans, rvdw, capi, capj, dgpair, dspair)
   real(wp), intent(in) :: rij(3), capi, capj, rvdw, kbc, trans(:, :)
   real(wp), intent(out) :: dgpair(3)
   real(wp), intent(out) :: dspair(3, 3)

   integer :: itr
   real(wp) :: r1, arg, dtmp, dgtmp(3), dstmp(3, 3), vec(3)

   dgpair(:) = 0.0_wp
   dspair(:, :) = 0.0_wp
   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      call get_dcpair(kbc, vec, rvdw, capi, capj, dgtmp, dstmp)
      dgpair(:) = dgpair + dgtmp
      dspair(:, :) = dspair + dstmp
   end do
end subroutine get_dcpair_dir

end module multicharge_model_eeqbc
