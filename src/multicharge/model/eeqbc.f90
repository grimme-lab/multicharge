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

!> Bond capacitor electronegativity equilibration charge model
module multicharge_model_eeqbc

   use mctc_env, only: error_type, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_ncoord, only: new_ncoord, cn_count
   use multicharge_wignerseitz, only: new_wignerseitz_cell, wignerseitz_cell_type
   use multicharge_model_type, only: mchrg_model_type, get_dir_trans
   use multicharge_blas, only: gemv, gemm
   use multicharge_model_cache, only: cache_container, model_cache
   implicit none
   private

   public :: eeqbc_model, new_eeqbc_model

   !> Cache for the EEQ-BC charge model
   type, extends(model_cache) :: eeqbc_cache
      !> Local charges
      real(wp), allocatable :: qloc(:)
      !> Local charge dr derivative
      real(wp), allocatable :: dqlocdr(:, :, :)
      !> Local charge dL derivative
      real(wp), allocatable :: dqlocdL(:, :, :)
      !> Full Maxwell capacitance matrix
      real(wp), allocatable :: cmat(:, :)
      !> Derivative of Maxwell capacitance matrix w.r.t positions
      real(wp), allocatable :: dcdr(:, :, :)
      !> Derivative of Maxwell capacitance matrix w.r.t lattice vectors
      real(wp), allocatable :: dcdL(:, :, :)
      !> Store tmp array from xvec calculation for reuse
      real(wp), allocatable :: xtmp(:)
   end type eeqbc_cache

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
      !> Calculate Coulomb matrix
      procedure :: get_coulomb_matrix
      !> Calculate derivatives of Coulomb matrix
      procedure :: get_coulomb_derivs
      !> Calculate right-hand side (electronegativity vector)
      procedure :: get_xvec
      !> Calculate derivatives of EN vector
      procedure :: get_xvec_derivs
      !> Calculate Coulomb matrix
      procedure :: get_amat_0d
      !> Calculate Coulomb matrix periodic
      procedure :: get_amat_3d
      !> Calculate Coulomb matrix derivative
      procedure :: get_damat_0d
      !> Calculate Coulomb matrix derivative (periodic)
      procedure :: get_damat_3d
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
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Default exponent of distance/CN normalization
   real(wp), parameter :: default_norm_exp = 1.0_wp

   !> Default exponent of error function in bond capacitance
   real(wp), parameter :: default_kbc = 0.65_wp
contains

   subroutine new_eeqbc_model(self, mol, error, chi, rad, &
      & eta, kcnchi, kqchi, kqeta, kcnrad, cap, avg_cn, &
      & kbc, cutoff, cn_exp, rcov, en, cn_max, norm_exp, rvdw)
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
      !> Van-der-Waals radii
      real(wp), intent(in), optional :: rvdw(:, :)

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
      ! Electronegativity weighted coordination number for local charge
      call new_ncoord(self%ncoord_en, mol, cn_count%erf_en, error, &
         & cutoff=cutoff, kcn=cn_exp, rcov=rcov, en=en, cut=cn_max, &
         & norm_exp=self%norm_exp)

   end subroutine new_eeqbc_model

   subroutine update(self, mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(in) :: cn(:)
      real(wp), intent(in), optional :: qloc(:)
      real(wp), intent(in), optional :: dcndr(:, :, :)
      real(wp), intent(in), optional :: dcndL(:, :, :)
      real(wp), intent(in), optional :: dqlocdr(:, :, :)
      real(wp), intent(in), optional :: dqlocdL(:, :, :)

      logical :: grad

      type(eeqbc_cache), pointer :: ptr

      call taint(cache, ptr)

      grad = present(dcndr) .and. present(dcndL) .and. present(dqlocdr) .and. present(dqlocdL)

      ! Refer CN and local charge arrays in cache
      ptr%cn = cn
      if (present(qloc)) then
         ptr%qloc = qloc
      else
         error stop "qloc required for eeqbc"
      end if

      if (grad) then
         ptr%dcndr = dcndr
         ptr%dcndL = dcndL
         ptr%dqlocdr = dqlocdr
         ptr%dqlocdL = dqlocdL
      end if

      ! Allocate (for get_xvec and xvec_derivs)
      if (.not. allocated(ptr%xtmp)) then
         allocate (ptr%xtmp(mol%nat + 1))
      end if

      ! Allocate cmat
      if (.not. allocated(ptr%cmat)) then
         allocate (ptr%cmat(mol%nat + 1, mol%nat + 1))
      end if

      if (any(mol%periodic)) then
         ! Create WSC
         call new_wignerseitz_cell(ptr%wsc, mol)

         ! Get full cmat sum over all WSC images (for get_xvec and xvec_derivs)
         call self%get_cmat_3d(mol, ptr%wsc, ptr%cmat)
         if (grad) then
            if (.not. allocated(ptr%dcdr) .and. .not. allocated(ptr%dcdL)) then
               allocate (ptr%dcdr(3, mol%nat, mol%nat + 1), ptr%dcdL(3, 3, mol%nat + 1))
            end if
            call self%get_dcmat_3d(mol, ptr%wsc, ptr%dcdr, ptr%dcdL)
         end if
      else
         call self%get_cmat_0d(mol, ptr%cmat)

         ! cmat gradients
         if (grad) then
            if (.not. allocated(ptr%dcdr) .and. .not. allocated(ptr%dcdL)) then
               allocate (ptr%dcdr(3, mol%nat, mol%nat + 1), ptr%dcdL(3, 3, mol%nat + 1))
            end if
            call self%get_dcmat_0d(mol, ptr%dcdr, ptr%dcdL)
         end if
      end if

   end subroutine update

   subroutine get_xvec(self, mol, cache, xvec)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out) :: xvec(:)

      type(eeqbc_cache), pointer :: ptr

      integer :: iat, izp, img
      real(wp) :: ctmp, vec(3), rvdw, capi, wsw
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private array for reduction
      real(wp), allocatable :: xvec_local(:)

      call view(cache, ptr)

      xvec(:) = 0.0_wp
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, ptr, xvec) &
      !$omp private(iat, izp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ptr%xtmp(iat) = -self%chi(izp) + self%kcnchi(izp)*ptr%cn(iat) &
            & + self%kqchi(izp)*ptr%qloc(iat)
      end do
      ptr%xtmp(mol%nat + 1) = mol%charge

      call gemv(ptr%cmat, ptr%xtmp, xvec)

      if (any(mol%periodic)) then
         call get_dir_trans(mol%lattice, dtrans)
         !$omp parallel default(none) &
         !$omp shared(mol, self, ptr, xvec, dtrans) private(iat, izp, img, wsw) &
         !$omp private(capi, vec, rvdw, ctmp, xvec_local)
         allocate (xvec_local, mold=xvec)
         xvec_local(:) = 0.0_wp
         !$omp do schedule(runtime)
         do iat = 1, mol%nat
            izp = mol%id(iat)
            capi = self%cap(izp)
            ! eliminate self-interaction (quasi off-diagonal)
            rvdw = self%rvdw(iat, iat)
            wsw = 1.0_wp/real(ptr%wsc%nimg(iat, iat), wp)
            do img = 1, ptr%wsc%nimg(iat, iat)
               vec = ptr%wsc%trans(:, ptr%wsc%tridx(img, iat, iat))

               call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
               xvec_local(iat) = xvec_local(iat) - wsw*ctmp*ptr%xtmp(iat)
            end do
         end do
         !$omp end do
         !$omp critical (get_xvec_)
         xvec(:) = xvec(:) + xvec_local(:)
         !$omp end critical (get_xvec_)
         deallocate (xvec_local)
         !$omp end parallel
      end if
   end subroutine get_xvec

   subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out), contiguous :: dxdr(:, :, :)
      real(wp), intent(out), contiguous :: dxdL(:, :, :)

      type(eeqbc_cache), pointer :: ptr

      integer :: iat, izp, jat, jzp, img
      real(wp) :: capi, capj, wsw, vec(3), ctmp, rvdw, dG(3), dS(3, 3)
      real(wp), allocatable :: dtmpdr(:, :, :), dtmpdL(:, :, :)
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dxdr_local(:, :, :), dxdL_local(:, :, :), dtmpdr_local(:, :, :), dtmpdL_local(:, :, :)

      call view(cache, ptr)
      allocate (dtmpdr(3, mol%nat, mol%nat + 1), dtmpdL(3, 3, mol%nat + 1))

      dxdr(:, :, :) = 0.0_wp
      dxdL(:, :, :) = 0.0_wp
      dtmpdr(:, :, :) = 0.0_wp
      dtmpdL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(mol, self, ptr, dtmpdr, dtmpdL) &
      !$omp private(iat, izp, dtmpdr_local, dtmpdL_local)
      allocate (dtmpdr_local, source=dtmpdr)
      allocate (dtmpdL_local, source=dtmpdL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! CN and effective charge derivative
         dtmpdr_local(:, :, iat) = self%kcnchi(izp)*ptr%dcndr(:, :, iat) + dtmpdr_local(:, :, iat)
         dtmpdL_local(:, :, iat) = self%kcnchi(izp)*ptr%dcndL(:, :, iat) + dtmpdL_local(:, :, iat)
         dtmpdr_local(:, :, iat) = self%kqchi(izp)*ptr%dqlocdr(:, :, iat) + dtmpdr_local(:, :, iat)
         dtmpdL_local(:, :, iat) = self%kqchi(izp)*ptr%dqlocdL(:, :, iat) + dtmpdL_local(:, :, iat)
      end do
      !$omp end do
      !$omp critical (get_xvec_derivs_)
      dtmpdr(:, :, :) = dtmpdr(:, :, :) + dtmpdr_local(:, :, :)
      dtmpdL(:, :, :) = dtmpdL(:, :, :) + dtmpdL_local(:, :, :)
      !$omp end critical (get_xvec_derivs_)
      deallocate (dtmpdL_local, dtmpdr_local)
      !$omp end parallel

      call gemm(dtmpdr, ptr%cmat, dxdr)
      call gemm(dtmpdL, ptr%cmat, dxdL)

      if (any(mol%periodic)) then
         call get_dir_trans(mol%lattice, dtrans)
         !$omp parallel default(none) &
         !$omp shared(mol, self, ptr, dxdr, dxdL, dtrans) &
         !$omp private(iat, izp, jat, jzp, img, wsw) &
         !$omp private(capi, capj, vec, rvdw, ctmp, dG, dS) &
         !$omp private(dxdr_local, dxdL_local)
         allocate (dxdr_local, mold=dxdr)
         allocate (dxdL_local, mold=dxdL)
         dxdr_local(:, :, :) = 0.0_wp
         dxdL_local(:, :, :) = 0.0_wp
         !$omp do schedule(runtime)
         do iat = 1, mol%nat
            izp = mol%id(iat)
            capi = self%cap(izp)
            do jat = 1, mol%nat
               ! Diagonal elements
               dxdr_local(:, iat, iat) = dxdr_local(:, iat, iat) + ptr%xtmp(jat)*ptr%dcdr(:, iat, jat)
               ! Derivative of capacitance matrix
               dxdr_local(:, iat, jat) = (ptr%xtmp(iat) - ptr%xtmp(jat))*ptr%dcdr(:, iat, jat) &
               & + dxdr_local(:, iat, jat)
               jzp = mol%id(jat)
               capj = self%cap(jzp)
               rvdw = self%rvdw(iat, jat)
               wsw = 1.0_wp/real(ptr%wsc%nimg(iat, jat), wp)
               do img = 1, ptr%wsc%nimg(iat, jat)
                  vec = mol%xyz(:, jat) - mol%xyz(:, iat) + ptr%wsc%trans(:, ptr%wsc%tridx(img, jat, iat))
                  call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
                  dxdL_local(:, :, iat) = dxdL_local(:, :, iat) + wsw*dS*ptr%xtmp(jat)
               end do
            end do
            dxdL_local(:, :, iat) = dxdL_local(:, :, iat) + ptr%xtmp(iat)*ptr%dcdL(:, :, iat)

            ! Capacitance terms for i = j, T != 0
            rvdw = self%rvdw(iat, iat)
            wsw = 1.0_wp/real(ptr%wsc%nimg(iat, iat), wp)
            do img = 1, ptr%wsc%nimg(iat, iat)
               vec = ptr%wsc%trans(:, ptr%wsc%tridx(img, iat, iat))

               call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, ctmp)
               ctmp = ctmp*wsw
               ! EN derivative
               dxdr_local(:, :, iat) = dxdr_local(:, :, iat) - ctmp*self%kcnchi(izp)*ptr%dcndr(:, :, iat)
               dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - ctmp*self%kcnchi(izp)*ptr%dcndL(:, :, iat)
               dxdr_local(:, :, iat) = dxdr_local(:, :, iat) - ctmp*self%kqchi(izp)*ptr%dqlocdr(:, :, iat)
               dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - ctmp*self%kqchi(izp)*ptr%dqlocdL(:, :, iat)
            end do
         end do
         !$omp end do
         !$omp critical (get_xvec_derivs_update)
         dxdr(:, :, :) = dxdr(:, :, :) + dxdr_local(:, :, :)
         dxdL(:, :, :) = dxdL(:, :, :) + dxdL_local(:, :, :)
         !$omp end critical (get_xvec_derivs_update)
         deallocate (dxdL_local, dxdr_local)
         !$omp end parallel
      else
         !$omp parallel default(none) &
         !$omp shared(mol, self, ptr, dxdr, dxdL) &
         !$omp private(iat, izp, jat, jzp, vec, dxdr_local, dxdL_local)
         allocate (dxdr_local, mold=dxdr)
         allocate (dxdL_local, mold=dxdL)
         dxdr_local(:, :, :) = 0.0_wp
         dxdL_local(:, :, :) = 0.0_wp
         !$omp do schedule(runtime)
         do iat = 1, mol%nat
            do jat = 1, mol%nat
               vec = mol%xyz(:, iat) - mol%xyz(:, jat)
               ! Derivative of capacitance matrix
               dxdr_local(:, iat, jat) = (ptr%xtmp(iat) - ptr%xtmp(jat))*ptr%dcdr(:, iat, jat) + dxdr_local(:, iat, jat)

               if (iat .eq. jat) cycle
               dxdL_local(:, :, iat) = dxdL_local(:, :, iat) + ptr%xtmp(jat)*spread(ptr%dcdr(:, iat, jat), 1, 3)*spread(vec, 2, 3)
               ! dxdL_local(:, :, iat) = dxdL_local(:, :, iat) - ptr%xtmp(iat)*spread(ptr%dcdr(:, iat, jat), 1, 3)*spread(vec, 2, 3) ! A
            end do
            do jat = 1, iat - 1
               ! Diagonal elements
               dxdr_local(:, iat, iat) = dxdr_local(:, iat, iat) + ptr%xtmp(jat)*ptr%dcdr(:, iat, jat)
               dxdr_local(:, jat, jat) = dxdr_local(:, jat, jat) + ptr%xtmp(iat)*ptr%dcdr(:, jat, iat)
            end do
            dxdr_local(:, iat, iat) = dxdr_local(:, iat, iat) + ptr%xtmp(iat)*ptr%dcdr(:, iat, iat)
            dxdL_local(:, :, iat) = dxdL_local(:, :, iat) + ptr%xtmp(iat)*ptr%dcdL(:, :, iat) ! remove if using A
         end do
         !$omp end do
         !$omp critical (get_xvec_derivs_)
         dxdr(:, :, :) = dxdr(:, :, :) + dxdr_local(:, :, :)
         dxdL(:, :, :) = dxdL(:, :, :) + dxdL_local(:, :, :)
         !$omp end critical (get_xvec_derivs_)
         deallocate (dxdL_local, dxdr_local)
         !$omp end parallel
      end if

   end subroutine get_xvec_derivs

   subroutine get_coulomb_matrix(self, mol, cache, amat)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out) :: amat(:, :)

      type(eeqbc_cache), pointer :: ptr
      call view(cache, ptr)

      if (any(mol%periodic)) then
         call self%get_amat_3d(mol, ptr%wsc, ptr%cn, ptr%qloc, ptr%cmat, amat)
      else
         call self%get_amat_0d(mol, ptr%cn, ptr%qloc, ptr%cmat, amat)
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

      ! Thread-private array for reduction
      real(wp), allocatable :: amat_local(:, :)

      amat(:, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(amat, mol, self, cn, qloc, cmat) &
      !$omp private(iat, izp, jat, jzp, gam2, vec, r2, tmp) &
      !$omp private(norm_cn, radi, radj, amat_local)
      allocate (amat_local, source=amat)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp/self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp)*(1.0_wp - self%kcnrad*cn(iat)*norm_cn)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! Effective charge width of j
            norm_cn = cn(jat)/self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp)*(1.0_wp - self%kcnrad*norm_cn)
            ! Coulomb interaction of Gaussian charges
            gam2 = 1.0_wp/(radi**2 + radj**2)
            tmp = erf(sqrt(r2*gam2))/sqrt(r2)*cmat(jat, iat)
            amat_local(jat, iat) = tmp
            amat_local(iat, jat) = tmp
         end do
         ! Effective hardness
         tmp = self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi
         amat_local(iat, iat) = amat_local(iat, iat) + tmp*cmat(iat, iat) + 1.0_wp
      end do
      !$omp end do
      !$omp critical (get_amat_0d_)
      amat(:, :) = amat(:, :) + amat_local(:, :)
      !$omp end critical (get_amat_0d_)
      deallocate (amat_local)
      !$omp end parallel

      amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
      amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
      amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

   end subroutine get_amat_0d

   subroutine get_amat_3d(self, mol, wsc, cn, qloc, cmat, amat)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(wignerseitz_cell_type), intent(in) :: wsc
      real(wp), intent(in) :: cn(:), qloc(:), cmat(:, :)
      real(wp), intent(out) :: amat(:, :)

      integer :: iat, jat, izp, jzp, img
      real(wp) :: vec(3), r1, gam, dtmp, ctmp, capi, capj, radi, radj, norm_cn, rvdw, wsw
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private array for reduction
      real(wp), allocatable :: amat_local(:, :)

      call get_dir_trans(mol%lattice, dtrans)

      amat(:, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(amat, cmat, mol, cn, qloc, self, wsc, dtrans)  &
      !$omp private(iat, izp, jat, jzp, gam, vec, dtmp, ctmp, norm_cn) &
      !$omp private(radi, radj, capi, capj, rvdw, r1, wsw, amat_local)
      allocate (amat_local, source=amat)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp/self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp)*(1.0_wp - self%kcnrad*cn(iat)*norm_cn)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            ! vdw distance in Angstrom (approximate factor 2)
            rvdw = self%rvdw(iat, jat)
            ! Effective charge width of j
            norm_cn = cn(jat)/self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp)*(1.0_wp - self%kcnrad*norm_cn)
            capj = self%cap(jzp)
            ! Coulomb interaction of Gaussian charges
            gam = 1.0_wp/sqrt(radi**2 + radj**2)
            wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))
               call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capj, dtmp)
               amat_local(jat, iat) = amat_local(jat, iat) + dtmp*wsw
               amat_local(iat, jat) = amat_local(iat, jat) + dtmp*wsw
            end do
         end do

         ! diagonal Coulomb interaction terms
         gam = 1.0_wp/sqrt(2.0_wp*radi**2)
         rvdw = self%rvdw(iat, iat)
         wsw = 1.0_wp/real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))
            call get_amat_dir_3d(vec, gam, dtrans, self%kbc, rvdw, capi, capi, dtmp)
            amat_local(iat, iat) = amat_local(iat, iat) + dtmp*wsw
         end do

         ! Effective hardness
         dtmp = self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi
         amat_local(iat, iat) = amat_local(iat, iat) + cmat(iat, iat)*dtmp + 1.0_wp
      end do
      !$omp end do
      !$omp critical (get_amat_3d_)
      amat(:, :) = amat(:, :) + amat_local(:, :)
      !$omp end critical (get_amat_3d_)
      deallocate (amat_local)
      !$omp end parallel

      amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
      amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
      amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

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
         tmp = -ctmp*erf(gam*r1)/r1
         amat = amat + tmp
      end do

   end subroutine get_amat_dir_3d

   subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: qvec(:)
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

      type(eeqbc_cache), pointer :: ptr
      call view(cache, ptr)

      if (any(mol%periodic)) then
         call self%get_damat_3d(mol, ptr%wsc, ptr%cn, &
         & ptr%qloc, qvec, ptr%dcndr, ptr%dcndL, ptr%dqlocdr, &
         & ptr%dqlocdL, ptr%cmat, ptr%dcdr, ptr%dcdL, dadr, dadL, atrace)

      else
         call self%get_damat_0d(mol, ptr%cn, &
         & ptr%qloc, qvec, ptr%dcndr, ptr%dcndL, ptr%dqlocdr, &
         & ptr%dqlocdL, ptr%cmat, ptr%dcdr, ptr%dcdL, dadr, dadL, atrace)
      end if
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

      allocate (dgamdr(3, mol%nat))

      atrace(:, :) = 0.0_wp
      dadr(:, :, :) = 0.0_wp
      dadL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(atrace, dadr, dadL, mol, self, cn, qloc, qvec) &
      !$omp shared(cmat, dcdr, dcdL, dcndr, dcndL, dqlocdr, dqlocdL) &
      !$omp private(iat, izp, jat, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
      !$omp private(radi, radj, dradi, dradj, dgamdr, dgamdL, dG, dS) &
      !$omp private(atrace_local, dadr_local, dadL_local)
      allocate (atrace_local, source=atrace)
      allocate (dadr_local, source=dadr)
      allocate (dadL_local, source=dadL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp/self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp)*(1.0_wp - self%kcnrad*cn(iat)*norm_cn)
         dradi = -self%rad(izp)*self%kcnrad*norm_cn
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            ! Effective charge width of j
            norm_cn = 1.0_wp/self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp)*(1.0_wp - self%kcnrad*cn(jat)*norm_cn)
            dradj = -self%rad(jzp)*self%kcnrad*norm_cn

            ! Coulomb interaction of Gaussian charges
            gam = 1.0_wp/sqrt(radi**2 + radj**2)
            dgamdr(:, :) = -(radi*dradi*dcndr(:, :, iat) + radj*dradj*dcndr(:, :, jat)) &
                          & *gam**3.0_wp
            dgamdL(:, :) = -(radi*dradi*dcndL(:, :, iat) + radj*dradj*dcndL(:, :, jat)) &
                          & *gam**3.0_wp

            ! Explicit derivative
            arg = gam*gam*r2
            dtmp = 2.0_wp*gam*exp(-arg)/(sqrtpi*r2) &
               & - erf(sqrt(arg))/(r2*sqrt(r2))
            dG(:) = dtmp*vec
            dS(:, :) = spread(dG, 1, 3)*spread(vec, 2, 3)
            atrace_local(:, iat) = -dG*qvec(jat)*cmat(jat, iat) + atrace_local(:, iat)
            atrace_local(:, jat) = +dG*qvec(iat)*cmat(iat, jat) + atrace_local(:, jat)
            dadr_local(:, iat, jat) = -dG*qvec(iat)*cmat(iat, jat) + dadr_local(:, iat, jat)
            dadr_local(:, jat, iat) = +dG*qvec(jat)*cmat(jat, iat) + dadr_local(:, jat, iat)
            dadL_local(:, :, iat) = +dS*qvec(jat)*cmat(jat, iat) + dadL_local(:, :, iat)
            dadL_local(:, :, jat) = +dS*qvec(iat)*cmat(iat, jat) + dadL_local(:, :, jat)

            ! Effective charge width derivative
            dtmp = 2.0_wp*exp(-arg)/(sqrtpi)
            atrace_local(:, iat) = -dtmp*qvec(jat)*dgamdr(:, jat)*cmat(jat, iat) + atrace_local(:, iat)
            atrace_local(:, jat) = -dtmp*qvec(iat)*dgamdr(:, iat)*cmat(iat, jat) + atrace_local(:, jat)
            dadr_local(:, iat, jat) = +dtmp*qvec(iat)*dgamdr(:, iat)*cmat(iat, jat) + dadr_local(:, iat, jat)
            dadr_local(:, jat, iat) = +dtmp*qvec(jat)*dgamdr(:, jat)*cmat(jat, iat) + dadr_local(:, jat, iat)
            dadL_local(:, :, iat) = +dtmp*qvec(jat)*dgamdL(:, :)*cmat(jat, iat) + dadL_local(:, :, iat)
            dadL_local(:, :, jat) = +dtmp*qvec(iat)*dgamdL(:, :)*cmat(iat, jat) + dadL_local(:, :, jat)

            ! Capacitance derivative off-diagonal
            dtmp = erf(sqrt(r2)*gam)/(sqrt(r2))
            ! potentially switch indices for dcdr
            atrace_local(:, iat) = -dtmp*qvec(jat)*dcdr(:, jat, iat) + atrace_local(:, iat)
            atrace_local(:, jat) = -dtmp*qvec(iat)*dcdr(:, iat, jat) + atrace_local(:, jat)
            dadr_local(:, iat, jat) = +dtmp*qvec(iat)*dcdr(:, iat, jat) + dadr_local(:, iat, jat)
            dadr_local(:, jat, iat) = +dtmp*qvec(jat)*dcdr(:, jat, iat) + dadr_local(:, jat, iat)
            dadL_local(:, :, iat) = dadL_local(:, :, iat) - dtmp*qvec(jat)*spread(dcdr(:, iat, jat), 2, 3)*spread(vec, 1, 3)
            dadL_local(:, :, jat) = dadL_local(:, :, jat) - dtmp*qvec(iat)*spread(dcdr(:, iat, jat), 2, 3)*spread(vec, 1, 3)

            ! Capacitance derivative diagonal
            dtmp = (self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi)*qvec(iat)
            dadr_local(:, jat, iat) = -dtmp*dcdr(:, jat, iat) + dadr_local(:, jat, iat)

            dtmp = (self%eta(jzp) + self%kqeta(jzp)*qloc(jat) + sqrt2pi/radj)*qvec(jat)
            dadr_local(:, iat, jat) = -dtmp*dcdr(:, iat, jat) + dadr_local(:, iat, jat)
         end do

         ! Hardness derivative
         dtmp = self%kqeta(izp)*qvec(iat)*cmat(iat, iat)
         dadr_local(:, :, iat) = +dtmp*dqlocdr(:, :, iat) + dadr_local(:, :, iat)
         dadL_local(:, :, iat) = +dtmp*dqlocdL(:, :, iat) + dadL_local(:, :, iat)

         ! Effective charge width derivative
         dtmp = -sqrt2pi*dradi/(radi**2)*qvec(iat)*cmat(iat, iat)
         dadr_local(:, :, iat) = +dtmp*dcndr(:, :, iat) + dadr_local(:, :, iat)
         dadL_local(:, :, iat) = +dtmp*dcndL(:, :, iat) + dadL_local(:, :, iat)

         ! Capacitance derivative
         dtmp = (self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi)*qvec(iat)
         dadr_local(:, iat, iat) = +dtmp*dcdr(:, iat, iat) + dadr_local(:, iat, iat)
         dadL_local(:, :, iat) = +dtmp*dcdL(:, :, iat) + dadL_local(:, :, iat)

      end do
      !$omp end do
      !$omp critical (get_damat_0d_)
      atrace(:, :) = atrace(:, :) + atrace_local(:, :)
      dadr(:, :, :) = dadr(:, :, :) + dadr_local(:, :, :)
      dadL(:, :, :) = dadL(:, :, :) + dadL_local(:, :, :)
      !$omp end critical (get_damat_0d_)
      deallocate (dadL_local, dadr_local, atrace_local)
      !$omp end parallel

   end subroutine get_damat_0d

   subroutine get_damat_3d(self, mol, wsc, cn, qloc, qvec, dcndr, dcndL, dqlocdr, &
      & dqlocdL, cmat, dcdr, dcdL, dadr, dadL, atrace)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(wignerseitz_cell_type), intent(in) :: wsc
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
      real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn, rvdw, wsw, dgam
      real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3)
      real(wp) :: dgamdL(3, 3), capi, capj
      real(wp), allocatable :: dgamdr(:, :), dtrans(:, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: atrace_local(:, :)
      real(wp), allocatable :: dadr_local(:, :, :), dadL_local(:, :, :)

      call get_dir_trans(mol%lattice, dtrans)

      allocate (dgamdr(3, mol%nat))

      atrace(:, :) = 0.0_wp
      dadr(:, :, :) = 0.0_wp
      dadL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(self, mol, cn, qloc, qvec, wsc, dadr, dadL, atrace) &
      !$omp shared (cmat, dcdr, dcdL, dcndr, dcndL, dqlocdr, dqlocdL, dtrans) &
      !$omp private(iat, izp, jat, jzp, img, gam, vec, r2, dtmp, norm_cn, arg, rvdw) &
      !$omp private(radi, radj, dradi, dradj, capi, capj, dgamdr, dgamdL, dG, dS, wsw) &
      !$omp private(dgam, dadr_local, dadL_local, atrace_local)
      allocate (atrace_local, source=atrace)
      allocate (dadr_local, source=dadr)
      allocate (dadL_local, source=dadL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! Effective charge width of i
         norm_cn = 1.0_wp/self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp)*(1.0_wp - self%kcnrad*cn(iat)*norm_cn)
         dradi = -self%rad(izp)*self%kcnrad*norm_cn
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(iat, jat)

            ! Effective charge width of j
            norm_cn = 1.0_wp/self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp)*(1.0_wp - self%kcnrad*cn(jat)*norm_cn)
            dradj = -self%rad(jzp)*self%kcnrad*norm_cn

            ! Coulomb interaction of Gaussian charges
            gam = 1.0_wp/sqrt(radi**2 + radj**2)
            dgamdr(:, :) = -(radi*dradi*dcndr(:, :, iat) + radj*dradj*dcndr(:, :, jat)) &
                          & *gam**3.0_wp
            dgamdL(:, :) = -(radi*dradi*dcndL(:, :, iat) + radj*dradj*dcndL(:, :, jat)) &
                          & *gam**3.0_wp

            wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))

               call get_damat_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS, dgam)
               dG = dG*wsw
               dS = dS*wsw
               dgam = dgam*wsw

               ! Explicit derivative
               atrace_local(:, iat) = -dG*qvec(jat) + atrace_local(:, iat)
               atrace_local(:, jat) = +dG*qvec(iat) + atrace_local(:, jat)
               dadr_local(:, iat, jat) = -dG*qvec(iat) + dadr_local(:, iat, jat)
               dadr_local(:, jat, iat) = +dG*qvec(jat) + dadr_local(:, jat, iat)
               dadL_local(:, :, jat) = +dS*qvec(iat) + dadL_local(:, :, jat)
               dadL_local(:, :, iat) = +dS*qvec(jat) + dadL_local(:, :, iat)

               ! Effective charge width derivative
               atrace_local(:, iat) = +dgam*qvec(jat)*dgamdr(:, jat) + atrace_local(:, iat)
               atrace_local(:, jat) = +dgam*qvec(iat)*dgamdr(:, iat) + atrace_local(:, jat)
               dadr_local(:, iat, jat) = -dgam*qvec(iat)*dgamdr(:, iat) + dadr_local(:, iat, jat)
               dadr_local(:, jat, iat) = -dgam*qvec(jat)*dgamdr(:, jat) + dadr_local(:, jat, iat)
               dadL_local(:, :, iat) = -dgam*qvec(jat)*dgamdL(:, :) + dadL_local(:, :, iat)
               dadL_local(:, :, jat) = -dgam*qvec(iat)*dgamdL(:, :) + dadL_local(:, :, jat)

               call get_damat_dc_dir(vec, dtrans, capi, capj, rvdw, self%kbc, gam, dG, dS)
               dG = dG*wsw
               dS = dS*wsw

               ! Capacitance derivative off-diagonal
               atrace_local(:, iat) = -qvec(jat)*dG(:) + atrace_local(:, iat)
               atrace_local(:, jat) = +qvec(iat)*dG(:) + atrace_local(:, jat)
               dadr_local(:, jat, iat) = +qvec(jat)*dG(:) + dadr_local(:, jat, iat)
               dadr_local(:, iat, jat) = -qvec(iat)*dG(:) + dadr_local(:, iat, jat)
               dadL_local(:, :, jat) = +qvec(iat)*dS(:, :) + dadL_local(:, :, jat)
               dadL_local(:, :, iat) = +qvec(jat)*dS(:, :) + dadL_local(:, :, iat)

               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)
               dG = dG*wsw

               ! Capacitance derivative diagonal
               dtmp = (self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi)*qvec(iat)
               dadr_local(:, jat, iat) = -dtmp*dG(:) + dadr_local(:, jat, iat)
               dtmp = (self%eta(jzp) + self%kqeta(jzp)*qloc(jat) + sqrt2pi/radj)*qvec(jat)
               dadr_local(:, iat, jat) = +dtmp*dG(:) + dadr_local(:, iat, jat)
            end do
         end do

         ! diagonal explicit, charge width, and capacitance derivative terms
         gam = 1.0_wp/sqrt(2.0_wp*radi**2)
         dtmp = -sqrt2pi*dradi/(radi**2)*qvec(iat)
         rvdw = self%rvdw(iat, iat)
         wsw = 1.0_wp/real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))
            call get_damat_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS, dgam)
            dgam = dgam*wsw

            ! Explicit derivative
            dadL_local(:, :, iat) = +dS*wsw*qvec(iat) + dadL_local(:, :, iat)

            ! Effective charge width derivative
            atrace_local(:, iat) = +dtmp*dcndr(:, iat, iat)*dgam + atrace_local(:, iat)
            dadr_local(:, iat, iat) = -dtmp*dcndr(:, iat, iat)*dgam + dadr_local(:, iat, iat)
            dadL_local(:, :, iat) = -dtmp*dcndL(:, :, iat)*dgam + dadL_local(:, :, iat)

            ! Capacitance derivative
            call get_damat_dc_dir(vec, dtrans, capi, capi, rvdw, self%kbc, gam, dG, dS)
            dadL_local(:, :, iat) = +qvec(iat)*dS*wsw + dadL_local(:, :, iat)
         end do

         ! Hardness derivative
         dtmp = self%kqeta(izp)*qvec(iat)*cmat(iat, iat)
         dadr_local(:, :, iat) = +dtmp*dqlocdr(:, :, iat) + dadr_local(:, :, iat)
         dadL_local(:, :, iat) = +dtmp*dqlocdL(:, :, iat) + dadL_local(:, :, iat)

         ! Effective charge width derivative
         dtmp = -sqrt2pi*dradi/(radi**2)*qvec(iat)*cmat(iat, iat)
         dadr_local(:, :, iat) = +dtmp*dcndr(:, :, iat) + dadr_local(:, :, iat)
         dadL_local(:, :, iat) = +dtmp*dcndL(:, :, iat) + dadL_local(:, :, iat)

         dtmp = (self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi)*qvec(iat)
         dadr_local(:, iat, iat) = +dtmp*dcdr(:, iat, iat) + dadr_local(:, iat, iat)
         dadL_local(:, :, iat) = +dtmp*dcdL(:, :, iat) + dadL_local(:, :, iat)

      end do
      !$omp end do
      !$omp critical (get_damat_3d_)
      atrace(:, :) = atrace(:, :) + atrace_local(:, :)
      dadr(:, :, :) = dadr(:, :, :) + dadr_local(:, :, :)
      dadL(:, :, :) = dadL(:, :, :) + dadL_local(:, :, :)
      !$omp end critical (get_damat_3d_)
      deallocate (dadL_local, dadr_local, atrace_local)
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

      gam2 = gam*gam

      do itr = 1, size(trans, 2)
         vec(:) = rij(:) + trans(:, itr)
         r1 = norm2(vec)
         if (r1 < eps) cycle
         r2 = r1*r1
         call get_cpair(kbc, cmat, r1, rvdw, capi, capj)
         gtmp = 2.0_wp*gam*exp(-r2*gam2)/(sqrtpi*r2) - erf(r1*gam)/(r2*r1)
         dG(:) = dG - cmat*gtmp*vec
         dS(:, :) = dS - cmat*gtmp*spread(vec, 1, 3)*spread(vec, 2, 3)
         dgam = dgam + cmat*2.0_wp*exp(-gam2*r2)/sqrtpi
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
         tmp = erf(gam*r1)/r1
         dG(:) = dG(:) + tmp*gtmp
         dS(:, :) = dS(:, :) + tmp*stmp
      end do

   end subroutine get_damat_dc_dir

   subroutine get_cmat_0d(self, mol, cmat)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(out) :: cmat(:, :)

      integer :: iat, jat, izp, jzp
      real(wp) :: vec(3), rvdw, tmp, capi, capj, r1

      ! Thread-private array for reduction
      real(wp), allocatable :: cmat_local(:, :)

      cmat(:, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(cmat, mol, self) &
      !$omp private(iat, izp, jat, jzp) &
      !$omp private(vec, r1, rvdw, tmp, capi, capj, cmat_local)
      allocate (cmat_local, source=cmat)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r1 = norm2(vec)
            rvdw = self%rvdw(iat, jat)
            capj = self%cap(jzp)

            call get_cpair(self%kbc, tmp, r1, rvdw, capi, capj)

            ! Off-diagonal elements
            cmat_local(jat, iat) = -tmp
            cmat_local(iat, jat) = -tmp
            ! Diagonal elements
            cmat_local(iat, iat) = cmat_local(iat, iat) + tmp
            cmat_local(jat, jat) = cmat_local(jat, jat) + tmp
         end do
      end do
      !$omp end do
      !$omp critical (get_cmat_0d_)
      cmat(:, :) = cmat(:, :) + cmat_local(:, :)
      !$omp end critical (get_cmat_0d_)
      deallocate (cmat_local)
      !$omp end parallel

      cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp

   end subroutine get_cmat_0d

   subroutine get_cmat_3d(self, mol, wsc, cmat)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(wignerseitz_cell_type), intent(in) :: wsc
      real(wp), intent(out) :: cmat(:, :)

      integer :: iat, jat, izp, jzp, img
      real(wp) :: vec(3), rvdw, tmp, capi, capj, wsw
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private array for reduction
      real(wp), allocatable :: cmat_local(:, :)

      call get_dir_trans(mol%lattice, dtrans)

      cmat(:, :) = 0.0_wp!1.0_wp/real(25, wp)

      !$omp parallel default(none) &
      !$omp shared(cmat, mol, self, wsc, dtrans) &
      !$omp private(iat, izp, jat, jzp, img) &
      !$omp private(vec, rvdw, tmp, capi, capj, wsw, cmat_local)
      allocate (cmat_local, source=cmat)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            rvdw = self%rvdw(iat, jat)
            capj = self%cap(jzp)
            wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))

               call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, tmp)

               ! Off-diagonal elements
               cmat_local(jat, iat) = cmat_local(jat, iat) - tmp*wsw
               cmat_local(iat, jat) = cmat_local(iat, jat) - tmp*wsw
               ! Diagonal elements
               !$omp atomic
               cmat_local(iat, iat) = cmat_local(iat, iat) + tmp*wsw
               !$omp atomic
               cmat_local(jat, jat) = cmat_local(jat, jat) + tmp*wsw
            end do
         end do

         ! diagonal capacitance (interaction with images)
         rvdw = self%rvdw(iat, iat)
         wsw = 1.0_wp/real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))
            call get_cpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, tmp)
            cmat_local(iat, iat) = cmat_local(iat, iat) + tmp*wsw
         end do
      end do
      !$omp end do
      !$omp critical (get_cmat_3d_)
      cmat(:, :) = cmat(:, :) + cmat_local(:, :)
      !$omp end critical (get_cmat_3d_)
      deallocate (cmat_local)
      !$omp end parallel
      !
      cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp

   end subroutine get_cmat_3d

   subroutine get_cpair(kbc, cpair, r1, rvdw, capi, capj)
      real(wp), intent(in) :: r1, capi, capj, rvdw, kbc
      real(wp), intent(out) :: cpair

      real(wp) :: arg

      ! Capacitance of bond between atom i and j
      arg = -kbc*(r1 - rvdw)/rvdw
      cpair = sqrt(capi*capj)*0.5_wp*(1.0_wp + erf(arg))
   end subroutine get_cpair

   subroutine get_cpair_dir(kbc, rij, trans, rvdw, capi, capj, cpair)
      real(wp), intent(in) :: rij(3), capi, capj, rvdw, kbc, trans(:, :)
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
      real(wp), intent(in) :: vec(3), capi, capj, rvdw, kbc
      real(wp), intent(out) :: dgpair(3)
      real(wp), intent(out) :: dspair(3, 3)

      real(wp) :: r1, arg, dtmp

      dgpair(:) = 0.0_wp
      dspair(:, :) = 0.0_wp

      r1 = norm2(vec)
      ! Capacitance of bond between atom i and j
      arg = -(kbc*(r1 - rvdw)/rvdw)**2
      dtmp = sqrt(capi*capj)*kbc*exp(arg)/(sqrtpi*rvdw)
      dgpair = dtmp*vec/r1
      dspair = spread(dgpair, 1, 3)*spread(vec, 2, 3)
   end subroutine get_dcpair

   subroutine get_dcmat_0d(self, mol, dcdr, dcdL)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(out) :: dcdr(:, :, :)
      real(wp), intent(out) :: dcdL(:, :, :)

      integer :: iat, jat, izp, jzp
      real(wp) :: vec(3), r2, rvdw, dtmp, arg, dG(3), dS(3, 3), capi, capj

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dcdr_local(:, :, :), dcdL_local(:, :, :)

      dcdr(:, :, :) = 0.0_wp
      dcdL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(dcdr, dcdL, mol, self) &
      !$omp private(iat, izp, jat, jzp, r2, vec, rvdw) &
      !$omp private(dG, dS, dtmp, arg, capi, capj) &
      !$omp private(dcdr_local, dcdL_local)
      allocate (dcdr_local, source=dcdr)
      allocate (dcdL_local, source=dcdL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(iat, jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)

            call get_dcpair(self%kbc, vec, rvdw, capi, capj, dG, dS)

            ! Negative off-diagonal elements
            dcdr_local(:, iat, jat) = -dG
            dcdr_local(:, jat, iat) = +dG
            ! Positive diagonal elements
            dcdr_local(:, iat, iat) = +dG + dcdr_local(:, iat, iat)
            dcdr_local(:, jat, jat) = -dG + dcdr_local(:, jat, jat)
            dcdL_local(:, :, iat) = -dS + dcdL_local(:, :, iat)
            dcdL_local(:, :, jat) = -dS + dcdL_local(:, :, jat)
         end do
      end do
      !$omp end do
      !$omp critical (get_dcmat_0d_)
      dcdr(:, :, :) = dcdr(:, :, :) + dcdr_local(:, :, :)
      dcdL(:, :, :) = dcdL(:, :, :) + dcdL_local(:, :, :)
      !$omp end critical (get_dcmat_0d_)
      deallocate (dcdL_local, dcdr_local)
      !$omp end parallel

   end subroutine get_dcmat_0d

   subroutine get_dcmat_3d(self, mol, wsc, dcdr, dcdL)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(wignerseitz_cell_type), intent(in) :: wsc
      real(wp), intent(out) :: dcdr(:, :, :)
      real(wp), intent(out) :: dcdL(:, :, :)

      integer :: iat, jat, izp, jzp, img
      real(wp) :: vec(3), r2, rvdw, dtmp, arg, dG(3), dS(3, 3), capi, capj, wsw
      real(wp), allocatable :: dtrans(:, :)

      ! Thread-private arrays for reduction
      real(wp), allocatable :: dcdr_local(:, :, :), dcdL_local(:, :, :)

      call get_dir_trans(mol%lattice, dtrans)

      dcdr(:, :, :) = 0.0_wp
      dcdL(:, :, :) = 0.0_wp

      !$omp parallel default(none) &
      !$omp shared(dcdr, dcdL, mol, self, dtrans, wsc) &
      !$omp private(iat, izp, jat, jzp, r2, vec, rvdw) &
      !$omp private(dG, dS, dtmp, arg, capi, capj, wsw) &
      !$omp private(dcdr_local, dcdL_local)
      allocate (dcdr_local, source=dcdr)
      allocate (dcdL_local, source=dcdL)
      !$omp do schedule(runtime)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            rvdw = self%rvdw(iat, jat)
            wsw = 1/real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, jat) - mol%xyz(:, iat) + wsc%trans(:, wsc%tridx(img, jat, iat))

               call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capj, dG, dS)

               ! Negative off-diagonal elements
               dcdr_local(:, iat, jat) = -dG*wsw + dcdr_local(:, iat, jat)
               dcdr_local(:, jat, iat) = +dG*wsw + dcdr_local(:, jat, iat)
               ! Positive diagonal elements
               dcdr_local(:, iat, iat) = +dG*wsw + dcdr_local(:, iat, iat)
               dcdr_local(:, jat, jat) = -dG*wsw + dcdr_local(:, jat, jat)
               dcdL_local(:, :, jat) = -dS*wsw + dcdL_local(:, :, jat)
               dcdL_local(:, :, iat) = -dS*wsw + dcdL_local(:, :, iat)
            end do
         end do

         rvdw = self%rvdw(iat, iat)
         wsw = 1/real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))

            call get_dcpair_dir(self%kbc, vec, dtrans, rvdw, capi, capi, dG, dS)

            ! Positive diagonal elements
            dcdL_local(:, :, iat) = -dS*wsw + dcdL_local(:, :, iat)
         end do
      end do
      !$omp end do
      !$omp critical (get_dcmat_3d_)
      dcdr(:, :, :) = dcdr(:, :, :) + dcdr_local(:, :, :)
      dcdL(:, :, :) = dcdL(:, :, :) + dcdL_local(:, :, :)
      !$omp end critical (get_dcmat_3d_)
      deallocate (dcdL_local, dcdr_local)
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

   ! NOTE: the following is basically identical to tblite versions of this pattern

   !> Inspect cache and reallocate it in case of type mismatch
   subroutine taint(cache, ptr)
      !> Instance of the cache
      type(cache_container), target, intent(inout) :: cache
      !> Reference to the cache
      type(eeqbc_cache), pointer, intent(out) :: ptr

      if (allocated(cache%raw)) then
         call view(cache, ptr)
         if (associated(ptr)) return
         deallocate (cache%raw)
      end if

      if (.not. allocated(cache%raw)) then
         block
            type(eeqbc_cache), allocatable :: tmp
            allocate (tmp)
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
      type(eeqbc_cache), pointer, intent(out) :: ptr
      nullify (ptr)
      select type (target => cache%raw)
      type is (eeqbc_cache)
         ptr => target
      end select
   end subroutine view

end module multicharge_model_eeqbc
