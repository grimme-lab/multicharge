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

   use iso_fortran_env, only: output_unit

   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_convert, only: autoaa
   use mctc_io_math, only: matdet_3x3
   use mctc_ncoord, only: new_ncoord
   use multicharge_wignerseitz, only: wignerseitz_cell_type
   use multicharge_model_type, only: mchrg_model_type, get_dir_trans, get_rec_trans
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
      !> Full constraint matrix for 0d case
      real(wp), allocatable :: cmat(:, :)
      !> Contributions for every WSC image for diagonal elements of constraint matrix
      real(wp), allocatable :: cmat_diag(:, :)
      !> Derivative of constraint matrix w.r.t positions
      real(wp), allocatable :: dcdr(:, :, :)
      !> Derivative of constraint matrix w.r.t lattice vectors
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
      !> Calculate Coulomb matrix derivative periodic
      procedure :: get_damat_3d
      !> Calculate constraint matrix (molecular case)
      procedure :: get_cmat_0d
      !> Calculate diagonal contributions (periodic case)
      procedure :: get_cmat_diag_3d
      !> Calculate constraint matrix derivatives (molecular)
      procedure :: get_dcmat_0d
      ! procedure :: get_dcmat_3d
   end type eeqbc_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

   !> Default exponent of distance/CN normalization
   real(wp), parameter :: default_norm_exp = 1.0_wp

   !> Default exponent of error function in bond capacitance
   real(wp), parameter :: default_kbc = 0.65_wp
contains

   subroutine new_eeqbc_model(self, mol, chi, rad, eta, kcnchi, kqchi, kqeta, &
      & kcnrad, cap, avg_cn, kbc, cutoff, cn_exp, rcov, en, cn_max, norm_exp, &
      & dielectric, rvdw)
      !> Bond capacitor electronegativity equilibration model
      type(eeqbc_model), intent(out) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
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
      !> Dielectric constant of the surrounding medium
      real(wp), intent(in), optional :: dielectric
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

      if (present(dielectric)) then
         self%dielectric = dielectric
      else
         self%dielectric = 1.0_wp
      end if

      ! Coordination number
      call new_ncoord(self%ncoord, mol, "erf", cutoff=cutoff, kcn=cn_exp, &
         & rcov=rcov, cut=cn_max, norm_exp=self%norm_exp)
      ! Electronegativity weighted coordination number for local charge
      call new_ncoord(self%ncoord_en, mol, "erf_en", cutoff=cutoff, kcn=cn_exp, &
         & rcov=rcov, en=en, cut=cn_max, norm_exp=self%norm_exp)

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
      call ptr%update(mol)

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

      if (any(mol%periodic)) then
         ! Allocate cmat diagonal WSC image contributions
         if (.not. allocated(ptr%cmat_diag)) then
            allocate (ptr%cmat_diag(mol%nat, ptr%wsc%nimg_max))
         end if
         ! Get cmat diagonal contributions for all WSC images
         call self%get_cmat_diag_3d(mol, ptr%wsc, ptr%cmat_diag)
         ! if (grad) then
         ! call self%get_dcmat_3d()
         ! end if
      else
         ! Allocate cmat
         if (.not. allocated(ptr%cmat)) then
            allocate (ptr%cmat(mol%nat + 1, mol%nat + 1))
         end if
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

      integer :: iat, izp

      call view(cache, ptr)

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, ptr) private(iat, izp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ptr%xtmp(iat) = -self%chi(izp) + self%kcnchi(izp)*ptr%cn(iat) &
            & + self%kqchi(izp)*ptr%qloc(iat)
      end do
      ptr%xtmp(mol%nat + 1) = mol%charge
      call gemv(ptr%cmat, ptr%xtmp, xvec)

   end subroutine get_xvec

   subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out) :: dxdr(:, :, :)
      real(wp), intent(out) :: dxdL(:, :, :)

      type(eeqbc_cache), pointer :: ptr

      integer :: iat, izp, jat
      real(wp) :: tmpdcn, tmpdqloc
      real(wp), allocatable :: dtmpdr(:, :, :), dtmpdL(:, :, :)

      call view(cache, ptr)
      allocate (dtmpdr(3, mol%nat, mol%nat + 1), dtmpdL(3, 3, mol%nat + 1))

      dxdr(:, :, :) = 0.0_wp
      dxdL(:, :, :) = 0.0_wp
      dtmpdr(:, :, :) = 0.0_wp
      dtmpdL(:, :, :) = 0.0_wp
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(ptr, self, mol, dtmpdr, dtmpdL) &
      !$omp private(iat, izp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         ! CN and effective charge derivative
         dtmpdr(:, :, iat) = self%kcnchi(izp)*ptr%dcndr(:, :, iat) + dtmpdr(:, :, iat)
         dtmpdL(:, :, iat) = self%kcnchi(izp)*ptr%dcndL(:, :, iat) + dtmpdL(:, :, iat)
         dtmpdr(:, :, iat) = self%kqchi(izp)*ptr%dqlocdr(:, :, iat) + dtmpdr(:, :, iat)
         dtmpdL(:, :, iat) = self%kqchi(izp)*ptr%dqlocdL(:, :, iat) + dtmpdL(:, :, iat)
      end do

      call gemm(dtmpdr(:, :, :mol%nat), ptr%cmat(:mol%nat, :mol%nat), dxdr)
      call gemm(dtmpdL(:, :, :mol%nat), ptr%cmat(:mol%nat, :mol%nat), dxdL)
      !call gemv(cache%dcdr(:, :, :mol%nat), tmp(:mol%nat), xvec)
      print'(a)', 'xtmp:'
      print'(3es21.14)', ptr%xtmp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:dxdr, dxdL) shared(self, mol, ptr) &
      !$omp private(iat, jat)
      do iat = 1, mol%nat
         do jat = 1, mol%nat
            dxdr(:, iat, iat) = ptr%xtmp(jat)*ptr%dcdr(:, iat, jat) + dxdr(:, iat, iat)
            dxdr(:, iat, jat) = (ptr%xtmp(iat) - ptr%xtmp(jat))*ptr%dcdr(:, iat, jat) &
               & + dxdr(:, iat, jat)
         end do
      end do
   end subroutine get_xvec_derivs

   subroutine get_coulomb_matrix(self, mol, cache, amat)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out) :: amat(:, :)

      type(eeqbc_cache), pointer :: ptr
      call view(cache, ptr)

      if (any(mol%periodic)) then
         call self%get_amat_3d(mol, ptr%wsc, ptr%alpha, ptr%cn, &
         & ptr%qloc, amat, ptr%cmat_diag)
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

      amat(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:amat) shared(mol, self, cn, qloc, cmat) &
      !$omp private(iat, izp, jat, jzp, gam2, vec, r2, tmp, norm_cn, radi, radj)
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
            tmp = erf(sqrt(r2*gam2))/(sqrt(r2)*self%dielectric)*cmat(jat, iat)
            amat(jat, iat) = amat(jat, iat) + tmp
            amat(iat, jat) = amat(iat, jat) + tmp
         end do
         ! Effective hardness
         tmp = self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi !
         amat(iat, iat) = amat(iat, iat) + tmp*cmat(iat, iat) + 1.0_wp
      end do

      amat(mol%nat + 1, 1:mol%nat + 1) = 1.0_wp
      amat(1:mol%nat + 1, mol%nat + 1) = 1.0_wp
      amat(mol%nat + 1, mol%nat + 1) = 0.0_wp

   end subroutine get_amat_0d

   subroutine get_amat_3d(self, mol, wsc, alpha, cn, qloc, amat, cmat_diag)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(wignerseitz_cell_type), intent(in) :: wsc
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: cn(:), qloc(:)
      real(wp), intent(out) :: amat(:, :)
      real(wp), intent(out) :: cmat_diag(:, :)

      integer :: iat, jat, isp, jsp, izp, jzp, img
      real(wp) :: vec(3), gam, wsw, dtmp, rtmp, vol, ctmp, capi, capj, radi, radj, norm_cn, rvdw
      real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

      amat(:, :) = 0.0_wp

      vol = abs(matdet_3x3(mol%lattice))
      call get_dir_trans(mol%lattice, dtrans)
      call get_rec_trans(mol%lattice, rtrans)

      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:amat) shared(mol, self, wsc, dtrans, rtrans, alpha, vol, cmat_diag) &
      !$omp shared(cn, qloc) &
      !$omp private(iat, izp, jat, jzp, gam, wsw, vec, dtmp, rtmp, ctmp, norm_cn) &
      !$omp private(isp, jsp, radi, radj, capi, capj, rvdw)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         isp = mol%num(izp)
         ! Effective charge width of i
         norm_cn = 1.0_wp/self%avg_cn(izp)**self%norm_exp
         radi = self%rad(izp)*(1.0_wp - self%kcnrad*cn(iat)*norm_cn)
         capi = self%cap(isp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            jsp = mol%num(jzp)
            ! vdw distance in Angstrom (approximate factor 2)
            rvdw = self%rvdw(iat, jat)
            ! Effective charge width of j
            norm_cn = cn(jat)/self%avg_cn(jzp)**self%norm_exp
            radj = self%rad(jzp)*(1.0_wp - self%kcnrad*norm_cn)
            capj = self%cap(jsp)
            ! Coulomb interaction of Gaussian charges
            gam = 1.0_wp/sqrt(radi**2 + radj**2)
            wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
               call get_cmat_pair(mol, self%kbc, ctmp, vec, rvdw, capi, capj)
               call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
               call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
               amat(jat, iat) = amat(jat, iat) + ctmp*(dtmp + rtmp)*wsw
               amat(iat, jat) = amat(iat, jat) + ctmp*(dtmp + rtmp)*wsw
            end do
         end do

         ! WSC image contributions
         gam = 1.0_wp/sqrt(2.0_wp*self%rad(izp)**2)
         wsw = 1.0_wp/real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))
            ctmp = cmat_diag(iat, img)
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
            call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
            amat(iat, iat) = amat(iat, iat) + ctmp*(dtmp + rtmp)*wsw
         end do
         ! Effective hardness
         dtmp = self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi
         amat(iat, iat) = amat(iat, iat) + cmat_diag(iat, 1)*dtmp + 1.0_wp
      end do

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
         tmp = erf(gam*r1)/r1 - erf(alp*r1)/r1
         amat = amat + tmp
      end do

   end subroutine get_amat_dir_3d

   subroutine get_amat_rec_3d(rij, vol, alp, trans, amat)
      real(wp), intent(in) :: rij(3)
      real(wp), intent(in) :: vol
      real(wp), intent(in) :: alp
      real(wp), intent(in) :: trans(:, :)
      real(wp), intent(out) :: amat

      integer :: itr
      real(wp) :: fac, vec(3), g2, tmp

      amat = 0.0_wp
      fac = 4*pi/vol

      do itr = 1, size(trans, 2)
         vec(:) = trans(:, itr)
         g2 = dot_product(vec, vec)
         if (g2 < eps) cycle
         tmp = cos(dot_product(rij, vec))*fac*exp(-0.25_wp*g2/(alp*alp))/g2
         amat = amat + tmp
      end do

   end subroutine get_amat_rec_3d

   subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: qvec(:)
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

      type(eeqbc_cache), pointer :: ptr
      call view(cache, ptr)

      if (any(mol%periodic)) then
         call self%get_damat_3d(mol, ptr%wsc, ptr%alpha, qvec, dadr, dadL, atrace)
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

      allocate (dgamdr(3, mol%nat))

      atrace(:, :) = 0.0_wp
      dadr(:, :, :) = 0.0_wp
      dadL(:, :, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:atrace, dadr, dadL) shared(self, mol, cn, qloc, qvec) &
      !$omp shared (cmat, dcdr, dcdL, dcndr, dcndL, dqlocdr, dqlocdL) &
      !$omp private(iat, izp, jat, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
      !$omp private(radi, radj, dradi, dradj, dgamdr, dgamdL, dG, dS)
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
            dtmp = 2.0_wp*gam*exp(-arg)/(sqrtpi*r2*self%dielectric) &
               & - erf(sqrt(arg))/(r2*sqrt(r2)*self%dielectric)
            dG(:) = -dtmp*vec ! questionable sign
            dS(:, :) = spread(dG, 1, 3)*spread(vec, 2, 3)
            atrace(:, iat) = +dG*qvec(jat)*cmat(jat, iat) + atrace(:, iat)
            atrace(:, jat) = -dG*qvec(iat)*cmat(iat, jat) + atrace(:, jat)
            dadr(:, iat, jat) = +dG*qvec(iat)*cmat(iat, jat) + dadr(:, iat, jat)
            dadr(:, jat, iat) = -dG*qvec(jat)*cmat(jat, iat) + dadr(:, jat, iat)
            dadL(:, :, jat) = +dS*qvec(iat)*cmat(iat, jat) + dadL(:, :, jat)
            dadL(:, :, iat) = +dS*qvec(jat)*cmat(jat, iat) + dadL(:, :, iat)

            ! Effective charge width derivative
            dtmp = 2.0_wp*exp(-arg)/(sqrtpi*self%dielectric)
            atrace(:, iat) = +dtmp*qvec(jat)*dgamdr(:, jat)*cmat(jat, iat) + atrace(:, iat)
            atrace(:, jat) = +dtmp*qvec(iat)*dgamdr(:, iat)*cmat(iat, jat) + atrace(:, jat)
            dadr(:, iat, jat) = +dtmp*qvec(iat)*dgamdr(:, iat)*cmat(iat, jat) + dadr(:, iat, jat)
            dadr(:, jat, iat) = +dtmp*qvec(jat)*dgamdr(:, jat)*cmat(jat, iat) + dadr(:, jat, iat)
            dadL(:, :, jat) = +dtmp*qvec(iat)*dgamdL(:, :)*cmat(iat, jat) + dadL(:, :, jat)
            dadL(:, :, iat) = +dtmp*qvec(jat)*dgamdL(:, :)*cmat(jat, iat) + dadL(:, :, iat)

            ! Capacitance derivative
            dtmp = erf(sqrt(r2)*gam)/(sqrt(r2)*self%dielectric)
            ! potentially switch indices for dcdr
            atrace(:, iat) = +dtmp*qvec(jat)*dcdr(:, jat, iat) + atrace(:, iat)
            atrace(:, jat) = +dtmp*qvec(iat)*dcdr(:, iat, jat) + atrace(:, jat)
            dadr(:, iat, jat) = +dtmp*qvec(iat)*dcdr(:, iat, jat) + dadr(:, iat, jat)
            dadr(:, jat, iat) = +dtmp*qvec(jat)*dcdr(:, jat, iat) + dadr(:, jat, iat)
            dadL(:, :, jat) = +dtmp*qvec(iat)*dcdL(:, :, iat) + dadL(:, :, jat)
            dadL(:, :, iat) = +dtmp*qvec(jat)*dcdL(:, :, jat) + dadL(:, :, iat)
         end do

         ! Hardness derivative
         dtmp = self%kqeta(izp)*qvec(iat)*cmat(iat, iat)
         !atrace(:, iat)    = -dtmp*dqlocdr(:, iat, iat) + atrace(:, iat)
         dadr(:, iat, iat) = +dtmp*dqlocdr(:, iat, iat) + dadr(:, iat, iat)
         dadL(:, :, iat) = +dtmp*dqlocdL(:, :, iat) + dadL(:, :, iat)

         ! Effective charge width derivative
         dtmp = -sqrt2pi*dradi/(radi**2)*qvec(iat)*cmat(iat, iat)
         !atrace(:, iat)    = -dtmp*dcndr(:, iat, iat) + atrace(:, iat)
         dadr(:, iat, iat) = +dtmp*dcndr(:, iat, iat) + dadr(:, iat, iat)
         dadL(:, :, iat) = +dtmp*dcndL(:, :, iat) + dadL(:, :, iat)

         ! Capacitance derivative
         dtmp = (self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi)*qvec(iat)
         !atrace(:, iat)    = -dtmp*dcdr(:, iat, iat) + atrace(:, iat)
         dadr(:, iat, iat) = +dtmp*dcdr(:, iat, iat) + dadr(:, iat, iat)
         dadL(:, :, iat) = +dtmp*dcdL(:, :, iat) + dadL(:, :, iat)

      end do

   end subroutine get_damat_0d

   subroutine get_damat_3d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(wignerseitz_cell_type), intent(in) :: wsc
      real(wp), intent(in) :: alpha
      real(wp), intent(in) :: qvec(:)
      real(wp), intent(out) :: dadr(:, :, :)
      real(wp), intent(out) :: dadL(:, :, :)
      real(wp), intent(out) :: atrace(:, :)

      integer :: iat, jat, izp, jzp, img
      real(wp) :: vol, gam, wsw, vec(3), dG(3), dS(3, 3)
      real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
      real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

      atrace(:, :) = 0.0_wp
      dadr(:, :, :) = 0.0_wp
      dadL(:, :, :) = 0.0_wp

      vol = abs(matdet_3x3(mol%lattice))
      call get_dir_trans(mol%lattice, dtrans)
      call get_rec_trans(mol%lattice, rtrans)

      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:atrace, dadr, dadL) &
      !$omp shared(mol, self, wsc, alpha, vol, dtrans, rtrans, qvec) &
      !$omp private(iat, izp, jat, jzp, img, gam, wsw, vec, dG, dS, &
      !$omp& dGr, dSr, dGd, dSd)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            dG(:) = 0.0_wp
            dS(:, :) = 0.0_wp
            gam = 1.0_wp/sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
            wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
               call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
               call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
               dG = dG + (dGd + dGr)*wsw
               dS = dS + (dSd + dSr)*wsw
            end do
            atrace(:, iat) = +dG*qvec(jat) + atrace(:, iat)
            atrace(:, jat) = -dG*qvec(iat) + atrace(:, jat)
            dadr(:, iat, jat) = +dG*qvec(iat) + dadr(:, iat, jat)
            dadr(:, jat, iat) = -dG*qvec(jat) + dadr(:, jat, iat)
            dadL(:, :, jat) = +dS*qvec(iat) + dadL(:, :, jat)
            dadL(:, :, iat) = +dS*qvec(jat) + dadL(:, :, iat)
         end do

         dS(:, :) = 0.0_wp
         gam = 1.0_wp/sqrt(2.0_wp*self%rad(izp)**2)
         wsw = 1.0_wp/real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))
            call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
            call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
            dS = dS + (dSd + dSr)*wsw
         end do
         dadL(:, :, iat) = +dS*qvec(iat) + dadL(:, :, iat)
      end do

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

      gam2 = gam*gam
      alp2 = alp*alp

      do itr = 1, size(trans, 2)
         vec(:) = rij + trans(:, itr)
         r1 = norm2(vec)
         if (r1 < eps) cycle
         r2 = r1*r1
         gtmp = +2*gam*exp(-r2*gam2)/(sqrtpi*r2) - erf(r1*gam)/(r2*r1)
         atmp = -2*alp*exp(-r2*alp2)/(sqrtpi*r2) + erf(r1*alp)/(r2*r1)
         dg(:) = dg + (gtmp + atmp)*vec
         ds(:, :) = ds + (gtmp + atmp)*spread(vec, 1, 3)*spread(vec, 2, 3)
      end do

   end subroutine get_damat_dir_3d

   subroutine get_damat_rec_3d(rij, vol, alp, trans, dg, ds)
      real(wp), intent(in) :: rij(3)
      real(wp), intent(in) :: vol
      real(wp), intent(in) :: alp
      real(wp), intent(in) :: trans(:, :)
      real(wp), intent(out) :: dg(3)
      real(wp), intent(out) :: ds(3, 3)

      integer :: itr
      real(wp) :: fac, vec(3), g2, gv, etmp, dtmp, alp2
      real(wp), parameter :: unity(3, 3) = reshape(&
         & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

      dg(:) = 0.0_wp
      ds(:, :) = 0.0_wp
      fac = 4*pi/vol
      alp2 = alp*alp

      do itr = 1, size(trans, 2)
         vec(:) = trans(:, itr)
         g2 = dot_product(vec, vec)
         if (g2 < eps) cycle
         gv = dot_product(rij, vec)
         etmp = fac*exp(-0.25_wp*g2/alp2)/g2
         dtmp = -sin(gv)*etmp
         dg(:) = dg + dtmp*vec
         ds(:, :) = ds + etmp*cos(gv) &
                   & *((2.0_wp/g2 + 0.5_wp/alp2)*spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
      end do

   end subroutine get_damat_rec_3d

   subroutine get_cmat_0d(self, mol, cmat)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(out) :: cmat(:, :)

      integer :: iat, jat, izp, jzp, isp, jsp
      real(wp) :: vec(3), rvdw, tmp, capi, capj

      cmat(:, :) = 0.0_wp
      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:cmat) shared(mol, self) &
      !$omp private(iat, izp, isp, jat, jzp, jsp) &
      !$omp private(vec, rvdw, tmp, capi, capj)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         isp = mol%num(izp)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            jsp = mol%num(jzp)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            rvdw = self%rvdw(iat, jat)
            capj = self%cap(jzp)
            call get_cmat_pair(mol, self%kbc, tmp, vec, rvdw, capi, capj)
            ! Off-diagonal elements
            cmat(jat, iat) = -tmp
            cmat(iat, jat) = -tmp
            ! Diagonal elements
            cmat(iat, iat) = cmat(iat, iat) + tmp
            cmat(jat, jat) = cmat(jat, jat) + tmp
         end do
      end do
      cmat(mol%nat + 1, mol%nat + 1) = 1.0_wp

   end subroutine get_cmat_0d

   subroutine get_cmat_pair(mol, kbc, cmat, vec, rvdw, capi, capj)
      type(structure_type), intent(in) :: mol
      real(wp), intent(in) :: vec(3), capi, capj, rvdw, kbc
      real(wp), intent(out) :: cmat

      real(wp) :: r2, arg

      r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
      ! Capacitance of bond between atom i and j
      arg = -kbc*(sqrt(r2) - rvdw)/rvdw
      cmat = sqrt(capi*capj)*0.5_wp*(1.0_wp + erf(arg))
   end subroutine get_cmat_pair

   subroutine get_cmat_diag_3d(self, mol, wsc, cmat)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(wignerseitz_cell_type), intent(in) :: wsc
      real(wp), intent(out) :: cmat(:, :)

      integer :: iat, jat, izp, jzp, isp, jsp, img
      real(wp) :: vec(3), rvdw, capi, capj, tmp

      cmat(:, :) = 0.0_wp
      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:cmat) shared(mol, self, wsc) &
      !$omp private(iat, izp, isp, jat, jzp, jsp, img) &
      !$omp private(vec, rvdw, tmp, capi, capj)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         isp = mol%num(izp)
         capi = self%cap(isp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            jsp = mol%num(jzp)
            rvdw = self%rvdw(iat, jat)
            capj = self%cap(jsp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
               call get_cmat_pair(mol, self%kbc, tmp, vec, rvdw, capi, capj)
               cmat(iat, img) = cmat(iat, img) + tmp
               cmat(jat, img) = cmat(jat, img) + tmp
            end do
         end do
      end do
   end subroutine get_cmat_diag_3d

   subroutine get_dcmat_0d(self, mol, dcdr, dcdL)
      class(eeqbc_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(out) :: dcdr(:, :, :)
      real(wp), intent(out) :: dcdL(:, :, :)

      integer :: iat, jat, izp, jzp
      real(wp) :: vec(3), r2, rvdw, dtmp, arg, dG(3), dS(3, 3), capi, capj

      dcdr(:, :, :) = 0.0_wp
      dcdL(:, :, :) = 0.0_wp
      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:dcdr, dcdL) shared(mol, self) &
      !$omp private(iat, izp, jat, jzp, r2) &
      !$omp private(vec, rvdw, dG, dS, dtmp, arg, capi, capj)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         capi = self%cap(izp)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            capj = self%cap(jzp)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            rvdw = self%rvdw(iat, jat)

            ! Capacitance of bond between atom i and j
            arg = -(self%kbc*(sqrt(r2) - rvdw)/rvdw)**2
            dtmp = sqrt(capi*capj)* &
               & self%kbc*exp(arg)/(sqrtpi*rvdw)
            dG = dtmp*vec/sqrt(r2)
            dS = spread(dG, 1, 3)*spread(vec, 2, 3)

            ! Negative off-diagonal elements
            dcdr(:, iat, jat) = -dG
            dcdr(:, jat, iat) = +dG
            ! Positive diagonal elements
            dcdr(:, iat, iat) = +dG + dcdr(:, iat, iat)
            dcdr(:, jat, jat) = -dG + dcdr(:, jat, jat)
            dcdL(:, :, jat) = +dS + dcdL(:, :, jat)
            dcdL(:, :, iat) = +dS + dcdL(:, :, iat)
         end do
      end do

   end subroutine get_dcmat_0d

   subroutine write_2d_matrix(matrix, name, unit, step)
      implicit none
      real(wp), intent(in) :: matrix(:, :)
      character(len=*), intent(in), optional :: name
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: step
      integer :: d1, d2
      integer :: i, j, k, l, istep, iunit

      d1 = size(matrix, dim=1)
      d2 = size(matrix, dim=2)

      if (present(unit)) then
         iunit = unit
      else
         iunit = output_unit
      end if

      if (present(step)) then
         istep = step
      else
         istep = 6
      end if

      if (present(name)) write (iunit, '(/,"matrix printed:",1x,a)') name

      do i = 1, d2, istep
         l = min(i + istep - 1, d2)
         write (iunit, '(/,6x)', advance='no')
         do k = i, l
            write (iunit, '(6x,i7,3x)', advance='no') k
         end do
         write (iunit, '(a)')
         do j = 1, d1
            write (iunit, '(i6)', advance='no') j
            do k = i, l
               write (iunit, '(1x,f15.8)', advance='no') matrix(j, k)
            end do
            write (iunit, '(a)')
         end do
      end do

   end subroutine write_2d_matrix

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
