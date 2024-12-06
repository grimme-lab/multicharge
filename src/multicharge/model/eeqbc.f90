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

   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_convert, only : autoaa
   use mctc_io_math, only : matdet_3x3
   use mctc_ncoord, only : new_ncoord
   use mctc_data, only : get_vdw_rad
   use multicharge_wignerseitz, only : wignerseitz_cell_type
   use multicharge_model_type, only : mchrg_model_type, get_dir_trans, get_rec_trans
   use multicharge_blas, only : gemv, gemm
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
   contains
      !> Calculate right-hand side (electronegativity)
      procedure :: get_vrhs
      !> Calculate Coulomb matrix 
      procedure :: get_amat_0d
      !> Calculate Coulomb matrix periodic
      procedure :: get_amat_3d
      !> Calculate Coulomb matrix derivative
      procedure :: get_damat_0d
      !> Calculate Coulomb matrix derivative periodic
      procedure :: get_damat_3d
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
   & dielectric)
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

   self%chi = chi
   self%rad = rad
   self%eta = eta
   self%kcnchi = kcnchi
   self%kqchi = kqchi
   self%kqeta = kqeta
   self%kcnrad = kcnrad
   self%cap = cap
   self%avg_cn = avg_cn
   
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

subroutine get_vrhs(self, mol, cn, qloc, xvec, dcndr, dcndL, &
   & dqlocdr, dqlocdL, dxdr, dxdL)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: qloc(:)
   real(wp), intent(out) :: xvec(:)
   real(wp), intent(in), optional :: dcndr(:, :, :)
   real(wp), intent(in), optional :: dcndL(:, :, :)
   real(wp), intent(in), optional :: dqlocdr(:, :, :)
   real(wp), intent(in), optional :: dqlocdL(:, :, :)
   real(wp), intent(out), optional :: dxdr(:, :, :)
   real(wp), intent(out), optional :: dxdL(:, :, :)
   
   integer :: iat, izp, jat, jzp
   real(wp) :: tmpdcn, tmpdqloc
   real(wp), allocatable :: tmp(:), dtmpdr(:, :, :), dtmpdL(:, :, :)

   allocate(tmp(mol%nat+1))
   if (present(dxdr) .and. present(dxdL) &
      & .and. present(dcndr) .and. present(dcndL) &
      & .and. present(dqlocdr) .and. present(dqlocdL)) then
      allocate(dtmpdr(3, mol%nat, mol%nat+1), dtmpdL(3, 3, mol%nat+1))

      dxdr(:, :, :) = 0.0_wp
      dxdL(:, :, :) = 0.0_wp
      dtmpdr(:, :, :) = 0.0_wp
      dtmpdL(:, :, :) = 0.0_wp
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(cn, dcndr, dcndL, qloc, dqlocdr, dqlocdL) &
      !$omp shared(self, mol, tmp, dtmpdr, dtmpdL) private(iat, izp) 
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp(iat) = -self%chi(izp) + self%kcnchi(izp)*cn(iat) &
            & + self%kqchi(izp)*qloc(iat)

         ! CN and effective charge derivative
         dtmpdr(:, :, iat) = self%kcnchi(izp) * dcndr(:, :, iat) + dtmpdr(:, :, iat)
         dtmpdL(:, :, iat) = self%kcnchi(izp) * dcndL(:, :, iat) + dtmpdL(:, :, iat)
         dtmpdr(:, :, iat) = self%kqchi(izp) * dqlocdr(:, :, iat) + dtmpdr(:, :, iat)
         dtmpdL(:, :, iat) = self%kqchi(izp) * dqlocdL(:, :, iat) + dtmpdL(:, :, iat)
      end do

      call gemm(dtmpdr(:, :, :mol%nat), cache%cmat(:mol%nat, :mol%nat), dxdr)
      call gemm(dtmpdL(:, :, :mol%nat), cache%cmat(:mol%nat, :mol%nat), dxdL)
      !call gemv(cache%dcdr(:, :, :mol%nat), tmp(:mol%nat), xvec)


      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:dxdr, dxdL) shared(self, mol, cache, tmp) &
      !$omp private(iat, jat)
      do iat = 1, mol%nat
         do jat = 1, mol%nat
            dxdr(:, iat, iat) = tmp(jat) * cache%dcdr(:, iat, jat) + dxdr(:, iat, iat)
            dxdr(:, iat, jat) = (tmp(iat) - tmp(jat)) * cache%dcdr(:, iat, jat) &
               & + dxdr(:, iat, jat)
         end do 
      end do
   else
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cn, qloc, tmp) private(iat, izp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp(iat) = -self%chi(izp) + self%kcnchi(izp)*cn(iat) &
            & + self%kqchi(izp)*qloc(iat)
      end do
   end if
   tmp(mol%nat+1) = mol%charge
   ! xvec = tmp
   call gemv(cache%cmat, tmp, xvec)

end subroutine get_vrhs


subroutine get_amat_0d(self, mol, cn, qloc, amat)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: qloc(:)
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp,  norm_cn, radi, radj

   amat(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:amat) shared(mol, self, cn, qloc, cache) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp, norm_cn, radi, radj)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
      radi = self%rad(izp) * (1.0_wp - self%kcnrad * cn(iat) * norm_cn)
      do jat = 1, iat-1
        jzp = mol%id(jat)
        vec = mol%xyz(:, jat) - mol%xyz(:, iat)
        r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
        ! Effective charge width of j
        norm_cn = cn(jat) / self%avg_cn(jzp)**self%norm_exp
        radj = self%rad(jzp) * (1.0_wp - self%kcnrad * norm_cn)
        ! Coulomb interaction of Gaussian charges
        gam = 1.0_wp / (radi**2 + radj**2)
        tmp = erf(sqrt(r2*gam))/(sqrt(r2)*self%dielectric) * cache%cmat(jat, iat)
        amat(jat, iat) = amat(jat, iat) + tmp
        amat(iat, jat) = amat(iat, jat) + tmp
      end do
      ! Effective hardness
      tmp = self%eta(izp) + self%kqeta(izp) * qloc(iat) + sqrt2pi / radi !
      amat(iat, iat) = amat(iat, iat) + tmp * cache%cmat(iat, iat) + 1.0_wp
   end do

   amat(mol%nat+1, 1:mol%nat+1) = 1.0_wp
   amat(1:mol%nat+1, mol%nat+1) = 1.0_wp
   amat(mol%nat+1, mol%nat+1) = 0.0_wp

end subroutine get_amat_0d

subroutine get_amat_3d(self, mol, wsc, alpha, amat)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), gam, wsw, dtmp, rtmp, vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   amat(:, :) = 0.0_wp

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:amat) shared(mol, self, wsc, dtrans, rtrans, alpha, vol) &
   !$omp private(iat, izp, jat, jzp, gam, wsw, vec, dtmp, rtmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
            call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
            amat(jat, iat) = amat(jat, iat) + (dtmp + rtmp) * wsw
            amat(iat, jat) = amat(iat, jat) + (dtmp + rtmp) * wsw
         end do
      end do

      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
         call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
         amat(iat, iat) = amat(iat, iat) + (dtmp + rtmp) * wsw
      end do

      dtmp = self%eta(izp) + sqrt2pi / self%rad(izp) - 2 * alpha / sqrtpi
      amat(iat, iat) = amat(iat, iat) + dtmp
   end do

   amat(mol%nat+1, 1:mol%nat+1) = 1.0_wp
   amat(1:mol%nat+1, mol%nat+1) = 1.0_wp
   amat(mol%nat+1, mol%nat+1) = 0.0_wp

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
      tmp = cos(dot_product(rij, vec)) * fac * exp(-0.25_wp*g2/(alp*alp))/g2
      amat = amat + tmp
   end do

end subroutine get_amat_rec_3d

subroutine get_damat_0d(self, mol, cn, qloc, qvec, dcndr, dcndL, &
      & dqlocdr, dqlocdL, dadr, dadL, atrace)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: cn(:)
   real(wp), intent(in) :: qloc(:)
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(in) :: dcndr(:, :, :)
   real(wp), intent(in) :: dcndL(:, :, :)
   real(wp), intent(in) :: dqlocdr(:, :, :)
   real(wp), intent(in) :: dqlocdL(:, :, :)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, arg, dtmp, norm_cn
   real(wp) :: radi, radj, dradi, dradj, dG(3), dS(3, 3), dgamdL(3, 3)
   real(wp), allocatable :: dgamdr(:, :)

   allocate(dgamdr(3, mol%nat))
   
   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace, dadr, dadL) shared(self, mol, cn, qloc, qvec) & 
   !$omp shared (cache, dcndr, dcndL, dqlocdr, dqlocdL) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, dtmp, norm_cn, arg) &
   !$omp private(radi, radj, dradi, dradj, dgamdr, dgamdL, dG, dS)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      ! Effective charge width of i
      norm_cn = 1.0_wp / self%avg_cn(izp)**self%norm_exp
      radi = self%rad(izp) * (1.0_wp - self%kcnrad * cn(iat) * norm_cn)
      dradi = -self%rad(izp) * self%kcnrad * norm_cn
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! Effective charge width of j
         norm_cn = 1.0_wp / self%avg_cn(jzp)**self%norm_exp
         radj = self%rad(jzp) * (1.0_wp - self%kcnrad * cn(jat) * norm_cn)
         dradj = -self%rad(jzp) * self%kcnrad * norm_cn
         
         ! Coulomb interaction of Gaussian charges
         gam = 1.0_wp / sqrt(radi**2 + radj**2)
         dgamdr(:, :) = - (radi*dradi*dcndr(:, :, iat) + radj*dradj*dcndr(:, :, jat)) &
            & * gam**3.0_wp
         dgamdL(:, :) = - (radi*dradi*dcndL(:, :, iat) + radj*dradj*dcndL(:, :, jat)) &
            & * gam**3.0_wp

         ! Explicit derivative
         arg = gam*gam*r2
         dtmp = 2.0_wp*gam*exp(-arg)/(sqrtpi*r2*self%dielectric) &
            & - erf(sqrt(arg))/(r2*sqrt(r2)*self%dielectric)
         dG(:) = -dtmp * vec ! questionable sign
         dS(:, :) = spread(dG, 1, 3) * spread(vec, 2, 3)
         atrace(:, iat)    = +dG*qvec(jat)*cache%cmat(jat, iat) + atrace(:, iat)
         atrace(:, jat)    = -dG*qvec(iat)*cache%cmat(iat, jat) + atrace(:, jat)
         dadr(:, iat, jat) = +dG*qvec(iat)*cache%cmat(iat, jat) + dadr(:, iat, jat)
         dadr(:, jat, iat) = -dG*qvec(jat)*cache%cmat(jat, iat) + dadr(:, jat, iat)
         dadL(:, :, jat)   = +dS*qvec(iat)*cache%cmat(iat, jat) + dadL(:, :, jat)
         dadL(:, :, iat)   = +dS*qvec(jat)*cache%cmat(jat, iat) + dadL(:, :, iat)

         ! Effective charge width derivative 
         dtmp = 2.0_wp*exp(-arg)/(sqrtpi*self%dielectric)
         atrace(:, iat)    = +dtmp*qvec(jat)*dgamdr(:, jat)*cache%cmat(jat, iat) + atrace(:, iat)
         atrace(:, jat)    = +dtmp*qvec(iat)*dgamdr(:, iat)*cache%cmat(iat, jat) + atrace(:, jat)
         dadr(:, iat, jat) = +dtmp*qvec(iat)*dgamdr(:, iat)*cache%cmat(iat, jat) + dadr(:, iat, jat)
         dadr(:, jat, iat) = +dtmp*qvec(jat)*dgamdr(:, jat)*cache%cmat(jat, iat) + dadr(:, jat, iat)
         dadL(:, :, jat)   = +dtmp*qvec(iat)*dgamdL(:, :)*cache%cmat(iat, jat)   + dadL(:, :, jat)
         dadL(:, :, iat)   = +dtmp*qvec(jat)*dgamdL(:, :)*cache%cmat(jat, iat)   + dadL(:, :, iat)

         ! Capacitance derivative 
         dtmp = erf(sqrt(r2)*gam)/(sqrt(r2)*self%dielectric)
         ! potentially switch indices for dcdr
         atrace(:, iat)    = +dtmp*qvec(jat)*cache%dcdr(:, jat, iat) + atrace(:, iat)
         atrace(:, jat)    = +dtmp*qvec(iat)*cache%dcdr(:, iat, jat) + atrace(:, jat)
         dadr(:, iat, jat) = +dtmp*qvec(iat)*cache%dcdr(:, iat, jat) + dadr(:, iat, jat)
         dadr(:, jat, iat) = +dtmp*qvec(jat)*cache%dcdr(:, jat, iat) + dadr(:, jat, iat)
         dadL(:, :, jat)   = +dtmp*qvec(iat)*cache%dcdL(:, :, iat)   + dadL(:, :, jat)
         dadL(:, :, iat)   = +dtmp*qvec(jat)*cache%dcdL(:, :, jat)   + dadL(:, :, iat)
      end do

      ! Hardness derivative
      dtmp = self%kqeta(izp) * qvec(iat) * cache%cmat(iat, iat)
      !atrace(:, iat)    = -dtmp*dqlocdr(:, iat, iat) + atrace(:, iat)
      dadr(:, iat, iat) = +dtmp*dqlocdr(:, iat, iat) + dadr(:, iat, iat) 
      dadL(:, :, iat)   = +dtmp*dqlocdL(:, :, iat)   + dadL(:, :, iat)

      ! Effective charge width derivative 
      dtmp = -sqrt2pi*dradi/(radi**2) * qvec(iat) * cache%cmat(iat, iat)
      !atrace(:, iat)    = -dtmp*dcndr(:, iat, iat) + atrace(:, iat)
      dadr(:, iat, iat) = +dtmp*dcndr(:, iat, iat) + dadr(:, iat, iat)
      dadL(:, :, iat)   = +dtmp*dcndL(:, :, iat)   + dadL(:, :, iat)

      ! Capacitance derivative 
      dtmp = (self%eta(izp) + self%kqeta(izp)*qloc(iat) + sqrt2pi/radi) * qvec(iat)
      !atrace(:, iat)    = -dtmp*cache%dcdr(:, iat, iat) + atrace(:, iat)
      dadr(:, iat, iat) = +dtmp*cache%dcdr(:, iat, iat) + dadr(:, iat, iat)
      dadL(:, :, iat)   = +dtmp*cache%dcdL(:, :, iat)   + dadL(:, :, iat)

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
      do jat = 1, iat-1
         jzp = mol%id(jat)
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
            call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         atrace(:, iat) = +dG*qvec(jat) + atrace(:, iat)
         atrace(:, jat) = -dG*qvec(iat) + atrace(:, jat)
         dadr(:, iat, jat) = +dG*qvec(iat) + dadr(:, iat, jat)
         dadr(:, jat, iat) = -dG*qvec(jat) + dadr(:, jat, iat)
         dadL(:, :, jat) = +dS*qvec(iat) + dadL(:, :, jat)
         dadL(:, :, iat) = +dS*qvec(jat) + dadL(:, :, iat)
      end do

      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(2.0_wp * self%rad(izp)**2)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_damat_dir_3d(vec, gam, alpha, dtrans, dGd, dSd)
         call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
         dS = dS + (dSd + dSr) * wsw
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
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * spread(vec, 1, 3) * spread(vec, 2, 3)
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
      etmp = fac * exp(-0.25_wp*g2/alp2)/g2
      dtmp = -sin(gv) * etmp
      dg(:) = dg + dtmp * vec
      ds(:, :) = ds + etmp * cos(gv) &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_rec_3d


subroutine get_cmat(self, mol, cmat)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: cmat(:, :)

   integer :: iat, jat, izp, jzp, isp, jsp
   real(wp) :: vec(3), r2, rvdw, tmp, arg

   cmat(:, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:cmat) shared(mol, self) &
   !$omp private(iat, izp, isp, jat, jzp, jsp, r2) & 
   !$omp private(vec, rvdw, tmp, arg)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      isp = mol%num(izp)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jsp = mol%num(jzp)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! vdw distance in Angstrom (approximate factor 2)
         rvdw = get_vdw_rad(isp, jsp) * autoaa
         ! Capacitance of bond between atom i and j
         arg = -self%kbc * (sqrt(r2) - rvdw) / rvdw
         tmp = sqrt(self%cap(izp) * self%cap(jzp)) * &
            & 0.5_wp * (1.0_wp + erf(arg))
         ! Off-diagonal elements
         cmat(jat, iat) = - tmp
         cmat(iat, jat) = - tmp
         ! Diagonal elements
         cmat(iat, iat) = cmat(iat, iat) + tmp
         cmat(jat, jat) = cmat(jat, jat) + tmp
      end do
   end do
   cmat(mol%nat+1, mol%nat+1) = 1.0_wp

end subroutine get_cmat

subroutine get_dcmat_0d(self, mol, dcdr, dcdL)
   class(eeqbc_model), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out), optional :: dcdr(:, :, :)
   real(wp), intent(out), optional :: dcdL(:, :, :)

   integer :: iat, jat, izp, jzp, isp, jsp
   real(wp) :: vec(3), r2, rvdw, dtmp, arg, dG(3), dS(3, 3)

   dcdr(:, :, :) = 0.0_wp
   dcdL(:, :, :) = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:dcdr, dcdL) shared(mol, self) &
   !$omp private(iat, izp, isp, jat, jzp, jsp, r2) & 
   !$omp private(vec, rvdw, dG, dS, dtmp, arg)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      isp = mol%num(izp)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         jsp = mol%num(jzp)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         ! vdw distance in Angstrom (approximate factor 2)
         rvdw = get_vdw_rad(isp, jsp) * autoaa

         ! Capacitance of bond between atom i and j
         arg = -(self%kbc * (sqrt(r2) - rvdw) / rvdw)**2
         dtmp = sqrt(self%cap(izp) * self%cap(jzp)) * &
            & self%kbc * exp(arg) / (sqrtpi * rvdw)
         dG = dtmp*vec/sqrt(r2)
         dS = spread(dG, 1, 3) * spread(vec, 2, 3)

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

end module multicharge_model_eeqbc
