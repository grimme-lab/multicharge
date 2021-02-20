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

module multicharge_model
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use multicharge_blas, only : gemv, symv, gemm
   use multicharge_lapack, only : sytrf, sytrs, sytri
   implicit none
   private

   public :: mchrg_model_type, new_mchrg_model


   type :: mchrg_model_type
      real(wp), allocatable :: rad(:)
      real(wp), allocatable :: chi(:)
      real(wp), allocatable :: eta(:)
      real(wp), allocatable :: kcn(:)
   contains
      procedure :: solve
   end type mchrg_model_type


   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)

contains


subroutine new_mchrg_model(self, chi, rad, eta, kcn)
   type(mchrg_model_type), intent(out) :: self
   real(wp), intent(in) :: rad(:)
   real(wp), intent(in) :: chi(:)
   real(wp), intent(in) :: eta(:)
   real(wp), intent(in) :: kcn(:)

   self%rad = rad
   self%chi = chi
   self%eta = eta
   self%kcn = kcn
end subroutine new_mchrg_model


subroutine get_vrhs(self, mol, cn, xvec, dxdcn)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: cn(:)
   real(wp), intent(out) :: xvec(:)
   real(wp), intent(out), optional :: dxdcn(:)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   if (present(dxdcn)) then
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp = self%kcn(izp) / sqrt(cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*cn(iat)
         dxdcn(iat) = 0.5_wp*tmp
      end do
      dxdcn(mol%nat+1) = 0.0_wp
   else
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp = self%kcn(izp) / sqrt(cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*cn(iat)
      end do
   end if
   xvec(mol%nat+1) = mol%charge

end subroutine get_vrhs


subroutine get_amat(self, mol, amat)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp / (self%rad(izp)**2 + self%rad(jzp)**2)
         tmp = erf(sqrt(r2*gam))/sqrt(r2)
         amat(jat, iat) = amat(jat, iat) + tmp
         amat(iat, jat) = amat(iat, jat) + tmp
      end do
      tmp = self%eta(izp) + sqrt2pi / self%rad(izp)
      amat(iat, iat) = amat(iat, iat) + tmp
   end do

   amat(mol%nat+1, 1:mol%nat+1) = 1.0_wp
   amat(1:mol%nat+1, mol%nat+1) = 1.0_wp
   amat(mol%nat+1, mol%nat+1) = 0.0_wp

end subroutine get_amat

subroutine get_damat(self, mol, qvec, dadr, dadL, atrace)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(out) :: dadr(:, :, :)
   real(wp), intent(out) :: dadL(:, :, :)
   real(wp), intent(out) :: atrace(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, arg, dtmp, dG(3), dS(3, 3)

   atrace(:, :) = 0.0_wp
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
         gam = 1.0_wp/sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
         arg = gam*gam*r2
         dtmp = 2.0_wp*gam*exp(-arg)/(sqrtpi*r2)-erf(sqrt(arg))/(r2*sqrt(r2))
         dG = dtmp*vec
         dS = spread(dG, 1, 3) * spread(vec, 2, 3)
         atrace(:, iat) = +dG*qvec(jat) + atrace(:, iat)
         atrace(:, jat) = -dG*qvec(iat) + atrace(:, jat)
         dadr(:, iat, jat) = +dG*qvec(iat)
         dadr(:, jat, iat) = -dG*qvec(jat)
         dadL(:, :, jat) = +dS*qvec(iat) + dadL(:, :, jat)
         dadL(:, :, iat) = +dS*qvec(jat) + dadL(:, :, iat)
      end do
   end do

end subroutine get_damat

subroutine solve(self, mol, cn, dcndr, dcndL, energy, gradient, sigma, qvec, dqdr, dqdL)
   class(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in), contiguous :: cn(:)
   real(wp), intent(in), contiguous, optional :: dcndr(:, :, :)
   real(wp), intent(in), contiguous, optional :: dcndL(:, :, :)
   real(wp), intent(out), contiguous, optional :: qvec(:)
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)
   real(wp), intent(inout), contiguous, optional :: energy(:)
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   integer :: ic, jc, iat, ndim, info
   real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :), ainv(:, :)
   real(wp), allocatable :: dxdcn(:), atrace(:, :), dadr(:, :, :), dadL(:, :, :)
   integer, allocatable :: ipiv(:)
   logical :: grad, cpq, dcn

   ndim = mol%nat + 1

   dcn = present(dcndr) .and. present(dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn
   cpq = present(dqdr) .and. present(dqdL) .and. dcn

   allocate(amat(ndim, ndim), xvec(ndim))
   allocate(ipiv(ndim))
   if (grad.or.cpq) then
      allocate(dxdcn(ndim))
   end if

   call get_vrhs(self, mol, cn, xvec, dxdcn)
   call get_amat(self, mol, amat)

   vrhs = xvec
   ainv = amat

   call sytrf(ainv, ipiv, info=info, uplo='l')

   if (info == 0) then
      if (cpq) then
         call sytri(ainv, ipiv, info=info, uplo='l')
         if (info == 0) then
            call symv(ainv, xvec, vrhs, uplo='l')
            do ic = 1, ndim
               do jc = ic+1, ndim
                  ainv(ic, jc) = ainv(jc, ic)
               end do
            end do
         end if
      else
         call sytrs(ainv, vrhs, ipiv, info=info, uplo='l')
      end if
   end if

   if (present(qvec)) then
      qvec(:) = vrhs(:mol%nat)
   end if

   if (present(energy)) then
      call symv(amat, vrhs, xvec, alpha=0.5_wp, beta=-1.0_wp, uplo='l')
      energy(:) = energy(:) + vrhs(:mol%nat) * xvec(:mol%nat)
   end if

   if (grad.or.cpq) then
      allocate(dadr(3, mol%nat, ndim), dadL(3, 3, ndim), atrace(3, mol%nat))
      call get_damat(self, mol, vrhs, dadr, dadL, atrace)
      xvec(:) = -dxdcn * vrhs
   end if

   if (grad) then
      call gemv(dadr, vrhs, gradient, beta=1.0_wp)
      call gemv(dcndr, xvec(:mol%nat), gradient, beta=1.0_wp)
      call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dcndL, xvec(:mol%nat), sigma, beta=1.0_wp)
   end if

   if (cpq) then
      do iat = 1, mol%nat
         dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
         dadr(:, :, iat) = -dcndr(:, :, iat) * dxdcn(iat) + dadr(:, :, iat)
         dadL(:, :, iat) = -dcndL(:, :, iat) * dxdcn(iat) + dadL(:, :, iat)
      end do

      call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
      call gemm(dadL, ainv(:, :mol%nat), dqdL, alpha=-1.0_wp)
   end if

end subroutine solve


end module multicharge_model
