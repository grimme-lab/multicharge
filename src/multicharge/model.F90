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

#ifndef IK
#define IK i4
#endif

module multicharge_model
   use mctc_env, only : error_type, wp, ik => IK
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use multicharge_blas, only : gemv, symv, gemm
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_ewald, only : get_alpha
   use multicharge_lapack, only : sytrf, sytrs, sytri
   use multicharge_wignerseitz, only : wignerseitz_cell_type, new_wignerseitz_cell
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


   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

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
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cn, xvec, dxdcn) private(iat, izp, tmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp = self%kcn(izp) / sqrt(cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*cn(iat)
         dxdcn(iat) = 0.5_wp*tmp
      end do
      dxdcn(mol%nat+1) = 0.0_wp
   else
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cn, xvec) private(iat, izp, tmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp = self%kcn(izp) / sqrt(cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*cn(iat)
      end do
   end if
   xvec(mol%nat+1) = mol%charge

end subroutine get_vrhs


subroutine get_dir_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2

   call get_lattice_points(lattice, rep, .true., trans)

end subroutine get_dir_trans

subroutine get_rec_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points(rec_lat, rep, .false., trans)

end subroutine get_rec_trans


subroutine get_amat_0d(self, mol, amat)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(out) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r2, gam, tmp

   amat(:, :) = 0.0_wp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:amat) shared(mol, self) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
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

end subroutine get_amat_0d

subroutine get_amat_3d(self, mol, wsc, alpha, amat)
   type(mchrg_model_type), intent(in) :: self
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

subroutine get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
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

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace, dadr, dadL) shared(mol, self, qvec) &
   !$omp private(iat, izp, jat, jzp, gam, r2, vec, dG, dS, dtmp, arg)
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

end subroutine get_damat_0d

subroutine get_damat_3d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace)
   type(mchrg_model_type), intent(in) :: self
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

   integer :: ic, jc, iat, ndim
   logical :: grad, cpq, dcn
   real(wp) :: alpha
   integer(ik) :: info
   integer(ik), allocatable :: ipiv(:)

   real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :), ainv(:, :)
   real(wp), allocatable :: dxdcn(:), atrace(:, :), dadr(:, :, :), dadL(:, :, :)
   type(wignerseitz_cell_type) :: wsc

   ndim = mol%nat + 1
   if (any(mol%periodic)) then
      call new_wignerseitz_cell(wsc, mol)
      call get_alpha(mol%lattice, alpha)
   end if

   dcn = present(dcndr) .and. present(dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn
   cpq = present(dqdr) .and. present(dqdL) .and. dcn

   allocate(amat(ndim, ndim), xvec(ndim))
   allocate(ipiv(ndim))
   if (grad.or.cpq) then
      allocate(dxdcn(ndim))
   end if

   call get_vrhs(self, mol, cn, xvec, dxdcn)
   if (any(mol%periodic)) then
      call get_amat_3d(self, mol, wsc, alpha, amat)
   else
      call get_amat_0d(self, mol, amat)
   end if

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
      if (any(mol%periodic)) then
         call get_damat_3d(self, mol, wsc, alpha, vrhs, dadr, dadL, atrace)
      else
         call get_damat_0d(self, mol, vrhs, dadr, dadL, atrace)
      end if
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
