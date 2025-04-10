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

!> Electronegativity equlibration charge model
module multicharge_model_eeq
   use mctc_env, only: wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matdet_3x3
   use mctc_ncoord, only: new_ncoord
   use multicharge_wignerseitz, only: wignerseitz_cell_type, new_wignerseitz_cell
   use multicharge_ewald, only: get_alpha
   use multicharge_model_type, only: mchrg_model_type, get_dir_trans, get_rec_trans
   use multicharge_model_cache, only: cache_container, model_cache
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
      !> Calculate Coulomb matrix
      procedure :: get_amat_0d
      !> Calculate Coulomb matrix periodic
      procedure :: get_amat_3d
      !> Calculate Coulomb matrix derivative
      procedure :: get_damat_0d
      !> Calculate Coulomb matrix derivative periodic
      procedure :: get_damat_3d
   end type eeq_model

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: sqrt2pi = sqrt(2.0_wp/pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains

   subroutine new_eeq_model(self, mol, chi, rad, eta, kcnchi, &
      & cutoff, cn_exp, rcov, cn_max, dielectric)
      !> Electronegativity equilibration model
      type(eeq_model), intent(out) :: self
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
      !> Cutoff radius for coordination number
      real(wp), intent(in), optional :: cutoff
      !> Steepness of the CN counting function
      real(wp), intent(in), optional :: cn_exp
      !> Covalent radii for CN
      real(wp), intent(in), optional :: rcov(:)
      !> Maximum CN cutoff for CN
      real(wp), intent(in), optional :: cn_max
      !> Dielectric constant of the surrounding medium
      real(wp), intent(in), optional :: dielectric

      self%chi = chi
      self%rad = rad
      self%eta = eta
      self%kcnchi = kcnchi

      if (present(dielectric)) then
         self%dielectric = dielectric
      else
         self%dielectric = 1.0_wp
      end if

      call new_ncoord(self%ncoord, mol, "erf", cutoff=cutoff, kcn=cn_exp, &
         & rcov=rcov, cut=cn_max)

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
         tmp = self%kcnchi(izp)/sqrt(ptr%cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*ptr%cn(iat)
      end do
      xvec(mol%nat + 1) = mol%charge

   end subroutine get_xvec

   subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
      class(eeq_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out) :: dxdr(:, :, :)
      real(wp), intent(out) :: dxdL(:, :, :)
      real(wp), parameter :: reg = 1.0e-14_wp

      integer :: iat, izp
      real(wp) :: tmp
      real(wp), allocatable :: dtrans(:, :), cn(:), dcndr(:, :, :), dcndL(:, :, :)

      type(eeq_cache), pointer :: ptr

      call view(cache, ptr)

      dxdr(:, :, :) = 0.0_wp
      dxdL(:, :, :) = 0.0_wp

      ! NOTE: just remove the mol%periodic branch to restore
      if (any(mol%periodic)) then
         !$omp parallel do default(none) schedule(runtime) &
         !$omp shared(mol, self, ptr, dxdr, dxdL, dtrans, cn) &
         !$omp private(iat, izp, tmp, wsw)
         do iat = 1, mol%nat
            izp = mol%id(iat)
            wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               call self%ncoord%get_coordination_number(mol, wsc%trans(:, wsc%tridx(img, jat, iat)), cn, dcndr, dcndL)
               tmp = self%kcnchi(izp)/sqrt(cn(iat) + reg)
               dxdr(:, :, iat) = 0.5_wp*tmp*dcndr(:, :, iat)*wsw + dxdr(:, :, iat)
               dxdL(:, :, iat) = 0.5_wp*tmp*dcndL(:, :, iat)*wsw + dxdL(:, :, iat)
            end do
         end do
      else
         !$omp parallel do default(none) schedule(runtime) &
         !$omp shared(mol, self, ptr, dxdr, dxdL) &
         !$omp private(iat, izp, tmp)
         do iat = 1, mol%nat
            izp = mol%id(iat)
            tmp = self%kcnchi(izp)/sqrt(ptr%cn(iat) + reg)
            dxdr(:, :, iat) = 0.5_wp*tmp*ptr%dcndr(:, :, iat) + dxdr(:, :, iat)
            dxdL(:, :, iat) = 0.5_wp*tmp*ptr%dcndL(:, :, iat) + dxdL(:, :, iat)
         end do
      end if
   end subroutine get_xvec_derivs

   subroutine get_coulomb_matrix(self, mol, cache, amat)
      class(eeq_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(out) :: amat(:, :)

      type(eeq_cache), pointer :: ptr

      call view(cache, ptr)

      if (any(mol%periodic)) then
         call self%get_amat_3d(mol, ptr%wsc, ptr%alpha, amat)
      else
         call self%get_amat_0d(mol, amat)
      end if
   end subroutine get_coulomb_matrix

   subroutine get_amat_0d(self, mol, amat)
      class(eeq_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(out) :: amat(:, :)

      integer :: iat, jat, izp, jzp
      real(wp) :: vec(3), r2, gam, tmp

      amat(:, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(amat, mol, self) &
      !$omp private(iat, izp, jat, jzp, gam, vec, r2, tmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            gam = 1.0_wp/(self%rad(izp)**2 + self%rad(jzp)**2)
            tmp = erf(sqrt(r2*gam))/(sqrt(r2)*self%dielectric)
            !$omp atomic
            amat(jat, iat) = amat(jat, iat) + tmp
            !$omp atomic
            amat(iat, jat) = amat(iat, jat) + tmp
         end do
         tmp = self%eta(izp) + sqrt2pi/self%rad(izp)
         !$omp atomic
         amat(iat, iat) = amat(iat, iat) + tmp
      end do

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
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            gam = 1.0_wp/sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
            wsw = 1.0_wp/real(wsc%nimg(jat, iat), wp)
            do img = 1, wsc%nimg(jat, iat)
               vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
               call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
               call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
               amat(jat, iat) = amat(jat, iat) + (dtmp + rtmp)*wsw
               amat(iat, jat) = amat(iat, jat) + (dtmp + rtmp)*wsw
            end do
         end do

         gam = 1.0_wp/sqrt(2.0_wp*self%rad(izp)**2)
         wsw = 1.0_wp/real(wsc%nimg(iat, iat), wp)
         do img = 1, wsc%nimg(iat, iat)
            vec = wsc%trans(:, wsc%tridx(img, iat, iat))
            call get_amat_dir_3d(vec, gam, alpha, dtrans, dtmp)
            call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
            amat(iat, iat) = amat(iat, iat) + (dtmp + rtmp)*wsw
         end do

         dtmp = self%eta(izp) + sqrt2pi/self%rad(izp) - 2*alpha/sqrtpi
         amat(iat, iat) = amat(iat, iat) + dtmp
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
      class(eeq_model), intent(in) :: self
      type(structure_type), intent(in) :: mol
      type(cache_container), intent(inout) :: cache
      real(wp), intent(in) :: qvec(:)
      real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

      type(eeq_cache), pointer :: ptr

      call view(cache, ptr)

      if (any(mol%periodic)) then
         call self%get_damat_3d(mol, ptr%wsc, ptr%alpha, qvec, dadr, dadL, atrace)
      else
         call self%get_damat_0d(mol, qvec, dadr, dadL, atrace)
      end if
   end subroutine get_coulomb_derivs

   subroutine get_damat_0d(self, mol, qvec, dadr, dadL, atrace)
      class(eeq_model), intent(in) :: self
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
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat)
            r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
            gam = 1.0_wp/sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
            arg = gam*gam*r2
            dtmp = 2.0_wp*gam*exp(-arg)/(sqrtpi*r2*self%dielectric) &
               & - erf(sqrt(arg))/(r2*sqrt(r2)*self%dielectric)
            dG = dtmp*vec
            dS = spread(dG, 1, 3)*spread(vec, 2, 3)
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
      class(eeq_model), intent(in) :: self
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

   ! NOTE: the following is basically identical to tblite versions of this pattern

   !> Inspect cache and reallocate it in case of type mismatch
   subroutine taint(cache, ptr)
      !> Instance of the cache
      type(cache_container), target, intent(inout) :: cache
      !> Reference to the cache
      type(eeq_cache), pointer, intent(out) :: ptr

      if (allocated(cache%raw)) then
         call view(cache, ptr)
         if (associated(ptr)) return
         deallocate (cache%raw)
      end if

      if (.not. allocated(cache%raw)) then
         block
            type(eeq_cache), allocatable :: tmp
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
      type(eeq_cache), pointer, intent(out) :: ptr
      nullify (ptr)
      select type (target => cache%raw)
      type is (eeq_cache)
         ptr => target
      end select
   end subroutine view

end module multicharge_model_eeq
