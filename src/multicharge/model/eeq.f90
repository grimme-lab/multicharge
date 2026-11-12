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
   use mctc_cutoff, only: get_lattice_points, get_periodic_axes
   use mctc_env, only: error_type, fatal_error, wp
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: crossprod, outerprod
   use mctc_ncoord, only: new_ncoord, cn_count
   use mctc_wignerseitz, only: wignerseitz_cell, new_wignerseitz_cell
   use multicharge_blas, only: gemm
   use multicharge_ewald, only: get_alpha, get_ewald_dir_trans, get_ewald_rec_trans, &
      & erfc_exp_2d, fast_f012_1d, fast_k012_1d
   use multicharge_model_cache, only: cache_container, model_cache
   use multicharge_model_hessian, only: coulomb_kernel, radial_hessian, &
      & add_pair_hessian, add_pair_block
   use multicharge_model_type, only: mchrg_model_type
   implicit none
   private

   public :: eeq_model, new_eeq_model

   !> Reusable data container for EEQ model evaluations
   type, extends(model_cache), public :: eeq_cache

      !> Number of periodic dimensions
      integer :: ndim = 0

      !> Cell length (1D), area (2D), or volume (3D)
      real(wp) :: vol = 0.0_wp

      !> Wire direction (1D) or slab normal (2D)
      real(wp) :: axis_vec(3) = 0.0_wp

      !> Direct-space Ewald translations
      real(wp), allocatable :: dtrans(:, :)

      !> Nonzero reciprocal-space Ewald vectors
      real(wp), allocatable :: rtrans(:, :)

      !> Reciprocal-space Ewald weights (3D)
      real(wp), allocatable :: rec_weight(:)

      !> Reciprocal-space Ewald strain tensors (3D)
      real(wp), allocatable :: rec_strain(:, :, :)

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
      if (.not. all(mol%periodic) .and. abs(mol%charge) > eps) then
         call fatal_error(error, &
            & "Charged 1D and 2D periodic systems are not supported by EEQ")
         return
      end if

      call new_wignerseitz_cell(ptr%wsc, mol, error)
      if (allocated(error)) return

      call setup_ewald(self, mol, ptr, error)
      if (allocated(error)) return
   end if

end subroutine update


!> Set up Ewald splitting parameter, cell geometry, and lattice translations
subroutine setup_ewald(self, mol, cache, error)

   !> Electronegativity equilibration model
   class(eeq_model), intent(in) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> EEQ model cache
   type(eeq_cache), intent(inout) :: cache

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: p_axes(3)
   real(wp) :: gam_min, max_rij, n_vec(3)

   call get_alpha(mol%lattice, cache%alpha, error, mol%periodic)
   if (allocated(error)) return

   cache%vol = 0.0_wp
   cache%axis_vec = 0.0_wp

   call get_periodic_axes(mol%periodic, cache%ndim, p_axes)
   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(self%rad(mol%id)))
   max_rij = get_max_wsc_disp(mol, cache%wsc)

   call get_ewald_dir_trans(mol%periodic, mol%lattice, cache%alpha, &
      & gam_min, cache%dtrans, error, max_rij)
   if (allocated(error)) return
   call get_ewald_rec_trans(mol%periodic, mol%lattice, cache%alpha, &
      & cache%rtrans, error, cache%vol)
   if (allocated(error)) return

   select case (cache%ndim)
   case (1)
      cache%axis_vec = mol%lattice(:, p_axes(1)) / cache%vol
   case (2)
      n_vec = crossprod(mol%lattice(:, p_axes(1)), mol%lattice(:, p_axes(2)))
      cache%axis_vec = n_vec / norm2(n_vec)
   case (3)
      call get_rec_weights(cache%vol, cache%alpha, cache%rtrans, cache%rec_weight, &
         & cache%rec_strain)
   case default
      continue
   end select

end subroutine setup_ewald


!> Compute maximum interatomic pair displacement within the Wigner-Seitz cell
pure function get_max_wsc_disp(mol, wsc) result(max_rij)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell instance
   type(wignerseitz_cell), intent(in) :: wsc

   !> Maximum pair displacement norm
   real(wp) :: max_rij

   integer :: iat, jat, img
   real(wp) :: vec(3)

   max_rij = 0.0_wp
   do iat = 1, mol%nat
      do jat = 1, iat - 1
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) &
               & + wsc%trans(:, wsc%tridx(img, jat, iat))
            max_rij = max(max_rij, norm2(vec))
         end do
      end do
   end do

end function get_max_wsc_disp


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
   real(wp), allocatable :: dEdcn(:), dcndr(:, :), work(:, :), trans(:, :)
   real(wp), allocatable :: diagonal(:, :, :)

   call view(cache, ptr)
   if (.not. associated(ptr)) then
      call fatal_error(error, "Incompatible EEQ response cache")
      return
   end if

   nat = mol%nat
   allocate(dEdcn(nat), work(3 * nat, nat))
   dcndr = reshape(ptr%dcndr, [3 * nat, nat])
   do iat = 1, nat
      izp = mol%id(iat)
      c = ptr%cn(iat)
      dEdcn(iat) = adj(iat) * self%kcnchi(izp) * (0.5_wp * c + reg) &
         & / (c + reg)**1.5_wp
      weight = -adj(iat) * self%kcnchi(izp) * (0.25_wp * c + reg) &
         & / (c + reg)**2.5_wp
      work(:, iat) = weight * dcndr(:, iat)
   end do
   call get_lattice_points(mol%periodic, mol%lattice, self%ncoord%cutoff, &
      & trans, error)
   if (allocated(error)) return

   call self%ncoord%add_coordination_number_hessian(mol, trans, dEdcn, hessian)
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
         call add_pair_hessian(iat, jat, vec, r, weight * dr, weight * drr, &
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

   ! Each pair owns its off-diagonal blocks. Only the diagonal blocks are
   ! reduced, avoiding a full (3N, 3N) Hessian copy for every OpenMP thread.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(self, mol, ptr, qvec, adj, hessian, npair) &
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
            if (r < eps) cycle

            call coulomb_kernel(r, gam, val, dr, drr)
            call coulomb_kernel(r, ptr%alpha, val, ar, arr)
            box = box + wsw * radial_hessian(vec, r, dr - ar, drr - arr)
         end do
         select case (ptr%ndim)
         case (1)
            call get_hess_rec_1d(rij, ptr%alpha, ptr%vol, ptr%axis_vec, &
               & ptr%rtrans, rbox)
         case (2)
            call get_hess_rec_2d(rij, ptr%alpha, ptr%vol, ptr%axis_vec, &
               & ptr%rtrans, rbox)
         case (3)
            call get_hess_rec_3d(rij, ptr%rtrans, ptr%rec_weight, rbox)
         case default
            rbox = 0.0_wp
         end select
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


!> Reciprocal-space and G=0 sheet pair Hessian block for 2D Ewald
pure subroutine get_hess_rec_2d(rij, alpha, area, n_hat, trans, rbox)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic surface area
   real(wp), intent(in) :: area

   !> Unit normal vector to the periodic plane
   real(wp), intent(in) :: n_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space Hessian block
   real(wp), intent(out) :: rbox(3, 3)

   integer :: itr
   real(wp) :: z, pi_area, inv_alp2, z2_alp2, g_vec(3), g2, g_norm, gv, cos_gv, sin_gv
   real(wp) :: ep, em, fp, fm, gauss, pref, gauss_z
   real(wp) :: gg_mat(3, 3), nn_mat(3, 3), ng_sym(3, 3)

   rbox = 0.0_wp
   z = dot_product(rij, n_hat)
   pi_area = pi / area
   inv_alp2 = 1.0_wp / (alpha * alpha)
   z2_alp2 = (alpha * z)**2
   nn_mat = outerprod(n_hat, n_hat)

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle

      g_norm = sqrt(g2)
      gv = dot_product(g_vec, rij)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      ep = erfc_exp_2d(g_norm, +z, alpha)
      em = erfc_exp_2d(g_norm, -z, alpha)
      fp = ep + em
      fm = ep - em
      gauss = exp(-0.25_wp * g2 * inv_alp2 - z2_alp2)
      pref = pi_area / g_norm
      gg_mat = outerprod(g_vec, g_vec)
      ng_sym = outerprod(g_vec, n_hat) + outerprod(n_hat, g_vec)
      rbox = rbox - pref * cos_gv * fp * gg_mat &
         & - pi_area * sin_gv * fm * ng_sym &
         & + pi_area * cos_gv &
         & * (g_norm * fp - (4.0_wp / sqrtpi) * alpha * gauss) * nn_mat
   end do

   gauss_z = exp(-z2_alp2)
   rbox = rbox - (4.0_wp * sqrtpi * alpha / area) * gauss_z * nn_mat

end subroutine get_hess_rec_2d


!> Reciprocal-space and G=0 wire pair Hessian block for 1D Ewald
pure subroutine get_hess_rec_1d(rij, alpha, len1, u_hat, trans, rbox)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic wire length
   real(wp), intent(in) :: len1

   !> Unit vector along the periodic wire axis
   real(wp), intent(in) :: u_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space Hessian block
   real(wp), intent(out) :: rbox(3, 3)

   integer :: itr
   real(wp) :: r_perp(3), rho2, alp2, inv_alp2, inv_len, x, y, g_vec(3), g2, gv
   real(wp) :: cos_gv, sin_gv, k0, k1, k2, f0, f1, f2
   real(wp) :: p_perp(3, 3), gg_mat(3, 3), rg_sym(3, 3), rr_mat(3, 3)

   rbox = 0.0_wp
   p_perp = unity - outerprod(u_hat, u_hat)
   r_perp = rij - dot_product(rij, u_hat) * u_hat
   rho2 = dot_product(r_perp, r_perp)
   alp2 = alpha * alpha
   inv_alp2 = 1.0_wp / alp2
   inv_len = 1.0_wp / len1
   y = alp2 * rho2
   rr_mat = outerprod(r_perp, r_perp)

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle

      x = 0.25_wp * g2 * inv_alp2
      gv = dot_product(g_vec, rij)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      call fast_k012_1d(x, y, k0, k1, k2)
      gg_mat = outerprod(g_vec, g_vec)
      rg_sym = outerprod(g_vec, r_perp) + outerprod(r_perp, g_vec)
      rbox = rbox + inv_len * (-cos_gv * k0 * gg_mat &
         & + 2.0_wp * alp2 * sin_gv * k1 * rg_sym &
         & + cos_gv * (-2.0_wp * alp2 * k1 * p_perp &
         & + 4.0_wp * alp2 * alp2 * k2 * rr_mat))
   end do

   call fast_f012_1d(y, f0, f1, f2)
   rbox = rbox + inv_len * (-2.0_wp * alp2 * f1 * p_perp &
      & - 4.0_wp * alp2 * alp2 * f2 * rr_mat)

end subroutine get_hess_rec_1d


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

   integer :: iat, jat, izp, jzp, img, ipair, npair
   real(wp) :: rij(3), vec(3), gam, wsw, dtmp, rtmp, tmp

   ! Every matrix element belongs to one pair; no thread-local matrix is needed.
   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel do default(none) schedule(runtime) &
   !$omp& shared(amat, mol, self, ptr, npair) &
   !$omp& private(iat, jat, izp, jzp, img, gam, wsw, rij, vec, dtmp, rtmp, tmp)
   do ipair = 1, npair
      iat = int(0.5_wp * (sqrt(8.0_wp * real(ipair, wp) + 1.0_wp) - 1.0_wp))
      if (iat * (iat + 1) / 2 < ipair) iat = iat + 1
      jat = ipair - iat * (iat - 1) / 2
      izp = mol%id(iat)
      jzp = mol%id(jat)
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      if (iat == jat) then
         vec = 0.0_wp
         call get_amat_dir(vec, gam, ptr%alpha, ptr%dtrans, dtmp)
         select case (ptr%ndim)
         case (1)
            call get_amat_rec_1d(vec, ptr%alpha, ptr%vol, ptr%axis_vec, &
               & ptr%rtrans, rtmp)
         case (2)
            call get_amat_rec_2d(vec, ptr%alpha, ptr%vol, ptr%axis_vec, &
               & ptr%rtrans, rtmp)
         case (3)
            call get_amat_rec_3d(vec, ptr%rtrans, ptr%rec_weight, rtmp)
         case default
            rtmp = 0.0_wp
         end select
         tmp = dtmp + rtmp + self%eta(izp) + sqrt2pi / self%rad(izp) &
            & - 2.0_wp * ptr%alpha / sqrtpi
      else
         rij = mol%xyz(:, jat) - mol%xyz(:, iat)
         wsw = 1.0_wp / real(ptr%wsc%nimg(jat, iat), wp)
         tmp = 0.0_wp
         do img = 1, ptr%wsc%nimg(jat, iat)
            vec = rij + ptr%wsc%trans(:, ptr%wsc%tridx(img, jat, iat))
            call get_amat_dir(vec, gam, ptr%alpha, ptr%dtrans, dtmp)
            select case (ptr%ndim)
            case (1)
               call get_amat_rec_1d(vec, ptr%alpha, ptr%vol, ptr%axis_vec, &
                  & ptr%rtrans, rtmp)
            case (2)
               call get_amat_rec_2d(vec, ptr%alpha, ptr%vol, ptr%axis_vec, &
                  & ptr%rtrans, rtmp)
            case (3)
               call get_amat_rec_3d(vec, ptr%rtrans, ptr%rec_weight, rtmp)
            case default
               rtmp = 0.0_wp
            end select
            tmp = tmp + dtmp + rtmp
         end do
         tmp = tmp * wsw
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

      tmp = (erf(gam * r1) - erf(alp * r1)) / r1
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


!> Reciprocal-space and G=0 sheet contribution for 2D Ewald
pure subroutine get_amat_rec_2d(rij, alpha, area, n_hat, trans, amat)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic surface area
   real(wp), intent(in) :: area

   !> Unit normal vector to the periodic plane
   real(wp), intent(in) :: n_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space Coulomb matrix element contribution
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: z, pi_area, g_vec(3), g_norm, gv, ep, em

   z = dot_product(rij, n_hat)
   pi_area = pi / area
   amat = 0.0_wp
   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g_norm = norm2(g_vec)
      if (g_norm < eps) cycle

      gv = dot_product(rij, g_vec)
      ep = erfc_exp_2d(g_norm, +z, alpha)
      em = erfc_exp_2d(g_norm, -z, alpha)
      amat = amat + (pi_area / g_norm) * cos(gv) * (ep + em)
   end do
   amat = amat - 2.0_wp * pi_area &
      & * (z * erf(alpha * z) + exp(-(alpha * z)**2) / (alpha * sqrtpi))

end subroutine get_amat_rec_2d


!> Reciprocal-space and G=0 wire contribution for 1D Ewald
pure subroutine get_amat_rec_1d(rij, alpha, len1, u_hat, trans, amat)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic wire length
   real(wp), intent(in) :: len1

   !> Unit vector along the periodic wire axis
   real(wp), intent(in) :: u_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space Coulomb matrix element contribution
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: r_perp(3), rho2, alp2, inv_alp2, x, y, g_vec(3), g2, gv
   real(wp) :: k0, k1, k2, f0, f1, f2

   amat = 0.0_wp
   r_perp = rij - dot_product(rij, u_hat) * u_hat
   rho2 = dot_product(r_perp, r_perp)
   alp2 = alpha * alpha
   inv_alp2 = 1.0_wp / alp2
   y = alp2 * rho2

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle

      x = 0.25_wp * g2 * inv_alp2
      gv = dot_product(rij, g_vec)
      call fast_k012_1d(x, y, k0, k1, k2)
      amat = amat + cos(gv) * k0
   end do
   call fast_f012_1d(y, f0, f1, f2)
   amat = (amat - f0 + 2.0_wp * log(alpha)) / len1

end subroutine get_amat_rec_1d


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
   real(wp) :: vec(3), r1, r2, inv_r, inv_r2, gam, arg, dtmp, dG(3), dS(3, 3)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   npair = mol%nat * (mol%nat + 1) / 2

   ! Pair derivatives have unique owners. Only atom-wise traces and strain
   ! derivatives need reduction, with O(N) storage per thread.
   !$omp parallel default(none) shared(dadr, mol, self, qvec, npair) &
   !$omp& private(iat, jat, izp, jzp, ipair, gam, r1, r2, inv_r, inv_r2, &
   !$omp& vec, dG, dS, dtmp, arg) &
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
      r2 = dot_product(vec, vec)
      r1 = sqrt(r2)
      inv_r = 1.0_wp / r1
      inv_r2 = inv_r * inv_r
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      arg = gam * r1
      dtmp = ((2.0_wp / sqrtpi) * gam * exp(-arg**2) - erf(arg) * inv_r) * inv_r2
      dG = dtmp * vec
      dS = outerprod(vec, dG)
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


!> Evaluate 1D, 2D, or 3D Ewald Coulomb matrix derivatives
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
   real(wp) :: gam, wsw, rij(3), vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)

   atrace(:, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp

   npair = mol%nat * (mol%nat + 1) / 2
   !$omp parallel default(none) &
   !$omp& shared(mol, self, ptr, qvec, dadr, npair) &
   !$omp& private(iat, jat, izp, jzp, img, ipair, gam, wsw, rij, vec, dG, dS, &
   !$omp& dGr, dSr, dGd, dSd) &
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
      dG(:) = 0.0_wp
      dS(:, :) = 0.0_wp
      gam = 1.0_wp / sqrt(self%rad(izp)**2 + self%rad(jzp)**2)
      if (iat == jat) then
         vec = 0.0_wp
         call get_damat_dir(vec, gam, ptr%alpha, ptr%dtrans, dGd, dSd)
         select case (ptr%ndim)
         case (1)
            call get_damat_rec_1d(vec, ptr%alpha, ptr%vol, ptr%axis_vec, &
               & ptr%rtrans, dGr, dSr)
         case (2)
            call get_damat_rec_2d(vec, ptr%alpha, ptr%vol, ptr%axis_vec, &
               & ptr%rtrans, dGr, dSr)
         case (3)
            call get_damat_rec_3d(vec, ptr%rtrans, ptr%rec_weight, &
               & ptr%rec_strain, dGr, dSr)
         case default
            dGr = 0.0_wp
            dSr = 0.0_wp
         end select
         dS = dSd + dSr
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(iat)
      else
         rij = mol%xyz(:, jat) - mol%xyz(:, iat)
         wsw = 1.0_wp / real(ptr%wsc%nimg(jat, iat), wp)
         do img = 1, ptr%wsc%nimg(jat, iat)
            vec = rij + ptr%wsc%trans(:, ptr%wsc%tridx(img, jat, iat))
            call get_damat_dir(vec, gam, ptr%alpha, ptr%dtrans, dGd, dSd)
            select case (ptr%ndim)
            case (1)
               call get_damat_rec_1d(vec, ptr%alpha, ptr%vol, ptr%axis_vec, &
                  & ptr%rtrans, dGr, dSr)
            case (2)
               call get_damat_rec_2d(vec, ptr%alpha, ptr%vol, ptr%axis_vec, &
                  & ptr%rtrans, dGr, dSr)
            case (3)
               call get_damat_rec_3d(vec, ptr%rtrans, ptr%rec_weight, &
                  & ptr%rec_strain, dGr, dSr)
            case default
               dGr = 0.0_wp
               dSr = 0.0_wp
            end select
            dG = dG + dGd + dGr
            dS = dS + dSd + dSr
         end do
         dG = dG * wsw
         dS = dS * wsw
         atrace(:, iat) = atrace(:, iat) - dG * qvec(jat)
         atrace(:, jat) = atrace(:, jat) + dG * qvec(iat)
         dadr(:, iat, jat) = -dG * qvec(iat)
         dadr(:, jat, iat) = +dG * qvec(jat)
         dadL(:, :, jat) = dadL(:, :, jat) + dS * qvec(iat)
         dadL(:, :, iat) = dadL(:, :, iat) + dS * qvec(jat)
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
   real(wp) :: vec(3), r1, r2, inv_r, inv_r2, gtmp, atmp, dtmp
   real(wp) :: gam2, alp2, fac_g, fac_a

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   gam2 = gam * gam
   alp2 = alp * alp
   fac_g = (2.0_wp / sqrtpi) * gam
   fac_a = (2.0_wp / sqrtpi) * alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle

      r2 = r1 * r1
      inv_r = 1.0_wp / r1
      inv_r2 = inv_r * inv_r
      gtmp = fac_g * exp(-r2 * gam2) - erf(r1 * gam) * inv_r
      atmp = erf(r1 * alp) * inv_r - fac_a * exp(-r2 * alp2)
      dtmp = (gtmp + atmp) * inv_r2
      dg(:) = dg + dtmp * vec
      ds(:, :) = ds + dtmp * outerprod(vec, vec)
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


!> Reciprocal-space and G=0 sheet derivatives for 2D Ewald
pure subroutine get_damat_rec_2d(rij, alpha, area, n_hat, trans, dg, ds)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic surface area
   real(wp), intent(in) :: area

   !> Unit normal vector to the periodic plane
   real(wp), intent(in) :: n_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Position derivative contribution
   real(wp), intent(out) :: dg(3)

   !> Strain derivative contribution
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: z, pi_area, inv_alp2, z2_alp2, fac_dg, r_perp(3), g_vec(3)
   real(wp) :: g_norm, inv_g, gv, cos_gv, sin_gv
   real(wp) :: ep, em, f_gz, df_dg, df_dz, pref, phi_0, dphi0_dz
   real(wp) :: i_par(3, 3), nn(3, 3), outer_gg(3, 3), outer_rg(3, 3)

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   z = dot_product(rij, n_hat)
   pi_area = pi / area
   inv_alp2 = 1.0_wp / alpha**2
   z2_alp2 = (alpha * z)**2
   fac_dg = 2.0_wp / (alpha * sqrtpi)
   r_perp = z * n_hat
   nn = outerprod(n_hat, n_hat)
   i_par = unity - nn

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g_norm = norm2(g_vec)
      if (g_norm < eps) cycle

      inv_g = 1.0_wp / g_norm
      gv = dot_product(rij, g_vec)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      ep = erfc_exp_2d(g_norm, +z, alpha)
      em = erfc_exp_2d(g_norm, -z, alpha)
      f_gz = ep + em
      df_dg = z * (ep - em) - fac_dg * exp(-0.25_wp * g_norm**2 * inv_alp2 - z2_alp2)
      df_dz = g_norm * (ep - em)
      pref = pi_area * inv_g

      dg(:) = dg + pref * (-sin_gv * f_gz * g_vec + cos_gv * df_dz * n_hat)

      outer_gg = outerprod(g_vec, g_vec)
      outer_rg = outerprod(r_perp, g_vec) + outerprod(g_vec, r_perp)
      ds(:, :) = ds + pref * cos_gv * (-f_gz * i_par &
         & + (f_gz * inv_g - df_dg) * inv_g * outer_gg + z * df_dz * nn) &
         & - pref * sin_gv * f_gz * outer_rg
   end do

   phi_0 = -2.0_wp * pi_area &
      & * (z * erf(alpha * z) + 0.5_wp * fac_dg * exp(-z2_alp2))
   dphi0_dz = -2.0_wp * pi_area * erf(alpha * z)
   dg(:) = dg + dphi0_dz * n_hat
   ds(:, :) = ds - phi_0 * i_par + dphi0_dz * z * nn

end subroutine get_damat_rec_2d


!> Reciprocal-space and G=0 wire derivatives for 1D Ewald
pure subroutine get_damat_rec_1d(rij, alpha, len1, u_hat, trans, dg, ds)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic wire length
   real(wp), intent(in) :: len1

   !> Unit vector along the periodic wire axis
   real(wp), intent(in) :: u_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Position derivative contribution
   real(wp), intent(out) :: dg(3)

   !> Strain derivative contribution
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: r_perp(3), rho2, alp2, inv_alp2, inv_len, x, y, g_vec(3), g2, gv
   real(wp) :: cos_gv, sin_gv, k0, k1, k2, f0, f1, f2, exp_xy
   real(wp) :: p_par(3, 3), outer_rr(3, 3), outer_rg(3, 3)

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   r_perp = rij - dot_product(rij, u_hat) * u_hat
   rho2 = dot_product(r_perp, r_perp)
   alp2 = alpha * alpha
   inv_alp2 = 1.0_wp / alp2
   inv_len = 1.0_wp / len1
   y = alp2 * rho2
   p_par = outerprod(u_hat, u_hat)
   outer_rr = outerprod(r_perp, r_perp)

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle

      x = 0.25_wp * g2 * inv_alp2
      gv = dot_product(rij, g_vec)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      call fast_k012_1d(x, y, k0, k1, k2)
      exp_xy = exp(-x - y)

      dg(:) = dg - inv_len * (sin_gv * k0 * g_vec &
         & + 2.0_wp * alp2 * cos_gv * k1 * r_perp)

      outer_rg = outerprod(g_vec, r_perp) + outerprod(r_perp, g_vec)
      ds(:, :) = ds + inv_len * (cos_gv * (2.0_wp * exp_xy + 2.0_wp * y * k1 - k0) &
         & * p_par - 2.0_wp * alp2 * cos_gv * k1 * outer_rr &
         & - sin_gv * k0 * outer_rg)
   end do

   call fast_f012_1d(y, f0, f1, f2)
   dg(:) = dg - (2.0_wp * alp2 * f1 * inv_len) * r_perp
   ds(:, :) = ds + inv_len * ((f0 - 2.0_wp * log(alpha)) * p_par &
      & - 2.0_wp * alp2 * f1 * outer_rr)

end subroutine get_damat_rec_1d


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
   real(wp) :: g2, inv_g2, inv_alp2, fac, vec(3)

   fac = 4.0_wp * pi / vol
   inv_alp2 = 1.0_wp / alpha**2
   allocate(weight(size(trans, 2)), source=0.0_wp)
   if (present(strain)) then
      allocate(strain(3, 3, size(trans, 2)), source=0.0_wp)
   end if

   do itr = 1, size(trans, 2)
      vec = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle

      inv_g2 = 1.0_wp / g2
      weight(itr) = fac * exp(-0.25_wp * g2 * inv_alp2) * inv_g2
      if (present(strain)) then
         strain(:, :, itr) = weight(itr) * ((2.0_wp * inv_g2 + 0.5_wp * inv_alp2) &
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
