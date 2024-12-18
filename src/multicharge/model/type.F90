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

!> @file multicharge/model/type.f90
!> Provides a general base class for the charge models

#ifndef IK
#define IK i4
#endif

!> General charge model
module multicharge_model_type

   use iso_fortran_env, only: output_unit

   use mctc_env, only: wp, ik => IK
   use mctc_io, only: structure_type
   use mctc_io_constants, only: pi
   use mctc_io_math, only: matinv_3x3
   use mctc_cutoff, only: get_lattice_points
   use mctc_ncoord, only: ncoord_type
   use multicharge_blas, only: gemv, symv, gemm
   use multicharge_ewald, only: get_alpha
   use multicharge_lapack, only: sytrf, sytrs, sytri
   use multicharge_wignerseitz, only: wignerseitz_cell_type, new_wignerseitz_cell
   use multicharge_model_cache, only: model_cache, cache_container
   implicit none
   private

   public :: mchrg_model_type, get_dir_trans, get_rec_trans

   !> Abstract multicharge model type
   type, abstract :: mchrg_model_type
      !> Electronegativity
      real(wp), allocatable :: chi(:)
      !> Charge width
      real(wp), allocatable :: rad(:)
      !> Chemical hardness
      real(wp), allocatable :: eta(:)
      !> CN scaling factor for electronegativity
      real(wp), allocatable :: kcnchi(:)
      !> Local charge scaling factor for electronegativity
      real(wp), allocatable :: kqchi(:)
      !> Local charge scaling factor for chemical hardness
      real(wp), allocatable :: kqeta(:)
      !> CN scaling factor for charge width
      real(wp), allocatable :: kcnrad
      !> Dielectric constant of the surrounding medium
      real(wp), allocatable :: dielectric
      !> Coordination number
      class(ncoord_type), allocatable :: ncoord
      !> Electronegativity weighted CN for local charge
      class(ncoord_type), allocatable :: ncoord_en
   contains
      !> Solve linear equations for the charge model
      procedure :: solve
      !> Calculate local charges from electronegativity weighted CN
      procedure :: local_charge
      !> Update cache
      procedure(update), deferred :: update
      !> Calculate right-hand side (electronegativity)
      procedure(get_xvec), deferred :: get_xvec
      !> Calculate xvec Gradients
      procedure(get_xvec_derivs), deferred :: get_xvec_derivs
      !> Calculate Coulomb matrix
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix
      !> Calculate Coulomb matrix derivatives
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs
   end type mchrg_model_type

   abstract interface
      subroutine update(self, mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(in) :: cn(:)
         real(wp), intent(in), optional :: qloc(:)
         real(wp), intent(in), optional :: dcndr(:, :, :)
         real(wp), intent(in), optional :: dcndL(:, :, :)
         real(wp), intent(in), optional :: dqlocdr(:, :, :)
         real(wp), intent(in), optional :: dqlocdL(:, :, :)
      end subroutine update

      subroutine get_coulomb_matrix(self, mol, cache, amat)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out) :: amat(:, :)
      end subroutine get_coulomb_matrix

      subroutine get_coulomb_derivs(self, mol, cache, vrhs, dadr, dadL, atrace)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(in) :: vrhs(:)
         real(wp), intent(out) :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
      end subroutine get_coulomb_derivs

      subroutine get_xvec(self, mol, cache, xvec)
         import :: mchrg_model_type, cache_container, structure_type, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out) :: xvec(:)
      end subroutine get_xvec

      subroutine get_xvec_derivs(self, mol, cache, dxdr, dxdL)
         import :: mchrg_model_type, structure_type, cache_container, wp
         class(mchrg_model_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(cache_container), intent(inout) :: cache
         real(wp), intent(out) :: dxdr(:, :, :)
         real(wp), intent(out) :: dxdL(:, :, :)
      end subroutine get_xvec_derivs

      !subroutine get_amat_0d(self, mol, amat, cn, qloc, cmat)
      !   import :: mchrg_model_type, cache_container, structure_type, wp
      !   class(mchrg_model_type), intent(in) :: self
      !   type(structure_type), intent(in) :: mol
      !   real(wp), intent(out) :: amat(:, :)
      !   real(wp), intent(in), optional :: cn(:)
      !   real(wp), intent(in), optional :: qloc(:)
      !   real(wp), intent(in), optional :: cmat(:, :)
      !end subroutine get_amat_0d

      !subroutine get_amat_3d(self, mol, cache, wsc, alpha, amat)
      !   import :: mchrg_model_type, cache_container, structure_type, &
      !      & wignerseitz_cell_type, wp
      !   class(mchrg_model_type), intent(in) :: self
      !   type(structure_type), intent(in) :: mol
      !   type(cache_container), intent(inout) :: cache
      !   type(wignerseitz_cell_type), intent(in) :: wsc
      !   real(wp), intent(in) :: alpha
      !   real(wp), intent(out) :: amat(:, :)
      !end subroutine get_amat_3d

      !subroutine get_damat_0d(self, mol, cn, qloc, qvec, dcndr, dcndL, &
      !   & dqlocdr, dqlocdL, dadr, dadL, atrace)
      !   import :: mchrg_model_type, structure_type, wp
      !   class(mchrg_model_type), intent(in) :: self
      !   type(structure_type), intent(in) :: mol
      !   real(wp), intent(in) :: qvec(:)
      !   real(wp), intent(out) :: dadr(:, :, :)
      !   real(wp), intent(out) :: dadL(:, :, :)
      !   real(wp), intent(out) :: atrace(:, :)
      !   real(wp), intent(in), optional :: cn(:)
      !   real(wp), intent(in), optional :: qloc(:)
      !   real(wp), intent(in), optional :: dcndr(:, :, :)
      !   real(wp), intent(in), optional :: dcndL(:, :, :)
      !   real(wp), intent(in), optional :: dqlocdr(:, :, :)
      !   real(wp), intent(in), optional :: dqlocdL(:, :, :)
      !end subroutine get_damat_0d

      !subroutine get_damat_3d(self, mol, wsc, alpha, qvec, dadr, dadL, atrace)
      !   import :: mchrg_model_type, structure_type, &
      !      & wignerseitz_cell_type, wp
      !   class(mchrg_model_type), intent(in) :: self
      !   type(structure_type), intent(in) :: mol
      !   type(wignerseitz_cell_type), intent(in) :: wsc
      !   real(wp), intent(in) :: alpha
      !   real(wp), intent(in) :: qvec(:)
      !   real(wp), intent(out) :: dadr(:, :, :)
      !   real(wp), intent(out) :: dadL(:, :, :)
      !   real(wp), intent(out) :: atrace(:, :)
      !end subroutine get_damat_3d

   end interface

   real(wp), parameter :: twopi = 2*pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains

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

   subroutine solve(self, mol, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, &
      & energy, gradient, sigma, qvec, dqdr, dqdL)
      class(mchrg_model_type), intent(in) :: self
      type(structure_type), intent(in) :: mol
      real(wp), intent(inout), contiguous :: cn(:)
      real(wp), intent(inout), contiguous :: qloc(:)
      real(wp), intent(in), contiguous, optional :: dcndr(:, :, :)
      real(wp), intent(in), contiguous, optional :: dcndL(:, :, :)
      real(wp), intent(in), contiguous, optional :: dqlocdr(:, :, :)
      real(wp), intent(in), contiguous, optional :: dqlocdL(:, :, :)
      real(wp), intent(out), contiguous, optional :: qvec(:)
      real(wp), intent(inout), contiguous, optional :: energy(:)
      real(wp), intent(inout), contiguous, optional :: gradient(:, :)
      real(wp), intent(inout), contiguous, optional :: sigma(:, :)
      real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)
      real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)

      integer :: ic, jc, iat, ndim
      logical :: grad, cpq, dcn
      integer(ik) :: info
      integer(ik), allocatable :: ipiv(:)

      !> Variables for solving ES equation
      real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :)
      real(wp), allocatable :: ainv(:, :)
      !> Gradients
      real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
      real(wp), allocatable :: dxdr(:, :, :), dxdL(:, :, :)
      type(cache_container), allocatable :: cache
      real(wp), allocatable :: trans(:, :)

      !> Calculate gradient if the respective arrays are present
      dcn = present(dcndr) .and. present(dcndL)
      grad = present(gradient) .and. present(sigma) .and. dcn
      cpq = present(dqdr) .and. present(dqdL) .and. dcn

      !> Update cache
      allocate (cache)
      call self%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)

      !> Get lattice points
      if (any(mol%periodic)) then
         call get_dir_trans(mol%lattice, trans)
      end if

      !> Get amat
      ndim = mol%nat + 1
      allocate (amat(ndim, ndim))
      call self%get_coulomb_matrix(mol, cache, amat)

      !> Get RHS of ES equation
      allocate (xvec(ndim))
      call self%get_xvec(mol, cache, xvec)

      vrhs = xvec
      ainv = amat

      allocate (ipiv(ndim))
      call sytrf(ainv, ipiv, info=info, uplo='l')

      if (info == 0) then
         if (cpq) then
            call sytri(ainv, ipiv, info=info, uplo='l')
            if (info == 0) then
               call symv(ainv, xvec, vrhs, uplo='l')
               do ic = 1, ndim
                  do jc = ic + 1, ndim
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
         call symv(amat(:mol%nat, :mol%nat), vrhs(:mol%nat), xvec(:mol%nat), alpha=0.5_wp, beta=-1.0_wp, uplo='l')
         energy(:) = energy(:) + vrhs(:mol%nat)*xvec(:mol%nat)
      end if

      !> Allocate and get amat derivatives
      if (grad .or. cpq) then
         allocate (dadr(3, mol%nat, ndim), dadL(3, 3, ndim), atrace(3, mol%nat))
         allocate (dxdr(3, mol%nat, ndim), dxdL(3, 3, ndim))
         call self%get_xvec_derivs(mol, cache, dxdr, dxdL)
         call self%get_coulomb_derivs(mol, cache, vrhs, dadr, dadL, atrace)
         print'(a)', "atrace:"
         print'(3es21.14)', atrace
      end if

      if (grad) then
         gradient = 0.0_wp
         call gemv(dadr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp)!, alpha=0.5_wp)
         call gemv(dxdr(:, :, :mol%nat), vrhs(:mol%nat), gradient, beta=1.0_wp, alpha=-1.0_wp)
         call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
         call gemv(dxdL, vrhs, sigma, beta=1.0_wp, alpha=-1.0_wp)
      end if

      if (cpq) then
         do iat = 1, mol%nat
            dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
            dadr(:, :, iat) = -dxdr(:, :, iat) + dadr(:, :, iat)
            dadL(:, :, iat) = -dxdL(:, :, iat) + dadL(:, :, iat)
         end do

         call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
         call gemm(dadL, ainv(:, :mol%nat), dqdL, alpha=-1.0_wp)
      end if
   end subroutine solve

   subroutine local_charge(self, mol, trans, qloc, dqlocdr, dqlocdL)
      !> Electronegativity equilibration model
      class(mchrg_model_type) :: self
      !> Molecular structure data
      type(structure_type), intent(in) :: mol
      !> Lattice points
      real(wp), intent(in) :: trans(:, :)
      !> Local atomic partial charges
      real(wp), intent(out) :: qloc(:)
      !> Optional derivative of local atomic partial charges w.r.t. atomic positions
      real(wp), intent(out), optional :: dqlocdr(3, mol%nat, mol%nat)
      !> Optional derivative of local atomic partial charges w.r.t. lattice vectors
      real(wp), intent(out), optional :: dqlocdL(3, 3, mol%nat)

      qloc = 0.0_wp
      if (present(dqlocdr) .and. present(dqlocdL)) then
         dqlocdr = 0.0_wp
         dqlocdL = 0.0_wp
      end if
      ! Get the electronegativity weighted CN for local charge
      ! Derivatives depend only in this CN
      if (allocated(self%ncoord_en)) then
         call self%ncoord_en%get_coordination_number(mol, trans, qloc, dqlocdr, dqlocdL)
      end if

      ! Distribute the total charge equally
      qloc = qloc + mol%charge/real(mol%nat, wp)

   end subroutine local_charge

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

end module multicharge_model_type
