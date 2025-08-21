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

module test_pbc
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only: structure_type
   use mctc_cutoff, only: get_lattice_points
   use mstore, only: get_structure
   use multicharge_model, only: mchrg_model_type
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   use multicharge_model_cache, only: cache_container
   implicit none
   private

   public :: collect_pbc

   real(wp), parameter :: thr = 1000 * epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_pbc(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("eeq-charges-cyanamide", test_eeq_q_cyanamide), &
      & new_unittest("eeq-energy-formamide", test_eeq_e_formamide), &
      & new_unittest("eeq-dbdr-co2", test_eeq_dbdr_co2), &
      & new_unittest("eeq-dbdL-co2", test_eeq_dbdL_co2), &
      & new_unittest("eeq-dadr-ice", test_eeq_dadr_ice), &
      & new_unittest("eeq-dadL-ice", test_eeq_dadL_ice), &
      & new_unittest("eeq-gradient-co2", test_eeq_g_co2), &
      & new_unittest("eeq-sigma-ice", test_eeq_s_ice), &
      & new_unittest("eeq-dqdr-urea", test_eeq_dqdr_urea), &
      & new_unittest("eeq-dqdL-oxacb", test_eeq_dqdL_oxacb), &
      & new_unittest("eeqbc-dbdr-co2", test_eeqbc_dbdr_co2), &
      & new_unittest("eeqbc-dbdL-co2", test_eeqbc_dbdL_co2), &
      & new_unittest("eeqbc-dadr-ice", test_eeqbc_dadr_ice), &
      & new_unittest("eeqbc-dadL-ice", test_eeqbc_dadL_ice), &
      & new_unittest("eeqbc-gradient-co2", test_eeqbc_g_co2), &
      & new_unittest("eeqbc-sigma-ice", test_eeqbc_s_ice), &
      & new_unittest("eeqbc-dqdr-urea", test_eeqbc_dqdr_urea), &
      & new_unittest("eeqbc-dqdL-oxacb", test_eeqbc_dqdL_oxacb) &
      & ]

end subroutine collect_pbc

subroutine gen_test(error, mol, model, qref, eref)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Reference charges
   real(wp), intent(in), optional :: qref(:)

   !> Reference energies
   real(wp), intent(in), optional :: eref(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), qloc(:), trans(:, :)
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   allocate(cn(mol%nat), qloc(mol%nat))

   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)

   if (present(eref)) then
      allocate(energy(mol%nat))
      energy(:) = 0.0_wp
   end if
   if (present(qref)) then
      allocate(qvec(mol%nat))
   end if

   call model%solve(mol, error, cn, qloc, energy=energy, qvec=qvec)
   if (allocated(error)) return

   if (present(qref)) then
      if (any(abs(qvec - qref) > thr)) then
         call test_failed(error, "Partial charges do not match")
         print'(a)', "Charges:"
         print'(3es21.14)', qvec
         print'("---")'
         print'(3es21.14)', qref
         print'("---")'
         print'(3es21.14)', qvec - qref
      end if
   end if
   if (allocated(error)) return

   if (present(eref)) then
      if (any(abs(energy - eref) > thr)) then
         call test_failed(error, "Energies do not match")
         print'(a)', "Energy:"
         print'(3es21.14)', energy
         print'("---")'
         print'(3es21.14)', eref
         print'("---")'
         print'(3es21.14)', energy - eref
      end if
   end if

end subroutine gen_test

subroutine test_numgrad(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :)
   real(wp) :: er, el

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), sigma(3, 3), numgrad(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         energy(:) = 0.0_wp
         er = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp * (er - el) / step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   energy(:) = 0.0_wp
   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, energy, gradient, sigma)
   if (allocated(error)) return

   if (any(abs(gradient(:, :) - numgrad(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print'(a)', "gradient:"
      print'(3es21.14)', gradient
      print'(a)', "numgrad:"
      print'(3es21.14)', numgrad
      print'(a)', "diff:"
      print'(3es21.14)', gradient - numgrad
   end if

end subroutine test_numgrad

subroutine test_numsigma(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, jc
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :)
   real(wp) :: er, el, eps(3, 3), numsigma(3, 3), sigma(3, 3), lattice(3, 3)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), xyz(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         numsigma(jc, ic) = 0.5_wp * (er - el) / step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   energy(:) = 0.0_wp
   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, energy, gradient, sigma)
   if (allocated(error)) return

   if (any(abs(sigma(:, :) - numsigma(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print'(a)', "sigma:"
      print'(3es21.14)', sigma
      print'(a)', "numgrad:"
      print'(3es21.14)', numsigma
      print'(a)', "diff:"
      print'(3es21.14)', sigma(:, :) - numsigma(:, :)
   end if

end subroutine test_numsigma

subroutine test_dbdr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dbdr(:, :, :), dbdL(:, :, :)
   real(wp), allocatable :: numgrad(:, :, :), xvecr(:), xvecl(:)
   type(cache_container), allocatable :: cache
   allocate(cache)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & xvecr(mol%nat + 1), xvecl(mol%nat + 1), numgrad(3, mol%nat, mol%nat + 1), &
      & dbdr(3, mol%nat, mol%nat + 1), dbdL(3, 3, mol%nat + 1))

   numgrad = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         ! Right-hand side
         xvecr(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat, :) = 0.5_wp * (xvecr(:) - xvecl(:)) / step
      end do
   end do lp

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   call model%get_xvec(mol, cache, xvecl) ! need to call this for xtmp in cache (eeqbc)
   call model%get_xvec_derivs(mol, cache, dbdr, dbdL)

   if (any(abs(dbdr(:, :, :) - numgrad(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the b vector does not match")
      print'(a)', "dbdr:"
      print'(3es21.14)', dbdr
      print'(a)', "numgrad:"
      print'(3es21.14)', numgrad
      print'(a)', "diff:"
      print'(3es21.14)', dbdr - numgrad
   end if

end subroutine test_dbdr

subroutine test_dbdL(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic, jc
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dbdr(:, :, :), dbdL(:, :, :)
   real(wp), allocatable :: numsigma(:, :, :), xvecr(:), xvecl(:)
   real(wp), allocatable :: xyz(:, :), lattr(:, :), trans(:, :)
   real(wp) :: lattice(3, 3)
   real(wp) :: eps(3, 3)
   type(cache_container), allocatable :: cache
   allocate(cache)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & xvecr(mol%nat + 1), xvecl(mol%nat + 1), numsigma(3, 3, mol%nat + 1), &
      & dbdr(3, mol%nat, mol%nat + 1), dbdL(3, 3, mol%nat + 1), xyz(3, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   numsigma = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         ! Right-hand side
         xvecr(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecl)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         do iat = 1, mol%nat
            numsigma(jc, ic, iat) = 0.5_wp * (xvecr(iat) - xvecl(iat)) / step
         end do
      end do
   end do lp

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   call model%get_xvec(mol, cache, xvecl) ! need to call this for xtmp in cache (eeqbc)
   call model%get_xvec_derivs(mol, cache, dbdr, dbdL)

   if (any(abs(dbdL(:, :, :) - numsigma(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the b vector does not match")
      print'(a)', "dbdL:"
      print'(3es21.14)', dbdL
      print'(a)', "numsigma:"
      print'(3es21.14)', numsigma
      print'(a)', "diff:"
      print'(3es21.14)', dbdL - numsigma
   end if

end subroutine test_dbdL

subroutine test_dadr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic, jat, kat
   real(wp) :: thr2_local
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: qloc(:)
   real(wp), allocatable :: trans(:, :)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: qvec(:), numgrad(:, :, :), amatr(:, :), amatl(:, :), numtrace(:, :)
   type(cache_container), allocatable :: cache
   allocate(cache)

   allocate(cn(mol%nat), qloc(mol%nat), amatr(mol%nat + 1, mol%nat + 1), amatl(mol%nat + 1, mol%nat + 1), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), dqlocdr(3, mol%nat, mol%nat), &
      & dqlocdL(3, 3, mol%nat), dadr(3, mol%nat, mol%nat + 1), dadL(3, 3, mol%nat + 1), &
      & atrace(3, mol%nat), numtrace(3, mol%nat), numgrad(3, mol%nat, mol%nat + 1), qvec(mol%nat))

   ! Set tolerance higher if testing eeqbc model
   select type(model)
   type is(eeqbc_model)
      thr2_local = 3.0_wp * thr2
   class default
      thr2_local = thr2
   end select

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   ! Obtain the vector of charges
   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, qvec=qvec)
   if (allocated(error)) return

   numgrad = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         ! Right-hand side
         amatr(:, :) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatr)

         ! Left-hand side
         amatl(:, :) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step

         do kat = 1, mol%nat
            do jat = 1, mol%nat
               ! Numerical gradient of the A matrix
               numgrad(ic, iat, kat) = 0.5_wp * qvec(jat) * (amatr(kat, jat) - amatl(kat, jat)) / step &
                  & + numgrad(ic, iat, kat)
            end do
         end do
      end do
   end do lp

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)

   ! Add trace of the A matrix
   do iat = 1, mol%nat
      dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
   end do

   ! higher tolerance for numerical gradient
   if (any(abs(dadr(:, :, :) - numgrad(:, :, :)) > thr2_local)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print'(a)', "dadr:"
      print'(3es21.14)', dadr
      print'(a)', "numgrad:"
      print'(3es21.14)', numgrad
      print'(a)', "diff:"
      print'(3es21.14)', dadr - numgrad
   end if

end subroutine test_dadr

subroutine test_dadL(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, jc, iat
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: xyz(:, :), lattr(:, :), trans(:, :)
   real(wp), allocatable :: qvec(:), numsigma(:, :, :), amatr(:, :), amatl(:, :)
   real(wp) :: lattice(3, 3)
   real(wp) :: eps(3, 3)
   type(cache_container), allocatable :: cache
   allocate(cache)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & amatr(mol%nat + 1, mol%nat + 1), amatl(mol%nat + 1, mol%nat + 1), &
      & dadr(3, mol%nat, mol%nat + 1), dadL(3, 3, mol%nat + 1), atrace(3, mol%nat), &
      & numsigma(3, 3, mol%nat + 1), qvec(mol%nat), xyz(3, mol%nat))

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, qvec=qvec)
   if (allocated(error)) return

   numsigma = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         amatr(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatr)
         if (allocated(error)) exit lp

         amatl(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         do iat = 1, mol%nat
            ! Numerical sigma of the a matrix
            numsigma(jc, ic, :) = 0.5_wp * qvec(iat) * (amatr(iat, :) - amatl(iat, :)) / step + numsigma(jc, ic, :)
         end do
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)

   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)
   if (allocated(error)) return

   if (any(abs(dadL(:, :, :) - numsigma(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print'(a)', "dadL:"
      print'(3es21.14)', dadL
      print'(a)', "numsigma:"
      print'(3es21.14)', numsigma
      print'(a)', "diff:"
      print'(3es21.14)', dadL - numsigma
   end if

end subroutine test_dadL

subroutine test_numdqdr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: ql(:), qr(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & ql(mol%nat), qr(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         qr = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=qr)
         if (allocated(error)) exit lp

         ql = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=ql)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp * (qr - ql) / step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdr(:, :, :) - numdr(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of charges does not match")
      print'(a)', "dqdr:"
      print'(3es21.14)', dqdr
      print'(a)', "numdr:"
      print'(3es21.14)', numdr
      print'(a)', "diff:"
      print'(3es21.14)', dqdr - numdr
   end if

end subroutine test_numdqdr

subroutine test_numdqdL(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, jc
   real(wp), parameter :: cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: qr(:), ql(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :), numdL(:, :, :)
   real(wp) :: eps(3, 3), lattice(3, 3)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & qr(mol%nat), ql(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & xyz(3, mol%nat), numdL(3, 3, mol%nat))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         qr = 0.0_wp
         ql = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%solve(mol, error, cn, qloc, qvec=qr)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%solve(mol, error, cn, qloc, qvec=ql)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp * (qr - ql) / step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdL(:, :, :) - numdL(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of charges does not match")
      print'(a)', "dqdL:"
      print'(3es21.14)', dqdL
      print'(a)', "numdL:"
      print'(3es21.14)', numdL
      print'(a)', "diff:"
      print'(3es21.14)', dqdL - numdL
   end if

end subroutine test_numdqdL

subroutine test_eeq_q_cyanamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(40) = [&
      & 3.47274007973765E-1_wp, 3.47273017259661E-1_wp, 3.47276806231462E-1_wp, &
      & 3.47267235121650E-1_wp, 3.47288751410308E-1_wp, 3.47296840650318E-1_wp, &
      & 3.47299868243747E-1_wp, 3.47299359754750E-1_wp, 3.45885710186856E-1_wp, &
      & 3.45874246032015E-1_wp, 3.45888242047875E-1_wp, 3.45877451600398E-1_wp, &
      & 3.45902365123333E-1_wp, 3.45902162041418E-1_wp, 3.45900336974539E-1_wp, &
      & 3.45903770121699E-1_wp, 3.58996736955583E-1_wp, 3.58987991603151E-1_wp, &
      & 3.59000859720791E-1_wp, 3.58990960169615E-1_wp, 3.58990422147650E-1_wp, &
      & 3.58997842028106E-1_wp, 3.58997124263551E-1_wp, 3.59001728635141E-1_wp, &
      &-5.86637659146611E-1_wp, -5.86601835343945E-1_wp, -5.86625078062709E-1_wp, &
      &-5.86579999312624E-1_wp, -5.86654394920376E-1_wp, -5.86669360667441E-1_wp, &
      &-5.86667481028248E-1_wp, -5.86674873278027E-1_wp, -4.65529092171338E-1_wp, &
      &-4.65527713392387E-1_wp, -4.65546450764424E-1_wp, -4.65540905286785E-1_wp, &
      &-4.65529860758785E-1_wp, -4.65533757413217E-1_wp, -4.65527683275370E-1_wp, &
      &-4.65527691475100E-1_wp]

   call get_structure(mol, "X23", "cyanamide")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeq_q_cyanamide

subroutine test_eeq_e_formamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(24) = [&
      & 4.01878116759118E-1_wp, 4.01884779867146E-1_wp, 4.01836430415694E-1_wp, &
      & 4.01859412981015E-1_wp, 4.30132679242368E-1_wp, 4.30105708698129E-1_wp, &
      & 4.30085668404029E-1_wp, 4.30108945050621E-1_wp, 1.90400193163051E-1_wp, &
      & 1.90401822781550E-1_wp, 1.90404130981992E-1_wp, 1.90407447306916E-1_wp, &
      & 3.33868335196572E-1_wp, 3.33848682655264E-1_wp, 3.33875977331594E-1_wp, &
      & 3.33893030511429E-1_wp, -8.35635584733599E-1_wp, -8.35614926719694E-1_wp, &
      &-8.35542684637248E-1_wp, -8.35589488474056E-1_wp, -6.31079122091240E-1_wp, &
      &-6.31059677948463E-1_wp, -6.31085206912995E-1_wp, -6.31081747027041E-1_wp]

   call get_structure(mol, "X23", "formamide")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeq_e_formamide

subroutine test_eeq_dbdr_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model)

end subroutine test_eeq_dbdr_co2

subroutine test_eeq_dbdL_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdL(error, mol, model)

end subroutine test_eeq_dbdL_co2

subroutine test_eeq_dadr_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeq_dadr_ice

subroutine test_eeq_dadL_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeq_dadL_ice

subroutine test_eeq_g_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_eeq_g_co2

subroutine test_eeq_s_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeq_s_ice

subroutine test_eeq_dqdr_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "urea")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_eeq_dqdr_urea

subroutine test_eeq_dqdL_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "oxacb")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine test_eeq_dqdL_oxacb

subroutine test_eeqbc_dbdr_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model)

end subroutine test_eeqbc_dbdr_co2

subroutine test_eeqbc_dbdL_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdL(error, mol, model)

end subroutine test_eeqbc_dbdL_co2

subroutine test_eeqbc_dadr_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeqbc_dadr_ice

subroutine test_eeqbc_dadL_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeqbc_dadL_ice

subroutine test_eeqbc_g_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_eeqbc_g_co2

subroutine test_eeqbc_s_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "ICE10", "vi")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeqbc_s_ice

subroutine test_eeqbc_dqdr_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "urea")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_eeqbc_dqdr_urea

subroutine test_eeqbc_dqdL_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "oxacb")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine test_eeqbc_dqdL_oxacb

end module test_pbc
