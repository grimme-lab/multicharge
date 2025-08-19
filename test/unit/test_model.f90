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

module test_model
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge_model, only: mchrg_model_type
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   use multicharge_model_cache, only: cache_container
   use multicharge_blas, only: gemv
   use multicharge_charge, only: get_charges, get_eeq_charges, get_eeqbc_charges
   implicit none
   private

   public :: collect_model

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

contains

!> Collect all exported unit tests
subroutine collect_model(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("eeq-dadr-mb01", test_eeq_dadr_mb01), &
      & new_unittest("eeq-dadL-mb01", test_eeq_dadL_mb01), &
      & new_unittest("eeq-dbdr-mb01", test_eeq_dbdr_mb01), &
      & new_unittest("eeq-dbdL-mb01", test_eeq_dbdL_mb01), &
      & new_unittest("eeq-charges-mb01", test_eeq_q_mb01), &
      & new_unittest("eeq-charges-mb02", test_eeq_q_mb02), &
      & new_unittest("eeq-charges-actinides", test_eeq_q_actinides), &
      & new_unittest("eeq-energy-mb03", test_eeq_e_mb03), &
      & new_unittest("eeq-energy-mb04", test_eeq_e_mb04), &
      & new_unittest("eeq-gradient-mb05", test_eeq_g_mb05), &
      & new_unittest("eeq-gradient-mb06", test_eeq_g_mb06), &
      & new_unittest("eeq-sigma-mb07", test_eeq_s_mb07), &
      & new_unittest("eeq-sigma-mb08", test_eeq_s_mb08), &
      & new_unittest("eeq-dqdr-mb09", test_eeq_dqdr_mb09), &
      & new_unittest("eeq-dqdr-mb10", test_eeq_dqdr_mb10), &
      & new_unittest("eeq-dqdL-mb11", test_eeq_dqdL_mb11), &
      & new_unittest("eeq-dqdL-mb12", test_eeq_dqdL_mb12), &
      & new_unittest("gradient-h2plus", test_g_h2plus), &
      & new_unittest("eeq-dadr-znooh", test_eeq_dadr_znooh), &
      & new_unittest("eeq-dbdr-znooh", test_eeq_dbdr_znooh), &
      & new_unittest("gradient-znooh", test_g_znooh), &
      & new_unittest("dqdr-znooh", test_dqdr_znooh), &
      & new_unittest("eeqbc-dadr-mb01", test_eeqbc_dadr_mb01), &
      & new_unittest("eeqbc-dadL-mb01", test_eeqbc_dadL_mb01), &
      & new_unittest("eeqbc-dbdr-mb01", test_eeqbc_dbdr_mb01), &
      & new_unittest("eeqbc-dbdL-mb01", test_eeqbc_dbdL_mb01), &
      & new_unittest("eeqbc-dadr-mb05", test_eeqbc_dadr_mb05), &
      & new_unittest("eeqbc-dadL-mb05", test_eeqbc_dadL_mb05), &
      & new_unittest("eeqbc-dbdr-mb05", test_eeqbc_dbdr_mb05), &
      & new_unittest("eeqbc-charges-mb01", test_eeqbc_q_mb01), &
      & new_unittest("eeqbc-charges-mb02", test_eeqbc_q_mb02), &
      & new_unittest("eeqbc-charges-actinides", test_eeqbc_q_actinides), &
      & new_unittest("eeqbc-energy-mb03", test_eeqbc_e_mb03), &
      & new_unittest("eeqbc-energy-mb04", test_eeqbc_e_mb04), &
      & new_unittest("eeqbc-gradient-mb05", test_eeqbc_g_mb05), &
      & new_unittest("eeqbc-gradient-mb06", test_eeqbc_g_mb06), &
      & new_unittest("eeqbc-sigma-mb07", test_eeqbc_s_mb07), &
      & new_unittest("eeqbc-sigma-mb08", test_eeqbc_s_mb08), &
      & new_unittest("eeqbc-dqdr-mb09", test_eeqbc_dqdr_mb09), &
      & new_unittest("eeqbc-dqdr-mb10", test_eeqbc_dqdr_mb10), &
      & new_unittest("eeqbc-dqdL-mb11", test_eeqbc_dqdL_mb11), &
      & new_unittest("eeqbc-dqdL-mb12", test_eeqbc_dqdL_mb12) &
      & ]

end subroutine collect_model

subroutine test_dadr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic, jat, kat
   real(wp) :: thr2_local
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: qloc(:)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: qvec(:), numgrad(:, :, :), amatr1(:, :), amatr2(:, :), amatl1(:, :), amatl2(:, :), numtrace(:, :)
   type(cache_container), allocatable :: cache
   allocate (cache)

   allocate (cn(mol%nat), qloc(mol%nat), amatr1(mol%nat + 1, mol%nat + 1), amatl1(mol%nat + 1, mol%nat + 1), &
      & amatr2(mol%nat + 1, mol%nat + 1), amatl2(mol%nat + 1, mol%nat + 1), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), dqlocdr(3, mol%nat, mol%nat), &
      & dqlocdL(3, 3, mol%nat), dadr(3, mol%nat, mol%nat + 1), dadL(3, 3, mol%nat + 1), &
      & atrace(3, mol%nat), numtrace(3, mol%nat), numgrad(3, mol%nat, mol%nat + 1), qvec(mol%nat))

   ! Set tolerance higher if testing eeqbc model
   select type (model)
   type is (eeqbc_model)
      thr2_local = 3.0_wp*thr2
   class default
      thr2_local = thr2
   end select

   ! Obtain the vector of charges
   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, qvec=qvec)
   if (allocated(error)) return

   numgrad = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         amatr1(:, :) = 0.0_wp
         amatr2(:, :) = 0.0_wp
         amatl1(:, :) = 0.0_wp
         amatl2(:, :) = 0.0_wp

         ! First right-hand side (x+h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatr1)

         ! Second right-hand side (x+2h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatr2)

         ! Return to original position before calculating left sides
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step

         ! First left-hand side (x-h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl1)

         ! Second left-hand side (x-2h)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl2)

         ! Return to original position
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + 2*step

         do kat = 1, mol%nat
            do jat = 1, mol%nat
               ! Numerical gradient using 4-step central difference formula
               ! f'(x) ≈ [f(x-2h) - 8f(x-h) + 8f(x+h) - f(x+2h)]/(12h)
               numgrad(ic, iat, kat) = numgrad(ic, iat, kat) + &
                  & qvec(jat)*(amatl2(kat, jat) - 8.0_wp*amatl1(kat, jat) + &
                  & 8.0_wp*amatr1(kat, jat) - amatr2(kat, jat))/(12.0_wp*step)
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

   if (any(abs(dadr(:, :, :) - numgrad(:, :, :)) > thr2_local)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print'(a)', "dadr:"
      print'(3es21.12)', dadr
      print'(a)', "numgrad:"
      print'(3es21.12)', numgrad
      print'(a)', "diff:"
      print'(3es21.12)', dadr - numgrad
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)
   real(wp), allocatable :: xyz(:, :)
   real(wp), allocatable :: qvec(:), numsigma(:, :, :), amatr(:, :), amatl(:, :)
   real(wp) :: eps(3, 3)
   type(cache_container), allocatable :: cache
   allocate (cache)

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & amatr(mol%nat + 1, mol%nat + 1), amatl(mol%nat + 1, mol%nat + 1), &
      & dadr(3, mol%nat, mol%nat + 1), dadL(3, 3, mol%nat + 1), atrace(3, mol%nat), &
      & numsigma(3, 3, mol%nat + 1), qvec(mol%nat), xyz(3, mol%nat))

   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, qvec=qvec)
   if (allocated(error)) return

   numsigma = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lp: do ic = 1, 3
      do jc = 1, 3
         amatr(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatr)
         if (allocated(error)) exit lp

         amatl(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_coulomb_matrix(mol, cache, amatl)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         do iat = 1, mol%nat
            ! Numerical sigma of the a matrix
            numsigma(jc, ic, :) = numsigma(jc, ic, :) + &
               & 0.5_wp*qvec(iat)*(amatr(iat, :) - amatl(iat, :))/step
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
      print'(3es21.12)', dadL
      print'(a)', "numsigma:"
      print'(3es21.12)', numsigma
      print'(a)', "diff:"
      print'(3es21.12)', dadL - numsigma
   end if

end subroutine test_dadL

subroutine test_dbdr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dbdr(:, :, :), dbdL(:, :, :)
   real(wp), allocatable :: numgrad(:, :, :), xvecr(:), xvecl(:)
   type(cache_container), allocatable :: cache
   allocate (cache)

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & xvecr(mol%nat + 1), xvecl(mol%nat + 1), numgrad(3, mol%nat, mol%nat + 1), &
      & dbdr(3, mol%nat, mol%nat + 1), dbdL(3, 3, mol%nat + 1))

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
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat, :) = 0.5_wp*(xvecr(:) - xvecl(:))/step
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: dbdr(:, :, :), dbdL(:, :, :)
   real(wp), allocatable :: numsigma(:, :, :), xvecr(:), xvecl(:)
   real(wp), allocatable :: xyz(:, :)
   real(wp) :: eps(3, 3)
   type(cache_container), allocatable :: cache
   allocate (cache)

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & xvecr(mol%nat + 1), xvecl(mol%nat + 1), numsigma(3, 3, mol%nat + 1), &
      & dbdr(3, mol%nat, mol%nat + 1), dbdL(3, 3, mol%nat + 1), xyz(3, mol%nat))

   numsigma = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lp: do ic = 1, 3
      do jc = 1, 3
         ! Right-hand side
         xvecr(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, cache, cn, qloc)
         call model%get_xvec(mol, cache, xvecl)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         do iat = 1, mol%nat
            numsigma(jc, ic, iat) = 0.5_wp*(xvecr(iat) - xvecl(iat))/step
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

   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: cn(:), qloc(:)
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)

   allocate (cn(mol%nat), qloc(mol%nat))

   call model%ncoord%get_coordination_number(mol, trans, cn)
   if (allocated(model%ncoord_en)) then
      call model%local_charge(mol, trans, qloc)
   end if

   if (present(eref)) then
      allocate (energy(mol%nat))
      energy(:) = 0.0_wp
   end if
   if (present(qref)) then
      allocate (qvec(mol%nat))
   end if

   call model%solve(mol, error, cn, qloc, energy=energy, qvec=qvec)
   if (allocated(error)) return

   if (present(qref)) then
      if (any(abs(qvec - qref) > thr)) then
         call test_failed(error, "Partial charges do not match")
         print'(a)', "Charges:"
         print'(3es21.14)', qvec
         print'(a)', "diff:"
         print'(3es21.14)', qvec - qref
      end if
   end if
   if (allocated(error)) return

   if (present(eref)) then
      if (any(abs(energy - eref) > thr)) then
         call test_failed(error, "Energies do not match")
         print'(a)', "Energy:"
         print'(3es21.14)', energy
         print'(a)', "diff:"
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :)
   real(wp) :: er, el

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), sigma(3, 3), numgrad(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         energy(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, gradient=gradient, sigma=sigma)
   if (allocated(error)) return

   if (any(abs(gradient(:, :) - numgrad(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
      print'(a)', "Energy gradient:"
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :)
   real(wp) :: er, el, eps(3, 3), numsigma(3, 3), sigma(3, 3)

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), xyz(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         lattr(:, :) = trans
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
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
      print'(a)', "Energy strain:"
      print'(3es21.14)', sigma
      print'(a)', "numsigma:"
      print'(3es21.14)', numsigma
      print'(a)', "diff:"
      print'(3es21.14)', sigma - numsigma
   end if

end subroutine test_numsigma

subroutine test_numdqdr(error, mol, model)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: ql(:), qr(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & ql(mol%nat), qr(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=qr)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=ql)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp*(qr - ql)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)

   call model%solve(mol, error, cn, qloc, dcndr, dcndL, &
      & dqlocdr, dqlocdL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdr(:, :, :) - numdr(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of charges does not match")
      print'(a)', "Charge gradient:"
      print'(3es21.14)', dqdr
      print'(a)', "numgrad:"
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: qr(:), ql(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :), numdL(:, :, :)
   real(wp) :: eps(3, 3)

   allocate (cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qloc(mol%nat), dqlocdr(3, mol%nat, mol%nat), dqlocdL(3, 3, mol%nat), &
      & qr(mol%nat), ql(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & xyz(3, mol%nat), numdL(3, 3, mol%nat))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=qr)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%solve(mol, error, cn, qloc, qvec=ql)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp*(qr - ql)/step
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
      print'(a)', "Charge gradient:"
      print'(3es21.14)', dqdL
      print'(a)', "numgrad:"
      print'(3es21.14)', numdL
      print'(a)', "diff:"
      print'(3es21.14)', dqdL - numdL
   end if

end subroutine test_numdqdL

subroutine test_eeq_dadr_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeq_dadr_mb01

subroutine test_eeq_dadL_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeq_dadL_mb01

subroutine test_eeq_dbdr_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model)

end subroutine test_eeq_dbdr_mb01

subroutine test_eeq_dbdL_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdL(error, mol, model)

end subroutine test_eeq_dbdL_mb01

subroutine test_eeq_q_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & 7.73347759615437E-1_wp, 1.07626897257271E-1_wp, -3.66999554268267E-1_wp, &
      & 4.92833775451616E-2_wp, -1.83332153007808E-1_wp, 2.33302084420314E-1_wp, &
      & 6.61837602813735E-2_wp, -5.43944147972069E-1_wp, -2.70264297953247E-1_wp, &
      & 2.66618970100409E-1_wp, 2.62725030332215E-1_wp, -7.15315061940473E-2_wp, &
      &-3.73300836681036E-1_wp, 3.84585142200261E-2_wp, -5.05851076468890E-1_wp, &
      & 5.17677178773158E-1_wp]

   real(wp), allocatable :: qvec(:)

   call get_structure(mol, "MB16-43", "01")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)
   if (allocated(error)) return

   ! Check wrapper functions
   allocate (qvec(mol%nat), source=0.0_wp)
   call get_charges(model, mol, error, qvec)
   if (allocated(error)) return

   if (any(abs(qvec - ref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print'(a)', "Charges:"
      print'(3es21.14)', qvec
      print'(a)', "diff:"
      print'(3es21.14)', qvec - ref
   end if
   if (allocated(error)) return

   qvec = 0.0_wp
   call get_eeq_charges(mol, error, qvec)
   if (allocated(error)) return

   if (any(abs(qvec - ref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print'(a)', "Charges:"
      print'(3es21.14)', qvec
      print'(a)', "diff:"
      print'(3es21.14)', qvec - ref
   end if
   if (allocated(error)) return

end subroutine test_eeq_q_mb01

subroutine test_eeq_q_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & 7.38394752482521E-2_wp, -1.68354859084778E-1_wp, -3.47642846358022E-1_wp, &
      &-7.05489251302223E-1_wp, 7.73548241620680E-1_wp, 2.30207580650128E-1_wp, &
      & 1.02748505731185E-1_wp, 9.47818154871089E-2_wp, 2.44259536057649E-2_wp, &
      & 2.34984928231320E-1_wp, -3.17839956573785E-1_wp, 6.67112952465234E-1_wp, &
      &-4.78119957747208E-1_wp, 6.57536208287042E-2_wp, 1.08259091466373E-1_wp, &
      &-3.58215294268738E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeq_q_mb02

subroutine test_eeq_q_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(17) = [&
      & 1.86904766283711E-02_wp, 2.89972818160259E-01_wp, 3.59298070941105E-02_wp, &
      &-4.61256458126589E-02_wp, -7.02605348653647E-02_wp, -7.42052215689073E-02_wp, &
      &-8.21938718945845E-02_wp, 1.64953118841151E-01_wp, 2.10381640633390E-01_wp, &
      &-6.65485355096282E-02_wp, -2.59873890255450E-01_wp, 1.33839147940414E-01_wp, &
      & 7.20768968601809E-02_wp, -3.36652347675997E-03_wp, -1.14546280789657E-01_wp, &
      &-8.55922398441004E-02_wp, -1.23131162140762E-01_wp]

   !> Molecular structure data
   mol%nat = 17
   mol%nid = 17
   mol%id = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, &
      & 12, 13, 14, 15, 16, 17]
   mol%num = [87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      & 98, 99, 100, 101, 102, 103]
   mol%xyz = reshape([ &
      & 0.98692316414074_wp, 6.12727238368797_wp, -6.67861597188102_wp, &
      & 3.63898862390869_wp, 5.12109301182962_wp, 3.01908613326278_wp, &
      & 5.14503571563551_wp, -3.97172984617710_wp, 3.82011791828867_wp, &
      & 6.71986847575494_wp, 1.71382138402812_wp, 3.92749159076307_wp, &
      & 4.13783589704826_wp, -2.10695793491818_wp, 0.19753203068899_wp, &
      & 8.97685097698326_wp, -3.08813636191844_wp, -4.45568615593938_wp, &
      & 12.5486412940776_wp, -1.77128765259458_wp, 0.59261498922861_wp, &
      & 7.82051475868325_wp, -3.97159756604558_wp, -0.53637703616916_wp, &
      &-0.43444574624893_wp, -1.69696511583960_wp, -1.65898182093050_wp, &
      &-4.71270645149099_wp, -0.11534827468942_wp, 2.84863373521297_wp, &
      &-2.52061680335614_wp, 1.82937752749537_wp, -2.10366982879172_wp, &
      & 0.13551154616576_wp, 7.99805359235043_wp, -1.55508522619903_wp, &
      & 3.91594542499717_wp, -1.72975169129597_wp, -5.07944366756113_wp, &
      &-1.03393930231679_wp, 4.69307230054046_wp, 0.02656940927472_wp, &
      & 6.20675384557240_wp, 4.24490721493632_wp, -0.71004195169885_wp, &
      & 7.04586341131562_wp, 5.20053667939076_wp, -7.51972863675876_wp, &
      & 2.01082807362334_wp, 1.34838807211157_wp, -4.70482633508447_wp],&
      & [3, 17])
   mol%periodic = [.false.]

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeq_q_actinides

subroutine test_eeq_e_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      &-2.18698345033562E-1_wp, -1.04793885931268E+0_wp, 4.78963353574572E-2_wp, &
      & 5.76566377591676E-1_wp, 7.37187470977927E-1_wp, 8.06020047053305E-2_wp, &
      &-4.19837955782898E-1_wp, 5.49627510550566E-2_wp, 8.01486728591565E-2_wp, &
      & 1.00618944521776E-1_wp, -6.61715169034150E-1_wp, -3.60531647289563E-1_wp, &
      &-4.87729666337974E-1_wp, 2.48257554279938E-1_wp, 6.96027176590956E-1_wp, &
      & 4.31679925875087E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeq_e_mb03

subroutine test_eeq_e_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & 1.13974214746111E-1_wp, -4.41735365367827E-1_wp, 8.99036489938394E-2_wp, &
      &-2.97539904703271E-1_wp, 8.05174117097006E-3_wp, 1.31105783760276E-1_wp, &
      & 1.54594451996644E-1_wp, 1.19929653841255E-1_wp, 1.26056586309101E-1_wp, &
      & 1.78439005754586E-1_wp, -1.98703462666082E-1_wp, 4.19630120027785E-1_wp, &
      & 7.05569220334930E-2_wp, -4.50925107441869E-1_wp, 1.39289602382354E-1_wp, &
      &-2.67853086061429E-1_wp]

   call get_structure(mol, "MB16-43", "04")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeq_e_mb04

subroutine test_eeq_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_eeq_g_mb05

subroutine test_eeq_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "06")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_eeq_g_mb06

subroutine test_eeq_s_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "07")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeq_s_mb07

subroutine test_eeq_s_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "08")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeq_s_mb08

subroutine test_eeq_dqdr_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "09")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_eeq_dqdr_mb09

subroutine test_eeq_dqdr_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "10")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_eeq_dqdr_mb10

subroutine test_eeq_dqdL_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "11")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine test_eeq_dqdL_mb11

subroutine test_eeq_dqdL_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "12")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine test_eeq_dqdL_mb12

subroutine test_g_h2plus(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   integer, parameter :: nat = 2
   real(wp), parameter :: charge = 1.0_wp
   integer, parameter :: num(nat) = [1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & +0.00000000000000_wp, +0.00000000000000_wp, +0.00000000000000_wp, &
      & +1.00000000000000_wp, +0.00000000000000_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_g_h2plus

subroutine test_eeq_dadr_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
      & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
      & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
      & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeq_dadr_znooh

subroutine test_eeq_dbdr_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
      & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
      & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
      & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model)

end subroutine test_eeq_dbdr_znooh

subroutine test_g_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
      & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
      & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
      & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_g_znooh

subroutine test_dqdr_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   integer, parameter :: nat = 4
   real(wp), parameter :: charge = -1.0_wp
   integer, parameter :: num(nat) = [30, 8, 8, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & -0.30631629283878_wp, -1.11507514203552_wp, +0.00000000000000_wp, &
      & -0.06543072660074_wp, -4.32862093666082_wp, +0.00000000000000_wp, &
      & -0.64012239724097_wp, +2.34966763895920_wp, +0.00000000000000_wp, &
      & +1.01186941668051_wp, +3.09402843973713_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_dqdr_znooh

subroutine test_eeqbc_dadr_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeqbc_dadr_mb01

subroutine test_eeqbc_dadL_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeqbc_dadL_mb01

subroutine test_eeqbc_dbdr_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model)

end subroutine test_eeqbc_dbdr_mb01

subroutine test_eeqbc_dbdL_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdL(error, mol, model)

end subroutine test_eeqbc_dbdL_mb01

subroutine test_eeqbc_dadr_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeqbc_dadr_mb05

subroutine test_eeqbc_dadL_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeqbc_dadL_mb05

subroutine test_eeqbc_dbdr_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dbdr(error, mol, model)

end subroutine test_eeqbc_dbdr_mb05

subroutine test_eeqbc_q_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      & 4.75783090912440E-1_wp, -4.26540500638442E-2_wp, -3.77871226005535E-1_wp, &
      &-9.67376090029522E-2_wp, -1.73364116997142E-1_wp, 1.08660101025683E-1_wp, &
      &-1.13628448410420E-1_wp, -3.17939699645693E-1_wp, -2.45655524697400E-1_wp, &
      & 1.76106572419156E-1_wp, 1.14510850652006E-1_wp, -1.22241025474265E-1_wp, &
      &-1.44595425453640E-2_wp, 2.57782082780412E-1_wp, -1.11777579535162E-1_wp, &
      & 4.83486124588080E-1_wp]

   real(wp), allocatable :: qvec(:)

   call get_structure(mol, "MB16-43", "01")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

   ! Check wrapper functions
   allocate (qvec(mol%nat), source=0.0_wp)
   call get_charges(model, mol, error, qvec)
   if (allocated(error)) return

   if (any(abs(qvec - ref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print'(a)', "Charges:"
      print'(3es21.14)', qvec
      print'(a)', "diff:"
      print'(3es21.14)', qvec - ref
   end if
   if (allocated(error)) return

   qvec = 0.0_wp
   call get_eeqbc_charges(mol, error, qvec)
   if (allocated(error)) return

   if (any(abs(qvec - ref) > thr)) then
      call test_failed(error, "Partial charges do not match")
      print'(a)', "Charges:"
      print'(3es21.14)', qvec
      print'(a)', "diff:"
      print'(3es21.14)', qvec - ref
   end if
   if (allocated(error)) return

end subroutine test_eeqbc_q_mb01

subroutine test_eeqbc_q_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      &-7.89571755894845E-2_wp, -1.84724587297173E-1_wp, -1.63060175795952E-2_wp, &
      &-2.36115890461711E-1_wp, 5.05729582512203E-1_wp, 1.37556939519704E-1_wp, &
      &-2.29048340967271E-2_wp, -4.31722346626804E-2_wp, 2.26466952977883E-1_wp, &
      & 1.25047857913714E-1_wp, 6.72899182661252E-3_wp, 3.08986208662492E-1_wp, &
      &-3.34344661086462E-1_wp, -3.16758668376149E-2_wp, -5.24170403450005E-2_wp, &
      &-3.09898225456160E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeqbc_q_mb02

subroutine test_eeqbc_q_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(17) = [&
      & 9.27195802124755E-2_wp, -2.78358027117801E-1_wp, 1.71815557281178E-1_wp, &
      & 7.85579953672371E-2_wp, -1.08186262417305E-2_wp, -4.81860290986309E-2_wp, &
      & 1.57794666483371E-1_wp, -1.61830258916072E-1_wp, -2.76569765724910E-1_wp, &
      & 2.99654899926371E-1_wp, -5.24433579322476E-1_wp, -1.99523360511699E-1_wp, &
      &-3.42285450387671E-2_wp, -3.15076271542101E-2_wp, 1.49700940990172E-1_wp, &
      & 1.45447393911445E-1_wp, 4.69764784954047E-1_wp]

   ! Molecular structure data
   mol%nat = 17
   mol%nid = 17
   mol%id = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, &
      & 12, 13, 14, 15, 16, 17]
   mol%num = [87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      & 98, 99, 100, 101, 102, 103]
   mol%xyz = reshape([ &
      & 0.98692316414074_wp, 6.12727238368797_wp, -6.67861597188102_wp, &
      & 3.63898862390869_wp, 5.12109301182962_wp, 3.01908613326278_wp, &
      & 5.14503571563551_wp, -3.97172984617710_wp, 3.82011791828867_wp, &
      & 6.71986847575494_wp, 1.71382138402812_wp, 3.92749159076307_wp, &
      & 4.13783589704826_wp, -2.10695793491818_wp, 0.19753203068899_wp, &
      & 8.97685097698326_wp, -3.08813636191844_wp, -4.45568615593938_wp, &
      & 12.5486412940776_wp, -1.77128765259458_wp, 0.59261498922861_wp, &
      & 7.82051475868325_wp, -3.97159756604558_wp, -0.53637703616916_wp, &
      &-0.43444574624893_wp, -1.69696511583960_wp, -1.65898182093050_wp, &
      &-4.71270645149099_wp, -0.11534827468942_wp, 2.84863373521297_wp, &
      &-2.52061680335614_wp, 1.82937752749537_wp, -2.10366982879172_wp, &
      & 0.13551154616576_wp, 7.99805359235043_wp, -1.55508522619903_wp, &
      & 3.91594542499717_wp, -1.72975169129597_wp, -5.07944366756113_wp, &
      &-1.03393930231679_wp, 4.69307230054046_wp, 0.02656940927472_wp, &
      & 6.20675384557240_wp, 4.24490721493632_wp, -0.71004195169885_wp, &
      & 7.04586341131562_wp, 5.20053667939076_wp, -7.51972863675876_wp, &
      & 2.01082807362334_wp, 1.34838807211157_wp, -4.70482633508447_wp],&
      & [3, 17])
   mol%periodic = [.false.]

   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, qref=ref)

end subroutine test_eeqbc_q_actinides

subroutine test_eeqbc_e_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      &-6.96992195046228E-2_wp, -1.62155815983893E+0_wp, -1.38060751929644E-3_wp, &
      &-9.06342279911342E-1_wp, -1.83281566961757E+0_wp, -1.20333262207652E-1_wp, &
      &-6.51187555181622E-1_wp, -3.27410111288548E-3_wp, -8.00565881078213E-3_wp, &
      &-2.60385867643294E-2_wp, -9.33285940415006E-1_wp, -1.48859947660327E-1_wp, &
      &-7.19456827995756E-1_wp, -9.58311834831915E-2_wp, -1.54672086637309E+0_wp, &
      &-1.03483694342593E-5_wp]

   call get_structure(mol, "MB16-43", "03")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeqbc_e_mb03

subroutine test_eeqbc_e_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(16) = [&
      &-3.91054587109712E-2_wp, -8.21933095021462E-4_wp, -1.28550631772418E-2_wp, &
      &-8.95571658260288E-2_wp, -4.94655224590082E-1_wp, -3.34598696522549E-2_wp, &
      &-3.75768676247744E-2_wp, -1.36087478076862E-2_wp, -2.07985587717960E-3_wp, &
      &-1.17711924662077E-2_wp, -2.68707428024071E-1_wp, -1.00650791933494E+0_wp, &
      &-5.64487253409848E-2_wp, -4.89693252471477E-1_wp, -3.74734977139679E-2_wp, &
      &-9.22642011641358E-3_wp]

   call get_structure(mol, "MB16-43", "04")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call gen_test(error, mol, model, eref=ref)

end subroutine test_eeqbc_e_mb04

subroutine test_eeqbc_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "05")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_eeqbc_g_mb05

subroutine test_eeqbc_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "06")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)

end subroutine test_eeqbc_g_mb06

subroutine test_eeqbc_s_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "07")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeqbc_s_mb07

subroutine test_eeqbc_s_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "08")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeqbc_s_mb08

subroutine test_eeqbc_dqdr_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "09")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_eeqbc_dqdr_mb09

subroutine test_eeqbc_dqdr_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "10")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)

end subroutine test_eeqbc_dqdr_mb10

subroutine test_eeqbc_dqdL_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "11")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine test_eeqbc_dqdL_mb11

subroutine test_eeqbc_dqdL_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "MB16-43", "12")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine test_eeqbc_dqdL_mb12

end module test_model
