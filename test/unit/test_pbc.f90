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
   use mctc_cutoff, only: get_lattice_points
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_convert, only: aatoau
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge, only: get_charges, add_charge_hessian
   use multicharge_ewald, only: get_alpha, get_ewald_trans
   use multicharge_ewald_cutoff, only: get_ewald_dir_trans, get_ewald_rec_trans
   use multicharge_ewald_special, only: fast_f012_1d, fast_k012_1d, erfc_exp_2d
   use multicharge_model, only: mchrg_model_type
   use multicharge_model_cache, only: cache_container
   use multicharge_model_eeq, only: eeq_cache
   use multicharge_model_eeqbc, only: eeqbc_model
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private

   public :: collect_pbc

   real(wp), parameter :: thr = 1000 * epsilon(1.0_wp)
   real(wp), parameter :: thr2 = 10 * sqrt(epsilon(1.0_wp))

   integer, parameter :: znoh_num(3) = [30, 8, 1]
   real(wp), parameter :: znoh_xyz(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp, 0.1_wp, 0.2_wp, &
      & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3])
   real(wp), parameter :: znoh_lat(3, 3) = reshape([ &
      & 5.1_wp, 0.0_wp, 0.0_wp, &
      & 0.4_wp, 6.2_wp, 0.0_wp, &
      & 0.2_wp, 0.7_wp, 7.1_wp], [3, 3])

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
      & new_unittest("eeq-dadr-znoh", test_eeq_dadr_znoh), &
      & new_unittest("eeq-dadL-znoh", test_eeq_dadL_znoh), &
      & new_unittest("eeq-gradient-co2", test_eeq_g_co2), &
      & new_unittest("eeq-sigma-co2", test_eeq_s_co2), &
      & new_unittest("eeq-dqdr-urea", test_eeq_dqdr_urea), &
      & new_unittest("eeq-dqdL-oxacb", test_eeq_dqdL_oxacb), &
      & new_unittest("eeqbc-dbdr-co2", test_eeqbc_dbdr_co2), &
      & new_unittest("eeqbc-dbdL-co2", test_eeqbc_dbdL_co2), &
      & new_unittest("eeqbc-dadr-znoh", test_eeqbc_dadr_znoh), &
      & new_unittest("eeqbc-response-znoh", test_eeqbc_response_znoh), &
      & new_unittest("eeqbc-dadL-znoh", test_eeqbc_dadL_znoh), &
      & new_unittest("eeqbc-gradient-co2", test_eeqbc_g_co2), &
      & new_unittest("eeqbc-sigma-co2", test_eeqbc_s_co2), &
      & new_unittest("eeqbc-dqdr-urea", test_eeqbc_dqdr_urea), &
      & new_unittest("eeqbc-dqdL-oxacb", test_eeqbc_dqdL_oxacb), &
      & new_unittest("eeq-1d", test_eeq_1d), &
      & new_unittest("eeq-2d", test_eeq_2d), &
      & new_unittest("eeqbc-1d", test_eeqbc_1d), &
      & new_unittest("eeqbc-2d", test_eeqbc_2d), &
      & new_unittest("ewald-setup", test_ewald_setup), &
      & new_unittest("ewald-alpha-1d", test_ewald_alpha_1d), &
      & new_unittest("ewald-alpha-2d", test_ewald_alpha_2d), &
      & new_unittest("ewald-lattice-invariance", test_ewald_lattice_invariance), &
      & new_unittest("eeq-slab-285", test_eeq_slab_285), &
      & new_unittest("eeq-slab-3d-limit", test_eeq_slab_3d_limit), &
      & new_unittest("supercell-1d", test_supercell_1d), &
      & new_unittest("supercell-2d", test_supercell_2d), &
      & new_unittest("supercell-3d", test_supercell_3d), &
      & new_unittest("special-functions", test_special_functions) &
      & ]

end subroutine collect_pbc

!> Compact skew cell probes nonlocal widths and their periodic self images.
subroutine test_eeqbc_response_znoh(error)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, [30, 8, 1], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
      & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), charge=1.0_wp, &
      & lattice=reshape([5.1_wp, 0.0_wp, 0.0_wp, 0.4_wp, 6.2_wp, 0.0_wp, &
      & 0.2_wp, 0.7_wp, 7.1_wp], [3, 3]), periodic=[.true., .true., .true.])
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)
   if (allocated(error)) return
   call test_dadL(error, mol, model)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)
end subroutine test_eeqbc_response_znoh

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

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
         print"(a)", "Charges:"
         print"(3es21.14)", qvec
         print'("---")'
         print"(3es21.14)", qref
         print'("---")'
         print"(3es21.14)", qvec - qref
      end if
   end if
   if (allocated(error)) return

   if (present(eref)) then
      if (any(abs(energy - eref) > thr)) then
         call test_failed(error, "Energies do not match")
         print"(a)", "Energy:"
         print"(3es21.14)", energy
         print'("---")'
         print"(3es21.14)", eref
         print'("---")'
         print"(3es21.14)", energy - eref
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

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
      print"(a)", "gradient:"
      print"(3es21.14)", gradient
      print"(a)", "numgrad:"
      print"(3es21.14)", numgrad
      print"(a)", "diff:"
      print"(3es21.14)", gradient - numgrad
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
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :)
   real(wp) :: er, el, eps(3, 3), numsigma(3, 3), sigma(3, 3), lattice(3, 3)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
      print"(a)", "sigma:"
      print"(3es21.14)", sigma
      print"(a)", "numgrad:"
      print"(3es21.14)", numsigma
      print"(a)", "diff:"
      print"(3es21.14)", sigma(:, :) - numsigma(:, :)
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

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
         call model%update(mol, error, cache, cn, qloc)
         if (allocated(error)) exit lp
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, error, cache, cn, qloc)
         if (allocated(error)) exit lp
         call model%get_xvec(mol, cache, xvecl)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat, :) = 0.5_wp * (xvecr(:) - xvecl(:)) / step
      end do
   end do lp
   if (allocated(error)) return

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, error, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(error)) return
   call model%get_xvec(mol, cache, xvecl) ! need to call this for xtmp in cache (eeqbc)
   call model%get_xvec_derivs(mol, cache, dbdr, dbdL)

   if (any(abs(dbdr(:, :, :) - numgrad(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the b vector does not match")
      print"(a)", "dbdr:"
      print"(3es21.14)", dbdr
      print"(a)", "numgrad:"
      print"(3es21.14)", numgrad
      print"(a)", "diff:"
      print"(3es21.14)", dbdr - numgrad
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
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
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

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
         call model%update(mol, error, cache, cn, qloc)
         if (allocated(error)) exit lp
         call model%get_xvec(mol, cache, xvecr)

         ! Left-hand side
         xvecl(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, error, cache, cn, qloc)
         if (allocated(error)) exit lp
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
   if (allocated(error)) return

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, error, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(error)) return
   call model%get_xvec(mol, cache, xvecl) ! need to call this for xtmp in cache (eeqbc)
   call model%get_xvec_derivs(mol, cache, dbdr, dbdL)

   if (any(abs(dbdL(:, :, :) - numsigma(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the b vector does not match")
      print"(a)", "dbdL:"
      print"(3es21.14)", dbdL
      print"(a)", "numsigma:"
      print"(3es21.14)", numsigma
      print"(a)", "diff:"
      print"(3es21.14)", dbdL - numsigma
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

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
         call model%update(mol, error, cache, cn, qloc)
         if (allocated(error)) exit lp
         call model%get_coulomb_matrix(mol, cache, amatr)

         ! Left-hand side
         amatl(:, :) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2 * step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%local_charge(mol, trans, qloc)
         call model%update(mol, error, cache, cn, qloc)
         if (allocated(error)) exit lp
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
   if (allocated(error)) return

   ! Analytical gradient
   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc, dqlocdr, dqlocdL)
   call model%update(mol, error, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(error)) return
   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)

   ! Add trace of the A matrix
   do iat = 1, mol%nat
      dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
   end do

   ! higher tolerance for numerical gradient
   if (any(abs(dadr(:, :, :) - numgrad(:, :, :)) > thr2_local)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print"(a)", "dadr:"
      print"(3es21.14)", dadr
      print"(a)", "numgrad:"
      print"(3es21.14)", numgrad
      print"(a)", "diff:"
      print"(3es21.14)", dadr - numgrad
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
   & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
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

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
         call model%update(mol, error, cache, cn, qloc)
         if (allocated(error)) exit lp
         call model%get_coulomb_matrix(mol, cache, amatr)

         amatl(:, :) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2 * step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%local_charge(mol, lattr, qloc)
         call model%update(mol, error, cache, cn, qloc)
         if (allocated(error)) exit lp
         call model%get_coulomb_matrix(mol, cache, amatl)

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
   call model%update(mol, error, cache, cn, qloc, dcndr, dcndL, dqlocdr, dqlocdL)
   if (allocated(error)) return

   call model%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)

   if (any(abs(dadL(:, :, :) - numsigma(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of the A matrix does not match")
      print"(a)", "dadL:"
      print"(3es21.14)", dadL
      print"(a)", "numsigma:"
      print"(3es21.14)", numsigma
      print"(a)", "diff:"
      print"(3es21.14)", dadL - numsigma
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

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
      print"(a)", "dqdr:"
      print"(3es21.14)", dqdr
      print"(a)", "numdr:"
      print"(3es21.14)", numdr
      print"(a)", "diff:"
      print"(3es21.14)", dqdr - numdr
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
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], [3, 3])
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), trans(:, :)
   real(wp), allocatable :: qloc(:), dqlocdr(:, :, :), dqlocdL(:, :, :)
   real(wp), allocatable :: qr(:), ql(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :), numdL(:, :, :)
   real(wp) :: eps(3, 3), lattice(3, 3)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans, error)
   if (allocated(error)) return

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
      print"(a)", "dqdL:"
      print"(3es21.14)", dqdL
      print"(a)", "numdL:"
      print"(3es21.14)", numdL
      print"(a)", "diff:"
      print"(3es21.14)", dqdL - numdL
   end if

end subroutine test_numdqdL

subroutine test_eeq_q_cyanamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: ref(40) = [&
      & 3.47273696479021E-1_wp, 3.47272705640270E-1_wp, 3.47276492903749E-1_wp, &
      & 3.47266922178364E-1_wp, 3.47288439460215E-1_wp, 3.47296528581096E-1_wp, &
      & 3.47299555433948E-1_wp, 3.47299047247046E-1_wp, 3.45885448217334E-1_wp, &
      & 3.45873984082476E-1_wp, 3.45887979518336E-1_wp, 3.45877189683493E-1_wp, &
      & 3.45902103370702E-1_wp, 3.45901899732038E-1_wp, 3.45900074770712E-1_wp, &
      & 3.45903507923248E-1_wp, 3.58994427905596E-1_wp, 3.58985684782718E-1_wp, &
      & 3.58998551488410E-1_wp, 3.58988652945645E-1_wp, 3.58988114752443E-1_wp, &
      & 3.58995533629608E-1_wp, 3.58994816857993E-1_wp, 3.58999419102342E-1_wp, &
      &-5.86636852376099E-1_wp, -5.86601030236713E-1_wp, -5.86624269538219E-1_wp, &
      &-5.86579192001313E-1_wp, -5.86653588719833E-1_wp, -5.86668553721606E-1_wp, &
      &-5.86666673919381E-1_wp, -5.86674064898699E-1_wp, -4.65527015747556E-1_wp, &
      &-4.65525637708548E-1_wp, -4.65544375954499E-1_wp, -4.65538830945934E-1_wp, &
      &-4.65527785578432E-1_wp, -4.65531681382053E-1_wp, -4.65525608429865E-1_wp, &
      &-4.65525615528052E-1_wp]

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
   ! Updated for accuracy-driven Ewald cutoffs and mctc-lib Wigner-Seitz tolerance.
   real(wp), parameter :: ref(24) = [&
      & 4.01820317641649E-1_wp, 4.01826970293715E-1_wp, 4.01778618929473E-1_wp, &
      & 4.01801620737465E-1_wp, 4.30133215169715E-1_wp, 4.30106241841601E-1_wp, &
      & 4.30086201257563E-1_wp, 4.30109482177696E-1_wp, 1.90391275024120E-1_wp, &
      & 1.90392904370207E-1_wp, 1.90395221896778E-1_wp, 1.90398527506007E-1_wp, &
      & 3.33798778194695E-1_wp, 3.33779134873330E-1_wp, 3.33806406791091E-1_wp, &
      & 3.33823483957305E-1_wp, -8.35501801552073E-1_wp, -8.35481140201098E-1_wp, &
      &-8.35408885793255E-1_wp, -8.35455735270739E-1_wp, -6.31069087329026E-1_wp, &
      &-6.31049654240391E-1_wp, -6.31075181828644E-1_wp, -6.31071707101098E-1_wp]

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

subroutine test_eeq_dadr_znoh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, znoh_num, znoh_xyz, lattice=znoh_lat, periodic=[.true., .true., .true.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeq_dadr_znoh

subroutine test_eeq_dadL_znoh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, znoh_num, znoh_xyz, lattice=znoh_lat, periodic=[.true., .true., .true.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeq_dadL_znoh

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

subroutine test_eeq_s_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeq_s_co2

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

subroutine test_eeqbc_dadr_znoh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, znoh_num, znoh_xyz, lattice=znoh_lat, periodic=[.true., .true., .true.])
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadr(error, mol, model)

end subroutine test_eeqbc_dadr_znoh

subroutine test_eeqbc_dadL_znoh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, znoh_num, znoh_xyz, lattice=znoh_lat, periodic=[.true., .true., .true.])
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_dadL(error, mol, model)

end subroutine test_eeqbc_dadL_znoh

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

subroutine test_eeqbc_s_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call get_structure(mol, "X23", "CO2")
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)

end subroutine test_eeqbc_s_co2

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


!> Run all finite-difference derivative checks for a given structure and model
subroutine check_all_derivs(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Electronegativity equilibration model
   class(mchrg_model_type), intent(in) :: model

   call test_dbdr(error, mol, model)
   if (allocated(error)) return
   call test_dbdL(error, mol, model)
   if (allocated(error)) return
   call test_dadr(error, mol, model)
   if (allocated(error)) return
   call test_dadL(error, mol, model)
   if (allocated(error)) return
   call test_numgrad(error, mol, model)
   if (allocated(error)) return
   call test_numsigma(error, mol, model)
   if (allocated(error)) return
   call test_numdqdr(error, mol, model)
   if (allocated(error)) return
   call test_numdqdL(error, mol, model)

end subroutine check_all_derivs


!> Test 1D periodic EEQ across all mask orientations and verify derivatives.
subroutine test_eeq_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(3) = [30, 8, 1]
   real(wp), parameter :: a_wire(3) = [5.3_wp, 0.6_wp, -0.5_wp]
   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.42_wp * 5.3_wp, 0.42_wp * 0.6_wp, -0.42_wp * 0.5_wp, &
      & 1.1_wp, 1.8_wp, 0.7_wp], [3, 3])
   logical, parameter :: masks(3, 3) = reshape([ &
      & .true., .false., .false., &
      & .false., .true., .false., &
      & .false., .false., .true.], [3, 3])

   type(error_type), allocatable :: expected
   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   real(wp) :: lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3), dqdr_ref(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   integer :: imask

   do imask = 1, 3
      ! Construct two lattices with identical wire vector and different nonperiodic columns
      lat1 = 0.0_wp
      lat1(:, imask) = a_wire
      lat1(:, mod(imask, 3) + 1) = [0.0_wp, 20.0_wp, 0.0_wp]
      lat1(:, mod(imask + 1, 3) + 1) = [0.0_wp, 0.0_wp, 20.0_wp]

      lat2 = 0.0_wp
      lat2(:, imask) = a_wire
      lat2(:, mod(imask, 3) + 1) = [15.0_wp, 200.0_wp, -30.0_wp]
      lat2(:, mod(imask + 1, 3) + 1) = [-40.0_wp, 10.0_wp, 250.0_wp]

      call new(mol, num, xyz, charge=0.0_wp, lattice=lat1, periodic=masks(:, imask))
      call new(mol_big, num, xyz, charge=0.0_wp, lattice=lat2, periodic=masks(:, imask))

      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return

      call get_charges(model, mol, error, q1, dqdr1, dqdL1)
      if (allocated(error)) return
      call get_charges(model, mol_big, error, q2, dqdr2, dqdL2)
      if (allocated(error)) return

      if (abs(sum(q1)) > thr) then
         call test_failed(error, "1D EEQ total charge constraint violated")
         return
      end if
      if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
         & .or. maxval(abs(dqdL1 - dqdL2)) > thr) then
         call test_failed(error, "1D EEQ depends on nonperiodic lattice vectors")
         return
      end if
      if (imask == 1) then
         q_ref = q1
         dqdr_ref = dqdr1
      else if (maxval(abs(q1 - q_ref)) > thr .or. maxval(abs(dqdr1 - dqdr_ref)) > thr) then
         call test_failed(error, "1D EEQ results depend on periodic mask index")
         return
      end if
   end do

   ! Validate analytical derivatives against finite differences once
   call check_all_derivs(error, mol, model)
   if (allocated(error)) return

   ! Charged 1D EEQ cell must be rejected
   mol%charge = 1.0_wp
   call get_charges(model, mol, expected, q1)
   if (.not. allocated(expected)) then
      call test_failed(error, "Charged 1D EEQ cell was not rejected")
      return
   end if

end subroutine test_eeq_1d


!> Test 2D periodic EEQ across all mask orientations and verify derivatives.
subroutine test_eeq_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(3) = [30, 8, 1]
   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp, -0.3_wp, &
      & 2.5_wp, 0.2_wp,  0.5_wp, &
      & 0.6_wp, 2.1_wp,  0.4_wp], [3, 3])
   real(wp), parameter :: a1(3) = [5.2_wp, 0.3_wp, 0.1_wp]
   real(wp), parameter :: a2(3) = [0.5_wp, 6.1_wp, -0.2_wp]
   logical, parameter :: masks(3, 3) = reshape([ &
      & .true., .true., .false., &
      & .true., .false., .true., &
      & .false., .true., .true.], [3, 3])
   integer, parameter :: p1_idx(3) = [1, 1, 2]
   integer, parameter :: p2_idx(3) = [2, 3, 3]
   integer, parameter :: np_idx(3) = [3, 2, 1]

   type(error_type), allocatable :: expected
   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   real(wp) :: lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3), dqdr_ref(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   integer :: imask

   do imask = 1, 3
      ! Construct two lattices with identical slab vectors and different nonperiodic column
      lat1 = 0.0_wp
      lat1(:, p1_idx(imask)) = a1
      lat1(:, p2_idx(imask)) = a2
      lat1(:, np_idx(imask)) = [0.0_wp, 0.0_wp, 20.0_wp]

      lat2 = 0.0_wp
      lat2(:, p1_idx(imask)) = a1
      lat2(:, p2_idx(imask)) = a2
      lat2(:, np_idx(imask)) = [35.0_wp, -70.0_wp, 200.0_wp]

      call new(mol, num, xyz, charge=0.0_wp, lattice=lat1, periodic=masks(:, imask))
      call new(mol_big, num, xyz, charge=0.0_wp, lattice=lat2, periodic=masks(:, imask))

      call new_eeq2019_model(mol, model, error)
      if (allocated(error)) return

      call get_charges(model, mol, error, q1, dqdr1, dqdL1)
      if (allocated(error)) return
      call get_charges(model, mol_big, error, q2, dqdr2, dqdL2)
      if (allocated(error)) return

      if (abs(sum(q1)) > thr) then
         call test_failed(error, "2D EEQ total charge constraint violated")
         return
      end if
      if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
         & .or. maxval(abs(dqdL1 - dqdL2)) > thr) then
         call test_failed(error, "2D EEQ depends on nonperiodic lattice vector")
         return
      end if
      if (imask == 1) then
         q_ref = q1
         dqdr_ref = dqdr1
      else if (maxval(abs(q1 - q_ref)) > thr .or. maxval(abs(dqdr1 - dqdr_ref)) > thr) then
         call test_failed(error, "2D EEQ results depend on periodic mask index")
         return
      end if
   end do

   ! Validate analytical derivatives against finite differences once
   call check_all_derivs(error, mol, model)
   if (allocated(error)) return

   ! Charged 2D EEQ cell must be rejected
   mol%charge = 1.0_wp
   call get_charges(model, mol, expected, q1)
   if (.not. allocated(expected)) then
      call test_failed(error, "Charged 2D EEQ cell was not rejected")
      return
   end if

end subroutine test_eeq_2d


!> Test 1D periodic EEQBC across all mask orientations (neutral and charged).
subroutine test_eeqbc_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(3) = [30, 8, 1]
   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp,  0.0_wp, &
      & 2.4_wp, 0.3_wp, -0.4_wp, &
      & 1.1_wp, 1.8_wp,  0.7_wp], [3, 3])
   real(wp), parameter :: a_wire(3) = [5.3_wp, 0.6_wp, -0.5_wp]
   logical, parameter :: masks(3, 3) = reshape([ &
      & .true., .false., .false., &
      & .false., .true., .false., &
      & .false., .false., .true.], [3, 3])
   real(wp), parameter :: charges(2) = [0.0_wp, 1.0_wp]

   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   real(wp) :: lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3), dqdr_ref(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   integer :: imask, iq

   do iq = 1, size(charges)
      do imask = 1, 3
         lat1 = 0.0_wp
         lat1(:, imask) = a_wire
         lat1(:, mod(imask, 3) + 1) = [0.0_wp, 20.0_wp, 0.0_wp]
         lat1(:, mod(imask + 1, 3) + 1) = [0.0_wp, 0.0_wp, 20.0_wp]

         lat2 = 0.0_wp
         lat2(:, imask) = a_wire
         lat2(:, mod(imask, 3) + 1) = [15.0_wp, 200.0_wp, -30.0_wp]
         lat2(:, mod(imask + 1, 3) + 1) = [-40.0_wp, 10.0_wp, 250.0_wp]

         call new(mol, num, xyz, charge=charges(iq), lattice=lat1, periodic=masks(:, imask))
         call new(mol_big, num, xyz, charge=charges(iq), lattice=lat2, periodic=masks(:, imask))

         call new_eeqbc2025_model(mol, model, error)
         if (allocated(error)) return

         call get_charges(model, mol, error, q1, dqdr1, dqdL1)
         if (allocated(error)) return
         call get_charges(model, mol_big, error, q2, dqdr2, dqdL2)
         if (allocated(error)) return

         if (abs(sum(q1) - charges(iq)) > thr) then
            call test_failed(error, "1D EEQBC total charge constraint violated")
            return
         end if
         if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
            & .or. maxval(abs(dqdL1 - dqdL2)) > thr) then
            call test_failed(error, "1D EEQBC depends on nonperiodic lattice vectors")
            return
         end if
         if (imask == 1) then
            q_ref = q1
            dqdr_ref = dqdr1
         else if (maxval(abs(q1 - q_ref)) > thr .or. maxval(abs(dqdr1 - dqdr_ref)) > thr) then
            call test_failed(error, "1D EEQBC results depend on periodic mask index")
            return
         end if
      end do

      ! Validate analytical derivatives against finite differences once per charge state
      call check_all_derivs(error, mol, model)
      if (allocated(error)) return
   end do

end subroutine test_eeqbc_1d


!> Test 2D periodic EEQBC across all mask orientations (neutral and charged).
subroutine test_eeqbc_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(3) = [30, 8, 1]
   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp, -0.3_wp, &
      & 2.5_wp, 0.2_wp,  0.5_wp, &
      & 0.6_wp, 2.1_wp,  0.4_wp], [3, 3])
   real(wp), parameter :: a1(3) = [5.2_wp, 0.3_wp, 0.1_wp]
   real(wp), parameter :: a2(3) = [0.5_wp, 6.1_wp, -0.2_wp]
   logical, parameter :: masks(3, 3) = reshape([ &
      & .true., .true., .false., &
      & .true., .false., .true., &
      & .false., .true., .true.], [3, 3])
   integer, parameter :: p1_idx(3) = [1, 1, 2]
   integer, parameter :: p2_idx(3) = [2, 3, 3]
   integer, parameter :: np_idx(3) = [3, 2, 1]
   real(wp), parameter :: charges(2) = [0.0_wp, -1.0_wp]

   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   real(wp) :: lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3), dqdr_ref(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   integer :: imask, iq

   do iq = 1, size(charges)
      do imask = 1, 3
         lat1 = 0.0_wp
         lat1(:, p1_idx(imask)) = a1
         lat1(:, p2_idx(imask)) = a2
         lat1(:, np_idx(imask)) = [0.0_wp, 0.0_wp, 20.0_wp]

         lat2 = 0.0_wp
         lat2(:, p1_idx(imask)) = a1
         lat2(:, p2_idx(imask)) = a2
         lat2(:, np_idx(imask)) = [35.0_wp, -70.0_wp, 200.0_wp]

         call new(mol, num, xyz, charge=charges(iq), lattice=lat1, periodic=masks(:, imask))
         call new(mol_big, num, xyz, charge=charges(iq), lattice=lat2, periodic=masks(:, imask))

         call new_eeqbc2025_model(mol, model, error)
         if (allocated(error)) return

         call get_charges(model, mol, error, q1, dqdr1, dqdL1)
         if (allocated(error)) return
         call get_charges(model, mol_big, error, q2, dqdr2, dqdL2)
         if (allocated(error)) return

         if (abs(sum(q1) - charges(iq)) > thr) then
            call test_failed(error, "2D EEQBC total charge constraint violated")
            return
         end if
         if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
            & .or. maxval(abs(dqdL1 - dqdL2)) > thr) then
            call test_failed(error, "2D EEQBC depends on nonperiodic lattice vector")
            return
         end if
         if (imask == 1) then
            q_ref = q1
            dqdr_ref = dqdr1
         else if (maxval(abs(q1 - q_ref)) > thr .or. maxval(abs(dqdr1 - dqdr_ref)) > thr) then
            call test_failed(error, "2D EEQBC results depend on periodic mask index")
            return
         end if
      end do

      ! Validate analytical derivatives against finite differences once per charge state
      call check_all_derivs(error, mol, model)
      if (allocated(error)) return
   end do

end subroutine test_eeqbc_2d


!> Verify one-element periodicity mask handling and degenerate-lattice error propagation.
subroutine test_ewald_setup(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: lat_ref(3, 3) = reshape([ &
      & 3.6_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 20.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 20.0_wp], [3, 3])

   type(error_type), allocatable :: expected
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :), trans(:, :)
   real(wp) :: lat_deg(3, 3), alpha_1, alpha_3

   call get_alpha(lat_ref, alpha_1, error, periodic=[.true.])
   if (allocated(error)) return

   call get_alpha(lat_ref, alpha_3, error, periodic=[.true., .true., .true.])
   if (allocated(error)) return

   if (abs(alpha_1 - alpha_3) > thr) then
      call test_failed(error, "get_alpha with 1-element mask differs from 3-element mask")
      return
   end if

   lat_deg = 0.0_wp
   call get_alpha(lat_deg, alpha_1, expected, periodic=[.true., .false., .false.])
   if (.not. allocated(expected)) then
      call test_failed(error, "get_alpha did not report error for degenerate lattice")
      return
   end if

   deallocate(expected)
   call get_ewald_trans([.true., .false., .false.], lat_deg, 0.25_wp, 0.2_wp, &
      & dtrans, rtrans, expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "get_ewald_trans did not report error for degenerate lattice")
      return
   end if

   deallocate(expected)
   call get_ewald_dir_trans([.true., .false., .false.], lat_deg, 0.25_wp, 0.2_wp, &
      & trans, expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "get_ewald_dir_trans did not report error for degenerate lattice")
      return
   end if

   deallocate(expected)
   call get_ewald_rec_trans([.true., .false., .false.], lat_deg, 0.25_wp, trans, expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "get_ewald_rec_trans did not report error for degenerate lattice")
      return
   end if

end subroutine test_ewald_setup


!> Verify Ewald splitting parameter invariance for 1D periodic EEQ.
subroutine test_ewald_alpha_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz_1d(3, 3) = reshape([ &
      & 0.0_wp,  0.0_wp,  0.0_wp, &
      & 1.1_wp,  0.8_wp, -0.3_wp, &
      & 0.5_wp, -0.9_wp,  0.6_wp], [3, 3])
   real(wp), parameter :: lat_1d(3, 3) = reshape([ &
      & 3.6_wp,  0.0_wp,  0.0_wp, &
      & 0.0_wp, 20.0_wp,  0.0_wp, &
      & 0.0_wp,  0.0_wp, 20.0_wp], [3, 3])

   type(structure_type) :: mol

   call new(mol, [8, 1, 6], xyz_1d, charge=0.0_wp, lattice=lat_1d, &
      & periodic=[.true., .false., .false.])
   call check_alpha_invariance(error, mol)

end subroutine test_ewald_alpha_1d


!> Verify Ewald splitting parameter invariance for 2D periodic EEQ.
subroutine test_ewald_alpha_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz_2d(3, 3) = reshape([ &
      & 0.1_wp, 0.2_wp, -0.1_wp, &
      & 1.7_wp, 0.6_wp,  0.3_wp, &
      & 0.3_wp, 1.8_wp,  0.5_wp], [3, 3])
   real(wp), parameter :: lat_2d(3, 3) = reshape([ &
      & 5.5_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 0.5_wp, 6.0_wp, 2.0_wp], [3, 3])

   type(structure_type) :: mol

   call new(mol, [8, 1, 1], xyz_2d, charge=0.0_wp, lattice=lat_2d, &
      & periodic=[.true., .false., .true.])
   call check_alpha_invariance(error, mol)

end subroutine test_ewald_alpha_2d


!> Check Ewald alpha-splitting invariance of Coulomb matrix, first derivatives, and Hessian
subroutine check_alpha_invariance(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Periodic structure (3 atoms)
   type(structure_type), intent(in) :: mol

   real(wp), parameter :: alphas(3) = [0.16_wp, 0.26_wp, 0.40_wp]
   real(wp), parameter :: qvec(4) = [0.35_wp, -0.55_wp, 0.20_wp, 0.05_wp]
   real(wp), parameter :: adj(4) = [-0.40_wp, 0.25_wp, 0.15_wp, 0.0_wp]

   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: cache
   real(wp), allocatable :: trans(:, :)
   real(wp) :: cn(3), qloc(3), gam_min
   real(wp) :: amat_ref(4, 4), amat_cur(4, 4)
   real(wp) :: dadr_ref(3, 3, 4), dadL_ref(3, 3, 4), atrace_ref(3, 3)
   real(wp) :: dadr_cur(3, 3, 4), dadL_cur(3, 3, 4), atrace_cur(3, 3)
   real(wp) :: hess_ref(9, 9), hess_cur(9, 9)
   integer :: ia

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   call get_lattice_points(mol%periodic, mol%lattice, 25.0_wp, trans, error)
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%update(mol, error, cache, cn, qloc)
   if (allocated(error)) return

   allocate(cache%qvec(4), source=qvec)
   select type (raw => cache%raw)
   type is (eeq_cache)
      allocate(raw%dcndr(3, 3, 3), source=0.0_wp)
   end select

   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(model%rad(mol%id)))

   do ia = 1, size(alphas)
      select type (raw => cache%raw)
      type is (eeq_cache)
         raw%alpha = alphas(ia)
         call get_ewald_trans(mol%periodic, mol%lattice, raw%alpha, gam_min, &
            & raw%dtrans, raw%rtrans, error)
         if (allocated(error)) return
      end select

      call model%get_coulomb_matrix(mol, cache, amat_cur)
      call model%get_coulomb_derivs(mol, cache, qvec, dadr_cur, dadL_cur, atrace_cur)
      hess_cur = 0.0_wp
      call model%add_model_hessian(mol, error, cache, adj, hess_cur)
      if (allocated(error)) return

      if (ia == 1) then
         amat_ref = amat_cur
         dadr_ref = dadr_cur
         dadL_ref = dadL_cur
         atrace_ref = atrace_cur
         hess_ref = hess_cur
      else
         if (maxval(abs(amat_cur - amat_ref)) > 1.0e-11_wp &
            & .or. maxval(abs(dadr_cur - dadr_ref)) > 1.0e-11_wp &
            & .or. maxval(abs(dadL_cur - dadL_ref)) > 1.0e-11_wp &
            & .or. maxval(abs(atrace_cur - atrace_ref)) > 1.0e-11_wp &
            & .or. maxval(abs(hess_cur - hess_ref)) > 1.0e-10_wp) then
            call test_failed(error, "Periodic Ewald sum is not invariant to alpha")
            return
         end if
      end if
   end do

end subroutine check_alpha_invariance


!> Verify invariance under unimodular basis shear (a2' = a2 + a1) and multi-cell shifts.
subroutine test_ewald_lattice_invariance(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(3) = [8, 1, 1]
   real(wp), parameter :: xyz_2d(3, 3) = reshape([ &
      & 0.1_wp, 0.2_wp, -0.1_wp, &
      & 1.7_wp, 0.6_wp,  0.3_wp, &
      & 0.3_wp, 1.8_wp,  0.5_wp], [3, 3])
   real(wp), parameter :: a1_2d(3) = [5.5_wp, 0.0_wp, 0.0_wp]
   real(wp), parameter :: a2_2d(3) = [0.5_wp, 6.0_wp, 2.0_wp]
   real(wp), parameter :: weights(3) = [0.6_wp, -0.4_wp, 0.25_wp]

   type(structure_type) :: mol, mol_shift
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: cache
   real(wp) :: lat_ref(3, 3), lat_skew(3, 3), xyz_shift(3, 3)
   real(wp) :: q_ref(3), q_shift(3), dqdr_ref(3, 3, 3), dqdr_shift(3, 3, 3)
   real(wp) :: dqdL_ref(3, 3, 3), dqdL_shift(3, 3, 3)
   real(wp) :: hess_ref(9, 9), hess_shift(9, 9)

   lat_ref = 0.0_wp
   lat_ref(:, 1) = a1_2d
   lat_ref(:, 3) = a2_2d

   lat_skew = 0.0_wp
   lat_skew(:, 1) = a1_2d
   lat_skew(:, 3) = a2_2d + a1_2d

   xyz_shift = xyz_2d
   xyz_shift(:, 2) = xyz_2d(:, 2) + a1_2d
   xyz_shift(:, 3) = xyz_2d(:, 3) - a2_2d

   call new(mol, num, xyz_2d, charge=0.0_wp, lattice=lat_ref, &
      & periodic=[.true., .false., .true.])
   call new(mol_shift, num, xyz_shift, charge=0.0_wp, lattice=lat_skew, &
      & periodic=[.true., .false., .true.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   call get_charges(model, mol, error, q_ref, dqdr_ref, dqdL_ref, cache=cache)
   if (allocated(error)) return

   hess_ref = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hess_ref, cache)
   if (allocated(error)) return

   call get_charges(model, mol_shift, error, q_shift, dqdr_shift, dqdL_shift, cache=cache)
   if (allocated(error)) return

   hess_shift = 0.0_wp
   call add_charge_hessian(model, mol_shift, error, weights, hess_shift, cache)
   if (allocated(error)) return

   if (maxval(abs(q_shift - q_ref)) > 1.0e-11_wp &
      & .or. maxval(abs(dqdr_shift - dqdr_ref)) > 1.0e-10_wp &
      & .or. maxval(abs(dqdL_shift - dqdL_ref)) > 1.0e-10_wp &
      & .or. maxval(abs(hess_shift - hess_ref)) > 1.0e-9_wp) then
      call test_failed(error, &
         & "2D EEQ not invariant under unimodular transform and multi-cell shift")
   end if

end subroutine test_ewald_lattice_invariance


!> Regression test for DFT-D4 issue #285: 2D slab reference charges across cell shapes.
subroutine test_eeq_slab_285(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(3) = [8, 1, 1]
   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & -0.3910542064024872_wp,  0.3976645623910747_wp, 2.977929077294636_wp, &
      & -0.3542054272959850_wp,  1.0567259314317760_wp, 2.610987477729398_wp, &
      & -0.9445722089116435_wp, -0.0304014548863456_wp, 3.370602014356701_wp], [3, 3]) * aatoau
   real(wp), parameter :: q2d_hex_ref(3) = [ &
      & -1.06091408_wp, 0.54516560_wp, 0.51574848_wp]
   real(wp), parameter :: q2d_ahex_ref(3) = [ &
      & -1.06069500_wp, 0.54509370_wp, 0.51560131_wp]
   real(wp), parameter :: q2d_orth_ref(3) = [ &
      & -1.06073102_wp, 0.54653780_wp, 0.51419322_wp]

   real(wp) :: lat(3, 3)

   lat = 0.0_wp
   lat(:, 1) = [2.794646093662685_wp, 0.0_wp, 0.0_wp] * aatoau
   lat(:, 2) = [-1.397323046831344_wp, 2.420234511698831_wp, 0.0_wp] * aatoau
   call check_slab_charges(error, num, xyz, lat, q2d_hex_ref, "hexagonal")
   if (allocated(error)) return

   lat(:, 1) = [2.8_wp, 0.0_wp, 0.0_wp] * aatoau
   call check_slab_charges(error, num, xyz, lat, q2d_ahex_ref, "almost-hexagonal")
   if (allocated(error)) return

   lat(:, 2) = [0.0_wp, 2.8_wp, 0.0_wp] * aatoau
   call check_slab_charges(error, num, xyz, lat, q2d_orth_ref, "orthorhombic")

end subroutine test_eeq_slab_285


!> Check 2D slab EEQ charges against a reference vector
subroutine check_slab_charges(error, num, xyz, lat, q_ref, label)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Atomic numbers
   integer, intent(in) :: num(:)

   !> Cartesian coordinates in Bohr
   real(wp), intent(in) :: xyz(:, :)

   !> Lattice vectors in Bohr
   real(wp), intent(in) :: lat(3, 3)

   !> Reference partial charges
   real(wp), intent(in) :: q_ref(:)

   !> Cell shape label
   character(len=*), intent(in) :: label

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp) :: q(size(num))

   call new(mol, num, xyz, lattice=lat, periodic=[.true., .true., .false.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   call get_charges(model, mol, error, q)
   if (allocated(error)) return

   if (maxval(abs(q - q_ref)) > 1.0e-7_wp) then
      call test_failed(error, "2D " // label // " #285 slab charges deviate from reference")
   end if

end subroutine check_slab_charges


!> Verify 3D large-vacuum slab convergence (1/Lz dipole-sheet scaling) toward 2D EEQ.
subroutine test_eeq_slab_3d_limit(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(3) = [8, 1, 1]
   real(wp), parameter :: xyz(3, 3) = reshape([ &
      & -0.3910542064024872_wp,  0.3976645623910747_wp, 2.977929077294636_wp, &
      & -0.3542054272959850_wp,  1.0567259314317760_wp, 2.610987477729398_wp, &
      & -0.9445722089116435_wp, -0.0304014548863456_wp, 3.370602014356701_wp], [3, 3]) * aatoau

   type(structure_type) :: mol2d, mol3d
   class(mchrg_model_type), allocatable :: model
   real(wp) :: lat(3, 3), q2d(3), q_z50(3), q_z100(3)

   lat = 0.0_wp
   lat(:, 1) = [2.794646093662685_wp, 0.0_wp, 0.0_wp] * aatoau
   lat(:, 2) = [-1.397323046831344_wp, 2.420234511698831_wp, 0.0_wp] * aatoau
   lat(:, 3) = [0.0_wp, 0.0_wp, 50.0_wp] * aatoau

   call new(mol2d, num, xyz, lattice=lat, periodic=[.true., .true., .false.])
   call new_eeq2019_model(mol2d, model, error)
   if (allocated(error)) return

   call get_charges(model, mol2d, error, q2d)
   if (allocated(error)) return

   call new(mol3d, num, xyz, lattice=lat, periodic=[.true., .true., .true.])
   call get_charges(model, mol3d, error, q_z50)
   if (allocated(error)) return

   lat(:, 3) = [0.0_wp, 0.0_wp, 100.0_wp] * aatoau
   call new(mol3d, num, xyz, lattice=lat, periodic=[.true., .true., .true.])
   call get_charges(model, mol3d, error, q_z100)
   if (allocated(error)) return

   if (maxval(abs(q_z50 - q2d)) > 3.0e-5_wp &
      & .or. maxval(abs(q_z100 - q2d)) > 1.5e-5_wp &
      & .or. maxval(abs(50.0_wp * (q_z50 - q2d) - 100.0_wp * (q_z100 - q2d))) > 1.0e-5_wp) then
      call test_failed(error, "3D slab charges failed to converge toward 2D EEQ limit")
   end if

end subroutine test_eeq_slab_3d_limit


!> Verify 1D and 2D Ewald special-function accuracy against 50-digit references.
subroutine test_special_functions(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: u_pts(6) = [ &
      & 0.0_wp, 1.0e-8_wp, 0.5_wp, 4.0_wp, 4.000001_wp, 25.0_wp]
   real(wp), parameter :: f0_ref(6) = [ &
      & 0.0_wp, 9.9999999750000000556e-9_wp, 0.44384207911774836294_wp, &
      & 1.9672893784312723859_wp, 1.9672896238523342756_wp, 3.7960914897702684998_wp]
   real(wp), parameter :: f1_ref(6) = [ &
      & 1.0_wp, 0.99999999500000001667_wp, 0.78693868057473315279_wp, &
      & 0.24542109027781645493_wp, 0.24542103350146551229_wp, 0.039999999999444482245_wp]
   real(wp), parameter :: f2_ref(6) = [ &
      & -0.5_wp, -0.49999999666666667917_wp, -0.36081604172419945838_wp, &
      & -0.056776362847270568658_wp, -0.056776339038005506227_wp, -0.0015999999994222615352_wp]

   real(wp), parameter :: x_pts(5) = [1.0e-8_wp, 1.0e-4_wp, 0.5_wp, 2.0_wp, 10.0_wp]
   real(wp), parameter :: y_pts(5) = [0.0_wp, 1.0e-2_wp, 0.5_wp, 0.1_wp, 25.0_wp]
   real(wp), parameter :: k0_ref(5) = [ &
      & 17.843465089050832587_wp, 8.6232592773087082697_wp, 0.42102443824070833334_wp, &
      & 0.045293668272415987242_wp, 8.1625283072617411727e-15_wp]
   real(wp), parameter :: k1_ref(5) = [ &
      & 0.99999981156534915949_wp, 0.99405430461170725993_wp, 0.23402778902579225314_wp, &
      & 0.034642560210863181845_wp, 5.2280069047533075135e-15_wp]
   real(wp), parameter :: k2_ref(5) = [ &
      & 0.49999999000000096717_wp, 0.49657968235608294261_wp, 0.15332113394940819643_wp, &
      & 0.027734685027132461108_wp, 3.4489111320542408121e-15_wp]

   real(wp), parameter :: g_2d(3) = [2.0_wp, 0.5_wp, 5.0_wp]
   real(wp), parameter :: s_2d(3) = [-1.0_wp, -5.0_wp, 20.0_wp]
   real(wp), parameter :: v2d_ref(3) = [ &
      & 0.00005300645451857405009_wp, &
      & 0.15243208859081749745_wp, &
      & 1.4504329021329886148e-45_wp]

   real(wp) :: f0, f1, f2, k0, k1, k2, v_2d
   integer :: i

   do i = 1, size(u_pts)
      call fast_f012_1d(u_pts(i), f0, f1, f2)
      if (abs(f0 - f0_ref(i)) > 1.0e-14_wp * max(1.0_wp, abs(f0_ref(i))) &
         & .or. abs(f1 - f1_ref(i)) > 1.0e-14_wp * max(1.0_wp, abs(f1_ref(i))) &
         & .or. abs(f2 - f2_ref(i)) > 1.0e-14_wp * max(1.0_wp, abs(f2_ref(i)))) then
         call test_failed(error, "fast_f012_1d deviates from 50-digit reference")
         return
      end if
   end do

   do i = 1, size(x_pts)
      call fast_k012_1d(x_pts(i), y_pts(i), k0, k1, k2)
      if (abs(k0 - k0_ref(i)) > 2.0e-14_wp * max(1.0e-15_wp, abs(k0_ref(i))) &
         & .or. abs(k1 - k1_ref(i)) > 2.0e-14_wp * max(1.0e-15_wp, abs(k1_ref(i))) &
         & .or. abs(k2 - k2_ref(i)) > 2.0e-14_wp * max(1.0e-15_wp, abs(k2_ref(i)))) then
         call test_failed(error, "fast_k012_1d deviates from 50-digit reference")
         return
      end if
   end do

   do i = 1, size(g_2d)
      v_2d = erfc_exp_2d(g_2d(i), s_2d(i), 0.35_wp)
      if (abs(v_2d - v2d_ref(i)) > 5.0e-14_wp * v2d_ref(i)) then
         call test_failed(error, "erfc_exp_2d deviates from 50-digit reference")
         return
      end if
   end do

end subroutine test_special_functions


!> Verify linear energy scaling and charge invariance for 1D elongated supercells.
subroutine test_supercell_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: lat_1d(3, 3) = reshape([ &
      & 5.3_wp,  0.6_wp, -0.5_wp, &
      & 0.0_wp, 20.0_wp,  0.0_wp, &
      & 0.0_wp,  0.0_wp, 20.0_wp], [3, 3])
   integer, parameter :: reps(3) = [2, 3, 4]

   type(structure_type) :: mol
   integer :: i

   call new(mol, znoh_num, znoh_xyz, charge=0.0_wp, lattice=lat_1d, &
      & periodic=[.true., .false., .false.])

   do i = 1, size(reps)
      call check_supercell_scaling(error, mol, .false., reps(i), 1, 1, thr)
      if (allocated(error)) return
   end do

   do i = 1, 2
      call check_supercell_scaling(error, mol, .true., reps(i), 1, 1, thr2)
      if (allocated(error)) return
   end do

   mol%charge = 1.0_wp
   do i = 1, 2
      call check_supercell_scaling(error, mol, .true., reps(i), 1, 1, thr2)
      if (allocated(error)) return
   end do

end subroutine test_supercell_1d


!> Verify linear energy scaling and charge invariance for 2D elongated supercells.
subroutine test_supercell_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: lat_2d(3, 3) = reshape([ &
      & 5.2_wp, 0.3_wp,  0.1_wp, &
      & 0.5_wp, 6.1_wp, -0.2_wp, &
      & 0.0_wp, 0.0_wp, 25.0_wp], [3, 3])
   integer, parameter :: n1_eeq(4) = [2, 1, 4, 2]
   integer, parameter :: n2_eeq(4) = [1, 3, 1, 3]

   type(structure_type) :: mol
   integer :: i

   call new(mol, znoh_num, znoh_xyz, charge=0.0_wp, lattice=lat_2d, &
      & periodic=[.true., .true., .false.])

   do i = 1, size(n1_eeq)
      call check_supercell_scaling(error, mol, .false., n1_eeq(i), n2_eeq(i), 1, thr)
      if (allocated(error)) return
   end do

   call check_supercell_scaling(error, mol, .true., 3, 1, 1, thr2)
   if (allocated(error)) return

   call check_supercell_scaling(error, mol, .true., 2, 2, 1, thr2)
   if (allocated(error)) return

   mol%charge = -1.0_wp
   call check_supercell_scaling(error, mol, .true., 2, 2, 1, thr2)

end subroutine test_supercell_2d


!> Verify linear energy scaling and charge invariance for 3D elongated supercells.
subroutine test_supercell_3d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n1_eeq(3) = [3, 1, 2]
   integer, parameter :: n2_eeq(3) = [1, 1, 3]
   integer, parameter :: n3_eeq(3) = [1, 4, 1]

   type(structure_type) :: mol
   integer :: i

   call new(mol, znoh_num, znoh_xyz, charge=0.0_wp, lattice=znoh_lat, &
      & periodic=[.true., .true., .true.])

   do i = 1, size(n1_eeq)
      call check_supercell_scaling(error, mol, .false., n1_eeq(i), n2_eeq(i), &
         & n3_eeq(i), thr)
      if (allocated(error)) return
   end do

   call check_supercell_scaling(error, mol, .true., 2, 1, 1, thr2)
   if (allocated(error)) return

   call check_supercell_scaling(error, mol, .true., 2, 2, 1, thr2)
   if (allocated(error)) return

   mol%charge = 1.0_wp
   call check_supercell_scaling(error, mol, .true., 2, 2, 1, thr2)

end subroutine test_supercell_3d


!> Check linear extensivity of total energy and invariance of atomic energies/charges
subroutine check_supercell_scaling(error, mol, use_eeqbc, n1, n2, n3, tol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Unit-cell structure
   type(structure_type), intent(in) :: mol

   !> Select EEQ-BC (.true.) or EEQ (.false.)
   logical, intent(in) :: use_eeqbc

   !> Supercell replication factor along axis 1
   integer, intent(in) :: n1

   !> Supercell replication factor along axis 2
   integer, intent(in) :: n2

   !> Supercell replication factor along axis 3
   integer, intent(in) :: n3

   !> Per-unit-cell comparison tolerance
   real(wp), intent(in) :: tol

   type(structure_type) :: sc
   real(wp), allocatable :: e_unit(:), q_unit(:), e_sc(:), q_sc(:)
   real(wp) :: etot_unit, etot_sc
   integer :: ntot, nat, irep, offset

   call solve_cell(mol, use_eeqbc, e_unit, q_unit, error)
   if (allocated(error)) return

   call make_supercell(mol, n1, n2, n3, sc)
   call solve_cell(sc, use_eeqbc, e_sc, q_sc, error)
   if (allocated(error)) return

   ntot = n1 * n2 * n3
   nat = mol%nat
   etot_unit = sum(e_unit)
   etot_sc = sum(e_sc)

   if (abs(etot_sc - real(ntot, wp) * etot_unit) > real(ntot, wp) * tol) then
      call test_failed(error, "Supercell total energy does not scale linearly")
      print "(a,3es21.14)", "E_sc, N*E_unit, diff:", etot_sc, &
         & real(ntot, wp) * etot_unit, etot_sc - real(ntot, wp) * etot_unit
      return
   end if

   do irep = 0, ntot - 1
      offset = irep * nat
      if (any(abs(e_sc(offset + 1:offset + nat) - e_unit) > tol) &
         & .or. any(abs(q_sc(offset + 1:offset + nat) - q_unit) > tol)) then
         call test_failed(error, "Supercell atomic energies or charges differ from unit cell")
         return
      end if
   end do

end subroutine check_supercell_scaling


!> Replicate a periodic structure into an (n1 x n2 x n3) supercell
subroutine make_supercell(mol, n1, n2, n3, sc)

   !> Unit-cell structure
   type(structure_type), intent(in) :: mol

   !> Replication factor along lattice vector 1
   integer, intent(in) :: n1

   !> Replication factor along lattice vector 2
   integer, intent(in) :: n2

   !> Replication factor along lattice vector 3
   integer, intent(in) :: n3

   !> Replicated supercell structure
   type(structure_type), intent(out) :: sc

   integer :: ntot, nat_sc, i1, i2, i3, ia, idx
   integer, allocatable :: num_sc(:)
   real(wp), allocatable :: xyz_sc(:, :)
   real(wp) :: lat_sc(3, 3), shift(3)

   ntot = n1 * n2 * n3
   nat_sc = mol%nat * ntot
   allocate(num_sc(nat_sc), xyz_sc(3, nat_sc))

   idx = 0
   do i1 = 0, n1 - 1
      do i2 = 0, n2 - 1
         do i3 = 0, n3 - 1
            shift = real(i1, wp) * mol%lattice(:, 1) &
               & + real(i2, wp) * mol%lattice(:, 2) &
               & + real(i3, wp) * mol%lattice(:, 3)
            do ia = 1, mol%nat
               idx = idx + 1
               num_sc(idx) = mol%num(mol%id(ia))
               xyz_sc(:, idx) = mol%xyz(:, ia) + shift
            end do
         end do
      end do
   end do

   lat_sc(:, 1) = real(n1, wp) * mol%lattice(:, 1)
   lat_sc(:, 2) = real(n2, wp) * mol%lattice(:, 2)
   lat_sc(:, 3) = real(n3, wp) * mol%lattice(:, 3)

   call new(sc, num_sc, xyz_sc, charge=real(ntot, wp) * mol%charge, &
      & lattice=lat_sc, periodic=mol%periodic)

end subroutine make_supercell


!> Evaluate atomic energies and partial charges for a structure
subroutine solve_cell(mol, use_eeqbc, energy, qvec, error)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Select EEQ-BC (.true.) or EEQ (.false.)
   logical, intent(in) :: use_eeqbc

   !> Atomic electrostatic energies
   real(wp), allocatable, intent(out) :: energy(:)

   !> Atomic partial charges
   real(wp), allocatable, intent(out) :: qvec(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   class(mchrg_model_type), allocatable :: model
   real(wp), allocatable :: trans(:, :), cn(:), qloc(:)

   if (use_eeqbc) then
      call new_eeqbc2025_model(mol, model, error)
   else
      call new_eeq2019_model(mol, model, error)
   end if
   if (allocated(error)) return

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, &
      & trans, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), qloc(mol%nat), energy(mol%nat), qvec(mol%nat))
   energy = 0.0_wp
   call model%ncoord%get_coordination_number(mol, trans, cn)
   call model%local_charge(mol, trans, qloc)
   call model%solve(mol, error, cn, qloc, energy=energy, qvec=qvec)

end subroutine solve_cell


end module test_pbc
