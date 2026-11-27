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
   use mctc_io_constants, only: pi
   use mctc_io_convert, only: aatoau
   use mctc_io_structure, only: structure_type, new
   use mstore, only: get_structure
   use multicharge, only: get_charges, add_charge_hessian
   use multicharge_ewald, only: get_alpha, get_ewald_setup
   use multicharge_ewald_cutoff, only: get_ewald_trans
   use multicharge_ewald_special, only: erfc_exp_2d, expint_e1, fast_f012_1d, &
      & fast_k012_1d
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
      & new_unittest("eeq-dadr-ice", test_eeq_dadr_ice), &
      & new_unittest("eeq-dadL-ice", test_eeq_dadL_ice), &
      & new_unittest("eeq-gradient-co2", test_eeq_g_co2), &
      & new_unittest("eeq-sigma-ice", test_eeq_s_ice), &
      & new_unittest("eeq-dqdr-urea", test_eeq_dqdr_urea), &
      & new_unittest("eeq-dqdL-oxacb", test_eeq_dqdL_oxacb), &
      & new_unittest("eeqbc-dbdr-co2", test_eeqbc_dbdr_co2), &
      & new_unittest("eeqbc-dbdL-co2", test_eeqbc_dbdL_co2), &
      & new_unittest("eeqbc-dadr-ice", test_eeqbc_dadr_ice), &
      & new_unittest("eeqbc-response-znoh", test_eeqbc_response_znoh), &
      & new_unittest("eeqbc-dadL-ice", test_eeqbc_dadL_ice), &
      & new_unittest("eeqbc-gradient-co2", test_eeqbc_g_co2), &
      & new_unittest("eeqbc-sigma-ice", test_eeqbc_s_ice), &
      & new_unittest("eeqbc-dqdr-urea", test_eeqbc_dqdr_urea), &
      & new_unittest("eeqbc-dqdL-oxacb", test_eeqbc_dqdL_oxacb), &
      & new_unittest("eeq-1d", test_eeq_1d), &
      & new_unittest("eeq-2d", test_eeq_2d), &
      & new_unittest("eeqbc-1d", test_eeqbc_1d), &
      & new_unittest("eeqbc-2d", test_eeqbc_2d), &
      & new_unittest("ewald-setup", test_ewald_setup), &
      & new_unittest("ewald-alpha-1d", test_ewald_alpha_1d), &
      & new_unittest("ewald-alpha-2d", test_ewald_alpha_2d), &
      & new_unittest("ewald-cutoff-convergence", test_ewald_cutoff_convergence), &
      & new_unittest("ewald-lattice-invariance", test_ewald_lattice_invariance), &
      & new_unittest("anisotropic-skewed-pbc", test_anisotropic_skewed), &
      & new_unittest("cache-invalidation", test_cache_invalidation), &
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
      & 3.47273696479028E-1_wp, 3.47272705640278E-1_wp, 3.47276492903757E-1_wp, &
      & 3.47266922178372E-1_wp, 3.47288439460222E-1_wp, 3.47296528581104E-1_wp, &
      & 3.47299555433954E-1_wp, 3.47299047247052E-1_wp, 3.45885448217350E-1_wp, &
      & 3.45873984082493E-1_wp, 3.45887979518354E-1_wp, 3.45877189683509E-1_wp, &
      & 3.45902103370717E-1_wp, 3.45901899732055E-1_wp, 3.45900074770728E-1_wp, &
      & 3.45903507923264E-1_wp, 3.58994427905607E-1_wp, 3.58985684782732E-1_wp, &
      & 3.58998551488422E-1_wp, 3.58988652945658E-1_wp, 3.58988114752453E-1_wp, &
      & 3.58995533629615E-1_wp, 3.58994816858005E-1_wp, 3.58999419102355E-1_wp, &
      &-5.86636852376127E-1_wp, -5.86601030236743E-1_wp, -5.86624269538247E-1_wp, &
      &-5.86579192001341E-1_wp, -5.86653588719858E-1_wp, -5.86668553721632E-1_wp, &
      &-5.86666673919406E-1_wp, -5.86674064898724E-1_wp, -4.65527015747563E-1_wp, &
      &-4.65525637708556E-1_wp, -4.65544375954505E-1_wp, -4.65538830945940E-1_wp, &
      &-4.65527785578439E-1_wp, -4.65531681382055E-1_wp, -4.65525608429872E-1_wp, &
      &-4.65525615528058E-1_wp]

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
      & 4.01820317641654E-1_wp, 4.01826970293720E-1_wp, 4.01778618929479E-1_wp, &
      & 4.01801620737470E-1_wp, 4.30133215169721E-1_wp, 4.30106241841606E-1_wp, &
      & 4.30086201257569E-1_wp, 4.30109482177700E-1_wp, 1.90391275024117E-1_wp, &
      & 1.90392904370205E-1_wp, 1.90395221896776E-1_wp, 1.90398527506004E-1_wp, &
      & 3.33798778194706E-1_wp, 3.33779134873342E-1_wp, 3.33806406791101E-1_wp, &
      & 3.33823483957318E-1_wp, -8.35501801552089E-1_wp, -8.35481140201114E-1_wp, &
      &-8.35408885793269E-1_wp, -8.35455735270754E-1_wp, -6.31069087329033E-1_wp, &
      &-6.31049654240399E-1_wp, -6.31075181828650E-1_wp, -6.31071707101107E-1_wp]

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


!> Test 1D periodic EEQ across all mask orientations, nonperiodic columns, and charge check
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
   real(wp), parameter :: weights(3) = [0.5_wp, -0.8_wp, 0.3_wp]

   type(error_type), allocatable :: expected
   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: c1, c2
   real(wp) :: lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3), dqdr_ref(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3), h1(9, 9), h2(9, 9)
   integer :: imask

   do imask = 1, 3
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

      call get_charges(model, mol, error, q1, dqdr1, dqdL1, cache=c1)
      if (allocated(error)) return

      call get_charges(model, mol_big, error, q2, dqdr2, dqdL2, cache=c2)
      if (allocated(error)) return

      h1 = 0.0_wp
      call add_charge_hessian(model, mol, error, weights, h1, c1)
      if (allocated(error)) return

      h2 = 0.0_wp
      call add_charge_hessian(model, mol_big, error, weights, h2, c2)
      if (allocated(error)) return

      if (abs(sum(q1)) > thr) then
         call test_failed(error, "1D EEQ total charge constraint violated")
         return
      end if
      if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
         & .or. maxval(abs(dqdL1 - dqdL2)) > thr &
         & .or. maxval(abs(h1 - h2)) > 1.0e-11_wp) then
         call test_failed(error, "1D EEQ depends on nonperiodic lattice vectors")
         return
      end if
      if (imask == 1) then
         q_ref = q1
         dqdr_ref = dqdr1
      else if (maxval(abs(q1 - q_ref)) > thr &
         & .or. maxval(abs(dqdr1 - dqdr_ref)) > thr) then
         call test_failed(error, "1D EEQ results depend on periodic mask index")
         return
      end if

      mol_big%charge = 1.0_wp
      if (allocated(expected)) deallocate(expected)
      call get_charges(model, mol_big, expected, q2)
      if (.not. allocated(expected)) then
         call test_failed(error, "Charged 1D EEQ cell was not rejected")
         return
      end if
   end do

   call check_all_derivs(error, mol, model)

end subroutine test_eeq_1d


!> Test 2D periodic EEQ across all mask orientations, nonperiodic column, and charge check
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
   real(wp), parameter :: weights(3) = [0.5_wp, -0.8_wp, 0.3_wp]

   type(error_type), allocatable :: expected
   type(structure_type) :: mol, mol_big
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: c1, c2
   real(wp) :: lat1(3, 3), lat2(3, 3)
   real(wp) :: q1(3), q2(3), q_ref(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3), dqdr_ref(3, 3, 3)
   real(wp) :: dqdL1(3, 3, 3), dqdL2(3, 3, 3), h1(9, 9), h2(9, 9)
   integer :: imask

   do imask = 1, 3
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

      call get_charges(model, mol, error, q1, dqdr1, dqdL1, cache=c1)
      if (allocated(error)) return

      call get_charges(model, mol_big, error, q2, dqdr2, dqdL2, cache=c2)
      if (allocated(error)) return

      h1 = 0.0_wp
      call add_charge_hessian(model, mol, error, weights, h1, c1)
      if (allocated(error)) return

      h2 = 0.0_wp
      call add_charge_hessian(model, mol_big, error, weights, h2, c2)
      if (allocated(error)) return

      if (abs(sum(q1)) > thr) then
         call test_failed(error, "2D EEQ total charge constraint violated")
         return
      end if
      if (maxval(abs(q1 - q2)) > thr .or. maxval(abs(dqdr1 - dqdr2)) > thr &
         & .or. maxval(abs(dqdL1 - dqdL2)) > thr &
         & .or. maxval(abs(h1 - h2)) > 1.0e-11_wp) then
         call test_failed(error, "2D EEQ depends on nonperiodic lattice vector")
         return
      end if
      if (imask == 1) then
         q_ref = q1
         dqdr_ref = dqdr1
      else if (maxval(abs(q1 - q_ref)) > thr &
         & .or. maxval(abs(dqdr1 - dqdr_ref)) > thr) then
         call test_failed(error, "2D EEQ results depend on periodic mask index")
         return
      end if

      mol_big%charge = -1.0_wp
      if (allocated(expected)) deallocate(expected)
      call get_charges(model, mol_big, expected, q2)
      if (.not. allocated(expected)) then
         call test_failed(error, "Charged 2D EEQ cell was not rejected")
         return
      end if
   end do

   call check_all_derivs(error, mol, model)

end subroutine test_eeq_2d


!> Test 1D periodic EEQBC across all mask orientations (neutral and charged)
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
         else if (maxval(abs(q1 - q_ref)) > thr &
            & .or. maxval(abs(dqdr1 - dqdr_ref)) > thr) then
            call test_failed(error, "1D EEQBC results depend on periodic mask index")
            return
         end if
      end do

      call check_all_derivs(error, mol, model)
      if (allocated(error)) return
   end do

end subroutine test_eeqbc_1d


!> Test 2D periodic EEQBC across all mask orientations (neutral and charged)
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
         else if (maxval(abs(q1 - q_ref)) > thr &
            & .or. maxval(abs(dqdr1 - dqdr_ref)) > thr) then
            call test_failed(error, "2D EEQBC results depend on periodic mask index")
            return
         end if
      end do

      call check_all_derivs(error, mol, model)
      if (allocated(error)) return
   end do

end subroutine test_eeqbc_2d


!> Verify closed-form alpha, subspace translation filtering, projected rmax, and error handling
subroutine test_ewald_setup(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: lat_ref(3, 3) = reshape([ &
      & 4.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 5.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 8.0_wp], [3, 3])
   integer, parameter :: num_2(2) = [5, 7]

   type(error_type), allocatable :: expected
   type(structure_type) :: mol_near, mol_far
   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: c_near, c_far
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)
   real(wp) :: lat_deg(3, 3), lat_empty(0, 0), xyz_2(3, 2), q(2)
   real(wp) :: a_def, a_1, a_3, a_1d, a_2d
   integer :: ndir_near, nrec_near, ndir_far, nrec_far

   ! Omitted mask, 1-element [.true.], and 3-element [.true., .true., .true.] agree
   call get_alpha(lat_ref, a_def, error)
   if (allocated(error)) return

   call get_alpha(lat_ref, a_1, error, periodic=[.true.])
   if (allocated(error)) return

   call get_alpha(lat_ref, a_3, error, periodic=[.true., .true., .true.])
   if (allocated(error)) return

   if (abs(a_def - a_3) > thr .or. abs(a_1 - a_3) > thr &
      & .or. abs(a_3 - 1.25_wp * sqrt(pi) / 160.0_wp**(1.0_wp / 3.0_wp)) > thr) then
      call test_failed(error, "3D closed-form alpha mismatch")
      return
   end if

   ! Closed-form 1D and 2D alpha and subspace-restricted translations
   call get_alpha(lat_ref, a_1d, error, periodic=[.false., .true., .false.])
   if (allocated(error)) return

   call get_alpha(lat_ref, a_2d, error, periodic=[.true., .false., .true.])
   if (allocated(error)) return

   if (abs(a_1d - 0.5_wp * sqrt(pi) / 5.0_wp) > thr &
      & .or. abs(a_2d - 1.0_wp * sqrt(pi) / sqrt(32.0_wp)) > thr) then
      call test_failed(error, "1D or 2D closed-form alpha mismatch")
      return
   end if

   call get_ewald_trans(lat_ref, a_2d, 0.25_wp, 2.0_wp, dtrans, rtrans, error, &
      & periodic=[.true., .false., .true.])
   if (allocated(error)) return

   if (maxval(abs(dtrans(2, :))) > 0.0_wp .or. maxval(abs(rtrans(2, :))) > 0.0_wp) then
      call test_failed(error, "Inactive axis has nonzero translations")
      return
   end if

   ! Projected rmax: large perpendicular separation must not expand 2D dtrans
   xyz_2 = reshape([0.0_wp, 0.0_wp, 0.0_wp, 1.4_wp, 0.8_wp, 1.5_wp], [3, 2])
   call new(mol_near, num_2, xyz_2, lattice=lat_ref, periodic=[.true., .true., .false.])
   xyz_2(3, 2) = 80.0_wp
   call new(mol_far, num_2, xyz_2, lattice=lat_ref, periodic=[.true., .true., .false.])
   call new_eeq2019_model(mol_near, model, error)
   if (allocated(error)) return

   call get_charges(model, mol_near, error, q, cache=c_near)
   if (allocated(error)) return

   call get_charges(model, mol_far, error, q, cache=c_far)
   if (allocated(error)) return

   select type (raw => c_near%raw)
   type is (eeq_cache)
      ndir_near = size(raw%dtrans, 2)
      nrec_near = size(raw%rtrans, 2)
   end select
   select type (raw => c_far%raw)
   type is (eeq_cache)
      ndir_far = size(raw%dtrans, 2)
      nrec_far = size(raw%rtrans, 2)
   end select
   if (ndir_near /= ndir_far .or. nrec_near /= nrec_far) then
      call test_failed(error, "Perpendicular pair separation inflated 2D Ewald translations")
      return
   end if

   ! Degenerate zero lattice
   lat_deg = 0.0_wp
   call get_alpha(lat_deg, a_1, expected, periodic=[.true., .true., .true.])
   if (.not. allocated(expected)) then
      call test_failed(error, "get_alpha did not report error for zero lattice")
      return
   end if

   ! Linearly dependent (collinear) 2D active lattice vectors
   deallocate(expected)
   lat_deg = 0.0_wp
   lat_deg(:, 1) = [4.0_wp, 1.0_wp, 0.0_wp]
   lat_deg(:, 3) = [8.0_wp, 2.0_wp, 0.0_wp]
   call get_alpha(lat_deg, a_1, expected, periodic=[.true., .false., .true.])
   if (.not. allocated(expected)) then
      call test_failed(error, "get_alpha did not reject collinear 2D lattice vectors")
      return
   end if

   deallocate(expected)
   call get_ewald_trans(lat_deg, 0.25_wp, 0.2_wp, 1.0_wp, dtrans, rtrans, expected, &
      & periodic=[.true., .false., .true.])
   if (.not. allocated(expected)) then
      call test_failed(error, "get_ewald_trans did not reject collinear 2D lattice vectors")
      return
   end if

   ! Nonperiodic mask
   deallocate(expected)
   call get_alpha(lat_ref, a_1, expected, periodic=[.false., .false., .false.])
   if (.not. allocated(expected)) then
      call test_failed(error, "get_alpha did not reject nonperiodic mask")
      return
   end if

   ! Undersized lattice array
   deallocate(expected)
   call get_ewald_trans(lat_empty, 0.25_wp, 0.2_wp, 0.0_wp, dtrans, rtrans, expected, &
      & periodic=[.true., .true., .true.])
   if (.not. allocated(expected)) then
      call test_failed(error, "get_ewald_trans did not reject empty lattice")
      return
   end if

   ! Invalid cutoff parameters and reciprocal-box integer overflow
   deallocate(expected)
   call get_ewald_trans(lat_ref, -0.2_wp, 0.2_wp, 0.0_wp, dtrans, rtrans, expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "get_ewald_trans did not reject nonpositive alpha")
      return
   end if

   deallocate(expected)
   call get_ewald_trans(lat_ref, 1.0e8_wp, 0.2_wp, 0.0_wp, dtrans, rtrans, expected)
   if (.not. allocated(expected)) then
      call test_failed(error, "get_ewald_trans did not reject reciprocal overflow")
      return
   end if

end subroutine test_ewald_setup


!> Verify Ewald splitting parameter invariance for 1D periodic EEQ
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


!> Verify Ewald splitting parameter invariance for 2D periodic EEQ
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

   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(model%rad))

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

   do ia = 1, size(alphas)
      select type (raw => cache%raw)
      type is (eeq_cache)
         call get_ewald_setup(mol, raw%wsc, gam_min, raw%ndim, raw%measure, &
            & raw%axis_vec, raw%alpha, raw%dtrans, raw%rtrans, error, &
            & alpha_in=alphas(ia))
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


!> Verify that potential-based cutoffs at machine epsilon converge 1D and 2D
!> charges, energies, gradients, strain derivatives, and Hessians against enlarged cutoffs
subroutine test_ewald_cutoff_convergence(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz_wire(3, 4) = reshape([ &
      &  0.00_wp,  0.00_wp,  0.00_wp, &
      &  0.35_wp,  1.95_wp, -0.25_wp, &
      & -0.60_wp,  3.10_wp,  1.45_wp, &
      &  0.80_wp, -0.45_wp, -1.30_wp], [3, 4])
   real(wp), parameter :: lat_wire(3, 3) = reshape([ &
      & 0.0_wp,  0.0_wp,  0.0_wp, &
      & 0.5_wp,  5.6_wp, -0.4_wp, &
      & 0.0_wp,  0.0_wp,  0.0_wp], [3, 3])
   real(wp), parameter :: xyz_slab(3, 4) = reshape([ &
      & 0.00_wp, -0.40_wp, 0.00_wp, &
      & 2.35_wp,  0.65_wp, 0.30_wp, &
      & 0.85_wp,  2.10_wp, 2.25_wp, &
      & 3.10_wp, -1.55_wp, 3.80_wp], [3, 4])
   real(wp), parameter :: lat_slab(3, 3) = reshape([ &
      & 5.10_wp,  0.20_wp, 0.35_wp, &
      & 0.00_wp,  0.00_wp, 0.00_wp, &
      & 1.40_wp, -0.30_wp, 6.40_wp], [3, 3])

   type(structure_type) :: mol

   call new(mol, [30, 8, 6, 1], xyz_wire, lattice=lat_wire, &
      & periodic=[.false., .true., .false.])
   call check_enlarged_cutoff(error, mol, "1D wire")
   if (allocated(error)) return

   call new(mol, [5, 7, 5, 7], xyz_slab, lattice=lat_slab, &
      & periodic=[.true., .false., .true.])
   call check_enlarged_cutoff(error, mol, "2D buckled slab")

end subroutine test_ewald_cutoff_convergence


!> Compare default machine-epsilon Ewald cutoffs against enlarged (1e-4 * eps) cutoffs
subroutine check_enlarged_cutoff(error, mol, label)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Periodic 4-atom structure
   type(structure_type), intent(in) :: mol

   !> Diagnostic label
   character(len=*), intent(in) :: label

   real(wp), parameter :: qvec(5) = [0.42_wp, -0.61_wp, 0.33_wp, -0.14_wp, 0.0_wp]
   real(wp), parameter :: adj(5) = [-0.30_wp, 0.50_wp, -0.45_wp, 0.25_wp, 0.0_wp]

   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: cache
   real(wp), allocatable :: trans(:, :)
   real(wp) :: cn(4), qloc(4), dcndr(3, 4, 4), dcndL(3, 3, 4), gam_min, alpha_def
   real(wp) :: amat_def(5, 5), amat_big(5, 5)
   real(wp) :: dadr_def(3, 4, 5), dadL_def(3, 3, 5), atrace_def(3, 4), hess_def(12, 12)
   real(wp) :: dadr_big(3, 4, 5), dadL_big(3, 3, 5), atrace_big(3, 4), hess_big(12, 12)
   integer :: ndir_def, nrec_def, ndir_big, nrec_big

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   gam_min = 1.0_wp / (sqrt(2.0_wp) * maxval(model%rad))

   call get_lattice_points(mol%periodic, mol%lattice, model%ncoord%cutoff, trans, error)
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)
   call model%local_charge(mol, trans, qloc)
   call model%update(mol, error, cache, cn, qloc, dcndr, dcndL)
   if (allocated(error)) return

   allocate(cache%qvec(5), source=qvec)
   call model%get_coulomb_matrix(mol, cache, amat_def)
   call model%get_coulomb_derivs(mol, cache, qvec, dadr_def, dadL_def, atrace_def)
   hess_def = 0.0_wp
   call model%add_model_hessian(mol, error, cache, adj, hess_def)
   if (allocated(error)) return

   select type (raw => cache%raw)
   type is (eeq_cache)
      alpha_def = raw%alpha
      ndir_def = size(raw%dtrans, 2)
      nrec_def = size(raw%rtrans, 2)
      call get_ewald_setup(mol, raw%wsc, gam_min, raw%ndim, raw%measure, &
         & raw%axis_vec, raw%alpha, raw%dtrans, raw%rtrans, error, &
         & alpha_in=alpha_def, tol_scale=1.0e-4_wp)
      if (allocated(error)) return
      ndir_big = size(raw%dtrans, 2)
      nrec_big = size(raw%rtrans, 2)
   end select

   if (ndir_big <= ndir_def .and. nrec_big <= nrec_def) then
      call test_failed(error, label // ": enlarged cutoff did not add lattice vectors")
      return
   end if

   call model%get_coulomb_matrix(mol, cache, amat_big)
   call model%get_coulomb_derivs(mol, cache, qvec, dadr_big, dadL_big, atrace_big)
   hess_big = 0.0_wp
   call model%add_model_hessian(mol, error, cache, adj, hess_big)
   if (allocated(error)) return

   if (maxval(abs(amat_def(:4, :4) - amat_big(:4, :4))) > 1.0e-13_wp &
      & .or. maxval(abs(dadr_def - dadr_big)) > 1.0e-13_wp &
      & .or. maxval(abs(atrace_def - atrace_big)) > 1.0e-13_wp &
      & .or. maxval(abs(dadL_def - dadL_big)) > 1.0e-12_wp &
      & .or. maxval(abs(hess_def - hess_big)) > 1.0e-12_wp) then
      call test_failed(error, label // ": potential-based Ewald cutoff is not converged")
   end if

end subroutine check_enlarged_cutoff


!> Verify invariance under unimodular basis transforms and multi-cell shifts in 1D, 2D, and 3D
subroutine test_ewald_lattice_invariance(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: num(3) = [8, 1, 1]
   real(wp), parameter :: xyz_ref(3, 3) = reshape([ &
      & 0.1_wp, 0.2_wp, -0.1_wp, &
      & 1.7_wp, 0.6_wp,  0.3_wp, &
      & 0.3_wp, 1.8_wp,  0.5_wp], [3, 3])
   real(wp), parameter :: a1_2d(3) = [5.5_wp, 0.0_wp, 0.0_wp]
   real(wp), parameter :: a2_2d(3) = [0.5_wp, 6.0_wp, 2.0_wp]

   type(structure_type) :: mol, mol_shift
   real(wp) :: lat_ref(3, 3), lat_skew(3, 3), xyz_shift(3, 3)

   ! 1D: axis inversion (unimodular U = [-1]) + multi-cell shift
   lat_ref = 0.0_wp
   lat_ref(:, 2) = [0.6_wp, 5.4_wp, -0.4_wp]
   lat_skew = 0.0_wp
   lat_skew(:, 2) = -lat_ref(:, 2)
   xyz_shift = xyz_ref
   xyz_shift(:, 1) = xyz_ref(:, 1) + 2.0_wp * lat_ref(:, 2)
   xyz_shift(:, 3) = xyz_ref(:, 3) - 3.0_wp * lat_ref(:, 2)
   call new(mol, num, xyz_ref, lattice=lat_ref, periodic=[.false., .true., .false.])
   call new(mol_shift, num, xyz_shift, lattice=lat_skew, periodic=[.false., .true., .false.])
   call check_lattice_equiv(error, mol, mol_shift, "1D")
   if (allocated(error)) return

   ! 2D: unimodular shear (a2' = a2 + 2*a1) + multi-cell shift
   lat_ref = 0.0_wp
   lat_ref(:, 1) = a1_2d
   lat_ref(:, 3) = a2_2d
   lat_skew = 0.0_wp
   lat_skew(:, 1) = a1_2d
   lat_skew(:, 3) = a2_2d + 2.0_wp * a1_2d
   xyz_shift = xyz_ref
   xyz_shift(:, 2) = xyz_ref(:, 2) + a1_2d
   xyz_shift(:, 3) = xyz_ref(:, 3) - a2_2d
   call new(mol, num, xyz_ref, lattice=lat_ref, periodic=[.true., .false., .true.])
   call new(mol_shift, num, xyz_shift, lattice=lat_skew, periodic=[.true., .false., .true.])
   call check_lattice_equiv(error, mol, mol_shift, "2D")
   if (allocated(error)) return

   ! 3D: unimodular SL(3, Z) shear (a1' = a1 + a2, a3' = a3 - a1, det U = 1) + shifts
   lat_ref = znoh_lat
   lat_skew(:, 1) = znoh_lat(:, 1) + znoh_lat(:, 2)
   lat_skew(:, 2) = znoh_lat(:, 2)
   lat_skew(:, 3) = znoh_lat(:, 3) - znoh_lat(:, 1)
   xyz_shift = znoh_xyz
   xyz_shift(:, 1) = znoh_xyz(:, 1) + znoh_lat(:, 1) - znoh_lat(:, 3)
   xyz_shift(:, 2) = znoh_xyz(:, 2) - 2.0_wp * znoh_lat(:, 2)
   xyz_shift(:, 3) = znoh_xyz(:, 3) + znoh_lat(:, 3)
   call new(mol, znoh_num, znoh_xyz, lattice=lat_ref, periodic=[.true., .true., .true.])
   call new(mol_shift, znoh_num, xyz_shift, lattice=lat_skew, periodic=[.true., .true., .true.])
   call check_lattice_equiv(error, mol, mol_shift, "3D")

end subroutine test_ewald_lattice_invariance


!> Compare energy, gradient, stress, charges, response, and Hessian of two equivalent cells
subroutine check_lattice_equiv(error, mol1, mol2, label)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Reference representation
   type(structure_type), intent(in) :: mol1

   !> Sheared/shifted equivalent representation
   type(structure_type), intent(in) :: mol2

   !> Dimensionality label
   character(len=*), intent(in) :: label

   real(wp), parameter :: weights(3) = [0.6_wp, -0.4_wp, 0.25_wp]

   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: c1, c2
   real(wp), allocatable :: tr1(:, :), tr2(:, :)
   real(wp) :: cn1(3), qloc1(3), dcndr1(3, 3, 3), dcndL1(3, 3, 3)
   real(wp) :: cn2(3), qloc2(3), dcndr2(3, 3, 3), dcndL2(3, 3, 3)
   real(wp) :: dqlr(3, 3, 3), dqlL(3, 3, 3)
   real(wp) :: e1(3), e2(3), g1(3, 3), g2(3, 3), s1(3, 3), s2(3, 3)
   real(wp) :: q1(3), q2(3), dqdr1(3, 3, 3), dqdr2(3, 3, 3), dqdL1(3, 3, 3), dqdL2(3, 3, 3)
   real(wp) :: h1(9, 9), h2(9, 9)

   call new_eeq2019_model(mol1, model, error)
   if (allocated(error)) return

   call get_lattice_points(mol1%periodic, mol1%lattice, model%ncoord%cutoff, tr1, error)
   if (allocated(error)) return

   call get_lattice_points(mol2%periodic, mol2%lattice, model%ncoord%cutoff, tr2, error)
   if (allocated(error)) return

   e1 = 0.0_wp
   g1 = 0.0_wp
   s1 = 0.0_wp
   call model%ncoord%get_coordination_number(mol1, tr1, cn1, dcndr1, dcndL1)
   call model%local_charge(mol1, tr1, qloc1, dqlr, dqlL)
   call model%solve(mol1, error, cn1, qloc1, dcndr1, dcndL1, dqlr, dqlL, &
      & energy=e1, gradient=g1, sigma=s1, qvec=q1, dqdr=dqdr1, dqdL=dqdL1, response=c1)
   if (allocated(error)) return

   e2 = 0.0_wp
   g2 = 0.0_wp
   s2 = 0.0_wp
   call model%ncoord%get_coordination_number(mol2, tr2, cn2, dcndr2, dcndL2)
   call model%local_charge(mol2, tr2, qloc2, dqlr, dqlL)
   call model%solve(mol2, error, cn2, qloc2, dcndr2, dcndL2, dqlr, dqlL, &
      & energy=e2, gradient=g2, sigma=s2, qvec=q2, dqdr=dqdr2, dqdL=dqdL2, response=c2)
   if (allocated(error)) return

   h1 = 0.0_wp
   call add_charge_hessian(model, mol1, error, weights, h1, c1)
   if (allocated(error)) return

   h2 = 0.0_wp
   call add_charge_hessian(model, mol2, error, weights, h2, c2)
   if (allocated(error)) return

   if (maxval(abs(e1 - e2)) > 1.0e-11_wp &
      & .or. maxval(abs(q1 - q2)) > 1.0e-11_wp &
      & .or. maxval(abs(g1 - g2)) > 1.0e-11_wp &
      & .or. maxval(abs(s1 - s2)) > 1.0e-10_wp &
      & .or. maxval(abs(dqdr1 - dqdr2)) > 1.0e-10_wp &
      & .or. maxval(abs(dqdL1 - dqdL2)) > 1.0e-10_wp &
      & .or. maxval(abs(h1 - h2)) > 1.0e-9_wp) then
      call test_failed(error, label // " EEQ not invariant under unimodular transform and shift")
   end if

end subroutine check_lattice_equiv


!> Validate analytical gradients, strain derivatives, and response on strongly
!> anisotropic and skewed cells in 1D, 2D, and 3D
subroutine test_anisotropic_skewed(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz_1d(3, 3) = reshape([ &
      & 0.0_wp,  0.0_wp,  0.0_wp, &
      & 3.8_wp, -2.5_wp,  1.2_wp, &
      & -2.1_wp, 3.4_wp,  2.1_wp], [3, 3])
   real(wp), parameter :: lat_1d(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp,  0.0_wp, &
      & 0.0_wp, 0.0_wp,  0.0_wp, &
      & 1.4_wp, -0.9_wp, 4.6_wp], [3, 3])
   real(wp), parameter :: xyz_2d(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp, -0.5_wp, &
      & 2.1_wp, 5.2_wp,  0.9_wp, &
      & 4.3_wp, 9.8_wp, -1.1_wp], [3, 3])
   real(wp), parameter :: lat_2d(3, 3) = reshape([ &
      & 4.5_wp,  0.2_wp,  0.1_wp, &
      & 6.8_wp, 17.5_wp, -0.6_wp, &
      & 0.0_wp,  0.0_wp,  0.0_wp], [3, 3])
   real(wp), parameter :: lat_3d(3, 3) = reshape([ &
      & 4.6_wp,  0.1_wp,  0.0_wp, &
      & 3.9_wp,  5.1_wp,  0.2_wp, &
      & 0.4_wp, -0.5_wp, 24.0_wp], [3, 3])

   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model

   call new(mol, znoh_num, xyz_1d, lattice=lat_1d, periodic=[.false., .false., .true.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   call test_numgrad(error, mol, model)
   if (allocated(error)) return

   call test_numsigma(error, mol, model)
   if (allocated(error)) return

   call new(mol, znoh_num, xyz_2d, lattice=lat_2d, periodic=[.true., .true., .false.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   call test_numgrad(error, mol, model)
   if (allocated(error)) return

   call test_numsigma(error, mol, model)
   if (allocated(error)) return

   call test_numdqdL(error, mol, model)
   if (allocated(error)) return

   call new(mol, znoh_num, znoh_xyz, lattice=lat_3d, periodic=[.true., .true., .true.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   call test_numgrad(error, mol, model)
   if (allocated(error)) return

   call test_numsigma(error, mol, model)

end subroutine test_anisotropic_skewed


!> Verify that cached calculations match fresh calculations and that repeated
!> geometry/dimensionality/model updates invalidate and regenerate the cache
subroutine test_cache_invalidation(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: lat_1d(3, 3) = reshape([ &
      & 5.3_wp, 0.6_wp, -0.5_wp, &
      & 0.0_wp, 0.0_wp,  0.0_wp, &
      & 0.0_wp, 0.0_wp,  0.0_wp], [3, 3])
   real(wp), parameter :: lat_2d(3, 3) = reshape([ &
      & 5.2_wp, 0.3_wp,  0.1_wp, &
      & 0.0_wp, 0.0_wp,  0.0_wp, &
      & 0.5_wp, 6.1_wp, -0.2_wp], [3, 3])

   type(structure_type) :: mol
   type(cache_container) :: shared_cache

   ! Step 1: 3D periodic EEQ
   call new(mol, znoh_num, znoh_xyz, lattice=znoh_lat, periodic=[.true., .true., .true.])
   call check_shared_vs_fresh(error, mol, .false., shared_cache)
   if (allocated(error)) return

   ! Step 2: 1D periodic EEQ on the same cache
   call new(mol, znoh_num, znoh_xyz, lattice=lat_1d, periodic=[.true., .false., .false.])
   call check_shared_vs_fresh(error, mol, .false., shared_cache)
   if (allocated(error)) return

   ! Step 3: 2D periodic EEQ with perturbed geometry on the same cache
   call new(mol, znoh_num, 1.08_wp * znoh_xyz, lattice=lat_2d, &
      & periodic=[.true., .false., .true.])
   call check_shared_vs_fresh(error, mol, .false., shared_cache)
   if (allocated(error)) return

   ! Step 4: 0D molecular EEQ with different atom count (2 atoms) on the same cache
   call new(mol, [8, 1], znoh_xyz(:, :2))
   call check_shared_vs_fresh(error, mol, .false., shared_cache)
   if (allocated(error)) return

   ! Step 5: Switch model type to charged 1D EEQBC on the same cache
   call new(mol, znoh_num, znoh_xyz, charge=1.0_wp, lattice=lat_1d, &
      & periodic=[.true., .false., .false.])
   call check_shared_vs_fresh(error, mol, .true., shared_cache)
   if (allocated(error)) return

   ! Step 6: Switch back to 2D EEQ with scaled lattice on the same cache
   call new(mol, znoh_num, znoh_xyz, lattice=1.25_wp * lat_2d, &
      & periodic=[.true., .false., .true.])
   call check_shared_vs_fresh(error, mol, .false., shared_cache)

end subroutine test_cache_invalidation


!> Compare a calculation reusing shared_cache against a freshly initialized cache
subroutine check_shared_vs_fresh(error, mol, use_eeqbc, shared_cache)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure
   type(structure_type), intent(in) :: mol

   !> Select EEQBC (.true.) or EEQ (.false.)
   logical, intent(in) :: use_eeqbc

   !> Persistent cache reused across geometry/model updates
   type(cache_container), intent(inout) :: shared_cache

   class(mchrg_model_type), allocatable :: model
   type(cache_container) :: fresh_cache
   real(wp), allocatable :: q_s(:), q_f(:), dqdr_s(:, :, :), dqdr_f(:, :, :)
   real(wp), allocatable :: dqdL_s(:, :, :), dqdL_f(:, :, :), w(:), h_s(:, :), h_f(:, :)
   integer :: n, i

   n = mol%nat
   allocate(q_s(n), q_f(n), dqdr_s(3, n, n), dqdr_f(3, n, n))
   allocate(dqdL_s(3, 3, n), dqdL_f(3, 3, n), w(n), h_s(3 * n, 3 * n), h_f(3 * n, 3 * n))
   w = [(sin(1.9_wp * real(i, wp)), i = 1, n)]

   if (use_eeqbc) then
      call new_eeqbc2025_model(mol, model, error)
   else
      call new_eeq2019_model(mol, model, error)
   end if
   if (allocated(error)) return

   call get_charges(model, mol, error, q_s, dqdr_s, dqdL_s, cache=shared_cache)
   if (allocated(error)) return

   call get_charges(model, mol, error, q_f, dqdr_f, dqdL_f, cache=fresh_cache)
   if (allocated(error)) return

   h_s = 0.0_wp
   call add_charge_hessian(model, mol, error, w, h_s, shared_cache)
   if (allocated(error)) return

   h_f = 0.0_wp
   call add_charge_hessian(model, mol, error, w, h_f, fresh_cache)
   if (allocated(error)) return

   if (maxval(abs(q_s - q_f)) > 1.0e-13_wp &
      & .or. maxval(abs(dqdr_s - dqdr_f)) > 1.0e-13_wp &
      & .or. maxval(abs(dqdL_s - dqdL_f)) > 1.0e-13_wp &
      & .or. maxval(abs(h_s - h_f)) > 1.0e-12_wp) then
      call test_failed(error, "Reused cache disagrees with freshly initialized cache")
   end if

end subroutine check_shared_vs_fresh


!> Regression test for DFT-D4 issue #285: 2D slab reference charges across cell shapes
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


!> Verify 3D large-vacuum slab convergence (1/Lz dipole-sheet scaling) toward 2D EEQ
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


!> Verify 1D and 2D Ewald special-function accuracy against 50-digit references
subroutine test_special_functions(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: e1_pts(7) = [ &
      & 1.0e-8_wp, 0.25_wp, 1.0_wp, 1.000001_wp, 4.0_wp, 4.000001_wp, 15.0_wp]
   real(wp), parameter :: e1_ref(7) = [ &
      & 17.843465089050832587_wp, 1.0442826344437381945_wp, &
      & 0.21938393439552027368_wp, 0.21938356651644698137_wp, &
      & 3.7793524098489064789e-3_wp, 3.7793478309420461127e-3_wp, &
      & 1.9186278921478669771e-8_wp]

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

   real(wp) :: e1, f0, f1, f2, k0, k1, k2, v_2d
   integer :: i

   do i = 1, size(e1_pts)
      e1 = expint_e1(e1_pts(i))
      if (abs(e1 - e1_ref(i)) > 1.0e-14_wp * e1_ref(i)) then
         call test_failed(error, "expint_e1 deviates from 50-digit reference")
         return
      end if
   end do

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


!> Verify linear energy scaling and charge invariance for 1D elongated supercells
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


!> Verify linear energy scaling and charge invariance for 2D elongated supercells
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


!> Verify linear energy scaling and charge invariance for 3D elongated supercells
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
