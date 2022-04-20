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
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_data, only : get_covalent_rad
   use multicharge_model
   use multicharge_ncoord, only : get_coordination_number, cut_coordination_number
   use multicharge_output, only : write_ascii_model, write_ascii_properties, &
      & write_ascii_results
   use multicharge_param, only : new_eeq2019_model
   implicit none
   private

   public :: collect_pbc

   real(wp), parameter :: thr = 1000*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_pbc(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("charges-cyanamide", test_q_cyanamide), &
      & new_unittest("energy-formamide", test_e_formamide), &
      & new_unittest("gradient-co2", test_g_co2), &
      & new_unittest("sigma-ice", test_s_ice), &
      & new_unittest("dqdr-urea", test_dqdr_urea), &
      & new_unittest("dqdL-oxacb", test_dqdL_oxacb) &
      & ]

end subroutine collect_pbc


subroutine gen_test(error, mol, qref, eref)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Reference charges
   real(wp), intent(in), optional :: qref(:)

   !> Reference energies
   real(wp), intent(in), optional :: eref(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(mchrg_model_type) :: model
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), rcov(:), trans(:, :)
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)

   call new_eeq2019_model(mol, model)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   allocate(cn(mol%nat))

   rcov = get_covalent_rad(mol%num)
   call get_coordination_number(mol, trans, cutoff, rcov, cn, cut=cn_max)

   if (present(eref)) then
      allocate(energy(mol%nat))
      energy(:) = 0.0_wp
   end if
   if (present(qref)) then
      allocate(qvec(mol%nat))
   end if

   call model%solve(mol, cn, energy=energy, qvec=qvec)
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


subroutine test_numgrad(error, mol)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   type(mchrg_model_type) :: model
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:), trans(:, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :)
   real(wp) :: er, el

   call new_eeq2019_model(mol, model)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   rcov = get_covalent_rad(mol%num)
   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), sigma(3, 3), numgrad(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         energy(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_coordination_number(mol, trans, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_coordination_number(mol, trans, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do lp
   if (allocated(error)) return

   call get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut=cn_max)

   energy(:) = 0.0_wp
   call model%solve(mol, cn, dcndr, dcndL, energy, gradient, sigma)
   if (allocated(error)) return

   if (any(abs(gradient(:, :) - numgrad(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, jc
   type(mchrg_model_type) :: model
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:), trans(:, :)
   real(wp), allocatable :: energy(:), gradient(:, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :)
   real(wp) :: er, el, eps(3, 3), numsigma(3, 3), sigma(3, 3), lattice(3, 3)

   call new_eeq2019_model(mol, model)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   rcov = get_covalent_rad(mol%num)
   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
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
         call get_coordination_number(mol, lattr, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do lp
   if (allocated(error)) return

   call get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut=cn_max)

   energy(:) = 0.0_wp
   call model%solve(mol, cn, dcndr, dcndL, energy, gradient, sigma)
   if (allocated(error)) return

   if (any(abs(sigma(:, :) - numsigma(:, :)) > thr2)) then
      call test_failed(error, "Derivative of energy does not match")
   end if

end subroutine test_numsigma


subroutine test_numdqdr(error, mol)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: iat, ic
   type(mchrg_model_type) :: model
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:), trans(:, :)
   real(wp), allocatable :: ql(:), qr(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)

   call new_eeq2019_model(mol, model)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   rcov = get_covalent_rad(mol%num)
   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & ql(mol%nat), qr(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_coordination_number(mol, trans, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, qvec=qr)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_coordination_number(mol, trans, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, qvec=ql)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp*(qr - ql)/step
      end do
   end do lp
   if (allocated(error)) return

   call get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut=cn_max)

   call model%solve(mol, cn, dcndr, dcndL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdr(:, :, :) - numdr(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of charges does not match")
   end if

end subroutine test_numdqdr


subroutine test_numdqdL(error, mol)

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer :: ic, jc
   type(mchrg_model_type) :: model
   real(wp), parameter :: cn_max = 8.0_wp, cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:), trans(:, :)
   real(wp), allocatable :: qr(:), ql(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :), numdL(:, :, :)
   real(wp) :: eps(3, 3), lattice(3, 3)

   call new_eeq2019_model(mol, model)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, trans)

   rcov = get_covalent_rad(mol%num)
   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & qr(mol%nat), ql(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & xyz(3, mol%nat), numdL(3, 3, mol%nat))

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   lattice(:, :) = mol%lattice
   lattr = trans
   lp: do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, qvec=qr)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         mol%lattice(:, :) = matmul(eps, lattice)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, qvec=ql)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         mol%lattice(:, :) = lattice
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp*(qr - ql)/step
      end do
   end do lp
   if (allocated(error)) return

   call get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut=cn_max)

   call model%solve(mol, cn, dcndr, dcndL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdL(:, :, :) - numdL(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of charges does not match")
   end if

end subroutine test_numdqdL


subroutine test_q_cyanamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(40) = [&
      & 3.47274007973765E-1_wp, 3.47273017259661E-1_wp, 3.47276806231462E-1_wp, &
      & 3.47267235121650E-1_wp, 3.47288751410308E-1_wp, 3.47296840650318E-1_wp, &
      & 3.47299868243747E-1_wp, 3.47299359754750E-1_wp, 3.45885710186856E-1_wp, &
      & 3.45874246032015E-1_wp, 3.45888242047875E-1_wp, 3.45877451600398E-1_wp, &
      & 3.45902365123333E-1_wp, 3.45902162041418E-1_wp, 3.45900336974539E-1_wp, &
      & 3.45903770121699E-1_wp, 3.58996736955583E-1_wp, 3.58987991603151E-1_wp, &
      & 3.59000859720791E-1_wp, 3.58990960169615E-1_wp, 3.58990422147650E-1_wp, &
      & 3.58997842028106E-1_wp, 3.58997124263551E-1_wp, 3.59001728635141E-1_wp, &
      &-5.86637659146611E-1_wp,-5.86601835343945E-1_wp,-5.86625078062709E-1_wp, &
      &-5.86579999312624E-1_wp,-5.86654394920376E-1_wp,-5.86669360667441E-1_wp, &
      &-5.86667481028248E-1_wp,-5.86674873278027E-1_wp,-4.65529092171338E-1_wp, &
      &-4.65527713392387E-1_wp,-4.65546450764424E-1_wp,-4.65540905286785E-1_wp, &
      &-4.65529860758785E-1_wp,-4.65533757413217E-1_wp,-4.65527683275370E-1_wp, &
      &-4.65527691475100E-1_wp]

   call get_structure(mol, "X23", "cyanamide")
   call gen_test(error, mol, qref=ref)

end subroutine test_q_cyanamide


subroutine test_e_formamide(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(24) = [&
      & 1.89757620852359E-1_wp, 1.89760760029133E-1_wp, 1.89737950594800E-1_wp, &
      & 1.89748797090202E-1_wp, 2.03099343220260E-1_wp, 2.03086622963944E-1_wp, &
      & 2.03077163423409E-1_wp, 2.03088145719017E-1_wp, 8.99074841863554E-2_wp, &
      & 8.99082451986937E-2_wp, 8.99093436347654E-2_wp, 8.99109153737500E-2_wp, &
      & 1.64940329265701E-1_wp, 1.64930611448213E-1_wp, 1.64944104143581E-1_wp, &
      & 1.64952535219774E-1_wp,-4.28537080109528E-1_wp,-4.28526505049323E-1_wp, &
      &-4.28489488988212E-1_wp,-4.28513466458521E-1_wp,-3.29599028782344E-1_wp, &
      &-3.29588887030536E-1_wp,-3.29602200135643E-1_wp,-3.29600393007701E-1_wp]


   call get_structure(mol, "X23", "formamide")
   call gen_test(error, mol, eref=ref)

end subroutine test_e_formamide


subroutine test_g_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "CO2")
   call test_numgrad(error, mol)

end subroutine test_g_co2


subroutine test_s_ice(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "ICE10", "vi")
   call test_numsigma(error, mol)

end subroutine test_s_ice


subroutine test_dqdr_urea(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "urea")
   call test_numdqdr(error, mol)

end subroutine test_dqdr_urea


subroutine test_dqdL_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "oxacb")
   call test_numdqdL(error, mol)

end subroutine test_dqdL_oxacb


end module test_pbc
