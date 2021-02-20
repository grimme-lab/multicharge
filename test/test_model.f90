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
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use multicharge_data, only : get_covalent_rad
   use multicharge_model
   use multicharge_ncoord, only : get_coordination_number, cut_coordination_number
   use multicharge_output, only : write_ascii_model, write_ascii_properties, &
      & write_ascii_results
   use multicharge_param, only : new_eeq2019_model
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
      & new_unittest("charges-mb01", test_q_mb01), &
      & new_unittest("charges-mb02", test_q_mb02), &
      & new_unittest("energy-mb03", test_e_mb03), &
      & new_unittest("energy-mb04", test_e_mb04), &
      & new_unittest("gradient-mb05", test_g_mb05), &
      & new_unittest("gradient-mb06", test_g_mb06), &
      & new_unittest("sigma-mb07", test_s_mb07), &
      & new_unittest("sigma-mb08", test_s_mb08), &
      & new_unittest("dqdr-mb09", test_dqdr_mb09), &
      & new_unittest("dqdr-mb10", test_dqdr_mb10), &
      & new_unittest("dqdL-mb11", test_dqdL_mb11), &
      & new_unittest("dqdL-mb12", test_dqdL_mb12) &
      & ]

end subroutine collect_model


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
   real(wp), parameter :: cn_max = 8.0_wp, trans(3, 1) = 0.0_wp, cutoff = 25.0_wp
   real(wp), allocatable :: cn(:), rcov(:)
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)

   call new_eeq2019_model(mol, model)

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
      end if
   end if
   if (allocated(error)) return

   if (present(eref)) then
      if (any(abs(energy - eref) > thr)) then
         call test_failed(error, "Energies do not match")
         print'(a)', "Energy:"
         print'(3es21.14)', energy
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
   real(wp), parameter :: cn_max = 8.0_wp, trans(3, 1) = 0.0_wp, cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :)
   real(wp) :: er, el

   call new_eeq2019_model(mol, model)

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
   real(wp), parameter :: cn_max = 8.0_wp, trans(3, 1) = 0.0_wp, cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:)
   real(wp), allocatable :: energy(:), gradient(:, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :)
   real(wp) :: er, el, eps(3, 3), numsigma(3, 3), sigma(3, 3)

   call new_eeq2019_model(mol, model)

   rcov = get_covalent_rad(mol%num)
   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
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
         call get_coordination_number(mol, lattr, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
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
   real(wp), parameter :: cn_max = 8.0_wp, trans(3, 1) = 0.0_wp, cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:)
   real(wp), allocatable :: ql(:), qr(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)
   real(wp) :: er, el

   call new_eeq2019_model(mol, model)

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
   real(wp), parameter :: cn_max = 8.0_wp, trans(3, 1) = 0.0_wp, cutoff = 25.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :), rcov(:)
   real(wp), allocatable :: qr(:), ql(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :), numdL(:, :, :)
   real(wp) :: eps(3, 3)

   call new_eeq2019_model(mol, model)

   rcov = get_covalent_rad(mol%num)
   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
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
         call get_coordination_number(mol, lattr, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, qvec=qr)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cn, cut=cn_max)
         call model%solve(mol, cn, qvec=ql)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
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


subroutine test_q_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = &
      &[ 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831011E-1_wp, &
      &  4.92833325937898E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp, &
      &  6.61837152062313E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp, &
      &  2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172569E-2_wp, &
      & -3.73300777019193E-1_wp, 3.84585237785616E-2_wp,-5.05851088366941E-1_wp, &
      &  5.17677238544190E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call gen_test(error, mol, qref=ref)

end subroutine test_q_mb01


subroutine test_q_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: ref(16) = &
      &[ 7.38394711236237E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp, &
      & -7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp, &
      &  1.02748501676354E-1_wp, 9.47818107467041E-2_wp, 2.44260351729181E-2_wp, &
      &  2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp, &
      & -4.78119977010489E-1_wp, 6.57536027459279E-2_wp, 1.08259054549882E-1_wp, &
      & -3.58215329983396E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call gen_test(error, mol, qref=ref)

end subroutine test_q_mb02


subroutine test_e_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = &
      &[-1.13826381452353E-1_wp,-5.62509093753757E-1_wp, 2.40314612017535E-2_wp, &
      &  2.34612404176049E-1_wp, 3.24513140675984E-1_wp, 4.02366546606962E-2_wp, &
      & -2.17529354309074E-1_wp, 2.75364780143586E-2_wp, 4.02137437159417E-2_wp, &
      &  5.04840431402129E-2_wp,-3.53634596002012E-1_wp,-1.87985781866652E-1_wp, &
      & -2.52739878326955E-1_wp, 1.24520644464529E-1_wp, 2.69468100612926E-1_wp, &
      &  2.15919302769455E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call gen_test(error, mol, eref=ref)

end subroutine test_e_mb03


subroutine test_e_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = &
      &[ 5.48650310058159E-2_wp,-2.25780895465614E-1_wp, 4.35281555165058E-2_wp, &
      & -1.57205784670962E-1_wp, 4.09838051265821E-3_wp, 6.31282510503942E-2_wp, &
      &  7.48306240708851E-2_wp, 5.87730318612068E-2_wp, 6.10308471954074E-2_wp, &
      &  8.63933900356271E-2_wp,-9.99483735030189E-2_wp, 2.02497866273279E-1_wp, &
      &  3.47529219386489E-2_wp,-2.37058815227503E-1_wp, 6.74225110397230E-2_wp, &
      & -1.36552336490113E-1_wp]

   call get_structure(mol, "MB16-43", "04")
   call gen_test(error, mol, eref=ref)

end subroutine test_e_mb04


subroutine test_g_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numgrad(error, mol)

end subroutine test_g_mb05


subroutine test_g_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_numgrad(error, mol)

end subroutine test_g_mb06


subroutine test_s_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_numsigma(error, mol)

end subroutine test_s_mb07


subroutine test_s_mb08(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "08")
   call test_numsigma(error, mol)

end subroutine test_s_mb08


subroutine test_dqdr_mb09(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "09")
   call test_numdqdr(error, mol)

end subroutine test_dqdr_mb09


subroutine test_dqdr_mb10(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "10")
   call test_numdqdr(error, mol)

end subroutine test_dqdr_mb10


subroutine test_dqdL_mb11(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "11")
   call test_numdqdL(error, mol)

end subroutine test_dqdL_mb11


subroutine test_dqdL_mb12(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "12")
   call test_numdqdL(error, mol)

end subroutine test_dqdL_mb12


end module test_model
