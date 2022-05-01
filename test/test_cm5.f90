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

module test_cm5
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_data, only : get_covalent_rad
   use multicharge_cm5
   implicit none
   private

   public :: collect_cm5

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = sqrt(10*epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_cm5(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("cm5-mb01", test_cm5_mb01), &
      & new_unittest("cm5-mb02", test_cm5_mb02), &
      & new_unittest("cm5-mb03", test_cm5_mb03), &
      & new_unittest("cm5-acetic", test_cm5_acetic), &
      & new_unittest("dcm5dr-mb04", test_dcm5dr_mb04), &
      & new_unittest("dcm5dr-mb05", test_dcm5dr_mb05), &
      & new_unittest("dcm5dr-ammonia", test_dcm5dr_ammonia), &
      & new_unittest("dcm5dL-mb06", test_dcm5dL_mb06), &
      & new_unittest("dcm5dL-mb07", test_dcm5dL_mb07), &
      & new_unittest("dcm5dL-antracene", test_dcm5dL_anthracene) &
      & ]

end subroutine collect_cm5


subroutine test_cm5_gen(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Reference cm5s
   real(wp), intent(in) :: ref(:)

   type(cm5_type) :: calc
   real(wp), allocatable :: cm5(:), rad(:)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)

   allocate(rad(mol%nid), cm5(mol%nat))
   rad(:) = get_covalent_rad(mol%num)
   call new_cm5(calc, mol, rad)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call calc%get_correction(mol, lattr, cutoff, cm5)

   if (any(abs(cm5 - ref) > thr)) then
      call test_failed(error, "Coordination numbers do not match")
      print'(3es21.14)', cm5
   end if

end subroutine test_cm5_gen


subroutine test_numgrad(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   type(cm5_type) :: calc
   integer :: iat, ic, mat
   real(wp), allocatable :: cm5(:), rad(:), cm5r(:), cm5l(:)
   real(wp), allocatable :: dcm5dr(:, :, :), dcm5dL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)
   real(wp), allocatable :: lattr(:, :)
   real(wp), parameter :: cutoff = 20.0_wp
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(rad(mol%nid), cm5(mol%nat), cm5r(mol%nat), cm5l(mol%nat), &
      & dcm5dr(3, mol%nat, mol%nat), dcm5dL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))
   rad(:) = get_covalent_rad(mol%num)
   call new_cm5(calc, mol, rad)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   if (any(mol%periodic)) then
      mat = min(mol%nat, 5)
   else
      mat = mol%nat
   end if

   do iat = 1, mat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call calc%get_correction(mol, lattr, cutoff, cm5r)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call calc%get_correction(mol, lattr, cutoff, cm5l)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp*(cm5r - cm5l)/step
      end do
   end do

   call calc%get_correction(mol, lattr, cutoff, cm5, dcm5dr, dcm5dL)

   if (any(abs(dcm5dr(:, :mat, :) - numdr(:, :mat, :)) > thr2)) then
      call test_failed(error, "Derivative of coordination number does not match")
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   type(cm5_type) :: calc
   integer :: ic, jc
   real(wp) :: eps(3, 3)
   real(wp), allocatable :: cm5(:), rad(:), cm5r(:), cm5l(:), xyz(:, :)
   real(wp), allocatable :: dcm5dr(:, :, :), dcm5dL(:, :, :)
   real(wp), allocatable :: numdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), trans(:, :)
   real(wp), parameter :: cutoff = 20.0_wp
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(rad(mol%nid), cm5(mol%nat), cm5r(mol%nat), cm5l(mol%nat), &
      & dcm5dr(3, mol%nat, mol%nat), dcm5dL(3, 3, mol%nat), xyz(3, mol%nat), &
      & numdL(3, 3, mol%nat))
   rad(:) = get_covalent_rad(mol%num)
   call new_cm5(calc, mol, rad)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   trans = lattr
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call calc%get_correction(mol, lattr, cutoff, cm5r)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call calc%get_correction(mol, lattr, cutoff, cm5l)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp*(cm5r - cm5l)/step
      end do
   end do

   call calc%get_correction(mol, lattr, cutoff, cm5, dcm5dr, dcm5dL)

   if (any(abs(dcm5dL - numdL) > thr3)) then
      call test_failed(error, "Derivative of coordination number does not match")
      print '(3es20.13)', dcm5dl - numdl
   end if

end subroutine test_numsigma


subroutine test_cm5_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 1.11895897328988E+0_wp, 2.16214959962204E-1_wp,-3.10468325172267E-1_wp, &
      & 1.96058146036230E-1_wp, 1.44402110173228E-1_wp, 4.62343021143816E-1_wp, &
      & 2.74877218633135E-1_wp,-3.96863410125448E-1_wp,-5.40293891904070E-1_wp, &
      & 3.82012512591119E-1_wp, 4.69916546063708E-1_wp, 2.31082727875368E-1_wp, &
      &-8.53576733849155E-1_wp,-4.54932092996253E-1_wp,-9.75254887961956E-1_wp, &
      & 3.55231262404571E-2_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_cm5_gen(error, mol, ref)

end subroutine test_cm5_mb01


subroutine test_cm5_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 2.11395850273080E-1_wp, 3.43105731316906E-2_wp,-6.57669452971921E-1_wp, &
      &-9.78417480427598E-1_wp, 1.07832086292190E+0_wp, 3.54023220587265E-1_wp, &
      & 2.18958685899068E-1_wp, 2.19511745916536E-1_wp,-5.14398778582881E-1_wp, &
      & 3.52856876925939E-1_wp,-7.35967645481895E-1_wp, 7.73145275809223E-1_wp, &
      &-6.62311788738417E-1_wp, 1.95480555279679E-1_wp, 2.93980077224277E-1_wp, &
      &-1.83218577765942E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_cm5_gen(error, mol, ref)

end subroutine test_cm5_mb02


subroutine test_cm5_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      &-1.34373509998640E-1_wp,-1.19534493468583E+0_wp, 1.26438437442476E-1_wp, &
      & 7.56508218833743E-1_wp, 1.07381780492030E+0_wp,-6.11505111992222E-1_wp, &
      &-3.33920024757557E-1_wp, 1.31612681861290E-1_wp, 1.30754366335334E-1_wp, &
      & 1.35179641962814E-1_wp,-6.84575253657633E-1_wp,-2.68387943164373E-1_wp, &
      &-4.21338979620731E-1_wp, 3.51291950392921E-1_wp, 8.26286829970604E-1_wp, &
      & 1.17555826157508E-1_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_cm5_gen(error, mol, ref)

end subroutine test_cm5_mb03


subroutine test_cm5_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(32) = [&
      & 4.19634302526966E-1_wp, 4.19686161643937E-1_wp, 4.19639315842737E-1_wp, &
      & 4.19679308743018E-1_wp, 1.51237758195438E-1_wp, 1.51228260058246E-1_wp, &
      & 1.51208370041931E-1_wp, 1.39678754693929E-1_wp, 1.39673463518620E-1_wp, &
      & 1.39667650445828E-1_wp, 1.39662715696985E-1_wp, 1.41331304306604E-1_wp, &
      & 1.41365692572454E-1_wp, 1.41352367081752E-1_wp, 1.41343242891877E-1_wp, &
      & 1.51233291612632E-1_wp, 1.43881363825290E-1_wp, 1.43898657679985E-1_wp, &
      & 1.43888970348442E-1_wp, 1.43875598498826E-1_wp,-3.33281760711736E-1_wp, &
      &-3.33296632639702E-1_wp,-3.33263567313097E-1_wp,-3.33272141246093E-1_wp, &
      &-4.30479914692791E-1_wp,-4.30532521040603E-1_wp,-4.30482005460891E-1_wp, &
      &-4.30531782858980E-1_wp,-2.31984922325696E-1_wp,-2.32020594433292E-1_wp, &
      &-2.32014740280446E-1_wp,-2.32005967222174E-1_wp]

   call get_structure(mol, "X23", "acetic")
   call test_cm5_gen(error, mol, ref)

end subroutine test_cm5_acetic


subroutine test_dcm5dr_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol)

end subroutine test_dcm5dr_mb04


subroutine test_dcm5dr_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numgrad(error, mol)

end subroutine test_dcm5dr_mb05


subroutine test_dcm5dr_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "ammonia")
   call test_numgrad(error, mol)

end subroutine test_dcm5dr_ammonia


subroutine test_dcm5dL_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol)

end subroutine test_dcm5dL_mb06


subroutine test_dcm5dL_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_numsigma(error, mol)

end subroutine test_dcm5dL_mb07


subroutine test_dcm5dL_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "anthracene")
   call test_numsigma(error, mol)

end subroutine test_dcm5dL_anthracene


end module test_cm5
