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

module test_ncoord
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, &
      & test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_data, only : get_covalent_rad
   use multicharge_ncoord
   implicit none
   private

   public :: collect_ncoord

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr3 = sqrt(10*epsilon(1.0_wp))


contains


!> Collect all exported unit tests
subroutine collect_ncoord(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("cn-mb01", test_cn_mb01), &
      & new_unittest("cn-mb02", test_cn_mb02), &
      & new_unittest("cn-mb03", test_cn_mb03), &
      & new_unittest("cn-acetic", test_cn_acetic), &
      & new_unittest("dcndr-mb04", test_dcndr_mb04), &
      & new_unittest("dcndr-mb05", test_dcndr_mb05), &
      & new_unittest("dcndr-ammonia", test_dcndr_ammonia), &
      & new_unittest("dcndL-mb06", test_dcndL_mb06), &
      & new_unittest("dcndL-mb07", test_dcndL_mb07), &
      & new_unittest("dcndL-antracene", test_dcndL_anthracene) &
      & ]

end subroutine collect_ncoord


subroutine test_cn_gen(error, mol, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type) :: mol

   !> Reference CNs
   real(wp), intent(in) :: ref(:)

   real(wp), allocatable :: cn(:), rcov(:)
   real(wp), parameter :: cutoff = 30.0_wp
   real(wp), allocatable :: lattr(:, :)

   allocate(rcov(mol%nid), cn(mol%nat))
   rcov(:) = get_covalent_rad(mol%num)

   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
   call get_coordination_number(mol, lattr, cutoff, rcov, cn)

   if (any(abs(cn - ref) > thr)) then
      call test_failed(error, "Coordination numbers do not match")
      print'(3es21.14)', cn
   end if

end subroutine test_cn_gen


subroutine test_numgrad(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: iat, ic, mat
   real(wp), allocatable :: cn(:), rcov(:), cnr(:), cnl(:)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)
   real(wp), allocatable :: lattr(:, :)
   real(wp), parameter :: cutoff = 20.0_wp
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(rcov(mol%nid), cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))
   rcov(:) = get_covalent_rad(mol%num)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   if (any(mol%periodic)) then
      mat = min(mol%nat, 5)
   else
      mat = mol%nat
   end if

   do iat = 1, mat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call get_coordination_number(mol, lattr, cutoff, rcov, cnr)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call get_coordination_number(mol, lattr, cutoff, rcov, cnl)
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp*(cnr - cnl)/step
      end do
   end do

   call get_coordination_number(mol, lattr, cutoff, rcov, cn, dcndr, dcndL)

   if (any(abs(dcndr(:, :mat, :) - numdr(:, :mat, :)) > thr2)) then
      call test_failed(error, "Derivative of coordination number does not match")
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   integer :: ic, jc
   real(wp) :: eps(3, 3)
   real(wp), allocatable :: cn(:), rcov(:), cnr(:), cnl(:), xyz(:, :)
   real(wp), allocatable :: dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: numdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), trans(:, :)
   real(wp), parameter :: cutoff = 20.0_wp
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(rcov(mol%nid), cn(mol%nat), cnr(mol%nat), cnl(mol%nat), &
      & dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), xyz(3, mol%nat), &
      & numdL(3, 3, mol%nat))
   rcov(:) = get_covalent_rad(mol%num)
   call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   trans = lattr
   do ic = 1, 3
      do jc = 1, 3
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cnr)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call get_coordination_number(mol, lattr, cutoff, rcov, cnl)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp*(cnr - cnl)/step
      end do
   end do

   call get_coordination_number(mol, lattr, cutoff, rcov, cn, dcndr, dcndL)

   if (any(abs(dcndL - numdL) > thr3)) then
      call test_failed(error, "Derivative of coordination number does not match")
   end if

end subroutine test_numsigma


subroutine test_cn_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 4.03670396918677E+0_wp, 9.72798721502297E-1_wp, 1.98698465669657E+0_wp, &
      & 1.47312608051590E+0_wp, 9.97552155866795E-1_wp, 9.96862039916965E-1_wp, &
      & 1.45188437942218E+0_wp, 1.99267278111197E+0_wp, 3.84566220624764E+0_wp, &
      & 1.00242959599510E+0_wp, 9.96715113655073E-1_wp, 1.92505296745902E+0_wp, &
      & 4.62015142034058E+0_wp, 3.81973465175781E+0_wp, 3.95710919750442E+0_wp, &
      & 5.33862698412205E+0_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb01


subroutine test_cn_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 9.61099101791137E-1_wp, 3.87581247819995E+0_wp, 3.80155140067831E+0_wp, &
      & 2.96990277678560E+0_wp, 5.43508021969867E+0_wp, 1.01156705157372E+0_wp, &
      & 9.70139042949472E-1_wp, 9.72142268717279E-1_wp, 4.98780441573354E+0_wp, &
      & 1.01084927946071E+0_wp, 3.92876025928151E+0_wp, 3.88754303198463E+0_wp, &
      & 1.99577129500205E+0_wp, 9.71947229716782E-1_wp, 1.66031989216595E+0_wp, &
      & 1.97969868901054E+0_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb02


subroutine test_cn_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 3.99833201310984E+0_wp, 2.98041466005736E+0_wp, 9.93008506739515E-1_wp, &
      & 4.67413521470136E+0_wp, 7.44812343610600E+0_wp, 4.69855063323806E+0_wp, &
      & 4.85780221547050E+0_wp, 1.13980049581933E+0_wp, 9.92951731064869E-1_wp, &
      & 9.93802296429801E-1_wp, 1.98628603868569E+0_wp, 2.95465591146130E+0_wp, &
      & 4.83084657196178E+0_wp, 1.02417949516932E+0_wp, 3.46452255374201E+0_wp, &
      & 1.31433623523357E+0_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb03


subroutine test_cn_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(32) = [&
      & 9.85484970254893E-1_wp, 9.85511619286810E-1_wp, 9.85484986365689E-1_wp, &
      & 9.85508121225053E-1_wp, 9.93479325369907E-1_wp, 9.93474399260879E-1_wp, &
      & 9.93464233927038E-1_wp, 9.92795267281333E-1_wp, 9.92791891712794E-1_wp, &
      & 9.92789986522430E-1_wp, 9.92786608672290E-1_wp, 9.92912009999909E-1_wp, &
      & 9.92928947391201E-1_wp, 9.92923985134838E-1_wp, 9.92917815220921E-1_wp, &
      & 9.93477003743712E-1_wp, 2.99430549384771E+0_wp, 2.99430629828867E+0_wp, &
      & 2.99430089621814E+0_wp, 2.99430574791506E+0_wp, 3.97661752162159E+0_wp, &
      & 3.97662838199445E+0_wp, 3.97660644822657E+0_wp, 3.97661392666413E+0_wp, &
      & 1.97939425433928E+0_wp, 1.97942728307424E+0_wp, 1.97940242602914E+0_wp, &
      & 1.97942686736739E+0_wp, 1.00629029792463E+0_wp, 1.00629206096689E+0_wp, &
      & 1.00629410070858E+0_wp, 1.00630136304720E+0_wp]

   call get_structure(mol, "X23", "acetic")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_acetic


subroutine test_dcndr_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol)

end subroutine test_dcndr_mb04


subroutine test_dcndr_mb05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "05")
   call test_numgrad(error, mol)

end subroutine test_dcndr_mb05


subroutine test_dcndr_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "ammonia")
   call test_numgrad(error, mol)

end subroutine test_dcndr_ammonia


subroutine test_dcndL_mb06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol)

end subroutine test_dcndL_mb06


subroutine test_dcndL_mb07(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "07")
   call test_numsigma(error, mol)

end subroutine test_dcndL_mb07


subroutine test_dcndL_anthracene(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "X23", "anthracene")
   call test_numsigma(error, mol)

end subroutine test_dcndL_anthracene


end module test_ncoord
