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
   real(wp), parameter :: ref(16) = &
      &[4.03670002889443E+0_wp, 9.72798602236489E-1_wp, 1.98698459289250E+0_wp, &
      & 1.47312500109475E+0_wp, 9.97552138855144E-1_wp, 9.96862023238879E-1_wp, &
      & 1.45188312357461E+0_wp, 1.99267274348098E+0_wp, 3.84566163867660E+0_wp, &
      & 1.00242950582398E+0_wp, 9.96715095665523E-1_wp, 1.92505259119407E+0_wp, &
      & 4.62014881800631E+0_wp, 3.81973402238472E+0_wp, 3.95710899267725E+0_wp, &
      & 5.33862393788989E+0_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb01


subroutine test_cn_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = &
      &[9.61098938644516E-1_wp, 3.87581069538993E+0_wp, 3.80155064682641E+0_wp, &
      & 2.96990263040762E+0_wp, 5.43507838905789E+0_wp, 1.01156688025885E+0_wp, &
      & 9.70138913446538E-1_wp, 9.72142146912119E-1_wp, 4.98780161714987E+0_wp, &
      & 1.01084909581035E+0_wp, 3.92875855863300E+0_wp, 3.88754217159268E+0_wp, &
      & 1.99577127238747E+0_wp, 9.71947107159613E-1_wp, 1.66031882708281E+0_wp, &
      & 1.97969858871118E+0_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb02


subroutine test_cn_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = &
      &[3.99833054569977E+0_wp, 2.98041456425212E+0_wp, 9.93008463291392E-1_wp, &
      & 4.67413241067993E+0_wp, 7.44811983163412E+0_wp, 4.69854800105065E+0_wp, &
      & 4.85780104966086E+0_wp, 1.13979986956202E+0_wp, 9.92951691462952E-1_wp, &
      & 9.93802264983146E-1_wp, 1.98628597194324E+0_wp, 2.95465508647572E+0_wp, &
      & 4.83084536412428E+0_wp, 1.02417925644388E+0_wp, 3.46452093869485E+0_wp, &
      & 1.31433531152694E+0_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_cn_gen(error, mol, ref)

end subroutine test_cn_mb03


subroutine test_cn_acetic(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(32) = &
      &[9.85484827652486E-1_wp, 9.85511476874469E-1_wp, 9.85484843883598E-1_wp, &
      & 9.85507978837200E-1_wp, 9.93479292804725E-1_wp, 9.93474366672860E-1_wp, &
      & 9.93464201291902E-1_wp, 9.92795231559757E-1_wp, 9.92791855975717E-1_wp, &
      & 9.92789950776612E-1_wp, 9.92786572910960E-1_wp, 9.92911974815435E-1_wp, &
      & 9.92928912284626E-1_wp, 9.92923950005444E-1_wp, 9.92917780063145E-1_wp, &
      & 9.93476971167766E-1_wp, 2.99430546284162E+0_wp, 2.99430626728671E+0_wp, &
      & 2.99430086518903E+0_wp, 2.99430571691001E+0_wp, 3.97661738191383E+0_wp, &
      & 3.97662824232731E+0_wp, 3.97660630846656E+0_wp, 3.97661378694019E+0_wp, &
      & 1.97939413938006E+0_wp, 1.97942716824100E+0_wp, 1.97940231110394E+0_wp, &
      & 1.97942675253856E+0_wp, 1.00629024113779E+0_wp, 1.00629200417302E+0_wp, &
      & 1.00629404390296E+0_wp, 1.00630130618483E+0_wp]

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
