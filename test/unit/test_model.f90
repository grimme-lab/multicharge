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
   use mctc_io_structure, only : structure_type, new
   use mstore, only : get_structure
   use multicharge_model, only: mchrg_model_type
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
      & new_unittest("charges-actinides", test_q_actinides), &
      & new_unittest("energy-mb03", test_e_mb03), &
      & new_unittest("energy-mb04", test_e_mb04), &
      & new_unittest("gradient-mb05", test_g_mb05), &
      & new_unittest("gradient-mb06", test_g_mb06), &
      & new_unittest("sigma-mb07", test_s_mb07), &
      & new_unittest("sigma-mb08", test_s_mb08), &
      & new_unittest("dqdr-mb09", test_dqdr_mb09), &
      & new_unittest("dqdr-mb10", test_dqdr_mb10), &
      & new_unittest("dqdL-mb11", test_dqdL_mb11), &
      & new_unittest("dqdL-mb12", test_dqdL_mb12), &
      & new_unittest("gradient-h2plus", test_g_h2plus), &
      & new_unittest("gradient-znooh", test_g_znooh), &
      & new_unittest("dqdr-znooh", test_dqdr_znooh) &
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), allocatable :: cn(:)
   real(wp), allocatable :: energy(:)
   real(wp), allocatable :: qvec(:)

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   allocate(cn(mol%nat))

   call model%ncoord%get_coordination_number(mol, trans, cn)

   if (present(eref)) then
      allocate(energy(mol%nat))
      energy(:) = 0.0_wp
   end if
   if (present(qref)) then
      allocate(qvec(mol%nat))
   end if

   call model%solve(mol, error, cn, energy=energy, qvec=qvec)
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :), sigma(:, :)
   real(wp), allocatable :: numgrad(:, :)
   real(wp) :: er, el

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & energy(mol%nat), gradient(3, mol%nat), sigma(3, 3), numgrad(3, mol%nat))
   energy(:) = 0.0_wp
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         energy(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%solve(mol, error, cn, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%solve(mol, error, cn, energy=energy)
         if (allocated(error)) exit lp
         el = sum(energy)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)

   energy(:) = 0.0_wp
   call model%solve(mol, error, cn, dcndr, dcndL, energy, gradient, sigma)
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: energy(:), gradient(:, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :)
   real(wp) :: er, el, eps(3, 3), numsigma(3, 3), sigma(3, 3)

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

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
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%solve(mol, error, cn, energy=energy)
         if (allocated(error)) exit lp
         er = sum(energy)

         energy(:) = 0.0_wp
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%solve(mol, error, cn, energy=energy)
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

   energy(:) = 0.0_wp
   call model%solve(mol, error, cn, dcndr, dcndL, energy, gradient, sigma)
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: ql(:), qr(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: numdr(:, :, :)

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   allocate(cn(mol%nat), dcndr(3, mol%nat, mol%nat), dcndL(3, 3, mol%nat), &
      & ql(mol%nat), qr(mol%nat), dqdr(3, mol%nat, mol%nat), dqdL(3, 3, mol%nat), &
      & numdr(3, mol%nat, mol%nat))

   lp: do iat = 1, mol%nat
      do ic = 1, 3
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%solve(mol, error, cn, qvec=qr)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call model%ncoord%get_coordination_number(mol, trans, cn)
         call model%solve(mol, error, cn, qvec=ql)
         if (allocated(error)) exit lp

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numdr(ic, iat, :) = 0.5_wp*(qr - ql)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)

   call model%solve(mol, error, cn, dcndr, dcndL, dqdr=dqdr, dqdL=dqdL)
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
   real(wp), parameter :: trans(3, 1) = 0.0_wp
   real(wp), parameter :: step = 1.0e-6_wp, unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), allocatable :: cn(:), dcndr(:, :, :), dcndL(:, :, :)
   real(wp), allocatable :: qr(:), ql(:), dqdr(:, :, :), dqdL(:, :, :)
   real(wp), allocatable :: lattr(:, :), xyz(:, :), numdL(:, :, :)
   real(wp) :: eps(3, 3)

   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

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
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%solve(mol, error, cn, qvec=qr)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         lattr(:, :) = matmul(eps, trans)
         call model%ncoord%get_coordination_number(mol, lattr, cn)
         call model%solve(mol, error, cn, qvec=ql)
         if (allocated(error)) exit lp

         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         lattr(:, :) = trans
         numdL(jc, ic, :) = 0.5_wp*(qr - ql)/step
      end do
   end do lp
   if (allocated(error)) return

   call model%ncoord%get_coordination_number(mol, trans, cn, dcndr, dcndL)

   call model%solve(mol, error, cn, dcndr, dcndL, dqdr=dqdr, dqdL=dqdL)
   if (allocated(error)) return

   if (any(abs(dqdL(:, :, :) - numdL(:, :, :)) > thr2)) then
      call test_failed(error, "Derivative of charges does not match")
   end if

end subroutine test_numdqdL


subroutine test_q_mb01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 7.73347759615437E-1_wp, 1.07626897257271E-1_wp,-3.66999554268267E-1_wp, &
      & 4.92833775451616E-2_wp,-1.83332153007808E-1_wp, 2.33302084420314E-1_wp, &
      & 6.61837602813735E-2_wp,-5.43944147972069E-1_wp,-2.70264297953247E-1_wp, &
      & 2.66618970100409E-1_wp, 2.62725030332215E-1_wp,-7.15315061940473E-2_wp, &
      &-3.73300836681036E-1_wp, 3.84585142200261E-2_wp,-5.05851076468890E-1_wp, &
      & 5.17677178773158E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call gen_test(error, mol, qref=ref)

end subroutine test_q_mb01


subroutine test_q_mb02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: ref(16) = [&
      & 7.38394752482521E-2_wp,-1.68354859084778E-1_wp,-3.47642846358022E-1_wp, &
      &-7.05489251302223E-1_wp, 7.73548241620680E-1_wp, 2.30207580650128E-1_wp, &
      & 1.02748505731185E-1_wp, 9.47818154871089E-2_wp, 2.44259536057649E-2_wp, &
      & 2.34984928231320E-1_wp,-3.17839956573785E-1_wp, 6.67112952465234E-1_wp, &
      &-4.78119957747208E-1_wp, 6.57536208287042E-2_wp, 1.08259091466373E-1_wp, &
      &-3.58215294268738E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call gen_test(error, mol, qref=ref)

end subroutine test_q_mb02


subroutine test_q_actinides(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   real(wp), parameter :: ref(17) = [&
      & 1.86904766283711E-02_wp, 2.89972818160259E-01_wp, 3.59298070941105E-02_wp, & 
      &-4.61256458126589E-02_wp,-7.02605348653647E-02_wp,-7.42052215689073E-02_wp, &
      &-8.21938718945845E-02_wp, 1.64953118841151E-01_wp, 2.10381640633390E-01_wp, &
      &-6.65485355096282E-02_wp,-2.59873890255450E-01_wp, 1.33839147940414E-01_wp, &
      & 7.20768968601809E-02_wp,-3.36652347675997E-03_wp,-1.14546280789657E-01_wp, &
      &-8.55922398441004E-02_wp,-1.23131162140762E-01_wp]

   !> Molecular structure data 
   mol%nat = 17
   mol%nid = 17
   mol%id = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, &
      & 12, 13, 14, 15, 16, 17]
   mol%num = [87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, &
      & 98, 99, 100, 101, 102, 103]
   mol%xyz = reshape([ &
      & 0.98692316414074_wp, 6.12727238368797_wp,-6.67861597188102_wp, &
      & 3.63898862390869_wp, 5.12109301182962_wp, 3.01908613326278_wp, &
      & 5.14503571563551_wp,-3.97172984617710_wp, 3.82011791828867_wp, &
      & 6.71986847575494_wp, 1.71382138402812_wp, 3.92749159076307_wp, &
      & 4.13783589704826_wp,-2.10695793491818_wp, 0.19753203068899_wp, &
      & 8.97685097698326_wp,-3.08813636191844_wp,-4.45568615593938_wp, &
      & 12.5486412940776_wp,-1.77128765259458_wp, 0.59261498922861_wp, &
      & 7.82051475868325_wp,-3.97159756604558_wp,-0.53637703616916_wp, &
      &-0.43444574624893_wp,-1.69696511583960_wp,-1.65898182093050_wp, &
      &-4.71270645149099_wp,-0.11534827468942_wp, 2.84863373521297_wp, &
      &-2.52061680335614_wp, 1.82937752749537_wp,-2.10366982879172_wp, &
      & 0.13551154616576_wp, 7.99805359235043_wp,-1.55508522619903_wp, &
      & 3.91594542499717_wp,-1.72975169129597_wp,-5.07944366756113_wp, &
      &-1.03393930231679_wp, 4.69307230054046_wp, 0.02656940927472_wp, &
      & 6.20675384557240_wp, 4.24490721493632_wp,-0.71004195169885_wp, &
      & 7.04586341131562_wp, 5.20053667939076_wp,-7.51972863675876_wp, &
      & 2.01082807362334_wp, 1.34838807211157_wp,-4.70482633508447_wp],&
      & [3, 17])
   mol%periodic = [.false.]

   call gen_test(error, mol, qref=ref)

end subroutine test_q_actinides


subroutine test_e_mb03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      &-2.18698345033562E-1_wp,-1.04793885931268E+0_wp, 4.78963353574572E-2_wp, &
      & 5.76566377591676E-1_wp, 7.37187470977927E-1_wp, 8.06020047053305E-2_wp, &
      &-4.19837955782898E-1_wp, 5.49627510550566E-2_wp, 8.01486728591565E-2_wp, &
      & 1.00618944521776E-1_wp,-6.61715169034150E-1_wp,-3.60531647289563E-1_wp, &
      &-4.87729666337974E-1_wp, 2.48257554279938E-1_wp, 6.96027176590956E-1_wp, &
      & 4.31679925875087E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call gen_test(error, mol, eref=ref)

end subroutine test_e_mb03


subroutine test_e_mb04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: ref(16) = [&
      & 1.13974214746111E-1_wp,-4.41735365367827E-1_wp, 8.99036489938394E-2_wp, &
      &-2.97539904703271E-1_wp, 8.05174117097006E-3_wp, 1.31105783760276E-1_wp, &
      & 1.54594451996644E-1_wp, 1.19929653841255E-1_wp, 1.26056586309101E-1_wp, &
      & 1.78439005754586E-1_wp,-1.98703462666082E-1_wp, 4.19630120027785E-1_wp, &
      & 7.05569220334930E-2_wp,-4.50925107441869E-1_wp, 1.39289602382354E-1_wp, &
      &-2.67853086061429E-1_wp]

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


subroutine test_g_h2plus(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   integer, parameter :: nat = 2
   real(wp), parameter :: charge = 1.0_wp
   integer, parameter :: num(nat) = [1, 1]
   real(wp), parameter :: xyz(3, nat) = reshape([ &
      & +0.00000000000000_wp, +0.00000000000000_wp, +0.00000000000000_wp, &
      & +1.00000000000000_wp, +0.00000000000000_wp, +0.00000000000000_wp],&
      & [3, nat])

   call new(mol, num, xyz, charge)
   call test_numgrad(error, mol)

end subroutine test_g_h2plus


subroutine test_g_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

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
   call test_numgrad(error, mol)

end subroutine test_g_znooh


subroutine test_dqdr_znooh(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

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
   call test_numdqdr(error, mol)

end subroutine test_dqdr_znooh


end module test_model
