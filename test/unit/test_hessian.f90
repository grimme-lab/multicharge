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

module test_hessian
   use mctc_env, only: wp
   use mctc_env_testing, only: new_unittest, unittest_type, error_type, test_failed
   use mctc_io, only: structure_type, new
   use mstore, only: get_structure
   use multicharge, only: get_charges, add_charge_hessian, cache_container
   use multicharge_model, only: mchrg_model_type
   use multicharge_param, only: new_eeq2019_model, new_eeqbc2025_model
   implicit none
   private

   public :: collect_hessian

   integer, parameter :: num_tri(3) = [30, 8, 1]
   integer, parameter :: num_di(2) = [8, 1]
   real(wp), parameter :: xyz_tri(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 3.0_wp, 0.1_wp, 0.2_wp, &
      & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3])
   real(wp), parameter :: xyz_wire(3, 3) = reshape([ &
      &  0.00_wp, 0.00_wp,  0.00_wp, &
      &  0.16_wp, 2.12_wp, -0.12_wp, &
      & -0.40_wp, 1.10_wp,  1.50_wp], [3, 3])
   real(wp), parameter :: xyz_slab(3, 3) = reshape([ &
      & 0.0_wp, -0.3_wp, 0.0_wp, &
      & 2.5_wp,  0.4_wp, 0.2_wp, &
      & 0.6_wp,  0.5_wp, 2.1_wp], [3, 3])
   real(wp), parameter :: xyz_half(3, 2) = reshape([ &
      & 0.0_wp, 0.0_wp, -0.2_wp, &
      & 3.0_wp, 0.2_wp,  0.3_wp], [3, 2])
   real(wp), parameter :: lat_skew(3, 3) = reshape([ &
      & 5.1_wp, 0.0_wp, 0.0_wp, &
      & 0.4_wp, 6.2_wp, 0.0_wp, &
      & 0.2_wp, 0.7_wp, 7.1_wp], [3, 3])
   real(wp), parameter :: lat_orth(3, 3) = reshape([ &
      & 6.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 7.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 8.0_wp], [3, 3])
   real(wp), parameter :: lat_wire(3, 3) = reshape([ &
      & 0.0_wp, 0.0_wp,  0.0_wp, &
      & 0.4_wp, 5.3_wp, -0.3_wp, &
      & 0.0_wp, 0.0_wp,  0.0_wp], [3, 3])
   real(wp), parameter :: lat_slab(3, 3) = reshape([ &
      & 5.2_wp,  0.1_wp, 0.3_wp, &
      & 0.0_wp,  0.0_wp, 0.0_wp, &
      & 0.5_wp, -0.2_wp, 6.1_wp], [3, 3])

contains


!> Collect all exported unit tests
subroutine collect_hessian(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("eeq-molecule", test_eeq), &
      & new_unittest("eeqbc-molecule", test_eeqbc), &
      & new_unittest("eeq-cn-cap", test_eeq_cap), &
      & new_unittest("eeqbc-cn-cap", test_eeqbc_cap), &
      & new_unittest("eeq-charged", test_eeq_charged), &
      & new_unittest("eeqbc-charged", test_eeqbc_charged), &
      & new_unittest("eeq-1d", test_eeq_1d), &
      & new_unittest("eeqbc-1d", test_eeqbc_1d), &
      & new_unittest("eeq-2d", test_eeq_2d), &
      & new_unittest("eeqbc-2d", test_eeqbc_2d), &
      & new_unittest("eeq-3d", test_eeq_3d), &
      & new_unittest("eeqbc-3d", test_eeqbc_3d), &
      & new_unittest("isolated-atom", test_atom), &
      & new_unittest("response-cache", test_cache), &
      & new_unittest("periodic-lattice-coincident", test_periodic_coincident), &
      & new_unittest("invalid-response", test_invalid) &
      & ]

end subroutine collect_hessian


subroutine test_eeq(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call check_model_response(error, mol, bond_capacity=.false.)

end subroutine test_eeq


subroutine test_eeqbc(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call check_model_response(error, mol, bond_capacity=.true.)

end subroutine test_eeqbc


subroutine test_eeq_cap(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call check_model_response(error, mol, bond_capacity=.false., cn_cut=2.5_wp)

end subroutine test_eeq_cap


subroutine test_eeqbc_cap(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call get_structure(mol, "MB16-43", "01")
   call check_model_response(error, mol, bond_capacity=.true., cn_cut=2.5_wp)

end subroutine test_eeqbc_cap


subroutine test_eeq_charged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call new(mol, num_tri, xyz_tri, charge=1.0_wp)
   call check_model_response(error, mol, bond_capacity=.false.)

end subroutine test_eeq_charged


subroutine test_eeqbc_charged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call new(mol, num_tri, xyz_tri, charge=1.0_wp)
   call check_model_response(error, mol, bond_capacity=.true.)

end subroutine test_eeqbc_charged


subroutine test_eeq_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call run_1d_cases(error, bond_capacity=.false.)

end subroutine test_eeq_1d


subroutine test_eeqbc_1d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call run_1d_cases(error, bond_capacity=.true.)

end subroutine test_eeqbc_1d


subroutine test_eeq_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call run_2d_cases(error, bond_capacity=.false.)

end subroutine test_eeq_2d


subroutine test_eeqbc_2d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call run_2d_cases(error, bond_capacity=.true.)

end subroutine test_eeqbc_2d


subroutine test_eeq_3d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call run_3d_cases(error, bond_capacity=.false.)

end subroutine test_eeq_3d


subroutine test_eeqbc_3d(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   call run_3d_cases(error, bond_capacity=.true.)

end subroutine test_eeqbc_3d


subroutine test_atom(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol

   call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp)
   call check_model_response(error, mol, bond_capacity=.false., label="EEQ isolated atom")
   if (allocated(error)) return

   call check_model_response(error, mol, bond_capacity=.true., label="EEQBC isolated atom")

end subroutine test_atom


!> Run 1D periodic Hessian checks (skew wire, half-cell tie, and charged EEQBC wire)
subroutine run_1d_cases(error, bond_capacity)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Select EEQBC (.true.) or EEQ (.false.)
   logical, intent(in) :: bond_capacity

   type(structure_type) :: mol

   ! Neutral 1D wire (FTF) with both on-axis and off-axis pairs
   call new(mol, num_tri, xyz_wire, lattice=lat_wire, &
      & periodic=[.false., .true., .false.])
   call check_model_response(error, mol, bond_capacity, label="1D wire")
   if (allocated(error)) return

   ! Half-cell displacement (TFF) with competing Wigner-Seitz images
   call new(mol, num_di, xyz_half, lattice=lat_orth, &
      & periodic=[.true., .false., .false.])
   call check_model_response(error, mol, bond_capacity, label="1D half-cell")
   if (allocated(error)) return

   if (bond_capacity) then
      call new(mol, num_tri, xyz_wire, charge=1.0_wp, lattice=lat_wire, &
         & periodic=[.false., .true., .false.])
      call check_model_response(error, mol, bond_capacity, label="1D charged wire")
   end if

end subroutine run_1d_cases


!> Run 2D periodic Hessian checks (skew slab, half-cell tie, and charged EEQBC slab)
subroutine run_2d_cases(error, bond_capacity)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Select EEQBC (.true.) or EEQ (.false.)
   logical, intent(in) :: bond_capacity

   type(structure_type) :: mol

   ! Neutral 2D slab (TFT)
   call new(mol, num_tri, xyz_slab, lattice=lat_slab, &
      & periodic=[.true., .false., .true.])
   call check_model_response(error, mol, bond_capacity, label="2D slab")
   if (allocated(error)) return

   ! Half-cell displacement (TTF) with competing Wigner-Seitz images
   call new(mol, num_di, xyz_half, lattice=lat_orth, &
      & periodic=[.true., .true., .false.])
   call check_model_response(error, mol, bond_capacity, label="2D half-cell")
   if (allocated(error)) return

   if (bond_capacity) then
      call new(mol, num_tri, xyz_slab, charge=-1.0_wp, lattice=lat_slab, &
         & periodic=[.true., .false., .true.])
      call check_model_response(error, mol, bond_capacity, label="2D charged slab")
   end if

end subroutine run_2d_cases


!> Run 3D periodic Hessian checks (charged skew cell, CN cap, crystal, atom, half-cell tie)
subroutine run_3d_cases(error, bond_capacity)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Select EEQBC (.true.) or EEQ (.false.)
   logical, intent(in) :: bond_capacity

   type(structure_type) :: mol

   call new(mol, num_tri, xyz_tri, charge=1.0_wp, lattice=lat_skew, &
      & periodic=[.true., .true., .true.])
   call check_model_response(error, mol, bond_capacity, label="3D charged skew")
   if (allocated(error)) return

   call check_model_response(error, mol, bond_capacity, cn_cut=2.5_wp, &
      & label="3D CN cap")
   if (allocated(error)) return

   call get_structure(mol, "X23", "urea")
   call check_model_response(error, mol, bond_capacity, label="3D urea")
   if (allocated(error)) return

   call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp, &
      & lattice=lat_skew, periodic=[.true., .true., .true.])
   call check_model_response(error, mol, bond_capacity, label="3D isolated atom")
   if (allocated(error)) return

   ! Half-cell displacement has two equivalent Wigner-Seitz images
   call new(mol, num_di, xyz_half, lattice=lat_orth, &
      & periodic=[.true., .true., .true.])
   call check_model_response(error, mol, bond_capacity, label="3D half-cell")

end subroutine run_3d_cases


!> Build the requested model and validate its weighted charge Hessian
subroutine check_model_response(error, mol, bond_capacity, cn_cut, label)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Select EEQBC (.true.) or EEQ (.false.)
   logical, intent(in) :: bond_capacity

   !> Optional smooth CN cutoff cap override
   real(wp), intent(in), optional :: cn_cut

   !> Optional diagnostic label prepended on failure
   character(len=*), intent(in), optional :: label

   class(mchrg_model_type), allocatable :: model

   if (bond_capacity) then
      call new_eeqbc2025_model(mol, model, error)
   else
      call new_eeq2019_model(mol, model, error)
   end if
   if (allocated(error)) return

   if (present(cn_cut)) model%ncoord%cut = cn_cut

   call check_response(error, mol, model)
   if (allocated(error) .and. present(label)) then
      error%message = label // ": " // error%message
   end if

end subroutine check_model_response


!> Independent reference: central differences of the weighted first response
subroutine check_response(error, mol, model)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Charge model instance
   class(mchrg_model_type), intent(in) :: model

   real(wp), parameter :: step = 1.0e-4_wp

   type(cache_container) :: cache
   real(wp) :: xyz, delta, scale
   real(wp), allocatable :: q(:), weights(:), hess(:, :), numeric(:, :), again(:, :)
   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :), grad(:, :), plus(:, :)
   integer :: n, i, j, a, col
   character(len=100) :: message

   n = mol%nat
   allocate(q(n), weights(n), hess(3 * n, 3 * n), numeric(3 * n, 3 * n), &
      & again(3 * n, 3 * n))
   allocate(dqdr(3, n, n), dqdL(3, 3, n), grad(3, n), plus(3, n))
   weights = [(sin(real(2 * i, wp)), i = 1, n)]

   call get_charges(model, mol, error, q, cache=cache)
   if (allocated(error)) return

   hess = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hess, cache)
   if (allocated(error)) return

   do j = 1, n
      do a = 1, 3
         col = 3 * (j - 1) + a
         xyz = mol%xyz(a, j)

         mol%xyz(a, j) = xyz + step
         call get_charges(model, mol, error, q, dqdr, dqdL)
         if (allocated(error)) return

         plus = 0.0_wp
         do i = 1, n
            plus = plus + weights(i) * dqdr(:, :, i)
         end do

         mol%xyz(a, j) = xyz - step
         call get_charges(model, mol, error, q, dqdr, dqdL)
         if (allocated(error)) return

         grad = 0.0_wp
         do i = 1, n
            grad = grad + weights(i) * dqdr(:, :, i)
         end do

         numeric(:, col) = reshape((plus - grad) / (2 * step), [3 * n])
         mol%xyz(a, j) = xyz
      end do
   end do

   delta = maxval(abs(hess - numeric))
   if (delta > 1.0e-6_wp) then
      write(message, "(a, es12.4)") &
         & "Weighted charge Hessian differs from finite differences: ", delta
      call test_failed(error, trim(message))
      return
   end if

   scale = max(1.0_wp, maxval(abs(hess)))
   if (maxval(abs(hess - transpose(hess))) > 1.0e-11_wp * scale) then
      call test_failed(error, "Charge Hessian is not symmetric")
      return
   end if

   do a = 1, 3
      if (maxval(abs(sum(hess(:, a:3 * n:3), dim=2))) > 1.0e-10_wp * scale) then
         call test_failed(error, "Charge Hessian violates translation invariance")
         return
      end if
   end do

   ! The same response cache may serve several weight vectors and accumulations
   do j = 1, 3 * n
      do i = 1, 3 * n
         numeric(i, j) = 0.125_wp + real(i + 2 * j, wp) * 1.0e-6_wp
      end do
   end do
   again = numeric
   call add_charge_hessian(model, mol, error, weights, again, cache)
   if (allocated(error)) return

   if (maxval(abs(again - hess - numeric)) > 1.0e-11_wp * scale) then
      call test_failed(error, "Charge Hessian does not accumulate")
      return
   end if

   again = 0.0_wp
   weights = 1.0_wp
   call add_charge_hessian(model, mol, error, weights, again, cache)
   if (allocated(error)) return

   if (maxval(abs(again)) > 1.0e-10_wp * scale) then
      call test_failed(error, "Second response violates total charge conservation")
   end if

end subroutine check_response


!> Refresh an existing cache for a different model and atom count, and exercise
!> the cache + Cartesian-response interface without a strain output.
subroutine test_cache(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   type(cache_container) :: cache, fresh
   class(mchrg_model_type), allocatable :: model
   real(wp) :: q1(1), q(3), refq(3), dqdr(3, 3, 3), refdr(3, 3, 3), dqdL(3, 3, 3)
   real(wp) :: hess(9, 9), refhess(9, 9), weights(3)

   call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   call get_charges(model, mol, error, q1, cache=cache)
   if (allocated(error)) return

   call new(mol, num_tri, xyz_tri, charge=1.0_wp)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return

   call get_charges(model, mol, error, q, dqdr=dqdr, cache=cache)
   if (allocated(error)) return

   call get_charges(model, mol, error, refq, refdr, dqdL, fresh)
   if (allocated(error)) return

   if (maxval(abs(q - refq)) > 1.0e-12_wp &
      & .or. maxval(abs(dqdr - refdr)) > 1.0e-12_wp) then
      call test_failed(error, "Refreshing the response cache changed charges or first response")
      return
   end if

   weights = [0.3_wp, -0.7_wp, 0.5_wp]
   hess = 0.0_wp
   refhess = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hess, cache)
   if (allocated(error)) return

   call add_charge_hessian(model, mol, error, weights, refhess, fresh)
   if (allocated(error)) return

   if (maxval(abs(hess - refhess)) > 1.0e-12_wp) then
      call test_failed(error, "Refreshing the response cache changed second response")
   end if

end subroutine test_cache


subroutine test_periodic_coincident(error)
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), parameter :: xyz(3, 2) = reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, &
      & 6.0_wp, 0.0_wp, 0.0_wp], [3, 2])

   type(structure_type) :: mol
   type(cache_container) :: cache
   class(mchrg_model_type), allocatable :: model
   real(wp) :: q(2), weights(2), hess(6, 6), scale
   integer :: a

   call new(mol, num_di, xyz, lattice=lat_orth, periodic=[.true., .true., .true.])
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   call get_charges(model, mol, error, q, cache=cache)
   if (allocated(error)) return

   weights = [0.4_wp, -0.9_wp]
   hess = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hess, cache)
   if (allocated(error)) return

   if (any(ieee_is_nan(hess))) then
      call test_failed(error, "Lattice-equivalent pair produced NaN in periodic Hessian")
      return
   end if

   scale = max(1.0_wp, maxval(abs(hess)))
   if (maxval(abs(hess - transpose(hess))) > 1.0e-11_wp * scale) then
      call test_failed(error, "Periodic coincident-image Hessian is not symmetric")
      return
   end if

   do a = 1, 3
      if (maxval(abs(sum(hess(:, a:6:3), dim=2))) > 1.0e-10_wp * scale) then
         call test_failed(error, "Periodic coincident-image Hessian violates translation invariance")
         return
      end if
   end do

end subroutine test_periodic_coincident


subroutine test_invalid(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(error_type), allocatable :: expected
   type(structure_type) :: mol
   type(cache_container) :: cache
   class(mchrg_model_type), allocatable :: model
   real(wp) :: hess(3, 3), q(1)

   call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]))
   call new_eeq2019_model(mol, model, error)
   if (allocated(error)) return

   hess = 0.125_wp
   call add_charge_hessian(model, mol, expected, [1.0_wp], hess, cache)
   if (.not. allocated(expected)) then
      call test_failed(error, "An unprepared response cache must be rejected")
      return
   end if

   call get_charges(model, mol, error, q, cache=cache)
   if (allocated(error)) return

   call add_charge_hessian(model, mol, expected, [1.0_wp, 2.0_wp], hess, cache)
   if (.not. allocated(expected)) then
      call test_failed(error, "Mismatched weight dimensions must be rejected")
      return
   end if

   if (any(hess /= 0.125_wp)) then
      call test_failed(error, "Rejected calls changed the Hessian")
   end if

end subroutine test_invalid


end module test_hessian
