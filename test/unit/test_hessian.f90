! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0

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
contains

subroutine collect_hessian(testsuite)
   type(unittest_type), allocatable, intent(out) :: testsuite(:)
   testsuite = [ &
      & new_unittest("eeq-molecule", test_eeq), &
      & new_unittest("eeqbc-molecule", test_eeqbc), &
      & new_unittest("eeq-cn-cap", test_eeq_cap), &
      & new_unittest("eeqbc-cn-cap", test_eeqbc_cap), &
      & new_unittest("eeq-charged", test_eeq_charged), &
      & new_unittest("eeqbc-charged", test_eeqbc_charged), &
      & new_unittest("eeq-periodic", test_eeq_periodic), &
      & new_unittest("eeqbc-periodic", test_eeqbc_periodic), &
      & new_unittest("isolated-atom", test_atom), &
      & new_unittest("response-cache", test_cache), &
      & new_unittest("invalid-response", test_invalid) &
      & ]
end subroutine collect_hessian

subroutine test_eeq(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .false., 1)
end subroutine test_eeq

subroutine test_eeqbc(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .true., 1)
end subroutine test_eeqbc

subroutine test_eeq_cap(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .false., 2)
end subroutine test_eeq_cap

subroutine test_eeqbc_cap(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .true., 2)
end subroutine test_eeqbc_cap

subroutine test_eeq_charged(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .false., 3)
end subroutine test_eeq_charged

subroutine test_eeqbc_charged(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .true., 3)
end subroutine test_eeqbc_charged

subroutine test_eeq_periodic(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: icase
   do icase = 5, 9
      call run_case(error, .false., icase)
      if (allocated(error)) return
   end do
end subroutine test_eeq_periodic

subroutine test_eeqbc_periodic(error)
   type(error_type), allocatable, intent(out) :: error
   integer :: icase
   do icase = 5, 9
      call run_case(error, .true., icase)
      if (allocated(error)) return
   end do
end subroutine test_eeqbc_periodic

subroutine test_atom(error)
   type(error_type), allocatable, intent(out) :: error
   call run_case(error, .false., 4)
   if (allocated(error)) return
   call run_case(error, .true., 4)
end subroutine test_atom

subroutine run_case(error, bond_capacity, icase)
   type(error_type), allocatable, intent(out) :: error
   logical, intent(in) :: bond_capacity
   integer, intent(in) :: icase
   type(structure_type) :: mol
   class(mchrg_model_type), allocatable :: model
   real(wp), parameter :: lattice(3, 3) = reshape([ &
      & 5.1_wp, 0.0_wp, 0.0_wp, 0.4_wp, 6.2_wp, 0.0_wp, &
      & 0.2_wp, 0.7_wp, 7.1_wp], [3, 3])
   character(len=30) :: label

   select case(icase)
   case(1, 2)
      call get_structure(mol, "MB16-43", "01")
   case(3)
      call new(mol, [30, 8, 1], reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
         & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), charge=1.0_wp)
   case(4)
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp)
   case(5, 6)
      call new(mol, [30, 8, 1], reshape([ &
         & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
         & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), charge=1.0_wp, &
         & lattice=lattice, periodic=[.true., .true., .true.])
   case(7)
      call get_structure(mol, "X23", "urea")
   case(8)
      call new(mol, [8], reshape([0.0_wp, 0.0_wp, 0.0_wp], [3, 1]), charge=-1.0_wp, &
         & lattice=lattice, periodic=[.true., .true., .true.])
   case(9)
      ! The half-cell displacement has two equivalent WS images.
      call new(mol, [8, 1], reshape([0.0_wp, 0.0_wp, 0.0_wp, &
         & 3.0_wp, 0.2_wp, 0.3_wp], [3, 2]), lattice=reshape([ &
         & 6.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 7.0_wp, 0.0_wp, &
         & 0.0_wp, 0.0_wp, 8.0_wp], [3, 3]), periodic=[.true., .true., .true.])
   case default
      call test_failed(error, "Unknown charge Hessian test case")
      return
   end select
   if (bond_capacity) then
      call new_eeqbc2025_model(mol, model, error)
   else
      call new_eeq2019_model(mol, model, error)
   end if
   if (allocated(error)) return
   ! Force an appreciable smooth-cap curvature in both models.
   if (icase == 2 .or. icase == 6) model%ncoord%cut = 2.5_wp
   call check_response(error, mol, model)
   if (allocated(error)) then
      write(label, "(a, i0, a)") "Case ", icase, ": "
      error%message = trim(label)//" "//error%message
   end if
end subroutine run_case

!> Independent reference: central differences of the weighted first response.
subroutine check_response(error, mol, model)
   type(error_type), allocatable, intent(out) :: error
   type(structure_type), intent(inout) :: mol
   class(mchrg_model_type), intent(in) :: model
   type(cache_container) :: cache
   real(wp), parameter :: step = 1.0e-4_wp
   real(wp) :: xyz, delta, scale
   real(wp), allocatable :: q(:), weights(:), hess(:, :), numeric(:, :), again(:, :)
   real(wp), allocatable :: dqdr(:, :, :), dqdL(:, :, :), grad(:, :), plus(:, :)
   integer :: n, i, j, a, col
   character(len=100) :: message

   n = mol%nat
   allocate(q(n), weights(n), hess(3*n, 3*n), numeric(3*n, 3*n), again(3*n, 3*n))
   allocate(dqdr(3, n, n), dqdL(3, 3, n), grad(3, n), plus(3, n))
   weights = [(sin(real(2*i, wp)), i=1, n)]
   call get_charges(model, mol, error, q, cache=cache)
   if (allocated(error)) return
   hess = 0.0_wp
   call add_charge_hessian(model, mol, error, weights, hess, cache)
   if (allocated(error)) return
   do j = 1, n
      do a = 1, 3
         col = 3*(j-1)+a
         xyz = mol%xyz(a, j)
         mol%xyz(a, j) = xyz+step
         call get_charges(model, mol, error, q, dqdr, dqdL)
         if (allocated(error)) return
         plus = 0.0_wp
         do i = 1, n
            plus = plus+weights(i)*dqdr(:, :, i)
         end do
         mol%xyz(a, j) = xyz-step
         call get_charges(model, mol, error, q, dqdr, dqdL)
         if (allocated(error)) return
         grad = 0.0_wp
         do i = 1, n
            grad = grad+weights(i)*dqdr(:, :, i)
         end do
         numeric(:, col) = reshape((plus-grad)/(2*step), [3*n])
         mol%xyz(a, j) = xyz
      end do
   end do
   delta = maxval(abs(hess-numeric))
   if (delta > 1.0e-6_wp) then
      write(message, "(a, es12.4)") "Weighted charge Hessian differs from finite differences: ", delta
      call test_failed(error, trim(message))
      return
   end if
   scale = max(1.0_wp, maxval(abs(hess)))
   if (maxval(abs(hess-transpose(hess))) > 1.0e-11_wp*scale) then
      call test_failed(error, "Charge Hessian is not symmetric")
      return
   end if
   do a = 1, 3
      if (maxval(abs(sum(hess(:, a:3*n:3), dim=2))) > 1.0e-10_wp*scale) then
         call test_failed(error, "Charge Hessian violates translation invariance")
         return
      end if
   end do
   ! The same response cache may serve several weight vectors and accumulations.
   again = 0.125_wp
   call add_charge_hessian(model, mol, error, weights, again, cache)
   if (allocated(error)) return
   if (maxval(abs(again-hess-0.125_wp)) > 1.0e-11_wp*scale) then
      call test_failed(error, "Charge Hessian does not accumulate")
      return
   end if
   again = 0.0_wp
   weights = 1.0_wp
   call add_charge_hessian(model, mol, error, weights, again, cache)
   if (allocated(error)) return
   if (maxval(abs(again)) > 1.0e-10_wp*scale) then
      call test_failed(error, "Second response violates total charge conservation")
   end if
end subroutine check_response

!> Refresh an existing cache for a different model and atom count, and exercise
!> the cache + Cartesian-response interface without a strain output.
subroutine test_cache(error)
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
   call new(mol, [30, 8, 1], reshape([ &
      & 0.0_wp, 0.0_wp, 0.0_wp, 3.0_wp, 0.1_wp, 0.2_wp, &
      & 0.5_wp, 2.0_wp, 0.3_wp], [3, 3]), charge=1.0_wp)
   call new_eeqbc2025_model(mol, model, error)
   if (allocated(error)) return
   call get_charges(model, mol, error, q, dqdr=dqdr, cache=cache)
   if (allocated(error)) return
   call get_charges(model, mol, error, refq, refdr, dqdL, fresh)
   if (allocated(error)) return
   if (maxval(abs(q-refq)) > 1.0e-12_wp .or. maxval(abs(dqdr-refdr)) > 1.0e-12_wp) then
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
   if (maxval(abs(hess-refhess)) > 1.0e-12_wp) then
     call test_failed(error, "Refreshing the response cache changed second response")
   end if
end subroutine test_cache

subroutine test_invalid(error)
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
   mol%periodic = [.true., .false., .false.]
   call add_charge_hessian(model, mol, expected, [1.0_wp], hess, cache)
   if (.not. allocated(expected)) then
      call test_failed(error, "Partial periodicity is not supported by the 3D response")
      return
   end if
   if (any(hess /= 0.125_wp)) call test_failed(error, "Rejected calls changed the Hessian")
end subroutine test_invalid

end module test_hessian
