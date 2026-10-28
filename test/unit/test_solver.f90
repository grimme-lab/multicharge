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

module test_solver
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, test_failed
   use multicharge_solver, only : new_mchrg_solver
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, &
      & direct_input
   use multicharge_solver_type, only : mchrg_solver_type, mchrg_solver_input
   implicit none
   private

   public :: collect_solver

   real(wp), parameter :: thr = 100 * epsilon(1.0_wp)
   real(wp), parameter :: thr1 = 1.0e5_wp*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))
   real(wp), parameter :: thr_rel = 1.0e-6_wp

   !> Solver input unknown to the solver factory
   type, extends(mchrg_solver_input) :: unknown_input
   end type unknown_input


contains


!> Collect all unit tests for the CG solver
subroutine collect_solver(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      & new_unittest("cg-identity-2x2", test_cg_identity_2x2), &
      & new_unittest("cg-diagonal-5x5", test_cg_diagonal_5x5), &
      & new_unittest("cg-spd-small", test_cg_spd_small), &
      & new_unittest("cg-spd-large", test_cg_spd_large), &
      & new_unittest("cg-ill-conditioned", test_cg_ill_conditioned), &
      & new_unittest("cg-zero-rhs", test_cg_zero_rhs), &
      & new_unittest("cg-random-spd", test_cg_random_spd), &
      & new_unittest("cg-invalid-input", test_cg_invalid_input), &
      & new_unittest("cg-not-converged", test_cg_not_converged), &
      & new_unittest("cg-breakdown", test_cg_breakdown), &
      & new_unittest("cg-verbosity", test_cg_verbosity), &
      & new_unittest("direct-inverse", test_direct_inverse), &
      & new_unittest("direct-singular", test_direct_singular), &
      & new_unittest("direct-verbosity", test_direct_verbosity), &
      & new_unittest("solver-factory", test_solver_factory), &
      & new_unittest("solver-factory-unknown", test_solver_factory_unknown) &
      & ]

end subroutine collect_solver


!> Check that a solution satisfies the linear system
subroutine check_residual(error, amat, xvec, vrhs, message)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in) :: amat(:, :)

   !> Right-hand side vector
   real(wp), intent(in) :: xvec(:)

   !> Solution vector
   real(wp), intent(in) :: vrhs(:)

   !> Failure message
   character(len=*), intent(in) :: message

   real(wp), allocatable :: res(:)

   res = matmul(amat, vrhs) - xvec
   if (any(abs(res) > thr)) then
      call test_failed(error, message)
      print"(a)", "Residual:"
      print"(3es21.14)", res
   end if

end subroutine check_residual


!> Test: Identity matrix 2x2
subroutine test_cg_identity_2x2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 2
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n)

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp), parameter :: tol = 1.0e-15_wp
   integer, parameter :: maxiter = 1000
   integer, parameter :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! Identity matrix
   amat = 0.0_wp
   amat(1,1) = 1.0_wp
   amat(2,2) = 1.0_wp

   ! RHS vector
   xvec = [1.0_wp, 2.0_wp]

   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp]

   ! Solve iteratively
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return
   expected = [0.0_wp, 0.0_wp]
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)

   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for identity matrix")
      print"(a)", "Solution:"
      print"(3es21.14)", vrhs
      print"(a)", "Expected:"
      print"(3es21.14)", expected
   end if

end subroutine test_cg_identity_2x2


!> Test: Diagonal matrix 5x5
subroutine test_cg_diagonal_5x5(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 5
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n), diag(n)
   integer :: i

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp), parameter :: tol = 1.0e-15_wp
   integer, parameter :: maxiter = 1000
   integer, parameter :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! Diagonal matrix with increasing values
   amat = 0.0_wp
   diag = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]
   do i = 1, n
      amat(i,i) = diag(i)
   end do

   ! RHS vector
   xvec = [1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp, 1.0_wp]

   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return
   expected = [0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp, 0.0_wp]
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)

   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for diagonal matrix")
      print"(a)", "Solution:"
      print"(3es21.14)", vrhs
      print"(a)", "Expected:"
      print"(3es21.14)", expected
   end if

end subroutine test_cg_diagonal_5x5


!> Test: Small SPD matrix
subroutine test_cg_spd_small(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 3
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n)

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp), parameter :: tol = 1.0e-15_wp
   integer, parameter :: maxiter = 1000
   integer, parameter :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! SPD matrix
   amat = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
   &1.0_wp, 3.0_wp, 2.0_wp, &
   &1.0_wp, 2.0_wp, 4.0_wp], [n, n])

   ! RHS vector
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]

   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp, 0.0_wp]

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return
   expected = [1.0_wp, 1.0_wp, 1.0_wp]
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)

   ! Check solution
   if (any(abs(vrhs - expected) > thr)) then
      call test_failed(error, "CG solver failed for small SPD matrix")
      print"(a)", "Solution:"
      print"(3es21.14)", vrhs
      print"(a)", "Expected:"
      print"(3es21.14)", expected
   end if

end subroutine test_cg_spd_small


!> Test: Large SPD matrix (1000x1000)
subroutine test_cg_spd_large(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 1000
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:)
   real(wp), allocatable :: expected(:), b(:)
   integer :: i, j

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp), parameter :: tol = 1.0e-15_wp
   integer, parameter :: maxiter = 1000
   integer, parameter :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   allocate(amat(n,n), xvec(n), vrhs(n), expected(n), b(n))

   ! Create a diagonally dominant SPD matrix
   amat = 0.0_wp
   do i = 1, n
      amat(i,i) = real(n, wp) + real(i, wp)
      do j = 1, n
         if (i /= j) then
            amat(i,j) = 1.0_wp / (abs(i-j) + 1.0_wp)
         end if
      end do
   end do

   ! Create a solution vector
   expected = [(sin(real(i, wp) * 0.1_wp), i=1, n)]

   ! Compute RHS: b = A * expected
   b = 0.0_wp
   do i = 1, n
      do j = 1, n
         b(i) = b(i) + amat(i,j) * expected(j)
      end do
   end do

   ! Initial guess
   vrhs = 0.0_wp
   xvec = b

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! Reference solution
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)

   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print"(a)", "Solution:"
      print"(3es21.14)", vrhs
      print"(a)", "Expected:"
      print"(3es21.14)", expected
   end if

end subroutine test_cg_spd_large


!> Test: Ill-conditioned matrix
subroutine test_cg_ill_conditioned(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 8
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n), b(n)
   integer :: i

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp), parameter :: tol = 1.0e-15_wp
   integer, parameter :: maxiter = 1000
   integer, parameter :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! Create an ill-conditioned diagonal matrix
   amat = 0.0_wp
   do i = 1, n
      amat(i,i) = 10.0_wp ** (-(i-1))
   end do

   ! Solution vector
   expected = [(1.0_wp, i=1, n)]

   ! Compute RHS
   b = 0.0_wp
   do i = 1, n
      b(i) = amat(i,i) * expected(i)
   end do

   ! Initial guess
   vrhs = 0.0_wp
   xvec = b

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! Reference solution
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)

   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print"(a)", "Solution:"
      print"(3es21.14)", vrhs
      print"(a)", "Expected:"
      print"(3es21.14)", expected
   end if

end subroutine test_cg_ill_conditioned


!> Test: Zero RHS vector
subroutine test_cg_zero_rhs(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 5
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   real(wp) :: expected(n)
   integer :: i

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp), parameter :: tol = 1.0e-15_wp
   integer, parameter :: maxiter = 1000
   integer, parameter :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! Create a simple SPD matrix
   amat = 0.0_wp
   do i = 1, n
      amat(i,i) = real(i, wp)
      if (i > 1) amat(i,i-1) = 0.1_wp
      if (i < n) amat(i,i+1) = 0.1_wp
   end do

   ! Zero RHS
   xvec = 0.0_wp

   ! Initial guess (non-zero to test convergence)
   vrhs = [1.0_wp, 2.0_wp, 3.0_wp, 4.0_wp, 5.0_wp]

   ! Expected solution (all zeros)
   expected = 0.0_wp

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! Reference solution
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)

   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print"(a)", "Solution:"
      print"(3es21.14)", vrhs
      print"(a)", "Expected:"
      print"(3es21.14)", expected
   end if

end subroutine test_cg_zero_rhs


!> Test: Random SPD matrix
subroutine test_cg_random_spd(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 100
   real(wp), allocatable :: amat(:,:), xvec(:), vrhs(:)
   real(wp), allocatable :: expected(:), b(:), temp(:,:)
   integer :: i, j, k, seed_size
   integer, allocatable :: seed(:)
   real(wp) :: r

   ! Solver variables
   class(mchrg_solver_type), allocatable :: solver
   class(mchrg_solver_input), allocatable :: solver_input
   real(wp), parameter :: tol = 1.0e-15_wp
   integer, parameter :: maxiter = 1000
   integer, parameter :: verbosity = 0

   allocate(cg_input :: solver_input)
   select type (solver_input)
   type is (cg_input)
      solver_input%cgtol = tol
      solver_input%cgmiter = maxiter
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   allocate(amat(n,n), xvec(n), vrhs(n), expected(n), b(n), temp(n,n))

   ! Initialize random seed
   call random_seed(size=seed_size)
   allocate(seed(seed_size))
   seed = 12345
   call random_seed(put=seed)

   ! Generate random matrix B
   temp = 0.0_wp
   do i = 1, n
      do j = 1, i
         call random_number(r)
         ! Random values in [-1, 1]
         temp(i,j) = 2.0_wp * r - 1.0_wp
         if (i /= j) then
            ! Make symmetric
            temp(j,i) = temp(i,j)
         end if
      end do
   end do

   ! Create SPD matrix: A = B^T * B + n*I (ensures positive definiteness)
   amat = 0.0_wp
   do i = 1, n
      do j = 1, n
         do k = 1, n
            amat(i,j) = amat(i,j) + temp(k,i) * temp(k,j)
         end do
      end do
      ! Add diagonal dominance
      amat(i,i) = amat(i,i) + real(n, wp)
   end do

   ! Generate random solution vector
   expected = 0.0_wp
   do i = 1, n
      call random_number(r)
      ! Random values in [-1, 1]
      expected(i) = 2.0_wp * r - 1.0_wp
   end do

   ! Compute RHS: b = A * expected
   b = 0.0_wp
   do i = 1, n
      do j = 1, n
         b(i) = b(i) + amat(i,j) * expected(j)
      end do
   end do

   ! Initial guess
   vrhs = 0.0_wp
   xvec = b

   ! Solve
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=error)
   if (allocated(error)) return

   ! Reference direct solver
   deallocate(solver_input)
   deallocate(solver)
   allocate(direct_input :: solver_input)
   select type(solver_input)
   type is (direct_input)
      solver_input%verbosity = verbosity
   end select
   call new_mchrg_solver(solver, solver_input, error)
   if (allocated(error)) return

   ! Reference solution
   call solver%solve(amat=amat, xvec=xvec, vrhs=expected, error=error)

   ! Check solution
   if (any(abs(vrhs - expected) / max(1.0_wp, abs(expected)) > thr_rel)) then
      call test_failed(error, "CG solver failed for medium SPD matrix")
      print"(a)", "Solution:"
      print"(3es21.14)", vrhs
      print"(a)", "Expected:"
      print"(3es21.14)", expected
   end if

end subroutine test_cg_random_spd


!> Test: Invalid CG settings and requests are rejected
subroutine test_cg_invalid_input(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 3
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   type(cg_input) :: input
   type(cg_solver) :: solver
   type(error_type), allocatable :: solve_error

   amat = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 3.0_wp, 2.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp], [n, n])
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]

   input%cgmiter = 0
   call new_cg_solver(solver, input)
   vrhs = 0.0_wp
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG accepted a non-positive iteration limit")
      return
   end if

   deallocate(input%cgmiter)
   input%cgtol = 0.0_wp
   call new_cg_solver(solver, input)
   vrhs = 0.0_wp
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG accepted a non-positive tolerance")
      return
   end if

   deallocate(input%cgtol)
   call new_cg_solver(solver, input)
   vrhs = 0.0_wp
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, ainv=ainv, &
      & error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG returned an inverse matrix")
      return
   end if

end subroutine test_cg_invalid_input


!> Test: CG reports an error when the iteration limit is exhausted
subroutine test_cg_not_converged(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 3
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   type(cg_input) :: input
   type(cg_solver) :: solver
   type(error_type), allocatable :: solve_error

   amat = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 3.0_wp, 2.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp], [n, n])
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]

   input%cgmiter = 1
   call new_cg_solver(solver, input)
   vrhs = 0.0_wp
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG did not report a missing convergence")
   end if

end subroutine test_cg_not_converged


!> Test: CG reports a breakdown for matrices that are not positive definite
subroutine test_cg_breakdown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 2
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   type(cg_input) :: input
   type(cg_solver) :: solver
   type(error_type), allocatable :: solve_error

   call new_cg_solver(solver, input)

   ! Singular matrix with a right-hand side outside of its range
   amat = reshape([1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp], [n, n])
   xvec = [1.0_wp, -1.0_wp]
   vrhs = 0.0_wp
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, &
         & "CG did not report a breakdown for a singular matrix")
      return
   end if

   ! Indefinite matrix with zero curvature along the first search direction
   amat = reshape([1.0_wp, 0.0_wp, &
      & 0.0_wp, -1.0_wp], [n, n])
   xvec = [1.0_wp, 1.0_wp]
   vrhs = 0.0_wp
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, &
         & "CG did not report a breakdown for an indefinite matrix")
   end if

end subroutine test_cg_breakdown


!> Test: CG solves the linear system at every verbosity level
subroutine test_cg_verbosity(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 3
   real(wp) :: amat(n, n), xvec(n), vrhs(n)
   integer :: verbosity, unit
   type(cg_input) :: input
   type(cg_solver) :: solver

   amat = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 3.0_wp, 2.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp], [n, n])
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]

   open(newunit=unit, status="scratch")
   do verbosity = 0, 2
      input%verbosity = verbosity
      call new_cg_solver(solver, input)
      vrhs = 0.0_wp
      call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, new_unit=unit, &
         & error=error)
      if (allocated(error)) exit
      call check_residual(error, amat, xvec, vrhs, &
         & "CG solution does not satisfy the linear system")
      if (allocated(error)) exit
   end do
   close(unit)
   if (allocated(error)) print '("Verbosity: ", i0)', verbosity

end subroutine test_cg_verbosity


!> Test: The direct solver returns the inverse only on request
subroutine test_direct_inverse(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 4
   real(wp) :: amat(n, n), xvec(n), vrhs(n), vref(n), ainv(n, n)
   real(wp) :: unity(n, n)
   integer :: i
   type(direct_input) :: input
   type(direct_solver) :: solver

   ! Constrained system with a zero diagonal element, symmetric indefinite
   amat = reshape([4.0_wp, 1.0_wp, 0.5_wp, 1.0_wp, &
      & 1.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, &
      & 0.5_wp, 2.0_wp, 4.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp], [n, n])
   xvec = [1.0_wp, -2.0_wp, 0.5_wp, 0.0_wp]

   call new_direct_solver(solver, input)
   call solver%solve(amat=amat, xvec=xvec, vrhs=vref, error=error)
   if (allocated(error)) return
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, ainv=ainv, error=error)
   if (allocated(error)) return

   unity = matmul(amat, ainv)
   do i = 1, n
      unity(i, i) = unity(i, i) - 1.0_wp
   end do
   if (any(abs(unity) > thr)) then
      call test_failed(error, "Direct solver did not return the inverse matrix")
      print"(a)", "A * ainv - 1:"
      print"(4es21.14)", unity
      return
   end if

   if (any(abs(vrhs - vref) > thr) .or. any(abs(matmul(amat, vref) - xvec) > thr)) then
      call test_failed(error, "Direct solutions with and without inverse differ")
      print"(a)", "Solution:"
      print"(4es21.14)", vrhs
      print"(a)", "Reference:"
      print"(4es21.14)", vref
   end if

end subroutine test_direct_inverse


!> Test: The direct solver reports a failed factorization of a singular matrix
subroutine test_direct_singular(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 2
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   type(direct_input) :: input
   type(direct_solver) :: solver
   type(error_type), allocatable :: solve_error

   amat = reshape([1.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp], [n, n])
   xvec = [1.0_wp, -1.0_wp]

   call new_direct_solver(solver, input)
   vrhs = 0.0_wp
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "Direct solver accepted a singular matrix")
      return
   end if

   vrhs = 0.0_wp
   call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, ainv=ainv, &
      & error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "Direct solver inverted a singular matrix")
   end if

end subroutine test_direct_singular


!> Test: The direct solver solves the linear system at every verbosity level
subroutine test_direct_verbosity(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 4
   real(wp) :: amat(n, n), xvec(n), vrhs(n), ainv(n, n)
   integer :: verbosity, unit
   type(direct_input) :: input
   type(direct_solver) :: solver

   ! Constrained system with a zero diagonal element, symmetric indefinite
   amat = reshape([4.0_wp, 1.0_wp, 0.5_wp, 1.0_wp, &
      & 1.0_wp, 3.0_wp, 2.0_wp, 1.0_wp, &
      & 0.5_wp, 2.0_wp, 4.0_wp, 1.0_wp, &
      & 1.0_wp, 1.0_wp, 1.0_wp, 0.0_wp], [n, n])
   xvec = [1.0_wp, -2.0_wp, 0.5_wp, 0.0_wp]

   open(newunit=unit, status="scratch")
   do verbosity = 0, 2
      input%verbosity = verbosity
      call new_direct_solver(solver, input)
      vrhs = 0.0_wp
      call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, new_unit=unit, &
         & error=error)
      if (allocated(error)) exit
      call check_residual(error, amat, xvec, vrhs, &
         & "Direct solution does not satisfy the linear system")
      if (allocated(error)) exit

      vrhs = 0.0_wp
      call solver%solve(amat=amat, xvec=xvec, vrhs=vrhs, ainv=ainv, &
         & new_unit=unit, error=error)
      if (allocated(error)) exit
      call check_residual(error, amat, xvec, vrhs, &
         & "Direct solution with inverse does not satisfy the linear system")
      if (allocated(error)) exit
   end do
   close(unit)
   if (allocated(error)) print '("Verbosity: ", i0)', verbosity

end subroutine test_direct_verbosity


!> Test: The solver factory constructs the solver matching its input
subroutine test_solver_factory(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(cg_input) :: cg_inp
   type(direct_input) :: direct_inp
   class(mchrg_solver_type), allocatable :: solver

   call new_mchrg_solver(solver, cg_inp, error)
   if (allocated(error)) return
   select type (solver)
   type is (cg_solver)
      if (.not. solver%need_pos_def) then
         call test_failed(error, &
            & "CG solver does not require a positive-definite matrix")
         return
      end if
   class default
      call test_failed(error, "Factory did not construct a CG solver")
      return
   end select

   call new_mchrg_solver(solver, direct_inp, error)
   if (allocated(error)) return
   select type (solver)
   type is (direct_solver)
      if (solver%need_pos_def) then
         call test_failed(error, &
            & "Direct solver requires a positive-definite matrix")
      end if
   class default
      call test_failed(error, "Factory did not construct a direct solver")
   end select

end subroutine test_solver_factory


!> Test: The solver factory rejects an unknown solver input
subroutine test_solver_factory_unknown(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(unknown_input) :: input
   class(mchrg_solver_type), allocatable :: solver
   type(error_type), allocatable :: factory_error

   call new_mchrg_solver(solver, input, factory_error)
   if (.not. allocated(factory_error)) then
      call test_failed(error, "Factory accepted an unknown solver input")
      return
   end if
   if (allocated(solver)) then
      call test_failed(error, "Factory constructed a solver for an unknown input")
   end if

end subroutine test_solver_factory_unknown


end module test_solver
