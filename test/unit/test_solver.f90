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
   use multicharge_model_cache, only : mchrg_cache
   use multicharge_solver, only : new_mchrg_solver
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, &
      & direct_input
   use multicharge_solver_trialprod, only : trialprod, dense_trialprod
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

   !> Matrix-free shifted 1D Laplacian, tridiagonal with 2 + shift and -1
   type, extends(trialprod) :: laplace_trialprod

      !> Diagonal shift keeping the operator positive definite
      real(wp) :: shift = 0.5_wp
   contains

      !> Apply the shifted Laplacian to a trial vector
      procedure :: multiply => laplace_multiply

   end type laplace_trialprod


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
      & new_unittest("cg-matrix-free", test_cg_matrix_free), &
      & new_unittest("cg-preconditioner", test_cg_preconditioner), &
      & new_unittest("direct-inverse", test_direct_inverse), &
      & new_unittest("direct-singular", test_direct_singular), &
      & new_unittest("direct-verbosity", test_direct_verbosity), &
      & new_unittest("solver-factory", test_solver_factory), &
      & new_unittest("solver-factory-unknown", test_solver_factory_unknown) &
      & ]

end subroutine collect_solver


!> Apply the shifted 1D Laplacian to a trial vector without storing a matrix
subroutine laplace_multiply(self, x, y)

   !> Matrix-free Laplacian operator
   class(laplace_trialprod), intent(inout) :: self

   !> Trial vector
   real(wp), intent(in) :: x(:)

   !> Resulting vector
   real(wp), intent(out) :: y(:)

   integer :: i, n

   n = size(x)
   do i = 1, n
      y(i) = (2.0_wp + self%shift) * x(i)
      if (i > 1) y(i) = y(i) - x(i - 1)
      if (i < n) y(i) = y(i) - x(i + 1)
   end do

end subroutine laplace_multiply


!> Compute the diagonal (Jacobi) preconditioner of a dense matrix
subroutine get_jacobi(amat, prec)

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in) :: amat(:, :)

   !> Inverse diagonal of the coefficient matrix
   real(wp), allocatable, intent(out) :: prec(:)

   integer :: i

   allocate(prec(size(amat, 1)))
   do i = 1, size(amat, 1)
      prec(i) = 1.0_wp / (amat(i, i) + tiny(1.0_wp))
   end do

end subroutine get_jacobi


!> Solve a dense linear system the same way the charge model does: wrap the
!> matrix in a cache and dense trial product, precompute the Jacobi
!> preconditioner and call the solver.
subroutine solve_dense(solver, amat, xvec, vrhs, error, ainv, new_unit)

   !> Solver instance
   class(mchrg_solver_type), intent(in) :: solver

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in) :: amat(:, :)

   !> Right-hand side vector
   real(wp), intent(in) :: xvec(:)

   !> On input: initial guess; on output: solution
   real(wp), intent(inout), contiguous :: vrhs(:)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Inverse matrix, requested by passing the argument
   real(wp), intent(out), optional :: ainv(:, :)

   !> Output unit
   integer, intent(in), optional :: new_unit

   type(mchrg_cache), target :: cache
   type(dense_trialprod) :: tprod
   real(wp), allocatable :: prec(:)

   cache%amat = amat
   tprod%cache => cache
   call get_jacobi(amat, prec)

   call solver%solve(tprod=tprod, xvec=xvec, vrhs=vrhs, prec=prec, &
      & ainv=ainv, new_unit=new_unit, error=error)

end subroutine solve_dense


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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=expected, error=error)

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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=expected, error=error)

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
      & 1.0_wp, 3.0_wp, 2.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp], [n, n])

   ! RHS vector
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]

   ! Initial guess
   vrhs = [0.0_wp, 0.0_wp, 0.0_wp]

   ! Solve
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=expected, error=error)

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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=expected, error=error)

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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=expected, error=error)

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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=expected, error=error)

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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=expected, error=error)

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
   type(mchrg_cache), target :: cache
   type(dense_trialprod) :: tprod
   real(wp), allocatable :: prec(:)

   amat = reshape([4.0_wp, 1.0_wp, 1.0_wp, &
      & 1.0_wp, 3.0_wp, 2.0_wp, &
      & 1.0_wp, 2.0_wp, 4.0_wp], [n, n])
   xvec = [6.0_wp, 6.0_wp, 7.0_wp]

   input%cgmiter = 0
   call new_cg_solver(solver, input)
   vrhs = 0.0_wp
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG accepted a non-positive iteration limit")
      return
   end if

   deallocate(input%cgmiter)
   input%cgtol = 0.0_wp
   call new_cg_solver(solver, input)
   vrhs = 0.0_wp
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG accepted a non-positive tolerance")
      return
   end if

   deallocate(input%cgtol)
   call new_cg_solver(solver, input)
   vrhs = 0.0_wp
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, ainv=ainv, &
      & error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG returned an inverse matrix")
      return
   end if

   ! The preconditioner has to be computed by the caller
   cache%amat = amat
   tprod%cache => cache
   vrhs = 0.0_wp
   call solver%solve(tprod=tprod, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG accepted a missing preconditioner")
      return
   end if

   allocate(prec(n - 1), source=1.0_wp)
   vrhs = 0.0_wp
   call solver%solve(tprod=tprod, xvec=xvec, vrhs=vrhs, prec=prec, &
      & error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "CG accepted a preconditioner of wrong size")
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
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
      call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, new_unit=unit, &
         & error=error)
      if (allocated(error)) exit
      call check_residual(error, amat, xvec, vrhs, &
         & "CG solution does not satisfy the linear system")
      if (allocated(error)) exit
   end do
   close(unit)
   if (allocated(error)) print '("Verbosity: ", i0)', verbosity

end subroutine test_cg_verbosity


!> Test: CG solves a matrix-free system given only by its trial product
subroutine test_cg_matrix_free(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 50
   real(wp), parameter :: shift = 0.5_wp
   real(wp) :: amat(n, n), xvec(n), vrhs(n), vref(n), prec(n)
   integer :: i
   type(cg_input) :: input
   type(cg_solver) :: cg
   type(direct_solver) :: direct
   type(laplace_trialprod) :: tprod
   type(error_type), allocatable :: solve_error

   tprod%shift = shift
   do i = 1, n
      xvec(i) = sin(real(i, wp))
   end do

   ! Diagonal preconditioner computed by the caller
   prec(:) = 1.0_wp / (2.0_wp + shift)

   call new_cg_solver(cg, input)
   vrhs = 0.0_wp
   call cg%solve(tprod=tprod, xvec=xvec, vrhs=vrhs, prec=prec, error=error)
   if (allocated(error)) return

   ! Reference from the direct solver on the equivalent dense matrix
   amat = 0.0_wp
   do i = 1, n
      amat(i, i) = 2.0_wp + shift
      if (i > 1) amat(i, i - 1) = -1.0_wp
      if (i < n) amat(i, i + 1) = -1.0_wp
   end do
   call new_direct_solver(direct, direct_input())
   vref = 0.0_wp
   call solve_dense(direct, amat=amat, xvec=xvec, vrhs=vref, error=error)
   if (allocated(error)) return

   call check_residual(error, amat, xvec, vrhs, &
      & "Matrix-free CG solution does not satisfy the linear system")
   if (allocated(error)) return

   if (any(abs(vrhs - vref) > thr1)) then
      call test_failed(error, "Matrix-free CG differs from the direct solution")
      print"(a)", "Difference:"
      print"(5es21.14)", vrhs - vref
      return
   end if

   ! The direct solver needs an explicitly stored matrix
   vref = 0.0_wp
   call direct%solve(tprod=tprod, xvec=xvec, vrhs=vref, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, &
         & "Direct solver accepted an operator without a dense matrix")
   end if

end subroutine test_cg_matrix_free


!> Test: The caller-provided diagonal preconditioner makes a badly scaled
!> system solvable within an iteration budget where plain CG fails
subroutine test_cg_preconditioner(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   integer, parameter :: n = 100
   integer, parameter :: maxiter = 80
   real(wp), parameter :: tol = 1.0e-12_wp
   real(wp) :: amat(n, n), xvec(n), vrhs(n), vref(n), scale(n)
   real(wp), allocatable :: prec(:)
   integer :: i
   type(cg_input) :: input
   type(cg_solver) :: cg
   type(direct_solver) :: direct
   type(mchrg_cache), target :: cache
   type(dense_trialprod) :: tprod
   type(error_type), allocatable :: solve_error

   ! Well-conditioned tridiagonal matrix with diagonal scaling over six
   ! orders of magnitude, amat = S * B * S
   do i = 1, n
      scale(i) = sqrt(10.0_wp**(-3.0_wp + 6.0_wp * real(i - 1, wp) &
         & / real(n - 1, wp)))
      xvec(i) = cos(real(i, wp))
   end do
   amat = 0.0_wp
   do i = 1, n
      amat(i, i) = 2.5_wp * scale(i)**2
      if (i > 1) amat(i, i - 1) = -scale(i) * scale(i - 1)
      if (i < n) amat(i, i + 1) = -scale(i) * scale(i + 1)
   end do
   cache%amat = amat
   tprod%cache => cache

   ! Reference from the direct solver
   call new_direct_solver(direct, direct_input())
   vref = 0.0_wp
   call solve_dense(direct, amat=amat, xvec=xvec, vrhs=vref, error=error)
   if (allocated(error)) return

   input%cgmiter = maxiter
   input%cgtol = tol
   call new_cg_solver(cg, input)

   ! Diagonal preconditioner computed before calling the solver
   call get_jacobi(amat, prec)
   do i = 1, n
      if (abs(prec(i) * amat(i, i) - 1.0_wp) > thr) then
         call test_failed(error, "Jacobi preconditioner is not the inverse diagonal")
         return
      end if
   end do

   vrhs = 0.0_wp
   call cg%solve(tprod=tprod, xvec=xvec, vrhs=vrhs, prec=prec, error=error)
   if (allocated(error)) return

   if (norm2(vrhs - vref) > thr_rel * norm2(vref)) then
      call test_failed(error, "Preconditioned CG differs from the direct solution")
      print '("Relative error: ", es12.4)', norm2(vrhs - vref) / norm2(vref)
      return
   end if

   ! Without preconditioning the same iteration budget is not sufficient
   prec(:) = 1.0_wp
   vrhs = 0.0_wp
   call cg%solve(tprod=tprod, xvec=xvec, vrhs=vrhs, prec=prec, &
      & error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, &
         & "Unpreconditioned CG converged within the reduced iteration budget")
   end if

end subroutine test_cg_preconditioner


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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vref, error=error)
   if (allocated(error)) return
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, ainv=ainv, error=error)
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
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, error=solve_error)
   if (.not. allocated(solve_error)) then
      call test_failed(error, "Direct solver accepted a singular matrix")
      return
   end if

   vrhs = 0.0_wp
   call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, ainv=ainv, &
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
      call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, new_unit=unit, &
         & error=error)
      if (allocated(error)) exit
      call check_residual(error, amat, xvec, vrhs, &
         & "Direct solution does not satisfy the linear system")
      if (allocated(error)) exit

      vrhs = 0.0_wp
      call solve_dense(solver, amat=amat, xvec=xvec, vrhs=vrhs, ainv=ainv, &
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
