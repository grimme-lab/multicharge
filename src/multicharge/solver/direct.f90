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

!> @file multicharge/solver/direct.f90
!> Provides implementation of the direct solver using LAPACK
!> for symmetric indefinite systems

!> Direct linear system solver
module multicharge_solver_direct
   use, intrinsic :: iso_fortran_env, only : output_unit
   use mctc_env, only : error_type, fatal_error, timer_type, wp
   use multicharge_blas, only : symv
   use multicharge_lapack, only : sytrf, sytri, sytrs
   use multicharge_output, only : print_direct_final, write_direct_solver
   use multicharge_solver_type, only : mchrg_solver_input, mchrg_solver_type
   implicit none
   private

   public :: direct_solver, direct_input, new_direct_solver


   !> Input configuration for the direct solver
   type, extends(mchrg_solver_input) :: direct_input

      !> Optional output verbosity
      integer, allocatable :: verbosity
   end type direct_input

   !> Direct LAPACK solver for symmetric indefinite systems
   type, extends(mchrg_solver_type) :: direct_solver

      !> Output verbosity
      integer, allocatable :: verbosity
   contains

      !> Solve the linear system directly
      procedure :: solve
   end type direct_solver

   !> Default verbosity level
   integer, parameter :: verbosity_def = 0


contains


!> Construct a direct solver from its input configuration
subroutine new_direct_solver(self, input)

   !> Direct solver instance
   class(direct_solver), intent(out) :: self

   !> Direct solver configuration
   type(direct_input), intent(in) :: input

   self%need_pos_def = .false.
   if (allocated(input%verbosity)) then
      self%verbosity = input%verbosity
   else
      self%verbosity = verbosity_def
   end if

end subroutine new_direct_solver


!> Solve a dense symmetric linear system using LAPACK
subroutine solve(self, amat, xvec, vrhs, ainv, new_unit, error)

   !> Direct solver instance
   class(direct_solver), intent(in) :: self

   !> Dense coefficient matrix of the linear system
   real(wp), intent(in) :: amat(:, :)

   !> Right-hand side vector
   real(wp), intent(in) :: xvec(:)

   !> On input: initial guess; on output: solution
   real(wp), intent(inout), contiguous :: vrhs(:)

   !> Inverse coefficient matrix, requested by passing the argument
   real(wp), intent(out), optional :: ainv(:, :)

   !> Output unit (optional)
   integer, intent(in), optional :: new_unit

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp), allocatable :: invmat(:, :)
   integer, allocatable :: ipiv(:)

   integer :: local_info, ndim, ic, jc, unit
   type(timer_type) :: timer

   if (self%verbosity > 1) call timer%push("total")

   if (present(new_unit)) then
      unit = new_unit
   else
      unit = output_unit
   end if

   if (self%verbosity > 0) then
      call write_direct_solver(unit)
   end if

   ! Dimensions match check
   ndim = size(xvec)

   if (size(amat, 1) /= ndim .or. size(amat, 2) /= ndim &
      & .or. size(vrhs) /= ndim) then
      call fatal_error(error, "solve_direct: dimension mismatch.")
      return
   end if
   if (present(ainv)) then
      if (size(ainv, 1) /= ndim .or. size(ainv, 2) /= ndim) then
         call fatal_error(error, "solve_direct: dimension mismatch of ainv.")
         return
      end if
   end if

   allocate(invmat(ndim, ndim))
   invmat = amat
   vrhs = xvec

   ! Factorize the Coulomb matrix
   allocate(ipiv(ndim))
   call sytrf(invmat, ipiv, info=local_info, uplo="l")
   if (local_info /= 0) then
      call fatal_error(error, "Bunch-Kaufman factorization failed.")
      return
   end if

   if (present(ainv)) then
      ! Inverted matrix is needed for coupled-perturbed equations
      call sytri(invmat, ipiv, info=local_info, uplo="l")
      if (local_info /= 0) then
         call fatal_error(error, "Inversion of factorized matrix failed.")
         return
      end if
      ! Solve the linear system
      call symv(invmat, xvec, vrhs, uplo="l")
      do ic = 1, ndim
         do jc = ic + 1, ndim
            invmat(ic, jc) = invmat(jc, ic)
         end do
      end do
      ainv = invmat
   else
      ! Solve the linear system
      call sytrs(invmat, vrhs, ipiv, info=local_info, uplo="l")
      if (local_info /= 0) then
         call fatal_error(error, "Solution of linear system failed.")
         return
      end if
   end if

   ! Stop the solve timer
   if (self%verbosity > 1) call timer%pop
   call print_direct_final(unit, timer, self%verbosity)

end subroutine solve


end module multicharge_solver_direct
