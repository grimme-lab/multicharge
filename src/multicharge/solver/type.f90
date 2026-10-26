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

!> @file multicharge/solver/type.f90
!> Provides a general base class for the linear system solvers

!> Base types for the linear system solvers
module multicharge_solver_type
   use mctc_env, only : error_type, wp
   implicit none
   private

   public :: mchrg_solver_type, mchrg_solver_input


   !> Abstract base type for multicharge solvers
   type, abstract :: mchrg_solver_type

      !> Whether the solver requires a positive-definite matrix
      logical, allocatable :: need_pos_def
   contains

      !> Solve the linear system
      procedure(solve), deferred :: solve

   end type mchrg_solver_type


   abstract interface
      !> Solve the linear system of equations
      subroutine solve(self, amat, xvec, vrhs, ainv, new_unit, error)
         import :: mchrg_solver_type, error_type, wp

         !> Solver instance
         class(mchrg_solver_type), intent(in) :: self

         !> Dense coefficient matrix of the linear system
         real(wp), intent(in) :: amat(:, :)

         !> Right-hand side vector
         real(wp), intent(in) :: xvec(:)

         !> On input: initial guess; on output: solution
         real(wp), intent(inout), contiguous :: vrhs(:)

         !> Inverse matrix, requested by passing the argument
         real(wp), intent(out), optional :: ainv(:, :)

         !> Output unit
         integer, intent(in), optional :: new_unit

         !> Error handling
         type(error_type), allocatable, intent(out) :: error

      end subroutine solve
   end interface


   !> Abstract base type for solver configuration
   type, abstract :: mchrg_solver_input
   end type mchrg_solver_input


end module multicharge_solver_type
