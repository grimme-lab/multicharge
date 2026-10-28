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

!> @dir multicharge/solver
!> Contains the implementation of the linear equations solvers

!> @file multicharge/solver.f90
!> Provides a reexport of the solvers implementations

!> Proxy module to reexport the solver implementations
module multicharge_solver
   use mctc_env, only : error_type, fatal_error
   use multicharge_solver_cg, only : cg_solver, new_cg_solver, cg_input
   use multicharge_solver_direct, only : direct_solver, new_direct_solver, &
      & direct_input
   use multicharge_solver_type, only : mchrg_solver_type, mchrg_solver_input
   implicit none
   private

   public :: mchrg_solver_type, mchrg_solver_input
   public :: direct_solver, new_direct_solver, direct_input
   public :: cg_solver, new_cg_solver, cg_input
   public :: new_mchrg_solver


contains


!> Construct a solver instance from its input configuration
subroutine new_mchrg_solver(solver, input, error)

   !> Solver instance to construct
   class(mchrg_solver_type), intent(out), allocatable :: solver

   !> Solver input configuration
   class(mchrg_solver_input), intent(in) :: input

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   select type (input)
   type is (cg_input)
      block
         class(cg_solver), allocatable :: tmp
         allocate(tmp)
         call new_cg_solver(tmp, input)
         call move_alloc(tmp, solver)
      end block
   type is (direct_input)
      block
         class(direct_solver), allocatable :: tmp
         allocate(tmp)
         call new_direct_solver(tmp, input)
         call move_alloc(tmp, solver)
      end block
   class default
      call fatal_error(error, "Unknown solver input type")
      return
   end select

end subroutine new_mchrg_solver


end module multicharge_solver
