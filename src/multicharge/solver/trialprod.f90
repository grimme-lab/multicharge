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

!> @file multicharge/solver/trialprod.f90
!> Provides a general base class for the multiplication between
!> a linear system matrix and a trial vector.

!> Abstract operator for the linear system matrix-vector product
module multicharge_solver_trialprod
   use mctc_env, only : wp
   use multicharge_blas, only : symv
   use multicharge_model_cache, only : mchrg_cache
   implicit none
   private

   public :: trialprod, dense_trialprod


   !> Abstract product of the linear system matrix with a trial vector
   type, abstract :: trialprod

      !> Cache holding the linear system, referenced but not owned
      type(mchrg_cache), pointer :: cache => null()
   contains

      !> Multiply the linear system matrix with a trial vector
      procedure(multiply), deferred :: multiply

   end type trialprod


   abstract interface
      !> Multiply the linear system matrix with a trial vector
      subroutine multiply(self, x, y)
         import :: trialprod, wp

         !> Matrix-vector multiplication type
         class(trialprod), intent(inout) :: self

         !> Trial vector
         real(wp), intent(in) :: x(:)

         !> Resulting vector
         real(wp), intent(out) :: y(:)

      end subroutine multiply
   end interface


   !> Trial product with the dense coefficient matrix stored in the cache
   type, extends(trialprod) :: dense_trialprod
   contains

      !> Perform dense matrix-vector multiplication
      procedure :: multiply => dense_multiply

   end type dense_trialprod


contains


!> Perform dense matrix-vector multiplication with the cached matrix
subroutine dense_multiply(self, x, y)

   !> Dense trial product operator
   class(dense_trialprod), intent(inout) :: self

   !> Trial vector
   real(wp), intent(in) :: x(:)

   !> Resulting vector
   real(wp), intent(out) :: y(:)

   call symv(self%cache%amat, x, y, alpha=1.0_wp, beta=0.0_wp)

end subroutine dense_multiply


end module multicharge_solver_trialprod
