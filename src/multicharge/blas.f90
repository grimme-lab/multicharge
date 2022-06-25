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

!> Interface to BLAS library for matrix-vector and matrix-matrix operations
module multicharge_blas
   use blas, only : symv => blas_symv, gemv => blas_gemv, gemm => blas_gemm
   use mctc_env, only : wp
   implicit none
   private

   public :: symv, gemv, gemm


   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface gemv
      module procedure :: mchrg_gemv312
      module procedure :: mchrg_gemv321
   end interface

   !> Performs one of the matrix-matrix operations
   !>
   !>    C := alpha*A*B + beta*C,
   !>
   !> or
   !>
   !>    C := alpha*B*A + beta*C,
   !>
   !> where alpha and beta are scalars,  A is a symmetric matrix and  B and
   !> C are  m by n matrices.
   interface gemm
      module procedure :: mchrg_gemm323
      module procedure :: mchrg_gemm233
      module procedure :: mchrg_gemm332
   end interface

contains


subroutine mchrg_gemv312(amat, xvec, yvec, alpha, beta, trans)
   real(wp), intent(in), contiguous, target :: amat(:, :, :)
   real(wp), intent(in) :: xvec(:)
   real(wp), intent(inout), contiguous, target :: yvec(:, :)
   real(wp), intent(in), optional :: alpha
   real(wp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(wp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine mchrg_gemv312


subroutine mchrg_gemv321(amat, xvec, yvec, alpha, beta, trans)
   real(wp), intent(in), contiguous, target :: amat(:, :, :)
   real(wp), intent(in), contiguous, target :: xvec(:, :)
   real(wp), intent(inout) :: yvec(:)
   real(wp), intent(in), optional :: alpha
   real(wp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(wp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine mchrg_gemv321


subroutine mchrg_gemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(wp), intent(in), contiguous, target :: amat(:, :, :)
   real(wp), intent(in) :: bmat(:, :)
   real(wp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(wp), intent(in), optional :: alpha
   real(wp), intent(in), optional :: beta
   real(wp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine mchrg_gemm323


subroutine mchrg_gemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(wp), intent(in) :: amat(:, :)
   real(wp), intent(in), contiguous, target :: bmat(:, :, :)
   real(wp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(wp), intent(in), optional :: alpha
   real(wp), intent(in), optional :: beta
   real(wp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   cptr(1:size(cmat, 1), 1:size(cmat, 2)*size(cmat, 3)) => cmat
   call gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine mchrg_gemm233


subroutine mchrg_gemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(wp), intent(in), contiguous, target :: amat(:, :, :)
   real(wp), intent(in), contiguous, target :: bmat(:, :, :)
   real(wp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(wp), intent(in), optional :: alpha
   real(wp), intent(in), optional :: beta
   real(wp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = 'n'
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = 'n'
   end if
   if (any(tra == ['n', 'N'])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ['n', 'N'])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine mchrg_gemm332


end module multicharge_blas
