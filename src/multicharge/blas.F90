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

#ifndef IK
#define IK i4
#endif

!> Interface to BLAS library for vector, matrix-vector and matrix-matrix operations
module multicharge_blas
   use mctc_env, only : sp, dp, ik => IK
   implicit none
   private

   public :: dot, scal, axpy, symv, gemv, gemm


   !> Performs the vector operation
   !>
   !>    y := alpha*x + y,
   !>
   !> where alpha is a scalar and x and y are vectors.
   interface axpy
      module procedure :: mchrg_saxpy
      module procedure :: mchrg_caxpy
      module procedure :: mchrg_daxpy
      module procedure :: mchrg_zaxpy
   end interface axpy

   !> Performs the vector operation
   !>
   !>    x := alpha*x,
   !>
   !> where alpha is a scalar and x is a vector.
   interface scal
      module procedure :: mchrg_sscal
      module procedure :: mchrg_cscal
      module procedure :: mchrg_csscal
      module procedure :: mchrg_dscal
      module procedure :: mchrg_zscal
      module procedure :: mchrg_zdscal
   end interface scal

   !> Forms the dot product of two vectors
   !>
   !>    x**T*y,
   !>
   !> where arrays of higher rank are treated as vectors.
   interface dot
      module procedure :: mchrg_sdot
      module procedure :: mchrg_ddot
      module procedure :: mchrg_sdot12
      module procedure :: mchrg_sdot21
      module procedure :: mchrg_sdot22
      module procedure :: mchrg_ddot12
      module procedure :: mchrg_ddot21
      module procedure :: mchrg_ddot22
   end interface dot

   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface gemv
      module procedure :: mchrg_sgemv
      module procedure :: mchrg_dgemv
      module procedure :: mchrg_sgemv312
      module procedure :: mchrg_sgemv321
      module procedure :: mchrg_dgemv312
      module procedure :: mchrg_dgemv321
   end interface gemv

   !> Performs the matrix-vector  operation
   !>
   !>    y := alpha*A*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are n element vectors and
   !> A is an n by n symmetric matrix.
   interface symv
      module procedure :: mchrg_ssymv
      module procedure :: mchrg_dsymv
   end interface symv

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
      module procedure :: mchrg_sgemm
      module procedure :: mchrg_dgemm
      module procedure :: mchrg_sgemm323
      module procedure :: mchrg_sgemm233
      module procedure :: mchrg_sgemm332
      module procedure :: mchrg_dgemm323
      module procedure :: mchrg_dgemm233
      module procedure :: mchrg_dgemm332
   end interface gemm


   !> Performs the vector operation
   !>
   !>    y := alpha*x + y,
   !>
   !> where alpha is a scalar and x and y are vectors.
   interface blas_axpy
      pure subroutine saxpy(n, alpha, x, incx, y, incy)
         import :: sp, ik
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: n
      end subroutine saxpy
      pure subroutine daxpy(n, alpha, x, incx, y, incy)
         import :: dp, ik
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: n
      end subroutine daxpy
      pure subroutine caxpy(n, alpha, x, incx, y, incy)
         import :: sp, ik
         complex(sp), intent(in) :: x(*)
         complex(sp), intent(inout) :: y(*)
         complex(sp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: n
      end subroutine caxpy
      pure subroutine zaxpy(n, alpha, x, incx, y, incy)
         import :: dp, ik
         complex(dp), intent(in) :: x(*)
         complex(dp), intent(inout) :: y(*)
         complex(dp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: n
      end subroutine zaxpy
   end interface blas_axpy

   !> Performs the vector operation
   !>
   !>    x := alpha*x,
   !>
   !> where alpha is a scalar and x is a vector.
   interface blas_scal
      pure subroutine sscal(n, alpha, x, incx)
         import :: sp, ik
         real(sp), intent(inout) :: x(*)
         real(sp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: n
      end subroutine sscal
      pure subroutine dscal(n, alpha, x, incx)
         import :: dp, ik
         real(dp), intent(inout) :: x(*)
         real(dp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: n
      end subroutine dscal
      pure subroutine cscal(n, alpha, x, incx)
         import :: sp, ik
         complex(sp), intent(inout) :: x(*)
         complex(sp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: n
      end subroutine cscal
      pure subroutine zscal(n, alpha, x, incx)
         import :: dp, ik
         complex(dp), intent(inout) :: x(*)
         complex(dp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: n
      end subroutine zscal
      pure subroutine csscal(n, alpha, x, incx)
         import :: sp, ik
         complex(sp), intent(inout) :: x(*)
         real(sp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: n
      end subroutine csscal
      pure subroutine zdscal(n, alpha, x, incx)
         import :: dp, ik
         complex(dp), intent(inout) :: x(*)
         real(dp), intent(in) :: alpha
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: n
      end subroutine zdscal
   end interface blas_scal

   !> Forms the dot product of two vectors
   !>
   !>    x**T*y,
   !>
   !> where x and y are vectors.
   interface blas_dot
      pure function sdot(n, x, incx, y, incy)
         import :: sp, ik
         real(sp) :: sdot
         real(sp), intent(in) :: x(*)
         real(sp), intent(in) :: y(*)
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: n
      end function sdot
      pure function ddot(n, x, incx, y, incy)
         import :: dp, ik
         real(dp) :: ddot
         real(dp), intent(in) :: x(*)
         real(dp), intent(in) :: y(*)
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: n
      end function ddot
   end interface blas_dot

   !> Performs one of the matrix-vector operations
   !>
   !>    y := alpha*A*x + beta*y,   or   y := alpha*A**T*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are vectors and A is an
   !> m by n matrix.
   interface blas_gemv
      pure subroutine sgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: m
         integer(ik), intent(in) :: n
      end subroutine sgemv
      pure subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         character(len=1), intent(in) :: trans
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: m
         integer(ik), intent(in) :: n
      end subroutine dgemv
   end interface blas_gemv

   !> Performs the matrix-vector  operation
   !>
   !>    y := alpha*A*x + beta*y,
   !>
   !> where alpha and beta are scalars, x and y are n element vectors and
   !> A is an n by n symmetric matrix.
   interface blas_symv
      pure subroutine ssymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: x(*)
         real(sp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: n
      end subroutine ssymv
      pure subroutine dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: x(*)
         real(dp), intent(inout) :: y(*)
         character(len=1), intent(in) :: uplo
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer(ik), intent(in) :: incx
         integer(ik), intent(in) :: incy
         integer(ik), intent(in) :: n
      end subroutine dsymv
   end interface blas_symv

   !> Performs one of the matrix-matrix operations
   !>
   !>    C := alpha*op( A )*op( B ) + beta*C,
   !>
   !> where  op( X ) is one of
   !>
   !>    op( X ) = X   or   op( X ) = X**T,
   !>
   !> alpha and beta are scalars, and A, B and C are matrices, with op( A )
   !> an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
   interface blas_gemm
      pure subroutine sgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: sp, ik
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         integer(ik), intent(in) :: ldc
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(in) :: b(ldb, *)
         real(sp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(sp), intent(in) :: alpha
         real(sp), intent(in) :: beta
         integer(ik), intent(in) :: m
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: k
      end subroutine sgemm
      pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, &
            & beta, c, ldc)
         import :: dp, ik
         integer(ik), intent(in) :: lda
         integer(ik), intent(in) :: ldb
         integer(ik), intent(in) :: ldc
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(in) :: b(ldb, *)
         real(dp), intent(inout) :: c(ldc, *)
         character(len=1), intent(in) :: transa
         character(len=1), intent(in) :: transb
         real(dp), intent(in) :: alpha
         real(dp), intent(in) :: beta
         integer(ik), intent(in) :: m
         integer(ik), intent(in) :: n
         integer(ik), intent(in) :: k
      end subroutine dgemm
   end interface blas_gemm


contains


pure subroutine mchrg_saxpy(xvec, yvec, alpha)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp) :: a
   integer(ik) :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_axpy(n, a, xvec, incx, yvec, incy)
end subroutine mchrg_saxpy


pure subroutine mchrg_caxpy(xvec, yvec, alpha)
   complex(sp), intent(in) :: xvec(:)
   complex(sp), intent(inout) :: yvec(:)
   complex(sp), intent(in), optional :: alpha
   complex(sp) :: a
   integer(ik) :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_axpy(n, a, xvec, incx, yvec, incy)
end subroutine mchrg_caxpy


pure subroutine mchrg_daxpy(xvec, yvec, alpha)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp) :: a
   integer(ik) :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_axpy(n, a, xvec, incx, yvec, incy)
end subroutine mchrg_daxpy


pure subroutine mchrg_zaxpy(xvec, yvec, alpha)
   complex(dp), intent(in) :: xvec(:)
   complex(dp), intent(inout) :: yvec(:)
   complex(dp), intent(in), optional :: alpha
   complex(dp) :: a
   integer(ik) :: incx, incy, n
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   incx = 1
   incy = 1
   n = size(xvec)
   call blas_axpy(n, a, xvec, incx, yvec, incy)
end subroutine mchrg_zaxpy


pure subroutine mchrg_sscal(alpha, xvec)
   real(sp), intent(in) :: alpha
   real(sp), intent(inout) :: xvec(:)
   integer(ik) :: incx, n
   incx = 1
   n = size(xvec)
   call blas_scal(n, alpha, xvec, incx)
end subroutine mchrg_sscal


pure subroutine mchrg_cscal(alpha, xvec)
   complex(sp), intent(in) :: alpha
   complex(sp), intent(inout) :: xvec(:)
   integer(ik) :: incx, n
   incx = 1
   n = size(xvec)
   call blas_scal(n, alpha, xvec, incx)
end subroutine mchrg_cscal


pure subroutine mchrg_csscal(alpha, xvec)
   real(sp), intent(in) :: alpha
   complex(sp), intent(inout) :: xvec(:)
   integer(ik) :: incx, n
   incx = 1
   n = size(xvec)
   call blas_scal(n, alpha, xvec, incx)
end subroutine mchrg_csscal


pure subroutine mchrg_dscal(alpha, xvec)
   real(dp), intent(in) :: alpha
   real(dp), intent(inout) :: xvec(:)
   integer(ik) :: incx, n
   incx = 1
   n = size(xvec)
   call blas_scal(n, alpha, xvec, incx)
end subroutine mchrg_dscal


pure subroutine mchrg_zscal(alpha, xvec)
   complex(dp), intent(in) :: alpha
   complex(dp), intent(inout) :: xvec(:)
   integer(ik) :: incx, n
   incx = 1
   n = size(xvec)
   call blas_scal(n, alpha, xvec, incx)
end subroutine mchrg_zscal


pure subroutine mchrg_zdscal(alpha, xvec)
   real(dp), intent(in) :: alpha
   complex(dp), intent(inout) :: xvec(:)
   integer(ik) :: incx, n
   incx = 1
   n = size(xvec)
   call blas_scal(n, alpha, xvec, incx)
end subroutine mchrg_zdscal


function mchrg_sdot(xvec, yvec) result(dot_prod)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in) :: yvec(:)
   real(sp) :: dot_prod
   integer(ik) :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   dot_prod = blas_dot(n, xvec, incx, yvec, incy)
end function mchrg_sdot


function mchrg_ddot(xvec, yvec) result(dot_prod)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in) :: yvec(:)
   real(dp) :: dot_prod
   integer(ik) :: incx, incy, n
   incx = 1
   incy = 1
   n = size(xvec)
   dot_prod = blas_dot(n, xvec, incx, yvec, incy)
end function mchrg_ddot


function mchrg_sdot12(xvec, yvec) result(dot_prod)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(in), contiguous, target :: yvec(:, :)
   real(sp) :: dot_prod
   real(sp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   dot_prod = dot(xvec, yptr)
end function mchrg_sdot12


function mchrg_sdot21(xvec, yvec) result(dot_prod)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(in) :: yvec(:)
   real(sp) :: dot_prod
   real(sp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   dot_prod = dot(xptr, yvec)
end function mchrg_sdot21


function mchrg_sdot22(xvec, yvec) result(dot_prod)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(in), contiguous, target :: yvec(:, :)
   real(sp) :: dot_prod
   real(sp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   dot_prod = dot(xptr, yptr)
end function mchrg_sdot22


function mchrg_ddot12(xvec, yvec) result(dot_prod)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(in), contiguous, target :: yvec(:, :)
   real(dp) :: dot_prod
   real(dp), pointer :: yptr(:)
   yptr(1:size(yvec)) => yvec
   dot_prod = dot(xvec, yptr)
end function mchrg_ddot12


function mchrg_ddot21(xvec, yvec) result(dot_prod)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(in) :: yvec(:)
   real(dp) :: dot_prod
   real(dp), pointer :: xptr(:)
   xptr(1:size(xvec)) => xvec
   dot_prod = dot(xptr, yvec)
end function mchrg_ddot21


function mchrg_ddot22(xvec, yvec) result(dot_prod)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(in), contiguous, target :: yvec(:, :)
   real(dp) :: dot_prod
   real(dp), pointer :: xptr(:), yptr(:)
   xptr(1:size(xvec)) => xvec
   yptr(1:size(yvec)) => yvec
   dot_prod = dot(xptr, yptr)
end function mchrg_ddot22


subroutine mchrg_sgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout), contiguous, target :: yvec(:, :)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = "n"
   end if
   if (any(tra == ["n", "N"])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine mchrg_sgemv312


subroutine mchrg_sgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: xvec(:, :)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = "n"
   end if
   if (any(tra == ["n", "N"])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine mchrg_sgemv321


subroutine mchrg_dgemv312(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout), contiguous, target :: yvec(:, :)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), yptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = "n"
   end if
   if (any(tra == ["n", "N"])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      yptr(1:size(yvec, 1)*size(yvec, 2)) => yvec
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      yptr(1:size(yvec, 1) * size(yvec, 2)) => yvec
   end if
   call gemv(aptr, xvec, yptr, alpha, beta, tra)
end subroutine mchrg_dgemv312


subroutine mchrg_dgemv321(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: xvec(:, :)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp), pointer :: aptr(:, :), xptr(:)
   character(len=1) :: tra
   if (present(trans)) then
      tra = trans
   else
      tra = "n"
   end if
   if (any(tra == ["n", "N"])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
      xptr(1:size(xvec, 1)*size(xvec, 2)) => xvec
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
      xptr(1:size(xvec, 1) * size(xvec, 2)) => xvec
   end if
   call gemv(aptr, xptr, yvec, alpha, beta, tra)
end subroutine mchrg_dgemv321


pure subroutine mchrg_sgemv(amat, xvec, yvec, alpha, beta, trans)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(sp) :: a, b
   character(len=1) :: tra
   integer(ik) :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = "n"
   end if
   incx = 1_ik
   incy = 1_ik
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine mchrg_sgemv


pure subroutine mchrg_dgemv(amat, xvec, yvec, alpha, beta, trans)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1), intent(in), optional :: trans
   real(dp) :: a, b
   character(len=1) :: tra
   integer(ik) :: incx, incy, m, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(trans)) then
      tra = trans
   else
      tra = "n"
   end if
   incx = 1_ik
   incy = 1_ik
   lda = max(1, size(amat, 1))
   m = size(amat, 1)
   n = size(amat, 2)
   call blas_gemv(tra, m, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine mchrg_dgemv


pure subroutine mchrg_ssymv(amat, xvec, yvec, uplo, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: xvec(:)
   real(sp), intent(inout) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: ula
   real(sp) :: a, b
   integer(ik) :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = "u"
   end if
   incx = 1_ik
   incy = 1_ik
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call blas_symv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine mchrg_ssymv


pure subroutine mchrg_dsymv(amat, xvec, yvec, uplo, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: xvec(:)
   real(dp), intent(inout) :: yvec(:)
   character(len=1), intent(in), optional :: uplo
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: ula
   real(dp) :: a, b
   integer(ik) :: incx, incy, n, lda
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0
   end if
   if (present(uplo)) then
      ula = uplo
   else
      ula = "u"
   end if
   incx = 1_ik
   incy = 1_ik
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   call blas_symv(ula, n, a, amat, lda, xvec, incx, b, yvec, incy)
end subroutine mchrg_dsymv


pure subroutine mchrg_sgemm(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   character(len=1) :: tra, trb
   real(sp) :: a, b
   integer(ik) :: m, n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_sp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_sp
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = "n"
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = "n"
   end if
   if ((tra=="n".or.tra=="N")) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
   call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
end subroutine mchrg_sgemm


pure subroutine mchrg_dgemm(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   character(len=1) :: tra, trb
   real(dp) :: a, b
   integer(ik) :: m, n, k, lda, ldb, ldc
   if (present(alpha)) then
      a = alpha
   else
      a = 1.0_dp
   end if
   if (present(beta)) then
      b = beta
   else
      b = 0.0_dp
   end if
   if (present(transa)) then
      tra = transa
   else
      tra = "n"
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = "n"
   end if
   if ((tra=="n".or.tra=="N")) then
      k = size(amat, 2)
   else
      k = size(amat, 1)
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   ldc = max(1, size(cmat, 1))
   m = size(cmat, 1)
   n = size(cmat, 2)
   call blas_gemm(tra, trb, m, n, k, a, amat, lda, bmat, ldb, b, cmat, ldc)
end subroutine mchrg_dgemm


subroutine mchrg_sgemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in) :: bmat(:, :)
   real(sp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra
   if (present(transa)) then
      tra = transa
   else
      tra = "n"
   end if
   if (any(tra == ["n", "N"])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine mchrg_sgemm323


subroutine mchrg_sgemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(in), contiguous, target :: bmat(:, :, :)
   real(sp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = "n"
   end if
   if (any(trb == ["n", "N"])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   cptr(1:size(cmat, 1), 1:size(cmat, 2)*size(cmat, 3)) => cmat
   call gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine mchrg_sgemm233


subroutine mchrg_sgemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(sp), intent(in), contiguous, target :: amat(:, :, :)
   real(sp), intent(in), contiguous, target :: bmat(:, :, :)
   real(sp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(sp), intent(in), optional :: alpha
   real(sp), intent(in), optional :: beta
   real(sp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = "n"
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = "n"
   end if
   if (any(tra == ["n", "N"])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ["n", "N"])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine mchrg_sgemm332


subroutine mchrg_dgemm323(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in) :: bmat(:, :)
   real(dp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: aptr(:, :), cptr(:, :)
   character(len=1) :: tra
   if (present(transa)) then
      tra = transa
   else
      tra = "n"
   end if
   if (any(tra == ["n", "N"])) then
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   end if
   cptr(1:size(cmat, 1)*size(cmat, 2), 1:size(cmat, 3)) => cmat
   call gemm(aptr, bmat, cptr, tra, transb, alpha, beta)
end subroutine mchrg_dgemm323


subroutine mchrg_dgemm233(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(in), contiguous, target :: bmat(:, :, :)
   real(dp), intent(inout), contiguous, target :: cmat(:, :, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: bptr(:, :), cptr(:, :)
   character(len=1) :: trb
   if (present(transb)) then
      trb = transb
   else
      trb = "n"
   end if
   if (any(trb == ["n", "N"])) then
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   end if
   cptr(1:size(cmat, 1), 1:size(cmat, 2)*size(cmat, 3)) => cmat
   call gemm(amat, bptr, cptr, transa, trb, alpha, beta)
end subroutine mchrg_dgemm233


subroutine mchrg_dgemm332(amat, bmat, cmat, transa, transb, alpha, beta)
   real(dp), intent(in), contiguous, target :: amat(:, :, :)
   real(dp), intent(in), contiguous, target :: bmat(:, :, :)
   real(dp), intent(inout) :: cmat(:, :)
   character(len=1), intent(in), optional :: transa
   character(len=1), intent(in), optional :: transb
   real(dp), intent(in), optional :: alpha
   real(dp), intent(in), optional :: beta
   real(dp), pointer :: aptr(:, :), bptr(:, :)
   character(len=1) :: tra, trb
   if (present(transa)) then
      tra = transa
   else
      tra = "n"
   end if
   if (present(transb)) then
      trb = transb
   else
      trb = "n"
   end if
   if (any(tra == ["n", "N"])) then
      aptr(1:size(amat, 1), 1:size(amat, 2)*size(amat, 3)) => amat
   else
      aptr(1:size(amat, 1)*size(amat, 2), 1:size(amat, 3)) => amat
   end if
   if (any(trb == ["n", "N"])) then
      bptr(1:size(bmat, 1)*size(bmat, 2), 1:size(bmat, 3)) => bmat
   else
      bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   end if
   call gemm(aptr, bptr, cmat, tra, trb, alpha, beta)
end subroutine mchrg_dgemm332


end module multicharge_blas
