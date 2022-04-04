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

module multicharge_lapack
   use mctc_env, only : sp, dp, wp
   use multicharge_blas, only : gemm
   implicit none
   private

   public :: sytrf, sytrs, sytri
   public :: syinv

   interface syinv
      module procedure :: mchrg_syinv
   end interface syinv

   interface sytrf
      module procedure :: mchrg_ssytrf
      module procedure :: mchrg_dsytrf
   end interface sytrf

   interface sytrs
      module procedure :: mchrg_ssytrs
      module procedure :: mchrg_ssytrs1
      module procedure :: mchrg_ssytrs3
      module procedure :: mchrg_dsytrs
      module procedure :: mchrg_dsytrs1
      module procedure :: mchrg_dsytrs3
   end interface sytrs

   interface sytri
      module procedure :: mchrg_ssytri
      module procedure :: mchrg_dsytri
   end interface sytri


   interface lapack_sytrf
      pure subroutine ssytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine ssytrf
      pure subroutine dsytrf(uplo, n, a, lda, ipiv, work, lwork, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: ipiv(*)
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(inout) :: work(*)
         integer, intent(in) :: lwork
      end subroutine dsytrf
   end interface lapack_sytrf

   interface lapack_sytrs
      pure subroutine ssytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: sp
         real(sp), intent(in) :: a(lda, *)
         real(sp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine ssytrs
      pure subroutine dsytrs(uplo, n, nrhs, a, lda, ipiv, b, ldb, info)
         import :: dp
         real(dp), intent(in) :: a(lda, *)
         real(dp), intent(inout) :: b(ldb, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: nrhs
         integer, intent(in) :: lda
         integer, intent(in) :: ldb
      end subroutine dsytrs
   end interface lapack_sytrs

   interface lapack_sytri
      pure subroutine ssytri(uplo, n, a, lda, ipiv, work, info)
         import :: sp
         real(sp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(sp), intent(in) :: work(*)
      end subroutine ssytri
      pure subroutine dsytri(uplo, n, a, lda, ipiv, work, info)
         import :: dp
         real(dp), intent(inout) :: a(lda, *)
         integer, intent(in) :: ipiv(*)
         character(len=1), intent(in) :: uplo
         integer, intent(out) :: info
         integer, intent(in) :: n
         integer, intent(in) :: lda
         real(dp), intent(in) :: work(*)
      end subroutine dsytri
   end interface lapack_sytri


contains


subroutine mchrg_ssytrf(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, lda, lwork, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   real(sp) :: test(1)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   lwork = -1
   call lapack_sytrf(ula, n, amat, lda, ipiv, test, lwork, stat)
   if (stat == 0) then
      lwork = nint(test(1))
      if (stat_alloc==0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc==0) then
         call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, stat)
      else
         stat = -1000
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytrf failed"
   end if
end subroutine mchrg_ssytrf


subroutine mchrg_dsytrf(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(out) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, lda, lwork, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   real(dp) :: test(1)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   lwork = -1
   call lapack_sytrf(ula, n, amat, lda, ipiv, test, lwork, stat)
   if (stat == 0) then
      lwork = nint(test(1))
      if (stat_alloc==0) then
         allocate(work(lwork), stat=stat_alloc)
      end if
      if (stat_alloc==0) then
         call lapack_sytrf(ula, n, amat, lda, ipiv, work, lwork, stat)
      else
         stat = -1000
      end if
      deallocate(work, stat=stat_dealloc)
   end if
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytrf failed"
   end if
end subroutine mchrg_dsytrf


subroutine mchrg_ssytrs(amat, bmat, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, nrhs, lda, ldb
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytrs failed"
   end if
end subroutine mchrg_ssytrs


subroutine mchrg_dsytrs(amat, bmat, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout) :: bmat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, nrhs, lda, ldb
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   ldb = max(1, size(bmat, 1))
   n = size(amat, 2)
   nrhs = size(bmat, 2)
   call lapack_sytrs(ula, n, nrhs, amat, lda, ipiv, bmat, ldb, stat)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytrs failed"
   end if
end subroutine mchrg_dsytrs


subroutine mchrg_ssytrs1(amat, bvec, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_ssytrs1


subroutine mchrg_ssytrs3(amat, bmat, ipiv, uplo, info)
   real(sp), intent(in) :: amat(:, :)
   real(sp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(sp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_ssytrs3


subroutine mchrg_dsytrs1(amat, bvec, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), target :: bvec(:)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_dsytrs1


subroutine mchrg_dsytrs3(amat, bmat, ipiv, uplo, info)
   real(dp), intent(in) :: amat(:, :)
   real(dp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   real(dp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_dsytrs3


subroutine mchrg_ssytri(amat, ipiv, uplo, info)
   real(sp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, lda, stat_alloc, stat_dealloc
   real(sp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc==0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, stat)
   else
      stat = -1000
   end if
   deallocate(work, stat=stat_dealloc)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] ssytri failed"
   end if
end subroutine mchrg_ssytri


subroutine mchrg_dsytri(amat, ipiv, uplo, info)
   real(dp), intent(inout) :: amat(:, :)
   integer, intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer, intent(out), optional :: info
   character(len=1) :: ula
   integer :: stat, n, lda, stat_alloc, stat_dealloc
   real(dp), allocatable :: work(:)
   if (present(uplo)) then
      ula = uplo
   else
      ula = 'u'
   end if
   lda = max(1, size(amat, 1))
   n = size(amat, 2)
   stat_alloc = 0
   allocate(work(n), stat=stat_alloc)
   if (stat_alloc==0) then
      call lapack_sytri(ula, n, amat, lda, ipiv, work, stat)
   else
      stat = -1000
   end if
   deallocate(work, stat=stat_dealloc)
   if (present(info)) then
      info = stat
   else
      if (stat /= 0) error stop "[multicharge_lapack] dsytri failed"
   end if
end subroutine mchrg_dsytri


subroutine mchrg_syinv(amat)
   real(wp), intent(inout) :: amat(:, :)

   real(wp), allocatable :: rmat(:, :), r2mat(:, :)
   real(wp), parameter :: thr = sqrt(epsilon(1.0_wp)) * 10.0_wp

   allocate(rmat(size(amat, 1), size(amat, 2)))
   allocate(r2mat(size(amat, 1), size(amat, 2)))
   rmat(:, :) = transpose(amat) / &
      (maxval(sum(abs(amat), 1))*maxval(sum(abs(amat), 2)))

   do
      call gemm(rmat, rmat, r2mat)
      call gemm(r2mat, amat, rmat, beta=2.0_wp, alpha=-1.0_wp)
      call set_eye(r2mat)
      call gemm(amat, rmat, r2mat, beta=-1.0_wp)
      if (all(abs(r2mat) < thr)) exit
   end do
   amat(:, :) = rmat
contains
   subroutine set_eye(mat)
      real(wp), intent(out) :: mat(:, :)
      integer :: i, j
      do concurrent(i = 1:size(mat, 1), j = 1:size(mat, 2))
         mat(i, j) = merge(1.0_wp, 0.0_wp, i == j)
      end do
   end subroutine set_eye
end subroutine mchrg_syinv


end module multicharge_lapack
