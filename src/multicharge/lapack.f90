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
   use lapack_kinds, only : ik
   use lapack, only : sytrf => lapack_sytrf, sytrs => lapack_sytrs, sytri => lapack_sytri
   use mctc_env, only : wp
   implicit none
   private

   public :: sytrf, sytrs, sytri, ik

   interface sytrs
      module procedure :: mchrg_sytrs1
      module procedure :: mchrg_sytrs3
   end interface


contains


subroutine mchrg_sytrs1(amat, bvec, ipiv, uplo, info)
   real(wp), intent(in) :: amat(:, :)
   real(wp), intent(inout), target :: bvec(:)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   real(wp), pointer :: bptr(:, :)
   bptr(1:size(bvec), 1:1) => bvec
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_sytrs1


subroutine mchrg_sytrs3(amat, bmat, ipiv, uplo, info)
   real(wp), intent(in) :: amat(:, :)
   real(wp), intent(inout), contiguous, target :: bmat(:, :, :)
   integer(ik), intent(in) :: ipiv(:)
   character(len=1), intent(in), optional :: uplo
   integer(ik), intent(out), optional :: info
   real(wp), pointer :: bptr(:, :)
   bptr(1:size(bmat, 1), 1:size(bmat, 2)*size(bmat, 3)) => bmat
   call sytrs(amat, bptr, ipiv, uplo, info)
end subroutine mchrg_sytrs3


end module multicharge_lapack
