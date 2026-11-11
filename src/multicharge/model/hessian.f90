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

!> @file multicharge/model/hessian.f90
!> Provides shared kernels for analytical charge-response Hessians

!> Small kernels shared by the charge-response implementations.
module multicharge_model_hessian
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use mctc_io_math, only: outerprod
   use multicharge_blas, only: gemm
   implicit none
   private

   public :: coulomb_kernel, radial_hessian, add_pair_hessian, add_pair_block, &
      & add_chain_hessian, add_symmetric_product


contains


!> Gaussian Coulomb kernel and partial derivatives in distance and exponent
pure subroutine coulomb_kernel(r, gam, val, dr, drr, dg, drg, dgg)

   !> Interatomic distance
   real(wp), intent(in) :: r

   !> Effective Gaussian charge width parameter
   real(wp), intent(in) :: gam

   !> Kernel value
   real(wp), intent(out) :: val

   !> First derivative w.r.t. distance
   real(wp), intent(out) :: dr

   !> Second derivative w.r.t. distance
   real(wp), intent(out) :: drr

   !> First derivative w.r.t. width parameter
   real(wp), intent(out), optional :: dg

   !> Mixed second derivative w.r.t. distance and width parameter
   real(wp), intent(out), optional :: drg

   !> Second derivative w.r.t. width parameter
   real(wp), intent(out), optional :: dgg

   real(wp) :: inv_r, erfval, ex

   inv_r = 1.0_wp / r
   erfval = erf(gam * r)
   ex = (2.0_wp / sqrt(pi)) * exp(-(gam * r)**2)
   val = erfval * inv_r
   dr = (gam * ex - val) * inv_r
   drr = -2.0_wp * (gam**3 * ex + dr * inv_r)
   if (present(dg)) dg = ex
   if (present(drg)) drg = -2.0_wp * gam**2 * r * ex
   if (present(dgg)) dgg = -2.0_wp * gam * r**2 * ex

end subroutine coulomb_kernel


!> Cartesian curvature of a radial kernel
pure function radial_hessian(vec, r, dr, drr) result(box)

   !> Interatomic displacement vector
   real(wp), intent(in) :: vec(3)

   !> Interatomic distance
   real(wp), intent(in) :: r

   !> First radial derivative
   real(wp), intent(in) :: dr

   !> Second radial derivative
   real(wp), intent(in) :: drr

   !> Cartesian curvature block
   real(wp) :: box(3, 3)

   real(wp) :: inv_r, d1, u(3)
   integer :: ic

   inv_r = 1.0_wp / r
   d1 = dr * inv_r
   u = vec * inv_r
   box = (drr - d1) * outerprod(u, u)
   do ic = 1, 3
      box(ic, ic) = box(ic, ic) + d1
   end do

end function radial_hessian


!> Add the Cartesian blocks of a radial pair term (vec = R_i - R_j - T)
subroutine add_pair_hessian(iat, jat, vec, r, dr, drr, hessian, diagonal)

   !> First atom index
   integer, intent(in) :: iat

   !> Second atom index
   integer, intent(in) :: jat

   !> Interatomic displacement vector
   real(wp), intent(in) :: vec(3)

   !> Interatomic distance
   real(wp), intent(in) :: r

   !> First radial derivative
   real(wp), intent(in) :: dr

   !> Second radial derivative
   real(wp), intent(in) :: drr

   !> Cartesian Hessian matrix
   real(wp), intent(inout) :: hessian(:, :)

   !> Optional diagonal block accumulator
   real(wp), intent(inout), optional :: diagonal(:, :, :)

   call add_pair_block(iat, jat, radial_hessian(vec, r, dr, drr), hessian, diagonal)

end subroutine add_pair_hessian


!> Scatter a pair curvature after summing its periodic images
subroutine add_pair_block(iat, jat, box, hessian, diagonal)

   !> First atom index
   integer, intent(in) :: iat

   !> Second atom index
   integer, intent(in) :: jat

   !> Pair curvature block
   real(wp), intent(in) :: box(3, 3)

   !> Cartesian Hessian matrix
   real(wp), intent(inout) :: hessian(:, :)

   !> Separate diagonal accumulator for pair-parallel callers. Each off-diagonal
   !> block must have a unique owner; only this O(N) array then needs reduction.
   real(wp), intent(inout), optional :: diagonal(:, :, :)

   integer :: ii, jj

   ii = 3 * (iat - 1) + 1
   jj = 3 * (jat - 1) + 1
   if (present(diagonal)) then
      diagonal(:, :, iat) = diagonal(:, :, iat) + box
      diagonal(:, :, jat) = diagonal(:, :, jat) + box
   else
      hessian(ii:ii + 2, ii:ii + 2) = hessian(ii:ii + 2, ii:ii + 2) + box
      hessian(jj:jj + 2, jj:jj + 2) = hessian(jj:jj + 2, jj:jj + 2) + box
   end if
   hessian(ii:ii + 2, jj:jj + 2) = hessian(ii:ii + 2, jj:jj + 2) - box
   hessian(jj:jj + 2, ii:ii + 2) = hessian(jj:jj + 2, ii:ii + 2) - box

end subroutine add_pair_block


!> Add F_Rz z_R + transpose and z_R^T F_zz z_R for an atomic variable z
subroutine add_chain_hessian(dzdr, mixed, hessian, second)

   !> Derivative of atomic variable w.r.t. Cartesian coordinates
   real(wp), contiguous, target, intent(in) :: dzdr(:, :, :)

   !> Mixed second derivative tensor
   real(wp), contiguous, target, intent(in) :: mixed(:, :, :)

   !> Cartesian Hessian matrix
   real(wp), intent(inout) :: hessian(:, :)

   !> Optional second derivative matrix w.r.t. atomic variable
   real(wp), intent(in), optional :: second(:, :)

   real(wp), pointer :: dr(:, :), cross(:, :)
   real(wp), allocatable :: work(:, :)
   integer :: ndim, nat

   nat = size(dzdr, 3)
   ndim = size(dzdr, 1) * size(dzdr, 2)
   dr(1:ndim, 1:nat) => dzdr
   cross(1:ndim, 1:nat) => mixed
   call add_symmetric_product(cross, dr, hessian)
   if (present(second)) then
      allocate(work(ndim, nat))
      call gemm(dr, second, work)
      call gemm(work, dr, hessian, transb="t", beta=1.0_wp)
   end if

end subroutine add_chain_hessian


!> Add alpha * (A B^T + B A^T), computing the transposed product only once.
!> Preserve the existing Hessian, including any asymmetric contribution.
subroutine add_symmetric_product(a, b, hessian, alpha)

   !> Left factor matrix
   real(wp), intent(in) :: a(:, :)

   !> Right factor matrix
   real(wp), intent(in) :: b(:, :)

   !> Cartesian Hessian matrix
   real(wp), intent(inout) :: hessian(:, :)

   !> Optional scaling factor
   real(wp), intent(in), optional :: alpha

   real(wp), allocatable :: work(:, :)
   integer :: ii, jj, ndim

   ndim = size(hessian, 1)
   allocate(work(ndim, ndim))
   call gemm(a, b, work, transb="t", alpha=alpha)
   !$omp parallel do collapse(2) default(none) schedule(static) &
   !$omp& shared(hessian, work, ndim) private(ii, jj)
   do jj = 1, ndim
      do ii = 1, ndim
         hessian(ii, jj) = hessian(ii, jj) + work(ii, jj) + work(jj, ii)
      end do
   end do
   !$omp end parallel do

end subroutine add_symmetric_product


end module multicharge_model_hessian
