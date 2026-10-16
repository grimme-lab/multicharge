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

!> Small kernels shared by the charge-response implementations.
module multicharge_model_hessian
   use mctc_env, only: wp
   use mctc_io_constants, only: pi
   use multicharge_blas, only: gemm
   implicit none
   private
   public :: coulomb_kernel, radial_hessian, add_pair_hessian, add_pair_block, add_chain_hessian

contains

!> Gaussian Coulomb kernel and partial derivatives in distance and exponent.
pure subroutine coulomb_kernel(r, gam, val, dr, drr, dg, drg, dgg)
   real(wp), intent(in) :: r, gam
   real(wp), intent(out) :: val, dr, drr
   real(wp), intent(out), optional :: dg, drg, dgg
   real(wp) :: erfval, ex

   erfval = erf(gam*r)
   ex = 2.0_wp/sqrt(pi)*exp(-(gam*r)**2)
   val = erfval/r
   dr = gam*ex/r - erfval/r**2
   drr = -2.0_wp*gam**3*ex - 2.0_wp*gam*ex/r**2 + 2.0_wp*erfval/r**3
   if (present(dg)) dg = ex
   if (present(drg)) drg = -2.0_wp*gam**2*r*ex
   if (present(dgg)) dgg = -2.0_wp*gam*r**2*ex
end subroutine coulomb_kernel


!> Cartesian curvature of a radial kernel.
pure function radial_hessian(vec, dr, drr) result(box)
   real(wp), intent(in) :: vec(3), dr, drr
   real(wp) :: r, box(3, 3)
   integer :: ic

   r = norm2(vec)
   box = (drr-dr/r)*spread(vec/r, 1, 3)*spread(vec/r, 2, 3)
   do ic = 1, 3
      box(ic, ic) = box(ic, ic) + dr/r
   end do
end function radial_hessian


!> Add the Cartesian blocks of a radial pair term (vec=R_i-R_j-T).
subroutine add_pair_hessian(iat, jat, vec, dr, drr, hessian)
   integer, intent(in) :: iat, jat
   real(wp), intent(in) :: vec(3), dr, drr
   real(wp), intent(inout) :: hessian(:, :)

   call add_pair_block(iat, jat, radial_hessian(vec, dr, drr), hessian)
end subroutine add_pair_hessian


!> Scatter a pair curvature after summing its periodic images.
subroutine add_pair_block(iat, jat, box, hessian)
   integer, intent(in) :: iat, jat
   real(wp), intent(in) :: box(3, 3)
   real(wp), intent(inout) :: hessian(:, :)
   integer :: ic, jc, ii, jj, ik, jk

   do jc = 1, 3
      jj = 3*(iat-1)+jc
      jk = 3*(jat-1)+jc
      do ic = 1, 3
         ii = 3*(iat-1)+ic
         ik = 3*(jat-1)+ic
         hessian(ii, jj) = hessian(ii, jj) + box(ic, jc)
         hessian(ik, jk) = hessian(ik, jk) + box(ic, jc)
         hessian(ii, jk) = hessian(ii, jk) - box(ic, jc)
         hessian(ik, jj) = hessian(ik, jj) - box(ic, jc)
      end do
   end do
end subroutine add_pair_block


!> Add F_Rz z_R + transpose and z_R^T F_zz z_R for an atomic variable z.
subroutine add_chain_hessian(dzdr, mixed, hessian, second)
   real(wp), contiguous, target, intent(in) :: dzdr(:, :, :), mixed(:, :, :)
   real(wp), intent(inout) :: hessian(:, :)
   real(wp), intent(in), optional :: second(:, :)
   real(wp), pointer :: dr(:, :), cross(:, :)
   real(wp), allocatable :: work(:, :)
   integer :: ndim, nat

   nat = size(dzdr, 3)
   ndim = size(dzdr, 1)*size(dzdr, 2)
   dr(1:ndim, 1:nat) => dzdr
   cross(1:ndim, 1:nat) => mixed
   call gemm(cross, dr, hessian, transb="t", beta=1.0_wp)
   call gemm(dr, cross, hessian, transb="t", beta=1.0_wp)
   if (present(second)) then
      allocate(work(ndim, nat))
      call gemm(dr, second, work)
      call gemm(work, dr, hessian, transb="t", beta=1.0_wp)
   end if
end subroutine add_chain_hessian

end module multicharge_model_hessian
