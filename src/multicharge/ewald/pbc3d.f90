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

!> @file multicharge/ewald/pbc3d.f90
!> Reciprocal-space kernels for 3D periodic Ewald summation

!> Reciprocal-space weights, potentials, derivatives, and Hessians for 3D Ewald
module multicharge_ewald_pbc3d
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : outerprod
   implicit none
   private

   public :: get_rec_weights, get_amat_rec_3d, get_damat_rec_3d, get_hess_rec_3d

   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: unity(3, 3) = reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])


contains


!> Reciprocal-space factors independent of the atom pair.
!> Prepare once for each matrix/response evaluation, including strain when needed.
pure subroutine get_rec_weights(vol, alpha, trans, weight, strain)

   !> Unit cell volume
   real(wp), intent(in) :: vol

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space weights
   real(wp), allocatable, intent(out) :: weight(:)

   !> Optional reciprocal-space strain factors
   real(wp), allocatable, intent(out), optional :: strain(:, :, :)

   integer :: itr
   real(wp) :: g2, inv_g2, inv_alp2, fac, vec(3)

   fac = 4.0_wp * pi / vol
   inv_alp2 = 1.0_wp / alpha**2
   allocate(weight(size(trans, 2)), source=0.0_wp)
   if (present(strain)) then
      allocate(strain(3, 3, size(trans, 2)), source=0.0_wp)
   end if

   do itr = 1, size(trans, 2)
      vec = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle

      inv_g2 = 1.0_wp / g2
      weight(itr) = fac * exp(-0.25_wp * g2 * inv_alp2) * inv_g2
      if (present(strain)) then
         strain(:, :, itr) = weight(itr) * ((2.0_wp * inv_g2 + 0.5_wp * inv_alp2) &
            & * outerprod(vec, vec) - unity)
      end if
   end do

end subroutine get_rec_weights


!> Reciprocal-space contribution to the Coulomb matrix for 3D Ewald
pure subroutine get_amat_rec_3d(rij, trans, weight, amat)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space weights
   real(wp), intent(in) :: weight(:)

   !> Reciprocal-space Coulomb matrix element contribution
   real(wp), intent(out) :: amat

   integer :: itr

   amat = 0.0_wp
   do itr = 1, size(trans, 2)
      amat = amat + weight(itr) * cos(dot_product(rij, trans(:, itr)))
   end do

end subroutine get_amat_rec_3d


!> Reciprocal-space Ewald derivatives of the Coulomb matrix for 3D periodicity
pure subroutine get_damat_rec_3d(rij, trans, weight, strain, dg, ds)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space weights
   real(wp), intent(in) :: weight(:)

   !> Reciprocal-space strain tensors
   real(wp), intent(in) :: strain(:, :, :)

   !> Position derivative contribution
   real(wp), intent(out) :: dg(3)

   !> Strain derivative contribution
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: gv

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   do itr = 1, size(trans, 2)
      gv = dot_product(rij, trans(:, itr))
      dg(:) = dg - weight(itr) * sin(gv) * trans(:, itr)
      ds(:, :) = ds + cos(gv) * strain(:, :, itr)
   end do

end subroutine get_damat_rec_3d


!> Reciprocal-space pair Hessian block for 3D Ewald
pure subroutine get_hess_rec_3d(rij, trans, weight, rbox)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space weights
   real(wp), intent(in) :: weight(:)

   !> Reciprocal-space Hessian block
   real(wp), intent(out) :: rbox(3, 3)

   integer :: itr
   real(wp) :: vec(3)

   rbox = 0.0_wp
   do itr = 1, size(trans, 2)
      vec = trans(:, itr)
      rbox = rbox - weight(itr) * cos(dot_product(rij, vec)) &
         & * outerprod(vec, vec)
   end do

end subroutine get_hess_rec_3d


end module multicharge_ewald_pbc3d
