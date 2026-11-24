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

!> @file multicharge/ewald/direct.f90
!> Direct-space short-range Ewald kernels

!> Short-range real-space Ewald potentials, derivatives, and Hessians
module multicharge_ewald_direct
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : outerprod
   use multicharge_model_hessian, only : coulomb_kernel, radial_hessian
   implicit none
   private

   public :: get_amat_dir, get_damat_dir, get_hess_dir

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))


contains


!> Direct-space Ewald contribution to the Coulomb matrix
pure subroutine get_amat_dir(rij, gam, alp, trans, amat)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Effective Gaussian charge width parameter
   real(wp), intent(in) :: gam

   !> Ewald splitting parameter
   real(wp), intent(in) :: alp

   !> Direct-space lattice translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Direct-space Coulomb matrix element contribution
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle

      tmp = (erf(gam * r1) - erf(alp * r1)) / r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir


!> Direct-space Ewald derivatives of the Coulomb matrix
pure subroutine get_damat_dir(rij, gam, alp, trans, dg, ds)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Effective Gaussian charge width parameter
   real(wp), intent(in) :: gam

   !> Ewald splitting parameter
   real(wp), intent(in) :: alp

   !> Direct-space lattice translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Position derivative contribution
   real(wp), intent(out) :: dg(3)

   !> Strain derivative contribution
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, inv_r, inv_r2, gtmp, atmp, dtmp
   real(wp) :: gam2, alp2, fac_g, fac_a

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   gam2 = gam * gam
   alp2 = alp * alp
   fac_g = (2.0_wp / sqrtpi) * gam
   fac_a = (2.0_wp / sqrtpi) * alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle

      r2 = r1 * r1
      inv_r = 1.0_wp / r1
      inv_r2 = inv_r * inv_r
      gtmp = fac_g * exp(-r2 * gam2) - erf(r1 * gam) * inv_r
      atmp = erf(r1 * alp) * inv_r - fac_a * exp(-r2 * alp2)
      dtmp = (gtmp + atmp) * inv_r2
      dg(:) = dg + dtmp * vec
      ds(:, :) = ds + dtmp * outerprod(vec, vec)
   end do

end subroutine get_damat_dir


!> Direct-space Ewald pair Hessian block
pure subroutine get_hess_dir(rij, gam, alp, trans, dbox)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Effective Gaussian charge width parameter
   real(wp), intent(in) :: gam

   !> Ewald splitting parameter
   real(wp), intent(in) :: alp

   !> Direct-space lattice translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Direct-space Hessian block
   real(wp), intent(out) :: dbox(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, val, dr_g, drr_g, dr_a, drr_a

   dbox = 0.0_wp

   do itr = 1, size(trans, 2)
      vec = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle

      call coulomb_kernel(r1, gam, val, dr_g, drr_g)
      call coulomb_kernel(r1, alp, val, dr_a, drr_a)
      dbox = dbox + radial_hessian(vec, dr_g - dr_a, drr_g - drr_a)
   end do

end subroutine get_hess_dir


end module multicharge_ewald_direct
