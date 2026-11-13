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

!> @file multicharge/ewald/pbc1d.f90
!> Reciprocal-space and G=0 wire kernels for 1D periodic Ewald summation

!> Reciprocal-space and zero-mode potentials, derivatives, and Hessians for 1D Ewald
module multicharge_ewald_pbc1d
   use mctc_env, only : wp
   use mctc_io_math, only : outerprod
   use multicharge_ewald_special, only : fast_f012_1d, fast_k012_1d
   implicit none
   private

   public :: get_amat_rec_1d, get_damat_rec_1d, get_hess_rec_1d

   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: unity(3, 3) = reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])


contains


!> Reciprocal-space and G=0 wire contribution for 1D Ewald.
!> Expects open half-space reciprocal vectors (one representative per +G/-G pair).
pure subroutine get_amat_rec_1d(rij, alpha, len1, u_hat, trans, amat)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic wire length
   real(wp), intent(in) :: len1

   !> Unit vector along the periodic wire axis
   real(wp), intent(in) :: u_hat(3)

   !> Half-space reciprocal lattice vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space Coulomb matrix element contribution
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: r_perp(3), rho2, alp2, inv_alp2, x, y, g_vec(3), g2, gv
   real(wp) :: k0, k1, k2, f0, f1, f2

   amat = 0.0_wp
   r_perp = rij - dot_product(rij, u_hat) * u_hat
   rho2 = dot_product(r_perp, r_perp)
   alp2 = alpha * alpha
   inv_alp2 = 1.0_wp / alp2
   y = alp2 * rho2

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle

      x = 0.25_wp * g2 * inv_alp2
      gv = dot_product(rij, g_vec)
      call fast_k012_1d(x, y, k0, k1, k2)
      amat = amat + cos(gv) * k0
   end do
   call fast_f012_1d(y, f0, f1, f2)
   amat = (2.0_wp * amat - f0 + 2.0_wp * log(alpha)) / len1

end subroutine get_amat_rec_1d


!> Reciprocal-space and G=0 wire derivatives for 1D Ewald.
!> Expects open half-space reciprocal vectors (one representative per +G/-G pair).
pure subroutine get_damat_rec_1d(rij, alpha, len1, u_hat, trans, dg, ds)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic wire length
   real(wp), intent(in) :: len1

   !> Unit vector along the periodic wire axis
   real(wp), intent(in) :: u_hat(3)

   !> Half-space reciprocal lattice vectors
   real(wp), intent(in) :: trans(:, :)

   !> Position derivative contribution
   real(wp), intent(out) :: dg(3)

   !> Strain derivative contribution
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: r_perp(3), rho2, alp2, inv_alp2, inv_len, x, y, g_vec(3), g2, gv
   real(wp) :: cos_gv, sin_gv, k0, k1, k2, f0, f1, f2, exp_xy
   real(wp) :: p_par(3, 3), outer_rr(3, 3), outer_rg(3, 3)

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   r_perp = rij - dot_product(rij, u_hat) * u_hat
   rho2 = dot_product(r_perp, r_perp)
   alp2 = alpha * alpha
   inv_alp2 = 1.0_wp / alp2
   inv_len = 1.0_wp / len1
   y = alp2 * rho2
   p_par = outerprod(u_hat, u_hat)
   outer_rr = outerprod(r_perp, r_perp)

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle

      x = 0.25_wp * g2 * inv_alp2
      gv = dot_product(rij, g_vec)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      call fast_k012_1d(x, y, k0, k1, k2)
      exp_xy = exp(-x - y)

      dg(:) = dg - (sin_gv * k0 * g_vec &
         & + 2.0_wp * alp2 * cos_gv * k1 * r_perp)

      outer_rg = outerprod(g_vec, r_perp) + outerprod(r_perp, g_vec)
      ds(:, :) = ds + cos_gv * (2.0_wp * exp_xy + 2.0_wp * y * k1 - k0) &
         & * p_par - 2.0_wp * alp2 * cos_gv * k1 * outer_rr &
         & - sin_gv * k0 * outer_rg
   end do

   call fast_f012_1d(y, f0, f1, f2)
   dg(:) = (2.0_wp * inv_len) * dg - (2.0_wp * alp2 * f1 * inv_len) * r_perp
   ds(:, :) = (2.0_wp * inv_len) * ds + inv_len * ((f0 - 2.0_wp * log(alpha)) &
      & * p_par - 2.0_wp * alp2 * f1 * outer_rr)

end subroutine get_damat_rec_1d


!> Reciprocal-space and G=0 wire pair Hessian block for 1D Ewald.
!> Expects open half-space reciprocal vectors (one representative per +G/-G pair).
pure subroutine get_hess_rec_1d(rij, alpha, len1, u_hat, trans, rbox)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic wire length
   real(wp), intent(in) :: len1

   !> Unit vector along the periodic wire axis
   real(wp), intent(in) :: u_hat(3)

   !> Half-space reciprocal lattice vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space Hessian block
   real(wp), intent(out) :: rbox(3, 3)

   integer :: itr
   real(wp) :: r_perp(3), rho2, alp2, inv_alp2, inv_len, x, y, g_vec(3), g2, gv
   real(wp) :: cos_gv, sin_gv, k0, k1, k2, f0, f1, f2
   real(wp) :: p_perp(3, 3), gg_mat(3, 3), rg_sym(3, 3), rr_mat(3, 3)

   rbox = 0.0_wp
   p_perp = unity - outerprod(u_hat, u_hat)
   r_perp = rij - dot_product(rij, u_hat) * u_hat
   rho2 = dot_product(r_perp, r_perp)
   alp2 = alpha * alpha
   inv_alp2 = 1.0_wp / alp2
   inv_len = 1.0_wp / len1
   y = alp2 * rho2
   rr_mat = outerprod(r_perp, r_perp)

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle

      x = 0.25_wp * g2 * inv_alp2
      gv = dot_product(g_vec, rij)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      call fast_k012_1d(x, y, k0, k1, k2)
      gg_mat = outerprod(g_vec, g_vec)
      rg_sym = outerprod(g_vec, r_perp) + outerprod(r_perp, g_vec)
      rbox = rbox - cos_gv * k0 * gg_mat &
         & + 2.0_wp * alp2 * sin_gv * k1 * rg_sym &
         & + cos_gv * (-2.0_wp * alp2 * k1 * p_perp &
         & + 4.0_wp * alp2 * alp2 * k2 * rr_mat)
   end do

   call fast_f012_1d(y, f0, f1, f2)
   rbox = (2.0_wp * inv_len) * rbox + inv_len * (-2.0_wp * alp2 * f1 * p_perp &
      & - 4.0_wp * alp2 * alp2 * f2 * rr_mat)

end subroutine get_hess_rec_1d


end module multicharge_ewald_pbc1d
