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

!> @file multicharge/ewald/pbc2d.f90
!> Reciprocal-space and G=0 sheet kernels for 2D periodic Ewald summation

!> Reciprocal-space and zero-mode potentials, derivatives, and Hessians for 2D Ewald
module multicharge_ewald_pbc2d
   use mctc_env, only : wp
   use mctc_io_constants, only : pi
   use mctc_io_math, only : outerprod
   use multicharge_ewald_special, only : erfc_exp_2d
   implicit none
   private

   public :: get_amat_rec_2d, get_damat_rec_2d, get_hess_rec_2d

   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
   real(wp), parameter :: unity(3, 3) = reshape([ &
      & 1.0_wp, 0.0_wp, 0.0_wp, &
      & 0.0_wp, 1.0_wp, 0.0_wp, &
      & 0.0_wp, 0.0_wp, 1.0_wp], [3, 3])


contains


!> Reciprocal-space and G=0 sheet contribution for 2D Ewald
pure subroutine get_amat_rec_2d(rij, alpha, area, n_hat, trans, amat)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic surface area
   real(wp), intent(in) :: area

   !> Unit normal vector to the periodic plane
   real(wp), intent(in) :: n_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space Coulomb matrix element contribution
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: z, pi_area, g_vec(3), g_norm, gv, ep, em

   z = dot_product(rij, n_hat)
   pi_area = pi / area
   amat = 0.0_wp
   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g_norm = norm2(g_vec)
      if (g_norm < eps) cycle

      gv = dot_product(rij, g_vec)
      ep = erfc_exp_2d(g_norm, +z, alpha)
      em = erfc_exp_2d(g_norm, -z, alpha)
      amat = amat + (pi_area / g_norm) * cos(gv) * (ep + em)
   end do
   amat = amat - 2.0_wp * pi_area &
      & * (z * erf(alpha * z) + exp(-(alpha * z)**2) / (alpha * sqrtpi))

end subroutine get_amat_rec_2d


!> Reciprocal-space and G=0 sheet derivatives for 2D Ewald
pure subroutine get_damat_rec_2d(rij, alpha, area, n_hat, trans, dg, ds)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic surface area
   real(wp), intent(in) :: area

   !> Unit normal vector to the periodic plane
   real(wp), intent(in) :: n_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Position derivative contribution
   real(wp), intent(out) :: dg(3)

   !> Strain derivative contribution
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: z, pi_area, inv_alp2, z2_alp2, fac_dg, r_perp(3), g_vec(3)
   real(wp) :: g_norm, inv_g, gv, cos_gv, sin_gv
   real(wp) :: ep, em, f_gz, df_dg, df_dz, pref, phi_0, dphi0_dz
   real(wp) :: i_par(3, 3), nn(3, 3), outer_gg(3, 3), outer_rg(3, 3)

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   z = dot_product(rij, n_hat)
   pi_area = pi / area
   inv_alp2 = 1.0_wp / alpha**2
   z2_alp2 = (alpha * z)**2
   fac_dg = 2.0_wp / (alpha * sqrtpi)
   r_perp = z * n_hat
   nn = outerprod(n_hat, n_hat)
   i_par = unity - nn

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g_norm = norm2(g_vec)
      if (g_norm < eps) cycle

      inv_g = 1.0_wp / g_norm
      gv = dot_product(rij, g_vec)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      ep = erfc_exp_2d(g_norm, +z, alpha)
      em = erfc_exp_2d(g_norm, -z, alpha)
      f_gz = ep + em
      df_dg = z * (ep - em) &
         & - fac_dg * exp(-0.25_wp * g_norm**2 * inv_alp2 - z2_alp2)
      df_dz = g_norm * (ep - em)
      pref = pi_area * inv_g

      dg(:) = dg + pref * (-sin_gv * f_gz * g_vec + cos_gv * df_dz * n_hat)

      outer_gg = outerprod(g_vec, g_vec)
      outer_rg = outerprod(r_perp, g_vec) + outerprod(g_vec, r_perp)
      ds(:, :) = ds + pref * cos_gv * (-f_gz * i_par &
         & + (f_gz * inv_g - df_dg) * inv_g * outer_gg + z * df_dz * nn) &
         & - pref * sin_gv * f_gz * outer_rg
   end do

   phi_0 = -2.0_wp * pi_area &
      & * (z * erf(alpha * z) + 0.5_wp * fac_dg * exp(-z2_alp2))
   dphi0_dz = -2.0_wp * pi_area * erf(alpha * z)
   dg(:) = dg + dphi0_dz * n_hat
   ds(:, :) = ds - phi_0 * i_par + dphi0_dz * z * nn

end subroutine get_damat_rec_2d


!> Reciprocal-space and G=0 sheet pair Hessian block for 2D Ewald
pure subroutine get_hess_rec_2d(rij, alpha, area, n_hat, trans, rbox)

   !> Interatomic displacement vector
   real(wp), intent(in) :: rij(3)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Periodic surface area
   real(wp), intent(in) :: area

   !> Unit normal vector to the periodic plane
   real(wp), intent(in) :: n_hat(3)

   !> Reciprocal-space translation vectors
   real(wp), intent(in) :: trans(:, :)

   !> Reciprocal-space Hessian block
   real(wp), intent(out) :: rbox(3, 3)

   integer :: itr
   real(wp) :: z, pi_area, inv_alp2, z2_alp2, g_vec(3), g2, g_norm, gv
   real(wp) :: cos_gv, sin_gv, ep, em, fp, fm, gauss, pref, gauss_z
   real(wp) :: gg_mat(3, 3), nn_mat(3, 3), ng_sym(3, 3)

   rbox = 0.0_wp
   z = dot_product(rij, n_hat)
   pi_area = pi / area
   inv_alp2 = 1.0_wp / (alpha * alpha)
   z2_alp2 = (alpha * z)**2
   nn_mat = outerprod(n_hat, n_hat)

   do itr = 1, size(trans, 2)
      g_vec = trans(:, itr)
      g2 = dot_product(g_vec, g_vec)
      if (g2 < eps) cycle

      g_norm = sqrt(g2)
      gv = dot_product(g_vec, rij)
      cos_gv = cos(gv)
      sin_gv = sin(gv)
      ep = erfc_exp_2d(g_norm, +z, alpha)
      em = erfc_exp_2d(g_norm, -z, alpha)
      fp = ep + em
      fm = ep - em
      gauss = exp(-0.25_wp * g2 * inv_alp2 - z2_alp2)
      pref = pi_area / g_norm
      gg_mat = outerprod(g_vec, g_vec)
      ng_sym = outerprod(g_vec, n_hat) + outerprod(n_hat, g_vec)
      rbox = rbox - pref * cos_gv * fp * gg_mat &
         & - pi_area * sin_gv * fm * ng_sym &
         & + pi_area * cos_gv &
         & * (g_norm * fp - (4.0_wp / sqrtpi) * alpha * gauss) * nn_mat
   end do

   gauss_z = exp(-z2_alp2)
   rbox = rbox - (4.0_wp * sqrtpi * alpha / area) * gauss_z * nn_mat

end subroutine get_hess_rec_2d


end module multicharge_ewald_pbc2d
