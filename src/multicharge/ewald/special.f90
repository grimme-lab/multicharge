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

!> @file multicharge/ewald/special.f90
!> Special functions for 1D and 2D periodic Ewald summation

!> Overflow-safe complementary error function products and 1D incomplete
!> Bessel/exponential integrals for mixed-periodicity Ewald summation.
module multicharge_ewald_special
   use mctc_env, only : wp
   implicit none
   private

   public :: erfc_exp_2d, fast_f012_1d, fast_k012_1d

   real(wp), parameter :: euler_gamma = 0.57721566490153286060_wp

   !> Positive nodes of 16-point Gauss-Legendre quadrature on [-1, 1]
   real(wp), parameter :: gl16_x(8) = [ &
      & 0.095012509837637440185_wp, &
      & 0.281603550779258913230_wp, &
      & 0.458016777657227386342_wp, &
      & 0.617876244402643748447_wp, &
      & 0.755404408355003033895_wp, &
      & 0.865631202387831743880_wp, &
      & 0.944575023073232576078_wp, &
      & 0.989400934991649932596_wp]

   !> Weights of 16-point Gauss-Legendre quadrature on [-1, 1]
   real(wp), parameter :: gl16_w(8) = [ &
      & 0.189450610455068496285_wp, &
      & 0.182603415044923588867_wp, &
      & 0.169156519395002538189_wp, &
      & 0.149595988816576732081_wp, &
      & 0.124628971255533872052_wp, &
      & 0.095158511682492784810_wp, &
      & 0.062253523938647892863_wp, &
      & 0.027152459411754094852_wp]


contains


!> Overflow-safe evaluation of exp(G * s) * erfc(G / (2 * alpha) + alpha * s)
pure function erfc_exp_2d(g_norm, s, alpha) result(val)

   !> Norm of the 2D reciprocal vector (G > 0)
   real(wp), intent(in) :: g_norm

   !> Signed out-of-plane coordinate (+z or -z)
   real(wp), intent(in) :: s

   !> Ewald splitting parameter (alpha > 0)
   real(wp), intent(in) :: alpha

   !> Value of exp(G * s) * erfc(G / (2 * alpha) + alpha * s)
   real(wp) :: val

   real(wp) :: u, exponent

   u = 0.5_wp * g_norm / alpha + alpha * s
   if (u > 0.0_wp) then
      exponent = -0.25_wp * (g_norm / alpha)**2 - (alpha * s)**2
      val = exp(exponent) * erfc_scaled(u)
   else
      val = exp(g_norm * s) * erfc(u)
   end if

end function erfc_exp_2d


!> Evaluate 1D zero-mode functions f0(u), f1(u), f2(u) for u = alpha^2 * rho^2 >= 0
pure subroutine fast_f012_1d(u, f0, f1, f2)

   !> Non-negative scalar argument u = alpha^2 * rho^2
   real(wp), intent(in) :: u

   !> Zero-mode potential function f0(u) = gamma_E + E1(u) + ln(u)
   real(wp), intent(out) :: f0

   !> First derivative f1(u) = f0'(u) = (1 - exp(-u)) / u
   real(wp), intent(out) :: f1

   !> Second derivative f2(u) = f0''(u) = ((u + 1)*exp(-u) - 1) / u^2
   real(wp), intent(out) :: f2

   real(wp) :: term, cf, e1_val, exp_u
   integer :: n, k

   if (u <= 0.0_wp) then
      f0 = 0.0_wp
      f1 = 1.0_wp
      f2 = -0.5_wp
      return
   end if

   ! Power series for small u, continued fraction for E1(u) at larger u
   if (u <= 4.0_wp) then
      f0 = 0.0_wp
      term = u
      do n = 1, 35
         f0 = f0 + term / real(n, wp)
         term = -term * u / real(n + 1, wp)
      end do
   else
      cf = 0.0_wp
      do k = 40, 1, -1
         cf = -real(k * k, wp) / (u + real(2 * k + 1, wp) + cf)
      end do
      e1_val = exp(-u) / (u + 1.0_wp + cf)
      f0 = euler_gamma + e1_val + log(u)
   end if

   ! Taylor expansion for f1 and f2 near u = 0 to avoid cancellation
   if (u <= 0.5_wp) then
      f1 = 0.0_wp
      f2 = 0.0_wp
      term = 1.0_wp
      do n = 0, 19
         f1 = f1 + term / real(n + 1, wp)
         f2 = f2 - term / real(n + 2, wp)
         term = -term * u / real(n + 1, wp)
      end do
   else
      exp_u = exp(-u)
      f1 = (1.0_wp - exp_u) / u
      f2 = ((u + 1.0_wp) * exp_u - 1.0_wp) / (u * u)
   end if

end subroutine fast_f012_1d


!> Evaluate incomplete Bessel/exponential integrals K0(x, y), K1(x, y), K2(x, y)
!> using composite 16-point Gauss-Legendre quadrature in w = ln(t).
pure subroutine fast_k012_1d(x, y, k0, k1, k2)

   !> Positive reciprocal argument x = G^2 / (4 * alpha^2)
   real(wp), intent(in) :: x

   !> Non-negative radial argument y = alpha^2 * rho^2
   real(wp), intent(in) :: y

   !> Integral K0(x, y)
   real(wp), intent(out) :: k0

   !> Integral K1(x, y) = -dK0/dy
   real(wp), intent(out) :: k1

   !> Integral K2(x, y) = d^2K0/dy^2
   real(wp), intent(out) :: k2

   real(wp) :: w0, psi_min, s_cut, disc, w_min, w_max, s0, s1, s2

   ! Locate the saddle point of the exponent psi(w) = x*exp(w) + y*exp(-w)
   if (y > x) then
      w0 = 0.5_wp * log(y / x)
      psi_min = 2.0_wp * sqrt(x * y)
   else
      w0 = 0.0_wp
      psi_min = x + y
   end if

   if (psi_min > 700.0_wp) then
      k0 = 0.0_wp
      k1 = 0.0_wp
      k2 = 0.0_wp
      return
   end if

   ! Truncate integration domain where integrand falls below machine precision
   s_cut = psi_min + 44.0_wp
   disc = sqrt(max(0.0_wp, s_cut * s_cut - 4.0_wp * x * y))
   w_max = log((s_cut + disc) / (2.0_wp * x))
   if (w_max <= 0.0_wp) then
      k0 = 0.0_wp
      k1 = 0.0_wp
      k2 = 0.0_wp
      return
   end if

   ! Split the integral at the minimum of the exponent
   if (w0 > 0.0_wp .and. w0 < w_max) then
      w_min = max(0.0_wp, log(y / (psi_min + 52.0_wp)))
      call integrate_k012_interval(x, y, w_min, w0, k0, k1, k2)
      call integrate_k012_interval(x, y, w0, w_max, s0, s1, s2)
      k0 = k0 + s0
      k1 = k1 + s1
      k2 = k2 + s2
   else
      call integrate_k012_interval(x, y, 0.0_wp, w_max, k0, k1, k2)
   end if

end subroutine fast_k012_1d


!> Composite 16-point Gauss-Legendre quadrature of exp(-x*exp(w) - y*exp(-w))
pure subroutine integrate_k012_interval(x, y, a_start, b_end, s0, s1, s2)

   !> Positive reciprocal argument x
   real(wp), intent(in) :: x

   !> Non-negative radial argument y
   real(wp), intent(in) :: y

   !> Lower integration bound in w = ln(t)
   real(wp), intent(in) :: a_start

   !> Upper integration bound in w = ln(t)
   real(wp), intent(in) :: b_end

   !> Accumulated K0 contribution
   real(wp), intent(out) :: s0

   !> Accumulated K1 contribution
   real(wp), intent(out) :: s1

   !> Accumulated K2 contribution
   real(wp), intent(out) :: s2

   integer :: n_sub, isub, q
   real(wp) :: h, a_sub, wt, w, ew, emw, val

   n_sub = max(4, ceiling((b_end - a_start) / 0.25_wp))
   h = (b_end - a_start) / real(n_sub, wp)
   s0 = 0.0_wp
   s1 = 0.0_wp
   s2 = 0.0_wp

   do isub = 0, n_sub - 1
      a_sub = a_start + real(isub, wp) * h
      do q = 1, 8
         wt = 0.5_wp * h * gl16_w(q)

         w = a_sub + 0.5_wp * h * (1.0_wp - gl16_x(q))
         ew = exp(w)
         emw = 1.0_wp / ew
         val = wt * exp(-x * ew - y * emw)
         s0 = s0 + val
         s1 = s1 + val * emw
         s2 = s2 + val * emw * emw

         w = a_sub + 0.5_wp * h * (1.0_wp + gl16_x(q))
         ew = exp(w)
         emw = 1.0_wp / ew
         val = wt * exp(-x * ew - y * emw)
         s0 = s0 + val
         s1 = s1 + val * emw
         s2 = s2 + val * emw * emw
      end do
   end do

end subroutine integrate_k012_interval


end module multicharge_ewald_special
