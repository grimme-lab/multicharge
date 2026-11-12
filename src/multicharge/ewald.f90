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

!> @file multicharge/ewald.f90
!> Ewald summation utilities for 1D, 2D, and 3D periodic boundary conditions.
!>
!> The 1D and 2D Ewald decompositions follow:
!> - F. Nestler and D. Potts, *Fast Ewald summation under 2d- and 1d-periodic
!>   boundary conditions based on NFFTs*, Proc. Sampling Theory and Applications
!>   (SampTA), 329-332 (2013).
!> - F. Nestler, M. Pippig, and D. Potts, *Fast Ewald summation based on NFFT with
!>   mixed periodicity*, J. Comput. Phys. **285**, 280-315 (2015),
!>   DOI: 10.1016/j.jcp.2014.12.052.
!> This implementation evaluates their 1D and 2D Ewald formulations by direct
!> lattice summation.
module multicharge_ewald
   use mctc_cutoff, only: get_inverse_metric, get_lattice_points, get_periodic_axes
   use mctc_env, only: error_type, wp
   use mctc_io_constants, only: pi
   use mctc_io_math, only: crossprod, matdet_3x3
   implicit none
   private

   public :: get_alpha, get_ewald_dir_trans, get_ewald_rec_trans
   public :: erfc_exp_2d, fast_f012_1d, fast_k012_1d

   real(wp), parameter :: twopi = 2.0_wp * pi
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))
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

   abstract interface
      !> Returns the maximum value of a reciprocal-space Ewald term for a given vector length
      pure function get_rec_term_gen(gg, alpha, vol) result(gTerm)
         import :: wp

         !> Length of the reciprocal space vector
         real(wp), intent(in) :: gg

         !> Parameter of the Ewald summation
         real(wp), intent(in) :: alpha

         !> Volume (or area/length) of the real space unit cell
         real(wp), intent(in) :: vol

         !> Reciprocal term
         real(wp) :: gTerm

      end function get_rec_term_gen
   end interface


contains


!> Determine the optimal Ewald splitting parameter alpha
subroutine get_alpha(lattice, alpha, error, periodic, vol, rec_lat)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Ewald splitting parameter
   real(wp), intent(out) :: alpha

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional periodicity mask (defaults to 3D periodic)
   logical, intent(in), optional :: periodic(:)

   !> Optional unit cell volume (or area/length)
   real(wp), intent(out), optional :: vol

   !> Optional reciprocal lattice vectors
   real(wp), intent(out), optional :: rec_lat(3, 3)

   integer :: ndim, p_axes(3)
   real(wp) :: volume, red_lat(3, 3), rlat(3, 3)

   alpha = 0.0_wp

   if (present(periodic)) then
      call prepare_periodic_lattice(periodic, lattice, ndim, p_axes, &
         & red_lat, rlat, volume, error)
   else
      call prepare_periodic_lattice([.true., .true., .true.], lattice, &
         & ndim, p_axes, red_lat, rlat, volume, error)
   end if
   if (allocated(error)) return

   ! Balance the decay rates of the direct- and reciprocal-space Ewald sums
   select case (ndim)
   case (0)
      alpha = 0.25_wp
   case (1)
      call search_alpha(red_lat(:, p_axes(:1)), rlat(:, p_axes(:1)), &
         & volume, eps, alpha, get_rec_term_1d)
   case (2)
      call search_alpha(red_lat(:, p_axes(:2)), rlat(:, p_axes(:2)), &
         & volume, eps, alpha, get_rec_term_2d)
   case (3)
      call search_alpha(red_lat, rlat, volume, eps, alpha, get_rec_term_3d)
      ! Prevent small alpha and excessive direct-space boxes in anisotropic 3D cells
      if (volume > eps) then
         alpha = max(alpha, sqrt(pi) / volume**(1.0_wp / 3.0_wp))
      end if
   case default
      alpha = 0.25_wp
   end select

   if (present(vol)) vol = volume
   if (present(rec_lat)) rec_lat = rlat

end subroutine get_alpha


!> Generate converged direct-space lattice translations for 1D, 2D, or 3D Ewald summation
subroutine get_ewald_dir_trans(periodic, lattice, alpha, gam_min, trans, error, max_rij)

   !> Periodicity mask along each lattice axis
   logical, intent(in) :: periodic(:)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Minimum Gaussian exponent min_ij(gamma_ij)
   real(wp), intent(in) :: gam_min

   !> Direct-space lattice translations (3, ntrans)
   real(wp), allocatable, intent(out) :: trans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional maximum pair displacement within the Wigner-Seitz cell
   real(wp), intent(in), optional :: max_rij

   integer :: ndim, p_axes(3), i
   real(wp) :: red_lat(3, 3), rec_lat(3, 3), vol
   real(wp) :: kappa, kappa_max, tol_val, u_cut, rcut, r_disp, r_dir

   call prepare_periodic_lattice(periodic, lattice, ndim, p_axes, red_lat, rec_lat, vol, error)
   if (allocated(error)) return
   if (ndim == 0) then
      allocate(trans(3, 1), source=0.0_wp)
      return
   end if

   ! Determine direct-space cutoff radius from slowest Gaussian/Ewald decay rate
   tol_val = epsilon(1.0_wp)
   kappa = min(alpha, gam_min)
   kappa_max = max(alpha, gam_min)
   u_cut = sqrt(max(42.25_wp, -log(tol_val) &
      & + log(max(1.0_wp, 4.0_wp * kappa_max**3)) + 4.0_wp))
   rcut = u_cut / kappa
   do while (get_dir_term(rcut, kappa) > tol_val)
      rcut = 1.1_wp * rcut
   end do

   ! The cutoff includes the maximum Wigner-Seitz pair displacement because sums are evaluated for shifted pairs
   r_disp = 0.0_wp
   do i = 1, ndim
      r_disp = r_disp + 0.5_wp * norm2(red_lat(:, p_axes(i)))
   end do
   if (present(max_rij)) r_disp = max(r_disp, max_rij)
   r_dir = rcut + r_disp

   call get_lattice_points(periodic, red_lat, r_dir, trans, error)
   if (allocated(error)) return

end subroutine get_ewald_dir_trans


!> Generate converged nonzero reciprocal-space lattice vectors for 1D, 2D, or 3D Ewald summation
subroutine get_ewald_rec_trans(periodic, lattice, alpha, trans, error, vol)

   !> Periodicity mask along each lattice axis
   logical, intent(in) :: periodic(:)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Ewald splitting parameter
   real(wp), intent(in) :: alpha

   !> Nonzero reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: trans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional generalized cell measure (length in 1D, area in 2D, volume in 3D)
   real(wp), intent(out), optional :: vol

   integer :: ndim, p_axes(3), rep(3), i, ax
   real(wp) :: red_lat(3, 3), rec_lat(3, 3), cell_vol
   real(wp) :: tol_val, u_cut, gcut

   call prepare_periodic_lattice(periodic, lattice, ndim, p_axes, red_lat, rec_lat, cell_vol, error)
   if (present(vol)) vol = cell_vol
   if (allocated(error)) return
   if (ndim == 0) then
      allocate(trans(3, 0))
      return
   end if

   ! Determine reciprocal-space cutoff radius for each periodic dimensionality
   tol_val = epsilon(1.0_wp)
   select case (ndim)
   case (1)
      u_cut = sqrt(max(42.25_wp, -log(tol_val) &
         & + log(max(1.0_wp, 4.0_wp * max(1.0_wp, alpha**2) / cell_vol)) + 4.0_wp))
      gcut = 2.0_wp * alpha * u_cut
      do while (get_rec_term_1d(gcut, alpha, cell_vol) > tol_val)
         gcut = 1.1_wp * gcut
      end do
   case (2)
      u_cut = sqrt(max(42.25_wp, -log(tol_val) &
         & + log(max(1.0_wp, 8.0_wp * sqrt(pi) * max(alpha, alpha**3) / cell_vol)) + 4.0_wp))
      gcut = 2.0_wp * alpha * u_cut
      do while (get_rec_term_2d(gcut, alpha, cell_vol) > tol_val)
         gcut = 1.1_wp * gcut
      end do
   case (3)
      u_cut = sqrt(max(42.25_wp, -log(tol_val) &
         & + log(max(1.0_wp, 16.0_wp * pi * max(1.0_wp, alpha**2) / cell_vol)) + 4.0_wp))
      gcut = 2.0_wp * alpha * u_cut
      do while (get_rec_term_3d(gcut, alpha, cell_vol) > tol_val)
         gcut = 1.1_wp * gcut
      end do
   case default
      gcut = 0.0_wp
   end select

   ! Project reciprocal cutoff sphere onto real-space basis: |m_i| <= gcut * |u_i| / (2*pi)
   rep = 0
   do i = 1, ndim
      ax = p_axes(i)
      rep(ax) = max(2, ceiling(gcut * norm2(red_lat(:, ax)) / twopi))
   end do

   call get_lattice_points(rec_lat, rep, .false., trans, error)
   if (allocated(error)) return

end subroutine get_ewald_rec_trans


!> Prepare reduced direct and reciprocal lattice vectors for a periodic subspace
subroutine prepare_periodic_lattice(periodic, lattice, ndim, p_axes, red_lat, rec_lat, vol, error)

   !> Periodicity mask along each lattice axis
   logical, intent(in) :: periodic(:)

   !> Real-space lattice vectors as columns (3, 3)
   real(wp), intent(in) :: lattice(:, :)

   !> Number of periodic dimensions
   integer, intent(out) :: ndim

   !> Active periodic axis indices
   integer, intent(out) :: p_axes(3)

   !> Size-reduced real-space lattice vectors (3, 3)
   real(wp), intent(out) :: red_lat(3, 3)

   !> Reciprocal-space lattice vectors (3, 3)
   real(wp), intent(out) :: rec_lat(3, 3)

   !> Generalized cell measure (length in 1D, area in 2D, volume in 3D)
   real(wp), intent(out) :: vol

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   real(wp) :: basis(3, 3), ginv(3, 3)

   red_lat = 0.0_wp
   rec_lat = 0.0_wp
   vol = 0.0_wp

   call get_periodic_axes(periodic, ndim, p_axes)

   call get_inverse_metric(periodic, lattice, basis, ginv, error)
   if (allocated(error)) return

   ! Reduce the periodic basis to shorten direct and reciprocal box bounds
   call reduce_basis(ndim, p_axes, basis)

   ! Compute inverse metric G^(-1) and reciprocal basis B = 2*pi * L * G^(-1)
   call get_inverse_metric(periodic, basis, red_lat, ginv, error)
   if (allocated(error)) return
   rec_lat = twopi * matmul(red_lat, ginv)

   select case (ndim)
   case (1)
      vol = norm2(red_lat(:, p_axes(1)))
   case (2)
      vol = norm2(crossprod(red_lat(:, p_axes(1)), red_lat(:, p_axes(2))))
   case (3)
      vol = abs(matdet_3x3(red_lat))
   case default
      vol = 0.0_wp
   end select

end subroutine prepare_periodic_lattice


!> Reduce periodic lattice basis by iterative pairwise size reduction
pure subroutine reduce_basis(ndim, p_axes, lat)

   !> Number of periodic dimensions
   integer, intent(in) :: ndim

   !> Active periodic axes
   integer, intent(in) :: p_axes(3)

   !> Lattice vectors, reduced in-place on periodic columns
   real(wp), intent(inout) :: lat(3, 3)

   integer :: iter, i, j, ai, aj
   real(wp) :: ni, nj, dot_ij, mu, tmp(3)
   logical :: changed

   if (ndim < 2) return

   ! Integer shears preserve the lattice while shortening the translation search
   do iter = 1, 100
      changed = .false.
      do i = 1, ndim - 1
         do j = i + 1, ndim
            ai = p_axes(i)
            aj = p_axes(j)
            ni = dot_product(lat(:, ai), lat(:, ai))
            nj = dot_product(lat(:, aj), lat(:, aj))
            if (ni > nj * (1.0_wp + 1.0e-12_wp)) then
               tmp = lat(:, ai)
               lat(:, ai) = lat(:, aj)
               lat(:, aj) = tmp
               ni = nj
               changed = .true.
            end if
            if (ni <= tiny(1.0_wp)) cycle

            dot_ij = dot_product(lat(:, ai), lat(:, aj))
            mu = anint(dot_ij / ni)
            if (abs(mu) >= 1.0_wp .and. abs(dot_ij) > 0.5_wp * ni * (1.0_wp + 1.0e-12_wp)) then
               lat(:, aj) = lat(:, aj) - mu * lat(:, ai)
               changed = .true.
            end if
         end do
      end do
      if (.not. changed) exit
   end do

end subroutine reduce_basis


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

   real(wp) :: w0, psi_min, s_cut, w_min, w_max, s0, s1, s2

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
   w_max = log((s_cut + sqrt(max(0.0_wp, s_cut * s_cut - 4.0_wp * x * y))) / (2.0_wp * x))
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


!> Composite 16-point Gauss-Legendre quadrature of exp(-x*exp(w) - y*exp(-w)) on [a_start, b_end]
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


!> Get optimal alpha-parameter for the Ewald summation by finding alpha where
!> decline of real and reciprocal part of Ewald are equal.
subroutine search_alpha(lattice, rec_lat, volume, tolerance, alpha, get_rec_term)

   !> Lattice vectors
   real(wp), intent(in) :: lattice(:, :)

   !> Reciprocal vectors
   real(wp), intent(in) :: rec_lat(:, :)

   !> Volume (or area/length) of the unit cell
   real(wp), intent(in) :: volume

   !> Tolerance for difference in real and rec. part
   real(wp), intent(in) :: tolerance

   !> Optimal alpha
   real(wp), intent(out) :: alpha

   !> Reciprocal term evaluator
   procedure(get_rec_term_gen) :: get_rec_term

   real(wp) :: alpl, alpr, rlen, dlen, diff
   real(wp), parameter :: alpha0 = 1.0e-8_wp
   integer, parameter :: niter = 30
   integer :: ibs, stat

   rlen = sqrt(minval(sum(rec_lat(:, :)**2, dim=1)))
   dlen = sqrt(minval(sum(lattice(:, :)**2, dim=1)))

   ! The search balances decay of real- and reciprocal-space contributions; failure uses a fixed fallback
   stat = 0
   alpha = alpha0
   diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
   do while (diff < -tolerance .and. alpha <= huge(1.0_wp))
      alpha = 2.0_wp * alpha
      diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
   end do
   if (alpha > huge(1.0_wp)) then
      stat = 1
   else if (alpha == alpha0) then
      stat = 2
   end if

   if (stat == 0) then
      alpl = 0.5_wp * alpha
      do while (diff < tolerance .and. alpha <= huge(1.0_wp))
         alpha = 2.0_wp * alpha
         diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
      end do
      if (alpha > huge(1.0_wp)) then
         stat = 3
      end if
   end if

   if (stat == 0) then
      alpr = alpha
      alpha = (alpl + alpr) * 0.5_wp
      ibs = 0
      diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
      do while (abs(diff) > tolerance .and. ibs <= niter)
         if (diff < 0) then
            alpl = alpha
         else
            alpr = alpha
         end if
         alpha = (alpl + alpr) * 0.5_wp
         diff = rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume)
         ibs = ibs + 1
      end do
      if (ibs > niter) then
         stat = 4
      end if
   end if

   if (stat /= 0) then
      alpha = 0.25_wp
   end if

end subroutine search_alpha


!> Returns the difference in the decrease of the real and reciprocal parts of the
!> Ewald sum. In order to make the real space part shorter than the reciprocal
!> space part, the values are taken at different distances for the real and the
!> reciprocal space parts.
pure function rec_dir_diff(alpha, get_rec_term, rlen, dlen, volume) result(diff)

   !> Parameter for the Ewald summation
   real(wp), intent(in) :: alpha

   !> Procedure pointer to reciprocal routine
   procedure(get_rec_term_gen) :: get_rec_term

   !> Length of the shortest reciprocal space vector in the sum
   real(wp), intent(in) :: rlen

   !> Length of the shortest real space vector in the sum
   real(wp), intent(in) :: dlen

   !> Volume of the real space unit cell
   real(wp), intent(in) :: volume

   !> Difference between changes in the two terms
   real(wp) :: diff

   diff = ((get_rec_term(4*rlen, alpha, volume) &
      & - get_rec_term(5*rlen, alpha, volume))) &
      & - (get_dir_term(2*dlen, alpha) - get_dir_term(3*dlen, alpha))

end function rec_dir_diff


!> Returns the maximum value of a term in the real space part of the Ewald summation
!> for a given vector length.
pure function get_dir_term(rr, alpha) result(dval)

   !> Length of the real space vector
   real(wp), intent(in) :: rr

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Real space term
   real(wp) :: dval

   dval = erfc(alpha*rr)/rr

end function get_dir_term


!> Returns the maximum value of a term in the 3D reciprocal space part of the Ewald
!> summation for a given vector length.
pure function get_rec_term_3d(gg, alpha, vol) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Volume of the real space unit cell
   real(wp), intent(in) :: vol

   !> Reciprocal term
   real(wp) :: rval

   rval = 4.0_wp*pi*(exp(-0.25_wp*gg*gg/(alpha**2))/(vol*gg*gg))

end function get_rec_term_3d


!> Returns the maximum value of a term in the 2D reciprocal space part of the Ewald
!> summation for a given vector length (at z = 0).
pure function get_rec_term_2d(gg, alpha, area) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Area of the 2D real space unit cell
   real(wp), intent(in) :: area

   !> Reciprocal term
   real(wp) :: rval

   rval = twopi * erfc(0.5_wp * gg / alpha) / (area * gg)

end function get_rec_term_2d


!> Returns the maximum value of a term in the 1D reciprocal space part of the Ewald
!> summation for a given vector length (at rho = 0).
pure function get_rec_term_1d(gg, alpha, len1) result(rval)

   !> Length of the reciprocal space vector
   real(wp), intent(in) :: gg

   !> Parameter of the Ewald summation
   real(wp), intent(in) :: alpha

   !> Length of the 1D real space unit cell
   real(wp), intent(in) :: len1

   !> Reciprocal term
   real(wp) :: rval
   real(wp) :: u

   u = 0.25_wp * gg * gg / (alpha * alpha)
   rval = 2.0_wp * exp(-u) / (len1 * (u + 1.0_wp))

end function get_rec_term_1d


end module multicharge_ewald
