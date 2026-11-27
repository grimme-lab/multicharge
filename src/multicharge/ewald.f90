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

!> @dir multicharge/ewald
!> Contains the 1D, 2D, and 3D periodic Ewald summation implementations
!>
!> The 1D and 2D Ewald decompositions follow:
!> - F. Nestler and D. Potts, *Fast Ewald summation under 2d- and 1d-periodic
!>   boundary conditions based on NFFTs*, Proc. Sampling Theory and Applications
!>   (SampTA), 329-332 (2013).
!> - F. Nestler, M. Pippig, and D. Potts, *Fast Ewald summation based on NFFT with
!>   mixed periodicity*, J. Comput. Phys. **285**, 280-315 (2015),
!>   DOI: [10.1016/j.jcp.2014.12.052](https://dx.doi.org/10.1016/j.jcp.2014.12.052).
!> This implementation evaluates their 1D and 2D Ewald formulations by direct
!> lattice summation within tolerance-derived cutoff spheres.

!> @file multicharge/ewald.f90
!> Provides the periodic Ewald summation interface

!> Ewald splitting parameter, lattice translations and reciprocal-space kernels
!> for 1D, 2D, and 3D periodic boundary conditions
module multicharge_ewald
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_wignerseitz, only : wignerseitz_cell
   use multicharge_ewald_cutoff, only : build_ewald_trans, compute_alpha, &
      & get_alpha, get_ewald_geometry
   use multicharge_ewald_pbc1d, only : get_amat_rec_1d, get_damat_rec_1d, &
      & get_hess_rec_1d
   use multicharge_ewald_pbc2d, only : get_amat_rec_2d, get_damat_rec_2d, &
      & get_hess_rec_2d
   use multicharge_ewald_pbc3d, only : get_rec_weights, get_amat_rec_3d, &
      & get_damat_rec_3d, get_hess_rec_3d
   implicit none
   private

   public :: get_alpha, get_ewald_setup
   public :: get_amat_rec_1d, get_damat_rec_1d, get_hess_rec_1d
   public :: get_amat_rec_2d, get_damat_rec_2d, get_hess_rec_2d
   public :: get_rec_weights, get_amat_rec_3d, get_damat_rec_3d, get_hess_rec_3d


contains


!> Prepare periodic Ewald geometry, splitting parameter, and converged lattice translations
subroutine get_ewald_setup(mol, wsc, gam, ndim, measure, axis_vec, alpha, &
   & dtrans, rtrans, error, alpha_in, tol_scale)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Wigner-Seitz cell images of all pairs
   type(wignerseitz_cell), intent(in) :: wsc

   !> Smallest Gaussian charge exponent of all pairs
   real(wp), intent(in) :: gam

   !> Number of active periodic dimensions (1, 2, or 3)
   integer, intent(out) :: ndim

   !> Generalized cell measure: length (1D), area (2D), or volume (3D)
   real(wp), intent(out) :: measure

   !> Unit wire direction (1D) or unit slab normal (2D); zero in 3D
   real(wp), intent(out) :: axis_vec(3)

   !> Ewald splitting parameter
   real(wp), intent(out) :: alpha

   !> Direct-space lattice translations (3, ndir)
   real(wp), allocatable, intent(out) :: dtrans(:, :)

   !> Half-space reciprocal lattice vectors (3, nrec)
   real(wp), allocatable, intent(out) :: rtrans(:, :)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Optional explicit Ewald splitting parameter override
   real(wp), intent(in), optional :: alpha_in

   !> Optional multiplicative scale on the per-term cutoff tolerance
   real(wp), intent(in), optional :: tol_scale

   logical :: mask(3), has_alpha_in
   integer :: iat, jat, img
   real(wp) :: red_lat(3, 3), rec_lat(3, 3), vec(3), v_par(3), r2, rmax2, alpha_override

   has_alpha_in = present(alpha_in)
   if (has_alpha_in) alpha_override = alpha_in
   alpha = 0.0_wp
   call get_ewald_geometry(mol%lattice, mask, ndim, measure, axis_vec, &
      & red_lat, rec_lat, error, mol%periodic)
   if (allocated(error)) return

   if (has_alpha_in) then
      alpha = alpha_override
   else
      alpha = compute_alpha(ndim, measure)
   end if

   ! Longest pair displacement projected onto the periodic subspace.
   ! Perpendicular separation is orthogonal to every lattice translation T,
   ! so |v + T| <= rcut implies |T| <= rcut + |v_par|.
   rmax2 = 0.0_wp
   !$omp parallel do default(none) schedule(runtime) &
   !$omp shared(mol, wsc, ndim, axis_vec) private(iat, jat, img, vec, v_par, r2) &
   !$omp reduction(max:rmax2)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, jat) - mol%xyz(:, iat) &
               & + wsc%trans(:, wsc%tridx(img, jat, iat))
            select case (ndim)
            case (1)
               r2 = dot_product(vec, axis_vec)**2
            case (2)
               v_par = vec - dot_product(vec, axis_vec) * axis_vec
               r2 = dot_product(v_par, v_par)
            case default
               r2 = dot_product(vec, vec)
            end select
            rmax2 = max(rmax2, r2)
         end do
      end do
   end do
   !$omp end parallel do

   call build_ewald_trans(mask, ndim, measure, red_lat, rec_lat, alpha, gam, &
      & sqrt(rmax2), dtrans, rtrans, error, tol_scale=tol_scale)

end subroutine get_ewald_setup


end module multicharge_ewald
