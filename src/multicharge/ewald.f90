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
!> Contains the 1D, 2D, and 3D Ewald summation implementations
!>
!> The 1D and 2D Ewald decompositions follow:
!> - F. Nestler and D. Potts, *Fast Ewald summation under 2d- and 1d-periodic
!>   boundary conditions based on NFFTs*, Proc. Sampling Theory and Applications
!>   (SampTA), 329-332 (2013).
!> - F. Nestler, M. Pippig, and D. Potts, *Fast Ewald summation based on NFFT with
!>   mixed periodicity*, J. Comput. Phys. **285**, 280-315 (2015),
!>   DOI: [10.1016/j.jcp.2014.12.052](https://dx.doi.org/10.1016/j.jcp.2014.12.052).
!> This implementation evaluates their 1D and 2D Ewald formulations by direct
!> lattice summation.

!> @file multicharge/ewald.f90
!> Provides a reexport of the Ewald summation implementations

!> Proxy module to reexport 1D, 2D, and 3D Ewald summation procedures
module multicharge_ewald
   use multicharge_ewald_cutoff, only : get_alpha, get_ewald_dir_trans, &
      & get_ewald_rec_trans
   use multicharge_ewald_direct, only : get_amat_dir, get_damat_dir, get_hess_dir
   use multicharge_ewald_pbc1d, only : get_amat_rec_1d, get_damat_rec_1d, &
      & get_hess_rec_1d
   use multicharge_ewald_pbc2d, only : get_amat_rec_2d, get_damat_rec_2d, &
      & get_hess_rec_2d
   use multicharge_ewald_pbc3d, only : get_rec_weights, get_amat_rec_3d, &
      & get_damat_rec_3d, get_hess_rec_3d
   use multicharge_ewald_special, only : erfc_exp_2d, fast_f012_1d, fast_k012_1d
   implicit none
   private

   public :: get_alpha, get_ewald_dir_trans, get_ewald_rec_trans
   public :: erfc_exp_2d, fast_f012_1d, fast_k012_1d
   public :: get_amat_dir, get_damat_dir, get_hess_dir
   public :: get_amat_rec_1d, get_damat_rec_1d, get_hess_rec_1d
   public :: get_amat_rec_2d, get_damat_rec_2d, get_hess_rec_2d
   public :: get_rec_weights, get_amat_rec_3d, get_damat_rec_3d, get_hess_rec_3d

end module multicharge_ewald
