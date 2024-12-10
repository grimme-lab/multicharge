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

module multicharge_wignerseitz
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_cutoff, only : get_lattice_points
   implicit none
   private

   public :: wignerseitz_cell_type, new_wignerseitz_cell

   type :: wignerseitz_cell_type
      integer :: nimg_max
      integer, allocatable :: nimg(:, :)
      integer, allocatable :: tridx(:, :, :)
      real(wp), allocatable :: trans(:, :)
   end type wignerseitz_cell_type


   !> Small cutoff threshold to create only closest cells
   real(wp), parameter :: thr = sqrt(epsilon(0.0_wp))

   !> Tolerance to consider equivalent images
   real(wp), parameter :: tol = 0.01_wp


contains


subroutine new_wignerseitz_cell(self, mol)

   !> Wigner-Seitz cell instance
   type(wignerseitz_cell_type), intent(out) :: self

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   integer :: iat, jat, ntr, nimg
   integer, allocatable :: tridx(:)
   real(wp) :: vec(3)
   real(wp), allocatable :: trans(:, :)

   call get_lattice_points(mol%periodic, mol%lattice, thr, trans)
   ntr = size(trans, 2)
   allocate(self%nimg(mol%nat, mol%nat), self%tridx(ntr, mol%nat, mol%nat), &
      & tridx(ntr))

   self%nimg_max = 0
   !$omp parallel do default(none) schedule(runtime) collapse(2) &
   !$omp shared(mol, trans, self) private(iat, jat, vec, nimg, tridx)
   do iat = 1, mol%nat
      do jat = 1, mol%nat
         vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat)
         call get_pairs(nimg, trans, vec, tridx)
         self%nimg(jat, iat) = nimg
         self%tridx(:, jat, iat) = tridx
         self%nimg_max = max(nimg, self%nimg_max)
      end do
   end do

   call move_alloc(trans, self%trans)
   
end subroutine new_wignerseitz_cell


subroutine get_pairs(iws, trans, rij, list)
   integer, intent(out) :: iws
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: trans(:, :)
   integer, intent(out) :: list(:)

   logical :: mask(size(list))
   real(wp) :: dist(size(list)), vec(3), r2
   integer :: itr, img, pos

   iws = 0
   img = 0
   list(:) = 0
   mask(:) = .true.

   do itr = 1, size(trans, 2)
      vec(:) = rij - trans(:, itr)
      r2 = vec(1)**2 + vec(2)**2 + vec(3)**2
      if (r2 < thr) cycle
      img = img + 1
      dist(img) = r2
   end do

   if (img == 0) return

   pos = minloc(dist(:img), dim=1)

   r2 = dist(pos)
   mask(pos) = .false.

   iws = 1
   list(iws) = pos
   if (img <= iws) return

   do
      pos = minloc(dist(:img), dim=1, mask=mask(:img))
      if (abs(dist(pos) - r2) > tol) exit
      mask(pos) = .false.
      iws = iws + 1
      list(iws) = pos
   end do

end subroutine get_pairs


end module multicharge_wignerseitz
