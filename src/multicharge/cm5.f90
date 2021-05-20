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

!> Implementation of the charge model 5
!>
!> A.V. Marenich, S.V. Jerome, C.J. Cramer, and D.G. Truhlar
!> J. Chem. Theory Comput. 2012, 8, 2, 527–541. DOI: 10.1021/ct200866d
module multicharge_cm5
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_convert, only : aatoau
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: cm5_type, new_cm5

   !> Charge model 5 data
   type :: cm5_type
      !> Global model parameter
      real(wp) :: alpha
      !> Atomic radii
      real(wp), allocatable :: rad(:)
      !> Pair parameters
      real(wp), allocatable :: par(:, :)
   contains
      procedure :: get_correction
   end type cm5_type


  !> Get pair parameters for two species (Dzz')
  interface get_pair_parameter
    module procedure :: get_pair_parameter_symbol
    module procedure :: get_pair_parameter_number
  end interface get_pair_parameter


  ! Charge model 5 atomwise parameters (Dz in paper)
  real(wp), parameter :: dz_parameter(118) = [ &
      & 0.0056_wp,-0.1543_wp, 0.0000_wp, 0.0333_wp,-0.1030_wp,-0.0446_wp, &
      &-0.1072_wp,-0.0802_wp,-0.0629_wp,-0.1088_wp, 0.0184_wp, 0.0000_wp, &
      &-0.0726_wp,-0.0790_wp,-0.0756_wp,-0.0565_wp,-0.0444_wp,-0.0767_wp, &
      & 0.0130_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      &-0.0512_wp,-0.0557_wp,-0.0533_wp,-0.0399_wp,-0.0313_wp,-0.0541_wp, &
      & 0.0092_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      &-0.0361_wp,-0.0393_wp,-0.0376_wp,-0.0281_wp,-0.0220_wp,-0.0381_wp, &
      & 0.0065_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp,-0.0255_wp,-0.0277_wp,-0.0265_wp,-0.0198_wp, &
      &-0.0155_wp,-0.0269_wp, 0.0046_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp, &
      & 0.0000_wp, 0.0000_wp, 0.0000_wp, 0.0000_wp,-0.0179_wp,-0.0195_wp, &
      &-0.0187_wp,-0.0140_wp,-0.0110_wp,-0.0189_wp]

   !> Global model parameter
   real(wp), parameter :: alpha_default = 2.474_wp/aatoau

contains

!> Constructor for charge model 5
subroutine new_cm5(self, mol, rad)
   !> Instance of the charge model
   type(cm5_type), intent(out) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Atomic radii
   real(wp), intent(in) :: rad(:)

   integer :: isp, jsp

   self%alpha = alpha_default
   self%rad = rad

    allocate(self%par(mol%nid, mol%nid))
    do isp = 1, mol%nid
      do jsp = 1, mol%nid
        self%par(isp, jsp) = get_pair_parameter(mol%num(isp), mol%num(jsp))
      end do
    end do
end subroutine new_cm5

!> Evaluate charge model correction
subroutine get_correction(self, mol, trans, cutoff, cm5, dcm5dr, dcm5dL)
   !> Instance of the charge model
   class(cm5_type), intent(in) :: self
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Real space cutoff
   real(wp), intent(in) :: cutoff
   !> Lattice translations
   real(wp), intent(in) :: trans(:, :)
   !> Charge model 5 correction
   real(wp), intent(out) :: cm5(:)
   !> Derivative of correction w.r.t. cartesian displacements
   real(wp), intent(out), optional :: dcm5dr(:, :, :)
   !> Derivative of correction w.r.t. strain deformations
   real(wp), intent(out), optional :: dcm5dL(:, :, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: vec(3), dist, pij, pji, dEr, dGr(3), dSr(3, 3)

   cm5(:) = 0.0_wp

   if (present(dcm5dr) .and. present(dcm5dL)) then
      dcm5dr(:, :, :) = 0.0_wp
      dcm5dL(:, :, :) = 0.0_wp

      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:cm5, dcm5dr, dcm5dL) shared(mol, self, cutoff, trans) &
      !$omp private(iat, jat, izp, jzp, itr, vec, dist, pij, pji, dEr, dGr, dSr)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            if (izp == jzp) cycle
            do itr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               dist = norm2(vec)
               if (dist > cutoff) cycle
               pij = self%par(izp, jzp)
               pji = self%par(jzp, izp)

               dEr = exp(-self%alpha*(dist-self%rad(izp)-self%rad(jzp)))
               dGr = dEr * self%alpha * vec/dist
               dSr = spread(dGr, 1, 3) * spread(vec, 2, 3)

               cm5(iat) = cm5(iat) + dEr * pij
               cm5(jat) = cm5(jat) + dEr * pji

               dcm5dr(:, iat, iat) = dcm5dr(:, iat, iat) - dGr * pij
               dcm5dr(:, jat, jat) = dcm5dr(:, jat, jat) + dGr * pji
               dcm5dr(:, iat, jat) = dcm5dr(:, iat, jat) - dGr * pji
               dcm5dr(:, jat, iat) = dcm5dr(:, jat, iat) + dGr * pij

               dcm5dL(:, :, iat) = dcm5dL(:, :, iat) - dSr * pij
               dcm5dL(:, :, jat) = dcm5dL(:, :, jat) - dSr * pji
            end do
         end do
      end do
   else
      !$omp parallel do default(none) schedule(runtime) &
      !$omp reduction(+:cm5) shared(mol, self, cutoff, trans) &
      !$omp private(iat, jat, izp, jzp, itr, vec, dist, pij, pji, dEr)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         do jat = 1, iat - 1
            jzp = mol%id(jat)
            if (izp == jzp) cycle
            do itr = 1, size(trans, 2)
               vec(:) = mol%xyz(:, iat) - mol%xyz(:, jat) - trans(:, itr)
               dist = norm2(vec)
               if (dist > cutoff) cycle
               pij = self%par(izp, jzp)
               pji = self%par(jzp, izp)

               dEr = exp(-self%alpha*(dist-self%rad(izp)-self%rad(jzp)))

               cm5(iat) = cm5(iat) + dEr * pij
               cm5(jat) = cm5(jat) + dEr * pji
            end do
         end do
      end do
   end if
end subroutine get_correction


!> Get pair parameter (Dzz') for species with a given symbols
elemental function get_pair_parameter_symbol(symbol1, symbol2) result(par)
   !> Element symbol
   character(len=*), intent(in) :: symbol1
   !> Element symbol
   character(len=*), intent(in) :: symbol2
   !> Pair parameter
   real(wp) :: par

   par = get_pair_parameter(to_number(symbol1), to_number(symbol2))

end function get_pair_parameter_symbol


!> Get pair parameter (Dzz') for species with a given atomic numbers
elemental function get_pair_parameter_number(number1, number2) result(par)
   !> Atomic number
   integer, intent(in) :: number1
   !> Atomic number
   integer, intent(in) :: number2
   !> Pair parameter
   real(wp) :: par

   if (number1 > 0 .and. number1 <= size(dz_parameter, dim=1) .and. &
      & number2 > 0 .and. number2 <= size(dz_parameter, dim=1)) then
      if (number1 == 1 .and. number2 == 6) then
         par = 0.0502_wp
      else if (number1 == 6 .and. number2 == 1) then
         par = -0.0502_wp
      else if (number1 == 1 .and. number2 == 7) then
         par = 0.1747_wp
      else if (number1 == 7 .and. number2 == 1) then
         par = -0.1747_wp
      else if (number1 == 1 .and. number2 == 8) then
         par = 0.1671_wp
      else if (number1 == 8 .and. number2 == 1) then
         par = -0.1671_wp
      else if (number1 == 6 .and. number2 == 7) then
         par = 0.0556_wp
      else if (number1 == 7 .and. number2 == 6) then
         par = -0.0556_wp
      else if (number1 == 6 .and. number2 == 8) then
         par = 0.0234_wp
      else if (number1 == 8 .and. number2 == 6) then
         par = -0.0234_wp
      else if (number1 == 7 .and. number2 == 8) then
         par = -0.0346_wp
      else if (number1 == 8 .and. number2 == 7) then
         par = 0.0346_wp
      else
         par = dz_parameter(number1) - dz_parameter(number2)
      end if
   else
      par = 0.0_wp
   end if

end function get_pair_parameter_number

end module multicharge_cm5
