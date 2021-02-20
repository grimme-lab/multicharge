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

module multicharge_ncoord
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   implicit none
   private

   public :: get_coordination_number, cut_coordination_number


   !> Steepness of counting function
   real(wp), parameter :: kcn = 7.5_wp


contains


!> Geometric fractional coordination number, supports error function counting.
subroutine get_coordination_number(mol, trans, cutoff, rcov, cn, dcndr, dcndL, cut)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out), optional :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out), optional :: dcndL(:, :, :)

   !> Cut coordination number
   real(wp), intent(in), optional :: cut

   if (present(dcndr) .and. present(dcndL)) then
      call ncoord_derf(mol, trans, cutoff, rcov, cn, dcndr, dcndL)
   else
      call ncoord_erf(mol, trans, cutoff, rcov, cn)
   end if

   if (present(cut)) then
      call cut_coordination_number(cut, cn, dcndr, dcndL)
   end if

end subroutine get_coordination_number


subroutine ncoord_erf(mol, trans, cutoff, rcov, cn)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r2, r1, rc, rij(3), countf, cutoff2

   cn(:) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do default(none) reduction(+:cn) &
   !$omp shared(mol, trans, cutoff2, rcov) &
   !$omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)

         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            rc = rcov(izp) + rcov(jzp)

            countf = erf_count(kcn, r1, rc)

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

         end do
      end do
   end do

end subroutine ncoord_erf


subroutine ncoord_derf(mol, trans, cutoff, rcov, cn, dcndr, dcndL)

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   !> Lattice points
   real(wp), intent(in) :: trans(:, :)

   !> Real space cutoff
   real(wp), intent(in) :: cutoff

   !> Covalent radius
   real(wp), intent(in) :: rcov(:)

   !> Error function coordination number.
   real(wp), intent(out) :: cn(:)

   !> Derivative of the CN with respect to the Cartesian coordinates.
   real(wp), intent(out) :: dcndr(:, :, :)

   !> Derivative of the CN with respect to strain deformations.
   real(wp), intent(out) :: dcndL(:, :, :)

   integer :: iat, jat, izp, jzp, itr
   real(wp) :: r2, r1, rc, rij(3), countf, countd(3), sigma(3, 3), cutoff2

   cn(:) = 0.0_wp
   dcndr(:, :, :) = 0.0_wp
   dcndL(:, :, :) = 0.0_wp
   cutoff2 = cutoff**2

   !$omp parallel do default(none) reduction(+:cn, dcndr, dcndL) &
   !$omp shared(mol, trans, cutoff2, rcov) &
   !$omp private(jat, itr, izp, jzp, r2, rij, r1, rc, countf, countd, sigma)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat
         jzp = mol%id(jat)

         do itr = 1, size(trans, dim=2)
            rij = mol%xyz(:, iat) - (mol%xyz(:, jat) + trans(:, itr))
            r2 = sum(rij**2)
            if (r2 > cutoff2 .or. r2 < 1.0e-12_wp) cycle
            r1 = sqrt(r2)

            rc = rcov(izp) + rcov(jzp)

            countf = erf_count(kcn, r1, rc)
            countd = derf_count(kcn, r1, rc) * rij/r1

            cn(iat) = cn(iat) + countf
            if (iat /= jat) then
               cn(jat) = cn(jat) + countf
            end if

            dcndr(:, iat, iat) = dcndr(:, iat, iat) + countd
            dcndr(:, jat, jat) = dcndr(:, jat, jat) - countd
            dcndr(:, iat, jat) = dcndr(:, iat, jat) + countd
            dcndr(:, jat, iat) = dcndr(:, jat, iat) - countd

            sigma = spread(countd, 1, 3) * spread(rij, 2, 3)

            dcndL(:, :, iat) = dcndL(:, :, iat) + sigma
            if (iat /= jat) then
               dcndL(:, :, jat) = dcndL(:, :, jat) + sigma
            end if

         end do
      end do
   end do

end subroutine ncoord_derf


!> Error function counting function for coordination number contributions.
pure function erf_count(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp) :: count

   count = 0.5_wp * (1.0_wp + erf(-k*(r-r0)/r0))

end function erf_count


!> Derivative of the counting function w.r.t. the distance.
pure function derf_count(k, r, r0) result(count)

   !> Steepness of the counting function.
   real(wp), intent(in) :: k

   !> Current distance.
   real(wp), intent(in) :: r

   !> Cutoff radius.
   real(wp), intent(in) :: r0

   real(wp), parameter :: sqrtpi = sqrt(pi)

   real(wp) :: count

   count = -k/sqrtpi/r0*exp(-k**2*(r-r0)**2/r0**2)

end function derf_count


!> Cutoff function for large coordination numbers
pure subroutine cut_coordination_number(cn_max, cn, dcndr, dcndL)

   !> Maximum CN (not strictly obeyed)
   real(wp), intent(in) :: cn_max

   !> On input coordination number, on output modified CN
   real(wp), intent(inout) :: cn(:)

   !> On input derivative of CN w.r.t. cartesian coordinates,
   !> on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndr(:, :, :)

   !> On input derivative of CN w.r.t. strain deformation,
   !> on output derivative of modified CN
   real(wp), intent(inout), optional :: dcndL(:, :, :)

   real(wp) :: dcnpdcn
   integer  :: iat

   if (present(dcndL)) then
      do iat = 1, size(cn)
         dcnpdcn = dlog_cn_cut(cn(iat), cn_max)
         dcndL(:, :, iat) = dcnpdcn*dcndL(:, :, iat)
      enddo
   endif

   if (present(dcndr)) then
      do iat = 1, size(cn)
         dcnpdcn = dlog_cn_cut(cn(iat), cn_max)
         dcndr(:, :, iat) = dcnpdcn*dcndr(:, :, iat)
      enddo
   endif

   do iat = 1, size(cn)
      cn(iat) = log_cn_cut(cn(iat), cn_max)
   enddo

end subroutine cut_coordination_number

elemental function log_cn_cut(cn, cnmax) result(cnp)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cnmax
   real(wp) :: cnp
   cnp = log(1.0_wp + exp(cnmax)) - log(1.0_wp + exp(cnmax - cn))
end function log_cn_cut

elemental function dlog_cn_cut(cn, cnmax) result(dcnpdcn)
   real(wp), intent(in) :: cn
   real(wp), intent(in) :: cnmax
   real(wp) :: dcnpdcn
   dcnpdcn = exp(cnmax)/(exp(cnmax) + exp(cn))
end function dlog_cn_cut


end module multicharge_ncoord
