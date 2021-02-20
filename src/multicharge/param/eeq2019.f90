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

!> Electronegativity equilibration charge model published in
!>
!> E. Caldeweyher, S. Ehlert, A. Hansen, H. Neugebauer, S. Spicher, C. Bannwarth
!> and S. Grimme, *J. Chem. Phys.*, **2019**, 150, 154122.
!> DOI: [10.1063/1.5090222](https://dx.doi.org/10.1063/1.5090222)
module multicharge_param_eeq2019
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_eeq_chi, get_eeq_eta, get_eeq_rad, get_eeq_kcn


   !> Element-specific electronegativity for the electronegativity equilibration charges.
   interface get_eeq_chi
      module procedure get_eeq_chi_sym
      module procedure get_eeq_chi_num
   end interface get_eeq_chi

   !> Element-specific chemical hardnesses for the electronegativity equilibration charges
   interface get_eeq_eta
      module procedure :: get_eeq_eta_sym
      module procedure :: get_eeq_eta_num
   end interface get_eeq_eta

   !> Element-specific CN scaling constant for the electronegativity equilibration charges.
   interface get_eeq_kcn
      module procedure :: get_eeq_kcn_sym
      module procedure :: get_eeq_kcn_num
   end interface get_eeq_kcn

   !> Element-specific charge widths for the electronegativity equilibration charges.
   interface get_eeq_rad
      module procedure :: get_eeq_rad_sym
      module procedure :: get_eeq_rad_num
   end interface get_eeq_rad


   !> Maximum atomic number allowed in EEQ calculations
   integer, parameter :: max_elem = 86


   !> Element-specific electronegativity for the electronegativity equilibration charges.
   real(wp), parameter :: eeq_chi(max_elem) = [&
      & 1.23695041_wp, 1.26590957_wp, 0.54341808_wp, 0.99666991_wp, 1.26691604_wp, &
      & 1.40028282_wp, 1.55819364_wp, 1.56866440_wp, 1.57540015_wp, 1.15056627_wp, &
      & 0.55936220_wp, 0.72373742_wp, 1.12910844_wp, 1.12306840_wp, 1.52672442_wp, &
      & 1.40768172_wp, 1.48154584_wp, 1.31062963_wp, 0.40374140_wp, 0.75442607_wp, &
      & 0.76482096_wp, 0.98457281_wp, 0.96702598_wp, 1.05266584_wp, 0.93274875_wp, &
      & 1.04025281_wp, 0.92738624_wp, 1.07419210_wp, 1.07900668_wp, 1.04712861_wp, &
      & 1.15018618_wp, 1.15388455_wp, 1.36313743_wp, 1.36485106_wp, 1.39801837_wp, &
      & 1.18695346_wp, 0.36273870_wp, 0.58797255_wp, 0.71961946_wp, 0.96158233_wp, &
      & 0.89585296_wp, 0.81360499_wp, 1.00794665_wp, 0.92613682_wp, 1.09152285_wp, &
      & 1.14907070_wp, 1.13508911_wp, 1.08853785_wp, 1.11005982_wp, 1.12452195_wp, &
      & 1.21642129_wp, 1.36507125_wp, 1.40340000_wp, 1.16653482_wp, 0.34125098_wp, &
      & 0.58884173_wp, 0.68441115_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
      & 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
      & 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, 0.56999999_wp, &
      & 0.56999999_wp, 0.87936784_wp, 1.02761808_wp, 0.93297476_wp, 1.10172128_wp, &
      & 0.97350071_wp, 1.16695666_wp, 1.23997927_wp, 1.18464453_wp, 1.14191734_wp, &
      & 1.12334192_wp, 1.01485321_wp, 1.12950808_wp, 1.30804834_wp, 1.33689961_wp, &
      & 1.27465977_wp]

   !> Element-specific chemical hardnesses for the electronegativity equilibration charges.
   real(wp), parameter :: eeq_eta(max_elem) = [&
      &-0.35015861_wp, 1.04121227_wp, 0.09281243_wp, 0.09412380_wp, 0.26629137_wp, &
      & 0.19408787_wp, 0.05317918_wp, 0.03151644_wp, 0.32275132_wp, 1.30996037_wp, &
      & 0.24206510_wp, 0.04147733_wp, 0.11634126_wp, 0.13155266_wp, 0.15350650_wp, &
      & 0.15250997_wp, 0.17523529_wp, 0.28774450_wp, 0.42937314_wp, 0.01896455_wp, &
      & 0.07179178_wp,-0.01121381_wp,-0.03093370_wp, 0.02716319_wp,-0.01843812_wp, &
      &-0.15270393_wp,-0.09192645_wp,-0.13418723_wp,-0.09861139_wp, 0.18338109_wp, &
      & 0.08299615_wp, 0.11370033_wp, 0.19005278_wp, 0.10980677_wp, 0.12327841_wp, &
      & 0.25345554_wp, 0.58615231_wp, 0.16093861_wp, 0.04548530_wp,-0.02478645_wp, &
      & 0.01909943_wp, 0.01402541_wp,-0.03595279_wp, 0.01137752_wp,-0.03697213_wp, &
      & 0.08009416_wp, 0.02274892_wp, 0.12801822_wp,-0.02078702_wp, 0.05284319_wp, &
      & 0.07581190_wp, 0.09663758_wp, 0.09547417_wp, 0.07803344_wp, 0.64913257_wp, &
      & 0.15348654_wp, 0.05054344_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
      & 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
      & 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, 0.11000000_wp, &
      & 0.11000000_wp,-0.02786741_wp, 0.01057858_wp,-0.03892226_wp,-0.04574364_wp, &
      &-0.03874080_wp,-0.03782372_wp,-0.07046855_wp, 0.09546597_wp, 0.21953269_wp, &
      & 0.02522348_wp, 0.15263050_wp, 0.08042611_wp, 0.01878626_wp, 0.08715453_wp, &
      & 0.10500484_wp]

   !> Element-specific CN scaling constant for the electronegativity equilibration charges.
   real(wp), parameter :: eeq_kcn(max_elem) = [&
      & 0.04916110_wp, 0.10937243_wp,-0.12349591_wp,-0.02665108_wp,-0.02631658_wp, &
      & 0.06005196_wp, 0.09279548_wp, 0.11689703_wp, 0.15704746_wp, 0.07987901_wp, &
      &-0.10002962_wp,-0.07712863_wp,-0.02170561_wp,-0.04964052_wp, 0.14250599_wp, &
      & 0.07126660_wp, 0.13682750_wp, 0.14877121_wp,-0.10219289_wp,-0.08979338_wp, &
      &-0.08273597_wp,-0.01754829_wp,-0.02765460_wp,-0.02558926_wp,-0.08010286_wp, &
      &-0.04163215_wp,-0.09369631_wp,-0.03774117_wp,-0.05759708_wp, 0.02431998_wp, &
      &-0.01056270_wp,-0.02692862_wp, 0.07657769_wp, 0.06561608_wp, 0.08006749_wp, &
      & 0.14139200_wp,-0.05351029_wp,-0.06701705_wp,-0.07377246_wp,-0.02927768_wp, &
      &-0.03867291_wp,-0.06929825_wp,-0.04485293_wp,-0.04800824_wp,-0.01484022_wp, &
      & 0.07917502_wp, 0.06619243_wp, 0.02434095_wp,-0.01505548_wp,-0.03030768_wp, &
      & 0.01418235_wp, 0.08953411_wp, 0.08967527_wp, 0.07277771_wp,-0.02129476_wp, &
      &-0.06188828_wp,-0.06568203_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
      &-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
      &-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp,-0.11000000_wp, &
      &-0.11000000_wp,-0.03585873_wp,-0.03132400_wp,-0.05902379_wp,-0.02827592_wp, &
      &-0.07606260_wp,-0.02123839_wp, 0.03814822_wp, 0.02146834_wp, 0.01580538_wp, &
      &-0.00894298_wp,-0.05864876_wp,-0.01817842_wp, 0.07721851_wp, 0.07936083_wp, &
      & 0.05849285_wp]

   !> Element-specific charge widths for the electronegativity equilibration charges.
   real(wp), parameter :: eeq_rad(max_elem) = [&
      & 0.55159092_wp, 0.66205886_wp, 0.90529132_wp, 1.51710827_wp, 2.86070364_wp, &
      & 1.88862966_wp, 1.32250290_wp, 1.23166285_wp, 1.77503721_wp, 1.11955204_wp, &
      & 1.28263182_wp, 1.22344336_wp, 1.70936266_wp, 1.54075036_wp, 1.38200579_wp, &
      & 2.18849322_wp, 1.36779065_wp, 1.27039703_wp, 1.64466502_wp, 1.58859404_wp, &
      & 1.65357953_wp, 1.50021521_wp, 1.30104175_wp, 1.46301827_wp, 1.32928147_wp, &
      & 1.02766713_wp, 1.02291377_wp, 0.94343886_wp, 1.14881311_wp, 1.47080755_wp, &
      & 1.76901636_wp, 1.98724061_wp, 2.41244711_wp, 2.26739524_wp, 2.95378999_wp, &
      & 1.20807752_wp, 1.65941046_wp, 1.62733880_wp, 1.61344972_wp, 1.63220728_wp, &
      & 1.60899928_wp, 1.43501286_wp, 1.54559205_wp, 1.32663678_wp, 1.37644152_wp, &
      & 1.36051851_wp, 1.23395526_wp, 1.65734544_wp, 1.53895240_wp, 1.97542736_wp, &
      & 1.97636542_wp, 2.05432381_wp, 3.80138135_wp, 1.43893803_wp, 1.75505957_wp, &
      & 1.59815118_wp, 1.76401732_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
      & 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
      & 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, 1.63999999_wp, &
      & 1.63999999_wp, 1.47055223_wp, 1.81127084_wp, 1.40189963_wp, 1.54015481_wp, &
      & 1.33721475_wp, 1.57165422_wp, 1.04815857_wp, 1.78342098_wp, 2.79106396_wp, &
      & 1.78160840_wp, 2.47588882_wp, 2.37670734_wp, 1.76613217_wp, 2.66172302_wp, &
      & 2.82773085_wp]

contains



!> Get electronegativity for species with a given symbol
elemental function get_eeq_chi_sym(symbol) result(chi)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> electronegativity
   real(wp) :: chi

   chi = get_eeq_chi(to_number(symbol))

end function get_eeq_chi_sym


!> Get electronegativity for species with a given atomic number
elemental function get_eeq_chi_num(number) result(chi)

   !> Atomic number
   integer, intent(in) :: number

   !> electronegativity
   real(wp) :: chi

   if (number > 0 .and. number <= size(eeq_chi, dim=1)) then
      chi = eeq_chi(number)
   else
      chi = -1.0_wp
   end if

end function get_eeq_chi_num


!> Get hardness for species with a given symbol
elemental function get_eeq_eta_sym(symbol) result(eta)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> hardness
   real(wp) :: eta

   eta = get_eeq_eta(to_number(symbol))

end function get_eeq_eta_sym


!> Get hardness for species with a given atomic number
elemental function get_eeq_eta_num(number) result(eta)

   !> Atomic number
   integer, intent(in) :: number

   !> hardness
   real(wp) :: eta

   if (number > 0 .and. number <= size(eeq_eta, dim=1)) then
      eta = eeq_eta(number)
   else
      eta = -1.0_wp
   end if

end function get_eeq_eta_num


!> Get CN scaling for species with a given symbol
elemental function get_eeq_kcn_sym(symbol) result(kcn)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> CN scaling
   real(wp) :: kcn

   kcn = get_eeq_kcn(to_number(symbol))

end function get_eeq_kcn_sym


!> Get CN scaling for species with a given atomic number
elemental function get_eeq_kcn_num(number) result(kcn)

   !> Atomic number
   integer, intent(in) :: number

   !> CN scaling
   real(wp) :: kcn

   if (number > 0 .and. number <= size(eeq_kcn, dim=1)) then
      kcn = eeq_kcn(number)
   else
      kcn = -1.0_wp
   end if

end function get_eeq_kcn_num


!> Get charge width for species with a given symbol
elemental function get_eeq_rad_sym(symbol) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> charge width
   real(wp) :: rad

   rad = get_eeq_rad(to_number(symbol))

end function get_eeq_rad_sym


!> Get charge width for species with a given atomic number
elemental function get_eeq_rad_num(number) result(rad)

   !> Atomic number
   integer, intent(in) :: number

   !> Charge width
   real(wp) :: rad

   if (number > 0 .and. number <= size(eeq_rad, dim=1)) then
      rad = eeq_rad(number)
   else
      rad = -1.0_wp
   end if

end function get_eeq_rad_num


end module multicharge_param_eeq2019
