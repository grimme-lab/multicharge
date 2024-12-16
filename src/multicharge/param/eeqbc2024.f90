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

!> Bond capacitor electronegativity equilibration charge model published in
!>
!> T. Froitzheim, M. Müller, A. Hansen, and S. Grimme
!> , *J. Chem. Phys.*, in preparation.
module multicharge_param_eeqbc2024
   use mctc_env, only : wp
   use mctc_io_symbols, only : to_number
   implicit none
   private

   public :: get_eeqbc_chi, get_eeqbc_eta, get_eeqbc_rad, get_eeqbc_kcnchi, &
      & get_eeqbc_kqchi, get_eeqbc_kqeta, get_eeqbc_cap, get_eeqbc_cov_radii, &
      & get_eeqbc_avg_cn


   !> Element-specific electronegativity for the EEQ_BC charges.
   interface get_eeqbc_chi
      module procedure get_eeqbc_chi_sym
      module procedure get_eeqbc_chi_num
   end interface get_eeqbc_chi

   !> Element-specific chemical hardnesses for the EEQ_BC charges.
   interface get_eeqbc_eta
      module procedure :: get_eeqbc_eta_sym
      module procedure :: get_eeqbc_eta_num
   end interface get_eeqbc_eta

   !> Element-specific charge widths for the EEQ_BC charges.
   interface get_eeqbc_rad
      module procedure :: get_eeqbc_rad_sym
      module procedure :: get_eeqbc_rad_num
   end interface get_eeqbc_rad

   !> Element-specific CN scaling of the electronegativity for the EEQ_BC charges.
   interface get_eeqbc_kcnchi
      module procedure :: get_eeqbc_kcnchi_sym
      module procedure :: get_eeqbc_kcnchi_num
   end interface get_eeqbc_kcnchi

   !> Element-specific local q scaling of the electronegativity for the EEQ_BC charges.
   interface get_eeqbc_kqchi
      module procedure :: get_eeqbc_kqchi_sym
      module procedure :: get_eeqbc_kqchi_num
   end interface get_eeqbc_kqchi

   !> Element-specific local q scaling of the chemical hardness for the EEQ_BC charges.
   interface get_eeqbc_kqeta
      module procedure :: get_eeqbc_kqeta_sym
      module procedure :: get_eeqbc_kqeta_num
   end interface get_eeqbc_kqeta

   !> Element-specific bond capacitance for the EEQ_BC charges.
   interface get_eeqbc_cap
      module procedure :: get_eeqbc_cap_sym
      module procedure :: get_eeqbc_cap_num
   end interface get_eeqbc_cap

   !> Element-specific covalent radii for the CN for the EEQ_BC charges.
   interface get_eeqbc_cov_radii
      module procedure :: get_eeqbc_cov_radii_sym
      module procedure :: get_eeqbc_cov_radii_num
   end interface get_eeqbc_cov_radii

   !> Element-specific average CN for the EEQ_BC charges.
   interface get_eeqbc_avg_cn
      module procedure :: get_eeqbc_avg_cn_sym
      module procedure :: get_eeqbc_avg_cn_num
   end interface get_eeqbc_avg_cn

   !> Maximum atomic number allowed in EEQ_BC calculations
   integer, parameter :: max_elem = 103


   !> Element-specific electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_chi(max_elem) = [&
      &  1.7500687479_wp,  0.7992983109_wp,  0.8817302909_wp,  1.2122559922_wp, & !1-4
      &  1.4042606312_wp,  1.7373300176_wp,  1.9224220861_wp,  2.0295674708_wp, & !5-8
      &  2.0914017724_wp,  0.2783743672_wp,  0.7909141712_wp,  0.9333749946_wp, & !9-12
      &  1.1280735350_wp,  1.3504642320_wp,  1.7084529806_wp,  1.9657999323_wp, & !13-16
      &  1.8796814465_wp,  0.8120477849_wp,  0.6229777212_wp,  0.8955669337_wp, & !17-20
      &  0.8887941055_wp,  0.9249293933_wp,  0.8910306356_wp,  0.8730274586_wp, & !21-24
      &  1.0692963783_wp,  1.1430792497_wp,  1.2352658732_wp,  1.2511161359_wp, & !25-28
      &  1.0995052580_wp,  1.0059572004_wp,  1.0390725738_wp,  1.2885924052_wp, & !29-32
      &  1.4638654613_wp,  1.7797597799_wp,  1.6400765990_wp,  0.8859889377_wp, & !33-36
      &  0.5142094052_wp,  0.8785352464_wp,  0.9716967887_wp,  0.8109573582_wp, & !37-40
      &  0.9361297862_wp,  0.9872048394_wp,  1.1290914832_wp,  1.0416755409_wp, & !41-44
      &  1.1579060755_wp,  1.1371382461_wp,  1.1490154759_wp,  1.0811257447_wp, & !45-48
      &  1.0201038561_wp,  1.2317318949_wp,  1.2546053590_wp,  1.6136334955_wp, & !49-52
      &  1.5949826440_wp,  0.8548714270_wp,  0.5116591821_wp,  0.8221154800_wp, & !53-56
      &  0.8992637384_wp,  0.7835477700_wp,  0.6865502434_wp,  0.6063027416_wp, & !57-60
      &  0.5428052646_wp,  0.4960578124_wp,  0.4660603851_wp,  0.4528129826_wp, & !61-64
      &  0.4563156049_wp,  0.4765682521_wp,  0.5135709241_wp,  0.5673236209_wp, & !65-68
      &  0.6378263425_wp,  0.7250790890_wp,  0.8290818603_wp,  0.8697550816_wp, & !69-72
      &  1.0442533196_wp,  1.1429836348_wp,  1.1622493128_wp,  1.2650483683_wp, & !73-76
      &  1.2650500943_wp,  1.3607929134_wp,  1.3186071563_wp,  1.0545750683_wp, & !77-80
      &  0.9074468503_wp,  1.0892548243_wp,  1.1983441731_wp,  1.3955974910_wp, & !81-84
      &  1.6266506350_wp,  0.9802627692_wp,  0.4952498716_wp,  0.7903508991_wp, & !85-88
      &  0.7482689572_wp,  0.8666000614_wp,  0.8153381406_wp,  0.7700731721_wp, & !89-92
      &  0.7308051560_wp,  0.6975340922_wp,  0.6702599807_wp,  0.6489828216_wp, & !93-96
      &  0.6337026148_wp,  0.6244193604_wp,  0.6211330583_wp,  0.6238437086_wp, & !97-100
      &  0.6325513112_wp,  0.6472558662_wp,  0.6679573735_wp] !101-103

   !> Element-specific chemical hardnesses for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_eta(max_elem) = [&
      &  0.3572813340_wp, 14.1713349136_wp, -0.0335574715_wp, -2.2617753890_wp, & !1-4
      & -2.9993990603_wp, -2.8456422314_wp, -2.2316836385_wp, -0.9048085573_wp, & !5-8
      & -3.3402942035_wp, 11.6677100883_wp,  0.0461110187_wp, -0.1623149426_wp, & !9-12
      & -0.1976009198_wp, -3.6156182254_wp, -4.8040123811_wp, -5.8989254120_wp, & !13-16
      & -1.7918672558_wp,  3.2077831067_wp,  0.4598658365_wp, -0.3196730368_wp, & !17-20
      & -0.0066012997_wp, -0.0650415781_wp,  0.0116105065_wp, -0.2020240365_wp, & !21-24
      & -0.0451985500_wp, -0.8983846024_wp, -0.5087624261_wp, -0.9360254729_wp, & !25-28
      & -0.3137611925_wp,  0.3714666864_wp, -0.5637510788_wp, -1.5811888792_wp, & !29-32
      & -2.5680164043_wp, -3.3791525742_wp, -0.9039263250_wp,  2.6191171553_wp, & !33-36
      &  0.4517188832_wp, -0.4737572247_wp, -0.3291918172_wp, -0.0641706161_wp, & !37-40
      & -0.4365721167_wp, -0.1388382729_wp,  0.0445179428_wp, -0.3077776724_wp, & !41-44
      & -0.1421769591_wp, -0.3718332953_wp, -0.9003899901_wp, -0.5034953355_wp, & !45-48
      & -0.3154724874_wp, -1.2061278491_wp, -1.0351395610_wp, -2.4727516433_wp, & !49-52
      & -0.5377076044_wp,  2.1647776210_wp,  0.3592585022_wp, -0.6373016543_wp, & !53-56
      & -0.1481956999_wp, -0.4595916155_wp, -0.6048435529_wp, -0.7208619618_wp, & !57-60
      & -0.8076468424_wp, -0.8651981945_wp, -0.8935160183_wp, -0.8926003136_wp, & !61-64
      & -0.8624510805_wp, -0.8030683191_wp, -0.7144520292_wp, -0.5966022109_wp, & !65-68
      & -0.4495188642_wp, -0.2732019891_wp, -0.0676515856_wp, -0.1339322663_wp, & !69-72
      & -0.7103642117_wp, -0.1700796179_wp, -0.1362891699_wp, -1.0705189016_wp, & !73-76
      & -0.8229572159_wp, -1.3207540081_wp, -2.0554362750_wp, -0.2654477885_wp, & !77-80
      & -0.0736143849_wp, -1.1221956034_wp, -0.1821999108_wp, -0.7727065022_wp, & !81-84
      & -0.4699768943_wp,  0.6377347433_wp,  0.4140010159_wp, -0.2353223377_wp, & !85-88
      & -0.1309097826_wp,  0.1881855179_wp,  0.2007222471_wp,  0.1912792246_wp, & !89-92
      &  0.1598564505_wp,  0.1064539248_wp,  0.0310716475_wp, -0.0662903814_wp, & !93-96
      & -0.1856321619_wp, -0.3269536941_wp, -0.4902549779_wp, -0.6755360133_wp, & !97-100
      & -0.8827968003_wp, -1.1120373389_wp, -1.3632576292_wp] !101-103

   !> Element-specific charge widths for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_rad(max_elem) = [&
      &  0.4537866920_wp,  0.8971879958_wp,  0.3987756594_wp,  0.2435934990_wp, & !1-4
      &  0.2119711703_wp,  0.2064066867_wp,  0.2398313485_wp,  0.3482853216_wp, & !5-8
      &  0.1479057386_wp,  1.4433940527_wp,  0.6317031456_wp,  0.7152255265_wp, & !9-12
      &  0.6920759433_wp,  0.1952261525_wp,  0.1478738486_wp,  0.1173410276_wp, & !13-16
      &  0.2188836429_wp,  0.7265491450_wp,  1.0062576628_wp,  0.6529550574_wp, & !17-20
      &  1.0787300626_wp,  1.0194369772_wp,  0.7673907688_wp,  0.8234907812_wp, & !21-24
      &  0.7956000862_wp,  0.4194926962_wp,  0.6577871621_wp,  0.4350022430_wp, & !25-28
      &  0.5436327263_wp,  1.2387687941_wp,  0.5125789654_wp,  0.3834386963_wp, & !29-32
      &  0.2781070074_wp,  0.2053677667_wp,  0.3191301456_wp,  3.4957602962_wp, & !33-36
      &  0.8847073217_wp,  0.6739335178_wp,  0.8092111775_wp,  0.8229663676_wp, & !37-40
      &  0.7341667740_wp,  0.8802988629_wp,  1.1234870897_wp,  0.5654595735_wp, & !41-44
      &  0.7749739189_wp,  0.6091511140_wp,  0.4788100227_wp,  0.6104947355_wp, & !45-48
      &  0.6518973596_wp,  0.4348284778_wp,  0.4885595700_wp,  0.2660054523_wp, & !49-52
      &  0.4274914591_wp,  2.3114324559_wp,  0.9734795056_wp,  0.6329900422_wp, & !53-56
      &  1.0109847900_wp,  0.6287499845_wp,  0.5401093486_wp,  0.4679527826_wp, & !57-60
      &  0.4122802864_wp,  0.3730918601_wp,  0.3503875036_wp,  0.3441672169_wp, & !61-64
      &  0.3544310001_wp,  0.3811788531_wp,  0.4244107759_wp,  0.4841267686_wp, & !65-68
      &  0.5603268311_wp,  0.6530109634_wp,  0.7621791656_wp,  1.0577606985_wp, & !69-72
      &  0.6844888492_wp,  0.9102124518_wp,  0.8550543040_wp,  0.4138761210_wp, & !73-76
      &  0.5593056202_wp,  0.3751752813_wp,  0.2949155601_wp,  0.6769971683_wp, & !77-80
      &  0.7124606732_wp,  0.4519163133_wp,  1.0405678353_wp,  0.6688421527_wp, & !81-84
      &  0.4838599292_wp,  0.9792188430_wp,  0.8793273061_wp,  0.8333325045_wp, & !85-88
      &  0.8202868436_wp,  1.7807640816_wp,  1.5641357264_wp,  1.3644976007_wp, & !89-92
      &  1.1818497047_wp,  1.0161920382_wp,  0.8675246014_wp,  0.7358473941_wp, & !93-96
      &  0.6211604164_wp,  0.5234636683_wp,  0.4427571498_wp,  0.3790408609_wp, & !97-100
      &  0.3323148016_wp,  0.3025789719_wp,  0.2898333718_wp] !101-103

   !> Element-specific CN scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kcnchi(max_elem) = [&
      &  1.3415783494_wp,  2.4226307746_wp,  0.0910702713_wp, -0.2802662922_wp, & !1-4
      & -0.0464303067_wp,  0.3049790613_wp,  0.5014914830_wp,  0.7131712513_wp, & !5-8
      &  1.5978006993_wp,  4.6934800245_wp, -0.2311835622_wp, -0.5722047540_wp, & !9-12
      & -0.1872404228_wp,  0.1355861183_wp,  0.5037598487_wp,  0.8257488249_wp, & !13-16
      &  1.5828922925_wp,  5.6324196990_wp, -1.3574661808_wp, -0.7114730764_wp, & !17-20
      & -0.8412840531_wp, -0.8100781799_wp, -0.7321477749_wp, -0.5690936866_wp, & !21-24
      & -0.7978421025_wp, -0.7081664947_wp, -0.5311094926_wp, -0.5561735098_wp, & !25-28
      &  0.1043470768_wp, -0.2459258932_wp, -0.2244771250_wp,  0.0378446029_wp, & !29-32
      &  0.2939641775_wp,  0.7336233202_wp,  1.1960377617_wp,  1.5974038323_wp, & !33-36
      & -0.5630850954_wp, -1.1059510466_wp, -0.7830773028_wp, -0.9114834757_wp, & !37-40
      & -0.4093603622_wp, -0.2717170095_wp, -0.4691579275_wp, -0.2257381361_wp, & !41-44
      & -0.1375984198_wp,  0.3330053570_wp,  0.0221109296_wp, -0.0920402467_wp, & !45-48
      & -0.3096506887_wp,  0.0088013637_wp,  0.0730363100_wp,  0.4356094483_wp, & !49-52
      &  1.0199146044_wp,  1.0092039203_wp, -0.7528024837_wp, -1.1365506475_wp, & !53-56
      & -0.9661197708_wp, -1.1514088354_wp, -1.1092964223_wp, -1.0718762355_wp, & !57-60
      & -1.0391482749_wp, -1.0111125406_wp, -0.9877690325_wp, -0.9691177507_wp, & !61-64
      & -0.9551586951_wp, -0.9458918658_wp, -0.9413172627_wp, -0.9414348859_wp, & !65-68
      & -0.9462447354_wp, -0.9557468111_wp, -0.9699411130_wp, -0.9467711075_wp, & !69-72
      & -0.5854657957_wp, -0.1956906192_wp, -0.3841246137_wp, -0.2184058724_wp, & !73-76
      & -0.2071244723_wp,  0.1769757167_wp,  0.5363613694_wp,  0.0342662426_wp, & !77-80
      & -0.5074824777_wp, -0.0048092213_wp, -0.0546120433_wp,  0.0560290491_wp, & !81-84
      &  0.8822097689_wp,  0.9546406691_wp, -1.8612818673_wp, -1.2559850201_wp, & !85-88
      & -0.8232940275_wp, -0.7432092987_wp, -0.9259469469_wp, -1.0588247895_wp, & !89-92
      & -1.1418428264_wp, -1.1750010577_wp, -1.1582994833_wp, -1.0917381033_wp, & !93-96
      & -0.9753169176_wp, -0.8090359263_wp, -0.5928951293_wp, -0.3268945267_wp, & !97-100
      & -0.0110341184_wp,  0.3546860955_wp,  0.7702661151_wp] !101-103

   !> Element-specific local q scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kqchi(max_elem) = [&
      &  0.7122604774_wp, -1.7351284097_wp,  3.0089829052_wp,  2.1166762050_wp, & !1-4
      &  1.5179774898_wp,  1.2180269092_wp,  1.0873609014_wp,  0.8994075937_wp, & !5-8
      &  0.1658248786_wp, -2.5747028940_wp,  3.1762170214_wp,  2.3987338612_wp, & !9-12
      &  2.2469063726_wp,  1.5639940746_wp,  1.2412557993_wp,  1.6283237163_wp, & !13-16
      &  1.5628790844_wp, -0.9249536928_wp,  3.0733040004_wp,  2.7596745507_wp, & !17-20
      &  2.9366708989_wp,  2.7004746183_wp,  2.2295030415_wp,  2.0304690076_wp, & !21-24
      &  1.9683561829_wp,  2.2302711526_wp,  1.8504904266_wp,  2.0575510119_wp, & !25-28
      &  2.2756603413_wp,  2.2094576537_wp,  2.1544064368_wp,  1.9327504630_wp, & !29-32
      &  1.4451438826_wp,  1.4813741556_wp,  2.0308095325_wp,  0.4032186085_wp, & !33-36
      &  3.6036894994_wp,  2.6513413398_wp,  2.6634586616_wp,  2.3940154835_wp, & !37-40
      &  2.3527262731_wp,  2.0735381213_wp,  1.7234564437_wp,  2.2302635382_wp, & !41-44
      &  2.1871313764_wp,  1.8061408427_wp,  1.9051691947_wp,  2.0424482278_wp, & !45-48
      &  2.8036578365_wp,  2.0783981020_wp,  2.0481231960_wp,  1.8544101088_wp, & !49-52
      &  2.1888387015_wp,  0.5779869189_wp,  3.2064625646_wp,  2.7406551784_wp, & !53-56
      &  2.5529621630_wp,  2.5391757608_wp,  2.4348800350_wp,  2.3484586230_wp, & !57-60
      &  2.2799115250_wp,  2.2292387408_wp,  2.1964402705_wp,  2.1815161141_wp, & !61-64
      &  2.1844662716_wp,  2.2052907430_wp,  2.2439895282_wp,  2.3005626274_wp, & !65-68
      &  2.3750100404_wp,  2.4673317674_wp,  2.5775278082_wp,  2.6463737671_wp, & !69-72
      &  2.3987259080_wp,  2.0862161326_wp,  1.8045334538_wp,  2.0382923920_wp, & !73-76
      &  1.6579982531_wp,  1.8353080915_wp,  1.8450710788_wp,  1.5696036105_wp, & !77-80
      &  2.8136219641_wp,  2.3784572290_wp,  1.9914691678_wp,  1.8625351100_wp, & !81-84
      &  2.1579257719_wp,  0.6206683275_wp,  3.5103871382_wp,  2.7327597379_wp, & !85-88
      &  2.7369312006_wp,  2.6004448612_wp,  2.7011486104_wp,  2.7879694953_wp, & !89-92
      &  2.8609075157_wp,  2.9199626718_wp,  2.9651349636_wp,  2.9964243909_wp, & !93-96
      &  3.0138309539_wp,  3.0173546524_wp,  3.0069954867_wp,  2.9827534565_wp, & !97-100
      &  2.9446285619_wp,  2.8926208030_wp,  2.8267301797_wp] !101-103

   !> Element-specific local q scaling of the chemical hardness for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kqeta(max_elem) = [&
      &  1.8222099473_wp, -0.2575679643_wp,  0.4393826724_wp,  1.1102162003_wp, & !1-4
      &  1.2310070946_wp,  0.9818102022_wp,  0.1502230497_wp,  0.4134119032_wp, & !5-8
      &  2.5030512016_wp, -0.4998596384_wp,  2.1023399046_wp,  1.1266337899_wp, & !9-12
      &  1.3785272689_wp,  0.9471745876_wp,  1.6601128471_wp, -0.0156796346_wp, & !13-16
      &  0.6525286877_wp, -2.8148211211_wp,  1.8730352397_wp,  0.4148795713_wp, & !17-20
      &  1.9811917137_wp,  1.3666346630_wp,  0.4773540249_wp,  0.6660383739_wp, & !21-24
      &  0.4949831426_wp,  0.9260098769_wp,  1.4071496248_wp,  0.7430722161_wp, & !25-28
      &  1.4792830405_wp,  1.4211880229_wp,  0.6613271421_wp,  1.3109487181_wp, & !29-32
      &  0.9539967321_wp,  0.0441858334_wp,  0.8506553360_wp, -0.7778128954_wp, & !33-36
      &  2.4456255294_wp,  0.6279760783_wp,  0.8504097502_wp,  0.1275277215_wp, & !37-40
      &  1.0244946467_wp,  0.3991961865_wp,  0.3007399180_wp,  0.8892405348_wp, & !41-44
      &  1.0358999274_wp,  0.5910349581_wp,  1.3306044793_wp,  1.0116510919_wp, & !45-48
      &  1.2017335753_wp,  1.0749481071_wp,  1.5278450966_wp,  0.3830852785_wp, & !49-52
      &  0.8039617911_wp, -1.6689377641_wp,  1.3153512507_wp,  0.6850807472_wp, & !53-56
      &  0.4068053082_wp,  0.2805275842_wp,  0.2612355874_wp,  0.2457254002_wp, & !57-60
      &  0.2339970224_wp,  0.2260504541_wp,  0.2218856952_wp,  0.2215027459_wp, & !61-64
      &  0.2249016061_wp,  0.2320822757_wp,  0.2430447548_wp,  0.2577890434_wp, & !65-68
      &  0.2763151415_wp,  0.2986230491_wp,  0.3247127662_wp,  0.9329386915_wp, & !69-72
      &  1.1124975975_wp,  0.3105056463_wp,  0.2119489274_wp,  0.3490965682_wp, & !73-76
      &  0.9303004996_wp,  0.6578893166_wp,  0.7625190003_wp,  0.6067448860_wp, & !77-80
      &  1.1098111282_wp,  0.9571986961_wp,  1.4674965889_wp,  0.7713149335_wp, & !81-84
      &  0.5513455799_wp, -0.7227615433_wp,  1.2895674764_wp,  0.5960416182_wp, & !85-88
      &  0.1671277145_wp,  0.1575313114_wp,  0.2863965715_wp,  0.4002506248_wp, & !89-92
      &  0.4990934714_wp,  0.5829251111_wp,  0.6517455441_wp,  0.7055547703_wp, & !93-96
      &  0.7443527897_wp,  0.7681396024_wp,  0.7769152082_wp,  0.7706796073_wp, & !97-100
      &  0.7494327996_wp,  0.7131747852_wp,  0.6619055639_wp] !101-103

   !> Element-specific bond capacitance for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cap(max_elem) = [&
      &  3.4358731613_wp,  0.2563012350_wp,  1.7336935111_wp,  1.4252599447_wp, & !1-4
      &  1.9821377790_wp,  7.9575330990_wp,  5.2650283958_wp,  5.3394223720_wp, & !5-8
      &  4.7702507597_wp,  0.5095753028_wp,  5.7961811482_wp,  2.8738819069_wp, & !9-12
      &  1.5730116016_wp,  0.7813507196_wp,  1.0337776163_wp,  1.4123845734_wp, & !13-16
      &  3.0340296817_wp,  0.5326667425_wp,  6.4794438076_wp,  4.1572236543_wp, & !17-20
      &  2.6197028418_wp,  1.9926557922_wp,  1.4258893003_wp,  3.4184301443_wp, & !21-24
      &  3.1337436912_wp,  4.5345735628_wp,  6.3426635435_wp,  4.8622181062_wp, & !25-28
      &  3.9658581319_wp,  2.4205042838_wp,  2.0153453160_wp,  1.3655709456_wp, & !29-32
      &  1.0879161652_wp,  0.8125045161_wp,  3.4331186365_wp,  1.1410555369_wp, & !33-36
      &  5.3302096260_wp,  8.9866820455_wp,  8.0879982654_wp,  1.3505819625_wp, & !37-40
      &  1.9761405818_wp,  4.8306789723_wp,  2.6167089975_wp,  4.9413659163_wp, & !41-44
      &  5.5889636514_wp,  3.7289038580_wp,  2.2978010245_wp,  2.9915912946_wp, & !45-48
      &  3.2084006372_wp,  2.4592286766_wp,  1.0482227697_wp,  1.4124670516_wp, & !49-52
      &  2.0699368746_wp,  2.3426022325_wp,  4.9766316345_wp,  4.7445931148_wp, & !53-56
      &  7.6556126582_wp,  2.2792827162_wp,  2.2265798615_wp,  2.2270872929_wp, & !57-60
      &  2.2808050104_wp,  2.3877330140_wp,  2.5478713036_wp,  2.7612198793_wp, & !61-64
      &  3.0277787411_wp,  3.3475478889_wp,  3.7205273228_wp,  4.1467170428_wp, & !65-68
      &  4.6261170489_wp,  5.1587273410_wp,  5.7445479192_wp,  1.9450532464_wp, & !69-72
      &  1.2082681633_wp,  5.4761913827_wp,  2.8688258387_wp,  3.4269533511_wp, & !73-76
      &  1.2827929585_wp,  4.2446334525_wp,  8.5466705292_wp,  2.7030553995_wp, & !77-80
      &  1.7482905639_wp,  4.5652515937_wp,  2.0750200204_wp,  2.1042278455_wp, & !81-84
      &  2.9249818593_wp,  1.1606670882_wp,  5.1339954989_wp,  5.4015367551_wp, & !85-88
      &  1.5278253705_wp,  0.7201439348_wp,  0.8778110607_wp,  1.0152634518_wp, & !89-92
      &  1.1325011080_wp,  1.2295240293_wp,  1.3063322158_wp,  1.3629256675_wp, & !93-96
      &  1.3993043843_wp,  1.4154683662_wp,  1.4114176133_wp,  1.3871521256_wp, & !97-100
      &  1.3426719030_wp,  1.2779769455_wp,  1.1930672532_wp] !101-103

   !> Element-specific covalent radii for the CN for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cov_radii(max_elem) = 0.5_wp * [&
      &  1.1980006149_wp,  2.2610217725_wp,  2.3787175190_wp,  2.4632164676_wp, & !1-4
      &  2.4613895807_wp,  2.6763007964_wp,  2.7655085211_wp,  2.6466398902_wp, & !5-8
      &  2.0647114131_wp,  2.2964278893_wp,  3.0473595746_wp,  3.3597126173_wp, & !9-12
      &  2.9413551863_wp,  3.3593150082_wp,  3.7124038217_wp,  3.7950496861_wp, & !13-16
      &  3.6218412465_wp,  2.3368507550_wp,  3.2678005729_wp,  2.6934639460_wp, & !17-20
      &  3.0942813806_wp,  3.1994190934_wp,  3.1865351525_wp,  3.0245247746_wp, & !21-24
      &  2.9516405455_wp,  2.7967091405_wp,  2.8624752847_wp,  2.9325871383_wp, & !25-28
      &  2.8750457420_wp,  3.2880254556_wp,  3.4129389757_wp,  3.5547538315_wp, & !29-32
      &  3.8824044195_wp,  4.1349986852_wp,  3.9489265588_wp,  3.6264077864_wp, & !33-36
      &  3.8936921777_wp,  3.5213571939_wp,  3.2417303558_wp,  3.5464510864_wp, & !37-40
      &  3.6460575764_wp,  3.4102131328_wp,  3.4469009914_wp,  3.3043609242_wp, & !41-44
      &  3.3574938698_wp,  3.4236287446_wp,  3.6341385880_wp,  3.8313216784_wp, & !45-48
      &  3.8271624151_wp,  4.1254263093_wp,  3.9895425056_wp,  4.5341418141_wp, & !49-52
      &  4.6758001321_wp,  4.3188262237_wp,  2.8314213049_wp,  4.8914010315_wp, & !53-56
      &  3.6004533315_wp,  2.8758092017_wp,  2.9967129499_wp,  3.1042023757_wp, & !57-60
      &  3.1982774792_wp,  3.2789382603_wp,  3.3461847191_wp,  3.4000168555_wp, & !61-64
      &  3.4404346695_wp,  3.4674381612_wp,  3.4810273305_wp,  3.4812021775_wp, & !65-68
      &  3.4679627021_wp,  3.4413089043_wp,  3.4012407842_wp,  3.5004027339_wp, & !69-72
      &  3.6576246465_wp,  3.2722427492_wp,  3.4847840299_wp,  3.3869572767_wp, & !73-76
      &  3.4600493844_wp,  3.5857257632_wp,  3.5138481825_wp,  4.0752898970_wp, & !77-80
      &  4.2705802544_wp,  4.3281934906_wp,  4.0616856521_wp,  4.3269140322_wp, & !81-84
      &  4.7950102995_wp,  4.0621301306_wp,  4.7045604278_wp,  4.3693314868_wp, & !85-88
      &  3.2349557337_wp,  2.0334056417_wp,  2.5551666814_wp,  3.0015806363_wp, & !89-92
      &  3.3726475065_wp,  3.6683672921_wp,  3.8887399930_wp,  4.0337656091_wp, & !93-96
      &  4.1034441406_wp,  4.0977755873_wp,  4.0167599494_wp,  3.8603972267_wp, & !97-100
      &  3.6286874194_wp,  3.3216305273_wp,  2.9392265506_wp] !101-103

   !> Element-specific averaged coordination number over the fitset for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_avg_cn(max_elem) = [&
      &  0.3921100000_wp,  0.0810600000_wp,  0.9910100000_wp,  0.7499500000_wp, & !1-4
      &  1.1543700000_wp,  1.6691400000_wp,  1.4250300000_wp,  0.8718100000_wp, & !5-8
      &  0.6334000000_wp,  0.0876700000_wp,  0.8740600000_wp,  0.8754800000_wp, & !9-12
      &  1.2147200000_wp,  1.1335000000_wp,  1.6890600000_wp,  1.0221600000_wp, & !13-16
      &  0.5386400000_wp,  0.0827800000_wp,  1.4096300000_wp,  1.1954700000_wp, & !17-20
      &  1.5142100000_wp,  1.7892000000_wp,  2.0646100000_wp,  1.6905600000_wp, & !21-24
      &  1.6563700000_wp,  1.5128400000_wp,  1.3179000000_wp,  0.9749800000_wp, & !25-28
      &  0.5334600000_wp,  0.6585000000_wp,  0.9696500000_wp,  1.0083100000_wp, & !29-32
      &  1.0871000000_wp,  0.8222200000_wp,  0.5449300000_wp,  0.1647100000_wp, & !33-36
      &  1.2490800000_wp,  1.2198700000_wp,  1.5657400000_wp,  1.8697600000_wp, & !37-40
      &  1.8947900000_wp,  1.7085000000_wp,  1.5521300000_wp,  1.4903300000_wp, & !41-44
      &  1.3177400000_wp,  0.6991700000_wp,  0.5528200000_wp,  0.6642200000_wp, & !45-48
      &  0.9069800000_wp,  1.0976200000_wp,  1.2183000000_wp,  0.7321900000_wp, & !49-52
      &  0.5498700000_wp,  0.2467100000_wp,  1.5680600000_wp,  1.1677300000_wp, & !53-56
      &  1.6642500000_wp,  1.6032600000_wp,  1.6032600000_wp,  1.6032600000_wp, & !57-60
      &  1.6032600000_wp,  1.6032600000_wp,  1.6032600000_wp,  1.6032600000_wp, & !61-64
      &  1.6032600000_wp,  1.6032600000_wp,  1.6032600000_wp,  1.6032600000_wp, & !65-68
      &  1.6032600000_wp,  1.6032600000_wp,  1.6032600000_wp,  1.8191000000_wp, & !69-72
      &  1.8175100000_wp,  1.6802300000_wp,  1.5224100000_wp,  1.4602600000_wp, & !73-76
      &  1.1110400000_wp,  0.9102600000_wp,  0.5218000000_wp,  1.4895900000_wp, & !77-80
      &  0.8441800000_wp,  0.9426900000_wp,  1.5171900000_wp,  0.7287100000_wp, & !81-84
      &  0.5137000000_wp,  0.2678200000_wp,  1.2122500000_wp,  1.5797100000_wp, & !85-88
      &  1.7549800000_wp,  1.7549800000_wp,  1.7549800000_wp,  1.7549800000_wp, & !89-92
      &  1.7549800000_wp,  1.7549800000_wp,  1.7549800000_wp,  1.7549800000_wp, & !93-96
      &  1.7549800000_wp,  1.7549800000_wp,  1.7549800000_wp,  1.7549800000_wp, & !97-100
      &  1.7549800000_wp,  1.7549800000_wp,  1.7549800000_wp] !101-103

contains


!> Get electronegativity for species with a given symbol
elemental function get_eeqbc_chi_sym(symbol) result(chi)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> electronegativity
   real(wp) :: chi

   chi = get_eeqbc_chi(to_number(symbol))

end function get_eeqbc_chi_sym


!> Get electronegativity for species with a given atomic number
elemental function get_eeqbc_chi_num(number) result(chi)

   !> Atomic number
   integer, intent(in) :: number

   !> electronegativity
   real(wp) :: chi

   if (number > 0 .and. number <= size(eeqbc_chi, dim=1)) then
      chi = eeqbc_chi(number)
   else
      chi = -1.0_wp
   end if

end function get_eeqbc_chi_num


!> Get hardness for species with a given symbol
elemental function get_eeqbc_eta_sym(symbol) result(eta)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> hardness
   real(wp) :: eta

   eta = get_eeqbc_eta(to_number(symbol))

end function get_eeqbc_eta_sym


!> Get hardness for species with a given atomic number
elemental function get_eeqbc_eta_num(number) result(eta)

   !> Atomic number
   integer, intent(in) :: number

   !> hardness
   real(wp) :: eta

   if (number > 0 .and. number <= size(eeqbc_eta, dim=1)) then
      eta = eeqbc_eta(number)
   else
      eta = -1.0_wp
   end if

end function get_eeqbc_eta_num


!> Get charge width for species with a given symbol
elemental function get_eeqbc_rad_sym(symbol) result(rad)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> charge width
   real(wp) :: rad

   rad = get_eeqbc_rad(to_number(symbol))

end function get_eeqbc_rad_sym


!> Get charge width for species with a given atomic number
elemental function get_eeqbc_rad_num(number) result(rad)

   !> Atomic number
   integer, intent(in) :: number

   !> charge width
   real(wp) :: rad

   if (number > 0 .and. number <= size(eeqbc_rad, dim=1)) then
      rad = eeqbc_rad(number)
   else
      rad = -1.0_wp
   end if

end function get_eeqbc_rad_num


!> Get CN scaling of the electronegativity for species with a given symbol
elemental function get_eeqbc_kcnchi_sym(symbol) result(kcnchi)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> CN scaling of EN
   real(wp) :: kcnchi

   kcnchi = get_eeqbc_kcnchi(to_number(symbol))

end function get_eeqbc_kcnchi_sym


!> Get CN scaling of the electronegativity for species with a given atomic number
elemental function get_eeqbc_kcnchi_num(number) result(kcnchi)

   !> Atomic number
   integer, intent(in) :: number

   !> CN scaling of EN
   real(wp) :: kcnchi

   if (number > 0 .and. number <= size(eeqbc_kcnchi, dim=1)) then
      kcnchi = eeqbc_kcnchi(number)
   else
      kcnchi = -1.0_wp
   end if

end function get_eeqbc_kcnchi_num


!> Get local q scaling of the electronegativity for species with a given symbol
elemental function get_eeqbc_kqchi_sym(symbol) result(kqchi)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> local q scaling of EN
   real(wp) :: kqchi

   kqchi = get_eeqbc_kqchi(to_number(symbol))

end function get_eeqbc_kqchi_sym


!> Get local q scaling of the electronegativity for species with a given atomic number
elemental function get_eeqbc_kqchi_num(number) result(kqchi)

   !> Atomic number
   integer, intent(in) :: number

   !> local q scaling of EN
   real(wp) :: kqchi

   if (number > 0 .and. number <= size(eeqbc_kqchi, dim=1)) then
      kqchi = eeqbc_kqchi(number)
   else
      kqchi = -1.0_wp
   end if

end function get_eeqbc_kqchi_num


!> Get local q scaling of the chemical hardness for species with a given symbol
elemental function get_eeqbc_kqeta_sym(symbol) result(kqeta)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> local q scaling of hardness
   real(wp) :: kqeta

   kqeta = get_eeqbc_kqeta(to_number(symbol))

end function get_eeqbc_kqeta_sym


!> Get local q scaling of the chemical hardness for species with a given atomic number
elemental function get_eeqbc_kqeta_num(number) result(kqeta)

   !> Atomic number
   integer, intent(in) :: number

   !> local q scaling of hardness
   real(wp) :: kqeta

   if (number > 0 .and. number <= size(eeqbc_kqeta, dim=1)) then
      kqeta = eeqbc_kqeta(number)
   else
      kqeta = -1.0_wp
   end if

end function get_eeqbc_kqeta_num


!> Get bond capacitance for species with a given symbol
elemental function get_eeqbc_cap_sym(symbol) result(cap)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> bond capacitance
   real(wp) :: cap

   cap = get_eeqbc_cap(to_number(symbol))

end function get_eeqbc_cap_sym


!> Get bond capacitance for species with a given atomic number
elemental function get_eeqbc_cap_num(number) result(cap)

   !> Atomic number
   integer, intent(in) :: number

   !> bond capacitance
   real(wp) :: cap

   if (number > 0 .and. number <= size(eeqbc_cap, dim=1)) then
      cap = eeqbc_cap(number)
   else
      cap = -1.0_wp
   end if

end function get_eeqbc_cap_num


!> Get covalent radius for species with a given symbol
elemental function get_eeqbc_cov_radii_sym(symbol) result(rcov)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> covalent radius
   real(wp) :: rcov

   rcov = get_eeqbc_cov_radii(to_number(symbol))

end function get_eeqbc_cov_radii_sym


!> Get covalent radius for species with a given atomic number
elemental function get_eeqbc_cov_radii_num(number) result(rcov)

   !> Atomic number
   integer, intent(in) :: number

   !> covalent radius
   real(wp) :: rcov

   if (number > 0 .and. number <= size(eeqbc_cov_radii, dim=1)) then
      rcov = eeqbc_cov_radii(number)
   else
      rcov = -1.0_wp
   end if

end function get_eeqbc_cov_radii_num


!> Get average CN for species with a given symbol
elemental function get_eeqbc_avg_cn_sym(symbol) result(avg_cn)

   !> Element symbol
   character(len=*), intent(in) :: symbol

   !> average CN
   real(wp) :: avg_cn

   avg_cn = get_eeqbc_avg_cn(to_number(symbol))

end function get_eeqbc_avg_cn_sym


!> Get average CN for species with a given atomic number
elemental function get_eeqbc_avg_cn_num(number) result(avg_cn)

   !> Atomic number
   integer, intent(in) :: number

   !> average CN
   real(wp) :: avg_cn

   if (number > 0 .and. number <= size(eeqbc_avg_cn, dim=1)) then
      avg_cn = eeqbc_avg_cn(number)
   else
      avg_cn = -1.0_wp
   end if

end function get_eeqbc_avg_cn_num


end module multicharge_param_eeqbc2024
