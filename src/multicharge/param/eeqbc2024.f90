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
      &  1.7500484721_wp,  1.0526984072_wp,  0.8799564774_wp,  1.1812492597_wp, & !1-4
      &  1.3212041667_wp,  1.7223828106_wp,  1.9243977269_wp,  2.0202199701_wp, & !5-8
      &  2.0665243348_wp,  0.4860725936_wp,  0.7996111235_wp,  0.9126524304_wp, & !9-12
      &  1.0598730902_wp,  1.3440091790_wp,  1.7452027676_wp,  1.9294302550_wp, & !13-16
      &  1.8883307190_wp,  0.9554966813_wp,  0.6235746313_wp,  0.9041626136_wp, & !17-20
      &  0.8833912883_wp,  0.9257481038_wp,  0.8686601408_wp,  0.8496202114_wp, & !21-24
      &  1.0551242931_wp,  1.1419198920_wp,  1.2477890267_wp,  1.2460708935_wp, & !25-28
      &  1.0960264160_wp,  1.0069253786_wp,  1.0302067991_wp,  1.2722661942_wp, & !29-32
      &  1.4516242676_wp,  1.7559534629_wp,  1.6029112325_wp,  0.8990990167_wp, & !33-36
      &  0.5368145278_wp,  0.8795114010_wp,  0.9766146097_wp,  0.8140575699_wp, & !37-40
      &  0.9200033125_wp,  0.9707908196_wp,  1.1121449676_wp,  1.0428359023_wp, & !41-44
      &  1.1647544472_wp,  1.1451621474_wp,  1.1431710957_wp,  1.0891346977_wp, & !45-48
      &  1.0230602557_wp,  1.2233014173_wp,  1.2220883831_wp,  1.6190775222_wp, & !49-52
      &  1.5987600838_wp,  0.8738325308_wp,  0.5201497778_wp,  0.8293959230_wp, & !53-56
      &  0.9045405345_wp,  0.7951792263_wp,  0.6930199522_wp,  0.6080724877_wp, & !57-60
      &  0.5403368329_wp,  0.4898129877_wp,  0.4565009522_wp,  0.4404007263_wp, & !61-64
      &  0.4415123101_wp,  0.4598357035_wp,  0.4953709066_wp,  0.5481179194_wp, & !65-68
      &  0.6180767418_wp,  0.7052473739_wp,  0.8096298156_wp,  0.8671796504_wp, & !69-72
      &  1.0584810052_wp,  1.1309552852_wp,  1.1345823102_wp,  1.2633956218_wp, & !73-76
      &  1.2767291962_wp,  1.3656714720_wp,  1.3344481533_wp,  1.0626915556_wp, & !77-80
      &  0.9125064985_wp,  1.0911338155_wp,  1.1920489483_wp,  1.3712452197_wp, & !81-84
      &  1.6133613549_wp,  0.9910744242_wp,  0.5096061529_wp,  0.8141916249_wp, & !85-88
      &  0.7444164650_wp,  0.9209739572_wp,  0.8712337946_wp,  0.8256069491_wp, & !89-92
      &  0.7840934209_wp,  0.7466932098_wp,  0.7134063158_wp,  0.6842327391_wp, & !93-96
      &  0.6591724795_wp,  0.6382255371_wp,  0.6213919119_wp,  0.6086716039_wp, & !97-100
      &  0.6000646130_wp,  0.5955709393_wp,  0.5951905828_wp] !101-103

   !> Element-specific chemical hardnesses for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_eta(max_elem) = [&
      &  0.4363893835_wp, 16.7215837203_wp, -0.0999410763_wp, -2.2380597504_wp, & !1-4
      & -2.6208170656_wp, -4.0267359584_wp, -2.2136385442_wp, -0.8302131282_wp, & !5-8
      & -3.2456070782_wp, 10.5191696081_wp,  0.1020059595_wp, -0.2576857625_wp, & !9-12
      & -0.1758932213_wp, -3.5735748164_wp, -4.9024526610_wp, -6.6231228053_wp, & !13-16
      & -1.2121677988_wp,  1.6023586513_wp,  0.4631331968_wp, -0.3176589560_wp, & !17-20
      & -0.0159627739_wp, -0.0826061164_wp, -0.0207700380_wp, -0.1210185666_wp, & !21-24
      & -0.0439102587_wp, -0.6811336736_wp, -0.3563000054_wp, -0.8436599135_wp, & !25-28
      & -0.1555211102_wp,  0.4423551731_wp, -0.5562475241_wp, -1.5299030336_wp, & !29-32
      & -2.6139753737_wp, -3.1879240007_wp, -0.4504079874_wp,  2.3751873858_wp, & !33-36
      &  0.4327734063_wp, -0.4919557851_wp, -0.4253814024_wp, -0.0826964082_wp, & !37-40
      & -0.4209715839_wp, -0.1366061663_wp, -0.0189035347_wp, -0.2321889200_wp, & !41-44
      & -0.1346968274_wp, -0.3191245451_wp, -0.7247952231_wp, -0.4128602594_wp, & !45-48
      & -0.3939995633_wp, -1.1587324481_wp, -0.9230230000_wp, -2.3355388900_wp, & !49-52
      & -0.4013922241_wp,  2.0242970416_wp,  0.1692007876_wp, -0.6500977825_wp, & !53-56
      & -0.1654800966_wp, -0.4573562397_wp, -0.6166625265_wp, -0.7438748869_wp, & !57-60
      & -0.8389933209_wp, -0.9020178284_wp, -0.9329484095_wp, -0.9317850642_wp, & !61-64
      & -0.8985277925_wp, -0.8331765943_wp, -0.7357314697_wp, -0.6061924187_wp, & !65-68
      & -0.4445594413_wp, -0.2508325375_wp, -0.0250117072_wp, -0.1581081093_wp, & !69-72
      & -0.7140164551_wp, -0.1957106066_wp, -0.1923000488_wp, -0.9160081969_wp, & !73-76
      & -0.8066420086_wp, -1.2598169305_wp, -1.8470963207_wp, -0.2372416688_wp, & !77-80
      & -0.1762113830_wp, -1.1073326781_wp, -0.0620055121_wp, -0.6414516463_wp, & !81-84
      & -0.2983776500_wp,  0.5073478774_wp,  0.4011656410_wp, -0.2566414046_wp, & !85-88
      & -0.2131530393_wp, -0.0137707685_wp, -0.0223759361_wp, -0.0323937506_wp, & !89-92
      & -0.0438242121_wp, -0.0566673204_wp, -0.0709230756_wp, -0.0865914777_wp, & !93-96
      & -0.1036725267_wp, -0.1221662226_wp, -0.1420725654_wp, -0.1633915552_wp, & !97-100
      & -0.1861231918_wp, -0.2102674753_wp, -0.2358244057_wp] !101-103

   !> Element-specific charge widths for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_rad(max_elem) = [&
      &  0.4948060787_wp,  3.8937905064_wp,  0.4207644693_wp,  0.2453021010_wp, & !1-4
      &  0.2229794675_wp,  0.1622556638_wp,  0.2388906321_wp,  0.3635933312_wp, & !5-8
      &  0.1522765621_wp,  1.3325071781_wp,  0.6337957149_wp,  0.7607282716_wp, & !9-12
      &  0.6486350301_wp,  0.2019805393_wp,  0.1521629677_wp,  0.1122745793_wp, & !13-16
      &  0.2621613704_wp,  0.4730017988_wp,  1.0034948675_wp,  0.6509535338_wp, & !17-20
      &  1.1361102335_wp,  1.0891087619_wp,  0.7940609626_wp,  0.8413545024_wp, & !21-24
      &  0.8096705687_wp,  0.4858701129_wp,  0.7992925448_wp,  0.4925566182_wp, & !25-28
      &  0.6066392150_wp,  1.5426746901_wp,  0.5148521159_wp,  0.4054401028_wp, & !29-32
      &  0.2856252228_wp,  0.2202603178_wp,  0.4427577812_wp,  3.0127830265_wp, & !33-36
      &  0.9025441363_wp,  0.6693672655_wp,  0.7625330059_wp,  0.8850560267_wp, & !37-40
      &  0.7790271504_wp,  0.9006173190_wp,  1.0215604229_wp,  0.6079642387_wp, & !41-44
      &  0.7962528734_wp,  0.6199263440_wp,  0.5349204670_wp,  0.6835104999_wp, & !45-48
      &  0.6365539496_wp,  0.4572444362_wp,  0.5398926656_wp,  0.2807259254_wp, & !49-52
      &  0.4695397417_wp,  2.4767895683_wp,  0.7809602772_wp,  0.6453086375_wp, & !53-56
      &  1.0363730007_wp,  0.6382084916_wp,  0.5444348120_wp,  0.4679023806_wp, & !57-60
      &  0.4086111973_wp,  0.3665612621_wp,  0.3417525750_wp,  0.3341851360_wp, & !61-64
      &  0.3438589451_wp,  0.3707740024_wp,  0.4149303077_wp,  0.4763278612_wp, & !65-68
      &  0.5549666628_wp,  0.6508467125_wp,  0.7639680103_wp,  1.1418438777_wp, & !69-72
      &  0.7373710500_wp,  0.8983413360_wp,  0.7933160956_wp,  0.4563374348_wp, & !73-76
      &  0.5992139646_wp,  0.4089588601_wp,  0.3213680881_wp,  0.6979534360_wp, & !77-80
      &  0.6808608432_wp,  0.4659807860_wp,  1.2050301680_wp,  0.8083552234_wp, & !81-84
      &  0.5557856918_wp,  0.8468384603_wp,  0.9072825504_wp,  0.8514470736_wp, & !85-88
      &  0.7751156726_wp,  1.3506770411_wp,  1.2357943143_wp,  1.1286142576_wp, & !89-92
      &  1.0291368711_wp,  0.9373621548_wp,  0.8532901087_wp,  0.7769207328_wp, & !93-96
      &  0.7082540270_wp,  0.6472899915_wp,  0.5940286261_wp,  0.5484699309_wp, & !97-100
      &  0.5106139058_wp,  0.4804605510_wp,  0.4580098663_wp] !101-103

   !> Element-specific CN scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kcnchi(max_elem) = [&
      &  1.2966834027_wp,  3.3917716211_wp,  0.2794844816_wp, -0.3296209845_wp, & !1-4
      & -0.1650917868_wp,  0.2871111972_wp,  0.4989642511_wp,  0.6918815394_wp, & !5-8
      &  1.5560916379_wp,  6.2196389442_wp, -0.2218580168_wp, -0.6677254542_wp, & !9-12
      & -0.3262811638_wp,  0.1096673450_wp,  0.5447945730_wp,  0.7984338829_wp, & !13-16
      &  1.4657291466_wp,  6.6547510887_wp,  0.1333155235_wp, -0.1019355232_wp, & !17-20
      & -0.8688857024_wp, -0.8177242012_wp, -0.7323344811_wp, -0.5855108129_wp, & !21-24
      & -0.7997652442_wp, -0.7167183641_wp, -0.4911358966_wp, -0.5511568681_wp, & !25-28
      &  0.0715026869_wp, -0.2836096896_wp, -0.2523661988_wp, -0.0008315805_wp, & !29-32
      &  0.2746325234_wp,  0.6971731027_wp,  1.0506766539_wp,  1.6825599754_wp, & !33-36
      &  0.0488026283_wp, -1.0742604257_wp, -0.8324237981_wp, -0.9400656318_wp, & !37-40
      & -0.4446941396_wp, -0.3082071229_wp, -0.4749748963_wp, -0.2250871447_wp, & !41-44
      & -0.1458031008_wp,  0.3304723756_wp, -0.0283067031_wp, -0.1068815225_wp, & !45-48
      & -0.3366529721_wp, -0.0077009414_wp, -0.0091824123_wp,  0.4301516555_wp, & !49-52
      &  1.0125302562_wp,  1.0407378421_wp, -0.6701933136_wp, -1.1368178059_wp, & !53-56
      & -0.9952692224_wp, -0.0808225036_wp, -0.2332924229_wp, -0.3718046190_wp, & !57-60
      & -0.4963590919_wp, -0.6069558416_wp, -0.7035948681_wp, -0.7862761714_wp, & !61-64
      & -0.8549997515_wp, -0.9097656085_wp, -0.9505737422_wp, -0.9774241528_wp, & !65-68
      & -0.9903168401_wp, -0.9892518043_wp, -0.9742290453_wp, -0.9760600486_wp, & !69-72
      & -0.5786116338_wp, -0.2214485477_wp, -0.4098471183_wp, -0.2374843373_wp, & !73-76
      & -0.2087017206_wp,  0.1463235529_wp,  0.5393082645_wp,  0.0299762559_wp, & !77-80
      & -0.5042005627_wp, -0.0087399131_wp, -0.0622529907_wp, -0.0066969174_wp, & !81-84
      &  0.8420991662_wp,  1.0860840951_wp, -1.9702741166_wp, -1.2652561265_wp, & !85-88
      &  0.0729251828_wp, -0.9405783943_wp, -1.0100873130_wp, -1.0570453678_wp, & !89-92
      & -1.0814525588_wp, -1.0833088857_wp, -1.0626143488_wp, -1.0193689480_wp, & !93-96
      & -0.9535726833_wp, -0.8652255546_wp, -0.7543275621_wp, -0.6208787056_wp, & !97-100
      & -0.4648789852_wp, -0.2863284009_wp, -0.0852269527_wp] !101-103

   !> Element-specific local q scaling of the electronegativity for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kqchi(max_elem) = [&
      &  0.7536628563_wp, -0.7928603399_wp,  3.0215141501_wp,  2.2481863628_wp, & !1-4
      &  1.9500288592_wp,  1.2853565733_wp,  1.1452707516_wp,  0.9138607237_wp, & !5-8
      &  0.2113663935_wp, -2.4560052154_wp,  3.0512780429_wp,  2.4361119597_wp, & !9-12
      &  2.5721090297_wp,  1.7057657722_wp,  1.1812315455_wp,  1.6019394481_wp, & !13-16
      &  1.4522992583_wp, -0.7878460491_wp,  3.0454094184_wp,  2.7226086331_wp, & !17-20
      &  2.9415302614_wp,  2.6632649943_wp,  2.2456126456_wp,  2.0358310909_wp, & !21-24
      &  1.9621123559_wp,  2.2889391653_wp,  1.7223126493_wp,  1.9738331072_wp, & !25-28
      &  2.3095815986_wp,  2.2565291827_wp,  2.2202438464_wp,  1.9227836388_wp, & !29-32
      &  1.4581719966_wp,  1.5439678882_wp,  1.9347027449_wp,  0.4701522608_wp, & !33-36
      &  3.5242173098_wp,  2.6363519218_wp,  2.6528722128_wp,  2.3945899061_wp, & !37-40
      &  2.3704975137_wp,  2.1110654023_wp,  1.6703242165_wp,  2.3192944186_wp, & !41-44
      &  2.2028748120_wp,  1.8548407102_wp,  1.9891660199_wp,  2.0492772247_wp, & !45-48
      &  2.7569259143_wp,  2.0988338945_wp,  2.1235183203_wp,  1.9558700254_wp, & !49-52
      &  2.2462818543_wp,  0.6239202114_wp,  3.2320067110_wp,  2.6615463162_wp, & !53-56
      &  2.4539903153_wp,  2.5535528745_wp,  2.4268879878_wp,  2.3194234111_wp, & !57-60
      &  2.2311591442_wp,  2.1620951873_wp,  2.1122315402_wp,  2.0815682031_wp, & !61-64
      &  2.0701051759_wp,  2.0778424585_wp,  2.1047800511_wp,  2.1509179535_wp, & !65-68
      &  2.2162561659_wp,  2.3007946881_wp,  2.4045335203_wp,  2.6806374754_wp, & !69-72
      &  2.4003713187_wp,  2.1327804161_wp,  1.8387035182_wp,  2.1402184648_wp, & !73-76
      &  1.6617101413_wp,  1.7672750906_wp,  1.9182018763_wp,  1.6743244331_wp, & !77-80
      &  2.7824646655_wp,  2.3974384763_wp,  2.0455202575_wp,  1.8808281365_wp, & !81-84
      &  2.3230291681_wp,  0.6322323373_wp,  3.4597909888_wp,  2.6583951533_wp, & !85-88
      &  2.8642150313_wp,  2.5882734739_wp,  2.6556922117_wp,  2.7152780818_wp, & !89-92
      &  2.7670310843_wp,  2.8109512191_wp,  2.8470384862_wp,  2.8752928857_wp, & !93-96
      &  2.8957144175_wp,  2.9083030817_wp,  2.9130588782_wp,  2.9099818071_wp, & !97-100
      &  2.8990718683_wp,  2.8803290618_wp,  2.8537533877_wp] !101-103

   !> Element-specific local q scaling of the chemical hardness for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_kqeta(max_elem) = [&
      &  2.1395929425_wp, -7.3617035539_wp,  0.6505518174_wp,  1.2370404121_wp, & !1-4
      &  1.3432971069_wp,  0.8816071592_wp,  0.2839088075_wp,  0.4693824536_wp, & !5-8
      &  2.6542988922_wp,  0.7609713722_wp,  1.5049131566_wp,  1.2481897948_wp, & !9-12
      &  1.2902774472_wp,  1.0726190048_wp,  1.8670630099_wp,  0.1359291998_wp, & !13-16
      &  0.1682685474_wp, -1.7321763019_wp,  1.7371190449_wp,  0.3244821964_wp, & !17-20
      &  1.9192639221_wp,  1.3343230052_wp,  0.6207563526_wp,  0.5178907968_wp, & !21-24
      &  0.4993882153_wp,  1.0689181095_wp,  1.1898819688_wp,  0.8250816694_wp, & !25-28
      &  1.3075645243_wp,  1.4572323840_wp,  0.6297367373_wp,  1.2555969543_wp, & !29-32
      &  1.0261216879_wp,  0.2238526662_wp,  1.0756350116_wp, -0.3481929742_wp, & !33-36
      &  2.2587549248_wp,  0.5761682168_wp,  0.8014079276_wp,  0.2259558595_wp, & !37-40
      &  1.0606549876_wp,  0.3764925563_wp,  0.1463005340_wp,  1.0814100757_wp, & !41-44
      &  0.9719588327_wp,  0.4414631456_wp,  1.2014651475_wp,  1.0414933308_wp, & !45-48
      &  1.0976715875_wp,  1.1141930862_wp,  1.5340540236_wp,  0.4659966479_wp, & !49-52
      &  0.8051128565_wp, -1.4352091023_wp,  1.3312677685_wp,  0.6124487966_wp, & !53-56
      &  0.3387688982_wp,  0.3137694073_wp,  0.2522053683_wp,  0.1967591511_wp, & !57-60
      &  0.1474307559_wp,  0.1042201825_wp,  0.0671274310_wp,  0.0361525014_wp, & !61-64
      &  0.0112953937_wp, -0.0074438921_wp, -0.0200653561_wp, -0.0265689981_wp, & !65-68
      & -0.0269548183_wp, -0.0212228166_wp, -0.0093729930_wp,  1.0045621095_wp, & !69-72
      &  1.1971115293_wp,  0.3412676640_wp,  0.1412297445_wp,  0.5402463837_wp, & !73-76
      &  0.9930279645_wp,  0.6799486848_wp,  0.7132512217_wp,  0.7135772074_wp, & !77-80
      &  1.0541292189_wp,  1.0177470373_wp,  1.5022288204_wp,  0.7733513665_wp, & !81-84
      &  0.6371862709_wp, -0.9286384593_wp,  1.0512977406_wp,  0.4478113543_wp, & !85-88
      &  0.2706215216_wp,  0.1211096926_wp,  0.2542354374_wp,  0.3674348447_wp, & !89-92
      &  0.4607079145_wp,  0.5340546468_wp,  0.5874750416_wp,  0.6209690989_wp, & !93-96
      &  0.6345368187_wp,  0.6281782009_wp,  0.6018932457_wp,  0.5556819530_wp, & !97-100
      &  0.4895443227_wp,  0.4034803550_wp,  0.2974900497_wp] !101-103

   !> Element-specific bond capacitance for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cap(max_elem) = [&
      &  3.2154265670_wp,  0.8290139573_wp,  2.1691150655_wp,  1.0944312636_wp, & !1-4
      &  1.4125792215_wp,  6.3369044258_wp,  4.5466291057_wp,  4.6537714106_wp, & !5-8
      &  4.5664793255_wp,  0.4857350415_wp,  3.9323133631_wp,  2.5876523882_wp, & !9-12
      &  1.1387203262_wp,  0.5834834124_wp,  0.5059679595_wp,  0.9373279723_wp, & !13-16
      &  1.7895689139_wp,  0.4444507597_wp,  5.9910132240_wp,  3.3204059133_wp, & !17-20
      &  2.4174028088_wp,  1.8244283865_wp,  1.2483865008_wp,  3.1680724947_wp, & !21-24
      &  2.5144743282_wp,  4.1583681619_wp,  3.6687900788_wp,  4.3025759444_wp, & !25-28
      &  3.7741322618_wp,  2.3040521845_wp,  1.7855042050_wp,  1.1771845352_wp, & !29-32
      &  0.8504814672_wp,  0.7947488555_wp,  4.3775869892_wp,  1.1417251132_wp, & !33-36
      &  5.1670808115_wp,  9.1394718584_wp,  6.8038459938_wp,  1.4491855197_wp, & !37-40
      &  1.7446574228_wp,  4.5152330241_wp,  2.2668311355_wp,  3.5320934917_wp, & !41-44
      &  5.0734525381_wp,  3.8797996216_wp,  2.2718530358_wp,  2.8178716034_wp, & !45-48
      &  3.1480373819_wp,  1.9656247358_wp,  1.0667891632_wp,  1.2453572293_wp, & !49-52
      &  1.9598464332_wp,  2.3223908517_wp,  3.9972831293_wp,  4.5707025621_wp, & !53-56
      &  6.3046306574_wp,  2.3576829597_wp,  2.2231545627_wp,  2.1539603437_wp, & !57-60
      &  2.1501003027_wp,  2.2115744398_wp,  2.3383827549_wp,  2.5305252481_wp, & !61-64
      &  2.7880019193_wp,  3.1108127686_wp,  3.4989577959_wp,  3.9524370013_wp, & !65-68
      &  4.4712503847_wp,  5.0553979461_wp,  5.7048796856_wp,  1.9109639168_wp, & !69-72
      &  1.0576347738_wp,  6.0070230852_wp,  2.2995013066_wp,  3.5105623066_wp, & !73-76
      &  1.1946216529_wp,  4.0247915184_wp,  8.0442460842_wp,  2.3756419385_wp, & !77-80
      &  1.5730903675_wp,  3.9340302707_wp,  1.9607007347_wp,  1.9751036203_wp, & !81-84
      &  4.5256278524_wp,  1.0804157381_wp,  4.7684322814_wp,  5.6592308767_wp, & !85-88
      &  1.5224734026_wp,  0.9287382385_wp,  1.1572187001_wp,  1.3496427135_wp, & !89-92
      &  1.5060102787_wp,  1.6263213957_wp,  1.7105760645_wp,  1.7587742851_wp, & !93-96
      &  1.7709160575_wp,  1.7470013818_wp,  1.6870302578_wp,  1.5910026857_wp, & !97-100
      &  1.4589186654_wp,  1.2907781969_wp,  1.0865812802_wp] !101-103

   !> Element-specific covalent radii for the CN for the EEQ_BC charges.
   real(wp), parameter :: eeqbc_cov_radii(max_elem) = 0.5_wp * [&
      &  1.2128733184_wp,  1.6286580002_wp,  2.1287453721_wp,  2.4216168539_wp, & !1-4
      &  2.3827439506_wp,  2.6808053318_wp,  2.7801790216_wp,  2.6691898814_wp, & !5-8
      &  2.0340698263_wp,  2.2228122744_wp,  3.1231211724_wp,  3.4622345110_wp, & !9-12
      &  3.0206576420_wp,  3.3230544359_wp,  3.5302182643_wp,  3.6655369036_wp, & !13-16
      &  3.6758848016_wp,  2.0677269305_wp,  2.8112907187_wp,  2.5461685765_wp, & !17-20
      &  3.2062127869_wp,  3.2436299544_wp,  3.2854214942_wp,  3.0190088169_wp, & !21-24
      &  2.9535273817_wp,  2.7895598404_wp,  2.8376327000_wp,  2.9404675083_wp, & !25-28
      &  2.8253153933_wp,  3.2908203354_wp,  3.3965111905_wp,  3.5308901294_wp, & !29-32
      &  3.9007783954_wp,  4.1124008156_wp,  3.8731396330_wp,  3.4442923110_wp, & !33-36
      &  3.5415477250_wp,  3.4887755843_wp,  3.3173607964_wp,  3.5992976332_wp, & !37-40
      &  3.7163289263_wp,  3.4067152378_wp,  3.4404124873_wp,  3.3251375471_wp, & !41-44
      &  3.3101656501_wp,  3.3806949974_wp,  3.5807347169_wp,  3.8147599730_wp, & !45-48
      &  3.7834157311_wp,  4.1335898647_wp,  3.9963439273_wp,  4.4909570543_wp, & !49-52
      &  4.6304944015_wp,  4.1412331913_wp,  1.7187274006_wp,  4.8533357994_wp, & !53-56
      &  3.7230258294_wp,  0.8597845508_wp,  1.3407180597_wp,  1.7740601478_wp, & !57-60
      &  2.1598108151_wp,  2.4979700616_wp,  2.7885378873_wp,  3.0315142922_wp, & !61-64
      &  3.2268992763_wp,  3.3746928397_wp,  3.4748949822_wp,  3.5275057040_wp, & !65-68
      &  3.5325250050_wp,  3.4899528852_wp,  3.3997893446_wp,  3.5821167198_wp, & !69-72
      &  3.7289425673_wp,  3.2800254061_wp,  3.5084712866_wp,  3.3566713050_wp, & !73-76
      &  3.4329474808_wp,  3.5189878745_wp,  3.4616764774_wp,  4.0585318982_wp, & !77-80
      &  4.2520114640_wp,  4.3051538750_wp,  4.0622599580_wp,  4.3336545861_wp, & !81-84
      &  4.6902082889_wp,  3.9747831153_wp,  4.5466756003_wp,  4.2951554725_wp, & !85-88
      &  2.6121533852_wp,  2.4938107984_wp,  2.8635632616_wp,  3.1787071560_wp, & !89-92
      &  3.4392424816_wp,  3.6451692385_wp,  3.7964874266_wp,  3.8931970460_wp, & !93-96
      &  3.9352980966_wp,  3.9227905784_wp,  3.8556744915_wp,  3.7339498358_wp, & !97-100
      &  3.5576166114_wp,  3.3266748182_wp,  3.0411244562_wp] !101-103

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
