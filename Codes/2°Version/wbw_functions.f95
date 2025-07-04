module wbw_functions

   implicit none

contains

   !====================================================================
   !Function to get the number of gray gases
   !====================================================================
   integer function wbw_number_gray_gases(model_id)
      
      implicit none
      character(*),intent(in) :: model_id
   
      selectcase(trim(model_id))
      case('p1-j1')
         wbw_number_gray_gases = 1
      case('p1-j2')
         wbw_number_gray_gases = 2
      case('p1-j3')
         wbw_number_gray_gases = 3
      case('p1-j4')
         wbw_number_gray_gases = 4

      case('p10-j1')
         wbw_number_gray_gases = 1
      case('p10-j2')
         wbw_number_gray_gases = 2
      case('p10-j3')
         wbw_number_gray_gases = 3
      case('p10-j4')
         wbw_number_gray_gases = 4

      case('p20-j1')
         wbw_number_gray_gases = 1
      case('p20-j2')
         wbw_number_gray_gases = 2
      case('p20-j3')
         wbw_number_gray_gases = 3
      case('p20-j4')
         wbw_number_gray_gases = 4
      endselect
   endfunction wbw_number_gray_gases

   !====================================================================
   !Function to get the total number of bands
   !====================================================================
   integer function wbw_number_bands(model_id)
      
      implicit none
      character(*),intent(in) :: model_id
      
      selectcase(trim(model_id))
      case('p1-j4')
         wbw_number_bands = 5
      case('p10-j4')
         wbw_number_bands = 5
      case('p20-j4')
         wbw_number_bands = 5
      endselect
      wbw_number_bands = 5

   endfunction wbw_number_bands
   
   !====================================================================
   !Subroutine to output the lower and upper bounds of the bands
   !====================================================================
   subroutine wbw_bounds(wbw_lbound,wbw_ubound,model_id,band_id)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use precision_parameters, only: dp
      implicit none
      character(*),intent(in) :: model_id
      integer,intent(in) :: band_id
      real(dp),intent(out) :: wbw_lbound,wbw_ubound
      real(dp) :: lbnd(10),ubnd(10)

      !-----------------------------------------------------------------
      !Define band bounds (in 1/cm)
      !-----------------------------------------------------------------
      selectcase(trim(model_id))
         case default
            lbnd(1:5) = (/ 0._dp   , 1000._dp, 2600._dp, 4400._dp,&
                           6000._dp  /)
            ubnd(1:5) = (/ 1000._dp, 2600._dp, 4400._dp, 6000._dp,&
                           10000._dp /)
      endselect

      !-----------------------------------------------------------------
      !Output
      !-----------------------------------------------------------------
      wbw_lbound = lbnd(band_id)*100._dp
      wbw_ubound = ubnd(band_id)*100._dp

   endsubroutine wbw_bounds

   !====================================================================
   !Function to compute the pressure-based absorption coefficient
   !====================================================================
   real(dp) function wbw_kappaj(model_id,band_id,gas_id)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use precision_parameters, only: dp
      implicit none
      character(*),intent(in) :: model_id
      integer,intent(in) :: band_id,gas_id
      real(dp) :: kcoeff(10)
      
      !-----------------------------------------------------------------
      !Special case: transparent window
      !-----------------------------------------------------------------
      if (gas_id.eq.5) then
         wbw_kappaj = 0._dp
         return
      endif     

      !-----------------------------------------------------------------
      !Define coefficients
      !-----------------------------------------------------------------
      selectcase(trim(model_id))
      case('p1-j1')
         selectcase(band_id)
         case(1)
            kcoeff(1:1) = (/ 4.24031_dp /)
         case(2)
            kcoeff(1:1) = (/ 3.00774_dp /)
         case(3)
            kcoeff(1:1) = (/ 0.96483_dp /)
         case(4)
            kcoeff(1:1) = (/ 0.28362_dp /)
         case(5)
            kcoeff(1:1) = (/ 0.19912_dp /)
         endselect
      case('p1-j2')
         selectcase(band_id)
         case(1)
            kcoeff(1:2) = (/ 0.74677_dp, 12.68970_dp /)
         case(2)
            kcoeff(1:2) = (/ 0.69680_dp, 16.24040_dp /)
         case(3)
            kcoeff(1:2) = (/ 0.30558_dp, 3.69346_dp /)
         case(4)
            kcoeff(1:2) = (/ 0.17626_dp, 1.52787_dp /)
         case(5)
            kcoeff(1:2) = (/ 0.12440_dp, 1.35349_dp /)
         endselect
      case('p1-j3')
         selectcase(band_id)
         case(1)
            kcoeff(1:3) = (/ 0.44617_dp, 4.36598_dp, 30.97350_dp /)
         case(2)
            kcoeff(1:3) = (/ 0.37450_dp, 3.12647_dp, 42.81430_dp /)
         case(3)
            kcoeff(1:3) = (/ 0.21039_dp, 1.68448_dp, 10.68390_dp /)
         case(4)
            kcoeff(1:3) = (/ 0.06952_dp, 0.43553_dp, 2.76981_dp /)
         case(5)
            kcoeff(1:3) = (/ 0.05345_dp, 0.05854_dp, 0.65831_dp /)
         endselect
      case('p1-j4')
         selectcase(band_id)
         case(1)
            kcoeff(1:4) = (/ 0.33344_dp, 2.09808_dp, 11.47687_dp, 76.98296_dp /)
         case(2)
            kcoeff(1:4) = (/ 0.27366_dp, 1.74412_dp, 11.12446_dp, 81.01230_dp /)
         case(3)
            kcoeff(1:4) = (/ 0.16483_dp, 1.07802_dp, 4.29001_dp, 25.40130_dp /)
         case(4)
            kcoeff(1:4) = (/ 0.11323_dp, 0.53958_dp, 2.03853_dp, 9.43548_dp /)
         case(5)
            kcoeff(1:4) = (/ 0.04016_dp, 0.27374_dp, 1.24082_dp, 6.13351_dp /)
         endselect

      case('p10-j1')
         selectcase(band_id)
         case(1)
            kcoeff(1:1) = (/ 4.12663_dp /)
         case(2)
            kcoeff(1:1) = (/ 1.65414_dp /)
         case(3)
            kcoeff(1:1) = (/ 0.64413_dp /)
         case(4)
            kcoeff(1:1) = (/ 0.17313_dp /)
         case(5)
            kcoeff(1:1) = (/ 0.06638_dp /)
         endselect
      case('p10-j2')
         selectcase(band_id)
         case(1)
            kcoeff(1:2) = (/ 0.63922_dp, 15.73136_dp /)
         case(2)
            kcoeff(1:2) = (/ 0.40978_dp, 9.20599_dp /)
         case(3)
            kcoeff(1:2) = (/ 0.08451_dp, 2.19982_dp /)
         case(4)
            kcoeff(1:2) = (/ 0.03631_dp, 0.39088_dp /)
         case(5)
            kcoeff(1:2) = (/ 0.02283_dp, 0.26064_dp /)
         endselect
      case('p10-j3')
         selectcase(band_id)
         case(1)
            kcoeff(1:3) = (/ 0.20142_dp, 2.02924_dp, 22.46205_dp /)
         case(2)
            kcoeff(1:3) = (/ 0.11471_dp, 1.79025_dp, 28.38313_dp /)
         case(3)
            kcoeff(1:3) = (/ 0.04986_dp, 0.41824_dp, 4.08123_dp /)
         case(4)
            kcoeff(1:3) = (/ 0.02557_dp, 0.19806_dp, 0.78262_dp /)
         case(5)
            kcoeff(1:3) = (/ 0.01564_dp, 0.13589_dp, 0.63747_dp /)
         endselect
      case('p10-j4')
         selectcase(band_id)
         case(1)
            kcoeff(1:4) = (/ 0.11155_dp, 0.78009_dp, 5.02382_dp, 28.70622_dp /)
         case(2)
            kcoeff(1:4) = (/ 0.07470_dp, 0.68312_dp, 4.09838_dp, 44.70016_dp /)
         case(3)
            kcoeff(1:4) = (/ 0.03792_dp, 0.20466_dp, 1.46723_dp, 7.34271_dp /)
         case(4)
            kcoeff(1:4) = (/ 0.01992_dp, 0.11064_dp, 0.38134_dp, 1.42667_dp /)
         case(5)
            kcoeff(1:4) = (/ 0.01258_dp, 0.09153_dp, 0.29208_dp, 1.25525_dp /)
         endselect
      
      case('p20-j1')
         selectcase(band_id)
         case(1)
            kcoeff(1:1) = (/ 3.69110_dp /)
         case(2)
            kcoeff(1:1) = (/ 1.94609_dp /)
         case(3)
            kcoeff(1:1) = (/ 0.33158_dp /)
         case(4)
            kcoeff(1:1) = (/ 0.15241_dp /)
         case(5)
            kcoeff(1:1) = (/ 0.04783_dp /)
         endselect
      case('p20-j2')
         selectcase(band_id)
         case(1)
            kcoeff(1:2) = (/ 0.43172_dp, 14.03305_dp /)
         case(2)
            kcoeff(1:2) = (/ 0.34838_dp, 10.81278_dp /)
         case(3)
            kcoeff(1:2) = (/ 0.07354_dp, 1.90469_dp /)
         case(4)
            kcoeff(1:2) = (/ 0.02293_dp, 0.33427_dp /)
         case(5)
            kcoeff(1:2) = (/ 0.01321_dp, 0.20655_dp /)
         endselect
      case('p20-j3')
         selectcase(band_id)
         case(1)
            kcoeff(1:3) = (/ 0.15501_dp, 1.46635_dp, 18.62135_dp /)
         case(2)
            kcoeff(1:3) = (/ 0.09687_dp, 1.32986_dp, 18.75608_dp /)
         case(3)
            kcoeff(1:3) = (/ 0.03588_dp, 0.26773_dp, 3.47165_dp /)
         case(4)
            kcoeff(1:3) = (/ 0.01501_dp, 0.11792_dp, 0.55932_dp /)
         case(5)
            kcoeff(1:3) = (/ 0.00967_dp, 0.08974_dp, 0.45268_dp /)
         endselect
      case('p20-j4')
         selectcase(band_id)
         case(1)
            kcoeff(1:4) = (/ 0.10414_dp, 0.85732_dp, 6.21553_dp, 36.19979_dp /)
         case(2)
            kcoeff(1:4) = (/ 0.07285_dp, 0.72279_dp, 4.47017_dp, 58.58844_dp /)
         case(3)
            kcoeff(1:4) = (/ 0.02663_dp, 0.14594_dp, 1.04175_dp, 6.35332_dp /)
         case(4)
            kcoeff(1:4) = (/ 0.01107_dp, 0.05950_dp, 0.29593_dp, 1.05692_dp /)
         case(5)
            kcoeff(1:4) = (/ 0.00779_dp, 0.05077_dp, 0.20746_dp, 0.86880_dp /)
         endselect
      endselect

      !-----------------------------------------------------------------
      !Obtain correct coefficient
      !-----------------------------------------------------------------
      wbw_kappaj = kcoeff(gas_id)
      
   endfunction wbw_kappaj

   !====================================================================
   !Function to compute the weighted coefficient
   !====================================================================
   real(dp) recursive function wbw_aj(model_id,band_id,gas_id,tmp) &
      result (aj)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use precision_parameters, only: dp
      implicit none
      character(*),intent(in) :: model_id
      integer,intent(in) :: band_id,gas_id
      integer :: j,k,nj,nk
      real(dp),intent(in) :: tmp
      real(dp) :: asum,bcoeff(10,10)
      
      !-----------------------------------------------------------------
      !Define coefficients
      !-----------------------------------------------------------------
      selectcase(trim(model_id))
      case('p1-j1')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 0.54966_dp, 0.31147_dp, 0.01976_dp, -0.06033_dp, 0.01262_dp, -0.00009_dp /)
         case(2)
            bcoeff(1,1:6) = (/ 0.00687_dp, 2.43177_dp, -3.06419_dp, 1.94727_dp, -0.60663_dp, 0.07358_dp /)
         case(3)
            bcoeff(1,1:6) = (/ -0.14087_dp, 1.26817_dp, -0.75569_dp, 0.22382_dp, -0.02534_dp, -0.00032_dp /)
         case(4)
            bcoeff(1,1:6) = (/ -0.25505_dp, 1.61567_dp, -1.63134_dp, 0.96136_dp, -0.28746_dp, 0.03371_dp /)
         case(5)
            bcoeff(1,1:6) = (/ -0.24766_dp, 1.11581_dp, -0.68516_dp, 0.00010_dp, 0.02870_dp, -0.00802_dp /)
         endselect
      case('p1-j2')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 0.19459_dp, 0.40261_dp, -0.11373_dp, -0.06001_dp, 0.02958_dp, -0.00343_dp /)
            bcoeff(2,1:6) = (/ 0.41241_dp, -0.26837_dp, 0.57448_dp, -0.40908_dp, 0.14357_dp, -0.01952_dp /)
         case(2)
            bcoeff(1,1:6) = (/ -0.16503_dp, 2.25381_dp, -2.83552_dp, 1.69209_dp, -0.48710_dp, 0.05504_dp /)
            bcoeff(2,1:6) = (/ 0.13939_dp, 0.46359_dp, -0.58348_dp, 0.45465_dp, -0.17243_dp, 0.02401_dp /)
         case(3)
            bcoeff(1,1:6) = (/ 0.27702_dp, -0.40105_dp, 1.01056_dp, -0.76321_dp, 0.25827_dp, -0.03262_dp /)
            bcoeff(2,1:6) = (/ -0.33554_dp, 1.50782_dp, -1.49656_dp, 0.80819_dp, -0.22734_dp, 0.02555_dp /)
         case(4)
            bcoeff(1,1:6) = (/ -0.18057_dp, 1.12125_dp, -0.95736_dp, 0.56563_dp, -0.17277_dp, 0.02063_dp /)
            bcoeff(2,1:6) = (/ -0.08597_dp, 0.55462_dp, -0.67012_dp, 0.39592_dp, -0.11660_dp, 0.01349_dp /)
         case(5)
            bcoeff(1,1:6) = (/ -0.09763_dp, 0.37126_dp, 0.38124_dp, -0.52417_dp, 0.20103_dp, -0.02634_dp /)
            bcoeff(2,1:6) = (/ -0.11348_dp, 0.54538_dp, -0.63699_dp, 0.32834_dp, -0.08096_dp, 0.00783_dp /)
         endselect
      case('p1-j3')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 0.27646_dp, -0.18734_dp, 0.80059_dp, -0.75205_dp, 0.25946_dp, -0.03139_dp /)
            bcoeff(2,1:6) = (/ -0.04876_dp, 1.09536_dp, -1.58484_dp, 1.17407_dp, -0.38095_dp, 0.04518_dp /)
            bcoeff(3,1:6) = (/ 0.41041_dp, -0.83660_dp, 1.39158_dp, -1.01681_dp, 0.33770_dp, -0.04200_dp /)
         case(2)
            bcoeff(1,1:6) = (/ -0.22879_dp, 2.21020_dp, -2.87387_dp, 1.66137_dp, -0.45030_dp, 0.04733_dp /)
            bcoeff(2,1:6) = (/ 0.10213_dp, 0.35760_dp, -0.35192_dp, 0.32667_dp, -0.14682_dp, 0.02297_dp /)
            bcoeff(3,1:6) = (/ 0.06635_dp, 0.42851_dp, -0.57170_dp, 0.37920_dp, -0.12178_dp, 0.01490_dp /)
         case(3)
            bcoeff(1,1:6) = (/ 0.27981_dp, -0.22975_dp, 0.70565_dp, -0.59515_dp, 0.21354_dp, -0.02786_dp /)
            bcoeff(2,1:6) = (/ -0.12986_dp, 0.42873_dp, -0.11541_dp, 0.01700_dp, -0.00046_dp, -0.00027_dp /)
            bcoeff(3,1:6) = (/ -0.16746_dp, 0.87736_dp, -0.98500_dp, 0.53201_dp, -0.14655_dp, 0.01617_dp /)
         case(4)
            bcoeff(1,1:6) = (/ -0.21407_dp, 1.38445_dp, -1.32870_dp, 0.77646_dp, -0.22679_dp, 0.02580_dp /)
            bcoeff(2,1:6) = (/ -0.11089_dp, 0.60364_dp, -0.61567_dp, 0.40403_dp, -0.12947_dp, 0.01540_dp /)
            bcoeff(3,1:6) = (/ -0.04684_dp, 0.34112_dp, -0.41200_dp, 0.22348_dp, -0.06068_dp, 0.00669_dp /)
         case(5)
            bcoeff(1,1:6) = (/ -1.05240_dp, 4.04570_dp, -3.36410_dp, 0.18963_dp, 0.60761_dp, -0.15121_dp /)
            bcoeff(2,1:6) = (/ 0.98290_dp, -3.89391_dp, 4.42341_dp, -1.28345_dp, -0.20043_dp, 0.09693_dp /)
            bcoeff(3,1:6) = (/ -0.17337_dp, 0.81430_dp, -0.86331_dp, 0.40318_dp, -0.08861_dp, 0.00739_dp /)
         endselect
      case('p1-j4')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 0.21575_dp, -0.16163_dp, 0.82673_dp, -0.76506_dp, 0.23618_dp, -0.02260_dp /)
            bcoeff(2,1:6) = (/ 0.06176_dp, 0.63157_dp, -1.17768_dp, 0.95994_dp, -0.30869_dp, 0.03310_dp /)
            bcoeff(3,1:6) = (/ 0.06790_dp, 0.23707_dp, 0.25197_dp, -0.38495_dp, 0.16130_dp, -0.02127_dp /)
            bcoeff(4,1:6) = (/ 0.29684_dp, -0.61917_dp, 0.71026_dp, -0.39851_dp, 0.11144_dp, -0.01285_dp /)
         case(2)
            bcoeff(1,1:6) = (/ -0.24092_dp, 2.23257_dp, -3.04608_dp, 1.80849_dp, -0.49786_dp, 0.05274_dp /)
            bcoeff(2,1:6) = (/ 0.06531_dp, 0.28328_dp, -0.05037_dp, 0.02206_dp, -0.02886_dp, 0.00704_dp /)
            bcoeff(3,1:6) = (/ 0.11022_dp, 0.31240_dp, -0.58593_dp, 0.46625_dp, -0.15992_dp, 0.01997_dp /)
            bcoeff(4,1:6) = (/ 0.05206_dp, 0.08593_dp, 0.01283_dp, -0.04399_dp, 0.01497_dp, -0.00156_dp /)
         case(3)
            bcoeff(1,1:6) = (/ 0.33469_dp, -0.41099_dp, 0.92450_dp, -0.73289_dp, 0.25442_dp,-0.03265_dp /)
            bcoeff(2,1:6) = (/ -0.13132_dp, 0.51104_dp, -0.35869_dp, 0.17965_dp, -0.04521_dp, 0.00465_dp /)
            bcoeff(3,1:6) = (/ -0.05698_dp, 0.27732_dp, -0.04617_dp, -0.03835_dp, 0.01564_dp, -0.00215_dp /)
            bcoeff(4,1:6) = (/ -0.12122_dp, 0.61735_dp, -0.77876_dp, 0.44666_dp, -0.12561_dp, 0.01408_dp /)
         case(4)
            bcoeff(1,1:6) = (/ -0.15529_dp, 1.00977_dp, -0.75017_dp, 0.38204_dp, -0.10700_dp, 0.01244_dp /)
            bcoeff(2,1:6) = (/ -0.08884_dp, 0.52746_dp, -0.66771_dp, 0.46662_dp, -0.15144_dp, 0.01826_dp /)
            bcoeff(3,1:6) = (/ -0.02385_dp, 0.14315_dp, -0.06660_dp, -0.01443_dp, 0.01228_dp, -0.00178_dp /)
            bcoeff(4,1:6) = (/ -0.01855_dp, 0.14137_dp, -0.22424_dp, 0.14877_dp, -0.04540_dp, 0.00530_dp /)
         case(5)
            bcoeff(1,1:6) = (/ 0.02224_dp, -0.00549_dp, 0.77540_dp, -0.62015_dp, 0.18667_dp, -0.02054_dp /)
            bcoeff(2,1:6) = (/ -0.10455_dp, 0.38567_dp, -0.08679_dp, -0.13265_dp, 0.07673_dp, -0.01193_dp /)
            bcoeff(3,1:6) = (/ -0.06684_dp, 0.32308_dp, -0.36758_dp, 0.18214_dp, -0.04285_dp, 0.00394_dp /)
            bcoeff(4,1:6) = (/ -0.02028_dp, 0.11407_dp, -0.17024_dp, 0.11116_dp, -0.03387_dp, 0.00395_dp /)
         endselect

      case('p10-j1')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 0.78542_dp, 0.34500_dp, -0.26514_dp, 0.11042_dp, -0.02216_dp, 0.00154_dp /)
         case(2)
            bcoeff(1,1:6) = (/ 0.40411_dp, 2.07458_dp, -2.87316_dp, 1.88622_dp, -0.59263_dp, 0.07180_dp /)
         case(3)
            bcoeff(1,1:6) = (/ 0.64413_dp, 0.16064_dp, 0.29526_dp, -0.32518_dp, 0.12058_dp, -0.01572_dp /)
         case(4)
            bcoeff(1,1:6) = (/ -0.26155_dp, 2.02261_dp, -1.72032_dp, 0.88130_dp, -0.24693_dp, 0.02853_dp /)
         case(5)
            bcoeff(1,1:6) = (/ -0.26645_dp, 1.48281_dp, -0.84276_dp, 0.17644_dp, 0.00848_dp, -0.00582_dp /)
         endselect
      case('p10-j2')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 0.05901_dp, 2.16814_dp, -3.34018_dp, 2.20694_dp, -0.69669_dp, 0.08536_dp /)
            bcoeff(2,1:6) = (/ 0.72017_dp, -1.57769_dp, 2.68756_dp, -1.83760_dp, 0.59227_dp, -0.07370_dp /)
         case(2)
            bcoeff(1,1:6) = (/ 0.65412_dp, -0.33768_dp, -0.00441_dp, 0.18120_dp, -0.09021_dp, 0.01397_dp /)
            bcoeff(2,1:6) = (/ 0.15978_dp, 0.67701_dp, 0.01049_dp, -0.50697_dp, 0.29024_dp, -0.04916_dp /)
         case(3)
            bcoeff(1,1:6) = (/ 0.80453_dp, -0.39958_dp, 0.00259_dp, 0.15349_dp, -0.07461_dp, 0.01122_dp /)
            bcoeff(2,1:6) = (/ -0.07505_dp, 0.64271_dp, -0.00270_dp, -0.22032_dp, 0.09986_dp, -0.01401_dp /)
         case(4)
            bcoeff(1,1:6) = (/ 0.05602_dp, 0.35256_dp, -0.00108_dp, -0.05686_dp, 0.00146_dp, 0.00352_dp /)
            bcoeff(2,1:6) = (/ -0.09346_dp, 0.72679_dp, 0.00502_dp, -0.41149_dp, 0.23239_dp, -0.03910_dp /)
         case(5)
            bcoeff(1,1:6) = (/ 0.06727_dp, 0.24622_dp, -0.00037_dp, 0.09784_dp, -0.08077_dp, 0.01577_dp /)
            bcoeff(2,1:6) = (/ -0.19168_dp, 0.73850_dp, 0.00083_dp, -0.49988_dp, 0.27528_dp, -0.04467_dp /)
         endselect
      case('p10-j3')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 0.46704_dp, -0.09185_dp, 0.01319_dp, -0.29651_dp, 0.20524_dp, -0.03818_dp /)
            bcoeff(2,1:6) = (/ -0.16962_dp, 0.51974_dp, -0.01357_dp, 0.01316_dp, -0.05893_dp, 0.01571_dp /)
            bcoeff(3,1:6) = (/ 0.59745_dp, -0.24034_dp, 0.00301_dp, 0.12721_dp, -0.05424_dp, 0.00670_dp /)
         case(2)
            bcoeff(1,1:6) = (/ 0.62279_dp, -0.67805_dp, -0.00771_dp, 0.48171_dp, -0.27159_dp, 0.04521_dp /)
            bcoeff(2,1:6) = (/ -0.15248_dp, 1.06367_dp, 0.01099_dp, -0.83223_dp, 0.48691_dp, -0.08262_dp /)
            bcoeff(3,1:6) = (/ 0.39982_dp, -0.12727_dp, 0.00189_dp, 0.08779_dp, -0.05184_dp, 0.00839_dp /)
         case(3)
            bcoeff(1,1:6) = (/ 0.45006_dp, 0.12947_dp, 0.01026_dp, -0.37502_dp, 0.23784_dp, -0.04235_dp /)
            bcoeff(2,1:6) = (/ 0.58928_dp, -0.81190_dp, -0.01349_dp, 0.86906_dp, -0.51699_dp, 0.08895_dp /)
            bcoeff(3,1:6) = (/ -0.29824_dp, 0.97403_dp, 0.00446_dp, -0.62829_dp, 0.34557_dp, -0.05658_dp /)
         case(4)
            bcoeff(1,1:6) = (/ 0.02923_dp, 0.54951_dp, -0.00141_dp, -0.31885_dp, 0.14642_dp, -0.01951_dp /)
            bcoeff(2,1:6) = (/ 0.00690_dp, 0.01939_dp, 0.00255_dp, 0.28182_dp, -0.15729_dp, 0.02443_dp /)
            bcoeff(3,1:6) = (/ -0.07555_dp, 0.60064_dp, 0.00259_dp, -0.50851_dp, 0.28440_dp, -0.04643_dp /)
         case(5)
            bcoeff(1,1:6) = (/ 0.07447_dp, 0.34322_dp, 0.00073_dp, -0.05547_dp, 0.00934_dp, 0.00081_dp /)
            bcoeff(2,1:6) = (/ -0.09669_dp, 0.29009_dp, -0.00143_dp, -0.00586_dp, -0.01196_dp, 0.00313_dp /)
            bcoeff(3,1:6) = (/ -0.09131_dp, 0.42554_dp, 0.00156_dp, -0.39227_dp, 0.22432_dp, -0.03704_dp /)
         endselect
      case('p10-j4')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/  0.08967_dp,  1.77489_dp, -3.90697_dp,  3.09948_dp, -1.08285_dp,  0.14087_dp /)
            bcoeff(2,1:6) = (/ -0.13981_dp, -0.01087_dp,  1.88687_dp, -2.12059_dp,  0.84764_dp, -0.11788_dp /)
            bcoeff(3,1:6) = (/  0.40656_dp, -0.82597_dp,  0.59384_dp,  0.07842_dp, -0.12698_dp,  0.02368_dp /)
            bcoeff(4,1:6) = (/  0.43170_dp, -0.14133_dp,  0.31748_dp, -0.32728_dp,  0.13201_dp, -0.01866_dp /)
         case(2)
            bcoeff(1,1:6) = (/  1.32032_dp, -4.07054_dp,  5.26909_dp, -3.28977_dp,  0.99784_dp, -0.11808_dp /)
            bcoeff(2,1:6) = (/ -1.06731_dp,  5.02513_dp, -6.67656_dp,  4.10672_dp, -1.20254_dp,  0.13684_dp /)
            bcoeff(3,1:6) = (/  0.32520_dp, -0.36950_dp,  0.69065_dp, -0.36776_dp,  0.06940_dp, -0.00293_dp /)
            bcoeff(4,1:6) = (/  0.10603_dp,  0.67413_dp, -1.10872_dp,  0.77599_dp, -0.25307_dp,  0.03135_dp /)
         case(3)
            bcoeff(1,1:6) = (/ -0.33537_dp,  3.07083_dp, -4.63043_dp,  3.07008_dp, -0.96408_dp,  0.11666_dp /)
            bcoeff(2,1:6) = (/  1.28877_dp, -3.60875_dp,  4.66999_dp, -2.86253_dp,  0.84860_dp, -0.09773_dp /)
            bcoeff(3,1:6) = (/  0.03292_dp, -0.32208_dp,  0.96967_dp, -0.72035_dp,  0.24622_dp, -0.03228_dp /)
            bcoeff(4,1:6) = (/ -0.34596_dp,  1.63897_dp, -1.76582_dp,  0.92040_dp, -0.24725_dp,  0.02695_dp /)
         case(4)
            bcoeff(1,1:6) = (/  0.18756_dp, -0.37186_dp,  1.71069_dp, -1.71029_dp,  0.65115_dp, -0.08725_dp /)
            bcoeff(2,1:6) = (/ -0.18896_dp,  1.38680_dp, -2.68364_dp,  2.32497_dp, -0.86381_dp,  0.11649_dp /)
            bcoeff(3,1:6) = (/  0.06175_dp, -0.64476_dp,  2.07855_dp, -1.78040_dp,  0.64339_dp, -0.08532_dp /)
            bcoeff(4,1:6) = (/ -0.23486_dp,  1.46247_dp, -2.12938_dp,  1.38233_dp, -0.42979_dp,  0.05183_dp /)
         case(5)
            bcoeff(1,1:6) = (/  0.18179_dp, -0.11823_dp,  0.75692_dp, -0.59841_dp,  0.18320_dp, -0.01987_dp /)
            bcoeff(2,1:6) = (/ -0.22801_dp,  1.03398_dp, -1.40703_dp,  1.01436_dp, -0.33414_dp,  0.04077_dp /)
            bcoeff(3,1:6) = (/  0.11290_dp, -0.66571_dp,  1.70661_dp, -1.42685_dp,  0.49793_dp, -0.06331_dp /)
            bcoeff(4,1:6) = (/ -0.21611_dp,  1.05955_dp, -1.45363_dp,  0.88362_dp, -0.25351_dp,  0.02811_dp /)
         endselect
      
      case('p20-j1')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 1.00226_dp, -0.54794_dp, 1.11402_dp, -0.87136_dp, 0.30495_dp, -0.03985_dp /)
         case(2)
            bcoeff(1,1:6) = (/ 0.29521_dp, 2.45625_dp, -3.37635_dp, 2.22156_dp, -0.70383_dp, 0.08622_dp /)
         case(3)
            bcoeff(1,1:6) = (/ 0.68393_dp, 0.23250_dp, 0.21845_dp, -0.32905_dp, 0.14019_dp, -0.02012_dp /)
         case(4)
            bcoeff(1,1:6) = (/ -0.18149_dp, 1.63319_dp, -0.97959_dp, 0.31879_dp, -0.05815_dp, 0.00500_dp /)
         case(5)
            bcoeff(1,1:6) = (/ -0.17348_dp, 1.15692_dp, -0.33619_dp, -0.14212_dp, 0.10151_dp, -0.01623_dp /)
         endselect
      case('p20-j2')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ 0.25452_dp, 1.60261_dp, -2.75830_dp, 1.86911_dp, -0.59024_dp, 0.07165_dp /)
            bcoeff(2,1:6) = (/ 0.75561_dp, -1.95510_dp, 3.53356_dp, -2.50613_dp, 0.81986_dp, -0.10223_dp /)
         case(2)
            bcoeff(1,1:6) = (/ 1.01188_dp, -1.87687_dp, 2.31965_dp, -1.43520_dp, 0.44013_dp, -0.05266_dp /)
            bcoeff(2,1:6) = (/ -0.60755_dp, 4.06616_dp, -5.36373_dp, 3.45090_dp, -1.08090_dp, 0.13132_dp /)
         case(3)
            bcoeff(1,1:6) = (/ 0.02030_dp, 2.39417_dp, -3.60051_dp, 2.35584_dp, -0.72910_dp, 0.08725_dp /)
            bcoeff(2,1:6) = (/ 0.61771_dp, -1.73668_dp, 3.15663_dp, -2.24189_dp, 0.72977_dp, -0.09043_dp /)
         case(4)
            bcoeff(1,1:6) = (/ 0.30426_dp, -0.79527_dp, 2.02423_dp, -1.65864_dp, 0.57113_dp, -0.07165_dp /)
            bcoeff(2,1:6) = (/ -0.40112_dp, 2.16098_dp, -2.43911_dp, 1.52903_dp, -0.47582_dp, 0.05738_dp /)
         case(5)
            bcoeff(1,1:6) = (/ 0.24812_dp, -0.41240_dp, 0.99030_dp, -0.60557_dp, 0.15432_dp, -0.01406_dp /)
            bcoeff(2,1:6) = (/ -0.32569_dp, 1.34728_dp, -0.96658_dp, 0.26463_dp, -0.00565_dp, -0.00611_dp /)
         endselect
      case('p20-j3')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ -0.39277_dp, 4.39230_dp, -7.83540_dp, 5.74737_dp, -1.92620_dp, 0.24377_dp /)
            bcoeff(2,1:6) = (/ 1.26780_dp, -6.30089_dp, 11.27670_dp, -8.43472_dp, 2.86472_dp, -0.36559_dp /)
            bcoeff(3,1:6) = (/ 0.08544_dp, 1.95703_dp, -3.39102_dp, 2.59159_dp, -0.89254_dp, 0.11468_dp /)
         case(2)
            bcoeff(1,1:6) = (/ 1.32785_dp, -3.64399_dp, 4.34470_dp, -2.53851_dp, 0.73120_dp, -0.08301_dp /)
            bcoeff(2,1:6) = (/ -1.16178_dp, 4.98994_dp, -5.82155_dp, 3.29881_dp, -0.91752_dp, 0.10107_dp /)
            bcoeff(3,1:6) = (/ 0.37050_dp, 0.43988_dp, -1.08827_dp, 0.97936_dp, -0.37636_dp, 0.05189_dp /)
         case(3)
            bcoeff(1,1:6) = (/ -0.56125_dp, 3.35700_dp, -4.25979_dp, 2.42773_dp, -0.66836_dp, 0.07218_dp /)
            bcoeff(2,1:6) = (/ 1.15870_dp, -2.43077_dp, 2.30665_dp, -0.90888_dp, 0.14498_dp, -0.00461_dp /)
            bcoeff(3,1:6) = (/ -0.01304_dp, -0.02726_dp, 1.23128_dp, -1.26372_dp, 0.48998_dp, -0.06759_dp /)
         case(4)
            bcoeff(1,1:6) = (/ 0.38122_dp, -1.25252_dp, 3.25545_dp, -2.91897_dp, 1.07649_dp, -0.14257_dp /)
            bcoeff(2,1:6) = (/ -0.20195_dp, 1.22843_dp, -2.43728_dp, 2.28988_dp, -0.89164_dp, 0.12340_dp /)
            bcoeff(3,1:6) = (/ -0.22505_dp, 1.19471_dp, -0.72005_dp, 0.01908_dp, 0.09303_dp, -0.01982_dp /)
         case(5)
            bcoeff(1,1:6) = (/ 0.26118_dp, -0.37301_dp, 1.11297_dp, -0.84384_dp, 0.26348_dp, -0.02972_dp /)
            bcoeff(2,1:6) = (/ -0.11537_dp, 0.42189_dp, -0.39991_dp, 0.35353_dp, -0.13642_dp, 0.01811_dp /)
            bcoeff(3,1:6) = (/ -0.19390_dp, 0.85931_dp, -0.54510_dp, 0.00288_dp, 0.07621_dp, -0.01557_dp /)
         endselect
      case('p20-j4')
         selectcase(band_id)
         case(1)
            bcoeff(1,1:6) = (/ -0.48870_dp, 4.61247_dp, -8.40442_dp, 6.30211_dp, -2.14461_dp, 0.274302_dp /)
            bcoeff(2,1:6) = (/ 1.29443_dp, -6.33832_dp, 11.7547_dp, -9.13175_dp, 3.17892_dp, -0.412438_dp /)
            bcoeff(3,1:6) = (/ -0.23020_dp, 2.37182_dp, -4.54194_dp, 3.78984_dp, -1.37441_dp, 0.182673_dp /)
            bcoeff(4,1:6) = (/ 0.37519_dp, -0.511232_dp, 1.07094_dp, -0.919529_dp, 0.337725_dp, -0.0453184_dp /)
         case(2)
            bcoeff(1,1:6) = (/ 1.43838_dp, -4.18105_dp, 5.00503_dp, -2.90524_dp, 0.82701_dp, -0.09273_dp /)
            bcoeff(2,1:6) = (/ -1.13189_dp, 4.60862_dp, -5.37978_dp, 2.95510_dp, -0.78221_dp, 0.08144_dp /)
            bcoeff(3,1:6) = (/ 0.12309_dp, 0.89357_dp, -1.37004_dp, 1.08494_dp, -0.40361_dp, 0.05550_dp /)
            bcoeff(4,1:6) = (/ 0.14060_dp, 0.35365_dp, -0.68262_dp, 0.52368_dp, -0.18067_dp, 0.02315_dp /)
         case(3)
            bcoeff(1,1:6) = (/ -0.16574_dp, 1.24119_dp, -1.08242_dp, 0.23636_dp, 0.05073_dp, -0.01827_dp /)
            bcoeff(2,1:6) = (/ -0.03331_dp, 2.39732_dp, -4.46452_dp, 3.50991_dp, -1.24834_dp, 0.16561_dp /)
            bcoeff(3,1:6) = (/ 1.52295_dp, -5.67984_dp, 8.20953_dp, -5.43816_dp, 1.73077_dp, -0.21284_dp /)
            bcoeff(4,1:6) = (/ -0.71640_dp, 2.84433_dp, -3.20869_dp, 1.80405_dp, -0.51408_dp, 0.05836_dp /)
         case(4)
            bcoeff(1,1:6) = (/ 0.53652_dp, -2.16383_dp, 4.88395_dp, -4.23173_dp, 1.55501_dp, -0.20708_dp /)
            bcoeff(2,1:6) = (/ -0.39671_dp, 2.42200_dp, -4.39692_dp, 3.63718_dp, -1.33951_dp, 0.18113_dp /)
            bcoeff(3,1:6) = (/ 0.24983_dp, -1.42629_dp, 3.05135_dp, -2.34328_dp, 0.81736_dp, -0.10729_dp /)
            bcoeff(4,1:6) = (/ -0.35723_dp, 1.97383_dp, -2.69302_dp, 1.69800_dp, -0.52415_dp, 0.06337_dp /)
         case(5)
            bcoeff(1,1:6) = (/ 0.40928_dp, -0.98407_dp, 2.06185_dp, -1.52370_dp, 0.48512_dp, -0.05665_dp /)
            bcoeff(2,1:6) = (/ -0.32679_dp, 1.44024_dp, -2.10747_dp, 1.54384_dp, -0.51536_dp, 0.06371_dp /)
            bcoeff(3,1:6) = (/ 0.20921_dp, -1.00012_dp, 2.10499_dp, -1.61442_dp, 0.54134_dp, -0.06732_dp /)
            bcoeff(4,1:6) = (/ -0.26490_dp, 1.19694_dp, -1.45433_dp, 0.79042_dp, -0.20559_dp, 0.02083_dp /)
         endselect
      endselect
      
      !-----------------------------------------------------------------
      !Compute the weighting coefficient for the transparent window
      !-----------------------------------------------------------------
      if (gas_id.eq.5) then
         nj = wbw_number_gray_gases(model_id)
         asum = 0._dp
         do j=1,nj
            asum = asum + wbw_aj(model_id,band_id,j,tmp)
         enddo
         aj = 1._dp - asum

      !-----------------------------------------------------------------
      !Compute the weighting coefficient for the gray gas
      !-----------------------------------------------------------------
      else
         asum = 0._dp; nk = 6       
         do k=1,nk
            asum = asum + bcoeff(gas_id,k)*((tmp/1000._dp)**(k-1))
         enddo
         aj = asum
      endif

   endfunction wbw_aj

endmodule wbw_functions
