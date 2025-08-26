!###############################################
!Functions and subroutines related to the 
!radiative transfer solution via the WSGG model
!###############################################
module wbm_functions
   
   use wbm_parameters
   implicit none

contains
   
   !====================================================================
   !Subroutine to load the WBM data from an external file
   !====================================================================
   subroutine load_wbm_data
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
                                get_file_unit
      use precision_parameters, only: dp
      implicit none
      character(200) :: file_name
      integer :: ierr
      integer :: inb,isp
      integer :: nnb,nsp,ntg,nxs
      integer,allocatable,dimension(:) :: file_unit,nnb_spc
      real(dp),allocatable,dimension(:,:) :: lnb_spc,unb_spc
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      nsp = wbm_data_nsp
      allocate(file_unit(nsp),stat=ierr)
      call CheckMemAlloc('file_unit',ierr)
      allocate(nnb_spc(nsp),stat=ierr)
      call CheckMemAlloc('nnb_spc',ierr)
      
      !-----------------------------------------------------------------
      !Initialize variables
      !-----------------------------------------------------------------
      wbm_data_ntg = 0
      wbm_data_nxs = 0
      
      !-----------------------------------------------------------------     
      !Read initial lines of data
      !-----------------------------------------------------------------
      reading_loop_1: do isp=1,nsp
         file_name = trim(wbm_prefix)//trim(wbm_data_file(isp))         !File name
         call CheckFileExists(file_name)                                !Check if file exists
         file_unit(isp) = get_file_unit()                               !Assign unit
         open(file=file_name,unit=file_unit(isp),action='read',&        !Open unit
              form='unformatted')
         read(file_unit(isp)) nnb_spc(isp),wbm_data_ntg(isp),&          !Read array sizes
            wbm_data_nxs(isp)
         if (isp.gt.1) call assert(nnb_spc(isp).eq.nnb_spc(isp-1))      !Check if the number of narrow bands is 
                                                                        !  the same across all data files
      enddo reading_loop_1
      number_wbm_bands = nnb_spc(1)

      !-----------------------------------------------------------------
      !Allocate data file arrays
      !-----------------------------------------------------------------
      nnb = number_wbm_bands
      ntg = maxval(wbm_data_ntg)
      nxs = maxval(wbm_data_nxs)
      if (allocated(wbm_data_tg)) deallocate(wbm_data_tg)
      allocate(wbm_data_tg(ntg,nsp),stat=ierr)
      call CheckMemAlloc('wbm_data_tg',ierr)
      if (allocated(wbm_data_xs)) deallocate(wbm_data_xs)
      allocate(wbm_data_xs(nxs,nsp),stat=ierr)
      call CheckMemAlloc('wbm_data_xs',ierr)
      if (allocated(wbm_kappa)) deallocate(wbm_kappa)
      allocate(wbm_kappa(ntg,nxs,nnb,nsp),stat=ierr)
      call CheckMemAlloc('wbm_kappa',ierr)
      if (allocated(wbm_lbound)) deallocate(wbm_lbound)
      allocate(wbm_lbound(nnb),stat=ierr)
      call CheckMemAlloc('wbm_lbound',ierr)
      if (allocated(wbm_ubound)) deallocate(wbm_ubound)
      allocate(wbm_ubound(nnb),stat=ierr)
      call CheckMemAlloc('wbm_ubound',ierr)

      !-----------------------------------------------------------------
      !Second reading loop: array values
      !-----------------------------------------------------------------
      !Allocate auxiliary arrays
      allocate(lnb_spc(nnb,nsp),stat=ierr)
      call CheckMemAlloc('lnb_spc',ierr)
      allocate(unb_spc(nnb,nsp),stat=ierr)
      call CheckMemAlloc('unb_spc',ierr)
      
      !Read data
      reading_loop_2: do isp=1,nsp
         read(file_unit(isp)) lnb_spc(:,isp),unb_spc(:,isp)
         read(file_unit(isp)) wbm_data_tg(:,isp)
         read(file_unit(isp)) wbm_data_xs(:,isp)
         read(file_unit(isp)) wbm_kappa(:,:,:,isp)
         if (isp.gt.1) then
            do inb=1,nnb
               call assert(dabs(lnb_spc(inb,isp)-lnb_spc(inb,isp-1))&
                           .lt.1e-8_dp)
               call assert(dabs(unb_spc(inb,isp)-unb_spc(inb,isp-1))&
                           .lt.1e-8_dp)
            enddo
         endif
         
         close(file_unit(isp))
      enddo reading_loop_2
      wbm_lbound = lnb_spc(:,1)
      wbm_ubound = unb_spc(:,1)

      !-----------------------------------------------------------------
      !Deallocate auxiliary arrays
      !-----------------------------------------------------------------
      deallocate(file_unit,nnb_spc)
      deallocate(lnb_spc,unb_spc)

   endsubroutine load_wbm_data

   !====================================================================
   !Function to define the number of bands fo the model
   !====================================================================
   subroutine wbm_prepare_model(model_id)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------      
      use comp_functions, only: CheckMemAlloc,shutdown
      use precision_parameters, only: dp
      character(*),intent(in) :: model_id
      integer :: ierr
      
      !-----------------------------------------------------------------
      !Special case: WBM parameters are provided in an external file
      !-----------------------------------------------------------------
      if (trim(model_id).eq.'external') then
         call load_wbm_data
         return
      endif
      
      !-----------------------------------------------------------------
      !Define number of bands according to the model
      !-----------------------------------------------------------------   
      selectcase(trim(model_id))
      case('bordbar2018')
         number_wbm_bands = 31
      case('bordbar2019')
         number_wbm_bands = 10
      case('paul')
         number_wbm_bands = 5
      case('johnson-baseline')
         number_wbm_bands = 4
      case('johnson-improved')
         number_wbm_bands = 6
      case default
         call shutdown('wbm_prepare_model: WBM not available')
      endselect
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      if (allocated(wbm_lbound)) deallocate(wbm_lbound)
      allocate(wbm_lbound(number_wbm_bands),stat=ierr)
      call CheckMemAlloc('wbm_lbound',ierr)
      if (allocated(wbm_ubound)) deallocate(wbm_ubound)
      allocate(wbm_ubound(number_wbm_bands),stat=ierr)
      call CheckMemAlloc('wbm_ubound',ierr)
      
      !-----------------------------------------------------------------
      !Mount the array with all band bounds
      !-----------------------------------------------------------------
      selectcase(trim(model_id))
      case('bordbar2018')                                               !Borbar & Hyppänen, 2018 (IJHMT)
         wbm_lbound(1:31) = (/  150.e2_dp, 351.6e2_dp,  555.e2_dp, &
                                712.e2_dp, 866.7e2_dp, 1104.e2_dp, &
                               1163.e2_dp, 1497.e2_dp, 1604.e2_dp, &
                               1850.e2_dp, 1996.e2_dp, 2065.e2_dp, &
                               2341.e2_dp, 2358.e2_dp, 2362.e2_dp, &
                               2496.e2_dp, 2659.e2_dp, 3124.e2_dp, &
                               3335.e2_dp, 3697.e2_dp, 3809.e2_dp, &
                               4011.e2_dp, 4210.e2_dp, 4531.e2_dp, &
                               5076.e2_dp, 5277.e2_dp, 5698.e2_dp, &
                               6477.e2_dp, 7136.e2_dp, 7400.e2_dp, &
                               7418.e2_dp /)
         wbm_ubound(1:31) = (/ 351.6e2_dp,  555.e2_dp,  712.e2_dp, &
                               866.7e2_dp, 1104.e2_dp, 1163.e2_dp, &
                               1497.e2_dp, 1604.e2_dp, 1850.e2_dp, &
                               1996.e2_dp, 2065.e2_dp, 2341.e2_dp, &
                               2358.e2_dp, 2362.e2_dp, 2496.e2_dp, &
                               2659.e2_dp, 3124.e2_dp, 3335.e2_dp, &
                               3697.e2_dp, 3809.e2_dp, 4011.e2_dp, &
                               4210.e2_dp, 4531.e2_dp, 5076.e2_dp, &
                               5277.e2_dp, 5698.e2_dp, 6477.e2_dp, &
                               7136.e2_dp, 7400.e2_dp, 7418.e2_dp, &
                               7472.e2_dp /)

      case('bordbar2019')                                               !Borbar et al., 2019 (AE)
         wbm_lbound(1:10) = (/ 1._dp/1.465e-6_dp, 1._dp/2.105e-6_dp, &
                               1._dp/2.632e-6_dp, 1._dp/2.878e-6_dp, &
                               1._dp/4.706e-6_dp, 1._dp/5.970e-6_dp, &
                               1._dp/7.843e-6_dp, 1._dp/12.90e-6_dp, &
                               1._dp/16.67e-6_dp, 1._dp/66.67e-6_dp /)
         wbm_ubound(1:10) = (/ 1._dp/1.338e-6_dp, 1._dp/1.770e-6_dp, &
                               1._dp/2.454e-6_dp, 1._dp/2.632e-6_dp, &
                               1._dp/4.167e-6_dp, 1._dp/5.263e-6_dp, &
                               1._dp/6.154e-6_dp, 1._dp/8.889e-6_dp, &
                               1._dp/12.90e-6_dp, 1._dp/16.67e-6_dp /)  
            
      case('paul')
         wbm_lbound(1:5) = (/ 20000._dp, 207500._dp, 240000._dp, &
                              315000._dp, 405000._dp /)
         wbm_ubound(1:5) = (/ 207500._dp, 240000._dp, 315000._dp, &
                              405000._dp, 1500000._dp /)

      case('johnson-baseline')   
         wbm_lbound(1:4) = (/ 472590._dp, 287520._dp, &
                              160000._dp, 100000._dp /)
         wbm_ubound(1:4) = (/ 754720._dp, 472590._dp, &
                              287520._dp, 160000._dp /)

      case('johnson-improved')
         wbm_lbound(1:6) = (/ 20000._dp, 95000._dp, 207500._dp, &
                              240000._dp, 315000._dp, 405000._dp /)
         wbm_ubound(1:6) = (/ 95000._dp, 207500._dp, 240000._dp, &
                              315000._dp, 405000._dp, 1000000._dp /)
      case default
         call shutdown('wbm_prepare_model: WBM not available')
      endselect

   endsubroutine wbm_prepare_model

   !====================================================================
   !Function to determine the absorption 
   !coefficient of the band in the box model
   !====================================================================
   real(dp) function box_kappa_func(model_id,iband,tmp,xs,tp)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use global_parameters, only: id_co,id_co2,id_h2o
      use precision_parameters, only: dp,small
      use math_functions, only: locate
      character(*),intent(in) :: model_id
      integer,intent(in) :: iband
      integer :: isp
      integer :: itg,ixs,ntg,nxs
      integer :: lxs,ltg,uxs,utg
      real(dp),intent(in) :: tmp,tp,xs(:)
      real(dp),dimension(1:6) :: aco,aco2,ah2o,bco,bco2,bh2o,&
         cco,cco2,ch2o,dco,dco2,dh2o
      real(dp),dimension(1:31) :: b0,b1,b2,b3,b4,b5,b6
      real(dp) :: kappa,kappap(size(xs))
      real(dp) :: Rtg,Rxs,Qtg,Qxs
   
      !-----------------------------------------------------------------
      !Mount the array with all band coefficients
      !-----------------------------------------------------------------
      kappap = 0
      selectcase(trim(model_id))
      case('bordbar2018')                                               !Borbar & Hyppänen, 2018 (IJHMT)
         !H2O
         b0(1:31) = (/ -3116.809_dp, 1358.922_dp,  -290.301_dp, &
                          37.857_dp,  -12.902_dp,  -42.214_dp,  &
                        -348.078_dp,   2537.687_dp,   1277.706_dp,  &
                         -63.611_dp,     -0.570_dp,      6.947_dp,  &
                          -4.880_dp,     -2.392_dp,     -3.330_dp,  &
                           0.603_dp,     -0.172_dp,     -9.537_dp,  &
                        -216.254_dp,   2690.634_dp,     22.888_dp,  &
                          -6.392_dp,      0.690_dp,      2.620_dp,  &
                          18.736_dp,    189.574_dp,      0.00689_dp,&
                          -6.789_dp,    118.614_dp,     26.724_dp,  &
                         -46.671_dp /)
         b1(1:31) = (/  28246.58_dp, -12072.860_dp,   2592.489_dp,  &
                        -328.578_dp,    114.698_dp,    369.576_dp,  &
                        3100.598_dp, -23106.540_dp, -11654.990_dp,  &
                         546.488_dp,     20.446_dp,    -63.901_dp,  &
                          40.749_dp,     22.697_dp,     27.110_dp,  &
                          -5.398_dp,      0.865_dp,     79.452_dp,  &
                        1913.582_dp, -24472.920_dp,   -294.920_dp,  &
                          52.296_dp,     -6.439_dp,    -24.360_dp,  &
                        -187.228_dp,  -1719.278_dp,     -0.00797_dp,&
                          62.988_dp,  -1082.151_dp,   -239.883_dp,  &
                         383.061_dp /)
         b2(1:31) = (/ -1020.932_dp,    420.809_dp,    -88.562_dp,  &
                          10.855_dp,     -3.909_dp,    -12.311_dp,  &
                        -110.389_dp,    851.226_dp,    430.618_dp,  &
                         -18.522_dp,     -1.469_dp,      2.342_dp,  &
                          -1.274_dp,     -0.825_dp,     -0.821_dp,  &
                           0.188_dp,     -0.023_dp,     -2.198_dp,  &
                         -67.437_dp,    900.064_dp,     15.377_dp,  &
                          -1.587_dp,      0.241_dp,      0.913_dp,  &
                           7.777_dp,     62.998_dp,     -0.00426_dp,&
                          -2.366_dp,     40.007_dp,      8.489_dp,  &
                         -11.757_dp /)
         b3(1:31) = (/   184.569_dp,    -70.960_dp,     14.000_dp,  &
                          -1.752_dp,      0.616_dp,      1.848_dp,  &
                          19.853_dp,   -162.507_dp,    -82.659_dp,  &
                           3.131_dp,      0.455_dp,     -0.426_dp,  &
                           0.176_dp,      0.139_dp,      0.108_dp,  &
                          -0.032_dp,      0.001_dp,      0.146_dp,  &
                          11.911_dp,   -171.371_dp,     -4.247_dp,  &
                           0.208_dp,     -0.045_dp,     -0.174_dp,  &
                          -1.743_dp,    -11.938_dp,      0.00171_dp,&
                           0.459_dp,     -7.690_dp,     -1.471_dp,  &
                           1.605_dp /)
         b4(1:31) = (/   -16.568_dp,      5.511_dp,     -0.983_dp,  &
                           0.163_dp,     -0.038_dp,     -0.111_dp,  &
                          -1.833_dp,     17.049_dp,      8.780_dp,  &
                          -0.275_dp,     -0.072_dp,      0.038_dp,  &
                          -0.009_dp,     -0.010_dp,     -0.005_dp,  &
                           0.003_dp,      0.002_dp,      0.025_dp,  &
                          -1.054_dp,     17.899_dp,      0.681_dp,  &
                          -0.010_dp,      0.004_dp,      0.017_dp,  &
                           0.228_dp,      1.240_dp,      4.0e-05_dp,&
                          -0.048_dp,      0.819_dp,      0.123_dp,  &
                          -0.088_dp  /)
         b5(1:31) = (/    50.776_dp,    -12.111_dp,      3.509_dp,  &
                          -0.631_dp,      0.110_dp,      0.328_dp,  &
                           6.711_dp,    -94.771_dp,    -50.537_dp,  &
                           1.121_dp,      0.560_dp,     -0.110_dp,  &
                           0.013_dp,      0.029_dp,      0.001_dp,  &
                          -0.018_dp,     -0.016_dp,     -0.210_dp,  &
                           3.170_dp,    -98.758_dp,     -6.504_dp,  &
                           0.015_dp,     -0.019_dp,     -0.053_dp,  &
                          -1.738_dp,     -6.829_dp,     -0.00067_dp,&
                           0.254_dp,     -4.724_dp,     -0.382_dp,  &
                           0.183_dp  /)
         b6(1:31) = (/  1694.838_dp,     89.408_dp,    -44.643_dp,  &
                           8.960_dp,     -0.780_dp,     -2.785_dp,  &
                         147.346_dp,   2288.787_dp,   1374.642_dp,  &
                           9.907_dp,     -9.270_dp,      1.230_dp,  &
                           0.164_dp,     -0.180_dp,      0.275_dp,  &
                           0.438_dp,      1.493_dp,      9.085_dp,  &
                         145.609_dp,   2357.223_dp,    342.373_dp,  &
                           1.265_dp,      0.382_dp,      0.618_dp,  &
                          68.624_dp,    169.872_dp,      0.07470_dp,&
                          -1.576_dp,    127.688_dp,      4.816_dp,  &
                          -0.169_dp  /)
         kappap(id_h2o) = b0(iband)*(tmp**6._dp)*1.e-20_dp + &
                b1(iband)*(tmp**5._dp)*1.e-17_dp + &
                b2(iband)*(tmp**4._dp)*1.e-12_dp + &
                b3(iband)*(tmp**3._dp)*1.e-8_dp  + &
                b4(iband)*(tmp**2._dp)*1.e-4_dp  + &
                b5(iband)*tmp*1.e-2_dp + b6(iband)*1.e-1_dp

         !CO2
         b0(1:31) = (/     0.315043_dp,  16.654_dp   , 2970.11_dp     ,&
                          23.029_dp   ,  11.094_dp   ,    0.357_dp    ,&
                           0.396_dp   ,   0.58_dp    ,    0.976_dp    ,&
                           2.469_dp   ,  48.089_dp   ,  667.07_dp     ,&
                        1444.158_dp   ,  93.862_dp   ,   32.461_dp    ,&
                           5.091_dp   ,   1.601_dp   ,    9.878_dp    ,&
                         -38.252_dp   ,  49.276_dp   ,    0.316_dp    ,&
                           0.147_dp   ,   0.252_dp   ,   -1.339_dp    ,&
                          -0.568_dp   ,   0.034_dp   ,   -0.005_dp    ,&
                          -0.045_dp   ,   0.015_dp   ,    0.019_dp    ,&
                           0.016_dp /)
         b1(1:31) = (/    -2.822_dp   , -127.239_dp  , -27227.68_dp   ,&
                        -210.581_dp   ,  -96.337_dp  ,     -3.234_dp  ,&
                          -3.593_dp   ,   -5.256_dp  ,     -8.741_dp  ,&
                         -23.019_dp   , -545.2_dp    ,  -8248.341_dp  ,&
                      -13826.88_dp    , -853.174_dp  ,   -293.978_dp  ,&
                         -46.144_dp   ,  -14.461_dp  ,    -78.409_dp  ,&
                         281.539_dp   , -389.89_dp   ,     -2.913_dp  ,&
                          -1.333_dp   ,   -2.039_dp  ,     11.645_dp  ,&
                           5.468_dp   ,   -0.311_dp  ,      0.047_dp  ,&
                           0.396_dp   ,   -0.133_dp  ,     -0.176_dp  ,&
                          -0.146_dp /)
         b2(1:31) = (/     0.102_dp   ,    3.456_dp  ,   1013.409_dp  ,&
                           7.766_dp   ,    3.298_dp  ,      0.118_dp  ,&
                           0.132_dp   ,    0.193_dp  ,      0.316_dp  ,&
                           0.843_dp   ,   21.501_dp  ,    423.586_dp  ,&
                         526.763_dp   ,   31.37_dp   ,     10.751_dp  ,&
                           1.69_dp    ,    0.527_dp  ,      2.285_dp  ,&
                          -6.682_dp   ,   11.294_dp  ,      0.11_dp   ,&
                           0.049_dp   ,    0.063_dp  ,     -0.391_dp  ,&
                          -0.214_dp   ,    0.012_dp  ,     -0.001_dp  ,&
                          -0.014_dp   ,    0.005_dp  ,      0.007_dp  ,&
                           0.005_dp /)
         b3(1:31) = (/    -0.019_dp   ,   -0.317_dp  ,   -196.781_dp  ,&
                          -1.46_dp    ,   -0.556_dp  ,     -0.022_dp  ,&
                          -0.025_dp   ,   -0.037_dp  ,     -0.059_dp  ,&
                          -0.155_dp   ,   -3.72_dp   ,   -117.595_dp  ,&
                         -99.871_dp   ,   -5.978_dp  ,     -2.031_dp  ,&
                          -0.32_dp    ,   -0.099_dp  ,     -0.301_dp  ,&
                           0.182_dp   ,   -1.312_dp  ,     -0.021_dp  ,&
                          -0.009_dp   ,   -0.01_dp   ,      0.061_dp  ,&
                           0.043_dp   ,   -0.002_dp  ,      0.0001_dp ,&
                           0.0023_dp  ,   -0.00091_dp,     -0.00125_dp,&
                          -0.00106_dp /)
         b4(1:31) = (/     0.002_dp   ,   -0.025_dp  ,     21.298_dp  ,&
                           0.139_dp   ,    0.045_dp  ,      0.002_dp  ,&
                           0.003_dp   ,    0.004_dp  ,      0.006_dp  ,&
                           0.016_dp   ,    0.316_dp  ,     19.158_dp  ,&
                           9.483_dp   ,    0.627_dp  ,      0.21_dp   ,&
                           0.033_dp   ,    0.01_dp   ,      0.02_dp   ,&
                           0.161_dp   ,    0.013_dp  ,      0.002_dp  ,&
                           0.001_dp   ,    0.001_dp  ,     -0.004_dp  ,&
                          -0.005_dp   ,    0.0003_dp ,      2.E-05_dp ,&
                           0.00018_dp ,    0.0001_dp ,      0.00013_dp,&
                           0.00011_dp /)
         b5(1:31) = (/    -0.011_dp   ,    0.655_dp  ,   -126.760_dp  ,&
                          -0.481_dp   ,   -0.108_dp  ,     -0.013_dp  ,&
                          -0.015_dp   ,   -0.021_dp  ,     -0.033_dp  ,&
                          -0.085_dp   ,   -1.287_dp  ,   -186.405_dp  ,&
                         -38.935_dp   ,   -3.502_dp  ,     -1.145_dp  ,&
                          -0.181_dp   ,   -0.055_dp  ,     -0.067_dp  ,&
                          -2.566_dp   ,    0.741_dp  ,     -0.012_dp  ,&
                          -0.005_dp   ,   -0.004_dp  ,      0.003_dp  ,&
                           0.026_dp   ,   -0.002_dp  ,      0.00035_dp,&
                           0.00074_dp ,   -0.00053_dp,     -0.00074_dp,&
                          -0.00063_dp /)
         b6(1:31) = (/     0.281_dp   ,  -12.998_dp  ,   3736.652_dp,&
                           5.87_dp    ,    0.713_dp  ,      0.302_dp,&
                           0.358_dp   ,    0.519_dp  ,      0.818_dp,&
                           2.202_dp   ,   21.959_dp  ,  10326.54_dp,&
                         591.792_dp   ,   86.076_dp  ,     26.99_dp,&
                           4.283_dp   ,    1.265_dp  ,      0.965_dp,&
                         177.407_dp   ,  -17.458_dp  ,      0.269_dp,&
                           0.125_dp   ,    0.089_dp  ,      1.423_dp,&
                          -0.287_dp   ,    0.038_dp  ,      0.027_dp,&
                           0.0331_dp  ,    0.013_dp  ,      0.017_dp,&
                           0.01499_dp /)
         kappap(id_co2) = b0(iband)*(tmp**6._dp)*1.e-20_dp + &
                b1(iband)*(tmp**5._dp)*1.e-17_dp + &
                b2(iband)*(tmp**4._dp)*1.e-12_dp + &
                b3(iband)*(tmp**3._dp)*1.e-8_dp  + &
                b4(iband)*(tmp**2._dp)*1.e-4_dp  + &
                b5(iband)*tmp*1.e-2_dp + b6(iband)*1.e-1_dp

!      case('bordbar2019')                                               !Borbar et al., 2019 (AE)
!         !Mount the array with the coefficients
!         selectcase(trim(wbm_mr_spec))
!         case('h2o')
!            a1(1:10) = (/  1.985829_dp,  1.209332_dp,  18.33313_dp, &
!                           19.79002_dp, -1.142476_dp,  20.84968_dp, &
!                           18.19244_dp,  14.54572_dp,  15.69315_dp, &
!                           11.07610_dp /)
!            a2(1:10) = (/ -0.169898_dp, -0.234023_dp, -3.842048_dp, &
!                          -3.751522_dp,  0.217618_dp, -4.122157_dp, &
!                          -3.654740_dp, -4.902579_dp, -6.305280_dp, &
!                          -2.645757_dp /)
!            a3(1:10) = (/ -0.648576_dp, -0.267457_dp, -2.470255_dp, &
!                          -3.361840_dp,  0.120564_dp, -3.257602_dp, &
!                          -2.814294_dp,  0.895729_dp,  3.299291_dp, &
!                          -0.735311_dp /)
!            a4(1:10) = (/ -0.007143_dp,  0.016826_dp,  0.213103_dp, &
!                           0.182403_dp, -0.003124_dp,  0.214057_dp, &
!                           0.216805_dp,  0.414010_dp,  0.621025_dp, &
!                           0.265270_dp /)
!            a5(1:10) = (/ -0.003844_dp, -0.000063_dp,  0.055931_dp, &
!                           0.080064_dp, -0.002225_dp,  0.083781_dp, &
!                           0.133542_dp,  0.001129_dp,  0.048814_dp, &
!                           0.382563_dp /)
!            a6(1:10) = (/  0.081536_dp,  0.026499_dp,  0.260911_dp, &
!                           0.356732_dp, -0.021402_dp,  0.339288_dp, &
!                           0.229317_dp, -0.145340_dp, -0.568111_dp, &
!                           -0.260229_dp /)
!         case('co2')
!            a1(1:10) = (/  0.060643_dp,  0.501492_dp, -0.068305_dp, &
!                           9.272950_dp,  61.62775_dp,  0.404750_dp, &
!                           0.015454_dp,  0.173689_dp,  21.84088_dp, &
!                           0.998090_dp /)
!            a2(1:10) = (/ -0.014286_dp, -0.098930_dp,  0.016185_dp, &
!                          -2.656498_dp, -24.72399_dp, -0.125082_dp, &
!                          -0.004939_dp, -0.379901_dp, -7.039825_dp, &
!                          -0.764430_dp /)
!            a3(1:10) = (/ -0.001099_dp, -0.071983_dp, -0.000042_dp, &
!                          -0.825930_dp,  12.61452_dp, -0.000306_dp, &
!                           0.000014_dp,  0.282935_dp,  0.275633_dp, &
!                           0.546698_dp /)
!            a4(1:10) = (/  0.000840_dp,  0.005336_dp, -0.000738_dp, &
!                           0.248210_dp,  2.734730_dp,  0.009740_dp, &
!                           0.000395_dp,  0.060332_dp,  0.707394_dp, &
!                           0.105802_dp /)
!            a5(1:10) = (/ -0.000002_dp, -0.000642_dp, -0.000048_dp, &
!                           0.046637_dp,  1.298749_dp, -0.000028_dp, &
!                           0._dp      , -0.003043_dp,  0.275312_dp, &
!                           0.009654_dp /)
!            a6(1:10) = (/  0.000142_dp,  0.009437_dp, -0.000023_dp, &
!                           0.008587_dp, -3.179420_dp,  0.000027_dp, &
!                          -0.000002_dp, -0.049357_dp, -0.404064_dp, &
!                          -0.106379_dp /)
!         case('1/8')
!            a1(1:10) = (/  0.199322_dp,  0.460116_dp,  2.567918_dp, &
!                           10.55948_dp,  56.53518_dp,  4.556402_dp, &
!                           2.413270_dp,  3.523753_dp,  23.16641_dp, &
!                           1.890447_dp /)
!            a2(1:10) = (/  0.002452_dp, -0.066706_dp, -0.485771_dp, &
!                          -2.613615_dp, -22.67026_dp, -0.988243_dp, &
!                          -0.407942_dp, -1.456466_dp, -7.684319_dp, &
!                          -0.937951_dp /)
!            a3(1:10) = (/ -0.081722_dp, -0.120368_dp, -0.460969_dp, &
!                          -1.550354_dp,  11.45638_dp, -0.605158_dp, &
!                          -0.546100_dp,  0.436200_dp,  0.816435_dp, &
!                           0.448627_dp /)
!            a4(1:10) = (/ -0.003181_dp,  0.002554_dp,  0.023309_dp, &
!                           0.219032_dp,  2.512294_dp,  0.056484_dp, &
!                           0.019668_dp,  0.146622_dp,  0.768879_dp, &
!                           0.131430_dp /)
!            a5(1:10) = (/ -0.000569_dp, -0.001168_dp,  0.003067_dp, &
!                           0.043157_dp,  1.195645_dp,  0.006011_dp, &
!                           0.010303_dp, -0.002814_dp,  0.270188_dp, &
!                           0.067787_dp /)
!            a6(1:10) = (/  0.010554_dp,  0.015310_dp,  0.054005_dp, &
!                           0.107953_dp, -2.908978_dp,  0.069641_dp, &
!                           0.055666_dp, -0.073248_dp, -0.473434_dp, &
!                          -0.151218_dp /)
!         case('1/4')
!            a1(1:10) = (/  0.331540_dp,  0.547611_dp,  4.245958_dp, &
!                           11.20230_dp,  52.42621_dp,  6.345440_dp, &
!                           4.112496_dp,  5.181861_dp,  22.65398_dp, &
!                           2.620030_dp /)
!            a2(1:10) = (/  0.001619_dp, -0.083501_dp, -0.842435_dp, &
!                          -2.576854_dp, -21.00060_dp, -1.337306_dp, &
!                          -0.764805_dp, -1.991517_dp, -7.623501_dp, &
!                          -1.075018_dp /)
!            a3(1:10) = (/ -0.144560_dp, -0.144379_dp, -0.680102_dp, &
!                          -1.919555_dp,  10.50136_dp, -0.893243_dp, &
!                          -0.781350_dp,  0.533868_dp,  1.042302_dp, &
!                           0.359544_dp /)
!            a4(1:10) = (/ -0.004914_dp,  0.003801_dp,  0.043555_dp, &
!                           0.201280_dp,  2.330265_dp,  0.074528_dp, &
!                           0.042099_dp,  0.189612_dp,  0.765537_dp, &
!                           0.145443_dp /)
!            a5(1:10) = (/ -0.000999_dp, -0.001309_dp,  0.008076_dp, &
!                           0.043926_dp,  1.110232_dp,  0.013760_dp, &
!                           0.022972_dp, -0.002467_dp,  0.263117_dp, &
!                           0.103480_dp /)
!            a6(1:10) = (/  0.018540_dp,  0.017803_dp,  0.076821_dp, &
!                           0.158788_dp, -2.685676_dp,  0.099064_dp, &
!                           0.073043_dp, -0.088538_dp, -0.494163_dp, &
!                          -0.170703_dp /)
!         case('1/2')
!            a1(1:10) = (/  0.579222_dp,  0.696623_dp,  6.660610_dp, &
!                           12.02059_dp,  45.96635_dp,  8.873031_dp, &
!                           6.543768_dp,  7.247368_dp,  21.39730_dp, &
!                           3.663677_dp /)
!            a2(1:10) = (/ -0.016965_dp, -0.116735_dp, -1.357929_dp, &
!                          -2.498614_dp, -18.37492_dp, -1.828428_dp, &
!                          -1.271468_dp, -2.653061_dp, -7.354324_dp, &
!                          -1.243549_dp /)
!            a3(1:10) = (/ -0.232978_dp, -0.171070_dp, -0.986873_dp, &
!                          -2.424063_dp,  9.002218_dp, -1.298180_dp, &
!                          -1.120598_dp,  0.647058_dp,  1.315303_dp, &
!                           0.178492_dp /)
!            a4(1:10) = (/ -0.006102_dp,  0.006394_dp,  0.072886_dp, &
!                           0.173683_dp,  2.043978_dp,  0.099676_dp, &
!                           0.073450_dp,  0.242334_dp,  0.743811_dp, &
!                           0.160816_dp /)
!            a5(1:10) = (/ -0.001557_dp, -0.001303_dp,  0.016014_dp, &
!                           0.045697_dp,  0.976053_dp,  0.025613_dp, &
!                           0.041902_dp, -0.001979_dp,  0.248349_dp, &
!                           0.152501_dp /)
!            a6(1:10) = (/  0.029676_dp,  0.020184_dp,  0.108338_dp, &
!                           0.228792_dp, -2.335062_dp,  0.140053_dp, &
!                           0.098248_dp, -0.106229_dp, -0.511823_dp, &
!                          -0.188653_dp /)
!         case('1')
!            a1(1:10) = (/  0.919122_dp,  0.877775_dp,  9.613162_dp, &
!                           12.89373_dp,  37.16659_dp,  11.92844_dp, &
!                           9.502207_dp,  9.472530_dp,  19.34950_dp, &
!                           4.952210_dp /)
!            a2(1:10) = (/ -0.050875_dp, -0.158920_dp, -1.987588_dp, &
!                          -2.377743_dp, -14.80333_dp, -2.417743_dp, &
!                          -1.881765_dp, -3.357394_dp, -6.853713_dp, &
!                          -1.433418_dp /)
!            a3(1:10) = (/ -0.339151_dp, -0.197953_dp, -1.360774_dp, &
!                          -2.992425_dp,  6.976153_dp, -1.792236_dp, &
!                          -1.542546_dp,  0.751424_dp,  1.623699_dp, &
!                          -0.075431_dp /)
!            a4(1:10) = (/ -0.006720_dp,  0.009788_dp,  0.108578_dp, &
!                           0.138231_dp,  1.655044_dp,  0.129435_dp, &
!                           0.110609_dp,  0.297754_dp,  0.701101_dp, &
!                           0.176454_dp /)
!            a5(1:10) = (/ -0.002178_dp, -0.001119_dp,  0.026032_dp, &
!                           0.049245_dp,  0.794430_dp,  0.040327_dp, &
!                           0.065159_dp, -0.001391_dp,  0.224327_dp, &
!                           0.210895_dp /)
!            a6(1:10) = (/  0.042970_dp,  0.022207_dp,  0.146709_dp, &
!                           0.307937_dp, -1.860952_dp,  0.190125_dp, &
!                           0.130306_dp, -0.122498_dp, -0.523703_dp, &
!                          -0.204262_dp /)
!         case('2')
!            a1(1:10) = (/  1.271382_dp,  1.043090_dp,  12.53369_dp, &
!                           13.72157_dp,  27.08033_dp,  14.92917_dp, &
!                           12.41867_dp,  11.43849_dp,  16.87678_dp, &
!                           6.291047_dp /)
!            a2(1:10) = (/ -0.089230_dp, -0.198014_dp, -2.609370_dp, &
!                          -2.273590_dp, -10.72507_dp, -2.993108_dp, &
!                          -2.479104_dp, -3.970542_dp, -6.218208_dp, &
!                          -1.629235_dp /)
!            a3(1:10) = (/ -0.443173_dp, -0.221590_dp, -1.731498_dp, &
!                          -3.464707_dp,  4.696161_dp, -2.282060_dp, &
!                          -1.965455_dp,  0.825736_dp,  1.945330_dp, &
!                          -0.335575_dp /)
!            a4(1:10) = (/ -0.006974_dp,  0.013038_dp,  0.143698_dp, &
!                           0.105218_dp,  1.212122_dp,  0.158181_dp, &
!                           0.146590_dp,  0.345208_dp,  0.645129_dp, &
!                           0.191338_dp /)
!            a5(1:10) = (/ -0.002755_dp, -0.000830_dp,  0.036029_dp, &
!                           0.055098_dp,  0.588806_dp,  0.054901_dp, &
!                           0.088115_dp, -0.000736_dp,  0.192917_dp, &
!                           0.268061_dp /)
!            a6(1:10) = (/  0.055954_dp,  0.023742_dp,  0.184815_dp, &
!                           0.373611_dp, -1.326135_dp,  0.239899_dp, &
!                           0.162958_dp, -0.134068_dp, -0.530026_dp, &
!                          -0.217563_dp /)
!         case default
!            call shutdown('wbm_mr_spec not available')
!         endselect
         
!         !Compute the absorption coefficient
!         pl = (xxco2 + xxh2o)*tp*path_length
!         if (pl.lt.0.1_dp) then
!            pl2 = 0.1_dp
!            k2 = a1(iband) + a2(iband)*log(tmp) + &
!               a3(iband)*log(pl2) + a4(iband)*(log(tmp))**2._dp + &
!               a5(iband)*(log(pl2))**2._dp + a6(iband)*log(tmp)*log(pl2)
!            box_kappa_func = k2*pl/pl2
!         else
!         !if (pl.gt.15_dp) pl = 15._dp
!            box_kappa_func = a1(iband) + a2(iband)*log(tmp) + &
!               a3(iband)*log(pl) + a4(iband)*(log(tmp))**2._dp + &
!               a5(iband)*(log(pl))**2._dp + a6(iband)*log(tmp)*log(pl)
!         endif

      case('paul')
         !Mount the array with the coefficients
         aco2(1:5) = (/ 0.6587_dp, 5.12_dp, 0.6359_dp, 0.07567_dp,&
            0.00023387_dp /)
         bco2(1:5) = (/ -0.004613_dp, -0.002712_dp, -0.006208_dp,&
            -0.004972_dp, 0._dp /)
         cco2(1:5) = (/ 0.02394_dp, 3.146_dp, 0.03651_dp, 0.03519_dp,&
            0._dp /)
         dco2(1:5) = (/ -0.000468_dp, -0.0005468_dp, -0.001691_dp,&
            -0.0003229_dp, 0._dp /)
       
         ah2o(1:5) = (/ 0.489_dp, 0.004_dp, 0.0018_dp, 0.09007_dp,&
            0.0027_dp /)
         bh2o(1:5) = (/ -0.003077_dp, 0._dp, 0._dp, -0.001306_dp,&
            0._dp /)
         ch2o(1:5) = (/ 0.1124_dp, 0._dp, 0._dp, 0.05733_dp, 0._dp /)
         dh2o(1:5) = (/ -0.0003568_dp, 0._dp, 0._dp, -0.000418_dp,&
            0._dp /)
         
         aco(1:5)  = (/ 0.0145_dp, 0.1128_dp, 5.3533e-5_dp,&
            3.3125e-5_dp, 0.00022833_dp /)
         bco(1:5)  = (/ 0._dp, 0._dp, 0._dp, 0._dp, 0._dp /)
         cco(1:5)  = (/ 0._dp, 0._dp, 0._dp, 0._dp, 0._dp /)
         dco(1:5)  = (/ 0._dp, 0._dp, 0._dp, 0._dp, 0._dp /)

         !Compute the pressure absorption coefficient of each species
         kappap(id_co2) = aco2(iband)*exp(bco2(iband)*tmp) + &
            cco2(iband)*exp(dco2(iband)*tmp)
         kappap(id_h2o) = ah2o(iband)*exp(bh2o(iband)*tmp) + &
            ch2o(iband)*exp(dh2o(iband)*tmp)
         kappap(id_co ) = aco(iband)*exp(bco(iband)*tmp) + &
            cco(iband)*exp(dco(iband)*tmp)
         
         !Convert to 1/m
         kappap = kappap*100._dp
         
      case('johnson-baseline')
         !Mount the array with the coefficients
         aco2(1:4) = (/ 6.50e-3_dp, 1.51e-2_dp, 1.14_dp , 1.07e-3_dp /)
         bco2(1:4) = (/ -4.33e-3_dp, -3.16e-3_dp, -1.15e-3_dp, 0._dp /)
         cco2(1:4) = (/ 7.47e-4_dp, 1.99e-2_dp, 2.76e-1_dp, 0._dp /)
         dco2(1:4) = (/ -6.12e-4_dp, -3.87e-4_dp, -2.65e-4_dp, 0._dp /)
         
         ah2o(1:4) = (/ 1.12e-2_dp, 7.26e-2_dp, 1.62e-2_dp, 1.65e-1_dp /)
         bh2o(1:4) = (/ -2.16e-3_dp, -7.92e-4_dp, 0._dp, -1.06e-3_dp /)
         ch2o(1:4) = (/ 4.12e-3_dp, 6.67e-4_dp, 0._dp, 3.83e-2_dp /)
         dh2o(1:4) = (/ -1.86e-4_dp, 5.99e-4_dp, 0._dp, -6.32e-5_dp /)
         
         aco(1:4) = (/ 2.18e-6_dp, 2.54e-4_dp, 4.09e-2_dp, 2.31e-5_dp /)
         bco(1:4) = (/ 0._dp, 0._dp, 0._dp, 0._dp /)
         cco(1:4) = (/ 0._dp, 0._dp, 0._dp, 0._dp /)
         dco(1:4) = (/ 0._dp, 0._dp, 0._dp, 0._dp /)

         !Compute the pressure absorption coefficient of each species
         kappap(id_co2) = aco2(iband)*exp(bco2(iband)*tmp) + &
            cco2(iband)*exp(dco2(iband)*tmp)
         kappap(id_h2o) = ah2o(iband)*exp(bh2o(iband)*tmp) + &
            ch2o(iband)*exp(dh2o(iband)*tmp)
         kappap(id_co ) = aco(iband)*exp(bco(iband)*tmp) + &
            cco(iband)*exp(dco(iband)*tmp)
         
         !Convert to 1/m
         kappap = kappap*100._dp
            
      case('johnson-improved')
         !Mount the array with the coefficients
         aco2(1:6) = (/ 5.8500e-1_dp, 1.5000e-3_dp, 5.1200_dp,&
            6.3590e-1_dp, 7.5670e-2_dp, 2.3387e-4_dp /)
         bco2(1:6) = (/ -3.7130e-3_dp, 0._dp, -2.7120e-3_dp,&
            -6.2080e-3_dp, -4.9720e-3_dp, 0._dp /)
         cco2(1:6) = (/ 1.6670e-1_dp, 0._dp, 3.1460_dp, 3.6510e-2_dp,&
            3.5190e-2_dp, 0._dp /)
         dco2(1:6) = (/ -5.6390e-4_dp, 0._dp, -5.4680e-4_dp,&
            -1.6910e-3_dp, -3.2290e-4_dp, 0._dp /)
         
         ah2o(1:6) = (/ 6.4650e-1_dp, 2.2070e-1_dp, 4.0000e-3_dp,&
            1.8000e-3_dp, 9.0070e-2_dp, 2.7000e-3_dp /)
         bh2o(1:6) = (/ -3.7420e-3_dp, -1.9990e-3_dp, 0._dp, 0._dp,&
            -1.3060e-3_dp, 0._dp /)
         ch2o(1:6) = (/ 2.1320e-1_dp, 6.6880e-2_dp, 0._dp, 0._dp,&
            5.7330e-2_dp, 0._dp /)
         dh2o(1:6) = (/ -1.2190e-4_dp, -3.3460e-4_dp, 0._dp, 0._dp,&
            -4.1800e-4_dp, 0._dp/)
         
         aco(1:6)  = (/ 6.5447e-6_dp, 1.7200e-2_dp, 1.1280e-1_dp,&
            5.3533e-5_dp, 3.3125e-5_dp, 2.2833e-4_dp /)
         bco(1:6)  = (/ 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp /)
         cco(1:6)  = (/ 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp /)
         dco(1:6)  = (/ 0._dp, 0._dp, 0._dp, 0._dp, 0._dp, 0._dp /)
      
         !Compute the pressure absorption coefficient of each species
         kappap(id_co2) = aco2(iband)*exp(bco2(iband)*tmp) + &
            cco2(iband)*exp(dco2(iband)*tmp)
         kappap(id_h2o) = ah2o(iband)*exp(bh2o(iband)*tmp) + &
            ch2o(iband)*exp(dh2o(iband)*tmp)
         kappap(id_co ) = aco(iband)*exp(bco(iband)*tmp) + &
            cco(iband)*exp(dco(iband)*tmp)
         
         !Convert to 1/m
         kappap = kappap*100._dp
      
      case('external')                                                  !Interpolate the data
         do isp=1,size(xs)
            !Surrogate names
            ntg = wbm_data_ntg(isp)
            nxs = wbm_data_nxs(isp)

            !Upper and lower temperature indexes
            ltg = locate(wbm_data_tg(:,isp),tmp,ntg)                    !Locate lower index
            ltg = max(1,min(ntg-1,ltg)); utg = min(ltg+1,ntg)           !Correct lower index, compute upper index
            if (tmp.lt.wbm_data_tg(1,isp))    utg = ltg                 !This is to prevent extrapolations
            if (tmp.gt.wbm_data_tg(ntg,isp))  ltg = utg                 !  (instead, simply take the value at the extreme)

            !!Upper and lower mole fraction indexes 
            lxs = locate(wbm_data_xs(:,isp),xs(isp),nxs)                !Locate lower index
            lxs = max(1,min(nxs-1,lxs)); uxs = min(lxs+1,nxs)           !Correct lower index, compute upper index
            if (xs(isp).lt.wbm_data_xs(1,isp))    uxs = lxs             !This is to prevent extrapolations
            if (xs(isp).gt.wbm_data_xs(nxs,isp))  lxs = uxs             !  (instead, simply take the value at the extreme)

            !Initial values for the interpolation on mole fraction
            Qxs = 0._dp
            Rxs = (xs(isp) - wbm_data_xs(lxs,isp))/&
                  (wbm_data_xs(uxs,isp) - wbm_data_xs(lxs,isp) + small)
            if (lxs.eq.uxs) Rxs = 0._dp                                 !If only one mole fraction value is provided,
                                                                        !  do not interpolate in mole fraction
            molfrac_loop: do ixs=lxs,uxs
               !Initial values for the interpolation on temperature
               Qtg = 0._dp
               Rtg = (tmp - wbm_data_tg(ltg,isp))/&
                  (wbm_data_tg(utg,isp) - wbm_data_tg(ltg,isp) + small)
               if (ltg.eq.utg) Rtg = 0._dp                              !If only one temperature value is provided, 
                                                                        !  do not interpolate in temperature
               !Interpolate on temperature
               tmp_loop: do itg=ltg,utg
                  Rtg = 1._dp - Rtg
                  Qtg = Qtg + Rtg*wbm_kappa(itg,ixs,iband,isp)
               enddo tmp_loop

               !Interpolate on mole fraction
               Rxs = 1._dp - Rxs
               Qxs = Qxs + Rxs*Qtg
            enddo molfrac_loop
            kappap(isp) = Qxs
         enddo
      
      case default
         call shutdown('box_kappa_func: WBM not available')
      endselect

      !-----------------------------------------------------------------
      !Compute the absorption coefficient of the mixture
      !-----------------------------------------------------------------
      kappa = 0._dp
      do isp=1,size(xs)
         if ((isp.ne.id_co2).and.(isp.ne.id_h2o).and.(isp.ne.id_co)) &
            cycle
         kappa = kappa + kappap(isp)*xs(isp)*tp
      enddo
      box_kappa_func = kappa

   endfunction box_kappa_func

endmodule wbm_functions
