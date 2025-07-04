module gg_functions

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp
   implicit none
   
contains
   
   !====================================================================
   !Function to compute the absorption coefficient 
   !of the gray medium via different correlations
   !====================================================================
   real(dp) function gg_kappa_func(ttmp,tp,xxs,wsgg_spec)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use gg_parameters, only: get_interpolated_kp,which_gg_model
      use global_parameters, only: id_ch4,id_co,id_co2,id_h2o,id_soot,&
                                   number_of_species,soot_constant
      implicit none
      character(*),optional :: wsgg_spec
      character(100) :: wsgg_spec_aux
      integer :: i,nsp
      real(dp),intent(in) :: ttmp,tp,xxs(:)
      real(dp) :: xxch4,xxco,xxco2,xxh2o,xxsoot
      real(dp) :: kch4,kco,kco2,kh2o,ksoot,kmix
      real(dp) :: c_ch4(5),c_co(5),c_co2(6),c_h2o(6),c_soot(3)
      
      !-----------------------------------------------------------------
      !Default parameters
      !-----------------------------------------------------------------
      if(.not.present(wsgg_spec)) then
         wsgg_spec_aux = 'null'
      else
         wsgg_spec_aux = wsgg_spec
      endif
      
      !-----------------------------------------------------------------
      !Concentrations of individual species
      !-----------------------------------------------------------------
      nsp = number_of_species
      xxch4  = 0._dp; if ((id_ch4.le.nsp).and.(id_ch4.gt.0))  xxch4  = xxs(id_ch4) 
      xxco   = 0._dp; if ((id_co.le.nsp).and.(id_co.gt.0))   xxco   = xxs(id_co)
      xxco2  = 0._dp; if ((id_co2.le.nsp).and.(id_co2.gt.0))  xxco2  = xxs(id_co2)
      xxh2o  = 0._dp; if ((id_h2o.le.nsp).and.(id_h2o.gt.0))  xxh2o  = xxs(id_h2o)
      xxsoot = 0._dp; if ((id_soot.le.nsp).and.(id_soot.gt.0)) xxsoot = xxs(id_soot)
      
      !-----------------------------------------------------------------
      !Selecting the model
      !-----------------------------------------------------------------
      selectcase(trim(which_gg_model))
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Barlow et al.(2002) (obtained from the TNF workshop website)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('TNF')
            !Constants of the model
            c_ch4(1:5) = (/ 6.6334_dp, -0.0035686_dp, 1.6682e-8_dp,&
                            2.5611e-10_dp, -2.6558e-14_dp /)
            if (ttmp.lt.750._dp) then
               c_co(1:5) = (/ 4.7869_dp, -0.06953_dp, 2.95775e-4_dp,&
                             -4.25732e-7_dp, 2.02894e-10_dp /)
            else
               c_co(1:5) = (/ 10.09_dp, -0.01183_dp, 4.7753e-6_dp,&
                             -5.87209e-10_dp, -2.5334e-14_dp /)
            endif
            c_co2(1:6) = (/ 1.8741e1_dp, -1.2131e2_dp, 2.7350e2_dp, &
                           -1.9405e2_dp, 5.6310e1_dp, -5.8169_dp /)
            c_h2o(1:6) = (/ -2.3093e-1_dp, 1.1239_dp, 9.4153_dp, &
                            -2.9988_dp, 5.1382e-1_dp, -1.8684e-5_dp /)
            
            !Absorption coefficient for CH4 and CO
            kch4 = 0._dp; kco = 0._dp
            do i=1,5
               kch4 = kch4 + c_ch4(i)*ttmp**(real(i-1,dp))
               kco  = kco  + c_co(i) *ttmp**(real(i-1,dp))
            enddo
            kch4 = tp*xxch4*kch4; kco = tp*xxco*kco
                          
            !Absorption coefficient for CO2 and H2O
            kco2 = 0._dp; kh2o = 0._dp
            do i=1,6
               kco2 = kco2 + c_co2(i)*((1000._dp/ttmp)**(real(i-1,dp)))
               kh2o = kh2o + c_h2o(i)*((1000._dp/ttmp)**(real(i-1,dp)))
            enddo
            kco2 = tp*xxco2*kco2; kh2o = tp*xxh2o*kh2o            
         
            !Mixture's absorption coefficient
            gg_kappa_func = kco2 + kh2o + kco + kch4

         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Kumar et al.(2002)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('kumar2002')
            !Constants of the model
            c_co2(1:5) = (/ 0.96986e3_dp, -0.58838e3_dp, 0.13289e3_dp,&
                           -0.13182e2_dp, 0.48396_dp /)
            c_h2o(1:5) = (/ 0.278713e3_dp, -0.153240e3_dp,&
                            0.321971e2_dp, -0.300870e1_dp, 0.104055_dp/)
   
            !Absorption coefficient for CO2 and H2O
            kco2 = 0._dp; kh2o = 0._dp
            do i=1,5
               kco2 = kco2 + c_co2(i)*((log(ttmp))**(real(i-1,dp)))
               kh2o = kh2o + c_h2o(i)*((log(ttmp))**(real(i-1,dp)))
            enddo
            kco2 = tp*xxco2*dexp(kco2); kh2o = tp*xxh2o*dexp(kh2o)
         
            !Absorption coefficient for soot
            ksoot = 1864.32_dp*soot_constant*xxsoot*ttmp/7._dp
            
            !Mixture's absorption coefficient
            gg_kappa_func = kco2 + kh2o + ksoot
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Cassol et al.(2015)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('cassol2015')
            !Constants of the model
            c_co2(1:6) = (/ -6.4750e-1_dp, 4.2895e-3_dp, -6.6089e-6_dp,&
                           4.4190e-9_dp, -1.3796e-12_dp, 1.6484e-16_dp/)
            c_h2o(1:6) = (/ 7.5702e-1_dp, -1.9716e-3_dp, 2.1998e-6_dp,&
                         -1.2492e-9_dp, 3.5385e-13_dp, -3.9663e-17_dp /)
            c_soot(1:3) = (/ -281.19_dp, 3.7161_dp, -6.7737e-4_dp /)
         
            !Absorption coefficient for CO2 and H2O
            kco2 = 0._dp; kh2o = 0._dp
            do i=1,6
               kco2 = kco2 + c_co2(i)*(ttmp**(real(i-1,dp)))
               kh2o = kh2o + c_h2o(i)*(ttmp**(real(i-1,dp)))
            enddo
            kco2 = tp*xxco2*kco2; kh2o = tp*xxh2o*kh2o
      
            ksoot = 0._dp
            do i=1,3
               ksoot = ksoot + c_soot(i)*(ttmp**(real(i-1,dp)))
            enddo
            ksoot = soot_constant*xxsoot*ksoot
      
            !Mixture's absorption coefficient
            kmix = kco2 + kh2o + ksoot
         
            !Converting from 1/cm to 1/m
            gg_kappa_func = kmix*100._dp
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Computed from data for the pressure
         !Planck-mean absorption coefficient
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('kp_data')
            kco = get_interpolated_kp(ttmp,'co')*tp*xxco
            kco2 = get_interpolated_kp(ttmp,'co2')*tp*xxco2
            kh2o = get_interpolated_kp(ttmp,'h2o')*tp*xxh2o
            ksoot = get_interpolated_kp(ttmp,'soot')*tp*xxsoot   
            gg_kappa_func = kco + kco2 + kh2o + ksoot
         
         case default
            call shutdown('gg_kappa_func: which_gg_model not defined')
      endselect

   endfunction gg_kappa_func

endmodule gg_functions
