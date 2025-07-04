!#######################################################################
!Module compiling WSGG parameters and correlations
!#######################################################################
module wsgg_parameters

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp
   
   !====================================================================
   !Declaration of variables
   !====================================================================
   implicit none
   character(80) :: soot_gg_correction_model                            !Model used to compute the gray soot absorption coefficient
                                                                        !  (to be added to the absorption coefficient of the gray
                                                                        !  gases of the gas phase, if soot_gg_correction = .true.)
   character(80) :: soot_wsgg_model                                     !If different soot correlations are provided in the same 
                                                                        !  work, defines which correlation to used
   character(80) :: which_wsgg_model                                    !Defines which general WSGG correlation to use                                           
   character(80) :: wsgg_ch4_spec,wsgg_co_spec,wsgg_co2_spec,&          !Defines which WSGG correlation to use for modeling
                    wsgg_h2o_spec,wsgg_soot_spec                        !  individual species (used in the superposition model)
   character(80) :: wsgg_species_spec(10)                               !Same as above, but in array form, with each position
                                                                        !  corresponding to a different species
   character(80) :: wsgg_mixture_spec                                   !Defines which WSGG correlation to use for modeling a
                                                                        !  mixture of gases (used in the superposition model)
   character(80) :: wsgg_print_file                                     !Name of the file where the gray gas properties
                                                                        !  are printed, if requested
   character(80) :: fraga_h2o_interpolation                             !In my 2021 H2O correlations, defines which approach to use
                                                                        !  to deal with varying H2O mole fractions
   integer :: degree_wsgg_polynomial                                    !For fixed-MR correlations, degree of the temperature
                                                                        !  polynomial that defines the weighting coefficient
                                                                        !  (used internally, not to be speficied by the user)
   integer :: number_wsgg_gases                                         !Number of gray gases in a correlation (used internally,
                                                                        !  not to be speficied by the user)
   integer :: number_wsgg_gases_soot                                    !For soot correlations proposed in a same work that were
                                                                        !  developed based on different number of gray gases,
                                                                        !  defines which of the correlations to use
   real(dp) :: b(1:10,1:10),kappa_p(1:10)                               !For fixed-MR correlations, array with the gray gas 
                                                                        !  polynomial and pressure absorption coefficients
                                                                        !  (used internally, not to be speficied by the user)
   real(dp) :: wsgg_Ib_lbound,wsgg_Ib_ubound                            !Bounds (in 1/m) for the Ib integration
   real(dp) :: bordbar_infty,smith_infty,yin_infty                      !For varying-MR correlation, defines the value of "infinity"
                                                                        !  use to interpolate the correlations from their upper MR
                                                                        !  bounds to a correlation for pure H2O
   real(dp) :: which_h2o_pressure                                       !For H2O correlations available for different mole fractions,
                                                                        !  defines which one to use (value, not string)
   real(dp) :: which_molar_ratio                                        !For fixed-MR correlations available for different MR values,
                                                                        !  defines defines which one to use (value, not string)
   real(dp) :: wsgg_superposition_theta                                 !Defines the value of the safety factor in the improved 
                                                                        !  superposition model
   real(dp) :: which_total_pressure                                     !For correlations available for different total pressures,
                                                                        !  defines which one to use (value, not string)
   real(dp) :: wsgg_reference_temperature                               !Reference temperature for fixed MR correlations
   logical :: soot_transparent_windows                                  !If .true., includes a clear gas for soot in the modeling
   logical :: soot_gg_correction                                        !If .true., incorporates gray soot to the modeling
   logical :: wsgg_print_properties                                     !If .true., print the properties of each gray gas into
                                                                        !  file. The name of the file is wsgg_print_file
   logical :: wsgg_bound_Ib                                             !If .true., the blackbody intensity is not computed over
                                                                        !  the entire spectrum, but rather between wsgg_Ib_lbound
                                                                        !  and wsgg_Ib_ubound
   logical :: bordbar_co2_interpolation,bordbar_h2o_interpolation       !If .true., turns on interpolation for the Bordbar 1atm
                                                                        !  correlations (these parameters can likely be removed
                                                                        !  in the future)
   logical :: wsgg_improved_superposition                               !If .true., use the improved superposition method instead
                                                                        !  of the standard method
   real(dp) :: wsgg_beam_length
   
contains   
   !====================================================================
   !Subroutine with the default parameters for the WSGG model
   !====================================================================
   !(this does not need to be included in the main 
   !code if the user knows what they are doing)
   subroutine set_default_wsgg_parameters
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use precision_parameters, only: big,small
      implicit none
      
      !-----------------------------------------------------------------
      !General parameters
      !-----------------------------------------------------------------
      which_total_pressure = 0._dp
      which_h2o_pressure = 0._dp
      which_molar_ratio = 0._dp
      which_wsgg_model = 'null'
      wsgg_reference_temperature = 1._dp
      wsgg_bound_Ib = .false.
      wsgg_Ib_lbound = small
      wsgg_Ib_ubound = big
      
      !-----------------------------------------------------------------
      !Parameters related to the superposition model
      !-----------------------------------------------------------------
      wsgg_improved_superposition = .false.
      wsgg_superposition_theta = 2._dp
      wsgg_species_spec = 'null'
      wsgg_ch4_spec     = 'null'
      wsgg_co_spec      = 'null'
      wsgg_co2_spec     = 'null'
      wsgg_h2o_spec     = 'null'
      wsgg_soot_spec    = 'null'
      wsgg_mixture_spec = 'null'
      
      !-----------------------------------------------------------------
      !Parameters related to soot
      !-----------------------------------------------------------------
      number_wsgg_gases_soot = 0
      soot_wsgg_model = 'null'
      soot_transparent_windows = .true.
      soot_gg_correction = .false.
      soot_gg_correction_model = 'null'
      
      !-----------------------------------------------------------------
      !Parameters specific to single correlations
      !-----------------------------------------------------------------
      fraga_h2o_interpolation = 'null'
      bordbar_infty = 100._dp
      smith_infty = 100._dp
      yin_infty = 100._dp     
      bordbar_co2_interpolation = .false.
      bordbar_h2o_interpolation = .false.

      !-----------------------------------------------------------------
      !Misc parameters
      !-----------------------------------------------------------------
      wsgg_beam_length      = 1._dp
      wsgg_print_properties = .false.
      wsgg_print_file       = 'null'

   endsubroutine set_default_wsgg_parameters

   !====================================================================
   !Subroutine compiling all WSGG correlations
   !====================================================================
   subroutine get_wsgg_correlations(wsgg_spec,species_spec)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      character(*),intent(in) :: wsgg_spec
      character(*),optional :: species_spec
      character(100) :: species_spec_aux
      character(80) :: msg
      
      !-----------------------------------------------------------------
      !Set defaults
      !-----------------------------------------------------------------
      species_spec_aux = species_spec
      if (.not.present(species_spec)) species_spec_aux = 'mixture'
      
      !-----------------------------------------------------------------
      !Defining WSGG parameters
      !-----------------------------------------------------------------
      select case(trim(wsgg_spec))           
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Smith et al. (1982)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('smith1982')         
            !Total number of gray gases
            number_wsgg_gases = 3
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 3
            
            if (trim(species_spec_aux).eq.'co2') then
               !Pressure-absorption coefficients
               kappa_p(1:3) = (/ 0.3966_dp, 15.64_dp, 394.3_dp /)

               !Polynomial coefficients
               b(1,1:4) = (/  0.4334e-1_dp,  2.6200e-4_dp,&
                             -1.5600e-7_dp,  2.5650e-11_dp  /)
               b(2,1:4) = (/ -0.4814e-1_dp,  2.8220e-4_dp,&
                             -1.7940e-7_dp,  3.2740e-11_dp  /)
               b(3,1:4) = (/  0.5492e-1_dp,  0.1087e-4_dp,&
                             -0.3500e-7_dp,  0.9123e-11_dp /)
 
            elseif (trim(species_spec_aux).eq.'h2o') then       !p_w -> 0 atm
               if (which_h2o_pressure.eq.0._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.4098_dp, 6.325_dp, 120.5_dp /)
                  
                  !Polynomial coefficients
                  b(1,1:4) = (/  5.977e-1_dp, -5.119e-4_dp,&
                                 3.042e-7_dp, -5.564e-11_dp /)
                  b(2,1:4) = (/  0.5677e-1_dp, 3.333e-4_dp,&
                                -1.9670e-7_dp, 2.718e-11_dp /)
                  b(3,1:4) = (/  1.800e-1_dp, -2.334e-4_dp,&
                                 1.008e-7_dp, -1.454e-11_dp /)
               
               elseif (which_h2o_pressure.eq.1._dp) then                !p_w = 1 atm
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.4496_dp, 7.113_dp, 119.7_dp /)
                  
                  !Polynomial coefficients
                  b(1,1:4) = (/  6.3240e-1_dp, -8.358e-4_dp,&
                                 6.135e-7_dp, -13.03e-11_dp /)
                  b(2,1:4) = (/ -0.2016e-1_dp,  7.145e-4_dp,&
                                -5.212e-7_dp,  9.868e-11_dp /)
                  b(3,1:4) = (/  3.5000e-1_dp, -5.040e-4_dp,&
                                 2.425e-7_dp, -3.888e-11_dp /)
               else
                  write(msg,'(a,f3.1,a,a)') 'H2O partial pressure ',&
                     which_h2o_pressure,' atm not available for &
                     &correlation ',trim(wsgg_spec)
                  call shutdown(msg)
               endif          
               
            else
               if (which_molar_ratio.eq.1._dp) then                     !p_w/p_c = 1
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.4303_dp, 7.055_dp, 178.1_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  5.150e-1_dp ,  7.749e-2_dp ,&
                                 1.907e-1_dp  /)                        !0th degree
                  b(1:3,2) = (/ -2.303e-4_dp ,  3.399e-4_dp ,&
                                -1.824e-4_dp  /)                        !1st degree
                  b(1:3,3) = (/  9.779e-8_dp , -2.297e-7_dp ,&
                                 5.608e-8_dp  /)                        !2nd degree
                  b(1:3,4) = (/ -1.494e-11_dp,  3.770e-11_dp,& 
                                -5.122e-12_dp /)                        !3rd degree
               
               elseif (which_molar_ratio.eq.2._dp) then                 !p_w/p_c = 2
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.4201_dp, 6.516_dp, 131.9_dp /)
         
                  !Polynomial coefficients
                  b(1:3,1) = (/  6.508e-1_dp, -0.2504e-1_dp,&
                                 2.718e-1_dp  /)                        !0th degree
                  b(1:3,2) = (/ -5.551e-4_dp,  6.112e-4_dp ,& 
                                -3.118e-4_dp  /)                        !1st degree
                  b(1:3,3) = (/  3.029e-7_dp, -3.882e-7_dp ,&  
                                 1.221e-7_dp  /)                        !2nd degree
                  b(1:3,4) = (/ -5.353e-11_dp, 6.528e-11_dp,& 
                                -1.612e-11_dp /)                        !3rd degree
                     
               else                                                     !Stop and print error message
                  write(msg,'(a,f3.1,a,a)') &
                     'Partial pressure ratio ',which_molar_ratio,&
                     ' not available for correlation ',trim(wsgg_spec)
                  call shutdown(msg)
               
               endif
            
            endif
            
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Soufiani & Djavdan (1994)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('soufiani1994')   
            !Total number of gray gases
            number_wsgg_gases = 3
         
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 5

            if (which_molar_ratio.eq.2._dp) then                        !p_w/p_c = 2
               !Pressure-absorption coefficients
               kappa_p(1:3) = (/ 1.2531_dp, 8.4258_dp, 87.064_dp /)

               !Polynomial coefficients
               b(1:3,1) = (/  1.6879e-1,   4.9577e-2,   2.7890e-1  /)   !0th degree
               b(1:3,2) = (/  2.5682e-4,   9.3954e-4,  -5.1265e-4  /)   !1st degree
               b(1:3,3) = (/  9.5161e-8,  -1.6416e-6,   6.7320e-7  /)   !2nd degree
               b(1:3,4) = (/ -3.1660e-10,  1.1478e-9,  -5.1488e-10 /)   !3rd degree
               b(1:3,5) = (/  1.4834e-13, -3.7600e-13,  1.8887e-13 /)   !4th degree
               b(1:3,6) = (/ -2.1560e-17,  4.7503e-17, -2.5856e-17 /)   !5th degree

            else                                                         !Stop and print error message
               write(msg,'(a,f3.1,a,a)') 'Partial pressure ratio ',&
                  which_molar_ratio,' not available for correlation ',&
                  trim(wsgg_spec)
               call shutdown(msg)
            endif

         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Bahador & Sunden (2008)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('bahador2008')
            !Total number of gray gases
            number_wsgg_gases = 3
         
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 3         
         
            if (which_molar_ratio.eq.1._dp) then
               if (which_total_pressure.eq.1._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.2924_dp, 4.5820_dp, 97.8332_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  5.2421e-1_dp ,  1.1647e-1_dp ,&
                                 2.2616e-1_dp  /)                       !0th degree
                  b(1:3,2) = (/ -1.5856e-4_dp ,  2.8420e-4_dp ,&  
                                -1.9445e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  6.8306e-8_dp , -2.1488e-7_dp ,& 
                                 5.6080e-8_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -1.4940e-11_dp,  3.7700e-11_dp,& 
                                -5.1220e-12_dp /)                       !3rd degree   

               elseif (which_total_pressure.eq.2.5_dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.2959_dp, 5.0180_dp, 115.7979_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  6.8211e-1_dp ,  1.6502e-1_dp ,&
                                 2.1804e-1_dp  /)                       !0th degree
                  b(1:3,2) = (/ -3.6623e-4_dp ,  2.8854e-4_dp ,&
                                -1.9193e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  1.3387e-7_dp , -2.2503e-7_dp ,&
                                 5.6080e-8_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -1.4940e-11_dp,  3.7700e-11_dp,&
                                -5.1220e-12_dp /)                       !3rd degree   
               
               elseif (which_total_pressure.eq.5._dp) then      
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.3102_dp, 5.3389_dp, 118.5489_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  3.5365e-1_dp ,  3.6378e-1_dp ,&
                                 2.1912e-1_dp  /)                       !0th degree
                  b(1:3,2) = (/  8.8373e-5_dp ,  5.8208e-5_dp ,&
                                -1.9273e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/ -4.5279e-10_dp, -1.6238e-7_dp ,& 
                                 5.6080e-8_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -1.4940e-11_dp,  3.7700e-11_dp,&
                                -5.1220e-12_dp /)                       !3rd degree   

               elseif (which_total_pressure.eq.7.5_dp) then   
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.3166_dp, 5.4404_dp, 121.2707_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  3.3797e-1_dp ,  4.0883e-1_dp ,&
                                 2.1783e-1_dp  /)                       !0th degree
                  b(1:3,2) = (/  1.1362e-4_dp ,  1.8948e-5_dp ,& 
                                -1.9233e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/ -7.5256e-9_dp , -1.5378e-7_dp ,&
                                 5.6080e-8_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -1.4940e-11_dp,  3.7700e-11_dp,&
                                -5.1220e-12_dp /)                       !3rd degree   
         
               elseif (which_total_pressure.eq.10._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.3208_dp, 5.4921_dp, 122.7539_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  3.2978e-1_dp ,  4.3950e-1_dp ,&
                                 2.1705e-1_dp  /)                       !0th degree
                  b(1:3,2) = (/  1.2866e-4_dp , -7.9118e-6_dp ,& 
                                -1.9207e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/ -1.1953e-8_dp , -1.4788e-7_dp ,&
                                 5.6080e-8_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -1.4940e-11_dp,  3.7700e-11_dp,&
                                -5.1220e-12_dp /)                       !3rd degree   
         
               elseif (which_total_pressure.eq.12.5_dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.3120_dp, 5.3784_dp, 122.5504_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  3.2238e-1_dp ,  4.6432e-1_dp ,&
                                 2.1759e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/  1.3753e-4_dp , -2.2944e-5_dp ,&
                                -1.9230e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/ -1.3636e-8_dp , -1.4558e-7_dp ,&
                                 5.6080e-8_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -1.4940e-11_dp,  3.7700e-11_dp,&
                                -5.1220e-12_dp /)                       !3rd degree   

               elseif (which_total_pressure.eq.15._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.3260_dp, 5.5372_dp, 124.2125_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  3.1885e-1_dp ,  4.8302e-1_dp ,&
                                 2.1621e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/  1.4928e-4_dp , -4.6337e-5_dp ,&
                                -1.9179e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/ -1.8156e-8_dp , -1.3937e-7_dp ,&
                                 5.6080e-8_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -1.4940e-11_dp,  3.7700e-11_dp,&
                                -5.1220e-12_dp /)                       !3rd degree   
         
               elseif (which_total_pressure.eq.20._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.3233_dp, 5.5180_dp, 125.4961_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  3.2718e-1_dp ,  5.0175e-1_dp ,&
                                 2.1537e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/  1.4189e-4_dp , -5.4062e-5_dp ,&
                                -1.9149e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/ -1.5774e-8_dp , -1.3935e-7_dp ,&
                                 5.6080e-8_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -1.4940e-11_dp,  3.7700e-11_dp,&
                                -5.1220e-12_dp /)                       !3rd degree   
         
               else                                                     !Stop and print error message
                  write(msg,'(a,f3.1,a,a)') &
                     'Total pressure ',which_total_pressure,&
                     ' not available for correlation ',trim(wsgg_spec)
                  call shutdown(msg)
               endif
            
            elseif (which_molar_ratio.eq.2._dp) then      
               if (which_total_pressure.eq.1._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.2926_dp, 4.4201_dp, 76.2197_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  4.6406e-1_dp ,  8.9179e-2_dp ,&
                                 3.0884e-01_dp /)                       !0th degre
                  b(1:3,2) = (/ -2.2005e-4_dp ,  4.5355e-4_dp ,&
                                -3.2474e-04_dp /)                       !1st degree   
                  b(1:3,3) = (/  2.0267e-7_dp , -3.4437e-7_dp ,&
                                 1.2210e-07_dp /)                       !2nd degree   
                  b(1:3,4) = (/ -5.3530e-11_dp,  6.5280e-11_dp,&
                                -1.6120e-11_dp /)                       !3rd degree   

               elseif (which_total_pressure.eq.2.5_dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.32164_dp, 5.19426_dp,&
                                    87.78467_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  4.0114e-1_dp ,  2.0600e-1_dp ,&
                                 3.0040e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/ -1.3702e-4_dp ,  3.5383e-4_dp ,&
                                -3.2221e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  1.8170e-7_dp , -3.2309e-7_dp ,&
                                 1.2210e-7_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -5.3530e-11_dp,  6.5280e-11_dp,&
                                -1.6120e-11_dp /)                       !3rd degree   
                     
               elseif (which_total_pressure.eq.5._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.3403_dp, 5.5908_dp, 93.8217_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  3.7498e-1_dp ,  2.7903e-1_dp ,&
                                 2.9564e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/ -9.9314e-5_dp ,  2.9151e-4_dp ,&
                                -3.2066e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  1.7188e-7_dp , -3.0986e-7_dp ,&
                                 1.2210e-7_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -5.3530e-11_dp,  6.5280e-11_dp,&
                                -1.6120e-11_dp /)                       !3rd degree   

               elseif (which_total_pressure.eq.7.5_dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.38088_dp, 6.46146_dp, 106.313_dp/)
             
                  !Polynomial coefficients
                  b(1:3,1) = (/  5.1790e-1_dp ,  2.4775e-1_dp ,&
                                 2.8137e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/ -2.6324e-4_dp ,  3.3442e-4_dp ,& 
                                -3.1509e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  2.1708e-7_dp , -3.2412e-7_dp ,&  
                                 1.2210e-7_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -5.3530e-11_dp,  6.5280e-11_dp,& 
                                -1.6120e-11_dp /)                       !3rd degree   
         
               elseif (which_total_pressure.eq.10._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.38127_dp, 6.42086_dp,&
                                    106.1716_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  5.2958e-1_dp ,  2.7559e-1_dp ,&
                                 2.8205e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/ -2.7863e-4_dp ,  3.1132e-4_dp ,&
                                -3.1542e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  2.2206e-7_dp , -3.1913e-7_dp ,&
                                 1.2210e-7_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -5.3530e-11_dp,  6.5280e-11_dp,&
                                -1.6120e-11_dp /)                       !3rd degree   

               elseif (which_total_pressure.eq.12.5_dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.3372_dp, 5.75196_dp,&
                                    105.2711_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  5.3511e-1_dp ,  2.9518e-1_dp ,&
                                 2.8441e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/ -3.2011e-4_dp ,  3.3023e-4_dp ,&
                                -3.1627e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  2.4052e-7_dp , -3.2917e-7_dp ,&
                                 1.2210e-7_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -5.3530e-11_dp,  6.5280e-11_dp,&
                                -1.6120e-11_dp /)                       !3rd degree   

               elseif (which_total_pressure.eq.15._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.35161_dp, 5.93231_dp,&
                                    106.3668_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  5.5112e-1_dp ,  2.9174e-1_dp ,&
                                 2.8327e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/ -3.2424e-4_dp ,  3.3059e-4_dp ,&
                                -3.1587e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  2.3895e-7_dp , -3.2907e-7_dp ,&
                                 1.2210e-7_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -5.3530e-11_dp,  6.5280e-11_dp,&
                                -1.6120e-11_dp /)                       !3rd degree   

               elseif (which_total_pressure.eq.20._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 0.28627_dp, 4.43817_dp,&
                                    98.55304_dp /)
               
                  !Polynomial coefficients
                  b(1:3,1) = (/  4.6612e-1_dp ,  3.7614e-1_dp ,&
                                 3.1282e-1_dp  /)                       !0th degre
                  b(1:3,2) = (/ -1.9094e-4_dp ,  2.1665e-4_dp ,&
                                -3.2598e-4_dp  /)                       !1st degree   
                  b(1:3,3) = (/  1.9243e-7_dp , -2.9508e-7_dp ,&
                                 1.2210e-7_dp  /)                       !2nd degree   
                  b(1:3,4) = (/ -5.3530e-11_dp,  6.5280e-11_dp,&
                                -1.6120e-11_dp /)                       !3rd degree   

               else                                                     !Stop and print error message
                  write(msg,'(a,f3.1,a,a)') &
                     'Total pressure ',which_total_pressure,&
                     ' not available for correlation ',trim(wsgg_spec)
                  call shutdown(msg)
               endif
               
            else                                                        !Stop and print error message
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)
            endif
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Johansson et al. (2010)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('johansson2010_4gg')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 2
            
            if (which_molar_ratio.eq.0.125_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.0408_dp, 0.4217_dp,&
                                 5.2010_dp, 122.48_dp /)

               !Polynomial coefficients
               b(1:4,1) = (/  0.2719_dp,  0.3677_dp,&
                              0.2324_dp,  0.1058_dp /)
               b(1:4,2) = (/  0.0896_dp, -0.1284_dp,&
                             -0.1214_dp, -0.0602_dp /)
               b(1:4,3) = (/ -0.0327_dp, -0.0030_dp,&
                              0.0170_dp,  0.0080_dp /)
            
            elseif (which_molar_ratio.eq.1._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.0668_dp, 0.6818_dp,&
                                 5.9261_dp, 86.014_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.2887_dp,  0.1330_dp,&
                              0.3610_dp,  0.1843_dp /)
               b(1:4,2) = (/ -0.0760_dp,  0.3030_dp,&
                             -0.2112_dp, -0.1545_dp /)
               b(1:4,3) = (/  0.0604_dp, -0.1426_dp,&
                              0.0332_dp,  0.0347_dp /)
            
            else                                                        !Stop and print error message
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)               
            endif
            
         case('johansson2010_3gg')
            !Total number of gray gases
            number_wsgg_gases = 3
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 2
            
            if (which_molar_ratio.eq.0.125_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:3) = (/  0.0992_dp, 2.6589_dp, 88.1078_dp /)

               !Polynomial coefficients
               b(1:3,1) = (/  0.4995_dp,  0.3418_dp,  0.1273_dp /)
               b(1:3,2) = (/ -0.0170_dp, -0.1701_dp, -0.0726_dp /)
               b(1:3,3) = (/ -0.0393_dp,  0.0196_dp,  0.0101_dp /)

            elseif (which_molar_ratio.eq.1._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:3) = (/  0.1281_dp, 2.4256_dp, 59.557_dp /)
               
               !Polynomial coefficients
               b(1:3,1) = (/  0.2621_dp,  0.3898_dp,  0.2603_dp /)
               b(1:3,2) = (/  0.1798_dp, -0.0232_dp, -0.2193_dp /)
               b(1:3,3) = (/ -0.0491_dp, -0.0523_dp,  0.0502_dp /)

            else                                                        !Stop and print error message
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)               
            endif
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Krishnamoorthy (2010a) - air/fuel conditions
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('krishnamoorthy2010a')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 1
            
            if (which_molar_ratio.eq.2._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.21051_dp, 1.33782_dp,&
                                 8.55495_dp, 99.75649_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  1.54129e-1_dp,  2.43637e-1_dp,&
                              2.84084e-1_dp,  8.57853e-2_dp /)
               b(1:4,2) = (/  1.07579e-4_dp, -3.09769e-5_dp,&
                             -1.13634e-4_dp, -3.43141e-5_dp  /)
            
            elseif (which_molar_ratio.eq.3._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.22606_dp, 1.42179_dp,&
                                 9.19411_dp, 99.99325_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  1.74045e-1_dp,  2.40128e-1_dp,&
                              2.98507e-1_dp,  7.08215e-2_dp /)
               b(1:4,2) = (/  9.87576e-5_dp, -3.08707e-5_dp,&
                             -1.19403e-4_dp, -2.83286e-5_dp /)
            
            else                                                        !Stop and print error message
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)               
            endif
            
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Krishnamoorthy (2010b) - oxy/fuel conditions
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('krishnamoorthy2010b')
            !Total number of gray gases
            number_wsgg_gases = 3
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 1
            
            if (which_molar_ratio.eq.1._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:3) = (/ 0.09227_dp, 1.56657_dp, 99.99909_dp /)
               
               !Polynomial coefficients
               b(1:3,1) = (/ 1.53735e-1_dp,  4.66878e-1_dp,&
                             1.77659e-1_dp /)
               b(1:3,2) = (/ 1.88981e-4_dp, -1.22246e-4_dp,&
                            -7.10637e-5_dp /)
            
            elseif (which_molar_ratio.eq.(1._dp/3._dp)) then
               !Pressure-absorption coefficients
               kappa_p(1:3) = (/ 0.09386_dp, 1.02754_dp, 100.00000_dp /)
               
               !Polynomial coefficients
               b(1:3,1) = (/  1.97795e-1_dp,  3.59000e-1_dp,&
                              2.61828e-1_dp /)
               b(1:3,2) = (/  8.18110e-5_dp, -7.08406e-5_dp,&
                             -1.04731e-4_dp /)
            
            elseif (which_molar_ratio.eq.(1._dp/9._dp)) then
               !Pressure-absorption coefficients
               kappa_p(1:3) = (/ 0.06288_dp, 1.02333_dp, 100.00000_dp /)
               
               !Polynomial coefficients
               b(1:3,1) = (/  2.63797e-1_dp,  2.81308e-1_dp,&
                              2.30675e-1_dp /)
               b(1:3,2) = (/  6.12094e-5_dp, -3.78803e-5_dp,& 
                             -9.22702e-5_dp /)
            
            else                                                        !Stop and print error message
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)               
            endif            
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Kangwanpongpan et al. (2012)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('kangwanpongpan2012_fixed')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 5
            
            if (which_molar_ratio.eq.0.125_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.049191_dp, 0.343891_dp,&
                                 3.710740_dp, 106.080000_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.4329_dp, -0.5151_dp,&
                              0.3461_dp,  0.1443_dp /)
               b(1:4,2) = (/ -0.7255_dp,  7.2071_dp,&
                             -1.1838_dp, -0.4216_dp /)
               b(1:4,3) = (/  1.1384_dp, -21.319_dp,&
                              3.0944_dp,  1.5071_dp /)
               b(1:4,4) = (/  0.2661_dp,  28.355_dp,&
                             -4.1157_dp, -2.7378_dp /)
               b(1:4,5) = (/ -1.1786_dp, -17.950_dp,&
                              2.5339_dp,  2.1386_dp /)
               b(1:4,6) = (/  0.4713_dp,  4.3778_dp,&
                             -0.5857_dp, -0.5994_dp /)
                  
            elseif (which_molar_ratio.eq.0.25_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.045079_dp, 0.386642_dp,&
                                 3.764160_dp, 96.034300_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.3607_dp, -0.4713_dp,&
                              0.2999_dp,  0.1755_dp /)
               b(1:4,2) = (/ -0.1849_dp,  6.3750_dp,&
                             -0.4354_dp, -0.4710_dp /)
               b(1:4,3) = (/ -1.0644_dp, -18.102_dp,&
                              0.8603_dp,  1.4103_dp /)
               b(1:4,4) = (/  4.0610_dp,  23.534_dp,&
                             -1.1531_dp, -2.4705_dp /)
               b(1:4,5) = (/ -3.9878_dp, -14.680_dp,&
                              0.6442_dp,  1.9323_dp /)
               b(1:4,6) = (/  1.2272_dp,  3.5389_dp,&
                             -0.1167_dp, -0.5459_dp /)

            elseif (which_molar_ratio.eq.0.5_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.049191_dp, 0.421502_dp,&
                                 3.852390_dp, 83.534000_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.3026_dp, -0.4235_dp,&
                              0.2614_dp,  0.1954_dp /)
               b(1:4,2) = (/  0.3705_dp,  5.4049_dp,&
                              0.1391_dp, -0.3374_dp /)
               b(1:4,3) = (/ -3.3957_dp, -14.682_dp,&
                             -0.5843_dp,  0.6829_dp /)
               b(1:4,4) = (/  7.9979_dp,  18.707_dp,&
                              0.5412_dp, -1.3044_dp /)
               b(1:4,5) = (/ -6.8847_dp, -11.509_dp,&
                             -0.3740_dp,  1.1421_dp /)
               b(1:4,6) = (/  2.0098_dp,  2.7377_dp,&
                              0.1325_dp, -0.3490_dp /)
                  
            elseif (which_molar_ratio.eq.0.75_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.050784_dp, 0.431878_dp,&
                                 3.908780_dp, 75.525500_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.2887_dp, -0.4066_dp,&
                              0.2529_dp,  0.1974_dp /)
               b(1:4,2) = (/  0.5335_dp,  4.9951_dp,&
                              0.2466_dp, -0.1456_dp /)
               b(1:4,3) = (/ -4.2296_dp, -13.339_dp,&
                             -0.6019_dp, -0.1062_dp /)
               b(1:4,4) = (/  9.4993_dp,  16.952_dp,&
                              0.3143_dp, -0.1181_dp /)
               b(1:4,5) = (/ -8.0354_dp, -10.410_dp,&
                             -0.1583_dp,  0.3554_dp /)
               b(1:4,6) = (/  2.3307_dp,  2.4659_dp,&
                              0.0746_dp, -0.1542_dp /)

            elseif (which_molar_ratio.eq.1._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.051446_dp, 0.436145_dp,&
                                 3.948270_dp, 69.781000_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.2853_dp, -0.3976_dp,&
                              0.2532_dp,  0.1941_dp /)
               b(1:4,2) = (/  0.5853_dp,  4.7799_dp,&
                              0.2151_dp,  0.0337_dp /)
               b(1:4,3) = (/ -4.5894_dp, -12.689_dp,&
                             -0.2515_dp, -0.7945_dp /)
               b(1:4,4) = (/  10.205_dp,  16.172_dp,&
                             -0.3589_dp,  0.8954_dp /)
               b(1:4,5) = (/ -8.6057_dp, -9.9506_dp,&
                              0.3293_dp, -0.3114_dp /)
               b(1:4,6) = (/  2.4960_dp,  2.3557_dp,&
                             -0.0498_dp,  0.0104_dp /)

            elseif (which_molar_ratio.eq.2._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.051832_dp, 0.440593_dp,&
                                 3.981020_dp, 56.081800_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.2872_dp, -0.3748_dp,&
                              0.2655_dp,  0.1750_dp /)
               b(1:4,2) = (/  0.5807_dp,  4.3960_dp,&
                             -0.0908_dp,  0.5481_dp /)
               b(1:4,3) = (/ -4.9092_dp, -11.671_dp,&
                              1.3035_dp, -2.6576_dp /)
               b(1:4,4) = (/  11.012_dp,  15.104_dp,&
                             -2.9227_dp,  3.5816_dp /)
               b(1:4,5) = (/ -9.3465_dp, -9.3820_dp,&
                              2.0978_dp, -2.0627_dp /)
               b(1:4,6) = (/  2.7294_dp,  2.2263_dp,&
                             -0.4950_dp,  0.4412_dp /)
         
            elseif (which_molar_ratio.eq.4._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.051687_dp, 0.444449_dp,&
                                 3.933740_dp, 44.750100_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.2911_dp, -0.3416_dp,&
                              0.2790_dp,  0.1471_dp /)
               b(1:4,2) = (/  0.5050_dp,  4.0487_dp,&
                             -0.4828_dp,  1.0690_dp /)
               b(1:4,3) = (/ -4.8039_dp, -10.799_dp,&
                              3.0103_dp, -4.4490_dp /)
               b(1:4,4) = (/  11.037_dp,  14.158_dp,&
                             -5.5730_dp,  6.1082_dp /)
               b(1:4,5) = (/ -9.4886_dp, -8.8559_dp,&
                              3.8801_dp, -3.6923_dp /)
               b(1:4,6) = (/  2.7969_dp,  2.1031_dp,&
                             -0.9391_dp,  0.8400_dp /)      
                                                                        
            else                                                        !Stop and print error message
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)               
            endif
            
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Sun et al. (2012)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('sun2012')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 3
            
            if (which_molar_ratio.eq.0.25_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.46095_dp, 2.87928_dp,&
                                 11.61799_dp, 98.9513_dp /)

               !Polynomial coefficients
               b(1:4,1) = (/  0.07191_dp,  0.01175_dp,&
                              0.10367_dp,  0.13361_dp /)
               b(1:4,2) = (/  0.70335_dp,  0.60358_dp,&
                             -0.01238_dp, -0.12734_dp /)
               b(1:4,3) = (/ -0.53734_dp, -0.92935_dp,&
                             -0.05355_dp,  0.11549_dp /)
               b(1:4,4) = (/  0.01492_dp,  0.41215_dp,&
                              0.00127_dp, -0.06533_dp /)
               
            elseif (which_molar_ratio.eq.0.5_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.17702_dp, 2.13774_dp,&
                                 10.00212_dp, 79.99998_dp /)

               !Polynomial coefficients
               b(1:4,1) = (/  0.35759_dp,  0.05724_dp,&
                              0.01599_dp,  0.18157_dp /)
               b(1:4,2) = (/ -0.06667_dp,  0.48880_dp,&
                              0.57736_dp, -0.21301_dp /)
               b(1:4,3) = (/  0.38491_dp, -0.49309_dp,&
                             -0.89519_dp,  0.15808_dp /)
               b(1:4,4) = (/ -0.30731_dp,  0.11829_dp,&
                              0.35941_dp, -0.06616_dp /)
               
            elseif (which_molar_ratio.eq.1._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.19469_dp, 2.42095_dp,&
                                 11.43826_dp, 99.97641_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.40156_dp, -0.01038_dp,&
                              0.00258_dp,  0.21823_dp /)
               b(1:4,2) = (/ -0.39783_dp,  0.74092_dp,&
                              0.73953_dp, -0.30511_dp /)
               b(1:4,3) = (/  0.47735_dp, -0.64727_dp,&
                             -1.13108_dp,  0.22613_dp /)
               b(1:4,4) = (/ -0.07228_dp,  0.10195_dp,&
                              0.46468_dp, -0.08484_dp /)
            
            elseif (which_molar_ratio.eq.2._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.12962_dp, 2.20030_dp,&
                                 9.99797_dp, 79.99998_dp /)

               !Polynomial coefficients
               b(1:4,1) = (/  0.28829_dp,  0.01281_dp,&
                             -0.00963_dp,  0.23181_dp /)
               b(1:4,2) = (/ -0.04761_dp,  0.64013_dp,&
                              0.79150_dp, -0.17721_dp /)
               b(1:4,3) = (/  0.42561_dp, -0.55587_dp,&
                             -1.03171_dp, -0.07489_dp /)
               b(1:4,4) = (/ -0.21707_dp,  0.11501_dp,&
                              0.34351_dp,  0.08105_dp /)

            elseif (which_molar_ratio.eq.6._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.62645_dp, 3.54451_dp,&
                                 12.08918_dp, 147.83565_dp /)

               !Polynomial coefficients
               b(1:4,1) = (/ 0.107034_dp, -0.06638_dp,&
                              0.06942_dp,  0.16753_dp /)
               b(1:4,2) = (/  0.35070_dp,  0.76871_dp,&
                              0.44534_dp,  0.08112_dp /)
               b(1:4,3) = (/ -0.25322_dp, -0.69694_dp,&
                             -0.64016_dp, -0.43438_dp /)
               b(1:4,4) = (/  0.06531_dp,  0.15164_dp,&
                              0.21836_dp,  0.23381_dp /)

            elseif (which_molar_ratio.eq.12._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.61931_dp, 3.63459_dp,&
                                 12.53632_dp, 93.95137_dp /)

               !Polynomial coefficients
               b(1:4,1) = (/  0.13802_dp, -0.07767_dp,&
                              0.06950_dp,  0.16094_dp /)
               b(1:4,2) = (/  0.30308_dp,  0.86072_dp,&
                              0.43418_dp,  0.07453_dp /)
               b(1:4,3) = (/ -0.15691_dp, -0.86728_dp,&
                             -0.60729_dp, -0.41931_dp /)
               b(1:4,4) = (/  0.00033_dp,  0.24841_dp,&
                              0.20469_dp,  0.21952_dp /)

            else
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)
            endif
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Ziemniczak et al. (2013)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('ziemniczak2013_3gg')
            !Total number of gray gases
            number_wsgg_gases = 3
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 3
         
            !Pressure-absorption coefficients
            kappa_p(1:3) = (/ 0.237_dp, 2.611_dp, 36.409_dp /)
            
            !Polynomial coefficients
            b(1,1:4) = (/  2.256e-1_dp,  1.629e-4_dp,&
                          -3.114e-8_dp, -1.699e-12_dp /)
            b(2,1:4) = (/  1.429e-1_dp,  3.509e-4_dp,&
                          -2.495e-7_dp,  4.184e-11_dp /)
            b(3,1:4) = (/  2.743e-1_dp, -9.669e-5_dp,&
                          -4.743e-8_dp,  1.858e-11_dp /)

         case('ziemniczak2013_4gg')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4
         
            !Pressure-absorption coefficients
            kappa_p(1:4) = (/ 0.183_dp, 1.533_dp, 9.339_dp, 100.340_dp/)
      
            !Polynomial coefficients
            b(1,1:5) = (/  5.305e-2_dp,  7.903e-4_dp, -8.783e-7_dp,&
                           4.393e-10_dp, -7.715e-14_dp /)
            b(2,1:5) = (/  1.296e-1_dp,  1.717e-4_dp,  1.027e-8_dp,&
                          -7.952e-11_dp,  1.923e-14_dp /)
            b(3,1:5) = (/  1.440e-1_dp,  2.648e-4_dp, -3.741e-7_dp,&
                           1.565e-10_dp, -2.237e-14_dp /)
            b(4,1:5) = (/  1.298e-1_dp, -1.203e-5_dp, -9.266e-8_dp,&
                           5.020e-11_dp, -7.629e-15_dp /)

         case('ziemniczak2013_5gg')
            !Total number of gray gases
            number_wsgg_gases = 5
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 5
         
            !Pressure-absorption coefficients
            kappa_p(1:5) = (/ 0.127_dp, 0.752_dp, 3.094_dp,&
                              14.828_dp, 124.977_dp /)

            !Polynomial coefficients
            b(1,1:6) = (/ -3.381e-1_dp,  2.714e-3_dp, -4.378e-6_dp,&
                           3.255e-9_dp, -1.109e-12_dp,  1.411e-16_dp /)
            b(2,1:6) = (/  1.479e-1_dp, -2.608e-4_dp,  8.186e-7_dp,& 
                          -7.516e-10_dp,  2.802e-13_dp, -3.771e-17_dp /)
            b(3,1:6) = (/  5.298e-2_dp,  5.607e-4_dp, -8.839e-7_dp,&
                           6.564e-10_dp, -2.406e-13_dp,  3.366e-17_dp /)
            b(4,1:6) = (/  1.051e-1_dp,  2.180e-4_dp, -2.962e-7_dp,&
                           1.047e-10_dp, -5.589e-15_dp, -1.992e-18_dp /)
            b(5,1:6) = (/  1.257e-1_dp, -9.273e-5_dp,  4.385e-8_dp,&
                          -4.233e-11_dp,  2.229e-14_dp, -3.816e-18_dp /)
            
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Dorigon et al. (2013)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('dorigon2013')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4

            if (which_molar_ratio.eq.1._dp) then                        !p_w/p_c = 1
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 1.873e-1_dp, 1.723e0_dp,&
                                 1.248e1_dp, 1.449e2_dp /)

               !Polynomial coefficients
               b(1:4,1) = (/  7.197e-2_dp ,  1.107e-1_dp ,&
                              2.091e-1_dp ,  7.092e-2_dp  /)            !0th degree
               b(1:4,2) = (/  8.724e-4_dp ,  3.397e-4_dp ,&
                             -6.423e-5_dp ,  6.586e-5_dp  /)            !1st degree
               b(1:4,3) = (/ -9.690e-7_dp , -2.467e-7_dp ,&
                             -3.200e-8_dp , -1.278e-7_dp  /)            !2nd degree
               b(1:4,4) = (/  4.651e-10_dp,  4.647e-11_dp,&
                              1.718e-11_dp,  5.577e-11_dp /)            !3rd degree
               b(1:4,5) = (/ -7.917e-14_dp, -1.039e-15_dp,&
                             -2.105e-15_dp, -7.709e-15_dp /)            !4th degree

            elseif (which_molar_ratio.eq.2._dp) then                    !p_w/p_c = 2
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 1.921e-1_dp, 1.719e0_dp,&
                                 1.137e1_dp, 1.110e2_dp /)

               !Polynomial coefficients
               b(1:4,1) = (/  5.617e-2_dp ,  1.426e-1_dp ,&
                              1.362e-1_dp ,  1.222e-1_dp  /)            !0th degree
               b(1:4,2) = (/  7.844e-4_dp ,  1.795e-4_dp ,&
                              2.574e-4_dp , -2.327e-5_dp  /)            !1st degree
               b(1:4,3) = (/ -8.563e-7_dp , -1.077e-8_dp ,&
                             -3.711e-7_dp , -7.492e-8_dp  /)            !2nd degree
               b(1:4,4) = (/  4.246e-10_dp, -6.971e-11_dp,&
                              1.575e-10_dp,  4.275e-11_dp /)            !3rd degree
               b(1:4,5) = (/ -7.440e-14_dp,  1.774e-14_dp,&
                             -2.267e-14_dp, -6.608e-15_dp /)            !4th degree
               
            else                                                        !Stop and print error message
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)
            endif
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Krishnamoorthy (2013)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case ('krishnamoorthy2013')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 1
            
            if (which_molar_ratio.eq.0.11_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.06592_dp, 0.99698_dp,&
                                 10.00038_dp, 100.00000_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  2.39641e-1_dp,  3.42342e-1_dp,&
                              1.37773e-1_dp, 4.40724e-2_dp /)
               b(1:4,2) = (/  7.85445e-5_dp, -9.47416e-5_dp,&
                             -5.51091e-5_dp, 7.26634e-6_dp /)

            elseif (which_molar_ratio.eq.0.5_dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.10411_dp, 1.00018_dp,&
                                 9.99994_dp, 100.00000_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  1.89029e-1_dp,  2.87021e-1_dp,&
                              2.54516e-1_dp,  6.54289e-2_dp /)
               b(1:4,2) = (/  9.33340e-5_dp, -5.32833e-5_dp,&
                             -1.01806e-4_dp, -2.25973e-5_dp /)

            elseif (which_molar_ratio.eq.1._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.20616_dp, 1.39587_dp,&
                                 8.56904_dp, 99.75698_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  1.91464e-1_dp,  2.34876e-1_dp,&
                              2.47320e-1_dp,  9.59426e-2_dp /)
               b(1:4,2) = (/  9.22363e-5_dp, -4.25444e-5_dp,&
                             -9.89282e-5_dp, -3.83770e-5_dp /)

            elseif (which_molar_ratio.eq.2._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.21051_dp, 1.33782_dp,&
                                 8.55495_dp, 99.75649_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  1.54129e-1_dp,  2.43637e-1_dp,&
                              2.84084e-1_dp,  8.57853e-2_dp /)
               b(1:4,2) = (/  1.07579e-4_dp, -3.09769e-5_dp,&
                             -1.13634e-4_dp, -3.43141e-5_dp /)

            else
               write(msg,'(a,f3.1,a,a)') &
                  'Partial pressure ratio ',which_molar_ratio,&
                  ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)
            endif
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Yin (2013)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('yin2013')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 3
            
            select case(trim(species_spec_aux))
               case('co2')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.163233_dp, 13.096584_dp,&
                                    175.474735_dp, 1310.847307_dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/  0.204623_dp, -0.020227_dp,&
                                 0.044221_dp,  0.039311_dp /)
                  b(1:4,2) = (/ -0.378060_dp,  0.256006_dp,&
                                 0.003850_dp, -0.054832_dp /)
                  b(1:4,3) = (/  0.666639_dp, -0.195201_dp,&
                                -0.020175_dp,  0.025370_dp /)
                  b(1:4,4) = (/ -0.203453_dp,  0.040493_dp,&
                                 0.004919_dp, -0.003891_dp /)

               case('h2o')
                  if (which_h2o_pressure.eq.0._dp) then
                     !Pressure-absorption coefficients
                     kappa_p(1:4) = (/ 0.085523_dp, 0.475777_dp,&
                                       8.549733_dp, 201.906503_dp /)

                     !Polynomial coefficients
                     b(1:4,1) = (/  0.966357_dp,  0.662059_dp,&
                                    0.060870_dp,  0.103568_dp /)
                     b(1:4,2) = (/ -0.790165_dp, -2.262877_dp,& 
                                    0.436788_dp, -0.153135_dp /)
                     b(1:4,3) = (/ -0.050144_dp,  2.309473_dp,&
                                   -0.395493_dp,  0.074910_dp /)
                     b(1:4,4) = (/  0.115202_dp, -0.572895_dp,&
                                    0.085146_dp, -0.012091_dp /)

                  elseif (which_h2o_pressure.eq.0.05_dp) then
                     !Pressure-absorption coefficients
                     kappa_p(1:4) = (/ 0.232724_dp, 2.134299_dp,&
                                       9.266065_dp, 134.988332_dp /)
                     
                     !Polynomial coefficients
                     b(1:4,1) = (/  0.340618_dp,  0.175818_dp,&
                                    0.044325_dp,  0.126628_dp /)
                     b(1:4,2) = (/ -0.105469_dp, -0.063466_dp,&
                                    0.288376_dp, -0.186480_dp /)
                     b(1:4,3) = (/  0.068051_dp,  0.086631_dp,&
                                   -0.258205_dp,  0.090755_dp /)
                     b(1:4,4) = (/ -0.017828_dp, -0.026581_dp,&
                                    0.054333_dp, -0.014569_dp /)

                  elseif (which_h2o_pressure.eq.1._dp) then
                     !Pressure-absorption coefficients
                     kappa_p(1:4) = (/ 0.065411_dp, 0.696552_dp,&
                                       4.862610_dp, 60.255980_dp /)

                     !Polynomial coefficients
                     b(1:4,1) = (/ -0.077336_dp,  0.506777_dp,&
                                   -0.079989_dp,  0.373898_dp /)
                     b(1:4,2) = (/  0.661776_dp, -0.758948_dp,&
                                    0.851078_dp, -0.540887_dp /)
                     b(1:4,3) = (/ -0.362515_dp,  0.516146_dp,&
                                   -0.604264_dp,  0.258923_dp /)
                     b(1:4,4) = (/  0.053534_dp, -0.102909_dp,&
                                    0.113500_dp, -0.040957_dp /)
                     
                  else
                     write(msg,'(a,f3.1,a,a)') 'H2O partial pressure ',&
                        which_h2o_pressure,' atm not available for &
                        &correlation ',trim(wsgg_spec)
                     call shutdown(msg)
                  endif
                  
               case('mixture')           
                  if (which_molar_ratio.eq.0.05_dp) then                !p_w/p_c = 0.05
                     !Pressure-absorption coefficients
                     kappa_p(1:4) = (/ 0.352505_dp, 8.210621_dp,&
                                       137.410012_dp, 1269.710976_dp /)

                     !Polynomial coefficients
                     b(1:4,1) = (/  0.315106_dp,  0.092474_dp,&
                                    0.031702_dp,  0.046138_dp /)        !0th degree
                     b(1:4,2) = (/  0.023475_dp,  0.109146_dp,&
                                    0.037396_dp, -0.061392_dp /)        !1st degree
                     b(1:4,3) = (/ -0.057930_dp, -0.121000_dp,&
                                   -0.040731_dp,  0.027164_dp /)        !2nd degree
                     b(1:4,4) = (/  0.008408_dp,  0.027145_dp,&
                                    0.008742_dp, -0.003996_dp /)        !3rd degree
            
                  elseif (which_molar_ratio.eq.1._dp) then              !p_w/p_c = 1
                     !Pressure-absorption coefficients
                     kappa_p(1:4) = (/ 0.261021_dp, 3.147817_dp,&
                                       54.265868_dp, 482.900353_dp /)

                     !Polynomial coefficients
                     b(1:4,1) = (/  0.500119_dp,  0.071592_dp,&
                                    0.155320_dp,  0.072615_dp /)        !0th degree
                     b(1:4,2) = (/ -0.447068_dp,  0.508252_dp,&
                                   -0.104294_dp, -0.100601_dp /)        !1st degree
                     b(1:4,3) = (/  0.286878_dp, -0.384253_dp,&
                                    0.014096_dp,  0.046681_dp /)        !2nd degree
                     b(1:4,4) = (/ -0.059165_dp,  0.073477_dp,&
                                    0.001643_dp, -0.007224_dp /)        !3rd degree

                  elseif (which_molar_ratio.eq.2._dp) then              !p_w/p_c = 2
                     !Pressure-absorption coefficients
                     kappa_p(1:4) = (/ 0.179160_dp, 2.388971_dp,&
                                       28.415805_dp, 253.059089_dp /)
      
                     !Polynomial coefficients
                     b(1:4,1) = (/  0.542458_dp,  0.101734_dp,&
                                    0.146066_dp,  0.129511_dp /)        !0th degree
                     b(1:4,2) = (/ -0.658411_dp,  0.518429_dp,&
                                   -0.008745_dp, -0.187993_dp /)        !1st degree
                     b(1:4,3) = (/  0.466444_dp, -0.386151_dp,&
                                   -0.058325_dp,  0.090709_dp /)        !2nd degree
                     b(1:4,4) = (/ -0.100186_dp,  0.073453_dp,&
                                    0.015984_dp, -0.014493_dp /)        !3rd degree

                  else                                                  !Stop and print error message
                     write(msg,'(a,f3.1,a,a)') &
                        'Partial pressure ratio ',which_molar_ratio,&
                        ' not available for correlation ',&
                        trim(wsgg_spec)
                     call shutdown(msg)
                  endif
            endselect
            
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Cassol et al. (2014)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('cassol2014')
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4
            
            select case(trim(species_spec_aux))
               case('co2')
                  !Number of gray gases
                  number_wsgg_gases = 4
                  
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.138_dp, 1.895_dp,&
                                    13.301_dp, 340.811_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/   0.09990_dp ,   0.00942_dp ,&
                                 0.14511_dp ,  -0.02915_dp  /)          !0th degree
                  b(1:4,2) = (/  64.41e-5_dp ,  10.36e-5_dp ,&
                                -30.73e-5_dp ,  25.23e-5_dp  /)         !1st degree
                  b(1:4,3) = (/ -89.94e-8_dp , -2.277e-8_dp ,&  
                                 37.65e-8_dp , -26.10e-8_dp  /)         !2nd degree
                  b(1:4,4) = (/  41.27e-11_dp, -2.134e-11_dp,& 
                                -18.41e-11_dp,  9.965e-11_dp /)         !3rd degree
                  b(1:4,5) = (/ -67.74e-15_dp,  6.497e-15_dp,&  
                                 30.16e-15_dp, -13.26e-15_dp /)         !4th degree
         
               case('h2o')
                  !Number of gray gases
                  number_wsgg_gases = 4
                  
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.171_dp, 1.551_dp,&
                                    5.562_dp, 49.159_dp /)
         
                  !Polynomial coefficients
                  b(1:4,1) = (/   0.06617_dp ,   0.11045_dp ,&
                                 -0.04915_dp ,   0.23675_dp  /)         !0th degree
                  b(1:4,2) = (/  55.48e-5_dp ,   0.576e-5_dp,&  
                                 70.63e-5_dp , -18.91e-5_dp  /)         !1st degree
                  b(1:4,3) = (/ -48.41e-8_dp ,  24.00e-8_dp ,&
                                -70.12e-8_dp , -0.907e-8_dp  /)         !2nd degree
                  b(1:4,4) = (/  22.27e-11_dp, -17.01e-11_dp,&
                                 26.07e-11_dp,  4.082e-11_dp /)         !3rd degree
                  b(1:4,5) = (/ -40.17e-15_dp,  30.96e-15_dp,&
                                -34.94e-15_dp, -8.778e-15_dp /)         !4th degree
         
               case('soot')
                  select case(number_wsgg_gases_soot)
                     case(2)
                        !Number of gray gases
                        number_wsgg_gases = 2
                  
                        !Pressure-absorption coefficients
                        kappa_p(1:2) = (/ 22313.49_dp, 466624.8_dp /)

                        !Polynomial coefficients
                        b(1:2,1) = (/  0.95552_dp  ,  0.08010_dp   /)   !0th degree
                        b(1:2,2) = (/ -1.431e-3_dp ,  1.290e-3_dp  /)   !1st degree
                        b(1:2,3) = (/  9.871e-7_dp , -7.874e-7_dp  /)   !2nd degree
                        b(1:2,4) = (/ -3.390e-10_dp,  2.322e-10_dp /)   !3rd degree
                        b(1:2,5) = (/  4.555e-14_dp, -3.084e-14_dp /)   !4th degree
                  
                     case(3)
                        !Number of gray gases
                        number_wsgg_gases = 3
                     
                        !Pressure-absorption coefficients
                        kappa_p(1:3) = (/ 1251.56_dp, 50470.6_dp,&
                                          460361.0_dp /)
         
                        !Polynomial coefficients
                        b(1:3,1) = (/  0.02812_dp   ,   1.25626_dp ,&
                                      -0.25179_dp   /)                  !0th degree
                        b(1:3,2) = (/ -1.271e-4_dp  , -18.74e-4_dp ,&  
                                       18.55e-4_dp  /)                  !1st degree
                        b(1:3,3) = (/  1.395e-7_dp  ,  12.09e-7_dp ,& 
                                      -11.39e-7_dp  /)                  !2nd degree
                        b(1:3,4) = (/ -5.672e-11_dp , -39.89e-11_dp,&  
                                      34.46e-11_dp /)                   !3rd degree
                        b(1:3,5) = (/  8.0868e-15_dp,  52.85e-15_dp,& 
                                     -45.58e-15_dp /)                   !4th degree
                  
                     case(4)
                        !Number of gray gases
                        number_wsgg_gases = 4
                     
                        !Pressure-absorption coefficients
                        kappa_p(1:4) = (/ 2857.86_dp, 39234.9_dp,&
                                          160748.0_dp, 495898.0_dp /)
   
                        !Polynomial coefficients
                        b(1:4,1) = (/  0.00129_dp   ,    1.26110_dp,&
                                         -0.25757_dp,   0.07980_dp  /)  !0th degree
                        b(1:4,2) = (/ -0.545e-5_dp  , -319.2e-5_dp ,&  
                                        362.1e-5_dp , -72.08e-5_dp  /)  !1st degree
                        b(1:4,3) = (/  0.123e-7_dp  ,   27.72e-7_dp,&  
                                        -40.12e-7_dp,  15.87e-7_dp  /)  !2nd degree
                        b(1:4,4) = (/ -0.847e-11_dp , -100.5e-11_dp,&  
                                        154.9e-11_dp, -70.89e-11_dp /)  !3rd degree
                        b(1:4,5) = (/  1.6807e-15_dp,  132.8e-15_dp,& 
                                       -207.8e-15_dp,  97.69e-15_dp /)  !4th degree
                                             
                     case default
                        !Number of gray gases
                        number_wsgg_gases = 0
                     
                        !Pressure-absorption coefficients
                        kappa_p(:) = 0._dp
                     
                        !Polynomial coefficients
                        b(:,:) = 0._dp
                     
                  endselect
                     
                  !Correct the number of gray gases if transparent
                  !windows for soot are to be neglected
                  if (.not.soot_transparent_windows) &
                     number_wsgg_gases = number_wsgg_gases - 1
                     
            end select
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Coelho & França (2017)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('coelho2017')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4
            
            if (trim(species_spec_aux).eq.'co2') then
               if (which_total_pressure.eq.1._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.261_dp, 3.441_dp,&
                                    25.691_dp, 309.289_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/  2.442e-2_dp ,  6.953e-2_dp ,&  
                                 1.102e-1_dp , -4.009e-2_dp  /)
                  b(1:4,2) = (/  6.183e-4_dp , -1.339e-4_dp ,& 
                                -2.044e-4_dp ,  2.737e-4_dp  /)
                  b(1:4,3) = (/ -7.486e-7_dp ,  2.582e-7_dp ,&  
                                 2.296e-7_dp , -2.829e-7_dp  /)
                  b(1:4,4) = (/  3.335e-10_dp, -1.481e-10_dp,& 
                                -1.105e-10_dp,  1.097e-10_dp /)
                  b(1:4,5) = (/ -5.198e-14_dp,  2.599e-14_dp,&  
                                 1.817e-14_dp, -1.491e-14_dp /)

               elseif (which_total_pressure.eq.10._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.045_dp, 0.374_dp,&
                                    5.217_dp, 151.891_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/  7.567e-1_dp ,  2.455e-2_dp ,&
                                 1.070e-1_dp ,  5.294e-2_dp  /)
                  b(1:4,2) = (/ -1.931e-3_dp ,  8.742e-4_dp ,&
                                -1.425e-4_dp ,  9.748e-5_dp  /)
                  b(1:4,3) = (/  2.233e-6_dp , -1.200e-6_dp ,&
                                 2.505e-7_dp , -1.124e-7_dp  /)
                  b(1:4,4) = (/ -1.030e-9_dp ,  5.665e-10_dp,&
                                -1.435e-10_dp,  3.459e-11_dp /)
                  b(1:4,5) = (/  1.632e-13_dp, -9.037e-14_dp,&
                                 2.525e-14_dp, -3.061e-15_dp /)

               elseif (which_total_pressure.eq.40._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.011_dp, 0.175_dp,&
                                    2.133_dp, 58.250_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/  5.185e-2_dp ,  5.762e-1_dp ,&
                                 1.122e-1_dp ,  1.244e-1_dp  /)
                  b(1:4,2) = (/ -1.082e-4_dp , -1.065e-3_dp ,&
                                 4.491e-4_dp , -4.333e-5_dp  /)
                  b(1:4,3) = (/  4.682e-7_dp ,  1.114e-6_dp ,&
                                -6.088e-7_dp ,  3.704e-8_dp  /)
                  b(1:4,4) = (/ -2.774e-10_dp, -5.083e-10_dp,&
                                 2.673e-10_dp, -3.141e-11_dp /)
                  b(1:4,5) = (/  4.658e-14_dp,  8.096e-14_dp,&
                                -3.985e-14_dp,  6.992e-15_dp /)
               
               else
                  write(msg,'(a,f3.1,a,a)') &
                     'Total pressure ',which_total_pressure,&
                     ' not available for correlation ',trim(wsgg_spec)
                  call shutdown(msg)
               endif
               
            elseif (trim(species_spec_aux).eq.'h2o') then
               if (which_total_pressure.eq.1._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.182_dp, 1.358_dp,&
                                    7.819_dp, 75.164_dp /)

                  !Polynomial coefficients  
                  b(1:4,1) = (/  8.998e-2_dp ,  6.831e-2_dp ,&
                                 3.448e-4_dp ,  2.203e-1_dp  /)
                  b(1:4,2) = (/  4.567e-4_dp ,  2.308e-4_dp ,&
                                 6.119e-4_dp , -3.131e-4_dp  /)
                  b(1:4,3) = (/ -3.369e-7_dp , -2.674e-9_dp ,&
                                -6.886e-7_dp ,  1.741e-7_dp  /)
                  b(1:4,4) = (/  1.510e-10_dp, -7.251e-11_dp,&
                                 2.789e-10_dp, -4.492e-11_dp /)
                  b(1:4,5) = (/ -2.723e-14_dp,  1.719e-14_dp,&
                                -3.980e-14_dp,  4.552e-15_dp /)

               elseif (which_total_pressure.eq.10._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.036_dp, 0.319_dp,&
                                    2.808_dp, 26.235_dp /)
                  
                  !Polynomial coefficients  
                  b(1:4,1) = (/  6.747e-1_dp , -2.711e-2_dp ,&
                                -7.584e-2_dp ,  4.024e-1_dp  /)
                  b(1:4,2) = (/ -1.591e-3_dp ,  8.393e-4_dp ,&
                                 5.586e-4_dp , -8.570e-6_dp  /)
                  b(1:4,3) = (/  1.750e-6_dp , -9.044e-7_dp ,&
                                -5.761e-8_dp , -4.818e-7_dp  /)
                  b(1:4,4) = (/ -7.523e-10_dp,  4.511e-10_dp,&
                                -1.575e-10_dp,  3.003e-10_dp /)
                  b(1:4,5) = (/  1.126e-13_dp, -7.964e-14_dp,&
                                 4.226e-14_dp, -5.284e-14_dp /)

               elseif (which_total_pressure.eq.40._dp) then
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.014_dp, 0.132_dp,&
                                    1.018_dp, 12.172_dp /)
                  
                  !Polynomial coefficients  
                  b(1:4,1) = (/  4.350e-1_dp ,  4.570e-1_dp ,&
                                -2.736e-1_dp ,  4.533e-1_dp  /)
                  b(1:4,2) = (/ -1.108e-3_dp , -5.952e-4_dp ,&
                                 9.947e-4_dp ,  3.093e-4_dp  /)
                  b(1:4,3) = (/  1.354e-6_dp ,  4.254e-7_dp ,&
                                -4.444e-7_dp , -7.827e-7_dp  /)
                  b(1:4,4) = (/ -6.439e-10_dp, -6.605e-11_dp,&
                                 3.448e-11_dp,  3.897e-10_dp /)
                  b(1:4,5) = (/  1.057e-13_dp, -6.940e-15_dp,&
                                 7.893e-15_dp, -6.140e-14_dp /)

               else
                  write(msg,'(a,f3.1,a,a)') &
                     'Total pressure ',which_total_pressure,&
                     ' not available for correlation ',trim(wsgg_spec)
                  call shutdown(msg)
               endif
               
            else
               write(msg,'(a,a,a,a)') &
                     'Species ',trim(species_spec_aux),&
                     ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)
            endif
            
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Centeno et al. (2018)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('centeno2018')                                            !p_w/p_c = 1.5
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4

            !Pressure-absorption coefficients
            kappa_p(1:4) = (/ 0.178_dp, 1.501_dp,&
                              9.510_dp, 110.348_dp /)

            !Polynomial coefficients
            b(1:4,1) = (/  6.345e-2_dp ,  1.042e-1_dp ,&
                           1.864e-1_dp ,  1.071e-1_dp  /)               !0th degree
            b(1:4,2) = (/  8.077e-4_dp ,  3.102e-4_dp ,&
                           4.624e-5_dp ,  2.209e-5_dp  /)               !1st degree
            b(1:4,3) = (/ -8.900e-7_dp , -1.904e-7_dp ,&
                          -1.318e-7_dp , -1.056e-7_dp  /)               !2nd degree
            b(1:4,4) = (/  4.337e-10_dp,  1.974e-11_dp,&
                           5.402e-11_dp,  4.986e-11_dp /)               !3rd degree
            b(1:4,5) = (/ -7.472e-14_dp,  3.091e-15_dp,&
                          -7.232e-15_dp, -7.007e-15_dp /)               !4th degree

         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Marzouk (2018)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('marzouk2018')
            !Total number of gray gases
            number_wsgg_gases = 4
         
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4
            
            if (trim(species_spec_aux).eq.'co2') then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.04089024_dp, 0.48865123_dp,&
                                 4.91049904_dp, 108.71403783_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.23633250_dp, -0.03365397_dp,&
                              0.10044540_dp,  0.03076686_dp /)
               b(1:4,2) = (/ -0.14137986_dp,  0.61079456_dp,&
                             -0.19075064_dp,  0.23689082_dp /)
               b(1:4,3) = (/ -0.09615643_dp, -0.98169775_dp,&
                              0.42907206_dp, -0.41816812_dp /)
               b(1:4,4) = (/  0.11093274_dp,  0.58100617_dp,&
                             -0.34556923_dp,  0.24019457_dp /)
               b(1:4,5) = (/ -0.02849616_dp, -0.11917163_dp,&
                              0.08630678_dp, -0.04679157_dp /)
            
            elseif (trim(species_spec_aux).eq.'h2o') then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.12995182_dp,  1.40673720_dp,&
                                 6.43286122_dp, 52.43600070_dp /)
               
               !Polynomial coefficients
               b(1:4,1) = (/  0.31147381_dp, -0.06295794_dp,&
                              0.06417476_dp,  0.34675301_dp /)
               b(1:4,2) = (/ -0.38154730_dp,  1.33734280_dp,&
                              0.50816204_dp, -0.78767144_dp /)
               b(1:4,3) = (/  1.19826122_dp, -1.96567794_dp,&
                             -0.63233721_dp,  0.66356908_dp /)
               b(1:4,4) = (/ -0.97093067_dp,  1.19703078_dp,&
                              0.20080577_dp, -0.24173161_dp /)
               b(1:4,5) = (/  0.24391581_dp, -0.26885101_dp,&
                             -0.00714479_dp,  0.03158120_dp /)

            else
               if (which_molar_ratio.eq.1._dp/20._dp) then    
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.03404788_dp,  0.40856120_dp,&
                                    4.98834612_dp, 110.10062265_dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/  0.34980614_dp,  0.13122203_dp,&
                                 0.22973752_dp,  0.04960046_dp /)
                  b(1:4,2) = (/  0.13570199_dp,  0.82233394_dp,&
                                -0.41343100_dp,  0.15949517_dp /)
                  b(1:4,3) = (/ -0.24465194_dp, -1.58968389_dp,&
                                 0.53075920_dp, -0.31369820_dp /)
                  b(1:4,4) = (/  0.12489443_dp,  0.98158637_dp,&
                                -0.33611838_dp,  0.18146967_dp /)
                  b(1:4,5) = (/ -0.02973515_dp, -0.20425005_dp,&
                                 0.07556505_dp, -0.03495927_dp /)
                  
               elseif (which_molar_ratio.eq.1._dp/8._dp) then    
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.04827306_dp,  0.57059748_dp,&
                                    5.88731430_dp, 112.02187600_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/  0.31104289_dp,  0.08873874_dp,&
                                 0.27104897_dp,  0.06689982_dp /)
                  b(1:4,2) = (/  0.16957294_dp,  0.98442535_dp,&
                                -0.43372133_dp,  0.09879398_dp /)
                  b(1:4,3) = (/ -0.25694276_dp, -1.65763972_dp,&
                                 0.44899506_dp, -0.24474902_dp /)
                  b(1:4,4) = (/  0.17010955_dp,  0.95399708_dp,&
                                -0.25495556_dp,  0.14826314_dp /)
                  b(1:4,5) = (/ -0.04841856_dp, -0.18988730_dp,&
                                 0.05489626_dp, -0.02910387_dp /)
                  
               elseif (which_molar_ratio.eq.1._dp/4._dp) then    
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.06181095_dp,  0.74689552_dp,&
                                    6.94567399_dp, 111.19582630_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/  0.30477949_dp,  0.05129759_dp,&
                                 0.29414139_dp,  0.08750500_dp /)
                  b(1:4,2) = (/  0.11241736_dp,  1.09457041_dp,&
                                -0.40524068_dp,  0.03091227_dp /)
                  b(1:4,3) = (/ -0.17679091_dp, -1.65423859_dp,&
                                 0.31944295_dp, -0.16985663_dp /)
                  b(1:4,4) = (/  0.16127789_dp,  0.88986523_dp,&
                                -0.15282118_dp,  0.11268674_dp /)
                  b(1:4,5) = (/ -0.05456268_dp, -0.16888241_dp,&
                                 0.03075060_dp, -0.02286324_dp /)
                  
               elseif (which_molar_ratio.eq.1._dp/2._dp) then    
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.07392381_dp,  0.87830788_dp,&
                                    7.62031942_dp, 103.11882353_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/  0.32418038_dp,  0.00555671_dp,&
                                 0.30468133_dp,  0.12105722_dp /)
                  b(1:4,2) = (/ -0.07613144_dp,  1.20372256_dp,&
                                -0.29414360_dp, -0.06682114_dp /)
                  b(1:4,3) = (/  0.06947528_dp, -1.62548826_dp,&
                                 0.08133544_dp, -0.06814394_dp /)
                  b(1:4,4) = (/  0.06388190_dp,  0.80603343_dp,&
                                 0.00596004_dp,  0.06582474_dp /)
                  b(1:4,5) = (/ -0.04344751_dp, -0.14345971_dp,&
                                -0.00407591_dp, -0.01475981_dp /)
                  
               elseif (which_molar_ratio.eq.1._dp) then    
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.08130729_dp,  0.91245096_dp,&
                                    7.57172165_dp, 88.53048852_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/  0.35860798_dp, -0.02738398_dp,&
                                 0.28289630_dp,  0.16995541_dp /)
                  b(1:4,2) = (/ -0.32518384_dp,  1.21396014_dp,&
                                -0.04872198_dp, -0.19251673_dp /)
                  b(1:4,3) = (/  0.37665224_dp, -1.45788036_dp,&
                                -0.30083097_dp,  0.04850927_dp /)
                  b(1:4,4) = (/ -0.06634910_dp,  0.65049879_dp,&
                                 0.22784510_dp,  0.01752049_dp /)
                  b(1:4,5) = (/ -0.02526599_dp, -0.10502236_dp,&
                                -0.04906112_dp, -0.00716657_dp /)
                  
               elseif (which_molar_ratio.eq.2._dp) then    
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.08516379_dp,  0.92208873_dp,&
                                    7.72024076_dp, 77.29996669_dp /)
                  
                  !Polynomial coefficients
                  b(1:4,1) = (/  0.38429579_dp, -0.01789618_dp,&
                                 0.22445124_dp,  0.21829730_dp /)
                  b(1:4,2) = (/ -0.51628961_dp,  1.09039550_dp,&
                                 0.27231656_dp, -0.32152661_dp /)
                  b(1:4,3) = (/  0.61663009_dp, -1.17058882_dp,&
                                -0.72218639_dp,  0.16134270_dp /)
                  b(1:4,4) = (/ -0.17595628_dp,  0.45398513_dp,&
                                 0.44900713_dp, -0.02406260_dp /)
                  b(1:4,5) = (/ -0.00790864_dp, -0.06194554_dp,&
                                -0.09093902_dp, -0.00160557_dp /)

               elseif (which_molar_ratio.eq.5._dp) then    
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.09039754_dp,  1.00421305_dp,&
                                    8.74756608_dp,  75.23553830_dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/  0.40181309_dp,  0.02902963_dp,&
                                 0.13754056_dp,  0.25784580_dp /)
                  b(1:4,2) = (/ -0.65335779_dp,  0.87186760_dp,&
                                 0.67684212_dp, -0.49478680_dp /)
                  b(1:4,3) = (/  0.84111554_dp, -0.82065876_dp,&
                                -1.23645727_dp,  0.35301430_dp /)
                  b(1:4,4) = (/ -0.31657749_dp,  0.25208835_dp,&
                                 0.70940728_dp, -0.10826644_dp /)
                  b(1:4,5) = (/  0.02285388_dp, -0.02258414_dp,&
                                -0.13804737_dp,  0.01158972_dp /)

               elseif (which_molar_ratio.eq.20._dp) then    
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 0.10276016_dp,  1.02692687_dp,&
                                    4.87450005_dp, 40.32952963_dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/  0.24440974_dp,  0.04520324_dp,&
                                -0.06133110_dp,  0.41485758_dp /)
                  b(1:4,2) = (/ -0.09020790_dp,  0.87741307_dp,&
                                 1.14799062_dp, -0.79332082_dp /)
                  b(1:4,3) = (/  0.29974144_dp, -1.36498366_dp,&
                                -1.33642854_dp,  0.52692383_dp /)
                  b(1:4,4) = (/ -0.12697120_dp,  0.87138940_dp,&
                                 0.51146296_dp, -0.12895066_dp /)
                  b(1:4,5) = (/  0.00335368_dp, -0.20292017_dp,&
                                -0.05776215_dp,  0.00563381_dp /)
               
               else
                  write(msg,'(a,f3.1,a,a)') &
                     'Partial pressure ratio ',which_molar_ratio,&
                     ' not available for correlation ',trim(wsgg_spec)
                  call shutdown(msg)
               endif   
            endif

         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Coelho & França (2019)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('coelho2019')
            !Total number of gray gases
            number_wsgg_gases = 4
         
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4
            
            if (which_total_pressure.eq.1._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.180_dp, 1.517_dp,&
                                 9.423_dp, 101.324_dp /)
   
               !Polynomial coefficients
               b(1:4,1) = (/  5.047e-2_dp ,  1.152e-1_dp ,&
                              1.763e-1_dp ,  1.251e-1_dp  /)            !0th degree
               b(1:4,2) = (/  7.887e-4_dp ,  2.581e-4_dp ,&  
                              1.266e-4_dp , -1.206e-5_dp  /)            !1st degree
               b(1:4,3) = (/ -8.564e-7_dp , -1.106e-7_dp ,& 
                             -2.262e-7_dp , -8.320e-8_dp  /)            !2nd degree
               b(1:4,4) = (/  4.203e-10_dp, -2.021e-11_dp,&  
                              9.513e-11_dp,  4.343e-11_dp /)            !3rd degree
               b(1:4,5) = (/ -7.299e-14_dp,  9.604e-15_dp,& 
                             -1.345e-14_dp, -6.309e-15_dp /)            !4th degree

            elseif (which_total_pressure.eq.2._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.109_dp, 0.852_dp,&
                                 5.533_dp, 68.011_dp /)
         
               !Polynomial coefficients
               b(1:4,1) = (/  1.847e-1_dp , -5.099e-2_dp ,&
                              2.225e-1_dp ,  2.147e-1_dp  /)            !0th degree
               b(1:4,2) = (/  3.984e-4_dp ,  6.629e-4_dp ,&
                              1.740e-4_dp , -7.606e-5_dp  /)            !1st degree
               b(1:4,3) = (/ -5.816e-7_dp , -4.744e-7_dp ,&
                             -2.281e-7_dp , -9.287e-8_dp  /)            !2nd degree
               b(1:4,4) = (/  3.585e-10_dp,  1.327e-10_dp,&
                              7.211e-11_dp,  6.141e-11_dp /)            !3rd degree
               b(1:4,5) = (/ -7.058e-14_dp, -1.382e-14_dp,&
                             -7.524e-15_dp, -1.005e-14_dp /)            !4th degree

            elseif (which_total_pressure.eq.5._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.055_dp, 0.402_dp,&
                                 3.399_dp, 42.955_dp /)
         
               !Polynomial coefficients
               b(1:4,1) = (/  6.267e-1_dp , -2.647e-1_dp ,&
                              1.699e-1_dp ,  3.519e-1_dp  /)            !0th degree
               b(1:4,2) = (/ -1.311e-3_dp ,  1.576e-3_dp ,&
                              2.546e-4_dp , -1.310e-4_dp  /)            !1st degree
               b(1:4,3) = (/  1.344e-6_dp , -1.664e-6_dp ,&
                             -5.565e-8_dp , -1.836e-7_dp  /)            !2nd degree
               b(1:4,4) = (/ -5.190e-10_dp,  7.462e-10_dp,&
                             -8.521e-11_dp,  1.313e-10_dp /)            !3rd degree
               b(1:4,5) = (/  6.939e-14_dp, -1.200e-13_dp,&
                              2.516e-14_dp, -2.339e-14_dp /)            !4th degree

            elseif (which_total_pressure.eq.10._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.039_dp, 0.328_dp,& 
                                 3.337_dp, 40.591_dp /)   
         
               !Polynomial coefficients
               b(1:4,1) = (/  7.523e-1_dp , -1.919e-1_dp ,&  
                              4.111e-2_dp ,  4.328e-1_dp  /)            !0th degree
               b(1:4,2) = (/ -2.035e-3_dp ,  1.549e-3_dp ,&  
                              6.564e-4_dp , -2.644e-4_dp  /)            !1st degree
               b(1:4,3) = (/  2.323e-6_dp , -1.773e-6_dp ,& 
                             -3.944e-7_dp , -1.120e-7_dp  /)            !2nd degree
               b(1:4,4) = (/ -1.016e-9_dp ,  8.482e-10_dp,&  
                              3.100e-11_dp,  1.177e-10_dp /)            !3rd degree
               b(1:4,5) = (/  1.542e-13_dp, -1.424e-13_dp,&  
                              1.076e-14_dp, -2.298e-14_dp /)            !4th degree

            elseif (which_total_pressure.eq.20._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.024_dp, 0.211_dp,& 
                                 2.118_dp, 25.103_dp /)
         
               !Polynomial coefficients
               b(1:4,1) = (/  5.776e-1_dp ,  1.820e-1_dp ,& 
                             -2.658e-1_dp ,  5.821e-1_dp  /)            !0th degree
               b(1:4,2) = (/ -1.834e-3_dp ,  4.348e-4_dp ,&  
                              1.517e-3_dp , -3.830e-4_dp  /)            !1st degree
               b(1:4,3) = (/  2.253e-6_dp , -6.735e-7_dp ,& 
                             -1.212e-6_dp , -6.810e-8_dp  /)            !2nd degree
               b(1:4,4) = (/ -1.032e-9_dp ,  4.041e-10_dp,&  
                              3.715e-10_dp,  1.055e-10_dp /)            !3rd degree
               b(1:4,5) = (/  1.620e-13_dp, -7.788e-14_dp,& 
                             -4.108e-14_dp, -2.112e-14_dp /)            !4th degree

            elseif (which_total_pressure.eq.40._dp) then
               !Pressure-absorption coefficients
               kappa_p(1:4) = (/ 0.014_dp, 0.138_dp,& 
                                 1.274_dp, 15.934_dp /)
         
               !Polynomial coefficients
               b(1:4,1) = (/  2.357e-1_dp ,  6.246e-1_dp ,& 
                             -4.705e-1_dp ,  6.633e-1_dp  /)            !0th degree
               b(1:4,2) = (/ -8.825e-4_dp , -1.178e-3_dp ,&  
                              2.150e-3_dp , -3.019e-4_dp  /)            !1st degree
               b(1:4,3) = (/  1.222e-6_dp ,  1.131e-6_dp ,& 
                             -1.888e-6_dp , -1.908e-7_dp  /)            !2nd degree
               b(1:4,4) = (/ -5.762e-10_dp, -4.045e-10_dp,&  
                              6.854e-10_dp,  1.514e-10_dp /)            !3rd degree
               b(1:4,5) = (/  9.134e-14_dp,  4.924e-14_dp,& 
                             -9.249e-14_dp, -2.674e-14_dp /)            !4th degree
               
            else                                                        !Stop and print error message
               write(msg,'(a,f3.1,a,a)') &
                     'Total pressure ',which_total_pressure,&
                     ' not available for correlation ',trim(wsgg_spec)
               call shutdown(msg)
            endif
            
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Selhorst et al. (2020)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('selhorst2020')
            !Total number of gray gases
            number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 6
         
            !Pressure-absorption coefficients
            kappa_p(1:4) = (/ 0.269348_dp, 2.28822390707116_dp,&
                             15.1959768410186_dp, 156.930802145608_dp /)
      
            !Polynomial coefficients
            b(1,1:7) = (/  9.8470e-2_dp,  6.85690e-1_dp,&
                          -5.8071e-1_dp,  1.66572e-1_dp,&
                           0.0000e+0_dp,  0.00000e+0_dp,&
                          -1.8795e-3_dp /)
            b(2,1:7) = (/  1.1001e-1_dp,  3.09624e-1_dp,&
                          -2.6010e-1_dp,  5.78303e-2_dp,&
                           0.0000e+0_dp,  0.00000e+0_dp,&
                          -3.3529e-4_dp /)
            b(3,1:7) = (/  1.8455e-1_dp, -7.44010e-2_dp,&
                          -1.0004e-2_dp,  6.33634e-3_dp,&
                           0.0000e+0_dp,  0.00000e+0_dp,&
                          -3.0838e-5_dp /)
            b(4,1:7) = (/  1.1558e-1_dp, -1.07150e-1_dp,&
                           5.5662e-2_dp, -1.58230e-2_dp,&
                           0.0000e+0_dp,  0.00000e+0_dp,&
                           2.4449e-4_dp /)
      
         !case('mr2')
            !Total number of gray gases
            !number_wsgg_gases = 4
      
            !Degree of the polynomial for the temperature coefficient
            !degree_wsgg_polynomial = 6
         
            !Pressure-absorption coefficients
            !kappa_p(1:4) = (/ 0.275757_dp, 2.249395_dp, 13.955648_dp, 125.40315_dp /)
      
            !Polynomial coefficients
            !b(1,1:6) = (/ 1.1411E-01_dp, 4.9428E-01_dp, -3.6579E-01_dp,   1.02590E-01_dp,   0.0000E+00_dp,   0.0000E+00_dp, -1.2598E-03_dp /)
            !b(2,1:6) = (/ 1.0832E-01_dp, 3.19134E-01_dp, -2.3526E-01_dp, 4.42865E-02_dp, 0.0000E+00_dp,   0.0000E+00_dp, -1.2684E-04_dp /)
            !b(3,1:6) = (/ 1.4118E-01_dp, 1.29516E-01_dp, -2.0227E-01_dp, 6.19913E-02_dp, 0.0000E+00_dp,   0.0000E+00_dp, -5.7748E-04_dp /)
            !b(4,1:6) = (/ 1.5930E-01_dp, -1.9197E-01_dp, 1.0956E-01_dp,   -2.7972E-02_dp,   0.0000E+00_dp,   0.0000E+00_dp, 3.2486E-04_dp /)
      
      !endselect
      
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Coelho et al. (RAD19)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('rad19')
            !Total number of gray gases
            number_wsgg_gases = 4
         
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4
            
            selectcase(trim(soot_wsgg_model))
               case('lee')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 168837._dp , 579945._dp ,&
                                    2162920._dp, 4074590._dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/  1.0834e+0_dp , -1.5854e-1_dp ,&
                                -6.8794e-1_dp ,  3.2998e-1_dp  /)       !0th degree
                  b(1:4,2) = (/ -2.4344e-3_dp ,  2.8795e-3_dp ,&
                                 1.7757e-3_dp , -1.1141e-3_dp  /)       !1st degree
                  b(1:4,3) = (/  2.0689e-6_dp , -3.6549e-6_dp ,&
                                -5.1673e-7_dp ,  1.0536e-6_dp  /)       !2nd degree
                  b(1:4,4) = (/ -7.7367e-10_dp,  1.6514e-09_dp,&
                                -1.3346e-10_dp, -3.1503e-10_dp /)       !3rd degree
                  b(1:4,5) = (/  1.0675e-13_dp, -2.5477e-13_dp,&
                                 5.1100e-14_dp,  3.3178e-14_dp /)       !4th degree

               case('dalzell')           
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 168984._dp , 629048._dp ,&
                                    2387410._dp, 5488390._dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/  9.8752e-1_dp , -5.5401e-2_dp ,&
                                -6.2116e-1_dp ,  2.5764e-1_dp  /)       !0th degree
                  b(1:4,2) = (/ -2.2909e-3_dp ,  2.6849e-3_dp ,&
                                 1.5792e-3_dp , -8.7130e-4_dp  /)       !1st degree
                  b(1:4,3) = (/  1.9861e-6_dp , -3.5019e-6_dp ,&
                                -3.1940e-7_dp ,  7.9086e-7_dp  /)       !2nd degree
                  b(1:4,4) = (/ -7.5272e-10_dp,  1.5942e-09_dp,&
                                -2.1026e-10_dp, -2.0382e-10_dp /)       !3rd degree
                  b(1:4,5) = (/  1.0482e-13_dp, -2.4639e-13_dp,&
                                 6.0896e-14_dp,  1.7213e-14_dp /)       !4th degree

               case('chang')          
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 135017._dp , 338321._dp ,&
                                    1401810._dp, 3570560._dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/  2.9789e-1_dp ,  9.0433e-1_dp ,&
                                -7.6923e-1_dp ,  1.3030e-1_dp  /)       !0th degree
                  b(1:4,2) = (/ -9.9361e-4_dp , -7.6433e-4_dp ,&
                                 3.7483e-3_dp , -8.8092e-4_dp  /)       !1st degree
                  b(1:4,3) = (/  1.0947e-6_dp , -1.7600e-7_dp ,&
                                -3.2012e-6_dp ,  1.2368e-06_dp /)       !2nd degree
                  b(1:4,4) = (/ -4.7997e-10_dp,  2.8188e-10_dp,&
                                 1.1178e-9_dp , -4.9359e-10_dp /)       !3rd degree
                  b(1:4,5) = (/  7.2890e-14_dp, -5.8550e-14_dp,&
                                -1.4805e-13_dp,  7.0586e-14_dp /)       !4th degree

               case('linear')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 226829._dp , 642573._dp ,&
                                    1636970._dp, 3349290._dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/  8.1024e-1_dp ,  1.2792e-1_dp ,&
                                -5.7641e-1_dp ,  2.0548e-1_dp  /)       !0th degree
                  b(1:4,2) = (/ -1.9068e-3_dp ,  2.4642e-3_dp ,&
                                 1.1414e-3_dp , -5.9331e-4_dp  /)       !1st degree
                  b(1:4,3) = (/  1.6735e-6_dp , -3.4210e-6_dp ,&
                                 2.7016e-7_dp ,  4.3060e-7_dp  /)       !2nd degree
                  b(1:4,4) = (/ -6.3827e-10_dp,  1.5851e-9_dp ,&
                                -4.9365e-10_dp, -2.5253e-11_dp /)       !3rd degree
                  b(1:4,5) = (/  8.8293e-14_dp, -2.4334e-13_dp,&
                                 1.0031e-13_dp, -8.7781e-15_dp /)       !4th degree
               
               case('linear_cassol_2')
                  !Number of gray gases
                  number_wsgg_gases = 2
                  
                  !Pressure-absorption coefficients
                  kappa_p(1:2) = (/ 22313.49_dp*4.1_dp,&
                                    466624.8_dp*4.1_dp /)

                  !Polynomial coefficients
                  b(1:2,1) = (/  0.95552_dp  ,  0.08010_dp   /)         !0th degree
                  b(1:2,2) = (/ -1.431e-3_dp ,  1.290e-3_dp  /)         !1st degree
                  b(1:2,3) = (/  9.871e-7_dp , -7.874e-7_dp  /)         !2nd degree
                  b(1:2,4) = (/ -3.390e-10_dp,  2.322e-10_dp /)         !3rd degree
                  b(1:2,5) = (/  4.555e-14_dp, -3.084e-14_dp /)         !4th degree
                  
               case('linear_cassol_3')
                  !Number of gray gases
                  number_wsgg_gases = 3
                     
                  !Pressure-absorption coefficients
                  kappa_p(1:3) = (/ 1251.56_dp*4.1_dp, &
                                    50470.6_dp*4.1_dp, &
                                    460361.0_dp*4.1_dp /)
         
                  !Polynomial coefficients
                  b(1:3,1) = (/  0.02812_dp   ,   1.25626_dp ,&
                                -0.25179_dp  /)                         !0th degree
                  b(1:3,2) = (/ -1.2710e-4_dp , -18.74e-4_dp ,&
                                 18.55e-4_dp  /)                        !1st degree
                  b(1:3,3) = (/  1.3950e-7_dp ,  12.09e-7_dp ,&
                                -11.39e-7_dp  /)                        !2nd degree
                  b(1:3,4) = (/ -5.6720e-11_dp, -39.89e-11_dp,&
                                 34.46e-11_dp /)                        !3rd degree
                  b(1:3,5) = (/  8.0868e-15_dp,  52.85e-15_dp,& 
                                -45.58e-15_dp /)                        !4th degree
                  
               case('linear_cassol_4')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 2857.86_dp*4.1_dp, &
                                    39234.9_dp*4.1_dp, &
                                    160748._dp*4.1_dp, &
                                    495898._dp*4.1_dp /)
   
                  !Polynomial coefficients
                  b(1:4,1) = (/  0.00129_dp   ,    1.26110_dp,&
                                   -0.25757_dp,   0.07980_dp  /)        !0th degree
                  b(1:4,2) = (/ -0.5450e-5_dp , -319.2e-5_dp ,&
                                  362.1e-5_dp , -72.08e-5_dp  /)        !1st degree
                  b(1:4,3) = (/  0.1230e-7_dp ,   27.72e-7_dp,&
                                  -40.12e-7_dp,  15.87e-7_dp  /)        !2nd degree
                  b(1:4,4) = (/ -0.8470e-11_dp, -100.5e-11_dp,&  
                                  154.9e-11_dp, -70.89e-11_dp /)        !3rd degree
                  b(1:4,5) = (/  1.6807e-15_dp,  132.8e-15_dp,&
                                 -207.8e-15_dp,  97.69e-15_dp /)        !4th degree

               case default
                  write(msg,'(3(a))') 'Soot model ',&
                     trim(soot_wsgg_model),' not available'
                  call shutdown(msg)
            endselect
            
            !Correct the number of gray gases if transparent
            !windows for soot are to be neglected
            if (.not.soot_transparent_windows) &
               number_wsgg_gases = number_wsgg_gases - 1

         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Hosein et al. (2021)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('hosein2021')
            !Total number of gray gases
            number_wsgg_gases = 4
         
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 5
            
            !Reference temperature
            wsgg_reference_temperature = 1400._dp

            select case(trim(species_spec_aux))
               case('soot')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.4530e5_dp, 7.5475e5_dp, &
                                    2.0836e6_dp, 4.2113e6_dp /)   
!                  kp = kp/pressure  
                  
                  !Polynomial coefficients
                  b(1,1:6) = (/  1.8613_dp   , -7.7857_dp   ,&
                                 1.2809e1_dp ,  -1.0158e1_dp,&
                                 3.8717_dp   , -5.6880e-1_dp /)
                  b(2,1:6) = (/ -1.1374_dp   ,  1.0625e1_dp ,&
                                -2.1665e1_dp ,  1.8750e1_dp ,&
                                -7.4578_dp   ,  1.1219_dp    /)
                  b(3,1:6) = (/  2.5975e-1_dp, -2.9708_dp   ,&
                                 9.6830_dp   , -9.7681_dp   ,&
                                 4.1073_dp   , -6.3057e-1_dp /)
                  b(4,1:6) = (/ -3.8367e-2_dp,  4.2159e-1_dp,&
                                -1.4101_dp   ,  1.7225_dp   ,&
                                 -7.5889e-1_dp, 1.1454e-1_dp /)
       
               case('co')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.7920e-1_dp, 1.7918_dp, &
                                    1.2953e1_dp, 1.2900e2_dp /) 
                  
                  !Polynomial coefficients
                  b(1,1:6) = (/ -5.3582e-3_dp, -1.4397e-3_dp,&
                                 4.0604e-1_dp, -5.7254e-1_dp,&
                                 2.8282e-1_dp, -4.7820e-2_dp /)
                  b(2,1:6) = (/ -6.7961e-2_dp,  4.2204e-1_dp,&
                                -5.4894e-1_dp,  2.8819e-1_dp,&
                                -6.2318e-2_dp,  3.7321e-3_dp /)
                  b(3,1:6) = (/ -5.7642e-2_dp,  4.2020e-1_dp,&
                                -7.6297e-1_dp,  6.0302e-1_dp,&
                                -2.2181e-1_dp, -3.1122e-2_dp /)
                  b(4,1:6) = (/ -1.6152e-2_dp,  1.2220e-1_dp,&
                                -2.2207e-1_dp,  1.7430e-1_dp,&
                                -6.3464e-2_dp,  8.8012e-3_dp /)

               case('ch4')
                  !Reference temperature
                  wsgg_reference_temperature = 2500._dp
            
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.749379e-1_dp, 2.571313e1_dp, &
                                    2.961692e3_dp, 2.810566_dp /)

                  !Polynomial coefficients
                  b(1:4,1) = (/ -2.637795e-1_dp, -1.291970e-1_dp,&
                                 1.311677e-3_dp, -2.549570e-1_dp /)
                  b(1:4,2) = (/  3.462304_dp   ,  2.129878_dp   ,&
                                 1.847697e-1_dp,  3.559554_dp    /)
                  b(1:4,3) = (/ -3.845975_dp   , -8.190559_dp   ,&
                                -1.011002_dp   , -1.100239e1_dp  /)
                  b(1:4,4) = (/ -7.555724_dp   ,  1.343738e1_dp ,&
                                 2.020896_dp   ,  1.429368e1_dp  /)
                  b(1:4,5) = (/  1.499125e1_dp , -1.023326e1_dp ,&
                                -1.770474_dp   , -8.453295_dp    /)
                  b(1:4,6) = (/ -6.783667_dp   ,  2.988867_dp   ,&
                                 5.756737e-1_dp,  1.855989_dp    /)
                               
            endselect

         case('fraga2021')
            !Total number of gray gases
            number_wsgg_gases = 4
         
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4
            
            !Reference temperature
            wsgg_reference_temperature = 1000._dp

            select case(trim(species_spec_aux))
               case('h2o')                                              !Use value for x = 0.1
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.919439e-1_dp, 1.498957_dp,&
                                    9.803681_dp, 1.109096e2_dp /)

                  !Polynomial coefficients
                  b(1,1:5) = (/ 1.590220e-1_dp,  2.179190e-1_dp,&
                               -3.412370e-2_dp,  2.162570e-3_dp,&
                               -2.333450e-3_dp /)
                  b(2,1:5) = (/ 4.463570e-2_dp,  3.852140e-1_dp,&
                               -2.165200e-1_dp,  2.950680e-2_dp,&
                                8.475130e-4_dp /)
                  b(3,1:5) = (/ 9.873750e-3_dp,  5.614650e-1_dp,&
                               -6.756940e-1_dp,  2.882870e-1_dp,&
                               -4.280770e-2_dp /)
                  b(4,1:5) = (/ 1.980150e-1_dp, -3.538890e-1_dp,&
                                2.627760e-1_dp, -9.273280e-2_dp,&
                                1.266630e-2_dp /)

               case('co2')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.584383e-1_dp, 2.670235e+0_dp, &
                                    2.394601e+1_dp, 2.912417e+2_dp /)  

                  !Polynomial coefficients
                  b(1,1:5) = (/  1.228740e-1_dp,  5.472120e-1_dp,&
                                -7.693430e-1_dp,  3.698630e-1_dp,&
                                -6.083820e-2_dp /)
                  b(2,1:5) = (/  5.223340e-2_dp, -4.409430e-2_dp,&
                                 1.666480e-1_dp, -1.128020e-1_dp,&
                                 2.128700e-2_dp /)
                  b(3,1:5) = (/  1.190670e-1_dp, -2.321120e-1_dp,&
                                 2.634100e-1_dp, -1.251010e-1_dp,&
                                 2.025260e-2_dp /)
                  b(4,1:5) = (/  2.285260e-3_dp,  1.334090e-1_dp,&
                                -1.254240e-1_dp,  3.722150e-2_dp,&
                                -3.179340e-3_dp /)
                  
               case('soot')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.631470e5_dp, 1.000450e6_dp, &
                                    2.953190e6_dp, 2.979880e6_dp /)

                  !Polynomial coefficients
                  b(1,1:5) = (/  1.450490_dp   , -2.907530_dp   ,&
                                 2.108700_dp   , -6.586220e-1_dp,&
                                 7.539150e-2_dp /)
                  b(2,1:5) = (/ -7.972750e-1_dp,  4.673440_dp   ,&
                                -4.797280_dp   ,  1.912160_dp   ,&
                                -2.740940e-1_dp /)
                  b(3,1:5) = (/  7.397080e-1_dp, -3.790640_dp   ,&
                                 5.861280_dp   , -3.089560_dp   ,&
                                 5.209490e-1_dp /)
                  b(4,1:5) = (/ -5.680900e-1_dp,  2.756410_dp   ,&
                                -4.109730_dp   ,  2.300270_dp   ,&
                                -4.005500e-1_dp /)

               case('co')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.910752e-1_dp, 1.664792e+0_dp, &
                                    1.238668e+1_dp, 1.882960e+2_dp /)

                  !Polynomial coefficients
                  b(1,1:5) = (/ -1.893390e-2_dp,  9.334080e-2_dp,&
                                 3.768110e-3_dp, -3.503150e-2_dp,&
                                 8.796540e-3_dp /)
                  b(2,1:5) = (/ -1.036400e-1_dp,  4.394090e-1_dp,&
                                -4.264890e-1_dp,  1.616730e-1_dp,&
                                -2.163050e-2_dp /)
                  b(3,1:5) = (/ -4.956480e-2_dp,  2.621080e-1_dp,&
                                -3.089880e-1_dp,  1.397920e-1_dp,&
                                -2.193810e-2_dp /)
                  b(4,1:5) = (/ -1.395890e-2_dp, 7.300760e-2_dp ,&
                                -8.310090e-2_dp,  3.628900e-2_dp,&
                                -5.522170e-3_dp /)

               case('ch4')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.611520e-1_dp, 1.152081e+0_dp, &
                                    7.220810e+0_dp, 5.986229e+1_dp /)

                  !Polynomial coefficients
                  b(1,1:5) = (/ -2.576940e-1_dp,  1.392900_dp   ,&
                                -1.107260_dp   ,  4.189890e-1_dp,&
                                -6.624220e-2_dp /)
                  b(2,1:5) = (/ -1.252560e-1_dp,  5.877710e-1_dp,&
                                -2.905960e-1_dp,  1.640850e-2_dp,&
                                 7.288570e-3_dp /)
                  b(3,1:5) = (/ -6.488220e-2_dp,  3.589360e-1_dp,&
                                -1.308920e-1_dp, -3.910640e-2_dp,&
                                 1.588830e-2_dp /)
                  b(4,1:5) = (/ -1.208290e-2_dp,  2.081460e-1_dp,&
                                -2.559860e-1_dp,  1.081910e-1_dp,&
                                -1.542060e-2_dp /)

            endselect

         case('fraga2021v2')
            !Total number of gray gases
            number_wsgg_gases = 4
         
            !Degree of the polynomial for the temperature coefficient
            degree_wsgg_polynomial = 4
            
            !Reference temperature
            wsgg_reference_temperature = 1000._dp

            select case(trim(species_spec_aux))
               case('h2o')                                              !Use value for x = 0.1
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.919439e-1_dp, 1.498957_dp,&
                                    9.803681_dp, 1.109096e2_dp /)

                  !Polynomial coefficients
                  b(1,1:5) = (/ 1.590220e-1_dp,  2.179190e-1_dp,&
                               -3.412370e-2_dp,  2.162570e-3_dp,&
                               -2.333450e-3_dp /)
                  b(2,1:5) = (/ 4.463570e-2_dp,  3.852140e-1_dp,&
                               -2.165200e-1_dp,  2.950680e-2_dp,&
                                8.475130e-4_dp /)
                  b(3,1:5) = (/ 9.873750e-3_dp,  5.614650e-1_dp,&
                               -6.756940e-1_dp,  2.882870e-1_dp,&
                               -4.280770e-2_dp /)
                  b(4,1:5) = (/ 1.980150e-1_dp, -3.538890e-1_dp,&
                                2.627760e-1_dp, -9.273280e-2_dp,&
                                1.266630e-2_dp /)

               case('co2')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 8.151690e-2_dp, 9.787520e-1_dp, &
                                    8.510160_dp, 1.911380e2_dp /) 

                  !Polynomial coefficients
                  b(1,1:5) = (/  2.0693e-1_dp,  5.2037e-1_dp,&
                                -8.5223e-1_dp,  4.3860e-1_dp,&
                                -7.5328e-2_dp /)
                  b(2,1:5) = (/ -8.8737e-3_dp,  2.0370e-1_dp,&
                                -1.4561e-1_dp,  3.4399e-2_dp,&
                                -2.1191e-3_dp /)
                  b(3,1:5) = (/  1.4750e-1_dp, -3.1930e-1_dp,&
                                 4.1047e-1_dp, -2.0347e-1_dp,&
                                 3.3415e-2_dp /)
                  b(4,1:5) = (/  2.8914e-2_dp,  9.9319e-2_dp,&
                                -8.9300e-2_dp,  1.8716e-2_dp,&
                                 6.0005e-5_dp /)

               case('soot')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.579500e5_dp, 9.851430e5_dp, &
                                    2.915820e6_dp, 3.005000e6_dp /)

                  !Polynomial coefficients
                  b(1,1:5) = (/  1.712830_dp   , -3.798550_dp   ,&
                                 3.132310_dp   , -1.143010_dp   ,&
                                 1.563750e-1_dp /)
                  b(2,1:5) = (/ -7.558450e-1_dp,  4.704080_dp   ,&
                                -4.976680_dp   ,  2.044010_dp   ,&
                                -3.015440e-1_dp /)
                  b(3,1:5) = (/  5.246250e-1_dp, -2.755590_dp   ,&
                                 4.346840_dp   , -2.236920_dp   ,&
                                 3.604750e-1_dp /)
                  b(4,1:5) = (/ -3.471090e-1_dp,  1.672340_dp   ,&
                                -2.483080_dp   ,  1.374970_dp   ,&
                                -2.259010e-1_dp /)

               case('co')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.563010e-1_dp, 9.291450e-1_dp, &
                                    4.548900_dp, 7.094230e1_dp /)

                  !Polynomial coefficients
                  b(1,1:5) = (/ -9.448530e-3_dp,  7.857400e-2_dp,&
                                -1.203580e-2_dp, -1.612550e-2_dp,&
                                 4.408760e-3_dp /)
                  b(2,1:5) = (/ -6.492880e-2_dp,  2.374770e-1_dp,&
                                -1.762510e-1_dp,  4.428570e-2_dp,&
                                -2.752670e-3_dp /)
                  b(3,1:5) = (/ -8.312190e-2_dp,  4.121860e-1_dp,&
                                -4.697440e-1_dp,  2.072590e-1_dp,&
                                -3.189440e-2_dp /)
                  b(4,1:5) = (/ -2.723380e-2_dp,  1.442960e-1_dp,&
                                -1.662020e-1_dp,  7.344460e-2_dp,&
                                -1.129700e-2_dp /)

               case('ch4')
                  !Pressure-absorption coefficients
                  kappa_p(1:4) = (/ 1.611520e-1_dp, 1.152081_dp, &
                                    7.220810_dp, 5.986229e1_dp /)

                  !Polynomial coefficients
                  b(1,1:5) = (/ -2.576940e-1_dp,  1.392900_dp   ,&
                                -1.107260_dp   ,  4.189890e-1_dp,&
                                -6.624220e-2_dp /)
                  b(2,1:5) = (/ -1.252560e-1_dp,  5.877710e-1_dp,&
                                -2.905960e-1_dp,  1.640850e-2_dp,&
                                 7.288570e-3_dp /)
                  b(3,1:5) = (/ -6.488220e-2_dp,  3.589360e-1_dp,&
                                -1.308920e-1_dp, -3.910640e-2_dp,&
                                 1.588830e-2_dp /)
                  b(4,1:5) = (/ -1.208290e-2_dp,  2.081460e-1_dp,&
                                -2.559860e-1_dp,  1.081910e-1_dp,&
                                -1.542060e-2_dp /)

            endselect
         case default
            write(msg,'(3(a))') 'Correlation ',trim(wsgg_spec),&
                  ' not available'
            call shutdown(msg)
      endselect
         
   endsubroutine get_wsgg_correlations   
   
endmodule wsgg_parameters
