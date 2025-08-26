!#######################################################################
!Functions and subroutines related to the 
!radiative transfer solution via the WSGG model
!#######################################################################
module wsgg_functions

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp,small
   use wsgg_parameters
   
   !====================================================================
   !Declaration of variables
   !====================================================================
   implicit none
   real(dp),allocatable,dimension(:,:) :: superposition_array

contains

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                      !
!                             MAIN FUNCTIONS                           !
!                                                                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   !====================================================================
   !Function to get the number of gray gases of the model
   !====================================================================
   integer recursive function wsgg_number_gray_gases(wsgg_spec,&
      species_spec) result(number_gray_gases)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      character(*),intent(in) :: wsgg_spec
      character(*),optional :: species_spec
      character(100) :: species_spec_aux
      integer :: ngg,ngg_co2,ngg_h2o,ngg_mix,ngg_soot
      character(200) :: msg
      
      !-----------------------------------------------------------------
      !Default parameters
      !-----------------------------------------------------------------
      species_spec_aux = species_spec
      if (.not.present(species_spec)) species_spec_aux = 'mixture'
      
      !-----------------------------------------------------------------
      !Getting the number of gray gases
      !-----------------------------------------------------------------
      selectcase(trim(wsgg_spec))
         case('smith1982')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('smith1982_stepwise')
            ngg = 3
            
         case('smith1982_linear')
            ngg = 3
         
         case('soufiani1994')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
            
         case('bahador2008')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('johansson2010_4gg')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('johansson2010_3gg')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
            
         case('krishnamoorthy2010a')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('krishnamoorthy2010b')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
            
         case('yin2010')
            ngg = 4
         
         case('yin2010_stepwise')
            ngg = 4 
         
         case('yin2010_linear')
            ngg = 4 
           
         case('johansson2011')
            ngg = 4
      
         case('kangwanpongpan2012_fixed')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('kangwanpongpan2012_variable')
            ngg = 4
         
         case('kangwanpongpan2012_stepwise')
            ngg = 4
         
         case('kangwanpongpan2012_linear')
            ngg = 4
         
         case('sun2012')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('ziemniczak2013_3gg')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('ziemniczak2013_4gg')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('ziemniczak2013_5gg')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
            
         case('dorigon2013')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('krishnamoorthy2013')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
         
         case('krishnamoorthy2013_stepwise')
            ngg = 4
         
         case('krishnamoorthy2013_linear')
            ngg = 4
            
         case('yin2013')
            ngg = 4
         
         case('yin2013_stepwise')
            ngg = 4
         
         case('yin2013_linear')
            ngg = 4
         
         case('bordbar2014')
            ngg = 4
   
         case('cassol2014')
            call get_wsgg_correlations(wsgg_spec,species_spec_aux)
            ngg = number_wsgg_gases

         case('guo2015')
            ngg = 4
         
         case('ge2017')
            ngg = 4
            
         case('coelho2017')
            call get_wsgg_correlations(wsgg_spec,species_spec_aux)
            ngg = number_wsgg_gases
            
         case('centeno2018')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases
      
         case('marzouk2018')
            call get_wsgg_correlations(wsgg_spec,species_spec_aux)
            ngg = number_wsgg_gases
      
         case('shan2018')
            ngg = 4
         
         case('shan2018_stepwise')
            ngg = 4
         
         case('shan2018_linear')
            ngg = 4
            
         case('coelho2019')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases

         case('selhorst2020')
            call get_wsgg_correlations(wsgg_spec)
            ngg = number_wsgg_gases

         case('rad19')
            call get_wsgg_correlations(wsgg_spec)
            ngg = (number_wsgg_gases+1)*5 - 1
         
         case('fatmir2024')
            ngg = 6-1
         
         case('superposition')
            if (trim(wsgg_mixture_spec).ne.'null') then                 !Gaseous mixture + soot
               ngg_mix = wsgg_number_gray_gases(wsgg_mixture_spec)
               if (trim(wsgg_soot_spec).ne.'null') then
                  ngg_soot = wsgg_number_gray_gases(wsgg_soot_spec,&
                                                    'soot')
               else
                  ngg_soot = 0
               endif   
               ngg = (ngg_mix + 1)*(ngg_soot + 1) - 1
            
            else                                                        !CO2 + H2O + soot
               ngg_co2 = wsgg_number_gray_gases(wsgg_co2_spec,&
                                                'co2')
               ngg_h2o = wsgg_number_gray_gases(wsgg_h2o_spec,&
                                                'h2o')
               if (trim(wsgg_soot_spec).ne.'null') then
                  ngg_soot = wsgg_number_gray_gases(wsgg_soot_spec,&
                                                    'soot')
               else
                  ngg_soot = 0
               endif   
               ngg = (ngg_co2 + 1)*(ngg_h2o + 1)*(ngg_soot + 1) - 1
            endif
          
         case('fraga2021')
            call get_wsgg_correlations(wsgg_spec,species_spec_aux)
            ngg = number_wsgg_gases
            
         case('fraga2021v2')
            call get_wsgg_correlations(wsgg_spec,species_spec_aux)
            ngg = number_wsgg_gases
            
         case default
            write(msg,'(a,a,a)') 'wsgg_number_gray_gases: Model ',&
                                 trim(wsgg_spec),' not available'
            call shutdown(msg)
      
      endselect
      number_gray_gases = ngg
      
   endfunction wsgg_number_gray_gases

!   integer function wsgg_species_ngg(wsgg_species_spec,species_name)
   
!      use comp_functions, only: shutdown
!      implicit none
!      character(*),intent(in) :: wsgg_species_spec,species_name
!      integer :: ngg
      
!      selectcase(trim(wsgg_species_spec))
!      case('null')
!         ngg = 0
!      case('cassol2014')
!         ngg = 4
!         if (trim(species_name).eq.'soot') ngg = number_wsgg_gases_soot
!      case('hosein2021')
!         ngg = 4
!      case('fraga2021')
!         selectcase(trim(species_name))
!         case('h2o')
!            ngg = 4
!         case('co2')
!            ngg = 4
!         case('soot')
!            ngg = 4
!         case('co')
!            ngg = 4
!         case('ch4')
!            ngg = 4
!         case default
!            call shutdown('wsgg_species_ngg: species_name &
!                          &not available')
!         endselect
!      case('fraga2021v2')
!         selectcase(trim(species_name))
!         case('co2')
!            ngg = 4
!         case('soot')
!            ngg = 4
!         case('co')
!            ngg = 4
!         case default
!            call shutdown('wsgg_species_ngg: species_name &
!                          &not available')
!         endselect
!      case default
!         call shutdown('wsgg_species_ngg: wsgg_species_spec &
!                       &not available')
!      endselect

!      wsgg_species_ngg = ngg

!   endfunction wsgg_species_ngg



   !====================================================================
   !Subroutine to compute the absorption and weighting coefficients
   !for a single gray gas (for either a single species or a mixture)
   !====================================================================
   !Inputs:
   !  j: gray gas (anything outside the range of gray gases for the
   !               correlation signifies a transparent window)
   !  Tg: local temperature
   !  xg: local species mole fractions (array)
   !  pg: local total pressure
   !  wsgg_spec: string identifying the correlation
   !  species_spec: optional parameter to indentify a specific species
   !                (default: 'mixture')
   !Outputs:
   !  kappa_j: absorption coefficient of gas j
   !  a_j: weighting coefficient of gas j
   subroutine wsgg_properties(kappa_j,a_j,j,Tg,xg,pg,wsgg_spec,&
                              species_spec)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use global_parameters, only: id_ch4,id_co,id_co2,id_h2o,id_soot,&
                                   number_of_species,soot_constant
      use precision_parameters, only: dp
      implicit none
      character(*),intent(in) :: wsgg_spec
      character(*),intent(in),optional :: species_spec
      character(200) :: msg,species_spec_aux
      integer :: nsp
      integer,intent(in) :: j
      real(dp),intent(in) :: pg,Tg,xg(:)
      real(dp),intent(out) :: a_j,kappa_j
      real(dp) :: mr,pa,pch4,pco,pco2,ph2o,Tref,xch4,xco,xco2,xh2o,xsoot

      !-----------------------------------------------------------------
      !Preperatory procedures 
      !-----------------------------------------------------------------
      !Set optional paramer
      species_spec_aux = 'mixture'
      if (present(species_spec)) species_spec_aux = species_spec
      
      !Mole fractions of the relevant species      
      nsp = number_of_species
      xco2 = 0._dp; xh2o = 0._dp; xsoot = 0._dp; xco = 0._dp; xch4 = 0._dp
      if ((id_co2.le.nsp) .and.(id_co2.gt.0) ) xco2  = xg(id_co2)
      if ((id_h2o.le.nsp) .and.(id_h2o.gt.0) ) xh2o  = xg(id_h2o)
      if ((id_soot.le.nsp).and.(id_soot.gt.0)) xsoot = xg(id_soot)
      if ((id_co.le.nsp)  .and.(id_co.gt.0)  ) xco   = xg(id_co)
      if ((id_ch4.le.nsp) .and.(id_ch4.gt.0) ) xch4  = xg(id_ch4)
      
      !Mole fraction ratio
      mr = xh2o/(xco2+small)
      
      !Partial pressure
      pco2 = pg*xco2; ph2o = pg*xh2o; pco = pg*xco; pch4 = pg*xch4
      pa = pco2 + ph2o + pco + pch4
      
      !Default reference temperature value
      Tref = 1._dp

      !-----------------------------------------------------------------
      !Select the appropriate function to compute 
      !the gray gas absorption coefficient
      !-----------------------------------------------------------------
      selectcase(trim(wsgg_spec))
         case('smith1982')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('smith1982_stepwise')
            kappa_j = smith_wsgg_kappa(pco2,ph2o,mr,j,'stepwise')
            a_j = smith_wsgg_a(Tg,pco2,ph2o,mr,j,'stepwise')
         
         case('smith1982_linear')
            kappa_j = smith_wsgg_kappa(pco2,ph2o,mr,j,'linear')
            a_j = smith_wsgg_a(Tg,pco2,ph2o,mr,j,'linear')
         
         case('soufiani1994')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('bahador2008')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('johansson2010_4gg')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            Tref = 1200._dp
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('johansson2010_3gg')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            Tref = 1200._dp
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('krishnamoorthy2010a')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('krishnamoorthy2010b')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
            
         case('yin2010')
            kappa_j = yin2010_wsgg_kappa(pco2,ph2o,which_molar_ratio,j,&
                                         'none')
            a_j = yin2010_wsgg_a(Tg,pco2,ph2o,which_molar_ratio,j,'none')
         
         case('yin2010_stepwise')
            kappa_j = yin2010_wsgg_kappa(pco2,ph2o,mr,j,'stepwise')
            a_j = yin2010_wsgg_a(Tg,pco2,ph2o,mr,j,'stepwise')
         
         case('yin2010_linear')
            kappa_j = yin2010_wsgg_kappa(pco2,ph2o,mr,j,'linear')
            a_j = yin2010_wsgg_a(Tg,pco2,ph2o,mr,j,'linear')
         
         case('johansson2011')
            kappa_j = johansson_wsgg_kappa(pa,mr,j)
            a_j = johansson_wsgg_a(Tg,mr,j)

         case('kangwanpongpan2012_fixed')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
            
         case('kangwanpongpan2012_variable')
            kappa_j = kangwanpongpan_variable_wsgg_kappa(pa,mr,j)
            a_j = kangwanpongpan_variable_wsgg_a(Tg,mr,j)
         
         case('kangwanpongpan2012_stepwise')
            kappa_j = kangwanpongpan_fixed_wsgg_kappa(pco2,ph2o,mr,j,&
                                                      'stepwise')
            a_j = kangwanpongpan_fixed_wsgg_a(Tg,pco2,ph2o,mr,j,'stepwise')
         
         case('kangwanpongpan2012_linear')
            kappa_j = kangwanpongpan_fixed_wsgg_kappa(pco2,ph2o,mr,j,&
                                                      'linear')
            a_j = kangwanpongpan_fixed_wsgg_a(Tg,pco2,ph2o,mr,j,'linear')
         
         case('sun2012')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            Tref = 1200._dp
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('ziemniczak2013_3gg')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('ziemniczak2013_4gg')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
            
         case('ziemniczak2013_5gg')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
            
         case('dorigon2013')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('krishnamoorthy2013')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('krishnamoorthy2013_stepwise')
            kappa_j = krishnamoorthy2013_wsgg_kappa(pco2,ph2o,mr,j,&
                                                    'stepwise')
            a_j = krishnamoorthy2013_wsgg_a(Tg,pco2,ph2o,mr,j,'stepwise')
         
         case('krishnamoorthy2013_linear')
            kappa_j = krishnamoorthy2013_wsgg_kappa(pco2,ph2o,mr,j,&
                                                    'linear')
            a_j = krishnamoorthy2013_wsgg_a(Tg,pco2,ph2o,mr,j,'linear')
         
         case('yin2013')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            Tref = 1200._dp
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('yin2013_stepwise')
            kappa_j = yin2013_wsgg_kappa(pco2,ph2o,mr,j,'stepwise')
            a_j = yin2013_wsgg_a(Tg,pco2,ph2o,mr,j,'stepwise')
         
         case('yin2013_linear')
            kappa_j = yin2013_wsgg_kappa(pco2,ph2o,mr,j,'linear')
            a_j = yin2013_wsgg_a(Tg,pco2,ph2o,mr,j,'linear')

         case('bordbar2014')
            kappa_j = bordbar_wsgg_kappa(pa,mr,j)
            a_j = bordbar_wsgg_a(Tg,mr,j)

         case('cassol2014')
            if (trim(species_spec_aux).eq.'co2') &
               kappa_j = fixed_wsgg_kappa(pco2,j,'cassol2014',&
                                          'co2')
            if (trim(species_spec_aux).eq.'h2o') &
               kappa_j = fixed_wsgg_kappa(ph2o,j,'cassol2014',&
                                          'h2o')
            if (trim(species_spec_aux).eq.'soot') &
               kappa_j = fixed_wsgg_kappa(xsoot*soot_constant,j,&
                                          'cassol2014','soot')
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec,species_spec_aux)
                  
         case('guo2015')
            kappa_j = guo_wsgg_kappa(Tg,pa,mr,j)
            a_j = guo_wsgg_a(Tg,j)
         
         case('ge2017')
            kappa_j = ge_wsgg_kappa(Tg,pa,mr,j)
            a_j = ge_wsgg_a(Tg,mr,j)
            
         case('coelho2017')
            if (trim(species_spec_aux).eq.'co2') &
               kappa_j = fixed_wsgg_kappa(pco2,j,'coelho2017',&
                                          'co2')
            if (trim(species_spec_aux).eq.'h2o') &
               kappa_j = fixed_wsgg_kappa(ph2o,j,'coelho2017',&
                                          'h2o')
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec,species_spec_aux)
               
         case('centeno2018')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('marzouk2018')
            if (trim(species_spec_aux).eq.'co2') &
               kappa_j = fixed_wsgg_kappa(pco2,j,'marzouk2018',&
                                          'co2')
            if (trim(species_spec_aux).eq.'h2o') &
               kappa_j = fixed_wsgg_kappa(ph2o,j,'marzouk2018',&
                                          'h2o')
            if (trim(species_spec_aux).eq.'mixture') &
               kappa_j = fixed_wsgg_kappa(pa,j,'marzouk2018','mixture')
               a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec,species_spec_aux)
         
         case('shan2018')
            kappa_j = shan_wsgg_kappa(pa,pg,which_molar_ratio,j,'none')
            a_j = shan_wsgg_a(Tg,pg,which_molar_ratio,j,'none')
         
         case('shan2018_stepwise')
            kappa_j = shan_wsgg_kappa(pa,pg,mr,j,'stepwise')
            a_j = shan_wsgg_a(Tg,pg,mr,j,'stepwise')
         
         case('shan2018_linear')
            kappa_j = shan_wsgg_kappa(pa,pg,mr,j,'linear')
            a_j = shan_wsgg_a(Tg,pg,mr,j,'linear')
            
         case('coelho2019')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
         
         case('selhorst2020')
            kappa_j = fixed_wsgg_kappa(pa,j,wsgg_spec)
            Tref = 1000._dp
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec)
          
         case('fraga2021')
            if (trim(species_spec_aux).eq.'co2') &
               kappa_j = fixed_wsgg_kappa(pco2,j,'fraga2021','co2')
            if (trim(species_spec_aux).eq.'h2o') &
               kappa_j = fixed_wsgg_kappa(ph2o,j,'fraga2021','h2o')
            if (trim(species_spec_aux).eq.'soot') &
               kappa_j = fixed_wsgg_kappa(xsoot,j,'fraga2021','soot')
            if (trim(species_spec_aux).eq.'co') &
               kappa_j = fixed_wsgg_kappa(pco,j,'fraga2021','co')
            if (trim(species_spec_aux).eq.'ch4') &
               kappa_j = fixed_wsgg_kappa(pch4,j,'fraga2021','ch4')
            Tref = 1000._dp
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec,species_spec_aux)
            
         case('fraga2021v2')
            if (trim(species_spec_aux).eq.'co2') &
               kappa_j = fixed_wsgg_kappa(pco2,j,'fraga2021v2','co2')
            if (trim(species_spec_aux).eq.'h2o') &
               kappa_j = fixed_wsgg_kappa(ph2o,j,'fraga2021v2','h2o')
            if (trim(species_spec_aux).eq.'soot') &
               kappa_j = fixed_wsgg_kappa(xsoot,j,'fraga2021v2','soot')
            if (trim(species_spec_aux).eq.'co') &
               kappa_j = fixed_wsgg_kappa(pco,j,'fraga2021v2','co')
            if (trim(species_spec_aux).eq.'ch4') &
               kappa_j = fixed_wsgg_kappa(pch4,j,'fraga2021v2','ch4')
            Tref = 1000._dp
            a_j = fixed_wsgg_a(Tg,Tref,j,wsgg_spec,species_spec_aux)
!            write(*,*) kappa_j,a_j,pco2,ph2o,xsoot,pco,pch4,trim(species_spec_aux)
         
!         case('rad19')
!            kappa_j = rad19_wsgg_kappa(pco2,ph2o,xsoot,j)
!            a_j = rad19_wsgg_a(Tg,mr,j)
         
!         case('fatmir2024')
!            call wsgg_asllanaj2024(kappa_j,a_j,j,tg,xg,pg)
         
!         case('superposition')
!            call find_superposition_indexes(j,jm,jc,jw,js)              !Find the indexes of each correlation      
!            if (trim(wsgg_mixture_spec).ne.'null') then                 !Gaseous mixture + soot
!               kappa_mix = wsgg_kappa_func(Tg,pg,xg,jm,&
!                                          wsgg_mixture_spec,'mixture')
!               if (trim(wsgg_soot_spec).ne.'null') then
!                  kappa_soot = wsgg_kappa_func(Tg,pg,xg,js,&
!                                          wsgg_soot_spec,'soot')
!               else
!                  kappa_soot = 0._dp
!               endif              
!               kappa = kappa_mix + kappa_soot
       
!            else                                                        !CO2 + H2O + soot
!               kappa_co2 = wsgg_kappa_func(ttmp,tp,xxs,jc,&
!                                       wsgg_co2_spec,'co2')
!               kappa_h2o = wsgg_kappa_func(ttmp,tp,xxs,jw,&
!                                          wsgg_h2o_spec,'h2o')
!               if (trim(wsgg_soot_spec).ne.'null') then
!                  kappa_soot = wsgg_kappa_func(ttmp,tp,xxs,js,&
!                                               wsgg_soot_spec,'soot')
!               else
!                  kappa_soot = 0._dp
!               endif              
!               kappa = kappa_co2 + kappa_h2o + kappa_soot
!            endif
            
         case default
            write(msg,'(a,a,a)') 'wsgg_properties: Model ',&
                                 trim(wsgg_spec),' not available'
            call shutdown(msg)      
      
      endselect
      
      !------------------------------------------------------------
      !Apply the correction for the soot contribution if necessary
      !------------------------------------------------------------
      if (soot_gg_correction) kappa_j = kappa_j + &
            gg_kappa_soot(Tg,xsoot,soot_gg_correction_model)

   endsubroutine wsgg_properties

   
!   !====================================================================
!   !Function to compute the absorption coefficient
!   !for WSGG models for a single species
!   !====================================================================
!   real(dp) function wsgg_species_kappa(mole_frac,pressure,jgas,&
!                                        wsgg_species_spec,species_name)

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use global_parameters, only: soot_constant
!      use comp_functions, only: shutdown
!      implicit none
!      character(*),intent(in) :: wsgg_species_spec,species_name
!      integer,intent(in) :: jgas
!      real(dp),intent(in) :: mole_frac,pressure
!      real(dp) :: kp(1:10)

!      !-----------------------------------------------------------------
!      !Defining the pressure absorption coefficient
!      !-----------------------------------------------------------------
!      selectcase(trim(wsgg_species_spec))
!      case('null')
!         kp = 0._dp
         
!      case('cassol2014')
!         selectcase(trim(species_name))
!         case('h2o')
!            kp(1:4) = (/ 0.171_dp, 1.551_dp, 5.562_dp, 49.159_dp /)
!         case('co2')
!            kp(1:4) = (/ 0.138_dp, 1.895_dp, 13.301_dp, 340.811_dp /)
!         case('soot')
!            select case(number_wsgg_gases_soot)
!            case(2)
!               kp(1:2) = (/ 22313.49_dp, 466624.8_dp /)
!            case(3)
!               kp(1:3) = (/ 1251.56_dp, 50470.6_dp, 460361.0_dp /)
!            case(4)
!               kp(1:4) = (/ 2857.86_dp, 39234.9_dp, &
!                            160748.0_dp, 495898.0_dp/)
!            case default
!               call shutdown('wsgg_species_kappa: &
!                             &number_wsgg_gases_soot not available')
!            endselect           
            
!            kp = kp*soot_constant/pressure
!         case default
!            call shutdown('wsgg_species_kappa: species_name &
!                          &not available')
!         endselect
         
!      case('hosein2021')
!         selectcase(trim(species_name))
!         case('soot')
!            kp(1:4) = (/ 1.4530e5_dp, 7.5475e5_dp, &
!                         2.0836e6_dp, 4.2113e6_dp /)   
!            kp = kp/pressure         
!         case('co')
!            kp(1:4) = (/ 1.7920e-1_dp, 1.7918_dp, &
!                         1.2953e1_dp, 1.2900e2_dp /)              
!         case('ch4')
!            kp(1:4) = (/ 1.749379e-1_dp, 2.571313e1_dp, &
!                         2.961692e3_dp, 2.810566_dp /)
!            !(/ 1.6352e-1_dp, 2.3385_dp, 1.7250e1_dp, 1.2668e2_dp /)
!         endselect
         
!      case('fraga2021')
!         selectcase(trim(species_name))
!         case('h2o')
!            kp = fraga_wsgg_kappa(mole_frac,jgas,&
!                                  fraga_h2o_interpolation)
!         case('co2')
!            kp(1:4) = (/ 1.584383e-1_dp, 2.670235_dp, &
!                         2.394601e1_dp, 2.912417e2_dp /)            
!         case('soot')
!            kp(1:4) = (/ 1.631470e5_dp, 1.000450e6_dp, &
!                         2.953190e6_dp, 2.979880e6_dp /)
!         case('co')
!            kp(1:4) = (/ 1.910752e-1_dp, 1.664792_dp, &
!                         1.238668e1_dp, 1.882960e2_dp /)
!         case('ch4')
!            kp(1:4) = (/ 1.611520e-1_dp, 1.152081_dp, &
!                         7.220810_dp, 5.986229e1_dp /)
!         case default
!            call shutdown('wsgg_species_kappa: species_name &
!                          &not available')
!         endselect
         
!      case('fraga2021v2')
!         selectcase(trim(species_name))
!         case('co2')
!            kp(1:4) = (/ 8.151690e-2_dp, 9.787520e-1_dp, &
!                         8.510160_dp, 1.911380e2_dp /)
!         case('soot')
!            kp(1:4) = (/ 1.579500e5_dp, 9.851430e5_dp, &
!                         2.915820e6_dp, 3.005000e6_dp /)
!         case('co')
!            kp(1:4) = (/ 1.563010e-1_dp, 9.291450e-1_dp, &
!                         4.548900_dp, 7.094230e1_dp /)
!         case default
!            call shutdown('wsgg_species_kappa: species_name &
!                          &not available')
!         endselect
         
!      case default
!         call shutdown('wsgg_species_kappa: wsgg_species_spec &
!                       &not available')
!      endselect

!      !-----------------------------------------------------------------
!      !Computing the absorption coefficient
!      !-----------------------------------------------------------------
!      if (jgas.eq.0) then
!         wsgg_species_kappa = 0._dp
!      else
!         wsgg_species_kappa = kp(jgas)*mole_frac*pressure
!      endif
      
!   endfunction wsgg_species_kappa
   
!   !====================================================================
!   !Function to compute the weighting coefficient
!   !for WSGG models for a single species
!   !====================================================================
!   real(dp) recursive function wsgg_species_a(tmp,mole_frac,jgas,&
!      wsgg_species_spec,species_name) result (a_wsgg)
      
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: shutdown
!      implicit none
!      character(*),intent(in) :: wsgg_species_spec,species_name
!      integer,intent(in) :: jgas
!      integer :: i,ndegree,ngg
!      real(dp),intent(in) :: mole_frac,tmp
!      real(dp) :: b(1:10,1:10),asum,tref
      
!      !Reference temperature (=1 by default)
!      tref = 1._dp
      
!      !-----------------------------------------------------------------
!      !Defining the polynomial coefficients
!      !-----------------------------------------------------------------
!      selectcase(trim(wsgg_species_spec))
!      case('null')
!         b = 0._dp

!      case('cassol2014')
!         ndegree = 4                                                    !Polynomial degree      
!         selectcase(trim(species_name))
!         case('h2o')
!            b(1:4,1) = (/  0.06617_dp  ,   0.11045_dp ,&
!                          -0.04915_dp  ,   0.23675_dp  /)               !0th degree
!            b(1:4,2) = (/  55.48e-5_dp ,   0.576e-5_dp,&
!                           70.63e-5_dp , -18.91e-5_dp  /)               !1st degree
!            b(1:4,3) = (/ -48.41e-8_dp ,  24.00e-8_dp ,&
!                          -70.12e-8_dp , -0.907e-8_dp  /)               !2nd degree
!            b(1:4,4) = (/  22.27e-11_dp, -17.01e-11_dp,&
!                           26.07e-11_dp,  4.082e-11_dp /)               !3rd degree
!            b(1:4,5) = (/ -40.17e-15_dp,  30.96e-15_dp,&
!                          -34.94e-15_dp, -8.778e-15_dp /)               !4th degree
!         case('co2')
!            b(1:4,1) = (/  0.09990_dp  ,  0.00942_dp  ,&
!                           0.14511_dp  , -0.02915_dp   /)               !0th degree
!            b(1:4,2) = (/  64.41e-5_dp ,  10.36e-5_dp ,&
!                          -30.73e-5_dp ,  25.23e-5_dp  /)               !1st degree
!            b(1:4,3) = (/ -89.94e-8_dp , -2.277e-8_dp ,&
!                           37.65e-8_dp , -26.10e-8_dp  /)               !2nd degree
!            b(1:4,4) = (/  41.27e-11_dp, -2.134e-11_dp,&
!                          -18.41e-11_dp,  9.965e-11_dp /)               !3rd degree
!            b(1:4,5) = (/ -67.74e-15_dp,  6.497e-15_dp,&
!                           30.16e-15_dp, -13.26e-15_dp /)               !4th degree
!         case('soot')
!            select case(number_wsgg_gases_soot)
!            case(2)
!               b(1:2,1) = (/  0.95552_dp  ,  0.08010_dp   /)            !0th degree
!               b(1:2,2) = (/ -1.431e-3_dp ,  1.290e-3_dp  /)            !1st degree
!               b(1:2,3) = (/  9.871e-7_dp , -7.874e-7_dp  /)            !2nd degree
!               b(1:2,4) = (/ -3.390e-10_dp,  2.322e-10_dp /)            !3rd degree
!               b(1:2,5) = (/  4.555e-14_dp, -3.084e-14_dp /)            !4th degree
!            case(3)
!               b(1:3,1) = (/  0.02812_dp   ,   1.25626_dp ,&
!                             -0.25179_dp  /)                            !0th degree
!               b(1:3,2) = (/ -1.271e-4_dp  , -18.74e-4_dp ,&
!                             18.55e-4_dp  /)                            !1st degree
!               b(1:3,3) = (/  1.395e-7_dp  ,  12.09e-7_dp ,&
!                            -11.39e-7_dp  /)                            !2nd degree
!               b(1:3,4) = (/ -5.672e-11_dp , -39.89e-11_dp,&
!                             34.46e-11_dp /)                            !3rd degree
!               b(1:3,5) = (/  8.0868e-15_dp,  52.85e-15_dp,&
!                            -45.58e-15_dp /)                            !4th degree
!            case(4)
!               b(1:4,1) = (/  0.00129_dp, 1.26110_dp,&
!                             -0.25757_dp, 0.07980_dp  /)                !0th degree
!               b(1:4,2) = (/ -0.545e-5_dp, -319.2e-5_dp,&
!                            362.1e-5_dp, -72.08e-5_dp  /)               !1st degree
!               b(1:4,3) = (/  0.123e-7_dp, 27.72e-7_dp,&
!                            -40.12e-7_dp, 15.87e-7_dp  /)               !2nd degree
!               b(1:4,4) = (/ -0.847e-11_dp, -100.5e-11_dp,&
!                            154.9e-11_dp, -70.89e-11_dp /)              !3rd degree
!               b(1:4,5) = (/  1.6807e-15_dp, 132.8e-15_dp,&
!                           -207.8e-15_dp,  97.69e-15_dp /)              !4th degree
!            case default
!               call shutdown('wsgg_species_a: number_wsgg_gases_soot &
!                             &not available')
!            endselect           
!         case default
!            call shutdown('wsgg_species_a: species_name &
!                          &not available')
!         endselect

!      case('hosein2021')
!         ndegree = 5                                                    !Polynomial degree 
!         tref =  1400._dp                                               !Reference temperature
!         selectcase(trim(species_name))
!         case('soot')
!            b(1,1:6) = (/  1.8613_dp   , -7.7857_dp   ,  1.2809e1_dp, &
!                          -1.0158e1_dp ,  3.8717_dp   , -5.6880e-1_dp /)
!            b(2,1:6) = (/ -1.1374_dp   ,  1.0625e1_dp , -2.1665e1_dp, &
!                           1.8750e1_dp , -7.4578_dp   ,  1.1219_dp    /)
!            b(3,1:6) = (/  2.5975e-1_dp, -2.9708_dp   ,  9.6830_dp,   &
!                          -9.7681_dp   ,  4.1073_dp   , -6.3057e-1_dp /)
!            b(4,1:6) = (/ -3.8367e-2_dp,  4.2159e-1_dp, -1.4101_dp,   &
!                           1.7225_dp   , -7.5889e-1_dp,  1.1454e-1_dp /)
!         case('co')
!            b(1,1:6) = (/ -5.3582e-3_dp, -1.4397e-3_dp,  4.0604e-1_dp,&
!                          -5.7254e-1_dp,  2.8282e-1_dp, -4.7820e-2_dp /)
!            b(2,1:6) = (/ -6.7961e-2_dp,  4.2204e-1_dp, -5.4894e-1_dp,&
!                           2.8819e-1_dp, -6.2318e-2_dp,  3.7321e-3_dp /)
!            b(3,1:6) = (/ -5.7642e-2_dp,  4.2020e-1_dp, -7.6297e-1_dp,&
!                           6.0302e-1_dp, -2.2181e-1_dp, -3.1122e-2_dp /)
!            b(4,1:6) = (/ -1.6152e-2_dp,  1.2220e-1_dp, -2.2207e-1_dp,&
!                           1.7430e-1_dp, -6.3464e-2_dp,  8.8012e-3_dp /)
!         case('ch4')
!            tref =  2500._dp
!            b(1:4,1) = (/ -2.637795e-1_dp, -1.291970e-1_dp,&
!                           1.311677e-3_dp, -2.549570e-1_dp /)
!            b(1:4,2) = (/  3.462304_dp   ,  2.129878_dp   ,&
!                           1.847697e-1_dp,  3.559554_dp   /)
!            b(1:4,3) = (/ -3.845975_dp   , -8.190559_dp   ,&
!                          -1.011002_dp   , -1.100239e1_dp /)
!            b(1:4,4) = (/ -7.555724_dp   ,  1.343738e1_dp ,&
!                           2.020896_dp   ,  1.429368e1_dp /)
!            b(1:4,5) = (/  1.499125e1_dp , -1.023326e1_dp ,&
!                          -1.770474_dp   , -8.453295_dp   /)
!            b(1:4,6) = (/ -6.783667_dp   ,  2.988867_dp   ,&
!                           5.756737e-1_dp,  1.855989_dp   /)
!!            b(1,1:6) = (/ -2.5429e-1_dp, 1.8623_dp   , -1.0442_dp   , -1.4615_dp  ,  1.5196_dp   , -3.7806e-1_dp /)
!!            b(2,1:6) = (/ -2.4021e-1_dp, 1.8795_dp   , -3.2156_dp   , 2.2964_dp   , -7.3711e-1_dp, 8.5605e-2_dp /)
!!            b(3,1:6) = (/ -1.4355e-1_dp, 1.2361_dp   , -2.5390_dp   , 2.2398_dp   , -9.2219e-1_dp, 1.4638e-1_dp /)
!!            b(4,1:6) = (/ -1.2161e-2_dp, 2.7405e-1_dp, -7.3582e-1_dp, 7.7714e-1_dp, -3.6778e-1_dp, 6.5290e-2_dp /)
!         case default
!            call shutdown('wsgg_species_a: species_name not available')
!         endselect

!      case('fraga2021')
!         ndegree = 4                                                    !Polynomial degree 
!         tref =  1000._dp                                               !Reference temperature
!         selectcase(trim(species_name))
!         case('h2o')
!            a_wsgg = fraga_wsgg_a(tmp,mole_frac,jgas,fraga_h2o_interpolation)
!            return
!         case('co2')
!            b(1,1:5) = (/ 1.228740e-1_dp,  5.472120e-1_dp,&
!               -7.693430e-1_dp,  3.698630e-1_dp, -6.083820e-2_dp /)
!            b(2,1:5) = (/ 5.223340e-2_dp, -4.409430e-2_dp,&
!                1.666480e-1_dp, -1.128020e-1_dp,  2.128700e-2_dp /)
!            b(3,1:5) = (/ 1.190670e-1_dp, -2.321120e-1_dp,&
!                2.634100e-1_dp, -1.251010e-1_dp,  2.025260e-2_dp /)
!            b(4,1:5) = (/ 2.285260e-3_dp,  1.334090e-1_dp,&
!               -1.254240e-1_dp,  3.722150e-2_dp, -3.179340e-3_dp /)
!         case('soot')
!            b(1,1:5) = (/  1.450490_dp   , -2.907530_dp,  2.108700_dp,&
!                          -6.586220e-1_dp,  7.539150e-2_dp /)
!            b(2,1:5) = (/ -7.972750e-1_dp,  4.673440_dp, -4.797280_dp,&
!                           1.912160_dp   , -2.740940e-1_dp /)
!            b(3,1:5) = (/  7.397080e-1_dp, -3.790640_dp,  5.861280_dp,&
!                          -3.089560_dp   ,  5.209490e-1_dp /)
!            b(4,1:5) = (/ -5.680900e-1_dp,  2.756410_dp, -4.109730_dp,&
!                           2.300270_dp   , -4.005500e-1_dp /)
!         case('co')
!            b(1,1:5) = (/ -1.893390e-2_dp, 9.334080e-2_dp,&
!                3.768110e-3_dp, -3.503150e-2_dp,  8.796540e-3_dp /)
!            b(2,1:5) = (/ -1.036400e-1_dp, 4.394090e-1_dp,&
!               -4.264890e-1_dp,  1.616730e-1_dp, -2.163050e-2_dp /)
!            b(3,1:5) = (/ -4.956480e-2_dp, 2.621080e-1_dp,&
!               -3.089880e-1_dp,  1.397920e-1_dp, -2.193810e-2_dp /)
!            b(4,1:5) = (/ -1.395890e-2_dp, 7.300760e-2_dp,&
!               -8.310090e-2_dp,  3.628900e-2_dp, -5.522170e-3_dp /)
!         case('ch4')
!            b(1,1:5) = (/ -2.576940e-1_dp, 1.392900_dp   ,&
!               -1.107260_dp   ,  4.189890e-1_dp, -6.624220e-2_dp /)
!            b(2,1:5) = (/ -1.252560e-1_dp, 5.877710e-1_dp,&
!               -2.905960e-1_dp,  1.640850e-2_dp,  7.288570e-3_dp /)
!            b(3,1:5) = (/ -6.488220e-2_dp, 3.589360e-1_dp,&
!               -1.308920e-1_dp, -3.910640e-2_dp,  1.588830e-2_dp /)
!            b(4,1:5) = (/ -1.208290e-2_dp, 2.081460e-1_dp,&
!               -2.559860e-1_dp,  1.081910e-1_dp, -1.542060e-2_dp /)
!         case default
!            call shutdown('wsgg_species_a: species_name not available')
!         endselect

!      case('fraga2021v2')
!         ndegree = 4                                                    !Polynomial degree 
!         tref =  1000._dp                                               !Reference temperature
!         selectcase(trim(species_name))
!         case('co2')
!            b(1,1:5) = (/  2.0693e-1_dp,  5.2037e-1_dp, -8.5223e-1_dp,&
!                           4.3860e-1_dp, -7.5328e-2_dp /)
!            b(2,1:5) = (/ -8.8737e-3_dp,  2.0370e-1_dp, -1.4561e-1_dp,&
!                           3.4399e-2_dp, -2.1191e-3_dp /)
!            b(3,1:5) = (/  1.4750e-1_dp, -3.1930e-1_dp,  4.1047e-1_dp,&
!                          -2.0347e-1_dp,  3.3415e-2_dp /)
!            b(4,1:5) = (/  2.8914e-2_dp,  9.9319e-2_dp, -8.9300e-2_dp,&
!                           1.8716e-2_dp,  6.0005e-5_dp /)
!         case('soot')
!            b(1,1:5) = (/  1.712830_dp   , -3.798550_dp,  3.132310_dp,&
!                          -1.143010_dp   ,  1.563750e-1_dp /)
!            b(2,1:5) = (/ -7.558450e-1_dp,  4.704080_dp, -4.976680_dp,&
!                           2.044010_dp   , -3.015440e-1_dp /)
!            b(3,1:5) = (/  5.246250e-1_dp, -2.755590_dp,  4.346840_dp,& 
!                          -2.236920_dp   ,  3.604750e-1_dp /)
!            b(4,1:5) = (/ -3.471090e-1_dp,  1.672340_dp, -2.483080_dp,&
!                           1.374970_dp   , -2.259010e-1_dp /)
!         case('co')
!            b(1,1:5) = (/ -9.448530e-3_dp, 7.857400e-2_dp,&
!               -1.203580e-2_dp, -1.612550e-2_dp,  4.408760e-3_dp /)
!            b(2,1:5) = (/ -6.492880e-2_dp, 2.374770e-1_dp,&
!               -1.762510e-1_dp,  4.428570e-2_dp, -2.752670e-3_dp /)
!            b(3,1:5) = (/ -8.312190e-2_dp, 4.121860e-1_dp,&
!               -4.697440e-1_dp,  2.072590e-1_dp, -3.189440e-2_dp /)
!            b(4,1:5) = (/ -2.723380e-2_dp, 1.442960e-1_dp,&
!               -1.662020e-1_dp,  7.344460e-2_dp, -1.129700e-2_dp /)
!         case default
!            call shutdown('wsgg_species_a: species_name not available')
!         endselect

!      case default
!         call shutdown('wsgg_species_a: wsgg_species_spec not available')
!      endselect    
      
!      !-----------------------------------------------------------------
!      !Compute the absorption coefficient
!      !-----------------------------------------------------------------
!      ngg = wsgg_species_ngg(wsgg_species_spec,species_name)
!      if (jgas.eq.0) then
!         asum = 0._dp
!         do i=1,ngg
!            asum = asum + wsgg_species_a(tmp,mole_frac,i,&
!                                         wsgg_species_spec,species_name)
!         enddo
!         a_wsgg = max(1._dp - asum,0._dp)
!      else
!         asum = 0._dp
!         do i=1,ndegree+1
!            asum = asum + b(jgas,i)*(tmp/tref)**real(i-1,dp)
!         enddo
!         a_wsgg = max(asum,0._dp)
!      endif
      
!   endfunction wsgg_species_a

!   !====================================================================
!   !Function to compute the absorption coefficient
!   !for my H2O WSGG model
!   !====================================================================
!   real(dp) recursive function fraga_wsgg_kappa(xw,jj,&
!      interpolation_spec) result (kappa_func_res)
   
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: shutdown
!      use math_functions, only: locate
!      implicit none
!      character(*),intent(in) :: interpolation_spec
!      character(8) :: xname(4)
!      character(80) :: msg
!      integer,intent(in) :: jj
!      integer :: i1,i2
!      real(dp),intent(in) :: xw
!      real(dp) :: kp(0:4),xfrac(4),x1,x2,y1,y2
   
!      selectcase(trim(interpolation_spec))
!         case('mr=0.00')
!            kp(0:4) = (/ 0._dp, 1.865631e-1_dp, 1.451070_dp,&
!                         9.828900_dp, 1.157114e2_dp /)
!            kappa_func_res = kp(jj)
!         case('mr=0.05')
!            kp(0:4) = (/ 0._dp, 1.896536e-1_dp, 1.476826_dp,&
!                         9.798522_dp, 1.127594e2_dp /)
!            kappa_func_res = kp(jj)
!         case('mr=0.10')
!            kp(0:4) = (/ 0._dp, 1.919439e-1_dp, 1.498957_dp,&
!                         9.803681_dp, 1.109096e2_dp /)
!            kappa_func_res = kp(jj)
!         case('mr=0.20')
!            kp(0:4) = (/ 0._dp, 1.954878e-1_dp, 1.536312_dp,&
!                         9.853210_dp, 1.085340e2_dp /)
!            kappa_func_res = kp(jj)
!         case('linear')
!            !Array with the available mole fractions
!            xfrac(1:4) = (/ 0._dp, 0.05_dp, 0.1_dp, 0.2_dp /)
!            xname(1:4) = (/ 'mr=0.00', 'mr=0.05', 'mr=0.10', 'mr=0.20'/)
            
!            !Find the upper and lower indexes
!            i1 = min(max(locate(xfrac,xw,4),1),3)
!            i2 = i1 + 1
            
!            !Interpolation bounds
!            x1 = xfrac(i1);   x2 = xfrac(i2)
!            y1 = fraga_wsgg_kappa(x1,jj,xname(i1))
!            y2 = fraga_wsgg_kappa(x2,jj,xname(i2))
         
!            !Interpolate
!            kappa_func_res = y1 + (y2-y1)*(xw-x1)/(x2-x1)
         
!      case default
!         msg = 'fraga_wsgg_kappa: interpolation_spec '//&
!            trim(interpolation_spec)//' not available for model &
!                                       &Fraga et al., 2021'
!         call shutdown(msg)
!      endselect
   
!   endfunction fraga_wsgg_kappa

!   !====================================================================
!   !Function to compute the weighting coefficient
!   !for my H2O WSGG model
!   !====================================================================
!   real(dp) recursive function fraga_wsgg_a(tmp,xw,jj,&
!      interpolation_spec) result (a_func_res)
   
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: shutdown
!      use math_functions, only: locate
!      implicit none
!      character(*),intent(in) :: interpolation_spec
!      character(8) :: xname(4)
!      character(80) :: msg
!      integer,intent(in) :: jj
!      integer :: i,i1,i2
!      real(dp),intent(in) :: tmp,xw
!      real(dp) :: asum,b(1:4,1:5),xfrac(4),x1,x2,y1,y2
   
!      selectcase(trim(interpolation_spec))
!         case('mr=0.00')
!            b(1,1:5) = (/ 1.497890e-1_dp,  2.432210e-1_dp,&
!               -4.425930e-2_dp, -4.850550e-4_dp, -1.108820e-3_dp /)
!            b(2,1:5) = (/ 1.905410e-2_dp,  5.028750e-1_dp,&
!               -3.950850e-1_dp,  1.217500e-1_dp, -1.465740e-2_dp /)
!            b(3,1:5) = (/ 4.667140e-2_dp,  3.727520e-1_dp,&
!               -4.692610e-1_dp,  2.013240e-1_dp, -2.990620e-2_dp /)
!            b(4,1:5) = (/ 1.666230e-1_dp, -2.902130e-1_dp,&
!                2.147880e-1_dp, -7.647290e-2_dp,  1.057280e-2_dp /)
!         case('mr=0.05')
!            b(1,1:5) = (/ 1.513460e-1_dp,  2.393470e-1_dp,&
!               -4.891120e-2_dp,  5.440640e-3_dp, -2.484350e-3_dp /)
!            b(2,1:5) = (/ 3.301640e-2_dp,  4.398470e-1_dp,&
!               -2.974750e-1_dp,  7.085780e-2_dp, -6.060610e-3_dp /)
!            b(3,1:5) = (/ 2.839470e-2_dp,  4.688340e-1_dp,&
!               -5.733440e-1_dp,  2.447280e-1_dp, -3.629550e-2_dp /)
!            b(4,1:5) = (/ 1.838120e-1_dp, -3.249590e-1_dp,&
!                2.411690e-1_dp, -8.554190e-2_dp,  1.176110e-2_dp /)
!         case('mr=0.10')
!            b(1,1:5) = (/ 1.590220e-1_dp,  2.179190e-1_dp,&
!               -3.412370e-2_dp,  2.162570e-3_dp, -2.333450e-3_dp /)
!            b(2,1:5) = (/ 4.463570e-2_dp,  3.852140e-1_dp,&
!               -2.165200e-1_dp,  2.950680e-2_dp,  8.475130e-4_dp /)
!            b(3,1:5) = (/ 9.873750e-3_dp,  5.614650e-1_dp,&
!               -6.756940e-1_dp,  2.882870e-1_dp, -4.280770e-2_dp /)
!            b(4,1:5) = (/ 1.980150e-1_dp, -3.538890e-1_dp,&
!                2.627760e-1_dp, -9.273280e-2_dp,  1.266630e-2_dp /)
!         case('mr=0.20')
!            b(1,1:5) = (/  1.744920e-1_dp,  1.811360e-1_dp,&
!               -1.209960e-2_dp, -1.335980e-3_dp, -2.460420e-3_dp /)
!            b(2,1:5) = (/  7.417510e-2_dp,  2.556630e-1_dp,&
!               -4.162530e-2_dp, -5.664160e-2_dp,  1.502710e-2_dp /)
!            b(3,1:5) = (/ -1.754590e-2_dp,  7.039350e-1_dp,&
!               -8.345180e-1_dp,  3.558690e-1_dp, -5.286660e-2_dp /)
!            b(4,1:5) = (/  2.221620e-1_dp, -4.051810e-1_dp,&
!                3.022710e-1_dp, -1.061560e-1_dp,  1.437860e-2_dp /)
!         case('linear')
!            !Array with the available mole fractions
!            xfrac(1:4) = (/ 0._dp, 0.05_dp, 0.1_dp, 0.2_dp /)
!            xname(1:4) = (/ 'mr=0.00', 'mr=0.05', 'mr=0.10', 'mr=0.20'/)
            
!            !Find the upper and lower indexes
!            i1 = min(max(locate(xfrac,xw,4),1),3)
!            i2 = i1 + 1
            
!            !Interpolation bounds
!            x1 = xfrac(i1);   x2 = xfrac(i2)
!            y1 = fraga_wsgg_a(tmp,x1,jj,xname(i1))
!            y2 = fraga_wsgg_a(tmp,x2,jj,xname(i2))
         
!            !Interpolate
!            a_func_res = y1 + (y2-y1)*(xw-x1)/(x2-x1)
!            return
         
!      case default
!         msg = 'fraga_wsgg_a: interpolation_spec '//&
!            trim(interpolation_spec)//' not available for model &
!                                       &Fraga et al., 2021'
!         call shutdown(msg)
!      endselect
      
!      if (jj.eq.0) then
!         asum = 0._dp
!         do i=1,4
!            asum = asum + fraga_wsgg_a(tmp,xw,i,interpolation_spec)
!         enddo
!         a_func_res = 1._dp - asum
!      else
!         asum = 0._dp
!         do i=1,5
!            asum = asum + b(jj,i)*(tmp/1000._dp)**real(i-1,dp)
!         enddo
!         a_func_res = asum
!      endif
   
!   endfunction fraga_wsgg_a







!   subroutine wsgg_asllanaj2024(kappa_j,a_j,j,tg,xg,pg)
      
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use constants, only: pi
!      use global_parameters, only: id_co2,id_h2o
!      use math_functions, only: linint
!      use precision_parameters, only: dp,small
!      implicit none
!      integer,intent(in) :: j
!      integer :: i
!      real(dp),intent(in) :: pg,tg,xg(:)
!      real(dp),intent(out) :: a_j,kappa_j
!      real(dp) :: a1_array(100),a0_array(6),k0_array(6),T1_array(100)
!      real(dp) :: a1,k1,L_ref,pa,pa_ref
      
!      !-----------------------------------------------------------------
!      !Reference values
!      !-----------------------------------------------------------------
!      pa_ref = 0.3_dp; L_ref = 1._dp
      
!      !-----------------------------------------------------------------
!      !Define arrays
!      !-----------------------------------------------------------------
!      k0_array(1:6) = (/ 0.6382089284_dp, 1.2344705560_dp, &
!                         3.7243386417_dp, 11.7017888278_dp, &
!                        12.1322089821_dp, 99.6564449070_dp /)
!      a0_array(1:6) = (/ 0.2233136596_dp, 0.3253242605_dp, &
!                         0.1556371594_dp, 0.1642889187_dp, &
!                         0.1433675242_dp, 0.1505063464_dp /)
!      a1_array(1:100) = (/ 0.4480865666900_dp, 0.445792521760_dp, &
!         0.442142171980_dp, 0.437230262520_dp, 0.431163124340_dp, &
!         0.424051219320_dp, 0.416007356270_dp, 0.407145731100_dp, &
!         0.397581110710_dp, 0.387427978930_dp, 0.376799568950_dp, &
!         0.365806741110_dp, 0.354556688240_dp, 0.343151482130_dp, &
!         0.331686546990_dp, 0.320249282410_dp, 0.308917824040_dp, &
!         0.297755954130_dp, 0.285998286370_dp, 0.274538324170_dp, &
!         0.263586614800_dp, 0.253192493890_dp, 0.243383028250_dp, &
!         0.234166387780_dp, 0.225533085140_dp, 0.217303355980_dp, &
!         0.209510313430_dp, 0.202360579970_dp, 0.195816457640_dp, &
!         0.189833747630_dp, 0.184361424280_dp, 0.179181838690_dp, &
!         0.174325997610_dp, 0.169905177160_dp, 0.165870965250_dp, &
!         0.162174656480_dp, 0.158744856740_dp, 0.155538590540_dp, &
!         0.152595287280_dp, 0.149875429380_dp, 0.147343909610_dp, &
!         0.144959447010_dp, 0.142723226150_dp, 0.140639583740_dp, &
!         0.138684589550_dp, 0.136835273350_dp, 0.135065015720_dp, &
!         0.133383869400_dp, 0.131780806430_dp, 0.130243772270_dp, &
!         0.128761438380_dp, 0.127339080520_dp, 0.125970400140_dp, &
!         0.124646792640_dp, 0.123361966820_dp, 0.122121591660_dp, &
!         0.120926976490_dp, 0.119771281760_dp, 0.118650619000_dp, &
!         0.117562404140_dp, 0.116514542490_dp, 0.115499878360_dp, &
!         0.114516082290_dp, 0.113561693210_dp, 0.112638045370_dp, &
!         0.111752947950_dp, 0.110895917890_dp, 0.110065768410_dp, &
!         0.109261923480_dp, 0.108484067650_dp, 0.107732100080_dp, &
!         0.107009690870_dp, 0.106316214420_dp, 0.105651163210_dp, &
!         0.105014157310_dp, 0.104405364080_dp, 0.103832217660_dp, &
!         0.103285131330_dp, 0.102763571040_dp, 0.102267170910_dp, &
!         0.101795652970_dp, 0.101348783750_dp, 0.100926342430_dp, &
!         0.100526726470_dp, 0.100151188190_dp, 0.099800653556_dp, &
!         0.099474815087_dp, 0.099173385527_dp, 0.098896102097_dp, &
!         0.098642726314_dp, 0.098413043007_dp, 0.098206859356_dp, &
!         0.098024004086_dp, 0.097864326784_dp, 0.097727697322_dp, &
!         0.097614005323_dp, 0.097523159680_dp, 0.097455088119_dp, &
!         0.097409736809_dp, 0.097387070049_dp /)
!      do i=1,100
!         T1_array(i) = 400._dp + &
!            1400._dp * (dsin(pi*real(i-1,dp)/real(198,dp)))**2
!      enddo

!      !-----------------------------------------------------------------
!      !Main calculation
!      !-----------------------------------------------------------------
!      !Partial pressure of the participating species
!      pa = pg*(xg(id_h2o) + xg(id_co2))
      
!      !Find weighting coefficient of the first gas
!      a1 = linint(T1_array,a1_array,Tg,100)
!      if (Tg.lt.T1_array(1))   a1 = a1_array(1)
!      if (Tg.gt.T1_array(100)) a1 = a1_array(100)
      
!      !Find pressure absorption coefficient of the first gas
!      k1 = k0_array(1) + (a0_array(1) - a1)/(pa_ref*L_ref + small)
      
!      !Find weighting coefficient of gas j
!      a_j = a1 + (a0_array(j) - a0_array(1))
      
!      !Find absorption coefficient of gas j
!      kappa_j = (k1 + (k0_array(j) - k0_array(1)))*pa
      
!   endsubroutine wsgg_asllanaj2024
      
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                      !
!           FUNCTIONS FOR THE ABSORPTION COEFFICIENT                   !
!                                                                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   !====================================================================
   !Function to compute the absorption coefficient for 
   !a single gray gas for fixed-ratio WSGG correlations
   !====================================================================
   real(dp) function fixed_wsgg_kappa(pp,jj,wsgg_spec,species_spec)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      character(*),intent(in) :: wsgg_spec
      character(*),optional :: species_spec
      character(100) :: species_spec_aux
      integer,intent(in) :: jj
      logical :: transparent_window
      real(dp),intent(in) :: pp
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      !Set optional parameter
      species_spec_aux = 'mixture'
      if (present(species_spec)) species_spec_aux = species_spec
      
      !Getting the appropriate WSGG correlation
      call get_wsgg_correlations(wsgg_spec,species_spec_aux)
      
      !Define if this gray gas is a transparent window
      transparent_window = .false.
      if (jj.lt.1) transparent_window = .true.
      if (jj.gt.number_wsgg_gases) transparent_window = .true.
!write(*,*) 'kappa',jj,kappa_p(jj),pp,trim(wsgg_spec),trim(species_spec_aux)
      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      if (transparent_window) then
         fixed_wsgg_kappa = 0._dp
      else
         fixed_wsgg_kappa = kappa_p(jj)*pp
      endif

   endfunction fixed_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a gray gas
   !following different versions of the model of Smith et al (1982)
   !====================================================================
   real(dp) recursive function smith_wsgg_kappa(ppc,ppw,mol_rat,&
      jj,interpolation_spec) result (kappa_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      real(dp),intent(in) :: ppc,ppw,mol_rat
      real(dp) :: kappap_smith,k1(1:4),k2(1:4),k3(1:4),k4(1:4),&
         k5(1:4),x1,x2,y1,y2
   
      !-----------------------------------------------------------------
      !Special case for transparent window
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.3)) then
         kappa_func_res = 0._dp
         return
      endif
   
      !-----------------------------------------------------------------
      !Mounting the ki's arrays
      !-----------------------------------------------------------------
      k1(1:3) = (/ 0.3966_dp, 15.64_dp, 394.3_dp /)                     !pc -> 0 atm
      k2(1:3) = (/ 0.4303_dp, 7.055_dp, 178.1_dp /)                     !p_w/p_c = 1
      k3(1:3) = (/ 0.4201_dp, 6.516_dp, 131.9_dp /)                     !p_w/p_c = 2
      k4(1:3) = (/ 0.4098_dp, 6.325_dp, 120.5_dp /)                     !pw -> 0 atm
      k5(1:3) = (/ 0.4496_dp, 7.113_dp, 119.7_dp /)                     !pw = 1 atm
      
      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      selectcase(trim(interpolation_spec))
         case('none')
            if (mol_rat.eq.0._dp) then
               kappap_smith = k1(jj)
            elseif (mol_rat.eq.1._dp) then
               kappap_smith = k2(jj)
            elseif (mol_rat.eq.2._dp) then
               kappap_smith = k3(jj)
            elseif (mol_rat.eq.smith_infty) then
               if (ppw.le.0.5_dp) then
                  kappap_smith = k4(jj)
               else
                  kappap_smith = k5(jj)
               endif
            else
               write(msg,'(a,f3.1,a)') &
                  'smith_wsgg_kappa: Partial pressure ratio ',&
                  mol_rat,' not available for correlation smith1982'
               call shutdown(msg)               
            endif
            kappa_func_res = (ppc+ppw)*kappap_smith
         case('stepwise')
            if (ppw.le.(0.5_dp*ppc)) then
               kappap_smith = k1(jj)
            elseif (ppw.le.(1.5_dp*ppc)) then
               kappap_smith = k2(jj)
            elseif (ppw.le.(2.5_dp*ppc)) then
               kappap_smith = k3(jj)
            elseif (ppw.le.0.5_dp) then
               kappap_smith = k4(jj)
            else
               kappap_smith = k5(jj)
            endif
            kappa_func_res = (ppc+ppw)*kappap_smith
         case('linear')
            if (mol_rat.le.1._dp) then
               x1 = 0._dp; x2 = 1._dp
               y1 = smith_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = smith_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2._dp) then
               x1 = 1._dp; x2 = 2._dp
               y1 = smith_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = smith_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.smith_infty) then
               x1 = 2._dp; x2 = smith_infty
               y1 = smith_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = smith_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               kappa_func_res = smith_wsgg_kappa(ppc,ppw,&
                                                 smith_infty,jj,'none')
            endif
      endselect
   
   endfunction smith_wsgg_kappa
   
   !====================================================================
   !Function to compute the absorption coefficient for a 
   !single gray gas following the formulation of Yin (2010) 
   !====================================================================
   real(dp) recursive function yin2010_wsgg_kappa(ppc,ppw,mol_rat,jj,&
      interpolation_spec) result (kappa_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      real(dp),intent(in) :: ppc,ppw,mol_rat
      real(dp) :: kappap_yin,k1(1:4),k2(1:4),k3(1:4),k4(1:4),k5(1:4),&
         k6(1:4),k7(1:4),k8(1:4),k9(1:4),k10(1:4)
      real(dp) :: x1,x2,y1,y2

      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         kappa_func_res = 0._dp
         return
      endif

      !-----------------------------------------------------------------
      !Mounting the ki's arrays
      !-----------------------------------------------------------------
      k1 (1:4) = (/ 0.009422_dp, 0.415646_dp, 11.617018_dp, 319.911168_dp /)   !pw -> 0, pc -> 0
      k2 (1:4) = (/ 0.256738_dp, 3.108033_dp, 52.585782_dp, 440.845718_dp /)   !pw = 0.1 atm, pc = 0.1 atm
      k3 (1:4) = (/ 0.132242_dp, 14.660767_dp, 1.750654_dp, 165.763926_dp /)   !pw = 0.3 atm, pc = 0.1 atm
      k4 (1:4) = (/ 0.051237_dp, 0.688383_dp, 13.763205_dp, 289.841885_dp /)   !pw/pc = 1/8
      k5 (1:4) = (/ 0.052694_dp, 0.752776_dp, 11.543306_dp, 252.938841_dp /)   !pw/pc = 1/4
      k6 (1:4) = (/ 0.052378_dp, 0.712283_dp, 8.067637_dp, 195.892573_dp /)    !pw/pc = 1/2
      k7 (1:4) = (/ 0.051639_dp, 0.617739_dp, 6.051770_dp, 150.875915_dp /)    !pw/pc = 3/4
      k8 (1:4) = (/ 0.051487_dp, 0.571797_dp, 5.398936_dp, 130.622859_dp /)    !pw/pc = 1
      k9 (1:4) = (/ 0.054480_dp, 0.555304_dp, 5.040174_dp, 100.372663_dp /)    !pw/pc = 2
      k10(1:4) = (/ 0.060800_dp, 5.608831_dp, 0.676040_dp, 84.540632_dp /)     !pw/pc = 4
      
      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      selectcase(trim(interpolation_spec))
         case('none')
            if (mol_rat.eq.0._dp) then
               kappap_yin = k1(jj)
            elseif (mol_rat.eq.0.125_dp) then
               kappap_yin = k4(jj)
            elseif (mol_rat.eq.0.25_dp) then
               kappap_yin = k5(jj)
            elseif (mol_rat.eq.0.5_dp) then
               kappap_yin = k6(jj)
            elseif (mol_rat.eq.0.75_dp) then
               kappap_yin = k7(jj)
            elseif (mol_rat.eq.1._dp) then
               kappap_yin = k8(jj)
            elseif (mol_rat.eq.2._dp) then
               kappap_yin = k9(jj)
            elseif (mol_rat.eq.4._dp) then
               kappap_yin = k10(jj)
            else
               write(msg,'(a,f3.1,a)') &
                  'yin2010_wsgg_kappa: Partial pressure ratio ',&
                  mol_rat,' not available for correlation yin2010'
               call shutdown(msg)               
            endif
            kappa_func_res = (ppc+ppw)*kappap_yin
               
         case('stepwise')
            if ((ppc+ppw).le.0.1_dp) then
               kappap_yin = k1(jj)
            elseif ((ppc+ppw).le.0.3_dp) then
               kappap_yin = k2(jj)
            elseif ((ppc+ppw).le.0.5_dp) then
               kappap_yin = k3(jj)
            elseif (ppw.le.(0.2_dp*ppc)) then            
               kappap_yin = k4(jj)
            elseif (ppw.le.(0.4_dp*ppc)) then
               kappap_yin = k5(jj)
            elseif (ppw.le.(0.6_dp*ppc)) then
               kappap_yin = k6(jj)
            elseif (ppw.le.(0.9_dp*ppc)) then
               kappap_yin = k7(jj)
            elseif (ppw.le.(1.1_dp*ppc)) then
               kappap_yin = k8(jj)
            elseif (ppw.le.(2.5_dp*ppc)) then
               kappap_yin = k9(jj)
            else
               kappap_yin = k10(jj)
            endif
            kappa_func_res = (ppc+ppw)*kappap_yin
               
         case('linear')
            if (mol_rat.le.0.125_dp) then
               kappa_func_res = yin2010_wsgg_kappa(ppc,ppw,0.125_dp,&
                                                   jj,'none')
            elseif (mol_rat.le.0.250_dp) then
               x1 = 0.125_dp; x2 = 0.250_dp
               y1 = yin2010_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = yin2010_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.500_dp) then
               x1 = 0.250_dp; x2 = 0.500_dp
               y1 = yin2010_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = yin2010_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.750_dp) then
               x1 = 0.500_dp; x2 = 0.750_dp
               y1 = yin2010_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = yin2010_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1.000_dp) then
               x1 = 0.750_dp; x2 = 1.000_dp
               y1 = yin2010_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = yin2010_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2.000_dp) then
               x1 = 1.000_dp; x2 = 2.000_dp
               y1 = yin2010_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = yin2010_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.4.000_dp) then
               x1 = 2.000_dp; x2 = 4.000_dp
               y1 = yin2010_wsgg_kappa(ppc,ppw,x1,jj,'none')
               y2 = yin2010_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + &
                  (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               kappa_func_res = yin2010_wsgg_kappa(ppc,ppw,4._dp,&
                                                   jj,'none')
            endif
      endselect   
   
   endfunction yin2010_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a 
   !single gray gas following the formulation of 
   !Kangwanpongpan et al (2012) for a variable molar ratio
   !====================================================================
   real(dp) function kangwanpongpan_variable_wsgg_kappa(pp,mol_rat,jj)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer,intent(in) :: jj
      real(dp),intent(in) :: mol_rat,pp
      real(dp) :: ck1(1:4),ck2(1:4),ck3(1:4)

      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         kangwanpongpan_variable_wsgg_kappa = 0._dp
         return
      endif

      !-----------------------------------------------------------------
      !Compute the absorption coefficient
      !-----------------------------------------------------------------
      !Mounting the ck1, ck2 and ck3 arrays
      ck1(1:4) = (/  0.0429_dp,  0.3647_dp,  3.7144_dp,  105.31_dp /)
      ck2(1:4) = (/  0.0093_dp,  0.0790_dp,  0.2565_dp, -39.265_dp /)
      ck3(1:4) = (/ -0.0018_dp, -0.0150_dp, -0.0509_dp,  6.0877_dp /)

      !Computing the absorption coefficient
      kangwanpongpan_variable_wsgg_kappa = &
         pp*( ck1(jj) + ck2(jj)*mol_rat + ck3(jj)*(mol_rat**2._dp) )

   endfunction kangwanpongpan_variable_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a single
   !gray gas following different versions of the formulation of 
   !Kangwanpongpan et al (2012) for a fixed molar ratio 
   !====================================================================
   real(dp) recursive function kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,&
      mol_rat,jj,interpolation_spec) result (kappa_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      real(dp),intent(in) :: ppc,ppw,mol_rat
      real(dp) :: kappap_kang,k1(1:4),k2(1:4),k3(1:4),k4(1:4),k5(1:4),&
         k6(1:4),k7(1:4)
      real(dp) :: x1,x2,y1,y2

      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         kappa_func_res = 0._dp
         return
      endif
      
      !-----------------------------------------------------------------
      !Mounting the ki's arrays
      !-----------------------------------------------------------------
      k1(1:4) = (/ 0.049191_dp, 0.343891_dp, 3.710740_dp, 106.080000_dp/)   !MR = 0.125
      k2(1:4) = (/ 0.045079_dp, 0.386642_dp, 3.764160_dp, 96.034300_dp /)   !MR = 0.25
      k3(1:4) = (/ 0.049191_dp, 0.421502_dp, 3.852390_dp, 83.534000_dp /)   !MR = 0.5
      k4(1:4) = (/ 0.050784_dp, 0.431878_dp, 3.908780_dp, 75.525500_dp /)   !MR = 0.75
      k5(1:4) = (/ 0.051446_dp, 0.436145_dp, 3.948270_dp, 69.781000_dp /)   !MR = 1.0
      k6(1:4) = (/ 0.051832_dp, 0.440593_dp, 3.981020_dp, 56.081800_dp /)   !MR = 2.0
      k7(1:4) = (/ 0.051687_dp, 0.444449_dp, 3.933740_dp, 44.750100_dp /)   !MR = 4.0

      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      selectcase(trim(interpolation_spec))
         case('none')
            if (mol_rat.eq.0.125_dp) then
               kappap_kang = k1(jj)
            elseif (mol_rat.eq.0.25_dp) then
               kappap_kang = k2(jj)
            elseif (mol_rat.eq.0.5_dp) then
               kappap_kang = k3(jj)
            elseif (mol_rat.eq.0.75_dp) then
               kappap_kang = k4(jj)
            elseif (mol_rat.eq.1._dp) then
               kappap_kang = k5(jj)
            elseif (mol_rat.eq.2._dp) then
               kappap_kang = k6(jj)
            elseif (mol_rat.eq.4._dp) then
               kappap_kang = k7(jj)
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',&
                  mol_rat,' not available for correlation&
                  & kangwanpongpan2012_fixed'
               call shutdown(msg)               
            endif
            kappa_func_res = (ppc+ppw)*kappap_kang
               
         case('stepwise')
            if (ppw.le.(0.1875_dp*ppc)) then            
               kappap_kang = k1(jj)
            elseif (ppw.le.(0.375_dp*ppc)) then
               kappap_kang = k2(jj)
            elseif (ppw.le.(0.625_dp*ppc)) then
               kappap_kang = k3(jj)
            elseif (ppw.le.(0.875_dp*ppc)) then
               kappap_kang = k4(jj)
            elseif (ppw.le.(1.5_dp*ppc)) then
               kappap_kang = k5(jj)
            elseif (ppw.le.(3._dp*ppc)) then
               kappap_kang = k6(jj)
            else
               kappap_kang = k7(jj)
            endif
            kappa_func_res = (ppc+ppw)*kappap_kang
               
         case('linear')
            if (mol_rat.le.0.125_dp) then
               kappa_func_res = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,0.125_dp,jj,'none')
            elseif (mol_rat.le.0.250_dp) then
               x1 = 0.125_dp  ; y1 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 0.250_dp  ; y2 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.500_dp) then
               x1 = 0.250_dp  ; y1 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 0.500_dp  ; y2 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.750_dp) then
               x1 = 0.500_dp  ; y1 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 0.750_dp  ; y2 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1.000_dp) then
               x1 = 0.750_dp  ; y1 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 1.000_dp  ; y2 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2.000_dp) then
               x1 = 1.000_dp  ; y1 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 2.000_dp  ; y2 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.4.000_dp) then
               x1 = 2.000_dp  ; y1 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 4.000_dp  ; y2 = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               kappa_func_res = kangwanpongpan_fixed_wsgg_kappa(ppc,ppw,4._dp,jj,'none')
            endif
      endselect   
   
   endfunction kangwanpongpan_fixed_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a single
   !gray gas following the formulation of Johansson et al (2011) 
   !====================================================================
   real(dp) function johansson_wsgg_kappa(pp,mol_rat,jj)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      real(dp),intent(in) :: mol_rat,pp
      real(dp) :: k1(1:4),k2(1:4)

      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         johansson_wsgg_kappa = 0._dp
         return
      endif
      
      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      !Mounting the k's array
      k1(1:4) = (/ 0.055_dp,  0.880_dp,  10._dp, 135._dp /)
      k2(1:4) = (/ 0.012_dp, -0.021_dp, -1.6_dp, -35._dp /)
      
      !Absorption coefficient
      johansson_wsgg_kappa = pp*( k1(jj) + k2(jj)*mol_rat )
      
   endfunction johansson_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a single
   !gray gas following the formulation of Krishnamoorthy (2013) 
   !====================================================================
   real(dp) recursive function krishnamoorthy2013_wsgg_kappa(ppc,ppw,&
      mol_rat,jj,interpolation_spec) result (kappa_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      real(dp),intent(in) :: ppc,ppw,mol_rat
      real(dp) :: kappap_krishna,k1(1:4),k2(1:4),k3(1:4),k4(1:4),&
         x1,x2,y1,y2

      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         kappa_func_res = 0._dp
         return
      endif
      
      !-----------------------------------------------------------------
      !Mounting the ki's arrays
      !-----------------------------------------------------------------
      k1(1:4) = (/ 0.06592_dp, 0.99698_dp, 10.00038_dp, 100.00000_dp /)   !MR = 0.11
      k2(1:4) = (/ 0.10411_dp, 1.00018_dp, 9.99994_dp,  100.00000_dp /)   !MR = 0.5
      k3(1:4) = (/ 0.20616_dp, 1.39587_dp, 8.56904_dp,  99.75698_dp /)   !MR = 1.0
      k4(1:4) = (/ 0.21051_dp, 1.33782_dp, 8.55495_dp,  99.75649_dp /)   !MR = 2.0
      
      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      selectcase(trim(interpolation_spec))
         case('none')
            if (mol_rat.eq.0.11_dp) then
               kappap_krishna = k1(jj)
            elseif (mol_rat.eq.0.5_dp) then
               kappap_krishna = k2(jj)
            elseif (mol_rat.eq.1._dp) then
               kappap_krishna = k3(jj)
            elseif (mol_rat.eq.2._dp) then
               kappap_krishna = k4(jj)
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation krishnamoorthy2013'
               call shutdown(msg)               
            endif
            kappa_func_res = (ppc+ppw)*kappap_krishna
         case('stepwise')
            if (ppw.le.(0.2_dp*ppc)) then
               kappap_krishna = k1(jj)
            elseif (ppw.le.(0.67_dp*ppc)) then
               kappap_krishna = k2(jj)
            elseif (ppw.le.(1.5_dp*ppc)) then
               kappap_krishna = k3(jj)
            else
               kappap_krishna = k4(jj)
            endif
            kappa_func_res = (ppc+ppw)*kappap_krishna
            
         case('linear')
            if (mol_rat.le.0.11_dp) then
               kappa_func_res = krishnamoorthy2013_wsgg_kappa(ppc,ppw,0.11_dp,jj,'none')
            elseif (mol_rat.le.0.5_dp) then
               x1 = 0.11_dp  ; y1 = krishnamoorthy2013_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 0.50_dp  ; y2 = krishnamoorthy2013_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1._dp) then
               x1 = 0.5_dp  ; y1 = krishnamoorthy2013_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 1.0_dp  ; y2 = krishnamoorthy2013_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2._dp) then
               x1 = 1._dp  ; y1 = krishnamoorthy2013_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 2._dp  ; y2 = krishnamoorthy2013_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               kappa_func_res = krishnamoorthy2013_wsgg_kappa(ppc,ppw,2._dp,jj,'none')
            endif
      endselect   
      
   endfunction krishnamoorthy2013_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a single
   !gray gas following the formulation of Yin (2013) 
   !====================================================================
   real(dp) recursive function yin2013_wsgg_kappa(ppc,ppw,mol_rat,&
      jj,interpolation_spec) result (kappa_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      real(dp),intent(in) :: ppc,ppw,mol_rat
      real(dp) :: kappap_yin,k1(1:4),k2(1:4),k3(1:4),k4(1:4),k5(1:4),&
         k6(1:4),k7(1:4),x1,x2,y1,y2

      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         kappa_func_res = 0._dp
         return
      endif
      
      !-----------------------------------------------------------------
      !Mounting the ki's arrays
      !-----------------------------------------------------------------
      k1(1:4) = (/ 0.163233_dp, 13.096584_dp, 175.474735_dp, 1310.847307_dp /) !pc -> 0 atm
      k2(1:4) = (/ 0.352505_dp, 8.210621_dp, 137.410012_dp, 1269.710976_dp /)  !p_w/p_c = 0.05
      k3(1:4) = (/ 0.261021_dp, 3.147817_dp, 54.265868_dp, 482.900353_dp /)    !p_w/p_c = 1
      k4(1:4) = (/ 0.179160_dp, 2.388971_dp, 28.415805_dp, 253.059089_dp /)    !p_w/p_c = 2
      k5(1:4) = (/ 0.085523_dp, 0.475777_dp, 8.549733_dp, 201.906503_dp /)     !pw -> 0 atm
      k6(1:4) = (/ 0.232724_dp, 2.134299_dp, 9.266065_dp, 134.988332_dp /)     !pw = 0.05 atm
      k7(1:4) = (/ 0.065411_dp, 0.696552_dp, 4.862610_dp, 60.255980_dp /)      !pw = 1 atm
      
      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      selectcase(trim(interpolation_spec))
         case('none')
            if (mol_rat.eq.0._dp) then
               kappap_yin = k1(jj)
            elseif (mol_rat.eq.0.05_dp) then
               kappap_yin = k2(jj)
            elseif (mol_rat.eq.1._dp) then
               kappap_yin = k3(jj)
            elseif (mol_rat.eq.2._dp) then
               kappap_yin = k4(jj)
            elseif (mol_rat.eq.yin_infty) then
               if (ppw.le.0.01) then
                  kappap_yin = k5(jj)
               elseif (ppw.le.0.2) then
                  kappap_yin = k6(jj)
               else
                  kappap_yin = k7(jj)
               endif
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation yin2013'
               call shutdown(msg)               
            endif
            kappa_func_res = (ppc+ppw)*kappap_yin
         case('stepwise')
            if (ppw.le.(0.01_dp*ppc)) then
               kappap_yin = k1(jj)
            elseif (ppw.le.(0.5_dp*ppc)) then
               kappap_yin = k2(jj)
            elseif (ppw.le.(1.5_dp*ppc)) then
               kappap_yin = k3(jj)
            elseif (ppw.le.(2.5_dp*ppc)) then
               kappap_yin = k4(jj)
            elseif (ppw.le.0.01_dp) then
               kappap_yin = k5(jj)
            elseif (ppw.le.0.2_dp) then
               kappap_yin = k6(jj)
            else
               kappap_yin = k7(jj)
            endif
            kappa_func_res = (ppc+ppw)*kappap_yin
            
         case('linear')
            if (mol_rat.le.0.05_dp) then
               x1 = 0.00_dp  ; y1 = yin2013_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 0.05_dp  ; y2 = yin2013_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1._dp) then
               x1 = 0.05_dp  ; y1 = yin2013_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 1.00_dp  ; y2 = yin2013_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2._dp) then
               x1 = 1._dp  ; y1 = yin2013_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = 2._dp  ; y2 = yin2013_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.yin_infty) then
               x1 = 2._dp     ; y1 = yin2013_wsgg_kappa(ppc,ppw,x1,jj,'none')
               x2 = yin_infty ; y2 = yin2013_wsgg_kappa(ppc,ppw,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               kappa_func_res = yin2013_wsgg_kappa(ppc,ppw,yin_infty,jj,'none')
            endif
      endselect
      
   endfunction yin2013_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a single
   !gray gas following the formulation of Bordbar et al (2014)
   !====================================================================
   real(dp) recursive function bordbar_wsgg_kappa(pp,mol_rat,jj) &
      result(kappa_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: mol_rat,pp
      real(dp) :: d(1:4,0:4),kappap(4),sum_kappa,x1,x2,y1,y2,mr

      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         kappa_func_res = 0._dp
         return
      endif
      
      !-----------------------------------------------------------------
      !Compute the absorption coefficient 
      !within five intervals of molar ratio
      !-----------------------------------------------------------------
      if ((mol_rat.eq.0._dp).and.bordbar_co2_interpolation) then        !Only CO2
         !Mounting the kappa_p array
         kappap(1:4) = (/ 3.388079e-2_dp, 4.544269e-1_dp,&
                          4.680226_dp, 1.038439e2_dp/)
         
         !Getting the absorption coefficient
         kappa_func_res = kappap(jj)*pp
         
      elseif ((mol_rat.ge.bordbar_infty)&
              .and.bordbar_h2o_interpolation) then                      !Only H2O (virtually)
         !Mounting the kappa_p array
         kappap(1:4) = (/ 7.703541e-2_dp, 8.242941e-1_dp,&
                          6.854761_dp, 6.593653e1_dp /)
         
         !Getting the absorption coefficient
         kappa_func_res = kappap(jj)*pp
         
      elseif ((mol_rat.lt.0.01_dp).and.bordbar_co2_interpolation) then  !0<Molar ratio<0.01
         !Bounds for the interpolation
         x1 = 0._dp  ; y1 = bordbar_wsgg_kappa(pp,x1,jj)
         x2 = 0.01_dp; y2 = bordbar_wsgg_kappa(pp,x2,jj)
         
         !Interpolating
         kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
         
      elseif ((mol_rat.gt.4._dp).and.bordbar_h2o_interpolation) then    !4<Molar ratio<bordbar_infty
         !Values for the interpolation
         x1 = 4._dp        ; y1 = bordbar_wsgg_kappa(pp,x1,jj)
         x2 = bordbar_infty; y2 = bordbar_wsgg_kappa(pp,x2,jj)
         
         !Interpolating
         kappa_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
         
      else                                                              !Normal CO2-H2O mixture
         mr = min(max(mol_rat,0.01_dp),4._dp)
         
         !Mounting the d's array
         d(1,0:4) = (/ 0.0340429_dp,  0.0652305_dp, -0.0463685_dp,&
                       0.0138684_dp, -0.0014450_dp /)
         d(2,0:4) = (/ 0.3509457_dp,  0.7465138_dp, -0.5293090_dp,&
                       0.1594423_dp, -0.0166326_dp /)
         d(3,0:4) = (/ 4.5707400_dp,  2.1680670_dp, -1.4989010_dp,&
                       0.4917165_dp, -0.0542999_dp /)
         d(4,0:4) = (/ 109.81690_dp, -50.923590_dp,  23.432360_dp,&
                      -5.1638920_dp,  0.4393889_dp /)
      
         !Computing the absorption coefficient
         sum_kappa = 0._dp
         do ii = 0,4
            sum_kappa = sum_kappa + d(jj,ii)*(mr**real(ii,dp))
         enddo
         kappa_func_res = pp*sum_kappa
      endif
      
   endfunction bordbar_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a single
   !gray gas following the formulation of Guo et al (2015)
   !====================================================================
   real(dp) function guo_wsgg_kappa(ttmp,pp,mol_rat,jj)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: mol_rat,pp,ttmp
      real(dp) :: c1(1:4,1:4),c2(1:4,1:4),c3(1:4,1:4),ck(1:4),&
                  sum_kappa,ttmp_ref

      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         guo_wsgg_kappa = 0._dp
         return
      endif

      !-----------------------------------------------------------------
      !Mounting the c1, c2 and c3 arrays
      !-----------------------------------------------------------------
      c1(1,1:4) = (/ -1.98431e-4_dp, -4.77095e-5_dp, -8.75575e-4_dp,  4.00925e-4_dp /)
      c1(2,1:4) = (/  2.59122e-3_dp, -1.10845e-2_dp, -1.09658e-3_dp, -1.85170e-3_dp /)
      c1(3,1:4) = (/  2.73794e-2_dp, -7.15296e-2_dp, -8.28257e-2_dp, -3.16580e-2_dp /)
      c1(4,1:4) = (/  0.373884_dp,    -3.16546_dp,     8.23483_dp,     -7.50258_dp   /)
      
      c2(1,1:4) = (/  2.10006e-4_dp, 3.03422e-3_dp, -9.03821e-4_dp, 2.83204e-4_dp /)
      c2(2,1:4) = (/ -1.28988e-2_dp, 5.54857e-2_dp, -9.26126e-3_dp, 7.78285e-3_dp /)
      c2(3,1:4) = (/ -6.77349e-2_dp, 9.74045e-2_dp,  0.545000_dp,   7.57066e-2_dp /)
      c2(4,1:4) = (/ -1.81510_dp,    12.9141_dp,    -30.7345_dp,    26.8828_dp    /)

      c3(1,1:4) = (/ -8.77357e-5_dp,  5.80970e-4_dp, -3.67320e-4_dp, 3.33187e-4_dp /)
      c3(2,1:4) = (/ -2.07199e-3_dp,  1.12076e-2_dp, -1.53483e-2_dp, 1.37598e-2_dp /)
      c3(3,1:4) = (/ -1.10101e-2_dp,  2.86980e-2_dp,  1.52726e-2_dp, 6.34926e-2_dp /)
      c3(4,1:4) = (/  0.729046_dp,   -4.53601_dp,       8.07119_dp,    3.17932_dp      /)

      !-----------------------------------------------------------------
      !Mounting the ck array
      !-----------------------------------------------------------------
      do ii=1,4
         ck(ii) = c1(jj,ii)*(mol_rat**2._dp) + c2(jj,ii)*mol_rat + c3(jj,ii)
      enddo

      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      ttmp_ref = 1500._dp; sum_kappa = 0._dp
      do ii=1,4
         sum_kappa = sum_kappa + ck(ii)*((ttmp/ttmp_ref)**real(4-ii,dp))
      enddo
      guo_wsgg_kappa = sum_kappa*pp

   endfunction guo_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a single
   !gray gas following the formulation of Ge et al (2017)
   !====================================================================
   real(dp) function ge_wsgg_kappa(ttmp,pp,mol_rat,jj)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: mol_rat,pp,ttmp
      real(dp) :: b_aux,b1(1:4,1:4),b2(1:4,1:4),b3(1:4,1:4),&
         sum_kappa,ttmp_ref
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         ge_wsgg_kappa = 0._dp
         return
      endif
      
      !-----------------------------------------------------------------
      !Mounting the bk's arrays
      !-----------------------------------------------------------------
      b1(1,1:4) = (/ 105.3082_dp, 0.017262_dp, -0.01525_dp, 0.003911_dp /)
      b1(2,1:4) = (/ 3.386233_dp, 1.693028_dp, -1.19358_dp, 0.266963_dp /)
      b1(3,1:4) = (/ -0.06371_dp, 1.173845_dp, -0.95320_dp, 0.232854_dp /)
      b1(4,1:4) = (/ 0.007763_dp, 0.028807_dp, -0.02902_dp, 0.007428_dp /)
      
      b2(1,1:4) = (/ -39.1954_dp, -0.36690_dp, 0.280155_dp, -0.06612_dp /)
      b2(2,1:4) = (/ -0.83532_dp, 4.250212_dp, -3.33618_dp, 0.799926_dp /)
      b2(3,1:4) = (/ 1.191368_dp, -4.86797_dp, 4.197609_dp, -1.05866_dp /)
      b2(4,1:4) = (/ 0.013514_dp, -0.02448_dp, 0.026231_dp, -0.00774_dp /)
      
      b3(1,1:4) = (/ 6.043568_dp, 0.237710_dp, -0.17932_dp, 0.041997_dp /)
      b3(2,1:4) = (/ 0.678994_dp, -3.00472_dp, 2.296622_dp, -0.54187_dp /)
      b3(3,1:4) = (/ -0.54334_dp, 2.702553_dp, -2.27965_dp, 0.568522_dp /)
      b3(4,1:4) = (/ -0.00993_dp, 0.030806_dp, -0.02717_dp, 0.006980_dp /)
         
      !-----------------------------------------------------------------
      !Computing the absorption coefficient
      !-----------------------------------------------------------------
      ttmp_ref = 1200._dp; sum_kappa = 0._dp
      do ii=1,4
         b_aux = b1(jj,ii) + b2(jj,ii)*mol_rat + b3(jj,ii)*mol_rat**2._dp
         sum_kappa = sum_kappa + b_aux*(ttmp/ttmp_ref)**real((ii-1),dp)
      enddo
      ge_wsgg_kappa = sum_kappa*pp
      
   endfunction ge_wsgg_kappa

   !====================================================================
   !Function to compute the absorption coefficient for a single
   !gray gas following the formulation of Shan et al (2018)
   !====================================================================
   real(dp) recursive function shan_wsgg_kappa(pp,tp,mol_rat,&
      jj,interpolation_spec) result(kappa_func_res)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      real(dp),intent(in) :: pp,tp,mol_rat
      real(dp) :: k0(1:4),k1(1:4),k2(1:4),&
         x1,x2,y1,y2
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         kappa_func_res = 0._dp
         return
      endif

      selectcase(trim(interpolation_spec))
         case('none')
            !-----------------------------------------------------------
            !Mounting the k's arrays
            !-----------------------------------------------------------
            if (mol_rat.eq.0.125_dp) then
               if (tp.le.2._dp) then
                  k0(1:4) = (/ -2.9167910000e-3_dp, -1.3571366200e-1_dp, -5.1576318000e-2_dp,  5.8876288600e+1_dp /)
                  k1(1:4) = (/  5.8045024000e-2_dp,  6.7388017200e-1_dp,  5.5070752990e+0_dp,  9.5668485800e+1_dp /)
                  k2(1:4) = (/ -1.5377880000e-2_dp, -1.9086652000e-1_dp, -1.5181246220e+0_dp, -2.5110022000e+1_dp /)
               elseif (tp.le.4._dp) then
                  k0(1:4) = (/ -1.6409414000e-2_dp, -9.7521380000e-3_dp,  3.0840317430e+0_dp,  1.3859741880e+2_dp /)
                  k1(1:4) = (/  4.7022294500e-2_dp,  3.1607257600e-1_dp,  1.2643348875e+0_dp,  9.6257067000e+0_dp /)
                  k2(1:4) = (/ -6.4933595000e-3_dp, -4.3453103000e-2_dp, -1.8065643150e-1_dp, -2.0189150000e+0_dp /)
               elseif (tp.le.6._dp) then
                  k0(1:4) = (/  5.8049010000e-2_dp,  5.0919353800e-1_dp,  5.3961325110e+0_dp,  1.6383510940e+2_dp /)
                  k1(1:4) = (/  4.3572685000e-3_dp,  1.7275237000e-2_dp, -6.1247200500e-2_dp, -5.5885669500e+0_dp /)
                  k2(1:4) = (/ -4.8075450000e-4_dp, -1.1878730000e-3_dp,  6.2327925000e-3_dp,  2.0729775000e-1_dp /)
               elseif (tp.le.10._dp) then
                  k0(1:4) = (/  7.9318545000e-2_dp,  5.7884122600e-1_dp,  5.2810289130e+0_dp,  1.5172829500e+2_dp /)
                  k1(1:4) = (/ -2.5535497500e-3_dp, -6.0209300000e-3_dp, -2.7338952500e-2_dp, -3.3358517250e+0_dp /)
                  k2(1:4) = (/  8.0228125000e-5_dp,  7.6016350000e-4_dp,  3.7787400000e-3_dp,  1.6814561250e-1_dp /)
               elseif (tp.le.20._dp) then
                  k0(1:4) = (/  3.0904496000e-2_dp,  7.9552120000e-2_dp,  3.2230159440e+0_dp,  1.3058027350e+2_dp /)
                  k1(1:4) = (/  4.3410448000e-3_dp,  6.8346068800e-2_dp,  2.8874469640e-1_dp,  3.3528269000e-1_dp /)
                  k2(1:4) = (/ -1.2509084000e-4_dp, -1.6836453200e-3_dp, -7.2494952000e-3_dp,  1.2512386000e-2_dp /)
               elseif (tp.le.30._dp) then
                  k0(1:4) = (/  6.7375840000e-2_dp,  5.1426632200e-1_dp,  5.1653596100e+0_dp,  1.2378543430e+2_dp /)
                  k1(1:4) = (/ -1.6779280000e-4_dp,  1.7448804700e-2_dp,  6.6467331900e-2_dp,  1.0250152900e+0_dp /)
                  k2(1:4) = (/  9.1726800000e-6_dp, -2.2556762000e-4_dp, -9.9148614000e-4_dp, -4.9871460000e-3_dp /)
               endif
            elseif (mol_rat.eq.0.25_dp) then
               if (tp.le.2._dp) then
                  k0(1:4) = (/  1.9706494000e-2_dp,  1.3494304200e-1_dp,  1.4139805140e+0_dp,  6.5741511300e+1_dp /)
                  k1(1:4) = (/  2.1098919000e-2_dp,  2.6016262800e-1_dp,  3.1658765890e+0_dp,  7.0810186300e+1_dp /)
                  k2(1:4) = (/  2.6014200000e-4_dp, -3.9464140000e-2_dp, -6.7431656200e-1_dp, -1.7895879800e+1_dp /)
               elseif (tp.le.4._dp) then
                  k0(1:4) = (/  9.4691834000e-2_dp,  5.4128123100e-1_dp,  4.4095695510e+0_dp,  1.3456946630e+2_dp /)
                  k1(1:4) = (/ -2.1398633000e-2_dp, -3.6425581500e-2_dp,  4.0497079950e-1_dp,  2.2604134000e+0_dp /)
                  k2(1:4) = (/  2.7625830000e-3_dp,  7.2454175000e-3_dp, -4.2760926500e-2_dp, -8.2798210000e-1_dp /)
               elseif (tp.le.6._dp) then
                  k0(1:4) = (/  6.4874500000e-2_dp,  4.5261486900e-1_dp,  4.8774798090e+0_dp,  1.5163382130e+2_dp /)
                  k1(1:4) = (/ -3.9824215000e-3_dp,  1.3514571000e-2_dp,  1.5755131300e-1_dp, -7.3798539500e+0_dp /)
                  k2(1:4) = (/  2.7211350000e-4_dp,  3.0202700000e-4_dp, -1.0150446000e-2_dp,  5.1556255000e-1_dp /)
               elseif (tp.le.10._dp) then
                  k0(1:4) = (/  3.2783110000e-3_dp,  4.8423850500e-1_dp,  6.1060145370e+0_dp,  1.4294459550e+2_dp /)
                  k1(1:4) = (/  9.8474045000e-3_dp, -1.4894994500e-2_dp, -2.6158086975e-1_dp, -4.6566358500e+0_dp /)
                  k2(1:4) = (/ -3.2185225000e-4_dp,  4.1585202500e-3_dp,  2.5578953125e-2_dp,  3.0306025000e-1_dp /)
               elseif (tp.le.20._dp) then
                  k0(1:4) = (/  5.7884953000e-2_dp,  3.6113113200e-1_dp,  5.4874661780e+0_dp,  1.2086133940e+2_dp /)
                  k1(1:4) = (/  1.4181754000e-3_dp,  5.2177217500e-2_dp,  8.0110110200e-2_dp,  3.3670318000e-1_dp /)
                  k2(1:4) = (/ -2.4995760000e-5_dp, -1.3176272200e-3_dp, -2.4046612800e-3_dp,  2.4558908000e-2_dp /)
               elseif (tp.le.30._dp) then
                  k0(1:4) = (/  4.1070755000e-2_dp,  3.4189475400e-1_dp,  5.8054307160e+0_dp,  1.0845788300e+2_dp /)
                  k1(1:4) = (/  2.7833601000e-3_dp,  4.7576588000e-2_dp,  3.0414480500e-2_dp,  1.6370367200e+0_dp /)
                  k2(1:4) = (/ -5.1219500000e-5_dp, -1.0395048000e-3_dp, -7.1479114000e-4_dp, -9.4491280000e-3_dp /)
               endif
            elseif (mol_rat.eq.0.5_dp) then
               if (tp.le.2._dp) then
                  k0(1:4) = (/  1.265864600e-2_dp,  4.551181290e-1_dp,  3.564692907e+0_dp,  7.213907520e+1_dp /)
                  k1(1:4) = (/  5.121773100e-2_dp, -2.640971800e-2_dp,  1.000638423e+0_dp,  4.635137980e+1_dp /)
                  k2(1:4) = (/ -1.207794200e-2_dp,  3.576782800e-2_dp, -7.221456200e-2_dp, -1.112042240e+1_dp /)
               elseif (tp.le.4._dp) then
                  k0(1:4) = (/  7.159707300e-2_dp,  5.417509180e-1_dp,  4.200229864e+0_dp,  1.176132419e+2_dp /)
                  k1(1:4) = (/ -2.695067500e-3_dp, -9.961335500e-3_dp,  6.854526175e-1_dp,  2.528834350e+0_dp /)
                  k2(1:4) = (/  1.438505000e-4_dp,  5.885439500e-3_dp, -7.350589850e-2_dp, -5.776913500e-1_dp /)
               elseif (tp.le.6._dp) then
                  k0(1:4) = (/  7.440038900e-2_dp,  4.421759900e-1_dp,  4.912116166e+0_dp,  1.302281541e+2_dp /)
                  k1(1:4) = (/ -4.409216500e-3_dp,  4.593814450e-2_dp,  3.280826280e-1_dp, -4.459285900e+0_dp /)
                  k2(1:4) = (/  3.971805000e-4_dp, -1.865997500e-3_dp, -2.865629500e-2_dp,  3.809067000e-1_dp /)
               elseif (tp.le.10._dp) then
                  k0(1:4) = (/  6.057653300e-2_dp,  4.225489880e-1_dp,  5.815796230e+0_dp,  1.226685405e+2_dp /)
                  k1(1:4) = (/ -4.207225000e-5_dp,  5.258522425e-2_dp,  3.202086450e-2_dp, -1.958664450e+0_dp /)
                  k2(1:4) = (/  5.331912500e-5_dp, -2.428649625e-3_dp, -4.414891750e-3_dp,  1.741257250e-1_dp /)
               elseif (tp.le.20._dp) then
                  k0(1:4) = (/  2.142314300e-2_dp,  4.873508580e-1_dp,  6.243964451e+0_dp,  1.055048727e+2_dp /)
                  k1(1:4) = (/  5.953597600e-3_dp,  3.250121460e-2_dp, -7.930961130e-2_dp,  1.380479860e+0_dp /)
                  k2(1:4) = (/ -1.547139600e-4_dp, -1.068267360e-3_dp,  2.436473620e-3_dp,  1.184797200e-2_dp /)
               elseif (tp.le.30._dp) then
                  k0(1:4) = (/  6.746009500e-2_dp,  8.894570960e-1_dp,  5.597251625e+0_dp,  9.638841450e+1_dp /)
                  k1(1:4) = (/  7.389980000e-4_dp, -1.248944330e-2_dp, -1.869969920e-2_dp,  2.364411290e+0_dp /)
                  k2(1:4) = (/ -9.076360000e-6_dp,  1.759999400e-4_dp,  1.022760080e-3_dp, -1.455745400e-2_dp /)
               endif
            elseif (mol_rat.eq.0.75_dp) then
               if (tp.le.2._dp) then
                  k0(1:4) = (/  1.037903300e-2_dp,  4.891076980e-1_dp,  3.610363064e+0_dp,  5.584487414e+1_dp /)
                  k1(1:4) = (/  5.721903100e-2_dp, -5.716888600e-2_dp,  9.979526690e-1_dp,  5.856955557e+1_dp /)
                  k2(1:4) = (/ -1.467670600e-2_dp,  4.423163600e-2_dp, -6.126358200e-2_dp, -1.627600622e+1_dp /)
               elseif (tp.le.4._dp) then
                  k0(1:4) = (/  7.469435700e-2_dp,  5.071092910e-1_dp,  4.008762698e+0_dp,  1.144740095e+2_dp /)
                  k1(1:4) = (/ -5.547543000e-3_dp,  1.656193850e-2_dp,  8.814215920e-1_dp, -6.076262450e+0_dp /)
                  k2(1:4) = (/  6.277500000e-4_dp,  2.865825500e-3_dp, -1.025979520e-1_dp,  1.389618950e+0_dp /)
               elseif (tp.le.6._dp) then
                  k0(1:4) = (/ -3.796061300e-2_dp, -5.056422810e-1_dp,  9.982265420e-1_dp,  1.114613439e+2_dp /)
                  k1(1:4) = (/  3.890298550e-2_dp,  4.166505195e-1_dp,  1.890247779e+0_dp, -2.904848500e-1_dp /)
                  k2(1:4) = (/ -3.443946500e-3_dp, -3.385934650e-2_dp, -1.666459890e-1_dp,  1.314661500e-1_dp /)
               elseif (tp.le.10._dp) then
                  k0(1:4) = (/  3.179900500e-2_dp, -4.788244200e-2_dp,  5.127685238e+0_dp,  1.100633439e+2_dp /)
                  k1(1:4) = (/  9.844493250e-3_dp,  2.097073562e-1_dp,  3.613406962e-1_dp,  2.846281750e-1_dp /)
                  k2(1:4) = (/ -5.386316250e-4_dp, -1.208437038e-2_dp, -2.653532788e-2_dp,  7.444731250e-2_dp /)
               elseif (tp.le.20._dp) then
                  k0(1:4) = (/  7.102111000e-2_dp,  1.384587938e+0_dp,  8.927256113e+0_dp,  1.038232912e+2_dp /)
                  k1(1:4) = (/  5.642975000e-4_dp, -7.380669090e-2_dp, -4.078991292e-1_dp,  1.472737430e+0_dp /)
                  k2(1:4) = (/ -2.833100000e-6_dp,  1.942330540e-3_dp,  1.239294592e-2_dp,  1.803691400e-2_dp /)
               elseif (tp.le.30._dp) then
                  k0(1:4) = (/  9.562258400e-2_dp,  1.604808182e+0_dp,  7.081416023e+0_dp,  9.700865080e+1_dp /)
                  k1(1:4) = (/ -1.580959000e-3_dp, -8.134269350e-2_dp, -1.612857771e-1_dp,  2.355046970e+0_dp /)
                  k2(1:4) = (/  4.292604000e-5_dp,  1.768580060e-3_dp,  4.676878540e-3_dp, -9.041962000e-3_dp /)
               endif
            elseif (mol_rat.eq.1._dp) then
               if (tp.le.2._dp) then
                  k0(1:4) = (/ -2.443296000e-3_dp,  2.184777200e-1_dp,  2.944584676e+0_dp,  5.957355714e+1_dp /)
                  k1(1:4) = (/  8.144882500e-2_dp,  3.664573060e-1_dp,  2.057737171e+0_dp,  4.234449532e+1_dp /)
                  k2(1:4) = (/ -2.379021400e-2_dp, -9.893208400e-2_dp, -4.100998820e-1_dp, -9.897833920e+0_dp /)
               elseif (tp.le.4._dp) then
                  k0(1:4) = (/  7.480039900e-2_dp,  4.823945990e-1_dp,  3.943178351e+0_dp,  9.414267060e+1_dp /)
                  k1(1:4) = (/ -6.482639500e-3_dp,  3.517206150e-2_dp,  9.738590265e-1_dp,  6.925639650e+0_dp /)
                  k2(1:4) = (/  8.645945000e-4_dp,  7.313185000e-4_dp, -1.178092285e-1_dp, -8.306844500e-1_dp /)
               elseif (tp.le.6._dp) then
                  k0(1:4) = (/  6.427893300e-2_dp,  3.998058210e-1_dp,  5.207776735e+0_dp,  1.074087334e+2_dp /)
                  k1(1:4) = (/ -1.346511000e-3_dp,  7.968439000e-2_dp,  3.244209105e-1_dp, -3.368132500e-1_dp /)
                  k2(1:4) = (/  2.381540000e-4_dp, -5.234965000e-3_dp, -3.448709850e-2_dp,  1.557998500e-1_dp /)
               elseif (tp.le.10._dp) then
                  k0(1:4) = (/  9.382017000e-2_dp,  7.351173870e-1_dp,  7.040362339e+0_dp,  1.049428687e+2_dp /)
                  k1(1:4) = (/ -9.479518000e-3_dp, -1.954280775e-2_dp, -2.675828095e-1_dp,  3.844726500e-1_dp /)
                  k2(1:4) = (/  7.730652500e-4_dp,  1.988691125e-3_dp,  1.327503250e-2_dp,  1.040817750e-1_dp /)
               elseif (tp.le.20._dp) then
                  k0(1:4) = (/  5.879357700e-2_dp,  6.971385550e-1_dp,  6.144833522e+0_dp,  9.309310980e+1_dp /)
                  k1(1:4) = (/  2.310565000e-3_dp,  9.334062900e-3_dp, -7.551552000e-2_dp,  2.746668550e+0_dp /)
                  k2(1:4) = (/ -5.567712000e-5_dp, -5.192076200e-4_dp,  3.023591720e-3_dp, -1.364022600e-2_dp /)
               elseif (tp.le.30._dp) then
                  k0(1:4) = (/  7.627320100e-2_dp,  8.330065730e-1_dp,  4.838028000e+0_dp,  8.876908040e+1_dp /)
                  k1(1:4) = (/  3.199630000e-4_dp, -1.230597200e-2_dp,  4.691371090e-2_dp,  3.200894780e+0_dp /)
                  k2(1:4) = (/  1.539200000e-7_dp,  2.231240800e-4_dp,  1.691439800e-4_dp, -2.554146400e-2_dp /)
               endif
            elseif (mol_rat.eq.2._dp) then
               if (tp.le.2._dp) then
                  k0(1:4) = (/  6.814734700e-2_dp,  6.504046660e-1_dp,  4.391868256e+0_dp,  7.508422843e+1_dp /)
                  k1(1:4) = (/ -1.037330400e-2_dp, -2.411952700e-1_dp, -2.176479900e-2_dp, -1.557011057e+1_dp /)
                  k2(1:4) = (/  4.183992000e-3_dp,  1.015320640e-1_dp,  2.864415900e-1_dp,  8.094834100e+0_dp /)
               elseif (tp.le.4._dp) then
                  k0(1:4) = (/  7.255491900e-2_dp,  4.614907180e-1_dp,  3.943379962e+0_dp,  6.707187336e+1_dp /)
                  k1(1:4) = (/ -6.308614000e-3_dp,  6.265385000e-2_dp,  1.047454190e+0_dp,  5.649844395e+0_dp /)
                  k2(1:4) = (/  1.049754000e-3_dp, -3.164009000e-3_dp, -1.360458310e-1_dp, -5.120546150e-1_dp /)
               elseif (tp.le.6._dp) then
                  k0(1:4) = (/  1.221046190e-1_dp,  7.520998160e-1_dp,  6.372755034e+0_dp,  7.667165036e+1_dp /)
                  k1(1:4) = (/ -2.741564700e-2_dp, -6.121696250e-2_dp, -1.544710020e-1_dp,  5.511971450e-1_dp /)
                  k2(1:4) = (/  3.229656000e-3_dp,  9.640625500e-3_dp,  1.259952500e-2_dp,  1.626211350e-1_dp /)
               elseif (tp.le.10._dp) then
                  k0(1:4) = (/  5.788464500e-2_dp,  6.103656260e-1_dp,  6.918604503e+0_dp,  7.345466762e+1_dp /)
                  k1(1:4) = (/  3.440489500e-3_dp,  3.560516400e-2_dp, -2.343829417e-1_dp,  1.913612077e+0_dp /)
                  k2(1:4) = (/ -1.291452500e-4_dp, -2.559334750e-3_dp,  1.075569638e-2_dp,  2.491261125e-2_dp /)
               elseif (tp.le.20._dp) then
                  k0(1:4) = (/  7.036769900e-2_dp,  8.969261050e-1_dp,  5.790460822e+0_dp, -5.747176148e+1_dp /)
                  k1(1:4) = (/  1.001521600e-3_dp, -2.692895300e-2_dp, -4.938834490e-2_dp,  2.030830140e+1_dp /)
                  k2(1:4) = (/ -1.007900000e-5_dp,  8.284721600e-4_dp,  3.537673500e-3_dp, -5.052920296e-1_dp /)
               elseif (tp.le.30._dp) then
                  k0(1:4) = (/  8.320686500e-2_dp,  2.432539870e-1_dp,  2.620395694e+0_dp,  7.978914580e+1_dp /)
                  k1(1:4) = (/ -3.579539000e-4_dp,  3.022676050e-2_dp,  2.412842339e-1_dp,  4.154631200e+0_dp /)
                  k2(1:4) = (/  2.579686000e-5_dp, -3.951332200e-4_dp, -3.070792620e-3_dp, -4.076078800e-2_dp /)
               endif
            elseif (mol_rat.eq.4._dp) then
               if (tp.eq.4._dp) then
                  k0(1:4) = (/  1.804554300e-2_dp,  3.948378160e-1_dp, -3.838082580e-1_dp, -2.719218350e+1_dp /)
                  k1(1:4) = (/  4.653235600e-2_dp,  6.547056600e-2_dp,  5.792356922e+0_dp,  1.006233201e+2_dp /)
                  k2(1:4) = (/ -1.188070000e-2_dp,  1.612535600e-2_dp, -1.416839776e+0_dp, -2.639996476e+1_dp /)
               elseif (tp.le.4._dp) then
                  k0(1:4) = (/  6.752071000e-2_dp,  4.474855270e-1_dp,  4.003185597e+0_dp,  5.242165954e+1_dp /)
                  k1(1:4) = (/ -3.475324500e-3_dp,  8.633064550e-2_dp,  1.048119719e+0_dp,  9.805418480e+0_dp /)
                  k2(1:4) = (/  7.543485000e-4_dp, -7.466611500e-3_dp, -1.414696385e-1_dp, -8.944747200e-1_dp /)
               elseif (tp.le.6._dp) then
                  k0(1:4) = (/ -3.325066600e-2_dp,  8.661400700e-2_dp,  5.391173361e+0_dp,  1.228572154e+2_dp /)
                  k1(1:4) = (/  3.779735550e-2_dp,  2.284274815e-1_dp,  2.482335145e-1_dp, -2.314126432e+1_dp /)
                  k2(1:4) = (/ -3.265610500e-3_dp, -2.043635050e-2_dp, -2.824732250e-2_dp,  2.939973740e+0_dp /)
               elseif (tp.le.10._dp) then
                  k0(1:4) = (/  5.996457500e-2_dp,  7.412526830e-1_dp,  6.645531657e+0_dp, -3.970579720e+1_dp /)
                  k1(1:4) = (/  3.537294750e-3_dp, -7.659060000e-4_dp, -1.957984047e-1_dp,  3.050830082e+1_dp /)
                  k2(1:4) = (/ -1.449126250e-4_dp, -4.218602500e-4_dp,  1.091471138e-2_dp, -1.485981212e+0_dp /)
               elseif (tp.le.20._dp) then
                  k0(1:4) = (/  7.105380900e-2_dp,  4.742799510e-1_dp,  4.471821161e+0_dp,  7.217389160e+1_dp /)
                  k1(1:4) = (/  8.778513000e-4_dp,  2.761937610e-2_dp,  1.481553722e-1_dp,  5.114723460e+0_dp /)
                  k2(1:4) = (/  1.013938000e-5_dp, -5.906611400e-4_dp, -1.743561360e-3_dp, -6.542036400e-2_dp /)
               elseif (tp.le.30._dp) then
                  k0(1:4) = (/  7.772109100e-2_dp,  7.931764870e-1_dp,  4.615561915e+0_dp,  8.114608580e+1_dp /)
                  k1(1:4) = (/  8.201060000e-4_dp, -1.613541900e-3_dp,  1.354406157e-1_dp,  4.207385630e+0_dp /)
                  k2(1:4) = (/ -3.641560000e-6_dp,  7.374342000e-5_dp, -1.467175420e-3_dp, -4.248395800e-2_dp /)  
               endif
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation shan2018'
               call shutdown(msg)   
            endif
         
            !-----------------------------------------------------------
            !Computing the absorption coefficient
            !-----------------------------------------------------------
            kappa_func_res = pp*(k2(jj)*tp**2._dp + k1(jj)*tp + k0(jj))
         
         case('stepwise')
            if (mol_rat.le.0.2_dp) then
               kappa_func_res = shan_wsgg_kappa(pp,tp,0.125_dp,jj,'none')
            elseif (mol_rat.le.0.4_dp) then
               kappa_func_res = shan_wsgg_kappa(pp,tp,0.25_dp,jj,'none')
            elseif (mol_rat.le.0.6_dp) then
               kappa_func_res = shan_wsgg_kappa(pp,tp,0.5_dp,jj,'none')
            elseif (mol_rat.le.0.9_dp) then
               kappa_func_res = shan_wsgg_kappa(pp,tp,0.75_dp,jj,'none')
            elseif (mol_rat.le.1.1_dp) then
               kappa_func_res = shan_wsgg_kappa(pp,tp,1._dp,jj,'none')
            elseif (mol_rat.le.2.5_dp) then
               kappa_func_res = shan_wsgg_kappa(pp,tp,2._dp,jj,'none')
            else
               kappa_func_res = shan_wsgg_kappa(pp,tp,4._dp,jj,'none')
            endif     
                
         case('linear')
            if (mol_rat.le.0.125_dp) then
               kappa_func_res = shan_wsgg_kappa(pp,tp,0.125_dp,jj,'none')
            elseif (mol_rat.le.0.25_dp) then
               x1 = 0.125_dp; y1 = shan_wsgg_kappa(pp,tp,x1,jj,'none')
               x2 = 0.250_dp; y2 = shan_wsgg_kappa(pp,tp,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.5_dp) then
               x1 = 0.250_dp; y1 = shan_wsgg_kappa(pp,tp,x1,jj,'none')
               x2 = 0.500_dp; y2 = shan_wsgg_kappa(pp,tp,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.75_dp) then
               x1 = 0.500_dp; y1 = shan_wsgg_kappa(pp,tp,x1,jj,'none')
               x2 = 0.750_dp; y2 = shan_wsgg_kappa(pp,tp,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1._dp) then
               x1 = 0.75_dp; y1 = shan_wsgg_kappa(pp,tp,x1,jj,'none')
               x2 = 1.00_dp; y2 = shan_wsgg_kappa(pp,tp,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2._dp) then
               x1 = 1._dp; y1 = shan_wsgg_kappa(pp,tp,x1,jj,'none')
               x2 = 2._dp; y2 = shan_wsgg_kappa(pp,tp,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.4._dp) then
               x1 = 2._dp; y1 = shan_wsgg_kappa(pp,tp,x1,jj,'none')
               x2 = 4._dp; y2 = shan_wsgg_kappa(pp,tp,x2,jj,'none')
               kappa_func_res =  y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               kappa_func_res = shan_wsgg_kappa(pp,tp,4._dp,jj,'none')
            endif   
      endselect
      
   endfunction shan_wsgg_kappa

!   !====================================================================
!   !Function to compute the absorption coefficient for a single
!   !gray gas for the RAD19 paper
!   !====================================================================
!   real(dp) function rad19_wsgg_kappa(ppc,ppw,pps,jj)
   
!      !-------------------------
!      !Declaration of variables
!      !-------------------------
!      integer,intent(in) :: jj
!      integer :: counter,jm,js,ngm,ngs
!      real(dp),intent(in) :: ppc,pps,ppw
!      real(dp) :: mol_rat
      
!      !---------------------------------------------------------
!      !Getting the indexes for the CO2-H2O mixture and for soot
!      !---------------------------------------------------------
!      ngm = wsgg_number_gray_gases('bordbar2014')
!      call get_wsgg_correlations('rad19')
!      ngs = number_wsgg_gases
!      counter = 1
!      main_loop: do jm=1,ngm+1
!         do js=1,ngs+1
!            if (jj.eq.counter) exit main_loop
!            counter = counter + 1
!         enddo
!      enddo main_loop
      
!      !----------------------------------------------------
!      !Computing the absorption coefficient of the mixture
!      !----------------------------------------------------
!      mol_rat = ppw/(ppc+small)                                       !Molar ratio
!      rad19_wsgg_kappa = bordbar_wsgg_kappa(ppc+ppw,mol_rat,jm) + &
!                          fixed_wsgg_kappa(pps,js,'rad19','soot')
   
!   endfunction rad19_wsgg_kappa

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                      !
!           FUNCTIONS FOR THE TEMPERATURE COEFFICIENT                  !
!                                                                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

   !====================================================================
   !Function to compute the temperature coefficient of
   !a single gray gas for fixed-ratio WSGG correlations
   !====================================================================
   real(dp) recursive function fixed_wsgg_a(ttmp,ttmp_ref,jj,&
      wsgg_spec,species_spec) result(a_func_res)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      character(*),intent(in) :: wsgg_spec
      character(*),optional :: species_spec
      character(100) :: species_spec_aux
      integer,intent(in) :: jj
      integer :: ii
      logical :: transparent_window
      real(dp),intent(in) :: ttmp,ttmp_ref
      real(dp) :: sum_a,sum_b
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      !Set default optional parameter
      species_spec_aux = 'mixture'
      if (present(species_spec)) species_spec_aux = species_spec
      
      !Get the appropriate WSGG correlation
      call get_wsgg_correlations(wsgg_spec,species_spec_aux)
      
      !Define if this gray gas is a transparent window
      transparent_window = .false.
      if (jj.lt.1) transparent_window = .true.
      if (jj.gt.number_wsgg_gases) transparent_window = .true.
      
      !-----------------------------------------------------------------
      !Computing the temperature coefficient of gray gas j
      !-----------------------------------------------------------------
      !Transparent windows
      if (transparent_window) then
         sum_a = 0._dp
         do ii=1,number_wsgg_gases
            sum_a = sum_a + &
               fixed_wsgg_a(ttmp,ttmp_ref,ii,wsgg_spec,species_spec_aux)
         enddo
         a_func_res = 1._dp - sum_a
         
      !Gray gases
      else
         sum_b = 0._dp
         do ii=1,degree_wsgg_polynomial+1
            sum_b = sum_b + b(jj,ii)*((ttmp/ttmp_ref)**(real(ii-1,dp)))
         enddo
         a_func_res = sum_b
      endif      
      
   endfunction fixed_wsgg_a

   !======================================================================
   !Function to compute the temperature coefficient for a single gray gas
   !following different versions of the formulation of Smith et al (1982)
   !======================================================================
   real(dp) recursive function smith_wsgg_a(ttmp,ppc,ppw,mol_rat,&
      jj,interpolation_spec) result(a_func_res)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: ppc,ppw,mol_rat,ttmp
      real(dp) :: bb(1:4),b1(1:3,1:4),b2(1:3,1:4),b3(1:3,1:4),&
         b4(1:3,1:4),b5(1:3,1:4),sum_a,sum_b,x1,x2,y1,y2
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + &
               smith_wsgg_a(ttmp,ppc,ppw,mol_rat,ii,interpolation_spec)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif
      
      !-----------------------------------------------------------------
      !Mounting the b's arrays
      !-----------------------------------------------------------------
      !pc -> 0 atm
      b1(1,1:4) = (/  0.4334e-1_dp,  2.620e-4_dp , -1.560e-7_dp ,  2.565e-11_dp  /)
      b1(2,1:4) = (/ -0.4814e-1_dp,  2.822e-4_dp , -1.794e-7_dp ,  3.274e-11_dp  /)
      b1(3,1:4) = (/  0.5492e-1_dp,  0.1087e-4_dp, -0.3500e-7_dp,  0.9123e-11_dp /)
      
      !p_w/p_c = 1
      b2(1:3,1) = (/ 5.150e-1_dp,   0.7749e-1_dp, 1.907e-1_dp    /)
      b2(1:3,2) = (/ -2.303e-4_dp,  3.399e-4_dp,  -1.824e-4_dp   /)
      b2(1:3,3) = (/ 0.9779e-7_dp,  -2.297e-7_dp, 0.5608e-7_dp   /)
      b2(1:3,4) = (/ -1.494e-11_dp, 3.770e-11_dp, -0.5122e-11_dp /)
      
      !p_w/p_c = 2
      b3(1:3,1) = (/  6.508e-1_dp,  -0.2504e-1_dp, 2.718e-1_dp    /)
      b3(1:3,2) = (/ -5.551e-4_dp,  6.112e-4_dp,   -3.118e-4_dp    /)
      b3(1:3,3) = (/ 3.029e-7_dp,   -3.882e-7_dp,  1.221e-7_dp    /)
      b3(1:3,4) = (/ -5.353e-11_dp, 6.528e-11_dp,  -1.612e-11_dp /)
      
      !pw -> 0 atm
      b4(1,1:4) = (/  5.977e-1_dp, -5.119e-4_dp,  3.042e-7_dp, -5.564e-11_dp /)
      b4(2,1:4) = (/ 0.5677e-1_dp,  3.333e-4_dp, -1.967e-7_dp,  2.718e-11_dp /)
      b4(3,1:4) = (/  1.800e-1_dp, -2.334e-4_dp,  1.008e-7_dp, -1.454e-11_dp /)
      
      !pw = 1 atm
      b5(1,1:4) = (/  6.324e-1_dp , -8.358e-4_dp,  6.135e-7_dp, -13.03e-11_dp /)
      b5(2,1:4) = (/ -0.2016e-1_dp,  7.145e-4_dp, -5.212e-7_dp,  9.868e-11_dp /)
      b5(3,1:4) = (/  3.500e-1_dp , -5.040e-4_dp,  2.425e-7_dp, -3.888e-11_dp /)
      
      !----------------------------------------------------------------------------
      !Computing the temperature coefficient according to the interpolation scheme
      !----------------------------------------------------------------------------
      selectcase(trim(interpolation_spec))
         case('none')
            !Defining the appropriate bb array
            if (mol_rat.eq.0._dp) then
               bb = b1(jj,:)
            elseif (mol_rat.eq.1._dp) then
               bb = b2(jj,:)
            elseif (mol_rat.eq.2._dp) then
               bb = b3(jj,:)
            elseif (mol_rat.eq.smith_infty) then
               if (ppw.le.0.5_dp) then
                  bb = b4(jj,:)
               else
                  bb = b5(jj,:)
               endif
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation smith1982'
               call shutdown(msg)               
            endif
               
            !Computing aj
            sum_b = 0._dp
            do ii=1,4
               sum_b = sum_b + bb(ii)*(ttmp**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
               
         case('stepwise')
            !Defining the appropriate bb array
            if (ppw.le.(0.5_dp*ppc)) then
               bb = b1(jj,:)
            elseif (ppw.le.(1.5_dp*ppc)) then
               bb = b2(jj,:)
            elseif (ppw.le.(2.5_dp*ppc)) then
               bb = b3(jj,:)
            elseif (ppw.le.0.5_dp) then
               bb = b4(jj,:)
            else
               bb = b5(jj,:)
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,4
               sum_b = sum_b + bb(ii)*(ttmp**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('linear')
            if (mol_rat.le.1._dp) then
               x1 = 0._dp; y1 = smith_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 1._dp; y2 = smith_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2._dp) then
               x1 = 1._dp; y1 = smith_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 2._dp; y2 = smith_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.smith_infty) then
               x1 = 2._dp        ; y1 = smith_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = smith_infty  ; y2 = smith_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               a_func_res = smith_wsgg_a(ttmp,ppc,ppw,smith_infty,jj,'none')
            endif
      endselect
         
   endfunction smith_wsgg_a   
   
   !====================================================================
   !Function to compute the temperature coefficient for a single 
   !gray gas following the formulation of Yin et al (2010) 
   !====================================================================
   real(dp) recursive function yin2010_wsgg_a(ttmp,ppc,ppw,mol_rat,&
      jj,interpolation_spec) result(a_func_res)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: ppc,ppw,mol_rat,ttmp
      real(dp) :: bb(1:4),b1(1:4,1:4),b2(1:4,1:4),b3(1:4,1:4),&
         b4(1:4,1:4),b5(1:4,1:4),b6(1:4,1:4),b7(1:4,1:4),&
         b8(1:4,1:4),b9(1:4,1:4),b10(1:4,1:4),sum_a,sum_b,tref
      real(dp) :: x1,x2,y1,y2      
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + &
               yin2010_wsgg_a(ttmp,ppc,ppw,mol_rat,ii,interpolation_spec)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif

      !-----------------------------------------------------------------
      !Mounting the b's arrays
      !-----------------------------------------------------------------
      !pw -> 0, pc -> 0
      b1(1:4,1) = (/  0.778969_dp, -0.011449_dp, -0.007627_dp,  0.080082_dp /)
      b1(1:4,2) = (/ -1.342848_dp,  0.343754_dp,  0.242233_dp, -0.049280_dp /)
      b1(1:4,3) = (/  0.964858_dp, -0.234886_dp, -0.173738_dp,  0.001861_dp /)
      b1(1:4,4) = (/ -0.195747_dp,  0.044008_dp,  0.033868_dp,  0.002232_dp /)
      
      !pw = 0.1 atm, pc = 0.1 atm
      b2(1:4,1) = (/  0.492304_dp,  0.082686_dp,  0.144385_dp,  0.079515_dp /)
      b2(1:4,2) = (/ -0.433789_dp,  0.486294_dp, -0.083662_dp, -0.110361_dp /)
      b2(1:4,3) = (/  0.279329_dp, -0.369752_dp,  0.002003_dp,  0.051379_dp /)
      b2(1:4,4) = (/ -0.057770_dp,  0.070509_dp,  0.003902_dp, -0.007983_dp /)
      
      !pw = 0.3 atm, pc = 0.1 atm
      b3(1:4,1) = (/  0.478371_dp,  0.101065_dp,  0.185155_dp,  0.191665_dp /)
      b3(1:4,2) = (/ -0.608643_dp,  0.204118_dp,  0.299794_dp, -0.277448_dp /)
      b3(1:4,3) = (/  0.475098_dp, -0.202202_dp, -0.240346_dp,  0.133514_dp /)
      b3(1:4,4) = (/ -0.109044_dp,  0.042771_dp,  0.046968_dp, -0.021280_dp /)
      
      !pw/pc = 1/8
      b4(1:4,1) = (/  0.515415_dp,  0.199807_dp,  0.138767_dp,  0.087511_dp /)
      b4(1:4,2) = (/ -0.618162_dp,  0.298581_dp, -0.001851_dp, -0.067295_dp /)
      b4(1:4,3) = (/  0.430921_dp, -0.265758_dp, -0.049353_dp,  0.013489_dp /)
      b4(1:4,4) = (/ -0.092082_dp,  0.052910_dp,  0.013012_dp, -5.540e-6_dp /)

      !pw/pc = 1/4
      b5(1:4,1) = (/  0.486247_dp,  0.213959_dp,  0.181991_dp,  0.106180_dp /)
      b5(1:4,2) = (/ -0.644137_dp,  0.306543_dp, -0.020460_dp, -0.096088_dp /)
      b5(1:4,3) = (/  0.485654_dp, -0.264417_dp, -0.053791_dp,  0.028114_dp /)
      b5(1:4,4) = (/ -0.107808_dp,  0.051889_dp,  0.015058_dp, -0.002443_dp /)
      
      !pw/pc = 1/2
      b6(1:4,1) = (/  0.383225_dp,  0.251481_dp,  0.208239_dp,  0.147259_dp /)
      b6(1:4,2) = (/ -0.510937_dp,  0.161562_dp,  0.070697_dp, -0.156339_dp /)
      b6(1:4,3) = (/  0.442201_dp, -0.150405_dp, -0.135668_dp,  0.057698_dp /)
      b6(1:4,4) = (/ -0.106398_dp,  0.028982_dp,  0.032090_dp, -0.007266_dp /)
      
      !pw/pc = 3/4
      b7(1:4,1) = (/  0.255953_dp,  0.340392_dp,  0.160253_dp,  0.201452_dp /)
      b7(1:4,2) = (/ -0.276222_dp, -0.126902_dp,  0.289548_dp, -0.233937_dp /)
      b7(1:4,3) = (/  0.311285_dp,  0.051357_dp, -0.284144_dp,  0.095159_dp /)
      b7(1:4,4) = (/ -0.084903_dp, -0.010259_dp,  0.060344_dp, -0.013302_dp /)
      
      !pw/pc = 1
      b8(1:4,1) = (/ 0.164048_dp,  0.412652_dp,  0.112364_dp,   0.238339_dp /)
      b8(1:4,2) = (/ -0.087793_dp, -0.339810_dp,  0.450929_dp, -0.288619_dp /)
      b8(1:4,3) = (/  0.195253_dp,  0.197886_dp, -0.388486_dp,  0.121962_dp /)
      b8(1:4,4) = (/ -0.063573_dp, -0.038963_dp,  0.079862_dp, -0.017651_dp /)
      
      !pw/pc = 2
      b9(1:4,1) = (/ -0.002188_dp,  0.546857_dp, -0.001911_dp,  0.317219_dp /)
      b9(1:4,2) = (/  0.286129_dp, -0.714799_dp,  0.764177_dp, -0.415470_dp /)
      b9(1:4,3) = (/ -0.048594_dp,  0.452812_dp, -0.581819_dp,  0.186570_dp /)
      b9(1:4,4) = (/ -0.016243_dp, -0.088841_dp,  0.115069_dp, -0.028335_dp /)

      !pw/pc = 4
      b10(1:4,1) = (/ -0.053999_dp, -0.094953_dp,  0.606525_dp,  0.369661_dp /)
      b10(1:4,2) = (/  0.434975_dp,  0.952010_dp, -0.853216_dp, -0.517493_dp /)
      b10(1:4,3) = (/ -0.152413_dp, -0.696161_dp,  0.545562_dp,  0.244011_dp /)
      b10(1:4,4) = (/  0.005094_dp,  0.136316_dp, -0.107328_dp, -0.038451_dp /)
      
      !-----------------------------------------------------------------
      !Computing the temperature coefficient 
      !according to the interpolation scheme
      !-----------------------------------------------------------------
      tref = 1200._dp                                                   !Reference temperature
      selectcase(trim(interpolation_spec))
         case('none')
            !Defining the appropriate bb array
            if (mol_rat.eq.0._dp) then
               bb = b1(jj,:)
            elseif (mol_rat.eq.0.125_dp) then
               bb = b4(jj,:)
            elseif (mol_rat.eq.0.25_dp) then
               bb = b5(jj,:)
            elseif (mol_rat.eq.0.5_dp) then
               bb = b6(jj,:)
            elseif (mol_rat.eq.0.75_dp) then
               bb = b7(jj,:)
            elseif (mol_rat.eq.1._dp) then
               bb = b8(jj,:)
            elseif (mol_rat.eq.2._dp) then
               bb = b9(jj,:)
            elseif (mol_rat.eq.4._dp) then
               bb = b10(jj,:)
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation yin2010'
               call shutdown(msg)               
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,4
               sum_b = sum_b + bb(ii)*((ttmp/tref)**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('stepwise')
            !Defining the appropriate bb array
            if ((ppc+ppw).le.0.1_dp) then
               bb = b1(jj,:)
            elseif ((ppc+ppw).le.0.3_dp) then
               bb = b2(jj,:)
            elseif ((ppc+ppw).le.0.5_dp) then
               bb = b3(jj,:)
            elseif (ppw.le.(0.2_dp*ppc)) then
               bb = b4(jj,:)
            elseif (ppw.le.(0.4_dp*ppc)) then
               bb = b5(jj,:)
            elseif (ppw.le.(0.6_dp*ppc)) then
               bb = b6(jj,:)
            elseif (ppw.le.(0.9_dp*ppc)) then
               bb = b7(jj,:)
            elseif (ppw.le.(1.1_dp*ppc)) then
               bb = b8(jj,:)
            elseif (ppw.le.(2.5_dp*ppc)) then
               bb = b9(jj,:)
            else
               bb = b10(jj,:)
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,4
               sum_b = sum_b + bb(ii)*((ttmp/tref)**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('linear')
            if (mol_rat.le.0.125_dp) then
               x1 = 0._dp     ; y1 = yin2010_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.125_dp  ; y2 = yin2010_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.250_dp) then
               x1 = 0.125_dp  ; y1 = yin2010_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.250_dp  ; y2 = yin2010_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.500_dp) then
               x1 = 0.250_dp  ; y1 = yin2010_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.500_dp  ; y2 = yin2010_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.750_dp) then
               x1 = 0.500_dp  ; y1 = yin2010_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.750_dp  ; y2 = yin2010_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1.000_dp) then
               x1 = 0.750_dp  ; y1 = yin2010_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 1.000_dp  ; y2 = yin2010_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2.000_dp) then
               x1 = 1.000_dp  ; y1 = yin2010_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 2.000_dp  ; y2 = yin2010_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.4.000_dp) then
               x1 = 2.000_dp  ; y1 = yin2010_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 4.000_dp  ; y2 = yin2010_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               a_func_res = yin2010_wsgg_a(ttmp,ppc,ppw,4._dp,jj,'none')
            endif
      endselect
      
   endfunction yin2010_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for a single gray
   !gas following the formulation of Kangwanpongpan et al (2012) 
   !====================================================================
   real(dp) recursive function kangwanpongpan_variable_wsgg_a(ttmp,mol_rat,jj) &
      result(a_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: mol_rat,ttmp
      real(dp) :: c(1:6),c1(1:4,1:6),c2(1:4,1:6),c3(1:4,1:6)
      real(dp) :: sum_a,sum_b,tref
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + &
               kangwanpongpan_variable_wsgg_a(ttmp,mol_rat,ii)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif
      
      !-----------------------------------------------------------------
      !Mounting the c1, c2 and c3 arrays
      !-----------------------------------------------------------------
      c1(1,1:6) = (/  0.3947_dp, -0.4512_dp,  0.1492_dp,  1.8824_dp, -2.3284_dp,  0.7698_dp /)
      c1(2,1:6) = (/ -0.4974_dp,  6.8986_dp, -19.988_dp,  26.208_dp, -16.440_dp,  3.9847_dp /)
      c1(3,1:6) = (/  0.3189_dp, -0.7222_dp,  1.5053_dp, -1.8378_dp,  1.0337_dp, -0.2107_dp /)
      c1(4,1:6) = (/  0.1648_dp, -0.6012_dp,  2.0308_dp, -3.4361_dp,  2.5803_dp, -0.7069_dp /)

      c2(1,1:6) = (/ -0.1214_dp,  1.1420_dp, -5.2222_dp,  9.1820_dp, -6.9298_dp,  1.9063_dp /)
      c2(2,1:6) = (/  0.1092_dp, -2.3198_dp,  8.0021_dp, -11.0070_dp, 7.1199_dp, -1.7876_dp /)
      c2(3,1:6) = (/ -0.0720_dp,  1.0304_dp, -1.9350_dp,  1.6332_dp, -0.7798_dp,  0.1782_dp /)
      c2(4,1:6) = (/  0.0329_dp,  0.6942_dp, -3.0960_dp,  4.7494_dp, -3.1714_dp,  0.7869_dp /)

      c3(1,1:6) = (/  0.0243_dp, -0.2296_dp,  1.0115_dp, -1.7493_dp,  1.3038_dp, -0.3549_dp /)
      c3(2,1:6) = (/ -0.0179_dp,  0.4077_dp, -1.4482_dp,  2.0311_dp, -1.3278_dp,  0.3349_dp /)
      c3(3,1:6) = (/  0.0158_dp, -0.2478_dp,  0.5931_dp, -0.6619_dp,  0.3857_dp, -0.0933_dp /)
      c3(4,1:6) = (/ -0.0095_dp, -0.0687_dp,  0.3691_dp, -0.5919_dp,  0.4017_dp, -0.1003_dp /)
      
      !-----------------------------------------------------------------
      !Mounting the c array
      !-----------------------------------------------------------------
      do ii=1,6
         c(ii) = c1(jj,ii) + c2(jj,ii)*mol_rat + c3(jj,ii)*(mol_rat**2._dp)
      enddo

      !-----------------------------------------------------------------
      !Computing the temperature coefficient of gray gas j
      !-----------------------------------------------------------------
      tref = 2000._dp; sum_b = 0._dp
      do ii=1,6
         sum_b = sum_b + c(ii)*((ttmp/tref)**real(ii-1,dp))
      enddo
      a_func_res = sum_b
  
   endfunction kangwanpongpan_variable_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for
   !a single gray gas following the formulation of
   !Kangwanpongpan et al (2012) for a fixed molar ratio 
   !====================================================================
   real(dp) recursive function kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,&
      mol_rat,jj,interpolation_spec) result(a_func_res)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: ppc,ppw,mol_rat,ttmp
      real(dp) :: bb(1:6),b1(1:4,1:6),b2(1:4,1:6),b3(1:4,1:6),&
         b4(1:4,1:6),b5(1:4,1:6),b6(1:4,1:6),b7(1:4,1:6),sum_a,sum_b,tref
      real(dp) :: x1,x2,y1,y2
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,&
               mol_rat,ii,interpolation_spec)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif

      !-----------------------------------------------------------------
      !Mounting the b's arrays
      !-----------------------------------------------------------------
      !MR = 0.125
      b1(1:4,1) = (/  0.4329_dp, -0.5151_dp,  0.3461_dp,  0.1443_dp /)
      b1(1:4,2) = (/ -0.7255_dp,  7.2071_dp, -1.1838_dp, -0.4216_dp /)
      b1(1:4,3) = (/  1.1384_dp, -21.319_dp,  3.0944_dp,  1.5071_dp /)
      b1(1:4,4) = (/  0.2661_dp,  28.355_dp, -4.1157_dp, -2.7378_dp /)
      b1(1:4,5) = (/ -1.1786_dp, -17.950_dp,  2.5339_dp,  2.1386_dp /)
      b1(1:4,6) = (/  0.4713_dp,  4.3778_dp, -0.5857_dp, -0.5994_dp /)
      
      !MR = 0.25
      b2(1:4,1) = (/  0.3607_dp, -0.4713_dp,  0.2999_dp,  0.1755_dp /)
      b2(1:4,2) = (/ -0.1849_dp,  6.3750_dp, -0.4354_dp, -0.4710_dp /)
      b2(1:4,3) = (/ -1.0644_dp, -18.102_dp,  0.8603_dp,  1.4103_dp /)
      b2(1:4,4) = (/  4.0610_dp,  23.534_dp, -1.1531_dp, -2.4705_dp /)
      b2(1:4,5) = (/ -3.9878_dp, -14.680_dp,  0.6442_dp,  1.9323_dp /)
      b2(1:4,6) = (/  1.2272_dp,  3.5389_dp, -0.1167_dp, -0.5459_dp /)
      
      !MR = 0.5
      b3(1:4,1) = (/  0.3026_dp, -0.4235_dp,  0.2614_dp,  0.1954_dp /)
      b3(1:4,2) = (/  0.3705_dp,  5.4049_dp,  0.1391_dp, -0.3374_dp /)
      b3(1:4,3) = (/ -3.3957_dp, -14.682_dp, -0.5843_dp,  0.6829_dp /)
      b3(1:4,4) = (/  7.9979_dp,  18.707_dp,  0.5412_dp, -1.3044_dp /)
      b3(1:4,5) = (/ -6.8847_dp, -11.509_dp, -0.3740_dp,  1.1421_dp /)
      b3(1:4,6) = (/  2.0098_dp,  2.7377_dp,  0.1325_dp, -0.3490_dp /)
      
      !MR = 0.75
      b4(1:4,1) = (/  0.2887_dp, -0.4066_dp,  0.2529_dp,  0.1974_dp /)
      b4(1:4,2) = (/  0.5335_dp,  4.9951_dp,  0.2466_dp, -0.1456_dp /)
      b4(1:4,3) = (/ -4.2296_dp, -13.339_dp, -0.6019_dp, -0.1062_dp /)
      b4(1:4,4) = (/  9.4993_dp,  16.952_dp,  0.3143_dp, -0.1181_dp /)
      b4(1:4,5) = (/ -8.0354_dp, -10.410_dp, -0.1583_dp,  0.3554_dp /)
      b4(1:4,6) = (/  2.3307_dp,  2.4659_dp,  0.0746_dp, -0.1542_dp /)
      
      !MR = 1.0
      b5(1:4,1) = (/  0.2853_dp, -0.3976_dp,  0.2532_dp,  0.1941_dp /)
      b5(1:4,2) = (/  0.5853_dp,  4.7799_dp,  0.2151_dp,  0.0337_dp /)
      b5(1:4,3) = (/ -4.5894_dp, -12.689_dp, -0.2515_dp, -0.7945_dp /)
      b5(1:4,4) = (/  10.205_dp,  16.172_dp, -0.3589_dp,  0.8954_dp /)
      b5(1:4,5) = (/ -8.6057_dp, -9.9506_dp,  0.3293_dp, -0.3114_dp /)
      b5(1:4,6) = (/  2.4960_dp,  2.3557_dp, -0.0498_dp,  0.0104_dp /)
      
      !MR = 2.0
      b6(1:4,1) = (/  0.2872_dp, -0.3748_dp,  0.2655_dp,  0.1750_dp /)
      b6(1:4,2) = (/  0.5807_dp,  4.3960_dp, -0.0908_dp,  0.5481_dp /)
      b6(1:4,3) = (/ -4.9092_dp, -11.671_dp,  1.3035_dp, -2.6576_dp /)
      b6(1:4,4) = (/  11.012_dp,  15.104_dp, -2.9227_dp,  3.5816_dp /)
      b6(1:4,5) = (/ -9.3465_dp, -9.3820_dp,  2.0978_dp, -2.0627_dp /)
      b6(1:4,6) = (/  2.7294_dp,  2.2263_dp, -0.4950_dp,  0.4412_dp /)
      
      !MR = 4.0
      b7(1:4,1) = (/  0.2911_dp, -0.3416_dp,  0.2790_dp,  0.1471_dp /)
      b7(1:4,2) = (/  0.5050_dp,  4.0487_dp, -0.4828_dp,  1.0690_dp /)
      b7(1:4,3) = (/ -4.8039_dp, -10.799_dp,  3.0103_dp, -4.4490_dp /)
      b7(1:4,4) = (/  11.037_dp,  14.158_dp, -5.5730_dp,  6.1082_dp /)
      b7(1:4,5) = (/ -9.4886_dp, -8.8559_dp,  3.8801_dp, -3.6923_dp /)
      b7(1:4,6) = (/  2.7969_dp,  2.1031_dp, -0.9391_dp,  0.8400_dp /)   
      
      !-----------------------------------------------------------------
      !Computing the temperature coefficient 
      !according to the interpolation scheme
      !-----------------------------------------------------------------
      tref = 2000._dp                                                   !Reference temperature
      selectcase(trim(interpolation_spec))
         case('none')
            !Defining the appropriate bb array
            if (mol_rat.eq.0.125_dp) then
               bb = b1(jj,:)
            elseif (mol_rat.eq.0.25_dp) then
               bb = b2(jj,:)
            elseif (mol_rat.eq.0.5_dp) then
               bb = b3(jj,:)
            elseif (mol_rat.eq.0.75_dp) then
               bb = b4(jj,:)
            elseif (mol_rat.eq.1._dp) then
               bb = b5(jj,:)
            elseif (mol_rat.eq.2._dp) then
               bb = b6(jj,:)
            elseif (mol_rat.eq.4._dp) then
               bb = b7(jj,:)
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation kangwanpongpan2012_fixed'
               call shutdown(msg)               
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,6
               sum_b = sum_b + bb(ii)*((ttmp/tref)**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('stepwise')
            !Defining the appropriate bb array
            if (ppw.le.(0.1875_dp*ppc)) then            
               bb = b1(jj,:)
            elseif (ppw.le.(0.375_dp*ppc)) then
               bb = b2(jj,:)
            elseif (ppw.le.(0.625_dp*ppc)) then
               bb = b3(jj,:)
            elseif (ppw.le.(0.875_dp*ppc)) then
               bb = b4(jj,:)
            elseif (ppw.le.(1.5_dp*ppc)) then
               bb = b5(jj,:)
            elseif (ppw.le.(3._dp*ppc)) then
               bb = b6(jj,:)
            else
               bb = b7(jj,:)
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,6
               sum_b = sum_b + bb(ii)*((ttmp/tref)**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('linear')
            if (mol_rat.le.0.125_dp) then
               x1 = 0._dp     ; y1 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.125_dp  ; y2 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.250_dp) then
               x1 = 0.125_dp  ; y1 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.250_dp  ; y2 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.500_dp) then
               x1 = 0.250_dp  ; y1 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.500_dp  ; y2 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.750_dp) then
               x1 = 0.500_dp  ; y1 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.750_dp  ; y2 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1.000_dp) then
               x1 = 0.750_dp  ; y1 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 1.000_dp  ; y2 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2.000_dp) then
               x1 = 1.000_dp  ; y1 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 2.000_dp  ; y2 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.4.000_dp) then
               x1 = 2.000_dp  ; y1 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 4.000_dp  ; y2 = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               a_func_res = kangwanpongpan_fixed_wsgg_a(ttmp,ppc,ppw,4._dp,jj,'none')
            endif            
      endselect
      
   endfunction kangwanpongpan_fixed_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for a single 
   !gray gas following the formulation of Johansson et al (2014) 
   !====================================================================
   real(dp) recursive function johansson_wsgg_a(ttmp,mol_rat,jj) &
      result(a_func_res)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: mol_rat,ttmp
      real(dp) :: c1(1:4,1:3),c2(1:4,1:3),c3(1:4,1:3),tref
      real(dp) :: sum_a,sum_b,sum_c
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + johansson_wsgg_a(ttmp,mol_rat,ii)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif
      
      !-----------------------------------------------------------------
      !Mounting the c's array
      !-----------------------------------------------------------------
      c1(1:4,1) = (/  0.3580_dp,  0.3920_dp,  0.1420_dp,  0.0798_dp /)
      c1(1:4,2) = (/  0.0731_dp, -0.2120_dp, -0.0831_dp, -0.0370_dp /)
      c1(1:4,3) = (/ -0.0466_dp,  0.0191_dp,  0.0148_dp,  0.0023_dp /)
      c2(1:4,1) = (/ -0.1650_dp, -0.2910_dp,  0.3480_dp,  0.0866_dp /)
      c2(1:4,2) = (/ -0.0554_dp,  0.6440_dp, -0.2940_dp, -0.1060_dp /)
      c2(1:4,3) = (/  0.0930_dp, -0.2090_dp,  0.0662_dp,  0.0305_dp /)
      c3(1:4,1) = (/  0.0598_dp,  0.0784_dp, -0.1220_dp, -0.0127_dp /)
      c3(1:4,2) = (/  0.0028_dp, -0.1970_dp,  0.1180_dp,  0.0169_dp /)
      c3(1:4,3) = (/ -0.0256_dp,  0.0662_dp, -0.0295_dp, -0.0051_dp /)

      !-----------------------------------------------------------------
      !Computing the temperature coefficient of gray gas j
      !-----------------------------------------------------------------
      sum_b = 0._dp; tref = 1200._dp
      do ii=1,3
         sum_c = c1(jj,ii) + c2(jj,ii)*mol_rat + c3(jj,ii)*mol_rat**2._dp
         sum_b = sum_b + sum_c*(ttmp/tref)**(real(ii-1,dp))
      enddo
      a_func_res = sum_b
      
   endfunction johansson_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for a single
   !gray gas following the formulation of Krishnamoorthy (2013) 
   !====================================================================
   real(dp) recursive function krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,&
      mol_rat,jj,interpolation_spec) result(a_func_res)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: ppc,ppw,mol_rat,ttmp
      real(dp) :: bb(1:2),b1(1:4,1:2),b2(1:4,1:2),b3(1:4,1:2),&
         b4(1:4,1:2),sum_a,sum_b,x1,x2,y1,y2
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,&
                                          mol_rat,ii,interpolation_spec)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif

      !-----------------------------------------------------------------
      !Mounting the b's arrays
      !-----------------------------------------------------------------
      !MR = 0.11
      b1(1:4,1) = (/ 2.39641e-01_dp,  3.42342e-01_dp,&
                     1.37773e-01_dp,  4.40724e-02_dp /)
      b1(1:4,2) = (/ 7.85445e-05_dp, -9.47416e-05_dp,& 
                    -5.51091e-05_dp,  7.26634e-06_dp /)
      
      !MR = 0.5
      b2(1:4,1) = (/ 1.89029e-01_dp,  2.87021e-01_dp,&
                     2.54516e-01_dp,  6.54289e-02_dp /)
      b2(1:4,2) = (/ 9.33340e-05_dp, -5.32833e-05_dp,&
                    -1.01806e-04_dp, -2.25973e-05_dp /)
      
      !MR = 1.0
      b3(1:4,1) = (/ 1.91464e-01_dp,  2.34876e-01_dp,&
                     2.47320e-01_dp,  9.59426e-02_dp /)
      b3(1:4,2) = (/ 9.22363e-05_dp, -4.25444e-05_dp,&
                    -9.89282e-05_dp, -3.83770e-05_dp /)
      
      !MR = 2.0
      b4(1:4,1) = (/ 1.54129e-01_dp,  2.43637e-01_dp,&
                     2.84084e-01_dp,  8.57853e-02_dp /)
      b4(1:4,2) = (/ 1.07579e-04_dp, -3.09769e-05_dp,&
                    -1.13634e-04_dp, -3.43141e-05_dp /)
      
      !-----------------------------------------------------------------
      !Computing the temperature coefficient 
      !according to the interpolation scheme
      !-----------------------------------------------------------------
      selectcase(trim(interpolation_spec))
         case('none')
            !Defining the appropriate bb array
            if (mol_rat.eq.0.11_dp) then
               bb = b1(jj,:)
            elseif (mol_rat.eq.0.5_dp) then
               bb = b2(jj,:)
            elseif (mol_rat.eq.1._dp) then
               bb = b3(jj,:)
            elseif (mol_rat.eq.2._dp) then
               bb = b4(jj,:)
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation krishnamoorthy2013'
               call shutdown(msg)               
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,2
               sum_b = sum_b + bb(ii)*(ttmp**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('stepwise')
            !Defining the appropriate bb array
            if (ppw.le.(0.2_dp*ppc)) then
               bb = b1(jj,:)
            elseif (ppw.le.(0.67_dp*ppc)) then
               bb = b2(jj,:)
            elseif (ppw.le.(1.5_dp*ppc)) then
               bb = b3(jj,:)
            else
               bb = b4(jj,:)
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,2
               sum_b = sum_b + bb(ii)*(ttmp**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('linear')
            if (mol_rat.le.0.11_dp) then
               a_func_res = krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,0.11_dp,jj,'none')
            elseif (mol_rat.le.0.5_dp) then
               x1 = 0.11_dp  ; y1 = krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 0.50_dp  ; y2 = krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1._dp) then
               x1 = 0.5_dp  ; y1 = krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 1.0_dp  ; y2 = krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2._dp) then
               x1 = 1._dp  ; y1 = krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
               x2 = 2._dp  ; y2 = krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               a_func_res = krishnamoorthy2013_wsgg_a(ttmp,ppc,ppw,2._dp,jj,'none')
            endif
      endselect
      
   endfunction krishnamoorthy2013_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for a 
   !single gray gas following the formulation of Yin (2013) 
   !====================================================================
   real(dp) recursive function yin2013_wsgg_a(ttmp,ppc,ppw,mol_rat,&
      jj,interpolation_spec) result(a_func_res)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*),intent(in) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: ppc,ppw,mol_rat,ttmp
      real(dp) :: bb(1:4),b1(1:4,1:4),b2(1:4,1:4),b3(1:4,1:4),&
         b4(1:4,1:4),b5(1:4,1:4),b6(1:4,1:4),b7(1:4,1:4),&
         sum_a,sum_b,tref,x1,x2,y1,y2
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + &
               yin2013_wsgg_a(ttmp,ppc,ppw,mol_rat,ii,interpolation_spec)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif

      !-----------------------------------------------------------------
      !Mounting the b's arrays
      !-----------------------------------------------------------------
      !pc -> 0 atm
      b1(1:4,1) = (/  0.204623_dp, -0.020227_dp,&
                      0.044221_dp,  0.039311_dp /)
      b1(1:4,2) = (/ -0.378060_dp,  0.256006_dp,&
                      0.003850_dp, -0.054832_dp /)
      b1(1:4,3) = (/  0.666639_dp, -0.195201_dp,&
                     -0.020175_dp,  0.025370_dp /)
      b1(1:4,4) = (/ -0.203453_dp,  0.040493_dp,&
                      0.004919_dp, -0.003891_dp /)
      
      !p_w/p_c = 0.05
      b2(1:4,1) = (/  0.315106_dp,  0.092474_dp,&
                      0.031702_dp,  0.046138_dp /)
      b2(1:4,2) = (/  0.023475_dp,  0.109146_dp,&
                      0.037396_dp, -0.061392_dp /)
      b2(1:4,3) = (/ -0.057930_dp, -0.121000_dp,&
                     -0.040731_dp,  0.027164_dp /)
      b2(1:4,4) = (/  0.008408_dp,  0.027145_dp,&
                      0.008742_dp, -0.003996_dp /)
      
      !p_w/p_c = 1
      b3(1:4,1) = (/  0.500119_dp,  0.071592_dp,&
        0.155320_dp,  0.072615_dp /)
      b3(1:4,2) = (/ -0.447068_dp,  0.508252_dp,&
       -0.104294_dp, -0.100601_dp /)
      b3(1:4,3) = (/  0.286878_dp, -0.384253_dp,&
        0.014096_dp,  0.046681_dp /)
      b3(1:4,4) = (/ -0.059165_dp,  0.073477_dp,&
        0.001643_dp, -0.007224_dp /)
      
      !p_w/p_c = 2
      b4(1:4,1) = (/  0.542458_dp,  0.101734_dp,&
                      0.146066_dp,  0.129511_dp /)
      b4(1:4,2) = (/ -0.658411_dp,  0.518429_dp,&
                     -0.008745_dp, -0.187993_dp /)
      b4(1:4,3) = (/  0.466444_dp, -0.386151_dp,&
                     -0.058325_dp,  0.090709_dp /)
      b4(1:4,4) = (/ -0.100186_dp,  0.073453_dp,&
                      0.015984_dp, -0.014493_dp /)
      
      !pw -> 0 atm
      b5(1:4,1) = (/  0.966357_dp,  0.662059_dp,&
                      0.060870_dp,  0.103568_dp /)
      b5(1:4,2) = (/ -0.790165_dp, -2.262877_dp,&
                      0.436788_dp, -0.153135_dp /)
      b5(1:4,3) = (/ -0.050144_dp,  2.309473_dp,&
                     -0.395493_dp,  0.074910_dp /)
      b5(1:4,4) = (/  0.115202_dp, -0.572895_dp,&
                      0.085146_dp, -0.012091_dp /)
      
      !pw = 0.05 atm
      b6(1:4,1) = (/  0.340618_dp,  0.175818_dp,&
                      0.044325_dp,  0.126628_dp /)
      b6(1:4,2) = (/ -0.105469_dp, -0.063466_dp,&
                      0.288376_dp, -0.186480_dp /)
      b6(1:4,3) = (/  0.068051_dp,  0.086631_dp,&
                     -0.258205_dp,  0.090755_dp /)
      b6(1:4,4) = (/ -0.017828_dp, -0.026581_dp,&
                      0.054333_dp, -0.014569_dp /)
      
      !pw = 1 atm
      b7(1:4,1) = (/ -0.077336_dp,  0.506777_dp,&
                     -0.079989_dp,  0.373898_dp /)
      b7(1:4,2) = (/  0.661776_dp, -0.758948_dp,&
                      0.851078_dp, -0.540887_dp /)
      b7(1:4,3) = (/ -0.362515_dp,  0.516146_dp,&
                     -0.604264_dp,  0.258923_dp /)
      b7(1:4,4) = (/  0.053534_dp, -0.102909_dp,&
                      0.113500_dp, -0.040957_dp /)
      
      !-----------------------------------------------------------------
      !Computing the temperature coefficient
      !according to the interpolation scheme
      !-----------------------------------------------------------------
      tref = 1200._dp                                                   !Reference temperature
      selectcase(trim(interpolation_spec))
         case('none')
            !Defining the appropriate bb array
            if (mol_rat.eq.0._dp) then
               bb = b1(jj,:)
            elseif (mol_rat.eq.0.05_dp) then
               bb = b2(jj,:)
            elseif (mol_rat.eq.1._dp) then
               bb = b3(jj,:)
            elseif (mol_rat.eq.2._dp) then
               bb = b4(jj,:)
            elseif (mol_rat.eq.yin_infty) then
               if (ppw.le.0.01) then
                  bb = b5(jj,:)
               elseif (ppw.le.0.2) then
                  bb = b6(jj,:)
               else
                  bb = b7(jj,:)
               endif
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation yin2013'
               call shutdown(msg)               
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,4
               sum_b = sum_b + bb(ii)*((ttmp/tref)**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('stepwise')
            !Defining the appropriate bb array
            if (ppw.le.(0.01_dp*ppc)) then
               bb = b1(jj,:)
            elseif (ppw.le.(0.5_dp*ppc)) then
               bb = b2(jj,:)
            elseif (ppw.le.(1.5_dp*ppc)) then
               bb = b3(jj,:)
            elseif (ppw.le.(2.5_dp*ppc)) then
               bb = b4(jj,:)
            elseif (ppw.le.0.01_dp) then
               bb = b5(jj,:)
            elseif (ppw.le.0.2_dp) then
               bb = b6(jj,:)
            else
               bb = b7(jj,:)
            endif
            
            !Computing aj
            sum_b = 0._dp
            do ii=1,4
               sum_b = sum_b + bb(ii)*((ttmp/tref)**(real(ii-1,dp)))
            enddo
            a_func_res = sum_b
            
         case('linear')
            if (mol_rat.le.0.05_dp) then
                  x1 = 0.00_dp  ; y1 = yin2013_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
                  x2 = 0.05_dp  ; y2 = yin2013_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
                  a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
               elseif (mol_rat.le.1._dp) then
                  x1 = 0.05_dp  ; y1 = yin2013_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
                  x2 = 1.00_dp  ; y2 = yin2013_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
                  a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
               elseif (mol_rat.le.2._dp) then
                  x1 = 1._dp  ; y1 = yin2013_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
                  x2 = 2._dp  ; y2 = yin2013_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
                  a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
               elseif (mol_rat.le.yin_infty) then
                  x1 = 2._dp     ; y1 = yin2013_wsgg_a(ttmp,ppc,ppw,x1,jj,'none')
                  x2 = yin_infty ; y2 = yin2013_wsgg_a(ttmp,ppc,ppw,x2,jj,'none')
                  a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
               else
                  a_func_res = yin2013_wsgg_a(ttmp,ppc,ppw,yin_infty,jj,'none')
               endif
      endselect
      
   endfunction yin2013_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for a single
   !gray gas following the formulation of Bordbar et al (2014) 
   !====================================================================
   real(dp) recursive function bordbar_wsgg_a(ttmp,mol_rat,jj) &
      result(a_func_res)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      integer :: ii,kk
      real(dp),intent(in) :: mol_rat,ttmp
      real(dp) :: bb(1:4,0:4),c(1:4,0:4,0:4),tref,sum_a,sum_b,sum_c,&
         x1,x2,y1,y2,mr
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + bordbar_wsgg_a(ttmp,mol_rat,ii)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif
      
      !-----------------------------------------------------------------
      !Computing the temperature coefficient for the gray gases 
      !(within five intervals of molar ratio)
      !-----------------------------------------------------------------
      tref = 1200._dp                                                   !Reference temperature
      if ((mol_rat.eq.0._dp).and.bordbar_co2_interpolation) then        !Only CO2
         !Mounting the b's array
         bb(1,0:4) = (/  8.425766e-1_dp, -1.442229e+0_dp,  1.286974e+0_dp, -5.202712e-1_dp,  7.581559e-2_dp /)
         bb(2,0:4) = (/ -3.023864e-2_dp,  5.264245e-1_dp, -6.209696e-1_dp,  2.704755e-1_dp, -4.090690e-2_dp /)
         bb(3,0:4) = (/  1.070243e-1_dp, -1.989596e-1_dp,  3.101602e-1_dp, -1.737230e-1_dp,  3.081180e-2_dp /)
         bb(4,0:4) = (/  3.108972e-2_dp,  1.981489e-1_dp, -2.543676e-1_dp,  1.061331e-1_dp, -1.498231e-2_dp /)
            
         !Computing the polynomial
         sum_b = 0._dp
         do ii=0,4
            sum_b = sum_b + bb(jj,ii)*(ttmp/tref)**(real(ii,dp))
         enddo
         a_func_res = sum_b
            
      elseif ((mol_rat.ge.bordbar_infty).and.&
              bordbar_h2o_interpolation) then                           !Only H2O (virtually)
         !Mounting the b's array
         bb(1,0:4) = (/  7.129509e-1_dp, -1.378353e+0_dp,  1.555028e+0_dp, -6.636291e-1_dp,  9.773674e-2_dp /)
         bb(2,0:4) = (/  1.589917e-1_dp,  5.635578e-2_dp,  2.666874e-1_dp, -2.040335e-1_dp,  3.742408e-2_dp /)
         bb(3,0:4) = (/ -1.196373e-1_dp,  1.349665e+0_dp, -1.544797e+0_dp,  6.397595e-1_dp, -9.153650e-2_dp /)
         bb(4,0:4) = (/  3.078250e-1_dp, -6.003555e-1_dp,  4.441261e-1_dp, -1.468813e-1_dp,  1.824702e-2_dp /)
         
         !Computing the polynomial
         sum_b = 0._dp
         do ii=0,4
            sum_b = sum_b + bb(jj,ii)*(ttmp/tref)**(real(ii,dp))
         enddo
         a_func_res = sum_b
            
      elseif ((mol_rat.lt.0.01_dp).and.bordbar_co2_interpolation) then  !0<Molar ratio<0.01
         !Bounds for the interpolation
         x1 = 0._dp      ; y1 = bordbar_wsgg_a(ttmp,x1,jj)
         x2 = 0.01_dp   ; y2 = bordbar_wsgg_a(ttmp,x2,jj)
      
         !Interpolating
         a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
      
      elseif ((mol_rat.gt.4._dp).and.bordbar_h2o_interpolation) then    !4<Molar ratio<bordbar_infty
         !Values for the interpolation
         x1 = 4._dp            ; y1 = bordbar_wsgg_a(ttmp,x1,jj)
         x2 = bordbar_infty   ; y2 = bordbar_wsgg_a(ttmp,x2,jj)
         
         !Interpolating
         a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
         
      else
         mr = min(max(mol_rat,0.01_dp),4._dp)
                                                                        !CO2-H2O mixture
         !Mounting the c's array
         c(1,0,0:4) = (/  0.7412956_dp, -0.5244441_dp,  0.5822860_dp, -0.2096994_dp,  0.0242031_dp /)
         c(1,1,0:4) = (/ -0.9412652_dp,  0.2799577_dp, -0.7672319_dp,  0.3204027_dp, -0.0391017_dp /)
         c(1,2,0:4) = (/  0.8531866_dp,  0.0823075_dp,  0.5289430_dp, -0.2468463_dp,  0.0310940_dp /)
         c(1,3,0:4) = (/ -0.3342806_dp,  0.1474987_dp, -0.4160689_dp,  0.1697627_dp, -0.0204066_dp /)
         c(1,4,0:4) = (/  0.0431436_dp, -0.0688622_dp,  0.1109773_dp, -0.0420861_dp,  0.0049188_dp /)
         c(2,0,0:4) = (/  0.1552073_dp, -0.4862117_dp,  0.3668088_dp, -0.1055508_dp,  0.0105857_dp /)
         c(2,1,0:4) = (/  0.6755648_dp,  1.4092710_dp, -1.3834490_dp,  0.4575210_dp, -0.0501976_dp /)
         c(2,2,0:4) = (/ -1.1253940_dp, -0.5913199_dp,  0.9085441_dp, -0.3334201_dp,  0.0384236_dp /)
         c(2,3,0:4) = (/  0.6040543_dp, -0.0553385_dp, -0.1733014_dp,  0.0791608_dp, -0.0098934_dp /)
         c(2,4,0:4) = (/ -0.1105453_dp,  0.0464663_dp, -0.0016129_dp, -0.0035398_dp,  0.0006121_dp /)
         c(3,0,0:4) = (/  0.2550242_dp,  0.3805403_dp, -0.4249709_dp,  0.1429446_dp, -0.0157408_dp /)
         c(3,1,0:4) = (/ -0.6065428_dp,  0.3494024_dp,  0.1853509_dp, -0.1013694_dp,  0.0130244_dp /)
         c(3,2,0:4) = (/  0.8123855_dp, -1.1020090_dp,  0.4046178_dp, -0.0811822_dp,  0.0062981_dp /)
         c(3,3,0:4) = (/ -0.4532290_dp,  0.6784475_dp, -0.3432603_dp,  0.0883088_dp, -0.0084152_dp /)
         c(3,4,0:4) = (/  0.0869309_dp, -0.1306996_dp,  0.0741446_dp, -0.0202929_dp,  0.0020110_dp /)
         c(4,0,0:4) = (/ -0.0345199_dp,  0.2656726_dp, -0.1225365_dp,  0.0300151_dp, -0.0028205_dp /)
         c(4,1,0:4) = (/  0.4112046_dp, -0.5728350_dp,  0.2924490_dp, -0.0798076_dp,  0.0079966_dp /)
         c(4,2,0:4) = (/ -0.5055995_dp,  0.4579559_dp, -0.2616436_dp,  0.0764841_dp, -0.0079084_dp /)
         c(4,3,0:4) = (/  0.2317509_dp, -0.1656759_dp,  0.1052608_dp, -0.0321935_dp,  0.0033870_dp /)
         c(4,4,0:4) = (/ -0.0375491_dp,  0.0229520_dp, -0.0160047_dp,  0.0050463_dp, -0.0005364_dp /)

         !Computing the temperature coefficient of gray gas j
         sum_b = 0._dp
         do ii=0,4
            sum_c = 0._dp
            do kk=0,4
               sum_c = sum_c + c(jj,ii,kk)*mr**(real(kk,dp))
            enddo
            sum_b = sum_b + sum_c*(ttmp/tref)**(real(ii,dp))
         enddo
         a_func_res = sum_b
         
      endif               
      
   endfunction bordbar_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for a single 
   !gray gas following the formulation of Guo et al (2015)
   !====================================================================
   real(dp) recursive function guo_wsgg_a(ttmp,jj) result (a_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: ttmp
      real(dp) :: c(1:4,1:4),sum_a,sum_b,tref
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + guo_wsgg_a(ttmp,ii)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif

      !-----------------------------------------------------------------
      !Mounting the c's array
      !-----------------------------------------------------------------
      c(1:4,1) = (/  0.00527774_dp,  0.0127396_dp,&
                     0.01060850_dp, -0.0286258_dp /)
      c(1:4,2) = (/  0.05877410_dp, -0.2107950_dp,&
                    -0.04988030_dp,  0.201901_dp  /)
      c(1:4,3) = (/ -0.00221317_dp,  0.6406980_dp,&
                    -0.07515800_dp, -0.563327_dp  /)
      c(1:4,4) = (/  0.03367980_dp, -0.2072980_dp,& 
                     0.51450900_dp,  0.659109_dp  /)

      !-----------------------------------------------------------------
      !Computing the temperature coefficient of gray gas j
      !-----------------------------------------------------------------
      sum_b = 0._dp; tref = 1500._dp
      do ii=1,4
         sum_b = sum_b + c(jj,ii)*((ttmp/tref)**real(4-ii,dp))
      enddo
      a_func_res = sum_b

   endfunction guo_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for a single 
   !gray gas following the formulation of Ge et al. (2017)
   !====================================================================
   real(dp) recursive function ge_wsgg_a(ttmp,mol_rat,jj) &
      result (a_func_res)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: ttmp,mol_rat
      real(dp) :: c_aux,c1(1:4,1:4),c2(1:4,1:4),c3(1:4,1:4),&
         sum_a,sum_b,tref
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + ge_wsgg_a(ttmp,mol_rat,ii)
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif

      !-----------------------------------------------------------------
      !Mounting the ck's arrays
      !-----------------------------------------------------------------
      c1(1,1:4) = (/  0.070113_dp, -0.082250_dp,&
                      0.039140_dp, -0.006610_dp /)
      c1(2,1:4) = (/  0.450440_dp, -0.520540_dp,&
                      0.259516_dp, -0.048550_dp /)
      c1(3,1:4) = (/  0.059645_dp,  1.099669_dp,&
                     -0.955080_dp,  0.221254_dp /)
      c1(4,1:4) = (/  0.419802_dp, -0.496880_dp,&
                      0.656425_dp, -0.166090_dp /)

      c2(1,1:4) = (/  0.254835_dp, -0.393990_dp,&
                      0.221961_dp, -0.042870_dp /)
      c2(2,1:4) = (/ -0.493390_dp,  2.391296_dp,&
                     -2.029220_dp,  0.489448_dp /)
      c2(3,1:4) = (/ -0.019610_dp, -0.945460_dp,&
                      1.064002_dp, -0.280390_dp /)
      c2(4,1:4) = (/  0.258166_dp, -1.051840_dp,&
                      0.743259_dp, -0.166190_dp /)

      c3(1,1:4) = (/ -0.057950_dp,  0.142473_dp,&
                     -0.104400_dp,  0.023894_dp /)
      c3(2,1:4) = (/  0.124878_dp, -0.892900_dp,&
                      0.770610_dp, -0.186270_dp /)
      c3(3,1:4) = (/  0.040412_dp,  0.296147_dp,&
                     -0.335840_dp,  0.086820_dp /)
      c3(4,1:4) = (/ -0.107340_dp,  0.454277_dp,&
                     -0.330370_dp,  0.075552_dp /)

      !-----------------------------------------------------------------
      !Computing the temperature coefficient of gray gas j
      !-----------------------------------------------------------------
      tref = 1200._dp; sum_b = 0._dp
      do ii=1,4
         c_aux = c1(jj,ii) + c2(jj,ii)*mol_rat + c3(jj,ii)*mol_rat**2._dp
         sum_b = sum_b + c_aux*((ttmp/tref)**real(ii-1,dp))
      enddo
      a_func_res = sum_b

   endfunction ge_wsgg_a

   !====================================================================
   !Function to compute the temperature coefficient for a single
   !gray gas following the formulation of Shan et al. (2018)
   !====================================================================
   real(dp) recursive function shan_wsgg_a(ttmp,tp,mol_rat,&
      jj,interpolation_spec) result (a_func_res)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      character(*) :: interpolation_spec
      character(80) :: msg
      integer,intent(in) :: jj
      integer :: ii
      real(dp),intent(in) :: ttmp,tp,mol_rat
      real(dp) :: c_aux,c0(1:4,0:4),c1(1:4,0:4),c2(1:4,0:4),tref,&
         sum_a,sum_b,x1,x2,y1,y2
      
      !-----------------------------------------------------------------
      !Special case for transparent windows
      !-----------------------------------------------------------------
      if ((jj.lt.1).or.(jj.gt.4)) then
         sum_a = 0._dp
         do ii=1,4
            sum_a = sum_a + shan_wsgg_a(ttmp,tp,mol_rat,ii,'none')
         enddo
         a_func_res = 1._dp - sum_a
         return
      endif

      selectcase(trim(interpolation_spec))
         case('none')
            !-----------------------------------------------------------
            !Mounting the c's array
            !-----------------------------------------------------------
            if (mol_rat.eq.0.125_dp) then
               if (tp.le.2._dp) then
                  c0(1,0:4) = (/ -1.735122558e+0_dp,  5.080530622e+0_dp, -5.125955165e+0_dp, &
                     1.946001892e+0_dp,  9.691470600e-2_dp /)
                  c0(2,0:4) = (/  9.065256600e-2_dp, -4.219255660e-1_dp,  6.707691410e-1_dp, &
                     -5.495152460e-1_dp,  3.315414210e-1_dp /)
                  c0(3,0:4) = (/  4.383250000e-4_dp,  1.680944520e-1_dp, -4.305723780e-1_dp, &
                     1.961745750e-1_dp,  1.605817450e-1_dp /)
                  c0(4,0:4) = (/  3.054463600e-2_dp,  6.285695400e-2_dp, -2.714205840e-1_dp, &
                      1.751395520e-1_dp,  2.109199800e-2_dp /)
                  c1(1,0:4) = (/  3.850902068e+0_dp, -1.206334535e+1_dp,  1.308073904e+1_dp, &
                     -5.480776293e+0_dp,  7.154054650e-1_dp /)
                  c1(2,0:4) = (/ -1.978296350e+0_dp,  6.781313752e+0_dp, -8.173162126e+0_dp, &
                      3.833825305e+0_dp, -4.655643420e-1_dp /)
                  c1(3,0:4) = (/  6.782227110e-1_dp, -2.345866871e+0_dp,  2.937262245e+0_dp, &
                     -1.503885041e+0_dp,  2.183161010e-1_dp /)
                  c1(4,0:4) = (/ -3.070343000e-1_dp,  8.641695280e-1_dp, -8.059927980e-1_dp, &
                      2.435556430e-1_dp,  1.397633100e-2_dp /)
                  c2(1,0:4) = (/ -1.165102140e+0_dp,  3.598783466e+0_dp, -3.838701366e+0_dp, &
                      1.610416686e+0_dp, -2.363320580e-1_dp /)
                  c2(2,0:4) = (/  7.512887240e-1_dp, -2.528108612e+0_dp,  2.969973660e+0_dp, &
                     -1.349620506e+0_dp,  1.705445320e-1_dp /)
                  c2(3,0:4) = (/ -3.689911060e-1_dp,  1.235267546e+0_dp, -1.472281114e+0_dp, &
                      6.937160820e-1_dp, -7.809705000e-2_dp /)
                  c2(4,0:4) = (/  1.418551520e-1_dp, -4.183546560e-1_dp,  4.295334800e-1_dp, &
                     -1.752349100e-1_dp,  2.283008600e-2_dp /)
               elseif (tp.le.4._dp) then
                  c0(1,0:4) = (/  1.174774809e+0_dp, -4.354172175e+0_dp,  5.188950281e+0_dp, &
                     -2.209437695e+0_dp,  5.540111190e-1_dp /)
                  c0(2,0:4) = (/ -8.732300940e-1_dp,  3.243279142e+0_dp, -3.940825688e+0_dp, &
                     1.503376939e+0_dp,  2.153985790e-1_dp /)
                  c0(3,0:4) = (/  6.769321930e-1_dp, -2.336147559e+0_dp,  2.758457272e+0_dp, &
                     -1.371826598e+0_dp,  3.834673460e-1_dp /)
                  c0(4,0:4) = (/ -4.136934350e-1_dp,  1.475026750e+0_dp, -1.829108199e+0_dp, &
                     8.578020500e-1_dp, -7.546576500e-2_dp /)
                  c1(1,0:4) = (/  2.057818835e-1_dp, -5.605291535e-1_dp,  7.255481620e-1_dp, &
                     -4.324113625e-1_dp,  6.090665350e-2_dp /)
                  c1(2,0:4) = (/ -2.338575940e-1_dp,  6.581085020e-1_dp, -8.371930645e-1_dp, &
                     5.690070295e-1_dp, -1.479133290e-1_dp /)
                  c1(3,0:4) = (/ -4.000178370e-1_dp,  1.389840662e+0_dp, -1.591283520e+0_dp, &
                     6.249676245e-1_dp, -3.370198250e-2_dp /)
                  c1(4,0:4) = (/  2.409388745e-1_dp, -8.181235640e-1_dp,  9.924706105e-1_dp, &
                     -5.238619460e-1_dp,  1.214680975e-1_dp /)
                  c2(1,0:4) = (/ -7.001638950e-2_dp,  2.060510675e-1_dp, -2.398322890e-1_dp, &
                     1.250941175e-1_dp, -2.335675550e-2_dp /)
                  c2(2,0:4) = (/  1.200400110e-1_dp, -3.828071640e-1_dp,  4.548878365e-1_dp, &
                     -2.304344145e-1_dp,  4.075473600e-2_dp /)
                  c2(3,0:4) = (/  1.005701000e-3_dp, -6.525717500e-3_dp, -5.265644000e-3_dp, &
                     2.129004250e-2_dp, -7.809408500e-3_dp /)
                  c2(4,0:4) = (/ -2.107191750e-2_dp,  6.974944100e-2_dp, -8.027632050e-2_dp, &
                     3.780826000e-2_dp, -6.776356500e-3_dp /)
               elseif (tp.le.6._dp) then
                  c0(1,0:4) = (/  1.800071297e+0_dp, -6.175275997e+0_dp,  7.371754937e+0_dp, &
                     -3.423063169e+0_dp,  8.205125090e-1_dp /)
                  c0(2,0:4) = (/ -2.783340948e+0_dp,  9.540784048e+0_dp, -1.155374041e+1_dp, &
                     5.333783871e+0_dp, -4.387648850e-1_dp /)
                  c0(3,0:4) = (/  1.573214305e+0_dp, -5.305114527e+0_dp,  6.420221700e+0_dp, &
                     -3.326028952e+0_dp,  7.158583820e-1_dp /)
                  c0(4,0:4) = (/ -5.180785330e-1_dp,  1.787681908e+0_dp, -2.166918793e+0_dp, &
                     1.009262706e+0_dp, -9.623209700e-2_dp /)
                  c1(1,0:4) = (/ -2.756448105e-1_dp,  8.512280240e-1_dp, -9.086338500e-1_dp, &
                     4.336373720e-1_dp, -1.231891120e-1_dp /)
                  c1(2,0:4) = (/  9.447244295e-1_dp, -3.189273463e+0_dp,  3.747513782e+0_dp, &
                     -1.699985396e+0_dp,  2.360590210e-1_dp /)
                  c1(3,0:4) = (/ -8.816535410e-1_dp,  2.977820712e+0_dp, -3.542008127e+0_dp, &
                     1.657856207e+0_dp, -2.028732295e-1_dp /)
                  c1(4,0:4) = (/  2.832596190e-1_dp, -9.446846755e-1_dp,  1.132371585e+0_dp, &
                     -5.901615100e-1_dp,  1.316689805e-1_dp /)
                  c2(1,0:4) = (/  1.125925350e-2_dp, -3.306923800e-2_dp,  3.228792300e-2_dp, &
                     -1.556647400e-2_dp,  6.010849000e-3_dp /)
                  c2(2,0:4) = (/ -5.522356650e-2_dp,  1.854442705e-1_dp, -2.154817050e-1_dp, &
                     9.741325850e-2_dp, -1.435313500e-2_dp /)
                  c2(3,0:4) = (/  6.539699500e-2_dp, -2.179602945e-1_dp,  2.535552310e-1_dp, &
                     -1.147944560e-1_dp,  1.370896350e-2_dp /)
                  c2(4,0:4) = (/ -2.512803500e-2_dp,  8.184877150e-2_dp, -9.413840200e-2_dp, &
                     4.491686000e-2_dp, -8.028681500e-3_dp /)
               elseif (tp.le.10._dp) then
                  c0(1,0:4) = (/  2.262824648e+0_dp, -7.631682442e+0_dp,  9.007340648e+0_dp, &
                     -4.061624548e+0_dp,  8.133642260e-1_dp /)
                  c0(2,0:4) = (/ -1.804339320e+0_dp,  6.207066814e+0_dp, -7.655158297e+0_dp, &
                     3.545988702e+0_dp, -1.584535520e-1_dp /)
                  c0(3,0:4) = (/ -2.729361230e-1_dp,  7.941796020e-1_dp, -5.802904080e-1_dp, &
                     -2.073040030e-1_dp,  3.388628090e-1_dp /)
                  c0(4,0:4) = (/  2.125043720e-1_dp, -5.654526080e-1_dp,  4.810434020e-1_dp, &
                     -1.994677350e-1_dp,  9.743007400e-2_dp /)
                  c1(1,0:4) = (/ -4.556264655e-1_dp,  1.417454420e+0_dp, -1.540539555e+0_dp, &
                     6.823673830e-1_dp, -1.254060900e-1_dp /)
                  c1(2,0:4) = (/  6.664823050e-1_dp, -2.223457898e+0_dp,  2.600137314e+0_dp, &
                     -1.166169427e+0_dp,  1.498578472e-1_dp /)
                  c1(3,0:4) = (/ -2.812309890e-1_dp,  9.954119580e-1_dp, -1.270425820e+0_dp, &
                     6.493560547e-1_dp, -8.195382250e-2_dp /)
                  c1(4,0:4) = (/  3.165156775e-2_dp, -1.394824507e-1_dp,  2.329435902e-1_dp, &
                     -1.834885925e-1_dp,  6.760615100e-2_dp /)
                  c2(1,0:4) = (/  2.840193625e-2_dp, -8.698456950e-2_dp,  9.217260412e-2_dp, &
                     -3.928365975e-2_dp,  6.578908750e-3_dp /)
                  c2(2,0:4) = (/ -3.604436875e-2_dp,  1.170782662e-1_dp, -1.325462412e-1_dp, &
                     5.810490737e-2_dp, -7.772698625e-3_dp /)
                  c2(3,0:4) = (/  1.660852600e-2_dp, -5.698367250e-2_dp,  6.941684950e-2_dp, &
                     -3.334234588e-2_dp,  4.027828250e-3_dp /)
                  c2(4,0:4) = (/ -3.487329375e-3_dp,  1.301324838e-2_dp, -1.778824163e-2_dp, &
                     1.071388600e-2_dp, -2.731048000e-3_dp /)
               elseif (tp.le.20._dp) then
                  c0(1,0:4) = (/  2.106149890e-1_dp, -1.166582095e+0_dp,  2.004424551e+0_dp, &
                     -1.068604819e+0_dp,  2.787529820e-1_dp /)
                  c0(2,0:4) = (/  8.407840400e-2_dp, -1.210015850e-1_dp, -3.713936060e-1_dp, &
                     4.813693770e-1_dp,  1.807259660e-1_dp /)
                  c0(3,0:4) = (/ -7.670632160e-1_dp,  2.617987680e+0_dp, -2.943634875e+0_dp, &
                     9.413713540e-1_dp,  3.062550040e-1_dp /)
                  c0(4,0:4) = (/  4.244838660e-1_dp, -1.435055924e+0_dp,  1.741123173e+0_dp, &
                     -1.012202532e+0_dp,  3.297726140e-1_dp /)
                  c1(1,0:4) = (/  4.449742950e-2_dp, -1.341929580e-1_dp,  1.169826815e-1_dp, &
                     -1.979829100e-2_dp, -6.348843500e-3_dp /)
                  c1(2,0:4) = (/  1.835499361e-1_dp, -6.450497940e-1_dp,  8.199874428e-1_dp, &
                     -4.125155734e-1_dp,  5.858829680e-2_dp /)
                  c1(3,0:4) = (/ -1.369392715e-1_dp,  4.866465364e-1_dp, -6.390364673e-1_dp, &
                     3.507629797e-1_dp, -6.278316250e-2_dp /)
                  c1(4,0:4) = (/ -1.232280940e-2_dp,  3.593659780e-2_dp, -1.573952770e-2_dp, &
                     -2.785550200e-2_dp,  2.499782860e-2_dp /)
                  c2(1,0:4) = (/ -1.088356660e-3_dp,  3.529164880e-3_dp, -3.550458580e-3_dp, &
                     1.002710360e-3_dp,  1.929654000e-5_dp /)
                  c2(2,0:4) = (/ -6.635309100e-3_dp,  2.251813984e-2_dp, -2.736890104e-2_dp, &
                     1.338571524e-2_dp, -2.037538760e-3_dp /)
                  c2(3,0:4) = (/  7.120625180e-3_dp, -2.434521112e-2_dp,  2.991135890e-2_dp, &
                     -1.496979194e-2_dp,  2.436840300e-3_dp /)
                  c2(4,0:4) = (/ -1.209686600e-3_dp,  4.167376680e-3_dp, -5.520727540e-3_dp, &
                     3.277924920e-3_dp, -7.936411600e-4_dp /)
               elseif (tp.le.30._dp) then
                  c0(1,0:4) = (/ -2.615673790e-1_dp,  8.228156500e-2_dp,  7.973469050e-1_dp, &
                     -4.957763230e-1_dp,  1.747398860e-1_dp /)
                  c0(2,0:4) = (/  2.508380620e+0_dp, -7.900418475e+0_dp,  8.594052674e+0_dp, &
                     -3.747412931e+0_dp,  8.283734940e-1_dp /)
                  c0(3,0:4) = (/ -3.575897160e+0_dp,  1.153986279e+1_dp, -1.321435518e+1_dp, &
                     5.811882354e+0_dp, -4.544166220e-1_dp /)
                  c0(4,0:4) = (/  1.357228718e+0_dp, -4.131319522e+0_dp,  4.658203015e+0_dp, &
                     -2.396772654e+0_dp,  5.795668940e-1_dp /)
                  c1(1,0:4) = (/  6.876590950e-2_dp, -1.915631090e-1_dp,  1.688316330e-1_dp, &
                     -4.967110620e-2_dp,  4.265841000e-4_dp /)
                  c1(2,0:4) = (/ -9.031462430e-2_dp,  2.406143237e-1_dp, -2.058903280e-1_dp, &
                     7.346660840e-2_dp, -1.642688600e-2_dp /)
                  c1(3,0:4) = (/  1.782746401e-1_dp, -5.279829744e-1_dp,  5.388217289e-1_dp, &
                     -2.091262271e-1_dp,  2.477042880e-2_dp /)
                  c1(4,0:4) = (/ -1.025942744e-1_dp,  3.037636941e-1_dp, -3.136801314e-1_dp, &
                     1.161623525e-1_dp, -7.819014000e-4_dp /)
                  c2(1,0:4) = (/ -1.121324740e-3_dp,  3.275513280e-3_dp, -3.125212040e-3_dp, &
                     1.064279880e-3_dp, -5.944210000e-5_dp /)
                  c2(2,0:4) = (/  9.971633800e-4_dp, -2.316523820e-3_dp,  1.511371800e-3_dp, &
                     -3.414380800e-4_dp,  9.410156000e-5_dp /)
                  c2(3,0:4) = (/ -1.617985540e-3_dp,  4.081576640e-3_dp, -3.304750140e-3_dp, &
                     8.483909000e-4_dp, -3.916020000e-5_dp /)
                  c2(4,0:4) = (/  9.720245200e-4_dp, -2.483319140e-3_dp,  2.083603040e-3_dp, &
                     -4.615425000e-4_dp, -1.291403600e-4_dp /)
               endif
            elseif (mol_rat.eq.0.25_dp) then
               if (tp.le.2._dp) then
                  c0(1,0:4) = (/ -3.502486507e+0_dp,  1.019417822e+1_dp, -1.011829757e+1_dp, &
                     3.899847524e+0_dp, -1.329567160e-1_dp /)
                  c0(2,0:4) = (/  1.043304609e+0_dp, -3.015865679e+0_dp,  3.032319635e+0_dp, &
                     -1.454377193e+0_dp,  5.064989540e-1_dp /)
                  c0(3,0:4) = (/ -4.765521750e-1_dp,  1.646393482e+0_dp, -2.017588060e+0_dp, &
                     9.209544210e-1_dp,  5.833325000e-3_dp /)
                  c0(4,0:4) = (/  1.447684160e-1_dp, -2.851927910e-1_dp,  8.021601900e-2_dp, &
                     3.867737400e-2_dp,  3.391765900e-2_dp /)
                  c1(1,0:4) = (/  6.064331778e+0_dp, -1.860904741e+1_dp,  1.973553983e+1_dp, &
                     -8.170297445e+0_dp,  9.824153130e-1_dp /)
                  c1(2,0:4) = (/ -3.000066738e+0_dp,  9.505174945e+0_dp, -1.073988929e+1_dp, &
                     5.058396691e+0_dp, -7.363398250e-1_dp /)
                  c1(3,0:4) = (/  1.046090132e+0_dp, -3.450386434e+0_dp,  4.112109282e+0_dp, &
                     -2.157159381e+0_dp,  4.739982290e-1_dp /)
                  c1(4,0:4) = (/ -3.989518410e-1_dp,  1.111209221e+0_dp, -9.740848610e-1_dp, &
                     2.347471340e-1_dp,  4.375512500e-2_dp /)
                  c2(1,0:4) = (/ -1.897774928e+0_dp,  5.789198448e+0_dp, -6.080963260e+0_dp, &
                     2.497085550e+0_dp, -3.134712100e-1_dp /)
                  c2(2,0:4) = (/  1.048111372e+0_dp, -3.319609574e+0_dp,  3.724120102e+0_dp, &
                     -1.713947402e+0_dp,  2.417904900e-1_dp /)
                  c2(3,0:4) = (/ -4.404642440e-1_dp,  1.464744372e+0_dp, -1.753974736e+0_dp, &
                     8.963297180e-1_dp, -1.680091140e-1_dp /)
                  c2(4,0:4) = (/  1.536864420e-1_dp, -4.534290740e-1_dp,  4.566455980e-1_dp, &
                     -1.803398440e-1_dp,  2.128399400e-2_dp /)
               elseif (tp.le.4._dp) then
                  c0(1,0:4) = (/ -1.291851864e+0_dp,  3.425695684e+0_dp, -3.048250836e+0_dp, &
                     1.030294444e+0_dp,  2.595333270e-1_dp /)
                  c0(2,0:4) = (/ -7.852438150e-1_dp,  3.551520309e+0_dp, -5.327250985e+0_dp, &
                     2.795665643e+0_dp, -1.939883550e-1_dp /)
                  c0(3,0:4) = (/  8.148779790e-1_dp, -2.981713633e+0_dp,  3.938067331e+0_dp, &
                     -2.261117783e+0_dp,  5.596218640e-1_dp /)
                  c0(4,0:4) = (/ -3.401872500e-1_dp,  1.222704318e+0_dp, -1.544584976e+0_dp, &
                     7.484493560e-1_dp, -7.164237200e-2_dp /)
                  c1(1,0:4) = (/  1.722102397e+0_dp, -5.405497313e+0_dp,  6.006012522e+0_dp, &
                     -2.619579555e+0_dp,  2.663822925e-1_dp /)
                  c1(2,0:4) = (/ -2.108701040e-1_dp,  1.933740290e-1_dp,  2.834855950e-1_dp, &
                     -2.235066150e-1_dp,  6.632782650e-2_dp /)
                  c1(3,0:4) = (/ -5.451055370e-1_dp,  2.060895535e+0_dp, -2.744070267e+0_dp, &
                     1.406221941e+0_dp, -1.468806215e-1_dp /)
                  c1(4,0:4) = (/  1.948612480e-1_dp, -7.069596685e-1_dp,  9.570391575e-1_dp, &
                     -5.978309990e-1_dp,  1.648046255e-1_dp /)
                  c2(1,0:4) = (/ -2.793188985e-1_dp,  8.795440360e-1_dp, -9.837112895e-1_dp, &
                     4.391148750e-1_dp, -5.357721050e-2_dp /)
                  c2(2,0:4) = (/  1.106501610e-1_dp, -3.055556130e-1_dp,  3.023253130e-1_dp, &
                     -1.355064580e-1_dp,  1.557849150e-2_dp /)
                  c2(3,0:4) = (/  3.227605200e-2_dp, -1.338698335e-1_dp,  1.852011905e-1_dp, &
                     -8.984289200e-2_dp,  3.983176500e-3_dp /)
                  c2(4,0:4) = (/ -2.198118600e-2_dp,  7.868109350e-2_dp, -1.027161625e-1_dp, &
                     5.850622700e-2_dp, -1.285074850e-2_dp /)
               elseif (tp.le.6._dp) then
                  c0(1,0:4) = (/  2.590983232e+0_dp, -8.674761112e+0_dp,  1.035999292e+1_dp, &
                     -4.907255998e+0_dp,  1.001810339e+0_dp /)
                  c0(2,0:4) = (/ -2.566324727e+0_dp,  8.691804773e+0_dp, -1.069355266e+1_dp, &
                     5.207936249e+0_dp, -4.735765830e-1_dp /)
                  c0(3,0:4) = (/  7.108171350e-1_dp, -2.327104121e+0_dp,  2.888663037e+0_dp, &
                     -1.772129583e+0_dp,  5.901654140e-1_dp /)
                  c0(4,0:4) = (/ -1.261581220e-1_dp,  4.661864720e-1_dp, -5.852972940e-1_dp, &
                     2.331816800e-1_dp,  2.668785200e-2_dp /)
                  c1(1,0:4) = (/ -5.101889085e-1_dp,  1.570180146e+0_dp, -1.737030073e+0_dp, &
                     8.240222055e-1_dp, -1.757601565e-1_dp /)
                  c1(2,0:4) = (/  9.212471680e-1_dp, -3.127513191e+0_dp,  3.801892279e+0_dp, &
                     -1.823426973e+0_dp,  2.670176515e-1_dp /)
                  c1(3,0:4) = (/ -5.699565300e-1_dp,  1.976911377e+0_dp, -2.499612539e+0_dp, &
                     1.304896567e+0_dp, -1.913324790e-1_dp /)
                  c1(4,0:4) = (/  1.003928780e-1_dp, -3.683899330e-1_dp,  5.225410790e-1_dp, &
                     -3.622648640e-1_dp,  1.196281335e-1_dp /)
                  c2(1,0:4) = (/  3.607673450e-2_dp, -1.080967790e-1_dp,  1.140341245e-1_dp, &
                     -5.068866250e-2_dp,  1.056608850e-2_dp /)
                  c2(2,0:4) = (/ -6.106160000e-2_dp,  2.033984130e-1_dp, -2.418825035e-1_dp, &
                     1.137067185e-1_dp, -1.711970050e-2_dp /)
                  c2(3,0:4) = (/  4.499260300e-2_dp, -1.537868885e-1_dp,  1.896745270e-1_dp, &
                     -9.507331100e-2_dp,  1.318716900e-2_dp /)
                  c2(4,0:4) = (/ -1.174091400e-2_dp,  4.132102500e-2_dp, -5.404712300e-2_dp, &
                     3.181892300e-2_dp, -7.702264500e-3_dp /)
               elseif (tp.le.10._dp) then
                  c0(1,0:4) = (/ -1.302983840e-1_dp, -1.335419893e+0_dp,  3.686421361e+0_dp, &
                     -2.419361536e+0_dp,  5.959311230e-1_dp /)
                  c0(2,0:4) = (/  8.998249360e-1_dp, -1.920510739e+0_dp,  8.167814290e-1_dp, &
                     2.481931400e-1_dp,  2.816900310e-1_dp /)
                  c0(3,0:4) = (/ -1.674243273e+0_dp,  5.332714522e+0_dp, -5.918545035e+0_dp, &
                     2.322103836e+0_dp,  3.411125900e-2_dp /)
                  c0(4,0:4) = (/  5.113994840e-1_dp, -1.643627230e+0_dp,  1.956637452e+0_dp, &
                     -1.072392265e+0_dp,  2.617423940e-1_dp /)
                  c1(1,0:4) = (/  2.297524087e-1_dp, -4.291798053e-1_dp,  8.340090400e-2_dp, &
                     1.443104520e-1_dp, -6.154530825e-2_dp /)
                  c1(2,0:4) = (/ -8.375163400e-2_dp, -3.401000700e-2_dp,  4.227547155e-1_dp, &
                     -3.529373250e-1_dp,  4.061255625e-2_dp /)
                  c1(3,0:4) = (/  1.672501485e-1_dp, -4.094783020e-1_dp,  2.670316205e-1_dp, &
                     8.693203250e-3_dp, -1.251905175e-2_dp /)
                  c1(4,0:4) = (/ -9.861322975e-2_dp,  2.944496310e-1_dp, -2.817258002e-1_dp, &
                     5.367108425e-2_dp,  4.426931950e-2_dp /)
                  c2(1,0:4) = (/ -1.165566238e-2_dp,  2.125929013e-2_dp, -3.994050500e-3_dp, &
                     -6.511549750e-3_dp,  2.804703125e-3_dp /)
                  c2(2,0:4) = (/  1.015626525e-2_dp, -1.739890900e-2_dp,  1.575588000e-3_dp, &
                     6.395752500e-3_dp, -3.651461250e-4_dp /)
                  c2(3,0:4) = (/ -1.162349875e-2_dp,  3.117198450e-2_dp, -2.678816425e-2_dp, &
                     7.231876875e-3_dp, -1.169120125e-3_dp /)
                  c2(4,0:4) = (/  3.716837125e-3_dp, -1.054629950e-2_dp,  9.388058375e-3_dp, &
                     -1.237792125e-3_dp, -1.671755000e-3_dp /)
               elseif (tp.le.20._dp) then
                  c0(1,0:4) = (/  9.069585790e-1_dp, -3.575065054e+0_dp,  4.811424145e+0_dp, &
                     -2.326662017e+0_dp,  4.652235370e-1_dp /)
                  c0(2,0:4) = (/  3.674748280e-1_dp, -1.135041672e+0_dp,  1.137203301e+0_dp, &
                     -4.944245850e-1_dp,  4.005858940e-1_dp /)
                  c0(3,0:4) = (/ -2.223807053e+0_dp,  7.292974644e+0_dp, -8.532992502e+0_dp, &
                     3.894901531e+0_dp, -2.518255460e-1_dp /)
                  c0(4,0:4) = (/  1.330598684e+0_dp, -4.161710142e+0_dp,  4.753010915e+0_dp, &
                     -2.508414561e+0_dp,  6.007174220e-1_dp /)
                  c1(1,0:4) = (/  1.854464690e-2_dp, -1.227298330e-2_dp, -6.180553220e-2_dp, &
                     7.729983620e-2_dp, -2.445311680e-2_dp /)
                  c1(2,0:4) = (/  1.008308087e-1_dp, -4.050656223e-1_dp,  5.719828485e-1_dp, &
                     -3.032398685e-1_dp,  3.840011970e-2_dp /)
                  c1(3,0:4) = (/  1.346446922e-1_dp, -3.467448860e-1_dp,  2.670523541e-1_dp, &
                     -5.301659810e-2_dp, -3.105552300e-3_dp /)
                  c1(4,0:4) = (/ -1.992046443e-1_dp,  5.962387736e-1_dp, -6.103457268e-1_dp, &
                     2.315798974e-1_dp, -7.236502700e-3_dp /)
                  c2(1,0:4) = (/ -9.074558200e-4_dp,  1.965059540e-3_dp, -7.234347200e-4_dp, &
                     -7.374833600e-4_dp,  4.025598400e-4_dp /)
                  c2(2,0:4) = (/ -2.978477940e-3_dp,  1.185196186e-2_dp, -1.655144402e-2_dp, &
                     8.852184100e-3_dp, -1.332861100e-3_dp /)
                  c2(3,0:4) = (/ -2.867315320e-3_dp,  5.296041680e-3_dp, -6.457629400e-4_dp, &
                     -2.325119940e-3_dp,  7.488979800e-4_dp /)
                  c2(4,0:4) = (/  5.583986580e-3_dp, -1.554438464e-2_dp,  1.428631640e-2_dp, &
                     -4.668450480e-3_dp,  8.907694000e-5_dp /)
               elseif (tp.le.30._dp) then
                  c0(1,0:4) = (/  9.548568100e-2_dp, -1.065999992e+0_dp,  1.855112595e+0_dp, &
                     -7.641756950e-1_dp,  1.449039090e-1_dp /)
                  c0(2,0:4) = (/  3.612973652e+0_dp, -1.157249382e+1_dp,  1.322977935e+1_dp, &
                     -6.216697403e+0_dp,  1.260250042e+0_dp /)
                  c0(3,0:4) = (/ -3.706726049e+0_dp,  1.236107518e+1_dp, -1.523713272e+1_dp, &
                     7.667090737e+0_dp, -9.395657700e-1_dp /)
                  c0(4,0:4) = (/  3.263972320e-1_dp, -1.326650302e+0_dp,  2.444439431e+0_dp, &
                     -2.082602037e+0_dp,  7.265590180e-1_dp /)
                  c1(1,0:4) = (/  8.753245460e-2_dp, -2.323543744e-1_dp,  2.067494481e-1_dp, &
                     -7.027488350e-2_dp,  6.908948600e-3_dp /)
                  c1(2,0:4) = (/ -1.830147733e-1_dp,  5.299094673e-1_dp, -5.373601921e-1_dp, &
                     2.316647168e-1_dp, -4.309028610e-2_dp /)
                  c1(3,0:4) = (/  2.147192052e-1_dp, -6.814264819e-1_dp,  7.920649339e-1_dp, &
                     -3.814025540e-1_dp,  6.065277730e-2_dp /)
                  c1(4,0:4) = (/ -5.791944290e-2_dp,  2.096362488e-1_dp, -2.923884942e-1_dp, &
                     1.599908592e-1_dp, -1.689690330e-2_dp /)
                  c2(1,0:4) = (/ -2.328163960e-3_dp,  6.696466440e-3_dp, -6.760404860e-3_dp, &
                     2.735036820e-3_dp, -3.647443600e-4_dp /)
                  c2(2,0:4) = (/  3.100054100e-3_dp, -8.803162260e-3_dp,  8.684267900e-3_dp, &
                     -3.587363120e-3_dp,  5.924988200e-4_dp /)
                  c2(3,0:4) = (/ -3.163743480e-3_dp,  9.359870140e-3_dp, -1.013604138e-2_dp, &
                     4.663704840e-3_dp, -7.196679400e-4_dp /)
                  c2(4,0:4) = (/  1.030230140e-3_dp, -3.301908000e-3_dp,  4.159883480e-3_dp, &
                     -2.153529880e-3_dp,  2.574929800e-4_dp /)
               endif
            elseif (mol_rat.eq.0.5_dp) then
               if (tp.le.2._dp) then
                  c0(1,0:4) = (/  4.957454130e-1_dp, -1.867044609e+0_dp,  2.168912647e+0_dp, &
                     -8.630056100e-1_dp,  4.859790230e-1_dp /)
                  c0(2,0:4) = (/ -1.796530149e+0_dp,  6.067914554e+0_dp, -7.149705326e+0_dp, &
                     3.146113604e+0_dp, -1.212621610e-1_dp /)
                  c0(3,0:4) = (/  1.453941704e+0_dp, -4.841806424e+0_dp,  5.750086165e+0_dp, &
                     -2.894957554e+0_dp,  5.683816020e-1_dp /)
                  c0(4,0:4) = (/ -2.904748340e-1_dp,  1.099037011e+0_dp, -1.477501611e+0_dp, &
                     7.747138100e-1_dp, -9.428591000e-2_dp /)
                  c1(1,0:4) = (/  2.121686280e-1_dp, -9.866596400e-1_dp,  1.900306802e+0_dp, &
                     -1.359462281e+0_dp,  1.231337550e-1_dp /)
                  c1(2,0:4) = (/  1.102518249e+0_dp, -3.436988220e+0_dp,  3.381831676e+0_dp, &
                     -9.695334830e-1_dp, -1.107827800e-2_dp /)
                  c1(3,0:4) = (/ -1.618862614e+0_dp,  5.464112511e+0_dp, -6.476488234e+0_dp, &
                     2.949540130e+0_dp, -2.234304080e-1_dp /)
                  c1(4,0:4) = (/  1.926686020e-1_dp, -8.134822010e-1_dp,  1.276439675e+0_dp, &
                     -8.989939380e-1_dp,  2.582713720e-1_dp /)
                  c2(1,0:4) = (/  7.210185200e-2_dp, -1.769387760e-1_dp,  3.702676000e-3_dp, &
                     1.453660300e-1_dp, -1.815193000e-2_dp /)
                  c2(2,0:4) = (/ -2.661180580e-1_dp,  7.881419280e-1_dp, -6.828004120e-1_dp, &
                     1.189458100e-1_dp,  3.198391600e-2_dp /)
                  c2(3,0:4) = (/  3.612951600e-1_dp, -1.193967318e+0_dp,  1.347066788e+0_dp, &
                     -5.492020240e-1_dp,  1.626363200e-2_dp /)
                  c2(4,0:4) = (/ -3.817637600e-2_dp,  1.659561060e-1_dp, -2.520249180e-1_dp, &
                     1.579479360e-1_dp, -3.501642000e-2_dp /)
               elseif (tp.le.4._dp) then
                  c0(1,0:4) = (/ -5.840553490e-1_dp,  7.139707970e-1_dp,  7.815621190e-1_dp, &
                     -1.141382811e+0_dp,  6.223972080e-1_dp /)
                  c0(2,0:4) = (/ -7.293710220e-1_dp,  3.814530105e+0_dp, -6.539187441e+0_dp, &
                     4.120598945e+0_dp, -5.335241220e-1_dp /)
                  c0(3,0:4) = (/  4.560403210e-1_dp, -2.003496167e+0_dp,  3.225363322e+0_dp, &
                     -2.337396937e+0_dp,  7.300650470e-1_dp /)
                  c0(4,0:4) = (/ -1.712332760e-1_dp,  6.846855630e-1_dp, -9.591471350e-1_dp, &
                     5.096686540e-1_dp, -4.890752200e-2_dp /)
                  c1(1,0:4) = (/  1.261179223e+0_dp, -3.730649005e+0_dp,  3.747861938e+0_dp, &
                     -1.391559555e+0_dp,  5.666090350e-2_dp /)
                  c1(2,0:4) = (/ -6.698147750e-2_dp, -5.609590045e-1_dp,  1.738476516e+0_dp, &
                     -1.343570303e+0_dp,  3.029815455e-1_dp /)
                  c1(3,0:4) = (/ -5.178551515e-1_dp,  2.146710819e+0_dp, -3.244253663e+0_dp, &
                     2.007261478e+0_dp, -3.418903795e-1_dp /)
                  c1(4,0:4) = (/  9.474098500e-2_dp, -4.276754150e-1_dp,  7.485861390e-1_dp, &
                     -6.132022060e-1_dp,  2.076505200e-1_dp /)
                  c2(1,0:4) = (/ -1.824532550e-1_dp,  5.498020550e-1_dp, -5.732372600e-1_dp, &
                     2.310089675e-1_dp, -1.902005050e-2_dp /)
                  c2(2,0:4) = (/  5.184202350e-2_dp, -8.652656750e-2_dp, -1.375230350e-2_dp, &
                     6.234288450e-2_dp, -2.198050550e-2_dp /)
                  c2(3,0:4) = (/  6.026677450e-2_dp, -2.448440365e-1_dp,  3.621302135e-1_dp, &
                     -2.174528525e-1_dp,  3.507275650e-2_dp /)
                  c2(4,0:4) = (/ -1.902295700e-2_dp,  7.664057500e-2_dp, -1.176867690e-1_dp, &
                     8.131335900e-2_dp, -2.105059100e-2_dp /)
               elseif (tp.le.6._dp) then
                  c0(1,0:4) = (/  2.287978977e+0_dp, -8.120156455e+0_dp,  1.031717796e+1_dp, &
                     -5.248748127e+0_dp,  1.085109552e+0_dp /)
                  c0(2,0:4) = (/ -2.175444392e+0_dp,  7.618304285e+0_dp, -9.737564919e+0_dp, &
                     5.055667613e+0_dp, -5.409624520e-1_dp /)
                  c0(3,0:4) = (/ -7.230738500e-2_dp,  2.500296950e-1_dp, -1.810497880e-1_dp, &
                     -3.050546570e-1_dp,  4.285935070e-1_dp /)
                  c0(4,0:4) = (/  6.750656800e-2_dp, -2.568304990e-1_dp,  4.223705030e-1_dp, &
                     -3.732257260e-1_dp,  1.492986000e-1_dp /)
                  c1(1,0:4) = (/ -2.650100885e-1_dp,  9.244847680e-1_dp, -1.223525872e+0_dp, &
                     7.198443535e-1_dp, -1.744146345e-1_dp /)
                  c1(2,0:4) = (/  6.317042950e-1_dp, -2.354068402e+0_dp,  3.167611044e+0_dp, &
                     -1.704408402e+0_dp,  2.838470660e-1_dp /)
                  c1(3,0:4) = (/ -1.789946950e-1_dp,  7.828564920e-1_dp, -1.265334252e+0_dp, &
                     8.560684605e-1_dp, -1.715202665e-1_dp /)
                  c1(4,0:4) = (/ -4.302132000e-2_dp,  1.007569585e-1_dp, -9.010188500e-3_dp, &
                     -1.377271470e-1_dp,  1.022894115e-1_dp /)
                  c2(1,0:4) = (/  1.959192750e-2_dp, -6.184843500e-2_dp,  7.363370250e-2_dp, &
                     -4.013167750e-2_dp,  9.829312500e-3_dp /)
                  c2(2,0:4) = (/ -3.244983400e-2_dp,  1.240148955e-1_dp, -1.711373430e-1_dp, &
                     9.411061750e-2_dp, -1.673199000e-2_dp /)
                  c2(3,0:4) = (/  8.573392000e-3_dp, -4.472582100e-2_dp,  8.030118000e-2_dp, &
                     -5.667599050e-2_dp,  1.132219950e-2_dp /)
                  c2(4,0:4) = (/  4.963790000e-4_dp,  3.377235500e-3_dp, -1.463253950e-2_dp, &
                     1.762549300e-2_dp, -7.098196500e-3_dp /)
               elseif (tp.le.10._dp) then
                  c0(1,0:4) = (/  2.280069366e+0_dp, -8.126285122e+0_dp,  1.031045446e+1_dp, &
                     -5.156571303e+0_dp,  9.818042280e-1_dp /)
                  c0(2,0:4) = (/ -2.776749308e+0_dp,  8.924126078e+0_dp, -1.043339224e+1_dp, &
                     5.083687178e+0_dp, -5.121760810e-1_dp /)
                  c0(3,0:4) = (/  8.724882010e-1_dp, -2.271947155e+0_dp,  1.940471414e+0_dp, &
                     -9.017155310e-1_dp,  5.164852150e-1_dp /)
                  c0(4,0:4) = (/ -1.089739640e-1_dp,  1.281598820e-1_dp,  3.133801190e-1_dp, &
                     -5.921205520e-1_dp,  2.657534610e-1_dp /)
                  c1(1,0:4) = (/ -1.756578077e-1_dp,  6.729084017e-1_dp, -9.602296487e-1_dp, &
                     5.799975635e-1_dp, -1.252153622e-1_dp /)
                  c1(2,0:4) = (/  7.513705042e-1_dp, -2.555663966e+0_dp,  3.159642079e+0_dp, &
                     -1.610590400e+0_dp,  2.587307387e-1_dp /)
                  c1(3,0:4) = (/ -4.463437972e-1_dp,  1.481354351e+0_dp, -1.821228055e+0_dp, &
                     9.869875772e-1_dp, -1.896657295e-1_dp /)
                  c1(4,0:4) = (/  7.625320500e-3_dp, -5.586986750e-3_dp,  7.641631250e-3_dp, &
                     -5.832752575e-2_dp,  6.300749700e-2_dp /)
                  c2(1,0:4) = (/  4.919592125e-3_dp, -1.974879988e-2_dp,  2.993776262e-2_dp, &
                     -1.938434650e-2_dp,  4.499026125e-3_dp /)
                  c2(2,0:4) = (/ -3.569128788e-2_dp,  1.213413287e-1_dp, -1.504806456e-1_dp, &
                     7.769596262e-2_dp, -1.334555688e-2_dp /)
                  c2(3,0:4) = (/  2.688725388e-2_dp, -9.108721837e-2_dp,  1.140190026e-1_dp, &
                     -6.192193013e-2_dp,  1.190500700e-2_dp /)
                  c2(4,0:4) = (/ -3.042490750e-3_dp,  1.040704913e-2_dp, -1.438033213e-2_dp, &
                     1.047263463e-2_dp, -3.786068000e-3_dp /)
               elseif (tp.le.20._dp) then
                  c0(1,0:4) = (/  7.172544870e-1_dp, -3.092478768e+0_dp,  4.360267975e+0_dp, &
                     -2.104196946e+0_dp,  3.912080220e-1_dp /)
                  c0(2,0:4) = (/  1.039409101e+0_dp, -3.306378970e+0_dp,  3.724474882e+0_dp, &
                     -1.730483647e+0_dp,  6.304385230e-1_dp /)
                  c0(3,0:4) = (/ -3.006622885e+0_dp,  9.636721619e+0_dp, -1.112955614e+1_dp, &
                     5.280086640e+0_dp, -5.605819330e-1_dp /)
                  c0(4,0:4) = (/  1.583595155e+0_dp, -4.846613208e+0_dp,  5.419557614e+0_dp, &
                     -2.934283646e+0_dp,  7.571373410e-1_dp /)
                  c1(1,0:4) = (/  6.724713620e-2_dp, -1.248919528e-1_dp,  1.112942480e-2_dp, &
                     6.414347040e-2_dp, -2.189147440e-2_dp /)
                  c1(2,0:4) = (/ -7.923041800e-3_dp, -9.169014460e-2_dp,  2.618977848e-1_dp, &
                     -1.882158498e-1_dp,  1.529575400e-2_dp /)
                  c1(3,0:4) = (/  3.508896907e-1_dp, -1.018218045e+0_dp,  1.000523135e+0_dp, &
                     -3.878736677e-1_dp,  5.541685070e-2_dp /)
                  c1(4,0:4) = (/ -3.180205331e-1_dp,  9.722680837e-1_dp, -1.031237201e+0_dp, &
                     4.374603581e-1_dp, -4.192363260e-2_dp /)
                  c2(1,0:4) = (/ -3.742753480e-3_dp,  9.693172040e-3_dp, -7.696279920e-3_dp, &
                     1.677319240e-3_dp,  7.259940000e-5_dp /)
                  c2(2,0:4) = (/  2.076482640e-3_dp, -2.751002960e-3_dp, -2.284887360e-3_dp, &
                     3.600215880e-3_dp, -4.282044400e-4_dp /)
                  c2(3,0:4) = (/ -1.404498406e-2_dp,  3.978333354e-2_dp, -3.745584078e-2_dp, &
                     1.374617266e-2_dp, -1.832579540e-3_dp /)
                  c2(4,0:4) = (/  1.259640342e-2_dp, -3.763072702e-2_dp,  3.844577610e-2_dp, &
                     -1.568452282e-2_dp,  1.793206160e-3_dp /)
               elseif (tp.le.30._dp) then
                  c0(1,0:4) = (/  1.427064935e+0_dp, -4.844598452e+0_dp,  5.595976783e+0_dp, &
                     -2.205812890e+0_dp,  3.194990120e-1_dp /)
                  c0(2,0:4) = (/  1.064968585e+0_dp, -4.712497206e+0_dp,  7.280725976e+0_dp, &
                     -4.422047783e+0_dp,  1.088719793e+0_dp /)
                  c0(3,0:4) = (/ -1.432353479e+0_dp,  6.228465747e+0_dp, -9.747258712e+0_dp, &
                     5.968191016e+0_dp, -8.339235670e-1_dp /)
                  c0(4,0:4) = (/ -2.844165700e-2_dp, -3.064334040e-1_dp,  1.248396584e+0_dp, &
                     -1.570586132e+0_dp,  7.276638290e-1_dp /)
                  c1(1,0:4) = (/ -5.031061420e-2_dp,  1.841710658e-1_dp, -2.447911100e-1_dp, &
                     1.260637340e-1_dp, -2.142069630e-2_dp /)
                  c1(2,0:4) = (/  4.953593760e-2_dp, -1.248775152e-1_dp,  8.069897530e-2_dp, &
                     8.603458600e-3_dp, -1.904155310e-2_dp /)
                  c1(3,0:4) = (/ -3.415986640e-2_dp,  1.626895500e-3_dp,  1.637354788e-1_dp, &
                     -1.605498449e-1_dp,  3.851654080e-2_dp /)
                  c1(4,0:4) = (/  2.927518830e-2_dp, -3.944191850e-2_dp, -3.984723820e-2_dp, &
                     5.948969960e-2_dp, -4.211589800e-3_dp /)
                  c2(1,0:4) = (/  3.606079200e-4_dp, -1.379679680e-3_dp,  2.010474800e-3_dp, &
                     -1.164654080e-3_dp,  2.283330200e-4_dp /)
                  c2(2,0:4) = (/ -8.603650400e-4_dp,  2.423661160e-3_dp, -2.115574620e-3_dp, &
                     4.881608000e-4_dp,  1.429577400e-4_dp /)
                  c2(3,0:4) = (/  1.271820280e-3_dp, -2.688273820e-3_dp,  9.277984400e-4_dp, &
                     6.597205800e-4_dp, -3.042099600e-4_dp /)
                  c2(4,0:4) = (/ -7.382906200e-4_dp,  1.604323580e-3_dp, -6.958194400e-4_dp, &
                     -1.952336800e-4_dp, -1.871220000e-5_dp /)
               endif
            elseif (mol_rat.eq.0.75_dp) then
               if (tp.le.2._dp) then
                  c0(1,0:4) = (/ -2.036176170e-1_dp,  3.638281070e-1_dp, -2.633099920e-1_dp, &
                     2.327527650e-1_dp,  3.054762650e-1_dp /)
                  c0(2,0:4) = (/ -1.227911341e+0_dp,  4.403336192e+0_dp, -5.571011018e+0_dp, &
                     2.620464287e+0_dp, -5.483454800e-2_dp /)
                  c0(3,0:4) = (/  1.370769812e+0_dp, -4.670575543e+0_dp,  5.731355735e+0_dp, &
                     -3.025640916e+0_dp,  6.266892880e-1_dp /)
                  c0(4,0:4) = (/ -2.896445200e-1_dp,  1.110638138e+0_dp, -1.503327308e+0_dp, &
                     7.923916430e-1_dp, -9.945274300e-2_dp /)
                  c1(1,0:4) = (/  1.075074122e+0_dp, -3.877877616e+0_dp,  5.291493432e+0_dp, &
                     -3.033539275e+0_dp,  4.039629480e-1_dp /)
                  c1(2,0:4) = (/  4.787196250e-1_dp, -1.600413924e+0_dp,  1.565443186e+0_dp, &
                     -2.272435680e-1_dp, -1.523600840e-1_dp /)
                  c1(3,0:4) = (/ -1.650060295e+0_dp,  5.722641236e+0_dp, -7.044800566e+0_dp, &
                     3.393261605e+0_dp, -3.076087990e-1_dp /)
                  c1(4,0:4) = (/  2.332788030e-1_dp, -9.821487640e-1_dp,  1.531563576e+0_dp, &
                     -1.073240047e+0_dp,  3.046350670e-1_dp /)
                  c2(1,0:4) = (/ -1.343749480e-1_dp,  5.497398920e-1_dp, -9.023865080e-1_dp, &
                     6.248322780e-1_dp, -1.035341000e-1_dp /)
                  c2(2,0:4) = (/ -1.136429340e-1_dp,  3.120071000e-1_dp, -1.536356680e-1_dp, &
                     -1.453775360e-1_dp,  9.023886400e-2_dp /)
                  c2(3,0:4) = (/  4.020327660e-1_dp, -1.353282952e+0_dp,  1.568086648e+0_dp, &
                     -6.671778220e-1_dp,  2.810760600e-2_dp /)
                  c2(4,0:4) = (/ -7.274267400e-2_dp,  2.836432800e-1_dp, -3.938771800e-1_dp, &
                     2.294990340e-1_dp, -4.725464200e-2_dp /)
               elseif (tp.le.4._dp) then
                  c0(1,0:4) = (/ -4.461859060e-1_dp,  5.436276000e-1_dp,  6.361122520e-1_dp, &
                     -8.935862330e-1_dp,  5.303711400e-1_dp /)
                  c0(2,0:4) = (/ -4.503103560e-1_dp,  2.577169601e+0_dp, -4.724907908e+0_dp, &
                     3.124460312e+0_dp, -3.625218390e-1_dp /)
                  c0(3,0:4) = (/  1.046668770e-1_dp, -6.318515350e-1_dp,  1.352692792e+0_dp, &
                     -1.340262695e+0_dp,  5.975845090e-1_dp /)
                  c0(4,0:4) = (/ -7.866739300e-2_dp,  3.086939020e-1_dp, -4.229995480e-1_dp, &
                     1.991051480e-1_dp,  8.318436000e-3_dp /)
                  c1(1,0:4) = (/  1.318724949e+0_dp, -4.127972131e+0_dp,  4.463101086e+0_dp, &
                     -1.866313160e+0_dp,  1.613538875e-1_dp /)
                  c1(2,0:4) = (/ -3.362418805e-1_dp,  4.684241905e-1_dp,  3.945445750e-1_dp, &
                     -6.471795715e-1_dp,  1.824325085e-1_dp /)
                  c1(3,0:4) = (/ -2.639760265e-1_dp,  1.245893332e+0_dp, -2.162567726e+0_dp, &
                     1.529461017e+0_dp, -2.943167505e-1_dp /)
                  c1(4,0:4) = (/ -2.196467500e-3_dp, -9.357646400e-2_dp,  3.613202160e-1_dp, &
                     -4.542719245e-1_dp,  1.984897665e-1_dp /)
                  c2(1,0:4) = (/ -1.955582895e-1_dp,  6.298372765e-1_dp, -7.130458960e-1_dp, &
                     3.228039700e-1_dp, -3.845328850e-2_dp /)
                  c2(2,0:4) = (/  9.943757250e-2_dp, -2.658703095e-1_dp,  2.202878600e-1_dp, &
                     -6.140854050e-2_dp, -2.356095000e-4_dp /)
                  c2(3,0:4) = (/  2.551636550e-2_dp, -1.245900020e-1_dp,  2.216359635e-1_dp, &
                     -1.566220835e-1_dp,  2.873777650e-2_dp /)
                  c2(4,0:4) = (/ -7.749320500e-3_dp,  3.984318900e-2_dp, -7.883744000e-2_dp, &
                     6.833659650e-2_dp, -2.112478650e-2_dp /)
               elseif (tp.le.6._dp) then
                  c0(1,0:4) = (/  3.088967402e+0_dp, -1.041470324e+1_dp,  1.270891572e+1_dp, &
                     -6.404938329e+0_dp,  1.180872848e+0_dp /)
                  c0(2,0:4) = (/ -1.357433308e+0_dp,  4.759582107e+0_dp, -6.386804662e+0_dp, &
                     3.811620702e+0_dp, -4.925820650e-1_dp /)
                  c0(3,0:4) = (/ -1.299569133e+0_dp,  3.510626021e+0_dp, -2.956901614e+0_dp, &
                     4.252437330e-1_dp,  5.765850330e-1_dp /)
                  c0(4,0:4) = (/  2.774296790e-1_dp, -7.095758120e-1_dp,  6.296422220e-1_dp, &
                     -3.336576920e-1_dp,  1.540018280e-1_dp /)
                  c1(1,0:4) = (/ -4.789507215e-1_dp,  1.516140858e+0_dp, -1.836812637e+0_dp, &
                     1.037307880e+0_dp, -1.963160595e-1_dp /)
                  c1(2,0:4) = (/  2.687641375e-1_dp, -1.138214290e+0_dp,  1.819563890e+0_dp, &
                     -1.238470683e+0_dp,  2.699086190e-1_dp /)
                  c1(3,0:4) = (/  3.668533940e-1_dp, -6.811994010e-1_dp, -3.463016200e-2_dp, &
                     5.730229225e-1_dp, -2.437996415e-1_dp /)
                  c1(4,0:4) = (/ -1.788025475e-1_dp,  4.436606665e-1_dp, -2.523358045e-1_dp, &
                     -1.131694145e-1_dp,  1.105946825e-1_dp /)
                  c2(1,0:4) = (/  3.291354650e-2_dp, -9.629529350e-2_dp,  1.073823180e-1_dp, &
                     -5.864178400e-2_dp,  1.030784150e-2_dp /)
                  c2(2,0:4) = (/  4.881252500e-3_dp, -6.114710000e-4_dp, -3.209842150e-2_dp, &
                     4.346671300e-2_dp, -1.397587300e-2_dp /)
                  c2(3,0:4) = (/ -4.442623900e-2_dp,  9.827833400e-2_dp, -4.099877700e-2_dp, &
                     -2.785671150e-2_dp,  1.742096650e-2_dp /)
                  c2(4,0:4) = (/  1.414613250e-2_dp, -3.082423650e-2_dp,  8.786454500e-3_dp, &
                     1.635864650e-2_dp, -8.256227500e-3_dp /)
               elseif (tp.le.10._dp) then
                  c0(1,0:4) = (/  5.213268050e-1_dp, -2.926887103e+0_dp,  4.702748964e+0_dp, &
                     -2.655326358e+0_dp,  6.027812180e-1_dp /)
                  c0(2,0:4) = (/ -1.184350252e+0_dp,  4.498019040e+0_dp, -5.995795069e+0_dp, &
                     3.257166924e+0_dp, -3.055356080e-1_dp /)
                  c0(3,0:4) = (/  8.406368220e-1_dp, -2.715431920e+0_dp,  2.825924075e+0_dp, &
                     -1.305728307e+0_dp,  6.041712390e-1_dp /)
                  c0(4,0:4) = (/ -3.505864030e-1_dp,  9.632536060e-1_dp, -6.010949110e-1_dp, &
                     -2.810872430e-1_dp,  2.678514410e-1_dp /)
                  c1(1,0:4) = (/  2.663980450e-1_dp, -6.403220515e-1_dp,  4.524026460e-1_dp, &
                     -4.029840550e-2_dp, -3.272240275e-2_dp /)
                  c1(2,0:4) = (/  3.059441602e-1_dp, -1.330167409e+0_dp,  1.978644304e+0_dp, &
                     -1.163997380e+0_dp,  2.158662782e-1_dp /)
                  c1(3,0:4) = (/ -3.577313240e-1_dp,  1.369466980e+0_dp, -1.861373156e+0_dp, &
                     1.063931472e+0_dp, -2.249090435e-1_dp /)
                  c1(4,0:4) = (/  2.544885425e-2_dp, -8.168991700e-2_dp,  1.006865058e-1_dp, &
                     -9.377288750e-2_dp,  6.605726825e-2_dp /)
                  c2(1,0:4) = (/ -1.998789800e-2_dp,  5.512029875e-2_dp, -5.176004150e-2_dp, &
                     1.680337550e-2_dp, -8.996671250e-4_dp /)
                  c2(2,0:4) = (/ -6.123280625e-3_dp,  3.864635625e-2_dp, -6.947320150e-2_dp, &
                     4.645598962e-2_dp, -1.016455113e-2_dp /)
                  c2(3,0:4) = (/  1.688771525e-2_dp, -7.055334225e-2_dp,  1.028243417e-1_dp, &
                     -6.159224650e-2_dp,  1.350625000e-2_dp /)
                  c2(4,0:4) = (/ -2.450876625e-3_dp,  1.026671025e-2_dp, -1.586345463e-2_dp, &
                     1.166560175e-2_dp, -3.995814375e-3_dp /)
               elseif (tp.le.20._dp) then
                  c0(1,0:4) = (/  3.022358992e+0_dp, -1.017853674e+1_dp,  1.197296449e+1_dp, &
                     -5.367844811e+0_dp,  8.740466640e-1_dp /)
                  c0(2,0:4) = (/ -7.005402620e-1_dp,  1.730747567e+0_dp, -1.128064219e+0_dp, &
                     1.184436600e-2_dp,  4.936155120e-1_dp /)
                  c0(3,0:4) = (/ -2.782412090e-1_dp,  1.569482919e+0_dp, -2.907214000e+0_dp, &
                     1.877570937e+0_dp, -1.549858500e-1_dp /)
                  c0(4,0:4) = (/ -5.107198620e-1_dp,  1.678877245e+0_dp, -1.800692423e+0_dp, &
                     3.752689190e-1_dp,  2.331113940e-1_dp /)
                  c1(1,0:4) = (/ -2.340754947e-1_dp,  8.110528699e-1_dp, -1.011823899e+0_dp, &
                     5.127747884e-1_dp, -8.941988260e-2_dp /)
                  c1(2,0:4) = (/  2.532031470e-1_dp, -8.758720168e-1_dp,  1.070007257e+0_dp, &
                     -5.175191603e-1_dp,  5.126662500e-2_dp /)
                  c1(3,0:4) = (/ -8.009074200e-2_dp,  2.689472688e-1_dp, -3.429177026e-1_dp, &
                     1.972682001e-1_dp, -2.545314460e-2_dp /)
                  c1(4,0:4) = (/ -2.068571900e-3_dp, -1.705488400e-3_dp,  3.550871210e-2_dp, &
                     -5.083604300e-2_dp,  3.939188960e-2_dp /)
                  c2(1,0:4) = (/  5.049134100e-3_dp, -1.750069706e-2_dp,  2.196045778e-2_dp, &
                     -1.137875936e-2_dp,  2.057426400e-3_dp /)
                  c2(2,0:4) = (/ -5.687279200e-3_dp,  2.088953176e-2_dp, -2.728680528e-2_dp, &
                     1.426139326e-2_dp, -1.696097000e-3_dp /)
                  c2(3,0:4) = (/  3.124373600e-4_dp, -3.350519520e-3_dp,  8.310177160e-3_dp, &
                     -6.758911700e-3_dp,  1.152231000e-3_dp /)
                  c2(4,0:4) = (/  1.902200580e-3_dp, -4.887969000e-3_dp,  2.650299860e-3_dp, &
                     8.083556800e-4_dp, -9.818760400e-4_dp /)
               elseif (tp.le.30._dp) then
                  c0(1,0:4) = (/  2.712437936e+0_dp, -8.924555694e+0_dp,  9.983880961e+0_dp, &
                     -3.966048071e+0_dp,  5.527217600e-1_dp /)
                  c0(2,0:4) = (/  1.213863640e-1_dp, -2.564678859e+0_dp,  5.882062259e+0_dp, &
                     -4.262546992e+0_dp,  1.114785748e+0_dp /)
                  c0(3,0:4) = (/ -1.345757265e+0_dp,  6.570709643e+0_dp, -1.077625801e+1_dp, &
                     6.676317627e+0_dp, -9.806843420e-1_dp /)
                  c0(4,0:4) = (/ -2.510596340e-1_dp,  4.352255690e-1_dp,  3.757755830e-1_dp, &
                     -1.169709457e+0_dp,  6.839876340e-1_dp /)
                  c1(1,0:4) = (/ -1.802842907e-1_dp,  6.084717566e-1_dp, -7.209982547e-1_dp, &
                     3.285773058e-1_dp, -4.980144700e-2_dp /)
                  c1(2,0:4) = (/  1.566372385e-1_dp, -3.969116591e-1_dp,  3.043737029e-1_dp, &
                     -5.502981440e-2_dp, -1.463963480e-2_dp /)
                  c1(3,0:4) = (/ -5.869835680e-2_dp,  1.990193580e-2_dp,  2.119915481e-1_dp, &
                     -2.064286012e-1_dp,  4.835757200e-2_dp /)
                  c1(4,0:4) = (/  5.656034150e-2_dp, -1.260707942e-1_dp,  5.189915460e-2_dp, &
                     2.229215180e-2_dp,  1.627852800e-3_dp /)
                  c2(1,0:4) = (/  3.134376540e-3_dp, -1.050659400e-2_dp,  1.239188438e-2_dp, &
                     -5.673377080e-3_dp,  8.798168800e-4_dp /)
                  c2(2,0:4) = (/ -2.913800340e-3_dp,  7.680079940e-3_dp, -6.530443780e-3_dp, &
                     1.822904360e-3_dp,  4.629040000e-5_dp /)
                  c2(3,0:4) = (/  1.911608240e-3_dp, -3.401319680e-3_dp,  2.373246600e-4_dp, &
                     1.429061640e-3_dp, -4.740586000e-4_dp /)
                  c2(4,0:4) = (/ -1.678395660e-3_dp,  4.439425480e-3_dp, -3.610392280e-3_dp, &
                     1.014391880e-3_dp, -2.208648000e-4_dp /)
               endif
            elseif (mol_rat.eq.1._dp) then
               if (tp.le.2._dp) then
                  c0(1,0:4) = (/  2.177563120e-1_dp, -4.435386560e-1_dp,  8.549412800e-2_dp, &
                     2.815893190e-1_dp,  2.182410770e-1_dp /)
                  c0(2,0:4) = (/ -4.542801710e-1_dp,  1.516783791e+0_dp, -1.920137248e+0_dp, &
                     8.818704880e-1_dp,  1.926097080e-1_dp /)
                  c0(3,0:4) = (/ -1.472897760e-1_dp,  1.947826850e-1_dp,  2.465187830e-1_dp, &
                     -5.764844610e-1_dp,  3.485357560e-1_dp /)
                  c0(4,0:4) = (/  4.914189000e-3_dp,  2.005788080e-1_dp, -4.993201110e-1_dp, &
                     3.401259100e-1_dp, -3.609668800e-2_dp /)
                  c1(1,0:4) = (/  3.700161210e-1_dp, -2.421659879e+0_dp,  4.521505523e+0_dp, &
                     -3.031216099e+0_dp,  5.190864470e-1_dp /)
                  c1(2,0:4) = (/ -5.800259640e-1_dp,  2.383235586e+0_dp, -3.545828770e+0_dp, &
                     2.293561238e+0_dp, -5.391689340e-1_dp /)
                  c1(3,0:4) = (/  5.566592110e-1_dp, -1.300978382e+0_dp,  7.966934480e-1_dp, &
                     -7.672881700e-2_dp,  9.653501400e-2_dp /)
                  c1(4,0:4) = (/ -2.188176180e-1_dp,  3.922541830e-1_dp,  6.766606500e-2_dp, &
                     -4.600643810e-1_dp,  2.329634610e-1_dp /)
                  c2(1,0:4) = (/  1.463884060e-1_dp, -8.312007400e-2_dp, -4.859813660e-1_dp, &
                     5.568403500e-1_dp, -1.323702780e-1_dp /)
                  c2(2,0:4) = (/  2.171186600e-1_dp, -9.647052000e-1_dp,  1.529733352e+0_dp, &
                     -1.002606020e+0_dp,  2.268459760e-1_dp /)
                  c2(3,0:4) = (/ -3.208093140e-1_dp,  9.587658640e-1_dp, -1.036536552e+0_dp, &
                     5.076341180e-1_dp, -1.158029000e-1_dp /)
                  c2(4,0:4) = (/  7.081147200e-2_dp, -1.525569460e-1_dp,  7.204694600e-2_dp, &
                     3.046821000e-2_dp, -2.140395000e-2_dp /)
               elseif (tp.le.4._dp) then
                  c0(1,0:4) = (/  2.395541310e-1_dp, -1.347389678e+0_dp,  2.406334374e+0_dp, &
                     -1.540539857e+0_dp,  5.932556470e-1_dp /)
                  c0(2,0:4) = (/ -6.838636250e-1_dp,  2.859544568e+0_dp, -4.486810116e+0_dp, &
                     2.781277460e+0_dp, -2.787039180e-1_dp /)
                  c0(3,0:4) = (/ -1.164767340e-1_dp,  4.285067420e-1_dp, -3.038126610e-1_dp, &
                     -3.516412760e-1_dp,  4.320905540e-1_dp /)
                  c0(4,0:4) = (/  2.339444600e-2_dp, -9.918080400e-2_dp,  1.534863730e-1_dp, &
                     -1.381925560e-1_dp,  7.398216100e-2_dp /)
                  c1(1,0:4) = (/  9.194568405e-1_dp, -3.083571614e+0_dp,  3.567226694e+0_dp, &
                     -1.589363583e+0_dp,  1.423237420e-1_dp /)
                  c1(2,0:4) = (/ -1.898515410e-1_dp,  2.775139385e-1_dp,  2.879162840e-1_dp, &
                     -4.532938380e-1_dp,  1.284061990e-1_dp /)
                  c1(3,0:4) = (/ -1.047106520e-1_dp,  5.416390005e-1_dp, -1.143489416e+0_dp, &
                     9.715977335e-1_dp, -2.056409830e-1_dp /)
                  c1(4,0:4) = (/ -9.663948750e-2_dp,  2.488673750e-1_dp, -6.880307500e-2_dp, &
                     -2.436707100e-1_dp,  1.699961695e-1_dp /)
                  c2(1,0:4) = (/ -1.337814085e-1_dp,  4.737985490e-1_dp, -5.890520130e-1_dp, &
                     2.914463860e-1_dp, -3.774256800e-2_dp /)
                  c2(2,0:4) = (/  7.942731200e-2_dp, -2.475345705e-1_dp,  2.545290420e-1_dp, &
                     -1.040302250e-1_dp,  1.088681600e-2_dp /)
                  c2(3,0:4) = (/  2.172357000e-3_dp, -2.097384150e-2_dp,  7.113774100e-2_dp, &
                     -7.273995350e-2_dp,  1.439639900e-2_dp /)
                  c2(4,0:4) = (/  5.102342500e-3_dp, -5.923639000e-3_dp, -2.292010500e-2_dp, &
                     4.185099100e-2_dp, -1.744001650e-2_dp /)
               elseif (tp.le.6._dp) then
                  c0(1,0:4) = (/  4.174270830e-1_dp, -2.875497350e+0_dp,  5.398362732e+0_dp, &
                     -3.523611609e+0_dp,  8.832882850e-1_dp /)
                  c0(2,0:4) = (/ -6.523797270e-1_dp,  3.035827478e+0_dp, -4.889508410e+0_dp, &
                     3.053939322e+0_dp, -2.862252140e-1_dp /)
                  c0(3,0:4) = (/ -1.201398232e+0_dp,  3.821318754e+0_dp, -4.371154983e+0_dp, &
                     1.819400304e+0_dp,  9.580552600e-2_dp /)
                  c0(4,0:4) = (/  3.113609840e-1_dp, -1.113036776e+0_dp,  1.593910447e+0_dp, &
                     -1.107109932e+0_dp,  3.130512150e-1_dp /)
                  c1(1,0:4) = (/  6.637746265e-1_dp, -1.733723964e+0_dp,  1.335731611e+0_dp, &
                     -2.169390690e-1_dp, -6.473365950e-2_dp /)
                  c1(2,0:4) = (/ -4.837839750e-2_dp, -3.817720870e-1_dp,  1.223858720e+0_dp, &
                     -9.820242415e-1_dp,  2.017111950e-1_dp /)
                  c1(3,0:4) = (/  3.669177285e-1_dp, -8.922033785e-1_dp,  5.492225905e-1_dp, &
                     7.122606250e-2_dp, -7.523623400e-2_dp /)
                  c1(4,0:4) = (/ -2.281127840e-1_dp,  7.105845480e-1_dp, -7.332810515e-1_dp, &
                     2.134610660e-1_dp,  5.497170000e-2_dp /)
                  c2(1,0:4) = (/ -8.097791450e-2_dp,  2.318433660e-1_dp, -2.181800145e-1_dp, &
                     7.228224200e-2_dp, -4.105257500e-3_dp /)
                  c2(2,0:4) = (/  4.209128250e-2_dp, -9.373074600e-2_dp,  4.571207650e-2_dp, &
                     1.111100950e-2_dp, -6.969352000e-3_dp /)
                  c2(3,0:4) = (/ -4.792714450e-2_dp,  1.254360025e-1_dp, -9.783136550e-2_dp, &
                     1.666286550e-2_dp,  2.813026000e-3_dp /)
                  c2(4,0:4) = (/  1.997275800e-2_dp, -5.798693400e-2_dp,  5.317288450e-2_dp, &
                     -1.187461700e-2_dp, -3.625715000e-3_dp /)
               elseif (tp.le.10._dp) then
                  c0(1,0:4) = (/  5.124685923e+0_dp, -1.664114300e+1_dp,  1.915680136e+1_dp, &
                     -8.867395968e+0_dp,  1.503262732e+0_dp /)
                  c0(2,0:4) = (/ -4.414367589e+0_dp,  1.339362032e+1_dp, -1.417379461e+1_dp, &
                     6.034368984e+0_dp, -5.215394090e-1_dp /)
                  c0(3,0:4) = (/  6.005700650e-1_dp, -5.176854090e-1_dp, -1.633042971e+0_dp, &
                     1.865377023e+0_dp, -1.617710090e-1_dp /)
                  c0(4,0:4) = (/  6.811882460e-1_dp, -2.735582834e+0_dp,  4.204337395e+0_dp, &
                     -2.931242880e+0_dp,  7.736108490e-1_dp /)
                  c1(1,0:4) = (/ -9.065811942e-1_dp,  2.878295743e+0_dp, -3.291830418e+0_dp, &
                     1.586894975e+0_dp, -2.739144047e-1_dp /)
                  c1(2,0:4) = (/  1.218702322e+0_dp, -3.889011426e+0_dp,  4.402233589e+0_dp, &
                     -2.028810925e+0_dp,  2.905794010e-1_dp /)
                  c1(3,0:4) = (/ -2.942329165e-1_dp,  7.636566492e-1_dp, -6.298035492e-1_dp, &
                     1.993508885e-1_dp, -1.681240700e-2_dp /)
                  c1(4,0:4) = (/ -2.921167767e-1_dp,  1.056976776e+0_dp, -1.374184225e+0_dp, &
                     7.074551325e-1_dp, -7.850087950e-2_dp /)
                  c2(1,0:4) = (/  4.999086562e-2_dp, -1.544475392e-1_dp,  1.709014730e-1_dp, &
                     -7.991831087e-2_dp,  1.353668763e-2_dp /)
                  c2(2,0:4) = (/ -6.458917462e-2_dp,  2.030926761e-1_dp, -2.261202295e-1_dp, &
                     1.027857439e-1_dp, -1.524421425e-2_dp /)
                  c2(3,0:4) = (/  1.220995475e-2_dp, -3.001277538e-2_dp,  2.261432413e-2_dp, &
                     -5.968403250e-3_dp,  2.306252500e-4_dp /)
                  c2(4,0:4) = (/  2.036711063e-2_dp, -7.064824825e-2_dp,  8.747822050e-2_dp, &
                     -4.353660175e-2_dp,  5.826391750e-3_dp /)
               elseif (tp.le.20._dp) then
                  c0(1,0:4) = (/  2.896240981e+0_dp, -9.418863550e+0_dp,  1.068699196e+1_dp, &
                     -4.601604253e+0_dp,  7.132835640e-1_dp /)
                  c0(2,0:4) = (/ -1.920853879e+0_dp,  5.263105328e+0_dp, -4.654557670e+0_dp, &
                     1.443070224e+0_dp,  1.848608210e-1_dp /)
                  c0(3,0:4) = (/  1.084708436e+0_dp, -2.523105561e+0_dp,  1.420921519e+0_dp, &
                     8.437155400e-2_dp,  1.810358500e-1_dp /)
                  c0(4,0:4) = (/ -4.717117760e-1_dp,  1.340388609e+0_dp, -1.213600829e+0_dp, &
                     -2.316892700e-2_dp,  3.577856110e-1_dp /)
                  c1(1,0:4) = (/ -2.333791930e-1_dp,  7.734324952e-1_dp, -9.300680472e-1_dp, &
                     4.597083271e-1_dp, -7.662365630e-2_dp /)
                  c1(2,0:4) = (/  4.318880086e-1_dp, -1.410836646e+0_dp,  1.630650863e+0_dp, &
                     -7.601859608e-1_dp,  9.973124130e-2_dp /)
                  c1(3,0:4) = (/ -2.993089903e-1_dp,  9.234432465e-1_dp, -1.031080724e+0_dp, &
                     4.888351883e-1_dp, -7.918430760e-2_dp /)
                  c1(4,0:4) = (/  2.798472390e-2_dp, -5.604002770e-2_dp,  5.283909330e-2_dp, &
                     -3.637550370e-2_dp,  3.096358820e-2_dp /)
                  c2(1,0:4) = (/  4.955114920e-3_dp, -1.618400900e-2_dp,  1.942332984e-2_dp, &
                     -9.857563260e-3_dp,  1.707404460e-3_dp /)
                  c2(2,0:4) = (/ -1.084288036e-2_dp,  3.658034794e-2_dp, -4.415432634e-2_dp, &
                     2.183623508e-2_dp, -3.223400580e-3_dp /)
                  c2(3,0:4) = (/  7.876178420e-3_dp, -2.593723358e-2_dp,  3.220239668e-2_dp, &
                     -1.710677854e-2_dp,  3.039746720e-3_dp /)
                  c2(4,0:4) = (/ -1.140392200e-4_dp, -1.062822600e-4_dp, -1.044729140e-3_dp, &
                     1.765722340e-3_dp, -9.618026400e-4_dp /)
               elseif (tp.le.30._dp) then
                  c0(1,0:4) = (/  1.299635767e+0_dp, -4.292228618e+0_dp,  4.649518347e+0_dp, &
                     -1.568112201e+0_dp,  1.862952920e-1_dp /)
                  c0(2,0:4) = (/  2.280998461e+0_dp, -8.780702192e+0_dp,  1.205129407e+1_dp, &
                     -6.696123616e+0_dp,  1.402424647e+0_dp /)
                  c0(3,0:4) = (/ -3.088976014e+0_dp,  1.127255500e+1_dp, -1.504320798e+1_dp, &
                     8.301229018e+0_dp, -1.194361692e+0_dp /)
                  c0(4,0:4) = (/  4.405017380e-1_dp, -1.593647249e+0_dp,  2.438381437e+0_dp, &
                     -2.060257637e+0_dp,  8.650529890e-1_dp /)
                  c1(1,0:4) = (/ -6.604851030e-2_dp,  2.379199886e-1_dp, -3.004928138e-1_dp, &
                     1.423630485e-1_dp, -2.136959950e-2_dp /)
                  c1(2,0:4) = (/  1.563143160e-2_dp,  4.163121000e-4_dp, -7.846767950e-2_dp, &
                     8.992141840e-2_dp, -3.037039520e-2_dp /)
                  c1(3,0:4) = (/  6.514612140e-2_dp, -3.060990111e-1_dp,  5.003598657e-1_dp, &
                     -3.150234513e-1_dp,  6.202354670e-2_dp /)
                  c1(4,0:4) = (/ -1.237053060e-2_dp,  7.555877280e-2_dp, -1.562202276e-1_dp, &
                     1.148727658e-1_dp, -1.576311310e-2_dp /)
                  c2(1,0:4) = (/  5.800938200e-4_dp, -2.224971000e-3_dp,  3.038252200e-3_dp, &
                     -1.574029460e-3_dp,  2.621723000e-4_dp /)
                  c2(2,0:4) = (/ -5.346823600e-4_dp,  1.127218860e-3_dp, -4.630285400e-4_dp, &
                     -3.211492800e-4_dp,  2.377716800e-4_dp /)
                  c2(3,0:4) = (/  8.763396000e-5_dp,  1.050727900e-3_dp, -3.209309060e-3_dp, &
                     2.544009780e-3_dp, -5.821521400e-4_dp /)
                  c2(4,0:4) = (/ -3.768102800e-4_dp,  6.488673600e-4_dp,  2.782812400e-4_dp, &
                     -7.039693600e-4_dp,  1.063639800e-4_dp /)
               endif
            elseif (mol_rat.eq.2._dp) then
               if (tp.le.2._dp) then
                  c0(1,0:4) = (/ -3.230723556e+0_dp,  9.803733010e+0_dp, -1.031072653e+1_dp, &
                     4.398299147e+0_dp, -2.513631170e-1_dp /)
                  c0(2,0:4) = (/  2.072909162e+0_dp, -5.451248065e+0_dp,  4.144747873e+0_dp, &
                     -7.813727140e-1_dp,  1.724324980e-1_dp /)
                  c0(3,0:4) = (/ -3.599706330e-1_dp,  7.966057440e-1_dp, -5.347968500e-2_dp, &
                     -8.688409880e-1_dp,  5.134539480e-1_dp /)
                  c0(4,0:4) = (/ -4.575668800e-1_dp,  1.371797813e+0_dp, -1.444010654e+0_dp, &
                     6.323323120e-1_dp, -1.020519080e-1_dp /)
                  c1(1,0:4) = (/  4.952590649e+0_dp, -1.603976646e+1_dp,  1.838733328e+1_dp, &
                     -8.583966607e+0_dp,  1.141064195e+0_dp /)
                  c1(2,0:4) = (/ -3.552667663e+0_dp,  1.026320718e+1_dp, -9.902423135e+0_dp, &
                     3.745651140e+0_dp, -4.280028480e-1_dp /)
                  c1(3,0:4) = (/  4.917494160e-1_dp, -9.335032770e-1_dp, -2.724823780e-1_dp, &
                     1.103029518e+0_dp, -2.358122410e-1_dp /)
                  c1(4,0:4) = (/  2.737353210e-1_dp, -7.469425280e-1_dp,  8.716022410e-1_dp, &
                     -7.135694430e-1_dp,  3.394745390e-1_dp /)
                  c2(1,0:4) = (/ -1.243515114e+0_dp,  4.024870484e+0_dp, -4.646282902e+0_dp, &
                     2.213411654e+0_dp, -3.121890980e-1_dp /)
                  c2(2,0:4) = (/  1.022651866e+0_dp, -3.034288688e+0_dp,  3.092021234e+0_dp, &
                     -1.283581888e+0_dp,  1.707012120e-1_dp /)
                  c2(3,0:4) = (/ -1.696526760e-1_dp,  4.426778980e-1_dp, -2.802467720e-1_dp, &
                     -2.807982400e-2_dp,  3.524402000e-3_dp /)
                  c2(4,0:4) = (/ -9.198781800e-2_dp,  2.210034400e-1_dp, -1.858343460e-1_dp, &
                     1.045450420e-1_dp, -5.065836600e-2_dp /)
               elseif (tp.le.4._dp) then
                  c0(1,0:4) = (/  1.124742640e-1_dp, -1.376094169e+0_dp,  3.128227747e+0_dp, &
                     -2.308273297e+0_dp,  7.671192040e-1_dp /)
                  c0(2,0:4) = (/ -1.314679965e+0_dp,  5.156706239e+0_dp, -7.497644722e+0_dp, &
                     4.523253504e+0_dp, -6.283518970e-1_dp /)
                  c0(3,0:4) = (/  1.693283319e+0_dp, -5.371458409e+0_dp,  6.221918426e+0_dp, &
                     -3.409042933e+0_dp,  9.624701410e-1_dp /)
                  c0(4,0:4) = (/ -1.270532531e+0_dp,  3.996331278e+0_dp, -4.417351799e+0_dp, &
                     1.958206249e+0_dp, -2.521044510e-1_dp /)
                  c1(1,0:4) = (/  1.125891931e+0_dp, -3.494118478e+0_dp,  3.613773630e+0_dp, &
                     -1.344376897e+0_dp,  6.498868550e-2_dp /)
                  c1(2,0:4) = (/  1.035729855e-1_dp, -9.440471660e-1_dp,  2.096348852e+0_dp, &
                     -1.583956059e+0_dp,  3.574345605e-1_dp /)
                  c1(3,0:4) = (/ -1.156649576e+0_dp,  3.994202978e+0_dp, -5.191630272e+0_dp, &
                     3.002659825e+0_dp, -5.790028065e-1_dp /)
                  c1(4,0:4) = (/  6.868308995e-1_dp, -2.237985249e+0_dp,  2.741791606e+0_dp, &
                     -1.586168718e+0_dp,  4.042942995e-1_dp /)
                  c2(1,0:4) = (/ -1.659652100e-1_dp,  5.470032875e-1_dp, -6.192416480e-1_dp, &
                     2.702599100e-1_dp, -2.877192350e-2_dp /)
                  c2(2,0:4) = (/  4.142882350e-2_dp, -8.265009100e-2_dp,  3.233389500e-3_dp, &
                     5.506515700e-2_dp, -2.182139350e-2_dp /)
                  c2(3,0:4) = (/  1.412333320e-1_dp, -4.791591915e-1_dp,  6.104776475e-1_dp, &
                     -3.428444915e-1_dp,  6.286563650e-2_dp /)
                  c2(4,0:4) = (/ -9.529419450e-2_dp,  3.103914345e-1_dp, -3.775937425e-1_dp, &
                     2.093761955e-1_dp, -4.555511050e-2_dp /)
               elseif (tp.le.6._dp) then
                  c0(1,0:4) = (/  3.299899722e+0_dp, -1.217842603e+1_dp,  1.572745821e+1_dp, &
                     -8.111319659e+0_dp,  1.560141458e+0_dp /)
                  c0(2,0:4) = (/ -3.192174961e+0_dp,  1.088884936e+1_dp, -1.323621587e+1_dp, &
                     6.636540992e+0_dp, -8.283248890e-1_dp /)
                  c0(3,0:4) = (/  3.103362390e-1_dp, -1.100461299e+0_dp,  1.342127302e+0_dp, &
                     -1.038311585e+0_dp,  6.369080530e-1_dp /)
                  c0(4,0:4) = (/ -3.789397150e-1_dp,  1.452765098e+0_dp, -1.831668777e+0_dp, &
                     8.346864990e-1_dp, -9.323379900e-2_dp /)
                  c1(1,0:4) = (/ -4.416096315e-1_dp,  1.863738417e+0_dp, -2.695273977e+0_dp, &
                     1.589277096e+0_dp, -3.350371480e-1_dp /)
                  c1(2,0:4) = (/  9.641418065e-1_dp, -3.566095004e+0_dp,  4.694839548e+0_dp, &
                     -2.513976583e+0_dp,  4.330150925e-1_dp /)
                  c1(3,0:4) = (/ -3.761563300e-1_dp,  1.563266795e+0_dp, -2.387283764e+0_dp, &
                     1.616689565e+0_dp, -3.790348725e-1_dp /)
                  c1(4,0:4) = (/  2.284233995e-1_dp, -9.078127645e-1_dp,  1.346845869e+0_dp, &
                     -9.466905110e-1_dp,  3.041477165e-1_dp /)
                  c2(1,0:4) = (/  2.669608950e-2_dp, -1.173151950e-1_dp,  1.705683500e-1_dp, &
                     -1.004631905e-1_dp,  2.167064400e-2_dp /)
                  c2(2,0:4) = (/ -5.636994450e-2_dp,  2.146029235e-1_dp, -2.877285875e-1_dp, &
                     1.554898200e-1_dp, -2.821821450e-2_dp /)
                  c2(3,0:4) = (/  3.254421300e-2_dp, -1.383624650e-1_dp,  2.143779655e-1_dp, &
                     -1.445226355e-1_dp,  3.322128350e-2_dp /)
                  c2(4,0:4) = (/ -3.641687050e-2_dp,  1.368211995e-1_dp, -1.904624970e-1_dp, &
                     1.197266280e-1_dp, -3.044788050e-2_dp /)
               elseif (tp.le.10._dp) then
                  c0(1,0:4) = (/  2.088941304e+0_dp, -7.924782496e+0_dp,  1.030490483e+1_dp, &
                     -5.100518141e+0_dp,  9.187205666e-1_dp /)
                  c0(2,0:4) = (/ -1.329388501e+0_dp,  4.038387442e+0_dp, -4.044771639e+0_dp, &
                     1.436905595e+0_dp,  2.242395530e-1_dp /)
                  c0(3,0:4) = (/ -1.602468165e+0_dp,  6.568598385e+0_dp, -1.001324806e+1_dp, &
                     6.266949190e+0_dp, -1.016486384e+0_dp /)
                  c0(4,0:4) = (/  1.132825919e+0_dp, -4.368530854e+0_dp,  6.451051206e+0_dp, &
                     -4.361685669e+0_dp,  1.132980309e+0_dp /)
                  c1(1,0:4) = (/  5.744987000e-3_dp,  2.576259427e-1_dp, -6.515053707e-1_dp, &
                     4.844191472e-1_dp, -1.074755637e-1_dp /)
                  c1(2,0:4) = (/  3.296019472e-1_dp, -1.244044355e+0_dp,  1.636750238e+0_dp, &
                     -8.405937747e-1_dp,  1.051187500e-1_dp /)
                  c1(3,0:4) = (/  2.401828972e-1_dp, -8.878341890e-1_dp,  1.206059370e+0_dp, &
                     -6.762224942e-1_dp,  1.326717450e-1_dp /)
                  c1(4,0:4) = (/ -3.482113447e-1_dp,  1.253366932e+0_dp, -1.640951928e+0_dp, &
                     8.741404425e-1_dp, -1.145385190e-1_dp /)
                  c2(1,0:4) = (/ -1.422527975e-2_dp,  3.221345263e-2_dp, -1.943326838e-2_dp, &
                     4.642537500e-5_dp,  1.560960263e-3_dp /)
                  c2(2,0:4) = (/ -2.357369625e-3_dp,  1.788509087e-2_dp, -3.336493137e-2_dp, &
                     2.102700188e-2_dp, -2.806725250e-3_dp /)
                  c2(3,0:4) = (/ -1.704553588e-2_dp,  5.712493000e-2_dp, -6.908546350e-2_dp, &
                     3.470546388e-2_dp, -6.135529500e-3_dp /)
                  c2(4,0:4) = (/  1.769543038e-2_dp, -6.167275125e-2_dp,  7.742824750e-2_dp, &
                     -3.940152625e-2_dp,  5.271655750e-3_dp /)
               elseif (tp.le.20._dp) then
                  c0(1,0:4) = (/  2.617542346e+0_dp, -8.445966029e+0_dp,  9.293738910e+0_dp, &
                     -3.726657912e+0_dp,  5.317658086e-1_dp /)
                  c0(2,0:4) = (/ -4.308786366e+0_dp,  1.177854649e+1_dp, -1.022277525e+1_dp, &
                     2.920733202e+0_dp,  9.278494500e-2_dp /)
                  c0(3,0:4) = (/  6.421483327e+0_dp, -1.857655868e+1_dp,  1.758816186e+1_dp, &
                     -5.944994268e+0_dp,  8.238601990e-1_dp /)
                  c0(4,0:4) = (/ -3.702958887e+0_dp,  1.148749491e+1_dp, -1.199795553e+1_dp, &
                     4.334881025e+0_dp, -1.179988370e-1_dp /)
                  c1(1,0:4) = (/ -2.420942901e-1_dp,  8.055588697e-1_dp, -9.482688548e-1_dp, &
                     4.444700667e-1_dp, -6.861447887e-2_dp /)
                  c1(2,0:4) = (/  8.519828833e-1_dp, -2.605336006e+0_dp,  2.742478323e+0_dp, &
                     -1.129153506e+0_dp,  1.365721097e-1_dp /)
                  c1(3,0:4) = (/ -1.069654245e+0_dp,  3.228997716e+0_dp, -3.343633889e+0_dp, &
                     1.361134376e+0_dp, -1.786099381e-1_dp /)
                  c1(4,0:4) = (/  4.585202772e-1_dp, -1.395777107e+0_dp,  1.455667037e+0_dp, &
                     -5.944289378e-1_dp,  9.583863530e-2_dp /)
                  c2(1,0:4) = (/  5.272637540e-3_dp, -1.736800474e-2_dp,  2.035473924e-2_dp, &
                     -9.697268860e-3_dp,  1.544399362e-3_dp /)
                  c2(2,0:4) = (/ -2.480148458e-2_dp,  7.661266546e-2_dp, -8.215770368e-2_dp, &
                     3.504469894e-2_dp, -4.637515140e-3_dp /)
                  c2(3,0:4) = (/  3.369866338e-2_dp, -1.031066899e-1_dp,  1.098697633e-1_dp, &
                     -4.691078854e-2_dp,  6.589172980e-3_dp /)
                  c2(4,0:4) = (/ -1.461988376e-2_dp,  4.468139494e-2_dp, -4.774358166e-2_dp, &
                     2.048974484e-2_dp, -3.256268220e-3_dp /)
               elseif (tp.le.30._dp) then
                  c0(1,0:4) = (/  8.075657160e-1_dp, -2.514634025e+0_dp,  2.428229084e+0_dp, &
                     -5.280674060e-1_dp,  3.081807400e-2_dp /)
                  c0(2,0:4) = (/  2.306312502e+0_dp, -8.803354850e+0_dp,  1.205569561e+1_dp, &
                     -6.739957044e+0_dp,  1.367970879e+0_dp /)
                  c0(3,0:4) = (/ -9.543679249e+0_dp,  3.246642159e+1_dp, -4.022038919e+1_dp, &
                     2.118599073e+1_dp, -3.543353951e+0_dp /)
                  c0(4,0:4) = (/  7.377143959e+0_dp, -2.458427511e+1_dp,  2.985978468e+1_dp, &
                     -1.592854058e+1_dp,  3.400402971e+0_dp /)
                  c1(1,0:4) = (/ -5.870442220e-2_dp,  2.084715983e-1_dp, -2.596823647e-1_dp, &
                     1.219467346e-1_dp, -1.766298770e-2_dp /)
                  c1(2,0:4) = (/  4.911390390e-2_dp, -1.110857541e-1_dp,  4.243262600e-2_dp, &
                     4.372017540e-2_dp, -2.183352740e-2_dp /)
                  c1(3,0:4) = (/  6.774854291e-1_dp, -2.311028430e+0_dp,  2.882857193e+0_dp, &
                     -1.530830950e+0_dp,  2.814940050e-1_dp /)
                  c1(4,0:4) = (/ -6.722330151e-1_dp,  2.246992034e+0_dp, -2.735096147e+0_dp, &
                     1.415510895e+0_dp, -2.501827727e-1_dp /)
                  c2(1,0:4) = (/  6.280857200e-4_dp, -2.341971180e-3_dp,  3.089189300e-3_dp, &
                     -1.567578520e-3_dp,  2.491941400e-4_dp /)
                  c2(2,0:4) = (/ -1.195782780e-3_dp,  3.354906220e-3_dp, -2.851596000e-3_dp, &
                     5.527404800e-4_dp,  9.480188000e-5_dp /)
                  c2(3,0:4) = (/ -1.374541386e-2_dp,  4.628716672e-2_dp, -5.693341318e-2_dp, &
                     2.986001528e-2_dp, -5.497988800e-3_dp /)
                  c2(4,0:4) = (/  1.421752374e-2_dp, -4.727763704e-2_dp,  5.715022700e-2_dp, &
                     -2.934869280e-2_dp,  5.248797660e-3_dp /)
               endif
            elseif (mol_rat.eq.4._dp) then
               if (tp.le.2._dp) then
                  c0(1,0:4) = (/ -5.790495710e-1_dp,  1.544643160e+0_dp, -1.185283270e+0_dp, &
                     2.934555060e-1_dp,  3.213097320e-1_dp /)
                  c0(2,0:4) = (/ -1.101283986e+0_dp,  3.744270311e+0_dp, -4.888539825e+0_dp, &
                     2.749736151e+0_dp, -2.698127540e-1_dp /)
                  c0(3,0:4) = (/  3.590446760e+0_dp, -1.131184997e+1_dp,  1.247918026e+1_dp, &
                     -5.708099259e+0_dp,  9.850551230e-1_dp /)
                  c0(4,0:4) = (/ -2.501760427e+0_dp,  8.040241285e+0_dp, -8.858495322e+0_dp, &
                     3.645660877e+0_dp, -2.854465300e-1_dp /)
                  c1(1,0:4) = (/  1.899480863e+0_dp, -6.574744127e+0_dp,  7.991516215e+0_dp, &
                     -3.947963937e+0_dp,  4.942508950e-1_dp /)
                  c1(2,0:4) = (/  2.818093530e-1_dp, -9.547989520e-1_dp,  1.271623902e+0_dp, &
                     -6.715700190e-1_dp,  1.258012960e-1_dp /)
                  c1(3,0:4) = (/ -4.276503486e+0_dp,  1.394028924e+1_dp, -1.610193081e+1_dp, &
                     7.533934773e+0_dp, -9.286427640e-1_dp /)
                  c1(4,0:4) = (/  2.720486059e+0_dp, -8.898774607e+0_dp,  1.020784647e+1_dp, &
                     -4.719207655e+0_dp,  6.603212660e-1_dp /)
                  c2(1,0:4) = (/ -3.613943060e-1_dp,  1.292774778e+0_dp, -1.650468142e+0_dp, &
                     8.790923580e-1_dp, -1.246076140e-1_dp /)
                  c2(2,0:4) = (/ -1.068696860e-1_dp,  2.752761880e-1_dp, -2.084109760e-1_dp, &
                     1.913329000e-2_dp,  6.765944000e-3_dp /)
                  c2(3,0:4) = (/  1.212447928e+0_dp, -3.896163398e+0_dp,  4.371944412e+0_dp, &
                     -1.930912478e+0_dp,  2.068702000e-1_dp /)
                  c2(4,0:4) = (/ -7.849292540e-1_dp,  2.556968686e+0_dp, -2.894585214e+0_dp, &
                     1.280560230e+0_dp, -1.443986600e-1_dp /)
               elseif (tp.le.4._dp) then
                  c0(1,0:4) = (/ -3.396509140e-1_dp, -1.979575790e-1_dp,  2.177164989e+0_dp, &
                     -2.076972979e+0_dp,  7.479018990e-1_dp /)
                  c0(2,0:4) = (/ -7.231721260e-1_dp,  3.354365269e+0_dp, -5.611571547e+0_dp, &
                     3.775703659e+0_dp, -5.338845970e-1_dp /)
                  c0(3,0:4) = (/  1.001655337e+0_dp, -3.154400510e+0_dp,  3.660562599e+0_dp, &
                     -2.147871972e+0_dp,  7.655493260e-1_dp /)
                  c0(4,0:4) = (/ -9.633283050e-1_dp,  2.990216208e+0_dp, -3.198465199e+0_dp, &
                     1.301464732e+0_dp, -1.286865290e-1_dp /)
                  c1(1,0:4) = (/  1.497825229e+0_dp, -4.520110356e+0_dp,  4.536728880e+0_dp, &
                     -1.647388273e+0_dp,  1.033679525e-1_dp /)
                  c1(2,0:4) = (/ -3.278305630e-1_dp,  3.168573590e-1_dp,  8.563799710e-1_dp, &
                     -1.123019209e+0_dp,  2.966516465e-1_dp /)
                  c1(3,0:4) = (/ -7.853171255e-1_dp,  2.885346145e+0_dp, -4.056848560e+0_dp, &
                     2.555029124e+0_dp, -5.354389725e-1_dp /)
                  c1(4,0:4) = (/  5.982219920e-1_dp, -1.975339157e+0_dp,  2.466137438e+0_dp, &
                     -1.472251605e+0_dp,  3.993750385e-1_dp /)
                  c2(1,0:4) = (/ -2.204161535e-1_dp,  7.011080775e-1_dp, -7.636865395e-1_dp, &
                     3.214116475e-1_dp, -3.581418450e-2_dp /)
                  c2(2,0:4) = (/  1.034223070e-1_dp, -2.630757070e-1_dp,  1.799689200e-1_dp, &
                     -1.163399200e-2_dp, -1.264127050e-2_dp /)
                  c2(3,0:4) = (/  1.140526035e-1_dp, -4.080542145e-1_dp,  5.540577035e-1_dp, &
                     -3.315164755e-1_dp,  6.514475350e-2_dp /)
                  c2(4,0:4) = (/ -1.084052510e-1_dp,  3.577572305e-1_dp, -4.387382285e-1_dp, &
                     2.431312415e-1_dp, -5.311554650e-2_dp /)
               elseif (tp.le.6._dp) then
                  c0(1,0:4) = (/  2.013580586e+0_dp, -7.028498975e+0_dp,  9.175156901e+0_dp, &
                     -5.005988089e+0_dp,  1.068223823e+0_dp /)
                  c0(2,0:4) = (/ -4.149980560e-1_dp,  1.505960585e+0_dp, -2.571111715e+0_dp, &
                     2.088209175e+0_dp, -2.032058810e-1_dp /)
                  c0(3,0:4) = (/ -5.199022105e+0_dp,  1.703356977e+1_dp, -1.969413637e+1_dp, &
                     9.021212572e+0_dp, -9.318756740e-1_dp /)
                  c0(4,0:4) = (/  4.738259175e+0_dp, -1.539713500e+1_dp,  1.787382777e+1_dp, &
                     -8.733021270e+0_dp,  1.493387883e+0_dp /)
                  c1(1,0:4) = (/  1.908931465e-1_dp, -5.768494475e-1_dp,  3.458001825e-1_dp, &
                     1.614528900e-1_dp, -1.104556005e-1_dp /)
                  c1(2,0:4) = (/ -3.103732185e-1_dp,  6.723622460e-1_dp, -4.462026700e-2_dp, &
                     -5.247886920e-1_dp,  1.611036755e-1_dp /)
                  c1(3,0:4) = (/  2.212440393e+0_dp, -6.864673228e+0_dp,  7.235065311e+0_dp, &
                     -2.860802940e+0_dp,  2.895017735e-1_dp /)
                  c1(4,0:4) = (/ -2.113834710e+0_dp,  6.750228428e+0_dp, -7.516743150e+0_dp, &
                     3.279010757e+0_dp, -3.689039365e-1_dp /)
                  c2(1,0:4) = (/ -4.076010150e-2_dp,  1.422016875e-1_dp, -1.533288595e-1_dp, &
                     5.226480100e-2_dp, -2.378416500e-3_dp /)
                  c2(2,0:4) = (/  7.979709150e-2_dp, -2.364266360e-1_dp,  2.151902400e-1_dp, &
                     -5.572321600e-2_dp,  5.783025000e-4_dp /)
                  c2(3,0:4) = (/ -2.478444360e-1_dp,  7.677024865e-1_dp, -8.092520790e-1_dp, &
                     3.243737565e-1_dp, -3.500137050e-2_dp /)
                  c2(4,0:4) = (/  2.132597070e-1_dp, -6.744252155e-1_dp,  7.399636080e-1_dp, &
                     -3.175289740e-1_dp,  3.757454650e-2_dp /)
               elseif (tp.le.10._dp) then
                  c0(1,0:4) = (/  3.956798381e+0_dp, -1.332905159e+1_dp,  1.576014222e+1_dp, &
                     -7.325211727e+0_dp,  1.224942086e+0_dp /)
                  c0(2,0:4) = (/ -1.303793839e+0_dp,  3.917296544e+0_dp, -3.903609475e+0_dp, &
                     1.404652596e+0_dp,  2.181800080e-1_dp /)
                  c0(3,0:4) = (/  3.290108588e+0_dp, -8.922602992e+0_dp,  6.661333007e+0_dp, &
                     -8.255746460e-1_dp, -9.785729600e-2_dp /)
                  c0(4,0:4) = (/ -4.695639852e+0_dp,  1.396082495e+1_dp, -1.321672981e+1_dp, &
                     4.003366167e+0_dp,  3.558673200e-2_dp /)
                  c1(1,0:4) = (/ -4.351210112e-1_dp,  1.548329084e+0_dp, -1.978043925e+0_dp, &
                     1.036115868e+0_dp, -1.846925000e-1_dp /)
                  c1(2,0:4) = (/  2.890560450e-1_dp, -1.137403168e+0_dp,  1.562949116e+0_dp, &
                     -8.410707285e-1_dp,  1.079880892e-1_dp /)
                  c1(3,0:4) = (/ -9.676929977e-1_dp,  2.919302005e+0_dp, -2.875353251e+0_dp, &
                     1.062013040e+0_dp, -9.216776300e-2_dp /)
                  c1(4,0:4) = (/  1.125850147e+0_dp, -3.351933535e+0_dp,  3.260708406e+0_dp, &
                     -1.191072373e+0_dp,  1.562649345e-1_dp /)
                  c2(1,0:4) = (/  9.597319375e-3_dp, -3.697938412e-2_dp,  5.106223287e-2_dp, &
                     -2.908948325e-2_dp,  5.641115000e-3_dp /)
                  c2(2,0:4) = (/  4.580986000e-3_dp, -1.780621375e-3_dp, -1.572416387e-2_dp, &
                     1.597813950e-2_dp, -2.274263375e-3_dp /)
                  c2(3,0:4) = (/  4.636860987e-2_dp, -1.419552534e-1_dp,  1.437213099e-1_dp, &
                     -5.590703962e-2_dp,  5.443041750e-3_dp /)
                  c2(4,0:4) = (/ -6.463501837e-2_dp,  1.937695576e-1_dp, -1.926517185e-1_dp, &
                     7.369634112e-2_dp, -9.459122250e-3_dp /)
               elseif (tp.le.20._dp) then
                  c0(1,0:4) = (/  2.472031235e+0_dp, -7.807661507e+0_dp,  8.288656062e+0_dp, &
                     -3.103947432e+0_dp,  4.008428130e-1_dp /)
                  c0(2,0:4) = (/ -2.282918347e+0_dp,  5.550600684e+0_dp, -3.504052730e+0_dp, &
                     2.952828000e-3_dp,  4.704619360e-1_dp /)
                  c0(3,0:4) = (/  5.510889090e-1_dp, -4.804535410e-1_dp, -1.906547168e+0_dp, &
                     2.517517992e+0_dp, -3.431984370e-1_dp /)
                  c0(4,0:4) = (/  2.061251350e-1_dp, -6.946222580e-1_dp,  1.293120049e+0_dp, &
                     -1.510069473e+0_dp,  7.403141600e-1_dp /)
                  c1(1,0:4) = (/ -2.436411471e-1_dp,  7.967991438e-1_dp, -9.121627499e-1_dp, &
                     4.084902035e-1_dp, -5.791596570e-2_dp /)
                  c1(2,0:4) = (/  5.954380034e-1_dp, -1.823487097e+0_dp,  1.905988933e+0_dp, &
                     -7.690517933e-1_dp,  8.976107330e-2_dp /)
                  c1(3,0:4) = (/ -3.497821537e-1_dp,  1.023484596e+0_dp, -9.640800466e-1_dp, &
                     3.256050806e-1_dp, -3.546285540e-2_dp /)
                  c1(4,0:4) = (/  8.060278900e-3_dp, -6.235992400e-3_dp, -6.427743080e-2_dp, &
                     8.139834240e-2_dp, -4.452237600e-3_dp /)
                  c2(1,0:4) = (/  5.297004420e-3_dp, -1.704029088e-2_dp,  1.918897694e-2_dp, &
                     -8.539559700e-3_dp,  1.204454300e-3_dp /)
                  c2(2,0:4) = (/ -1.626596476e-2_dp,  5.049473008e-2_dp, -5.402371300e-2_dp, &
                     2.279324366e-2_dp, -2.974381060e-3_dp /)
                  c2(3,0:4) = (/  1.196772226e-2_dp, -3.679500702e-2_dp,  3.827279116e-2_dp, &
                     -1.569717004e-2_dp,  2.225962400e-3_dp /)
                  c2(4,0:4) = (/ -1.873681460e-3_dp,  5.754275400e-3_dp, -5.251633400e-3_dp, &
                     1.583625960e-3_dp, -4.346793200e-4_dp /)
               elseif (tp.le.30._dp) then
                  c0(1,0:4) = (/  7.441587330e-1_dp, -2.161982325e+0_dp,  1.716931614e+0_dp, &
                     4.642501200e-2_dp, -9.731476300e-2_dp /)
                  c0(2,0:4) = (/  2.731857021e+0_dp, -1.027177454e+1_dp,  1.384094089e+1_dp, &
                     -7.617635500e+0_dp,  1.534623680e+0_dp /)
                  c0(3,0:4) = (/  2.082318902e+1_dp, -5.984336250e+1_dp,  5.992378613e+1_dp, &
                     -2.504523761e+1_dp,  4.088875107e+0_dp /)
                  c0(4,0:4) = (/ -2.321515909e+1_dp,  6.854329471e+1_dp, -7.128236945e+1_dp, &
                     3.075433215e+1_dp, -4.321817070e+0_dp /)
                  c1(1,0:4) = (/ -6.052405840e-2_dp,  2.085410595e-1_dp, -2.441886051e-1_dp, &
                     9.919437570e-2_dp, -1.103283010e-2_dp /)
                  c1(2,0:4) = (/  3.979336700e-2_dp, -8.340100600e-2_dp,  1.687770400e-2_dp, &
                     5.052498390e-2_dp, -2.326505310e-2_dp /)
                  c1(3,0:4) = (/ -1.812676790e+0_dp,  5.278545817e+0_dp, -5.379549880e+0_dp, &
                     2.298152732e+0_dp, -3.541594518e-1_dp /)
                  c1(4,0:4) = (/  1.819215925e+0_dp, -5.348286015e+0_dp,  5.530199727e+0_dp, &
                     -2.407115612e+0_dp,  3.840676371e-1_dp /)
                  c2(1,0:4) = (/  4.608312400e-4_dp, -1.741584620e-3_dp,  2.219580820e-3_dp, &
                     -9.506994200e-4_dp,  1.056914600e-4_dp /)
                  c2(2,0:4) = (/ -1.020671360e-3_dp,  3.046363600e-3_dp, -2.930635600e-3_dp, &
                     8.658756200e-4_dp,  1.652090000e-5_dp /)
                  c2(3,0:4) = (/  3.443220378e-2_dp, -1.011407957e-1_dp,  1.044704496e-1_dp, &
                     -4.541766362e-2_dp,  7.080608360e-3_dp /)
                  c2(4,0:4) = (/ -3.387825318e-2_dp,  9.976198410e-2_dp, -1.035367675e-1_dp, &
                     4.534831962e-2_dp, -7.205344980e-3_dp /)
               endif
            else
               write(msg,'(a,f3.1,a)') 'Partial pressure ratio ',mol_rat,&
                  ' not available for correlation shan2018'
               call shutdown(msg)   
            endif
         
            !-----------------------------------------------------------
            !Computing the temperature coefficient of gray gas j
            !-----------------------------------------------------------
            tref = 2000._dp; sum_b = 0._dp
            do ii=0,4
               c_aux = c2(jj,ii)*tp**2._dp + c1(jj,ii)*tp + c0(jj,ii)
               sum_b = sum_b + c_aux*((ttmp/tref)**real(4-ii,dp))
            enddo
            a_func_res = sum_b
         
         case('stepwise')
            if (mol_rat.le.0.2_dp) then
               a_func_res = shan_wsgg_a(ttmp,tp,0.125_dp,jj,'none')
            elseif (mol_rat.le.0.4_dp) then
               a_func_res = shan_wsgg_a(ttmp,tp,0.25_dp,jj,'none')
            elseif (mol_rat.le.0.6_dp) then
               a_func_res = shan_wsgg_a(ttmp,tp,0.5_dp,jj,'none')
            elseif (mol_rat.le.0.9_dp) then
               a_func_res = shan_wsgg_a(ttmp,tp,0.75_dp,jj,'none')
            elseif (mol_rat.le.1.1_dp) then
               a_func_res = shan_wsgg_a(ttmp,tp,1._dp,jj,'none')
            elseif (mol_rat.le.2.5_dp) then
               a_func_res = shan_wsgg_a(ttmp,tp,2._dp,jj,'none')
            else
               a_func_res = shan_wsgg_a(ttmp,tp,4._dp,jj,'none')
            endif     
                
         case('linear')
            if (mol_rat.le.0.125_dp) then
               a_func_res = shan_wsgg_a(ttmp,tp,0.125_dp,jj,'none')
            elseif (mol_rat.le.0.25_dp) then
               x1 = 0.125_dp; y1 = shan_wsgg_a(ttmp,tp,x1,jj,'none')
               x2 = 0.250_dp; y2 = shan_wsgg_a(ttmp,tp,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.5_dp) then
               x1 = 0.250_dp; y1 = shan_wsgg_a(ttmp,tp,x1,jj,'none')
               x2 = 0.500_dp; y2 = shan_wsgg_a(ttmp,tp,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.0.75_dp) then
               x1 = 0.500_dp; y1 = shan_wsgg_a(ttmp,tp,x1,jj,'none')
               x2 = 0.750_dp; y2 = shan_wsgg_a(ttmp,tp,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.1._dp) then
               x1 = 0.75_dp; y1 = shan_wsgg_a(ttmp,tp,x1,jj,'none')
               x2 = 1.00_dp; y2 = shan_wsgg_a(ttmp,tp,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.2._dp) then
               x1 = 1._dp; y1 = shan_wsgg_a(ttmp,tp,x1,jj,'none')
               x2 = 2._dp; y2 = shan_wsgg_a(ttmp,tp,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            elseif (mol_rat.le.4._dp) then
               x1 = 2._dp; y1 = shan_wsgg_a(ttmp,tp,x1,jj,'none')
               x2 = 4._dp; y2 = shan_wsgg_a(ttmp,tp,x2,jj,'none')
               a_func_res = y1 + (y2 - y1)*(mol_rat - x1)/(x2 - x1 + small)
            else
               a_func_res = shan_wsgg_a(ttmp,tp,4._dp,jj,'none')
            endif   
      endselect
    
   endfunction shan_wsgg_a

!   !====================================================================
!   !Function to compute the temperature coefficient for a single
!   !gray gas for the RAD19 paper
!   !====================================================================
!   real(dp) function rad19_wsgg_a(ttmp,mol_rat,jj)
      
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      integer,intent(in) :: jj
!      integer :: counter,jm,js,ngm,ngs
!      real(dp),intent(in) :: mol_rat,ttmp
      
!      !-----------------------------------------------------------------
!      !Getting the indexes for the CO2-H2O mixture and for soot
!      !-----------------------------------------------------------------
!      ngm = wsgg_number_gray_gases('bordbar2014')
!      call get_wsgg_correlations('rad19')
!      ngs = number_wsgg_gases
!      counter = 1
!      main_loop: do jm=1,ngm+1
!         do js=1,ngs+1
!            if (jj.eq.counter) exit main_loop
!            counter = counter + 1
!         enddo
!      enddo main_loop
      
!      !-----------------------------------------------------------------
!      !Computing the absorption coefficient of the mixture
!      !-----------------------------------------------------------------
!      rad19_wsgg_a = bordbar_wsgg_a(ttmp,mol_rat,jm)*&
!                      fixed_wsgg_a(ttmp,1._dp,js,'rad19','soot')
!   endfunction rad19_wsgg_a

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!                                                                      !
!                          MISC FUNCTIONS                              !
!                                                                      !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
   
   !=======================================
   !Function to compute the gas emissivity 
   !=======================================
!   real(dp) function wsgg_gas_emissivity(ttmp,tp,xxc,xxw,xxs,pl,wsgg_spec)

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      character(*),intent(in) :: wsgg_spec
!      integer :: jj,ngg
!      real(dp),intent(in) :: pl,tp,ttmp,xxc,xxs,xxw
!      real(dp) :: a,kappa,sum_emi
   
!      !-----------------------------------------------------------------
!      !Some parameters
!      !-----------------------------------------------------------------
!      ngg = wsgg_number_gray_gases(wsgg_spec)
      
!      sum_emi = 0._dp
!      do jj=1,ngg
!         !--------------------------------------------------------------
!         !Computing local quantities
!         !--------------------------------------------------------------
!         kappa = wsgg_kappa_func(ttmp,tp,xxc,xxw,xxs,jj,wsgg_spec)
!         a = wsgg_a_func(ttmp,tp,xxc,xxw,xxs,jj,wsgg_spec)
         
!         !--------------------------------------------------------------
!         !Gray gas emissivity
!         !--------------------------------------------------------------
!         sum_emi = sum_emi + a*(1._dp-dexp(-kappa*pl))      
!      enddo
      
!      wsgg_gas_emissivity = sum_emi

!   endfunction wsgg_gas_emissivity

   !=====================================================================
   !Function to compute the gray gas absorption coefficient of soot only
   !=====================================================================
   real(dp) function gg_kappa_soot(ttmp,xxs,gg_spec)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use global_parameters, only: path_length,soot_constant,soot_density
      character(*),intent(in) :: gg_spec
      integer :: ii,jj,ngs,nps
      real(dp),intent(in) :: ttmp,xxs
      real(dp) :: a,c_s(3),dummy_sum,kappa,kappa_s
      
      !-----------------------------------------------------------------
      !Selecting the appropriate model
      !-----------------------------------------------------------------
      selectcase(trim(gg_spec))
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Fluent model
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('fluent')
            kappa_s = 1232.4_dp*soot_density*xxs*&
                      (1._dp+4.8e-4_dp*(ttmp-2000._dp))
            
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Kumar et al.(2002)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('kumar2002')
            kappa_s = 1864.32_dp*soot_constant*xxs*ttmp/7._dp
      
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Cassol et al.(2014)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('cassol2014')
            !Getting the soot WSGG correlations
            call get_wsgg_correlations('cassol2014','soot')
            ngs = number_wsgg_gases                                     !Number of gray gases
            nps = degree_wsgg_polynomial                                !Degree of the temperature coefficient polynomial
            
            dummy_sum = 0._dp
            do jj=1,ngs+1
               !Computing the absorption coefficient
               kappa = kappa_p(jj)*xxs*soot_constant
               
               !Computing the temperature coefficient
               a = 0._dp
               do ii=1,nps+1
                  a = a + b(jj,ii)*(ttmp)**(real(ii-1,dp))
               enddo
               
               !Computing the gas emissivity
               dummy_sum = dummy_sum + &
                           a*(1._dp - dexp(-kappa*path_length))
            enddo

            !Computing the absorption coefficient
            kappa_s =-log(1._dp-dummy_sum)/(path_length+small)
      
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Cassol et al.(2015)
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         case('cassol2015')
            !Constants of the model
            c_s(1:3) = (/ -281.19_dp, 3.7161_dp, -6.7737e-4_dp /)
            
            !Computing the polynomial
            dummy_sum = 0._dp
            do ii=1,3
               dummy_sum = dummy_sum + c_s(ii)*(ttmp**(real(ii-1,dp)))
            enddo
            kappa_s = dummy_sum*soot_constant*xxs*100._dp
      
      endselect
      gg_kappa_soot = kappa_s      
      
   endfunction gg_kappa_soot

   !====================================================================
   !Function to compute the Planck-mean absorption coefficient 
   !in any arbitrary WSGG model
   !====================================================================
   real(dp) function wsgg_kappa_planck(tmp,tp,xs)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      integer :: jj,ngg
      real(dp),intent(in) :: tmp,tp,xs(:)
      real(dp) :: aj,kj,sum_kp
      
      !-----------------------------------------------------------------
      !Computing
      !-----------------------------------------------------------------
      ngg = wsgg_number_gray_gases(which_wsgg_model)                    !Number of gray gases
      sum_kp = 0._dp
      do jj=1,ngg
         call wsgg_properties(kj,aj,jj,tmp,xs,tp,which_wsgg_model)
         sum_kp = sum_kp + kj*aj
      enddo      
      wsgg_kappa_planck = sum_kp
   
   endfunction wsgg_kappa_planck


endmodule wsgg_functions   
