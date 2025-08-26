!#######################################################################
!Module concerning with validation procedures
!#######################################################################
module validation

   !====================================================================
   !Modules & misc
   !====================================================================
   use precision_parameters, only: dp
   implicit none
   character(200) :: validation_case
   integer :: validation_subcase

contains

   !====================================================================
   !Subroutine to setup validation cases
   !====================================================================
   subroutine validation_setup
   
      use comp_functions, only: shutdown
      use constants, only: pi,sigrpi,twopi
      use global_parameters
      use mesh
      implicit none
      integer :: i,ii,j,jj,k
      real(dp) :: alpha,radius,y0
      real(dp) :: asum,bsum,xch4,xco,xco2,xh2o,xsoot
      real(dp) :: bb(1:10,0:10),cc(1:10)
      
      mesh_practice = 'B'
      one_d = .false.; two_d = .false.; three_d = .false.
      
      selectcase(trim(validation_case))

            case('los')
            !Participating species
            number_of_species = 1
            id_co2 = 1
         
            !Computational domain and mesh
            one_d = .true.
            xmin = 0._dp; xmax = 1._dp
            xcells = 50
            call initialize_mesh_variables
            call build_spatial_mesh
            
            !Composition profiles
            p = 1._dp; T = 0._dp; xs = 0._dp
            j = 1; k = 1
            do i=1,xpoints
               selectcase(validation_subcase)
                  case(1)
                     alpha = dsin(pi*(x(i) - xmin)/(xmax - xmin))
                     T(i,j,k) = 400._dp + 1400._dp*alpha**2._dp
                     xs(id_co2,i,j,k) = 0.2_dp*alpha**2._dp
               endselect
            enddo
            
         case('dorigon2013')
            !Participating species
            number_of_species = 2
            id_h2o = 1; id_co2 = 2
         
            !Computational domain and mesh
            one_d = .true.
            xmin = 0._dp; xmax = 1._dp
            xcells = 50
            call initialize_mesh_variables
            call build_spatial_mesh
            
            !Composition profiles
            p = 1._dp; T = 0._dp; xs = 0._dp
            j = 1; k = 1
            do i=1,xpoints
               selectcase(validation_subcase)
                  case(1)
                     alpha = dsin(pi*(x(i) - xmin)/(xmax - xmin))
                     T(i,j,k) = 400._dp + 1400._dp*alpha**2._dp
                     xs(id_co2,i,j,k) = 0.2_dp*alpha**2._dp
                     xs(id_h2o,i,j,k) = 0.4_dp*alpha**2._dp
                  case(2)
                     alpha = dsin(twopi*(x(i) - xmin)/(xmax - xmin))
                     T(i,j,k) = 400._dp + 1400._dp*alpha**2._dp
                     xs(id_co2,i,j,k) = 0.2_dp*alpha**2._dp
                     xs(id_h2o,i,j,k) = 0.4_dp*alpha**2._dp
                  case(3)
                     if (((x(i) - xmin)/(xmax - xmin)).le.0.25_dp) then
                        alpha = dsin(twopi*(x(i) - xmin)/(xmax - xmin))
                        T(i,j,k) = 880._dp + 920._dp*alpha**2._dp
                        xs(id_co2,i,j,k) = 0.25_dp*alpha**2._dp
                        xs(id_h2o,i,j,k) = 0.50_dp*alpha**2._dp
                     else
                        alpha = pi*((x(i) - xmin)/(xmax - xmin) - 0.25_dp)
                        T(i,j,k) = 400._dp + 1400._dp*&
                             (1._dp - (dsin(2._dp*alpha/3._dp))**1.5_dp)
                        xs(id_co2,i,j,k) = 0.25_dp*&
                           (1._dp - dsin(2._dp*alpha/3._dp))
                        xs(id_h2o,i,j,k) = 0.50_dp*&
                           (1._dp - dsin(2._dp*alpha/3._dp))
                     endif
                  case default
                     call shutdown('validation_setup: &
                           &validation_subcase incorrectly specified')
               endselect
            enddo
         
         case('goutiere2000')
            !Participating species
            number_of_species = 1!2
            !id_h2o = 1;
            id_co2 = 1!2
         
            !Computational domain and mesh
            two_d = .true.
            xmin = 0._dp; xmax = 1._dp
            ymin = 0._dp; ymax = 0.5_dp
            xcells = 61; ycells = 31
            call initialize_mesh_variables
            call build_spatial_mesh
            
            !Composition profiles
            p = 1._dp; T = 0._dp; xs = 0._dp
            k = 1
            do i=2,xpoints-1
               do j=2,ypoints-1
                  selectcase(validation_subcase)
                     case(1)
                        T(i,j,k) = 1000._dp
                        xs(id_co2,i,j,k) = 0.1_dp
                     case(2)
                        alpha = (1._dp - 2._dp*dabs(x(i) - 0.5_dp))*&
                           (1._dp - 4._dp*dabs(y(j) - 0.25_dp))
                        T(i,j,k) = 1200._dp*(alpha/3._dp + 1._dp)
                        xs(id_co2,i,j,k) = 0.02_dp*(alpha/3._dp + 1._dp)
                     case(3)
                        T(i,j,k) = 1000._dp
                        xs(id_h2o,i,j,k) = 0.2_dp
                     case(4)
                        alpha = (1._dp - 2._dp*dabs(x(i) - 0.5_dp))*&
                           (1._dp - 4._dp*dabs(y(j) - 0.25_dp))
                        T(i,j,k) = 1200._dp*(alpha/3._dp + 1._dp)
                        xs(id_h2o,i,j,k) = 0.04_dp*(alpha/3._dp + 1._dp)
                     case(5)
                        y0 = dabs(0.25_dp - y(j))/0.25_dp
                        alpha = 1._dp - 3._dp*y0**2._dp + &
                                2._dp*y0**3._dp
!                       xs(id_h2o,:,:,:) = 0.2_dp
                        xs(id_co2,:,:,:) = 0.1_dp
                        if (x(i).le.0.1_dp) &
                           T(i,j,k) = (14000._dp*x(i) - 400._dp)*&
                                    alpha + 800._dp
                        if (x(i).gt.0.1_dp) &
                           T(i,j,k) = (10000._dp/9._dp)*(1._dp - x(i))*&
                                    alpha + 800._dp

                     case default
                        call shutdown('validation_setup: &
                           &validation_subcase incorrectly specified')
                  endselect
               enddo
            enddo

         case('deshmukh2009')         
            !Participating species
            number_of_species = 2
            id_h2o = 1                                                  !Surrogate for kappa
            id_co2 = 2
            
            !Computational domain and mesh
            two_d = .true.
            xmin = -1._dp; xmax = 1._dp
            ymin = -1._dp; ymax = 1._dp
            xcells = 41; ycells = 41
            call initialize_mesh_variables
            call build_spatial_mesh
            
            !Composition profiles
            p = 1._dp; T = 0._dp; xs = 0._dp
            k = 1
            do i=1,xpoints
               do j=1,ypoints
                  alpha = 0.5_dp*(x(i)**2._dp + y(i)**2._dp)
                  T(i,j,k) = (1._dp + 20._dp*alpha*(1._dp - alpha))
                  T(i,j,k) = (T(i,j,k)/sigrpi)**(1._dp/4._dp)
                  xs(id_h2o,i,j,k) = (1._dp + 15._dp*(1._dp - alpha)**2._dp)
               enddo
            enddo

         case('modest2012')         
            !Participating species
            number_of_species = 2
            id_h2o = 1
            id_co2 = 2
            
            !Computational domain and mesh
            two_d = .true.
            xmin = -1._dp; xmax = 1._dp
            ymin = -1._dp; ymax = 1._dp
            xcells = 41; ycells = 41
            call initialize_mesh_variables
            call build_spatial_mesh
            
            !Composition profiles
            p = 1._dp; T = 0._dp; xs = 1._dp
            k = 1
            do i=2,xpoints-1
               do j=2,ypoints-1
                  alpha = x(i)**2._dp + y(i)**2._dp
                  T(i,j,k) = (1._dp + 5._dp*alpha*(2._dp - alpha))
                  T(i,j,k) = (T(i,j,k)/sigrpi)**(1._dp/4._dp)
                  xs(id_h2o,i,j,k) = (1._dp + 5._dp*(2._dp - alpha)**2._dp)
               enddo
            enddo
      
         case('step-2d')
            !Participating species
            number_of_species = 1; id_co2 = 1

            !Computational domain and mesh
            two_d = .true.
            xmin = 0._dp; xmax = 1._dp
            ymin = 0._dp; ymax = 1._dp
            xcells = 41; ycells = 41
            call initialize_mesh_variables
            call build_spatial_mesh
            
            !Composition profiles
            p = 1._dp; T = 0._dp; xs = 0._dp
            k = 1
            do i=2,xpoints-1
               do j=2,ypoints-1
                  !xs(id_h2o,i,j,k) = 0.2_dp
                  xs(id_co2,i,j,k) = 0.1_dp
                  if (x(i).lt.0.2_dp) T(i,j,k) = 2000._dp
                  if (x(i).ge.0.2_dp) T(i,j,k) = 500._dp
               enddo
            enddo
            
         case('step-1d')
            !Participating species
            number_of_species = 1; id_co2 = 1

            !Computational domain and mesh
            one_d = .true.
            xmin = 0._dp; xmax = 1._dp
            xcells = 200
            call initialize_mesh_variables
            call build_spatial_mesh
            
            !Composition profiles
            p = 1._dp; T = 0._dp; xs = 0._dp
            j = 1; k = 1
            do i=2,xpoints-1
               xs(id_co2,i,j,k) = 0.1_dp
               if (x(i).lt.0.2_dp) T(i,j,k) = 2000._dp
               if (x(i).ge.0.2_dp) T(i,j,k) = 500._dp
            enddo
      
         case('2021-wsgg-case1')
            !Computational domain and mesh
            one_d = .true.
            xmin = 0._dp; xmax = 1._dp
            xcells = 200            
            
            selectcase(validation_subcase)
               
               case(2)
                  !Species
                  number_of_species = 2
                  id_h2o = 1; id_co2 = 2
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  xs(id_h2o,:,:,:) = .1_dp
                  xs(id_co2,:,:,:) = .05_dp
                  
                  xs(id_h2o,1,:,:) = .0_dp; xs(id_h2o,xpoints,:,:) = .0_dp
                  xs(id_co2,1,:,:) = .0_dp; xs(id_co2,xpoints,:,:) = .0_dp
               
               case(3)
                  !Species
                  number_of_species = 3
                  id_h2o = 1; id_co2 = 2; id_soot = 3
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  xs(id_h2o,:,:,:) = .1_dp
                  xs(id_co2,:,:,:) = .05_dp
                  xs(id_soot,:,:,:) = 1.e-7_dp
                  
                  xs(id_h2o,1,:,:) = .0_dp; xs(id_h2o,xpoints,:,:) = .0_dp
                  xs(id_co2,1,:,:) = .0_dp; xs(id_co2,xpoints,:,:) = .0_dp
                  xs(id_soot,1,:,:) = .0_dp; xs(id_soot,xpoints,:,:) = .0_dp
               
               case(4)
                  !Species
                  number_of_species = 4
                  id_h2o = 1; id_co2 = 2; id_soot = 3; id_co = 4
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  xs(id_h2o,:,:,:) = .1_dp
                  xs(id_co2,:,:,:) = .05_dp
                  xs(id_soot,:,:,:) = 1.e-7_dp
                  xs(id_co,:,:,:) = .04_dp
                  
                  xs(id_h2o,1,:,:) = .0_dp; xs(id_h2o,xpoints,:,:) = .0_dp
                  xs(id_co2,1,:,:) = .0_dp; xs(id_co2,xpoints,:,:) = .0_dp
                  xs(id_soot,1,:,:) = .0_dp; xs(id_soot,xpoints,:,:) = .0_dp
                  xs(id_co,1,:,:) = .0_dp; xs(id_co,xpoints,:,:) = .0_dp
                  
               case(5)
                  !Species
                  number_of_species = 5
                  id_h2o = 1; id_co2 = 2; id_soot = 3; id_co = 4; id_ch4 = 5
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  xs(id_h2o,:,:,:) = .1_dp
                  xs(id_co2,:,:,:) = .05_dp
                  xs(id_soot,:,:,:) = 1.e-7_dp
                  xs(id_co,:,:,:) = .04_dp
                  xs(id_ch4,:,:,:) = .04_dp
                  
                  xs(id_h2o,1,:,:) = .0_dp; xs(id_h2o,xpoints,:,:) = .0_dp
                  xs(id_co2,1,:,:) = .0_dp; xs(id_co2,xpoints,:,:) = .0_dp
                  xs(id_soot,1,:,:) = .0_dp; xs(id_soot,xpoints,:,:) = .0_dp
                  xs(id_co,1,:,:) = .0_dp; xs(id_co,xpoints,:,:) = .0_dp
                  xs(id_ch4,1,:,:) = .0_dp; xs(id_ch4,xpoints,:,:) = .0_dp
                  
            endselect
            
            !Pressure and temperature profiles
            p = 1._dp; T = 1800._dp
            T(1,:,:) = 300._dp
            T(xpoints,:,:) = 300._dp

         case('2021-wsgg-case2')
            !Computational domain and mesh
            one_d = .true.
            xmin = 0._dp; xmax = 1._dp
            xcells = 200            
            
            selectcase(validation_subcase)
               case(2)
                  !Species
                  number_of_species = 2
                  id_h2o = 1; id_co2 = 2
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  do i=1,xpoints
                     xs(id_h2o,i,:,:) = &
                        .1_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co2,i,:,:) = &
                        .05_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                  enddo
                  
               case(3)
                  !Species
                  number_of_species = 3
                  id_h2o = 1; id_co2 = 2; id_soot = 3
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  do i=1,xpoints
                     xs(id_h2o,i,:,:) = &
                        .1_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co2,i,:,:) = &
                        .05_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_soot,i,:,:) = &
                        1.e-7_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                  enddo
               
               case(4)
                  !Species
                  number_of_species = 4
                  id_h2o = 1; id_co2 = 2; id_soot = 3; id_co = 4
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  do i=1,xpoints
                     xs(id_h2o,i,:,:) = &
                        .1_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co2,i,:,:) = &
                        .05_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_soot,i,:,:) = &
                        1.e-7_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co,i,:,:) = &
                        .04_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                  enddo
                  
               case(5)
                  !Species
                  number_of_species = 5
                  id_h2o = 1; id_co2 = 2; id_soot = 3; id_co = 4; id_ch4 = 5
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  do i=1,xpoints
                     xs(id_h2o,i,:,:) = &
                        .1_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co2,i,:,:) = &
                        .05_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_soot,i,:,:) = &
                        1.e-7_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co,i,:,:) = &
                        .04_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_ch4,i,:,:) = &
                        .04_dp*(dsin(1._dp*pi*x(i)/(xmax - xmin)))**2
                  enddo
                  
            endselect
            
            !Pressure and temperature profiles
            p = 1._dp
            do i=1,xpoints
               T(i,:,:) = 300._dp + &
                  1500._dp*(dsin(pi*x(i)/(xmax - xmin)))**2
            enddo
            T(1,:,:) = 300._dp
            T(xpoints,:,:) = 300._dp

         case('2021-wsgg-case3')
            !Computational domain and mesh
            two_d = .true.
            xmin = -1.2_dp; xmax = 1.2_dp
            ymin = 0.25_dp; ymax = 5._dp
            xcells = 48; ycells = 96            
            
            selectcase(validation_subcase)
               case(2)
                  !Species
                  number_of_species = 2
                  id_h2o = 1; id_co2 = 2
                  
               case(3)
                  !Species
                  number_of_species = 3
                  id_h2o = 1; id_co2 = 2; id_soot = 3
                  
               case(4)
                  !Species
                  number_of_species = 4
                  id_h2o = 1; id_co2 = 2; id_soot = 3; id_co = 4
                  
               case(5)
                  !Species
                  number_of_species = 5
                  id_h2o = 1; id_co2 = 2; id_soot = 3; id_co = 4; id_ch4 = 5
                  
            endselect
            
            !Mesh arrays
            call initialize_mesh_variables
            call build_spatial_mesh
            
            do i=1,xpoints
               do j=1,ypoints
                  radius = min(abs(x(i)),1.2_dp)

                  !Temperature profile
                  cc(1:4) = (/ 4664.74122541374_dp, &
                               0.210072431323855_dp,&
                               6.67190142493502_dp, &
                               3.58716736913414_dp /)
                  bb(1,0:4) = (/ -3165.10473950106_dp,&
                                  67392.6543905640_dp,&
                                 -301363.056289393_dp,&
                                  383102.403315692_dp,&
                                 -147288.257468992_dp /)
                  bb(2,0:4) = (/ 1203.24171735589_dp,&
                                -9798.76868179257_dp,&
                                 31348.4664493848_dp,& 
                                -33779.9331570168_dp,&
                                 11464.5415568891_dp /)
                  bb(3,0:4) = (/ 6149.68226144286_dp,&
                                -140568.767911782_dp,&
                                 613743.909559609_dp,&
                                -766110.222649197_dp,&
                                 290489.935933031_dp /)
                  bb(4,0:4) = (/ -2544.45353359079_dp,&
                                  79750.6876530367_dp,&
                                 -335957.473611559_dp,&
                                  408178.833603675_dp,&
                                 -151668.920493191_dp /)
                  bsum = 0._dp
                  do ii=1,4
                     asum = 0._dp
                     do jj=0,4
                        asum = asum + bb(ii,jj)*(radius**jj)
                     enddo
                     bsum = bsum + asum*(1._dp - dexp(-cc(ii)*y(j)))
                  enddo
                  T(i,j,:) = max(bsum,300._dp)
                  
                  !CO2 mole fraction profile
                  cc(1:4) = (/ 3.01642949876554_dp, &
                               6.60007062717780_dp, &
                               0.240590806935068_dp,&
                               68.8225811897622_dp /)
                  bb(1,0:4) = (/ -0.135581740056667_dp,&
                                  4.26445903362457_dp, &
                                 -18.2804049746972_dp, &
                                  22.5688816339567_dp, &
                                 -8.52118831824259_dp  /)
                  bb(2,0:4) = (/  0.376390250806248_dp,&
                                 -7.56414169858682_dp, &
                                  34.4449052764868_dp, &
                                 -44.2003158516455_dp, &
                                  17.1252768277067_dp  /)
                  bb(3,0:4) = (/  0.100566720864732_dp, &
                                 -0.865473263633969_dp, &
                                  2.83939793322812_dp,  &
                                 -3.13723531175224_dp,  &
                                  1.09515160897571_dp /)
                  bb(4,0:4) = (/ -0.237657521292915_dp, &
                                  3.90942933824436_dp,  &
                                 -18.3723313134367_dp,  &
                                  24.0495672462375_dp,  &
                                 -9.44033431732454_dp  /)
                  bsum = 0._dp
                  do ii=1,4
                     asum = 0._dp
                     do jj=0,4
                        asum = asum + bb(ii,jj)*(radius**jj)
                     enddo
                     bsum = bsum + asum*(1._dp - dexp(-cc(ii)*y(j)))
                  enddo
                  xco2 = max(bsum,0._dp)
                  
                  !H2O mole fraction profile
                  cc(1:4) = (/ 4.83369950888859_dp,&
                               89.7406474065145_dp, &
                               0.130414300395918_dp,&
                               4.81003154380251_dp /)
                  bb(1,0:4) = (/ -45.4183903077146_dp,&
                                  207.520501667280_dp,&
                                  1308.78631024042_dp,&
                                 -2533.18943227628_dp,&
                                  1100.83993541154_dp /)
                  bb(2,0:4) = (/  0.104358171997472_dp,&
                                 -2.35208352362674_dp, &
                                 -2.73241407938704_dp, &
                                  9.26286100548538_dp, &
                                 -4.48315723382295_dp /)
                  bb(3,0:4) = (/  0.0614102996237949_dp,&
                                 -0.473062633853121_dp, &
                                  2.14671723789815_dp,  &
                                 -2.59947995625785_dp,  &
                                  0.918234940884834_dp /)
                  bb(4,0:4) = (/  45.3966481340503_dp,&
                                 -204.863438902823_dp,&
                                 -1307.71315501106_dp,&
                                  2525.89052803827_dp,&
                                 -1097.05291610239_dp /)
                  bsum = 0._dp
                  do ii=1,4
                     asum = 0._dp
                     do jj=0,4
                        asum = asum + bb(ii,jj)*(radius**jj)
                     enddo
                     bsum = bsum + asum*(1._dp - dexp(-cc(ii)*y(j)))
                  enddo
                  xh2o = max(bsum,0._dp)
                  
                  !Soot volume fraction profile
                  cc(1:4) = (/  10.0886242515777_dp, &
                                6.24714960817725_dp, &
                                0.237969648946843_dp,&
                                10.9098689880697_dp /)
                  bb(1,0:4) = (/  0.142494337125892_dp,&
                                 -15.6524948185568_dp, &
                                  145.649117980623_dp, &
                                 -229.514002622469_dp, &
                                  99.4362746008028_dp /)
                  bb(2,0:4) = (/  0.0776390626162987_dp,&
                                  1.73009250508413_dp,  &
                                 -16.3278233789723_dp,  &
                                  25.4950111965039_dp,  &
                                 -10.9826389113942_dp  /)
                  bb(3,0:4) = (/ -0.00782016291190991_dp,&
                                 -0.0404185091766290_dp, &
                                  0.282886357757614_dp,  &
                                 -0.395644171037235_dp,  &
                                  0.161036374034997_dp /)
                  bb(4,0:4) = (/ -0.205405259329755_dp,&
                                  13.9450566000409_dp, &
                                 -129.553203762301_dp, &
                                  204.346418475683_dp, &
                                 -88.5861831755553_dp /)
                  bsum = 0._dp
                  do ii=1,4
                     asum = 0._dp
                     do jj=0,4
                        asum = asum + bb(ii,jj)*(radius**jj)
                     enddo
                     bsum = bsum + asum*(1._dp - dexp(-cc(ii)*y(j)))
                  enddo
                  xsoot = max(bsum*1.0613e-4_dp,0._dp)

                  !CH4 mole fraction profile
                  cc(1:4) = (/  1.86058084508000_dp,&
                                112.260137050129_dp,&
                                1.84369700946972_dp,&
                                9.56090722635469_dp /)
                  bb(1,0:4) = (/  135.063514701867_dp,&
                                 -776.330001668283_dp,&
                                  905.206563454137_dp,&
                                  121.860239924745_dp,&
                                 -357.630286928144_dp /)
                  bb(2,0:4) = (/  11.8979702941785_dp,&
                                 -46.7911966283841_dp,&
                                  129.532630536039_dp,&
                                 -100.561909218538_dp,&
                                  13.5172956800249_dp  /)
                  bb(3,0:4) = (/ -134.234541007791_dp,&
                                  766.621955329639_dp,&
                                 -887.087615303635_dp,&
                                 -129.832382267588_dp,&
                                  356.766111260115_dp /)
                  bb(4,0:4) = (/ -12.6677246261884_dp,&
                                  56.3701986764812_dp,&
                                 -147.581137161693_dp,&
                                  108.521762103360_dp,&
                                 -12.6445438771843_dp /)
                  bsum = 0._dp
                  do ii=1,4
                     asum = 0._dp
                     do jj=0,4
                        asum = asum + bb(ii,jj)*(radius**jj)
                     enddo
                     bsum = bsum + asum*(1._dp - dexp(-cc(ii)*y(j)))
                  enddo
                  xch4 = min(max(bsum*0.01_dp,0._dp),1._dp)
            
                  !CO mole fraction profile
                  cc(1:4) = (/  2.15945078998510_dp,&
                                88.7106994277095_dp,&
                                11.7153268528063_dp,&
                                2.08759304304498_dp /)
                  bb(1,0:4) = (/  1.25933956001295_dp,&
                                  10.3157225154377_dp,&
                                 -75.7165741465910_dp,&
                                  110.110481773630_dp,&
                                 -45.9137423105262_dp /)
                  bb(2,0:4) = (/ -0.0942385586002653_dp,&
                                  2.12179674764034_dp,  &
                                 -19.8126963383394_dp,  &
                                  31.5789812447823_dp,  &
                                 -13.7845706446015_dp /)
                  bb(3,0:4) = (/  0.0666961363663398_dp,&
                                 -2.52202234198902_dp,  &
                                  23.2876983615641_dp,  &
                                 -36.9045303560429_dp,  &
                                  16.0619976371957_dp /)
                  bb(4,0:4) = (/ -1.21530855918155_dp,&
                                 -9.92817568763522_dp,&
                                  72.1692529922646_dp,&
                                 -104.666968510003_dp,&
                                  43.5873481186482_dp /)
                  bsum = 0._dp
                  do ii=1,4
                     asum = 0._dp
                     do jj=0,4
                        asum = asum + bb(ii,jj)*(radius**jj)
                     enddo
                     bsum = bsum + asum*(1._dp - dexp(-cc(ii)*y(j)))
                  enddo
                  xco = max(bsum,0._dp)
                  
                  !Boundary values
                  if ((i.eq.1).or.(i.eq.xpoints).or.&
                      (j.eq.1).or.(j.eq.ypoints)) then
                     T(i,j,:) = 0._dp
                     xco2 = 0._dp
                     xh2o = 0._dp
                     xsoot = 0._dp
                     xch4 = 0._dp
                     xco = 0._dp
                  endif
                  
                  selectcase(validation_subcase)
                     case(2)
                        xs(id_h2o,i,j,:) = xh2o
                        xs(id_co2,i,j,:) = xco2
                     case(3)
                        xs(id_h2o,i,j,:) = xh2o
                        xs(id_co2,i,j,:) = xco2
                        xs(id_soot,i,j,:) = xsoot
                     case(4)
                        xs(id_h2o,i,j,:) = xh2o
                        xs(id_co2,i,j,:) = xco2
                        xs(id_soot,i,j,:) = xsoot
                        xs(id_co,i,j,:) = xco
                     case(5)
                        xs(id_h2o,i,j,:) = xh2o
                        xs(id_co2,i,j,:) = xco2
                        xs(id_soot,i,j,:) = xsoot
                        xs(id_co,i,j,:) = xco
                        xs(id_ch4,i,j,:) = xch4
                  endselect
                  
               enddo
            enddo
            T(:,ypoints,:) = 473._dp
            
            !Pressure
            p = 1._dp

         case('solovjov2000')                                           !Sec 4.3, 10.1016/S0022-4073(99)00133-8
            !Computational domain and mesh
            one_d = .true.
            xmin = 0._dp; xmax = 1._dp
            xcells = 200            
            
            selectcase(validation_subcase)
               case(2)
                  !Species
                  number_of_species = 2
                  id_h2o = 1; id_co2 = 2
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  do i=1,xpoints
                     xs(id_h2o,i,:,:) = &
                        .1_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co2,i,:,:) = &
                        .05_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                  enddo
                  
               case(3)
                  !Species
                  number_of_species = 3
                  id_h2o = 1; id_co2 = 2; id_soot = 3
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  do i=1,xpoints
                     xs(id_h2o,i,:,:) = &
                        .1_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co2,i,:,:) = &
                        .05_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_soot,i,:,:) = &
                        1.e-7_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                  enddo
               
               case(4)
                  !Species
                  number_of_species = 4
                  id_h2o = 1; id_co2 = 2; id_soot = 3; id_co = 4
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  do i=1,xpoints
                     xs(id_h2o,i,:,:) = &
                        .1_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co2,i,:,:) = &
                        .05_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_soot,i,:,:) = &
                        1.e-7_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co,i,:,:) = &
                        .04_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                  enddo
                  
               case(5)
                  !Species
                  number_of_species = 5
                  id_h2o = 1; id_co2 = 2; id_soot = 3; id_co = 4; id_ch4 = 5
                  
                  !Mesh arrays
                  call initialize_mesh_variables
                  call build_spatial_mesh
                  
                  !Species concentrations
                  do i=1,xpoints
                     xs(id_h2o,i,:,:) = &
                        .1_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co2,i,:,:) = &
                        .05_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_soot,i,:,:) = &
                        1.e-7_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_co,i,:,:) = &
                        .04_dp*(dsin(2._dp*pi*x(i)/(xmax - xmin)))**2
                     xs(id_ch4,i,:,:) = &
                        .04_dp*(dsin(1._dp*pi*x(i)/(xmax - xmin)))**2
                  enddo
                  
            endselect
            
            !Pressure and temperature profiles
            p = 1._dp
            do i=1,xpoints
               T(i,:,:) = 300._dp + &
                  1500._dp*(dsin(pi*x(i)/(xmax - xmin)))**2
            enddo
            T(1,:,:) = 300._dp
            T(xpoints,:,:) = 300._dp

         case('solovjovJHT2011')
            !Participating species
            number_of_species = 2
            id_h2o = 1; id_co2 = 2
            
            !Set domain dimensions based on case
            select case(validation_subcase)
               case(1)  
                  !Computational domain and mesh
                  one_d = .true.
                  xcells = 50
                  xmin = 0._dp; xmax = 1._dp
                  number_of_surface_bands = 1
                     
                  !Initialize mesh after setting dimensions
                  call initialize_mesh_variables
                  call build_spatial_mesh
                     
                  !Emissivities
                  xmin_emissivity = 1._dp
                  xmax_emissivity = 1._dp
                     
                  !Composition profiles
                  p = 1._dp; T = 10._dp; xs = 0._dp        
                  j = 1; k = 1
                  do i=2,xpoints-1
                     T(i,j,k) = 1000._dp
                     xs(id_h2o,i,j,k) = 0.2_dp
                     xs(id_co2,i,j,k) = 0.1_dp
                  enddo
               
               case(3)
                  !Computational domain and mesh
                  one_d = .true.
                  xcells = 300
                  xmin = 0._dp; xmax = 3._dp
                  number_of_surface_bands = 1
                     
                  !Initialize mesh after setting dimensions
                  call initialize_mesh_variables
                  call build_spatial_mesh
                     
                  !Emissivities
                  xmin_emissivity = 0.8_dp
                  xmax_emissivity = 0.8_dp
                     
                  !Composition profiles
                  p = 1._dp
                  j = 1; k = 1
                  do i=1,xpoints
                     alpha = dcos(pi*x(i)/(xmax - xmin))
                     T(i,j,k) = 700._dp - 300._dp*alpha                   
                     xs(id_h2o,i,j,k) = 0.2_dp - 0.15_dp*alpha
                     xs(id_co2,i,j,k) = 2._dp*xs(id_h2o,i,j,k)/3._dp
                  enddo
                  T(1,j,k) = 300._dp
                  T(xpoints,j,k) = 300._dp
               
               case(4) 
                  !Computational domain and mesh
                  one_d = .true.
                  xcells = 500
                  xmin = 0._dp; xmax = 10._dp
                  number_of_surface_bands = 1
                     
                  !Initialize mesh after setting dimensions
                  call initialize_mesh_variables
                  call build_spatial_mesh
                     
                  !Emissivities
                  xmin_emissivity = 1._dp
                  xmax_emissivity = 1._dp
                     
                  !Composition profiles
                  p = 1._dp; T = 500._dp; xs = 0._dp
                  j = 1; k = 1
                  do i=2,xpoints-1
                     alpha = dcos(pi*x(i)/(xmax - xmin))
                     T(i,j,k) = 1000._dp + 500._dp*alpha                   
                     xs(id_h2o,i,j,k) = 0.5_dp - 0.5_dp*alpha
                     xs(id_co2,i,j,k) = 0.5_dp + 0.5_dp*alpha
                  enddo
               case default
                  call shutdown('validation_setup: &
                     &validation_subcase incorrectly specified')
            end select        

         case('solovjovJQSRT2011')
            !Participating species
            number_of_species = 1  
            id_h2o = 1
            
            selectcase(validation_subcase)
               case(2)
                  !Computational domain setup
                  one_d = .true.
                  xcells = 100
                  xmin = 0._dp; xmax = 2._dp  
                  number_of_surface_bands = 1
                  
                  !Initialize mesh after setting dimensions
                  call initialize_mesh_variables
                  call build_spatial_mesh
                     
                  !Emissivities
                  xmin_emissivity = 1._dp
                  xmax_emissivity = 1._dp
                  
                  !Composition profiles
                  p = 1._dp
                  j = 1; k = 1
                  do i=1,xpoints
                     alpha = dcos(pi*x(i)/(xmax - xmin))
                     T(i,j,k) = 1000._dp + 250._dp*alpha
                     xs(id_h2o,i,j,k) = 0.1_dp
                  enddo
               
               case(3)
                  !Computational domain setup
                  one_d = .true.
                  xcells = 500
                  xmin = 0._dp; xmax = 10._dp  
                  number_of_surface_bands = 1
                  
                  !Initialize mesh after setting dimensions
                  call initialize_mesh_variables
                  call build_spatial_mesh
                     
                  !Emissivities
                  xmin_emissivity = 1._dp
                  xmax_emissivity = 1._dp
                  
                  !Composition profiles
                  p = 1._dp
                  j = 1; k = 1
                  do i=1,xpoints
                     alpha = dsin(pi*x(i)/(xmax - xmin))
                     T(i,j,k) = 500._dp + 1000._dp*alpha
                     xs(id_h2o,i,j,k) = alpha
                  enddo
               
               case default
                  call shutdown('validation_setup: &
                     &validation_subcase incorrectly specified')
            endselect
                     
      endselect

   
   endsubroutine validation_setup

endmodule validation

