program main_slw

   !====================================================================
   !INITIALIZATION
   !====================================================================
   !--------------------------------------------------------------------
   !Modules
   !--------------------------------------------------------------------
   use precision_parameters
   use fvm_parameters
   use mesh
   use comp_functions, only: get_file_unit,print_to_prompt
   use lbl_parameters
   use lbl_routines
   use slw_parameters
   use slw_functions
   use slw_routines
   use fvm_parameters
   use global_parameters
   use validation

   !--------------------------------------------------------------------
   !Declaration of variables
   !--------------------------------------------------------------------
   implicit none     
   character(200) :: local_file
   integer :: i,isp,j,k,ntmp,nx,ny,nz
   integer :: local_unit
   real(dp),allocatable,dimension(:) :: T_albdf,xs_albdf
   real(dp),allocatable,dimension(:,:,:) :: source_lbl,source_slw,&
                                            source_slw1
   real(dp),allocatable,dimension(:,:,:,:) :: flux_lbl,flux_slw,&
                                              flux_slw1
   logical :: convert_lbl,generate_albdf,solve_lbl,solve_slw,solve_slw1

   call print_to_prompt('PROGRAM STARTED')
   
   !--------------------------------------------------------------------
   !Input parameters
   !--------------------------------------------------------------------
   call print_to_prompt('Set up the input parameters',3)
   
   !Flags
   convert_lbl = .false.
   generate_albdf = .false.
   solve_lbl = .false.
   solve_slw = .true.
   solve_slw1 = .true.

   !--------------------------------------------------------------------
   !Convert formatted LBL data to unformatted data
   !--------------------------------------------------------------------
   if (convert_lbl) then
      call print_to_prompt('Convert formatted LBL data to &
                           &unformatted data',3)

      number_of_species = 2; id_h2o = 1; id_co2 = 2
      call set_default_lbl_parameters
      lbl_data_cm = .true.; lbl_ceta_input = .true.
      lbl_xceta_input = .false.
      lbl_data_prefix = 'Databases/'
      lbl_data_ext = '.dat'
      lbl_data_nxs(id_h2o) = 9; lbl_data_nxs(id_co2) = 3
      lbl_data_xs(1,id_h2o) = 0._dp
      lbl_data_xs(2,id_h2o) = 0.05_dp
      lbl_data_xs(3,id_h2o) = 0.1_dp
      lbl_data_xs(4,id_h2o) = 0.2_dp
      lbl_data_xs(5,id_h2o) = 0.4_dp
      lbl_data_xs(6,id_h2o) = 0.5_dp
      lbl_data_xs(7,id_h2o) = 0.6_dp
      lbl_data_xs(8,id_h2o) = 0.8_dp
      lbl_data_xs(9,id_h2o) = 1._dp
      lbl_data_xs(1,id_co2) = 0._dp
      lbl_data_xs(2,id_co2) = 0.5_dp
      lbl_data_xs(3,id_co2) = 1._dp
      lbl_data_file(1,id_h2o) = 'h2o-fraga-p1-y0'
      lbl_data_file(2,id_h2o) = 'h2o-fraga-p1-y005'
      lbl_data_file(3,id_h2o) = 'h2o-fraga-p1-y01'
      lbl_data_file(4,id_h2o) = 'h2o-fraga-p1-y02'
      lbl_data_file(5,id_h2o) = 'h2o-fraga-p1-y04'
      lbl_data_file(6,id_h2o) = 'h2o-fraga-p1-y05'
      lbl_data_file(7,id_h2o) = 'h2o-fraga-p1-y06'
      lbl_data_file(8,id_h2o) = 'h2o-fraga-p1-y08' 
      lbl_data_file(9,id_h2o) = 'h2o-fraga-p1-y1'
      lbl_data_file(1,id_co2) = 'co2-fraga-p1-y0'
      lbl_data_file(2,id_co2) = 'co2-fraga-p1-y05'
      lbl_data_file(3,id_co2) = 'co2-fraga-p1-y1'
      lbl_data_ntg = 23
      do i=1,23
         lbl_data_tg(i,:,:) = 300._dp + 100._dp*real(i-1,dp)
      enddo
      call convert_lbl_data

   endif

   !--------------------------------------------------------------------
   !Set up case
   !--------------------------------------------------------------------
   call print_to_prompt('Set up case',3)
   validation_case = 'Pattern Search'
   validation_subcase = 3
   call validation_setup                                                !Não mexer

   !--------------------------------------------------------------------
   !SLW parameters
   !--------------------------------------------------------------------
   call print_to_prompt('Defining SLW parameters',3)
   call set_default_slw_parameters                                      !Set default parameters for the model
   if (id_h2o.gt.0) then
      albdf_file(id_h2o) = 'ALBDFs/albdfH2O.data'
      albdf_info_file(id_h2o) = 'ALBDFs/infoH2O.data'
   endif
   if (id_co2.gt.0) then
      albdf_file(id_co2) = 'ALBDFs/albdfCO2.data'
      albdf_info_file(id_co2) = 'ALBDFs/infoCO2.data'
   endif
   slw_nonuniform_method = 'rank_correlated'                            !Não mexer
   slw_mixture_method = 'one_species'                                !'one_species', 'multiplication'
   slw1_approach = 'Q-Q'
   slw1_length(1) = 10.0_dp
   slw1_length(2) = 10.0_dp
   slw1_position(1) = 0.5_dp
   slw1_position(2) = 1.5_dp
   slw_Fmin=0._dp; slw_Fmax=1._dp 

   !--------------------------------------------------------------------
   !Generate the ALBDF
   !--------------------------------------------------------------------
   if (generate_albdf) then
      call print_to_prompt('Generate the ALBDF',3)
      
      !LBL parameters
      number_of_species = 2; id_h2o = 1; id_co2 = 2
      call set_default_lbl_parameters
      lbl_binary_data = .true.
      lbl_data_prefix = 'Databases/'; lbl_data_ext = '.data'
      lbl_data_nxs(id_h2o) = 9; lbl_data_nxs(id_co2) = 3
      lbl_data_xs(1,id_h2o) = 0._dp
      lbl_data_xs(2,id_h2o) = 0.05_dp
      lbl_data_xs(3,id_h2o) = 0.1_dp
      lbl_data_xs(4,id_h2o) = 0.2_dp
      lbl_data_xs(5,id_h2o) = 0.4_dp
      lbl_data_xs(6,id_h2o) = 0.5_dp
      lbl_data_xs(7,id_h2o) = 0.6_dp
      lbl_data_xs(8,id_h2o) = 0.8_dp
      lbl_data_xs(9,id_h2o) = 1._dp
      lbl_data_xs(1,id_co2) = 0._dp
      lbl_data_xs(2,id_co2) = 0.5_dp
      lbl_data_xs(3,id_co2) = 1._dp
      lbl_data_file(1,id_h2o) = 'h2o-fraga-p1-y0'
      lbl_data_file(2,id_h2o) = 'h2o-fraga-p1-y005'
      lbl_data_file(3,id_h2o) = 'h2o-fraga-p1-y01'
      lbl_data_file(4,id_h2o) = 'h2o-fraga-p1-y02'
      lbl_data_file(5,id_h2o) = 'h2o-fraga-p1-y04'
      lbl_data_file(6,id_h2o) = 'h2o-fraga-p1-y05'
      lbl_data_file(7,id_h2o) = 'h2o-fraga-p1-y06'
      lbl_data_file(8,id_h2o) = 'h2o-fraga-p1-y08' 
      lbl_data_file(9,id_h2o) = 'h2o-fraga-p1-y1'
      lbl_data_file(1,id_co2) = 'co2-fraga-p1-y0'
      lbl_data_file(2,id_co2) = 'co2-fraga-p1-y05'
      lbl_data_file(3,id_co2) = 'co2-fraga-p1-y1'
      lbl_data_ntg = 23
      do i=1,23
         lbl_data_tg(i,:,:) = 300._dp + 100._dp*real(i-1,dp)
      enddo
            
      !Define discrete values for which to generate the ALBDF
      ntmp = 23; allocate(T_albdf(ntmp))
      do i=1,ntmp
         T_albdf(i) = 300._dp + 100._dp*real(i-1,dp)
      enddo
      allocate(xs_albdf(9))
      xs_albdf(1:9) = (/ 0._dp, 0.05_dp, 0.1_dp, 0.2_dp, 0.4_dp, &
                         0.5_dp, 0.6_dp, 0.8_dp, 1._dp /)
      debug_mode = .false.; print_xeta = .false.
      
      !Generate the ALBDF
      do isp=1,number_of_species
         call generate_fsck(T_albdf,xs_albdf,T_albdf,isp,&
            albdf_file(isp),albdf_info_file(isp),albdf=.true.,&
            pressure_based=.true.,tabulate_k=.true.)
      enddo

   endif

   !--------------------------------------------------------------------
   !Set up FVM method
   !--------------------------------------------------------------------
   call print_to_prompt('Set up FVM method parameters',3)
   call set_default_fvm_parameters                                      !Não mexer
   rte_solution_method = 'FVM'                                          !Não mexer
   fvm_angles = 100                                                     !Não mexer
   call initalize_fvm_parameters                                        !Não mexer
   call build_fvm_angular_grid                                          !Não mexer
   
   !--------------------------------------------------------------------
   !Allocate arrays
   !--------------------------------------------------------------------
   call print_to_prompt('Allocating arrays',3)
   nx = xpoints; ny = ypoints; nz = zpoints                             !Não mexer
   allocate(flux_lbl(3,nx,ny,nz))                                       !Não mexer
   allocate(flux_slw(3,nx,ny,nz))                                       !Não mexer
   allocate(flux_slw1(3,nx,ny,nz))                                      !Não mexer
   allocate(source_lbl(nx,ny,nz))                                       !Não mexer
   allocate(source_slw(nx,ny,nz))                                       !Não mexer
   allocate(source_slw1(nx,ny,nz))                                      !Não mexer
   
   !--------------------------------------------------------------------
   !Solve the radiation field
   !--------------------------------------------------------------------
   if (solve_slw) then
      call print_to_prompt('Solving the radiation field: SLW model',3)
      call slw_gray_solution(flux_slw,source_slw,input_ngas=20)
   endif
   
   if (solve_slw1) then
      call print_to_prompt('Solving the radiation field: SLW-1 model',3)
      call slw1_gray_solution(flux_slw1,source_slw1)
   endif
   
   if (solve_lbl) then
      call print_to_prompt('Define LBL parameters',3)
      call set_default_lbl_parameters
      lbl_binary_data = .true.
      lbl_data_prefix = 'Databases/'; lbl_data_ext = '.data'
      if (id_h2o.gt.0) then
         lbl_data_nxs(id_h2o) = 9
         lbl_data_xs(1,id_h2o) = 0._dp
         lbl_data_xs(2,id_h2o) = 0.05_dp
         lbl_data_xs(3,id_h2o) = 0.1_dp
         lbl_data_xs(4,id_h2o) = 0.2_dp
         lbl_data_xs(5,id_h2o) = 0.4_dp
         lbl_data_xs(6,id_h2o) = 0.5_dp
         lbl_data_xs(7,id_h2o) = 0.6_dp
         lbl_data_xs(8,id_h2o) = 0.8_dp
         lbl_data_xs(9,id_h2o) = 1._dp
         lbl_data_file(1,id_h2o) = 'h2o-fraga-p1-y0'
         lbl_data_file(2,id_h2o) = 'h2o-fraga-p1-y005'
         lbl_data_file(3,id_h2o) = 'h2o-fraga-p1-y01'
         lbl_data_file(4,id_h2o) = 'h2o-fraga-p1-y02'
         lbl_data_file(5,id_h2o) = 'h2o-fraga-p1-y04'
         lbl_data_file(6,id_h2o) = 'h2o-fraga-p1-y05'
         lbl_data_file(7,id_h2o) = 'h2o-fraga-p1-y06'
         lbl_data_file(8,id_h2o) = 'h2o-fraga-p1-y08' 
         lbl_data_file(9,id_h2o) = 'h2o-fraga-p1-y1'
      endif
      if (id_co2.gt.0) then
         lbl_data_nxs(id_co2) = 3
         lbl_data_xs(1,id_co2) = 0._dp
         lbl_data_xs(2,id_co2) = 0.5_dp
         lbl_data_xs(3,id_co2) = 1._dp
         lbl_data_file(1,id_co2) = 'co2-fraga-p1-y0'
         lbl_data_file(2,id_co2) = 'co2-fraga-p1-y05'
         lbl_data_file(3,id_co2) = 'co2-fraga-p1-y1'
      endif
      
      lbl_data_ntg = 23
      do i=1,23
         lbl_data_tg(i,:,:) = 300._dp + 100._dp*real(i-1,dp)
      enddo
         
      call print_to_prompt('Solving the radiation field: LBL model',3)
      call lbl_nongray_solution(flux_lbl,source_lbl)
   endif
   
   
   !   call slw1_black_solution(flux_slw,source_slw)                        !Não mexer
   
   !--------------------------------------------------------------------
   !Dump local data
   !--------------------------------------------------------------------
   !Output files
   write(local_file, '(A,A,A,I0,A,A,A)') 'Results/Simulation_', &
    trim(validation_case), '_Subcase', validation_subcase, '_', &
    trim(slw1_approach), '.csv'
   
   call print_to_prompt('Dumping local data',3)
   local_unit = get_file_unit()
   open(unit=local_unit,file=trim(local_file))
   j = 1; k = 1
   if (id_co2.le.0) then
      write(local_unit,'(1000(a,:,","))') '#x','T','xh2o',&
         'S_LBL','S_SLW','S_SLW1','q_LBL','q_SLW','q_SLW1'
      write(local_unit,'(1000(a,:,","))') '#[m]','[K]','[]',&
         '[kW/m3]','[kW/m2]'
      do i=1,nx
         write(local_unit,'(1000(e26.15e3,:,","))') &
            x(i),T(i,j,k),xs(id_h2o,i,j,k),&
            source_lbl(i,j,k)*1.e-3_dp,source_slw(i,j,k)*1.e-3_dp,&
            source_slw1(i,j,k)*1.e-3_dp,&
            flux_lbl(1,i,j,k)*1.e-3_dp,flux_slw(1,i,j,k)*1.e-3_dp,&
            flux_slw1(1,i,j,k)*1.e-3_dp
      enddo

   elseif (id_h2o.le.0) then
      write(local_unit,'(1000(a,:,","))') '#x','T','xco2',&
         'S_LBL','S_SLW','S_SLW1','q_LBL','q_SLW','q_SLW1'
      write(local_unit,'(1000(a,:,","))') '#[m]','[K]','[]',&
         '[kW/m3]','[kW/m2]'
      do i=1,nx
         write(local_unit,'(1000(e26.15e3,:,","))') &
            x(i),T(i,j,k),xs(id_co2,i,j,k),&
            source_lbl(i,j,k)*1.e-3_dp,source_slw(i,j,k)*1.e-3_dp,&
            source_slw1(i,j,k)*1.e-3_dp,&
            flux_lbl(1,i,j,k)*1.e-3_dp,flux_slw(1,i,j,k)*1.e-3_dp,&
            flux_slw1(1,i,j,k)*1.e-3_dp
      enddo

   else
      write(local_unit,'(1000(a,:,","))') '#x','T','xh2o','xco2',&
         'S_LBL','S_SLW','S_SLW1','q_LBL','q_SLW','q_SLW1'
      write(local_unit,'(1000(a,:,","))') '#[m]','[K]','[]','[]',&
         '[kW/m3]','[kW/m2]'
      do i=1,nx
         write(local_unit,'(1000(e26.15e3,:,","))') &
            x(i),T(i,j,k),xs(id_h2o,i,j,k),xs(id_co2,i,j,k),&
            source_lbl(i,j,k)*1.e-3_dp,source_slw(i,j,k)*1.e-3_dp,&
            source_slw1(i,j,k)*1.e-3_dp,&
            flux_lbl(1,i,j,k)*1.e-3_dp,flux_slw(1,i,j,k)*1.e-3_dp,&
            flux_slw1(1,i,j,k)*1.e-3_dp
      enddo
   endif
   close(local_unit)
   
   call print_to_prompt('PROGRAM ENDED')
   
end program main_slw

