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
   use slw_parameters
   use slw_functions
   use slw_routines
   use global_parameters
   use validation

   !--------------------------------------------------------------------
   !Declaration of variables
   !--------------------------------------------------------------------
   implicit none     
   character(200) :: local_file
   integer :: i,j,k,nx,ny,nz
   integer :: local_unit
   real(dp),allocatable,dimension(:,:,:) :: source_slw
   real(dp),allocatable,dimension(:,:,:,:) :: flux_slw

   call print_to_prompt('PROGRAM STARTED')
   
   !--------------------------------------------------------------------
   !File names
   !--------------------------------------------------------------------
   call print_to_prompt('Defining file names',3)
   local_file = 'slw1-results.csv'

   !--------------------------------------------------------------------
   !Set up case
   !--------------------------------------------------------------------
   call print_to_prompt('Set up case',3)
   validation_case = 'solovjovJQSRT2011'
   validation_subcase = 2
   call validation_setup                                                !Não mexer

   !--------------------------------------------------------------------
   !Set up FVM method
   !--------------------------------------------------------------------
   call print_to_prompt('Set up FVM method parameters',3)
   call set_default_fvm_parameters                                      !Não mexer
   rte_solution_method = 'FVM'                                          !Não mexer
   fvm_angles = 48                                                      !Não mexer
   call initalize_fvm_parameters                                        !Não mexer
   call build_fvm_angular_grid                                          !Não mexer
   
   !--------------------------------------------------------------------
   !SLW parameters
   !--------------------------------------------------------------------
   call print_to_prompt('Defining SLW parameters',3)
   call set_default_slw_parameters                                      !Set default parameters for the model
   albdf_file(id_h2o) = 'ALBDFs/albdfH2O.data'                          !Não mexer
   !albdf_file(id_co2) = 'ALBDFs/albdfCO2.data'                          !Não mexer
   albdf_info_file(id_h2o) = 'ALBDFs/infoH2O.data'                      !Não mexer
   !albdf_info_file(id_co2) = 'ALBDFs/infoCO2.data'                      !Não mexer
   slw_nonuniform_method = 'rank_correlated'                            !Não mexer
   slw_mixture_method = 'multiplication'                                !Não mexer
   slw1_approach = 'kp-epsilon'
   slw1_length(1) = 1._dp
   slw1_length(2) = 2._dp
   slw_Fmin=0._dp; slw_Fmax=1._dp 
!   slw_print_properties = .false.

   !--------------------------------------------------------------------
   !Allocate arrays
   !--------------------------------------------------------------------
   call print_to_prompt('Allocating arrays',3)
   nx = xpoints; ny = ypoints; nz = zpoints                             !Não mexer
   allocate(flux_slw(3,nx,ny,nz))                                       !Não mexer
   allocate(source_slw(nx,ny,nz))                                       !Não mexer
   
   !--------------------------------------------------------------------
   !Solve the radiation field
   !--------------------------------------------------------------------
   call print_to_prompt('Solving the radiation field: SLW model',3)
   call slw1_black_solution(flux_slw,source_slw)                        !Não mexer
   
   !--------------------------------------------------------------------
   !Dump local data
   !--------------------------------------------------------------------
   call print_to_prompt('Dumping local data',3)
   local_unit = get_file_unit()
   open(unit=local_unit,file=trim(local_file))
   write(local_unit,'(1000(a,:,","))') '#x','T','xh2o','S','qx'          !,'xco2'
   write(local_unit,'(1000(a,:,","))') '#[m]','[K]','[]',&               !,'[]'
      '[kW/m3]','[kW/m2]'
   j = 1; k = 1
   do i=1,nx
      write(local_unit,'(1000(e26.15e3,:,","))') &
         x(i),T(i,j,k),xs(id_h2o,i,j,k),&                            !,xs(id_co2,i,j,k)
         source_slw(i,j,k)*1.e-3_dp,flux_slw(1,i,j,k)*1.e-3_dp
   enddo
   close(local_unit)
   
   call print_to_prompt('PROGRAM ENDED')
   
end program main_slw
