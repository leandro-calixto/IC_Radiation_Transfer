module wbm_parameters

   use precision_parameters, only: dp
   implicit none
   integer,parameter :: wbm_max_nsp = 5                                !Maximum size for arrays related to the WBM input data
   character(80) :: wbm_prefix
   character(80) :: wbm_print_file
   character(200) :: wbm_data_file(wbm_max_nsp)
   integer :: number_wbm_bands
   integer :: wbm_data_nsp
   integer :: wbm_data_ntg(wbm_max_nsp)
   integer :: wbm_data_nxs(wbm_max_nsp)
   logical :: wbm_print_properties
   real(dp),allocatable,dimension(:) :: wbm_lbound,wbm_ubound
   real(dp),allocatable,dimension(:,:) :: wbm_data_tg
   real(dp),allocatable,dimension(:,:) :: wbm_data_xs
   real(dp),allocatable,dimension(:,:,:,:) :: wbm_kappa

contains

   !==================================================
   !Subroutine with the default parameters of the WBM
   !==================================================
   !(this does not need to be included in the main 
   !code if the user knows what they are doing)
   subroutine set_default_wbm_parameters

      wbm_print_file = 'null'
      wbm_print_properties = .false.
      
   endsubroutine set_default_wbm_parameters
   
endmodule wbm_parameters
