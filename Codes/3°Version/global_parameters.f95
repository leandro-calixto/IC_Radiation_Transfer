module global_parameters 
 
   use precision_parameters, only: dp
 
   implicit none 
   character(200) :: soot_correlation_file
   character(20) :: rte_solution_method
   real(dp) :: path_length
   real(dp) ::soot_constant,soot_density
   real(dp) :: total_bundles
   
   integer :: number_of_species,number_of_surface_bands
   integer :: id_ch4=-1,id_co=-1,id_co2=-1,id_h2o=-1,id_soot=-1
   character(10) :: spec_name(10)
   
   !File extensions
   character(4),parameter :: dat_ext = '.dat'
   character(5),parameter :: bin_ext = '.data'
   
   !TRI parameters
   logical :: tri_correction
   
   !Parameters for the get_random_gauss function
   real(dp) :: next_gauss_rnd_num
   logical :: gauss_rnd_num_ready
   
   logical :: debug_mode = .false.

end module global_parameters 
