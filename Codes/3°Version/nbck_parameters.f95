!#######################################################################
!Module compiling LBL parameters
!#######################################################################
module nbck_parameters

   use precision_parameters, only: dp
   implicit none
   integer,parameter :: nbck_max_nsp = 5                                !Maximum size for arrays related to the NBCK input data
   character(200) :: nbck_prefix                                        !Prefix of all NBCK external data files
   character(200) :: nbck_info_file                                     !Name of the file to which all information about the NBCK
                                                                        !  data files will be outputted (for verification purposes)
   character(200) :: nbck_data_file(nbck_max_nsp)                       !Array with the name of the external NBCK data files 
                                                                        !  (one for each species)
   character(200) :: nbck_print_file
   character(200) :: nbck_premix_file                                   !File storing the premixed NBCK
   character(20) :: nbck_mixing_method                                  !String to specify the method for mixing the NBCK distributions
                                                                        !  of single species
   character(20) :: nbck_quadrature                                     !String to specify which quadrature is to be used 
                                                                        !  for the NBCK solution
   integer :: nbck_data_ntg(nbck_max_nsp)                               !Array with the number of temperature points in the NBCK
                                                                        !  data files (it can be different for each species)
   integer :: nbck_data_nxs(nbck_max_nsp)                               !Array with the number of mole fraction points in the NBCK
                                                                        !  data files (it can be different for each species)
   integer :: nbck_data_nkg(nbck_max_nsp)                               !Array with the number of g points in the NBCK
                                                                        !  data files (it can be different for each species)
   integer :: nbck_data_nsp                                             !Number of species for which the NBCK data file is available
   integer :: number_nbck_bands                                         !Number of narrow bands
   integer :: max_nbck_points
   logical :: nbck_print_properties                                     !If .true., dump local properties used in the RTE solution
                                                                        !  to an external file
   logical :: nbck_dump_info                                            !If .true., write information on the NBCK data files
                                                                        !  to an external file
   logical :: nbck_loaded                                               !Flag to indicate if the NBCK has already been loaded
   logical :: nbck_is_premixed                                          !Flag to indicate if the NBCK has already been premixed
   logical :: scaled_nbck                                               !If .true., use the scaled NBCK model
   real(dp) :: nbck_T0
   real(dp),allocatable,dimension(:,:) :: nbck_g                        !Array with the g values in the NBCK data files
   real(dp),allocatable,dimension(:,:) :: nbck_data_tg                  !Array with the temperature values in the NBCK data files
   real(dp),allocatable,dimension(:,:) :: nbck_data_xs                  !Array with the mole fraction values in the NBCK data files
   real(dp),allocatable,dimension(:) :: nbck_xs0
   real(dp),allocatable,dimension(:) :: nbck_lbound,nbck_ubound         !Arrays with the lower and upper bounds of the narrow bands
   real(dp),allocatable,dimension(:,:,:,:,:) :: nbck_k                  !Array with the k values in the NBCK data files

   !Variables for the premixed NBCK database
   integer :: nbck_mix_nkg                                              !Number of k and g points in the premixed NBCK database
   integer :: nbck_mix_ntg                                              !Number of temperature values in the premixed NBCK database
   integer :: nbck_mix_nxc                                              !Number of CO2 mole fraction values in the premixed NBCK database
   integer :: nbck_mix_nxw                                              !Number of H2O mole fraction valuesin the premixed NBCK database
   integer :: nbck_mix_nfv                                              !Number of soot volume fraction values in the premixed NBCK database
   logical :: nbck_mix_from_file                                        !If .true., read the premixed NBCK from an external file
   real(dp),allocatable,dimension(:) :: nbck_mix_g                      !Array with the g values in the premixed NBCK database
   real(dp),allocatable,dimension(:) :: nbck_mix_tg                     !Array with the temperature values in the premixed NBCK database
   real(dp),allocatable,dimension(:) :: nbck_mix_xc                     !Array with the CO2 mole fraction values in the premixed NBCK database
   real(dp),allocatable,dimension(:) :: nbck_mix_xw                     !Array with the H2O mole fraction values in the premixed NBCK database
   real(dp),allocatable,dimension(:) :: nbck_mix_fv                     !Array with the soot volume fraction values in the premixed NBCK database
   real(dp),allocatable,dimension(:,:,:,:,:,:) :: nbck_mix_k            !Array with the k values in the premixed NBCK database   
   
contains

   !====================================================================
   !Subroutine to set the default NBCK parameters
   !====================================================================
   subroutine set_default_nbck_parameters

      max_nbck_points = 10000
      nbck_prefix = ''
      nbck_info_file = 'nbck_info'
      nbck_quadrature = 'GC-even'
      nbck_mixing_method = 'riazzi'
      nbck_print_file = 'null'
      nbck_premix_file = 'null'
      scaled_nbck = .false.
      nbck_loaded = .false.
      nbck_dump_info = .false.
      nbck_mix_from_file = .false.
   
   endsubroutine set_default_nbck_parameters


endmodule nbck_parameters
