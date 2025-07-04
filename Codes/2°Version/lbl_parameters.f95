!#######################################################################
!Module compiling LBL parameters
!#######################################################################
module lbl_parameters

   !====================================================================
   !Declaration of variables
   !====================================================================
   use precision_parameters, only: dp
   implicit none
   integer,parameter :: lbl_max_nsp=5,lbl_max_nxs=10,lbl_max_ntg=30     !These parameters give the maximum sizes for arrays related to the
                                                                        !  input LBL spectral data. (ns -> number of species, nx -> number
                                                                        !  of mole fraction values, nt -> number of temperature values)
   character(80) :: which_lbl_soot                                      !String that specifies the soot treatment   
   character(80) :: lbl_data_file(lbl_max_nxs,lbl_max_nsp)              !Array with the name of the external files with the LBL spectral data
   character(80) :: lbl_data_prefix
   character(80) :: lbl_data_ext
   character(80) :: lbl_deta_file                                       !External file containing the size of each LBL line
   character(80) :: lbl_data_averaging                                  !String that specifies how to average the LBL spectral data to obtain
                                                                        !  wavenumber positions and spectral absorption coefficients at a given line
                                                                        !  ('arithmetic', 'geometric', 'upwind', 'none')
   integer :: lbl_nlines                                                !Total number of lines in the LBL data files (must be equal for all files)
   integer :: lbl_nsp                                                   !Actual number of participating species
   integer :: lbl_nxs                                                   !Actual (maximum) number of discrete mole fraction values
   integer :: lbl_ntg                                                   !Actual (maximum) number of discrete temperature values
   integer :: lbl_data_nxs(lbl_max_nsp)                                 !Array with the number of discrete mole fraction values for which the LBL data
                                                                        !  are available for each species
   integer :: lbl_data_ntg(lbl_max_nxs,lbl_max_nsp)                     !Array with the number of dicsrete temperature values for which the LBL data
                                                                        !  are available for each species and each mole fraction value
   integer,allocatable,dimension(:) :: lbl_data_unit                    !Array with the units for reading of each LBL data file (the code mounts and
                                                                        !  handles this array alone; it is not to be specified by the user)
   real(dp) :: lbl_data_xs(lbl_max_nxs,lbl_max_nsp)                     !Array with the values of discrete mole fraction values for which the LBL data
                                                                        !  are available for each species
   real(dp) :: lbl_data_tg(lbl_max_ntg,lbl_max_nxs,lbl_max_nsp)         !Array with the values of discrete temperature values for which the LBL data
                                                                        !  are available for each species and each mole fraction value
   real(dp),allocatable,dimension(:) :: lbl_deta,lbl_xeta               !Arrays containing all the center line positions and widths pertaining to the
                                                                        !  external LBL data files
   real(dp),allocatable,dimension(:,:,:,:) :: lbl_data                  !Array containing the input data read from the external LBL data files
   real(dp) :: lbl_eta_min,lbl_eta_max                                  !Minimum and maximum values of eta considered in the LBL solution (1/m)
   logical :: lbl_ceta_input                                            !If .true., the code expects the LBL data to be in the ceta format (default)
   logical :: lbl_xceta_input                                           !If .true., the code expects the LBL data to be in the x*ceta format
   logical :: lbl_kappa_input                                           !If .true., the code expects the LBL data to be in kappa_p format
   logical :: lbl_data_cm                                               !If .true., the code expects the LBL data to be in 1/cm rather than 1/m
   logical :: lbl_data_ready                                            !If .true., the external LBL data files are assumed to already be opened
                                                                        !  (so prepare_lbl_data does not need to be called)
   logical :: lbl_constant_ib                                           !If .true., assume that the Planck function is constant for each wavenumber interval
   logical :: lbl_binary_data                                           !If .true., expects the external LBL data to be unformatted (binary);
                                                                        !  the data should be formatted (e.g., in csv or dat file)
   logical :: print_xeta                                                !If .true., prints the xeta values during the RTE solution
   
contains   
   
   !====================================================================
   !Subroutine with the default parameters of the LBL method
   !====================================================================
   !(this does not need to be included in the main 
   !code if the user knows what they are doing)
   subroutine set_default_lbl_parameters
   
      !Set parameters of the LBL data files 
      !(so as to skip the reading of the files for all species)
      lbl_data_file = 'null'
      lbl_deta_file = 'null'
      lbl_data_prefix = ''
      lbl_data_ext = ''
      lbl_data_averaging = 'none'
      lbl_data_nxs = 0
      lbl_data_ntg = 0
      lbl_data_xs = 0._dp
      lbl_data_tg = 0._dp
      
      !Set flags
      lbl_ceta_input    = .false.
      lbl_xceta_input   = .false.
      lbl_kappa_input   = .false.
      lbl_data_cm       = .false.
      lbl_data_ready    = .false.
      print_xeta        = .false.
      lbl_constant_ib   = .false.
      lbl_binary_data   = .false.
      
      !Set soot parameters
      which_lbl_soot = 'chang1990'
      
      !Set bounding eta intervals such that no eta is skipped
      lbl_eta_min = -1._dp
      lbl_eta_max = 1.e20_dp
   
   endsubroutine set_default_lbl_parameters   
   
endmodule lbl_parameters
