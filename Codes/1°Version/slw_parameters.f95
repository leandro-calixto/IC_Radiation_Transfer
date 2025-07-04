!#######################################################################
!Module compiling SLW parameters
!#######################################################################
module slw_parameters

   use precision_parameters, only: dp
   implicit none
   integer,parameter :: slw_ns = 10, slw_nx = 10                        !These parameters give the maximum allowed sizes for arrays 
                                                                        !  related to the input ALBDF data. (ns -> number of species, 
                                                                        !  nx -> number of discrete mole fraction values)
   character(80) :: slw_mixture_method                                  !Specification of the mixture method
   character(80) :: slw_nonuniform_method                               !Specification of the method for treating
                                                                        !  non-uniform media
   character(80) :: slw_quadrature                                      !Speficiation of the quadrature used to define
                                                                        !the F values when needed
!   character(80) :: slw1_method                                         !Method to determine the gray gas coefficients
                                                                        !  in the SLW1 model
   character(80) :: slw_print_file                                      !Name of the file where the gray gas properties
                                                                        !  are printed, if requested
   character(80) :: albdf_file(slw_ns)                                  !Array with the name of the ALBDF file for each species
   character(80) :: albdf_info_file(slw_ns)                             !Array with the name of the file containing information
                                                                        !  on the ALBDF for each species
   character(80) :: slw1_approach                                       !Name of the approach to be used for computing the
                                                                        !  reference absorption and weighting coefficients in
                                                                        !  the SLW-1 model
   integer :: albdf_nx(slw_ns)                                          !Array containing the number of discrete mole 
                                                                        !  fractions in the ALBDF data file for each species
   integer :: albdf_ntg,albdf_ntb                                       !Number of discrete gas and source temperatures
                                                                        !  in the ALBDF data files
   integer :: albdf_nCj                                                 !Number of discrete absorption cross-sections in the
                                                                        !  ALBDF data files
   integer :: slw_ngas                                                  !Number of gray gases to be used in the RTE solution (this 
                                                                        !  parameter can be overuled when calling the subroutine)
   integer :: slw_ngas_per_species(10)                                  !Number of gray gases per individual species to be used
                                                                        !  in the RTE solution (this parameter is only used when
                                                                        !  mixing is done via multiple integration)
   integer :: albdf_nbands                                              !Number of discrete bands in the ALBDF data files
   integer :: slw1_ngases                                               !Number of gray gases used to compute the reference
                                                                        !  solution for determining the single-gray gas
                                                                        !  coefficients in the SLW1 model
   real(dp) :: albdf_cmin,albdf_cmax                                    !Minimum and maximum absorption cross-section values
                                                                        !  in the ALBDF data files
   real(dp) :: slw_cmin,slw_cmax                                        !Minimum and maximum absorption cross-section values
                                                                        !  for the application of the SLW model
   real(dp) :: slw_Fmin,slw_Fmax                                        !Minimum and maximum F values
   real(dp) :: slw_Ib_lbound,slw_Ib_ubound                              !Bounds (in 1/m) for the Ib integration
!   real(dp) :: slw1_length1,slw1_length2                                !Characteristic lengths used in the SLW1 model
   real(dp) :: albdf_x(slw_nx,slw_ns)                                   !Array with the mole fraction values for which the
                                                                        !  ALBDF data is available for each species
   real(dp) :: slw_pref,slw_Tref,slw_xsref(slw_ns)                      !Reference pressure, temperature and mole fraction of
                                                                        !  the participating species (needed for non-uniform SLW)
   real(dp) :: slw_superposition_theta                                  !Defines the value of the safety factor in the improved 
                                                                        !  multiple integration model
   real(dp) :: slw1_length(2)
   real(dp),allocatable,dimension(:) :: albdf_Tg,albdf_Tb               !Array with the gas and source temperature values 
                                                                        !  for which the ALBDF data is available
   real(dp),allocatable,dimension(:) :: albdf_cj,albdf_iarr             !Array with the absorption cross-section values for
                                                                        !  which the ALBDF data is available
   real(dp),allocatable,dimension(:) :: bslw_lbound,bslw_ubound         !Array with the lower and upper wavelength
                                                                        !  bounds (in m) of each band in the BLSW model
   real(dp),allocatable,dimension(:,:,:,:,:) :: ijk_albdf               !Array storing the mixed ALBDF (C and F values) at each
                                                                        !  grid point
   real(dp),allocatable,dimension(:,:,:,:,:,:) :: albdf_darr            !Arrays storing all the ALBDF data for all species
   logical :: slw_print_properties                                      !If .true., print the properties of each gray gas into
                                                                        !  file. The name of the file is slw_print_file
   logical :: slw_bound_Ib                                              !If .true., the blackbody intensity is not computed over
                                                                        !  the entire spectrum, but rather between slw_Ib_lbound
                                                                        !  and slw_Ib_ubound
   logical :: albdf_unformatted                                         !If .true., the ALBDF files are expected to be in
                                                                        !  binary (unformatted)
   logical :: albdf_loaded                                              !If .true., the ALBDF is already loaded, and does not need 
                                                                        !  to be read again
   logical :: albdf_inverted                                            !If .true., the ALBDF external data file comes as C(F); 
                                                                        !  if .false., F(C)
   contains
   !====================================================================
   !Subroutine with the default parameters of the SLW model
   !====================================================================
   !(this does not need to be included in the main 
   !code if the user knows what they are doing)
   subroutine set_default_slw_parameters
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use precision_parameters, only: big,small
      implicit none
   
      !-----------------------------------------------------------------
      !ALBDF file names
      !-----------------------------------------------------------------
      albdf_file = 'null'
      albdf_info_file = 'null'
   
      !-----------------------------------------------------------------
      !Maximum and minimum values for the C discretization
      !-----------------------------------------------------------------
      slw_cmin=0.0001_dp;     slw_cmax=1000._dp
      slw_Fmin=0._dp;         slw_Fmax=1._dp
      
      !-----------------------------------------------------------------
      !Reference values
      !-----------------------------------------------------------------
      slw_Tref = -1.e5_dp
      slw_pref = -1.e5_dp
      slw_xsref = -1.e5_dp
      
      !-----------------------------------------------------------------
      !SLW model parameters
      !-----------------------------------------------------------------
      slw_ngas                = 10
      slw_ngas_per_species    = 8
      slw_mixture_method      = 'multiplication'
      slw_nonuniform_method   = 'rank_correlated' 
      slw_quadrature          = 'Gauss-Legendre'
      slw_bound_Ib            = .false.
      slw_Ib_lbound           = small
      slw_Ib_ubound           = big
      slw_superposition_theta = big
      
      !-----------------------------------------------------------------
      !SLW1 model parameters
      !-----------------------------------------------------------------
      slw1_approach  = 'epsilon-epsilon'
      slw1_ngases  = 25
      slw1_length(1) = 1._dp
      slw1_length(2) = 2._dp
      
      !-----------------------------------------------------------------
      !Misc parameters
      !-----------------------------------------------------------------
      albdf_loaded         = .false.
      albdf_unformatted    = .true.
      albdf_inverted       = .false.
      slw_print_properties = .false.
      slw_print_file       = 'null'
      
   endsubroutine set_default_slw_parameters 
   
endmodule slw_parameters
