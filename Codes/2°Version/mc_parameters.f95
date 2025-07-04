module mc_parameters

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp,i8b
   use lbl_parameters, only: lbl_max_nsp,lbl_max_nxs,lbl_max_ntg
   implicit none
   character(100) :: which_spectral_mc                                  !Defines the spectral model to be used for the Monte Carlo solution
   character(100) :: which_mc_random                                    !Defines which random number generator to use
   character(100) :: which_mc_absorption,which_mc_emission              !Defines the schemes for absorption and emission
   integer(i8b) :: mc_total_rays,mc_rays_per_point                      !Defines the number of rays (either total or per grid point)
   integer :: mc_total_runs                                             !Defines the number of consecutive runs for a given problem
                                                                        !  (results are then averaged over the runs)
   integer :: mc_seed                                                   !Seed for the random number generation
   logical :: check_mc_energy                                           !If .true., checks for conservation of total energy in
                                                                        !  the MC method calculation
   logical :: randomize_cell_emission                                   !If .true., randomize position of emission within the cell

   !Parameters of the LBL random number relations
   character(200) :: lblrnr_file(lbl_max_nxs,lbl_max_nsp)               !Array with the names of the (unformatted) external 
                                                                        !  files containing all the data
   integer :: lblrnr_nsp                                                !Number of species
   integer :: lblrnr_ntg(lbl_max_nxs,lbl_max_nsp)                       !Number of discrete temperature values
   integer :: lblrnr_nxs(lbl_max_nsp)                                   !Number of discrete mole fraction values
   integer :: lblrnr_nln(lbl_max_nxs,lbl_max_nsp)                       !Number of discrete spectral intervals
   integer,parameter :: lblrnr_max_iterations = 1000                    !Maximum number of iterations for all inversions
   real(dp) :: lblrnr_etamax,lblrnr_etamin
   real(dp) :: lblrnr_xs(lbl_max_nxs,lbl_max_nsp)                       !Array with discrete mole fraction values for
                                                                        !  each species
   real(dp),allocatable,dimension(:,:,:) :: lblrnr_eta                  !Array with the spectral positions
   real(dp),allocatable,dimension(:,:,:) :: lblrnr_tg                   !Array with discrete temperature values

   real(dp),allocatable,dimension(:,:,:) :: lblrnr_kppi                 !Array with discrete values for the pressure-based
                                                                        !  Planck-mean absorption coefficient of each species
   real(dp),allocatable,dimension(:,:,:,:) :: lblrnr_kpi                !Array with discrete values for the pressure-based
                                                                        !  spectral absorption coefficient of each species
   real(dp),allocatable,dimension(:,:,:,:) :: lblrnr_Retai              !Array with discrete values for the spectral random number
                                                                        !  for each species
   real(dp),parameter :: lblrnr_tolerance = 1.e-4_dp                    !Tolerance for all inversions
   
   
contains
   !========================================================
   !Subroutine with the default parameters of the MC method
   !========================================================
   !(this does not need to be included in the main 
   !code if the user knows what they are doing)
   subroutine set_default_mc_parameters

      implicit none
      which_spectral_mc = 'GG'
      which_mc_random = 'numrep'
      which_mc_absorption = 'energy_partitioning'
      which_mc_emission = 'distributed'
      mc_total_rays = 1000
      mc_total_runs = 1
      mc_rays_per_point = 10
      mc_seed = 1234
      check_mc_energy = .true.
      randomize_cell_emission = .false.
      
      !Parameters of the data files with the random number relations
      lblrnr_file = 'null'
      lblrnr_nln = 0
      lblrnr_nxs = 0
      lblrnr_ntg = 0
      lblrnr_xs = 0._dp
      
   endsubroutine set_default_mc_parameters
   
endmodule mc_parameters
