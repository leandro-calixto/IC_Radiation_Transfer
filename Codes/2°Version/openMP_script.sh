#! /bin/bash
program_name=program
fortran_version=gfortran
main_name=main_slw1.f95
run=y
numthreads=1

#Colors
RED='\033[0;31m'

#Start
echo Script has started

# Define number of OpenMP threads
echo Defining number of threads
export OMP_NUM_THREADS=$numthreads

# Compiling modules
echo Compiling modules
$fortran_version -c precision_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c constants.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c global_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c mesh.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c lbl_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c nbck_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c slw_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c wsgg_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c wbm_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c gg_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c mc_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c fvm_parameters.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c lbl_functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c nbck_functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c slw_functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c wsgg_functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c gg_functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c mc_functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c wbw_functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c wbm_functions.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c exact_solution.f95 -g -fcheck=all -Wall -fbacktrace -frecursive -fopenmp
$fortran_version -c fvm_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c mc_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c lbl_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c nbck_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c slw_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c wsgg_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c wbw_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c wbm_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c gg_routines.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp
$fortran_version -c validation.f95 -g -fcheck=all -Wall -fbacktrace -fopenmp

#Compiling the main program
echo Compiling main program
$fortran_version -c $main_name -g -fcheck=all -frecursive -Wall -fopenmp

# Generating the executable
echo Generating executable
$fortran_version *.o -o $program_name -fopenmp

# Running
rm *.o *.mod
if [ "$run" == "y" ]; then
   echo Running the executable
   ./$program_name
fi

# End
rm $program_name
echo Script has ended
