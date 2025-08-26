!#######################################################################
!Functions related to the computational process
!#######################################################################
module comp_functions

   !====================================================================
   !Declaring variables for the module
   !====================================================================
   implicit none
   integer :: current_unit                                              !These parameters are necessary
   logical :: first_time=.true.                                         !   for the get_file_unit() function
   
contains
   !====================================================================
   !Subroutine to check if a logical argument is true; if it
   !is not, the program stops and an error message is printed
   !====================================================================
   subroutine assert(istrue,name_of_condition)
   
      implicit none
      logical,intent(in) :: istrue
      character(*),optional :: name_of_condition
      character(200) :: msg
      
      if ((present(name_of_condition)).and.(.not.istrue)) then
         msg = 'Condition '//trim(name_of_condition)//' not true'
         call shutdown(msg)
      endif
      if (.not.istrue) call shutdown('Assert condition not true')
         
   endsubroutine assert

   !====================================================================
   !Function to make a text string msg into bold face
   !====================================================================
   character(400) function bftext(msg)
   
      character(*),intent(in) :: msg
      character(3) :: cstr,rstr      
      rstr = '[0m'; cstr = '[1m'
      
      bftext = achar(27)//cstr//trim(msg)//achar(27)//rstr
   
   endfunction bftext

   !====================================================================
   !Subroutine for checking if an external file exists
   !====================================================================
   subroutine CheckFileExists(file_name)

      implicit none
      character(*),intent(in) :: file_name
      character(180) :: msg
      logical :: file_exists
      
      inquire(file=trim(file_name),exist=file_exists)
      if (.not.file_exists) then
         msg = 'File '//trim(file_name)//' does not exist'
         call shutdown(msg) 
      endif

   endsubroutine CheckFileExists

   !====================================================================
   !Subroutine for checking for memory allocation error
   !====================================================================
   subroutine CheckMemAlloc(array_name,error_value)
      
      implicit none
      character(*),intent(in) :: array_name
      character(180) :: msg
      integer,intent(in) :: error_value

      if (error_value.ne.0) then
         write(msg,'(a,a,a,i4)') 'Allocation error for array ',&
            trim(array_name),'. Error code:',error_value
         call shutdown(msg)
      endif
   
   endsubroutine CheckMemAlloc

   !====================================================================
   !Function to convert a string msg into a color
   !textcolor for outputs to the terminal
   !====================================================================
   character(400) function clrtext(msg,textcolor)
   
      character(*),intent(in) :: msg,textcolor
      character(3) :: rstr
      character(4) :: cstr
      
      rstr = '[0m'
      if (trim(textcolor).eq.'black')   cstr = '[30m'
      if (trim(textcolor).eq.'red')     cstr = '[31m'
      if (trim(textcolor).eq.'green')   cstr = '[32m'
      if (trim(textcolor).eq.'yellow')  cstr = '[33m'
      if (trim(textcolor).eq.'blue')    cstr = '[34m'
      if (trim(textcolor).eq.'magenta') cstr = '[35m'
      if (trim(textcolor).eq.'cyan')    cstr = '[36m'
      if (trim(textcolor).eq.'white')   cstr = '[37m'
      if (trim(textcolor).eq.'pink')    cstr = '[95m'
      
      clrtext = achar(27)//cstr//trim(msg)//achar(27)//rstr
   
   endfunction clrtext

   !====================================================================
   !Subroutine for printing to prompt when in debug mode
   !====================================================================
   subroutine dprint(msg,whitespace)
   
      use global_parameters, only: debug_mode
      implicit none
      character(*),intent(in) :: msg
      integer,optional :: whitespace
      integer :: wsp
      
      wsp = 0; if (present(whitespace)) wsp = whitespace
      if (debug_mode) call print_to_prompt(msg,wsp)

   endsubroutine dprint
   
   !====================================================================
   !Function for giving out file units
   !====================================================================
   integer function get_file_unit()

      if (first_time) then                                              !Set unit file of the first file as 20
         first_time = .false.
         current_unit = 20
      else
         current_unit = current_unit + 1                                !Ensure that other files have different units
      endif
      get_file_unit = current_unit

   endfunction get_file_unit

   !====================================================================
   !Function to determine the exact wall clock time
   !====================================================================
   real(dp) function get_wall_time()
      
      use precision_parameters, only: dp
      implicit none
      integer :: cnt,cnt_rate
      call system_clock(cnt,cnt_rate)
      get_wall_time = real(cnt,dp)/real(cnt_rate,dp)
      
   endfunction get_wall_time

   !====================================================================
   !Subroutine to print to prompt
   !====================================================================
   subroutine print_to_prompt(message,whitespace)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      character(*),intent(in) :: message
      character(10) :: out_format
      integer,optional :: whitespace
      integer :: char_length,w_length

      !-----------------------------------------------------------------
      !The default whitespace is zero
      !-----------------------------------------------------------------
      if (.not.present(whitespace)) then
         w_length=0
      else
         w_length = whitespace
      endif

      !-----------------------------------------------------------------
      !Define the format of the output to account for whitespace
      !-----------------------------------------------------------------
      char_length = len_trim(message)                                   !Length of the message
      char_length = char_length + w_length                              !Length of the message + whitespace
      if (char_length.eq.0) char_length = 1                             !For the total length is zero, set it to 1
      write(out_format,'(i5)') char_length                              !Write the total length as a string
      out_format = trim(adjustl(out_format))                            !Left align and trim the total length
      out_format = '(a'//trim(out_format)//')'                          !Define the final format for the output

      !-----------------------------------------------------------------
      !Print the message
      !-----------------------------------------------------------------
      write(*,out_format) trim(message)

      !-----------------------------------------------------------------
      !Flush the data
      !-----------------------------------------------------------------
      call flush(6)
      call flush(0)
      
   endsubroutine print_to_prompt
   
   !====================================================================
   !Subroutine to stop print an error message and stop the execution
   !====================================================================
   subroutine shutdown(message)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      character(*),intent(in),optional :: message   
      character(400) :: msg
      
      !-----------------------------------------------------------------
      !Set up message to be printed
      !-----------------------------------------------------------------
      if (present(message)) then
         msg = 'ERROR: '//trim(message)
      else
         msg = ''
      endif
      
      !-----------------------------------------------------------------
      !Print message and stop
      !-----------------------------------------------------------------
      call print_to_prompt(msg)
      stop
   
   endsubroutine shutdown

endmodule comp_functions

!#######################################################################
!Miscellaneous functions
!#######################################################################
module misc_functions

contains
   !====================================================================
   !Subroutine to swap two double precision elements
   !====================================================================
   subroutine dpswap(a,b,mask)
   
      use precision_parameters, only: dp
      implicit none
      real(dp),intent(inout) :: a,b
      real(dp) :: dum
      logical,intent(in),optional :: mask
      
      if (present(mask)) then
         if (.not.mask) return
      endif
      dum = a; a = b; b = dum
   
   endsubroutine dpswap
   
   !====================================================================
   !Subroutine to find the index in an array corresponding  
   !to the imin-to-last minimum value of the array
   !====================================================================
   integer recursive function find_imin(array,istart,iend,imin) &
      result (min_pos)

      use precision_parameters, only: dp
      implicit none
      integer,intent(in) :: istart,iend,imin
      integer :: ii,jj
      integer,allocatable,dimension(:) :: iprev
      real(dp),intent(in) :: array(:)
      real(dp) :: min_val
      allocate(iprev(1:imin))

      !Initial (very large) minimum value
      min_val = 1.e100_dp

      !Find the index corresponding to the absolute minimum
      if (imin.eq.1) then
         do jj=istart,iend
            if (array(jj).lt.min_val) then
               min_val = array(jj)
               min_pos = jj
            endif
         enddo      
      else
         !Find indexes of the previous minima
         do ii=imin-1,1,-1
            iprev(ii) = find_imin(array,istart,iend,ii)
         enddo
         
         !Find the correct index for the minimum
         min_loop: do jj=istart,iend
            if ((array(jj).lt.min_val)) then
               do ii=1,imin-1
                  if (jj.eq.iprev(ii)) cycle min_loop
               enddo      
               min_val = array(jj)
               min_pos = jj   
            endif
         enddo min_loop
      endif
      
      deallocate(iprev)
      
   endfunction find_imin

   !====================================================================
   !Subroutine to compute the average value of a given 3d array
   !====================================================================
   subroutine get_3d_average(array,x,y,z,nx,ny,nz,avg_value)
      
      use precision_parameters, only: dp
      implicit none
      integer,intent(in) :: nx,ny,nz
      integer :: i,j,jj,k,kk,nny,nnz
      real(dp),intent(in) :: array(0:nx+1,0:ny+1,0:nz+1),x(0:nx+1),&
         y(0:ny+1),z(0:nz+1)
      real(dp),intent(out) :: avg_value
      real(dp) :: dx,dy,dz,sum_array,volume
      
      !Correct ny and nz if necessary
      nny = ny; if (ny.lt.0) nny = 1
      nnz = nz; if (nz.lt.0) nnz = 1    
      
      sum_array = 0._dp
      volume = 0._dp      
      do k=1,nnz
         !Cell size in z direction
         if (nnz.eq.1) then                                              !Correction for 1d and 2d arrays
            kk = 0
            dz = 1._dp
         else
            kk = k
            if (kk.eq.1) then
               dz = 0.5_dp*(z(kk+1) + z(kk)) - z(kk-1)
            elseif (kk.eq.nz) then
               dz = z(kk+1) - 0.5_dp*(z(kk) + z(kk-1))
            else
               dz = 0.5_dp*(z(kk+1) - z(kk-1))
            endif
         endif
      
         do j=1,nny
            !Cell size in y direction
            if (nny.eq.1) then                                           !Correction for 1d arrays
               jj = 0
               dy = 1._dp
            else
               jj = j
               if (jj.eq.1) then
                  dy = 0.5_dp*(y(jj+1) + y(jj)) - y(jj-1)
               elseif (jj.eq.ny) then
                  dy = y(jj+1) - 0.5_dp*(y(jj) + y(jj-1))
               else
                  dy = 0.5_dp*(y(jj+1) - y(jj-1))
               endif
            endif

            do i=1,nx
               !Cell size in x direction
               if (i.eq.1) then
                  dx = 0.5_dp*(x(i+1) + x(i)) - x(i-1)
               elseif (i.eq.nx) then
                  dx = x(i+1) - 0.5_dp*(x(i) + x(i-1))
               else
                  dx = 0.5_dp*(x(i+1) - x(i-1))
               endif
               
               !Integrate the array in space
               sum_array = sum_array + array(i,jj,kk)*dx*dy*dz
               
               !Compute the total volume
               volume = volume + dx*dy*dz
            enddo
         enddo
      enddo
      
      avg_value = sum_array/volume
   
   endsubroutine get_3d_average
   
   !====================================================================
   !Subroutine to compute differences between two arrays
   !====================================================================
   subroutine get_differences(array_a,array_b,max_dif_rel_a,&
                              max_dif_norm_max_a,max_dif_norm_avg_a,&
                              avg_dif_rel_a,avg_dif_norm_max_a,&
                              avg_dif_norm_avg_a,te_a)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use precision_parameters, only: dp,small
      implicit none
      real(dp),intent(in) :: array_a(:),array_b(:)
      real(dp) :: avg_dif_rel_a,avg_dif_norm_avg_a,avg_dif_norm_avg_b,&
         avg_dif_norm_max_a,avg_dif_norm_max_b,avg_dif_rel_b,&
         max_dif_norm_avg_a,max_dif_norm_avg_b,max_dif_norm_max_a,&
         max_dif_norm_max_b,max_dif_rel_a,max_dif_rel_b,te_a,te_b                  
      integer :: i,n,n_a,n_b
      real(dp) :: aavg,amax,bavg,bmax
      real(dp),allocatable,dimension(:) :: aabs,babs,dif_abs,&
         dif_rel_a,dif_rel_b
      
      !-----------------------------------------------------------------
      !Preparing the arrays
      !-----------------------------------------------------------------
      !Finding the sizes of the arrays
      n_a = size(array_a)
      n_b = size(array_b)
      n   = n_a
      
      !If the arrays are not of the same size, stop
      if (n_a.ne.n_b) then
         write(*,*) "ERROR: Arrays are not of the same size to &
                     &compute differences"
         stop
      endif      
      
      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      allocate(aabs     (1:n))
      allocate(babs     (1:n))
      allocate(dif_abs  (1:n))
      allocate(dif_rel_a(1:n))
      allocate(dif_rel_b(1:n))      
      
      !-----------------------------------------------------------------
      !Extracting the maxima
      !-----------------------------------------------------------------
      !Computing absolute values of the arrays
      do i=1,n
         aabs(i) = dabs(array_a(i))
         babs(i) = dabs(array_b(i))
      enddo
      
      !Computing the maxima
      amax = maxval(aabs)
      bmax = maxval(babs)
      
      !-----------------------------------------------------------------
      !Extracting the averages
      !-----------------------------------------------------------------
      aavg = sum(aabs)/real(n,dp)
      bavg = sum(babs)/real(n,dp)
            
      !-----------------------------------------------------------------
      !Computing local differences
      !-----------------------------------------------------------------
      do i=2,n_a-1
         !Absolute difference
         dif_abs(i) = dabs(array_a(i) - array_b(i))
         
         !Relative differences
         dif_rel_a(i) = dif_abs(i)/(dabs(array_a(i))+small)
         dif_rel_b(i) = dif_abs(i)/(dabs(array_b(i))+small)
         
      enddo
      
      !-----------------------------------------------------------------
      !Finding maxima differences
      !-----------------------------------------------------------------
      max_dif_rel_a = maxval(dif_rel_a)*100._dp
      max_dif_rel_b = maxval(dif_rel_b)*100._dp
      
      max_dif_norm_max_a = maxval(dif_abs)*100._dp/amax
      max_dif_norm_max_b = maxval(dif_abs)*100._dp/bmax
      
      max_dif_norm_avg_a = maxval(dif_abs)*100._dp/aavg
      max_dif_norm_avg_b = maxval(dif_abs)*100._dp/bavg
   
      !-----------------------------------------------------------------
      !Finding averaged differences
      !-----------------------------------------------------------------
      avg_dif_rel_a = sum(dif_rel_a)*100._dp/real(n_a,dp)
      avg_dif_rel_b = sum(dif_rel_b)*100._dp/real(n_b,dp)
      
      avg_dif_norm_max_a = sum(dif_abs)*100._dp/(amax*real(n,dp))
      avg_dif_norm_max_b = sum(dif_abs)*100._dp/(bmax*real(n,dp))
      
      avg_dif_norm_avg_a = sum(dif_abs)*100._dp/(aavg*real(n,dp))
      avg_dif_norm_avg_b = sum(dif_abs)*100._dp/(bavg*real(n,dp))
      
      !-----------------------------------------------------------------
      !Total error (proposed by Webb)
      !-----------------------------------------------------------------
      te_a = sum(dif_abs)*100._dp/sum(aabs)
      te_b = sum(dif_abs)*100._dp/sum(babs)

   endsubroutine get_differences
   
   !====================================================================
   !Subroutine to get a random number
   !====================================================================
   subroutine get_random_number(rnd_num)
      
      use precision_parameters, only: dp
      implicit none
      real(dp),intent(out) :: rnd_num
      call random_number(rnd_num)
      
   endsubroutine get_random_number
   
   !====================================================================
   !Function to find the position of the
   !last non-negative value in an array
   !===================================================================
   integer function last_non_negative_position(x_array)
   
      use precision_parameters, only: dp
      implicit none
      integer :: i,darray,uarray
      real(dp),intent(in) :: x_array(:)
      
      !Find the dimension bounds of the array
      uarray = ubound(x_array,1)
      darray = lbound(x_array,1)
      
      !Check the array position by position
      do i=darray,uarray
         if (x_array(i).lt.0._dp) exit
      enddo
      
      last_non_negative_position = i - 1
   
   endfunction last_non_negative_position
   
   !====================================================================
   !Function to find the index of the last nonzero in an array
   !====================================================================
   integer function nonzero_ubound(lookup)
   
      use precision_parameters, only: dp
      implicit none
      integer :: i,darray,uarray
      real(dp),intent(in) :: lookup(:)
      
      !Find the dimension bounds of the array
      uarray = ubound(lookup,1)
      darray = lbound(lookup,1)
      
      !Check the array position by position
      do i=darray,uarray
         if (lookup(i).eq.0._dp) exit
      enddo
      nonzero_ubound = i - 1
   
   endfunction nonzero_ubound

   !====================================================================
   !Subroutine to check if all members in the array 
   !have the same value within a given tolerance
   !====================================================================
   subroutine same_array_members(array,npts,tol,okay)
   
      use precision_parameters, only: dp
      implicit none
      integer,intent(in) :: npts
      integer :: i
      real(dp),intent(in) :: array(1:npts),tol
      real(dp) :: val1,val2
      logical,intent(out) :: okay
      
      okay = .true.
      val1 = array(1)
      do i=2,npts
         val2 = array(i)
         if (abs((val1-val2)/val1).gt.tol) then
            okay = .false.
            exit
         endif
         val1 = val2
      enddo
   
   endsubroutine same_array_members

   !====================================================================
   !Subroutine to define the random number seed
   !====================================================================
   subroutine set_random_seed(seed_value)
      
      implicit none
      integer,intent(in) :: seed_value
      integer :: n
      integer,allocatable :: seed(:)
      
      call random_seed(size=n)
      allocate(seed(n))
      seed = seed_value
      call random_seed(put=seed)
      deallocate(seed)
   
   endsubroutine set_random_seed
   
endmodule misc_functions

!#######################################################################
!Physical functions
!#######################################################################
module physical_functions

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp,small
   implicit none

contains

   !====================================================================
   !Function to assess a property at an interface 
   !between two grid cells, given that only the values 
   !of the property at the grid points are known
   !Inputs:
   !  p1_val, p2_val: values of the property at the grid points 
   !                  of the cells that neighbor the face
   !  dx1f: distance from point 1 (to the left of the face) to the face
   !  dx12: distance between grid points
   !  inteface_approach: either 'linear' or 'harmonic'
   !====================================================================
   real(dp) function get_interface_property(p1_val,p2_val,dx1f,dx12,&
                                            interface_approach)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      character(12),intent(in) :: interface_approach
      real(dp),intent(in) :: dx12,dx1f,p1_val,p2_val
      real(dp) :: f,r1,r2
      
      !-----------------------------------------------------------------
      !Define interface value according to interface_approach string
      !-----------------------------------------------------------------
      f = dx1f/(dx12 + small)
      selectcase(trim(interface_approach))
      case('linear')
         get_interface_property = f*p2_val + (1._dp - f)*p1_val
      case('harmonic')
         r1 = 1._dp/(p1_val + small); r2 = 1._dp/(p2_val + small)
         get_interface_property = 1._dp/(f*r1 + (1._dp - f)*r2 + small)
      case default
         call shutdown('get_interface_property: &
                        &interface_approach undefined')
      endselect
      
   endfunction get_interface_property

   !====================================================================
   !Function to compute the spectral blackbody radiation intensity 
   !at a given wavenumber value (for a refractive index of 1)
   !Input parameters:
   !  tmp: temperature
   !  x_eta: wavenumber [1/m]
   !  emissive_power: if true, compute the spectral blackbody emissive
   !                  power instead (default: false)
   !  wavelength: if true, assume the spectral input is in wavelength
   !              [m], and output in the Planck function in the
   !              corresponding units
   !Output:
   !  Planck_function: in W/m^3 or W/m
   !====================================================================
   real(dp) function Planck_function(tmp,x_eta,emissive_power,&
      wavelength)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: c_planck1,c_planck2,pi
      implicit none
      logical :: compute_E,lambda_input
      logical,intent(in),optional :: emissive_power,wavelength
      real(dp),intent(in) :: tmp,x_eta
      real(dp) :: mult_factor
      
      !-----------------------------------------------------------------
      !Optional parameters
      !-----------------------------------------------------------------
      compute_E = .false.
      if (present(emissive_power)) compute_E = emissive_power
      
      lambda_input = .false.
      if (present(wavelength)) lambda_input = wavelength
      
      !-----------------------------------------------------------------
      !Computing the function
      !-----------------------------------------------------------------
      mult_factor = 1._dp; if (compute_E) mult_factor = pi
      if (tmp.eq.0._dp) then
         !Assure the function to be zero for null temperatures
         Planck_function = 0._dp
         
      else
         if (lambda_input) then
            Planck_function = 2._dp*mult_factor*c_planck1/&
               ((x_eta**5)*(DEXP(c_planck2/(tmp*x_eta)) - 1._dp + small))
         
         else
            !Eq. (1.19) [Howell et al., 2010]
            Planck_function = 2._dp*mult_factor*c_planck1*(x_eta**3._dp)/&
                        (DEXP(c_planck2*x_eta/tmp) - 1._dp + small)
         endif
      endif

   endfunction Planck_function


   !====================================================================
   !Function to compute the total blackbody radiation intensity 
   !at a given temperature (for a refractive index of 1)
   !Input parameters:
   !  tmp: temperature [K]
   !  emissive_power: if true, compute the blackbody emissive power 
   !                  instead (default: false)
   !Output:
   !  Planck_function: in W/m^2
   !====================================================================
   real(dp) function Ib_function(tmp,emissive_power)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: pi,sigrpi
      implicit none
      logical :: compute_E
      logical,intent(in),optional :: emissive_power
      real(dp),intent(in) :: tmp
      real(dp) :: mult_factor
      
      !-----------------------------------------------------------------
      !Optional parameters
      !-----------------------------------------------------------------
      compute_E = .false.
      if (present(emissive_power)) compute_E = emissive_power
      
      !-----------------------------------------------------------------
      !Computing the function
      !-----------------------------------------------------------------
      mult_factor = 1._dp; if (compute_E) mult_factor = pi
      Ib_function = sigrpi*mult_factor*(tmp**4)

   endfunction Ib_function

   !====================================================================
   !Function to compute the blackbody radiation emission between zero 
   !and a given wavelength (if eta_in = .true., the input is in 
   !wavenumber [1/m]; otherwise, it is in wavelength [m])
   !====================================================================
   real(dp) function bb_emission_frac(wvqt,tmp,eta_in)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use constants, only: c_planck2,pi
      use precision_parameters, only: small
      implicit none
      real(dp),intent(in) :: wvqt,tmp
      real(dp) :: F_bnd,F_bnd0,lambda,m,sum_dif,sum_tol,zeta
      logical,optional :: eta_in
      logical :: eta_to_lambda
      
      !-----------------------------------------------------------------
      !Halt if the temperature is too small
      !-----------------------------------------------------------------
      if (tmp.le.small) call shutdown('bb_emission_frac: T = 0')
   
      !-----------------------------------------------------------------
      !Set up logical parameter
      !-----------------------------------------------------------------
      eta_to_lambda = .false.
      if (present(eta_in)) eta_to_lambda = eta_in

      !-----------------------------------------------------------------
      !Initializing parameter
      !-----------------------------------------------------------------
      sum_tol = 1.e-10_dp                                               !Tolerance for converging the posterior sum
      sum_dif = 10000._dp                                               !Initializing the difference between consecutive terms of the posterior sum
      m       = 0._dp                                                   !Initializing the counting term of the posterior sum
      F_bnd   = 0._dp                                                   !Initializing the current function value as zero (necessary for the sum)
      F_bnd0  = 1000._dp                                                !Initializing the old function value as a large value (necessary to compute the first difference)
      
      !-----------------------------------------------------------------
      !Computing the blackbody emission fraction
      !following Eq. (1.25) [Modest, 2013]
      !-----------------------------------------------------------------
      lambda = wvqt; if (eta_to_lambda) lambda = 1._dp/wvqt
      zeta = c_planck2/(lambda*tmp+small)
      do while (sum_dif.ge.sum_tol)
         m = m + 1._dp
         F_bnd = F_bnd + (dexp(-m*zeta)/(m**4._dp))*(6._dp + &
                  6._dp*m*zeta + 3._dp*(m*zeta)**2._dp + &
                  (m*zeta)**3._dp)
         sum_dif = dabs((F_bnd - F_bnd0)/(F_bnd + small))               !Difference for testing the convergence of the sum
         F_bnd0 = F_bnd                                                   !Updating the old function value
      enddo
         bb_emission_frac = (15._dp/(pi**4._dp))*F_bnd
         
   endfunction bb_emission_frac
   
   !====================================================================
   !Function to compute the blackbody radiation emission 
   !fraction within a given wavelength
   !====================================================================
   real(dp) function F_band(lambda,lambda_array,tmp)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use misc_functions, only: nonzero_ubound
      implicit none
      integer :: i,n_array,n_max,n_min
      real(dp),intent(in) :: tmp,lambda_array(:),lambda
      real(dp) :: F
      
      !-----------------------------------------------------------------
      !Finding bounds of the spectral array
      !-----------------------------------------------------------------
      n_min      = lbound(lambda_array,1)                               !Lower bound of the wavelength array
      n_max      = nonzero_ubound(lambda_array)                         !Upper bound of the wavelength array      
      n_array    = n_max-n_min+1                                        !Size of the wavelength array

      !-----------------------------------------------------------------
      !Computing the blackbody radiation fraction
      !-----------------------------------------------------------------
      if (n_array.eq.0) then                                            !For a single interval, ensure that F = 1
         F = 1._dp
      
      elseif (lambda.eq.lambda_array(n_min)) then
         F = bb_emission_frac(lambda,tmp)                               !Defining F for the first spectral interval
         
      elseif (lambda.eq.lambda_array(n_max+1)) then
         F = 1._dp - bb_emission_frac(lambda_array(n_max),tmp)          !Defining F for the last spectral interval
         
      else
         !Finding the index of the current variable inside the array
         do i=n_min+1,n_max-1   
            if (lambda.eq.lambda_array(i)) exit
         enddo
         
         !Defining F for the remaining spectral intervals
         F = bb_emission_frac(lambda,tmp) - &
             bb_emission_frac(lambda_array(i-1),tmp)             
      endif
      F_band = F
      
   endfunction F_band

   !====================================================================
   !Function to compute the total wall emissivity or absorptivity
   !====================================================================
   real(dp) function get_total_emissivity(eps_array,lmbd_array,tmp)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use misc_functions, only: nonzero_ubound
      implicit none
      integer :: i,eps_down,eps_up
      real(dp),intent(in) :: eps_array(:),lmbd_array(:),tmp
      real(dp) :: sum_eps
      
      !-----------------------------------------------------------------
      !Array operations
      !-----------------------------------------------------------------
      eps_down = lbound(eps_array,1)
      eps_up = nonzero_ubound(eps_array)

      !-----------------------------------------------------------------
      !Computing the total emissivity/absorptivity
      !-----------------------------------------------------------------
      sum_eps = 0._dp
      do i=eps_down,eps_up
         sum_eps = sum_eps + &
                  eps_array(i)*F_band(lmbd_array(i),lmbd_array,tmp)
      enddo
      get_total_emissivity = sum_eps
   
   endfunction get_total_emissivity

   !====================================================================
   !Function to compute the emissivity variance in wavelength space
   !====================================================================
   real(dp) function get_emissivity_variance(eps_array,lmbd_array,tmp,&
                                             eps_gray)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use misc_functions, only: nonzero_ubound
      implicit none
      integer :: i,eps_down,eps_up
      real(dp),intent(in) :: eps_array(:),lmbd_array(:),tmp,eps_gray
      real(dp) :: sum_eps
      
      !-----------------------------------------------------------------
      !Array operations
      !-----------------------------------------------------------------
      eps_down = lbound(eps_array,1)
      eps_up = nonzero_ubound(eps_array)

      !-----------------------------------------------------------------
      !Computing the total emissivity/absorptivity
      !-----------------------------------------------------------------
      sum_eps = 0._dp
      do i=eps_down,eps_up
         sum_eps = sum_eps + ((eps_array(i) - eps_gray)**2._dp)*&
            F_band(lmbd_array(i),lmbd_array,tmp)
      enddo
      get_emissivity_variance = sum_eps
   
   endfunction get_emissivity_variance

   !====================================================================
   !Function to compute the spectral blackbody radiative 
   !intensity integrated from eta_minus to eta_plus
   !====================================================================
   real(dp) function integrate_Ib(ttmp,eta_minus,eta_plus)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      real(dp),intent(in) :: ttmp,eta_minus,eta_plus
      real(dp) :: deta,ib_diff,curr_eta,old_ib_sum,ib_sum
      real(dp) :: ib_tol
   
      !-----------------------------------------------------------------
      !Check input parameters
      !-----------------------------------------------------------------
      if (eta_plus.lt.eta_minus) call shutdown('integrate_Ib: &
         &Ib integration with eta_plus < eta_minus')
      
      !-----------------------------------------------------------------
      !Define initial values
      !-----------------------------------------------------------------
      ib_tol = 1.e-4_dp                                                 !Tolerance for the integral
      old_ib_sum = -10000._dp
      deta = (eta_plus - eta_minus)/10._dp

      int_ib_loop: do
         curr_eta = eta_minus
         ib_sum = 0._dp
         do while(curr_eta.le.eta_plus)
            ib_sum = ib_sum + Planck_function(ttmp,curr_eta)*deta       !Compute the Ib integration
            curr_eta = curr_eta + deta                                  !Update curr_eta
         enddo
   
         !Check for convergence
         ib_diff = abs(ib_sum - old_ib_sum)/(ib_sum + small)
         if (ib_diff.lt.ib_tol) then
            exit int_ib_loop
         else
            old_ib_sum = ib_sum
            deta = 0.5_dp*deta
            cycle int_ib_loop
         endif
      enddo int_ib_loop
      integrate_Ib = ib_sum

   endfunction integrate_Ib

   !====================================================================
   !Function to convert from spectral absorption cross-section [m2/mol] 
   !to pressure spectral absorption coefficient [1/(m.atm)]
   !(if ktoC = .true., convert the other way around)
   !====================================================================
   real(dp) function Ceta_to_keta(input,T,ktoC)
   
      use constants, only: atm_pa,Ru
      logical,intent(in), optional :: ktoC
      logical :: Ctok
      real(dp),intent(in) :: input,T
      real(dp) :: mol_density
      
      Ctok = .true.; if (present(ktoC)) Ctok = .not.ktoC
      mol_density = atm_pa/(T*Ru)
      if (Ctok) then
         Ceta_to_keta = input*mol_density
      else
         Ceta_to_keta = input/mol_density
      endif
   
   endfunction Ceta_to_keta

   !====================================================================
   !Subroutine to compute the real and imaginary part of the complex 
   !index of refraction of soot according to the Chang and 
   !Charalampopoulos (1990) correlation
   !x_eta is expected in 1/m
   !===================================================================
   subroutine chang_index_of_refraction(x_eta,nn,kk)
   
      !-----------------------------------------------------------------
		!Declaration of variables
		!-----------------------------------------------------------------
      implicit none
      real(dp),intent(in) :: x_eta
      real(dp),intent(out) :: nn,kk
      real(dp) :: x_lambda,lnlambda
   
      !Convert from m-1 to mum
      x_lambda = 1.e6_dp/x_eta
      
      !-----------------------------------------------------------------
      !Main calculation
      !-----------------------------------------------------------------
      if ((x_lambda.lt.0.4_dp).or.(x_lambda.gt.30_dp)) then             !Restrict to the bounds of the correlation
         nn = 0._dp
         kk = 0._dp
      else
         !Ln(x_lambda)
         lnlambda = log(x_lambda)
      
         !Real part
         nn = 1.811_dp + 0.1263_dp*lnlambda + 0.027_dp*lnlambda**2._dp + &
               0.0417_dp*lnlambda**3._dp   
      
         !Imaginary part
         kk = 0.5821_dp + 0.1213_dp*lnlambda + 0.2309_dp*lnlambda**2._dp - &
               0.01_dp*lnlambda**3._dp
      endif
   
   endsubroutine chang_index_of_refraction

   !====================================================================
   !Function that outputs the spectral absorption coefficient of soot
   !at a given wavenumber x_eta (1/m) and for a volume fraction fv, 
   !following Rayleigh theory. By default, the Chang and 
   !Charalampopoulos (1990) correlation is used.
   !====================================================================
   real(dp) function kappa_eta_soot(x_eta,fv)

      use constants, only: pi
      use precision_parameters, only: small
      implicit none
      real(dp),intent(in) :: fv,x_eta
      real(dp) :: denum,nsoot,nsoot2,num,ksoot,ksoot2

      call chang_index_of_refraction(x_eta,nsoot,ksoot)
      nsoot2 = nsoot**2; ksoot2 = ksoot**2
      
      num = 36._dp*pi*fv*nsoot*ksoot*x_eta
      denum = (nsoot2 - ksoot2 + 2._dp)**2 + 4._dp*nsoot2*ksoot2
      
      kappa_eta_soot = num/(denum + small)

   endfunction kappa_eta_soot
   
endmodule physical_functions


!#######################################################################
!Mathematical functions
!#######################################################################
module math_functions

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp,small
   implicit none

contains

   !====================================================================
   !Sort array arr in ascending order using Heapsort
   !(adapted from Numerical Recipes)
   !====================================================================
   subroutine sort_heap(arr)
      use misc_functions, only: dpswap
      implicit none
      real(dp),dimension(:),intent(inout) :: arr
      integer :: i,n
   
      n=size(arr)
      do i=n/2,1,-1
         call sift_down(i,n)
      end do
      do i=n,2,-1
         call dpswap(arr(1),arr(i))
         call sift_down(1,i-1)

      end do
   
   contains
      subroutine sift_down(l,r)
      integer, intent(in) :: l,r
      integer :: j,jold
      real(dp) :: a

      a=arr(l)
      jold=l
      j=l+l
      do
         if (j > r) exit
         if (j < r) then
            if (arr(j) < arr(j+1)) j=j+1
         end if
         if (a >= arr(j)) exit
         arr(jold)=arr(j)
         jold=j
         j=j+j
      end do
      arr(jold)=a
   endsubroutine sift_down
   endsubroutine sort_heap

   !====================================================================
   !Factorial function
   !====================================================================
   real(dp) function factorial(n)
      
      use comp_functions, only: shutdown
      implicit none
      integer, intent(in) :: n
      integer :: i

      if (n.lt.0) call shutdown('factorial:&
                                &undefined for negative integers')
      factorial = 1._dp
      if ((n.eq.0).or.n.eq.1) return
      do i=2,n
         factorial = factorial*real(i,dp)
      enddo
   
   end function factorial
   
   !====================================================================
   !Subroutine that performs a gauss elimination to matrix 
   !a_matrix, giving the final solution in vector sol_vector
   !====================================================================
   subroutine gauss_elimination(a_matrix,sol_vector)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      real(dp) :: a_matrix(:,:),sol_vector(:)
      real(dp) :: f
      integer :: ii,jj,hh,kk,mm,nn
      
      !-----------------------------------------------------------------
      !Pre-processing
      !-----------------------------------------------------------------
      !Find number of rows and columns of the a_matrix
      mm = size(a_matrix,1)
      nn = size(a_matrix,2)
   
      !-----------------------------------------------------------------
      !Constructing the row-echelon matrix
      !-----------------------------------------------------------------
      !Initialization (part I)
      hh=1                                                              !Pivot row
      kk=1                                                              !Pivot column
      do while ((hh.le.mm).and.(kk.le.nn))
         if (a_matrix(hh,kk).eq.0._dp) then
            !No pivot in the column, pass to the next
            kk = kk + 1
         else
            !Divide whole pivot row by pivot value
            f = 1._dp/a_matrix(hh,kk)
            do jj=kk,nn
               a_matrix(hh,jj) = a_matrix(hh,jj)*f
            enddo
         
            !Loop over all rows
            do ii=hh+1,mm
               f = a_matrix(ii,kk)/a_matrix(hh,kk)
            
               !Loop over all columns
               do jj=kk,nn
                  a_matrix(ii,jj) = a_matrix(ii,jj) - a_matrix(hh,jj)*f
               enddo
            enddo
               
            !Update pivot row and column
            hh = hh + 1
            kk = kk + 1
         
         endif
      enddo
   
      !-----------------------------------------------------------------
      !Finishing up the Gauss elimination
      !-----------------------------------------------------------------
      !Initialization (part II)
      hh=mm                                                             !Pivot row
      kk=nn-1                                                           !Pivot column
      do while((hh.ge.1).and.(kk.ge.1))
         !Loop over all rows (backwards)
         do ii=hh-1,1,-1
            f = a_matrix(ii,kk)/a_matrix(hh,kk)
         
            !Loop over all columns (backwards)
            do jj=1,nn
               a_matrix(ii,jj) = a_matrix(ii,jj) - a_matrix(hh,jj)*f
            enddo   
         enddo
      
         !Update pivot row and column
         hh = hh - 1
         kk = kk - 1
      enddo
      
      !Mounting the solution vector
      do ii=1,mm
         sol_vector(ii) = a_matrix(ii,nn)
      enddo
   
   endsubroutine gauss_elimination   

   !====================================================================
   !Double precision Heaviside function
   !====================================================================
   real(dp) function heaviside_dp(xval)
   
      implicit none
      real(dp),intent(in) :: xval
      heaviside_dp = 0._dp; if (xval.gt.0._dp) heaviside_dp = 1._dp
   
   endfunction heaviside_dp
   
   !====================================================================
   !Kronecker delta function
   !====================================================================
   integer function kron_delta(ii,jj)
   
      implicit none
      integer,intent(in) :: ii,jj
      kron_delta = 1; if (ii.ne.jj) kron_delta = 0
   
   endfunction kron_delta
     
   !====================================================================
   !Subroutine to update the maximum value of scalar 
   !max_val given a new scalar value new_val
   !====================================================================
   subroutine update_max(new_val,max_val)
   
      implicit none
      real(dp),intent(in) :: new_val
      real(dp),intent(inout) :: max_val
      if (new_val.gt.max_val) max_val = new_val
      
   endsubroutine update_max

   !====================================================================
   !Subroutine to update the minimum value of scalar 
   !min_val given a new scalar value new_val
   !====================================================================
   subroutine update_min(new_val,min_val)
   
      implicit none
      real(dp),intent(in) :: new_val
      real(dp),intent(inout) :: min_val
      if (new_val.lt.min_val) min_val = new_val
      
   endsubroutine update_min
 
   !====================================================================
   !Function to locate a value within an array with a bisection method
   !====================================================================
   !Returns a value j such that xval is between xarray(j) and 
   !xarray(j+1); npts is the number of points in the array.
   !Returns j=0 or j=npts if xval is out of bounds.
   !(adapted from Numerical Recipes)
   integer function locate(xarray,xval,npts)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer :: npts
      integer :: jl,ju,jm
      real(dp),intent(in) :: xarray(npts),xval
      logical :: ascnd

      !-----------------------------------------------------------------
      !Find if the array is in ascending or descending order
      !-----------------------------------------------------------------
      if (xarray(npts).ge.xarray(1)) ascnd = .true.
      if (xarray(npts).lt.xarray(1)) ascnd = .false.
      
      !-----------------------------------------------------------------
      !Main loop
      !-----------------------------------------------------------------
      jl = 0                                                            !Initialize lower
      ju = npts + 1                                                     !and upper limits
      do
         if ((ju-jl).le.1) exit
         jm=(ju+jl)/2                                                   !Compute midpoint
         if (ascnd.eqv.(xval.ge.xarray(jm))) then                       !and replace the 
            jl=jm                                                       !upper or lower
         else                                                           !limits accordingly
            ju=jm
         end if
      end do
      
      !-----------------------------------------------------------------
      !Set the output
      !-----------------------------------------------------------------
      if (xval.eq.xarray(1)) then
         locate=1
      elseif (xval.eq.xarray(npts)) then
         locate=npts-1
      else
         locate=jl
      endif
   
   endfunction locate

   !====================================================================
   !Function to locate the position of the array minimum 
   !====================================================================
   integer function iminloc(arr)
   
      implicit none
      real(dp),intent(in) :: arr(:)
      integer :: imin(1)
      imin = minloc(arr(:))
      iminloc = imin(1)
   
   endfunction iminloc

   !====================================================================
   !Function to perform a linear interpolation
   !====================================================================
   real(dp) function linint(xarray,yarray,xtarget,npts)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer,intent(in) :: npts
      integer :: idown,iup
      real(dp),intent(in) :: xarray(:),yarray(:),xtarget
      real(dp) :: x1,x2,y1,y2
      
      !-----------------------------------------------------------------
      !Particularize the case where the 
      !arrays have only a single member
      !-----------------------------------------------------------------
      if (npts.eq.1) then
         linint = yarray(1)
      else
      !-----------------------------------------------------------------
      !For npts>1, carry out the interpolation
      !-----------------------------------------------------------------
         !Find upper and lower indexes
         idown = max(1,min(npts-1,locate(xarray,xtarget,npts)))
         iup = idown + 1
      
         !Interpolating for ytarget
         x1 = xarray(idown); x2 = xarray(iup)
         y1 = yarray(idown); y2 = yarray(iup)
         if (abs(x2-x1).le.spacing(x1)) then 
            linint = y1
         else
            linint = y1 + (y2 - y1)*(xtarget - x1)/(x2 - x1 + small)
         endif
      endif
      
   endfunction

!   subroutine polint(xa,ya,x,y,dy)

!      implicit none
!      real(dp), dimension(:), intent(in) :: xa,ya
!      real(dp), intent(in) :: x
!      real(dp), intent(out) :: y,dy
!      integer :: m,n,ns
!      real(sp), dimension(size(xa)) :: c,d,den,ho
!n=assert_eq(size(xa),size(ya),’polint’)
!c=ya initialize the tableau of c’s and d’s.
!d=ya
!ho=xa-x
!ns=iminloc(abs(x-xa)) find index ns of closest table entry.
!y=ya(ns) this is the initial approximation to y.
!ns=ns-1
!do m=1,n-1 for each column of the tableau,
!den(1:n-m)=ho(1:n-m)-ho(1+m:n) we loop over the current c’s and d’s and up-
!date them.if (any(den(1:n-m) == 0.0)) &
!call nrerror(’polint: calculation failure’)
!this error can occur only if two input xa’s are (to within roundoff) identical.
!den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
!d(1:n-m)=ho(1+m:n)*den(1:n-m) here the c’s and d’s are updated.
!c(1:n-m)=ho(1:n-m)*den(1:n-m)
!if (2*ns < n-m) then after each column in the tableau is completed, we decide
!which correction, c or d, we want to add to our accu-
!mulating value of y, i.e., which path to take through
!the tableau—forking up or down. we do this in such a
!way as to take the most “straight line” route through the
!tableau to its apex, updating ns accordingly to keep track
!of where we are. this route keeps the partial approxima-
!tions centered (insofar as possible) on the target x. the
!last dy added is thus the error indication.
!dy=c(ns+1)
!else
!dy=d(ns)
!ns=ns-1
!end if
!y=y+dy
!end do
!end subroutine polint

   !====================================================================
   !Subroutine for polynomial interpolation
   !====================================================================
   !Given arrays xarray and yarray, and given a value x, this subroutine
   !returns a value y and an error estimate (the latter optional),
   !which is obtained via polynomial interpolation using a polynomial of
   !degree deg. Adapted from Numerical Recipes
   subroutine polint(xarray,yarray,x,y,deg,error)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,shutdown
      implicit none
      integer,intent(in) :: deg
      integer :: il,iu,j,m,n,ns
      real(dp),dimension(:), intent(in) :: xarray,yarray
      real(dp),intent(in) :: x
      real(dp),intent(out) :: y
      real(dp),intent(out),optional :: error
      real(dp) :: dy
      real(dp),dimension(deg+1) :: c,d,den,ho,xa,ya

      call assert(size(xarray).eq.size(yarray),'polint')
      n = size(xarray)
!write(*,*) 'n',n,xarray(1),xarray(n),yarray(1),yarray(n)
!      if (n.eq.1) then
!         y = yarray(1)
!         if (present(error)) error = 0._dp
!         return
!      endif

      if ((deg+1).lt.n) then
         j = locate(xarray,x,n)
         il = min(max(j-ceiling(real(deg-1,dp)/2._dp),1),n-deg)
         iu = il + deg
         do m=il,iu
            xa(m-il+1) = xarray(m); ya(m-il+1) = yarray(m)
         enddo
         n = deg + 1
      else
         xa = xarray; ya = yarray
      endif

      !Initialize the tableau of c’s and d’s.
      c=ya; d=ya
      ho=xa-x
      
      !Find index ns of closest table entry.
      ns=iminloc(abs(x-xa))
      y=ya(ns)                                                          !Initial approximation to y.
      ns=ns-1

      do m=1,n-1                                                        !For each column of the tableau,
         den(1:n-m)=ho(1:n-m)-ho(1+m:n)                                 !we loop over the current c’s and d’s and up-
         if (any(den(1:n-m) == 0.0)) &                                  !date them.
            call shutdown('polint: calculation failure')                !this error can occur only if two input xa’s are (to within roundoff) identical.
         den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
         d(1:n-m)=ho(1+m:n)*den(1:n-m)                                  !Here the c’s and d’s are updated.
         c(1:n-m)=ho(1:n-m)*den(1:n-m)
         if (2*ns < n-m) then                                           !After each column in the tableau is completed, we decide
            dy=c(ns+1)                                                  !which correction, c or d, we want to add to our accu-
         else                                                           !mulating value of y, i.e., which path to take through      
            dy=d(ns)                                                    !the tableau—forking up or down. we do this in such a
            ns=ns-1                                                     !way as to take the most “straight line” route through the
         endif                                                          !tableau to its apex, updating ns accordingly to keep track
         y=y+dy                                                         !of where we are. this route keeps the partial approxima-
      enddo                                                             !tions centered (insofar as possible) on the target x. the
                                                                        !last dy added is thus the error indication.
      if (present(error)) error = dy
   
   endsubroutine polint

   !====================================================================
   !Subroutine to compute the average value of a given array
   !====================================================================
   subroutine get_average(array,avg_value)
   
      implicit none
      real(dp),intent(in) :: array(:)
      real(dp),intent(out) :: avg_value
      integer :: n_array
      
      n_array = size(array)
      avg_value = sum(array)/real(n_array,dp)      
   
   endsubroutine get_average

   !====================================================================
   !Subroutine to define the quadrature points and 
   !weights according to the Gauss-Legendre quadrature
   !====================================================================
   !(Adapted from the gauleg subroutine in Numerical Recipes)
   !x1 and x2 are the upper and lower limits of the quadrature
   subroutine quad_gauss_legendre(x1,x2,npts,xquad,wquad)
   
      use constants, only: pi
      implicit none
      integer,intent(in) :: npts
      real(dp),intent(in) :: x1,x2
      real(dp),intent(out) :: xquad(npts),wquad(npts)
      real(dp), parameter :: eps=3.0e-14_dp
      integer :: i,j,m
      real(dp) :: xl,xm
      real(dp) :: p1,p2,p3,pp,z,z1
      
      m=(npts+1)/2                                                      !The roots are symmetrical, so only 
                                                                        !half of them needs to be found
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do i=1,m                                                          !Loop over the roots
         z = cos(pi*(real(i,dp)-0.25_dp)/(real(npts,dp)+0.5_dp))        !Approximation for the ith root
         do
            p1 = 1._dp
            p2 = 0._dp
            do j=1,npts                                                 !Loop up the recurrence relation
               p3 = p2                                                  !to get the Legendre polynomial
               p2 = p1                                                  !evaluated at z
               p1=(real(2*j-1,dp)*z*p2-real(j-1,dp)*p3)/real(j,dp)      !Desired Legendre polynomial 
            enddo
            pp=real(npts,dp)*(z*p1-p2)/(z*z-1._dp)
            z1 = z
            z = z1 - p1/pp                                              !Newton's method
            if (abs(z-z1).lt.eps) exit
         enddo
         xquad(i) = xm-xl*z                                             !Scale the root to the desired interval,
         xquad(npts+1-i) = xm+xl*z                                      !and put in its symmetric counterpart.
         wquad(i) = 2._dp*xl/((1._dp-z*z)*pp*pp)                        !Compute the weight
         wquad(npts+1-i) = wquad(i)                                     !and its symmetric counterpar
      
      enddo

   endsubroutine quad_gauss_legendre
   
   !====================================================================
   !Subroutine to define the quadrature points and weights according to 
   !the Gauss-Chebyshev quadrature (see Wang & Modest, 2005, JQSRT)
   !Obs: I have made a modification relative to the JQSRT paper; there,
   !the weights are made to sum up to 2, which I do not think is 
   !correct, so I have divided the weights in that formulation by 2
   !====================================================================
   subroutine quad_gauss_chebyshev(npoints,which_quad,xquad,wquad)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use constants, only: pi
      implicit none
      character(*),intent(in) :: which_quad
      integer,intent(in) :: npoints
      integer :: ii,kk,ll
      real(dp),intent(out) :: xquad(npoints),wquad(npoints)
      real(dp) :: thetak,wsum
            
      !-----------------------------------------------------------------
      !Compute the quadrature points and weights
      !-----------------------------------------------------------------
      if (trim(which_quad).eq.'natural') then
         do kk=1,npoints
            ii = npoints - kk + 1                                       !Flip the resulting arrays
            thetak = real(kk,dp)*pi/real(npoints + 1,dp)
            xquad(ii) = dcos(thetak)
            wsum = 0._dp
            do ll=1,(npoints+1)/2
               wsum = wsum + dsin(real(2*ll-1,dp)*thetak)/&
                  real(2*ll-1,dp)
            enddo
            wquad(ii) = 4._dp*dsin(thetak)*wsum/real(npoints+1,dp)
         enddo
         
      elseif (trim(which_quad).eq.'even-rank') then
         do kk=1,npoints
            ii = npoints - kk + 1                                       !Flip the resulting arrays
            thetak = real(kk,dp)*pi/real(2*npoints + 1,dp)
            xquad(ii) = dcos(thetak)
            wsum = 0._dp
            do ll=1,npoints
               wsum = wsum + dsin(real(2*ll-1,dp)*thetak)/&
                  real(2*ll-1,dp)
            enddo
            wquad(ii) = 8._dp*dsin(thetak)*wsum/real(2*npoints+1,dp)
         enddo     
         wquad = wquad/2._dp                                            !My addition
      
      elseif (trim(which_quad).eq.'odd-rank') then
         do kk=1,npoints
            ii = npoints - kk + 1                                       !Flip the resulting arrays
            thetak = real(kk,dp)*pi/real(2*npoints,dp)
            xquad(ii) = dcos(thetak)
            wsum = 0._dp
            do ll=1,npoints
               wsum = wsum + dsin(real(2*ll-1,dp)*thetak)/&
                  real(2*ll-1,dp)
            enddo
            wquad(ii) = 2._dp*real(2-kron_delta(kk,npoints),dp)*&
                        dsin(thetak)*wsum/npoints
         enddo      
         wquad = wquad/2._dp                                            !My addition
      else
         call shutdown('quad_gauss_chebyshev: quadrature not specified')
      endif  
      
   endsubroutine quad_gauss_chebyshev

   function dlegendre(n,x) result(dleg)
      
      implicit none
      integer, intent(in) :: n
      integer :: i 
      real(dp), intent(in) :: x
      real(dp) :: dleg
      real(dp) :: leg_down1, leg_down2, leg
      real(dp) :: dleg_down1, dleg_down2

      selectcase(n)
         case(0)
            dleg = 0
         case(1)
            dleg = 1
         case default                   
            leg_down1  = x; dleg_down1 = 1
            leg_down2  = 1; dleg_down2 = 0
            do i=2,n
               leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i
               dleg = dleg_down2 + (2*i-1)*leg_down1
               leg_down2 = leg_down1
               leg_down1 = leg
               dleg_down2 = dleg_down1
               dleg_down1 = dleg
            enddo
        endselect
    endfunction dlegendre

   function legendre(n,x) result(leg)
      
      implicit none
      integer, intent(in) :: n
      integer :: i
      real(dp), intent(in) :: x
      real(dp) :: leg,leg_down1,leg_down2
      
      selectcase(n)
         case(0)
            leg  = 1
         case(1)
            leg  = x
         case default
            
            leg_down1  = x; leg_down2  = 1._dp
            do i=2,n
               leg = (2*i-1)*x*leg_down1/i - (i-1)*leg_down2/i
               leg_down2 = leg_down1
               leg_down1 = leg
            enddo
        endselect
    endfunction legendre

   subroutine gauss_legendre_lobatto(x,w,interval)
      
      use constants, only: pi
      implicit none
      integer :: i,j,n
      integer,parameter :: newton_iters = 10000
      real(dp),intent(out) :: x(:), w(:)
      real(dp),intent(in),optional :: interval(2)
      real(dp) :: a,b,leg,dleg,delta
      real(dp),parameter :: tolerance = 0.0001_dp

      n = size(x)-1
      selectcase(n)
         case(1)
            x(1) = -1
            x(2) =  1
            w = 1
         
         case default
            x(1)   = -1._dp
            x(n+1) =  1._dp
            w(1)   =  2._dp/(n*(n+1._dp))
            w(n+1) =  2._dp/(n*(n+1._dp))

            do i = 1, (n+1)/2 - 1
               !Initial guess from an approximate form
               !given by SV Parter (1999)
               x(i+1) = -dcos( (i+0.25_dp)*pi/n  - &
                                3/(8*n*pi*(i+0.25_dp)))
               do j=1,newton_iters
                  leg  = legendre(n+1,x(i+1)) - legendre(n-1,x(i+1))
                  dleg = dlegendre(n+1,x(i+1)) - dlegendre(n-1,x(i+1))
                  delta = -leg/dleg
                  x(i+1) = x(i+1) + delta
                  if ( abs(delta) <= tolerance * abs(x(i+1)) )  exit
               enddo
               x(n-i+1) = -x(i+1)

               leg = legendre(n, x(i+1))
               w(i+1)   = 2._dp/(n*(n+1._dp)*leg**2) 
               w(n-i+1) = w(i+1)
            enddo

            if (mod(n,2).eq.0) then
               x(n/2+1) = 0
               leg = legendre(n, 0.0_dp)
               w(n/2+1)   = 2._dp/(n*(n+1._dp)*leg**2) 
            endif
        endselect
        
        if (present(interval)) then
            a = interval(1); b = interval(2)
            x = 0.5_dp*(b-a)*x+0.5_dp*(b+a)
            x(1)       = interval(1)
            x(size(x)) = interval(2)
            w = 0.5_dp*(b-a)*w
        endif
        
    endsubroutine

   !====================================================================
   !Univariate Gaussian distribution function
   !with mean mu and standard deviation std_dev
   !====================================================================
   real(dp) function univariate_gauss(xx,mu,std_dev)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: twopi
      implicit none
      real(dp),intent(in) :: xx,mu,std_dev
      real(dp) :: exp_term
      
      !-----------------------------------------------------------------
      !Calculate the function
      !-----------------------------------------------------------------
      !Exponential term
      exp_term = exp(-0.5_dp*((xx - mu)/std_dev)**2._dp)
      
      !Gauss distribution
      univariate_gauss = exp_term/(std_dev*sqrt(twopi))
   
   endfunction univariate_gauss
   
   !====================================================================
   !Univariate, truncated Gaussian distribution function with mean mu, 
   !standard deviation sigma, and truncation points a and b (with a<b)
   !(note: mu and sigma refer to the non-clipped Gaussian distribution)
   !====================================================================
   real(dp) function truncated_gauss(xx,mu,sigma,a,b)
   
      use constants, only: sqtwo,twopi
      implicit none
      real(dp),intent(in) :: a,b,mu,sigma,xx
      real(dp) :: alpha,beta,phi,xi,z
   
      xi = (xx - mu)/sigma
      alpha = (a - mu)/sigma
      beta = (b - mu)/sigma
      z = 0.5_dp*(erf_function(beta/sqtwo) - erf_function(alpha/sqtwo))
      phi = dexp(-0.5_dp*xi*xi)/dsqrt(twopi)
      truncated_gauss = phi/(sigma*z)
   
   endfunction truncated_gauss
   
   !====================================================================
   !Function to generate a uniform random deviate between 0.0 and 1.0,
   !adapted from  the "ran" function in Numerical Recipes
   !(call with idum a negative integer to initialize; 
   !thereafter, do not alter idum except to reinitialize)
   !====================================================================
   real(dp) function get_numrep_random(idum)

      use precision_parameters, only: i4b
      implicit none
      integer(i4b), intent(inout) :: idum      
      integer(i4b), parameter :: ia=16807,im=2147483647,iq=127773,ir=2836
      real(dp),save :: am
      integer(i4b),save :: ix=-1,iy=-1,k
      
      if (idum <= 0 .or. iy < 0) then
         am=nearest(1.0,-1.0)/real(im,dp)
         iy=ior(ieor(888889999,abs(idum)),1)
         ix=ieor(777755555,abs(idum))
         idum=abs(idum)+1
      end if
      ix=ieor(ix,ishft(ix,13))
      ix=ieor(ix,ishft(ix,-17))
      ix=ieor(ix,ishft(ix,5))
      k=iy/iq
      iy=ia*(iy-k*iq)-ir*k
      if (iy < 0) iy=iy+im
      get_numrep_random=am*ior(iand(im,ieor(ix,iy)),1)
      
   end function get_numrep_random   
   
   !====================================================================
   !Function to get a random number from a uniform distribution 
   !within the range 0 and 1 including 1 and excluding 0 
   !(added also arbitrary mean and standard deviations)
   !====================================================================
   real(dp) function get_random_uniform(mu,std_dev,idum)
   
      use misc_functions, only: get_random_number
      implicit none
      integer,intent(inout) :: idum
      real(dp),intent(in) :: mu,std_dev
      real(dp) :: r1
      
      !Get random number
      r1 = get_numrep_random(idum)
      get_random_uniform = (1._dp - r1)*std_dev + mu
   
   endfunction get_random_uniform
   
   !====================================================================
   !Function to define a random number from a set that follows a 
   !Gaussian distribution with mean mu and standard deviation std_dev
   !====================================================================
   real(dp) function get_random_gauss(mu,std_dev,idum)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: pi
      implicit none
      integer,intent(inout) :: idum
      real(dp),intent(in) :: mu,std_dev
      real(dp) :: r1,r2,r0
   
      !-----------------------------------------------------------------
      !For null standard deviation, the function reduces to a Dirac 
      !delta function and returns the mean value
      !-----------------------------------------------------------------
      if (std_dev.lt.small) then
         get_random_gauss = mu
         return
      endif
      
      !-----------------------------------------------------------------
      !Pick two random numbers with zero
      !mean and unitary standard deviation
      !-----------------------------------------------------------------
      r1 = get_numrep_random(idum)
      r2 = get_numrep_random(idum)
   
      !-----------------------------------------------------------------
      !Apply the method of Murison, 2000
      !-----------------------------------------------------------------
      r0 = sqrt(-2._dp*log(r1))*cos(2._dp*pi*r2)
      get_random_gauss = mu + std_dev*r0
   
   endfunction get_random_gauss

   !====================================================================
   !Function that generates a random number that follows a gamma 
   !distribution with pdf x^(k-1) exp(-x)/Gamma(k); that is, with shape
   !parameter k and scale parameter equal to unity. We use here the 
   !approach of Marsaglia & Tang, 2000
   !====================================================================
   real(dp) recursive function get_random_gamma(shape_parameter,idum) &
      result (ran_num)

      implicit none
      integer,intent(inout) :: idum
      real(dp),intent(in) :: shape_parameter
      real(dp) :: c,d,u,v,x
      
      if (shape_parameter.lt.1._dp) then
         u = get_numrep_random(idum)
         x = get_random_gamma(1._dp + shape_parameter,idum)
         ran_num = x*(u**(1._dp/shape_parameter))
      else
         do
            v = -1._dp
            d = shape_parameter - 1._dp/3._dp; c = 1._dp/dsqrt(9._dp*d)
            do while (v.le.0._dp)
               x = get_random_gauss(0._dp,1._dp,idum)
               v = (1._dp + c*x)**3._dp
            enddo
            u = get_numrep_random(idum)
            if (u.lt.(1._dp - 0.0331_dp*x**4._dp)) exit
            if (log(u).lt.(0.5_dp*x*x + d*(1._dp - v + log(v)))) exit
         enddo
         ran_num = d*v
      endif
      
   endfunction get_random_gamma
   
   !====================================================================
   !Function to generate a random number that follows a beta 
   !distribution with parameters alpha and beta
   !====================================================================
   real(dp) function get_random_beta(alpha,beta,idum)
      
      implicit none
      integer,intent(inout) :: idum
      real(dp),intent(in) :: alpha,beta
      real(dp) :: ga,gb
      
      ga = get_random_gamma(alpha,idum)
      gb = get_random_gamma(beta,idum)
      get_random_beta = min(max(0._dp,ga/(ga + gb + small)),1._dp)
      
   endfunction get_random_beta

   !====================================================================
   !Function to generate a random number that follows a clipped Gaussian
   !distribution with mean mu, standard deviation sigma, and clipped
   !at values a and b (with a<b). 
   !(note: mu and sigma refer to the non-clipped Gaussian distribution)
   !====================================================================
   real(dp) function get_random_clipped_gauss(mu,sigma,a,b,idum)

      implicit none
      integer,intent(inout) :: idum
      real(dp),intent(in) :: a,b,mu,sigma
      real(dp) :: z
      
      do
         z = get_random_gauss(mu,sigma,idum)
         if ((z.ge.a).and.(z.le.b)) exit
      enddo
      get_random_clipped_gauss = z
!      real(dp) :: rho,u,z
      
      !Method of Robert, 1995 (does not work)
!      do
!         z = a + (b-a)*get_numrep_random(idum)
!         if ((a.le.0._dp).and.(b.ge.0._dp)) then
!            rho = dexp(-z*z/2._dp)
!         elseif (b.lt.0._dp) then
!            rho = dexp((b*b - z*z)/2._dp)
!         elseif (a.gt.0._dp) then
!            rho = dexp((a*a - z*z)/2._dp)
!         endif
!         u = get_numrep_random(idum)
!         if (u.le.rho) exit
!      enddo
!      get_random_clipped_gauss = z

   endfunction get_random_clipped_gauss

   !====================================================================
   !Subroutine to apply the 1D TDMA method
   !====================================================================
   !The routine solves the following discretized equation for x_i:
   !       a_i x_i = b_i x_{i+1} + c_i x_{i-1} + d_i
   !where a_i, b_i, c_i and d_i are coefficients and i = ilow,...,iup
   subroutine tdma_1d(x,a,b,c,d,ilow,iup)
      
      !-----------------------------------------------------------------
      !Declarations of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      integer,intent(in) :: ilow,iup
      integer :: ii
      real(dp),intent(in) :: a(ilow:iup),b(ilow:iup),&
                             c(ilow:iup),d(ilow:iup)
      real(dp),intent(out) :: x(ilow:iup)
      real(dp) :: p(ilow:iup),q(ilow:iup)
      
      !-----------------------------------------------------------------
      !Calculate p and q for the first point
      !-----------------------------------------------------------------
      if (a(ilow).eq.0._dp) call shutdown('tdma_1d: null a_1')
      p(ilow) = b(ilow)/(a(ilow) + small)
      q(ilow) = d(ilow)/(a(ilow) + small)
      
      !-----------------------------------------------------------------
      !Obtain p_i and q_i for the remaining points
      !-----------------------------------------------------------------
      do ii=ilow+1,iup
         if (dabs((a(ii) - c(ii)*p(ii-1))).lt.small) &
            call shutdown('tdma_1d: null p_i and q_i denominators')
         p(ii) = b(ii)/(a(ii) - c(ii)*p(ii-1) + small)
         q(ii) = (d(ii) + c(ii)*q(ii-1))/(a(ii) - c(ii)*p(ii-1) + small)
      enddo
      
      !-----------------------------------------------------------------
      !Find x_i
      !-----------------------------------------------------------------
      x(iup) = q(iup)
      do ii=iup-1,ilow,-1
         x(ii) = p(ii)*x(ii+1) + q(ii)
      enddo
      
   endsubroutine tdma_1d

   !====================================================================
   !Subroutine to apply the 2D TDMA method
   !====================================================================
   !The routine solves the following discretized equation for x_i:
   !  a_{i,j} x_{i,j} = b_{i,j} x_{i+1,j} + c_{i,j} x_{i-1,j} 
   !                  + d_{i,j} x_{i,j+1} + e_{i,j} x_{i,j-1} + f_{i,j}
   !where a_{i,j}, b_{i,j}, c_{i,j}, d_{i,j}, e_{i,j} and f_{i,j} are
   !coefficients and i = ilow,...,iup, j = jlow,...jup
   subroutine tdma_2d(x,a,b,c,d,e,f,ilow,iup,jlow,jup)
   
      !-----------------------------------------------------------------
      !Declarations of variables
      !-----------------------------------------------------------------
      implicit none
      integer,intent(in) :: ilow,iup,jlow,jup
      integer :: i,j
      integer :: im,ip,jm,jp
      real(dp),intent(in) :: a(ilow:iup,jlow:jup),b(ilow:iup,jlow:jup),&
                             c(ilow:iup,jlow:jup),d(ilow:iup,jlow:jup),&
                             e(ilow:iup,jlow:jup),f(ilow:iup,jlow:jup)
      real(dp),intent(out) :: x(ilow:iup,jlow:jup)
      real(dp) :: ff(ilow:iup,jlow:jup),x0(ilow:iup,jlow:jup)
      real(dp) :: diff,diffMax,tol
      
      !-----------------------------------------------------------------
      !Setting up parameters for the iteration
      !-----------------------------------------------------------------
      tol = 1.e-6_dp                                                    !Tolerance for the iterative process
      x = 0._dp                                                         !Initial values for the x vector
      x0 = 0._dp                                                        !Initial values for the x vector at
                                                                        !  the previous iteration
      diffMax = 1000._dp                                                !Initial value of the residual
      
      !-----------------------------------------------------------------
      !Iteration loop
      !-----------------------------------------------------------------
      do while (diffMax.gt.tol)
         do i=ilow,iup
            ip = min(i+1,iup); im = max(i-1,ilow)
            ff(i,:) = f(i,:) + b(i,:)*x(ip,:) + c(i,:)*x(im,:)
            call tdma_1d(x(i,:),a(i,:),d(i,:),e(i,:),ff(i,:),jlow,jup)
         enddo

         do j=jlow,jup
            jp = min(j+1,jup); jm = max(j-1,jlow)
            ff(:,j) = f(:,j) + d(:,j)*x(:,jp) + e(:,j)*x(:,jm)
            call tdma_1d(x(:,j),a(:,j),b(:,j),c(:,j),ff(:,j),ilow,iup)         
         enddo
          
         diffMax = -1._dp
         do i=ilow,iup
            do j=jlow,jup
               diff = abs((x(i,j) - x0(i,j))/(x0(i,j) + small))*100._dp
               if (diff.gt.diffMax) diffMax = diff
               x0(i,j) = x(i,j)
            enddo
         enddo
      enddo
   
   endsubroutine tdma_2d
   
   !====================================================================
   !Subroutine to apply the 3D TDMA method
   !====================================================================
   !The routine solves the following discretized equation for x_i:
   !  a_{i,j,k} x_{i,j,k} = 
   !     b_{i,j,k} x_{i+1,j,k} + c_{i,j,k} x_{i-1,j,k} 
   !   + d_{i,j,k} x_{i,j+1,k} + e_{i,j,k} x_{i,j-1,k} 
   !   + f_{i,j,k} x_{i,j,k+1} + g_{i,j,k} x_{i,j,k-1} + h_{i,j,k}
   !where a_{i,j,k}, b_{i,j,k}, c_{i,j,k}, d_{i,j,k}, e_{i,j,k}, 
   !f_{i,j,k}, g_{i,j,k} and h_{i,j,k} are coefficients and 
   !i = ilow,...,iup, j = jlow,...,jup, k = klow,...,kup
   subroutine tdma_3d(x,a,b,c,d,e,f,g,h,ilow,iup,jlow,jup,klow,kup)
   
      !-----------------------------------------------------------------
      !Declarations of variables
      !-----------------------------------------------------------------
      implicit none
      integer,intent(in) :: ilow,iup,jlow,jup,klow,kup
      integer :: i,j,k
      integer :: jp,jm,kp,km
      real(dp),intent(in) :: a(ilow:iup,jlow:jup,klow:kup),&
         b(ilow:iup,jlow:jup,klow:kup),c(ilow:iup,jlow:jup,klow:kup),&
         d(ilow:iup,jlow:jup,klow:kup),e(ilow:iup,jlow:jup,klow:kup),&
         f(ilow:iup,jlow:jup,klow:kup),g(ilow:iup,jlow:jup,klow:kup),&
         h(ilow:iup,jlow:jup,klow:kup)
      real(dp),intent(out) :: x(ilow:iup,jlow:jup,klow:kup)
      real(dp) :: hh(ilow:iup),x0(ilow:iup,jlow:jup,klow:kup)
      real(dp) :: diff,diffMax,tol
      
      !-----------------------------------------------------------------
      !Setting up parameters for the iteration
      !-----------------------------------------------------------------
      tol = 1.e-6_dp                                                    !Tolerance for the iterative process
      x = 0._dp                                                         !Initial values for the x vector
      x0 = 0._dp                                                        !Initial values for the x vector at
                                                                        !  the previous iteration
      diffMax = 1000._dp                                                !Initial value of the residual
      
      !-----------------------------------------------------------------
      !Iteration loop
      !-----------------------------------------------------------------
      do while (diffMax.gt.tol)
         do j=jlow,jup
            do k=klow,kup
               jp = min(j+1,jup); jm = max(j-1,jlow)
               kp = min(k+1,kup); km = max(k-1,klow)
               do i=ilow,iup
                  hh(i) = h(i,j,k) + d(i,j,k)*x(i,jp,k) + &
                          e(i,j,k)*x(i,jm,k) + f(i,j,k)*x(i,j,kp) + &
                          g(i,j,k)*x(i,j,km)
               enddo
               call tdma_1d(x(:,j,k),a(:,j,k),b(:,j,k),&
                            c(:,j,k),hh,ilow,iup)
            enddo
         enddo
         
         diffMax = -1._dp
         do i=ilow,iup
            do j=jlow,jup
               do k=klow,kup
                  diff = abs((x(i,j,k) - x0(i,j,k))/&
                             (x0(i,j,k) + small))*100._dp
                  if (diff.gt.diffMax) diffMax = diff
                  x0(i,j,k) = x(i,j,k)
               enddo
            enddo
         enddo
      enddo
   
   endsubroutine tdma_3d
   
   subroutine tridag(a,b,c,r,u)

      use comp_functions, only: shutdown
      implicit none
      real(dp), dimension(:), intent(in) :: a,b,c,r
      real(dp), dimension(:), intent(out) :: u
      real(dp), dimension(size(b)) :: gam
      integer :: n,j
      real(dp) :: bet
      n=size(u)
      bet=b(1)
      if (bet == 0.0) call shutdown('tridag: error at code stage 1')

      u(1)=r(1)/bet
      do j=2,n   
         gam(j)=c(j-1)/bet
         bet=b(j)-a(j-1)*gam(j)
         if (bet == 0.0) &
         call shutdown('tridag: error at code stage 2')
         u(j)=(r(j)-a(j-1)*u(j-1))/bet
      end do
      do j=n-1,1,-1
         u(j)=u(j)-gam(j+1)*u(j+1)
      end do
   endsubroutine tridag
   
   !====================================================================
   !Subroutine to compute the second-order derivative
   !of a function using a cubic spline 
   !====================================================================
   !Inputs: yarray, xarray: arrays with data for a function y(x)
   !        npt: number of points in the above arrays
   !        ybc_1, ybc_n: first derivatives of y(x) at the first and
   !                      last points
   !
   !Output: y2array: array with the second-order derivative of y
   !(adapted from Numerical Recipes)   
   subroutine spline(xarray,yarray,npts,ybc_1,ybc_n,y2array)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer,intent(in) :: npts
      integer :: ii
      real(dp),intent(in) :: xarray(1:npts),yarray(1:npts),ybc_1,ybc_n
      real(dp),intent(out) :: y2array(1:npts)
      real(dp) :: a(npts),b(npts),c(npts),d(npts)
      
      !-----------------------------------------------------------------
      !Determine coefficients for the core of points
      !-----------------------------------------------------------------
      do ii=2,npts-1
         a(ii) = 2._dp*(xarray(ii-1) - xarray(ii+1))
         b(ii) = (xarray(ii+1) - xarray(ii))
         c(ii) = (xarray(ii) - xarray(ii-1))
         d(ii) = 6._dp*(yarray(ii) - yarray(ii-1))/(c(ii) + small) - &
                 6._dp*(yarray(ii+1) - yarray(ii))/(b(ii) + small)    
!write(*,*) ii,a(ii),b(ii),c(ii),d(ii)
      enddo
      
      !-----------------------------------------------------------------
      !Determine coefficient at the boundaries
      !-----------------------------------------------------------------
      !Left boundary
      if (ybc_1.gt.0.99e30_dp) then
         a(1) = 1._dp
         b(1) = 0._dp
         c(1) = 0._dp
         d(1) = 0._dp
      else
         a(1) = 1._dp
         b(1) = -0.5_dp
         c(1) = 0._dp
         d(1) = (3._dp/(xarray(2)-xarray(1)))*&
                ((yarray(2)-yarray(1))/(xarray(2)-xarray(1)) - ybc_1)
      endif
      
      !Right boundary
      if (ybc_n.gt.0.99e30_dp) then
         a(npts) = 1._dp
         b(npts) = 0._dp
         c(npts) = 0._dp
         d(npts) = 0._dp
      else
         a(npts) = 1._dp
         b(npts) = 0._dp
         c(npts) = -0.5_dp
         d(npts) = -(3._dp/(xarray(npts)-xarray(npts-1)))*&
                    ((yarray(npts)-yarray(npts-1))/&
                     (xarray(npts)-xarray(npts-1)) - ybc_n)
      endif
      
      !-----------------------------------------------------------------
      !Solve the resulting system of equations
      !-----------------------------------------------------------------
      call tdma_1d(y2array,a,b,c,d,1,npts)
      
   endsubroutine spline
   
   !====================================================================
   !Function to compute a spline interpolation
   !====================================================================
   !Inputs: yarray, xarray: arrays with data for a function y(x)
   !        npt: number of points in the above arrays
   !        y2array: array with the second-order derivative of y
   !                 (output of the spline subroutine)
   !        xval: target value of x
   !(adapted from Numerical Recipes)
   real(dp) function splint(xarray,yarray,y2array,xval,npts)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer,intent(in) :: npts
      integer :: klo,khi
      real(dp),intent(in) :: xarray(npts),yarray(npts),y2array(npts),&
         xval
      real(dp) :: a,b,h
      
      !-----------------------------------------------------------------
      !Find the correct indexes for the interpolation
      !-----------------------------------------------------------------
      klo = max(min(locate(xarray,xval,npts),npts-1),1)                 !Lower index
      khi = klo + 1                                                     !Upper index

      !-----------------------------------------------------------------
      !Evaluate the interpolation
      !-----------------------------------------------------------------
      h = xarray(khi)-xarray(klo)
      a =(xarray(khi)-xval)/h
      b =(xval-xarray(klo))/h
      splint = a*yarray(klo) + b*yarray(khi) + &
               ((a**3-a)*y2array(klo)+(b**3-b)*y2array(khi))*(h**2)/6._dp
      
   endfunction splint
   
   !====================================================================
   !Error function
   !====================================================================
   !This implementation relies on the erf function, which is only
   !available in Fortran 2008 and up
   real(dp) function erf_function(xx)
      
      real(dp), intent(in) :: xx
      erf_function = derf(xx)
   
   endfunction erf_function

   !====================================================================
   !Integer function to return an arithmetic progression (for integers)
   !====================================================================
   !Inputs:
   !  ifirst: fisrt value in the sequence
   !  iinc: increment
   !  nvals: number of values in the sequence
   function arth_i(ifirst,iinc,nvals)
   
      use precision_parameters, only: i4b
      integer(i4b),dimension(nvals) :: arth_i
      integer(i4b),intent(in) :: ifirst,iinc,nvals
      integer(i4b) :: kk
      
      if (nvals.gt.0) arth_i(1) = ifirst
      do kk=2,nvals
         arth_i(kk) = arth_i(kk-1) + iinc
      enddo
      
   endfunction arth_i

   !====================================================================
   !Function to evaluate the exponential integral E_n(x), 
   !cf. Eq. (14.31) of Modest (adapted from Numerical Recipes)
   !====================================================================
   real(dp) function expint(n,x)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,shutdown
      use constants, only: euler
      use precision_parameters, only: big,eps,i4b
      implicit none
      integer, intent(in) :: n
      real(dp), intent(in) :: x
      integer(i4b), parameter :: maxit=1000
      integer(i4b) :: i,nm1
      real(dp) :: a,b,c,d,del,fact,h

      call assert(n.ge.0,'n.ge.0')
      call assert(x.ge.0._dp,'x.ge.0._dp')
      call assert((x.gt.0._dp .or. n.gt.1),'(x.gt.0._dp .or. n.gt.1)')
      
      !-----------------------------------------------------------------
      !Special case for n=0
      !-----------------------------------------------------------------
      if (n.eq.0) then
         expint=exp(-x)/x
         return
      end if
      
      nm1=n-1
      if (x.le.small) then
         !--------------------------------------------------------------
         !Special case for x = 0
         !--------------------------------------------------------------
         expint = 1._dp/real(nm1,dp)
      elseif (x.gt.1._dp) then
         !--------------------------------------------------------------
         !Use Lentz' algorithm
         !--------------------------------------------------------------
         b=x+n
         c=big
         d=1._dp/b
         h=d
         do i=1,maxit
            a=-i*(nm1+i)
            b=b+2._dp
            d=1._dp/(a*d+b)                                             !Denominators cannot be zero
            c=b+a/c
            del=c*d
            h=h*del
            if (abs(del-1._dp).le.eps) exit
         end do
         if (i.ge.maxit) call shutdown('expint function: &
            &continued fraction failed')
         expint=h*exp(-x)
      else
         !--------------------------------------------------------------
         !Evaluate series
         !--------------------------------------------------------------
         if (nm1.ne.0) then                                             !Set first term
            expint=1.0_dp/nm1
         else
            expint=-log(x)-euler
         end if
         fact=1._dp
         do i=1,maxit
            fact=-fact*x/i
            if (i.ne.nm1) then
               del=-fact/(i-nm1)
            else
               del=fact*(-log(x)-euler+sum(1._dp/arth_i(1,1,nm1)))
            end if
            expint=expint+del
            if (abs(del).le.abs(expint)*eps) exit
         end do
         if (i.ge.maxit) call shutdown('expint function: series failed')
      end if
   
   endfunction expint

   !====================================================================
   !Ln(gamma) function. Adapted from Numerical Recipes
   !====================================================================
   real(dp) function gammaln(xx)

      use constants, only: twopi
      use comp_functions, only: assert
      implicit none
      integer :: ii
      real(dp), intent(in) :: xx
      real(dp) :: denum,dsum,sqrtp,temp
      real(dp), dimension(6) :: coef

      call assert(xx>0._dp)
      
      !Some parameters
      sqrtp = dsqrt(twopi)
      coef = (/ 76.18009172947146_dp, -86.50532032941677_dp,&
                24.01409824083091_dp, -1.231739572450155_dp,&
                0.1208650973866179e-2_dp, -0.5395239384953e-5_dp /)
      
      !Compute ln(gamma)
      temp = xx + 5.5_dp; temp = (xx + 0.5_dp)*log(temp) - temp
      dsum = 1.000000000190015_dp                                       !c_0
      denum = xx
      do ii=1,6
         denum = denum + 1._dp
         dsum = dsum + coef(ii)/denum
      enddo
      gammaln = temp + log(sqrtp*dsum)
      
   endfunction gammaln
   
   !====================================================================
   !Function to compute the beta PDF at a value xx. It takes two values
   !as input (in1 and in2), whose meaning depend on the value of the
   !string which_data:
   !  If which_data = 'alpha/beta': in1 = alpha, in2 = beta
   !  If which_data = 'mu/sigma': in1 = mu; in2 = sigma
   !====================================================================
   real(dp) function beta_pdf(xx,in1,in2,which_data)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      character(*),intent(in) :: which_data
      real(dp),intent(in) :: in1,in2,xx
      real(dp) :: alpha,beta,mu,sigma
      real(dp) :: denum,lnga,lngab,lngb,rms2
      
      !-----------------------------------------------------------------
      !Define shape parameters from the input values
      !-----------------------------------------------------------------
      if (trim(which_data).eq.('alpha/beta')) then
         alpha = in1; beta = in2
      elseif (trim(which_data).eq.('mu/sigma')) then
         mu = in1; sigma = in2
         rms2 = (mu/(sigma + small))**2
         alpha = rms2*(1._dp - mu) - mu
         beta = alpha*(1._dp/(mu + small) - 1._dp)
      else
         call shutdown('beta_pdf: type of input &
                        &parameters not specified correctly')
      endif
      
      !-----------------------------------------------------------------
      !Compute gamma function at alpha, beta and alpha + beta
      !-----------------------------------------------------------------
      lnga = dlgama(alpha); lngb = dlgama(beta)
      lngab = dlgama(alpha + beta)
      denum = dexp(lnga + lngb - lngab)

      !-----------------------------------------------------------------
      !Compute the pdf
      !-----------------------------------------------------------------
      beta_pdf = (xx**(alpha - 1._dp))*((1._dp - xx)**(beta - 1._dp))/&
                 (denum + small)
      
   endfunction beta_pdf
   
   !====================================================================
   !Subroutine to determine the mean and standard deviation of a
   !truncated Gaussian distribution obtained from a non-truncated 
   !Gaussian distribution (or vice-versa). IMPORTANT: here it is assumed
   !a "true" truncated distribution; i.e., the distribution is only
   !defined within the bounds [a,b].
   !Input:
   !  mu_in       Mean value of the non-truncated distribution
   !  sigma_in    Standard deviation of the non-truncated distribution
   !  a,b         Lower and upper bounds that define the truncated
   !              Gaussian distribution (a<b)
   !  reverse     If .true., reverse the calculation (i.e., mu_in and
   !              sigma_in become the parameters of the truncated
   !              distribution, and mu_out and sigma_out, the parameters
   !              of the non-truncated distribution)
   !Output:
   !  mu_out      Mean value of the truncated distribution
   !  sigma_out   Standard deviation of the truncated distribution
   !====================================================================   
   subroutine Gauss_trunc_conversion(mu_in,sigma_in,a,b,&
                                     mu_out,sigma_out,reverse)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: sqtwo,twopi
      use precision_parameters, only: big
      implicit none
      real(dp),intent(in) :: a,b,mu_in,sigma_in
      real(dp),intent(out) :: mu_out,sigma_out
      real(dp) :: alpha,beta,phi_alpha,phi_beta,z
      real(dp) :: diff,sigma_0,tol
      logical,optional :: reverse
      logical :: clip_to_unclip,unclip_to_clip
      
      !-----------------------------------------------------------------
      !Set logical flags
      !-----------------------------------------------------------------
      clip_to_unclip = .true.
      if (present(reverse)) clip_to_unclip = .not.reverse
      unclip_to_clip = .not.clip_to_unclip
      
      !-----------------------------------------------------------------
      !Convert from non-truncated to truncated
      !-----------------------------------------------------------------
      if (clip_to_unclip) then
         !Auxiliary terms
         alpha = (a - mu_in)/sigma_in
         beta = (b - mu_in)/sigma_in
         z = 0.5_dp*(erf_function(beta/sqtwo) - &
                     erf_function(alpha/sqtwo))
         
         !phi(alpha) and phi(beta), see Wikipedia page
         phi_alpha = dexp(-0.5_dp*alpha*alpha)/dsqrt(twopi)
         phi_beta = dexp(-0.5_dp*beta*beta)/dsqrt(twopi)
         
         !Resulting mean and standard deviation
         mu_out = mu_in + (phi_alpha - phi_beta)*sigma_in/z
         sigma_out = sigma_in*&
            dsqrt(1._dp + (alpha*phi_alpha - beta*phi_beta)/z + &
                  ((phi_alpha - phi_beta)/z)**2._dp)
      endif

      !-----------------------------------------------------------------
      !Convert from truncated to non-truncated
      !(this requires an iterative process)
      !-----------------------------------------------------------------      
      if (unclip_to_clip) then

         !Parameters for the iterative loop
         tol = 1.e-5_dp                                                 !Tolerance for the iterative process
         diff = 10._dp*tol                                              !Value of the relative difference to ensure the iterative loop starts
         sigma_0 = sigma_in                                             !Initial guess value for sigma

         do while(diff.gt.tol)
            !Auxiliary terms
            alpha = (a - mu_in)/sigma_0
            beta = (b - mu_in)/sigma_0
            z = 0.5_dp*(erf_function(beta/sqtwo) - &
                        erf_function(alpha/sqtwo))
         
            !phi(alpha) and phi(beta), see Wikipedia page
            phi_alpha = dexp(-0.5_dp*alpha*alpha)/dsqrt(twopi)
            phi_beta = dexp(-0.5_dp*beta*beta)/dsqrt(twopi)
         
            !Iterate on sigma
            !(this simple iterative procedure should work fine here)
            sigma_out = sigma_0*&
               dsqrt(1._dp + (alpha*phi_alpha - beta*phi_beta)/z + &
                     ((phi_alpha - phi_beta)/z)**2._dp)
            
            !Checking for convergence
            diff = abs((sigma_out - sigma_0)/sigma_0)*100._dp
            
            !Prepare for next iteration
            sigma_0 = sigma_out
         enddo
         
         !Compute mu from the converged sigma
         mu_out = mu_in - (phi_alpha - phi_beta)*sigma_out/z
      endif          
      
   endsubroutine Gauss_trunc_conversion

   !====================================================================
   !Subroutine to sort an array arr_in into the ascending order arry
   !arr_out using straight insertion (adapted from Numerical Recipes)
   !SLOW, BUT WORKS
   !====================================================================
   subroutine sort_si(arr_in,arr_out)

      use precision_parameters, only: dp
      implicit none
      real(dp), dimension(:), intent(in) :: arr_in
      real(dp), dimension(:), intent(out) :: arr_out
      integer :: i,j,n
      real(dp) :: a
      
      arr_out = arr_in
      n=size(arr_in)
      do j=2,n
         a=arr_out(j)
         do i=j-1,1,-1
            if (arr_out(i) <= a) exit
            arr_out(i+1)=arr_out(i)
         enddo
         arr_out(i+1)=a
      enddo
   endsubroutine sort_si

   !====================================================================
   !Subroutine to sort an array arr_in into the ascending order arry
   !arr_out using the Quicksort algorithm
   !(Adapted from Numerical Recipes)
   !For some reason, it is not working right
   !====================================================================
   subroutine quicksort(arr_in,arr_out)
      
      use comp_functions, only: shutdown
      use misc_functions, only: dpswap
      implicit none
      integer, parameter :: nn=15                                       !Size of the subarrays sorted by straight insertion
      integer, parameter :: nstack=50                                   !Required auxiliary storage
      integer :: n,k,i,j,jstack,l,r
      integer,dimension(nstack) :: istack
      real(dp),dimension(:),intent(in) :: arr_in
      real(dp),dimension(:),intent(out) :: arr_out
      real(dp) :: a

      arr_out = arr_in
      n=size(arr_in)
      jstack=0
      l=1
      r=n
      do
         if (r-l < nn) then                                             !Insertion sort when subarray small enough
            do j=l+1,r
               a=arr_out(j)
               do i=j-1,l,-1
                  if (arr_out(i) <= a) exit
                  arr_out(i+1)=arr_out(i)
               enddo
               arr_out(i+1)=a
            enddo
            if (jstack == 0) return
            r=istack(jstack)                                            !Pop stack and begin a new
            l=istack(jstack-1)                                          !   round of partitioning
            jstack=jstack-2
         else                                                           !Choose median of left, center, and right elements
            k=(l+r)/2                                                   !  as partitioning element a. Also rearrange so
            call dpswap(arr_out(k),arr_out(l+1))                        !  that a(l).le.a(l+1).le.a(r).
            call dpswap(arr_out(l),arr_out(r),arr_out(l)>arr_out(r))
            call dpswap(arr_out(l+1),arr_out(r),arr_out(l+1)>arr_out(r))
            call dpswap(arr_out(l),arr_out(l+1),arr_out(l)>arr_out(l+1))
            i=l+1                                                       !Initialize pointers for partitioning
            j=r
            a=arr_out(l+1)                                              !Partitioning element
            do                                                          !Here is the meat
               do                                                       !Scan up to ﬁnd element >= a.
                  i=i+1
                  if (arr_out(i) >= a) exit
               enddo
               do                                                       !Scan down to ﬁnd element <= a
                  j=j-1
                  if (arr_out(j) <= a) exit
               enddo
               if (j < i) exit                                          !Pointers crossed. Exit with partitioning complete
               call dpswap(arr_out(i),arr_out(j))                                 !Exchange elements
            enddo
            arr_out(l+1)=arr_out(j)                                             !Insert partitioning element
            arr_out(j)=a
            jstack=jstack+2                                             !Push pointers to larger subarray on stack; 
                                                                        !  process smaller subarray immediately
            if (jstack > nstack) &
               call shutdown('quicksort: nstack too small')
            if (r-i+1 >= j-l) then
               istack(jstack)=r
               istack(jstack-1)=i
               r=j-1
            else
               istack(jstack)=j-1
               istack(jstack-1)=l
               l=i
            endif
         endif
      enddo
   endsubroutine quicksort

   !====================================================================
   !Routine that generates a distribution (into out_array) within the
   !bounds of min_val and max_val, with a total of npts points, 
   !following the distribution type dist_type
   !====================================================================
   subroutine generate_distribution(dist_type,min_val,max_val,npts,&
      out_array,beta)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      implicit none
      character(*),intent(in) :: dist_type
      integer,intent(in) :: npts
      integer :: ii
      real(dp),intent(in) :: max_val,min_val
      real(dp),intent(in),optional :: beta
      real(dp),intent(out) :: out_array(npts)
      real(dp) :: bbeta,dx,invbbeta,dummy_arr(npts)
      
      bbeta = 0.1_dp; if (present(beta)) bbeta = beta
      if (npts.eq.1) then
         out_array = (max_val + min_val)/2._dp
         return
      endif
      
      selectcase(trim(dist_type))
      case('lin')
         dx = (max_val - min_val)/real(npts-1,dp)
         do ii=1,npts
            out_array(ii) = min_val + real(ii-1,dp)*dx
         enddo

      case('lin-offset')
         dx = (max_val - min_val)/real(npts,dp)
         do ii=1,npts
            out_array(ii) = min_val + 0.5_dp*real(ii,dp)*dx
         enddo

      case('wang')
         invbbeta = 1._dp/bbeta
         dx = (max_val**bbeta - min_val**bbeta)/real(npts-1,dp)
         do ii=1,npts
            out_array(ii) = (min_val**bbeta + &
               real(ii-1,dp)*dx)**invbbeta
         enddo
         out_array(1) = min_val                                         !This is to make sure that the distribution 
         out_array(npts) = max_val                                      !   includes the maximum and minimum values specified
         
      case('webb')                                                      !Standard ALBDF discretization
         do ii=1,npts
            out_array(ii) = (min_val+small)*(max_val/(min_val+small))**&
                           (real(ii-1,dp)/real(npts-1,dp))
         enddo

      case('gauss-legendre')
         call quad_gauss_legendre(0._dp,1._dp,npts,out_array,dummy_arr)
      
      case('gauss-chebyshev-even')
         call quad_gauss_chebyshev(npts,'even-rank',out_array,dummy_arr)
         
      case('gauss-chebyshev-odd')
         call quad_gauss_chebyshev(npts,'odd-rank',out_array,dummy_arr)

      case default
         call shutdown('generate_distribution: dist_type unspecified')

      endselect

   endsubroutine generate_distribution

endmodule math_functions




!   !====================================================================
!   !Function to get the total number of 
!   !emissivity bands between two walls
!   !====================================================================
!   integer function get_n_epsilon()

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      implicit none
!      integer :: i,i1,i2,n_lmbd,n_lmbd_right,n_lmbd_left
      
!      !-----------------------------------------------------------------
!      !Finding input array sizes
!      !-----------------------------------------------------------------
!      n_lmbd_right   = number_epsilon_left_intervals - 1
!      n_lmbd_left    = number_epsilon_right_intervals - 1
   
!      !-----------------------------------------------------------------
!      !Computing output array sizes
!      !-----------------------------------------------------------------
!      i1 = 1
!      i2 = 2
!      i  = 1
!      do
!         if ((i1.ge.n_lmbd_right).and.(i2.ge.n_lmbd_left)) then
!            exit
!         elseif (lambda_right_array(i1).eq.lambda_left_array(i2)) then
!            i1 = i1 + 1
!            i2 = i2 + 1
!            i  = i  + 1
!         elseif (lambda_right_array(i1).lt.lambda_left_array(i2)) then
!            if (i1.ge.n_lmbd_right) then
!               i2 = i2 + 1
!               i  = i  + 1
!            else
!               i1 = i1 + 1
!               i  = i  + 1
!            endif
!         elseif (lambda_left_array(i2).lt.lambda_right_array(i1)) then
!            if (i2.ge.n_lmbd_left) then
!               i1 = i1 + 1
!               i  = i  + 1
!            else
!               i2 = i2 + 1
!               i  = i  + 1
!            endif      
!         endif
!      enddo
!      n_lmbd = i                                                        !(Initial) size of the output lambda array   
!      get_n_epsilon = n_lmbd + 1                                        !Size of the output epsilon arrays
      
!   endfunction get_n_epsilon

!   !====================================================================
!   !Subroutine to order the emissivity 
!   !function for both left and right walls
!   !====================================================================
!   subroutine order_epsilon2(eps_right_out,eps_left_out,lmbd_out,n_eps)

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use precision_parameters, only: dp
!      implicit none
!      integer,intent(in) :: n_eps
!      integer :: i,i1,i2,n_eps_right,n_eps_left,n_lmbd_right,n_lmbd_left
!      real(dp),intent(out) :: eps_right_out(n_eps),eps_left_out(n_eps),&
!         lmbd_out(n_eps)

!      !-----------------------------------------------------------------
!      !Finding input array sizes
!      !-----------------------------------------------------------------  
!      !Emissivities
!      n_eps_right = number_epsilon_left_intervals
!      n_eps_left  = number_epsilon_right_intervals
   
!      !Wavenumbers
!      n_lmbd_right   = n_eps_right - 1
!      n_lmbd_left    = n_eps_left - 1

!      !-----------------------------------------------------------------
!      !Mounting the output wavenumbers array
!      !-----------------------------------------------------------------
!      i1 = 1
!      i2 = 2
!      i  = 1
!      do
!         if ((i1.ge.n_lmbd_right).and.(i2.ge.n_lmbd_left)) then
!            if (n_lmbd_right.ge.n_lmbd_left) then
!               lmbd_out(i) = lambda_right_array(n_lmbd_right)
!            else
!               lmbd_out(i) = lambda_left_array(n_lmbd_left)
!            endif
!            exit
!         elseif (lambda_right_array(i1).eq.lambda_left_array(i2)) then
!            lmbd_out(i) = lambda_right_array(i1)
!            i1 = i1 + 1
!            i2 = i2 + 1
!            i  = i  + 1
!         elseif (lambda_right_array(i1).lt.lambda_left_array(i2)) then
!            if (i1.ge.n_lmbd_right) then
!               lmbd_out(i) = lambda_left_array(i2)
!               i2 = i2 + 1
!               i  = i  + 1
!            else
!               lmbd_out(i) = lambda_right_array(i1)
!               i1 = i1 + 1
!               i  = i  + 1
!            endif
!         elseif (lambda_left_array(i2).lt.lambda_right_array(i1)) then
!            if (i2.ge.n_lmbd_left) then
!               lmbd_out(i) = lambda_right_array(i1)
!               i1 = i1 + 1
!               i  = i  + 1
!            else
!               lmbd_out(i) = lambda_left_array(i2)
!               i2 = i2 + 1
!               i  = i  + 1
!            endif      
!         endif
!      enddo
!      !For convenience, set the last position of the
!      !spectral variable array as zero (its true value 
!      !could also be infinity for wavelengths)
!      lmbd_out(n_eps) = 0._dp
   
!      !-----------------------------------------------------------------
!      !Mounting the output emissivities array
!      !-----------------------------------------------------------------
!      do i=1,n_eps-1
!         i1_loop: do i1=1,n_eps_right-1
!            if (lmbd_out(i).le.lambda_right_array(i1)) then
!               eps_right_out(i) = epsilon_right_array(i1)
!               exit i1_loop
!            else
!               eps_right_out(i) = epsilon_right_array(n_eps_right)
!            endif
!         enddo i1_loop
!         i2_loop: do i2=1,n_eps_left-1
!            if (lmbd_out(i).le.lambda_left_array(i2)) then
!               eps_left_out(i) = epsilon_left_array(i2)
!               exit i2_loop
!            else
!               eps_left_out(i) = epsilon_left_array(n_eps_left)
!            endif
!         enddo i2_loop
!      enddo
!      eps_right_out(n_eps) = epsilon_right_array(n_eps_right)
!      eps_left_out(n_eps)  = epsilon_left_array(n_eps_left)
   
!   endsubroutine order_epsilon2

!   !====================================================================
!   !Subroutine to order the emissivity 
!   !function for both left and right walls
!   !====================================================================
!   subroutine order_epsilon(eps_right_out,eps_left_out,lmbd_out,n_eps)

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use precision_parameters, only: dp
!      implicit none
!      integer,intent(out) :: n_eps
!      integer :: i,i1,i2,n_eps_right,n_eps_left,n_lmbd,&
!         n_lmbd_right,n_lmbd_left
!      real(dp) :: eps_right(100),eps_left(100),lmbd_right(99),&
!         lmbd_left(99)
!      real(dp),allocatable,dimension(:) :: eps_right_out,eps_left_out,&
!         lmbd_out

!      !-----------------------------------------------------------------
!      !Setting the arrays
!      !-----------------------------------------------------------------
!      eps_right = epsilon_right_array
!      eps_left = epsilon_left_array

!      lmbd_right = lambda_right_array
!      lmbd_left = lambda_left_array

!      !-----------------------------------------------------------------
!      !Finding input array sizes
!      !-----------------------------------------------------------------
!      !Emissivities
!      n_eps_right = number_epsilon_left_intervals
!      n_eps_left  = number_epsilon_right_intervals
   
!      !Wavenumbers
!      n_lmbd_right = n_eps_right - 1
!      n_lmbd_left  = n_eps_left - 1
   
!      !-----------------------------------------------------------------
!      !Computing output array sizes
!      !-----------------------------------------------------------------
!      i1 = 1
!      i2 = 2
!      i  = 1
!      do
!         if ((i1.ge.n_lmbd_right).and.(i2.ge.n_lmbd_left)) then
!            exit
!         elseif (lmbd_right(i1).eq.lmbd_left(i2)) then
!            i1 = i1 + 1
!            i2 = i2 + 1
!            i  = i  + 1
!         elseif (lmbd_right(i1).lt.lmbd_left(i2)) then
!            if (i1.ge.n_lmbd_right) then
!               i2 = i2 + 1
!               i  = i  + 1
!            else
!               i1 = i1 + 1
!               i  = i  + 1
!            endif
!         elseif (lmbd_left(i2).lt.lmbd_right(i1)) then
!            if (i2.ge.n_lmbd_left) then
!               i1 = i1 + 1
!               i  = i  + 1
!            else
!               i2 = i2 + 1
!               i  = i  + 1
!            endif      
!         endif
!      enddo
!      n_lmbd = i                                                        !(Initial) size of the output lambda array   
!      n_eps = n_lmbd + 1                                                !Size of the output epsilon arrays
   
!      !-----------------------------------------------------------------
!      !Allocating arrays
!      !-----------------------------------------------------------------
!      allocate(eps_right_out  (1:n_eps))
!      allocate(eps_left_out   (1:n_eps))
!      allocate(lmbd_out       (1:n_eps))
   
!      !-----------------------------------------------------------------
!      !Mounting the output wavenumbers array
!      !-----------------------------------------------------------------
!      i1 = 1
!      i2 = 2
!      i  = 1
!      do
!         if ((i1.ge.n_lmbd_right).and.(i2.ge.n_lmbd_left)) then
!            if (n_lmbd_right.ge.n_lmbd_left) then
!               lmbd_out(i) = lmbd_right(n_lmbd_right)
!            else
!               lmbd_out(i) = lmbd_left(n_lmbd_left)
!            endif
!            exit
!         elseif (lmbd_right(i1).eq.lmbd_left(i2)) then
!            lmbd_out(i) = lmbd_right(i1)
!            i1 = i1 + 1
!            i2 = i2 + 1
!            i  = i  + 1
!         elseif (lmbd_right(i1).lt.lmbd_left(i2)) then
!            if (i1.ge.n_lmbd_right) then
!               lmbd_out(i) = lmbd_left(i2)
!               i2 = i2 + 1
!               i  = i  + 1
!            else
!               lmbd_out(i) = lmbd_right(i1)
!               i1 = i1 + 1
!               i  = i  + 1
!            endif
!         elseif (lmbd_left(i2).lt.lmbd_right(i1)) then
!            if (i2.ge.n_lmbd_left) then
!               lmbd_out(i) = lmbd_right(i1)
!               i1 = i1 + 1
!               i  = i  + 1
!            else
!               lmbd_out(i) = lmbd_left(i2)
!               i2 = i2 + 1
!               i  = i  + 1
!            endif      
!         endif
!      enddo
!      !For convenience, set the last position of the
!      !spectral variable array as zero (its true value 
!      !could also be infinity for wavelengths)
!      lmbd_out(n_eps) = 0._dp
   
!      !-----------------------------------------------------------------
!      !Mounting the output emissivities array
!      !-----------------------------------------------------------------
!      do i=1,n_eps-1
!         i1_loop: do i1=1,n_eps_right-1
!            if (lmbd_out(i).le.lmbd_right(i1)) then
!               eps_right_out(i) = eps_right(i1)
!               exit i1_loop
!            else
!               eps_right_out(i) = eps_right(n_eps_right)
!            endif
!         enddo i1_loop
!         i2_loop: do i2=1,n_eps_left-1
!            if (lmbd_out(i).le.lmbd_left(i2)) then
!               eps_left_out(i) = eps_left(i2)
!               exit i2_loop
!            else
!               eps_left_out(i) = eps_left(n_eps_left)
!            endif
!         enddo i2_loop
!      enddo
!      eps_right_out(n_eps) = eps_right(n_eps_right)
!      eps_left_out(n_eps) = eps_left(n_eps_left)
   
!   endsubroutine order_epsilon
