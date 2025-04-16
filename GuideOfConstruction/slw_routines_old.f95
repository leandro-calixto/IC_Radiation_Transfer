!#######################################################################
!Routines to solve the radiative heat transfer 
!via the SLW model for different scenarios
!#######################################################################
module slw_routines

   !===============
   !Modules & Misc
   !===============
   use precision_parameters
   use constants, only: fourpi,pi,sigma
   use global_parameters
   use mesh
   use slw_parameters, only: slw_cmin,slw_Fmin
   use slw_functions
   use fvm_parameters
   use fvm_routines, only: fvm_solution
   implicit none
   
contains
   
   !=================================================================
   !Subroutine for the radiative heat transfer solution via the
   !SLW- model for a non-scattering medium bounded by black surfaces
   !=================================================================
   subroutine slw_black_solution(n_gases,flux,source,proc_time,&
                                 total_time,absorption,emission,&
                                 kappaPlanck,albdf_ready)
   
      !-------------------------
      !Declaration of variables
      !-------------------------
      use comp_functions, only: CheckMemAlloc,get_file_unit,shutdown
      use slw_parameters, only: slw_bound_Ib,slw_Ib_lbound,slw_Ib_ubound
      use physical_functions, only: bb_emission_frac
      integer,intent(in) :: n_gases
      integer :: nx,ny,nz,nl
      integer :: nch4,nco,nco2,nh2o,nsoot
      integer :: i,j,jgas,k
      integer :: ierr,uout
      real(dp),intent(out) :: proc_time,total_time,&
         flux(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
         source(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp),intent(out),optional :: &
         absorption(0:xcells+1,0:ycells+1,0:zcells+1),&
         emission(0:xcells+1,0:ycells+1,0:zcells+1),&
         kappaPlanck(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp) :: flux_j(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
                  source_j(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp) :: a,cj_loc,F
      real(dp) :: p_ref,T_ref,xh2o_ref,xco2_ref,xsoot_ref,xco_ref,&
                  xch4_ref
      real(dp) :: xmix
      real(dp) :: kp_loc,kp_ref
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref,&
                                           Fj_ref,Fj_sup_ref
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: aIb,Ib,kappa,kappa_sct,&
                                             kIb,u_loc,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: cj_sup_loc,Fj_sup_loc,&
                                                 a_j,kappa_j
      logical,optional :: albdf_ready
      logical :: albdf_read,transparent_window
      logical :: compute_absorption=.false.,compute_emission=.false.,&
                 compute_planck=.false.
   
      !Starting the total computing time counter
      call cpu_time(start_total_time)
      
      !Zeroing out processing time counter 
      !(necessary for the posterior sum over 
      !all instances of the RTE solution)
      proc_time = 0._dp
      
      !Set flags for optional output arguments
      if (present(absorption))   compute_absorption = .true.
      if (present(emission))     compute_emission = .true.
      if (present(kappaPlanck))  compute_planck = .true.
      
      !---------------------------
      !Set up optional parameters
      !---------------------------
      if (present(albdf_ready)) then
         albdf_read = albdf_ready
      else
         albdf_read = .false.                                           !By default, assume that the ALBDF data
      endif                                                             !  has not been read yet
      
      !----------------
      !Surrogate names
      !----------------
      nx = xcells                                                       !Number of cells along direction x
      ny = ycells                                                       !Number of cells along direction y
      nz = zcells                                                       !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization

      nch4  = albdf_ch4_nx                                              !Number of ALBDF mole fraction values for CH4
      nco   = albdf_co_nx                                               !Number of ALBDF mole fraction values for CO
      nco2  = albdf_co2_nx                                              !Number of ALBDF mole fraction values for CO2
      nh2o  = albdf_h2o_nx                                              !Number of ALBDF mole fraction values for H2O
      nsoot = albdf_soot_nx                                             !Number of ALBDF mole fraction values for soot
      
      !---------------------------------------------
      !Prepare for dumping properties, if requested
      !---------------------------------------------
      if (slw_print_properties) then
         !Allocate extra arrays
         allocate(a_j(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
         call CheckMemAlloc('a_j',ierr)
         allocate(kappa_j(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
         call CheckMemAlloc('kappa_j',ierr)
         
         !Prepare output unit
         if (trim(slw_print_file).eq.'null') &                          !Check if a file name was given                          
            call shutdown('slw_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(slw_print_file))                      !Open file
         write(uout,'(3(a,:,","),100(i3,:,","))') ' ',' ','j =',&       !Write header
                                             (jgas,jgas,jgas=0,n_gases)
         write(uout,'(100(a,:,","))') '#x','y','z',&                    !Write subheader
                                    ('kappa_j','a_j',jgas=0,n_gases)
         write(uout,'(100(a,:,","))') '[m]','[m]','[m]',&               !Write subsubheader
                                    ('[1/m]','[-]',jgas=0,n_gases)
      endif
      
      !------------------
      !Allocating arrays
      !------------------
      allocate(aIb(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(cj_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('cj_ref',ierr)
      allocate(cj_sup_loc(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
      call CheckMemAlloc('cj_sup_loc',ierr)
      allocate(cj_sup_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('cj_sup_ref',ierr)
      allocate(Fj_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('Fj_ref',ierr)
      allocate(Fj_sup_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('Fj_sup_ref',ierr)
      allocate(Fj_sup_loc(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
      call CheckMemAlloc('Fj_sup_loc',ierr)
      allocate(Ib(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(phi(1:nl,1:nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(u_loc(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('u_loc',ierr)
      allocate(x_eps(1:2,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(1:2,0:nx+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(1:2,0:nx+1,0:ny+1),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !--------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !--------------------------------------------------------
      !Blackbody intensity
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1
               if (slw_bound_Ib) then
                  F = bb_emission_frac(slw_Ib_lbound,T(i,j,k),.true.) -&
                      bb_emission_frac(slw_Ib_ubound,T(i,j,k),.true.)
               else
                  F = 1._dp
               endif
               Ib(i,j,k) = F*sigma*(T(i,j,k)**4._dp)/pi
            enddo
         enddo
      enddo

      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !Set wall emissivities to unit for all wall cells
      x_eps = 1._dp
      y_eps = 1._dp
      z_eps = 1._dp
      
      !------------------------------------------------------------------------------
      !Defining absorption cross-sections and supplementar absorption cross-sections
      !------------------------------------------------------------------------------
      if (trim(slw_nonuniform_method).eq.'rank_correlated') then
         call get_slw_Fj(n_gases,Fj_ref(1:n_gases),Fj_sup_ref(1:n_gases))
         Fj_sup_ref(0) = slw_Fmin  

      else
         !For other SLWs, divide the reference Cj
         call get_slw_cj(slw_cmin,slw_cmax,n_gases,cj_ref(1:n_gases),&
            cj_sup_ref(1:n_gases))
         cj_sup_ref(0) = slw_cmin
      endif

      !-----------------------------------
      !If needed, read the ALBDF database
      !-----------------------------------
      if (trim(slw_mixture_method).eq.'albdf_precombined') then
         if (.not.albdf_read) call read_aldbf_mix
      else
         if (.not.albdf_read) call read_aldbf
      endif

      !-----------------------------------------
      !Define the reference state, if necessary
      !-----------------------------------------
      if (trim(slw_nonuniform_method).ne.'uniform') &
         call get_slw_reference_state(T_ref,p_ref,xh2o_ref,xco2_ref,&
                                      xsoot_ref,xco_ref,xch4_ref)

      !----------------------------------------
      !Compute scaling parameter, if necessary
      !----------------------------------------
      if (trim(slw_nonuniform_method).eq.'scaled') then
         do i=0,nx+1
            do j=0,ny+1
               do k=0,nz+1
                  !Planck-mean absorption coefficients at
                  !local and reference conditions
                  kp_loc = get_slw_kp(T(i,j,k),p(i,j,k),xh2o(i,j,k),&
                     xco2(i,j,k),xsoot(i,j,k),xco(i,j,k),xch4(i,j,k),&
                     T(i,j,k),n_gases)
                  kp_ref = get_slw_kp(T_ref,p_ref,xh2o_ref,xco2_ref,&
                     xsoot_ref,xco_ref,xch4_ref,T(i,j,k),n_gases)

                  !Scaling parameter
                  u_loc(i,j,k) = kp_loc/(kp_ref+small)
               enddo
            enddo
         enddo
      endif
      
      !----------------------     
      !Initialize sum arrays
      !----------------------
      flux = 0._dp
      source = 0._dp
      kIb = 0._dp
      
      gas_loop: do jgas=0,n_gases
         
         !Check if the gas is a transparent window
         transparent_window = .false.
         if (jgas.eq.0) transparent_window = .true.
         
         !-------------------------------------
         !Computing local radiative properties
         !-------------------------------------
         do i=0,nx+1
            do j=0,ny+1
               do k=0,nz+1
                  !Mole fraction of the mixture of participating species
                  xmix = 0._dp
                  if ((nch4.gt.0).and.(albdf_ch4_name.ne.'null')) &
                     xmix = xmix + xch4(i,j,k)
                  if ((nco.gt.0).and.(albdf_co_name.ne.'null')) &
                     xmix = xmix + xco(i,j,k)   
                  if ((nco2.gt.0).and.(albdf_co2_name.ne.'null')) &
                     xmix = xmix + xco2(i,j,k)    
                  if ((nh2o.gt.0).and.(albdf_h2o_name.ne.'null')) &
                     xmix = xmix + xh2o(i,j,k) 
                  if ((nsoot.gt.0).and.(albdf_soot_name.ne.'null')) &
                     xmix = xmix + xsoot(i,j,k) 

                  selectcase(trim(slw_nonuniform_method))
                  case('uniform')
                     !+++++++++++++++
                     !Uniform medium
                     !+++++++++++++++
                     !(for uniform media, cj_ref = cj_loc, 
                     !so use cj_ref always)
                     if (transparent_window) then
                        kappa(i,j,k) = 0._dp                            !Absorption coefficient of the
                                                                        !  transparent window
                        a = slw_a_func(T(i,j,k),xh2o(i,j,k),&           !Weighting coefficient of the
                           xco2(i,j,k),xsoot(i,j,k),xco(i,j,k),&        !  transparent window
                           xch4(i,j,k),T(i,j,k),cj_sup_ref(jgas))
                     
                     else
                        kappa(i,j,k) = slw_kappa_func(T(i,j,k),&        !Absorption coefficient of the 
                           p(i,j,k),xmix,cj_ref(jgas))                  !  gray gas
                        a = slw_a_func(T(i,j,k),xh2o(i,j,k),&           !Weighting coefficient of the
                           xco2(i,j,k),xsoot(i,j,k),xco(i,j,k),&        !  gray gas
                           xch4(i,j,k),T(i,j,k),cj_sup_ref(jgas),&
                           cj_sup_ref(jgas-1))
                     endif               
            
                  case('scaled')
                     !+++++++++++
                     !Scaled-SLW
                     !+++++++++++
                     if (transparent_window) then
                        kappa(i,j,k) = 0._dp                            !Absorption coefficient of the
                                                                        !  transparent window
                        a = slw_a_func(T_ref,xh2o_ref,&                 !Weighting coefficient of the
                           xco2_ref,xsoot_ref,xco_ref,&                 !  transparent window
                           xch4_ref,T(i,j,k),cj_sup_ref(jgas))
                     else
                        kappa(i,j,k) = u_loc(i,j,k)*slw_kappa_func(&    !Absorption coefficient of the 
                           T(i,j,k),p(i,j,k),xmix,cj_ref(jgas))         !  gray gas
                        a = slw_a_func(T_ref,xh2o_ref,xco2_ref,&        !Weighting coefficient of the
                           xsoot_ref,xco_ref,xch4_ref,T(i,j,k),&        !  gray gas
                           cj_sup_ref(jgas),cj_sup_ref(jgas-1))                           
                     endif
               
                  case('reference_approach')
                     !+++++++
                     !RA-SLW
                     !+++++++
                     Fj_sup_ref(jgas) = get_albdf(T_ref,xh2o_ref,&      !ALBDF evaluated at the
                        xco2_ref,xsoot_ref,xco_ref,xch4_ref,T_ref,&     !  reference state
                        cj_sup_ref(jgas))
                     cj_sup_loc(i,j,k,jgas) = &                         !Get the local supplementar cross-section
                        inverse_albdf(Fj_sup_ref(jgas),T(i,j,k),&       !  by solving the implicit equation
                        xh2o(i,j,k),xco2(i,j,k),xsoot(i,j,k),&
                        xco(i,j,k),xch4(i,j,k),T_ref)
                     if (transparent_window) then
                        kappa(i,j,k) = 0._dp                            !Absorption coefficient of the
                                                                        !  transparent window
                        a = slw_a_func(T_ref,xh2o_ref,xco2_ref,&        !Weighting coefficient of the
                           xsoot_ref,xco_ref,xch4_ref,T(i,j,k),&        !  transparent window
                           cj_sup_ref(jgas))  
                     else
                        cj_loc = &                                      !Local cross-section determined
                           get_slw_single_cj(cj_sup_loc(i,j,k,jgas-1),& !  from the supplementar ones
                                             cj_sup_loc(i,j,k,jgas))   
                        kappa(i,j,k) = slw_kappa_func(&                 !Absorption coefficient of the 
                           T(i,j,k),p(i,j,k),xmix,cj_loc)               !  gray gas
                        a = slw_a_func(T_ref,xh2o_ref,xco2_ref,&        !Weighting coefficient of the
                           xsoot_ref,xco_ref,xch4_ref,T(i,j,k),&        !  gray gas
                           cj_sup_ref(jgas),cj_sup_ref(jgas-1))                                 
                     endif
               
                  case('rank_correlated')
                     !+++++++
                     !RC-SLW
                     !+++++++
                     cj_sup_loc(i,j,k,jgas) = &                         !Get the local supplementar cross-section
                        inverse_albdf(Fj_sup_ref(jgas),T(i,j,k),&       !  by solving the implicit equation using
                           xh2o(i,j,k),xco2(i,j,k),xsoot(i,j,k),&       !  the previously divided supplementar Fs
                           xco(i,j,k),xch4(i,j,k),T_ref)   
                     if (transparent_window) then
                        kappa(i,j,k) = 0._dp                            !Absorption coefficient of the
                                                                        !  transparent window
                        a = slw_a_func(T(i,j,k),xh2o(i,j,k),&           !Weighting coefficient of the
                           xco2(i,j,k),xsoot(i,j,k),xco(i,j,k),&        !  transparent window
                           xch4(i,j,k),T(i,j,k),cj_sup_loc(i,j,k,jgas)) 
                     else
                        cj_loc = inverse_albdf(Fj_ref(jgas),T(i,j,k),&  !Get the local cross-section by solving 
                           xh2o(i,j,k),xco2(i,j,k),xsoot(i,j,k),&       !  the implicit equation using the
                           xco(i,j,k),xch4(i,j,k),T_ref)                !  the previously defined Fs
                        kappa(i,j,k) = slw_kappa_func(T(i,j,k),&        !Absorption coefficient of the 
                           p(i,j,k),xmix,cj_loc)                        !  gray gas
                        a = slw_a_func(T(i,j,k),xh2o(i,j,k),&           !Weighting coefficient of the
                           xco2(i,j,k),xsoot(i,j,k),xco(i,j,k),&        !  gray gas
                           xch4(i,j,k),T(i,j,k),cj_sup_loc(i,j,k,jgas),&
                           cj_sup_loc(i,j,k,jgas-1))
                     endif
                     
                  case('locally_correlated')
                     !+++++++
                     !LC-SLW
                     !+++++++
                     Fj_sup_loc(i,j,k,jgas) = get_albdf(T_ref,xh2o_ref,&!Local ALBDF evaluated using
                        xco2_ref,xsoot_ref,xco_ref,xch4_ref,T(i,j,k),&  !  reference state
                        cj_sup_ref(jgas))
                     cj_sup_loc(i,j,k,jgas) = &                         !Get the local supplementar cross-section
                           inverse_albdf(Fj_sup_loc(i,j,k,jgas),&       !  by solving the implicit equation
                           T(i,j,k),xh2o(i,j,k),xco2(i,j,k),&           
                           xsoot(i,j,k),xco(i,j,k),xch4(i,j,k),T(i,j,k))
                     if (transparent_window) then
                        kappa(i,j,k) = 0._dp                            !Absorption coefficient of the
                                                                        !  transparent window
                        a = Fj_sup_loc(i,j,k,jgas)                      !Weighting coefficient of the
                                                                        !  transparent window
                     else
                        cj_loc = &                                      !Local cross-section determined
                           get_slw_single_cj(cj_sup_loc(i,j,k,jgas-1),& !  from the supplementar ones
                                             cj_sup_loc(i,j,k,jgas))
                        kappa(i,j,k) = slw_kappa_func(&                 !Absorption coefficient of the 
                           T(i,j,k),p(i,j,k),xmix,cj_loc)               !  gray gas
                        a = Fj_sup_loc(i,j,k,jgas) - &                  !Weighting coefficient of the gray gas
                           Fj_sup_loc(i,j,k,jgas-1)
                     endif
               
                  case default
                     call shutdown('Problem in the specification of &
                        &slw_nonuniform_method')
                  endselect
                  aIb(i,j,k) = a*Ib(i,j,k)                              !Emission term
                  kIb(i,j,k) = kIb(i,j,k) + kappa(i,j,k)*aIb(i,j,k)     !Total RTE emission term
                  
                  !Store properties of each gray gas if requested
                  if (slw_print_properties) then
                     kappa_j(i,j,k,jgas) = kappa(i,j,k)
                     a_j(i,j,k,jgas) = a
                  endif
               enddo
            enddo
         enddo 

         !--------------------------
         !Solve the radiation field
         !--------------------------
         !Starting the processing time counter
         call cpu_time(start_proc_time)
   
         !Solving the radiation field for gas jgas
         call fvm_solution(kappa,kappa,kappa_sct,aIb,phi,x_eps,y_eps,&
                           z_eps,flux_j,source_j,.true.)
                           
         !Updating the heat flux and heat source
         flux = flux + flux_j
         source = source + source_j
      
         !Stopping the processing time counter and adding to the total
         call cpu_time(end_proc_time)
         proc_time = proc_time + end_proc_time - start_proc_time
      
      enddo gas_loop
      
      !--------------------------------------------
      !Manage dump of gas properties, if necessary
      !--------------------------------------------
      if (slw_print_properties) then
         !Dump data
         do i=0,nx+1
            do j=0,ny+1
               do k=0,nz+1
                  write(uout,'(1000(e26.15e3,:,","))')&
                     x(i),y(j),z(k),(kappa_j(i,j,k,jgas),&
                        a_j(i,j,k,jgas),jgas=0,n_gases)
               enddo
            enddo
         enddo
         
         !Close output unit
         if (slw_print_properties) close(uout)
      
         !Deallocate extra arrays
         deallocate(kappa_j)
         deallocate(a_j)
      endif
      
      !-----------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1
               if (compute_emission) &
                  emission(i,j,k) = fourpi*kIb(i,j,k)                   !Total emission
               if (compute_absorption) &
                  absorption(i,j,k) = source(i,j,k) + fourpi*kIb(i,j,k) !Total absorption
               if (compute_planck) &
                  kappaPlanck(i,j,k) = kIb(i,j,k)/Ib(i,j,k)             !Planck-mean absorption coefficient
            enddo
         enddo
      enddo
      
      !------------------
      !Deallocate arrays
      !------------------
      deallocate(aIb)
      deallocate(cj_ref)
      deallocate(cj_sup_loc)
      deallocate(cj_sup_ref)
      deallocate(Fj_ref)
      deallocate(Fj_sup_ref)
      deallocate(Fj_sup_loc)
      deallocate(Ib)
      deallocate(kappa)
      deallocate(kappa_sct)
      deallocate(kIb)
      deallocate(phi)
      deallocate(u_loc)
      deallocate(x_eps)
      deallocate(y_eps)
      deallocate(z_eps)
      
		!---------------------
      !Total computing time
      !---------------------
      call cpu_time(end_total_time)
      total_time = end_total_time - start_total_time
      
   endsubroutine slw_black_solution   

   !==============================================================
   !Subroutine for the radiative heat transfer solution
   !via the multiple integration SLW model for a non-
   !scattering medium bounded by black walls
   !==============================================================
   subroutine slw_multint_black_solution(n_gases,flux,source,proc_time,&
                                         total_time,albdf_ready)
	
		!-------------------------
		!Declaration of variables
		!-------------------------
      use comp_functions, only: CheckMemAlloc,shutdown
      real(dp),intent(out) :: proc_time,total_time,&
         flux(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
         source(0:xcells+1,0:ycells+1,0:zcells+1)
      integer,intent(in) :: n_gases
      integer :: i,j,jg,k,m,n
      integer :: counter_species,ierr
      integer :: ich4,ico,ico2,ih2o,isoot
      integer :: nx,ny,nz,nl,ns,ng_total
      integer :: nch4,nco,nco2,nh2o,nsoot
      integer,allocatable,dimension(:,:) :: jn
      real(dp) :: flux_j(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
                  source_j(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp) :: xxh2o,xxco2,xxsoot,xxco,xxch4
      real(dp) :: p_ref,T_ref,xh2o_ref,xco2_ref,xsoot_ref,xco_ref,&
                  xch4_ref
      real(dp) :: xxh2o_ref,xxco2_ref,xxsoot_ref,xxco_ref,xxch4_ref
      real(dp) :: cj_loc,xmix
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref,&
         Fj_ref,Fj_sup_ref
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: aIb,Ib,kappa,kappa_sct,&
         x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: an,kappan
      real(dp),allocatable,dimension(:,:,:,:,:) :: cj_sup_loc,Fj_sup_loc
      logical,optional :: albdf_ready
      logical :: albdf_read,next_gas,transparent_window
      
      !Starting the total computing time counter
      call cpu_time(start_total_time)
      
      !Zeroing out processing time counter 
      !(necessary for the posterior sum over 
      !all instances of the RTE solution)
      proc_time = 0._dp
      
      !---------------------------
      !Set up optional parameters
      !---------------------------
      if (present(albdf_ready)) then
         albdf_read = albdf_ready
      else
         albdf_read = .false.                                           !By default, assume that the ALBDF data
      endif                                                             !  has not been read yet
      
      !----------------
      !Surrogate names
      !----------------
      nx = xcells                                                       !Number of cells along direction x
      ny = ycells                                                       !Number of cells along direction y
      nz = zcells                                                       !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization

      nch4  = albdf_ch4_nx                                              !Number of ALBDF mole fraction values for CH4
      nco   = albdf_co_nx                                               !Number of ALBDF mole fraction values for CO
      nco2  = albdf_co2_nx                                              !Number of ALBDF mole fraction values for CO2
      nh2o  = albdf_h2o_nx                                              !Number of ALBDF mole fraction values for H2O
      nsoot = albdf_soot_nx                                             !Number of ALBDF mole fraction values for soot
     
      !-------------------------------------
      !Parameters of the individual species
      !-------------------------------------
      !Counting the species
      counter_species = 0
      if ((nh2o.gt.0).and.(albdf_h2o_name.ne.'null')) then
         counter_species = counter_species + 1
         ih2o = counter_species                                         !Set the species identifier
     endif
      if ((nco2.gt.0).and.(albdf_co2_name.ne.'null')) then
         counter_species = counter_species + 1
         ico2 = counter_species                    
      endif
      if ((nsoot.gt.0).and.(albdf_soot_name.ne.'null')) then
         counter_species = counter_species + 1
         isoot = counter_species
      endif
      if ((nco.gt.0).and.(albdf_co_name.ne.'null')) then
         counter_species = counter_species + 1
         ico = counter_species
      endif
      if ((nch4.gt.0).and.(albdf_ch4_name.ne.'null')) then
         counter_species = counter_species + 1
         ich4 = counter_species
      endif
      
      !Total number of species
      ns = counter_species

		!------------------------
		!Defining SLW parameters
		!------------------------
		!Total number of gray gases
      ng_total = (n_gases+1)**ns

      !Allocate arrays
      allocate(jn(1:ng_total,1:ns),stat=ierr)
      call CheckMemAlloc('jn',ierr)
      
      !Mount the jn array 
      !(array storing the gray gas index of each individual
      !species n that corresponds to the overal gray gas jg)
      jg = 1
      jn(1,:) = 0
      do n=1,ns
         do jg=2,ng_total
            if (n.eq.1) then                                            !Particularize n=1
               jn(jg,n) = jn(jg-1,n) + 1
               if (jn(jg,n).gt.n_gases) jn(jg,n) = 0                    !Return to the transparent windows
            else
               next_gas = .true.                                        !This flag indicats if all
               do m=1,n-1                                               !previous species have reached
                  if (jn(jg-1,m).ne.n_gases) next_gas = .false.         !their maximum number of gray gases
               enddo                                                    !for this jg
               if (next_gas) then
                  jn(jg,n) = jn(jg-1,n) + 1                             !Update gas index for this species
                  if (jn(jg,n).gt.n_gases) jn(jg,n) = 0                 !Return to the transparent windows
               else
                  jn(jg,n) = jn(jg-1,n)
               endif
            endif
         enddo
      enddo
      
      !----------------------------
      !Allocating all other arrays
      !----------------------------     
      allocate(an(0:nx+1,0:ny+1,0:nz+1,1:ns),stat=ierr)
      call CheckMemAlloc('an',ierr)
      allocate(aIb(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(cj_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('cj_ref',ierr)
      allocate(cj_sup_loc(0:nx+1,0:ny+1,0:nz+1,0:n_gases,1:ns),stat=ierr)
      call CheckMemAlloc('cj_sup_loc',ierr)
      allocate(cj_sup_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('cj_sup_ref',ierr)
      allocate(Fj_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('Fj_ref',ierr)
      allocate(Fj_sup_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('Fj_sup_ref',ierr)
      allocate(Fj_sup_loc(0:nx+1,0:ny+1,0:nz+1,0:n_gases,1:ns),stat=ierr)
      call CheckMemAlloc('Fj_sup_loc',ierr)
      allocate(Ib(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappan(0:nx+1,0:ny+1,0:nz+1,1:ns),stat=ierr)
      call CheckMemAlloc('kappan',ierr)
      allocate(kappa_sct(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(phi(1:nl,1:nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(x_eps(1:2,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(1:2,0:nx+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(1:2,0:nx+1,0:ny+1),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !------------------------------------------------------------------------------
      !Defining absorption cross-sections and supplementar absorption cross-sections
      !------------------------------------------------------------------------------
      if (trim(slw_nonuniform_method).eq.'rank_correlated') then
         call get_slw_Fj(n_gases,Fj_ref(1:n_gases),Fj_sup_ref(1:n_gases))
         Fj_sup_ref(0) = slw_Fmin  
      else
         !For other SLWs, divide the reference Cj
         call get_slw_cj(slw_cmin,slw_cmax,n_gases,cj_ref(1:n_gases),&
            cj_sup_ref(1:n_gases))
         cj_sup_ref(0) = slw_cmin
      endif

      !-----------------------------------
      !If needed, read the ALBDF database
      !-----------------------------------
      if (.not.albdf_read) call read_aldbf

      !-----------------------------------------
      !Define the reference state, if necessary
      !-----------------------------------------
      if (trim(slw_nonuniform_method).ne.'uniform') &
         call get_slw_reference_state(T_ref,p_ref,xh2o_ref,xco2_ref,&
                                      xsoot_ref,xco_ref,xch4_ref)
                                      
      !----------------------------------
      !Initializing arrays and variables
      !----------------------------------
      flux = 0._dp
      source = 0._dp

      !--------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !--------------------------------------------------------
      !Blackbody intensity
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1
               Ib(i,j,k) = sigma*(T(i,j,k)**4._dp)/pi
            enddo
         enddo
      enddo
      
      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !Set wall emissivities to unit for all wall cells
      x_eps = 1._dp
      y_eps = 1._dp
      z_eps = 1._dp

      all_gases_loop: do jg=1,ng_total
         !-----------------------------------------------------
         !Computing the properties that depend on the gray gas
         !-----------------------------------------------------
         do i=0,nx+1
            do j=0,ny+1
               do k=0,nz+1                   
                  do n=1,ns
                     if (jg.ge.2) then                                  !Avoid computing the properties
                        if (jn(jg,n).eq.(jn(jg-1,n))) cycle             !  of a repeated individual species
                     endif                                              !  gray gas
                     
                     !Check if the gas is a transparent window
                     transparent_window = .false.
                     if (jn(jg,n).eq.0) transparent_window = .true.
         
                     !Define the local mole fraction
                     xxh2o  = 0._dp
                     xxco2  = 0._dp
                     xxco   = 0._dp
                     xxch4  = 0._dp
                     xxsoot = 0._dp
                     if (n.eq.ih2o)  xxh2o  = xh2o(i,j,k)
                     if (n.eq.ico2)  xxco2  = xco2(i,j,k)
                     if (n.eq.ico)   xxco   = xco(i,j,k)
                     if (n.eq.ich4)  xxch4  = xch4(i,j,k)
                     if (n.eq.isoot) xxsoot = xsoot(i,j,k)
                     xmix = max(xxh2o,xxco2,xxco,xxch4,xxsoot)
                     if (n.eq.isoot) xmix = 1._dp                       !Correct for the soot
                                                                        !absorption coefficient calculation
               
                     !Define the reference mole fraction values
                     xxh2o_ref  = 0._dp
                     xxco2_ref  = 0._dp
                     xxco_ref   = 0._dp
                     xxch4_ref  = 0._dp
                     xxsoot_ref = 0._dp
                     if (n.eq.ih2o)  xxh2o_ref  = xh2o_ref
                     if (n.eq.ico2)  xxco2_ref  = xco2_ref
                     if (n.eq.ico)   xxco_ref   = xco_ref
                     if (n.eq.ich4)  xxch4_ref  = xch4_ref
                     if (n.eq.isoot) xxsoot_ref = xsoot_ref
                     
                     selectcase(trim(slw_nonuniform_method))
                     case('uniform')
                        !+++++++++++++++
                        !Uniform medium
                        !+++++++++++++++
                        !(for uniform media, cj_ref = cj_loc, 
                        !so use cj_ref always)
                        if (transparent_window) then
                           kappan(i,j,k,n) = 0._dp                            !Absorption coefficient of the
                                                                        !  transparent window
                           an(i,j,k,n) = slw_a_func(T(i,j,k),xxh2o,xxco2,xxsoot,&!Weighting coefficient of the
                              xxco,xxch4,T(i,j,k),cj_sup_ref(jn(jg,n))) !  transparent window
                        
                        else
                           kappan(i,j,k,n) = slw_kappa_func(T(i,j,k),p(i,j,k),&  !Absorption coefficient of the 
                              xmix,cj_ref(jn(jg,n)))                    !  gray gas
                           an(i,j,k,n) = slw_a_func(T(i,j,k),xxh2o,xxco2,xxsoot,&!Weighting coefficient of the
                              xxco,xxch4,T(i,j,k),cj_sup_ref(jn(jg,n)),&!  gray gas
                              cj_sup_ref(jn(jg,n)-1))
                        endif               
            
                     case('scaled')
                        !+++++++++++
                        !Scaled-SLW
                        !+++++++++++
                        call shutdown('Scaled model not available for now')
                  
                     case('reference_approach')
                        !+++++++
                        !RA-SLW
                        !+++++++
                        if (transparent_window) then
                           kappan(i,j,k,n) = 0._dp                      !Absorption coefficient of the
                                                                        !  transparent window
                           an(i,j,k,n) = slw_a_func(T_ref,xxh2o_ref,xxco2_ref,&  !Weighting coefficient of the
                              xxsoot_ref,xxco_ref,xxch4_ref,T(i,j,k),&  !  transparent window
                              cj_sup_ref(jn(jg,n)))  
                        else
                           Fj_sup_ref(jn(jg,n)) = get_albdf(T_ref,&     !ALBDF evaluated at the
                              xxh2o_ref,xxco2_ref,xxsoot_ref,xxco_ref,& !  reference state
                              xxch4_ref,T_ref,cj_sup_ref(jn(jg,n)))
                           cj_sup_loc(i,j,k,jn(jg,n),n) = &             !Get the local supplementar cross-section
                              inverse_albdf(Fj_sup_ref(jn(jg,n)),&      !  by solving the implicit equation
                              T(i,j,k),xxh2o,xxco2,xxsoot,xxco,xxch4,&
                              T_ref)
                           cj_loc = &                                   !Local cross-section determined
                              get_slw_single_cj(&                       !  from the supplementar ones
                                 cj_sup_loc(i,j,k,jn(jg,n)-1,n),& 
                                 cj_sup_loc(i,j,k,jn(jg,n),n))   
                           kappan(i,j,k,n) = slw_kappa_func(T(i,j,k),p(i,j,k),&  !Absorption coefficient of the 
                              xmix,cj_loc)                              !  gray gas
                           an(i,j,k,n) = slw_a_func(T_ref,xxh2o_ref,xxco2_ref,&  !Weighting coefficient of the
                              xxsoot_ref,xxco_ref,xxch4_ref,T(i,j,k),&  !  gray gas
                              cj_sup_ref(jn(jg,n)),&
                              cj_sup_ref(jn(jg,n)-1))    
                        endif
               
                     case('rank_correlated')
                        !+++++++
                        !RC-SLW
                        !+++++++
                        cj_sup_loc(i,j,k,jn(jg,n),n) = &                !Get the local supplementar cross-section
                           inverse_albdf(Fj_sup_ref(jn(jg,n)),T(i,j,k),&!  by solving the implicit equation using
                              xxh2o,xxco2,xxsoot,xxco,xxch4,T_ref)      !  the previously divided supplementar Fs
                        if (transparent_window) then
                           kappan(i,j,k,n) = 0._dp                               !Absorption coefficient of the
                                                                        !  transparent window
                           an(i,j,k,n) = slw_a_func(T(i,j,k),xxh2o,xxco2,xxsoot,&!Weighting coefficient of the
                              xxco,xxch4,T(i,j,k),&                     !  transparent window
                              cj_sup_loc(i,j,k,jn(jg,n),n))  
                        else
                           cj_loc = inverse_albdf(Fj_ref(jn(jg,n)),&    !Get the local cross-section by solving 
                              T(i,j,k),xxh2o,xxco2,xxsoot,xxco,xxch4,&  !  the implicit equation using the
                              T_ref)                                    !  the previously defined Fs
                           kappan(i,j,k,n) = slw_kappa_func(T(i,j,k),p(i,j,k),&  !Absorption coefficient of the 
                              xmix,cj_loc)                              !  gray gas
                           an(i,j,k,n) = slw_a_func(T(i,j,k),xxh2o,xxco2,xxsoot,&!Weighting coefficient of the
                              xxco,xxch4,T(i,j,k),&                     !  gray gas
                              cj_sup_loc(i,j,k,jn(jg,n),n),&
                              cj_sup_loc(i,j,k,jn(jg,n)-1,n))     
                        endif
         
                     case('locally_correlated')
                        !+++++++
                        !LC-SLW
                        !+++++++
                        Fj_sup_loc(i,j,k,jn(jg,n),n) = get_albdf(T_ref,&!Local ALBDF evaluated using
                           xxh2o_ref,xxco2_ref,xxsoot_ref,xxco_ref,&    !  reference state
                           xxch4_ref,T(i,j,k),cj_sup_ref(jn(jg,n)))
                        if (transparent_window) then
                           kappan(i,j,k,n) = 0._dp                               !Absorption coefficient of the
                                                                        !  transparent window
                           an(i,j,k,n) = Fj_sup_loc(i,j,k,jn(jg,n),n)            !Weighting coefficient of the
                                                                        !  transparent window
                        else
                           cj_sup_loc(i,j,k,jn(jg,n),n) = &             !Get the local supplementar cross-section
                              inverse_albdf(&                           !  by solving the implicit equation
                              Fj_sup_loc(i,j,k,jn(jg,n),n),T(i,j,k),&
                              xxh2o,xxco2,xxsoot,xxco,xxch4,T(i,j,k))
                           cj_loc = &                                   !Local cross-section determined
                              get_slw_single_cj(&                       !  from the supplementar ones
                                 cj_sup_loc(i,j,k,jn(jg,n)-1,n),&
                                 cj_sup_loc(i,j,k,jn(jg,n),n))
                           kappan(i,j,k,n) = slw_kappa_func(&                    !Absorption coefficient of the 
                              T(i,j,k),p(i,j,k),xmix,cj_loc)            !  gray gas
                           an(i,j,k,n) = Fj_sup_loc(i,j,k,jn(jg,n),n) - &        !Weighting coefficient of the gray gas
                              Fj_sup_loc(i,j,k,jn(jg,n)-1,n)
                        endif
                  
                     case default
                        call shutdown('Problem in the specification of &
                           &slw_nonuniform_method')
                     endselect                  
                  enddo        
                            
                  !Computing the global kappa and a
                  kappa(i,j,k) = sum(kappan(i,j,k,:))
                  aIb(i,j,k) = product(an(i,j,k,:))*Ib(i,j,k)
!                  if (i.eq.nx/2) write(*,*) jg,jn(jg,:),&
!                     kappa(i,j,k),product(an(i,j,k,:))
               enddo
            enddo
         enddo
          
         !--------------------------
         !Solve the radiation field
         !--------------------------
         write(*,*) jn(jg,:)
         !Starting the processing time counter
         call cpu_time(start_proc_time)
   
         !Solving the radiation field for gas jgas
         call fvm_solution(kappa,kappa,kappa_sct,aIb,phi,x_eps,y_eps,&
                           z_eps,flux_j,source_j,.true.)
                           
         !Updating the heat flux and heat source
         flux = flux + flux_j
         source = source + source_j

         !Stopping the processing time counter and adding to the total
         call cpu_time(end_proc_time)
         proc_time = proc_time + end_proc_time - start_proc_time
         
      enddo all_gases_loop
      
      !------------------
      !Deallocate arrays
      !------------------
      deallocate(an)
      deallocate(aIb)
      deallocate(cj_ref)
      deallocate(cj_sup_loc)
      deallocate(cj_sup_ref)
      deallocate(Fj_ref)
      deallocate(Fj_sup_ref)
      deallocate(Fj_sup_loc)
      deallocate(Ib)
      deallocate(jn)
      deallocate(kappa)
      deallocate(kappan)
      deallocate(kappa_sct)
      deallocate(phi)
      deallocate(x_eps)
      deallocate(y_eps)
      deallocate(z_eps)
      
		!---------------------
      !Total computing time
      !---------------------
      call cpu_time(end_total_time)
      total_time = end_total_time - start_total_time

   endsubroutine slw_multint_black_solution

   !==================================================================
   !Subroutine for the radiative heat transfer solution via the
   !SLW-1 model for a non-scattering medium bounded by black surfaces
   !==================================================================
   subroutine slw1_black_solution(flux,source,proc_time,total_time,&
                                  absorption,emission,kappaPlanck,&
                                  albdf_ready)
   
      !-------------------------
      !Declaration of variables
      !-------------------------
      use comp_functions, only: CheckFileExists,CheckMemAlloc,&
                                get_file_unit,shutdown
      integer :: ng_ref,nl,nx,ny,nz
      integer :: nch4,nco,nco2,nh2o,nsoot
      integer :: i,j,jgas,k
      integer :: ierr,uout
      real(dp),intent(out) :: proc_time,total_time,&
         flux(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
         source(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp),intent(out),optional :: &
         absorption(0:xcells+1,0:ycells+1,0:zcells+1),&
         emission(0:xcells+1,0:ycells+1,0:zcells+1),&
         kappaPlanck(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp) :: flux_j(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
                  source_j(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp) :: a_ref(0:1),kappa_ref(0:1)
      real(dp) :: C_ref,C_lcl,F_ref
      real(dp) :: p_ref,T_ref,xh2o_ref,xco2_ref,xsoot_ref,xco_ref,&
                  xch4_ref,xmix_ref
      real(dp) :: xmix
      real(dp) :: len1,len2,val1,val2
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: a,aIb,Ib,kappa,&
         kappa_sct,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: a_j,kappa_j
      logical,optional :: albdf_ready
      logical :: albdf_read,transparent_window
      logical :: compute_absorption=.false.,compute_emission=.false.,&
                 compute_planck=.false.
   
      !Starting the total computing time counter
      call cpu_time(start_total_time)
      
      !Zeroing out processing time counter 
      !(necessary for the posterior sum over 
      !all instances of the RTE solution)
      proc_time = 0._dp
      
      !Set flags for optional output arguments
      if (present(absorption))   compute_absorption = .true.
      if (present(emission))     compute_emission = .true.
      if (present(kappaPlanck))  compute_planck = .true.
      
      !---------------------------
      !Set up optional parameters
      !---------------------------
      if (present(albdf_ready)) then
         albdf_read = albdf_ready
      else
         albdf_read = .false.                                           !By default, assume that the ALBDF data
      endif                                                             !  has not been read yet
      
      !----------------
      !Surrogate names
      !----------------
      nx = xcells                                                       !Number of cells along direction x
      ny = ycells                                                       !Number of cells along direction y
      nz = zcells                                                       !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization

      nch4  = albdf_ch4_nx                                              !Number of ALBDF mole fraction values for CH4
      nco   = albdf_co_nx                                               !Number of ALBDF mole fraction values for CO
      nco2  = albdf_co2_nx                                              !Number of ALBDF mole fraction values for CO2
      nh2o  = albdf_h2o_nx                                              !Number of ALBDF mole fraction values for H2O
      nsoot = albdf_soot_nx                                             !Number of ALBDF mole fraction values for soot
      
      len1 = slw1_length1                                               !First characteristic length
      len2 = slw1_length2                                               !Second characteristic length
      ng_ref = slw1_ngases                                              !Number of gray gases for the reference SLW solution

      !---------------------------------------------
      !Prepare for dumping properties, if requested
      !---------------------------------------------
      if (slw_print_properties) then
         !Allocate extra arrays
         allocate(a_j(0:nx+1,0:ny+1,0:nz+1,0:1),stat=ierr)
         call CheckMemAlloc('a_j',ierr)
         allocate(kappa_j(0:nx+1,0:ny+1,0:nz+1,0:1),stat=ierr)
         call CheckMemAlloc('kappa_j',ierr)
         
         !Prepare output unit
         if (trim(slw_print_file).eq.'null') &                          !Check if a file name was given                          
            call shutdown('slw_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(slw_print_file))                      !Open file
         write(uout,'(100(a,:,","))') '#x','y','z','kappa_0','a_0',&    !Write header
                                      'kappa_1','a_1'
         write(uout,'(100(a,:,","))') '[m]','[m]','[m]','[1/m]','[-]',& !Write subheader
                                      '[1/m]','[-]'
      endif
      
      !------------------
      !Allocating arrays
      !------------------
      allocate(a(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('a',ierr)
      allocate(aIb(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(Ib(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(phi(1:nl,1:nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(x_eps(1:2,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(1:2,0:nx+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(1:2,0:nx+1,0:ny+1),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)
      
      !----------------------------------
      !Initializing arrays and variables
      !----------------------------------
      flux = 0._dp
      source = 0._dp

      !--------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !--------------------------------------------------------
      !Blackbody intensity
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1
               Ib(i,j,k) = sigma*(T(i,j,k)**4._dp)/pi
            enddo
         enddo
      enddo
      
      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !Set wall emissivities to unit for all wall cells
      x_eps = 1._dp
      y_eps = 1._dp
      z_eps = 1._dp
      
      !-----------------------------------
      !If needed, read the ALBDF database
      !-----------------------------------
      if (trim(slw_mixture_method).eq.'albdf_precombined') then
         if (.not.albdf_read) call read_aldbf_mix
      else
         if (.not.albdf_read) call read_aldbf
      endif

      !---------------------------
      !Define the reference state
      !---------------------------
      call get_slw_reference_state(T_ref,p_ref,xh2o_ref,xco2_ref,&
                                   xsoot_ref,xco_ref,xch4_ref)
      xmix_ref = xh2o_ref + xco2_ref + xsoot_ref + xco_ref + xch4_ref   !Reference mixture mole fraction

      !----------------------------------
      !Compute the two target quantities
      !----------------------------------
      selectcase(trim(slw1_method))
      case('kappa-epsilon')
         val1 = get_slw_kp(T_ref,p_ref,xh2o_ref,xco2_ref,xsoot_ref,&
                           xco_ref,xch4_ref,T_ref,ng_ref)
         val2 = get_slw_emissivity(T_ref,p_ref,xh2o_ref,xco2_ref,&
                           xsoot_ref,xco_ref,xch4_ref,T_ref,ng_ref,len1)
      case('epsilon-epsilon')
         val1 = get_slw_emissivity(T_ref,p_ref,xh2o_ref,xco2_ref,&
                           xsoot_ref,xco_ref,xch4_ref,T_ref,ng_ref,len1)
         val2 = get_slw_emissivity(T_ref,p_ref,xh2o_ref,xco2_ref,&
                           xsoot_ref,xco_ref,xch4_ref,T_ref,ng_ref,len2)
      case default
         call shutdown('SLW-1 method undefined')
      endselect
      
      !---------------------------------------
      !Solve for the absorption and weighting
      !coefficients at the reference state
      !---------------------------------------
      call solve_slw1_parameters(slw1_method,val1,val2,len1,len2,&      !Absorption and weighting coefficients
                                 a_ref(1),kappa_ref(1))                 !  of the gray gas
      kappa_ref(0) = 0._dp                                              !Absorption coefficient of the transparent window
      a_ref(0) = 1._dp - a_ref(1)                                       !Weighting coefficient of the transparent window

      gas_loop: do jgas=0,1
         !Check if the gas is a transparent window
         transparent_window = .false.
         if (jgas.eq.0) transparent_window = .true.

         !-------------------------------------
         !Computing local radiative properties
         !-------------------------------------           
         do i=0,nx+1
            do j=0,ny+1
               do k=0,nz+1
                  !Mole fraction of the mixture of participating species
                  xmix = 0._dp
                  if ((nch4.gt.0).and.(albdf_ch4_name.ne.'null')) &
                     xmix = xmix + xch4(i,j,k)
                  if ((nco.gt.0).and.(albdf_co_name.ne.'null')) &
                     xmix = xmix + xco(i,j,k)   
                  if ((nco2.gt.0).and.(albdf_co2_name.ne.'null')) &
                     xmix = xmix + xco2(i,j,k)    
                  if ((nh2o.gt.0).and.(albdf_h2o_name.ne.'null')) &
                     xmix = xmix + xh2o(i,j,k) 
                  if ((nsoot.gt.0).and.(albdf_soot_name.ne.'null')) &
                     xmix = xmix + xsoot(i,j,k) 
                  
                  selectcase(trim(slw_nonuniform_method))
                  case('uniform')                                       !Uniform medium
                     kappa(i,j,k) = kappa_ref(jgas)
                     a(i,j,k) = a_ref(jgas)
                  
                  case('reference_approach1')                           !Solovjov et al., 2011 (JHT)
                     if (transparent_window) then
                        C_ref = inverse_albdf(a_ref(0),T_ref,xh2o_ref,& !Supplemental absorption cross-section
                           xco2_ref,xsoot_ref,xco_ref,xch4_ref,T_ref)   !  at the reference state       
                        a(i,j,k) = get_albdf(T_ref,xh2o_ref,xco2_ref,&  !Weighting coefficient of the
                           xsoot_ref,xco_ref,xch4_ref,T(i,j,k),C_ref)   !  transparent window
                        kappa(i,j,k) = 0._dp
                     else
                        C_ref = slw_kappa_func(T_ref,p_ref,xmix_ref,&   !Absorption cross-section at 
                                               kappa_ref(1),.true.)     !  the reference state
                        F_ref = get_albdf(T_ref,xh2o_ref,xco2_ref,&     !ALBDF at the reference state
                           xsoot_ref,xco_ref,xch4_ref,T(i,j,k),C_ref)
                        C_lcl = inverse_albdf(F_ref,T(i,j,k),&          !Absorption cross-section
                           xh2o(i,j,k),xco2(i,j,k),xsoot(i,j,k),&       !  at the local state
                           xco(i,j,k),xch4(i,j,k),T_ref)
                        kappa(i,j,k) = slw_kappa_func(T(i,j,k),&        !Absorption coefficient
                           p(i,j,k),xmix,C_lcl)
                        a(i,j,k) = 1._dp - a(i,j,k)                     !Weighting coefficient
                     endif
                  case('reference_approach2')                           !Alternative approach mentioned in 
                                                                        !  Solovjov et al., 2011 (JQSRT)
                     if (transparent_window) then
                        C_ref = inverse_albdf(a_ref(0),T_ref,xh2o_ref,& !Supplemental absorption cross-section
                           xco2_ref,xsoot_ref,xco_ref,xch4_ref,T_ref)   !  at the reference state
                        a(i,j,k) = get_albdf(T(i,j,k),xh2o(i,j,k),&     !Weighting coefficient of the
                           xco2(i,j,k),xsoot(i,j,k),xco(i,j,k),&        !  transparent window
                           xch4(i,j,k),T(i,j,k),C_ref)   
                        kappa(i,j,k) = 0._dp
                     else
                        C_ref = slw_kappa_func(T_ref,p_ref,xmix_ref,&   !Absorption cross-section at 
                                               kappa_ref(1),.true.)     !  the reference state
                        F_ref = get_albdf(T_ref,xh2o_ref,xco2_ref,&     !ALBDF at the reference state
                           xsoot_ref,xco_ref,xch4_ref,T(i,j,k),C_ref)
                        C_lcl = inverse_albdf(F_ref,T(i,j,k),&          !Absorption cross-section
                           xh2o(i,j,k),xco2(i,j,k),xsoot(i,j,k),&       !  at the local state
                           xco(i,j,k),xch4(i,j,k),T_ref)
                        kappa(i,j,k) = slw_kappa_func(T(i,j,k),&        !Absorption coefficient
                           p(i,j,k),xmix,C_lcl)
                        a(i,j,k) = 1._dp - a(i,j,k)                     !Weighting coefficient
                     endif
               
                  case default
                     call shutdown('slw_nonuniform_method unspecified')
                  endselect
                  aIb(i,j,k) = a(i,j,k)*Ib(i,j,k)                       !Emission term
                  
                  !Store properties of each gray gas if requested
                  if (slw_print_properties) then
                     kappa_j(i,j,k,jgas) = kappa(i,j,k)
                     a_j(i,j,k,jgas) = a(i,j,k)
                  endif
               enddo
            enddo
         enddo 

         !--------------------------
         !Solve the radiation field
         !--------------------------
         !Starting the processing time counter
         call cpu_time(start_proc_time)
   
         !Solving the radiation field for gas jgas
         call fvm_solution(kappa,kappa,kappa_sct,aIb,phi,x_eps,y_eps,&
                           z_eps,flux_j,source_j,.true.)
                           
         !Updating the heat flux and heat source
         flux = flux + flux_j
         source = source + source_j
      
         !Stopping the processing time counter and adding to the total
         call cpu_time(end_proc_time)
         proc_time = proc_time + end_proc_time - start_proc_time
      
      enddo gas_loop
      
      !--------------------------------------------
      !Manage dump of gas properties, if necessary
      !--------------------------------------------
      if (slw_print_properties) then
         !Dump data
         do i=0,nx+1
            do j=0,ny+1
               do k=0,nz+1
                  write(uout,'(1000(e26.15e3,:,","))')&
                     x(i),y(j),z(k),(kappa_j(i,j,k,jgas),&
                        a_j(i,j,k,jgas),jgas=0,1)
               enddo
            enddo
         enddo
         
         !Close output unit
         if (slw_print_properties) close(uout)
      
         !Deallocate extra arrays
         deallocate(kappa_j)
         deallocate(a_j)
      endif
      
      !-----------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1
               if (compute_emission) &
                  emission(i,j,k) = fourpi*kappa(i,j,k)*aIb(i,j,k)      !Total emission
               if (compute_absorption) &
                  absorption(i,j,k) = source(i,j,k) + &
                                      fourpi*kappa(i,j,k)*aIb(i,j,k)    !Total absorption
               if (compute_planck) &
                  kappaPlanck(i,j,k) = kappa(i,j,k)*a(i,j,k)            !Planck-mean absorption coefficient
            enddo
         enddo
      enddo
      
      !------------------
      !Deallocate arrays
      !------------------
      deallocate(a)
      deallocate(aIb)
      deallocate(Ib)
      deallocate(kappa)
      deallocate(kappa_sct)
      deallocate(phi)
      deallocate(x_eps)
      deallocate(y_eps)
      deallocate(z_eps)
      
		!---------------------
      !Total computing time
      !---------------------
      call cpu_time(end_total_time)
      total_time = end_total_time - start_total_time
      
   endsubroutine slw1_black_solution

   !==================================================================
   !Subroutine for the radiative heat transfer solution via the
   !SLW-1 model for a non-scattering medium bounded by black surfaces
   !==================================================================
   subroutine bslw1_black_solution(flux,source,proc_time,total_time,&
                                   absorption,emission,kappaPlanck,&
                                   albdf_ready)
   
      !-------------------------
      !Declaration of variables
      !-------------------------
      use comp_functions, only: CheckFileExists,CheckMemAlloc,&
                                get_file_unit,shutdown
      use physical_functions, only: bb_emission_frac
      integer :: n_bands,ng_ref,nl,nx,ny,nz
      integer :: nch4,nco,nco2,nh2o,nsoot
      integer :: i,iband,j,jgas,k
      integer :: ierr,uout
      real(dp),intent(out) :: proc_time,total_time,&
         flux(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
         source(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp),intent(out),optional :: &
         absorption(0:xcells+1,0:ycells+1,0:zcells+1),&
         emission(0:xcells+1,0:ycells+1,0:zcells+1),&
         kappaPlanck(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp) :: flux_j(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
                  source_j(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp) :: a_ref(0:1),kappa_ref(0:1)
      real(dp) :: C_ref,C_lcl,F_ref
      real(dp) :: p_ref,T_ref,xh2o_ref,xco2_ref,xsoot_ref,xco_ref,&
                  xch4_ref,xmix_ref
      real(dp) :: xmix
      real(dp) :: len1,len2,val1,val2
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: a,F,FaIb,Ib,kappa,&
         kappa_sct,kIb,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: a_j,kappa_j
      logical,optional :: albdf_ready
      logical :: albdf_read,transparent_window
      logical :: compute_absorption=.false.,compute_emission=.false.,&
                 compute_planck=.false.
   
      !Starting the total computing time counter
      call cpu_time(start_total_time)
      
      !Zeroing out processing time counter 
      !(necessary for the posterior sum over 
      !all instances of the RTE solution)
      proc_time = 0._dp
      
      !Set flags for optional output arguments
      if (present(absorption))   compute_absorption = .true.
      if (present(emission))     compute_emission = .true.
      if (present(kappaPlanck))  compute_planck = .true.
      
      !---------------------------
      !Set up optional parameters
      !---------------------------
      if (present(albdf_ready)) then
         albdf_read = albdf_ready
      else
         albdf_read = .false.                                           !By default, assume that the ALBDF data
      endif                                                             !  has not been read yet
      
      !----------------
      !Surrogate names
      !----------------
      nx = xcells                                                       !Number of cells along direction x
      ny = ycells                                                       !Number of cells along direction y
      nz = zcells                                                       !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization

      nch4  = albdf_ch4_nx                                              !Number of ALBDF mole fraction values for CH4
      nco   = albdf_co_nx                                               !Number of ALBDF mole fraction values for CO
      nco2  = albdf_co2_nx                                              !Number of ALBDF mole fraction values for CO2
      nh2o  = albdf_h2o_nx                                              !Number of ALBDF mole fraction values for H2O
      nsoot = albdf_soot_nx                                             !Number of ALBDF mole fraction values for soot
      
      len1 = slw1_length1                                               !First characteristic length
      len2 = slw1_length2                                               !Second characteristic length
      ng_ref = slw1_ngases                                              !Number of gray gases for the reference SLW solution

      n_bands = albdf_nbands                                            !Number of wide bands in the solution

      !---------------------------------------------
      !Prepare for dumping properties, if requested
      !---------------------------------------------
      if (slw_print_properties) then
         !Allocate extra arrays
         allocate(a_j(0:nx+1,0:ny+1,0:nz+1,0:1),stat=ierr)
         call CheckMemAlloc('a_j',ierr)
         allocate(kappa_j(0:nx+1,0:ny+1,0:nz+1,0:1),stat=ierr)
         call CheckMemAlloc('kappa_j',ierr)
         
         !Prepare output unit
         if (trim(slw_print_file).eq.'null') &                          !Check if a file name was given                          
            call shutdown('slw_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(slw_print_file))                      !Open file
         write(uout,'(100(a,:,","))') 'lambda-','lambda+',&             !Write header
            'x','y','z','F','kappa_0','a_0','kappa_1','a_1'
         write(uout,'(100(a,:,","))') '[mum]','[mum]','[m]','[m]',&     !Write subheader
                                 '[m]','[-]','[1/m]','[-]','[1/m]','[-]'
      endif
      
      !------------------
      !Allocating arrays
      !------------------
      allocate(a(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('a',ierr)
      allocate(F(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('F',ierr)
      allocate(Ib(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(faIb(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('faIb',ierr)
      allocate(kappa(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(phi(1:nl,1:nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(x_eps(1:2,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(1:2,0:nx+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(1:2,0:nx+1,0:ny+1),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !--------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !--------------------------------------------------------
      !Blackbody intensity
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1
               Ib(i,j,k) = sigma*(T(i,j,k)**4._dp)/pi
            enddo
         enddo
      enddo
      
      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !Set wall emissivities to unit for all wall cells
      x_eps = 1._dp
      y_eps = 1._dp
      z_eps = 1._dp
      
      !-----------------------------------
      !If needed, read the ALBDF database
      !-----------------------------------
      if (trim(slw_mixture_method).eq.'albdf_precombined') then
         if (.not.albdf_read) call read_aldbf_mix
      else
         if (.not.albdf_read) call read_aldbf
      endif

      !---------------------------
      !Define the reference state
      !---------------------------
      call get_slw_reference_state(T_ref,p_ref,xh2o_ref,xco2_ref,&
                                   xsoot_ref,xco_ref,xch4_ref)
      xmix_ref = xh2o_ref + xco2_ref + xsoot_ref + xco_ref + xch4_ref   !Reference mixture mole fraction

      !----------------------------------
      !Initializing arrays and variables
      !----------------------------------
      flux = 0._dp
      source = 0._dp
      kIb = 0._dp

      band_loop: do iband=1,n_bands
         !---------------------------------------------------------
         !Compute the fraction of blackbody energy within the band
         !---------------------------------------------------------
         do i=0,nx+1
            do j=0,ny+1
               do k=0,nz+1
                  F(i,j,k) = &
                     bb_emission_frac(bslw_ubound(iband),T(i,j,k)) - &
                     bb_emission_frac(bslw_lbound(iband),T(i,j,k))
               enddo
            enddo
         enddo
            
         !----------------------------------
         !Compute the two target quantities
         !----------------------------------
         selectcase(trim(slw1_method))
         case('kappa-epsilon')
            val1 = get_slw_kp(T_ref,p_ref,xh2o_ref,xco2_ref,xsoot_ref,&
                              xco_ref,xch4_ref,T_ref,ng_ref,.false.,&
                              iband)
            val2 = get_slw_emissivity(T_ref,p_ref,xh2o_ref,xco2_ref,&
                                      xsoot_ref,xco_ref,xch4_ref,T_ref,&
                                      ng_ref,len1,.false.,iband)
         case('epsilon-epsilon')
            val1 = get_slw_emissivity(T_ref,p_ref,xh2o_ref,xco2_ref,&
                                      xsoot_ref,xco_ref,xch4_ref,T_ref,&
                                      ng_ref,len1,.false.,iband)
            val2 = get_slw_emissivity(T_ref,p_ref,xh2o_ref,xco2_ref,&
                                      xsoot_ref,xco_ref,xch4_ref,T_ref,&
                                      ng_ref,len2,.false.,iband)
         case default
            call shutdown('SLW-1 method undefined')
         endselect
      
         !---------------------------------------
         !Solve for the absorption and weighting
         !coefficients at the reference state
         !---------------------------------------
         write(*,*) iband,val1,val2
         call solve_slw1_parameters(slw1_method,val1,val2,len1,len2,&   !Absorption and weighting coefficients
                                    a_ref(1),kappa_ref(1))              !  of the gray gas
         kappa_ref(0) = 0._dp                                           !Absorption coefficient of the transparent window
         a_ref(0) = 1._dp - a_ref(1)                                    !Weighting coefficient of the transparent window

         gas_loop: do jgas=0,1
            !Check if the gas is a transparent window
            transparent_window = .false.
            if (jgas.eq.0) transparent_window = .true.

            !-------------------------------------
            !Computing local radiative properties
            !-------------------------------------           
            do i=0,nx+1
               do j=0,ny+1
                  do k=0,nz+1
                     !Mole fraction of the mixture 
                     !of participating species
                     xmix = 0._dp
                     if ((nch4.gt.0).and.(albdf_ch4_name.ne.'null')) &
                        xmix = xmix + xch4(i,j,k)
                     if ((nco.gt.0).and.(albdf_co_name.ne.'null')) &
                        xmix = xmix + xco(i,j,k)   
                     if ((nco2.gt.0).and.(albdf_co2_name.ne.'null')) &
                        xmix = xmix + xco2(i,j,k)    
                     if ((nh2o.gt.0).and.(albdf_h2o_name.ne.'null')) &
                        xmix = xmix + xh2o(i,j,k) 
                     if ((nsoot.gt.0).and.(albdf_soot_name.ne.'null')) &
                        xmix = xmix + xsoot(i,j,k) 
                  
                     selectcase(trim(slw_nonuniform_method))
                     case('uniform')                                    !Uniform medium
                        kappa(i,j,k) = kappa_ref(jgas)
                        a(i,j,k) = a_ref(jgas)
                     
                     case('reference_approach1')                        !Solovjov et al., 2011 (JHT)
                        if (transparent_window) then
                           C_ref = inverse_albdf(a_ref(0),T_ref,&       !Supplemental absorption cross-section
                              xh2o_ref,xco2_ref,xsoot_ref,xco_ref,&     !  at the reference state       
                              xch4_ref,T_ref,iband)
                           a(i,j,k) = get_albdf(T_ref,xh2o_ref,&        !Weighting coefficient of the
                              xco2_ref,xsoot_ref,xco_ref,xch4_ref,&     !  transparent window
                              T(i,j,k),C_ref,iband)
                           kappa(i,j,k) = 0._dp
                        else
                           C_ref = slw_kappa_func(T_ref,p_ref,xmix_ref,&!Absorption cross-section at 
                                                  kappa_ref(1),.true.)  !  the reference state
                           F_ref = get_albdf(T_ref,xh2o_ref,xco2_ref,&  !ALBDF at the reference state
                              xsoot_ref,xco_ref,xch4_ref,T(i,j,k),&
                              C_ref,iband)
                           C_lcl = inverse_albdf(F_ref,T(i,j,k),&       !Absorption cross-section
                              xh2o(i,j,k),xco2(i,j,k),xsoot(i,j,k),&    !  at the local state
                              xco(i,j,k),xch4(i,j,k),T_ref,iband)
                           kappa(i,j,k) = slw_kappa_func(T(i,j,k),&     !Absorption coefficient
                              p(i,j,k),xmix,C_lcl)
                           a(i,j,k) = 1._dp - a(i,j,k)                  !Weighting coefficient
                        endif
                     case('reference_approach2')                        !Alternative approach mentioned in 
                                                                        !  Solovjov et al., 2011 (JQSRT)
                        if (transparent_window) then
                           C_ref = inverse_albdf(a_ref(0),T_ref,&       !Supplemental absorption cross-section
                              xh2o_ref,xco2_ref,xsoot_ref,xco_ref,&     !  at the reference state
                              xch4_ref,T_ref)
                           a(i,j,k) = get_albdf(T(i,j,k),xh2o(i,j,k),&  !Weighting coefficient of the
                              xco2(i,j,k),xsoot(i,j,k),xco(i,j,k),&     !  transparent window
                              xch4(i,j,k),T(i,j,k),C_ref,iband)   
                           kappa(i,j,k) = 0._dp
                        else
                           C_ref = slw_kappa_func(T_ref,p_ref,xmix_ref,&!Absorption cross-section at 
                                                kappa_ref(1),.true.)    !  the reference state
                           F_ref = get_albdf(T_ref,xh2o_ref,xco2_ref,&  !ALBDF at the reference state
                              xsoot_ref,xco_ref,xch4_ref,T(i,j,k),&
                              C_ref,iband)
                           C_lcl = inverse_albdf(F_ref,T(i,j,k),&       !Absorption cross-section
                              xh2o(i,j,k),xco2(i,j,k),xsoot(i,j,k),&    !  at the local state
                              xco(i,j,k),xch4(i,j,k),T_ref,iband)
                           kappa(i,j,k) = slw_kappa_func(T(i,j,k),&     !Absorption coefficient
                              p(i,j,k),xmix,C_lcl)
                           a(i,j,k) = 1._dp - a(i,j,k)                  !Weighting coefficient
                        endif

                     case default
                        call shutdown('slw_nonuniform_method unspecified')
                     endselect
                     faIb(i,j,k) = F(i,j,k)*a(i,j,k)*Ib(i,j,k)          !Emission term
                     kIb(i,j,k) = kIb(i,j,k) + kappa(i,j,k)*faIb(i,j,k) !Total RTE emission term
                     
                     !Store properties of each gray gas if requested
                     if (slw_print_properties) then
                        kappa_j(i,j,k,jgas) = kappa(i,j,k)
                        a_j(i,j,k,jgas) = a(i,j,k)
                     endif
                  enddo
               enddo
            enddo 

            !--------------------------
            !Solve the radiation field
            !--------------------------
            !Starting the processing time counter
            call cpu_time(start_proc_time)

            !Solving the radiation field for gas jgas
            call fvm_solution(kappa,kappa,kappa_sct,faIb,phi,x_eps,&
                              y_eps,z_eps,flux_j,source_j,.true.)
                              
            !Updating the heat flux and heat source
            flux = flux + flux_j
            source = source + source_j

            !Stopping the processing time counter and adding to the total
            call cpu_time(end_proc_time)
            proc_time = proc_time + end_proc_time - start_proc_time

         enddo gas_loop

         !-------------------------------------
         !Dump of gas properties, if necessary
         !-------------------------------------
         if (slw_print_properties) then
            do i=0,nx+1
               do j=0,ny+1
                  do k=0,nz+1
                     write(uout,'(1000(e26.15e3,:,","))')&
                        bslw_lbound(iband)*1.e6_dp,&
                        bslw_ubound(iband)*1.e6_dp,&
                        x(i),y(j),z(k),F(i,j,k),&
                        (kappa_j(i,j,k,jgas),a_j(i,j,k,jgas),jgas=0,1)
                  enddo
               enddo
            enddo
         endif
      
      enddo band_loop
      
      !-------------------------------------
      !Finished up with the dump properties
      !-------------------------------------
      if (slw_print_properties) then
         !Close output unit
         if (slw_print_properties) close(uout)
      
         !Deallocate extra arrays
         deallocate(kappa_j)
         deallocate(a_j)
      endif
      
      !-----------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1
               if (compute_emission) &
                  emission(i,j,k) = fourpi*kIb(i,j,k)                   !Total emission
               if (compute_absorption) &
                  absorption(i,j,k) = source(i,j,k) + fourpi*kIb(i,j,k) !Total absorption
               if (compute_planck) &
                  kappaPlanck(i,j,k) = kIb(i,j,k)/Ib(i,j,k)             !Planck-mean absorption coefficient
            enddo
         enddo
      enddo
      
      !------------------
      !Deallocate arrays
      !------------------
      deallocate(a)
      deallocate(F)
      deallocate(FaIb)
      deallocate(Ib)
      deallocate(kappa)
      deallocate(kappa_sct)
      deallocate(kIb)
      deallocate(phi)
      deallocate(x_eps)
      deallocate(y_eps)
      deallocate(z_eps)
      
		!---------------------
      !Total computing time
      !---------------------
      call cpu_time(end_total_time)
      total_time = end_total_time - start_total_time
      
   endsubroutine bslw1_black_solution












   !==============================================================
   !Subroutine for the radiative heat transfer solution
   !via the multiple integration SLW model for a non-
   !scattering medium bounded by black walls
   !(only for testing purposes)
   !==============================================================
   subroutine slw_multint_black_solution2(n_gases,flux,source,proc_time,&
                                         total_time,albdf_ready)
	
		!-------------------------
		!Declaration of variables
		!-------------------------
      use comp_functions, only: CheckMemAlloc,shutdown
      real(dp),intent(out) :: proc_time,total_time,&
         flux(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
         source(0:xcells+1,0:ycells+1,0:zcells+1)
      integer,intent(in) :: n_gases
      integer :: i,j,k
      integer :: jco2,jh2o
      integer :: ierr
      integer :: nx,ny,nz,nl
      integer :: nch4,nco,nco2,nh2o,nsoot
      real(dp) :: flux_j(1:3,0:xcells+1,0:ycells+1,0:zcells+1),&
                  source_j(0:xcells+1,0:ycells+1,0:zcells+1)
      real(dp) :: p_ref,T_ref,xh2o_ref,xco2_ref,xsoot_ref,xco_ref,&
                  xch4_ref
      real(dp) :: cj_loc
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref,&
         Fj_ref,Fj_sup_ref
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: aIb,Ib,kappa,kappa_sct,&
         x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: cj_sup_loc,Fj_sup_loc,&
         ah2o,aco2,kappah2o,kappaco2
      logical,optional :: albdf_ready
      logical :: albdf_read,transparent_window
      
      !Starting the total computing time counter
      call cpu_time(start_total_time)
      
      !Zeroing out processing time counter 
      !(necessary for the posterior sum over 
      !all instances of the RTE solution)
      proc_time = 0._dp
      
      !---------------------------
      !Set up optional parameters
      !---------------------------
      if (present(albdf_ready)) then
         albdf_read = albdf_ready
      else
         albdf_read = .false.                                           !By default, assume that the ALBDF data
      endif                                                             !  has not been read yet
      
      !----------------
      !Surrogate names
      !----------------
      nx = xcells                                                       !Number of cells along direction x
      ny = ycells                                                       !Number of cells along direction y
      nz = zcells                                                       !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization

      nch4  = albdf_ch4_nx                                              !Number of ALBDF mole fraction values for CH4
      nco   = albdf_co_nx                                               !Number of ALBDF mole fraction values for CO
      nco2  = albdf_co2_nx                                              !Number of ALBDF mole fraction values for CO2
      nh2o  = albdf_h2o_nx                                              !Number of ALBDF mole fraction values for H2O
      nsoot = albdf_soot_nx                                             !Number of ALBDF mole fraction values for soot

      !----------------------------
      !Allocating all other arrays
      !----------------------------     
      allocate(ah2o(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
      call CheckMemAlloc('ah2o',ierr)
      allocate(aco2(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
      call CheckMemAlloc('aco2',ierr)
      allocate(aIb(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(cj_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('cj_ref',ierr)
      allocate(cj_sup_loc(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
      call CheckMemAlloc('cj_sup_loc',ierr)
      allocate(cj_sup_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('cj_sup_ref',ierr)
      allocate(Fj_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('Fj_ref',ierr)
      allocate(Fj_sup_ref(0:n_gases),stat=ierr)
      call CheckMemAlloc('Fj_sup_ref',ierr)
      allocate(Fj_sup_loc(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
      call CheckMemAlloc('Fj_sup_loc',ierr)
      allocate(Ib(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappah2o(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
      call CheckMemAlloc('kappah2o',ierr)
      allocate(kappaco2(0:nx+1,0:ny+1,0:nz+1,0:n_gases),stat=ierr)
      call CheckMemAlloc('kappaco2',ierr)
      allocate(kappa_sct(0:nx+1,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(phi(1:nl,1:nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(x_eps(1:2,0:ny+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(1:2,0:nx+1,0:nz+1),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(1:2,0:nx+1,0:ny+1),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !------------------------------------------------------------------------------
      !Defining absorption cross-sections and supplementar absorption cross-sections
      !------------------------------------------------------------------------------
      if (trim(slw_nonuniform_method).eq.'rank_correlated') then
         call get_slw_Fj(n_gases,Fj_ref(1:n_gases),Fj_sup_ref(1:n_gases))
         Fj_sup_ref(0) = slw_Fmin  
      else
         !For other SLWs, divide the reference Cj
         call get_slw_cj(slw_cmin,slw_cmax,n_gases,cj_ref(1:n_gases),&
            cj_sup_ref(1:n_gases))
         cj_sup_ref(0) = slw_cmin
      endif

      !-----------------------------------
      !If needed, read the ALBDF database
      !-----------------------------------
      if (.not.albdf_read) call read_aldbf

      !-----------------------------------------
      !Define the reference state, if necessary
      !-----------------------------------------
      if (trim(slw_nonuniform_method).ne.'uniform') &
         call get_slw_reference_state(T_ref,p_ref,xh2o_ref,xco2_ref,&
                                      xsoot_ref,xco_ref,xch4_ref)
                                      
      !----------------------------------
      !Initializing arrays and variables
      !----------------------------------
      flux = 0._dp
      source = 0._dp

      !--------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !--------------------------------------------------------
      !Blackbody intensity
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1
               Ib(i,j,k) = sigma*(T(i,j,k)**4._dp)/pi
            enddo
         enddo
      enddo
      
      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !Set wall emissivities to unit for all wall cells
      x_eps = 1._dp
      y_eps = 1._dp
      z_eps = 1._dp

		!-----------------------------
		!Defining gray gas properties
		!-----------------------------
      do i=0,nx+1
         do j=0,ny+1
            do k=0,nz+1    
               do jh2o = 0,n_gases
                  !Check if the gas is a transparent window
                  transparent_window = .false.
                  if (jh2o.eq.0) transparent_window = .true.
                  
                  selectcase(trim(slw_nonuniform_method))
                  case('uniform')
                     !+++++++++++++++
                     !Uniform medium
                     !+++++++++++++++
                     !(for uniform media, cj_ref = cj_loc, 
                     !so use cj_ref always)
                     if (transparent_window) then
                        kappah2o(i,j,k,jh2o) = 0._dp                    !Absorption coefficient of the
                                                                        !  transparent window
                        ah2o(i,j,k,jh2o) = slw_a_func(T(i,j,k),&        !Weighting coefficient of the
                           xh2o(i,j,k),0._dp,0._dp,0._dp,0._dp,&        !  transparent window
                           T(i,j,k),cj_sup_ref(jh2o))
                        
                     else
                        kappah2o(i,j,k,jh2o) = slw_kappa_func(T(i,j,k),&!Absorption coefficient of the 
                           p(i,j,k),xh2o(i,j,k),cj_ref(jh2o))           !  gray gas
                        ah2o(i,j,k,jh2o) = slw_a_func(T(i,j,k),&        !Weighting coefficient of the
                           xh2o(i,j,k),0._dp,0._dp,0._dp,0._dp,&        !  gray gas
                           T(i,j,k),cj_sup_ref(jh2o),cj_sup_ref(jh2o-1))
                     endif               
            
                  case('scaled')
                     !+++++++++++
                     !Scaled-SLW
                     !+++++++++++
                     call shutdown('Scaled model not available for now')
                  
                  case('reference_approach')
                     !+++++++
                     !RA-SLW
                     !+++++++
                     if (transparent_window) then
                        kappah2o(i,j,k,jh2o) = 0._dp                         !Absorption coefficient of the
                                                                        !  transparent window
                        ah2o(i,j,k,jh2o) = slw_a_func(T_ref,xh2o_ref,0._dp,& !Weighting coefficient of the
                           0._dp,0._dp,0._dp,T(i,j,k),cj_sup_ref(jh2o)) !  transparent window
                    else
                        Fj_sup_ref(jh2o) = get_albdf(T_ref,xh2o_ref,&   !ALBDF evaluated at the
                           0._dp,0._dp,0._dp,0._dp,T_ref,&              !  reference state
                           cj_sup_ref(jh2o))
                        cj_sup_loc(i,j,k,jh2o) = inverse_albdf(&        !Get the local supplementar cross-section
                           Fj_sup_ref(jh2o),T(i,j,k),xh2o(i,j,k),&      !  by solving the implicit equation
                           0._dp,0._dp,0._dp,0._dp,T_ref)
                        cj_loc = get_slw_single_cj(&                    !Local cross-section determined
                           cj_sup_loc(i,j,k,jh2o-1),&                   !  from the supplementar ones
                              cj_sup_loc(i,j,k,jh2o))                        
                        kappah2o(i,j,k,jh2o) = slw_kappa_func(T(i,j,k),&     !Absorption coefficient of the 
                           p(i,j,k),xh2o(i,j,k),cj_loc)                 !  gray gas
                        ah2o(i,j,k,jh2o) = slw_a_func(T_ref,xh2o_ref,0._dp,& !Weighting coefficient of the
                           0._dp,0._dp,0._dp,T(i,j,k),cj_sup_ref(jh2o),&!  gray gas
                           cj_sup_ref(jh2o-1))    
                     endif
               
                  case('rank_correlated')
                     !+++++++
                     !RC-SLW
                     !+++++++
                     cj_sup_loc(i,j,k,jh2o) = inverse_albdf(&           !Get the local supplementar cross-section
                        Fj_sup_ref(jh2o),T(i,j,k),xh2o(i,j,k),0._dp,&   !  by solving the implicit equation using
                           0._dp,0._dp,0._dp,T_ref)                     !  the previously divided supplementar Fs
                     if (transparent_window) then
                        kappah2o(i,j,k,jh2o) = 0._dp                         !Absorption coefficient of the
                                                                        !  transparent window
                        ah2o(i,j,k,jh2o) = slw_a_func(T(i,j,k),xh2o(i,j,k),& !Weighting coefficient of the
                           0._dp,0._dp,0._dp,0._dp,T(i,j,k),&           !  transparent window
                           cj_sup_loc(i,j,k,jh2o))
                      else
                         cj_loc = inverse_albdf(Fj_ref(jh2o),T(i,j,k),& !Get the local cross-section by solving 
                            xh2o(i,j,k),0._dp,0._dp,0._dp,0._dp,T_ref)  !  the implicit equation using the
                                                                        !  the previously defined Fs
                        kappah2o(i,j,k,jh2o) = slw_kappa_func(T(i,j,k),&     !Absorption coefficient of the 
                           p(i,j,k),xh2o(i,j,k),cj_loc)                 !  gray gas
                        ah2o(i,j,k,jh2o) = slw_a_func(T(i,j,k),xh2o(i,j,k),& !Weighting coefficient of the
                           0._dp,0._dp,0._dp,0._dp,T(i,j,k),&           !  gray gas
                           cj_sup_loc(i,j,k,jh2o),cj_sup_loc(i,j,k,jh2o-1))
                     endif
         
                  case('locally_correlated')
                     !+++++++
                     !LC-SLW
                     !+++++++
                     Fj_sup_loc(i,j,k,jh2o) = get_albdf(T_ref,&         !Local ALBDF evaluated using
                        xh2o_ref,0._dp,0._dp,0._dp,0._dp,T(i,j,k),&     !  reference state
                        cj_sup_ref(jh2o))
                     if (transparent_window) then
                        kappah2o(i,j,k,jh2o) = 0._dp                         !Absorption coefficient of the
                                                                        !  transparent window
                        ah2o(i,j,k,jh2o) = Fj_sup_loc(i,j,k,jh2o)            !Weighting coefficient of the
                                                                        !  transparent window
                     else
                        cj_sup_loc(i,j,k,jh2o) = inverse_albdf(&        !Get the local supplementar cross-section
                           Fj_sup_loc(i,j,k,jh2o),T(i,j,k),xh2o(i,j,k),&!  by solving the implicit equation
                           0._dp,0._dp,0._dp,0._dp,T(i,j,k))
                           
                        cj_loc = get_slw_single_cj(&                    !Local cross-section determined        
                           cj_sup_loc(i,j,k,jh2o-1),&                   !  from the supplementar ones
                           cj_sup_loc(i,j,k,jh2o))
                        kappah2o(i,j,k,jh2o) = slw_kappa_func(&              !Absorption coefficient of the 
                           T(i,j,k),p(i,j,k),xh2o(i,j,k),cj_loc)        !  gray gas
                        ah2o(i,j,k,jh2o) = Fj_sup_loc(i,j,k,jh2o) - &        !Weighting coefficient of the gray gas
                           Fj_sup_loc(i,j,k,jh2o-1)
                     endif
                  
                  case default
                     call shutdown('Problem in the specification of &
                        &slw_nonuniform_method')
                  endselect        
                  
               enddo

               do jco2 = 0,n_gases
                  !Check if the gas is a transparent window
                  transparent_window = .false.
                  if (jco2.eq.0) transparent_window = .true.
                  
                  selectcase(trim(slw_nonuniform_method))
                  case('uniform')
                     !+++++++++++++++
                     !Uniform medium
                     !+++++++++++++++
                     !(for uniform media, cj_ref = cj_loc, 
                     !so use cj_ref always)
                     if (transparent_window) then
                        kappaco2(i,j,k,jco2) = 0._dp                         !Absorption coefficient of the
                                                                        !  transparent window
                        aco2(i,j,k,jco2) = slw_a_func(T(i,j,k),0._dp,&       !Weighting coefficient of the
                           xco2(i,j,k),0._dp,0._dp,0._dp,T(i,j,k),&     !  transparent window
                           cj_sup_ref(jco2))
                        
                     else
                        kappaco2(i,j,k,jco2) = slw_kappa_func(T(i,j,k),&     !Absorption coefficient of the 
                           p(i,j,k),xco2(i,j,k),cj_ref(jco2))           !  gray gas
                        aco2(i,j,k,jco2) = slw_a_func(T(i,j,k),0._dp,&       !Weighting coefficient of the
                           xco2(i,j,k),0._dp,0._dp,0._dp,T(i,j,k),&     !  gray gas
                           cj_sup_ref(jco2),cj_sup_ref(jco2-1))
                     endif               

                  case('scaled')
                     !+++++++++++
                     !Scaled-SLW
                     !+++++++++++
                     call shutdown('Scaled model not available for now')
                  
                  case('reference_approach')
                     !+++++++
                     !RA-SLW
                     !+++++++
                     if (transparent_window) then
                        kappaco2(i,j,k,jco2) = 0._dp                         !Absorption coefficient of the
                                                                        !  transparent window
                        aco2(i,j,k,jco2) = slw_a_func(T_ref,0._dp,xco2_ref,& !Weighting coefficient of the
                           0._dp,0._dp,0._dp,T(i,j,k),cj_sup_ref(jco2)) !  transparent window
                    else
                        Fj_sup_ref(jco2) = get_albdf(T_ref,0._dp,&      !ALBDF evaluated at the
                           xco2_ref,0._dp,0._dp,0._dp,T_ref,&           !  reference state
                           cj_sup_ref(jco2))
                        cj_sup_loc(i,j,k,jco2) = inverse_albdf(&        !Get the local supplementar cross-section
                           Fj_sup_ref(jco2),T(i,j,k),0._dp,xco2(i,j,k),&!  by solving the implicit equation
                           0._dp,0._dp,0._dp,T_ref)
                        cj_loc = get_slw_single_cj(&                    !Local cross-section determined
                           cj_sup_loc(i,j,k,jco2-1),&                   !  from the supplementar ones
                           cj_sup_loc(i,j,k,jco2))                        
                        kappaco2(i,j,k,jco2) = slw_kappa_func(T(i,j,k),&     !Absorption coefficient of the 
                           p(i,j,k),xco2(i,j,k),cj_loc)                 !  gray gas
                        aco2(i,j,k,jco2) = slw_a_func(T_ref,0._dp,xco2_ref,& !Weighting coefficient of the
                           0._dp,0._dp,0._dp,T(i,j,k),cj_sup_ref(jco2),&!  gray gas
                           cj_sup_ref(jco2-1))
                     endif
               
                  case('rank_correlated')
                     !+++++++
                     !RC-SLW
                     !+++++++
                     cj_sup_loc(i,j,k,jco2) = inverse_albdf(&           !Get the local supplementar cross-section
                        Fj_sup_ref(jco2),T(i,j,k),0._dp,xco2(i,j,k),&   !  by solving the implicit equation using
                           0._dp,0._dp,0._dp,T_ref)                     !  the previously divided supplementar Fs
                     if (transparent_window) then
                        kappaco2(i,j,k,jco2) = 0._dp                         !Absorption coefficient of the
                                                                        !  transparent window
                        aco2(i,j,k,jco2) = slw_a_func(T(i,j,k),0._dp,&       !Weighting coefficient of the
                           xco2(i,j,k),0._dp,0._dp,0._dp,T(i,j,k),&     !  transparent window
                           cj_sup_loc(i,j,k,jco2))
                      else
                         cj_loc = inverse_albdf(Fj_ref(jco2),T(i,j,k),& !Get the local cross-section by solving 
                            0._dp,xco2(i,j,k),0._dp,0._dp,0._dp,T_ref)  !  the implicit equation using the
                                                                        !  the previously defined Fs
                        kappaco2(i,j,k,jco2) = slw_kappa_func(T(i,j,k),&     !Absorption coefficient of the 
                           p(i,j,k),xco2(i,j,k),cj_loc)                 !  gray gas
                        aco2(i,j,k,jco2) = slw_a_func(T(i,j,k),0._dp,&       !Weighting coefficient of the
                           xco2(i,j,k),0._dp,0._dp,0._dp,T(i,j,k),&     !  gray gas
                           cj_sup_loc(i,j,k,jco2),cj_sup_loc(i,j,k,jco2-1))
                     endif
         
                  case('locally_correlated')
                     !+++++++
                     !LC-SLW
                     !+++++++
                     Fj_sup_loc(i,j,k,jco2) = get_albdf(T_ref,&         !Local ALBDF evaluated using
                        0._dp,xco2_ref,0._dp,0._dp,0._dp,T(i,j,k),&     !  reference state
                        cj_sup_ref(jco2))
                     if (transparent_window) then
                        kappaco2(i,j,k,jco2) = 0._dp                         !Absorption coefficient of the
                                                                        !  transparent window
                        aco2(i,j,k,jco2) = Fj_sup_loc(i,j,k,jco2)            !Weighting coefficient of the
                                                                        !  transparent window
                     else
                        cj_sup_loc(i,j,k,jco2) = inverse_albdf(&        !Get the local supplementar cross-section
                           Fj_sup_loc(i,j,k,jco2),T(i,j,k),0._dp,&      !  by solving the implicit equation
                           xco2(i,j,k),0._dp,0._dp,0._dp,T(i,j,k))
                           
                        cj_loc = get_slw_single_cj(&                    !Local cross-section determined        
                           cj_sup_loc(i,j,k,jco2-1),&                   !  from the supplementar ones
                           cj_sup_loc(i,j,k,jco2))
                        kappaco2(i,j,k,jco2) = slw_kappa_func(&              !Absorption coefficient of the 
                           T(i,j,k),p(i,j,k),xco2(i,j,k),cj_loc)        !  gray gas
                        aco2(i,j,k,jco2) = Fj_sup_loc(i,j,k,jco2) - &        !Weighting coefficient of the gray gas
                           Fj_sup_loc(i,j,k,jco2-1)
                     endif
                  
                  case default
                     call shutdown('Problem in the specification of &
                        &slw_nonuniform_method')
                  endselect        
                  
               enddo

            enddo
         enddo
      enddo
      
      !--------------------------
      !Solve the radiation field
      !--------------------------
      do jh2o = 0,n_gases
         do jco2 = 0,n_gases
            do i=0,nx+1
               do j=0,ny+1
                  do k=0,nz+1 
                     kappa(i,j,k) = kappah2o(i,j,k,jh2o) + kappaco2(i,j,k,jco2)
                     aIb(i,j,k) = ah2o(i,j,k,jh2o)*aco2(i,j,k,jco2)*Ib(i,j,k)
                     if (i.eq.nx/2) write(*,*) jh2o,jco2,kappa(i,j,k),&
                        ah2o(i,j,k,jh2o)*aco2(i,j,k,jco2)
                  enddo
               enddo
            enddo
            
            !Starting the processing time counter
            call cpu_time(start_proc_time)
   
            !Solving the radiation field for gas jgas
            call fvm_solution(kappa,kappa,kappa_sct,aIb,phi,x_eps,y_eps,&
                              z_eps,flux_j,source_j,.true.)
                           
            !Updating the heat flux and heat source
            flux = flux + flux_j
            source = source + source_j
         
            !Stopping the processing time counter and adding to the total
            call cpu_time(end_proc_time)
            proc_time = proc_time + end_proc_time - start_proc_time
         
         enddo
      enddo
      
      !------------------
      !Deallocate arrays
      !------------------
      
		!---------------------
      !Total computing time
      !---------------------
      call cpu_time(end_total_time)
      total_time = end_total_time - start_total_time

   endsubroutine slw_multint_black_solution2



endmodule slw_routines
