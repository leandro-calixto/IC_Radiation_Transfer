!##############################################
!Routines to solve the radiative heat transfer 
!via the SLW model for different scenarios
!##############################################
module slw_routines

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp,small
   use constants, only: fourpi,invpi,sigma
   use slw_parameters, only: slw_pref,slw_Tref,slw_xsref
   use slw_functions
   implicit none
   
contains

   !====================================================================
   !Subroutine for the solution of the RTE along a line-of-sight for a
   !non-scattering medium. The code solves the RTE for a single gray
   !gas only, and outputs the corresponding intensity
   !====================================================================
!   subroutine slw_sg_los_solution(x,T,p,xs,Fsup0,Fsup,F_ref,&
!      I0,Irad,emission_correction,absorption_correction,&
!      processing_time,total_time,absorption,emission,kappaPlanck,&
!      Planck,points,solve_RTE,quick_RTE)
   
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: assert,CheckMemAlloc,get_file_unit,&
!                                get_wall_time,shutdown
!      use constants, only: sigma,invpi
!      use exact_solution, only: exact_los_nonsct
!      use global_parameters, only: number_of_species
!      implicit none
!      integer,intent(in),optional :: points
!      integer :: i,ierr,n,ns,nx
!      real(dp),intent(in) :: Fsup,Fsup0,F_ref,I0,p(:),T(:),x(:),xs(:,:)
!      real(dp),intent(in),optional :: absorption_correction(:),&
!                                      emission_correction(:)
!      real(dp),intent(out),optional :: absorption(:),emission(:),&
!         kappaPlanck(:),Planck(:),processing_time,total_time
!      real(dp),parameter :: sigrpi = invpi*sigma
!      real(dp) :: a,csup,csup0,cloc,xmix
!      real(dp) :: end_proc_time,start_proc_time,end_total_time,&
!         start_total_time
!      real(dp),intent(out) :: Irad(:)
!      real(dp),allocatable,dimension(:) :: acor,ecor,kappa,Ib
!      logical,optional :: solve_RTE,quick_RTE
!      logical :: compute_absorption,compute_emission,compute_kplanck,&
!         compute_planck,compute_radiation,simplified

!      !-----------------------------------------------------------------
!      !Preparatory procedures
!      !-----------------------------------------------------------------
!      !Starting the total computing time counter
!      start_total_time = get_wall_time()
      
!      !Check if all arrays have the same size
!      call assert(size(x).eq.size(T))
!      call assert(size(T).eq.size(p))
!      call assert(size(p).eq.size(xs,1))
      
!      !Setting default input parameters
!      nx = size(x); if (present(points)) nx = points      
!      compute_radiation = .true.; if (present(solve_RTE)) &
!                                          compute_radiation = solve_RTE
!      simplified = .false.; if (present(quick_RTE)) &
!                                                   simplified = quick_RTE
      
!      !Set flags for optional output arguments
!      compute_absorption = .false.; if (present(absorption)) &
!                                             compute_absorption = .true.
!      compute_emission = .false.; if (present(emission)) &
!                                             compute_emission = .true.
!      compute_kplanck = .false.; if (present(kappaPlanck)) &
!                                             compute_kplanck = .true.
!      compute_planck = .false.; if (present(Planck)) &
!                                             compute_planck = .true.

!      !Surrogate names
!      ns = number_of_species

!      !-----------------------------------------------------------------
!      !Allocating arrays
!      !-----------------------------------------------------------------
!      allocate(acor(nx),stat=ierr)
!      call CheckMemAlloc('acor',ierr)
!      allocate(ecor(nx),stat=ierr)
!      call CheckMemAlloc('ecor',ierr)
!      allocate(kappa(nx),stat=ierr)
!      call CheckMemAlloc('kappa',ierr)
!      allocate(Ib(nx),stat=ierr)
!      call CheckMemAlloc('Ib',ierr)

!      !-----------------------------------------------------------------
!      !Set up emission and absorption correction factors
!      !-----------------------------------------------------------------
!      ecor = 1._dp; acor = 1._dp
!      if (present(emission_correction)) ecor = emission_correction
!      if (present(absorption_correction)) acor = absorption_correction

!      !-----------------------------------------------------------------
!      !Precomputing local radiative properties
!      !-----------------------------------------------------------------
!      prop_loop: do i=1,nx     
!         !Mole fraction of the mixture of participating species
!         xmix = 0._dp
!         do n=1,ns
!            if ((albdf_nx(n).gt.0).and.(albdf_file(n).ne.'null')) &
!               xmix = xmix + xs(i,n)
!         enddo
      
!         !Get the local supplementar cross-sections by solving the 
!         !implicit equation using the previously divided supplementar Fs
!         csup0 = inverse_albdf(Fsup0,T(i),xs(i,:),SLW_Tref)
!         csup = inverse_albdf(Fsup,T(i),xs(i,:),SLW_Tref)

!         !Get the local cross-section by solving the implicit equation 
!         !using the the previously defined Fs
!         cloc = inverse_albdf(F_ref,T(i),xs(i,:),SLW_Tref)
         
!         !Absorption coefficient of the gray gas
!         kappa(i) = acor(i)*slw_kappa_func(T(i),p(i),xmix,cloc)        

!         !Weighting coefficient of the gray gas
!         a = slw_a_func(T(i),xs(i,:),T(i),csup,csup0)

!         !RTE source term
!         Ib(i) = (ecor(i)/(acor(i)+small))*a*sigrpi*(T(i)**4._dp)
      
!      enddo prop_loop

!      !-----------------------------------------------------------------
!      !Solving radiation intensities
!      !-----------------------------------------------------------------
!      !Starting the processing time counter
!      start_proc_time = get_wall_time()

!      !Solve the RTE
!      if (compute_radiation) &
!         call exact_los_nonsct(x,kappa,Ib,I0,nx,Irad,simplified)

!      !Stopping the processing time counter
!      end_proc_time = get_wall_time()

!      !-----------------------------------------------------------------
!      !Final loop over all grid cells to finish up
!      !the calculation of the the optional properties
!      !-----------------------------------------------------------------
!      do i=1,nx
!         if (compute_emission)      emission(i) = kappa(i)*Ib(i)        !Emission                              
!         if (compute_absorption)    absorption(i) = kappa(i)*Irad(i)    !Absorption
!         if (compute_kplanck)       kappaPlanck(i) = kappa(i)           !Planck-mean absorption coefficient
!         if (compute_planck)        Planck(i) = Ib(i)                   !Planck function
!      enddo

!      !-----------------------------------------------------------------
!      !Deallocate arrays
!      !-----------------------------------------------------------------
!      deallocate(acor,ecor)
!      deallocate(Ib,kappa)

!      !-----------------------------------------------------------------
!      !Total times
!      !-----------------------------------------------------------------
!      !Processing
!      if (present(processing_time)) &                                   
!         processing_time = end_proc_time - start_proc_time
      
!      !Computing (total)
!      end_total_time = get_wall_time()
!      if (present(total_time)) &
!         total_time = end_total_time - start_total_time
      
!   endsubroutine slw_sg_los_solution

   !====================================================================
   !Subroutine for the exact solution of the RTE for a one-dimensional
   !cylindrical, non-scattering medium
   !====================================================================
   subroutine slw_los_solution(x,T,p,xs,I0,ngas,Irad,processing_time,&
      total_time,absorption,emission,kappaPlanck,points,quick_RTE)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckMemAlloc,dprint,&
         get_file_unit,get_wall_time,print_to_prompt,shutdown
      use exact_solution, only: exact_los_nonsct
      use global_parameters, only: debug_mode,number_of_species
      use omp_lib 
      use physical_functions, only: bb_emission_frac
      use slw_parameters, only: slw_bound_Ib,slw_Ib_lbound,slw_Ib_ubound
      implicit none
      character(10) :: gas_name
      character(200) :: msg
      integer,intent(in) :: ngas
      integer :: nx,ns
      integer :: i,jgas,n
      integer :: ierr,uout
      integer,intent(in),optional :: points
      logical,intent(in),optional :: quick_RTE
      logical :: compute_absorption,compute_emission,compute_planck,&
                 transparent_window,simplified
      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
      real(dp),intent(out) :: Irad(:)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                absorption(:),emission(:),kappaPlanck(:)
      real(dp) :: a,F,kp_loc,kp_ref
      real(dp) :: end_proc_time,end_total_time,start_proc_time,&
                  start_total_time
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref,&
                                        cj_sup_loc,F0,Fj_ref,Fj_sup_ref
      real(dp),allocatable,dimension(:) :: Ib,kI,kIb,u_loc
      real(dp),allocatable,dimension(:,:) :: aIb,Irad_j,kappa

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Check if all arrays have the same size
      call assert(size(x).eq.size(T))
      call assert(size(T).eq.size(p))
      call assert(size(p).eq.size(xs,2))
      
      !Set flags for optional output arguments
      simplified = .false.
      if (present(quick_RTE)) simplified = quick_RTE
      compute_absorption=.false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission=.false.
      if (present(emission))    compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: Surrogate names')
      nx = size(x); if (present(points)) nx = points                    !Number of grid points
      ns = number_of_species                                            !Number of participating species
      
      !-----------------------------------------------------------------
      !Prepare for dumping properties, if requested
      !-----------------------------------------------------------------
      if (slw_print_properties) then
         call dprint('slw_los_solution: Prepare for dumping properties')
         
         !Prepare output unit
         if (trim(slw_print_file).eq.'null') &                          !Check if a file name was given                          
            call shutdown('slw_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(slw_print_file),&                     !Open file
              form='formatted',action='write')
         write(uout,'(a,(500(i3,:,",")))') '#0,',(i,i=1,3*(ngas+1))     !Write first line
         write(uout,'(3(a,:,","),100(i3,:,","))') '#j =',&              !Write header
                                            (jgas,jgas,jgas,jgas=0,ngas)
         write(uout,'(500(a,:,","))') '#x',&                            !Write subheader
                                    ('kappa_j','a_j','I_j',jgas=0,ngas)
         write(uout,'(500(a,:,","))') '#[m]',&                          !Write subsubheader
                                    ('[1/m]','[-]','[W/m3]',jgas=0,ngas)
      endif
      
      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('slw_los_solution: &
                                          &allocating arrays')
      allocate(aIb(nx,0:ngas),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(cj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_ref',ierr)
      allocate(cj_sup_loc(nx),stat=ierr)
      call CheckMemAlloc('cj_sup_loc',ierr)
      allocate(cj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_sup_ref',ierr)
      allocate(F0(nx),stat=ierr)
      call CheckMemAlloc('F0',ierr)
      allocate(Fj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_ref',ierr)
      allocate(Fj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_sup_ref',ierr)
      allocate(Ib(nx),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(Irad_j(nx,0:ngas),stat=ierr)
      call CheckMemAlloc('Irad_j',ierr)
      allocate(kappa(nx,0:ngas),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kI(nx),stat=ierr)
      call CheckMemAlloc('kI',ierr)
      allocate(kIb(nx),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(u_loc(nx),stat=ierr)
      call CheckMemAlloc('u_loc',ierr)

      !-----------------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !-----------------------------------------------------------------
      !Blackbody intensity
      call dprint('slw_los_solution: Compute blackbody intensity')
      do i=1,nx
         if (slw_bound_Ib) then
            F = bb_emission_frac(slw_Ib_lbound,T(i),.true.) -&
                bb_emission_frac(slw_Ib_ubound,T(i),.true.)
         else
            F = 1._dp
         endif
            Ib(i) = F*sigma*invpi*(T(i)**4._dp)
      enddo
      
      !-----------------------------------------------------------------
      !Defining absorption cross-sections and 
      !supplementar absorption cross-sections
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: Define absorption cross-sections')
      if (trim(slw_nonuniform_method).eq.'rank_correlated') then
         call get_slw_Fj(ngas,Fj_ref(1:ngas),Fj_sup_ref(1:ngas))
         Fj_sup_ref(0) = slw_Fmin  

      else
         !For other SLWs, divide the reference Cj
         call get_slw_cj(slw_cmin,slw_cmax,ngas,cj_ref(1:ngas),&
                         cj_sup_ref(1:ngas))
         cj_sup_ref(0) = slw_cmin
      endif

      !-----------------------------------------------------------------
      !Load the ALBDF database
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: Load the ALBDF')
      call load_albdf

      !-----------------------------------------------------------------
      !Define the reference state, if necessary
      !-----------------------------------------------------------------
      if (trim(slw_nonuniform_method).ne.'uniform') then
         call dprint('slw_los_solution: Define the reference state')
         slw_Tref = maxval(T)
         slw_pref = maxval(p)
         do n=1,ns
           ! slw_xsref(n) = maxval(xs(n,:,:,:))
         enddo
      endif
      
      !-----------------------------------------------------------------
      !Compute scaling parameter, if necessary
      !-----------------------------------------------------------------
      if (trim(slw_nonuniform_method).eq.'scaled') then
         call dprint('slw_los_solution: Compute scaling parameter')
         do i=1,nx
            !Planck-mean absorption coefficients at
            !local and reference conditions
            kp_loc = get_slw_kp(T(i),p(i),xs(:,i),T(i),ngas)
            kp_ref = get_slw_kp(slw_Tref,slw_pref,slw_xsref,T(i),ngas)

            !Scaling parameter
            u_loc(i) = kp_loc/(kp_ref+small)
         enddo
      endif

      !-----------------------------------------------------------------
      !Compute scaling parameter, if necessary
      !-----------------------------------------------------------------
      u_loc = 1._dp
      if (trim(slw_nonuniform_method).eq.'scaled') then
         call dprint('slw_los_solution: Compute scaling parameter')
         do i=1,nx
            kp_loc = get_slw_kp(T(i),p(i),xs(:,i),T(i),ngas)            !Planck-mean absorption coefficient
                                                                        !  at the local thermodynamic state                  
            kp_ref = get_slw_kp(slw_Tref,slw_pref,slw_xsref,T(i),ngas)  !Planck-mean absorption coefficient
                                                                        !  at the reference thermodynamic state
            u_loc(i) = kp_loc/(kp_ref+small)                            !Scaling parameter
         enddo
      endif

      !-----------------------------------------------------------------
      !Compute radiative properties at each grid point and gray gas
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: Compute local properties')
      kIb = 0._dp                                                       !Zero out array with the total RTE emission term
      do jgas=0,ngas
         transparent_window = .false.                                   !Set transparent window flag
         if (jgas.eq.0) transparent_window = .true.                     ! for jgas = 0
         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP& PRIVATE(a) 
         do i=1,nx
            call compute_slw_gas_parameters(kappa(i,jgas),a,&           !Subroutine responsible for computing
               cj_sup_ref(jgas),cj_sup_loc(i),Fj_ref(jgas),&            !  the absorption and weighting coefficients
               Fj_sup_ref(jgas),F0(i),T(i),xs(:,i),transparent_window=& !  of the gray gas
               transparent_window,scaling_coeff=u_loc(i))
            aIb(i,jgas) = a*Ib(i)                                       !RTE emission term
         enddo
         !$OMP END PARALLEL DO
      enddo

      !-----------------------------------------------------------------
      !Main radiative transfer solution loop
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: &
                  &Main radiative transfer solution loop started')
      start_proc_time = get_wall_time()                                 !Start the processing time counter

      !$OMP PARALLEL DO DEFAULT(SHARED)
      gas_loop: do jgas=0,ngas
         !Longer debug information
         write(gas_name,'(i5)') jgas
         msg = 'slw_los_solution: gas_loop for gas '//&
            trim(adjustl(gas_name))//' started'
         call dprint(msg)

         !Solve the RTE for gas jgas
         call exact_los_nonsct(x,kappa(:,jgas),aIb(:,jgas),I0,nx,&
                               Irad_j(:,jgas),simplified)      
      enddo gas_loop
      !$OMP END PARALLEL DO
      end_proc_time = get_wall_time()                                   !Stop the processing time counter

      !-----------------------------------------------------------------
      !Dump gas properties, if necessary
      !-----------------------------------------------------------------
      if (slw_print_properties) then
         call dprint('slw_los_solution: Dump gas properties')
         do i=1,nx
            write(uout,'(1000(e26.15e3,:,","))') x(i),(kappa(i,jgas),&
               aIb(i,jgas)/Ib(i),Irad_j(i,jgas),jgas=0,ngas)
         enddo
         close(uout)                                                    !Close output unit
      endif
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of total quantities
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: &
                               &Finish calculation of total quantities')
      Irad = 0._dp; kI = 0._dp; kIb = 0._dp                             !Initialize sum arrays
      do i=1,nx
         do jgas=0,ngas
            Irad(i) = Irad(i) + Irad_j(i,jgas)
            kI(i) = kI(i) + kappa(i,jgas)*Irad_j(i,jgas)
            kIb(i) = kIb(i) + kappa(i,jgas)*aIb(i,jgas)
         enddo
         
         !Optional quantities
         if (compute_emission)   emission(i) = kIb(i)                   !Total emission
         if (compute_absorption) absorption(i) = kI(i)                  !Total absorption
         if (compute_planck)     kappaPlanck(i) = kIb(i)/Ib(i)          !Planck-mean absorption coefficient
      enddo
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: Deallocate arrays')
      deallocate(aIb,Ib,Irad_j,kappa,kI,kIb)
      deallocate(cj_ref,cj_sup_loc,cj_sup_ref)
      deallocate(Fj_ref,Fj_sup_ref)
      deallocate(u_loc)
      
      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('slw_los_solution: Compute computational times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()                                  !Stop the total time counter
      if (present(total_time)) &                                        !Total
         total_time = end_total_time - start_total_time
      
   endsubroutine slw_los_solution

   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !SLW- model for a non-scattering medium bounded by black surfaces
   !====================================================================
   subroutine slw_black_solution(flux,source,emission_correction,&
      absorption_correction,processing_time,total_time,absorption,&
      emission,kappaPlanck,solve_RTE,input_ngas,output_ngas)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,print_to_prompt,shutdown
      use constants, only: sigrpi
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: number_of_species 
      use mesh
      use omp_lib
      use physical_functions, only: bb_emission_frac
      use slw_parameters, only: slw_bound_Ib,slw_Ib_lbound,slw_Ib_ubound
      implicit none
      character(10) :: gas_name
      character(200) :: msg
      integer,optional,intent(in) :: input_ngas
      integer,optional,intent(out) :: output_ngas
      integer :: ncj,ndata,ngas,nx,ny,nz,nl,nspc
      integer :: i,idata,idsp,iispc,ispc,j,jgas,jjgas,k
      integer :: gas_counter,gas_lbound,ierr!,uout
      integer,allocatable,dimension(:) :: &
         dominant_spc,original_gas_index
      integer,allocatable,dimension(:,:) :: data_array,spc_gas_index
      real(dp),intent(in),optional :: &
         absorption_correction(:,:,:,:),&
         emission_correction(:,:,:,:)
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints)
      real(dp) :: a,aa,cj_sup_loc,F,F0,kappa_comp,kappa_mix,kp_loc,kp_ref
      real(dp) :: p_avg,T_avg,theta
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref,&
                                           Fj_ref,Fj_sup_ref
      real(dp),allocatable,dimension(:) :: xs_avg,xs_cur
      real(dp),allocatable,dimension(:,:) :: a_spc,kappa_spc,phi
      real(dp),allocatable,dimension(:,:) :: &
         cj_ref_spc,cj_sup_ref_spc,Fj_ref_spc,Fj_sup_ref_spc
      real(dp),allocatable,dimension(:,:,:) :: Ib,kappa_sct,kIb,u_loc,&
         source_j,source_omp,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: aIb,kappa
      real(dp),allocatable,dimension(:,:,:,:) :: flux_j,flux_omp
      real(dp),allocatable,dimension(:,:,:,:) :: acor,ecor
      logical,optional :: solve_RTE
      logical :: advance_index,improved_multiple_integration,&
                 multiple_integration,transparent_window
      logical :: compute_absorption,compute_emission,compute_planck,&
                 compute_radiation
      logical,allocatable,dimension(:) :: skip_gas

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Set flags for optional output arguments
      compute_radiation = .true.
      if (present(solve_RTE)) compute_radiation = solve_RTE
      compute_absorption=.false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission=.false.
      if (present(emission))    compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.
      
      !Surrogate names
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nspc = number_of_species                                          !Number of participating species
      ndata = nx*ny*nz                                                  !Size of the first rank of data_array
      theta = slw_superposition_theta
      
      !-----------------------------------------------------------------
      !Define number of gases 
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Define number of gases')
      if (trim(slw_mixture_method).eq.'multiple_integration') then
         !Set flags and general parameters
         multiple_integration = .true.
         improved_multiple_integration = .true.
         if (theta.ge.1.e30_dp) improved_multiple_integration = .false.
         gas_lbound = 1
         
         !For the multiple integration approach, determine the (max) 
         !total number of gray gases from the provided number of gray
         !gases per species
         ngas = 1
         do ispc=1,nspc
            ngas = ngas*(1 + slw_ngas_per_species(ispc))
         enddo
      
      else      
         !For all other mixing approaches, use 
         !the given number of gray gases
         ngas = slw_ngas
         if (present(input_ngas)) ngas = input_ngas
         
         !Set flags and general parameters
         multiple_integration = .false.
         gas_lbound = 0
      endif
      
!      !-----------------------------------------------------------------
!      !Prepare for dumping properties, if requested
!      !-----------------------------------------------------------------
!      if (slw_print_properties) then
!         call dprint('slw_black_solution: &
!                     &Prepare for dumping properties')
         
!         !Prepare output unit
!         if (trim(slw_print_file).eq.'null') &                          !Check if a file name was given                          
!            call shutdown('slw_print_file undefined')
!         uout = get_file_unit()                                         !Get unit
!         open(unit=uout,file=trim(slw_print_file),&                     !Open file
!              form='formatted',action='write')
!         write(uout,'(a,(100(i3,:,",")))') '#0,',(i,i=1,2+2*(n_gases+1))!Write first line
!         write(uout,'(3(a,:,","),100(i3,:,","))') '#',' ','j =',&       !Write header
!                                             (jgas,jgas,jgas=0,n_gases)
!         write(uout,'(100(a,:,","))') '#x','y','z',&                    !Write subheader
!                                    ('kappa_j','a_j',jgas=0,n_gases)
!         write(uout,'(100(a,:,","))') '#[m]','[m]','[m]',&              !Write subsubheader
!                                    ('[1/m]','[-]',jgas=0,n_gases)
!      endif
      
      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Allocating arrays')
      allocate(aIb(nx,ny,nz,gas_lbound:ngas),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(acor(nx,ny,nz,gas_lbound:ngas),stat=ierr)
      call CheckMemAlloc('acor',ierr)
      allocate(data_array(ndata,3),stat=ierr)
      call CheckMemAlloc('data_array',ierr)
      allocate(ecor(nx,ny,nz,gas_lbound:ngas),stat=ierr)
      call CheckMemAlloc('ecor',ierr)      
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,gas_lbound:ngas),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(u_loc(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('u_loc',ierr)
      allocate(source_j(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(source_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_omp',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)
      if (multiple_integration) then
         allocate(a_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('a_spc',ierr)
         allocate(cj_ref_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('cj_ref_spc',ierr)
         allocate(cj_sup_ref_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('cj_sup_ref_spc',ierr)
         allocate(dominant_spc(1:ngas),stat=ierr)
         call CheckMemAlloc('dominant_spc',ierr)
         allocate(Fj_ref_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('Fj_ref_spc',ierr)
         allocate(Fj_sup_ref_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('Fj_sup_ref_spc',ierr)
         allocate(kappa_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('kappa_spc',ierr)
         allocate(original_gas_index(1:ngas),stat=ierr)
         call CheckMemAlloc('original_gas_index',ierr)
         allocate(skip_gas(1:ngas),stat=ierr)
         call CheckMemAlloc('skip_gas',ierr)
         allocate(spc_gas_index(1:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('spc_gas_index',ierr)
         allocate(xs_avg(1:nspc),stat=ierr)
         call CheckMemAlloc('xs_avg',ierr)
         allocate(xs_cur(1:nspc),stat=ierr)
         call CheckMemAlloc('xs_cur',ierr)
      else
         allocate(cj_ref(0:ngas),stat=ierr)
         call CheckMemAlloc('cj_ref',ierr)
         allocate(cj_sup_ref(0:ngas),stat=ierr)
         call CheckMemAlloc('cj_sup_ref',ierr)
         allocate(Fj_ref(0:ngas),stat=ierr)
         call CheckMemAlloc('Fj_ref',ierr)
         allocate(Fj_sup_ref(0:ngas),stat=ierr)
         call CheckMemAlloc('Fj_sup_ref',ierr)
      endif

      !-----------------------------------------------------------------
      !Set up emission and absorption correction factors
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: &
                  &Set up emission and absorption correction factors')
      ecor = 1._dp; acor = 1._dp
      if (present(emission_correction))   ecor = emission_correction
      if (present(absorption_correction)) acor = absorption_correction

      !-----------------------------------------------------------------
      !Mount data_array
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Mount data_array')
      idata = 0
      do i=1,nx
         do j=1,ny
            do k=1,nz
               idata = idata + 1
               data_array(idata,1:3) = (/ i,j,k /)
            enddo
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !-----------------------------------------------------------------
      !Blackbody intensity
      call dprint('slw_black_solution: Compute blackbody intensity')
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (slw_bound_Ib) then
                  F = bb_emission_frac(slw_Ib_lbound,T(i,j,k),.true.) -&
                      bb_emission_frac(slw_Ib_ubound,T(i,j,k),.true.)
               else
                  F = 1._dp
               endif
               Ib(i,j,k) = F*sigrpi*(T(i,j,k)**4._dp)
            enddo
         enddo
      enddo

      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !Set wall emissivities to unit for all wall cells
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp

      !-----------------------------------------------------------------
      !Load the ALBDF database
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Load the ALBDF')
!      if (trim(cell_albdf).ne.'null') then        
!         !Allocate point-by-point ALBDF and Cj arrays
!         if (allocated(ijk_albdf)) deallocate(ijk_albdf)
!         allocate(ijk_albdf(5*n_gases,nx,ny,nz),stat=ierr)
!         call CheckMemAlloc('ijk_albdf',ierr)
!         if (allocated(ijk_albdf_cj)) deallocate(ijk_albdf)
!         allocate(ijk_albdf_cj(5*n_gases),stat=ierr)
!         call CheckMemAlloc('ijk_albdf_cj',ierr)

!         !Mount Cj array (same for all grid points)
         
!         if (cell_albdf.eq.'from-nbck') then
!            do i=1,nx
!               do j=1,ny
!                  do k=1,nz
!                     call albdf_from_nbck(ijk_albdf_cj,&
!                        ijk_albdf(:,i,j,k),T(i,j,k),p(i,j,k),xgas,Tsource)
!                  enddo
!               enddo
!            enddo
         
!         endif 
      
!      else
         call load_albdf
!      endif

      !-----------------------------------------------------------------
      !Define the reference state, if necessary
      !-----------------------------------------------------------------
      if (trim(slw_nonuniform_method).ne.'uniform') then
         call dprint('slw_black_solution: Define the reference state')
         slw_Tref = maxval(T)
         slw_pref = maxval(p)
         do ispc=1,nspc
            slw_xsref(ispc) = maxval(xs(ispc,:,:,:))
         enddo
      endif

      !-----------------------------------------------------------------
      !Compute scaling parameter, if necessary
      !-----------------------------------------------------------------
      u_loc = 1._dp
      if (trim(slw_nonuniform_method).eq.'scaled') then
         call dprint('slw_black_solution: Compute scaling parameter')
         do idata=1,ndata
            i = data_array(idata,1)                                     !Define i,
            j = data_array(idata,2)                                     !  j and
            k = data_array(idata,3)                                     !  k indexes
            kp_loc = get_slw_kp(T(i,j,k),p(i,j,k),xs(:,i,j,k),&         !Planck-mean absorption coefficient
                                T(i,j,k),ngas)                          !  at the local thermodynamic state                  
            kp_ref = get_slw_kp(slw_Tref,slw_pref,slw_xsref,&           !Planck-mean absorption coefficient
                                T(i,j,k),ngas)                          !  at the reference thermodynamic state
            u_loc(i,j,k) = kp_loc/(kp_ref+small)                        !Scaling parameter
         enddo
      endif

      !-----------------------------------------------------------------
      !Prepare superposition arrays
      !-----------------------------------------------------------------
      if (multiple_integration) then
          call dprint('slw_black_solution: &
                      &Prepare superposition arrays')
      
         !Mount the spc_gas_index array 
         !(array storing the gray gas index of each individual
         !species ispc that corresponds to the overal gray gas jgas)
         spc_gas_index = 0                                              !For jgas = 1, all species' gray gases are
                                                                        !  transparent windows
         do jgas=2,ngas
            do ispc=1,nspc        
               if (ispc.eq.1) then                                         !The first species is treated normally,
                  spc_gas_index(jgas,ispc) = &                             !  with each new jgas corresponding to a
                     mod((spc_gas_index(jgas-1,ispc)+1),&                  !  new j for the species, and the counter
                        (slw_ngas_per_species(ispc)+1))                           !  resetting once j = J
               else                                                        !For all subsequent species, a new step in
                  spc_gas_index(jgas,ispc) = spc_gas_index(jgas-1,ispc)    !  j is only taken (advance_index = .true)
                  advance_index = .true.                                   !  if j = J @ jgas-1 for all previous species
                  do iispc=1,ispc-1                                        !This loops checks if j = J @ jgas-1 for all
                     if (spc_gas_index(jgas-1,iispc).ne.&                  !  previous species. If that is the case, then
                        slw_ngas_per_species(iispc)) advance_index = .false.      !  advance_index = .true.
                  enddo                                                    !And, if that is the case, than the j stepping
                  if (advance_index) spc_gas_index(jgas,ispc) = &          !  is done, once again resetting if j = J for
                     mod((spc_gas_index(jgas,ispc)+1),&                    !  that species
                        (slw_ngas_per_species(ispc)+1)) 
               endif
            enddo 
         enddo

         !Define the Fj and Cj divisions for each species separately
         do ispc=1,nspc
            ncj = slw_ngas_per_species(ispc)
            if (trim(slw_nonuniform_method).eq.'rank_correlated') then
               !For the Rank-SLW, divide the reference Fj
               call get_slw_Fj(ncj,Fj_ref_spc(1:ncj,ispc),&
                  Fj_sup_ref_spc(1:ncj,ispc))
               Fj_sup_ref_spc(0,ispc) = slw_Fmin  
            else
               !For other SLWs, divide the reference Cj
               call get_slw_cj(slw_cmin,slw_cmax,ngas,&
                  cj_ref_spc(1:ncj,ispc),cj_sup_ref_spc(1:ncj,ispc))
               cj_sup_ref_spc(0,ispc) = slw_cmin
            endif
         enddo

         !Compute the average composition
         if (improved_multiple_integration) &
            call get_average_composition(T_avg,xs_avg,p_avg,'T4')

         !Mount the dominant_spc array
         !(this array informs which species is the dominant one, in the
         !framework of the improved superposition method; a value of -1
         !indicates that no species dominates)
         dominant_spc(1:ngas) = -1                                      !Initially, assume that no species dominate
         if (improved_multiple_integration) then
            !First compute kappa for each species and each species'
            !gray gas separately
            do ispc=1,nspc
               xs_cur = 0._dp; xs_cur(ispc) = xs_avg(ispc)
               do jgas=0,slw_ngas_per_species(ispc)
                  transparent_window = .false.
                  if (jgas.eq.0) transparent_window = .true.
                  call compute_slw_gas_parameters(kappa_spc(jgas,ispc),&
                     a_spc(jgas,ispc),cj_sup_ref_spc(jgas,ispc),&
                     cj_sup_loc,Fj_ref_spc(jgas,ispc),&
                     Fj_sup_ref_spc(jgas,ispc),F0,T_avg,&
                     xs_cur,transparent_window=transparent_window)!,&
!                     scaling_coeff=u_loc(i,j,k))
               enddo
            enddo
            
            do jgas=1,ngas                                              !Now loop over all mixture gray gases
               !Compute reference values of kappa
               kappa_mix = 0._dp
               do ispc=1,nspc
                  kappa_mix = kappa_mix + &
                     kappa_spc(spc_gas_index(jgas,ispc),ispc)
               enddo
               do ispc=1,nspc                                           !For each species, compute kappa of the mixture minus
                  kappa_comp = kappa_mix - &                            !  the current species and compare to kappa of the 
                     kappa_spc(spc_gas_index(jgas,ispc),ispc)           !  species, taking into account the safety factor theta
                  if (kappa_spc(spc_gas_index(jgas,ispc),ispc).gt.&     !If the comparison is successful, flag that as the
                      kappa_comp*theta) then                            !  dominant species for this jgas
                     dominant_spc(jgas) = ispc                             
                     exit                                                  
                  endif
               enddo
            enddo
         endif

         !Mount the skip_gas array
         !(this array informs if a jgas gas should be skipped when 
         !solving the RTE in the improved superposition method to 
         !prevent repeated calculations. This is done by sweeping 
         !through all gray gases larger than the current one and 
         !checking if the main condition of the improved method is not 
         !verified for any other combination of species' gray gases that
         !include the current gas for the dominating species species
         do jgas=1,ngas
            idsp = dominant_spc(jgas)
            if (idsp.le.0) cycle
            skip_gas(jgas) = .false.                                    !In principle, do not skip this gas
            gas_check_loop: do jjgas=jgas+1,ngas                        !Loop over all remaining gray gases, but consider only
               if (spc_gas_index(jjgas,idsp).ne.&                       !  jgas correspondig to the same gray gas index as that
                  spc_gas_index(jgas,idsp)) cycle gas_check_loop        !  of the current gray gas index of the dominant species
               if (dominant_spc(jjgas).eq.idsp) then                    !If, for any of the remaining gases, the current dominant
                  skip_gas(jgas) = .true.                               !  species dominates once again, we should flag jgas to
                  exit gas_check_loop                                   !  be skipped
               endif
            enddo gas_check_loop
         enddo      
      endif

      !-----------------------------------------------------------------
      !Compute radiative properties at each grid point and gray gas
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Compute local properties',3)
      kIb = 0._dp                                                       !Zero out array with the total RTE emission term
      if (multiple_integration) then
         do idata=1,ndata
            !Define i, j and k indexes
            i = data_array(idata,1)
            j = data_array(idata,2)
            k = data_array(idata,3)
            
            !First compute kappa for each species and each species'
            !gray gas separately
            do ispc=1,nspc
               xs_cur = 0._dp; xs_cur(ispc) = xs(ispc,i,j,k)
               do jgas=0,slw_ngas_per_species(ispc)
                  transparent_window = .false.
                  if (jgas.eq.0) transparent_window = .true.
                  call compute_slw_gas_parameters(kappa_spc(jgas,ispc),&
                     a_spc(jgas,ispc),cj_sup_ref_spc(jgas,ispc),&
                     cj_sup_loc,Fj_ref_spc(jgas,ispc),&
                     Fj_sup_ref_spc(jgas,ispc),F0,T(i,j,k),&
                     xs_cur,transparent_window=transparent_window,&
                     scaling_coeff=u_loc(i,j,k))
               enddo
            enddo
            
            gas_counter = 0
            do jgas=1,ngas
               idsp = dominant_spc(jgas)
               if (idsp.lt.0) then                                      !Standard superposition
                  gas_counter = gas_counter + 1
                  original_gas_index(gas_counter) = jgas
                  kappa(i,j,k,gas_counter) = 0._dp; a = 1._dp
                  do ispc=1,nspc
                     kappa(i,j,k,gas_counter) = &
                        kappa(i,j,k,gas_counter) + &
                           kappa_spc(spc_gas_index(jgas,ispc),ispc)
                     a = a*a_spc(spc_gas_index(jgas,ispc),ispc)
                  enddo
               else                                                     !Improved superposition
                  if (skip_gas(jgas)) cycle                             !Skip the jgas if the jgas is to be skipped
                  gas_counter = gas_counter + 1
                  original_gas_index(gas_counter) = jgas
                  kappa(i,j,k,gas_counter) = &                          !kappa for the gas is equal to kappa of the
                     kappa_spc(spc_gas_index(jgas,idsp),idsp)           !   dominant species
                  a = 1._dp
                  do ispc=1,nspc                                        !Initialize a as the value of a
                     a = a*a_spc(spc_gas_index(jgas,ispc),ispc)
                  enddo
                  do jjgas=1,jgas-1                                     !To compute a for the gas, loop over all previous jgas
                     if (dominant_spc(jjgas).ne.idsp) cycle             !  skipping jgas for which the dominant species is not 
                                                                        !  equal to the current dominant species, and
                     if (spc_gas_index(jjgas,idsp).ne.&                 !  jgas for which the dominant species' gray gas is 
                         spc_gas_index(jgas,idsp)) cycle                !  not equal to the current dominant species' gray gas
                     aa = 1._dp
                     do ispc=1,nspc
                        aa = aa*a_spc(spc_gas_index(jjgas,ispc),ispc)
                     enddo
                     a = a + aa
                  enddo
               endif  
               aIb(i,j,k,gas_counter) = a*Ib(i,j,k)                     !RTE emission term
               kIb(i,j,k) = kIb(i,j,k) + &                              !Total RTE emission term
                     kappa(i,j,k,gas_counter)*aIb(i,j,k,gas_counter)
            enddo
         enddo
      
         !Update the actual total number of gray gases
         ngas = gas_counter
         if (present(output_ngas)) output_ngas = ngas
      
      else

         !For the Rank-SLW, divide the reference Fj; for all other SLWs,
         !divide the reference Cj. If we are using multiple 
         !integration, this is done for each species one by one; 
         !otherwise, this is done only once, for the mixture
         if (trim(slw_nonuniform_method).eq.'rank_correlated') then
            call get_slw_Fj(ngas,Fj_ref(1:ngas),Fj_sup_ref(1:ngas))
               Fj_sup_ref(0) = slw_Fmin  
         else
            call get_slw_cj(slw_cmin,slw_cmax,ngas,cj_ref(1:ngas),&
                            cj_sup_ref(1:ngas))
            cj_sup_ref(0) = slw_cmin         
         endif

         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP& PRIVATE(a,cj_sup_loc,F0,i,j,jgas,k,transparent_window) 
         do idata=1,ndata
            !Define i, j and k indexes
            i = data_array(idata,1)
            j = data_array(idata,2)
            k = data_array(idata,3)

            !Compute kappa and a for the mixture
            do jgas=0,ngas
               transparent_window = .false.                             !Set transparent window flag
               if (jgas.eq.0) transparent_window = .true.               ! for jgas = 0
               call compute_slw_gas_parameters(kappa(i,j,k,jgas),a,&
                  cj_sup_ref(jgas),cj_sup_loc,Fj_ref(jgas),&
                  Fj_sup_ref(jgas),F0,T(i,j,k),&
                  xs(:,i,j,k),transparent_window=transparent_window,&
                  scaling_coeff=u_loc(i,j,k))
               aIb(i,j,k,jgas) = a*Ib(i,j,k)                            !RTE emission term
               kIb(i,j,k) = kIb(i,j,k) + &                              !Total RTE emission term
                  kappa(i,j,k,jgas)*aIb(i,j,k,jgas)
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

      !-----------------------------------------------------------------
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Initialize sum arrays')
      flux = 0._dp; source = 0._dp
      start_proc_time = get_wall_time()                                 !Start the processing time counter
      
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_j,flux_omp,source_j,source_omp)
      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      flux_omp = 0._dp; source_omp = 0._dp
      
      !-----------------------------------------------------------------
      !Solve the radiation field
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Radiative transfer loop started')
      !$OMP DO
      gas_loop: do jgas=gas_lbound,ngas
         !Set up longer debug information
         write(gas_name,'(i5)') jgas
         msg = 'slw_black_solution: gas_loop for gas '//&
               trim(adjustl(gas_name))//' started'
         call dprint(msg)

         !Solve the radiation field for gas jgas
         call dprint('slw_black_solution: RTE solution',3)
         if (compute_radiation) &
            call fvm_solution(kappa(:,:,:,jgas),kappa(:,:,:,jgas),&
                              kappa_sct,aIb(:,:,:,jgas),phi,x_eps,&
                              y_eps,z_eps,flux_j,source_j,.true.)
                           
         !Updating the heat flux and heat source
         call dprint('slw_black_solution: &
                     &Update the heat flux and heat source',3)
         flux_omp = flux_omp + flux_j
         source_omp = source_omp + source_j
      
      enddo gas_loop
      !$OMP ENDDO
      
      !-----------------------------------------------------------------
      !Finish the calculation of the total radiative quantities
      !-----------------------------------------------------------------
      !$OMP CRITICAL
      call dprint('slw_black_solution: Finish the &
                  &calculation of the total radiative quantities')
      flux = flux + flux_omp
      source = source + source_omp
      !$OMP ENDCRITICAL
      !$OMP ENDPARALLEL
      
      !Stop processing time counter
      call dprint('slw_black_solution: Radiative transfer loop ended')
      end_proc_time = get_wall_time()
      
!      !-----------------------------------------------------------------
!      !Dump of gas properties, if necessary
!      !-----------------------------------------------------------------
!      if (slw_print_properties) then
!         call dprint('slw_black_solution: Dump of gas properties')
!         do i=1,nx
!            do j=1,ny
!               do k=1,nz
!                  write(uout,'(1000(e26.15e3,:,","))')&
!                     x(i),y(j),z(k),(kappa(i,j,k,jgas),&
!                        aIb(i,j,k,jgas)/Ib(i,j,k),jgas=0,ngas)
!               enddo
!            enddo
!         enddo
         
!         !Close output unit
!         close(uout)

!      endif
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: &
                            &Finish calculation of optional properties')
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (compute_emission) &
                  emission(i,j,k) = fourpi*kIb(i,j,k)                   !Total emission
               if (compute_absorption) &
                  absorption(i,j,k) = source(i,j,k) + fourpi*kIb(i,j,k) !Total absorption
               if (compute_planck) &
                  kappaPlanck(i,j,k) = kIb(i,j,k)/Ib(i,j,k)             !Planck-mean absorption coefficient
            enddo
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Deallocate arrays')
      deallocate(aIb,Ib,kIb)
      deallocate(acor,ecor)
      deallocate(data_array)
      deallocate(flux_j,flux_omp,source_j,source_omp)
      deallocate(kappa,kappa_sct,phi)
      deallocate(u_loc)
      deallocate(x_eps,y_eps,z_eps)
      if (multiple_integration) then
         deallocate(a_spc,kappa_spc)
         deallocate(cj_ref_spc,cj_sup_ref_spc)
         deallocate(Fj_ref_spc,Fj_sup_ref_spc)
         deallocate(dominant_spc,skip_gas,spc_gas_index)
         deallocate(original_gas_index)
         deallocate(xs_avg,xs_cur)
      else
         deallocate(cj_ref,cj_sup_ref)
         deallocate(Fj_ref,Fj_sup_ref)
      endif
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time
      
   endsubroutine slw_black_solution

   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !SLW- model for a non-scattering medium bounded by gray surfaces
   !====================================================================
   subroutine slw_gray_solution(flux,source,emission_correction,&
      absorption_correction,processing_time,total_time,absorption,&
      emission,kappaPlanck,solve_RTE,input_ngas,output_ngas)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,print_to_prompt,shutdown
      use constants, only: sigrpi
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: number_of_species 
      use mesh
      use omp_lib
      use physical_functions, only: bb_emission_frac
      use slw_parameters, only: slw_bound_Ib,slw_Ib_lbound,slw_Ib_ubound
      implicit none
      character(10) :: gas_name
      character(200) :: msg
      integer,optional,intent(in) :: input_ngas
      integer,optional,intent(out) :: output_ngas
      integer :: ncj,ndata,ngas,nx,ny,nz,nl,nspc
      integer :: i,idata,idsp,ieps,iispc,ispc,j,jgas,jjgas,k
      integer :: gas_counter,gas_lbound,ierr!,uout
      integer,allocatable,dimension(:) :: &
         dominant_spc,original_gas_index
      integer,allocatable,dimension(:,:) :: data_array,spc_gas_index
      real(dp),intent(in),optional :: &
         absorption_correction(:,:,:,:),&
         emission_correction(:,:,:,:)
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints)
      real(dp) :: a,aa,cj_sup_loc,F,F0,kappa_comp,kappa_mix,kp_loc,kp_ref
      real(dp) :: p_avg,T_avg,theta
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref,&
                                           Fj_ref,Fj_sup_ref
      real(dp),allocatable,dimension(:) :: xs_avg,xs_cur
      real(dp),allocatable,dimension(:,:) :: a_spc,kappa_spc,phi
      real(dp),allocatable,dimension(:,:) :: &
         cj_ref_spc,cj_sup_ref_spc,Fj_ref_spc,Fj_sup_ref_spc
      real(dp),allocatable,dimension(:,:,:) :: Ib,kappa_sct,kIb,u_loc,&
         source_j,source_omp,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: aIb,kappa
      real(dp),allocatable,dimension(:,:,:,:) :: flux_j,flux_omp
      real(dp),allocatable,dimension(:,:,:,:) :: acor,ecor
      logical,optional :: solve_RTE
      logical :: advance_index,improved_multiple_integration,&
                 multiple_integration,transparent_window
      logical :: compute_absorption,compute_emission,compute_planck,&
                 compute_radiation
      logical,allocatable,dimension(:) :: skip_gas

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Set flags for optional output arguments
      compute_radiation = .true.
      if (present(solve_RTE)) compute_radiation = solve_RTE
      compute_absorption=.false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission=.false.
      if (present(emission))    compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.
      
      !Surrogate names
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nspc = number_of_species                                          !Number of participating species
      ndata = nx*ny*nz                                                  !Size of the first rank of data_array
      theta = slw_superposition_theta
      
      !-----------------------------------------------------------------
      !Define number of gases 
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Define number of gases')
      if (trim(slw_mixture_method).eq.'multiple_integration') then
         !Set flags and general parameters
         multiple_integration = .true.
         improved_multiple_integration = .true.
         if (theta.ge.1.e30_dp) improved_multiple_integration = .false.
         gas_lbound = 1
         
         !For the multiple integration approach, determine the (max) 
         !total number of gray gases from the provided number of gray
         !gases per species
         ngas = 1
         do ispc=1,nspc
            ngas = ngas*(1 + slw_ngas_per_species(ispc))
         enddo
      
      else      
         !For all other mixing approaches, use 
         !the given number of gray gases
         ngas = slw_ngas
         if (present(input_ngas)) ngas = input_ngas
         
         !Set flags and general parameters
         multiple_integration = .false.
         gas_lbound = 0
      endif
      
!      !-----------------------------------------------------------------
!      !Prepare for dumping properties, if requested
!      !-----------------------------------------------------------------
!      if (slw_print_properties) then
!         call dprint('slw_black_solution: &
!                     &Prepare for dumping properties')
         
!         !Prepare output unit
!         if (trim(slw_print_file).eq.'null') &                          !Check if a file name was given                          
!            call shutdown('slw_print_file undefined')
!         uout = get_file_unit()                                         !Get unit
!         open(unit=uout,file=trim(slw_print_file),&                     !Open file
!              form='formatted',action='write')
!         write(uout,'(a,(100(i3,:,",")))') '#0,',(i,i=1,2+2*(n_gases+1))!Write first line
!         write(uout,'(3(a,:,","),100(i3,:,","))') '#',' ','j =',&       !Write header
!                                             (jgas,jgas,jgas=0,n_gases)
!         write(uout,'(100(a,:,","))') '#x','y','z',&                    !Write subheader
!                                    ('kappa_j','a_j',jgas=0,n_gases)
!         write(uout,'(100(a,:,","))') '#[m]','[m]','[m]',&              !Write subsubheader
!                                    ('[1/m]','[-]',jgas=0,n_gases)
!      endif
      
      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('slw_black_solution: Allocating arrays')
      allocate(aIb(nx,ny,nz,gas_lbound:ngas),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(acor(nx,ny,nz,gas_lbound:ngas),stat=ierr)
      call CheckMemAlloc('acor',ierr)
      allocate(data_array(ndata,3),stat=ierr)
      call CheckMemAlloc('data_array',ierr)
      allocate(ecor(nx,ny,nz,gas_lbound:ngas),stat=ierr)
      call CheckMemAlloc('ecor',ierr)      
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,gas_lbound:ngas),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(u_loc(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('u_loc',ierr)
      allocate(source_j(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(source_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_omp',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)
      if (multiple_integration) then
         allocate(a_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('a_spc',ierr)
         allocate(cj_ref_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('cj_ref_spc',ierr)
         allocate(cj_sup_ref_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('cj_sup_ref_spc',ierr)
         allocate(dominant_spc(1:ngas),stat=ierr)
         call CheckMemAlloc('dominant_spc',ierr)
         allocate(Fj_ref_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('Fj_ref_spc',ierr)
         allocate(Fj_sup_ref_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('Fj_sup_ref_spc',ierr)
         allocate(kappa_spc(0:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('kappa_spc',ierr)
         allocate(original_gas_index(1:ngas),stat=ierr)
         call CheckMemAlloc('original_gas_index',ierr)
         allocate(skip_gas(1:ngas),stat=ierr)
         call CheckMemAlloc('skip_gas',ierr)
         allocate(spc_gas_index(1:ngas,1:nspc),stat=ierr)
         call CheckMemAlloc('spc_gas_index',ierr)
         allocate(xs_avg(1:nspc),stat=ierr)
         call CheckMemAlloc('xs_avg',ierr)
         allocate(xs_cur(1:nspc),stat=ierr)
         call CheckMemAlloc('xs_cur',ierr)
      else
         allocate(cj_ref(0:ngas),stat=ierr)
         call CheckMemAlloc('cj_ref',ierr)
         allocate(cj_sup_ref(0:ngas),stat=ierr)
         call CheckMemAlloc('cj_sup_ref',ierr)
         allocate(Fj_ref(0:ngas),stat=ierr)
         call CheckMemAlloc('Fj_ref',ierr)
         allocate(Fj_sup_ref(0:ngas),stat=ierr)
         call CheckMemAlloc('Fj_sup_ref',ierr)
      endif

      !-----------------------------------------------------------------
      !Set up emission and absorption correction factors
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: &
                  &Set up emission and absorption correction factors')
      ecor = 1._dp; acor = 1._dp
      if (present(emission_correction))   ecor = emission_correction
      if (present(absorption_correction)) acor = absorption_correction

      !-----------------------------------------------------------------
      !Mount data_array
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Mount data_array')
      idata = 0
      do i=1,nx
         do j=1,ny
            do k=1,nz
               idata = idata + 1
               data_array(idata,1:3) = (/ i,j,k /)
            enddo
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !-----------------------------------------------------------------
      !Blackbody intensity
      call dprint('slw_gray_solution: Compute blackbody intensity')
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (slw_bound_Ib) then
                  F = bb_emission_frac(slw_Ib_lbound,T(i,j,k),.true.) -&
                      bb_emission_frac(slw_Ib_ubound,T(i,j,k),.true.)
               else
                  F = 1._dp
               endif
               Ib(i,j,k) = F*sigrpi*(T(i,j,k)**4._dp)
            enddo
         enddo
      enddo

      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !-----------------------------------------------------------------
      !Define wall emissivities
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Define wall emissivities')
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Initialize arrays
      ieps = 1                                                          !Assume a single band
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (i.eq.1)  x_eps(1,j,k) = xmin_emissivity(j,k,ieps)
               if (i.eq.nx) x_eps(2,j,k) = xmax_emissivity(j,k,ieps)
                  
               if (j.eq.1)  y_eps(1,i,k) = ymin_emissivity(i,k,ieps)
               if (j.eq.ny) y_eps(2,i,k) = ymax_emissivity(i,k,ieps)
                  
               if (k.eq.1)  z_eps(1,i,j) = zmin_emissivity(i,j,ieps)
               if (k.eq.nz) z_eps(2,i,j) = zmax_emissivity(i,j,ieps)
            enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      !Load the ALBDF database
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Load the ALBDF')
      call load_albdf

      !-----------------------------------------------------------------
      !Define the reference state, if necessary
      !-----------------------------------------------------------------
      if (trim(slw_nonuniform_method).ne.'uniform') then
         call dprint('slw_gray_solution: Define the reference state')
         slw_Tref = maxval(T)
         slw_pref = maxval(p)
         do ispc=1,nspc
            slw_xsref(ispc) = maxval(xs(ispc,:,:,:))
         enddo
      endif

      !-----------------------------------------------------------------
      !Compute scaling parameter, if necessary
      !-----------------------------------------------------------------
      u_loc = 1._dp
      if (trim(slw_nonuniform_method).eq.'scaled') then
         call dprint('slw_gray_solution: Compute scaling parameter')
         do idata=1,ndata
            i = data_array(idata,1)                                     !Define i,
            j = data_array(idata,2)                                     !  j and
            k = data_array(idata,3)                                     !  k indexes
            kp_loc = get_slw_kp(T(i,j,k),p(i,j,k),xs(:,i,j,k),&         !Planck-mean absorption coefficient
                                T(i,j,k),ngas)                          !  at the local thermodynamic state                  
            kp_ref = get_slw_kp(slw_Tref,slw_pref,slw_xsref,&           !Planck-mean absorption coefficient
                                T(i,j,k),ngas)                          !  at the reference thermodynamic state
            u_loc(i,j,k) = kp_loc/(kp_ref+small)                        !Scaling parameter
         enddo
      endif

      !-----------------------------------------------------------------
      !Prepare superposition arrays
      !-----------------------------------------------------------------
      if (multiple_integration) then
          call dprint('slw_gray_solution: &
                      &Prepare superposition arrays')
      
         !Mount the spc_gas_index array 
         !(array storing the gray gas index of each individual
         !species ispc that corresponds to the overal gray gas jgas)
         spc_gas_index = 0                                              !For jgas = 1, all species' gray gases are
                                                                        !  transparent windows
         do jgas=2,ngas
            do ispc=1,nspc        
               if (ispc.eq.1) then                                         !The first species is treated normally,
                  spc_gas_index(jgas,ispc) = &                             !  with each new jgas corresponding to a
                     mod((spc_gas_index(jgas-1,ispc)+1),&                  !  new j for the species, and the counter
                        (slw_ngas_per_species(ispc)+1))                           !  resetting once j = J
               else                                                        !For all subsequent species, a new step in
                  spc_gas_index(jgas,ispc) = spc_gas_index(jgas-1,ispc)    !  j is only taken (advance_index = .true)
                  advance_index = .true.                                   !  if j = J @ jgas-1 for all previous species
                  do iispc=1,ispc-1                                        !This loops checks if j = J @ jgas-1 for all
                     if (spc_gas_index(jgas-1,iispc).ne.&                  !  previous species. If that is the case, then
                        slw_ngas_per_species(iispc)) advance_index = .false.      !  advance_index = .true.
                  enddo                                                    !And, if that is the case, than the j stepping
                  if (advance_index) spc_gas_index(jgas,ispc) = &          !  is done, once again resetting if j = J for
                     mod((spc_gas_index(jgas,ispc)+1),&                    !  that species
                        (slw_ngas_per_species(ispc)+1)) 
               endif
            enddo 
         enddo

         !Define the Fj and Cj divisions for each species separately
         do ispc=1,nspc
            ncj = slw_ngas_per_species(ispc)
            if (trim(slw_nonuniform_method).eq.'rank_correlated') then
               !For the Rank-SLW, divide the reference Fj
               call get_slw_Fj(ncj,Fj_ref_spc(1:ncj,ispc),&
                  Fj_sup_ref_spc(1:ncj,ispc))
               Fj_sup_ref_spc(0,ispc) = slw_Fmin  
            else
               !For other SLWs, divide the reference Cj
               call get_slw_cj(slw_cmin,slw_cmax,ngas,&
                  cj_ref_spc(1:ncj,ispc),cj_sup_ref_spc(1:ncj,ispc))
               cj_sup_ref_spc(0,ispc) = slw_cmin
            endif
         enddo

         !Compute the average composition
         if (improved_multiple_integration) &
            call get_average_composition(T_avg,xs_avg,p_avg,'T4')

         !Mount the dominant_spc array
         !(this array informs which species is the dominant one, in the
         !framework of the improved superposition method; a value of -1
         !indicates that no species dominates)
         dominant_spc(1:ngas) = -1                                      !Initially, assume that no species dominate
         if (improved_multiple_integration) then
            !First compute kappa for each species and each species'
            !gray gas separately
            do ispc=1,nspc
               xs_cur = 0._dp; xs_cur(ispc) = xs_avg(ispc)
               do jgas=0,slw_ngas_per_species(ispc)
                  transparent_window = .false.
                  if (jgas.eq.0) transparent_window = .true.
                  call compute_slw_gas_parameters(kappa_spc(jgas,ispc),&
                     a_spc(jgas,ispc),cj_sup_ref_spc(jgas,ispc),&
                     cj_sup_loc,Fj_ref_spc(jgas,ispc),&
                     Fj_sup_ref_spc(jgas,ispc),F0,T_avg,&
                     xs_cur,transparent_window=transparent_window)!,&
!                     scaling_coeff=u_loc(i,j,k))
               enddo
            enddo
            
            do jgas=1,ngas                                              !Now loop over all mixture gray gases
               !Compute reference values of kappa
               kappa_mix = 0._dp
               do ispc=1,nspc
                  kappa_mix = kappa_mix + &
                     kappa_spc(spc_gas_index(jgas,ispc),ispc)
               enddo
               do ispc=1,nspc                                           !For each species, compute kappa of the mixture minus
                  kappa_comp = kappa_mix - &                            !  the current species and compare to kappa of the 
                     kappa_spc(spc_gas_index(jgas,ispc),ispc)           !  species, taking into account the safety factor theta
                  if (kappa_spc(spc_gas_index(jgas,ispc),ispc).gt.&     !If the comparison is successful, flag that as the
                      kappa_comp*theta) then                            !  dominant species for this jgas
                     dominant_spc(jgas) = ispc                             
                     exit                                                  
                  endif
               enddo
            enddo
         endif

         !Mount the skip_gas array
         !(this array informs if a jgas gas should be skipped when 
         !solving the RTE in the improved superposition method to 
         !prevent repeated calculations. This is done by sweeping 
         !through all gray gases larger than the current one and 
         !checking if the main condition of the improved method is not 
         !verified for any other combination of species' gray gases that
         !include the current gas for the dominating species species
         do jgas=1,ngas
            idsp = dominant_spc(jgas)
            if (idsp.le.0) cycle
            skip_gas(jgas) = .false.                                    !In principle, do not skip this gas
            gas_check_loop: do jjgas=jgas+1,ngas                        !Loop over all remaining gray gases, but consider only
               if (spc_gas_index(jjgas,idsp).ne.&                       !  jgas correspondig to the same gray gas index as that
                  spc_gas_index(jgas,idsp)) cycle gas_check_loop        !  of the current gray gas index of the dominant species
               if (dominant_spc(jjgas).eq.idsp) then                    !If, for any of the remaining gases, the current dominant
                  skip_gas(jgas) = .true.                               !  species dominates once again, we should flag jgas to
                  exit gas_check_loop                                   !  be skipped
               endif
            enddo gas_check_loop
         enddo      
      endif

      !-----------------------------------------------------------------
      !Compute radiative properties at each grid point and gray gas
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Compute local properties',3)
      kIb = 0._dp                                                       !Zero out array with the total RTE emission term
      if (multiple_integration) then
         do idata=1,ndata
            !Define i, j and k indexes
            i = data_array(idata,1)
            j = data_array(idata,2)
            k = data_array(idata,3)
            
            !First compute kappa for each species and each species'
            !gray gas separately
            do ispc=1,nspc
               xs_cur = 0._dp; xs_cur(ispc) = xs(ispc,i,j,k)
               do jgas=0,slw_ngas_per_species(ispc)
                  transparent_window = .false.
                  if (jgas.eq.0) transparent_window = .true.
                  call compute_slw_gas_parameters(kappa_spc(jgas,ispc),&
                     a_spc(jgas,ispc),cj_sup_ref_spc(jgas,ispc),&
                     cj_sup_loc,Fj_ref_spc(jgas,ispc),&
                     Fj_sup_ref_spc(jgas,ispc),F0,T(i,j,k),&
                     xs_cur,transparent_window=transparent_window,&
                     scaling_coeff=u_loc(i,j,k))
               enddo
            enddo
            
            gas_counter = 0
            do jgas=1,ngas
               idsp = dominant_spc(jgas)
               if (idsp.lt.0) then                                      !Standard superposition
                  gas_counter = gas_counter + 1
                  original_gas_index(gas_counter) = jgas
                  kappa(i,j,k,gas_counter) = 0._dp; a = 1._dp
                  do ispc=1,nspc
                     kappa(i,j,k,gas_counter) = &
                        kappa(i,j,k,gas_counter) + &
                           kappa_spc(spc_gas_index(jgas,ispc),ispc)
                     a = a*a_spc(spc_gas_index(jgas,ispc),ispc)
                  enddo
               else                                                     !Improved superposition
                  if (skip_gas(jgas)) cycle                             !Skip the jgas if the jgas is to be skipped
                  gas_counter = gas_counter + 1
                  original_gas_index(gas_counter) = jgas
                  kappa(i,j,k,gas_counter) = &                          !kappa for the gas is equal to kappa of the
                     kappa_spc(spc_gas_index(jgas,idsp),idsp)           !   dominant species
                  a = 1._dp
                  do ispc=1,nspc                                        !Initialize a as the value of a
                     a = a*a_spc(spc_gas_index(jgas,ispc),ispc)
                  enddo
                  do jjgas=1,jgas-1                                     !To compute a for the gas, loop over all previous jgas
                     if (dominant_spc(jjgas).ne.idsp) cycle             !  skipping jgas for which the dominant species is not 
                                                                        !  equal to the current dominant species, and
                     if (spc_gas_index(jjgas,idsp).ne.&                 !  jgas for which the dominant species' gray gas is 
                         spc_gas_index(jgas,idsp)) cycle                !  not equal to the current dominant species' gray gas
                     aa = 1._dp
                     do ispc=1,nspc
                        aa = aa*a_spc(spc_gas_index(jjgas,ispc),ispc)
                     enddo
                     a = a + aa
                  enddo
               endif  
               aIb(i,j,k,gas_counter) = a*Ib(i,j,k)                     !RTE emission term
               kIb(i,j,k) = kIb(i,j,k) + &                              !Total RTE emission term
                     kappa(i,j,k,gas_counter)*aIb(i,j,k,gas_counter)
            enddo
         enddo
      
         !Update the actual total number of gray gases
         ngas = gas_counter
         if (present(output_ngas)) output_ngas = ngas
      
      else

         !For the Rank-SLW, divide the reference Fj; for all other SLWs,
         !divide the reference Cj. If we are using multiple 
         !integration, this is done for each species one by one; 
         !otherwise, this is done only once, for the mixture
         if (trim(slw_nonuniform_method).eq.'rank_correlated') then
            call get_slw_Fj(ngas,Fj_ref(1:ngas),Fj_sup_ref(1:ngas))
               Fj_sup_ref(0) = slw_Fmin  
         else
            call get_slw_cj(slw_cmin,slw_cmax,ngas,cj_ref(1:ngas),&
                            cj_sup_ref(1:ngas))
            cj_sup_ref(0) = slw_cmin         
         endif

         !$OMP PARALLEL DO DEFAULT(SHARED) &
         !$OMP& PRIVATE(a,cj_sup_loc,F0,i,j,jgas,k,transparent_window) 
         do idata=1,ndata
            !Define i, j and k indexes
            i = data_array(idata,1)
            j = data_array(idata,2)
            k = data_array(idata,3)

            !Compute kappa and a for the mixture
            do jgas=0,ngas
               transparent_window = .false.                             !Set transparent window flag
               if (jgas.eq.0) transparent_window = .true.               ! for jgas = 0
               call compute_slw_gas_parameters(kappa(i,j,k,jgas),a,&
                  cj_sup_ref(jgas),cj_sup_loc,Fj_ref(jgas),&
                  Fj_sup_ref(jgas),F0,T(i,j,k),&
                  xs(:,i,j,k),transparent_window=transparent_window,&
                  scaling_coeff=u_loc(i,j,k))
               aIb(i,j,k,jgas) = a*Ib(i,j,k)                            !RTE emission term
               kIb(i,j,k) = kIb(i,j,k) + &                              !Total RTE emission term
                  kappa(i,j,k,jgas)*aIb(i,j,k,jgas)
            enddo
         enddo
         !$OMP END PARALLEL DO
      endif

      !-----------------------------------------------------------------
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Initialize sum arrays')
      flux = 0._dp; source = 0._dp
      start_proc_time = get_wall_time()                                 !Start the processing time counter
      
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_j,flux_omp,source_j,source_omp)
      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      flux_omp = 0._dp; source_omp = 0._dp
      
      !-----------------------------------------------------------------
      !Solve the radiation field
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Radiative transfer loop started')
      !$OMP DO
      gas_loop: do jgas=gas_lbound,ngas
         !Set up longer debug information
         write(gas_name,'(i5)') jgas
         msg = 'slw_black_solution: gas_loop for gas '//&
               trim(adjustl(gas_name))//' started'
         call dprint(msg)

         !Solve the radiation field for gas jgas
         call dprint('slw_black_solution: RTE solution',3)
         if (compute_radiation) &
            call fvm_solution(kappa(:,:,:,jgas),kappa(:,:,:,jgas),&
                              kappa_sct,aIb(:,:,:,jgas),phi,x_eps,&
                              y_eps,z_eps,flux_j,source_j,.true.)
                           
         !Updating the heat flux and heat source
         call dprint('slw_gray_solution: &
                     &Update the heat flux and heat source',3)
         flux_omp = flux_omp + flux_j
         source_omp = source_omp + source_j
      
      enddo gas_loop
      !$OMP ENDDO
      
      !-----------------------------------------------------------------
      !Finish the calculation of the total radiative quantities
      !-----------------------------------------------------------------
      !$OMP CRITICAL
      call dprint('slw_gray_solution: Finish the &
                  &calculation of the total radiative quantities')
      flux = flux + flux_omp
      source = source + source_omp
      !$OMP ENDCRITICAL
      !$OMP ENDPARALLEL
      
      !Stop processing time counter
      call dprint('slw_gray_solution: Radiative transfer loop ended')
      end_proc_time = get_wall_time()
      
!      !-----------------------------------------------------------------
!      !Dump of gas properties, if necessary
!      !-----------------------------------------------------------------
!      if (slw_print_properties) then
!         call dprint('slw_black_solution: Dump of gas properties')
!         do i=1,nx
!            do j=1,ny
!               do k=1,nz
!                  write(uout,'(1000(e26.15e3,:,","))')&
!                     x(i),y(j),z(k),(kappa(i,j,k,jgas),&
!                        aIb(i,j,k,jgas)/Ib(i,j,k),jgas=0,ngas)
!               enddo
!            enddo
!         enddo
         
!         !Close output unit
!         close(uout)

!      endif
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: &
                            &Finish calculation of optional properties')
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (compute_emission) &
                  emission(i,j,k) = fourpi*kIb(i,j,k)                   !Total emission
               if (compute_absorption) &
                  absorption(i,j,k) = source(i,j,k) + fourpi*kIb(i,j,k) !Total absorption
               if (compute_planck) &
                  kappaPlanck(i,j,k) = kIb(i,j,k)/Ib(i,j,k)             !Planck-mean absorption coefficient
            enddo
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Deallocate arrays')
      deallocate(aIb,Ib,kIb)
      deallocate(acor,ecor)
      deallocate(data_array)
      deallocate(flux_j,flux_omp,source_j,source_omp)
      deallocate(kappa,kappa_sct,phi)
      deallocate(u_loc)
      deallocate(x_eps,y_eps,z_eps)
      if (multiple_integration) then
         deallocate(a_spc,kappa_spc)
         deallocate(cj_ref_spc,cj_sup_ref_spc)
         deallocate(Fj_ref_spc,Fj_sup_ref_spc)
         deallocate(dominant_spc,skip_gas,spc_gas_index)
         deallocate(original_gas_index)
         deallocate(xs_avg,xs_cur)
      else
         deallocate(cj_ref,cj_sup_ref)
         deallocate(Fj_ref,Fj_sup_ref)
      endif
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      call dprint('slw_gray_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time
      
   endsubroutine slw_gray_solution

   !==================================================================
   !Subroutine for the radiative heat transfer solution via the
   !SLW-1 model for a non-scattering medium bounded by black surfaces
   !==================================================================
   subroutine slw1_black_solution(flux,source)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,print_to_prompt,shutdown
      use constants, only: sigrpi
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: number_of_species 
      use mesh
      use omp_lib
      implicit none
      integer :: i,ierr,ispc,j,jgas,k
      integer :: nl,nspc,nx,ny,nz
      real(dp),intent(out) :: &
         flux(3,xpoints,ypoints,zpoints),&
         source(xpoints,ypoints,zpoints)
      real(dp) :: p_ref,T_ref,xs_mix_ref
      real(dp) :: xs_mix
      real(dp) :: a_loc(0:1),a_ref(0:1),aux,c_loc_1,c_ref(0:1),&
                  kappa_ref(0:1)
      real(dp),allocatable,dimension(:) :: xs_ref
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: Ib,source_j
      real(dp),allocatable,dimension(:,:,:) :: kappa_sct,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: aIb,kappa
      real(dp),allocatable,dimension(:,:,:,:) :: flux_j

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Preparatory procedures')
      
      !Surrogate names
      nspc = number_of_species                                          !Number of species
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      
      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Allocating arrays')
      allocate(aIb(nx,ny,nz,0:1),stat=ierr)
      call CheckMemAlloc('a_jIb',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,0:1),stat=ierr)
      call CheckMemAlloc('kappa_j',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(source_j(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(xs_ref(nspc),stat=ierr)
      call CheckMemAlloc('xs_ref',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !-----------------------------------------------------------------
      !Define the reference state
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Define the reference state')
      T_ref = maxval(T)
      p_ref = maxval(p)
      do ispc=1,nspc
         xs_ref(ispc) = maxval(xs(ispc,:,:,:))
      enddo
      
      !-----------------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Computing properties that do &
                  &not depend on the gray gas')
      !Blackbody intensity
      do i=1,nx
         do j=1,ny
            do k=1,nz
               Ib(i,j,k) = sigrpi*(T(i,j,k)**4._dp)
            enddo
         enddo
      enddo

      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !Set wall emissivities to unit for all wall cells
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp

      !-----------------------------------------------------------------
      !Read the ALBDF
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Read the ALBDF')
      call load_albdf

      !-----------------------------------------------------------------
      !Computing gray gas reference properties
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Computing reference properties')
 
      !Compute overall mixture mole fraction
      xs_mix_ref = 0._dp
      do ispc=1,nspc
         if (is_slw_species(ispc)) &
            xs_mix_ref = xs_mix_ref + xs_ref(ispc)
      enddo
 
      !Properties for gas 1
      call slw1_compute_ref(T_ref,xs_ref,p_ref,kappa_ref(1),a_ref(1))
      c_ref(1) = slw_kappa_func(T_ref,p_ref,xs_mix_ref,kappa_ref(1),&
                                compute_cj=.true.,mixing_method='SLW1') !(163)
      aux = albdf_mix(T_ref,xs_ref,T_ref,c_ref(1),invert=.false.)

      !Properties for gas 0
      kappa_ref(0) = 0._dp
      a_ref(0) = 1._dp - a_ref(1)
      c_ref(0) = albdf_mix(T_ref,xs_ref,T_ref,a_ref(0),invert=.true.)   !(162)

      !-----------------------------------------------------------------
      !Computing local properties
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Computing local properties')
      do i=1,nx
         do j=1,ny
            do k=1,nz
               !Compute mixture overall mole fraction
               xs_mix = 0._dp
               do ispc=1,nspc
                  if (is_slw_species(ispc)) xs_mix = xs(ispc,i,j,k)
               enddo
               
               !Preliminary calculations
               c_loc_1 = albdf_mix(T(i,j,k),xs(:,i,j,k),T_ref,aux,&
                                   invert=.true.)                       !(164)
            
               !Weighting coefficient
               a_loc(0) = albdf_mix(T_ref,xs_ref,T(i,j,k),&
                                    c_ref(0),invert=.false.)            !(165)
               a_loc(1) = 1._dp - a_loc(0)                  
            
               !Absorption coefficient
               kappa(i,j,k,0) = 0._dp
               kappa(i,j,k,1) = slw_kappa_func(T(i,j,k),p(i,j,k),&
                                               xs_mix,c_loc_1,&
                                               mixing_method='SLW1')
            
               !Weighted blackbody function
               aIb(i,j,k,0) = a_loc(0)*Ib(i,j,k)
               aIb(i,j,k,1) = a_loc(1)*Ib(i,j,k)
            enddo
         enddo
      enddo 

      !-----------------------------------------------------------------
      !Solve the radiation field
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Solve the radiation field')
      
      !Initialize sum variables
      flux = 0; source = 0
      
      gas_loop: do jgas = 0,1
         
         !Solving the radiation field for gas jgas
         call fvm_solution(kappa(:,:,:,jgas),kappa(:,:,:,jgas),&
                     kappa_sct,aIb(:,:,:,jgas),phi,&
                     x_eps,y_eps,z_eps,flux_j,source_j)
                           
         !Updating the heat flux and heat source
         flux = flux + flux_j
         source = source + source_j
      
      enddo gas_loop
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('slw1_black_solution: Deallocate arrays')
      deallocate(aIb,Ib,kappa)
      deallocate(flux_j,source_j)
      deallocate(kappa_sct,phi,x_eps,y_eps,z_eps)
      
   endsubroutine slw1_black_solution

   !==================================================================
   !Subroutine for the radiative heat transfer solution via the
   !SLW-1 model for a non-scattering medium bounded by gray surfaces
   !==================================================================
   subroutine slw1_gray_solution(flux,source)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,print_to_prompt,shutdown
      use constants, only: sigrpi
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: number_of_species 
      use mesh
      use omp_lib
      implicit none
      integer :: i,ieps,ierr,ispc,j,jgas,k
      integer :: nl,nspc,nx,ny,nz
      real(dp),intent(out) :: &
         flux(3,xpoints,ypoints,zpoints),&
         source(xpoints,ypoints,zpoints)
      real(dp) :: p_ref,T_ref,xs_mix_ref
      real(dp) :: xs_mix
      real(dp) :: a_loc(0:1),a_ref(0:1),aux,c_loc_1,c_ref(0:1),&
                  kappa_ref(0:1)
      real(dp),allocatable,dimension(:) :: xs_ref
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: Ib,source_j
      real(dp),allocatable,dimension(:,:,:) :: kappa_sct,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: aIb,kappa
      real(dp),allocatable,dimension(:,:,:,:) :: flux_j

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Preparatory procedures')
      
      !Surrogate names
      nspc = number_of_species                                          !Number of species
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      
      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Allocating arrays')
      allocate(aIb(nx,ny,nz,0:1),stat=ierr)
      call CheckMemAlloc('a_jIb',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,0:1),stat=ierr)
      call CheckMemAlloc('kappa_j',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(source_j(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(xs_ref(nspc),stat=ierr)
      call CheckMemAlloc('xs_ref',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !-----------------------------------------------------------------
      !Define the reference state
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Define the reference state')
      T_ref = maxval(T)
      p_ref = maxval(p)
      do ispc=1,nspc
         xs_ref(ispc) = maxval(xs(ispc,:,:,:))
      enddo
      
      !-----------------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Computing properties that do &
                  &not depend on the gray gas')
      !Blackbody intensity
      do i=1,nx
         do j=1,ny
            do k=1,nz
               Ib(i,j,k) = sigrpi*(T(i,j,k)**4._dp)
            enddo
         enddo
      enddo

      !Scattering coefficient is null
      kappa_sct = 0._dp
      
      !Phase function (isotropic)
      phi = 1._dp
      
      !-----------------------------------------------------------------
      !Define wall emissivities
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Define wall emissivities')
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Initialize arrays
      ieps = 1                                                          !Assume a single band
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (i.eq.1)  x_eps(1,j,k) = xmin_emissivity(j,k,ieps)
               if (i.eq.nx) x_eps(2,j,k) = xmax_emissivity(j,k,ieps)
                  
               if (j.eq.1)  y_eps(1,i,k) = ymin_emissivity(i,k,ieps)
               if (j.eq.ny) y_eps(2,i,k) = ymax_emissivity(i,k,ieps)
                  
               if (k.eq.1)  z_eps(1,i,j) = zmin_emissivity(i,j,ieps)
               if (k.eq.nz) z_eps(2,i,j) = zmax_emissivity(i,j,ieps)
            enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      !Read the ALBDF
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Read the ALBDF')
      call load_albdf

      !-----------------------------------------------------------------
      !Computing gray gas reference properties
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Computing reference properties')
 
      !Compute overall mixture mole fraction
      xs_mix_ref = 0._dp
      do ispc=1,nspc
         if (is_slw_species(ispc)) &
            xs_mix_ref = xs_mix_ref + xs_ref(ispc)
      enddo
 
      !Properties for gas 1
      call slw1_compute_ref(T_ref,xs_ref,p_ref,kappa_ref(1),a_ref(1))
      c_ref(1) = slw_kappa_func(T_ref,p_ref,xs_mix_ref,kappa_ref(1),&
                                compute_cj=.true.,mixing_method='SLW1') !(163)
      aux = albdf_mix(T_ref,xs_ref,T_ref,c_ref(1),invert=.false.)

      !Properties for gas 0
      kappa_ref(0) = 0._dp
      a_ref(0) = 1._dp - a_ref(1)
      c_ref(0) = albdf_mix(T_ref,xs_ref,T_ref,a_ref(0),invert=.true.)   !(162)

      !-----------------------------------------------------------------
      !Computing local properties
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Computing local properties')
      do i=1,nx
         do j=1,ny
            do k=1,nz
               !Compute mixture overall mole fraction
               xs_mix = 0._dp
               do ispc=1,nspc
                  if (is_slw_species(ispc)) xs_mix = xs(ispc,i,j,k)

               enddo
               
               !Preliminary calculations
               c_loc_1 = albdf_mix(T(i,j,k),xs(:,i,j,k),T_ref,aux,&
                                   invert=.true.)                       !(164)
            
               !Weighting coefficient
               a_loc(0) = albdf_mix(T_ref,xs_ref,T(i,j,k),&
                                    c_ref(0),invert=.false.)            !(165)
               a_loc(1) = 1._dp - a_loc(0)      
!               print *, i, a_loc(0), a_loc(1)            
            
               !Absorption coefficient
               kappa(i,j,k,0) = 0._dp
               kappa(i,j,k,1) = slw_kappa_func(T(i,j,k),p(i,j,k),&
                                               xs_mix,c_loc_1,&
                                               mixing_method='SLW1')
!               print *, kappa(i,j,k,1)
               
               !Weighted blackbody function
               aIb(i,j,k,0) = a_loc(0)*Ib(i,j,k)
               aIb(i,j,k,1) = a_loc(1)*Ib(i,j,k)
            enddo
         enddo
      enddo 

      !-----------------------------------------------------------------
      !Solve the radiation field
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Solve the radiation field')
      
      !Initialize sum variables
      flux = 0; source = 0
      
      gas_loop: do jgas = 0,1
         
         !Solving the radiation field for gas jgas
         call fvm_solution(kappa(:,:,:,jgas),kappa(:,:,:,jgas),&
                     kappa_sct,aIb(:,:,:,jgas),phi,&
                     x_eps,y_eps,z_eps,flux_j,source_j)
                           
         !Updating the heat flux and heat source
         flux = flux + flux_j
         source = source + source_j
      
      enddo gas_loop
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('slw1_gray_solution: Deallocate arrays')
      deallocate(aIb,Ib,kappa)
      deallocate(flux_j,source_j)
      deallocate(kappa_sct,phi,x_eps,y_eps,z_eps)
      
   endsubroutine slw1_gray_solution

   !====================================================================
   !Subroutine for the exact solution of the RTE for a one-dimensional
   !cylindrical, non-scattering medium
   !====================================================================
!   subroutine slw_cyl_solution(x,T,p,xs,I0,n_gases,flux,source,&
!      processing_time,total_time,absorption,emission,kappaPlanck,&
!      solve_RTE,albdf_ready,points)
   
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: assert,CheckMemAlloc,get_file_unit,&
!                                get_wall_time,print_to_prompt,shutdown
!      use exact_solution, only: exact_1d_cyl
!      use global_parameters, only: debug_mode,number_of_species 
!      use physical_functions, only: bb_emission_frac
!      use slw_parameters, only: slw_bound_Ib,slw_Ib_lbound,slw_Ib_ubound
!      implicit none
!      character(10) :: gas_name
!      character(200) :: msg
!      integer,intent(in) :: n_gases
!      integer :: nx,ns
!      integer :: i,jgas,n
!      integer :: ierr,uout
!      integer,intent(in),optional :: points
!      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
!      real(dp),intent(out) :: flux(:),source(:)
!      real(dp),intent(out),optional :: processing_time,total_time,&
!                                absorption(:),emission(:),kappaPlanck(:)
!      real(dp) :: a,cj_loc,F
!      real(dp) :: xmix
!      real(dp) :: kp_loc,kp_ref
!      real(dp) :: end_proc_time,start_proc_time,&
!                  end_total_time,start_total_time,time_sum
!      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref,&
!                                           Fj_ref,Fj_sup_ref
!      real(dp),allocatable,dimension(:) :: aIb,Ib,kappa,kIb,u_loc,&
!         flux_j,source_j
!      real(dp),allocatable,dimension(:,:) :: cj_sup_loc,Fj_sup_loc
!      real(dp),allocatable,dimension(:,:) :: a_j,kappa_j
!      logical,optional :: albdf_ready,solve_RTE
!      logical :: albdf_read,transparent_window
!      logical :: compute_absorption,compute_emission,compute_planck,&
!                 compute_radiation
!start_total_time = I0 !Remove this
!      !-----------------------------------------------------------------
!      !Preparatory procedures
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                          &preparatory procedures')
      
!      !Starting the total computing time counter
!      start_total_time = get_wall_time()
      
!      !Check if all arrays have the same size
!      call assert(size(x).eq.size(T))
!      call assert(size(T).eq.size(p))
!      call assert(size(p).eq.size(xs,2))
      
!      !Zeroing out processing time counter 
!      !(necessary for the posterior sum over 
!      !all instances of the RTE solution)
!      time_sum = 0._dp
      
!      !Set flags for optional output arguments
!      compute_radiation = .true.
!      if (present(solve_RTE)) compute_radiation = solve_RTE
!      compute_absorption=.false.
!      if (present(absorption))  compute_absorption = .true.
!      compute_emission=.false.
!      if (present(emission))    compute_emission = .true.
!      compute_planck=.false.
!      if (present(kappaPlanck)) compute_planck = .true.
!      albdf_read = .false.                                              !By default, assume that the ALBDF data
!      if (present(albdf_ready)) albdf_read = albdf_ready                !  has not been loaded yet      

!      !-----------------------------------------------------------------
!      !Surrogate names
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                           &surrogate names')
!      nx = size(x); if (present(points)) nx = points                    !Number of grid points
!      ns = number_of_species                                            !Number of participating species
      
!      !-----------------------------------------------------------------
!      !Prepare for dumping properties, if requested
!      !-----------------------------------------------------------------
!      if (slw_print_properties) then
!         if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                       &prepare for dumping properties')
      
!         !Allocate extra arrays
!         allocate(a_j(nx,0:n_gases),stat=ierr)
!         call CheckMemAlloc('a_j',ierr)
!         allocate(kappa_j(nx,0:n_gases),stat=ierr)
!         call CheckMemAlloc('kappa_j',ierr)
         
!         !Prepare output unit
!         if (trim(slw_print_file).eq.'null') &                          !Check if a file name was given                          
!            call shutdown('slw_print_file undefined')
!         uout = get_file_unit()                                         !Get unit
!         open(unit=uout,file=trim(slw_print_file),&                     !Open file
!              form='formatted',action='write')
!         write(uout,'(a,(100(i3,:,",")))') '#0,',(i,i=1,2*(n_gases+1))  !Write first line
!         write(uout,'(3(a,:,","),100(i3,:,","))') '#j =',&              !Write header
!                                             (jgas,jgas,jgas=0,n_gases)
!         write(uout,'(100(a,:,","))') '#x',&                            !Write subheader
!                                    ('kappa_j','a_j',jgas=0,n_gases)
!         write(uout,'(100(a,:,","))') '#[m]',&                          !Write subsubheader
!                                    ('[1/m]','[-]',jgas=0,n_gases)
!      endif
      
!      !-----------------------------------------------------------------
!      !Allocating arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                          &allocating arrays')
!      allocate(aIb(nx),stat=ierr)
!      call CheckMemAlloc('aIb',ierr)
!      allocate(cj_ref(0:n_gases),stat=ierr)
!      call CheckMemAlloc('cj_ref',ierr)
!      allocate(cj_sup_loc(nx,0:n_gases),stat=ierr)
!      call CheckMemAlloc('cj_sup_loc',ierr)
!      allocate(cj_sup_ref(0:n_gases),stat=ierr)
!      call CheckMemAlloc('cj_sup_ref',ierr)
!      allocate(Fj_ref(0:n_gases),stat=ierr)
!      call CheckMemAlloc('Fj_ref',ierr)
!      allocate(Fj_sup_ref(0:n_gases),stat=ierr)
!      call CheckMemAlloc('Fj_sup_ref',ierr)
!      allocate(Fj_sup_loc(nx,0:n_gases),stat=ierr)
!      call CheckMemAlloc('Fj_sup_loc',ierr)
!      allocate(flux_j(nx),stat=ierr)
!      call CheckMemAlloc('flux_j',ierr)
!      allocate(Ib(nx),stat=ierr)
!      call CheckMemAlloc('Ib',ierr)
!      allocate(kappa(nx),stat=ierr)
!      call CheckMemAlloc('kappa',ierr)
!      allocate(kIb(nx),stat=ierr)
!      call CheckMemAlloc('kIb',ierr)
!      allocate(u_loc(nx),stat=ierr)
!      call CheckMemAlloc('u_loc',ierr)
!      allocate(source_j(nx),stat=ierr)
!      call CheckMemAlloc('source_j',ierr)

!      !-----------------------------------------------------------------
!      !Computing properties that do not depend on the gray gas
!      !-----------------------------------------------------------------
!      !Blackbody intensity
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                          &compute blackbody intensity')
!      do i=1,nx
!         if (slw_bound_Ib) then
!            F = bb_emission_frac(slw_Ib_lbound,T(i),.true.) -&
!                bb_emission_frac(slw_Ib_ubound,T(i),.true.)
!         else
!            F = 1._dp
!         endif
!            Ib(i) = F*sigma*invpi*(T(i)**4._dp)
!      enddo
      
!      !-----------------------------------------------------------------
!      !Defining absorption cross-sections and 
!      !supplementar absorption cross-sections
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                     &define absorption cross-sections')
!      if (trim(slw_nonuniform_method).eq.'rank_correlated') then
!         call get_slw_Fj(n_gases,Fj_ref(1:n_gases),&
!                         Fj_sup_ref(1:n_gases))
!         Fj_sup_ref(0) = slw_Fmin  

!      else
!         !For other SLWs, divide the reference Cj
!         call get_slw_cj(slw_cmin,slw_cmax,n_gases,cj_ref(1:n_gases),&
!                         cj_sup_ref(1:n_gases))
!         cj_sup_ref(0) = slw_cmin
!      endif

!      !-----------------------------------------------------------------
!      !If needed, read the ALBDF database
!      !-----------------------------------------------------------------
!      if (.not.albdf_read) then
!         if (debug_mode) &
!            call print_to_prompt('slw_cyl_solution: &
!                                 &read the ALBDF database')
!         call read_albdf_arrays_size
!         call allocate_albdf_parameters
!         call read_albdf_aux_arrays
!         call read_aldbf
!      endif

!      !-----------------------------------------------------------------
!      !Define the reference state, if necessary
!      !-----------------------------------------------------------------
!      if (trim(slw_nonuniform_method).ne.'uniform') then
!         if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                           &define the reference state')
!         slw_Tref = maxval(T)
!         slw_pref = maxval(p)
!         do n=1,ns
!           ! slw_xsref(n) = maxval(xs(n,:,:,:))
!         enddo
!      endif
      
!      !-----------------------------------------------------------------
!      !Compute scaling parameter, if necessary
!      !-----------------------------------------------------------------
!      if (trim(slw_nonuniform_method).eq.'scaled') then
!         if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                            &compute scaling parameter')
!         do i=1,nx
!            !Planck-mean absorption coefficients at
!            !local and reference conditions
!            kp_loc = get_slw_kp(T(i),p(i),xs(:,i),T(i),n_gases)
!            kp_ref = get_slw_kp(slw_Tref,slw_pref,slw_xsref,T(i),&
!                                n_gases)

!            !Scaling parameter
!            u_loc(i) = kp_loc/(kp_ref+small)
!         enddo
!      endif
       
!      !-----------------------------------------------------------------
!      !Initialize sum arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                           &initialize sum arrays')
!      flux = 0._dp
!      source = 0._dp
!      kIb = 0._dp
      
!      gas_loop: do jgas=0,n_gases

!         !--------------------------------------------------------------
!         !Preparatory procedures for the new gas
!         !--------------------------------------------------------------
!         !Longer debug information
!         if (debug_mode) then
!            write(gas_name,'(i5)') jgas
!            msg = 'slw_cyl_solution: gas_loop for gas '//&
!               trim(adjustl(gas_name))//' started'
!            call print_to_prompt(msg)
!         endif

!         !Check if the gas is a transparent window
!         transparent_window = .false.
!         if (jgas.eq.0) transparent_window = .true.

!         !--------------------------------------------------------------
!         !Computing local radiative properties
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                           &compute local properties',3)
!         do i=1,nx
!            !Mole fraction of the mixture of participating species
!            xmix = 0._dp
!            do n=1,ns
!               if ((albdf_nx(n).gt.0).and.(albdf_file(n).ne.'null')) &
!                  xmix = xmix + xs(n,i)
!            enddo

!            selectcase(trim(slw_nonuniform_method))
!            case('uniform')
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               !Uniform medium
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               !(for uniform media, cj_ref = cj_loc, 
!               !so use cj_ref always)
!               if (transparent_window) then
!                  kappa(i) = 0._dp                                      !Absorption coefficient of the transparent window
!                  a = slw_a_func(T(i),xs(:,i),T(i),cj_sup_ref(jgas))    !Weighting coefficient of the transparent window
!               else
!                  kappa(i) = slw_kappa_func(T(i),p(i),xmix,cj_ref(jgas))!Absorption coefficient of the gray gas
!                  a = slw_a_func(T(i),xs(:,i),T(i),cj_sup_ref(jgas),&   !Weighting coefficient of the gray gas
!                                 cj_sup_ref(jgas-1))
!               endif               
            
!            case('scaled')
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               !Scaled-SLW
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               if (transparent_window) then
!                  kappa(i) = 0._dp                                      !Absorption coefficient of the transparent window
!                  a = slw_a_func(slw_Tref,slw_xsref,T(i),&              !Weighting coefficient of the transparent window
!                                 cj_sup_ref(jgas))
!               else
!                  kappa(i) = u_loc(i)*&                                 !Absorption coefficient of the gray gas
!                     slw_kappa_func(T(i),p(i),xmix,cj_ref(jgas))
!                  a = slw_a_func(slw_Tref,slw_xsref,T(i),&              !Weighting coefficient of the gray gas
!                                 cj_sup_ref(jgas),cj_sup_ref(jgas-1))                           
!               endif
               
!            case('reference_approach')
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               !RA-SLW
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               Fj_sup_ref(jgas) = get_albdf(slw_Tref,slw_xsref,&        !ALBDF evaluated at the reference state
!                                             slw_Tref,cj_sup_ref(jgas))
!               cj_sup_loc(i,jgas) = inverse_albdf(Fj_sup_ref(jgas),&    !Get the local supplementar cross-section
!                                                  T(i),xs(:,i),slw_Tref)!  by solving the implicit equation                                      
!               if (transparent_window) then
!                  kappa(i) = 0._dp                                      !Absorption coefficient of the transparent window
!                  a = slw_a_func(slw_Tref,slw_xsref,T(i),&              !Weighting coefficient of the transparent window
!                                 cj_sup_ref(jgas)) 
!               else
!                  cj_loc = get_slw_single_cj(cj_sup_loc(i,jgas-1),&     !Local cross-section determined
!                                             cj_sup_loc(i,jgas))        !  from the supplementar ones
!                  kappa(i) = slw_kappa_func(T(i),p(i),xmix,cj_loc)      !Absorption coefficient of the gray gas
!                  a = slw_a_func(slw_Tref,slw_xsref,T(i),&              !Weighting coefficient of the gray gas
!                                 cj_sup_ref(jgas),cj_sup_ref(jgas-1))                                 
!               endif
               
!            case('rank_correlated')
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               !RC-SLW
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               cj_sup_loc(i,jgas) = inverse_albdf(Fj_sup_ref(jgas),&    !Get the local supplementar cross-section
!                                                  T(i),xs(:,i),slw_Tref)!  by solving the implicit equation using
!                                                                        !  the previously divided supplementar Fs
!               if (transparent_window) then
!                  kappa(i) = 0._dp                                      !Absorption coefficient of the transparent window
!                  a = slw_a_func(T(i),xs(:,i),T(i),cj_sup_loc(i,jgas))  !Weighting coefficient of the transparent window
!               else
!                  cj_loc = inverse_albdf(Fj_ref(jgas),T(i),xs(:,i),&    !Get the local cross-section by solving the implicit
!                                         slw_Tref)                      !   equation using the previously defined Fs
!                  kappa(i) = slw_kappa_func(T(i),p(i),xmix,cj_loc)      !Absorption coefficient of the gray gas
!                  a = slw_a_func(T(i),xs(:,i),T(i),cj_sup_loc(i,jgas),&  !Weighting coefficient of the gray gas
!                                 cj_sup_loc(i,jgas-1))
!               endif
                     
!            case('locally_correlated')
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               !LC-SLW
!               !++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!               Fj_sup_loc(i,jgas) = get_albdf(slw_Tref,slw_xsref,T(i),& !Local ALBDF evaluated using reference state                        
!                                              cj_sup_ref(jgas))
!               cj_sup_loc(i,jgas) = inverse_albdf(Fj_sup_loc(i,jgas),&  !Get the local supplementar cross-section
!                                                  T(i),xs(:,i),T(i))    !  by solving the implicit equation
!               if (transparent_window) then
!                  kappa(i) = 0._dp                                      !Absorption coefficient of the transparent window
!                  a = Fj_sup_loc(i,jgas)                                !Weighting coefficient of the transparent window
!               else
!                  cj_loc = get_slw_single_cj(cj_sup_loc(i,jgas-1),&     !Local cross-section determined
!                                             cj_sup_loc(i,jgas))        !  from the supplementar ones
!                  kappa(i) = slw_kappa_func(T(i),p(i),xmix,cj_loc)      !Absorption coefficient of the  gray gas
!                  a = Fj_sup_loc(i,jgas) - Fj_sup_loc(i,jgas-1)         !Weighting coefficient of the gray gas                           
!               endif
               
!            case default
!               call shutdown('slw_cyl_solution: Problem in the &
!                  &specification of slw_nonuniform_method')
!            endselect
!            aIb(i) = a*Ib(i)                                            !Emission term
!            kIb(i) = kIb(i) + kappa(i)*aIb(i)                           !Total RTE emission term
                  
!            !Store properties of each gray gas if requested
!            if (slw_print_properties) then
!               kappa_j(i,jgas) = kappa(i)
!               a_j(i,jgas) = a
!            endif
!         enddo

!         !--------------------------------------------------------------
!         !Solve the radiation field
!         !--------------------------------------------------------------
!         !Starting the processing time counter
!         start_proc_time = get_wall_time()
   
!         !Solving the radiation field for gas jgas
!         if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                              &solve radiation field',3)
!         if (compute_radiation) &
!            call exact_1d_cyl(x,kappa,aIb,0._dp,nx,source_j)            !Temporary
                           
!         !Updating the heat flux and heat source
!         if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                               &update the heat flux and heat source',3)
!         flux = flux + flux_j
!         source = source + source_j
      
!         !Stopping the processing time counter and adding to the total
!         end_proc_time = get_wall_time()
!         time_sum = time_sum + end_proc_time - start_proc_time
      
!      enddo gas_loop
      
!      !-----------------------------------------------------------------
!      !Dump of gas properties, if necessary
!      !-----------------------------------------------------------------
!      if (slw_print_properties) then
!         if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                               &dump of gas properties')
!         do i=1,nx
!            write(uout,'(1000(e26.15e3,:,","))') x(i),(kappa_j(i,jgas),&
!               a_j(i,jgas),jgas=0,n_gases)
!         enddo
         
!         !Close output unit
!         close(uout)
      
!         !Deallocate additional arrays
!         deallocate(a_j,kappa_j)
!      endif
      
!      !-----------------------------------------------------------------
!      !Final loop over all grid cells to finish up
!      !the calculation of the the optional properties
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                            &finish calculation of optional properties')
!      do i=1,nx
!         if (compute_emission) emission(i) = fourpi*kIb(i)              !Total emission
!         if (compute_absorption) &
!            absorption(i) = source(i) + fourpi*kIb(i)                   !Total absorption
!         if (compute_planck) kappaPlanck(i) = kIb(i)/Ib(i)              !Planck-mean absorption coefficient

!      enddo
      
!      !-----------------------------------------------------------------
!      !Deallocate arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                           &deallocate arrays')
!      deallocate(aIb,Ib,kIb)
!      deallocate(cj_ref,cj_sup_loc,cj_sup_ref)
!      deallocate(Fj_ref,Fj_sup_ref,Fj_sup_loc)
!      deallocate(flux_j,source_j)
!      deallocate(kappa)
!      deallocate(u_loc)
      
!      !-----------------------------------------------------------------
!      !Total times
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('slw_cyl_solution: &
!                                           &compute total times')
!      if (present(processing_time)) processing_time = time_sum          !Processing
!      end_total_time = get_wall_time(); if (present(total_time)) &      !Computing (total)
!                        total_time = end_total_time - start_total_time
      
!   endsubroutine slw_cyl_solution

endmodule slw_routines
