module gg_routines

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp
   use gg_functions, only: gg_kappa_func
   use gg_parameters, only: gg_print_file,gg_print_properties,&
      which_gg_model
   implicit none
   
contains
   
   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !GG model for a non-scattering medium bounded by black walls
   !====================================================================
   subroutine gg_black_solution(flux,source,processing_time,total_time,&
               wsgg_spec,absorption,emission,kappaPlanck,fixed_kappa)
                                 
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,get_file_unit,&
                                get_wall_time,print_to_prompt,shutdown
      use constants, only: fourpi,sigrpi
      use exact_solution, only: exact_1d_pp
      use fvm_parameters
      use fvm_routines, only: fvm_solution
      use global_parameters, only: rte_solution_method,debug_mode
      use mesh   
      integer :: i,j,k,nl,nx,ny,nz
      integer :: ierr,uout
      character(*),optional :: wsgg_spec
      character(100) :: wsgg_spec_aux
      real(dp),intent(in),optional :: &
         fixed_kappa(xpoints,ypoints,zpoints)
      real(dp),intent(out) :: flux(1:3,xpoints,ypoints,zpoints),&
         source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: &
         absorption(xpoints,ypoints,zpoints),&
         emission(xpoints,ypoints,zpoints),&
         kappaPlanck(xpoints,ypoints,zpoints),processing_time,total_time
      real(dp) :: end_proc_time,start_proc_time,end_total_time,&
         start_total_time
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: kappa,kappa_sct,Ib,&
         x_eps,y_eps,z_eps
      logical :: compute_absorption=.false.,compute_emission=.false.,&
         compute_planck=.false.
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_black_solution: &
                                           &preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Setting default parameter
      wsgg_spec_aux = 'null'
      if(present(wsgg_spec)) wsgg_spec_aux = wsgg_spec
      
      !Set flags for optional output arguments
      if (present(absorption)) compute_absorption = .true.
      if (present(emission)) compute_emission = .true.
      if (present(kappaPlanck)) compute_planck = .true.
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_black_solution: &
                                           &surrogate names')
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization

      !-----------------------------------------------------------------
      !Prepare output file if requested
      !-----------------------------------------------------------------
      if (gg_print_properties) then
         if (debug_mode) call print_to_prompt('gg_black_solution: &
                                     &prepare output file if requested')
         if (trim(gg_print_file).eq.'null') &                           !Check if a file name was given                          
            call shutdown('gg_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(gg_print_file))                       !Open file
         write(uout,'(4(a,:,","))') 'x','y','z','kappa'                 !Write header
         write(uout,'(4(a,:,","))') '[m]','[m]','[m]','[1/m]'           !Write subheader
      endif

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_black_solution: &
                                           &allocate arrays')
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(kappa(xpoints,ypoints,zpoints),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(xpoints,ypoints,zpoints),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(Ib(xpoints,ypoints,zpoints),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(x_eps(1:2,ypoints,zpoints),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(1:2,xpoints,zpoints),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(1:2,xpoints,ypoints),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)
      
      !-----------------------------------------------------------------
      !Precompute local radiative properties
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_black_solution: &
                                &precompute local radiative properties')
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (trim(which_gg_model).eq.'fixed') then                !If a fixed absorption coefficient is imposed,
                  if (.not.present(fixed_kappa)) call shutdown(&        !  check if a value for the absorption coefficient 
                     'gg_black_solution: fixed_kappa unspecified')      !  was given (through fixed_kappa),
                  kappa(i,j,k) = fixed_kappa(i,j,k)                     !  and assign the value to the absorption coefficient
               else                                                     !Otherwise,
                  kappa(i,j,k) = gg_kappa_func(T(i,j,k),p(i,j,k),&      !  compute the local mixture absorption coefficient
                                              xs(:,i,j,k),wsgg_spec_aux)!  using the dedicated function                 
               endif
               Ib(i,j,k) = sigrpi*(T(i,j,k)**4._dp)                     !Local blackbody intensity
               if (gg_print_properties) &                               !Print local properties if requested
                  write(uout,'(4(e26.15e3,:,","))') &
                     x(i),y(j),z(k),kappa(i,j,k)
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
      
      !-----------------------------------------------------------------
      !Solve the radiation field
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_black_solution: &
                                           &solve the radiation field')
      !Starting the processing time counter
      start_proc_time = get_wall_time()
      
      !Apply chosen RTE solution method
      selectcase(rte_solution_method)
      case('exact-1d')
         call exact_1d_pp(x,kappa(:,1,1),Ib(:,1,1),nx,source(:,1,1),&
                          flux(1,:,1,1))
      case('FVM')                                                       !Finite volume method solution
         call fvm_solution(kappa,kappa,kappa_sct,Ib,phi,&
                           x_eps,y_eps,z_eps,flux,source,.true.)
!      case('P1')
!         call p1_solution(kappa,Ib,flux,source)                         !P1 method
!      case('MC')                                                        !Monte Carlo method
!         call mc_solution(kappa,kappa,Ib,flux,source)      
      endselect
      
      !Stop the processing time counter
      end_proc_time = get_wall_time()

      !-----------------------------------------------------------------
      !Close output unit if it was used
      !-----------------------------------------------------------------
      if (gg_print_properties) close(uout)

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_black_solution: &
                    &finish calculation of the the optional properties')
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (compute_emission) &
                  emission(i,j,k) = fourpi*kappa(i,j,k)*Ib(i,j,k)       !Total emission
               if (compute_absorption) &
                  absorption(i,j,k) = source(i,j,k) + &
                     fourpi*kappa(i,j,k)*Ib(i,j,k)                      !Total absorption
               if (compute_planck) &
                  kappaPlanck(i,j,k) = kappa(i,j,k)                     !Planck-mean absorption coefficient
            enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_black_solution: &
                                           &deallocate arrays')
      deallocate(phi)
      deallocate(kappa,kappa_sct)
      deallocate(Ib)
      deallocate(x_eps,y_eps,z_eps)
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_black_solution: &
                                           &compute total times')
      !Processing
      if (present(processing_time)) &
         processing_time = end_proc_time - start_proc_time

      !Computing (total)
      end_total_time = get_wall_time()
      if (present(total_time)) &
         total_time = end_total_time - start_total_time

   endsubroutine gg_black_solution

!   !=============================================
!   !Subroutine for the solution of the RTE along 
!   !a line-of-sight for a non-scattering medium
!   !=============================================
!   subroutine gg_los_solution(gg_spec,I0,Irad,proc_time,total_time,&
!      wsgg_spec,absorption,emission,kappaPlanck,fixed_kappa)
   
!      !-------------------------
!      !Declaration of variables
!      !-------------------------
!      use comp_functions, only: CheckMemAlloc,get_file_unit,shutdown
!      use constants, only: sigma,invpi
!      use exact_solution, only: exact_los_nonsct
!      use mesh     
!      integer :: i,j,k,nx
!      integer :: ierr,uout
!      character(*),intent(in) :: gg_spec
!      character(*),optional :: wsgg_spec
!      character(100) :: wsgg_spec_aux
!      real(dp),intent(in) :: I0
!      real(dp),intent(in),optional :: fixed_kappa
!      real(dp),intent(out) :: proc_time,total_time,Irad(xpoints)
!      real(dp),intent(out),optional :: absorption(xpoints),&
!         emission(xpoints),kappaPlanck(xpoints)
!      real(dp) :: end_proc_time,start_proc_time,end_total_time,&
!         start_total_time
!      real(dp),allocatable,dimension(:) :: kappa,Ib
!      logical :: compute_absorption=.false.,compute_emission=.false.,&
!         compute_planck=.false.
      
!      !-----------------------
!      !Preparatory procedures
!      !-----------------------
!      !Starting the total computing time counter
!      call cpu_time(start_total_time)
      
!      !Setting default parameter
!      if(.not.present(wsgg_spec)) then
!         wsgg_spec_aux = 'null'
!      else
!         wsgg_spec_aux = wsgg_spec
!      endif
      
!      !Set flags for optional output arguments
!      if (present(absorption)) compute_absorption = .true.
!      if (present(emission)) compute_emission = .true.
!      if (present(kappaPlanck)) compute_planck = .true.
      
!      !----------------
!      !Surrogate names
!      !----------------
!      nx = xpoints                                                      !Number of cells along the path

!      !---------------------------------
!      !Prepare output unit if requested
!      !---------------------------------
!      if (gg_print_properties) then
!         if (trim(gg_print_file).eq.'null') &                           !Check if a file name was given                          
!            call shutdown('gg_print_file undefined')
!         uout = get_file_unit()                                         !Get unit
!         open(unit=uout,file=trim(gg_print_file))                       !Open file
!         write(uout,'(3(a,:,","))') 's','kappa','Ib'                    !Write header
!         write(uout,'(3(a,:,","))') '[m]','[1/m]','[W/m²]'              !Write subheader
!      endif

!      !------------------
!      !Allocating arrays
!      !------------------
!      allocate(kappa(xpoints),stat=ierr)
!      call CheckMemAlloc('kappa',ierr)
!      allocate(Ib(xpoints),stat=ierr)
!      call CheckMemAlloc('Ib',ierr)
      
!      !----------------------------------------
!      !Precomputing local radiative properties
!      !----------------------------------------
!      j = 1; k = 1                                                      !Keep j and k fixed (1d case)
!      do i=1,nx
!         if (trim(gg_spec).eq.'fixed') then                             !If a fixed absorption coefficient is imposed,
!            if (.not.present(fixed_kappa)) call shutdown(&              !  check if a value for the absorption coefficient 
!               'gg_black_solution: const_kappa unspecified')            !  was given (through fixed_kappa),
!            kappa(i) = fixed_kappa                                      !  and assign the value to the absorption coefficient
!         else                                                           !Otherwise,
!            kappa(i) = gg_kappa_func(gg_spec,T(i,j,k),p(i,j,k),&        !  compute the local mixture absorption coefficient
!                     xh2o(i,j,k),xco2(i,j,k),xco(i,j,k),xch4(i,j,k),&  !  using the dedicated function
!                     xsoot(i,j,k),wsgg_spec_aux)                        
!         endif
!         Ib(i) = invpi*sigma*(T(i,j,k)**4._dp)                          !Local blackbody intensity
!         if (gg_print_properties) &                                     !Print local properties if requested
!            write(uout,'(4(e26.15e3,:,","))') x(i),kappa(i),Ib(i)
!      enddo
      
!      !------------------------------
!      !Solving radiation intensities
!      !------------------------------
!      !Starting the processing time counter
!      call cpu_time(start_proc_time)
      
!      !Solve the RTE
!      call exact_los_nonsct(x,kappa,Ib,I0,nx,Irad)
      
!      !Stopping the processing time counter
!      call cpu_time(end_proc_time)

!      !---------------------------------
!      !Close output unit if it was used
!      !---------------------------------
!      if (gg_print_properties) close(uout)

!      !-----------------------------------------------
!      !Final loop over all grid cells to finish up
!      !the calculation of the the optional properties
!      !-----------------------------------------------
!      do i=1,nx
!         if (compute_emission) emission(i) = kappa(i)*Ib(i)             !Emission                              
!         if (compute_absorption) absorption(i) = kappa(i)*Irad(i)       !Absorption
!         if (compute_planck) kappaPlanck(i) = kappa(i)                  !Planck-mean absorption coefficient
!      enddo

!      !------------------
!      !Deallocate arrays
!      !------------------
!      deallocate(kappa)
!      deallocate(Ib)
      
!      !------------
!      !Total times
!      !------------
!      !Processing
!      proc_time = end_proc_time - start_proc_time

!      !Computing (total)
!      call cpu_time(end_total_time)
!      total_time = end_total_time - start_total_time

!   endsubroutine gg_los_solution

   !====================================================================
   !Subroutine for the exact solution of the RTE 
   !along for cylindrical, non-scattering medium
   !====================================================================
   subroutine gg_cyl_solution(x,T,p,xs,I0,source,flux,processing_time,&
                           total_time,absorption,emission,kappaPlanck,&
                           fixed_kappa,points,wsgg_spec)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
                  get_file_unit,get_wall_time,print_to_prompt,shutdown
      use constants, only: fourpi,invpi,sigma
      use exact_solution, only: exact_1d_cyl!,exact_1d_cyl2
      use global_parameters, only: debug_mode,number_of_species
      integer,intent(in),optional :: points
      integer :: i,nsp,nx
      integer :: ierr,uout
      character(*),optional :: wsgg_spec
      character(100) :: wsgg_spec_aux
      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
      real(dp),intent(in),optional :: fixed_kappa(:)
      real(dp),intent(out) :: flux(:),source(:)
      real(dp),intent(out),optional :: absorption(:),emission(:),&
         kappaPlanck(:),processing_time,total_time
      real(dp) :: end_proc_time,start_proc_time,end_total_time,&
         start_total_time
      real(dp),allocatable,dimension(:) :: kappa,Ib
      logical :: compute_absorption=.false.,compute_emission=.false.,&
         compute_planck=.false.
start_total_time = I0 !Remove this
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                                           &preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Check if all arrays have the same size
      call assert(size(x).eq.size(T))
      call assert(size(T).eq.size(p))
      call assert(size(p).eq.size(xs,2))
      
      !Setting default parameter
      wsgg_spec_aux = 'null'
      if(present(wsgg_spec)) wsgg_spec_aux = wsgg_spec
      
      !Set flags for optional output arguments
      compute_absorption=.false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission=.false.
      if (present(emission))    compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                                           &surrogate names')
      nx = size(x); if (present(points)) nx = points                    !Number of grid points
      nsp = number_of_species                                           !Number of participating species

      !-----------------------------------------------------------------
      !Prepare output unit if requested
      !-----------------------------------------------------------------
      if (gg_print_properties) then
         if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                                              &prepare output unit')
         if (trim(gg_print_file).eq.'null') &                           !Check if a file name was given                          
            call shutdown('gg_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(gg_print_file))                       !Open file
         write(uout,'(3(a,:,","))') 'r','kappa','Ib'                    !Write header                  
         write(uout,'(3(a,:,","))') '[m]','[1/m]','[kW/m²]'             !Write subheader
      endif

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                                           &allocating arrays')
      allocate(kappa(nx),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(Ib(nx),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      
      !-----------------------------------------------------------------
      !Precomputing local radiative properties
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                                      &precompute radiative properties')
      do i=1,nx
         if (trim(which_gg_model).eq.'fixed') then                      !If a fixed absorption coefficient is imposed,
            if (.not.present(fixed_kappa)) call shutdown(&              !  check if a value for the absorption coefficient 
               'gg_cyl_solution: fixed_kappa array unspecified')        !  was given (through fixed_kappa),
            kappa(i) = fixed_kappa(i)                                   !  and assign the value to the absorption coefficient
         else                                                           !Otherwise,
             kappa(i) = gg_kappa_func(T(i),p(i),xs(:,i),wsgg_spec_aux)  !  compute the local mixture absorption coefficient
         endif
         Ib(i) = invpi*sigma*(T(i)**4._dp)                              !Local blackbody intensity
         if (gg_print_properties) &                                     !Print local properties if requested
            write(uout,'(4(e26.15e3,:,","))') x(i),kappa(i),&
               Ib(i)/1000._dp
      enddo
      
      !-----------------------------------------------------------------
      !Solving radiation intensities
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                                           &solve the radiation field')
      
      !Starting the processing time counter
      start_proc_time = get_wall_time()
      
      !Solve the RTE
      flux = 0._dp
      !call exact_1d_cyl2(x,kappa,Ib,0._dp,nx,source)
      call exact_1d_cyl(x,kappa,Ib,0._dp,nx,source)                     !Temporary
      
      !Stopping the processing time counter
      end_proc_time = get_wall_time()

      !-----------------------------------------------------------------
      !Close output unit if it was used
      !-----------------------------------------------------------------
      if (gg_print_properties) close(uout)

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                        &finish calculation of the optional properties')
      do i=1,nx
         if (compute_emission) emission(i) = fourpi*kappa(i)*Ib(i)      !Emission                              
         if (compute_absorption) &
            absorption(i) = source(i) + fourpi*kappa(i)*Ib(i)           !Total absorption
         if (compute_planck) kappaPlanck(i) = kappa(i)                  !Planck-mean absorption coefficient
      enddo

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                                           &deallocate arrays')
      deallocate(kappa)
      deallocate(Ib)
      
      !------------
      !Total times
      !------------
      if (debug_mode) call print_to_prompt('gg_cyl_solution: &
                                           &total times')
      !Processing
      if (present(processing_time)) &
         processing_time = end_proc_time - start_proc_time

      !Computing (total)
      call cpu_time(end_total_time)
      if (present(total_time)) &
         total_time = end_total_time - start_total_time
   
   endsubroutine gg_cyl_solution

   !====================================================================
   !Subroutine for the exact solution of the RTE 
   !for a plane-parallel, non-scattering medium
   !====================================================================
   subroutine gg_pp_solution(x,T,p,xs,source,flux,processing_time,&
                           total_time,absorption,emission,kappaPlanck,&
                           fixed_kappa,points,wsgg_spec)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
                  get_file_unit,get_wall_time,print_to_prompt,shutdown
      use constants, only: fourpi,invpi,sigma
      use exact_solution, only: exact_1d_pp
      use global_parameters, only: debug_mode,number_of_species
      integer,intent(in),optional :: points
      integer :: i,nsp,nx
      integer :: ierr,uout
      character(*),optional :: wsgg_spec
      character(100) :: wsgg_spec_aux
      real(dp),intent(in) :: p(:),T(:),x(:),xs(:,:)
      real(dp),intent(in),optional :: fixed_kappa(:)
      real(dp),intent(out) :: flux(:),source(:)
      real(dp),intent(out),optional :: absorption(:),emission(:),&
         kappaPlanck(:),processing_time,total_time
      real(dp) :: end_proc_time,start_proc_time,end_total_time,&
         start_total_time
      real(dp),allocatable,dimension(:) :: kappa,Ib
      logical :: compute_absorption=.false.,compute_emission=.false.,&
         compute_planck=.false.
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_pp_solution: &
                                           &preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Check if all arrays have the same size
      call assert(size(x).eq.size(T))
      call assert(size(T).eq.size(p))
      call assert(size(p).eq.size(xs,2))
      
      !Setting default parameter
      wsgg_spec_aux = 'null'
      if(present(wsgg_spec)) wsgg_spec_aux = wsgg_spec
      
      !Set flags for optional output arguments
      compute_absorption=.false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission=.false.
      if (present(emission))    compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_pp_solution: &
                                           &surrogate names')
      nx = size(x); if (present(points)) nx = points                    !Number of grid points
      nsp = number_of_species                                           !Number of participating species

      !-----------------------------------------------------------------
      !Prepare output unit if requested
      !-----------------------------------------------------------------
      if (gg_print_properties) then
         if (debug_mode) call print_to_prompt('gg_pp_solution: &
                                              &prepare output unit')
         if (trim(gg_print_file).eq.'null') &                           !Check if a file name was given                          
            call shutdown('gg_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(gg_print_file))                       !Open file
         write(uout,'(3(a,:,","))') 'r','kappa','Ib'                    !Write header                  
         write(uout,'(3(a,:,","))') '[m]','[1/m]','[kW/m²]'             !Write subheader
      endif

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_pp_solution: &
                                           &allocating arrays')
      allocate(kappa(nx),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(Ib(nx),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      
      !-----------------------------------------------------------------
      !Precomputing local radiative properties
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_pp_solution: &
                                      &precompute radiative properties')
      do i=1,nx
         if (trim(which_gg_model).eq.'fixed') then                      !If a fixed absorption coefficient is imposed,
            if (.not.present(fixed_kappa)) call shutdown(&              !  check if a value for the absorption coefficient 
               'gg_pp_solution: fixed_kappa array unspecified')         !  was given (through fixed_kappa),
            kappa(i) = fixed_kappa(i)                                   !  and assign the value to the absorption coefficient
         else                                                           !Otherwise,
             kappa(i) = gg_kappa_func(T(i),p(i),xs(:,i),wsgg_spec_aux)  !  compute the local mixture absorption coefficient
         endif
         Ib(i) = invpi*sigma*(T(i)**4._dp)                              !Local blackbody intensity
         if (gg_print_properties) &                                     !Print local properties if requested
            write(uout,'(4(e26.15e3,:,","))') x(i),kappa(i),&
               Ib(i)/1000._dp
!write(*,*) i,kappa(i),Ib(i)
      enddo
      
      !-----------------------------------------------------------------
      !Solving radiation intensities
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_pp_solution: &
                                           &solve the radiation field')
      
      !Starting the processing time counter
      start_proc_time = get_wall_time()
      
      !Solve the RTE
      call exact_1d_pp(x,kappa,Ib,nx,source,flux)
      
      !Stopping the processing time counter
      end_proc_time = get_wall_time()

      !-----------------------------------------------------------------
      !Close output unit if it was used
      !-----------------------------------------------------------------
      if (gg_print_properties) close(uout)

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_pp_solution: &
                        &finish calculation of the optional properties')
      do i=1,nx
         if (compute_emission) emission(i) = fourpi*kappa(i)*Ib(i)      !Emission                              
         if (compute_absorption) &
            absorption(i) = source(i) + fourpi*kappa(i)*Ib(i)           !Total absorption
         if (compute_planck) kappaPlanck(i) = kappa(i)                  !Planck-mean absorption coefficient
      enddo

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('gg_pp_solution: &
                                           &deallocate arrays')
      deallocate(kappa)
      deallocate(Ib)
      
      !------------
      !Total times
      !------------
      if (debug_mode) call print_to_prompt('gg_pp_solution: &
                                           &total times')
      !Processing
      if (present(processing_time)) &
         processing_time = end_proc_time - start_proc_time

      !Computing (total)
      call cpu_time(end_total_time)
      if (present(total_time)) &
         total_time = end_total_time - start_total_time
   
   endsubroutine gg_pp_solution

endmodule gg_routines
