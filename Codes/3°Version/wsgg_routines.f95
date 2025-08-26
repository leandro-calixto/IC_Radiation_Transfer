!#######################################################################
!Routines to solve the radiative heat transfer 
!via the WSGG model for different scenarios
!#######################################################################
module wsgg_routines

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp,small
   use constants, only: invpi,fourpi,pi,sigma
   use global_parameters
   use omp_lib
   implicit none
   
contains
   
   !====================================================================
   !Subroutine for the solution of the RTE along a line-of-sight for a
   !non-scattering medium. The code solves the RTE for a single gray
   !gas only, and outputs the corresponding intensity
   !====================================================================
   subroutine wsgg_sg_los_solution(x,T,p,xs,jgas,I0,Irad,wsgg_spec,&
      emission_correction,absorption_correction,processing_time,&
      total_time,absorption,emission,kappaPlanck,Planck,points,&
      solve_RTE,quick_RTE)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckMemAlloc,dprint,&
                                get_file_unit,get_wall_time,shutdown
      use constants, only: sigrpi
      use exact_solution, only: exact_los_nonsct
      use wsgg_functions, only: wsgg_properties
      implicit none
      character(*) :: wsgg_spec
      integer,intent(in) :: jgas
      integer,intent(in),optional :: points
      integer :: i,ierr,ns,nx
      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
      real(dp),intent(in),optional :: absorption_correction(:),&
                                      emission_correction(:)
      real(dp),intent(out),optional :: absorption(:),&
         emission(:),kappaPlanck(:),Planck(:),processing_time,total_time
      real(dp) :: a
      real(dp) :: end_proc_time,start_proc_time,end_total_time,&
         start_total_time
      real(dp),intent(out) :: Irad(:)
      real(dp),allocatable,dimension(:) :: acor,ecor,kappa,Ib
      logical,optional :: solve_RTE,quick_RTE
      logical :: compute_absorption,compute_emission,compute_kplanck,&
         compute_planck,compute_radiation,simplified

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('wsgg_sg_los_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Check if all arrays have the same size
      call assert(size(x).eq.size(T))
      call assert(size(T).eq.size(p))
      call assert(size(p).eq.size(xs,2))
      
      !Setting default input parameters
      nx = size(x); if (present(points)) nx = points      
      compute_radiation = .true.
      if (present(solve_RTE)) compute_radiation = solve_RTE
      simplified = .false.
      if (present(quick_RTE)) simplified = quick_RTE
      
      !Set flags for optional output arguments
      compute_absorption = .false.
      if (present(absorption))   compute_absorption = .true.
      compute_emission = .false.
      if (present(emission))     compute_emission = .true.
      compute_kplanck = .false.
      if (present(kappaPlanck))  compute_kplanck = .true.
      compute_planck = .false.
      if (present(Planck))       compute_planck = .true.

      !Surrogate names
      ns = number_of_species

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_sg_los_solution: Allocate arrays')
      allocate(acor(nx),stat=ierr)
      call CheckMemAlloc('acor',ierr)
      allocate(ecor(nx),stat=ierr)
      call CheckMemAlloc('ecor',ierr)
      allocate(kappa(nx),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(Ib(nx),stat=ierr)
      call CheckMemAlloc('Ib',ierr)

      !-----------------------------------------------------------------
      !Set up emission and absorption correction factors
      !-----------------------------------------------------------------
      call dprint('wsgg_sg_los_solution: &
                  &Set up emission and absorption correction factors')
      ecor = 1._dp; acor = 1._dp
      if (present(emission_correction)) ecor = emission_correction
      if (present(absorption_correction)) acor = absorption_correction

      !-----------------------------------------------------------------
      !Precompute local radiative properties
      !-----------------------------------------------------------------
      call dprint('wsgg_sg_los_solution: & 
                  &Precompute local radiative properties')
      prop_loop: do i=1,nx
         call wsgg_properties(kappa(i),a,jgas,T(i),xs(:,i),p(i),&
                              wsgg_spec)
      
         !Absorption coefficient of the gray gas
         kappa(i) = acor(i)*kappa(i)

         !RTE source term
         Ib(i) = (ecor(i)/(acor(i)+small))*a*sigrpi*(T(i)**4._dp)
         
      enddo prop_loop

      !-----------------------------------------------------------------
      !Solve radiation intensities
      !-----------------------------------------------------------------
      call dprint('wsgg_sg_los_solution: Solve radiation intensities')
      
      !Starting the processing time counter
      start_proc_time = get_wall_time()

      !Solve the RTE
      if (compute_radiation) &
         call exact_los_nonsct(x,kappa,Ib,I0,nx,Irad,simplified)

      !Stopping the processing time counter
      call cpu_time(end_proc_time)

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the optional properties
      !-----------------------------------------------------------------
      call dprint('wsgg_sg_los_solution: &
                  &Finish computing optional properties')
      do i=1,nx
         if (compute_emission)      emission(i) = kappa(i)*Ib(i)        !Emission                              
         if (compute_absorption)    absorption(i) = kappa(i)*Irad(i)    !Absorption
         if (compute_kplanck)       kappaPlanck(i) = kappa(i)           !Planck-mean absorption coefficient
         if (compute_planck)        Planck(i) = Ib(i)                   !Planck function
      enddo

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_sg_los_solution: Deallocate arrays')
      deallocate(acor,ecor,Ib,kappa)

      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('wsgg_sg_los_solution: Computational times')
      
      !Processing
      if (present(processing_time)) &
         processing_time = end_proc_time - start_proc_time

      !Computing (total)
      call cpu_time(end_total_time)
      if (present(total_time)) &
         total_time = end_total_time - start_total_time
   
   endsubroutine wsgg_sg_los_solution

   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !WSGG model for a non-scattering medium bounded by black walls
   !====================================================================
   subroutine wsgg_black_solution(wsgg_spec,flux,source,&
      processing_time,total_time,absorption,emission,kappaPlanck,&
      solve_RTE,total_loss)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,shutdown
      use constants, only: sigrpi
      use exact_solution, only: exact_1d_pp
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use physical_functions, only: bb_emission_frac
      use wsgg_functions, only: wsgg_number_gray_gases,wsgg_properties
      use wsgg_parameters, only : wsgg_bound_Ib,wsgg_Ib_lbound,&
         wsgg_Ib_ubound,wsgg_print_file,wsgg_print_properties
      character(*),intent(in) :: wsgg_spec
      character(5) :: gas_name
      character(200) :: msg
      integer :: i,j,jgas,k
      integer :: ngas,nl,nx,ny,nz
      integer :: ierr,uout
      logical,intent(in),optional :: solve_RTE
      logical :: compute_absorption,compute_emission,compute_planck,&
                 compute_radiation
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints),&
                                 total_loss
      real(dp) :: a,F,loss,loss_j,loss_omp
      real(dp) :: end_proc_time,end_total_time,start_proc_time,&
                  start_total_time
      real(dp),allocatable,dimension(:) :: flux_1d
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: kappa_sct,kIb,Ib,&
         source_j,source_omp,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: aIb,flux_j,flux_omp,&
                                                 kappa
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
            
      !Set flags for optional output arguments
      compute_radiation = .true.
      if (present(solve_RTE)) compute_radiation = solve_RTE
      compute_absorption=.false.
      if (present(absorption)) compute_absorption = .true.
      compute_emission=.false.
      if (present(emission)) compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.        

      !Surrogate names
      ngas = wsgg_number_gray_gases(wsgg_spec)                          !Number of gray gases of the model
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      
      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: Allocating arrays')
      allocate(aIb(nx,ny,nz,ngas+1),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(flux_1d(nx),stat=ierr)
      call CheckMemAlloc('flux_1d',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,ngas+1),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(source_j(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(source_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !-----------------------------------------------------------------
      !Prepare for dumping properties, if requested
      !-----------------------------------------------------------------
      if (wsgg_print_properties) then
         call dprint('wsgg_black_solution: &
                     &Prepare for dumping properties')
         
         !Prepare output unit
         if (trim(wsgg_print_file).eq.'null') &                         !Check if a file name was given                          
            call shutdown('wsgg_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(wsgg_print_file),&                    !Open file
              form='formatted',action='write')
         write(uout,'(a,(100(i3,:,",")))') '#0,',(i,i=1,2+2*(ngas+1))   !Write first line
         write(uout,'(3(a,:,","),100(i3,:,","))') '#',' ','j =',&       !Write header
                                             (jgas,jgas,jgas=1,ngas+1)
         write(uout,'(100(a,:,","))') '#x','y','z',&                    !Write subheader
                                    ('kappa_j','a_j',jgas=1,ngas+1)
         write(uout,'(100(a,:,","))') '#[m]','[m]','[m]',&              !Write subsubheader
                                    ('[1/m]','[-]',jgas=1,ngas+1)
      endif

      !-----------------------------------------------------------------
      !Computing properties
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: Compute properties')
      kIb = 0._dp                                                       !Zero out sum array
      do i=1,nx
         do j=1,ny
            do k=1,nz
               !Blackbody intensity
               if (wsgg_bound_Ib) then
                  F = bb_emission_frac(wsgg_Ib_lbound,T(i,j,k),.true.)-&
                      bb_emission_frac(wsgg_Ib_ubound,T(i,j,k),.true.)
               else
                  F = 1._dp
               endif
               Ib(i,j,k) = F*sigrpi*(T(i,j,k)**4._dp)
               
               !Absorption coefficient and weighting coefficient
               do jgas=1,ngas+1
                  call wsgg_properties(kappa(i,j,k,jgas),a,jgas,&
                     T(i,j,k),xs(:,i,j,k),p(i,j,k),wsgg_spec)
                  aIb(i,j,k,jgas) = a*Ib(i,j,k)                         !Emission term
                  kIb(i,j,k) = kIb(i,j,k) + &                           !Total RTE emission term
                               kappa(i,j,k,jgas)*aIb(i,j,k,jgas)
               enddo
            enddo
         enddo
      enddo
      
      !Misc parameters needed for the RTE solution subroutine
      kappa_sct = 0._dp                                                 !Scattering coefficient (null)
      phi = 1._dp                                                       !Phase function (isotropic)
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Set wall emissivities to unit for all wall cells
                         
      !-----------------------------------------------------------------
      !Main radiative transfer solution loop
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: Main solution loop started')
      flux = 0._dp; source = 0._dp; loss = 0._dp                        !Zero out sum arrays
      start_proc_time = get_wall_time()                                 !Start the processing time counter

      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_1d,flux_j,flux_omp,gas_name,jgas,loss_j,&
      !$OMP&           loss_omp,msg,source_j,source_omp)
      flux_omp = 0._dp; source_omp = 0._dp; loss_omp = 0._dp            !Zero out OMP sum arrays

      !$OMP DO
      do jgas=1,ngas+1
         !Longer debug information
         write(gas_name,'(i5)') jgas
         msg = 'wsgg_black_solution: gas_loop for gas '//&
               trim(adjustl(gas_name))
         call dprint(msg,3)

         !Solve the radiation field for gas jgas
         if (compute_radiation) then
            selectcase(rte_solution_method)
            case('exact-1d')
               call exact_1d_pp(x,kappa(:,1,1,jgas),aIb(:,1,1,jgas),nx,&
                                source_j(:,1,1),flux_1d)
               flux_j(1,:,1,1) = flux_1d
            case('FVM')                                                 !Finite volume method solution
               call fvm_solution(kappa(:,:,:,jgas),kappa(:,:,:,jgas),&
                  kappa_sct,aIb(:,:,:,jgas),phi,x_eps,y_eps,z_eps,&
                  flux_j,source_j,dont_iterate=.true.,total_loss=loss_j)

            case('PMC')                                                 !In this approach, the PMC method is executed
               which_spectral_mc = 'GG'                                 !  by-passing the spectral random number selection;
               call mc_solution(flux_j,source_j,fixed_Ib=&              !  i.e., the spectral integration is made in the
                  aIb(:,:,:,jgas),fixed_kappa=kappa(:,:,:,jgas))        !  same way as in the FVM

!           case('P1')
!              call p1_solution(kappa,aIb,flux_j,source_j)              !P1 method
            endselect
         endif
         
         !Update the heat flux, heat source and total loss
         flux_omp = flux_omp + flux_j
         source_omp = source_omp + source_j
         loss_omp = loss_omp + loss_j

      enddo
      !$OMP ENDDO
      
      !Finish the calculation of the total radiative quantities
      !$OMP CRITICAL
      flux = flux + flux_omp
      source = source + source_omp
      loss = loss + loss_omp
      !$OMP ENDCRITICAL
      
      !$OMP ENDPARALLEL
      end_proc_time = get_wall_time()                                   !End the processing time counter
      call dprint('wsgg_black_solution: Main solution loop ended')

      !-----------------------------------------------------------------
      !Dump of gas properties, if necessary
      !-----------------------------------------------------------------
      if (wsgg_print_properties) then
         call dprint('wsgg_black_solution: Dump of gas properties')
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  write(uout,'(1000(e26.15e3,:,","))')&
                     x(i),y(j),z(k),(kappa(i,j,k,jgas),&
                        aIb(i,j,k,jgas)/Ib(i,j,k),jgas=1,ngas+1)
               enddo
            enddo
         enddo
      endif

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: &
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
      if (present(total_loss)) total_loss = loss

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: Deallocating arrays')
      deallocate(aIb,Ib,kappa,kIb)
      deallocate(kappa_sct,phi,x_eps,y_eps,z_eps)
      
      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()       
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time
      
   endsubroutine wsgg_black_solution

   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !WSGG model for a non-scattering medium bounded by black walls
   !====================================================================
   subroutine wsgg_nongray_solution(wsgg_spec,flux,source,&
      processing_time,total_time,absorption,emission,kappaPlanck,&
      solve_RTE,total_loss)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,shutdown
      use constants, only: sigrpi
      use exact_solution, only: exact_1d_pp
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: number_of_surface_bands
!      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use physical_functions, only: bb_emission_frac
      use wsgg_functions, only: wsgg_number_gray_gases,wsgg_properties
      use wsgg_parameters, only : wsgg_print_file,wsgg_print_properties
      character(*),intent(in) :: wsgg_spec
      character(200) :: msg
      integer :: i,ieps,iomp,j,jgas,k
      integer :: neps,ngas,nomp,nl,nx,ny,nz
      integer :: ierr,uout
      integer,allocatable,dimension(:,:) :: iomp_array
      logical,intent(in),optional :: solve_RTE
      logical :: compute_absorption,compute_emission,compute_planck,&
                 compute_radiation
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints),&
                                 total_loss
      real(dp) :: a,F,loss,loss_j,loss_omp
      real(dp) :: end_proc_time,end_total_time,start_proc_time,&
                  start_total_time
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: kappa_sct,kIb,Ib,&
                                               source_j,source_omp
      real(dp),allocatable,dimension(:,:,:,:) :: flux_j,flux_omp,kappa,&
                                                 x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:,:) :: aIb
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
            
      !Set flags for optional output arguments
      compute_radiation = .true.
      if (present(solve_RTE)) compute_radiation = solve_RTE
      compute_absorption=.false.
      if (present(absorption)) compute_absorption = .true.
      compute_emission=.false.
      if (present(emission)) compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.        

      !Surrogate names
      neps = number_of_surface_bands                                    !Number of wall emissivity bands
      ngas = wsgg_number_gray_gases(wsgg_spec)                          !Number of gray gases of the model
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_nongray_solution: Allocating arrays')
      allocate(aIb(nx,ny,nz,ngas+1,neps),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,ngas+1),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(source_j(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(source_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(x_eps(2,ny,nz,neps),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz,neps),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny,neps),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !-----------------------------------------------------------------
      !Prepare for dumping properties, if requested
      !-----------------------------------------------------------------
      if (wsgg_print_properties) then
         call dprint('wsgg_nongray_solution: &
                     &Prepare for dumping properties')
         
         !Prepare output unit
         if (trim(wsgg_print_file).eq.'null') &                         !Check if a file name was given                          
            call shutdown('wsgg_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(wsgg_print_file),&                    !Open file
              form='formatted',action='write')
         write(uout,'(a,(100(i3,:,",")))') '#0,',(i,i=1,2+2*(ngas+1))   !Write first line
         write(uout,'(3(a,:,","),100(i3,:,","))') '#',' ','j =',&       !Write header
                                             (jgas,jgas,jgas=1,ngas+1)
         write(uout,'(100(a,:,","))') '#x','y','z',&                    !Write subheader
                                    ('kappa_j','a_j',jgas=1,ngas+1)
         write(uout,'(100(a,:,","))') '#[m]','[m]','[m]',&              !Write subsubheader
                                    ('[1/m]','[-]',jgas=1,ngas+1)
      endif

      !-----------------------------------------------------------------
      !Computing properties
      !-----------------------------------------------------------------
      call dprint('wsgg_nongray_solution: Compute properties')
      kIb = 0._dp                                                       !Zero out sum array
      do i=1,nx
         do j=1,ny
            do k=1,nz
               !Blackbody intensity
               Ib(i,j,k) = sigrpi*(T(i,j,k)**4._dp)
               
               !Absorption coefficient and weighting coefficient
               do jgas=1,ngas+1
                  call wsgg_properties(kappa(i,j,k,jgas),a,jgas,&
                     T(i,j,k),xs(:,i,j,k),p(i,j,k),wsgg_spec)
                  kIb(i,j,k) = kIb(i,j,k) + &                           !Total RTE emission term
                               kappa(i,j,k,jgas)*a*Ib(i,j,k)
                  do ieps=1,neps
                     F = bb_emission_frac(surface_bands(ieps,1),&       !Blackbody fraction within the band
                                          T(i,j,k),.true.) - &
                         bb_emission_frac(surface_bands(ieps,2),&
                                          T(i,j,k),.true.)
                     aIb(i,j,k,jgas,ieps) = a*F*Ib(i,j,k)               !Emission term
                  enddo
               enddo
            enddo
         enddo
      enddo
      
      !Misc parameters needed for the RTE solution subroutine
      kappa_sct = 0._dp                                                 !Scattering coefficient (null)
      phi = 1._dp                                                       !Phase function (isotropic)

      !-----------------------------------------------------------------
      !Define wall emissivities
      !-----------------------------------------------------------------
      call dprint('wsgg_nongray_solution: Define wall emissivities')
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Initialize arrays
      do ieps=1,neps
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  if (i.eq.1)  x_eps(1,j,k,ieps) = &
                     xmin_emissivity(j,k,ieps)
                  if (i.eq.nx) x_eps(2,j,k,ieps) = &
                     xmax_emissivity(j,k,ieps)
                  
                  if (j.eq.1)  y_eps(1,i,k,ieps) = &
                     ymin_emissivity(i,k,ieps)
                  if (j.eq.ny) y_eps(2,i,k,ieps) = &
                     ymax_emissivity(i,k,ieps)
                  
                  if (k.eq.1)  z_eps(1,i,j,ieps) = &
                     zmin_emissivity(i,j,ieps)
                  if (k.eq.nz) z_eps(2,i,j,ieps) = &
                     zmax_emissivity(i,j,ieps)
               enddo
            enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      !Mount iomp array
      !-----------------------------------------------------------------
      nomp = (ngas+1)*neps; allocate(iomp_array(nomp,2))
      iomp = 0
      do ieps=1,neps
         do jgas=1,ngas+1
            iomp = iomp + 1
            iomp_array(iomp,1:2) = (/ ieps, jgas /)
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Main radiative transfer solution loop
      !-----------------------------------------------------------------
      call dprint('wsgg_nongray_solution: Main solution loop started')
      flux = 0._dp; source = 0._dp; loss = 0._dp                        !Zero out sum arrays
      start_proc_time = get_wall_time()                                 !Start the processing time counter
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_j,flux_omp,ieps,jgas,loss_j,loss_omp,msg,&
      !$OMP&           source_j,source_omp)
      
      flux_omp = 0._dp; source_omp = 0._dp; loss_omp = 0._dp            !Zero out OMP sum arrays
      !$OMP DO
      do iomp=1,nomp
         !Define indexes
         ieps = iomp_array(iomp,1); jgas = iomp_array(iomp,2)

         !Solve the radiation field for gas jgas
         if (compute_radiation) then
            selectcase(rte_solution_method)
            case('FVM')                                                 !Finite volume method solution
               call fvm_solution(kappa(:,:,:,jgas),kappa(:,:,:,jgas),&
                  kappa_sct,aIb(:,:,:,jgas,ieps),phi,x_eps(:,:,:,ieps),&
                  y_eps(:,:,:,ieps),z_eps(:,:,:,ieps),flux_j,source_j,&
                  total_loss=loss_j)

!            case('PMC')                                                 !In this approach, the PMC method is executed
!               which_spectral_mc = 'GG'                                 !  by-passing the spectral random number selection;
!               call mc_solution(flux_j,source_j,fixed_Ib=&              !  i.e., the spectral integration is made in the
!                  aIb(:,:,:,jgas),fixed_kappa=kappa(:,:,:,jgas))        !  same way as in the FVM

!           case('P1')
!              call p1_solution(kappa,aIb,flux_j,source_j)              !P1 method
            endselect
         endif
         
         !Update the heat flux, heat source and total loss
         flux_omp = flux_omp + flux_j
         source_omp = source_omp + source_j
         loss_omp = loss_omp + loss_j

      enddo
      !$OMP ENDDO
      
      !Finish the calculation of the total radiative quantities
      !$OMP CRITICAL
      flux = flux + flux_omp
      source = source + source_omp
      loss = loss + loss_omp
      !$OMP ENDCRITICAL
      
      !$OMP ENDPARALLEL
      end_proc_time = get_wall_time()                                   !End the processing time counter
      call dprint('wsgg_nongray_solution: Main solution loop ended')

      !-----------------------------------------------------------------
      !Dump of gas properties, if necessary
      !-----------------------------------------------------------------
      if (wsgg_print_properties) then
         call dprint('wsgg_nongray_solution: Dump of gas properties')
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  write(uout,'(1000(e26.15e3,:,","))')&
                     x(i),y(j),z(k),(kappa(i,j,k,jgas),&
                        aIb(i,j,k,jgas,ieps)/Ib(i,j,k),jgas=1,ngas+1)
               enddo
            enddo
         enddo
      endif

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('wsgg_nongray_solution: &
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
      if (present(total_loss)) total_loss = loss

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_nongray_solution: Deallocating arrays')
      deallocate(aIb,Ib,kappa,kIb)
      deallocate(kappa_sct,phi,x_eps,y_eps,z_eps)
      
      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('wsgg_nongray_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()       
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time
      
   endsubroutine wsgg_nongray_solution   

   !====================================================================
   !Subroutine for the exact solution of the RTE for a one-dimensional
   !cylindrical, non-scattering medium
   !====================================================================
   subroutine wsgg_cyl_solution(x,T,p,xs,I0,source,flux,&
      processing_time,total_time,absorption,emission,kappaPlanck,points)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckMemAlloc,get_wall_time
      use wsgg_functions, only: wsgg_number_gray_gases,wsgg_properties
      use wsgg_parameters, only: which_wsgg_model
      use exact_solution, only: exact_1d_cyl
      integer :: i,ierr,j,jgas,k,ng,nx
      integer,intent(in),optional :: points
      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
      real(dp),intent(out) :: flux(:),source(:)
      real(dp),intent(out),optional :: processing_time,total_time,&
         absorption(:),emission(:),kappaPlanck(:)
      real(dp) :: a
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time,time_sum
      real(dp),allocatable,dimension(:) :: flux_j,Ib,kIb,source_j
      real(dp),allocatable,dimension(:,:) :: aIb,kappa
      logical :: compute_absorption=.false.,compute_emission=.false.,&
         compute_planck=.false.
start_total_time = I0 !Remove this
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Check if all arrays have the same size
      call assert(size(x).eq.size(T))
      call assert(size(T).eq.size(p))
      call assert(size(p).eq.size(xs,2))
      
      !Set flags for optional output arguments
      if (present(absorption)) compute_absorption = .true.
      if (present(emission)) compute_emission = .true.
      if (present(kappaPlanck)) compute_planck = .true.
      
      !Zeroing out processing time counter 
      !(necessary for the posterior sum over 
      !all instances of the RTE solution)
      time_sum = 0._dp
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      nx = size(x); if (present(points)) nx = points                    !Number of grid points
      
      !-----------------------------------------------------------------
      !Defining WSGG parameters
      !-----------------------------------------------------------------
      ng = wsgg_number_gray_gases(which_wsgg_model)                     !Number of gray gases of the model

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      allocate(aIb(nx,ng+1),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(flux_j(nx),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(Ib(nx),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kIb(nx),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(kappa(nx,ng+1),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(source_j(nx),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      
      !-----------------------------------------------------------------
      !Precomputing local radiative properties
      !-----------------------------------------------------------------
      j = 1; k = 1
      x_prop_loop: do i=1,nx
         Ib(i) = invpi*sigma*T(i)**4                                    !Local blackbody intensity
         kIb(i) = 0._dp
         gas_prop_loop: do jgas=1,ng+1
            call wsgg_properties(kappa(i,jgas),a,jgas,T(i),xs(:,i),&
                                 p(i),which_wsgg_model)
            aIb(i,jgas) = a*Ib(i)                                       !Local weighted blackbody
                                                                        !  intensity of gas jgas
            kIb(i) = kIb(i) + kappa(i,jgas)*aIb(i,jgas)                 !Total RTE emission term
         enddo gas_prop_loop
      enddo x_prop_loop
      
      !-----------------------------------------------------------------
      !Main solution loop   
      !-----------------------------------------------------------------
      !Starting the processing time counter
      start_proc_time = get_wall_time()

      !Prepare flux and source arrays
      flux = 0._dp
      source = 0._dp
      
      gas_solu_loop: do jgas=1,ng+1
         flux_j = 0._dp
         !Solving the radiation field for gas jgas
         call exact_1d_cyl(x,kappa(:,jgas),aIb(:,jgas),0._dp,nx,&       !Temporary
            source_j)               
      
         !Updating the heat flux and heat source
         flux = flux + flux_j
         source = source + source_j
      enddo gas_solu_loop
      
      !Stopping the processing time counter and adding to the total
      end_proc_time = get_wall_time()
      time_sum = time_sum + end_proc_time - start_proc_time
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      do i=1,nx
         if (compute_emission) emission(i) = fourpi*kIb(i)              !Total emission
         if (compute_absorption) absorption(i) = source(i) + &          !Total absorption
                                                 fourpi*kIb(i)
         if (compute_planck) kappaPlanck(i) = kIb(i)/Ib(i)              !Planck-mean absorption coefficient
      enddo
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      deallocate(aIb,Ib,kIb,kappa)
      deallocate(flux_j,source_j)
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      if (present(processing_time)) processing_time = time_sum          !Processing
      end_total_time = get_wall_time(); if (present(total_time)) &      !Computing (total)
                        total_time = end_total_time - start_total_time
      
   endsubroutine wsgg_cyl_solution   


























   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !superposition WSGG model for a non-scattering medium bounded by
   !black walls
   !====================================================================
   subroutine wsgg_sup_black_solution(flux,source,processing_time,&
      total_time,absorption,emission,kappaPlanck,solve_RTE,total_loss,&
      number_of_gases)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,shutdown
      use constants, only: sigrpi
      use exact_solution, only: exact_1d_pp
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: spec_name
      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use physical_functions, only: bb_emission_frac
      use wsgg_functions, only: wsgg_number_gray_gases,wsgg_properties
      use wsgg_parameters, only : wsgg_bound_Ib,wsgg_Ib_lbound,&
         wsgg_Ib_ubound,wsgg_print_file,wsgg_print_properties,&
         wsgg_improved_superposition,wsgg_superposition_theta,&
         wsgg_species_spec
      character(5) :: gas_name
      character(200) :: msg
      integer :: i,idsp,iispc,ispc,j,jjmix,jmix,k
      integer :: gas_counter,ngas,nl,nspc,nx,ny,nz
      integer :: ierr,uout
      integer,intent(out),optional :: number_of_gases
      integer,allocatable,dimension(:) :: dominant_spc,ngas_per_spc,&
                                          original_gas_index
      integer,allocatable,dimension(:,:) :: spc_gas_index
      logical,intent(in),optional :: solve_RTE
      logical :: advance_index,compute_absorption,compute_emission,&
                 compute_planck,compute_radiation
      logical,allocatable,dimension(:) :: skip_gas
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints),&
                                 total_loss
      real(dp) :: a,F,loss,loss_j,loss_omp
      real(dp) :: kappa_comp,kappa_mix,p_avg,T_avg,theta
      real(dp) :: end_proc_time,end_total_time,start_proc_time,&
                  start_total_time
      real(dp) :: dummy_dp
      real(dp),allocatable,dimension(:) :: flux_dummy,kappa_ref,xs_avg
      real(dp),allocatable,dimension(:,:) :: a_spc,kappa_spc,phi
      real(dp),allocatable,dimension(:,:,:) :: kappa_sct,kIb,Ib,&
         source_j,source_omp,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: aIb,flux_j,flux_omp,&
                                                 kappa
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('wsgg_sup_black_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
            
      !Set flags for optional output arguments
      compute_radiation = .true.
      if (present(solve_RTE)) compute_radiation = solve_RTE
      compute_absorption=.false.
      if (present(absorption)) compute_absorption = .true.
      compute_emission=.false.
      if (present(emission)) compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.        

      !Surrogate names
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nspc = number_of_species                                          !Number of participating species
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      theta = wsgg_superposition_theta                                  !Correction factor for the superposition method

      !Total number of gray gases
      allocate(ngas_per_spc(nspc),stat=ierr)
      call CheckMemAlloc('ngas_per_spc',ierr)
      ngas = 1
      do ispc=1,nspc
         ngas_per_spc(ispc) = &                                         !Number of gases in the correlation
            wsgg_number_gray_gases(wsgg_species_spec(ispc),&            !  for the individual species
                                   spec_name(ispc))
         ngas = ngas*(ngas_per_spc(ispc) + 1)                           !Total number of gray gases (standard superposition)
                                                                        !  or maximum number of gray gases (improved method)
      enddo

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_sup_black_solution: Allocating arrays')
      allocate(a_spc(1:ngas,1:nspc),stat=ierr)
      call CheckMemAlloc('a_spc',ierr)
      allocate(aIb(nx,ny,nz,ngas),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(dominant_spc(ngas),stat=ierr)
      call CheckMemAlloc('dominant_spc',ierr)
      allocate(flux_dummy(nx),stat=ierr)
      call CheckMemAlloc('flux_dummy',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,ngas),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_ref(1:nspc),stat=ierr)
      call CheckMemAlloc('kappa_ref',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kappa_spc(1:ngas,1:nspc),stat=ierr)
      call CheckMemAlloc('kappa_spc',ierr)
      allocate(kIb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(original_gas_index(ngas),stat=ierr)
      call CheckMemAlloc('original_gas_index',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(skip_gas(1:ngas),stat=ierr)
      call CheckMemAlloc('skip_gas',ierr)
      allocate(spc_gas_index(1:ngas,1:nspc),stat=ierr)
      call CheckMemAlloc('spc_gas_index',ierr)
      allocate(source_j(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_j',ierr)
      allocate(source_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_omp',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(xs_avg(1:nspc),stat=ierr)
      call CheckMemAlloc('xs_avg',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)
      
      !Compute the average composition
      if (wsgg_improved_superposition) &
         call get_average_composition(T_avg,xs_avg,p_avg,'T4')

      !-----------------------------------------------------------------
      !Prepare superposition WSGG arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_sup_black_solution: &
                  &Prepare superposition WSGG arrays')
      
      !Mount the spc_gas_index array 
      !(array storing the gray gas index of each individual
      !species ispc that corresponds to the overal gray gas jmix)
      spc_gas_index = 0                                                 !For jmix = 1, all species' gray gases are
                                                                        !  transparent windows
      do jmix=2,ngas
         do ispc=1,nspc        
            if (ispc.eq.1) then                                         !The first species is treated normally,
               spc_gas_index(jmix,ispc) = &                             !  with each new jmix corresponding to a
                  mod((spc_gas_index(jmix-1,ispc)+1),&                  !  new j for the species, and the counter
                      (ngas_per_spc(ispc)+1))                           !  resetting once j = J
            else                                                        !For all subsequent species, a new step in
               spc_gas_index(jmix,ispc) = spc_gas_index(jmix-1,ispc)    !  j is only taken (advance_index = .true)
               advance_index = .true.                                   !  if j = J @ jmix-1 for all previous species
               do iispc=1,ispc-1                                        !This loops checks if j = J @ jmix-1 for all
                  if (spc_gas_index(jmix-1,iispc).ne.&                  !  previous species. If that is the case, then
                      ngas_per_spc(iispc)) advance_index = .false.      !  advance_index = .true.
               enddo                                                    !And, if that is the case, than the j stepping
               if (advance_index) spc_gas_index(jmix,ispc) = &          !  is done, once again resetting if j = J for
                  mod((spc_gas_index(jmix,ispc)+1),&                    !  that species
                      (ngas_per_spc(ispc)+1)) 
            endif
         enddo 
      enddo
      
      !Mount the dominant_spc array
      !(this array informs which species is the dominant one, in the
      !framework of the improved superposition method; a value of -1
      !indicates that no species dominates)
      dominant_spc(1:ngas) = -1                                         !Initially, assume that no species dominate
      if (wsgg_improved_superposition) then
         do jmix=1,ngas
            !Compute reference values of kappa
            kappa_mix = 0._dp
            do ispc=1,nspc
               call wsgg_properties(kappa_ref(ispc),dummy_dp,&          !Compute the reference kappa for each species for
                  spc_gas_index(jmix,ispc),T_avg,xs_avg,p_avg,&         !  the current jmix
                  wsgg_species_spec(ispc),spec_name(ispc))
               kappa_mix = kappa_mix + kappa_ref(ispc)
!write(*,*) ispc,spc_gas_index(jmix,ispc),kappa_ref(ispc),dummy_dp
            enddo
            do ispc=1,nspc                                              !For each species, compute kappa of the mixture minus
               kappa_comp = kappa_mix - kappa_ref(ispc)                 !  the current species and compare to kappa of the 
               if (kappa_ref(ispc).gt.kappa_comp*theta) then            !  species, taking into account the safety factor theta
                  dominant_spc(jmix) = ispc                             !If the comparison is successful, flag that as the
                  exit                                                  !  dominant species for this jmix
               endif
            enddo
         enddo
      endif
!stop
      !Mount the skip_gas array
      !(this array informs if a jmix gas should be skipped when solving
      !the RTE in the improved superposition method to prevent repeated
      !calculations. This is done by sweeping through all gray gases 
      !larger than the current one and checking if the main condition 
      !of the improved method is not verified for any other combination 
      !of species' gray gases that include the current gas for the 
      !dominating species species
      do jmix=1,ngas
         idsp = dominant_spc(jmix)
         if (idsp.le.0) cycle
         skip_gas(jmix) = .false.                                       !In principle, do not skip this gas
         gas_check_loop: do jjmix=jmix+1,ngas                           !Loop over all remaining gray gases, but consider only
            if (spc_gas_index(jjmix,idsp).ne.&                          !  jmix correspondig to the same gray gas index as that
                spc_gas_index(jmix,idsp)) cycle gas_check_loop          !  of the current gray gas index of the dominant species
            if (dominant_spc(jjmix).eq.idsp) then                       !If, for any of the remaining gases, the current dominant
               skip_gas(jmix) = .true.                                  !  species dominates once again, we should flag jmix to
               exit gas_check_loop                                      !  be skipped
            endif
         enddo gas_check_loop
      enddo

      !-----------------------------------------------------------------
      !Computing properties
      !-----------------------------------------------------------------
      call dprint('wsgg_sup_black_solution: Compute properties')
      kIb = 0._dp                                                       !Zero out sum array
      do i=1,nx
         do j=1,ny
            do k=1,nz
               !Blackbody intensity
               if (wsgg_bound_Ib) then
                  F = bb_emission_frac(wsgg_Ib_lbound,T(i,j,k),.true.)-&
                      bb_emission_frac(wsgg_Ib_ubound,T(i,j,k),.true.)
               else
                  F = 1._dp
               endif
               Ib(i,j,k) = F*sigrpi*(T(i,j,k)**4._dp)

               !Absorption coefficient and weighting coefficient of
               !each individual species' gray gas
               do jmix=1,ngas
                  do ispc=1,nspc
                     call wsgg_properties(kappa_spc(jmix,ispc),&
                        a_spc(jmix,ispc),spc_gas_index(jmix,ispc),&
                        T(i,j,k),xs(:,i,j,k),p(i,j,k),&
                        wsgg_species_spec(ispc),spec_name(ispc))
                  enddo
               enddo

               !Actual absorption coefficient and weighting coefficient
               gas_counter = 0
               do jmix=1,ngas
                  idsp = dominant_spc(jmix)
                  if (idsp.lt.0) then                                   !Standard superposition
                     gas_counter = gas_counter + 1
                     original_gas_index(gas_counter) = jmix
                     kappa(i,j,k,gas_counter) = sum(kappa_spc(jmix,:))
                     a = product(a_spc(jmix,:))
                  else                                                  !Improved superposition
                     if (skip_gas(jmix)) cycle                          !Skip the jmix if the jmix is to be skipped
                     gas_counter = gas_counter + 1
                     original_gas_index(gas_counter) = jmix
                     kappa(i,j,k,gas_counter) = kappa_spc(jmix,idsp)    !kappa for the gas is equal to kappa of the
                                                                        !   dominant species
                     a = product(a_spc(jmix,:))                         !Initialize a as the value of a                                  
                     do jjmix=1,jmix-1                                  !To compute a for the gas, loop over all previous jmix
                        if (dominant_spc(jjmix).ne.idsp) cycle          !  skipping jmix for which the dominant species is not 
                                                                        !  equal to the current dominant species, and
                        if (spc_gas_index(jjmix,idsp).ne.&              !  jmix for which the dominant species' gray gas is 
                           spc_gas_index(jmix,idsp)) cycle              !  not equal to the current dominant species' gray gas
                        a = a + product(a_spc(jjmix,:))
                     enddo
                  endif                     
                  aIb(i,j,k,gas_counter) = a*Ib(i,j,k)                  !Emission term
                  kIb(i,j,k) = kIb(i,j,k) + kappa(i,j,k,gas_counter)*&  !Total RTE emission term
                               aIb(i,j,k,gas_counter)
               enddo
            enddo
         enddo
      enddo

      !Update the actual total number of gray gases
      ngas = gas_counter
      if (present(number_of_gases)) number_of_gases = ngas

      !Misc parameters needed for the RTE solution subroutine
      kappa_sct = 0._dp                                                 !Scattering coefficient (null)
      phi = 1._dp                                                       !Phase function (isotropic)
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Set wall emissivities to unit for all wall cells

      !-----------------------------------------------------------------
      !Main radiative transfer solution loop
      !-----------------------------------------------------------------
      call dprint('wsgg_sup_black_solution: Main solution loop started')
      flux = 0._dp; source = 0._dp; loss = 0._dp                        !Zero out sum arrays
      start_proc_time = get_wall_time()                                 !Start the processing time counter
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_dummy,flux_j,flux_omp,jmix,loss_j,loss_omp,&
      !$OMP&           msg,source_j,source_omp)

      flux_omp = 0._dp; source_omp = 0._dp; loss_omp = 0._dp            !Zero out OMP sum arrays
      !$OMP DO
      do jmix=1,ngas
         !Longer debug information
         write(gas_name,'(i5)') original_gas_index(jmix)
         msg = 'wsgg_sup_black_solution: gas_loop for gas '//&
               trim(adjustl(gas_name))
         call dprint(msg,3)

         !Solve the radiation field for gas jgas
         if (compute_radiation) then
            selectcase(rte_solution_method)
            case('exact-1d')
               call exact_1d_pp(x,kappa(:,1,1,jmix),aIb(:,1,1,jmix),nx,&
                                source_j(:,1,1),flux_dummy)
               flux_j(1,:,1,1) = flux_dummy
            case('FVM')                                                 !Finite volume method solution
               call fvm_solution(kappa(:,:,:,jmix),kappa(:,:,:,jmix),&
                  kappa_sct,aIb(:,:,:,jmix),phi,x_eps,y_eps,z_eps,&
                  flux_j,source_j,dont_iterate=.true.,total_loss=loss_j)
            case('PMC')                                                 !In this approach, the PMC method is executed
               which_spectral_mc = 'GG'                                 !  by-passing the spectral random number selection;
               call mc_solution(flux_j,source_j,fixed_Ib=&              !  i.e., the spectral integration is made in the
                  aIb(:,:,:,jmix),fixed_kappa=kappa(:,:,:,jmix))        !  same way as in the FVM

!           case('P1')
!              call p1_solution(kappa,aIb,flux_j,source_j)              !P1 method
            endselect
         endif

         !Update the heat flux, heat source and total loss
         flux_omp = flux_omp + flux_j
         source_omp = source_omp + source_j
         loss_omp = loss_omp + loss_j

      enddo
      !$OMP ENDDO

      !Finish the calculation of the total radiative quantities
      !$OMP CRITICAL
      flux = flux + flux_omp
      source = source + source_omp
      loss = loss + loss_omp
      !$OMP ENDCRITICAL

      !$OMP ENDPARALLEL
      end_proc_time = get_wall_time()                                   !End the processing time counter
      call dprint('wsgg_sup_black_solution: Main solution loop ended')

      !-----------------------------------------------------------------
      !Dump of gas properties, if necessary
      !-----------------------------------------------------------------
      if (wsgg_print_properties) then
         call dprint('wsgg_sup_black_solution: Dump of gas properties')

         !Prepare output unit
         if (trim(wsgg_print_file).eq.'null') &                         !Check if a file name was given                          
            call shutdown('wsgg_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(wsgg_print_file),&                    !Open file
              form='formatted',action='write')
         write(uout,'(a,(100(i3,:,",")))') '#0,',(i,i=1,2+2*(ngas))     !Write first line
         write(uout,'(3(a,:,","),100(i3,:,","))') '#',' ','j =',&       !Write header
            (original_gas_index(jmix),original_gas_index(jmix),&
            jmix=1,ngas)
         write(uout,'(100(a,:,","))') '#x','y','z',&                    !Write subheader
                                    ('kappa_j','a_j',jmix=1,ngas)
         write(uout,'(100(a,:,","))') '#[m]','[m]','[m]',&              !Write subsubheader
                                    ('[1/m]','[-]',jmix=1,ngas)

         !Main dump
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  write(uout,'(1000(e26.15e3,:,","))')&
                     x(i),y(j),z(k),(kappa(i,j,k,jmix),&
                        aIb(i,j,k,jmix)/Ib(i,j,k),jmix=1,ngas)
               enddo
            enddo
         enddo

         !Close unit
         close(uout)
      endif

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('wsgg_sup_black_solution: &
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
      if (present(total_loss)) total_loss = loss

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('wsgg_sup_black_solution: Deallocating arrays')
      deallocate(a_spc,aIb,Ib,kIb)
      deallocate(dominant_spc,original_gas_index,skip_gas,spc_gas_index)
      deallocate(flux_dummy,flux_j,flux_omp)
      deallocate(kappa,kappa_ref,kappa_sct,kappa_spc)
      deallocate(phi,x_eps,y_eps,z_eps)
      deallocate(source_j,source_omp)
      deallocate(xs_avg)

      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('wsgg_sup_black_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()       
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time

   endsubroutine wsgg_sup_black_solution

endmodule wsgg_routines
