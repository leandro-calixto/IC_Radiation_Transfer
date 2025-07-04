!##############################################
!Routines to solve the radiative heat transfer 
!via the NBCK method for different scenarios
!##############################################
module nbck_routines

   !====================================================================
   !Modules & Misc
   !====================================================================
   use nbck_parameters
   use nbck_functions
   implicit none
   
contains

   !====================================================================
   !Subroutine for the solution of the RTE along a line-of-sight for a
   !non-scattering medium. The code solves the RTE for a single narrow
   !band (of index inb) and single quadrature point (with CK value g),
   !and outputs the corresponding intensity
   !====================================================================
   subroutine nbck_qd_los_solution(x,T,p,xs,g,inb,I0,Irad,&
               absorption,emission,k,Planck,points,solve_RTE,quick_RTE,&
               emission_correction,absorption_correction)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckMemAlloc,dprint,&
         get_file_unit,get_wall_time,shutdown
      use constants, only: sigrpi
      use exact_solution, only: exact_los_nonsct 
      use physical_functions, only: bb_emission_frac
      use precision_parameters, only: dp,small
      implicit none
      integer,intent(in) :: inb
      integer :: i,ierr,nx
      integer,intent(in),optional :: points
      real(dp),intent(in) :: g,I0,p(:),T(:),x(:),xs(:,:)
      real(dp),intent(in),optional :: absorption_correction(:),&
                                      emission_correction(:)
      real(dp),intent(out) :: Irad(:)
      real(dp),intent(out),optional :: absorption(:),emission(:),k(:),&
                                       Planck(:)
      real(dp) :: F,nb_lbound,nb_ubound
      real(dp),allocatable,dimension(:) :: acor,ecor,Ib,kk
      logical,intent(in),optional :: quick_RTE,solve_RTE
      logical :: compute_absorption,compute_emission,compute_radiation,&
         compute_k,compute_planck,simplified

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('nbck_qd_los_solution: Preparatory procedures')

      !Define narrow band bounds
      nb_lbound = nbck_lbound(inb)
      nb_ubound = nbck_ubound(inb)

      !Set flags for optional output arguments
      simplified = .false.
      if (present(quick_RTE)) simplified = quick_RTE
      compute_radiation = .true.
      if (present(solve_RTE)) compute_radiation = solve_RTE
      compute_absorption = .false.
      if (present(absorption)) compute_absorption = .true.
      compute_emission = .false.
      if (present(emission)) compute_emission = .true.
      compute_k = .false.
      if (present(k)) compute_k = .true.
      compute_planck = .false.
      if (present(planck)) compute_planck = .true.

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('nbck_qd_los_solution: Surrogate names')
      nx = size(x); if  (present(points)) nx = points                   !Number of grid points
  
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('nbck_qd_los_solution: Allocate arrays')
      allocate(acor(nx),stat=ierr)
      call CheckMemAlloc('acor',ierr)
      allocate(ecor(nx),stat=ierr)
      call CheckMemAlloc('ecor',ierr)
      allocate(Ib(nx),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kk(nx),stat=ierr)
      call CheckMemAlloc('kk',ierr)

      !-----------------------------------------------------------------
      !Set up emission and absorption correction factors
      !-----------------------------------------------------------------
      call dprint('nbck_qd_los_solution: Set up correction factors')
      ecor = 1._dp; acor = 1._dp
      if (present(emission_correction))   ecor = emission_correction
      if (present(absorption_correction)) acor = absorption_correction

      !--------------------------------------------------------------
      !Computing local radiative properties
      !--------------------------------------------------------------
      call dprint('nbck_qd_los_solution: Compute local properties')
      grid_loop: do i=1,nx
         !Compute the local k
         kk(i) = acor(i)*nbck_mix_kg(inb,T(i),xs(i,:),p(i),g,'k')

         !Blackbody intensity
         F = bb_emission_frac(nb_lbound,T(i),.true.) - &
             bb_emission_frac(nb_ubound,T(i),.true.)
         Ib(i) = (ecor(i)/(acor(i)+small))*F*sigrpi*(T(i)**4._dp)

      enddo grid_loop

      !-----------------------------------------------------------------
      !Solve the radiation field
      !-----------------------------------------------------------------
      !Solving the radiation field for gas jgas
      call dprint('nbck_qd_los_solution: Solve radiation field')
      if (compute_radiation) &
         call exact_los_nonsct(x,kk,Ib,I0,nx,Irad,simplified)
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional, band-averaged properties
      !-----------------------------------------------------------------
      call dprint('nbck_qd_los_solution: &
                            &Finish calculation of optional properties')
      do i=1,nx
         if (compute_emission) emission(i) = kk(i)*Ib(i)                !Emission
         if (compute_absorption) absorption(i) = kk(i)*Irad(i)          !Absorption
         if (compute_planck) Planck(i) = Ib(i)                          !Planck function 
         if (compute_k) k(i) = kk(i)
      enddo
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('nbck_qd_los_solution: Deallocate arrays')
      deallocate(Ib,kk)

   endsubroutine nbck_qd_los_solution

   !====================================================================
   !Subroutine for the solution of the RTE along a line-of-sight for a
   !non-scattering medium. The code solves the RTE for a single narrow
   !band (of index inb), and outputs the corresponding intensity
   !====================================================================
   subroutine nbck_sb_los_solution(x,T,p,xs,I0,n_quad,inb,Irad,&
      processing_time,total_time,absorption,emission,kabs,kemi,Planck,&
      points,solve_RTE,quick_RTE)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckMemAlloc,dprint,&
                                get_file_unit,get_wall_time,shutdown
      use constants, only: sigrpi
      use exact_solution, only: exact_los_nonsct
      use physical_functions, only: bb_emission_frac
      use precision_parameters, only: dp,small
      implicit none
      character(10) :: qname
      character(200) :: msg
      integer,intent(in) :: inb,n_quad
      integer :: nx,nqd
      integer :: i,iqd
      integer :: ierr,uout
      integer,intent(in),optional :: points
      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
      real(dp),intent(out) :: Irad(:)
      real(dp),intent(out),optional :: processing_time,total_time,&
                     absorption(:),emission(:),kabs(:),kemi(:),Planck(:)
      real(dp) :: F,nb_lbound,nb_ubound
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time,time_sum
      real(dp),allocatable,dimension(:) :: Ib,Irad_qd,k,kI,kIb
      real(dp),allocatable,dimension(:) :: xquad,wquad
      real(dp),allocatable,dimension(:,:) :: k_qd
      logical,intent(in),optional :: quick_RTE,solve_RTE
      logical :: compute_absorption,compute_emission,compute_radiation,&
         compute_kabs,compute_kemi,compute_planck,simplified

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('nbck_sb_los_solution: Preparatory procedures')

      !Starting the total computing time counter
      start_total_time = get_wall_time()

      !Check if all arrays have the same size
      call assert(size(x).eq.size(T),'size(x)=size(T)')
      call assert(size(T).eq.size(p),'size(x)=size(p)')
      call assert(size(p).eq.size(xs,2),'size(x)=size(xs,2)')

      !Define narrow band bounds
      nb_lbound = nbck_lbound(inb)
      nb_ubound = nbck_ubound(inb)

      !Set flags for optional output arguments
      simplified = .false.
      if (present(quick_RTE)) simplified = quick_RTE
      compute_radiation = .true.
      if (present(solve_RTE)) compute_radiation = solve_RTE
      compute_absorption = .false.
      if (present(absorption)) compute_absorption = .true.
      compute_emission = .false.
      if (present(emission)) compute_emission = .true.
      compute_kabs = .false.
      if (present(kabs)) compute_kabs = .true.
      compute_kemi = .false.
      if (present(kemi)) compute_kemi = .true.
      compute_planck = .false.
      if (present(planck)) compute_planck = .true.

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('nbck_sb_los_solution: Surrogate names')
      nx = size(x); if  (present(points)) nx = points                   !Number of grid points
      nqd = n_quad                                                      !Number of quadrature points

      !-----------------------------------------------------------------
      !Prepare for dumping properties, if requested
      !-----------------------------------------------------------------
      if (nbck_print_properties) then
         call dprint('nbck_sb_los_solution: &
                                       &Prepare for dumping properties')

         !Allocate extra arrays
         allocate(k_qd(nx,nqd),stat=ierr)
         call CheckMemAlloc('k_qd',ierr)

         !Prepare output unit
         if (trim(nbck_print_file).eq.'null') &                         !Check if a file name was given                          
            call shutdown('nbck_sb_los_solution: &
                          &nbck_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(nbck_print_file),&                    !Open file
              form='formatted',action='write')
         write(uout,'(a,(100(i3,:,",")))') '#0,',(iqd,iqd=1,nqd)        !Write first line
         write(uout,'(3(a,:,","),100(i3,:,","))') '#i =',(iqd,iqd=0,nqd)!Write header
         write(uout,'(100(a,:,","))') '#x',('k_i',iqd=0,nqd)            !Write subheader
         write(uout,'(100(a,:,","))') '#[m]',('[1/m]',iqd=1,nqd)        !Write subsubheader
                                    
      endif
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('nbck_sb_los_solution: Allocate arrays')
      allocate(Ib(nx),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(Irad_qd(nx),stat=ierr)
      call CheckMemAlloc('Irad_qd',ierr)
      allocate(k(nx),stat=ierr)
      call CheckMemAlloc('k',ierr)
      allocate(kI(nx),stat=ierr)
      call CheckMemAlloc('kI',ierr)
      allocate(kIb(nx),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(xquad(nqd),stat=ierr)
      call CheckMemAlloc('xquad',ierr)
      allocate(wquad(nqd),stat=ierr)
      call CheckMemAlloc('wquad',ierr)

      !-----------------------------------------------------------------
      !Computing properties that do not depend on the gray gas
      !-----------------------------------------------------------------
      !Blackbody intensity
      call dprint('nbck_sb_los_solution: Compute blackbody intensity')
      do i=1,nx
         F = bb_emission_frac(nb_lbound,T(i),.true.) -&
             bb_emission_frac(nb_ubound,T(i),.true.)
         Ib(i) = F*sigrpi*(T(i)**4._dp)
      enddo

      !Quadrature points and weights
      call dprint('nbck_sb_los_solution: Define quadrature')
      call get_ck_quadrature(nbck_quadrature,xquad,wquad)

      !-----------------------------------------------------------------
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('nbck_sb_los_solution: Initialize sum arrays')
      Irad = 0._dp; kI = 0._dp; kIb = 0._dp
      
      quad_loop: do iqd=1,nqd
         !--------------------------------------------------------------
         !Preparatory procedures for the new gas
         !--------------------------------------------------------------
         !Longer debug information
         write(qname,'(i5)') iqd
         msg = 'nbck_sb_los_solution: quad_loop for quadrature '//&
            trim(adjustl(qname))//' started'
         call dprint(msg)

         !--------------------------------------------------------------
         !Computing local radiative properties
         !--------------------------------------------------------------
         call dprint('nbck_sb_los_solution: Compute local properties')
         grid_loop: do i=1,nx
            !Compute the local k
            k(i) = nbck_mix_kg(inb,T(i),xs(:,i),p(i),xquad(iqd),'k')

            !Store properties of each gray gas if requested
            if (nbck_print_properties) k_qd(i,iqd) = k(i)
         enddo grid_loop

         !--------------------------------------------------------------
         !Solve the radiation field
         !--------------------------------------------------------------
         !Starting the processing time counter
         start_proc_time = get_wall_time()
  
         !Solving the radiation field for gas jgas
         call dprint('nbck_sb_los_solution: Solve radiation field')
         if (compute_radiation) &
            call exact_los_nonsct(x,k,Ib,I0,nx,Irad_qd,simplified)
      
         !Stopping the processing time counter and adding to the total
         end_proc_time = get_wall_time()
         time_sum = time_sum + end_proc_time - start_proc_time

         !Updating the heat flux and heat source
         call dprint('nbck_sb_los_solution: &
                                        &Integrate spectral quantities')
         Irad = Irad + Irad_qd*wquad(iqd)
         kIb = kIb + k*Ib*wquad(iqd)
         kI = kI + k*Irad_qd*wquad(iqd)

      enddo quad_loop
      
      !-----------------------------------------------------------------
      !Dump gas properties, if necessary
      !-----------------------------------------------------------------
      if (nbck_print_properties) then
         call dprint('nbck_sb_los_solution: Dump gas properties')
         do i=1,nx
            write(uout,'(1000(e26.15e3,:,","))') &
               x(i),(k_qd(i,iqd),iqd=0,nqd)
         enddo
         
         !Close output unit
         close(uout)
      
         !Deallocate additional arrays
         deallocate(k_qd)
      endif
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional, band-averaged properties
      !-----------------------------------------------------------------
      call dprint('nbck_sb_los_solution: &
                            &Finish calculation of optional properties')
      do i=1,nx
         if (compute_emission) emission(i) = kIb(i)                     !Emission
         if (compute_absorption) absorption(i) = kI(i)                  !Absorption
         if (compute_planck) Planck(i) = Ib(i)                          !Planck function 
         if (compute_kabs) kabs(i) = kI(i)/(Irad(i) + small)            !Incident-mean absorption coefficient
         if (compute_kemi) kemi(i) = kIb(i)/(Ib(i) + small)             !Planck-mean absorption coefficient
      enddo
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('nbck_sb_los_solution: Deallocate arrays')
      deallocate(Ib,Irad_qd,k,kI,kIb)
      deallocate(xquad,wquad)
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      call dprint('nbck_sb_los_solution: Compute total times')
      if (present(processing_time)) processing_time = time_sum          !Processing
      end_total_time = get_wall_time(); if (present(total_time)) &      !Computing (total)
                        total_time = end_total_time - start_total_time
      
   endsubroutine nbck_sb_los_solution

   !====================================================================
   !Subroutine for the solution of the RTE along a line-of-sight for a
   !non-scattering medium. The code solves the RTE for the entire 
   !spectrum, and outputs the corresponding total intensity
   !====================================================================
   subroutine nbck_los_solution(x,T,p,xs,I0,n_quad,Irad,&
      processing_time,total_time,absorption,emission,kappaPlanck,&
      points,quick_RTE,read_premixed_nbck,write_premixed_nbck)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckMemAlloc,dprint,&
                                get_file_unit,get_wall_time,shutdown,print_to_prompt
      use constants, only: sigrpi
      use exact_solution, only: exact_los_nonsct
      use global_parameters, only: id_co2,id_h2o,number_of_species
      use omp_lib
      use physical_functions, only: bb_emission_frac
      use precision_parameters, only: dp,small
      implicit none
      character(10) :: bname,cname
      character(200) :: msg,premixed_file
      integer,intent(in) :: n_quad
      integer :: counter,i,ierr,inb,nnb,nx
      integer,intent(in),optional :: points
      logical,intent(in),optional :: quick_RTE,read_premixed_nbck,&
                                     write_premixed_nbck
      logical :: compute_absorption,compute_emission,compute_planck,&
                 read_premix,simplified,write_premix
      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
      real(dp),intent(out) :: Irad(:)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                absorption(:),emission(:),kappaPlanck(:)
      real(dp) :: Ib,Tmax,Tmin
      real(dp) :: end_proc_time,end_total_time,&
                  start_proc_time,start_total_time
      real(dp),allocatable,dimension(:) :: Irad_eta,Irad_omp,kI,kI_eta,&
         kI_omp,kIb,kIb_eta,kIb_omp,xsmax,xsmin

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('nbck_los_solution: Preparatory procedures')

      !Starting the total computing time counter
      start_total_time = get_wall_time()

      !Check if all arrays have the same size
      call assert(size(x).eq.size(T),'size(x)=size(T)')
      call assert(size(T).eq.size(p),'size(x)=size(p)')
      call assert(size(p).eq.size(xs,2),'size(x)=size(xs,2)')

      !Set flags for optional output arguments
      simplified = .false.
      if (present(quick_RTE)) simplified = quick_RTE
      compute_absorption=.false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission=.false.
      if (present(emission))    compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.
      read_premix = .false.
      if (present(read_premixed_nbck)) read_premix = read_premixed_nbck
      write_premix = .false.
      if (present(write_premixed_nbck)) &
         write_premix = write_premixed_nbck

      !-----------------------------------------------------------------
      !Load NBCK data
      !-----------------------------------------------------------------
      call dprint('nbck_los_solution: Load NBCK data')
      call load_nbck_data

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('nbck_los_solution: Surrogate names')
      nx = size(x); if  (present(points)) nx = points                   !Number of grid points
      nnb = number_nbck_bands                                           !Number of narrow bands

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('nbck_los_solution: Allocating arrays')
      allocate(Irad_eta(nx),stat=ierr)
      call CheckMemAlloc('Irad_eta',ierr)
      allocate(Irad_omp(nx),stat=ierr)
      call CheckMemAlloc('Irad_omp',ierr)
      allocate(kI(nx),stat=ierr)
      call CheckMemAlloc('kI',ierr)
      allocate(kI_eta(nx),stat=ierr)
      call CheckMemAlloc('kI_eta',ierr)
      allocate(kI_omp(nx),stat=ierr)
      call CheckMemAlloc('kI_omp',ierr)
      allocate(kIb(nx),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(kIb_eta(nx),stat=ierr)
      call CheckMemAlloc('kIb_eta',ierr)
      allocate(kIb_omp(nx),stat=ierr)
      call CheckMemAlloc('kIb_omp',ierr)
      
      !-----------------------------------------------------------------
      !Prepare the premixing of the NBCK distributions, if need be
      !-----------------------------------------------------------------
      if (trim(nbck_mixing_method).eq.'premixed') then
         call dprint('nbck_los_solution: NBCK premixing')

         if (nbck_mix_from_file) then
            call load_premix_nbck_data
         else
            !Allocate auxiliary arrays
            allocate(xsmax(number_of_species),stat=ierr)
            call CheckMemAlloc('xsmax',ierr)
            allocate(xsmin(number_of_species),stat=ierr)
            call CheckMemAlloc('xsmin',ierr)

            !Compute maxima and minima
            Tmax = maxval(T); Tmin = minval(T)
            xsmax(id_h2o) = maxval(xs(id_h2o,:))
            xsmin(id_h2o) = minval(xs(id_h2o,:))
            xsmax(id_co2) = maxval(xs(id_co2,:))
            xsmin(id_co2) = minval(xs(id_co2,:))
   
            !Prepare for the premixing
            call prepare_premixed_nbck(Tmin,Tmax,xsmin,xsmax)

            !Premixing
            premixed_file = trim(nbck_prefix)//trim(nbck_premix_file)
            if (write_premix) then
               call premix_nbck_data(write_to_file=premixed_file)
            elseif (read_premix) then
               call premix_nbck_data(read_from_file=premixed_file)
            else
               call premix_nbck_data
            endif

            !Deallocate auxiliary arrays
            deallocate(xsmax,xsmin)
         endif
      endif

      !-----------------------------------------------------------------
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('nbck_los_solution: Initialize sum arrays')
      Irad = 0._dp; kI = 0._dp; kIb = 0._dp; counter = 0
      
      !Start the processing time counter 
      start_proc_time = get_wall_time()  

      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(bname,Irad_eta,Irad_omp,kI_eta,kI_omp,kIb_eta,&
      !$OMP&           kIb_omp,msg)

      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      Irad_omp = 0._dp; kI_omp = 0._dp; kIb_omp = 0._dp

      !$OMP DO
      band_loop: do inb=1,nnb
         !--------------------------------------------------------------
         !Preparatory procedures for the new band
         !--------------------------------------------------------------
         !Longer debug information
         !$OMP CRITICAL
         counter = counter + 1
         write(bname,'(i5)') inb
         write(cname,'(f8.3)') real(counter,dp)*100._dp/real(nnb,dp)
         msg = 'nbck_los_solution: solving for narrow band '//&
            trim(adjustl(bname))//' ('//trim(adjustl(cname))//&
            '% concluded)'
         call print_to_prompt(msg,3)
         !$OMP ENDCRITICAL

         !--------------------------------------------------------------
         !Solve RTE for the band
         !--------------------------------------------------------------
         call dprint('nbck_los_solution: Solve the RTE')
         call nbck_sb_los_solution(x,T,p,xs,I0,n_quad,inb,Irad_eta,&
            absorption=kI_eta,emission=kIb_eta,points=nx,&
            quick_RTE=simplified)

         !Updating the heat flux and heat source
         call dprint('nbck_sb_los_solution: &
                                        &Integrate spectral quantities')
         Irad_omp = Irad_omp + Irad_eta
         kIb_omp = kIb_omp + kIb_eta
         kI_omp = kI_omp + kI_eta

      enddo band_loop
      !$OMP ENDDO
      
      !$OMP CRITICAL
      Irad = Irad + Irad_omp
      kIb = kIb + kIb_omp
      kI = kI + kI_omp
      !$OMP ENDCRITICAL
      !$OMP ENDPARALLEL
      
      !Stop the processing time counter
      end_proc_time = get_wall_time()
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('nbck_los_solution: &
                            &Finish calculation of optional properties')
      do i=1,nx
         Ib = sigrpi*T(i)**4
         if (compute_emission) emission(i) = kIb(i)                     !Total emission
         if (compute_absorption) absorption(i) = kI(i)                  !Total absorption
         if (compute_planck) kappaPlanck(i) = kIb(i)/Ib                 !Planck-mean absorption coefficient
      enddo
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('nbck_los_solution: Deallocate arrays')
      deallocate(Irad_eta,Irad_omp,kI,kI_eta,kI_omp,kIb,kIb_eta,kIb_omp)
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      call dprint('nbck_los_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time(); if (present(total_time)) &      !Computing (total)
                        total_time = end_total_time - start_total_time

   endsubroutine nbck_los_solution

   !====================================================================
   !Subroutine for the solution of the RTE via the NBCK method for a
   !non-scattering medium bounded by black walls
   !====================================================================
   subroutine nbck_black_solution(n_quad,flux,source,processing_time,&
                  total_time,absorption,emission,kappaPlanck,&
                  total_loss,read_premixed_nbck,write_premixed_nbck)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,CheckMemAlloc,dprint,&
                                get_wall_time,shutdown
      use constants, only: fourpi,sigrpi
      use exact_solution, only: exact_1d_pp,exact_los_nonsct
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: id_co2,id_h2o,number_of_species,&
                                   rte_solution_method
      use mesh
      use omp_lib                                                       !Parallelization
      use physical_functions, only: bb_emission_frac
      implicit none
      character(5) :: bname,qname
      character(200) :: cname,msg,premixed_file
      integer,intent(in) :: n_quad
      integer :: nl,nnb,nqd,nx,ny,nz
      integer :: counter,i,ierr,inb,iqd,j,k
      logical,optional,intent(in) :: &
         read_premixed_nbck,write_premixed_nbck
      logical :: compute_absorption,compute_emission,compute_planck,&
                 read_premix,write_premix
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
         total_loss,absorption(xpoints,ypoints,zpoints),&
         emission(xpoints,ypoints,zpoints),&
         kappaPlanck(xpoints,ypoints,zpoints)
      real(dp) :: F,Ib_total,nb_lbound,nb_ubound,Tmax,Tmin
      real(dp) :: loss,loss_omp,loss_qd
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp),allocatable,dimension(:) :: xquad,xsmax,xsmin,wquad
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: Ib,kappa,kappa_sct,&
         kIb,kIb_nb,kIb_omp,source_omp,source_qd,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:) :: source_nb
      real(dp),allocatable,dimension(:,:,:,:) :: flux_omp,flux_qd

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: Preparatory procedures')

      !Starting the total computing time counter
      start_total_time = get_wall_time()

      !Set up optional flags
      compute_absorption = .false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission = .false.
      if (present(emission))    compute_emission = .true.
      compute_planck = .false.
      if (present(kappaPlanck)) compute_planck = .true.
      read_premix = .false.
      if (present(read_premixed_nbck)) read_premix = read_premixed_nbck
      write_premix = .false.
      if (present(write_premixed_nbck)) &
         write_premix = write_premixed_nbck

      !-----------------------------------------------------------------
      !Load NBCK data
      !-----------------------------------------------------------------
      if (.not.nbck_mix_from_file) then
         call dprint('nbck_black_solution: Load NBCK data')
         call load_nbck_data
      endif
      
      !-----------------------------------------------------------------
      !Prepare the premixing of the NBCK distributions, if need be
      !-----------------------------------------------------------------
      if (trim(nbck_mixing_method).eq.'premixed') then
         call dprint('nbck_black_solution: NBCK premixing')
         
         if (nbck_mix_from_file) then
            call load_premix_nbck_data
         else
            !Allocate auxiliary arrays
            allocate(xsmax(number_of_species),stat=ierr)
            call CheckMemAlloc('xsmax',ierr)
            allocate(xsmin(number_of_species),stat=ierr)
            call CheckMemAlloc('xsmin',ierr)

            !Compute maxima and minima
            Tmax = maxval(T(:,:,:)); Tmin = minval(T(:,:,:))
            xsmax(id_h2o) = maxval(xs(id_h2o,:,:,:))
            xsmin(id_h2o) = minval(xs(id_h2o,:,:,:))
            xsmax(id_co2) = maxval(xs(id_co2,:,:,:))
            xsmin(id_co2) = minval(xs(id_co2,:,:,:))
   
            !Prepare for the premixing
            call prepare_premixed_nbck(Tmin,Tmax,xsmin,xsmax)

            !Premixing
            premixed_file = trim(nbck_prefix)//trim(nbck_premix_file)
            if (write_premix) then
               call premix_nbck_data(write_to_file=premixed_file)
            elseif (read_premix) then
               call premix_nbck_data(read_from_file=premixed_file)
            else
               call premix_nbck_data
            endif

            !Deallocate auxiliary arrays
            deallocate(xsmax,xsmin)
         endif
      endif

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: Surrogate names')
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nnb = number_nbck_bands                                           !Number of narrow bands
      nqd = n_quad                                                      !Number of quadrature points

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: Allocate arrays')
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(flux_qd(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_qd',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(kIb_nb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb_nb',ierr)
      allocate(kIb_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb_omp',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(source_nb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_nb',ierr)
      allocate(source_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_omp',ierr)
      allocate(source_qd(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_qd',ierr)
      allocate(xquad(nqd),stat=ierr)
      call CheckMemAlloc('xquad',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)
      allocate(wquad(nqd),stat=ierr)
      call CheckMemAlloc('wquad',ierr)

      !-----------------------------------------------------------------
      !Fixed properties
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: Fixed properties')
      kappa_sct = 0._dp                                                 !Scattering coefficient is null
      phi = 1._dp                                                       !Phase function (isotropic)
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Set wall emissivities to unit for all wall cells

      !-----------------------------------------------------------------
      !Quadrature points and weights
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: Define quadrature')
      call get_ck_quadrature(nbck_quadrature,xquad,wquad)

      !-----------------------------------------------------------------
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: Initialize sum arrays')
      flux = 0._dp; source = 0._dp; kIb = 0._dp; loss = 0._dp
      counter = 0
      
      !Start the processing time counter 
      start_proc_time = get_wall_time()  
      
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(bname,cname,end_proc_time,F,flux_omp,flux_qd,Ib,&
      !$OMP&           kappa,kIb_nb,kIb_omp,loss_omp,loss_qd,msg,&
      !$OMP&           nb_lbound,nb_ubound,qname,source_nb,source_omp,&
      !$OMP&           source_qd)

      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      flux_omp = 0._dp; source_omp = 0._dp
      kIb_omp = 0._dp; loss_omp = 0._dp
      kIb_nb = 0._dp; source_nb = 0._dp

      !$OMP DO
      band_loop: do inb=1,nnb
         !--------------------------------------------------------------
         !Preparatory procedures for the new band
         !--------------------------------------------------------------
         !Longer debug information
         !$OMP CRITICAL
         counter = counter + 1
         write(bname,'(i5)') inb
         write(cname,'(f8.3)') real(counter,dp)*100._dp/real(nnb,dp)
         msg = 'nbck_black_solution: solving for narrow band '//&
            trim(adjustl(bname))//' ('//trim(adjustl(cname))//&
            '% concluded)'
         call dprint(msg,3)
         !$OMP ENDCRITICAL
         
         !Narrow band bounds
         nb_lbound = nbck_lbound(inb)
         nb_ubound = nbck_ubound(inb)

         !--------------------------------------------------------------
         !Compute the local band-integrated blackbody intensity
         !--------------------------------------------------------------
         call dprint('nbck_black_solution: &
                      &Compute the band-integrated blackbody intensity')
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  F = bb_emission_frac(nb_lbound,T(i,j,k),.true.) -&
                  bb_emission_frac(nb_ubound,T(i,j,k),.true.)
                  Ib(i,j,k) = F*sigrpi*(T(i,j,k)**4._dp)
               enddo
            enddo
         enddo

         quad_loop: do iqd=1,nqd
            !-----------------------------------------------------------
            !Preparatory procedures for the new quadrature point
            !-----------------------------------------------------------
            !Longer debug information
            write(qname,'(i5)') iqd
            msg = 'nbck_black_solution: solving RTE for quadrature &
               & point'//trim(adjustl(qname))
            call dprint(msg)

            !-----------------------------------------------------------
            !Compute local k
            !-----------------------------------------------------------
            call dprint('nbck_black_solution: Compute local k')
            do i=1,nx
               do j=1,ny
                  do k=1,nz
!if (i.ne.(nx/2)) cycle
                     kappa(i,j,k) = nbck_mix_kg(inb,T(i,j,k),&
                                    xs(:,i,j,k),p(i,j,k),xquad(iqd),'k')
!write(*,*) i,kappa(i,j,k)
                  enddo
               enddo
            enddo

            !-----------------------------------------------------------
            !Solve the radiation field
            !-----------------------------------------------------------
            !Starting the processing time counter
            start_proc_time = get_wall_time()
  
            !Solving the radiation field for gas jgas
            call dprint('nbck_black_solution: &
                                            &Solve the radiation field')
            selectcase(rte_solution_method)
            case('exact-los')                                              !Solve for intensity (stored in source)
               call exact_los_nonsct(x,kappa(:,1,1),Ib(:,1,1),0._dp,nx,&
                                     source_qd(:,1,1))
            case('exact-1d')
               call exact_1d_pp(x,kappa(:,1,1),Ib(:,1,1),nx,&
                                source_qd(:,1,1),flux_qd(1,:,1,1))
            case('FVM')
               call fvm_solution(kappa,kappa,kappa_sct,Ib,phi,x_eps,&
                        y_eps,z_eps,flux_qd,source_qd,dont_iterate=&
                        .true.,total_loss=loss_qd)
            case default
               call shutdown('nbck_black_solution: rte_solution_method &
                             &not specified correctly')
            endselect

            !Updating the heat flux and heat source
            call dprint('nbck_black_solution: &
                                        &Integrate spectral quantities')
            kIb_nb = kIb_nb + kappa*Ib*wquad(iqd)
            source_nb = source_nb + source_qd*wquad(iqd)
            source_omp = source_omp + source_qd*wquad(iqd)
            flux_omp = flux_omp + flux_qd*wquad(iqd)
            loss_omp = loss_omp + loss_qd*wquad(iqd)
            kIb_omp = kIb_omp + kappa*Ib*wquad(iqd)

         enddo quad_loop
         
      enddo band_loop
      !$OMP ENDDO
      
      !-----------------------------------------------------------------
      !Finish the determination of the total radiative quantities
      !-----------------------------------------------------------------
      !$OMP CRITICAL
      source = source + source_omp
      flux = flux + flux_omp
      loss = loss + loss_omp
      kIb = kIb + kIb_omp
      !$OMP ENDCRITICAL
      !$OMP ENDPARALLEL
      
      !Stop the processing time counter
      end_proc_time = get_wall_time()
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: &
                            &Finish calculation of optional properties')
      if (trim(rte_solution_method).eq.'exact-los') kIb = kIb/fourpi
      do i=1,nx
         do j=1,ny
            do k=1,nz
               Ib_total = sigrpi*T(i,j,k)**4._dp                        !Total blackbody intensity
               if (compute_emission) &
                  emission(i,j,k) = fourpi*kIb(i,j,k)                   !Total emission
               if (compute_absorption) &
                  absorption(i,j,k) = source(i,j,k) + fourpi*kIb(i,j,k) !Total absorption
               if (compute_planck) &
                  kappaPlanck(i,j,k) = kIb(i,j,k)/Ib_total              !Planck-mean absorption coefficient
            enddo
         enddo
      enddo
      if (present(total_loss)) total_loss = loss
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: Deallocate arrays')
      deallocate(flux_qd,source_nb,source_qd)
      deallocate(Ib,kIb,kIb_nb)
      deallocate(kappa,kappa_sct,phi)
      deallocate(x_eps,y_eps,z_eps)
      deallocate(xquad,wquad)
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      call dprint('nbck_black_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time(); if (present(total_time)) &      !Computing (total)
                        total_time = end_total_time - start_total_time

   endsubroutine nbck_black_solution

endmodule nbck_routines
