!#######################################################################
!Routines to solve the radiative heat transfer 
!via LBL integration for different scenarios
!#######################################################################
module lbl_routines

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp,small
   use physical_functions, only: ib_function
   use lbl_functions
   use lbl_parameters
   implicit none
   
   contains

   !====================================================================
   !Subroutine for the solution of the RTE along a line-of-sight for a
   !non-scattering medium
   !====================================================================
   subroutine lbl_los_solution(x,T,p,xs,I0,Irad,processing_time,&
      total_time,absorption,emission,kappaPlanck,points,quick_RTE)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
         dprint,get_wall_time,print_to_prompt,shutdown
      use constants, only: sigrpi
      use exact_solution, only: exact_los_nonsct
      use global_parameters, only: id_soot,number_of_species
      use physical_functions, only: kappa_eta_soot
      implicit none
      character(200) :: msg
      integer,intent(in),optional :: points
      integer :: i,ierr,ilb,isp,nlb,nsp,nx
      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
      real(dp),intent(out) :: Irad(:)
      real(dp),intent(out),optional :: absorption(:),emission(:),&
         kappaPlanck(:),processing_time,total_time
      real(dp) :: Ib_total
      real(dp) :: end_proc_time,end_total_time,&
                  start_proc_time,start_total_time
      real(dp),allocatable,dimension(:) :: Ib,Irad_eta,Irad_omp,kappa,&
                                           kI,kI_omp,kIb,kIb_omp
      logical,optional :: quick_RTE
      logical :: compute_absorption,compute_emission,compute_kplanck,&
         simplified

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Check if all arrays have the same size
      call assert(size(x).eq.size(T),'size(x) = size(T)')
      call assert(size(T).eq.size(p),'size(T) = size(p)')
      call assert(size(p).eq.size(xs,1),'size(p) = size(xs)')
      
      !Setting default input parameters
      nx = size(x); if (present(points)) nx = points
      simplified = .false.
      if (present(quick_RTE)) simplified = quick_RTE
      
      !Set flags for optional output arguments
      compute_absorption = .false.
      if (present(absorption))   compute_absorption = .true.
      compute_emission = .false.
      if (present(emission))     compute_emission = .true.
      compute_kplanck = .false.
      if (present(kappaPlanck))  compute_kplanck = .true.

      !-----------------------------------------------------------------
      !Load the LBL data
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: Read the LBL data')
      call load_lbl_data

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: Surrogate names')
      nsp = number_of_species                                           !Number of species
      nlb = lbl_nlines                                                  !Number of lines in the LBL data files

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: Allocate arrays')
      allocate(Ib(nx),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(Irad_eta(nx),stat=ierr)
      call CheckMemAlloc('Irad_eta',ierr)
      allocate(Irad_omp(nx),stat=ierr)
      call CheckMemAlloc('Irad_omp',ierr)
      allocate(kappa(nx),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kI(nx),stat=ierr)
      call CheckMemAlloc('kI',ierr)
      allocate(kI_omp(nx),stat=ierr)
      call CheckMemAlloc('kI_omp',ierr)
      allocate(kIb(nx),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(kIb_omp(nx),stat=ierr)
      call CheckMemAlloc('kIb_omp',ierr)

      !-----------------------------------------------------------------
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: Initialize sum arrays')
      Irad = 0._dp; kI = 0._dp; kIb = 0._dp

      !Start the processing time counter 
      start_proc_time = get_wall_time()  

      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(Ib,Irad_eta,Irad_omp,kappa,kI_omp,kIb_omp,msg)
      
      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      Irad_omp = 0._dp; kI_omp = 0._dp; kIb_omp = 0._dp

      !$OMP DO
      lines_loop: do ilb=1,nlb
         !--------------------------------------------------------------
         !Prepare the new spectral line
         !--------------------------------------------------------------
         !Check if this eta is to be solved for
         if ((lbl_xeta(ilb).lt.lbl_eta_min).or.&
             (lbl_xeta(ilb).gt.lbl_eta_max)) cycle lines_loop

         !Print xeta if requested
         if (print_xeta) then
            write(msg,'(f12.4)') lbl_xeta(ilb)/100._dp
            msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
            call print_to_prompt(msg,3)
         endif

         !--------------------------------------------------------------
         !Determine the properties for this wavenumber
         !--------------------------------------------------------------
         call dprint('lbl_los_solution: Determine line properties',6)         
         do i=1,nx
            !Absorption coefficient
            kappa(i) = 0._dp
            do isp=1,nsp
               if (isp.eq.id_soot) then
                  kappa(i) = kappa(i) + &
                     kappa_eta_soot(lbl_xeta(ilb),xs(i,isp))
               else
                  kappa(i) = kappa(i) + &
                     lbl_kappa_func(xs(i,isp),T(i),p(i),isp,ilb)
               endif
            enddo

            !Spectral blackbory radiative intensity
            Ib(i) = Ib_function(T(i),lbl_xeta(ilb))

            !Integrate RTE emission term
            kIb_omp(i) = kIb_omp(i) + kappa(i)*Ib(i)*lbl_deta(ilb)

         enddo
         
         !--------------------------------------------------------------
         !Solve the radiation field
         !--------------------------------------------------------------
         call dprint('lbl_los_solution: Solve the radiation field',6)  
         call exact_los_nonsct(x,kappa,Ib,I0,nx,Irad_eta,simplified)
                           
         !Update the radiative intensity and absorption
         Irad_omp = Irad_omp + Irad_eta*lbl_deta(ilb)
         kI_omp = kI_omp + kappa*Irad_eta*lbl_deta(ilb)
 
      enddo lines_loop
      !$OMP ENDDO

      !-----------------------------------------------------------------
      !Finish the calculation of the total radiative quantities
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: &
           &Finish the calculation of the total radiative quantities',6)  
      !$OMP CRITICAL
      Irad = Irad + Irad_omp
      kI = kI + kI_omp
      kIb = kIb + kIb_omp
      !$OMP ENDCRITICAL
      !$OMP ENDPARALLEL
      
      !Stop the processing time counter
      end_proc_time = get_wall_time()

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: &
                        &Finish calculation of the optional properties')
      do i=1,nx
         Ib_total = sigrpi*T(i)**4._dp                                  !Total blackbody intensity
         if (compute_emission) emission(i) = kIb(i)                     !Total emission
         if (compute_absorption) absorption(i) = kI(i)                  !Total absorption
         if (compute_kplanck) kappaPlanck(i) = kIb(i)/Ib_total          !Planck-mean absorption coefficient
      enddo

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: Deallocate arrays')
      deallocate(Irad_eta,Irad_omp,Ib,kappa,kI,kI_omp,kIb,kIb_omp)

      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      call dprint('lbl_los_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()
      if (present(total_time)) &                                        !Computing (total)
                        total_time = end_total_time - start_total_time

   endsubroutine lbl_los_solution

   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !LBL method for a non-scattering medium bounded by black walls
   !====================================================================
   subroutine lbl_black_solution(flux,source,processing_time,&
                  total_time,absorption,emission,kappaPlanck,&
                  total_loss,sri_output,solve_RTE)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,CheckMemAlloc,dprint,&
                                get_wall_time,print_to_prompt,shutdown
      use constants, only: fourpi,invpi,sigma
      use exact_solution, only: exact_1d_pp,exact_los_nonsct
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: id_soot,number_of_species,&
                                   rte_solution_method
      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use omp_lib
      use physical_functions, only: kappa_eta_soot
      implicit none
      character(100) :: msg
      integer :: i,ierr,ilbl,isp,j,k
      integer :: nl,nlbl,nsp,nx,ny,nz
      integer :: icm,icp,jcm,jcp,kcm,kcp
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
         total_loss,absorption(xpoints,ypoints,zpoints),&
         emission(xpoints,ypoints,zpoints),&
         kappaPlanck(xpoints,ypoints,zpoints),sri_output(:,:)
      real(dp) :: dV,dx,dy,dz,intV(13),loss,loss_eta,loss_omp,magflux
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp) :: kappa_spc,Ib_total
      real(dp),allocatable,dimension(:) :: flux_dummy
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: Ib,kappa,kappa_sct,&
         kIb,kIb_omp,source_eta,source_omp,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: flux_eta,flux_omp
      logical,intent(in),optional :: solve_RTE
      logical :: compute_absorption,compute_emission,compute_planck,&
                 compute_radiation

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('lbl_black_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Check for repeated flags that should be unique
      call check_lbl_parameters
      
      !Set flags for optional output arguments
      compute_absorption=.false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission=.false.
      if (present(emission))    compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.
      compute_radiation=.true.
      if (present(solve_RTE))   compute_radiation = solve_RTE

      !-----------------------------------------------------------------
      !Load the LBL data
      !-----------------------------------------------------------------
      call dprint('lbl_black_solution: Read the LBL data')
      call load_lbl_data

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('lbl_black_solution: Surrogate names')
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nlbl = lbl_nlines                                                 !Number of lines in the LBL data files
      nsp = number_of_species                                           !Number of participating species
      
      !Upper and lower cell indexes
      if (present(sri_output)) then
         icm = 2; icp = nx-1
         jcm = 2; jcp = ny-1
         kcm = 2; kcp = nz-1
         if (one_d) jcm = 1; if (one_d) jcp = 1
         if (one_d.or.two_d) kcm = 1; if (one_d.or.two_d) kcp = 1
      endif
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('lbl_black_solution: Allocate arrays')
      allocate(flux_dummy(nx),stat=ierr)
      call CheckMemAlloc('flux_dummy',ierr)
      allocate(flux_eta(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_eta',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(kIb_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb_omp',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(source_eta(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_eta',ierr)
      allocate(source_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_omp',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !-----------------------------------------------------------------
      !Define fixed properties
      !-----------------------------------------------------------------
      call dprint('lbl_black_solution: Define fixed properties')
      kappa_sct = 0._dp                                                 !Scattering coefficient is null
      phi = 1._dp                                                       !Phase function (isotropic)
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Set wall emissivities to unit for all wall cells
      
      !-----------------------------------------------------------------
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('lbl_black_solution: Initialize sum arrays')
      flux = 0._dp; source = 0._dp; kIb = 0._dp; loss = 0._dp
      
      !Start the processing time counter 
      start_proc_time = get_wall_time()  

      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(dV,dx,dy,dz,flux_dummy,flux_eta,flux_omp,Ib,&
      !$OMP&           kappa,kappa_spc,intV,kIb_omp,loss_eta,loss_omp,&
      !$OMP&           magflux,source_eta,source_omp)
      
      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      flux_omp = 0._dp; source_omp = 0._dp
      kIb_omp = 0._dp; loss_omp = 0._dp
      
      !$OMP DO
      lines_loop: do ilbl=1,nlbl
         !--------------------------------------------------------------
         !Prepare the new spectral line
         !--------------------------------------------------------------
         !Check if this eta is to be solved for
         if ((lbl_xeta(ilbl).lt.lbl_eta_min).or.&
             (lbl_xeta(ilbl).gt.lbl_eta_max)) cycle lines_loop

         !Print xeta if requested
         if (print_xeta) then
            write(msg,'(f12.4)') lbl_xeta(ilbl)/100._dp
            msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
            call print_to_prompt(msg,3)
         endif

         !--------------------------------------------------------------
         !Determine the properties for this wavenumber
         !--------------------------------------------------------------
         call dprint('lbl_black_solution: Determine line properties',3)
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  !Absorption coefficient
                  kappa(i,j,k) = 0._dp
                  do isp=1,nsp
                     if (isp.eq.id_soot) then
                        kappa_spc = &
                           kappa_eta_soot(lbl_xeta(ilbl),xs(isp,i,j,k))
                     else                  
                        kappa_spc = lbl_kappa_func(xs(isp,i,j,k),&
                                             T(i,j,k),p(i,j,k),isp,ilbl)
                           
                     endif
                     kappa(i,j,k) = kappa(i,j,k) + kappa_spc
                  enddo

                  !Spectral blackbory radiative intensity
                  Ib(i,j,k) = Ib_function(T(i,j,k),lbl_xeta(ilbl))

                  !Integrate RTE emission term
                  kIb_omp(i,j,k) = kIb_omp(i,j,k) + &
                     kappa(i,j,k)*Ib(i,j,k)*lbl_deta(ilbl)

               enddo
            enddo
         enddo

         !--------------------------------------------------------------
         !Solve the radiation field
         !--------------------------------------------------------------
         call dprint('lbl_black_solution: Solve the radiation field',3)
         if (compute_radiation) then
            selectcase(rte_solution_method)
            case('exact-los')                                           !Solve for intensity (stored in source)
               call exact_los_nonsct(x,kappa(:,1,1),Ib(:,1,1),0._dp,nx,&
                                    source_eta(:,1,1))
            case('exact-1d')
               call exact_1d_pp(x,kappa(:,1,1),Ib(:,1,1),nx,&
                                source_eta(:,1,1),flux_dummy)
               flux_eta(1,:,1,1) = flux_dummy
            case('FVM')
               call fvm_solution(kappa,kappa,kappa_sct,Ib,phi,x_eps,&
                        y_eps,z_eps,flux_eta,source_eta,dont_iterate=&
                        .true.,total_loss=loss_eta)
            case('PMC')                                                 !In this approach, the PMC method is executed
               which_spectral_mc = 'GG'                                 !  by-passing the spectral random number selection;
               call mc_solution(flux_eta,source_eta,fixed_Ib=Ib,&       !  i.e., the spectral integration is made in the
                                fixed_kappa=kappa,total_loss=loss_eta)  !  same way as in the FVM
            case default
               call shutdown('lbl_black_solution: rte_solution_method &
                             &not specified correctly')
            endselect
         endif

         !Update the heat flux and heat source
         call dprint('lbl_black_solution: &
                     &Update the heat flux and heat source',3)
         flux_omp = flux_omp + flux_eta*lbl_deta(ilbl)
         source_omp = source_omp + source_eta*lbl_deta(ilbl)
         loss_omp = loss_omp + loss_eta*lbl_deta(ilbl)
 
         !--------------------------------------------------------------
         !Compute and store quantities for the SRI analyses, 
         !if requested
         !--------------------------------------------------------------
         !$OMP CRITICAL
         if (present(sri_output)) then     
            call dprint('lbl_black_solution: &
                        &Compute quantities for the SRI analyses',3)
            intV = 0._dp; intV(1) = lbl_xeta(ilbl)
            do i=icm,icp
               do j=jcm,jcp
                  do k=kcm,kcp
                     dx = xf(i) - xf(i-1)                               !Cell size: x direction
                     dy = 1._dp; if (.not.one_d) dy = yf(j) - yf(j-1)   !Cell size: y direction
                     dz = 1._dp; if (three_d)    dz = zf(k) - zf(k-1)   !Cell size: z direction
                     dV = dx*dy*dz                                      !Cell volume
                     intV(2) = intV(2) + dabs(source_eta(i,j,k))*dV
                     intV(3) = intV(3) + source_eta(i,j,k)*dV
                     intV(4) = intV(4) + dabs(flux_eta(1,i,j,k))*dV
                     intV(5) = intV(5) + flux_eta(1,i,j,k)*dV
                     intV(6) = intV(6) + dabs(flux_eta(2,i,j,k))*dV
                     intV(7) = intV(7) + flux_eta(2,i,j,k)*dV
                     intV(8) = intV(8) + dabs(flux_eta(3,i,j,k))*dV
                     intV(9) = intV(9) + flux_eta(3,i,j,k)*dV
                     magflux = dsqrt(flux_eta(1,i,j,k)**2 + &
                                     flux_eta(2,i,j,k)**2 + &
                                     flux_eta(3,i,j,k)**2)
                     intV(10) = intV(10) + magflux*dV
                     intV(11) = intV(11) + Ib(i,j,k)*dV
                     intV(12) = intV(12) + &
                        fourpi*kappa(i,j,k)*Ib(i,j,k)*dV
                     intV(13) = intV(13) + (source(i,j,k) + &
                        fourpi*kappa(i,j,k)*Ib(i,j,k))*dV
                  enddo
               enddo
            enddo
            sri_output(ilbl,:) = intV
         endif
         !$OMP ENDCRITICAL

      enddo lines_loop
      !$OMP ENDDO

      !-----------------------------------------------------------------
      !Finish the calculation of the total radiative quantities
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
      call dprint('lbl_black_solution: &
                  &Finish calculation of the optional properties')
      if (trim(rte_solution_method).eq.'exact-los') kIb = kIb/fourpi
      do i=1,nx
         do j=1,ny
            do k=1,nz
               Ib_total = sigma*invpi*T(i,j,k)**4._dp                   !Total blackbody intensity
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
      call dprint('lbl_black_solution: Deallocate arrays')
      deallocate(flux_eta,source_eta)
      deallocate(flux_omp,source_omp)
      deallocate(Ib,kIb,kIb_omp)
      deallocate(kappa,kappa_sct,phi)
      deallocate(x_eps,y_eps,z_eps)
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      call dprint('lbl_black_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()
      if (present(total_time)) &                                        !Computing (total)
                        total_time = end_total_time - start_total_time
      
   endsubroutine lbl_black_solution

   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !LBL method for a non-scattering medium bounded by black walls
   !====================================================================
   subroutine lbl_nongray_solution(flux,source,processing_time,&
                  total_time,absorption,emission,kappaPlanck,&
                  total_loss,sri_output)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,CheckMemAlloc,dprint,&
                                get_wall_time,print_to_prompt,shutdown
      use constants, only: fourpi,invpi,sigma
      use exact_solution, only: exact_1d_pp,exact_los_nonsct
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: number_of_species,&
         number_of_surface_bands,rte_solution_method
!      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use omp_lib
      implicit none
      character(100) :: msg
      integer :: i,ieps,ierr,ilbl,isp,j,k
      integer :: neps,nl,nlbl,nsp,nx,ny,nz
      integer :: icm,icp,jcm,jcp,kcm,kcp
      logical :: compute_absorption,compute_emission,compute_planck
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
         total_loss,absorption(xpoints,ypoints,zpoints),&
         emission(xpoints,ypoints,zpoints),&
         kappaPlanck(xpoints,ypoints,zpoints),sri_output(:,:)
      real(dp) :: dV,dx,dy,dz,intV(13),loss,loss_eta,loss_omp,magflux
      real(dp) :: end_proc_time,start_proc_time,&
                  end_total_time,start_total_time
      real(dp) :: Ib_total
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: Ib,kappa,kappa_sct,&
         kIb,kIb_omp,source_eta,source_omp
      real(dp),allocatable,dimension(:,:,:,:) :: flux_eta,flux_omp,&
         x_eps,y_eps,z_eps

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('lbl_nongray_solution: Preparatory procedures')
      
      !Starting the total computing time counter
      start_total_time = get_wall_time()
      
      !Check for repeated flags that should be unique
      call check_lbl_parameters
      
      !Set flags for optional output arguments
      compute_absorption=.false.
      if (present(absorption))  compute_absorption = .true.
      compute_emission=.false.
      if (present(emission))    compute_emission = .true.
      compute_planck=.false.
      if (present(kappaPlanck)) compute_planck = .true.

      !-----------------------------------------------------------------
      !Load the LBL data
      !-----------------------------------------------------------------
      call dprint('lbl_nongray_solution: Read the LBL data')
      call load_lbl_data

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('lbl_nongray_solution: Surrogate names')
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nlbl = lbl_nlines                                                 !Number of lines in the LBL data files
      nsp = number_of_species                                           !Number of participating species
      neps = number_of_surface_bands                                    !Number of wall emissivity bands

      !Upper and lower cell indexes
      if (present(sri_output)) then
         icm = 2; icp = nx-1
         jcm = 2; jcp = ny-1
         kcm = 2; kcp = nz-1
         if (one_d) jcm = 1; if (one_d) jcp = 1
         if (one_d.or.two_d) kcm = 1; if (one_d.or.two_d) kcp = 1
      endif

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('lbl_nongray_solution: Allocate arrays')
      allocate(flux_eta(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_eta',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa',ierr)
      allocate(kappa_sct(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kappa_sct',ierr)
      allocate(kIb(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb',ierr)
      allocate(kIb_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('kIb_omp',ierr)
      allocate(phi(nl,nl),stat=ierr)
      call CheckMemAlloc('phi',ierr)
      allocate(source_eta(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_eta',ierr)
      allocate(source_omp(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('source_omp',ierr)
      allocate(x_eps(2,ny,nz,neps),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz,neps),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny,neps),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)

      !-----------------------------------------------------------------
      !Define fixed properties
      !-----------------------------------------------------------------
      call dprint('lbl_nongray_solution: Define fixed properties')
      kappa_sct = 0._dp                                                 !Scattering coefficient is null
      phi = 1._dp                                                       !Phase function (isotropic)

      !-----------------------------------------------------------------
      !Define wall emissivities
      !-----------------------------------------------------------------
      call dprint('lbl_nongray_solution: Define wall emissivities')
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
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('lbl_nongray_solution: Initialize sum arrays')
      flux = 0._dp; source = 0._dp; kIb = 0._dp; loss = 0._dp
      
      !Start the processing time counter 
      start_proc_time = get_wall_time()  

      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(dV,dx,dy,dz,flux_eta,flux_omp,Ib,ieps,kappa,intV,&
      !$OMP&           kIb_omp,loss_eta,loss_omp,magflux,source_eta,&
      !$OMP&           source_omp)
      
      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      flux_omp = 0._dp; source_omp = 0._dp
      kIb_omp = 0._dp; loss_omp = 0._dp
      
      !$OMP DO
      lines_loop: do ilbl=1,nlbl
         !--------------------------------------------------------------
         !Prepare the new spectral line
         !--------------------------------------------------------------
         !Check if this eta is to be solved for
         if ((lbl_xeta(ilbl).lt.lbl_eta_min).or.&
             (lbl_xeta(ilbl).gt.lbl_eta_max)) cycle lines_loop

         !Print xeta if requested
         if (print_xeta) then
            write(msg,'(f12.4)') lbl_xeta(ilbl)/100._dp
            msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
            call print_to_prompt(msg,3)
         endif
         
         !Define which wall emissivity band to use
         do ieps=1,neps
            if ((lbl_xeta(ilbl).ge.surface_bands(ieps,1)).and.&
                (lbl_xeta(ilbl).le.surface_bands(ieps,2))) exit
         enddo
         if (ieps.gt.neps) cycle lines_loop

         !--------------------------------------------------------------
         !Determine the properties for this wavenumber
         !--------------------------------------------------------------
         call dprint('lbl_nongray_solution: &
                     &Determine line properties',3)
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  !Absorption coefficient
                  kappa(i,j,k) = 0._dp
                  do isp=1,nsp
                     kappa(i,j,k) = kappa(i,j,k) + &
                        lbl_kappa_func(xs(isp,i,j,k),T(i,j,k),&
                                       p(i,j,k),isp,ilbl)
                           
                  enddo

                  !Spectral blackbory radiative intensity
                  Ib(i,j,k) = Ib_function(T(i,j,k),lbl_xeta(ilbl))

                  !Integrate RTE emission term
                  kIb_omp(i,j,k) = kIb_omp(i,j,k) + &
                     kappa(i,j,k)*Ib(i,j,k)*lbl_deta(ilbl)

               enddo
            enddo
         enddo

         !--------------------------------------------------------------
         !Solve the radiation field
         !--------------------------------------------------------------
         call dprint('lbl_nongray_solution: &
                     &Solve the radiation field',3)
         selectcase(rte_solution_method)
!         case('exact-los')                                              !Solve for intensity (stored in source)
!            call exact_los_nonsct(x,kappa(:,1,1),Ib(:,1,1),0._dp,nx,&
!                                  source_eta(:,1,1))
!         case('exact-1d')
!            call exact_1d_pp(x,kappa(:,1,1),Ib(:,1,1),nx,&
!                             source_eta(:,1,1),flux_eta(1,:,1,1))
         case('FVM')
            call fvm_solution(kappa,kappa,kappa_sct,Ib,phi,&
               x_eps(:,:,:,ieps),y_eps(:,:,:,ieps),z_eps(:,:,:,ieps),&
               flux_eta,source_eta,total_loss=loss_eta)
!         case('PMC')                                                    !In this approach, the PMC method is executed
!            which_spectral_mc = 'GG'                                    !  by-passing the spectral random number selection;
!            call mc_solution(flux_eta,source_eta,fixed_Ib=Ib,&          !  i.e., the spectral integration is made in the
!                             fixed_kappa=kappa,total_loss=loss_eta)     !  same way as in the FVM
         case default
            call shutdown('lbl_nongray_solution: rte_solution_method &
                          &not specified correctly')
         endselect

         !Update the heat flux and heat source
         call dprint('lbl_nongray_solution: &
                     &Update the heat flux and heat source',3)
         flux_omp = flux_omp + flux_eta*lbl_deta(ilbl)
         source_omp = source_omp + source_eta*lbl_deta(ilbl)
         loss_omp = loss_omp + loss_eta*lbl_deta(ilbl)
 
         !--------------------------------------------------------------
         !Compute and store quantities for the SRI analyses, 
         !if requested
         !--------------------------------------------------------------
         !$OMP CRITICAL
         if (present(sri_output)) then     
            call dprint('lbl_nongray_solution: &
                        &Compute quantities for the SRI analyses',3)
            intV = 0._dp; intV(1) = lbl_xeta(ilbl)
            do i=icm,icp
               do j=jcm,jcp
                  do k=kcm,kcp
                     dx = xf(i) - xf(i-1)                               !Cell size: x direction
                     dy = 1._dp; if (.not.one_d) dy = yf(j) - yf(j-1)   !Cell size: y direction
                     dz = 1._dp; if (three_d)    dz = zf(k) - zf(k-1)   !Cell size: z direction
                     dV = dx*dy*dz                                      !Cell volume
                     intV(2) = intV(2) + dabs(source_eta(i,j,k))*dV
                     intV(3) = intV(3) + source_eta(i,j,k)*dV
                     intV(4) = intV(4) + dabs(flux_eta(1,i,j,k))*dV
                     intV(5) = intV(5) + flux_eta(1,i,j,k)*dV
                     intV(6) = intV(6) + dabs(flux_eta(2,i,j,k))*dV
                     intV(7) = intV(7) + flux_eta(2,i,j,k)*dV
                     intV(8) = intV(8) + dabs(flux_eta(3,i,j,k))*dV
                     intV(9) = intV(9) + flux_eta(3,i,j,k)*dV
                     magflux = dsqrt(flux_eta(1,i,j,k)**2 + &
                                     flux_eta(2,i,j,k)**2 + &
                                     flux_eta(3,i,j,k)**2)
                     intV(10) = intV(10) + magflux*dV
                     intV(11) = intV(11) + Ib(i,j,k)*dV
                     intV(12) = intV(12) + &
                        fourpi*kappa(i,j,k)*Ib(i,j,k)*dV
                     intV(13) = intV(13) + (source(i,j,k) + &
                        fourpi*kappa(i,j,k)*Ib(i,j,k))*dV
                  enddo
               enddo
            enddo
            sri_output(ilbl,:) = intV
         endif
         !$OMP ENDCRITICAL

      enddo lines_loop
      !$OMP ENDDO

      !-----------------------------------------------------------------
      !Finish the calculation of the total radiative quantities
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
      call dprint('lbl_nongray_solution: &
                  &Finish calculation of the optional properties')
      if (trim(rte_solution_method).eq.'exact-los') kIb = kIb/fourpi
      do i=1,nx
         do j=1,ny
            do k=1,nz
               Ib_total = sigma*invpi*T(i,j,k)**4._dp                   !Total blackbody intensity
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
      call dprint('lbl_nongray_solution: Deallocate arrays')
      deallocate(flux_eta,source_eta)
      deallocate(flux_omp,source_omp)
      deallocate(Ib,kIb,kIb_omp)
      deallocate(kappa,kappa_sct,phi)
      deallocate(x_eps,y_eps,z_eps)
      
      !-----------------------------------------------------------------
      !Total times
      !-----------------------------------------------------------------
      call dprint('lbl_nongray_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()
      if (present(total_time)) &                                        !Computing (total)
                        total_time = end_total_time - start_total_time
      
   endsubroutine lbl_nongray_solution

!   !====================================================================
!   !Subroutine for the exact RTE solution via the LBL method for a 
!   !non-scattering, one-dimensional cylindrical medium
!   !====================================================================
!   subroutine lbl_cyl_solution(x,T,p,xs,I0,source,flux,processing_time,&
!                     total_time,absorption,emission,kappaPlanck,points)
   
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
!                                get_wall_time,print_to_prompt,shutdown
!      use constants, only: fourpi,invpi,sigma
!      use exact_solution, only: exact_1d_cyl
!      use global_parameters, only: debug_mode,number_of_species
!      implicit none
!      character(100) :: msg
!      integer :: i,ierr,isp
!      integer :: n,nx,nlines,nsp
!      integer,intent(in),optional :: points
!      real(dp),intent(in) :: I0,p(:),T(:),x(:),xs(:,:)
!      real(dp),intent(out) :: flux(:),source(:)
!      real(dp),intent(out),optional :: absorption(:),emission(:),&
!         kappaPlanck(:),processing_time,total_time
!      real(dp) :: end_proc_time,start_proc_time,&
!                  end_total_time,start_total_time,time_sum
!      real(dp) :: deta,eta_down,eta_up,xeta
!      real(dp) :: Ib_total
!      real(dp),allocatable,dimension(:) :: flux_eta,source_eta
!      real(dp),allocatable,dimension(:) :: Ib,kappa,kIb
!      real(dp),allocatable,dimension(:,:,:) :: acs_down,acs_up
!      logical :: compute_absorption=.false.,&
!         compute_emission=.false.,compute_planck=.false.
!start_total_time = I0 !Remove this
!      !-----------------------------------------------------------------
!      !Preparatory procedures
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                           &preparatory procedures')
      
!      !Starting the total computing time counter
!      start_total_time = get_wall_time()

!      !Check if all arrays have the same size
!      call assert(size(x).eq.size(T))
!      call assert(size(T).eq.size(p))
!      call assert(size(p).eq.size(xs,2))

!      !Zeroing out processing time counter (necessary for the posterior
!      !sum over all instances of the RTE solution)
!      time_sum = 0._dp
      
!      !Check for repeated flags that should be unique
!      call check_lbl_parameters
      
!      !Set flags for optional output arguments
!      compute_absorption=.false.
!      if (present(absorption))  compute_absorption = .true.
!      compute_emission=.false.
!      if (present(emission))    compute_emission = .true.
!      compute_planck=.false.
!      if (present(kappaPlanck)) compute_planck = .true.

!      !-----------------------------------------------------------------
!      !Surrogate names
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                           &surrogate names')
!      nx = size(x); if (present(points)) nx = points                    !Number of grid points
!      nsp = number_of_species                                           !Number of participating species

!      !-----------------------------------------------------------------
!      !Allocating arrays
!      !-----------------------------------------------------------------
!      allocate(acs_down(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_down',ierr)
!      allocate(acs_up(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_up',ierr)
!      allocate(flux_eta(nx),stat=ierr)
!      call CheckMemAlloc('flux_eta',ierr)
!      allocate(Ib(nx),stat=ierr)
!      call CheckMemAlloc('Ib',ierr)
!      allocate(kappa(nx),stat=ierr)
!      call CheckMemAlloc('kappa',ierr)
!      allocate(kIb(nx),stat=ierr)
!      call CheckMemAlloc('kIb',ierr)
!      allocate(source_eta(nx),stat=ierr)
!      call CheckMemAlloc('source_eta',ierr)

!      !-----------------------------------------------------------------
!      !Define, open units and read the first line of data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                           &prepare LBL data files')
!      if (.not.lbl_data_ready.or.(lbl_nlines.le.0)) &                   !Either prepare the LBL data to be read
!         call prepare_lbl_data                                          !  (i.e., get units and open the external files)
!      if (lbl_data_ready.and.(lbl_nlines.gt.0)) call rewind_lbl_data    !  or rewind the already prepared LBL data
!      if (lbl_data_averaging.ne.'none') &
!         call real_lbl_data(eta_down,acs_down)                          !Read the first line of the data
!      nlines = lbl_nlines-1                                             !Total number of lines in the data files
!      if (lbl_data_averaging.eq.'none') nlines = nlines + 1

!      !-----------------------------------------------------------------
!      !Initialize sum arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                           &initialize sum arrays')
!      flux = 0._dp
!      source = 0._dp
!      kIb = 0._dp

!      solu_loop: do n=1,nlines
!         !--------------------------------------------------------------
!         !Prepare the new spectral line
!         !--------------------------------------------------------------
!         call real_lbl_data(eta_up,acs_up,deta)                         !Read a new line of the LBL data files
                 
!         !Adjusting spectral bands for the analysis
!         if (lbl_data_averaging.ne.'none') &                            !Interval between bands
!            deta = (eta_up - eta_down)                                 
!         if ((lbl_data_averaging.eq.'arithmetic').or.&                  !Wavenumber position
!             (lbl_data_averaging.eq.'geometric')) &
!               xeta = 0.5_dp*(eta_up + eta_down)              
!         if ((lbl_data_averaging.eq.'upwind').or.&                      !Wavenumber position (no mean: use the
!             (lbl_data_averaging.eq.'none')) xeta = eta_up              !  current wavenumber value)

!         !Converting from 1/cm to 1/m
!         if (lbl_data_cm) deta = deta*100._dp
!         if (lbl_data_cm) xeta = xeta*100._dp
         
!         !Print xeta if requested
!         if (print_xeta.or.debug_mode) then
!            write(msg,'(f12.4)') xeta/100._dp
!            msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
!            call print_to_prompt(msg)
!         endif
         
!         !Check if this eta is to be solved for
!         if (xeta.lt.lbl_eta_min) then
!            eta_down = eta_up
!            acs_down = acs_up
!            cycle solu_loop
!         endif
!         if (xeta.gt.lbl_eta_max) exit solu_loop

!         !--------------------------------------------------------------
!         !Determine the properties for this wavenumber
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                       &determine the properties for this wavenumber',3)
                                            
!         !Current absorption cross-section array
!         if (lbl_data_averaging.eq.'arithmetic') &
!            lbl_data = 0.5_dp*(acs_up + acs_down)
!         if (lbl_data_averaging.eq.'geometric') &
!            lbl_data = dsqrt(acs_up*acs_down)
!         if ((lbl_data_averaging.eq.'upwind').or.&
!             (lbl_data_averaging.eq.'none')) lbl_data = acs_up
         
!         do i=1,nx
!            !Absorption coefficient
!            kappa(i) = 0._dp
!            do isp=1,nsp
!               kappa(i) = kappa(i) + &
!                  lbl_kappa_func(xs(isp,i),T(i),p(i),isp)
!            enddo
                     
!            !Spectral blackbory radiative intensity
!            Ib(i) = Ib_function(T(i),xeta)
                  
!            !Integrate RTE emission term
!            kIb(i) = kIb(i) + kappa(i)*Ib(i)*deta

!         enddo
         
!         !--------------------------------------------------------------
!         !Solve the radiation field
!         !--------------------------------------------------------------
!         !Start the processing time counter
!         start_proc_time = get_wall_time()
         
!         !Solve the radiation field for the current wavenumber
!         if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                          &solve the radiation field',3)
!         call exact_1d_cyl(x,kappa,Ib,0._dp,nx,source_eta)
         
!         !Update the heat flux and heat source
!         if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                               &update the heat flux and heat source',3)
!         flux = flux + flux_eta*deta
!         source = source + source_eta*deta
         
!         !Stop the processing time counter and adding to the total
!         end_proc_time = get_wall_time()
!         time_sum = time_sum + end_proc_time - start_proc_time
 
!         !--------------------------------------------------------------
!         !Update spectral values
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                             &update spectral values',3)
!         eta_down = eta_up                                              !Wavenumber
!         acs_down = acs_up                                              !Absorption cross-section

!      enddo solu_loop
      
!      !-----------------------------------------------------------------
!      !Close units
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                            &close units')
!      call close_lbl_data

!      !-----------------------------------------------------------------
!      !Final loop over all grid cells to finish up
!      !the calculation of the the optional properties
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                        &finish calculation of the optional properties')
!      do i=1,nx
!         Ib_total = sigma*invpi*T(i)**4._dp                             !Total blackbody intensity
!         if (compute_emission) emission(i) = fourpi*kIb(i)              !Total emission
!         if (compute_absorption) absorption(i) = source(i) +&           !Total absorption
!                                                 fourpi*kIb(i)
!         if (compute_planck) kappaPlanck(i) = kIb(i)/Ib_total           !Planck-mean absorption coefficient
!      enddo
         
!      !-----------------------------------------------------------------
!      !Deallocate arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                           &deallocate arrays')
!      deallocate(acs_down,acs_up)
!      deallocate(flux_eta,source_eta)
!      deallocate(Ib,kappa,kIb)
      
!      !-----------------------------------------------------------------
!      !Total times
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('lbl_cyl_solution: &
!                                           &total times')

!      !Processing
!      if (present(processing_time)) &
!         processing_time = end_proc_time - start_proc_time

!      !Computing (total)
!      call cpu_time(end_total_time)
!      if (present(total_time)) &
!         total_time = end_total_time - start_total_time
      
!   endsubroutine lbl_cyl_solution

endmodule lbl_routines
