module wbm_routines

   implicit none
   
   contains
   
   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !box model for a non-scattering medium bounded by black walls
   !====================================================================
   subroutine box_black_solution(model_id,flux,source,processing_time,&
      total_time,absorption,emission,kappaPlanck,solve_RTE,total_loss,&
      band_id)
                                 
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: fourpi,sigrpi
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
         get_wall_time,shutdown
      use exact_solution, only: exact_1d_pp
      use fvm_parameters
      use fvm_routines, only: fvm_solution
      use global_parameters, only: rte_solution_method
      use physical_functions, only: bb_emission_frac
      use precision_parameters, only: dp
      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use wbm_parameters
      use wbm_functions
      character(*),intent(in) :: model_id
      integer,optional,intent(in) :: band_id
      integer :: i,j,k,nl,nx,ny,nz,ierr,uout
      integer :: iband,lband,nband,uband
      logical,intent(in),optional :: solve_RTE
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints),&
                                 total_loss
      real(dp) :: F,loss,loss_j,loss_omp
      real(dp) :: end_proc_time,start_proc_time,end_total_time,&
         start_total_time
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: Ib,kappa_sct,kIb,&
         source_j,source_omp,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: fIb,flux_j,flux_omp,&
         kappa
      logical :: compute_absorption,compute_emission,compute_planck,&
         compute_radiation

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('box_black_solution: Preparatory procedures')
      
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

      !Prepare the WBM data
      call wbm_prepare_model(model_id)

      !Define number of bands
      if (present(band_id)) then
         lband = band_id; uband = band_id; nband = 1
      else
         nband = number_wbm_bands
         lband = 1; uband = nband
      endif

      !Surrogate names
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      
      !-----------------------------------------------------------------
      !Prepare output unit if requested
      !-----------------------------------------------------------------
      if (wbm_print_properties) then
         if (trim(wbm_print_file).eq.'null') &                          !Check if a file name was given                          
            call shutdown('wbm_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(wbm_print_file))                      !Open file
         write(uout,'(8(a,:,","))') 'eta-','eta+','x','y','z',&         !Write header
            'kappa','f','fIb'
         write(uout,'(8(a,:,","))') '[1/m]','[1/m]','[m]','[m]','[m]',& !Write subheader
            '[1/m]','[-]','[W/m2]'
      endif

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('box_black_solution: Allocating arrays')
      allocate(fIb(nx,ny,nz,lband:uband),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,lband:uband),stat=ierr)
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
      call CheckMemAlloc('source_omp',ierr)
      allocate(x_eps(2,ny,nz),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)
      
      !-----------------------------------------------------------------
      !Computing properties
      !-----------------------------------------------------------------
      call dprint('box_black_solution: Compute properties')
      kIb = 0._dp                                                       !Zero out sum array
      do iband=lband,uband
         do i=1,nx
            do j=1,ny
               do k=1,nz               
                  !Blackbody intensity
                  Ib(i,j,k) = sigrpi*(T(i,j,k)**4._dp)
                  
                  !Blackbody weighting fraction
                  F = bb_emission_frac(wbm_lbound(iband),T(i,j,k),.true.)-&
                      bb_emission_frac(wbm_ubound(iband),T(i,j,k),.true.)
                  
                  !Absorption coefficient and weighting coefficient
                  kappa(i,j,k,iband) = box_kappa_func(model_id,iband,&
                                          T(i,j,k),xs(:,i,j,k),p(i,j,k))
                  fIb(i,j,k,iband) = F*Ib(i,j,k)                        !Emission term
                  kIb(i,j,k) = kIb(i,j,k) + &                           !Total RTE emission term
                     kappa(i,j,k,iband)*fIb(i,j,k,iband)
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
      call dprint('box_black_solution: Main solution loop started')
      flux = 0._dp; source = 0._dp; loss = 0._dp                        !Zero out sum arrays
      start_proc_time = get_wall_time()                                 !Start the processing time counter
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_j,flux_omp,loss_j,loss_omp,source_j,&
      !$OMP&           source_omp)
      
      flux_omp = 0._dp; source_omp = 0._dp; loss_omp = 0._dp            !Zero out OMP sum arrays
      !$OMP DO
      band_loop: do iband=lband,uband
         
         !Solve the radiation field for gas jgas, band iband
         if (compute_radiation) then
            selectcase(rte_solution_method)
            case('exact-1d')
               call exact_1d_pp(x,kappa(:,1,1,iband),fIb(:,1,1,iband),&
                                nx,source_j(:,1,1),flux_j(1,:,1,1))
                                
            case('FVM')                                                 !Finite volume method solution
               call fvm_solution(kappa(:,:,:,iband),kappa(:,:,:,iband),&
                  kappa_sct,fIb(:,:,:,iband),phi,x_eps,y_eps,z_eps,&
                  flux_j,source_j,dont_iterate=.true.,total_loss=loss_j)

            case('PMC')                                                 !In this approach, the PMC method is executed
               which_spectral_mc = 'GG'                                 !  by-passing the spectral random number selection;
               call mc_solution(flux_j,source_j,fixed_Ib=&              !  i.e., the spectral integration is made in the
                                fIb(:,:,:,iband),&                      !  same way as in the FVM
                                fixed_kappa=kappa(:,:,:,iband))
!           case('P1')
!              call p1_solution(kappa,aIb,flux_j,source_j)              !P1 method
            endselect
         endif
         
         !Update the heat flux, heat source and total loss
         flux_omp = flux_omp + flux_j
         source_omp = source_omp + source_j
         loss_omp = loss_omp + loss_j
      
      enddo band_loop
      !$OMP ENDDO
      
      !Finish the calculation of the total radiative quantities
      !$OMP CRITICAL
      flux = flux + flux_omp
      source = source + source_omp
      loss = loss + loss_omp
      !$OMP ENDCRITICAL
      
      !$OMP ENDPARALLEL
      end_proc_time = get_wall_time()                                   !End the processing time counter
      call dprint('box_black_solution: Main solution loop ended')
      
      !-----------------------------------------------------------------
      !Close output unit if it was used
      !-----------------------------------------------------------------
      if (wbm_print_properties) close(uout)
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('box_black_solution: &
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
      call dprint('box_black_solution: Deallocating arrays')
      deallocate(fIb,Ib,kappa,kIb)
      deallocate(kappa_sct,phi,x_eps,y_eps,z_eps)
      deallocate(flux_j,flux_omp,source_j,source_omp)
      
      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('box_black_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()       
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time

   endsubroutine box_black_solution
   
   !====================================================================
   !Subroutine for the radiative heat transfer solution via the
   !box model for a non-scattering medium bounded by non-gray walls
   !====================================================================
   subroutine box_nongray_solution(model_id,flux,source,&
      processing_time,total_time,absorption,emission,kappaPlanck,&
      solve_RTE,total_loss,band_id)
                                 
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: fourpi,sigrpi
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
         get_wall_time,shutdown
      use exact_solution, only: exact_1d_pp
      use fvm_parameters
      use fvm_routines, only: fvm_solution
      use global_parameters, only: rte_solution_method
      use physical_functions, only: bb_emission_frac
      use precision_parameters, only: dp
!      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use wbm_parameters
      use wbm_functions
      character(*),intent(in) :: model_id
      integer,optional,intent(in) :: band_id
      integer :: i,j,k,nl,nx,ny,nz,ierr,uout
      integer :: iband,lband,nband,uband
      logical,intent(in),optional :: solve_RTE
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints),&
                                 total_loss
      real(dp) :: F,loss,loss_j,loss_omp
      real(dp) :: end_proc_time,start_proc_time,end_total_time,&
         start_total_time
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: Ib,kappa_sct,kIb,&
                                               source_j,source_omp
      real(dp),allocatable,dimension(:,:,:,:) :: fIb,flux_j,flux_omp,&
                                                 kappa,x_eps,y_eps,z_eps
      logical :: compute_absorption,compute_emission,compute_planck,&
         compute_radiation

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('box_nongray_solution: Preparatory procedures')
      
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

      !Prepare the WBM data
      call wbm_prepare_model(model_id)

      !Define number of bands
      if (present(band_id)) then
         lband = band_id; uband = band_id; nband = 1
      else
         nband = number_wbm_bands
         lband = 1; uband = nband
      endif

      !Surrogate names
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      
      !-----------------------------------------------------------------
      !Prepare output unit if requested
      !-----------------------------------------------------------------
      if (wbm_print_properties) then
         if (trim(wbm_print_file).eq.'null') &                          !Check if a file name was given                          
            call shutdown('wbm_print_file undefined')
         uout = get_file_unit()                                         !Get unit
         open(unit=uout,file=trim(wbm_print_file))                      !Open file
         write(uout,'(8(a,:,","))') 'eta-','eta+','x','y','z',&         !Write header
            'kappa','f','fIb'
         write(uout,'(8(a,:,","))') '[1/m]','[1/m]','[m]','[m]','[m]',& !Write subheader
            '[1/m]','[-]','[W/m2]'
      endif

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('box_nongray_solution: Allocating arrays')
      allocate(fIb(nx,ny,nz,lband:uband),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,lband:uband),stat=ierr)
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
      call CheckMemAlloc('source_omp',ierr)
      allocate(x_eps(2,ny,nz,lband:uband),stat=ierr)
      call CheckMemAlloc('x_eps',ierr)
      allocate(y_eps(2,nx,nz,lband:uband),stat=ierr)
      call CheckMemAlloc('y_eps',ierr)
      allocate(z_eps(2,nx,ny,lband:uband),stat=ierr)
      call CheckMemAlloc('z_eps',ierr)
      
      !-----------------------------------------------------------------
      !Computing properties
      !-----------------------------------------------------------------
      call dprint('box_nongray_solution: Compute properties')
      kIb = 0._dp                                                       !Zero out sum array
      do iband=lband,uband
         do i=1,nx
            do j=1,ny
               do k=1,nz               
                  !Blackbody intensity
                  Ib(i,j,k) = sigrpi*(T(i,j,k)**4._dp)
                  
                  !Blackbody weighting fraction
                  F = bb_emission_frac(wbm_lbound(iband),T(i,j,k),.true.)-&
                      bb_emission_frac(wbm_ubound(iband),T(i,j,k),.true.)
                  
                  !Absorption coefficient and weighting coefficient
                  kappa(i,j,k,iband) = box_kappa_func(model_id,iband,&
                                          T(i,j,k),xs(:,i,j,k),p(i,j,k))
                  fIb(i,j,k,iband) = F*Ib(i,j,k)                        !Emission term
                  kIb(i,j,k) = kIb(i,j,k) + &                           !Total RTE emission term
                     kappa(i,j,k,iband)*fIb(i,j,k,iband)
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
      call dprint('box_nongray_solution: Define wall emissivities')
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Initialize arrays
      do iband=lband,uband
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  if (i.eq.1)  x_eps(1,j,k,iband) = &
                     xmin_emissivity(j,k,iband)
                  if (i.eq.nx) x_eps(2,j,k,iband) = &
                     xmax_emissivity(j,k,iband)
                  
                  if (j.eq.1)  y_eps(1,i,k,iband) = &
                     ymin_emissivity(i,k,iband)
                  if (j.eq.ny) y_eps(2,i,k,iband) = &
                     ymax_emissivity(i,k,iband)
                  
                  if (k.eq.1)  z_eps(1,i,j,iband) = &
                     zmin_emissivity(i,j,iband)
                  if (k.eq.nz) z_eps(2,i,j,iband) = &
                     zmax_emissivity(i,j,iband)
               enddo
            enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      !Main radiative transfer solution loop
      !-----------------------------------------------------------------
      call dprint('box_nongray_solution: Main solution loop started')
      flux = 0._dp; source = 0._dp; loss = 0._dp                        !Zero out sum arrays
      start_proc_time = get_wall_time()                                 !Start the processing time counter
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_j,flux_omp,loss_j,loss_omp,source_j,&
      !$OMP&           source_omp)
      
      flux_omp = 0._dp; source_omp = 0._dp; loss_omp = 0._dp            !Zero out OMP sum arrays
      !$OMP DO
      band_loop: do iband=lband,uband
         
         !Solve the radiation field for gas jgas, band iband
         if (compute_radiation) then
            selectcase(rte_solution_method)
            case('exact-1d')
               call exact_1d_pp(x,kappa(:,1,1,iband),fIb(:,1,1,iband),&
                                nx,source_j(:,1,1),flux_j(1,:,1,1))
                                
            case('FVM')                                                 !Finite volume method solution
               call fvm_solution(kappa(:,:,:,iband),kappa(:,:,:,iband),&
                  kappa_sct,fIb(:,:,:,iband),phi,x_eps(:,:,:,iband),&
                  y_eps(:,:,:,iband),z_eps(:,:,:,iband),flux_j,&
                  source_j,total_loss=loss_j)

!            case('PMC')                                                 !In this approach, the PMC method is executed
!               which_spectral_mc = 'GG'                                 !  by-passing the spectral random number selection;
!               call mc_solution(flux_j,source_j,fixed_Ib=&              !  i.e., the spectral integration is made in the
!                                fIb(:,:,:,iband),&                      !  same way as in the FVM
!                                fixed_kappa=kappa(:,:,:,iband))
!           case('P1')
!              call p1_solution(kappa,aIb,flux_j,source_j)              !P1 method
            endselect
         endif
         
         !Update the heat flux, heat source and total loss
         flux_omp = flux_omp + flux_j
         source_omp = source_omp + source_j
         loss_omp = loss_omp + loss_j
      
      enddo band_loop
      !$OMP ENDDO
      
      !Finish the calculation of the total radiative quantities
      !$OMP CRITICAL
      flux = flux + flux_omp
      source = source + source_omp
      loss = loss + loss_omp
      !$OMP ENDCRITICAL
      
      !$OMP ENDPARALLEL
      end_proc_time = get_wall_time()                                   !End the processing time counter
      call dprint('box_nongray_solution: Main solution loop ended')
      
      !-----------------------------------------------------------------
      !Close output unit if it was used
      !-----------------------------------------------------------------
      if (wbm_print_properties) close(uout)
      
      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('box_nongray_solution: &
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
      call dprint('box_nongray_solution: Deallocating arrays')
      deallocate(fIb,Ib,kappa,kIb)
      deallocate(kappa_sct,phi,x_eps,y_eps,z_eps)
      deallocate(flux_j,flux_omp,source_j,source_omp)
      
      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('box_nongray_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()       
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time

   endsubroutine box_nongray_solution

endmodule wbm_routines
