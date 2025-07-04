module wbw_routines

   implicit none

contains

   !====================================================================
   !Subroutine for the radiative heat transfer solution via the WBW 
   !model for a non-scattering medium bounded by black walls
   !====================================================================
   subroutine wbw_black_solution(model_id,flux,source,&
      processing_time,total_time,absorption,emission,kappaPlanck,&
      solve_RTE,total_loss,band_id)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,shutdown
      use constants, only: fourpi,sigrpi
      use exact_solution, only: exact_1d_pp
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: rte_solution_method
      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use physical_functions, only: bb_emission_frac
      use precision_parameters, only: dp
      use wbw_functions
      character(*),intent(in) :: model_id
      integer,optional,intent(in) :: band_id
      integer :: i,iband,j,jgas,k,lband,uband
      integer :: nband,ngas,nl,nk,nx,ny,nz
      integer :: ierr
      integer,allocatable,dimension(:,:) :: mj_array
      logical,intent(in),optional :: solve_RTE
      logical :: compute_absorption,compute_emission,compute_planck,&
                 compute_radiation
      real(dp) :: wbw_lbound,wbw_ubound
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints),&
                                 total_loss
      real(dp) :: a,F,kp,loss,loss_j,loss_omp,pa
      real(dp) :: end_proc_time,end_total_time,start_proc_time,&
                  start_total_time
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: kappa_sct,kIb,Ib,&
         source_j,source_omp,x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:) :: flux_j,flux_omp
      real(dp),allocatable,dimension(:,:,:,:,:) :: aIb,kappa
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('wbw_black_solution: Preparatory procedures')
      
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

      !Define number of bands
      if (present(band_id)) then
         lband = band_id; uband = band_id; nband = 1
      else
         nband = wbw_number_bands(model_id)
         lband = 1; uband = nband
      endif

      !Surrogate names
      ngas = wbw_number_gray_gases(model_id)                            !Number of gray gases of the model
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z
      
      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('wbw_black_solution: Allocating arrays')
      allocate(aIb(nx,ny,nz,ngas+1,lband:uband),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,ngas+1,lband:uband),stat=ierr)
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
      call dprint('wbw_black_solution: Compute properties')
      kIb = 0._dp                                                       !Zero out sum array
      do iband=lband,uband
         call wbw_bounds(wbw_lbound,wbw_ubound,model_id,iband)          !Obtain the band bounds
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  !Local partial pressure
                  pa = p(i,j,k)*sum(xs(:,i,j,k))
               
                  !Blackbody intensity
                  Ib(i,j,k) = sigrpi*(T(i,j,k)**4._dp)
                  
                  !Blackbody weighting fraction
                  F = bb_emission_frac(wbw_lbound,T(i,j,k),.true.)-&
                      bb_emission_frac(wbw_ubound,T(i,j,k),.true.)
if (isnan(F)) then
write(*,*) wbw_lbound,wbw_ubound,T(i,j,k)
stop
endif
                  !Absorption coefficient and weighting coefficient
                  do jgas=1,ngas+1
                     kp = wbw_kappaj(model_id,iband,jgas)               !Local pressure-based absorption coefficient
                     kappa(i,j,k,jgas,iband) = kp*pa                    !Local absorption coefficient         
                     a = wbw_aj(model_id,iband,jgas,T(i,j,k))           !Local weighting coefficient
                     aIb(i,j,k,jgas,iband) = a*F*Ib(i,j,k)              !Emission term
                     kIb(i,j,k) = kIb(i,j,k) + &                        !Total RTE emission term
                        kappa(i,j,k,jgas,iband)*aIb(i,j,k,jgas,iband)

                  enddo
               enddo
            enddo
         enddo
      enddo
      
      !Misc parameters needed for the RTE solution subroutine
      kappa_sct = 0._dp                                                 !Scattering coefficient (null)
      phi = 1._dp                                                       !Phase function (isotropic)
      x_eps = 1._dp; y_eps = 1._dp; z_eps = 1._dp                       !Set wall emissivities to unit for all wall cells
      
      !-----------------------------------------------------------------
      !Mount mj_array
      !-----------------------------------------------------------------
      call dprint('wbw_black_solution: Mount mj_array')
      
      !Prepare the array      
      nk = nband*(ngas+1)                                               !Array size
      allocate(mj_array(nk,2),stat=ierr)                                !Allocate array
      call CheckMemAlloc('mj_array',ierr)
      
      !Fill in the array
      k = 0
      do iband=lband,uband
         do jgas=1,ngas+1
            k = k + 1
            mj_array(k,1) = jgas; mj_array(k,2) = iband
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Main radiative transfer solution loop
      !-----------------------------------------------------------------
      call dprint('wbw_black_solution: Main solution loop started')
      flux = 0._dp; source = 0._dp; loss = 0._dp                        !Zero out sum arrays
      start_proc_time = get_wall_time()                                 !Start the processing time counter
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_j,flux_omp,iband,jgas,k,loss_j,loss_omp,&
      !$OMP&           source_j,source_omp)
      
      flux_omp = 0._dp; source_omp = 0._dp; loss_omp = 0._dp            !Zero out OMP sum arrays
      !$OMP DO
      do k=1,nk
         !Define gas and band indexes         
         jgas = mj_array(k,1); iband = mj_array(k,2)
         
         !Solve the radiation field for gas jgas, band iband
         if (compute_radiation) then
            selectcase(rte_solution_method)
            case('exact-1d')
               call exact_1d_pp(x,kappa(:,1,1,jgas,iband),&
                                aIb(:,1,1,jgas,iband),nx,&
                                source_j(:,1,1),flux_j(1,:,1,1))
            case('FVM')                                                 !Finite volume method solution
               call fvm_solution(kappa(:,:,:,jgas,iband),&
                  kappa(:,:,:,jgas,iband),kappa_sct,&
                  aIb(:,:,:,jgas,iband),phi,x_eps,y_eps,z_eps,&
                  flux_j,source_j,dont_iterate=.true.,total_loss=loss_j)

            case('PMC')                                                 !In this approach, the PMC method is executed
               which_spectral_mc = 'GG'                                 !  by-passing the spectral random number selection;
               call mc_solution(flux_j,source_j,fixed_Ib=&              !  i.e., the spectral integration is made in the
                                aIb(:,:,:,jgas,iband),&                 !  same way as in the FVM
                                fixed_kappa=kappa(:,:,:,jgas,iband))
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
      call dprint('wbw_black_solution: Main solution loop ended')

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('wbw_black_solution: &
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
      deallocate(flux_j,flux_omp,source_j,source_omp)
      
      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('wsgg_black_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()       
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time
      
   endsubroutine wbw_black_solution  

   !====================================================================
   !Subroutine for the radiative heat transfer solution via the WBW 
   !model for a non-scattering medium bounded by black walls
   !====================================================================
   subroutine wbw_nongray_solution(model_id,flux,source,&
      processing_time,total_time,absorption,emission,kappaPlanck,&
      solve_RTE,total_loss,band_id)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                get_wall_time,shutdown
      use constants, only: fourpi,sigrpi
      use exact_solution, only: exact_1d_pp
      use fvm_parameters, only: fvm_angles
      use fvm_routines, only: fvm_solution
      use global_parameters, only: rte_solution_method
!      use mc_parameters, only: which_spectral_mc
      use mc_routines
      use mesh
      use physical_functions, only: bb_emission_frac
      use precision_parameters, only: dp
      use wbw_functions
      character(*),intent(in) :: model_id
      integer,optional,intent(in) :: band_id
      integer :: i,iband,j,jgas,k,lband,uband
      integer :: nband,ngas,nl,nk,nx,ny,nz
      integer :: ierr
      integer,allocatable,dimension(:,:) :: mj_array
      logical,intent(in),optional :: solve_RTE
      logical :: compute_absorption,compute_emission,compute_planck,&
                 compute_radiation
      real(dp) :: wbw_lbound,wbw_ubound
      real(dp),intent(out) :: flux(3,xpoints,ypoints,zpoints),&
                              source(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: processing_time,total_time,&
                                 absorption(xpoints,ypoints,zpoints),&
                                 emission(xpoints,ypoints,zpoints),&
                                 kappaPlanck(xpoints,ypoints,zpoints),&
                                 total_loss
      real(dp) :: a,F,kp,loss,loss_j,loss_omp,pa
      real(dp) :: end_proc_time,end_total_time,start_proc_time,&
                  start_total_time
      real(dp),allocatable,dimension(:,:) :: phi
      real(dp),allocatable,dimension(:,:,:) :: kappa_sct,kIb,Ib,&
                                               source_j,source_omp
      real(dp),allocatable,dimension(:,:,:,:) :: flux_j,flux_omp,&
                                                 x_eps,y_eps,z_eps
      real(dp),allocatable,dimension(:,:,:,:,:) :: aIb,kappa
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('wbw_nongray_solution: Preparatory procedures')
      
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

      !Define number of bands
      if (present(band_id)) then
         lband = band_id; uband = band_id; nband = 1
      else
         nband = wbw_number_bands(model_id)
         lband = 1; uband = nband
      endif

      !Surrogate names
      ngas = wbw_number_gray_gases(model_id)                            !Number of gray gases of the model
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nx = xpoints                                                      !Number of cells along direction x
      ny = ypoints                                                      !Number of cells along direction y
      nz = zpoints                                                      !Number of cells along direction z

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('wbw_nongray_solution: Allocating arrays')
      allocate(aIb(nx,ny,nz,ngas+1,lband:uband),stat=ierr)
      call CheckMemAlloc('aIb',ierr)
      allocate(flux_j(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_j',ierr)
      allocate(flux_omp(3,nx,ny,nz),stat=ierr)
      call CheckMemAlloc('flux_omp',ierr)
      allocate(Ib(nx,ny,nz),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(kappa(nx,ny,nz,ngas+1,lband:uband),stat=ierr)
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
      call dprint('wbw_nongray_solution: Compute properties')
      kIb = 0._dp                                                       !Zero out sum array
      do iband=lband,uband
         call wbw_bounds(wbw_lbound,wbw_ubound,model_id,iband)          !Obtain the band bounds
         do i=1,nx
            do j=1,ny
               do k=1,nz
                  !Local partial pressure
                  pa = p(i,j,k)*sum(xs(:,i,j,k))
               
                  !Blackbody intensity
                  Ib(i,j,k) = sigrpi*(T(i,j,k)**4._dp)
                  
                  !Blackbody weighting fraction
                  F = bb_emission_frac(wbw_lbound,T(i,j,k),.true.)-&
                      bb_emission_frac(wbw_ubound,T(i,j,k),.true.)
                  
                  !Absorption coefficient and weighting coefficient
                  do jgas=1,ngas+1
                     kp = wbw_kappaj(model_id,iband,jgas)               !Local pressure-based absorption coefficient
                     kappa(i,j,k,jgas,iband) = kp*pa                    !Local absorption coefficient         
                     a = wbw_aj(model_id,iband,jgas,T(i,j,k))           !Local weighting coefficient
                     aIb(i,j,k,jgas,iband) = a*F*Ib(i,j,k)              !Emission term
                     kIb(i,j,k) = kIb(i,j,k) + &                        !Total RTE emission term
                        kappa(i,j,k,jgas,iband)*aIb(i,j,k,jgas,iband)
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
      call dprint('wbw_nongray_solution: Define wall emissivities')
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
      !Mount mj_array
      !-----------------------------------------------------------------
      call dprint('wbw_nongray_solution: Mount mj_array')
      
      !Prepare the array      
      nk = nband*(ngas+1)                                               !Array size
      allocate(mj_array(nk,2),stat=ierr)                                !Allocate array
      call CheckMemAlloc('mj_array',ierr)
      
      !Fill in the array
      k = 0
      do iband=lband,uband
         do jgas=1,ngas+1
            k = k + 1
            mj_array(k,1) = jgas; mj_array(k,2) = iband
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Main radiative transfer solution loop
      !-----------------------------------------------------------------
      call dprint('wbw_nongray_solution: Main solution loop started')
      flux = 0._dp; source = 0._dp; loss = 0._dp                        !Zero out sum arrays
      start_proc_time = get_wall_time()                                 !Start the processing time counter
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(flux_j,flux_omp,iband,jgas,k,loss_j,loss_omp,&
      !$OMP&           source_j,source_omp)
      
      flux_omp = 0._dp; source_omp = 0._dp; loss_omp = 0._dp            !Zero out OMP sum arrays
      !$OMP DO
      do k=1,nk
         !Define gas and band indexes         
         jgas = mj_array(k,1); iband = mj_array(k,2)
         
         !Solve the radiation field for gas jgas, band iband
         if (compute_radiation) then
            selectcase(rte_solution_method)
!            case('exact-1d')
!               call exact_1d_pp(x,kappa(:,1,1,jgas,iband),&
!                                aIb(:,1,1,jgas,iband),nx,&
!                                source_j(:,1,1),flux_j(1,:,1,1))
            case('FVM')                                                 !Finite volume method solution
               call fvm_solution(kappa(:,:,:,jgas,iband),&
                  kappa(:,:,:,jgas,iband),kappa_sct,&
                  aIb(:,:,:,jgas,iband),phi,x_eps(:,:,:,iband),&
                  y_eps(:,:,:,iband),z_eps(:,:,:,iband),flux_j,&
                  source_j,total_loss=loss_j)

!            case('PMC')                                                 !In this approach, the PMC method is executed
!               which_spectral_mc = 'GG'                                 !  by-passing the spectral random number selection;
!               call mc_solution(flux_j,source_j,fixed_Ib=&              !  i.e., the spectral integration is made in the
!                                aIb(:,:,:,jgas,iband),&                 !  same way as in the FVM
!                                fixed_kappa=kappa(:,:,:,jgas,iband))
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
      call dprint('wbw_nongray_solution: Main solution loop ended')

      !-----------------------------------------------------------------
      !Final loop over all grid cells to finish up
      !the calculation of the the optional properties
      !-----------------------------------------------------------------
      call dprint('wbw_nongray_solution: &
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
      call dprint('wbw_nongray_solution: Deallocating arrays')
      deallocate(aIb,Ib,kappa,kIb)
      deallocate(kappa_sct,phi,x_eps,y_eps,z_eps)
      deallocate(flux_j,flux_omp,source_j,source_omp)
      
      !-----------------------------------------------------------------
      !Computational times
      !-----------------------------------------------------------------
      call dprint('wbw_nongray_solution: Compute total times')
      if (present(processing_time)) &                                   !Processing
         processing_time = end_proc_time - start_proc_time
      end_total_time = get_wall_time()       
      if (present(total_time)) &                                        !Computing (total)
         total_time = end_total_time - start_total_time
      
   endsubroutine wbw_nongray_solution

endmodule wbw_routines
