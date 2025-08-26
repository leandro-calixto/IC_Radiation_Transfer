module mc_routines

   !====================================================================
   !Modules & Misc
   !====================================================================
   use mc_parameters
   use mc_functions
   use mesh
   implicit none
   
contains
   
   !====================================================================
   !Subroutine to solve the RTE via the Monte Carlo method
   !(currently, only for non-scattering media bounded by black walls)
   !====================================================================
   subroutine mc_solution(qrad,Srad,fixed_Ib,fixed_kappa,kappa_planck,&
      absorption,emission,total_absorption,total_emission,total_loss)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,print_to_prompt,shutdown
      use constants, only: fourpi,pi,pio2,thrpio2,sigrpi,twopi
      use gg_functions, only: gg_kappa_func
      use global_parameters, only: debug_mode
      use math_functions, only: locate
      use misc_functions, only: set_random_seed,find_imin
      use precision_parameters, only: dp,small,big
      use wsgg_functions, only: wsgg_kappa_planck,wsgg_properties
      use wsgg_parameters, only: which_wsgg_model
      implicit none
      character(20) :: istr,jstr,kstr,nstr,rstr
      character(200) :: msg
      integer :: ierr,idum,irun
      integer :: i,j,k
      integer :: np,nrun,nx,ny,nz
      integer :: icm,icp,iwm,iwp,jcm,jcp,jwm,jwp,kcm,kcp,kwm,kwp
      integer :: idir,jdir,kdir,next_dir
      integer :: ipath,jpath,kpath,path_counter,path_max,path_size
      integer :: jgas
      integer(i8b) :: n,nr
      integer,allocatable,dimension(:,:) :: path_ijk
      integer(i8b),allocatable,dimension(:,:,:) :: n_emi
      integer,allocatable,dimension(:,:,:) :: cell_type
      real(dp) :: psi,theta,R_abs,R_eta,R_psi,R_theta,R_x,R_y,R_z
      real(dp) :: abs0,bundle_emi,ckds,tkIb,ds,dface(3),s1,s2,s3
      real(dp) :: xpos,xfpos,ypos,yfpos,zpos,zfpos
      real(dp) :: dummy_dp,eta,kappa
      real(dp) :: net_run,total_absorption_run,total_emission_run
      real(dp) :: avcoeff,qrad_lcl,Srad_lcl
      real(dp),intent(out) :: qrad(3,xpoints,ypoints,zpoints),&
                              Srad(xpoints,ypoints,zpoints)
      real(dp),intent(in),optional :: &
         fixed_Ib(xpoints,ypoints,zpoints),&
         fixed_kappa(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: total_absorption,total_emission,&
         total_loss
      real(dp),intent(out),optional :: &
         kappa_planck(xpoints,ypoints,zpoints),&
         emission(xpoints,ypoints,zpoints),&
         absorption(xpoints,ypoints,zpoints)
      real(dp),allocatable,dimension(:) :: dx,dy,dz,ckIb,tau
      real(dp),allocatable,dimension(:) :: tkIb_y
      real(dp),allocatable,dimension(:,:) :: ckIb_y,tkIb_z
      real(dp),allocatable,dimension(:,:,:) :: Ib,kplanck
      real(dp),allocatable,dimension(:,:,:) :: bundle_abs,ckIb_z,kIb,dV
      logical :: compute_absorption,compute_emission
      logical :: ggmc,lblmc,wsggmc
      logical :: losmc = .false.                                        !For debugging only

      if (mesh_practice.eq.'A') &
         call shutdown('fvm_solution: Practice A not yet implemented')
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      if (debug_mode) &
         call print_to_prompt('mc_solution: surrogate names',3)
      nx = xpoints                                                      !Number of grid points along direction x
      ny = ypoints                                                      !Number of grid points along direction y
      nz = zpoints                                                      !Number of grid points along direction z
      np = nx*ny*nz                                                     !Total number of grid points
      nr = max(mc_total_rays,mc_rays_per_point*int(np,i8b))             !Total number of rays in the MC method
      nrun = mc_total_runs                                              !Number of consecutive runs

      !-----------------------------------------------------------------
      !Set up flags
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('mc_solution: &
                                           &set up flags',3)
      !Spectral model
      ggmc = .false.
      if (trim(which_spectral_mc).eq.'GG') ggmc = .true.
      lblmc = .false.
      if (trim(which_spectral_mc).eq.'LBL') lblmc = .true.
      wsggmc = .false.
      if (trim(which_spectral_mc).eq.'WSGG') wsggmc = .true.

      !Optional arguments
      compute_absorption = .false.
      if (present(absorption)) compute_absorption = .true.
      compute_emission = .false.
      if (present(emission)) compute_emission = .true.
      
      !Misc
      if (losmc.and.(.not.one_d)) &
         call shutdown('mc_solution: set one_d = .true. for losmc')

      !-----------------------------------------------------------------
      !Set up limiting grid indexes
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('mc_solution: set up &
                                           &limiting grid indexes',3)
      iwm = 1                                                           !x-index of the first wall point
      iwp = nx                                                          !x-index of the last wall point
      icm = iwm + 1                                                     !x-index of the first grid point in the medium
      icp = iwp - 1                                                     !x-index of the last grid point in the medium
      
      jwm = 1; jcm = jwm; jcp = jcm; jwp = jcp                          !Default values for the y-indexes (for 1D calculations)
      if (two_d.or.three_d) then
         jwm = 1                                                        !y-index of the first wall point
         jwp = ny                                                       !y-index of the last wall point
         jcm = jwm + 1                                                  !y-index of the first grid point in the medium
         jcp = jwp - 1                                                  !y-index of the last grid point in the medium
      endif
      
      kwm = 1; kcm = kwm; kcp = kcm; kwp = kcp                          !Default values for the z-indexes (for 1D and 2D calculations)
      if (three_d) then
         kwm = 1                                                        !z-index of the first point cell
         kwp = nz                                                       !z-index of the last point cell
         kcm = kwm + 1                                                  !z-index of the first point cell in the medium
         kcp = kwp - 1                                                  !z-index of the last point cell in the medium
      endif
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('mc_solution: &
                                           &allocate arrays',3)
!      path_max = ceiling(sqrt((iwp - iwm + 1)**2._dp + &                    !Maximum number of cells in a single optical path
!                          (jwp - jwm + 1)**2._dp + &                    !   (assume the path passes through 
!                          (kwp - kwm + 1)**2._dp))                      !    the diagonal of the domain)
      path_max = np                                                     !Maximum number of cells in a single optical path
      allocate(n_emi(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)                !Number of bundles emitted at each grid point
      call CheckMemAlloc('n_emi',ierr)
      allocate(bundle_abs(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)           !Total energy (in a single run) absorbed by each grid point
      call CheckMemAlloc('bundle_abs',ierr)
      allocate(cell_type(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)            !Gives the type of grid point (if located within the medium,
      call CheckMemAlloc('cell_type',ierr)                              !  at a wall, or in a corner)
      allocate(ckIb(iwm:iwp),stat=ierr)                                 !Cumulative emitted energy at a grid point index along the x direction
      call CheckMemAlloc('ckIb',ierr)                                   !  (computed from the left, integrated across the cross-section area)
      allocate(ckIb_y(iwm:iwp,jwm:jwp),stat=ierr)                       !Cumulative emitted energy at a grid point indexes along the y direction
      call CheckMemAlloc('ckIb_y',ierr)                                 !  (for a given x index, computed from the bottom corner, integrated 
                                                                        !  across the z direction)
      allocate(ckIb_z(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)               !Cumulative emitted energy at a grid point, computed with reference to the
      call CheckMemAlloc('ckIb_z',ierr)                                 !  z_min position (at given x and y indexes)
      allocate(dV(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)                   !Cell volume
      call CheckMemAlloc('dV',ierr)
      allocate(dx(iwm:iwp),stat=ierr)                                   !Cell size along the x direction
      call CheckMemAlloc('dx',ierr)
      allocate(dy(jwm:jwp),stat=ierr)                                   !Cell size along the y direction
      call CheckMemAlloc('dy',ierr)
      allocate(dz(kwm:kwp),stat=ierr)                                   !Cell size along the z direction
      call CheckMemAlloc('dz',ierr)
      allocate(Ib(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)                   !Local total blackbody intensity
      call CheckMemAlloc('Ib',ierr)
      allocate(kIb(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)                  !Local emitted energy per unit volume
      call CheckMemAlloc('kIb',ierr)
      allocate(kplanck(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)              !Local Planck-mean absorption coefficient
      call CheckMemAlloc('kplanck',ierr)
      allocate(path_ijk(path_max,1:3),stat=ierr)                        !Array with the indexes of the cells along the path traveled
      call CheckMemAlloc('path_ijk',ierr)                               !  by the emitted bundle, ordered by distance from the emitting point
      allocate(tau(path_max),stat=ierr)                                 !Optical thickness of a cell along the path traveled
      call CheckMemAlloc('tau',ierr)                                    !  by the emitted bundle
      allocate(tkIb_y(iwm:iwp),stat=ierr)                               !Total emitted energy integrated along the yz-plane at a given
      call CheckMemAlloc('tkIb_y',ierr)                                 !  cell x-index [W/m]
      allocate(tkIb_z(iwm:iwp,jwm:jwp),stat=ierr)                       !Total emitted energy integrated along the z direction at a given
      call CheckMemAlloc('tkIb_z',ierr)                                 !  cell x- and y-indexes [W/m2]
      
      !-----------------------------------------------------------------
      !Compute cell dimensions
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('mc_solution: &
                                           &compute cell dimensions',3)
      dx = 1._dp                                                        !Default x grid cell dimension 
      do i=icm,icp                                                      !  (keep equal to 1 for x-normal walls)
         dx(i) = xf(i) - xf(i-1)                                        !x grid cell dimension
      enddo

      dy = 1._dp                                                        !Default y grid cell dimension;
      if (two_d.or.three_d) then                                        !  only correct it for 2d and 3d geometries
         do j=jcm,jcp
            dy(j) = yf(j) - yf(j-1)                                     !y grid cell dimension
         enddo
      endif

      dz = 1._dp                                                        !Default y grid cell dimension;
      if (three_d) then                                                 !  only correct it for 3d geometries
         do k=kcm,kcp
            dz(k) = zf(k) - zf(k-1)                                     !z grid cell dimension
         enddo
      endif

      dface = big                                                       !Default cell face size
      
      !-----------------------------------------------------------------
      !Define if the cell is within the medium (0), in a wall (1), or
      !in a corner (2)
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('mc_solution: &
                                      &define cell indentifier index',3)
      cell_type = 0                                                     !By default, all grid points correspond to
                                                                        !  cells that are within the medium
      cell_type(iwm,:,:) = 1; cell_type(iwp,:,:) = 1                    !East and west walls
      if (two_d.or.three_d) then
         cell_type(:,jwm,:) = 1; cell_type(:,jwp,:) = 1                 !North and south walls
         cell_type(iwm,jwm,:) = 2; cell_type(iwp,jwm,:) = 2             !SW and SE vertices 
         cell_type(iwm,jwp,:) = 2; cell_type(iwp,jwp,:) = 2             !NW and NE vertices
      endif
      if (three_d) then
         cell_type(:,:,kwm) = 1; cell_type(:,:,kwp) = 1                 !Top and bottom walls           
         cell_type(iwm,jwm,:) = 2; cell_type(iwp,jwm,:) = 2             !SW and SE edges 
         cell_type(iwm,jwp,:) = 2; cell_type(iwp,jwp,:) = 2             !NW and NE edges
         cell_type(iwm,:,kwm) = 2; cell_type(iwp,:,kwm) = 2             !BW and BE edges 
         cell_type(iwm,:,kwp) = 2; cell_type(iwp,:,kwp) = 2             !TW and TE edges
         cell_type(:,jwm,kwm) = 2; cell_type(:,jwp,kwm) = 2             !BS and BN edges 
         cell_type(:,jwm,kwp) = 2; cell_type(:,jwp,kwp) = 2             !TS and TN edges
      endif   
      
      !-----------------------------------------------------------------
      !Set up the GG Monte Carlo method
      !-----------------------------------------------------------------
      if (ggmc) then
         if (debug_mode) call print_to_prompt('mc_solution: &
                                            &set up the GG/MC method',3)
         if (present(fixed_kappa)) then                                 !If the fixed_kappa optional argument is present,
            kplanck = fixed_kappa                                       !  use the provided array for the local Planck-
                                                                        !  mean absorption coefficient values
         else
            kplanck = 0._dp                                             !The Planck-mean absorption coefficient of grid points
            do k=kcm,kcp                                                !  not within the medium will be zero
               do j=jcm,jcp
                  do i=icm,icp
                     kplanck(i,j,k) = gg_kappa_func(T(i,j,k),p(i,j,k),& !Compute the local Planck-mean absorption
                                                    xs(:,i,j,k))        !   coefficient for grid points within the medium
                  enddo
               enddo
            enddo
         endif
      endif

      !-----------------------------------------------------------------
      !Set up the WSGG Monte Carlo method
      !-----------------------------------------------------------------
      if (wsggmc) then
         call shutdown('WSGG MC not ready')
         if (debug_mode) call print_to_prompt('mc_solution: &
                                          &set up the WSGG/MC method',3)
         
         !Compute the Planck-mean absorption coefficient
         kplanck = 0._dp                                                !The Planck-mean absorption coefficient of grid points
         do k=kcm,kcp                                                   !  not within the medium will be zero
               do j=jcm,jcp
                  do i=icm,icp              
                  kplanck(i,j,k) = wsgg_kappa_planck(T(i,j,k),p(i,j,k),&!Compute the local Planck-mean absorption
                                                     xs(:,i,j,k))       !   coefficient for grid points within the medium
               enddo
            enddo
         enddo 
      endif

      !-----------------------------------------------------------------
      !Set up the LBL Monte Carlo method
      !-----------------------------------------------------------------
      if (lblmc) then
         if (debug_mode) call print_to_prompt('mc_solution: &
                                           &set up the LBL/MC method',3)
         
         !Prepare the spectral random number data
         call prepare_lblrnr
         
         !Compute the Planck-mean absorption coefficient
         kplanck = 0._dp                                                !The Planck-mean absorption coefficient of grid points
            do k=kcm,kcp                                                !  not within the medium will be zero
               do j=jcm,jcp
                  do i=icm,icp
                  kplanck(i,j,k) = p(i,j,k)*get_mc_kplanck(T(i,j,k),&   !Compute the local Planck-mean absorption
                                                           xs(:,i,j,k)) !   coefficient for grid points within the medium
               enddo
            enddo
         enddo 
      endif
      
      !-----------------------------------------------------------------
      !Precomputing cell volume and direction-integrated emission
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('mc_solution: &
         &precompute cell volume and direction-integrated emission',3)
      kIb = 0._dp; dV = 0._dp                                           !Initialize grid cell volume and emission as zero
      do k=kwm,kwp                                                      !  (this ensures that corner points have
         do j=jwm,jwp                                                   !  zero emission and volume)
            do i=iwm,iwp               
               if (present(kappa_planck)) &                             !If the optional argument kappa_planck is present, store
                  kappa_planck(i,j,k) = kplanck(i,j,k)                  !  the local Planck-mean absorption coefficient for output
               if (ggmc.and.present(fixed_Ib)) then                     !For the GG model, if the optional argument fixed_Ib
                  Ib(i,j,k) = fixed_Ib(i,j,k)                           !  is present, use the provided array as the local 
                                                                        !  total blackbody intensity values
               else                                                     !Otherwise, compute the local total blackbody intensity
                  Ib(i,j,k) = sigrpi*T(i,j,k)**4._dp                    !  from the Stefan-Boltzmann law
               endif
               if (cell_type(i,j,k).eq.2) cycle                         !Skip corner points
               dV(i,j,k) = dx(i)*dy(j)*dz(k)                            !Compute local cell volume
               if (cell_type(i,j,k).eq.1) kIb(i,j,k) = pi*Ib(i,j,k)     !Emission at the wall points
               if (cell_type(i,j,k).eq.0) &
                  kIb(i,j,k) = fourpi*kplanck(i,j,k)*Ib(i,j,k)          !Emission at the points within the medium
               
               if ((cell_type(i,j,k).eq.1).and.losmc) &
                  kIb(i,j,k) = Ib(i,j,k)
               if ((cell_type(i,j,k).eq.0).and.losmc) &
                  kIb(i,j,k) = kplanck(i,j,k)*Ib(i,j,k)
            enddo
         enddo
      enddo
      if (compute_emission) emission = kIb                              !Store the local emission if requested
           
      !-----------------------------------------------------------------
      !Compute cumulative cell-integrated emission
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('mc_solution: &
                        &compute cumulative cell-integrated emission',3)
      ckIb = 0._dp; ckIb_y = 0._dp; ckIb_z = 0._dp                      !Initialize values for the summations that follow
      do i=iwm,iwp        
         do j=jwm,jwp               
            do k=kwm,kwp
               if (cell_type(i,j,k).eq.2) cycle                         !Skip corner cells
               ckIb_z(i,j,k) = ckIb_z(i,j,max(k-1,kwm)) + &             !Cumulative cell-integrated emission per unit area
                               kIb(i,j,k)*dz(k)                         !  [max(k-1,kwm) is used here to avoid accessing 
            enddo                                                       !  the k=kwm-1 position of the ckIb_z array]
            tkIb_z(i,j) = ckIb_z(i,j,kwp)                               !Total cell-integrated emission per unit area at i,j
            ckIb_y(i,j) = ckIb_y(i,max(j-1,jwm)) + tkIb_z(i,j)*dy(j)    !Cumulative cell-itegrated emission per unit length
         enddo
         tkIb_y(i) = ckIb_y(i,jwp)                                      !Total cell-integrated emission per unit length at i
         ckIb(i) = ckIb(max(i-1,iwm)) + tkIb_y(i)*dx(i)                 !Cumulative cell-integrated emission
      enddo
      tkIb = ckIb(iwp)                                                  !Domain-integrated emission [W]

      !-----------------------------------------------------------------
      !Zero out quantities
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('mc_solution: &
                                           &zero out quantities',3)
      Srad = 0._dp; qrad = 0._dp                                        !Radiative heat source and flux
      if (present(total_absorption)) total_absorption = 0._dp           !Domain-integrated absorption
      if (present(total_emission)) total_emission = 0._dp               !Domain-integrated emission
      if (present(total_loss)) total_loss = 0._dp                       !Domain-integrated net radiative loss
      
      run_loop: do irun=1,nrun                                          !Start loop over all statistical runs
         !--------------------------------------------------------------
         !Print loop status
         !--------------------------------------------------------------
         write(msg,'(i5)') irun
         msg = 'mc_solution: realization '//trim(adjustl(msg))//' began'
         if (debug_mode) call print_to_prompt(msg,3)
         
         !--------------------------------------------------------------
         !Define the seed (starts with the predefined value   
         !                 and changes with each new irun)
         !--------------------------------------------------------------
         if (debug_mode) call print_to_prompt('mc_solution: &
                                              &define seed',6)
         if (trim(which_mc_random).eq.'fortran') &
            call set_random_seed(irun*mc_seed)                               
         if (trim(which_mc_random).eq.'numrep') &
            idum = min(irun*mc_seed,-irun*mc_seed)                      !idum must be negative
         
         !--------------------------------------------------------------
         !Defining how many bundles each cell emits
         !(different treatments available)
         !--------------------------------------------------------------
         if (debug_mode) call print_to_prompt('mc_solution: &
                                &define bundles emitted by each cell',6)
         selectcase(trim(which_mc_emission))
         case ('fixed')                                                 !All cells emit the same number of rays
            n_emi = nr/((iwp-iwm+1)*(jwp-jwm+1)*(kwp-kwm+1))
         
         case('random')                                                 !Number of rays per cell dependent on the random number
            n_emi = 0                                                   !Zero out the array
            do n=1,nr
               R_x = get_mc_ran(idum)                                   !Define random number for x-position
               R_y = get_mc_ran(idum)                                   !Define random number for y-position
               R_z = get_mc_ran(idum)                                   !Define random number for z-position
               i = max(min(locate(ckIb/tkIb,R_x,iwp),iwp),iwm)          !Define emission x-position
               j = max(min(locate(ckIb_y(i,:)/&                         !Define emission y-position
                           tkIb_y(i),R_y,jwp),jwp),jwm)
               k = max(min(locate(ckIb_z(i,j,:)/&                       !Define emission z-position
                           tkIb_z(i,j),R_z,kwp),kwp),kwm)
               n_emi(i,j,k) = n_emi(i,j,k) + 1                          !Increase the counter for the number of bundles
            enddo                                                       !  for that emission poistion
            
         case('distributed')                                            !Number of bundles per cell distributed such that each
                                                                        !  has exactly the same energy
            n_emi = 0                                                   !Zero out the array
            do i=iwm,iwp
               do j=jwm,jwp
                  do k=kwm,kwp
                     n_emi(i,j,k) = ceiling(real(nr,dp)*kIb(i,j,k)*&    !Compute number of bundles for cell ijk (rounding up)
                                            dV(i,j,k)/tkIb,i8b)
                        
                  enddo
               enddo
            enddo
            
         case default
            call shutdown('mc_routines: which_mc_emission undefined')
         endselect

         !--------------------------------------------------------------
         !Computing absorption
         !--------------------------------------------------------------
         if (debug_mode) call print_to_prompt('mc_solution: &
                                              &compute absorption',6)
         bundle_abs = 0._dp                                             !Zero out absorbed energy for the run
         total_emission_run = 0._dp                                     !Zero out total emitted energy for the run
         x_grid_loop: do i=iwm,iwp
            y_grid_loop: do j=jwm,jwp
               z_grid_loop: do k=kwm,kwp
                  if (cell_type(i,j,k).eq.2) cycle z_grid_loop          !Cycle corner cells
                  if (n_emi(i,j,k).eq.0) cycle z_grid_loop              !Skip the loop if the cell does not emit
                  bundle_emi = kIb(i,j,k)*dv(i,j,k)/&                   !Compute the energy of the bundles emitted by the cell
                                 real(n_emi(i,j,k),dp)                  
                  emi_loop: do n=1,n_emi(i,j,k)                         !Loop over the bundles emitted by the cell
                     if (debug_mode) then
                        write(rstr,'(i5)') irun                         !Write into string the indexes of the run
                        write(istr,'(i5)') i                            !  the cell in the x-,
                        write(jstr,'(i5)') j                            !  y-,
                        write(kstr,'(i5)') k                            !  and z-directions,
                        write(nstr,'(i8)') n                            !  and the bundle count
                        write(msg,'(a,a,a,a,a,a,a,a,a,a,a)') &          !  Print status to screen
                           'mc_solution emi_loop: ',&
                           'run= ',trim(adjustl(rstr)),&
                           '; i= ',trim(adjustl(istr)),&
                           '; j= ',trim(adjustl(jstr)),&
                           '; k= ',trim(adjustl(kstr)),&
                           '; n= ',trim(adjustl(nstr))
                        call print_to_prompt(msg,9)
                     endif
                     
                     !Compute domain-integrated emission for the run
                     total_emission_run = total_emission_run + &
                                          bundle_emi

                     !For the LBL model, determine the
                     !wavelength of absorption
                     if (lblmc) then
                        R_eta = get_mc_ran(idum)                        !Random number for spectral position
                        if (cell_type(i,j,k).eq.0) eta = &              !Wavelength (medium)
                           get_mc_wavenumber(R_eta,T(i,j,k),xs(:,i,j,k))
                        if (cell_type(i,j,k).eq.1) eta = &              !Wavelength (black walls)
                           get_mc_wavenumber(R_eta,T(i,j,k),xs(:,i,j,k),&
                                             wall=.true.)
                     endif
             
                     !For the WSGG model, determine the
                     !gas of absorption
                     if (wsggmc) then
                        call shutdown('WSGG MC not ready')
                        R_eta = get_mc_ran(idum)                        !Random number for gray gas position
                        if (cell_type(i,j,k).eq.0) jgas = &             !Gray gas selection (medium)
                           get_mc_WSGGj(R_eta,T(i,j,k),p(i,j,k),&
                                        xs(:,i,j,k))
                        if (cell_type(i,j,k).eq.1) jgas = &             !Gray gas selection (black walls)
                           get_mc_WSGGj(R_eta,T(i,j,k),p(i,j,k),&
                                        xs(:,i,j,k),wall=.true.)
                     endif

                     !Define components of the direction vector ŝ
                     !(ŝ = s1 ê1 + s2 ê2 + s3 ê3)
                     R_theta = get_mc_ran(idum)                         !Random number for polar direction
                     R_psi = get_mc_ran(idum)                           !Random number for azimuthal direction
                     if (cell_type(i,j,k).eq.1) then                    !Vector components for wall cells (assumed isotropic and diffuse)
                        theta = dasin(dsqrt(R_theta))                   !Polar angle (measured from the normal to the wall)
                        psi = twopi*R_psi                               !Azimuthal angle
                        if ((i.eq.iwm).or.(i.eq.iwp)) then              !Direction vector components for grid points
                           s1 = dcos(theta); if (i.eq.iwp) s1 = -s1     !   at the east and west walls
                           s2 = dsin(theta)*dcos(psi)
                           s3 = dsin(theta)*dsin(psi)
                        endif
                        if (((j.eq.jwm).or.(j.eq.jwp))&
                             .and.(two_d.or.three_d)) then              !Direction vector components for grid points
                           s1 = dsin(theta)*dsin(psi)                   !   at the north and south walls
                           s2 = dcos(theta); if (j.eq.jwp) s2 = -s2
                           s3 = dsin(theta)*dcos(psi)
                        endif                   
                        if (((k.eq.kwm).or.(k.eq.kwp))&
                             .and.three_d) then                         !Direction vector components for grid points
                           s1 = dsin(theta)*dcos(psi)                   !   at the top and bottom walls
                           s2 = dsin(theta)*dsin(psi)
                           s3 = dcos(theta); if (k.eq.kwp) s3 = -s3
                        endif
                     else                                               !Vector components for cells within the medium
                        theta = dacos(1._dp - 2._dp*R_theta)            !Polar angle
                        psi = twopi*R_psi                               !Azimuthal angle
                        if (one_d) then                                 !Components of the direction vector
                           s1 = dcos(theta)                             !  for 1D geometries
                        else
                           s1 = dsin(theta)*dcos(psi)                   !Components of the direction vector
                           s2 = dsin(theta)*dsin(psi)                   !  for 2D and 3D geometries
                           s3 = dcos(theta)
                        endif
                     endif     
                     if (losmc)          s1 = 1._dp
                     if (one_d)          s2 = small                     !Correct s2 for 1D geometries
                     if (one_d.or.two_d) s3 = small                     !Correct s3 for 1D and 2D geometries              
   
                     !Define whether the direction vector is >1 (1)
                     !or <1 (-1) in each main orthogonal direction
                     idir = 1; if (s1.lt.0) idir = -1                   !x-direction
                     jdir = 1; if (s2.lt.0) jdir = -1                   !y-direction
                     kdir = 1; if (s3.lt.0) kdir = -1                   !z-direction
                     if (one_d)          jdir = 0                       !Correct jdir for 1d geometries
                     if (one_d.or.two_d) kdir = 0                       !Correct kdir for 1d and 2d geometries
               
                     !Follow the path of the bundle until it hits a wall
                     ipath = max(min(i,icp),icm)                        !Indexes of the first cell crossed by the bundle within
                     jpath = max(min(j,jcp),jcm)                        !  the medium (either the cell of emission, if within the
                     kpath = max(min(k,kcp),kcm)                        !  medium, or the cell adjacent to the wall cell of emission)
                     if ((cell_type(i,j,k).eq.0).and.&                  !Randomize the position of emission within
                        randomize_cell_emission) then                   !  the medium cell, if requested
                        R_x = get_mc_ran(idum)                          !  (in this case, assume that the local
                        R_y = get_mc_ran(idum)                          !  properties are uniform within the cell)
                        R_z = get_mc_ran(idum)
                        xpos = xf(i-1) + R_x*dx(i)                      !Emission position in the x direction
                        if (two_d.or.three_d) ypos = yf(j-1) + R_y*dy(j)!Emission position in the y direction
                        if (three_d) zpos = zf(k-1) + R_z*dz(k)         !Emission position in the z direction
                     else                                               !Otherwise, the originating position of the bundle
                        xpos = x(i); ypos = y(j); zpos = z(k)           !  is the cell center
                     endif
                     path_ijk = -10
                     path_counter = 1                                   !Set initial counter for the number of cells crossed by the bundle
                     path_ijk(path_counter,:) = (/ ipath, jpath, kpath/)!Cell indexes for the first cell in the path
                     follow_ray_loop: do while &                        !Loop over the path until a wall point is reached
                                    (cell_type(ipath,jpath,kpath).eq.0)
                        
                        !Compute optical thickness
                        xfpos = xf(ipath + min(idir,0))                 !Position of the downwind face in the x-, y- and z-directions
                        yfpos = yf(jpath + min(jdir,0))                 ! (note that the cell face position is indexed such that the  
                        zfpos = zf(kpath + min(kdir,0))                 !  grid point has the same index to its downwind face)
                        dface(1) = (xfpos - xpos)/s1                    !Compute distances from the incoming point and the faces normal 
                        if (two_d.or.three_d) &                         !  to the x-, y- and z-directions (1 = x, 2 = y, 3 = z)
                           dface(2) = (yfpos - ypos)/s2
                        if (three_d) dface(3) = (zfpos - zpos)/s3
                        ds = minval(dface(1:3))                         !The actual path length is the smallest of these distances
                        if (ggmc) kappa = kplanck(ipath,jpath,kpath)    !Local absorption coefficient for the GG model
                        if (wsggmc) &
                           call wsgg_properties(kappa,dummy_dp,jgas,&   !Local absorption coefficient for the WSGG model
                              T(ipath,jpath,kpath),&
                              xs(:,ipath,jpath,kpath),&
                              p(ipath,jpath,kpath),which_wsgg_model)
                        if (lblmc) kappa = &                            !Local absorption coefficient for the LBL method
                           get_mc_keta(eta,T(ipath,jpath,kpath),&
                                       xs(:,ipath,jpath,kpath))                                 
                        tau(path_counter) = kappa*ds                    !Optical thickness of the cell

                        !Find next cell intersected by the ray
                        next_dir = find_imin(dface,1,3,1)               !Find which one the downwind cell faces
                        if (next_dir.eq.1) ipath = ipath + idir         !  is intersected by the ray 
                        if (next_dir.eq.2) jpath = jpath + jdir         !  and update the path indexes
                        if (next_dir.eq.3) kpath = kpath + kdir         !  accordingly
                        xpos = xpos + ds*s1                             !Update position of the ray within the cell in the x,
                        ypos = ypos + ds*s2                             !  y and
                        zpos = zpos + ds*s3                             !  z-directions
                        path_counter = path_counter + 1                 !Add to the counter for the number of cells in the path
                        path_ijk(path_counter,1:3) = &                  !Store the current cell indexes to the path_ijk array
                           (/ ipath, jpath, kpath/)

                     enddo follow_ray_loop
                     path_size = path_counter                           !Total number of cells in the path

                     !Absorb the bundle
                     selectcase(trim(which_mc_absorption))
                     case('standard')                                   !Apply the MC method as in Modest's book
                        R_abs = get_mc_ran(idum)                        !Random number for absorption
                        ckds = 0._dp                                    !Zero out sum for the cumulative optical thickness
                        do path_counter=1,path_size-1                   !Find out where bundle is absorbed 
                           ckds = ckds + tau(path_counter)              !Cumulative optical thickness
                           if (ckds.gt.log(1._dp/R_abs)) exit           !Eq. (21.19), Modest 2013
                        enddo
                        ipath = path_ijk(path_counter,1)                !x-Index of the grid cell where the bundle is absorbed
                        jpath = path_ijk(path_counter,2)                !y-Index of the grid cell where the bundle is absorbed
                        kpath = path_ijk(path_counter,3)                !z-Index of the grid cell where the bundle is absorbed
                        bundle_abs(ipath,jpath,kpath) = &               !Absorb the bundle
                           bundle_abs(ipath,jpath,kpath) + bundle_emi                     

                     case('energy_partitioning')                        !Energy-partitioning method
                        abs0 = bundle_emi                               !The energy at the start of the path is the bundle energy
                        do path_counter=1,path_size-1                   !Trace the bundle from the first cell all the way to the
                                                                        !  last cell within the medium
                           ipath = path_ijk(path_counter,1)             !x-index of the grid cell
                           jpath = path_ijk(path_counter,2)             !y-index of the grid cell
                           kpath = path_ijk(path_counter,3)             !z-index of the grid cell
                           bundle_abs(ipath,jpath,kpath) = &            !Add to the cell absorption
                              bundle_abs(ipath,jpath,kpath) + &
                                 abs0*(1._dp - dexp(-tau(path_counter)))
                           abs0 = abs0*dexp(-tau(path_counter))         !Update the incoming energy at the next cell
                        enddo
                        ipath = path_ijk(path_size,1)                   !x-index of the last cell (wall) along the path
                        jpath = path_ijk(path_size,2)                   !y-index of the last cell (wall) along the path
                        kpath = path_ijk(path_size,3)                   !z-index of the last cell (wall) along the path
                        bundle_abs(ipath,jpath,kpath) = &               !All energy not absorbed in the path is absorbed
                           bundle_abs(ipath,jpath,kpath) + abs0         !  by the wall point at the end of the path

                     case default
                        call shutdown('mc_routines: &
                                      &which_mc_absorption undefied')                                 
                     endselect
                  enddo emi_loop
               enddo z_grid_loop
            enddo y_grid_loop
         enddo x_grid_loop
   
         !--------------------------------------------------------------
         !Compute heat source and heat flux at each cell
         !--------------------------------------------------------------
         if (debug_mode) call print_to_prompt('mc_solution: &
                                 &compute local heat source and flux',6)
         
         total_absorption_run = 0._dp                                   !Zero out the domain-integrated absorption for the run
         avcoeff = 1._dp/real(irun,dp)                                  !Coefficient for computing the averages
         
         !Compute x-normal wall radiative heat fluxes
         do j=jcm,jcp
            do k=kcm,kcp
               i = iwm                                                  !West wall
               total_absorption_run = total_absorption_run + &          !Add to domain-integrated absorption
                                       bundle_abs(i,j,k)                !  for the run
               qrad_lcl = kIb(i,j,k) - bundle_abs(i,j,k)/dV(i,j,k)      !Radiative heat flux for this run
               qrad(1,i,j,k) = (1._dp - avcoeff)*qrad(1,i,j,k) + &      !Update the average radiative heat flux over all runs
                                 avcoeff*qrad_lcl
                                 
               i = iwp                                                  !East wall
               total_absorption_run = total_absorption_run + &          !Add to domain-integrated absorption
                                       bundle_abs(i,j,k)                !  for the run
               qrad_lcl = bundle_abs(i,j,k)/dV(i,j,k) - kIb(i,j,k)      !Radiative heat flux for this run
               qrad(1,i,j,k) = (1._dp - avcoeff)*qrad(1,i,j,k) + &      !Update the average radiative heat flux over all runs
                                 avcoeff*qrad_lcl                       
            enddo
         enddo

         !Compute y-normal wall radiative heat fluxes
         if (two_d.or.three_d) then
            do i=icm,icp
               do k=kcm,kcp
                  j = jwm                                               !South wall
                  total_absorption_run = total_absorption_run + &       !Add to domain-integrated absorption
                                          bundle_abs(i,j,k)             !  for the run
                  qrad_lcl = kIb(i,j,k) - bundle_abs(i,j,k)/dV(i,j,k)   !Radiative heat flux for this run
                  qrad(2,i,j,k) = (1._dp - avcoeff)*qrad(2,i,j,k) + &   !Update the verage radiative heat flux over all runs
                                    avcoeff*qrad_lcl
                                    
                  j = jwp                                               !North wall
                  total_absorption_run = total_absorption_run + &       !Add to domain-integrated absorption
                                          bundle_abs(i,j,k)             !  for the run
                  qrad_lcl = bundle_abs(i,j,k)/dV(i,j,k) - kIb(i,j,k)   !Radiative heat flux for this run
                  qrad(2,i,j,k) = (1._dp - avcoeff)*qrad(2,i,j,k) + &   !Update the verage radiative heat flux over all runs
                                    avcoeff*qrad_lcl                       
               enddo
            enddo
         endif
   
         !Compute y-normal wall radiative heat fluxes
         if (three_d) then
            do i=icm,icp
               do j=jcm,jcp
                  k = kwm                                               !Bottom wall
                  total_absorption_run = total_absorption_run + &       !Add to domain-integrated absorption
                                          bundle_abs(i,j,k)             !  for the run
                  qrad_lcl = kIb(i,j,k) - bundle_abs(i,j,k)/dV(i,j,k)   !Radiative heat flux for this run
                  qrad(3,i,j,k) = (1._dp - avcoeff)*qrad(3,i,j,k) + &   !Update the verage radiative heat flux over all runs
                                    avcoeff*qrad_lcl
                                    
                  k = kwp                                               !Top wall
                  total_absorption_run = total_absorption_run + &       !Add to domain-integrated absorption
                                          bundle_abs(i,j,k)             !  for the run
                  qrad_lcl = bundle_abs(i,j,k)/dV(i,j,k) - kIb(i,j,k)   !Radiative heat flux for this run
                  qrad(3,i,j,k) = (1._dp - avcoeff)*qrad(3,i,j,k) + &   !Update the verage radiative heat flux over all runs
                                    avcoeff*qrad_lcl                       
               enddo
            enddo
         endif
 
         !Compute radiative source in the medium
         do i=icm,icp                                                   !Loop over grid points within the medium
            do j=jcm,jcp
               do k=kcm,kcp
                  total_absorption_run = total_absorption_run + &       !Add to domain-integrated absorption
                                          bundle_abs(i,j,k)             !  for the run
                  Srad_lcl = bundle_abs(i,j,k)/dV(i,j,k) - kIb(i,j,k)   !Radiative heat source for this irun
                  Srad(i,j,k) = (1._dp - avcoeff)*Srad(i,j,k) + &       !Update the average radiative heat source
                                 avcoeff*Srad_lcl
               enddo
            enddo
         enddo

         !Update domain-integrated quantities
         if (present(total_absorption)) &
            total_absorption = (1._dp - avcoeff)*total_absorption + &
                               avcoeff*total_absorption_run
         if (present(total_emission)) &
            total_emission = (1._dp - avcoeff)*total_emission + &
                               avcoeff*total_emission_run
                               
         !--------------------------------------------------------------
         !Check energy conservation, if requested
         !--------------------------------------------------------------
         if (debug_mode) then
            call print_to_prompt('mc_solution: &
                                 &check energy conservation',6)
            write(msg,'(a,f12.4)') &
               'Total emission [kW]: ',total_emission_run/1000._dp
            call print_to_prompt(msg,9)   
            write(msg,'(a,f12.4)') &
               'Total absorption [kW]: ',total_absorption_run/1000._dp
            call print_to_prompt(msg,9)
         endif  
         net_run = total_emission_run - total_absorption_run        
         if (check_mc_energy.and.&
            (dabs(net_run)/total_emission_run.gt.1.e-10)) &
               call shutdown ('mc_solution: Total energy not conserved')

      enddo run_loop

      !-----------------------------------------------------------------
      !Compute the domain-integrated net-radiation loss
      !-----------------------------------------------------------------
      if (present(total_loss)) then
         if (debug_mode) call print_to_prompt('mc_solution: &
                                              &compute total loss',3)
         total_loss = 0._dp
         do i=icm,icp                                                   !Loop over grid points within the medium
            do j=jcm,jcp
               do k=kcm,kcp
                  total_loss = total_loss - Srad(i,j,k)*dV(i,j,k) 
               enddo
            enddo
         enddo
      endif

      !-----------------------------------------------------------------
      !Compute the local absorption
      !-----------------------------------------------------------------
      if (compute_absorption) then
         if (debug_mode) &
            call print_to_prompt('mc_solution: &
                                 &compute local absorption',3)
         absorption = 0._dp
         do i=icm,icp                                                   !Loop over grid points within the medium
            do j=jcm,jcp
               do k=kcm,kcp
                  absorption(i,j,k) = Srad(i,j,k) + kIb(i,j,k)
               enddo
            enddo
         enddo
      endif

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      if (debug_mode) &
         call print_to_prompt('mc_solution: deallocate arrays',3)
      deallocate(n_emi,bundle_abs)
      deallocate(cell_type,path_ijk)
      deallocate(ckIb,ckIb_y,ckIb_z)
      deallocate(dV,dx,dy,dz)
      deallocate(Ib,kplanck,kIb,tau)
      deallocate(tkIb_y,tkIb_z)

   endsubroutine mc_solution

endmodule mc_routines
