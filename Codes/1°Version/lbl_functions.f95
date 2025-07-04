!#######################################################################
!Functions and subroutines related to the 
!radiative transfer solution via LBL integration
!#######################################################################
module lbl_functions

   !====================================================================
   !Modules & Misc
   !====================================================================
   use lbl_parameters
   implicit none

contains

   !====================================================================
   !This subroutine checks if there is any conflict between
   !the specified parameters related to the LBL method
   !====================================================================
   subroutine check_lbl_parameters

      use comp_functions, only: shutdown
      implicit none
      if (lbl_kappa_input.and.lbl_xceta_input) &
         call shutdown('True lbl_kappa_input and lbl_xceta_input')
      if (lbl_kappa_input.and.lbl_ceta_input) &
         call shutdown('True lbl_kappa_input and lbl_ceta_input')
      if (lbl_ceta_input.and.lbl_xceta_input) &
         call shutdown('True lbl_ceta_input and lbl_xceta_input')

   endsubroutine check_lbl_parameters

   !====================================================================
   !Subroutine to open each external LBL data file
   !====================================================================
   subroutine open_lbl_data

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,get_file_unit,shutdown
      use global_parameters, only: number_of_species
      implicit none
      character(20) :: form_string
      character(200) :: file_name
      integer :: ict,iit,isp,ixs
      integer :: nsp,nxs

      !-----------------------------------------------------------------
      !Set up form string
      !-----------------------------------------------------------------
      form_string = 'formatted'
      if (lbl_binary_data) form_string = 'unformatted'
      
      !-----------------------------------------------------------------
      !Assign an unit value for each of the files (stored in the 
      !lbl_data_unit array) and open each of them
      !-----------------------------------------------------------------
      if (allocated(lbl_data_unit)) deallocate(lbl_data_unit)           !Deallocate lbl_data_unit if it is already allocated
      nsp = number_of_species                                           !Surrogate name for the total number of species
      do iit=1,2                                                        !iit=1 -> Allocation; iit=2 -> Get and open unit
         ict = 0                                                        !Initialize unit counter
         do isp=1,nsp                                                   !Loop over all species
            nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle                !Skip species for which there is no LBL data file
            do ixs=1,nxs                                                !Loop over all discrete mole fraction values
               if (lbl_data_file(ixs,isp).eq.'null') cycle              !Skip mole fraction values for wehich there is no LBL data file
               ict = ict + 1                                            !Add to unit counter
               if (iit.eq.1) cycle                                      !Skip getting and opening the unit for allocation step
               file_name = trim(lbl_data_prefix)//&                     !Define the full file name
                  trim(lbl_data_file(ixs,isp))//trim(lbl_data_ext)
!               file_name = adjustl(trim(file_name))
               call CheckFileExists(file_name)                          !Check if the file exists
               lbl_data_unit(ict) = get_file_unit()                     !Get file unit
               open(unit=lbl_data_unit(ict),file=file_name,&            !Open the unit
                    form=trim(form_string),action='read')
            enddo
         enddo
         deta_if: if (trim(lbl_deta_file).ne.'null') then               !Prepare and read the wavenumber spacing data from the 
                                                                        !  extenal file, if one has been provided
            ict = ict + 1                                               !Add again to the unit counter
            if (iit.eq.1) exit deta_if                                     !Skip getting and opening the unit for allocation step
            lbl_data_unit(ict) = get_file_unit()                        !Get file unit
            call CheckFileExists(lbl_deta_file)                         !Check if the file exists
            open(unit=lbl_data_unit(ict),file=lbl_deta_file,&           !Open the unit
                 form=trim(form_string),action='read')
         endif deta_if
         if (iit.eq.1) allocate(lbl_data_unit(ict))                     !Allocate array with all the units for reading the LBL data files
      enddo
   
   endsubroutine open_lbl_data

   !====================================================================
   !Subroutine to finish the reading of the external LBL data files
   !====================================================================
   subroutine close_lbl_data
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use global_parameters, only: number_of_species
      implicit none
      integer :: ict,isp,ixs
      integer :: nsp,nxs
      
      !-----------------------------------------------------------------
      !Loop for closing all units
      !-----------------------------------------------------------------
      nsp = number_of_species                                           !Surrogate name for the total number of species
      ict = 0                                                           !Initialize unit counter
      do isp=1,nsp                                                      !Loop over all species
         nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle                   !Skip species for which there is no LBL data file
         do ixs=1,nxs                                                   !Loop over all discrete mole fraction values
            if (lbl_data_file(ixs,isp).eq.'null') cycle                 !Skip mole fraction values for wehich there is no LBL data file
            ict = ict + 1                                               !Update unit counter
            close(lbl_data_unit(ict))                                   !Close current unit
         enddo
      enddo
      if (trim(lbl_deta_file).ne.'null') then                           !Also close the unit corresponding to the data file with the
         ict = ict + 1                                                  !  wavenumber interval size, if this is to be used
         close(lbl_data_unit(ict))
      endif

   endsubroutine close_lbl_data

   !====================================================================
   !Subroutine to rewind all external LBL data files
   !====================================================================
   subroutine rewind_lbl_data

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use global_parameters, only: number_of_species
      implicit none
      integer :: ict,isp,ixs
      integer :: nsp,nxs
      
      !-----------------------------------------------------------------
      !Loop for rewinding all units
      !-----------------------------------------------------------------
      nsp = number_of_species                                           !Surrogate name for the total number of species
      ict = 0                                                           !Initialize unit counter
      do isp=1,nsp                                                      !Loop over all species
         nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle                   !Skip species for which there is no LBL data file
         do ixs=1,nxs                                                   !Loop over all discrete mole fraction values
            if (lbl_data_file(ixs,isp).eq.'null') cycle                 !Skip mole fraction values for wehich there is no LBL data file
            ict = ict + 1                                               !Add to the unit counter
            rewind(lbl_data_unit(ict))                                  !Rewind current unit
         enddo
      enddo
      if (trim(lbl_deta_file).ne.'null') then                           !Rewind the unit corresponding to the data file with the
         ict = ict + 1                                                  !  wavenumber interval size, if this is to be used
         rewind(lbl_data_unit(ict))
      endif

   endsubroutine rewind_lbl_data

   !====================================================================
   !Subroutine to read the one line (i.e., for one discrete wavenumber)
   !of external LBL data files
   !====================================================================
   subroutine read_lbl_data

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,shutdown
      use global_parameters, only: number_of_species
      use precision_parameters, only: dp,small
      implicit none
      character(300) :: msg
      integer :: dummy_int
      integer :: ict,ierr,ilbl,isp,itg,ixs
      integer :: nlbl,nsp,ntg,nxs
      logical :: skip_check
      real(dp),allocatable,dimension(:) :: deta_0,deta_1,xeta_0,xeta_1
      real(dp),allocatable,dimension(:,:,:) :: acs_d,acs_u
      real(dp) :: cm_to_m,eta_d,eta_u,tol,x0,x1,xeta_e,xeta_p,xeta_w
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      !Set parameters
      tol = 0.01_dp                                                     !Tolerance for defining if wavenumber positions
                                                                        !  match between different input files
      cm_to_m = 1._dp; if (lbl_data_cm) cm_to_m = 100._dp
      
      !Surrogate names
      nlbl = lbl_nlines                                                 !Surrogate name for the total number of LBL lines
      nsp = number_of_species                                           !Surrogate name for the total number of species
      
      !-----------------------------------------------------------------
      !Main loop for binary data
      !-----------------------------------------------------------------
      if (lbl_binary_data) then
         !Allocate arrays
         allocate(xeta_0(nlbl),xeta_1(nlbl))
         allocate(deta_0(nlbl),deta_1(nlbl))
      
         !Reading loop
         ict = 0                                                        !Zero out unit counter
         spec_loop_bin: do isp=1,nsp                                    !Loop over all species
            nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle                !Skip species for which there is no LBL data file
            xspec_loop_bin: do ixs=1,nxs                                !Loop over all discrete mole fraction values
               if (lbl_data_file(ixs,isp).eq.'null') cycle              !Skip mole fraction values for which there is no LBL data file
               ict = ict + 1                                            !Add to the unit counter
               ntg = lbl_data_ntg(ixs,isp)                              !Surrogate name for number of discrete temperature values      
               read(lbl_data_unit(ict)) dummy_int,skip_check            !Skip first line of data
               read(lbl_data_unit(ict)) xeta_0                          !Read wavenumber positions
               read(lbl_data_unit(ict)) deta_0                          !Read wavenumber intervals
               read(lbl_data_unit(ict)) lbl_data(:,:,ixs,isp)           !Read data

               !Check if the wavenumber positions and intervals match
               !from file to file
               if (skip_check) cycle xspec_loop_bin 
               if (ict.eq.1) then
                  xeta_1 = xeta_0
                  deta_1 = deta_0
                  cycle xspec_loop_bin
               endif
               line_loop_bin: do ilbl=1,nlbl
                  if (dabs(xeta_1(ilbl)-xeta_0(ilbl))/&
                     (xeta_1(ilbl)+small).gt.tol) then
                     msg = 'Mismatch for LBL input data file '//&
                           trim(lbl_data_file(ixs,isp))
                     call shutdown(msg)
                  endif 
               enddo line_loop_bin
               xeta_1 = xeta_0; deta_1 = deta_0
            enddo xspec_loop_bin
         enddo spec_loop_bin
         lbl_xeta = xeta_0; lbl_deta = deta_0
         
         !Deallocate arrays
         deallocate(deta_0,deta_1,xeta_0,xeta_1)
      endif

      !-----------------------------------------------------------------
      !Main loop for non-binary data
      !-----------------------------------------------------------------
      if (.not.lbl_binary_data) then
         !Allocate arrays
         allocate(acs_d(lbl_ntg,lbl_nxs,lbl_nsp),stat=ierr)
         call CheckMemAlloc('acs_d',ierr)
         allocate(acs_u(lbl_ntg,lbl_nxs,lbl_nsp),stat=ierr)
         call CheckMemAlloc('acs_u',ierr)
      
         !Correct bounds for the reading of the LBL data files
         if (trim(lbl_data_averaging).eq.'arithmetic') nlbl = nlbl + 1
         if (trim(lbl_data_averaging).eq.'geometric')  nlbl = nlbl + 1
         if (trim(lbl_data_averaging).eq.'upwind')     nlbl = nlbl + 1
         if (trim(lbl_data_averaging).eq.'none')       nlbl = nlbl         
         
         !Reading loop
         acs_u = 0._dp                                                  !Initialze input array
         line_loop: do ilbl=1,nlbl                                      !Loop over all LBL lines
            ict = 0                                                     !Zero out unit counter
            spec_loop: do isp=1,nsp                                     !Loop over all species
               nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle             !Skip species for which there is no LBL data file
               xspec_loop: do ixs=1,nxs                                 !Loop over all discrete mole fraction values
                  if (lbl_data_file(ixs,isp).eq.'null') cycle           !Skip mole fraction values for which there is no LBL data file
                  ict = ict + 1                                         !Add to the unit counter
                  ntg = lbl_data_ntg(ixs,isp)                           !Surrogate name for number of discrete temperature values      
                  read(lbl_data_unit(ict),*) &                          !Read one line of the LBL data file for species index isp and
                     x0,(acs_u(itg,ixs,isp),itg=1,ntg)                  !  mole fraction index ixs
                  if (ict.eq.1) x1 = x0                                 !Skip comparison between wavenumber positions for the first file
                  if (dabs(x1-x0)/(x1+small).gt.tol) then               !Check if the wavenumber position read between two consecutive
                     msg = 'Mismatch for LBL input data file '//&       !  files is the same; if not, stop the calculation and print
                           trim(lbl_data_file(ixs,isp))                 !  an error message
                     call shutdown(msg)
                  endif                              
               enddo xspec_loop
            enddo spec_loop  
            
            !Convert to 1/m, if needed
            x0 = x0*cm_to_m
            
            !Read wavenumber interval, if requested            
            if (trim(lbl_data_averaging).eq.'none') then
               if (trim(lbl_deta_file).ne.'null') then
                  ict = ict + 1                                         !Add to the unit counter
                  read(lbl_data_unit(ict),*) lbl_deta(ilbl)             !Read one line of the LBL data file
               else
                  xeta_w = xeta_p; xeta_p = xeta_e; xeta_e = x0                  
                  if (ilbl.eq.2) &
                     lbl_deta(ilbl-1) = (xeta_e - xeta_p)
                  if (ilbl.gt.2) &
                     lbl_deta(ilbl-1) = 0.5_dp*(xeta_e - xeta_w)
                  if (ilbl.eq.nlbl) &
                     lbl_deta(ilbl) = (xeta_e - xeta_p)
               endif
            endif
            
            if (ilbl.eq.1) then
               eta_d = x0; acs_d = acs_u; cycle line_loop
            endif
   
            !Define absorption cross-section for the line
            do isp=1,nsp
               do ixs=1,lbl_data_nxs(isp)
                  do itg=1,lbl_data_ntg(ixs,isp)
                     if (trim(lbl_data_averaging).eq.'arithmetic') &
                        lbl_data(ilbl-1,itg,ixs,isp) = &
                           0.5_dp*(acs_u(itg,ixs,isp) + &
                              acs_d(itg,ixs,isp))
                     if (trim(lbl_data_averaging).eq.'geometric') &
                        lbl_data(ilbl-1,itg,ixs,isp) = &
                           dsqrt(acs_u(itg,ixs,isp)*acs_d(itg,ixs,isp))
                     if (trim(lbl_data_averaging).eq.'upwind') &
                        lbl_data(ilbl-1,itg,ixs,isp) = &
                           acs_u(itg,ixs,isp)
                     if (trim(lbl_data_averaging).eq.'none') &
                        lbl_data(ilbl,itg,ixs,isp) = &
                           acs_u(itg,ixs,isp)
                  enddo
               enddo
            enddo
            acs_d = acs_u

            !Defining the spectral position
            eta_d = eta_u; eta_u = x0
            if (trim(lbl_data_averaging).eq.'arithmetic') &
               lbl_xeta(ilbl-1) = 0.5_dp*(eta_u + eta_d)
            if (trim(lbl_data_averaging).eq.'geometric') &
               lbl_xeta(ilbl-1) = 0.5_dp*(eta_u + eta_d)
            if (trim(lbl_data_averaging).eq.'upwind') &  
               lbl_xeta(ilbl-1) = eta_u
            if (trim(lbl_data_averaging).eq.'none') &  
               lbl_xeta(ilbl) = eta_u
         enddo line_loop

         !Deallocate arrays
         deallocate(acs_d,acs_u)
      endif

   endsubroutine read_lbl_data

   !====================================================================
   !Subroutine to load all the needed data for the LBL calculations
   !from external files
   !====================================================================
   subroutine load_lbl_data
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,CheckMemAlloc,&
                                dprint,get_file_unit,shutdown
      use constants, only: atm_pa,avogadro,Ru
      use global_parameters, only: number_of_species
      implicit none
      character(300) :: msg
      integer :: ierr
      integer :: ict,isp,itg,ixs
      integer :: nsp,nxs
      integer,allocatable,dimension(:) :: nlines
      
      !-----------------------------------------------------------------
      !Check if the data is already loaded
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Check if the data is already loaded')
      if (lbl_data_ready) then
         call rewind_lbl_data
         return
      endif
      lbl_data_ready = .true.
   
      !-----------------------------------------------------------------
      !Check the LBL parameters for inconsistencies
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Check the LBL parameters')
      call check_lbl_parameters
   
      !-----------------------------------------------------------------
      !Open LBL units
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Open LBL units')
      call open_lbl_data

      !-----------------------------------------------------------------
      !Determine the total number of spectral lines in the data files
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Determine number of spectral lines')
      nsp = number_of_species                                           !Surrogate name for the total number of species
      allocate(nlines(size(lbl_data_unit)))                             !Allocate array with the number of lines of each LBL data file
      ict = 0                                                           !Reinitialize unit counter
      do isp=1,nsp                                                      !Loop over all species
         nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle                   !Skip species for which there is no LBL data file
         do ixs=1,nxs                                                   !Loop over all discrete mole fraction values
            if (lbl_data_file(ixs,isp).eq.'null') cycle                 !Skip mole fraction values for wehich there is no LBL data file
            ict = ict + 1                                               !Update unit counter
            if (lbl_binary_data) then
               read(lbl_data_unit(ict)) nlines(ict)
            else
               nlines(ict) = 0                                          !Initialize line counter
               do                                                       !Count number of lines in the LBL data file
                  read(lbl_data_unit(ict),*,iostat=ierr)
                  if (ierr.ne.0) exit; nlines(ict) = nlines(ict) + 1
               enddo
            endif
            if (ict.eq.1) cycle                                         !Skip comparison between number of lines for the first file
            if (nlines(ict).ne.nlines(ict-1)) then                      !Check if two consecutive LBL data files have the same 
               msg = 'Number of lines mismatch for LBL '//&             !  number of lines; if not, stop the calculation and 
                     'input data file '//trim(lbl_data_file(ixs,isp))   !  print an error message
               call shutdown(msg)
            endif            
         enddo
      enddo
      
      !Assign final total number of lines for all files
      if (trim(lbl_data_averaging).eq.'arithmetic') &
         lbl_nlines = nlines(1) - 1
      if (trim(lbl_data_averaging).eq.'geometric') &
         lbl_nlines = nlines(1) - 1
      if (trim(lbl_data_averaging).eq.'upwind') &
         lbl_nlines = nlines(1) - 1
      if (trim(lbl_data_averaging).eq.'none') &
         lbl_nlines = nlines(1)
      deallocate(nlines)                                                !Deallocate nlines array

      !-----------------------------------------------------------------
      !Define arrays sizes
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Define arrays sizes')
      lbl_nsp = nsp
      lbl_nxs = maxval(lbl_data_nxs)
      lbl_ntg = maxval(lbl_data_ntg)

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Allocate arrays')
      if (allocated(lbl_data)) deallocate(lbl_data)
      allocate(lbl_data(lbl_nlines,lbl_ntg,lbl_nxs,lbl_nsp),stat=ierr)
      call CheckMemAlloc('lbl_data',ierr)
      if (allocated(lbl_deta)) deallocate(lbl_deta)
      allocate(lbl_deta(lbl_nlines),stat=ierr)
      call CheckMemAlloc('lbl_deta',ierr)
      if (allocated(lbl_xeta)) deallocate(lbl_xeta)
      allocate(lbl_xeta(lbl_nlines),stat=ierr)
      call CheckMemAlloc('lbl_xeta',ierr)
      lbl_data = 0._dp; lbl_deta = 0._dp; lbl_xeta = 0._dp

      !-----------------------------------------------------------------
      !Rewind the LBL data files
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Rewind the LBL data files')
      call rewind_lbl_data
      
      !-----------------------------------------------------------------
      !Read the LBL data files
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Read the LBL data files')
      call read_lbl_data
      
      !-----------------------------------------------------------------
      !Convert the LBL data into data for
      !pressure absorption coefficient
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Convert the LBL data into &
                  &data for pressure absorption coefficient')
      if (lbl_xceta_input.or.lbl_ceta_input) then
         do isp=1,lbl_nsp
            do ixs=1,lbl_data_nxs(isp)                                  !Loop over all discrete mole fraction values
               if (lbl_xceta_input) lbl_data(:,:,ixs,isp) = &           !Convert from x*ceta to ceta
                  lbl_data(:,:,ixs,isp)/lbl_data_xs(ixs,isp)            !  if necessary
               do itg=1,lbl_data_ntg(ixs,isp)                           !Loop over all discrete temperature values
                  lbl_data(:,itg,ixs,isp) = &                           !Convert from absorption cross-section
                     avogadro*lbl_data(:,itg,ixs,isp)*atm_pa/&          !  (ceta) to absortpion coefficient,
                        (Ru*lbl_data_tg(itg,ixs,isp)*1.e4_dp)           !  if necessary
               enddo
            enddo            
         enddo
      endif
      
      !-----------------------------------------------------------------
      !Close the LBL data files
      !-----------------------------------------------------------------
      call dprint('load_lbl_data: Close the LBL data files')
      call close_lbl_data

   endsubroutine load_lbl_data

   subroutine deallocate_lbl_data
      if (allocated(lbl_data)) deallocate(lbl_data)
      if (allocated(lbl_deta)) deallocate(lbl_deta)
      if (allocated(lbl_xeta)) deallocate(lbl_xeta)
   
   endsubroutine deallocate_lbl_data

   !====================================================================
   !Function to compute the spectral absorption coefficient
   !of species with index isp at a given temperature tmp and mole
   !fraction xsp. If kappa_p = .true., compute the pressure-based
   !coefficient instead (by default, kappa_p = .false.)
   !====================================================================
   real(dp) function lbl_kappa_func(xsp,tmp,pres,isp,ilbl,kappa_p)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,shutdown
      use math_functions, only: locate
      use precision_parameters, only: dp,small
      integer,intent(in) :: ilbl,isp
      integer :: itg,ixs,ntg,nxs
      integer :: lxs,ltg,uxs,utg
      real(dp),intent(in) :: pres,tmp,xsp
      real(dp) :: Rtg,Rxs,Qtg,Qxs
      logical,intent(in),optional :: kappa_p
      logical :: pbased
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
!      call dprint('lbl_kappa_func: Preparatory procedures')
      pbased = .false.; if (present(kappa_p)) pbased = kappa_p
      nxs = lbl_data_nxs(isp)
      
      !-----------------------------------------------------------------
      !Interpolation
      !-----------------------------------------------------------------   
!      call dprint('lbl_kappa_func: Interpolation')
      
      !Upper and lower indexes for mole fraction
      lxs = locate(lbl_data_xs(:,isp),xsp,nxs)                          !Locate lower index
      lxs = max(1,min(nxs-1,lxs)); uxs = min(lxs+1,nxs)                 !Correct lower index, compute upper index
      if (xsp.lt.lbl_data_xs(1,isp))    uxs = lxs                       !This is to prevent extrapolations
      if (xsp.gt.lbl_data_xs(nxs,isp))  lxs = uxs                       !  (instead, simply take the value at the extreme)

      !Initial values for the interpolation on mole fraction
      Qxs = 0._dp
      Rxs = (xsp - lbl_data_xs(lxs,isp))/&
         (lbl_data_xs(uxs,isp) - lbl_data_xs(lxs,isp) + small)
      if (lxs.eq.uxs) Rxs = 0._dp                                       !If only one mole fraction value is provided,
                                                                        !  do not interpolate in mole fraction
                                                                        !  (recall: we are interpolating for kp, not k)

      do ixs=lxs,uxs
         !Upper and lower indexes for temperature
         ntg = lbl_data_ntg(ixs,isp)
         ltg = locate(lbl_data_tg(:,ixs,isp),tmp,ntg)
         ltg = max(1,min(ntg-1,ltg)); utg = min(ltg+1,ntg)
         if (tmp.lt.lbl_data_tg(1,ixs,isp))    utg = ltg
         if (tmp.gt.lbl_data_tg(ntg,ixs,isp))  ltg = utg

         !Initial values for the interpolation on temperature
         Qtg = 0._dp
         Rtg = (tmp - lbl_data_tg(ltg,ixs,isp))/&
            (lbl_data_tg(utg,ixs,isp) - lbl_data_tg(ltg,ixs,isp) + small)
         if (ltg.eq.utg) Rtg = 0._dp                                    !If only one temperature value is provided, 
                                                                        !  do not interpolate in temperature

         do itg=ltg,utg
            !Interpolate on temperature
            Rtg = 1._dp - Rtg
            Qtg = Qtg + Rtg*lbl_data(ilbl,itg,ixs,isp)
         enddo
         
         !Interpolate on mole fraction
         Rxs = 1._dp - Rxs
         Qxs = Qxs + Rxs*Qtg
      enddo

      !-----------------------------------------------------------------
      !Final value
      !-----------------------------------------------------------------
!      call dprint('lbl_kappa_func: Final value')
      if (.not.pbased) Qxs = Qxs*xsp*pres
      lbl_kappa_func = Qxs

   endfunction lbl_kappa_func
   
!   !====================================================================
!   !Subroutine dumps to external file kfile the values of the spectral 
!   !absorption coefficient (or the pressure-based absorption 
!   !coefficient, if kappa_p = .true.) at specific thermodynamic states
!   !(temperature T and mole fractions xs) and within the spectral bounds
!   !leta < eta < ueta
!   !====================================================================
!   subroutine dump_kappa_eta(kfile,T,xs,leta,ueta,kappa_p)
      
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: assert,get_file_unit,print_to_prompt
!      use global_parameters, only: debug_mode
!      use precision_parameters, only: dp
!      implicit none
!      character(*) :: kfile
!      integer :: kunit,isp,its,nsp,nts
!      logical,intent(in),optional :: kappa_p
!      logical :: pbased
!      real(dp),intent(in) :: leta,T(:),ueta,xs(:,:)
!      real(dp) :: deta,eta_down,eta_factor,eta_up,xeta
!      real(dp) :: acs_down(lbl_ns,lbl_nx,lbl_nt),&
!                  acs_up(lbl_ns,lbl_nx,lbl_nt),kappa(size(T))

!      !-----------------------------------------------------------------
!      !Initial procedures
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('get_local_nbck: &
!                                           &initial procedures')
      
!      !Check for consistency in the array sizes
!      call assert(size(T).eq.size(xs,1))
      
!      !Surrogate names
!      nsp = size(xs,2)
!      nts = size(T)

!      !Set up optional parameter
!      pbased = .false.; if (present(kappa_p)) pbased = kappa_p
      
!      !-----------------------------------------------------------------
!      !Prepare the LBL data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('dump_kappa_eta: &
!                                           &prepare LBL data files')
!      call open_lbl_data
!      if (lbl_data_ready) call rewind_lbl_data
!      if (lbl_data_averaging.ne.'none') &
!         call real_lbl_data(eta_down,acs_down)                          !Read the first line of the data
!      eta_factor = 1._dp; if (lbl_data_cm) eta_factor = 100._dp         !Factor to convert from 1/cm to 1/m, if needed

!      !-----------------------------------------------------------------
!      !Prepare output file
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('dump_kappa_eta: &
!                                           &prepare output file')
!      kunit = get_file_unit()
!      open(file=kfile,unit=kunit,action='write',form='formatted')
!      write(kunit,'(a,1000(i5,:,","))') '#xeta [1/cm],',(its,its=1,nts)
                                           
!      !-----------------------------------------------------------------
!      !Main calculation
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('dump_kappa_eta: &
!                                           &Main calculation')
!      lbl_loop: do
         
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Read a line of the LBL databases
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         call real_lbl_data(eta_up,acs_up,deta)                         !Read a new line of the LBL data files

!         !Adjusting spectral bands for the analysis
!         if (lbl_data_averaging.ne.'none') &                            !Interval between bands
!            deta = (eta_up - eta_down)*eta_factor                                 
!         if ((lbl_data_averaging.eq.'arithmetic').or.&                  !Wavenumber position
!             (lbl_data_averaging.eq.'geometric')) &
!               xeta = 0.5_dp*(eta_up + eta_down)*eta_factor            
!         if ((lbl_data_averaging.eq.'upwind').or.&                      !Wavenumber position (no mean: use the
!             (lbl_data_averaging.eq.'none')) &                          !  current wavenumber value)
!               xeta = eta_up*eta_factor   
         
!         !Only consider lines that fall within the bounds of the
!         !narrow-band under consideration
!         if (xeta.lt.leta) then
!            eta_down = eta_up; acs_down = acs_up
!            cycle lbl_loop
!         elseif (xeta.gt.ueta) then
!            exit lbl_loop
!         endif
                  
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Define current absorption cross-section array
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         if (lbl_data_averaging.eq.'arithmetic') &
!            lbl_data = 0.5_dp*(acs_up + acs_down)
!         if (lbl_data_averaging.eq.'geometric') &
!            lbl_data = dsqrt(acs_up*acs_down)
!         if ((lbl_data_averaging.eq.'upwind').or.&
!             (lbl_data_averaging.eq.'none')) lbl_data = acs_up
         
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Compute the spectral absorption coefficient
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         state_loop: do its=1,nts
!            kappa(its) = 0._dp
!            do isp=1,nsp
!               kappa(its) = kappa(its) + &
!                  lbl_kappa_func(xs(its,isp),T(its),1._dp,isp,&
!                                 kappa_p=pbased)
!            enddo         
!         enddo state_loop
             
!         !Update spectral values
!         eta_down = eta_up                                              !Wavenumber
!         acs_down = acs_up                                              !Absorption cross-section

!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Dump the data
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         write(kunit,'(1000(e26.15e3,:,","))') &
!            xeta/100._dp,(kappa(its),its=1,nts)

!      enddo lbl_loop
      
!      !-----------------------------------------------------------------
!      !Close units
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('dump_kappa_eta: &
!                                           &close units')
!      close(kunit)
!      call close_lbl_data                                               !Close LBL units

!   endsubroutine dump_kappa_eta

   !====================================================================
   !Function to compute the band-averaged transmissivity of a gas
   !column with non-uniform temperature T and mole fraction xs, for a
   !band defined by the bounds leta < eta < ueta
   !====================================================================
   real(dp) function lbl_band_transmissivity(x,T,xs,leta,ueta)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,dprint,get_file_unit
      use precision_parameters, only: dp
      use omp_lib
      implicit none
      integer :: i,ilbl,isp,nsp,nx
      real(dp),intent(in) :: leta,T(:),ueta,x(:),xs(:,:)
      real(dp) :: dx,kappa,ot,sum_deta,sum_omp,tau,tau_omp
      
      !-----------------------------------------------------------------
      !Initial procedures
      !-----------------------------------------------------------------
      call dprint('lbl_band_transmissivity: Initial procedures')
      
      !Check for consistency in the array sizes
      call assert(size(x).eq.size(T))
      call assert(size(T).eq.size(xs,2))
      
      !Surrogate names
      nx = size(x)                                                      !Number of grid points
      nsp = size(xs,1)                                                  !Number of species
      
      !-----------------------------------------------------------------
      !Load the LBL data
      !-----------------------------------------------------------------
      call dprint('lbl_band_transmissivity: Read the LBL data')
      call load_lbl_data
      
      tau = 0._dp; sum_deta = 0._dp
      
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(dx,kappa,ot,tau_omp,sum_omp)
      
      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      tau_omp = 0._dp; sum_omp = 0._dp
      
      !$OMP DO
      lines_loop: do ilbl=1,lbl_nlines
         !--------------------------------------------------------------
         !Prepare the new spectral line
         !--------------------------------------------------------------
         !Check if this eta is to be solved for
         if ((lbl_xeta(ilbl).lt.leta).or.&
             (lbl_xeta(ilbl).gt.ueta)) cycle lines_loop

         !--------------------------------------------------------------
         !Compute the transmissivity
         !--------------------------------------------------------------
         ot = 0._dp
         grid_loop: do i=1,nx
            !Absorption coefficient
            kappa = 0._dp
            do isp=1,nsp
               kappa = kappa + &
                  lbl_kappa_func(xs(isp,i),T(i),1._dp,isp,ilbl)
            enddo
!write(*,*) lbl_xeta(ilbl),kappa
            !Cell size
            if (nx.eq.1) then
               dx = x(1)
            elseif (i.eq.1) then
               dx = 0.5_dp*(x(i+1) - x(i))
            elseif (i.eq.nx) then
               dx = 0.5_dp*(x(i) - x(i-1))
            else
               dx = 0.5_dp*(x(i+1) - x(i-1))
            endif
            
            ot = ot + kappa*dx
         enddo grid_loop
         tau_omp = tau_omp + dexp(-ot)*lbl_deta(ilbl)
         sum_omp = sum_omp + lbl_deta(ilbl)

      enddo lines_loop
      !$OMP ENDDO

      !-----------------------------------------------------------------
      !Finish the calculation
      !-----------------------------------------------------------------
      !$OMP CRITICAL
      tau = tau + tau_omp
      sum_deta = sum_deta + sum_omp
      !$OMP ENDCRITICAL
      !$OMP ENDPARALLEL

      lbl_band_transmissivity = tau/sum_deta

   endfunction lbl_band_transmissivity


   !====================================================================
   !Function to compute the emissivity of a gas column of width W, with
   !temperature T, total pressure p, and mole fraction xs
   !====================================================================
   real(dp) function lbl_emissivity(T,p,xs,w)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: dprint
      use physical_functions, only: Ib_function
      use precision_parameters, only: dp
      use omp_lib
      implicit none
      integer :: ilbl,isp,nsp
      real(dp),intent(in) :: p,T,xs(:),W
      real(dp) :: denum,denum_omp,expot,Ib,kappa,num,num_omp
      
      !-----------------------------------------------------------------
      !Initial procedures
      !-----------------------------------------------------------------
      call dprint('lbl_emissivity: Initial procedures')
      
      !Surrogate names
      nsp = size(xs)                                                    !Number of species
      
      !-----------------------------------------------------------------
      !Load the LBL data
      !-----------------------------------------------------------------
      call dprint('lbl_emissivity: Read the LBL data')
      call load_lbl_data
      
      num = 0._dp; denum = 0._dp
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(expot,Ib,kappa,num_omp,denum_omp)
      
      !-----------------------------------------------------------------
      !Initialize OMP sum arrays
      !-----------------------------------------------------------------
      num_omp = 0._dp; denum_omp = 0._dp
      
      !$OMP DO
      lines_loop: do ilbl=1,lbl_nlines
         !--------------------------------------------------------------
         !Prepare the new spectral line
         !--------------------------------------------------------------
         !Check if this eta is to be solved for
         if ((lbl_xeta(ilbl).lt.lbl_eta_min).or.&
             (lbl_xeta(ilbl).gt.lbl_eta_max)) cycle lines_loop

         !--------------------------------------------------------------
         !Compute the transmissivity
         !--------------------------------------------------------------
         !Absorption coefficient
         kappa = 0._dp
         do isp=1,nsp
            kappa = kappa + &
               lbl_kappa_func(xs(isp),T,p,isp,ilbl)
         enddo
         
         !Planck function
         Ib = Ib_function(T,lbl_xeta(ilbl))

         !Compute numerator and denominator
         expot = dexp(-kappa*W)
         num_omp = num_omp + Ib*(1._dp - expot)*lbl_deta(ilbl)
         denum_omp = denum_omp + Ib*lbl_deta(ilbl)
      enddo lines_loop
      !$OMP ENDDO

      !-----------------------------------------------------------------
      !Finish the calculation
      !-----------------------------------------------------------------
      !$OMP CRITICAL
      num = num + num_omp
      denum = denum + denum_omp
      !$OMP ENDCRITICAL
      !$OMP ENDPARALLEL

      lbl_emissivity = num/denum

   endfunction lbl_emissivity


   !====================================================================
   !Subroutine to uniformly reduce the LBL spectrum
   !====================================================================
   !Input:
   !  n_red: number of wavernumber positions that the base data will be
   !         averaged over for the analysis
   !  which_avg: integer indicating how the absorption cross-section for
   !             a reduced wavenumber position is to be determined
   !             (1: arithmetic mean, 2: log mean, 3: no average)
   subroutine uniformly_reduce_spectrum(n_red,which_avg,append_name)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit
      use global_parameters, only: number_of_species
      implicit none
      character(*),intent(in) :: append_name
      character(200) :: new_name,original_name
      integer,intent(in) :: n_red,which_avg
      integer :: ierr,ounit
      integer :: ilbl,iilbl,isp,ixs
      integer :: nlbl,nlbl_red,nsp,nxs
      real(dp) :: int_deta,int_eta
      real(dp),allocatable,dimension(:) :: deta_red,int_data,xeta_red
      real(dp),allocatable,dimension(:,:) :: data_red

      !-----------------------------------------------------------------
      !Load LBL data
      !-----------------------------------------------------------------
      call dprint('uniformly_reduce_spectrum: Load LBL data')
      call load_lbl_data
   
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('uniformly_reduce_spectrum: Allocate arrays')
      nlbl_red = ceiling(real(lbl_nlines,dp)/real(n_red,dp))
      allocate(data_red(nlbl_red,lbl_ntg),stat=ierr)
      call CheckMemAlloc('data_red',ierr)
      allocate(int_data(lbl_ntg),stat=ierr)
      call CheckMemAlloc('int_data',ierr)
      allocate(deta_red(nlbl_red),stat=ierr)
      call CheckMemAlloc('deta_red',ierr)
      allocate(xeta_red(nlbl_red),stat=ierr)
      call CheckMemAlloc('xeta_red',ierr)
      
      !-----------------------------------------------------------------
      !Main calculation loop
      !-----------------------------------------------------------------
      call dprint('uniformly_reduce_spectrum: Main calculation loop')
      nlbl = lbl_nlines
      nsp = number_of_species
      do isp=1,nsp
         nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle
         do ixs=1,nxs
            if (lbl_data_file(ixs,isp).eq.'null') cycle
            
            
            iilbl = 0
            int_data = 0._dp; int_eta = 0._dp; int_deta = 0._dp
            do ilbl=1,nlbl
               int_eta = int_eta + lbl_xeta(ilbl)*lbl_deta(ilbl)
               int_deta = int_deta + lbl_deta(ilbl)
               
               if (which_avg.eq.1) int_data = int_data + &
                  lbl_data(ilbl,:,ixs,isp)*lbl_deta(ilbl)
               if (which_avg.eq.2) int_data = &
                  int_data*(lbl_data(ilbl,:,ixs,isp)**lbl_deta(ilbl))
               if (which_avg.eq.3) int_data = lbl_data(ilbl,:,ixs,isp)
               
               if (mod(ilbl,n_red).eq.0) then
                  iilbl = iilbl + 1
                  if (which_avg.eq.1) data_red(iilbl,:) = &
                      int_data/int_deta
                  if (which_avg.eq.2) data_red(iilbl,:) = &
                     int_data**(1._dp/int_deta)
                  if (which_avg.eq.3) data_red(iilbl,:) = &
                     int_data
                  xeta_red(iilbl) = int_eta/int_deta
                  deta_red(iilbl) = int_deta
                  int_data = 0._dp; int_eta = 0._dp; int_deta = 0._dp
               endif
            enddo            
            
            !Define the file name
            original_name = trim(lbl_data_prefix)//&
                            trim(lbl_data_file(ixs,isp))
            new_name = trim(original_name)//'_'//trim(append_name)//&
                       trim(lbl_data_ext)
            
            !Open unit
            ounit = get_file_unit()
            open(file=new_name,unit=ounit,form='unformatted',&
                 action='write')
            
            !Dump data
            write(ounit) iilbl,.true.
            write(ounit) xeta_red
            write(ounit) deta_red
            write(ounit) data_red
write(*,*) nlbl,iilbl,nlbl_red
            !Close unit
            close(ounit)
         enddo
      enddo

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('uniformly_reduce_spectrum: Deallocate arrays')
      deallocate(data_red,deta_red,xeta_red,int_data)
      
   endsubroutine uniformly_reduce_spectrum

   !====================================================================
   !Subroutine to uniformly reduce the LBL spectrum
   !====================================================================
!   subroutine uniformly_reduce_spectrum(n_red,which_avg,red_nlines)
!   !Input:
!   !  n_red: number of wavernumber positions that the base data will be
!   !         averaged over for the analysis
!   !  which_avg: integer indicating how the absorption cross-section for
!   !             a reduced wavenumber position is to be determined
!   !             (1: arithmetic mean, 2: geometric mean, 
!   !              3: log mean, 4: no average)
!   !Output:
!   !  red_niles: optional integer with the total number of lines of the
!   !             reduced spectrum

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: assert,CheckMemAlloc,get_file_unit,&
!                                print_to_prompt
!      use global_parameters, only: debug_mode,number_of_species
!      use math_functions, only: locate
!      implicit none
!      character(8) :: sred
!      character(200) :: base_file,deta_file,msg,red_file
!      integer,intent(in) :: n_red
!      integer,intent(in),optional :: which_avg
!      integer,intent(out),optional :: red_nlines
!      integer :: ilbl,nlbl,nlbl_out,nsp,ntg,nxs
!      integer :: deta_unit,ierr
!      integer :: ict,iit,isp,itg,ixs
!      integer :: ind_avg
!      integer,allocatable,dimension(:) :: red_unit
!      real(dp) :: deta,deta_denum,eta_down,eta_up,xeta
!      real(dp) :: int_eta,eta_factor,eta_red
!      real(dp),allocatable,dimension(:,:,:) :: acs,acs_down,acs_red,&
!         acs_up,int_acs

!      !-----------------------------------------------------------------
!      !Preparatory procedures
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('uniformly_reduce_spectrum: &
!                                           &preparatory procedures')

!      !Optional input parameters
!      ind_avg = 4; if (present(which_avg)) ind_avg = which_avg

!      !Check for repeated flags that should be unique
!      call check_lbl_parameters

!      !Set up unit conversion factor
!      eta_factor = 100._dp; if (lbl_data_cm) eta_factor = 1._dp

!      !Set up string
!      write(sred,'(i8)') n_red

!      !-----------------------------------------------------------------
!      !Give a name, assign an unit value and open each output file
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('uniformly_reduce_spectrum: &
!                                           &name output files')
!      nsp = number_of_species                                           !Surrogate name for the total number of species
!      do iit=1,2                                                        !iit=1 -> Allocation; iit=2 -> Get and open unit
!         ict = 0                                                        !Initialize unit counter
!         do isp=1,nsp                                                   !Loop over all species
!            nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle                !Skip species for which there is no LBL data file
!            do ixs=1,nxs                                                !Loop over all discrete mole fraction values
!               if (lbl_data_file(isp,ixs).eq.'null') cycle              !Skip mole fraction values for wehich there is no LBL data file
!               ict = ict + 1                                            !Update unit counter
!               if (iit.eq.1) cycle                                      !Skip getting and opening the unit for allocation step
!               red_unit(ict) = get_file_unit()                          !Get file unit
!               base_file = lbl_data_file(isp,ixs)
!               base_file = base_file(1:len_trim(base_file)-4)
!               red_file = trim(base_file)//'_unired_'//&
!                          trim(adjustl(sred))//'.csv'               
!               open(unit=red_unit(ict),file=red_file,&                  !Open the unit
!                    form='formatted',action='write')
!            enddo
!         enddo
!         if (iit.eq.1) allocate(red_unit(ict))                          !Allocate array with all the units for reading the LBL data file
!      enddo

!      deta_file = 'lbldeta_unired_'//trim(adjustl(sred))//'.csv'
!      deta_unit = get_file_unit()
!      open(unit=deta_unit,file=deta_file,&
!           form='formatted',action='write')

!      !-----------------------------------------------------------------
!      !Allocate arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('uniformly_reduce_spectrum: &
!                                           &allocate arrays')
!      allocate(acs(lbl_nsp,lbl_nxs,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs',ierr)
!      allocate(acs_down(lbl_nsp,lbl_nxs,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_down',ierr)
!      allocate(acs_red(lbl_nsp,lbl_nxs,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_red',ierr)
!      allocate(acs_up(lbl_nsp,lbl_nxs,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_up',ierr)
!      allocate(int_acs(lbl_nsp,lbl_nxs,lbl_nt),stat=ierr)
!      call CheckMemAlloc('int_acs',ierr)

!      !-----------------------------------------------------------------
!      !Define, open units and read the first line of data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('uniformly_reduce_spectrum: &
!                                           &prepare LBL data files')
!      if (.not.lbl_data_ready.or.(lbl_nlines.le.0)) &                   !Either prepare the LBL data to be read
!         call prepare_lbl_data                                          !  (i.e., get units and open the external files)
!      if (lbl_data_ready.and.(lbl_nlines.gt.0)) call rewind_lbl_data    !  or rewind the already prepared LBL data
!      if (lbl_data_averaging.ne.'none') &
!         call real_lbl_data(eta_down,acs_down)                          !Read the first line of the data
!      nlbl = lbl_nlines-1                                               !Total number of lines in the data files
!      if (lbl_data_averaging.eq.'none') nlbl = nlbl + 1

!      !-----------------------------------------------------------------
!      !Main loop for reducing the spectral data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('uniformly_reduce_spectrum: &
!                                           &main loop')
!      nlbl_out = 0
!      int_acs = 0._dp; if (ind_avg.eq.2) int_acs = 1._dp
!      int_eta = 0._dp; deta_denum = 0._dp
!      main_loop: do ilbl=1,nlbl
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Prepare the new spectral line
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         call read_lbl_data(eta_up,acs_up)                              !Read a new line of the LBL data files

!         !Adjusting spectral bands for the analysis
!         deta = (eta_up - eta_down)                                     !Interval between bands
!         if (lbl_data_averaging.eq.'arithmetic') then                   !Arithmetic mean
!            xeta = 0.5_dp*(eta_up + eta_down)        
!            acs = 0.5_dp*(acs_up + acs_down)                  
!         elseif (lbl_data_averaging.eq.'geometric') then                !Geometric mean
!            xeta = 0.5_dp*(eta_up + eta_down)                           
!            acs = dsqrt(acs_up*acs_down)
!         elseif (lbl_data_averaging.eq.'none') then                     !No mean: use the current wavenumber value
!            xeta = eta_up
!            acs = acs_up
!         endif

!         !Check if this eta is to be solved for
!         if (xeta.lt.lbl_eta_min/100._dp) then
!            eta_down = eta_up
!            acs_down = acs_up
!            cycle main_loop
!         endif
!         if (xeta.gt.lbl_eta_max/100._dp) exit main_loop

!         !Print xeta if requested
!         if (print_xeta.or.debug_mode) then
!            write(msg,'(f12.4)') xeta/eta_factor
!            msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
!            call print_to_prompt(msg)
!         endif

!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Update spectral values
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         eta_down = eta_up                                              !Wavenumber
!         acs_down = acs_up                                              !Absorption cross-section

!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Apply the reduction to the absorption cross-section data
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         if (ind_avg.eq.1) int_acs = int_acs + acs*deta
!         if (ind_avg.eq.2) int_acs = int_acs*(acs**deta)
!         if (ind_avg.eq.4) int_acs = acs
!         int_eta = int_eta + xeta*deta
!         deta_denum = deta_denum + deta

!         if (mod(ilbl,n_red).eq.0) then
!            nlbl_out = nlbl_out + 1
!            if (ind_avg.eq.1) acs_red = int_acs/deta_denum
!            if (ind_avg.eq.2) acs_red = int_acs**(1._dp/deta_denum)
!            if (ind_avg.eq.4) acs_red = int_acs
!            eta_red = int_eta/deta_denum

!            !Write the data to file
!            ict = 0
!            do isp=1,nsp                                                !Loop over all species
!               nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle             !Skip species for which there is no LBL data file
!               do ixs=1,nxs                                             !Loop over all discrete mole fraction values
!                  if (lbl_data_file(isp,ixs).eq.'null') cycle           !Skip mole fraction values for wehich there is no LBL data file
!                  ict = ict + 1                                         !Update unit counter
!                  ntg = lbl_data_ntg(isp,ixs)                           !Surrogate name for number of discrete temperature values
!                  write(red_unit(ict),'(100(e26.15e3,:,","))') &        !Write one line of the LBL data file for species index isp and
!                     eta_red,(acs_red(isp,ixs,itg),itg=1,ntg)           !  mole fraction index ixs
!               enddo
!            enddo
!            write(deta_unit,'(e26.15e3)') deta_denum
!            int_acs = 0._dp; int_eta = 0._dp; deta_denum = 0._dp

!         endif
!      enddo main_loop

!      !Record total number of lines of the reduced spectrum
!      if (present(red_nlines)) red_nlines = nlbl_out

!      !-----------------------------------------------------------------
!      !Close the units
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('uniformly_reduce_spectrum: &
!                                           &close units')
!      do ict=1,size(red_unit)
!         close(red_unit(ict))
!      enddo
!      close(deta_unit)
!      call close_lbl_data

!      !-----------------------------------------------------------------
!      !Deallocate arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('uniformly_reduce_spectrum: &
!                                           &deallocate arrays')
!      deallocate(acs,acs_down,acs_red,acs_up,int_acs)

!   endsubroutine uniformly_reduce_spectrum

   !====================================================================
   !Subroutine to reduce the LBL spectrum according to the Spectrally
   !Reduced Integration (SRI) method of Coelho & França, 2021
   !====================================================================
   !Input:
   !  eta_in: array containing the wavenumber positions of the spectral
   !          data that will be the basis for the spectrum reduction
   !  deta_in: array containing the wavenumber interval of the spectral
   !           data that will be the basis for the spectrum reduction
   !  data_in: array containing the spectral data that will be the basis
   !           for the spectrum reduction
   !  n_red: number of wavernumber positions that the base data will be
   !         averaged over for the analysis
   !  zeta_array: array containing the downwind positions used for
   !              defining the different refinement levels
   !  red_array: array contaning the number of wavenumber positions to
   !             be skipped for each refinement level
   !  which_avg: integer indicating how the absorption cross-section for
   !             a reduced wavenumber position is to be determined
   !             (1: arithmetic mean, 2: geometric mean, 3: log mean)
   !Output:
   !  red_niles: optional integer with the total number of lines of the
   !             reduced spectrum
!   subroutine reduce_spectrum(eta_in,deta_in,data_in,n_red,&
!                              zeta_array,red_array,which_avg,&
!                              red_nlines)

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: assert,CheckMemAlloc,get_file_unit,&
!                                print_to_prompt
!      use global_parameters, only: debug_mode,number_of_species
!      use math_functions, only: locate
!      implicit none
!      character(200) :: base_file,deta_file,msg,red_file
!      integer,intent(in) :: red_array(:),n_red
!      integer,intent(in),optional :: which_avg
!      integer,intent(out),optional :: red_nlines
!      integer :: n,nlines,nlines_lbl,nlines_out,nsp,ntg,nxs,nzeta
!      integer :: deta_unit,ierr,red_counter
!      integer :: ict,idg,iit,isp,itg,ixs
!      integer :: ind_avg
!      integer,allocatable,dimension(:) :: red_factor,red_unit
!      real(dp),intent(in) :: data_in(:),deta_in(:),eta_in(:),&
!         zeta_array(:)
!      real(dp) :: deta,deta_denum,eta_down,eta_up,xeta
!      real(dp) :: int_data,int_eta,eta_red
!      real(dp) :: data_max,data_min
!      real(dp),allocatable,dimension(:,:,:) :: acs,acs_down,acs_red,&
!         acs_up,int_acs
      
!      !-----------------------------------------------------------------
!      !Preparatory procedures
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('reduce_spectrum: &
!                                           &preparatory procedures')
      
!      !Check if the input arrays have the correct sizes
!      call assert(size(eta_in).eq.size(deta_in))
!      call assert(size(deta_in).eq.size(data_in))
!      call assert(size(zeta_array).eq.size(red_array))
      
!      !Array sizes
!      nlines = size(deta_in)
!      nzeta = size(zeta_array)
      
!      !Optional input parameters
!      ind_avg = 1; if (present(which_avg)) ind_avg = which_avg
      
!      !Check for repeated flags that should be unique
!      call check_lbl_parameters

!      !Maximum and minimum values for the data
!      data_max = maxval(data_in)
!      data_min = minval(data_in)

!      !-----------------------------------------------------------------
!      !Give a name, assign an unit value and open each output file
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('reduce_spectrum: &
!                                           &name output files')
!      nsp = number_of_species                                           !Surrogate name for the total number of species
!      do iit=1,2                                                        !iit=1 -> Allocation; iit=2 -> Get and open unit
!         ict = 0                                                        !Initialize unit counter
!         do isp=1,nsp                                                   !Loop over all species
!            nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle                !Skip species for which there is no LBL data file
!            do ixs=1,nxs                                                !Loop over all discrete mole fraction values
!               if (lbl_data_file(isp,ixs).eq.'null') cycle              !Skip mole fraction values for wehich there is no LBL data file
!               ict = ict + 1                                            !Update unit counter
!               if (iit.eq.1) cycle                                      !Skip getting and opening the unit for allocation step
!               red_unit(ict) = get_file_unit()                          !Get file unit
!               base_file = lbl_data_file(isp,ixs)
!               base_file = base_file(1:len_trim(base_file)-4)
!               red_file = trim(base_file)//'_red.csv'               
!               open(unit=red_unit(ict),file=red_file,&                  !Open the unit
!                    form='formatted',action='write')
!            enddo
!         enddo
!         if (iit.eq.1) allocate(red_unit(ict))                          !Allocate array with all the units for reading the LBL data file
!      enddo
      
!      deta_file = 'lbldeta_red.csv'
!      deta_unit = get_file_unit()
!      open(unit=deta_unit,file=deta_file,&
!           form='formatted',action='write')

!      !-----------------------------------------------------------------
!      !Allocate arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('reduce_spectrum: &
!                                           &allocate arrays')
!      allocate(acs(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs',ierr)
!      allocate(acs_down(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_down',ierr)
!      allocate(acs_red(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_red',ierr)
!      allocate(acs_up(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_up',ierr)
!      allocate(int_acs(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('int_acs',ierr)
!      allocate(red_factor(nlines),stat=ierr)
!      call CheckMemAlloc('red_factor',ierr)
!      red_factor = -1

!      !-----------------------------------------------------------------
!      !Reduce the input spectral data by integrating across the 
!      !spectral interval deta_red, following Eqs. (18) and (19) of
!      !Coelho & França, 2021
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('reduce_spectrum: &
!                                           &reduce input spectral data')
!      red_counter = 0
!      int_data = 0._dp; int_eta = 0._dp
!      reduction_loop: do n=1,nlines
!         int_data = int_data + &
!            (data_in(n) - data_min)/(data_max - data_min)*deta_in(n)
!         int_eta = int_eta + eta_in(n)*deta_in(n)
!         deta_denum = deta_denum + deta_in(n)
!         if (mod(n,n_red).eq.0) then
!            red_counter = red_counter + 1
!            int_data = int_data/deta_denum
!            idg = min(max(locate(zeta_array,int_data,nzeta)+1,1),nzeta)
!            red_factor(red_counter) = red_array(idg)
!            int_data = 0._dp; deta_denum = 0._dp; int_eta = 0._dp
!         endif
!      enddo reduction_loop

!      !-----------------------------------------------------------------
!      !Define, open units and read the first line of data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('reduce_spectrum: &
!                                           &prepare LBL data files')
!      if (.not.lbl_data_ready.or.(lbl_nlines.le.0)) &                   !Either prepare the LBL data to be read
!         call prepare_lbl_data                                          !  (i.e., get units and open the external files)
!      if (lbl_data_ready.and.(lbl_nlines.gt.0)) call rewind_lbl_data    !  or rewind the already prepared LBL data
!      if (lbl_data_averaging.ne.'none') &
!         call real_lbl_data(eta_down,acs_down)                          !Read the first line of the data
!      nlines_lbl = lbl_nlines-1                                         !Total number of lines in the data files
!      if (lbl_data_averaging.eq.'none') nlines_lbl = nlines_lbl + 1

!      !-----------------------------------------------------------------
!      !Main loop for reducing the spectral data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('reduce_spectrum: &
!                                           &main loop')
!      red_counter = 1; nlines_out = 0
!      int_acs = 0._dp; if (ind_avg.eq.2) int_acs = 1._dp
!      int_eta = 0._dp; deta_denum = 0._dp
!      main_loop: do n=1,min(nlines,nlines_lbl)
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Prepare the new spectral line
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         call real_lbl_data(eta_up,acs_up)                              !Read a new line of the LBL data files

!         !Adjusting spectral bands for the analysis
!         deta = (eta_up - eta_down)                                     !Interval between bands
!         if (lbl_data_averaging.eq.'arithmetic') then                   !Arithmetic mean
!            xeta = 0.5_dp*(eta_up + eta_down)        
!            acs = 0.5_dp*(acs_up + acs_down)                  
!         elseif (lbl_data_averaging.eq.'geometric') then                !Geometric mean
!            xeta = dsqrt(eta_up*eta_down)                           
!            acs = dsqrt(acs_up*acs_down)
!         elseif (lbl_data_averaging.eq.'none') then                     !No mean: use the current wavenumber value
!            xeta = eta_up
!            acs = acs_up
!         endif

!         !Check is the spectral positions match
!         call assert(dabs(xeta-eta_in(n)).lt.1e-5_dp)

!         !Converting from 1/cm to 1/m
!         if (lbl_data_cm) deta = deta*100._dp
!         if (lbl_data_cm) xeta = xeta*100._dp

!         !Check if this eta is to be solved for
!         if (xeta.lt.lbl_eta_min) then
!            eta_down = eta_up
!            acs_down = acs_up
!            cycle main_loop
!         endif
!         if (xeta.gt.lbl_eta_max) exit main_loop

!         !Print xeta if requested
!         if (print_xeta.or.debug_mode) then
!            write(msg,'(f12.4)') xeta/100._dp
!            msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
!            call print_to_prompt(msg)
!         endif

!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Update spectral values
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         eta_down = eta_up                                              !Wavenumber
!         acs_down = acs_up                                              !Absorption cross-section

!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Select appropriate correction factor
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         if (mod(n,n_red).eq.0) red_counter = red_counter + 1
!         if (red_factor(red_counter).lt.0) exit main_loop

!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         !Apply the reduction to the absorption cross-section data
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         if (ind_avg.eq.1) int_acs = int_acs + acs*deta
!         if (ind_avg.eq.2) int_acs = int_acs*(acs**deta)
!         int_eta = int_eta + xeta*deta
!         deta_denum = deta_denum + deta

!         if (mod(n,red_factor(red_counter)).eq.0) then
!            nlines_out = nlines_out + 1
!            if (ind_avg.eq.1) acs_red = int_acs/deta_denum
!            if (ind_avg.eq.2) acs_red = int_acs**(1._dp/deta_denum)
!            eta_red = int_eta/deta_denum
         
!            !Write the data to file
!            ict = 0
!            do isp=1,nsp                                                !Loop over all species
!               nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle             !Skip species for which there is no LBL data file
!               do ixs=1,nxs                                             !Loop over all discrete mole fraction values
!                  if (lbl_data_file(isp,ixs).eq.'null') cycle           !Skip mole fraction values for wehich there is no LBL data file
!                  ict = ict + 1                                         !Update unit counter
!                  ntg = lbl_data_ntg(isp,ixs)                           !Surrogate name for number of discrete temperature values
!                  write(red_unit(ict),'(100(e26.15e3,:,","))') &        !Write one line of the LBL data file for species index isp and
!                     eta_red,(acs_red(isp,ixs,itg),itg=1,ntg)           !  mole fraction index ixs
!               enddo
!            enddo
!            write(deta_unit,'(e26.15e3)') deta_denum
            
!            int_acs = 0._dp; int_eta = 0._dp; deta_denum = 0._dp
            
!         endif
!      enddo main_loop
   
!      !Record total number of lines of the reduced spectrum
!      if (present(red_nlines)) red_nlines = nlines_out
      
!      !-----------------------------------------------------------------
!      !Close the units
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('reduce_spectrum: &
!                                           &close units')
!      do ict=1,size(red_unit)
!         close(red_unit(ict))
!      enddo
!      close(deta_unit)
      
!      !-----------------------------------------------------------------
!      !Deallocate arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('reduce_spectrum: &
!                                           &deallocate arrays')
!      deallocate(acs,acs_down,acs_red,acs_up,int_acs)
!      deallocate(red_factor)
      
!   endsubroutine reduce_spectrum

   !====================================================================
   !Subroutine to generate the random number relation for a single 
   !species, to be used in the LBL/Monte Carlo method
   !====================================================================
   subroutine get_lbl_rnr(T_in,p_in,xs_in,ofile,binary_output,&
      spec_index,pressure_based)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,get_file_unit,&
                                print_to_prompt
      use global_parameters, only: id_soot,number_of_species
      use physical_functions, only: Ib_function,kappa_eta_soot
      use precision_parameters, only: dp,small
      implicit none
      character(*),intent(in) :: ofile
      character(200) :: msg
      integer,intent(in),optional :: spec_index
      integer :: ierr,ounit
      integer :: nlbl,nspc
      integer :: ilbl,ispc
      real(dp) :: Ib,Ib_total,kappa,kappa_spc,kPi
      real(dp),intent(in) :: p_in,T_in,xs_in(:)
      real(dp),allocatable,dimension(:) :: ki,Rk
      logical,intent(in),optional :: binary_output,pressure_based
      logical :: binary,single_species,pbased

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('get_lbl_rnr: preparatory procedures')
      
      !Set optional flags
      binary = .false.
      if (present(binary_output)) binary = binary_output
      
      single_species = .false.
      if (present(spec_index)) single_species = .true.
      if (present(spec_index)) ispc = spec_index
      
      pbased = .false.
      if (present(pressure_based)) pbased = pressure_based
      
      !-----------------------------------------------------------------
      !Load the LBL data
      !-----------------------------------------------------------------
      call dprint('get_lbl_rnr: Read the LBL data')
      call load_lbl_data
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('get_lbl_rnr: Surrogate names')
      nlbl = lbl_nlines                                                 !Number of lines in the LBL data files
      nspc = number_of_species                                          !Number of participating species
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      allocate(ki(nlbl),stat=ierr)
      call CheckMemAlloc('ki',ierr)
      allocate(Rk(nlbl),stat=ierr)
      call CheckMemAlloc('Rk',ierr)
      
      !-----------------------------------------------------------------
      !Zero-out quantities
      !-----------------------------------------------------------------
      call dprint('get_lbl_rnr: zero-out quantities')
      Rk = 0._dp; Ib_total = 0._dp
      
      lines_loop: do ilbl=1,nlbl                                        !Main loop
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
         !Compute properties
         !--------------------------------------------------------------
         !Absorption coefficient
         kappa = 0._dp
         if (single_species) then
            if (ispc.eq.id_soot) then
               kappa = kappa_eta_soot(lbl_xeta(ilbl),xs_in(ispc))
            else                  
               kappa = lbl_kappa_func(xs_in(ispc),T_in,p_in,&
                                          ispc,ilbl,kappa_p=pbased)
            endif
         else
            do ispc=1,nspc
               if (ispc.eq.id_soot) then
                  kappa_spc = kappa_eta_soot(lbl_xeta(ilbl),xs_in(ispc))
               else                  
                  kappa_spc = lbl_kappa_func(xs_in(ispc),T_in,p_in,ispc,ilbl)
               endif
               kappa = kappa + kappa_spc
            enddo
            if (pbased) kappa = kappa/(sum(xs_in)*p_in)
         endif

         !Spectral blackbory radiative intensity
         Ib = Ib_function(T_in,lbl_xeta(ilbl))

         !--------------------------------------------------------------
         !Store quantities for later
         !--------------------------------------------------------------
         ki(ilbl) = kappa                                               !Store the current pressure-based spectral 
                                                                        !  absorption coefficient into array ki 
                                                                        !  (which will be writted to file later)
         Rk(ilbl) = Rk(max(ilbl-1,1)) + kappa*Ib*lbl_deta(ilbl)         !Compute the cumulative emitted energy
         Ib_total = Ib_total + Ib*lbl_deta(ilbl)                        !Integrate for the total blackbody function

      enddo lines_loop
      
      !-----------------------------------------------------------------
      !Final calculations
      !-----------------------------------------------------------------
      call dprint('get_lbl_rnr: final calculations',3)
      
      !Compute Planck-mean absorption coefficient
      kPi = Rk(nlbl)/Ib_total
      
      !Finish computing the random number
      do ilbl=1,nlbl
         Rk(ilbl) = Rk(ilbl)/(Rk(nlbl) + small)
      enddo
      
      !-----------------------------------------------------------------
      !Dump data
      !-----------------------------------------------------------------
      call dprint('get_lbl_rnr: manage output file',3)
      if (binary) then
         !Prepare output unit
         ounit = get_file_unit()                                        !Get unit
         open(unit=ounit,file=ofile,form='unformatted',action='write')  !Open unit
         
         !Dump initial lines
         write(ounit) nlbl
         write(ounit) T_in,p_in,xs_in
         
         !Dump data
         write(ounit) lbl_xeta
         write(ounit) Rk
         write(ounit) ki
         write(ounit) kPi
         
         !Close unit
         close(ounit)
      
      else
         !Prepare output unit
         ounit = get_file_unit()                                        !Get unit
         open(unit=ounit,file=ofile,form='formatted',action='write')    !Open unit
         
         !Write header
         write(ounit,'(100(a,:,","))') '#eta','Rk_eta','kappa_eta'
         
         !Dump data
         do ilbl=1,nlbl
            write(ounit,'(100(e22.15e3,:,","))') lbl_xeta(ilbl),&
               Rk(ilbl),ki(ilbl)
         enddo
         
         !Close unit
         close(ounit)
      
      endif
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('get_lbl_rnr: deallocate arrays')
      deallocate(ki,Rk)

   endsubroutine get_lbl_rnr

!   !====================================================================
!   !Subroutine to generate the random number relation for a single 
!   !species, to be used in the LBL/Monte Carlo method
!   !====================================================================
!   subroutine precompute_lblrnr(file_in,T_in,xs_in)

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: CheckMemAlloc,get_file_unit,&
!                                print_to_prompt
!      use constants, only: invpi,sigma
!      use global_parameters, only: debug_mode,number_of_species
!      use physical_functions, only: Ib_function
!      use precision_parameters, only: dp,small
!      implicit none
!      character(*),intent(in) :: file_in(:)
!      character(200) :: msg,ofile
!      integer :: dpos,ierr,ounit
!      integer :: nfl,nln,ntg
!      integer :: ifl,iln,itg
!      integer :: number_of_species_mem,lbl_data_nxs_mem(lbl_ns)
!      real(dp),parameter :: sigrpi = sigma*invpi
!      real(dp) :: deta,eta_down,eta_factor,eta_up,xeta
!      real(dp) :: Ib,kappa
!      real(dp) :: lbl_data_xs_mem(lbl_nx,lbl_ns)
!      real(dp),intent(in) :: T_in(:),xs_in(:)
!      real(dp),allocatable,dimension(:) :: eta,kPi
!      real(dp),allocatable,dimension(:,:) :: ki,Rk
!      real(dp),allocatable,dimension(:,:,:) :: acs_down,acs_up

!      !-----------------------------------------------------------------
!      !Preparatory procedures
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                           &preparatory procedures')
!      !Surrogate names
!      ntg = size(T_in)                                                  !Number of discrete temperature values
!      nfl = size(file_in)                                               !Number of input files

!      !Set up factor to convert from 1/cm to 1/m, if needed
!      eta_factor = 1._dp; if (lbl_data_cm) eta_factor = 100._dp

!      !-----------------------------------------------------------------
!      !Redefine parameters of the LBL database so that 
!      !only a single datafile is expected
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                     &redefine LBL database parameters')
      
!      !Record original values of the LBL database parameters
!      number_of_species_mem = number_of_species
!      lbl_data_nxs_mem = lbl_data_nxs
!      lbl_data_xs_mem = lbl_data_xs
      
!      !Update values of the LBL database parameters
!      number_of_species = 1                                             !Each input file will be processed separately, being
!      lbl_data_nxs(1) = 1                                               !  treated as if it was the only one available for the
!                                                                        !  calculation of the spectral absorption coefficient

!      !-----------------------------------------------------------------
!      !Allocate arrays (pt I)
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                           &allocate arrays (pt I)')
!      allocate(acs_down(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_down',ierr)
!      allocate(acs_up(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_up',ierr)

!      file_loop: do ifl=1,nfl                                           !Loop over the LBL data files
!         !--------------------------------------------------------------
!         !Prepare the current LBL database to be read
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                   &prepare the current LBL database',3)
!         lbl_data_file(1,1) = file_in(ifl)                              !Database name
!         call check_lbl_parameters                                      !Check for repeated flags that should be unique
!         call prepare_lbl_data                                          !Prepare the data
!         if (lbl_data_averaging.ne.'none') &
!            call real_lbl_data(eta_down,acs_down)                       !Read the first line of the data
!         nln = lbl_nlines-1                                             !Total number of lines in the data files
!         if (lbl_data_averaging.eq.'none') nln = nln + 1
!         lbl_data_xs(1,1) = xs_in(ifl)                                  !Mole fraction of the current LBL data file
         
!         !--------------------------------------------------------------
!         !Allocate arrays (pt II)
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                             &allocate arrays (pt III)')
!         allocate(eta(nln),stat=ierr)
!         call CheckMemAlloc('eta',ierr)
!         allocate(ki(nln,ntg),stat=ierr)
!         call CheckMemAlloc('ki',ierr)
!         allocate(kPi(ntg),stat=ierr)
!         call CheckMemAlloc('kPi',ierr)
!         allocate(Rk(nln,ntg),stat=ierr)
!         call CheckMemAlloc('Rk',ierr)

!         !--------------------------------------------------------------
!         !Manage output file
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                              &manage output file',3)
!         ounit = get_file_unit()                                        !Get unit
!         ofile = trim(file_in(ifl))                                     !Set up the file name, which is defined by appending
!         dpos = index(trim(ofile),'.',.true.) - 1                       !  '_rnr' to the end of the name of the original LBL
!         ofile = ofile(1:dpos)//'_rnr.data'                             !  data file
!         open(unit=ounit,file=ofile,form='unformatted',action='write')  !Open unit
!         write(ounit) ntg, nln; write(ounit) T_in                       !Dump initial lines
      
!         !--------------------------------------------------------------
!         !Main spectral loop
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                              &spectral loop',3)
!         ki = 0._dp; Rk = 0._dp
!         spectral_loop: do iln=1,nln
!            !Read a new line of the LBL data files
!            call real_lbl_data(eta_up,acs_up,deta)

!            !Adjusting spectral bands for the analysis
!            if (lbl_data_averaging.ne.'none') &                         !Interval between bands
!               deta = (eta_up - eta_down)                                 
!            if ((lbl_data_averaging.eq.'arithmetic').or.&               !Wavenumber position
!               (lbl_data_averaging.eq.'geometric')) &
!                  xeta = 0.5_dp*(eta_up + eta_down)              
!            if ((lbl_data_averaging.eq.'upwind').or.&                   !Wavenumber position (no mean: use the
!               (lbl_data_averaging.eq.'none')) xeta = eta_up            !  current wavenumber value)

!            !Converting from 1/cm to 1/m
!            deta = deta*eta_factor                                      !Wavenumber spacing
!            xeta = xeta*eta_factor                                      !Wavenumber position
!            eta(iln) = xeta                                             !Store the current wavenumber position into array eta
!                                                                        !  (which will be writted to file later)
!            !Print xeta if requested
!            if (print_xeta.or.debug_mode) then
!               write(msg,'(f12.4)') xeta/100._dp
!               msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
!               call print_to_prompt(msg,6)
!            endif

!            !Current absorption cross-section array
!            if (lbl_data_averaging.eq.'arithmetic') &
!               lbl_data = 0.5_dp*(acs_up + acs_down)
!            if (lbl_data_averaging.eq.'geometric') &
!               lbl_data = dsqrt(acs_up*acs_down)
!            if ((lbl_data_averaging.eq.'upwind').or.&
!               (lbl_data_averaging.eq.'none')) lbl_data = acs_up

!            !Add to the integrated quantities
!            do itg=1,ntg
!               Ib = Ib_function(T_in(itg),xeta)                         !Spectral blackbory radiative intensity
!               kappa = lbl_kappa_func(xs_in(ifl),T_in(itg),1._dp,1,.true.)!Pressure-based spectral absorption coefficient
!               ki(iln,itg) = kappa                                      !Store the current pressure-based spectral 
!                                                                        !  absorption coefficient into array ki 
!                                                                        !  (which will be writted to file later)
!               Rk(iln,itg) = Rk(max(iln-1,1),itg) + kappa*Ib*deta       !Compute the cumulative emitted energy
!            enddo

!            !Update spectral values
!            eta_down = eta_up
!            acs_down = acs_up

!         enddo spectral_loop

!         !--------------------------------------------------------------
!         !Compute the pressure-based Planck-mean  
!         !absorption coefficient for the species
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                     &compute the Planck-mean absorption coefficient',3)
!         do itg=1,ntg
!            Ib = sigrpi*T_in(itg)**4._dp                                !Total blackbody intensity
!            kPi(itg) = Rk(nln,itg)/Ib                                   !Pressure-based Planck-mean absorption coefficient
!         enddo

!         !--------------------------------------------------------------
!         !Finish the calculation of the random number relations
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!              &finish the calculation of the random number relations',3)
!         do iln=1,nln
!            do itg=1,ntg
!               Rk(iln,itg) = Rk(iln,itg)/(Rk(nln,itg) + small)
!            enddo
!         enddo

!         !--------------------------------------------------------------
!         !Dump the data
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                              &dump the data',3)
!         write(ounit) eta
!         write(ounit) Rk
!         write(ounit) ki
!         write(ounit) kPi

!         !--------------------------------------------------------------
!         !Prepare for the next step of the file loop
!         !--------------------------------------------------------------
!         if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                         &prepare for the next step of the file loop',3)
!         close(ounit)                                                   !Close current output unit
!         call close_lbl_data                                            !Close current LBL data file
!         deallocate(eta,ki,kPi,Rk)                                      !Deallocate arrays that will be allocated again

!      enddo file_loop

!      !-----------------------------------------------------------------
!      !Dealocate arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                                           &deallocate arrays')
!      deallocate(acs_down,acs_up)

!      !-----------------------------------------------------------------
!      !Restore the LBL database parameters to their original values
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('precompute_lblrnr: &
!                             &restore original LBL database parameters')
!      number_of_species = number_of_species_mem
!      lbl_data_nxs = lbl_data_nxs_mem
!      lbl_data_xs = lbl_data_xs_mem

!   endsubroutine precompute_lblrnr

   !====================================================================
   !Subroutine to generate an artificial absorption spectrum
   !(for testing purposes)
   !====================================================================
   subroutine generate_artificial_spectrum(file_name,nlines,ntmp,&
                                           eta_min,eta_max,scale_factor)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: get_file_unit,get_wall_time
      use math_functions, only: get_numrep_random
      implicit none
      character(*),intent(in) :: file_name
      integer,intent(in) :: nlines,ntmp
      integer :: idum,iln,itg,nln,ntg,ounit
      real(dp),intent(in),optional :: eta_min,eta_max,scale_factor
      real(dp) :: deta,emin,emax,eta,sf,idum_dp
      real(dp),allocatable,dimension(:) :: acs

      !-----------------------------------------------------------------
      !Set up optional parameters
      !-----------------------------------------------------------------
      emin = 0._dp; if (present(eta_min)) emin = eta_min
      emax = 1._dp; if (present(eta_max)) emax = eta_max
      sf = 1.e-20_dp; if (present(scale_factor)) sf = scale_factor

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      nln = nlines; ntg = ntmp

      !-----------------------------------------------------------------
      !Set up output
      !-----------------------------------------------------------------
      ounit = get_file_unit()
      open(file=file_name,unit=ounit,form='formatted',action='write')

      !-----------------------------------------------------------------
      !Set a random seed
      !-----------------------------------------------------------------
      idum_dp = get_wall_time()
      idum_dp = (idum_dp - (floor(idum_dp)))*999._dp
      idum = -nint(idum_dp)

      !-----------------------------------------------------------------
      !Main loop
      !-----------------------------------------------------------------
      allocate(acs(ntg))
      eta = emin
      deta = (emax - emin)/real(nlines-1,dp)
      do iln=1,nln
!         if (iln.eq.1) then
!            acs = (/ 0.25_dp, 1._dp /)*sf
!         else
!            acs = (/ 2._dp, 0.25_dp /)*sf
!         endif
!         if (iln.eq.1) eta = emin
!         if (iln.eq.2) eta = 0.2_dp*(emin + emax)
!         if (iln.eq.3) eta = emax
         
         do itg=1,ntg
            acs(itg) = get_numrep_random(idum)*sf
         enddo
         write(ounit,'(100(e26.15,:,","))') eta,(acs(itg),itg=1,ntg)
         eta = eta + deta
      enddo
      deallocate(acs)
      close(ounit)

   endsubroutine generate_artificial_spectrum

!   !====================================================================
!   !Subroutine to generate the random number relation for a single 
!   !species, to be used in the LBL/Monte Carlo method
!   !====================================================================
!   subroutine LBLMC_precompute(T_in,xs_in,output_file)

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: CheckMemAlloc,get_file_unit,&
!                                print_to_prompt
!      use constants, only: invpi,sigma
!      use global_parameters, only: debug_mode
!      use physical_functions, only: Ib_function
!      use precision_parameters, only: dp,small
!      implicit none
!      character(*),intent(in) :: output_file
!      character(200) :: msg
!      integer :: ierr,ounit
!      integer :: nln,nsp,ntg,nxp
!      integer :: iln,isp,itg,ixp
!      real(dp),parameter :: sigrpi = sigma*invpi
!      real(dp) :: deta,eta_down,eta_factor,eta_up,xeta
!      real(dp) :: Ib,kappa
!      real(dp),intent(in) :: T_in(:),xs_in(:,:)
!      real(dp),allocatable,dimension(:) :: eta
!      real(dp),allocatable,dimension(:,:,:) :: acs_down,acs_up,kPi
!      real(dp),allocatable,dimension(:,:,:,:) :: ki,Rk

!      !-----------------------------------------------------------------
!      !Preparatory procedures
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                                           &preparatory procedures')
!      !Surrogate names
!      ntg = size(T_in)
!      nxp = size(xs_in,1)
!      nsp = size(xs_in,2)

!      !Factor to convert from 1/cm to 1/m, if needed
!      eta_factor = 1._dp; if (lbl_data_cm) eta_factor = 100._dp
      
!      !-----------------------------------------------------------------
!      !Allocate arrays (pt I)
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                                           &allocate arrays (pt I)')
!      allocate(acs_down(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_down',ierr)
!      allocate(acs_up(lbl_ns,lbl_nx,lbl_nt),stat=ierr)
!      call CheckMemAlloc('acs_up',ierr)
      
!      !-----------------------------------------------------------------
!      !Define, open units and read the first line of data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                                           &prepare LBL data files')
!      call check_lbl_parameters                                         !Check for repeated flags that should be unique
!      if (.not.lbl_data_ready.or.(lbl_nlines.le.0)) &                   !Either prepare the LBL data to be read
!         call prepare_lbl_data                                          !  (i.e., get units and open the external files)
!      if (lbl_data_ready.and.(lbl_nlines.gt.0)) call rewind_lbl_data    !  or rewind the already prepared LBL data
!      if (lbl_data_averaging.ne.'none') &
!         call real_lbl_data(eta_down,acs_down)                          !Read the first line of the data
!      nln = lbl_nlines-1                                                !Total number of lines in the data files
!      if (lbl_data_averaging.eq.'none') nln = nln + 1

!      !-----------------------------------------------------------------
!      !Allocate arrays (pt II)
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                                           &allocate arrays (pt II)')
!      allocate(eta(nln),stat=ierr)
!      call CheckMemAlloc('eta',ierr)
!      allocate(ki(nsp,nxp,ntg,nln),stat=ierr)
!      call CheckMemAlloc('ki',ierr)
!      allocate(kPi(nsp,nxp,ntg),stat=ierr)
!      call CheckMemAlloc('kPi',ierr)
!      allocate(Rk(nsp,nxp,ntg,nln),stat=ierr)
!      call CheckMemAlloc('Rk',ierr)

!      !-----------------------------------------------------------------
!      !Manage output file
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                                           &manage output file')
!      ounit = get_file_unit()                                           !Get unit
!      open(unit=ounit,file=output_file,&                                !Open unit
!         form='unformatted',action='write')
!      write(ounit) ntg,nsp,nxp,nln; write(ounit) T_in,xs_in             !Dump initial lines
      
!      !-----------------------------------------------------------------
!      !Main calculation loop
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                                           &main loop')
!      ki = 0._dp; Rk = 0._dp
!      main_loop: do iln=1,nln
         
!         !Read a new line of the LBL data files
!         call real_lbl_data(eta_up,acs_up,deta)

!         !Adjusting spectral bands for the analysis
!         if (lbl_data_averaging.ne.'none') &                            !Interval between bands
!            deta = (eta_up - eta_down)                                 
!         if ((lbl_data_averaging.eq.'arithmetic').or.&                  !Wavenumber position
!             (lbl_data_averaging.eq.'geometric')) &
!               xeta = 0.5_dp*(eta_up + eta_down)              
!         if ((lbl_data_averaging.eq.'upwind').or.&                      !Wavenumber position (no mean: use the
!             (lbl_data_averaging.eq.'none')) xeta = eta_up              !  current wavenumber value)

!         !Converting from 1/cm to 1/m
!         deta = deta*eta_factor
!         xeta = xeta*eta_factor
!         eta(iln) = xeta
         
!         !Print xeta if requested
!         if (print_xeta.or.debug_mode) then
!            write(msg,'(f12.4)') xeta/100._dp
!            msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
!            call print_to_prompt(msg,3)
!         endif
                                                     
!         !Current absorption cross-section array
!         if (lbl_data_averaging.eq.'arithmetic') &
!            lbl_data = 0.5_dp*(acs_up + acs_down)
!         if (lbl_data_averaging.eq.'geometric') &
!            lbl_data = dsqrt(acs_up*acs_down)
!         if ((lbl_data_averaging.eq.'upwind').or.&
!             (lbl_data_averaging.eq.'none')) lbl_data = acs_up

!         !Add to the integrated quantities
!         do itg=1,ntg
!            Ib = Ib_function(T_in(itg),xeta)                            !Spectral blackbory radiative intensity
!            do isp=1,nsp
!               do ixp=1,nxp
!                  kappa = lbl_kappa_func(xs_in(ixp,isp),T_in(itg),&
!                                         1._dp,isp)
!                  ki(isp,ixp,itg,iln) = kappa
!                  Rk(isp,ixp,itg,iln) = Rk(isp,ixp,itg,max(iln-1,1)) + &
!                     kappa*Ib*deta
!               enddo
!            enddo
!         enddo
!      enddo main_loop

!      !-----------------------------------------------------------------
!      !Compute the Planck-mean absorption coefficient for each species
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                     &compute the Planck-mean absorption coefficient')
!      do itg=1,ntg
!         Ib = sigrpi*T_in(itg)**4._dp
!         do isp=1,nsp
!            do ixp=1,nxp
!               kPi(isp,ixp,itg) = Rk(isp,ixp,itg,nln)/Ib
!            enddo
!         enddo
!      enddo

!      !-----------------------------------------------------------------
!      !Finish the calculation of the random number relations
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!               &finish the calculation of the random number relations')
!      do iln=1,nln
!         do itg=1,ntg
!            do isp=1,nsp
!               do ixp=1,nxp
!                  Rk(isp,ixp,itg,iln) = Rk(isp,ixp,itg,iln)/&
!                     (Rk(isp,ixp,itg,nln) + small)
!               enddo
!            enddo
!         enddo
!      enddo
      
!      !-----------------------------------------------------------------
!      !Dump the data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                                           &dump the data')
!      write(ounit) eta
!      write(ounit) Rk
!      write(ounit) ki
!      write(ounit) kPi
!      close(ounit)
      
!      !-----------------------------------------------------------------
!      !Dealocate arrays
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('LBLMC_precompute: &
!                                           &deallocate arrays')
!      deallocate(acs_down,acs_up)
!      deallocate(eta,ki,kPi,Rk)
      
!   endsubroutine LBLMC_precompute
   
   
   subroutine generate_fsck(&
      T_array,xs_array,Tb_array,id_spec,g_file,info_file,&              !Mandatory parameters
      k_array,k_max,k_min,k_points,k_distribution,g_points,&
      formatted_output,albdf,pressure_based,tabulate_k,&
      wb_lbound,wb_ubound)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckMemAlloc,dprint,&
                  get_file_unit,get_wall_time,print_to_prompt,shutdown
      use constants, only: atm_pa,Ru
      use global_parameters, only: id_soot
      use math_functions, only: generate_distribution,linint,locate,&
                                spline,splint
      use omp_lib
      use physical_functions, only: Ib_function,kappa_eta_soot
      use precision_parameters
      implicit none
      character(*),intent(in) :: g_file,info_file
      character(*),intent(in),optional :: k_distribution
      character(200) :: msg,which_distribution
      integer,intent(in) :: id_spec
      integer,intent(in),optional :: k_points,g_points
      integer :: ierr,g_unit,info_unit,npts
      integer :: ngv,nlbl,nTg,nTb,ntp,nwb,nxs,nxv
      integer :: igv,ilbl,iTg,iTb,itp,iwb,ixs,ixv
      logical,intent(in),optional :: albdf,formatted_output,&
         pressure_based,tabulate_k
      logical :: binary_output,compute_albdf,compute_fsckd
      logical :: full_k_range,k_from_dist,k_is_provided
      logical :: k_is_x,g_is_x
      logical :: pbased
      real(dp),intent(in) :: T_array(:),Tb_array(:),xs_array(:)
      real(dp) :: kappa,ktoC,xloc
      real(dp),intent(in),optional :: k_array(:),k_max,k_min,&
                                      wb_lbound(:),wb_ubound(:)
      real(dp),allocatable,dimension(:) :: Ibeta,wb_low,wb_up
      real(dp),allocatable,dimension(:) :: kaux,gx
      real(dp),allocatable,dimension(:) :: p_array
      real(dp),allocatable,dimension(:,:) :: Ib,Ib_omp
      real(dp),allocatable,dimension(:,:,:) :: k_realmax,k_realmin,&
                                               omp_max,omp_min
      real(dp),allocatable,dimension(:,:,:,:) :: x_array
      real(dp),allocatable,dimension(:,:,:,:,:) :: ky,g,g_omp
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Preparatory procedures')
      
      !Define if we are calculating the the 
      !FSCK distribution (default) or the ALBDF
      compute_albdf = .false.; if (present(albdf)) compute_albdf = albdf
      compute_fsckd = .not.compute_albdf
      
      !Set up output format (by default, all output files are in binary)
      binary_output = .true.
      if (present(formatted_output)) binary_output = .not.formatted_output
      
      !Define wheter k or g are to be tabulated
      g_is_x = .false.; if (present(tabulate_k)) g_is_x = tabulate_k
      k_is_x = .not.g_is_x
      
      !Set up how the discrete k (or C) values are defined
      k_is_provided = .false.
      k_from_dist = .false.
      full_k_range = .false.
      if (g_is_x) then
         full_k_range = .true.
      elseif (present(k_array)) then
         k_is_provided = .true.
      elseif (present(k_max).and.present(k_min).and.&
             present(k_points)) then
         k_from_dist = .true.
      else
         full_k_range = .true.
      endif
    
      !Set up default parameters for the k/C discrete values that
      !will be used to generate  the FSKd or ALBDF
      nxv = 5000; if (present(k_points)) nxv = k_points                 !Number of k or C points
      ngv = 64; if (present(g_points)) ngv = g_points                   !Number of g or F points (if tabulate_k = .true.)
      which_distribution = 'wang'                                       !Function that describes how the k or C points
      if (present(k_distribution)) &                                    !  are distributed
         which_distribution = trim(k_distribution)
      
      !Set up default parameters for the wide bands
      nwb = 1
      if (present(wb_lbound).and.present(wb_ubound)) then
         call assert(size(wb_lbound).eq.size(wb_ubound))
         nwb = size(wb_lbound)
      endif
      
      !Set up default flag
      pbased = .false.
      if (present(pressure_based)) pbased = pressure_based
      if (compute_albdf) pbased = .true.
      
      !Number of k/C points in the tabulated data
      npts = nxv; if (g_is_x) npts = ngv
      
      !-----------------------------------------------------------------
      !Load the LBL data
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Load the LBL data')
      call load_lbl_data
      
      !-----------------------------------------------------------------
      !Define array sizes
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Define array sizes')
      nlbl = lbl_nlines
      nTg = size(T_array)
      nTb = size(Tb_array)
      nxs = size(xs_array)
      
      !-----------------------------------------------------------------
      !Set up artificial pressure array (this is temporary)
      !-----------------------------------------------------------------
      ntp = 1
      allocate(p_array(ntp))
      p_array = 1._dp
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Allocate arrays')
      allocate(Ib(nTb,nwb),stat=ierr)
      call CheckMemAlloc('Ib',ierr)
      allocate(Ib_omp(nTb,nwb),stat=ierr)
      call CheckMemAlloc('Ib_omp',ierr)
      allocate(Ibeta(nTb),stat=ierr)
      call CheckMemAlloc('Ibeta',ierr)
      allocate(g(nxv,nTg,nxs,nTb,nwb),stat=ierr)
      call CheckMemAlloc('g',ierr)
      allocate(g_omp(nxv,nTg,nxs,nTb,nwb),stat=ierr)
      call CheckMemAlloc('g_omp',ierr)
      allocate(k_realmax(nTg,nxs,nwb),stat=ierr)
      call CheckMemAlloc('k_realmax',ierr)
      allocate(k_realmin(nTg,nxs,nwb),stat=ierr)
      call CheckMemAlloc('k_realmin',ierr)
      allocate(omp_max(nTg,nxs,nwb),stat=ierr)
      call CheckMemAlloc('omp_max',ierr)
      allocate(omp_min(nTg,nxs,nwb),stat=ierr)
      call CheckMemAlloc('omp_min',ierr)
      allocate(x_array(nxv,nTg,nxs,nwb),stat=ierr)
      call CheckMemAlloc('x_array',ierr)
      allocate(wb_low(nwb),stat=ierr)
      call CheckMemAlloc('wb_low',ierr)
      allocate(wb_up(nwb),stat=ierr)
      call CheckMemAlloc('wb_up',ierr)
!      allocate(test(nxv),stat=ierr)
!      call CheckMemAlloc('',ierr)

      !-----------------------------------------------------------------
      !Set up wide band bounds
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Set up wide band bounds')
      wb_low = small; wb_up = big                                       !Default values (consider the full spectrum)
      if (present(wb_lbound)) wb_low = wb_lbound
      if (present(wb_ubound)) wb_up = wb_ubound

      !-----------------------------------------------------------------
      !Compute the maximum and minimum values of k (or C), if needed
      !-----------------------------------------------------------------
      if (full_k_range) then
         call dprint('generate_fsck: Compute the maximum and minimum &
                     &values of k')
         k_realmax = -big; k_realmin = big
         
         !$OMP PARALLEL DEFAULT(SHARED) &
         !$OMP&   PRIVATE(itg,itp,ixs,iwb,kappa,ktoC,omp_max,omp_min,&
         !$OMP&           xloc)
         omp_max = -big; omp_min = big
         !$OMP DO
         do ilbl=1,nlbl
         
            !Find the appropriate wide-band index
            if (lbl_xeta(ilbl).lt.(wb_low(1))) cycle
            if (lbl_xeta(ilbl).gt.(wb_up(nwb))) cycle
            do iwb=1,nwb
               if ((lbl_xeta(ilbl).ge.wb_low(iwb)).and.&
                  (lbl_xeta(ilbl).lt.wb_up(iwb))) exit
            enddo
         
            do itg=1,ntg
               do ixs=1,nxs
                  do itp=1,ntp

                     !Absorption coefficient
                     if (id_spec.eq.id_soot) then
                        kappa = kappa_eta_soot(lbl_xeta(ilbl),&
                           xs_array(ixs))
                        if (pbased) kappa = kappa/xs_array(ixs)
                     else
                        kappa = lbl_kappa_func(xs_array(ixs),&
                           T_array(iTg),p_array(itp),id_spec,ilbl,&
                           kappa_p=pbased)
                     endif
                     if (compute_fsckd) ktoC = 1._dp
                     if (compute_albdf) ktoC = atm_pa/(T_array(iTg)*Ru)
                     xloc = kappa/ktoC

                     if (xloc.gt.omp_max(itg,ixs,iwb)) &
                        omp_max(itg,ixs,iwb) = xloc
                     if (xloc.lt.omp_min(itg,ixs,iwb)) &
                        omp_min(itg,ixs,iwb) = xloc
                  enddo
               enddo
            enddo   
         enddo
         !$OMP ENDDO
           
         !$OMP CRITICAL
         do itg=1,ntg
            do ixs=1,nxs
               do iwb=1,nwb
                  if (omp_max(itg,ixs,iwb).gt.k_realmax(itg,ixs,iwb)) &
                     k_realmax(itg,ixs,iwb) = omp_max(itg,ixs,iwb)
                  if (omp_min(itg,ixs,iwb).lt.k_realmin(itg,ixs,iwb)) &
                     k_realmin(itg,ixs,iwb) = omp_min(itg,ixs,iwb)
               enddo
            enddo
         enddo
         !$OMP ENDCRITICAL
         !$OMP ENDPARALLEL
         
         call rewind_lbl_data
      endif

      !-----------------------------------------------------------------
      !Define the discrete k (or C) values to be used in the generation
      !of the full-spectrum CK distributions or ALBDF
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Define the discrete k values')
      do itg=1,ntg
         do ixs=1,nxs
            do iwb=1,nwb
               !Option 1: k_array is provided as input
               if (k_is_provided) x_array(:,iTg,ixs,iwb) = k_array
               
               !Option 2: the parameters of the distribution (maximum 
               !and minimum values and number of points) are provided
               !as input
               if (k_from_dist) &
                  call generate_distribution(which_distribution,k_min,&
                                       k_max,nxv,x_array(:,iTg,ixs,iwb))
               
               !Option 3: the full range of possible k values at each 
               !of the speficied thermodynamic states is used
               if (full_k_range) &
                  call generate_distribution(which_distribution,&
                     k_realmin(iTg,ixs,iwb),k_realmax(iTg,ixs,iwb),nxv,&
                     x_array(:,iTg,ixs,iwb))
            enddo
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Initialize sum arrays
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Initialize sum arrays')
      g = 0._dp; Ib = 0._dp
      
      !-----------------------------------------------------------------
      !Main calculation loop
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Main calculation loop')
      !$OMP PARALLEL DEFAULT(SHARED) &
      !$OMP&   PRIVATE(g_omp,Ib_omp,Ibeta,itb,itg,itp,ixs,ixv,iwb,&
      !$OMP&           kappa,ktoC,xloc)
      
      !Initialize omp sum arrays
      g_omp = 0; Ib_omp = 0._dp
      
      !$OMP DO
      lbl_loop: do ilbl=1,nlbl
      
         !Print xeta if requested
         if (print_xeta) then
            write(msg,'(f12.4)') lbl_xeta(ilbl)/100._dp
            msg = 'Wavenumber position [1/cm]: '//trim(adjustl(msg))
            call print_to_prompt(msg,3)
         endif
      
         !Find the appropriate wide-band index
         if (lbl_xeta(ilbl).lt.(wb_low(1))) cycle lbl_loop
         if (lbl_xeta(ilbl).gt.(wb_up(nwb))) cycle lbl_loop
         do iwb=1,nwb
            if ((lbl_xeta(ilbl).ge.wb_low(iwb)).and.&
               (lbl_xeta(ilbl).lt.wb_up(iwb))) exit
         enddo
      
         !Compute the blackbody intensity
         tb_loop_0: do itb=1,ntb
            Ibeta(itb) = &
               Ib_function(Tb_array(itb),lbl_xeta(ilbl))*lbl_deta(ilbl)
            Ib_omp(itb,iwb) = Ib_omp(itb,iwb) + Ibeta(itb)
         enddo tb_loop_0

         tg_loop_1: do itg=1,ntg
            xs_loop_1: do ixs=1,nxs
               tp_loop_1: do itp=1,ntp
                  !Compute the absorption coefficient
                  if (id_spec.eq.id_soot) then
                     kappa = kappa_eta_soot(lbl_xeta(ilbl),&
                           xs_array(ixs))
                     if (pbased) kappa = kappa/xs_array(ixs)
                  else
                     kappa = lbl_kappa_func(xs_array(ixs),T_array(iTg),&
                               p_array(itp),id_spec,ilbl,kappa_p=pbased)
                  endif

                  if (compute_fsckd) ktoC = 1._dp
                  if (compute_albdf) ktoC = atm_pa/(T_array(iTg)*Ru)
                  xloc = kappa/ktoC

                  !Find where to add to the Ib*deta integration
                  ixv = locate(x_array(:,itg,ixs,iwb),xloc,nxv) + 1
                  if (ixv.gt.nxv) cycle tp_loop_1
!                  if (xloc.lt.x_array(1,itg,ixs,iwb)) then
!                     ixv = 1
!                  elseif (xloc.gt.x_array(nxv,itg,ixs,itp)) then
!                     i=nxv
!                  else
!                  endif 
                  
                  !Add to the Ib*deta integration
                  tb_loop_1: do itb=1,ntb
                     g_omp(ixv,itg,ixs,itb,iwb) = &
                        g_omp(ixv,itg,ixs,itb,iwb) + Ibeta(itb)
                  enddo tb_loop_1
               enddo tp_loop_1
            enddo xs_loop_1
         enddo tg_loop_1
      enddo lbl_loop
      !$OMP ENDDO
      
      !-----------------------------------------------------------------
      !Finish the calculation of the total blackbody intensity
      !-----------------------------------------------------------------
      !$OMP CRITICAL
      Ib = Ib + Ib_omp
      g = g + g_omp
      !$OMP ENDCRITICAL
      
      !$OMP ENDPARALLEL

      !-----------------------------------------------------------------
      !Finish calculation of g
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Finish calculation of g')
      tg_loop_2: do itg=1,ntg
         xs_loop_2: do ixs=1,nxs
            tp_loop_2: do itp=1,ntp
               tb_loop_2: do itb=1,ntb
                  wb_loop: do iwb=1,nwb
                     x_loop: do ixv=2,nxv
                        g(ixv,itg,ixs,itb,iwb) = &
                           g(ixv,itg,ixs,itb,iwb) + &
                              g(ixv-1,itg,ixs,itb,iwb)
                     enddo x_loop
                     g(:,itg,ixs,itb,iwb) = &
                        g(:,itg,ixs,itb,iwb)/Ib(itb,iwb)
                  enddo wb_loop
               enddo tb_loop_2
            enddo tp_loop_2
         enddo xs_loop_2
      enddo tg_loop_2
      
      !-----------------------------------------------------------------
      !If needed, invert the data as to make k(g), rather than g(k)
      !-----------------------------------------------------------------
      if (g_is_x) then
         call dprint('generate_fsck: Invert the data as to make k(g), &
                     &rather than g(k)')

         !Allocate additional arrays
         allocate(kaux(nxv),stat=ierr)
         call CheckMemAlloc('kaux',ierr)
         allocate(ky(ngv,ntg,nxs,ntb,nwb),stat=ierr)
         call CheckMemAlloc('ky',ierr)
         allocate(gx(ngv),stat=ierr)
         call CheckMemAlloc('gx',ierr)
         
         !Define the discrete g positions
         call generate_distribution('gauss-chebyshev-even',&
                                    0._dp,1._dp,ngv,gx)
         
         !Determine k from g throught spline interpolations      
         do iTg=1,nTg
            do ixs=1,nxs
               do iTb=1,nTb
                  do iwb=1,nwb
!                     call spline(g(:,itg,ixs,itb,iwb),&
!                        x_array(:,itg,ixs,iwb),nxv,big,big,kaux)
!do ixv=1,nxv
!   write(*,*) ixv,g(ixv,itg,ixs,itb,iwb),x_array(ixv,itg,ixs,iwb)
!enddo
!write(*,*) T_array(iTg),Tb_array(itb)
                     do igv=1,ngv
                        ky(igv,itg,ixs,itb,iwb) = &
                           linint(g(:,itg,ixs,itb,iwb),&
                              x_array(:,itg,ixs,iwb),gx(igv),nxv)
!write(*,*) ky(igv,itg,ixs,itb,iwb)
!                        ky(igv,itg,ixs,itb,iwb) = &
!                           splint(g(:,itg,ixs,itb,iwb),&
!                                x_array(:,itg,ixs,iwb),kaux,gx(igv),nxv)

                     enddo
                  enddo
               enddo
            enddo
         enddo
      endif

      !-----------------------------------------------------------------
      !Dump the FSCKd/ALBDF data
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Dump data')
      if (binary_output) then
         !Dump g data
         g_unit = get_file_unit()
         open(file=g_file,unit=g_unit,action='write',form='unformatted')
         if (k_is_x) write(g_unit) g
         if (g_is_x) write(g_unit) ky
         close(g_unit)
         
         !Dump info
         info_unit = get_file_unit()
         open(file=info_file,unit=info_unit,&
            action='write',form='unformatted')
         write(info_unit) g_is_x
         if (compute_albdf) write(info_unit) nxs,nTg,nTb,npts,nwb
         if (.not.compute_albdf) &
            write(info_unit) nxs,nTg,nTb,nxv,nwb,pbased
         write(info_unit) xs_array
         if (k_is_x) then
            if (full_k_range) write(info_unit) x_array
            if (k_is_provided.or.k_from_dist) &
               write(info_unit) x_array(:,1,1,1)
         endif
         if (g_is_x) write(info_unit) gx
         write(info_unit) Tb_array
         write(info_unit) T_array
         write(info_unit) wb_low
         write(info_unit) wb_up
         close(info_unit)
         
      else 
         call shutdown('generate_fsck: formatted output not available')
      endif
   
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('generate_fsck: Deallocate arrays')
      deallocate(Ib,Ib_omp,Ibeta)
      deallocate(g)
      deallocate(k_realmax,k_realmin)
      deallocate(omp_max,omp_min)
      deallocate(wb_low,wb_up)
      if (g_is_x) deallocate(kaux,ky,gx)
      
   endsubroutine generate_fsck

   !====================================================================
   !Subroutine that generates narrow-band cumulative-k distributions 
   !from a LBL database and dumps to an external file
   !====================================================================
   subroutine generate_nbck(&
      T_array,xs_array,lnb_array,unb_array,id_spec,g_file,&
      error_file,g_points,g_distribution,path_lengths,pressure,run_time)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
         dprint,get_file_unit,get_wall_time,shutdown
      use math_functions, only: generate_distribution,linint,locate,&
                                polint,sort_heap
      use lbl_parameters
      use omp_lib
      use precision_parameters, only: big,dp,small
      implicit none
      character(*),intent(in) :: g_file
      character(*),optional,intent(in) :: error_file,g_distribution
      character(15) :: Tstr,xstr
      character(200) :: msg,which_distribution
      integer,intent(in) :: id_spec
      integer,optional,intent(in) :: g_points
      integer :: error_unit,g_unit,ierr
      integer :: i,ilb,inb,ipl,itg,ixs
      integer :: nlb,nnb,npl,ntg,nxs
      integer :: llb,ulb
      integer :: g_total,k_total
      logical :: dump_error
      real(dp),intent(in) :: lnb_array(:),T_array(:),unb_array(:),&
         xs_array(:)
      real(dp),optional,intent(in) :: path_lengths(:),pressure
      real(dp),optional,intent(out) :: run_time
      real(dp) :: nb_deta,pres
      real(dp) :: dg,kappa
      real(dp) :: start_time,end_time
      real(dp),allocatable,dimension(:) :: g_array,g_out,k_array,plength
      real(dp),allocatable,dimension(:,:,:) :: ekappa,kappa_ck,kappa_lbl
      real(dp),allocatable,dimension(:,:,:,:) :: k_out
      real(dp),allocatable,dimension(:,:,:,:) :: etau,tau_ck,tau_lbl
      
      !-----------------------------------------------------------------
      !Set up initial parameters
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Set up initial parameters')
      start_time = get_wall_time()                                      !Start time counter
      
      !Set up optional parameters
      which_distribution = 'gauss-chebyshev-even'                       !Type of g distribution in the
      if (present(g_distribution)) which_distribution = g_distribution  !  output file
      g_total = 100; if (present(g_points)) g_total = g_points          !Number of points in the output file
      pres = 1._dp; if (present(pressure)) pres = pressure              !Total pressure
      
      !Set up flags
      dump_error = .false.; if (present(error_file)) dump_error = .true.
      
      !-----------------------------------------------------------------
      !Load LBL data files
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Load LBL data files')
      call load_lbl_data

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Define surrogate names')
      nnb = size(lnb_array)                                             !Total number of narrow bands
      call assert(nnb.eq.size(unb_array),&
                  'size(lnb_array)=size(unb_array)')                    !Check if the sizes of leta_in and ueta_in match
      ntg = size(T_array)                                               !Number of discrete temperature values
      nxs = size(xs_array)                                              !Number of discrete mole fraction values
                                                                        !  (same for all species)
      nlb = lbl_nlines                                                  !Total number of lines in the external LBL data files

      !-----------------------------------------------------------------
      !Set up path lengths for computing the transmittance
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Set up path lengths')
      if (present(path_lengths)) then
         npl = size(path_lengths); allocate(plength(npl))
         plength = path_lengths
      else
         npl = 5; allocate(plength(npl))
         plength(1) = 0.01_dp
         do ipl=2,npl
            plength(ipl) = plength(ipl-1)*10._dp
         enddo
      endif

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Allocate arrays')
      allocate(ekappa(nnb,ntg,nxs),stat=ierr)
      call CheckMemAlloc('ekappa',ierr)
      allocate(etau(nnb,ntg,nxs,npl),stat=ierr)
      call CheckMemAlloc('etau',ierr)
      allocate(g_out(g_total),stat=ierr)
      call CheckMemAlloc('g_out',ierr)
      allocate(k_out(g_total,ntg,nxs,nnb),stat=ierr)
      call CheckMemAlloc('k_out',ierr)
      allocate(kappa_ck(nnb,ntg,nxs),stat=ierr)
      call CheckMemAlloc('kappa_ck',ierr)
      allocate(kappa_lbl(nnb,ntg,nxs),stat=ierr)
      call CheckMemAlloc('kappa_lbl',ierr)
      allocate(tau_ck(nnb,ntg,nxs,npl),stat=ierr)
      call CheckMemAlloc('tau_ck',ierr)
      allocate(tau_lbl(nnb,ntg,nxs,npl),stat=ierr)
      call CheckMemAlloc('tau_lbl',ierr)

      !-----------------------------------------------------------------
      !Generate reference g distribution
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Generate reference g distribution')
      call generate_distribution(which_distribution,0._dp,1._dp,&
                                 g_total,g_out)

      mole_fraction_loop: do ixs=1,nxs
         write(xstr,'(f6.4)') xs_array(ixs)
         xstr = adjustl(trim(xstr))
            
         temperature_loop: do itg=1,ntg
            write(Tstr,'(f12.2)') T_array(itg)
            Tstr = adjustl(trim(Tstr))

            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP&   PRIVATE(dg,g_array,i,ilb,ipl,k_array,k_total,&
            !$OMP&           kappa,llb,msg,nb_deta,ulb) 
            narrow_band_loop: do inb=1,nnb
               !Printing status
               write(msg,'(a,a,a,a,a,i3)') &
                  'generate_nbck: Main loop started. xs = ',trim(xstr),&
                  '; T = ',trim(Tstr),' K; band ',inb
               call dprint(msg)
               
               !--------------------------------------------------------
               !Preparatory procedures
               !--------------------------------------------------------
               call dprint('generate_nbck: Preparatory procedures',3)
       
               !Indexes of the LBL data corresponding to the 
               !positions of the  lower and upper bounds of the
               !current narrow band
               llb = locate(lbl_xeta,lnb_array(inb),nlb) + 1
               ulb = locate(lbl_xeta,unb_array(inb),nlb)
               k_total = ulb - llb + 1

               allocate(k_array(k_total))
               allocate(g_array(k_total))

               !--------------------------------------------------------
               !Find appropriate maximum and minimum k
               !--------------------------------------------------------
               call dprint('generate_nbck: &
                           &Find maximum and minimum k for the band',3)
               
               kappa_lbl(inb,itg,ixs) = 0._dp
               tau_lbl(inb,itg,ixs,:) = 0._dp
               nb_deta = 0._dp
               do ilb=llb,ulb
                  !Absorption coefficient
                  kappa = lbl_kappa_func(xs_array(ixs),T_array(itg),&
                                        pres,id_spec,ilb,kappa_p=.true.)
                  k_array(ilb-llb+1) = kappa
                  kappa_lbl(inb,itg,ixs) = kappa_lbl(inb,itg,ixs) + &
                     kappa*lbl_deta(ilb)

                  !Compute band width
                  nb_deta = nb_deta + lbl_deta(ilb)

                  !Compute the LBL band-average transmittance
                  do ipl=1,npl
                     tau_lbl(inb,itg,ixs,ipl) = &
                        tau_lbl(inb,itg,ixs,ipl) + &
                           dexp(-kappa*plength(ipl))*lbl_deta(ilb)
                  enddo
               enddo
               
               !Finish computing the absorption coefficient and the
               !transmittance for the LBL method
               kappa_lbl(inb,itg,ixs) = kappa_lbl(inb,itg,ixs)/nb_deta
               tau_lbl(inb,itg,ixs,:) = tau_lbl(inb,itg,ixs,:)/nb_deta
               
               !--------------------------------------------------------
               !Sort the spectral absorption coefficient data for the
               !narrow band in increasing order
               !--------------------------------------------------------
               call sort_heap(k_array)

               !--------------------------------------------------------
               !Compute g
               !--------------------------------------------------------
               !Zero out the sum arrays
               g_array = 0._dp
            
               !Check where to add to the eta integration
               do ilb=llb,ulb
                  i = ilb - llb + 1
                  g_array(i) = g_array(i) + lbl_deta(ilb)
               enddo
               g_array = g_array/nb_deta
                     
               !Finish up the g calculation
               do i=2,k_total
                  g_array(i) = min(g_array(i-1) + g_array(i),1._dp)
               enddo

               !--------------------------------------------------------
               !Invert k(g) dependence
               !--------------------------------------------------------
               do i=1,g_total
                  call polint(g_array,k_array,g_out(i),&
                              k_out(i,itg,ixs,inb),2)
!if (inb.eq.83) write(*,*) xs_array(ixs),T_array(itg),g_out(i),k_out(i,itg,ixs,inb)
!write(*,*) i,g_out(i),k_out(i,itg,ixs,inb)
               enddo
               deallocate(g_array,k_array)
               
               !--------------------------------------------------------
               !Compute the absorption coefficient and 
               !the transmittance from the k-distributions
               !--------------------------------------------------------
               kappa_ck(inb,itg,ixs) = 0._dp
               tau_ck(inb,itg,ixs,:) = 0._dp
               do i=1,g_total
                  !Absorption coefficient
                  kappa = k_out(i,itg,ixs,inb)
                        
                  !Quadrature weight
                  if (i.eq.1) then
                     dg = g_out(i)
                  elseif (i.eq.g_total) then
                     dg = 1._dp - g_out(i)
                  else
                     dg = 0.5_dp*(g_out(i+1) - g_out(i-1))
                  endif
                  
                  !Band-average absorption coefficient   
                  kappa_ck(inb,itg,ixs) = &
                     kappa_ck(inb,itg,ixs) + kappa*dg
                  
                  !Band-average transmittance
                  do ipl=1,npl
                     tau_ck(inb,itg,ixs,ipl) = &
                        tau_ck(inb,itg,ixs,ipl) + &
                           dexp(-kappa*plength(ipl))*dg
                  enddo
               enddo

               !--------------------------------------------------------
               !Compute errors
               !--------------------------------------------------------
               ekappa(inb,itg,ixs) = dabs(kappa_ck(inb,itg,ixs) - kappa_lbl(inb,itg,ixs))/&
                                               (kappa_lbl(inb,itg,ixs)+small)*100._dp
               do ipl=1,npl
                  etau(inb,itg,ixs,ipl) = &
                     dabs(tau_ck(inb,itg,ixs,ipl) - tau_lbl(inb,itg,ixs,ipl))/&
                        (tau_lbl(inb,itg,ixs,ipl) + small)*100._dp
               enddo
            enddo narrow_band_loop
            !$OMP END PARALLEL DO
         enddo temperature_loop
      enddo mole_fraction_loop

      !-----------------------------------------------------------------
      !Dump data
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Dump data')
      g_unit = get_file_unit()
      open(file=g_file,unit=g_unit,action='write',form='unformatted')
      write(g_unit) g_total,nnb,ntg,nxs
      write(g_unit) lnb_array,unb_array
      write(g_unit) T_array
      write(g_unit) xs_array
      write(g_unit) g_out
      write(g_unit) k_out
      close(g_unit)

      !-----------------------------------------------------------------
      !Dump errors
      !-----------------------------------------------------------------
      if (dump_error) then
         error_unit = get_file_unit()
         open(file=error_file,unit=error_unit,action='write',&
              form='formatted')
         write(error_unit,'(100(a,:,","))') 'T','xs','lbound','ubound',&
            'kappa','kappa','kappa',('tau',ipl=1,3*npl)
         write(error_unit,'(100(a,:,","))') '','','','','LBL','CK',&
            'Error',('LBL','CK','Error',ipl=1,npl)   
         write(error_unit,'(7(a,:,","),100(e26.15e3,:,","))') &
            '','','','','','','',(plength(ipl),plength(ipl),&
            plength(ipl),ipl=1,npl)
         do itg=1,ntg
            do ixs=1,nxs
               do inb=1,nnb
                  write(error_unit,'(100(e26.15e3,:,","))') &
                     T_array(itg),xs_array(ixs),lnb_array(inb),&
                     unb_array(inb),kappa_lbl(inb,itg,ixs),&
                     kappa_ck(inb,itg,ixs),ekappa(inb,itg,ixs),&
                     (tau_lbl(inb,itg,ixs,ipl),tau_ck(inb,itg,ixs,ipl),&
                     etau(inb,itg,ixs,ipl),ipl=1,npl)
               enddo
            enddo
         enddo
         close(error_unit)
      endif

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Deallocate arrays')
      deallocate(ekappa,etau)
      deallocate(g_out,k_out)
      deallocate(kappa_ck,kappa_lbl)
      deallocate(tau_ck,tau_lbl)

      !-----------------------------------------------------------------
      !Compute total run time
      !-----------------------------------------------------------------
      end_time = get_wall_time()                                        !Stop time counter
      if (present(run_time)) run_time = end_time - start_time           !Compute run time, if requested

   endsubroutine generate_nbck

   !====================================================================
   !Subroutine that generates narrow-band cumulative-k distributions 
   !from a LBL database for a H2O-CO2 mixture and dumps to an external 
   !file
   !====================================================================
   subroutine generate_premix_nbck(&
      T_array,xs_array,lnb_array,unb_array,g_file,&
      error_file,g_points,g_distribution,path_lengths,pressure,run_time)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
         dprint,get_file_unit,get_wall_time,shutdown
      use global_parameters, only: id_co2,id_h2o,id_soot
      use math_functions, only: generate_distribution,linint,locate,&
                                polint,sort_heap
      use lbl_parameters
      use omp_lib
      use physical_functions, only: kappa_eta_soot
      use precision_parameters, only: big,dp,small
      implicit none
      character(*),intent(in) :: g_file
      character(*),optional,intent(in) :: error_file,g_distribution
      character(15) :: fvstr,Tstr,xcstr,xwstr
      character(200) :: msg,which_distribution
      integer,optional,intent(in) :: g_points
      integer :: error_unit,g_unit,ierr
      integer :: i,ifv,ilb,inb,ipl,itg,ixc,ixw
      integer :: nfv,nlb,nnb,npl,ntg,nxc,nxw
      integer :: llb,ulb
      integer :: g_total,k_total
      logical :: dump_error
      real(dp),intent(in) :: lnb_array(:),T_array(:),unb_array(:),&
         xs_array(:,:)
      real(dp),optional,intent(in) :: path_lengths(:),pressure
      real(dp),optional,intent(out) :: run_time
      real(dp) :: nb_deta,pres
      real(dp) :: dg,kappa,kappa_c,kappa_s,kappa_w
      real(dp) :: start_time,end_time
      real(dp),allocatable,dimension(:) :: g_array,g_out,k_array,plength
      real(dp),allocatable,dimension(:,:,:,:,:) :: ekappa,kappa_ck,kappa_lbl
      real(dp),allocatable,dimension(:,:,:,:,:,:) :: k_out
      real(dp),allocatable,dimension(:,:,:,:,:,:) :: etau,tau_ck,tau_lbl
      
      !-----------------------------------------------------------------
      !Set up initial parameters
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Set up initial parameters')
      start_time = get_wall_time()                                      !Start time counter
      
      !Set up optional parameters
      which_distribution = 'gauss-chebyshev-even'                       !Type of g distribution in the
      if (present(g_distribution)) which_distribution = g_distribution  !  output file
      g_total = 100; if (present(g_points)) g_total = g_points          !Number of points in the output file
      pres = 1._dp; if (present(pressure)) pres = pressure              !Total pressure
      
      !Set up flags
      dump_error = .false.; if (present(error_file)) dump_error = .true.
      
      !-----------------------------------------------------------------
      !Load LBL data files
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Load LBL data files')
      call load_lbl_data

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Define surrogate names')
      nnb = size(lnb_array)                                             !Total number of narrow bands
      call assert(nnb.eq.size(unb_array),&
                  'size(lnb_array)=size(unb_array)')                    !Check if the sizes of leta_in and ueta_in match
      ntg = size(T_array)                                               !Number of discrete temperature values
      nxc = size(xs_array(:,id_co2))                                    !Number of discrete mole fraction values for CO2
      nxw = size(xs_array(:,id_h2o))                                    !Number of discrete mole fraction values for H2O
      nfv = size(xs_array(:,id_soot))                                   !Number of discrete soot volume fraction values
      nlb = lbl_nlines                                                  !Total number of lines in the external LBL data files

      !-----------------------------------------------------------------
      !Set up path lengths for computing the transmittance
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Set up path lengths')
      if (present(path_lengths)) then
         npl = size(path_lengths); allocate(plength(npl))
         plength = path_lengths
      else
         npl = 5; allocate(plength(npl))
         plength(1) = 0.01_dp
         do ipl=2,npl
            plength(ipl) = plength(ipl-1)*10._dp
         enddo
      endif

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Allocate arrays')
      allocate(ekappa(nnb,ntg,nxw,nxc,nfv),stat=ierr)
      call CheckMemAlloc('ekappa',ierr)
      allocate(etau(nnb,ntg,nxw,nxc,nfv,npl),stat=ierr)
      call CheckMemAlloc('etau',ierr)
      allocate(g_out(g_total),stat=ierr)
      call CheckMemAlloc('g_out',ierr)
      allocate(k_out(g_total,ntg,nxw,nxc,nfv,nnb),stat=ierr)
      call CheckMemAlloc('k_out',ierr)
      allocate(kappa_ck(nnb,ntg,nxw,nxc,nfv),stat=ierr)
      call CheckMemAlloc('kappa_ck',ierr)
      allocate(kappa_lbl(nnb,ntg,nxw,nxc,nfv),stat=ierr)
      call CheckMemAlloc('kappa_lbl',ierr)
      allocate(tau_ck(nnb,ntg,nxw,nxc,nfv,npl),stat=ierr)
      call CheckMemAlloc('tau_ck',ierr)
      allocate(tau_lbl(nnb,ntg,nxw,nxc,nfv,npl),stat=ierr)
      call CheckMemAlloc('tau_lbl',ierr)

      !-----------------------------------------------------------------
      !Generate reference g distribution
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Generate reference g distribution')
      call generate_distribution(which_distribution,0._dp,1._dp,&
                                 g_total,g_out)

      h2o_loop: do ixw=1,nxw
         write(xwstr,'(f6.4)') xs_array(ixw,id_h2o)
         xwstr = adjustl(trim(xwstr))
         
         co2_loop: do ixc=1,nxc
            write(xcstr,'(f6.4)') xs_array(ixc,id_co2)
            xcstr = adjustl(trim(xcstr))
            
            soot_loop: do ifv=1,nfv
               write(fvstr,'(f6.4)') xs_array(ifv,id_soot)
               fvstr = adjustl(trim(fvstr))       
               temperature_loop: do itg=1,ntg
                  write(Tstr,'(f12.2)') T_array(itg)
                  Tstr = adjustl(trim(Tstr))
write(*,*) ixw,ixc,ifv,itg 

                  !$OMP PARALLEL DO DEFAULT(SHARED) &
                  !$OMP&   PRIVATE(dg,g_array,i,ilb,ipl,k_array,&
                  !$OMP&           k_total,kappa,kappa_c,kappa_s,&
                  !$OMP&           kappa_w,llb,msg,nb_deta,ulb) 
                  narrow_band_loop: do inb=1,nnb
                     !Printing status
                     write(msg,'(a,a,a,a,a,a,a,a,a,i3)') &
                        'generate_nbck: Main loop started. xw = ',&
                        trim(xwstr),'; xc = ',trim(xcstr),'; fv = ',&
                        trim(fvstr),'; T = ',trim(Tstr),' K; band ',inb
                     call dprint(msg)
               
                     !--------------------------------------------------
                     !Preparatory procedures
                     !--------------------------------------------------
                     call dprint('generate_premix_nbck: &
                                 &Preparatory procedures',3)
       
                     !Indexes of the LBL data corresponding to the 
                     !positions of the  lower and upper bounds of the
                     !current narrow band
                     llb = locate(lbl_xeta,lnb_array(inb),nlb) + 1
                     ulb = locate(lbl_xeta,unb_array(inb),nlb)
                     k_total = ulb - llb + 1

                     allocate(k_array(k_total))
                     allocate(g_array(k_total))
   
                     !--------------------------------------------------
                     !Build input kappa array for the band
                     !--------------------------------------------------
                     call dprint('generate_premix_nbck: &
                        &Build input kappa array for the band',3)
               
                     kappa_lbl(inb,itg,ixw,ixc,ifv) = 0._dp
                     tau_lbl(inb,itg,ixw,ixc,ifv,:) = 0._dp
                     nb_deta = 0._dp
                     do ilb=llb,ulb                     
                        !Absorption coefficient
                        kappa_w = lbl_kappa_func(xs_array(ixw,id_h2o),&
                           T_array(itg),pres,id_h2o,ilb)
                        kappa_c = lbl_kappa_func(xs_array(ixc,id_co2),&
                           T_array(itg),pres,id_co2,ilb)
                        kappa_s = kappa_eta_soot(lbl_xeta(ilb),&
                                                 xs_array(ifv,id_soot))
                        kappa = kappa_w + kappa_c + kappa_s
                        k_array(ilb-llb+1) = kappa
                        kappa_lbl(inb,itg,ixw,ixc,ifv) = &
                           kappa_lbl(inb,itg,ixw,ixc,ifv) + &
                              kappa*lbl_deta(ilb)

                        !Compute band width
                        nb_deta = nb_deta + lbl_deta(ilb)

                        !Compute the LBL band-average transmittance
                        do ipl=1,npl
                           tau_lbl(inb,itg,ixw,ixc,ifv,ipl) = &
                              tau_lbl(inb,itg,ixw,ixc,ifv,ipl) + &
                                 dexp(-kappa*plength(ipl))*lbl_deta(ilb)
                        enddo
                     enddo
               
                     !Finish computing the absorption coefficient and 
                     !the transmittance for the LBL method
                     kappa_lbl(inb,itg,ixw,ixc,ifv) = &
                        kappa_lbl(inb,itg,ixw,ixc,ifv)/nb_deta
                     tau_lbl(inb,itg,ixw,ixc,ifv,:) = &
                        tau_lbl(inb,itg,ixw,ixc,ifv,:)/nb_deta
               
                     !Sort the spectral absorption coefficient data for
                     !the narrow band in increasing order
                     call sort_heap(k_array)

                     !--------------------------------------------------
                     !Compute g
                     !--------------------------------------------------
                     !Zero out the sum arrays
                     g_array = 0._dp
            
                     !Check where to add to the eta integration
                     do ilb=llb,ulb
                        i = ilb - llb + 1
                        g_array(i) = g_array(i) + lbl_deta(ilb)
                     enddo
                     g_array = g_array/nb_deta
                     
                     !Finish up the g calculation
                     do i=2,k_total
                        g_array(i) = min(g_array(i-1) + g_array(i),1._dp)
                     enddo

                     !--------------------------------------------------
                     !Invert k(g) dependence
                     !--------------------------------------------------
                     if (k_total.eq.1) then
                        k_out(i,itg,ixw,ixc,ifv,inb) = k_array(k_total)
                     else
                        do i=1,g_total
!write(*,*) k_total,lnb_array(inb),unb_array(inb),inb
!write(*,*) 'g',k_total,g_array(1),g_array(k_total),k_array(1),k_array(k_total)
                           call polint(g_array,k_array,g_out(i),&
                                       k_out(i,itg,ixw,ixc,ifv,inb),2)
                        enddo
                     endif
                     deallocate(g_array,k_array)
               
                     !--------------------------------------------------
                     !Compute the absorption coefficient and 
                     !the transmittance from the k-distributions
                     !--------------------------------------------------
                     kappa_ck(inb,itg,ixw,ixc,ifv) = 0._dp
                     tau_ck(inb,itg,ixw,ixc,ifv,:) = 0._dp
                     do i=1,g_total
                        !Absorption coefficient
                        kappa = k_out(i,itg,ixw,ixc,ifv,inb)
                           
                        !Quadrature weight
                        if (i.eq.1) then
                           dg = g_out(i)
                        elseif (i.eq.g_total) then
                           dg = 1._dp - g_out(i)
                        else
                           dg = 0.5_dp*(g_out(i+1) - g_out(i-1))
                        endif
                     
                        !Band-average absorption coefficient   
                        kappa_ck(inb,itg,ixw,ixc,ifv) = &
                           kappa_ck(inb,itg,ixw,ixc,ifv) + kappa*dg

                        !Band-average transmittance
                        do ipl=1,npl
                           tau_ck(inb,itg,ixw,ixc,ifv,ipl) = &
                              tau_ck(inb,itg,ixw,ixc,ifv,ipl) + &
                                 dexp(-kappa*plength(ipl))*dg
                        enddo
                     enddo
                     !--------------------------------------------------
                     !Compute errors
                     !--------------------------------------------------
                     ekappa(inb,itg,ixw,ixc,ifv) = &
                        dabs(kappa_ck(inb,itg,ixw,ixc,ifv) - &
                           kappa_lbl(inb,itg,ixw,ixc,ifv))/&
                              (kappa_lbl(inb,itg,ixw,ixc,ifv) + small)*100._dp
                     do ipl=1,npl
                        etau(inb,itg,ixw,ixc,ifv,ipl) = &
                           dabs(tau_ck(inb,itg,ixw,ixc,ifv,ipl) - &
                              tau_lbl(inb,itg,ixw,ixc,ifv,ipl))/&
                                 (tau_lbl(inb,itg,ixw,ixc,ifv,ipl) +&
                                    small)*100._dp
                     enddo
                  enddo narrow_band_loop
                  !$OMP END PARALLEL DO
               enddo temperature_loop
            enddo soot_loop
         enddo co2_loop
      enddo h2o_loop

      !-----------------------------------------------------------------
      !Dump data
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Dump data')
      g_unit = get_file_unit()
      open(file=g_file,unit=g_unit,action='write',form='unformatted')
      write(g_unit) g_total,nnb,ntg,nxc,nxw,nfv
      write(g_unit) lnb_array,unb_array
      write(g_unit) T_array
      write(g_unit) xs_array(:,id_co2)
      write(g_unit) xs_array(:,id_h2o)
      write(g_unit) xs_array(:,id_soot)
      write(g_unit) g_out
      write(g_unit) k_out
      close(g_unit)

      !-----------------------------------------------------------------
      !Dump errors
      !-----------------------------------------------------------------
      if (dump_error) then
         call dprint('generate_premix_nbck: Dump errors')
         error_unit = get_file_unit()
         open(file=error_file,unit=error_unit,action='write',&
              form='formatted')
         write(error_unit,'(100(a,:,","))') 'T','xc','xw','fv',&
            'lbound','ubound','kappa','kappa','kappa',&
            ('tau',ipl=1,3*npl)
         write(error_unit,'(100(a,:,","))') '','','','','','','LBL',&
            'CK','Error',('LBL','CK','Error',ipl=1,npl)   
         write(error_unit,'(9(a,:,","),100(e26.15e3,:,","))') &
            '','','','','','','','','',(plength(ipl),plength(ipl),&
            plength(ipl),ipl=1,npl)
         do itg=1,ntg
            do ixc=1,nxc
               do ixw=1,nxw
                  do ifv=1,nfv
                     do inb=1,nnb
!write(*,*) inb,kappa_lbl(inb,itg,ixw,ixc,ifv),kappa_ck(inb,itg,ixw,ixc,ifv)
                        write(error_unit,'(100(e26.15e3,:,","))') &
                           T_array(itg),xs_array(ixc,id_co2),&
                           xs_array(ixw,id_h2o),xs_array(ifv,id_soot),&
                           lnb_array(inb),unb_array(inb),&
                           kappa_lbl(inb,itg,ixw,ixc,ifv),&
                           kappa_ck(inb,itg,ixw,ixc,ifv),&
                           ekappa(inb,itg,ixw,ixc,ifv),&
                           (tau_lbl(inb,itg,ixw,ixc,ifv,ipl),&
                           tau_ck(inb,itg,ixw,ixc,ifv,ipl),&
                           etau(inb,itg,ixw,ixc,ifv,ipl),ipl=1,npl)
                     enddo
                  enddo
               enddo
            enddo
         enddo
         close(error_unit)
      endif

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Deallocate arrays')
      deallocate(ekappa,etau)
      deallocate(g_out,k_out)
      deallocate(kappa_ck,kappa_lbl)
      deallocate(tau_ck,tau_lbl)

      !-----------------------------------------------------------------
      !Compute total run time
      !-----------------------------------------------------------------
      call dprint('generate_premix_nbck: Compute total run time')
      end_time = get_wall_time()                                        !Stop time counter
      if (present(run_time)) run_time = end_time - start_time           !Compute run time, if requested

   endsubroutine generate_premix_nbck

   !====================================================================
   !Subroutine that generates data for the Planck-mean absorption 
   !coefficient from a LBL database and dumps to an external file
   !====================================================================
   subroutine generate_box_kappa(&
      T_array,xs_array,lnb_array,unb_array,id_spec,file_name,&
      pressure,run_time)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
         dprint,get_file_unit,get_wall_time,shutdown
      use math_functions, only: locate
      use lbl_parameters
      use omp_lib
      use physical_functions, only: Ib_function
      use precision_parameters, only: dp,small
      implicit none
      character(*),intent(in) :: file_name
      integer,intent(in) :: id_spec
      integer :: ounit,ierr
      integer :: ilb,inb,itg,its,ixs
      integer :: nlb,nnb,ntg,nts,nxs
      integer :: llb,ulb
      integer,allocatable,dimension(:,:) :: its_array
      real(dp),intent(in) :: lnb_array(:),T_array(:),unb_array(:),&
         xs_array(:)
      real(dp),optional,intent(in) :: pressure
      real(dp),optional,intent(out) :: run_time
      real(dp) :: denum,Ib,kappa,num,pres
      real(dp) :: start_time,end_time
      real(dp),allocatable,dimension(:,:,:) :: kp
      
      !-----------------------------------------------------------------
      !Set up initial parameters
      !-----------------------------------------------------------------
      call dprint('generate_box_kappa: Set up initial parameters')
      start_time = get_wall_time()                                      !Start time counter
      
      !Set up optional parameters
      pres = 1._dp; if (present(pressure)) pres = pressure              !Total pressure

      !-----------------------------------------------------------------
      !Load LBL data files
      !-----------------------------------------------------------------
      call dprint('generate_nbck: Load LBL data files')
      call load_lbl_data

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('generate_box_kappa: Define surrogate names')
      nnb = size(lnb_array)                                             !Total number of narrow bands
      call assert(nnb.eq.size(unb_array),&
                  'size(lnb_array)=size(unb_array)')                    !Check if the sizes of leta_in and ueta_in match
      ntg = size(T_array)                                               !Number of discrete temperature values
      nxs = size(xs_array)                                              !Number of discrete mole fraction values
                                                                        !  (same for all species)
      nlb = lbl_nlines                                                  !Total number of lines in the external LBL data files

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('generate_box_kappa: Allocate arrays')
      allocate(kp(ntg,nxs,nnb),stat=ierr)
      call CheckMemAlloc('kp',ierr)

      !-----------------------------------------------------------------
      !Mount iTS array
      !-----------------------------------------------------------------
      nts = ntg*nxs*nnb
      allocate(its_array(nts,3))
      its = 0
      do itg=1,ntg
         do ixs=1,nxs
            do inb=1,nnb
               its = its + 1
               its_array(its,1:3) = (/ itg, ixs, inb /)
            enddo
         enddo
      enddo
      
      !$OMP PARALLEL DO DEFAULT(SHARED) &
      !$OMP&   PRIVATE(itg,ixs,inb,its,llb,ulb,ilb,kappa,Ib,num,denum) 
      ts_loop: do its=1,nts
      
         !Get current indexes
         itg = its_array(its,1)
         ixs = its_array(its,2)
         inb = its_array(its,3)
         
         !Indexes of the LBL data corresponding to the 
         !positions of the  lower and upper bounds of the
         !current narrow band
         llb = locate(lbl_xeta,lnb_array(inb),nlb) + 1
         ulb = locate(lbl_xeta,unb_array(inb),nlb)

         num = 0._dp; denum = 0._dp
         do ilb=llb,ulb
            !Absorption coefficient
            kappa = lbl_kappa_func(xs_array(ixs),T_array(itg),&
                                   pres,id_spec,ilb,kappa_p=.true.)
               
            !Planck function
            Ib = Ib_function(T_array(itg),lbl_xeta(ilb))
   
            !Integrate
            num = num + kappa*Ib*lbl_deta(ilb)
            denum = denum + Ib*lbl_deta(ilb)
            
         enddo
         
         !Finish computing the Planck-mean absorption coefficient
         kp(itg,ixs,inb) = num/denum

      enddo ts_loop
      !$OMP END PARALLEL DO

      !-----------------------------------------------------------------
      !Dump data
      !-----------------------------------------------------------------
      call dprint('generate_box_kappa: Dump data')
      ounit = get_file_unit()
      open(file=file_name,unit=ounit,action='write',form='unformatted')
      write(ounit) nnb,ntg,nxs
      write(ounit) lnb_array,unb_array
      write(ounit) T_array
      write(ounit) xs_array
      write(ounit) kp
      close(ounit)

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('generate_box_kappa: Deallocate arrays')
      deallocate(kp)

      !-----------------------------------------------------------------
      !Compute total run time
      !-----------------------------------------------------------------
      end_time = get_wall_time()                                        !Stop time counter
      if (present(run_time)) run_time = end_time - start_time           !Compute run time, if requested

   endsubroutine generate_box_kappa

   !====================================================================
   !This subroutine converts formatted LBL data files into unformatted
   !(binary) files, which are quickly to read; the data files are also
   !properly post-processed before outputting
   !====================================================================
   subroutine convert_lbl_data
   
      use comp_functions, only: get_file_unit
      use global_parameters, only: number_of_species
      implicit none
      character(5),parameter :: new_ext = '.data'
      character(200) :: new_name,original_name
      integer :: ounit
      integer :: isp,ixs,nsp,nxs
      
      !-----------------------------------------------------------------
      !Load LBL data
      !-----------------------------------------------------------------
      call load_lbl_data
   
      !-----------------------------------------------------------------
      !Dump binary data
      !-----------------------------------------------------------------
      nsp = number_of_species
      do isp=1,nsp
         nxs = lbl_data_nxs(isp); if (nxs.le.0) cycle
         do ixs=1,nxs
            if (lbl_data_file(ixs,isp).eq.'null') cycle
            
            !Define the file name
            original_name = trim(lbl_data_prefix)//&
                            trim(lbl_data_file(ixs,isp))
            new_name = trim(original_name)//trim(new_ext)
            
            !Open unit
            ounit= get_file_unit()
            open(file=new_name,unit=ounit,form='unformatted',&
                 action='write')
            
            !Dump data
            write(ounit) lbl_nlines,.true.
            write(ounit) lbl_xeta
            write(ounit) lbl_deta
            write(ounit) lbl_data(:,:,ixs,isp)
            
            !Close unit
            close(ounit)
         enddo
      enddo
      
   endsubroutine convert_lbl_data
   
   !====================================================================
   !Subroutine to dump line properties to an external file name
   !====================================================================
   subroutine dump_lbl_properties(T,p,xs,file_name)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: dprint,get_file_unit
      use omp_lib
      implicit none
      character(*),intent(in) :: file_name
      integer :: ounit
      integer :: ilb,its,isp
      integer :: nlb,nts,nsp
      real(dp),intent(in) :: T(:),p(:),xs(:,:)
      real(dp),allocatable,dimension(:,:) :: kappa

      !-----------------------------------------------------------------
      !Load the LBL data
      !-----------------------------------------------------------------
      call dprint('dump_lbl_properties: Read the LBL data')
      call load_lbl_data

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('dump_lbl_properties: Surrogate names')
      nts = size(T)                                                     !Number of thermodynamic states
      nsp = size(xs,1)                                                  !Number of species
      nlb = lbl_nlines                                                  !Number of lines in the LBL data files

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      allocate(kappa(nts,nlb))
      
      !-----------------------------------------------------------------
      !Determine spectral absorption coefficient
      !-----------------------------------------------------------------
      call dprint('dump_lbl_properties: &
                  &Determine spectral absorption coefficient')
      !$OMP PARALLEL DEFAULT(SHARED)
      !$OMP DO
      lines_loop: do ilb=1,nlb

         !Check if this eta is to be solved for
         if ((lbl_xeta(ilb).lt.lbl_eta_min).or.&
             (lbl_xeta(ilb).gt.lbl_eta_max)) cycle lines_loop

         !Determine the properties for this wavenumber
         do its=1,nts
            do isp=1,nsp
               kappa(its,ilb) = kappa(its,ilb) + &
                  lbl_kappa_func(xs(isp,its),T(its),p(its),isp,ilb)
            enddo
         enddo
      enddo lines_loop
      !$OMP ENDDO
      !$OMP ENDPARALLEL
   
      !-----------------------------------------------------------------
      !Dump data
      !-----------------------------------------------------------------
      call dprint('dump_lbl_properties: Dump data')

      !Prepare unit
      ounit = get_file_unit()
      open(file=file_name,unit=ounit,action='write',form='formatted')

      !Dump data
      do ilb=1,nlb
         if ((lbl_xeta(ilb).lt.lbl_eta_min).or.&
             (lbl_xeta(ilb).gt.lbl_eta_max)) cycle
if (mod(ilb,4).ne.0) cycle
         write(ounit,'(100(e26.15e3,:,","))') &
            lbl_xeta(ilb),(kappa(its,ilb),its=1,nts)
      enddo
      
      !Close unit
      close(ounit)
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      deallocate(kappa)
      
   endsubroutine dump_lbl_properties
   
endmodule lbl_functions
