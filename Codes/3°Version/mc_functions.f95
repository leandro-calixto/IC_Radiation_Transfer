!#######################################################################
!Module containing functions related to the Monte Carlo solution
!#######################################################################
module mc_functions

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp
   use mc_parameters
   implicit none
   
contains

   !====================================================================
   !Function to get a random number (selects according to 
   !which_mc_random from one of the random number generators in the 
   !math_functions module)
   !====================================================================
   real(dp) function get_mc_ran(idum)
      
      use comp_functions, only: shutdown
      use math_functions, only: get_numrep_random,get_random_uniform
      use precision_parameters, only: i4b
      implicit none
      integer(i4b), intent(inout) :: idum
      
   
      selectcase(trim(which_mc_random))
!      case('fortran')
!         get_mc_ran = get_random_uniform(0._dp,1._dp)
      case('numrep')
         get_mc_ran = get_numrep_random(idum)
      case default
         call shutdown('get_mc_ran: which_mc_random undefined')
      
      endselect
   
   endfunction get_mc_ran

   !====================================================================
   !Subroutine to read the external file containing the random number
   !relations (RNR) needed for the LBL/MC method
   !====================================================================
   subroutine prepare_lblrnr
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,CheckFileExists,&
                                get_file_unit
      use precision_parameters, only: big,small
      implicit none
      integer :: ict,iit,isp,ixs,ierr
      integer :: nln,nsp,ntg,nxs
      integer,allocatable,dimension(:) :: lblrnr_unit
      
      !-----------------------------------------------------------------
      !Preparing the input file
      !-----------------------------------------------------------------
      nsp = lblrnr_nsp                                                  !Number of species
      do iit=1,2                                                        !iit=1 -> Allocation; iit=2 -> Get and open unit
         ict = 0                                                        !Initialize unit counter
         do isp=1,nsp                                                   !Loop over all species
            nxs = lblrnr_nxs(isp); if (nxs.le.0) cycle                  !Skip species for which there is no RNR data file
            do ixs=1,nxs                                                !Loop over all discrete mole fraction values
               if (lblrnr_file(ixs,isp).eq.'null') cycle                !Skip mole fraction values for which there is no RNR data file
               ict = ict + 1                                            !Add to unit counter
               if (iit.eq.1) cycle                                      !Skip getting and opening the unit for allocation step
               lblrnr_unit(ict) = get_file_unit()                       !Get file unit
               call CheckFileExists(lblrnr_file(ixs,isp))               !Check if the file exists
               open(unit=lblrnr_unit(ict),file=lblrnr_file(ixs,isp),&   !Open the unit
                    form='unformatted',action='read')
               read(lblrnr_unit(ict)) lblrnr_ntg(ixs,isp),&             !Read first line
                                      lblrnr_nln(ixs,isp)
            enddo
         enddo
         if (iit.eq.1) allocate(lblrnr_unit(ict))                       !Allocate array with input units
      enddo 

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      !Array sizes at each dimension
      ntg = maxval(lblrnr_ntg)                                          !Temperature
      nln = maxval(lblrnr_nln)                                          !Spectral intervals
      nxs = maxval(lblrnr_nxs)                                          !Species mole fractions

      !Deallocate, allocate and check for memory errors
      if (allocated(lblrnr_eta)) deallocate(lblrnr_eta)
      allocate(lblrnr_eta(nln,nxs,nsp),stat=ierr)
      call CheckMemAlloc('lblrnr_eta',ierr)
      if (allocated(lblrnr_kpi)) deallocate(lblrnr_kpi)
      allocate(lblrnr_kpi(nln,ntg,nxs,nsp),stat=ierr)
      call CheckMemAlloc('lblrnr_kpi',ierr)
      if (allocated(lblrnr_kppi)) deallocate(lblrnr_kppi)
      allocate(lblrnr_kppi(ntg,nxs,nsp),stat=ierr)
      call CheckMemAlloc('lblrnr_kpplancki',ierr)
      if (allocated(lblrnr_Retai)) deallocate(lblrnr_Retai)
      allocate(lblrnr_Retai(nln,ntg,nxs,nsp),stat=ierr)
      call CheckMemAlloc('lblrnr_Retai',ierr)
      if (allocated(lblrnr_tg)) deallocate(lblrnr_tg)
      allocate(lblrnr_tg(ntg,nxs,nsp),stat=ierr)
      call CheckMemAlloc('lblrnr_tg',ierr)
      
      !-----------------------------------------------------------------
      !Read the data
      !-----------------------------------------------------------------
      ict = 0                                                           !Reset unit counter
      lblrnr_etamax = small; lblrnr_etamin = big                        !Initial values for the maximum and minimum 
                                                                        !   wavenumber values in the RNR data
      do isp=1,nsp                                                      !Loop over all species
         nxs = lblrnr_nxs(isp); if (nxs.le.0) cycle                     !Skip species for which there is no RNR data file
         do ixs=1,nxs                                                   !Loop over all discrete mole fraction values
            if (lblrnr_file(ixs,isp).eq.'null') cycle                   !Skip mole fraction values for which there is no RNR data file
            ict = ict + 1                                               !Add to unit counter
            read(lblrnr_unit(ict)) lblrnr_tg(:,ixs,isp)                 !Read array with the discrete temperature values
            read(lblrnr_unit(ict)) lblrnr_eta(:,ixs,isp)                !Read array with the discrete wavenumber values
            read(lblrnr_unit(ict)) lblrnr_Retai(:,:,ixs,isp)            !Read the random number array
            read(lblrnr_unit(ict)) lblrnr_kpi(:,:,ixs,isp)              !Read the pressure-based spectral absorption coefficient array
            read(lblrnr_unit(ict)) lblrnr_kppi(:,ixs,isp)               !Read the pressure-based Planck-mean absorption coefficient array
            lblrnr_etamax = max(lblrnr_etamax,&                         !Compute maximum wavenumber value
                                maxval(lblrnr_eta(:,ixs,isp)))
            lblrnr_etamin = min(lblrnr_etamin,&                         !Compute minimum wavenumber value
                                minval(lblrnr_eta(:,ixs,isp)))
            close(lblrnr_unit(ict))                                     !Close the unit
         enddo
      enddo

   endsubroutine prepare_lblrnr
   
   !====================================================================
   !Function to interpolate for the pressure-based Planck-mean 
   !absorption coefficient of a single species 
   !(to be used in the framework of the LBL/MC method)
   !====================================================================
   real(dp) function get_kpplancki(tmp,xsp,isp)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use math_functions, only: locate
      use precision_parameters, only: small
      implicit none
      integer,intent(in) :: isp
      integer :: itg,ixs,ntg,nxs
      integer :: ltg,lxs,utg,uxs
      real(dp),intent(in) :: tmp,xsp
      real(dp) :: Qtg,Qxs,Rtg,Rxs
      
      !-----------------------------------------------------------------
      !Interpolation
      !-----------------------------------------------------------------   
      !Upper and lower indexes for mole fraction
      nxs = lblrnr_nxs(isp)
      lxs = locate(lblrnr_xs(:,isp),xsp,nxs)
      lxs = max(1,min(nxs-1,lxs)); uxs = min(lxs+1,nxs)
      if (xsp.lt.lblrnr_xs(1,isp))    uxs = lxs                         !These lines prevent extrapolations
      if (xsp.gt.lblrnr_xs(nxs,isp))  lxs = uxs                         !  in the mole fraction space

      !Initial values for the interpolation on mole fraction
      Qxs = 0._dp
      Rxs = (xsp - lblrnr_xs(lxs,isp))/&
         (lblrnr_xs(uxs,isp) - lblrnr_xs(lxs,isp) + small)
      if (lxs.eq.uxs) Rxs = 0._dp

      do ixs=lxs,uxs
        !Upper and lower indexes for temperature
         ntg = lblrnr_ntg(ixs,isp)
         ltg = locate(lblrnr_tg(:,ixs,isp),tmp,ntg)
         ltg = max(1,min(ntg-1,ltg)); utg = min(ltg + 1,ntg)
         if (tmp.lt.lblrnr_tg(1,ixs,isp))    utg = ltg                  !These lines prevent extrapolations
         if (tmp.gt.lblrnr_tg(ntg,ixs,isp))  ltg = utg                  !  in the temperature space

         !Initial values for the interpolation on temperature
         Qtg = 0._dp
         Rtg = (tmp - lblrnr_tg(ltg,ixs,isp))/&
            (lblrnr_tg(utg,ixs,isp) - lblrnr_tg(ltg,ixs,isp) + small)
         if (ltg.eq.utg) Rtg = 0._dp
         
         do itg=ltg,utg
            !Interpolate on temperature
            Rtg = 1._dp - Rtg
            Qtg = Qtg + Rtg*lblrnr_kppi(itg,ixs,isp)
         enddo
         
         !Interpolate on mole fraction
         Rxs = 1._dp - Rxs
         Qxs = Qxs + Rxs*Qtg
      enddo

      !Final value
      get_kpplancki = Qxs
   
   endfunction get_kpplancki

   !====================================================================
   !Function to interpolate for the spectral random number (inum=0) or 
   !the pressure-based spectral absorption coefficient (inum=1) of a
   !single species (to be used in the framework of the LBL/MC method)
   !====================================================================
   real(dp) function get_R_or_ki(eta,tmp,xsp,isp,inum)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use math_functions, only: locate
      use precision_parameters, only: small
      implicit none
      integer,intent(in) :: isp,inum
      integer :: iln,itg,ixs,nln,ntg,nxs
      integer :: lln,ltg,lxs,uln,utg,uxs
      real(dp),intent(in) :: eta,tmp,xsp
      real(dp) :: Qln,Qtg,Qval,Qxs,Rln,Rtg,Rxs
      
      !-----------------------------------------------------------------
      !Avoid errors with the definition of inum
      !-----------------------------------------------------------------
      if ((inum.ne.0).and.(inum.ne.1)) &
         call shutdown('get_R_or_ki: inum defined incorrectly')
      
      !-----------------------------------------------------------------
      !Interpolation
      !-----------------------------------------------------------------   
      !Upper and lower indexes for mole fraction
      nxs = lblrnr_nxs(isp)
      lxs = locate(lblrnr_xs(:,isp),xsp,nxs)
      lxs = max(1,min(nxs-1,lxs)); uxs = min(lxs+1,nxs)
      if (xsp.lt.lblrnr_xs(1,isp))    uxs = lxs
      if (xsp.gt.lblrnr_xs(nxs,isp))  lxs = uxs
      
      !Initial values for the interpolation on mole fraction
      Qxs = 0._dp
      Rxs = (xsp - lblrnr_xs(lxs,isp))/&
         (lblrnr_xs(uxs,isp) - lblrnr_xs(lxs,isp) + small)
      if (lxs.eq.uxs) Rxs = 0._dp

      do ixs=lxs,uxs
         !Upper and lower indexes for temperature
         ntg = lblrnr_ntg(ixs,isp)
         ltg = locate(lblrnr_tg(:,ixs,isp),tmp,ntg)
         ltg = max(1,min(ntg-1,ltg)); utg = min(ltg + 1,ntg)
         if (tmp.lt.lblrnr_tg(1,ixs,isp))    utg = ltg
         if (tmp.gt.lblrnr_tg(ntg,ixs,isp))  ltg = utg
         
         !Initial values for the interpolation on temperature
         Qtg = 0._dp
         Rtg = (tmp - lblrnr_tg(ltg,ixs,isp))/&
            (lblrnr_tg(utg,ixs,isp) - lblrnr_tg(ltg,ixs,isp) + small)
         if (ltg.eq.utg) Rtg = 0._dp
         
         do itg=ltg,utg         
            !Upper and lower indexes for spectral interval
            nln = lblrnr_nln(ixs,isp)
            lln = locate(lblrnr_eta(:,ixs,isp),eta,nln)
            lln = max(1,min(nln-1,lln)); uln = min(lln + 1,nln) 
            if (eta.le.lblrnr_eta(1,ixs,isp))   uln = lln
            if (eta.ge.lblrnr_eta(nln,ixs,isp)) lln = uln

            !Initial values for the interpolation on temperature
            Qln = 0._dp
            Rln = (eta - lblrnr_eta(lln,ixs,isp))/&
               (lblrnr_eta(uln,ixs,isp) - lblrnr_eta(lln,ixs,isp) +&
                  small)
            if (inum.eq.1) Rln = 1._dp
            if (lln.eq.uln) Rln = 0._dp
!            if ((inum.eq.1).and.(Rln.le.0.5_dp)) Rln = 0._dp            !When interpolating for the spectral absorption coefficient,
!            if ((inum.eq.1).and.(Rln.gt.0.5_dp)) Rln = 1._dp            !  skip the interpolation in the wavenumber space 
            
            do iln=lln,uln
               !Interpolate on the spectral interval
               if (inum.eq.0) Qval = lblrnr_Retai(iln,itg,ixs,isp)
               if (inum.eq.1) Qval = lblrnr_kpi(iln,itg,ixs,isp)
!write(*,*) ixs,itg,iln,inum,Qval
               Rln = 1._dp - Rln; Qln = Qln + Rln*Qval
            enddo
            
            !Interpolate on temperature
            Rtg = 1._dp - Rtg; Qtg = Qtg + Rtg*Qln
         enddo
         
         !Interpolate on mole fraction
         Rxs = 1._dp - Rxs; Qxs = Qxs + Rxs*Qtg
      enddo
      get_R_or_ki = Qxs

   endfunction get_R_or_ki

   !====================================================================
   !Function to determine the mixture spectral random number
   !====================================================================
   real(dp) function get_Reta(eta,tmp,xs)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc
      use precision_parameters, only: small
      implicit none
      integer :: ierr
      integer :: isp,nsp
      real(dp),intent(in) :: eta,tmp,xs(:)
      real(dp) :: kpl_m,Reta
      real(dp),allocatable,dimension(:) :: kpl_i

      !-----------------------------------------------------------------
      !Compute the Planck-mean absorption coefficient 
      !for each single species and for the mixture
      !-----------------------------------------------------------------
      nsp = size(xs)
      allocate(kpl_i(nsp),stat=ierr)
      call CheckMemAlloc('kpl_i',ierr)
      kpl_m = 0._dp
      do isp=1,nsp
         kpl_i(isp) = xs(isp)*get_kpplancki(tmp,xs(isp),isp)
         kpl_m = kpl_m + kpl_i(isp)
      enddo

      !-----------------------------------------------------------------
      !Compute Reta
      !-----------------------------------------------------------------
      Reta = 0._dp
      do isp=1,nsp
         Reta = Reta + &
            kpl_i(isp)*get_R_or_ki(eta,tmp,xs(isp),isp,0)/(kpl_m+small)
      enddo
      get_Reta = Reta
   
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      deallocate(kpl_i)
   
   endfunction get_Reta

   !====================================================================
   !Function to determine the mixture spectral absorption coefficient
   !====================================================================
   real(dp) function get_mc_keta(eta,tmp,xs)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer :: isp,nsp
      real(dp),intent(in) :: eta,tmp,xs(:)
      real(dp) :: kmix
      
      !-----------------------------------------------------------------
      !Compute the spectral absorption coefficient
      !-----------------------------------------------------------------
      nsp = size(xs)                                                    !Number of species
      kmix = 0._dp                                                      !Initialize the variable for the sum
      do isp=1,nsp                                                      !Loop over species
         kmix = kmix + xs(isp)*get_R_or_ki(eta,tmp,xs(isp),isp,1)       !Sum over all species absorption coefficients
      enddo                                                             !  (assuming p = 1 atm)
      get_mc_keta = kmix
   
   endfunction get_mc_keta

   !====================================================================
   !Function to determine the mixture Planck mean-absorption coefficient
   !====================================================================
   real(dp) function get_mc_kplanck(tmp,xs)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer :: isp,nsp
      real(dp),intent(in) :: tmp,xs(:)
      real(dp) :: kmix
      
      !-----------------------------------------------------------------
      !Compute the spectral absorption coefficient
      !-----------------------------------------------------------------
      nsp = size(xs)                                                    !Number of species
      kmix = 0._dp                                                      !Initialize the variable for the sum
      do isp=1,nsp                                                      !Loop over species
         kmix = kmix + xs(isp)*get_kpplancki(tmp,xs(isp),isp)           !Sum over all species Planck-mean absorption
      enddo                                                             !   coefficients (assuming p = 1 atm)
      get_mc_kplanck = kmix
      
   endfunction get_mc_kplanck
   
   !====================================================================
   !Function to determine the correct wavenumber from a given spectral
   !random number in the framework of the LBL/MC method
   !====================================================================
   real(dp) function get_mc_wavenumber(Rtarget,tmp,xs,wall)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,shutdown
      use physical_functions, only: bb_emission_frac
      use precision_parameters, only: small
      implicit none
      integer :: counter
      real(dp),intent(in) :: Rtarget,tmp,xs(:)
      real(dp) :: xleft,xright,xguess
      real(dp) :: yleft,yright,yguess
      real(dp) :: tol,yfx
      logical,intent(in),optional :: wall
      logical :: is_wall
      is_wall = .false.; if (present(wall)) is_wall = .true.            !Set up the optional parameter
      
      !-----------------------------------------------------------------
      !Initial boundary values and guess
      !-----------------------------------------------------------------
      if (is_wall) then                                                 
         xleft = 99999999999._dp; xright = small                        !Boundaries in x
         yleft = bb_emission_frac(xleft,tmp,eta_in=.true.) - Rtarget    !Boundaries in y
         yright = bb_emission_frac(xright,tmp,eta_in=.true.) - Rtarget
      else
         xleft = lblrnr_etamin; xright = lblrnr_etamax                  !Boundaries in x
         yleft = get_Reta(xleft,tmp,xs) - Rtarget                       !Boundaries in y
         yright = get_Reta(xright,tmp,xs) - Rtarget               
      endif
      xguess = xleft + (xright - xleft)*Rtarget                         !Initial guess

      !-----------------------------------------------------------------
      !This serves to avoid entering the iteration process if the target
      !is outside the bounds of the data (this should never occur for
      !properly constructed data)
      !-----------------------------------------------------------------
      if (yleft.gt.0) then
         get_mc_wavenumber = xleft
         return
      elseif (yright.lt.0) then
         get_mc_wavenumber = xright
         return
      endif      

      !-----------------------------------------------------------------
      !Root finding
      !-----------------------------------------------------------------
      counter = 0                                                       !Initialize iteration counter
      tol = lblrnr_tolerance                                            !Tolerance for the iterative process
      bisection_loop: do
         
         !Compute y (tentative root)
         if (is_wall) yfx = bb_emission_frac(xguess,tmp,eta_in=.true.)
         if (.not.is_wall) yfx = get_Reta(xguess,tmp,xs)
         yguess = yfx - Rtarget

         !Check for convergence
         if (dabs(yguess).lt.tol) exit

         !Update boundary values
         if (yguess*yleft.lt.0) then
            xright = xguess; yright = yguess
         elseif (yguess*yright.lt.0) then
            xleft = xguess; yleft = yguess
         endif
         
         !Update x guess
         xguess = 0.5_dp*(xleft + xright)
      
         !Update counter
         counter = counter + 1
         if (counter.gt.lblrnr_max_iterations) &
            call shutdown('get_mc_wavenumber: maximum number of &
                          &iterations exceeded')
      
      enddo bisection_loop
      get_mc_wavenumber = xguess

   endfunction get_mc_wavenumber
   
   !====================================================================
   !Function to determine the correct gray gas from a given spectral
   !random number in the framework of the WSGG/MC method
   !====================================================================
   integer function get_mc_WSGGj(Rtarget,tmp,tp,xs,wall)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use wsgg_parameters, only: which_wsgg_model
      use wsgg_functions, only: wsgg_kappa_planck,&
         wsgg_number_gray_gases,wsgg_properties
      implicit none
      integer :: j,ngg
      real(dp),intent(in) :: Rtarget,tmp,tp,xs(:)
      real(dp) :: aj,kj,kplanck,rsum
      logical,intent(in),optional :: wall
      logical :: is_wall
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      is_wall = .false.; if (present(wall)) is_wall = wall              !Set up optional parameter      
      ngg = wsgg_number_gray_gases(which_wsgg_model)                    !Number of gray gases
      
      !-----------------------------------------------------------------
      !Main calculation
      !-----------------------------------------------------------------
      if (is_wall) then
         !For a wall cell, the gray gas selection is made on the basis
         !of the weighting coefficient only
         rsum = 0._dp                                                   !Initialize the sum array
         do j=1,ngg+1                                                   !Loop over all gray gases and transparent windows
            call wsgg_properties(kj,aj,j,tmp,xs,tp,which_wsgg_model)    !Weighting coefficient
            rsum = rsum + aj                                            !Add to the weghting coefficient sum
            if (rsum.gt.Rtarget) exit
         enddo      
      else
         !Within the medium, the gray gas selection is made on the basis
         !of the fraction of radiative emission
         kplanck = wsgg_kappa_planck(tmp,tp,xs)                         !Planck-mean absorption coefficient
         rsum = 0._dp                                                   !Initialize the sum array
         do j=1,ngg                                                     !Loop over gray gases only
            call wsgg_properties(kj,aj,j,tmp,xs,tp,which_wsgg_model)    !Absorption and weighting coefficients
            rsum = rsum + kj*aj/kplanck                                 !Add to the cumulative fraction of radiative emission
            if (rsum.gt.Rtarget) exit
         enddo
      endif         
      get_mc_WSGGj = j                                                  !Final value
      
   endfunction get_mc_WSGGj
      
endmodule mc_functions
