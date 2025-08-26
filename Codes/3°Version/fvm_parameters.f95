!#######################################################################
!Module that declares and sets up all parameters related to the
!Finite Volume Method for the solution of the radiation field
!#######################################################################
module fvm_parameters

   !====================================================================
   !Modules & Misc
   !====================================================================
   use mesh, only: one_d,two_d,three_d
   use precision_parameters, only: dp
   implicit none
   character(20) :: which_fvm_scheme,which_fvm_gamma
   integer :: fvm_angles,fvm_iterations
   real(dp) :: fvm_tolerance
   real(dp),allocatable,dimension(:) :: dlx,dly,dlz,domega
   
contains
   !====================================================================
   !Subroutine to correct the number of FVM angles
   !according to which_fvm_scheme
   !====================================================================
   subroutine correct_fvm_angles

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown,CheckMemAlloc
      use constants, only: pi,twopi
      implicit none
      integer :: itheta,ll,ntheta,nphi
      real(dp) :: theta_up,theta_low

      !-----------------------------------------------------------------
      !Correct the number of angles
      !-----------------------------------------------------------------
      if (trim(which_fvm_scheme).eq.'FT-FVM') then
         !1D: must be even
         if (one_d) then
            fvm_angles = ceiling(real(fvm_angles,dp)/2._dp)*2
         
         !2D: must be a multiple of four
         elseif (two_d) then
            fvm_angles = ceiling(real(fvm_angles,dp)/4._dp)*4
      
         !3D: ntheta(ntheta+1) correction (with ntheta even)
         else if (three_d) then
            ntheta = ceiling(sqrt(real(1+fvm_angles,dp)) - 1._dp)       !Preliminary ntheta calculation
            ntheta = ceiling(real(ntheta,dp)/2._dp)*2                   !Correct odd ntheta
            fvm_angles = ntheta*(ntheta + 2)                            !Update total number of anlges
         endif
      elseif (trim(which_fvm_scheme).eq.'FDS') then
         !1D: must be even
         if (one_d) then
            fvm_angles = ceiling(real(fvm_angles,dp)/2._dp)*2
         
         !2D: must be a multiple of four
         elseif (two_d) then
            fvm_angles = ceiling(real(fvm_angles,dp)/4._dp)*4
      
         !3D: use FDS's model
         else if (three_d) then
            ntheta = 2*nint(0.5_dp*1.17*&
                               real(fvm_angles)**(1._dp/2.26_dp))       !Number of theta divisions (even)
            ll = 0                                                      !Initializing the counter for the corrected total number of angles
            do itheta=1,ntheta
               theta_low = pi*real(itheta-1)/real(ntheta)
               theta_up  = pi*real(itheta)/real(ntheta)
               nphi = nint(0.5_dp*real(fvm_angles,dp)*&                 !Number of phi divisions for this
                                         (cos(theta_low)-cos(theta_up)))!  theta angle
               nphi = max(4,nphi)                                       !Minimum is 4
               nphi = 4*nint(0.25_dp*real(nphi,dp))                     !Rounding up to the nearest integer
               ll = ll + nphi                                           !Updating the counter                             
            enddo
            fvm_angles = ll                                             !Update total number of angles
         endif
      else
         call shutdown('correct_fvm_angles: which_fvm_scheme &
                       &not specified')
      endif
   
   endsubroutine correct_fvm_angles
   
   !====================================================================
   !Subroutine to allocate and initialize all FVM parameters
   !====================================================================
   subroutine initalize_fvm_parameters
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,shutdown
      implicit none
      integer :: ierr,nl
      
      !-----------------------------------------------------------------
      !First, correct the number of FVM angles
      !-----------------------------------------------------------------
      call correct_fvm_angles

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      nl = fvm_angles
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      if (allocated(domega)) deallocate(domega)
      if (allocated(dlx)) deallocate(dlx)
      if (allocated(dly)) deallocate(dly)
      if (allocated(dlz)) deallocate(dlz)
      
      allocate(domega(1:nl),stat=ierr)
      call CheckMemAlloc('domega',ierr)
      allocate(dlx(1:nl),stat=ierr)
      call CheckMemAlloc('dlx',ierr)
      allocate(dly(1:nl),stat=ierr)
      call CheckMemAlloc('dly',ierr)
      allocate(dlz(1:nl),stat=ierr)
      call CheckMemAlloc('dlz',ierr)
      
   endsubroutine initalize_fvm_parameters
   
   !====================================================================
   !Subroutine with the default parameters of the FVM
   !====================================================================
   !(this does not need to be included in the main 
   !code if the user knows what they are doing)
   subroutine set_default_fvm_parameters

      which_fvm_scheme = 'FT-FVM'
      which_fvm_gamma = 'step'
      fvm_angles = 400
      fvm_iterations = 1000
      fvm_tolerance = 1.e-6

   endsubroutine set_default_fvm_parameters

   !====================================================================
   !Subroutine to build the FVM angular grid
   !====================================================================
   subroutine build_fvm_angular_grid
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown,CheckMemAlloc
      use constants, only: pi,twopi
      implicit none
      integer :: iphi,itheta,ierr,ll,nphi,ntheta
      real(dp) :: dtheta,dphi,f_theta
      real(dp) :: theta_up,theta_low,phi_up,phi_low
      real(dp),allocatable,dimension(:,:) :: thetaf,phif

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      allocate(thetaf(1:2,1:fvm_angles),stat=ierr)
      call CheckMemAlloc('thetaf',ierr)
      allocate(phif(1:2,1:fvm_angles),stat=ierr)
      call CheckMemAlloc('phif',ierr)
   
      !-----------------------------------------------------------------
      !Define the theta and phi bounds
      !-----------------------------------------------------------------
      !1D: only theta varies
      if (one_d) then
         dtheta = pi/real(fvm_angles,dp)                                !(uniform) theta spacing
         do ll=1,fvm_angles
            thetaf(1,ll) = dtheta*real(ll-1,dp)                         !Lower theta bound
            thetaf(2,ll) = dtheta*real(ll,dp)                           !Upper theta bound
         enddo
         phif(1,:) = 0._dp                                              !Lower and upper phi bounds are
         phif(2,:) = twopi                                              !0 and 2pi for all control angles
            
      !2D: only phi varies
      elseif (two_d) then
         dphi = twopi/real(fvm_angles,dp)                               !(uniform) phi spacing
         do ll=1,fvm_angles
            phif(1,ll) = dphi*real(ll-1,dp)                             !Lower phi bound
            phif(2,ll) = dphi*real(ll,dp)                               !Upper phi bound
         enddo
         thetaf(1,:) = 0._dp                                            !Lower and upper theta bounds are
         thetaf(2,:) = pi                                               !0 and pi for all control angles

      !3D:
      else if (three_d) then
         if (trim(which_fvm_scheme).eq.'FT-FVM') then         
            ntheta = ceiling(sqrt(real(1+fvm_angles,dp)) - 1._dp)       !Number of theta divisions
            dtheta = pi/real(ntheta,dp)                                 !(uniform) theta spacing
            ll = 1                                                      !Initializing angle counter
            do itheta=1,ntheta
               if (itheta.le.ntheta/2) nphi = 4*itheta                  !Number of phi divisions 
               if (itheta.gt.ntheta/2) nphi = 4*(ntheta - itheta + 1)   !for this theta angle                      
               dphi = twopi/real(nphi,dp)                               !(uniform) phi spacing (for this theta)
               do iphi=1,nphi
                  thetaf(1,ll) = dtheta*real(itheta-1,dp)               !Lower theta bound
                  thetaf(2,ll) = dtheta*real(itheta,dp)                 !Upper theta bound             
                  phif(1,ll) = dphi*real(iphi-1,dp)                     !Lower phi bound
                  phif(2,ll) = dphi*real(iphi,dp)                       !Upper phi bound
                  ll = ll + 1                                           !Update angle counter
               enddo
            enddo
         elseif (trim(which_fvm_scheme).eq.'FDS') then
            ntheta = nint(1.17*real(fvm_angles)**(1._dp/2.26_dp))       !Number of theta divisions
            dtheta = pi/real(ntheta,dp)                                 !(uniform) theta spacing
            ll = 1                                                      !Initializing angle counter
            do itheta=1,ntheta                   
               theta_low = dtheta*real(itheta-1,dp)
               theta_up = dtheta*itheta 
               nphi = max(4,4*nint(0.25_dp*0.5_dp*real(fvm_angles,dp)*& !Number of phi divisions for this
                                (cos(theta_low)-cos(theta_up))))        !  theta angle               
               dphi = twopi/real(nphi,dp)                               !(uniform) phi spacing (for this theta)
               do iphi=1,nphi
                  thetaf(1,ll) = theta_low                              !Lower theta bound
                  thetaf(2,ll) = theta_up                               !Upper theta bound             
                  phif(1,ll) = dphi*real(iphi-1,dp)                     !Lower phi bound
                  phif(2,ll) = dphi*real(iphi,dp)                       !Upper phi bound
                  ll = ll + 1                                           !Update angle counter
               enddo
            enddo
         endif
      endif
      
      !-----------------------------------------------------------------
      !Compute the Dlm coefficients and angle sizes
      !-----------------------------------------------------------------
      dlx = 0._dp; dly = 0._dp; dlz = 0._dp
      if (one_d) then
         do ll=1,fvm_angles
            theta_low = thetaf(1,ll)
            theta_up  = thetaf(2,ll)
            dlx(ll) = pi*((sin(theta_up))**2._dp - &
                          (sin(theta_low))**2._dp)
            domega(ll) = 2._dp*pi*(cos(theta_low) - cos(theta_up))
         enddo 
         
      elseif (two_d.or.three_d) then
         do ll=1,fvm_angles
            theta_low = thetaf(1,ll); phi_low = phif(1,ll)
            theta_up  = thetaf(2,ll); phi_up  = phif(2,ll)
            f_theta = (theta_up - theta_low) - &
                      (sin(theta_up)*cos(theta_up) - &
                       sin(theta_low)*cos(theta_low))
            dlx(ll) = 0.5_dp*(sin(phi_up) - sin(phi_low))*f_theta
            dly(ll) = 0.5_dp*(cos(phi_up) - cos(phi_low))*f_theta
            if (three_d) dlz(ll) = 0.5_dp*(phi_up - phi_low)*&
                     ((sin(theta_up))**2._dp - (sin(theta_low))**2._dp)
            domega(ll) = (phi_low - phi_up)*&
               (cos(theta_up) - cos(theta_low))
         enddo
      endif

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      deallocate(thetaf)
      deallocate(phif)
      
   endsubroutine build_fvm_angular_grid

   !====================================================================
   !Subroutine to define the value of coefficient gamma 
   !(that relates the face intensities to the cell intensity) 
   !according to the which_fvm_gamma specification
   !====================================================================
   subroutine define_gamma(gval)

      use comp_functions, only: shutdown
      implicit none
      real(dp) :: gval
      
      selectcase(trim(which_fvm_gamma))
         case('diamond')
            gval = 0.5_dp
         case('step')
            gval = 1._dp
         case default
            call shutdown('define_gamma: which_fvm_gamma not specified')
      endselect
   
   endsubroutine define_gamma   
   
endmodule fvm_parameters
