!#######################################################################
!Module containing routines related to the exact RTE solution
!#######################################################################
module exact_solution

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp,small
   implicit none
   character(40) :: exact_interpolation
   real(dp) :: exact_tol
   
contains

   !====================================================================
   !Subroutine to compute the exact RTE solution along
   !a line-of-sight for a non-scattering medium
   !====================================================================
   subroutine exact_los_nonsct(x,kappa,Ib,I0,nx,Irad,simplified)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc
      use math_functions, only: linint,splint,spline
      use precision_parameters, only: big
      implicit none
      integer,intent(in) :: nx
      integer :: counter,i,ierr,intmin
      real(dp),intent(in) :: kappa(nx),Ib(nx),I0,x(nx)
      real(dp),intent(out) :: Irad(1:nx)
      real(dp) :: diff,ds,dtau,s,sint,sint0,tau,taul,tol
      real(dp) :: kappa_s,S_tau
      real(dp),allocatable,dimension(:) :: k2array,Ib2array
      real(dp),allocatable,dimension(:) :: emi,tau_cell
      logical,optional :: simplified
      logical :: assume_iso,lin_int,spl_int
      
      !-----------------------------------------------------------------
      !Setup parameters
      !-----------------------------------------------------------------
      assume_iso = .false.
      if (present(simplified)) assume_iso = simplified
      
      lin_int = .false.; spl_int = .false.
      if (trim(exact_interpolation).eq.'linear') lin_int = .true.       !Set flag for linear interpolation
      if (trim(exact_interpolation).eq.'spline') spl_int = .true.       !Set flag for spline interpolation
      tol = exact_tol                                                   !Tolerance for the integration
      intmin = 3                                                        !Minimum number of iterations for the
                                                                        !  numerical integrations

      if (assume_iso) then
         !--------------------------------------------------------------
         !Compute LOS intensity following a simplified, faster approach
         !--------------------------------------------------------------
         Irad(1) = I0
         los_loop: do i=2,nx
            tau = dexp(-kappa(i)*(x(i) - x(i-1)))
            Irad(i) = Irad(i-1)*tau + Ib(i)*(1._dp - tau)
         enddo los_loop
      else
         !--------------------------------------------------------------
         !Allocate arrays
         !--------------------------------------------------------------
         allocate(emi(1:nx),stat=ierr)
         call CheckMemAlloc('emi',ierr)
         allocate(k2array(1:nx),stat=ierr)
         call CheckMemAlloc('k2array',ierr)
         allocate(Ib2array(1:nx),stat=ierr)
         call CheckMemAlloc('Ib2array',ierr)
         allocate(tau_cell(1:nx),stat=ierr)
         call CheckMemAlloc('tau_cell',ierr)
   
        !--------------------------------------------------------------
        !Compute cell-centered optical thickness
        !--------------------------------------------------------------
        if (spl_int) call spline(x,kappa,nx,big,big,k2array)           !Prepare spline for the subsequent interpolation
        tau_cell(1) = 0._dp                                            !Optical tickness at the first point is zero
        los_loop_1: do i=2,nx         
           diff = 2._dp*tol; sint0 = 0._dp; counter = 0                !Initial values to start the integration loop
           ds = (x(i) - x(i-1))/4._dp
           do while (diff.gt.tol)
              counter = counter + 1
              s = x(i-1)
              sint = sint0/(2._dp*ds)
              do while(s.le.x(i))
                 if (lin_int) kappa_s = linint(x,kappa,s,nx)
                 if (spl_int) kappa_s = splint(x,kappa,k2array,s,nx)
                 sint = sint + kappa_s
                 s = s + ds
                 if (sint0.ne.0._dp) s = s + ds
              enddo
              sint = sint*ds
              if (counter.gt.intmin) &
                 diff = abs((sint - sint0)/(sint0 + small))*100._dp
              ds = ds/2._dp
              sint0 = sint
           enddo
           tau_cell(i) = tau_cell(i-1) + sint
        enddo los_loop_1
           
        if (spl_int) call spline(tau_cell,Ib,nx,big,big,Ib2array)      !Prepare spline for the subsequent interpolation
        taul = tau_cell(1)
        Irad(1) = I0                                                   !Intensity at the first grid point is equal 
                                                                       !  to the intensity at the boundary
        los_loop_2: do i=2,nx
            !-----------------------------------------------------------
            !Compute source term integral
            !-----------------------------------------------------------
            diff = 2._dp*tol; sint0 = 0._dp; counter = 0     
            dtau = (tau_cell(i) - tau_cell(i-1))/4._dp
            do while (diff.gt.tol)
               counter = counter + 1
               taul = tau_cell(i-1)
               sint = sint0/(2._dp*dtau)
               do while(taul.le.tau_cell(i))
                  if (lin_int) S_tau = linint(tau_cell,Ib,taul,nx)
                  if (spl_int) S_tau = &
                     splint(tau_cell,Ib,Ib2array,taul,nx)
                  sint = sint + S_tau*dexp(taul - tau_cell(i))
                  taul = taul + dtau
                  if (sint0.ne.0._dp) taul = taul + dtau
               enddo
               sint = sint*dtau
               if (counter.gt.intmin) &
                  diff = abs((sint - sint0)/(sint0 + small))*100._dp
               dtau = dtau/2._dp
               sint0 = sint
            enddo
           
            !-----------------------------------------------------------
            !Compute the cell-centered intensity
            !-----------------------------------------------------------
            Irad(i) = Irad(i-1)*dexp(-(tau_cell(i)-(tau_cell(i-1)))) + &
                      sint

         enddo los_loop_2
         
         !--------------------------------------------------------------
         !Deallocate arrays
         !--------------------------------------------------------------
         deallocate(emi)
         deallocate(k2array)
         deallocate(Ib2array)
         deallocate(tau_cell)
      
      endif
   
   endsubroutine exact_los_nonsct
      
   !====================================================================
   !Subroutine to compute the exact RTE solution for a  
   !plane-parallel, non-scattering, gray, 1d medium
   !====================================================================
   subroutine exact_1d_pp(x,kappa,Ib,nx,source,flux,absorption,emission)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: fourpi,twopi
      use comp_functions, only: CheckMemAlloc
      use math_functions, only: expint,spline,splint,polint
      use precision_parameters, only: big
      implicit none
      integer,intent(in) :: nx
      integer :: i,ierr
      logical :: transparent_medium
      real(dp),intent(in) :: x(nx),kappa(nx),Ib(nx)
      real(dp),intent(out) :: flux(nx),source(nx)
      real(dp),intent(out),optional :: absorption(nx),emission(nx)
      real(dp) :: diff,ds,dtau,int0,s,taul,tol
      real(dp) :: expnum,kappa_s,Ib_tau,tau_inc,tau_L
      real(dp),allocatable,dimension(:) :: radabs,rademi
      real(dp),allocatable,dimension(:) :: k2array,Ib2array,tau_cell
      real(dp),allocatable,dimension(:) :: IbE1_b,IbE1_f,IbE2_b,IbE2_f
      
      !-----------------------------------------------------------------
      !Set parameters
      !-----------------------------------------------------------------
      tol = exact_tol                                                   !Tolerance for the integration
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      allocate(k2array(2:nx-1),stat=ierr)
      call CheckMemAlloc('k2array',ierr)
      allocate(Ib2array(2:nx-1),stat=ierr)
      call CheckMemAlloc('Ib2array',ierr)
      allocate(IbE1_b(1:nx),stat=ierr)
      call CheckMemAlloc('IbE1_b',ierr)
      allocate(IbE1_f(1:nx),stat=ierr)
      call CheckMemAlloc('IbE1_f',ierr)
      allocate(IbE2_b(1:nx),stat=ierr)
      call CheckMemAlloc('IbE2_b',ierr)
      allocate(IbE2_f(1:nx),stat=ierr)
      call CheckMemAlloc('IbE2_f',ierr)
      allocate(radabs(1:nx),stat=ierr)
      call CheckMemAlloc('radabs',ierr)
      allocate(rademi(1:nx),stat=ierr)
      call CheckMemAlloc('rademi',ierr)
      allocate(tau_cell(1:nx),stat=ierr)
      call CheckMemAlloc('tau_cell',ierr)

      !-----------------------------------------------------------------
      !Special case for fully transparent media
      !-----------------------------------------------------------------
      transparent_medium = .true.
      do i=2,nx-1
         if (kappa(i).gt.small) then
            transparent_medium = .false.
            exit
         endif
      enddo

      if (transparent_medium) then
         tau_cell = 0._dp
         IbE1_f = 0._dp; IbE1_b = 0._dp
         IbE2_f = 1._dp; IbE2_b = 1._dp      
         
      else

      !-----------------------------------------------------------------
      !Compute cell-centered optical thickness
      !-----------------------------------------------------------------
      call spline(x(2:nx-1),kappa(2:nx-1),nx-2,big,big,k2array(2:nx-1)) !Prepare spline for the subsequent interpolation
      tau_cell(1) = 0._dp                                               !Optical tickness at the first point is zero
      tau_loop: do i=2,nx         
         diff = 10._dp*tol                                              !Initial difference between two successive iterations
         ds = (x(i) - x(i-1))/10._dp                                    !Initial finite increment for the integration
         int0 = 0._dp; tau_inc = 0._dp                                  !Initial "old" integration value 
         do while (diff.gt.tol)
            s = x(i-1)                                                  !Initial optical path coordinate value (previous grid point position)
            tau_inc = 0._dp                                             !kappa*ds integration is only performed between the previous and the
                                                                        !  current grid point, i.e., we are computing the increment in tau
                                                                        !  between these two points
            do while(s.le.x(i))
               if (i.eq.2) then
                  kappa_s = kappa(i)
               elseif (i.eq.nx) then
                  kappa_s = kappa(nx-1)
               else
                  kappa_s = splint(x(2:nx-1),kappa(2:nx-1),&
                                   k2array(2:nx-1),s,nx-2)              !Local absorption coefficient (interpolated)
               endif
               tau_inc = tau_inc + kappa_s*ds                           !Add to the kappa*ds integration
               s = s + ds                                               !Update the optical path coordinate value
            enddo
            diff = dabs((tau_inc - int0)/(int0 + small))*100._dp        !Difference relative to the previous integration
            ds = ds/2._dp                                               !Reduce integration increment for the next iteration step
            int0 = tau_inc                                              !Update "old" integration value
         enddo
         tau_cell(i) = tau_cell(i-1) + tau_inc
      enddo tau_loop

      tau_L = tau_cell(nx)                                              !Optical thickness at the right-hand boundary

      !-----------------------------------------------------------------
      !Compute the IbE1 and IbE2 integrals
      !-----------------------------------------------------------------
      call spline(tau_cell(2:nx-1),Ib(2:nx-1),nx-2,big,big,&
                  Ib2array(2:nx-1))                                     !Prepare spline for the subsequent interpolation

      IbE1_f(1) = 0._dp                                                 !Set the "forward" integral values at the left-hand
      IbE2_f(1) = 0._dp                                                 !  boundary to zero
      IbE_loop_1: do i=2,nx 

         !int_{0}^{tau_{i}} Ib(tau')*E1(tau-tau') dtau'
         diff = 10._dp*tol                                              !Initial difference between two successive iterations
         dtau = (tau_cell(i) - tau_cell(1))/10._dp                      !Initial finite increment for the integration
         int0 = 0._dp                                                   !Initial previous integration value 
         do while (diff.gt.tol)
            taul = tau_cell(1)                                          !Integration starts at tau' = 0
            IbE1_f(i) = 0._dp                                           !Initialize the sum variable to 0
            do while(taul.lt.tau_cell(i))
!               if (taul.lt.tau_cell(2)) then
!                  Ib_tau = Ib(2)
!               elseif (taul.gt.tau_cell(nx-1)) then
!                  Ib_tau = Ib(nx-1)
!               else
                  Ib_tau = splint(tau_cell(2:nx-1),Ib(2:nx-1),&         !Local blackbody intensity (interpolated)
                                  Ib2array(2:nx-1),taul,nx-2)
!               endif
               expnum = tau_cell(i) - taul                              !Exponent of the E_n function (should always be greater than 0)
               IbE1_f(i) = IbE1_f(i) + Ib_tau*expint(1,expnum)*dtau     !Add to the integration
               taul = taul + dtau                                       !Update tau'
            enddo
            diff = dabs((IbE1_f(i) - int0)/(int0 + small))*100._dp      !Difference relative to the previous integration
            dtau = dtau/2._dp                                           !Reduce tau' integration step for next iteration
            int0 = IbE1_f(i)                                            !Update "old" integration value
         enddo

         !int_{0}^{tau} Ib(tau')*E2(tau-tau') dtau'
         diff = 10._dp*tol
         dtau = (tau_cell(i) - tau_cell(1))/10._dp
         int0 = 0._dp
         do while (diff.gt.tol)
            taul = tau_cell(1)
            IbE2_f(i) = 0._dp
            do while(taul.lt.tau_cell(i))
!               if (taul.lt.tau_cell(2)) then
!                  Ib_tau = Ib(2)
!               elseif (taul.gt.tau_cell(nx-1)) then
!                  Ib_tau = Ib(nx-1)
!               else
                  Ib_tau = splint(tau_cell(2:nx-1),Ib(2:nx-1),&
                                  Ib2array(2:nx-1),taul,nx-2)
!               endif
               expnum = max(tau_cell(i) - taul,small)
               IbE2_f(i) = IbE2_f(i) + Ib_tau*expint(2,expnum)*dtau
               taul = taul + dtau
            enddo
            diff = dabs((IbE2_f(i) - int0)/(int0 + small))*100._dp
            dtau = dtau/2._dp
            int0 = IbE2_f(i)
         enddo
      enddo IbE_loop_1

      IbE1_b(nx) = 0._dp
      IbE2_b(nx) = 0._dp
      IbE_loop_2: do i=nx-1,1,-1

         !int_{tau}^{tau_L} Ib(tau')*E1(tau'-tau) dtau'
         diff = 10._dp*tol
         dtau = (tau_cell(nx) - tau_cell(i))/10._dp
         int0 = 0._dp
         do while (diff.gt.tol)
            taul = tau_cell(nx)
            IbE1_b(i) = 0._dp
            do while (taul.gt.tau_cell(i))
!               if (taul.lt.tau_cell(2)) then
!                  Ib_tau = Ib(2)
!               elseif (taul.gt.tau_cell(nx-1)) then
!                  Ib_tau = Ib(nx-1)
!               else
                  Ib_tau = splint(tau_cell(2:nx-1),Ib(2:nx-1),&
                                  Ib2array(2:nx-1),taul,nx-2)
!               endif
               expnum = max(taul - tau_cell(i),small)
               IbE1_b(i) = IbE1_b(i) + Ib_tau*expint(1,expnum)*dtau
               taul = taul - dtau
            enddo
            diff = dabs((IbE1_b(i) - int0)/(int0 + small))*100._dp
            dtau = dtau/2._dp
            int0 = IbE1_b(i)
         enddo                        

         !int_{tau}^{tau_L} Ib(tau')*E2(tau'-tau) dtau'
         diff = 10._dp*tol
         dtau = (tau_cell(nx) - tau_cell(i))/10._dp
         int0 = 0._dp
         do while (diff.gt.tol)
            taul = tau_cell(nx)
            IbE2_b(i) = 0._dp
            do while (taul.gt.tau_cell(i))
!               if (taul.lt.tau_cell(2)) then
!                  Ib_tau = Ib(2)
!               elseif (taul.gt.tau_cell(nx-1)) then
!                  Ib_tau = Ib(nx-1)
!               else
                  Ib_tau = splint(tau_cell(2:nx-1),Ib(2:nx-1),&
                                  Ib2array(2:nx-1),taul,nx-2)
!               endif
               expnum = max(taul - tau_cell(i),small)
               IbE2_b(i) = IbE2_b(i) + Ib_tau*expint(2,expnum)*dtau
               taul = taul - dtau
            enddo
            diff = abs((IbE2_b(i) - int0)/(int0 + small))*100._dp
            dtau = dtau/2._dp
            int0 = IbE2_b(i)
         enddo
      enddo IbE_loop_2

      endif

      !-----------------------------------------------------------------
      !Compute radiative absorption, emission, source and flux
      !-----------------------------------------------------------------
      final_loop: do i=1,nx
         rademi(i) = fourpi*kappa(i)*Ib(i)                              !Emission
         radabs(i) = twopi*(Ib(1)*expint(2,tau_cell(i)) + &             !Absorption
                            Ib(nx)*expint(2,tau_L - tau_cell(i)) + &
                            IbE1_f(i) + IbE1_b(i) )*kappa(i)
                            
         source(i) = radabs(i) - rademi(i)                              !Source
         flux(i) = twopi*(Ib(1)*expint(3,tau_cell(i)) - &               !Flux
                            Ib(nx)*expint(3,tau_L - tau_cell(i)) + &
                            IbE2_f(i) - IbE2_b(i) )
      enddo final_loop
      if (present(absorption)) absorption = radabs
      if (present(emission)) emission = rademi
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      deallocate(k2array,Ib2array)
      deallocate(IbE1_b,IbE1_f,IbE2_b,IbE2_f)
      deallocate(radabs,rademi)
      deallocate(tau_cell)
      
   endsubroutine exact_1d_pp



   subroutine exact_1d_pp_old(x,kappa,Ib,nx,source,flux,absorption,emission)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: fourpi,twopi
      use comp_functions, only: CheckMemAlloc
      use math_functions, only: expint,spline,splint
      use precision_parameters, only: big
      implicit none
      integer,intent(in) :: nx
      integer :: i,ierr
      real(dp),intent(in) :: x(nx),kappa(nx),Ib(nx)
      real(dp),intent(out) :: flux(nx),source(nx)
      real(dp),intent(out),optional :: absorption(nx),emission(nx)
      real(dp) :: diff,ds,dtau,int0,s,taul,tol
      real(dp) :: expnum,kappa_s,Ib_tau,tau_L
      real(dp),allocatable,dimension(:) :: radabs,rademi
      real(dp),allocatable,dimension(:) :: k2array,Ib2array,tau_cell
      real(dp),allocatable,dimension(:) :: IbE1_b,IbE1_f,IbE2_b,IbE2_f
      
      !-----------------------------------------------------------------
      !Set parameters
      !-----------------------------------------------------------------
      tol = exact_tol                                                   !Tolerance for the integration
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      allocate(k2array(1:nx),stat=ierr)
      call CheckMemAlloc('k2array',ierr)
      allocate(Ib2array(1:nx),stat=ierr)
      call CheckMemAlloc('Ib2array',ierr)
      allocate(IbE1_b(1:nx),stat=ierr)
      call CheckMemAlloc('IbE1_b',ierr)
      allocate(IbE1_f(1:nx),stat=ierr)
      call CheckMemAlloc('IbE1_f',ierr)
      allocate(IbE2_b(1:nx),stat=ierr)
      call CheckMemAlloc('IbE2_b',ierr)
      allocate(IbE2_f(1:nx),stat=ierr)
      call CheckMemAlloc('IbE2_f',ierr)
      allocate(radabs(1:nx),stat=ierr)
      call CheckMemAlloc('radabs',ierr)
      allocate(rademi(1:nx),stat=ierr)
      call CheckMemAlloc('rademi',ierr)
      allocate(tau_cell(1:nx),stat=ierr)
      call CheckMemAlloc('tau_cell',ierr)
      
      !-----------------------------------------------------------------
      !Compute cell-centered optical thickness
      !-----------------------------------------------------------------
      call spline(x,kappa,nx,big,big,k2array)                           !Prepare spline for the subsequent interpolation
      tau_cell(1) = 0._dp                                               !Optical tickness at the first point is zero
      tau_loop: do i=2,nx         

         diff = 1000._dp                                                !Initial difference between two successive iterations
         ds = (x(i) - x(i-1))/10._dp                                    !Initial finite increment for the integration
         int0 = tau_cell(i-1)                                           !Initial "old" integration value 
         do while (diff.gt.tol)
            s = x(i-1)                                                  !Initial optical path coordinate value (previous grid point position)
            tau_cell(i) = tau_cell(i-1)                                 !kappa*ds integration is only performed between the previous and the
                                                                        !  current grid point, so always add to the integration the optical
                                                                        !  thickness computed up to the previous grid point
            do while(s.le.x(i))
               kappa_s = splint(x,kappa,k2array,s,nx)                   !Local absorption coefficient (interpolated)
               tau_cell(i) = tau_cell(i) + kappa_s*ds                   !Add to the kappa*ds integration
               s = s + ds                                               !Update the optical path coordinate value
            enddo
            diff = abs((tau_cell(i) - int0)/(int0 + small))*100._dp     !Difference relative to the previous integration
!write(*,*) i,tau_cell(i),int0,diff
            ds = ds/2._dp                                               !Reduce integration increment for the next iteration step
            int0 = tau_cell(i)                                          !Update "old" integration value
         enddo

      enddo tau_loop
      tau_L = tau_cell(nx)                                              !Optical thickness at the right-hand boundary
      
      !-----------------------------------------------------------------
      !Compute the IbE1 and IbE2 integrals
      !-----------------------------------------------------------------
      call spline(x,Ib,nx,big,big,Ib2array)                             !Prepare spline for the subsequent interpolation
      IbE1_f(1) = 0._dp                                                 !Set the "forward" integral values at the left-hand
      IbE2_f(1) = 0._dp                                                 !  boundary to zero
      IbE_loop_1: do i=2,nx 

         !int_{0}^{tau} Ib(tau')*E1(tau-tau') dtau'
         diff = 1000._dp                                                !Initial difference between two successive iterations
         dtau = (tau_cell(i) - tau_cell(1))/10._dp                      !Initial finite increment for the integration
         int0 = 0._dp                                                   !Initial previous integration value 
         do while (diff.gt.tol)
!            if (dtau.le.small) exit
            taul = 0._dp                                                !Integration starts at tau' = 0
            IbE1_f(i) = 0._dp                                           !Initialize the sum variable to 0
            do while(taul.le.tau_cell(i))
               Ib_tau = splint(tau_cell,Ib,Ib2array,taul,nx)            !Local blackbody intensity (interpolated)
               expnum = max(tau_cell(i) - taul,small)                   !Exponent of the E_n function (should always be greater than 0)
               IbE1_f(i) = IbE1_f(i) + Ib_tau*expint(1,expnum)*dtau     !Add to the integration
               taul = taul + dtau                                       !Update tau'
write(*,*) i,taul,tau_cell(i),Ib_tau
            enddo
            diff = abs((IbE1_f(i) - int0)/(int0 + small))*100._dp       !Difference relative to the previous integration
            dtau = dtau/2._dp                                           !Reduce tau' integration step for next iteration
            int0 = IbE1_f(i)                                            !Update "old" integration value
         enddo

         !int_{0}^{tau} Ib(tau')*E2(tau-tau') dtau'
         diff = 1000._dp
         dtau = (tau_cell(i) - tau_cell(1))/10._dp
         int0 = 0._dp
         do while (diff.gt.tol)
!            if (dtau.le.small) exit
            taul = 0._dp
            IbE2_f(i) = 0._dp
            do while(taul.le.tau_cell(i))
               Ib_tau = splint(tau_cell,Ib,Ib2array,taul,nx)
               expnum = max(tau_cell(i) - taul,small)
               IbE2_f(i) = IbE2_f(i) + Ib_tau*expint(2,expnum)*dtau
               taul = taul + dtau
            enddo
            diff = abs((IbE2_f(i) - int0)/(int0 + small))*100._dp
            dtau = dtau/2._dp
            int0 = IbE2_f(i)
         enddo
      enddo IbE_loop_1
      
      IbE1_b(nx) = 0._dp
      IbE2_b(nx) = 0._dp
      IbE_loop_2: do i=1,nx-1

         !int_{tau}^{tau_L} Ib(tau')*E1(tau'-tau) dtau'
         diff = 1000._dp
         dtau = (tau_L - tau_cell(i))/10._dp
         int0 = 0._dp
         do while (diff.gt.tol)
!            if (dtau.le.small) exit
            taul = tau_cell(i)                                          !Integration starts from the optical thickness of the grid cell center
            IbE1_b(i) = 0._dp
            do while(taul.le.tau_L)
               Ib_tau = splint(tau_cell,Ib,Ib2array,taul,nx)
               expnum = max(taul - tau_cell(i),small)
               IbE1_b(i) = IbE1_b(i) + Ib_tau*expint(1,expnum)*dtau
               taul = taul + dtau
            enddo
            diff = abs((IbE1_b(i) - int0)/(int0 + small))*100._dp
            dtau = dtau/2._dp
            int0 = IbE1_b(i)
         enddo
      
         !int_{tau}^{tau_L} Ib(tau')*E2(tau'-tau) dtau'
         diff = 1000._dp
         dtau = (tau_L - tau_cell(i))/10._dp
         int0 = 0._dp
         do while (diff.gt.tol)
!            if (dtau.le.small) exit
            taul = tau_cell(i)
            IbE2_b(i) = 0._dp
            do while(taul.le.tau_L)
               Ib_tau = splint(tau_cell,Ib,Ib2array,taul,nx) 
               expnum = max(taul - tau_cell(i),small)
               IbE2_b(i) = IbE2_b(i) + Ib_tau*expint(2,expnum)*dtau
               taul = taul + dtau
            enddo
            diff = abs((IbE2_b(i) - int0)/(int0 + small))*100._dp
            dtau = dtau/2._dp
            int0 = IbE2_b(i)
         enddo
      enddo IbE_loop_2
      
      !-----------------------------------------------------------------
      !Compute radiative absorption, emission, source and flux
      !-----------------------------------------------------------------
      final_loop: do i=1,nx
         rademi(i) = fourpi*kappa(i)*Ib(i)                              !Emission
         radabs(i) = twopi*(Ib(1)*expint(2,tau_cell(i)) + &             !Absorption
                            Ib(nx)*expint(2,tau_L - tau_cell(i)) + &
                            IbE1_f(i) + IbE1_b(i) )*kappa(i)
                            !
                            
         source(i) = radabs(i) - rademi(i)                              !Source
         flux(i) = twopi*(Ib(1)*expint(3,tau_cell(i)) + &               !Flux
                            Ib(nx)*expint(3,tau_L - tau_cell(i)) + &
                            IbE2_f(i) - IbE2_b(i) )
      enddo final_loop
      if (present(absorption)) absorption = radabs
      if (present(emission)) emission = rademi
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      deallocate(k2array,Ib2array)
      deallocate(IbE1_b,IbE1_f,IbE2_b,IbE2_f)
      deallocate(radabs,rademi)
      deallocate(tau_cell)
      
   endsubroutine exact_1d_pp_old


   !====================================================================
   !Subroutine to compute the exact RTE solution for a cylindrical, 
   !non-scattering 1d medium (temporary code)   
   !====================================================================
   subroutine exact_1d_cyl(r,a,Sb,I0,nx,source)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: fourpi,pi,pio2
      use comp_functions, only: CheckMemAlloc
      use math_functions, only: expint,spline,splint
      use omp_lib
      implicit none
      integer,intent(in) :: nx
      integer :: i,ierr
      integer :: ipsi,itheta,npsi,ntau,ntheta
      real(dp),intent(in) :: I0,Sb(nx),a(nx),r(nx)
      real(dp),intent(out) :: source(nx)
      real(dp) :: dpsi,dr,dtheta,g1,g2,Gint,tau0
      real(dp),allocatable,dimension(:) :: tau,theta,psi,Rad_abs,Rad_emi
dpsi = I0 !delete
      !-----------------------------------------------------------------
      ! Discretization parameters (part 1)
      !-----------------------------------------------------------------
      ntheta = 41
      npsi = 41
      ntau = 41

      !-----------------------------------------------------------------
      ! Allocating arrays
      !-----------------------------------------------------------------
      allocate(tau(1:nx),stat=ierr)
      call CheckMemAlloc('tau',ierr)
      allocate(theta(1:ntheta),stat=ierr)
      call CheckMemAlloc('theta',ierr)
      allocate(psi(1:npsi),stat=ierr)
      call CheckMemAlloc('psi',ierr)
      allocate(Rad_abs(1:nx),stat=ierr)
      call CheckMemAlloc('Rad_abs',ierr)
      allocate(Rad_emi(1:nx),stat=ierr)
      call CheckMemAlloc('Rad_emi',ierr)

      !-----------------------------------------------------------------
      ! Discretization parameters (part 2)
      !-----------------------------------------------------------------
      dtheta = pi/real(ntheta-1,dp)
      theta(1) = 0._dp
      do i=2,ntheta
         theta(i) = theta(i-1) + dtheta
      enddo

      psi(1) = -0.5_dp*pi
      dpsi = 2*pi/real(npsi-1,dp)
      do i=2,npsi
         psi(i) = psi(i-1) + dpsi
      enddo
   
      !-----------------------------------------------------------------
      !Computing predefined optical thickness
      !-----------------------------------------------------------------
      tau = 0._dp
      do i=2,nx
         dr = r(i)-r(i-1)
         tau(i) = tau(i-1) + 0.5_dp*(a(i) + a(i-1))*dr
      enddo
      tau0 = tau(nx)

      !-----------------------------------------------------------------
      !Computing all integrals needed for G
      !-----------------------------------------------------------------
      !$OMP PARALLEL DO PRIVATE(g1,G2,Gint,itheta,ipsi)
      do i=1,nx
         if (a(i).lt.0.0000001_dp) cycle
         Gint = 0._dp
         do itheta=1,ntheta
            do ipsi=1,npsi
               if (psi(ipsi).le.(pi/2._dp)) then
                  g1 = g_integral(i,psi(ipsi),theta(itheta),&
                                  tau(i)*abs(sin(psi(ipsi))),tau0,&
                                  tau,Sb,nx,ntau,.true.)
                  g2 = g_integral(i,psi(ipsi),theta(itheta),&
                                  tau(i)*abs(sin(psi(ipsi))),tau(i),&
                                  tau,Sb,nx,ntau,.false.)
               else
                  g1 = g_integral(i,psi(ipsi),theta(itheta),&
                                  tau(i),tau0,tau,Sb,nx,ntau,.true.)
                  g2 = 0._dp
               endif
               Gint = Gint + (g1+g2)*sin(theta(itheta))*dtheta*dpsi
            enddo
         enddo
         Rad_abs(i) = a(i)*Gint
         Rad_emi(i) = 4._dp*pi*a(i)*Sb(i)
         source(i) = Rad_abs(i) - Rad_emi(i)
      enddo
      !$OMP END PARALLEL DO 
      
   endsubroutine exact_1d_cyl
   
   !====================================================================
   !Subroutine to compute the incident radiation integral for the exact
   !one-dimensional cylindrical solution
   !====================================================================
   real(dp) function g_integral(i_in,psi_in,theta_in,mintau,maxtau,&
                                tau,Sb,nx,ntau,plus)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use math_functions, only: linint
      use precision_parameters, only: big
      implicit none
      integer,intent(in) :: i_in,nx,ntau
      integer :: ii,nintpoints
      logical,intent(in) :: plus
      real(dp),intent(in) :: psi_in,theta_in,mintau,maxtau,tau(nx),Sb(nx)
      real(dp) :: gsum,dtau,taufrac,expterm,sqrtterm,taull,Sbll
      real(dp),allocatable,dimension(:) :: Sbl,taul
      logical :: isnumber,notinfty

      !-----------------------------------------------------------------
      !Define integration parameters and arrays
      !-----------------------------------------------------------------
      nintpoints = ntau
      dtau = (maxtau - mintau)/real(nintpoints-1,dp)
      allocate(Sbl(1:nintpoints))
      allocate(taul(1:nintpoints))

      !-----------------------------------------------------------------
      !Define interpolated arrays
      !-----------------------------------------------------------------
      taul(1) = mintau
      Sbl(1) = linint(tau,Sb,taul(1),nx)
      do ii=2,nintpoints
         taul(ii) = taul(ii-1) + dtau
         Sbl(ii) = linint(tau,Sb,taul(ii),nx)
      enddo
 
      !-----------------------------------------------------------------
      !Compute the integral    
      !-----------------------------------------------------------------
      gsum = 0._dp
      do ii=1,nintpoints-1
         taull = (taul(ii) + taul(ii+1))/2._dp
         Sbll = ( Sbl(ii) + Sbl(ii+1) )/2._dp
         sqrtterm = sqrt( taull**2._dp - (tau(i_in)*sin(psi_in))**2._dp )
         taufrac = taull / (sin(theta_in)* sqrtterm )

         if (plus) then
            expterm = - (tau(i_in)*cos(psi_in) + sqrtterm )/sin(theta_in)
         else
            expterm = - (tau(i_in)*cos(psi_in) - sqrtterm )/sin(theta_in)
         endif
       	
         if (isnan(Sbll*taufrac*exp(expterm)*dtau)) then
            isnumber = .false.
         else
            isnumber = .true.
         endif
      
         if (dabs(Sbll*taufrac*exp(expterm)*dtau).gt.big) then
            notinfty = .false.
         else
            notinfty = .true.
         endif
      
         if (isnumber.and.notinfty) gsum = gsum + Sbll*taufrac*exp(expterm)*dtau
            
      enddo
      g_integral = gsum

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      deallocate(taul)
      deallocate(Sbl)

   endfunction g_integral
      
!   !=============================================
!   !Subroutine to compute the exact RTE solution
!   !for a cylindrical, non-scattering 1d medium
!   !=============================================
!   subroutine exact_1d_cyl(r,kappa,Ib,I0,nx,source,flux)
   
!      !-------------------------
!      !Declaration of variables
!      !-------------------------
!      use constants, only: big,fourpi,pi,pio2
!      use comp_functions, only: CheckMemAlloc
!      use math_functions, only: expint,spline,splint
!      integer,intent(in) :: nx
!      integer :: i,ierr
!      real(dp),intent(in) :: I0,Ib(nx),kappa(nx),r(nx)
!      real(dp),intent(out) :: flux(nx),source(nx)
!      real(dp) :: diff_tau,diff_max,dpsi,ds,dtau,dtheta,taul,tol
!      real(dp) :: exp1,exp2,expw,frac,Iw,sqr,tau,tt0,ttl,s
!      real(dp) :: int_0,int_1,int_2,int_term
!      real(dp) :: cpsi,psi,theta,spsi,stheta
!      real(dp) :: G_0,G_diff,G_psi,G_theta,q_0,q_diff,q_psi,q_theta
!      real(dp) :: kappa_s,Ib_tau,tau_0
!      real(dp),allocatable,dimension(:) :: k2array,Ib2array,tau_cell
      
!      !---------------
!      !Set parameters
!      !---------------
!      tol = 1.e-1_dp                                                    !Tolerance for the integration
      
!      !----------------
!      !Allocate arrays
!      !----------------
!      allocate(k2array(nx),stat=ierr)
!      call CheckMemAlloc('k2array',ierr)
!      allocate(Ib2array(nx),stat=ierr)
!      call CheckMemAlloc('Ib2array',ierr)
!      allocate(tau_cell(nx),stat=ierr)
!      call CheckMemAlloc('tau_cell',ierr)
      
!      !----------------------------------------
!      !Compute cell-centered optical thickness
!      !----------------------------------------
!      call spline(r,kappa,nx,big,big,k2array)                           !Prepare spline for the subsequent interpolation
!      tau_cell(1) = 0._dp                                               !Optical tickness at the first point is zero
!      tau_loop: do i=2,nx         
!         diff_tau = 1000._dp                                            !Initial difference between two successive iterations
!         ds = (r(i) - r(i-1))/10._dp                                    !Initial finite increment for the integration
!         int_0 = tau_cell(i-1)                                          !Initial "old" integration value 
!         do while (diff_tau.gt.tol)
!            s = r(i-1) + ds/2._dp                                       !Initial optical path coordinate value (previous grid point position)
!            tau_cell(i) = tau_cell(i-1)                                 !kappa*ds integration is only performed between the previous and the
!                                                                        !  current grid point, so always add to the integration the optical
!                                                                        !  thickness computed up to the previous grid point
!            do while(s.lt.r(i))
!               kappa_s = splint(r,kappa,k2array,s,nx)                   !Local absorption coefficient (interpolated)
!               tau_cell(i) = tau_cell(i) + kappa_s*ds                   !Add to the kappa*ds integration
!               if (i.eq.102) write(*,*) kappa_s,ds,s,r(i-1),r(i)
!               s = s + ds                                               !Update the optical path coordinate value
!            enddo
!            diff_tau = abs((tau_cell(i) - int_0)/&
!                           (int_0 + small))*100._dp                     !Difference relative to the previous integration
!            ds = ds/2._dp                                               !Reduce integration increment for the next iteration step
!            int_0 = tau_cell(i)                                          !Update "old" integration value
!         enddo
!      enddo tau_loop
!      tau_0 = tau_cell(nx)                                              !Optical thickness at the right-hand boundary
      
!!      do i=1,nx
!!         write(*,*) r(i),tau_cell(i),splint(r,kappa,k2array,r(i),nx)
!!      enddo
!!      stop
!      !-----------------------
!      !Main calculation loops
!      !-----------------------
!      Iw = I0                                                           !Boundary condition
!      r_loop: do i=1,nx
         
!         diff_max = 1000._dp
!         q_0 = 0._dp
!         G_0 = 0._dp
!         tau = max(tau_cell(i),small)                                   !Local optical thickness
!         dtheta = pi/10._dp                                             !Integration increment in theta
!         dpsi = pio2/10._dp                                             !Integration increment in psi
!         do while (diff_max.gt.tol)
!            theta = dtheta/2._dp                                        !Initial theta value
!            q_theta = 0._dp; G_theta = 0._dp
!            theta_loop: do while (theta.le.pi)
!               stheta = max(dsin(theta),small) !Store sin(theta)
!               psi = dpsi/2._dp       !Initial psi value
!               q_psi = 0._dp; G_psi = 0._dp
!               psi_loop: do while (psi.le.pio2)
!                  spsi = max(dsin(psi),small)  !Store sin(psi)
!                  cpsi = max(dcos(psi),small)  !Store cos(psi)
!                  tt0 = tau_0/tau
!                  sqr = sqrt(tt0**2._dp - spsi**2._dp)
!                  expw = dexp(-tau*(cpsi + sqr)/stheta)
                  
!                  diff_tau = 1000._dp
!                  dtau = (tau - tau*abs(spsi))/10._dp
!                  int_0 = 0._dp
!                  do while (diff_tau.gt.tol)
!                     taul = tau*abs(spsi) + dtau/2._dp
!                     int_1 = 0._dp
!                     if (dabs(tau - tau*abs(spsi)).lt.small) exit
!                     tau_loop_1: do while (taul.le.tau)
!                        Ib_tau = splint(tau_cell,Ib,Ib2array,taul,nx)
!                        ttl = taul/tau
!                        sqr = sqrt(ttl**2._dp - spsi**2._dp)
!                        if (isnan(sqr)) cycle!sqr = small
!                        frac = ttl/(stheta*sqr + small)
!                        exp1 = dexp(-tau*(cpsi + sqr)/stheta)
!                        exp2 = dexp(-tau*(cpsi - sqr)/stheta)
!                        if (isnan(exp1).or.isnan(exp2)) cycle
!                        int_1 = int_1 + Ib_tau*frac*(exp1 + exp2)*dtau
!                        taul = taul + dtau
!                     enddo tau_loop_1
!                     diff_tau = dabs((int_1 - int_0)/(int_0 + small))*100._dp
!                     dtau = dtau/2._dp
!                     int_0 = int_1
!                  enddo
   
!                  diff_tau = 1000._dp
!                  dtau = (tau_0 - tau)/10._dp
!                  int_0 = 0._dp
!                  do while (diff_tau.gt.tol)
!                     taul = tau + dtau/2._dp
!                     int_2 = 0._dp
!                     if (dabs(tau_0 - tau).lt.small) exit
!                     tau_loop_2: do while (taul.le.tau_0)
!                        Ib_tau = splint(tau_cell,Ib,Ib2array,taul,nx)
!                        ttl = taul/tau 
!                        sqr = sqrt(ttl**2._dp - spsi**2._dp)
!                        if (isnan(sqr)) cycle!sqr = small
!                        frac = ttl/(stheta*sqr + small)
!                        exp1 = dexp(-tau*(cpsi + sqr)/stheta)
!                        if (isnan(exp1)) cycle
!                        int_2 = int_2 + Ib_tau*frac*exp1*dtau
!                        taul = taul + dtau
!                     enddo tau_loop_2
!                     diff_tau = dabs((int_2 - int_0)/(int_0 + small))*100._dp
!                     dtau = dtau/2._dp
!                     int_0 = int_2
!                  enddo
!                  int_term = Iw*expw + 4._dp*(int_1 + int_2)
!                  q_psi = q_psi + int_term*cpsi*dpsi
!                  G_psi = G_psi + int_term*dpsi
!                  psi = psi + dpsi
!               enddo psi_loop
!               q_theta = q_theta + q_psi*stheta*stheta*dtheta
!               G_theta = G_theta + G_psi*stheta*dtheta
!               theta = theta + dtheta
!            enddo theta_loop
            
!            q_diff = dabs((q_theta - q_0)/q_0)*100._dp
!            G_diff = dabs((G_theta - G_0)/q_0)*100._dp
!            diff_max = max(q_diff,G_diff)
!            dtheta = dtheta/2._dp
!            dpsi = dpsi/2._dp
!            q_0 = q_theta
!            G_0 = G_theta
!         enddo
!         flux(i) = q_theta
!         source(i) = kappa(i)*(G_theta - fourpi*Ib(i))
!write(*,*) i,source(i),flux(i)
!      enddo r_loop

!      !------------------
!      !Deallocate arrays
!      !------------------
!      deallocate(k2array)
!      deallocate(Ib2array)
!      deallocate(tau_cell)
   
!   endsubroutine exact_1d_cyl

   
   
   
   
!   subroutine exact_1d_cyl2(r,a,Sb,I0,nx,source)
   
!      use constants, only: fourpi,pi,pio2
!      use comp_functions, only: CheckMemAlloc
!      use math_functions, only: expint,spline,splint
!      implicit none
!      integer,intent(in) :: nx
!      integer :: i
!      integer :: ipsi,itheta,npsi,ntau,ntheta
!      real(dp),intent(in) :: I0,Sb(nx),a(nx),r(nx)
!      real(dp),intent(out) :: source(nx)
!      real(dp) :: dpsi,dr,dtheta,g1,g2,Gint,tau0
!      real(dp),allocatable,dimension(:) :: tau,theta,psi,Rad_abs,Rad_emi
      
      
!      !------------------------------------
!      ! Discretization parameters (part 1)
!      !------------------------------------
!      ntheta = 41
!      npsi = 41
!      ntau = 41

!      !-------------------
!      ! Allocating arrays
!      !-------------------
!      allocate(tau(1:nx))
!      allocate(theta(1:ntheta))
!      allocate(psi(1:npsi))
!      allocate(Rad_abs(1:nx))
!      allocate(Rad_emi(1:nx))

!      !------------------------------------
!      ! Discretization parameters (part 2)
!      !------------------------------------
!      dtheta = pi/real(ntheta-1,dp)
!      theta(1) = 0._dp
!      do i=2,ntheta
!         theta(i) = theta(i-1) + dtheta
!      enddo

!      psi(1) = -0.5_dp*pi
!      dpsi = 2*pi/real(npsi-1,dp)
!      do i=2,npsi
!         psi(i) = psi(i-1) + dpsi
!      enddo
   
!   	!---------------------------------------
!      !Computing predefined optical thickness
!      !---------------------------------------
!      tau = 0._dp
!      do i=2,nx
!         dr = r(i)-r(i-1)
!         tau(i) = tau(i-1) + 0.5_dp*(a(i) + a(i-1))*dr
!!write(*,*) i,tau(i)
!      enddo
!      tau0 = tau(nx)                        !Optical thickness tau0
   
!      !-------------------------------------
!      !Computing all integrals needed for G
!      !-------------------------------------
!      do i=1,nx
!         Gint = 0._dp
!         do itheta=1,ntheta
!            do ipsi=1,npsi
!write(*,*) i,itheta,ipsi
!               if (psi(ipsi).le.(pi/2._dp)) then
!                  g1 = g_integral(i,psi(ipsi),theta(itheta),&
!                        tau(i)*abs(sin(psi(ipsi))),tau0,tau,Sb,nx,ntau,.true.)
!                  g2 = g_integral(i,psi(ipsi),theta(itheta),&
!                        tau(i)*abs(sin(psi(ipsi))),tau(i),tau,Sb,nx,ntau,.false.)
!               else
!                  g1 = g_integral(i,psi(ipsi),theta(itheta),&
!                        tau(i),tau0,tau,Sb,nx,ntau,.true.)
!                  g2 = 0._dp
!               endif
!               Gint = Gint + (g1+g2)*sin(theta(itheta))*dtheta*dpsi
!            enddo
!         enddo
!         Rad_abs(i) = a(i)*Gint
!         Rad_emi(i) = 4._dp*pi*a(i)*Sb(i)
!         source(i) = Rad_abs(i) - Rad_emi(i)
!!write(*,*) i,a(i),Sb(i),Rad_abs(i),Rad_emi(i),source(i)
!      enddo
!   endsubroutine exact_1d_cyl2
   
   
!   real(dp) function g_integral(i_in,psi_in,theta_in,mintau,maxtau,&
!                                tau,Sb,nx,ntau,plus)

!      use precision_parameters, only: big
!      implicit none
!      integer,intent(in) :: i_in,nx,ntau
!      integer :: ii,nintpoints
!      logical,intent(in) :: plus
!      real(dp),intent(in) :: psi_in,theta_in,mintau,maxtau,tau(nx),Sb(nx)
!      real(dp) :: gsum,dtau,taufrac,expterm,sqrtterm,taull,Sbll
!      real(dp),allocatable,dimension(:) :: Sbl,taul
!      logical :: isnumber,notinfty

!      !Define integration parameters and arrays
!      nintpoints = ntau
!      allocate(Sbl(1:nintpoints))
!      allocate(taul(1:nintpoints))

!      !Define integration bounds
!      dtau = (maxtau - mintau)/real(nintpoints-1,dp)

!      !Define interpolated arrays
!      taul(1) = mintau
!      Sbl(1) = LinInterp(taul(1),tau,Sb)
!      do ii=2,nintpoints
!         taul(ii) = taul(ii-1) + dtau
!         Sbl(ii) = LinInterp(taul(ii),tau,Sb)
!      enddo
 
!      !Compute the integral    
!      gsum = 0._dp
!      do ii=1,nintpoints-1
!         taull = (taul(ii) + taul(ii+1))/2._dp
!         Sbll = ( Sbl(ii) + Sbl(ii+1) )/2._dp
      
!         sqrtterm = sqrt( taull**2._dp - (tau(i_in)*sin(psi_in))**2._dp )
!         taufrac = taull / (sin(theta_in)* sqrtterm )
            
!         if (plus) then
!            expterm = - (tau(i_in)*cos(psi_in) + sqrtterm )/sin(theta_in)
!         else
!            expterm = - (tau(i_in)*cos(psi_in) - sqrtterm )/sin(theta_in)
!         endif
       	
!         if (isnan(Sbll*taufrac*exp(expterm)*dtau)) then
!            isnumber = .false.
!         else
!            isnumber = .true.
!         endif
      
!         if (dabs(Sbll*taufrac*exp(expterm)*dtau).gt.big) then
!            notinfty = .false.
!         else
!            notinfty = .true.
!         endif
      
!         if (isnumber.and.notinfty) gsum = gsum + Sbll*taufrac*exp(expterm)*dtau
            
!      enddo
!      g_integral = gsum

!      !Deallocate arrays
!      deallocate(taul)
!      deallocate(Sbl)

!   endfunction g_integral    
   
   
   real(dp) function LinInterp(xtarget,xarray,yarray)

      integer :: ii,ntotal
      real(dp),intent(in) :: xtarget,xarray(:),yarray(:)
      real(dp) :: x1,y1,rxy
      
      ntotal = size(xarray)
   
      !Find bounds
      if (xtarget.le.xarray(1)) then
         x1 = 0._dp
         y1 = 0._dp
         rxy = yarray(1)/(xarray(1)+small)
      elseif(xtarget.ge.xarray(ntotal)) then
         x1 = xarray(ntotal)
         y1 = yarray(ntotal)
         rxy = (yarray(ntotal) - yarray(ntotal-1))/(xarray(ntotal) - xarray(ntotal-1) + small)
      else
         do ii=2,ntotal
            if ((xtarget.gt.xarray(ii-1)).and.(xtarget.le.xarray(ii))) then
               x1 = xarray(ii-1)
               y1 = yarray(ii-1)
               rxy = (yarray(ii) - yarray(ii-1)) / (xarray(ii) - xarray(ii-1) + small)
            endif
         enddo
      endif

      !Interpolate
      LinInterp = y1 + rxy*(xtarget - x1)
   
   endfunction LinInterp

endmodule exact_solution

