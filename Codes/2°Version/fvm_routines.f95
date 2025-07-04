module fvm_routines

   !====================================================================
   !Modules & Misc
   !====================================================================
   use mesh
   use fvm_parameters
   use precision_parameters, only: dp,small
   implicit none
   
contains
   
   !====================================================================
   !Main subroutine to solve the RTE with the finite volume method
   !====================================================================
   subroutine fvm_solution(k_abs,k_emi,k_sct,Ib,phi,x_eps,y_eps,z_eps,&
                           qrad,Srad,dont_iterate,total_loss)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: pi,fourpi
      use comp_functions, only: CheckMemAlloc,print_to_prompt,shutdown
      use global_parameters, only: debug_mode
      implicit none
      integer :: i,j,k,l,ll,n,nx,ny,nz,nl,nn
      integer :: icm,icp,iwm,iwp
      integer :: jcm,jcp,jwm,jwp
      integer :: kcm,kcp,kwm,kwp
      integer :: iend,istart,istep,iwall
      integer :: jend,jstart,jstep,jwall
      integer :: kend,kstart,kstep,kwall
      integer :: ierr
      real(dp),intent(in) :: k_abs(xpoints,ypoints,zpoints),&
                             k_emi(xpoints,ypoints,zpoints),&
                             k_sct(xpoints,ypoints,zpoints),&
                             Ib(xpoints,ypoints,zpoints),&
                             phi(fvm_angles,fvm_angles),&
                             x_eps(2,ypoints,zpoints),&
                             y_eps(2,xpoints,zpoints),&
                             z_eps(2,xpoints,ypoints)
      real(dp),intent(out) :: qrad(3,xpoints,ypoints,zpoints),&
                              Srad(xpoints,ypoints,zpoints)
      real(dp),intent(out),optional :: total_loss
      real(dp) :: Ax,Ay,Az,dx,dy,dz,V
      real(dp) :: A_ind,A_ijk,A_xu,A_yu,A_zu,AI_xu,AI_yu,AI_zu,B
      real(dp) :: gamma_x,gamma_y,gamma_z
      real(dp) :: eps_w,I_w,Ib_w,Iphi
      real(dp) :: Idiff,Idiff_max,iter_tol
      real(dp),allocatable,dimension(:,:,:) :: x_inc,y_inc,z_inc,&
                                               I_xu,I_yu,I_zu
      real(dp),allocatable,dimension(:,:,:,:) :: Irad,Irad0
      logical,intent(in),optional :: dont_iterate
      logical :: skip_iteration
      
      if (mesh_practice.eq.'A') &
         call shutdown('fvm_solution: Practice A not yet implemented')
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      if (debug_mode) &
         call print_to_prompt('fvm_solution: surrogate names')
      nx = xpoints; ny = ypoints; nz = zpoints                          !Number of grid points along directions x, y and z
      nl = fvm_angles                                                   !Number of angles in the FVM discretization
      nn = fvm_iterations                                               !Maximum number of iterations in the solution
      iter_tol = fvm_tolerance                                          !Tolerance for the iterative process
           
      !-----------------------------------------------------------------
      !Configure optional parameter
      !-----------------------------------------------------------------
      skip_iteration = .false.
      if (present(dont_iterate)) skip_iteration = dont_iterate
      
      !-----------------------------------------------------------------
      !Set up limiting grid indexes (in the current implementation, a
      !half-grid cell neighbors all domain borders)
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('fvm_solution: &
                                         &set up limiting grid indexes')
      iwm = 1                                                           !x-index of the first wall cell
      iwp = nx                                                          !x-index of the last wall cell
      icm = 2                                                           !x-index of the first grid cell in the medium
      icp = nx-1                                                        !x-index of the last grid cell in the medium

      jwm = 1                                                           !y-index of the first wall cell
      jwp = ny                                                          !y-index of the first wall cell
      jcm = 2                                                           !y-index of the first grid cell in the medium
      jcp = ny-1                                                        !y-index of the last grid cell in the medium
      if (one_d) then                                                   !Correct the cell y-indexes for the case of
         jcm = 1; jcp = 1                                               !  one-dimensional geometry
      endif

      kwm = 1                                                           !z-index of the first wall cell
      kwp = nz                                                          !z-index of the first wall cell
      kcm = 2                                                           !z-index of the first grid point in the medium
      kcp = nz-1                                                        !z-index of the last grid cell in the medium
      if (one_d.or.two_d) then                                          !Correct the cell z-indexes for the case of
         kcm = 1; kcp = 1                                               !  one- or two-dimensional geometries
      endif

      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('fvm_solution: &
                                           &allocate arrays')
      allocate(x_inc(2,jwm:jwp,kwm:kwp),stat=ierr)
      call CheckMemAlloc('x_inc',ierr)
      allocate(y_inc(2,iwm:iwp,kwm:kwp),stat=ierr)
      call CheckMemAlloc('y_inc',ierr)
      allocate(z_inc(2,iwm:iwp,jwm:jwp),stat=ierr)
      call CheckMemAlloc('z_inc',ierr)
      allocate(I_xu(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)
      call CheckMemAlloc('I_xu',ierr)
      allocate(I_yu(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)
      call CheckMemAlloc('I_yu',ierr)
      allocate(I_zu(iwm:iwp,jwm:jwp,kwm:kwp),stat=ierr)
      call CheckMemAlloc('I_zu',ierr)
      allocate(Irad(iwm:iwp,jwm:jwp,kwm:kwp,nl),stat=ierr)
      call CheckMemAlloc('Irad',ierr)
      allocate(Irad0(iwm:iwp,jwm:jwp,kwm:kwp,nl),stat=ierr)
      call CheckMemAlloc('Irad0',ierr)
      
      !-----------------------------------------------------------------
      !Determine the values of gamma
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('fvm_solution: &
                                           &define gamma')
      call define_gamma(gamma_x)
      call define_gamma(gamma_y)
      call define_gamma(gamma_z)
      
      !-----------------------------------------------------------------
      !Initialize geometry parameters
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('fvm_solution: &
                                       &initialize geometry parameters')
      dy = 1._dp; dz = 1._dp                                            !Cell sizes in the y and z directions
      Ay = 0._dp; Az = 0._dp                                            !Cell areas in the y and z directions
      AI_yu = 0._dp; AI_zu = 0._dp                                      !Contribution of the upstream intensities
                                                                        !  in the y and z directions to the local intensity
      
      !-----------------------------------------------------------------
      !Initialize parameters for the iterative procedure
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('fvm_solution: &
                    &initialize parameters for the iterative procedure')
      n = 0                                                             !Initialize iteration counter
      Irad = 0._dp                                                      !Initialize current radiative intensity
      Irad0 = 0._dp                                                     !Initialize previous radiative intensity
      
      iter_loop: do
         !--------------------------------------------------------------
         !Check if maximum number of iterations has been reached
         !--------------------------------------------------------------
         n = n + 1                                                      !Update iteration counter
         if (n.gt.nn) &
            call shutdown('fvm_solution: &
                          &maximum number of FVM iteractions reached')
         
         !--------------------------------------------------------------
         !Pre-compute the incident intensity at the boundaries 
         !--------------------------------------------------------------
         x_inc = 0._dp; y_inc = 0._dp; z_inc = 0._dp
         bc_loop: do l=1,nl
            !x faces
            do j=jcm,jcp
               do k=kcm,kcp
                  if (dlx(l).lt.0) then                                 !Incoming at the left wall
                     iwall = 1; I_w = Irad0(iwm,j,k,l)
                  else                                                  !Incoming at the right wall
                     iwall = 2; I_w = Irad0(iwp,j,k,l)
                  endif 
               
                  !Integrate the incident intensity
                  x_inc(iwall,j,k) = x_inc(iwall,j,k) +&
                     I_w*dabs(dlx(l))*domega(l)
               enddo
            enddo
            if (one_d) cycle bc_loop                                    !Skip BC for the cells in the y and z directions
                                                                        !  for one-dimensional geometries
            !y faces
            do i=icm,icp
               do k=kcm,kcp
                  if (dly(l).lt.0) then                                 !Incoming at the bottom wall
                     jwall = 1; I_w = Irad0(i,jwm,k,l)
                  else                                                  !Incoming at the top wall
                     jwall = 2; I_w = Irad0(i,jwp,k,l)
                  endif 
            
                  !Integrate the incident intensity
                  y_inc(jwall,i,k) = y_inc(jwall,i,k) +&
                     I_w*dabs(dly(l))*domega(l)               
               enddo
            enddo
            if (two_d) cycle bc_loop                                    !Skip BC for the cells in the z direction
                                                                        !  for two-dimensional geometries

            !z faces
            do i=icm,icp
               do j=jcm,jcp
                  if (dlz(l).lt.0) then                                 !Incoming at the back wall
                     kwall = 1; I_w = Irad0(i,j,kwm,l)
                  else                                                  !Incoming at the front wall
                     kwall = 2; I_w = Irad0(i,j,kwp,l)
                  endif 
            
                  !Integrate the incident intensity
                  z_inc(kwall,i,j) = z_inc(kwall,i,j) +&
                     I_w*dabs(dlz(l))*domega(l)               
               enddo
            enddo
            
         enddo bc_loop
         
         solution_loop: do l=1,nl
            !-----------------------------------------------------------
            !Define parameters of the line-of-sight loop
            !-----------------------------------------------------------
            !Direction x
            istart = icm; iend = icp; istep = 1
            iwall    = 1                                                !1: left wall; 2: right wall
            if (dlx(l).lt.0) then
               istart = icp; iend = icm; istep = -1
               iwall = 2
            endif
               
            !Direction y
            jstart = jcm; jend = jcp; jstep = 1
            jwall = 1                                                   !1: bottom wall; 2: top wall
            if (dly(l).lt.0) then
               jstart = jcp; jend = jcm; jstep = -1
               jwall = 2
            endif
               
            !Direction z
            kstart = kcm; kend = kcp; kstep = 1
            kwall = 1                                                   !1: back wall; 2: front wall
            if (dlz(l).lt.0) then
               kstart = kcp; kend = kcm; kstep = -1
               kwall = 2
            endif

            los_xloop: do i=istart,iend,istep
               los_yloop: do j=jstart,jend,jstep
                  los_zloop: do k=kstart,kend,kstep
                     !--------------------------------------------------
                     !Determine the intensity at the first cell
                     !--------------------------------------------------
                     !x direction
                     if (i.eq.istart) then
                        eps_w = x_eps(iwall,j,k)                        !Wall emissivity
                        Ib_w = Ib(istart-istep,j,k)                     !Blackbody intensity at the wall
                        I_w = x_inc(iwall,j,k)                          !Incident intensity at the wall
                        I_xu(istart-istep,j,k) = eps_w*Ib_w + &         !Outgoing intensity at the wall
                           ((1._dp-eps_w)/pi)*I_w
                        Irad(istart-istep,j,k,l) = &                    !Intensity at the first grid point
                           I_xu(istart-istep,j,k)                       !  equals the one at the wall
                     endif

                     !y direction (only for 2D and 3D geometries)
                     if ((.not.one_d).and.(j.eq.jstart)) then
                        eps_w = y_eps(jwall,i,k)                        !Wall emissivity
                        Ib_w = Ib(i,jstart-jstep,k)                     !Blackbody intensity at the wall
                        I_w = y_inc(jwall,i,k)                          !Incident intensity at the wall
                        I_yu(i,jstart-jstep,k) = eps_w*Ib_w + &         !Outgoing intensity at the wall
                           ((1._dp-eps_w)/pi)*I_w
                        Irad(i,jstart-jstep,k,l) = &                    !Intensity at the first grid point
                           I_yu(i,jstart-jstep,k)                       !  equals the one at the wall
                     endif 

                     !z direction (only for 3D geometries)
                     if ((three_d).and.(k.eq.kstart)) then
                        eps_w = z_eps(kwall,i,j)                        !Wall emissivity
                        Ib_w = Ib(i,j,kstart-kstep)                     !Blackbody intensity at the wall
                        I_w = z_inc(kwall,i,j)                          !Incident intensity at the wall
                        I_zu(i,j,kstart-kstep) = eps_w*Ib_w + &         !Outgoing intensity at the wall
                           ((1._dp-eps_w)/pi)*I_w
                        Irad(i,j,kstart-kstep,l) = &                    !Intensity at the first grid point
                           I_zu(i,j,kstart-kstep)                       !  equals the one at the wall
                     endif                      

                     !Volume size
                     dx = xf(i) - xf(i-1)                               !Cell size: x direction
                     if (.not.one_d) dy = yf(j) - yf(j-1)               !Cell size: y direction
                     if (three_d)    dz = zf(k) - zf(k-1)               !Cell size: z direction
                     V = dx*dy*dz

                     !Area of the cell faces
                     Ax = dy*dz                                         !Face normal to the x direction
                     if (.not.one_d) Ay = dx*dz                         !Face normal to the x direction
                     if (three_d)    Az = dx*dy                         !Face normal to the x direction

                     !Compute summation of the incoming intensity times  
                     !scattering phase function over all incoming angles
                     Iphi = 0._dp
                     do ll = 1,nl
                        if (ll.eq.l) cycle                              !Skip the current angle
                        Iphi = Iphi + Irad0(i,j,k,ll)*phi(l,ll)
                     enddo

                     !Parameters of the discretized equation
                     A_xu  = Ax*dabs(dlx(l))/(gamma_x + small)
                     A_yu  = Ay*dabs(dly(l))/(gamma_y + small)          !A_yu = 0 for 1D geometries
                     A_zu  = Az*dabs(dlz(l))/(gamma_z + small)          !A_zu = 0 for 1D and 2D geometries
                     A_ind = (k_abs(i,j,k) + k_sct(i,j,k) - &
                              k_sct(i,j,k)*phi(l,l)/fourpi)*V*domega(l)
                     A_ijk = A_xu + A_yu + A_zu + A_ind
                     B     = k_emi(i,j,k)*Ib(i,j,k)*V*domega(l) + &
                              k_sct(i,j,k)*Iphi*V*domega(l)/fourpi

                     !Solving the grid cell intensity
                     AI_xu = A_xu*I_xu(i-istep,j,k)
                     if (.not.one_d) AI_yu = A_yu*I_yu(i,j-jstep,k)
                     if (three_d)    AI_zu = A_zu*I_zu(i,j,k-kstep)
                     Irad(i,j,k,l) = (AI_xu + AI_yu + AI_zu + B)/&
                                     (A_ijk + small)

                     !Solving the intensity at the next face
                     I_xu(i,j,k) = (Irad(i,j,k,l) - (1._dp-gamma_x)*&
                                    I_xu(i-istep,j,k))/gamma_x
                     if (.not.one_d) I_yu(i,j,k) = &
                        (Irad(i,j,k,l) - (1._dp-gamma_y)*&
                           I_yu(i,j-jstep,k))/gamma_y
                     if (three_d) I_zu(i,j,k) = &
                        (Irad(i,j,k,l) - (1._dp-gamma_z)*&
                           I_zu(i,j,k-kstep))/gamma_z

                     !Finding the intensity at the last volume 
                     !(coincident with the last face)
                     if (i.eq.iend) &                                   !x direction
                        Irad(iend+istep,j,k,l) = I_xu(i,j,k)
                     if ((.not.one_d).and.(j.eq.jend)) &                !y direction (only for 2D and 3D geometries)
                        Irad(i,jend+jstep,k,l) = I_yu(i,j,k)
                     if ((three_d).and.(k.eq.kend)) &                   !z direction (only for 3D geometries)
                        Irad(i,j,kend+kstep,l) = I_zu(i,j,k)

                  enddo los_zloop
               enddo los_yloop
            enddo los_xloop

         enddo solution_loop
         
         !--------------------------------------------------------------
         !Check convergence
         !--------------------------------------------------------------
         Idiff_max = 0._dp
         do i=iwm,iwp
            do j=jwm,jwp
               do k=kwm,kwp
                  do l=1,nl
                     !Compute relative difference 
                     !to previous iteration
                     Idiff = (Irad(i,j,k,l) - Irad0(i,j,k,l))/&
                             (Irad(i,j,k,l) + small)
                  
                     !Store the maximum difference
                     if (abs(Idiff).gt.Idiff_max) &
                        Idiff_max = abs(Idiff)

                     !Update radiative intensity value
                     Irad0(i,j,k,l) = Irad(i,j,k,l)
                  enddo
               enddo
            enddo
         enddo
            
         !If the convergence criterion has been reached, 
         !escape the loop
         if (Idiff_max.lt.iter_tol) exit iter_loop
         
         !--------------------------------------------------------------
         !Skip the iteration loop if requested
         !--------------------------------------------------------------
         if (skip_iteration) exit iter_loop
         
      enddo iter_loop

      !-----------------------------------------------------------------
      !Compute the radiative heat flux and radiative heat source
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('fvm_solution: &
                                &compute the radiative flux and source')
      qrad = 0._dp; Srad = 0._dp
      do i=iwm,iwp
         do j=jwm,jwp
            do k=kwm,kwp
               do l=1,nl
                  !Heat fluxes
                  qrad(1,i,j,k) = qrad(1,i,j,k) + Irad(i,j,k,l)*dlx(l)  !x direction
                  if (.not.one_d) qrad(2,i,j,k) = qrad(2,i,j,k) + &     !y direction (only for 2D and 3D geometries)
                                                Irad(i,j,k,l)*dly(l)
                  if (three_d) qrad(3,i,j,k) = qrad(3,i,j,k) + &        !z direction (only for 3D geometries)
                                             Irad(i,j,k,l)*dlz(l)
            
                  !Heat source
                  Srad(i,j,k) = Srad(i,j,k) + &
                     (k_abs(i,j,k)*Irad(i,j,k,l) - &
                        k_emi(i,j,k)*Ib(i,j,k))*domega(l)
               enddo
            enddo
         enddo
      enddo

      !-----------------------------------------------------------------
      !Compute the domain-integrated net loss, if requested
      !-----------------------------------------------------------------
      if (present(total_loss)) then
         if (debug_mode) call print_to_prompt('fvm_solution: &
                                     &compute the total radiative loss')
         total_loss = 0._dp
         do i=icm,icp
            do j=jcm,jcp
               do k=kcm,kcp
                  dx = xf(i) - xf(i-1)                                  !Cell size: x direction
                  if (.not.one_d) dy = yf(j) - yf(j-1)                  !Cell size: y direction
                  if (three_d)    dz = zf(k) - zf(k-1)                  !Cell size: z direction
                  V = dx*dy*dz                                          !Cell volume
                  total_loss = total_loss - Srad(i,j,k)*V
               enddo
            enddo
         enddo
      endif

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      if (debug_mode) call print_to_prompt('fvm_solution: &
                                           &deallocate arrays')
      deallocate(x_inc,y_inc,z_inc)
      deallocate(I_xu,I_yu,I_zu)
      deallocate(Irad,Irad0)
      
   endsubroutine fvm_solution
   
endmodule fvm_routines
