!#######################################################################
!Module containing all mesh-related 
!parameters, subroutines and functions
!#######################################################################
module mesh

   !====================================================================
   !Declaration of variables
   !====================================================================
   use precision_parameters, only: dp
   implicit none
   character(1) :: mesh_practice
   integer :: xcells,xfaces,xpoints
   integer :: ycells,yfaces,ypoints
   integer :: zcells,zfaces,zpoints
   real(dp) :: xmin,xmax,ymin,ymax,zmin,zmax
   real(dp),allocatable,dimension(:) :: x,y,z,xf,yf,zf
   real(dp),allocatable,dimension(:,:) :: surface_bands
   real(dp),allocatable,dimension(:,:,:) :: T,p
   real(dp),allocatable,dimension(:,:,:) :: xmin_emissivity,xmax_emissivity
   real(dp),allocatable,dimension(:,:,:) :: ymin_emissivity,ymax_emissivity
   real(dp),allocatable,dimension(:,:,:) :: zmin_emissivity,zmax_emissivity
   real(dp),allocatable,dimension(:,:,:,:) :: xs
   logical :: one_d,two_d,three_d

contains
   !====================================================================
   !Subroutine to build the spatial grid
   !====================================================================
   subroutine build_spatial_mesh
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,shutdown
      implicit none
      integer :: i,ierr
      integer :: nx,nxf,ny,nyf,nz,nzf
      real(dp) :: dx,dy,dz      
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      nx = xpoints; ny = ypoints; nz = zpoints
      nxf = xfaces; nyf = yfaces; nzf = zfaces
      
      !-----------------------------------------------------------------
      !Allocate arrays with grid point and face positions
      !-----------------------------------------------------------------
      if (allocated(x)) deallocate(x); if (allocated(xf)) deallocate(xf)
      if (allocated(y)) deallocate(y); if (allocated(yf)) deallocate(yf)
      if (allocated(z)) deallocate(z); if (allocated(zf)) deallocate(zf)
      
      allocate(x(1:nx),stat=ierr); call CheckMemAlloc('x',ierr)
      allocate(y(1:ny),stat=ierr); call CheckMemAlloc('y',ierr)
      allocate(z(1:nz),stat=ierr); call CheckMemAlloc('z',ierr)
      allocate(xf(1:nxf),stat=ierr); call CheckMemAlloc('xf',ierr)
      allocate(yf(1:nyf),stat=ierr); call CheckMemAlloc('yf',ierr)
      allocate(zf(1:nzf),stat=ierr); call CheckMemAlloc('zf',ierr)
      
      !-----------------------------------------------------------------
      !Cell sizes (uniform for now)
      !-----------------------------------------------------------------
      dx = (xmax - xmin)/real(xcells,dp)                                !x direction
      dy = (ymax - ymin)/real(ycells,dp)                                !y direction
      dz = (zmax - zmin)/real(zcells,dp)                                !z direction
      selectcase(mesh_practice)
      case('A')
         !--------------------------------------------------------------
         !Practice A: faces located midway between grid points
         !--------------------------------------------------------------
         !First define the position of the grid points
         do i=1,nx
            x(i) = xmin + real(i-1,dp)*dx                               !Grid point position, x direction
         enddo
         do i=1,ny
            y(i) = ymin + real(i-1,dp)*dy                               !Grid point position, y direction
         enddo
         do i=1,nz
            z(i) = zmin + real(i-1,dp)*dz                               !Grid point position, z direction
         enddo
         
         !Faces are located midway between grid points
         xf(1) = x(1); xf(nxf) = x(nx)                                  !First and last face, x direction
         yf(1) = y(1); yf(nyf) = y(ny)                                  !First and last face, y direction
         zf(1) = z(1); zf(nzf) = z(nz)                                  !First and last face, z direction
         do i=2,nxf-1
            xf(i) = (x(i) + x(i-1))/2._dp                               !Face position, x direction
         enddo
         do i=2,nyf-1
            yf(i) = (y(i) + y(i-1))/2._dp                               !Face position, y direction
         enddo
         do i=2,nzf-1
            zf(i) = (z(i) + z(i-1))/2._dp                               !Face position, z direction
         enddo
      
      case('B')
         !--------------------------------------------------------------
         !Practice B: grid points located midway between faces
         !--------------------------------------------------------------
         !First define the position of the cell faces
         do i=1,nxf
            xf(i) = xmin + real(i-1,dp)*dx                              !Grid point position, x direction
         enddo
         do i=1,nyf
            yf(i) = ymin + real(i-1,dp)*dy                              !Grid point position, y direction
         enddo
         do i=1,nzf
            zf(i) = zmin + real(i-1,dp)*dz                              !Grid point position, z direction
         enddo
         
         !Grid points are located midway between the faces
         x(1) = xf(1); x(nx) = xf(nxf)                                  !First and last point, x direction
         y(1) = yf(1); y(ny) = yf(nyf)                                  !First and last point, y direction
         z(1) = zf(1); z(nz) = zf(nzf)                                  !First and last point, z direction
         do i=2,nx-1
            x(i) = (xf(i) + xf(i-1))/2._dp                              !Grid point position, x direction
         enddo
         do i=2,ny-1
            y(i) = (yf(i) + yf(i-1))/2._dp                              !Grid point position, y direction
         enddo
         do i=2,nz-1
            z(i) = (zf(i) + zf(i-1))/2._dp                              !Grid point position, z direction
         enddo
      case default
         call shutdown('build_spatial_mesh: mesh_practice unspecified')
      endselect
      
   endsubroutine build_spatial_mesh
   
   !====================================================================
   !Subroutine to initialize all mesh variables
   !====================================================================
   subroutine initialize_mesh_variables
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,shutdown
      use global_parameters, only: number_of_species,&
         number_of_surface_bands
      implicit none
      integer :: imin,imax,jmin,jmax,kmin,kmax
      integer :: ierr,padd,nb,ns
      logical :: dim_error
      
      !-----------------------------------------------------------------
      !Preliminary cheks of dimensionality
      !-----------------------------------------------------------------
      !Check if no dimensionality is specified
      if ((.not.one_d).and.(.not.two_d).and.(.not.three_d)) &
         call shutdown('Dimensionality of the problem not defined')
      
      !Check if more than one dimensionality is specified   
      dim_error = .false.
      if (one_d.and.two_d)    dim_error = .true.
      if (one_d.and.three_d)  dim_error = .true.
      if (two_d.and.three_d)  dim_error = .true.
      if (dim_error) &
         call shutdown('Problem with dimensionality definition')
      
      !-----------------------------------------------------------------
      !For 1D and 2D solutions, the non-used dimensions 
      !only have one grid cell (with index 1)
      !-----------------------------------------------------------------
      if (one_d) then
         ycells = 1; zcells = 1
      elseif (two_d) then
         zcells = 1
      endif
      
      !-----------------------------------------------------------------
      !Define number of grid points
      !(equal to the number of cells for Practice A;
      !equal to 2 + number of cells for Practice B)
      !-----------------------------------------------------------------
      if (mesh_practice.eq.'A') padd = 0
      if (mesh_practice.eq.'B') padd = 2
      xpoints = xcells + padd
      ypoints = ycells + padd; if (one_d)          ypoints = ycells
      zpoints = zcells + padd; if (one_d.or.two_d) zpoints = zcells
      
      !-----------------------------------------------------------------
      !Define number of faces
      !-----------------------------------------------------------------
      xfaces = xcells + 1
      yfaces = ycells + 1; if (one_d) yfaces = ycells
      zfaces = zcells + 1; if (one_d.or.two_d) zfaces = zcells
      
      !-----------------------------------------------------------------
      !Set bounds for allocation
      !-----------------------------------------------------------------
      imin = 1; imax = xpoints
      jmin = 1; jmax = ypoints
      kmin = 1; kmax = zpoints
      
      !-----------------------------------------------------------------
      !Allocating
      !-----------------------------------------------------------------
      if (allocated(T))                deallocate(T)
      if (allocated(xs))               deallocate(xs)
      if (allocated(p))                deallocate(p)
      if (allocated(surface_bands))    deallocate(surface_bands)
      if (allocated(xmin_emissivity))  deallocate(xmin_emissivity)
      if (allocated(xmax_emissivity))  deallocate(xmax_emissivity)
      if (allocated(ymin_emissivity))  deallocate(ymin_emissivity)
      if (allocated(ymax_emissivity))  deallocate(ymax_emissivity)
      if (allocated(zmin_emissivity))  deallocate(zmin_emissivity)
      if (allocated(zmax_emissivity))  deallocate(zmax_emissivity)
      
      ns = number_of_species; nb = number_of_surface_bands
      allocate(T(imin:imax,jmin:jmax,kmin:kmax),stat=ierr)
      call CheckMemAlloc('T',ierr)
      allocate(xs(ns,imin:imax,jmin:jmax,kmin:kmax),stat=ierr)
      call CheckMemAlloc('xs',ierr)
      allocate(p(imin:imax,jmin:jmax,kmin:kmax),stat=ierr)
      call CheckMemAlloc('p',ierr)
      allocate(surface_bands(nb,2),stat=ierr)
      call CheckMemAlloc('xmin_emissivity',ierr)
      allocate(xmin_emissivity(jmin:jmax,kmin:kmax,nb),stat=ierr)
      call CheckMemAlloc('xmin_emissivity',ierr)
      allocate(xmax_emissivity(jmin:jmax,kmin:kmax,nb),stat=ierr)
      call CheckMemAlloc('xmax_emissivity',ierr)
      allocate(ymin_emissivity(imin:imax,kmin:kmax,nb),stat=ierr)
      call CheckMemAlloc('ymin_emissivity',ierr)
      allocate(ymax_emissivity(imin:imax,kmin:kmax,nb),stat=ierr)
      call CheckMemAlloc('ymax_emissivity',ierr)
      allocate(zmin_emissivity(imin:imax,jmin:jmax,nb),stat=ierr)
      call CheckMemAlloc('zmin_emissivity',ierr)
      allocate(zmax_emissivity(imin:imax,jmin:jmax,nb),stat=ierr)
      call CheckMemAlloc('zmax_emissivity',ierr)
   
      !-----------------------------------------------------------------
      !Giving an initial value to the arrays
      !-----------------------------------------------------------------
      T = 0._dp; xs = 0._dp; p = 1._dp
      xmin_emissivity = 1._dp; xmax_emissivity = 1._dp
      ymin_emissivity = 1._dp; ymax_emissivity = 1._dp
      zmin_emissivity = 1._dp; zmax_emissivity = 1._dp
      
   endsubroutine initialize_mesh_variables

   !====================================================================
   !Subroutine to compute the average composition
   !====================================================================
   subroutine get_average_composition(T_avg,xs_avg,p_avg,which_average)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc
      implicit none
      character,intent(in),optional :: which_average
      integer :: ierr
      integer :: i,j,k
      integer :: nspc,nx,ny,nz
      logical :: T4_avg
      real(dp),intent(out) :: p_avg,T_avg,xs_avg(:)
      real(dp) :: dV,dx,dy,dz,vol
      real(dp) :: pow_T,p_sum,T_sum
      real(dp),allocatable,dimension(:) :: xs_sum
      
      !-----------------------------------------------------------------
      !Set up optional parameters
      !-----------------------------------------------------------------
      T4_avg = .false.
      if (present(which_average)) then
         if (trim(which_average).eq.'T4') T4_avg = .true.
      endif
      
      !-----------------------------------------------------------------
      !Mesh parameters
      !-----------------------------------------------------------------
      !Array sizes
      nspc = size(xs_avg)                                               !Number of species
      nx = xcells
      ny = ycells
      nz = zcells
   
      !Grid cell size (uniform)
      dx = (xmax - xmin)/real(nx,dp)                                    !x direction
      dy = (ymax - ymin)/real(ny,dp); if (one_d) dy = 1._dp             !y direction
      dz = (zmax - zmin)/real(nz,dp); if (one_d.or.two_d) dz = 1._dp    !z direction
      dV = dx*dy*dz
      
      !Volume of the domain
      vol = (xmax - xmin)
      if (two_d.or.three_d) vol = vol*(ymax - ymin)
      if (three_d) vol = vol*(zmax - zmin)
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      allocate(xs_sum(1:nspc),stat=ierr)
      call CheckMemAlloc('xs_sum',ierr)
      
      !-----------------------------------------------------------------
      !Sum loop
      !-----------------------------------------------------------------
      T_sum = 0._dp; p_sum = 0._dp; xs_sum = 0._dp
      do i=1,nx
         do j=1,ny
            do k=1,nz
               pow_T = 1._dp; if (T4_avg) pow_T = 4._dp
               T_sum = T_sum + (T(i,j,k)**pow_T)*dV
               p_sum = p_sum + p(i,j,k)*dV
               xs_sum(:) = xs_sum(:) + xs(:,i,j,k)*dV
            enddo
         enddo
      enddo
      
      !Finish calculation of the averages
      T_avg = (T_sum/vol)**(1._dp/pow_T)
      p_avg = p_sum/vol
      xs_avg(:) = xs_sum(:)/vol   
   
   endsubroutine get_average_composition

endmodule mesh
