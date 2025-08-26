module nbck_functions

   use nbck_parameters

contains

   !====================================================================
   !Subroutine to read the NBCK data from external files
   !====================================================================
   subroutine load_nbck_data
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,CheckFileExists,CheckMemAlloc,&
                                dprint,get_file_unit
      use global_parameters, only: bin_ext,dat_ext,spec_name
      implicit none
      character(10) :: str1,str2,str3
      character(200) :: file_name,info_file
      integer :: ierr,info_unit
      integer :: ikg,inb,isp,itg,ixs
      integer :: nkg,nnb,nsp,ntg,nxs
      integer,allocatable,dimension(:) :: file_unit,nnb_spc
      real(dp),allocatable,dimension(:,:) :: lnb_spc,unb_spc
      
      !-----------------------------------------------------------------
      !If the premixed NBCK has already been loaded, skip everything
      !-----------------------------------------------------------------
      if (nbck_loaded) return
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('load_nbck_data: Allocate arrays')
      nsp = nbck_data_nsp
      allocate(file_unit(nsp),stat=ierr)
      call CheckMemAlloc('file_unit',ierr)
      allocate(nnb_spc(nsp),stat=ierr)
      call CheckMemAlloc('nnb_spc',ierr)
      
      !-----------------------------------------------------------------
      !Initialize variables
      !-----------------------------------------------------------------
      call dprint('load_nbck_data: Initialize variables')
      nbck_data_nkg = 0
      nbck_data_ntg = 0
      nbck_data_nxs = 0
      
      !-----------------------------------------------------------------
      !First reading loop: array sizes
      !-----------------------------------------------------------------
      call dprint('load_nbck_data: First reading loop')
      reading_loop_1: do isp=1,nsp
         file_name = trim(nbck_prefix)//trim(nbck_data_file(isp))       !File name
         call CheckFileExists(file_name)                                !Check if file exists
         file_unit(isp) = get_file_unit()                               !Assign unit
         open(file=file_name,unit=file_unit(isp),action='read',&        !Open unit
              form='unformatted')
         read(file_unit(isp)) nbck_data_nkg(isp),nnb_spc(isp),&         !Read array sizes
                         nbck_data_ntg(isp),nbck_data_nxs(isp)
         if (isp.gt.1) call assert(nnb_spc(isp).eq.nnb_spc(isp-1))      !Check if the number of narrow bands is 
                                                                        !  the same across all data files
      enddo reading_loop_1
      number_nbck_bands = nnb_spc(1)

      !-----------------------------------------------------------------
      !Allocate NBCK data file arrays
      !-----------------------------------------------------------------
      call dprint('load_nbck_data: Allocate NBCK data file arrays')
      nkg = maxval(nbck_data_nkg)
      nnb = number_nbck_bands
      ntg = maxval(nbck_data_ntg)
      nxs = maxval(nbck_data_nxs)
      if (allocated(nbck_data_tg)) deallocate(nbck_data_tg)
      allocate(nbck_data_tg(ntg,nsp),stat=ierr)
      call CheckMemAlloc('nbck_data_tg',ierr)
      if (allocated(nbck_data_xs)) deallocate(nbck_data_xs)
      allocate(nbck_data_xs(nxs,nsp),stat=ierr)
      call CheckMemAlloc('nbck_data_xs',ierr)
      if (allocated(nbck_g)) deallocate(nbck_g)
      allocate(nbck_g(nkg,nsp),stat=ierr)
      call CheckMemAlloc('nbck_g',ierr)
      if (allocated(nbck_k)) deallocate(nbck_k)
      allocate(nbck_k(nkg,ntg,nxs,nnb,nsp),stat=ierr)
      call CheckMemAlloc('nbck_k',ierr)
      if (allocated(nbck_lbound)) deallocate(nbck_lbound)
      allocate(nbck_lbound(nnb),stat=ierr)
      call CheckMemAlloc('nbck_lbound',ierr)
      if (allocated(nbck_ubound)) deallocate(nbck_ubound)
      allocate(nbck_ubound(nnb),stat=ierr)
      call CheckMemAlloc('nbck_ubound',ierr)
      
      !-----------------------------------------------------------------
      !Second reading loop: array values
      !-----------------------------------------------------------------
      call dprint('load_nbck_data: Second reading loop')
      
      !Allocate auxiliary arrays
      allocate(lnb_spc(nnb,nsp),stat=ierr)
      call CheckMemAlloc('lnb_spc',ierr)
      allocate(unb_spc(nnb,nsp),stat=ierr)
      call CheckMemAlloc('unb_spc',ierr)
      
      !Read data
      reading_loop_2: do isp=1,nsp
         read(file_unit(isp)) lnb_spc(:,isp),unb_spc(:,isp)
         read(file_unit(isp)) nbck_data_tg(:,isp)
         read(file_unit(isp)) nbck_data_xs(:,isp)
         read(file_unit(isp)) nbck_g(:,isp)
         read(file_unit(isp)) nbck_k(:,:,:,:,isp)
         if (isp.gt.1) then
            do inb=1,nnb
               call assert(dabs(lnb_spc(inb,isp)-lnb_spc(inb,isp-1))&
                           .lt.1e-8_dp)
               call assert(dabs(unb_spc(inb,isp)-unb_spc(inb,isp-1))&
                           .lt.1e-8_dp)
            enddo
         endif
         
         close(file_unit(isp))
      enddo reading_loop_2
      nbck_lbound = lnb_spc(:,1)
      nbck_ubound = unb_spc(:,1)

      !-----------------------------------------------------------------
      !Deallocate auxiliary arrays
      !-----------------------------------------------------------------
      call dprint('load_nbck_data: Deallocate auxiliary arrays')
      deallocate(file_unit,nnb_spc)
      deallocate(lnb_spc,unb_spc)
      
      !-----------------------------------------------------------------
      !Dump info data
      !-----------------------------------------------------------------
      call dprint('load_nbck_data: Dump info data')
      if (nbck_dump_info) then
         !Open unit
         info_unit = get_file_unit()
         info_file = trim(nbck_prefix)//trim(nbck_info_file)
         open(file=info_file,unit=info_unit,&
           form='formatted',action='write')

         !Dump information
         write(str1,'(i4)') nnb
         write(info_unit,'(a,a)') &
            'Number of narrow bands: ',trim(adjustl(str1))
         write(info_unit,'(a3,2x,a10,2x,a10)') &
            'n','lbound [1/cm]','ubound [1/cm]'
         do inb=1,nnb
            write(str1,'(i3)') inb
            write(str2,'(f10.4)') nbck_lbound(inb)/100._dp
            write(str3,'(f10.4)') nbck_ubound(inb)/100._dp
            write(info_unit,'(a3,2x,a10,2x,a10)') trim(adjustl(str1)),&
               trim(adjustl(str2)),trim(adjustl(str3))
         enddo
         write(info_unit,*)
         
         write(str1,'(i3)') nsp
         write(info_unit,'(a,a)') &
            'Number of species: ',trim(adjustl(str1))       
         write(info_unit,*)
         
         do isp=1,nsp
            write(info_unit,'(a,a)') &
               'Species: ',trim(spec_name(isp))
            write(info_unit,'(6x,a3,2x,a)') 'i','Mole fractions'
            do ixs=1,nbck_data_nxs(isp)
               write(str1,'(i3)') ixs
               write(str2,'(f10.4)') nbck_data_xs(ixs,isp)
               write(info_unit,'(6x,a3,2x,a10)') &
                  trim(adjustl(str1)),trim(adjustl(str2))
            enddo
            write(info_unit,*)
            write(info_unit,'(6x,a3,2x,a)') 'i','Temperatures [K]'
            do itg=1,nbck_data_ntg(isp)
               write(str1,'(i3)') itg
               write(str2,'(f10.4)') nbck_data_tg(itg,isp)
               write(info_unit,'(6x,a3,2x,a10)') &
                  trim(adjustl(str1)),trim(adjustl(str2))
            enddo
            write(info_unit,*)
            write(info_unit,'(6x,a3,2x,a)') 'i','g'
            do ikg=1,nbck_data_nkg(isp)
               write(str1,'(i3)') ikg
               write(str2,'(f10.4)') nbck_g(ikg,isp)
               write(info_unit,'(6x,a3,2x,a10)') &
                  trim(adjustl(str1)),trim(adjustl(str2))
            enddo
            write(info_unit,*)
         enddo
      endif
      close(info_unit)

      !-----------------------------------------------------------------
      !Set flag to indicate that the NBCK has been loaded
      !-----------------------------------------------------------------
      nbck_loaded = .true.

   endsubroutine load_nbck_data

   !====================================================================
   !Subroutine to read the premixed NBCK data from external files
   !====================================================================
   subroutine load_premix_nbck_data
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,CheckMemAlloc,&
                                dprint,get_file_unit
      use global_parameters, only: bin_ext,dat_ext
      implicit none
      character(200) :: file_name
      integer :: ierr,file_unit
      integer :: nfv,nkg,nnb,ntg,nxc,nxw
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('load_premix_nbck_data: Allocate arrays')
      
      !-----------------------------------------------------------------
      !Initialize variables
      !-----------------------------------------------------------------
      call dprint('load_premix_nbck_data: Initialize variables')
      nbck_mix_nkg = 0
      nbck_mix_ntg = 0
      nbck_mix_nxc = 0
      nbck_mix_nxw = 0
      nbck_mix_nfv = 0
      
      !-----------------------------------------------------------------
      !First reading loop: array sizes
      !-----------------------------------------------------------------
      call dprint('load_premix_nbck_data: First reading loop')
      file_name = trim(nbck_prefix)//trim(nbck_premix_file)             !File name
      call CheckFileExists(file_name)                                   !Check if file exists
      file_unit = get_file_unit()                                       !Assign unit
      open(file=file_name,unit=file_unit,action='read',&                !Open unit
            form='unformatted')
      read(file_unit) nbck_mix_nkg,number_nbck_bands,nbck_mix_ntg,&     !Read array sizes
                      nbck_mix_nxc,nbck_mix_nxw,nbck_mix_nfv

      !-----------------------------------------------------------------
      !Allocate NBCK data file arrays
      !-----------------------------------------------------------------
      call dprint('load_premix_nbck_data: &
                  &Allocate NBCK data file arrays')
      nkg = nbck_mix_nkg
      nnb = number_nbck_bands
      ntg = nbck_mix_ntg
      nxc = nbck_mix_nxc
      nxw = nbck_mix_nxw
      nfv = nbck_mix_nfv
      if (allocated(nbck_mix_tg)) deallocate(nbck_mix_tg)
      allocate(nbck_mix_tg(ntg),stat=ierr)
      call CheckMemAlloc('nbck_mix_tg',ierr)
      if (allocated(nbck_mix_xc)) deallocate(nbck_mix_xc)
      allocate(nbck_mix_xc(nxc),stat=ierr)
      call CheckMemAlloc('nbck_mix_xc',ierr)
      if (allocated(nbck_mix_xw)) deallocate(nbck_mix_xw)
      allocate(nbck_mix_xw(nxw),stat=ierr)
      call CheckMemAlloc('nbck_mix_xw',ierr)
      if (allocated(nbck_mix_fv)) deallocate(nbck_mix_fv)
      allocate(nbck_mix_fv(nfv),stat=ierr)
      call CheckMemAlloc('nbck_mix_fv',ierr)
      if (allocated(nbck_mix_g)) deallocate(nbck_mix_g)
      allocate(nbck_mix_g(nkg),stat=ierr)
      call CheckMemAlloc('nbck_mix_g',ierr)
      if (allocated(nbck_mix_k)) deallocate(nbck_mix_k)
      allocate(nbck_mix_k(nkg,ntg,nxw,nxc,nfv,nnb),stat=ierr)
      call CheckMemAlloc('nbck_mix_k',ierr)
      if (allocated(nbck_lbound)) deallocate(nbck_lbound)
      allocate(nbck_lbound(nnb),stat=ierr)
      call CheckMemAlloc('nbck_lbound',ierr)
      if (allocated(nbck_ubound)) deallocate(nbck_ubound)
      allocate(nbck_ubound(nnb),stat=ierr)
      call CheckMemAlloc('nbck_ubound',ierr)
      
      !-----------------------------------------------------------------
      !Second reading loop: array values
      !-----------------------------------------------------------------
      call dprint('load_premix_nbck_data: Second reading loop')
      
      !Read data
      read(file_unit) nbck_lbound,nbck_ubound
      read(file_unit) nbck_mix_tg
      read(file_unit) nbck_mix_xc
      read(file_unit) nbck_mix_xw
      read(file_unit) nbck_mix_fv
      read(file_unit) nbck_mix_g
      read(file_unit) nbck_mix_k
      close(file_unit)

   endsubroutine load_premix_nbck_data

   !====================================================================
   !Subroutine to initialize NBCK parameters
   !====================================================================
   subroutine prepare_premixed_nbck(Tg_min,Tg_max,xs_min,xs_max)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc
      use global_parameters, only: id_h2o,id_co2,id_soot
      use math_functions, only: generate_distribution
      implicit none
      integer :: ierr
      integer :: nfv,nkg,nnb,ntg,nxc,nxw
      real(dp),intent(in) :: Tg_min,Tg_max,xs_min(:),xs_max(:)
      real(dp),allocatable,dimension(:) :: dummy_arr

      !-----------------------------------------------------------------
      !If the premixed NBCK has already been computed, skip everything
      !-----------------------------------------------------------------
!     if (nbck_is_premixed) return  
!Commented because I don't know if it is a good idea to have this flag

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      nnb = number_nbck_bands
      nkg = nbck_mix_nkg
      ntg = nbck_mix_ntg
      nxc = nbck_mix_nxc
      nxw = nbck_mix_nxw
      nfv = nbck_mix_nxw

      !-----------------------------------------------------------------
      !Allocate arrays needed if the NBCK is to be premixed
      !-----------------------------------------------------------------
      if (allocated(nbck_mix_g)) deallocate(nbck_mix_g)
      allocate(nbck_mix_g(nkg),stat=ierr)
      call CheckMemAlloc('nbck_mix_g',ierr)
      if (allocated(nbck_mix_tg)) deallocate(nbck_mix_tg)
      allocate(nbck_mix_tg(ntg),stat=ierr)
      call CheckMemAlloc('nbck_mix_tg',ierr)
      if (allocated(nbck_mix_xc)) deallocate(nbck_mix_xc)
      allocate(nbck_mix_xc(nxc),stat=ierr)
      call CheckMemAlloc('nbck_mix_xc',ierr)
      if (allocated(nbck_mix_xw)) deallocate(nbck_mix_xw)
      allocate(nbck_mix_xw(nxw),stat=ierr)
      call CheckMemAlloc('nbck_mix_xw',ierr)
      if (allocated(nbck_mix_fv)) deallocate(nbck_mix_fv)
      allocate(nbck_mix_xw(nfv),stat=ierr)
      call CheckMemAlloc('nbck_mix_fv',ierr)
      if (allocated(nbck_mix_k)) deallocate(nbck_mix_k)
      allocate(nbck_mix_k(nkg,ntg,nxw,nxc,nfv,nnb),stat=ierr)
      call CheckMemAlloc('nbck_mix_k',ierr)
    
      !-----------------------------------------------------------------
      !Fill in the arrays
      !-----------------------------------------------------------------
      !Temperature
      call generate_distribution('lin',Tg_min,Tg_max,ntg,&
                                 nbck_mix_tg)
      
      !H2O mole fraction
      call generate_distribution('lin',xs_min(id_h2o),xs_max(id_h2o),&
                                 nxw,nbck_mix_xw)
      
      !CO2 mole fraction
      call generate_distribution('lin',xs_min(id_co2),xs_max(id_co2),&
                                 nxc,nbck_mix_xc)
      
      !Soot volume fraction
      call generate_distribution('lin',xs_min(id_soot),xs_max(id_soot),&
                                 nxc,nbck_mix_xc)
      
      !g points
      allocate(dummy_arr(nkg))
      call get_ck_quadrature(nbck_quadrature,nbck_mix_g,dummy_arr) 
      deallocate(dummy_arr)

   endsubroutine prepare_premixed_nbck

   !====================================================================
   !Subroutine to carry out a pre-mixing of the NBCK data from the
   !single species data that is already loaded (for now, only for a
   !CO2-H2O mixture)
   !====================================================================
   subroutine premix_nbck_data(band_index,read_from_file,write_to_file)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,dprint,get_file_unit
      use global_parameters, only: id_h2o,id_co2,id_soot
      use omp_lib
      implicit none
      character(*),optional,intent(in) :: read_from_file,write_to_file
      integer,intent(in),optional :: band_index
      integer :: nbck_unit
      logical :: rff,single_band,wtf
      integer :: ifv,ikg,inb,itg,its,ixc,ixw
      integer :: nfv,nkg,nnb,ntg,nts,nxc,nxw
      integer,allocatable,dimension(:,:) :: its_array
      real(dp) :: g,tg,xs(3)
      real(dp),allocatable,dimension(:,:) :: ts_array
      
      !-----------------------------------------------------------------
      !If the premixed NBCK has already been computed, skip everything
      !-----------------------------------------------------------------
!      if (nbck_is_premixed) return
!Commented because I don't know if it is a good idea to have this flag
      
      !-----------------------------------------------------------------
      !Set up optional parameters
      !-----------------------------------------------------------------
      call dprint('premix_nbck_data: Set up optional parameters')
      single_band = .false.; if (present(band_index)) single_band = .true.
      rff = .false.; if (present(read_from_file)) rff = .true.
      wtf = .false.; if (present(write_to_file))  wtf = .true.
      if (rff.and.wtf) rff = .false.
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('premix_nbck_data: Surrogate names')
      nkg = nbck_mix_nkg
      ntg = nbck_mix_ntg
      nxw = nbck_mix_nxw
      nxc = nbck_mix_nxc
      nfv = nbck_mix_nfv
      nnb = number_nbck_bands
      
      !-----------------------------------------------------------------
      !Set up thermodynamic states arrays
      !-----------------------------------------------------------------
      nts = ntg*nxw*nxc*nfv
      allocate(its_array(nts,4))
      allocate(ts_array(nts,4))
      its = 0
      do itg=1,ntg
         tg = nbck_mix_tg(itg)
         do ixw=1,nxw
            xs(id_h2o) = nbck_mix_xw(ixw)
            do ixc=1,nxc
               xs(id_co2) = nbck_mix_xc(ixc)
               do ifv=1,nfv
                  xs(id_soot) = nbck_mix_fv(ifv)
                  its = its + 1
                  its_array(its,1:4) = (/ itg, ixw, ixc, ifv /)
                  ts_array(its,1:4) = (/ tg, xs(id_h2o), xs(id_co2), &
                                         xs(id_soot) /)
               enddo
            enddo
         enddo
      enddo
      
      !-----------------------------------------------------------------
      !Main calculation
      !-----------------------------------------------------------------
      call dprint('premix_nbck_data: Main calculation')
      if (single_band) then
         if (rff) then
            nbck_unit = get_file_unit()
            open(unit=nbck_unit,file=read_from_file,&
               action='read',form='unformatted')
            read(nbck_unit) nbck_mix_k(:,:,:,:,:,band_index)
            close(nbck_unit)
         else
            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP&         PRIVATE(tg,xs,g,itg,ixw,ixc,ikg)
            do its=1,nts
               tg = ts_array(its,1); itg = its_array(its,1)
               xs(id_h2o) = ts_array(its,2); ixw = its_array(its,2)
               xs(id_co2) = ts_array(its,3); ixc = its_array(its,3)
               xs(id_soot) = ts_array(its,4); ifv = its_array(its,4)
               do ikg=1,nkg
                  g = nbck_mix_g(ikg)
                  nbck_mix_k(ikg,itg,ixw,ixc,ifv,band_index) = &
                     nbck_mix_kg(band_index,Tg,xs,1._dp,g,'k','riazzi')
               enddo
            enddo
            !$OMP ENDPARALLEL DO
         
            !Write to file
            if (wtf) then
               nbck_unit = get_file_unit()
               open(unit=nbck_unit,file=write_to_file,&
                  action='write',form='unformatted')
               write(nbck_unit) nbck_mix_k(:,:,:,:,:,band_index)
            endif
            close(nbck_unit)
         endif
      
      else
         if (rff) then
            nbck_unit = get_file_unit()
            open(unit=nbck_unit,file=read_from_file,&
               action='read',form='unformatted')
            read(nbck_unit) nbck_mix_k
            close(nbck_unit)
         else
            !$OMP PARALLEL DO DEFAULT(SHARED) &
            !$OMP&         PRIVATE(tg,xs,g,itg,ixw,ixc,ikg,inb)
            do its=1,nts
               tg = ts_array(its,1); itg = its_array(its,1)
               xs(id_h2o) = ts_array(its,2); ixw = its_array(its,2)
               xs(id_co2) = ts_array(its,3); ixc = its_array(its,3)
               xs(id_soot) = ts_array(its,4); ifv = its_array(its,3)
               do inb=1,nnb
                  do ikg=1,nkg
                     g = nbck_mix_g(ikg)
                     nbck_mix_k(ikg,itg,ixw,ixc,ifv,inb) = &
                        nbck_mix_kg(inb,Tg,xs,1._dp,g,'k','riazzi')
!if (inb.eq.83) & 
!   write(*,*) inb,tg,xs(id_h2o),xs(id_co2),g,nbck_mix_k(ikg,itg,ixw,ixc,inb)
                  enddo
               enddo
            enddo
            !$OMP ENDPARALLEL DO
!         stop
            !Write to file
            if (wtf) then
               nbck_unit = get_file_unit()
               open(unit=nbck_unit,file=write_to_file,&
                  action='write',form='unformatted')
               write(nbck_unit) nbck_mix_k
            endif
            close(nbck_unit)
         endif
      endif

      !-----------------------------------------------------------------
      !Set flag to indicate that the NBCK has been premixed
      !-----------------------------------------------------------------
      nbck_is_premixed = .true.
      
   endsubroutine premix_nbck_data

   !====================================================================
   !Subroutine to generate the quadrature points and weights
   !====================================================================
   subroutine get_ck_quadrature(which_quad,x,w)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use math_functions, only: quad_gauss_chebyshev
      implicit none
      character(*),intent(in) :: which_quad
      real(dp),intent(out) :: x(:),w(:)
      
      !-----------------------------------------------------------------
      !Select the appropriate quadrature
      !-----------------------------------------------------------------
      selectcase(trim(which_quad))
         case('single-node')
            x(1) = 0.5_dp; w(1) = 1._dp
         case('GC-even')
            call quad_gauss_chebyshev(size(x),'even-rank',x,w)
         case('GC-odd')
            call quad_gauss_chebyshev(size(x),'odd-rank',x,w)
         case default
            call shutdown('get_ck_quadrature: problem in the definition&
                           & of the quadrature type')
      endselect
      
   endsubroutine get_ck_quadrature

   !====================================================================
   !Function to interpolate the NBCK distribution for a single species
   !at a specific thermodynamic state. The interpolation can be made
   !for g from k (out_id = 'g',default) or for k from g (out_id = 'k')
   !====================================================================
   real(dp) function nbck_ss_kg(inb,isp,tmp,molfrac,pres,&
                                input_val,out_id)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use math_functions, only: locate
      use precision_parameters, only: small
      implicit none
      character(*),intent(in),optional :: out_id
      integer,intent(in) :: inb,isp
      integer :: ikg,itg,ixs,nkg,ntg,nxs
      integer :: lkg,lxs,ltg,ukg,uxs,utg
      logical :: interpolate_k,interpolate_g
      real(dp),intent(in) :: input_val,molfrac,pres,tmp
      real(dp) :: Rkg,Rtg,Rxs,Qkg,Qtg,Qxs
      real(dp) :: xl,xu,xval

      !-----------------------------------------------------------------
      !Set optional flag
      !-----------------------------------------------------------------
      interpolate_k = .false.; interpolate_g = .false.
      if (present(out_id)) then
         if (out_id.eq.'k') interpolate_k = .true.
         if (out_id.eq.'g') interpolate_g = .true.
      else
         interpolate_g = .true.
      endif

      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      ntg = nbck_data_ntg(isp)
      nxs = nbck_data_nxs(isp)

      !-----------------------------------------------------------------
      !Upper and lower indexes
      !-----------------------------------------------------------------
      !For temperature
      ltg = locate(nbck_data_tg(:,isp),tmp,ntg)                         !Locate lower index
      ltg = max(1,min(ntg-1,ltg)); utg = min(ltg+1,ntg)                 !Correct lower index, compute upper index
      if (tmp.lt.nbck_data_tg(1,isp))    utg = ltg                      !This is to prevent extrapolations
      if (tmp.gt.nbck_data_tg(ntg,isp))  ltg = utg                      !  (instead, simply take the value at the extreme)

      !For mole fraction
      lxs = locate(nbck_data_xs(:,isp),molfrac,nxs)                     !Locate lower index
      lxs = max(1,min(nxs-1,lxs)); uxs = min(lxs+1,nxs)                 !Correct lower index, compute upper index
      if (molfrac.lt.nbck_data_xs(1,isp))    uxs = lxs                  !This is to prevent extrapolations
      if (molfrac.gt.nbck_data_xs(nxs,isp))  lxs = uxs                  !  (instead, simply take the value at the extreme)

      !-----------------------------------------------------------------
      !Interpolate
      !-----------------------------------------------------------------
      !Initial values for the interpolation on mole fraction
      Qxs = 0._dp
      Rxs = (molfrac - nbck_data_xs(lxs,isp))/&
         (nbck_data_xs(uxs,isp) - nbck_data_xs(lxs,isp) + small)
      if (lxs.eq.uxs) Rxs = 0._dp                                       !If only one mole fraction value is provided,
                                                                        !  do not interpolate in mole fraction
      molfrac_loop: do ixs=lxs,uxs
         !Initial values for the interpolation on temperature
         Qtg = 0._dp
         Rtg = (tmp - nbck_data_tg(ltg,isp))/&
            (nbck_data_tg(utg,isp) - nbck_data_tg(ltg,isp) + small)
         if (ltg.eq.utg) Rtg = 0._dp                                    !If only one temperature value is provided, 
                                                                        !  do not interpolate in temperature

         tmp_loop: do itg=ltg,utg
            !Upper and lower indexes for k at 
            !the current thermodynamic state
            nkg = nbck_data_nkg(isp)
            if (interpolate_k) then
               xval = input_val
               lkg = locate(nbck_g(:,isp),xval,nkg)
               lkg = max(1,min(nkg-1,lkg)); ukg = min(lkg+1,nkg)
               if (xval.lt.nbck_g(1,isp))    ukg = lkg
               if (xval.gt.nbck_g(nkg,isp))  lkg = ukg
               xl = nbck_g(lkg,isp)
               xu = nbck_g(ukg,isp)
            endif
            if (interpolate_g) then
               xval = input_val/(molfrac*pres + small)
               lkg = locate(nbck_k(:,itg,ixs,inb,isp),xval,nkg)
               lkg = max(1,min(nkg-1,lkg)); ukg = min(lkg+1,nkg)
               if (xval.lt.nbck_k(1,itg,ixs,inb,isp))    ukg = lkg
               if (xval.gt.nbck_k(nkg,itg,ixs,inb,isp))  lkg = ukg
               xl = nbck_k(lkg,itg,ixs,inb,isp)
               xu = nbck_k(ukg,itg,ixs,inb,isp)
            endif

            !Initial values for the interpolation on k
            Qkg = 0._dp; Rkg = (xval - xl)/(xu - xl + small)
            if (lkg.eq.ukg) Rkg = 0._dp                                 !If only one k value is provided, 
                                                                        !  do not interpolate in k
!write(*,*) nbck_data_xs(ixs,isp),nbck_data_tg(itg,isp)
!do ikg=1,nkg
!write(*,*) isp,ikg,nbck_g(ikg,isp),nbck_k(ikg,itg,ixs,inb,isp)
!enddo
            kg_loop: do ikg=lkg,ukg
               !Interpolate on k or g
               Rkg = 1._dp - Rkg
               if (interpolate_k) &
                  Qkg = Qkg + Rkg*nbck_k(ikg,itg,ixs,inb,isp)
               if (interpolate_g) &
                  Qkg = Qkg + Rkg*nbck_g(ikg,isp)
            enddo kg_loop
            
            !Interpolate on temperature
            Rtg = 1._dp - Rtg
            Qtg = Qtg + Rtg*Qkg
         enddo tmp_loop

         !Interpolate on mole fraction
         Rxs = 1._dp - Rxs
         Qxs = Qxs + Rxs*Qtg
      enddo molfrac_loop

      !-----------------------------------------------------------------
      !Final value
      !-----------------------------------------------------------------
      if (interpolate_k) Qxs = Qxs*molfrac*pres
      nbck_ss_kg = Qxs

   endfunction nbck_ss_kg

   !====================================================================
   !Function to compute the mixture CK distribution
   !====================================================================
   real(dp) recursive function nbck_mix_kg(bindex,TT,xxs,pp,kg_in,&
                                        out_id,mixing_method) result(kg)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use global_parameters, only: id_h2o,id_co2,id_soot
      use math_functions, only: heaviside_dp,locate
      use precision_parameters, only: big,dp,small
      implicit none
      character(*),intent(in) :: out_id
      character(*),intent(in),optional :: mixing_method
      character(200) :: which_mixing
      integer,intent(in) :: bindex
      integer,parameter :: max_iter = 1000                              !Maximum number of iterations
      integer :: counter,iqd,jqd,nqd,isp,nsp,single_spec_id
      integer :: ifv,ikg,itg,ixc,ixw,nfv,nkg,ntg,nxc,nxw
      integer :: lfv,lkg,ltg,lxc,lxw,ufv,ukg,utg,uxc,uxw
      real(dp),intent(in) :: kg_in,pp,TT,xxs(:)
      real(dp),parameter :: ktol = 1.e-3_dp,gtol = 1.e-3_dp             !Tolerances      
      real(dp) :: diff,gmix,gmix0,H,ki,kj,spec_cutoff
      real(dp) :: kleft,kright,ktry
      real(dp) :: gleft,gright,gtry
      real(dp) :: fv,xc,xw
      real(dp) :: Qfv,Qkg,Qtg,Qxc,Qxw,Rfv,Rkg,Rtg,Rxc,Rxw
      real(dp),allocatable,dimension(:) :: xquad,wquad
      kg = 0._dp                                                        !This is here to avoid compilation warnings

      !Set up optional parameter
      which_mixing = nbck_mixing_method
      if (present(mixing_method)) which_mixing = mixing_method      

      !-----------------------------------------------------------------
      !Special case for premixed NBCK databases (for now, only works for
      !CO2-H2O mixtures and for obtaining k, not g)
      !-----------------------------------------------------------------
   
      if (trim(which_mixing).eq.'premixed') then
         !Upper and lower temperature indexes
         ntg = nbck_mix_ntg
         ltg = locate(nbck_mix_tg,TT,ntg)
         ltg = max(1,min(ntg-1,ltg)); utg = min(ltg+1,ntg)
         if (TT.lt.nbck_mix_tg(1))    utg = ltg
         if (TT.gt.nbck_mix_tg(ntg))  ltg = utg
!do itg=1,nbck_mix_ntg
!write(*,*) 'T',itg,nbck_mix_tg(itg)
!enddo
!write(*,*) 'T',TT,nbck_mix_tg(ltg),nbck_mix_tg(utg)

         !Upper and lower H2O mole fraction indexes
         xw = xxs(id_h2o)
         nxw = nbck_mix_nxw
         lxw = locate(nbck_mix_xw,xw,nxw)                               !Locate lower index
         lxw = max(1,min(nxw-1,lxw)); uxw = min(lxw+1,nxw)              !Correct lower index, compute upper index
         if (xw.lt.nbck_mix_xw(1))   uxw = lxw                          !This is to prevent extrapolations
         if (xw.gt.nbck_mix_xw(nxw)) lxw = uxw                          !  (instead, simply take the value at the extreme)

!do ixw=1,nbck_mix_nxw
!write(*,*) 'xw',ixw,nbck_mix_xw(ixw)
!enddo
!write(*,*) 'xw',xw,nbck_mix_xw(lxw),nbck_mix_xw(uxw)

         !Upper and lower CO2 mole fraction indexes
         xc = xxs(id_co2)
         nxc = nbck_mix_nxc
         lxc = locate(nbck_mix_xc,xc,nxc)                               !Locate lower index
         lxc = max(1,min(nxc-1,lxc)); uxc = min(lxc+1,nxc)              !Correct lower index, compute upper index
         if (xc.lt.nbck_mix_xc(1))    uxc = lxc                         !This is to prevent extrapolations
         if (xc.gt.nbck_mix_xc(nxc))  lxc = uxc                         !  (instead, simply take the value at the extreme)


!do ixc=1,nbck_mix_nxc
!write(*,*) 'xc',ixc,nbck_mix_xw(ixc)
!enddo
!write(*,*) 'xc',xc,nbck_mix_xc(lxc),nbck_mix_xc(uxc)

         !Upper and lower soot volume fraction indexes
         fv = xxs(id_soot)
         nfv = nbck_mix_nfv
         lfv = locate(nbck_mix_fv,fv,nfv)                               !Locate lower index
         lfv = max(1,min(nfv-1,lfv)); ufv = min(lfv+1,nfv)              !Correct lower index, compute upper index
         if (fv.lt.nbck_mix_fv(1))    ufv = lfv                         !This is to prevent extrapolations
         if (fv.gt.nbck_mix_fv(nfv))  lfv = ufv                         !  (instead, simply take the value at the extreme)

!do ifv=1,nbck_mix_nfv
!write(*,*) 'fv',ifv,nbck_mix_fv(ifv)
!enddo
!write(*,*) 'fv',fv,nbck_mix_fv(lfv),nbck_mix_fv(ufv)

         !Upper and lower g indexes
         nkg = nbck_mix_nkg
         lkg = locate(nbck_mix_g,kg_in,nkg)
         lkg = max(1,min(nkg-1,lkg)); ukg = min(lkg+1,nkg)
         if (kg_in.lt.nbck_mix_g(1))    ukg = lkg
         if (kg_in.gt.nbck_mix_g(nkg))  lkg = ukg


!do ikg=1,nbck_mix_nkg
!write(*,*) 'g',ikg,nbck_mix_g(ikg)
!enddo
!write(*,*) 'g',kg_in,nbck_mix_g(lkg),nbck_mix_g(ukg)

         !Interpolation on temperature
         Qtg = 0._dp
         Rtg = (TT - nbck_mix_tg(ltg))/&
            (nbck_mix_tg(utg) - nbck_mix_tg(ltg) + small)
         if (ltg.eq.utg) Rtg = 0._dp  
         tg_loop: do itg=ltg,utg

            !Interpolation in CO2
            Qxc = 0._dp
            Rxc = (xc - nbck_mix_xc(lxc))/&
               (nbck_mix_xc(uxc) - nbck_mix_xc(lxc) + small)
            if (lxc.eq.uxc) Rxc = 0._dp
            xc_loop: do ixc=lxc,uxc
 
               !Interpolation in H2O
               Qxw = 0._dp
               Rxw = (xw - nbck_mix_xw(lxw))/&
                  (nbck_mix_xw(uxw) - nbck_mix_xw(lxw) + small)
               if (lxw.eq.uxw) Rxw = 0._dp
               xw_loop: do ixw=lxw,uxw

                  !Interpolation in soot
                  Qfv = 0._dp
                  Rfv = (fv - nbck_mix_fv(lfv))/&
                     (nbck_mix_fv(ufv) - nbck_mix_fv(lfv) + small)
                  if (lfv.eq.ufv) Rfv = 0._dp

                  fv_loop: do ifv=lfv,ufv
                     !Interpolation in g 
                     Qkg = 0._dp
                     Rkg = (kg_in - nbck_mix_g(lkg))/&
                        (nbck_mix_g(ukg) - nbck_mix_g(lkg) + small)
                     if (lkg.eq.ukg) Rkg = 0._dp
                     kg_loop: do ikg=lkg,ukg
                        Rkg = 1._dp - Rkg
                        Qkg = Qkg + Rkg*&
                           nbck_mix_k(ikg,itg,ixw,ixc,ifv,bindex)
!write(*,*) itg,ixc,ixw,ifv,ikg,nbck_mix_k(ikg,itg,ixw,ixc,ifv,bindex)
                     enddo kg_loop
!write(*,*) 'g',Qkg
                     Rfv = 1._dp - Rfv
                     Qfv = Qfv + Rfv*Qkg
                  enddo fv_loop
!write(*,*) 'fv',Qfv
                  Rxw = 1._dp - Rxw
                  Qxw = Qxw + Rxw*Qfv                     
               enddo xw_loop
!write(*,*) 'xc',Qxw
               Rxc = 1._dp - Rxc
               Qxc = Qxc + Rxc*Qxw
            enddo xc_loop
!write(*,*) 'xw',Qxc
            Rtg = 1._dp - Rtg
            Qtg = Qtg + Rtg*Qxc
         enddo tg_loop
         kg = Qtg!*pp
         return
      endif

      !-----------------------------------------------------------------
      !Special case for a medium with no or only 
      !one participating species
      !-----------------------------------------------------------------
      !Determine how many participating species with 
      !non-negligible concentrations exist
      spec_cutoff = 0._dp                                               !Cutoff: any mole fraction below this value 
                                                                        !  is considered negligible
      nsp = size(xxs); single_spec_id = -1; counter = 0                 !Some initial parameters
      do isp=1,nsp
         if (xxs(isp).gt.spec_cutoff) then
            single_spec_id = isp; counter = counter + 1
         endif
      enddo

      !No participating species
      if (counter.eq.0) then
         if (out_id.eq.'k') kg = 0._dp
         if (out_id.eq.'g') kg = 1._dp
         return
      endif

      !Only one participating species
      if (counter.eq.1) then
         if (out_id.eq.'k') &
            kg = nbck_ss_kg(bindex,single_spec_id,TT,&
                            xxs(single_spec_id),pp,kg_in,'k')
         if (out_id.eq.'g') &
            kg = nbck_ss_kg(bindex,single_spec_id,TT,&
                            xxs(single_spec_id),pp,kg_in,'g')
         return
      endif

      !-----------------------------------------------------------------
      !Compute g for the mixture
      !-----------------------------------------------------------------
      if (trim(out_id).eq.'g') then      
         selectcase(which_mixing)              
         case('riazzi')
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !Use the mixing model of Modest and Riazzi, 2005
            !(only two species are implemented so far; more can probably
            !be implemented recursively)
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            nqd = 32                                                    !Number of quadrature points for the integration
            allocate(xquad(nqd)); allocate(wquad(nqd))                  !Allocate arrays with quadrature points and weights
            gmix = 0._dp                                                !Zero out gmix for the integration
            call get_ck_quadrature('GC-even',xquad,wquad)               !Get the quadrature points and weights
            do iqd=1,nqd
               ki = nbck_ss_kg(bindex,1,TT,xxs(1),pp,xquad(iqd),'k')    !k value for species 1
               do jqd=1,nqd
                  kj = nbck_ss_kg(bindex,2,TT,xxs(2),pp,xquad(jqd),'k') !k value for species 2
                  H = heaviside_dp(kg_in - ki - kj)                     !Apply the heaviside function              
                  gmix = gmix + H*wquad(iqd)*wquad(jqd)                 !Integrate for gmix
               enddo
            enddo
            deallocate(xquad,wquad)

         case('riazzi-iterative')
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !Use the mixing model of Modest and Riazzi, 2005, with the
            !addition of an iterative scheme for the integrations, in
            !order to increase the accuracy (increases cost as well)
            !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            nqd = 16; gmix0 = 0._dp; diff = 2._dp*gtol                  !Initial parameters
            do while (diff.gt.gtol)
               gmix = 0._dp                                             !Zero out gmix for the integration
               allocate(xquad(nqd))                                     !Allocate arrays with quadrature points and
               allocate(wquad(nqd))                                     !  weights according to the current number
                                                                        !  of quadrature points
               call get_ck_quadrature('GC-even',xquad,wquad)            !Get the quadrature points and weights
               do iqd=1,nqd
                  ki = nbck_ss_kg(bindex,1,TT,xxs(1),pp,xquad(iqd),'k') !k value for species 1
                  do jqd=1,nqd
                     kj = nbck_ss_kg(bindex,2,TT,xxs(2),pp,&            !k value for species 2
                                     xquad(jqd),'k')
                     H = heaviside_dp(kg_in - ki - kj)                  !Apply the heaviside function              
                     gmix = gmix + H*wquad(iqd)*wquad(jqd)              !Integrate for gmix
                  enddo
               enddo
               diff = dabs((gmix - gmix0)/(gmix + small))               !Compute difference relative to the previous interation
               gmix0 = gmix; nqd = 2*nqd                                !Update parameters for the next iteration
               deallocate(xquad,wquad)                                  !Deallocate arrays that will be reallocated    
            enddo
         case default
            call shutdown('nbck_mix_kg: problem with the specification &
                           &of nbck_mixture_method')
         endselect
         
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Final value
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         kg = gmix
         return
      
      !-----------------------------------------------------------------
      !Compute k for the mixture 
      !(for a mixture, this is done iteratively)
      !-----------------------------------------------------------------
      elseif (out_id.eq.'k') then 

         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Initial parameters for the bisection method
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !k and g to the left and to the right
         gleft = 0._dp; gright = 1._dp
         kleft = big; kright = -big
         do isp=1,nsp
            kleft = min(kleft,&
                        nbck_ss_kg(bindex,isp,TT,xxs(isp),pp,gleft,'k'))
            kright = max(kright,&
                         nbck_ss_kg(bindex,isp,TT,xxs(isp),pp,gright,'k'))
         enddo
 
         !Counter for the number of iterations
         counter = 0
      
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Apply the bisection method to find k
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         inversion_loop: do
            !Update counter
            counter = counter + 1
      
            !Define ktry
            ktry = 0.5_dp*(kleft+kright)
         
            !Compute gtry from ktry
            gtry = nbck_mix_kg(bindex,TT,xxs,pp,ktry,'g',which_mixing)

            !Check convergence
            diff = (kg_in - gtry)/(kg_in + small)
            if ((dabs(diff).le.gtol).or.&
               (dabs(kg_in - gtry).le.gtol)) exit inversion_loop
         
            !Escape if max iterations exceeded
            if (counter.gt.max_iter) then
               if ((dabs(kright-kleft).lt.ktol)) exit inversion_loop    !Temporary        
               call shutdown('nbck_mix_kg: Counter exceeded')
            endif
         
            !Update boundary values
            if ((gtry - kg_in)/(gleft - kg_in).lt.0) then
               kright = ktry; gright = gtry
            elseif ((gtry - kg_in)/(gright - kg_in).lt.0) then
               kleft = ktry; gleft = gtry
            else
               call shutdown('nbck_mix_kg: problem in the inversion')
            endif
         enddo inversion_loop     
      
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         !Final value
         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         kg = ktry
         return
      
      else
         call shutdown('nbck_mix_kg: problem with the specification &
                       &of out_id')
      endif
      
   endfunction nbck_mix_kg

   !====================================================================
   !Subroutine to compute the transmissivity according to the NBCK
   !model for narrow band of index band_index
   !====================================================================
   real(dp) function nbck_band_transmissivity(x,T,xs,p,inb,nqd)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: assert,dprint
      use global_parameters, only: number_of_species
      use omp_lib
      implicit none
      integer,intent(in) :: inb,nqd
      integer :: i,iqd,nx,nsp
      real(dp),intent(in):: x(:),p(:),T(:),xs(:,:)
      real(dp) :: xquad(nqd),wquad(nqd)
      real(dp) :: dx,int_omp,int_quad,k,g,ot
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('nbck_band_transmissivity: Preparatory procedures')
      
      !Check the consistency of the input arrays
      call assert(size(x).eq.size(T),'size(x) = size(T)')
      call assert(size(x).eq.size(T),'size(x) = size(p)')
      call assert(size(x).eq.size(xs,2),'size(x) = size(xs,2)')
      
      !Surrogate names
      nx = size(x)
      nsp = number_of_species

      !-----------------------------------------------------------------
      !Load NBCK data
      !-----------------------------------------------------------------
      call dprint('nbck_band_transmissivity: Load NBCK data')
      call load_nbck_data

      !-----------------------------------------------------------------
      !Generate quadrature points and weights
      !-----------------------------------------------------------------
      call dprint('nbck_band_transmissivity: &
                               &Generate quadrature points and weights')
      call get_ck_quadrature(nbck_quadrature,xquad,wquad)
      
      !-----------------------------------------------------------------
      !Main calculation
      !-----------------------------------------------------------------
      call dprint('nbck_band_transmissivity: Main calculation')
      int_quad = 0._dp
      !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(dx,g,i,int_omp,k,ot)
      int_omp = 0._dp
      
      !$OMP DO   
      quadrature_loop: do iqd=1,nqd
         g = xquad(iqd)
         ot = 0._dp
         do i=1,nx
            !Compute k
            k = nbck_mix_kg(inb,T(i),xs(:,i),p(i),g,'k')
!write(*,*) iqd,g,k
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
            
            !Integrate for the optical thickness
            ot = ot + k*dx
         enddo
         
         !Integrate for the transmissivity
         int_omp = int_omp + dexp(-ot)*wquad(iqd)
      enddo quadrature_loop
      !$OMP ENDDO
      
      !-----------------------------------------------------------------
      !Finish the integration by summing up over all OMP threads
      !-----------------------------------------------------------------
      !$OMP CRITICAL
      int_quad = int_quad + int_omp
      !$OMP ENDCRITICAL
      !$OMP ENDPARALLEL

      !-----------------------------------------------------------------
      !Final value
      !-----------------------------------------------------------------
      nbck_band_transmissivity = int_quad

   endfunction nbck_band_transmissivity
   
   !====================================================================
   !Function to compute the scaling function in the NBCK model, 
   !following the formulation of Modest & Zhang, 2002 (JHT)
   !====================================================================
!   real(dp) function nbck_scaling_function(TT,xxs,LL,g,w,inb)
   
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use precision_parameters, only: dp,small
!      use comp_functions, only: shutdown
!      implicit none
!      integer :: counter,iqd,nqd
!      integer,parameter :: max_iter = 10000
!      integer,intent(in) :: inb
!      real(dp),intent(in) :: g(:),LL,TT,xxs(:),w(:)
!      real(dp),parameter :: tol = 1.e-3_dp
!      real(dp) :: k_0(size(g)),k_l(size(g)),tau_0,tau_l
!      real(dp) :: root,root_m,root_p,u,u_m,u_p
      
!      !-----------------------------------------------------------------
!      !Preparatory procedures   
!      !-----------------------------------------------------------------
!      !Number of quadrature points
!      nqd = size(g)
      
!      !Compute k for each quadrature point at 
!      !the local and reference conditions
!      do iqd=1,nqd
!         k_l(iqd) = nbck_mix_kg(inb,TT,xxs,1._dp,g(iqd),'k')
!         k_0(iqd) = nbck_mix_kg(inb,nbck_T0,nbck_xs0,1._dp,g(iqd),'k')
!      enddo
   
!      !Compute gas column transmissitivy at the local conditions
!      tau_l = 0._dp
!      do iqd=1,nqd
!         tau_l = tau_l + dexp(-k_l(iqd)*LL)*w(iqd)
!      enddo
      
!      !-----------------------------------------------------------------
!      !Prepare for the iteration for u
!      !-----------------------------------------------------------------
!      !u bound and root to the left
!      u_m = 0._dp; root_m = 1._dp
      
!      !Determine u bound and root to the right
!      root_p = 1._dp; u = 1._dp
!      do while (root_p.gt.1._dp)
!         u = u*1.5_dp; tau_0 = 0._dp        
!         do iqd=1,nqd
!            tau_0 = tau_0 + dexp(-u*k_0(iqd)*LL)*w(iqd)
!         enddo
!         root_p = (tau_l - tau_0)/(tau_l + small)
!      enddo
!      u_p = u
      
!      !-----------------------------------------------------------------
!      !Determine u using the bisection method
!      !-----------------------------------------------------------------
!      counter = 0
!      bisection_loop: do
!         !Update counter
!         counter = counter + 1
      
!         !Define try value for u
!         u = 0.5_dp*(u_p + u_m)
         
!         !Compute transmissivity at the reference condition
!         tau_0 = 0._dp        
!         do iqd=1,nqd
!            tau_0 = tau_0 + dexp(-u*k_0(iqd)*LL)*w(iqd)
!         enddo
         
!         !Compute root
!         root = (tau_l - tau_0)/(tau_l + small)
         
!         !Check convergence
!         if (dabs(root).le.tol) exit bisection_loop
         
!         !Escape if max iterations exceeded
!         if (counter.gt.max_iter) then
!            if ((dabs(u_p-u_m).lt.tol)) exit bisection_loop             !Temporary        
!            call shutdown('nbck_scaling_function: Counter exceeded')
!         endif
         
!         !Update boundary values
!         if ((root/root_p).gt.0) then
!            u_p = u; root_p = root
!         elseif ((root/root_m).gt.0) then
!            u_m = u; root_m = root
!         else
!            call shutdown('nbck_scaling_function: &
!                          &problem in the inversion')
!         endif
!      enddo bisection_loop   
      
!      !Final value      
!      nbck_scaling_function = u

!   endfunction nbck_scaling_function


   !====================================================================
   !Subroutine that generates narrow-band cumulative-k distributions 
   !from a LBL database for a set of thermodynamic states, specified
   !in terms of a temperature (array T_in) and chemical species
   !concentrations (array xs_in); the pressure is assumed to be 1 atm.
   !The NBCK is generated for the discrete k values given in array k_in.
   !The resulting NBCK values is saved into array g_out (rank 2)
   !====================================================================
!   subroutine get_local_nbck(T_in,xs_in,k_in,lower_bound,upper_bound,&
!                             g_out)
   
!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: assert,CheckMemAlloc,print_to_prompt,&
!                                shutdown
!      use global_parameters, only: debug_mode
!      use lbl_functions, only: close_lbl_data,lbl_kappa_func,&
!                               open_lbl_data,real_lbl_data,&
!                               rewind_lbl_data
!      use lbl_parameters, only: lbl_data,lbl_data_averaging,&
!         lbl_data_cm,lbl_data_ready,lbl_ns,lbl_nx,lbl_nt
!      use nbck_parameters
!      use precision_parameters, only: dp
!      implicit none
!      integer :: ikg,isp,its
!      integer :: nkg,nsp,nts
!      real(dp),intent(in) :: lower_bound,upper_bound
!      real(dp),intent(in) :: k_in(:),T_in(:),xs_in(:,:)
!      real(dp),intent(out) :: g_out(:,:)
!      real(dp) :: deta,eta_down,eta_factor,eta_up,nb_deta,xeta
!      real(dp) :: kappa
!      real(dp) :: acs_down(lbl_ns,lbl_nx,lbl_nt),&
!                  acs_up(lbl_ns,lbl_nx,lbl_nt)

!      !-----------------------------------------------------------------
!      !Set up initial parameters
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('get_local_nbck: &
!                                            &set up initial parameters')
      
!      !Check for consistency in the array sizes
!      call assert(size(T_in).eq.size(xs_in,1))
!      call assert(size(T_in).eq.size(g_out,2))
!      call assert(size(k_in).eq.size(g_out,1))
      
!      !-----------------------------------------------------------------
!      !Surrogate names
!      !-----------------------------------------------------------------
!      nkg = size(k_in)
!      nsp = size(xs_in,2)
!      nts = size(T_in)
      
!      !-----------------------------------------------------------------
!      !Prepare the LBL data
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('get_local_nbck: &
!                                           &prepare LBL data files')
!      call open_lbl_data
!      if (lbl_data_ready) call rewind_lbl_data
!      if (lbl_data_averaging.ne.'none') &
!         call real_lbl_data(eta_down,acs_down)                          !Read the first line of the data
!      eta_factor = 1._dp; if (lbl_data_cm) eta_factor = 100._dp         !Factor to convert from 1/cm to 1/m, if needed

!      !-----------------------------------------------------------------
!      !Main calculation
!      !-----------------------------------------------------------------
!      if (debug_mode) call print_to_prompt('get_local_nbck: &
!                                           &Main calculation')
!      g_out = 0._dp; nb_deta = 0._dp
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
!         if (xeta.lt.lower_bound) then
!            eta_down = eta_up; acs_down = acs_up
!            cycle lbl_loop
!         elseif (xeta.gt.upper_bound) then
!            exit lbl_loop
!         endif
         
!         !Add to the narrow band width
!         nb_deta = nb_deta + deta
         
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
!         !Compute the NBCKs
!         !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!         state_loop: do its=1,nts
!            !Compute the absorption coefficient
!            kappa = 0._dp
!            do isp=1,nsp
!               kappa = kappa + &
!                  lbl_kappa_func(xs_in(its,isp),T_in(its),1._dp,isp,&
!                                 kappa_p=.false.)
!            enddo
!write(*,*) 'k',kappa
!            !Sweep through all k values and determine g for each of them
!            do ikg=1,nkg
!               if (kappa.le.k_in(ikg)) &
!                  g_out(ikg,its) = g_out(ikg,its) + deta
!            enddo
         
!         enddo state_loop
             
!         !Update spectral values
!         eta_down = eta_up                                              !Wavenumber
!         acs_down = acs_up                                              !Absorption cross-section

!      enddo lbl_loop
!      call close_lbl_data                                               !Close LBL units

!      !Finish the calculation of g
!      g_out = g_out/nb_deta

!   endsubroutine get_local_nbck

endmodule nbck_functions
