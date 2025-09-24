module slw_functions

   !====================================================================
   !Modules & Misc
   !====================================================================
   use precision_parameters, only: dp
   use slw_parameters
   implicit none
   
   contains   
   
   !====================================================================
   !Subroutine to read the sizes of the ALBDF-related arrays from an
   !external information file
   !====================================================================
   subroutine read_albdf_info_1
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,dprint,get_file_unit,&
                                shutdown
      use global_parameters, only: number_of_species
      implicit none
      integer :: in_unit,isp
      integer :: nTg(number_of_species),nTb(number_of_species),&
         nCj(number_of_species),npb(number_of_species)
      logical :: albdf_invert_flag(number_of_species)
      
      !-----------------------------------------------------------------
      !Read the parameters that define the ALBDF array sizes
      !and check if they are the same for all species
      !-----------------------------------------------------------------
      call dprint('read_albdf_info_1: Read the parameters that &
                  &define the ALBDF array sizes')
      do isp=1,number_of_species
         if (albdf_info_file(isp).eq.'null') cycle      
         in_unit = get_file_unit()
         call CheckFileExists(albdf_info_file(isp))
         open(unit=in_unit,file=albdf_info_file(isp),&
              form='unformatted',action='read')
         read(in_unit) albdf_invert_flag(isp)
         read(in_unit) albdf_nx(isp),nTg(isp),nTb(isp),nCj(isp),&
                       npb(isp)
      
         if (isp.gt.1) then
            if (nTg(isp).ne.(nTg(isp-1))) &
               call shutdown('read_albdf_info_1: &
                             &ALBDFs with different Tg sizes')
            if (nTb(isp).ne.(nTb(isp-1))) &
               call shutdown('read_albdf_info_1: &
                             &ALBDFs with different Tb sizes')
            if (nCj(isp).ne.(nCj(isp-1))) &
               call shutdown('read_albdf_info_1: &
                             &ALBDFs with different Cj sizes')
            if (npb(isp).ne.(npb(isp-1))) &
               call shutdown('read_albdf_info_1: &
                             &ALBDFs with different Pb sizes')
            if (xor(albdf_invert_flag(isp),albdf_invert_flag(isp-1))) &
               call shutdown('read_albdf_info_1: &
                             &ALBDFs invert flag mismatch')
         endif
         
         close(in_unit)
      enddo
      
      !-----------------------------------------------------------------
      !If no error occurred, assign the sizes of the ALBDF array
      !-----------------------------------------------------------------
      call dprint('read_albdf_info_1: &
                  &Assign the sizes of the ALBDF array')
      albdf_nTg = nTg(1)
      albdf_nTb = nTb(1)
      albdf_nCj = nCj(1)
      albdf_nbands = npb(1)
      albdf_inverted = albdf_invert_flag(1)
      
   endsubroutine read_albdf_info_1
   
   !====================================================================
   !Subroutine to allocate all arrays related to the ALBDF
   !====================================================================
   subroutine allocate_albdf_parameters
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint
      use global_parameters, only: number_of_species
      implicit none
      integer :: ns,nTg,nTb,nCj,npb
      integer :: ierr
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('allocate_albdf_parameters: Surrogate names')
      ns  = number_of_species
      nTg = albdf_nTg
      nTb = albdf_nTb
      nCj = albdf_nCj
      npb = albdf_nbands
   
      !-----------------------------------------------------------------
      !Deallocate arrays if they are already allocated
      !-----------------------------------------------------------------
      call dprint('allocate_albdf_parameters: &
                  &Deallocate arrays if they are already allocated')
      if (allocated(albdf_darr))    deallocate(albdf_darr)
      if (allocated(albdf_iarr))    deallocate(albdf_iarr)
      if (allocated(albdf_Tb))      deallocate(albdf_Tb)
      if (allocated(albdf_Tg))      deallocate(albdf_Tg)
      if (allocated(bslw_lbound))   deallocate(bslw_lbound)
      if (allocated(bslw_ubound))   deallocate(bslw_ubound)      
      
      !-----------------------------------------------------------------
      !Allocate arrays
      !-----------------------------------------------------------------
      call dprint('allocate_albdf_parameters: Allocate arrays')
      allocate(albdf_Tb(nTb),stat=ierr)
      call CheckMemAlloc('albdf_Tb',ierr)
      allocate(albdf_Tg(nTg),stat=ierr)
      call CheckMemAlloc('albdf_Tg',ierr)
      allocate(albdf_iarr(nCj),stat=ierr)
      call CheckMemAlloc('albdf_iarr',ierr)
      allocate(bslw_lbound(npb),stat=ierr)
      call CheckMemAlloc('bslw_lbound',ierr)
      allocate(bslw_ubound(npb),stat=ierr)
      call CheckMemAlloc('bslw_ubound',ierr)
      if (trim(slw_mixture_method).ne.'albdf_precombined') then
         allocate(albdf_darr(nCj,nTg,slw_nx,ns,nTb,npb),stat=ierr)
         call CheckMemAlloc('albdf_darr',ierr)
      else

      endif
      
   endsubroutine allocate_albdf_parameters

   !====================================================================
   !Subroutine to read the arrays needed for ALBDF interpolations
   !from an external information file
   !====================================================================
   subroutine read_albdf_info_2
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,CheckMemAlloc,dprint,&
                                get_file_unit,shutdown
      use global_parameters, only: number_of_species
      implicit none
      integer :: in_unit,ierr,ii,isp
      integer :: nTg,nTb,nCj,npb,nsp
      real(dp),allocatable,dimension(:) :: Tg_aux,Tb_aux,iarr_aux,&
                                           lpb_aux,upb_aux
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('read_albdf_info_2: Surrogate names')
      nsp = number_of_species
      nTg = albdf_nTg
      nTb = albdf_nTb
      nCj = albdf_nCj
      npb = albdf_nbands

      !-----------------------------------------------------------------
      !Allocate auxiliary arrays
      !-----------------------------------------------------------------
      call dprint('read_albdf_info_2: Allocate auxiliary arrays')
      allocate(iarr_aux(nCj),stat=ierr)
      call CheckMemAlloc('Cj_aux',ierr)
      allocate(lpb_aux(npb),stat=ierr)
      call CheckMemAlloc('lpb_aux',ierr)
      allocate(upb_aux(npb),stat=ierr)
      call CheckMemAlloc('upb_aux',ierr)
      allocate(Tb_aux(nTb),stat=ierr)
      call CheckMemAlloc('Tb_aux',ierr)
      allocate(Tg_aux(nTg),stat=ierr)
      call CheckMemAlloc('Tg_aux',ierr)
            
      !-----------------------------------------------------------------
      !Read the arrays needed for ALBDF interpolations
      !and check if they are the same for all species
      !-----------------------------------------------------------------
      call dprint('read_albdf_info_2: Read the arrays needed &
                  &for ALBDF interpolations')
      do isp=1,nsp
         if (albdf_info_file(isp).eq.'null') cycle
         in_unit = get_file_unit()
         call CheckFileExists(albdf_info_file(isp))
         open(unit=in_unit,file=albdf_info_file(isp),&
              form='unformatted',action='read')
         read(in_unit)
         read(in_unit) 
         read(in_unit) albdf_x(1:albdf_nx(isp),isp)

         if (isp.eq.1) then
            read(in_unit) albdf_iarr
            read(in_unit) albdf_Tb
            read(in_unit) albdf_Tg
            read(in_unit) bslw_lbound
            read(in_unit) bslw_ubound
         else
            read(in_unit) iarr_aux
            read(in_unit) Tb_aux
            read(in_unit) Tg_aux
            read(in_unit) lpb_aux
            read(in_unit) upb_aux
      
            !Check if Cj array is equal for all species
            do ii=1,nCj
               if (albdf_iarr(ii).ne.iarr_aux(ii)) &
                  call shutdown('read_albdf_info_2: &
                        &albdf_iarr mismatch between different ALBDFs')
            enddo
            
            !Check if Tb array is equal for all species
            do ii=1,nTb
               if (albdf_Tb(ii).ne.Tb_aux(ii)) &
                  call shutdown('read_albdf_info_2: &
                        &albdf_Tb mismatch between different ALBDFs')
            enddo
      
            !Check if Tg array is equal for all species
            do ii=1,nTg
               if (albdf_Tg(ii).ne.Tg_aux(ii)) &
                  call shutdown('read_albdf_info_2: &
                        &albdf_Tg mismatch between different ALBDFs')
            enddo
            
            !Check if bslw_lbound array is equal for all species
            if (npb.gt.1) then
               do ii=1,npb
                  if (bslw_lbound(ii).ne.lpb_aux(ii)) &
                  call shutdown('read_albdf_info_2: &
                        &bslw_lbound mismatch between different ALBDFs')
                  if (bslw_ubound(ii).ne.upb_aux(ii)) &
                  call shutdown('read_albdf_info_2: &
                        &bslw_ubound mismatch between different ALBDFs')
               enddo
            endif
         endif
         close(in_unit)
      enddo

   endsubroutine read_albdf_info_2
   
   !====================================================================
   !Subroutine to read the ALBDF from an external file
   !====================================================================
   subroutine read_aldbf
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckFileExists,dprint,get_file_unit,&
                                shutdown
      use global_parameters, only: number_of_species
      implicit none
      integer :: in_unit,isp,ipb,ixp,iTg,iTb,iCj,ierr
      integer :: nsp,nTg,nTb,nCj,npb,nxp
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      call dprint('read_aldbf: Surrogate names')
      nsp = number_of_species
      nTg = albdf_nTg
      nTb = albdf_nTb
      nCj = albdf_nCj
      npb = albdf_nbands
      
      !-----------------------------------------------------------------
      !Read the data
      !-----------------------------------------------------------------
      call dprint('read_aldbf: Read the data')
      do isp=1,nsp
         !Surrogate name for the number of species mole fractions
         !in the ALBDF of species isp
         nxp = albdf_nx(isp)
      
         !Skip species for which no ALBDF is available
         if (nxp.le.0) cycle; if (albdf_file(isp).eq.'null') cycle
         
         !Prepare input
         call CheckFileExists(albdf_file(isp))                          !Check if data file exists
         in_unit = get_file_unit()                                      !Get unit
         
         if (albdf_unformatted) then
            open(unit=in_unit,file=trim(albdf_file(isp)),&              !Open the unit to be read
               form='unformatted',action='read')
            read(in_unit) albdf_darr(1:nCj,1:nTg,1:nxp,isp,1:nTb,1:npb) !Read the ALBDF all at once
            
         else
            open(unit=in_unit,file=trim(albdf_file(isp)),&              !Open the unit to be read
               form='formatted',action='read')
            do ipb=1,npb
               do ixp=1,nxp
                  do iTg=1,nTg
                     do iTb=1,nTb
                        do iCj=1,nCj
                           read(in_unit,*,iostat=ierr) &                !Read de ALBDF line by line
                              albdf_darr(iCj,iTg,ixp,isp,iTb,ipb)
                           if (ierr.lt.0) &                             !Stop if the file ends suddenly
                              call shutdown('read_aldbf: &              
                                             &file missing data')
                        enddo
                     enddo
                  enddo
               enddo   
            enddo
         endif
         
         !Close the unit
         close(in_unit)
         
      enddo

   endsubroutine read_aldbf
   
   !====================================================================
   !Subroutine to fully load the ALBDF
   !====================================================================
   subroutine load_albdf
   
      if (albdf_loaded) return
      call read_albdf_info_1
      call allocate_albdf_parameters
      call read_albdf_info_2
      call read_aldbf
      albdf_loaded = .true.
   
   endsubroutine load_albdf

   !====================================================================
   !Function to interpolate the ALBDF database
   !====================================================================
   real(dp) function albdf_ss(Tloc,xloc,Tsrc,CFin,id_spec,&
                              bslw_band_index,invert)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint,shutdown
      use math_functions, only: locate
      use precision_parameters, only: dp,small
      implicit none
      integer,intent(in) :: id_spec
      integer,optional :: bslw_band_index
      integer :: ipb
      integer :: iTg,iTb,iCF,ixs
      integer :: nTg,nTb,nCF,nxs
      integer :: lTg,lTb,lCF,lxs
      integer :: uTg,uTb,uCF,uxs
      logical,intent(in),optional :: invert
      logical :: compute_F,compute_C,get_darr,get_iarr
      real(dp),intent(in) :: Tloc,Tsrc,xloc,CFin
      real(dp) :: Qtg,Qtb,Qcf,Qxs,Rtg,Rtb,Rcf,Rxs,xl,xu,xval
      
      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('interpolate_albdf: Preparatory procedures')
      
      !Surrogate names
      nCF = albdf_nCj
      nTb = albdf_nTb
      nTg = albdf_nTg
      nxs = albdf_nx(id_spec)
      
      !Set default bslw_band_index value
      call dprint('interpolate_albdf: &
                  &Set default bslw_band_index value')
      ipb = 1; if (present(bslw_band_index)) ipb = bslw_band_index

      !Set interpolation flags
      compute_C = .false.; if (present(invert)) compute_C = invert
      compute_F = .not.compute_C
      get_iarr = .false.; get_darr = .false.
      if (compute_C.and.(.not.albdf_inverted)) get_darr = .true.
      if (compute_F.and.(.not.albdf_inverted)) get_iarr = .true.
      if (compute_C.and.albdf_inverted)        get_iarr = .true.
      if (compute_F.and.albdf_inverted)        get_darr = .true.
      
      !Set x value for the interpolation
      xval = CFin

      !-----------------------------------------------------------------
      !Find upper and lower indexes
      !-----------------------------------------------------------------
      call dprint('interpolate_albdf: Find upper and lower indexes')
      !Gas temperature
      lTg = locate(albdf_Tg,Tloc,nTg)
      lTg = max(1,min(lTg,nTg-1)); uTg = lTg + 1
      if (Tloc.lt.albdf_Tg(1))   uTg = lTg
      if (Tloc.gt.albdf_Tg(nTg)) lTg = uTg

      !Source temperature
      lTb = locate(albdf_Tb,Tsrc,nTb)
      lTb = max(1,min(nTb-1,lTb)); uTb = lTb + 1
      if (Tsrc.lt.albdf_Tb(1))   uTb = lTb
      if (Tsrc.gt.albdf_Tb(nTb)) lTb = uTb

      !For mole fraction
      lxs = locate(albdf_x(:,id_spec),xloc,nxs)                         !Locate lower index
      lxs = max(1,min(nxs-1,lxs)); uxs = min(lxs+1,nxs)                 !Correct lower index, compute upper index
      if (xloc.lt.albdf_x(1,id_spec))    uxs = lxs                      !This is to prevent extrapolations
      if (xloc.gt.albdf_x(nxs,id_spec))  lxs = uxs                      !  (instead, simply take the value at the extreme)

      !Absorption cross-section
      if (get_iarr) then
         lCF = locate(albdf_iarr,xval,nCF)
         lCF = max(1,min(nCF-1,lCF)); uCF = min(lCF+1,nCF)
         if (xval.lt.albdf_iarr(1))   uCF = lCF
         if (xval.gt.albdf_iarr(nCF)) lCF = uCF
         xl = albdf_iarr(lCF); xu = albdf_iarr(uCF)
      endif
      
      !-----------------------------------------------------------------
      !Interpolation
      !-----------------------------------------------------------------
      !Initial values for the interpolation on mole fraction
      Qxs = 0._dp
      Rxs = (xloc - albdf_x(lxs,id_spec))/&
         (albdf_x(uxs,id_spec) - albdf_x(lxs,id_spec) + small)
      if (lxs.eq.uxs) Rxs = 0._dp                                       !If only one mole fraction value is provided,
                                                                        !  do not interpolate in mole fraction
      xloc_loop: do ixs=lxs,uxs
         !Initial values for the interpolation on local temperature
         Qtg = 0._dp
         Rtg = (Tloc - albdf_Tg(ltg))/&
            (albdf_Tg(utg) - albdf_Tg(ltg) + small)
         if (ltg.eq.utg) Rtg = 0._dp                                    !If only one temperature value is provided, 
                     
         Tloc_loop: do itg=ltg,utg
            !Initial values for the interpolation on source temperature
            Qtb = 0._dp
            Rtb = (Tsrc - albdf_Tb(ltb))/&
               (albdf_Tb(utb) - albdf_Tb(ltb) + small)
            if (ltb.eq.utb) Rtb = 0._dp
      
            Tsrc_loop: do itb=ltb,utb
               if (get_darr) then
                  lCF = &
                     locate(albdf_darr(:,iTg,ixs,id_spec,iTb,ipb),xval,nCF)
                  lCF = max(1,min(nCF-1,lCF)); uCF = min(lCF+1,nCF)
                  if (xval.lt.albdf_darr(1,iTg,ixs,id_spec,iTb,ipb)) &
                     uCF = lCF 
                  if (xval.gt.albdf_darr(nCF,iTg,ixs,id_spec,iTb,ipb)) &
                     lCF = uCF
                  xl = albdf_darr(lCF,iTg,ixs,id_spec,iTb,ipb)
                  xu = albdf_darr(uCF,iTg,ixs,id_spec,iTb,ipb)
               endif
               
               !Initial values for the interpolation on C/F
               Qcf = 0._dp; Rcf = (xval - xl)/(xu - xl + small)
               if (lCF.eq.uCF) Rcf = 0._dp
               
               CF_loop: do iCF=lCF,uCF
                  !Interpolate on C or F
                  Rcf = 1._dp - Rcf
                  if (get_iarr) &
                     Qcf = Qcf + Rcf*albdf_darr(iCF,iTg,ixs,id_spec,iTb,ipb)
                  if (get_darr) &
                     Qcf = Qcf + Rcf*albdf_iarr(iCF)
               enddo CF_loop
      
               !Interpolate on source temperature
               Rtb = 1._dp - Rtb
               Qtb = Qtb + Rtb*Qcf
            enddo Tsrc_loop
            
            !Interpolate on local temperature
            Rtg = 1._dp - Rtg
            Qtg = Qtg + Rtg*Qtb
         enddo Tloc_loop

         !Interpolate on mole fraction
         Rxs = 1._dp - Rxs
         Qxs = Qxs + Rxs*Qtg
      enddo xloc_loop

      if (compute_F) Qxs = max(min(1._dp,Qxs),0._dp)
      albdf_ss = Qxs

   endfunction albdf_ss   

   !====================================================================
   !Function to get the ALBDF
   !====================================================================
   real(dp) recursive function albdf_mix(Tloc,xloc,Tsrc,CFin,&
      bslw_band_index,invert,imesh,jmesh,kmesh) result (CF)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: dprint,shutdown
      use global_parameters, only: id_soot,number_of_species
      use math_functions, only: locate
      use precision_parameters, only: small
      implicit none
      integer,optional :: bslw_band_index
      integer :: id_spc,it_counter,max_iter,nsp,spc_counter
      integer :: ipb,isp,single_spec_id
      integer :: lcf,ncf,ucf,xind,yind
      integer :: iimesh,jjmesh,kkmesh
      integer,intent(in),optional :: imesh,jmesh,kmesh
      logical,intent(in),optional :: invert
      logical :: albdf_in_mesh,compute_C,compute_F
      real(dp),intent(in) :: CFin,Tloc,Tsrc,xloc(:)
      real(dp) :: denum,Fj,Rcf
      real(dp) :: spc_cutoff,xmax
      real(dp) :: Cleft,Cright,Ctol,Ctry,xval
      real(dp) :: Fdiff,Fleft,Fright,Ftarget,Ftol,Ftry
      CF = 0._dp                                                        !This is here to avoid compilation warnings

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      call dprint('albdf_mix: Preparatory procedures')

      !Set default bslw_band_index value
      ipb = 1; if (present(bslw_band_index)) ipb = bslw_band_index

      !Set interpolation flags
      compute_C = .false.; if (present(invert)) compute_C = invert
      compute_F = .not.compute_C

      !Set flag for local ALBDF
      iimesh = -1; if (present(imesh)) iimesh = imesh
      jjmesh = -1; if (present(jmesh)) jjmesh = jmesh
      kkmesh = -1; if (present(kmesh)) kkmesh = kmesh
      albdf_in_mesh = .true.
      if (iimesh.le.0) albdf_in_mesh = .false.
      if (jjmesh.le.0) albdf_in_mesh = .false.
      if (kkmesh.le.0) albdf_in_mesh = .false.

      !Set x value for the interpolation
      xval = CFin

      !Surrogate names
      nsp = number_of_species

      !-----------------------------------------------------------------
      !Special case for premixed ALBDFs stored for each grid cell
      !-----------------------------------------------------------------
      if (albdf_in_mesh) then
         !Set up the adequate x and y indexes for the interpolation
         if (compute_F) xind = 1
         if (compute_C) xind = 2
         yind = mod(2,xind)
         
         !Interpolate
         ncf = size(ijk_albdf(:,xind,imesh,jmesh,kmesh))                !Number of C/F discrete values
         lcf = locate(ijk_albdf(:,xind,imesh,jmesh,kmesh),CFin,ncf)     !Locate lower index
         lcf = max(1,min(ncf-1,lcf)); ucf = min(lcf+1,ncf)              !Correct lower index, compute upper index
         Rcf = (CFin - ijk_albdf(lcf,xind,imesh,jmesh,kmesh))/&
                  (ijk_albdf(ucf,xind,imesh,jmesh,kmesh) - &
                     ijk_albdf(lcf,xind,imesh,jmesh,kmesh) + small)
         CF = (1._dp - Rcf)*ijk_albdf(lcf,yind,imesh,jmesh,kmesh) + &
               Rcf*ijk_albdf(ucf,yind,imesh,jmesh,kmesh)
      
         return
      endif

      !-----------------------------------------------------------------
      !Special case for a medium with no or only 
      !one participating species
      !-----------------------------------------------------------------
      call dprint('albdf_mix: Special case for single or no species')

      !Determine how many participating species with 
      !non-negligible concentrations exist
      spc_cutoff = small                                                !Cutoff: any mole fraction below this value 
      single_spec_id = -1; spc_counter = 0                              !  is considered negligible
      do isp=1,nsp
         if (xloc(isp).gt.spc_cutoff) then
            single_spec_id = isp; spc_counter = spc_counter + 1
         endif
      enddo

      !No participating species
      if (spc_counter.eq.0) then
         if (compute_C) CF = 0._dp
         if (compute_F) CF = 1._dp
         return
      endif

      !Only one participating species
      if (spc_counter.eq.1) then
         CF = albdf_ss(Tloc,xloc(single_spec_id),Tsrc,xval,&
            single_spec_id,ipb,compute_C)
         return
      endif

      !-----------------------------------------------------------------
      !Compute F from C
      !-----------------------------------------------------------------
      if (compute_F) then
         call dprint('albdf_mix: Compute F from C')
         select case(trim(slw_mixture_method))  
            case('one_species')   
            !Find the dominating species
            xmax = xloc(1); id_spc = 1
            do isp=2,nsp
               if ((is_slw_species(isp)).and.(xloc(isp).gt.xmax)) then
                  xmax = xloc(isp); id_spc = isp
               endif
            enddo
            Fj = albdf_ss(Tloc,xmax,Tsrc,xval,id_spc,ipb)               !Interpolate the ABLDF for this species
          
            case('multiplication')
            !Only consider species with non-negligible mole fractions
            Fj = 1._dp
            do isp=1,nsp
               if (is_slw_species(isp).and.(xloc(isp).gt.spc_cutoff)) &
                  then
                  denum = xloc(isp)
                  if (isp.eq.id_soot) denum = denum!*1.e5_dp             !1e5 is a scaling factor used in the soot ALBDF 
                  Fj = Fj*albdf_ss(Tloc,xloc(isp),Tsrc,xval/denum,&     !  construction to produce data from suitable
                                   isp,ipb)                             !  values of Cmin and Cmax
               endif
            enddo
        
            case('superposition')
            !Only consider species with non-negligible mole fractions
            spc_counter = 0; Fj = 0._dp
            do isp=1,nsp
               if (is_slw_species(isp).and.(xloc(isp).gt.spc_cutoff)) &
                  then
                  denum = xloc(isp) 
                  if (isp.eq.id_soot) denum = 1._dp
                  Fj = Fj*albdf_ss(Tloc,xloc(isp),Tsrc,xval/denum,&
                                   isp,ipb)
                  spc_counter = spc_counter + 1
               endif
            enddo
            Fj = real(1-spc_counter,dp) + Fj
      
            case('multiple_integration')
            Fj = 1._dp                                                  !Set default value, in case there
            do isp=1,nsp                                                !  are no participating species
               denum = 1._dp
               if (isp.eq.id_soot) denum = xloc(isp)*1.e5_dp
               if (xloc(isp).gt.small) &
               Fj = albdf_ss(Tloc,xloc(isp),Tsrc,xval/denum,isp,ipb)
            enddo

            case default
               call shutdown('albdf_mix: incorrect slw_mixture_method')
      
         endselect
         CF = Fj
         return
      endif

      !-----------------------------------------------------------------
      !Compute C from F (bisection method)
      !-----------------------------------------------------------------
      if (compute_C) then
         call dprint('albdf_mix: Compute C from F')

         !Set initial parameters for the iterative process
         Ftarget = xval                                                 !Target value
         Ftol = 1.e-7_dp; Ctol = 1.e-5_dp                               !Tolerances for the iterative process
         it_counter = 0                                                 !Counter for the number of iterations
         max_iter = 100                                                 !Maximum number of iterations
         Cleft = 1.e-18_dp;   Fleft = 0._dp                             !C and F to the left 
         Cright = 10000._dp;  Fright = 1._dp                            !  and to the right
      
         !Check for extreme values of Ftarget
         if (dabs(Ftarget-Fleft).le.small) then
            CF = Cleft
            return
         elseif (dabs(Ftarget-Fright).le.small) then
            CF = Cright
            return
         endif

         !Inversion loop
         inversion_loop: do
            it_counter = it_counter + 1                                 !Update counter
            Ctry = dexp(0.5_dp*(log(Cleft)+log(Cright)))                !Define Ctry from Cleft and Cright
            Ftry = albdf_mix(Tloc,xloc,Tsrc,Ctry,ipb,imesh=imesh,&      !Compute Ftry from Ctry
                             jmesh=jmesh,kmesh=kmesh)                   
            Fdiff = (Ftarget - Ftry)/(Ftarget + small)                  !Error
            if ((dabs(Fdiff).le.Ftol).or.&                              !Check convergence
                (dabs(Ftarget - Ftry).le.Ftol)) exit inversion_loop
            if (dabs(Cleft-Cright).lt.Ctol) exit inversion_loop         !Escape if the C interval is too small
            if (it_counter.gt.max_iter) &                               !Halt if max iterations exceeded
               call shutdown('albdf_mix: Counter exceeded')         
            if ((Ftry - Ftarget)/(Fleft - Ftarget).lt.0) then           !Update boundary values
               Cright = Ctry; Fright = Ftry
            elseif ((Ftry - Ftarget)/(Fright - Ftarget).lt.0) then
               Cleft = Ctry; Fleft = Ftry
            else

               call shutdown('albdf_mix: Problem in ALBDF inversion')
            endif
         enddo inversion_loop

         !Final value
         CF = Ctry
         return
      endif

   endfunction albdf_mix
 
   !====================================================================
   !Function to define the absorption cross-section 
   !of a gas from its supplementar cross-sections
   !====================================================================
   real(dp) function get_slw_single_cj(ccj0,ccj1)
   
      implicit none
      real(dp),intent(in) :: ccj0,ccj1
      get_slw_single_cj = sqrt(ccj0*ccj1)                               !For now, only this approach
                                                                        !has been implemented   
   endfunction get_slw_single_cj
   
   !====================================================================
   !Subroutine to define the main and supplementar
   !absorption cross-sections for all gray gases
   !====================================================================
   subroutine get_slw_cj(cmin,cmax,ngg,ccj,ccj_sup)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer,intent(in) :: ngg
      integer :: ii
      real(dp),dimension(:),intent(out) :: ccj,ccj_sup
      real(dp),intent(in) :: cmin,cmax
   
      !-----------------------------------------------------------------
      !Computing the suplementar cross-sections
      !-----------------------------------------------------------------
      do ii=1,ngg
         ccj_sup(ii) = cmin*(cmax/cmin)**&
                                 (real(ii,dp)/real(ngg,dp))
      enddo
      
      !-----------------------------------------------------------------
      !Computing the gray gas absorption cross-sections
      !-----------------------------------------------------------------
      ccj(1) = sqrt(ccj_sup(1)*slw_cmin)
      do ii=2,ngg
         ccj(ii) = get_slw_single_cj(ccj_sup(ii-1),ccj_sup(ii))
      enddo
      
   endsubroutine get_slw_cj
   
   !====================================================================
   !Function to define the ALDBFs divisions (for RC-SLW)
   !====================================================================
   subroutine get_slw_Fj(ngg,ffj,ffj_sup)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use math_functions, only: quad_gauss_legendre,&
                                quad_gauss_chebyshev
      implicit none
      integer,intent(in) :: ngg
      integer :: ii,kk
      real(dp) :: Fmax,Fmin,sum_w
      real(dp),dimension(:) :: ffj,ffj_sup
      real(dp) :: x_i(1:ngg),w_i(1:ngg)
   
      !-----------------------------------------------------------------
      !Setting bounds for the supplementar ALBDF
      !-----------------------------------------------------------------
      Fmin = slw_Fmin
      Fmax = slw_Fmax
   
      !-----------------------------------------------------------------
      !Get the quadrature with standard bounds (-1 to 1)
      !-----------------------------------------------------------------
      if (trim(slw_quadrature).eq.'Gauss-Legendre') then
         call quad_gauss_legendre(-1._dp,1._dp,ngg,x_i,w_i)
      elseif (trim(slw_quadrature).eq.'Gauss-Chebyshev') then                !By default, use the 
         call quad_gauss_chebyshev(ngg,'even-rank',x_i,w_i)             !even GC quadrature
      elseif (trim(slw_quadrature).eq.'Gauss-Chebyshev-even') then
         call quad_gauss_chebyshev(ngg,'even-rank',x_i,w_i)
      elseif (trim(slw_quadrature).eq.'Gauss-Chebyshev-odd') then
         call quad_gauss_chebyshev(ngg,'odd-rank',x_i,w_i)
      else
         call shutdown('Problem in the specification of slw_quadrature')
      endif
      
      !-----------------------------------------------------------------
      !Adjusting the quadrature limits
      !-----------------------------------------------------------------
      do ii=1,ngg
         x_i(ii) = 0.5_dp*x_i(ii) + 0.5_dp
         w_i(ii) = 0.5_dp*w_i(ii)
      enddo
      
      !-----------------------------------------------------------------
      !Computing the ALBDFs
      !-----------------------------------------------------------------
      do ii=1,ngg
         sum_w = 0._dp
         do kk=1,ii
            sum_w = sum_w + w_i(kk)
         enddo
         ffj_sup(ii) = Fmin + (Fmax - Fmin)*sum_w                        !Supplementar
         ffj(ii) = Fmin + x_i(ii)*(Fmax - Fmin)
      enddo

   endsubroutine get_slw_Fj

   !====================================================================
   !Function to compute the gray gas 
   !absorption coefficient in the SLW model
   !====================================================================
   real(dp) function slw_kappa_func(ttmp,press,mole_frac,ccj,&
                                    compute_cj,mixing_method)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use constants, only: atm_pa,Ru
      implicit none
      character(*),optional,intent(in) :: mixing_method
      character(200) :: mixing_method_aux
      real(dp),intent(in) :: ttmp,press,ccj,mole_frac
      real(dp) :: mol_density,molm
      logical,optional :: compute_cj
      logical :: reverse

      !-----------------------------------------------------------------
      !Set up parameters
      !-----------------------------------------------------------------
      !Set up optional parameter for the reverse calculation
      reverse = .false.; if (present(compute_cj)) reverse = compute_cj  
  
      !Set up optional parameter for the mixing method
      mixing_method_aux = trim(slw_mixture_method)
      if (present(mixing_method)) mixing_method_aux = trim(mixing_method)
      mol_density = atm_pa*press/(ttmp*Ru)                              !Mole density

      !-----------------------------------------------------------------
      !Particularize the calculation depending on the mixture method
      !-----------------------------------------------------------------
      select case(trim(mixing_method_aux))   
         case('one_species')
            molm = mol_density*mole_frac
            
         case('multiplication')
            molm = mol_density
   
         case('superposition')
            molm = mol_density

         case('multiple_integration')
            molm = mol_density*mole_frac
         
         case('albdf_precombined')
            molm = mol_density
            
         case('SLW1')
            molm = mol_density*mole_frac
            
         case default
            call shutdown('Problem in the specification of&
                          & slw_mixture_method')
      endselect

      !-----------------------------------------------------------------
      !Finish up the calculation
      !-----------------------------------------------------------------
      if (reverse) then                                                 !Compute cj from kappaj = ccj
         slw_kappa_func = ccj/molm
      else                                                              !Compute kappaj from cj = ccj (default)
         slw_kappa_func = molm*ccj
      endif
         
   endfunction slw_kappa_func

   !====================================================================
   !Function to compute the gray gas 
   !weighting coefficient in the SLW model
   !====================================================================
   real(dp) function slw_a_func(ttmp,xxs,ttsource,ccj_sup1,ccj_sup0,&
                                bslw_band_index)
      
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      implicit none
      integer,optional :: bslw_band_index
      integer :: ipb
      real(dp),intent(in) :: ttmp,ttsource,xxs(:),ccj_sup1
      real(dp),optional :: ccj_sup0
      real(dp) :: F1,F0
      logical :: slw_twindow
      
      !-----------------------------------------------------------------
      !Set up flags
      !-----------------------------------------------------------------
      ipb = 1; if (present(bslw_band_index)) ipb = bslw_band_index      !Set up optional parameter
      slw_twindow = .false.; if(.not.present(ccj_sup0)) &               !Check if the gas is a transparent window
                                                   slw_twindow = .true.

      !-----------------------------------------------------------------
      !Compute ALBDFs
      !-----------------------------------------------------------------
      F1 = albdf_mix(ttmp,xxs,ttsource,ccj_sup1,ipb)
      if (slw_twindow) then
         F0 = 0._dp
      else
         F0 = albdf_mix(ttmp,xxs,ttsource,ccj_sup0,ipb)
      endif

      !-----------------------------------------------------------------
      !Apply the definition of the weighting coefficient
      !-----------------------------------------------------------------
      slw_a_func =  F1 - F0
            
   endfunction slw_a_func
   
   !====================================================================
   !Function to check if species isp is a participating 
   !species  in the framework of the SLW model
   !====================================================================
   logical function is_slw_species(isp)

      integer,intent(in) :: isp
      is_slw_species = (albdf_nx(isp).gt.0).and.&
                       (albdf_file(isp).ne.'null')

   endfunction is_slw_species

   !====================================================================
   !Function to compute the Planck-mean
   !absorption coefficient with the SLW model
   !====================================================================
   real(dp) function get_slw_kp(Tgas,Pgas,Xgas,Tsource,ngas,&
                                prepare_albdf,bslw_band_index)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint
      implicit none
      integer :: ierr,jgas,ipb
      integer,intent(in) :: ngas
      integer,optional :: bslw_band_index
      logical,optional :: prepare_albdf
      logical :: go_albdf,transparent_window
      real(dp),intent(in) :: Pgas,Tgas,Tsource,Xgas(:)
      real(dp) :: cj_sup_loc,F0
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref
      real(dp),allocatable,dimension(:) :: Fj_ref,Fj_sup_ref
      real(dp) :: a_j,kappa_j,sum_kp
a_j = Pgas  !Added just to avoid a compilation warning
      !-----------------------------------------------------------------
      !Set up the optional parameters
      !-----------------------------------------------------------------
      call dprint('get_slw_kp: Set up the optional parameters')
      
      !Flag for loading the ALBDF
      go_albdf = .false.
      if (present(prepare_albdf)) go_albdf = prepare_albdf
      
      !Band index
      ipb = 1; if (present(bslw_band_index)) ipb = bslw_band_index

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('get_slw_kp: Allocating arrays')
      allocate(cj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_ref',ierr)
      allocate(cj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_sup_ref',ierr)
      allocate(Fj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_ref',ierr)
      allocate(Fj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_sup_ref',ierr)

      !-----------------------------------------------------------------
      !Load the ALBDF, if requested
      !-----------------------------------------------------------------
      if (go_albdf) then
         call dprint('get_slw_kp: Load the ALBDF')
         call load_albdf
      endif
      
      !-----------------------------------------------------------------
      !Define cross-sections
      !-----------------------------------------------------------------
      call dprint('get_slw_kp: Define cross-sections')
      if (trim(slw_nonuniform_method).eq.'rank_correlated') then
         call get_slw_Fj(ngas,Fj_ref(1:ngas),Fj_sup_ref(1:ngas))
         Fj_sup_ref(0) = slw_Fmin  
      else
         call get_slw_cj(slw_cmin,slw_cmax,ngas,cj_ref(1:ngas),&
                         cj_sup_ref(1:ngas))
         cj_sup_ref(0) = slw_cmin         
      endif

      !-----------------------------------------------------------------
      !Main loop
      !-----------------------------------------------------------------
      call dprint('get_slw_kp: Main loop')
      sum_kp = 0._dp
      gas_loop: do jgas=0,ngas
         !Check if the gas is a transparent window
         transparent_window = .false.                                   !The transparent window should be included in 
         if (jgas.eq.0) transparent_window = .true.                     !  the gas_loop loop to initialize F0
      
         !Compute properties for the gas
         call compute_slw_gas_parameters(kappa_j,a_j,cj_sup_ref(jgas),&
            cj_sup_loc,Fj_ref(jgas),Fj_sup_ref(jgas),F0,Tgas,Xgas,&
            reference_T=Tsource,reference_xs=Xgas,transparent_window=&
            transparent_window)

         !Summing for the Planck-mean absorption coefficient
         sum_kp = sum_kp + kappa_j*a_j
            
      enddo gas_loop
      get_slw_kp = sum_kp
      
      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('get_slw_kp: Deallocate arrays')
      deallocate(cj_ref,cj_sup_ref)
      deallocate(Fj_ref,Fj_sup_ref)

   endfunction get_slw_kp 
   
   !====================================================================
   !Function to compute the gas emissivity with the SLW model
   !====================================================================
   real(dp) function get_slw_emissivity(Tgas,Pgas,Xgas,Tsource,ngas,&
                                   length,prepare_albdf,bslw_band_index)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint
      implicit none
      integer :: ierr,jgas,ipb
      integer,intent(in) :: ngas
      integer,optional :: bslw_band_index
      logical,optional :: prepare_albdf
      logical :: go_albdf,transparent_window
      real(dp),intent(in) :: length,Pgas,Tgas,Tsource,Xgas(:)
      real(dp) :: a_j,kappa_j,sum_emi
      real(dp) :: cj_sup_loc,F0
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref
      real(dp),allocatable,dimension(:) :: Fj_ref,Fj_sup_ref      
a_j = Pgas  !Added just to avoid a compilation warning
      !-----------------------------------------------------------------
      !Set up the optional parameters
      !-----------------------------------------------------------------
      call dprint('get_slw_emissivity: Set up the optional parameters')
      
      !Flag for loading the ALBDF
      go_albdf = .false.
      if (present(prepare_albdf)) go_albdf = prepare_albdf
      
      !Band index
      ipb = 1; if (present(bslw_band_index)) ipb = bslw_band_index

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('get_slw_emissivity: Allocating arrays')
      allocate(cj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_ref',ierr)
      allocate(cj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_sup_ref',ierr)
      allocate(Fj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_ref',ierr)
      allocate(Fj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_sup_ref',ierr)

      !-----------------------------------------------------------------
      !Load the ALBDF, if requested
      !-----------------------------------------------------------------
      if (go_albdf) then
         call dprint('get_slw_emissivity: Load the ALBDF')
         call load_albdf
      endif
      
      !-----------------------------------------------------------------
      !Define cross-sections
      !-----------------------------------------------------------------
      call dprint('get_slw_emissivity: Define cross-sections')
      if (trim(slw_nonuniform_method).eq.'rank_correlated') then
         call get_slw_Fj(ngas,Fj_ref(1:ngas),Fj_sup_ref(1:ngas))
         Fj_sup_ref(0) = slw_Fmin  
      else
         call get_slw_cj(slw_cmin,slw_cmax,ngas,cj_ref(1:ngas),&
                         cj_sup_ref(1:ngas))
         cj_sup_ref(0) = slw_cmin         
      endif

      !-----------------------------------------------------------------
      !Main loop
      !-----------------------------------------------------------------
      call dprint('get_slw_emissivity: Main loop')
      sum_emi = 0._dp
      gas_loop: do jgas=0,ngas
         !Check if the gas is a transparent window
         transparent_window = .false.                                   !The transparent window should be included in 
         if (jgas.eq.0) transparent_window = .true.                     !  the gas_loop loop to initialize F0
      
         !Compute properties for the gas
         call compute_slw_gas_parameters(kappa_j,a_j,cj_sup_ref(jgas),&
            cj_sup_loc,Fj_ref(jgas),Fj_sup_ref(jgas),F0,Tgas,Xgas,&
            reference_T=Tsource,reference_xs=Xgas,transparent_window=&
            transparent_window)

         !Summing the emissivity
         sum_emi = sum_emi + a_j*(1._dp - dexp(-kappa_j*length))
      enddo gas_loop
      get_slw_emissivity = sum_emi

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('get_slw_emissivity: Deallocate arrays')
      deallocate(cj_ref,cj_sup_ref)
      deallocate(Fj_ref,Fj_sup_ref)

   endfunction get_slw_emissivity

   !====================================================================
   !Function to compute the heat flux for a uniform one-dimensional
   !gas column using the SLW model
   !====================================================================
   real(dp) function get_slw_uniform_flux(Tgas,Pgas,Xgas,Tsource,ngas,&
                                   xpos,prepare_albdf,bslw_band_index)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint
      use constants, only: pi,sigrpi
      use math_functions, only: expint
      implicit none
      integer :: ierr,jgas,ipb
      integer,intent(in) :: ngas
      integer,optional :: bslw_band_index
      logical,optional :: prepare_albdf
      logical :: go_albdf,transparent_window
      real(dp),intent(in) :: Pgas,Tgas,Tsource,Xgas(:),xpos
      real(dp) :: a_j,Ib,kappa_j,sum_F
      real(dp) :: cj_sup_loc,F0
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref
      real(dp),allocatable,dimension(:) :: Fj_ref,Fj_sup_ref      
a_j = Pgas  !Added just to avoid a compilation warning
      !-----------------------------------------------------------------
      !Set up the optional parameters
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_flux: Set up the optional parameters')
      
      !Flag for loading the ALBDF
      go_albdf = .false.
      if (present(prepare_albdf)) go_albdf = prepare_albdf
      
      !Band index
      ipb = 1; if (present(bslw_band_index)) ipb = bslw_band_index

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_flux: Allocating arrays')
      allocate(cj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_ref',ierr)
      allocate(cj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_sup_ref',ierr)
      allocate(Fj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_ref',ierr)
      allocate(Fj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_sup_ref',ierr)

      !-----------------------------------------------------------------
      !Load the ALBDF, if requested
      !-----------------------------------------------------------------
      if (go_albdf) then
         call dprint('get_slw_uniform_flux: Load the ALBDF')
         call load_albdf
      endif
      
      !-----------------------------------------------------------------
      !Define cross-sections
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_flux: Define cross-sections')
      call get_slw_cj(slw_cmin,slw_cmax,ngas,cj_ref(1:ngas),&
                      cj_sup_ref(1:ngas))
      cj_sup_ref(0) = slw_cmin

      !-----------------------------------------------------------------
      !Main loop
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_flux: Main loop')
      sum_F = 0._dp
      gas_loop: do jgas=0,ngas
         !Check if the gas is a transparent window
         transparent_window = .false.                                   !The transparent window should be included in 
         if (jgas.eq.0) transparent_window = .true.                     !  the gas_loop loop to initialize F0
      
         !Compute properties for the gas
         call compute_slw_gas_parameters(kappa_j,a_j,cj_sup_ref(jgas),&
            cj_sup_loc,Fj_ref(jgas),Fj_sup_ref(jgas),F0,Tgas,Xgas,&
            reference_T=Tsource,reference_xs=Xgas,transparent_window=&
            transparent_window,nonuniform_method='uniform')

         !Summing the emissivity
         Ib = sigrpi*(Tsource**4)
         sum_F = sum_F + &
                 pi*a_j*Ib*(1._dp - 2._dp*expint(3,kappa_j*xpos))
         
      enddo gas_loop
      get_slw_uniform_flux = sum_F

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_flux: Deallocate arrays')
      deallocate(cj_ref,cj_sup_ref)
      deallocate(Fj_ref,Fj_sup_ref)

   endfunction get_slw_uniform_flux

   !====================================================================
   !Function to compute the heat flux for a uniform one-dimensional
   !gas column using the SLW model
   !====================================================================
   real(dp) function get_slw_uniform_source(Tgas,Pgas,Xgas,Tsource,&
                         ngas,xpos,length,prepare_albdf,bslw_band_index)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: CheckMemAlloc,dprint
      use constants, only: sigrpi,twopi
      use math_functions, only: expint
      implicit none
      integer :: ierr,jgas,ipb
      integer,intent(in) :: ngas
      integer,optional :: bslw_band_index
      logical,optional :: prepare_albdf
      logical :: go_albdf,transparent_window
      real(dp),intent(in) :: length,Pgas,Tgas,Tsource,Xgas(:),xpos
      real(dp) :: a_j,Ib,kappa_j,sum_Q
      real(dp) :: cj_sup_loc,F0
      real(dp),allocatable,dimension(:) :: cj_ref,cj_sup_ref
      real(dp),allocatable,dimension(:) :: Fj_ref,Fj_sup_ref      
a_j = Pgas  !Added just to avoid a compilation warning
      !-----------------------------------------------------------------
      !Set up the optional parameters
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_source: &
                  &Set up the optional parameters')
      
      !Flag for loading the ALBDF
      go_albdf = .false.
      if (present(prepare_albdf)) go_albdf = prepare_albdf
      
      !Band index
      ipb = 1; if (present(bslw_band_index)) ipb = bslw_band_index

      !-----------------------------------------------------------------
      !Allocating arrays
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_source: Allocating arrays')
      allocate(cj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_ref',ierr)
      allocate(cj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('cj_sup_ref',ierr)
      allocate(Fj_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_ref',ierr)
      allocate(Fj_sup_ref(0:ngas),stat=ierr)
      call CheckMemAlloc('Fj_sup_ref',ierr)

      !-----------------------------------------------------------------
      !Load the ALBDF, if requested
      !-----------------------------------------------------------------
      if (go_albdf) then
         call dprint('get_slw_uniform_source: Load the ALBDF')
         call load_albdf
      endif
      
      !-----------------------------------------------------------------
      !Define cross-sections
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_source: Define cross-sections')
      call get_slw_cj(slw_cmin,slw_cmax,ngas,cj_ref(1:ngas),&
                      cj_sup_ref(1:ngas))
      cj_sup_ref(0) = slw_cmin

      !-----------------------------------------------------------------
      !Main loop
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_source: Main loop')
      sum_Q = 0._dp
      gas_loop: do jgas=0,ngas
         !Check if the gas is a transparent window
         transparent_window = .false.                                   !The transparent window should be included in 
         if (jgas.eq.0) transparent_window = .true.                     !  the gas_loop loop to initialize F0
      
         !Compute properties for the gas
         call compute_slw_gas_parameters(kappa_j,a_j,cj_sup_ref(jgas),&
            cj_sup_loc,Fj_ref(jgas),Fj_sup_ref(jgas),F0,Tgas,Xgas,&
            reference_T=Tsource,reference_xs=Xgas,transparent_window=&
            transparent_window,nonuniform_method='uniform')

         !Summing the emissivity
         Ib = sigrpi*(Tsource**4)
         sum_Q = sum_Q - twopi*a_j*kappa_j*Ib*&
            (expint(2,kappa_j*xpos) + (expint(2,kappa_j*(length-xpos))))
         
      enddo gas_loop
      get_slw_uniform_source = sum_Q

      !-----------------------------------------------------------------
      !Deallocate arrays
      !-----------------------------------------------------------------
      call dprint('get_slw_uniform_source: Deallocate arrays')
      deallocate(cj_ref,cj_sup_ref)
      deallocate(Fj_ref,Fj_sup_ref)

   endfunction get_slw_uniform_source

   !====================================================================
   !Routine to determine the absorption coefficient 
   !and the weighting coefficient of a gray gas
   !====================================================================
   !Inputs:
   !  Tloc -> local temperature
   !  xsloc -> array with the local species mole fraction
   !  csup_ref -> reference supplemental cross-section, 
   !     \tilde{C}_{j}^{ref}
   !  csup_loc0 -> local supplemental cross-section @ j-1, 
   !     \tilde{C}_{j-1}^{ref}
   !  Fref -> ALBDF @ reference cross-section, F_{ref}
   !  Fsup_ref -> ALBDF @ reference supplemental cross-section,
   !     \tilde{F}_{j}^{ref}
   !  F0 -> lower bound of blackbody energy fraction to be used for the
   !     calculation of the weighting coefficient; this value is updated
   !     to the one corresponding to the next supplemental cross-section
   !  reference_T, reference_xs -> reference temperature and species
   !     mole fractions. Optional variables; if not provided, the values
   !     of slw_Tref and slw_xsref are used
   !  nonuniform_method -> string specifying which nonuniform treatment
   !     to be used. Optional variable; if not provided, the value of
   !     slw_nonuniform_method is used
   !  transparent_window -> if .true., compute the gas parameters for
   !     a transparent window. Optional variable; default = .false.
   !  scaling_coeff -> local scaling coefficient, necessary for the
   !     scaled SLW method. Optional parameter
   !Outputs:
   !  kj -> absorption coefficient
   !  aj -> weighting coefficient
   !  F0, csup_loc0* -> updated for the calculation with the next gray
   !     gas (*depending on the method)
   subroutine compute_slw_gas_parameters(kj,aj,csup_ref,csup_loc0,&
      Fref,Fsup_ref,F0,Tloc,xsloc,reference_T,reference_xs,&
      nonuniform_method,transparent_window,scaling_coeff,&
      imesh,jmesh,kmesh)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use comp_functions, only: shutdown
      use global_parameters, only: number_of_species
      implicit none
      character(*),intent(in),optional :: nonuniform_method
      character(200) :: nuni_method
      integer,intent(in),optional :: imesh,jmesh,kmesh
      integer :: isp
      integer :: iimesh,jjmesh,kkmesh
      logical,intent(in),optional :: transparent_window
      logical :: twindow
      real(dp),intent(in) :: Fref,Tloc,xsloc(:)
      real(dp),intent(in),optional :: reference_T,reference_xs(:),&
         scaling_coeff
      real(dp),intent(out) :: aj,kj
      real(dp),intent(inout) :: csup_loc0,csup_ref,Fsup_ref,F0
      real(dp) :: p,Tref,uloc,xmix,xsref(size(xsloc))
      real(dp) :: cloc,cref,csup_loc
      real(dp) :: F1,Floc,Fsup_loc

      !-----------------------------------------------------------------
      !Preparatory procedures
      !-----------------------------------------------------------------
      !For now, only cases with p = 1 atm are supported
      p = 1._dp                                                         

      !Compute total mole fraction of participating species
      xmix = 0._dp
      do isp=1,number_of_species
         if (is_slw_species(isp)) xmix = xmix + xsloc(isp)
      enddo

      !-----------------------------------------------------------------
      !Set up optional parameters
      !-----------------------------------------------------------------
      Tref = slw_Tref
      if (present(reference_T)) Tref = reference_T

      xsref = slw_xsref(1:size(xsloc))
      if (present(reference_xs)) xsref = reference_xs

      nuni_method = slw_nonuniform_method
      if (present(nonuniform_method)) nuni_method = nonuniform_method

      twindow = .false.
      if (present(transparent_window)) twindow = transparent_window

      uloc = 1._dp
      if (present(scaling_coeff)) uloc = scaling_coeff
      
      iimesh = -1; if (present(imesh)) iimesh = imesh
      jjmesh = -1; if (present(jmesh)) jjmesh = jmesh
      kkmesh = -1; if (present(kmesh)) kkmesh = kmesh

      !-----------------------------------------------------------------
      !Set up surrogate names for nonuniform approaches
      !-----------------------------------------------------------------
      if (trim(nuni_method).eq.'RA-SLW') nuni_method = 'I.1.1'
      if (trim(nuni_method).eq.'RC-SLW') nuni_method = 'I.2.2'
      if (trim(nuni_method).eq.'LC-SLW') nuni_method = 'II.1.1'
      if (trim(nuni_method).eq.'reference_approach') &
         nuni_method = 'I.1.1'
      if (trim(nuni_method).eq.'rank_correlated') &
         nuni_method = 'I.2.2'
      if (trim(nuni_method).eq.'locally_correlated') &
         nuni_method = 'II.1.1'
      if (trim(nuni_method).eq.'II.1.2') nuni_method = 'II.1.1'
      if (trim(nuni_method).eq.'II.2.2') nuni_method = 'II.2.1'
      if ((trim(nuni_method).eq.'scaled').and.&
         (.not.present(scaling_coeff))) &
            call shutdown('compute_slw_gas_parameters: scaling_coeff &
                          &must be specified for scaled-SLW method')

      !-----------------------------------------------------------------
      !Compute parameters, following the recipes in Solovjov et al. 
      !(2017) and Webb et al. (2018). The nomenclature of the I and II
      !methods follows the former paper
      !-----------------------------------------------------------------
      selectcase(trim(nuni_method))
      case('uniform')
         !Inputs: csup_loc and csup_loc0
         cloc = get_slw_single_cj(csup_loc0,csup_ref)                   !Local cross-section determined from the supplementar ones
         F1 = albdf_mix(Tloc,xsloc,Tloc,csup_ref)
         csup_loc0 = csup_ref         

      case('scaled')
         !Inputs csup_ref and csup_ref0
         cloc = get_slw_single_cj(csup_loc0,csup_ref)                   !Local cross-section determined from the supplementar ones
         F1 = albdf_mix(Tref,xsref,Tloc,csup_ref,&
                        imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         csup_loc0 = csup_ref

      case('I.1.1')
         !Inputs: csup_ref and csup_loc0
         Fsup_ref = albdf_mix(Tref,xsref,Tref,csup_ref,&                !ALBDF evaluated at the reference state, reference Tb
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         csup_loc = albdf_mix(Tloc,xsloc,Tref,Fsup_ref,invert=.true.,&  !Get the local supplementar cross-section
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         cloc = get_slw_single_cj(csup_loc0,csup_loc)                   !Local cross-section determined from the supplementar ones
         F1 = albdf_mix(Tref,xsref,Tloc,csup_ref,&                      !ALBDF evaluated at the reference state, local Tb
                        imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         csup_loc0 = csup_loc                                           !Update local supplemental cross-section for the calculation
                                                                        !  with the next gray gas
      case('I.1.2')
         !Inputs: csup_ref and csup_loc0
         Fsup_ref = albdf_mix(Tref,xsref,Tref,csup_ref,&                !ALBDF evaluated at the reference state, reference Tb
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         csup_loc = albdf_mix(Tloc,xsloc,Tref,Fsup_ref,invert=.true.,&  !Get the local supplementar cross-section
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         cloc = get_slw_single_cj(csup_loc0,csup_loc)                   !Local cross-section determined from the supplementar ones
         F1 = albdf_mix(Tloc,xsloc,Tloc,csup_loc,&                      !ALBDF evaluated at the reference state, local Tb
                        imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         csup_loc0 = csup_loc                                           !Update local supplemental cross-section for the calculation
                                                                        !  with the next gray gas
      case('I.2.1')
         !Inputs: Fref and Fsup_ref
         cloc = albdf_mix(Tloc,xsloc,Tref,Fref,invert=.true.,&          !Local cross-section
                          imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         csup_ref = albdf_mix(Tref,xsref,Tref,Fsup_ref,invert=.true.,&  !Reference supplemental cross-section
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         F1 = albdf_mix(Tref,xsref,Tloc,csup_ref)                       !ALBDF evaluated at the reference state, local Tb

      case('I.2.2')
         !Inputs: Fref and Fsup_ref
         if (.not.twindow) cloc = albdf_mix(Tloc,xsloc,Tref,Fref,&      !Local cross-section
            invert=.true.,imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         csup_loc = albdf_mix(Tloc,xsloc,Tref,Fsup_ref,invert=.true.,&  !Local supplemental cross-section
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         F1 = albdf_mix(Tloc,xsloc,Tloc,csup_loc,&                      !ALBDF evaluated at the reference state, local Tb
                        imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)

      case('II.1.1')
         !Inputs: csup_ref and csup_loc0
         Fsup_loc = albdf_mix(Tref,xsref,Tloc,csup_ref,&                !ALBDF evaluated at the reference state, local Tb
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         csup_loc = albdf_mix(Tloc,xsloc,Tloc,Fsup_loc,invert=.true.,&  !Local supplemental cross-section
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         cloc = get_slw_single_cj(csup_loc0,csup_loc)                   !Local cross-section determined from the supplementar ones
         F1 = Fsup_loc                                                  !ALBDF evaluated at the reference state, local Tb
         csup_loc0 = csup_loc                                           !Update local supplemental cross-section for the calculation
                                                                        !  with the next gray gas
      case ('II.2.1')
         !Inputs: Fref and Fsup_ref
         csup_ref = albdf_mix(Tref,xsref,Tref,Fsup_ref,invert=.true.,&  !Reference supplemental cross-section
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         cref = albdf_mix(Tref,xsref,Tref,Fref,invert=.true.,&          !Reference cross-section
                          imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         Fsup_loc = albdf_mix(Tref,xsref,Tloc,csup_ref,&                !ALBDF @ reference state, local Tb
                              imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         Floc = albdf_mix(Tref,xsref,Tloc,cref,&                        !ALBDF @ reference state, local Tb
                          imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         cloc = albdf_mix(Tloc,xsloc,Tref,Floc,invert=.true.,&          !Local cross-section
                          imesh=iimesh,jmesh=jjmesh,kmesh=kkmesh)
         F1 = Fsup_loc                                                  !ALBDF @ the reference state, local Tb

      case default
         call shutdown('compute_slw_gas_parameters: nonuniform method &
                       &unspecified')
      endselect

      !-----------------------------------------------------------------
      !Finish the calculation
      !-----------------------------------------------------------------
      if (twindow) F0 = 0._dp
      if (twindow) cloc = 0._dp
      kj = uloc*slw_kappa_func(Tloc,p,xmix,cloc)                        !Absorption coefficient
      aj = F1 - F0                                                      !Weighting coefficient
      F0 = F1                                                           !Update F0 value for the calculation with
                                                                        !  the next gray gas
   endsubroutine compute_slw_gas_parameters

!   !====================================================================
!   !Subroutine to compute the reference temperature to the used
!   !in the framework of rank-correlated SLW models
!   !====================================================================
!   subroutine compute_Tref(xx,TT,npts)
   
!      implicit none
!      integer :: npts
!      real(dp),intent(in) :: xx(npts),TT(npts)

!      slw_Tref = maxval(TT)
   
!   endsubroutine compute_Tref


   !====================================================================
   !Subroutine to assemble ALBDFs from NBCKs
   !====================================================================
   subroutine albdf_from_nbck(Carray,Farray,Tgas,pgas,xgas,Tsource)
   
      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: sigrpi
      use nbck_functions, only: nbck_mix_kg
      use nbck_parameters
      use physical_functions, only: bb_emission_frac,Ceta_to_keta
      implicit none
      integer :: icj,inb,ncj,nnb
      real(dp),intent(in) :: Carray(:),pgas,Tgas,Tsource,xgas(:)
      real(dp),intent(out) :: Farray(:)
      real(dp) :: fIb,k
      
      !-----------------------------------------------------------------
      !Surrogate names
      !-----------------------------------------------------------------
      ncj = size(Carray)
      nnb = number_nbck_bands
      
      !-----------------------------------------------------------------
      !Main loop
      !-----------------------------------------------------------------
      Farray = 0._dp
      band_loop: do inb=1,nnb
         fIb = bb_emission_frac(nbck_lbound(inb),Tsource,.true.) - &
                  bb_emission_frac(nbck_ubound(inb),Tsource,.true.)
         c_loop: do icj=1,ncj
            k = Ceta_to_keta(Carray(icj),Tgas)                          !Convert from C to k   
            Farray(icj) = Farray(icj) + &
               fIb*nbck_mix_kg(inb,Tgas,xgas,pgas,k,'g')
         enddo c_loop
      enddo band_loop
   
   endsubroutine albdf_from_nbck


!   !====================================================================
!   !Subroutine to solve for the single gray gas absorption
!   !and weighting coefficientsin the SLW-1 model
!   !====================================================================
!   subroutine solve_slw1_parameters(slw1_spec,val1,val2,len1,len2,&
!                                    a_res,kappa_res) 

!      !-----------------------------------------------------------------
!      !Declaration of variables
!      !-----------------------------------------------------------------
!      use comp_functions, only: shutdown
!      character(*),intent(in) :: slw1_spec
!      integer :: max_iter,counter
!      real(dp),intent(in) :: len1,len2,val1,val2
!      real(dp),intent(out) :: a_res,kappa_res
!      real(dp) :: a_new,diff,iter_tol,k_new,k_old

!      !---------------------
!      !Iteration parameters
!      !---------------------
!      max_iter = 100000
!      iter_tol = 1.e-2_dp
   
!      !----------------------
!      !Initialize parameters
!      !----------------------
!      diff = 1000._dp
!      counter = 0
!      k_old = 1._dp
   
!      !------------------------
!      !Iterate to find k and a
!      !------------------------
!      do while (diff.gt.iter_tol)
!         !Update counter
!         counter = counter + 1
      
!         !Stop the program if the maximum number of iterations is reached
!         if (counter.gt.max_iter) &
!            call shutdown('Maximum number of iterations in &
!                           &solve_slw1_parameters reached')
      
!         selectcase(trim(slw1_spec))
!         case('kappa-epsilon')                                          !Value1: Kp, Value2: epsilon @ len1
!            k_new = (val1/val2)*(1._dp - dexp(-k_old*len1))
!            a_new = val1/k_new
         
!         case('epsilon-epsilon')
!            if (val1.lt.val2) then
!               k_new = -dlog(1._dp - (val1/val2)*&
!                  (1._dp - dexp(-k_old*len2)))/len1
!               a_new = val2/(1._dp - dexp(-k_new*len2))
!            elseif (val2.lt.val1) then
!               k_new = -dlog(1._dp - (val2/val1)*&
!                  (1._dp - dexp(-k_old*len1)))/len2
!               a_new = val1/(1._dp - dexp(-k_new*len1))
!            else
!               call shutdown('Emissivities with the same value in &
!                              &solve_slw1_parameters')
!            endif
      
!         case default
!            call shutdown('SLW-1 method unspecified')
!         endselect
      
!         !Compute difference
!         diff = dabs((k_new-k_old)/(k_new+small))*100._dp

!         !Update k
!         k_old = k_new
!      enddo
   
!      !--------------------
!      !Assign final values
!      !--------------------
!      a_res = a_new
!      kappa_res = k_new
   
!   endsubroutine solve_slw1_parameters



   !====================================================================
   !Subroutine to compute the reference parameters (kappa and a) for
   !the SLW-1 model
   !Input:
   !   Tref -> Reference temperature
   !   xsref -> Reference mole fractions of the species
   !   pref -> Reference pressure
   !Output:
   !   kappa_out -> Reference absorption coefficient for gas 1
   !   a_out -> Reference weighting coefficient for gas 1
   !====================================================================
   subroutine slw1_compute_ref(Tref,xsref,pref,kappa_out,a_out)

      !-----------------------------------------------------------------
      !Declaration of variables
      !-----------------------------------------------------------------
      use constants, only: pi
      use math_functions, only: expint
      use comp_functions, only: shutdown
      use precision_parameters, only: small
      use physical_functions, only: Ib_function
      implicit none
      integer,parameter :: max_iter = 1000
      integer :: counter,i,j,x,y,nx
      real(dp),intent(in) :: pref,Tref,xsref(:)
      real(dp),intent(out) :: a_out,kappa_out
      real(dp),parameter :: iter_tol=1.e-6_dp
      real(dp) :: diff_grid(-1:1, -1:1),temp_grid(-1:1, -1:1)
      real(dp) :: eps,eps_1,eps_2,kp,length
      real(dp) :: denom_l,denom_r,denom_c
      real(dp) :: numer_l,numer_r,numer_c
      real(dp) :: diff, diff_test 
      real(dp) :: k_new,k_old,factor_x
      real(dp) :: step_a,step_kappa,dx
      real(dp) :: a_test,kappa_test 
      real(dp) :: f_1,f_2,q_1,q_2
      real(dp) :: x_1,x_2,x_c,x_l,x_r
      real(dp) :: y_c,y_l,y_r
      
      selectcase(trim(slw1_approach))
         case('kp-epsilon')
            !Compute target values
            length = maxval(slw1_length)  
            kp = get_slw_kp(Tref,pref,xsref,Tref,slw1_ngases)
            eps = get_slw_emissivity(Tref,pref,xsref,Tref,slw1_ngases,&
                                     length)  

            !Compute kappa iteratively
            counter = 0; k_old = 0.1_dp*kp; diff = 2._dp*iter_tol
            do while (diff.gt.iter_tol)
               counter = counter + 1
               if (counter.gt.max_iter) &
                  call shutdown('slw1_compute_ref: Maximum number of &
                                 &iterations exceeded')
               k_new = (kp/eps)*(1._dp - dexp(-k_old*length))
               diff = dabs((k_new - k_old)/(k_old + small))
               k_old = k_new
            enddo
          
            !Finish computing kappa and a
            kappa_out = k_new
            a_out = kp/k_new

         case('epsilon-epsilon')
            !Compute target values
            eps_1 = get_slw_emissivity(Tref,pref,xsref,Tref,&
                                       slw1_ngases,slw1_length(1))
            eps_2 = get_slw_emissivity(Tref,pref,xsref,Tref,&
                                       slw1_ngases,slw1_length(2))
            
            !Begin bisection method
            x_l = small; x_r = 1000._dp
            y_l = eps_1/eps_2 - (1._dp - dexp(-x_l*slw1_length(1)))/&
                                (1._dp - dexp(-x_l*slw1_length(2)))
            y_r = eps_1/eps_2 - (1._dp - dexp(-x_r*slw1_length(1)))/&
                                (1._dp - dexp(-x_r*slw1_length(2)))
            if (y_l*y_r.gt.0) &
                  call shutdown('slw1_compute_ref: Problem with &
                                &bisection method')
            do while(diff.gt.iter_tol)
               counter = counter + 1
               if (counter.gt.max_iter) &
                  call shutdown('slw1_compute_ref: Maximum number of &
                                 &iterations exceeded')
               
               x_c = (x_l + x_r)/2._dp
               y_c = eps_1/eps_2 - (1._dp - dexp(-x_c*slw1_length(1)))/&
                                   (1._dp - dexp(-x_c*slw1_length(2)))
               if (y_c*y_l.gt.0) then
                  x_l = x_c
                  y_l = y_c
               endif
               if (y_c*y_r.gt.0) then
                  x_r = x_c 
                  y_r = y_c 
               endif
               
               diff = dabs((y_l - y_r)/(y_c + small))
               
            enddo

            !Finish computing kappa and a
            kappa_out = x_c
            a_out = eps_1/(1._dp - dexp(-x_c*slw1_length(1)))

         case('F-F')
            !Surrogate names
            x_1 = slw1_position(1)
            x_2 = slw1_position(2)

            !Reference values
            f_1 = get_slw_uniform_flux(Tref,pref,xsref,Tref,&
                                       slw1_ngases,x_1)  
            f_2 = get_slw_uniform_flux(Tref,pref,xsref,Tref,&
                                       slw1_ngases,x_2)
                                       
            !Begin bisection method
            counter = 0 
            diff = 2._dp
            factor_x = 0.5_dp
            x_r = 10000._dp; x_l = x_r
            y_r = f_1/f_2 - (1._dp - 2._dp*expint(3, x_r*x_1))/&
                            (1._dp - 2._dp*expint(3, x_r*x_2))
            y_l = y_r
            do while (y_r*y_l.ge.0._dp)
               x_r = x_l; y_r = y_l
               x_l = x_r*factor_x
               y_l = f_1/f_2 - (1._dp - 2._dp*expint(3, x_l*x_1))/&
                                (1._dp - 2._dp*expint(3, x_l*x_2))
            enddo
     
            if (y_l*y_r.gt.0) &
               call shutdown('slw1_compute_ref: Problem with &
                                &bisection method')
                                
            do while(diff .gt. iter_tol)
               counter = counter + 1
               if (counter.gt.max_iter) &
                  call shutdown('slw1_compute_ref: Maximum number of &
                                 &iterations exceeded')
               
               x_c = (x_l + x_r)/2._dp
               y_c = f_1/f_2 - (1._dp - 2._dp*expint(3, x_c*x_1))/&
                               (1._dp - 2._dp*expint(3, x_c*x_2))
               if (y_c * y_l .lt. 0) then
                  x_r = x_c 
                  y_r = y_c
               endif
               if (y_c * y_r .lt. 0) then
                  x_l = x_c
                  y_l = y_c 
               endif
           
               diff = dabs((y_r - y_l)/(y_c + 1.e-6_dp))
               
            enddo

            !Finish computing kappa and a
            kappa_out = x_c
            a_out = f_1/(pi*Ib_function(Tref)*&
                    (1._dp - 2._dp*expint(3, x_c*x_1)))
            
            f_1 = pi*a_out*Ib_function(Tref)*(1._dp - 2._dp*expint(3, x_c*x_1))
            f_2 = pi*a_out*Ib_function(Tref)*(1._dp - 2._dp*expint(3, x_c*x_2))
            

         case('Q-Q')
            !Surrogate names
            x_1 = slw1_position(1)
            x_2 = slw1_position(2)
            length = maxval(slw1_length)
            x_l = small; x_r = 1000._dp

            !Reference values
            q_1 = get_slw_uniform_source(Tref,pref,xsref,Tref,&
                                         slw1_ngases,x_1,length)  
            q_2 = get_slw_uniform_source(Tref,pref,xsref,Tref,&
                                         slw1_ngases,x_2,length)
                                         
            !Begin bisection method
            counter = 0 
            diff = 2._dp
            factor_x = 0.5_dp
            x_r = 10000._dp; x_l = x_r
            numer_r = expint(2, x_r*x_1) + expint(2, x_r*(length - x_1))
            denom_r = expint(2, x_r*x_2) + expint(2, x_r*(length - x_2))
            y_r = q_1*denom_r - q_2*numer_r
            y_l = y_r

            do while (y_r*y_l.ge.0._dp)

               x_r = x_l; y_r = y_l
               x_l = x_r*factor_x
               numer_l = expint(2, x_l*x_1) + expint(2, x_l*(length - x_1))
               denom_l = expint(2, x_l*x_2) + expint(2, x_l*(length - x_2))
               y_l = q_1*denom_l - q_2*numer_l
            enddo
            
            if (y_l*y_r.gt.0) &
                  call shutdown('slw1_compute_ref: Problem with &
                                &bisection method')
                                
            do while(diff .gt. iter_tol)
               counter = counter + 1
               if (counter.gt.max_iter) &
                  call shutdown('slw1_compute_ref: Maximum number of &
                                 &iterations exceeded')
               
               x_c = (x_l + x_r)/2._dp
               numer_c = expint(2, x_c*x_1) + expint(2, x_c*(length - x_1))
               denom_c = expint(2, x_c*x_2) + expint(2, x_c*(length - x_2))
               y_c = q_1*denom_c - q_2*numer_c

               if (y_c * y_l .lt. 0) then
                  x_r = x_c 
                  y_r = y_c
               endif
               
               if (y_c * y_r .lt. 0) then
                  x_l = x_c
                  y_l = y_c 
               endif
       
               diff = dabs((x_r - x_l)/(x_c + small))
               if (diff.gt.(0.1_dp*iter_tol)) &
                  diff = dabs((y_r - y_l)/(y_c + small))
               
            enddo
            
            !Finish computing kappa and a
            kappa_out = x_c
            denom_c = (2._dp*pi*kappa_out*Ib_function(Tref)*&
              (expint(2, x_c*x_1) + expint(2, x_c*(length - x_1))))
            a_out = -q_1/(denom_c + small)
            
            q_1 = -(2._dp*pi*a_out*kappa_out*Ib_function(Tref)*&
                   (expint(2, x_c*x_1) + expint(2, x_c*(length - x_1))))
            q_2 = -(2._dp*pi*a_out*kappa_out*Ib_function(Tref)*&
                   (expint(2, x_c*x_2) + expint(2, x_c*(length - x_2))))
                   
                   
         case('Pattern Search')
             
            length = maxval(slw1_length)
         
            ! Pré varredura para encontrar um bom chute 
            diff_grid = huge(1.0_dp)
            step_a = 0.1_dp
            step_kappa = 0.1_dp
            nx = 50                ! Número de pontos da malha
            dx = length / (nx - 1._dp)

            do i = 1, 50
               a_test = 0.01_dp + i * step_a
               do j = 1, 50
                  kappa_test = 0.01_dp + j * step_kappa
                  
                  diff_test = 0._dp
                  do x = 1, nx
                     x_c = real(x-1,dp) * dx
                     eps_1 = get_slw_emissivity(Tref,pref,xsref,Tref,&
                                                slw1_ngases,x_c)
                     eps_2 = a_test * (1.0_dp - exp(-kappa_test * x_c))
                     diff_test = diff_test + (eps_1 - eps_2)**2 * dx
                  end do
                  
!write(*,*) a_test,kappa_test,diff_test
                  
                  if (diff_test < diff_grid(0,0)) then
                     diff_grid(0,0) = diff_test   ! Ponto central para comparação
                     a_out = a_test
                     kappa_out = kappa_test
                  endif
               end do
            end do
            
            ! Passos para a busca
            step_a = 0.01_dp
            step_kappa = 0.01_dp
                                         
            ! Começar a busca direta
            counter = 0 
            nx = 100
            dx = length / (nx - 1._dp)

            do while (diff_grid(0,0).gt.1.65e-2_dp)
               counter = counter + 1
               
               ! Testar 8 pontos ao redor do ponto atual (um quadrado)
               do i = -1, 1
                  do j = -1, 1
                  
                     if (i.eq.0 .and. j.eq.0) cycle ! Pula o ponto central
                     if (diff_grid(i,j).lt.huge(1.0_dp)) cycle !Pula este também, já foi calculado
                     
                     ! Atualizando os parâmetros
                     a_test = a_out + i * step_a
                     kappa_test = kappa_out + j * step_kappa
                     
                     diff_test = 0._dp
                     do x = 1, nx
                        x_c = real(x-1,dp) * dx
                        eps_1 = get_slw_emissivity(Tref,pref,xsref,Tref,&
                                                   slw1_ngases,x_c)
                        eps_2 = a_test*(1._dp - dexp(-kappa_test*x_c))
                        diff_test = diff_test + (eps_1 - eps_2)**2 * dx
                     end do
                     
write(*,*) a_test,kappa_test,diff_test                     

                     diff_grid(i,j) = diff_test
                  end do
               end do
               
               ! Encontrar o ponto com menor diferença
               do i = -1, 1
                  do j = -1, 1
                     if (i.eq.0 .and. j.eq.0) cycle ! Pular o ponto central
                     
                     ! Comparar todos e guardar para atualizar depois
                     if (diff_grid(i, j).lt.diff_grid(0,0)) then
                        x = i
                        y = j
                        diff_test = diff_grid(i, j)
                        a_test = a_out + i * step_a
                        kappa_test = kappa_out + j * step_kappa
                     end if
                     
                  end do
               end do
               
               !Shiftar a matriz e atualizar os parâmetros  
               temp_grid = huge(1.0_dp)

               do i = -1, 1
                  do j = -1, 1
                     ! nova posição = ponto antigo menos o deslocamento (x,y)
                     if ( (i-x >= -1) .and. (i-x <= 1) .and. &
                        (j-y >= -1) .and. (j-y <= 1) ) then
                        temp_grid(i-x, j-y) = diff_grid(i, j)
                     end if
                  end do
               end do

               diff_grid = temp_grid
               a_out = a_test
               kappa_out = kappa_test
               
               if (counter.gt.max_iter) &
                  call shutdown('slw1_compute_ref: Maximum number of &
                                 &iterations exceeded')
                                 
            enddo
            

         case default
            call shutdown('slw1_compute_ref: no valid option for &
                          &slw1_approach')
      endselect
   endsubroutine slw1_compute_ref

endmodule slw_functions



