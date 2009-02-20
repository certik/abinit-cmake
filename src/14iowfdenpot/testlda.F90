!{\src2tex{textfont=tt}}
!!****f* ABINIT/testlda
!! NAME
!! testlda
!!
!! FUNCTION
!! Test QPLDA or ABINIT LDA or KSS type file
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  Dtset<type(dataset_type)>=all input variables for this dataset, variables used :
!!     %accesswff=define the access mode
!!     %localrdwf=1 if each machine can access the file, 0 if only master
!!     %nsppol=only for checking purpose 
!!     %awtr=in case of time-reversal we need to know the number of occupied states
!!     %fixmom=to calculate new occ factors enforcing fixmom as magnetization  
!!  Dtfil<type(datafiles_type)>=datatype containing filenames
!!     %filkss=name of the KSS file
!!     %unkss=unit number associated to the KSS file
!!  MPI_enreg<type(MPI_type)>=datatype gathering information about the MPI parallelisation
!!
!! OUTPUT
!!  ibocc(Dtset%nsppol)=for each spin, the band index after which the occupations numbers 
!!    are smaller than tol8 (for screening calculations based on the Adler-wiser expression 
!!    with time-reversal or spectral method)
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotential
!!  nbnds_kss=number of bands contained in the KSS file
!!  ng_kss=number of plane waves in KSS file
!!  nsym=number of symmetries
!!  Hdr<hdr_type>=The abinit header.
!!
!! SIDE EFFECTS
!!  gvec_p(3,ng_kss)=
!!   In input pointer to integers, supposed to be not associated.
!!   In output it contains all the G vectors reported in the KSS file.
!!
!! NOTES
!!  Starting version 5.6, KSS files in single precision are not supported anymore.
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      hdr_clean,hdr_io,hdr_io_netcdf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine testlda(Dtset,Dtfil,nsym,nbnds_kss,ng_kss,mpsang,gvec_p,Hdr,MPI_enreg,&
 ibocc) ! Optional

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit, get_unit
 use m_errors, only : assert
#if defined HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_ETSF_IO
 use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
 use interfaces_13ionetcdf
 use interfaces_15gw
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: mpsang,nbnds_kss,ng_kss,nsym
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(in) :: Dtset
 type(Hdr_type),intent(out) :: Hdr
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 integer,pointer :: gvec_p(:,:)
 integer,intent(inout),optional :: ibocc(Dtset%nsppol)

!Local variables-------------------------------
#if defined HAVE_ETSF_IO
 integer,allocatable,target :: kg_k(:,:)
 type(ETSF_io_low_error) :: Error_data
 type(ETSF_electrons),target :: Electrons_folder
 type(ETSF_dims) :: Dims
 type(ETSF_basisdata),target :: Wave_folder
 logical :: lstat
 character(len=etsf_io_low_error_len) :: errmess
 real(dp),allocatable,target :: eigen(:),occ_vec(:)
#endif                                  
!scalars
 integer :: fform
 integer :: iatom,ib,ierr,ig,ii,ik,il,imin,ios,ispinor,isppol,itypat,master
 integer :: nel,nsh,nspinor,nsym_kss,rank,rdwr,spaceComm,untkss
 real(dp) :: efermi
 logical,parameter :: DEBUG=.FALSE.
 logical :: ltest
 character(len=500) :: msg
 character(len=fnlen) :: fname
!arrays
 real(dp),allocatable :: energyd(:),temp_en(:,:,:),temp_oc(:,:,:)
 character(len=80) :: tit(2)

! *************************************************************************

#if defined DEBUG_MODE 
 write(msg,'(a)')' testlda : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

#if !defined HAVE_ETSF_IO 
 if (Dtset%accesswff==3) then 
  write(msg,'(6a)')ch10,&
&  ' testlda: BUG - ',ch10,&
&  ' when accesswff==3, support for the ETSF I/O library ',ch10,&
&  ' must be compiled. Use --enable-etsf-io when configuring '
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end if
#endif
 !
 ! === Initialize MPI related quantities ===
 call xcomm_init   (MPI_enreg,spaceComm) 
 call xme_init     (MPI_enreg,rank     )          
 call xmaster_init (MPI_enreg,master   )  

 untkss=Dtfil%unkss

 if (rank==master.or.Dtset%localrdwf==1) then
  fname=Dtfil%filkss ; if (Dtset%accesswff==3) fname=TRIM(Dtfil%filkss)//'-etsf.nc'
  ! 
  ! === Try KSS file format by reading the header of the GS wavefunction file ===
  rdwr=1 
  select case (Dtset%accesswff)
  case (0)
   ! * Fortran unformatted file
   write(msg,'(3a)')ch10,&
&   ' testlda : testing Fortran KSS structure file ',TRIM(Dtfil%filkss)
   call wrtout(std_out,msg,'PERS')
   open(untkss,file=fname,status='old',form='unformatted',iostat=ios)
   if (ios/=0) then 
    write(msg,'(6a)')ch10,&
&    ' testlda : ERROR- ',ch10,&
&    ' opening file: ',TRIM(fname),' as old '
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
   call hdr_io(fform,Hdr,rdwr,untkss)
#if defined HAVE_ETSF_IO
   case (3)
    ! === NETCDF-ETSF file format ===
    write(msg,'(3a)')ch10,&
&    ' testlda : testing NETCDF-ETSF KSS file ',TRIM(fname)
    call wrtout(std_out,msg,'COLL')
    call etsf_io_low_open_read(untkss,fname,lstat,Error_data=Error_data)
    if (.not.lstat) then
     call etsf_io_low_error_to_str(errmess,Error_data)
     write(msg,'(4a)')ch10,' testlda: ERROR -',ch10,errmess(1:min(475,len(errmess)))
     call wrtout(std_out, msg, 'COLL') 
     call leave_new('COLL')
    end if
    call hdr_io_etsf(fform,Hdr,rdwr,untkss)
#endif
  case default
   write(msg,'(4a,i4)')ch10,&
&   ' testlda : BUG - ',ch10,&
&   ' wrong value for accesswff = ',Dtset%accesswff 
   call wrtout(std_out, msg, 'COLL') 
   call leave_new('COLL')
  end select

  if (fform==602) then
   write(msg,'(6a)')ch10,&
&   ' testlda : ERROR - ',ch10,&
&   ' starting v5.6, KSS files in single precision are not supported anymore,',ch10,&
&   ' Please, use an older version of abinit.'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if

  if (fform>=1.and.fform<=2) then
   ! STA or the QPLDA form (OBSOLETE)
    write(msg,'(6a,i4)')ch10,&
&   ' testlda : ERROR ',ch10,&
&   ' STA|QPLDA format not supported anymore ',ch10,&
&   ' fform = ',fform
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
  end if 

  if (fform/=502) then
    write(msg,'(6a,i4)')ch10,&
&   ' testlda : ERROR ',ch10,&
&   ' Found unknown file format ',ch10,&
&   ' fform = ',fform
   call wrtout(std_out,msg,'COLL')
   call leave_new('COLL')
  end if 

  ! === KSS file successfully read ===
  ! * Note that, in case of Fortran file, nsym is read from the second record
  nsym=Hdr%nsym
  write(msg,'(1x,47a)')('-',ii=1,47) 
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a)')' KSS abinit double precision form'
  call wrtout(std_out,msg,'COLL')
  write(msg,'(2a,a6,a,i3)')ch10,' generated by ABINIT ',Hdr%codvsn,' header version ',Hdr%headform
  call wrtout(std_out,msg,'COLL')
  if ( Hdr%headform/=23 .and. Hdr%headform/=34 .and. &
&      Hdr%headform/=40 .and. Hdr%headform/=41 .and. &
&      Hdr%headform/=42 .and. Hdr%headform/=44 .and. &
&      Hdr%headform/=53 .and. Hdr%headform/=56       &
&    ) then 
   write(msg,'(4a,i4)')ch10,&
&   ' testlda : ERROR - ',ch10,&
&   ' unknown  header version = ',Hdr%headform
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if 

  if (Dtset%accesswff==0) then 
   read(untkss) tit(1)
   read(untkss) tit(2)
   write(msg,'(a)')' title of file:'
   call wrtout(std_out,msg,'COLL')
   write(msg,'(1x,a79,a,1x,a79,a)')tit(1)(:79),ch10,tit(2)(:79),ch10
   call wrtout(std_out,msg,'COLL')
   read(untkss)nsym_kss,nbnds_kss,ng_kss,nsh,mpsang
   read(untkss) !(((symrel2(jj,ii,isym),ii=1,3),jj=1,3),isym=1,nsym_kss)
   read(untkss) !((tnons(i,isym),i=1,3),isym=1,nsym_kss)
   allocate(gvec_p(3,ng_kss))
   read(untkss)((gvec_p(ii,ig),ii=1,3),ig=1,ng_kss)
   ltest=ALL(gvec_p(:,1)==0)
   call assert(ltest,'First G-vector should be Gamma',__FILE__,__LINE__)
   nsym=nsym_kss
   if (Hdr%nsym/=nsym_kss) then 
    write(msg,'(4a)')ch10,&
&    ' testlda : WARNING - ',ch10,&
&    ' code does not use the original set of symmetries : Hdr%nsym/=nsym_kss '
    call wrtout(std_out,msg,'COLL')
   end if 
#if defined HAVE_ETSF_IO
  else if (Dtset%accesswff==3) then 
    ! TODO spin-orbit not treated, number of projectors not treated
   call etsf_io_dims_get(untkss,Dims,lstat,Error_data)
   nsym_kss =Dims%number_of_symmetry_operations
   nbnds_kss=Dims%max_number_of_states
   ng_kss=Dims%max_number_of_coefficients
   mpsang=Dims%max_number_of_angular_momenta
   allocate(gvec_p(3,ng_kss),kg_k(3,ng_kss))
   Wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k(:,:)
   call etsf_io_basisdata_get(untkss,Wave_folder,lstat,Error_data)
   gvec_p(:,:)=kg_k(:,:) ; deallocate(kg_k)
   ltest=ALL(gvec_p(:,1)==0)
   call assert(ltest,'First G-vector should be Gamma',__FILE__,__LINE__)
   ! nshells is not defined in the ETSF spefications but it is not used 
   nsh=0
#endif 
  end if 
  !  
  ! === Output important dimensions on the log file ===
  write(msg,'(a,i8)')' number of atomic species       ',Hdr%ntypat
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of atoms                ',Hdr%natom
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' Highest angular component +1   ',mpsang
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of symmetry operations  ',nsym
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of plane waves          ',ng_kss
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of irred k-points       ',Hdr%nkpt
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of bands                ',nbnds_kss
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of spinorial components ',Hdr%nspinor
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of spin polarisations   ',Hdr%nsppol
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of density components   ',Hdr%nspden
  call wrtout(std_out,msg,'COLL')
  write(msg,'(1x,47a)')('-',ii=1,47) 
  call wrtout(std_out,msg,'COLL')
  !  
  ! === Check the value of some variables ===
! MG TODO this will lead to array bound mismatch in some calculation, 
!    this is due to the awful treatment of symmetries done in outkss
  nspinor=Hdr%nspinor
  imin=MIN(nsym,nsym_kss)
  if (Hdr%nsym/=nsym_kss) then 
   write(msg,'(4a)')ch10,&
&   ' testlda : COMMENT - ',ch10,&
&   ' code does not use the original set of symmetries : Hdr%nsym/=nsym_kss'
   call wrtout(std_out,msg,'COLL') 
  end if 
  if (ANY(Hdr%symrel(:,:,1:imin)/=Dtset%symrel(:,:,1:imin))) then 
   write(msg,'(4a)')ch10,&
&   ' testlda : COMMENT - ',ch10,&
&   ' GW symmetries and symmetries inferred from input file do not agree '
   call wrtout(std_out,msg,'COLL') 
  end if 
  if (ANY(Hdr%tnons(:,1:imin)/=Dtset%tnons(:,1:imin))) then 
   write(msg,'(4a)')ch10,&
&   ' testlda : COMMENT - ',ch10,&
&   ' GW tnons and tnons inferred from input file do not agree '
   call wrtout(std_out,msg,'COLL') 
  end if 
  if (ANY(Hdr%symafm(1:imin)/=Dtset%symafm(1:imin))) then 
   write(msg,'(4a)')ch10,&
&   ' testlda : COMMENT - ',ch10,&
&   ' GW symafm and symafm inferred from input file do not agree '
   call wrtout(std_out,msg,'COLL') 
  end if 
  if (Hdr%nsppol/=Dtset%nsppol) then 
   write(msg,'(6a)')ch10,&
&   ' testlda ERROR - ',ch10,&
&   ' the value of nsspol read in the KSS file does not agree with the value ',ch10,&
&   ' specified in the input file, please modify your input file'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if 
  !  
  ! In case of parallelism over bands or Adler-Wiser with timereversal find the band 
  ! index separating the occupied and partially occupied from the empty states (for each spin)
  ! Each processor will store in memory the occupied states while the conduction 
  ! states will be shared between different processors
  !  
  if (PRESENT(ibocc)) then 
   write(msg,'(2a)')ch10,' testlda: reading occupation numbers to calculate ibocc '
   call wrtout(std_out,msg,'COLL')
   ! NOTE : the number of bands defined in the header and in the first 
   ! section of the KSS file might differ if the KSS file has been generated 
   ! using kssform 3 with nbandkss<nband.
   allocate(temp_oc(Hdr%nkpt,nbnds_kss,Hdr%nsppol))
   allocate(temp_en(Hdr%nkpt,nbnds_kss,Hdr%nsppol))
   temp_oc(:,:,:)=zero ; temp_en(:,:,:)=zero
   do isppol=1,Hdr%nsppol
    do ik=1,Hdr%nkpt
     do ib=1,Hdr%nband(ik+Hdr%nkpt*(isppol-1))
      ii = (ik-1)*Hdr%nband(ik) + (isppol-1)*Hdr%nkpt*Hdr%nband(ik) + ib 
      if (ib<=nbnds_kss) temp_oc(ik,ib,isppol)=Hdr%occ(ii)
     end do
    end do
   end do 
   allocate(energyd(1:nbnds_kss))

   if (Dtset%accesswff==0) then 
    ! Read eigenvalues from the KSS file in the FORTRAN format
    ! read(untkss) !tit(1) 
    ! read(untkss) !tit(2)
    ! read(untkss) !nsym_kss,nbandkss,npwkss,nsh,mpsang
    !read(untkss)   !(((symrel2(jj,ii,isym),ii=1,3),jj=1,3),isym=1,nsym_kss)
    !read(untkss)   !((tnons(i,isym),i=1,3),isym=1,nsym_kss)
    !read(untkss)   !((gvec(i,ig),i=1,3),ig=1,ngx)
    read(untkss)   !(shlim(i),i=1,nsh)
    read(untkss)   !((vkbsignd(il,is),il=1,mpsang),is=1,Hdr%ntypat)

    do isppol=1,Hdr%nsppol
     do ik=1,Hdr%nkpt

      do itypat=1,Hdr%ntypat
       do il=1,mpsang
        read(untkss) !vkbdb(:,itypat,il)
        read(untkss) !vkbdd(:,itypat,il)
       end do 
      end do 

      read(untkss) energyd(1:nbnds_kss)
      temp_en(ik,1:nbnds_kss,isppol)=energyd(1:nbnds_kss)

      do ib=1,nbnds_kss
       read(untkss) !wfgd(npwkss*nspinor)
       if (Hdr%usepaw==1) then
        do ispinor=1,nspinor
         do iatom=1,Hdr%natom
          read(untkss) !(cprjnk_k(ia,ibsp)%cp(:,1:cprjnk_k(ia,ibsp)%nlmn))
         end do
        end do
       end if
      end do

     end do !ik 
    end do !isppol
#if defined HAVE_ETSF_IO
   else if (Dtset%accesswff==3) then 
    allocate(eigen(nbnds_kss))
    eigen(:)=zero
    ! allocate(occ_vec(nbnds_kss))
    do isppol=1,Hdr%nsppol
     do ik=1,Hdr%nkpt
      Electrons_folder%eigenvalues%data1D => eigen
      Electrons_folder%eigenvalues__kpoint_access = ik
      Electrons_folder%eigenvalues__spin_access = isppol
!     Note : occupation numbers have been read from Hdr
!     Electrons_folder%occupations%data1D => occ_vec
!     Electrons_folder%occupations__kpoint_access = ik
!     Electrons_folder%occupations__spin_access = isppol
      call etsf_io_electrons_get(untkss,Electrons_folder,lstat,Error_data)
      if (.not.lstat) then
       call etsf_io_low_error_to_str(errmess,Error_data)
       write(msg,'(4a)')ch10,' testlda: ERROR -',ch10,errmess(1:min(475, len(errmess)))
       call wrtout(std_out,msg,'COLL')
      end if
      temp_en(ik,1:nbnds_kss,isppol)=eigen(1:nbnds_kss)
      if (DEBUG) write(*,*)isppol,ik,eigen(:)*Ha_eV 
     end do 
    end do 
    nullify(Electrons_folder%eigenvalues%data1D)
    deallocate(eigen) 
!   nullify(Electrons_folder%occupations%data1D) 
!   deallocate(occ_vec)
#endif
   end if
   !   
   ! Call fermi to calculate ibocc in case of semiconductors and occupation numbers in metals. 
   ! Note that occupations are zero if the KSS has been generated through a NSCF calculation 
   nel=zero
   call fermi(Hdr,nbnds_kss,Hdr%nkpt,Dtset%fixmom,Hdr%nsppol,Hdr%wtk,temp_en,temp_oc,nel,ibocc,efermi)

   if (Hdr%occopt>=3.and.Hdr%occopt<=7) then ! Metallic systems
    write(*,*)"Found metallic system"
    ibocc(:)=-HUGE(0)
    do isppol=1,Hdr%nsppol
     do ik=1,Hdr%nkpt
      do ib=1,nbnds_kss
       ! abs is used to avoid problems with negative occupation numbers 
       ! tol8 is chosen to be consistent with crho.F90 and density.F90
       ! TODO Should define a globlal variable (global to GW)
       if (ABS(temp_oc(ik,ib,isppol))<tol8) then 
        ibocc(isppol)=MAX(ib,ibocc(isppol)) ; EXIT 
       end if 
      end do 
     end do 
    end do 
   end if 
   if (rank==master) write(*,*)' Max band index for fully/partially occupied states : ',ibocc(:)
   deallocate(energyd)
   deallocate(temp_oc,temp_en)
  end if ! Present ibocc

  if (Dtset%accesswff==0) then 
   close(untkss)
#if defined HAVE_ETSF_IO
  else if (Dtset%accesswff==3) then 
   call etsf_io_low_close(untkss,lstat,Error_data)
#endif 
  end if 

 end if ! (rank==master.or.Dtset%localrdwf==1)
 !
 ! ============================================
 ! === Transmit data for parallel execution ===
 ! ============================================
 if (MPI_enreg%nproc>1.and.Dtset%localrdwf==0) then
  write(msg,'(2a)')ch10,' testlda : master is sending data...'
  call wrtout(std_out,msg,'COLL')
  call xcast_mpi(mpsang,master,spaceComm,ierr)
  call xcast_mpi(nbnds_kss,master,spaceComm,ierr)
  call xcast_mpi(ng_kss,master,spaceComm,ierr)
  call xcast_mpi(nsym,master,spaceComm,ierr)
  if (PRESENT(ibocc)) call xcast_mpi(ibocc,master,spaceComm,ierr)
  if (rank/=master.and.Dtset%localrdwf==0) allocate(gvec_p(3,ng_kss)) ! this proc has not read.
  call xcast_mpi(gvec_p,master,spaceComm,ierr)
#if defined MPI
  !TODO avoid preprocessor skipping hdr_comm inside the procedure if !MPI
  call hdr_comm(Hdr,master,rank,spaceComm)
#endif
  call leave_test(MPI_enreg) 
 end if

#if defined DEBUG_MODE 
 write(msg,'(a)')' testlda : exit'
 call wrtout(std_out,msg,'PERS') 
 call flush_unit(std_out)
 call leave_test(MPI_enreg) 
#endif

end subroutine testlda
!!***
