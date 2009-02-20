!{\src2tex{textfont=tt}}
!!****f* ABINIT/rdkss
!! NAME
!! rdkss
!!
!! FUNCTION
!! Read a _KSS file
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MM, XG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! Cprj_ibz(natom,nspinor*mband*mkmem*nsppol*psps%usepaw)<type(cprj_type)>=
!!  projected input wave functions <Proj_i|Cnk> with all NL projectors (only for PAW)
!! Dtfil<type(datafiles_type)>=datatype containing filenames, variables used
!!    %filkss=name of the KSS file 
!!    %unkss=unit number associated to the KSS file.
!! Dtset <type(dataset_type)>=input variables for this dataset, used variables
!!    %accesswff=option definig the file format of the KSS file
!!    %rprimd_orig(3,3)=value of rprimd coming from the input file
!!    %localrdwf=input variable (for parallel case) if 1, the KSS file is local to each machine
!!    %accesswff=1 for normal Fortran IO, 3 for ETSF-IO mode
!!    %prtvol=input variable defining the verbosity of the output
!!    %usepaw=1 in case of PAW calculations
!!    %charge=Additional charge to be added to the system 
!! (FIXME MG: this is not safe!!! one should calculate the number of electrons from the occupations numbers
!! mpsang= 1+maximum angular momentum for nonlocal pseudopotential
!! natom=Number of atoms
!! nspinor=number of spinorial components
!! nbndsA=Number of bands required (<=nbandkss)
!! npwvec=Max Number of G vectors required
!! nkibz=Number of irreducible k-points
!! nsym_gw=maximum number of symmetry operations
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! ntypat=maximum number of types of atoms
!! npwwfn=maximum number of planewaves
!! my_minb,my_maxb=index of the min and max band stored and treated by this proc
!! nbvw= maximum number of fully and partially occupied states over spin  (used for the new algorithm)
!! Pawtab(psps%ntypat*psps%usepaw) <type(pawtab_type)>=paw tabulated starting data obtained from PAW datasets
!! MPI_enreg=informations about MPI parallelization
!!
!! OUTPUT
!! en(nkibz,nbndsA,nsppol)=KS energies for each k point and each band
!! gvec(3,npwvec)=integer coordinates of G vectors
!! nel=number of electrons
!! occ(nkibz,nbndsA,nsppol)=occupation numbers for each k point, each band and each spin
!! symrec(3,3,nsym_gw)=symmetry operations in reciprocal space
!! tit(2)(len=80)= titles in the KSS file
!! tnons(3,nsym_gw)=fractional translations
!! vkb(npwwfn,ntypat,mpsang,nkibz)=KB projector function
!! vkbd(npwwfn,ntypat,mpsang,nkibz)=derivative of the KB projector function in reciprocal space
!! vkbsign(mpsang,ntypat)=sign of each KB dyadic product
!! wf(npwwfn,my_minb:my_maxb,nkibz,nsppol)=wavefunctions in G space, for each band, each k point and spin
!! wf_val(npwwfn,nbvw,nkibz,nsppol)=(OPTIONAL) array containing all the states whose occupation is less than tol8 
!!  used only if called from screening
!! 
!! NOTES 
!! For historical reasons, in case of Fortran files, symmetry operations and fractional translations 
!! are not read from the  Hdr but directly from the KSS file. 
!! This behavior should be changed but most of the automatic tests should be updated
!!
!! PARENTS
!!      rdm,screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rdkss(Dtfil,Dtset,Pawtab,nsym_gw,nbndsA,nbvw,nkibz,npwvec,nspinor,nsppol,npwwfn,tit,symrec,gvec,&
& en,occ,wf,Cprj_ibz,ntypat,natom,mpsang,tnons,vkbsign,vkb,vkbd,nel,MPI_enreg,my_minb,my_maxb,&
& wf_val) ! Optional arguments 

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit
#if defined HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_ETSF_IO
 use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13io_mpi
 use interfaces_13ionetcdf
 use interfaces_14iowfdenpot, except_this_one => rdkss
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang,my_maxb,my_minb,natom,nbndsA,nbvw,nkibz,npwvec
 integer,intent(in) :: npwwfn,nspinor,nsppol,nsym_gw,ntypat
 integer,intent(out) :: nel
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(in) :: Dtset
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 integer,intent(out) :: gvec(3,npwvec),symrec(3,3,nsym_gw)
 real(dp),intent(out) :: en(nkibz,nbndsA,nsppol),occ(nkibz,nbndsA,nsppol)
 real(dp),intent(out) :: tnons(3,nsym_gw),vkb(npwwfn,ntypat,mpsang,nkibz)
 real(dp),intent(out) :: vkbd(npwwfn,ntypat,mpsang,nkibz)
 real(dp),intent(out) :: vkbsign(mpsang,ntypat)
 complex(gwpc),intent(inout),optional,target :: wf_val(npwwfn*nspinor,nbvw,nkibz,nsppol)
 complex(gwpc),intent(inout),target :: wf(npwwfn*nspinor,my_minb:my_maxb,nkibz,nsppol)
 character(len=80),intent(out) :: tit(2)
 type(Cprj_type),intent(out) :: Cprj_ibz(natom,nspinor*nbndsA*nkibz*nsppol*Dtset%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Dtset%ntypat*Dtset%usepaw)

!Local variables-------------------------------
#if defined HAVE_ETSF_IO
 ! TODO check GWPC in case of ETSF_IO
 integer :: icg,option,formeig,tim_rwwf,nband_disk,optkg
 integer :: mband,mcg
 type(ETSF_main),target :: Main_folder
 type(ETSF_io_low_error) :: Error_data
 type(ETSF_dims) :: Dims
 type(ETSF_basisdata),target :: Wave_folder
 type(ETSF_gwdata),target :: GW_folder
 type(ETSF_electrons),target :: Electrons_folder
 type(Wffile_type) :: Wff
 logical :: lstat
 character(len=etsf_io_low_error_len) :: errmess
 character(len=fnlen) :: file_etsf
 integer,allocatable,target :: kg_k(:,:)
 integer,allocatable,target :: vkbsign_int(:,:)
 real(dp),allocatable,target :: vkb_tgt(:,:,:),vkbd_tgt(:,:,:)
 real(dp),allocatable,target :: cg(:,:),eigen(:),occ_vec(:)
#endif                                  
!scalars
 integer :: accesswff,enough,fform,ia,iatom,ib,ibg,ibp,ibsp,ibsp1,ibsp2,ierr,ig
 integer :: ii,ik,il,ilmn,ios,ipa,ipsp,is,ish,ispinor,isppol,isym,it,itypat
 integer :: j0lmn,jj,jlmn,klmn,localrdwf,master,mpsang_,nbandkss,nela,nprocs
 integer :: npwkss,nsh,nsym2,osl,pad1,pad2,prtvol,rank,rdwr,spaceComm,spad_kss
 integer :: spad_wfn,spinor_shift,untkss,usepaw
 real(dp) :: cinf,cinf1,csup,csup1,dum,einf,einf1,esup,esup1,sij
 complex(dpc) :: cdum
 logical,parameter :: DEBUG=.FALSE.
 logical :: i_read,master_casts,prior_to_55=.FALSE.,read_val,test
 character(len=100) :: frmt
 character(len=500) :: msg
 type(Hdr_type) :: Hdr
!arrays
 integer,allocatable :: shlim(:),symrec2(:,:,:),symrel2(:,:,:)
 real(dp) :: paw_ovlp(2)
 real(dp),allocatable :: energyd(:),temp_en(:,:),tnons2(:,:),vkbdb(:,:,:)
 real(dp),allocatable :: vkbdd(:,:,:),vkbsignd(:,:)
 complex(dpc),allocatable :: wfgd(:)
 complex(gwpc),pointer :: wfg1(:),wfg2(:)
 character(len=6) :: tag_spin(2)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' rdkss : enter '
 call wrtout(std_out,msg,'PERS') 
 call flush_unit(std_out)
#endif

 prtvol   = Dtset%prtvol
 accesswff= Dtset%accesswff 
 usepaw   = Dtset%usepaw
 localrdwf= Dtset%localrdwf ! 1 if each machine has access, 0 if only master

 untkss   = Dtfil%unkss     

#if !defined HAVE_ETSF_IO 
 if (accesswff==3) then 
  write(msg,'(6a)')ch10,&
&  ' rdkss : BUG - ',ch10,&
&  ' when accesswff==3, support for the ETSF I/O library ',ch10,&
&  ' must be compiled. Use --enable-etsf-io when configuring '
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end if 
#endif
!
!=== Initialize MPI related quantities ===
 call xcomm_init   (MPI_enreg,spaceComm) 
 call xme_init     (MPI_enreg,rank     )          
 call xmaster_init (MPI_enreg,master   )  
 call xproc_max    (nprocs,ierr)
!
!=== Check if valence and conduction are stored in two different arrays ===
 read_val=.FALSE.
 if (PRESENT(wf_val)) then 
  read_val=.TRUE.
  write(msg,'(4a)')ch10,&
&  ' rdkss : taking advantage of time reversal symmetry ',ch10,&
&  ' occupied and unoccupied states stored in two different arrays '
  call wrtout(std_out,msg,'COLL')
 end if 
!
!=== Define I/O for parallel execution ===
 i_read       = (localrdwf==1.or.rank==master)
 master_casts = (localrdwf==0.and.nprocs>1)
!
!=== Read file according to the fileformat ===
 if (i_read) then
  rdwr=1
  if (accesswff==0) then 
!  * Formatted Fortran File
   write(msg,'(3a)')ch10,&
&   ' rdkss : reading FORTRAN Kohn-Sham structure file ',TRIM(Dtfil%filkss)
   call wrtout(std_out,msg,'COLL')
   open(unit=untkss,file=Dtfil%filkss,form='unformatted',status='old',iostat=ios)
   if (ios/=0) then 
    write(msg,'(6a)')ch10,&
&    ' rdkss: ERROR- ',ch10,  &
&    ' opening file: ',TRIM(Dtfil%filkss),' as old'
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
   call hdr_io(fform,Hdr,rdwr,untkss)
#if defined HAVE_ETSF_IO
  else if (Dtset%accesswff==3) then 
!  * NETCDF-ETSF file format 
   if (usepaw==1) call not_implemented_rdkss('PAW+ETSF-IO')
   file_etsf=TRIM(Dtfil%filkss)//'-etsf.nc'
   write(msg,'(3a)')ch10,&
&   ' rdkss : reading NETCDF-ETSF Kohn-Sham structure file ',TRIM(file_etsf)
   call wrtout(std_out,msg,'COLL')
   call etsf_io_low_open_read(untkss,file_etsf,lstat,Error_data=Error_data)
   if (.not.lstat) then
    call etsf_io_low_error_to_str(errmess,Error_data)
    write(msg,'(4a)')ch10,' rdkss : ERROR -',ch10,errmess(1:min(475,len(errmess)))
    call wrtout(std_out,msg,'COLL') 
    call leave_new('COLL')
   end if
   call hdr_io_etsf(fform,Hdr,rdwr,untkss)
   if (fform==602) call not_implemented_rdkss('single precision + ETSF-NETCDF')
#endif
  end if 

  if (fform/=502.and.fform/=602) then
   write(msg,'(6a,i8)')ch10,&
&   ' rdkss : ERROR-  ',ch10,&
&   ' file ',TRIM(Dtfil%filkss),' has fform/=502,602 ',fform
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if
  if (fform==602) then
   write(msg,'(6a)')ch10,&
&   ' rdkss : ERROR - ',ch10,&
&   ' starting v5.6, KSS files in single precision are not supported anymore,',ch10,&
&   ' Please, use an older version of abinit.'
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
  end if

  if (fform==502) write(msg,'(2a)')' Kohn Sham Structure double precision file',ch10
  call wrtout(std_out,msg,'COLL')
! 
! === Output the header of the GS wavefunction file ===
  rdwr=4
  if (accesswff==0) then 
   call hdr_io(fform,Hdr,rdwr,std_out)
#if defined HAVE_ETSF_IO
  else if (accesswff==3) then 
   call hdr_io_etsf(fform,Hdr,rdwr,std_out)
#endif
  end if 
! 
! === Perform some check ===
  if (nspinor/=Hdr%nspinor.or.ntypat/=Hdr%ntypat.or.Hdr%natom/=natom.or.Hdr%nkpt/=nkibz.or.Hdr%nsppol/=nsppol) then 
   write(msg,'(3a)')ch10,' rdkss : BUG - inconsistent data1',ch10
   call wrtout(std_out,msg,'COLL')
   write(std_out,*)Hdr%nspinor,nspinor,Hdr%ntypat,ntypat,Hdr%natom,natom,Hdr%nkpt,nkibz,Hdr%nsppol,nsppol 
   call leave_new('COLL')
  end if 
  if (Hdr%usepaw/=Dtset%usepaw) stop "Hdr%usepaw/=Dtset%usepaw"
  tag_spin(:)=(/'      ','      '/) ; if (nsppol==2) tag_spin(:)=(/' UP   ',' DOWN '/) 
! 
! === Save occupation numbers ===
! NOTE: Hdr%nband(nkibz*nsppol). Hdr%occ(bantot) where bantot is the total number 
! of bands (i.e sum of nband on all kpts and spins). In the GW part the number 
! of bands must be the same for each k-point and spin
  do isppol=1,nsppol
   do ik=1,nkibz
    do ib=1,Hdr%nband(ik+nkibz*(isppol-1))
     ii= (ik-1)*Hdr%nband(ik) +(isppol-1)*nkibz*Hdr%nband(ik) + ib 
     if (ib<=nbndsA) occ(ik,ib,isppol)=Hdr%occ(ii)
    end do
   end do
  end do 
! === Test spin-orbit characteristic ===
  if (Hdr%headform<53) then 
!  * Old format, previous to version 5.5, now pspo is obsolete and has been substituted by so_psp.
   test=ALL(Hdr%pspso(1:ntypat)==1)
   call assert(test,'pspso/=1 value not programmed',__FILE__,__LINE__)
  else 
!  * New format containing so_psp
   test=ALL(Hdr%so_psp(1:Hdr%npsp)==1)
   call assert(test,'so_psp/=1 value not programmed',__FILE__,__LINE__)
  end if

  nel=0
  do ia=1,natom
   it=Hdr%typat(ia)
   nela=Hdr%znucltypat(it)  ! if it is not a pseudoatom, add all electron
   do ipa=1,Hdr%npsp
    if (Hdr%znuclpsp(ipa)==Hdr%znucltypat(it)) nela=Hdr%zionpsp(ipa)
   end do
   nel=nel+nela
  end do
  nel=nel-dtset%charge  !FIXME MG: this part has to be done in a safer way!
 end if !i_read
!
!=== End of treatment of the header ===
!* In case of parallel execution broadcast data
 if (master_casts) then 
  call xcast_mpi(occ,master,spaceComm,ierr)
  call xcast_mpi(nel,master,spaceComm,ierr)
  call xcast_mpi(fform,master,spaceComm,ierr)
! call leave_test(MPI_enreg)
 end if

 tit(1)='NoNe' ; tit(2)='NoNe'
 if (i_read) then 
  if (accesswff==0) then 
   read(untkss) tit(1)
   read(untkss) tit(2)
   write(msg,'(a,1x,a79,a,1x,a79)')ch10,tit(1)(:79),ch10,tit(2)(:79)
   call wrtout(std_out,msg,'COLL')
   read(untkss) nsym2,nbandkss,npwkss,nsh,mpsang_
#if defined HAVE_ETSF_IO
  else if (accesswff==3) then ! TODO spin-orbit not treated, number of projectors not treated
   call etsf_io_dims_get(untkss,Dims,lstat,Error_data)
   nsym2   =Dims%number_of_symmetry_operations
   nbandkss=Dims%max_number_of_states                                             
   npwkss  =Dims%max_number_of_coefficients
   mpsang_ =Dims%max_number_of_angular_momenta
!  nshells not defined in NETCDF-ETSF-IO, anyway not used
   nsh=0
#endif
  end if 

  write(msg,'(a,i8)')' number of electrons                    ',nel
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of symmetries without inversion ',nsym2
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of bands                        ',nbandkss
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of plane waves                  ',npwkss
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8)')' number of shells                       ',nsh
  call wrtout(std_out,msg,'COLL')
  write(msg,'(a,i8,a)')' maximum angular momentum for PSPs+1    ',mpsang_,ch10
  call wrtout(std_out,msg,'COLL')
! Enforcing equal nsym would be better, but then we cant change on the fly the operations 
  test=(mpsang_/=mpsang .or. nsym_gw>nsym2 .or. npwvec>npwkss .or. npwwfn>npwkss .or. nbndsA>nbandkss)
  call assert(.not.test,'Inconsistent data',__FILE__,__LINE__)
! 
! === Here we deal with the reading of symmetry operations ===
! * Get symmetry operations in reciprocal space (inverse-transpose)
  allocate(symrel2(3,3,nsym2),tnons2(3,nsym2),symrec2(3,3,nsym2))
  if (accesswff==0) then 
   read(untkss) (((symrel2(jj,ii,isym),jj=1,3),ii=1,3),isym=1,nsym2)
   read(untkss) ((tnons2(ii,isym),ii=1,3),isym=1,nsym2)
   do isym=1,nsym2
    call mati3inv(symrel2(:,:,isym),symrec2(:,:,isym))
   end do 
   symrec(:,:,1:nsym_gw)=symrec2(:,:,1:nsym_gw)
   tnons(:,1:nsym_gw)=tnons2(:,1:nsym_gw)
#if defined HAVE_ETSF_IO
  else if (accesswff==3) then 
!  Here we read all the operations from the Hdr
   allocate(symrec2(3,3,nsym2),symrel2(3,3,nsym2),tnons2(3,nsym2))
   do isym=1,nsym2
    symrel2(:,:,isym)=Hdr%symrel(:,:,isym)
    call mati3inv(Hdr%symrel(:,:,isym),symrec2(:,:,isym))
   end do 
   symrec(:,:,1:nsym_gw)=symrec2(:,:,1:nsym_gw)
   tnons(:,1:nsym_gw)=Hdr%tnons(:,1:nsym_gw)
#endif
  end if
  deallocate(symrec2,tnons2)

  write(msg,'(2a)')ch10,' Rotations                           Translations '
  call wrtout(std_out,msg,'COLL')
  do isym=1,nsym_gw
   write(msg,'(1x,(3(3i3,1x),4x,3(f11.7,1x)))')symrel2(:,:,isym),tnons(:,isym)
   call wrtout(std_out,msg,'COLL')
  end do 
 end if !i_read

 if (master_casts) then
  call xcast_mpi(nsym2,   master,spaceComm,ierr)
  call xcast_mpi(nbandkss,master,spaceComm,ierr)
  call xcast_mpi(npwkss,  master,spaceComm,ierr)
  call xcast_mpi(nsh,     master,spaceComm,ierr)
  call xcast_mpi(symrec,  master,spaceComm,ierr)
  call xcast_mpi(tnons,   master,spaceComm,ierr)
! call leave_test(MPI_enreg) 
 end if
!
!=== Read G-vectors and shell ===
 if (i_read) then
  if (accesswff==0) then 
   allocate(shlim(nsh))
   read(untkss) ((gvec(ii,ig),ii=1,3),ig=1,npwvec)
   read(untkss) (shlim(ii),ii=1,nsh)
   test=ALL(gvec(:,1)==0)
   call assert(test,'First G-vector should be Gamma',__FILE__,__LINE__)
   write(msg,'(3a)')ch10,'  shell limit G-vectors [reduced coordinates] ',ch10
   call wrtout(std_out,msg,'COLL')
   if (rank==master) then
    enough=0 ; osl=1
    do ish=1,nsh
     enough=enough+1 
     if (enough>20.and.prtvol==0) then 
      write(msg,'(a)')' rdkss : prtvol=0, stop printing more G-vector information'
      call wrtout(std_out,msg,'COLL') ; EXIT
     end if
     if (shlim(ish)>npwvec) osl= shlim(ish)+1
     write(*,'(i6,i6,4(i5," (",3i3,")"),20(/12x,4(i5," (",3i3,")")))')&
&     ish,shlim(ish),(ig,(gvec(ii,ig),ii=1,3),ig=osl,shlim(ish))
     osl=shlim(ish)+1
    end do
    write(*,*)
   end if
!  
!  === Read Kleynmann-Bylander sign and matrix elements ===
   allocate(vkbsignd(mpsang,ntypat))
   read(untkss) ((vkbsignd(il,is),il=1,mpsang),is=1,ntypat)
   vkbsign(1:mpsang,1:ntypat)=vkbsignd(1:mpsang,1:ntypat)
   deallocate(vkbsignd)
#if defined HAVE_ETSF_IO
  else if (accesswff==3) then 
   allocate(kg_k(3,npwkss))
!  Read reduced_coordinates_of plane waves (k-point independent)
   Wave_folder%reduced_coordinates_of_plane_waves%data2D => kg_k(:,:)
!  Wave_folder%red_coord_pw__kpoint_access = ikpt
   call etsf_io_basisdata_get(untkss,Wave_folder,lstat,Error_data)
   gvec(:,:)=kg_k(:,:)
   nullify(Wave_folder%reduced_coordinates_of_plane_waves%data2D)
   deallocate(kg_k)
   if (fform==602) call not_implemented_rdkss(' single precision + ETST_IO not implemented')
!  Read Kleynmann-Bylander sign and matrix elements (stored in double precision variables)
!  Note that this is integer, while output is real !
   allocate(vkbsign_int(mpsang,ntypat)) 
   vkbsign_int(:,:)=0
   GW_folder%kb_formfactor_sign%data2D => vkbsign_int
   call etsf_io_gwdata_get(untkss,GW_folder,lstat,Error_data)
   vkbsign(1:mpsang,1:ntypat)=REAL(vkbsign_int(1:mpsang,1:ntypat),dp)
   nullify(GW_folder%kb_formfactor_sign%data2D)
   deallocate(vkbsign_int)
#endif
  end if 
 end if !i_read

 if (master_casts) then
  call xcast_mpi(vkbsign,master,spaceComm,ierr)
  call xcast_mpi(gvec,master,spaceComm,ierr)
! call leave_test(MPI_enreg) 
 end if

 write(frmt,*)' (2a,',mpsang*ntypat,'(1x,f4.1)) '
 write(msg,frmt)ch10,' vkbsign : ',vkbsign(:,:)
 call wrtout(std_out,msg,'COLL')
!
!=== Allocate arrays for KB form factors and wavefunctions ===
 allocate(vkbdb(npwkss,ntypat,mpsang))
 allocate(vkbdd(npwkss,ntypat,mpsang))
 allocate(energyd(nbandkss))
 allocate(wfgd(npwkss*nspinor))
!
!Here it is possible to reduce the number of BCAST
!Moreover vkb and vkbd do not depend on the spin
 write(msg,'(2a)')ch10,' k       eigenvalues [eV]'
 call wrtout(ab_out,msg,'COLL') 
 call wrtout(std_out,msg,'COLL')

 ibg=0
 do isppol=1,nsppol
  do ik=1,nkibz

!  === Read Kleynmann-Bylander matrix elements and derivatives for this k-point ===
   if (accesswff==0) then 
    do is=1,ntypat
     do il=1,mpsang
      if (i_read) then
       read(untkss) vkbdb(:,is,il)
       read(untkss) vkbdd(:,is,il)
      end if
      if (master_casts) then
       call xcast_mpi(vkbdb(:,is,il),master,spaceComm,ierr)
       call xcast_mpi(vkbdd(:,is,il),master,spaceComm,ierr)
      end if
      vkb (1:npwwfn,is,il,ik)=vkbdb(1:npwwfn,is,il)
      vkbd(1:npwwfn,is,il,ik)=vkbdd(1:npwwfn,is,il)
     end do !il
    end do !is
#if defined HAVE_ETSF_IO
   else if (accesswff==3) then 
    if (i_read) then
     allocate(vkb_tgt (npwkss,mpsang,ntypat))
     allocate(vkbd_tgt(npwkss,mpsang,ntypat))
     GW_folder%kb_coeff__kpoint_access    =ik
     GW_folder%kb_coeff_der__kpoint_access=ik
     GW_folder%kb_formfactors%data3D => vkb_tgt
     GW_folder%kb_formfactor_derivative%data3D => vkbd_tgt
     call etsf_io_gwdata_get(untkss,GW_folder,lstat,Error_data)
!    FIXME note that mpsang and ntypat are inverted wrt to vkbd
!    should always use this ordering since it should be faster
     do is=1,ntypat
      do il=1,mpsang
       vkbdb(:,is,il)=vkb_tgt (:,il,is)
       vkbdd(:,is,il)=vkbd_tgt(:,il,is)
      end do
     end do
     nullify(GW_folder%kb_formfactors%data3D)
     deallocate(vkb_tgt)
     nullify(GW_folder%kb_formfactor_derivative%data3D)
     deallocate(vkbd_tgt)
    end if
    if (master_casts) then
     call xcast_mpi(vkbdb(:,:,:),master,spaceComm,ierr)
     call xcast_mpi(vkbdd(:,:,:),master,spaceComm,ierr)
    end if
    vkb (1:npwwfn,:,:,ik)=vkbdb(1:npwwfn,:,:)
    vkbd(1:npwwfn,:,:,ik)=vkbdd(1:npwwfn,:,:)
#endif         
   end if 
   if (DEBUG) write(*,*)' vkb  : ik = ',ik,vkb (:,:,:,ik)
   if (DEBUG) write(*,*)' vkbd : ik = ',ik,vkbd(:,:,:,ik)
!  
!  === Each proc stores only the sub-block of eigenvalues required by the user ===
   if (accesswff==0) then 
    if (i_read) read(untkss) energyd(1:nbandkss)
    if (master_casts) call xcast_mpi(energyd,master,spaceComm,ierr)
    en(ik,1:nbndsA,isppol)=energyd(1:nbndsA)
    if (rank==master) then 
     if (nsppol==2) then 
      write(std_out,'(i3,a,10f7.2/50(10x,10f7.2/))') ik,tag_spin(isppol),(en(ik,ib,isppol)*Ha_eV,ib=1,nbndsA)
      write(ab_out,'(i3,a,10f7.2/50(10x,10f7.2/))') ik,tag_spin(isppol),(en(ik,ib,isppol)*Ha_eV,ib=1,nbndsA)
     else 
      write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik,(Ha_eV*en(ik,ib,isppol),ib=1,nbndsA)
      write(ab_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik,(Ha_eV*en(ik,ib,isppol),ib=1,nbndsA)
     end if
    end if
!   === Each proc stores only the sub-block of wavefunction required by user ===
    ibsp=0

    do ib=1,nbndsA
     if (i_read) then 
      read(untkss) wfgd(1:npwkss*nspinor)
!     write(78,'(2f8.4)')  wfgd(1:npwkss*nspinor)
      if (usepaw==1) then 
!      FIXME here I should used dimlmn and check ordering used in ctocprj
       do ispinor=1,nspinor
        ibsp=ibsp+1
        do ia=1,natom
         ii=Cprj_ibz(ia,ibg+ibsp)%nlmn 
         read(untkss) Cprj_ibz(ia,ibg+ibsp)%cp(:,1:ii)
!        write(78,'(f8.4)')(Cprj_ibz(ia,ibg+ibsp)%cp(:,1:ii))
        end do
       end do
      end if
     end if
     if (master_casts) then 
!     MG this part is awful, should be rewritten!
      call xcast_mpi(wfgd,master,spaceComm,ierr)
      if (usepaw==1) then 
       do ispinor=1,nspinor
        jj=ispinor-1
        do ia=1,natom
         ii=Cprj_ibz(ia,ibg+ib+jj)%nlmn
         call xcast_mpi(Cprj_ibz(ia,ibg+ib+jj)%cp(1:2,1:ii),master,spaceComm,ierr)
        end do
       end do
      end if
     end if

     if (ib>=my_minb.and.ib<=my_maxb) then 
!     Get a slice of each spinor.
      do ispinor=1,nspinor
       spad_kss=(ispinor-1)*npwkss
       spad_wfn=(ispinor-1)*npwwfn
       wf(spad_wfn+1:spad_wfn+npwwfn,ib,ik,isppol)=wfgd(spad_kss+1:spad_kss+npwwfn)
      end do
     end if
!    In the new algorithm each proc stores the valence states in a separate array
!    FIXME should optimize the communication should check the localrdwf==0 case
     if (read_val.and.ib<=nbvw) then 
      do ispinor=1,nspinor
       spad_kss=(ispinor-1)*npwkss
       spad_wfn=(ispinor-1)*npwwfn
       wf_val(spad_wfn+1:spad_wfn+npwwfn,ib,ik,isppol)=wfgd(spad_kss+1:spad_kss+npwwfn)
      end do
     end if
    end do

    if (i_read) then ! Skip the remaining part
     do ib=nbndsA+1,nbandkss 
      read(untkss) wfgd(1:npwkss*nspinor)
!     write(78,*)   wfgd(1:npwkss*nspinor)
      if (usepaw==1) then 
       do ispinor=1,nspinor
        do ia=1,natom 
         read(untkss) !(Cprj_ibz(ia,ibsp)%cp
        end do
       end do
      end if
     end do
    end if
#if defined HAVE_ETSF_IO
   else if (Dtset%accesswff==3) then 
    if (nspinor==2) call not_implemented_rdkss("spinor==2")
!   icg=0 ; option=1 
    formeig=0 !; tim_rwwf=0 ; optkg=0
!   Wff%master=0 ; Wff%me=MPI_enreg%me ; Wff%unwff=untkss ; Wff%accesswff=3
    mcg=nbandkss*npwkss ; mband=nbandkss
!   !this leads to a crash in rwwf
!   nband_disk=mband
!   !cannot reshape since cg is real while wfg is complex 
!   ! but it is possibile to reshape occ_vec
    allocate(cg(2,mcg),eigen((2*mband)**formeig*mband),occ_vec(mband))
!   allocate(kg_k(3,optkg*npwkss))
    if (i_read) then
!    
!    We get eigenvalues, occupations and coefficients 
!    
!    TODO here it is better to use rwwf but it seems there is a problem in wffwritenpwrec 
!    since in case of ETSF_IO the npw record is not correctly written on the nc file, ask Damien 
!    For the time being Im using low level subroutines, they should be replaced
!    by rwwf but we have to fix the problem in wffwritenpwrec
!    
!    call rwwf(cg,eigen,formeig,Hdr%headform,icg,ik,isppol,kg_k,mband,mcg,&
!    & nbandkss,nband_disk,npwkss,Hdr%nspinor,occ_vec,option,optkg,tim_rwwf,Wff)
!    
!    BEGIN LOW LEVEL
     Electrons_folder%eigenvalues%data1D => eigen
     Electrons_folder%eigenvalues__kpoint_access = ik
     Electrons_folder%eigenvalues__spin_access = isppol
     Electrons_folder%occupations%data1D => occ_vec
     Electrons_folder%occupations__kpoint_access = ik
     Electrons_folder%occupations__spin_access = isppol
     call etsf_io_electrons_get(untkss,Electrons_folder,lstat,Error_data)
     if (.not.lstat) then
      call etsf_io_low_error_to_str(errmess,Error_data)
      write(msg,'(4a)')ch10,' rdkss: ERROR -',ch10,errmess(1:min(475, len(errmess)))
      call wrtout(std_out,msg,'COLL')
     end if
!    === Get the coefficients_of_wavefunctions ===
     Main_folder%coefficients_of_wavefunctions%data2D => cg(:,:)
     Main_folder%wfs_coeff__kpoint_access = ik
     Main_folder%wfs_coeff__spin_access   = isppol
     call etsf_io_main_get(untkss,Main_folder,lstat,Error_data)
     if (.not. lstat) then
      call etsf_io_low_error_to_str(errmess,Error_data)
      write(msg,'(4a)')ch10,' rdkss: ERROR -', ch10, errmess(1:min(475, len(errmess)))
      call wrtout(std_out, msg, 'COLL')
     end if
!    END LOW LEVEL
     if (DEBUG) write(*,*)' occ_vec    ik = ',ik,occ_vec
     if (DEBUG) write(*,*)' eigen [eV] ik = ',ik,eigen*Ha_eV
     if (DEBUG) write(*,*)' cg         ik = ',ik,cg
    end if 
    if (master_casts) then 
     call xcast_mpi(eigen,master,spaceComm,ierr)
     call xcast_mpi(cg,master,spaceComm,ierr)
    end if
    en(ik,1:nbndsA,isppol)=eigen(1:nbndsA)
!   Convert double precision real to complex quantities i.e cg ==> wfgd
    do ib=1,nbndsA
     ig=npwkss*(ib-1)
     wfgd(1:npwwfn)=CMPLX(cg(1,ig+1:ig+npwwfn),cg(2,ig+1:ig+npwwfn))
     if (ib>=my_minb .and. ib<=my_maxb) wf(1:npwwfn,ib,ik,isppol)=wfgd(1:npwwfn)
!    In the new algorithm each proc stores the valence states in a separate array
!    FIXME should optimize the communication should check the localrdwf==0 case
     if (read_val .and. ib<=nbvw) wf_val(1:npwwfn,ib,ik,isppol)=wfgd(1:npwwfn)
    end do 
    nullify(Main_folder%coefficients_of_wavefunctions%data2D)
    nullify(Electrons_folder%occupations%data1D)
    nullify(Electrons_folder%eigenvalues%data1D)
    deallocate(cg,eigen,occ_vec)
#endif 
   end if 

   ibg=ibg+nspinor*nbndsA ! index for Cprj_ibz 
  end do !ik
 end do !is
!Reading completed 
!
!=== Test on the orthonormalization of wavefunctions ===
 ibg=0
 do isppol=1,nsppol
  einf1=greatest_real ; esup1=zero
  cinf1=greatest_real ; csup1=zero
  do ik=1,nkibz
!  
!  * Test on the normalization
   ibsp = my_minb-1
   do ib=my_minb,my_maxb
    cdum=czero
    do ispinor=1,nspinor
     ibsp=ibsp+1
     spinor_shift=(ispinor-1)*npwwfn
     wfg1 => wf(1+spinor_shift:npwwfn+spinor_shift,ib,ik,isppol)
     cdum= cdum + overlap_cmplx(wfg1,wfg1,usepaw,Cprj_ibz(:,ibg+ibsp),Cprj_ibz(:,ibg+ibsp),Hdr%typat,Pawtab) 
    end do
    if (REAL(cdum)<einf1) einf1=REAL(cdum)
    if (REAL(cdum)>esup1) esup1=REAL(cdum)
   end do
   if (nprocs>1) then
    call xmin_mpi(einf1,einf,spaceComm,ierr)
    call xmax_mpi(esup1,esup,spaceComm,ierr)
   else
    einf=einf1 ; esup=esup1
   end if
!  === Test on the orthogonality of wavefunctions ===
   do ib=my_minb,my_maxb
    pad1=(ib-1)*nspinor
    do ibp=ib+1,my_maxb
     pad2=(ibp-1)*nspinor
     cdum=czero
     do ispinor=1,nspinor
      ibsp1=pad1+ispinor
      ibsp2=pad2+ispinor
      spinor_shift=(ispinor-1)*npwwfn
      wfg1 => wf(1+spinor_shift:npwwfn+spinor_shift,ib, ik,isppol)
      wfg2 => wf(1+spinor_shift:npwwfn+spinor_shift,ibp,ik,isppol)
      cdum = cdum + overlap_cmplx(wfg1,wfg2,usepaw,Cprj_ibz(:,ibg+ibsp1),Cprj_ibz(:,ibg+ibsp2),Hdr%typat,Pawtab) 
     end do
     if (ABS(cdum)<cinf1) cinf1=ABS(cdum)
     if (ABS(cdum)>csup1) csup1=ABS(cdum)
    end do !ibp
   end do !ib
   if (nprocs>1) then
    call xmin_mpi(cinf1,cinf,spaceComm,ierr)
    call xmax_mpi(csup1,csup,spaceComm,ierr)
   else
    csup=csup1 ; cinf=cinf1
   end if
   ibg=ibg+nspinor*nbndsA ! index for Cprj_ibz 
  end do ! ik
! 
! === Output results for this spin ===
  write(msg,'(2a)')ch10,&
&  ' test on the normalization of the wavefunctions'
  if (nsppol==2) write(msg,'(3a)')ch10,&
&  ' test on the normalization of the wavefunctions with spin ',tag_spin(isppol)
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,f9.6,a,a,f9.6)')&
&  ' min sum_G |a(n,k,G)| = ',einf,ch10,&
&  ' max sum_G |a(n,k,G)| = ',esup
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a)')' test on the orthogonalization of the wavefunctions'
  if (nsppol==2) write(msg,'(2a)')&
&  ' test on the orthogonalization of the wavefunctions with spin ',tag_spin(isppol)
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
  write(msg,'(a,f9.6,a,a,f9.6,a)')&
&  ' min sum_G a(n,k,G)* a(n",k,G) = ',cinf,ch10,&
&  ' max sum_G a(n,k,G)* a(n",k,G) = ',csup,ch10
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
 end do ! isppol
!
!In case of parallelization over bands we are evaluating the scalar product only 
!between bands belonging to the same sub-block, write a COMMENT on ab_out !
 if (MPI_enreg%gwpara==2) then 
  write(msg,'(6a)')&
&  ' rdkss : COMMENT -',ch10,&
&  '  Note that the test on the orthogonalization is not complete ',ch10,&
&  '  since bands are spread among different processors',ch10
  call wrtout(std_out,msg,'COLL') 
  call wrtout(ab_out,msg,'COLL')
 end if
!
!=== Close files ===
 if (i_read.and.accesswff==0) close(untkss)
#if defined HAVE_ETSF_IO
 if (i_read.and.accesswff==3) then 
  call etsf_io_low_close(untkss,lstat,Error_data)
 end if 
#endif
!
!=== Free memory ===
 if (allocated(symrel2)) deallocate(symrel2)
 if (allocated(symrec2)) deallocate(symrec2)
 if (allocated(shlim  )) deallocate(shlim  )
 if (allocated(vkbdb  )) deallocate(vkbdb  )
 if (allocated(vkbdb  )) deallocate(vkbdb  )
 if (allocated(vkbdd  )) deallocate(vkbdd  )
 if (allocated(energyd)) deallocate(energyd)
 if (allocated(wfgd   )) deallocate(wfgd   )

 call hdr_clean(Hdr)

#if defined DEBUG_MODE
 call leave_test(MPI_enreg) 
 write(msg,'(a)')' rdkss : exit '
 call wrtout(std_out,msg,'PERS')
 call flush_unit(std_out)
#endif

 CONTAINS  !===========================================================

  subroutine not_implemented_rdkss(string)



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

  implicit none

! Arguments ------------------------------------
  character(len=*),intent(in) :: string

! Local variables-------------------------------
  character(len=500) :: msg
! *********************************************************************

  write(msg,'(a)')TRIM(string)//', not yet implemented'
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')

 end subroutine not_implemented_rdkss

end subroutine rdkss
!!***
