!{\src2tex{textfont=tt}}
!!****p* ABINIT/newsp
!! NAME
!! newsp
!!
!! FUNCTION
!! This program is used to take wavefunction data at a set of k-points, a
!! lattice constant, and an energy cutoff and write them out at a
!! new set of k-points, lattice constant, and energy cutoff.  Any
!! or all of these parameters is allowed to vary.  Even if no
!! parameters vary, this program may be used to re-order a data
!! set (i.e. change the order of ikpt and kptns), or to change the
!! wavefunction storage mode (istwfk).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, ZL)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!! ABINIT coding rules are NOT followed strictly in this program.
!! Most of NEWSP functionalities are now included in ABINIT (from version 1.9).
!!
!! PARENTS
!!
!! CHILDREN
!!      bstruct_clean,bstruct_init,getmpw,getng,handle_ncerr,hdr_clean,hdr_init
!!      hdr_io,hdr_io_netcdf,hdr_update,herald,importcml,ini_wf_netcdf,instrng
!!      intagm,inupper,invars1,invars2,isfile,kpgsph,leave_new,listkk,matr3inv
!!      metric,newkpt,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program newsp

 use defs_basis
 use defs_datatypes
 use defs_infos
#if defined HAVE_NETCDF
 use netcdf
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_12parser
 use interfaces_13io_mpi
 use interfaces_13ionetcdf
 use interfaces_13iovars
 use interfaces_13recipspace
 use interfaces_13xml
 use interfaces_14iowfdenpot
 use interfaces_14wfs
 use interfaces_15common
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
!
! Set maximum array dimensions **********************************************
!mpsang=highest angular momentum allowed for consideration plus one
!msppol=maximum number of spin-polarization (default=2)
!msym=maximum number of symmetry operations in space group
 integer,parameter :: mpsang=4,msppol=2,msym=48
!Fortran units for disk files
 integer,parameter :: iout=6
 integer :: wfinp=10,wfout=11
 integer :: fform=2
 integer :: bant0,bantot,bantot2,ceksp2,debug,doorth,fform1,fill,formeig,gscase
 integer :: headform,iband,ii,ikpt,ikpt1,ikpt2,ipsp,ireadwf,isppol,isym,itypat
 integer :: jdtset,lenstr,marr,mband1,mband2,mband_rd,mband_upper,mcg2,me_fft
 integer :: mgfft2
 integer :: mkmem,mkmem2,mpw,mpw2,mpw2_tmp,mu,natom,natom2,nbks1,nbks2
 integer :: nfft2
 integer :: nkpt,nkpt2,nproc_fft,npsp,npsp2,nspden,nspinor,nsppol,nsppol2,nsym
 integer :: nsym2
 integer :: ntypat,ntypat2,nu,option,optorth,paral_fft,prtvol,rdwr,restart,sppoldbl
 integer :: tread,unkg2,usepaw
 integer :: ncid_hdr_in,ncid_hdr_out,ncerr,ndim,nvar,natt,uid
 integer :: bravais2(11),identity(3,3,1),intarr(1),ngfft(18),ngfft2(18)
 integer,allocatable :: indkk(:,:),istwfk(:),istwfk2(:),kg2(:,:),nband(:)
 integer,allocatable :: nband2(:),npwarr(:),npwarr2(:),symafm(:),symafmtmp(:)
 integer,allocatable :: symrel(:,:,:),symreltmp(:,:,:)
 real(dp) :: boxcutmin2,dksqmax,ecut,ecut2,ucvol2,zion_max2
 real(dp) :: acell(3),apwsph(3),dprarr(1),gmet(3,3),gmet2(3,3),gprim2(3,3)
 real(dp) :: gprimd(3,3),gprimd2(3,3),rmet2(3,3),rprimd(3,3),rprimd2(3,3)
 real(dp),allocatable :: cg2(:,:),doccde(:),eigen(:),kptns(:,:),occ(:),occ2(:)
 real(dp),allocatable :: tnonstmp(:,:),xred(:,:),zionpsp2(:)
 character(len=24) :: codename
 character(len=fnlen) :: f1,f2,file2
 character(len=strlen) :: string,string_raw
 character(len=30) :: token
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg
 type(dataset_type) :: dtset2
 type(hdr_type) :: hdr,hdr2
 type(bandstructure_type) :: bstruct
 type(pseudopotential_type) :: psps
 type(wffile_type) :: wffinp,wffout
 type(pawrhoij_type), allocatable :: pawrhoij(:)
 type(pawtab_type), allocatable :: pawtab(:,:)
 type(pspheader_type),allocatable :: pspheads(:)

! ************************************************************************
!BEGIN EXECUTABLE SECTION

 codename='NEWSP '//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

!Set up for debugging and set up for orthonormalized output
 debug=1
 doorth=1

!************************************************************************
!Take care of wf input file

 write(06, '(a)' ) ' Enter name of unformatted wf data input file'
 read (05, '(a)' ) f1
 write(06, '(a,a)' ) ' wfinp=',trim(f1)

!DEBUG
!write(6,*)' newsp : enter '
!stop
!ENDDEBUG

!Check that old wf file exists
 call isfile(f1,'old')

 wffinp%accesswff = 0
 wffout%accesswff = 0

#if defined HAVE_NETCDF
!test if file is a netcdf file
 ncerr = nf90_open(path=f1, mode=NF90_NOWRITE, ncid=ncid_hdr_in)
 if (ncerr == NF90_NOERR) then
! NOTE: the choice could be made here to use any combination of
! netcdf or binary input or output.
  wffinp%accesswff = 2
  wffout%accesswff = 2
! here we modify wfinp to save the netCDF identifier in it
  write (*,*) ' newsp : open a netCDF file ', trim(f1), ncid_hdr_in
  wfinp = ncid_hdr_in
  ncerr = nf90_Inquire(ncid=ncid_hdr_in,nDimensions=ndim,nVariables=nvar,&
&  nAttributes=natt,unlimitedDimId=uid)
  call handle_ncerr(ncerr, " general Inquire ")
  write (*,*) 'newsp : input found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
 end if
#endif

 wffinp%unwff=wfinp

!Open wf file for input wf
 if (wffinp%accesswff /= 2) then
  open (unit=wfinp,file=f1,form='unformatted',status='old')
 end if

!Read the header from the file
 rdwr=1
 if (wffinp%accesswff/=2) then
  call hdr_io(fform1,hdr,rdwr,wfinp)
#if defined HAVE_NETCDF
 else if (wffinp%accesswff==2) then
  call hdr_io_netcdf(fform1,hdr,rdwr,wfinp)
 else if (wffinp%accesswff == 3) then
  write (std_out,*) "FIXME: ETSF I/O support in newsp"
#endif
 end if

!Echo the header of the file
 rdwr=4
 call hdr_io(fform1,hdr,rdwr,6)

!Format 1 or 2 are allowed for input wavefunction
 if ( (fform1+1)/2 /= (fform+1)/2 ) then
  write(message, '(a,a,a,a,i10,a,i10,a,a)' ) ch10,&
&  ' newsp: ERROR -',ch10,&
&  '  Input wf file fform=',fform1,' should equal coded fform=',fform,ch10,&
&  '  Action : check that wf file is correct.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 headform=hdr%headform
 bantot=hdr%bantot
 natom=hdr%natom
 ngfft(1:3)=hdr%ngfft(1:3)
 nkpt=hdr%nkpt
 nspden=hdr%nspden
 nspinor=hdr%nspinor
 nsppol=hdr%nsppol
 nsym=hdr%nsym
 ntypat=hdr%ntypat
 npsp=ntypat
 acell(:)=one
 ecut=hdr%ecut_eff
 rprimd(:,:)=hdr%rprimd(:,:)

!Make sure nsppol does not exceed msppol
 if (nsppol>msppol) then
  write(message, '(a,a,a,a,i4,a,a,a)' ) ch10,&
&  ' newsp: ERROR -',ch10,&
&  '  Wf file nsppol=',nsppol,' exceeds msppol=2.',ch10,&
&  '  Action : the wf file must be corrupted, change it.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Make sure nsym does not exceed msym
 if (nsym>msym) then
  write(message, '(a,a,a,a,i8,a,i8,a,a)' ) ch10,&
&  ' newsp: BUG -',ch10,&
&  '  Wf file nsym=',nsym,' exceeds dimensioned msym=',msym,ch10,&
&  '  msym must be raised in newsp, and the code recompiled.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 allocate ( istwfk(nkpt),kptns(3,nkpt),nband(nkpt*msppol))
 allocate ( npwarr(nkpt),occ(bantot))
 allocate ( symafm(msym),symrel(3,3,msym))
 allocate ( xred(3,natom))

 nband(1:nkpt*nsppol)=hdr%nband(1:nkpt*nsppol)
 npwarr(1:nkpt)=hdr%npwarr(1:nkpt)
 symrel(:,:,1:nsym)=hdr%symrel(:,:,1:nsym)
 istwfk(1:nkpt)=hdr%istwfk(1:nkpt)
 kptns(:,1:nkpt)=hdr%kptns(:,1:nkpt)
 occ(1:bantot)=hdr%occ(1:bantot)
 xred(:,:)=hdr%xred(:,:)

!Get mpw, as the maximum value of npwarr(:)
 mpw=maxval(npwarr(:))

!DEBUG
!write(6,*)' newsp : after read, stop'
!stop
!ENDDEBUG

!************************************************************************
!Take care of input file for new data set

 write(06, '(a)' ) ' Enter name of formatted input file for new data set'
 read (05, '(a)' ) file2
 write(06, '(a,a)' ) ' input=',trim(file2)
!Check that second input file exists
 call isfile(file2,'old')

!Read the data from input file,
!Really need: ceksp2, ecut2, istwfk2
!nband2, ngfft2, ntypat2, occ2,
!Note : all variables in the next three calls and the allocate statement
!have 2 as index, EXCEPT : abinit_version,iout,lenstr,string

!strlen from defs_basis module
 option=1
 call instrng (file2,lenstr,option,strlen,string)

!Copy original file, without change of case
 string_raw=string

!To make case-insensitive, map characters of string to upper case:
 call inupper(string(1:lenstr))

!Might import data from CML file(s)
 call importcml(lenstr,string_raw,string,strlen)

 ntypat2=1 ; marr=1
 token = 'ntypat'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) ntypat2=intarr(1)
!Check that ntypat2 is greater than 0
 if (ntypat2<=0) then
  write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&  ' newsp : ERROR -',ch10,&
&  '  Input ntypat must be > 0, but was ',ntypat2,ch10,&
&  '  This is not allowed.  ',ch10,&
&  '  Action : modify ntypat in the input file.'
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

!Note that the default npsp is ntypat
 npsp2=ntypat2 ; marr=1
 token = 'npsp'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 if(tread==1) npsp2=intarr(1)
!Check that npsp is greater than 0
 if (npsp2<=0) then
  write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
&  ' abinit : ERROR -',ch10,&
&  '  Input npsp must be > 0, but was ',npsp2,ch10,&
&  '  This is not allowed.  ',ch10,&
&  '  Action : modify npsp in the input file.'
  call wrtout(06,  message,'COLL')
  call leave_new('COLL')
 end if

 allocate(zionpsp2(npsp2))
 zionpsp2(1:npsp2)=zero
 zion_max2=zero

!In newsp, msym is a parameter.
 allocate(symafmtmp(msym),symreltmp(3,3,msym),tnonstmp(3,msym))

 token = 'natom'
 call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
!Might initialize natom from CML file
 if(tread==0)then
  token = '_natom'
  call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
 end if
 if(tread==1)then
  natom2=intarr(1)
 else
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' abinit: ERROR -',ch10,&
&  '  Input natom must be defined, but was absent.',ch10,&
&  '  Action : check the input file.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 dtset2%natom=natom2
 dtset2%npsp=npsp2
 dtset2%ntypat=ntypat2

 allocate(dtset2%amu(ntypat2),dtset2%densty(ntypat2,4),dtset2%typat(natom2))
 allocate(dtset2%iatfix(3,natom2),dtset2%so_psp(npsp2))
 allocate(dtset2%spinat(3,natom2),dtset2%shiftk(3,8))
 allocate(dtset2%kberry(3,20),dtset2%vel_orig(3,natom2),dtset2%xred_orig(3,natom2))
 allocate(dtset2%ziontypat(ntypat2),dtset2%znucl(npsp2))
 allocate(dtset2%algalch(ntypat2),dtset2%mixalch(npsp2,ntypat2))
 allocate(dtset2%ptcharge(ntypat2),dtset2%quadmom(ntypat2))
 allocate(dtset2%lpawu(ntypat2))
 allocate(dtset2%lexexch(ntypat2))

!Here set defaults
 dtset2%acell_orig(:)=zero
 dtset2%densty(:,:)=zero
 dtset2%iatfix(:,:)=0
 dtset2%spinat(:,:)=zero
 dtset2%typat(:)=0
 dtset2%nsym=0
 dtset2%so_psp(:)=1
 dtset2%rprim_orig(:,:)=zero
 dtset2%rprim_orig(1,1)=one
 dtset2%rprim_orig(2,2)=one
 dtset2%rprim_orig(3,3)=one
 symafmtmp(:)=1
 symreltmp(:,:,:)=0
 do ii=1,3
  symreltmp(ii,ii,:)=1
 end do
 tnonstmp(:,:)=0.0d0
 dtset2%vel_orig(:,:)=zero
 dtset2%xred_orig(:,:)=zero

!Note that this data is transferred directly from header of disk file
!It requires npsp=npsp2
 dtset2%znucl(:)=hdr%znuclpsp(:)

 jdtset=0
!Defaults
 dtset2%mkmem=-1
 dtset2%mkqmem=-1
 dtset2%mk1mem=-1
 dtset2%natpawu=0
 dtset2%natvshift=0
 dtset2%ncenter=0
 dtset2%nconeq=0
 dtset2%nkpt=1
 dtset2%norb=0
 dtset2%nspinor=1
 dtset2%pawspnorb=0
 dtset2%nsppol=1
 dtset2%nkptgw=0
 dtset2%kptopt=0
 dtset2%nshiftk=1

 mpi_enreg%paral_compil_kpt=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%nproc=0
 mpi_enreg%me=0
 mpi_enreg%mode_para="n"
 call invars1(bravais2,dtset2,iout,jdtset,lenstr,&
& mband_upper,mpi_enreg,msym,string,&
& symafmtmp,symreltmp,tnonstmp,zion_max2)

 nsppol2=dtset2%nsppol
 nsym2=dtset2%nsym
 nkpt2=dtset2%nkpt
 allocate(istwfk2(nkpt2),nband2(nkpt2*nsppol2) )

 allocate(dtset2%istwfk(nkpt2),dtset2%occ_orig(mband_upper*nkpt2*nsppol2))
 allocate(dtset2%kpt(3,nkpt2),dtset2%kptns(3,nkpt2),dtset2%wtk(nkpt2))
 allocate(dtset2%nband(nkpt2*nsppol2),dtset2%bdgw(2,dtset2%nkptgw))
 allocate(dtset2%kptgw(3,dtset2%nkptgw))

 allocate(dtset2%symafm(nsym2),dtset2%symrel(3,3,nsym2),dtset2%tnons(3,nsym2))
 dtset2%symafm(:)=symafmtmp(1:nsym2)
 dtset2%symrel(:,:,:)=symreltmp(:,:,1:nsym2)
 dtset2%tnons(:,:)=tnonstmp(:,1:nsym2)

!Really need: ceksp2, ecut2
!ngfft2, ntypat2, occ2

!Here set defaults
 ceksp2=0
 dtset2%kptrlen=zero
 ecut2=-1.0d0
 dtset2%istwfk(:)=0
 dtset2%kpt(:,:)=zero
 dtset2%kptnrm=one
 dtset2%nband(:)=0
 dtset2%nspden=1
 dtset2%nspinor=1
 dtset2%occ_orig(:)=zero
 dtset2%wtk(:)=zero

!These default values should come from indefo.f
!except for iscf, that is set to 0 here

 dtset2%ceksph=0
 dtset2%dilatmx=one
 dtset2%enunit=0
 dtset2%exchn2n3d=0
 dtset2%ionmov=0
 dtset2%intxc=0
 dtset2%iprcch=2
 dtset2%iprcel=0
 dtset2%iprcfc=0
 dtset2%irdwfk=0
 dtset2%iscf=0
 dtset2%isecur=0
 dtset2%ixc=1
 dtset2%nqpt=0
 dtset2%restartxf=0
 dtset2%optcell=0
 dtset2%irdwfq=0
 dtset2%ird1wf=0
 dtset2%irdddk=0
 dtset2%kptopt=0
 dtset2%chkexit=2
 dtset2%ikhxc=0
 dtset2%nbdbuf=0
 dtset2%localrdwf=1
 dtset2%nberry=1
 dtset2%bdberry(1)=0
 dtset2%bdberry(2)=0
 dtset2%bdberry(3)=0
 dtset2%bdberry(4)=0
 dtset2%delayperm=0
 dtset2%signperm=1
 dtset2%nbandkss=0
 dtset2%npwkss=-1
 dtset2%ntypalch=0
 dtset2%berryopt=0
 dtset2%wfoptalg=0
 dtset2%nbdblock=1
 dtset2%kssform=0
 dtset2%usedmatpu=0
 dtset2%usepawu=0
 dtset2%useylm=0
 dtset2%td_mexcit=0
 dtset2%prtdensph=0

 dtset2%rfasr=0
 dtset2%rfatpol(1)=1
 dtset2%rfatpol(2)=1
 dtset2%rfdir(1)=0
 dtset2%rfdir(2)=0
 dtset2%rfdir(3)=0
 dtset2%rfelfd=0
 dtset2%rfmeth=1
 dtset2%rfmgfd=0
 dtset2%rfphon=0
 dtset2%rfstrs=0
 dtset2%rfthrd=0
 dtset2%rfuser=0
 dtset2%rf1elfd=0
 dtset2%rf1phon=0
 dtset2%rf2elfd=0
 dtset2%rf2phon=0
 dtset2%rf3elfd=0
 dtset2%rf3phon=0

 dtset2%qpt(:)=zero
 dtset2%qptnrm=one

!DEBUG
!write(6,*)' newsp : string '
!write(6,*)string(1:lenstr)
!ENDDEBUG
!DEBUG
!write(6,*)' newsp : before invars2'
!write(6,*)' nkpt2=',nkpt2
!stop
!write(6,*)' newsp : before invars2 '
!write(6,*)' newsp : dtset2%nband=',dtset2%nband
!ENDDEBUG

 usepaw=0
!Default values for FFT sequential execution
 paral_fft=0
 me_fft=0
 nproc_fft=1

 allocate(pspheads(npsp2))
 call invars2(bravais2,dtset2,iout,jdtset,lenstr,mband_upper,msym,&
& npsp2,pspheads,string,usepaw,zionpsp2)
 deallocate(pspheads)

 istwfk2(:)=dtset2%istwfk(:)
 ecut2=dtset2%ecut
 boxcutmin2=dtset2%boxcutmin
 ngfft2(:)=dtset2%ngfft(:)
 ceksp2=dtset2%ceksph

!DEBUG
!write(6,*)' newsp : after invars2 '
!write(6,*)' newsp : nkpt2=',nkpt2
!write(6,*)' newsp : dtset2%nband=',dtset2%nband
!ENDDEBUG

 deallocate(symafmtmp,symreltmp,tnonstmp)

!Compute mgfft2,mpw2_tmp,nfft2 for this data set
 rprimd2(:,:)=dtset2%rprimd_orig(:,:)
 call metric(gmet2,gprimd2,-1,rmet2,rprimd2,ucvol2)
 call getng(boxcutmin2,ecut2,gmet2,me_fft,mgfft2,nfft2,ngfft2,nproc_fft,nsym2,1,paral_fft,dtset2%symrel)
 dtset2%ngfft(:)=ngfft2(:)
 mpi_enreg%nproc_fft=nproc_fft; ngfft(10)=nproc_fft
 mpi_enreg%me_fft=me_fft;ngfft(11)=me_fft
 mpi_enreg%fft_option_lob=1
 call getmpw(ecut2,dtset2%exchn2n3d,gmet2,istwfk2,dtset2%kptns,&
& mpi_enreg,mpw2_tmp,nkpt2,ucvol2)
 nband2(1:nkpt2*nsppol2)=dtset2%nband(1:nkpt2*nsppol2)
 mband2=maxval(nband2(1:nkpt2*nsppol2))
 dtset2%mband=mband2

 allocate(occ2(mband2*nkpt2*nsppol2))
 occ2(1:mband2*nkpt2*nsppol2)=dtset2%occ_orig(1:mband2*nkpt2*nsppol2)

!NOTE: At this point, much more checking could be done to
!guarantee some consistency between the new input file and the
!old wf file (e.g. same number of atoms, etc.).  I am NOT doing
!this checking right now (23 Jun 1993).  Most of the data from
!the original wf file (e.g. natom) is being passed along to the newly
!created file.

!Check same choice of spin-polarization
!(But I think it should work ...)
 if (nsppol/=nsppol2) then
  write(message, '(a,a,a,a,i6,a,i6,a,a,a,a)' ) ch10,&
&  ' newsp: ERROR -',ch10,&
&  '  Wf disk file nsppol=',nsppol,&
&  '  not equal input file nsppol=',nsppol2,ch10,&
&  '  These must be same for newsp to run.',ch10,&
&  '  Action : correct one of these files.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!************************************************************************
!At last, need name of wf output file, then initialize its header.

 write(06, '(a)' ) ' Enter name of unformatted wf data output file'
 read(05, '(a)' ) f2
 write(06, '(a,a)' ) ' wfout=',trim(f2)
!Check that new file does NOT already exist
 call isfile(f2,'new')

 write(06, '(a,3i6)' ) ' wf    file ng',ngfft(1:3)
 write(06, '(a,3i6)' ) ' input file ng',ngfft2(1:3)

 write(06, '(a,1p,e12.4)' ) ' wf    file ecut',ecut
 write(06, '(a,1p,e12.4)' ) ' input file ecut',ecut2

 write(06, '(a,1p,3e15.7)' ) ' wf    file acell',acell(1:3)
 write(06, '(a,1p,3e15.7)' ) ' input file acell',dtset2%acell_orig(1:3)

!Compute reciprocal space metric gmet from original unit cell
 call matr3inv(rprimd,gprimd)
 do nu=1,3
  do mu=1,3
   gmet(mu,nu)=gprimd(1,mu)*gprimd(1,nu)+&
&   gprimd(2,mu)*gprimd(2,nu)+&
&   gprimd(3,mu)*gprimd(3,nu)
  end do
 end do

 sppoldbl=1
 if(minval(dtset2%symafm(:))==-1)then
  if(nsppol==1 .and. nsppol2==2)sppoldbl=2
 end if
 allocate (indkk(nkpt2*sppoldbl,6))

!Find nearest old k point to each new k point
 call listkk(dksqmax,gmet,indkk,kptns,dtset2%kptns,nkpt,nkpt2,nsym2,sppoldbl,&
& dtset2%symafm,dtset2%symrel,1)

 if (debug>0 .and. sppoldbl==1) then
! For new wf file, compute overall total number of bands, bantot
  write(6,*)' original bantot=',bantot
  bantot2=0
  do isppol=1,nsppol2
   do ikpt2=1,nkpt2
    bantot2=bantot2+nband2(indkk(ikpt2,1)+(isppol-1)*nkpt2)
   end do
  end do
  write(6,*)' final    bantot=',bantot2
 end if

!DEBUG
!write(6,*)' newsp : before kpgsph, nkpt2=',nkpt2
!stop
!ENDDEBUG

!Compute npw for each k point of new wf file, then compute actual mpw2
 allocate(kg2(3,mpw2_tmp),npwarr2(nkpt2))
 do ikpt2=1,nkpt2
  call kpgsph(ecut2,dtset2%exchn2n3d,gmet2,0,ikpt2,istwfk2(ikpt2),kg2,dtset2%kptns(:,ikpt2),&
&  1,mpi_enreg,mpw2_tmp,npwarr2(ikpt2))
 end do
 mpw2=maxval(npwarr2(:))
 deallocate(kg2)
 allocate(kg2(3,mpw2))

!DEBUG
!write(6,*)' newsp : after kpgsph '
!write(6,*)' nsppol,nkpt,nkpt2',nsppol,nkpt,nkpt2
!write(6,*)' nband(1:nkpt)',nband(1:nkpt)
!write(6,*)' nband2(1:nkpt2)',nband2(1:nkpt2)
!write(6,*)' occ(1:nband(1)*nkpt)',occ(1:nband(1)*nkpt)
!stop
!ENDDEBUG

!Copy over occupancies, band by band, for new k points
!(also gives value of bantot2)
!--NOTE that OLD occupancies are used here, NOT the presumably
!desired new occupancies as requested in the new input file

 bantot2=0
!Warning : should be adapted to sppoldbl=2
 do isppol=1,nsppol
  do ikpt2=1,nkpt2
!  Count over skipped bands from original set
   bant0=0
   if (indkk(ikpt2,1)>1) then
    do ikpt1=1,indkk(ikpt2,1)-1
     bant0=bant0+nband(ikpt1+(isppol-1)*nkpt)
    end do
   end if
!  Determine nband2 for each k, spin--this may not exceed
!  nband1 for the same spin and the associated k
!  nbks1 is the nband for the same spin and associated k
!  nbks2 is the requested nband from formatted input file
   nbks1=nband(indkk(ikpt2,1)+(isppol-1)*nkpt)
   nbks2=nband2(ikpt2+(isppol-1)*nkpt2)
!  If number of bands is being increased, print warning and
!  reset nband2 to only available number
!  WARNING : does not take into account nspinor, while it should
   if (nbks2>nbks1) then
    write(06, '(a,i2,i5,2i6,a,i6)' ) &
&    ' newsp: isppol,ikpt2,nband1,nband2=',&
&    isppol,ikpt2,nbks1,nbks2,'=> lower nband2 to',nbks1
    nband2(ikpt2+(isppol-1)*nkpt2)=nbks1
    nbks2=nbks1
   end if
!  Define occupancies for new wf file using appropriate
!  nband for each k,spin
   do iband=1,nbks2
    occ2(iband+bantot2)=occ(iband+bant0)
   end do
   bantot2=bantot2+nbks2
  end do
 end do

!After all checks out, write first part of new header for new wf file

!The following are from the original wf file, NOT from the new
!input file:
!abinit_version, fform
!natom,nsppol,ntypat
!occ(shifted)
!All data from psp headers (of course), and xred.
!The following are modified from original wf file:
!bantot (now reflects new set of k points),ngfft,nkpt,
!nband(shifted and possibly made smaller), ecut,rprimd,kptns,
!Note especially that no account is taken of any consistency in
!the use of symmetry.

!Set up the psps part of the header (many other infos form
!the psps datastructure, but only these matters here)
 psps%ntypat=ntypat
 psps%npsp=npsp
 allocate(psps%pspcod(npsp))
 allocate(psps%pspdat(npsp))
 allocate(psps%pspso(npsp))
 allocate(psps%pspxc(npsp))
 allocate(psps%title(npsp))
 allocate(psps%znuclpsp(npsp))
 allocate(psps%znucltypat(ntypat))
 allocate(psps%zionpsp(npsp))
 do ipsp=1,npsp
  psps%pspcod(ipsp)=hdr%pspcod(ipsp)
  psps%pspdat(ipsp)=hdr%pspdat(ipsp)
  psps%pspso(ipsp)=hdr%pspso(ipsp)
  psps%pspxc(ipsp)=hdr%pspxc(ipsp)
  psps%title(ipsp)=hdr%title(ipsp)
  psps%znuclpsp(ipsp)=hdr%znuclpsp(ipsp)
  psps%zionpsp(ipsp)=hdr%zionpsp(ipsp)
 end do
 do itypat=1,ntypat
  psps%znucltypat(itypat)=hdr%znucltypat(itypat)
 end do

!Initialize band structure datatype
 allocate(doccde(bantot2),eigen(bantot2))
 doccde(:)=zero ; eigen(:)=zero
 call bstruct_init(bantot2,bstruct,doccde,eigen,dtset2%istwfk,dtset2%kptns,&
& nband2,nkpt2,npwarr2,nsppol2,occ2,dtset2%wtk)
 deallocate(doccde,eigen)

!DEBUG
!write(6,*)' newsp: before hdr_init'
!ENDDEBUG

!Set up the band structure information in hdr2
 gscase=0
 call hdr_init(bstruct,abinit_version,dtset2,hdr2,pawtab,gscase,psps)

!DEBUG
!write(6,*)' newsp: after hdr_init'
!ENDDEBUG


!Clean band structure datatype
 call bstruct_clean(bstruct)

!Update header, with evolving variables
!WARNING : The number of atom in the disk file and the input file
!cannot change.
 allocate(pawrhoij(natom*usepaw))
 call hdr_update(bantot2,hdr%etot,hdr%fermie,hdr2,natom,&
& hdr%residm,rprimd2,occ2,pawrhoij,usepaw,xred)
 deallocate(pawrhoij)

#if defined HAVE_NETCDF
 if(wffout%accesswff==2) then
! Create empty netCDF file
  ncerr = nf90_create(path=f2, cmode=NF90_CLOBBER, ncid=ncid_hdr_out)
  call handle_ncerr(ncerr," create netcdf wavefunction file")
  ncerr = nf90_close(ncid_hdr_out)
  call handle_ncerr(ncerr," close netcdf wavefunction file")
 else if (wffout%accesswff == 3) then
  write (std_out,*) "FIXME: ETSF I/O support in newsp"
 end if
#endif

!Open wf file for output wf
 if (wffout%accesswff/=2) then
  open (unit=wfout,file=f2,form='unformatted',status='unknown')
#if defined HAVE_NETCDF
 else if (wffout%accesswff==2) then
! here we modify wfout to save the netCDF identifier in it
  ncerr = nf90_open(path=f2, mode=NF90_WRITE, ncid=ncid_hdr_out)
  call handle_ncerr(ncerr," newsp : open netcdf wavefunction file")
  write (*,*) ' newsp : open a netCDF file ', trim(f2), ncid_hdr_out
  wfout = ncid_hdr_out
 else if (wffout%accesswff == 3) then
  write (std_out,*) "FIXME: ETSF I/O support in newsp"
#endif
 end if

!Only now can we initialize wffout%unwff, in the NetCDF case.
 wffout%unwff=wfout

 rdwr=2
 if (wffout%accesswff/=2) then
  call hdr_io(fform,hdr2,rdwr,wfout)
#if defined HAVE_NETCDF
 else if (wffout%accesswff==2) then
  call hdr_io_netcdf(fform,hdr2,rdwr,wfout)

  call ini_wf_netcdf(mpw2,ncid_hdr_out,gscase)
  ncerr = nf90_Inquire(ncid=ncid_hdr_out,nDimensions=ndim,nVariables=nvar,&
&  nAttributes=natt,unlimitedDimId=uid)
  call handle_ncerr(ncerr, " general Inquire ")
  write (*,*) 'newsp : output found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
 else if (wffout%accesswff == 3) then
  write (std_out,*) "FIXME: ETSF I/O support in newsp"
#endif
 end if

!Echo header of future disk file
 write(6,*)' Echo the header of the new file, ',trim(f2)

 rdwr=4
 call hdr_io(fform,hdr2,rdwr,6)

 call hdr_clean(hdr2)
 deallocate(psps%pspcod,psps%pspdat,psps%pspso,psps%pspxc)
 deallocate(psps%title,psps%znuclpsp,psps%zionpsp,psps%znucltypat)

!write(unit=wfout) hdr%etot
!Header data is now complete for new file
!For old file, one is ready to read the first wf block

!nband  tells how many bands are to be read (how many expected)
!nband2 tells how many are to be written

!************************************************************************

 if (debug>0) then
  write(6,*)'newsp: Calling newkpt with :'
  write(6,*)' mband2,mpw,mpw2,nkpt,nkpt2='
  write(6, '(6i5)' ) mband2,mpw,mpw2,nkpt,nkpt2
  write(6,*)' nspinor,nspinor2,nsppol,nsppol2='
  write(6, '(6i5)' ) nspinor,dtset2%nspinor,nsppol,nsppol2
  write(6, '(15i4)')istwfk(1:nkpt)
  write(6, '(15i4)')istwfk2(1:nkpt2)
 end if

 fill=0 ; formeig=0 ; ireadwf=1 ; mkmem=0 ; mkmem2=0
 prtvol=dtset2%prtvol ; mpi_enreg%paralbd=0 ; restart=2

 mband1=maxval(nband(1:nkpt*nsppol))
 mband_rd=min(mband1,(mband2/dtset2%nspinor)*nspinor)
 mcg2=max(mpw*nspinor*mband_rd,mpw2*mband2*dtset2%nspinor)
 allocate(cg2(2,mcg2),eigen(mband2*nkpt2*nsppol))

!DEBUG
!write(6,*)' newsp : before newkpt, stop'
!stop
!ENDDEBUG
 optorth=1
 call newkpt (ceksp2,cg2,debug,doorth,ecut,ecut2,ecut2,eigen,dtset2%exchn2n3d,fill,formeig,&
& gmet,gmet2,headform,indkk,iout,ireadwf,istwfk,istwfk2,kg2,kptns,dtset2%kptns,&
& mband2,mcg2,mkmem,mkmem2,mpi_enreg,mpw,mpw2,&
& nband,nband2,ngfft,nkpt,nkpt2,npwarr,npwarr2,nspinor,dtset2%nspinor,&
& nsppol,nsppol2,nsym2,occ2,optorth,prtvol,restart,rprimd2,sppoldbl,&
& dtset2%symafm,dtset2%symrel,dtset2%tnons,unkg2,wffinp,wffout)

!-------------------------------------------------------------------

 deallocate(istwfk,kptns,nband)
 deallocate(npwarr,occ)
 deallocate(symafm,symrel)
 deallocate(xred,zionpsp2)
 deallocate(dtset2%amu,dtset2%bdgw,dtset2%kptgw,dtset2%densty,dtset2%typat,dtset2%istwfk)
 deallocate(dtset2%kpt,dtset2%kptns)
 deallocate(dtset2%iatfix,dtset2%so_psp,dtset2%spinat,dtset2%shiftk,dtset2%kberry)
 deallocate(dtset2%symafm,dtset2%symrel,dtset2%tnons,dtset2%wtk,dtset2%ziontypat)
 deallocate(dtset2%vel_orig,dtset2%xred_orig,dtset2%occ_orig,dtset2%nband,dtset2%znucl)
 deallocate(dtset2%algalch,dtset2%mixalch)
 deallocate(dtset2%lpawu,dtset2%ptcharge,dtset2%quadmom)
 deallocate(dtset2%lexexch)
 deallocate(istwfk2,nband2)
 deallocate(occ2,indkk,npwarr2,kg2,cg2,eigen)

 call hdr_clean(hdr)

#if defined HAVE_NETCDF
 if (wffinp%accesswff==2) then
  ncerr = nf90_close (ncid=ncid_hdr_in)
  call handle_ncerr(ncerr," newsp : final close netcdf input wavefunction file")
 else if (wffinp%accesswff == 3) then
  write (std_out,*) "FIXME: ETSF I/O support in newsp"
 end if
 if (wffout%accesswff==2) then
  ncerr = nf90_close (ncid=ncid_hdr_out)
  call handle_ncerr(ncerr," newsp : final close netcdf output wavefunction file")
 else if (wffout%accesswff == 3) then
  write (std_out,*) "FIXME: ETSF I/O support in newsp"
 end if
#endif

 write(6,*)' newsp:  program ends normally'

 end program newsp
!!***
