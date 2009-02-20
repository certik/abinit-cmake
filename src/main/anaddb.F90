!{\src2tex{textfont=tt}}
!!****p* ABINIT/anaddb
!! NAME
!! anaddb
!!
!! FUNCTION
!! Main routine for analysis of the interatomic force constants and associated
!! properties.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG,DCA,JCC,CL,XW)
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
!! WARNING
!!
!! PARENTS
!!
!! CHILDREN
!!      asria9,asrprs,bigbx9,chneu9,diel9,dtchi,dtech9,elast9,electrooptic
!!      elphon,gtblk9,gtdyn9,herald,init9,inprep8,instr9,instrng,inupper
!!      invars9,isfile,leave_new,mkifc9,mkphdos,mkrdim,outlwf9,outvars9,phfrq3
!!      piezo9,prtph3,ramansus,rdddb9,relaxpol,sortph,symph3,thm9,thmeig,timein
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program anaddb

 use defs_basis
 use defs_datatypes
 use defs_infos


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_12parser
 use interfaces_16response
 use interfaces_17ddb
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
! Set array dimensions
!  msym =maximum number of symmetry elements in space group
!  (not really worth to eliminate msym)
 integer,parameter :: msym=192
!Define input and output unit numbers (some are defined in defs_basis -all should be there ...):
 integer,parameter :: ddbun=2,iodyn=8,mdun=3
!Define unit number for the files that can be analysed with band2eps
 integer,parameter :: udispl=19,ufreq=18
 integer :: choice,dimekb,dims
 integer :: iatom,iblok,iblok_stress,ibloknl,idir,ii,index,index1,iphl1,iphl2,jj,lenstr
 integer :: mband,mblktyp,mpert,msize,natom,nblok
 integer :: nkpt,nph2l,nrpt,nsym,ntypat
 integer :: nunit,nversn,occopt,option,rftyp
 integer :: usepaw,vrsddb
 integer :: elphflag,unitgkk
 integer :: rfelfd(4),rfphon(4),rfstrs(4),symq(4,2,msym)
 integer :: themflag
 integer :: phfrq_unit
 integer :: symrec(3,3,msym),symrel(3,3,msym)
 integer,allocatable :: blkflg(:,:),blktyp(:),d2flg(:),indsym(:,:),typat(:)
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: etotal,tcpu,tcpui,twall
 real(dp) :: twalli,ucvol
 real(dp) :: acell(3),dielt(3,3),dielt_rlx(3,3),elast(6,6),epsinf(3,3),gmet(3,3),gprim(3,3)
 real(dp) :: pel(3)
 real(dp) :: piezo(6,3),qphnrm(3),qphon(3,3),qpt(3),rmet(3,3),rprim(3,3)
 real(dp) :: rprimd(3,3),strten(6),tnons(3,msym)
 real(dp),allocatable :: amu(:),atmfrc(:,:,:,:,:,:),blknrm(:,:),blkqpt(:,:)
 real(dp),allocatable :: blkval(:,:,:),blkval2(:,:,:,:,:),d2asr(:,:,:,:),d2cart(:,:),dchide(:,:,:)
 real(dp),allocatable :: dchidt(:,:,:,:),displ(:),dyewq0(:),eigval(:,:)
 real(dp),allocatable :: eigvec(:,:,:,:,:),fact_oscstr(:,:,:),instrain(:,:),kpnt(:,:,:)
 real(dp),allocatable :: fred(:,:),lst(:),phfreq(:,:),phfrq(:),eigvect(:,:,:,:,:)
 real(dp),allocatable :: rcan(:,:),rpt(:,:),rsus(:,:,:),trans(:,:)
 real(dp),allocatable :: singular(:),uinvers(:,:), vtinvers(:,:)
 real(dp),allocatable :: wghatm(:,:,:),xcart(:),xred(:,:),zeff(:,:,:),zion(:)
 integer :: irpt,irpt_new,nrpt_new
 real(dp),allocatable :: atmfrc_tmp(:,:,:,:,:,:),rpt_tmp(:,:),wghatm_tmp(:,:,:)
 real(dp),allocatable :: minvers(:,:)
 character(len=24) :: codename
 character(len=strlen) :: string
 character(len=fnlen) :: filnam(7),elph_base_name,ddkfilename
 character(len=fnlen) :: phonon_freq_filename,outfilename_radix
 character(len=500) :: message
 type(anaddb_dataset_type) :: anaddb_dtset
 ! here only since needed for call to rwwf
 type(MPI_type) :: mpi_enreg

!******************************************************************
!BEGIN EXECUTABLE SECTION

!Initialisation of the timing
 call timein(tcpui,twalli)

 codename='ANADDB'//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

!Initialise the code : write heading, and read names of files.
 call init9(filnam)

!******************************************************************

 write(message, '(a,a,a,a)' )&
& ch10,ch10,' Read the input file',ch10
 call wrtout(6,message,'COLL')

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a)' )&
& '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
 call wrtout(6,message,'COLL')

!Must read natom from the DDB before being able to allocate
!some arrays needed for invars9
 vrsddb=010929
 call inprep8(dimekb,filnam(3),mband,mblktyp,msym,&
& natom,nblok,nkpt,ntypat,ddbun,usepaw,vrsddb)

 mpert=natom+6
 msize=3*mpert*3*mpert
 if(mblktyp==3)msize=msize*3*mpert

!Read the input file, and store the information in a long string of characters
!strlen from defs_basis module
 option=1
 call instrng (filnam(1),lenstr,option,strlen,string)

!To make case-insensitive, map characters to upper case:
 call inupper(string(1:lenstr))

 write(*,*) 'will read the inputs completely'

!Read the inputs
 call invars9 (anaddb_dtset,lenstr,natom,tmp_unit,qtol,string)

 nph2l=anaddb_dtset%nph2l
 allocate(lst(nph2l))

 write(*,*) 'read the inputs completely'

!Echo the inputs to console
 nunit=6
 call outvars9 (anaddb_dtset,nunit)

!save output file name
 outfilename_radix = filnam(2)

!Open output files iodyn and ab_out (might change its name if needed)
 call isfile(filnam(2),'new')
 open (unit=ab_out,file=filnam(2),form='formatted',status='new')
 rewind (unit=ab_out)
 call herald(codename,abinit_version,ab_out)
 if (anaddb_dtset%eivec==3) then
  call isfile(filnam(4),'new')
  open (unit=iodyn,file=filnam(4),form='formatted',status='new')
  rewind (unit=iodyn)
 end if

!Echo the inputs to long printout
 call outvars9 (anaddb_dtset,ab_out)

!******************************************************************

!Read the DDB information,
!also perform some checks, and symmetrize
!partially the DDB

 write(message, '(a,a)' ) &
& ' read the DDB information and perform some checks',ch10
 call wrtout(6,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 call timein(tcpu,twall)
 write(message, '(a,f11.3,a,f11.3,a,a)' )&
& '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec',ch10
 call wrtout(6,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 allocate(amu(ntypat),d2cart(2,msize))
 allocate(blkflg(msize,nblok),blknrm(3,nblok),blkqpt(9,nblok))
 allocate(blktyp(nblok),blkval(2,msize,nblok))
 allocate(indsym(4,msym*natom))
 allocate(typat(natom),xcart(3*natom),xred(3,natom),zion(ntypat))
 allocate(instrain(3*natom,6))

 call rdddb9(acell,anaddb_dtset%atifc,amu,blkflg,blknrm,blkqpt,&
& blktyp,blkval,ddbun,dimekb,filnam(3),gmet,gprim,indsym,ab_out,&
& mband,mpert,msize,msym,&
& anaddb_dtset%natifc,natom,nblok,nkpt,nsym,ntypat,&
& occopt,rmet,rprim,symq,symrec,symrel,1,&
& tnons,typat,ucvol,usepaw,xcart,xred,zion)

 call mkrdim(acell,rprim,rprimd)


!Now the whole DDB is in central memory, contained in the
!array blkval(2,msize,nblok).
!The information on it is contained in the four arrays
!blkflg(msize,nblok) : blok flag for each element
!blkqpt(9,nblok)  : blok wavevector (unnormalized)
!blknrm(3,nblok)  : blok wavevector normalization
!blktyp(nblok)    : blok type

 allocate(displ(2*3*natom*3*natom),dyewq0(3*3*natom),d2asr(2,3,3,natom))
 allocate(eigval(3,natom),eigvec(2,3,natom,3,natom))
 allocate(phfrq(3*natom),rcan(3,natom),trans(3,natom))
 allocate(zeff(3,3,natom))

!**********************************************************************
!**********************************************************************

!Acoustic Sum Rule

!In case the interatomic forces are not calculated, the
!ASR-correction (d2asr) has to be determined here from
!the Dynamical matrix at Gamma.
 if(anaddb_dtset%ifcflag==0)then

! Find the Gamma block in the DDB (no need for E-field entries)
  qphon(:,1)=0.0d0
  qphnrm(1)=0.0d0
  rfphon(1:2)=1
  rfelfd(:)=0
  rfstrs(:)=0
  rftyp=anaddb_dtset%rfmeth

  call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&  qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

  if(anaddb_dtset%asr==1 .or. anaddb_dtset%asr==2) then
   if (iblok <= nblok) then
    call asria9(anaddb_dtset%asr,1,d2asr,blkval(1,1,iblok),mpert,natom)
   else
    d2asr(:,:,:,:)=0.0d0
   end if
  end if

! Rotational invariance for 1D and 0D systems
  
  if(anaddb_dtset%asr==3 .or. anaddb_dtset%asr==4) then
   dims=3*natom*(3*natom-1)/2
   allocate(uinvers(1:dims,1:dims),vtinvers(1:dims,1:dims))
   allocate(singular(1:dims))
   uinvers=0d0
   vtinvers=0d0
   singular=0d0
   if (iblok <= nblok) then
    call asrprs(anaddb_dtset%asr,1,3,uinvers,vtinvers,singular,blkval(1,1,iblok),mpert,natom,rprim,xcart)
   end if
  end if  
 end if

!**********************************************************************

!Dielectric Tensor and Effective Charges

!Look after the Gamma Blok in the DDB
 qphon(:,1)=0.0d0
 qphnrm(1)=0.0d0
 rfphon(1:2)=1
 rfelfd(1:2)=2
 rfstrs(1:2)=0
 rftyp=anaddb_dtset%rfmeth

 call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
& qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

!Compute effective charges and dielectric tensor only if the
!Gamma-blok was found in the DDB or the occupation is metallic
!In case it was not found, iblok = nblok + 1

 if ((iblok <= nblok).or.(3<=occopt.and.occopt<=7)) then

  write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&  ' Dielectric Tensor and Effective Charges ',ch10
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call timein(tcpu,twall)
  write(message, '(a,f11.3,a,f11.3,a)' )&
&  '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  write(message, '(a,i6)' )' The Gamma block is : ',iblok
  call wrtout(6,message,'COLL')

  if (0<=occopt .and. occopt<=2) then

!  Make the imaginary part of the Gamma block vanish
   write(message, '(a,a,a,a,a)'  ) ch10,&
&   ' anaddb : Zero the imaginary part of the Dynamical Matrix at Gamma,',ch10,&
&   '   and impose the ASR on the effective charges ',ch10
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  Impose the charge neutrality on the effective charges
!  and eventually select some parts of the effective charges
   call chneu9(anaddb_dtset%chneut,blkval(1,1,iblok),mpert,natom,ntypat,&
&   anaddb_dtset%selectz,typat,zion)

!  Extraction of the dielectric tensor and the effective charges
   call dtech9(blkval,dielt,iblok,mpert,natom,nblok,zeff)

  else if (3<=occopt.and.occopt<=7) then

   write(message, '(a,a)' ) ch10,&
&   ' Metallic case : effective charges are set to 0'
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   zeff(:,:,:)=0.0d0
   dielt(:,:)=0.0d0
   dielt(1,1)=1.0d0 ; dielt(2,2)=1.0d0 ; dielt(3,3)=1.0d0

  else

   write(message, '(a,a,a,a,i6,a,a,a,a,a)' )ch10,&
&   ' anaddb : ERROR -',ch10,&
&   '  The input DDB value of occopt is',occopt,',',ch10,&
&   '  while it should be between 0 and 7.',ch10,&
&   '  Action : check the value of occopt in your DDB.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')

  end if

 end if    ! iblok < nblok

!**********************************************************************

!Structural response at fixed polarization
 if (anaddb_dtset%polflag == 1) then

  allocate(d2flg(msize))

  if(iblok<=nblok)then

!  Save the second-order derivatives
   d2cart(1:2,1:msize) = blkval(1:2,1:msize,iblok)
   d2flg(1:msize) = blkflg(1:msize,iblok)

  else ! the gamma blok has not been found

   if(anaddb_dtset%relaxat==0 .and. anaddb_dtset%relaxstr==0)then

!   The gamma blok is not needed
    d2cart(1:2,1:msize)=zero
    d2flg(1:msize)=1

   else ! There is a problem !

    write(message, '(10a)' ) ch10,&
&    ' anaddb : ERROR -',ch10,&
&    '  The dynamical matrix at Gamma is needed, in order to perform ',ch10,&
&    "  relaxation at constant polarisation (Na Sai's method)",ch10,&
&    '  However, this was not found in the DDB.',ch10,&
&    '  Action : complete your DDB with the dynamical matrix at Gamma.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')

   end if

  end if ! iblok <= nblok

! Read the block with the total energy
  qphon(:,:) = zero
  qphnrm(:) = zero
  rfphon(:) = 0
  rfelfd(:) = 0
  rfstrs(:) = 0
  rftyp = 0
  call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&  qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
  etotal = blkval(1,1,iblok)

! Read the block with the gradients
  allocate(fred(3,natom))
  rftyp = 4
  rfelfd(:) = 2
  if (anaddb_dtset%relaxat == 1) rfphon(:) = 1
  if (anaddb_dtset%relaxstr == 1) rfstrs(:) = 3
  call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&  qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

  if (anaddb_dtset%relaxat == 1) then
   index = 0
   do iatom = 1, natom
    do idir = 1, 3
     index = index + 1
     fred(idir,iatom) = blkval(1,index,iblok)
    end do
   end do
  end if

  pel(1:3) = blkval(1,3*natom+4:3*natom+6,iblok)

  if (anaddb_dtset%relaxstr == 1) then
   index = 3*natom + 6
   do ii = 1, 6
    index = index + 1
    strten(ii) = blkval(1,index,iblok)
   end do
  end if

  call relaxpol(d2flg,d2cart,etotal,fred,anaddb_dtset%iatfix,&
&  indsym,ab_out,anaddb_dtset%istrfix,&
&  mpert,msize,msym,anaddb_dtset%natfix,natom,&
&  anaddb_dtset%nstrfix,nsym,ntypat,pel,&
&  anaddb_dtset%relaxat,anaddb_dtset%relaxstr,&
&  rprimd,strten,symrel,anaddb_dtset%targetpol,typat,ucvol,xcart,xred,zion)

  deallocate(fred,d2flg)

 end if

!***************************************************************************

!Compute non-linear optical susceptibilities and
!First-order change in the linear dielectric susceptibility
!induced by an atomic displacement

 if (anaddb_dtset%nlflag > 0) then

  qphon(:,:) = 0_dp
  qphnrm(:)  = 1_dp
  rfphon(1)  = 1 ; rfphon(2:3) = 0
  rfelfd(:)  = 2
  rfstrs(:)  = 0
  rftyp = 3

  call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&  qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

  ibloknl = iblok
  allocate(dchide(3,3,3),dchidt(natom,3,3,3))

  call dtchi(blkval(:,:,ibloknl),dchide,dchidt,mpert,natom,anaddb_dtset%ramansr)

 end if ! nlflag

!**********************************************************************
!**********************************************************************

!Interatomic Forces Calculation

!DEBUG
!write(6,*)' anaddb : before ifcflag check, ifcflg=',anaddb_dtset%ifcflag,anaddb_dtset%thmflag
!stop
!ENDDEBUG

 if (anaddb_dtset%ifcflag/=0 ) then

  write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&  ' Calculation of the interatomic forces ',ch10
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call timein(tcpu,twall)
  write(message, '(a,f11.3,a,f11.3,a)' )&
&  '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! Compute the number of points (cells) in real space
  choice=0
  allocate(rpt(3,1))
  call bigbx9(anaddb_dtset%brav,choice,1,anaddb_dtset%ngqpt,&
&  anaddb_dtset%nqshft,nrpt,rprim,rpt)
  deallocate(rpt)

  allocate(atmfrc_tmp(2,3,natom,3,natom,nrpt),rpt_tmp(3,nrpt))
  allocate(wghatm_tmp(natom,natom,nrpt))

  call mkifc9(acell,amu,anaddb_dtset,atmfrc_tmp,blkflg,blknrm,blkqpt,blktyp,blkval,&
&  dielt,displ,dyewq0,d2cart,eigval,eigvec,gmet,gprim,indsym,ab_out,&
&  mpert,msym,natom,nblok,nrpt,nsym,ntypat,phfrq,rcan,rmet,rprim,&
&  rpt_tmp,symrec,symrel,tcpui,tnons,trans,twalli,typat,&
&  ucvol,wghatm_tmp,xred,zeff)

! Only conserve the necessary points in rpt: in the FT algorithm the
! the order of the points is unimportant
  nrpt_new = 0
  do irpt=1,nrpt
   if (sum(wghatm_tmp(:,:,irpt)) /= 0) then
    nrpt_new = nrpt_new+1
   end if
  end do

  allocate (atmfrc(2,3,natom,3,natom,nrpt_new),rpt(3,nrpt_new))
  allocate(wghatm(natom,natom,nrpt_new))

  irpt_new = 1
  do irpt=1,nrpt
   if (sum(wghatm_tmp(:,:,irpt)) /= 0) then
    atmfrc(:,:,:,:,:,irpt_new) = atmfrc_tmp(:,:,:,:,:,irpt)
    rpt(:,irpt_new) = rpt_tmp(:,irpt)
    wghatm(:,:,irpt_new) = wghatm_tmp(:,:,irpt)
    irpt_new = irpt_new+1
   end if
  end do

  nrpt = nrpt_new
  deallocate (atmfrc_tmp,rpt_tmp,wghatm_tmp)

  write(message, '(a)' )' anaddb    : end of the IFC section '
  call wrtout(6,message,'COLL')

 end if

!DEBUG
!write(6,*)' anaddb : after ifcflag check, ifcflg=',anaddb_dtset%ifcflag,anaddb_dtset%thmflag
!stop
!ENDDEBUG

!**********************************************************************
!**********************************************************************

!Electron-phonon section
 if (anaddb_dtset%elphflag == 1) then

  write(message, '(a,a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&  ' Properties based on electron-phonon coupling ',ch10
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  unitgkk = 90
  open (unit=unitgkk,file=filnam(5),form='unformatted',status='old')
  elph_base_name=trim(filnam(6))
  ddkfilename=trim(filnam(7))

  call elphon(anaddb_dtset%a2fsmear,acell,amu,atmfrc,blkflg,blkqpt,blknrm,blktyp,blkval,&
&  anaddb_dtset%brav,ddkfilename,dielt,anaddb_dtset%dipdip,&
&  dyewq0,anaddb_dtset%elphsmear,anaddb_dtset%elph_fermie,&
&  elph_base_name,anaddb_dtset%enunit,anaddb_dtset%ep_b_min,anaddb_dtset%ep_b_max,&
&  anaddb_dtset%gkk2exist,anaddb_dtset%gkk2write,anaddb_dtset%gkk_rptexist,&
&  anaddb_dtset%gkk_rptwrite,anaddb_dtset%gkqexist,anaddb_dtset%gkqwrite,&
&  anaddb_dtset%phfrqexist,anaddb_dtset%phfrqwrite,anaddb_dtset%prtfsurf,anaddb_dtset%prtnest,gmet,&
&  gprim,anaddb_dtset%ifcflag,anaddb_dtset%ifltransport,&
&  indsym,anaddb_dtset%kptrlatt,&
&  mpert,mpi_enreg,msym,anaddb_dtset%mustar,&
&  natom,nblok,anaddb_dtset%ngqpt(1:3),anaddb_dtset%nqpath,&
&  anaddb_dtset%nqshft,nrpt,nsym,ntypat,anaddb_dtset%qpath,&
&  anaddb_dtset%q1shft,rcan,rmet,rprim,rpt,&
&  symrec,symrel,anaddb_dtset%telphint,anaddb_dtset%tkeepbands,&
&  anaddb_dtset%doscalprod,tnons,trans,typat,&
&  ucvol,unitgkk,wghatm,xred,zeff)

 end if

!**********************************************************************
!**********************************************************************

!Phonon density of states calculation, Start if interatomic forces have been calculated
 if (anaddb_dtset%prtdos==1 .and. anaddb_dtset%ifcflag==1) then
  write(message,'(a,(80a),4a)')ch10,('=',ii=1,80),ch10,ch10,&
&  ' Calculation of phonon density of states with gaussian method, ',ch10
  call wrtout(6,message,'COLL') ; !call wrtout(ab_out,message,'COLL')

  call mkphdos(acell,amu,anaddb_dtset,atmfrc,dielt,dyewq0,filnam(2),gmet,gprim,indsym,&
&  mpert,msym,natom,nrpt,nsym,ntypat,rcan,rmet,rprim,rpt,symrec,symrel,tcpui,&
&  trans,twalli,typat,ucvol,wghatm,xred,zeff) 
 end if

!Phonon density of states and thermodynamical properties calculation 
!Start if interatomic forces and thermal flags are on
 if(anaddb_dtset%ifcflag==1 .and. anaddb_dtset%thmflag/=0) then

  write(message, '(a,a,a)' ) ch10,&
&  ' anaddb   : start phonon density of states calculation',ch10
  call wrtout(6,message,'COLL')
  write(message, '(a,(80a),a,a,a,a,a,a,a,a)' ) ch10,('=',ii=1,80),ch10,ch10,&
&  ' Calculation of phonon density of states, ',ch10,&
&  '    thermodynamical properties, ',ch10,&
&  '    and Debye-Waller factors.',ch10
  call wrtout(ab_out,message,'COLL')

  call timein(tcpu,twall)
  write(message, '(a,f11.3,a,f11.3,a)' )&
&  '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  if(anaddb_dtset%thmflag==1) then
   call thm9(acell,amu,anaddb_dtset,atmfrc,dielt,displ,dyewq0,d2cart,&
&   eigval,eigvec,gmet,gprim,indsym,ab_out,mpert,msym,natom,&
&   nrpt,nsym,ntypat,phfrq,rcan,rmet,rprim,rpt,symrec,symrel,tcpui,&
&   trans,twalli,typat,ucvol,wghatm,xred,zeff)
  else if (anaddb_dtset%thmflag==2) then
   write(message, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&   ch10,' Entering modified part ',ch10
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   themflag = anaddb_dtset%thmflag
   call thm9(acell,amu,anaddb_dtset,atmfrc,dielt,displ,dyewq0,d2cart,&
&   eigval,eigvec,gmet,gprim,indsym,ab_out,mpert,msym,natom,&
&   nrpt,nsym,ntypat,phfrq,rcan,rmet,rprim,rpt,symrec,symrel,tcpui,&
&   trans,twalli,typat,ucvol,wghatm,xred,zeff, themflag, filnam(4), udispl, ufreq)

  end if

 end if

!**********************************************************************
!**********************************************************************

!Now treat the first list of vectors (without non-analyticities)
!MJV NOTE 28/7/2008: the whole of this should be in a subroutine
!for phonon interpolation on 1st list of q vectors
 if(anaddb_dtset%nph1l/=0)then

  write(message, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&  ch10,' Treat the first list of vectors ',ch10
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call timein(tcpu,twall)
  write(message, '(a,f11.3,a,f11.3,a)' )&
&  '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! Write some information in the lwf-formatted file
  if (anaddb_dtset%eivec==3) then
   call outlwf9(acell,iodyn,msym,natom,anaddb_dtset%nph1l,nsym,ntypat,rprim,symrel,typat,xred)
  end if

! open and write header for phonon frequency file
  phfrq_unit=300
  phonon_freq_filename = trim(outfilename_radix)//"_PHFRQ"
  open (unit=phfrq_unit,file=phonon_freq_filename)
  write (phfrq_unit,*) '#'
  write (phfrq_unit,*) '# phonon frequencies (in Ha) on qph1l list of qpoints'
  write (phfrq_unit,*) '#'

  do iphl1=1,anaddb_dtset%nph1l

!  Initialisation of the phonon wavevector
   qphon(:,1)=anaddb_dtset%qph1l(:,iphl1)
   qphnrm(1)=anaddb_dtset%qnrml1(iphl1)

!  Generation of the dynamical matrix in cartesian coordinates
   if(anaddb_dtset%ifcflag==1)then

!   Get d2cart using the interatomic forces and the
!   long-range coulomb interaction through Ewald summation
    write(message, '(a)' )' anaddb    : enter gtdyn9 '
    call wrtout(6,message,'COLL')

    call timein(tcpu,twall)
    write(message, '(a,f11.3,a,f11.3,a)' )&
&    '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
    call wrtout(6,message,'COLL')
    call gtdyn9(acell,atmfrc,dielt,anaddb_dtset%dipdip,&
&    dyewq0,d2cart,gmet,gprim,mpert,natom,&
&    nrpt,qphnrm(1),qphon,rcan,rmet,rprim,rpt,&
&    trans,ucvol,wghatm,xred,zeff)

   else if(anaddb_dtset%ifcflag==0)then

!   Look after the information in the DDB
    rfphon(1:2)=1
    rfelfd(1:2)=0
    rfstrs(1:2)=0
    rftyp=anaddb_dtset%rfmeth
    call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&    qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

!   Copy the dynamical matrix in d2cart
    d2cart(:,1:msize)=blkval(:,:,iblok)

!   Eventually impose the acoustic sum rule
    if (anaddb_dtset%asr==1 .or. anaddb_dtset%asr==2) then
     call asria9(anaddb_dtset%asr,2,d2asr,d2cart,mpert,natom)
    end if

!   Impose acoustic sum rule plus rotational symmetry for 0D and 1D systems    
    if (anaddb_dtset%asr==3 .or. anaddb_dtset%asr==4) then
     call asrprs(anaddb_dtset%asr,2,3,uinvers,vtinvers,singular,d2cart,mpert,natom,rprim,xcart)
    end if
   end if

!  Calculation of the eigenvectors and eigenvalues
!  of the dynamical matrix
   write(message, '(a)' )' anaddb    : enter phfrq3 '
   call wrtout(6,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(6,message,'COLL')
   write(*,*) ' in anaddb with start phfrq3'
   call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&   mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol,xred)

   if(anaddb_dtset%thmflag==3) then
    if(.not.allocated(phfreq)) then
     allocate(phfreq(3*natom,anaddb_dtset%nph1l))
     allocate(eigvect(2,3,natom,3*natom,anaddb_dtset%nph1l))
    end if
    phfreq(:,iphl1) = phfrq(:)
    index1 = 0
    do ii=1,natom
     do jj=1,3
      index1=index1+1
      eigvect(:,:,:,index1,iphl1) = eigvec(:,:,:,jj,ii) 
     end do
    end do
   end if
!  In case eivec == 4, write output files for band2eps
!  (vizualization of phonon band structures)
   if (anaddb_dtset%eivec == 4) then
    call sortph(displ,filnam(4),&
&    ab_out,natom,phfrq,qphnrm(1),qphon,udispl,ufreq)
   end if

!  Write the phonon frequencies
   write(message, '(a)' )' anaddb    : enter prtph3 '
   call wrtout(6,message,'COLL')

   call timein(tcpu,twall)
   write(message, '(a,f11.3,a,f11.3,a)' )&
&   '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
   call wrtout(6,message,'COLL')
!  if (anaddb_dtset%eivec==3) then
   call prtph3(displ,anaddb_dtset%eivec,anaddb_dtset%enunit,iodyn,ab_out,natom,phfrq,qphnrm(1),qphon)
!  else
!  call prtph3(displ,anaddb_dtset%eivec,anaddb_dtset%enunit,-1,ab_out,natom,phfrq,qphnrm,qphon)
!  end if

!  write to file - present version of prtph3 is not atomic enough to do this and
!  has lots of other junk etc...
   write (phfrq_unit,'(i6)',ADVANCE='NO') iphl1
   do ii=1,3*natom
    write (phfrq_unit,'(E16.8,2x)',ADVANCE='NO') phfrq(ii)
   end do
   write (phfrq_unit,*)

!  Determine the symmetries of the phonon mode at Gamma
   if(sum(abs(qphon(:,1)))<qtol)then
    call symph3(ab_out,acell,eigvec,indsym,natom,nsym,phfrq,rprim,symrel)
   end if

  end do
! close unit for phonon frequencies
  close (phfrq_unit)

 end if

!***********************************************************************
!***********************************************************************
!Test thmeig
 if(anaddb_dtset%thmflag==3) then
  call inprep8(dimekb,filnam(5),mband,mblktyp,msym,&
&  natom,nblok,nkpt,ntypat,ddbun,usepaw,vrsddb)

  mpert=natom
  msize=mpert*3*3*mpert

  allocate(blkval2(2,msize,mband,nkpt,nblok),kpnt(3,nkpt,nblok))
  
  call rdddb9(acell,anaddb_dtset%atifc,amu,blkflg,blknrm,blkqpt,&
&  blktyp,blkval,ddbun,dimekb,filnam(5),gmet,gprim,indsym,ab_out,&
&  mband,mpert,msize,msym,&
&  anaddb_dtset%natifc,natom,nblok,nkpt,nsym,ntypat,&
&  occopt,rmet,rprim,symq,symrec,symrel,anaddb_dtset%thmflag,&
&  tnons,typat,ucvol,usepaw,xcart,xred,zion,blkval2,kpnt)

  write(*,*)'Entering thmeig: '
  write(*,*) 'blkval2', blkval2(1,1,1,1,1)
  call thmeig(acell,amu,blkval2,eigvect,filnam(6),kpnt,mband,msize,natom,nkpt,nblok,anaddb_dtset%ntemper,&
&  ntypat,phfreq,blkqpt,rprim,anaddb_dtset%temperinc,anaddb_dtset%tempermin,typat,blknrm,xred)

 end if

!**********************************************************************

!DEBUG
!write(6,*)' anaddb : nph2l,dieflag=',nph2l,anaddb_dtset%dieflag
!stop
!ENDDEBUG
!Now treat the second list of vectors (only at the Gamma point,
!but can include non-analyticities), as well as the
!frequency-dependent dielectric tensor

 if (anaddb_dtset%nlflag > 0) allocate(rsus(3*natom,3,3))
 allocate(fact_oscstr(2,3,3*natom))

 if( nph2l/=0 .or. anaddb_dtset%dieflag==1 )then

  write(message, '(a,(80a),a,a,a,a)' ) ch10,('=',ii=1,80),ch10,&
&  ch10,' Treat the second list of vectors ',ch10
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  call timein(tcpu,twall)
  write(message, '(a,f11.3,a,f11.3,a)' )&
&  '-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! Before examining every direction or the dielectric tensor,
! generates the dynamical matrix at gamma
  qphon(:,1)=0.0d0
  qphnrm(1)=0.0d0

! Generation of the dynamical matrix in cartesian coordinates
  if(anaddb_dtset%ifcflag==1)then

!  Get d2cart using the interatomic forces and the
!  long-range coulomb interaction through Ewald summation
   call gtdyn9(acell,atmfrc,dielt,anaddb_dtset%dipdip,&
&   dyewq0,d2cart,gmet,gprim,mpert,natom,&
&   nrpt,qphnrm(1),qphon,rcan,rmet,rprim,rpt,&
&   trans,ucvol,wghatm,xred,zeff)

  else if(anaddb_dtset%ifcflag==0)then

!  Look after the information in the DDB
   rfphon(1:2)=1
   rfelfd(1:2)=2
   rfstrs(1:2)=0
   rftyp=anaddb_dtset%rfmeth
   call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

!  Copy the dynamical matrix in d2cart
   d2cart(:,1:msize)=blkval(:,:,iblok)

!  Eventually impose the acoustic sum rule
   call asria9(anaddb_dtset%asr,2,d2asr,d2cart,mpert,natom)

!  end of the generation of the dynamical matrix at gamma.
  end if

  if(nph2l/=0)then

!  Examine every wavevector of this list
   do iphl2=1,nph2l

!   Initialisation of the phonon wavevector
    qphon(:,1)=anaddb_dtset%qph2l(:,iphl2)
    qphnrm(1)=anaddb_dtset%qnrml2(iphl2)

!   Calculation of the eigenvectors and eigenvalues
!   of the dynamical matrix
    call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&    mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,rprimd,anaddb_dtset%symdynmat,&
&    symrel,typat,ucvol,xred)

!   Write the phonon frequencies
    call prtph3(displ,anaddb_dtset%eivec,anaddb_dtset%enunit,-1,ab_out,natom,phfrq,qphnrm(1),qphon)

!   Determine the symmetries of the phonon modes at Gamma
    if(sum(abs(qphon(:,1)))<qtol)then
     call symph3(ab_out,acell,eigvec,indsym,natom,nsym,phfrq,rprim,symrel)
    end if

!   Write Raman susceptibilities
    if (anaddb_dtset%nlflag == 1) then
     call ramansus(d2cart,dchide,dchidt,displ,mpert,&
&     natom,phfrq,qphon,qphnrm(1),rsus,ucvol)
    end if

!   Prepare the evaluation of the Lyddane-Sachs-Teller relation
    if(anaddb_dtset%dieflag==1 .and. natom>1)then
     lst(iphl2)=zero
!    The fourth mode should have positive frequency, otherwise,
!    there is an unstability, and the LST relationship should not
!    be evaluated
     if(phfrq(4)>tol6)then
      do ii=4,3*natom
       lst(iphl2)=lst(iphl2)+2*log(phfrq(ii))
      end do
     end if
    end if

   end do ! iphl2
  end if ! nph2l/=0

! The frequency-dependent dielectric tensor (and oscillator strength).
  if (anaddb_dtset%dieflag==1)then

   write(message, '(a,a,a,a,a,a)' )&
&   ' anaddb : the frequency-dependent dielectric tensor (and also once more',&
&   ch10,' the phonons at gamma - without non-analytic part )',ch10,&
&   ch10,' The frequency-dependent dielectric tensor'
   call wrtout(6,message,'COLL')

!  Initialisation of the phonon wavevector
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0

!  Calculation of the eigenvectors and eigenvalues
!  of the dynamical matrix
   call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&   mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol,xred)

!  Write the phonon frequencies (not to ab_out, however)
   call prtph3(displ,0,anaddb_dtset%enunit,-1,-1,natom,phfrq,qphnrm(1),qphon)

!  Evaluation of the oscillator strengths and frequency-dependent
!  dielectric tensor.
   call diel9(amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
&   ab_out,lst,mpert,natom,nph2l,ntypat,phfrq,qtol,typat,ucvol)

!  DEBUG
!  write(6,*)' anaddb : after diel9, dielt_rlx(:,:)=',dielt_rlx(:,:)
!  ENDDEBUG

  end if

! If the electronic dielectric tensor only is needed...
  if (anaddb_dtset%dieflag==2.or.anaddb_dtset%dieflag==3&
&  .or. anaddb_dtset%dieflag==4)then

!  Everything is already in place...
   call diel9(amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
&   ab_out,lst,mpert,natom,nph2l,ntypat,phfrq,qtol,typat,ucvol)

  end if

! End the condition of either nph2l/=0  or  dieflag==1
 end if

!**********************************************************************

!In case nph2l was equal to 0, the electronic dielectric tensor
!has to be computed independently.

 if( anaddb_dtset%dieflag==2 .and. anaddb_dtset%nph2l==0 )then

  write(message, '(a)' )&
&  ' anaddb : nph2l=0, so compute the electronic dielectric tensor independently'
  call wrtout(6,message,'COLL')

! Look after the second derivative matrix at gamma in the DDB
! Note that the information on the dielectric tensor is completely
! independent of the interatomic force constant calculation
  qphon(:,1)=0.0d0
  qphnrm(1)=0.0d0
  rfphon(1:2)=0
  rfelfd(1:2)=2
  rfstrs(1:2)=0
  rftyp=anaddb_dtset%rfmeth
  call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&  qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)

  d2cart(:,1:msize)=blkval(:,:,iblok)

! Print the electronic dielectric tensor
  call diel9(amu,anaddb_dtset,dielt_rlx,displ,d2cart,epsinf,fact_oscstr,&
&  ab_out,lst,mpert,natom,nph2l,ntypat,phfrq,qtol,typat,ucvol)

! DEBUG
! write(6,*)' anaddb : after third diel9, dielt_rlx(:,:)=',dielt_rlx(:,:)
! ENDDEBUG

 end if

!**********************************************************************

!Compute the electrooptic tensor

 if (anaddb_dtset%nlflag == 1) then

! In case dieflag = 2, recompute phonon frequencies and
! eigenvectors without non-analyticity
  if (anaddb_dtset%dieflag == 2) then
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0
   call phfrq3(amu,displ,d2cart,eigval,eigvec,indsym,&
&   mpert,msym,natom,nsym,ntypat,phfrq,qphnrm(1),qphon,rprimd,anaddb_dtset%symdynmat,symrel,typat,ucvol,xred)
  end if

  rsus(:,:,:) = 0_dp
  call ramansus(d2cart,dchide,dchidt,displ,mpert,&
&  natom,phfrq(1),qphon,qphnrm(1),rsus,ucvol)

  call electrooptic(dchide,anaddb_dtset%dieflag,epsinf,&
&  fact_oscstr,natom,phfrq,anaddb_dtset%prtmbm,rsus,ucvol)

 end if  ! condition on nlflag and dieflag

 deallocate(fact_oscstr)
 if (anaddb_dtset%nlflag > 0) deallocate(dchide,dchidt,rsus)

!**********************************************************************

!here treating the internal strain tensors at Gamma point
 if(anaddb_dtset%instrflag/=0)then

  write(message, '(a,a,(80a),a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&  ' Calculation of the internal-strain  tensor',ch10
  call wrtout(6,message,'coll')
  call wrtout(ab_out,message,'coll')

  call timein(tcpu,twall)
  write(message,'(a,f11.3,a,f11.3,a)')&
&  '-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'

  call wrtout(6,message,'coll')
  call wrtout(ab_out,message,'coll')

! allocate(instrain(3*natom,6))
  if(anaddb_dtset%instrflag==1)then
   write(message,'(a)' )&
&   ' anaddb : instrflag=1, so extract the internal strain constant from the 2DTE'
   call wrtout(6,message,'coll')

!  looking after the no. of blok that contians internal strain tensor
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0
   rfphon(1:2)=0
   rfelfd(1:2)=0
   rfstrs(1:2)=3
   rftyp=anaddb_dtset%rfmeth
   call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
!  then print the internal stain tensor
!  write(ab_out,'(/,a,i6)')'iblok is',iblok
   call instr9(blkval,iblok,instrain,ab_out,mpert,natom,nblok)
  end if
 end if
!end the part for internal strain

!**********************************************************************

!here treating the elastic tensors at Gamma Point
 if(anaddb_dtset%elaflag/=0)then
  write(message, '(a,a,(80a),a,a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&  ' Calculation of the elastic and compliances tensor (Voigt notation)',ch10
  call wrtout(6,message,'coll')
  call wrtout(ab_out,message,'coll')

  call timein(tcpu,twall)
  write(message,'(a,f11.3,a,f11.3,a)')&
&  '-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'

  call wrtout(6,message,'coll')
  call wrtout(ab_out,message,'coll')


  if(anaddb_dtset%elaflag==1 .or.anaddb_dtset%elaflag==2&
&  .or. anaddb_dtset%elaflag==3 .or.anaddb_dtset%elaflag==4&
&  .or. anaddb_dtset%elaflag==5)then
   write(message,'(a)' )&
&   ' anaddb : so extract the elastic constant from the 2DTE'
   call wrtout(6,message,'coll')

!  look after the blok no. that contains the stress tensor
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0
   rfphon(1:2)=0
   rfelfd(1:2)=0
   rfstrs(1:2)=0
   rftyp=4
   call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
   iblok_stress=iblok

!  DEBUG
!  check the iblok number containing first order derivative
!  write(6,'(/,a,/)')'iblok_stress number'
!  write(6,'(i)')iblok_stress
!  ENDDEBUG

!  look after the blok no.iblok that contains the elastic tensor
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0
   rfphon(1:2)=0
   rfelfd(1:2)=0
   rfstrs(1:2)=3
!  for both diagonal and shear parts
   rftyp=anaddb_dtset%rfmeth
   call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
!  write(ab_out,'(/,a,i6)')'iblok is',iblok
!  print the elastic tensor
   call elast9(anaddb_dtset,blkval,elast,iblok,iblok_stress,instrain,ab_out,mpert,&
&   natom,nblok,ucvol)
  end if
 end if
!ending the part for elastic tensors

!**********************************************************************

!here treating the piezoelectric tensor at Gamma Point
 if(anaddb_dtset%piezoflag/=0 .or. anaddb_dtset%dieflag==4&
& .or. anaddb_dtset%elaflag==4)then
  write(message, '(a,a,(80a),a,a,a,a,a)') ch10,('=',ii=1,80),ch10,ch10,&
&  ' Calculation of the tensor related to piezoelectric effetc',ch10,&
&  '  (Elastic indices in Voigt notation)',ch10
  call wrtout(6,message,'coll')
  call wrtout(ab_out,message,'coll')

  call timein(tcpu,twall)
  write(message,'(a,f11.3,a,f11.3,a)')&
&  '-begin at tcpu',tcpu-tcpui,'   and twall',twall-twalli,'sec'

  call wrtout(6,message,'coll')
  call wrtout(ab_out,message,'coll')

  if(anaddb_dtset%piezoflag==1 .or.anaddb_dtset%piezoflag==2&
&  .or.anaddb_dtset%piezoflag==3 .or. anaddb_dtset%piezoflag==4&
&  .or.anaddb_dtset%piezoflag==5 .or. anaddb_dtset%piezoflag==6&
&  .or.anaddb_dtset%piezoflag==7 .or. anaddb_dtset%dieflag==4&
&  .or.anaddb_dtset%elaflag==4)then
   write(message,'(a)' )&
&   ' anaddb : extract the piezoelectric constant from the 2DTE'
   call wrtout(6,message,'coll')
!  looking for the gamma point block
   qphon(:,1)=0.0d0
   qphnrm(1)=0.0d0
   rfphon(1:2)=0
   rfelfd(1:2)=0
   rfstrs(1:2)=3
!  for both diagonal and shear parts
   rftyp=anaddb_dtset%rfmeth

   call gtblk9(blkflg,blknrm,blkqpt,blktyp,iblok,mpert,msize,natom,nblok,&
&   qphon,qphnrm,qtol,rfphon,rfelfd,rfstrs,rftyp)
!  write(ab_out,'(/,a,i6)')'iblok is',iblok
!  then print out the piezoelectric constants

!  DEBUG
!  write(6,*)' anaddb : before piezo9, dielt_rlx(:,:)=',dielt_rlx(:,:)
!  ENDDEBUG

   call piezo9(anaddb_dtset,blkval,dielt_rlx,elast,iblok,instrain,ab_out,mpert,&
&   natom,nblok,piezo,ucvol)
  end if
 end if

!**********************************************************************

 deallocate(anaddb_dtset%atifc)
 deallocate(anaddb_dtset%iatfix)
 deallocate(amu,blkflg,blknrm,blkqpt,blktyp,blkval)
 deallocate(displ,dyewq0,d2asr,d2cart,eigval,eigvec,indsym,instrain)
 deallocate(lst,phfrq,rcan,trans,typat,xcart,zeff,zion)
 if(anaddb_dtset%ifcflag/=0)deallocate(atmfrc,rpt,wghatm)

 deallocate(anaddb_dtset%qnrml1)
 deallocate(anaddb_dtset%qnrml2)
 deallocate(anaddb_dtset%qph1l)
 deallocate(anaddb_dtset%qph2l)

 call timein(tcpu,twall)
 write(message, '(a,(80a),a,a,a,f11.3,a,f11.3,a,a,a,a)' ) ch10,&
& ('=',ii=1,80),ch10,ch10,&
& '+Total cpu time',tcpu-tcpui,&
& '  and wall time',twall-twalli,' sec',ch10,ch10,&
& ' anaddb : the run completed succesfully.'
 call wrtout(6,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 close(ab_out)

 end program anaddb
!!***
