!{\src2tex{textfont=tt}}
!!****p* ABINIT/mrgddb
!! NAME
!! mrgddb
!!
!! FUNCTION
!! This code merges the derivative databases.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
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
!! The heading of the constituted database is read,
!! then the heading of the temporary database to be added is read,
!! the code check their compatibility, and create a new
!! database that mixes the old and the temporary ones.
!! This process can be iterated.
!! The whole database will be stored in
!! central memory. One could introduce a third mode in which
!! only the temporary DDB is in central memory, while the
!! input DDB is read twice : first to make a table of blocks,
!! counting the final number of blocks, and second to merge
!! the two DDBs. This would save memory.
!!
!! PARENTS
!!
!! CHILDREN
!!      blok8,cmpar8,herald,init8,inprep8,ioddb8,leave_new,psddb8,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program mrgddb

 use defs_basis
 use defs_infos


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_16response
 use interfaces_17ddb
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
!
! Set array dimensions
!  mddb=maximum number of databases (cannot be made dynamic)
!  msym=maximum number of symmetry elements in space group
 integer,parameter :: mddb=500,msym=192
!Define input and output unit numbers:
 integer,parameter :: ddbun=2,unit00=1
 integer :: choice,dimekb,dimekb_tmp,fullinit,fullinit8,iblok,iblok1,iblok2
 integer :: iddb,ii,intxc,intxc8,iscf,iscf8,ixc,ixc8,lmnmax,lnmax,matom,mband
 integer :: mband_tmp,mblktyp,mblok,mkpt,mpert,msize,msize_tmp,mtypat,natom
 integer :: natom8,nblok,nblok8,nblokt,nddb,nkpt,nkpt8,nline,nq,nspden,nspden8
 integer :: nspinor,nspinor8,nsppo8,nsppol,nsym,nsym8,ntypat,ntypat8,nunit
 integer :: occop8,occopt,telphon,tmerge,usepaw,usepaw_tmp,useylm,vrsddb
 integer :: eatpol(2),ngfft(18),ngfft8(18),rfdir(3),symafm(msym),symafm8(msym)
 integer :: symre8(3,3,msym),symrel(3,3,msym)
 integer,allocatable :: blkflg(:,:),blktyp(:),indlmn(:,:,:),lloc(:),mgblok(:)
 integer,allocatable :: nband(:),nband8(:),pspso(:),typat(:),typat8(:)
 real(dp),parameter :: qtol=2.0d-8
 real(dp) :: diff,dilatmx,dilatmx8,ecut,ecut8,ecutsm,ecutsm8,kptnr8,kptnrm
 real(dp) :: sciss,sciss8,tolwf8,tolwfr,tphysel,tphysel8,tsmear,tsmear8
 real(dp) :: acell(3),acell8(3),qpt(3),rprim(3,3),rprim8(3,3),tnons(3,msym)
 real(dp) :: tnons8(3,msym)
 real(dp),allocatable :: amu(:),amu8(:),blknrm(:,:),blkqpt(:,:),blkval(:,:,:)
 real(dp),allocatable :: blkval2(:,:,:,:,:)
 real(dp),allocatable :: ekb(:,:),ekb8(:,:),kpt(:,:),kpt8(:,:),kpnt(:,:,:),occ(:),occ8(:)
 real(dp),allocatable :: spinat(:,:),spinat8(:,:),vel(:,:),wtk(:),wtk8(:)
 real(dp),allocatable :: xcart(:,:),xred(:,:),xred8(:,:),zion(:),zion8(:)
 real(dp),allocatable :: znucl(:),znucl8(:)
 character(len=24) :: codename
 character(len=fnlen) :: dscrpt,dummy
 character(len=fnlen) :: filnam(mddb+1)
 character(len=strlen) :: string
 character(len=500) :: message

!******************************************************************
!BEGIN EXECUTABLE SECTION

 codename='MRGDDB'//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

!Initialise the code : write heading,
!read names of files, operating mode (also check its value),
!and short description of new database.
 call init8(dscrpt,filnam,mddb,nddb)

!Set the ddb version
 vrsddb=010929

!Evaluate the maximal dimensions of arrays
 dimekb=0 ; matom=0 ; mband=0  ; mblok=0 ; mkpt=0
 msize=0  ;mtypat=0 ; usepaw=0
 do iddb=1,nddb
  call inprep8(dimekb_tmp,filnam(iddb+1),mband_tmp,mblktyp,&
&  msym,natom,nblok,nkpt,ntypat,ddbun,usepaw_tmp,vrsddb)
  dimekb=max(dimekb,dimekb_tmp)
  matom=max(matom,natom)
  mband=max(mband,mband_tmp)
  mblok=mblok+nblok
  mkpt=max(mkpt,nkpt)
  mtypat=max(mtypat,ntypat)
  usepaw=max(usepaw,usepaw_tmp)
 end do

!Check the value of usepaw
 if (usepaw==1) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' mrgddb: BUG -',ch10,&
&  '  Paw not yet allowed !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 mpert=matom+6
 msize=3*mpert*3*mpert
 if(mblktyp==3)msize=msize*3*mpert

!Allocate arrays
 allocate(blkflg(msize,mblok),blktyp(mblok))
 allocate(lloc(mtypat),mgblok(mblok))
 allocate(nband(mkpt),nband8(mkpt))
 allocate(typat(matom),typat8(matom))
 allocate(amu(mtypat),amu8(mtypat))
 allocate(blkqpt(9,mblok),blknrm(3,mblok))
 allocate(blkval(2,msize,mblok),ekb(dimekb,mtypat))
 allocate(ekb8(dimekb,mtypat),kpt(3,mkpt),kpt8(3,mkpt))
 allocate(indlmn(6,dimekb,mtypat),pspso(mtypat))
 allocate(occ(mband*mkpt),occ8(mband*mkpt))
 allocate(spinat(3,matom),spinat8(3,matom))
 allocate(vel(3,matom),wtk(mkpt),wtk8(mkpt),xcart(3,matom))
 allocate(xred8(3,matom),xred(3,matom))
 allocate(znucl(mtypat),znucl8(mtypat))
 allocate(zion(mtypat),zion8(mtypat))

!This is needed to read the DDBs in the old format
 symafm(:)=1 ; symafm8(:)=1
 if(mtypat>=1)then
  pspso(:)=0
  znucl(:)=zero ; znucl8(:)=zero
  ekb(:,:)=zero ; ekb8(:,:)=zero
 end if
 if(matom>=1)then
  spinat(:,:)=zero ; spinat8(:,:)=zero
 end if

!**********************************************************************

!Read the first database

 if(nddb==1)then

  write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&  ' mrgddb : ERROR -',ch10,&
&  '  The initialisation mode of MRGDDB, that uses nddb=1,',&
&  '  has been disabled in version 2.1 of ABINIT.',&
&  '  Action : you should use DDBs that include the symmetry',&
&  '  information (and that can be used and merged without',&
&  '  initialisation), or you should use ABINITv2.0.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')

 else if(nddb>=2)then

! Open the first derivative database file
! and read the preliminary information
  write(6,*)' read the input derivative database information'
  choice=1
  nunit=ddbun
  call ioddb8 (choice,dummy,filnam(2),matom,mband,&
&  mkpt,msym,mtypat,nunit,vrsddb,&
&  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,&
&  tphysel,tsmear,typat,wtk,xred,zion,znucl)

  telphon=0
  if(telphon==1)then
   allocate(blkval2(2,msize,nband(1),mkpt,mblok),kpnt(3,mkpt,mblok))
  end if 
! DEBUG
! write(6,*)' rprim=',rprim
! ENDDEBUG

! Read the psp information of the input DDB
  lmnmax=dimekb;lnmax=dimekb;useylm=0
  call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,lnmax,&
&  nblok,ntypat,nunit,pspso,usepaw,useylm,vrsddb)

  if(nblok>=1)then
!  Read the blocks from the input database.
   write(message, '(a,i5,a)' ) ' read ',nblok,' blocks from the input DDB '
   call wrtout(6,message,'COLL')
   choice=1
   nunit=ddbun
   do iblok=1,nblok
    if(telphon==1)then
     call blok8(blkflg(1,iblok),blknrm(1,iblok),blkqpt(1,iblok), &
&     blktyp(iblok),blkval(1,1,iblok),choice,nband(1),mpert,&
&     msize,natom,nkpt,nunit,blkval2(:,:,:,:,iblok),kpnt(:,:,iblok))
    else
     call blok8(blkflg(1,iblok),blknrm(1,iblok),blkqpt(1,iblok), &
&     blktyp(iblok),blkval(1,1,iblok),choice,nband(1),mpert,&
&     msize,natom,nkpt,nunit)
    end if
!   Setup merged indicator
    mgblok(iblok)=0
   end do
  else
   write(message, '(a)' )' No bloks in the first ddb '
   call wrtout(6,message,'COLL')
  end if
! Close the first ddb
  close(ddbun)

 end if

!*********************************************

!In case of merging of DDBs, iterate the reading
 do iddb=2,nddb

! Open the corresponding input DDB,
! and read the database file informations
  write(message, '(a,a,i6)' )ch10,&
&  ' read the input derivative database number',iddb
  call wrtout(6,message,'COLL')
  choice=1
  nunit=ddbun
  call ioddb8 (choice,dummy,filnam(iddb+1),matom,mband,&
&  mkpt,msym,mtypat,nunit,vrsddb,&
&  acell8,amu8,dilatmx8,ecut8,ecutsm8,intxc8,iscf8,ixc8,kpt8,kptnr8,&
&  natom8,nband8,ngfft8,nkpt8,nspden8,nspinor8,nsppo8,nsym8,ntypat8,occ8,occop8,&
&  rprim8,sciss8,spinat8,symafm8,symre8,tnons8,tolwf8,&
&  tphysel8,tsmear8,typat8,wtk8,xred8,zion8,znucl8)

! Read the psp information of the input DDB
  call psddb8 (choice,dimekb,ekb8,fullinit8,indlmn,lmnmax,lnmax,&
&  nblok8,ntypat8,nunit,pspso,usepaw,useylm,vrsddb)

! Compare the current DDB and input DDB information.
! In case of an inconsistency, halt the execution.
  write(message, '(a)' )' compare the current and input DDB information'
  call wrtout(6,message,'COLL')
! DEBUG
! write(6,*)' occopt=',occopt
! ENDDEBUG

! Should also compare indlmn and pspso ... but suppose that
! the checking of ekb is enough for the psps.
! Should also compare many other variables ... this is still
! to be done ...
  call cmpar8 (acell,acell8,amu,amu8,dimekb,ecut,ecut8,ekb,ekb8,&
&  fullinit,fullinit8,iscf,iscf8,ixc,ixc8,kpt,kpt8,kptnrm,kptnr8,&
&  natom,natom8,nband,nband8,ngfft,ngfft8,nkpt,nkpt8,&
&  nsppol,nsppo8,nsym,nsym8,ntypat,ntypat8,occ,occ8,&
&  occopt,occop8,rprim,rprim8,sciss,sciss8,symrel,symre8,&
&  tnons,tnons8,tolwfr,tolwf8,typat,typat8,&
&  usepaw,wtk,wtk8,xred,xred8,zion,zion8)
! DEBUG
! write(6,*)' occopt=',occopt
! ENDDEBUG

  write(message, '(a,a)' )' Will try to merge this input DDB with',&
&  ' the current one.'
  call wrtout(6,message,'COLL')

! First estimate of the total number of bloks, and error
! message if too large
  write(message, '(a,i5)' ) ' Current number of bloks =',nblok
  call wrtout(6,message,'COLL')
  write(message, '(a,i5,a)' )&
&  ' Will read ',nblok8,' blocks from the input DDB '
  call wrtout(6,message,'COLL')
  nblokt=nblok+nblok8
  if(nblokt>mblok)then
   write(message, '(a,a,a,i5,a,a,a,i5,a)' )&
&   ' mrgddb : ERROR -',ch10,&
&   '  The expected number of blocks',nblokt,' is larger than',ch10,&
&   'the maximum number of blocks',mblok,'.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Read the bloks from the temporary database, and close it.
! Also setup the merging indicator
  choice=1
  nunit=ddbun
  do iblok=nblok+1,nblokt
   if(telphon==1)then
    call blok8(blkflg(1,iblok),blknrm(1,iblok),blkqpt(1,iblok),&
&    blktyp(iblok),blkval(1,1,iblok),choice,nband(1),mpert,&
&    msize,natom,nkpt8,nunit,blkval2(:,:,:,:,iblok),kpnt(:,:,iblok))
   else
    call blok8(blkflg(1,iblok),blknrm(1,iblok),blkqpt(1,iblok),&
&    blktyp(iblok),blkval(1,1,iblok),choice,nband(1),mpert,&
&    msize,natom,nkpt8,nunit)
   end if
   mgblok(iblok)=0
  end do
  close(ddbun)

  nblok=nblokt
  write(message, '(a,i5)' ) ' Now, current number of bloks =',nblok
  call wrtout(6,message,'COLL')

 end do

 write(message, '(a)' )' All DDBs have been read '
 call wrtout(6,message,'COLL')

!*********************************************************

!Check the equality of blocks, and eventually merge them

 if(nblok>=1)then
  write(message, '(a)' )' check the equality of blocks, and eventually merge '
  call wrtout(6,message,'COLL')
  do iblok2=2,nblok
   do iblok1=1,iblok2-1
    tmerge=0

!   Check the block type identity
    if(blktyp(iblok1)==blktyp(iblok2))then

!    Check the wavevector identities
     tmerge=1
     if(blktyp(iblok1)==1.or.blktyp(iblok1)==2)then
      nq=1
     else if(blktyp(iblok1)==3)then
!     Note : do not merge permutation related elements ....
      nq=3
     else if(blktyp(iblok1)==4 .or. blktyp(iblok1)==0)then
      nq=0
     else if(blktyp(iblok1)==5)then
      nq=0
      tmerge=0
     end if
     if(nq/=0)then
      do ii=1,nq
       diff=blkqpt(1+3*(ii-1),iblok1)/blknrm(ii,iblok1)&
&       -blkqpt(1+3*(ii-1),iblok2)/blknrm(ii,iblok2)
       if(abs(diff)>qtol)tmerge=0
       diff=blkqpt(2+3*(ii-1),iblok1)/blknrm(ii,iblok1)&
&       -blkqpt(2+3*(ii-1),iblok2)/blknrm(ii,iblok2)
       if(abs(diff)>qtol)tmerge=0
       diff=blkqpt(3+3*(ii-1),iblok1)/blknrm(ii,iblok1)&
&       -blkqpt(3+3*(ii-1),iblok2)/blknrm(ii,iblok2)
       if(abs(diff)>qtol)tmerge=0
      end do ! ii
     end if

!    Now merges,
     if(tmerge==1)then
      write(message, '(a,i5,a,i5)' )&
&      ' merge block #',iblok2,' to block #',iblok1
      call wrtout(6,message,'COLL')
      mgblok(iblok2)=1
      do ii=1,msize
       if(blkflg(ii,iblok2)==1)then
        blkflg(ii,iblok1)=1
        blkval(1,ii,iblok1)=blkval(1,ii,iblok2)
        blkval(2,ii,iblok1)=blkval(2,ii,iblok2)
       end if
      end do
     end if

    end if
   end do
  end do

! Count the final number of bloks
  tmerge=0
  do ii=1,nblok
   if(mgblok(ii)==1)tmerge=tmerge+1
  end do
  nblok=nblok-tmerge

! Summarize the merging phase
  write(message, '(i6,a,i6,a)' )&
&  tmerge,' blocks are merged; the new DDB will have ',nblok,&
&  ' blocks.'
  call wrtout(6,message,'COLL')

! End the condition on existence of more than one blok in current DDB
 end if

!**********************************************************************

!Open the output database, then
!Write the preliminary informations
 write(message, '(a,a)' )' open the output database, write the',&
& ' preliminary information '
 call wrtout(6,message,'COLL')

!DEBUG
!write(6,*)' occopt=',occopt
!ENDDEBUG

 choice=2
 nunit=ddbun
 call ioddb8 (choice,dscrpt,filnam(1),matom,mband,&
& mkpt,msym,mtypat,nunit,vrsddb,&
& acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
& natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
& rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,&
& tphysel,tsmear,typat,wtk,xred,zion,znucl)

!DEBUG
!write(6,*)' mrgddb : after ioddb8 '
!stop
!ENDDEBUG

!Write the psp information in the output DDB
!as well as the value of the number of blocks.
 write(message, '(a)' )' write the psp information '
 call wrtout(6,message,'COLL')
 fullinit=1
 call psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,lnmax,&
& nblok,ntypat,nunit,pspso,usepaw,useylm,vrsddb)

!DEBUG
!write(6,*)' mrgddb : after psddb8, nddb= ',nddb
!stop
!ENDDEBUG

 if(nddb>1)then

! Write the whole database
  write(message, '(a)' )' write the DDB '
  call wrtout(6,message,'COLL')
  choice=2
  nunit=ddbun
  do iblok=1,nblok+tmerge
   if(mgblok(iblok)==0)then
    write(6, '(a,i4)' ) ' Write bloc number',iblok
    if(telphon==1)then
     call blok8(blkflg(1,iblok),blknrm(1,iblok),blkqpt(1,iblok), &
&     blktyp(iblok),blkval(1,1,iblok),choice,nband(1),mpert,&
&     msize,natom,nkpt,nunit,blkval2(:,:,:,:,iblok),kpnt(:,:,iblok))
    else
     call blok8(blkflg(1,iblok),blknrm(1,iblok),blkqpt(1,iblok), &
&     blktyp(iblok),blkval(1,1,iblok),choice,nband(1),mpert,&
&     msize,natom,nkpt,nunit)
    end if
   else
    write(message, '(a,i4,a)' )&
&    ' Bloc number',iblok,' was merged, so do not write it'
    call wrtout(6,message,'COLL')
   end if
  end do

! Also write summary of bloks at the end
  write(ddbun, '(/,a)' )' List of bloks and their characteristics '
  choice=3
  nunit=ddbun
  do iblok=1,nblok+tmerge
   if(mgblok(iblok)==0)then
    if(telphon==1)then
     call blok8(blkflg(1,iblok),blknrm(1,iblok),blkqpt(1,iblok),&
&     blktyp(iblok),blkval(1,1,iblok),choice,nband(1),mpert,&
&     msize,natom,nkpt,nunit,blkval2(:,:,:,:,iblok),kpnt(:,:,iblok))
    else
     call blok8(blkflg(1,iblok),blknrm(1,iblok),blkqpt(1,iblok),&
&     blktyp(iblok),blkval(1,1,iblok),choice,nband(1),mpert,&
&     msize,natom,nkpt,nunit)
    end if
   end if
  end do

 end if

 close (ddbun)

!*********************************************************************

 write(message, '(a)' )'+mrgddb : the run completed successfully '
 call wrtout(6,message,'COLL')

 end program
!!***
