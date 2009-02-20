!{\src2tex{textfont=tt}}
!!****f* ABINIT/inprep8
!!
!! NAME
!! inprep8
!!
!! FUNCTION
!! Open Derivative DataBase, then reads the variables that
!! must be known in order to dimension the arrays before complete reading
!!
!! Note : only one processor read or write the DDB.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! character(len=fnlen) filnam: name of input or output file
!! msym=maximum number of symmetries
!! unddb=unit number for input or output
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.  Have to agree with vrsinddb of this routine.
!!
!! OUTPUT
!! dimekb=dimension of ekb (for the time being, only for norm-
!!                          conserving psps)
!! mband=maximum number of bands
!! mblktyp=largest block type
!! natom=number of atoms
!! nblok=number of bloks in the DDB
!! nkpt=number of k points
!! ntypat=number of atom types
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! PARENTS
!!      anaddb,mrgddb
!!
!! CHILDREN
!!      chknm8,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine inprep8 (dimekb,filnam,mband,mblktyp,msym,natom,nblok,nkpt,&
& ntypat,unddb,usepaw,vrsddb)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: msym,unddb,vrsddb
 integer,intent(out) :: dimekb,mband,mblktyp,natom,nblok,nkpt,ntypat,usepaw
 character(len=fnlen),intent(in) :: filnam

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer,parameter :: vrsio8=010929,vrsio8_old=990527
 integer,save :: count=0
 integer :: bantot,blktyp,ddbvrs,iband,iblok,iekb,ii,ikpt,iline,im,ios,iproj
 integer :: itypat,itypat0,jekb,lmnmax,mproj,mpsang,nekb,nelmts,nsppol,nsym
 integer :: occopt,pspso0
 logical :: testn,testv
 character(len=12) :: string
 character(len=32) :: blkname
 character(len=500) :: message
 character(len=6) :: name_old
 character(len=80) :: rdstring
!arrays
 integer,allocatable :: nband(:)
 character(len=9) :: name(9)

! *********************************************************************

!DEBUG
!write(6,*)' inprep8 : enter'
!count=count+1
!write(6,*)' count=',count
!ENDDEBUG

!Check inprep8 version number (vrsio8) against mkddb version number
!(vrsddb)
 if (vrsio8/=vrsddb) then
  write(message, '(a,a,a,i10,a,a,i10,a)' )&
&  ' inprep8: BUG -',ch10,&
&  '  The input/output DDB version number=',vrsio8,ch10,&
&  '  is not equal to the DDB version number=',vrsddb,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Open the input derivative database.
 write(6, '(a,a)') ' inprep8 : open file ',trim(filnam)
 open (unit=unddb,file=filnam,status='old',form='formatted')

!Check the compatibility of the input DDB with the DDB code
 read (unddb,*)
 read (unddb,*)
 read (unddb, '(20x,i10)' )ddbvrs
 if(ddbvrs/=vrsio8 .and. ddbvrs/=vrsio8_old)then
  write(message, '(3a,i10,3a,i10,a,i10,a)' )&
&  ' inprep8 : BUG - ',ch10,&
&  '  The input DDB version number=',ddbvrs,' does not agree',ch10,&
&  '  with the allowed code DDB version numbers,',vrsio8,' and ',vrsio8_old,' .'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Read the 4 n-integers, also testing the names of data,
!and checking that their value is acceptable.
!This is important to insure that any array has a sufficient
!dimension.
 read (unddb,*)
 read (unddb,*)
 read (unddb,*)
 testn=.true.
 testv=.true.
!1. natom
 if(ddbvrs==vrsio8)then
  read (unddb, '(1x,a9,i10)' )name(1),natom
 else
  read (unddb, '(1x,a6,i10)' )name_old,natom ; name(1)='   '//name_old
 end if
 if(name(1)/='    natom')testn=.false.
 if(natom<=0)testv=.false.
!2. nkpt
 if(ddbvrs==vrsio8)then
  read (unddb, '(1x,a9,i10)' )name(2),nkpt
 else
  read (unddb, '(1x,a6,i10)' )name_old,nkpt ; name(2)='   '//name_old
 end if
 if(name(2)/='     nkpt')testn=.false.
 if(nkpt <=0)testv=.false.
!3. nsppol
 if(ddbvrs==vrsio8)then
  read (unddb, '(1x,a9,i10)' )name(3),nsppol
 else
  read (unddb, '(1x,a6,i10)' )name_old,nsppol ; name(3)='   '//name_old
 end if
 if(name(3)/='   nsppol')testn=.false.
 if(nsppol<=0.or.nsppol>2)testv=.false.
!4. nsym
 if(ddbvrs==vrsio8)then
  read (unddb, '(1x,a9,i10)' )name(4),nsym
 else
  read (unddb, '(1x,a6,i10)' )name_old,nsym ; name(4)='   '//name_old
 end if
 if(name(4)/='     nsym')testn=.false.
 if(nsym <=0.or.nsym >msym )testv=.false.
!5. ntypat
 if(ddbvrs==vrsio8)then
  read (unddb, '(1x,a9,i10)' )name(5),ntypat
 else
  read (unddb, '(1x,a6,i10)' )name_old,ntypat ; name(5)='   '//name_old
 end if
 if(ntypat<=0)testv=.false.
!6. occopt
!Before reading nband, the last parameters that define
!the dimension of some array, need to know what is their
!representation, given by occopt
 if(ddbvrs==vrsio8)then
  read (unddb, '(1x,a9,i10)' )name(6),occopt
 else
  read (unddb, '(1x,a6,i10)' )name_old,occopt ; name(6)='   '//name_old
 end if
 if(name(6)/='   occopt')testn=.false.
 if(occopt<0.or.occopt>7)testv=.false.
!Message if the names or values are not right
 if (.not.testn.or..not.testv) then
  write(message, '(a,a)' )' inprep8 : An error has been found in the',&
&  ' positive n-integers contained in the DDB : '
  call wrtout(6,message,'COLL')
  write(message, '(a)' )   '     Expected                      Found     '
  call wrtout(6,message,'COLL')
  write(message, '(a,i9,a,a,a,i10)' )&
&  '    natom , larger than',0,'    ',trim(name(1)),' =',natom
  call wrtout(6,message,'COLL')
  write(message, '(a,i9,a,a,a,i10)' )&
&  '    nkpt  , larger than',0,'    ',trim(name(2)),' =',nkpt
  call wrtout(6,message,'COLL')
  write(message, '(a,i1,a,a,a,i10)' )&
&  '    nsppol, either    1 or     ',2,'    ',trim(name(3)),' =',nsppol
  call wrtout(6,message,'COLL')
  write(message, '(a,i10,a,a,a,i10)' )&
&  '    nsym  , lower than',msym,'    ',trim(name(4)),' =',nsym
  call wrtout(6,message,'COLL')
  write(message, '(a,i9,a,a,a,i10)' )&
&  '    ntypat , larger than',0,'    ',trim(name(5)),' =',ntypat
  call wrtout(6,message,'COLL')
  write(message, '(a,a,a,i10)' )&
&  '    occopt,     equal to 0,1 or 2   ',trim(name(6)),' =',occopt
  call wrtout(6,message,'COLL')
  write(message, '(a,a,a)' )&
&  ' inprep8 : ERROR -',ch10,&
&  '  See the error message above.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!One more set of parameters define the dimensions of the
!array : nband. Morever, it depends on occopt and
!nkpt, and has to be
!tested after the test on nkpt is performed.
!7. nband
!DEBUG
!write(6,*)' inprep8 : nkpt=',nkpt
!if(count==2)stop
!ENDDEBUG
 allocate(nband(nkpt))
 if(occopt==2)then
  im=12
  do iline=1,(nkpt+11)/12
   if(iline==(nkpt+11)/12)im=nkpt-12*(iline-1)
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,5x,12i5)' )name(1),&
&    (nband((iline-1)*12+ii),ii=1,im)
   else
    read (unddb, '(1x,a6,5x,12i5)' )name_old,&
&    (nband((iline-1)*12+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
    call chknm8(name(1),'    nband')
   else
    call chknm8(name(1),'         ')
   end if
  end do
 else
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),nband(1)
  else
   read (unddb, '(1x,a6,i10)' )name_old,nband(1) ; name(1)='   '//name_old
  end if
  call chknm8(name(1),'    nband')
  if(nkpt>1)then
   do ikpt=2,nkpt
    nband(ikpt)=nband(1)
   end do
  end if
 end if

!Check all nband values, and sum them
 bantot=0
 do ikpt=1,nkpt
  if(nband(ikpt)<0)then
   write(message, '(a,a,a,i4,a,i4,3a)' )&
&   ' inprep8 : ERROR -',ch10,&
&   '  For ikpt = ',ikpt,'  nband = ',nband(ikpt),' is negative.',ch10,&
&   '  Action : correct your DDB.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  bantot=bantot+nband(ikpt)
 end do

 mband=maxval(nband(:))

!Skip the rest of variables
!8. acell
 read (unddb,*)
!9. amu
 do iline=1,(ntypat+2)/3
  read (unddb,*)
 end do
!10. dilatmx
 if(ddbvrs==vrsio8)then
  read (unddb,*)
 end if
!11. ecut
 read (unddb,*)
!12. ecutsm
 if(ddbvrs==vrsio8)then
  read (unddb,*)
 end if
!13. intxc
 if(ddbvrs==vrsio8)then
  read (unddb,*)
 end if
!14. iscf
 read (unddb,*)
!15. ixc
 read (unddb,*)
!16. kpt
 do iline=1,nkpt
  read (unddb,*)
 end do
!17. kptnrm
 read (unddb,*)
!18. ngfft
 read (unddb,*)
!19. nspden
 if(ddbvrs==vrsio8)then
  read (unddb,*)
 end if
!20. nspinor
 if(ddbvrs==vrsio8)then
  read (unddb,*)
 end if
!21. occ
 if(occopt==2)then
  do iline=1,(bantot+2)/3
   read (unddb,*)
  end do
 else
  write(6,*)' inprep8 : nband(1)=',nband(1)
  do iline=1,(nband(1)+2)/3
   read (unddb,'(a80)')rdstring
   write(6,*)trim(rdstring)
  end do
 end if
!22. rprim
 do iline=1,3
  read (unddb,*)
 end do
!23. sciss
 read (unddb,*)
!24. spinat
 if(ddbvrs==vrsio8)then
  do iline=1,natom
   read (unddb,*)
  end do
 end if
!25. symafm
 if(ddbvrs==vrsio8)then
  do iline=1,(nsym+11)/12
   read (unddb,*)
  end do
 end if
!26. symrel
 do iline=1,nsym
  read (unddb,*)
 end do
!27old. xred
 if(ddbvrs/=vrsio8)then
  do iline=1,natom
   read (unddb,*)
  end do
 end if
!27. tnons
 do iline=1,nsym
  read (unddb,*)
 end do
!28. tolwfr
 if(ddbvrs==vrsio8)then
  read (unddb,*)
 end if
!29. tphysel
 if(ddbvrs==vrsio8)then
  read (unddb,*)
 end if
!30. tsmear
 if(ddbvrs==vrsio8)then
  read (unddb,*)
 end if
!31. type
 do iline=1,(natom+11)/12
  read (unddb,*)
 end do
!31old. tolwfr
 if(ddbvrs/=vrsio8)then
  read (unddb,*)
 end if
!32. wtk
 do iline=1,(nkpt+2)/3
  read (unddb,*)
 end do
!33. xred
 if(ddbvrs==vrsio8)then
  do iline=1,natom
   read (unddb,*)
  end do
 end if
!34. znucl
 if(ddbvrs==vrsio8)then
  do iline=1,(ntypat+2)/3
   read (unddb,*)
  end do
 end if
!35. zion
 do iline=1,(ntypat+2)/3
  read (unddb,*)
 end do

 read (unddb,*)

!Now, take care of the pseudopotentials
 read(unddb, '(a12)' )string

 if(string=='  Descriptio')then

  read (unddb,*)
  read (unddb, '(10x,i3,14x,i3,11x,i3)', iostat=ios )dimekb,lmnmax,usepaw
  if(ios/=0)then
   backspace(unddb)
   read (unddb, '(10x,i3,14x,i3)')dimekb,lmnmax
   usepaw=0
  end if
  do itypat=1,ntypat
   read(unddb, '(13x,i4,9x,i3,8x,i4)' )itypat0,pspso0,nekb
   read(unddb,*)
   do iekb=1,nekb
    do jekb=1,nekb,4
     read(unddb,*)
    end do
   end do
  end do

 else if(string==' Description')then

  read (unddb, '(10x,i3,10x,i3)' )mproj,mpsang
  dimekb=mproj*mpsang
  usepaw=0
  do itypat=1,ntypat
   read (unddb,*)
!  For f-electrons, one more line has been written
   do iproj=1,mproj*max(1,(mpsang+2)/3)
    read (unddb,*)
   end do
  end do

 else if(string==' No informat')then

  dimekb=0
  usepaw=0

 else
  write(message, '(a,a,a,a,a,a)' )&
&  ' inprep8 : BUG -',ch10,&
&  '  Error when reading the psp information',ch10,&
&  '  String=',string
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Check the value of usepaw
 if (usepaw==1) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' psddb8: BUG -',ch10,&
&  '  Paw not yet allowed !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Now, the number of blocks
 read(unddb,*)
 read(unddb,*)
 read(unddb, '(24x,i4)' )nblok

!Now, the type of each blok, in turn
 mblktyp=1
 if(nblok>=1)then
  do iblok=1,nblok

   read(unddb,*)
   read(unddb, '(a32,12x,i8)' )blkname,nelmts
   if(blkname==' 2rd derivatives (non-stat.)  - ')then
    blktyp=1
   else if(blkname==' 2rd derivatives (stationary) - ')then
    blktyp=2
   else if(blkname==' 3rd derivatives              - ')then
    blktyp=3
   else if(blkname==' Total energy                 - ')then
    blktyp=0
   else if(blkname==' 1st derivatives              - ')then
    blktyp=4
   else if(blkname==' 2rd eigenvalue derivatives   - ')then
    blktyp=5
   else
    write(message, '(a,a,a,a,a,a,a,a)' )&
&    ' inprep8 : ERROR -',ch10,&
&    '  The following string appears in the DDB in place of',&
&    ' the block type description :',ch10,name,ch10,&
&    '  Action : check your DDB.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   if(blktyp==1.or.blktyp==2)then
!   Read the phonon wavevector
    read(unddb,*)
   else if(blktyp==3)then
!   Read the perturbation wavevectors
    read(unddb,*)
    read(unddb,*)
    read(unddb,*)
    mblktyp=3
   else if(blktyp==5)then
    read(unddb,*)
   end if

!  Read every element
   if(blktyp==5)then
    do ikpt=1,nkpt
     read(unddb,*)
     do iband=1,nband(ikpt)
      read(unddb,*)
      do ii=1,nelmts
       read(unddb,*)
      end do
     end do
    end do
   else
    do ii=1,nelmts
     read(unddb,*)
    end do
   end if

  end do
 end if

 deallocate(nband)

!Close the DDB
 close(unddb)

end subroutine inprep8
!!***
