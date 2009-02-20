!{\src2tex{textfont=tt}}
!!****f* ABINIT/ioddb8
!!
!! NAME
!! ioddb8
!!
!! FUNCTION
!! Open Derivative DataBase, then
!! reads or write Derivative DataBase preliminary information.
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
!! choice=(1=> input only),(2=> output only)
!! character(len=fnlen) dscrpt:string that describe the output database (not needed
!!   if choice=1)
!! character(len=fnlen) filnam: name of input or output file
!! matom=maximum number of atoms
!! mband=maximum number of bands
!! mkpt=maximum number of special points
!! msym=maximum number of symetries
!! mtypat=maximum number of atom types
!! unddb=unit number for input or output
!! vrsddb=6 digit integer giving date, in form yymmdd for month=mm(1-12),
!!  day=dd(1-31), and year=yy(90-99 for 1990 to 1999,00-89 for 2000 to 2089),
!!  of current DDB version.  Have to agree with vrsinddb of this routine.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! The following arguments are outputs or inputs, according to the
!! option choice, and are either read or written to the Derivative
!! DataBase:
!! acell(3)=length scales of primitive translations (bohr)
!! amu(mtypat)=mass of the atoms (atomic mass unit)
!! dilatmx=the maximal dilatation factor
!! ecut=kinetic energy planewave cutoff (hartree)
!! ecutsm=smearing energy for plane wave kinetic energy (Ha)
!! intxc=control xc quadrature
!! iscf=parameter controlling scf or non-scf choice
!! ixc=exchange-correlation choice parameter
!! kpt(3,mkpt)=k point set (reduced coordinates)
!! kptnrm=normalisation of k points
!! natom=number of atoms in the unit cell
!! nband(mkpt)=number of bands at each k point, for each polarization
!! ngfft(18)=contain all needed information about 3D FFT,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkpt=number of k points
!! nspden=number of spin-density components
!! nspinor=number of spinorial components of the wavefunctions
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetry elements in space group
!! ntypat=number of atom types
!! occ(mband*mkpt)=occupation number for each band and k
!! occopt=option for occupancies
!! rprim(3,3)=dimensionless primitive translations in real space
!! sciss=scissor shift (Ha)
!! spinat(3,matom)=initial spin of each atom, in unit of hbar/2
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym)=symmetry operations in real space
!! tnons(3,msym)=nonsymmorphic translations for symmetry operations
!! tolwfr=tolerance on largest wf residual
!! tphysel="physical" electronic temperature with FD occupations
!! tsmear=smearing width (or temperature) in Hartree
!! typat(matom)=type of each atom
!! wtk(mkpt)=weight assigned to each k point
!! xred(3,matom)=reduced atomic coordinates
!! zion(mtypat)=valence charge of each type of atom
!! znucl(mtypat)=atomic number of atom type
!!
!! TODO
!!
!! PARENTS
!!      gstate,loper3,mrgddb,nonlinear,rdddb9,respfn
!!
!! CHILDREN
!!      chknm8,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ioddb8 (choice,dscrpt,filnam,matom,mband,&
&  mkpt,msym,mtypat,unddb,vrsddb,&
&  acell,amu,dilatmx,ecut,ecutsm,intxc,iscf,ixc,kpt,kptnrm,&
&  natom,nband,ngfft,nkpt,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  rprim,sciss,spinat,symafm,symrel,tnons,tolwfr,tphysel,tsmear,typat,&
&  wtk,xred,zion,znucl)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,matom,mband,mkpt,msym,mtypat,unddb,vrsddb
 integer,intent(inout) :: intxc,iscf,ixc,natom,nkpt,nspden,nspinor,nsppol,nsym
 integer,intent(inout) :: ntypat,occopt
 real(dp),intent(inout) :: dilatmx,ecut,ecutsm,kptnrm,sciss,tolwfr,tphysel
 real(dp),intent(inout) :: tsmear
 character(len=fnlen),intent(in) :: dscrpt,filnam
!arrays
 integer,intent(inout) :: nband(mkpt),ngfft(18),symafm(msym),symrel(3,3,msym)
 integer,intent(inout) :: typat(matom)
 real(dp),intent(inout) :: acell(3),amu(mtypat),kpt(3,mkpt),occ(mband*mkpt)
 real(dp),intent(inout) :: rprim(3,3),spinat(3,matom),tnons(3,msym),wtk(mkpt)
 real(dp),intent(inout) :: xred(3,matom),zion(mtypat),znucl(mtypat)

!Local variables -------------------------
!Set routine version number here:
!scalars
 integer,parameter :: vrsio8=010929,vrsio8_old=990527
 integer :: bantot,ddbvrs,iband,ii,ij,ikpt,iline,im
 logical :: testn,testv
 character(len=500) :: message
 character(len=6) :: name_old
!arrays
 character(len=9) :: name(9)

! *********************************************************************

!DEBUG
!write(6,*)' ioddb8 : enter'
!write(6,*)' ioddb8 : mtypat=',mtypat
!ENDDEBUG

!Check ioddb8 version number (vrsio8) against mkddb version number
!(vrsddb)
 if (vrsio8/=vrsddb) then
  write(message, '(a,a,a,i10,a,a,i10,a)' )&
&  ' ioddb8: BUG -',ch10,&
&  '  The input/output DDB version number=',vrsio8,ch10,&
&  '  is not equal to the DDB version number=',vrsddb,'.'
  call wrtout(6,message,'COLL')
! call leave_new('COLL')
 end if

!Check the value of choice :
 if (choice<=0.or.choice>=3) then
  write(message, '(a,a,a,a,i4,a)' )&
&  ' ioddb8: BUG -',ch10,&
&  '  The permitted values for choice are 1 or 2.',&
&  '  The calling routine asks choice =',choice,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (choice==1) then

! Open the input derivative database.
  open (unit=unddb,file=filnam,status='old',form='formatted')

! Check the compatibility of the input DDB with the DDB code
  read (unddb,*)
  read (unddb,*)
  read (unddb, '(20x,i10)' )ddbvrs
  if(ddbvrs/=vrsio8 .and. ddbvrs/=vrsio8_old)then
   write(message, '(a,a,a,a,a)' )&
&   ' ioddb8 : BUG - ',ch10,&
&   '  The ioddb8 DDB version number do not agree',ch10,&
&   '  with the calling code DDB version number.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Read the 4 n-integers, also testing the names of data,
! and checking that their value is acceptable.
! This is important to insure that any array has a sufficient
! dimension.
  read (unddb,*)
  read (unddb,*)
  read (unddb,*)
  testn=.true.
  testv=.true.
! 1. natom
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),natom
  else
   read (unddb, '(1x,a6,i10)' )name_old,natom ; name(1)='   '//name_old
  end if
  if(name(1)/='    natom')testn=.false.
  if(natom<=0.or.natom>matom)testv=.false.
! 2. nkpt
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(2),nkpt
  else
   read (unddb, '(1x,a6,i10)' )name_old,nkpt ; name(2)='   '//name_old
  end if
  if(name(2)/='     nkpt')testn=.false.
  if(nkpt <=0.or.nkpt >mkpt )testv=.false.
! 3. nsppol
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(3),nsppol
  else
   read (unddb, '(1x,a6,i10)' )name_old,nsppol ; name(3)='   '//name_old
  end if
  if(name(3)/='   nsppol')testn=.false.
  if(nsppol<=0.or.nsppol>2)testv=.false.
! 4. nsym
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(4),nsym
  else
   read (unddb, '(1x,a6,i10)' )name_old,nsym ; name(4)='   '//name_old
  end if
  if(name(4)/='     nsym')testn=.false.
  if(nsym <=0.or.nsym >msym )testv=.false.
! 5. ntypat
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(5),ntypat
  else
   read (unddb, '(1x,a6,i10)' )name_old,ntypat ; name(5)='   '//name_old
  end if
  if(name(5)/='   ntypat' .and. name(5)/='    ntype')testn=.false.
  if(ntypat<=0.or.ntypat>mtypat)testv=.false.
! 6. occopt
! Before reading nband, the last parameters that define
! the dimension of some array, need to know what is their
! representation, given by occopt
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(6),occopt
  else
   read (unddb, '(1x,a6,i10)' )name_old,occopt ; name(6)='   '//name_old
  end if
  if(name(6)/='   occopt')testn=.false.
  if(occopt<0.or.occopt>7)testv=.false.
! Message if the names or values are not right
  if (.not.testn.or..not.testv) then
   write(message, '(a,a,a)' )' ioddb8 : An error has been found in one',ch10,&
&   ' of the positive n-integers contained in the DDB : '
   call wrtout(6,message,'COLL')
   write(message, '(a)' )&
&   '               Expected                      Found     '
   call wrtout(6,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    natom , lower than',matom+1,'    ',trim(name(1)),' =',natom
   call wrtout(6,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nkpt  , lower than',mkpt+1 ,'    ',trim(name(2)),' =',nkpt
   call wrtout(6,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nsppol, lower than',3      ,'    ',trim(name(3)),' =',nsppol
   call wrtout(6,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    nsym  , lower than',msym+1 ,'    ',trim(name(4)),' =',nsym
   call wrtout(6,message,'COLL')
   write(message, '(a,i10,a,a,a,i10)' )&
&   '    ntypat, lower than',mtypat+1,'    ',trim(name(5)),' =',ntypat
   call wrtout(6,message,'COLL')
   write(message, '(a,a,a,i10)' )&
&   '    occopt,  between 0 and 7        ',trim(name(6)),' =',occopt
   call wrtout(6,message,'COLL')
   write(message, '(a,a,a)' )&
&   ' ioddb8 : ERROR -',ch10,&
&   '  See the error message above.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! One more set of parameters define the dimensions of the
! array : nband. Morever, it depends on occopt and
! nkpt, and has to be
! tested after the test on nkpt is performed.
! 7. nband
  if(occopt==2)then
   im=12
   do iline=1,(nkpt+11)/12
    if(iline==(nkpt+11)/12)im=nkpt-12*(iline-1)
    if(ddbvrs==vrsio8)then
     read (unddb, '(1x,a9,5x,12i5)' )name(1),&
&     (nband((iline-1)*12+ii),ii=1,im)
    else
     read (unddb, '(1x,a6,5x,12i5)' )name_old,&
&     (nband((iline-1)*12+ii),ii=1,im) ; name(1)='   '//name_old
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

! check all nband values, and sum them
  bantot=0
  do ikpt=1,nkpt
   if(nband(ikpt)<0)then
    write(message, '(3a,i4,a,i4,3a)' )&
&    ' ioddb8 : ERROR -',ch10,&
&    '  For ikpt = ',ikpt,'  nband = ',nband(ikpt),' is negative.',ch10,&
&    '  Action : correct your DDB.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   else if(nband(ikpt)>mband)then
    write(message, '(3a,i4,a,i4,a,a,i4,3a)' )&
&    ' ioddb8 : ERROR -',ch10,&
&    ' For ikpt = ',ikpt,', nband = ',nband(ikpt),ch10,&
&    ' is larger than mband = ',mband,'.',ch10,&
&    ' Action : recompile the calling code with a larger mband.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   bantot=bantot+nband(ikpt)
  end do

! Read the rest of variables, with check of the names
! 8. acell
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,3d22.14)' )name(1),acell
  else
   read (unddb, '(1x,a6,3d22.14)' )name_old,acell ; name(1)='   '//name_old
  end if
  call chknm8(name(1),'    acell')
! 9. amu
  im=3
  do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (amu((iline-1)*3+ii),ii=1,im)
   else
    read (unddb, '(1x,a6,3d22.14)' )name_old,&
&    (amu((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
    call chknm8(name(1),'      amu')
   else
    call chknm8(name(1),'         ')
   end if
  end do
! 10. dilatmx
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,d22.14)' )name(1),dilatmx
   call chknm8(name(1),'  dilatmx')
  else
   dilatmx=one
  end if
! 11. ecut
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,d22.14)' )name(1),ecut
  else
   read (unddb, '(1x,a6,d22.14)' )name_old,ecut ; name(1)='   '//name_old
  end if
  call chknm8(name(1),'     ecut')
! 12. ecutsm
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,d22.14)' )name(1),ecutsm
   call chknm8(name(1),'   ecutsm')
  else
   ecutsm=zero
  end if
! 13. intxc
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),intxc
   call chknm8(name(1),'    intxc')
  else
   intxc=1
  end if
! 14. iscf
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),iscf
  else
   read (unddb, '(1x,a6,i10)' )name_old,iscf ; name(1)='   '//name_old
  end if
  call chknm8(name(1),'     iscf')
! 15. ixc
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),ixc
  else
   read (unddb, '(1x,a6,i10)' )name_old,ixc ; name(1)='   '//name_old
  end if
  call chknm8(name(1),'      ixc')
! 16. kpt
  do iline=1,nkpt
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (kpt(ii,iline),ii=1,3)
   else
    read (unddb, '(1x,a6,3d22.14)' )name_old,&
&    (kpt(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
    call chknm8(name(1),'      kpt')
   else
    call chknm8(name(1),'         ')
   end if
  end do
! 17. kptnrm
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,d22.14)' )name(1),kptnrm
  else
   read (unddb, '(1x,a6,d22.14)' )name_old,kptnrm ; name(1)='   '//name_old
  end if
  call chknm8(name(1),'   kptnrm')
! 18. ngfft
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,5x,3i5)' )name(1),ngfft(1:3)
  else
   read (unddb, '(1x,a6,5x,3i5)' )name_old,ngfft(1:3) ; name(1)='   '//name_old
  end if
! For the time being, do not check the validity of the name,
! in order to accept both ng and ngfft
! 19. nspden
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),nspden
   call chknm8(name(1),'   nspden')
  else
   nspden=0
  end if
! 20. nspinor
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,i10)' )name(1),nspinor
   call chknm8(name(1),'  nspinor')
  else
   nspinor=0
  end if
! 21. occ
  if(occopt==2)then
   im=3
   do iline=1,(bantot+2)/3
    if(iline==(bantot+2)/3)im=bantot-3*(iline-1)
    if(ddbvrs==vrsio8)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (occ((iline-1)*3+ii),ii=1,im)
    else
     read (unddb, '(1x,a6,3d22.14)' )name_old,&
&     (occ((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
    end if
    if (iline==1) then
     call chknm8(name(1),'      occ')
    else
     call chknm8(name(1),'         ')
    end if
   end do
  else
   im=3
   do iline=1,(nband(1)+2)/3
    if(iline==(nband(1)+2)/3)im=nband(1)-3*(iline-1)
    if(ddbvrs==vrsio8)then
     read (unddb, '(1x,a9,3d22.14)' )name(1),&
&     (occ((iline-1)*3+ii),ii=1,im)
    else
     read (unddb, '(1x,a6,3d22.14)' )name_old,&
&     (occ((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
    end if
    if (iline==1) then
     call chknm8(name(1),'      occ')
    else
     call chknm8(name(1),'         ')
    end if
   end do
   if(nkpt>1)then
    do ikpt=2,nkpt
     do iband=1,nband(1)
      occ(iband+nband(1)*(ikpt-1))=occ(iband)
     end do
    end do
   end if
  end if
! 22. rprim
  do iline=1,3
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (rprim(ii,iline),ii=1,3)
   else
    read (unddb, '(1x,a6,3d22.14)' )name_old,&
&    (rprim(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
    call chknm8(name(1),'    rprim')
   else
    call chknm8(name(1),'         ')
   end if
  end do
! 23. sciss
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,d22.14)' )name(1),sciss
  else
   read (unddb, '(1x,a6,d22.14)' )name_old,sciss ; name(1)='   '//name_old
  end if
  call chknm8(name(1),'    sciss')
! 24. spinat
  if(ddbvrs==vrsio8)then
   do iline=1,natom
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (spinat(ii,iline),ii=1,3)
    if (iline==1) then
     call chknm8(name(1),'   spinat')
    else
     call chknm8(name(1),'         ')
    end if
   end do
  else
!  spinat is set to zero by default in mrgddb.f
!  spinat(:,1:natom)=zero
  end if
! 25. symafm
  if(ddbvrs==vrsio8)then
   im=12
   do iline=1,(nsym+11)/12
    if(iline==(nsym+11)/12)im=nsym-12*(iline-1)
    read (unddb, '(1x,a9,5x,12i5)' )name(1),&
&    (symafm((iline-1)*12+ii),ii=1,im)
    if (iline==1) then
     call chknm8(name(1),'   symafm')
    else
     call chknm8(name(1),'         ')
    end if
   end do
  else
!  symafm is set to 1 by default in mrgddb.f
!  symafm(1:nsym)=1
  end if
! 26. symrel
  do iline=1,nsym
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,5x,9i5)' )name(1),&
&    ((symrel(ii,ij,iline),ii=1,3),ij=1,3)
   else
    read (unddb, '(1x,a6,5x,9i5)' )name_old,&
&    ((symrel(ii,ij,iline),ii=1,3),ij=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
    call chknm8(name(1),'   symrel')
   else
    call chknm8(name(1),'         ')
   end if
  end do
! 27old. xred
  if(ddbvrs/=vrsio8)then
   do iline=1,natom
    read (unddb, '(1x,a6,3d22.14)' )name(1),&
&    (xred(ii,iline),ii=1,3)
   end do
!  No check of name, to allow the old tn
  end if
! 27. tnons
  do iline=1,nsym
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (tnons(ii,iline),ii=1,3)
   else
    read (unddb, '(1x,a6,3d22.14)' )name_old,&
&    (tnons(ii,iline),ii=1,3) ; name(1)='   '//name_old
   end if
   if (iline==1) then
    call chknm8(name(1),'    tnons')
   else
    call chknm8(name(1),'         ')
   end if
  end do
! 28. tolwfr
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tolwfr
  end if
! Do not check the name, in order to allow both tolwfr and wftol
! 29. tphysel
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tphysel
   call chknm8(name(1),'  tphysel')
  else
   tphysel=zero
  end if
! 30. tsmear
  if(ddbvrs==vrsio8)then
   read (unddb, '(1x,a9,d22.14)' )name(1),tsmear
   call chknm8(name(1),'   tsmear')
  else
   tsmear=zero
  end if
! 31. typat
  im=12
  do iline=1,(natom+11)/12
   if(iline==(natom+11)/12)im=natom-12*(iline-1)
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,5x,12i5)' )name(1),&
&    (typat((iline-1)*12+ii),ii=1,im)
   else
    read (unddb, '(1x,a6,5x,12i5)' )name_old,&
&    (typat((iline-1)*12+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
!   Both type and typat are allowed => no check
!   call chknm8(name(1),'    typat')
   else
    call chknm8(name(1),'         ')
   end if
  end do
! 31old. tolwfr
  if(ddbvrs/=vrsio8)then
   read (unddb, '(1x,a6,d22.14)' )name(1),tolwfr
  end if
! Do not check the name, in order to allow both tolwfr and wftol
! 32. wtk
  im=3
  do iline=1,(nkpt+2)/3
   if(iline==(nkpt+2)/3)im=nkpt-3*(iline-1)
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (wtk((iline-1)*3+ii),ii=1,im)
   else
    read (unddb, '(1x,a6,3d22.14)' )name_old,&
&    (wtk((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
    call chknm8(name(1),'      wtk')
   else
    call chknm8(name(1),'         ')
   end if
  end do
! 33. xred
  if(ddbvrs==vrsio8)then
   do iline=1,natom
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (xred(ii,iline),ii=1,3)
    if (iline==1) then
     call chknm8(name(1),'     xred')
    else
     call chknm8(name(1),'         ')
    end if
   end do
  end if
! 34. znucl
  if(ddbvrs==vrsio8)then
   im=3
   do iline=1,(ntypat+2)/3
    if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (znucl((iline-1)*3+ii),ii=1,im)
    if (iline==1) then
     call chknm8(name(1),'    znucl')
    else
     call chknm8(name(1),'         ')
    end if
   end do
  else
!  znucl is set to zero by default in mrgddb.f
!  znucl(:)=zero
  end if
! 35. zion
  im=3
  do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   if(ddbvrs==vrsio8)then
    read (unddb, '(1x,a9,3d22.14)' )name(1),&
&    (zion((iline-1)*3+ii),ii=1,im)
   else
    read (unddb, '(1x,a6,3d22.14)' )name_old,&
&    (zion((iline-1)*3+ii),ii=1,im) ; name(1)='   '//name_old
   end if
   if (iline==1) then
!   Do not check the names, to allow both zion and znucl - the latter for 990527 format
!   call chknm8(name(1),'     zion')
   else
    call chknm8(name(1),'         ')
   end if
  end do

! *****************

 else if (choice==2) then
! 
! Open the output derivative database.
! (version 2.1. : changed because of a bug in a Perl script
! should set up a name checking procedure, with change of name
! like for the output file)
! open (unit=unddb,file=filnam,status='new',form='formatted')
  open (unit=unddb,file=filnam,status='unknown',form='formatted')

! Write the heading
  write(unddb, '(/,a,/,a,i10,/,/,a,a,/)' ) &
&  ' **** DERIVATIVE DATABASE ****    ',&
&  '+DDB, Version number',vrsddb,' ',dscrpt

! Write the descriptive data
! 1. natom
  write(unddb, '(1x,a9,i10)' )'    natom',natom
! 2. nkpt
  write(unddb, '(1x,a9,i10)' )'     nkpt',nkpt
! 3. nsppol
  write(unddb, '(1x,a9,i10)' )'   nsppol',nsppol
! 4. nsym
  write(unddb, '(1x,a9,i10)' )'     nsym',nsym
! 5. ntypat
  write(unddb, '(1x,a9,i10)' )'   ntypat',ntypat
! 6. occopt
  write(unddb, '(1x,a9,i10)' )'   occopt',occopt
! 7. nband
  if(occopt==2)then
   im=12
   name(1)='    nband'
   do iline=1,(nkpt+11)/12
    if(iline==(nkpt+11)/12)im=nkpt-12*(iline-1)
    write(unddb, '(1x,a9,5x,12i5)' )name(1),&
&    (nband((iline-1)*12+ii),ii=1,im)
    name(1)='         '
   end do
   bantot=0
   do ikpt=1,nkpt
    bantot=bantot+nband(ikpt)
   end do
  else
   write(unddb, '(1x,a9,i10)' )'    nband',nband(1)
   bantot=nkpt*nband(1)
  end if

! 8. acell
  write(unddb, '(1x,a9,3d22.14)' )'    acell',acell
! 9. amu
  im=3
  name(1)='      amu'
  do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write (unddb, '(1x,a9,3d22.14)' )name(1),&
&   (amu((iline-1)*3+ii),ii=1,im)
   name(1)='         '
  end do
! 10. dilatmx
  write(unddb, '(1x,a9,d22.14)' )'  dilatmx',dilatmx
! 11. ecut
  write(unddb, '(1x,a9,d22.14)' )'     ecut',ecut
! 12. ecutsm
  write(unddb, '(1x,a9,d22.14)' )'   ecutsm',ecutsm
! 13. intxc
  write(unddb, '(1x,a9,i10)' )'    intxc',intxc
! 14. iscf
  write(unddb, '(1x,a9,i10)' )'     iscf',iscf
! 15. ixc
  write(unddb, '(1x,a9,i10)' )'      ixc',ixc
! 16. kpt
  name(1)='      kpt'
  do iline=1,nkpt
   write (unddb, '(1x,a9,3d22.14)' )name(1),&
&   (kpt(ii,iline),ii=1,3)
   name(1)='      '
  end do
! 17. kptnrm
  write(unddb, '(1x,a9,d22.14)' )'   kptnrm',kptnrm
! 18. ngfft
  write(unddb, '(1x,a9,5x,3i5)' )'    ngfft',ngfft(1:3)
! 19. nspden
  write(unddb, '(1x,a9,i10)' )'   nspden',nspden
! 20. nspinor
  write(unddb, '(1x,a9,i10)' )'  nspinor',nspinor
! 21. occ
  if(occopt==2)then
   im=3
   name(1)='      occ'
   do iline=1,(bantot+2)/3
    if(iline==(bantot+2)/3)im=bantot-3*(iline-1)
    write(unddb, '(1x,a9,3d22.14)' )name(1),&
&    (occ((iline-1)*3+ii),ii=1,im)
    name(1)='         '
   end do
  else
   im=3
   name(1)='      occ'
   do iline=1,(nband(1)+2)/3
    if(iline==(nband(1)+2)/3)im=nband(1)-3*(iline-1)
    write(unddb, '(1x,a9,3d22.14)' )name(1),&
&    (occ((iline-1)*3+ii),ii=1,im)
    name(1)='         '
   end do
  end if
! 22. rprim
  name(1)='    rprim'
  do iline=1,3
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (rprim(ii,iline),ii=1,3)
   name(1)='      '
  end do
! 23. sciss
  write(unddb, '(1x,a9,d22.14)' )'    sciss',sciss
! 24. spinat
  name(1)='   spinat'
  do iline=1,natom
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (spinat(ii,iline),ii=1,3)
   name(1)='         '
  end do
! 25. symafm
! DEBUG
! write(6,*)' ioddb8 : symafm before write=',symafm(:)
! ENDDEBUG
  im=12
  name(1)='   symafm'
  do iline=1,(nsym+11)/12
   if(iline==(nsym+11)/12)im=nsym-12*(iline-1)
   write(unddb, '(1x,a9,5x,12i5)' )name(1),&
&   (symafm((iline-1)*12+ii),ii=1,im)
   name(1)='         '
  end do
! 26. symrel
  name(1)='   symrel'
  do iline=1,nsym
   write(unddb, '(1x,a9,5x,9i5)' )name(1),&
&   ((symrel(ii,ij,iline),ii=1,3),ij=1,3)
   name(1)='         '
  end do
! 27. tnons
  name(1)='    tnons'
  do iline=1,nsym
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (tnons(ii,iline),ii=1,3)
   name(1)='         '
  end do
! 28. tolwfr
  write(unddb, '(1x,a9,d22.14)' )'   tolwfr',tolwfr
! 29. tphysel
  write(unddb, '(1x,a9,d22.14)' )'  tphysel',tphysel
! 30. tsmear
  write(unddb, '(1x,a9,d22.14)' )'   tsmear',tsmear
! 31. typat
  im=12
  name(1)='    typat'
  do iline=1,(natom+11)/12
   if(iline==(natom+11)/12)im=natom-12*(iline-1)
   write(unddb, '(1x,a9,5x,12i5)' )name(1),&
&   (typat((iline-1)*12+ii),ii=1,im)
   name(1)='         '
  end do
! 32. wtk
  name(1)='      wtk'
  im=3
  do iline=1,(nkpt+2)/3
   if(iline==(nkpt+2)/3)im=nkpt-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (wtk((iline-1)*3+ii),ii=1,im)
   name(1)='         '
  end do
! 33. xred
  name(1)='     xred'
  do iline=1,natom
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (xred(ii,iline),ii=1,3)
   name(1)='         '
  end do
! 34. znucl
  name(1)='    znucl'
  im=3
  do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (znucl((iline-1)*3+ii),ii=1,im)
   name(1)='         '
  end do
! 35. zion
  name(1)='     zion'
  im=3
  do iline=1,(ntypat+2)/3
   if(iline==(ntypat+2)/3)im=ntypat-3*(iline-1)
   write(unddb, '(1x,a9,3d22.14)' )name(1),&
&   (zion((iline-1)*3+ii),ii=1,im)
   name(1)='         '
  end do

 end if

end subroutine ioddb8
!!***
