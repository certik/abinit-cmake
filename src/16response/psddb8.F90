!{\src2tex{textfont=tt}}
!!****f* ABINIT/psddb8
!!
!! NAME
!! psddb8
!!
!! FUNCTION
!! Take care of the i/o of pseudopotentials for the
!! Derivative DataBase, and also the number of data blocks.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  choice=(1 => read), (2=> write)
!!  dimekb=dimension of ekb (for the time being, only for norm-
!!                           conserving psps)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  nunit=unit number for the Derivative DataBase.
!!  ntypat=number of atom types
!!  pspso(ntypat)=For each type of psp, 1 if no spin-orbit component is taken
!!     into account, 2 if a spin-orbit component is used
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  vrsddb=Derivative Database version, for check of compatibility.
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                                 or i=lmn (if useylm=1)
!!  ekb(dimekb,ntypat)= (for the time being, norm-conserving psps only)
!!    (Real) Kleinman-Bylander energies (hartree)
!!    Presently the only characteristics of the psp
!!  fullinit=0 if the ekb are not available, at input as well as at output
!!  nblok=number of blocks
!!
!! NOTES
!! Only executed by one processor
!! This routine is in a transient state, since the
!! treatment of pseudopotentials is changing (introduction
!! of the paw method)
!!
!! PARENTS
!!      gstate,mrgddb,nonlinear,rdddb9,respfn
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine psddb8 (choice,dimekb,ekb,fullinit,indlmn,lmnmax,&
&      lnmax,nblok,ntypat,nunit,pspso,usepaw,useylm,vrsddb)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,dimekb,lmnmax,lnmax,ntypat,nunit,usepaw,useylm
 integer,intent(in) :: vrsddb
 integer,intent(inout) :: fullinit,nblok
!arrays
 integer,intent(in) :: pspso(ntypat)
 integer,intent(inout) :: indlmn(6,lmnmax,ntypat)
 real(dp),intent(inout) :: ekb(dimekb,ntypat)

!Local variables -------------------------
!Set the version number
!scalars
 integer,parameter :: vrsio8=010929,vrsio8_old=990527
 integer :: dimekb0,iekb,ii,ij,ilmn,iln,iln0,im,ios,iproj,iproj0,itypat,itypat0
 integer :: jekb,jj,jlmn,jln,kk,lmnmax0,lpsang,nekb,nproj,npsang,pspso0,usepaw0
 integer :: vrspsp8
 character(len=12) :: string
 character(len=500) :: message
!arrays
 real(dp) :: ekb0(dimekb,dimekb)

! *********************************************************************


!Check psddb8 version number (vrsio8) against DDB version number
!(vrsddb)
 if (vrsio8/=vrsddb) then
  write(message, '(a,a,a,a,i10,a,a,i10,a)' ) ch10,&
&  ' psddb8 : BUG -',ch10,&
&  '  the psddb8 DDB version number=',vrsio8,ch10,&
&  '    is not equal to the calling code DDB version number=',vrsddb,'.'
  call wrtout(6,message,'COLL')
! call leave_new('COLL')
 end if

!Check the value of choice
 if (choice<=0.or.choice>=3) then
  write(message, '(a,a,a,a,a,a,i10,a)' ) ch10,&
&  ' psddb8: BUG -',ch10,&
&  '  The permitted values for choice are 1 or 2.',ch10,&
&  '  The calling routine asks ',choice,'.'
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

!Check the compatibility of lnmax and dimekb
 if (lnmax>dimekb) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' psddb8: BUG -',ch10,&
&  '  For the time being lnmax must <= dimekb'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (choice==1) then
! Read option

  read(nunit,*)
  read(nunit, '(a12)' )string

  fullinit=1 ; indlmn(:,:,:)=0
  if (string=='  Descriptio')then
!  This is the new format
   read(nunit, '(32x,i6)' )vrspsp8
   usepaw0=0
   read(nunit, '(10x,i3,14x,i3,11x,i3)', iostat=ios)dimekb0,lmnmax0,usepaw0
   if(ios/=0)then
    backspace(nunit)
    read (nunit, '(10x,i3,14x,i3)' )dimekb0,lmnmax0
    usepaw0=0
   end if
   ekb(:,:)=zero
   do itypat=1,ntypat

    read(nunit, '(13x,i4,9x,i3,8x,i4)' )itypat0,pspso0,nekb

!   Check the compatibility with the main code dimensioning
    if(nekb>dimekb)then
     write(message, '(a,a,a,a,i8,a,a,a,i3,a)' ) ch10,&
&     ' psddb8 : BUG -',ch10,&
&     '  ',nekb,' components of ekb are announced',ch10,&
&     '  but dimekb=',dimekb,'.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
    if(usepaw0/=usepaw)then
     write(message, '(a,a,a,a)' ) ch10,&
&     ' psddb8: BUG -',ch10,&
&     '  Read usepaw0 must equal usepaw !'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if

    read(nunit,*)
    ilmn=0;iproj0=0
    do iekb=1,nekb
     read(nunit, '(3i6,3x,8d15.7)' ) iln,lpsang,iproj,&
&     (ekb0(ii,iekb),ii=1,min(nekb,4))
     if(nekb>4)then
      do jekb=5,nekb,4
       read(nunit, '(21x,8d15.7)' )&
&       (ekb0(ii,iekb),ii=jekb,min(nekb,jekb+3))
      end do
     end if
     if (lpsang==0.and.iproj>iproj0) iproj0=iproj
     if (useylm==1) then
      do im=-lpsang,lpsang
       ilmn=ilmn+1
       indlmn(1,ilmn,itypat)=lpsang
       indlmn(2,ilmn,itypat)=im
       indlmn(3,ilmn,itypat)=iproj
       indlmn(4,ilmn,itypat)=lpsang**2+lpsang+1+im
       indlmn(5,ilmn,itypat)=iln
       indlmn(6,ilmn,itypat)=1
       if (pspso0/=1.and.iln>(nekb-iproj0)/2) indlmn(6,ilmn,itypat)=2
      end do
     else
      ilmn=ilmn+1
      indlmn(1,ilmn,itypat)=lpsang
      indlmn(2,ilmn,itypat)=lpsang
      indlmn(3,ilmn,itypat)=iproj
      indlmn(4,ilmn,itypat)=lpsang**2+lpsang+1
      indlmn(5,ilmn,itypat)=iln
      indlmn(6,ilmn,itypat)=1
      if (pspso0/=1.and.iln>(nekb-iproj0)/2) indlmn(6,ilmn,itypat)=2
     end if
!    For the time being, only diagonal ekb are treated in abinit v3
     ekb(iekb,itypat)=ekb0(iekb,iekb)
!    For non-diagonal ekb, one could use:
!    do jekb=iekb to nekb
!    ekb(jekb+iekb*(iekb-1)/2,itypat)=ekb0(jekb,iekb)
!    end do
    end do
   end do

  else if (string==' Description')then
   read (nunit, '(10x,i3,10x,i3)' )nproj,npsang
   nekb=nproj*npsang

!  Check the compatibility with the main code dimensioning
   if(nekb>dimekb)then
    write(message, '(a,a,a,a,i8,a,a,a,i3,a)' ) ch10,&
&    ' psddb8 : BUG -',ch10,&
&    '  ',nekb,' components of ekb are announced',ch10,&
&    '  but the maximum is dimekb=',dimekb,'.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if(useylm/=0)then
    write(message, '(a,a,a,a)' ) ch10,&
&    ' psddb8: BUG -',ch10,&
&    '  useylm must be 0 !'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

!  Read the data
   ekb0(:,:)=zero
   do itypat=1,ntypat
    read (nunit, '(13x,i4)' )ij
    do iproj=1,nproj
     read (nunit, '(6x,3d22.14)' )&
&     (ekb0(iproj+nproj*(ii-1),iproj+nproj*(ii-1)),ii=1,min(npsang,3))
     if(npsang>3)read (nunit, '(6x,3d22.14)' )&
&     (ekb0(iproj+nproj*(ii-1),iproj+nproj*(ii-1)),ii=4,npsang)
     do ii=1,npsang
      iekb=iproj+nproj*(ii-1)
      indlmn(1,iekb,itypat)=ii-1
      indlmn(2,iekb,itypat)=ii-1
      indlmn(3,iekb,itypat)=iproj
      indlmn(4,iekb,itypat)=ii**2-ii+1
      indlmn(5,iekb,itypat)=iekb
      indlmn(6,iekb,itypat)=1
!     For the time being, only diagonal ekb are treated in abinit v3
      ekb(iekb,itypat)=ekb0(iekb,iekb)
     end do
    end do
   end do

  else if(string==' No informat')then
   fullinit=0
  else
   write(message, '(a,a,a)' )&
&   ' psddb8 : BUG -',ch10,&
&   '  Error when reading the psp information'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Now, the number of blocks
  read(nunit,*)
  read(nunit,*)
  read(nunit, '(24x,i4)' )nblok

 else if(choice==2)then

! Write option

  write(nunit, '(a)' )' '

  if (fullinit/=0)then

   write(nunit, '(a)' )&
&   '  Description of the potentials (KB energies)'
   write(nunit, '(a,i6)' )&
&   '  vrsio8 (for pseudopotentials)=',vrsio8
   write(nunit, '(a,i3,a,i3,a,i3)' )&
&   '  dimekb =',dimekb,'       lmnmax=',lmnmax,&
&   '       usepaw=',usepaw
   do itypat=1,ntypat

!   Compute nekb
    nekb=0
    do jlmn=1,lmnmax
     jln=indlmn(5,jlmn,itypat)
     if(jln>nekb)then
      nekb=jln
     end if
    end do

    write(nunit, '(a,i4,a,i3,a,i4)' ) &
&    '  Atom type= ',itypat,'   pspso=',pspso(itypat),'   nekb=',nekb
    write(nunit, '(a)' ) '  iln lpsang iproj  ekb(:)'
    iln0=0
    ekb0(:,:)=zero
    do ilmn=1,lmnmax
     iln =indlmn(5,ilmn,itypat)
     if (iln>iln0) then
      iln0=iln
      lpsang=indlmn(1,ilmn,itypat)
      iproj=indlmn(3,ilmn,itypat)
!     For the time being, only diagonal ekb are treated in abinit v3
      ekb0(iln,iln)=ekb(iln,itypat)
!     For non-diagonal ekb, one could use:
!     do ii=iln to nekb
!     ekb0(ii,iln)=ekb(ii+iln*(iln-1)/2,itypat)
!     end do
      write(nunit, '(3i6,3x,4es15.7)' ) iln,lpsang,iproj,&
&      (ekb0(ii,iln),ii=1,min(nekb,4))
      if(nekb>4)then
       do iekb=5,nekb,4
        write(nunit, '(21x,4es15.7)' )&
&        (ekb0(ii,iekb),ii=iekb,min(nekb,iekb+3))
       end do
      end if
     end if

    end do
   end do

  else if(fullinit==0)then

!  This possibility is used when the DDB is initialized,
!  and the ekb s are not available from the GS input file...
   write(nunit, '(a)' )&
&   ' No information on the potentials yet '

  end if

! Now, write the number of blocks
  write(nunit, '(a)' )' '
  write(nunit, '(a)' )' **** Database of total energy derivatives ****'
  write(nunit, '(a,i4)' ) ' Number of data blocks= ',nblok

 end if

end subroutine psddb8
!!***
