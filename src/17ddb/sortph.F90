!{\src2tex{textfont=tt}}
!!****f* ABINIT/sortph
!! NAME
!! sortph
!!
!! FUNCTION
!! Sort the energies in order to have fine phonon
!! dispersion curves
!! It is best not to include the gamma point in the list
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (FDortu,MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  displ(2,3*natom,3*natom)= contain
!!   the displacements of atoms in cartesian coordinates.
!!   The first index means either the real or the imaginary part,
!!   The second index runs on the direction and the atoms displaced
!!   The third index runs on the modes.
!!  filnam=name of output files
!!   hacmm1,hartev,harthz,xkb= different conversion factors
!!  iout= unit for long print (if negative, the routine only print
!!        on unit 6, and in Hartree only).
!!  natom= number of atom
!!  phfrq(3*natom)= phonon frequencies in Hartree
!!  qphnrm=phonon wavevector normalisation factor
!!  qphon(3)=phonon wavevector
!!  udispl=unit number for output of phonon eigendisplacements
!!  ufreq=unit number for output of phonon frequencies
!!
!! OUTPUT
!!  (only writing ?)
!!
!! NOTES
!! Called by one processor only
!!
!! PARENTS
!!      anaddb,thm9
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine sortph(displ,filnam,&
     &  iout,natom,phfrq,qphnrm,qphon,udispl,ufreq)

 use defs_basis
   use defs_datatypes

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: iout,natom,udispl,ufreq
 real(dp),intent(in) :: qphnrm
 character(len=fnlen),intent(in) :: filnam
!arrays
 real(dp),intent(in) :: displ(2,3*natom,3*natom),phfrq(3*natom),qphon(3)

!Local variables-------------------------------
  !displNew : hold the sorted dislpacements
  !phfrqNew : hold the sorted frequencies
!scalars
 integer,save :: nappels=0
 integer :: bandmax,bandmin,i,iatom,idir,ii,iimode,imode,modeNum,whichqm1
 real(dp) :: deltaE,deltaqm1,deltaqm1prev,displmodq,displmodqm1,displmodqm1prev
 character(len=500) :: message
 character(len=fnlen) :: file_displ,file_freq
!arrays
 integer,allocatable,save :: displLocqm1(:)
 integer :: actualMinPos(3*natom,natom),displLoc(3*natom)
 real(dp),allocatable,save :: displqm1(:,:,:)
 real(dp) :: actualMin(3*natom,natom),displNew(2,3*natom,3*natom),mean(3*natom)
 real(dp) :: phfrqNew(3*natom),variance(3*natom)

  ! *********************************************************************

 write(6, '(a)' )' sortph : enter '

!which band must be considered (in Hartree) :
!band are swaped only if they are not to far from each other
!deltaE=0.0003        !This should be set by the user (fine for ABO3)
 deltaE=0.00003       !fine for BaO
 actualMin(:,:)=1000   !Something arbitrary large enough
 nappels=nappels+1


!*********************************************************************
!We calculate the difference between the displacement at current q and
!the displacement at previous q

 if(allocated(displqm1)) then

  do imode=1,3*natom
   do iatom=1,natom


!   Only compare the eigenmodes whose the energy is in
!   the range 'deltaE' higher or lower
!   But almost one band lower and higher are considered

    bandmin=-1
    bandmax=1

    if(imode+bandmin>=1) then
     do while(abs(phfrq(imode)-phfrq(imode+bandmin))<deltaE)
      bandmin=bandmin-1
      if(imode+bandmin<1) exit
     end do
    end if


    if(imode+bandmax<=3*natom) then
     do while(abs(phfrq(imode)-phfrq(imode+bandmax))<deltaE)
      bandmax=bandmax+1
      if(imode+bandmax>3*natom) exit
     end do
    end if



!   Calculate the (squared) modulus of the displacement of the current atome (iatom)
!   corresponding to imode at current q
!   'displmodq' stands for 'displacement modulus for the current mode (and atom) at q'
    displmodq=displ(1,3*(iatom-1)+1,imode)**2 &
&    +displ(2,3*(iatom-1)+1,imode)**2 &
&    +displ(1,3*(iatom-1)+2,imode)**2 &
&    +displ(2,3*(iatom-1)+2,imode)**2 &
&    +displ(1,3*(iatom-1)+3,imode)**2 &
&    +displ(2,3*(iatom-1)+3,imode)**2



!   Calculate the (squared) modulus of the displacement of the current atome (iatom)
!   corresponding to imode+bandmin at q-1
!   bandmin and bandmax are relative index to the current band
!   ex: bandmin=-1 and bandmax=1
!   'displmodqm1prev' stands for
!   'displacement modulus for the previous mode at q-1'

    if(imode+bandmin >= 1) then
     displmodqm1prev=displqm1(1,3*(iatom-1)+1,imode+bandmin)**2 &
&     +displqm1(2,3*(iatom-1)+1,imode+bandmin)**2 &
&     +displqm1(1,3*(iatom-1)+2,imode+bandmin)**2 &
&     +displqm1(2,3*(iatom-1)+2,imode+bandmin)**2 &
&     +displqm1(1,3*(iatom-1)+3,imode+bandmin)**2 &
&     +displqm1(2,3*(iatom-1)+3,imode+bandmin)**2
    else
     displmodqm1prev=1000   !something arbitrary large enough
    end if

!   Find which of the band at q-1 is closer than the one at q

    do whichqm1=bandmin+1,bandmax
     if((imode+whichqm1 <= 3*natom).and.(imode+whichqm1 >= 1)) then
      displmodqm1=displqm1(1,3*(iatom-1)+1,imode+whichqm1)**2 &
&      +displqm1(2,3*(iatom-1)+1,imode+whichqm1)**2 &
&      +displqm1(1,3*(iatom-1)+2,imode+whichqm1)**2 &
&      +displqm1(2,3*(iatom-1)+2,imode+whichqm1)**2 &
&      +displqm1(1,3*(iatom-1)+3,imode+whichqm1)**2 &
&      +displqm1(2,3*(iatom-1)+3,imode+whichqm1)**2
     else
      displmodqm1=1000     !something arbitrary large enough
     end if



     deltaqm1=abs(displmodq-displmodqm1)
     deltaqm1prev=abs(displmodq-displmodqm1prev)


     if(deltaqm1<deltaqm1prev) then

      if(deltaqm1<actualMin(imode,iatom)) then
       actualMinPos(imode,iatom)=whichqm1
       actualMin(imode,iatom)=deltaqm1

      end if
     else
      if(deltaqm1prev<actualMin(imode,iatom)) then

       actualMinPos(imode,iatom)=whichqm1-1
       actualMin(imode,iatom)=deltaqm1prev

      end if
     end if




     displmodqm1prev=displmodqm1
    end do        !loop whichqm1


   end do       !big loop iatom

  end do       !big loop imode

! At the end of the loop, the only thing that is of interest is
! actualMinPos(imode,iatom) that gives for the current q point
! which band at previous q is of the same 'family'
! ex : actualMinPos(5,3) can give -1, 0 or +1 (or -2,-1,0,1,2 if more neboors are
! explored


! Determine what must swap according 'actualMinPos'


  mean(:)=0
  variance(:)=0
  do imode=1,3*natom

!  We calculate the mean and the variance
   do iatom=1,natom
    mean(imode)=mean(imode)+actualMinPos(imode,iatom)
   end do

   mean(imode)=mean(imode)/natom

   do iatom=1,natom
    variance(imode)=variance(imode)+(actualMinPos(imode,iatom)-mean(imode))**2
   end do
   variance(imode)=variance(imode)/natom

   mean(imode)=nint(mean(imode))

!  The 1 below is somewhat arbitrary
   if(variance(imode)>1) then
    mean(imode)=0
   end if


!  if the difference between two band is higher than deltaE, ignore the swap
!  This should have been set earlier in the code but it was problematic...
   if((imode-1) > 0) then
    if(abs(phfrq(imode)-phfrq(imode-1)) > deltaE) then
     mean(imode)=0
    end if
   end if
   if((imode+1)<3*natom) then
    if(abs(phfrq(imode)-phfrq(imode+1)) > deltaE) then
     mean(imode)=0
    end if
   end if



  end do


! Determine if coherent
! ex : if a band want to go up, the band just above must want to go down

  do imode=1,3*natom
   if(-int(mean(imode+int(mean(imode))))/=int(mean(imode))) then
    mean(imode)=0
    write(6, '(a)' )' sortph :a ambiguous swap -> do nothing  '
    write(6, '(a,3f10.4)' )' sortph :a     qpt= ', qphon(:)
    write(6, '(a,i5)' )' sortph :a     imode= ',imode
   end if

   if(int(mean(imode))/=0) then
    write(6, '(a)' )' sortph :s  swap '
    write(6, '(a,3f10.4)' )' sortph :s     qpt ', qphon(:)
    write(6, '(a,i4)' )' sortph :s     mode ', imode
    write(6, '(a,i4)' )' sortph :s     swap with ', int(mean(imode))
   end if
  end do



! Swap the eigenvalues

  do imode=1,3*natom
   modeNum=imode+mean(imode)
   do iimode=1,3*natom
!   Find the position(iimode) of imode
    if(allocated(displLocqm1)) then
     if(displLocqm1(iimode)==modeNum) then
      exit
     end if
    end if
   end do

   displNew(:,:,iimode)=displ(:,:,imode)
   phfrqNew(iimode)=phfrq(imode)
   displLoc(iimode)=imode

  end do


! Write frequencies in a file
! CE FORMAT DOIT ETRE REECRIT CAR CAS PARTICULIER 3*natom=15
  write(ufreq,'(15(5d14.6))')  (phfrqNew(ii),ii=1,3*natom)


! write displacements in a file
  do imode=1,3*natom

   do iatom=1,natom

    write(udispl,'(d14.6)') sqrt(displNew(1,3*(iatom-1)+1,imode)**2 &
&    +displNew(2,3*(iatom-1)+1,imode)**2 &
&    +displNew(1,3*(iatom-1)+2,imode)**2 &
&    +displNew(2,3*(iatom-1)+2,imode)**2 &
&    +displNew(1,3*(iatom-1)+3,imode)**2 &
&    +displNew(2,3*(iatom-1)+3,imode)**2)

   end do
  end do

  displLocqm1=displLoc


 else
! This part of the code is only executed once, at the first call

  file_freq = trim(filnam)//".freq"
  file_displ = trim(filnam)//".displ"

  write(6, '(a)' )' sortph : allocating displqm1 '
  allocate(displqm1(2,3*natom,3*natom))
  write(6, '(a)' )' sortph : displqm1 allocated '

  write(6, '(a,a)' )' sortph : opening file ',trim(file_freq)
  open(ufreq,FILE=trim(file_freq),&
&  STATUS='replace',ACCESS='sequential',ACTION='write')
  write(6, '(a,a,a)' )' sortph : file ',trim(file_freq),' opened '

  write(6, '(a,a)' )' sortph : opening file ',trim(file_displ)
  open(udispl,FILE=trim(file_displ),&
&  STATUS='replace',ACCESS='sequential',ACTION='write')
  write(6, '(a,a,a)' )' sortph : file ',trim(file_displ),' opened '

  write(6, '(a)' )' sortph : allocating displLocqm1 '
  allocate(displLocqm1(3*natom))
  write(6, '(a)' )' sortph : displLocqm1 allocated '
  displNew(:,:,:)=displ(:,:,:)

  do imode=1,3*natom
   displLocqm1(imode)=imode

  end do


! Write frequencies and diplacements unchanged
  write(ufreq,'(15(5d14.6))') (phfrq(ii),ii=1,3*natom)
  do imode=1,3*natom

   do iatom=1,natom

    write(udispl,'(d14.6)') sqrt(displ(1,3*(iatom-1)+1,imode)**2 &
&    +displ(2,3*(iatom-1)+1,imode)**2 &
&    +displ(1,3*(iatom-1)+2,imode)**2 &
&    +displ(2,3*(iatom-1)+2,imode)**2 &
&    +displ(1,3*(iatom-1)+3,imode)**2 &
&    +displ(2,3*(iatom-1)+3,imode)**2)

   end do
  end do

  actualMinPos(:,:)=0    !for the first q, the modes are in the right order(by definition)
 end if


 displqm1=displ

!**********************************************************************************



 write(6, '(a)' )' sortph : exit '

end subroutine sortph
!!***
