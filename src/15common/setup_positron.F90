!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_positron
!! NAME
!! setup_positron
!!
!! FUNCTION
!! Do various setup calculations for the positron lifetime calculation
!!
!! NOTE
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!! vhae(nfft,2)
!! vhap(nfft,2)
!! rhorp(nfft,nspden)
!! rhore(nfft,nspden)
!! rhocore(nfft)
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      ioarr,mkcore,status
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setup_positron(dtset, fildensin,filstat,filvhain,mpi_enreg,n1,n2,n3,&
&                         hdr,iexit,level,nfft,n1xccc,ntypat,&
&                         xcccrc,xccc1d,etotal,rhocore,rhore,rhorp,rprimd,ucvol,&
&                         vhae,vhap,xred,ngfft)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_13xc
 use interfaces_14iowfdenpot
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iexit,level,n1,n1xccc,n2,n3,nfft,ntypat
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: etotal
 character(len=fnlen),intent(in) :: fildensin,filstat,filvhain
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(hdr_type),intent(inout) :: hdr
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rprimd(3,3),xccc1d(n1xccc,6,ntypat),xcccrc(ntypat)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(inout) :: rhocore(nfft),rhore(nfft,2),rhorp(nfft,2)
 real(dp),intent(inout) :: vhae(nfft,2),vhap(nfft,2)

!Local variables-------------------------------
!scalars
 integer :: fform,option,rdwr
!arrays
 real(dp) :: dummy2(6)
 real(dp),allocatable :: dummy3(:,:),dyfrx2(:,:,:)
 type(pawrhoij_type),allocatable :: rhoij_dum(:)

! *************************************************************************

 if (dtset%positron==1) then
! DEBUG
! print*,'Reading the electronic density and VHartree'
! ENDDEBUG

  call status(0,filstat,iexit,level,'call ioarr    ')
! Read rho(r) from a disk file
  rdwr=1
  fform=52
  call ioarr(0,rhore, dtset, etotal,fform,fildensin,hdr, mpi_enreg, &
&  nfft,rhoij_dum,rdwr,0,ngfft)

  call status(0,filstat,iexit,level,'call ioarr    ')
! Read vha(r) from a disk file
  rdwr=1
  fform=102
  call ioarr(0,vhae, dtset, etotal,fform,filvhain,hdr, mpi_enreg, &
&  nfft,rhoij_dum,rdwr,0,ngfft)

! Compute core electron density rhocore
  allocate(dummy3(3,dtset%natom))
  if (n1xccc/=0) then
   option=1
   allocate(dyfrx2(3,3,dtset%natom))
   call status(0,filstat,iexit,level,'call mkcore   ')
   call mkcore(dummy2,dyfrx2,dummy3,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&   n1,n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,dummy2,&
&   xcccrc,xccc1d,rhocore,xred)
   deallocate(dyfrx2)
  end if
  deallocate(dummy3)



 end if

 if (dtset%positron==2) then
! DEBUG
! print*,'Reading the electronic density, the positron density and the VHartree(e-p) potential'
! ENDDEBUG

  call status(0,filstat,iexit,level,'call ioarr    ')
! Read rho(r) from a disk file
  rdwr=1
  fform=52
  call ioarr(0,rhorp, dtset, etotal,fform,fildensin,hdr, mpi_enreg, &
&  nfft,rhoij_dum,rdwr,0,ngfft)

  call status(0,filstat,iexit,level,'call ioarr    ')
! Read vha(r) from a disk file
  rdwr=1
  fform=102
  call ioarr(0,vhap, dtset, etotal,fform,filvhain,hdr, mpi_enreg, &
&  nfft,rhoij_dum,rdwr,0,ngfft)

! Compute core electron density rhocore
  allocate(dummy3(3,dtset%natom))
  if (n1xccc/=0) then
   option=1
   allocate(dyfrx2(3,3,dtset%natom))
   call status(0,filstat,iexit,level,'call mkcore   ')
   call mkcore(dummy2,dyfrx2,dummy3,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&   n1,n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,dummy2,&
&   xcccrc,xccc1d,rhocore,xred)
   deallocate(dyfrx2)
  end if
  deallocate(dummy3)

 end if
end subroutine setup_positron
!!***
