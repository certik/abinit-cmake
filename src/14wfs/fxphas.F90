!{\src2tex{textfont=tt}}
!!****f* ABINIT/fxphas
!!
!! NAME
!! fxphas
!!
!! FUNCTION
!! Fix phase of all bands. Keep normalization but maximize real part
!! (minimize imag part). Also fix the sign of real part
!! by setting the first non-zero element to be positive.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mcg)= contains the wavefunction |c> coefficients.
!!  gsc(2,mgsc)= if useoverlap==1, contains the S|c> coefficients,
!!               where S is an overlap matrix.
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  istwfk=input option parameter that describes the storage of wfs
!!    (set to 1 if usual complex vectors)
!!  mcg=size of second dimension of cg
!!  mgsc=size of second dimension of gsc
!!  mpi_enreg=informations about MPI parallelization
!!  nband_k=number of bands
!!  npw_k=number of planewaves
!!  useoverlap=describe the overlap of wavefunctions:
!!               0: no overlap (S=Identi0,ty_matrix)
!!               1: wavefunctions are overlapping
!!
!! OUTPUT
!!  cg(2,mcg)=same array with altered phase.
!!  gsc(2,mgsc)= same array with altered phase.
!!
!! NOTES
!! When the sign of the real part was fixed (modif v3.1.3g.6), the
!! test Tv3#5 , dataset 5, behaved differently than previously.
!! This should be cleared up.
!!
!! PARENTS
!!      phfrq3,vtowfk
!!
!! CHILDREN
!!      leave_new,timab,wrtout,xcomm_init,xmaster_init_fft,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fxphas(cg,gsc,icg,igsc,istwfk,mcg,mgsc,mpi_enreg,nband_k,npw_k,useoverlap)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,igsc,istwfk,mcg,mgsc,nband_k,npw_k,useoverlap
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc*useoverlap)

!Local variables-------------------------------
!scalars
 integer :: iband,ierr,ii,index,master,old_paral_level,spacecomm
 real(dp) :: cim,cre,gscim,gscre,norm,quotient,root1,root2,saa,sab,sbb,theta
 real(dp) :: thppi,xx,yy
 character(len=500) :: message
!arrays
 real(dp) :: buffer2(nband_k,2),buffer3(nband_k,3),tsec(2)
 real(dp),allocatable :: cimb(:),creb(:),saab(:),sabb(:),sbbb(:)

! *************************************************************************

!DEBUG
!write(6,*)' fxphas : enter '
!ENDDEBUG

!The general case, where a complex phase indeterminacy is present
 if(istwfk==1)then

  allocate(cimb(nband_k),creb(nband_k),saab(nband_k),sabb(nband_k),sbbb(nband_k))
  cimb(:)=zero ; creb(:)=zero
! Loop over bands
  do iband=1,nband_k

   index=icg+(iband-1)*npw_k

!  Compute several sums over Re, Im parts of c
   saa=0.0_dp ; sbb=0.0_dp ; sab=0.0_dp
!  $OMP PARALLEL DO PRIVATE(ii) REDUCTION(+:saa,sbb,sab) &
!  $OMP&SHARED(cg,index,npw_k)
   do ii=1+index,npw_k+index
    saa=saa+cg(1,ii)*cg(1,ii)
    sbb=sbb+cg(2,ii)*cg(2,ii)
    sab=sab+cg(1,ii)*cg(2,ii)
   end do
!  $OMP END PARALLEL DO
   saab(iband)=saa
   sbbb(iband)=sbb
   sabb(iband)=sab

  end do ! iband

! XG030513 : MPIWF : should transmit saab,sbbb,sabb from the procs
! of the WF group to the master processor of the WF group
  if (mpi_enreg%paral_fft == 1) then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   buffer3(:,1)=saab(:)
   buffer3(:,2)=sbbb(:)
   buffer3(:,3)=sabb(:)
   call timab(48,1,tsec)
   call xsum_mpi(buffer3,spaceComm,ierr)
!  call xsum_mpi(saab,spaceComm,ierr)
!  call xsum_mpi(sbbb,spaceComm,ierr)
!  call xsum_mpi(sabb,spaceComm,ierr)
   call timab(48,2,tsec)
   saab(:)=buffer3(:,1)
   sbbb(:)=buffer3(:,2)
   sabb(:)=buffer3(:,3)
   mpi_enreg%paral_level=old_paral_level
  end if

! XG030513 : MPIWF this loop should only be executed by the
! master of the WF group
  call xmaster_init_fft(mpi_enreg,master)
  if (mpi_enreg%me==master) then

   do iband=1,nband_k

    index=icg+(iband-1)*npw_k

    saa=saab(iband)
    sbb=sbbb(iband)
    sab=sabb(iband)

!   Get phase angle theta
    if (sbb+saa>tol8)then
     if (abs(sbb-saa)>tol8*abs(sab)) then
      quotient=sab/(sbb-saa)
      theta=0.5_dp*atan(2.0_dp*quotient)
     else
!     Taylor expansion of the atan in terms of inverse of its argument. Correct up to 1/x2, included.
      theta=0.25_dp*(pi-(sbb-saa)/sab)
     end if
    else
     write(message,'(a,a,a,i4,a)' )&
&     ' fxphas : BUG -',ch10,&
&     '  The eigenvector number',iband,' has zero norm.'
     call wrtout(6,message,'PERS')
     call leave_new('PERS')
    end if

!   Check roots to get theta for max Re part
    root1=cos(theta)**2*saa+sin(theta)**2*sbb-&
&    2.0_dp*cos(theta)*sin(theta)*sab
    thppi=theta+0.5_dp*pi
    root2=cos(thppi)**2*saa+sin(thppi)**2*sbb-&
&    2.0_dp*cos(thppi)*sin(thppi)*sab
    if (root2>root1) theta=thppi

    xx=cos(theta)
    yy=sin(theta)

!   Here, set the first non-zero element to be positive
!   Comment the next nine lines to recover the behaviour of pre v3.1.3g
    do ii=1+index,npw_k+index
     cre=cg(1,ii)
     cim=cg(2,ii)
     cre=xx*cre-yy*cim
     if(abs(cre)>tol8)exit
    end do
    if(cre<zero)then
     xx=-xx ; yy=-yy
    end if

    creb(iband)=xx
    cimb(iband)=yy

   end do
  end if
! XG030513 : MPIWF : should transmit creb(:),cimb(:) of the master
! processor of the WF group to the others procs of the WF group
  if (mpi_enreg%paral_fft == 1) then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   call timab(48,1,tsec)
   buffer2(:,1)=creb(:)
   buffer2(:,2)=cimb(:)
   call xsum_mpi(buffer2,spaceComm,ierr)
!  call xsum_mpi(creb,spaceComm,ierr)
!  call xsum_mpi(cimb,spaceComm,ierr)
   call timab(48,2,tsec)
   creb(:)=buffer2(:,1)
   cimb(:)=buffer2(:,2)
   mpi_enreg%paral_level=old_paral_level
  end if

  do iband=1,nband_k

   index=icg+(iband-1)*npw_k

   xx=creb(iband)
   yy=cimb(iband)
!  Alter phase of array |cg>
!  $OMP PARALLEL DO PRIVATE(cim,cre,ii) &
!  $OMP&SHARED(cg,index,npw_k,xx,yy)
   do ii=1+index,npw_k+index
    cre=cg(1,ii)
    cim=cg(2,ii)
    cg(1,ii)=xx*cre-yy*cim
    cg(2,ii)=xx*cim+yy*cre
   end do
!  $OMP END PARALLEL DO

!  Alter phase of array S|cg>
   if (useoverlap==1) then
    index=igsc+(iband-1)*npw_k
!   $OMP PARALLEL DO PRIVATE(gscim,gscre,ii) &
!   $OMP&SHARED(gsc,index,npw_k,xx,yy)
    do ii=1+index,npw_k+index
     gscre=gsc(1,ii)
     gscim=gsc(2,ii)
     gsc(1,ii)=xx*gscre-yy*gscim
     gsc(2,ii)=xx*gscim+yy*gscre
    end do
!   $OMP END PARALLEL DO
   end if

  end do ! iband

  deallocate(cimb,creb,saab,sabb,sbbb)

! ====================================================================

! Storages that take into account the time-reversal symmetry :
! the freedom is only a sign freedom
 else  ! if istwfk/=1

  allocate(creb(nband_k))
  creb(:)=zero
! XG030513 : MPIWF : this loop should be done only by the master
! processor of the WF group

  call xmaster_init_fft(mpi_enreg,master)
  if (mpi_enreg%me==master) then

!  Loop over bands
   do iband=1,nband_k

    index=icg+(iband-1)*npw_k

!   Here, set the first non-zero real element to be positive
    do ii=1+index,npw_k+index
     cre=cg(1,ii)
     if(abs(cre)>tol8)exit
    end do
    creb(iband)=cre

   end do ! iband

  end if
! XG030513 : MPIWF : should transmit cre(:) of the master
! processor of the WF group to the others
  if (mpi_enreg%paral_fft == 1) then
   old_paral_level=mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   call timab(48,1,tsec)
   call xsum_mpi(creb,spaceComm,ierr)
   call timab(48,2,tsec)
   mpi_enreg%paral_level=old_paral_level
  end if
  do iband=1,nband_k


   cre=creb(iband)
   if(cre<zero)then
    index=icg+(iband-1)*npw_k
    do ii=1+index,npw_k+index
     cg(1,ii)=-cg(1,ii)
     cg(2,ii)=-cg(2,ii)
    end do
    if(useoverlap==1)then
     index=igsc+(iband-1)*npw_k
     do ii=1+index,npw_k+index
      gsc(1,ii)=-gsc(1,ii)
      gsc(2,ii)=-gsc(2,ii)
     end do
    end if
   end if

  end do ! iband

  deallocate(creb)

 end if ! istwfk

!DEBUG
!write(6,*)' fxphas : exit'
!do iband=1,nband_k
!index=icg+(iband-1)*npw_k
!do ii=1+index,npw_k+index
!write(6,*)ii,cg(1,ii),cg(2,ii)
!end do
!end do
!ENDDEBUG

end subroutine fxphas
!!***
