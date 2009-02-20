!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfeig
!!
!! NAME
!! scfeig
!!
!! FUNCTION
!! Compute the largest eigenvalue and eigenvector of the SCF cycle.
!! A brute force algorithm is presently used.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG,GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  istep= number of the step in the SCF cycle
!!  i_vresid1=index of the residual potential in the array f_fftgr.
!!  i_vrespc1=index of the preconditioned residual pot. in the array f_fftgr.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nspden=number of spin-density components
!!  n_fftgr=third dimension of the array f_fftgr (should be equal to 5 here)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  vtrial(nfft,nspden)= at input, it is the trial potential that gave vresid .
!!       at output, it is an updated trial potential
!!  f_fftgr(nfft,nspden,n_fftgr)=different functions defined on the fft grid :
!!    the fixed point potential is kept in f_fftgr(:,:,1),
!!    the input residual potential (vresid) is in f_fftgr(:,:,i_vresid1),
!!    the input preconditioned residual potential (vrespc)
!!                                          is in f_fftgr(:,:,i_vrespc1),
!!    work space is provided in f_fftgr(:,:,4:5)
!!
!! PARENTS
!!      newrho,newvtr
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine scfeig(f_fftgr,istep,i_vresid1,i_vrespc1,nfft,nspden,n_fftgr,vtrial)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: i_vresid1,i_vrespc1,istep,n_fftgr,nfft,nspden
!arrays
 real(dp),intent(inout) :: f_fftgr(nfft,nspden,n_fftgr),vtrial(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ifft,isp
 real(dp) :: arg,eigen_scf,factor,fix_resid,norm,resid_new,resid_old
 character(len=500) :: message

! *************************************************************************

 if(nspden==4)then
  write(6,*)' scfeig : does not work yet for nspden=4'
  stop
 end if

!Set a fixed residual square for normalization of eigenvectors
 fix_resid=1.0d-4

!A few initialisations for the first istep
 if(istep==1)then

  write(message, '(a,es12.4,a,a,a,a,a,a,a)' )&
&  ' scfeig : fixed PC_residual square =',fix_resid,ch10,&
&  '    Note that fixed resid should always be much larger',ch10,&
&  '    than initial PC resid square, still sufficiently',ch10,&
&  '    small to reduce anharmonic effects ',ch10
  call wrtout(6,message,'COLL')

! Compute the preconditioned residual
! Also transfer vtrial in vtrial_old
  resid_old=0.0_dp
  do isp=1,nspden
   do ifft=1,nfft
    resid_old=resid_old+f_fftgr(ifft,isp,i_vrespc1)**2
   end do
  end do

  f_fftgr(:,:,1)=vtrial(:,:)
  write(message, '(a,es12.4)' )&
&  ' scfeig : initial PC_residual square =',resid_old
  call wrtout(6,message,'COLL')

  if(resid_old>1.0d-8)then
   write(message,'(a,a,a,a,a,a,a,a,a,a)') ch10,&
&   ' scfeig : ERROR -',ch10,&
&   '  This value is not good enough to allow',ch10,&
&   '  the computation of the eigenvectors of the SCF cycle.',ch10,&
&   '  It should be better than 1.0d-8 .',ch10,&
&   '  Action : improve the accuracy of your starting wavefunctions.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! In order to start the search for eigenvectors,
! use the tiny residual vector, renormalized
  factor=sqrt(fix_resid/resid_old)
  f_fftgr(:,:,4)=f_fftgr(:,:,i_vrespc1)*factor
  vtrial(:,:)=f_fftgr(:,:,1)+f_fftgr(:,:,4)

! If istep is not equal to 1
 else if(istep>=2)then
! 
! Compute the corresponding operator expectation value
! And put the residual vector minus the difference
! between vtrial and vtrial_old
! (this is actually the action of the operator !) in vect(*,2)
  eigen_scf=0.0_dp
  do isp=1,nspden
   do ifft=1,nfft
    eigen_scf=eigen_scf+&
&    f_fftgr(ifft,isp,4) * f_fftgr(ifft,isp,i_vrespc1)
   end do
  end do

  do isp=1,nspden
   do ifft=1,nfft
    f_fftgr(ifft,isp,i_vrespc1)=f_fftgr(ifft,isp,i_vrespc1)&
&    +vtrial(ifft,isp)-f_fftgr(ifft,isp,1)
    f_fftgr(ifft,isp,5)=f_fftgr(ifft,isp,i_vrespc1)
   end do
  end do
  eigen_scf=eigen_scf/fix_resid
  write(message, '(a,es12.4,a)' ) &
&  ' scfeig : Operator expectation value ',eigen_scf,' (extremal eigenvalue * diemix)'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')
! 
! Compute residual of vect(*,2)
  resid_new=zero
  do isp=1,min(nspden,2)
   do ifft=1,nfft
    resid_new=resid_new+ f_fftgr(ifft,isp,5) ** 2
   end do
  end do
  if (nspden==4) then
   do ifft=1,nfft
    resid_new=resid_new+two*(f_fftgr(ifft,3,5)**2+f_fftgr(ifft,4,5)**2)
   end do
  end if
  factor=sqrt(fix_resid/resid_new)
  if(eigen_scf<zero) then
   factor=-factor ! the new vector MAY be oposite to the old one
!  if(factor<-one) factor=-factor ! the new vector is not opposed to the old one
  end if
  write(message, '(a,es12.4)' ) &
&  ' scfeig : Inverse of renormalization factor ',one/factor
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')
  write(message, '(a,es12.4)' ) &
&  ' scfeig : Convergence criterion value (->0 at convergency) ',one/factor-eigen_scf-one
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

  f_fftgr(:,:,4)=f_fftgr(:,:,5)*factor
  vtrial(:,:)=f_fftgr(:,:,1)+f_fftgr(:,:,4)
! 
! End the different istep cases
 end if

end subroutine scfeig
!!***
