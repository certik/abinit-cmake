!{\src2tex{textfont=tt}}
!!****f* ABINIT/cppm1par
!! NAME
!! cppm1par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters big-omega-twiddle-squared and omega-twiddle from
!! epsilon-twiddle^-1 calculated for nomega (usually 2) frequencies omega=0 and omega=iE0.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  epsm1(npwc,npwc,nomega,nqibz)=dielectric matrix at nomega frequencies, and nqibz wavevectors
!!  npwc=number of plane waves
!!  nomega=number of frequencies (usually 2)
!!  nqibz=number of irreducible q-points
!!  omega(nomega)=frequencies
!!  omegaplasma=input variable or Drude plasma frequency
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc,nqibz)=parameter of the plasmon-pole model (see gwa.pdf file)
!!  omegatw(npwc,npwc,nqibz)=parameter of the plasmon-pole model (see gwa.pdf file)
!!
!! NOTES 
!!  Note the use of intent(inout) since elements of data types are supposed 
!!  to be passed this routine
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cppm1par(npwc,nqibz,nomega,epsm1,omega,bigomegatwsq,omegatw,omegaplasma)

 use defs_basis
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,nqibz
 real(dp),intent(in) :: omegaplasma
!arrays
 complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega,nqibz),omega(nomega)
 complex(gwpc),intent(inout) :: bigomegatwsq(npwc,npwc,nqibz)
 complex(gwpc),intent(inout) :: omegatw(npwc,npwc,nqibz)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,io,io0,ioe0,iq
 real(dp) :: e0,minomega
 logical :: ltest
 character(len=500) :: msg
!arrays
 real(dp) :: q(3,nqibz),qplusg(npwc)
 complex(gwpc),allocatable :: AA(:,:),omegatwsq(:,:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(2a)')ch10,' cppm1par : enter '
 call wrtout(std_out,msg,'COLL')
#endif

 allocate(AA(npwc,npwc),omegatwsq(npwc,npwc))
 !
 ! === Find omega=0 and omega=imag (closest to omegaplasma) where to fit ppm parameters ===
 minomega=1.0d-3 ; io0=0
 do io=1,nomega
  if (ABS(omega(io))<minomega) then
   io0=io ; minomega=ABS(omega(io))
  end if
 end do
 ltest=(io0/=0) 
 call assert(ltest,'omega=0 not found',__FILE__,__LINE__)

 minomega=1.0d-3 ; e0=200.0 ; ioe0=0
 do io=1,nomega
  if (REAL(omega(io))<minomega.and.AIMAG(omega(io))>minomega) then
   if (ABS(AIMAG(omega(io))-omegaplasma)<ABS(e0-omegaplasma)) then
    ioe0=io ; e0=AIMAG(omega(io))
   end if
  end if
 end do
 ltest=(ioe0/=0)
 call assert(ltest,'imaginary omega not found',__FILE__,__LINE__)

 do iq=1,nqibz
  ! === Calculate plasmon-pole A parameter A=epsilon^-1(0)-delta ===
  AA(:,:)=epsm1(:,:,io0,iq)
  do ig=1,npwc
   AA(ig,ig)=AA(ig,ig)-one
  end do
  ! === Calculate plasmon-pole omega-twiddle-square parameter ===
  omegatwsq(:,:)=(AA(:,:)/(epsm1(:,:,io0,iq)-epsm1(:,:,ioe0,iq))-one)*e0**2
  !
  ! If omega-twiddle-squared is negative,set omega-twiddle-squared to 1.0 (a reasonable way of treating
  ! such terms, in which epsilon**-1 was originally increasing along this part of the imaginary axis)
  ! (note: originally these terms were ignored in Sigma; this was changed on 6 March 1990.)
  WHERE (REAL(omegatwsq)<=zero) 
   omegatwsq=one
  END WHERE
  !do igp=1,npwc
  ! do ig=1,npwc
  !  if (REAL(omegatwsq(ig,igp))<=zero) omegatwsq(ig,igp)=one
  ! end do
  !end do
  !
  ! === Get omega-twiddle ===
  ! * Neglect the imag part (if one) in omega-twiddle-squared
  omegatw(:,:,iq)=SQRT(REAL(omegatwsq(:,:)))
  !
  ! === Get big-omega-twiddle-squared=-omega-twiddle-squared AA ===
  bigomegatwsq(:,:,iq)=-AA(:,:)*omegatw(:,:,iq)**2
 end do
 deallocate(AA,omegatwsq)

 write(msg,'(2a,f15.12)')ch10,&
& ' cppm1par : omega twiddle minval [eV] = ',MINVAL(ABS(omegatw(:,:,:)))*Ha_eV
 call wrtout(std_out,msg,'COLL')
 write(msg,'(2a,3i5,a)')ch10,&
& '           omega twiddle min location',MINLOC(ABS(omegatw(:,:,:))),ch10
 call wrtout(std_out,msg,'COLL')

#if defined DEBUG_MODE
 write(msg,'(a)')' cppm1par : exit '
 call wrtout(std_out,msg,'COLL')
#endif

end subroutine cppm1par
!!***
