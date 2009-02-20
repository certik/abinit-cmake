!{\src2tex{textfont=tt}}
!!****f* ABINIT/ylmcd
!! NAME
!! ylmcd
!!
!! FUNCTION
!!
!!  This subroutine computes dth and dphi, the first derivatives of
!!  (complex) Ylm as a function of th and phi (the angles of the spherical coordinates)
!!  It works for all spherical harmonics with il <=2
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (FB, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  il=angular quantum number
!!  im=magnetic quantum number
!!  kcart=cartesian coordinates of the vector where the first derivatives of Ylm are evaluated
!!
!! OUTPUT
!!  dth=derivative of Y_lm with respect to \theta
!!  dphi=derivative of Y_lm with respect to \phi
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Using double precision complex even thou it does not conform the ABINIT
!!  convention
!!  case l=3 must be implemented
!!
!! PARENTS
!!      ccgradvnl_ylm
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ylmcd(il,im,kcart,dth,dphi)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,im
 complex(dpc),intent(out) :: dphi,dth
!arrays
 real(dp),intent(in) :: kcart(3)

!Local variables-------------------------------
!scalars
 integer,parameter :: lx=3
 real(dp),parameter :: ppad=tol8
 real(dp) :: cosphi,costh,costhreephi,costwophi,r,rxy,sinphi,sinth,sinthreephi
 real(dp) :: sintwophi
 character(len=500) :: msg

! *************************************************************************

 if (ABS(im)>ABS(il))then
  write(msg,'(4a,i6,2a,i6,a,i6)') ch10,&
&  ' ylmcd: ERROR -',ch10,&
&  '  m is,',im,ch10,&
&  '  however it should be between ',-il,' and ',il
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 r=SQRT(kcart(1)**2+kcart(2)**2+kcart(3)**2)
 if (r<ppad) r=r+ppad
 rxy=SQRT(kcart(1)**2+kcart(2)**2)
 if (rxy<ppad) rxy=r+ppad

!(th,ph) spherical representation.
 costh= kcart(3)/r
 sinth= rxy/r
 cosphi= kcart(1)/rxy
 sinphi= kcart(2)/rxy
 costwophi= two*cosphi**2 - one
 sintwophi= two*sinphi*cosphi
 costhreephi=cosphi*costwophi-sinphi*sintwophi
 sinthreephi=cosphi*sintwophi+sinphi*costwophi

 select case (il)
  case (0)
   dth  = (0.d0,0.d0)
   dphi = (0.d0,0.d0)
  case (1)
   if (ABS(im)==0) then
    dth= -SQRT(three/(four_pi))*sinth
    dphi= (0.d0,0.d0)
   else if (abs(im)==1) then
    dth= -SQRT(3.d0/(8.d0*pi))*costh*CMPLX(cosphi,sinphi)
    dphi=-SQRT(3.d0/(8.d0*pi))*sinth*CMPLX(-sinphi,cosphi)
   end if
  case (2)
   if (ABS(im)==0) then
    dth= -SQRT(5.d0/(16.d0*pi))*6.d0*costh*sinth
    dphi= (0.d0,0.d0)
   else if (ABS(im)==1) then
    dth=  -SQRT(15.d0/(8.d0*pi))*(costh**2-sinth**2)*CMPLX(cosphi,sinphi)
    dphi= -SQRT(15.d0/(8.d0*pi))*costh*sinth*(0.d0,1.d0)*CMPLX(cosphi,sinphi)
   elseif (abs(im)==2) then
    dth  = SQRT(15.d0/(32.d0*pi))*2.d0*costh*sinth*CMPLX(costwophi,sintwophi)
    dphi = SQRT(15.d0/(32.d0*pi))*sinth**2*(0.d0,2.d0)*CMPLX(costwophi,sintwophi)
   end if
  case (3)
   if (ABS(im)==0) then
    dth = SQRT(7.d0/(16*pi))*(-15.d0*costh**2*sinth + 3.d0**sinth)
    dphi= (0.d0,0.d0)
   else if (ABS(im)==1) then
    dth= -SQRT(21.d0/(64.d0*pi))*CMPLX(cosphi,sinphi)*(5.d0*costh**3-costh-10.d0*sinth**2*costh)
    dphi=-SQRT(21.d0/(64.d0*pi))*sinth*(5.d0*costh**2-1)*(0.d0,1.d0)*CMPLX(cosphi,sinphi)
   else if (ABS(im)==2) then
    dth =SQRT(105.d0/(32.d0*pi))*(2.d0*sinth*costh**2-sinth**3)*CMPLX(costwophi,sintwophi)
    dphi=SQRT(105.d0/(32*pi))*sinth**2*costh*(0.d0,2.d0)*CMPLX(costwophi,sintwophi)
   else if (abs(im)==3) then
    dth =-SQRT(35.d0/(64.d0*pi))*3.d0*sinth**2*costh*CMPLX(costhreephi,sinthreephi)
    dphi= SQRT(35.d0/(64.d0*pi))*sinth**3*(0.d0,3.d0)*CMPLX(costhreephi,sinthreephi)
   end if
   case default
   write(msg,'(4a,i6,2a,i6)')ch10,&
&   ' ylmcd: ERROR -',ch10,&
&   '  The maximum allowed value for l is,',lx,ch10,&
&   '  however, l=',il
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
 end select
!
!Treat the case im < 0
 if (im<0) then
  dth = (-one)**(im)*CONJG(dth)
  dphi= (-one)**(im)*CONJG(dphi)
 end if

end subroutine ylmcd
!!***
