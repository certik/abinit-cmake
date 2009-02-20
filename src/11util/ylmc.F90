!{\src2tex{textfont=tt}}
!!****f* ABINIT/ylmc
!! NAME
!! ylmc
!!
!! FUNCTION
!!  Calculate all (complex) spherical harmonics for il<=3
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
!!  kcart=vector in cartesian coordinates defining the value of \theta and \psi
!!   where calculate the spherical harmonic
!!
!! OUTPUT
!!  ylm= spherical harmonic
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Using double precision complex even thou it does not conform the ABINIT
!!  convention
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 function ylmc(il,im,kcart)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: il,im
 complex(dpc) :: ylmc
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

 if (ABS(im)>ABS(il)) then
  write(msg,'(4a,i6,2a,i6,a,i6)')ch10,&
&  ' ylmc: ERROR -',ch10,&
&  '  m is,',im,ch10,&
&  '  however it should be between ',-il,' and ',il
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 r=SQRT(kcart(1)**2+kcart(2)**2+kcart(3)**2)
 if (r<ppad) r=r+ppad
 rxy=SQRT(kcart(1)**2+kcart(2)**2)
 if (rxy<ppad)rxy=r+ppad
!
!(th,phi) spherical coordinates
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
   ylmc= one/SQRT(four_pi)
  case (1)
   if (ABS(im)==0) then
    ylmc = SQRT(three/(four_pi))*costh
   else if (ABS(im)==1) then
    ylmc = -SQRT(three/(eight*pi))*sinth*CMPLX(cosphi,sinphi)
   end if
  case (2)
   if (ABS(im)==0) then
    ylmc = SQRT(5.d0/(16.d0*pi))*(three*costh**2-one)
   else if (ABS(im)==1) then
    ylmc = -SQRT(15.d0/(8.d0*pi))*sinth*costh*cmplx(cosphi,sinphi)
   else if (ABS(im)==2) then
    ylmc = SQRT(15.d0/(32.d0*pi))*(sinth)**2*CMPLX(costwophi,sintwophi)
   end if
  case (3)
   if (ABS(im)==0) then
    ylmc= SQRT(7.d0/(16.d0*pi))*(5.d0*costh**3 -3.d0*costh)
   else if (ABS(im)==1) then
    ylmc= -SQRT(21.d0/(64.d0*pi))*sinth*(5.d0*costh**2-one)*CMPLX(cosphi,sinphi)
   else if (ABS(im)==2) then
    ylmc= SQRT(105.d0/(32.d0*pi))*sinth**2*costh*CMPLX(costwophi,sintwophi)
   else if (ABS(im)==3) then
    ylmc=-SQRT(35.d0/(64.d0*pi))*sinth**3*CMPLX(costhreephi,sinthreephi)
   end if
   case default
   write(msg,'(4a,i6,2a,i6)')ch10,&
&   ' ylmc: ERROR -',ch10,&
&   '  The maximum allowed value for l is,',lx,ch10,&
&   '  however, l=',il
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end select
!
!== Treat the case im < 0 ==
 if (im < 0) then
  ylmc=(-one)**(im)*CONJG(ylmc)
 end if

end function ylmc
!!***
