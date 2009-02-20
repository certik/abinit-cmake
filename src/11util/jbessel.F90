!{\src2tex{textfont=tt}}
!!****f* ABINIT/jbessel
!! NAME
!! jbessel
!!
!! FUNCTION
!! Compute spherical Bessel function j_l(x) and derivative(s)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ll=l-order of the Bessel function
!!  order=1 if first derivative is requested
!!        2 if first and second derivatives are requested
!!  xx=where to compute j_l
!!
!! OUTPUT
!!  bes= Bessel function j_l at xx
!!  besp= first derivative of j_l at xx (only if order>=1)
!!  bespp= second derivative of j_l at xx (only if order=2)
!!
!! PARENTS
!!      cutoff_cylinder,m_special_funcs,pawgylm,pawshpfun,shapebes,solvbes
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine jbessel(bes,besp,bespp,ll,order,xx)

 use defs_basis

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll,order
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: bes,besp,bespp

!Local variables ---------------------------------------
!scalars
 integer,parameter :: imax=40
 integer :: ii,il
 real(dp),parameter :: prec=1.d-15
 real(dp) :: besp1,fact,factp,factpp,jn,jnp,jnpp,jr,xx2,xxinv

! *********************************************************************

 if (order>2) stop "Wrong order in jbessel !"

 if (abs(xx)<prec) then
  bes=0.d0;if (ll==0) bes=1.d0
  if (order>=1) then
   besp=0.d0;if (ll==1) besp=1.d0/3.d0
  end if
  if (order==2) then
   bespp=0.d0
   if (ll==0) bespp=-1.d0/3.d0
   if (ll==2) bespp=2.d0/15.d0
  end if
  return
 end if

 xxinv=1.d0/xx

 if (xx<1.d0) then
  xx2=0.5d0*xx*xx
  fact=one
  do il=1,ll
   fact=fact*xx/dble(2*il+1)
  end do
  jn=1.D0;jr=1.D0;ii=0
  do while(abs(jr)>=prec.and.ii<imax)
   ii=ii+1;jr=-jr*xx2/dble(ii*(2*(ll+ii)+1))
   jn=jn+jr
  end do
  bes=jn*fact
  if (abs(jr)>prec) stop 'Error: Bessel function did not converge !'
  if (order>=1) then
   factp=fact*xx/dble(2*ll+3)
   jnp=1.D0;jr=1.D0;ii=0
   do while(abs(jr)>=prec.AND.ii<imax)
    ii=ii+1;jr=-jr*xx2/dble(ii*(2*(ll+ii)+3))
    jnp=jnp+jr
   end do
   besp=-jnp*factp+jn*fact*xxinv*dble(ll)
   if (abs(jr)>prec) stop 'Error: 1st der. of Bessel function did not converge !'
  end if
  if (order==2) then
   factpp=factp*xx/dble(2*ll+5)
   jnpp=1.D0;jr=1.D0;ii=0
   do while(abs(jr)>=prec.AND.ii<imax)
    ii=ii+1;jr=-jr*xx2/dble(ii*(2*(ll+ii)+5))
    jnpp=jnpp+jr
   end do
   besp1=-jnpp*factpp+jnp*factp*xxinv*dble(ll+1)
   if (abs(jr)>prec) stop 'Error: 2nd der. of Bessel function did not converge !'
  end if
 else
  jn =sin(xx)*xxinv
  jnp=(-cos(xx)+jn)*xxinv
  do il=2,ll+1
   jr=-jn+dble(2*il-1)*jnp*xxinv
   jn=jnp;jnp=jr
  end do
  bes=jn
  if (order>=1) besp =-jnp+jn *xxinv*dble(ll)
  if (order==2) besp1= jn -jnp*xxinv*dble(ll+2)
 end if

 if (order==2) bespp=-besp1+besp*ll*xxinv-bes*ll*xxinv*xxinv

end subroutine jbessel
!!***
