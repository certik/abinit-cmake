!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_special_funcs
!! NAME
!! m_special_funcs
!!
!! FUNCTION
!! This module contains routines and functions used to 
!! evaluate special functions frequently occuring in Abinit
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!@ ABIDOC
module m_special_funcs

 use defs_basis
 use m_errors

 implicit none

! === List of available public routines and functions ===
 public ::           &
&  jbessel_4spline      ! Compute spherical Bessel functions and derivatives employing a polynomial approximation for q->0
!@END ABIDOC

CONTAINS  !===========================================================

!!***

!------------------------------------------------------------------------

!!****f* m_special_funcs/jbessel_4spline
!! NAME
!!  jbessel_4spline
!!
!! FUNCTION
!!  Compute spherical Bessel functions and derivatives. 
!!  A polynomial approximation is employed for q-->0.
!!  
!! INPUT
!!  ll=l-order of the Bessel function
!!  tol=tolerance below which a Polynomial approximation is employed
!!   both for jl and its derivative (if required)
!!  order=1 if only first derivative is requested
!!        2 if first and second derivatives are requested
!!  xx=where to compute j_l
!!
!! OUTPUT
!!  bes=Spherical Bessel function j_l at xx
!!  besp= first derivative of j_l at xx (only if order>=1)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      paw_mkrhox_spl,psp7nl
!!
!! CHILDREN
!!      jbessel,wrtout
!!
!! SOURCE

subroutine jbessel_4spline(bes,besp,ll,order,xx,tol)
!Arguments ---------------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 integer,intent(in) :: ll,order
 real(dp),intent(in) :: xx,tol
 real(dp),intent(out) :: bes,besp

!Local variables ---------------------------------------
!scalars
 real(dp) :: bespp
 character(len=500) :: msg
 real(dp) :: arg,bes0a,bes0ap,bes0b,bes0bp,bes1a,bes1ap,bes1b,bes1bp
 real(dp) :: bes2a,bes2ap,bes2b,bes2bp,bes3a,bes3ap,bes3b,bes3bp
! *********************************************************************

!=== l=0,1,2 and 3 spherical Bessel functions (and derivatives) ===
 bes0a(arg)=1.0_dp-arg**2/6.0_dp*(1.0_dp-arg**2/20.0_dp)
 bes0b(arg)=sin(arg)/arg
 bes1a(arg)=(10.0_dp-arg*arg)*arg/30.0_dp
 bes1b(arg)=(sin(arg)-arg*cos(arg))/arg**2
 bes2a(arg)=arg*arg/15.0_dp-arg**4/210.0_dp
 bes2b(arg)=((3.0_dp-arg**2)*sin(arg)-3.0_dp*arg*cos(arg))/arg**3
 bes3a(arg)=arg*arg*arg/105.0_dp-arg**5/1890.0_dp+arg**7/83160.0_dp
 bes3b(arg)=(15.0_dp*sin(arg)-15.0_dp*arg*cos(arg)-6.0_dp*arg**2*sin(arg)+arg**3*cos(arg))/arg**4
 bes0ap(arg)=(-10.0_dp+arg*arg)*arg/30.0_dp
 bes0bp(arg)=-(sin(arg)-arg*cos(arg))/arg**2
 bes1ap(arg)=(10.0_dp-3.0_dp*arg*arg)/30.0_dp
 bes1bp(arg)=((arg*arg-2.0_dp)*sin(arg)+2.0_dp*arg*cos(arg))/arg**3
 bes2ap(arg)=(1.0_dp-arg*arg/7.0_dp)*2.0_dp*arg/15.0_dp
 bes2bp(arg)=((4.0_dp*arg*arg-9.0_dp)*sin(arg)+(9.0_dp-arg*arg)*arg*cos(arg))/arg**4
 bes3ap(arg)=(1.0_dp/35-arg*arg/378.0_dp+arg**4/11880.0_dp)*arg*arg
 bes3bp(arg)=((-60.0_dp+27.0_dp*arg*arg-arg**4)*sin(arg)+(60.0_dp*arg-7.0_dp*arg**3)*cos(arg))/arg**5

 if (order>2) stop "Wrong order in jbessel !"

 select case (ll)

 case (0)
  if (xx<TOL) then
   bes=bes0a(xx)
   if (order>=1) besp=bes0ap(xx)
  else
   bes=bes0b(xx)
   if (order>=1) besp=bes0bp(xx)
  end if

 case (1)
  if (xx<TOL) then
   bes=bes1a(xx)
   if (order>=1) besp=bes1ap(xx)
  else
   bes=bes1b(xx)
   if (order>=1) besp=bes1bp(xx)
  end if

 case (2)
  if (xx<TOL) then
   bes=bes2a(xx)
   if (order>=1) besp=bes2ap(xx)
  else
   bes=bes2b(xx)
   if (order>=1) besp=bes2bp(xx)
  end if

 case (3)
  if (xx<TOL) then
   bes=bes3a(xx)
   if (order>=1) besp=bes3ap(xx)
  else
   bes=bes3b(xx)
   if (order>=1) besp=bes3bp(xx)
  end if

 case (4:)
  call jbessel(bes,besp,bespp,ll,order,xx)

 case default
  write(msg,'(4a,i4)')ch10,&
&  ' jbessel_4spline : BUG - ',ch10,&
&  ' wrong value for ll = ',ll
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end select

end subroutine jbessel_4spline

end module m_special_funcs
!!***
