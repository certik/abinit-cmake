!{\src2tex{textfont=tt}}
!!****f* ABINIT/plm_coeff
!! NAME
!! plm_coeff
!!
!! FUNCTION
!! Compute coefficients depending on Plm and its derivatives where P_lm is a legendre polynome.
!! They are used to compute the second derivatives of spherical harmonics
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mpsang=1+ maximum l quantum number
!!  xx= input value
!!
!! OUTPUT
!!  blm(5,mpsang*mpsang)=coefficients depending on Plm and its derivatives where P_lm is a legendre polynome
!!
!! NOTES
!!
!!
!! PARENTS
!!      initylmg
!!
!! CHILDREN
!!      leave_new,pl_deriv,plm_d2theta,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine plm_coeff(blm,mpsang,xx)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util, except_this_one => plm_coeff
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpsang
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: blm(5,mpsang*mpsang)

!Local variables ---------------------------------------
!scalars
 integer :: il,ilm,ilm1,im
 real(dp) :: sqrx
 character(len=500) :: message
!arrays
 real(dp) :: pl_d2(mpsang),plm_d2t(mpsang*mpsang)

!************************************************************************
 if (abs(xx).gt.1.d0) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' plm_d2theta : ERROR -',ch10,&
&  '   xx > 1 !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

 call plm_d2theta(mpsang,plm_d2t,xx)

 blm=zero
 sqrx=sqrt(abs((one-xx)*(one+xx)))

 do il=0,mpsang-1
  do im=0,il
   ilm=il*il+il+1+im
   ilm1=il*il+il+1-im
   blm(1,ilm)=two*xx*sqrx*(-1)**im*plm_dtheta(il,im,xx)+sqrx*sqrx*plm_d2t(ilm)
   blm(1,ilm1)=blm(1,ilm)
   blm(2,ilm)=(one-two*xx*xx)*(-1)**im*plm_dtheta(il,im,xx)-xx*sqrx*plm_d2t(ilm)
   blm(2,ilm1)=blm(2,ilm)
   blm(3,ilm)=il*(il+1)*(-1)**im*ass_leg_pol(il,im,xx)+plm_d2t(ilm)
   blm(3,ilm1)=blm(3,ilm)
   blm(4,ilm)=-two*xx*sqrx*(-1)**im*plm_dtheta(il,im,xx)+xx*xx*plm_d2t(ilm)
   blm(4,ilm1)=blm(4,ilm)
  end do
 end do

 if(abs(abs(xx)-one)>tol12) then
  do il=1,mpsang-1
   do im=0,il
    ilm=il*il+il+1+im
    if(im==0) then
     blm(5,ilm)=zero
    else
     ilm1=il*il+il+1-im
     blm(5,ilm)=(-1)**im*ass_leg_pol(il,im,xx)/(sqrx*sqrx)-(-1)**im*plm_dtheta(il,im,xx)*xx/sqrx
     blm(5,ilm1)=blm(5,ilm)
    end if
   end do
  end do
 else
  call pl_deriv(mpsang,pl_d2,one)
  do il=1,mpsang-1
   if(il>0) then
    ilm=il*il+il+2
    ilm1=il*il+il
    blm(5,ilm)=-il*(il+1)*ass_leg_pol(il,1,xx)+plm_d2t(ilm)
    blm(5,ilm1)=blm(5,ilm)
   end if
   if(il>1) then
    ilm=il*il+il+3
    ilm1=il*il+il-1
    blm(5,ilm)=-pl_d2(il+1)
    blm(5,ilm1)=blm(5,ilm)
   end if
  end do
 end if
 end subroutine plm_coeff
!!***
