!{\src2tex{textfont=tt}}
!!****f* ABINIT/pl_deriv
!! NAME
!! pl_deriv
!!
!! FUNCTION
!! Compute d2(Pl (x)))/d(x)2  where P_l is a legendre polynomial
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
!!  pl_d2(mpsang*mpsang)
!!
!! NOTES
!!
!!
!! PARENTS
!!      plm_coeff
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pl_deriv(mpsang,pl_d2,xx)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mpsang
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: pl_d2(mpsang)

!Local variables ---------------------------------------
!scalars
 integer :: il,ilm
 character(len=500) :: message
!arrays
 real(dp) :: pl(mpsang),pl_d1(mpsang)

! *********************************************************************

 if (abs(xx).gt.1.d0) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pl_deriv : ERROR -',ch10,&
&  '   xx > 1 !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if


 pl_d2=zero; pl_d1=zero; pl=zero
 pl(1)=one; pl(2)=xx
 pl_d1(1)=zero; pl_d1(2)=one
 pl_d2(1)=zero; pl_d2(2)=zero
 if (mpsang>1) then
! terme Pl(x)
  do il=1,mpsang-1
   ilm=il+1
   pl(ilm)=((2*il-1)*xx*pl(ilm-1)-(il-1)*pl(ilm-2))/dfloat(il)
   pl_d1(ilm)=((2*il-1)*(xx*pl_d1(ilm-1)+pl(ilm-1))-(il-1)*pl_d1(ilm-2))/dfloat(il)
  end do
! terme d(Plm)/dx
  do il=1,mpsang-1
   ilm=il+1
   pl_d1(ilm)=((2*il-1)*(xx*pl_d1(ilm-1)+pl(ilm-1))-(il-1)*pl_d1(ilm-2))/dfloat(il)
  end do
! terme d2(Plm)/dx2
  do il=1,mpsang-1
   ilm=il+1
   pl_d2(ilm)=((2*il-1)*(xx*pl_d2(ilm-1)+2.d0*pl_d1(ilm-1))-(il-1)*pl_d2(ilm-2))/dfloat(il)
  end do
 end if

 end subroutine pl_deriv
!!***
