!{\src2tex{textfont=tt}}
!!****f* ABINIT/sbf8
!! NAME
!! sbf8
!!
!! FUNCTION
!! Computes set of spherical bessel functions using accurate algorithm
!! based on downward recursion in order and normalization sum.
!! Power series used at small arguments.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nm=maximum angular momentum wanted + one
!!  xx=argument of sbf
!!
!! OUTPUT
!!  sb_out(nm)=values of spherical bessel functions for l=0,nm-1
!!
!! PARENTS
!!      psp8nl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine sbf8(nm,xx,sb_out)

 use defs_basis

 implicit none

!Arguments----------------------------------------------------------
!scalars
 integer,intent(in) :: nm
 real(dp),intent(in) :: xx
!arrays
 real(dp),intent(out) :: sb_out(nm)

!Local variables-------------------------------
!scalars
 integer :: nlim,nn
 real(dp) :: fn,sn,xi,xn,xs
!arrays
 real(dp),allocatable :: sb(:)

! *************************************************************************

 if(xx<= 1.0e-36_dp) then
! zero argument section
  sb_out(:)=zero
  sb_out(1)=one
 else if(xx<1.e-3_dp) then
! small argument section
  xn=one
  xs=half*xx**2
  do nn=1,nm
   sb_out(nn)=xn*(one - xs*(one - xs/(4*nn+6))/(2*nn+1))
   xn=xx*xn/(2*nn+1)
  end do
 else
! recursion method
  if(xx<one) then
   nlim=nm+int(15.0e0_dp*xx)+1
  else
   nlim=nm+int(1.36e0_dp*xx)+15
  end if
  allocate(sb(nlim+1))
  nn=nlim
  xi=one/xx
  sb(nn+1)=zero
  sb(nn)=1.e-18_dp
  sn=dble(2*nn-1)*1.e-36_dp
  do nn=nlim-1,1,-1
   sb(nn)=dble(2*nn+1)*xi*sb(nn+1) - sb(nn+2)
  end do
  do nn=1,nlim-1
   sn=sn + dble(2*nn-1)*sb(nn)*sb(nn)
  end do
  fn=1.d0/sqrt(sn)
  sb_out(:)=fn*sb(1:nm)
  deallocate(sb)
 end if
end subroutine sbf8
!!***
