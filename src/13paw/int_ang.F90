!{\src2tex{textfont=tt}}
!!****f* ABINIT/int_ang
!! NAME
!! int_ang
!!
!! FUNCTION
!! routine for evaluation of angular part for <phi_i|nabla|phi_j> and <tphi_i|nabla|tphi_j>
!! COPYRIGHT
!! Copyright (C) 2005-2006 ABINIT group (VR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  mpsang=1+ max. angular momentum
!!
!!
!! OUTPUT
!!  ang_phipphj :: angular part for <phi_i|nabla|phi_j> and <tphi_i|nabla|tphi_j>
!!  ang_phipphj(i,j,1)=\int sin\theta cos\phi Si Sj d\omega
!!  ang_phipphj(i,j,2)=\int cos\theta cos\phi Si \frac{d}{d\theta}Sj d\Omega
!!  ang_phipphj(i,j,3)=\int -sin\phi  Si \frac{d}{d\phi}Sj d\Omega
!!  ang_phipphj(i,j,4)=\int sin\theta sin\phi Si Sj d\Omega
!!  ang_phipphj(i,j,5)=\int cos\theta sin\phi Si \frac{d}{d\theta}Sj d\Omega
!!  ang_phipphj(i,j,6)=\int cos\phi Si \frac{d}{d\phi}Sj d\Omega
!!  ang_phipphj(i,j,7)=\int cos\theta  Si Sj d\Omega
!!  ang_phipphj(i,j,8)=\int -sin\theta Si \frac{d}{d\theta}Sj d\Omega
!!
!! SIDE EFFECTS
!!  (only writing, printing)
!!
!! NOTES
!!
!! PARENTS
!!      optics_paw,pawnabla_init
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine int_ang(ang_phipphj,mpsang)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang
!arrays
 real(dp),intent(out) :: ang_phipphj(mpsang**2,mpsang**2,8)

!Local variables-------------------------------
 character(len=500) :: message
 real(dp) :: ang_phipphj_tmp(9,9,8)

! ************************************************************************
!
!DEBUG
! write(*,*)' int_ang : enter'
!ENDDEBUG
 if (mpsang>3) then
  write(message, '(4a)' )ch10,&
&  ' int_ang :  -',ch10,&
&  '  Not designed for angular momentum greater than 2 !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 ang_phipphj_tmp=zero
!
 ang_phipphj_tmp(1,2,4)=one/sqrt(three)
 ang_phipphj_tmp(1,3,7)=one/sqrt(three)
 ang_phipphj_tmp(1,4,1)=one/sqrt(three)
 ang_phipphj_tmp(2,1,4)=one/sqrt(three)
 ang_phipphj_tmp(2,5,1)=one/sqrt(five)
 ang_phipphj_tmp(2,6,7)=one/sqrt(five)
 ang_phipphj_tmp(2,7,4)=-one/sqrt(15._dp)
 ang_phipphj_tmp(2,9,4)=-one/sqrt(five)
 ang_phipphj_tmp(3,1,7)=one/sqrt(three)
 ang_phipphj_tmp(3,6,4)=one/sqrt(five)
 ang_phipphj_tmp(3,7,7)=two/sqrt(15._dp)
 ang_phipphj_tmp(3,8,1)=one/sqrt(five)
 ang_phipphj_tmp(4,1,1)=one/sqrt(three)
 ang_phipphj_tmp(4,5,4)=one/sqrt(five)
 ang_phipphj_tmp(4,7,1)=-one/sqrt(15._dp)
 ang_phipphj_tmp(4,8,7)=one/sqrt(five)
 ang_phipphj_tmp(4,9,1)=one/sqrt(five)
 ang_phipphj_tmp(5,2,1)=one/sqrt(five)
 ang_phipphj_tmp(5,4,4)=one/sqrt(five)
 ang_phipphj_tmp(6,2,7)=one/sqrt(five)
 ang_phipphj_tmp(6,3,4)=one/sqrt(five)
 ang_phipphj_tmp(7,2,4)=-one/sqrt(15._dp)
 ang_phipphj_tmp(7,3,7)=two/sqrt(15._dp)
 ang_phipphj_tmp(7,4,1)=-one/sqrt(15._dp)
 ang_phipphj_tmp(8,3,1)=one/sqrt(five)
 ang_phipphj_tmp(8,4,7)=one/sqrt(five)
 ang_phipphj_tmp(9,2,4)=-one/sqrt(five)
 ang_phipphj_tmp(9,4,1)=one/sqrt(five)
!
 ang_phipphj_tmp(1,2,5)=one/two/sqrt(three)
 ang_phipphj_tmp(1,3,8)=two/sqrt(three)
 ang_phipphj_tmp(1,4,2)=one/two/sqrt(three)
 ang_phipphj_tmp(2,5,2)=one/two/sqrt(five)
 ang_phipphj_tmp(2,6,8)=three/sqrt(five)
 ang_phipphj_tmp(2,7,5)=-sqrt(three/five)
 ang_phipphj_tmp(2,9,5)=-one/two/sqrt(five)
 ang_phipphj_tmp(3,6,5)=one/two/sqrt(five)
 ang_phipphj_tmp(3,7,8)=two/sqrt(three)
 ang_phipphj_tmp(3,8,2)=one/two/sqrt(five)
 ang_phipphj_tmp(4,5,5)=one/two/sqrt(five)
 ang_phipphj_tmp(4,7,2)=-sqrt(three/five)
 ang_phipphj_tmp(4,8,8)=three/sqrt(five)
 ang_phipphj_tmp(4,9,2)=one/two/sqrt(five)
 ang_phipphj_tmp(5,2,2)=one/four/sqrt(five)
 ang_phipphj_tmp(5,4,5)=one/four/sqrt(five)
 ang_phipphj_tmp(6,2,8)=-one/sqrt(five)
 ang_phipphj_tmp(6,3,5)=-one/sqrt(five)
 ang_phipphj_tmp(7,2,5)=one/sqrt(15.0_dp)
 ang_phipphj_tmp(7,3,8)=-two/sqrt(15.0_dp)
 ang_phipphj_tmp(7,4,2)=one/sqrt(15.0_dp)
 ang_phipphj_tmp(8,3,2)=-one/sqrt(five)
 ang_phipphj_tmp(8,4,8)=one/sqrt(five)
 ang_phipphj_tmp(9,2,5)=-one/four/sqrt(five)
 ang_phipphj_tmp(9,4,2)=one/four/sqrt(five)
!
 ang_phipphj_tmp(1,2,6)=sqrt(three)/two
 ang_phipphj_tmp(1,4,3)=sqrt(three)/two
 ang_phipphj_tmp(2,5,3)=sqrt(five)/two
 ang_phipphj_tmp(2,9,6)=-sqrt(five)/two
 ang_phipphj_tmp(3,6,6)=sqrt(five)/two
 ang_phipphj_tmp(3,8,3)=sqrt(five)/two
 ang_phipphj_tmp(4,5,6)=sqrt(five)/two
 ang_phipphj_tmp(4,9,3)=sqrt(five)/two
 ang_phipphj_tmp(5,2,3)=-sqrt(five)/four
 ang_phipphj_tmp(5,4,6)=-sqrt(five)/four
 ang_phipphj_tmp(9,2,6)=sqrt(five)/four
 ang_phipphj_tmp(9,4,3)=-sqrt(five)/four

 ang_phipphj(:,:,:)=ang_phipphj_tmp(1:mpsang**2,1:mpsang**2,:)

!DEBUG
!write(6,*)' int_ang : exit '
!stop
!ENDDEBUG

 end subroutine int_ang
!!***
