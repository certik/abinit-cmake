!{\src2tex{textfont=tt}}
!!****f* ABINIT/findk
!! NAME
!! findk
!!
!! FUNCTION
!! Check whether the k-point is in the set of the kbz
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kbz(3,nkbz)=coordinates of k points in the BZ
!!  nkcalc= number of k points for GW calculation (input variable)
!!  nkbz=number of k points in Brillouin zone
!!  xkcalc(3,nkcalc)= coordinates of the k points
!!  umklp_opt=0 if no umklapp vector is admitted, 1 otherwise
!!
!! OUTPUT
!!  kcalc=index of the k points inside kbz
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

subroutine findk(nkcalc,nkbz,xkcalc,kbz,kcalc,umklp_opt)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkbz,nkcalc,umklp_opt
!arrays
 integer,intent(out) :: kcalc(nkcalc)
 real(dp),intent(in) :: kbz(3,nkbz),xkcalc(3,nkcalc)

!Local variables-------------------------------
!scalars
 integer :: ik,jj
 real(dp) :: shift
 character(len=500) :: message
!arrays
 real(dp) :: dummy(3),ktest(3)

! *************************************************************************

 write(message,'(2a)')ch10,' findk : check if the k-points for sigma are in the set of BZ'
 call wrtout(06,message,'COLL')

 if (umklp_opt==0) then 
  do jj=1,nkcalc
   kcalc(jj)=0
   do ik=1,nkbz
    if (all(abs(xkcalc(:,jj)-kbz(:,ik))<1.e-3)) kcalc(jj)=ik
   end do
   if (kcalc(jj)==0) then
    write(message,'(4a,3(f6.3,1x),a)')ch10,' findk : ERROR -',&
&    ch10,' k-point ',xkcalc(:,jj),' not in the set of kbz'
    call wrtout(06,message,'COLL') ; call leave_new('COLL')
   end if
  end do
 else if (umklp_opt==1) then 
  do jj=1,nkcalc
   kcalc(jj)=0
   do ik=1,nkbz
    dummy=xkcalc(:,jj)-kbz(:,ik)
    call canon9(dummy(1),ktest(1),shift)
    call canon9(dummy(2),ktest(2),shift)
    call canon9(dummy(3),ktest(3),shift)
    if (all(abs(ktest)<1.e-3)) kcalc(jj)=ik
   end do
   if(kcalc(jj)==0) then
    write(message,'(4a,3(f6.3,1x),3a)')ch10,' findk : ERROR -',&
&    ch10,' k-point ',xkcalc(:,jj),' not in the set of kbz',&
&    ch10,' even though umklapp G0 vectors are allowed '
    call wrtout(06,message,'COLL') ; call leave_new('COLL')
   end if
  end do
 else 
  write(message,'(4a)')ch10,' findk : BUG-',ch10,' wrong value for uklp'
  call wrtout(06,message,'COLL') ; call leave_new('COLL')
 end if 

 write(message,'(2a)')' check ok !',ch10
 call wrtout(06,message,'COLL')

end subroutine findk
!!***
