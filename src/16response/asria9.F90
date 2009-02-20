!{\src2tex{textfont=tt}}
!!****f* ABINIT/asria9
!! NAME
!! asria9
!!
!! FUNCTION
!! Imposition of the Acoustic sum rule on the InterAtomic Forces
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! asr=(0 => no ASR, 1 or 2=> the diagonal element is modified to give the ASR)
!! asrflg=(1 => the correction to enforce asr is computed from
!!           d2cart, but NOT applied;
!!         2 => one uses the previously determined correction)
!! mpert =maximum number of ipert
!! natom=number of atom
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!! d2cart=matrix of second derivatives of total energy, in cartesian coordinates
!! d2asr=matrix used to store the correction needed to fulfill
!! the acoustic sum rule.
!!
!! PARENTS
!!      anaddb,gath3
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine asria9(asr,asrflg,d2asr,d2cart,mpert,natom)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: asr,asrflg,mpert,natom
!arrays
 real(dp),intent(inout) :: d2asr(2,3,3,natom),d2cart(2,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: idir1,idir2,ii,ipert1,ipert2
 character(len=500) :: message

! *********************************************************************

 if(asrflg/=1 .and. asrflg/=2)then
  write(message, '(a,a,a,a,a,i4,a,a,a)' )&
&  ' asria9 : ERROR -',ch10,&
&  '  asrflg should be equal to 1 or 2 while',ch10,&
&  '  it is equal to ',asrflg,'.',ch10,&
&  '  Action : change asrflg in your input file.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (asr==1.or.asr==2)then
  write(message, '(a)' )&
&  ' asria9 : imposition of the ASR for the interatomic forces.'
  call wrtout(6,message,'COLL')
  do ii=1,2
   do ipert1=1,natom
    do idir1=1,3
     do idir2=1,3

!     Compute d2asr
      if(asrflg==1)then
       d2asr(ii,idir1,idir2,ipert1)=0.0_dp
       do ipert2=1,natom
        d2asr(ii,idir1,idir2,ipert1)=&
&        d2asr(ii,idir1,idir2,ipert1)+&
&        d2cart(ii,idir1,ipert1,idir2,ipert2)
       end do
      end if

!     Apply d2asr
      if(asrflg==2)then
       d2cart(ii,idir1,ipert1,idir2,ipert1)=&
&       d2cart(ii,idir1,ipert1,idir2,ipert1)-&
&       d2asr(ii,idir1,idir2,ipert1)
      end if

     end do
    end do
   end do
  end do
 end if

end subroutine asria9
!!***
