!{\src2tex{textfont=tt}}
!!****f* ABINIT/order_fs_kpts
!!
!! NAME
!! order_fs_kpts
!!
!! FUNCTION
!! This routine re-orders the kpoints on the standard grid which belong
!!  to the Fermi surface: put them in increasing z, then y,  then x
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   FSkptflag = flag for inclusion of kpoint in fermi-surface set
!!   hdr = GS header containing dimensions
!!   nFSkptirred = number of irreducible FS kpoints
!!
!! OUTPUT
!!   FSirredtoGS = mapping of irreducible kpoints to GS set
!!   FSkptirrank = rank of kpoint based on its coordinates
!!   FSkptirred = irreducible FS kpoint coordinates
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      canon9
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine order_fs_kpts(FSkptflag,FSkptirrank,FSkptirred,FSirredtoGS,hdr,nFSkptirred)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nFSkptirred
 type(hdr_type),intent(in) :: hdr
!arrays
 integer,intent(in) :: FSkptflag(hdr%nkpt)
 integer,intent(out) :: FSirredtoGS(nFSkptirred),FSkptirrank(nFSkptirred)
 real(dp),intent(inout) :: FSkptirred(3,nFSkptirred)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,ikpt,jFSkpt,kFSkpt,new
 real(dp) :: res
!arrays
 integer :: kptrank(hdr%nkpt)
 real(dp) :: kpt(3)

! *************************************************************************

!rank is used to order kpoints
!write (*,*) 'order_fs_kpt : ikpt,kpt,kptrank = '
 kptrank(:) = 0
 do ikpt=1,hdr%nkpt
  call canon9(hdr%kptns(1,ikpt),kpt(1),res)
  call canon9(hdr%kptns(2,ikpt),kpt(2),res)
  call canon9(hdr%kptns(3,ikpt),kpt(3),res)

  kptrank(ikpt) = 100000000.0_dp*(kpt(1)+one) + &
&  100000.0_dp*(kpt(2)+one) + &
&  100.0_dp*(kpt(3)+one)
! kptrank(ikpt) = 2*hdr%nkpt - ikpt + 1

! DEBUG
! write (*,*) ikpt,kpt,kptrank(ikpt)
! ENDDEBUG
 end do
!DEBUG
!write (*,*) 'order_fs_kpt : kptrank = '
!write (*,'(6I10)') kptrank
!ENDDEBUG

 iFSkpt=1
 do ikpt=1,hdr%nkpt
! is ikpt in fermi surface?
  if (FSkptflag(ikpt) == 1) then
!  add kpt to FS kpts, in order, increasing z, then y, then x !
   new = 1
!  look for position to insert kpt ikpt among irredkpts already found
   do jFSkpt=1,iFSkpt-1
    if (FSkptirrank(jFSkpt) > kptrank(ikpt)) then
!    shift all the others up
     do kFSkpt=iFSkpt-1,jFSkpt,-1
      FSkptirred(:,kFSkpt+1) = FSkptirred(:,kFSkpt)
      FSkptirrank(kFSkpt+1) = FSkptirrank(kFSkpt)
      FSirredtoGS(kFSkpt+1) = FSirredtoGS(kFSkpt)
     end do
!    insert kpoint ikpt
     call canon9(hdr%kptns(1,ikpt),FSkptirred(1,jFSkpt),res)
     call canon9(hdr%kptns(2,ikpt),FSkptirred(2,jFSkpt),res)
     call canon9(hdr%kptns(3,ikpt),FSkptirred(3,jFSkpt),res)

     FSkptirrank(jFSkpt) = kptrank(ikpt)
     FSirredtoGS(jFSkpt) = ikpt
     new=0
     exit
    end if
   end do
!  ikpt in FS, but not counted yet and higher rank than all previous
   if (new == 1) then
    call canon9(hdr%kptns(1,ikpt),FSkptirred(1,iFSkpt),res)
    call canon9(hdr%kptns(2,ikpt),FSkptirred(2,iFSkpt),res)
    call canon9(hdr%kptns(3,ikpt),FSkptirred(3,iFSkpt),res)
    FSkptirrank(iFSkpt) = kptrank(ikpt)
    FSirredtoGS(iFSkpt) = ikpt
   end if
   iFSkpt=iFSkpt+1
  end if
 end do

!DEBUG
!do ikpt=1,hdr%nkpt
!call canon9(hdr%kptns(1,ikpt),FSkptirred(1,ikpt),res)
!call canon9(hdr%kptns(2,ikpt),FSkptirred(2,ikpt),res)
!call canon9(hdr%kptns(3,ikpt),FSkptirred(3,ikpt),res)
!FSirredtoGS(ikpt) = ikpt
!end do
!ENDDEBUG

end subroutine order_fs_kpts
!!***
