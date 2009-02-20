!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkfsqgrid
!!
!! NAME
!! mkfsqgrid
!!
!! FUNCTION
!! This routine sets up the qpoints between the full FS kpt grid points
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  FSkpt = kpoints on the FS
!!  nFSkpt = number of kpts on the FS
!!
!! OUTPUT
!!  nFSqpt = full number of qpoints between points on the FS
!!  FStoqpt = qpoint index for each pair of kpt on the FS
!!  tmpFSqpt = temp array with coordinates of the
!!    qpoints between points on the FS
!!
!! NOTES
!!   might need the inverse indexing qpt -> (kpt1,kpt2) (kpt3,kpt4) ...
!!
!! PARENTS
!!
!! CHILDREN
!!      canon9
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkfsqgrid(FSkpt,FStoqpt,nFSkpt,nFSqpt,tmpFSqpt)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nFSkpt
 integer,intent(out) :: nFSqpt
!arrays
 integer,intent(out) :: FStoqpt(nFSkpt,nFSkpt)
 real(dp),intent(in) :: FSkpt(3,nFSkpt)
 real(dp),intent(out) :: tmpFSqpt(3,nFSkpt*nFSkpt)

!Local variables-------------------------------
!scalars
 integer :: iFSqpt,ikpt1,ikpt2,new
 real(dp) :: shift,ss
!arrays
 real(dp) :: k1(3),kpt(3)

! *************************************************************************

 nFSqpt=0
 tmpFSqpt(:,:)=zero
 do ikpt1=1,nFSkpt
  do ikpt2=1,nFSkpt
   k1(:) = FSkpt(:,ikpt1) - FSkpt(:,ikpt2)
   call canon9(k1(1),kpt(1),shift)
   call canon9(k1(2),kpt(2),shift)
   call canon9(k1(3),kpt(3),shift)

   new=1
!  is kpt among the FS qpts found already?
   do iFSqpt=1,nFSqpt
    ss=(kpt(1)-tmpFSqpt(1,iFSqpt))**2 + &
&    (kpt(2)-tmpFSqpt(2,iFSqpt))**2 + &
&    (kpt(3)-tmpFSqpt(3,iFSqpt))**2
    if (ss < tol6) then
     FStoqpt(ikpt1,ikpt2) = iFSqpt
     new=0
     exit
    end if
   end do
   if (new == 1) then
    nFSqpt=nFSqpt+1
    tmpFSqpt(:,nFSqpt) = kpt(:)
    FStoqpt(ikpt1,ikpt2) = nFSqpt
   end if

  end do
 end do

!got nFSqpt,tmpFSqpt,FStoqpt

end subroutine mkfsqgrid
!!***
