!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkkptrank
!!
!! NAME
!! mkkptrank
!!
!! FUNCTION
!! This routine sets up the kpt ranks for comparing kpts
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt = number of kpoints
!!  kpt = coordinates of kpoints
!!
!! OUTPUT
!!  rank = rank of each kpoint
!!  invrank = inverse of rank of each kpoint (recuperate ikpt)
!!
!! NOTES
!!  WARNING: supposes kpt grid is coarser than 100 points/direction
!!    in the BZ
!!
!! PARENTS
!!      bfactor,mkfskgrid,mkqptequiv
!!
!! CHILDREN
!!      canon9
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkkptrank (kpt,nkpt,rank,invrank)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt
!arrays
 integer,intent(out) :: invrank(16000000),rank(nkpt)
 real(dp),intent(in) :: kpt(3,nkpt)

!Local variables -------------------------
!scalars
 integer :: ikpt
 real(dp) :: res
!arrays
 real(dp) :: redkpt(3)

! *********************************************************************

 invrank(:) = -1

!Ensure kpt(i)+one is positive, and the smallest
!difference between kpts should be larger than 1/100
!ie ngkpt < 100.
 do ikpt=1,nkpt
  call canon9(kpt(1,ikpt),redkpt(1),res)
  call canon9(kpt(2,ikpt),redkpt(2),res)
  call canon9(kpt(3,ikpt),redkpt(3),res)
  rank(ikpt) = int(8000000.0_dp*(redkpt(1)+half+tol8) + &
&  40000.0_dp*(redkpt(2)+half+tol8) + &
&  200.0_dp*(redkpt(3)+half+tol8))
  if (rank(ikpt) > 16000000) then
   write (*,*) ' mkkptrank : error : rank should be inferior to ', 2000000
   stop
  end if
  invrank(rank(ikpt)) = ikpt
 end do


end subroutine mkkptrank
!!***
