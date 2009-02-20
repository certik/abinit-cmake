!{\src2tex{textfont=tt}}
!!****f* ABINIT/prmat
!!
!! NAME
!! prmat
!!
!! FUNCTION
!! This subroutine prints real*8 matrices in an attractive format.
!!
!! COPYRIGHT
!! Copyright (C) 1987-2008 ABINIT group (ZL, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mat(mi,nj)= matrix to be printed
!!  mi
!!  ni
!!  nj
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      newkpt
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prmat (mat, ni, nj, mi)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mi,ni,nj
!arrays
 real(dp),intent(in) :: mat(mi,nj)

!Local variables-------------------------------
!scalars
 integer,parameter :: nline=10
 integer :: ii,jj,jstart,jstop

! *************************************************************************

 do  jstart = 1, nj, nline
  jstop = min(nj, jstart+nline-1)
  write (6, '(3x,10(i4,8x))' ) (jj,jj=jstart,jstop)
 end do

 do ii = 1,ni
  do jstart= 1, nj, nline
   jstop = min(nj, jstart+nline-1)
   if (jstart==1) then
    write (6, '(i3,1p,10e12.4)' ) ii, (mat(ii,jj),jj=jstart,jstop)
   else
    write (6, '(3x,1p,10e12.4)' )    (mat(ii,jj),jj=jstart,jstop)
   end if
  end do
 end do
!
end subroutine prmat
!!***
