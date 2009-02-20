!{\src2tex{textfont=tt}}
!!****f* ABINIT/mvrecord
!! NAME
!! mvrecord
!!
!! FUNCTION
!! This subroutine moves forward or backward in a Fortran binary file by nn records.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (ZL,DCA,XG,GMR,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nrec=number of records
!! unitfile= file unit number
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!! For the future : one should treat the possible errors of backspace
!!
!! PARENTS
!!      WffReadSkipRec
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mvrecord(ierr,nrec,unitfile)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nrec,unitfile
 integer,intent(out) :: ierr

!Local variables-------------------------------
!scalars
 integer :: irec

! *************************************************************************

 if ( nrec > 0) then
! Move forward   nrec records
  do irec=1,nrec
   read(unitfile,iostat=ierr)
  end do
 else
! Move backward   nrec records
  do irec=1,-nrec
   backspace (unit=unitfile,iostat=ierr)
  end do
 end if

end subroutine mvrecord
!!***
