!{\src2tex{textfont=tt}}
!!****f* ABINIT/leave_myproc
!! NAME
!! leave_myproc
!!
!! FUNCTION
!! Routine for clean exit of f90 code by one processor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, NCJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (nothing)
!!
!! OUTPUT
!!  (only writing, then stop)
!!
!! NOTES
!! By default, it uses "call exit(1)", that is not completely
!! portable.
!!
!! PARENTS
!!      leave_new,leave_test
!!
!! CHILDREN
!!      exit
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine leave_myproc

 use defs_basis
!no_abirules
#if defined FC_NAG
 use f90_unix
#endif

 implicit none

!Arguments ------------------------------------

!Local variables-------------------------------
!scalars
 logical :: testopen

! **********************************************************************

 inquire(ab_out,OPENED=testopen)
 if(testopen)close(ab_out)

#if defined FC_NAG
 call exit(-1)
#elif defined HAVE_FORTRAN_EXIT
 call exit(1)
#else
 stop 1
#endif

end subroutine leave_myproc
!!***
