!{\src2tex{textfont=tt}}
!!****f* ABINIT/mk_hdr_check_fmt
!! NAME
!! mk_hdr_check_fmt
!!
!! FUNCTION
!! make a format needed in hdr_check, for arrays of nint integers
!! each of format i3
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (the_author missing ...)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nelm=number of elements to be printed
!!
!! OUTPUT
!!  character(len=26), typfmt= format needed
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      hdr_check
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mk_hdr_check_fmt(nelm,typfmt)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nelm
 character(len=26),intent(out) :: typfmt

!Local variables-------------------------------
 character(len=1), parameter :: number(0:10)=(/'0','1','2','3','4','5','6','7','8','9',' '/)
 character(len=26), parameter :: templatefmt='(2x,  i3,t41   ,a,2x,  i3)'
!character(len=500) :: message                          ! to be uncommented, eventually
!scalars
 integer :: ii,istop,mu

! *************************************************************************

!DEBUG
!write(6,*)' mk_hdr_check_fmt : enter, nelm= ',nelm
!ENDDEBUG

!DEBUG                                             ! to be uncommented, eventually
!write(message,'(a,a,a,a,a,a,i6)') ch10,&
!& ' mk_hdr_check_fmt: BUG -',ch10,&
!& '  The argument sizein should be a positive number,',ch10,&
!& '  however, sizein=',sizein
!call wrtout(06,message,'COLL')
!call leave_new('COLL')
!ENDDEBUG

!Initialize the format
 typfmt=templatefmt

!Generate the type format specifier
 ii=nelm/10
 if ( ii /= 0 ) then
  typfmt(5:5) = number(ii)
  typfmt(22:22) = number(ii)
 else
  typfmt(5:5) = ' '
  typfmt(22:22) = ' '
 end if
 ii = nelm - 10 * (nelm/10)
 typfmt(6:6) = number(ii)
 typfmt(23:23) = number(ii)

!DEBUG
!write(6,*)' mk_hdr_check_fmt : exit'
!write(6,*)' mk_hdr_check_fmt : typfmt="',typfmt,'"'
!stop
!ENDDEBUG

end subroutine mk_hdr_check_fmt
!!***
