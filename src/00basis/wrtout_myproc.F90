!{\src2tex{textfont=tt}}
!!****f* ABINIT/wrtout_myproc
!! NAME
!! wrtout_myproc
!!
!! FUNCTION
!! Do the output for one proc. For parallel or sequential output use wrtout()
!! instead. Also allows to treat correctly the write operations for
!! Unix (+DOS) and MacOS.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  unit=unit number for writing
!!
!! OUTPUT
!!  (only writing)
!!
!! SIDE EFFECTS
!!  message=(character(len=500)) message to be written
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wrtout_myproc(unit,message)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=500),intent(inout) :: message

!Local variables-------------------------------
!scalars
 integer,save :: iexit=0,ncomment=0,nwarning=0
 integer :: lenmessage,rtnpos
 character(len=500) :: messtmp

!******************************************************************
!BEGIN EXECUTABLE SECTION

 if(message/=' ') then
  messtmp=message
  lenmessage=len(message)
! Here, split the message, according to the char(10)
! characters (carriage return). This technique is
! portable accross different OS.
  rtnpos=index(messtmp,ch10)
  do while(rtnpos/=0)
   write(unit, '(a)' ) trim(messtmp(1:rtnpos-1))
   messtmp=messtmp(rtnpos+1:lenmessage)
   lenmessage=lenmessage-rtnpos
   rtnpos=index(messtmp,ch10)
  end do
  write(unit, '(a)' ) trim(messtmp)
 else
  write(unit,*)
 end if

 if( index(trim(message),'BUG') /= 0 )then
  write(unit, '(a)' ) '  Action : contact ABINIT group.'
  write(unit,*)
 end if

 if( index(trim(message),'BUG') /= 0   .or. &
& index(trim(message),'Calculation completed') /= 0 )then
  if(nwarning<10000 .and. ncomment<1000)then
   write(unit, '(a,i5,a,i4,a)' ) &
&   '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
  else
   write(unit, '(a,i6,a,i6,a)' ) &
&   '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
  end if
  if(iexit/=0)then
   write(unit, '(a)' ) ' Note : exit requested by the user.'
  end if
 end if

 if( index(trim(message),'Exit') /= 0 )then
  iexit=1
 end if

!Count the number of warnings and comments. Only take into
!account unit 6, in order not to duplicate these numbers.
 if( index(trim(message),'WARNING') /= 0 .and. unit==6 )then
  nwarning=nwarning+1
 end if
 if( index(trim(message),'COMMENT') /= 0 .and. unit==6 )then
  ncomment=ncomment+1
 end if

end subroutine wrtout_myproc
!!***
