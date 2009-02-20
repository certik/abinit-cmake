!{\src2tex{textfont=tt}}
!!****f* ABINIT/isfile
!!
!! NAME
!! isfile
!!
!! FUNCTION
!! Inquire Status of FILE
!! Checks that for status =
!! 'old': file already exists
!! 'new': file does not exist; if file exists,
!! filnam is modified to filnam.A or filnam.B,....
!!    (except on OpenVMS where the file version numbers are used)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! filnam=character string to specify filename
!! status='old' or 'new'
!!
!! OUTPUT
!! stops processing if old file does not exist; changes name
!! and returns new name in redefined filnam if new file
!! already exists
!!
!! PARENTS
!!      anaddb,crho,iofn1,mkphdos,newsp,paw_symcprj,screening
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine isfile(filnam,status)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=3),intent(in) :: status
 character(len=fnlen),intent(inout) :: filnam

!Local variables-------------------------------
!scalars
 integer :: ii,ios
 logical :: ex
 character(len=500) :: message
 character(len=fnlen) :: trialnam
!arrays
 character(len=1) :: alpha(26)

! *************************************************************************

 alpha(1:26)=(/'A','B','C','D','E','F','G','H','I','J','K','L','M','N',&
& 'O','P','Q','R','S','T','U','V','W','X','Y','Z'/)

 if (status=='old') then

! Check that old file exists
  inquire(file=trim(filnam),iostat=ios,exist=ex)
  if (ios/=0) then
   write(message, '(a,a,a,a,a,a,a,i8,a,a)' ) ch10,&
&   ' isfile : ERROR -',ch10,&
&   '  Checks for existence of file  ',trim(filnam),ch10,&
&   '  but INQUIRE statement returns error code',ios,ch10,&
&   '  Action : identify which problem appears with this file.'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  else if (.not.ex) then
   write(message, '(a,a,a,a,a,a,a,a)' ) ch10,&
&   ' isfile : ERROR -',ch10,&
&   '  Checks for existence of file  ',trim(filnam),ch10,&
&   '  but INQUIRE finds file does not exist.',&
&   '  Action : check file name and re-run.'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  end if

 else if (status=='new') then

  trialnam=trim(filnam)

  do ii=1,27

!  Check that new output file does NOT exist
   inquire(file=trim(trialnam),iostat=ios,exist=ex)

#ifdef VMS
!  For VMS, use the file version numbers, so there is no need to rename
!  a file if it already exists. With the following, the loop
!  will be exited gracefully, and the name will be preserved.
   ex=.false.
#endif

   if (ios/=0) then

!   There is a problem => stop
    write(message, '(a,a,a,a,a,a,a,i8,a,a)' ) ch10,&
&    ' isfile : ERROR -',ch10,&
&    '  Checks for existence of file  ',trim(trialnam),ch10,&
&    '  but INQUIRE statement returns error code',ios,ch10,&
&    '  Action : identify which problem appears with this file.'
    call wrtout(6,message,'PERS')
    call leave_new('PERS')

   else if (ex) then

    write(message, '(a,a,a,a,a,a,a)' ) ch10,&
&    ' isfile : WARNING -',ch10,&
&    '  Finds that output file ',trim(trialnam),&
&    ch10,' already exists.'
    call wrtout(6,message,'PERS')
!   'New' file already exists; define a new file name
    if (ii<=26) then
     trialnam=trim(filnam)//alpha(ii)
     write(message, '(a,a,a)' ) ' new name assigned:',trim(trialnam),ch10
     call wrtout(6,message,'PERS')
     cycle
    else
     write(message, '(a,a,a,a,a,a)' ) ch10,&
&     ' isfile : ERROR -',ch10,&
&     '  Have used up all names of the form filename.[A-Z]',ch10,&
&     '  Action : clean up your directory and start over.'
     call wrtout(6,message,'PERS')
     call leave_new('PERS')
    end if

   else

!   The name (or the new name) is correct
    exit

   end if

!  End loop on ii : scan the alphabet.
!  There is a "cycle" and an "exit" in the loop.
  end do

  filnam=trim(trialnam)

 else

! status not recognized
  write(message, '(a,a,a,a,a,a)' ) ch10,&
&  ' isfile : BUG -',ch10,&
&  '  Input status= ',status,' not recognized.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')

 end if

end subroutine isfile
!!***
