!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_errors
!! NAME
!!  m_errors
!!
!! FUNCTION
!!  This module contains low-level procedure to 
!!  check assertions and to handle errors.
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!@ ABIDOC
MODULE m_errors

 use defs_basis
 !£ use defs_infos

 use interfaces_01manage_mpi     ! abilint fails in recognizing wrtout and leave_new

 implicit none

 private

! === List of available public routines and functions ===
 public ::       &
&  assert_eq,    &   ! Report and die gracefully if integers not all equal (used for size checking).
&  assert,       &   ! Report and die if any logical is false (used for arg range checking).
&  msg_hndl,     &   ! Basic Error handlers 
&  io_hndl,      &   ! Error handler for IO operations on external files
&  wrap_wrtout        ! Simple wrapper around wrtout to pass assuems shape strings
!@END ABIDOC

 interface assert_eq  
  module procedure assert_eq2,assert_eq3,assert_eq4
  module procedure assert_eqn
 end interface

 interface assert 
  module procedure assert1,assert2,assert3,assert4
  module procedure assert_v
 end interface

CONTAINS  !===========================================================

!!***

!!****f* m_errors/assert_eq
!! NAME
!!  assert_eq
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
!!
!! INPUT 
!!  l1,l2,.. Integers to be checked (array version is also provided)
!!  message(len=*)=tag with additiona information
!!
!! SOURCE

function assert_eq2(l1,l2,message,f90name,line)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,intent(in) :: l1,l2,line
 integer :: assert_eq2
 character(len=*),intent(in) :: message,f90name
! *************************************************************************

 if (l1==l2) then
  assert_eq2=l1
 else
  call msg_hndl(message,f90name,line,'BUG')
 end if

end function assert_eq2

function assert_eq3(l1,l2,l3,message,f90name,line)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,intent(in) :: l1,l2,l3,line
 integer :: assert_eq3
 character(len=*),intent(in) :: message,f90name
! *************************************************************************

 if (l1==l2.and.l2==l3) then
  assert_eq3=l1
 else
  call msg_hndl(message,f90name,line,'BUG')
 end if

end function assert_eq3

function assert_eq4(l1,l2,l3,l4,message,f90name,line)

!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,intent(in) :: l1,l2,l3,l4,line
 integer :: assert_eq4
 character(len=*),intent(in) :: message,f90name
! *************************************************************************

 if (l1==l2.and.l2==l3.and.l3==l4) then
  assert_eq4=l1
 else
  call msg_hndl(message,f90name,line,'BUG')
 end if

end function assert_eq4

function assert_eqn(nn,message,f90name,line)

!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,intent(in) :: line
 integer :: assert_eqn
 character(len=*),intent(in) :: message,f90name
!arrays
 integer,intent(in) :: nn(:)
! *************************************************************************

 if (ALL(nn(2:)==nn(1))) then
  assert_eqn=nn(1)
 else
  call msg_hndl(message,f90name,line,'BUG')
 end if

end function assert_eqn
!!***

!!****f* m_errors/assert
!! NAME
!!  assert
!!
!! FUNCTION
!!  Routines for argument checking and error handling.
!!  Report and die if any logical is false (used for arg range checking).
!!
!! INPUT 
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additiona information
!!
!! SOURCE

subroutine assert1(l1,message,file,line)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Unknown'
! *************************************************************************

 if (.not.l1) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name=TRIM(file)
  call msg_hndl(message,f90name,f90line,'BUG')
 end if

end subroutine assert1

subroutine assert2(l1,l2,message,file,line)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Unknown'
! *************************************************************************

 if (.not.(l1.and.l2)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name=TRIM(file)
  call msg_hndl(message,f90name,f90line,'BUG')
 end if

end subroutine assert2

subroutine assert3(l1,l2,l3,message,file,line)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2,l3

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Unknown'
! *************************************************************************

 if (.not.(l1.and.l2.and.l3)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name=TRIM(file)
  call msg_hndl(message,f90name,f90line,'BUG')
 end if

end subroutine assert3

subroutine assert4(l1,l2,l3,l4,message,file,line)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2,l3,l4

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Unknown'
! *************************************************************************

 if (.not.(l1.and.l2.and.l3.and.l4)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name=TRIM(file)
  call msg_hndl(message,f90name,f90line,'BUG')
 end if

end subroutine assert4

subroutine assert_v(n,message,file,line)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: n(:)

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Unknown'
! *************************************************************************

 if (.not.ALL(n)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name=TRIM(file)
  call msg_hndl(message,f90name,f90line,'BUG')
 end if

end subroutine assert_v
!!***

!!****f* m_errors/msg_hndl
!! NAME
!!  msg_hndl
!!
!! FUNCTION
!!  Basic error handler for abinit. This routine is usually  interfaced though some macro 
!!  defined in FIXME
!!
!! INPUT 
!!  message=string containing additional information on the nature of the problem
!!  level=string defining the type of problem. Possible values are
!!   COMMENT
!!   WARNING
!!   ERROR
!!   BUG
!!  line=line number of the file where problem occurred
!!  f90name=name of the f90 file containing the caller
!!
!! PARENTS
!!      m_errors
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine msg_hndl(message,f90name,line,level)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 integer,intent(in) :: line
 character(len=*),intent(in) :: level,message,f90name

!Local variables-------------------------------
 character(len=500) :: msg
 character(len=10) :: lnum
! *********************************************************************

 !FIXME problem with dependencies !!
 !£ call int2char(line,lnum)
 msg=TRIM(f90name)//":"//TRIM(lnum)

 select case (level)

 !FIXME format might be improved, problem if message is more than one line
 case ('COMMENT')
  write(msg,'(a,2x,3a,2x,a)')ch10,&
   TRIM(msg),' COMMENT- ',ch10,TRIM(message)
  call wrtout(std_out,msg,'COLL') 

 case ('WARNING')
  write(msg,'(a,2x,3a,2x,a)')ch10,&
   TRIM(msg),' WARNING- ',ch10,TRIM(message)
  call wrtout(std_out,msg,'COLL') 

 case ('ERROR')
  write(msg,'(a,2x,3a,2x,a)')ch10,&
   TRIM(msg),' ERROR- ',ch10,TRIM(message)
  !£ call print_defs_infos()
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')

 case ('BUG')
  write(msg,'(a,2x,3a,2x,a)')ch10,&
   TRIM(msg),' BUG- ',ch10,TRIM(message)
  !£ call print_defs_infos() 
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')

 case default 
  write(msg,'(4a)')ch10,&
&  ' msg_hndl: BUG - ',ch10,&
&  ' Wrong value for level '
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end select

end subroutine  msg_hndl
!!***

!!****f* m_errors/io_hndl
!! NAME
!!  io_hndl
!!
!! FUNCTION
!!  Basic error handler for I/O operations on external files. 
!!  This routine is usually interfaced though some macro 
!!
!! INPUT 
!!  unit=Fortran unit number
!!  ios=IO status
!!  f90name=name of the subroutine where error occurs
!!  line=line number in the f90name file
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine io_hndl(ios,unit,f90name,line)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 integer,intent(in) :: unit,ios,line
 character(len=*),intent(in) :: f90name

!Local variables-------------------------------
 character(len=fnlen) :: fname
 character(len=10) :: lnum,s_ios,s_unt
 character(len=500) :: msg
 logical :: lexist,lopened,lnamed
! *********************************************************************

 !FIXME
 !£ call int2char(line,lnum) ; msg=TRIM(f90name)//':'//TRIM(lnum)
 !£ call int2char(ios,s_ios) 
 !£ call int2char(unit,s_unt)

 write(msg,'(8a)')ch10,&
& TRIM(msg),' I/O ERROR- ',ch10,&
& ' while operating on unit ',TRIM(s_unt),', iostat = ',TRIM(s_ios)
 call wrtout(std_out,msg,'COLL') 

 inquire(unit=unit,exist=lexist,named=lnamed,opened=lopened)
 fname='None' ; if (lnamed) inquire(unit=unit,name=fname)

 write(msg,'(2a,2(a,l7,a),2a)')&
& ' Inquire reports : ',ch10,&
& '  exist  = ',lexist,ch10,&
& '  opened = ',lopened,ch10,&
& '  name   = ',TRIM(fname)
 call wrtout(std_out,msg,'COLL') 
 
 call leave_new('COLL')

end subroutine io_hndl
!!***

!Simple wrapper around wrout to pass assumed shape strings
subroutine wrap_wrtout(unit,message,mode_paral)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 integer,intent(in) :: unit
 character(len=4),intent(in) :: mode_paral
 character(len=*),intent(in) :: message

!Local variables-------------------------------
 character(len=500) :: msg
! *********************************************************************

 msg=message
 call wrtout(unit,msg,mode_paral)

end subroutine wrap_wrtout
!!***

END MODULE m_errors
!!***
