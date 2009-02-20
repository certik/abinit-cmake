!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_io_tools
!! NAME
!!  m_io_tools
!!
!! FUNCTION
!!  This module contains basic tools to deal with Fortran IO
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
module m_io_tools

 use defs_basis
 use m_fstrings
 use m_errors, only : assert

 implicit none

! === List of available public routines and functions ===
 public ::         &
&  get_unit,       &  ! Get a free unit if no argument is specified or report the unit associated to a file name
&  open_file,      &  ! Open a file after having performed some bacis tests
&  is_open,        &  ! .TRUE. if file is open
&  is_connected,   &  ! .TRUE. if file is connected to a logical unit number
&  prompt,         &  ! Simple prompt
&  read_line,      &  ! Read line from unit ignoring blank lines and deleting comments beginning with !
&  flush_unit         ! Wrapper to the intrinsic flush routine, not implemented by every compiler 
!@END ABIDOC

 interface get_unit
  module procedure get_free_unit
  module procedure get_unit_from_fname
 end interface

 interface is_open
  module procedure is_open_unit
  module procedure is_open_fname
 end interface

 interface prompt 
  module procedure prompt_int_0D
  module procedure prompt_rdp_0D
  module procedure prompt_string
 end interface

 !PRIVATE
  integer,parameter :: STDIN=std_in
  integer,parameter :: STDOUT=std_out
  integer,parameter :: MIN_UNIT_NUMBER=10  ! Fortran does not define the range for logical unit numbers (they not be negative)
  integer,parameter :: MAX_UNIT_NUMBER=99  ! The following values should be safe
  integer,parameter :: IO_MAX_LEN=500
  character(len=1),parameter :: BLANK=' ' 

  ! === For Interactive sessions ===
  integer,parameter :: IO_EOT=-1           ! End of transmission i.e CTRL+D 
  character(len=1),parameter :: COMMENT='!'
  character(len=4),parameter :: PS1='>>> '
  character(len=4),parameter :: PS2='??? '

  ! === Built in IO exceptions, negative identifiers are used ===
  integer,parameter :: ERROR_UNKNOWN=-2  ! No units are available for Fortran I/O
  integer,parameter :: IO_NO_AVAILABLE_UNITS=-3  ! No units are available for Fortran I/O
  integer,parameter :: IO_FILE_EXISTS=-4         ! File already exists
  integer,parameter :: IO_FILE_DOES_NOT_EXIST=-5 ! File does not already exist
  integer,parameter :: IO_FILE_IS_OPEN=-6        ! File is already open 
  integer,parameter :: IO_FILE_NOT_ASSOCIATED=-7 ! File is not associated with any unit 
  !integer,parameter :: IO_END_OF_FILE=-8        ! End of file reached 

CONTAINS  !===========================================================

!!***

!!****f* m_io_tools/get_unit
!! NAME
!!  get_unit
!!
!! FUNCTION
!!  Obtain a logical Fortran unit. 
!!  A free unit is reported if no argument is specified.
!!  If the file name is supplied, the function reports the unit number 
!!  associated to the file
!!
!! INPUTS
!!
!! OUTPUT
!!  The unit number (free unit or unit associated to the file)
!!  Raises:
!!   IO_NO_AVAILABLE_UNITS if no logical unit is free (!)
!!   IO_FILE_NOT_ASSOCIATED if the file is not linked to a logical unit
!!
!! SOURCE

integer function get_free_unit() 

!Local variables-------------------------------

  integer :: iunt
  logical :: is_open
! *********************************************************************

 do iunt=MIN_UNIT_NUMBER,MAX_UNIT_NUMBER
  inquire(unit=iunt,opened=is_open)
  if (.not.is_open) then
   get_free_unit=iunt ; RETURN
  end if
 end do
 get_free_unit=IO_NO_AVAILABLE_UNITS

end function get_free_unit
!!***

integer function get_unit_from_fname(fname)

!Arguments ------------------------------------

 character(len=*),intent(in) :: fname

!Local variables-------------------------------
 integer :: unit
! *********************************************************************

 inquire(file=fname,number=unit)

 get_unit_from_fname=unit
 if (unit==-1) get_unit_from_fname=IO_FILE_NOT_ASSOCIATED

end function get_unit_from_fname
!!***

!!****f* m_io_tools/open_file
!! NAME
!!  open_file
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!  
!!
!! SOURCE

integer function open_file(fname,status,form,access,unit) result(unt)

!!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,optional,intent(in) :: unit
 character(len=*),intent(in) :: fname,status,form,access

!!Local variables-------------------------------
 integer :: ios
 logical :: exists,is_open,ltest
! *********************************************************************

 !ltest=(ANY(toupper(status)==(/'NEW','OLD','UNKNOWN'/)))
 !call assert(ltest,'Unknown status',__FILE__,__LINE__)
 !ltest=(ANY(toupper(form)==(/'FORMATTED','UNFORMATTED'/)))
 !call assert(ltest,'Unknown form',__FILE__,__LINE__)
 !ltest=(ANY(toupper(access)==(/'SEQUENTIAL','DIRECT'/)))
 !call assert(ltest,'Unknown access',__FILE__,__LINE__)

 inquire(file=TRIM(fname),exist=exists,opened=is_open)
 if (is_open) then
  unt=IO_FILE_IS_OPEN ; RETURN
 end if

 if (exists.and.toupper(status)=='NEW') then
  unt=IO_FILE_EXISTS ; RETURN
 end if

 if (.not.exists.and.toupper(status)=='OLD') then
  unt=IO_FILE_DOES_NOT_EXIST ; RETURN
 end if

 if (PRESENT(unit)) then
  unt=unit
 else
  unt=get_unit() ; if (unt==IO_NO_AVAILABLE_UNITS) RETURN
 end if
 !
 ! === Now we can open the file ===
 open(unit=unt,file=fname,status=status,form=form,iostat=ios)
 if (ios/=0) unt=ERROR_UNKNOWN

end function open_file
!!***

!!****f* m_io_tools/is_connected
!! NAME
!!  is_connected
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!  
!!
!! SOURCE

logical function is_connected(unit,fname)

!!Arguments ------------------------------------

 integer,intent(in) :: unit
 character(len=*),intent(in) :: fname

!Local variables-------------------------------
 integer :: unt_found
 logical :: is_open
! *********************************************************************

 inquire(file=fname,number=unt_found,opened=is_open)
 is_connected=(is_open.and.(unt_found==unit))

end function is_connected
!!***

!!****f* m_io_tools/is_open
!! NAME
!!  is_open
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!  
!!
!! SOURCE

logical function is_open_unit(unit)

!Arguments ------------------------------------

 integer, intent(in) :: unit
! *********************************************************************

 inquire(unit=unit,opened=is_open_unit)

end function is_open_unit
!!***

logical function is_open_fname(fname)

!Arguments ------------------------------------

 character(len=*),intent(in) :: fname
! *********************************************************************

 inquire(file=fname,opened=is_open_fname)

end function is_open_fname
!!***

!!****f* m_io_tools/check_unit
!! NAME
!!  check_unit
!!
!! FUNCTION
!!  Test the result of oper file and, if something went wrong, print an error 
!!  message with useful information of the nature of the problem.
!!
!! INPUTS
!!
!! OUTPUT
!!  
!!
!! PARENTS
!!
!! CHILDREN
!!      flush
!!
!! SOURCE

subroutine check_unit(unit,fatal)

!Arguments ------------------------------------

 integer,intent(in) :: unit
 logical,optional,intent(in) :: fatal

!Local variables-------------------------------
 integer :: ierr=-1
 logical :: is_fatal
! *********************************************************************

 is_fatal=.TRUE. ; if (PRESENT(fatal)) is_fatal=fatal

 select case (unit)
 case (ERROR_UNKNOWN) 
  write(*,*)' Error unknown'
 case (IO_NO_AVAILABLE_UNITS) 
  write(*,*)' Not able to find a free Fortran unit! '
 case (IO_FILE_EXISTS) 
  write(*,*)' File already exists'
 case (IO_FILE_DOES_NOT_EXIST) 
  write(*,*)'File does not exist'
 case (IO_FILE_IS_OPEN)
  write(*,*)'File is already open'
 case (IO_FILE_NOT_ASSOCIATED) 
  write(*,*)'File is not associated with any unit'
 case default
  ierr=0
  continue
 end select

 if (ierr/=0.and.is_fatal) STOP

end subroutine check_unit
!!***


!!****f* m_io_tools/prompt
!! NAME
!!  prompt
!!
!! FUNCTION
!!  Simple prompt
!!
!! INPUTS
!!
!! OUTPUT
!!  
!!
!! SOURCE

subroutine prompt_int_0D(msg,ivalue)

!Arguments ------------------------------------

 character(len=*),intent(in) :: msg
 integer,intent(out) :: ivalue
 character(len=4) :: PS
!Local variables-------------------------------
 integer :: ios
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
  write(STDOUT,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
  read(STDIN,*,IOSTAT=ios)ivalue
  if (ios==IO_EOT) call prompt_exit()
  PS=PS2
 end do

end subroutine prompt_int_0D
!!***

subroutine prompt_rdp_0D(msg,rvalue)

!Arguments ------------------------------------

 character(len=*),intent(in) :: msg
 real(dp),intent(out) :: rvalue
 character(len=4) :: PS
!Local variables-------------------------------
 integer :: ios
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
  write(STDOUT,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
  read(STDIN,*,IOSTAT=ios)rvalue
  if (ios==IO_EOT) call prompt_exit()
  PS=PS2
 end do

end subroutine prompt_rdp_0D
!!***


subroutine prompt_string(msg,string)

!Arguments ------------------------------------

 character(len=*),intent(in) :: msg
 character(len=*),intent(out) :: string
!Local variables-------------------------------
 integer :: ios
 character(len=4) :: PS
! *********************************************************************

 ios=-1 ; PS=PS1
 do while (ios/=0)
  write(STDOUT,'(a)',ADVANCE='NO')PS//TRIM(msg)//BLANK
  read(STDIN,*,IOSTAT=ios)string
  if (ios==IO_EOT) call prompt_exit()
  PS=PS2
 end do

end subroutine prompt_string
!!***

subroutine prompt_exit()

!Arguments ------------------------------------

 integer :: ios
 character(len=IO_MAX_LEN) :: ans
! *********************************************************************

 write(STDOUT,*)
 ios=-1
 do while (ios/=0.or.(ans/='y'.or.ans/='n'))
  write(STDOUT,'(a)')' Do you really want to exit (y/n)? '
  read(STDIN,*,IOSTAT=ios)ans
  if (ans=='y') STOP
  if (ans=='n') RETURN
 end do

end subroutine prompt_exit
!!***

!!***
!! NAME
!!  read_line
!!
!! FUNCTION
!!  Reads line from unit=std_in_ or unit if specified, ignoring blank lines
!!  and deleting comments beginning with !

subroutine read_line(line,ios,unit)

!Arguments ------------------------------------

 character(len=*),intent(out):: line
 integer,optional,intent(in) :: unit
 integer,intent(out) :: ios

!Local variables-------------------------------
 integer :: ipos,unt
! *********************************************************************

 unt=STDIN ; if (PRESENT(unit)) unt=unit

 do  
  read(unt,'(a)',iostat=ios) line  ! read input line
  if (ios/=0) RETURN
  line=ADJUSTL(line) ; ipos=INDEX(line,COMMENT)
  if (ipos==1) CYCLE
  if (ipos/=0) line=line(:ipos-1)
  if (LEN_TRIM(line)/=0) EXIT
 end do

end subroutine read_line
!!***

!!****f* m_io_tools/flush_unit
!! NAME
!! flush_unit
!!
!! FUNCTION
!! Wrapper to the standard flush_unit routine 
!! There might be problems with XLF and ABSOFT with the underscore
!!
!! INPUTS
!!  unit(optional)=the unit number to be flushed (if not specified ALL open units are flushed
!!
!! NOTES
!!  Only available if the compiler implements this intrinsic procedure
!!
!! PARENTS
!!
!! CHILDREN
!!      flush
!!
!! SOURCE

subroutine flush_unit(unit)
    
!Arguments ------------------------------------

 integer,optional,intent(in) :: unit

!Local variables-------------------------------
 integer :: unt
 logical :: is_open
!************************************************************************

#if defined HAVE_FLUSH
 if (PRESENT(unit)) then 
 unt=unit
 inquire(unit=unt,opened=is_open)
 if (is_open) call flush(unt)
 else 
  call flush()
 end if
#endif

end subroutine flush_unit
!!***

end module m_io_tools
!!***
