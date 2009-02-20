!!****m* ABINIT/interfaces_12parser
!! NAME
!! interfaces_12parser
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/12parser
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_12parser

 implicit none

interface
 subroutine inarray(b1,cs,dprarr,intarr,marr,narr,string,typevarphys)
  use defs_basis
  implicit none
  integer,intent(inout) :: b1
  integer,intent(in) :: marr
  integer,intent(in) :: narr
  character(len=fnlen),intent(inout) :: cs
  character(len=*),intent(in) :: string
  character(len=3),intent(in) :: typevarphys
  real(dp),intent(out) :: dprarr(marr)
  integer,intent(out) :: intarr(marr)
 end subroutine inarray
end interface

interface
 subroutine incomprs(string,length)
  implicit none
  integer,intent(out) :: length
  character(len=*),intent(inout) :: string
 end subroutine incomprs
end interface

interface
 subroutine inread(string,ndig,typevarphys,outi,outr,errcod)
  use defs_basis
  implicit none
  integer,intent(out) :: errcod
  integer,intent(in) :: ndig
  integer,intent(out) :: outi
  real(dp),intent(out) :: outr
  character(len=*),intent(in) :: string
  character(len=3),intent(in) :: typevarphys
 end subroutine inread
end interface

interface
 subroutine inreplsp(string)
  implicit none
  character(len=*),intent(inout) :: string
 end subroutine inreplsp
end interface

interface
 subroutine instrng (filnam,lenstr,option,strln,string)
  implicit none
  integer,intent(out) :: lenstr
  integer,intent(in) :: option
  integer,intent(in) :: strln
  character(len=*),intent(in) :: filnam
  character(len=*),intent(out) :: string
 end subroutine instrng
end interface

interface
 subroutine intagm(dprarr,intarr,jdtset,marr,narr,string,token,tread,typevarphys)
  use defs_basis
  implicit none
  integer,intent(in) :: jdtset
  integer,intent(in) :: marr
  integer,intent(in) :: narr
  integer,intent(out) :: tread
  character(len=*),intent(in) :: string
  character(len=*),intent(inout) :: token
  character(len=3),intent(in) :: typevarphys
  real(dp),intent(out) :: dprarr(marr)
  integer,intent(out) :: intarr(marr)
 end subroutine intagm
end interface

end module interfaces_12parser
!!***
