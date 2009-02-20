!!****m* ABINIT/interfaces_00basis
!! NAME
!! interfaces_00basis
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/00basis
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

module interfaces_00basis

 implicit none

interface
 subroutine initpapichoice(string, lenstr, papichoice)
  implicit none
  integer,intent(in) :: lenstr
  integer,intent(out) :: papichoice
  character(len=*),intent(in) :: string
 end subroutine initpapichoice
end interface

interface
 subroutine leave_myproc
  implicit none
 end subroutine leave_myproc
end interface

interface
 subroutine papi_init()
  implicit none
 end subroutine papi_init
end interface

interface
 subroutine timab(nn,option,tottim)
  use defs_basis
  implicit none
  integer,intent(in) :: nn
  integer,intent(in) :: option
  real(dp),intent(out) :: tottim(2)
 end subroutine timab
end interface

interface
 subroutine time_accu(nn,return_ncount,tottim, totflops, totftimes)
  use defs_basis
  implicit none
  integer,intent(in) :: nn
  integer,intent(out) :: return_ncount
  real(dp),intent(out) :: totflops
  real(dp),intent(out) :: totftimes(2)
  real(dp),intent(out) :: tottim(2)
 end subroutine time_accu
end interface

interface
 subroutine timein(cpu,wall)
  use defs_basis
  implicit none
  real(dp),intent(out) :: cpu
  real(dp),intent(out) :: wall
 end subroutine timein
end interface

interface
 subroutine wrtout_myproc(unit,message)
  implicit none
  integer,intent(in) :: unit
  character(len=500),intent(inout) :: message
 end subroutine wrtout_myproc
end interface

end module interfaces_00basis
!!***
