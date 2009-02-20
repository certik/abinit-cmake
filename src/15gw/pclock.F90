!{\src2tex{textfont=tt}}
!!****f* ABINIT/pclock
!! NAME
!! pclock
!!
!! FUNCTION
!! Print the timing at point number itimpt
!! if itimpt=0 then start the clock
!!             else print the cpu time elapsed since clock started
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  itimpt=see function description
!!
!! OUTPUT
!!  (only printing)
!!
!! PARENTS
!!      cchi0q0,rdkss,rdscr,screening,sigma,testlda,testscr
!!
!! CHILDREN
!!      timein
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pclock(itimpt)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itimpt

!Local variables-------------------------------
!scalars
 real(dp),save :: cstart,wstart
 real(dp) :: cpu,wall
 character(len=500) :: msg

! *************************************************************************

 if (itimpt==0) then
  ! === Start clock ===
  call timein(cstart,wstart)
  cpu=zero ; wall=zero
 else
  call timein(cpu,wall)
  cpu=cpu-cstart
  wall=wall-wstart
  if (itimpt==9999) then
   write(msg,'(a,f13.1,a,f13.1)')&
&   '+Overall time at end (sec) : cpu=',cpu,'  wall=',wall
   call wrtout(std_out,msg,'PERS')
  else
   write(msg,'(2a,i4,a,a,f10.2,3a,f10.2,2a)')ch10,&
&   ' timing point number ',itimpt,ch10,&
&   '-cpu time  = ',cpu,' seconds',ch10,&
&   '-real time = ',wall,' seconds',ch10
   call wrtout(std_out,msg,'PERS')
  end if
 end if

end subroutine pclock
!!***
