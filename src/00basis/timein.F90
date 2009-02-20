!{\src2tex{textfont=tt}}
!!****f* ABINIT/timein
!! NAME
!! timein
!!
!! FUNCTION
!! Timing routine.
!! Returns cpu and wall clock time in seconds since some arbitrary start.
!!
!! For wall clock time, call the F90 intrinsic date_and_time .
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, LSI, MM, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (no inputs)
!!
!! OUTPUT
!!  cpu= cpu time in seconds
!!  wall= wall clock time in seconds
!!
!! NOTES
!! For CPU time, contains machine-dependent code (choice will be selected
!! by c preprocessor).
!! Note that all supported machines are listed explicitly below; there
!! is no "else" which covers "other".  The C preprocessor will place
!! a spurious line of code (see below) into the fortran source unless
!! preprocessed with -Dflag where flag refers to one of the supported machines.
!!
!! WARNING: the following list is no more accurate (YP 20060530)
!!
!! Presently supported flags: "ibm", "hp", "P6", "dec_alpha", "sgi",
!!    "T3E", "vpp", "sun", "mac", "nec", "sr8k" , "VMS".
!! Previously supported flags:  "ultrix". Might still work !
!!
!! Calls machine-dependent "mclock" for "ibm" .
!! Calls machine-dependent "second" for "T3E"
!! Calls ANSI C subroutine "cclock" for "hp" and "sgi".
!! Calls machine-dependent "etime" for "P6", "mac", "dec_alpha", "sun", "nec" .
!! Calls machine-dependent "clock" for "vpp"
!! Calls machine-dependent "xclock" for "sr8k"
!! Calls F95 intrinsic "cpu_time" for "VMS"
!!
!! PARENTS
!!      abinit,aim,aim_follow,anaddb,chkexi,cpdrv,drvaim,mkifc9,pclock,rdddb9
!!      rsiaf9,rsurf,surf,thm9,timab
!!
!! CHILDREN
!!      cclock,clock,cpu_time,date_and_time,leave_new,system_clock,wrtout
!!      xclock
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine timein(cpu,wall)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: cpu,wall

!Local variables-------------------------------
!no_abirules
 integer, parameter :: nday(24)=(/31,28,31,30,31,30,31,31,30,31,30,31,&
&                                 31,28,31,30,31,30,31,31,30,31,30,31/)
 integer, save :: month_init,month_now,start=1,year_init
 integer :: count,count_max,count_rate,months
 integer :: values(8)
!Machine-dependent declarations
#if defined USE_CCLOCK
 real(dp) :: tmp
#elif defined FC_IBM
 integer :: mclock
#elif defined T3E
 real(dp) :: second
#elif defined FC_SUN || defined FC_NEC
 real :: tmp(2)
 real :: etime
#elif defined i386 || defined FC_COMPAQ || defined OS_MACOSX
 real :: tmp(2)           !real array only needed by etime
 real(dp) :: etime
#elif defined FC_NAG
 real :: second
#endif
 character(len=10) :: date,time,zone
 character(len=500) :: message

! *************************************************************************

!CPU time _______________________________________

!It is possible to suppress the call to an external routine, and leave
!the cpu time to 0.0d0, provided the timing of the timer is suppressed in timana.f
!(simply set the loop counter maximum value to 1 in that routine)
 cpu = 0.0d0

!Machine-dependent timers
#if defined USE_CCLOCK

 call cclock(tmp)
 cpu = tmp

#elif defined FC_IBM

 cpu = mclock()*0.01d0

#elif defined T3E

 cpu = second()

#elif defined FC_HP

 call cclock(cpu)

#elif defined FC_MIPSPRO

!call cclock(cpu)  ! XG041104 This timing routine caused some problem on the SGI "Spinoza" machine
!This is the Fortran standard subroutine.
 call system_clock(count,count_rate,count_max)
 cpu=dble(count)/dble(count_rate)

#elif defined i386 || defined OS_MACOSX

#if defined FC_NAG
!Here, the f95 intrinsic is used
 call cpu_time(second)
 cpu = second
#else
 cpu = etime(tmp)
#endif

#elif defined FC_COMPAQ || defined FC_SUN || defined FC_NEC

 cpu = etime(tmp)

#elif defined FC_FUJITSU

 call clock(cpu,0,2)

#elif defined FC_HITACHI

 call xclock(cpu,5)

#elif defined VMS

 call cpu_time( cpu )

#else

!This is the Fortran standard subroutine,
!might not always be sufficiently accurate
!Prior to v5, each machine type had its own timer.
 call system_clock(count,count_rate,count_max)
 cpu=dble(count)/dble(count_rate)

#endif

!Wallclock time ________________________________

!write(std_out,*)' timein : before if start '
!wall = 0.0d0

!The following section of code is standard F90, but
!it is useful only if the intrinsics
!date_and_time is accurate at the 0.01 sec level,
!which is not the case for a P6 with the pghpf compiler ...
!Year and month initialisation
 if(start==1)then
  start=0
  call date_and_time(date,time,zone,values)
  year_init=values(1)
  month_init=values(2)
 end if

!write(std_out,*)' timein : before date_and_time '

!Uses intrinsic F90 subroutine Date_and_time for
!wall clock (not correct when a change of year happen)
 call date_and_time(date,time,zone,values)

!Compute first the number of seconds from the beginning of the month
 wall=(values(3)*24.0d0+values(5))*3600.0d0&
& +values(6)*60.0d0+values(7)+values(8)*0.001d0

!If the month has changed, compute the number of seconds
!to be added. This fails if the program ran one year !!
 month_now=values(2)
 if(month_now/=month_init)then
  if(year_init+1==values(1))then
   month_now=month_now+12
  end if
  if(month_now<=month_init)then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' timein : BUG -',ch10,&
&   '  Problem with month and year numbers.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
  end if
  do months=month_init,month_now-1
   wall=wall+86400.0d0*nday(months)
  end do
 end if

!Now take into account bissextile years (I think 2000
!is bissextile, but I am not sure ...)
 if(mod(year_init,4)==0 .and. month_init<=2 .and. month_now>2)&
& wall=wall+3600.0d0
 if(mod(values(1),4)==0 .and. month_init<=14 .and. month_now>14)&
& wall=wall+3600.0d0

!write(std_out,*)' timein : wall at exit ',wall

end subroutine timein
!!***
