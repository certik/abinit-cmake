!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_time
!! NAME
!! defs_time
!!
!! FUNCTION
!! This module contains accumulators for the timer.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!! Include the name of all routines: better modularity
!!
!! PARENTS
!!    timab,time_accu
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_time

 use defs_basis

 implicit none

!mtim determines the maximum number of "timing slots" available
 integer,parameter :: mtim=599

!timeopt is a flag which indicates the suppression or not of the timing.
 integer :: timopt=1

!papiopt is a flag which indicates if there is or not an analysis of speed execution  
!imade by by defaut the analysis is not done 
 integer :: papiopt=0

!!initpapiopt is a flag which permit initialisation of overall timingi by papi 
!! and speed of run
 integer :: initpapiopt=1 

!Number of times that the routine has been called
 integer :: ncount(mtim)

!Accumulating cpu time (1) and wall to wall time (2) for each "timing slots"
 real(dp)  :: acctim(2,mtim),tzero(2,mtim)

!Accumulating number of floating point operation and cpu time (1) and wall to wall times (2)
!for each "performance slot"
 real(dp) :: papi_accflops(mtim), papi_acctim(2,mtim)

!reference value for number of floating point operation and time ( cpu and wall)! for each performance slot
 real(dp) :: flops(mtim) , papi_tzero(2,mtim)


!Ellapsed time and ellased number of floating point operation since a reference
real(dp) :: papi_tottim(2,mtim), papi_totflops(mtim)

end module defs_time
!!***
