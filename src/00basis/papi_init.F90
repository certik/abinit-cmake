!{\src2tex{textfont=tt}}
!!****f* ABINIT/papi_init
!! NAME
!! papi_init
!!
!! FUNCTION
!!
!! This function initialize  papi hight level inteface , set up counters
!! to monitor PAPI_FP_OPS and PAPI_TOT_CYC events and start the counters
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine papi_init()


implicit none

!Arguments ------------------------------------

!Local variables-------------------------------

#ifdef HAVE_PAPI
#include "f90papi.h"

#define  C_FLOAT REAL
#define C_LONG_LONG INTEGER*8
#define C_INT INTEGER
character*(PAPI_MAX_STR_LEN) papi_errstr

C_INT :: retval
C_FLOAT :: unused1, unused2, unused4 
C_LONG_LONG :: unused3 

#endif

! *************************************************************************

#ifdef HAVE_PAPI

 retval = PAPI_VER_CURRENT
 call PAPIf_library_init(retval)
 if ( retval.NE.PAPI_VER_CURRENT) then
  PRINT *,'Problem library PAPI'
  endif

! First pass. Initializing counter


  call PAPIf_flops(unused1, unused2, unused3, unused4, retval)
  if (retval.NE.PAPI_OK) then
   PRINT *, 'Problem to initialize papi high level inteface'
   call papif_perror(retval,papi_errstr,retval)
   PRINT *, 'Error code', papi_errstr
  end if ! DEBUG


! call PAPIf_query_event(PAPI_FP_INS, retval) 
! if (retval .NE. PAPI_OK) then
! PRINT *, 'Problem query event'
! endif

! inializing papi hight level inteface , set up counters
! to monitor PAPI_FP_OPS and PAPI_TOT_CYC events and start the counters
! Subsequent calls will read the counters and return total
! real time, total process time, total floting point instructions
! or operations since the start of the mesurement and the Mflop/s rate
! since latests call to PAPI_flops
! call PAPIf_flops(real_time, proc_time, flpops, mflops, retval)
! if (retval.NE.PAPI_OK) then
! PRINT *, 'Problem to initialize papi high level inteface'
! endif

#endif

end subroutine papi_init
!!***
