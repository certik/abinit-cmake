!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_GWdefs
!! NAME
!! m_GWdefs
!!
!! FUNCTION
!! This module contains definitions for a number of named constants used in the GW part
!! of abinit
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module m_gwdefs

 use defs_basis

 implicit none

! Unit number for formatted files produced by GW calculations.
! These files are not supposed to be read by abinit therefore
! their names and unit numbers are not defined in the Dtfil% structure.

 integer,parameter :: unt_gw  = 21  ! GW corrections
 integer,parameter :: unt_sig = 22  ! Self-energy as a function of frequency
 integer,parameter :: unt_sgr = 23  ! Derivative wrt omega of the Self-energy
 integer,parameter :: unt_sgm = 20  ! Sigma on the Matsubara axis

!nb of bytes related to GW arrays, that can be tuned from sp to dp independently
!of other variables in ABINIT. Presently single precision is the default.
!#if defined HAVE_GW_DPC
! integer, parameter :: gwp=kind(1.0_dp)
! integer, parameter :: gwpc=kind((1.0_dp,1.0_dp))
!#else
! integer, parameter :: gwp=kind(1.0)
! integer, parameter :: gwpc=kind((1.0,1.0))
!#endif

 real(dp),parameter :: GW_TOLQ =0.0001_dp  
 ! Tolerance below which two BZ points are considered equal within a RL vector:
 ! for each red. direct. the abs value of the difference btw the two coord must be smaller that tolq.

 real(dp),parameter :: GW_TOLQ0=0.001_dp   
 ! Tolerance below which a q-point is treated as zero (long wavelength limit)

 real(dp),parameter :: GW_TOL_DOCC=0.01_dp 
 ! Tolerance on the difference between two occupation numbers.
 ! below this value, the contribution of the transition is neglected in the evaluation of chi0

 real(dp),parameter :: GW_TOL_W0=0.001_dp 
 ! Tolerance on the real part of the frequency appearing in the denominator of the 
 ! non-interacting Green function G0. Above this value, a small purely imaginary 
 ! complex shift is added to the denominator during the evaluation of chi0.

 complex(gwpc),parameter :: czero_gw=(0.0_gwp,0.0_gwp)
 complex(gwpc),parameter :: cone_gw =(1.0_gwp,0.0_gwp)
 complex(gwpc),parameter :: j_gw    =(0.0_gwp,1.0_gwp)
 complex(dpc) ,parameter :: j_dpc   =(0.0_dp ,1.0_dp )

CONTAINS  !===========================================================

 subroutine print_kinds(unit)

 !Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 integer,intent(in) :: unit

 !Local variables-------------------------------
 character(len=500) :: msg
 ! *********************************************************************

 write(msg,'(2a)')ch10,' DATA TYPE INFORMATION: '
 call wrtout(unit,msg,'COLL') 

 write(msg,'(4a,2(a,i6,a),2(a,e14.8,a),e14.8,a)')ch10,&
& ' REAL: Data type name:',' DP ',ch10,&
& '       Kind value: ',KIND(0.0_dp),ch10,&
& '       Precision: ',PRECISION(0.0_dp),ch10,&
& '       Smallest nonnegligible quantity relative to 1:',EPSILON(0.0_dp),ch10,&
& '       Smallest positive number:',TINY(0.0_dp),ch10,&
& '       Largest representable number:',HUGE(0.0_dp)
 call wrtout(unit,msg,'COLL') 

 write(msg,'(3a,2(a,i20,a),a,i20)')&
  ' INTEGER: Data type name:',' (default) ',ch10,&
& '          Kind value:',KIND(0),ch10,&
& '          Bit size:',BIT_SIZE(0),ch10,&
  '          Largest representable number:',HUGE(0)
 call wrtout(unit,msg,'COLL') 

 write(msg,'(3a,a,i6)')&
& ' LOGICAL: Data type name:',' (default) ',ch10,&
& '         Kind value:',KIND(.TRUE.)
 call wrtout(unit,msg,'COLL') 

 write(msg,'(2a,a,i6)')&
& ' CHARACTER: Data type name:',' (default) ',ch10,&
& '            Kind value:',KIND('C')
 call wrtout(unit,msg,'COLL') 

 end subroutine print_kinds

end module m_gwdefs
!!***
