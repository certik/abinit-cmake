!{\src2tex{textfont=tt}}
!!****f* ABINIT/contract_int_ge_val
!! NAME
!! contract_int_ge_val
!!
!! FUNCTION
!!  "Design by contract" routine, ensuring that the
!!  target integer argument is greater or equal to some limiting value.
!!  Might be used to test composite quantities :
!!  integer_name might be "arg1-arg2*arg3", with corresponding input value
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  character(len=*) calling_routine = name of the routine of which the argument is checked
!!  character(len=*) integer_name = name of the integer variable that is checked
!!  integer_value = value that is checked
!!  limit = lower admitted limit for the integer
!!
!! OUTPUT
!!  (only checking, writing and stopping)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  "Design by contract" routines are simpler than the routines
!!  that check the input variables. Indeed, the error message
!!  can be much primitive.
!!
!! PARENTS
!!      dotprod_g,dotprod_v,dotprod_vn,dotprodm_v,dotprodm_vn,matrixelmt_g
!!      mean_fftr,meanvalue_g,sqnorm_g,sqnorm_v,sqnormm_v
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine contract_int_ge_val(calling_routine,integer_name,integer_value,limit)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: integer_value,limit
 character(len=*),intent(in) :: calling_routine,integer_name

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(6,*)' contract_int_ge_val : enter '
!ENDDEBUG

 if(integer_value<limit)then
  write(message,'(10a,i6,5a,i6,a)') ch10,&
&  ' contract_int_ge_val: BUG -',ch10,&
&  '  For the routine ',trim(calling_routine),',',ch10,&
&  '  the value of "',trim(integer_name),'" should be greater or equal to ',limit,'.',ch10,&
&  '  However, ',trim(integer_name),' = ',integer_value,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!DEBUG
!write(6,*)' contract_int_ge_val : exit'
!stop
!ENDDEBUG

end subroutine contract_int_ge_val
!!***
