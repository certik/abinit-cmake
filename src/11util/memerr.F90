!{\src2tex{textfont=tt}}
!!****f* ABINIT/memerr
!! NAME
!! memerr
!!
!! FUNCTION
!!  This function deals with the failure of a memory allocation.
!!  It writes a message reporting the name of the subroutine where 
!!  the allocation has failed and the dimension of the array.
!!  It uses leave_now to stop the code smoothly
!! 
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  sub_name=name of the subroutine where the allocation has failed
!!  array_name=name of the array that has not been allocated
!!  nelements=number of elements in the array
!!  kindp=string defining the kind of the array 
!!   sp  for single precision real
!!   spc for single precision complex
!!   dp  for double-precision real
!!   dpc for double-precision complex 
!!     
!! OUTPUT
!!  Only write, then stop the code 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      anascr,cchi0,cchi0q0,cppm2par,csigme,fftwfn,findggp,findnq
!!      make_epsm1_driver,mkphdos,mrgscr,outkss,rdm,rdscr,screening,sigma
!!      wf_info
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine memerr(sub_name,array_name,nelements,kindp)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nelements
 character(len=*),intent(in) :: array_name,kindp,sub_name

!Local variables-------------------------------
!scalars
 integer :: ib
 real(dp) :: sizeMb
 character(len=500) :: msg

! *************************************************************************
 
 select case (TRIM(kindp))
  case ('sp') 
   ib=KIND(1.0) 
  case ('spc')
   ib=KIND((1.0d0,1.0d0))
  case ('dp') 
   ib=dp
  case ('dpc') 
   ib=2*dpc 
  case ('gwpc') 
   ib=2*gwpc
  case ('i1b')
   ib=i1b
  case ('i2b')
   ib=i2b
  case ('i4b')
   ib=i4b
   case default
   write(msg,'(5a)')ch10,&
&   ' memerr: ERROR - ',ch10,&
&   ' called with a wrong value of kindp: ',TRIM(kindp)
   call wrtout(std_out,msg,'PERS') 
   call leave_new('COLL')
 end select
 write(*,*)ib,nelements
 sizeMb=(DBLE(nelements)/1024.0_dp**2)
 sizeMb=ib*sizeMb 

 write(msg,'(1x,8a,f12.2,a)')ch10,&
& TRIM(sub_name),': ERROR - ',ch10,&
& ' out of memory in ',TRIM(array_name),ch10,&
& ' requiring ',sizeMb,' Mb.'
 call wrtout(std_out,msg,'PERS') 
 call wrtout(ab_out,msg,'PERS')
 call leave_new('COLL')

end subroutine memerr
!!***
