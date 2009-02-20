!{\src2tex{textfont=tt}}
!!****f* ABINIT/normv
!! NAME
!! normv
!!
!! FUNCTION
!! Compute the norm of a vector expressed in reduce coordinates
!! The result is multiplied by 2pi in case of a vector in reciprocal space
!! to take into account the correct normalisation of the reciprocal lattice vectors
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  v(3)=vector in reduced coordinates
!!  met(3,3)=metric tensor
!!  space=character defining whether we are working in real (r) or reciprocal space (g)
!!
!! OUTPUT
!!  normv=norm of v 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  
!!
!! CHILDREN
!! 
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function normv(v,met,space)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: normv
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3),v(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 normv= (        v(1)*met(1,1)*v(1) + v(2)*met(2,2)*v(2) + v(3)*met(3,3)*v(3)   &
& +two * (v(1)*met(1,2)*v(2) + v(1)*met(1,3)*v(3) + v(2)*met(2,3)*v(3))  )

 select case (TRIM(space)) 
  case ('r','R')
   normv=SQRT(normv)
  case ('g','G')
   normv=two_pi*SQRT(normv)
   case default
   write(msg,'(5a)')ch10,&
&   ' normv : BUG -',ch10,   &
&   ' wrong value for space = ',space 
   call wrtout(std_out,msg,'COLL') 
   call leave_new('COLL')
 end select
 
end function normv
!!***
