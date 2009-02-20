!{\src2tex{textfont=tt}}
!!****f* ABINIT/vdotw
!! NAME
!! vdotw
!!
!! FUNCTION
!! Compute the scalar product of two vectors expressed in reduced coordinates
!! The result is multiplied by (2pi)**2 in case of vectors in reciprocal space
!! to take into account the correct normalisation of the reciprocal lattice vectors
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  v(3),w(3)=vectors in reduced coordinates
!!  met(3,3)=metric tensor
!!  space=character defining whether we are working in real (r) or reciprocal space (g)
!! 
!!
!! OUTPUT
!!  vdotw=scalar product of v and w  
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

function vdotw(v,w,met,space)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: vdotw
 character(len=1),intent(in) :: space
!arrays
 real(dp),intent(in) :: met(3,3),v(3),w(3)

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************
 
 vdotw = (  met(1,1)* v(1)*w(1)              &
& +met(2,2)* v(1)*w(2)              &
& +met(3,3)* v(3)*w(3)              &
& +met(1,2)*(v(1)*w(2) + v(2)*w(1)) &
& +met(1,3)*(v(1)*w(3) + v(3)*w(1)) &
& +met(2,3)*(v(2)*w(3) + v(3)*w(2)) ) 

 select case (space)
  case ('r')
   return 
  case ('g')
   vdotw= two_pi**2 * vdotw
   case default
   write(message,'(5a)')ch10,&
&   ' vdotw : BUG -',ch10,   &
&   ' wrong value for space = ',space 
   call wrtout(06,message,'COLL') 
   call leave_new('COLL')
 end select

end function vdotw
!!***
