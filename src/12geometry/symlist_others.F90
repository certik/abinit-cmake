!{\src2tex{textfont=tt}}
!!****f* ABINIT/symlist_others
!! NAME
!! symlist_others
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!! Non primitive, non BCC, non FCC case : rhombohedral or face centered
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! additional_info=information that is needed beyond n_axes, in order
!!  to discriminate between specific space groups
!! brvltt=Bravais lattice type
!! nsym=actual number of symmetries
!! n_axes(31)=array containing the number of all the possible symmetry operations
!!
!! OUTPUT
!! spgroup=space group number ; returns 0 if not found
!!
!! NOTES
!!
!! The list of symmetry operations is for the conventional cell
!!
!! TODO
!! For the time being there are several groups where uncertainties still exist
!! This will be solved in the very next ABINIT version
!!
!! PARENTS
!!      symspgr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symlist_others(additional_info,brvltt,nsym,n_axes,spgroup)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: additional_info,brvltt,nsym
 integer,intent(out) :: spgroup
!arrays
 integer,intent(in) :: n_axes(31)

!Local variables-------------------------------
!character(len=500) :: message
!arrays
 integer :: n_axest(31)

!**************************************************************************

!DEBUG
!write(6,*) ' symlist_others : enter '
!write(6,*) ' nsym = ', nsym
!write(6,*) ' brvltt = ',brvltt
!write(6, '(a,10i3)' ) ' n_axes(1:10) =',n_axes(1:10)
!write(6, '(a,10i3)' ) ' n_axes(11:20)=',n_axes(11:20)
!write(6, '(a,11i3)' ) ' n_axes(21:31)=',n_axes(21:31)
!ENDDEBUG

 spgroup=0

 if(brvltt==4 .or. brvltt==5 .or. brvltt==6)then       !  A, B, C face centered

  select case(nsym)

   case(4)

    n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,0,0,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=5
    n_axest=(/0,0,0,0,0,0,1,1,0,0,  0,0,0,0,1,1,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=8
    n_axest=(/0,0,0,0,0,0,1,1,0,0,  0,0,0,0,0,1,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=9

   case(8)

    n_axest=(/0,0,0,0,0,0,1,1,0,0,  0,0,0,0,1,2,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=36

    n_axest=(/0,0,0,0,2,0,1,1,1,0,  0,0,0,0,1,1,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=12
    n_axest=(/0,0,0,0,2,0,1,1,1,0,  0,0,0,0,0,1,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=15
    n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,2,1,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=38
    n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,1,3,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=39
    n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,1,1,0,2,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=40
    n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,0,3,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=41

    n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,0,0,0,0,0,4,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=20
    n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,2,2,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=35
    n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,0,2,0,2,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=37

    n_axest=(/0,0,0,0,0,0,1,1,4,0,  0,0,0,0,0,0,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=21

   case(16)

    n_axest=(/0,0,0,0,2,0,1,1,2,0,  0,0,0,0,2,2,0,2,0,4,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=63
    n_axest=(/0,0,0,0,2,0,1,1,2,0,  0,0,0,0,1,4,0,1,0,4,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=64

    n_axest=(/0,0,0,0,2,0,1,1,4,0,  0,0,0,0,3,2,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=65
    n_axest=(/0,0,0,0,2,0,1,1,4,0,  0,0,0,0,2,4,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=67
    n_axest=(/0,0,0,0,2,0,1,1,4,0,  0,0,0,0,0,4,0,2,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=68
    n_axest=(/0,0,0,0,2,0,1,1,4,0,  0,0,0,0,1,2,0,3,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=66

  end select

 else if(brvltt==7)then                ! Rhombohedral lattice

  select case(nsym)

   case(3)

    n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=146

   case(6)

    n_axest=(/0,0,2,0,1,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=148
    n_axest=(/0,0,0,0,0,0,0,1,3,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=155
    n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,3,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=160
    n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,3,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=161

   case(12)

    n_axest=(/0,0,2,0,1,0,0,1,3,2,  0,0,0,0,3,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=166
    n_axest=(/0,0,2,0,1,0,0,1,3,2,  0,0,0,0,0,3,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
    if(sum((n_axes-n_axest)**2)==0) spgroup=167

  end select

 end if

end subroutine symlist_others
!!***
