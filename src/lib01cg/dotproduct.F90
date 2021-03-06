!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotproduct
!! NAME
!! dotproduct
!!
!! FUNCTION
!! scalar product of two vectors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! v1 and v2: two real(dp) vectors
!!
!! OUTPUT
!! scalar product of the two vectors
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! vector size is not checked
!!
!! NOTES
!! I've benchmarked this to be speedier than the intrinsic dot_product even on
!! big vectors. The point is that less check is performed.
!!
!! PARENTS
!! cgpr,brent
!!
!! CHILDREN
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function dotproduct(nv1,nv2,v1,v2)

 use defs_basis
  use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nv1,nv2
 real(dp) :: dotproduct
!arrays
 real(dp),intent(in) :: v1(nv1,nv2),v2(nv1,nv2)

!Local variables-------------------------------
!scalars
 integer :: i,j

! *************************************************************************
 dotproduct=zero
 do j=1,nv2
  do i=1,nv1
   dotproduct=dotproduct+v1(i,j)*v2(i,j)
  end do
 end do
end function dotproduct

!!***
