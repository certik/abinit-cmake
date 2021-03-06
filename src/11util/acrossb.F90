!{\src2tex{textfont=tt}}
!!****f* ABINIT/acrossb
!! NAME
!! acrossb
!!
!! FUNCTION
!! Calculates the cross product of two 3-vectors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   a(3): real(dp) vector
!!   b(3): real(dp) vector
!!
!! OUTPUT
!!   c(3): real(dp) vector = a X b
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine acrossb(a,b,c)

 use defs_basis

 implicit none

!Arguments ---------------------------------------------
!arrays
 real(dp),intent(in) :: a(3),b(3)
 real(dp),intent(out) :: c(3)

!Local variables ---------------------------------------

! *********************************************************************
 
 c(1) =  a(2)*b(3) - a(3)*b(2)
 c(2) = -a(1)*b(3) + a(3)*b(1)
 c(3) =  a(1)*b(2) - b(1)*a(2)

end subroutine acrossb

!!***
