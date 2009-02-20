!{\src2tex{textfont=tt}}
!!****f* ABINIT/factorial
!! NAME
!! factorial
!!
!! FUNCTION
!! Calculates N! . Returns a (dp) real.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   nn=number to use
!!
!! OUTPUT
!!   factorial= n! (real)
!!
!! PARENTS
!!      gaunt,setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function factorial(nn)

 use defs_basis

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: nn
 real(dp) :: factorial

!Local variables ---------------------------------------
!scalars
 integer :: ii
 real(dp) :: ff

! *********************************************************************

 ff=one
 do ii=2,nn
  ff=ff*ii
 end do

 factorial=ff

 end function factorial

!!***
