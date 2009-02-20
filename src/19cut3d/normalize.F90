!{\src2tex{textfont=tt}}
!!****f* ABINIT/normalize
!! NAME
!! normalize
!!
!! FUNCTION
!! Normalizes the value of v
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (see side effects)
!!
!! OUTPUT
!!  (see side effects)

!! SIDE EFFECTS
!!   v=value to be normalized
!!
!! PARENTS
!!      lineint,planeint,volumeint
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine normalize(v)

 use defs_basis

 implicit none

!Arguments-------------------------------------------------------------
!arrays
 real(dp),intent(inout) :: v(3)

!Local variables--------------------------------------------------------
!scalars
 integer :: idir
 real(dp) :: norm

! *************************************************************************

 norm=0.0
 do idir=1,3
  norm=norm+v(idir)**2
 end do
 norm=sqrt(norm)

 do idir=1,3
  v(idir)=v(idir)/norm
 end do

 return
end subroutine normalize
!!***
