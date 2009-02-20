!{\src2tex{textfont=tt}}
!!****f* ABINIT/bstruct_clean
!! NAME
!! bstruct_clean
!!
!! FUNCTION
!! This subroutine deallocates the components of the bandstructure structured datatype
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  (only deallocate)
!!
!! OUTPUT
!!  (only deallocate)
!!
!! SIDE EFFECTS
!!  bstruct<type(bstruct_type)>=the band structure
!!
!! PARENTS
!!      gstate,loper3,newsp,nonlinear,respfn
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine bstruct_clean(bstruct)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(bandstructure_type),intent(inout) :: bstruct

!Local variables-------------------------------

! *************************************************************************

!DEBUG
!write(6,*)' bstruct_clean : enter'
!stop
!ENDDEBUG

!Deallocate all components of bstruct
 deallocate(bstruct%istwfk)
 deallocate(bstruct%nband)
 deallocate(bstruct%npwarr)
 deallocate(bstruct%kptns)
 deallocate(bstruct%eig)
 deallocate(bstruct%occ)
 deallocate(bstruct%doccde)
 deallocate(bstruct%wtk)

end subroutine bstruct_clean
!!***
