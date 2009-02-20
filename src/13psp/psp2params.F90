!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp2params_init
!! NAME
!! psp2params_init
!!
!! FUNCTION
!! Allocate and initialise the data structure holding parameters for the GTH
!! pseudo-potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npsp=number of true pseudo used (not alchemy).
!!
!! OUTPUT
!!  gth_params <type (pseudopotential_gth_type)>=the values to allocate and initialise.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine psp2params_init(gth_params, npsp)

 use defs_basis
  use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npsp
 type(pseudopotential_gth_type),intent(out) :: gth_params

!Local variables-------------------------------

! *********************************************************************
!Check array, no params are currently set.
 allocate(gth_params%set(npsp))
 gth_params%set(:) = .false.

!Check array, have geometric informations been filled?
 allocate(gth_params%hasGeometry(npsp))
 gth_params%hasGeometry(:) = .false.

!Coefficients for local part and projectors
 allocate(gth_params%psppar(0:4, 0:6, npsp))
 gth_params%psppar = real(0, dp)

!Different radii
 allocate(gth_params%radii_cov(npsp))
 allocate(gth_params%radii_cf(npsp, 2))
end subroutine psp2params_init
!!***

!!****f* ABINIT/psp2params_free
!! NAME
!! psp2params_free
!!
!! FUNCTION
!! Deallocate a previously allocated data structure for storage of GTH parameters.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! SIDE EFFECTS
!!  gth_params <type (pseudopotential_gth_type)>=the values to deallocate.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine psp2params_free(gth_params)

 use defs_basis
  use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(pseudopotential_gth_type),intent(inout) :: gth_params

!Local variables-------------------------------

! *********************************************************************

!Check arrays.
 deallocate(gth_params%set)
 deallocate(gth_params%hasGeometry)

!Coefficients for local part and projectors
 deallocate(gth_params%psppar)

!Different radii
 deallocate(gth_params%radii_cov)
 deallocate(gth_params%radii_cf)
end subroutine psp2params_free
!!***
