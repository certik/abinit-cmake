!{\src2tex{textfont=tt}}
!!****f* ABINIT/energies_init
!!
!! NAME
!! energies_init
!!
!! FUNCTION
!! Set zero in all values of a type(energies_type) object.
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
!! OUTPUT
!!   energies <type(energies_type)>=values to initialise
!!
!! PARENTS
!!      driver,gstate,scfcv,sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine energies_init(energies)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(energies_type),intent(out) :: energies

!Local variables-------------------------------

! *************************************************************************

 energies%entropy       = zero 
 energies%e_entropy     = zero 
 energies%e_fermie      = zero
 energies%e_paw         = zero
 energies%e_pawdc       = zero
 energies%e_kinetic     = zero
 energies%e_localpsp    = zero
 energies%e_nonlocalpsp = zero
 energies%e_eigenvalues = zero
 energies%e_hartree     = zero
 energies%e_ewald       = zero
 energies%e_xc          = zero
 energies%e_vxc         = zero
 energies%e_xcdc        = zero
 energies%e_elecfield   = zero
 energies%e_corepsp     = zero
 energies%e_pulay       = zero

end subroutine energies_init
!!***
