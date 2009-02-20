!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_free_type_wfs
!!
!! NAME
!! wvl_free_type_wfs
!!
!! FUNCTION
!! Freeing routine.
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
!!
!! SIDE EFFECTS
!!  wfs <type(wvl_wf_type)>=wavefunctions informations in a wavelet basis.
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_free_type_wfs(wfs)

 use defs_basis
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

  !Arguments ------------------------------------
  !scalars
  type(wvl_wf_type),intent(inout) :: wfs
  !Local variables -------------------------
  character(len=500)       :: message
  ! *********************************************************************
#if defined HAVE_BIGDFT
 if (associated(wfs%keys%keyg) .and. associated(wfs%keys%keyv)) then
  call deallocate_wfd(wfs%keys, "wvl_free_type_wfs")
 end if
 call deallocate_bounds(wfs%bounds, "wvl_free_type_wfs")
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_free_type_wfs: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif
 if (associated(wfs%psi)) deallocate(wfs%psi)
 if (associated(wfs%hpsi)) deallocate(wfs%hpsi)
 if (associated(wfs%eval)) deallocate(wfs%eval)
 deallocate(wfs%spinar)
 if (associated(wfs%psit))deallocate(wfs%psit)
 if (associated(wfs%psidst)) deallocate(wfs%psidst)
 if (associated(wfs%hpsidst)) deallocate(wfs%hpsidst)
 if (associated(wfs%ads)) deallocate(wfs%ads)
end subroutine wvl_free_type_wfs
!!***

!!****f* ABINIT/wvl_free_type_proj
!!
!! NAME
!! wvl_free_type_proj
!!
!! FUNCTION
!! Freeing routine.
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
!!
!! SIDE EFFECTS
!!  proj <type(wvl_projectors_type)>=projectors informations in a wavelet basis.
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_free_type_proj(proj)

 use defs_basis
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

  !Arguments ------------------------------------
  !scalars
  type(wvl_projectors_type),intent(inout) :: proj
  !Local variables -------------------------
  character(len=500)       :: message
  ! *********************************************************************
#if defined HAVE_BIGDFT
 deallocate(proj%keys%nvctr_p)
 deallocate(proj%keys%nseg_p)
 deallocate(proj%keys%keyv_p)
 deallocate(proj%keys%keyg_p)
 deallocate(proj%keys%nboxp_c, proj%keys%nboxp_f)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_free_type_proj: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif
 deallocate(proj%proj)
 if (associated(proj%der)) deallocate(proj%der)
end subroutine wvl_free_type_proj
!!***
