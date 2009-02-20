!{\src2tex{textfont=tt}}
!!****f* ABINIT/dtsetfree
!! NAME
!! dtsetfree
!!
!! FUNCTION
!! Free a dataset after use.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MF, GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SIDE EFFECTS
!!  dtset <type(dataset_type)>=free all associated pointers.
!!
!! PARENTS
!!      abinit,afterscfloop,cvxclda,driver,kxc_alda,outkss,wvl_memory,xc_kernel
!!
!! CHILDREN
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dtsetFree(dtset)

 use defs_basis
  use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(inout) :: dtset

!Local variables-------------------------------

! *************************************************************************

 if (associated(dtset%algalch)) then
  deallocate(dtset%algalch)
 end if
 if (associated(dtset%atvshift)) then
  deallocate(dtset%atvshift)
 end if
 if (associated(dtset%bdgw)) then
  deallocate(dtset%bdgw)
 end if
 if (associated(dtset%corecs)) then
  deallocate(dtset%corecs)
 end if
 if (associated(dtset%iatfix)) then
  deallocate(dtset%iatfix)
 end if
 if (associated(dtset%iatsph)) then
  deallocate(dtset%iatsph)
 end if
 if (associated(dtset%istwfk)) then
  deallocate(dtset%istwfk)
 end if
 if (associated(dtset%kberry)) then
  deallocate(dtset%kberry)
 end if
 if (associated(dtset%lexexch)) then
  deallocate(dtset%lexexch)
 end if
 if (associated(dtset%lpawu)) then
  deallocate(dtset%lpawu)
 end if
 if (associated(dtset%ltypeorb)) then
  deallocate(dtset%ltypeorb)
 end if
 if (associated(dtset%nband)) then
  deallocate(dtset%nband)
 end if
 if (associated(dtset%numorb)) then
  deallocate(dtset%numorb)
 end if
 if (associated(dtset%qmass)) then
  deallocate(dtset%qmass)
 end if
 if (associated(dtset%so_psp)) then
  deallocate(dtset%so_psp)
 end if
 if (associated(dtset%symafm)) then
  deallocate(dtset%symafm)
 end if
 if (associated(dtset%symrel)) then
  deallocate(dtset%symrel)
 end if
 if (associated(dtset%typat)) then
  deallocate(dtset%typat)
 end if
 if (associated(dtset%amu)) then
  deallocate(dtset%amu)
 end if
 if (associated(dtset%densty)) then
  deallocate(dtset%densty)
 end if
 if (associated(dtset%dmatpawu)) then
  deallocate(dtset%dmatpawu)
 end if
 if (associated(dtset%jpawu)) then
  deallocate(dtset%jpawu)
 end if
 if (associated(dtset%kpt)) then
  deallocate(dtset%kpt)
 end if
 if (associated(dtset%kptgw)) then
  deallocate(dtset%kptgw)
 end if
 if (associated(dtset%kptns)) then
  deallocate(dtset%kptns)
 end if
 if (associated(dtset%mixalch)) then
  deallocate(dtset%mixalch)
 end if
 if (associated(dtset%occ_orig)) then
  deallocate(dtset%occ_orig)
 end if
 if (associated(dtset%ptcharge)) then
  deallocate(dtset%ptcharge)
 end if
 if (associated(dtset%quadmom)) then
  deallocate(dtset%quadmom)
 end if
 if (associated(dtset%qptdm)) then
  deallocate(dtset%qptdm)
 end if
 if (associated(dtset%ratsph)) then
  deallocate(dtset%ratsph)
 end if
 if (associated(dtset%rcoord)) then
  deallocate(dtset%rcoord)
 end if
 if (associated(dtset%rtheta)) then
  deallocate(dtset%rtheta)
 end if
 if (associated(dtset%shiftk)) then
  deallocate(dtset%shiftk)
 end if
 if (associated(dtset%spinat)) then
  deallocate(dtset%spinat)
 end if
 if (associated(dtset%tnons)) then
  deallocate(dtset%tnons)
 end if
 if (associated(dtset%upawu)) then
  deallocate(dtset%upawu)
 end if
 if (associated(dtset%vel_orig)) then
  deallocate(dtset%vel_orig)
 end if
 if (associated(dtset%wtatcon)) then
  deallocate(dtset%wtatcon)
 end if
 if (associated(dtset%wtk)) then
  deallocate(dtset%wtk)
 end if
 if (associated(dtset%w90lplot)) then
  deallocate(dtset%w90lplot)
 end if
 if (associated(dtset%xred_orig)) then
  deallocate(dtset%xred_orig)
 end if
 if (associated(dtset%ziontypat)) then
  deallocate(dtset%ziontypat)
 end if
 if (associated(dtset%znucl)) then
  deallocate(dtset%znucl)
 end if

end subroutine dtsetFree
!!***
