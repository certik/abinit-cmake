!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_clean
!! NAME
!! hdr_clean
!!
!! FUNCTION
!! This subroutine deallocates the components of the header structured datatype
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! hdr <type(hdr_type)>=the header
!!
!! OUTPUT
!!  (only deallocate)
!!
!! PARENTS
!!      anascr,compare_interpol,conducti_nc,conducti_paw,cut3d,elphon,getgsc
!!      gstate,initaim,inwffil,ioarr,loper3,macroave,mrggkk,mrgscr,newsp
!!      nonlinear,optic,rdscr,read_el_veloc,read_gkk,respfn,rhoij_free,screening,sigma
!!      suscep,testlda,testscr,wannier
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine hdr_clean(hdr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(hdr_type),intent(inout) :: hdr

!Local variables-------------------------------

! *************************************************************************

!DEBUG
!write(6,*)' hdr_clean : enter'
!stop
!ENDDEBUG

!Deallocate all components of hdr
 deallocate(hdr%istwfk)
 deallocate(hdr%kptns)
 deallocate(hdr%lmn_size)
 deallocate(hdr%nband)
 deallocate(hdr%npwarr)
 deallocate(hdr%occ)
 deallocate(hdr%pspcod)
 deallocate(hdr%pspdat)
 deallocate(hdr%pspso)
 deallocate(hdr%pspxc)
 deallocate(hdr%so_psp)
 deallocate(hdr%symafm)
 deallocate(hdr%symrel)
 deallocate(hdr%title)
 deallocate(hdr%tnons)
 deallocate(hdr%typat)
 deallocate(hdr%wtk)
 deallocate(hdr%xred)
 deallocate(hdr%zionpsp)
 deallocate(hdr%znuclpsp)
 deallocate(hdr%znucltypat)
 if (hdr%usepaw==1) then
  call rhoij_free(hdr%pawrhoij)
  deallocate(hdr%pawrhoij)
 end if

!DEBUG
!write(6,*)' hdr_clean : exit'
!ENDDEBUG

end subroutine hdr_clean
!!***
