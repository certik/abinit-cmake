!{\src2tex{textfont=tt}}
!!****f* ABINIT/test_ftgkk
!! NAME
!! test_ftgkk
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVerstra)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with matrix elements
!!   gprim = reciprocal lattice vectors
!!   natom = number of atoms
!!   nrpt = number of real space points for FT interpolation
!!   rpt = coordinates of real space points for FT interpolation
!!   spqpt = qpoint coordinates
!!   wghatm = weights for pairs of atoms in FT interpolation
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!  MJV 18/5/2008 reverted to old syntax/use for ftgkk, with all ft being done
!!   in a batch. Might come back to 5.5 version with atomic FT in ftgkk, but later.
!!
!! PARENTS
!!
!! CHILDREN
!!      ftgkk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine test_ftgkk(elph_ds,gprim,natom,nrpt,rpt,spqpt,wghatm)

 use defs_basis
  use defs_datatypes
  use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_17ddb, except_this_one => test_ftgkk
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),spqpt(3,elph_ds%nqpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,iqpt,isppol,qtor
!arrays
 real(dp),allocatable :: gkq_disk(:,:,:,:,:),tmp_gkq(:,:,:,:,:)

! *************************************************************************

!for each qpt do FT to recuperate original values

 isppol = 1
 qtor = 0
 allocate (gkq_disk(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol))
 allocate (tmp_gkq(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)) 

 do iqpt=1,elph_ds%nqpt
  tmp_gkq(:,:,:,:,:) = zero

  call ftgkk (wghatm,tmp_gkq,elph_ds%gkk_rpt,elph_ds%gkqwrite,&
&  elph_ds%gkk_rptwrite,gprim,1,natom,&
&  elph_ds%nFSkpt,elph_ds%ngkkband,elph_ds%nFSkpt,1,&
&  nrpt,elph_ds%nsppol,qtor,rpt,spqpt,elph_ds%unit_gkk_rpt,elph_ds%unitgkq)

  if (elph_ds%gkqwrite == 0) then
   do iFSkpt=1,10
    write (93,*) tmp_gkq(:,:,:,iFSkpt,isppol)-elph_ds%gkk_qpt(:,:,:,iFSkpt,isppol,iqpt)
   end do
  else
   read (elph_ds%unitgkq,REC=iqpt) gkq_disk
   do iFSkpt=1,10
    write (93,*) tmp_gkq(:,:,:,iFSkpt,isppol)-gkq_disk(:,:,:,iFSkpt,isppol)
   end do
  end if
 end do

end subroutine test_ftgkk
!!***
