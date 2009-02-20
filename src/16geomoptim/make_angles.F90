!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_angles
!! NAME
!! make_angles
!!
!! FUNCTION
!!  (to be completed)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      make_prim_internals
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine make_angles(angs,bonds,icenter,nang,nbond,dtset)

 use defs_basis
  use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter,nbond
 integer,intent(out) :: nang
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,pointer :: angs(:,:,:),bonds(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ia1,ia2,iang,iatom,ibond,ii,is1,is2,ishift,ja1,ja2,jatom
 integer :: jbond,js1,js2
!arrays
 integer,allocatable :: angs_tmp(:,:,:)
 integer,allocatable,target :: angs_tmp2(:,:,:)

! *************************************************************************

!DEBUG
!write (*,*) 'make_angs: enter'
!ENDDEBUG


!tentative first allocation: < 6 angles per bond.
 allocate (angs_tmp(2,3,72*dtset%natom))

 nang = 0

 do ibond=1,nbond
  ia1 = bonds(1,1,ibond)
  is1 = bonds(2,1,ibond)
  ia2 = bonds(1,2,ibond)
  is2 = bonds(2,2,ibond)
  do jbond=ibond+1,nbond
   ja1 = bonds(1,1,jbond)
   ja2 = bonds(1,2,jbond)
   do ishift=-(icenter-1),+(icenter-1)
    js1 = bonds(2,1,jbond)+ishift
    js2 = bonds(2,2,jbond)+ishift

    if      (ia1==ja1 .and. is1==js1) then
     nang = nang+1
     angs_tmp(:,1,nang) = (/ia2,is2/)
     angs_tmp(:,2,nang) = (/ia1,is1/)
     angs_tmp(:,3,nang) = (/ja2,js2/)

    else if (ia1==ja2 .and. is1==js2) then
     nang = nang+1
     angs_tmp(:,1,nang) = (/ia2,is2/)
     angs_tmp(:,2,nang) = (/ia1,is1/)
     angs_tmp(:,3,nang) = (/ja1,js1/)

    else if (ia2==ja2 .and. is2==js2) then
     nang = nang+1
     angs_tmp(:,1,nang) = (/ia1,is1/)
     angs_tmp(:,2,nang) = (/ia2,is2/)
     angs_tmp(:,3,nang) = (/ja1,js1/)

    else if (ia2==ja1 .and. is2==js1) then
     nang = nang+1
     angs_tmp(:,1,nang) = (/ia1,is1/)
     angs_tmp(:,2,nang) = (/ia2,is2/)
     angs_tmp(:,3,nang) = (/ja2,js2/)

    end if
    if (nang > 72*dtset%natom) then
     write (*,*) 'make_angles : too many angles found > 72*dtset%natom'
     stop
    end if
   end do
  end do
! end jbond do
 end do
!end ibond do

 allocate (angs(2,3,nang))
 do iang=1,nang
  angs(:,:,iang) = angs_tmp(:,:,iang)
 end do
 deallocate (angs_tmp)
!angs => angs_tmp2

end subroutine make_angles
!!***
