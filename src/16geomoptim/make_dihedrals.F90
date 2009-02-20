!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_dihedrals
!! NAME
!! make_dihedrals
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
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine make_dihedrals(angs,badangles,bonds,dihedrals,icenter,nbond,nang,ndihed,nrshift,dtset)

 use defs_basis
  use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter,nang,nbond,nrshift
 integer,intent(out) :: ndihed
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: badangles(nang)
 integer,pointer :: angs(:,:,:),bonds(:,:,:),dihedrals(:,:,:)

!Local variables-------------------------------
!scalars
 integer :: chkdihed,i1,i2,i3,ia1,ia2,ia3,iang,iatom,ibond,idihed,ii,is1,is2
 integer :: is3,ishift,ja1,ja2,ja3,jang,jatom,jbond,js1,js2,js3,maxshift
 integer :: minshift
!arrays
 integer,allocatable :: diheds_tmp(:,:,:)
 integer,allocatable,target :: diheds_tmp2(:,:,:)

! *************************************************************************

!DEBUG
!write (*,*) 'make_dihedrals: enter'
!ENDDEBUG


!tentative first allocation: < 6 dihedrals per angle.
 allocate (diheds_tmp(2,4,6*nang))

 ndihed = 0
 diheds_tmp(:,:,:) = 0

 do iang=1,nang
  if (badangles(iang) == 1) cycle
  ia1 = angs(1,1,iang)
  is1 = angs(2,1,iang)
  ia2 = angs(1,2,iang)
  is2 = angs(2,2,iang)
  ia3 = angs(1,3,iang)
  is3 = angs(2,3,iang)

  do jang=iang+1,nang
   if (badangles(jang) == 1) cycle
   ja1 = angs(1,1,jang)
   ja2 = angs(1,2,jang)
   ja3 = angs(1,3,jang)
   do ishift=-(icenter-1),(icenter-1)
    js1 = angs(2,1,jang)+ishift
    js2 = angs(2,2,jang)+ishift
    js3 = angs(2,3,jang)+ishift

    chkdihed=0
    if (ia2==ja1 .and. is2==js1) then
     if (ia1==ja2 .and. is1==js2) then
      ndihed = ndihed+1
      diheds_tmp(:,1,ndihed) = (/ia3,is3/)
      diheds_tmp(:,2,ndihed) = (/ia2,is2/)
      diheds_tmp(:,3,ndihed) = (/ja2,js2/)
      diheds_tmp(:,4,ndihed) = (/ja3,js3/)
      chkdihed=1
     else if (ia3==ja2 .and. is3==js2) then
      ndihed = ndihed+1
      diheds_tmp(:,1,ndihed) = (/ia1,is1/)
      diheds_tmp(:,2,ndihed) = (/ia2,is2/)
      diheds_tmp(:,3,ndihed) = (/ja2,js2/)
      diheds_tmp(:,4,ndihed) = (/ja3,js3/)
      chkdihed=1
     end if
    else if (ia2==ja3 .and. is2==js3) then
     if (ia1==ja2 .and. is1==js2) then
      ndihed = ndihed+1
      diheds_tmp(:,1,ndihed) = (/ia3,is3/)
      diheds_tmp(:,2,ndihed) = (/ia2,is2/)
      diheds_tmp(:,3,ndihed) = (/ja2,js2/)
      diheds_tmp(:,4,ndihed) = (/ja1,js1/)
      chkdihed=1
     else if (ia3==ja2 .and. is3==js2) then
      ndihed = ndihed+1
      diheds_tmp(:,1,ndihed) = (/ia1,is1/)
      diheds_tmp(:,2,ndihed) = (/ia2,is2/)
      diheds_tmp(:,3,ndihed) = (/ja2,js2/)
      diheds_tmp(:,4,ndihed) = (/ja1,js1/)
      chkdihed=1
     end if
    end if
    if (ndihed > 6*nang) then
     write (*,*) 'make_dihedrals : too many dihedrals found > 6*nang'
     stop
    end if
    if (chkdihed == 1) then
     if (   diheds_tmp(1,4,ndihed) == diheds_tmp(1,1,ndihed) .and.&
&     diheds_tmp(2,4,ndihed) == diheds_tmp(2,1,ndihed) ) then
      write (*,*) 'make_dihedrals : Bad dihedral was found: atom1 == atom4. Discarding.'
      diheds_tmp(:,:,ndihed) = 0
      ndihed = ndihed-1
     end if
    end if
   end do
  end do
! end jang do
 end do
!end iang do



!write (*,*) ' copy ndihed dihedrals'

 allocate (dihedrals(2,4,ndihed))
 do idihed=1,ndihed
  dihedrals(:,:,idihed) = diheds_tmp(:,:,idihed)
! minshift = minval(diheds_tmp(2,:,idihed))
! if (minshift <= 0) then
! dihedrals(2,:,idihed) = dihedrals(2,:,idihed)+minshift+1
! end if
! maxshift = maxval(diheds_tmp(2,:,idihed))
! if (maxshift > nrshift) then
! dihedrals(2,:,idihed) = dihedrals(2,:,idihed)-maxshift
! end if
! 
  minshift = minval(diheds_tmp(2,:,idihed))
  maxshift = maxval(diheds_tmp(2,:,idihed))
  if (minshift <= 0 .or. maxshift > nrshift) then
   write (*,*) ' make_dihedrals : Error : dihedral extends beyond '
   write (*,*) '  first neighboring unit cells ! '
   stop
  end if
 end do
 deallocate (diheds_tmp)
!dihedrals => diheds_tmp2

end subroutine make_dihedrals
!!***
