!{\src2tex{textfont=tt}}
!!****f* ABINIT/getwtk
!! NAME
!! getwtk
!!
!! FUNCTION
!! Routine called by the program optic
!! Presumes kpts are the irreducible ones of a good uniform grid
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (SSharma,MVer,VRecoules,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!      canon9,dgemv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine getwtk(kpt,nkpt,nsym,symrel,wtk)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments -----------------------------------------------
! in
! out
!scalars
 integer,intent(in) :: nkpt,nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: kpt(3,nkpt)
 real(dp),intent(out) :: wtk(nkpt)

!Local variables -----------------------------------------
!scalars
 integer :: ii,ikpt,istar,isym,itim,new,nkpt_tot
 real(dp) :: shift,timsign,tmp
!arrays
 integer :: nstar(nkpt)
 real(dp) :: dkpt(3),kptstar(3,2*nkpt*nsym),rsymrel(3,3,nsym),symkpt(3)
 real(dp) :: tsymkpt(3)

! *************************************************************************

 do isym=1,nsym
  rsymrel(:,:,isym) = dble(symrel(:,:,isym))
 end do

!for each kpt find star and accumulate nkpts
 do ikpt=1,nkpt
  write (*,*) ' getwtk : ikpt = ', ikpt
  nstar(ikpt) = 0
  kptstar(:,:) = zero
  do isym=1,nsym

   call dgemv('N',3,3,one,rsymrel(:,:,isym),3,kpt(:,ikpt),1,zero,symkpt,1)

!  is symkpt already in star?
   do itim=0,1
    timsign=one-itim*two
    tsymkpt(:) = timsign*symkpt(:)
    call canon9(tsymkpt(1),tmp,shift) ;  tsymkpt(1) = tmp
    call canon9(tsymkpt(2),tmp,shift) ;  tsymkpt(2) = tmp
    call canon9(tsymkpt(3),tmp,shift) ;  tsymkpt(3) = tmp
    new=1
    do istar=1,nstar(ikpt)
     dkpt(:) = abs(tsymkpt(:)-kptstar(:,istar))
     if ( sum(dkpt) < 1.0d-6) then
      new=0
      exit
     end if
    end do
    if (new==1) then
     nstar(ikpt) = nstar(ikpt)+1
     kptstar(:,nstar(ikpt)) = tsymkpt(:)
    end if
   end do

  end do
! end do nsym
! DEBUG
! write (*,*) ' getwtk : nstar = ', nstar(ikpt)
! write (*,*) ' getwtk : star = '
! write (*,*)  kptstar(:,1:nstar(ikpt))
! ENDDEBUG
 end do
!end do nkpt

 nkpt_tot = sum(nstar)
!write (*,*) ' getwtk : nkpt_tot = ', nkpt_tot
 do ikpt=1,nkpt
  wtk(ikpt) = dble(nstar(ikpt))/dble(nkpt_tot)
 end do

end subroutine getwtk
!!***
