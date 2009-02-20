!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkqptequiv
!!
!! NAME
!! mkqptequiv
!!
!! FUNCTION
!! This routine determines the equivalence between 
!!   1) qpoints and fermi surface kpoints
!!   2) qpoints under symmetry operations
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   FSkpt = fermi surface kpoints
!!   FSkptirred = coordinates of irreducible FS kpoints
!!   nFSkpt = number of kpoints in the full FS set
!!   nFSkptirred = number of kpoints in the irreducible FS set
!!   nqpt = number of qpoints
!!   nsym = number of symmetries
!!   spqpt = qpoint coordinates
!!   symrec = reciprocal space symops
!!
!! OUTPUT
!!   FSfullpqtofull = mapping of k + q onto k' for k and k' in full BZ
!!   qpttoqpt(itim,isym,iqpt) = qpoint index which
!!     transforms to iqpt under isym and with time reversal itim.
!!
!! NOTES
!!   REMOVED 3/6/2008: much too large matrix, and not used at present
!!       FStoqpt = mapping of kpoint pairs (1 irreducible and 1 full) to qpoints
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      canon9,mkkptrank,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkqptequiv (FSfullpqtofull,FSkpt,FSkptirred,nFSkpt,nFSkptirred,nqpt,nsym,&
&   qpttoqpt,spqpt,symrec)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_17ddb, except_this_one => mkqptequiv
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!  integer,intent(out) :: FStoqpt(nFSkpt,nFSkptirred)
!scalars
 integer,intent(in) :: nFSkpt,nFSkptirred,nqpt,nsym
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(out) :: FSfullpqtofull(nFSkpt,nqpt),qpttoqpt(2,nsym,nqpt)
 real(dp),intent(in) :: FSkpt(3,nFSkpt),FSkptirred(3,nFSkptirred),spqpt(3,nqpt)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,iFSqpt,iqpt,isym,jFSkpt,symrankFSkpt
 real(dp) :: res
 character(len=500) :: message
!arrays
 integer,allocatable :: invrankFSkpt(:),rankFSkpt(:)
 real(dp) :: kpt(3),tmpkpt(3)

! *************************************************************************

 allocate (rankFSkpt(nFSkpt),invrankFSkpt(16000000))
 
 write (message,'(a)')' mkqptequiv : making rankFSkpt and invrankFSkpt'
 call wrtout(06,message,'COLL')

 call mkkptrank (FSkpt,nFSkpt,rankFSkpt,invrankFSkpt)

!write (*,*) 'mkqptequiv : FStoqpt'
!
!FStoqpt(:,:) = -999
!do iFSkpt=1,nFSkptirred
!! ============= NEW VERSION ==================
!do iqpt=1,nqpt
!!  tmpkpt = jkpt = ikpt + qpt
!tmpkpt(:) = FSkptirred(:,iFSkpt)+spqpt(:,iqpt)
!call canon9(tmpkpt(1),kpt(1),res)
!call canon9(tmpkpt(2),kpt(2),res)
!call canon9(tmpkpt(3),kpt(3),res)
!!  which kpt is it among the full FS kpts
!symrankFSkpt = int(8000000.0_dp*(kpt(1)+half+tol8) + &
!&   40000.0_dp*(kpt(2)+half+tol8) + &
!&   200.0_dp*(kpt(3)+half+tol8))
!if (symrankFSkpt > 16000000) then
!write (*,*) ' mkfskgrid : error : rank should be inferior to ', 16000000
!stop
!end if
!jFSkpt = invrankFSkpt(symrankFSkpt)
!if (jFSkpt == -1) then
!write (*,*) ' mkqptequiv : Error : looks like no kpoint equiv to k+q !!!'
!stop
!end if
!
!FStoqpt(jFSkpt,iFSkpt) = iqpt
!end do
!
!============= OLD VERSION ==================
!do jFSkpt=1,nFSkpt
!!   tmpkpt = qpt = jkpt - ikpt
!tmpkpt(:) = FSkpt(:,jFSkpt) - FSkptirred(:,iFSkpt)
!call canon9(tmpkpt(1),kpt(1),res)
!call canon9(tmpkpt(2),kpt(2),res)
!call canon9(tmpkpt(3),kpt(3),res)
!do iqpt=1,nqpt
!if (  abs(kpt(1)-spqpt(1,iqpt)) + &
!&abs(kpt(2)-spqpt(2,iqpt)) + &
!&abs(kpt(3)-spqpt(3,iqpt)) < tol6) then
!FStoqpt(jFSkpt,iFSkpt) = iqpt
!!DEBUG
!!write (*,*) 'mkqptequiv : FStoqpt'
!!write (*,*) iFSkpt,jFSkpt,FStoqpt(jFSkpt,iFSkpt)
!!ENDDEBUG
!exit
!end if
!end do
!end do
!============= END VERSIONS ==================
!end do
!
!write (message,'(a)')' mkqptequiv : FStoqpt made. Do FSfullpqtofull'
!call wrtout(06,message,'COLL')

 FSfullpqtofull(:,:) = -999

 do iFSkpt=1,nFSkpt
  do iqpt=1,nqpt
!  tmpkpt = jkpt = ikpt + qpt
   tmpkpt(:) = FSkpt(:,iFSkpt) + spqpt(:,iqpt)
   call canon9(tmpkpt(1),kpt(1),res)
   call canon9(tmpkpt(2),kpt(2),res)
   call canon9(tmpkpt(3),kpt(3),res)

!  ============= NEW VERSION ==================
!  which kpt is it among the full FS kpts
   symrankFSkpt = int(8000000.0_dp*(kpt(1)+half+tol8) + &
&   40000.0_dp*(kpt(2)+half+tol8) + &
&   200.0_dp*(kpt(3)+half+tol8))
   if (symrankFSkpt > 16000000) then
    write (*,*) ' mkfskgrid : error : rank should be inferior to ', 16000000
    stop
   end if
   jFSkpt = invrankFSkpt(symrankFSkpt)
   if (jFSkpt == -1) then
    write (*,*) ' mkqptequiv : Error : looks like no kpoint equiv to k+q !!!'
    stop
   end if
   FSfullpqtofull(iFSkpt,iqpt) = jFSkpt

!  ============= OLD VERSION ==================
!  do jFSkpt=1,nFSkpt
!  if (  abs(kpt(1)-FSkpt(1,jFSkpt)) + &
!  &abs(kpt(2)-FSkpt(2,jFSkpt)) + &
!  &abs(kpt(3)-FSkpt(3,jFSkpt)) < tol6) then
!  FSfullpqtofull(iFSkpt,iqpt) = jFSkpt
!  exit
!  end if
!  end do
!  ============= END VERSIONS ==================
  end do
 end do

 deallocate (rankFSkpt,invrankFSkpt)

 write(message,'(a)')' mkqptequiv : FSfullpqtofull made. Do qpttoqpt'
 call wrtout(06,message,'COLL')

 qpttoqpt(:,:,:) = 0
 do iFSqpt=1,nqpt
  do isym=1,nsym
   tmpkpt(:) = symrec(:,1,isym)*spqpt(1,iFSqpt) &
&   + symrec(:,2,isym)*spqpt(2,iFSqpt) &
&   + symrec(:,3,isym)*spqpt(3,iFSqpt)
   call canon9(tmpkpt(1),kpt(1),res)
   call canon9(tmpkpt(2),kpt(2),res)
   call canon9(tmpkpt(3),kpt(3),res)
   do iqpt=1,nqpt
    if ( abs(kpt(1)-spqpt(1,iqpt)) + &
&    abs(kpt(2)-spqpt(2,iqpt)) + &
&    abs(kpt(3)-spqpt(3,iqpt)) < tol6) then
     qpttoqpt(1,isym,iqpt) = iFSqpt
!    DEBUG
!    write (*,*) 0,isym,iFSqpt,iqpt
!    ENDDEBUG
     exit
    end if
   end do
   tmpkpt = -tmpkpt
   call canon9(tmpkpt(1),kpt(1),res)
   call canon9(tmpkpt(2),kpt(2),res)
   call canon9(tmpkpt(3),kpt(3),res)
   do iqpt=1,nqpt
    if (  abs(kpt(1)-spqpt(1,iqpt)) + &
&    abs(kpt(2)-spqpt(2,iqpt)) + &
&    abs(kpt(3)-spqpt(3,iqpt)) < tol6) then
     qpttoqpt(2,isym,iqpt) = iFSqpt
!    DEBUG
!    write (*,*) 1,isym,iFSqpt,iqpt
!    ENDDEBUG
     exit
    end if
   end do

  end do
 end do

!DEBUG
!write (124,*) FStoqpt
!write (124,*) FSfullpqtofull
!write (124,*) qpttoqpt
!ENDDEBUG

end subroutine mkqptequiv
!!***
