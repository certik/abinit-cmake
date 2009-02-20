!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkfskgrid
!!
!! NAME
!! mkfskgrid
!!
!! FUNCTION
!! This routine sets up the full FS kpt grid by symmetry
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nFSkptirred = number of irreducible kpoints close to the FS
!!  nsym        = number of symmetries for the full system
!!  symrec      = reciprocal space symmetries (those for the kpts)
!!  FSkptirred  = coordinates of the irreducible kpoints close to the FS
!!  timrev      = 1 if time reversal symmetry is to be used
!!
!! OUTPUT
!!  nFSkpt           = full number of kpoints close to the FS
!!  tmpFSfulltoirred = temp array with correspondence between
!!                     full and irred kpts close to the FS
!!  tmpFSfulltofull  = temp array with correspondence between
!!                     full kpts close to the FS under different symmetries
!!  tmpFSkpt         = temp array with coordinates of the full
!!                     kpoints close to the FS
!!  FSirredwtk       = weights of the irreducible kpoints
!!  FSirredtofull    = indices of irred kpoints in full array
!!
!! NOTES
!!  WARNING: supposes kpt grid has full symmetry!! Not always true!!!
!!    but should be for Monkhorst-Pack, efficient grids.
!!    otherwise you get an error message in interpolate_gkk because
!!    an FS kpt can not be found in the gkk file.
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      canon9,leave_new,mkkptrank,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkFSkgrid (FSirredwtk,FSirredtofull,tmpFSfulltoirred,&
&          tmpFSfulltofull,tmpFSkpt,FSkptirred,nFSkptirred,nFSkpt,&
&          nsym,symrec,timrev)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_17ddb, except_this_one => mkFSkgrid
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nFSkptirred,nsym,timrev
 integer,intent(out) :: nFSkpt
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 integer,intent(out) :: FSirredtofull(nFSkptirred)
 integer,intent(out) :: tmpFSfulltofull(2,nsym,2*nFSkptirred*nsym)
 integer,intent(out) :: tmpFSfulltoirred(3,2*nFSkptirred*nsym)
 real(dp),intent(in) :: FSkptirred(3,nFSkptirred)
 real(dp),intent(out) :: FSirredwtk(nFSkptirred),tmpFSkpt(3,2*nFSkptirred*nsym)

!Local variables-------------------------------
!scalars
 integer :: ikpt1,ikpt2,isym,itim,new,symrankFSkpt
 real(dp) :: res,ss,timsign
 character(len=500) :: message
!arrays
 integer,allocatable :: invrankFSkpt(:),rankFSkpt(:)
 real(dp) :: kpt(3),redkpt(3)

! *************************************************************************

 if(timrev /= 1 .and. timrev /= 0)then
  write (message,'(4a)')ch10,&
&  ' mkfskgrid : BUG-',ch10,' timrev must be 1 or 0'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 FSirredtofull(:) = 0
 FSirredwtk(:) = zero
 nFSkpt=0 !zero k-points found

 do isym=1,nsym
  do itim=0,1

   timsign = one-two*itim

   do ikpt1=1,nFSkptirred
!   generate symmetrics of kpt ikpt1
    kpt(:) = timsign*(symrec(:,1,isym)*FSkptirred(1,ikpt1) + &
&    symrec(:,2,isym)*FSkptirred(2,ikpt1) + &
&    symrec(:,3,isym)*FSkptirred(3,ikpt1))
    call canon9(kpt(1),redkpt(1),res)
    call canon9(kpt(2),redkpt(2),res)
    call canon9(kpt(3),redkpt(3),res)
    new=1
!   is the kpt on the full grid (may have lower symmetry than full spgroup)
!   is kpt among the full FS kpts found already?
    do ikpt2=1,nFSkpt
     ss= (redkpt(1)-tmpFSkpt(1,ikpt2))**2 + &
&     (redkpt(2)-tmpFSkpt(2,ikpt2))**2 + &
&     (redkpt(3)-tmpFSkpt(3,ikpt2))**2
     if (ss < tol6) then
      new=0
      exit
     end if
    end do

    if (new == 1) then
     FSirredwtk(ikpt1)=FSirredwtk(ikpt1)+1
     nFSkpt=nFSkpt+1
     tmpFSkpt(:,nFSkpt) = redkpt(:)
     tmpFSfulltoirred(1,nFSkpt) = ikpt1
!    save sym that sends irred kpt ikpt1 onto full kpt
     tmpFSfulltoirred(2,nFSkpt) = isym
     tmpFSfulltoirred(3,nFSkpt) = itim
    end if

   end do !end loop over irred k points
  end do !end loop over timrev
 end do !end loop over symmetry

 write (*,*)'mkfskgrid: after first evaluation, nFSkpt = ', nFSkpt

 FSirredwtk(:) = FSirredwtk(:) / nFSkpt

!find correspondence table between irred FS kpoints and a full one
 do ikpt1=1,nFSkptirred
  do ikpt2=1,nFSkpt
   ss= (FSkptirred(1,ikpt1)-tmpFSkpt(1,ikpt2))**2 + &
&   (FSkptirred(2,ikpt1)-tmpFSkpt(2,ikpt2))**2 + &
&   (FSkptirred(3,ikpt1)-tmpFSkpt(3,ikpt2))**2
   if (ss < tol6) then
    FSirredtofull(ikpt1) = ikpt2
    exit 
   end if
  end do
 end do

!write (*,*) 'mkfskgrid: after determination of FSirredtofull '

 allocate (rankFSkpt(nFSkpt),invrankFSkpt(16000000))
 call mkkptrank (tmpFSkpt,nFSkpt,rankFSkpt,invrankFSkpt)

!DEBUG
!write (*,*) 'mkfskgrid: rankFSkpt = '
!write (*,'(6i14)') rankFSkpt
!write (*,*) 'mkfskgrid: invrankFSkpt = '
!write (*,'(6i14)') invrankFSkpt
!ENDDEBUG

!find correspondence table between FS kpoints

 tmpFSfulltofull(:,:,:) = -999
 do ikpt1=1,nFSkpt
! generate symmetrics of kpt ikpt1
  do isym=1,nsym
   do itim=0,timrev
    timsign = one-two*itim
    kpt(:) = timsign*(symrec(:,1,isym)*tmpFSkpt(1,ikpt1) + &
&    symrec(:,2,isym)*tmpFSkpt(2,ikpt1) + &
&    symrec(:,3,isym)*tmpFSkpt(3,ikpt1))
    call canon9(kpt(1),redkpt(1),res)
    call canon9(kpt(2),redkpt(2),res)
    call canon9(kpt(3),redkpt(3),res)

    new=1
!   which kpt is it among the full FS kpts
    symrankFSkpt = int(8000000.0_dp*(redkpt(1)+half+tol8) + &
&    40000.0_dp*(redkpt(2)+half+tol8) + &
&    200.0_dp*(redkpt(3)+half+tol8))
    if (symrankFSkpt > 16000000) then
     write (*,*) ' mkfskgrid : error : rank should be inferior to ', 2000000
     stop
    end if
    ikpt2 = invrankFSkpt(symrankFSkpt)
    if (ikpt2 /= -1) then
     tmpFSfulltofull(itim+1,isym,ikpt2) = ikpt1
     new = 0
    end if

!   ============= OLD VERSION ==================
!   ! which kpt is it among the full FS kpts
!   do ikpt2=1,nFSkpt
!   ss=(redkpt(1)-tmpFSkpt(1,ikpt2))**2 + &
!   & (redkpt(2)-tmpFSkpt(2,ikpt2))**2 + &
!   & (redkpt(3)-tmpFSkpt(3,ikpt2))**2
!   if (ss < tol6) then
!   new=0
!   !DEBUG
!   !if (tmpFSfulltofull(itim+1,isym,ikpt2) /= -999) then
!   !    write (*,*) 'Error : refilling element',tmpFSfulltofull(itim+1,isym,ikpt2)
!   !    write (121,*) 'Error : refilling element',tmpFSfulltofull(itim+1,isym,ikpt2)
!   !end if
!   !ENDDEBUG
!   
!   !tmpFSfulltofull(itim+1,isym,ikpt1) = ikpt2
!   tmpFSfulltofull(itim+1,isym,ikpt2) = ikpt1
!   
!   !DEBUG
!   !                write (121,'(4I5)') itim+1,isym,ikpt2,ikpt1
!   !                write (121,'(3(E16.6,2x))') tmpFSkpt(:,ikpt1),tmpFSkpt(:,ikpt2)
!   !                write (121,'(3(3(I3),2x))') symrec(:,:,isym)
!   !ENDDEBUG
!   exit
!   end if
!   end do
    if (new == 1) then
     write (*,*) 'mkfskgrid Error: FS kpt ',ikpt1,&
&     ' has no symmetric under sym', isym,&
&     ' with itim ',itim
     write (*,*) ' redkpt = ', redkpt
     write (*,*) 'symrankFSkpt,ikpt2 = ', symrankFSkpt,ikpt2
     stop
    end if
   end do
  end do
 end do

 deallocate (rankFSkpt,invrankFSkpt)

 write (*,*) 'mkfskgrid: after determination of FSfulltofull '

!DEBUG
!write (123,*) 'mkFSkgrid : FSfulltofull = '
!do ikpt1=1,nFSkpt
!write (123,*) 'ikpt1 = ',ikpt1
!write (123,'(4(2i5,4x))') (tmpFSfulltofull(:,isym,ikpt1),isym=1,nsym)
!end do
!ENDDEBUG


!got nFSkpt, tmpFSkpt, tmpFSfulltoirred, tmpFSfulltofull, and FSirredwtk

end subroutine mkFSkgrid
!!***
