!{\src2tex{textfont=tt}}
!!****f* ABINIT/findnq
!! NAME
!! findnq
!!
!! FUNCTION
!! Identify the number of q-points in the IBZ by which the k-points in BZ differ
!! (count the q points in the k-point difference set)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MT, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kbz(3,nkbz)=coordinates of k points in BZ
!!  timrev=2 if time-reversal symmetry is used, 1 otherwise
!!  nkbz=number of k points in Brillouin zone
!!  nsym=number of symmetry operations
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!
!! OUTPUT
!!  nqibz=number of q points
!!
!! PARENTS
!!      mrgscr,screening
!!
!! CHILDREN
!!      dosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine findnq(nkbz,kbz,nsym,symrec,nqibz,timrev)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => findnq
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: timrev,nkbz,nsym
 integer,intent(out) :: nqibz
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: kbz(3,nkbz)

!Local variables ------------------------------
!scalars
 integer :: ifound,ik,isym,iq,istat,memory_exhausted,nqall,nqallm,itim
 character(len=500) :: msg
!arrays
 integer :: gtemp(3),g0(3)
 real(dp) :: qposs(3),qrot(3)
 real(dp),allocatable :: qall(:,:)
!************************************************************************

 ! === Infinite do-loop to be able to allocate sufficient memory ===
 nqallm=1000 ; memory_exhausted=0
 do 
  allocate(qall(3,nqallm),stat=istat)
  if (istat/=0) call memerr('findnq','qall',3*nqallm,'dp')
  nqall=0
  ! === Loop over all k-points in BZ, forming k-k1 ===
  do ik=1,nkbz
   qposs(:)=kbz(:,ik)-kbz(:,1)
   ! === Check whether this q (or its equivalent) has already been found ===
   ifound=0
   do iq=1,nqall
    do itim=1,timrev
     do isym=1,nsym
      !FIXME this is for g95
      call dosym(REAL(symrec(:,:,isym),dp),itim,qall(:,iq),qrot)
      if (is_samek(qrot,qposs,g0)) ifound=ifound+1
     end do
    end do
   end do

   if (ifound==0) then
    nqall=nqall+1
    !
    ! === If not yet found, check that the allocation is big enough ===
    if (nqall>nqallm) then
     memory_exhausted=1 ; deallocate(qall)
     nqallm=nqallm*2    ; EXIT ! Exit the do ik=1 loop
    end if
    ! === Add to the list ===
    qall(:,nqall)=qposs(:)
   end if
  end do

  if (memory_exhausted==0) EXIT
 end do !infinite loop

 deallocate(qall)
 nqibz=nqall

 write(msg,'(a,i8)')' number of q-points found ',nqibz
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

end subroutine findnq
!!***
