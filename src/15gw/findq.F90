!{\src2tex{textfont=tt}}
!!****f* ABINIT/findq
!! NAME
!! findq
!!
!! FUNCTION
!! Identify the q-points by which the k-points in BZ differ
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MT, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  kbz(3,nkbz)=coordinates of k points in BZ
!!  timrev=2 if time-reversal symmetry is used, 1 otherwise
!!  nkbz=number of k points in Brillouin zone
!!  nsym=number of symmetry operations
!!  nqibz=number of q points in the IBZ by which k points differ (computed in findnq)
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!
!! OUTPUT
!!  qibz(3,nqibz)=coordinates of q points by which k points differ
!!
!! PARENTS
!!      mrgscr,rdm,setup_qmesh
!!
!! CHILDREN
!!      bz1,dosym,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine findq(nkbz,kbz,nsym,symrec,gprimd,nqibz,qibz,timrev,avoid_zero)

 use defs_basis
 use m_numeric_tools, only : is_zero


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => findq
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkbz,nqibz,nsym,timrev
 logical,intent(in) :: avoid_zero
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),kbz(3,nkbz)
 real(dp),intent(out) :: qibz(3,nqibz)

!Local variables ------------------------------
!scalars
 integer :: ii,ik,iq,iqp,isym,itim,jj,nq0
 real(dp) :: shift,tolq0
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: g0(3),gtemp(3)
 real(dp) :: gmet(3,3),qposs(3),qrot(3)

!************************************************************************

 write(msg,'(a)')' find q-points q = k - k1 and translate in first BZ'
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 ! Compute reciprocal space metrics
 do ii=1,3
  gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&            gprimd(2,ii)*gprimd(2,:)+&
&            gprimd(3,ii)*gprimd(3,:)
 end do

 !b1=two_pi*gprimd(:,1)
 !b2=two_pi*gprimd(:,2)
 !b3=two_pi*gprimd(:,3)

 tolq0=0.001_dp !old behaviour
 !
 ! === Loop over all k-points in BZ, forming k-k1 ===
 ! iq is the no. of q-points found, zero at the beginning
 iq=0
 do ik=1,nkbz
  qposs(:)=kbz(:,ik)-kbz(:,1)
  ! === Check whether this q (or its equivalent) has already been found ===
  found=.FALSE.
  do iqp=1,iq
   do itim=1,timrev
    do isym=1,nsym
     !FIXME this is for g95
     call dosym(REAL(symrec(:,:,isym),dp),itim,qibz(:,iqp),qrot)
     if (is_samek(qrot,qposs,g0)) found=.TRUE.
    end do
   end do
  end do
  if (.not.found) then
   iq=iq+1
   if (iq>nqibz) then 
    write(msg,'(2a)')ch10,&
&    ' findq : BUG in findnq: iq>nqibz '
    call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
   end if 
   qibz(:,iq)=qposs(:)
  end if
 end do

 if (iq/=nqibz) then 
  write(msg,'(2a)')ch10,&
&  ' findq : BUG in findnq: iq/=nqibz '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 
 !
 ! Translate q-points (calculated for info) to 1st BZ
 !
!MG it seems that bz1 sometimes does not work, 
!in SiO2 for example I got 0.333333   -0.666667    0.000000
!it is better to use canon9 also because if an irred qpoint 
!lies outside the 1BZ then most probably we have to consider
!an umklapp G0 vector to reconstruct the full BZ. The 
!correct treatment of this case is not yet implemented yet, see csigme.F90
!Anyway I should check this new method because likely mrscr will complain
!The best idea consists in  writing a new subroutine canon10 
!wich reduce the q point in the interval [-1/2,1/2[ which is
!supposed to be the scope of bz1. Just to obtain the same results as 
!the automatic tests
!FIXME for the moment use old version, easy for debugging

 do iq=1,nqibz
  call bz1(qibz(:,iq),gtemp,gmet)
 end do

!DEBUG
!do iq=1,nqibz
!do ii=1,3
!call canon9(qibz(ii,iq),qibz(ii,iq),shift)
!end do
!end do 
!ENDDEBUG
 write(msg,'(a)')' q-points [reduced coordinates]'
 call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')

 nq0=0
 do jj=1,nqibz
  if (is_zero(qibz(:,jj),tolq0).and.avoid_zero) then
   qibz(1,jj)=0.000010
   qibz(2,jj)=0.000020
   qibz(3,jj)=0.000030
   nq0=nq0+1
  end if
  write(msg,'(3f12.6)') (qibz(ii,jj),ii=1,3)
  call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 end do
 
 !write(msg,'(a)')ch10
 !call wrtout(std_out,msg,'COLL') ; call wrtout(ab_out,msg,'COLL')
 write(std_out,*) ; write(ab_out,*)

 if (nq0/=1) then 
  write(msg,'(4a,i2,5a)')ch10,&
&  ' findq : ERROR ',ch10,&
&  ' Found ',nq0,' "small" qpoints ',ch10,&
&  ' Check the q-mesh and, if it is correct, decrease the tolerance value ',ch10,&
&  ' below which two points are considered equivalent. '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

end subroutine findq
!!***
