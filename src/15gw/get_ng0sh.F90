!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_ng0sh
!! NAME
!! get_ng0sh
!!
!! FUNCTION
!!  This routine calculates the difference k_BZ-q_IBZ for each k point in the full BZ and each q point in the IBZ.
!!  For each difference it finds the umklapp G0 vector required to restore k-q in the first Brillouin zone.
!!  The optimal value of G0 shells is returned, namely the smallest box around Gamma that is sufficient to 
!!  treat all the possible umklapp processes.
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gmet(3,3)=reciprocal space metric ($\textrm{bohr}^{-2}$).
!!  nfold=number of points in the array kfold
!!  nk1, nk2=number of points in the arrays kbz1, kbz2
!!  kbz1(3,nk1)=reduced coordinates of the first set of points
!!  kbz2(3,nk2)=reduced coordinates of the second set of points
!!  kfold(3,nkfol)=reduced coordinated of the points in the BZ
!!  tolq0=tolerance below which a q-point is treated as zero
!!  mg0sh=Integer defining the Max number of shells to be considered
!!
!! OUTPUT
!!  opt_ng0(3)=Minimum reduced compoments of the G0 vectors that must to be considered 
!!   such as k-q+G0 belong to the first BZ for each considered k- and q-point. 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  
!!
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_ng0sh(nk1,kbz1,nk2,kbz2,nkfold,kfold,gmet,tolq0,mg0sh,opt_ng0)

 use defs_basis
 use m_gwdefs, only : GW_TOLQ
 use m_numeric_tools, only : is_zero


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mg0sh,nk1,nk2,nkfold
 real(dp),intent(in) :: tolq0
!arrays
 integer,intent(out) :: opt_ng0(3)
 real(dp),intent(in) :: gmet(3,3),kbz1(3,nk1),kbz2(3,nk2),kfold(3,nkfold)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ig,ig1,ig2,ig3,ii,ikf,itr,ng0
 logical :: found
 character(len=500) :: msg
!arrays
 integer :: gtrial(3)
 integer,allocatable :: iperm(:)
 real(dp) :: k1mk2(3)
 real(dp),allocatable :: norm(:)
!no_abirules
 real(dp),allocatable :: g0(:,:) !TODO should be integer but I have to overload

!************************************************************************
 
 ng0=(2*mg0sh+1)**3
 allocate(g0(3,ng0),norm(ng0),iperm(ng0))

 ii=1
 do ig3=-mg0sh,mg0sh
  do ig2=-mg0sh,mg0sh
   do ig1=-mg0sh,mg0sh
    g0(1,ii)=ig1
    g0(2,ii)=ig2
    g0(3,ii)=ig3
    !TODO overload normv
    norm(ii)=normv(g0(:,ii),gmet,'g')
    iperm(ii)=ii
    ii=ii+1
   end do
  end do
 end do
 !
 ! === Order g0-vectors in order of increasing module ===
 ! * Should speed up the search in the while statement below
 call sort_dp(ng0,norm,iperm,tol14)

 opt_ng0(:)=0 
 do i2=1,nk2
  ! This is used in case of screening calculation
  ! If q is small treat it as zero. In this case, indeed, 
  ! we use q=0 to calculate the oscillator matrix elements. 
  if (is_zero(kbz2(:,i2),tolq0)) CYCLE
  do i1=1,nk1 
   k1mk2(:)=kbz1(:,i1)-kbz2(:,i2)
   itr=1 ; found=.FALSE.
   do while (itr<=ng0.and..not.found)
    do ikf=1,nkfold
     if (ALL(ABS(k1mk2(:)-kfold(:,ikf)-g0(:,iperm(itr))) < GW_TOLQ)) then 
      found=.TRUE. 
      gtrial(:)=ABS(g0(:,iperm(itr)))
      opt_ng0(1)=MAX(opt_ng0(1),gtrial(1))
      opt_ng0(2)=MAX(opt_ng0(2),gtrial(2))
      opt_ng0(3)=MAX(opt_ng0(3),gtrial(3))
      EXIT
     end if
    end do
    itr=itr+1
   end do
   if (.not.found) then 
    write(msg,'(4a,i5,3a,i2,2(2a,i4,3es16.8),a)')ch10,&
&    ' get_ng0sh : ERROR -',ch10,&
&    '  Not able to found umklapp G0 vector among ',ng0,' vectors',ch10,&
&    '  increase mg0sh such as k1-k2 = kf+G0, present value = ',mg0sh,ch10,&
&    '  point1 = ',i1,kbz1(:,i1),ch10,&
&    '  point2 = ',i2,kbz2(:,i2),ch10
    call wrtout(std_out,msg,'COLL') ;  write(*,*)kfold
    call leave_new('COLL')
   end if
  end do
 end do 

 write(msg,'(a,3i2)')' optimal value for ng0sh = ',opt_ng0(:)
 call wrtout(std_out,msg,'COLL')

 deallocate(g0,norm,iperm)

#if defined DEBUG_MODE
! do i1=1,nk1
!  do i2=1,nk2
!   k1mk2(:)=kbz1(:,i1)-kbz2(:,i2)
!   if (.not.findkp(nkfold,kfold,ikf,k1mk2,opt_ng0,gtrial)) then
!    write(msg,'(4a)')ch10,' get_ng0sh : BUG -',ch10,' kf=k1-k2-G0 not found'
!    call wrtout(std_out,msg,'COLL')
!    write(msg,'(2(2a,i4,a,3f12.6),2a,3f12.6)')ch10,&
!&    ' i1 = ',i1,' k1 = ',(kbz1(ii,i1),ii=1,3),ch10,&
!&    ' i2 = ',i2,' k2 = ',(kbz2(ii,i2),ii=1,3),ch10,&
!&    ' k1-k2 = ',(k1mk2(ii),ii=1,3)
!    call wrtout(std_out,msg,'COLL') ; call leave_new('COLL') ! the weight 1.0/nkbz
!   end if 
!  end do
! end do
! write(msg,'(a)')' get_ng0sh: exit'
! call wrtout(std_out,msg,'PERS')
! call flush_unit(std_out)
#endif

end subroutine get_ng0sh
!!***
