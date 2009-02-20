!{\src2tex{textfont=tt}}
!!****f* ABINIT/findqg0
!! NAME
!! findqg0
!!
!! FUNCTION
!! Identify q + g0 = k - kp
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kmkp(3)= k - kp input vector
!!  nqbz=number of q points
!!  qbz(3,nqbz)=coordinates of q points
!!  mG0(3)= For each reduced direction gives the max G0 component to account for umklapp processes
!!
!! OUTPUT
!!  iq=index of q vector in BZ table
!!  g0(3)=reciprocal space vector, to be used in igfft
!!
!! PARENTS
!!      csigme
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine findqg0(iq,g0,kmkp,nqbz,qbz,mG0)

 use defs_basis
 use m_gwdefs, only : GW_TOLQ


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqbz
 integer,intent(out) :: iq
!arrays
 integer,intent(in) :: mG0(3)
 integer,intent(out) :: g0(3)
 real(dp),intent(in) :: kmkp(3),qbz(3,nqbz)

!Local variables-------------------------------
!FIXME if I use 1.0d-4 the jobs crash, should understand why
!scalars
 integer :: iqbz,jg01,jg02,jg03,nq0found
 real(dp) :: tolq0=1.0D-3
 character(len=500) :: msg
!arrays
 real(dp) :: kp(3),qpg0(3)

! *************************************************************************

 iq=0
 if (ALL(ABS(kmkp(:))<epsilon(one))) then
  ! === Find q close to 0 ===
  nq0found=0
  do iqbz=1,nqbz
   if (ALL(ABS(qbz(:,iqbz))<tolq0)) then 
    iq=iqbz
    nq0found=nq0found+1
   end if 
  end do
  if (iq==0) then 
   write(msg,'(5a)')ch10,&
&  ' findqg0 BUG- ',ch10,&
   ' wrong SCR file: q=0 not existing',ch10
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if
  if (nq0found/=1) then 
   write(msg,'(6a)')ch10,&
&   ' findqg0 BUG- ',ch10,&
&   ' It seems that the SCR file contains more than one small q-point',ch10,&
&   ' Try to reduce the value of tolq0'
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if 
  g0(:)= 0 ; RETURN
 else
  ! === q is not zero, find q such as k-kp=q+G0 ===
  do iqbz=1,nqbz
   do jg01=-mG0(1),mG0(1)
    do jg02=-mG0(2),mG0(2)
     do jg03=-mG0(3),mG0(3)
      ! === Form q-G0 and check if it is the one ===
      qpg0(1)= qbz(1,iqbz)+jg01
      qpg0(2)= qbz(2,iqbz)+jg02
      qpg0(3)= qbz(3,iqbz)+jg03
      if (ALL(ABS(qpg0(:)-kmkp(:))<GW_TOLQ)) then
       if (iq/=0) then
        write(msg,'(4a)')ch10,&
&        ' findqg0 : ERROR - ',ch10,&
&        ' duplicated q+g0.'
        call wrtout(std_out,msg,'COLL')
        write(std_out,*)iqbz,qbz(:,iqbz),jg01,jg02,jg03
        write(std_out,*)iq,qbz(:,iq),g0
        call leave_new('COLL')
       end if
       iq=iqbz
       g0(1)=jg01
       g0(2)=jg02
       g0(3)=jg03
      end if
     end do
    end do
   end do
  end do
  if (iq==0) then
   write(msg,'(4a)')ch10,&
&   ' findqg0 : BUG - ',ch10,&
&   ' q = k - kp + g0 not found.'
   call wrtout(std_out,msg,'COLL') ; write(*,*)' kmkp  = ',kmkp(:) ; call leave_new('COLL')
  end if
 end if

end subroutine findqg0
!!***
