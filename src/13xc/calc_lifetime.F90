!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_lifetime
!! NAME
!! calc_lifetime
!!
!! FUNCTION
!! Calculate the positron lifetime
!!
!! NOTE
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      invcb,lifetime_bn,lifetime_rpa
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_lifetime(ixcpositron,lifetime,nfft,nspden,positron,rhocore,rhore,rhototp,ucvol)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13xc, except_this_one => calc_lifetime
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixcpositron,nfft,nspden,positron
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: lifetime
!arrays
 real(dp),intent(in) :: rhocore(nfft),rhore(nfft,nspden),rhototp(nfft)

!Local variables-------------------------------
!scalars
 integer :: ifft
 real(dp),parameter :: rsfac=0.6203504908994000d0
!arrays
 real(dp),allocatable :: rhotote(:),rhovalre(:),rsepts(:),rseptstot(:)
 real(dp),allocatable :: rsppts(:)

! *************************************************************************

 allocate(rsepts(nfft),rseptstot(nfft),rsppts(nfft),rhotote(nfft),rhovalre(nfft))

 call invcb(rhototp,rsppts,nfft)
 rsppts(:)=rsfac*rsppts(:)


 if (positron==2) then
  do ifft=1,nfft
   rhotote(ifft)=rhore(ifft,1)+rhocore(ifft)
  end do
  if(nspden==2) then
   do ifft=1,nfft
    rhotote(ifft)=rhotote(ifft)+rhore(ifft,2)+0.5d0*rhocore(ifft)
   end do
  end if
  call invcb(rhotote,rseptstot,nfft)
  rseptstot(:)=rsfac*rseptstot(:)
  call lifetime_bn(lifetime,nfft,rhotote,rseptstot,rhototp,rsppts,ucvol)
! call lifetime_rpa(lifetime,nfft,rhor,rhocore,rsepts,rhototp,ucvol)
 else
! call lifetime_bn(lifetime,nfft,rhotote,rsepts,rhototp,rsppts,ucvol)
  rhovalre(:)=rhore(:,1)
  call invcb(rhovalre,rsepts,nfft)
  rsepts(:)=rsfac*rsepts(:)
! allocate(rseptstot(nfft))
  do ifft=1,nfft
   rhotote(ifft)=rhore(ifft,1)+rhocore(ifft)
  end do
  if(nspden==2) then
   do ifft=1,nfft
    rhotote(ifft)=rhotote(ifft)+rhore(ifft,2)+0.5d0*rhocore(ifft)
   end do
  end if
  call invcb(rhotote,rseptstot,nfft)
  rseptstot(:)=rsfac*rseptstot(:)
! call lifetime_bn(lifetime,nfft,rhotote,rseptstot,rhototp,rsppts,ucvol)
  call lifetime_rpa(lifetime,nfft,rhovalre,rhocore,rsepts,rseptstot,rhototp,ucvol)
 end if
 deallocate(rsepts,rseptstot,rsppts,rhotote,rhovalre)
end subroutine calc_lifetime
!!***
