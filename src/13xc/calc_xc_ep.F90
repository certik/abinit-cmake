!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_xc_ep
!! NAME
!! calc_xc_ep
!!
!! FUNCTION
!! Calculate the electrons/positron correlation term
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
!!  ixcpositron=
!!  nfft=
!!  nspden=
!!  option=
!!  positron=
!!
!! OUTPUT
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      invcb,xce_ap,xcepsn_tcdft,xcp_ap,xcppsn_tcdft
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_xc_ep(excapn,ixcpositron,positron,nfft,nspden,option,rhocore,rhor,rhore,rhorp,rhotote,rhototp,rsepts,rsppts,&
&                      vhae,vhap,vpsp,vtrial,vxcapn)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13xc, except_this_one => calc_xc_ep
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixcpositron,nfft,nspden,option,positron
!arrays
 real(dp),intent(in) :: rhocore(nfft),rhor(nfft,nspden),rhore(nfft,nspden)
 real(dp),intent(in) :: rhorp(nfft,nspden),vhae(nfft,2),vhap(nfft,2),vpsp(nfft)
 real(dp),intent(inout) :: rhotote(nfft),rhototp(nfft),rsepts(nfft)
 real(dp),intent(inout) :: rsppts(nfft),vtrial(nfft,nspden),vxcapn(nfft)
 real(dp),intent(out) :: excapn(nfft)

!Local variables-------------------------------
!scalars
 integer :: ifft
 real(dp),parameter :: rsfac=0.6203504908994000d0
!arrays
 real(dp),allocatable :: rseptstot(:)

! *************************************************************************

 if (positron==1) then
  if(nspden==1)then
   rhototp(:)=rhor(:,1)
  else
   rhototp(:)=rhor(:,1)+rhor(:,2)
  end if

  call invcb(rhototp,rsppts,nfft)
  rsppts(:)=rsfac*rsppts(:)

  do ifft=1,nfft
   rhotote(ifft)=rhore(ifft,1)+rhocore(ifft)
  end do
  if(nspden==2) then
   do ifft=1,nfft
    rhotote(ifft)=rhotote(ifft)+rhore(ifft,2)
   end do
  end if

  call invcb(rhotote,rsepts,nfft)
  rsepts(:)=rsfac*rsepts(:)
  if (ixcpositron==1) then
   call xcp_ap(excapn,nfft,rhotote,rsepts,vxcapn)
  else if (ixcpositron==2) then
   call xcppsn_tcdft(nfft,rhotote,rsepts,rhototp,rsppts,vxcapn)
  else
   write(6,*) 'This type of electron-positron correlation is not taken into account !'
   stop
  end if

  if (option==0) vtrial(:,1)=-vpsp(:)-vhae(:,1)+vxcapn(:)

 elseif (positron==2) then

  if(nspden==1)then
   rhototp(:)=rhorp(:,1)
  else
   rhototp(:)=rhorp(:,1)+rhorp(:,2)
  end if

  call invcb(rhototp,rsppts,nfft)
  rsppts(:)=rsfac*rsppts(:)

  do ifft=1,nfft
   rhotote(ifft)=rhor(ifft,1)+rhocore(ifft)
  end do
  if(nspden==2) then
   do ifft=1,nfft
    rhotote(ifft)=rhotote(ifft)+rhor(ifft,2)+0.5d0*rhocore(ifft)
   end do
  end if
  allocate(rseptstot(nfft))
  call invcb(rhotote,rseptstot,nfft)
  rseptstot(:)=rsfac*rseptstot(:)

  if (ixcpositron==1) then
   call xce_ap(excapn,nfft,rhotote,rsepts,rhototp,vxcapn)
  else if (ixcpositron==2) then
   call xcepsn_tcdft(nfft,rhotote,rseptstot,rhototp,rsppts,vxcapn)
  else
   write(6,*) 'This type of electron-positron correlation is not taken into account !'
   stop
  end if

  if (option==0) vtrial(:,1)= vtrial(:,1)-vhap(:,1)+vxcapn(:)
  deallocate(rseptstot)
 end if
end subroutine calc_xc_ep
!!***
