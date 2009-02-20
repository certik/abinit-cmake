!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcppsn_tcdft
!! NAME
!! xcppsn_tcdft
!!
!! FUNCTION
!! Compute electron-positron correlation potential for positron.
!! Returns vxc from input rhopr and rhoer for
!! positron and electrons.
!!
!! NOTE
!! Puska, Seitsonen and Nieminen parameterization of Arponen and Pajanne
!! electron-positron  gas energy data--
!! J. Arponen and E. Pajanne, Ann. Phys. (N.Y.) 121, 343 (1979).
!! E. Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986).
!! M.J. Puska, A.P. Seitsonen, and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  rhoer(npt)=electron number density (bohr^-3)
!!  rhopr(npt)=positron number density (bohr^-3)
!!  rsepts(npt)=corresponding Wigner-Seitz radii, precomputed
!!  rsppts(npt)=corresponding Wigner-Seitz radii, precomputed
!!
!! OUTPUT
!!  vxc(npt)=xc potential (d($\rho$*exc)/d($\rho$)) (hartree)
!!
!! PARENTS
!!      calc_xc_ep
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xcppsn_tcdft(npt,rhoer,rsepts,rhopr,rsppts,vxc)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt
!arrays
 real(dp),intent(in) :: rhoer(npt),rhopr(npt),rsepts(npt),rsppts(npt)
 real(dp),intent(out) :: vxc(npt)

!Local variables-------------------------------
!Perdew-Zunger parameters a, b, b1, b2, c, d, gamma
!scalars
 integer :: ipt
 real(dp),parameter :: aa=69.7029_dp,apa=-0.6298_dp,apb=-13.15111_dp
 real(dp),parameter :: apc=2.8655_dp,ba=-107.4927_dp,bb=141.8458_dp
 real(dp),parameter :: ca=23.7182_dp,cb=-33.6472_dp,cc=5.21152_dp
 real(dp) :: arse,brse,crse,deape,deapp,eape,eapp,ec,fac,invrhoe,invrhop,invrse
 real(dp) :: invrsp,rhoe,rhop,rse,rsesq,rsp

! *************************************************************************

!Compute fac=4*Pi/3
 fac=4.0_dp*pi/3.0_dp

!Loop over grid points
 do ipt=1,npt
  rse=rsepts(ipt)
  rsp=rsppts(ipt)
  invrse=1.0_dp/rse
  invrsp=1.0_dp/rsp
  rhop=rhopr(ipt)
  invrhop=1.0_dp/rhop
  rhoe=rhoer(ipt)
  invrhoe=1.0_dp/rhoe

  rsesq=rse**2.0_dp
  arse=aa+ba*rse+ca*rsesq
  brse=ba+bb*rse+cb*rsesq
  crse=ca+cb*rse+cc*rsesq

  if (rse<0.302_dp) then
   eape=0.5_dp*(-1.56_dp/rse**0.5_dp+(0.051_dp*log(rse)-&
&   0.081_dp)*log(rse)+1.14_dp)
   deape=0.5_dp*(-invrhoe*rse/3._dp)&
&   *(-0.5_dp*1.56*rse**(-0.5_dp)+1/rse*(0.051*log(rse)-0.081)&
&   +0.051/rse*log(rse))
  else if (rse>=0.302_dp .and. rse<=0.56_dp) then
   eape=0.5_dp*(-0.92305_dp-0.05459_dp/rse**2.0_dp)
   deape=0.5_dp*(-invrhoe*rse/3._dp)*2._dp*0.05459/rse**3
  else if (rse>0.56_dp .and. rse<=8.0_dp) then
   eape=0.5_dp*(apa+apb/(rse+2.5)**2.0_dp+apc/(rse+2.5))
   deape=0.5_dp*(-invrhoe*rse/3._dp)*&
&   (2._dp*apa/(rse+2.5)**3.0_dp-apb/(rse+2.5)**2.0_dp)
  else
   eape=0.5_dp*(-179856.2768_dp*rhoe**2.0_dp+186.4207_dp*rhoe-0.524_dp)
   deape=0.5_dp*(-2._dp*179856.2768_dp*rhoe+186.4207_dp)
  end if

  if (rsp<0.302_dp) then
   eapp=0.5_dp*(-1.56_dp/rsp**0.5_dp+(0.051_dp*log(rsp)-&
&   0.081_dp)*log(rsp)+1.14_dp)
   deapp=0.5_dp*(-invrhop*rsp/3._dp)&
&   *(-0.5_dp*1.56*rsp**(-0.5_dp)+1/rsp*(0.051*log(rsp)-0.081)&
&   +0.051/rsp*log(rsp))
  else if (rsp>=0.302_dp .and. rsp<=0.56_dp) then
   eapp=0.5_dp*(-0.92305_dp-0.05459_dp/rsp**2.0_dp)
   deapp=0.5_dp*(-invrhop*rsp/3._dp)*2._dp*0.05459/rsp**3
  else if (rsp>0.56_dp .and. rsp<=8.0_dp) then
   eapp=0.5_dp*(apa+apb/(rsp+2.5)**2.0_dp+apc/(rsp+2.5))
   deapp=0.5_dp*(-invrhop*rsp/3._dp)*&
&   (2._dp*apa/(rsp+2.5)**3.0_dp-apb/(rsp+2.5)**2.0_dp)
  else
   eapp=0.5_dp*(-179856.2768_dp*rhop**2.0_dp+186.4207_dp*rhop-0.524_dp)
   deapp=0.5_dp*(-2._dp*179856.2768_dp*rhop+186.4207_dp)
  end if

! eape=apa+apb/(rse+2.5)**2.0_dp+apc/(rse+2.5)
! eapp=apa+apb/(rsp+2.5)**2.0_dp+apc/(rsp+2.5)
  ec=1.0_dp/(arse+brse*rsp+crse*rsp**2.0_dp+fac*rsp**3.0_dp/eape+&
&  fac*rse**3.0_dp/eapp)
! vxc given in rydberg units
  vxc(ipt)=-ec**2.0_dp*(-brse*rsp*invrhop/3.0_dp &
&  -crse*2.0_dp/3.0_dp*rsp**2.0_dp*invrhop &
&  -invrhop**2/eape &
&  -invrhoe/eapp**2*deapp&
&  )
! conversion in hartree units
  vxc(ipt)=0.5_dp*vxc(ipt)
 end do
!
end subroutine xcppsn_tcdft
!!***
