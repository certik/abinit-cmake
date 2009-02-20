!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcepsn_tcdft
!! NAME
!! xcepsn_tcdft
!!
!! FUNCTION
!! Compute electron-positron correlation potential for electron.
!! Returns vxc from input rhopr and rhoer for positron and electrons.
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

subroutine xcepsn_tcdft(npt,rhoer,rsepts,rhopr,rsppts,vxc)

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
  invrse=one/rse
  invrsp=one/rsp
  rhop=rhopr(ipt)
  invrhop=one/rhop
  rhoe=rhoer(ipt)
  invrhoe=one/rhoe


  rsesq=rse**2.0_dp
  arse=aa+ba*rse+ca*rsesq
  brse=ba+bb*rse+cb*rsesq
  crse=ca+cb*rse+cc*rsesq

  if (rse<0.302_dp) then
   eape=half*(-1.56_dp/rse**half+(0.051_dp*log(rse)-&
&   0.081_dp)*log(rse)+1.14_dp)
   deape=half*(-invrhoe*rse/3._dp)&
&   *(-half*1.56*rse**(-half)+1/rse*(0.051*log(rse)-0.081)&
&   +0.051/rse*log(rse))
  else if (rse>=0.302_dp .and. rse<=0.56_dp) then
   eape=half*(-0.92305_dp-0.05459_dp/rse**2.0_dp)
   deape=half*(-invrhoe*rse/3._dp)*2._dp*0.05459/rse**3
  else if (rse>0.56_dp .and. rse<=8.0_dp) then
   eape=half*(apa+apb/(rse+2.5)**2.0_dp+apc/(rse+2.5))
   deape=half*(-invrhoe*rse/3._dp)*&
&   (2._dp*apa/(rse+2.5)**3.0_dp-apb/(rse+2.5)**2.0_dp)
  else
   eape=half*(-179856.2768_dp*rhoe**2.0_dp+186.4207_dp*rhoe-0.524_dp)
   deape=half*(-2._dp*179856.2768_dp*rhoe+186.4207_dp)
  end if

  if (rsp<0.302_dp) then
   eapp=half*(-1.56_dp/rsp**half+(0.051_dp*log(rsp)-&
&   0.081_dp)*log(rsp)+1.14_dp)
   deapp=half*(-invrhop*rsp/3._dp)&
&   *(-half*1.56*rsp**(-half)+1/rsp*(0.051*log(rsp)-0.081)&
&   +0.051/rsp*log(rsp))
  else if (rsp>=0.302_dp .and. rsp<=0.56_dp) then
   eapp=half*(-0.92305_dp-0.05459_dp/rsp**2.0_dp)
   deapp=half*(-invrhop*rsp/3._dp)*2._dp*0.05459/rsp**3
  else if (rsp>0.56_dp .and. rsp<=8.0_dp) then
   eapp=half*(apa+apb/(rsp+2.5)**2.0_dp+apc/(rsp+2.5))
   deapp=half*(-invrhop*rsp/3._dp)*&
&   (2._dp*apa/(rsp+2.5)**3.0_dp-apb/(rsp+2.5)**2.0_dp)
  else
   eapp=half*(-179856.2768_dp*rhop**2.0_dp+186.4207_dp*rhop-0.524_dp)
   deapp=half*(-2._dp*179856.2768_dp*rhop+186.4207_dp)
  end if

! eape=apa+apb/(rse+2.5)**2.0_dp+apc/(rse+2.5)
! eapp=apa+apb/(rsp+2.5)**2.0_dp+apc/(rsp+2.5)
  ec=one/(arse+brse*rsp+crse*rsp**2.0_dp+fac*rsp**3.0_dp/eape+&
&  fac*rse**3.0_dp/eapp)
! vxc given in rydberg units
  vxc(ipt)=-ec**2.0_dp*(-ba*rse*invrhoe/3._dp-ca*2._dp*rse**2*invrhoe/3._dp &
&  -rsp*(bb*rse*invrhoe/3._dp+cb*2._dp*rse**2*invrhoe/3._dp) &
&  -rsp**2*(cb*rse*invrhoe/3._dp+cc*2._dp*rse**2*invrhoe/3._dp) &
&  -invrhop/eape**2*deape &
&  -invrhoe**2/eapp &
&  )
! conversion in hartree units
  vxc(ipt)=half*vxc(ipt)
 end do
!
end subroutine xcepsn_tcdft
!!***
