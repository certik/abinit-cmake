!{\src2tex{textfont=tt}}
!!****f* ABINIT/xce_ap
!! NAME
!! xce_ap
!!
!! FUNCTION
!! Compute electron-positron correlation potential for electron, in the zero
!! positron density limit.
!! Returns exc and vxc from input rhopr and rhoer for positron and electrons.
!!
!! NOTE
!! Boronski and Nieminen parameterization of Arponen and Pajanne
!! electron-positron  gas energy data--
!! J. Arponen and E. Pajanne, Ann. Phys. (N.Y.) 121, 343 (1979).
!! E. Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986).
!! M.J. Puska, A.P. Seitsonen, and R.M. Nieminen, Phys. Rev. B 52, 10947 (1994).
!!
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
!!  rhoer(npt)=electron density (bohr^-3)
!!  rhopr(npt)=positron density (bohr^-3)
!!  rsepts(npt)=corresponding Wigner-Seitz radii, precomputed
!!
!! OUTPUT
!!  exc(npt)=exchange-correlation energy density (hartree)
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

subroutine xce_ap(exc,npt,rhoer,rsepts,rhopr,vxc)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt
!arrays
 real(dp),intent(in) :: rhoer(npt),rhopr(npt),rsepts(npt)
 real(dp),intent(out) :: exc(npt),vxc(npt)

!Local variables-------------------------------
!Perdew-Zunger parameters a, b, b1, b2, c, d, gamma
!scalars
 integer :: ipt
 real(dp),parameter :: apa=-0.6298_dp,apb=-13.15111_dp,apc=2.8655_dp
 real(dp) :: drse,invrse,rhoe,rhop,rse

! *************************************************************************

!Loop over grid points
 do ipt=1,npt
  rse=rsepts(ipt)
  invrse=one/rse
  rhoe=rhoer(ipt)
  rhop=rhopr(ipt)
! Consider four regimes: rse<0.302
! 0.302<=rse<=0.56
! 0.56<rse<=8.0
! rse>8.0

  drse = -one/3.0_dp*rse/rhoe
  if (rse<0.302_dp) then
!  compute energy density exc (hartree)
   exc(ipt)=-drse*half*(half*1.56_dp/rse**1.5_dp &
&   +2_dp*0.051*log(rse)/rse-0.081/rse)

!  compute potential vxc=d(rhop*exc)/d(rhop) (hartree)
   vxc(ipt)=rhop*exc(ipt)
  else if (rse>=0.302_dp .and. rse<=0.56_dp) then
!  compute energy density exc (hartree)
   exc(ipt)=drse*half*(two*0.05459_dp/rse**3)
!  compute potential vxc=d(rhop*exc)/d(rhop) (hartree)
   vxc(ipt)=rhop*exc(ipt)
  else if (rse>0.56_dp .and. rse<=8.0_dp) then
!  compute energy density exc (hartree)
   exc(ipt)=drse*half*(-2._dp*apb/(rse+2.5_dp)**3-apc/(rse+2.5_dp)**2)
!  compute potential vxc=d(rhop*exc)/d(rhop) (hartree)
   vxc(ipt)=rhop*exc(ipt)
  else
!  compute energy density exc (hartree)
   exc(ipt)=half*(-two*179856.2768_dp*rhoe+186.4207_dp)
!  compute potential vxc=d(rhop*exc)/d(rhop) (hartree)
   vxc(ipt)=rhop*exc(ipt)
  end if
 end do
!
end subroutine xce_ap
!!***
