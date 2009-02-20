!{\src2tex{textfont=tt}}
!!****f* ABINIT/lifetime_rpa
!! NAME
!! lifetime_rpa
!!
!! FUNCTION
!! Compute positron lifetime in the case of a single positron in a
!! homogeneous electrons gas for the zero positron density limit.
!! Returns lifetime.
!!
!! NOTE
!! Boronski and Nieminen parameterization of Arponen and Pajanne
!! of the contact density g0(rse).
!!
!! J. Arponen and E. Pajanne, Ann. Phys. (N.Y.) 121, 343 (1979).
!! E. Boronski and R.M. Nieminen, Phys. Rev. B 34, 3820 (1986).
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npt=number of real space points on which density is provided
!!  stepint=radial step used for radial integrations
!!  rhoer(npt)=electron number density (bohr^-3)
!!  rsepts(npt)=corresponding Wigner-Seitz radii, precomputed
!!  rhopr(npt)=positron number density (bohr^-3)
!!
!! OUTPUT
!!  lifetime=1/lambda (picoseconds).
!!  lambda=Pi*r0^2*c*Int(rhoe(r)*rhop(r)*g0(0;rhoe;rhop)*dr)
!!
!! PARENTS
!!      calc_lifetime
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine lifetime_rpa(lifetimeEq1,npt,rhoer,xccc3d,rsepts,rseptstot,rhopr,ucvol)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: npt
 real(dp) :: lifetimeEq1,ucvol
!arrays
 real(dp) :: rhoer(npt),rhopr(npt),rsepts(npt),rseptstot(npt),xccc3d(npt)

!Local variables-------------------------------
!cc is light celerity in angstrom/second
!r0sq is Bohr radius squared
!scalars
 integer :: ipt
 real(dp),parameter :: cc=3.0_dp*10.0_dp**18/0.529177
 real(dp),parameter :: r0sq=2.817940285*10.0_dp**(-5)/0.529177
 real(dp) :: g0,g0tot,lifetimeEq5,nbe,nbp,nbval,ratecoreEq1,ratecoreEq5
 real(dp) :: ratevalEq1,ratevalEq5,rhoe,rhoecore,rhop,rse,rsetot
 character(len=500) :: message
!arrays
 real(dp) :: ff(npt)

! *************************************************************************

 lifetimeEq1 = 0._dp
 lifetimeEq5 = 0._dp
 ratevalEq1 = 0._dp
 ratevalEq5 = 0._dp
 ratecoreEq1 = 0._dp
 ratecoreEq5 = 0._dp
 nbe =  0._dp
 nbp = 0._dp
 nbval = 0._dp
!Loop over grid points
 do ipt=1,npt
  rse=rsepts(ipt)
  rsetot=rseptstot(ipt)
  rhoe=rhoer(ipt)
  rhop=rhopr(ipt)
  rhoecore=xccc3d(ipt)

  g0=1.0_dp+1.23_dp*rse+0.8295*rse**(3.0_dp/2.0_dp)-1.26_dp*rse**2.0_dp+&
&  0.3286_dp*rse**(5.0_dp/2.0_dp)+1.0_dp/6.0_dp*rse**3._dp
  g0tot=1.0_dp+1.23_dp*rsetot+0.8295*rsetot**(3.0_dp/2.0_dp)-1.26_dp*rsetot**2.0_dp+&
&  0.3286_dp*rsetot**(5.0_dp/2.0_dp)+1.0_dp/6.0_dp*rsetot**3._dp
! modif to use the same formula as Sterne and Kaiser
! g0=1.0_dp+0.1512_dp*rse+2.414_dp*rse**(3.0_dp/2.0_dp)-2.01_dp*rse**2._dp+&
! & 0.4466*rse**(2.5_dp)+0.1667*rse**(3.0_dp)
  ff(ipt)= rhop*(rhoe*g0+rhoecore)
  nbval = nbval + rhoe
  nbe = nbe + rhoecore
  nbp = nbp + rhop
  lifetimeEq5 = lifetimeEq5 + ff(ipt)
  lifetimeEq1 = lifetimeEq1 + rhop*(rhoe+rhoecore)*g0tot
  ratevalEq1 = ratevalEq1 + rhop*rhoe*g0tot
  ratevalEq5 = ratevalEq5 + rhop*rhoe*g0
  ratecoreEq1 = ratecoreEq1 + rhop*rhoecore*g0tot
  ratecoreEq5 = ratecoreEq5 + rhop*rhoecore
 end do

 lifetimeEq1 = lifetimeEq1 * ucvol / npt
 lifetimeEq5 = lifetimeEq5 * ucvol / npt
 lifetimeEq1 = 1.0_dp/(pi*r0sq*r0sq*cc*lifetimeEq1)*1.e12_dp
 lifetimeEq5 = 1.0_dp/(pi*r0sq*r0sq*cc*lifetimeEq5)*1.e12_dp

 ratevalEq1 = ratevalEq1 * ucvol / npt * pi*r0sq*r0sq*cc * 1.e-9_dp
 ratevalEq5 = ratevalEq5 * ucvol / npt * pi*r0sq*r0sq*cc * 1.e-9_dp
 ratecoreEq1 = ratecoreEq1 * ucvol / npt * pi*r0sq*r0sq*cc * 1.e-9_dp
 ratecoreEq5 = ratecoreEq5  * ucvol / npt * pi*r0sq*r0sq*cc * 1.e-9_dp


 nbe = nbe * ucvol / npt
 nbp = nbp * ucvol / npt
 nbval = nbval * ucvol / npt
!lifetime = lifetime + lifetime_core
!lifetime = 1.0_dp/(pi*r0sq*r0sq*cc*lifetime)
 write(message, '(a,a,es16.8,a,a,es16.8,a,a,es16.8,a,a,es16.8,a,a,es16.8,a,a,es16.8)' ) ch10,&
& ' Lifetime EQ 1 (in psec)       =',lifetimeEq1,ch10,&
& ' Lifetime EQ 5 (in psec)       =',lifetimeEq5,ch10,&
& ' Annihilation rate valence EQ1 =',ratevalEq1,ch10,&
& ' Annihilation rate core    EQ1 =',ratecoreEq1,ch10,&
& ' Annihilation rate valence EQ5 =',ratevalEq5,ch10,&
& ' Annihilation rate core    EQ5 =',ratecoreEq5
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')
 write(message, '(a,es16.8,a,a,es16.8,a,a,es16.8,a,a,i6)' ) &
& ' Number of core electrons      =',nbe,ch10,&
& ' Number of valence electrons   =',nbval,ch10,&
& ' Number of positrons           =',nbp,ch10,&
& ' Number of real space points on which density is provided = ',npt
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

end subroutine lifetime_rpa
!!***
