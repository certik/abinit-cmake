!{\src2tex{textfont=tt}}
!!****f* ABINIT/lifetime_psn
!! NAME
!! lifetime_psn
!!
!! FUNCTION
!! Compute positron lifetime in the case of a single positron in a
!! homogeneous electrons gas.
!! Returns lifetime.
!!
!! NOTE
!! Puska Seitsonen and Nieminen parameterization
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
!!  stepint=radial step used for radial integrations
!!  rhoer(npt)=electron number density (bohr^-3)
!!  rsepts(npt)=corresponding Wigner-Seitz radii, precomputed
!!  rhopr(npt)=positron number density (bohr^-3)
!!  ucvol=volume of the unit cell
!!
!! OUTPUT
!!  lifetime=1/lambda (picoseconds).
!!           lambda=Pi*r0^2*c*Int(rhoe(r)*rhop(r)*g0(0;rhoe;rhop)*dr)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine lifetime_psn(lifetime,npt,rhoer,rsepts,rhopr,rsppts,ucvol)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npt
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: lifetime
!arrays
 real(dp),intent(in) :: rhoer(npt),rhopr(npt),rsepts(npt),rsppts(npt)

!Local variables-------------------------------
!cc is light celerity in angstrom/second
!r0sq is Bohr radius squared
!scalars
 integer :: ipt
 real(dp),parameter :: cc=3.0d0*10.0d0**18/0.529177
 real(dp),parameter :: r0sq=2.817940285*10.0d0**(-5)/0.529177
 real(dp) :: g0,g0e,g0p,g1e,g1p,g2e,g2p,ke,kp,nbe,nbp,rhoe,rhop,rse,rsp
!arrays
 real(dp) :: ff(npt)

! *************************************************************************

 lifetime = 0.d0
 nbe =  0.d0
 nbp = 0.d0
!Loop over grid points
 do ipt=1,npt
  rse=rsepts(ipt)
  rsp=rsppts(ipt)
  rhoe=rhoer(ipt)
  rhop=rhopr(ipt)

  if (rhoe > rhop) then
   ke = -1.d0/6.d0*rse*(0.51+2.d0*0.65*rse &
&   -2.5*0.51*rse**(1.d0/2.d0) &
&   +3.d0*0.176*rse**2)
   g0e = 1.d0+1.23*rse+0.8295*rse**(1.5d0)-1.26*rse**2 &
&   +0.3286*rse**(2.5d0)+rse**3/6.d0
   g1e = 1.d0+0.51*rse+0.65*rse**2-0.51*rse**(2.5d0)+0.176*rse**3
   g2e = 1.d0+0.6*rse+0.63*rse**2-0.48*rse**(2.5d0)+0.167*rse**3


   g0 = 1.d0/(rhoe**3)*(2.d0*ke-6.d0*g1e+8.d0*g2e-2.d0*g0e)*rhop**3 &
&   +1.d0/(rhoe**2)*(-3.d0*ke+11.d0*g1e-16.d0*g2e+5.d0*g0e)*rhop**2 &
&   +1.d0/rhoe*(ke-4.d0*g1e+8.d0*g2e-4.d0*g0e)*rhop &
&   +g0e
  else
   kp = -1.d0/6.d0*rsp*(0.51+2.d0*0.65*rsp &
&   -2.5*0.51*rsp**(1.d0/2.d0) &
&   +3.d0*0.176*rsp**2)
   g0p = 1.d0+1.23*rsp+0.8295*rsp**(1.5d0)-1.26*rsp**2 &
&   +0.3286*rsp**(2.5d0)+rsp**3/6.d0
   g1p = 1.d0+0.51*rsp+0.65*rsp**2-0.51*rsp**(2.5d0)+0.176*rsp**3
   g2p = 1.d0+0.6*rsp+0.63*rsp**2-0.48*rsp**(2.5d0)+0.167*rsp**3


   g0 = 1.d0/(rhop**3)*(2.d0*kp-6.d0*g1p+8.d0*g2p-2.d0*g0p) * rhoe**3 &
&   +1.d0/(rhop**2)*(-3.d0*kp+11.d0*g1p-16.d0*g2p+5.d0*g0p)*rhoe**2 &
&   +1.d0/rhop*(kp-4.d0*g1p+8.d0*g2p- 4.d0*g0p) * rhoe &
&   +g0p
  end if

  ff(ipt)= rhoe*rhop*g0
  nbe = nbe + rhoe
  nbp = nbp + rhop
  lifetime = lifetime + ff(ipt)

 end do

 lifetime = lifetime * ucvol / npt
 nbe = nbe * ucvol / npt
 nbp = nbp * ucvol / npt

 lifetime = 1.0d0/(pi*r0sq*r0sq*cc*lifetime)

!
end subroutine lifetime_psn
!!***
