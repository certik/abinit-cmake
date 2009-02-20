!{\src2tex{textfont=tt}}
!!****f* ABINIT/geteexc_cc
!! NAME
!! geteexc_cc
!!
!! FUNCTION
!! For the calculation of the exchange-correlation energy using
!! the fluctuation-dissipation theorem:
!! Calculate the potential energy associated with the susceptibility
!! matrix in reciprocal space, with a real-space cutoff Coulomb interaction.
!! Unlike for the true, infinite-range Coulomb interaction there are
!! contributions from $\vec G=0$ only if the $\vec G=0$-element of the
!! susd array is nonzero.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (MF).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gsq(npwdiel)=squares of G vectors
!! npwdiel=number of plane waves
!! npw_tiny=tiny number of plane waves (dimension of energy)
!! optextrap=0 no extrapolation of G=0 term
!!          =1 extrapolate the G=0 term
!! rcut_coulomb=cutoff radius for Coulomb interaction (bohr)
!! susd(npwdiel)=the (real) diagonal of the susceptibility matrix
!!
!! OUTPUT
!! energy(npw_tiny)=trace of susceptibility matrix times coulomb interaction
!! energy_raw=dto, but without G=0 contribution
!!
!! PARENTS
!!      acfd_dyson,xcacfd
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine geteexc_cc(energy,energy_raw,gsq,npwdiel,npw_tiny,optextrap,rcut_coulomb,susd)

 use defs_basis

 implicit none

!Arguments-------------------------------------
!scalars
 integer,intent(in) :: npw_tiny,npwdiel,optextrap
 real(dp),intent(in) :: rcut_coulomb
 real(dp),intent(out) :: energy_raw
!arrays
 real(dp),intent(in) :: gsq(npwdiel),susd(npwdiel)
 real(dp),intent(out) :: energy(npw_tiny)

!Local variables-------------------------------
!scalars
 integer :: iorder,ipw,ipwnull
 real(dp) :: trace

! *************************************************************************

!Compute trace without G=0 term
 trace=0._dp
 do ipw=1,npwdiel
  if(gsq(ipw) > 1.d-12) then
   trace=trace+susd(ipw)*(1._dp-cos(sqrt(gsq(ipw))*rcut_coulomb))/gsq(ipw)
  else
   ipwnull=ipw
  end if
 end do

!DEBUG
!Note: optextrap is used to check whether the extrapolation is consistent with
!what is intended for. Technically, zero susd needs no extrapolation,
!while nonzero susd does.
!if(optextrap==0) then
!if(susd(ipwnull)/=0._dp) then
!write(6,'(1x,a)') '%geteexc_cc: BUG: no extrapolation to zero G terms assumed, but'
!write(6,'(3x,a,i5,a,es14.7,a)') 'for ipwnull=',ipwnull,' have susd(ipwnull)= ',susd(ipwnull),&
!&    ' /= 0.'
!write(6,'(3x,a)') 'Action: Check that this makes sense !'
!end if
!else if(optextrap==1) then
!if(susd(ipwnull)==0._dp) then
!write(6,'(1x,a)') '%geteexc_cc: BUG: extrapolation to zero G terms assumed, but'
!write(6,'(3x,a,i5,a,es14.7,a)') 'for ipwnull=',ipwnull,' have susd(ipwnull)= ',susd(ipwnull),&
!&    ' == 0.'
!write(6,'(3x,a)') 'Action: Check that this makes sense !'
!end if
!else
!write(6,'(1x,a,i3,a)') '%geteexc_cc: BUG: using undefined optextrap=',optextrap,'. Stopping.'
!stop
!end if
!ENDDEBUG

!Add G=0 term
 energy(1)=four_pi*(trace+0.5_dp*rcut_coulomb**2*susd(ipwnull))
 energy_raw=four_pi*trace

end subroutine geteexc_cc
!!***
