!{\src2tex{textfont=tt}}
!!****f* ABINIT/fkin
!! NAME
!! fkin
!!
!! FUNCTION
!! Calculate the integral over frequency (up to prefactors) appearing
!! in the kinetic contribution to the bandgap energy in the GW approximation.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (YMN)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ff = (integer !) occupation number for the KS state with energy epsilon0.
!!  omegame0 = omega-epsilon0 (omega is the frequency at which the self-energy
!!    is calculated).
!!  otw1, otw2 = two plasmon poles.
!!  zcut = +/-i zcut is added to the denominators of kincontrib if their absolute
!!    values are below zcut (avoids spurious divergences).
!!
!! OUTPUT
!!  kincontrib = the integral over frequency (up to prefactors) appearing in the
!!   kinetic contribution to the bandgap energy in the GW approximation.
!!
!! NOTES
!!
!! TODO
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fkin(ff,kincontrib,omegame0,otw1,otw2,zcut)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ff
 real(dp),intent(in) :: omegame0,otw1,otw2,zcut
 real(dp),intent(out) :: kincontrib

!Local variables ------------------------------
!scalars
 real(dp) :: fact,x
 complex :: den1,den2

!*************************************************************************

 if (ff == 1) then

  x = omegame0+otw1
  if (abs(x) < zcut) then
   den1 = cmplx(x,-zcut)
  else
   den1 = x
  end if
  x = omegame0+otw2
  if (abs(x) < zcut) then
   den2 = cmplx(x,-zcut)
  else
   den2 = x
  end if

  fact = omegame0/(((otw1+otw2)+one)*(otw1*otw2))

 else

  x = omegame0-otw1
  if (abs(x) < zcut) then
   den1 = cmplx(x,zcut)
  else
   den1 = x
  end if
  x = omegame0-otw2
  if (abs(x) < zcut) then
   den2 = cmplx(x,zcut)
  else
   den2 = x
  end if

  fact = omegame0/(((otw1+otw2)-one)*(otw1*otw2))

 end if

 kincontrib = fact/(den1*den2)

end subroutine fkin

!!***
