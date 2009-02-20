!{\src2tex{textfont=tt}}
!!****f* ABINIT/eli_z_1d
!!
!! NAME
!! eli_z_1d
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   lambda_1d = coupling constant as a function of frequency
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!
!! OUTPUT
!!   z_1d = renormalizing Z as a function of frequency
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      eliashberg_1d
!!
!! CHILDREN
!!
!! NOTES
!!  Because Z only depends on lambda(n-n´), and lambda(omega)
!!   is an even function, Z is symmetrical in n and -n
!!   hence only calculate for n>0 and complete the rest
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine eli_z_1d (lambda_1d,nmatsu,tc,z_1d)

 use defs_basis
 use defs_datatypes
 use defs_elphon

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
 real(dp),intent(in) :: tc
!arrays
 real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu)
 real(dp),intent(out) :: z_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,jmatsu
 real(dp) :: nu_matsu,nu_matsu2,omega

! *********************************************************************


 do imatsu=0,nmatsu

  z_1d(imatsu) = zero
! count $\mathrm{sign}(omega_{Matsubara})$
  do jmatsu=-nmatsu+imatsu,-1
   z_1d(imatsu) = z_1d(imatsu) - lambda_1d(imatsu-jmatsu)
  end do
  do jmatsu=0,nmatsu
   z_1d(imatsu) = z_1d(imatsu) + lambda_1d(imatsu-jmatsu)
  end do

! NOTE: the pi*Tc factor in Z cancels the one in the Matsubara frequency.
  z_1d(imatsu) = one + z_1d(imatsu) / (two*imatsu+one)
  z_1d(-imatsu) = z_1d(imatsu)
 end do

end subroutine eli_z_1d
!!***
