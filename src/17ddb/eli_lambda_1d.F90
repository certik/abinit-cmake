!{\src2tex{textfont=tt}}
!!****f* ABINIT/eli_lambda_1d
!!
!! NAME
!! eli_lambda_1d
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
!!   a2f_1d = 1D alpha2F function
!!   elph_ds = elphon dataset
!!   nmatsu = number of Matsubara frequencies
!!   tc = guess for critical temperature
!!
!! OUTPUT
!!   lambda_1d = coupling constant as a function of frequency
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      eliashberg_1d
!!
!! CHILDREN
!!      simpson_int
!!
!! NOTES
!!  lambda is used at points which are differences of Matsubara freqs,
!!  and hence is tabulated on points going through 0.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine eli_lambda_1d (a2f_1d,elph_ds,lambda_1d,nmatsu,tc)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_14occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
 real(dp),intent(in) :: tc
 type(elph_type),intent(in) :: elph_ds
!arrays
 real(dp),intent(in) :: a2f_1d(elph_ds%na2f)
 real(dp),intent(out) :: lambda_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: imatsu,iomega
 real(dp) :: nu_matsu,nu_matsu2,omega
!arrays
 real(dp) :: lambda_int(elph_ds%na2f),tmplambda(elph_ds%na2f)

! *********************************************************************

 do imatsu=-nmatsu,nmatsu

  nu_matsu = (two*imatsu)*pi*tc
  nu_matsu2 = nu_matsu*nu_matsu

  tmplambda(:) = zero
  omega=elph_ds%domega
  do iomega=2,elph_ds%na2f
   tmplambda(iomega) = a2f_1d(iomega) * two * omega / (nu_matsu2 + omega*omega)
   omega=omega+elph_ds%domega
  end do
  call simpson_int(elph_ds%na2f,elph_ds%domega,tmplambda,lambda_int)

  lambda_1d(imatsu) = lambda_int(elph_ds%na2f)

 end do


end subroutine eli_lambda_1d
!!***
