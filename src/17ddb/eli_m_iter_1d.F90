!{\src2tex{textfont=tt}}
!!****f* ABINIT/eli_m_iter_1d
!!
!! NAME
!! eli_m_iter_1d
!!
!! FUNCTION
!!  Find largest eigenvalue of M matrix, to deduce superconducting Tc
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
!!   z_1d = renormalization Z as a function of frequency
!!
!! OUTPUT
!!   maxeigval = estimation for maximum eigenvalue of M
!!
!! SIDE EFFECTS
!!   delta_1d = imaginary gap function as a function of frequency
!!
!! PARENTS
!!
!! CHILDREN
!!      eli_app_m_1d
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine eli_m_iter_1d (delta_1d,lambda_1d,maxeigval,nmatsu,tc,z_1d)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_17ddb, except_this_one => eli_m_iter_1d
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nmatsu
 real(dp),intent(in) :: tc
 real(dp),intent(out) :: maxeigval
!arrays
 real(dp),intent(in) :: lambda_1d(-nmatsu:nmatsu),z_1d(-nmatsu:nmatsu)
 real(dp),intent(inout) :: delta_1d(-nmatsu:nmatsu)

!Local variables-------------------------------
!scalars
 integer :: iiterm,imatsu,jmatsu,nfilter,ngeteig
 real(dp) :: dnewnorm,dnorm,fact,nu_matsu,nu_matsu2,omega
!arrays
 real(dp) :: delta_old(-nmatsu:nmatsu),delta_tmp(-nmatsu:nmatsu)

! *********************************************************************

 nfilter = 10
 ngeteig = 10


!
!1) apply M matrix enough times to filter out largest eigenvalue
!
 do iiterm=1,nfilter
  call eli_app_m_1d (delta_1d,lambda_1d,nmatsu,tc,z_1d)

! DEBUG
! dnorm=zero
! do imatsu=-nmatsu,nmatsu
! dnorm = dnorm + delta_1d(imatsu)*delta_1d(imatsu)/(two*imatsu+one)
! end do
! dnorm = sqrt(dnorm)
! write (*,*) 'eli_m_iter_1d : dnorm ', dnorm
! ENDDEBUG
 end do

!
!2) calculate norm
!
 dnorm=zero
 do imatsu=-nmatsu,nmatsu
  dnorm = dnorm + delta_1d(imatsu)*delta_1d(imatsu)/abs(two*imatsu+one)
 end do
 dnorm = sqrt(dnorm)

!normalize delta_1d
 delta_1d(:) = delta_1d(:) / dnorm

 delta_old = delta_1d

!DEBUG
!dnewnorm=zero
!do imatsu=-nmatsu,nmatsu
!dnewnorm = dnewnorm + delta_1d(imatsu)*delta_1d(imatsu)/abs(two*imatsu+one)
!end do
!dnewnorm = sqrt(dnewnorm)
!write (*,*) 'eli_m_iter_1d : dnewnorm1 ', dnewnorm
!ENDDEBUG

!
!3) re-apply M matrix ngeteig times
!
 do iiterm=1,ngeteig
  call eli_app_m_1d (delta_1d,lambda_1d,nmatsu,tc,z_1d)
! DEBUG
! dnewnorm=zero
! do imatsu=-nmatsu,nmatsu
! dnewnorm = dnewnorm + delta_1d(imatsu)*delta_1d(imatsu)/abs(two*imatsu+one)
! end do
! dnewnorm = sqrt(dnewnorm)
! write (*,*) 'eli_m_iter_1d : dnewnorm ', dnewnorm
! 
! do imatsu=-nmatsu,nmatsu
! write (112,*) imatsu,delta_1d(imatsu)/delta_old(imatsu)
! end do
! write (112,*)
! delta_old = delta_1d
! ENDDEBUG
 end do

!
!4) calculate new norm and estimate eigenvalue
!
 dnewnorm=zero
 do imatsu=-nmatsu,nmatsu
  dnewnorm = dnewnorm + delta_1d(imatsu)*delta_1d(imatsu)/abs(two*imatsu+one)
 end do
 dnewnorm = sqrt(dnewnorm)

!maxeigval = exp ( log(dnewnorm/dnorm) / ngeteig )
 maxeigval = exp ( log(dnewnorm) / ngeteig )

 write (*,*) 'eli_m_iter_1d : maxeigval =' , maxeigval
!fact = exp(-log(maxeigval) * (ngeteig+nfilter))
 fact = exp(-log(maxeigval) * (ngeteig))
 do imatsu=-nmatsu,nmatsu
  delta_1d(imatsu) = delta_1d(imatsu) * fact
 end do

end subroutine eli_m_iter_1d
!!***
