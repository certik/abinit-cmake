!{\src2tex{textfont=tt}}
!!****f* ABINIT/init_bess_spl
!! NAME
!! init_bess_spl
!!
!! FUNCTION
!! Pre-calculate the j_v(y) for recip_ylm on regular grid
!!     NOTE: spherical Bessel function small j!
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mbess=  max number of points on grid for integral
!!  bessargmax= max point to which we will integrate
!!  bessint_delta = space between integral arguments
!!  mlang=  max angular momentum
!!
!! OUTPUT
!!  bess_spl=array of integrals
!!  bess_spl_der=array of derivatives of integrals
!!  x_bess=coordinates of points belonging to the grid
!!
!! NOTES
!!
!! PARENTS
!!      partial_dos_fractions,wffile
!!
!! CHILDREN
!!      besjm,spline
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine init_bess_spl(mbess,bessargmax,bessint_delta,mlang,&
&    bess_spl,bess_spl_der,x_bess)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mbess,mlang
 real(dp),intent(in) :: bessargmax,bessint_delta
!arrays
 real(dp),intent(out) :: bess_spl(mbess,mlang),bess_spl_der(mbess,mlang)
 real(dp),intent(out) :: x_bess(mbess)

!Local variables -------------------------
! function calls (NumRec)
!scalars
 integer :: ii,iintarg,ll
 real(dp) :: cosdelta_bess,sindelta_bess,yp1,ypn
!arrays
 real(dp) :: work(mbess)
 real(dp),allocatable :: cosbessx(:),sinbessx(:)

! *********************************************************************


!DEBUG
!write (*,*) 'init_bess_ylm enter '
!ENDDEBUG

!-----------------------------------------------------------------
!Bessel function into array
!-----------------------------------------------------------------

!integration grid is nfiner times finer than the interpolation grid
 allocate (sinbessx(mbess))
 allocate (cosbessx(mbess))

 sindelta_bess = sin(bessint_delta)
 cosdelta_bess = cos(bessint_delta)
!
!could be done by chain rule for cos sin (is it worth it?) but
!precision problems as numerical errors are propagated.
!
 do iintarg=1,mbess
  x_bess(iintarg) = (iintarg-1)*bessint_delta
  sinbessx(iintarg) = sin(x_bess(iintarg))
  cosbessx(iintarg) = cos(x_bess(iintarg))

! x_bess(iintarg) = x_bess(iintarg-1)+bessint_delta
! !  get sin and cos of x_bess arguments
! sinbessx(iintarg) = sinbessx(iintarg-1)*cosdelta_bess &
! &                     + cosbessx(iintarg-1)*sindelta_bess
! cosbessx(iintarg) = cosbessx(iintarg-1)*cosdelta_bess &
! &                     - sinbessx(iintarg-1)*sindelta_bess
 end do

!write (*,*) 'x_bess = ', x_bess
!
!fill bess_spl array
!
 do ll=0,mlang-1

  call besjm(one,bess_spl(:,ll+1),cosbessx,   &
&  ll,mbess,sinbessx,x_bess)
! 
! call spline to get 2nd derivative (reuse in splint later)
! 
  yp1 = zero
  ypn = zero
  call spline (x_bess, bess_spl(:,ll+1), mbess, yp1, ypn, &
&  bess_spl_der(:,ll+1),work)
 end do

!DEBUG
!write (*,*) ' bess funct  0   1   2   3   4'
!do iintarg=1,mbess
!write (*,*) x_bess(iintarg), (bess_spl(iintarg,ll),ll=1,mlang)
!end do
!ENDDEBUG

 deallocate (sinbessx,cosbessx)

end subroutine init_bess_spl
!!***
