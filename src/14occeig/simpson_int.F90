!{\src2tex{textfont=tt}}
!!****f* ABINIT/simpson_int
!! NAME
!! simpson_int
!!
!! FUNCTION
!!   Simpson integral of input function
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mbessint=max number of points on grid for integral
!!  bessargmax=max point to which we will integrate
!!  bessint_delta = space between integral arguments
!!  mlang= max angular momentum
!!
!! OUTPUT
!!  spl_bessint=array of integrals
!!
!! NOTES
!!
!! PARENTS
!!      eli_lambda_1d,mka2f,mka2fQgrid,mka2f_tr,recip_ylm
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine simpson_int(nsimpson,simp_delta,simp_funct,simp_res)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsimpson
 real(dp),intent(in) :: simp_delta
!arrays
 real(dp),intent(in) :: simp_funct(nsimpson)
 real(dp),intent(out) :: simp_res(nsimpson)

!Local variables -------------------------
!scalars
 integer :: ii,iintarg
 real(dp) :: coef1,coef2,coef3

! *********************************************************************

 if (nsimpson < 6) then
  write (*,*) ' simpson_int : number of points in function must be >=6'
  stop
 end if

!-----------------------------------------------------------------
!Simpson integral of input function
!-----------------------------------------------------------------

 coef1 = 9.0_dp  / 24.0_dp
 coef2 = 28.0_dp / 24.0_dp
 coef3 = 23.0_dp / 24.0_dp

!first point is 0: don t store it

!do integration equivalent to Simpson O(1/N^4)
!from NumRec in C p 134  NumRec in Fortran p 128
 simp_res(1) =               coef1*simp_funct(1)
 simp_res(2) = simp_res(1) + coef2*simp_funct(2)
 simp_res(3) = simp_res(2) + coef3*simp_funct(3)

 do ii=4, nsimpson-3
  simp_res(ii) = simp_res(ii-1) + simp_funct(ii)
 end do

 simp_res(nsimpson-2) = simp_res(nsimpson-3) + coef3*simp_funct(nsimpson-2)
 simp_res(nsimpson-1) = simp_res(nsimpson-2) + coef2*simp_funct(nsimpson-1)
 simp_res(nsimpson  ) = simp_res(nsimpson-1) + coef1*simp_funct(nsimpson  )

!DEBUG
!write (*,*) ' simpson int ', simp_res(nsimpson, &
!&           simp_res(nsimpson) * bessint_delta
!ENDDEBUG

 simp_res(:) = simp_res(:) * simp_delta

end subroutine simpson_int
!!***
