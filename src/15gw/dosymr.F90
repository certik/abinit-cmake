!{\src2tex{textfont=tt}}
!!****f* ABINIT/dosymr
!! NAME
!! dosymr
!!
!! FUNCTION
!! Perform an inverse symmetry operation on a real-space vector in FFT units
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iinv= if 2 invert the vector, if 1, do not invert it
!!  ngfft(3)=1D dimensions of FFT grid
!!  op(3,3)=symmetry operation in reciprocal space
!!  r1(3)=input real space point
!!
!! OUTPUT
!!  r2(3)=symmetric/inverted real space point
!!
!! NOTES
!! If the OP matrices are in reciprocal-lattice units then
!! since the inverse-transpose of the OP matrix is the
!! corresponding matrix in real-space-lattice units, NFFT(I)/NFFT(J)
!! times the transpose is the matrix corresponding to the inverse
!! operation in FFT units.  Thus calling DOSYMR with OP and r1 performs
!! the operation r2=R**-1 r1 if OP corresponds to R.
!! IINV is 1 for non-inversion, 2 for inversion
!!
!! PARENTS
!!      surot
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dosymr(op,iinv,r1,ngfft,r2)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iinv
!arrays
 integer,intent(in) :: ngfft(3)
 real(dp),intent(in) :: op(3,3),r1(3)
 real(dp),intent(out) :: r2(3)

!Local variables-------------------------------
!scalars
 integer :: ii,imult,jj

! *************************************************************************

 r2(:)=zero

 imult=3-2*iinv
 do jj=1,3
  do ii=1,3
   r2(ii)=r2(ii)+imult*op(jj,ii)*ngfft(ii)*r1(jj)/float(ngfft(jj))
  end do
 end do

end subroutine dosymr
!!***
