!{\src2tex{textfont=tt}}
!!****f* ABINIT/dosym
!! NAME
!! dosym
!!
!! FUNCTION
!! Perform a symmetry operation k2=RI k1
!! If the OP matrices are in reciprocal-lattice units then the vectors
!! must be too.
!! iinv is 1 for non-inversion, 2 for inversion
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
!!  k1=input k vector
!!  op=symmetry matrix
!!
!! OUTPUT
!!  k2= symmetric/inverted k vector
!!
!! PARENTS
!!      assemblychi0q0_sym,cchi0q0,findnq,findq,identk,identq
!!      setup_little_group,surot
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dosym(op,iinv,k1,k2)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iinv
!arrays
 real(dp),intent(in) :: k1(3),op(3,3)
 real(dp),intent(out) :: k2(3)

!Local variables-------------------------------
!scalars
 integer :: ii,imult,jj

! *************************************************************************

 k2(:)=zero

 imult=3-2*iinv

 do jj=1,3
  do ii=1,3
   k2(ii)=k2(ii)+imult*op(ii,jj)*k1(jj)
  end do
 end do

end subroutine dosym
!!***
