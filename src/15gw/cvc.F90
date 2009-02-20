!{\src2tex{textfont=tt}}
!!****f* ABINIT/cvc
!! NAME
!! cvc
!!
!! FUNCTION
!! Set up table of length of |q+G| for coulombian potential
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gvec(3,npwvec)=coordinates of G vectors
!! gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!! iq=specific q point in BZ
!! npwvec=number of planewaves
!! nq=number of q points in BZ
!! q=coordinates of q points in BZ
!!
!! OUTPUT
!! qplusg(npwvec)=norm of q+G vector
!!
!! PARENTS
!!      cppm4par,screening
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cvc(nq,iq,q,npwvec,gvec,gprimd,qplusg)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq,npwvec,nq
!arrays
 integer,intent(in) :: gvec(3,npwvec)
 real(dp),intent(in) :: gprimd(3,3),q(3,nq)
 real(dp),intent(out) :: qplusg(npwvec)

!Local variables ------------------------------
!scalars
 integer :: ig,ii
!arrays
 real(dp) :: gmet(3,3),gpq(3)

!************************************************************************

 ! Compute reciprocal space metrics
 do ii=1,3
  gmet(ii,:)=gprimd(1,ii)*gprimd(1,:)+&
&            gprimd(2,ii)*gprimd(2,:)+&
&            gprimd(3,ii)*gprimd(3,:)
 end do

 if (ALL(ABS(q(:,iq))<1.e-3)) then
  qplusg(1)=two_pi*SQRT(DOT_PRODUCT(q(:,iq),MATMUL(gmet,q(:,iq))))
  do ig=2,npwvec
   gpq(:)=gvec(:,ig)
   qplusg(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
  end do
 else
  do ig=1,npwvec
   gpq(:)=gvec(:,ig)+q(:,iq)
   qplusg(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
  end do
 end if

end subroutine cvc
!!***
