!{\src2tex{textfont=tt}}
!!****f* ABINIT/bz1
!! NAME
!! bz1
!!
!! FUNCTION
!! Given a k-point not necessarily in the 1st Brillouin Zone (but quite
!! close to it) find kBZ and G where kBZ is in the first BZ and G
!! a RL vector such that k=kBZ+G.  G is assumed to have all components
!! either -1, 0, 1, though in fact -2 to 2 are checked and if -2 or 2
!! is found as a component of G an error is flagged.
!! K, G and kBZ are expressed in RL coords.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gmet(3,3==reciprocal space metrics
!!
!! OUTPUT
!!  g(3)=integer vector to bring initial k vector in BZ
!!
!! SIDE EFFECTS
!!  k(3)=in/out k vector
!!
!!
!! PARENTS
!!      findq
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine bz1(k,g,gmet)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(out) :: g(3)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(inout) :: k(3)

!Local variables ------------------------------
!scalars
 integer :: ix,iy,iz
 real(dp) :: x,xmin
!arrays
 real(dp) :: kmg(3)

!************************************************************************
 xmin=99999999.99

 do ix=-2,2
  do iy=-2,2
   do iz=-2,2

    kmg(1)=k(1)-ix
    kmg(2)=k(2)-iy
    kmg(3)=k(3)-iz
    x = two_pi*SQRT(DOT_PRODUCT(kmg,MATMUL(gmet,kmg)))

!   Heuristic method for choosing the right G when two Gs will do
!   (avoid the G such that kBZ has all components with same sign
!   otherwise k1BZ+k2BZ may need a G0 with components other than
!   -1,0,1 in RL coordinates to translate it back to the 1st BZ)
!   so add a small "penalty" to x if all 3 components of G have
!   the same sign
    x=x+0.001*abs(sign(one,kmg(1))+sign(one,kmg(2))+sign(one,kmg(3)))

    if(xmin-x>1e-12) then
     g(1)=ix
     g(2)=iy
     g(3)=iz
     xmin=x
    end if
   end do  ! iz
  end do ! iy
 end do ! iz

!Check that G does not have a component of-2 or 2
 if(abs(g(1))==2) stop
 if(abs(g(2))==2) stop
 if(abs(g(3))==2) stop

!Subtract G from k
 k(1)=k(1)-g(1)
 k(2)=k(2)-g(2)
 k(3)=k(3)-g(3)

end subroutine bz1
!!***
