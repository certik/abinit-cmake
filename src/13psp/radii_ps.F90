!{\src2tex{textfont=tt}}
!!****f* ABINIT/radii_ps
!! NAME
!! radii_ps
!!
!! FUNCTION
!!  (to be completed)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      smoothvlocal
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine radii_ps( vps, rofi, zval, nrval, lmxkb, nrgauss, rgauss, rgauss2)

 use defs_basis
!     This routine returns the maximum radius for the
!     Kleinman-Bylander projectors with a standard choice
!     of the local potential.
!     Check also at which radius the asymptotic 2*Zval/r
!     behaviour is achieved.
!     D. Sanchez-Portal, Aug. 1998

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmxkb,nrval
 integer,intent(out) :: nrgauss
 real(dp),intent(in) :: zval
 real(dp),intent(out) :: rgauss,rgauss2
!arrays
 real(dp),intent(in) :: rofi(:),vps(:,0:)

!Local variables ------------------------------
!scalars
 integer :: ir,l,nrgauss2
 real(dp),parameter :: eps=1.0d-4
 real(dp) :: dincv,r

! *************************************************************************

!**   Iterate over the possible local potentials**
 rgauss   = 0.0d0
 rgauss2  = 0.0d0
 nrgauss  = 0
 nrgauss2 = 0

 do l = 0, lmxkb-1
  do ir = nrval, 2, -1
   dincv = abs( vps(ir,l) - vps(ir,lmxkb) )
   if( dincv .gt. eps ) then
    exit
   end if
  end do
  rgauss  = max( rofi(ir), rgauss )
  nrgauss = max( ir, nrgauss )
 end do
!
!New: Use all potentials, not just l=0, since
!potentials with larger l can converge later...
!
 do l = 0, lmxkb
  do ir = nrval, 2, -1
   r = rofi(ir)
   dincv = abs( vps(ir,l)*r + 2.0d0*zval )
   if( dincv .gt. eps ) then
    exit
   end if
  end do
! write(6,'(a,i1,a,f8.4)') &
! &    'V l=', l,' = -2*Zval/r beyond r=', rofi(ir)
  rgauss2  = max( rofi(ir), rgauss2 )
  nrgauss2 = max( ir, nrgauss2 )
 end do

 if( lmxkb .eq. 0 ) then
  rgauss  = rgauss2
  nrgauss = nrgauss2
 end if

!write(6,'(a,f8.4)') 'All V_l potentials equal beyond r=', rgauss
!write(6,'(a)') &
!&  'This should be close to max(r_c) in ps generation'
!write(6,'(a,f8.4)') &
!&  'All pots = -2*Zval/r beyond r=', rgauss2

end subroutine radii_ps
!---
!!***
