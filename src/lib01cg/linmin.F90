!{\src2tex{textfont=tt}}
!!****f* ABINIT/linmin
!! NAME
!! linmin
!!
!! FUNCTION
!! minimizes a function along a gradient line:
!! first bracket the minimum then perform the minimization
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! dp_dum_vdp: function  to be minimized (return a dp from a vector of dp)
!! vdp_dum_vdp: derivative of f
!! itmax: number of iterations allowed
!!
!! OUTPUT
!! fmin: minimun value reached for dp_dum_vdp
!!
!! SIDE EFFECTS
!! grad: the gradient line along which the minimization is performed (not changed)
!! v: the starting and then ending point of the minimization
!!
!! WARNINGS
!!
!! NOTES
!!
!! PARENTS
!! cgpr
!!
!! CHILDREN
!! bracketing,brent
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine linmin(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,itmax,v,grad,fmin)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01cg, except_this_one => linmin
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 include "dummy_functions.inc"
!scalars
 integer,intent(in) :: itmax,nv1,nv2
 real(dp),intent(out) :: fmin
!arrays
 real(dp),intent(inout) :: grad(nv1,nv2),v(nv1,nv2)

!Local variables-------------------------------
!scalars
 real(dp),parameter :: maglimit=10000.0_dp,tol=tol8*tol8*tol8
 real(dp) :: a,b,fa,fb,fx,x,xmin
!no_abirules

!************************************************************************
 a=zero
 x=ninth*real(1e-4,dp)
 call bracketing (nv1,nv2,dp_dum_v2dp,v,grad,a,x,b,fa,fx,fb)
!DEBUG
!write(6,*) 'linmin (01cg) : linmin after bracketing'
!write(6,*) 'linmin (01cg) : point',a,x,b,'value',fa,fx,fb
!ENDDEBUG
 fmin =brent(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,6,v,grad,a,x,b,tol,xmin)

 end subroutine linmin
!!***
