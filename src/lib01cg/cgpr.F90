!{\src2tex{textfont=tt}}
!!****f* ABINIT/cgpr
!! NAME
!! cgpr
!!
!! FUNCTION
!! perform Polak-Ribiere conjugate gradient on a function f
!! implementation based on the cg recipe of "numerical recipe"
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
!! dtol: precision precision required for the minimization
!! itmax: number of iterations allowed (each linmin will be done with at max 10 times
!! this number
!!
!! OUTPUT
!! fmin: value of f at the minimum
!! lastdelta: absolute value of the last delta between steps
!! SIDE EFFECTS
!! v: vector on which minimization is to be performed, starting point
!! and resulting min
!!
!! WARNINGS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!! linmin
!! dotproduct
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cgpr(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,dtol,itmax,v,fmin,delta)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01cg, except_this_one => cgpr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 include "dummy_functions.inc"
!scalars
 integer,intent(in) :: itmax,nv1,nv2
 real(dp),intent(in) :: dtol
 real(dp),intent(out) :: delta,fmin
!arrays
 real(dp),intent(inout) :: v(nv1,nv2)

!Local variables-------------------------------
!scalars
 integer :: iiter
 real(dp) :: fv,gam,gscal,gscal2,sto
!arrays
 real(dp) :: grad0(nv1,nv2),grad1(nv1,nv2),grad2(nv1,nv2),grad3(nv1,nv2)
!no_abirules

!************************************************************************
 fv = dp_dum_v2dp(nv1,nv2,v(:,:))
 grad0(:,:) = -v2dp_dum_v2dp(nv1,nv2,v(:,:))
 grad1(:,:) = grad0(:,:)
 grad2(:,:) = grad0(:,:)
 do iiter=1,itmax
  call linmin(nv1,nv2,dp_dum_v2dp,v2dp_dum_v2dp,sub_dum_dp_v2dp_v2dp,itmax,v,grad0,fmin)
! return if the min is reached
  sto=dtol*(abs(fmin)+abs(fv)+tol14)
  delta=abs(fv-fmin)
  delta=abs(delta)
  if((delta.lt.sto).or.(iiter==itmax)) then
!  DEBUG
!  write(6,*) 'cgpr (01cg) : stop cond for cgpr:',sto,'delta:',delta,'fv:',fv,'fmin:',fmin
!  ENDDEBUG
   return
  end if
! a new step
  fv=fmin
  grad0(:,:)=v2dp_dum_v2dp(nv1,nv2,v(:,:))
  gscal=dotproduct(nv1,nv2,grad1(:,:),grad1(:,:))
  grad3(:,:)=grad0(:,:)+grad1(:,:)
  gscal2=dotproduct(nv1,nv2,grad3(:,:),grad0(:,:))
  gam=gscal2/gscal
  grad1(:,:)=-grad0(:,:)
  grad2(:,:)=grad1(:,:)+gam*grad2(:,:)
  grad0(:,:)=grad2(:,:)
! DEBUG
! write(6,*) 'cgpr (01cg) :================================================================================='
! write(6,*) 'cgpr (01cg) : step',iiter,'delta:',delta ,'fv',fv,'fmin',fmin
! write(6,*) 'cgpr (01cg) :================================================================================='
! ENDDEBUG
 end do
end subroutine cgpr

!!***
