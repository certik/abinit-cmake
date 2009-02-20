!{\src2tex{textfont=tt}}
!!****f* ABINIT/cc_derivatives
!! NAME
!! cc_derivatives
!!
!! FUNCTION
!! subroutine to spline the core charge and get derivatives
!!   extracted from previous version of psp6cc_drh
!! input on log grid, and splined to regular grid between
!!  0 and rchrg
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (AF,DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rchrg=cut-off radius for the core density
!!  rad=radial grid points
!!  ff=core charge at points in rad
!!  ff1=first derivative of ff on log grid
!!  ff2=second derivative of ff on log grid
!!
!!
!! OUTPUT
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! PARENTS
!!      psp6in
!!
!! CHILDREN
!!      spline,splint
!!
!! NOTES
!! Test version by DRH - requires very smooth model core charge
!!
!! SOURCE
  
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cc_derivatives(rad,ff,ff1,ff2,mmax,n1xccc,rchrg,xccc1d)


 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
! scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(in) :: rchrg
!arrays
 real(dp),intent(in) :: rad(mmax),ff(mmax),ff1(mmax),ff2(mmax)
 real(dp),intent(out) :: xccc1d(n1xccc,6)

!Local variables-------------------------------
! scalars
 integer :: i1xccc
 real(dp) :: der1,dern
!arrays
 real(dp),allocatable :: ff3(:),ff4(:),gg(:),gg1(:),gg2(:)
 real(dp),allocatable :: gg3(:),gg4(:),work(:),xx(:)
  
! *************************************************************************
 allocate(ff3(mmax),ff4(mmax))
 allocate(gg(n1xccc),gg1(n1xccc),gg2(n1xccc),gg3(n1xccc),gg4(n1xccc))
 allocate(work(mmax),xx(n1xccc))

 write (*,*) 'cc_derivatives : enter'

!calculate third derivative ff3 on logarithmic grid
 der1=ff2(1)
 dern=ff2(mmax)
 call spline(rad,ff1,mmax,der1,dern,ff3,work)

!calculate fourth derivative ff4 on logarithmic grid
 der1=0.d0
 dern=0.d0
 call spline(rad,ff2,mmax,der1,dern,ff4,work)

!generate uniform mesh xx in the box cut by rchrg:

 do i1xccc=1,n1xccc
  xx(i1xccc)=(i1xccc-1)* rchrg/dble(n1xccc-1)
 end do
!
!now interpolate core charge and derivatives on the uniform grid
!
!core charge, input=ff,  output=gg
 call splint(mmax,rad,ff,ff2,n1xccc,xx,gg)

!first derivative input=ff1, output=gg1
 call splint(mmax,rad,ff1,ff3,n1xccc,xx,gg1)

!normalize gg1
!gg1(:)=gg1(:)*rchrg

!second derivative input=ff2, output=gg2
 call splint(mmax,rad,ff2,ff4,n1xccc,xx,gg2)

!normalize gg2
!gg2(:)=gg2(:)*rchrg**2

!reallocate work otherwise the calls to spline crash
 deallocate (work)
 allocate (work(n1xccc))

!recalculate 3rd derivative consistent with spline fit to first derivative
!on linear grid
 der1=gg2(1)
 dern=gg2(n1xccc)
 call spline(xx,gg1,n1xccc,der1,dern,gg3,work)

!calculate 4th derivative consistent with spline fit to second derivative
!on linear grid
 der1=0.0d0
 dern=0.0d0
 call spline(xx,gg2,n1xccc,der1,dern,gg4,work)


!now calculate second to fourth derivative by forward differences
!to avoid numerical noise uses a smoothing function
!
!call smooth(gg1,n1xccc,10)

!gg2(n1xccc)=0.0
!do i1xccc=1,n1xccc-1
!gg2(i1xccc)=(gg1(i1xccc+1)-gg1(i1xccc))*dble(n1xccc-1)
!end do

!call smooth(gg2,n1xccc,10)

!gg3(n1xccc)=0.0
!do i1xccc=1,n1xccc-1
!gg3(i1xccc)=(gg2(i1xccc+1)-gg2(i1xccc))*dble(n1xccc-1)
!end do

!call smooth(gg3,n1xccc,10)

!gg4(n1xccc)=0.0
!do i1xccc=1,n1xccc-1
!gg4(i1xccc)=(gg3(i1xccc+1)-gg3(i1xccc))*dble(n1xccc-1)
!end do

!call smooth(gg4,n1xccc,10)

!write on xcc1d
!normalize to unit range usage later in program
 xccc1d(:,1)=gg(:)
 xccc1d(:,2)=gg1(:)*rchrg
 xccc1d(:,3)=gg2(:)*rchrg**2
 xccc1d(:,4)=gg3(:)*rchrg**3
 xccc1d(:,5)=gg4(:)*rchrg**4
!***drh test
!write(6,'(a,2i6)') 'drh:psp6cc_drh - mmax,n1xccc',mmax,n1xccc
!***end drh test


!DEBUG
!note: the normalization condition is the following:
!4pi rchrg /dble(n1xccc-1) sum xx^2 xccc1d(:,1) = qchrg
!
!norm=0.d0
!do i1xccc=1,n1xccc
!norm = norm + 4.d0*pi*rchrg/dble(n1xccc-1)*&
!&             xx(i1xccc)**2*xccc1d(i1xccc,1)
!end do
!write(6,*) ' norm=',norm
!
!write(6,*)' psp6cc_drh : output of core charge density and derivatives '
!write(6,*)'   xx          gg           gg1  '
!do i1xccc=1,n1xccc
!write(10, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,1),xccc1d(i1xccc,2)
!end do
!write(6,*)'   xx          gg2          gg3  '
!do i1xccc=1,n1xccc
!write(11, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,3),xccc1d(i1xccc,4)
!end do
!write(6,*)'   xx          gg4          gg5  '
!do i1xccc=1,n1xccc
!write(12, '(3es14.6)' ) xx(i1xccc),xccc1d(i1xccc,5),xccc1d(i1xccc,6)
!end do
!write(6,*)' psp1cc : debug done, stop '
!stop
!ENDDEBUG

 deallocate(ff3,ff4,gg,gg1,gg2,gg3,gg4,work,xx)

end subroutine cc_derivatives

!!***
