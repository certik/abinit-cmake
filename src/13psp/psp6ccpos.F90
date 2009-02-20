!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp6ccpos
!! NAME
!! psp6ccpos
!!
!! FUNCTION
!! Compute the core charge density, for use in the positron lifetime calculation
!! following the function definition valid for the format 6 of pseudopotentials. !!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (AF,GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mmax=maximum number of points in real space grid in the psp file
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  rchrg=cut-off radius for the core density
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  vhtnzc(mmax) = hartree potential induced by density tnzc (pseudo core density + nucleus)
!!  xccc1d(n1xccc,6)= 1D core charge function and its five first derivatives
!!
!! PARENTS
!!      psp6in
!!
!! CHILDREN
!!      calc_psden_log,calc_vhtnzc,smooth,spline,splint
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine psp6ccpos(mmax,n1xccc,rchrg,xccc1d,vhtnzc,znucl)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax,n1xccc
 real(dp),intent(in) :: rchrg,znucl
!arrays
 real(dp),intent(out) :: vhtnzc(mmax),xccc1d(n1xccc,6)

!Local variables-------------------------------
!scalars
 integer :: i1xccc,irad
 real(dp) :: der1,dern,rcut
!arrays
 real(dp),allocatable :: ff(:),ff1(:),ff2(:),ff3(:),gg(:),gg1(:),gg2(:),gg3(:)
 real(dp),allocatable :: gg4(:),rad(:),tnc(:),work(:),xx(:)

!**********************************************************************

 allocate(ff(mmax),ff1(mmax),ff2(mmax),ff3(mmax),rad(mmax))
 allocate(gg(n1xccc),gg1(n1xccc),gg2(n1xccc),gg3(n1xccc),gg4(n1xccc),&
& work(n1xccc),xx(n1xccc))

!
!read from pp file the model core charge (ff) and first (ff1) and
!second (ff2) derivative on logarithmic mesh mmax; rad is the radial grid
!the input functions contain the 4pi factor, it must be rescaled.

 do irad=1,mmax
  read(tmp_unit,*) rad(irad),ff(irad),ff1(irad),ff2(irad)
  ff(irad)=ff(irad)/4.d0/pi
  ff1(irad)=ff1(irad)/4.d0/pi
  ff2(irad)=ff2(irad)/4.d0/pi
 end do

 allocate(tnc(mmax))
 rcut=rchrg
 call calc_psden_log(tnc,mmax,ff,rcut,rad)
 call calc_vhtnzc(tnc,rcut,vhtnzc,mmax,rad,znucl)
 deallocate(tnc)



 rad(1)=0.d0


!calculate third derivative ff3 on logarithmic grid
 der1=ff2(1)
 dern=ff2(mmax)
 call spline(rad,ff1,mmax,der1,dern,ff3,work)

!generate uniform mesh xx in the box cut by rchrg:

 do i1xccc=1,n1xccc
  xx(i1xccc)=(i1xccc-1)* rchrg/dble(n1xccc-1)
 end do
!
!now interpolate core charge and derivatives on the uniform grid
!
!core charge, input=ff,  output=gg
 call splint(mmax,rad,ff,ff2 ,n1xccc,xx,gg)

!first derivative input=ff1, output=gg1
 call splint(mmax,rad,ff1,ff3,n1xccc,xx,gg1)

!normalize gg1
 gg1(:)=gg1(:)*rchrg

!now calculate second to fourth derivative by forward differences
!to avoid numerical noise uses a smoothing function

 call smooth(gg1,n1xccc,10)

 gg2(n1xccc)=0.0
 do i1xccc=1,n1xccc-1
  gg2(i1xccc)=(gg1(i1xccc+1)-gg1(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg2,n1xccc,10)

 gg3(n1xccc)=0.0
 do i1xccc=1,n1xccc-1
  gg3(i1xccc)=(gg2(i1xccc+1)-gg2(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg3,n1xccc,10)

 gg4(n1xccc)=0.0
 do i1xccc=1,n1xccc-1
  gg4(i1xccc)=(gg3(i1xccc+1)-gg3(i1xccc))*dble(n1xccc-1)
 end do

 call smooth(gg4,n1xccc,10)

!write on xcc1d
 xccc1d(:,1)=gg(:)
 xccc1d(:,2)=gg1(:)
 xccc1d(:,3)=gg2(:)
 xccc1d(:,4)=gg3(:)
 xccc1d(:,5)=gg4(:)

!WARNING : fifth derivative not yet computed
 xccc1d(:,6)=zero

!DEBUG
!note: the normalization condition is the following:
!4pi rchrg /dble(n1xccc-1) sum xx^2 xccc1d(:,1) = qchrg
!
!norm=0.d0
!do i1xccc=1,n1xccc
!norm = norm + 4.d0*pi*rchrg/dble(n1xccc-1)*&
!&             xx(i1xccc)**2*xccc1d(i1xccc,1)
!end do
!write(6,*) ' norm core=',norm
!
!write(6,*)' psp1cc : output of core charge density and derivatives '
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

 deallocate(ff,ff1,ff2,ff3,gg,gg1,gg2,gg3,gg4,rad,work,xx)

end subroutine psp6ccpos


!!***
