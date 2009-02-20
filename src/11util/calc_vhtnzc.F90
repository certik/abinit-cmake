!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_vhtnzc
!! NAME
!! calc_vhtnzc
!!
!! FUNCTION
!! Compute vh(tnZc)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GJ, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mesh= dimension of nc
!!  nc= density to be pseudized
!!  rad(mesh)=radial mesh
!!  rc= cut-off radius
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  vhtnzc(mesh) = hartree potential induced by density tnzc (pseudo core density + nucleus)
!!
!! SIDE EFFECT
!!  rc= cut-off radius    (XG050317 : divided by 5 at output ?! strange)
!!
!! NOTES
!!
!! PARENTS
!!      psp6ccpos
!!
!! CHILDREN
!!      ctrap
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_vhtnzc(nc,rc,vhtnzc,mesh,rad,znucl)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util, except_this_one => calc_vhtnzc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mesh
 real(dp),intent(in) :: znucl
 real(dp),intent(inout) :: rc
!arrays
 real(dp),intent(in) :: nc(mesh),rad(mesh)
 real(dp),intent(out) :: vhtnzc(mesh)

!Local variables-------------------------------
!scalars
 integer :: ir,irad,nc1
 real(dp) :: amesh,gnorm,intg,rc1,step,yp1,yp2,yp3
!arrays
 real(dp) :: shapefunc(mesh)
 real(dp),allocatable :: den1(:),den2(:),den3(:),den4(:),ff(:),nzc(:),rvhn(:)

! *************************************************************************

 allocate(nzc(mesh))

 shapefunc(1)=1.d0
 amesh=rad(2)/rad(1)
 rc=rc/5.0
 nc1=int(log(rc/rad(1))/log(amesh))+1
 rc1=rad(nc1)
 do irad=2,nc1
  shapefunc(irad)=(sin(pi*rad(irad)/rc1)/(pi*rad(irad)/rc1))**2
 end do
 do irad=nc1+1,mesh
  shapefunc(irad)=0.d0
 end do

 allocate(ff(mesh))
 ff(1:mesh)=4.d0*pi*shapefunc(1:mesh)*rad(1:mesh)**3
 call ctrap(mesh,ff,log(amesh),intg)
 gnorm =1.d0/intg
 deallocate(ff)

!print*,'NZC'
 do irad=1,mesh
  nzc(irad)=4.d0*pi*nc(irad)*rad(irad)**2.d0-&
&  4.d0*pi*shapefunc(irad)*rad(irad)**2*znucl*gnorm
! nzc(irad)=4.d0*pi*nc(irad)*rad(irad)**2.d0
! print*,rad(irad),nzc(irad),shapefunc(irad)*znucl*gnorm
 end do

 allocate(rvhn(mesh))
 rvhn(1)=zero

 allocate (den1(mesh),den2(mesh),&
& den3(mesh),den4(mesh))

 den1(1)=zero;den2(1)=zero
 den3(1)=zero;den4(1)=zero

 do ir=2,mesh
  den1(ir)= rad(ir)*nzc(ir)
  den2(ir)= den1(ir)/rad(ir)
 end do

 step=log(amesh)
!for first few points do stupid integral
 den3(1) = zero
 den4(1) = zero
 do ir=2,mesh
  call ctrap(ir,den1(1:ir),step,den3(ir))
  call ctrap(ir,den2(1:ir),step,den4(ir))
 end do

 do ir=1,mesh
  rvhn(ir)=den3(ir)+rad(ir)*(den4(mesh)-den4(ir))
 end do


!print*, 'VHTNZC'
 do ir=2,mesh
  vhtnzc(ir)=rvhn(ir)/rad(ir)
 end do
 yp2=(vhtnzc(3)-vhtnzc(2))/(rad(3)-rad(2))
 yp3=(vhtnzc(4)-vhtnzc(3))/(rad(4)-rad(3))
 yp1=yp2+(yp2-yp3)*rad(2)/(rad(3)-rad(2))
 vhtnzc(1)=vhtnzc(2)-(yp1+yp2)*rad(2)

 deallocate(nzc,rvhn)


end subroutine calc_vhtnzc
!!***
