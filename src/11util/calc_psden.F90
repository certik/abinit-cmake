!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_psden
!! NAME
!! calc_psden
!!
!! FUNCTION
!! Calculate a pseudo-density from an original density on a radial grid
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
!!  rc= cut-off radius
!!  step= step of the linear radial grid
!!
!! OUTPUT
!!  ff= pseudized density
!!
!! NOTES
!! ff=bb*sin(aa*r)/r
!!
!! PARENTS
!!      psp6in
!!
!! CHILDREN
!!      ctrap
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_psden(ff,mesh,nc,rc,step)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util, except_this_one => calc_psden
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mesh
 real(dp),intent(in) :: rc,step
!arrays
 real(dp),intent(in) :: nc(mesh)
 real(dp),intent(out) :: ff(mesh)

!Local variables-------------------------------
!scalars
 integer :: ii,nc1
 real(dp) :: aa,aa1,aa2,bb,c1,c3,cc,f0,f0p,norm1,norm2,rc1,rj

! *************************************************************************

 rc1=rc/5.d0
 nc1=int(rc1/step)+1
 rc1=(nc1-1)*step
 do ii=1,nc1
  ff(ii)=four_pi*nc(ii)*(step*dble(ii-1))**2
 end do
 call ctrap(nc1,ff(1:nc1),step,c3)
 f0=nc(nc1)
 f0p=(nc(nc1+1)-nc(nc1-1))/(2.d0*step)
 c1=-log(f0)

 aa1=zero
 aa=c1-aa1*rc1**4+rc1*(f0p/f0+4.d0*aa1*rc1**3)/2.d0
 bb=-(f0p/f0+4*aa1*rc1**3)/(2.d0*rc1)

 rj=step
 ff(1)=zero
 do ii=2,nc1
  ff(ii) = four_pi*rj**2*exp(-aa-bb*rj**2-aa1*rj**4)
  rj=rj+step
 end do
 call ctrap (nc1,ff(1:nc1),step,norm1)

 aa2=-10.4d0
 aa=c1-aa2*rc1**4+rc1*(f0p/f0+4.d0*aa2*rc1**3)/2.d0
 bb=-(f0p/f0+4*aa2*rc1**3)/(2.d0*rc1)


 rj=step
 ff(1)=zero
 do ii=2,nc1
! if((aa+bb*rj**2+aa1*rj**4)>
  ff(ii) =four_pi*rj**2*exp(-aa-bb*rj**2-aa2*rj**4)
  rj=rj+step
 end do
 call ctrap (nc1,ff(1:nc1),step,norm2)

 print*,c3,norm1,norm2
 if ((norm1-c3)*(norm2-c3)>0.d0) then
  print *, "Big pb in calc_psdens"
  stop
 end if

 do while (abs(norm2-c3)>tol6)
  cc=(aa1+aa2)/2.d0
  aa=c1-cc*rc1**4+rc1*(f0p/f0+4.d0*cc*rc1**3)/2.d0
  bb=-(f0p/f0+4*cc*rc1**3)/(2.d0*rc1)


  rj=step
  ff(1)=zero
  do ii=2,nc1
   ff(ii) =four_pi*rj**2*exp(-aa-bb*rj**2-cc*rj**4)
   rj=rj+step
  end do
  call ctrap (nc1,ff(1:nc1),step,norm2)


  if ((norm1-c3)*(norm2-c3)>0.d0) then
   aa1=cc
   norm1=norm2
  else
   aa2=cc
  end if
 end do

 call ctrap(nc1,ff(1:mesh),step,norm1)
 print*, norm1
 rj=step
 ff(1)=exp(-aa)
 do ii=2,nc1
  ff(ii)=ff(ii)/four_pi/rj**2
  rj=rj+step
! write (6,*) rj,ff(ii)
 end do
 do ii=nc1+1,mesh
  ff(ii) = nc(ii)
  rj=rj+step
! write (6,*) rj,ff(ii)
 end do

 print*,'aa=',aa
 print*,'bb=',bb
 print*,'cc=',cc

!stop

end subroutine calc_psden
!!***
