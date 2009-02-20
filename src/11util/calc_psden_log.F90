!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_psden_log
!! NAME
!! calc_psden_log
!!
!! FUNCTION
!! Calculate a pseudo-density from an original density on a log grid
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GJ, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mmax= dimension of nc
!!  nc= density to be pseudized
!!  rc= cut-off radius
!!  rad(mmax) = radial mesh
!!
!! OUTPUT
!!  ff= pseudized density
!!
!! NOTES
!! ff=bb*sin(aa*r)/r
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

subroutine calc_psden_log(ff,mmax,nc,rc,rad)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util, except_this_one => calc_psden_log
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mmax
 real(dp),intent(in) :: rc
!arrays
 real(dp),intent(in) :: nc(mmax),rad(mmax)
 real(dp),intent(out) :: ff(mmax)

!Local variables-------------------------------
!scalars
 integer :: ii,nc1
 real(dp) :: aa,aa1,aa2,amesh,bb,c1,c3,cc,f0,f0p,norm1,norm2,rc1
!arrays
 real(dp) :: gg(mmax)

! *************************************************************************

 amesh=rad(2)/rad(1)

 rc1=rc/5.d0
 nc1=int(log(rc1/rad(1))/log(amesh))+1
 rc1=rad(nc1)
 do ii=1,nc1
  ff(ii)=four_pi*nc(ii)*rad(ii)**3
 end do
 call ctrap(nc1,ff(1:nc1),log(amesh),c3)
 c3=c3+ff(1)/2.d0
 f0=nc(nc1)
 f0p=(nc(nc1+1)-nc(nc1-1))/(2.d0*log(amesh))
 c1=-log(f0)

 aa1=zero
 aa=c1-aa1*rc1**4+rc1*(f0p/f0+4.d0*aa1*rc1**3)/2.d0
 bb=-(f0p/f0+4*aa1*rc1**3)/(2.d0*rc1)

 do ii=1,nc1
  ff(ii) = four_pi*rad(ii)**3*exp(-aa-bb*rad(ii)**2-aa1*rad(ii)**4)
 end do
 call ctrap (nc1,ff(1:nc1),log(amesh),norm1)
 norm1=norm1+ff(1)/2.d0
!Pu aa2=-3.0d0
 aa2=-10.0d0
 aa=c1-aa2*rc1**4+rc1*(f0p/f0+4.d0*aa2*rc1**3)/2.d0
 bb=-(f0p/f0+4*aa2*rc1**3)/(2.d0*rc1)

 ff(1)=zero
 do ii=1,nc1
  ff(ii) =four_pi*rad(ii)**3*exp(-aa-bb*rad(ii)**2-aa2*rad(ii)**4)
 end do
 call ctrap (nc1,ff(1:nc1),log(amesh),norm2)
 norm2=norm2+ff(1)/2.d0

 print*,c3,norm1,norm2
 if ((norm1-c3)*(norm2-c3)>0.d0) then
  print *, "Big pb in calc_psdens"
  stop
 end if

 do while (abs(norm2-c3)>tol8)
  print*,abs(norm2-c3),norm1,norm2,c3
  cc=(aa1+aa2)/2.d0
  aa=c1-cc*rc1**4+rc1*(f0p/f0+4.d0*cc*rc1**3)/2.d0
  bb=-(f0p/f0+4*cc*rc1**3)/(2.d0*rc1)


  ff(1)=zero
  do ii=1,nc1
   ff(ii) =four_pi*rad(ii)**3*exp(-aa-bb*rad(ii)**2-cc*rad(ii)**4)
  end do
  call ctrap (nc1,ff(1:nc1),log(amesh),norm2)
  norm2=norm2+ff(1)/2.d0

  if ((norm1-c3)*(norm2-c3)>0.d0) then
   aa1=cc
   norm1=norm2
  else
   aa2=cc
  end if
 end do

 call ctrap(nc1,ff(1:nc1),log(amesh),norm1)
 norm1=norm1+ff(1)/2.d0
 print*, norm1

 do ii=1,nc1
  ff(ii)=ff(ii)/four_pi/rad(ii)**3
! write (6,*) rad(ii),ff(ii)
 end do
 do ii=nc1+1,mmax
  ff(ii) = nc(ii)
! write (6,*) rad(ii),ff(ii)
 end do




 gg(1)=zero
 do ii=1,mmax
  gg(ii)=four_pi*rad(ii)**3*ff(ii)
 end do
 call ctrap(mmax,gg(1:mmax),log(amesh),norm1)
 norm1=norm1+gg(1)/2.d0
 print*,'Inregrale de tnc sur la grille log=',norm1,gg(1)/2.D0

 print*,'aa=',aa
 print*,'bb=',bb
 print*,'cc=',cc

!stop

end subroutine calc_psden_log
!!***
