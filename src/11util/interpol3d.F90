!{\src2tex{textfont=tt}}
!!****f* ABINIT/interpol3d
!! NAME
!! interpol3d
!!
!! FUNCTION
!! Computes the density at any point r by linear interpolation
!! inside the eight vertices of the surrounding cube
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! r(3)=point coordinate
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! grid(nr1,nr2,nr3)=grid matrix
!!
!! OUTPUT
!! denval=density value
!!
!! PARENTS
!!      lineint,mknesting,planeint,pointint,volumeint
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine interpol3d(r,nr1,nr2,nr3,denval,grid)

 use defs_basis

 implicit none

!Arguments-------------------------------------------------------------
!scalars
 integer,intent(in) :: nr1,nr2,nr3
 real(dp),intent(out) :: denval
!arrays
 real(dp),intent(in) :: grid(nr1,nr2,nr3),r(3)

!Local variables--------------------------------------------------------
!scalars
 integer :: ir1,ir2,ir3,pr1,pr2,pr3
 real(dp) :: d1,d2,d3,iri1,iri2,iri3,x1,x2,x3

! *************************************************************************

!grid density
 d1=1.0/nr1
 d2=1.0/nr2
 d3=1.0/nr3

!lower left
 ir1=int(r(1)/d1)+1
 ir2=int(r(2)/d2)+1
 ir3=int(r(3)/d3)+1

!upper right
 pr1=mod(ir1+1,nr1)
 pr2=mod(ir2+1,nr2)
 pr3=mod(ir3+1,nr3)

!weight
 x1=1.0+(r(1)-ir1*d1)/d1
 x2=1.0+(r(2)-ir2*d2)/d2
 x3=1.0+(r(3)-ir3*d3)/d3

 if(ir1.eq.0) ir1=nr1
 if(ir2.eq.0) ir2=nr2
 if(ir3.eq.0) ir3=nr3

 if(pr1.eq.0) pr1=nr1
 if(pr2.eq.0) pr2=nr2
 if(pr3.eq.0) pr3=nr3

!calculation of the density value
 denval=0.0d0
 denval=denval+grid(ir1,ir2,ir3)*(1.0-x1)*(1.0-x2)*(1.0-x3)
 denval=denval+grid(pr1,ir2,ir3)*x1*(1.0-x2)*(1.0-x3)
 denval=denval+grid(ir1,pr2,ir3)*(1.0-x1)*x2*(1.0-x3)
 denval=denval+grid(ir1,ir2,pr3)*(1.0-x1)*(1.0-x2)*x3
 denval=denval+grid(pr1,pr2,ir3)*x1*x2*(1.0-x3)
 denval=denval+grid(ir1,pr2,pr3)*(1.0-x1)*x2*x3
 denval=denval+grid(pr1,ir2,pr3)*x1*(1.0-x2)*x3
 denval=denval+grid(pr1,pr2,pr3)*x1*x2*x3

end subroutine interpol3d
!!***
