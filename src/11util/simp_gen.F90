!{\src2tex{textfont=tt}}
!!****f* ABINIT/simp_gen
!! NAME
!! simp_gen
!!
!! FUNCTION
!! Do integral on a given (generalized) grid using Simpson rule
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!  func(:)=integrand values
!!
!! OUTPUT
!!  intg=resulting integral by Simpson rule
!!
!! PARENTS
!!      optics_paw,pawdenpot,pawdij,pawinit,pawpuxinit,pawxc,pawxcm,poisson
!!      psp7cg,psp7in,psp7lo,psp7nl
!!
!! CHILDREN
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine simp_gen(intg,func,radmesh)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: intg
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: func(radmesh%int_meshsz)

!Local variables-------------------------------
!scalars
 integer :: ii,nn
 real(dp) :: resid,simp

! *************************************************************************

 nn=radmesh%int_meshsz

 simp=zero
 do ii=1,nn
  simp=simp+func(ii)*radmesh%simfact(ii)
 end do

 resid=zero
 if (radmesh%mesh_type==3) then
  resid=half*(func(2)+func(1))*(radmesh%rad(2)-radmesh%rad(1))
  if (mod(nn,2)==1) resid=resid+radmesh%stepint/3.d0*(1.25d0*func(2)*radmesh%radfact(2) &
&  +2.d0*func(3)*radmesh%radfact(3)-0.25d0*func(4)*radmesh%radfact(4))
 else if (mod(nn,2)==0) then
  resid=radmesh%stepint/3.d0*(1.25d0*func(1)*radmesh%radfact(1)+2.d0*func(2)*radmesh%radfact(2) &
&  -0.25d0*func(3)*radmesh%radfact(3))
 end if

 intg=simp+resid

end subroutine simp_gen
!!***
