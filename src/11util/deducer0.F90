!{\src2tex{textfont=tt}}
!!****f* ABINIT/deducer0
!! NAME
!! deducer0
!!
!! FUNCTION
!! Extrapolate r=0 value of a function from values near r=0
!! using a 3 points formula
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  funcsz=size of array func
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! SIDE EFFECTS
!!  func(funcsz)=array containing values of function to extrapolate
!!
!! PARENTS
!!      nderiv_gen,nhatgrid,pawdens,pawxcgrad,poisson,psp7in
!!
!! CHILDREN
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine deducer0(func,funcsz,radmesh)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: funcsz
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp) :: func(funcsz)

!Local variables-------------------------------

! *************************************************************************

 if (radmesh%mesh_type==1.or.radmesh%mesh_type==2.or.radmesh%mesh_type==4) then
  func(1)=func(4)+3*(func(2)-func(3))
 else if (radmesh%mesh_type==3) then
  func(1)=func(4)+exp(two*radmesh%lstep)/(exp(radmesh%lstep)-one)*(func(2)-func(3))
 end if

end subroutine deducer0
!!***
