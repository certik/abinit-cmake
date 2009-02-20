!{\src2tex{textfont=tt}}
!!****f* ABINIT/ifromr
!! NAME
!! ifromr
!!
!! FUNCTION
!! Retreive Index FROM a given R value in a radial grid
!! Grid can be regular or logarithimc
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  rr=given input r value
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! OUTPUT
!!  ifromr=index of rr in radial grid
!!
!! PARENTS
!!      compmesh,psp7in
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

function ifromr(radmesh,rr)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: ifromr
 real(dp),intent(in) :: rr
 type(pawrad_type),intent(in) :: radmesh

!Local variables-------------------------------

! *************************************************************************

 if (radmesh%mesh_type==1) then
  ifromr=int(tol8+rr/radmesh%rstep)+1
 else if (radmesh%mesh_type==2) then
  ifromr=int(tol8+log(1.d0+rr/radmesh%rstep)/radmesh%lstep)+1
 else if (radmesh%mesh_type==3) then
  if (rr<radmesh%rstep) then
   ifromr=1
  else
   ifromr=int(tol8+log(rr/radmesh%rstep)/radmesh%lstep)+2
  end if
 else if (radmesh%mesh_type==4) then
  ifromr=int(tol8+(1.d0-exp(-rr/radmesh%rstep))/radmesh%lstep)+1
 end if

 end function ifromr
!!***
