!{\src2tex{textfont=tt}}
!!****f* ABINIT/copymesh
!! NAME
!! copymesh
!!
!! FUNCTION
!! Copy one radial mesh (in a generalized format) to another
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mesh1 <type(pawrad_type)>=data containing radial grid information of input mesh
!!
!! OUTPUT
!!  mesh2 <type(pawrad_type)>=data containing radial grid information of output mesh
!!
!! PARENTS
!!      paw_mkrhox_spl,psp7in,psp7nl
!!
!! CHILDREN
!!
!! NOTES
!!  Possible mesh types (mesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine copymesh(mesh1,mesh2)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(pawrad_type),intent(in) :: mesh1
 type(pawrad_type),intent(out) :: mesh2

!Local variables-------------------------------
!scalars
 integer :: ierr,ir

! *************************************************************************
 mesh2%mesh_type =mesh1%mesh_type
 mesh2%mesh_size =mesh1%mesh_size
 mesh2%int_meshsz=mesh1%int_meshsz
 mesh2%lstep     =mesh1%lstep
 mesh2%rstep     =mesh1%rstep
 mesh2%stepint   =mesh1%stepint
 mesh2%rmax      =mesh1%rmax
!those lines are not useful and the behavior is undefined. => problems with g95 (PMA)
!if (associated(mesh2%rad)) deallocate(mesh2%rad,STAT=ierr)
!if (associated(mesh2%radfact)) deallocate(mesh2%radfact,STAT=ierr)
!if (associated(mesh2%simfact)) deallocate(mesh2%simfact,STAT=ierr)
 allocate(mesh2%rad(mesh1%mesh_size))
 allocate(mesh2%radfact(mesh1%mesh_size))
 allocate(mesh2%simfact(mesh1%mesh_size))
 do ir=1,mesh1%mesh_size
  mesh2%rad(ir)    =mesh1%rad(ir)
  mesh2%radfact(ir)=mesh1%radfact(ir)
  mesh2%simfact(ir)=mesh1%simfact(ir)
 end do
end subroutine copymesh
!!***
