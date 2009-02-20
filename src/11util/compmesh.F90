!{\src2tex{textfont=tt}}
!!****f* ABINIT/compmesh
!! NAME
!! compmesh
!!
!! FUNCTION
!! Compute all points of a radial mesh
!! Grid can be regular or logarithimc
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SIDE EFFECTS
!!  mesh <type(pawrad_type)>=data containing radial grid information
!!
!! PARENTS
!!      paw_mkrhox_spl,psp7in
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
!! INPUTS
!! mesh=data containing radial grid information
!! r_for_intg=upper boundary for (future) integration over the radial grid
!!            (can be negative for an integration over the whole grid)
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine compmesh(mesh,r_for_intg)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util, except_this_one => compmesh
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: r_for_intg
 type(pawrad_type),intent(inout) :: mesh

!Local variables-------------------------------
!scalars
 integer :: ir,ir_last,isim
 real(dp) :: hh

!************************************************************************
 if (mesh%mesh_type==1) then
  isim=3
  mesh%stepint=mesh%rstep
  mesh%rad(1)=zero;mesh%radfact(1)=one
  do ir=2,mesh%mesh_size
   mesh%rad(ir) =mesh%rstep*dble(ir-1)
   mesh%radfact(ir)=one
  end do
 else if (mesh%mesh_type==2) then
  isim=3
  mesh%stepint=mesh%lstep
  mesh%rad(1)=zero;mesh%radfact(1)=mesh%rstep
  do ir=2,mesh%mesh_size
   mesh%rad(ir) =mesh%rstep*(exp(mesh%lstep*dble(ir-1))-one)
   mesh%radfact(ir)=mesh%rad(ir)+mesh%rstep
  end do
 else if (mesh%mesh_type==3) then
  isim=4
  mesh%stepint=mesh%lstep
  mesh%rad(1)=zero;mesh%radfact(1)=zero
  do ir=2,mesh%mesh_size
   mesh%rad(ir) =mesh%rstep*exp(mesh%lstep*dble(ir-2))
   mesh%radfact(ir)=mesh%rad(ir)
  end do
 else if (mesh%mesh_type==4) then
  isim=3
  mesh%lstep=one/dble(mesh%mesh_size)
  mesh%stepint=mesh%lstep
  mesh%rad(1)=zero;mesh%radfact(1)=mesh%rstep
  do ir=2,mesh%mesh_size
   mesh%rad(ir) =-mesh%rstep*log(one-mesh%lstep*dble(ir-1))
   mesh%radfact(ir)=mesh%rstep/(one-mesh%lstep*dble(ir-1))
  end do
 end if

 mesh%int_meshsz=mesh%mesh_size
 if (r_for_intg>0.d0) then
  ir=min(ifromr(mesh,r_for_intg),mesh%mesh_size)
  if (ir<mesh%mesh_size) then
   if (abs(mesh%rad(ir+1)-r_for_intg)<abs(mesh%rad(ir)-r_for_intg)) ir=ir+1
  end if
  if (ir>1) then
   if (abs(mesh%rad(ir-1)-r_for_intg)<abs(mesh%rad(ir)-r_for_intg)) ir=ir-1
  end if
  mesh%int_meshsz=ir
 end if

 hh=mesh%stepint/3.d0
 mesh%simfact(mesh%int_meshsz)=hh*mesh%radfact(mesh%int_meshsz)
 mesh%simfact(1:isim-2)=zero
 ir_last=1
 do ir=mesh%int_meshsz,isim,-2
  mesh%simfact(ir-1)=4.d0*hh*mesh%radfact(ir-1)
  mesh%simfact(ir-2)=2.d0*hh*mesh%radfact(ir-2)
  ir_last=ir-2
 end do
 mesh%simfact(ir_last)=half*mesh%simfact(ir_last)
 if (mesh%int_meshsz<mesh%mesh_size) mesh%simfact(mesh%int_meshsz+1:mesh%mesh_size)=zero

 mesh%rmax=mesh%rad(mesh%mesh_size)

end subroutine compmesh
!!***
