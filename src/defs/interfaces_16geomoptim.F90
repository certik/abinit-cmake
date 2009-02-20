!!****m* ABINIT/interfaces_16geomoptim
!! NAME
!! interfaces_16geomoptim
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/16geomoptim
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_16geomoptim

 implicit none

interface
 subroutine calc_b_matrix(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&  
  &  dtset,nrshift,rprimd,rshift,xcart,b_matrix)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nang
  integer,intent(in) :: nbond
  integer,intent(in) :: ncart
  integer,intent(in) :: ndihed
  integer,intent(in) :: ninternal
  integer,intent(in) :: nrshift
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: angs(:,:,:)
  integer,pointer :: bonds(:,:,:)
  integer,pointer :: carts(:,:)
  integer,pointer :: dihedrals(:,:,:)
  real(dp),intent(out) :: b_matrix(ninternal,3*dtset%natom)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rshift(3,nrshift)
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine calc_b_matrix
end interface

interface
 subroutine dbond_length_d1(r1,r2,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
 end subroutine dbond_length_d1
end interface

interface
 subroutine dang_d1(r1,r2,r3,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
 end subroutine dang_d1
end interface

interface
 subroutine dang_d2(r1,r2,r3,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
 end subroutine dang_d2
end interface

interface
 subroutine ddihedral_d1(r1,r2,r3,r4,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
  real(dp),intent(in) :: r4(3)
 end subroutine ddihedral_d1
end interface

interface
 subroutine ddihedral_d2(r1,r2,r3,r4,bb)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bb(3)
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
  real(dp),intent(in) :: r4(3)
 end subroutine ddihedral_d2
end interface

interface
 subroutine calc_prim_int(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&  
  &  dtset,nrshift,rprimd,rshift,xcart,prim_int)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nang
  integer,intent(in) :: nbond
  integer,intent(in) :: ncart
  integer,intent(in) :: ndihed
  integer,intent(in) :: ninternal
  integer,intent(in) :: nrshift
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: angs(:,:,:)
  integer,pointer :: bonds(:,:,:)
  integer,pointer :: carts(:,:)
  integer,pointer :: dihedrals(:,:,:)
  real(dp),intent(out) :: prim_int(ninternal)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rshift(3,nrshift)
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine calc_prim_int
end interface

interface
 function bond_length(r1,r2)
  use defs_basis
  implicit none
  real(dp) :: bond_length
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
 end function bond_length
end interface

interface
 function angle_ang(r1,r2,r3)
  use defs_basis
  implicit none
  real(dp) :: angle_ang
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
 end function angle_ang
end interface

interface
 function angle_dihedral(r1,r2,r3,r4)
  use defs_basis
  implicit none
  real(dp) :: angle_dihedral
  real(dp),intent(in) :: r1(3)
  real(dp),intent(in) :: r2(3)
  real(dp),intent(in) :: r3(3)
  real(dp),intent(in) :: r4(3)
 end function angle_dihedral
end interface

interface
 subroutine deloc2xcart(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,nrshift,&  
  &  dtset,rprimd,rshift,xcart,&  
  &  deloc_int,btinv,u_matrix)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nang
  integer,intent(in) :: nbond
  integer,intent(in) :: ncart
  integer,intent(in) :: ndihed
  integer,intent(in) :: ninternal
  integer,intent(in) :: nrshift
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: angs(:,:,:)
  integer,pointer :: bonds(:,:,:)
  integer,pointer :: carts(:,:)
  integer,pointer :: dihedrals(:,:,:)
  real(dp),intent(out) :: btinv(3*(dtset%natom-1),3*dtset%natom)
  real(dp),intent(in) :: deloc_int(3*(dtset%natom-1))
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rshift(3,nrshift)
  real(dp),intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
  real(dp),intent(inout) :: xcart(3,dtset%natom)
 end subroutine deloc2xcart
end interface

interface
 function get_mix(x1,y1,x2,y2,x3,y3)
  use defs_basis
  implicit none
  real(dp) :: get_mix
  real(dp),intent(in) :: x1
  real(dp),intent(in) :: x2
  real(dp),intent(in) :: x3
  real(dp),intent(in) :: y1
  real(dp),intent(in) :: y2
  real(dp),intent(in) :: y3
 end function get_mix
end interface

interface
 subroutine fred2fdeloc(btinv,deloc_force,fred,dtset,rprimd)
  use defs_basis
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: btinv(3*(dtset%natom-1),3*dtset%natom)
  real(dp),intent(out) :: deloc_force(3*(dtset%natom-1))
  real(dp),intent(out) :: fred(3,dtset%natom)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine fred2fdeloc
end interface

interface
 subroutine hessinit(dtfil, dtset, hessin, init_matrix, ndim, ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ndim
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: hessin(ndim,ndim)
  real(dp),intent(in) :: init_matrix(3,3)
 end subroutine hessinit
end interface

interface
 subroutine hessupdt(hessin,iatfix,natom,ndim,vin,vin_prev,vout,vout_prev)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  real(dp),intent(inout) :: hessin(ndim,ndim)
  integer,intent(in) :: iatfix(3,natom)
  real(dp),intent(in) :: vin(ndim)
  real(dp),intent(in) :: vin_prev(ndim)
  real(dp),intent(in) :: vout(ndim)
  real(dp),intent(in) :: vout_prev(ndim)
 end subroutine hessupdt
end interface

interface
 subroutine make_angles(angs,bonds,icenter,nang,nbond,dtset)
  use defs_datatypes
  implicit none
  integer,intent(in) :: icenter
  integer,intent(out) :: nang
  integer,intent(in) :: nbond
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: angs(:,:,:)
  integer,pointer :: bonds(:,:,:)
 end subroutine make_angles
end interface

interface
 subroutine make_bonds(bonds,nbond,dtset,icenter,nrshift,rprimd,rshift,xcart)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: icenter
  integer,intent(out) :: nbond
  integer,intent(in) :: nrshift
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: bonds(:,:,:)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rshift(3,nrshift)
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine make_bonds
end interface

interface
 subroutine make_dihedrals(angs,badangles,bonds,dihedrals,icenter,nbond,nang,ndihed,nrshift,dtset)
  use defs_datatypes
  implicit none
  integer,intent(in) :: icenter
  integer,intent(in) :: nang
  integer,intent(in) :: nbond
  integer,intent(out) :: ndihed
  integer,intent(in) :: nrshift
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: angs(:,:,:)
  integer,pointer :: bonds(:,:,:)
  integer,pointer :: dihedrals(:,:,:)
  integer,intent(in) :: badangles(nang)
 end subroutine make_dihedrals
end interface

interface
 subroutine make_prim_internals(angs,bonds,carts,dihedrals,icenter,nbond,nang,ndihed,ncart,&  
  &  dtset,nrshift,rprimd,rshift,xcart)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: icenter
  integer,intent(out) :: nang
  integer,intent(out) :: nbond
  integer,intent(out) :: ncart
  integer,intent(out) :: ndihed
  integer,intent(in) :: nrshift
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: angs(:,:,:)
  integer,pointer :: bonds(:,:,:)
  integer,pointer :: carts(:,:)
  integer,pointer :: dihedrals(:,:,:)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rshift(3,nrshift)
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine make_prim_internals
end interface

interface
 subroutine xcart2deloc(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&  
  &  dtset,nrshift,rprimd,rshift,xcart,&  
  &  bt_inv_matrix,u_matrix,deloc_int,prim_int)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nang
  integer,intent(in) :: nbond
  integer,intent(in) :: ncart
  integer,intent(in) :: ndihed
  integer,intent(in) :: ninternal
  integer,intent(in) :: nrshift
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: angs(:,:,:)
  integer,pointer :: bonds(:,:,:)
  integer,pointer :: carts(:,:)
  integer,pointer :: dihedrals(:,:,:)
  real(dp),intent(out) :: bt_inv_matrix(3*(dtset%natom-1),3*dtset%natom)
  real(dp),intent(out) :: deloc_int(3*(dtset%natom-1))
  real(dp),intent(out) :: prim_int(ninternal)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rshift(3,nrshift)
  real(dp),intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine xcart2deloc
end interface

interface
 subroutine calc_btinv_matrix(b_matrix,dtset,nang,nbond,ndihed,ncart,ninternal,&  
  &  bt_inv_matrix,u_matrix)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nang
  integer,intent(in) :: nbond
  integer,intent(in) :: ncart
  integer,intent(in) :: ndihed
  integer,intent(in) :: ninternal
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: b_matrix(ninternal,3*dtset%natom)
  real(dp),intent(out) :: bt_inv_matrix(3*(dtset%natom-1),3*dtset%natom)
  real(dp),intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
 end subroutine calc_btinv_matrix
end interface

interface
 subroutine align_u_matrices(dtset,ninternal,u_matrix,u_matrix_old,s_matrix,f_eigs)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: ninternal
  type(dataset_type), intent(in) :: dtset
  real(dp), intent(inout) :: f_eigs(3*dtset%natom)
  real(dp), intent(inout) :: s_matrix(3*dtset%natom,3*dtset%natom)
  real(dp), intent(inout) :: u_matrix(ninternal,3*(dtset%natom-1))
  real(dp), intent(in) :: u_matrix_old(ninternal,3*(dtset%natom-1))
 end subroutine align_u_matrices
end interface

interface
 subroutine cross_product(r1,r2,cp)
  use defs_basis
  implicit none
  real(dp), intent(out) :: cp(3)
  real(dp), intent(in) :: r1(3)
  real(dp), intent(in) :: r2(3)
 end subroutine cross_product
end interface

interface
 subroutine xcart2deloc_fixb(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&  
  &  dtset,nrshift,rprimd,rshift,xcart,&  
  &  u_matrix,deloc_int,prim_int)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nang
  integer,intent(in) :: nbond
  integer,intent(in) :: ncart
  integer,intent(in) :: ndihed
  integer,intent(in) :: ninternal
  integer,intent(in) :: nrshift
  type(dataset_type),intent(in) :: dtset
  integer,pointer :: angs(:,:,:)
  integer,pointer :: bonds(:,:,:)
  integer,pointer :: carts(:,:)
  integer,pointer :: dihedrals(:,:,:)
  real(dp),intent(out) :: deloc_int(3*(dtset%natom-1))
  real(dp),intent(out) :: prim_int(ninternal)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: rshift(3,nrshift)
  real(dp),intent(in) :: u_matrix(ninternal,3*(dtset%natom-1))
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine xcart2deloc_fixb
end interface

interface
 subroutine xfpack(acell,acell0,fred,natom,ndim,nsym,optcell,option,rprim,rprimd0,&  
  &  strtarget,strten,symrel,ucvol,ucvol0,vin,vout,xred)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ndim
  integer,intent(in) :: nsym
  integer,intent(in) :: optcell
  integer,intent(in) :: option
  real(dp),intent(inout) :: ucvol
  real(dp),intent(in) :: ucvol0
  integer,intent(in) :: symrel(3,3)
  real(dp),intent(inout) :: acell(3)
  real(dp),intent(in) :: acell0(3)
  real(dp),intent(inout) :: fred(3,natom)
  real(dp),intent(inout) :: rprim(3,3)
  real(dp),intent(in) :: rprimd0(3,3)
  real(dp),intent(in) :: strtarget(6)
  real(dp),intent(inout) :: strten(6)
  real(dp),intent(inout) :: vin(ndim)
  real(dp),intent(inout) :: vout(ndim)
  real(dp),intent(inout) :: xred(3,natom)
 end subroutine xfpack
end interface

end module interfaces_16geomoptim
!!***
