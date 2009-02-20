!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_prim_internals
!! NAME
!! make_prim_internals
!!
!! FUNCTION
!!  determine the bonds, angles and dihedrals for a starting
!!  geometry, based on covalent radii for the atoms.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables for this dataset
!! icenter= index of the center of the number of shifts
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! angs= number of angles
!! bonds(2,2,nbond)=for a bond between iatom and jatom
!!              bonds(1,1,nbond) = iatom
!!              bonds(2,1,nbond) = icenter
!!              bonds(1,2,nbond) = jatom
!!              bonds(2,2,nbond) = irshift
!! carts(2,ncart)= index of total primitive internal, and atom (carts(2,:))
!! dihedrals(2,4,ndihed)=indexes to characterize dihedrals
!! nang(2,3,nang)=indexes to characterize angles
!! nbond=number of bonds
!! ncart=number of cartesian coordinates included (for constraints)
!! ndihed= number of dihedrals
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!!   adds cartesian coordinates if the number of internals with a given atom is < 4
!!   the chosen coordinate could be optimized to be less dependent of the internals
!!   already incorporated.
!!
!! PARENTS
!!      delocint
!!
!! CHILDREN
!!      make_angles,make_bonds,make_dihedrals
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine make_prim_internals(angs,bonds,carts,dihedrals,icenter,nbond,nang,ndihed,ncart,&
&   dtset,nrshift,rprimd,rshift,xcart)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => make_prim_internals
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter,nrshift
 integer,intent(out) :: nang,nbond,ncart,ndihed
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,pointer :: angs(:,:,:),bonds(:,:,:),carts(:,:),dihedrals(:,:,:)
 real(dp),intent(in) :: rprimd(3,3),rshift(3,nrshift),xcart(3,dtset%natom)

!Local variables ------------------------------
! function
!scalars
 integer :: iang,iatom,ibond,icart,idihed,ii
 real(dp) :: sp
!arrays
 integer :: particip_atom(dtset%natom)
 integer,allocatable :: badangles(:)
 real(dp) :: rpt1(3),rpt2(3),rpt3(3)
!no_abirules

!************************************************************************

!DEBUG
!write (*,*) 'make_deloc_internals: enter'
!ENDDEBUG
 particip_atom(:) = 0

 call make_bonds(bonds,nbond,dtset,icenter,nrshift,rprimd,rshift,xcart)

!DEBUG
!write(*,*) size(bonds,1)
!write(*,*) size(bonds,2)
!write(*,*) size(bonds,3)
!ENDDEBUG

 do ibond=1,nbond
  write (*,'(a,i4,2(2i5,2x))') 'bond ', ibond, bonds(:,:,ibond)
  particip_atom(bonds(1,1,ibond)) = particip_atom(bonds(1,1,ibond))+1
  particip_atom(bonds(1,2,ibond)) = particip_atom(bonds(1,2,ibond))+1
 end do

!DEBUG
!return
!ENDDEBUG

 call make_angles(angs,bonds,icenter,nang,nbond,dtset)
 allocate(badangles(nang))
 badangles(:) = 0
 do iang=1,nang
  write (*,'(a,i4,3(2i5,2x))') 'angle ', iang, angs(:,:,iang)
  particip_atom(angs(1,1,iang)) = particip_atom(angs(1,1,iang))+1
  particip_atom(angs(1,2,iang)) = particip_atom(angs(1,2,iang))+1
  particip_atom(angs(1,3,iang)) = particip_atom(angs(1,3,iang))+1

! DEBUG
! rpt1(:) = xcart(:,angs(1,1,iang)) &
! & + rshift(1,angs(2,1,iang))*rprimd(:,1) &
! & + rshift(2,angs(2,1,iang))*rprimd(:,2) &
! & + rshift(3,angs(2,1,iang))*rprimd(:,3)
! rpt2(:) = xcart(:,angs(1,2,iang)) &
! & + rshift(1,angs(2,2,iang))*rprimd(:,1) &
! & + rshift(2,angs(2,2,iang))*rprimd(:,2) &
! & + rshift(3,angs(2,2,iang))*rprimd(:,3)
! rpt3(:) = xcart(:,angs(1,3,iang)) &
! & + rshift(1,angs(2,3,iang))*rprimd(:,1) &
! & + rshift(2,angs(2,3,iang))*rprimd(:,2) &
! & + rshift(3,angs(2,3,iang))*rprimd(:,3)
! write (*,*) rpt1,rpt2,rpt3,bond_length(rpt1,rpt2),bond_length(rpt2,rpt3)
! ENDDEBUG

! check if angles are 180 degrees: discard the dihedrals in that case.
  rpt1(:) = xcart(:,angs(1,1,iang)) &
&  + rshift(1,angs(2,1,iang))*rprimd(:,1) &
&  + rshift(2,angs(2,1,iang))*rprimd(:,2) &
&  + rshift(3,angs(2,1,iang))*rprimd(:,3) &
&  - xcart(:,angs(1,2,iang)) &
&  - rshift(1,angs(2,2,iang))*rprimd(:,1) &
&  - rshift(2,angs(2,2,iang))*rprimd(:,2) &
&  - rshift(3,angs(2,2,iang))*rprimd(:,3)

  rpt3(:) = xcart(:,angs(1,3,iang)) &
&  + rshift(1,angs(2,3,iang))*rprimd(:,1) &
&  + rshift(2,angs(2,3,iang))*rprimd(:,2) &
&  + rshift(3,angs(2,3,iang))*rprimd(:,3) &
&  - xcart(:,angs(1,2,iang)) &
&  - rshift(1,angs(2,2,iang))*rprimd(:,1) &
&  - rshift(2,angs(2,2,iang))*rprimd(:,2) &
&  - rshift(3,angs(2,2,iang))*rprimd(:,3)
  sp = (rpt1(1)*rpt3(1)+rpt1(2)*rpt3(2)+rpt1(3)*rpt3(3))&
&  / sqrt(rpt1(1)*rpt1(1)+rpt1(2)*rpt1(2)+rpt1(3)*rpt1(3)) &
&  / sqrt(rpt3(1)*rpt3(1)+rpt3(2)*rpt3(2)+rpt3(3)*rpt3(3))
  if (abs(abs(sp) - one) < tol6) then
   write (*,*) 'make_prim_internals : an angle is too close to 180 degrees:'
   write (*,*) '   will discard dihedrals using it '
   badangles(iang) = 1
  end if
 end do


 call make_dihedrals(angs,badangles,bonds,dihedrals,icenter,nbond,nang,ndihed,nrshift,dtset)
 do idihed=1,ndihed
  write (*,'(a,i4,4(2i5,2x))') 'dihedral ', idihed, dihedrals(:,:,idihed)
  particip_atom(dihedrals(1,1,idihed)) = particip_atom(dihedrals(1,1,idihed))+1
  particip_atom(dihedrals(1,2,idihed)) = particip_atom(dihedrals(1,2,idihed))+1
  particip_atom(dihedrals(1,3,idihed)) = particip_atom(dihedrals(1,3,idihed))+1
  particip_atom(dihedrals(1,4,idihed)) = particip_atom(dihedrals(1,4,idihed))+1

! DEBUG
! do ii=1,4
! write (*,'((3E16.6,2x))') xcart(:,dihedrals(1,ii,idihed)) + &
! &  rshift(1,dihedrals(2,ii,idihed))*rprimd(:,1)   + &
! &  rshift(2,dihedrals(2,ii,idihed))*rprimd(:,2)   + &
! &  rshift(2,dihedrals(2,ii,idihed))*rprimd(:,3)
! end do
! ENDDEBUG

 end do

 write (*,*) 'make_deloc_internals: nbond,nang,ndihed = ', nbond,nang,ndihed

!Check all atoms participate in at least 4 primitives. Otherwise, we should
!probably add cartesian coordinates to the internal ones.
 ncart = 0
 do iatom=1,dtset%natom
  if (particip_atom(iatom) < 4) then
   write (*,*) ' make_prim_internals : Warning : atom ', iatom, &
&   ' does not belong to enough primitives to determine its'
   write (*,*) ' position uniquely ! instead : ', particip_atom(iatom)
   write (*,*) ' Will add cartesian coordinates to set of internals.'
!  write (*,*) ' Not done yet.'
!  stop
   ncart = ncart + 4-particip_atom(iatom)
  end if
 end do
 allocate (carts (2,ncart))
 icart = 0
 do iatom=1,dtset%natom
  if (particip_atom(iatom) < 4) then
!  
!  kind of arbitrary : include first few directions for the atom:
!  x, then y then z
!  
   do ii=1,4-particip_atom(iatom)
    icart = icart+1
    carts(1,icart) = ii
    carts(2,icart) = iatom
   end do
  end if
 end do


end subroutine make_prim_internals

!-----------------------------------------------------------------------------
!!***
