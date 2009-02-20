!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_prim_int
!! NAME
!! calc_prim_int
!!
!! FUNCTION
!!  calculate values of primitive internal coordinates as a function of
!!  cartesian ones.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! angs= number of angles
!! bonds(2,2,nbond)=for a bond between iatom and jatom
!!              bonds(1,1,nbond) = iatom
!!              bonds(2,1,nbond) = icenter
!!              bonds(1,2,nbond) = jatom
!!              bonds(2,2,nbond) = irshift
!! carts(2,ncart)= index of total primitive internal, and atom (carts(2,:))
!! dihedrals(2,4,ndihed)=indexes to characterize dihedrals
!! dtset <type(dataset_type)>=all input variables for this dataset
!! nang(2,3,nang)=indexes to characterize angles
!! nbond=number of bonds
!! ncart=number of cartesian coordinates used
!! ndihed= number of dihedrals
!! ninternal=nbond+nang+ndihed+ncart: number of internal coordinates
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! prim_int(ninternal)=values of primitive internal coordinates
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      delocint,xcart2deloc,xcart2deloc_fixb
!!
!! CHILDREN
!!      cross_product
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_prim_int(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
& dtset,nrshift,rprimd,rshift,xcart,prim_int)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => calc_prim_int
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nang,nbond,ncart,ndihed,ninternal,nrshift
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,pointer :: angs(:,:,:),bonds(:,:,:),carts(:,:),dihedrals(:,:,:)
 real(dp),intent(in) :: rprimd(3,3),rshift(3,nrshift),xcart(3,dtset%natom)
 real(dp),intent(out) :: prim_int(ninternal)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,i4,iang,ibond,icart,idihed,iprim,s1,s2,s3,s4
!arrays
 real(dp) :: r1(3),r2(3),r3(3),r4(3)

!************************************************************************

!DEBUG
!write (*,*) ' calc_prim_int : enter'
!write (*,*) shape(angs),shape(bonds),shape(dihedrals)
!do ibond=1,nbond
!do i1=1,2
!write (*,'(2I5)') bonds(:,i1,ibond)
!end do
!end do
!ENDDEBUG
 iprim=1
!first: bond values
 do ibond=1,nbond
  i1 = bonds(1,1,ibond)
  s1 = bonds(2,1,ibond)
  r1(:) = xcart(:,i1)+rshift(1,s1)*rprimd(:,1)&
&  +rshift(2,s1)*rprimd(:,2)&
&  +rshift(3,s1)*rprimd(:,3)
  i2 = bonds(1,2,ibond)
  s2 = bonds(2,2,ibond)
  r2(:) = xcart(:,i2)+rshift(1,s2)*rprimd(:,1)&
&  +rshift(2,s2)*rprimd(:,2)&
&  +rshift(3,s2)*rprimd(:,3)
  prim_int(iprim) = bond_length(r1,r2)
  iprim=iprim+1
 end do

!second: angle values (ang)
 do iang=1,nang
  i1 = angs(1,1,iang)
  s1 = angs(2,1,iang)
  r1(:) = xcart(:,i1)+rshift(1,s1)*rprimd(:,1)&
&  +rshift(2,s1)*rprimd(:,2)&
&  +rshift(3,s1)*rprimd(:,3)
  i2 = angs(1,2,iang)
  s2 = angs(2,2,iang)
  r2(:) = xcart(:,i2)+rshift(1,s2)*rprimd(:,1)&
&  +rshift(2,s2)*rprimd(:,2)&
&  +rshift(3,s2)*rprimd(:,3)
  i3 = angs(1,3,iang)
  s3 = angs(2,3,iang)
  r3(:) = xcart(:,i3)+rshift(1,s3)*rprimd(:,1)&
&  +rshift(2,s3)*rprimd(:,2)&
&  +rshift(3,s3)*rprimd(:,3)
  prim_int(iprim) = angle_ang(r1,r2,r3)
  iprim=iprim+1
 end do

!third: dihedral values
 do idihed=1,ndihed
  i1 = dihedrals(1,1,idihed)
  s1 = dihedrals(2,1,idihed)
  r1(:) = xcart(:,i1)+rshift(1,s1)*rprimd(:,1)&
&  +rshift(2,s1)*rprimd(:,2)&
&  +rshift(3,s1)*rprimd(:,3)
  i2 = dihedrals(1,2,idihed)
  s2 = dihedrals(2,2,idihed)
  r2(:) = xcart(:,i2)+rshift(1,s2)*rprimd(:,1)&
&  +rshift(2,s2)*rprimd(:,2)&
&  +rshift(3,s2)*rprimd(:,3)
  i3 = dihedrals(1,3,idihed)
  s3 = dihedrals(2,3,idihed)
  r3(:) = xcart(:,i3)+rshift(1,s3)*rprimd(:,1)&
&  +rshift(2,s3)*rprimd(:,2)&
&  +rshift(3,s3)*rprimd(:,3)
  i4 = dihedrals(1,4,idihed)
  s4 = dihedrals(2,4,idihed)
  r4(:) = xcart(:,i4)+rshift(1,s4)*rprimd(:,1)&
&  +rshift(2,s4)*rprimd(:,2)&
&  +rshift(3,s4)*rprimd(:,3)
  prim_int(iprim) = angle_dihedral(r1,r2,r3,r4)
  iprim=iprim+1
 end do

 do icart=1,ncart
  prim_int(iprim) = xcart(carts(1,icart),carts(2,icart))
  iprim=iprim+1
 end do

!DEBUG
!write (*,*) 'Primitive internal coordinate values:'
!do iprim=1,ninternal
!if (iprim <= nbond) then
!write (*,*) iprim, prim_int(iprim)
!else
!write (*,*) iprim, prim_int(iprim), prim_int(iprim)/pi*180.0_dp
!end if
!end do
!ENDDEBUG

end subroutine calc_prim_int
!!***


!!****if* ABINIT/bond_length
!! NAME
!! bond_length
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 function bond_length(r1,r2)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: bond_length
!arrays
 real(dp),intent(in) :: r1(3),r2(3)

!Local variables ------------------------------------
!arrays
 real(dp) :: rpt(3)

!******************************************************************
 rpt(:) = r1(:)-r2(:)
 bond_length = sqrt(rpt(1)**2+rpt(2)**2+rpt(3)**2)

end function bond_length
!!***

!!****if* ABINIT/angle_ang
!! NAME
!! angle_ang
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 function angle_ang(r1,r2,r3)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => angle_ang
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: angle_ang
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3)

!Local variables ------------------------------
!scalars
 real(dp) :: cos_ang,n1,n2
!arrays
 real(dp) :: rpt12(3),rpt32(3)

!******************************************************************
 n1=bond_length(r1,r2)
 n2=bond_length(r3,r2)

 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)

 cos_ang = (rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3))/n1/n2

 if (cos_ang > one - epsilon(one)*two) then
  cos_ang = one
 else if(cos_ang < -one + epsilon(one)*two) then
  cos_ang = -one
 end if

 angle_ang=acos(cos_ang)

end function angle_ang
!!***

!!****if* ABINIT/angle_dihedral
!! NAME
!! angle_dihedral
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 function angle_dihedral(r1,r2,r3,r4)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => angle_dihedral
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: angle_dihedral
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3),r4(3)

!Local variables------------------------------------
!scalars
 real(dp) :: cos_dihedral,dih_sign,n1,n2,sin_dihedral
!arrays
 real(dp) :: cp1232(3),cp3432(3),cpcp(3),rpt12(3),rpt32(3),rpt34(3)

!******************************************************************

 rpt12(:) = r1(:)-r2(:)
 rpt32(:) = r3(:)-r2(:)
 rpt34(:) = r3(:)-r4(:)

 call cross_product(rpt12,rpt32,cp1232)
 call cross_product(rpt34,rpt32,cp3432)

!DEBUG
!write (*,*) ' cos_dihedral : cp1232 = ', cp1232
!write (*,*) ' cos_dihedral : cp3432 = ', cp3432
!ENDDEBUG

 n1 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 n2 = sqrt(cp3432(1)**2+cp3432(2)**2+cp3432(3)**2)

 cos_dihedral = (cp1232(1)*cp3432(1)+cp1232(2)*cp3432(2)+cp1232(3)*cp3432(3))/n1/n2
!we use complementary of standard angle, so
 cos_dihedral = -cos_dihedral

 call cross_product(cp1232,cp3432,cpcp)
 cpcp(:) = cpcp(:)/n1/n2
 sin_dihedral = -(cpcp(1)*rpt32(1)+cpcp(2)*rpt32(2)+cpcp(3)*rpt32(3))&
& /sqrt(rpt32(1)**2+rpt32(2)**2+rpt32(3)**2)
 dih_sign = one
!if (abs(sin_dihedral) > tol12) then
!dih_sign = sin_dihedral/abs(sin_dihedral)
!end if
 if (sin_dihedral < -tol12) then
  dih_sign = -one
 end if

!DEBUG
!write (*,'(a,3E20.10)') 'angle_dihedral : cos sin dih_sign= ',&
!&    cos_dihedral,sin_dihedral,dih_sign
!ENDDEBUG

 if (cos_dihedral > one - epsilon(one)*two) then
  cos_dihedral = one
 else if(cos_dihedral < -one + epsilon(one)*two) then
  cos_dihedral = -one
 end if

 angle_dihedral = dih_sign*acos(cos_dihedral)

end function angle_dihedral
!!***
