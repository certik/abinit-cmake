!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_b_matrix
!! NAME
!! calc_b_matrix
!!
!! FUNCTION
!!  calculate values of derivatives of internal coordinates as a function of
!!  cartesian ones =  B matrix
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
!! ncart=number of auxiliary cartesian atom coordinates (used for constraints)
!! ndihed= number of dihedrals
!! ninternal=nbond+nang+ndihed+ncart: number of internal coordinates
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! b_matrix(ninternal,3*natom)=matrix of derivatives of internal coordinates
!!   wrt cartesians
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      xcart2deloc
!!
!! CHILDREN
!!      cross_product
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_b_matrix(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
& dtset,nrshift,rprimd,rshift,xcart,b_matrix)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => calc_b_matrix
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nang,nbond,ncart,ndihed,ninternal,nrshift
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,pointer :: angs(:,:,:),bonds(:,:,:),carts(:,:),dihedrals(:,:,:)
 real(dp),intent(in) :: rprimd(3,3),rshift(3,nrshift),xcart(3,dtset%natom)
 real(dp),intent(out) :: b_matrix(ninternal,3*dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,i4,iang,ibond,icart,idihed,iprim,s1,s2,s3,s4
!arrays
 real(dp) :: bb(3),r1(3),r2(3),r3(3),r4(3)

! *************************************************************************

!DEBUG
!write (*,*) 'calc_b_matrix : enter'
!ENDDEBUG

 iprim=0
 b_matrix(:,:) = zero

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
  iprim=iprim+1
  call dbond_length_d1(r1,r2,bb)
  b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
  call dbond_length_d1(r2,r1,bb)
  b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
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
  iprim=iprim+1
  call dang_d1(r1,r2,r3,bb)
  b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
  call dang_d2(r1,r2,r3,bb)
  b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
  call dang_d1(r3,r2,r1,bb)
  b_matrix(iprim,3*(i3-1)+1:3*i3) = b_matrix(iprim,3*(i3-1)+1:3*i3) + bb(:)
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
! write (*,*) 'dihed ',idihed
! write (*,*) r1
! write (*,*) r2
! write (*,*) r3
! write (*,*) r4

  iprim=iprim+1
  call ddihedral_d1(r1,r2,r3,r4,bb)
  b_matrix(iprim,3*(i1-1)+1:3*i1) = b_matrix(iprim,3*(i1-1)+1:3*i1) + bb(:)
  call ddihedral_d2(r1,r2,r3,r4,bb)
  b_matrix(iprim,3*(i2-1)+1:3*i2) = b_matrix(iprim,3*(i2-1)+1:3*i2) + bb(:)
  call ddihedral_d2(r4,r3,r2,r1,bb)
  b_matrix(iprim,3*(i3-1)+1:3*i3) = b_matrix(iprim,3*(i3-1)+1:3*i3) + bb(:)
  call ddihedral_d1(r4,r3,r2,r1,bb)
  b_matrix(iprim,3*(i4-1)+1:3*i4) = b_matrix(iprim,3*(i4-1)+1:3*i4) + bb(:)
 end do

 do icart=1,ncart
  iprim=iprim+1
  b_matrix(iprim,3*(carts(2,icart)-1)+carts(1,icart)) = &
&  b_matrix(iprim,3*(carts(2,icart)-1)+carts(1,icart)) + one
 end do

!DEBUG
!write (200,*) 'calc_b_matrix : b_matrix = '
!do iprim=1,ninternal
!write (200,'(6E16.6)') b_matrix(iprim,:)
!write (200,*)
!end do
!ENDDEBUG


end subroutine calc_b_matrix
!!***


!!****f* ABINIT/dbond_length_d1
!! NAME
!! dbond_length_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      cross_product
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group ( ).
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
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dbond_length_d1(r1,r2,bb)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => dbond_length_d1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!arrays
 real(dp) :: rpt(3)

!************************************************************************
 rpt(:) = r1(:)-r2(:)
 bb(:) = rpt(:)/bond_length(r1,r2)

end subroutine dbond_length_d1
!!***


!!****f* ABINIT/dang_d1
!! NAME
!! dang_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      cross_product
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
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dang_d1(r1,r2,r3,bb)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => dang_d1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!scalars
 real(dp) :: cos_ang,n1,n1232,n2,tmp
!arrays
 real(dp) :: cp1232(3),rpt(3),rpt12(3),rpt32(3)

!************************************************************************
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

 rpt(:) = rpt32(:)/n1/n2 - rpt12(:)*cos_ang/n1/n1

 tmp = sqrt(one-cos_ang**2)
 bb(:) = zero
 if (tmp > epsilon(one)) then
  bb(:) = rpt(:) * (-one)/tmp
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson
 call cross_product(rpt12,rpt32,cp1232)
 n1232 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 rpt(:) = (cos_ang*rpt12(:)*n2/n1 - rpt32(:))/n1232
 if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3)) > tol10) then
  write (*,*) 'Compare bb ang 1 : '
  write (*,*) bb(:), rpt(:), bb(:)-rpt(:)
 end if
 bb(:) = rpt(:)

end subroutine dang_d1
!!***


!!****f* ABINIT/dang_d2
!! NAME
!! dang_d2
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      cross_product
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dang_d2(r1,r2,r3,bb)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => dang_d2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------
!scalars
 real(dp) :: cos_ang,n1,n1232,n2,tmp
!arrays
 real(dp) :: cp1232(3),rpt(3),rpt12(3),rpt32(3)

!************************************************************************
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

 rpt(:) = -rpt32(:)/n1/n2 - rpt12(:)/n1/n2 &
& + rpt12(:)*cos_ang/n1/n1 + rpt32(:)*cos_ang/n2/n2

 tmp = sqrt(one-cos_ang**2)
 bb(:) = zero
 if (tmp > tol12) then
  bb(:) = rpt(:) * (-one)/tmp
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson
 call cross_product(rpt12,rpt32,cp1232)
 n1232 = sqrt(cp1232(1)**2+cp1232(2)**2+cp1232(3)**2)
 rpt(:) = ((n1-n2*cos_ang)*rpt12(:)/n1 + (n2-n1*cos_ang)*rpt32(:)/n2) / n1232
 if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
  write (*,*) 'Compare bb ang 2 : '
  write (*,*) bb(:), rpt(:), bb(:)-rpt(:)
 end if
 bb(:) = rpt(:)


end subroutine dang_d2
!!***


!!****f* ABINIT/ddihedral_d1
!! NAME
!! ddihedral_d1
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      cross_product
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ddihedral_d1(r1,r2,r3,r4,bb)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => ddihedral_d1
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3),r4(3)
 real(dp),intent(out) :: bb(3)

!Local variables ------------------------------------
!scalars
 real(dp) :: cos_dihedral,dih_sign,n1,n2,n23,sin_dihedral,tmp
!arrays
 real(dp) :: cp1232(3),cp32_1232(3),cp32_3432(3),cp3432(3),cpcp(3),rpt(3)
 real(dp) :: rpt12(3),rpt32(3),rpt34(3)

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
 if (cos_dihedral > one - epsilon(one)*two) then
  cos_dihedral = one
 else if(cos_dihedral < -one + epsilon(one)*two) then
  cos_dihedral = -one
 end if
!we use complementary of standard angle, so
!cos_dihedral = -cos_dihedral

 call cross_product(cp1232,cp3432,cpcp)
 cpcp(:) = cpcp(:)/n1/n2
!we use complementary of standard angle, but sin is invariant
 sin_dihedral = -(cpcp(1)*rpt32(1)+cpcp(2)*rpt32(2)+cpcp(3)*rpt32(3))&
& /sqrt(rpt32(1)**2+rpt32(2)**2+rpt32(3)**2)
 dih_sign = one
 if (sin_dihedral < -epsilon(one)) then
  dih_sign = -one
 end if

!DEBUG
!write (*,'(a,3E16.6)') 'ddihedral_d1 : cos abs(sin) dih_sign= ',&
!&    cos_dihedral,sin_dihedral,dih_sign
!ENDDEBUG

!ddihedral_d1 = dih_sign* acos(cos_dihedral)
 call cross_product(rpt32,cp1232,cp32_1232)
 call cross_product(rpt32,cp3432,cp32_3432)

 rpt(:) = cp32_3432(:)/n1/n2 - cp32_1232(:)/n1/n1 * cos_dihedral
 bb(:) = zero

!DEBUG
!write (*,*) 'ddihedral_d1 cp1232 cp3432 = ',cp1232,cp3432,rpt32
!write (*,*) 'ddihedral_d1 cp32_1232 cp32_3432 = ',cp32_1232,cp32_3432,cos_dihedral,n1,n2
!write (*,*) 'ddihedral_d1 rpt = ',rpt
!ENDDEBUG

 tmp = sqrt(one-cos_dihedral**2)
 if (tmp > tol12) then
! we use complementary of standard angle, so cosine in acos has - sign,
! and it appears for the derivative
  bb(:) = -dih_sign * rpt(:) * (-one) / tmp
 else
  bb(:) = dih_sign * cp32_3432(:) / n1 / n2 / &
&  sqrt(cp32_3432(1)**2+cp32_3432(2)**2+cp32_3432(3)**2)
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson

 n23 = sqrt(rpt32(1)*rpt32(1)+rpt32(2)*rpt32(2)+rpt32(3)*rpt32(3))
 rpt(:) = cp1232(:)*n23/n1/n1
!if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
!write (*,*) 'Compare bb1 : '
!write (*,*) bb(:), rpt(:), bb(:)-rpt(:)
!end if
 bb(:) = rpt(:)

end subroutine ddihedral_d1
!!***


!!****f* ABINIT/ddihedral_d2
!! NAME
!! ddihedral_d2
!!
!! FUNCTION
!!
!! PARENTS
!!      calc_b_matrix
!!
!! CHILDREN
!!      cross_product
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ddihedral_d2(r1,r2,r3,r4,bb)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => ddihedral_d2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: r1(3),r2(3),r3(3),r4(3)
 real(dp),intent(out) :: bb(3)

!Local variables
!scalars
 real(dp) :: cos_dihedral,dih_sign,n1,n2,n23,sin_dihedral,sp1232,sp3432,tmp
!arrays
 real(dp) :: cp1232(3),cp1232_12(3),cp1232_34(3),cp32_1232(3),cp32_3432(3)
 real(dp) :: cp3432(3),cp3432_12(3),cp3432_34(3),cpcp(3),rpt(3),rpt12(3)
 real(dp) :: rpt32(3),rpt34(3)

! *************************************************************************
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
 if (cos_dihedral > one - epsilon(one)*two) then
  cos_dihedral = one
 else if(cos_dihedral < -one + epsilon(one)*two) then
  cos_dihedral = -one
 end if
!we use complementary of standard angle, so
!cos_dihedral = -cos_dihedral

 call cross_product(cp1232,cp3432,cpcp)
 cpcp(:) = cpcp(:)/n1/n2
!we use complementary of standard angle, but sin is invariant
 sin_dihedral = -(cpcp(1)*rpt32(1)+cpcp(2)*rpt32(2)+cpcp(3)*rpt32(3))&
& /sqrt(rpt32(1)**2+rpt32(2)**2+rpt32(3)**2)
 dih_sign = one
 if (sin_dihedral <  -tol12) then
  dih_sign = -one
 end if

!DEBUG
!write (*,'(a,3E16.6)') 'ddihedral_d2 : cos abs(sin) dih_sign= ',&
!&    cos_dihedral,sin_dihedral,dih_sign
!ENDDEBUG

!ddihedral_d2 = dih_sign* acos(cos_dihedral)
 call cross_product(rpt32,cp3432,cp32_3432)
 call cross_product(cp3432,rpt12,cp3432_12)
 call cross_product(cp1232,rpt34,cp1232_34)

 call cross_product(rpt32,cp1232,cp32_1232)
 call cross_product(cp1232,rpt12,cp1232_12)
 call cross_product(cp3432,rpt34,cp3432_34)

 rpt(:) = -(cp32_3432(:) + cp3432_12(:) + cp1232_34(:))/n1/n2 &
& +cos_dihedral*(cp32_1232(:)/n1/n1 + cp1232_12(:)/n1/n1 + cp3432_34(:)/n2/n2)
 bb(:) = zero
 tmp = sqrt(one-cos_dihedral**2)
 if (tmp > tol12) then
! we use complementary of standard angle, so cosine in acos has - sign,
! and it appears for derivative
  bb(:) = -dih_sign * rpt(:) * (-one) / tmp
 else
  bb(:) = dih_sign * cos_dihedral * &
&  ( cp32_1232(:)/n1/n1/sqrt(cp32_1232(1)**2+cp32_1232(2)**2+cp32_1232(3)**2) &
&  +cp1232_12(:)/n1/n1/sqrt(cp1232_12(1)**2+cp1232_12(2)**2+cp1232_12(3)**2) &
&  +cp3432_34(:)/n2/n2/sqrt(cp3432_34(1)**2+cp3432_34(2)**2+cp3432_34(3)**2) )
 end if

!TEST: version from MOLECULAR VIBRATIONS EB Wilson p. 61
 n23 = sqrt(rpt32(1)*rpt32(1)+rpt32(2)*rpt32(2)+rpt32(3)*rpt32(3))
 sp1232 = rpt12(1)*rpt32(1)+rpt12(2)*rpt32(2)+rpt12(3)*rpt32(3)
 sp3432 = rpt34(1)*rpt32(1)+rpt34(2)*rpt32(2)+rpt34(3)*rpt32(3)

 rpt(:) = -cp1232(:)*(n23-sp1232/n23)/n1/n1 - cp3432(:)*sp3432/n23/n2/n2
!if (abs(bb(1)-rpt(1))+abs(bb(2)-rpt(2))+abs(bb(3)-rpt(3))  > tol10) then
!write (*,*) 'Compare bb2 : '
!write (*,*) bb(:), rpt(:), bb(:)-rpt(:)
!write (*,*) -cp1232(:)*(n23-sp1232/n23)/n1/n1, -cp3432(:)*sp3432/n23/n2/n2
!end if
 bb(:) = rpt(:)


end subroutine ddihedral_d2
!!***
