!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcart2deloc_fixb
!! NAME
!! xcart2deloc_fixb
!!
!! FUNCTION
!!  calculate values of delocalized coordinates as a function of
!!  cartesian ones. First primitive internals, then B matrix, then F, then U
!!  then delocalized internals.
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
!! ncart=number of cartesian coordinates included
!! ndihed= number of dihedrals
!! ninternal=nbond+nang+ndihed+ncart: number of internal coordinates
!! nrshift= dimension of rshift
!! rprimd(3,3)=dimensional real space primitive translations (bohr)
!! rshift(3,nrshift)=shift in xred that must be done to find all neighbors of
!!                   a given atom within a given number of neighboring shells
!! u_matrix(ninternal,3*(natom-1))=eigenvectors of BB^{T}
!! xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!! deloc_int(3*(natom-1))=delocalized internal coordinates
!! prim_int(ninternal)=primitive internal coordinates
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      calc_prim_int,dgemv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xcart2deloc_fixb(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
& dtset,nrshift,rprimd,rshift,xcart,&
& u_matrix,deloc_int,prim_int)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16geomoptim, except_this_one => xcart2deloc_fixb
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nang,nbond,ncart,ndihed,ninternal,nrshift
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,pointer :: angs(:,:,:),bonds(:,:,:),carts(:,:),dihedrals(:,:,:)
 real(dp),intent(in) :: rprimd(3,3),rshift(3,nrshift)
 real(dp),intent(in) :: u_matrix(ninternal,3*(dtset%natom-1))
 real(dp),intent(in) :: xcart(3,dtset%natom)
 real(dp),intent(out) :: deloc_int(3*(dtset%natom-1)),prim_int(ninternal)

!Local variables-------------------------------
!scalars
 integer :: i1,ibond
!no_abirules

! *************************************************************************

!DEBUG
!write (*,*) 'xcart2deloc_fixb : enter'
!do ibond=1,nbond
!do i1=1,2
!write (*,'(2I5)') bonds(:,i1,ibond)
!end do
!end do
!ENDDEBUG


 call calc_prim_int(angs,bonds,carts,dihedrals,nbond,nang,ndihed,ncart,ninternal,&
& dtset,nrshift,rprimd,rshift,xcart,prim_int)

!calculate value of delocalized internals with the given U matrix


 call dgemv('T',ninternal,3*(dtset%natom-1),one,&
& u_matrix,ninternal,prim_int,1,zero,deloc_int,1)

!write (*,'(a,6E16.6)') 'xcart2deloc_fixb : deloc_int = ',  deloc_int


end subroutine xcart2deloc_fixb
!!***
