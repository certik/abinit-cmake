!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_update
!! NAME
!! hdr_update
!!
!! FUNCTION
!! This subroutine update the header structured datatype.
!! Most of its records had been initialized correctly, but some corresponds
!! to evolving variables, or change with the context (like fform),
!! This routine is to be called before writing the header
!! to a file, in order to have up-to-date information.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! bantot=total number of bands
!! etot=total energy (Hartree)
!! fermie=Fermi energy (Hartree)
!! natom=number of atoms
!! residm=maximal residual
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! occ(bantot)=occupancies for each band and k point
!! pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! xred(3,natom)= relative coords of atoms in unit cell (dimensionless)
!!
!! OUTPUT
!! hdr <type(hdr_type)>=the header, initialized, and for most part of
!!   it, contain its definite values, except for evolving variables
!!
!! PARENTS
!!      afterscfloop,gstate,loper3,newsp,nonlinear,respfn,vtorho
!!
!! CHILDREN
!!      rhoij_copy
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine hdr_update(bantot,etot,fermie,hdr,natom,&
&                     residm,rprimd,occ,pawrhoij,usepaw,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot,natom,usepaw
 real(dp),intent(in) :: etot,fermie,residm
 type(hdr_type),intent(out) :: hdr
!arrays
 real(dp),intent(in) :: occ(bantot),rprimd(3,3),xred(3,natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(:)

!Local variables-------------------------------
!scalars
 integer :: irhoij,nrhoij

! *************************************************************************

!Update of the "evolving" data
 hdr%etot     =etot
 hdr%fermie   =fermie
 hdr%residm   =residm
 hdr%rprimd(:,:)=rprimd(:,:)
 hdr%occ(:)   =occ(:)
 hdr%xred(:,:)=xred(:,:)
 if (usepaw==1) then
  nrhoij=size(pawrhoij)
  if (nrhoij==natom) then
   call rhoij_copy(pawrhoij,hdr%pawrhoij)
  else if (nrhoij<natom) then
   call rhoij_copy(pawrhoij(1:nrhoij),hdr%pawrhoij(1:nrhoij))
   do irhoij=nrhoij+1,natom
    hdr%pawrhoij(irhoij)%nrhoijsel=0
    hdr%pawrhoij(irhoij)%rhoijselect(:)=0
    hdr%pawrhoij(irhoij)%rhoijp(:,:)=zero
   end do
  else
   call rhoij_copy(pawrhoij(1:natom),hdr%pawrhoij(1:natom))
  end if
 end if

end subroutine hdr_update
!!***
