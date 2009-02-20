!{\src2tex{textfont=tt}}
!!****f* ABINIT/gipaw_j_dia_aug
!! NAME
!! gipaw_j_dia_aug
!!
!! FUNCTION
!! Compute the current due to the diamagnetic augemtation term in GIPAW
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! integer natom,nfft,ntypat: number of atoms in cell, number of grid points, number of atom types
!! integer typat(natom): array of atom types
!! type(gipaw_type) gipaw_aug(natom): data structure giving dia and para current strengths around each atom
!! type(pawfgrtab_type) pawfgrtab(natom): data on fine grid in each paw sphere
!! type(pawrhoij_type) pawrhoij(natom): paw density in each sphere
!! type(pawtab_type) pawtab(ntypat): paw wavefunctions around each atom type
!!
!! OUTPUT
!! real(dp) jdia(B_idir=1..3,jdia_i=1..3,nfft): vector current field at each point on the grid for each
!!                                              direction of the external field B
!! NOTES
!! In GIPAW the diagmagnetic augmentation current at position r' is determined from
!! $-\frac{\mathbf{B\wedge(r'-R)}}{2c}\sum_{ij}\rho_{ij}\phi_j(r')\phi_i(r')-\tilde{\phi}_j(r')\tilde{\phi}_i(r')$,
!! summed over atomic sites $R$ (see Yates, Pickard, Mauri PRB 76, 024401 (2007) Eq. 19). Notice that in that
!! paper the approximation is made that the response at a given atomic site is obtained just from the sphere
!! around that atom; here we do not make that approximation and instead compute the response for this term
!! over the entire grid and then compute the induced fields just as all the other terms are handled.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine gipaw_j_dia_aug(gipaw_aug,jdia,natom,nfft,ntypat,pawfgrtab,pawrhoij,pawtab,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,ntypat
!arrays
 integer,intent(in) :: typat(natom)
 type(gipaw_type),intent(in) :: gipaw_aug(natom)
 real(dp),intent(out) :: jdia(3,3,nfft)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,ifgd,ifftsph,irhoij,ispden,itypat,klmn
 real(dp) :: rhoij_jdij
!arrays
 real(dp) :: B(3,3),Bxr(3),rvec(3)

! ************************************************************************

!DEBUG
!write(*,*)' gipaw_j_dia_aug : enter'
!ENDDEBUG

!make sure current field is set to zero everywhere
 jdia(:,:,:) = zero

!define the directions of the external B field
 B = zero
 B(1,1) = 1.0; B(2,2) = 1.0; B(3,3) = 1.0;

!loop over atoms in cell
 do iatom = 1, natom
  itypat = typat(iatom)

! loop over spin components
  do ispden=1,pawrhoij(iatom)%nspden

   do irhoij=1,pawrhoij(iatom)%nrhoijsel ! Loop over non-zero elements of rhoij
    klmn=pawrhoij(iatom)%rhoijselect(irhoij)
    do ifgd=1, pawfgrtab(iatom)%nfgd ! loop over fine grid points in current PAW sphere
     ifftsph = pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid
!    here is rhoij*[phi_i(r')phi_j(r')-tphi_i(r')tphi_j(r')]/2c; the minus sign makes it charge current, not density 
     rhoij_jdij = -pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*gipaw_aug(iatom)%dia(ifgd,klmn)/(2.0*Sp_Lt)
     rvec(:) = pawfgrtab(iatom)%rfgd(:,ifgd)
!    now accumulate -B x (r-R)*rhoij*jdij into the vector fields jdia
     do idir = 1, 3
      call acrossb(B(idir,:),rvec,Bxr)
      jdia(idir,:,ifftsph) = jdia(idir,:,ifftsph) - rhoij_jdij*Bxr(:)
     end do ! end loop over B field directions
    end do ! end loop on nfgd points around atom
   end do ! end loop on non-zero elements of rhoij

  end do ! end loop over nspden components
 end do ! end loop over atoms in cell

!DEBUG
!write(6,*)' gipaw_j_dia_aug : exit '
!stop
!ENDDEBUG

 end subroutine gipaw_j_dia_aug
!!***
