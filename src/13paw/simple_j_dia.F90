!{\src2tex{textfont=tt}}
!!****f* ABINIT/gipaw_j_dia
!! NAME
!! simple_j_dia
!!
!! FUNCTION
!! simple test current for H atoms
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT

!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      acrossb
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine simple_j_dia(jdia,natom,nfft,pawfgrtab)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft
!arrays
 real(dp),intent(out) :: jdia(3,3,nfft)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,ifgd,ifftsph
 real(dp) :: Bvecsum,nrm,nrmmax,psi2,scale
!arrays
 real(dp) :: axb(3),B(3,3),rvec(3)
 real(dp),allocatable :: Bvec(:,:)

! ************************************************************************

!DEBUG
!write(*,*)' gipaw_j_dia : enter'
!ENDDEBUG

!make sure current field is set to zero everywhere
 jdia(:,:,:) = zero

!define the directions of the external B field
 B = zero
 B(1,1) = 1.0; B(2,2) = 1.0; B(3,3) = 1.0;

!loop over atoms in cell
 scale = 0.5
 nrmmax = 0.0
 Bvecsum = 0.0
 do iatom = 1, natom
  allocate(Bvec(3,pawfgrtab(iatom)%nfgd))
  do ifgd=1, pawfgrtab(iatom)%nfgd
   ifftsph = pawfgrtab(iatom)%ifftsph(ifgd)
   rvec(:) = pawfgrtab(iatom)%rfgd(:,ifgd)
   nrm = sqrt(dot_product(rvec,rvec))
   if (nrm > nrmmax) nrmmax = nrm
   psi2=exp(-2.0*nrm/scale)/(pi*scale*scale*scale)
   do idir=1, 3
    call acrossb(B(idir,:),rvec,axb)
    jdia(idir,:,ifftsph) = -psi2*axb(:)/(2.0*Sp_Lt)
   end do ! end loop over idir for jdia
   if (nrm > zero) then
    call acrossb(rvec,jdia(3,:,ifftsph),axb)
    Bvec(:,ifgd) = axb(:)/(Sp_Lt*nrm*nrm*nrm)
   end if
   Bvecsum = Bvecsum + Bvec(3,ifgd)
  end do ! end loop over nfgd points in sphere
  write(6,'(i8,f12.8,f12.8)')pawfgrtab(iatom)%nfgd,Bvecsum,Bvecsum*(4.0*pi*nrm**3/3.0)/pawfgrtab(iatom)%nfgd
  deallocate(Bvec)
 end do ! end loop over atoms in cell

!DEBUG
!write(6,*)' simple_j_dia : exit '
!stop
!ENDDEBUG

 end subroutine simple_j_dia
!!***
