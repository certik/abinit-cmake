!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_dos_1band
!! NAME
!! get_dos_1band
!!
!! FUNCTION
!! calculate DOS from tetrahedron method for 1 band and 1 sppol
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dos_fractions=fractional DOS at each irred kpoint
!! enemin=minimal energy for DOS
!! enemax=maximal energy for DOS
!! nene=number of energies for DOS
!! nkpt=number of irreducible kpoints
!! ndosfraction=number of different fractional DOSs
!! tweight=sum of tetrahedron weights for each irred kpoint
!! dtweightde=energy derivative of tweight
!!
!! OUTPUT
!!  partial_dos(nene,ndosfraction)=partial DOS, for each different channel
!!  integ_dos(nene,ndosfraction)=integrated DOS, for each different channel
!!
!! PARENTS
!!      tetrahedron
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_dos_1band (dos_fractions,enemin,enemax,&
&            integ_dos,nene,nkpt,ndosfraction,&
&            partial_dos,tweight,dtweightde)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndosfraction,nene,nkpt
 real(dp),intent(in) :: enemax,enemin
!arrays
 real(dp),intent(in) :: dos_fractions(nkpt,ndosfraction),dtweightde(nkpt,nene)
 real(dp),intent(in) :: tweight(nkpt,nene)
 real(dp),intent(out) :: integ_dos(nene,ndosfraction)
 real(dp),intent(out) :: partial_dos(nene,ndosfraction)

!Local variables-------------------------------
!scalars
 integer :: ieps,ifract,ikpt,itetra
 real(dp) :: deltaene,eigmid,eps,max_occ,prefac_tetr,tmp,tmp1,tmp2
!arrays
 real(dp) :: eigen_in(nkpt)

! *********************************************************************

!DEBUG
!write (*,*) ' get_dos_1band : enter'
!ENDDEBUG

 partial_dos(:,:) = zero
 integ_dos(:,:) = zero

 deltaene = (enemax-enemin) / (nene-1)
!
!Calculate parameters of DOS at each point eps in [epsmin,epsmax]
!
 eps=enemin
 do ieps=1,nene

  do ifract=1,ndosfraction
   do ikpt=1,nkpt
    partial_dos(ieps,ifract) = partial_dos(ieps,ifract)+&
&    dtweightde(ikpt,ieps)*dos_fractions(ikpt,ifract)
    integ_dos(ieps,ifract) = integ_dos(ieps,ifract)+&
&    tweight(ikpt,ieps)*dos_fractions(ikpt,ifract)
   end do
  end do
! DEBUG
! tmp = zero
! do ikpt=1,nkpt
! tmp = tmp + dtweightde(ikpt,ieps)
! end do
! tmp2 = zero
! do ikpt=1,nkpt
! tmp2 = tmp2 + tweight(ikpt,ieps)
! end do
! write (*,*) ' pDOStDOSiDOS', ieps, eps, partial_dos(ieps,1), tmp, tmp2
! ENDDEBUG
  eps = eps + deltaene
 end do
!end ieps do



!DEBUG
!write(6,*)' get_dos_1band : exit '
!ENDDEBUG

end subroutine get_dos_1band
!!***
