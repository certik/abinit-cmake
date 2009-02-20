!{\src2tex{textfont=tt}}
!!****f* ABINIT/setsym
!! NAME
!! setsym
!!
!! FUNCTION
!! Set up irreducible zone in  G space by direct calculation.
!! Do not call this routine if nsym=1 (only identity symmetry).
!! Only indsym and symrec get returned if iscf=0.
!! symrec needed to symmetrize coordinate gradients in sygrad.
!! (symrec is redundant and could be removed later in favor of symrel)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! densymop <type(dens_sym_operator_type)>=the density symmetrization operator
!! iscf=(<= 0  =>non-SCF), >0 => SCF
!! natom=number of atoms in unit cell
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nspden=number of spin-density components
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetries in space group (at least 1)
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in terms of real space
!! primitive translations
!! tnons(3,nsym)=nonsymmorphic translations of space group in terms
!! of real space primitive translations (may be 0)
!! typat(natom)=atom type (integer) for each atom
!! xred(3,natom)=atomic coordinates in terms of real space primitive
!! translations
!!
!! OUTPUT
!! indsym(4,nsym,natom)=indirect indexing of atom labels--see subroutine
!!   symatm for definition (if nsym>1)
!! irrzon(nfft,2,nspden/nsppol)=irreducible zone data
!! phnons(2,nfft,nspden/nsppol)=nonsymmorphic translation phases
!! symrec(3,3,nsym)=symmetry operations in terms of reciprocal
!!   space primitive translations (if nsym>1)
!!
!! NOTES
!! nsppol and nspden are needed in case of (anti)ferromagnetic symmetry operations
!!
!! PARENTS
!!      gstate,loper3,nonlinear,respfn,scfcv,suscep
!!
!! CHILDREN
!!      chkgrp,irrzg,mati3inv,symatm,symdet,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setsym(densymop,indsym,irrzon,iscf,natom,&
& nfft,ngfft,nspden,nsppol,nsym,phnons,&
& symafm,symrec,symrel,tnons,typat,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_13recipspace, except_this_one => setsym
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iscf,natom,nfft,nspden,nsppol,nsym
 type(dens_sym_operator_type),intent(in) :: densymop
!arrays
 integer,intent(in) :: ngfft(18),symafm(nsym),symrel(3,3,nsym),typat(natom)
 integer,intent(out) :: indsym(4,nsym,natom),irrzon(nfft,2,nspden/nsppol)
 integer,intent(out) :: symrec(3,3,nsym)
 real(dp),intent(in) :: tnons(3,nsym),xred(3,natom)
 real(dp),intent(out) :: phnons(2,nfft,nspden/nsppol)

!Local variables-------------------------------
!scalars
 integer,save :: count=0
 integer :: isym
!arrays
 integer,allocatable :: determinant(:)
 real(dp) :: tsec(2)

! *************************************************************************

!DEBUG
!write(6,*)' setsym : enter '
!count=count+1
!write(6,*)' count=',count
!if(count==3)stop
!ENDDEBUG

 call timab(6,1,tsec)

!Check that symmetries have unity determinant
 allocate(determinant(nsym))
 call symdet(determinant,nsym,symrel)
 deallocate(determinant)

!DEBUG
!write(6,*)' setsym : after symdet '
!stop
!ENDDEBUG

!Get the symmetry matrices in terms of reciprocal basis
 do isym=1,nsym
  call mati3inv(symrel(:,:,isym),symrec(:,:,isym))
 end do

!DEBUG
!write(6,*)' setsym : after invt '
!stop
!ENDDEBUG

!Check for group closure
 call chkgrp(nsym,symafm,symrel)
 call chkgrp(nsym,symafm,symrec)

!DEBUG
!write(6,*)' setsym : after chkgrp '
!stop
!ENDDEBUG

!Obtain a list of rotated atom labels:
 call symatm(indsym,natom,nsym,symrec,tnons,typat,xred)

!DEBUG
!write(6,*)' setsym : before irrzg '
!write(6,*)' nspden,nsppol,nsym,ngfft(1),ngfft(2),ngfft(3),nfft',&
!&            nspden,nsppol,nsym,ngfft(1),ngfft(2),ngfft(3),nfft
!do isym=1,nsym
!write(6,*)' isym,symafm(isym)=',isym,symafm(isym)
!write(6,*)'   tnons(1:3,isym)',tnons(1:3,isym)
!write(6,*)'   symrel(1:3,1,isym)',symrel(1:3,1,isym)
!write(6,*)'   symrel(1:3,2,isym)',symrel(1:3,2,isym)
!write(6,*)'   symrel(1:3,3,isym)',symrel(1:3,3,isym)
!end do
!irrzon(:,:,:)=0
!phnons(:,:,:)=zero
!if(count==3)stop
!ENDDEBUG

!If non-SCF calculation, or nsym==1, do not need IBZ data
 if ( (iscf>0 .or. iscf==-3) .and. nsym>1 ) then
! Locate irreducible zone in reciprocal space for symmetrization:
  call irrzg(densymop,irrzon,nspden,nsppol,nsym,ngfft(1),ngfft(2),ngfft(3),phnons,&
&  symafm,symrel,tnons)
 end if

 call timab(6,2,tsec)

!DEBUG
!if(count==3)stop
!write(6,*)' setsym : exit '
!ENDDEBUG

end subroutine setsym
!!***
