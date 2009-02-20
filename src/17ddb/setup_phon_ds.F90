!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_phon_ds
!!
!! NAME
!! setup_phon_ds
!!
!! FUNCTION
!! This routine copies scalars and arrays into structure phon_ds
!!  to be passed later to phonon interpolation routines.
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!!
!! INPUTS
!!   acell =  input length scales of cell (bohr)
!!   amu = mass of the atoms (atomic mass unit)
!!   atmfrc = inter-atomic force constants from anaddb
!!   dielt = dielectric tensor
!!   dipdip = dipole dipole interaction flag
!!   dyewq0 = atomic self-interaction correction to the
!!        dynamical matrix (only when dipdip=1)
!!   gmet = metric in reciprocal space (telphint=1)
!!   gprim =dimensionless basis vectors of reciprocal space
!!   indsym = mapping of atoms btw themselves under symmetry
!!   mpert = maximum number of ipert
!!   natom = number of atoms in cell
!!   nsym = number of space group symmetries
!!   ntypat = number of types of atoms
!!   nrpt =number of real space points used to integrate IFC (for
!!        interpolation of dynamical matrices)
!!   rcan = canonical positions of atoms
!!   rmet = metric tensor in real space (bohr^2)
!!   rprim =  primitive translation vectors (normalized)
!!   rprimd =  primitive translation vectors (dimensional)
!!   rpt = canonical positions of R points in the unit cell
!!   symrel = 3x3 matrices of the group symmetries (real space)
!!   trans = Atomic translations : xred = rcan + trans
!!   typat = type integer for each atom in cell
!!   ucvol = unit cell volume in bohr**3
!!   wghatm = Weight for the pair of atoms and the R vector
!!   xred = fractional dimensionless atomic coordinates
!!   zeff = effective charge on each atom, versus electric
!!        field and atomic displacement
!!
!! OUTPUT
!!   phon_ds = data structure for phonon interpolation - filled and allocated
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setup_phon_ds(phon_ds,dipdip,mpert,nsym,natom,ntypat,nrpt,&
&     ucvol,indsym,symrel,typat,acell,amu,atmfrc,dielt,dyewq0,gprim,gmet,&
&     xred,zeff,rcan,rmet,rprim,rprimd,rpt,trans,wghatm)

 use defs_basis
 use defs_datatypes
 use defs_elphon

 implicit none

 !Arguments ------------------------------------
 !scalars
 integer,intent(in) :: dipdip,mpert,nsym,natom,ntypat,nrpt
 real(dp),intent(in) :: ucvol
 type(phon_type) :: phon_ds
 !arrays
 integer,intent(in) :: indsym(4,nsym,natom)
 integer,intent(in) :: symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: acell(3),amu(ntypat),atmfrc(2,3,natom,3,natom,nrpt)
 real(dp),intent(in) :: dielt(3,3),dyewq0(3,3,natom),gmet(3,3),gprim(3,3)
 real(dp),intent(in) :: xred(3,natom),zeff(3,3,natom),rcan(3,natom)
 real(dp),intent(in) :: rmet(3,3),rprim(3,3),rprimd(3,3),rpt(3,nrpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt)
 !Local variables-------------------------------
! *************************************************************************
 allocate(phon_ds%indsym(4,nsym,natom))
 allocate(phon_ds%symrel(3,3,nsym))
 allocate(phon_ds%typat(natom))
 allocate(phon_ds%acell(3))
 allocate(phon_ds%amu(ntypat))
 allocate(phon_ds%atmfrc(2,3,natom,3,natom,nrpt))
 allocate(phon_ds%dielt(3,3))
 allocate(phon_ds%dyewq0(3,3,natom))
 allocate(phon_ds%gprim(3,3))
 allocate(phon_ds%gmet(3,3))
 allocate(phon_ds%xred(3,natom))
 allocate(phon_ds%zeff(3,3,natom))
 allocate(phon_ds%rcan(3,natom))
 allocate(phon_ds%rmet(3,3))
 allocate(phon_ds%rprim(3,3))
 allocate(phon_ds%rprimd(3,3))
 allocate(phon_ds%rpt(3,nrpt))
 allocate(phon_ds%trans(3,natom))
 allocate(phon_ds%wghatm(natom,natom,nrpt))
 phon_ds%dipdip = dipdip
 phon_ds%mpert = mpert
 phon_ds%nsym = nsym
 phon_ds%natom = natom
 phon_ds%ntypat = ntypat
 phon_ds%nrpt = nrpt
 phon_ds%ucvol = ucvol
 phon_ds%indsym(:,:,:) = indsym(:,:,:)
 phon_ds%symrel(:,:,:) = symrel(:,:,:)
 phon_ds%typat(:) = typat(:)
 phon_ds%acell(:) = acell(:)
 phon_ds%amu(:) = amu(:)
 phon_ds%atmfrc(:,:,:,:,:,:) = atmfrc(:,:,:,:,:,:)
 phon_ds%dielt(:,:) = dielt(:,:)
 phon_ds%dyewq0(:,:,:) = dyewq0(:,:,:)
 phon_ds%gprim(:,:) = gprim(:,:)
 phon_ds%gmet(:,:) = gmet(:,:)
 phon_ds%xred(:,:) = xred(:,:)
 phon_ds%zeff(:,:,:) = zeff(:,:,:)
 phon_ds%rcan(:,:) = rcan(:,:)
 phon_ds%rmet(:,:) = rmet(:,:)
 phon_ds%rprim(:,:) = rprim(:,:)
 phon_ds%rprimd(:,:) = rprimd(:,:)
 phon_ds%rpt(:,:) = rpt(:,:)
 phon_ds%trans(:,:) = trans(:,:)
 phon_ds%wghatm(:,:,:) = wghatm(:,:,:)
end subroutine setup_phon_ds
!!***

