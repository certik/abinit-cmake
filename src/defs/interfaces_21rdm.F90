!!****m* ABINIT/interfaces_21rdm
!! NAME
!! interfaces_21rdm
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/21rdm
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_21rdm

 implicit none

interface
 subroutine crho(paral_kgb,ngfft,irottb,nbnds,timrev,nkbz,nkibz,nsym,symrel,tnons,symafm,&  
  &  nfftot,nspden,nsppol,occ,omegaplasma,rho,ucvol,wfr,wtk,mpi_enreg,my_minb,my_maxb)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: nbnds
  integer,intent(in) :: nfftot
  integer,intent(in) :: nkbz
  integer,intent(in) :: nkibz
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: timrev
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp),intent(out) :: omegaplasma
  real(dp),intent(in) :: ucvol
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: irottb(nfftot,nsym)
  real(dp),intent(in) :: occ(nkibz,nbnds,nsppol)
  real(dp),intent(out) :: rho(nfftot,nsppol)
  integer,intent(in) :: symafm(nsym)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  complex,intent(in) :: wfr(nfftot,my_minb:my_maxb,nkibz,nsppol)
  real(dp),intent(in) :: wtk(nkibz)
 end subroutine crho
end interface

interface
 subroutine cvxclda(dtset,ixc,mpi_enreg,ngfft,nfftot,nsppol,rho,rprimd,vxclda)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: ixc
  integer,intent(in) :: nfftot
  integer,intent(in) :: nsppol
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: rho(nfftot,nsppol)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vxclda(nfftot,nsppol)
 end subroutine cvxclda
end interface

interface
 subroutine fftwfn(paral_kgb,npwwfn,my_minb,my_maxb,nkibz,nsppol,wfg,wfr,igfft,ngfft,tim_fourdp,mpi_enreg)
  use defs_datatypes
  implicit none
  integer,intent(in) :: my_maxb
  integer,intent(in) :: my_minb
  integer,intent(in) :: nkibz
  integer,intent(in) :: npwwfn
  integer,intent(in) :: nsppol
  integer,intent(in) :: paral_kgb
  integer,intent(in) :: tim_fourdp
  type(mpi_type),intent(inout) :: mpi_enreg
  integer,intent(in) :: ngfft(18)
  integer,intent(in) :: igfft(npwwfn)
  complex,intent(in) :: wfg(npwwfn,my_minb:my_maxb,nkibz,nsppol)
  complex,intent(out) :: wfr(ngfft(1)*ngfft(2)*ngfft(3),my_minb:my_maxb,nkibz,nsppol)
 end subroutine fftwfn
end interface

interface
 subroutine lattice(a1,a2,a3,b1,b2,b3,ucvol,bzvol)
  use defs_basis
  implicit none
  real(dp),intent(out) :: bzvol
  real(dp),intent(out) :: ucvol
  real(dp),intent(in) :: a1(3)
  real(dp),intent(in) :: a2(3)
  real(dp),intent(in) :: a3(3)
  real(dp),intent(out) :: b1(3)
  real(dp),intent(out) :: b2(3)
  real(dp),intent(out) :: b3(3)
 end subroutine lattice
end interface

interface
 subroutine occred(rdocc,nsppol,nk,nb)
  use defs_basis
  implicit none
  integer,intent(in) :: nb
  integer,intent(in) :: nk
  integer,intent(in) :: nsppol
  real(dp),intent(inout) :: rdocc(nk,nb,nsppol)
 end subroutine occred
end interface

interface
 subroutine old_setmesh(gmet,gvec,ngfft,npwvec,npwsigx,npwwfn,nfftot,method,mG0,nsym,symrel,tnons,enforce_sym)
  use defs_basis
  implicit none
  integer,intent(in) :: enforce_sym
  integer,intent(in) :: method
  integer,intent(out) :: nfftot
  integer,intent(in) :: npwsigx
  integer,intent(in) :: npwvec
  integer,intent(in) :: npwwfn
  integer,intent(in) :: nsym
  integer,intent(in) :: mG0(3)
  integer,intent(inout) :: ngfft(18)
  real(dp),intent(in) :: gmet(3,3)
  integer,intent(in) :: gvec(3,npwvec)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
 end subroutine old_setmesh
end interface

interface
 subroutine rdm(acell,dtfil,dtset,pawtab,mpi_enreg,rprim)
  use defs_basis
  use defs_datatypes
  implicit none
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: acell(3)
  type(pawtab_type),intent(inout) :: pawtab(dtset%ntypat*dtset%usepaw)
  real(dp),intent(in) :: rprim(3,3)
 end subroutine rdm
end interface

interface
 subroutine setup_G_rotation_old(only_one_kpt,nsym,symrec,timrev,npw,gvec,grottb,grottbm1)
  implicit none
  integer,intent(in) :: npw
  integer,intent(in) :: nsym
  integer,intent(in) :: timrev
  logical,intent(in) :: only_one_kpt
  integer,intent(inout) :: grottb(npw,timrev,nsym)
  integer,intent(inout) :: grottbm1(npw,timrev,nsym)
  integer,intent(in) :: gvec(3,npw)
  integer,intent(in) :: symrec(3,3,nsym)
 end subroutine setup_G_rotation_old
end interface

end module interfaces_21rdm
!!***
