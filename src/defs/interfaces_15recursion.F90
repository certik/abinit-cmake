!!****m* ABINIT/interfaces_15recursion
!! NAME
!! interfaces_15recursion
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/15recursion
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

module interfaces_15recursion

 implicit none

interface
 subroutine entropyrec(an,bn2,nrec,trotter,ent_out, mult,&  
  &  prtvol,n_pt_integ,xmax,&  
  &  ent_out1,ent_out2,ent_out3,ent_out4)
  use defs_basis
  implicit none
  integer,intent(in) :: n_pt_integ
  integer,intent(in) :: nrec
  integer,intent(in) :: prtvol
  integer,intent(in) :: trotter
  real(dp),intent(out) :: ent_out
  real(dp),intent(out) :: ent_out1
  real(dp),intent(out) :: ent_out2
  real(dp),intent(out) :: ent_out3
  real(dp),intent(out) :: ent_out4
  real(dp),intent(in) :: mult
  real(dp),intent(in) :: xmax
  real(dp),intent(in) :: an(0:nrec)
  real(dp),intent(in) :: bn2(0:nrec)
 end subroutine entropyrec
end interface

interface
 subroutine fermisolverec(fermie,rho,a,b2,nb_rec,&  
  &  temperature,trotter,nelect,&  
  &  acc, max_it,&  
  &  longueur_tranche,mpi_enreg,rang,&  
  &  rmet,inf_ucvol,tim_fourdp)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: longueur_tranche
  integer,intent(in) :: max_it
  integer,intent(in) :: nb_rec
  integer,intent(in) :: rang
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: trotter
  real(dp),intent(in) :: acc
  real(dp),intent(inout) :: fermie
  real(dp),intent(in) :: inf_ucvol
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: nelect
  real(dp),intent(in) :: temperature
  real(dp),intent(inout) :: a(0:nb_rec,1:longueur_tranche)
  real(dp),intent(inout) :: b2(0:nb_rec,1:longueur_tranche)
  real(dp),intent(inout) :: rho(1:longueur_tranche)
  real(dp),intent(in) :: rmet(3,3)
 end subroutine fermisolverec
end interface

interface
 subroutine free_energyrec(an,bn2,nrec,trotter,ene_out, mult,&  
  &  prtvol,n_pt_integ,xmax,&  
  &  ene_out1,ene_out2,ene_out3,ene_out4)
  use defs_basis
  implicit none
  integer,intent(in) :: n_pt_integ
  integer,intent(in) :: nrec
  integer,intent(in) :: prtvol
  integer,intent(in) :: trotter
  real(dp),intent(out) :: ene_out
  real(dp),intent(out) :: ene_out1
  real(dp),intent(out) :: ene_out2
  real(dp),intent(out) :: ene_out3
  real(dp),intent(out) :: ene_out4
  real(dp),intent(in) :: mult
  real(dp),intent(in) :: xmax
  real(dp),intent(in) :: an(0:nrec)
  real(dp),intent(in) :: bn2(0:nrec)
 end subroutine free_energyrec
end interface

interface
 subroutine getngrec(ngfft,inf_rmet,ngfftrec,nfftrec,rtroncat)
  use defs_basis
  implicit none
  integer,intent(out) :: nfftrec
  real(dp),intent(in) :: rtroncat
  integer,intent(in) :: ngfft(18)
  integer,intent(out) :: ngfftrec(18)
  real(dp),intent(in) :: inf_rmet(3,3)
 end subroutine getngrec
end interface

interface
 subroutine green_kernel(ZT_p,inf_rmet,inf_ucvol,mult,mpi_enreg,ngfft,nfft,prtvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: nfft
  integer,intent(in) :: prtvol
  real(dp),intent(in) :: inf_ucvol
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: mult
  integer,intent(in) :: ngfft(18)
  real(dp),intent(out) :: ZT_p(1:2,0:nfft-1)
  real(dp),intent(in) :: inf_rmet(3,3)
 end subroutine green_kernel
end interface

interface
 subroutine recursion(exppot,coordx,coordy,coordz,an,bn2,&  
  &  rho_out,&  
  &  nrec,fermie,tsmear,trotter,&  
  &  ZT_p, tol,&  
  &  get_rec_coef,prtvol,&  
  &  mpi_enreg,nfft,ngfft,rmet,inf_ucvol,tim_fourdp)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: coordx
  integer,intent(in) :: coordy
  integer,intent(in) :: coordz
  integer,intent(in) :: get_rec_coef
  integer,intent(in) :: nfft
  integer,intent(in) :: nrec
  integer,intent(in) :: prtvol
  integer,intent(in) :: tim_fourdp
  integer,intent(in) :: trotter
  real(dp),intent(in) :: fermie
  real(dp),intent(in) :: inf_ucvol
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(out) :: rho_out
  real(dp),intent(in) :: tol
  real(dp),intent(in) :: tsmear
  integer,intent(in) :: ngfft(18)
  real(dp),intent(in) :: ZT_p(1:2,0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1)
  real(dp),intent(inout) :: an(0:nrec)
  real(dp),intent(inout) :: bn2(0:nrec)
  real(dp),intent(in) :: exppot(0:ngfft(1)-1,0:ngfft(2)-1,0:ngfft(3)-1)
  real(dp),intent(in) :: rmet(3,3)
 end subroutine recursion
end interface

interface
 subroutine vtorhorec(densymop_gs,dtfil,dtset,&  
  &  ek,enl,entropy,e_eigenvalues,fermie,&  
  &  grnl,irrzon,mpi_enreg,natom,nfftf,nspden,nsppol,nsym,phnons,&  
  &  rhog, rhor,ucvol, vtrial, rmet, quit,get_ek,get_entropy)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: get_ek
  integer,intent(in) :: get_entropy
  integer,intent(in) :: natom
  integer,intent(in) :: nfftf
  integer,intent(in) :: nspden
  integer,intent(in) :: nsppol
  integer,intent(in) :: nsym
  integer,intent(inout) :: quit
  type(dens_sym_operator_type),intent(in) :: densymop_gs
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(out) :: e_eigenvalues
  real(dp),intent(out) :: ek
  real(dp),intent(out) :: enl
  real(dp),intent(out) :: entropy
  real(dp),intent(out) :: fermie
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: ucvol
  real(dp),intent(out) :: grnl(3*natom)
  integer,intent(in) :: irrzon((dtset%nfft)**(1-1/nsym),2,nspden/nsppol)
  real(dp),intent(in) :: phnons(2,(dtset%nfft)**(1-1/nsym),nspden/nsppol)
  real(dp),intent(inout) :: rhog(2,nfftf)
  real(dp),intent(out) :: rhor(nfftf,dtset%nspden)
  real(dp),intent(in) :: rmet(3,3)
  real(dp),intent(in) :: vtrial(nfftf,nspden)
 end subroutine vtorhorec
end interface

end module interfaces_15recursion
!!***
