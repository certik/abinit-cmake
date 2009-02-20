!!****m* ABINIT/interfaces_21paral_md
!! NAME
!! interfaces_21paral_md
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/21paral_md
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

module interfaces_21paral_md

 implicit none

interface
 subroutine pstate(acell,codvsn,cpui,dtfil,dtset,iexit,mband,mgfft,mkmem,mpi_enreg,&  
  &  mpw,natom,nfft,nkpt,npwtot,nspden,nspinor,nsppol,nsym,occ,&  
  &  pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(inout) :: iexit
  integer, intent(in) :: mband
  integer, intent(in) :: mgfft
  integer, intent(in) :: mkmem
  integer, intent(in) :: mpw
  integer, intent(inout) :: natom
  integer, intent(in) :: nfft
  integer, intent(inout) :: nkpt
  integer, intent(inout) :: nspden
  integer, intent(inout) :: nspinor
  integer, intent(inout) :: nsppol
  integer, intent(inout) :: nsym
  character(len=6), intent(in) :: codvsn
  real(dp), intent(in) :: cpui
  type(datafiles_type), intent(inout) :: dtfil
  type(dataset_type), intent(inout) :: dtset
  type(mpi_type), intent(inout) :: mpi_enreg
  type(pseudopotential_type), intent(inout) :: psps
  type(results_gs_type), intent(out) :: results_gs
  real(dp), intent(in) :: walli
  real(dp), intent(inout) :: acell(3)
  integer, intent(out) :: npwtot(nkpt)
  real(dp), intent(inout) :: occ(mband*nkpt*nsppol)
  type(pawrad_type), intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
  type(pawtab_type), intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
  real(dp), intent(inout) :: rprim(3,3)
  real(dp), intent(inout) :: vel(3,natom)
  real(dp) :: xred(3,natom)
 end subroutine pstate
end interface

end module interfaces_21paral_md
!!***
