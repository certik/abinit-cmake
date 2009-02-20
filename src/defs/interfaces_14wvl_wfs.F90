!!****m* ABINIT/interfaces_14wvl_wfs
!! NAME
!! interfaces_14wvl_wfs
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/14wvl_wfs
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

module interfaces_14wvl_wfs

 implicit none

interface
 subroutine wvl_free_type_wfs(wfs)
  use defs_wvltypes
  implicit none
  type(wvl_wf_type),intent(inout) :: wfs
 end subroutine wvl_free_type_wfs
end interface

interface
 subroutine wvl_free_type_proj(proj)
  use defs_wvltypes
  implicit none
  type(wvl_projectors_type),intent(inout) :: proj
 end subroutine wvl_free_type_proj
end interface

interface
 subroutine wvl_init_type_proj(dtset, mpi_enreg, proj, psps, rprimd, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  type(wvl_projectors_type),intent(inout) :: proj
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine wvl_init_type_proj
end interface

interface
 subroutine wvl_init_type_wfs(dtset, mpi_enreg, psps, rprimd, wfs, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_wf_type),intent(inout) :: wfs
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(inout) :: xred(3,dtset%natom)
 end subroutine wvl_init_type_wfs
end interface

interface
 subroutine wvl_nl_gradient(dtset, grnl, mpi_enreg, occ, psps, rprimd, wvl, xcart)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(inout) :: grnl(3,dtset%natom)
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine wvl_nl_gradient
end interface

interface
 subroutine wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(hdr_type), intent(in) :: hdr
  type(hdr_type), intent(in) :: hdr0
  type(mpi_type), intent(in) :: mpi_enreg
  type(wffile_type),intent(in) :: wff
  type(wvl_wf_type), intent(inout) :: wfs
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(in) :: xred(3, dtset%natom)
 end subroutine wvl_read
end interface

interface
 subroutine wvl_write(dtset, eigen, mpi_enreg, option, rprimd, wff, wfs, xred)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer, intent(in) :: option
  type(dataset_type), intent(in) :: dtset
  type(mpi_type), intent(in) :: mpi_enreg
  type(wffile_type),intent(in) :: wff
  type(wvl_wf_type), intent(in) :: wfs
  real(dp), intent(in), target :: eigen(dtset%mband)
  real(dp), intent(in) :: rprimd(3, 3)
  real(dp), intent(inout) :: xred(3, dtset%natom)
 end subroutine wvl_write
end interface

interface
 subroutine wvl_setngfft(dtset, mpi_enreg)
  use defs_datatypes
  implicit none
  type(dataset_type),intent(inout) :: dtset
  type(mpi_type),intent(inout) :: mpi_enreg
 end subroutine wvl_setngfft
end interface

interface
 subroutine wvl_tail_corrections(dtset, energies, etotal, mpi_enreg, occ, psps,&  
  &  vtrial, wvl, xcart)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  type(dataset_type),intent(inout) :: dtset
  type(energies_type),intent(inout) :: energies
  real(dp),intent(out) :: etotal
  type(mpi_type),intent(inout) :: mpi_enreg
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_data),intent(inout) :: wvl
  real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
  real(dp),intent(in),target :: vtrial(dtset%nfft)
  real(dp),intent(in) :: xcart(3,dtset%natom)
 end subroutine wvl_tail_corrections
end interface

end module interfaces_14wvl_wfs
!!***
