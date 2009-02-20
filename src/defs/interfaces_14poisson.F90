!!****m* ABINIT/interfaces_14poisson
!! NAME
!! interfaces_14poisson
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/14poisson
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

module interfaces_14poisson

 implicit none

interface
 subroutine Psolver_hartree(dtset, enhartr, mpi_enreg, rhor, rprimd, vhartr)
  use defs_basis
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  real(dp), intent(out) :: enhartr
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp),intent(in) :: rhor(dtset%nfft,dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vhartr(dtset%nfft)
 end subroutine Psolver_hartree
end interface

interface
 subroutine PSolver_kernel(dtset, iaction, kernel, mpi_enreg, rprimd)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: iaction
  type(dataset_type),intent(in) :: dtset
  type(mpi_type),intent(in) :: mpi_enreg
  real(dp), pointer :: kernel(:)
  real(dp),intent(in) :: rprimd(3,3)
 end subroutine PSolver_kernel
end interface

interface
 subroutine Psolver_rhohxc(dtset, enhartr, enxc, envxc, mpi_enreg, rhor, rprimd,&  
  &  vhartr, vxc, vxcavg)
  use defs_basis
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  real(dp), intent(out) :: enhartr
  real(dp), intent(out) :: envxc
  real(dp), intent(out) :: enxc
  type(mpi_type),intent(inout) :: mpi_enreg
  real(dp), intent(out) :: vxcavg
  real(dp),intent(inout) :: rhor(dtset%nfft, dtset%nspden)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(out) :: vhartr(dtset%nfft)
  real(dp),intent(out) :: vxc(dtset%nfft, dtset%nspden)
 end subroutine Psolver_rhohxc
end interface

end module interfaces_14poisson
!!***
