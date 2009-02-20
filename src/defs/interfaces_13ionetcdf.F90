!!****m* ABINIT/interfaces_13ionetcdf
!! NAME
!! interfaces_13ionetcdf
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/13ionetcdf
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

module interfaces_13ionetcdf

 implicit none

interface
 subroutine abi_etsf_electrons_put(dtset, filapp)
  use defs_basis
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  character(len=fnlen),intent(in) :: filapp
 end subroutine abi_etsf_electrons_put
end interface

interface
 subroutine abi_etsf_geo_put(dtset, filapp, psps, rprimd, xred)
  use defs_basis
  use defs_datatypes
  implicit none
  type(dataset_type),intent(in) :: dtset
  character(len=fnlen),intent(in) :: filapp
  type(pseudopotential_type),intent(in) :: psps
  real(dp),intent(in),target :: rprimd(3,3)
  real(dp),intent(in),target :: xred(3,dtset%natom)
 end subroutine abi_etsf_geo_put
end interface

interface
 subroutine abi_etsf_init(dtset, filapp, itype, kdep, lmn_size, psps, wfs)
  use defs_basis
  use defs_datatypes
  use defs_wvltypes
  implicit none
  integer,intent(in) :: itype
  type(dataset_type),intent(in) :: dtset
  character(len=fnlen),intent(in) :: filapp
  logical,intent(in) :: kdep
  type(pseudopotential_type),intent(in) :: psps
  type(wvl_wf_type),intent(in) :: wfs
  integer,intent(in) :: lmn_size(psps%npsp)
 end subroutine abi_etsf_init
end interface

interface
 subroutine handle_err_netcdf(status)
  implicit none
  integer,intent(in) :: status
 end subroutine handle_err_netcdf
end interface

interface
 subroutine hdr_io_etsf(fform,hdr,rdwr,unitwff)
  use defs_datatypes
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  integer,intent(in) :: unitwff
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_io_etsf
end interface

interface
 subroutine ini_wf_etsf(dtset, lmn_size, npsp, ntypat, unwff)
  use defs_datatypes
  implicit none
  integer,intent(in) :: npsp
  integer,intent(in) :: ntypat
  integer,intent(in) :: unwff
  type(dataset_type),intent(in) :: dtset
  integer,intent(in) :: lmn_size(npsp)
 end subroutine ini_wf_etsf
end interface

interface
 subroutine ini_wf_netcdf(mpw,ncid_hdr,response)
  implicit none
  integer,intent(in) :: mpw
  integer,intent(in) :: ncid_hdr
  integer,intent(in) :: response
 end subroutine ini_wf_netcdf
end interface

interface
 subroutine write_header_moldynnetcdf(dtfil, dtset, natom, ncoord, nelt_strten)
  use defs_datatypes
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: ncoord
  integer,intent(in) :: nelt_strten
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
 end subroutine write_header_moldynnetcdf
end interface

interface
 subroutine write_moldynvaluenetcdf(amass, itime1, dtfil, dtset, Epot, Ekin,  nbat,&  
  &  nbdir, nb1, pos, cel, stress, rprimd, ucvol)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: itime1
  integer,intent(in) :: nb1
  integer,intent(in) :: nbat
  integer,intent(in) :: nbdir
  real(dp),intent(in) :: Ekin
  real(dp),intent(in) :: Epot
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  real(dp),intent(in) :: ucvol
  real(dp),intent(in) :: amass(nbat)
  real(dp),intent(in) :: cel(nbdir,nbat)
  real(dp),intent(in) :: pos(nbdir,nbat)
  real(dp),intent(in) :: rprimd(3,3)
  real(dp),intent(in) :: stress(nb1)
 end subroutine write_moldynvaluenetcdf
end interface

end module interfaces_13ionetcdf
!!***
