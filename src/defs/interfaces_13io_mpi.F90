!!****m* ABINIT/interfaces_13io_mpi
!! NAME
!! interfaces_13io_mpi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/13io_mpi
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

module interfaces_13io_mpi

 implicit none

interface
 subroutine chkexi(cpus,filnam,iexit,iout,mpi_enreg,openexit)
  use defs_basis
  use defs_datatypes
  implicit none
  integer          ,intent(out) :: iexit
  integer          ,intent(in) :: iout
  integer          ,intent(in) :: openexit
  real(dp)         ,intent(in) :: cpus
  character(len=fnlen),intent(in) :: filnam
  type(mpi_type)   ,intent(inout) :: mpi_enreg
 end subroutine chkexi
end interface

interface
 subroutine handle_ncerr(ncerr,message)
  implicit none
  integer         ,intent(in) :: ncerr
  character(len=*),intent(in) :: message
 end subroutine handle_ncerr
end interface

interface
 subroutine hdr_comm(hdr,master,me,spaceComm)
  use defs_datatypes
  implicit none
  integer, intent(in) :: master
  integer, intent(in) :: me
  integer, intent(in) :: spaceComm
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_comm
end interface


!Generic interface of the routines hdr_io
interface hdr_io
 subroutine hdr_io_wfftype(fform,hdr,rdwr,wff)
  use defs_datatypes
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  type(hdr_type),intent(inout) :: hdr
  type(wffile_type),intent(inout) :: wff
 end subroutine hdr_io_wfftype
 subroutine hdr_io_int(fform,hdr,rdwr,unitfi)
  use defs_datatypes
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  integer,intent(in) :: unitfi
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_io_int
end interface
!End of the generic interface of hdr_io


!Generic interface of the routines hdr_io_netcdf
interface hdr_io_netcdf
 subroutine hdr_io_netcdf_wfftype(fform,hdr,rdwr,wff)
  use defs_datatypes
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  type(hdr_type),intent(inout) :: hdr
  type(wffile_type),intent(inout) :: wff
 end subroutine hdr_io_netcdf_wfftype
 subroutine hdr_io_netcdf_int(fform,hdr,rdwr,unitfi)
  use defs_datatypes
  implicit none
  integer,intent(inout) :: fform
  integer,intent(in) :: rdwr
  integer,intent(in) :: unitfi
  type(hdr_type),intent(inout) :: hdr
 end subroutine hdr_io_netcdf_int
end interface
!End of the generic interface of hdr_io_netcdf


!Generic interface of the routines hdr_skip
interface hdr_skip
 subroutine hdr_skip_int(unitfi,ierr)
  implicit none
  integer, intent(out) :: ierr
  integer, intent(in) :: unitfi
 end subroutine hdr_skip_int
 subroutine hdr_skip_wfftype(wff,ierr)
  use defs_datatypes
  implicit none
  integer, intent(out) :: ierr
  type(wffile_type),intent(inout) :: wff
 end subroutine hdr_skip_wfftype
end interface
!End of the generic interface of hdr_skip

interface
 subroutine outxfhist(nxfh,natom,mxfh,xfhist,option,wff2,ios)
  use defs_basis
  use defs_datatypes
  implicit none
  integer          ,intent(out) :: ios
  integer          ,intent(in) :: mxfh
  integer          ,intent(in) :: natom
  integer          ,intent(inout) :: nxfh
  integer          ,intent(in) :: option
  type(wffile_type),intent(inout) :: wff2
  real(dp)         ,intent(inout) :: xfhist(3,natom+4,2,mxfh)
 end subroutine outxfhist
end interface

interface
 subroutine rwwf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,nband,&  
  &  nband_disk,npw,nspinor,occ,option,optkg,tim_rwwf,wff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(inout) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: optkg
  integer,intent(in) :: tim_rwwf
  type(mpi_type), intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(inout), target :: cg(2,mcg)
  real(dp),intent(inout), target :: eigen((2*mband)**formeig*mband)
  integer,intent(inout), target :: kg_k(3,optkg*npw)
  real(dp),intent(inout), target :: occ(mband)
 end subroutine rwwf
end interface

interface
 subroutine writewf(cg,eigen,formeig,headform,icg,ikpt,isppol,kg_k,mband,mcg,mpi_enreg,nband,nband_disk,&  
  &  npw,nspinor,occ,option,optkg,wff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: formeig
  integer,intent(in) :: headform
  integer,intent(in) :: icg
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(in) :: mband
  integer,intent(in) :: mcg
  integer,intent(in) :: nband
  integer,intent(in) :: nband_disk
  integer,intent(in) :: npw
  integer,intent(in) :: nspinor
  integer,intent(in) :: option
  integer,intent(in) :: optkg
  type(mpi_type),intent(inout) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in), target :: cg(2,mcg)
  real(dp),intent(in), target :: eigen((2*mband)**formeig*mband)
  integer,intent(in), target :: kg_k(3,optkg*npw)
  real(dp),intent(in), target :: occ(mband)
 end subroutine writewf
end interface

interface
 subroutine WffClose(wff,ier)
  use defs_datatypes
  implicit none
  integer, intent(out) :: ier
  type(wffile_type), intent(inout) :: wff
 end subroutine WffClose
end interface

interface
 subroutine WffDelete(wff,ier)
  use defs_datatypes
  implicit none
  integer, intent(out) :: ier
  type(wffile_type),intent(inout) :: wff
 end subroutine WffDelete
end interface

interface
 subroutine WffKg(wff,optkg)
  use defs_datatypes
  implicit none
  integer          ,intent(in) :: optkg
  type(wffile_type),intent(inout) :: wff
 end subroutine WffKg
end interface

interface
 subroutine WffOffset(wff,sender,spaceComm,ier)
  use defs_datatypes
  implicit none
  integer          ,intent(out) :: ier
  integer          ,intent(inout) :: sender
  integer          ,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
 end subroutine WffOffset
end interface

interface
 subroutine WffOpen(accesswff,spaceComm,filename,ier,wff,master,me,unwff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(in) :: accesswff
  integer, intent(out) :: ier
  integer, intent(in) :: master
  integer, intent(in) :: me
  integer, intent(in) :: spaceComm
  integer, intent(in) :: unwff
  character(len=fnlen), intent(in) :: filename
  type(wffile_type), intent(out) :: wff
 end subroutine WffOpen
end interface

interface
 subroutine WffReadDataRec(dparray,ierr,ndp,wff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(out) :: ierr
  integer, intent(in) :: ndp
  type(wffile_type), intent(inout) :: wff
  real(dp), intent(out) :: dparray(ndp)
 end subroutine WffReadDataRec
end interface

interface
 subroutine WffReadNpwRec(ierr,ikpt,isppol,nband_disk,npw,nspinor,wff)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: ikpt
  integer,intent(in) :: isppol
  integer,intent(out) :: nband_disk
  integer,intent(out) :: npw
  integer,intent(out) :: nspinor
  type(wffile_type),intent(inout) :: wff
 end subroutine WffReadNpwRec
end interface

interface
 subroutine WffReadSkipRec(ierr,nrec,wff)
  use defs_datatypes
  implicit none
  integer          ,intent(out) :: ierr
  integer          ,intent(in) :: nrec
  type(wffile_type),intent(inout) :: wff
 end subroutine WffReadSkipRec
end interface

interface
 subroutine wffwritecg(wff,cg,mcg, icg,nband_disk,npwso, spaceComm, ierr1)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(in) :: icg
  integer,intent(out) :: ierr1
  integer,intent(in) :: mcg
  integer,intent(in) :: nband_disk
  integer,intent(in) :: npwso
  integer, intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: cg(2,mcg)
 end subroutine wffwritecg
end interface

interface
 subroutine WffWriteDataRec(dparray,ierr,ndp,wff)
  use defs_basis
  use defs_datatypes
  implicit none
  integer, intent(out) :: ierr
  integer, intent(in) :: ndp
  type(wffile_type), intent(inout) :: wff
  real(dp), intent(in) :: dparray(ndp)
 end subroutine WffWriteDataRec
end interface

interface
 subroutine WffWriteDataRecInt(intarray,ierr,nn,wff)
  use defs_datatypes
  implicit none
  integer, intent(out) :: ierr
  integer, intent(in) :: nn
  type(wffile_type), intent(inout) :: wff
  integer, intent(in) :: intarray(nn)
 end subroutine WffWriteDataRecInt
end interface

interface
 subroutine WffWriteNpwRec(ierr,nband_disk,npw,nspinor,wff)
  use defs_datatypes
  implicit none
  integer, intent(out) :: ierr
  integer, intent(in) :: nband_disk
  integer, intent(in) :: npw
  integer, intent(in) :: nspinor
  type(wffile_type), intent(inout) :: wff
 end subroutine WffWriteNpwRec
end interface

interface
 subroutine WffWriteNpwRec_cs(ierr,mpi_enreg,nband_disk,npw,nspinor,wff)
  use defs_datatypes
  implicit none
  integer, intent(out) :: ierr
  integer, intent(in) :: nband_disk
  integer, intent(in) :: npw
  integer, intent(in) :: nspinor
  type(mpi_type), intent(inout) :: mpi_enreg
  type(wffile_type), intent(inout) :: wff
 end subroutine WffWriteNpwRec_cs
end interface

end module interfaces_13io_mpi
!!***
