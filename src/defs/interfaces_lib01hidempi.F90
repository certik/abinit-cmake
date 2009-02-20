!!****m* ABINIT/interfaces_lib01hidempi
!! NAME
!! interfaces_lib01hidempi
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/lib01hidempi
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

module interfaces_lib01hidempi

 implicit none


!Generic interface of the routines xallgather_mpi
interface xallgather_mpi
 subroutine xallgather_mpi_int(xval,recvcounts,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval
  integer,intent(inout) :: recvcounts(:)
 end subroutine xallgather_mpi_int
 subroutine xallgather_mpi_dp1d(xval,nelem,recvcounts,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: nelem
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvcounts(:)
  real(dp),intent(in) :: xval(:)
 end subroutine xallgather_mpi_dp1d
 subroutine xallgather_mpi_dp4d(xval,nelem,recvcounts,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: nelem
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvcounts(:,:,:,:)
  real(dp),intent(in) :: xval(:,:,:,:)
 end subroutine xallgather_mpi_dp4d
end interface
!End of the generic interface of xallgather_mpi


!Generic interface of the routines xallgatherv_mpi
interface xallgatherv_mpi
 subroutine xallgatherv_mpi_int2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(inout) :: recvbuf(:,:)
  integer,intent(in) :: recvcounts(:)
  integer,intent(in) :: xval(:,:)
 end subroutine xallgatherv_mpi_int2d
 subroutine xallgatherv_mpi_int(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(inout) :: recvbuf(:)
  integer,intent(in) :: recvcounts(:)
  integer,intent(in) :: xval(:)
 end subroutine xallgatherv_mpi_int
 subroutine xallgatherv_mpi_dp(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:)
  real(dp),intent(in) :: xval(:)
 end subroutine xallgatherv_mpi_dp
 subroutine xallgatherv_mpi_dp2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xallgatherv_mpi_dp2d
 subroutine xallgatherv_mpi_dp3d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:,:)
  real(dp),intent(in) :: xval(:,:,:)
 end subroutine xallgatherv_mpi_dp3d
 subroutine xallgatherv_mpi_dp4d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nelem
  integer,intent(in) :: spaceComm
  integer,intent(in) :: displs(:)
  integer,intent(in) :: recvcounts(:)
  real(dp),intent(inout) :: recvbuf(:,:,:,:)
  real(dp),intent(in) :: xval(:,:,:,:)
 end subroutine xallgatherv_mpi_dp4d
end interface
!End of the generic interface of xallgatherv_mpi


!Generic interface of the routines xalltoall_mpi
interface xalltoall_mpi
 subroutine xalltoall_mpi_dp2d(xval, sendsize, recvbuf, recvsize, spaceComm, ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: recvsize
  integer ,intent(in) :: sendsize
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: recvbuf(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xalltoall_mpi_dp2d
end interface
!End of the generic interface of xalltoall_mpi


!Generic interface of the routines xalltoallv_mpi
interface xalltoallv_mpi
 subroutine xalltoallv_mpi_dp2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: rdispls(:)
  integer ,intent(in) :: recvcnts(:)
  integer ,intent(in) :: sdispls(:)
  integer ,intent(in) :: sendcnts(:)
  real(dp),intent(inout) :: recvbuf(:,:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xalltoallv_mpi_dp2d
 subroutine xalltoallv_mpi_int2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(in) :: rdispls(:)
  integer,intent(inout) :: recvbuf(:,:)
  integer ,intent(in) :: recvcnts(:)
  integer ,intent(in) :: sdispls(:)
  integer ,intent(in) :: sendcnts(:)
  integer,intent(in) :: xval(:,:)
 end subroutine xalltoallv_mpi_int2d
end interface
!End of the generic interface of xalltoallv_mpi


!Generic interface of the routines xcast_mpi
interface xcast_mpi
 subroutine xcast_mpi_intv(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval
 end subroutine xcast_mpi_intv
 subroutine xcast_mpi_int1d(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:)
 end subroutine xcast_mpi_int1d
 subroutine xcast_mpi_int2d(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_int2d
 subroutine xcast_mpi_int3d(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_int3d
 subroutine xcast_mpi_dpv(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval
 end subroutine xcast_mpi_dpv
 subroutine xcast_mpi_dp1d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:)
 end subroutine xcast_mpi_dp1d
 subroutine xcast_mpi_dp2d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_dp2d
 subroutine xcast_mpi_dp3d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_dp3d
 subroutine xcast_mpi_dp4d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:)
 end subroutine xcast_mpi_dp4d
 subroutine xcast_mpi_spv(xval,master,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: master
  integer,intent(in) :: spaceComm
  real,intent(inout) :: xval
 end subroutine xcast_mpi_spv
 subroutine xcast_mpi_sp1d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real,intent(inout) :: xval(:)
 end subroutine xcast_mpi_sp1d
 subroutine xcast_mpi_sp2d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real,intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_sp2d
 subroutine xcast_mpi_sp3d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real,intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_sp3d
 subroutine xcast_mpi_sp4d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real,intent(inout) :: xval(:,:,:,:)
 end subroutine xcast_mpi_sp4d
 subroutine xcast_mpi_cplxv(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval
 end subroutine xcast_mpi_cplxv
 subroutine xcast_mpi_cplx1d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:)
 end subroutine xcast_mpi_cplx1d
 subroutine xcast_mpi_cplx2d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_cplx2d
 subroutine xcast_mpi_cplx3d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_cplx3d
 subroutine xcast_mpi_cplx4d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:,:)
 end subroutine xcast_mpi_cplx4d
 subroutine xcast_mpi_dcv(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dp),intent(inout) :: xval
 end subroutine xcast_mpi_dcv
 subroutine xcast_mpi_dc1d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dp),intent(inout) :: xval(:)
 end subroutine xcast_mpi_dc1d
 subroutine xcast_mpi_dc2d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dp),intent(inout) :: xval(:,:)
 end subroutine xcast_mpi_dc2d
 subroutine xcast_mpi_dc3d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dp),intent(inout) :: xval(:,:,:)
 end subroutine xcast_mpi_dc3d
end interface
!End of the generic interface of xcast_mpi

interface
 subroutine xcomm_world(mpi_enreg,spaceComm)
  use defs_datatypes
  implicit none
  integer,intent(out) :: spaceComm
  type(mpi_type) :: mpi_enreg
 end subroutine xcomm_world
end interface

interface
 subroutine xcomm_init(mpi_enreg,spaceComm)
  use defs_datatypes
  implicit none
  integer,intent(out) :: spaceComm
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xcomm_init
end interface

interface
 subroutine xmaster_init(mpi_enreg,master)
  use defs_datatypes
  implicit none
  integer,intent(out) :: master
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xmaster_init
end interface

interface
 subroutine xmaster_init_fft(mpi_enreg,master)
  use defs_datatypes
  implicit none
  integer,intent(out) :: master
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xmaster_init_fft
end interface

interface
 subroutine xme_init(mpi_enreg,me)
  use defs_datatypes
  implicit none
  integer,intent(out) :: me
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xme_init
end interface

interface
 subroutine xproc_init(mpi_enreg,nproc_max)
  use defs_datatypes
  implicit none
  integer,intent(out) :: nproc_max
  type(mpi_type),intent(in) :: mpi_enreg
 end subroutine xproc_init
end interface

interface
 subroutine xme_whoiam(me)
  implicit none
  integer,intent(out) :: me
 end subroutine xme_whoiam
end interface

interface
 subroutine xproc_max(nproc,ierr)
  implicit none
  integer,intent(out) :: ierr
  integer,intent(out) :: nproc
 end subroutine xproc_max
end interface

interface
 subroutine xdefineOff(formeig,wff,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)
  use defs_datatypes
  implicit none
  integer, intent(in) :: formeig
  integer, intent(in) :: nkpt
  integer, intent(in) :: nspinor
  integer, intent(in) :: nsppol
  type(mpi_type),intent(in) :: mpi_enreg
  type(wffile_type),intent(inout) :: wff
  integer, intent(in) :: nband(nkpt*nsppol)
  integer, intent(in) :: npwarr(nkpt)
 end subroutine xdefineOff
end interface

interface
 subroutine xderiveWRecEnd(wff,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveWRecEnd
end interface

interface
 subroutine xderiveWRecEnd_cs(wff,ierr,me)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: me
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveWRecEnd_cs
end interface

interface
 subroutine xderiveWRecInit(wff,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveWRecInit
end interface

interface
 subroutine xderiveWRecInit_cs(wff,ierr,me)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: me
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveWRecInit_cs
end interface

interface
 subroutine xderiveRRecEnd(wff,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveRRecEnd
end interface

interface
 subroutine xderiveRRecInit(wff,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveRRecInit
end interface


!Generic interface of the routines xderiveread
interface xderiveread
 subroutine xderiveRead_int(wff,xval,n1,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  type(wffile_type),intent(inout) :: wff
  integer,intent(out) :: xval(:)
 end subroutine xderiveRead_int
 subroutine xderiveRead_int_mpio(wff,xval,n1,ierr,spaceComm)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(out) :: xval(:)
 end subroutine xderiveRead_int_mpio
 subroutine xderiveRead_int2d(wff,xval,n1,n2,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  integer,intent(out) :: xval(:,:)
 end subroutine xderiveRead_int2d
 subroutine xderiveRead_int2d_mpio(wff,xval,n1,n2,ierr,spaceComm)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(out) :: xval(:,:)
 end subroutine xderiveRead_int2d_mpio
 subroutine xderiveRead_dp(wff,xval,n1,ierr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: xval(:)
 end subroutine xderiveRead_dp
 subroutine xderiveRead_dp_mpio(wff,xval,n1,ierr,spaceComm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: xval(:)
 end subroutine xderiveRead_dp_mpio
 subroutine xderiveRead_dp2d(wff,xval,n1,n2,ierr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: xval(:,:)
 end subroutine xderiveRead_dp2d
 subroutine xderiveRead_dp2d_mpio(wff,xval,n1,n2,ierr,spaceComm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: xval(:,:)
 end subroutine xderiveRead_dp2d_mpio
end interface
!End of the generic interface of xderiveread


!Generic interface of the routines xderivereadval
interface xderivereadval
 subroutine xderiveReadVal_dp(wff,xval)
  use defs_basis
  use defs_datatypes
  implicit none
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(out) :: xval
 end subroutine xderiveReadVal_dp
 subroutine xderiveReadVal_int(wff,xval)
  use defs_datatypes
  implicit none
  integer,intent(out) :: xval
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveReadVal_int
 subroutine xderiveReadVal_char(wff,xval,n)
  use defs_datatypes
  implicit none
  integer,intent(in) :: n
  type(wffile_type),intent(inout) :: wff
  character(len=*),intent(out) :: xval
 end subroutine xderiveReadVal_char
end interface
!End of the generic interface of xderivereadval


!Generic interface of the routines xderivewrite
interface xderivewrite
 subroutine xderiveWrite_int(wff,xval,n1,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: xval(:)
 end subroutine xderiveWrite_int
 subroutine xderiveWrite_int_mpio(wff,xval,n1,ierr,spaceComm)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: xval(:)
 end subroutine xderiveWrite_int_mpio
 subroutine xderiveWrite_int2d(wff,xval,n1,n2,ierr)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: xval(:,:)
 end subroutine xderiveWrite_int2d
 subroutine xderiveWrite_int2d_mpio(wff,xval,n1,n2,ierr,spaceComm)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: xval(:,:)
 end subroutine xderiveWrite_int2d_mpio
 subroutine xderiveWrite_int2d_mpio_arr(wff,xval,n1,n2,ierr,local_offset)
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: local_offset(:)
  integer,intent(in) :: xval(:,:)
 end subroutine xderiveWrite_int2d_mpio_arr
 subroutine xderiveWrite_dp(wff,xval,n1,ierr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: xval(:)
 end subroutine xderiveWrite_dp
 subroutine xderiveWrite_dp_mpio(wff,xval,n1,ierr,spaceComm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: xval(:)
 end subroutine xderiveWrite_dp_mpio
 subroutine xderiveWrite_dp2d_mpio_arr(wff,xval,n1,n2,ierr,local_offset)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  integer,intent(in) :: local_offset(:)
  real(dp),intent(in) :: xval(:,:)
 end subroutine xderiveWrite_dp2d_mpio_arr
 subroutine xderiveWrite_dp2d(wff,xval,n1,n2,ierr)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: xval(:,:)
 end subroutine xderiveWrite_dp2d
 subroutine xderiveWrite_dp2d_mpio(wff,xval,n1,n2,ierr,spaceComm)
  use defs_basis
  use defs_datatypes
  implicit none
  integer,intent(out) :: ierr
  integer,intent(in) :: n1
  integer,intent(in) :: n2
  integer,intent(in) :: spaceComm
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: xval(:,:)
 end subroutine xderiveWrite_dp2d_mpio
end interface
!End of the generic interface of xderivewrite


!Generic interface of the routines xderivewriteval
interface xderivewriteval
 subroutine xderiveWriteVal_dp(wff,xval)
  use defs_basis
  use defs_datatypes
  implicit none
  type(wffile_type),intent(inout) :: wff
  real(dp),intent(in) :: xval
 end subroutine xderiveWriteVal_dp
 subroutine xderiveWriteVal_int(wff,xval)
  use defs_datatypes
  implicit none
  integer,intent(in) :: xval
  type(wffile_type),intent(inout) :: wff
 end subroutine xderiveWriteVal_int
 subroutine xderiveWriteVal_char(wff,xval,n)
  use defs_datatypes
  implicit none
  integer,intent(in) :: n
  type(wffile_type),intent(inout) :: wff
  character(len=*),intent(in) :: xval
 end subroutine xderiveWriteVal_char
end interface
!End of the generic interface of xderivewriteval


!Generic interface of the routines xexch_mpi
interface xexch_mpi
 subroutine xexch_mpi_intn(vsend,n1,sender,vrecv,recever,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: n1
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: vrecv(:)
  integer,intent(in) :: vsend(:)
 end subroutine xexch_mpi_intn
 subroutine xexch_mpi_int2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nt
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: vrecv(:,:)
  integer,intent(in) :: vsend(:,:)
 end subroutine xexch_mpi_int2d
 subroutine xexch_mpi_dpn(vsend,n1,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: n1
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: vrecv(:)
  real(dp),intent(in) :: vsend(:)
 end subroutine xexch_mpi_dpn
 subroutine xexch_mpi_dp2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nt
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: vrecv(:,:)
  real(dp),intent(in) :: vsend(:,:)
 end subroutine xexch_mpi_dp2d
 subroutine xexch_mpi_dp3d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: nt
  integer,intent(in) :: recever
  integer,intent(in) :: sender
  integer,intent(in) :: spaceComm
  real(dp),intent(inout) :: vrecv(:,:,:)
  real(dp),intent(in) :: vsend(:,:,:)
 end subroutine xexch_mpi_dp3d
end interface
!End of the generic interface of xexch_mpi


!Generic interface of the routines xmax_mpi
interface xmax_mpi
 subroutine xmax_mpi_intv(xval,xmax,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(inout) :: xmax
  integer ,intent(inout) :: xval
 end subroutine xmax_mpi_intv
 subroutine xmax_mpi_dpv(xval,xmax,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xmax
  real(dp),intent(inout) :: xval
 end subroutine xmax_mpi_dpv
end interface
!End of the generic interface of xmax_mpi


!Generic interface of the routines xmin_mpi
interface xmin_mpi
 subroutine xmin_mpi_intv(xval,xmin,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  integer ,intent(inout) :: xmin
  integer ,intent(inout) :: xval
 end subroutine xmin_mpi_intv
 subroutine xmin_mpi_dpv(xval,xmin,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xmin
  real(dp),intent(inout) :: xval
 end subroutine xmin_mpi_dpv
end interface
!End of the generic interface of xmin_mpi


!Generic interface of the routines xsum_master
interface xsum_master
 subroutine xsum_master_dp1d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:)
 end subroutine xsum_master_dp1d
 subroutine xsum_master_dp2d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:)
 end subroutine xsum_master_dp2d
 subroutine xsum_master_dp3d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xsum_master_dp3d
 subroutine xsum_master_dp4d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_master_dp4d
 subroutine xsum_master_dp5d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:)
 end subroutine xsum_master_dp5d
 subroutine xsum_master_dp6d(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:,:)
 end subroutine xsum_master_dp6d
 subroutine xsum_master_int4d(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  integer ,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_master_int4d
 subroutine xsum_master_c2cplx(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex ,intent(inout) :: xval(:,:)
 end subroutine xsum_master_c2cplx
 subroutine xsum_master_c3cplx(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex ,intent(inout) :: xval(:,:,:)
 end subroutine xsum_master_c3cplx
 subroutine xsum_master_c4cplx(xval,master,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex ,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_master_c4cplx
 subroutine xsum_master_c2dpc(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc) ,intent(inout) :: xval(:,:)
 end subroutine xsum_master_c2dpc
 subroutine xsum_master_c3dpc(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc) ,intent(inout) :: xval(:,:,:)
 end subroutine xsum_master_c3dpc
 subroutine xsum_master_c4dpc(xval,master,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: master
  integer ,intent(in) :: spaceComm
  complex(dpc) ,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_master_c4dpc
end interface
!End of the generic interface of xsum_master


!Generic interface of the routines xsum_mpi
interface xsum_mpi
 subroutine xsum_mpi_int(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:)
 end subroutine xsum_mpi_int
 subroutine xsum_mpi_intv(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval
 end subroutine xsum_mpi_intv
 subroutine xsum_mpi_intv2(xval,xsum,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xsum
  integer,intent(inout) :: xval
 end subroutine xsum_mpi_intv2
 subroutine xsum_mpi_intn(xval,n1,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: n1
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:)
 end subroutine xsum_mpi_intn
 subroutine xsum_mpi_int2t(xval,xsum,n1,spaceComm,ier)
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  integer ,intent(inout) :: xsum(:)
  integer ,intent(inout) :: xval(:)
 end subroutine xsum_mpi_int2t
 subroutine xsum_mpi_int2d(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_int2d
 subroutine xsum_mpi_int3d(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_int3d
 subroutine xsum_mpi_int4d(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  integer,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_int4d
 subroutine xsum_mpi_dp(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:)
 end subroutine xsum_mpi_dp
 subroutine xsum_mpi_dpvt(xval,xsum,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(out) :: xsum
  real(dp),intent(in) :: xval
 end subroutine xsum_mpi_dpvt
 subroutine xsum_mpi_dpv(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval
 end subroutine xsum_mpi_dpv
 subroutine xsum_mpi_dpn(xval,n1,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:)
 end subroutine xsum_mpi_dpn
 subroutine xsum_mpi_dp2d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_dp2d
 subroutine xsum_mpi_dp3d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_dp3d
 subroutine xsum_mpi_dp4d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_dp4d
 subroutine xsum_mpi_dp5d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:)
 end subroutine xsum_mpi_dp5d
 subroutine xsum_mpi_dp6d(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xval(:,:,:,:,:,:)
 end subroutine xsum_mpi_dp6d
 subroutine xsum_mpi_dp2t(xval,xsum,n1,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xsum(:)
  real(dp),intent(inout) :: xval(:)
 end subroutine xsum_mpi_dp2t
 subroutine xsum_mpi_dp3d2t(xval,xsum,n1,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xsum(:,:,:)
  real(dp),intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_dp3d2t
 subroutine xsum_mpi_dp4d2t(xval,xsum,n1,spaceComm,ier)
  use defs_basis
  implicit none
  integer ,intent(out) :: ier
  integer ,intent(in) :: n1
  integer ,intent(in) :: spaceComm
  real(dp),intent(inout) :: xsum(:,:,:,:)
  real(dp),intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_dp4d2t
 subroutine xsum_mpi_c2dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_c2dc
 subroutine xsum_mpi_c3dc(xval,spaceComm,ier)
  use defs_basis
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex(dpc),intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_c3dc
 subroutine xsum_mpi_c1cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:)
 end subroutine xsum_mpi_c1cplx
 subroutine xsum_mpi_c2cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:)
 end subroutine xsum_mpi_c2cplx
 subroutine xsum_mpi_c3cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:)
 end subroutine xsum_mpi_c3cplx
 subroutine xsum_mpi_c4cplx(xval,spaceComm,ier)
  implicit none
  integer,intent(out) :: ier
  integer,intent(in) :: spaceComm
  complex,intent(inout) :: xval(:,:,:,:)
 end subroutine xsum_mpi_c4cplx
end interface
!End of the generic interface of xsum_mpi

end module interfaces_lib01hidempi
!!***
