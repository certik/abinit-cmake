!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_xfuncmpi
!! NAME
!! defs_xfuncmpi
!!
!! FUNCTION
!! This module contains interfaces for using several procedures with the
!! same generic name, linked with MPI.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (MBoulet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 module defs_xfuncmpi

 implicit none

!Generic interface of the routines xsum_mpi
 interface xsum_mpi

  subroutine xsum_mpi_int(xval,spaceComm,ier)
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   integer,intent(inout) :: xval(:)
  end subroutine xsum_mpi_int

  subroutine xsum_mpi_intv(xval,spaceComm,ier)
   integer,intent(inout) :: xval
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
  end subroutine xsum_mpi_intv

  subroutine xsum_mpi_intv2(xval,xsum,spaceComm,ier)
   integer,intent(inout) :: xval
   integer,intent(inout) :: xsum
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
  end subroutine xsum_mpi_intv2

  subroutine xsum_mpi_intn(xval,n1,spaceComm,ier)
   integer,intent(in) :: n1
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   integer,intent(inout) :: xval(:)
  end subroutine xsum_mpi_intn

  subroutine xsum_mpi_int2d(xval,spaceComm,ier)
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   integer,intent(inout) :: xval(:,:)
  end subroutine xsum_mpi_int2d

  subroutine xsum_mpi_int3d(xval,spaceComm,ier)
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   integer,intent(inout) :: xval(:,:,:)
  end subroutine xsum_mpi_int3d

  subroutine xsum_mpi_dp(xval,spaceComm,ier)
   use defs_basis
   real(dp),intent(inout) :: xval(:)
   integer ,intent(in) :: spaceComm
   integer ,intent(out)   :: ier
  end subroutine xsum_mpi_dp

  subroutine xsum_mpi_dpv(xval,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval
  end subroutine xsum_mpi_dpv

  subroutine xsum_mpi_dpn(xval,n1,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: n1
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:)
  end subroutine xsum_mpi_dpn

  subroutine xsum_mpi_dp2d(xval,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:)
  end subroutine xsum_mpi_dp2d

  subroutine xsum_mpi_dp3d(xval,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:)
  end subroutine xsum_mpi_dp3d

  subroutine xsum_mpi_dp4d(xval,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:,:)
  end subroutine xsum_mpi_dp4d

  subroutine xsum_mpi_dp5d(xval,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:,:,:)
  end subroutine xsum_mpi_dp5d

  subroutine xsum_mpi_dp6d(xval,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:,:,:,:)
  end subroutine xsum_mpi_dp6d

  subroutine xsum_mpi_dp2t(xval,xsum,n1,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: n1
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:)
   real(dp),intent(inout) :: xsum(:)
  end subroutine xsum_mpi_dp2t

  subroutine xsum_mpi_dp3d2t(xval,xsum,n1,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: n1
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:)
   real(dp),intent(inout) :: xsum(:,:,:)
  end subroutine xsum_mpi_dp3d2t

  subroutine xsum_mpi_dp4d2t(xval,xsum,n1,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: n1
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:,:)
   real(dp),intent(inout) :: xsum(:,:,:,:)
  end subroutine xsum_mpi_dp4d2t

  subroutine xsum_mpi_c2dc(xval,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   complex(dpc),intent(inout) :: xval(:,:)
  end subroutine xsum_mpi_c2dc

  subroutine xsum_mpi_c3dc(xval,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   complex(dpc),intent(inout) :: xval(:,:,:)
  end subroutine xsum_mpi_c3dc

  subroutine xsum_mpi_c2cplx(xval,spaceComm,ier)
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   complex,intent(inout) :: xval(:,:)
  end subroutine xsum_mpi_c2cplx

  subroutine xsum_mpi_c3cplx(xval,spaceComm,ier)
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   complex,intent(inout) :: xval(:,:,:)
  end subroutine xsum_mpi_c3cplx

 end interface
!End of the generic interface of xsum_mpi


!Generic interface of the routines xsum_master
 interface xsum_master
 
  subroutine xsum_master_int4d(xval,master,spaceComm,ier)
   integer ,intent(in) :: master
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   integer,intent(inout) :: xval(:,:,:,:)
  end subroutine xsum_master_int4d

  subroutine xsum_master_dp2d(xval,master,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: master
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:)
  end subroutine xsum_master_dp2d

  subroutine xsum_master_dp3d(xval,master,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: master
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:)
  end subroutine xsum_master_dp3d

  subroutine xsum_master_dp4d(xval,master,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: master
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:,:)
  end subroutine xsum_master_dp4d

  subroutine xsum_master_dp5d(xval,master,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: master
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval(:,:,:,:,:)
  end subroutine xsum_master_dp5d

  subroutine xsum_master_c3cplx(xval,master,spaceComm,ier)
   integer ,intent(in) :: master
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   complex ,intent(inout) :: xval(:,:,:)
  end subroutine xsum_master_c3cplx
 end interface
!End of the generic interface of xsum_master



!Generic interface of the routines xmax_mpi
 interface xmax_mpi

  subroutine xmax_mpi_intv(xval,xmax,spaceComm,ier)
   integer ,intent(inout) :: xval
   integer ,intent(inout) :: xmax
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
  end subroutine xmax_mpi_intv

  subroutine xmax_mpi_dpv(xval,xmax,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval
   real(dp),intent(inout) :: xmax
  end subroutine xmax_mpi_dpv

 end interface
!End of the generic interface of xmax_mpi


!Generic interface of the routines xmin_mpi
 interface xmin_mpi

  subroutine xmin_mpi_intv(xval,xmin,spaceComm,ier)
   integer ,intent(inout) :: xval
   integer ,intent(inout) :: xmin
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
  end subroutine xmin_mpi_intv

  subroutine xmin_mpi_dpv(xval,xmin,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(inout) :: xval
   real(dp),intent(inout) :: xmin
  end subroutine xmin_mpi_dpv

 end interface
!End of the generic interface of xmin_mpi

!Generic interface of the routines xallgather_mpi
 interface xallgather_mpi

  subroutine xallgather_mpi_int(xval,recvcounts,spaceComm,ier)
   integer,intent(inout) :: xval
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   integer,intent(inout) :: recvcounts(:)
  end subroutine xallgather_mpi_int

  subroutine xallgather_mpi_dp4d(xval,nelem,recvcounts,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: nelem
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(in) :: xval(:,:,:,:)
   real(dp),intent(inout) :: recvcounts(:,:,:,:)
  end subroutine xallgather_mpi_dp4d

 end interface
!End of the generic interface of xallgather_mpi


!Generic interface of the routines xallgatherv_mpi
 interface xallgatherv_mpi

  subroutine xallgatherv_mpi_int2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
   integer,intent(in) :: nelem
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   integer,intent(in) :: xval(:,:)
   integer,intent(inout) :: recvbuf(:,:)
   integer,intent(in) :: recvcounts(:)
   integer,intent(in) :: displs(:)
  end subroutine xallgatherv_mpi_int2d

  subroutine xallgatherv_mpi_int(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
   integer,intent(in) :: nelem
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   integer,intent(in) :: xval(:)
   integer,intent(inout) :: recvbuf(:)
   integer,intent(in) :: recvcounts(:)
   integer,intent(in) :: displs(:)
  end subroutine xallgatherv_mpi_int

  subroutine xallgatherv_mpi_dp(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: nelem
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   real(dp),intent(in) :: xval(:)
   real(dp),intent(inout) :: recvbuf(:)
   integer,intent(in) :: recvcounts(:)
   integer,intent(in) :: displs(:)
  end subroutine xallgatherv_mpi_dp

  subroutine xallgatherv_mpi_dp2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: nelem
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   real(dp),intent(in) :: xval(:,:)
   real(dp),intent(inout) :: recvbuf(:,:)
   integer,intent(in) :: recvcounts(:)
   integer,intent(in) :: displs(:)
  end subroutine xallgatherv_mpi_dp2d

  subroutine xallgatherv_mpi_dp3d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: nelem
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   real(dp),intent(in) :: xval(:,:,:)
   real(dp),intent(inout) :: recvbuf(:,:,:)
   integer,intent(in) :: recvcounts(:)
   integer,intent(in) :: displs(:)
  end subroutine xallgatherv_mpi_dp3d

  subroutine xallgatherv_mpi_dp4d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: nelem
   integer,intent(in) :: spaceComm
   integer,intent(out) :: ier
   real(dp),intent(in) :: xval(:,:,:,:)
   real(dp),intent(inout) :: recvbuf(:,:,:,:)
   integer,intent(in) :: recvcounts(:)
   integer,intent(in) :: displs(:)
  end subroutine xallgatherv_mpi_dp4d

 end interface
!End of the generic interface of xallgatherv_mpi

!Generic interface of the routines xalltoallv_mpi
 interface xalltoallv_mpi

  subroutine xalltoallv_mpi_dp2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(in) :: xval(:,:)
   integer ,intent(in) :: sendcnts(:)
   integer ,intent(in) :: sdispls(:)
   real(dp),intent(inout) :: recvbuf(:,:)
   integer ,intent(in) :: recvcnts(:)
   integer ,intent(in) :: rdispls(:)
  end subroutine xalltoallv_mpi_dp2d

 end interface
!End of the generic interface of xalltoallv_mpi

!Generic interface of the routines xalltoall_mpi
 interface xalltoall_mpi

  subroutine xalltoall_mpi_dp2d(xval,sendsize,recvbuf,recvsize,spaceComm,ier)
   use defs_basis
   integer ,intent(in) :: spaceComm
   integer ,intent(out) :: ier
   real(dp),intent(in) :: xval(:,:)
   integer ,intent(in) :: sendsize
   real(dp),intent(inout) :: recvbuf(:,:)
   integer ,intent(in) :: recvsize
  end subroutine xalltoall_mpi_dp2d

 end interface
!End of the generic interface of xalltoall_mpi

!Generic interface of the routines xexch_mpi
 interface xexch_mpi

  subroutine xexch_mpi_intn(vsend,n1,sender,vrecv,recever,spaceComm,ier)
   integer,intent(in) :: n1
   integer,intent(in) :: vsend(:)
   integer,intent(inout) :: vrecv(:)
   integer,intent(in) :: sender,recever,spaceComm
   integer,intent(out)   :: ier
  end subroutine xexch_mpi_intn

  subroutine xexch_mpi_int2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
   integer,intent(in) :: nt
   integer,intent(in) :: vsend(:,:)
   integer,intent(inout) :: vrecv(:,:)
   integer,intent(in) :: sender,recever,spaceComm
   integer,intent(out)   :: ier
  end subroutine xexch_mpi_int2d

  subroutine xexch_mpi_dpn(vsend,n1,sender,vrecv,recever,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: n1
   real(dp),intent(in) :: vsend(:)
   real(dp),intent(inout) :: vrecv(:)
   integer,intent(in) :: sender,recever,spaceComm
   integer,intent(out)   :: ier
  end subroutine xexch_mpi_dpn

  subroutine xexch_mpi_dp2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: nt
   real(dp),intent(in) :: vsend(:,:)
   real(dp),intent(inout) :: vrecv(:,:)
   integer,intent(in) :: sender,recever,spaceComm
   integer,intent(out)   :: ier
  end subroutine xexch_mpi_dp2d

  subroutine xexch_mpi_dp3d(vsend,nt,sender,vrecv,recever,spaceComm,ier)
   use defs_basis
   integer,intent(in) :: nt
   real(dp),intent(in) :: vsend(:,:,:)
   real(dp),intent(inout) :: vrecv(:,:,:)
   integer,intent(in) :: sender,recever,spaceComm
   integer,intent(out)   :: ier
  end subroutine xexch_mpi_dp3d

 end interface

!End of the generic interface of xexch_mpi

!Generic interface of the routines xcast_mpi

 interface xcast_mpi

!   subroutine xcast_mpi_intoff(xval,master,spaceComm,ier)
!     use defs_basis
!     integer(abinit_offset), intent(inout) :: xval
!     integer,intent(in) :: spaceComm,master
!     integer,intent(out) :: ier
!   end subroutine xcast_mpi_intoff

  subroutine xcast_mpi_intv(xval,master,spaceComm,ier)
   use defs_basis
   integer,intent(inout) :: xval
   integer,intent(in) :: spaceComm,master
   integer,intent(out)   :: ier
  end subroutine xcast_mpi_intv

  subroutine xcast_mpi_int1d(xval,master,spaceComm,ier)
   use defs_basis
   integer,intent(inout) :: xval(:)
   integer,intent(in) :: spaceComm,master
   integer,intent(out)   :: ier
  end subroutine xcast_mpi_int1d

  subroutine xcast_mpi_int2d(xval,master,spaceComm,ier)
   use defs_basis
   integer,intent(inout) :: xval(:,:)
   integer,intent(in) :: spaceComm,master
   integer,intent(out)   :: ier
  end subroutine xcast_mpi_int2d

  subroutine xcast_mpi_int3d(xval,master,spaceComm,ier)
   use defs_basis
   integer,intent(inout) :: xval(:,:,:)
   integer,intent(in) :: spaceComm,master
   integer,intent(out)   :: ier
  end subroutine xcast_mpi_int3d

  subroutine xcast_mpi_dpv(xval,master,spaceComm,ier)
   use defs_basis
   real(dp),intent(inout) :: xval
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_dpv

  subroutine xcast_mpi_dp1d(xval,master,spaceComm,ier)
   use defs_basis
   real(dp),intent(inout) :: xval(:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_dp1d

  subroutine xcast_mpi_dp2d(xval,master,spaceComm,ier)
   use defs_basis
   real(dp),intent(inout) :: xval(:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_dp2d

  subroutine xcast_mpi_dp3d(xval,master,spaceComm,ier)
   use defs_basis
   real(dp),intent(inout) :: xval(:,:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_dp3d

  subroutine xcast_mpi_dp4d(xval,master,spaceComm,ier)
   use defs_basis
   real(dp),intent(inout) :: xval(:,:,:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
   integer::n1,n2,n3,n4
  end subroutine xcast_mpi_dp4d

  subroutine xcast_mpi_spv(xval,master,spaceComm,ier)
   use defs_basis
   real,intent(inout) :: xval
   integer,intent(in) :: spaceComm,master
   integer,intent(out) :: ier
  end subroutine xcast_mpi_spv

  subroutine xcast_mpi_sp1d(xval,master,spaceComm,ier)
   use defs_basis
   real,intent(inout) :: xval(:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_sp1d

  subroutine xcast_mpi_sp2d(xval,master,spaceComm,ier)
   use defs_basis
   real,intent(inout) :: xval(:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_sp2d

  subroutine xcast_mpi_sp3d(xval,master,spaceComm,ier)
   use defs_basis
   real,intent(inout) :: xval(:,:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_sp3d

  subroutine xcast_mpi_sp4d(xval,master,spaceComm,ier)
   use defs_basis
   real,intent(inout) :: xval(:,:,:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_sp4d

  subroutine xcast_mpi_cplxv(xval,master,spaceComm,ier)
   use defs_basis
   complex,intent(inout) :: xval
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_cplxv

  subroutine xcast_mpi_cplx1d(xval,master,spaceComm,ier)
   use defs_basis
   complex,intent(inout) :: xval(:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_cplx1d

  subroutine xcast_mpi_cplx2d(xval,master,spaceComm,ier)
   use defs_basis
   complex,intent(inout) :: xval(:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_cplx2d

  subroutine xcast_mpi_cplx3d(xval,master,spaceComm,ier)
   use defs_basis
   complex,intent(inout) :: xval(:,:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_cplx3d

  subroutine xcast_mpi_dcv(xval,master,spaceComm,ier)
   use defs_basis
   complex(dpc),intent(inout):: xval
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_dcv

  subroutine xcast_mpi_dc1d(xval,master,spaceComm,ier)
   use defs_basis
   complex(dpc),intent(inout):: xval(:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_dc1d

  subroutine xcast_mpi_dc2d(xval,master,spaceComm,ier)
   use defs_basis
   complex(dpc),intent(inout):: xval(:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_dc2d

  subroutine xcast_mpi_dc3d(xval,master,spaceComm,ier)
   use defs_basis
   complex(dpc),intent(inout):: xval(:,:,:)
   integer ,intent(in) :: spaceComm,master
   integer ,intent(out)   :: ier
  end subroutine xcast_mpi_dc3d

 end interface
!End of the generic interface of xcast_mpi

 end module defs_xfuncmpi
!!***
