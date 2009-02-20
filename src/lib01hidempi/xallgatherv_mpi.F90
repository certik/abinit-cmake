!{\src2tex{textfont=tt}}
!!****f* ABINIT/xallgatherv_mpi
!! NAME
!! xallgatherv_mpi
!!
!! FUNCTION
!! This module contains functions that calls MPI routine,
!! if we compile the code using the MPI CPP flags.
!! xallgatherv_mpi is the generic function.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (AR,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xallgatherv_mpi_int2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
           include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(in) :: xval(:,:)
 integer,intent(inout) :: recvbuf(:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

 ier=0
#if defined MPI
 if (spaceComm /=  MPI_COMM_SELF) then
          !allgather xval on all proc. in spaceComm
  call MPI_ALLGATHERV(xval,nelem,MPI_INTEGER,recvbuf,recvcounts,displs,&
         &  MPI_INTEGER,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_int2d


!--------------------------------------------------------------------
subroutine xallgatherv_mpi_int(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(in) :: xval(:)
 integer,intent(inout)   :: recvbuf(:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
 if (spaceComm /=  MPI_COMM_SELF) then
!allgather xval on all proc. in spaceComm
 call MPI_ALLGATHERV(xval,nelem,MPI_INTEGER,recvbuf,recvcounts,displs,&
          &  MPI_INTEGER,spaceComm,ier)
end if
#endif
end subroutine xallgatherv_mpi_int

!--------------------------------------------------------------------

subroutine xallgatherv_mpi_dp(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:)
 real(dp),intent(inout)   :: recvbuf(:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

 ier=0
#if defined MPI
 if (spaceComm /=  MPI_COMM_SELF) then
!allgather xval on all proc. in spaceComm
  call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
          &  MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_dp

!--------------------------------------------------------------------

subroutine xallgatherv_mpi_dp2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:)
 real(dp),intent(inout) :: recvbuf(:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

 ier=0
#if defined MPI
 if (spaceComm /=  MPI_COMM_SELF) then
!allgather xval on all proc. in spaceComm
  call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
          &  MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_dp2d

!--------------------------------------------------------------------

subroutine xallgatherv_mpi_dp3d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:,:)
 real(dp),intent(inout)   :: recvbuf(:,:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

 ier=0
#if defined MPI
 if (spaceComm /=  MPI_COMM_SELF) then
!allgather xval on all proc. in spaceComm
  call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
          &  MPI_DOUBLE_PRECISION,spaceComm,ier)
end if
#endif
end subroutine xallgatherv_mpi_dp3d

!--------------------------------------------------------------------

subroutine xallgatherv_mpi_dp4d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in) :: xval(:,:,:,:)
 real(dp),intent(inout)   :: recvbuf(:,:,:,:)
 integer,intent(in) :: recvcounts(:),displs(:)
 integer,intent(in) :: nelem,spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
if (spaceComm /=  MPI_COMM_SELF) then
!allgather xval on all proc. in spaceComm
  call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
          &  MPI_DOUBLE_PRECISION,spaceComm,ier)
 end if
#endif
end subroutine xallgatherv_mpi_dp4d
!!***
