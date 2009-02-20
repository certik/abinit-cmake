!{\src2tex{textfont=tt}}
!!****f* ABINIT/xalltoallv_mpi
!! NAME
!! xalltoallv_mpi
!!
!! FUNCTION
!! This module contains functions that calls MPI routine,
!! if we compile the code using the   CPP flags.
!! xalltoallv_mpi is the generic function.
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

subroutine xalltoallv_mpi_dp2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)

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
 integer ,intent(in) :: sendcnts(:),sdispls(:),rdispls(:),recvcnts(:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
 if (spaceComm /=  MPI_COMM_SELF) then
!allgather xval on all proc. in spaceComm
  call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_DOUBLE_PRECISION,recvbuf,&
&                    recvcnts,rdispls,MPI_DOUBLE_PRECISION,spaceComm,ier)
end if
#endif
end subroutine xalltoallv_mpi_dp2d


subroutine xalltoallv_mpi_int2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,spaceComm,ier)

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
 integer ,intent(in) :: sendcnts(:),sdispls(:),rdispls(:),recvcnts(:)
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI 
 if (spaceComm /=  MPI_COMM_SELF) then
!allgather xval on all proc. in spaceComm
  call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_INTEGER,recvbuf,&
&                    recvcnts,rdispls,MPI_INTEGER,spaceComm,ier)
end if
#endif
end subroutine xalltoallv_mpi_int2d
!!***
