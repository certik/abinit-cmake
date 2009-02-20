!{\src2tex{textfont=tt}}
!!****f* ABINIT/xalltoall_mpi
!! NAME
!! xalltoall_mpi
!!
!! FUNCTION
!! This module contains functions that calls MPI routine,
!! if we compile the code using the MPI CPP flags.
!! xalltoall_mpi is the generic function.
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

subroutine xalltoall_mpi_dp2d(xval, sendsize, recvbuf, recvsize, spaceComm, ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(in)    :: xval(:,:)
 real(dp),intent(inout) :: recvbuf(:,:)
 integer ,intent(in)    :: sendsize, recvsize
 integer ,intent(in)    :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
 if (spaceComm /=  MPI_COMM_SELF) then
    !allgather xval on all proc. in spaceComm
    call MPI_ALLTOALL(xval, sendsize, MPI_DOUBLE_PRECISION, recvbuf, &
         & recvsize, MPI_DOUBLE_PRECISION, spaceComm, ier)
 end if
#endif
end subroutine xalltoall_mpi_dp2d
!!***
