!{\src2tex{textfont=tt}}
!!****f* ABINIT/xallgather_mpi
!! NAME
!! xallgather_mpi
!!
!! FUNCTION
!! This module contains functions that calls MPI routine,
!! if we compile the code using the  MPI CPP flags.
!! xallgather_mpi is the generic function.
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

subroutine xallgather_mpi_int(xval,recvcounts,spaceComm,ier)

 use defs_basis
#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(inout) :: xval
 integer,intent(inout) :: recvcounts(:)
 integer,intent(in) :: spaceComm
 integer,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
                if (spaceComm /=  MPI_COMM_SELF) then
          !allgather xval on all proc. in spaceComm
                call MPI_ALLGATHER(xval,1,MPI_INTEGER,recvcounts,1,MPI_INTEGER,&
          &  spaceComm,ier)
                end if
#endif
end subroutine xallgather_mpi_int

!--------------------------------------------------------------------

subroutine xallgather_mpi_dp1d(xval,nelem,recvcounts,spaceComm,ier)

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
   real(dp),intent(inout) :: recvcounts(:)
   integer ,intent(in) :: nelem,spaceComm
   integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
                if (spaceComm /=  MPI_COMM_SELF) then
          !allgather xval on all proc. in spaceComm
                call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvcounts,nelem,&
          &  MPI_DOUBLE_PRECISION,spaceComm,ier)
                end if
#endif
end subroutine xallgather_mpi_dp1d

!--------------------------------------------------------------------

subroutine xallgather_mpi_dp4d(xval,nelem,recvcounts,spaceComm,ier)

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
   real(dp),intent(inout) :: recvcounts(:,:,:,:)
   integer ,intent(in) :: nelem,spaceComm
   integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
                if (spaceComm /=  MPI_COMM_SELF) then
          !allgather xval on all proc. in spaceComm
                call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvcounts,nelem,&
          &  MPI_DOUBLE_PRECISION,spaceComm,ier)
                end if
#endif
end subroutine xallgather_mpi_dp4d
!!***
