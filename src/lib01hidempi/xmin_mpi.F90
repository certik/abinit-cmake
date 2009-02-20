!{\src2tex{textfont=tt}}
!!****f* ABINIT/xmin_mpi
!! NAME
!! xmin_mpi
!!
!! FUNCTION
!! This module contains functions that calls MPI routine,
!! if we compile the code using the MPI  CPP flags.
!! xmin_mpi is the generic function.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (AR,XG,MB)
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

subroutine xmin_mpi_intv(xval,xmin,spaceComm,ier)

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer ,intent(inout) :: xval,xmin
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
            call MPI_ALLREDUCE(xval,xmin,1,MPI_INTEGER,&
            &  MPI_MIN,spaceComm,ier)
#endif
end subroutine xmin_mpi_intv

!--------------------------------------------------------------------

subroutine xmin_mpi_dpv(xval,xmin,spaceComm,ier)
 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 real(dp),intent(inout) :: xval,xmin
 integer ,intent(in) :: spaceComm
 integer ,intent(out)   :: ier

!Local variables-------------------

 ier=0
#if defined MPI
           call MPI_ALLREDUCE(xval,xmin,1,MPI_DOUBLE_PRECISION,&
           &  MPI_MIN,spaceComm,ier)
#endif
end subroutine xmin_mpi_dpv
!!***
