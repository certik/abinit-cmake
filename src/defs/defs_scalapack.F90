!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_scalapack
!! NAME
!! defs_scalapack
!!
!! FUNCTION
!! This module contains descriptions of the variables used in the ScaLAPACK routines. 
!! The ScaLAPACK routines and functions are defined in 01managempi/scalapack.F90
!! and used in 08seqpar/subdiago.F90.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (CS,GZ,FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!! To be translated
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_scalapack

  use defs_basis

  implicit none


!Parameters of the array of description of ScaLAPACK matrices
  integer, parameter :: DLEN_ = 9    ! length
  integer, parameter :: Dtype_ = 1   ! type
  integer, parameter :: CTXT_ = 2    ! BLACS context
  integer, parameter :: M_ = 3       ! nb global lines
  integer, parameter :: N_ = 4       ! nb global columns
  integer, parameter :: MB_ = 5      ! nb lines of a block
  integer, parameter :: NB_ = 6      ! nb columns of a bloc
  integer, parameter :: RSRC_ = 7    ! line of processors at the beginning 
  integer, parameter :: CSRC_ = 8    ! column of processors at the beginning 
  integer, parameter :: LLD_ = 9     ! local number of lines

! Grid of ScaLAPACK processors
  type grid_scalapack

     integer   :: nbprocs ! total number of processors
     integer,dimension(1:2)   :: dims

     integer   :: ictxt   ! blacs context

  end type grid_scalapack

! One processor in the grid
  type processor_scalapack

     integer    :: myproc ! number of the processor
     integer    :: comm   ! MPI communicator MPI underlying the grid BLACS
     integer,dimension(1:2)   :: coords

     type(grid_scalapack) :: grid   ! the grid to which the processor is associated 

  end type processor_scalapack

! Description of a ScaLAPACK matrix 
  type descript_scalapack
!     SEQUENCE
     integer, dimension(DLEN_) :: tab

  end type descript_scalapack

! The local part of a ScaLAPACK matrix 
  type matrix_scalapack

     complex(dpc),dimension(:,:),pointer :: buffer     ! buffer of the local part of the matrix 

     integer, dimension(:), pointer     :: ipiv

     integer,dimension(1:2)             :: sizeb_global

     type(processor_scalapack),pointer :: processor

     integer,dimension(1:2)             :: sizeb_blocs ! size of the block of consecutive data

     integer,dimension(1:2)             :: sizeb_local ! dimensions of the buffer

     type(descript_scalapack)         :: descript

  end type matrix_scalapack

end module defs_scalapack
!!***
