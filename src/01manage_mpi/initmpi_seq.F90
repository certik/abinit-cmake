!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_seq
!! NAME
!! initmpi_seq
!!
!! FUNCTION
!! Initialize the mpi informations for a sequential use of other routines
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!  mpi_enreg=informations about MPI parallelization
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! PARENTS
!!      phfrq3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initmpi_seq(mpi_enreg)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 type(MPI_type),intent(out) :: mpi_enreg

!Local variables-------------------------------

! ***********************************************************************

!DEBUG
!write(6,*)' initmpi_seq : enter'
!stop
!ENDDEBUG

 mpi_enreg%gwpara=0
 mpi_enreg%paral_compil_kpt=0
 mpi_enreg%paral_compil_fft=0
 mpi_enreg%paral_compil_mpio=0
 mpi_enreg%paral_level=0
 mpi_enreg%paralbd=0
 mpi_enreg%me=0
 mpi_enreg%nproc=0
 mpi_enreg%me_group=0
 mpi_enreg%nproc_group=0
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=0
 mpi_enreg%paral_fft=0
 mpi_enreg%me_g0=0
 mpi_enreg%num_group_fft=0
 mpi_enreg%num_group=0
 mpi_enreg%nproc_per_kpt=0
 mpi_enreg%world_group=0
 mpi_enreg%parareel=0
 mpi_enreg%has_band_comm=0

!DEBUG
!write(6,*)' initmpi_seq : exit '
!stop
!ENDDEBUG

end subroutine initmpi_seq
!!***
