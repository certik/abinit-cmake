!!****f* ABINIT/pre_scatter
!!
!! NAME
!! pre_scatter
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .

!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pre_scatter(array,array_allgather,n1,n2,n3,mpi_enreg,option)
 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

 type(mpi_type) :: mpi_enreg
 integer :: spaceComm
 integer ::  old_paral_level,ier,n1,n2,n3
 real(dp) :: array(n1,n2,n3/mpi_enreg%nproc_fft,1)
 real(dp) :: array_allgather(n1,n2,n3,1)
 character(*) :: option



 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm)
 if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
!Gather the array on all procs
 if(option=='gather') then
 call xallgather_mpi(array,n1*n2*n3/mpi_enreg%nproc_fft,array_allgather,spaceComm,ier)
! Perform the reverse operation
 elseif(option=='scatter') then
 array(:,:,:,:)=array_allgather(:,:,&
&               n3/mpi_enreg%nproc_fft*mpi_enreg%me_fft+1:n3/mpi_enreg%nproc_fft*(mpi_enreg%me_fft+1),:)
 endif
 mpi_enreg%paral_level=old_paral_level
end subroutine pre_scatter
!!***
