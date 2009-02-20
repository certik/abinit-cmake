!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_band
!! NAME
!! clnmpi_band
!!
!! FUNCTION
!! Clean up the mpi informations for the BAND parallelism (paralbd=1)
!!
!! COPYRIGHT
!! Copyright (C) 2008-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  nkpt= number of k-points
!!  nsppol= 1 for unpolarized, 2 for polarized
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!    mpi_enreg%band_comm(nkpt*nsppol)=comm array of BAND set
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      mpi_comm_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine clnmpi_band(nkpt,nsppol,mpi_enreg)

 use defs_basis
 use defs_datatypes
 use defs_infos

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: nkpt,nsppol
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
#if defined MPI
         integer :: ierr,ii,isppol,ikpt
#endif

! ***********************************************************************

!DEBUG
! write(6,*)' clnmpi_band : enter'
!ENDDEBUG

#if defined MPI
 if (mpi_enreg%has_band_comm==1) then

  do isppol=1,nsppol
   do ikpt=1,nkpt
    ii=ikpt+(isppol-1)*nkpt
    if (mpi_enreg%band_comm(ii)/=MPI_COMM_SELF) then
     call MPI_COMM_FREE(mpi_enreg%band_comm(ii),ierr)
    end if
   end do
  end do

  deallocate(mpi_enreg%band_comm)

 end if
#endif

!DEBUG
! write(6,*)' clnmpi_band : end'
!ENDDEBUG

end subroutine clnmpi_band
!!***
