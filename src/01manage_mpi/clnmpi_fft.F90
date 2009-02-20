!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_fft
!! NAME
!! clnmpi_fft
!!
!! FUNCTION
!! Clean up the mpi informations for FFT or BAND-FFT parallelism
!! (mostly deallocate parts of mpi_enreg)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (AR, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!   If the cpp option MPI is activated, deallocate
!!   mpi_enreg%proc_distrb
!!   mpi_enreg%band_comm
!!   mpi_enreg%kpt_group
!!
!! TODO
!!
!! PARENTS
!!      gstate,invars2m
!!
!! CHILDREN
!!      mpi_comm_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine clnmpi_fft(dtset, mpi_enreg)

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

 type(dataset_type), intent(in) :: dtset
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
#if defined MPI
         integer :: iisppol,nsppol,iikpt,nkpt,group,ierr,result
#endif

! ***********************************************************************

!DEBUG
! write(6,*)' clnmpi_fft : enter'
!ENDDEBUG

#if defined MPI
          nkpt=dtset%nkpt
          nsppol=dtset%nsppol
#endif

#if defined MPI
          if ( mpi_enreg%fft_master_comm /= MPI_COMM_NULL .and. &
&              mpi_enreg%fft_master_comm /= MPI_COMM_SELF      ) then
#if !defined FC_MIPSPRO
           call MPI_COMM_FREE(mpi_enreg%fft_master_comm,ierr)
#endif
          end if
#endif

#if defined MPI
         if (mpi_enreg%paral_fft==1) then
          do iisppol=1,nsppol
           do iikpt=1,nkpt
            group=iikpt+(iisppol-1)*nkpt
            if ( mpi_enreg%fft_comm(group) /= MPI_COMM_NULL .and. &
&                mpi_enreg%fft_comm(group) /= MPI_COMM_SELF      ) then
#if !defined FC_MIPSPRO
             call MPI_COMM_FREE(mpi_enreg%fft_comm(group),ierr)
#endif
            end if
           end do
          end do
         end if
#endif

#if defined MPI
          deallocate(mpi_enreg%fft_comm)
          deallocate(mpi_enreg%fft_group)
#endif

 if (mpi_enreg%paral_compil_fft==1) then
  deallocate(mpi_enreg%nplanes_fft)
  deallocate(mpi_enreg%ind_fft_planes)
 end if

!DEBUG
! write(6,*)' clnmpi_fft : end'
!ENDDEBUG

end subroutine clnmpi_fft
!!***
