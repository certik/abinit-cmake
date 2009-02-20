!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnmpi_gs
!! NAME
!! clnmpi_gs
!!
!! FUNCTION
!! Clean up the mpi informations for the ground-state datasets
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
!!   mpi_enreg%kpt_comm
!!   mpi_enreg%kpt_group
!!
!! TODO
!!
!! PARENTS
!!      gstate,pstate
!!
!! CHILDREN
!!      mpi_comm_free,mpi_group_free
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine clnmpi_gs(dtset, mpi_enreg)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------

 type(dataset_type) :: dtset
 type(MPI_type) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
#if defined MPI
          !Variables introduced for MPI version
           integer :: ierr,nkpt,iikpt,nsppol,iisppol,group,result
#endif

! ***********************************************************************

!DEBUG
! write(6,*)' clnmpi_gs : enter'
!stop
!ENDDEBUG

#if defined MPI
           nsppol=dtset%nsppol
           nkpt=dtset%nkpt
           if (mpi_enreg%parareel == 0) then
              deallocate(mpi_enreg%proc_distrb)
           end if
           if (mpi_enreg%paralbd >= 1) then
            do iisppol=1,nsppol
             do iikpt=1,nkpt
              group=iikpt+(iisppol-1)*nkpt
              if (mpi_enreg%kpt_comm(group) /= MPI_COMM_NULL) then
               call MPI_COMM_FREE(mpi_enreg%kpt_comm(group),ierr)
              end if
             end do
            end do
            call MPI_GROUP_FREE(mpi_enreg%world_group,ierr)
            deallocate(mpi_enreg%kpt_comm,mpi_enreg%kpt_group)
           end if
           if (mpi_enreg%parareel == 1) then
            do iisppol=1,nsppol
             do iikpt=1,nkpt
              group=iikpt+(iisppol-1)*nkpt
              if (mpi_enreg%kpt_comm(group) /= MPI_COMM_NULL) then
               call MPI_COMM_FREE(mpi_enreg%kpt_comm(group),ierr)
              end if
             end do
            end do
            call MPI_GROUP_FREE(mpi_enreg%world_group,ierr)
            deallocate(mpi_enreg%proc_distrb_para)
            deallocate(mpi_enreg%kpt_comm_para,mpi_enreg%kpt_group_para)
           end if
#endif

end subroutine clnmpi_gs
!!***
