!{\src2tex{textfont=tt}}
!!****f* ABINIT/xdef_comm
!! NAME
!! xdef_comm
!!
!! FUNCTION
!! Defines communicator and tools for MPI.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! Should become a module.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! xcomm_init  definition

!Used to distinguish serial execution when performing reduction operations.

subroutine xcomm_world(mpi_enreg,spaceComm)

!BEGIN TF_CHANGES
use defs_datatypes
!END TF_CHANGES

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
!BEGIN TF_CHANGES
 type(MPI_type) :: mpi_enreg
!END TF_CHANGES
 integer,intent(out) :: spaceComm

!Local variables-------------------

#if defined MPI
!BEGIN TF_CHANGES
                if (mpi_enreg%paral_compil_respfn == 1) then
                  spaceComm = mpi_enreg%spaceComm
                else
                    spaceComm = MPI_COMM_WORLD
                end if
!END TF_CHANGES
#else
                spaceComm = abinit_comm_serial
#endif
end subroutine xcomm_world
!!***


subroutine xcomm_init(mpi_enreg,spaceComm)
 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: spaceComm

!Local variables-------------------

#if defined MPI
!BEGIN TF_CHANGES
!        write(6,*) 'xcomm_init: enter'
!        write(6,*) 'paral_level=',mpi_enreg%paral_level
! comment
! init paral_level  in gstate and driver.f90
! if parareel == 0 ==> paral_level = 2
! if parareel == 1 paral_level = 1
!
       if (mpi_enreg%paral_compil_respfn == 1) then
               if (mpi_enreg%paral_level == 2) then

                        spaceComm = mpi_enreg%spaceComm
               else
                        if (mpi_enreg%num_group_fft /= 0) then
                                spaceComm =  mpi_enreg%fft_comm(mpi_enreg%num_group_fft)
                                else
                                spaceComm = MPI_COMM_SELF
                        end if
               end if
!END TF_CHANGES
       elseif (mpi_enreg%paral_level > 1) then
                if (mpi_enreg%paral_level == 2) then
                   spaceComm = MPI_COMM_WORLD
                else
                   if (mpi_enreg%num_group_fft /= 0) then
                      spaceComm =  mpi_enreg%fft_comm(mpi_enreg%num_group_fft)
                   else
                      spaceComm = MPI_COMM_SELF
                   end if
                end if
       else if (associated(mpi_enreg%kpt_comm_para)) then
                spaceComm = mpi_enreg%kpt_comm_para(mpi_enreg%ipara)
       else
          spaceComm = MPI_COMM_WORLD
       end if
#else
                spaceComm = abinit_comm_serial
#endif
end subroutine xcomm_init


!Define master
subroutine xmaster_init(mpi_enreg,master)
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: master

!Local variables-------------------

#if defined MPI
            if (mpi_enreg%parareel == 0) then
                master = 0
                else
                master = mpi_enreg%master_group_para
                end if
#else
                master = 0
#endif
end subroutine xmaster_init

!Define master_fft
subroutine xmaster_init_fft(mpi_enreg,master)
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: master

!Local variables-------------------

#if defined MPI
            if (mpi_enreg%paral_fft == 0) then
                master = mpi_enreg%me
                else
                master = mpi_enreg%master_fft
                end if
#else
                master = 0
#endif
end subroutine xmaster_init_fft

!Define me
subroutine xme_init(mpi_enreg,me)
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 type(MPI_type),intent(in) :: mpi_enreg
 integer,intent(out) :: me

!Local variables-------------------

!BEGIN TF_CHANGES
#if defined MPI
            if(mpi_enreg%paral_compil_respfn == 1 .AND. mpi_enreg%me_respfn/=-1) then
                me = mpi_enreg%me_respfn
            else
                me=mpi_enreg%me
            end if
#else
        me = 0
#endif
!END TF_CHANGES
end subroutine xme_init

!Define ntot proc
subroutine xproc_init(mpi_enreg,nproc_max)
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(out) :: nproc_max
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------

#if defined MPI
        if (mpi_enreg%parareel == 0) then
!BEGIN TF_CHANGES
            if (mpi_enreg%paral_compil_respfn == 1) then
              nproc_max=mpi_enreg%nproc_respfn
            else
              nproc_max=mpi_enreg%nproc
            end if
!END TF_CHANGES
        else
                 nproc_max=mpi_enreg%nproc_group_para
        end if

#else
        nproc_max = 1
#endif
end subroutine xproc_init

!Define who i am
subroutine xme_whoiam(me)


#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(out) :: me

!Local variables-------------------

#if defined MPI
           integer :: ierr
#endif
 me=0
#if defined MPI
        call MPI_COMM_RANK(MPI_COMM_WORLD,me,ierr)
#endif
end subroutine xme_whoiam

!Define ntotmax proc
subroutine xproc_max(nproc,ierr)


#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments-------------------------
 integer,intent(out) :: nproc,ierr

!Local variables-------------------

 ierr=0
 nproc = 1
#if defined MPI
        call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
#endif
end subroutine xproc_max
