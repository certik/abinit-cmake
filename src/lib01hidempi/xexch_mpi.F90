!{\src2tex{textfont=tt}}
!!****f* ABINIT/xexch_mpi
!! NAME
!! xexch_mpi  
!!
!! FUNCTION
!!   exchange data
!! This module contains functions that calls MPI routine,
!! if we compile the code using the MPI CPP flags.
!! xexch_mpi is the generic function.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (MB)
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

!--------------------------------------------------------------------

subroutine xexch_mpi_intn(vsend,n1,sender,vrecv,recever,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 integer,intent(in) :: n1
 integer,intent(in) :: vsend(:)
 integer,intent(inout) :: vrecv(:)
 integer,intent(in) :: sender,recever,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined MPI
           integer :: statux(MPI_STATUS_SIZE)
           integer :: tag,me
#endif

 ier=0

#if defined MPI
           if ( sender == recever ) return
           call MPI_COMM_RANK(spaceComm,me,ier)
           tag=n1
           if ( recever == me ) then
            call MPI_RECV(vrecv,n1,MPI_INTEGER,sender,tag  &
            &                ,spaceComm,statux,ier)
           end if
           if ( sender == me ) then
            call MPI_SEND(vsend,n1,MPI_INTEGER,recever,tag  &
            &                ,spaceComm,statux,ier)
           end if
                
#endif

 ier=0

end subroutine xexch_mpi_intn

!--------------------------------------------------------------------

subroutine xexch_mpi_int2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 integer,intent(in) :: nt
 integer,intent(in) :: vsend(:,:)
 integer,intent(inout) :: vrecv(:,:)
 integer,intent(in) :: sender,recever,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined MPI
           integer :: statux(MPI_STATUS_SIZE)
           integer :: tag,me
#endif

 ier=0

#if defined MPI
           if ( sender == recever ) return
           call MPI_COMM_RANK(spaceComm,me,ier)
           tag=nt
           if ( recever == me ) then
            call MPI_RECV(vrecv,nt,MPI_INTEGER,sender,tag  &
            &                ,spaceComm,statux,ier)
           end if
           if ( sender == me ) then
            call MPI_SEND(vsend,nt,MPI_INTEGER,recever,tag  &
            &                ,spaceComm,statux,ier)
           end if
                
#endif

 ier=0

end subroutine xexch_mpi_int2d


!--------------------------------------------------------------------
subroutine xexch_mpi_dpn(vsend,n1,sender,vrecv,recever,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 integer,intent(in) :: n1
 real(dp),intent(in) :: vsend(:)
 real(dp),intent(inout) :: vrecv(:)
 integer,intent(in) :: sender,recever,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined MPI
           integer :: statux(MPI_STATUS_SIZE)
           integer :: tag,me
#endif

 ier=0

#if defined MPI
           if ( sender == recever ) return
           call MPI_COMM_RANK(spaceComm,me,ier)
           tag=n1
           if ( recever == me ) then
            call MPI_RECV(vrecv,n1,MPI_DOUBLE_PRECISION,sender,tag  &
            &                ,spaceComm,statux,ier)
           end if
           if ( sender == me ) then
            call MPI_SEND(vsend,n1,MPI_DOUBLE_PRECISION,recever,tag  &
            &                ,spaceComm,statux,ier)
           end if
#endif

 ier=0

end subroutine xexch_mpi_dpn

!--------------------------------------------------------------------

subroutine xexch_mpi_dp2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 integer,intent(in) :: nt
 real(dp),intent(in) :: vsend(:,:)
 real(dp),intent(inout) :: vrecv(:,:)
 integer,intent(in) :: sender,recever,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined MPI
           integer :: statux(MPI_STATUS_SIZE)
           integer :: tag,me
#endif

 ier=0

#if defined MPI
           if ( sender == recever ) return

           call MPI_COMM_RANK(spaceComm,me,ier)
           tag=nt
           if ( recever == me ) then
            call MPI_RECV(vrecv,nt,MPI_DOUBLE_PRECISION,sender,tag  &
            &                ,spaceComm,statux,ier)
           end if
           if ( sender == me ) then
            call MPI_SEND(vsend,nt,MPI_DOUBLE_PRECISION,recever,tag  &
            &                ,spaceComm,statux,ier)
           end if
               
#endif

 ier=0

end subroutine xexch_mpi_dp2d
!--------------------------------------------------------------------

subroutine xexch_mpi_dp3d(vsend,nt,sender,vrecv,recever,spaceComm,ier)

 use defs_basis

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments----------------
 integer,intent(in) :: nt
 real(dp),intent(in) :: vsend(:,:,:)
 real(dp),intent(inout) :: vrecv(:,:,:)
 integer,intent(in) :: sender,recever,spaceComm
 integer,intent(out)   :: ier

!Local variables--------------
#if defined MPI
           integer :: statux(MPI_STATUS_SIZE)
           integer :: tag,me
#endif

 ier=0

#if defined MPI
           if ( sender == recever ) return
           call MPI_COMM_RANK(spaceComm,me,ier)
           tag=nt
           if ( recever == me ) then
            call MPI_RECV(vrecv,nt,MPI_DOUBLE_PRECISION,sender,tag  &
            &                ,spaceComm,statux,ier)
           end if
           if ( sender == me ) then
            call MPI_SEND(vsend,nt,MPI_DOUBLE_PRECISION,recever,tag  &
            &                ,spaceComm,statux,ier)
           end if
               
#endif

 ier=0

end subroutine xexch_mpi_dp3d


!--------------------------------------------------------------------

!!***
