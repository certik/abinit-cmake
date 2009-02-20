!{\src2tex{textfont=tt}}
!!****f* ABINIT/xderiveWrite
!! NAME
!! xderiveWrite
!!
!! FUNCTION
!! Generic routine to read/write wf files.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!   we use several procedures with the same generic name
!!   xderiveWrite  contains
!!               xderiveWrite_int  :  write integer  value
!!               xderiveWrite_int2d  :  write integer  array 2d
!!               xderiveWrite_dp   :  write double precision value
!!               xderiveWrite_dp2d   : write double precision array 2d
!!
!! PARENTS
!!      outxfhist,rwwf
!!
!! CHILDREN
!!      MPI_FILE_WRITE_AT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xderiveWrite_int(wff,xval,n1,ierr)

 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: xval(:),n1
 integer,intent(out) :: ierr

 ierr = 0
#if defined MPI_IO
            call MPI_FILE_WRITE_AT(wff%fhwff,  wff%offwff,xval,n1  &
            & , MPI_INTEGER , statux, ierr)

            wff%offwff = wff%offwff + wff%nbOct_int * n1
#endif

end subroutine xderiveWrite_int

subroutine xderiveWrite_int_mpio(wff,xval,n1,ierr,spaceComm)

 use defs_basis
 use defs_datatypes
 
#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: xval(:)

!Local variables
 !integer :: iproc,last_size,me,n3,numproc
 !integer, allocatable:: local_offset(:)
 integer(abinit_offset) :: nboct,dispoct,totoct
 integer(abinit_offset) :: posit
 ierr=0
#if defined MPI_IO
 nboct = n1*wff%nbOct_int
 posit = wff%offwff

 ! dispoct = sum (nboct, rank=0..me)
 call MPI_SCAN(nboct,dispoct,1,MPI_INTEGER8,MPI_SUM,spaceComm,ierr)
 posit = posit+dispoct-nboct
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1,MPI_INTEGER,spaceComm,ierr)
 posit = posit+nboct
 ! gather the bigest offset
 call MPI_ALLREDUCE(dispoct,totoct,1,MPI_INTEGER8,MPI_MAX,spaceComm,ierr)
 wff%offwff = wff%offwff+totoct

 ! Disable old code
#ifdef DEADCODE
 call MPI_COMM_SIZE(spaceComm,numproc,ierr)
 call MPI_COMM_RANK(spaceComm,me,ierr)
 allocate(local_offset(numproc))
 n3=n1
 !call xallgather_mpi(n3,local_offset,spaceComm,ierr)
 call flush(6)
 call MPI_ALLGATHER(n3,1,MPI_INTEGER,local_offset,1,MPI_INTEGER,&
      &  spaceComm,ierr)
 last_size=local_offset(1)
 !local_offset(1)=0
 do iproc=2,numproc
    local_offset(iproc)=local_offset(iproc-1)+local_offset(iproc)
 enddo
 local_offset(:)=local_offset(:)-last_size
 call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff+local_offset(me+1),xval,n1 &
      &  , MPI_INTEGER , statux, ierr)
 wff%offwff = wff%offwff + wff%nbOct_int * (local_offset(numproc)+last_size)
 deallocate(local_offset)
#endif
#endif
end subroutine xderiveWrite_int_mpio

subroutine xderiveWrite_int2d(wff,xval,n1,n2,ierr)

 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif

 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: xval(:,:),n1,n2
 integer,intent(out) :: ierr

 ierr = 0
#if defined MPI_IO
           call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff,xval,n1*n2  &
           & , MPI_INTEGER , statux, ierr)

           wff%offwff = wff%offwff + wff%nbOct_int * n1 *n2
#endif
end subroutine xderiveWrite_int2d

subroutine xderiveWrite_int2d_mpio(wff,xval,n1,n2,ierr,spaceComm)

 use defs_basis
 use defs_datatypes 
 
#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 integer,intent(in):: xval(:,:)

!Local variables
 !integer :: iproc,last_size,me,n3,numproc,dispoct,posit
 !integer, allocatable:: local_offset(:)
 integer(abinit_offset) :: nboct,dispoct,totoct
 integer(abinit_offset) :: posit
 ierr=0
#if defined MPI_IO
 nboct = n1*n2*wff%nbOct_int
 posit = wff%offwff

 ! dispoct = sum(nboct, rank=0..me)
 call MPI_SCAN(nboct,dispoct,1,MPI_INTEGER8,MPI_SUM,spaceComm,ierr)
 posit = posit + dispoct-nboct
 call mpi_file_write_at(wff%fhwff,posit,xval,n1*n2,MPI_INTEGER,statux,ierr)
 posit = posit + nboct
 ! gather the biggest offset
 call MPI_ALLREDUCE(dispoct,totoct,1,MPI_INTEGER8,MPI_MAX,spaceComm,ierr)
 wff%offwff = wff%offwff + totoct

 ! Disable old code
#ifdef DEADCODE
 call MPI_COMM_SIZE(spaceComm,numproc,ierr)
 call MPI_COMM_RANK(spaceComm,me,ierr)
 allocate(local_offset(numproc))
 n3=n2*n1
 !call xallgather_mpi(n3,local_offset,spaceComm,ierr)
 call flush(6)
 call MPI_ALLGATHER(n3,1,MPI_INTEGER,local_offset,1,MPI_INTEGER,&
      &  spaceComm,ierr)
 last_size=local_offset(1)
 !local_offset(1)=0
 do iproc=2,numproc
    local_offset(iproc)=local_offset(iproc-1)+local_offset(iproc)
 enddo
 local_offset(:)=local_offset(:)-last_size
 call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff+local_offset(me+1),xval,n1*n2 &
      &  , MPI_INTEGER , statux, ierr)
 wff%offwff = wff%offwff + wff%nbOct_int * (local_offset(numproc)+last_size)
 deallocate(local_offset)
#endif
#endif
end subroutine xderiveWrite_int2d_mpio

subroutine xderiveWrite_int2d_mpio_arr(wff,xval,n1,n2,ierr,local_offset)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,local_offset(:)
 integer,intent(out) :: ierr
 integer,intent(in):: xval(:,:)

!Local variables
 integer :: iproc,last_size,me,n3,numproc
 integer, allocatable :: loc_arr(:)
 ierr=0
#if defined MPI_IO
           numproc=size(local_offset)
           allocate(loc_arr(numproc))
	   loc_arr(:)=local_offset(:)
           n3=n2*n1
	   last_size=loc_arr(1)
	   do iproc=2,numproc
	    loc_arr(iproc)=loc_arr(iproc-1)+loc_arr(iproc)
	   enddo
	   loc_arr(:)=loc_arr(:)-last_size
           call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff+loc_arr(me+1),xval,n1*n2 &
           &  , MPI_INTEGER , statux, ierr)
           wff%offwff = wff%offwff + wff%nbOct_int * (loc_arr(numproc)+last_size)
	   deallocate(loc_arr)
#endif
end subroutine xderiveWrite_int2d_mpio_arr

subroutine xderiveWrite_dp(wff,xval,n1,ierr)
 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
 integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1
 integer,intent(out) :: ierr
 real(dp),intent(in) :: xval(:)

 ierr=0
#if defined MPI_IO
           call MPI_FILE_WRITE_AT(wff%fhwff,wff%offwff,xval,n1 &
           & , MPI_DOUBLE_PRECISION , statux, ierr)
           wff%offwff = wff%offwff + wff%nbOct_dp * n1

#endif
end subroutine xderiveWrite_dp

subroutine xderiveWrite_dp_mpio(wff,xval,n1,ierr,spaceComm)

 use defs_basis
 use defs_datatypes 

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:)

!Local variables
 integer :: iproc,last_size,me,n3,numproc
 integer, allocatable:: local_offset(:)
 integer(abinit_offset) :: nboct,dispoct,totoct,posit
!call leave_new("COLL")
 ierr=0
#if defined MPI_IO
 nboct = n1*wff%nbOct_dp
 posit = wff%offwff
 ! dispoct = sum (nboct, rank = 0..me)
 call MPI_SCAN(nboct,dispoct,1,MPI_INTEGER8,MPI_SUM,spaceComm,ierr)
 posit = posit + dispoct - nboct
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1 &
      & , MPI_DOUBLE_PRECISION , statux, ierr)
 posit = posit + nboct
 ! Gather the biggest offset
 call MPI_ALLREDUCE(dispoct,totoct,1,MPI_INTEGER8,MPI_MAX,spaceComm,ierr)
 wff%offwff = wff%offwff + totoct

 ! Disable old code
#ifdef DEADCODE
 call MPI_COMM_SIZE(spaceComm,numproc,ierr)
 call MPI_COMM_RANK(spaceComm,me,ierr)
 allocate(local_offset(numproc))
 n3=n1
 !call xallgather_mpi(n3,local_offset,spaceComm,ierr)
 call flush(6)
 call MPI_ALLGATHER(n3,1,MPI_INTEGER,local_offset,1,MPI_INTEGER,&
      &  spaceComm,ierr)
 last_size=local_offset(1)
 !local_offset(1)=0
 do iproc=2,numproc
    local_offset(iproc)=local_offset(iproc-1)+local_offset(iproc)
 enddo
 local_offset(:)=local_offset(:)-last_size
 call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff+local_offset(me+1),xval,n3 &
      &  , MPI_DOUBLE_PRECISION , statux, ierr)
 wff%offwff = wff%offwff + wff%nbOct_dp * (local_offset(numproc)+last_size)
 deallocate(local_offset)
#endif
#endif
end subroutine xderiveWrite_dp_mpio

subroutine xderiveWrite_dp2d_mpio_arr(wff,xval,n1,n2,ierr,local_offset)

 use defs_basis
 use defs_datatypes 
 
#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,local_offset(:)
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:,:)

!Local variables
 integer :: iproc,last_size,me,n3,numproc
 integer, allocatable :: loc_arr(:)
 ierr=0
#if defined MPI_IO
           numproc=size(local_offset)
           allocate(loc_arr(numproc))
	   loc_arr(:)=local_offset(:)
           n3=n2*n1
	   last_size=loc_arr(1)
	   do iproc=2,numproc
	    loc_arr(iproc)=loc_arr(iproc-1)+loc_arr(iproc)
	   enddo
	   loc_arr(:)=loc_arr(:)-last_size
           call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff+loc_arr(me+1),xval,n1*n2 &
           &  , MPI_DOUBLE_PRECISION , statux, ierr)
           wff%offwff = wff%offwff + wff%nbOct_dp * (loc_arr(numproc)+last_size)
	   deallocate(loc_arr)
#endif
end subroutine xderiveWrite_dp2d_mpio_arr

subroutine xderiveWrite_dp2d(wff,xval,n1,n2,ierr)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:,:)

 ierr=0
#if defined MPI_IO
           call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff,xval,n1*n2 &
           &  , MPI_DOUBLE_PRECISION , statux, ierr)
           wff%offwff = wff%offwff + wff%nbOct_dp * n1*n2
#endif
end subroutine xderiveWrite_dp2d


subroutine xderiveWrite_dp2d_mpio(wff,xval,n1,n2,ierr,spaceComm)
 use defs_basis
 use defs_datatypes 

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(in):: xval(:,:)

!Local variables
 integer :: iproc,last_size,me,n3,numproc
 integer, allocatable:: local_offset(:)
 integer(abinit_offset) :: nboct,dispoct,totoct,posit
 ierr=0
#if defined MPI_IO
 nboct = n1*n2*wff%nbOct_dp
 posit = wff%offwff
 ! dispoct = sum(nboct, rank=0..me)
 call MPI_SCAN(nboct,dispoct,1,MPI_INTEGER8,MPI_SUM,spaceComm,ierr)
 posit = posit+dispoct-nboct
 call MPI_FILE_WRITE_AT(wff%fhwff,posit,xval,n1*n2,MPI_DOUBLE_PRECISION,statux,ierr)
 posit = posit+nboct
 ! gather the biggest offset
 call MPI_ALLREDUCE(dispoct,totoct,1,MPI_INTEGER8,MPI_MAX,spaceComm,ierr)
 wff%offwff = wff%offwff+totoct
 ! old code
#ifdef DEADCODE
 call MPI_COMM_SIZE(spaceComm,numproc,ierr)
 call MPI_COMM_RANK(spaceComm,me,ierr)
 allocate(local_offset(numproc))
 n3=n2*n1
 !call xallgather_mpi(n3,local_offset,spaceComm,ierr)
 call flush(6)
 call MPI_ALLGATHER(n3,1,MPI_INTEGER,local_offset,1,MPI_INTEGER,&
      &  spaceComm,ierr)
 last_size=local_offset(1)
 !local_offset(1)=0
 do iproc=2,numproc
    local_offset(iproc)=local_offset(iproc-1)+local_offset(iproc)
 enddo
 local_offset(:)=local_offset(:)-last_size
 call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff+local_offset(me+1),xval,n1*n2 &
      &  , MPI_DOUBLE_PRECISION , statux, ierr)
 wff%offwff = wff%offwff + wff%nbOct_dp * (local_offset(numproc)+last_size)
 deallocate(local_offset)
#endif
#endif
end subroutine xderiveWrite_dp2d_mpio
!!***
