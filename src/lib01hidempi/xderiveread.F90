!{\src2tex{textfont=tt}}
!!****f* ABINIT/xderiveRead
!! NAME
!! xderiveRead
!!
!! FUNCTION
!! Generic subroutines to read wf files.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!   we use several procedures with the same generic name
!!   xderiveRead  contains
!!               xderiveRead_int  :  read integer  value
!!               xderiveRead_int2d  :  read integer  array 2d
!!               xderiveRead_dp   :  read double precision value
!!               xderiveRead_dp2d   :  read double precision array 2d
!!
!! PARENTS
!!      outxfhist,rwwf
!!
!! CHILDREN
!!      MPI_FILE_READ_AT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xderiveRead_int(wff,xval,n1,ierr)

 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
         integer(abinit_offset)  :: posit
         integer  :: statux(MPI_STATUS_SIZE)
         integer :: delim_record,nboct
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval(:)
 integer,intent(in) :: n1
 integer,intent(out) :: ierr

 xval(:)=0 ; ierr=0 ! Initialization, for the compiler
#if defined MPI_IO
           nboct = wff%nbOct_int * n1
           posit = wff%offwff
           delim_record = posit - wff%off_recs   &
           &            + wff%lght_recs - wff%nbOct_int

           if ( delim_record >= nboct ) then
            call MPI_FILE_READ_AT(wff%fhwff, posit,xval,n1   &
            & , MPI_INTEGER , statux, ierr)
            posit = posit + nboct
           else
            ierr = 1
            nboct =0
           end if

           ! new offset
           wff%offwff=wff%offwff + nboct
#endif
end subroutine xderiveRead_int


subroutine xderiveRead_int_mpio(wff,xval,n1,ierr,spaceComm)
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
 integer(abinit_offset)  :: posit,nboct,dispoct,totoct
 integer  :: statux(MPI_STATUS_SIZE)
 integer :: delim_record
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval(:)
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr

 xval(:)=0 ; ierr=0 ! Initialization, for the compiler
#if defined MPI_IO
 nboct = wff%nbOct_int * n1
 posit = wff%offwff
 delim_record = posit - wff%off_recs   &
      &            + wff%lght_recs - wff%nbOct_int
 
 if ( delim_record >= nboct ) then
    ! Compute offset for local part
    ! dispoct = sum (nboct, rank=0..me)
    call MPI_SCAN(nboct,dispoct,1,MPI_INTEGER8,MPI_SUM,spaceComm,ierr)
    posit = posit+dispoct-nboct
    call MPI_FILE_READ_AT(wff%fhwff, posit,xval,n1   &
         & , MPI_INTEGER , statux, ierr)
    posit = posit + nboct

    ! get the total number of bits wrote by processors
    call MPI_ALLREDUCE(dispoct,totoct,1,MPI_INTEGER8,MPI_MAX,spaceComm,ierr)
 else
    ierr = 1
    nboct =0
    totoct = 0
 end if
 
 ! new offset
 !wff%offwff=wff%offwff + nboct
 wff%offwff = wff%offwff + totoct
#endif
end subroutine xderiveRead_int_mpio

subroutine xderiveRead_int2d(wff,xval,n1,n2,ierr)

 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer(abinit_offset)  :: posit
           integer  :: statux(MPI_STATUS_SIZE)
           integer :: delim_record,nboct
#endif

 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval(:,:)
 integer,intent(in) :: n1,n2
 integer,intent(out) :: ierr

 xval(:,:)=0 ; ierr=0 ! Initialization, for the compiler
#if defined MPI_IO

           nboct = wff%nbOct_int * n1 *n2
           posit = wff%offwff
           delim_record = posit - wff%off_recs   &
           &            + wff%lght_recs - wff%nbOct_int

           if ( delim_record >= nboct ) then
            call MPI_FILE_READ_AT(wff%fhwff, posit,xval,n1*n2   &
            & , MPI_INTEGER , statux, ierr)
            posit = posit + nboct
           else
            ierr = 1
            nboct =0
           end if

           ! new offset
           wff%offwff=wff%offwff + nboct
#endif

end subroutine xderiveRead_int2d

subroutine xderiveRead_int2d_mpio(wff,xval,n1,n2,ierr,spaceComm)
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
 integer(abinit_offset)  :: posit
 integer  :: statux(MPI_STATUS_SIZE)
 integer(abinit_offset) :: delim_record,nboct,dispoct,totoct
#endif

 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval(:,:)
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr

 xval(:,:)=0 ; ierr=0 ! Initialization, for the compiler
#if defined MPI_IO
 nboct = wff%nbOct_int * n1 * n2
 posit = wff%offwff
 delim_record = posit - wff%off_recs   &
      &            + wff%lght_recs - wff%nbOct_int

 if ( delim_record >= nboct ) then
    ! Compute offset for local part
    ! dispoct = sum (nboct, rank=0..me)
    call MPI_SCAN(nboct,dispoct,1,MPI_INTEGER8,MPI_SUM,spaceComm,ierr)
    posit = posit + dispoct - nboct
    call MPI_FILE_READ_AT(wff%fhwff, posit,xval,n1*n2   &
         & , MPI_INTEGER , statux, ierr)
    posit = posit + nboct

    ! get the total number of bits wrote by processors
    call MPI_ALLREDUCE(dispoct,totoct,1,MPI_INTEGER8,MPI_MAX,spaceComm,ierr)
 else
    ierr = 1
    nboct =0
    totoct = 0
 end if

 ! new offset
 !wff%offwff=wff%offwff + nboct
 wff%offwff=wff%offwff + totoct
#endif
end subroutine xderiveRead_int2d_mpio

subroutine xderiveRead_dp(wff,xval,n1,ierr)

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
           integer(abinit_offset)  :: posit
           integer  :: statux(MPI_STATUS_SIZE)
           integer :: delim_record,nboct
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval(:)

 integer :: i

 xval(:)=zero ; ierr=0 ! Initialization, for the compiler

#if defined MPI_IO
           nboct = wff%nbOct_dp * n1
           posit = wff%offwff
           delim_record = posit - wff%off_recs   &
           &          + wff%lght_recs - wff%nbOct_int

           if ( delim_record >= nboct ) then
            call MPI_FILE_READ_AT(wff%fhwff, posit,xval,n1   &
            &     , MPI_DOUBLE_PRECISION , statux, ierr)
           else
            ierr = 1
            nboct =0
           end if

           ! new offset
           wff%offwff=wff%offwff + nboct
#endif
    end subroutine xderiveRead_dp

subroutine xderiveRead_dp_mpio(wff,xval,n1,ierr,spaceComm)
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
 integer(abinit_offset)  :: posit,nboct,dispoct,totoct
 integer  :: statux(MPI_STATUS_SIZE)
 integer :: delim_record
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval(:)
 
 integer :: i
 
 xval(:)=zero ; ierr=0 ! Initialization, for the compiler

#if defined MPI_IO
 nboct = wff%nbOct_dp * n1
 posit = wff%offwff
 delim_record = posit - wff%off_recs   &
      &          + wff%lght_recs - wff%nbOct_int
 
 if ( delim_record >= nboct ) then
    ! Compute offset for local part
    ! dispoct = sum (nboct, rank=0..me)
    call MPI_SCAN(nboct,dispoct,1,MPI_INTEGER8,MPI_SUM,spaceComm,ierr)
    posit = posit + dispoct - nboct
    
    call MPI_FILE_READ_AT(wff%fhwff, posit,xval,n1   &
         &     , MPI_DOUBLE_PRECISION , statux, ierr)

    posit = posit + nboct

    ! get the total number of bits wrote by processors
    call MPI_ALLREDUCE(dispoct,totoct,1,MPI_INTEGER8,MPI_MAX,spaceComm,ierr)
 else
    ierr = 1
    nboct =0
    totoct = 0
 end if
 
 ! new offset
 !wff%offwff=wff%offwff + nboct
 wff%offwff=wff%offwff + totoct
#endif
end subroutine xderiveRead_dp_mpio


subroutine xderiveRead_dp2d(wff,xval,n1,n2,ierr)

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
           integer(abinit_offset)  :: posit
           integer  :: statux(MPI_STATUS_SIZE)
           integer :: delim_record,nboct
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval(:,:)

 integer :: i

 xval(:,:)=zero ; ierr=0 ! Initialization, for the compiler
#if defined MPI_IO

           nboct = wff%nbOct_dp * n1 *n2
           posit = wff%offwff
           delim_record = posit - wff%off_recs   &
           &          + wff%lght_recs - wff%nbOct_int

           if ( delim_record >= nboct ) then
            call MPI_FILE_READ_AT(wff%fhwff, posit,xval,n1*n2   &
            &     , MPI_DOUBLE_PRECISION , statux, ierr)
           else
            ierr = 1
            nboct =0
           end if

           ! new offset
           wff%offwff=wff%offwff + nboct
#endif

end subroutine xderiveRead_dp2d

subroutine xderiveRead_dp2d_mpio(wff,xval,n1,n2,ierr,spaceComm)

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
 integer(abinit_offset)  :: posit
 integer  :: statux(MPI_STATUS_SIZE)
 integer(abinit_offset) :: delim_record,nboct,dispoct,totoct
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: n1,n2,spaceComm
 integer,intent(out) :: ierr
 real(dp),intent(out) :: xval(:,:)

 integer :: i

 xval(:,:)=zero ; ierr=0 ! Initialization, for the compiler
#if defined MPI_IO

 nboct = wff%nbOct_dp * n1 *n2
 posit = wff%offwff
 delim_record = posit - wff%off_recs   &
      &          + wff%lght_recs - wff%nbOct_int
 if ( delim_record >= nboct ) then
    ! Compute offset for local part
    ! dispoct = sum (nboct, rank=0..me)
    call MPI_SCAN(nboct,dispoct,1,MPI_INTEGER8,MPI_SUM,spaceComm,ierr)
    posit = posit + dispoct - nboct
    call MPI_FILE_READ_AT(wff%fhwff, posit,xval,n1*n2   &
         &     , MPI_DOUBLE_PRECISION , statux, ierr)

    posit = posit + nboct

    ! get the total number of bits wrote by processors
    call MPI_ALLREDUCE(dispoct,totoct,1,MPI_INTEGER8,MPI_MAX,spaceComm,ierr)
 else
    ierr = 1
    nboct =0
    totoct = 0
 end if

 ! new offset
 !wff%offwff=wff%offwff + nboct
 wff%offwff=wff%offwff + totoct
#endif
end subroutine xderiveRead_dp2d_mpio
!!***
