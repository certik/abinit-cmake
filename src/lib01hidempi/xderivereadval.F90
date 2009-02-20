!{\src2tex{textfont=tt}}
!!****f* ABINIT/xderiveReadVal
!! NAME
!! xderiveReadVal
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
!!   xderiveReadVal  contains
!!               xderiveReadVal_dp  :  read double precision array
!!               xderiveReadVal_int :  read integer array
!!               xderiveReadVal_char:  read character array
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

subroutine xderiveReadVal_dp(wff,xval)

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
           integer ::  ierr
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 real(dp),intent(out) :: xval

 xval=zero ! Initialization, for the compiler
#if defined MPI_IO
	   ierr = 0
           call MPI_FILE_READ_AT(wff%fhwff, wff%offwff,xval,1 &
           & , MPI_DOUBLE_PRECISION , statux, ierr)
           wff%offwff = wff%offwff +wff%nbOct_dp
#endif

end subroutine xderiveReadVal_dp

subroutine xderiveReadVal_int(wff,xval)

 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

#if defined MPI_IO
           integer ::  ierr
           integer  :: statux(MPI_STATUS_SIZE)
#endif
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) :: xval

 xval=0 ! Initialization, for the compiler
#if defined MPI_IO
           ierr = 0
           call MPI_FILE_READ_AT(wff%fhwff, wff%offwff,xval,1 &
           & , MPI_INTEGER , statux, ierr)

           wff%offwff = wff%offwff +wff%nbOct_int
#endif

end subroutine xderiveReadVal_int


subroutine xderiveReadVal_char(wff,xval,n)

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
 character(len=*),intent(out) :: xval
 integer,intent(in) :: n

 integer :: ierr

 xval=' ' ! Initialization, for the compiler

#if defined MPI_IO
           ierr = 0
           call MPI_FILE_READ_AT(wff%fhwff, wff%offwff,xval,n &
           & , MPI_CHARACTER , statux, ierr)

           wff%offwff = wff%offwff + wff%nbOct_ch * n
#endif

end subroutine xderiveReadVal_char
!!***
