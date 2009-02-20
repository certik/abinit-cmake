!{\src2tex{textfont=tt}}
!!****f* ABINIT/xderiveWriteVal
!! NAME
!! xderiveWriteVal
!!
!! FUNCTION
!! Generic subroutine functions to read/write wf files.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!   we use several procedures with the same generic name
!!   xderiveWriteVal  contains
!!               xderiveWriteVal_dp  :  write double precision array
!!               xderiveWriteVal_int :  write integer array
!!               xderiveWriteVal_char:  write character array
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

subroutine xderiveWriteVal_dp(wff,xval)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

 type(wffile_type),intent(inout) :: wff
 real(dp),intent(in) :: xval

#if defined MPI_IO
           integer ::  ierr
           integer(abinit_offset)  ::posit
           integer  :: statux(MPI_STATUS_SIZE)

           posit = wff%offwff
           call MPI_FILE_WRITE_AT(wff%fhwff, posit,xval,1 &
           &  , MPI_DOUBLE_PRECISION , statux, ierr)
           wff%offwff = wff%offwff +wff%nbOct_dp
#endif
end subroutine xderiveWriteVal_dp

subroutine xderiveWriteVal_int(wff,xval)

 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

 type(wffile_type),intent(inout) :: wff
 integer,intent(in) :: xval

#if defined MPI_IO
           integer ::  ierr
           integer(abinit_offset)  ::posit
           integer  :: statux(MPI_STATUS_SIZE)

           posit = wff%offwff
           call MPI_FILE_WRITE_AT(wff%fhwff, posit,xval,1 &
           & , MPI_INTEGER , statux, ierr)

           wff%offwff = wff%offwff +wff%nbOct_int
#endif
end subroutine xderiveWriteVal_int


subroutine xderiveWriteVal_char(wff,xval,n)

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
 character(len=*),intent(in) :: xval
 integer,intent(in) ::  n

 integer ::  ierr

 ierr=0

#if defined MPI_IO
           call MPI_FILE_WRITE_AT(wff%fhwff, wff%offwff,xval,n &
           & , MPI_CHARACTER , statux, ierr)
           wff%offwff = wff%offwff + wff%nbOct_ch * n
#endif

end subroutine xderiveWriteVal_char
!!***
