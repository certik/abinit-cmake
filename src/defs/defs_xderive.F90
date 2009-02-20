!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_xderive
!! NAME
!! defs_xderive
!!
!! FUNCTION
!! This module contains generic interfaces to read/write wf files.
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
!!               xderiveRead_dp   :  read real(dp) value
!!               xderiveRead_dp2d   :  read real(dp) array 2d
!!   xderiveReadVal  contains
!!               xderiveReadVal_dp  :  read real(dp) array
!!               xderiveReadVal_int :  read integer array
!!               xderiveReadVal_char:  read character array
!!   xderiveWrite  contains
!!               xderiveWrite_int  :  write integer  value
!!               xderiveWrite_int2d  :  write integer  array 2d
!!               xderiveWrite_dp   :  write real(dp) value
!!               xderiveWrite_dp2d   : write real(dp) array 2d
!!   xderiveWriteVal  contains
!!               xderiveWriteVal_dp  :  write real(dp) array
!!               xderiveWriteVal_int :  write integer array
!!               xderiveWriteVal_char:  write character array
!!
!! PARENTS
!!      outxfhist,rwwf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_xderive

 implicit none

!Generic interface of the routines xderiveRead
 interface xderiveRead

  subroutine xderiveRead_int(wff,xval,n1,ierr)
   use defs_datatypes
   implicit none
   integer,intent(in) :: n1
   integer,intent(out) :: ierr
   type(wffile_type),intent(inout) :: wff
   integer,intent(out) :: xval(:)
  end subroutine xderiveRead_int

  subroutine xderiveRead_int2d(wff,xval,n1,n2,ierr)
   use defs_datatypes
   implicit none
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   integer,intent(out) :: ierr
   type(wffile_type),intent(inout) :: wff
   integer,intent(out) :: xval(:,:)
  end subroutine xderiveRead_int2d

  subroutine xderiveRead_dp(wff,xval,n1,ierr)
   use defs_basis
   use defs_datatypes
   implicit none
   integer,intent(in) :: n1
   integer,intent(out) :: ierr
   type(wffile_type),intent(inout) :: wff
   real(dp),intent(out) :: xval(:)
  end subroutine xderiveRead_dp

  subroutine xderiveRead_dp2d(wff,xval,n1,n2,ierr)
   use defs_basis
   use defs_datatypes
   implicit none
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   integer,intent(out) :: ierr
   type(wffile_type),intent(inout) :: wff
   real(dp),intent(out) :: xval(:,:)
  end subroutine xderiveRead_dp2d

 end interface
!End of the generic interface of xderiveRead



!Generic interface of the routines xderiveReadVal
 interface xderiveReadVal

  subroutine xderiveReadVal_dp(wff,xval)
   use defs_basis
   use defs_datatypes
   implicit none
   type(wffile_type),intent(inout) :: wff
   real(dp),intent(out) :: xval
  end subroutine xderiveReadVal_dp

  subroutine xderiveReadVal_int(wff,xval)
   use defs_datatypes
   implicit none
   integer,intent(out) :: xval
   type(wffile_type),intent(inout) :: wff
  end subroutine xderiveReadVal_int

  subroutine xderiveReadVal_char(wff,xval,n)
   use defs_datatypes
   implicit none
   integer,intent(in) :: n
   type(wffile_type),intent(inout) :: wff
   character(len=*),intent(out) :: xval
  end subroutine xderiveReadVal_char

 end interface
!End of the generic interface of xderiveReadVal



!Generic interface of the routines xderiveWrite
 interface xderiveWrite

  subroutine xderiveWrite_int(wff,xval,n1,ierr)
   use defs_datatypes
   implicit none
   integer,intent(in) :: n1
   integer,intent(out) :: ierr
   type(wffile_type),intent(inout) :: wff
   integer,intent(in) :: xval(:)
  end subroutine xderiveWrite_int

  subroutine xderiveWrite_int_mpio(wff,xval,n1,ierr,spaceComm)
   use defs_datatypes
   implicit none
   integer,intent(out) :: ierr
   integer,intent(in) :: n1
   integer,intent(in) :: spaceComm
   type(wffile_type),intent(inout) :: wff
   integer,intent(in) :: xval(:)
  end subroutine xderiveWrite_int_mpio

  subroutine xderiveWrite_int2d(wff,xval,n1,n2,ierr)
   use defs_datatypes
   implicit none
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   integer,intent(out) :: ierr
   type(wffile_type),intent(inout) :: wff
   integer,intent(in) :: xval(:,:)
  end subroutine xderiveWrite_int2d

  subroutine xderiveWrite_int2d_mpio(wff,xval,n1,n2,ierr,spaceComm)
   use defs_datatypes
   implicit none
   integer,intent(out) :: ierr
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   integer,intent(in) :: spaceComm
   type(wffile_type),intent(inout) :: wff
   integer,intent(in) :: xval(:,:)
  end subroutine xderiveWrite_int2d_mpio

  subroutine xderiveWrite_int2d_mpio_arr(wff,xval,n1,n2,ierr,local_offset)
   use defs_datatypes
   implicit none
   integer,intent(out) :: ierr
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   type(wffile_type),intent(inout) :: wff
   integer,intent(in) :: local_offset(:)
   integer,intent(in) :: xval(:,:)
  end subroutine xderiveWrite_int2d_mpio_arr

  subroutine xderiveWrite_dp(wff,xval,n1,ierr)
   use defs_basis
   use defs_datatypes
   implicit none
   integer,intent(in) :: n1
   integer,intent(out) :: ierr
   type(wffile_type),intent(inout) :: wff
   real(dp),intent(in) :: xval(:)
  end subroutine xderiveWrite_dp

  subroutine xderiveWrite_dp_mpio(wff,xval,n1,ierr,spaceComm)
   use defs_basis
   use defs_datatypes
   implicit none
   integer,intent(out) :: ierr
   integer,intent(in) :: n1
   integer,intent(in) :: spaceComm
   type(wffile_type),intent(inout) :: wff
   real(dp),intent(in) :: xval(:)
  end subroutine xderiveWrite_dp_mpio

  subroutine xderiveWrite_dp2d_mpio_arr(wff,xval,n1,n2,ierr,local_offset)
   use defs_basis
   use defs_datatypes
   implicit none
   integer,intent(out) :: ierr
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   type(wffile_type),intent(inout) :: wff
   integer,intent(in) :: local_offset(:)
   real(dp),intent(in) :: xval(:,:)
  end subroutine xderiveWrite_dp2d_mpio_arr

  subroutine xderiveWrite_dp2d(wff,xval,n1,n2,ierr)
   use defs_basis
   use defs_datatypes
   implicit none
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   integer,intent(out) :: ierr
   type(wffile_type),intent(inout) :: wff
   real(dp),intent(in) :: xval(:,:)
  end subroutine xderiveWrite_dp2d

  subroutine xderiveWrite_dp2d_mpio(wff,xval,n1,n2,ierr,spaceComm)
   use defs_basis
   use defs_datatypes
   implicit none
   integer,intent(out) :: ierr
   integer,intent(in) :: n1
   integer,intent(in) :: n2
   integer,intent(in) :: spaceComm
   type(wffile_type),intent(inout) :: wff
   real(dp),intent(in) :: xval(:,:)
  end subroutine xderiveWrite_dp2d_mpio

 end interface
!End of the generic interface of xderiveWrite



!Generic interface of the routines xderiveWriteVal
 interface xderiveWriteVal

  subroutine xderiveWriteVal_dp(wff,xval)
   use defs_basis
   use defs_datatypes
   implicit none
   type(wffile_type),intent(inout) :: wff
   real(dp),intent(in) :: xval
  end subroutine xderiveWriteVal_dp

  subroutine xderiveWriteVal_int(wff,xval)
   use defs_datatypes
   implicit none
   integer,intent(in) :: xval
   type(wffile_type),intent(inout) :: wff
  end subroutine xderiveWriteVal_int

  subroutine xderiveWriteVal_char(wff,xval,n)
   use defs_datatypes
   implicit none
   integer,intent(in) :: n
   type(wffile_type),intent(inout) :: wff
   character(len=*),intent(in) :: xval
  end subroutine xderiveWriteVal_char

 end interface
!End of the generic interface of xderiveWriteVal

 interface WffWriteNpwRec
    subroutine WffWriteNpwRec(ierr,nband_disk,npw,nspinor,wff)
      use defs_datatypes
      implicit none
      integer, intent(out) :: ierr
      integer, intent(in) :: nband_disk
      integer, intent(in) :: npw
      integer, intent(in) :: nspinor
      type(wffile_type), intent(inout) :: wff
    end subroutine WffWriteNpwRec

    subroutine WffWriteNpwRec_mpio(ierr,mpi_enreg,nband_disk,npw,nspinor,wff)
      use defs_datatypes
      implicit none
      integer, intent(out) :: ierr
      integer, intent(in) :: nband_disk
      integer, intent(in) :: npw
      integer, intent(in) :: nspinor
      type(MPI_type), intent(inout) :: mpi_enreg
      type(wffile_type), intent(inout) :: wff
    end subroutine WffWriteNpwRec_mpio
 end interface
end module defs_xderive
!!***
