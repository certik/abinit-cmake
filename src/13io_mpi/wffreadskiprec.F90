!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffReadSkipRec
!! NAME
!! WffReadSkipRec
!!
!! FUNCTION
!! This subroutine move forward or backward in a Wf file by nrec records.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (ZL,DCA,XG,GMR,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nrec=number of records
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!! For the future : one should treat the possible errors of backspace
!!
!! PARENTS
!!      gstate,nstdy3,nstwf3,nstwf4,randac,rwwf,vtowfk3,wfkfermi3
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_type_size,mvrecord
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffReadSkipRec(ierr,nrec,wff)

 use defs_basis
 use defs_datatypes
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer          ,intent(out) :: ierr
 integer          ,intent(in)  :: nrec
 type(wffile_type),intent(inout)  :: wff

!Local variables-------------------------------
 integer :: irec
!no_abirules
#if defined MPI_IO
           integer :: nbOct_int
           integer :: statux(MPI_STATUS_SIZE)
           integer(abinit_offset)  :: offset,posit
           integer :: delim_record
#endif

! *************************************************************************

 ierr=0
 if( wff%accesswff==0 .or.                         &
&   (wff%accesswff==-1 .and. wff%master==wff%me)  )then

  call mvrecord(ierr,nrec,wff%unwff)

#if defined MPI_IO
           else if(wff%accesswff==1)then
            call MPI_Type_size(MPI_INTEGER,wff%nbOct_int,ierr)
            call MPI_Type_size(MPI_DOUBLE_PRECISION,wff%nbOct_dp,ierr)

            if ( nrec > 0) then
!            Move forward   nrec records
             do irec=1,nrec
              offset = wff%offwff
              wff%off_recs = wff%offwff
              call MPI_FILE_READ_AT(wff%fhwff, offset,delim_record,1 &
          &    , MPI_INTEGER , statux, ierr)
              wff%lght_recs = delim_record
              wff%offwff = offset + delim_record + 2 * wff%nbOct_int
             end do
            else
!            Move backward   -nrec records
             do irec=1,-nrec
              offset = wff%offwff-wff%nbOct_int
              call MPI_FILE_READ_AT(wff%fhwff, offset,delim_record,1 &
          &       , MPI_INTEGER , statux, ierr)
              wff%lght_recs = delim_record
              wff%offwff = wff%offwff - delim_record - 2 * wff%nbOct_int
              wff%off_recs = wff%offwff
             end do
            end if
#endif

 end if ! wff%accesswff==0,1 or -1

end subroutine WffReadSkipRec
!!***
