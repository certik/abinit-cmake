!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffDelete
!! NAME
!! WffDelete
!!
!! FUNCTION
!! This subroutine closes a Wf file, and delete it.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! PARENTS
!!      gstate,loper3,outwf,respfn
!!
!! CHILDREN
!!      mpi_file_close
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffDelete(wff,ier)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer, intent(out) :: ier

!Local variables-------------------------------

! *************************************************************************

 ier=0
 if(wff%accesswff==0)then

! All processors see a local file
  close(unit=wff%unwff,status='delete')

 else if(wff%accesswff==-1)then

! Only the master processor see a local file
  if(wff%master==wff%me)then
   close (unit=wff%unwff,status='delete')
  end if

#if defined MPI_IO
           else if(wff%accesswff==1)then
            if ( wff%fhwff /= -1 )then
             call MPI_FILE_CLOSE(wff%fhwff,ier)
            end if
            if (wff%master==wff%me ) then
             close(unit=wff%unwff,status='delete')
             wff%fhwff = -1
            end if
#endif

 end if

end subroutine WffDelete
!!***
