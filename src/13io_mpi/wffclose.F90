!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffClose
!! NAME
!! WffClose
!!
!! FUNCTION
!! This subroutine closes a Wf file.
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
!!      conducti,gstate,loop3dte,loper3,nonlinear,nstdy3,optic,respfn,scfcv3
!!      suscep,uderiv,wannier
!!
!! CHILDREN
!!      mpi_file_close
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffClose(wff,ier)

 use defs_basis
 use defs_datatypes
#if defined HAVE_ETSF_IO
 use etsf_io
#endif
#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type), intent(inout) :: wff
 integer, intent(out) :: ier

!Local ------------------------------------
 character(len=500) :: message
#if defined HAVE_ETSF_IO
  type(etsf_io_low_error) :: error
  logical                 :: lstat
  character(len = etsf_io_low_error_len)   :: errmess
#endif

! *************************************************************************

 ier=0
 if(wff%accesswff==0)then

! All processors see a local file
  close(unit=wff%unwff)

#if defined HAVE_ETSF_IO
 else if(wff%accesswff == 3)then

   call etsf_io_low_close(wff%unwff, lstat, error_data = error)
   if (.not. lstat) then
      call etsf_io_low_error_to_str(errmess, error)
      write(message, "(A,A,A,A)") ch10, " WffClose: ERROR -", ch10, &
                                & errmess(1:min(475, len(errmess)))
      call wrtout(std_out, message, 'COLL')
      call leave_new('COLL')
   end if

#endif
 else if(wff%accesswff==-1)then

! Only the master processor see a local file
  if(wff%master==wff%me)then
   close (unit=wff%unwff)
  end if

#if defined MPI_IO
           else if(wff%accesswff==1)then
            call MPI_FILE_CLOSE(wff%fhwff,ier)
            if (wff%master==wff%me ) then
             close(unit=wff%unwff)
            end if
#endif

 end if

end subroutine WffClose
!!***
