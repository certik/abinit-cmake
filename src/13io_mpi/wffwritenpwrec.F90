!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffWriteNpwRec
!! NAME
!! WffWriteNpwRec
!!
!! FUNCTION
!! This subroutine writes the npw record of a wavefunction file
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (XG,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!! nband_disk=number of bands
!! npw=number of plane waves
!! nspinor=number of spinorial components of the wavefunctions
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      rwwf,vtowfk3
!!
!! CHILDREN
!!      xderivewrecend,xderivewrecinit,xderivewriteval
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffWriteNpwRec(ierr,nband_disk,npw,nspinor,wff)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wffile_type), intent(inout) :: wff
 integer, intent(in) :: nband_disk,npw,nspinor
 integer, intent(out) :: ierr

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if( wff%accesswff == 0   .or.                     &
&   (wff%accesswff ==-1 .and. wff%master==wff%me) ) then
  write(wff%unwff,iostat=ierr) npw,nspinor,nband_disk
#if defined MPI_IO
           else if(wff%accesswff==1)then
            call xderiveWRecInit(wff,ierr)
            call xderiveWriteVal(wff,npw)
            call xderiveWriteVal(wff,nspinor)
            call xderiveWriteVal(wff,nband_disk)
            call xderiveWRecEnd(wff,ierr)
#endif
 end if

end subroutine WffWriteNpwRec
!!***


subroutine WffWriteNpwRec_cs(ierr,mpi_enreg,nband_disk,npw,nspinor,wff)

 use defs_basis
 use defs_datatypes
#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type), intent(inout) :: wff
 integer, intent(in) :: nband_disk,npw,nspinor
 integer, intent(out) :: ierr
 type(MPI_type), intent(inout) :: mpi_enreg

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if( wff%accesswff == 0   .or.                     &
      &   (wff%accesswff ==-1 .and. wff%master==wff%me) ) then
    write(wff%unwff,iostat=ierr) npw,nspinor,nband_disk
#if defined MPI_IO
 else if(wff%accesswff==1)then
    ! Only one mpi task do this writing
    if (mpi_enreg%me_cart_2d == 0) then
       call xderiveWRecInit(wff,ierr)
       call xderiveWriteVal(wff,npw)
       call xderiveWriteVal(wff,nspinor)
       call xderiveWriteVal(wff,nband_disk)
       call xderiveWRecEnd(wff,ierr)
    end if
    !call xcast_mpi(wff%offwff,0,mpi_enreg%commcart,ierr)
    call MPI_BCAST(wff%offwff,1,MPI_INTEGER8,0,mpi_enreg%commcart,ierr)
#endif
 end if

end subroutine WffWriteNpwRec_cs
