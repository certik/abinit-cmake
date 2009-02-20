!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffOffset
!! NAME
!! WffOffset
!!
!! FUNCTION
!! Tool to manage WF file in the MPI/IO case : broadcast the offset of
!! the first k-point data block
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  wff <type(wffile_type)> = structured info about the wavefunction file
!!  sender = id of the sender
!!  spaceComm = id of the space communicator handler
!!
!! OUTPUT
!!  ier = error code returned by the MPI call
!!
!! PARENTS
!!      outwf
!!
!! CHILDREN
!!      mpi_bcast,xmax_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffOffset(wff,sender,spaceComm,ier)

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
 type(wffile_type),intent(inout) :: wff
 integer          ,intent(inout) :: sender
 integer          ,intent(in)    :: spaceComm
 integer          ,intent(out)   :: ier

!Local variables ------------------------------
 integer          :: icom,ima
#if defined MPI_IO
           if ( wff%accesswff == 1) then
           call xmax_mpi(sender,icom,spaceComm,ier)
           if ( icom.ge.0)then
            ima=wff%offwff
            call MPI_BCAST(ima,1,MPI_INTEGER,icom,spaceComm,ier)
            wff%offwff=ima
           end if
          end if ! accesswff
#else
 ier = 0
#endif

end subroutine WffOffset

!!***
