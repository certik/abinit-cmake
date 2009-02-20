!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffWriteDataRecInt
!! NAME
!! WffWriteDataRecInt
!!
!! FUNCTION
!! This subroutine writes integer data in
!! one record of a wavefunction file
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (XG,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! intarray=array of integer numbers
!! nn=size of the array to be written
!! wff= structured info about the wavefunction file
!!
!! OUTPUT
!! ierr=error code
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      vtowfk3
!!
!! CHILDREN
!!      xderivewrecend,xderivewrecinit,xderivewrite
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffWriteDataRecInt(intarray,ierr,nn,wff)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wffile_type), intent(inout) :: wff
 integer, intent(in) ::  nn
 integer, intent(out) ::  ierr
 integer, intent(in) :: intarray(nn)

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if( wff%accesswff == 0   .or.                     &
&   (wff%accesswff ==-1 .and. wff%master==wff%me) ) then

  write(wff%unwff,iostat=ierr)intarray(1:nn)
#if defined MPI_IO
           else if(wff%accesswff==1)then
            call xderiveWRecInit(wff,ierr)
            call xderiveWrite(wff,intarray,nn,ierr)
            call xderiveWRecEnd(wff,ierr)
            
#endif
 end if

end subroutine WffWriteDataRecInt
!!***
