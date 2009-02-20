!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffWriteDataRec
!! NAME
!! WffWriteDataRec
!!
!! FUNCTION
!! This subroutine writes double precision data in
!! one record of a wavefunction file
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (XG,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! dparray=array of double precision numbers
!! ndp=size of the double precision array to be read
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

subroutine WffWriteDataRec(dparray,ierr,ndp,wff)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wffile_type), intent(inout) :: wff
 integer, intent(in) ::  ndp
 integer, intent(out) ::  ierr
 real(dp), intent(in) :: dparray(ndp)

!Local variables-------------------------------

! *************************************************************************

 ierr=0
 if( wff%accesswff == 0   .or.                     &
&   (wff%accesswff ==-1 .and. wff%master==wff%me) ) then

  write(wff%unwff,iostat=ierr)dparray(1:ndp)
#if defined MPI_IO
           else if(wff%accesswff==1)then
            call xderiveWRecInit(wff,ierr)
            call xderiveWrite(wff,dparray,ndp,ierr)
            call xderiveWRecEnd(wff,ierr)
#endif
 end if

end subroutine WffWriteDataRec
!!***
