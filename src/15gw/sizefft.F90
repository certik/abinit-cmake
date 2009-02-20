!{\src2tex{textfont=tt}}
!!****f* ABINIT/sizefft
!! NAME
!! sizefft
!!
!! FUNCTION
!! Calculate N, the size of the smallest allowed FFT that will encompass
!! the 2M+1 points-M..+M. (Ideally, N=2M+1, but this may not be an allowable FFT.)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  m=input M value
!!
!! OUTPUT
!!  n=output N value
!!
!! NOTES
!!  See defs_fftdata for a list of allowable sizes of FFT. 
!!
!! PARENTS
!!      setmesh
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine sizefft(m,n)

 use defs_basis
 use defs_fftdata


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: m
 integer,intent(out) :: n

!Local variables-------------------------------
!scalars
 integer,parameter :: nn=ndata
 integer :: ii,nbest
 character(len=500) :: msg

! *************************************************************************

 nbest=2*m+1

 if (nbest<2) then
  write(msg,'(4a,i8)')ch10,&
&  ' sizefft : BUG-',ch10,&
&  ' nbest = ',nbest 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 if (nbest>ifftsizes(nn)) then
  write(msg,'(4a,i8,2a)')ch10,&
&  ' sizefft : ERROR-',ch10,&
&  ' nbest = ',nbest,ch10,&
&  ' is larger than any allowable FFT'
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 do ii=nn,1,-1
  if (ifftsizes(ii)>=nbest) n=ifftsizes(ii)
 end do

end subroutine sizefft
!!***
