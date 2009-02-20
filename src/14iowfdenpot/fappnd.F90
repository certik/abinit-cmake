!{\src2tex{textfont=tt}}
!!****f* ABINIT/fappnd
!!
!! NAME
!! fappnd
!!
!! FUNCTION
!! Create the modified root name to be used for output of density, potential,
!! and geometry files. See the description of the iapp input variable.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! filnam= generic output root name
!! iapp=indicates the eventual suffix to be appended to the generic output root
!!        if 0 : no suffix to be appended (called directly from gstate)
!!        if positive : append "_TIM//iapp" (called from move or brdmin)
!!        if -1 : append "_TIM0" (called from brdmin)
!!        if -2, -3, -4, -5: append "_TIMA", ... ,"_TIMD", (called from move)
!! OUTPUT
!! filapp= filename with appended string
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      leave_new
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fappnd(filapp,filnam,iapp)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iapp
 character(len=fnlen),intent(in) :: filnam
 character(len=fnlen),intent(out) :: filapp

!Local variables-------------------------------
!scalars
 integer :: ndig
 character(len=8) :: nchar

! *************************************************************************

 if(iapp==0)then
  filapp=trim(filnam)
 else if(iapp>0)then
! Create character string for filename. Determine the number of digits in iapp.
  ndig=int(log10(dble(iapp)+0.5_dp))+1
! Make integer format field of exact size (internal write)
! for assumed nchar string of 8 characters
  write(nchar, '(i8)' ) iapp
  if (ndig>8) then
   write(06, '(/,a,/,a,/,a,/,a,i12)' ) &
&   ' fappnd : ERROR - ',&
&   '  Requested file name extension has more than the allowed 8 digits.',&
&   '  Action : resubmit the job with smaller value for ntime.',&
&   '  Value computed here was ndig=',ndig
   call leave_new('COLL')
  end if
! Concatenate into character string, picking off exact number of digits
! The potential or density label will be appended in ioarr
  filapp=trim(filnam)//'_TIM'//nchar(9-ndig:8)
 else if(iapp==-1)then
  filapp=trim(filnam)//'_TIM0'
 else if(iapp==-2)then
  filapp=trim(filnam)//'_TIMA'
 else if(iapp==-3)then
  filapp=trim(filnam)//'_TIMB'
 else if(iapp==-4)then
  filapp=trim(filnam)//'_TIMC'
 else if(iapp==-5)then
  filapp=trim(filnam)//'_TIMD'
 end if
!DEBUG
!write(6, '(a)' )filapp
!ENDDEBUG

end subroutine fappnd
!!***
