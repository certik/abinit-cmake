!{\src2tex{textfont=tt}}
!!****f* ABINIT/clsopn
!! NAME clsopn
!! clsopn
!!
!!
!! FUNCTION
!! Close wavefunction file (provided its access is standard F90 IO),
!! then reopen the same.
!! Uses fortran inquire statement to reopen with same
!! characteristics.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  wff=number of unit to which on which file is already
!!  opened.
!!
!! OUTPUT
!!
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      nselt3,nstdy3,optics_paw,outkss,outwant,partial_dos_fractions,rhofermi3
!!      vtorho,vtorho3,wffile
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine clsopn(wff)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wffile_type),intent(inout) :: wff

!Local variables-------------------------------
!scalars
 integer :: ios,unit
 logical :: nmd,od
 character(len=11) :: fm
 character(len=500) :: message
 character(len=fnlen) :: filnam

! *************************************************************************

 if(wff%accesswff<=0)then

  unit=wff%unwff
  inquire (unit=unit,iostat=ios,opened=od,name=filnam,&
&  form=fm,named=nmd)

! ios is a status specifier.  If an error condition exists,
! ios is assigned a processor-dependent value > 0.
  if (ios/=0) then
   write(message, '(/,a,/,a,i8,a,i8,/,a,/,a,/,a)' ) &
&   ' clsopn : ERROR -',&
&   '  Attempt to inquire about unit=',unit,&
&   '  indicates error condition iostat=',ios,&
&   '  May be due to temporary problem with file, disks or network.',&
&   '  Action : check whether there might be some external problem,',&
&   '  then resubmit.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')

!  od is a logical variable which is set to true if the specified
!  unit is connected to a file; otherwise it is set to false.
#if !defined FC_HITACHI
  else if (.not.od) then
   write(message, '(/,a,/,a,i8,/,a,/,a,/,a,/,a)' ) &
&   ' clsopn : ERROR -',&
&   '  Tried to inquire about unit',unit,&
&   '  and found it not connected to a file.',&
&   '  May be due to temporary problem with file, disks or network.',&
&   '  Action : check whether there might be some external problem,',&
&   '  then resubmit.'
   call wrtout(std_out,message,'COLL')
   call leave_new('COLL')
#endif

!  nmd is a logical variable assigned the value true if the file
!  has a name; otherwise false.  A scratch file is not named.
  else if (.not.nmd) then

!  No action for the time being. Possibility to debug.

  else

!  May now close the file and then reopen it
!  (file is already opened according to above checks)

#if defined FC_HITACHI
   if (.not.od) then
    write(message, '(/,a,/,a,i8,/,a,/,a,/,a)' ) &
&    ' clsopn : ERROR - (it might be a bug on SR8k sytem)',&
&    '  Tried to inquire about unit',unit,&
&    '  and found it not connected to a file.',&
&    '  May be due to temporary problem with file, disks or network.',&
&    '  Action : disregard this error and continue the process anyway.'
    call wrtout(std_out,message,'COLL')
   end if
#endif
   close (unit=unit)
   open (unit=unit,file=filnam,form=fm,status='old')

  end if

 else if (wff%accesswff == 1) then
  wff%offwff = 0
 else if (wff%accesswff == 3) then
! We do nothing, ETSF access already not being sequential.
 end if

end subroutine clsopn
!!***
