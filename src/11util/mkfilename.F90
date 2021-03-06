!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkfilename
!!
!! NAME
!! mkfilename
!!
!! FUNCTION
!! From the root (input or output) file names, produce
!! a real file name.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! character(len=fnlen):: filnam(5)=the root file names
!!  (only filnam(3) and filnam(4) are really needed)
!! get=input 'get variable', if 1, must get the file from another dataset
!! idtset=number of the dataset
!! ird=input 'iread variable', if 1, must get the file from the input root
!! jdtset_(0:ndtset)=actual index of the dataset
!! ndtset=number of datasets
!! character(len=4) :: stringfil=the string of characters to be appended
!!  e.g. '_WFK' or '_DEN'
!! character(len=9) :: stringvar=the string of characters to be appended
!!  that defines the 'get' or 'ird' variables, e.g. 'wfk' or 'ddk'
!!
!! OUTPUT
!! character(len=fnlen):: filnam_out=the new file name
!! will_read=1 if the file must be read ; 0 otherwise (ird and get were zero)
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkfilename(filnam,filnam_out,get,idtset,&
& ird,jdtset_,ndtset,stringfil,stringvar,will_read)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: get,idtset,ird,ndtset
 integer,intent(out) :: will_read
 character(len=4),intent(in) :: stringfil
 character(len=9),intent(in) :: stringvar
 character(len=fnlen),intent(out) :: filnam_out
!arrays
 integer,intent(in) :: jdtset_(0:ndtset)
 character(len=fnlen),intent(in) :: filnam(5)

!Local variables-------------------------------
!scalars
 integer :: jdtset,jget
 character(len=2) :: appen
 character(len=500) :: message
 character(len=fnlen) :: filnam_appen

! *************************************************************************

!Here, defaults if no get variable
 will_read=ird

 filnam_appen=trim(filnam(3))
 if(ndtset>0)then
  jdtset=jdtset_(idtset)
  if(jdtset<10)write(appen,'(i1)')jdtset
  if(jdtset>=10)write(appen,'(i2)')jdtset
  filnam_appen=trim(filnam_appen)//'_DS'//appen
 end if
 filnam_out=trim(filnam_appen)//trim(stringfil)

!Treatment of the multi-dataset case  (get is not relevant otherwise)
 if(ndtset/=0)then

  if(idtset+get<0)then
   write(message, '(a,a,a,a,a,a,a,a,i3,a,a,a,i3,a,a,a,a)' )ch10,&
&   ' mkfilename : ERROR -',ch10,&
&   '  The sum of idtset and get',trim(stringvar),' cannot be negative,',ch10,&
&   '  while they are idtset=',idtset,', and get',trim(stringvar),'=',get,ch10,&
&   '  Action : correct get',trim(stringvar),' in your input file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

  if(get>0 .or. (get<0 .and. idtset+get>0) )then

   if(ird/=0 .and. get/=0)then
    write(message, '(a,a,a,a,a,a,a,a,a,a,i3,a,a,a,i3,a,i3,a,a,a,a,a,a,a)' )ch10,&
&    ' driver : ERROR -',ch10,&
&    '  The input variables ird',trim(stringvar),&
&    ' and get',trim(stringvar),' cannot be',ch10,&
&    '  simultaneously non-zero, while for idtset=',idtset,',',ch10,&
&    '  they are ',ird,', and ',get,'.',ch10,&
&    '  Action : correct ird',trim(stringvar),&
&    ' or get',trim(stringvar),' in your input file.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if

   will_read=1

!  Compute the dataset from which to take the file, and the corresponding index
   if(get<0 .and. idtset+get>0) jget=jdtset_(idtset+get)
   if(get>0) jget=get
   if(jget<10)write(appen,'(i1)')jget
   if(jget>=10)write(appen,'(i2)')jget

!  Note use of output filename (filnam(4))
   filnam_out=trim(filnam(4))//'_DS'//trim(appen)//trim(stringfil)

   write(message, '(a,a,a,a,a,i3,a,a)' )&
&   ' mkfilename : get',trim(stringvar) ,'/=0, take file ',trim(stringfil),&
&   ' from output of DATASET ',jget,'.',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
  end if ! conditions on get and idtset

 end if ! ndtset/=0

end subroutine mkfilename
!!***
