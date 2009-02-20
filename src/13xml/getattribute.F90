!{\src2tex{textfont=tt}}
!!****f* ABINIT/getattribute
!! NAME
!! getattribute
!!
!! FUNCTION
!! Given a XML mark-up, identifies the value of one of its attributes
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  attribute=(string) the attribute whose value must be found
!!  attributelen=length of the attribute
!!  indices_markup(3)=the three indices : position of the '<' in '<markup',
!!   position of the '>' in '<markup ...>', and
!!   position of the '<' in '</markup'
!!  strln=maximal number of characters of string
!!  string_xml*(strln)=string of characters to be searched
!!
!! OUTPUT
!!  valattrib=(string) value of the attribute
!!
!! NOTES
!!
!! PARENTS
!!      append_cml2
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine getattribute(attribute,attributelen,indices_markup,strln,string_xml,valattrib)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: attributelen,strln
 character(len=*),intent(out) :: valattrib
 character(len=attributelen),intent(in) :: attribute
 character(len=strln),intent(in) :: string_xml
!arrays
 integer,intent(in) :: indices_markup(3)

!Local variables-------------------------------
 character(len=1), parameter :: quote='"'
!scalars
 integer :: index_attribute,index_first,index_second
 character(len=500) :: message

!************************************************************************

!DEBUG
!write(6,*)' getattribute : enter'
!write(6,*)' string_xml(indices_markup(1):indices_markup(2))='
!write(6,*)'"',string_xml(indices_markup(1):indices_markup(2)),'"'
!write(6,*)' trim(attribute)='
!write(6,*)'"',trim(attribute),'"'
!stop
!ENDDEBUG

!Find the place where the attribute string appears
 index_attribute=index(string_xml(indices_markup(1):indices_markup(2)),trim(attribute)) &
& +indices_markup(1)-1

 if(index_attribute==0)then
  write(message,'(4a)')ch10,&
&  ' getattribute: BUG -',ch10,&
&  '  index_attribute == 0 '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!Find the indices of the "..."
 index_first=index(string_xml(index_attribute:indices_markup(2)),quote) &
& +index_attribute-1
 index_second=index(string_xml(index_first+1:indices_markup(2)),quote) &
& +index_first

!Check that there is something inside "..."
 if(index_second-index_first < 2)then
  write(message,'(6a,i6)') ch10,&
&  ' getattribute : BUG -',ch10,&
&  '  The difference between index_first and index_second must be 2 at least,',ch10,&
&  '  however, index_second-index_first=',index_second-index_first
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 valattrib=string_xml(index_first+1:index_second-1)

!DEBUG
!write(6,*)' getattribute : exit '
!write(6,*)' attribute = "',attribute,'"'
!ENDDEBUG

end subroutine getattribute
!!***
