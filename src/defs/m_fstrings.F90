!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fstring
!! NAME
!!  m_fstrings
!!
!! FUNCTION
!!  This module contains basic tools to operate on Fortran strings
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!@ ABIDOC
MODULE m_fstrings

 use defs_basis

 implicit none
 
 public ::      &
&  is_letter,    &  ! Returns .TRUE. if ch is a letter and .FALSE. otherwise
&  is_digit,     &  ! Returns .TRUE. if ch is a digit (0,1,...,9) and .FALSE. otherwise
&  upper,        &  ! Convert lower case letters to UPPER CASE
&  toupper,      &  ! Convert lower case letters to UPPER CASE (function version)
&  lower,        &  ! Convert UPPER CASE letters to lower case 
&  tolower,      &  ! Convert UPPER CASE letters to lower case  (function version)
&  compact,      &  ! Converts multiple spaces and tabs to single spaces; deletes control characters and initial spaces
&  removesp,     &  ! Removes spaces, tabs, and control characters in string str
&  tovalue,      &  ! Converts number string to a number
&  lstrip,       &  ! Remove leading spaces from string
&  write_num,    &  ! Writes a number to a string using format fmt
&  trimzero,     &  ! Deletes nonsignificant trailing zeroes from a number string. 
&  writeq,       &  ! Writes a string of the form <name> = value to unit
&  match,        &  ! Find the position of a delimiter in a string
&  OPERATOR(.jn.)   ! Concatenation operator for Fortran strings
 !TODO method to center a string
!@END ABIDOC

 interface tovalue
  module procedure value_int
  module procedure value_dp
 end interface 

 interface write_num
  module procedure write_rdp_0D
  module procedure write_int_0D
 end interface

 interface writeq
  module procedure writeq_rdp_0D
  module procedure writeq_int_0D
 end interface

 interface is_digit
  module procedure is_digit_0D
 end interface

 interface operator (.jn.)    
  module procedure str_conct
 end interface 

 PRIVATE
  character(len=1),parameter :: BLANK=' ' 
  integer,parameter :: ASCII_A=ICHAR('A')
  integer,parameter :: ASCII_Z=ICHAR('Z')
  integer,parameter :: ASCII_aa=ICHAR('a')
  integer,parameter :: ASCII_zz=ICHAR('z')
  integer,parameter :: SHIFT=ASCII_A-ASCII_aa
  integer,parameter :: ASCII_0=ICHAR('0')
  integer,parameter :: ASCII_9=ICHAR('9')

CONTAINS  !===========================================================


!!***
!! NAME
!!  is_letter 
!!
!! FUNCTION
!!  Returns .TRUE. if ch is a letter and .FALSE. otherwise

function is_letter(ch) result(ans)

!Arguments ------------------------------------

 character(len=1),intent(in) :: ch
 logical :: ans
! *********************************************************************

 select case (ICHAR(ch))
 case (ASCII_A:ASCII_Z,ASCII_aa:ASCII_zz)
  ans=.TRUE.
 case DEFAULT
  ans=.FALSE.
 end select

end function is_letter
!!***

!!***
!! NAME
!!  is_digit
!!
!! FUNCTION
!!  Returns .TRUE. if ch is a digit (0,1,...,9) and .FALSE. otherwise

function is_digit_0D(ch) result(ans)

!Arguments ------------------------------------

 character(len=1),intent(in) :: ch
 logical :: ans
! *********************************************************************

 select case (ICHAR(ch))
 case(ASCII_0:ASCII_9)
  ans=.TRUE.
 case default
  ans=.FALSE.
 end select

end function is_digit_0D
!!***

!!***
!! NAME
!!  upper 
!!
!! FUNCTION
!!  Convert lower case letters to UPPER CASE 
!!

subroutine upper(str)

!Arguments ------------------------------------
!scalars

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: ic,iasc
! ********************************************************************* 

 do ic=1,LEN_TRIM(str)
  iasc=IACHAR(str(ic:ic))
  if (iasc>=ASCII_A.and.iasc<=ASCII_Z) str(ic:ic)=ACHAR(iasc-SHIFT)
 end do

end subroutine upper
!!***

!!***
!! NAME
!!  toupper 
!!
!! FUNCTION
!!  Convert lower case letters to UPPER CASE (function version)
!!

function toupper(str_in) result(str_out)

!Arguments ------------------------------------
!scalars

 character(len=*),intent(in) :: str_in
 character(len=LEN_TRIM(str_in)) :: str_out

!Local variables-------------------------------
 integer :: ic,iasc
! ********************************************************************* 

 do ic=1,LEN_TRIM(str_in)
  iasc=IACHAR(str_in(ic:ic))
  if (iasc>=ASCII_A.and.iasc<=ASCII_Z) then 
   str_out(ic:ic)=ACHAR(iasc-SHIFT)
  else 
   str_out(ic:ic)=str_in(ic:ic)
  end if
 end do

end function toupper
!!***


!!***
!! NAME
!!  lower 
!!
!! FUNCTION
!!  Convert UPPER CASE letters to lower case
!!  

subroutine lower(str)

!Arguments ------------------------------------

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: ic,iasc
! *********************************************************************

 do ic=1,LEN_TRIM(str)
  iasc=IACHAR(str(ic:ic))
  if (iasc>=ASCII_aa.and.iasc<=ASCII_zz) str(ic:ic)=ACHAR(iasc+SHIFT)
 end do

end subroutine lower
!!***

!!***
!! NAME
!!  tolower 
!!
!! FUNCTION
!!  Convert UPPER CASE letters to lower case (function version)
!!  

function tolower(str_in) result (str_out)

!Arguments ------------------------------------

 character(len=*),intent(in) :: str_in
 character(len=LEN_TRIM(str_in)) :: str_out

!Local variables-------------------------------
 integer :: ic,iasc
! *********************************************************************

 do ic=1,LEN_TRIM(str_in)
  iasc=IACHAR(str_in(ic:ic))
  if (iasc>=ASCII_aa.and.iasc<=ASCII_zz) then 
   str_out(ic:ic)=ACHAR(iasc+SHIFT)
  else 
   str_out(ic:ic)=str_in(ic:ic)
  end if
 end do

end function tolower
!!***


!!***
!! NAME
!!  compact  
!!
!! FUNCTION
!! Converts multiple spaces and tabs to single spaces; 
!! deletes control characters; removes initial spaces.
!!  

subroutine compact(str)

!Arguments ------------------------------------

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: isp,i,k,lenstr,ich
 character(len=1):: ch
 character(len=LEN_TRIM(str)):: outstr
! *********************************************************************

 str=ADJUSTL(str) ; lenstr=LEN_TRIM(str)

 outstr=BLANK ; isp=0 ; k=0
 do i=1,lenstr
  ch=str(i:i) ; ich=IACHAR(ch)
   
  select case(ich)
  case(9,32)     ! space or tab character
   if (isp==0) then
    k=k+1
    outstr(k:k)=' '
   end if
   isp=1
  case(33:)      ! not a space, quote, or control character
   k=k+1
   outstr(k:k)=ch
   isp=0
  end select
 end do

 str=ADJUSTL(outstr)

end subroutine compact
!!***


!!***
!! NAME
!!  removesp 
!!
!! FUNCTION
!!  Removes spaces, tabs, and control characters in string str
!!  

subroutine removesp(str)

!Arguments ------------------------------------

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: i,k,lenstr,ich
 character(len=1):: ch
 character(len=LEN_TRIM(str)):: outstr
! *********************************************************************

 str=ADJUSTL(str) ; lenstr=LEN_TRIM(str)

 outstr=BLANK ; k=0
 do i=1,lenstr
  ch=str(i:i)
  ich=IACHAR(ch)
  select case(ich)    
  case(0:32)  ! space, tab, or control character
   CYCLE       
  case(33:)  
   k=k+1
   outstr(k:k)=ch
  end select
 end do

 str=ADJUSTL(outstr)

end subroutine removesp
!!***


!!***
!! NAME
!!  tovalue 
!!
!! FUNCTION
!!  Convert number string to a number
!!  

subroutine value_dp(str,rnum,ios)

!Arguments ------------------------------------

 character(len=*),intent(in) ::str
 integer,intent(out) :: ios
 real(dp),intent(out) :: rnum

!Local variables-------------------------------
 integer :: ilen,ipos
! *********************************************************************

 ilen=LEN_TRIM(str) ; ipos=SCAN(str,'Ee')
 if (.not.is_digit(str(ilen:ilen)).and.ipos/=0) then
  ios=3
  RETURN
 end if
 read(str,*,iostat=ios) rnum

end subroutine value_dp


subroutine value_int(str,inum,ios)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 character(len=*),intent(in) ::str
 integer,intent(out) :: inum,ios

!Local variables-------------------------------
 real(dp) :: rnum
! *********************************************************************

 call value_dp(str,rnum,ios)
 if (ABS(rnum)>HUGE(inum)) then
  ios=1
  RETURN
 end if
 inum=NINT(rnum)

end subroutine value_int 
!!***


!!***
!! NAME
!!  lstrip 
!!
!! FUNCTION
!!  Removes leading spaces from string
!!  

subroutine lstrip(str)

!Arguments ------------------------------------

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: ii,jj,lg
! *********************************************************************

 lg=LEN(str)
 do ii=1,lg
  if (str(ii:ii)/=BLANK) EXIT
 end do

 do jj=1,lg-ii+1
  str(jj:jj)=str(ii:ii) ; ii=ii+1
 end do
 str(ii:lg)=BLANK

end subroutine lstrip

!!***
!! NAME
!!  write_num 
!!
!! FUNCTION
!!  Writes a number to a string using format fmt
!

subroutine write_rdp_0d(rnum,str,fmt)

!Arguments ------------------------------------

 real(dp),intent(in) :: rnum
 character(len=*),intent(in) :: fmt
 character(len=*),intent(out) :: str

!Local variables-------------------------------
 character(len=LEN(fmt)+2) :: formt
! *********************************************************************

 formt='('//TRIM(fmt)//')'
 write(str,formt)rnum
 str=ADJUSTL(str)

end subroutine write_rdp_0D

subroutine write_int_0D(inum,str,fmt)

!Arguments ------------------------------------

 integer,intent(in) :: inum
 character(len=*),intent(in) :: fmt
 character(len=*),intent(out) :: str

!Local variables-------------------------------
 character(len=LEN(fmt)+2) :: formt
! *********************************************************************

 formt='('//TRIM(fmt)//')'
 write(str,formt) inum
 str=ADJUSTL(str)

end subroutine write_int_0D


!!***
!! NAME
!!  trimzero 
!!
!! FUNCTION
!! Deletes nonsignificant trailing zeroes from number string str. If number
!! string ends in a decimal point, one trailing zero is added.

! NOT sure it will work

subroutine trimzero(str)

!Arguments ------------------------------------

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: i,ipos,lstr
 character :: ch
 character(len=10) :: sexp
! *********************************************************************

 ipos=SCAN(str,'eE')
 if (ipos>0) then
  sexp=str(ipos:)
  str=str(1:ipos-1)
 end if
 lstr=LEN_TRIM(str)
 do i=lstr,1,-1
  ch=str(i:i)
  if (ch=='0') CYCLE          
  if (ch=='.') then
   str=str(1:i)//'0'
   if (ipos>0) str=TRIM(str)//TRIM(sexp)
   EXIT
  end if
  str=str(1:i)
  EXIT
 end do

 if (ipos>0) str=TRIM(str)//TRIM(sexp)

end subroutine trimzero
!!***


!!***
!! NAME
!!  writeq 
!!
!! FUNCTION
!!  Writes a string of the form <name> = value to unit

subroutine writeq_rdp_0D(unit,namestr,value,fmt)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 real(dp),intent(in) :: value
 integer,intent(in) :: unit
 character(len=*),intent(in) :: fmt
 character(len=*),intent(in) :: namestr

!Local variables-------------------------------
 character(len=32) :: tempstr
! *********************************************************************

 call write_num(value,tempstr,fmt)
 call trimzero(tempstr)
 write(unit,*)TRIM(namestr)//' = '//TRIM(tempstr)

end subroutine writeq_rdp_0D


subroutine writeq_int_0D(unit,namestr,ivalue,fmt)

!Arguments ------------------------------------


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,intent(in) :: ivalue
 integer,intent(in) :: unit
 character(len=*),intent(in) :: namestr
 character(len=*),intent(in) :: fmt

!Local variables-------------------------------
 character(len=32) :: tempstr
! *********************************************************************

 call write_num(ivalue,tempstr,fmt)
 call trimzero(tempstr)
 write(unit,*)TRIM(namestr)//' = '//TRIM(tempstr)

end subroutine writeq_int_0D
!!***

!!***
!! NAME
!!  match 
!!
!! FUNCTION
!!  Sets imatch to the position in string of the delimiter matching the delimiter
!!  in position ipos. Allowable delimiters are (), [], {}, <>.

subroutine match(str,ipos,imatch,ios)

!Arguments ------------------------------------

 character(len=*),intent(in) :: str
 integer,intent(in) :: ipos
 integer,intent(out) :: ios,imatch

!Local variables-------------------------------
 integer :: i,isum,lenstr,idelim2,istart,iend,inc
 character :: delim1,delim2,ch
! *********************************************************************

 lenstr=LEN_TRIM(str)
 delim1=str(ipos:ipos)
 select case(delim1)
  case('(')
   idelim2=IACHAR(delim1)+1
   istart=ipos+1
   iend=lenstr
   inc=1
  case(')')
   idelim2=IACHAR(delim1)-1
   istart=ipos-1
   iend=1
   inc=-1
  case('[','{','<')
   idelim2=IACHAR(delim1)+2
   istart=ipos+1
   iend=lenstr
   inc=1
  case(']','}','>')
   idelim2=IACHAR(delim1)-2
   istart=ipos-1
   iend=1
   inc=-1
  case default
   write(*,*)delim1,' is not a valid delimiter'
   RETURN
 end select

 if (istart<1 .or. istart>lenstr) then
  write(*,*) delim1,' has no matching delimiter'
  RETURN
 end if

 delim2=ACHAR(idelim2) ! matching delimiter

 isum=1
 do i=istart,iend,inc
  ch=str(i:i)
  if (ch/=delim1 .and. ch/=delim2) CYCLE
  if (ch==delim1) isum=isum+1
  if (ch==delim2) isum=isum-1
  if (isum==0) EXIT
 end do

 if (isum/=0) then
  write(*,*)delim1,' has no matching delimiter'
  ios=1
  RETURN
 end if   
 ios=0 ; imatch=i

end subroutine match
!!***


function str_conct(str1,str2) result(cnct)

!Arguments ------------------------------------

 character(len=*),intent(in) :: str1,str2
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)) :: cnct

!Local variables-------------------------------
 character(len=LEN_TRIM(str1)) :: tmp1
 character(len=LEN_TRIM(str2)) :: tmp2
! *********************************************************************

 tmp1=TRIM(str1)
 tmp2=TRIM(str2)

 cnct=str1//str2

end function str_conct 
!!***

END MODULE m_fstrings
!!***
