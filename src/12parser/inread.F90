!{\src2tex{textfont=tt}}
!!****f* ABINIT/inread
!! NAME
!! inread
!!
!! FUNCTION
!! Carry out internal read from input character string, starting
!! at first character in string, reading ndig digits (including possible
!! sign, decimal, and exponent) by computing the appropriate format and
!! performing a formatted read (list-directed read would be perfect for
!! this application but is inconsistent with internal read according to
!! Fortran90 standard).
!! In case of a real number, this routine
!! is also able to read SQRT(number): return the square root of the number.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string=character string.
!!  ndig=length of field to be read (including signs, decimals, and exponents).
!!  typevarphys=variable type (might indicate the physical meaning of
!!   for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR','LEN','ENE'=>real(dp) (no special treatment)
!!   'LOG'=>integer, but read logical variable T,F,.true., or .false.
!!   'KEY'=>character, returned in token
!!
!! OUTPUT
!!  outi or outr (integer or real respectively)
!!  errcod, =0 for success, 1,2 for ini, inr failure resp.
!!
!! PARENTS
!!      adini,inarray
!!
!! CHILDREN
!!      leave_new
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine inread(string,ndig,typevarphys,outi,outr,errcod)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndig
 integer,intent(out) :: errcod,outi
 real(dp),intent(out) :: outr
 character(len=*),intent(in) :: string
 character(len=3),intent(in) :: typevarphys

!Local variables-------------------------------
!scalars
 integer :: done,idig,index_slash
 real(dp) :: den,num
 logical :: logi

! *************************************************************************

!DEBUG
!write(6,*)' inread : enter '
!write(6,*)'  string(1:ndig)=',string(1:ndig)
!ENDDEBUG

 if (typevarphys=='INT') then

! integer input section
  read (unit=string(1:ndig),fmt=*,iostat=errcod) outi
  if(errcod/=0)then
!  integer reading error
   write(6, '(/,a,/,a,i6,a)' ) &
&   ' inread : ERROR -',&
&   '  Attempted to read ndig=',ndig,' integer digits,'
   write(06, '(a,a,a)' ) '   from string(1:ndig)= ',string(1:ndig),&
&   ', to initialize an integer variable'
   errcod=1
  end if

 else if (typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE') then

! real(dp) input section

! Special treatment of SQRT(xxx) chain of characters
  done=0
  if (ndig>5) then
   if(string(1:5)=='SQRT(' .and. string(ndig:ndig)==')')then
    done=1
    read (unit=string(6:ndig-1),fmt=*,iostat=errcod) outr
    if(outr<-tol12)then
     errcod=1
    else
     outr=sqrt(outr)
    end if
   end if
  end if

! Special treatment of -SQRT(xxx) chain of characters
  if(done==0 .and. ndig > 6 )then
   if(string(1:6)=='-SQRT(' .and. string(ndig:ndig)==')')then
    done=1
    read (unit=string(7:ndig-1),fmt=*,iostat=errcod) outr
    if(outr<-tol12)then
     errcod=1
    else
     outr=-sqrt(outr)
    end if
   end if
  end if

! Special treatment of fractions
  if(done==0)then
   index_slash=index(string(1:ndig),'/')
   if(index_slash/=0)then
    done=1
    read (unit=string(1:index_slash-1),fmt=*,iostat=errcod) num
    if(errcod==0)then
     read (unit=string(index_slash+1:ndig),fmt=*,iostat=errcod) den
     if(errcod==0)then
      if(abs(den)<tol12)then
       errcod=1
      else
       outr=num/den
      end if
     end if
    end if
   end if
  end if

! Normal treatment of floats
  if(done==0)then ! Normal treatment of float numbers
   read (unit=string(1:ndig),fmt=*,iostat=errcod) outr
  end if

! Treatment of errors
  if(errcod/=0)then
!  real(dp) data reading error
   write(06, '(/,a,/,a,i6,a)' ) &
&   ' inread : ERROR -',&
&   '  Attempted to read ndig=',ndig,' floating point digits,'
   write(06, '(a,a,a)' ) '   from string(1:ndig) ',string(1:ndig),&
&   ', to initialize a floating variable.'
   errcod=2
  end if

 else if (typevarphys=='LOG') then

  read (unit=string(1:ndig),fmt=*,iostat=errcod) logi
  if(errcod/=0)then
!  integer reading error
   write(6, '(/,a,/,a,i6,a)' ) &
&   ' inread : ERROR -',&
&   '  Attempted to read ndig=',ndig,' integer digits,'
   write(06, '(a,a,a)' ) '   from string(1:ndig)= ',string(1:ndig),&
&   ', to initialize a logical variable.'
   errcod=3
  end if
  if(logi)outi=1
  if(.not.logi)outi=0

 else

  write(06, '(/,a,/,a,a,a,a,a)' ) &
&  ' inread: BUG -',&
&  '  Argument typevarphys must be INT,DPR,LEN,ENE or LOG ',ch10,&
&  '  but input value was',typevarphys,'.'
  call leave_new('COLL')

 end if

 if(errcod /= 0)then
  do idig=1,ndig
   if( string(idig:idig) == 'O' )then
    write(06, '(/,a,/,a,a,a)' ) &
&    '  inread : WARNING -',&
&    '   Note that this string contains the letter O. ',ch10,&
&    '   It is likely that this letter should be replaced by the number 0.'
    exit
   end if
  end do
 end if

end subroutine inread
!!***
