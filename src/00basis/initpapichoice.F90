!{\src2tex{textfont=tt}}
!!****f* ABINIT/initpapichoice
!! NAME
!! initpapichoice
!!
!! FUNCTION
!! This function  :
!!   - read in the string string the value associated with the keyword 
!! papiopt 
!!  - return the value 
!! 
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string : string where to find the keyword
!!  lenstr : length of the string string
!! 
!! OUTPUT
!!  papichoice : value of the keyword
!!  
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine initpapichoice(string, lenstr, papichoice)

 use defs_basis
 use defs_datatypes
 use defs_time


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12parser
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!string containing all input file
! papi choice for speed and time analysis  
!scalars
 integer,intent(in) :: lenstr
 integer,intent(out) :: papichoice
 character(len=*),intent(in) :: string

!Local variables-------------------------------
!scalars
 integer :: tread
 character(len=30) :: token
!arrays
 integer :: intarr(1)
 real(dp) :: dprarr(1)

! *************************************************************************

!Read papiopt 
 papiopt=0

#ifdef HAVE_PAPI
 token = 'papiopt'
 call intagm(dprarr,intarr,0,1,1,string(1:lenstr),token,tread,'INT')
 write(6,*) "papiopt lu = ", intarr(1) 
 if(tread==1) papiopt=intarr(1)
#else
 papiopt=0
#endif

 papichoice=papiopt
!write(6,*) "papiopt = ", papiopt

end subroutine initpapichoice
!!***
