#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine parsefile(filnamin, lenstr, ndtset, string)

  use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12parser
 use interfaces_13xml
!End of the abilint section

  implicit none

  character(len = *), intent(in) :: filnamin
  integer, intent(out) :: ndtset, lenstr
  character(len = strlen), intent(out) :: string

  integer :: option, marr, tread
  character(len = strlen) :: string_raw
  character(len = 30) :: token
  integer :: intarr(1)
  real(dp) :: dprarr(1)
  character(len=500) :: message




  ! Read the input file, and store the information in a long string of characters
  ! Note : this could be done only by me=0, and then string would be BCAST
!!$ !!!!!!!!!!!!!!!!!!!!!!
!!$ DC 2007-02-06: WARNING, this preprocessing option has been removed since
!!$ the variable netcdf_input is not initialized anymore (previously done
!!$ in iofn1_ncdf() that had disapeared).
!!$
!!$#if defined HAVE_NETCDF
!!$ !If the input is in a NetCDF file, a string should be rebuilt to have the
!!$ !parser believe it was plain text.
!!$ if ( netcdf_input ) then
!!$  write(message,'(a,a,a,a)') ch10,&
!!$&  ' abinit : ERROR -',ch10,&
!!$&  ' Restarting from a NetCDF file has not been implemented yet.'
!!$  call wrtout(std_out,message,'COLL')
!!$  call leave_new('COLL')
!!$ else
!!$#endif
!!$ !!!!!!!!!!!!!!!!!!!!!!

  !strlen from defs_basis module
  option=1
  call instrng (filnamin,lenstr,option,strlen,string)

  !Copy original file, without change of case
  string_raw=string

  !To make case-insensitive, map characters of string to upper case:
  call inupper(string(1:lenstr))

  !Might import data from CML file(s) into string
  !Need string_raw to deal properly with CML filenames
  call importcml(lenstr,string_raw,string,strlen)
!!$ !!!!!!!!!!!!!!!!!!!!!!
!!$#if defined HAVE_NETCDF
!!$ end if
!!$#endif
!!$ !!!!!!!!!!!!!!!!!!!!!!
  
  !6) Take ndtset from the input string, then allocate
  !the arrays whose dimensions depends only on ndtset and msym.
  
  ndtset=0 ; marr=1
  token = 'ndtset'
  call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),token,tread,'INT')
  if(tread==1) ndtset=intarr(1)
  !Check that ndtset is not negative
  if (ndtset<0 .or. ndtset>99) then
     write(message, '(a,a,a,a,i12,a,a,a,a)' ) ch10,&
          &  ' abinit : ERROR -',ch10,&
          &  '  Input ndtset must be non-negative and < 100, but was ',ndtset,ch10,&
          &  '  This is not allowed.  ',ch10,&
          &  '  Action : modify ndtset in the input file.'
     call wrtout(06,  message,'COLL')
     call leave_new('COLL')
  end if

end subroutine parsefile
