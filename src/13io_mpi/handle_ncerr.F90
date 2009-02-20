!{\src2tex{textfont=tt}}
!!****f* ABINIT/handle_ncerr
!! NAME
!! handle_ncerr
!!
!! FUNCTION
!! Rudimentary Error catching for NetCDF using routines: prints out a string
!! describing what was being done when the error code showed up.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (Mver)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!   ncerr = error code. If different from the preset nf90_noerr, then error
!!   message = subroutine specified message (description of action which
!!    produced the error)
!!
!! OUTPUT
!!  (only writes)
!!
!! PARENTS
!!      cut3d,gstate,hdr_io_netcdf,ini_wf_netcdf,inwffil,ioarr,loper3,newsp
!!      outwf,outxfhist,rwwf,uderiv,wffile,wffopen,wffreadnpwrec
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine handle_ncerr(ncerr,message)

#if defined HAVE_NETCDF
 use netcdf
#endif

 implicit none

!Arguments --------------------------
 integer         ,intent(in) :: ncerr
 character(len=*),intent(in) :: message

!Local variables --------------------

#if defined HAVE_NETCDF
  if (ncerr /= nf90_noerr) then
   write (*,*) 'Error in netcdf call while : ', trim(message)
   write (*,*) trim(nf90_strerror(ncerr))
   stop
  end if
#else
  write (*,*) ' handle_ncerr : Error : NETCDF not defined at compile time. handle_ncerr should not be called.'
  stop
#endif

end subroutine handle_ncerr
!!***
