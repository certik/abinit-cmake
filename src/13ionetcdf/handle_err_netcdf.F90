!{\src2tex{textfont=tt}}
!!****f* ABINIT/handle_err_netcdf
!! NAME
!! handle_err_netcdf
!!
!! FUNCTION
!! handle netcdf error
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  status  = netcdf error status
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!    None
!! TODO
!!
!! PARENTS
!!      write_header_moldynnetcdf,write_moldynvaluenetcdf
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine handle_err_netcdf(status)

 use defs_basis
#if defined HAVE_NETCDF
 use netcdf
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: status

!Local variables-------------------------------
#if defined HAVE_NETCDF
 character(len=500) :: message
#endif

! *************************************************************************

#if defined HAVE_NETCDF
 if ( status /= nf90_NoErr) then
  message = trim(nf90_strerror(status))
  call wrtout(6, message, 'PERS')
  call leave_new('PERS')
 end if
#endif

end subroutine  handle_err_netcdf
!!***
