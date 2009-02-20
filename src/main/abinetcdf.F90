!{\src2tex{textfont=tt}}
!!****p* ABINIT/abinetcdf
!! NAME
!! abinetcdf
!!
!! FUNCTION
!! Tests if NetCDF support is working properly.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (JPMinet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! OUTPUT
!!  (print all)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program abinetcdf

 use defs_basis
#if defined HAVE_NETCDF
 use netcdf
#endif

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!scalars
 integer :: EnID,dim_len,i,kxID,ncid,status
 character(len=8) :: att_value,dim_name
!arrays
 real(dp) :: energy(10),energy_read(10)

! *************************************************************************

 call random_number(energy)

#if defined HAVE_NETCDF
 write (*,'(2a)') "We are using NetCDF version ", NF90_INQ_LIBVERS()

!Creating dataset with one variable of dimenssion 10, with one attribute

 status=NF90_CREATE(path="test.nc", cmode=NF90_CLOBBER, ncid=ncid)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_DEF_DIM(ncid, "kx", 10, kxID)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_DEF_VAR(ncid, "En", NF90_DOUBLE, kxID, EnID)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_PUT_ATT(ncid, EnID, "units", "anything")
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_ENDDEF(ncid)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_PUT_VAR(ncid, EnID, energy)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_CLOSE(ncid)
 if (status /= NF90_NOERR) call handle_err(status)

!Reading dataset and comparing with original

 status=NF90_OPEN(path="test.nc", mode=NF90_NOWRITE , ncid=ncid)
 if (status /= NF90_NOERR) call handle_err(status)

 status=NF90_Inquire_Dimension(ncid, 1, dim_name, dim_len)
 if (status /= NF90_NOERR) call handle_err(status)
 if (dim_name /= "kx") call problem("Inquire_dimension")
 if (dim_len /= 10) call problem("Inquire_dimension")

 status=NF90_get_att(ncid, 1,"units",att_value)
 if (status /= NF90_NOERR) call handle_err(status)
 if (att_value /= "anything") call problem("get_att")

 status=NF90_get_var(ncid=ncid, varid=1, values=energy_read)
 if (status /= NF90_NOERR) call handle_err(status)
 do i=1,10
  if (energy_read(i) /= energy(i)) call problem("get var")
 end do

 status=NF90_close(ncid)
 if (status /= NF90_NOERR) call handle_err(status)

 write (*,'(a)') "Test OK: NetCDF correctly integrated in ABINIT tree :-)"
#else
 write (*,'(a)') "This build of ABINIT does not provide NetCDF support"
#endif

 end program abinetcdf
!!***

!!****f* ABINIT/handle_err
!! NAME
!! handle_err
!!
!! FUNCTION
!!
!! PARENTS
!!      abinetcdf
!!
!! CHILDREN
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine handle_err(status)

 use defs_basis
#if defined HAVE_NETCDF
  use netcdf
#endif

 implicit none

!Arguments -----------------------------------
!scalars
 integer,intent(in) :: status

!Local variables-------------------------------

! *************************************************************************

#if defined HAVE_NETCDF
  if (status /= NF90_NOERR) then
   write (*,'(a)') trim(NF90_STRERROR(status))
   stop "Stopped"
  end if
#endif
 end subroutine handle_err
!!***

!!****f* ABINIT/problem
!! NAME
!! problem
!!
!! FUNCTION
!!
!! PARENTS
!!      abinetcdf
!!
!! CHILDREN
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine problem(msg)

 use defs_basis

 implicit none

!Arguments -----------------------------------
!scalars
 character(len=*),intent(in) :: msg

!Local variables-------------------------------

! *************************************************************************

  write (*,'(2a)') "Problem with NetCDF function ", msg
  stop
 end subroutine problem
!!***
