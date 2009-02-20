!{\src2tex{textfont=tt}}
!!****f* ABINIT/write_moldynvaluenetcdf
!! NAME
!! write_moldynvaluenetcdf
!!
!! FUNCTION
!! Write NetCdf file for moldyn output
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! itime1 = position vs time in netcdf files
!! dtfil <type(datafiles_type)>=variables related to files
!! jdtset = dataset number
!! Epot = potential energy
!! Ekin = kinetic energy
!! nbat = number of atoms in unit cell
!! nb1 = size of stress tensor
!! pos = atoms positions
!! cel = atoms Celerity
!! stress = stress tensor
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!    None
!! TODO
!!
!! PARENTS
!!      moldyn
!!
!! CHILDREN
!!      handle_err_netcdf
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine write_moldynvaluenetcdf(amass, itime1, dtfil, dtset, Epot, Ekin,  nbat, &
&                 nbdir, nb1, pos, cel, stress, rprimd, ucvol)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itime1,nb1,nbat,nbdir
 real(dp),intent(in) :: Ekin,Epot,ucvol
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: amass(nbat),cel(nbdir,nbat),pos(nbdir,nbat),rprimd(3,3)
 real(dp),intent(in) :: stress(nb1)

!Local variables-------------------------------
#if defined HAVE_NETCDF
 integer :: Time_stepId,MassId,Atomic_NumberId
 integer :: E_potId, E_kinId, StressId, PosId, CelId
 integer :: CellVolumeId
 integer :: jdtset
 character(fnlen) ficname
 integer :: ncid, status, kk, idim, optcell
 character(len=2) :: appen, chain1
 character(len=100) :: chain
 integer :: PrimVectId(3)
 real(dp) :: vect(3)
#endif

! *************************************************************************

#if defined HAVE_NETCDF
 jdtset = dtset%jdtset
 if (jdtset < 10) write(appen,'(i1)')jdtset
 if (jdtset >= 10) write(appen,'(i2)')jdtset

 ficname = trim(dtfil%filnam_ds(2))//'_moldyn'//trim(appen)//'.nc'

!Opening file netcdf
 status = nf90_open(ficname, nf90_write, ncid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Time step
 status = nf90_inq_varid(ncid, "Time_step", Time_stepId)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_var(ncid, Time_stepId, dtset%dtion)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Ionic masses
 status = nf90_inq_varid(ncid, "Ionic_Mass", MassId)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_var(ncid, MassId, amass, start = (/ 1 /), count = (/ nbat /) )
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Ionic atomic numbers
 status = nf90_inq_varid(ncid, "Ionic_Atomic_Number", Atomic_NumberId)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_var(ncid, Atomic_NumberId, dtset%znucl(dtset%typat(:)), &
 	&start = (/ 1 /), count = (/ nbat /) )
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Epot
 status = nf90_inq_varid(ncid, "E_pot", E_potId)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_var(ncid, E_potId, (/ Epot /), start = (/ itime1 /), count = (/ 1 /) )
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Ekin
 status = nf90_inq_varid(ncid, "E_kin", E_kinId)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_var(ncid, E_kinId, (/ Ekin /), start = (/ itime1 /), count = (/ 1 /) )
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Stress
 status = nf90_inq_varid(ncid, "Stress", StressId)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_var(ncid, StressId, stress, start = (/1,  itime1 /), count = (/ nb1 /)  )
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Pos
 status = nf90_inq_varid(ncid, "Position", PosId)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_var(ncid, PosId, pos, start = (/1, 1, itime1/), &
& count = (/  nbdir, nbat, 1 /)  )
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Cel
 status = nf90_inq_varid(ncid, "Celerity", CelId)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_var(ncid, CelId, cel, start = (/1, 1, itime1/), &
& count = (/ nbdir, nbat, 1 /)  )

 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Primitive vector
!In case of volume cell constant and first step in time
 optcell = dtset%optcell

 if ( (optcell == 0) .and. (itime1 == 1)) then
  do idim = 1, 3
   write(chain1,'(i1)')idim
   chain = "PrimitiveVector"//trim(chain1)

   status = nf90_inq_varid(ncid, trim(chain), PrimVectId(idim) )
   if ( status /= nf90_NoErr) call handle_err_netcdf(status)

   do kk = 1, 3
    vect(kk) = rprimd(kk, idim)
   end do
   status = nf90_put_var(ncid, PrimVectId(idim), vect )
   if ( status /= nf90_NoErr) call handle_err_netcdf(status)
  end do

! Cell Volume
  status = nf90_inq_varid(ncid, "Cell_Volume" , CellVolumeId )
  if ( status /= nf90_NoErr) call handle_err_netcdf(status)

  status = nf90_put_var(ncid, CellVolumeId, ucvol )
  if ( status /= nf90_NoErr) call handle_err_netcdf(status)

 end if

 status = nf90_close(ncid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
#endif

end subroutine write_moldynvaluenetcdf
!!***
