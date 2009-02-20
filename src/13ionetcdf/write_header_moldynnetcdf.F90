!{\src2tex{textfont=tt}}
!!****f* ABINIT/write_header_moldynnetcdf
!! NAME
!! write_header_moldynnetcdf
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
!!  dtfil <type(datafiles_type)>=variables related to files
!!  jdtset = dataset number
!!  natom=number of atoms in unit cell
!!  ncoord = number of coordinates
!!  nelt_strten = number of elements in stress tensor
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

subroutine write_header_moldynnetcdf(dtfil, dtset, natom, ncoord, nelt_strten)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ncoord,nelt_strten
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset

!Local variables-------------------------------
#if defined HAVE_NETCDF
 character(fnlen) :: ficname
 integer :: ncid, status
 integer :: timeDimid, TensorSymDimid, DimVectorid
 integer :: DimScalarid, DimCoordid
 integer :: NbAtomsid, E_potDimid, E_kinDimid
 integer :: StressDimid,  Posid, Celid, CellVolumeId
 integer :: Time_stepDimid,MassDimid,Atomic_NumberDimid
 integer :: idim, optcell
 integer :: PrimVectId(3)
 integer :: jdtset
 character(len=2) :: appen, chain1
 character(len=100) :: chain
#endif

! *************************************************************************

#if defined HAVE_NETCDF
 jdtset = dtset%jdtset
 if (jdtset < 10) write(appen,'(i1)')jdtset
 if (jdtset >= 10) write(appen,'(i2)')jdtset

 optcell = dtset%optcell

 ficname = trim(dtfil%filnam_ds(2))//'_moldyn'//trim(appen)//'.nc'

!Creating file netcdf
 status = nf90_create(ficname, NF90_CLOBBER , ncid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Defining dimension
!Dimension time for netcdf (time dim is unlimited)
 status = nf90_def_dim(ncid, "time", nf90_unlimited, timeDimid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!symetric Tensor Dimension
 status = nf90_def_dim(ncid, "DimTensor", nelt_strten, TensorSymDimid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)


!Coordinates Dimension
 status = nf90_def_dim(ncid, "DimCoord", ncoord, DimCoordid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Atoms Dimensions
 status = nf90_def_dim(ncid, "NbAtoms", natom, NbAtomsid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Vector Dimension
 status = nf90_def_dim(ncid, "DimVector", 3 , DimVectorid )
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Scalar Dimension
 status = nf90_def_dim(ncid, "DimScalar", 1 , DimScalarid )
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Defining variables
!Time step
 status = nf90_def_var(ncid, "Time_step", nf90_double , &
& DimScalarid, Time_stepDimid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

 status = nf90_put_att(ncid, Time_stepDimid, "units", "atomic time unit")
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Ionic masses
 status = nf90_def_var(ncid, "Ionic_Mass", nf90_double , &
& NbAtomsid, MassDimid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

 status = nf90_put_att(ncid, MassDimid, "units", "atomic mass unit")
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Ionic atomic numbers
 status = nf90_def_var(ncid, "Ionic_Atomic_Number", nf90_double , &
& NbAtomsid, Atomic_NumberDimid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!E_pot
 status = nf90_def_var(ncid, "E_pot", nf90_double , &
& timeDimid, E_potDimid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

 status = nf90_put_att(ncid, E_potDimid, "units", "hartree")
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!E_kin
 status = nf90_def_var(ncid, "E_kin", nf90_double , &
& timeDimid, E_kinDimid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_att(ncid, E_kinDimid, "units", "hartree")
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Stress
 status = nf90_def_var(ncid, "Stress", nf90_double , &
& (/  TensorSymDimid , timeDimid  /), StressDimid)

 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_att(ncid, StressDimid, "units", "hartree/bohr^3")
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Position
 status = nf90_def_var(ncid, "Position", nf90_double , &
& (/ DimCoordid, NbAtomsid, timeDimid /), Posid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_att(ncid, Posid, "units", "bohr")
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!Celerity
 status = nf90_def_var(ncid, "Celerity", nf90_double , &
& (/ DimCoordid, NbAtomsid, timeDimid /), Celid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
 status = nf90_put_att(ncid, Celid, "units", "bohr/(atomic time unit)")
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

!In case of volume cell constant
!Primitive vectors
 if ( optcell == 0) then
  do idim = 1, 3
   write(chain1,'(i1)')idim
   chain = "PrimitiveVector"//trim(chain1)
   status = nf90_def_var(ncid, trim(chain), nf90_double , DimVectorid &
&   , PrimVectId(idim))
   if ( status /= nf90_NoErr) call handle_err_netcdf(status)
  end do

! Cell Volume
  status = nf90_def_var(ncid, "Cell_Volume", nf90_double , &
&  DimScalarid, CellVolumeId)
  if ( status /= nf90_NoErr) call handle_err_netcdf(status)
  status = nf90_put_att(ncid, CellVolumeId, "units", "bohr^3")
  if ( status /= nf90_NoErr) call handle_err_netcdf(status)

 end if

!Leaving define mode
 status = nf90_enddef(ncid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)

 status = nf90_close(ncid)
 if ( status /= nf90_NoErr) call handle_err_netcdf(status)
#endif

end subroutine write_header_moldynnetcdf
!!***

