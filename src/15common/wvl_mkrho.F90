!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_mkrho
!! NAME
!! wvl_mkrho
!!
!! FUNCTION
!! This method is just a wrapper around the BigDFT routine to compute the
!! density from the wavefunctions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=input variables.
!!  mpi_enreg=informations about MPI parallelization
!!  occ(dtset%mband)=occupation numbers.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!!  rhor(dtset%nfft)=electron density in r space
!!
!! SIDE EFFECTS
!!  proj <type(wvl_projector_type)>=projectors informations for wavelets.
!!   | proj(OUT)=computed projectors.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_mkrho(dtset, mpi_enreg, occ, rhor, wfs)

 use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(wvl_wf_type),intent(in) :: wfs
!arrays
 real(dp),intent(in) :: occ(dtset%mband)
 real(dp),intent(inout) :: rhor(dtset%nfft,dtset%nspden)

!Local variables-------------------------------
!scalars
 logical :: parallel
 character(len=500) :: message

! *************************************************************************

#if defined HAVE_BIGDFT
 parallel = (mpi_enreg%nproc > 1)

 call sumrho(parallel, mpi_enreg%me, mpi_enreg%nproc, wfs%nstates, wfs%mbandp, &
& dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
& dtset%wvl_internal%nSize(3), dtset%wvl_hgrid, occ, &
& wfs%keys, wfs%psi, rhor, &
& dtset%wvl_internal%dpSize(1) * dtset%wvl_internal%dpSize(2) * &
& mpi_enreg%nscatterarr(mpi_enreg%me, 1), mpi_enreg%nscatterarr, &
& dtset%nsppol, wfs%spinar, &
& dtset%wvl_internal%fineGrid(1, 1), dtset%wvl_internal%fineGrid(2, 1), &
& dtset%wvl_internal%fineGrid(1, 2), dtset%wvl_internal%fineGrid(2, 2), &
& dtset%wvl_internal%fineGrid(1, 3), dtset%wvl_internal%fineGrid(2, 3), &
& wfs%bounds)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_mkrho : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_mkrho
!!***
