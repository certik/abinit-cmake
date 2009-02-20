!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_nl_gradient
!! NAME
!! wvl_nl_gradient
!!
!! FUNCTION
!! Compute the non local part of the wavefunction gradient.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      afterscfloop,vtorho
!!
!! CHILDREN
!!      leave_new,nonlocal_forces,projectors_derivatives,wrtout,xcomm_world
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_nl_gradient(dtset, grnl, mpi_enreg, occ, psps, rprimd, wvl, xcart)

 use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data),intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),rprimd(3,3)
 real(dp),intent(in) :: xcart(3,dtset%natom)
 real(dp),intent(inout) :: grnl(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ierr,igeo,spaceComm
 character(len=500) :: message
!arrays
 real(dp),allocatable :: gxyz(:,:)

! *************************************************************************

!Compute forces
 write(message, '(a,a)' ) 'wvl_nl_gradient(): compute non-local part to gradient.'
 call wrtout(6,message,'COLL')

!Nullify output arrays.
 grnl(:, :) = zero
 
 allocate(gxyz(3, dtset%natom))
 gxyz(:,:) = zero
!Add the nonlocal part of the forces to grtn (BigDFT routine)
#if defined HAVE_BIGDFT
!If derivatives of projectors have not been computed yet, then we
!do it now.
 if (.not.associated(wvl%projectors%der)) then
  allocate(wvl%projectors%der(3 * wvl%projectors%keys%nprojel))
! the calculation of the derivatives of the projectors has been decoupled
! from the one of nonlocal forces, in this way forces can be calculated
! diring the wavefunction minimization if needed.
  call projectors_derivatives(mpi_enreg%me, dtset%wvl_internal%nSize(1), &
&  dtset%wvl_internal%nSize(2), dtset%wvl_internal%nSize(3), &
&  dtset%ntypat, dtset%natom, wvl%wfs%nstates, dtset%typat, &
&  psps%gth_params%psppar, wvl%projectors%keys, wvl%projectors%proj, &
&  xcart, psps%gth_params%radii_cf, &
&  dtset%wvl_cpmult, dtset%wvl_fpmult, dtset%wvl_hgrid, wvl%projectors%der)
 end if

 call nonlocal_forces(mpi_enreg%me, dtset%ntypat, dtset%natom, wvl%wfs%nstates, &
& wvl%wfs%mbandp, dtset%typat, psps%gth_params%psppar, psps%pspcod, occ, &
& wvl%projectors%keys, wvl%projectors%proj, wvl%projectors%der, wvl%wfs%keys, &
& wvl%wfs%psi, gxyz)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_nl_gradient: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif

 if (mpi_enreg%nproc > 1) then
  call xcomm_world(mpi_enreg, spaceComm)
  call xsum_mpi(gxyz, spaceComm, ierr)
 end if

!Forces should be in reduced coordinates.
 do ia = 1, dtset%natom, 1
  do igeo = 1, 3, 1
   grnl(igeo, ia) = - rprimd(1, igeo) * gxyz(1, ia) - &
&   rprimd(2, igeo) * gxyz(2, ia) - &
&   rprimd(3, igeo) * gxyz(3, ia)
  end do
 end do
 deallocate(gxyz)

end subroutine wvl_nl_gradient
!!***
