!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_tail_corrections
!! NAME
!! wvl_tail_corrections
!!
!! FUNCTION
!! Perform a minimization on the wavefunctions (especially the treatment
!! of the kinetic operator) with exponentialy decreasing functions on
!! boundaries.
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
!!      afterscfloop
!!
!! CHILDREN
!!      atmdata,calculatetailcorrection,leave_new,wrtout,xallgatherv_mpi_dp
!!      xcomm_world
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_tail_corrections(dtset, energies, etotal, mpi_enreg, occ, psps, &
     & vtrial, wvl, xcart)

 use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: etotal
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data),intent(inout) :: wvl
!arrays
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: xcart(3,dtset%natom)
 real(dp),intent(in),target :: vtrial(dtset%nfft)

!Local variables-------------------------------
!scalars
 integer :: iatom,ierr,nbuf,nsize,spaceComm,vtrial_shift
 real(dp) :: amu,ekin_sum,epot_sum,eproj_sum,rcov
 logical :: parallel
 character(len=500) :: message
!arrays
 integer :: ntails(3)
 real(dp) :: atails(3)
 real(dp),pointer :: vtotal(:)
 character(len=20) :: atomnames(100)

! *************************************************************************
 
 parallel = (mpi_enreg%nproc > 1)

!Write a message with the total energy before tail corrections.
 etotal = energies%e_kinetic + energies%e_localpsp + energies%e_nonlocalpsp + &
& energies%e_hartree + energies%e_xc - energies%e_vxc + &
& energies%e_ewald + energies%e_corepsp
 write(message,'(a,2x,e19.12)') ' Total energy before tail correction', etotal
 call wrtout(06, message, 'COLL')

!Calculate kinetic energy correction due to boundary conditions
 nbuf = nint(dtset%tl_radius / dtset%wvl_hgrid)
 ntails = dtset%wvl_internal%nSize + 2 * nbuf
 atails = real(ntails, dp) * dtset%wvl_hgrid
 write(message,'(a,a,i6,a,A,A,3F12.6,A,A,3I12,A)') ch10,&
& ' Tail requires ',nbuf,' additional grid points around cell.', ch10, &
& '  | new acell:', atails, ch10, &
& '  | new box size for wavelets:', ntails, ch10
 call wrtout(6,message,'COLL')
 call wrtout(ab_out,message,'COLL')

!---reformat potential
 if (parallel) then
  allocate(vtotal(product(dtset%wvl_internal%dpSize)))
  nsize = dtset%wvl_internal%dpSize(1) * dtset%wvl_internal%dpSize(2)
  vtrial_shift = 1 + dtset%wvl_internal%dpSize(1) * &
&  dtset%wvl_internal%dpSize(2) * mpi_enreg%nscatterarr(mpi_enreg%me, 4)
  call xcomm_world(mpi_enreg, spaceComm)
  call xallgatherv_mpi_dp(vtrial(vtrial_shift:dtset%nfft), &
&  nsize * mpi_enreg%nscatterarr(mpi_enreg%me, 2), &
&  vtotal,  nsize * mpi_enreg%nscatterarr(:,2), &
&  nsize * mpi_enreg%nscatterarr(:,3), spaceComm, ierr)
 else
  vtotal => vtrial
 end if

!Create atomnames
 do iatom = 1, dtset%ntypat, 1
  write(atomnames(iatom), "(A)") repeat(" ", 20)
  call atmdata(amu, rcov, atomnames(iatom), dtset%znucl(iatom))
 end do

#if defined HAVE_BIGDFT
 call CalculateTailCorrection(mpi_enreg%me, mpi_enreg%nproc, &
& dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
& dtset%wvl_internal%nSize(3), dtset%tl_radius, &
& wvl%wfs%nstates, wvl%wfs%mbandp, dtset%natom, dtset%ntypat, &
& dtset%wvl_internal%fineGrid(1, 1), dtset%wvl_internal%fineGrid(2, 1), &
& dtset%wvl_internal%fineGrid(1, 2), dtset%wvl_internal%fineGrid(2, 2), &
& dtset%wvl_internal%fineGrid(1, 3), dtset%wvl_internal%fineGrid(2, 3), &
& wvl%wfs%keys, wvl%projectors%keys, dtset%tl_nprccg, &
& psps%gth_params%psppar, psps%pspcod, wvl%wfs%eval, &
& vtotal, dtset%wvl_hgrid, xcart, psps%gth_params%radii_cf, &
& dtset%wvl_crmult, dtset%wvl_frmult, dtset%typat, atomnames, dtset%nsppol, &
& wvl%wfs%spinar, wvl%projectors%proj, wvl%wfs%psi, occ, &
& .false., parallel, ekin_sum, epot_sum, eproj_sum)
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_tail_corrections: BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif

 if (parallel) deallocate(vtotal)

 energies%e_kinetic = ekin_sum
 energies%e_localpsp = epot_sum - real(2., dp) * energies%e_hartree
 energies%e_nonlocalpsp = eproj_sum
 energies%e_corepsp = real(0., dp)
 etotal = energies%e_kinetic + energies%e_localpsp + energies%e_nonlocalpsp + &
& energies%e_hartree + energies%e_xc - energies%e_vxc + &
& energies%e_ewald + energies%e_corepsp

 write(message,'(a,3(1x,e18.11))') ' ekin_sum,epot_sum,eproj_sum',  & 
 ekin_sum,epot_sum,eproj_sum
 call wrtout(06, message, 'COLL')
 write(message,'(a,3(1x,e18.11))') ' ehart,eexcu,vexcu', &
& energies%e_hartree,energies%e_xc,energies%e_vxc
 call wrtout(06, message, 'COLL')
 write(message,'(a,2x,e19.12)') ' Total energy with tail correction', etotal
 call wrtout(06, message, 'COLL')

!--- End if of tail calculation
end subroutine wvl_tail_corrections
!!***
