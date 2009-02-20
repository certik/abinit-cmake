!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_vtorho
!! NAME
!! wvl_vtorho
!!
!! FUNCTION
!! Heart of the wavelet resolution, compute new wavefunctions mixed witf previous
!! by computing the gradient of the wavefunctions knowing the external potential.
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
!!  istep=id of the current iteration (first is 1).
!!  mpi_enreg=informations about MPI parallelization
!!  occ(dtset%mband * dtset%nsppol)=occupation numbers.
!!  proj <type(wvl_projector_type)>=projectors informations for wavelets.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  vtrial(dtset%nfft)=external potential.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_kinetic(OUT)=kinetic energy part of total energy
!!   | e_localpsp(OUT)=local pseudopotential part of total energy
!!   | e_nonlocalpsp(OUT)=nonlocal pseudopotential part of total energy
!!  residm=max value for gradient in the minimisation process.
!!  rhor(dtset%nfft)=electron density in r space
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      applylocpotkinall,applyprojectorsall,daxpy,diisstp,leave_new
!!      mpi_allreduce,orthoconstraint,orthoconstraint_p,orthon,orthon_p
!!      preconditionall,solveks,sumrho,transallwaves,untransallwaves,wrtout
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_vtorho(dtset, energies, istep, mpi_enreg, &
     & occ, proj, psps, residm, rhor, vtrial, wfs)

  use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15common
!End of the abilint section

  implicit none

!Arguments -------------------------------
  type(dataset_type), intent(inout)         :: dtset
  type(energies_type), intent(inout)        :: energies
  integer, intent(in)                       :: istep
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wvl_projectors_type), intent(in)     :: proj
  type(pseudopotential_type), intent(in)    :: psps
  real(dp), intent(inout)                   :: residm
  type(wvl_wf_type), intent(inout)          :: wfs
  real(dp), intent(in)                      :: occ(dtset%mband * dtset%nsppol)
  real(dp), intent(inout)                   :: rhor(dtset%nfft)
  real(dp), intent(in)                      :: vtrial(dtset%nfft * dtset%nspden)

!Local variables-------------------------------
  real(dp), save        :: alpha, etotal_local
  character(len = 500)  :: message
  real(dp)              :: epot_sum, ekin_sum, eproj_sum, etotal
  real(dp)              :: scprsum
  integer               :: mids, vtrial_shift, ia, igeo
  logical               :: parallel
  real(dp), allocatable :: xcart(:, :), gxyz(:, :)
  character(len = 1), save :: datacode
  
  parallel = (mpi_enreg%nproc > 1)

#if defined HAVE_BIGDFT

  write(message, '(a,a)' ) ch10,&
    &  ' wvl_vtorho: compute the new density from the trial potential.'
  call wrtout(6,message,'COLL')

  ! Initialisation of mixing parameter alpha
  if (istep == 1) then
    alpha        = real(1., dp)
    etotal_local = real(1.d100, dp)

    ! We allocate the DIIS arrays if necessary.
    if (dtset%nwfshist > 0) then
       allocate(wfs%psidst(wfs%mvctrp, wfs%mbandp * mpi_enreg%nproc, dtset%nwfshist))
       allocate(wfs%hpsidst(wfs%mvctrp, wfs%mbandp * mpi_enreg%nproc, dtset%nwfshist))
       allocate(wfs%ads(dtset%nwfshist + 1, dtset%nwfshist + 1, 3))
       wfs%ads = zero
    end if

    if (parallel) then
       datacode = 'D'
    else
       datacode = 'G'
    end if
  end if 
  ! Compute the current position in the DIIS array.
  if (dtset%nwfshist > 0) mids = modulo(istep - 1, dtset%nwfshist) + 1

  ! Apply vtrial and the projectors to the wavefubctions, computing HPsi.
  vtrial_shift = 1 + dtset%wvl_internal%dpSize(1) * dtset%wvl_internal%dpSize(2) * &
       & mpi_enreg%nscatterarr(mpi_enreg%me, 4)
  call HamiltonianApplication(parallel, datacode, mpi_enreg%me, mpi_enreg%nproc, &
       & dtset%natom, dtset%ntypat, dtset%typat, dtset%wvl_hgrid,&
       & psps%gth_params%psppar, psps%pspcod, wfs%nstates, wfs%mbandp, occ, &
       & dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
       & dtset%wvl_internal%nSize(3), dtset%wvl_internal%fineGrid(1, 1), &
       & dtset%wvl_internal%fineGrid(2, 1), dtset%wvl_internal%fineGrid(1, 2), &
       & dtset%wvl_internal%fineGrid(2, 2), dtset%wvl_internal%fineGrid(1, 3), &
       & dtset%wvl_internal%fineGrid(2, 3), wfs%keys, wfs%bounds, &
       & proj%keys, proj%proj, mpi_enreg%ngatherarr, &
       & mpi_enreg%nscatterarr(mpi_enreg%me, 2), vtrial(vtrial_shift), wfs%psi, &
       & wfs%hpsi, ekin_sum, epot_sum, eproj_sum, dtset%nsppol, wfs%spinar)

  ! WARNING! e_hartree is taken from the previous iteration as e_xc
  ! Update physical values
  energies%e_kinetic = ekin_sum
  energies%e_localpsp = epot_sum - real(2., dp) * energies%e_hartree
  energies%e_nonlocalpsp = eproj_sum
  energies%e_corepsp = real(0., dp)
  etotal = energies%e_kinetic + energies%e_localpsp + energies%e_nonlocalpsp + &
         & energies%e_hartree + energies%e_xc - energies%e_vxc + &
         & energies%e_ewald + energies%e_corepsp

  ! Precondition, minimise (DIIS or steepest descent) and ortho.
  ! Compute also the norm of the gradient.
  call hpsitopsi(istep, parallel, mpi_enreg%me, mpi_enreg%nproc, wfs%nstates, &
       & wfs%mbandp, occ, dtset%wvl_hgrid, dtset%wvl_internal%nSize(1), &
       & dtset%wvl_internal%nSize(2), dtset%wvl_internal%nSize(3), &
       & dtset%wvl_internal%fineGrid(1, 1), dtset%wvl_internal%fineGrid(2, 1), &
       & dtset%wvl_internal%fineGrid(1, 2), dtset%wvl_internal%fineGrid(2, 2), &
       & dtset%wvl_internal%fineGrid(1, 3), dtset%wvl_internal%fineGrid(2, 3), &
       & wfs%mvctrp, wfs%keys, wfs%bounds%kb, wfs%eval, dtset%wvl_nprccg, &
       & mids, dtset%nwfshist, wfs%ads, etotal, etotal_local, alpha, residm, &
       & scprsum, wfs%psi, wfs%psit, wfs%hpsi, wfs%psidst, wfs%hpsidst, &
       & dtset%nsppol, wfs%spinar)
  etotal_local = etotal

  ! Density from wavefunctions.
  call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wfs)
  
  ! Debugging messages
  write(message, '(1x,a,3(1x,1pe18.11))') 'ekin_sum,epot_sum,eproj_sum', & 
       & ekin_sum,epot_sum,eproj_sum
  call wrtout(06, message, 'COLL')
  write(message, '(1x,a,3(1x,1pe18.11))') '   ehart,   eexcu,    vexcu', &
       & energies%e_hartree,energies%e_xc,energies%e_vxc
  call wrtout(06, message, 'COLL')
  write(message, '(1x,a,i6,2x,1pe19.12,1x,1pe9.2)') 'iter,total energy,gnrm', &
       & istep,etotal,residm
  call wrtout(06, message, 'COLL')

#else
  write(message, '(a,a,a,a)' ) ch10,&
    &  ' wvl_vtorho : BigDFT library is not compiled.', ch10, &
    &  '   Action, used the flag --enable-bigdft when configuring.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
#endif
end subroutine wvl_vtorho
!!***
