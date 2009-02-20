!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfsinp_scratch
!! NAME
!! wvl_wfsinp_scratch
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from input guess.
!! See wvl_wfsinp_disk() or wvl_wfsinp_reformat() from other initialisation
!! routines.
!! 
!! When initialised from scratch or from disk, wvl%wfs%[h]psi comes unallocated
!! and will be allocated inside this routine.
!! When initialised from memory (reformating), wvl%wfs%[h]psi will be reallocated.
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
!!  hdr0 <type(hdr_type)>=the header of wf, den and pot files (read from restart)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  ireadwf=1 for reading from file, 0 otherwise.
!!  mpi_enreg=informations about MPI parallelization
!!  option=1 for reading a file following ABINIT format, -1 for a BigDFT format.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>= structure with informations on wf file.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wvl <type(wvl_data)>=wavefunctions & projectors informations for wavelets.
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_wfsinp_scratch(dtset, mpi_enreg, psps, rprimd, wvl, xred)

  use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_14poisson
!End of the abilint section

  implicit none

!Arguments -------------------------------
  !scalars
  type(dataset_type), intent(in)            :: dtset
  type(MPI_type), intent(in)                :: mpi_enreg
  type(pseudopotential_type),intent(in)     :: psps
  type(wvl_data), intent(inout)             :: wvl
  !arrays
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(inout)                   :: xred(3, dtset%natom)

!Local variables-------------------------------
  character(len = 500)  :: message
  logical               :: parallel
  integer               :: i, iorb, itype, natsc, iat, ierr, spaceComm
  integer, parameter    :: nmax = 6, lmax = 3
  integer               :: neleconf(nmax, 0:lmax)
  real(dp)              :: accurex, eion
  real(dp)              :: rcov, rprb, ehomo, atomMass
  real(dp), allocatable :: xcart(:,:)
  integer, allocatable  :: iasctype(:)
  real(dp), allocatable :: rhor(:,:), vpsp(:)
  real(dp), pointer     :: kernel(:)

  character(len = 20)   :: atomnames(100)
  character(len = 2)    :: atomSymbol
  character(len = 1)    :: datacode

  write(message, '(a,a)' ) ch10,&
    &  ' wvl_wfsinp_scratch: wavefunction initialisation.'
  call wrtout(6,message,'COLL')

  parallel = (mpi_enreg%nproc > 1)

#if defined HAVE_BIGDFT
  ! Test if our pseudos have semi-core using eleconf from BigDFT...
  ! Temporary array.
  allocate(iasctype(dtset%ntypat))

  do itype = 1, dtset%ntypat, 1
     !see whether the atom is semicore or not
     call eleconf(int(psps%znucltypat(itype)), int(psps%ziontypat(itype)), &
          & atomSymbol, rcov, rprb, ehomo, neleconf, iasctype(itype))        
     write(atomnames(itype), "(A,A)") trim(adjustl(atomSymbol)), "_lda"
  end do
  natsc = 0
  do iat = 1, dtset%natom, 1
     itype = dtset%typat(iat)
     if (iasctype(itype) /= 0) natsc = natsc + 1
  enddo

  ! Store xcart for each atom
  allocate(xcart(3, dtset%natom))
  call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

  ! We grep the kernel for the Poisson solver.
  call PSolver_kernel(dtset, 2, kernel, mpi_enreg, rprimd)

  ! We allocate temporary arrays for rho and vpsp.
  !allocate ionic potential
  if (mpi_enreg%ngfft3_ionic > 0) then
     allocate(vpsp(dtset%wvl_internal%dpSize(1) * &
          & dtset%wvl_internal%dpSize(2) * mpi_enreg%ngfft3_ionic))
  else
     allocate(vpsp(1))
  end if
  call createIonicPotential(mpi_enreg%me, mpi_enreg%nproc, dtset%natom, &
       & dtset%ntypat, dtset%typat, psps%gth_params%psppar, &
       & int(psps%ziontypat), xcart, dtset%wvl_hgrid, 0.d0, &
       & dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
       & dtset%wvl_internal%nSize(3), mpi_enreg%ngfft3_ionic, &
       & mpi_enreg%nscatterarr(mpi_enreg%me, 3) + 1, kernel, vpsp, eion)

  !Allocate Charge density, Potential in real space
  if (mpi_enreg%nscatterarr(mpi_enreg%me, 1) > 0) then
     allocate(rhor(dtset%wvl_internal%dpSize(1) * &
          & dtset%wvl_internal%dpSize(2) * &
          & mpi_enreg%nscatterarr(mpi_enreg%me, 1), 1))
  else
     allocate(rhor(1, 1))
  end if

  if (parallel) then
     datacode = 'D'
  else
     datacode = 'G'
  end if
  call input_wf_diag(parallel, mpi_enreg%me, mpi_enreg%nproc, &
       & dtset%wvl_internal%fineGrid(1, 1), dtset%wvl_internal%fineGrid(2, 1), &
       & dtset%wvl_internal%fineGrid(1, 2), dtset%wvl_internal%fineGrid(2, 2), &
       & dtset%wvl_internal%fineGrid(1, 3), dtset%wvl_internal%fineGrid(2, 3), &
       & dtset%natom, natsc, wvl%wfs%nstates, wvl%wfs%mbandp, &
       & dtset%wvl_internal%nSize(1), dtset%wvl_internal%nSize(2), &
       & dtset%wvl_internal%nSize(3), wvl%wfs%mvctrp, dtset%wvl_hgrid, xcart, rhor, &
       & vpsp, wvl%wfs%keys, wvl%wfs%bounds, wvl%projectors%keys, &
       & wvl%projectors%proj, atomnames, &
       & dtset%ntypat, dtset%typat, iasctype, kernel, int(psps%znucltypat), &
       & int(psps%ziontypat), psps%gth_params%psppar, psps%pspcod, &
       & dtset%ixc, wvl%wfs%psi, wvl%wfs%psit, wvl%wfs%eval, accurex, datacode, &
       & mpi_enreg%nscatterarr, mpi_enreg%ngatherarr, dtset%nsppol, wvl%wfs%spinar)

  write(message, '(a)' ) &
       &  '  | wavefunctions have been calculated.'
  call wrtout(6,message,'COLL')

  deallocate(iasctype)
  deallocate(xcart)
  deallocate(vpsp)
  deallocate(rhor)

  !allocate hpsi array (used also as transposed)
  !allocated in the transposed way (in parallel case only) such as 
  !it can also be used as the transposed hpsi
  if (parallel) then
     allocate(wvl%wfs%hpsi(wvl%wfs%mvctrp, wvl%wfs%mbandp * mpi_enreg%nproc))
  else
     allocate(wvl%wfs%hpsi(wvl%wfs%keys%nvctr_c + 7 * wvl%wfs%keys%nvctr_f, &
          & wvl%wfs%mbandp))
  endif

#else
  write(message, '(a,a,a,a)' ) ch10,&
    &  ' wvl_wfsinp_scratch: BigDFT library is not compiled.', ch10, &
    &  '   Action, used the flag --enable-bigdft when configuring.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
#endif
end subroutine wvl_wfsinp_scratch
!!***

