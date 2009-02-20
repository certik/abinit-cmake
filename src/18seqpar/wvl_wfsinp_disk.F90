!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfsinp_disk
!! NAME
!! wvl_wfsinp_disk
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from disk.
!! See wvl_wfsinp_scratch() or wvl_wfsinp_reformat() from other initialisation
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
!!  mpi_enreg=informations about MPI parallelization
!!  option=1 for reading a file following ABINIT format, -1 for a BigDFT format.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wff <type(wffile_type)>= structure with informations on wf file.
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
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

subroutine wvl_wfsinp_disk(dtset, hdr0, hdr, mpi_enreg, option, &
     & rprimd, wff, wfs, xred)

  use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_14wvl_wfs
!End of the abilint section

  implicit none

!Arguments -------------------------------
  !scalars
  integer, intent(in)                       :: option
  type(dataset_type), intent(in)            :: dtset
  type(hdr_type), intent(in)                :: hdr0
  type(hdr_type), intent(in)                :: hdr
  type(MPI_type), intent(in)                :: mpi_enreg
  type(wffile_type), intent(in)             :: wff
  type(wvl_wf_type), intent(inout)          :: wfs
  !arrays
  real(dp), intent(in)                      :: rprimd(3, 3)
  real(dp), intent(inout)                   :: xred(3, dtset%natom)

!Local variables-------------------------------
  character(len = 500)  :: message
  logical               :: parallel
  integer               :: spaceComm, ierr

  write(message, '(a,a)' ) ch10,&
    &  ' wvl_wfsinp_disk: wavefunction initialisation.'
  call wrtout(6,message,'COLL')

#if defined HAVE_BIGDFT
  parallel = (mpi_enreg%nproc > 1)

  ! We allocate psi.
  if (parallel) then
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(wfs%psi(wfs%mvctrp, wfs%mbandp * mpi_enreg%nproc))
  else
     allocate(wfs%psi(wfs%keys%nvctr_c + 7 * wfs%keys%nvctr_f, wfs%mbandp))
  end if
  write(message, '(a,a,a,a,I0)' ) ch10, &
    &  ' wvl_wfsinp_disk: allocate wavefunctions,', ch10, &
    &  '  size of the compressed array per proc: ', &
    & product(shape(wfs%psi))
  call wrtout(6,message,'COLL')

  call wvl_read(dtset, hdr0, hdr, mpi_enreg, option, rprimd, wff, wfs, xred)

  ! We orthogonalise.
  call first_orthon(mpi_enreg%me, mpi_enreg%nproc, parallel, wfs%nstates_up, &
       & wfs%nstates_dn, wfs%nstates, wfs%mbandp, wfs%keys%nvctr_c, &
       & wfs%keys%nvctr_f, wfs%mvctrp, dtset%nsppol, wfs%psi, wfs%hpsi, wfs%psit)
#else
  write(message, '(a,a,a,a)' ) ch10,&
    &  ' wvl_wfs_inp: BigDFT library is not compiled.', ch10, &
    &  '   Action, used the flag --enable-bigdft when configuring.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
#endif
end subroutine wvl_wfsinp_disk
!!***

