!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfsinp_reformat
!! NAME
!! wvl_wfsinp_reformat
!!
!! FUNCTION
!! This method allocates and initialises wavefunctions with values from disk.
!! See wvl_wfsinp_scratch() or wvl_wfsinp_reformat() from other initialisation
!! routines.
!! 
!! When initialised from scratch or from disk, wvl%wfs%[h]psi comes unallocated
!! and will be allocated inside this routine.
!! When initialised from memory (reformating), wvl%wfs%[h]psi will be reallocated.
!! The projectors are also recomputed.
!!
!! The scalar arrays should be reallocated using dtset%nfft after a call to
!! this routine.
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
!!      scfcv
!!
!! CHILDREN
!!      copy_old_wavefunctions,deallocate_wfd,first_orthon,leave_new
!!      reformatmywaves,wrtout,wvl_free_type_proj,wvl_free_type_wfs
!!      wvl_init_type_proj,wvl_init_type_wfs,wvl_setboxgeometry,wvl_setngfft
!!      xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_wfsinp_reformat(acell, dtset, mpi_enreg, occ, psps, &
     & rprimd, wvl, xred, xred_old)

  use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_14wvl_wfs
!End of the abilint section

  implicit none

  !Arguments ------------------------------------
  type(dataset_type), intent(inout)      :: dtset
  type(MPI_type), intent(inout)          :: mpi_enreg
  type(pseudopotential_type), intent(in) :: psps
  type(wvl_data), intent(inout)          :: wvl
  real(dp), intent(inout)                :: acell(3), rprimd(3,3)
  real(dp), intent(in)                   :: occ(dtset%mband * dtset%nkpt * &
       & dtset%nsppol)
  real(dp), intent(inout)                :: xred_old(3, dtset%natom)
  real(dp), intent(inout)                :: xred(3, dtset%natom)

  !Local variables-------------------------------
  integer                  :: nSize_old(3)
  real(dp)                 :: wvl_hgrid_old
  real(dp), allocatable    :: xcart(:,:), xcart_old(:,:)
  real(dp), pointer        :: psi_old(:,:), eigen_old(:)
#if defined HAVE_BIGDFT
  type(wavefunctions_descriptors) :: keys_old
#endif
  character(len=20)        :: atomnames(100)
  character(len=500)       :: message
  
  write(message, '(a,a)' ) ch10,&
       & ' wvl_wfsinp_reformat: reformat the wavefunctions.'
  call wrtout(06, message, 'COLL')

#if defined HAVE_BIGDFT
  ! Convert input xred_old (reduced coordinates) to xcart_old (cartesian)
  allocate(xcart_old(3, dtset%natom))
  call xredxcart(dtset%natom, 1, rprimd, xcart_old, xred_old)

  ! Copy current to old.
  call copy_old_wavefunctions(mpi_enreg%me, mpi_enreg%nproc, wvl%wfs%nstates, &
       & wvl%wfs%mbandp, dtset%wvl_hgrid, dtset%wvl_internal%nSize(1), &
       & dtset%wvl_internal%nSize(2), dtset%wvl_internal%nSize(3), wvl%wfs%eval, &
       & wvl%wfs%keys, wvl%wfs%psi, wvl_hgrid_old, nSize_old(1), nSize_old(2), &
       & nSize_old(3), eigen_old, keys_old, psi_old)

  ! We deallocate the previous projectors.
  call wvl_free_type_proj(wvl%projectors)

  ! Deallocate old wavefunctions
  call wvl_free_type_wfs(wvl%wfs)

  ! We change the box geometry.
  call wvl_setBoxGeometry(acell, dtset, mpi_enreg, psps%gth_params%radii_cf, &
       & rprimd, xred)
  call wvl_setngfft(dtset, mpi_enreg)

  ! Reallocate them with new size.
  call wvl_init_type_wfs(dtset, mpi_enreg, psps, rprimd, wvl%wfs, xred)

  ! Recopy old eval for precond.
  wvl%wfs%eval = eigen_old
  deallocate(eigen_old)

  ! We allocate psi.
  if (mpi_enreg%nproc > 1) then
     !allocated in the transposed way such as 
     !it can also be used as a work array for transposition
     allocate(wvl%wfs%psi(wvl%wfs%mvctrp, wvl%wfs%mbandp * mpi_enreg%nproc))
  else
     allocate(wvl%wfs%psi(wvl%wfs%keys%nvctr_c + 7 * wvl%wfs%keys%nvctr_f, &
          & wvl%wfs%mbandp))
  end if
  write(message, '(a,a,a,a,I0)' ) ch10, &
    &  ' wvl_wfsinp_reformat: allocate wavefunctions,', ch10, &
    &  '  size of the compressed array per proc: ', &
    & product(shape(wvl%wfs%psi))
  call wrtout(6,message,'COLL')

  ! Convert input xred (reduced coordinates) to xcart (cartesian)
  allocate(xcart(3, dtset%natom))
  call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

  ! We transfer the old wavefunctions to the new ones.
  call reformatmywaves(mpi_enreg%me, wvl%wfs%nstates, wvl%wfs%mbandp, dtset%natom, &
       & wvl_hgrid_old, nSize_old(1), nSize_old(2), nSize_old(3), xcart_old, &
       & keys_old, psi_old, dtset%wvl_hgrid, dtset%wvl_internal%nSize(1), &
       & dtset%wvl_internal%nSize(2), dtset%wvl_internal%nSize(3), xcart, &
       & wvl%wfs%keys, wvl%wfs%psi)
  deallocate(xcart, xcart_old)

  ! We free the old descriptors and arrays.
  deallocate(psi_old)
  call deallocate_wfd(keys_old, "wvl_wfsinp_reformat")

  ! Reallocate projectors for the new positions.
  call wvl_init_type_proj(dtset, mpi_enreg, wvl%projectors, psps, rprimd, xred)

  ! Orthogonilise new wavefunctions.
  call first_orthon(mpi_enreg%me, mpi_enreg%nproc, (mpi_enreg%nproc > 1), &
       & wvl%wfs%nstates_up, wvl%wfs%nstates_dn, wvl%wfs%nstates, &
       & wvl%wfs%mbandp, wvl%wfs%keys%nvctr_c, wvl%wfs%keys%nvctr_f, wvl%wfs%mvctrp, &
       & dtset%nsppol, wvl%wfs%psi, wvl%wfs%hpsi, wvl%wfs%psit)

#else
  write(message, '(a,a,a,a)' ) ch10,&
       & ' wvl_wfsinp_reformat: BUG -',ch10,&
       & '  BigDFT is not compile. Use --enable-bigdft during configure.'
  call wrtout(06, message, 'COLL')
  call leave_new('COLL')
#endif

end subroutine wvl_wfsinp_reformat
!!***
