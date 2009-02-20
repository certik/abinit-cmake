!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_wvltypes
!! NAME
!! defs_wvltypes
!!
!! FUNCTION
!! This module contains definitions of all structured datatypes for the
!! wavelet part of the ABINIT package.
!!
!! List of datatypes :
!! * wvl_projectors_type : structure to store projectors for wavelets calculations.
!! * wvl_wf_type : structure to store wavefunctions for wavelets calculations.
!! * wvl_data : container for all required wavelets data.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_wvltypes

 use defs_basis
#if defined HAVE_BIGDFT
 use BigDFT_API
#endif

 implicit none

!!***

 !!****t* defs_wvltypes/wvl_projectors_type
 !! NAME
 !! wvl_projectors_type
 !!
 !! FUNCTION
 !! This type constains all physical data associated to
 !! projectors of a given system (Z and atomic positions).
 !!
 !! SOURCE
 type wvl_projectors_type
    ! These arrays are compacted arrays. One should use a nonlocal_psp_descriptors
    ! object to read these arrays and grep values on real space grid.
#if defined HAVE_BIGDFT
    type(nonlocal_psp_descriptors) :: keys
#endif

    ! Data for projectors.
    !  size (%nprojel).
    real(dp), pointer :: proj(:)
    ! Derivatives of projectors (may be unassociated).
    !  size (%nprojel * 3).
    real(dp), pointer :: der(:)

 end type wvl_projectors_type

!!***

 !!****t* defs_wvltypes/wvl_wf_type
 !! NAME
 !! wvl_wf_type
 !!
 !! FUNCTION
 !! This type constains all physical data associated to
 !! wavefunctions of a given system (Z and atomic positions).
 !!
 !! SOURCE
 type wvl_wf_type
    ! This is the total number of orbitals as used in the BigDFT code.
    ! It is equal to mband when nsppol == 1 and is equal to nelect
    ! when nsppol == 2.
    ! Number of states per polarisation is not used in non-polarised case.
    integer           :: nstates, nstates_up, nstates_dn
    ! Number of band allocated for this processor
    integer           :: mbandp
    ! Number of elements allocated for this processor
    ! Usualy keys%nseg(0) + 7 * keys%nseg(1) when only one proc.
    integer           :: mvctrp

#if defined HAVE_BIGDFT
    type(wavefunctions_descriptors) :: keys
    ! Localisation informations.
    type(convolutions_bounds)       :: bounds
#endif

    ! wavefunctions, size (mvctrp, mbandp)
    real(dp), pointer :: psi(:, :)
    ! wavefunction gradients, size (mvctrp, mbandp)
    real(dp), pointer :: hpsi(:, :)

    ! Eigenvalues, used with the preconditionner (nstates).
    ! Should be removed and replaced by the eigen vector of abinit.
    real(dp), pointer :: eval(:)

    ! Spin projection on z axis, used in the collinear spin case (nstates).
    real(dp), pointer :: spinar(:)

    ! Temporary wavefunction storage when several proc are used.
    real(dp), pointer :: psit(:, :)

    ! Temporary storage when DIIS is used to minimise the wavefunctions.
    real(dp), pointer :: psidst(:,:,:), hpsidst(:,:,:)
    real(dp), pointer :: ads(:,:,:)
 end type wvl_wf_type

!!***


 !!****t* defs_wvltypes/wvl_data
 !! NAME
 !! wvl_data
 !!
 !! FUNCTION
 !! This type is a container to limit the number of arguments in
 !! ABINIT, it should not contains attributes other than types
 !! created for wavelets.
 !!
 !! SOURCE
 type wvl_data
    ! The data associated to projectors (computed from pseudo)
    type(wvl_projectors_type) :: projectors
    ! The data associated to the wavefunctions
    type(wvl_wf_type) :: wfs
 end type wvl_data

end module defs_wvltypes
!!***
