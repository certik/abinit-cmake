!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_setngfft
!! NAME
!! wvl_setngfft
!!
!! FUNCTION
!! When wavelets are used, the FFT grid is used to store potentials and
!! density. The size of the grid takes into account the two resolution in wavelet
!! description and also the distribution over processor in the parallel case.
!!
!! The FFT grid is not in strict terms an FFT grid but rather a real space grid.
!! Its dimensions are not directly compatible with FFTs. This is not relevant
!! when using the wavelet part of the code and in the Poisson solver the arrays
!! are extended to match FFT dimensions internally. But for other parts of the
!! code, this must be taken into account.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization (description of the
!!            density and potentials scatterring is allocated and updated).
!!  dtset <type(dataset_type)>=the FFT grid is changed.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      leave_new,ps_dim4allocation,wrtout
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_setngfft(dtset, mpi_enreg)

 use defs_basis
  use defs_datatypes
#if defined HAVE_BIGDFT
  use Poisson_Solver
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset

!Local variables-------------------------------
!scalars
 integer :: density_start,jproc,ngfft3_density,ngfft3_ionic,ngfft3_potential
 integer :: potential_shift
 character(len=1) :: datacode
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(6,*)' wvl_setngfft : enter '
!write(6,*)' associated(mpi_enreg%nscatterarr)=',associated(mpi_enreg%nscatterarr)
!stop
!ENDDEBUG

#if defined HAVE_BIGDFT
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_setngfft : Changing the FFT grid definition.'
 call wrtout(6,message,'COLL')

!Change nfft and ngfft
!We use the routine given in the poisson solver part.

 if (mpi_enreg%nproc > 1) then
  datacode = 'D'
 else
  datacode = 'G'
 end if

 do jproc = 0, mpi_enreg%nproc - 1, 1
! Call the data distribution for the free ('F') boundary counditions.
! We distribute ('D') the data among processors
  call PS_dim4allocation('F', datacode, jproc, mpi_enreg%nproc, &
&  dtset%wvl_internal%dpSize(1), dtset%wvl_internal%dpSize(2), &
&  dtset%wvl_internal%dpSize(3), dtset%ixc, ngfft3_density, &
&  ngfft3_potential, ngfft3_ionic, potential_shift, density_start)
  if (jproc == mpi_enreg%me) then
   mpi_enreg%ngfft3_ionic = ngfft3_ionic
  end if
! number of planes for the density
  mpi_enreg%nscatterarr(jproc, 1) = ngfft3_density
! number of planes for the potential
  mpi_enreg%nscatterarr(jproc, 2) = ngfft3_potential
! starting offset for the potential
  mpi_enreg%nscatterarr(jproc, 3) = density_start + potential_shift - 1
! GGA XC shift between density and potential
  mpi_enreg%nscatterarr(jproc, 4) = potential_shift
 end do

 mpi_enreg%ngatherarr(:, 1) = dtset%wvl_internal%dpSize(1) * &
& dtset%wvl_internal%dpSize(2) * mpi_enreg%nscatterarr(:, 2)
 mpi_enreg%ngatherarr(:, 2) = dtset%wvl_internal%dpSize(1) * &
& dtset%wvl_internal%dpSize(2) * mpi_enreg%nscatterarr(:, 3)

!Now ngfft will use the density definition (since the potential size
!is always smaller than the density one).
 dtset%ngfft(1) = dtset%wvl_internal%dpSize(1)
 dtset%ngfft(2) = dtset%wvl_internal%dpSize(2)
 dtset%ngfft(3) = mpi_enreg%nscatterarr(mpi_enreg%me, 1)

 dtset%nfft = product(dtset%ngfft(1:3))
!Set up fft array dimensions ngfft(4,5,6) to avoid cache conflicts
!Code paste from getng()
 dtset%ngfft(4) = 2 * (dtset%ngfft(1) / 2) + 1
 dtset%ngfft(5) = 2 * (dtset%ngfft(2) / 2) + 1
 dtset%ngfft(6) = dtset%ngfft(3)
 if (mpi_enreg%nproc == 0) then
  dtset%ngfft(9)  = 0    ! paral_fft
  dtset%ngfft(10) = 1    ! nproc_fft
  dtset%ngfft(11) = 0    ! me_fft
  dtset%ngfft(12) = 0    ! n2proc
  dtset%ngfft(13) = 0    ! n3proc
 else
  dtset%ngfft(9)  = 1    ! paral_fft
  dtset%ngfft(10) = mpi_enreg%nproc_fft
  dtset%ngfft(11) = mpi_enreg%me_fft
  dtset%ngfft(12) = dtset%ngfft(2)
  dtset%ngfft(13) = dtset%ngfft(3)
 end if
 write(message, '(a,3I12)' ) &
& '  | ngfft(1:3) is now:    ', dtset%ngfft(1:3)
 call wrtout(6,message,'COLL')
 write(message, '(a,3I12)' ) &
& '  | ngfft(4:6) is now:    ', dtset%ngfft(4:6)
 call wrtout(6,message,'COLL')

!Set mgfft
 dtset%mgfft= max(dtset%ngfft(1), dtset%ngfft(2), dtset%ngfft(3))  
 
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_setngfft : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_setngfft
!!***
