!{\src2tex{textfont=tt}}
!!****f* ABINIT/fft_onewfn
!! NAME
!! fft_onewfn
!!
!! FUNCTION
!! Calculate ONE wavefunction in real space using FFT
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf, FB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nspinor=number of spinorial components
!! igfft(npwwfn)=index of each plane wave in FFT grid
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nkibz=number of k points
!! npwwfn=number of plane waves
!! nsppol=number of independent spin polarizations 
!! tim_fourdp=4 if called from within screening ; =5 if called from within sigma
!! wfg(npwwfn,my_minb:my_maxb,nkibz,nsppol)=wavefunctions in reciprocal space treated by this processor.
!! my_minb,my_maxb = min and max band treated by this processor
!! MPI_enreg= datatype containing information on parallelism to be passed to fourdp
!!
!! OUTPUT
!!  wfr(ngfft(1)*ngfft(2)*ngfft(3)*nspinor)=wavefunctions in real space.
!!
!! PARENTS
!!      calc_density,calc_vHxc_braket,debug_tools,wf_info
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fft_onewfn(paral_kgb,nspinor,npwwfn,nfftot,wfg,wfr,igfft,ngfft,tim_fourdp,MPI_enreg)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: paral_kgb,npwwfn,nfftot,tim_fourdp,nspinor
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: igfft(npwwfn),ngfft(18)
 complex(gwpc),intent(in) :: wfg(npwwfn*nspinor)
 complex(gwpc),intent(out) :: wfr(ngfft(1)*ngfft(2)*ngfft(3)*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ispinor,ig,master,me,spaceComm,rspad,gspad
!arrays
 real(dp),allocatable :: wfg_dp(:,:),wfr_dp(:,:)

! *************************************************************************

 !call xcomm_init  (MPI_enreg,spaceComm) 
 !call xme_init    (MPI_enreg,me)         
 !call xmaster_init(MPI_enreg,master)  
 !MPI_enreg%me_fft=0 ; MPI_enreg%nproc_fft=1

 allocate(wfg_dp(2,nfftot),wfr_dp(2,nfftot))

 do ispinor=1,nspinor
  gspad=(ispinor-1)*npwwfn
  rspad=(ispinor-1)*nfftot
  !
  ! === Fill FFT array from PW array ===
  wfg_dp(:,:)=zero
  do ig=1,npwwfn
   wfg_dp(1,igfft(ig))=REAL (wfg(ig+gspad))
   wfg_dp(2,igfft(ig))=AIMAG(wfg(ig+gspad))
  end do
  !
  ! === Take FFT to give wfn in real space ===
  ! here wfr_dp doesnt has same shape as fofr
  call fourdp(2,wfg_dp(:,:),wfr_dp(:,:),+1,MPI_enreg,nfftot,ngfft,paral_kgb,tim_fourdp)
  
  wfr(1+rspad:nfftot+rspad)=CMPLX(wfr_dp(1,:),wfr_dp(2,:),kind=gwpc)
 end do

 deallocate(wfg_dp,wfr_dp)

end subroutine fft_onewfn
!!***
