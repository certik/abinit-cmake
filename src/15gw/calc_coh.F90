!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_coh
!! NAME
!! calc_coh
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (FB,MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! iqibz=index of the irreducible q-point in the array qibz, point which is 
!!  related by a symmetry operation to the point q summed over (see csigme). 
!!  This index is also used to treat the integrable coulombian singularity at q=0
!! jb,kb=left and righ band indeces definining the left and right states where the 
!!  partial contribution to the matrix element of $\Sigma_{COH}$ is evaluated
!! ngfft_gw(18)=contain all needed information about 3D FFT for GW wavefuntions,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nsig_ab=Number of components in the self-energy operator (1 for collinear magnetism) 
!! npwc=number of plane waves in $\tilde epsilon^{-1}$
!! npwx=number of G vectors in the arrays gvec and vc_sqrt
!! nspinor=Number of spinorial components.
!! nfftot=number of points in real space
!! i_sz=contribution arising from the integrable coulombian singularity at q==0 
!! (see csigme for the method used), note that in case of 3-D systems the factor 
!! 4pi in the coulombian potential is included in the definition of i_sz 
!! gvec(3,npwx)=G vectors in reduced coordinates 
!! vc_sqrt(npwx)= square root of the coulombian matrix elements for this q-point
!! epsm1q_o(npwc,npwc)= contains $\tilde epsilon^{-1}(q,w=0) - \delta_{G Gp}$ for 
!!  the particular q-point considered in the sum
!! wfg2_jk(nfftot)= Fourier Transform of $\u_{jb k}^*(r) u_{kb k}$ 
!! spinrot_k(4)=components of the spinor rotation matrix.
!!
!! OUTPUT
!! sigcohme=partial contribution to the matrix element of 
!!     $<jb k \sigma|\Sigma_{COH} | kb k \sigma>$ 
!!  coming from this single q-point
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      csigme
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_coh(paral_kgb,nspinor,nsig_ab,nfftot,ngfft_gw,tim_fourdp,MPI_enreg,ktabi_k,ktabr_k,spinrot_k,&
& wfr_jb,wfr_kb,npwx,npwc,gvec,epsm1q_o,vc_sqrt,i_sz,iqibz,same_band,sigcohme)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : czero_gw, cone_gw


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => calc_coh
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqibz,ktabi_k,nfftot,npwc,npwx,nsig_ab,nspinor,paral_kgb
 integer,intent(in) :: tim_fourdp
 real(dp),intent(in) :: i_sz
 logical,intent(in) :: same_band
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: gvec(3,npwx),ktabr_k(nfftot),ngfft_gw(18)
 real(dp),intent(in) :: spinrot_k(4)
 complex(gwpc),intent(in) :: epsm1q_o(npwc,npwc),vc_sqrt(npwx),wfr_jb(nfftot*nspinor)
 complex(gwpc),intent(in) :: wfr_kb(nfftot*nspinor)
 complex(gwpc),intent(out) :: sigcohme(nsig_ab)

!Local variables-------------------------------
!scalars
 integer :: ig,ig4,ig4x,ig4y,ig4z,igp,ispinor,map2sphere,ngfft1,ngfft2,ngfft3
 integer :: spad
!arrays
 integer,allocatable :: igfftg0_dummy(:)
 complex(gwpc),allocatable :: epsm1q_qpg2(:,:),test(:),wfg2_jk(:)

! *************************************************************************

!DEBUG
!write(6,*)' calc_coh : enter '
!ENDDEBUG

 allocate(epsm1q_qpg2(npwc,npwc))
 allocate(wfg2_jk(nfftot*nspinor))
 !
 ! === Calculate the product of two wfr and Fourier-transform it in recip. space ===
 !call calc_wfwfg(mpi_enreg,paral_kgb,tim_fourdp,ktabr_k,ktabi_k,nfftot,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)
 !allocate(test(nfftot*nspinor)) ; test=wfg2_jk

 ! TODO for PAW:
 ! 1) add onsite contributions to calc_wfwfg on the FFT box, see paw_rho_tw_g
 ! non-symmporphic phases not needed as wfs are evaluated at the same k-point

 map2sphere=0
 allocate(igfftg0_dummy(map2sphere))
 call rho_tw_g(paral_kgb,nspinor,nfftot,nfftot,ngfft_gw,map2sphere,igfftg0_dummy,&
& wfr_jb,ktabi_k,ktabr_k,cone_gw,spinrot_k,& !ph_mkt  ,& 
& wfr_kb,ktabi_k,ktabr_k,cone_gw,spinrot_k,& !ph_mkgwt,& 
& nspinor,wfg2_jk,tim_fourdp,MPI_enreg)
 deallocate(igfftg0_dummy)

 !if (allocated(test)) then
 ! write(std_out,*)' == COH test == ',MAXVAL(ABS(test-wfg2_jk))
 ! deallocate(test)
 !end if
 !
 ! === Set up of \epsilon^{-1} vc_sqrt(q,G) vc_sqrt(q,Gp) ===
 ! * vc_sqrt contains sqrt(vc) i.e 4\pi/|q+G| in 3-D systems
 do igp=1,npwc
  do ig=1,npwc
   epsm1q_qpg2(ig,igp) = epsm1q_o(ig,igp)*vc_sqrt(ig)*vc_sqrt(igp)
  end do
 end do
 !
 ! === Treat the case q --> 0 adequately ===
 ! TODO Better treatment of wings     
 !      check cutoff in the coulombian interaction.
 if (iqibz==1) then
  if (same_band) then
   epsm1q_qpg2(1,1)=epsm1q_o(1,1)*i_sz
   epsm1q_qpg2(2:npwc,1)=czero_gw
   epsm1q_qpg2(1,2:npwc)=czero_gw
  else
   epsm1q_qpg2(:,1)=czero_gw
   epsm1q_qpg2(1,:)=czero_gw
  end if
 end if
 !
 ! === Partial contribution to the matrix element of Sigma_c ===
 ! * For nspinor==2, the closure relation reads: 
 !  $\sum_s \psi_a^*(1)\psi_b(2) = \delta_{ab} \delta(1-2)$
 !  where a,b are the spinor components. As a consequence, Sigma_{COH} is always 
 !  diagonal in spin-space and only diagonal matrix elements have to be calculated.
 ! MG  TODO wfg2_jk should be calculated on an augmented FFT box to avoid spurious wrapping of G1-G2.
 !
 ngfft1=ngfft_gw(1) 
 ngfft2=ngfft_gw(2) 
 ngfft3=ngfft_gw(3)
 sigcohme(:)=czero_gw

 do ispinor=1,nspinor
  spad=(ispinor-1)*nfftot

  do igp=1,npwc
   do ig=1,npwc
    ig4x=MODULO(gvec(1,igp)-gvec(1,ig),ngfft1)
    ig4y=MODULO(gvec(2,igp)-gvec(2,ig),ngfft2)
    ig4z=MODULO(gvec(3,igp)-gvec(3,ig),ngfft3)
    ig4= 1+ig4x+ig4y*ngfft1+ig4z*ngfft1*ngfft2
    sigcohme(ispinor) = sigcohme(ispinor) + half*wfg2_jk(spad+ig4)*epsm1q_qpg2(ig,igp)
   end do
  end do

 end do !ispinor

 deallocate(epsm1q_qpg2)
 deallocate(wfg2_jk)

end subroutine calc_coh
!!***


!!****if* ABINIT/calc_wfwfg
!! NAME
!! calc_wfwfg
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_wfwfg(MPI_enreg,paral_kgb,tim_fourdp,ktabr_k,ktabi_k,nfftot,ngfft_gw,wfr_jb,wfr_kb,wfg2_jk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ktabi_k,nfftot,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: ktabr_k(nfftot),ngfft_gw(18)
 complex(gwpc),intent(in) :: wfr_jb(nfftot),wfr_kb(nfftot)
 complex(gwpc),intent(out) :: wfg2_jk(nfftot)

!Local variables-------------------------------
!arrays
 real(dp),allocatable :: wfg2_tmp(:,:),wfr2_tmp(:,:)

! *************************************************************************

 allocate(wfg2_tmp(2,nfftot),wfr2_tmp(2,nfftot))
 ! There is no need to take into account phases arising from non-symmorphic
 ! operations since the wavefunctions are evaluated at the same k-point.

 wfr2_tmp(1,:)=   real(wfr_jb(ktabr_k(:))) *  real(wfr_kb(ktabr_k(:)))&
&               +aimag(wfr_jb(ktabr_k(:))) * aimag(wfr_kb(ktabr_k(:)))

 wfr2_tmp(2,:)=  real(wfr_jb(ktabr_k(:))) * aimag(wfr_kb(ktabr_k(:)))&
&              -aimag(wfr_jb(ktabr_k(:))) *  real(wfr_kb(ktabr_k(:)))

 ! Conjugate the product if time-reversal is used to reconstruct this k-point
 if (ktabi_k==-1) wfr2_tmp(2,:)=-wfr2_tmp(2,:)

 call fourdp(2,wfg2_tmp,wfr2_tmp,-1,MPI_enreg,nfftot,ngfft_gw,paral_kgb,tim_fourdp)

 wfg2_jk(:)= wfg2_tmp(1,:)+(0.,1.)*wfg2_tmp(2,:)

 deallocate(wfg2_tmp,wfr2_tmp)

end subroutine calc_wfwfg
!!***
