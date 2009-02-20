!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_sig_noppm
!! NAME
!! calc_sig_noppm
!!
!! FUNCTION
!! Calculate contributions to self-energy operator without a plasmon pole model.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (FB, GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nomega= total number of frequencies where evaluate $\Sigma_c$ matrix elements
!!  nomegae=number of frequencies where $\epsilon^{-1}$ has been evaluated 
!!  nomegaei= number of imaginary frequencies for $\epsilon^{-1}$ (non zero)
!!  nomegaer= number of real frequencies for $\epsilon^{-1}$
!!  npwc=number of G vectors for correlation part.
!!  npwx=number of G vectors in rhotwgp for each spinorial component.
!!  nspinor=Number of spinorial components.
!!  theta_mu_minus_e0i=1 if e0i is occupied, 0 otherwise. Fractional occupancy in case of metals. 
!!  omegame0i(nomega)= contains $\omega-\epsilon_{k-q,b1,\sigma}$
!!  epsm1q(npwc,npwc,nomegae)=symmetrized inverse dielectric matrix (exchange part is subtracted).
!!  omega(nomegae)=frequencies for $\epsilon^{-1}$
!!  rhotwgp(npwx*nspinor)=oscillator matrix elements: $<k-q,b1,\sigma|e^{-i(q+G)r} |k,b2,\sigma>*vc_sqrt$
!!
!! OUTPUT
!! ket(npwc,nomega)=Contain \Sigma_c(\omega)|\phi> in reciprocal space. 
!!
!! NOTES
!!
!! PARENTS
!!      csigme
!!
!! CHILDREN
!!      cgemv,spline,splint
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_sig_noppm(npwc,npwx,nspinor,nomega,nomegae,nomegaer,nomegaei,rhotwgp,&
& omega,epsm1q,omegame0i,theta_mu_minus_e0i,ket)

 use defs_basis
 use m_gwdefs, only : czero_gw, cone_gw


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib00numeric
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,nomegae,nomegaei,nomegaer,npwc,npwx,nspinor
 real(dp),intent(in) :: theta_mu_minus_e0i
!arrays
 real(dp),intent(in) :: omegame0i(nomega)
 complex(gwpc),intent(in) :: epsm1q(npwc,npwc,nomegae),omega(nomegae)
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(npwc*nspinor,nomega)

!Local variables-------------------------------
#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'CGEMV' :: cgemv
#endif
!scalars
 integer :: ig,io,ios,ispinor,spadc,spadx
 real(dp) :: omg1,omg2,rt_imag,rt_real
 complex(gwpc) :: ct,domegaleft,domegaright,fact
!arrays
 real(dp) :: omegame0i_tmp(nomega),rtmp(nomegaer),tmp_x(2),tmp_y(2)
 real(dp) :: work(nomegaer)
 complex(gwpc) :: epsrho(npwc,nomegae),epsrho_imag(npwc,nomegaei+1)
 complex(gwpc) :: omega_imag(nomegaei+1)

!*************************************************************************

!DEBUG
!write(6,*)' calc_sig_noppm : enter '
!ENDDEBUG

 ! Avoid divergences
 omegame0i_tmp(:)=omegame0i(:)
 do ios=1,nomega
  if (ABS(omegame0i_tmp(ios))<1.d-6) omegame0i_tmp(ios)=1.d-6
 end do

 do ispinor=1,nspinor
  spadx=(ispinor-1)*npwx
  spadc=(ispinor-1)*npwc
  !
  ! Calculate $\sum_{Gp} (\epsilon^{-1}_{G Gp}(\omega)-\delta_{G Gp}\rhotwgp(Gp)$
  do io=1,nomegae
#if defined HAVE_GW_DPC
   call ZGEMV('N',npwc,npwc,cone_gw,epsm1q(:,:,io),npwc,rhotwgp(1+spadx),1,czero_gw,epsrho(:,io),1)
#else
   call CGEMV('N',npwc,npwc,cone_gw,epsm1q(:,:,io),npwc,rhotwgp(1+spadx),1,czero_gw,epsrho(:,io),1)
#endif
  end do

  epsrho_imag(:,1)=epsrho(:,1)
  epsrho_imag(:,2:nomegaei+1)=epsrho(:,nomegaer+1:nomegae)
  omega_imag(1)=omega(1)
  omega_imag(2:nomegaei+1)=omega(nomegaer+1:nomegae)
  !
  ! === Perform integration along the imaginary axis ===
  do io=1,nomegaei+1
   if (io==1) then
    domegaleft  = omega_imag(io)
    domegaright =(omega_imag(io+1)-omega_imag(io  ))*half
   else if (io==nomegaei+1) then
    domegaleft  =(omega_imag(io  )-omega_imag(io-1))*half
    domegaright =(omega_imag(io  )-omega_imag(io-1))*half
   else
    domegaleft  =(omega_imag(io  )-omega_imag(io-1))*half
    domegaright =(omega_imag(io+1)-omega_imag(io  ))*half
   end if

   do ios=1,nomega
    omg2 = -AIMAG(omega_imag(io)+domegaright)/REAL(omegame0i_tmp(ios))
    omg1 = -AIMAG(omega_imag(io)-domegaleft )/REAL(omegame0i_tmp(ios))
    fact = ATAN(omg2)-ATAN(omg1)
    ket(spadc+1:spadc+npwc,ios)=ket(spadc+1:spadc+npwc,ios)+epsrho_imag(:,io)*fact
   end do
  end do !io

  ket(spadc+1:spadc+npwc,:)=ket(spadc+1:spadc+npwc,:)/pi

  ! === Add contribution coming from poles ===
  do ios=1,nomega
   do ig=1,npwc
    ! === Interpolate real and imaginary part of epsrho at |omegame0i_tmp| ===
    call spline(DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),nomegaer,zero,zero,rtmp,work)
    tmp_x(1) = ABS(omegame0i_tmp(ios))
    call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(epsrho(ig,1:nomegaer)),rtmp,1,tmp_x,tmp_y)
    rt_real = tmp_y(1)

    call spline(DBLE(omega(1:nomegaer)),DBLE(AIMAG(epsrho(ig,1:nomegaer))),nomegaer,zero,zero,rtmp,work)
    tmp_x(1) = ABS(omegame0i_tmp(ios))
    call splint(nomegaer,DBLE(omega(1:nomegaer)),DBLE(AIMAG(epsrho(ig,1:nomegaer))),rtmp,1,tmp_x,tmp_y)
    rt_imag = tmp_y(1)

    ct=CMPLX(rt_real,rt_imag)

    if (omegame0i_tmp(ios)>tol12) then
     ket(spadc+ig,ios)=ket(spadc+ig,ios)+ct*(one-theta_mu_minus_e0i)
    end if
    if (omegame0i_tmp(ios)<-tol12) then
     ket(spadc+ig,ios)=ket(spadc+ig,ios)-ct*theta_mu_minus_e0i
    end if
   end do !ig
  end do !ios

 end do !ispinor

!DEBUG
!write(*,*) npwc,ket(:)
!stop
!ENDDEBUG

end subroutine calc_sig_noppm
!!***
