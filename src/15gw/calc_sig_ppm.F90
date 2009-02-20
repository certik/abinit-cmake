!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_sig_ppm
!!
!! NAME
!! calc_sig_ppm
!!
!! FUNCTION
!!  Calculate the contribution to self-energy operator using a plasmon-pole model
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (FB, GMR, VO, LR, RWG, RShaltaf, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nomega=number of frequencies
!!  nspinor=Number of spinorial components.
!!  npwc= number of G vectors in the plasmon pole 
!!  npwx=number of G vectors in rhotwgp, i.e no. of G"s for the exchange part
!!  theta_mu_minus_e0i= $\theta(\mu-\epsilon_{k-q,b1,s}), defines if the state is occupied or not 
!!  zcut=small imaginary part to avoid the divergence. (see related input variable)
!!  omegame0i(nomega)=frequencies where evaluate \Sigma_c ($\omega$ - $\epsilon_i$ 
!!  otq(npwc,dm2_otq)=plasmon pole parameters for this q-point
!!  PPm<PPmodel_type>=structure gathering info on the Plasmon-pole technique.
!!     %model=type plasmon pole model
!!     %dm2_botsq= 1 if model==3, =npwc if model== 4, 1 for all the other cases
!!     %dm2_otq= 1 if model==3, =1    if model== 4, 1 for all the other cases
!!     %dm_eig=npwc if model=3, 0 otherwise
!!  botsq(npwc,dm2_botsq)=plasmon pole parameters for this q-point
!!  eig(dm_eig,dm_eig)=the eigvectors of the symmetrized inverse dielectric matrix for this q point
!!   (first index for G, second index for bands)
!!  rhotwgp(npwx)=oscillator matrix elements divided by |q+G| i.e 
!!    $\frac{\langle b1 k-q s | e^{-i(q+G)r | b2 k s \rangle}{|q+G|}$ 
!!
!! OUTPUT
!!  sigcme(nomega) (to be described), only relevant if ppm3 or ppm4
!!
!!  ket(npwc,nomega): 
!!
!!  In case of model==1,2 it contains
!!
!!   ket(G,omega) = Sum_G2       conjg(rhotw(G)) * Omega(G,G2) * rhotw(G2)
!!                          ---------------------------------------------------
!!                            2 omegatw(G,G2) (omega-E_i + omegatw(G,G2)(2f-1))
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

subroutine calc_sig_ppm(PPm,nspinor,npwc,nomega,rhotwgp,botsq,otq,&
& omegame0i,zcut,theta_mu_minus_e0i,eig,npwx,ket,sigcme)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,npwx,nspinor
 real(dp),intent(in) :: theta_mu_minus_e0i,zcut
 type(PPmodel_type),intent(in) :: PPm
!arrays
 real(dp),intent(in) :: omegame0i(nomega),otq(npwc,PPm%dm2_otq)
 complex(gwpc),intent(in) :: botsq(npwc,PPm%dm2_botsq),eig(PPm%dm_eig,PPm%dm_eig)
 complex(gwpc),intent(in) :: rhotwgp(npwx*nspinor)
 complex(gwpc),intent(inout) :: ket(npwc*nspinor,nomega)
 complex(gwpc),intent(out) :: sigcme(nomega)

!Local variables-------------------------------
!scalars
 integer :: ig,igp,ii,ios,ispinor,spadc,spadx
 real(dp) :: den,ff,inv_den,omegame0i_io,otw,twofm1,twofm1_zcut
 complex(gwpc) :: ct,num,numf,rhotwgdp_igp
 logical :: fully_occupied,totally_empty
 character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: rhotwgdpcc(:)

!*************************************************************************

!DEBUG
!write(6,*)' calc_sig_ppm : enter '
!write(*,*)npwc,zcut
!write(*,*)rhotwgp(1),botsq(1,1),otq(1,1)
!ENDDEBUG

 SELECT CASE (PPm%model)

 CASE (1,2)

  fully_occupied =(ABS(theta_mu_minus_e0i-one)<0.001)
  totally_empty  =(ABS(theta_mu_minus_e0i    )<0.001)

  do ispinor=1,nspinor

   spadx=(ispinor-1)*npwx
   spadc=(ispinor-1)*npwc

   if (.not.totally_empty) then  
    ! \Bomega^2_{G1G2}/\omegat_{G1G2} M_{G1,G2}. \theta(\mu-e_s) / (\omega+\omegat_{G1G2}-e_s-i\delta)
    twofm1_zcut=zcut
    do ios=1,nomega
     omegame0i_io=omegame0i(ios)
     do igp=1,npwc
      rhotwgdp_igp=rhotwgp(spadx+igp)
      do ig=1,npwc
       otw=otq(ig,igp) !in principle otw -> otw - ieta
       num = botsq(ig,igp)*rhotwgdp_igp
       den = omegame0i_io+otw
       if (den**2>zcut**2)then
        ket(spadc+ig,ios)=ket(spadc+ig,ios) + num/(den*otw)*theta_mu_minus_e0i
       else
        ket(spadc+ig,ios)=ket(spadc+ig,ios) + num*CMPLX(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)*theta_mu_minus_e0i
       end if
      end do !ig
     end do !igp
    end do !ios
   end if !not totally empty

   if (.not.(fully_occupied)) then 
    ! \Bomega^2_{G1G2}/\omegat_{G1G2} M_{G1,G2}. \theta(e_s-\mu) / (\omega-\omegat_{G1G2}-e_s+i\delta)
    twofm1_zcut=-zcut
    do ios=1,nomega
     omegame0i_io=omegame0i(ios)
     do igp=1,npwc
      rhotwgdp_igp=rhotwgp(spadx+igp)
      do ig=1,npwc
       otw=otq(ig,igp) !in principle otw -> otw - ieta
       num = botsq(ig,igp)*rhotwgdp_igp
       den=omegame0i_io-otw
       if (den**2>zcut**2) then
        ket(spadc+ig,ios)=ket(spadc+ig,ios) + num/(den*otw)*(one-theta_mu_minus_e0i)
       else
        ket(spadc+ig,ios)=ket(spadc+ig,ios) + &
&        num*CMPLX(den,twofm1_zcut)/((den**2+twofm1_zcut**2)*otw)*(one-theta_mu_minus_e0i)
       end if
      end do !ig
     end do !igp
    end do !ios
   end if !not fully occupied

  end do !ispinor

  ket(:,:)=ket(:,:)*half

 CASE (3,4)
  if (nspinor==2) call assert(.FALSE.,'nspinor==2 not allowed',&
&  __FILE__,__LINE__)
  !
  ! * rho-twiddle(G) is formed, introduce rhotwgdpcc, for speed reason
  allocate(rhotwgdpcc(npwx))

  ff=theta_mu_minus_e0i      ! occupation number f (include poles if ...)
  twofm1=two*ff-one          ! 2f-1
  twofm1_zcut=twofm1*zcut
  rhotwgdpcc(:)=CONJG(rhotwgp(:))

  do ios=1,nomega
   omegame0i_io=omegame0i(ios)
   ct=(0.0_gwp,0.0_gwp)
   do ii=1,npwc   ! DM bands
    num=(0.0_gwp,0.0_gwp)

    SELECT CASE (PPm%model)
    CASE (3) 
     ! * Calculate \beta (eq. 106 pag 47) 
     do ig=1,npwc
      num=num+rhotwgdpcc(ig)*eig(ig,ii)
     end do
     numf=num*CONJG(num) !MG this means that we cannot do SCGW 
     numf=numf*botsq(ii,1)
    CASE (4)
     do ig=1,npwc
      num=num+rhotwgdpcc(ig)*botsq(ig,ii)
     end do
     numf=num*CONJG(num) !MG this means that we cannot do SCGW 
    CASE DEFAULT 
     call assert(.FALSE.,'wrong PPm%model',__FILE__,__LINE__)
    END SELECT

    !numf=num*CONJG(num) !MG this means that we cannot do SCGW 
    !if (PPm%model==3) numf=numf*botsq(ii,1)

    otw=DBLE(otq(ii,1)) ! in principle otw -> otw - ieta
    den=omegame0i_io+otw*twofm1

    if (den**2>zcut**2) then
     inv_den=one/den
     ct=ct+numf*inv_den
    else 
     inv_den=one/((den**2+twofm1_zcut**2))
     ct=ct+numf*CMPLX(den,twofm1_zcut)*inv_den
    end if

   end do !ii DM bands
   sigcme(ios)=ct*half
  end do !ios
  deallocate(rhotwgdpcc)

 CASE DEFAULT
  write(msg,'(4a,i3)')ch10,&
&  ' calc_sig_ppm : BUG ',ch10,&
&  '  Wrong value for PPm%model = ',PPm%model
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT

end subroutine calc_sig_ppm
!!***
