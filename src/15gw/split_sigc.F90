!{\src2tex{textfont=tt}}
!!****f* ABINIT/split_sigc
!! NAME
!! split_sigc
!!
!! FUNCTION
!!  Split the correlation part in Coulomb hole and screened Exchange contributions
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!!
!! NOTES
!! MG has created this routine while trying to clean a bit csigme.F90.
!! I just did a cut&paste of the piece of code present in pristine version 5.6
!! Use at your own risk since this part is not supported anymore and no 
!! automatic test has been provided.
!!
!! PARENTS
!!  
!!
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine split_sigc(sp,sr,npwc1,npwc2,jb,isppol,io,ioe0j,theta_mu_minus_e0i,&
& rhotwg,omegame0i,otq,botsq,sigccoh,sigcsex)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!io has no meaning in the present status of the code ??
!scalars
 integer,intent(in) :: io,ioe0j,isppol,jb,npwc1,npwc2
 real(dp),intent(in) :: theta_mu_minus_e0i
 type(sigma_parameters),intent(in) :: sp
 type(sigma_results),intent(in) :: sr
!arrays
 real(dp),intent(in) :: omegame0i(sp%nomegasr+sp%nomegasrd),otq(sp%npwc,npwc2)
 complex(gwpc) :: botsq(sp%npwc,npwc1)
 complex(gwpc),intent(in) :: rhotwg(sp%npwx)
 complex(gwpc),intent(inout) :: sigccoh(sp%nbnds,sp%nsppol)
 complex(gwpc),intent(inout) :: sigcsex(sp%nbnds,sp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: ff,ig,igp,twofm1
 real(dp) :: den_coh,den_sex,kincontrib,otw,reomegame0i,twofm1_zcut
 complex(gwpc) :: ccoh,csex,ct,dct,dsigc,idelta,num,omegame0i2_ac,omegame0i_ac
 complex(gwpc) :: rhotwgdp_igp,sigc,sigcohme,sigxme,twofm1_idelta,zz
 logical :: cohsex
 character(len=500) :: msg
!arrays
 complex(gwpc) :: rhotwgdpcc(sp%npwx)

! *************************************************************************

 ! Decomposition Sigma_c into Coulomb-hole (coh) and screened-exchange (sex) (GMR)
 ! MG case nsppol==2 is not tested Moreover there is no automatic test
 cohsex=.TRUE.

 ff=0
 if (theta_mu_minus_e0i>0.5) ff=1
 twofm1=2*ff-1 
 twofm1_zcut=twofm1*sp%zcut 
 idelta=CMPLX(0,sp%zcut)
 twofm1_idelta=twofm1*idelta

 !Now, introduce rhotwgdpcc, for speed reasons
 rhotwgdpcc(:)=conjg(rhotwg(:))
 reomegame0i=omegame0i(ioe0j)
 if ((cohsex) .and. (io==sr%nomega+ioe0j)) then !MG this won t work for sure but it was present also in v5.6
  ccoh=0 ; csex=0
  do igp =1,sp%npwc
   rhotwgdp_igp= rhotwg(igp)
   do ig =1,sp%npwc
    if (sp%ppmodel==3 .and. ig/=igp) cycle
    otw= dble(otq(ig,igp)) ! in principle otw -> otw - ieta
    num= rhotwgdpcc(ig)*botsq(ig,igp)*rhotwgdp_igp
    den_coh= reomegame0i-otw
    if (den_coh**2>sp%zcut**2) then
     ccoh= ccoh + num/(den_coh*otw)
    else ! if den is small
     ccoh= ccoh + num*cmplx(den_coh,twofm1_zcut)/ &
&     ((den_coh**2+twofm1_zcut**2)*otw)
    end if
    den_sex= reomegame0i**2-otw**2
    if (den_sex**2>sp%zcut**2) then
     csex= csex - 2*ff*num/den_sex
    else ! if den is small
    csex= csex -2*ff*num*cmplx(den_sex,twofm1_zcut)/ &
&    ((den_sex**2+twofm1_zcut**2)*otw)
    end if
   end do !ig
  end do !igp
  ! Must be multiplied later by 1/(ucvol*nkbz). 
  ! Note that 4pi is contained in vc_sqrt
  sigccoh(jb,isppol)= sigccoh(jb,isppol) + ccoh/two
  sigcsex(jb,isppol)= sigcsex(jb,isppol) + csex/two
 end if ! cohsex
 !
 ! Calculation of the kinetic contribution to the bandgap energy. (YMN on 07/04/04)
 ! if ((kinetic).and.(io == sr%nomega+ioe0j)) then
 ! print *,' omega = ',sr%omegasrd(jb,jkibz,io-sp%nomegasr)
 ! do igp = 1, sp%npwc
 ! do ig = 1, sp%npwc
 ! dct = zero
 ! do ig1 = 1, sp%npwc
 ! call fkin(ff,kincontrib,reomegame0i,otq(ig,ig1),otq(ig1,igp),sp%zcut)
 ! dct = dct+botsq(ig,ig1)*botsq(ig1,igp)*kincontrib
 ! end do
 ! ct = ct-rhotwgdpcc(ig)*dct*rhotwg(igp)
 ! end do
 ! end do
 ! sigckin(jb) = sigckin(jb)-ct*half  !Must be multiplied later by 4*pi/(ucvol*nkbz).
 ! end if ! kinetic

end subroutine split_sigc
!!***
