!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_rho_tw_g
!! NAME
!! paw_rho_tw_g
!!
!! FUNCTION
!! Evaluate localized contribution to oscillator matrix elements in case of PAW calculations. Namely:
!!  sum_{a,i,j} <\tilde\psi_{k-q,b1}|\tilde p_i^a> <\tilde p_j^a|\tilde\psi_{k,b2}>*
!!   \[<\phi_i^a|e^{-i(q+G).r}|\phi_j^a> - <\tilde\phi_i^a|e^{-i(q+G).r}|\tilde\phi_j^a>\]
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! dim_rtwg=Define the size of the array rhotwg
!!   === for nspinor==1 ===
!!    dim_rtwg=1
!!   === for nspinor==2 ===
!!    dim_rtwg=2 if only <up|up>, <dwn|dwn> matrix elements are required
!!    dim_rtwg=4 for <up|up>, <dwn|dwn>, <up|dwn> and <dwn|up>.
!!  nspinor=number of spinorial components.
!!  npw=number of plane waves for oscillator matrix elements
!!  natom=number of atoms
!!  lmnmax=Max number of (l,m,n) channels 
!!  Cprj_k1(natom),Cprj_k2(natom) <type(cprj_type)>=projected input wave functions <Proj_i|Cnk> with all NL projectors
!!   k1 corresponds to k-q, k2 to k.
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  rhotwg(npw*dim_rtwg)=Updated oscillator matrix elements containing on-site PAW contributions
!!
!! NOTES
!!
!! TODO 
!! Add option map2sphere. paw_rhox has to be changed.
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

subroutine paw_rho_tw_g(npw,dim_rtwg,nspinor,natom,lmnmax,dimlmn,Cprj_k1,Cprj_k2,paw_rhox,rhotwg)
    
 use defs_basis
 use defs_datatypes

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,npw,lmnmax,nspinor,dim_rtwg
!arrays
 integer,intent(in) :: dimlmn(natom)
 real(dp),intent(in) :: paw_rhox(2,npw,lmnmax*(lmnmax+1)/2,natom)
 complex(gwpc),intent(inout) :: rhotwg(npw*dim_rtwg)
 type(Cprj_type),intent(in) :: Cprj_k1(natom,nspinor),Cprj_k2(natom,nspinor)

!Local variables-------------------------------
!scalars
 integer :: ig,iat,nlmn,ilmn,jlmn,k0lmn,klmn,iab,isp1,isp2,spad
 real(dp) :: fij,re_p,im_p
!arrays
 integer,parameter :: spinor_idxs(2,4)=RESHAPE((/1,1,2,2,1,2,2,1/),(/2,4/))
 real(dp) :: tmp(2)
 
! *************************************************************************

 ! === Loop over the four spinorial combinations ===
 do iab=1,dim_rtwg
  isp1=spinor_idxs(1,iab)
  isp2=spinor_idxs(2,iab)
  spad=npw*(iab-1)

  do ig=1,npw
   tmp(:)=zero

   do iat=1,natom
    nlmn=dimlmn(iat)
    ! === Loop on [(jl,jm,jn),(il,im,in)] channels. packed storage mode ===
    do jlmn=1,nlmn 
     k0lmn=jlmn*(jlmn-1)/2 
     do ilmn=1,jlmn

      re_p=   Cprj_k1(iat,isp1)%cp(1,ilmn) * Cprj_k2(iat,isp2)%cp(1,jlmn) &
&            +Cprj_k1(iat,isp1)%cp(2,ilmn) * Cprj_k2(iat,isp2)%cp(2,jlmn) &
&            +Cprj_k1(iat,isp1)%cp(1,jlmn) * Cprj_k2(iat,isp2)%cp(1,ilmn) &
&            +Cprj_k1(iat,isp1)%cp(2,jlmn) * Cprj_k2(iat,isp2)%cp(2,ilmn) 

      im_p=   Cprj_k1(iat,isp1)%cp(1,ilmn) * Cprj_k2(iat,isp2)%cp(2,jlmn) &
&            -Cprj_k1(iat,isp1)%cp(2,ilmn) * Cprj_k2(iat,isp2)%cp(1,jlmn) &
&            +Cprj_k1(iat,isp1)%cp(1,jlmn) * Cprj_k2(iat,isp2)%cp(2,ilmn) &
&            -Cprj_k1(iat,isp1)%cp(2,jlmn) * Cprj_k2(iat,isp2)%cp(1,ilmn)

      klmn=k0lmn+ilmn ; fij=one ; if (jlmn==ilmn) fij=half
      tmp(1)=tmp(1)+ fij*(paw_rhox(1,ig,klmn,iat)*re_p - paw_rhox(2,ig,klmn,iat)*im_p) 
      tmp(2)=tmp(2)+ fij*(paw_rhox(1,ig,klmn,iat)*im_p + paw_rhox(2,ig,klmn,iat)*re_p)  

     end do !ilmn
    end do !jlmn
   end do !iat
   !
   ! === Update input date using the appropriate index ===
   rhotwg(ig+spad) = rhotwg(ig+spad) + CMPLX(tmp(1),tmp(2),kind=gwpc)
  end do !ig

 end do !dim_rtwg

end subroutine paw_rho_tw_g
!!***
