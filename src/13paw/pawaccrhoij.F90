!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawaccrhoij
!!
!! NAME
!! pawaccrhoij
!!
!! FUNCTION
!! Accumulate the PAW quantities rhoij (augmentation occupancies)
!! or their 1-st order change or their gradient vs r
!! Add the contribution of a given k-point and band
!! Remember: for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  cwaveprj(dimpaw,nspinor)= wave function at given n,k
!!                         projected with non-local projectors: cwaveprj=<p_i|Cnk>
!!  cwaveprj1(dimpaw,nspinor)= 1st-order wave function at n,k,q
!!                          projected with non-local projectors: cwaveprj1=<p_i|C1nk,q>
!!                          USED for RF only - can be set to cwaveprj in case of GS
!!  dimpaw=size of cwaveprj, cprjwave1, pawrhoij (can be 1 or natom)
!!  ipert=index of perturbation (RF only, i.e. option=2)
!!  isppol=index of current spin component
!!  natom=number of atoms in cell
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components
!!  nsppol=number of channels for spin-polarization (1 or 2)
!!  occ_k=occupation number for current band n,k
!!  option: choice of calculation:
!!          1: update rhoij (Ground-State)
!!          2: update 1-st order rhoij (Response Function) according to ipert
!!          3: update gradients of rhoij with respect to r
!!  wtk_k=weight assigned to current k-point
!!
!! SIDE EFFECTS
!!  pawrhoij(dimpaw) <type(pawrhoij_type)>= GS: paw rhoij occupancies and related data
!!                                          RF: 1-st order paw rhoij occupancies and related data
!!  On output, has been updated with the contribution of current n,k
!!    === option=1:
!!        pawrhoij(:)%rhoij_(lmn2_size,nspden)=      (non symetrized)
!!            Sum_{n,k} {occ(n,k)*conjugate[cprj_nk(ii)].cprj_nk(jj)}
!!    === option=2:
!!        pawrhoij(:)%rhoij_(lmn2_size,nspden)=      (non symetrized)
!!            Sum_{n,k} {occ(n,k)*(conjugate[cprj_nk(ii)].cprj1_nk,q(jj)
!!                                 conjugate[cprj_nk(jj)].cprj1_nk,q(ii)}
!!          + Sum_{n,k} {occ(n,k)*(conjugate[dcprj_nk(ii)/dlambda].cprj_nk(jj)
!!                                +conjugate[cprj_nk(ii)].dcprj_nk(jj)/dlambda)}
!!    === option=3:
!!        pawrhoij(:)%grhoij(lmn2_size,mu,nspden)=   (non symetrized)
!!            Sum_{n,k} {occ(n,k)*(conjugate[dcprj_nk(ii)/dr_mu].cprj_nk(jj)
!!                                +conjugate[cprj_nk(ii)].dcprj_nk(jj)/dr_mu)}
!!
!! PARENTS
!!      dyfnl3,energy,pawmkrhoij,pawmkrhoij3,vtowfk3
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj1,dimpaw,ipert,isppol,natom,nspden,&
&                       nspinor,nsppol,occ_k,option,pawrhoij,wtk_k)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,dimpaw,ipert,isppol,natom,nspden,nspinor,nsppol,option
 real(dp) :: occ_k,wtk_k
!arrays
 integer,intent(in) :: atindx1(natom)
 type(cprj_type),intent(in) :: cwaveprj(dimpaw,nspinor),cwaveprj1(dimpaw,nspinor)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dimpaw)

!Local variables ---------------------------------------
!scalars
 integer :: iatm,iatom,ilmn,iplex,j0lmn,jlmn,klmn,klmn_im,klmn_re,mu,ncpgr
 real(dp) :: grhoij_tmp,rhoij_tmp,ro11_im,ro11_re,ro12_im,ro12_re,ro21_im,ro21_re,ro22_im,ro22_re,weight
 character(len=500) :: message
!arrays
 real(dp) :: cpi0(2,nspinor),cpi1(2,nspinor),cpj0(2,nspinor),cpj1(2,nspinor)
 real(dp) :: dcpi0(2,nspinor,3),dcpj0(2,nspinor,3)

!************************************************************************

 ncpgr=0
 if (option==2.and.ipert<=natom) ncpgr=1
 if (option==3) ncpgr=3

!Tests
 if(option==2.and.ipert>natom) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawaccrhoij : ERROR -',ch10,&
&  '  ipert>natom not yet implemented !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(option==2.and.cwaveprj1(1,1)%ncpgr/=ncpgr) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawaccrhoij : BUG -',ch10,&
&  '  Error on cwaveprj1 factors derivatives !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(option==3.and.cwaveprj(1,1)%ncpgr/=ncpgr) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawaccrhoij : BUG -',ch10,&
&  '  Error on cwaveprj factors derivatives !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 weight=wtk_k*occ_k
 if (nspden==2.and.nsppol==1.and.nspinor==1) weight=half*weight

 if (option==1) then

! ==================================================================
! === OPTION 1: Accumulate (n,k) contribution to rhoij =============
! ==================================================================

  if (nspinor==1) then
   do iatm=1,natom
    iatom=atindx1(iatm)
    do jlmn=1,pawrhoij(iatom)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     cpj0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,jlmn)
     do ilmn=1,jlmn
      klmn=j0lmn+ilmn
      cpi0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,ilmn)
      rhoij_tmp=zero
      do iplex=1,cplex
       rhoij_tmp=rhoij_tmp+cpi0(iplex,1)*cpj0(iplex,1)
      end do
      pawrhoij(iatom)%rhoij_(klmn,isppol)=pawrhoij(iatom)%rhoij_(klmn,isppol)+weight*rhoij_tmp
     end do
    end do
   end do
  else ! nspinor=2
   do iatm=1,natom
    iatom=atindx1(iatm)
    do jlmn=1,pawrhoij(iatom)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     cpj0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,jlmn)
     cpj0(1:cplex,2)=cwaveprj(iatm,2)%cp(1:cplex,jlmn)
     do ilmn=1,jlmn
      klmn=j0lmn+ilmn
      cpi0(1:cplex,1)=cwaveprj(iatm,1)%cp(1:cplex,ilmn)
      cpi0(1:cplex,2)=cwaveprj(iatm,2)%cp(1:cplex,ilmn)
      ro11_re=cpi0(1,1)*cpj0(1,1);ro22_re=cpi0(1,2)*cpj0(1,2)
      pawrhoij(iatom)%rhoij_(klmn,1)=pawrhoij(iatom)%rhoij_(klmn,1)+weight*(ro11_re+ro22_re)
      if (nspden>1) then
       ro12_re=cpi0(1,1)*cpj0(1,2);ro21_re=cpi0(1,2)*cpj0(1,1)
       pawrhoij(iatom)%rhoij_(klmn,4)=pawrhoij(iatom)%rhoij_(klmn,4)+weight*(ro11_re-ro22_re)
       pawrhoij(iatom)%rhoij_(klmn,2)=pawrhoij(iatom)%rhoij_(klmn,2)+weight*(ro12_re+ro21_re)
      end if
      if (cplex==2) then
       ro11_re=cpi0(2,1)*cpj0(2,1);ro22_re=cpi0(2,2)*cpj0(2,2)
       pawrhoij(iatom)%rhoij_(klmn,1)=pawrhoij(iatom)%rhoij_(klmn,1)+weight*(ro11_re+ro22_re)
       if (nspden>1) then
        ro12_re=cpi0(2,1)*cpj0(2,2);ro21_re=cpi0(2,2)*cpj0(2,1)
        ro12_im=cpi0(1,1)*cpj0(2,2)-cpi0(2,1)*cpj0(1,2)
        ro21_im=cpi0(1,2)*cpj0(2,1)-cpi0(2,2)*cpj0(1,1)
        pawrhoij(iatom)%rhoij_(klmn,4)=pawrhoij(iatom)%rhoij_(klmn,4)+weight*(ro11_re-ro22_re)
        pawrhoij(iatom)%rhoij_(klmn,2)=pawrhoij(iatom)%rhoij_(klmn,2)+weight*(ro12_re+ro21_re)
        pawrhoij(iatom)%rhoij_(klmn,3)=pawrhoij(iatom)%rhoij_(klmn,3)+weight*(ro21_im-ro12_im)
       end if
      end if
     end do
    end do
   end do
  end if

 else if (option==2) then

! ==================================================================
! === OPTION 2: Accumulate (n,k) contribution to 1st-order rhoij ===
! ==================================================================

! Accumulate (n,k) contribution to rhoij1
! due to derivative of wave-function
  if (nspinor==1) then
   do iatm=1,dimpaw
    iatom=min(atindx1(iatm),dimpaw)
    do jlmn=1,pawrhoij(iatom)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     cpj0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,jlmn)
     cpj1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,jlmn)
     do ilmn=1,jlmn
      klmn=j0lmn+ilmn
      klmn_re=cplex*(klmn-1)+1
      cpi0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,ilmn)
      cpi1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,ilmn)
      pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol) &
&      +weight*(cpi0(1,1)*cpj1(1,1)+cpj0(1,1)*cpi1(1,1)+cpi0(2,1)*cpj1(2,1)+cpj0(2,1)*cpi1(2,1))
      if (cplex==2) then
       klmn_im=klmn_re+1
       pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol) &
&       +weight*(cpi0(1,1)*cpj1(2,1)-cpi0(2,1)*cpj1(1,1)+cpj0(1,1)*cpi1(2,1)-cpj0(2,1)*cpi1(1,1))
      end if
     end do
    end do
   end do
  else ! nspinor=2
   do iatm=1,dimpaw
    iatom=min(atindx1(iatm),dimpaw)
    do jlmn=1,pawrhoij(iatom)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     cpj0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,jlmn)
     cpj0(1:2,2)=cwaveprj (iatm,2)%cp(1:2,jlmn)
     cpj1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,jlmn)
     cpj1(1:2,2)=cwaveprj1(iatm,2)%cp(1:2,jlmn)
     do ilmn=1,jlmn
      klmn=j0lmn+ilmn
      klmn_re=cplex*(klmn-1)+1
      cpi0(1:2,1)=cwaveprj (iatm,1)%cp(1:2,ilmn)
      cpi0(1:2,2)=cwaveprj (iatm,2)%cp(1:2,ilmn)
      cpi1(1:2,1)=cwaveprj1(iatm,1)%cp(1:2,ilmn)
      cpi1(1:2,2)=cwaveprj1(iatm,2)%cp(1:2,ilmn)
      ro11_re=cpj0(1,1)*cpi1(1,1)+cpj0(2,1)*cpi1(2,1)+cpi0(1,1)*cpj1(1,1)+cpi0(2,1)*cpj1(2,1)
      ro22_re=cpj0(1,2)*cpi1(1,2)+cpj0(2,2)*cpi1(2,2)+cpi0(1,2)*cpj1(1,2)+cpi0(2,2)*cpj1(2,2)
      pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight*(ro11_re+ro22_re)
      if (nspden>1) then
       ro12_re=cpj0(1,2)*cpi1(1,1)+cpj0(2,2)*cpi1(2,1)+cpi0(1,1)*cpj1(1,2)+cpi0(2,1)*cpj1(2,2)
       ro12_im=cpj0(1,2)*cpi1(2,1)-cpi1(1,2)*cpj0(2,1)+cpi0(1,1)*cpj1(2,2)-cpj1(1,1)*cpi0(2,2)
       ro21_re=cpj0(1,1)*cpi1(1,2)+cpj0(2,1)*cpi1(2,2)+cpi0(1,2)*cpj1(1,1)+cpi0(2,2)*cpj1(2,1)
       ro21_im=cpj0(1,1)*cpi1(2,2)-cpi1(1,1)*cpj0(2,2)+cpi0(1,2)*cpj1(2,1)-cpj1(1,2)*cpi0(2,1)
       pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight*(ro11_re-ro22_re)
       pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_re+ro21_re)
       pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro21_im-ro12_im)
      end if
      if (cplex==2) then
       klmn_im=klmn_re+1
       ro11_im=cpj0(1,1)*cpi1(2,1)-cpi1(1,1)*cpj0(2,1)+cpi0(1,1)*cpj1(2,1)-cpj1(1,1)*cpi0(2,1)
       ro22_im=cpj0(1,2)*cpi1(2,2)-cpi1(1,2)*cpj0(2,2)+cpi0(1,2)*cpj1(2,2)-cpj1(1,2)*cpi0(2,2)
       pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight*(ro11_im+ro22_im)
       if (nspden>1) then
        pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight*(ro11_im-ro22_im)
        pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_im+ro21_im)
        pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro12_re-ro21_re)
       end if
      end if
     end do
    end do
   end do
  end if

! Accumulate (n,k) contribution to rhoij1
! due to derivative of projectors
  if (ipert<=natom) then
   if (nspinor==1) then
    do iatm=1,dimpaw
     iatom=min(atindx1(iatm),dimpaw)
     do jlmn=1,pawrhoij(iatom)%lmn_size
      j0lmn=jlmn*(jlmn-1)/2
      cpj0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,jlmn)
      dcpj0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,jlmn)
      do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       klmn_re=cplex*(klmn-1)+1
       cpi0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,ilmn)
       dcpi0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,ilmn)
       pawrhoij(iatom)%rhoij_(klmn_re,isppol)=pawrhoij(iatom)%rhoij_(klmn_re,isppol) &
&       +weight*(dcpi0(1,1,1)*cpj0(1,1) +dcpi0(2,1,1)*cpj0(2,1) &
&       +cpi0(1,1) *dcpj0(1,1,1)+cpi0(2,1) *dcpj0(2,1,1))
       if (cplex==2) then
        klmn_im=klmn_re+1
        pawrhoij(iatom)%rhoij_(klmn_im,isppol)=pawrhoij(iatom)%rhoij_(klmn_im,isppol) &
&        +weight*(dcpi0(1,1,1)*cpj0(2,1)-dcpi0(2,1,1)*cpj0(1,1)+cpi0(1,1)*dcpj0(2,1,1)-cpi0(2,1)*dcpj0(1,1,1))
       end if
      end do
     end do
    end do
   else ! nspinor=2
    do iatm=1,dimpaw
     iatom=min(atindx1(iatm),dimpaw)
     do jlmn=1,pawrhoij(iatom)%lmn_size
      j0lmn=jlmn*(jlmn-1)/2
      cpj0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,jlmn)
      dcpj0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,jlmn)
      cpj0 (1:2,2)  =cwaveprj(iatm,2)%cp (1:2  ,jlmn)
      dcpj0(1:2,2,1)=cwaveprj(iatm,2)%dcp(1:2,1,jlmn)
      do ilmn=1,jlmn
       klmn=j0lmn+ilmn
       klmn_re=cplex*(klmn-1)+1
       cpi0 (1:2,1)  =cwaveprj(iatm,1)%cp (1:2  ,ilmn)
       dcpi0(1:2,1,1)=cwaveprj(iatm,1)%dcp(1:2,1,ilmn)
       cpi0 (1:2,2)  =cwaveprj(iatm,2)%cp (1:2  ,ilmn)
       dcpi0(1:2,2,1)=cwaveprj(iatm,2)%dcp(1:2,1,ilmn)
       ro11_re=dcpi0(1,1,1)*cpj0(1,1)+dcpi0(2,1,1)*cpj0(2,1)+cpi0(1,1)*dcpj0(1,1,1)+cpi0(2,1)*dcpj0(2,1,1)
       ro22_re=dcpi0(1,2,1)*cpj0(1,2)+dcpi0(2,2,1)*cpj0(2,2)+cpi0(1,2)*dcpj0(1,2,1)+cpi0(2,2)*dcpj0(2,2,1)
       pawrhoij(iatom)%rhoij_(klmn_re,1)=pawrhoij(iatom)%rhoij_(klmn_re,1)+weight*(ro11_re+ro22_re)
       if (nspden>1) then
        ro12_re=dcpi0(1,1,1)*cpj0(1,2)+dcpi0(2,1,1)*cpj0(2,2)+cpi0(1,1)*dcpj0(1,2,1)+cpi0(2,1)*dcpj0(2,2,1)
        ro12_im=dcpi0(1,1,1)*cpj0(2,2)-dcpi0(2,1,1)*cpj0(1,2)+cpi0(1,1)*dcpj0(2,2,1)-cpi0(2,1)*dcpj0(1,2,1)
        ro21_re=dcpi0(1,2,1)*cpj0(1,1)+dcpi0(2,2,1)*cpj0(2,1)+cpi0(1,2)*dcpj0(1,1,1)+cpi0(2,2)*dcpj0(2,1,1)
        ro21_im=dcpi0(1,2,1)*cpj0(2,1)-dcpi0(2,2,1)*cpj0(1,1)+cpi0(1,2)*dcpj0(2,1,1)-cpi0(2,2)*dcpj0(1,1,1)
        pawrhoij(iatom)%rhoij_(klmn_re,4)=pawrhoij(iatom)%rhoij_(klmn_re,4)+weight*(ro11_re-ro22_re)
        pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_re+ro21_re)
        pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro21_im-ro12_im)
       end if
       if (cplex==2) then
        klmn_im=klmn_re+1
        ro11_im=dcpi0(1,1,1)*cpj0(2,1)-dcpi0(2,1,1)*cpj0(1,1)+cpi0(1,1)*dcpj0(2,1,1)-cpi0(2,1)*dcpj0(1,1,1)
        ro22_im=dcpi0(1,2,1)*cpj0(2,2)-dcpi0(2,2,1)*cpj0(1,2)+cpi0(1,2)*dcpj0(2,2,1)-cpi0(2,2)*dcpj0(1,2,1)
        pawrhoij(iatom)%rhoij_(klmn_im,1)=pawrhoij(iatom)%rhoij_(klmn_im,1)+weight*(ro11_im+ro22_im)
        if (nspden>1) then
         pawrhoij(iatom)%rhoij_(klmn_im,4)=pawrhoij(iatom)%rhoij_(klmn_im,4)+weight*(ro11_im-ro22_im)
         pawrhoij(iatom)%rhoij_(klmn_re,2)=pawrhoij(iatom)%rhoij_(klmn_re,2)+weight*(ro12_im+ro21_im)
         pawrhoij(iatom)%rhoij_(klmn_re,3)=pawrhoij(iatom)%rhoij_(klmn_re,3)+weight*(ro12_re-ro21_re)
        end if
       end if
      end do
     end do
    end do
   end if
  end if

 else if (option==3) then

! ==================================================================
! === OPTION 3: Accumulate (n,k) contribution to drhoij/dr =========
! ==================================================================

  if (nspinor==1) then
   do iatm=1,natom
    iatom=atindx1(iatm)
    do jlmn=1,pawrhoij(iatom)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     cpj0(1:cplex,1)     =cwaveprj(iatm,1)%cp (1:cplex,jlmn)
     dcpj0(1:cplex,1,1:3)=cwaveprj(iatm,1)%dcp(1:cplex,1:3,jlmn)
     do ilmn=1,jlmn
      klmn=j0lmn+ilmn
      cpi0(1:cplex,1)     =cwaveprj(iatm,1)%cp (1:cplex,ilmn)
      dcpi0(1:cplex,1,1:3)=cwaveprj(iatm,1)%dcp(1:cplex,1:3,ilmn)
      do mu=1,3
       grhoij_tmp=zero
       do iplex=1,cplex
        grhoij_tmp=grhoij_tmp+dcpi0(iplex,1,mu)*cpj0(iplex,1)+cpi0(iplex,1)*dcpj0(iplex,1,mu)
       end do
       pawrhoij(iatom)%grhoij(mu,klmn,isppol)=pawrhoij(iatom)%grhoij(mu,klmn,isppol)+weight*grhoij_tmp
      end do
     end do
    end do
   end do
  else ! nspinor=2
   do iatm=1,natom
    iatom=atindx1(iatm)
    do jlmn=1,pawrhoij(iatom)%lmn_size
     j0lmn=jlmn*(jlmn-1)/2
     cpj0(1:cplex,1)     =cwaveprj(iatm,1)%cp (1:cplex,jlmn)
     cpj0(1:cplex,2)     =cwaveprj(iatm,2)%cp (1:cplex,jlmn)
     dcpj0(1:cplex,1,1:3)=cwaveprj(iatm,1)%dcp(1:cplex,1:3,jlmn)
     dcpj0(1:cplex,2,1:3)=cwaveprj(iatm,2)%dcp(1:cplex,1:3,jlmn)
     do ilmn=1,jlmn
      klmn=j0lmn+ilmn
      cpi0(1:cplex,1)     =cwaveprj(iatm,1)%cp (1:cplex,ilmn)
      cpi0(1:cplex,2)     =cwaveprj(iatm,2)%cp (1:cplex,ilmn)
      dcpi0(1:cplex,1,1:3)=cwaveprj(iatm,1)%dcp(1:cplex,1:3,ilmn)
      dcpi0(1:cplex,2,1:3)=cwaveprj(iatm,2)%dcp(1:cplex,1:3,ilmn)
      do mu=1,3
       ro11_re=dcpi0(1,1,mu)*cpj0(1,1)+cpi0(1,1)*dcpj0(1,1,mu)
       ro22_re=dcpi0(1,2,mu)*cpj0(1,2)+cpi0(1,2)*dcpj0(1,2,mu)
       pawrhoij(iatom)%grhoij(mu,klmn,1)=pawrhoij(iatom)%grhoij(mu,klmn,1)+weight*(ro11_re+ro22_re)
       if (nspden>1) then
        ro12_re=dcpi0(1,1,mu)*cpj0(1,2)+cpi0(1,1)*dcpj0(1,2,mu)
        ro21_re=dcpi0(1,2,mu)*cpj0(1,1)+cpi0(1,2)*dcpj0(1,1,mu)
        pawrhoij(iatom)%grhoij(mu,klmn,4)=pawrhoij(iatom)%grhoij(mu,klmn,4)+weight*(ro11_re-ro22_re)
        pawrhoij(iatom)%grhoij(mu,klmn,2)=pawrhoij(iatom)%grhoij(mu,klmn,2)+weight*(ro12_re+ro21_re)
       end if
       if (cplex==2) then
        ro11_im=dcpi0(2,1,mu)*cpj0(2,1)+cpi0(2,1)*dcpj0(2,1,mu)
        ro22_im=dcpi0(2,2,mu)*cpj0(2,2)+cpi0(2,2)*dcpj0(2,2,mu)
        pawrhoij(iatom)%grhoij(mu,klmn,1)=pawrhoij(iatom)%grhoij(mu,klmn,1)+weight*(ro11_im+ro22_im)
        if (nspden>1) then
         ro12_re=dcpi0(2,1,mu)*cpj0(2,2)+cpi0(2,1)*dcpj0(2,2,mu)
         ro21_re=dcpi0(2,2,mu)*cpj0(2,1)+cpi0(2,2)*dcpj0(2,1,mu)
         ro12_im=dcpi0(1,1,mu)*cpj0(2,2)+cpi0(1,1)*dcpj0(2,2,mu)-dcpi0(2,1,mu)*cpj0(1,2)-cpi0(2,1)*dcpj0(1,2,mu)
         ro21_im=dcpi0(1,2,mu)*cpj0(2,1)+cpi0(1,2)*dcpj0(2,1,mu)-dcpi0(2,2,mu)*cpj0(1,1)-cpi0(2,2)*dcpj0(1,1,mu)
         pawrhoij(iatom)%grhoij(mu,klmn,4)=pawrhoij(iatom)%grhoij(mu,klmn,4)+weight*(ro11_im-ro22_im)
         pawrhoij(iatom)%grhoij(mu,klmn,2)=pawrhoij(iatom)%grhoij(mu,klmn,2)+weight*(ro12_re+ro21_re)
         pawrhoij(iatom)%grhoij(mu,klmn,3)=pawrhoij(iatom)%grhoij(mu,klmn,3)+weight*(ro21_im-ro12_im)
        end if
       end if
      end do
     end do
    end do
   end do
  end if

! End
 end if ! option

end subroutine pawaccrhoij
!!***
