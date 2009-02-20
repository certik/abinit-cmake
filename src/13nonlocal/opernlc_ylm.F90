!{\src2tex{textfont=tt}}
!!****f* ABINIT/opernlc_ylm
!! NAME
!! opernlc_ylm
!!
!! FUNCTION
!! * Operate with the non-local part of the hamiltonian,
!!   in order to reduce projected scalars
!! * Operate with the non-local projectors and the overlap matrix,
!!   in order to reduce projected scalars
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms (gives the absolute index of
!!                 an atom from its rank in a block of atoms)
!!  cplex=1 if <p_lmn|c> scalars are real (equivalent to istwfk>1)
!!        2 if <p_lmn|c> scalars are complex
!!  cplex_enl=1 if enl factors are real, 2 if they are complex
!!  cplex_fac=1 if gxfac scalars are real, 2 if gxfac scalars are complex
!!  dgxdt(cplex,ndgxdt,nlmn,nincat)=grads of projected scalars (only if optder>0)
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  enl(cplex_enl*dimenl1,dimenl2,nspinor**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1, 2 or 4 ====
!!                      (Real or complex, hermitian) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  gx(cplex,nlmn,nincat*abs(enl_opt))= projected scalars
!!  iatm=absolute rank of first atom of the current block of atoms
!!  indlmn(6,nlmn)= array giving l,m,n,lm,ln,s for i=lmn
!!  itypat=type of atoms
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!  natom=number of atoms in cell
!!  ndgxdt=second dimension of dgxdt
!!  ndgxdtfac=second dimension of dgxdtfac
!!  nincat=number of atoms in the subset here treated
!!  nlmn=number of (l,m,n) numbers for current type of atom
!!  nspinor=number of spinorial components of the wavefunctions
!!  optder=0=only gxfac is computed, 1=both gxfac and dgxdtfac are computed
!!  paw_opt= define the nonlocal operator concerned with:
!!           paw_opt=0 : Norm-conserving Vnl (use of Kleinman-Bylander ener.)
!!           paw_opt=1 : PAW nonlocal part of H (use of Dij coeffs)
!!           paw_opt=2 : PAW: (Vnl-lambda.Sij) (Sij=overlap matrix)
!!           paw_opt=3 : PAW overlap matrix (Sij)
!!           paw_opt=4 : both PAW nonlocal part of H (Dij) and overlap matrix (Sij)
!!  sij(nlm*(nlmn+1)/2)=overlap matrix components (only if paw_opt=2, 3 or 4)
!!
!! OUTPUT
!!  if (paw_opt=0, 1, 2 or 4)
!!    gxfac(cplex_fac,nlmn,nincat,nspinor)= reduced projected scalars related to Vnl (NL operator)
!!  if (paw_opt=3 or 4)
!!    gxfac_sij(cplex,nlmn,nincat,nspinor)= reduced projected scalars related to Sij (overlap)
!!  if (optder==1.and.paw_opt=0, 1, 2 or 4)
!!    dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Vnl (NL operator)
!!  if (optder==1.and.paw_opt=3 or 4)
!!    dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor)= gradients of gxfac related to Sij (overlap)
!!
!! NOTES
!! This routine operates for one type of atom, and within this given type of atom,
!! for a subset of at most nincat atoms.
!!
!! PARENTS
!!      nonlop_ylm
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine opernlc_ylm(atindx1,cplex,cplex_enl,cplex_fac,dgxdt,dgxdtfac,dgxdtfac_sij,&
&                      dimenl1,dimenl2,enl,gx,gxfac,gxfac_sij,iatm,indlmn,itypat,&
&                      lambda,natom,ndgxdt,ndgxdtfac,nincat,nlmn,nspinor,optder,paw_opt,sij)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,cplex_enl,cplex_fac,dimenl1,dimenl2,iatm,itypat,natom,ndgxdt,ndgxdtfac
 integer,intent(in) :: nincat,nlmn,nspinor,optder,paw_opt
 real(dp) :: lambda
!arrays
 integer,intent(in) :: atindx1(natom),indlmn(6,nlmn)
 real(dp),intent(in) :: dgxdt(cplex,ndgxdt,nlmn,nincat,nspinor)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinor**2)
 real(dp),intent(in) :: gx(cplex,nlmn,nincat,nspinor)
 real(dp),intent(in) :: sij(((paw_opt+1)/3)*nlmn*(nlmn+1)/2)
 real(dp),intent(out) :: dgxdtfac(cplex_fac,ndgxdtfac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: dgxdtfac_sij(cplex,ndgxdtfac,nlmn,nincat,nspinor*(paw_opt/3))
 real(dp),intent(out) :: gxfac(cplex_fac,nlmn,nincat,nspinor)
 real(dp),intent(out) :: gxfac_sij(cplex,nlmn,nincat,nspinor*(paw_opt/3))

!Local variables-------------------------------
!Arrays
!scalars
 integer :: ia,ijlmn,ilm,ilmn,iln,index_enl,iplex,jplex,ispinor,j0lmn,jjlmn,jlm,jlmn,jspinor,mu
 real(dp) :: sijr
!arrays
 real(dp), parameter :: facti(2)=(/1.d0,-1.d0/)
 real(dp) :: enl_(2),gxfi(2),gxi(cplex),gxj(cplex)
 real(dp),allocatable :: gxfj(:,:)

! *************************************************************************

!Accumulate gxfac related to non-local operator (Norm-conserving)
!-------------------------------------------------------------------
 if (paw_opt==0) then                  ! Enl is E(Kleinman-Bylander)
  if (cplex_enl==2    ) stop "opernlc_ylm: BUG - invalid cplex_enl=2 !"
  if (cplex_fac/=cplex) stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex) !"
  do ispinor=1,nspinor
   do ia=1,nincat
    do ilmn=1,nlmn
     iln=indlmn(5,ilmn)
     enl_(1)=enl(iln,itypat,ispinor)
     gxfac(1:cplex,ilmn,ia,ispinor)=enl_(1)*gx(1:cplex,ilmn,ia,ispinor)
    end do
   end do
  end do
 end if

!Accumulate gxfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
 if (paw_opt==1.or.paw_opt==2.or.paw_opt==4) then        ! Enl is psp strength Dij
  gxfac(1:cplex_fac,1:nlmn,1:nincat,1:nspinor)=zero      ! or (Dij-lambda.Sij)
! === Diagonal term(s) (up-up, down-down)
! 1-Enl is real
  if (cplex_enl==1) then
   do ispinor=1,nspinor
    do ia=1,nincat
     index_enl=atindx1(iatm+ia)
     do jlmn=1,nlmn
      j0lmn=jlmn*(jlmn-1)/2
      jjlmn=j0lmn+jlmn
      enl_(1)=enl(jjlmn,index_enl,ispinor)
      if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
      gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
      gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
      do ilmn=1,jlmn-1
       ijlmn=j0lmn+ilmn
       enl_(1)=enl(ijlmn,index_enl,ispinor)
       if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
       gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
       gxfac(1:cplex,ilmn,ia,ispinor)=gxfac(1:cplex,ilmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
       gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxi(1:cplex)
      end do
     end do
    end do
   end do
!  2-Enl is complex
  else
   if (cplex_fac/=cplex_enl) stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex_enl !"
   do ispinor=1,nspinor
    do ia=1,nincat
     index_enl=atindx1(iatm+ia)
     do jlmn=1,nlmn
      j0lmn=jlmn*(jlmn-1)/2
      jjlmn=j0lmn+jlmn
      enl_(1)=enl(2*jjlmn-1,index_enl,ispinor)
      if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
      gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
      gxfac(1:cplex,jlmn,ia,ispinor)=gxfac(1:cplex,jlmn,ia,ispinor)+enl_(1)*gxj(1:cplex)
      do ilmn=1,jlmn-1
       ijlmn=j0lmn+ilmn
       enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor)
       if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
       gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
       do iplex=1,cplex
        jplex=3-iplex
        gxfac(iplex,ilmn,ia,ispinor)=gxfac(iplex,ilmn,ia,ispinor)+             enl_(1)*gxj(iplex)
        gxfac(iplex,jlmn,ia,ispinor)=gxfac(iplex,jlmn,ia,ispinor)+             enl_(1)*gxi(iplex)
        gxfac(jplex,ilmn,ia,ispinor)=gxfac(jplex,ilmn,ia,ispinor)+facti(iplex)*enl_(2)*gxj(iplex)
        gxfac(jplex,jlmn,ia,ispinor)=gxfac(jplex,jlmn,ia,ispinor)-facti(iplex)*enl_(2)*gxi(iplex)
       end do
      end do
     end do
    end do
   end do
  end if

! === Off-diagonal term(s) (up-down, down-up)
  if (nspinor==2) then
   if (cplex_enl/=2) stop "opernlc_ylm: BUG - invalid cplex_enl/=2 !"
   if (cplex_fac/=2) stop "opernlc_ylm: BUG - invalid cplex_fac/=2 !"
   do ispinor=1,nspinor
    jspinor=3-ispinor
    do ia=1,nincat
     index_enl=atindx1(iatm+ia)
     do jlmn=1,nlmn
      j0lmn=jlmn*(jlmn-1)/2
      jjlmn=j0lmn+jlmn
      enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor)
      gxi(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
      gxj(1:cplex)=gx(1:cplex,jlmn,ia,jspinor)
      do iplex=1,cplex
       jplex=3-iplex
       gxfac(iplex,jlmn,ia,jspinor)=gxfac(iplex,jlmn,ia,jspinor)+             enl_(1)*gxi(iplex)
       gxfac(jplex,jlmn,ia,jspinor)=gxfac(jplex,jlmn,ia,jspinor)-facti(iplex)*enl_(2)*gxi(iplex)
      end do
      do ilmn=1,jlmn-1
       ijlmn=j0lmn+ilmn
       enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
       gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
       do iplex=1,cplex
        jplex=3-iplex
        gxfac(iplex,ilmn,ia,ispinor)=gxfac(iplex,ilmn,ia,ispinor)+             enl_(1)*gxj(iplex)
        gxfac(iplex,jlmn,ia,jspinor)=gxfac(iplex,jlmn,ia,jspinor)+             enl_(1)*gxi(iplex)
        gxfac(jplex,ilmn,ia,ispinor)=gxfac(jplex,ilmn,ia,ispinor)+facti(iplex)*enl_(2)*gxj(iplex)
        gxfac(jplex,jlmn,ia,jspinor)=gxfac(jplex,jlmn,ia,jspinor)-facti(iplex)*enl_(2)*gxi(iplex)
       end do
      end do
     end do
    end do
   end do
  end if
 end if

!Accumulate gxfac related to overlap (Sij) (PAW)
!------------------------------------------- ------------------------
 if (paw_opt==3.or.paw_opt==4) then                    ! Use Sij, overlap contribution
  gxfac_sij(1:cplex,1:nlmn,1:nincat,1:nspinor)=zero
  do ispinor=1,nspinor
   do ia=1,nincat
    do jlmn=1,nlmn
     j0lmn=jlmn*(jlmn-1)/2
     jjlmn=j0lmn+jlmn
     jlm=indlmn(4,jlmn)
     sijr=sij(jjlmn);gxj(1:cplex)=gx(1:cplex,jlmn,ia,ispinor)
     gxfac_sij(1:cplex,jlmn,ia,ispinor)=gxfac_sij(1:cplex,jlmn,ia,ispinor)+sijr*gxj(1:cplex)
     do ilmn=1,jlmn-1
      ilm=indlmn(4,ilmn)
      if (ilm==jlm) then
       ijlmn=j0lmn+ilmn
       sijr=sij(ijlmn)
       gxi(1:cplex)=gx(1:cplex,ilmn,ia,ispinor)
       gxfac_sij(1:cplex,ilmn,ia,ispinor)=gxfac_sij(1:cplex,ilmn,ia,ispinor)+sijr*gxj(1:cplex)
       gxfac_sij(1:cplex,jlmn,ia,ispinor)=gxfac_sij(1:cplex,jlmn,ia,ispinor)+sijr*gxi(1:cplex)
      end if
     end do
    end do
   end do
  end do
 end if

!Accumulate dgxdtfac related to nonlocal operator (Norm-conserving)
!-------------------------------------------------------------------
 if (optder==1.and.paw_opt==0) then    ! Enl is E(Kleinman-Bylander)
  if (cplex_enl/=1    ) stop "opernlc_ylm: BUG - invalid cplex_enl=2 !"
  if (cplex_fac/=cplex) stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex !"
  do ispinor=1,nspinor
   do ia=1,nincat
    do ilmn=1,nlmn
     iln=indlmn(5,ilmn)
     enl_(1)=enl(iln,itypat,ispinor)
     do mu=1,ndgxdtfac
      dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)=enl_(1)*dgxdt(1:cplex,mu,ilmn,ia,ispinor)
     end do
    end do
   end do
  end do
 end if

!Accumulate dgxdtfac related to nonlocal operator (PAW)
!-------------------------------------------------------------------
 if (optder==1.and.(paw_opt==1.or.paw_opt==2.or.paw_opt==4)) then  ! Enl is psp strength Dij
  allocate(gxfj(cplex,ndgxdtfac))                                  ! or (Dij-lambda.Sij)
  dgxdtfac(1:cplex,1:ndgxdtfac,1:nlmn,1:nincat,1:nspinor)=zero
! === Diagonal term(s) (up-up, down-down)
! 1-Enl is real
  if (cplex_enl==1) then
   do ispinor=1,nspinor
    do ia=1,nincat
     index_enl=atindx1(iatm+ia)
     do jlmn=1,nlmn
      j0lmn=jlmn*(jlmn-1)/2
      jjlmn=j0lmn+jlmn
      enl_(1)=enl(jjlmn,index_enl,ispinor)
      if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
      do mu=1,ndgxdtfac
       gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
       dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
      end do
      do ilmn=1,jlmn-1
       ijlmn=j0lmn+ilmn
       enl_(1)=enl(ijlmn,index_enl,ispinor)
       if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
       do mu=1,ndgxdtfac
        gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
        dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)=dgxdtfac(1:cplex,mu,ilmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
        dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfi(1:cplex)
       end do
      end do
     end do
    end do
   end do
!  2-Enl is complex
  else
   if (cplex_fac/=cplex_enl) stop "opernlc_ylm: BUG - invalid cplex_fac/=cplex_enl !"
   do ispinor=1,nspinor
    do ia=1,nincat
     index_enl=atindx1(iatm+ia)
     do jlmn=1,nlmn
      j0lmn=jlmn*(jlmn-1)/2
      jjlmn=j0lmn+jlmn
      enl_(1)=enl(2*jjlmn-1,index_enl,ispinor)
      if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(jjlmn)
      do mu=1,ndgxdtfac
       gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
       dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac(1:cplex,mu,jlmn,ia,ispinor)+enl_(1)*gxfj(1:cplex,mu)
      end do
      do ilmn=1,jlmn-1
       ijlmn=j0lmn+ilmn
       enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,ispinor)
       if (paw_opt==2) enl_(1)=enl_(1)-lambda*sij(ijlmn)
       do mu=1,ndgxdtfac
        gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
        do iplex=1,cplex
         jplex=3-iplex
         dgxdtfac(iplex,mu,ilmn,ia,ispinor)=dgxdtfac(iplex,mu,ilmn,ia,ispinor)+             enl_(1)*gxfj(iplex,mu)
         dgxdtfac(iplex,mu,jlmn,ia,ispinor)=dgxdtfac(iplex,mu,jlmn,ia,ispinor)+             enl_(1)*gxfi(iplex)
         dgxdtfac(jplex,mu,ilmn,ia,ispinor)=dgxdtfac(jplex,mu,ilmn,ia,ispinor)+facti(iplex)*enl_(2)*gxfj(iplex,mu)
         dgxdtfac(jplex,mu,jlmn,ia,ispinor)=dgxdtfac(jplex,mu,jlmn,ia,ispinor)-facti(iplex)*enl_(2)*gxfi(iplex)
        end do
       end do
      end do
     end do
    end do
   end do
  end if
! === Off-diagonal term(s) (up-down, down-up)
  if (nspinor==2) then
   if (cplex_enl/=2) stop "opernlc_ylm: BUG - invalid cplex_enl/=2 !"
   if (cplex_fac/=2) stop "opernlc_ylm: BUG - invalid cplex_fac/=2 !"
   do ispinor=1,nspinor
    jspinor=3-ispinor
    do ia=1,nincat
     index_enl=atindx1(iatm+ia)
     do jlmn=1,nlmn
      j0lmn=jlmn*(jlmn-1)/2
      jjlmn=j0lmn+jlmn
      enl_(1:2)=enl(2*jjlmn-1:2*jjlmn,index_enl,2+ispinor)
      do mu=1,ndgxdtfac
       gxfi(1:cplex)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
       gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,jspinor)
       do iplex=1,cplex
        jplex=3-iplex
        dgxdtfac(iplex,mu,jlmn,ia,jspinor)=dgxdtfac(iplex,mu,jlmn,ia,jspinor)+             enl_(1)*gxfi(iplex)
        dgxdtfac(jplex,mu,jlmn,ia,jspinor)=dgxdtfac(jplex,mu,jlmn,ia,jspinor)-facti(iplex)*enl_(2)*gxfi(iplex)
       end do
      end do
      do ilmn=1,jlmn-1
       ijlmn=j0lmn+ilmn
       enl_(1:2)=enl(2*ijlmn-1:2*ijlmn,index_enl,2+ispinor)
       do mu=1,ndgxdtfac
        gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
        do iplex=1,cplex
         jplex=3-iplex
         dgxdtfac(iplex,mu,ilmn,ia,ispinor)=dgxdtfac(iplex,mu,ilmn,ia,ispinor)+             enl_(1)*gxfj(iplex,mu)
         dgxdtfac(iplex,mu,jlmn,ia,jspinor)=dgxdtfac(iplex,mu,jlmn,ia,jspinor)+             enl_(1)*gxfi(iplex)
         dgxdtfac(jplex,mu,ilmn,ia,ispinor)=dgxdtfac(jplex,mu,ilmn,ia,ispinor)+facti(iplex)*enl_(2)*gxfj(iplex,mu)
         dgxdtfac(jplex,mu,jlmn,ia,jspinor)=dgxdtfac(jplex,mu,jlmn,ia,jspinor)-facti(iplex)*enl_(2)*gxfi(iplex)
        end do
       end do
      end do
     end do
    end do
   end do
  end if
  deallocate(gxfj)
 end if

!Accumulate dgxdtfac related to overlap (Sij) (PAW)
!-------------------------------------------------------------------
 if (optder==1.and.(paw_opt==3.or.paw_opt==4)) then  ! Use Sij, overlap contribution
  allocate(gxfj(cplex,ndgxdtfac))
  dgxdtfac_sij(1:cplex,1:ndgxdtfac,1:nlmn,1:nincat,1:nspinor)=zero
  do ispinor=1,nspinor
   do ia=1,nincat
    do jlmn=1,nlmn
     j0lmn=jlmn*(jlmn-1)/2
     jjlmn=j0lmn+jlmn
     sijr=sij(jjlmn)
     do mu=1,ndgxdtfac
      gxfj(1:cplex,mu)=dgxdt(1:cplex,mu,jlmn,ia,ispinor)
      dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
     end do
     do ilmn=1,jlmn-1
      ijlmn=j0lmn+ilmn
      sijr=sij(ijlmn)
      do mu=1,ndgxdtfac
       gxfi(1:cplex)=dgxdt(1:cplex,mu,ilmn,ia,ispinor)
       dgxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,ilmn,ia,ispinor)+sijr*gxfj(1:cplex,mu)
       dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)=dgxdtfac_sij(1:cplex,mu,jlmn,ia,ispinor)+sijr*gxfi(1:cplex)
      end do
     end do
    end do
   end do
  end do
  deallocate(gxfj)
 end if

end subroutine opernlc_ylm
!!***
