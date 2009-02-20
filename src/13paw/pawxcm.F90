!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxcm
!! NAME
!! pawxcm
!!
!! FUNCTION
!! PAW only
!! Start from the density or spin-density, and compute xc correlation
!! potential and energies inside a paw sphere.
!! LDA+GGA - USE A DEVELOPMENT OF THE DENSITY OVER (L,M) MOMENTS
!! Driver of XC functionals.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT, GJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc_coll
!!
!! INPUTS
!!  corexc(pawrad%mesh_size)=core density on radial grid
!!  exexch= choice of local exact exchange. Active if exexch=3
!!  ixc= choice of exchange-correlation scheme (see above and below)
!!  lm_size=size of density array rhor (see below)
!!  lmselect(lm_size)=select the non-zero LM-moments of input density rhor
!!  nhat(pawrad%mesh_size,lm_size,nspden)=compensation density
!!                                        (total in 1st half and spin-up in 2nd half if nspden=2)
!!  nspden=number of spin-density components
!!  option=0 compute both XC energies (direct+double-counting) and potential
!!         1 compute only XC potential
!!         2 compute only XC energies (direct+double-counting)
!!         3 compute only XC energy by direct scheme
!!         4 compute only XC energy by direct scheme for spherical part of the density
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  pawxcdev=order of Vxc development
!!  rhor(pawrad%mesh_size,lm_size,nspden)=electron density in real space in electrons/bohr**3
!!                                       (total in 1st half and spin-up in 2nd half if nspden=2)
!!  usecore= 1 if core density has to be used in Exc/Vxc ; 0 otherwise
!!  usexcnhat= 0 if compensation density does not have to be used
!!             1 if compensation density has to be used in double counting energy term only
!!             2 if compensation density (nhat) has to be used in Exc/Vxc and double counting energy term
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  == if option==0, 2, 3, or 4 ==
!!    enxc=returned exchange and correlation energy (hartree)
!!  == if option==0 or 2 ==
!!    enxcdc=returned exchange-cor. contribution to double-counting energy
!!  == if option==0 or 1 ==
!!    vxc(pawrad%mesh_size,lm_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!
!! PARENTS
!!      pawdenpot,psp7in
!!
!! CHILDREN
!!      mkdenpos,pawxcsph,simp_gen,size_dvxc,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pawxcm(corexc,enxc,enxcdc,exexch,ixc,lm_size,lmselect,nhat,nspden,option,&
&                  pawang,pawrad,pawxcdev,rhor,usecore,usexcnhat,vxc,xclevel)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13paw, except_this_one => pawxcm
 use interfaces_13xc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exexch,ixc,lm_size,nspden,option,pawxcdev,usecore
 integer,intent(in) :: usexcnhat,xclevel
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc(pawrad%mesh_size)
 real(dp),intent(in) :: nhat(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(in) :: rhor(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(out) :: vxc(pawrad%mesh_size,lm_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ilm1,ilm2,ir,ir1,ir2,isel,ispden,iwarn,jr,ndvxc,ngr2,ngrad,nrad
 integer :: nspden_updn,nspgrad,nvxcdgr,order
 real(dp),parameter :: delta=1.d-4
 real(dp) :: dvxc1,dvxc2,dvxc3,dvxc4,dvxca,dvxcb,dvxcc,dvxcd
 real(dp) :: fact,invsqfpi,invsqfpi2,m_norm,sqfpi,sqfpi2
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: d1vxc(:,:),d2vxc(:,:),exc_(:),exci(:),ff(:),gg(:),m_norm_inv(:)
 real(dp),allocatable :: rho_(:,:),rho_updn(:,:,:),rhoinv(:,:),rhosph(:,:)
 real(dp),allocatable :: v0sum(:,:),v1sum(:,:),v2sum(:,:,:),vxc1(:,:),vxc2(:,:)
 real(dp),allocatable :: vxcdn1(:,:),vxcdn2(:,:),vxci(:,:)
!************************************************************************

!DEBUG
!write(6,*)' pawxcm : enter with option, nspden ',option,nspden
!ENDDEBUG
 call timab(81,1,tsec)

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

!Arrays dimensions and constants
 iwarn=0
 order=1
 nspden_updn=min(nspden,2)
 ngrad=1;if (xclevel==2) ngrad=2 ! ngrad=1 is for LDAs or LSDs; ngrad=2 is for GGAs
 nspgrad=nspden_updn*ngrad;if(nspden_updn==2.and.ngrad==2) nspgrad=5
 nrad=pawrad%mesh_size
 sqfpi=sqrt(four_pi);sqfpi2=half*sqfpi
 invsqfpi=one/sqfpi;invsqfpi2=half*invsqfpi

!Compute sizes of arrays
 call size_dvxc(ixc,ndvxc,ngr2,nspden_updn,nvxcdgr,order)

!Initializations of output arrays
 if (option/=1) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option<3) vxc(:,:,:)=zero

 if (xclevel==0) then ! No xc at all is applied (usually for testing)
  write(message, '(a,a,a,a)' ) ch10,&
&  ' pawxcm : WARNING -',ch10,&
&  '  Note that no xc is applied (ixc=0).'
  call wrtout(06,message,'COLL')
  return
 end if

!----------------------------------------------------------------------
!----- Build several densities
!----------------------------------------------------------------------

!rho_updn contains the effective density used for XC
!with core density and/or compensation density eventually included
!-----------------------------------------------------------------
 allocate(rho_updn(nrad,lm_size,nspden))
 rho_updn(:,:,:)=rhor(:,:,:)
 if (usexcnhat==2) rho_updn(:,:,:)=rho_updn(:,:,:)+nhat(:,:,:)
 if (usecore==1) then
  if (nspden==1.or.nspden==4) then
   rho_updn(:,1,1)=rho_updn(:,1,1)+sqfpi*corexc(:)
  else if (nspden==2) then
   rho_updn(:,1,1)=rho_updn(:,1,1)+sqfpi*corexc(:)
   rho_updn(:,1,2)=rho_updn(:,1,2)+sqfpi2*corexc(:)
  end if
 end if

!In case of collinear magnetism, separate up and down contributions
 if (nspden==2) then
  allocate(ff(nrad))
  do ilm=1,lm_size
   ff(:)=rho_updn(:,ilm,2)
   rho_updn(:,ilm,2)=rho_updn(:,ilm,1)-ff(:)
   rho_updn(:,ilm,1)=ff(:)
  end do
  deallocate(ff)
 end if

!rhoSPH contains the spherical part of effective density
!(including Y00 spherical harmonic)
!-----------------------------------------------------------------
 allocate(rhosph(nrad,nspden_updn))

!Non-magnetic system: rhoSPH(;,1)=(1/2).rhoSPH_total
 if (nspden==1) then
  rhosph(:,1)=rho_updn(:,1,1)*invsqfpi2

! Collinear magnetism: rhoSPH = (rhoSPH_up, rhoSPH_dn)
 else if (nspden==2) then
  rhosph(:,1:2)=rho_updn(:,1,1:2)*invsqfpi

! Non-collinear magnetism: rhoSPH = (rhoSPH_up, rhoSPH_dn)
! obtained by rotating rho_updn
 else if (nspden==4) then
  allocate(m_norm_inv(nrad)) ! Store 1/norm of zero-order moment of magnetization
  do ir=1,nrad
   m_norm=sqrt(rho_updn(ir,1,2)**2+rho_updn(ir,1,3)**2+rho_updn(ir,1,4)**2)
   rhosph(ir,1)=(rho_updn(ir,1,1)+m_norm)*invsqfpi2
   rhosph(ir,2)=(rho_updn(ir,1,1)-m_norm)*invsqfpi2
   if (m_norm>abs(rho_updn(ir,1,1))*tol10+tol14) then
    m_norm_inv(ir)=one/m_norm
   else
    m_norm_inv(ir)=zero
   end if
  end do
 end if

!Make spherical density positive
 call mkdenpos(iwarn,nrad,nspden_updn,0,rhosph)

!----------------------------------------------------------------------
!----- Compute Exc(rhoSPH) and Vxc(rhoSPH)
!----------------------------------------------------------------------

 allocate(exci(nrad),vxci(nrad,nspden_updn))
 call pawxcsph(exci,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
& order,pawrad,rhosph,vxci,xclevel)

!----------------------------------------------------------------------
!----- Compute numerical derivatives of Vxc (by finite difference scheme)
!----------------------------------------------------------------------

 if (option/=4) then
  allocate(exc_(nrad),rho_(nrad,nspden_updn))

  if (nspden_updn==2) rho_(:,2)=rhosph(:,2)

! Compute Exc, Vxc for rho+delta_rho
  allocate(vxc1(nrad,nspden_updn))
  rho_(:,1)=(one+delta)*rhosph(:,1)
  call pawxcsph(exc_,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
&  order,pawrad,rho_,vxc1,xclevel)

! Compute Exc, Vxc for rho-delta_rho
  allocate(vxc2(nrad,nspden_updn))
  rho_(:,1)=(one-delta)*rhosph(:,1)
  call pawxcsph(exc_,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
&  order,pawrad,rho_,vxc2,xclevel)

! Additional terms for spin-polarized systems
  if (nspden_updn==2) then
   rho_(:,1)=rhosph(:,1)

!  Compute Exc, Vxc for rho+delta_rho_down
   allocate(vxcdn1(nrad,nspden_updn))
   rho_(:,2)=(one+delta)*rhosph(:,2)
   call pawxcsph(exc_,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
&   order,pawrad,rho_,vxcdn1,xclevel)

!  Compute Exc, Vxc for rho-delta_rho_down
   allocate(vxcdn2(nrad,nspden_updn))
   rho_(:,2)=(one-delta)*rhosph(:,2)
   call pawxcsph(exc_,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden_updn,nspgrad,nvxcdgr,&
&   order,pawrad,rho_,vxcdn2,xclevel)

  end if !nspden_updn==2
  deallocate(exc_,rho_)

! Store inverse of density finite step
  allocate(rhoinv(nrad,nspden_updn))
  fact=one/delta;if (nspden_updn==1) fact=half*fact
  do ispden=1,nspden_updn
   do ir=1,nrad
    if (rhosph(ir,ispden)>tol14) then
     rhoinv(ir,ispden)=fact/rhosph(ir,ispden)
    else
     rhoinv(ir,ispden)=zero
    end if
   end do
  end do

! Compute numerical first derivatives of Vxc (by finite difference scheme)
  allocate(d1vxc(nrad,2*nspden_updn-1))
! Non-magnetic system: compute dVxc/dn
  if (nspden==1) then
   d1vxc(1:nrad,1)=(vxc1(1:nrad,1)-vxc2(1:nrad,1))*half*rhoinv(1:nrad,1)
!  Collinear magnetism: compute dVxc_up/dn_up,dVxc_dn/dn_up,dVxc_dn/dn_dn
  else if (nspden==2) then
   d1vxc(1:nrad,1)=(vxc1(1:nrad,1)-vxc2(1:nrad,1))*half*rhoinv(1:nrad,1)
   d1vxc(1:nrad,2)=(vxc1(1:nrad,2)-vxc2(1:nrad,2))*half*rhoinv(1:nrad,1)
   d1vxc(1:nrad,3)=(vxcdn1(1:nrad,2)-vxcdn2(1:nrad,2))*half*rhoinv(1:nrad,2)
!  Non-collinear magnetism: compute 1/2 d(Vxc_up+Vxc_dn)/dn,1/2 d(Vxc_up-Vxc_dn)/dn
!  1/2 d(Vxc_up-Vxc_dn)/dm
  else if (nspden==4) then
   do ir=1,nrad
    fact=half*rhoinv(ir,1)
    dvxc1=(vxc1  (ir,1)-vxc2  (ir,1))*fact !dVxc_up/dn_up
    dvxc2=(vxc1  (ir,2)-vxc2  (ir,2))*fact !dVxc_dn/dn_up
    dvxc3=(vxcdn1(ir,2)-vxcdn2(ir,2))*fact !dVxc_dn/dn_dn
    dvxca=dvxc1+dvxc3;dvxcb=dvxc1-dvxc3;dvxcc=two*dvxc1 !Temporary terms
    d1vxc(ir,1)=quarter*(dvxca+dvxcc)  ! 1/2 d(Vxc_up+Vxc_dn)/dn
    d1vxc(ir,2)=quarter* dvxcb         ! 1/2 d(Vxc_up-Vxc_dn)/dn
    d1vxc(ir,3)=quarter*(dvxca-dvxcc)  ! 1/2 d(Vxc_up-Vxc_dn)/dm
   end do
  end if

! Compute numerical second derivatives of Vxc (by finite difference scheme)
  if (option<3.or.pawxcdev>1) then
   allocate(d2vxc(nrad,3*nspden_updn-2))
!  Non-magnetic system: compute d2Vxc/dn2
   if (nspden==1) then
    d2vxc(1:nrad,1)=(vxc1(1:nrad,1)+vxc2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,1)**2
!   Collinear magnetism: compute d2Vxc_up/dn_up2,d2Vxc_dn/dn_up2,d2Vxc_up/dn_dn2,d2Vxc_dn/dn_dn2
   else if (nspden==2) then
    d2vxc(1:nrad,1)=(vxc1(1:nrad,1)+vxc2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,1)**2
    d2vxc(1:nrad,2)=(vxc1(1:nrad,2)+vxc2(1:nrad,2)-two*vxci(1:nrad,2))*rhoinv(1:nrad,1)**2
    d2vxc(1:nrad,3)=(vxcdn1(1:nrad,1)+vxcdn2(1:nrad,1)-two*vxci(1:nrad,1))*rhoinv(1:nrad,2)**2
    d2vxc(1:nrad,4)=(vxcdn1(1:nrad,2)+vxcdn2(1:nrad,2)-two*vxci(1:nrad,2))*rhoinv(1:nrad,2)**2
!   Non-collinear magnetism: compute 1/2 d(2Vxc_up+Vxc_dn)/dn2,1/2 d2(Vxc_up-Vxc_dn)/dn2
!   1/2 d2(Vxc_up+Vxc_dn)/dm2,1/2 d2(Vxc_up-Vxc_dn)/dm2
   else if (nspden==4) then
    do ir=1,nrad
     fact=rhoinv(ir,1)**2
     dvxc1=(vxc1  (ir,1)+vxc2  (ir,1)-two*vxci(ir,1))*fact !d2Vxc_up/dn_up2
     dvxc2=(vxc1  (ir,2)+vxc2  (ir,2)-two*vxci(ir,2))*fact !d2Vxc_dn/dn_up2
     dvxc3=(vxcdn1(ir,1)+vxcdn2(ir,1)-two*vxci(ir,1))*fact !d2Vxc_up/dn_dn2
     dvxc4=(vxcdn1(ir,2)+vxcdn2(ir,2)-two*vxci(ir,2))*fact !d2Vxc_dn/dn_dn2
     dvxca=dvxc1+dvxc4;dvxcb=dvxc1-dvxc4 !Temporary terms
     dvxcc=dvxc2+dvxc3;dvxcd=dvxc2-dvxc3 !Temporary terms
     d2vxc(ir,1)=eighth*(dvxca+three*dvxcc)  ! 1/2 d(2Vxc_up+Vxc_dn)/dn2
     d2vxc(ir,2)=eighth*(dvxcb-three*dvxcd)  ! 1/2 d2(Vxc_up-Vxc_dn)/dn2
     d2vxc(ir,3)=eighth*(dvxcb+dvxcd)        ! 1/2 d2(Vxc_up+Vxc_dn)/dm2
     d2vxc(ir,4)=eighth*(dvxca-dvxcc)        ! 1/2 d2(Vxc_up-Vxc_dn)/dm2
    end do
   end if
  end if

! If non-collinear magnetism, store 1/2(Vxc_up+Vxc_dn) and 1/2(Vxc_up-Vxc_dn)
  if (nspden==4) then
   vxci(:,1)=half*(vxci(:,1)+vxci(:,2))
   vxci(:,2)=vxci(:,1)-vxci(:,2)
  end if

  deallocate(rhoinv,vxc1,vxc2);if (nspden_updn==2) deallocate(vxcdn1,vxcdn2)

 end if ! option/=4

 deallocate(rhosph)

!----------------------------------------------------------------------
!----- Compute useful sums of densities
!----------------------------------------------------------------------

 if (option<3.or.option/=1) then

! Non-collinear magnetism: V0SUM=(m_0.m_L)/|m_0|
! --------------------------------------------------
  if (nspden==4) then
   allocate(v0sum(nrad,lm_size));v0sum(:,1)=zero
   do ilm=2,lm_size
    v0sum(1:nrad,ilm)=(rho_updn(1:nrad,1,2)*rho_updn(1:nrad,ilm,2) &
&    +rho_updn(1:nrad,1,3)*rho_updn(1:nrad,ilm,3) &
&    +rho_updn(1:nrad,1,4)*rho_updn(1:nrad,ilm,4))*m_norm_inv(1:nrad)
   end do
  end if

! V1SUM(r,sig1,sig2)= Sum_L{n^sig1_L(r)*n^sig2_L(r)}
! where (sig1,sig2) are (up,dn) or (n,m)
! --------------------------------------------------
  if (pawxcdev>=1) then
   allocate(v1sum(nrad,2*nspden_updn-1));v1sum=zero

!  Non-magnetic system: compute V1SUM1(r)=Sum_L{n_L(r)^2}
   if (nspden==1) then
    do ilm=2,lm_size
     if (lmselect(ilm)) then
      v1sum(1:nrad,1)=v1sum(1:nrad,1)+rho_updn(1:nrad,ilm,1)**2
     end if
    end do
!   Collinear magnetism: compute V1SUM1(r)=Sum_L{n^up_L(r)^2}
!   V1SUM2(r)=Sum_L{n^up_L(r)*n^dn_L(r)}
!   V1SUM3(r)=Sum_L{n^dn_L(r)^2}
   else if (nspden==2) then
    do ilm=2,lm_size
     if (lmselect(ilm)) then
      v1sum(1:nrad,1)=v1sum(1:nrad,1)+rho_updn(1:nrad,ilm,1)**2
      v1sum(1:nrad,2)=v1sum(1:nrad,2)+rho_updn(1:nrad,ilm,1)*rho_updn(1:nrad,ilm,2)
      v1sum(1:nrad,3)=v1sum(1:nrad,3)+rho_updn(1:nrad,ilm,2)**2
     end if
    end do
!   Non-collinear magnetism: compute V1SUM1(r)=Sum_L{n_L(r)^2}
!   V1SUM2(r)=Sum_L{n_L(r) (m_0.m_L)}/|m_0|
!   V1SUM3(r)=Sum_L{(m_0.m_L)^2}/|m_0|^2
   else if (nspden==4) then
    do ilm=2,lm_size
     if (lmselect(ilm)) then
      v1sum(1:nrad,1)=v1sum(1:nrad,1)+rho_updn(1:nrad,ilm,1)**2
      v1sum(1:nrad,2)=v1sum(1:nrad,2)+rho_updn(1:nrad,ilm,1)*v0sum(1:nrad,ilm)
      v1sum(1:nrad,3)=v1sum(1:nrad,3)+v0sum(1:nrad,ilm)**2
     end if
    end do
   end if
  end if !pawxcdev

! V2SUM(r,L,sig1,sig2)= Sum_L1_L2{n^sig1_L1(r)*n^sig2_L2(r)*Gaunt(L,L1,L2)}
! where (sig1,sig2) are (up,dn) or (n,m)
! --------------------------------------------------
  if (pawxcdev>=2) then
   allocate(v2sum(nrad,lm_size,2*nspden_updn-1));v2sum=zero

!  Non-magnetic system: compute V2SUM1(r,L)=Sum_L1_L2{n_L1(r)*n_L2(r)*Gaunt(L,L1,L2)}
   if (nspden==1) then
    do ilm=2,lm_size
     do ilm1=2,lm_size
      if (lmselect(ilm1)) then
       do ilm2=2,ilm1
        if (lmselect(ilm2)) then
         isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
         if (isel>0) then
          fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
          v2sum(1:nrad,ilm,1)=v2sum(1:nrad,ilm,1)+fact*rho_updn(1:nrad,ilm1,1)*rho_updn(1:nrad,ilm2,1)
         end if
        end if
       end do
      end if
     end do
    end do
!   Collinear magnetism: compute V2SUM1(r,L)=Sum_L1_L2{n^up_L1(r)*n^up_L2(r)*Gaunt(L,L1,L2)}
!   V2SUM2(r,L)=Sum_L1_L2{n^up_L1(r)*n^dn_L2(r)*Gaunt(L,L1,L2)}
!   V2SUM3(r,L)=Sum_L1_L2{n^dn_L1(r)*n^dn_L2(r)*Gaunt(L,L1,L2)}
   else if (nspden==2) then
    do ilm=2,lm_size
     do ilm1=2,lm_size
      if (lmselect(ilm1)) then
       do ilm2=2,ilm1
        if (lmselect(ilm2)) then
         isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
         if (isel>0) then
          fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
          v2sum(1:nrad,ilm,1)=v2sum(1:nrad,ilm,1)+fact*rho_updn(1:nrad,ilm1,1)*rho_updn(1:nrad,ilm2,1)
          v2sum(1:nrad,ilm,2)=v2sum(1:nrad,ilm,2)+fact*rho_updn(1:nrad,ilm1,1)*rho_updn(1:nrad,ilm2,2)
          v2sum(1:nrad,ilm,3)=v2sum(1:nrad,ilm,3)+fact*rho_updn(1:nrad,ilm1,2)*rho_updn(1:nrad,ilm2,2)
         end if
        end if
       end do
      end if
     end do
    end do
!   Non-collinear magnetism: compute V2SUM1(r,L)=Sum_L1_L2{n_L1(r)*n_L2(r)*Gaunt(L,L1,L2)}
!   V2SUM2(r,L)=Sum_L1_L2{n_L1(r) (m_0.m_L2)*Gaunt(L,L1,L2)}/|m_0|
!   V2SUM3(r,L)=Sum_L1_L2{(m_0.m_L1)*(m_0.m_L2)*Gaunt(L,L1,L2)}/|m_0|^2
   else if (nspden==4) then
    do ilm=2,lm_size
     do ilm1=2,lm_size
      if (lmselect(ilm1)) then
       do ilm2=2,ilm1
        if (lmselect(ilm2)) then
         isel=pawang%gntselect(ilm,ilm2+ilm1*(ilm1-1)/2)
         if (isel>0) then
          fact=pawang%realgnt(isel);if (ilm1/=ilm2) fact=two*fact
          v2sum(1:nrad,ilm,1)=v2sum(1:nrad,ilm,1)+fact*rho_updn(1:nrad,ilm1,1)*rho_updn(1:nrad,ilm2,1)
          v2sum(1:nrad,ilm,2)=v2sum(1:nrad,ilm,2)+fact*rho_updn(1:nrad,ilm1,1)*v0sum(1:nrad,ilm2)
          v2sum(1:nrad,ilm,3)=v2sum(1:nrad,ilm,3)+fact*v0sum(1:nrad,ilm1)     *v0sum(1:nrad,ilm2)
         end if
        end if
       end do
      end if
     end do
    end do
   end if
  end if !pawxcdev

 end if !option

!----------------------------------------------------------------------
!----- Accumulate and store XC potential
!----------------------------------------------------------------------

 if (option<3) then

! === First order development
! ---------------------------
  if (pawxcdev>=1) then

!  Non-magnetic system:
   if (nspden_updn==1) then
    vxc(1:nrad,1,1)=v1sum(1:nrad,1)*d2vxc(1:nrad,1)*invsqfpi2+vxci(1:nrad,1)*sqfpi
    do ilm=2,lm_size
     if (lmselect(ilm)) then
      vxc(1:nrad,ilm,1)=d1vxc(1:nrad,1)*rho_updn(1:nrad,ilm,1)
     end if
    end do

!   Magnetic system:
   else if (nspden_updn==2) then
    vxc(1:nrad,1,1)=vxci(1:nrad,1)*sqfpi+invsqfpi2*(v1sum(1:nrad,1)*d2vxc(1:nrad,1) &
&    +two*v1sum(1:nrad,2)*d2vxc(1:nrad,2)+v1sum(1:nrad,3)*d2vxc(1:nrad,3))
    vxc(1:nrad,1,2)=vxci(1:nrad,2)*sqfpi+invsqfpi2*(v1sum(1:nrad,1)*d2vxc(1:nrad,2) &
&    +two*v1sum(1:nrad,2)*d2vxc(1:nrad,3)+v1sum(1:nrad,3)*d2vxc(1:nrad,4))
    if (nspden==2) then
     do ilm=2,lm_size
      if (lmselect(ilm)) then
       vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,2)+d1vxc(1:nrad,1)*rho_updn(1:nrad,ilm,1) &
&       +d1vxc(1:nrad,2)*rho_updn(1:nrad,ilm,2)
       vxc(1:nrad,ilm,2)=vxc(1:nrad,ilm,2)+d1vxc(1:nrad,2)*rho_updn(1:nrad,ilm,1) &
&       +d1vxc(1:nrad,3)*rho_updn(1:nrad,ilm,2)
      end if
     end do
    else if (nspden==4) then
     do ilm=2,lm_size
      if (lmselect(ilm)) then
       vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,2)+d1vxc(1:nrad,1)*rho_updn(1:nrad,ilm,1) &
&       +d1vxc(1:nrad,2)*v0sum(1:nrad,ilm)
       vxc(1:nrad,ilm,2)=vxc(1:nrad,ilm,2)+d1vxc(1:nrad,2)*rho_updn(1:nrad,ilm,1) &
&       +d1vxc(1:nrad,3)*v0sum(1:nrad,ilm)
      end if
     end do
    end if
   end if
  end if ! pawxcdev>=1

! == 2nd order development
! ---------------------------
  if (pawxcdev>=2) then

!  Non-magnetic system:
   if (nspden_updn==1) then
    do ilm=2,lm_size
     vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1)+half*d2vxc(1:nrad,1)*v2sum(1:nrad,ilm,1)
    end do

!   Magnetic system:
   else if (nspden_updn==2) then
    do ilm=2,lm_size
     vxc(1:nrad,ilm,1)=vxc(1:nrad,ilm,1)+d2vxc(1:nrad,2)*v2sum(1:nrad,ilm,2) &
&     +half*(d2vxc(1:nrad,1)*v2sum(1:nrad,ilm,1)+d2vxc(1:nrad,3)*v2sum(1:nrad,ilm,3))
     vxc(1:nrad,ilm,2)=vxc(1:nrad,ilm,2)+d2vxc(1:nrad,3)*v2sum(1:nrad,ilm,2) &
&     +half*(d2vxc(1:nrad,2)*v2sum(1:nrad,ilm,1)+d2vxc(1:nrad,4)*v2sum(1:nrad,ilm,3))
    end do
   end if
  end if !pawxcdev=2

! === Pathological case: if rho(r) is negative, interpolate Vxc
! -------------------------------------------------------------
  if (lmselect(1)) then
   do ispden=1,nspden_updn
    ir1=0;ir2=0
    do ir=1,nrad
     if (rho_updn(ir,1,ispden)<tol14) then
      if (ir1==0) ir1=ir-1
      ir2=ir+1
     else if (ir1>0) then
      if (ir1>1.or.ir2<nrad) then
       fact=(vxc(ir2,1,ispden)-vxc(ir1,1,ispden))/(pawrad%rad(ir2)-pawrad%rad(ir1))
       do jr=ir1+1,ir2-1
        vxc(jr,1,ispden)=vxc(ir1,1,ispden)+fact*(pawrad%rad(jr)-pawrad%rad(ir1))
       end do
      end if
      ir1=0;ir2=0
     end if
    end do
   end do
  end if

! === Non-collinear magnetism: "rotate" back the XC potential
! ------- ---------------------------------------------------
  if (nspden==4) then
   do ilm=1,lm_size
    do ir=1,nrad
     dvxca=vxc(ir,ilm,1)
     fact=vxc(ir,ilm,2)*m_norm_inv(ir)
     vxc(ir,ilm,1)=dvxca+fact*rho_updn(ir,1,4)
     vxc(ir,ilm,2)=dvxca-fact*rho_updn(ir,1,4)
     vxc(ir,ilm,3)=      fact*rho_updn(ir,1,2)
     vxc(ir,ilm,4)=     -fact*rho_updn(ir,1,3)
    end do
   end do
  end if

 end if !option<3

 if (nspden==4) deallocate(m_norm_inv)

!----------------------------------------------------------------------
!----- Accumulate and store XC energies
!----------------------------------------------------------------------

!----- Calculate Exc (direct scheme) term
!----------------------------------------
 if (option/=1) then
  allocate(ff(nrad))

! Contribution from spherical part of rho
  if (nspden==1.or.nspden==4) then
   ff(1:nrad)=rho_updn(1:nrad,1,1)*exci(1:nrad)*sqfpi
  else if (nspden==2) then
   ff(1:nrad)=(rho_updn(1:nrad,1,1)+rho_updn(1:nrad,1,2))*exci(1:nrad)*sqfpi
  end if

! Contribution from aspherical part of rho
  if (option/=4) then

!  First order development
   if (pawxcdev>=1) then
    if (nspden_updn==1) then
     ff(1:nrad)=ff(1:nrad)+half*v1sum(1:nrad,1)*d1vxc(1:nrad,1)
    else if (nspden_updn==2) then
     ff(1:nrad)=ff(1:nrad)+v1sum(1:nrad,2)*d1vxc(1:nrad,2) &
&     +half*(v1sum(1:nrad,1)*d1vxc(1:nrad,1)+v1sum(1:nrad,3)*d1vxc(1:nrad,3))
    end if
   end if

!  Second order development
   if (pawxcdev>=2) then
    allocate(gg(nrad))

    gg=zero
    do ilm=2,lm_size
     if (lmselect(ilm)) then
      gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,1)*rho_updn(1:nrad,ilm,1)
     end if
    end do
    ff(1:nrad)=ff(1:nrad)+sixth*gg(1:nrad)*d2vxc(1:nrad,1)

    if (nspden_updn==2) then
     gg=zero
     if (nspden==2) then
      do ilm=2,lm_size
       if (lmselect(ilm)) then
        gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,3)*rho_updn(1:nrad,ilm,2)
       end if
      end do
     else if (nspden==4) then
      do ilm=2,lm_size
       if (lmselect(ilm)) then
        gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,3)*v0sum(1:nrad,ilm)
       end if
      end do
     end if
     ff(1:nrad)=ff(1:nrad)+sixth*gg(1:nrad)*d2vxc(1:nrad,4)
     gg=zero
     do ilm=2,lm_size
      if (lmselect(ilm)) then
       gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,2)*rho_updn(1:nrad,ilm,1)
      end if
     end do
     ff(1:nrad)=ff(1:nrad)+half*gg(1:nrad)*d2vxc(1:nrad,2)
     gg=zero
     do ilm=2,lm_size
      if (lmselect(ilm)) then
       gg(1:nrad)=gg(1:nrad)+v2sum(1:nrad,ilm,3)*rho_updn(1:nrad,ilm,1)
      end if
     end do
     ff(1:nrad)=ff(1:nrad)+half*gg(1:nrad)*d2vxc(1:nrad,3)

     deallocate(gg)
    end if
   end if

  end if ! option/=4

  ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
  call simp_gen(enxc,ff,pawrad)
  deallocate(ff)
 end if ! option/=1

 deallocate(exci,vxci)
 if (nspden==4  .and.(option<3.or.option/=1)) deallocate(v0sum)
 if (pawxcdev>=1.and.(option<3.or.option/=1)) deallocate(v1sum)
 if (pawxcdev>=2.and.(option<3.or.option/=1)) deallocate(v2sum)
 if (option<3.or.(option/=4.and.pawxcdev>1)) deallocate(d2vxc)
 if (option/=4) deallocate(d1vxc)

!----- Calculate Excdc double counting term
!------------------------------------------
 if (option==0.or.option==2) then

! Build appropriate density
  if (usexcnhat==1) then
   if (nspden==1.or.nspden==4) then
    rho_updn(:,:,:)=rho_updn(:,:,:)+nhat(:,:,:)
   else if (nspden==2) then
    rho_updn(:,:,1)=rho_updn(:,:,1)+nhat(:,:,2)
    rho_updn(:,:,2)=rho_updn(:,:,2)+nhat(:,:,1)-nhat(:,:,2)
   end if
  end if
  if (usecore==1) then
   if (nspden==1.or.nspden==4) then
    rho_updn(:,1,1)=rho_updn(:,1,1)-sqfpi*corexc(:)
   else if (nspden==2) then
    rho_updn(:,1,1)=rho_updn(:,1,1)-sqfpi2*corexc(:)
    rho_updn(:,1,2)=rho_updn(:,1,2)-sqfpi2*corexc(:)
   end if
  end if

  allocate(ff(nrad));ff(1:nrad)=zero

! Non magnetic or collinear magnetic system:
  if (nspden/=4) then
   do ispden=1,nspden_updn
    do ilm=1,lm_size
     if (lmselect(ilm)) ff(1:nrad)=ff(1:nrad)+vxc(1:nrad,ilm,ispden)*rho_updn(1:nrad,ilm,ispden)
    end do
   end do
  else
!  Non-collinear magnetic system:
   do ilm=1,lm_size
    if (lmselect(ilm)) then
     do ir=1,nrad
      dvxca=vxc(ir,ilm,1)+vxc(ir,ilm,2);dvxcb=vxc(ir,ilm,1)-vxc(ir,ilm,2)
      ff(ir)=ff(ir)+half*(dvxca*rho_updn(ir,ilm,1)+dvxcb*rho_updn(ir,ilm,4)) &
&      +vxc(ir,ilm,3)*rho_updn(ir,ilm,2)-vxc(ir,ilm,4)*rho_updn(ir,ilm,3)
     end do
    end if
   end do
  end if

  ff(1:nrad)=ff(1:nrad)*pawrad%rad(1:nrad)**2
  call simp_gen(enxcdc,ff,pawrad)
  deallocate(ff)

 end if ! option

 deallocate(rho_updn)

!----- End of routine
 call timab(81,2,tsec)

!DEBUG
!write(6,*)' pawxcm : exit '
!ENDDEBUG

 end subroutine pawxcm
!!***
