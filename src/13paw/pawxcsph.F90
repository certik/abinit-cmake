!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxcsph
!! NAME
!! pawxcsph
!!
!! FUNCTION
!! PAW only
!! Compute XC energy and potential for a spherical density rho(r) given as (up,dn)
!! Driver of XC functionals. Only treat collinear spins. LDA and GGA
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc_coll
!!
!! INPUTS
!!  exexch= choice of local exact exchange. Active if exexch>0
!!  ixc= choice of exchange-correlation scheme (see above and below)
!!  ndvxc= size of dvxc(npts,ndvxc)
!!  ngr2= size of grho2_updn(npts,ngr2)
!!  ngrad : =1, only compute the density ; =2 also compute the gradient
!!  nrad= dimension of the radial mesh
!!  nspden=number of spin-density components
!!  nspgrad=number of spin-density and spin-density-gradient components
!!  nvxcdgr=second dimension of dvxcdgr
!!  order=gives the maximal derivative of Exc computed
!!  pawrad <type(pawrad_type)>=paw radial mesh and related data
!!  rho_updn(nrad,lm_size,nspden)=electron density in real space
!!             up (ispden=1) and down (ispden=2) parts
!!             If nspden=1, rho_updn(:,:,1) contains (1/2).rho_total
!!  xclevel= XC functional level
!!
!! OUTPUT
!!  exc= XC energy
!!  vxc((nrad,nspden)= XC potential
!!
!! PARENTS
!!      pawxcm
!!
!! CHILDREN
!!      deducer0,drivexc,leave_new,nderiv_gen,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pawxcsph(exc,exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden,&
 &          nspgrad,nvxcdgr,order,pawrad,rho_updn,vxc,xclevel)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13xc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exexch,ixc,ndvxc,ngr2,ngrad,nrad,nspden,nspgrad,nvxcdgr
 integer,intent(in) :: order,xclevel
 type(pawrad_type),intent(in) :: pawrad
!arrays
 real(dp),intent(in) :: rho_updn(nrad,nspden)
 real(dp),intent(out) :: exc(nrad),vxc(nrad,nspden)

!Local variables-------------------------------
!scalars
 integer :: ir,ispden
 real(dp),parameter :: tol24=tol12*tol12
 real(dp) :: coeff,grho_tot,grho_up
 character(len=500) :: message
!arrays
 real(dp),allocatable :: d2vxcar(:),dff(:),dnexcdn(:,:),dvxcdgr(:,:),dvxci(:,:)
 real(dp),allocatable :: grho2(:,:),grho_updn(:,:)

! *************************************************************************

 if(nspden>2)then
  write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&  ' pawxcsph :  BUG -',ch10,&
&  '  Only non-spin-polarised or collinear spin-densities are allowed,',ch10,&
&  '  while the argument nspden=',nspden
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(nrad/=pawrad%mesh_size)then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' pawxcsph :  BUG -',ch10,&
&  '  nrad is not equal to radial mesh size !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!--------------------------------------------------------------------------
!-------------- GGA: computation of the gradient of the density
!--------------------------------------------------------------------------

 if (ngrad==2) then

! grho_updn contains the gradient of the radial part
! grho2(:,1:3) contains the squared norm of this gradient (up, dn and total)
  allocate(grho_updn(nrad,nspden),grho2(nrad,ngr2))

! Gradient of radial part of density
  allocate(dff(nrad))
  do ispden=1,nspden
   call nderiv_gen(dff,rho_updn(:,ispden),1,pawrad)
   grho_updn(:,ispden)=dff(:)
  end do
  deallocate(dff)

! Squared norm of the gradient
  grho2(:,1)=grho_updn(:,1)**2
  if (nspden==2) then
   grho2(:,2)=grho_updn(:,2)**2
   grho2(:,3)=(grho_updn(:,1)+grho_updn(:,2))**2
  end if

 end if

!--------------------------------------------------------------------------
!-------------- Computation of Exc and Vxc
!--------------------------------------------------------------------------

!Used variable
 if (nvxcdgr/=0) allocate(dvxcdgr(nrad,nvxcdgr))
!Not used variables
 if (ndvxc/=0) allocate(dvxci(nrad,ndvxc))
 if ((ixc==3.or.(7<=ixc.and.ixc<=9).or.xclevel==2).and.order==3) allocate(d2vxcar(nrad))

!Cases with gradient
 if (xclevel==2)then
  if (order**2 <= 1 .or. ixc == 16) then
   if (ixc /= 13) then
    call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
&    exexch=exexch,grho2_updn=grho2,vxcgr=dvxcdgr)
   else
    call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
&    grho2_updn=grho2)
   end if
  else if (order /= 3) then
   if (ixc /= 13) then
    call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
&    dvxc=dvxci,grho2_updn=grho2,vxcgr=dvxcdgr)
   else
    call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
&    dvxc=dvxci,grho2_updn=grho2)
   end if
  else if (order == 3) then
   if (ixc /= 13) then
    call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
&    dvxc=dvxci,d2vxc=d2vxcar,grho2_updn=grho2,vxcgr=dvxcdgr)
   else
    call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
&    dvxc=dvxci,d2vxc=d2vxcar,grho2_updn=grho2)
   end if
  end if

! Cases without gradient
 else
  if (order**2 <=1 .or. ixc >= 31 .and. ixc<=34) then
   call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr)
  else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
   call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
&   dvxc=dvxci,d2vxc=d2vxcar)
  else
   call drivexc(exc,ixc,nrad,nspden,order,rho_updn,vxc,ndvxc,ngr2,nvxcdgr,&
&   dvxc=dvxci)
  end if
 end if

 if (ndvxc/=0) deallocate(dvxci)
 if ((ixc==3.or.(7<=ixc.and.ixc<=9).or.xclevel==2).and.order==3) deallocate(d2vxcar)

!--------------------------------------------------------------------------
!-------------- GGA: gardient corrections
!--------------------------------------------------------------------------

 if (ngrad==2) then

! Compute the derivative of Exc with respect to the (spin-)density,
! or to the norm of the gradient of the (spin-)density,
! Further divided by the norm of the gradient of the (spin-)density
! The different components of dnexcdn will be
! for nspden=1,         dnexcdn(:,1)=d(n.exc)/d(n)
! and if ngrad=2, dnexcdn(:,2)=1/2*1/|grad n_up|*d(n.exc)/d(|grad n_up|)
! +   1/|grad n|*d(n.exc)/d(|grad n|)
! (do not forget : |grad n| /= |grad n_up| + |grad n_down|
! for nspden=2,         dnexcdn(:,1)=d(n.exc)/d(n_up)
! dnexcdn(:,2)=d(n.exc)/d(n_down)
! and if ngrad=2, dnexcdn(:,3)=1/|grad n_up|*d(n.exc)/d(|grad n_up|)
! dnexcdn(:,4)=1/|grad n_down|*d(n.exc)/d(|grad n_down|)
! dnexcdn(:,5)=1/|grad n|*d(n.exc)/d(|grad n|)
  allocate(dnexcdn(nrad,nspgrad))
! LDA term
  dnexcdn(:,1:nspden)=vxc(:,1:nspden)
! Additional GGA terms
  do ir=1,nrad
   do ispden=1,3  ! spin_up, spin_down and total spin density
    if (nspden==1.and.ispden>=2) exit
!   If the norm of the gradient vanishes, then the different terms
!   vanishes, but the inverse of the gradient diverges,
!   so skip the update.
    if(grho2(ir,ispden)<tol24) then
     dnexcdn(ir,ispden+nspden)=zero;cycle
    end if
!   Compute the derivative of n.e_xc wrt the spin up, spin down,
!   or total density. In the non-spin-polarized case take the coeff.
!   that will be multiplied by the gradient of the total density.
    if (nvxcdgr/=0) then
     if (nspden==1) then
!     Definition of dvxcdgr changed in v3.3
      if (nvxcdgr==3) then
       coeff=half*dvxcdgr(ir,1)+dvxcdgr(ir,3)
      else
       coeff=half*dvxcdgr(ir,1)
      end if
     else if (nspden==2)then
      if (nvxcdgr==3) then
       coeff=dvxcdgr(ir,ispden)
      else if (ispden/=3) then
       coeff=dvxcdgr(ir,ispden)
      else if (ispden==3) then
       coeff=zero
      end if
     end if
    end if
    dnexcdn(ir,ispden+nspden)=coeff
   end do
  end do

! Calculate grad(rho)*dnexcdn and put it in rho(:,:,2)
  if (nvxcdgr/=0) then
   if(nspden==1)then
    grho_updn(:,1)=grho_updn(:,1)*dnexcdn(:,2)
   else
    do ir=1,nrad
     grho_up=grho_updn(ir,1);grho_tot=grho_up+grho_updn(ir,2)
     grho_updn(ir,1)=grho_up*dnexcdn(ir,3)+grho_tot*dnexcdn(ir,5)
     grho_updn(ir,2)=(grho_tot-grho_up)*dnexcdn(ir,4)+grho_tot*dnexcdn(ir,5)
    end do
   end if
  end if
  deallocate(dnexcdn)

! Compute Vxc
  allocate(dff(nrad))
  if (nspden==1) then
   call nderiv_gen(dff,grho_updn(:,1),1,pawrad)
   vxc(2:nrad,1)=vxc(2:nrad,1)-two*(dff(2:nrad)+two*grho_updn(2:nrad,1)/pawrad%rad(2:nrad))
   call deducer0(vxc(:,1),nrad,pawrad)
  else if (nspden==2) then
   do ispden=1,nspden
    call nderiv_gen(dff,grho_updn(:,ispden),1,pawrad)
    vxc(2:nrad,ispden)=vxc(2:nrad,ispden)-(dff(2:nrad)+two*grho_updn(2:nrad,ispden)/pawrad%rad(2:nrad))
    call deducer0(vxc(:,ispden),nrad,pawrad)
   end do
  end if
  deallocate(dff)

 end if ! ngrad==2

!--------------------------------------------------------------------------
!-------------- Deallocations
!--------------------------------------------------------------------------

 if (ngrad==2) deallocate(grho_updn,grho2)
 if (nvxcdgr/=0) deallocate(dvxcdgr)

 end subroutine pawxcsph
!!***
