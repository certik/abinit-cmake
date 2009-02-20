!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawxc
!! NAME
!! pawxc
!!
!! FUNCTION
!! PAW only
!! Start from the density or spin-density, and compute xc correlation
!! potential and energies inside a paw sphere.
!! LDA ONLY - USE THE DENSITY OVER A WHOLE SPHERICAL GRID (r,theta,phi)
!! Driver of XC functionals.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!! This routine has been written from rhohxc_coll
!!
!! INPUTS
!!  corexc(pawrad%mesh_size)=core density on radial grid
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
!!    vxc(pawrad%mesh_size,pawang%angl_size,nspden)=xc potential
!!       (spin up in 1st half and spin-down in 2nd half if nspden=2)
!!
!! PARENTS
!!      pawdenpot,psp7in
!!
!! CHILDREN
!!      drivexc,leave_new,mkdenpos,simp_gen,size_dvxc,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawxc(corexc,enxc,enxcdc,ixc,lm_size,lmselect,nhat,nspden,option,&
&                 pawang,pawrad,rhor,usecore,usexcnhat,vxc,xclevel)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13xc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,lm_size,nspden,option,usecore,usexcnhat,xclevel
 real(dp),intent(out) :: enxc,enxcdc
 type(pawang_type),intent(in) :: pawang
 type(pawrad_type),intent(in) :: pawrad
!arrays
 logical,intent(in) :: lmselect(lm_size)
 real(dp),intent(in) :: corexc(pawrad%mesh_size)
 real(dp),intent(in) :: nhat(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(in) :: rhor(pawrad%mesh_size,lm_size,nspden)
 real(dp),intent(out) :: vxc(pawrad%mesh_size,pawang%angl_size,nspden)

!Local variables-------------------------------
!scalars
 integer :: ilm,ipts,ir,ispden,iwarn,ndvxc,ngr2,ngrad,npts,nrad,nspden_updn,nspgrad,nvxcdgr,order
 real(dp) :: dvdn,dvdz,enxcr,factor,vxcrho
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: d2vxcar(:),dnexcdn(:,:),dvxcdgr(:,:),dvxci(:,:)
 real(dp),allocatable :: exci(:),ff(:),grho2_updn(:,:),m_norm(:),rho_updn(:,:)
 real(dp),allocatable :: rhoarr(:,:),rhoarrdc(:,:),rhohat(:,:,:),rhonow(:,:,:)
 real(dp),allocatable :: vxci(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' pawxc : enter with option, nspden ',option,nspden
!ENDDEBUG

 call timab(81,1,tsec)

!----------------------------------------------------------------------
!----- Check options
!----------------------------------------------------------------------

 if(xclevel==2) then
  write(message, '(4a)' ) ch10,&
&  ' pawxc : ERROR -',ch10,&
&  '  GGA is not implemented !'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!----------------------------------------------------------------------
!----- Initializations
!----------------------------------------------------------------------

 iwarn=0
 nspden_updn=min(nspden,2)
 nrad=pawrad%mesh_size
 npts=pawang%angl_size
 ngrad=1 ! ngrad=1 is for LDAs or LSDs
 nspgrad=nspden_updn*ngrad
 if (option/=1) enxc=zero
 if (option==0.or.option==2) enxcdc=zero
 if (option<3) vxc(:,:,:)=zero

 if (ixc==0) then ! No xc at all is applied (usually for testing)
  write(message, '(a,a,a,a)' ) ch10,&
&  ' pawxc : WARNING -',ch10,&
&  '  Note that no xc is applied (ixc=0).'
  call wrtout(06,message,'COLL')

 else if (1<=ixc .and. ixc<=10)  then

  allocate(rhonow(nrad,nspden,ngrad*ngrad),rhoarr(nrad,nspden))
  if ((usexcnhat==1.or.usecore==1).and.(option==0.or.option==2)) allocate(rhoarrdc(nrad,nspden))
  if (usexcnhat>0) then
   allocate(rhohat(nrad,lm_size,nspden))
   rhohat(:,:,:)=rhor(:,:,:)+nhat(:,:,:)
  end if
  if (nspden==4) allocate(m_norm(nrad))

! ----------------------------------------------------------------------
! ----- Loop on the angular part and inits
! ----------------------------------------------------------------------

! Do loop on the angular part
  do ipts=1,npts

!  Copy the input density for this (theta,phi)
   rhoarr(:,:)=zero
   if (usexcnhat==1.and.(option==0.or.option==2)) rhoarrdc(:,:)=zero
   do ispden=1,nspden
    do ilm=1,lm_size
     if (lmselect(ilm)) then
      if (usexcnhat==0) then
       do ir=1,nrad
        rhoarr(ir,ispden)=rhoarr(ir,ispden)+rhor(ir,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
      else if (usexcnhat==1.and.(option==0.or.option==2)) then
       do ir=1,nrad
        rhoarr  (ir,ispden)=rhoarr  (ir,ispden)+rhor  (ir,ilm,ispden)*pawang%ylmr(ilm,ipts)
        rhoarrdc(ir,ispden)=rhoarrdc(ir,ispden)+rhohat(ir,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
      else
       do ir=1,nrad
        rhoarr(ir,ispden)=rhoarr(ir,ispden)+rhohat(ir,ilm,ispden)*pawang%ylmr(ilm,ipts)
       end do
      end if
     end if
    end do
   end do

   if ((option==0.or.option==2).and.usecore==1.and.usexcnhat/=1) rhoarrdc(:,:)=rhoarr(:,:)
   if (usecore==1) then
    rhoarr(:,1)=rhoarr(:,1)+corexc(:)
    if (nspden==2) rhoarr(:,2)=rhoarr(:,2)+half*corexc(:)
   end if
   rhonow(1:nrad,1:nspden,1)=rhoarr(1:nrad,1:nspden)

!  The variable order indicates to which derivative of the energy
!  the computation must be done. Here, no derivative.
   order=1

!  Allocation of mandatory arguments of drivexc
   allocate(exci(nrad),vxci(nrad,nspden_updn),rho_updn(nrad,nspden_updn))
!  Allocation of optional arguments
   call size_dvxc(ixc,ndvxc,ngr2,nspden_updn,nvxcdgr,order)
   if (ndvxc/=0) allocate(dvxci(nrad,ndvxc))
   if (nvxcdgr/=0) allocate(dvxcdgr(nrad,nvxcdgr))
   if ((ixc==3 .or. xclevel==2) .and. order==3) allocate(d2vxcar(nrad))
   if (ngrad == 2) allocate(grho2_updn(nrad,ngr2))

   if (nspden==1) then
    rho_updn(1:nrad,1)=rhonow(1:nrad,1,1)*half
   else if (nspden==2) then
    rho_updn(1:nrad,1)=rhonow(1:nrad,2,1)
    rho_updn(1:nrad,2)=rhonow(1:nrad,1,1)-rhonow(1:nrad,2,1)
   else if (nspden==4) then
    do ir=1,nrad
     m_norm(ir)=sqrt(rhonow(ir,2,1)**2+rhonow(ir,3,1)**2+rhonow(ir,4,1)**2)
     rho_updn(ir,1)=(rhonow(ir,1,1)+m_norm(ir))*half
     rho_updn(ir,2)=(rhonow(ir,1,1)-m_norm(ir))*half
    end do
   end if

!  Make the density positive everywhere (but do not care about gradients)
   call mkdenpos(iwarn,nrad,nspden_updn,0,rho_updn)

!  Cases with gradient
   if (xclevel==2)then
    if (order**2 <= 1 .or. ixc == 16) then
     if (ixc /= 13) then
      call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,&
&      grho2_updn=grho2_updn,vxcgr=dvxcdgr)
     else
      call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,&
&      grho2_updn=grho2_updn)
     end if
    else if (order /= 3) then
     if (ixc /= 13) then
      call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,&
&      dvxc=dvxci,grho2_updn=grho2_updn,vxcgr=dvxcdgr)
     else
      call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,&
&      dvxc=dvxci,grho2_updn=grho2_updn)
     end if
    else if (order == 3) then
     if (ixc /= 13) then
      call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,&
&      dvxc=dvxci,d2vxc=d2vxcar,grho2_updn=grho2_updn,vxcgr=dvxcdgr)
     else
      call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,&
&      dvxc=dvxci,d2vxc=d2vxcar,grho2_updn=grho2_updn)
     end if
    end if
!   Cases without gradient
   else
    if (order**2 <=1 .or. ixc >= 31 .and. ixc<=34) then
     call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr)
    else if (order==3 .and. (ixc==3 .or. ixc>=7 .and. ixc<=10)) then
     call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,&
&     dvxc=dvxci,d2vxc=d2vxcar)
    else
     call drivexc(exci,ixc,nrad,nspden_updn,order,rho_updn,vxci,ndvxc,ngr2,nvxcdgr,&
&     dvxc=dvxci)
    end if
   end if

!  ----------------------------------------------------------------------
!  ----- Accumulate and store XC potential
!  ----------------------------------------------------------------------
   if (option<3) then
    if (nspden/=4) then
     do ispden=1,nspden
      do ir=1,nrad
       vxc(ir,ipts,ispden)=vxci(ir,ispden)
      end do
     end do
    else
     do ir=1,nrad
      dvdn=(vxci(ir,1)+vxci(ir,2))*half
      dvdz=(vxci(ir,1)-vxci(ir,2))*half
      if(m_norm(ir)>rhonow(ir,1,1)*tol10+tol14)then
       factor=dvdz/m_norm(ir)
       vxc(ir,ipts,1)=dvdn+rhonow(ir,4,1)*factor
       vxc(ir,ipts,2)=dvdn-rhonow(ir,4,1)*factor
       vxc(ir,ipts,3)= rhonow(ir,2,1)*factor
       vxc(ir,ipts,4)=-rhonow(ir,3,1)*factor
      else
       vxc(ir,ipts,1:2)=dvdn
       vxc(ir,ipts,3:4)=zero
      end if
     end do
    end if
   end if !option

!  ----------------------------------------------------------------------
!  ----- Accumulate and store XC energies
!  ----------------------------------------------------------------------

   if (option/=1) allocate(ff(nrad))

!  ----- Calculate Exc term
   if (option/=1) then
    ff(:)=rhoarr(:,1)*exci(:)*pawrad%rad(:)**2
    call simp_gen(enxcr,ff,pawrad)
    if (option/=4) enxc=enxc+enxcr*pawang%angwgth(ipts)
    if (option==4) enxc=enxc+enxcr
   end if

!  ----- Calculate Excdc double counting term
   if (option==0.or.option==2) then
    if (nspden/=4) then
     if (usexcnhat==1.or.usecore==1) then
      ff(:)=vxc(:,ipts,1)*rhoarrdc(:,nspden)
      if (nspden==2) ff(:)=ff(:)+vxc(:,ipts,2)*(rhoarrdc(:,1)-rhoarrdc(:,2))
     else
      ff(:)=vxc(:,ipts,1)*rhoarr(:,nspden)
      if (nspden==2) ff(:)=ff(:)+vxc(:,ipts,2)*(rhoarr(:,1)-rhoarr(:,2))
     end if
    else
     if (usexcnhat==1.or.usecore==1) then
      ff(:)=half*(vxc(:,ipts,1)*(rhoarrdc(:,1)+rhoarrdc(:,4))+vxc(:,ipts,2)*(rhoarrdc(:,1)-rhoarrdc(:,4))) &
&      +vxc(:,ipts,3)*rhoarrdc(:,2)-vxc(:,ipts,4)*rhoarrdc(:,3)
     else
      ff(:)=half*(vxc(:,ipts,1)*(rhoarr(:,1)+rhoarr(:,4))+vxc(:,ipts,2)*(rhoarr(:,1)-rhoarr(:,4))) &
&      +vxc(:,ipts,3)*rhoarr(:,2)-vxc(:,ipts,4)*rhoarr(:,3)
     end if
    end if
    ff(:)=ff(:)*pawrad%rad(:)**2
    call simp_gen(vxcrho,ff,pawrad)
    enxcdc=enxcdc+vxcrho*pawang%angwgth(ipts)
   end if

   if (option/=1) deallocate(ff)
!  ---------------------------------------------------
   deallocate(exci,rho_updn,vxci)
   if (allocated(dvxci)) deallocate(dvxci)
   if (allocated(dvxcdgr)) deallocate(dvxcdgr)
   if (allocated(d2vxcar)) deallocate(d2vxcar)
   if (allocated(grho2_updn)) deallocate(grho2_updn)

!  ----- End of the loop on npts (angular part)
  end do

! Add the four*pi factor of the angular integration
  if (option/=1) enxc=enxc*four_pi
  if (option==0.or.option==2) enxcdc=enxcdc*four_pi

  deallocate(rhonow,rhoarr)
  if ((usexcnhat==1.or.usecore==1).and.(option==0.or.option==2)) deallocate(rhoarrdc)
  if (usexcnhat>0) deallocate(rhohat)
  if (nspden==4) deallocate(m_norm)

! ------------------------------------
! End IF a xc part has to be computed
 end if
 call timab(81,2,tsec)

!DEBUG
!write(6,*)' pawxc : exit '
!ENDDEBUG

 end subroutine pawxc
!!***
