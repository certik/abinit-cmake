!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawgylm
!! NAME
!! pawgylm
!!
!! FUNCTION
!! Compute g_l(r)*Y_lm(r) (and derivatives) on the fine (rectangular) grid
!! around one atom (g_l=radial shape function).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  iatom=index of current atom
!!  ifftsph(nfgd)=FFT index (fine grid) of points in the paw sphere around current atom
!!  itypat=type of current atom
!!  lm_size=number of lm components to be calculated
!!  nfgd= number of (fine grid) FFT points in the paw sphere around current atom
!!  optgr0= 1 if g_l(r)*Y_lm(r) are computed
!!  optgr1= 1 if first derivatives of g_l(r)*Y_lm(r) are computed
!!  optgr2= 1 if second derivatives of g_l(r)*Y_lm(r) are computed
!!  pawtab <type(pawtab_type)>=paw tabulated starting data for current atom
!!  rfgd(3,nfgd)= coordinates of r-r_atom on the fine rect. grid around current atom
!!  rfgd_allocated= >0 if rfgd array has been allocated (for testing purpose only)
!!
!! OUTPUT
!!  if (optgr0==1)
!!    gylm(nfgd,lm_size)= g_l(r)*Y_lm(r) around current atom
!!  if (optgr1==1)
!!    gylmgr(3,nfgd,lm_size)= derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!  if (optgr2==1)
!!    gylmgr2(6,nfgd,lm_size)= second derivatives of g_l(r)*Y_lm(r) wrt cart. coordinates
!!
!! PARENTS
!!      nhatgrid,pawdij,pawgrnl,pawmknhat
!!
!! CHILDREN
!!      initylmr,jbessel,leave_new,nderiv_gen,sort_dp,splint,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawgylm(gylm,gylmgr,gylmgr2,iatom,ifftsph,itypat,lm_size,nfgd,optgr0,optgr1,optgr2,&
&                  pawtab,rfgd,rfgd_allocated)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13paw, except_this_one => pawgylm
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: iatom,itypat,lm_size,nfgd,optgr0,optgr1,optgr2
 integer,intent(in) :: rfgd_allocated
 type(pawtab_type),intent(in) :: pawtab
!arrays
 integer,intent(in) :: ifftsph(nfgd)
 real(dp),intent(in) :: rfgd(3,nfgd)
 real(dp),intent(out) :: gylm(nfgd,optgr0*lm_size)
 real(dp),intent(out) :: gylmgr(3,nfgd,optgr1*lm_size)
 real(dp),intent(out) :: gylmgr2(6,nfgd,optgr2*lm_size)

!Local variables ------------------------------
!scalars
 integer :: argl,ic,ilm,izero,l0,l_size,lambda,ll,mm,normchoice,option,shape_class
 integer :: shape_type
 real(dp) :: arg,cc,d2shpfunc1,d2shpfunc2,d2shpfunc3,dshpfunc1,dshpfunc2
 real(dp) :: dshpfunc3,g_tmp,jbes1,jbes2,jbesp1,jbesp2,jbespp1,jbespp2
 real(dp) :: pi_over_rshp,r2,s1,s2,s3,shapefunc1,shapefunc2,shapefunc3,sigma
 logical :: compute_gr0,compute_gr1,compute_gr2
 character(len=500) :: message
!arrays
 integer,allocatable :: isort(:)
 real(dp) :: ss(4)
 real(dp),allocatable :: alpha(:,:),d2gfact(:,:),d2shpfuncnum(:,:),dgfact(:,:)
 real(dp),allocatable :: dshpfunc(:,:),dshpfuncnum(:,:),gfact(:,:)
 real(dp),allocatable :: qq(:,:),rnrm(:),rnrm_inv(:),rnrm_sort(:)
 real(dp),allocatable :: shpfuncnum(:,:),work(:),ylmr(:,:),ylmrgr(:,:,:)
!no_abirules
!Statement functions -----------------------------------
 shapefunc1(arg)= exp(-(arg/sigma)**lambda)
 shapefunc2(arg)= (sin(pi_over_rshp*arg)/(pi_over_rshp*arg))**2
 shapefunc3(jbes1,jbes2,argl)= alpha(1,1+argl)*jbes1+alpha(2,1+argl)*jbes2
 dshpfunc1(arg)=-lambda/sigma*(arg/sigma)**(lambda-1)*exp(-(arg/sigma)**lambda)
 dshpfunc2(arg)=two*pi_over_rshp*sin(pi_over_rshp*arg)/(pi_over_rshp*arg)**3&
&              *(pi_over_rshp*arg*cos(pi_over_rshp*arg)-sin(pi_over_rshp*arg))
 dshpfunc3(jbesp1,jbesp2,argl)= alpha(1,1+argl)*qq(1,1+argl)*jbesp1 &
&                              +alpha(2,1+argl)*qq(2,1+argl)*jbesp2
 d2shpfunc1(arg)=lambda/(sigma**2)*(lambda*(arg/sigma)**(2*lambda-2) &
&                -(lambda-1)*(arg/sigma)**(lambda-2))*exp(-(arg/sigma)**lambda)
 d2shpfunc2(arg)=two/(pi_over_rshp**2*arg**4)* &
&               (pi_over_rshp**2*arg**2*(cos(pi_over_rshp*arg))**2 &
&               +(three-pi_over_rshp**2*arg**2)*(sin(pi_over_rshp*arg))**2 &
&               -four*pi_over_rshp*arg*cos(pi_over_rshp*arg)*sin(pi_over_rshp*arg))
 d2shpfunc3(jbespp1,jbespp2,argl)= alpha(1,1+argl)*(qq(1,1+argl)**2)*jbespp1 &
&                                 +alpha(2,1+argl)*(qq(2,1+argl)**2)*jbespp2

! *************************************************************************

 if (optgr0==0.and.optgr1==0.and.optgr2==0) return
 if (nfgd==0) return

!Compatibility test
!==========================================================
 if (rfgd_allocated==0) then
  write(message, '(a,a,a,a,i3,a)' )ch10,&
&  ' pawgylm : BUG -',ch10,&
&  '  rfgd array must be allocated (atom ',iatom, ') !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Initializations
!==========================================================
!Options for computation
 compute_gr0=(optgr0==1.or.optgr1==1.or.optgr2==1)
 compute_gr1=(optgr1==1.or.optgr2==1)
 compute_gr2=(optgr2==1)
 l_size=pawtab%lcut_size

!Norms of vectors around the atom
 allocate(rnrm(nfgd));izero=-1
 do ic=1,nfgd
  rnrm(ic)=sqrt(rfgd(1,ic)**2+rfgd(2,ic)**2+rfgd(3,ic)**2)
  if (rnrm(ic)<=ten*epsilon(one)) izero=ic
 end do

!Initializations
 if (optgr0==1) gylm=zero
 if (optgr1==1) gylmgr=zero
 if (optgr2==1) gylmgr2=zero

!Some definitions concerning shape function g_l(r)
 shape_type=pawtab%shape_type
 shape_class=0;if (shape_type==1.or.shape_type==2) shape_class=1 !Shape_class=1 if gl(r)=N.k(r).r^l
 if (shape_type==1) then
  sigma=pawtab%shape_sigma;lambda=pawtab%shape_lambda
 end if
 if (shape_type==2) pi_over_rshp=pi/pawtab%rshp
 if (shape_type==3) then
  allocate(alpha(2,l_size),qq(2,l_size))
  do ll=1,l_size
   alpha(1:2,ll)=pawtab%shape_alpha(1:2,ll)
   qq(1:2,ll)=pawtab%shape_q(1:2,ll)
  end do
 end if

!If needed, sort selected radii by increasing norm
 if (shape_type==-1) then
  allocate(isort(nfgd),rnrm_sort(nfgd))
  do ic=1,nfgd
   isort(ic)=ic
  end do
  rnrm_sort(1:nfgd)=rnrm(1:nfgd)
  call sort_dp(nfgd,rnrm_sort,isort,tol16)
 end if

!If shape function is "numeric", spline it onto selected radii
 if (shape_type==-1) then
  allocate(work(nfgd))
  if (compute_gr0) then
   allocate(shpfuncnum(nfgd,l_size))
   do ll=1,l_size
    call splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%shapefunc(:,ll),&
&    pawtab%dshpfunc(:,ll,2),nfgd,rnrm_sort,work)
    do ic=1,nfgd
     shpfuncnum(isort(ic),ll)=work(ic)
    end do
   end do
  end if
  if(compute_gr1) then
   allocate(dshpfuncnum(nfgd,l_size))
   do ll=1,l_size
    call splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%dshpfunc(:,ll,1),&
&    pawtab%dshpfunc(:,ll,3),nfgd,rnrm_sort,work)
    do ic=1,nfgd
     dshpfuncnum(isort(ic),ll)=work(ic)
    end do
   end do
  end if
  if(compute_gr2) then
   allocate(d2shpfuncnum(nfgd,l_size))
   do ll=1,l_size
    call splint(pawtab%mesh_size,pawtab%rad_for_spline,pawtab%dshpfunc(:,ll,2),&
&    pawtab%dshpfunc(:,ll,4),nfgd,rnrm_sort,work)
    do ic=1,nfgd
     d2shpfuncnum(isort(ic),ll)=work(ic)
    end do
   end do
  end if
  deallocate(work)
 end if

 if (shape_type==-1) deallocate(isort,rnrm_sort)

!Y_lm(r) calculation
!==========================================================
 normchoice=1 ; option=max(optgr0,2*optgr1,3*optgr2)
 if(compute_gr0) allocate(ylmr(l_size**2,nfgd))
 if(compute_gr1.and.(.not.compute_gr2)) allocate(ylmrgr(3,l_size**2,nfgd))
 if(compute_gr2) allocate(ylmrgr(9,l_size**2,nfgd))
 call initylmr(l_size,normchoice,nfgd,rnrm,option,rfgd,ylmr,ylmrgr)

!g_l(r) calculation (and factors for derivatives)
!==========================================================
 if (compute_gr0) allocate(gfact(nfgd,0:l_size-1))
 if (compute_gr1) allocate(dgfact(nfgd,0:l_size-1))
 if (compute_gr2) allocate(d2gfact(nfgd,0:l_size-1))
 if(compute_gr1) then
  allocate(rnrm_inv(nfgd))
  do ic=1,nfgd
   if (ic/=izero) rnrm_inv(ic)=1._dp/rnrm(ic)
  end do
  if (izero>0) rnrm_inv(izero)=zero
 end if
!----- type -1 -----
 if (shape_type==-1) then
  if (compute_gr0) then
   do ll=0,l_size-1
    gfact(1:nfgd,ll)=shpfuncnum(1:nfgd,ll+1)
   end do
  end if
  if (compute_gr1) then
   do ll=0,l_size-1
    dgfact(1:nfgd,ll)=dshpfuncnum(1:nfgd,ll+1)*rnrm_inv(1:nfgd)
   end do
  end if
  if(compute_gr2) then
   do ll=0,l_size-1
    d2gfact(1:nfgd,ll)=(d2shpfuncnum(1:nfgd,ll+1)-dgfact(1:nfgd,ll))*rnrm_inv(1:nfgd)**2
   end do
  end if
 end if
!----- type 1 or 2 -----
 if (shape_type==1.or.shape_type==2) then
  if (compute_gr0) then
   if (shape_type==1) then
    do ic=1,nfgd
     if (ic/=izero) gfact(ic,0)=shapefunc1(rnrm(ic))
    end do
    if (izero>0) gfact(izero,0)=one
   else if (shape_type==2) then
    do ic=1,nfgd
     if (ic/=izero) gfact(ic,0)=shapefunc2(rnrm(ic))
    end do
    if (izero>0) gfact(izero,0)=one
   end if
   if (l_size>1) then
    do ll=1,l_size-1
     do ic=1,nfgd
      gfact(ic,ll)=pawtab%gnorm(ll+1)*gfact(ic,0)*(rnrm(ic)**ll)
     end do
    end do
   end if
  end if
  if (compute_gr1) then
   if (shape_type==1) then
    do ic=1,nfgd
     if (ic/=izero) dgfact(ic,0)=dshpfunc1(rnrm(ic))*rnrm_inv(ic)
    end do
    if (izero>0) then
     if (lambda/=2) dgfact(izero,0)=zero  ! This not true for lambda=1
     if (lambda==2) dgfact(izero,0)=-two/lambda**2
    end if
   else if (shape_type==2) then
    do ic=1,nfgd
     if (ic/=izero) dgfact(ic,0)=dshpfunc2(rnrm(ic))*rnrm_inv(ic)
    end do
    if (izero>0) dgfact(izero,0)=-two*pi_over_rshp**2/three
   end if
   if (l_size>1) then
    do ic=1,nfgd
     dgfact(ic,1)=pawtab%gnorm(2)*(dgfact(ic,0)*rnrm(ic)+rnrm_inv(ic)*gfact(ic,0))
    end do
   end if
   if (l_size>2) then
    do ic=1,nfgd
     dgfact(ic,2)=pawtab%gnorm(3)*(dgfact(ic,0)*rnrm(ic)**2+two*gfact(ic,0))
    end do
   end if
   if (l_size>3) then
    do ll=3,l_size-1
     do ic=1,nfgd
      dgfact(ic,ll)=pawtab%gnorm(ll+1)*(dgfact(ic,0)*rnrm(ic)**ll+ll*gfact(ic,0)*rnrm(ic)**(ll-2))
     end do
    end do
   end if
  end if
  if (compute_gr2) then
   if (shape_type==1) then
    do ic=1,nfgd
     if (ic/=izero) d2gfact(ic,0)=rnrm_inv(ic)**2*(d2shpfunc1(rnrm(ic))-dgfact(ic,0))
    end do
    if (izero>0) then
     if (lambda/=2.and.lambda/=4) d2gfact(izero,0)=zero  ! This not true for lambda=1 or 3
     if (lambda==2) d2gfact(izero,0)=four/lambda**4
     if (lambda==4) d2gfact(izero,0)=-8._dp/lambda**4
    end if
   else if (shape_type==2) then
    do ic=1,nfgd
     if (ic/=izero) d2gfact(ic,0)=rnrm_inv(ic)**2*(d2shpfunc2(rnrm(ic))-dgfact(ic,0))
    end do
    if (izero>0) d2gfact(izero,0)=16._dp*pi_over_rshp**4/45._dp
   end if
   if (l_size>1) then
    do ic=1,nfgd
     d2gfact(ic,1)=pawtab%gnorm(2)*rnrm_inv(ic)**3*(rnrm(ic)**4*d2gfact(ic,0) &
&     +two*rnrm(ic)**2*dgfact(ic,0)-gfact(ic,0))
    end do
   end if
   if (l_size>2) then
    do ic=1,nfgd
     dgfact(ic,2)=pawtab%gnorm(3)*(rnrm(ic)**2*d2gfact(ic,0)+four*dgfact(ic,0))
    end do
   end if
   if (l_size>3) then
    do ll=3,l_size-1
     do ic=1,nfgd
      d2gfact(ic,ll)=pawtab%gnorm(ll+1)*rnrm_inv(ic)**2*(rnrm(ic)**(ll+2)*d2gfact(ic,0) &
&      +two*ll*rnrm(ic)**ll*dgfact(ic,0)+ll*(ll-2)*rnrm(ic)**(ll-2)*gfact(ic,0))
     end do
    end do
   end if
  end if

  if (compute_gr0) gfact(:,0)=gfact(:,0)*pawtab%gnorm(1)
  if (compute_gr1) dgfact(:,0)=dgfact(:,0)*pawtab%gnorm(1)
  if (compute_gr2) d2gfact(:,0)=d2gfact(:,0)*pawtab%gnorm(1)

! ----- type 3 -----
 else if (shape_type==3) then
  if (optgr0==1.and.optgr1==0.and.optgr2==0) then
   do ll=0,l_size-1
    do ic=1,nfgd
     call jbessel(jbes1,jbesp1,jbespp1,ll,0,qq(1,1+ll)*rnrm(ic))
     call jbessel(jbes2,jbesp2,jbespp2,ll,0,qq(2,1+ll)*rnrm(ic))
     gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
    end do
   end do
  else if (optgr1==1.and.optgr2==0) then
   do ll=0,l_size-1
    do ic=1,nfgd
     call jbessel(jbes1,jbesp1,jbespp1,ll,1,qq(1,1+ll)*rnrm(ic))
     call jbessel(jbes2,jbesp2,jbespp2,ll,1,qq(2,1+ll)*rnrm(ic))
     gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
     dgfact(ic,ll)=dshpfunc3(jbesp1,jbesp2,ll)*rnrm_inv(ic)
    end do
   end do
  else if (optgr2==1) then
   do ll=0,l_size-1
    do ic=1,nfgd
     call jbessel(jbes1,jbesp1,jbespp1,ll,2,qq(1,1+ll)*rnrm(ic))
     call jbessel(jbes2,jbesp2,jbespp2,ll,2,qq(2,1+ll)*rnrm(ic))
     gfact(ic,ll)=shapefunc3(jbes1,jbes2,ll)
     dgfact(ic,ll)=dshpfunc3(jbesp1,jbesp2,ll)*rnrm_inv(ic)
     d2gfact(ic,ll)=(d2shpfunc3(jbespp1,jbespp2,ll)-dgfact(ic,ll))*rnrm_inv(ic)**2
    end do
   end do
  end if
 end if

!g_l(r)*Y_lm(r) calculation
!==========================================================
 if (optgr0==1) then

  do ll=0,l_size-1
   l0=ll**2+ll+1
   do mm=-ll,ll
    ilm=l0+mm
    if (ilm<=lm_size) then
     do ic=1,nfgd
      if (ic==izero.and.shape_class==1) cycle
      gylm(ic,ilm)=gfact(ic,ll)*ylmr(ilm,ic)
     end do
    end if
   end do
  end do

! Special value at r=0  (supposing shapefunc(r)->C.r**l when r->0)
  if (izero>0.and.shape_class==1) gylm(izero,1)=ylmr(1,izero)*pawtab%gnorm(1)

 end if

!d/dr{g_l(r)*Y_lm(r)} calculation
!==========================================================
 if(optgr1==1) then

  do ll=0,l_size-1
   l0=ll**2+ll+1
   do mm=-ll,ll
    ilm=l0+mm
    if (ilm<=lm_size) then
     do ic=1,nfgd
      if (ic==izero) cycle
      gylmgr(1:3,ic,ilm)=gfact(ic,ll)*ylmrgr(1:3,ilm,ic)&
&      +dgfact(ic,ll)*rfgd(1:3,ic)*ylmr(ilm,ic)
     end do
    end if
   end do
  end do

! Special values at r=0  (supposing shapefunc(r)->C.r**l when r->0)
  if (izero>0.and.l_size>1) then
   if (shape_type==-1) then
    ss(2)=pawtab%shapefunc(2,2)/pawtab%rad_for_spline(2)
    ss(3)=pawtab%shapefunc(3,2)/pawtab%rad_for_spline(3)
    ss(4)=pawtab%shapefunc(4,2)/pawtab%rad_for_spline(4)
    ss(1)=ss(4)+(ss(2)-ss(3))*(pawtab%rad_for_spline(4)-pawtab%rad_for_spline(1))&
&    /(pawtab%rad_for_spline(3)-pawtab%rad_for_spline(2))
    cc=ss(1)*pawtab%gnorm(2)
   else if (shape_type==1.or.shape_type==2) then
    cc=pawtab%gnorm(2)
   else if (shape_type==3) then
    cc=(alpha(1,2)*qq(1,2)+alpha(2,2)*qq(2,2))/three
   end if
   g_tmp=cc*sqrt(three/four_pi)
   if (lm_size>=2) gylmgr(2,izero,2)=g_tmp
   if (lm_size>=3) gylmgr(3,izero,3)=g_tmp
   if (lm_size>=4) gylmgr(1,izero,4)=g_tmp
  end if

 end if

!d2/dridrj{g_l(r)*Y_lm(r)} calculation
!==========================================================
 if(optgr2==1) then

  do ll=0,l_size-1
   l0=ll**2+ll+1
   do mm=-ll,ll
    ilm=l0+mm
    if (ilm<=lm_size) then
     do ic=1,nfgd
      if (ic==izero) cycle
      gylmgr2(1,ic,ilm)=gfact(ic,ll)*ylmrgr(4,ilm,ic) &
&      +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(1,ic)*ylmrgr(1,ilm,ic)) &
&      +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(1,ic)*rfgd(1,ic)
      gylmgr2(2,ic,ilm)=gfact(ic,ll)*ylmrgr(5,ilm,ic) &
&      +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(2,ic)*ylmrgr(2,ilm,ic)) &
&      +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(2,ic)*rfgd(2,ic)
      gylmgr2(3,ic,ilm)=gfact(ic,ll)*ylmrgr(6,ilm,ic) &
&      +dgfact(ic,ll)*(ylmr(ilm,ic)+two*rfgd(3,ic)*ylmrgr(3,ilm,ic)) &
&      +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(3,ic)
      gylmgr2(4,ic,ilm)=gfact(ic,ll)*ylmrgr(7,ilm,ic) &
&      +dgfact(ic,ll)*(rfgd(3,ic)*ylmrgr(2,ilm,ic)+rfgd(2,ic)*ylmrgr(3,ilm,ic)) &
&      +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(2,ic)
      gylmgr2(5,ic,ilm)=gfact(ic,ll)*ylmrgr(8,ilm,ic) &
&      +dgfact(ic,ll)*(rfgd(3,ic)*ylmrgr(1,ilm,ic)+rfgd(1,ic)*ylmrgr(3,ilm,ic)) &
&      +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(3,ic)*rfgd(1,ic)
      gylmgr2(6,ic,ilm)=gfact(ic,ll)*ylmrgr(9,ilm,ic) &
&      +dgfact(ic,ll)*(rfgd(1,ic)*ylmrgr(2,ilm,ic)+rfgd(2,ic)*ylmrgr(1,ilm,ic)) &
&      +d2gfact(ic,ll)*ylmr(ilm,ic)*rfgd(1,ic)*rfgd(2,ic)
     end do
    end if
   end do
  end do

! Special values at r=0  (supposing shapefunc(r)->C.r**l when r->0)
  if (izero>0.and.l_size>2) then
   if (shape_type==-1) then
    ss(2)=pawtab%shapefunc(2,2)/pawtab%rad_for_spline(2)**2
    ss(3)=pawtab%shapefunc(3,2)/pawtab%rad_for_spline(3)**2
    ss(4)=pawtab%shapefunc(4,2)/pawtab%rad_for_spline(4)**2
    ss(1)=ss(4)+(ss(2)-ss(3))*(pawtab%rad_for_spline(4)-pawtab%rad_for_spline(1))&
&    /(pawtab%rad_for_spline(3)-pawtab%rad_for_spline(2))
    cc=ss(1)*pawtab%gnorm(3)
   else if (shape_type==1.or.shape_type==2) then
    cc=pawtab%gnorm(3)
   else if (shape_type==3) then
    cc=(alpha(1,2)*qq(1,2)**2+alpha(2,2)*qq(2,2)**2)/15._dp
   end if
   g_tmp=cc*sqrt(5._dp/four_pi)
   if (lm_size>=5) gylmgr2(6,izero,5)=sqrt3*g_tmp
   if (lm_size>=6) gylmgr2(4,izero,6)=sqrt3*g_tmp
   if (lm_size>=7) then
    gylmgr2(1,izero,7)=-g_tmp
    gylmgr2(2,izero,7)=-g_tmp
    gylmgr2(3,izero,7)=two*g_tmp
   end if
   if (lm_size>=8) gylmgr2(5,izero,8)=sqrt3*g_tmp
   if (lm_size>=9) then
    gylmgr2(1,izero,9)=sqrt3*g_tmp
    gylmgr2(2,izero,9)=-sqrt3*g_tmp
   end if
  end if

 end if

!Memory deallocation
!==========================================================
 deallocate(rnrm)
 if (compute_gr0) deallocate(gfact)
 if (compute_gr1) deallocate(dgfact)
 if (compute_gr2) deallocate(d2gfact)
 if (compute_gr1) deallocate(rnrm_inv)
 if (shape_type==3) deallocate(alpha,qq)
 if (compute_gr0) deallocate (ylmr)
 if (compute_gr1) deallocate(ylmrgr)
 if (shape_type==-1) then
  if (compute_gr0) deallocate(shpfuncnum)
  if (compute_gr1) deallocate(dshpfuncnum)
  if (compute_gr2) deallocate(d2shpfuncnum)
 end if

end subroutine pawgylm
!!***
