!{\src2tex{textfont=tt}}
!!****f* ABINIT/psp6in
!! NAME
!! psp6in
!!
!! FUNCTION
!! Initialize pspcod=6 (Pseudopotentials from the fhi98pp code):
!! continue to read the corresponding file, then compute the
!! local and non-local potentials.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, AF)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  lloc=angular momentum choice of local pseudopotential
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  lnmax=max. number of (l,n) components over all type of psps
!!  mmax=maximum number of points in real space grid in the psp file
!!   angular momentum of nonlocal pseudopotential
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid=dimension of q (or G) grid for arrays.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  optnlxccc=option for nl XC core correction (input variable)
!!  positron=0 if electron GS calculation
!!          1 if positron GS calculation
!!          2 if electron GS calculation in presence of the positron
!!  qgrid(mqgrid)=values of q (or |G|) on grid from 0 to qmax
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  zion=nominal valence of atom as specified in psp file
!!  znucl=nuclear number of atom as specified in psp file
!!
!! OUTPUT
!!  ekb(lnmax)=Kleinman-Bylander energy,
!!             {{\ \begin{equation}
!!               \frac{\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))^2 dr]}
!!             {\int_0^\infty [Rl(r)^2 (Vl(r)-Vloc(r))   dr]}
!!              \end{equation} }}
!!             for each (l,n)
!!  epsatm=$ (4\pi)\int_0^\infty [r^2 (V(r)+\frac{Zv}{r} dr]$ (hartree)
!!  ffspl(mqgrid,2,lnmax)=Kleinman-Bylander form factor f_l(q) and
!!   second derivative from spline fit for each angular momentum and
!!   each projector
!!  indlmn(6,i)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                           or i=lmn (if useylm=1)
!!  nproj(mpsang)=number of projection functions for each angular momentum
!!  qchrg is not used, and could be suppressed later
!!  vlspl(mqgrid,2)=q^2 Vloc(q) and second derivatives from spline fit
!!  xcccrc=XC core correction cutoff radius (bohr)
!!  xccc1d(n1xccc,6)=1D core charge function and five derivatives, from psp file
!!
!! PARENTS
!!      pspatm
!!
!! CHILDREN
!!      calc_psden,psp5lo,psp5nl,psp6cc,psp6cc_drh,psp6ccpos,spline,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine psp6in(ekb,epsatm,ffspl,indlmn,lloc,lmax,lmnmax,lnmax,&
&                  mmax,mpsang,mqgrid,nproj,n1xccc,optnlxccc,positron,qchrg,qgrid,&
&                  useylm,vlspl,xcccrc,xccc1d,zion,znucl)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13psp, except_this_one => psp6in
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lloc,lmax,lmnmax,lnmax,mmax,mpsang,mqgrid,n1xccc
 integer,intent(in) :: optnlxccc,positron,useylm
 real(dp),intent(in) :: zion,znucl
 real(dp),intent(out) :: epsatm,qchrg,xcccrc
!arrays
 integer,intent(out) :: indlmn(6,lmnmax),nproj(mpsang)
 real(dp),intent(in) :: qgrid(mqgrid)
 real(dp),intent(out) :: ekb(lnmax),ffspl(mqgrid,2,lnmax),vlspl(mqgrid,2)
 real(dp),intent(out) :: xccc1d(n1xccc,6)

!Local variables-------------------------------
!scalars
 integer :: ii,index,ipsang,irad,jj,jpsang,mm,mmax2
 real(dp) :: al,al_announced,amesh,fchrg,ratio,rc,rchrg,step,yp1,ypn
 character(len=3) :: testxc
 character(len=500) :: message
!arrays
 real(dp) :: vhtnzc(mmax)
 real(dp),allocatable :: ekb_tmp(:),ffspl_tmp(:,:,:),rad(:),radbis(:),vloc(:)
 real(dp),allocatable :: vpspll(:,:),wfll(:,:),work_space(:),work_spl(:)
 real(dp),allocatable :: xccc1d_wk(:)

! ***************************************************************************

!DEBUG
!write(6,*)' psp6in : enter'
!return
!stop
!ENDDEBUG

!File format of formatted fhi psp input, as adapted for use
!by the ABINIT code (the 3 first lines
!have already been read in calling -pspatm- routine) :

!(1) title (character) line
!(2) znucl,zion,pspdat
!(3) pspcod,pspxc,lmax,lloc,mmax,r2well
!(4) rchrg,fchrg,qchrg
!Note : prior to version 2.2, this 4th line started with  4--  ,
!and no core-correction was available.
!(5)-(18) -empty-
!(19) mmax, amesh ( mesh increment r(m+1)/r(m) )
!Then, for ll=0,lmax :
!for  irad=1,mmax  : irad, r(irad), upsp(irad,ll), vpsp(irad,ll)

 read (tmp_unit, '(a3)') testxc
 if(testxc/='4--')then
  backspace(tmp_unit)
  read (tmp_unit,*) rchrg,fchrg,qchrg
  write(message, '(3f20.14,t64,a)' ) rchrg,fchrg,qchrg,&
&  'rchrg,fchrg,qchrg'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 else
  write(message, '(a)' ) '  No XC core correction.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
  rchrg=0.0d0 ; fchrg=0.0d0 ; qchrg=0.0d0
 end if
 do ii=5,18
  read(tmp_unit,*)
 end do

!--------------------------------------------------------------------
!Will now proceed at the reading of pots and wfs

!rad(:)=radial grid r(i)
!vpspll(:,1),...,vpspll(:,4)=nonlocal pseudopotentials
!wfll(:,1),...,wfll(:,4)=reference config. wavefunctions

 allocate(rad(mmax),vpspll(mmax,mpsang),wfll(mmax,mpsang))

!Read atomic pseudopotential for each l, filling up arrays vpspll
!and wfll. Also set up rad array (actually read more than once)
!Note: put each l into vpspll(:,l+1)
 do ipsang=1,lmax+1
  nproj(ipsang)=1
  read(tmp_unit,*)mmax2,amesh
  if(ipsang==1)then
   write(message, '(f10.6,t20,a)' ) amesh,' amesh (Hamman grid)'
   al_announced=log(amesh)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if
  do irad=1,mmax
   read(tmp_unit,*)jj,rad(irad),wfll(irad,ipsang),vpspll(irad,ipsang)
!  DEBUG
!  Maybe the normalization is different
!  wfll(irad,ipsang)=wfll(irad,ipsang)/rad(irad)
!  ENDDEBUG
  end do
 end do


!Generate core charge function and derivatives, if needed
 if(fchrg>1.0d-15)then
  if (positron == 1) then
   call psp6ccpos(mmax,n1xccc,rchrg,xccc1d,vhtnzc,znucl)
   vpspll(:,1)=vhtnzc(:)
  else
   if(optnlxccc==1)then
    call psp6cc(mmax,n1xccc,rchrg,xccc1d)
   else if(optnlxccc==2)then
    call psp6cc_drh(mmax,n1xccc,rchrg,xccc1d)
   end if
  end if
! The core charge function for pspcod=6
! becomes zero beyond rchrg. Thus xcccrc must be set
! equal to rchrg .
  xcccrc=rchrg

  if (positron /= 0) then
   allocate(xccc1d_wk(n1xccc))
   step = xcccrc/dble(n1xccc-1)
   rc = xcccrc
   call calc_psden(xccc1d_wk,n1xccc,xccc1d(:,1),rc,step)
   xccc1d(:,1)=xccc1d_wk(:)
   deallocate(xccc1d_wk)
  end if
 else
  xccc1d(:,:)=0.0d0
  xcccrc=0.0d0
 end if


!Compute in real(dp) al : the announced amesh is inaccurate.
 ratio=rad(mmax)/rad(1)
 al=log(ratio)/dble(mmax-1)

!DEBUG
!write(6,*)' psp6in : al ; al_announced =',al,al_announced
!allocate(radbis(mmax))
!write(6,*)' psp6in : lloc  ',lloc
!do ipsang=1,lmax+1
!write(6,*)' psp6in : ipsang  ',ipsang
!do irad=1,mmax
!write(6,*)irad,rad(irad),wfll(irad,ipsang),vpspll(irad,ipsang)
!end do
!end do
!deallocate(radbis)
!ENDDEBUG

!vloc(:)=Vlocal(r), lloc=0, 1, or 2 or -1 for avg.
 allocate(vloc(mmax))
!Copy appropriate nonlocal psp for use as local one
 vloc( 1:mmax ) = vpspll( 1:mmax , lloc+1 )

!--------------------------------------------------------------------
!Carry out calculations for local (lloc) pseudopotential.
!Obtain Fourier transform (1-d sine transform)
!to get q^2 V(q).

 call psp5lo(al,epsatm,mmax,mqgrid,qgrid,&
& vlspl(:,1),rad,vloc,yp1,ypn,zion)

!Fit spline to q^2 V(q) (Numerical Recipes subroutine)
 allocate(work_space(mqgrid),work_spl(mqgrid))
 call spline (qgrid,vlspl(:,1),mqgrid,yp1,ypn,work_spl,work_space)
 vlspl(:,2)=work_spl(:)
 deallocate(work_space,work_spl)

!--------------------------------------------------------------------
!Take care of non-local part

 allocate(ekb_tmp(mpsang),ffspl_tmp(mqgrid,2,mpsang))

!Zero out all Kleinman-Bylander energies to initialize
 ekb_tmp(:)=zero
 ekb(:)=zero

!Allow for option of no nonlocal corrections (lloc=lmax=0)
 if (lloc==0.and.lmax==0) then
  write(message, '(a,f5.1)' ) ' Note: local psp for atom with Z=',znucl
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,  message,'COLL')
 else

! ----------------------------------------------------------------------
! Compute KB form factors and fit splines

  call psp5nl(al,ekb_tmp,ffspl_tmp,lmax,mmax,mpsang,mqgrid,qgrid,rad,vloc,&
&  vpspll,wfll)

 end if

 jj=0;index=0;indlmn(:,:)=0
 do ipsang=1,lmax+1
! nproj had been set at 1, by default
  if(abs(ekb_tmp(ipsang))<tol10)then
   nproj(ipsang)=0
  end if
! Possible values for nproj in this routine : 0 or 1.
  if(nproj(ipsang)==1)then
   if (useylm==1) then
    jj=jj+1
    do mm=1,2*ipsang-1
     index=index+1
     indlmn(1,index)=ipsang-1
     indlmn(2,index)=mm-ipsang
     indlmn(3,index)=1
     indlmn(4,index)=mm+(ipsang-1)*(ipsang-1)
     indlmn(5,index)=jj
     indlmn(6,index)=1
    end do
   else
    jj=jj+1
    index=index+1
    indlmn(1,index)=ipsang-1
    indlmn(2,index)=0
    indlmn(3,index)=1
    indlmn(4,index)=ipsang+(ipsang-1)*(ipsang-1)
    indlmn(5,index)=jj
    indlmn(6,index)=1
   end if
  end if
 end do
!Transfer ekb and ffspl to their definitive location
 jpsang=1
 do ipsang=1,lmax+1
  if(nproj(ipsang)/=0)then
   ekb(jpsang)=ekb_tmp(ipsang)
   ffspl(:,:,jpsang)=ffspl_tmp(:,:,ipsang)
   jpsang=jpsang+1
   if(jpsang>lnmax)then
    write(message,'(6a,2i6)') ch10,&
&    ' psp6in : BUG -',ch10,&
&    '  Problem with the dimension of the ekb and ffspl arrays.',ch10,&
&    '  ipsang,lnmax=',ipsang,lnmax
   end if
  end if
 end do

 deallocate(ekb_tmp,ffspl_tmp)

!DEBUG
!write(6,*)' psp6in : exit '
!write(6,*)' psp6in : indlmn(1:6,jj)'
!do jj=1,lmnmax
!write(6,*)indlmn(1:6,jj)
!end do
!ENDDEBUG

 deallocate(vpspll,rad,vloc,wfll)

end subroutine psp6in
!!***
