!{\src2tex{textfont=tt}}
!!****f* ABINIT/getgh1c
!!
!! NAME
!! getgh1c
!!
!! FUNCTION
!! Compute <G|H(1)|C> for input vector |C> expressed in reciprocal space
!! (H(1) is the 1st-order pertubed Hamiltonian)
!! Result is put in array gh1c.
!! <G|K(1)+Vnonlocal(1)|C> is also returned in gvnl1c.
!! if required, <G|S(1)|C> is returned in gs1c (S=overlap - PAW only)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  berryopt=option for Berry phase
!!  cplex=1 if vlocal1 is real, 2 if vlocal1 is complex
!!  cwave(2,npw*nspinor)=input wavefunction, in reciprocal space
!!  cwaveprj(ncprj,nspinor*usecprj)=<p_lmn|C> coefficients for wavefunction |C> (and 1st derivatives)
!!  dimekb=first dimension of ekb (see ekb_typ)
!!  dimffnlk=second dimension of ffnl (1+number of derivatives)
!!  dimffnl1=second dimension of ffnl1 and ffnlkq (1+number of derivatives)
!!  dkinpw(npw)=derivative of the (modified) kinetic energy for each plane wave at k (Hartree)
!!  ekb_typ(dimekb,1,nspinor**2)=
!!  ->Norm conserving : (Real) Kleinman-Bylander energies (hartree) for the displaced atom
!!          for number of basis functions (l,n) (lnmax)
!!          dimekb=lnmax
!!    ->PAW : (Real, symmetric) Frozen part of Dij coefficients
!!                               to connect projectors for the displaced atom
!!          for number of basis functions (l,m,n) (lmnmax)
!!          dimekb=lmnmax*(lmnmax+1)/2
!!          These are complex numbers in particular cases (spin-orbit)
!!               ekb_typ(:,:,1) contains Dij^up-up
!!               ekb_typ(:,:,2) contains Dij^dn-dn
!!               ekb_typ(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!               ekb_typ(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  ekb1_typ(dimekb,1,useekb1*nspinor**2)=1st derivative of ekb_typ for the current pertubation (see above)
!!  ffnlk(npw,dimffnlk,lmnmax,1)=nonloc form factors at k, for the displaced atom.
!!  ffnlkq(npw1,dimffnl1,lmnmax,1)=nonloc form fact at k+q for the displaced atom
!!  ffnl1(npw1,dimffnl1,lmnmax,ntypat)=nonloc form factors at k+q
!!  filstat=name of the status file
!!  gbound(2*mgfft+8,2)=G sphere boundary at k
!!  grad_berry(2,npw1*nspinor)= the gradient of the Berry phase term
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  idir=direction of the perturbation
!!  indlmn_typ(6,lmnmax,1)=indlmn info for the displaced atom
!!  ipert=type of the perturbation
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere at k.
!!  kg1_k(3,npw1)=coordinates of planewaves in basis sphere at k+q.
!!  kinpw1(npw1)=(modified) kinetic energy for each plane wave at k+q (Hartree)
!!  kpg_k(npw,nkpg)= (k+G) components at k (only if useylm=1)
!!  kpg1_k(npw1,nkpg1)= (k+G) components at k+q (only if useylm=1)
!!  kpt(3)=coordinates of k point.
!!  lmnmax= max number of (l,n) comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in unit cell.
!!  ncprj= dimension of cwaveprj array (1 for atom. displ. perturbation, else natom)
!!  nkpg,nkpg1=second dimensions of kpg_k and kpg1_k (0 if useylm=0)
!!  npw=number of planewaves in basis sphere at given k.
!!  npw1=number of planewaves in basis sphere at k+q
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell.
!!  n4,n5,n6 dimensions of vlocal and vlocal1
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  ph3d(2,npw1,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  pspso_typ(1)=spin-orbit info for the displaced atom
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H1|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S1|C> have to be computed in gs1 in addition to gh1
!!                       if -1, matrix elements <G|H1-lambda.S1|C> have to be computed in gh1 (gs1 not used)
!!  sij_typ(dimekb,1)=overlap matrix components (only if sij_opt/=0)
!!  tim_getgh1c=timing code of the calling subroutine (can be set to 0 if not attributed)
!!  usecprj=1 if cwaveprj coefficients (and 1st derivatives) are already in memory (PAW only)
!!  useekb1=1 if ekb derivatives (ekb1) exist
!!  vlocal1(cplex*n4,n5,n6)= 1st-order local pot in real space, on the augmented fft grid
!!
!! OUTPUT
!! gh1(2,npw*nspinor)=matrix elements <G|H(1)|C>.
!! gvnl1(2,npw1*nspinor)=matrix elements <G|Vnl(1)|C>
!! if (sij_opt=1)
!!  gs1(2,npw*nspinor)=matrix elements <G|S(1)|C> (S=overlap).
!!
!! SIDE EFFECTS
!! wfraug(2,n4,n5,n6)=is a dummy array
!!
!! PARENTS
!!      cgwf3
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,fourwf,leave_new,nonlop,status,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine getgh1c(berryopt,cplex,cwave,cwaveprj,dimekb,dimffnlk,dimffnl1,dkinpw,ekb_typ,ekb1_typ,ffnlk,ffnlkq,ffnl1,&
&                  filstat,gbound,gh1,grad_berry,gs1,gs_hamkq,gvnl1,idir,indlmn_typ,ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,&
&                  kpt,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,ncprj,nkpg,nkpg1,npw,npw1,nspinor,&
&                  ntypat,n4,n5,n6,paral_kgb,ph3d,prtvol,pspso_typ,sij_opt,sij_typ,tim_getgh1c,usecprj,useekb1,vlocal1,wfraug)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_13nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,cplex,dimekb,dimffnl1,dimffnlk,idir,ipert
 integer,intent(in) :: lmnmax,matblk,mgfft,mpsang,mpssoang,n4,n5,n6,natom,ncprj
 integer,intent(in) :: nkpg,nkpg1,npw,npw1,nspinor,ntypat,paral_kgb,prtvol
 integer,intent(in) :: sij_opt,tim_getgh1c,usecprj,useekb1
 real(dp),intent(in) :: lambda
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),indlmn_typ(6,lmnmax,1),kg1_k(3,npw1)
 integer,intent(in) :: kg_k(3,npw),pspso_typ(1)
 real(dp),intent(in) :: dkinpw(npw),ekb1_typ(dimekb,1,useekb1*nspinor**2)
 real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
 real(dp),intent(in) :: ffnl1(npw1,dimffnl1,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlk(npw,dimffnlk,lmnmax,1)
 real(dp),intent(in) :: ffnlkq(npw1,dimffnl1,lmnmax,1)
 real(dp),intent(in) :: grad_berry(2,npw1*nspinor*(berryopt/4)),kinpw1(npw1)
 real(dp),intent(in) :: kpg1_k(npw1,nkpg1),kpg_k(npw,nkpg),kpt(3)
 real(dp),intent(in) :: sij_typ(dimekb,gs_hamkq%usepaw)
 real(dp),intent(inout) :: cwave(2,npw*nspinor),ph3d(2,npw1,matblk)
 real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6),wfraug(2,n4,n5,n6)
 real(dp),intent(out) :: gh1(2,npw1*nspinor)
 real(dp),intent(out) :: gs1(2,npw1*nspinor*((sij_opt+1)/2))
 real(dp),intent(out) :: gvnl1(2,npw1*nspinor)
 type(cprj_type),intent(inout) :: cwaveprj(ncprj,nspinor*usecprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=20
 integer :: choice,cpopt,dimekb2_der,iexit,ipw,ipws,ispinor,istr,matblk_der,n1
 integer :: n2,n3,natom_der,nnlout,ntypat_der,paw_opt,shift1,shift2,shift3
 integer :: signs,tim_fourwf,tim_nonlop
 real(dp) :: arg,dum,weight
 character(len=500) :: message
!arrays
 integer :: atindx1_der(1),atindx_der(1),dimlmn_typ(1),nattyp_der(1)
 integer :: nloalg_der(5)
 real(dp) :: nonlop_dum(1,1),phkxredin(2,1),phkxredout(2,1),tsec(2),xred_der(3)
 real(dp),allocatable :: cwave_sp(:,:),enlout(:),gh1_sp(:,:),gvnl2(:,:)
 real(dp),allocatable :: ph1d_der(:,:),ph3din(:,:,:),ph3dout(:,:,:)
 type(cprj_type),allocatable :: cprj_tmp(:,:)

! *********************************************************************

!Keep track of total time spent in getgh1c
 call timab(206+tim_getgh1c,1,tsec)

 if(prtvol<0)then
  call status(0,filstat,iexit,level,'enter         ')
 end if

!Compatibility tests
 if(gs_hamkq%usepaw==0.and.usecprj==1)then
  write(message,'(4a)') ch10,&
&  ' getgh1c : BUG -',ch10,&
&  '   usecprj==1 not allowed for NC psps !'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 if(gs_hamkq%usepaw==1.and.ipert<=natom.and.useekb1==0)then
  write(message,'(4a)') ch10,&
&  ' getgh1c : BUG -',ch10,&
&  '   ekb derivatives must be allocated !'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!======================================================================
!== Apply the 1st-order local potential to the wavefunction
!======================================================================

!Phonon perturbation
!or Electric field perturbation
!or Strain perturbation
!-------------------------------------------
 if(ipert<=natom+4.and.ipert/=natom+1) then

  weight=one ; tim_fourwf=4
  call fourwf(cplex,vlocal1,cwave,gh1,wfraug,gbound,gs_hamkq%gbound,&
&  gs_hamkq%istwf_k,kg_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&  npw,npw1,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
  if(nspinor==2)then
   allocate(cwave_sp(2,npw),gh1_sp(2,npw1))
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cwave,cwave_sp,npw)
   do ipw=1,npw
    cwave_sp(1,ipw)=cwave(1,ipw+npw)
    cwave_sp(2,ipw)=cwave(2,ipw+npw)
   end do
!  $OMP END PARALLEL DO
   call fourwf(cplex,vlocal1,cwave_sp,gh1_sp,wfraug,gbound,gs_hamkq%gbound,&
&   gs_hamkq%istwf_k,kg_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&   npw,npw1,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(gh1,gh1_sp,npw1)
   do ipw=1,npw1
    gh1(1,ipw+npw1)=gh1_sp(1,ipw)
    gh1(2,ipw+npw1)=gh1_sp(2,ipw)
   end do
!  $OMP END PARALLEL DO
   deallocate(cwave_sp,gh1_sp)
  end if

! k-point perturbation
! -------------------------------------------
 else if(ipert==natom+1)then

! In the case of ddk operator, no local contribution
! (also because no self-consistency)
! $OMP PARALLEL DO PRIVATE(ipw) &
! $OMP&SHARED(gh1,npw1,nspinor)
  do ipw=1,npw1*nspinor
   gh1(:,ipw)=zero
  end do
! $OMP END PARALLEL DO
! Magnetic field perturbation (currently mimics electric field for testing)
! -------------------------------------------
 else if(ipert == natom+5) then

  weight=one ; tim_fourwf=4
  call fourwf(cplex,vlocal1,cwave,gh1,wfraug,gbound,gs_hamkq%gbound,&
&  gs_hamkq%istwf_k,kg_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&  npw,npw1,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
  if(nspinor==2)then
   allocate(cwave_sp(2,npw),gh1_sp(2,npw1))
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cwave,cwave_sp,npw)
   do ipw=1,npw
    cwave_sp(1,ipw)=cwave(1,ipw+npw)
    cwave_sp(2,ipw)=cwave(2,ipw+npw)
   end do
!  $OMP END PARALLEL DO
   call fourwf(cplex,vlocal1,cwave_sp,gh1_sp,wfraug,gbound,gs_hamkq%gbound,&
&   gs_hamkq%istwf_k,kg_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&   npw,npw1,n4,n5,n6,2,paral_kgb,tim_fourwf,weight)
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(gh1,gh1_sp,npw1)
   do ipw=1,npw1
    gh1(1,ipw+npw1)=gh1_sp(1,ipw)
    gh1(2,ipw+npw1)=gh1_sp(2,ipw)
   end do
!  $OMP END PARALLEL DO
   deallocate(cwave_sp,gh1_sp)
  end if

 end if

!======================================================================
!== Apply the 1st-order non-local potential to the wavefunction
!======================================================================

!Phonon perturbation
!-------------------------------------------
 if(ipert<=natom) then

  signs=2 ; nnlout=3 ; natom_der=1 ; nattyp_der(1)=1 ; ntypat_der=1
  dimekb2_der=1 ; allocate(enlout(nnlout))
  matblk_der=1 ; tim_nonlop=7
  xred_der(:)=gs_hamkq%xred(:,ipert)
  atindx_der(1)=1 ; atindx1_der(1)=1
  n1=gs_hamkq%ngfft(1) ; n2=gs_hamkq%ngfft(2) ; n3=gs_hamkq%ngfft(3)

! Store at the right place the 1d phases
  allocate(ph1d_der(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
  shift1=(gs_hamkq%atindx(ipert)-1)*(2*n1+1)
  ph1d_der(:,1:2*n1+1)=gs_hamkq%ph1d(:,1+shift1:2*n1+1+shift1)
  shift2=(gs_hamkq%atindx(ipert)-1)*(2*n2+1)+natom*(2*n1+1)
  ph1d_der(:,1+2*n1+1:2*n2+1+2*n1+1)=gs_hamkq%ph1d(:,1+shift2:2*n2+1+shift2)
  shift3=(gs_hamkq%atindx(ipert)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
  ph1d_der(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=&
&  gs_hamkq%ph1d(:,1+shift3:2*n3+1+shift3)

! Will compute the 3D phase factors inside nonlop
  allocate(ph3din(2,npw,1),ph3dout(2,npw1,1))
  nloalg_der(:)=gs_hamkq%nloalg(:)
  nloalg_der(1)=-abs(gs_hamkq%nloalg(1))
  nloalg_der(4)=1

! Compute here phkxred for kpt and kpq
  arg=two_pi*(kpt(1)*gs_hamkq%xred(1,ipert) &
&  +kpt(2)*gs_hamkq%xred(2,ipert) &
&  +kpt(3)*gs_hamkq%xred(3,ipert))
  phkxredin(1,1)=cos(arg)  ; phkxredin(2,1)=sin(arg)
  arg=two_pi*(gs_hamkq%kpoint(1)*gs_hamkq%xred(1,ipert) &
&  +gs_hamkq%kpoint(2)*gs_hamkq%xred(2,ipert) &
&  +gs_hamkq%kpoint(3)*gs_hamkq%xred(3,ipert))
  phkxredout(1,1)=cos(arg) ; phkxredout(2,1)=sin(arg)

! Application of 1st-order nl operator

! PAW:
  if (gs_hamkq%usepaw==1) then
   allocate(gvnl2(2,npw1*nspinor))
   if (usecprj==1) then
!   A- The <p_i|c> (cprj) factors (and derivatives) are already in memory
!   A1- Compute derivatives due to projectors |p_i>
    cpopt=4 ; choice=2
    paw_opt=1;if (sij_opt/=0) paw_opt=sij_opt+3
    call nonlop(atindx1_der,choice,cpopt,cwaveprj,dimekb,dimekb2_der,dimffnlk,dimffnl1,&
&    ekb_typ,enlout,ffnlk,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&    indlmn_typ,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&    lambda,lmnmax,matblk_der,mgfft,mpi_enreg,mpsang,mpssoang,natom_der,nattyp_der,&
&    gs_hamkq%ngfft,nkpg,nkpg1,nloalg_der,nnlout,npw,npw1,nspinor,ntypat_der,&
&    0,paw_opt,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,pspso_typ,signs,&
&    sij_typ,gs1,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1)
!   A2- Compute derivatives due to Ekb (D_ij)
    cpopt=2 ; choice=1 ; paw_opt=1
    call nonlop(atindx1_der,choice,cpopt,cwaveprj,dimekb,dimekb2_der,dimffnlk,dimffnl1,&
&    ekb1_typ,enlout,ffnlk,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&    indlmn_typ,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&    lambda,lmnmax,matblk_der,mgfft,mpi_enreg,mpsang,mpssoang,natom_der,nattyp_der,&
&    gs_hamkq%ngfft,nkpg,nkpg1,nloalg_der,nnlout,npw,npw1,nspinor,ntypat_der,&
&    0,paw_opt,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,pspso_typ,signs,&
&    nonlop_dum,nonlop_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl2)
   else
!   B- The <p_i|c> (cprj) factors (and derivatives) have to be computed
    allocate(cprj_tmp(1,nspinor))
    dimlmn_typ(1)=count(indlmn_typ(3,:,1)>0)
    call cprj_alloc(cprj_tmp,0,dimlmn_typ)
!   B1- Compute derivatives due to projectors |p_i>
    cpopt=0 ; choice=2
    paw_opt=1;if (sij_opt/=0) paw_opt=sij_opt+3
    call nonlop(atindx1_der,choice,cpopt,cprj_tmp,dimekb,dimekb2_der,dimffnlk,dimffnl1,&
&    ekb_typ,enlout,ffnlk,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&    indlmn_typ,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&    lambda,lmnmax,matblk_der,mgfft,mpi_enreg,mpsang,mpssoang,natom_der,nattyp_der,&
&    gs_hamkq%ngfft,nkpg,nkpg1,nloalg_der,nnlout,npw,npw1,nspinor,ntypat_der,&
&    0,paw_opt,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,pspso_typ,signs,&
&    sij_typ,gs1,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1)
!   B2- Compute derivatives due to Ekb (D_ij)
    cpopt=2 ; choice=1 ; paw_opt=1
    call nonlop(atindx1_der,choice,cpopt,cprj_tmp,dimekb,dimekb2_der,dimffnlk,dimffnl1,&
&    ekb1_typ,enlout,ffnlk,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&    indlmn_typ,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&    lambda,lmnmax,matblk_der,mgfft,mpi_enreg,mpsang,mpssoang,natom_der,nattyp_der,&
&    gs_hamkq%ngfft,nkpg,nkpg1,nloalg_der,nnlout,npw,npw1,nspinor,ntypat_der,&
&    0,paw_opt,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,pspso_typ,signs,&
&    nonlop_dum,nonlop_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl2)
    call cprj_free(cprj_tmp)
    deallocate(cprj_tmp)
   end if
!  Sum up both contributions to Vnl(1):
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(gvnl1,gvnl2,npw1)
   do ipw=1,npw1
    gvnl1(1,ipw)=gvnl1(1,ipw)+gvnl2(1,ipw)
    gvnl1(2,ipw)=gvnl1(2,ipw)+gvnl2(2,ipw)
   end do
!  $OMP END PARALLEL DO

!  Norm-conserving psps:
  else
!  Compute only derivatives due to projectors |p_i>
   cpopt=-1 ; choice=2 ; paw_opt=0
   call nonlop(atindx1_der,choice,cpopt,cwaveprj,dimekb,dimekb2_der,dimffnlk,dimffnl1,&
&   ekb_typ,enlout,ffnlk,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,idir,&
&   indlmn_typ,gs_hamkq%istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,gs_hamkq%kpoint,&
&   lambda,lmnmax,matblk_der,mgfft,mpi_enreg,mpsang,mpssoang,natom_der,nattyp_der,&
&   gs_hamkq%ngfft,nkpg,nkpg1,nloalg_der,nnlout,npw,npw1,nspinor,ntypat_der,&
&   0,paw_opt,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,pspso_typ,signs,&
&   nonlop_dum,nonlop_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1)
  end if

  deallocate(enlout,ph1d_der,ph3din,ph3dout)

! k-point perturbation
! -------------------------------------------
 else if(ipert==natom+1)then

! Remember, q=0, so can take all RF data...
  choice=5 ; nnlout=1 ; signs=2 ; tim_nonlop=8
  cpopt=-1+4*usecprj ; paw_opt=0 ; allocate(enlout(nnlout))
  call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj,gs_hamkq%dimekb1,gs_hamkq%dimekb2,&
&  dimffnl1,dimffnl1,gs_hamkq%ekb,enlout,ffnl1,ffnl1,gs_hamkq%gmet,&
&  gs_hamkq%gprimd,idir,gs_hamkq%indlmn,gs_hamkq%istwf_k,kg1_k,kg1_k,&
&  kpg1_k,kpg1_k,gs_hamkq%kpoint,gs_hamkq%kpoint,lambda,lmnmax,matblk,&
&  mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&  nkpg1,nkpg1,gs_hamkq%nloalg,nnlout,npw1,npw1,nspinor,ntypat,paw_opt,&
&  0,gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,gs_hamkq%pspso,&
&  signs,nonlop_dum,nonlop_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1)
  deallocate(enlout)

! Electric field perturbation without Berry phase
! -------------------------------------------
 else if(ipert==natom+2.and.berryopt/=4)then
! gvnl1 was already initialized in the calling routine, by reading a ddk file

! Magnetic field perturbation (mimics electric field for testing)
! -------------------------------------------
 else if(ipert==natom+5) then
! gvnl1 was already initialized in the calling routine, by reading a ddk file

! Electric field perturbation with Berry phase
! -------------------------------------------
 else if(ipert==natom+2.and.berryopt==4)then

  do ipw=1,npw1*nspinor
   gvnl1(1,ipw)=-grad_berry(2,ipw)
   gvnl1(2,ipw)= grad_berry(1,ipw)
  end do

! Strain perturbation
! -------------------------------------------
 else if(ipert==natom+3 .or. ipert==natom+4)then

  if(ipert==natom+3) then
   istr=idir
  else
   istr = idir+3
  end if

! Remember, q=0, so can take all RF data)
! Copied from d/dk above; changes may be needed; tim_nonlop may need changing
  signs=2 ; choice=3 ; nnlout=6 ; tim_nonlop=8
  cpopt=-1+4*usecprj ; paw_opt=0 ; allocate(enlout(nnlout))
  call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj,gs_hamkq%dimekb1,gs_hamkq%dimekb2,&
&  dimffnl1,dimffnl1,gs_hamkq%ekb,enlout,ffnl1,ffnl1,gs_hamkq%gmet,&
&  gs_hamkq%gprimd,istr,gs_hamkq%indlmn,gs_hamkq%istwf_k,kg1_k,kg1_k,&
&  kpg1_k,kpg1_k,gs_hamkq%kpoint,gs_hamkq%kpoint,lambda,lmnmax,matblk,&
&  mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&  nkpg1,nkpg1,gs_hamkq%nloalg,nnlout,npw1,npw1,nspinor,ntypat,paw_opt,&
&  0,gs_hamkq%phkxred,gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,gs_hamkq%pspso,&
&  signs,nonlop_dum,nonlop_dum,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave,gvnl1)
  deallocate(enlout)

 end if

!======================================================================
!== Apply the 1st-order kinetic operator to the wavefunction
!== (add it to nl contribution)
!======================================================================

!Phonon perturbation
!or Electric field perturbation
!-------------------------------------------
!No kinetic contribution

!k-point perturbation
!or Strain perturbation
!-------------------------------------------
 if(ipert==natom+1 .or. &
& ipert==natom+3 .or. ipert==natom+4) then

! Remember that npw=npw1 for ddk perturbation
  do ispinor=1,nspinor
!  $OMP PARALLEL DO PRIVATE(ipw,ipws) &
!  $OMP&SHARED(cwave,ispinor,gvnl1,dkinpw,kinpw1,npw,nspinor)
   do ipw=1,npw
    ipws=ipw+npw*(ispinor-1)
    if(kinpw1(ipw)<huge(zero)*1.d-11)then
     gvnl1(1,ipws)=gvnl1(1,ipws)+dkinpw(ipw)*cwave(1,ipws)
     gvnl1(2,ipws)=gvnl1(2,ipws)+dkinpw(ipw)*cwave(2,ipws)
    else
     gvnl1(1,ipws)=zero
     gvnl1(2,ipws)=zero
    end if
   end do
!  $OMP END PARALLEL DO
  end do

 end if

!======================================================================
!== Sum contributions to get the application of H(1) to the wf
!======================================================================
!Also filter the wavefunctions for large modified kinetic energy

 do ispinor=1,nspinor
  ipws=(ispinor-1)*npw1
! $OMP PARALLEL DO PRIVATE(ipw) &
! $OMP&SHARED(gh1,gvnl1,kinpw1,ipws,npw1)
  do ipw=1+ipws,npw1+ipws
   if(kinpw1(ipw-ipws)<huge(zero)*1.d-11)then
    gh1(1,ipw)=gh1(1,ipw)+gvnl1(1,ipw)
    gh1(2,ipw)=gh1(2,ipw)+gvnl1(2,ipw)
   else
    gh1(1,ipw)=zero
    gh1(2,ipw)=zero
   end if
  end do
! $OMP END PARALLEL DO
 end do

 if(prtvol<0)then
  call status(0,filstat,iexit,level,'exit          ')
 end if

 call timab(206+tim_getgh1c,1,tsec)

end subroutine getgh1c
!!***
