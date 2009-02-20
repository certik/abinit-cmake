!{\src2tex{textfont=tt}}
!!****f* ABINIT/cgwf3
!! NAME
!! cgwf3
!!
!!
!! FUNCTION
!! Update one single wavefunction (cwavef), non self-consistently.
!! Uses a conjugate-gradient algorithm.
!! Try to keep close to the formulas in PRB55, 10337 (1997), for the
!! non-self-consistent case, except that we are computing here
!! the second-derivative of the total energy, and not E(2). There
!! is a factor of 2 between the two quantities ...
!! The wavefunction that is generated is always orthogonal to cgq .
!! It is orthogonal to the active Hilbert space, and will be complemented
!! by contributions from the active space in the calling routine, if needed.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG,DRH,XW,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  band=which particular band we are converging.
!!  berryopt=option for Berry phase
!!  cgq(2,mcgq)=wavefunction coefficients for ALL bands.
!!  cplex=1 if vlocal1 is real, 2 if vlocal1 is complex
!!  cwave0(2,npw*nspinor)=GS wavefunction at k, in reciprocal space
!!  cwaveprj0(ncprj,nspinor*usecprj)= wave functions at k projected with nl projectors
!!  dimekb=first dimension of ekb (see ekb_typ)
!!  dimffnlk=second dimension of ffnl (1+number of derivatives)
!!  dimffnl1=second dimension of ffnl1 and ffnlkq (1+number of derivatives)
!!  dkinpw(npw)=derivative of the (modified) kinetic energy for each plane wave at k (Hartree)
!!  eig0nk=0-order eigenvalue for the present wavefunction at k
!!  eig0_kq(nband)=GS eigenvalues at k+Q (hartree)
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
!!  grad_berry(2,mpw1,dtefield%nband_occ) = the gradient of the Berry phase term
!!  gbound(2*mgfft+8,2)=G sphere boundary at k
!!  gscq(2,mgscq)=<g|S|Cnk+q> coefficients for ALL bands (PAW)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  icgq=shift to be applied on the location of data in the array cgq
!!  igscq=shift to be applied on the location of data in the array gscq
!!  idir=direction of the perturbation
!!  indlmn_typ(6,lmnmax,1)=indlmn info for the displaced atom
!!  ipert=type of the perturbation
!!  kg_k(3,npw)=coordinates of planewaves in basis sphere at k.
!!  kg1_k(3,npw1)=coordinates of planewaves in basis sphere at k+q.
!!  kinpw1(npw1)=(modified) kinetic energy for each plane wave at k+q (Hartree)
!!  kpg_k(npw,nkpg)= (k+G) components at k (only if useylm=1)
!!  kpg1_k(npw1,nkpg1)= (k+G) components at k+q (only if useylm=1)
!!  kpt(3)=coordinates of k point.
!!  lmnmax= max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mcgq=second dimension of the cgq array
!!  mgfft=maximum size of 1D FFTs
!!  mgscq=second dimension of gscq
!!  mpi_enreg=informations about MPI parallelization
!!  mpw1=maximum number of planewave for first-order wavefunctions
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband=number of bands.
!!  nbdbuf=number of buffer bands for the minimisation
!!  ncprj= dimension of cwaveprj0 array (1 for atom. displ. perturbation, else natom)
!!  nkpg,nkpg1=second dimensions of kpg_k and kpg1_k (0 if useylm=0)
!!  nline=number of line minimizations per band.
!!  npw=number of planewaves in basis sphere at given k.
!!  npw1=number of planewaves in basis sphere at k+q
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell.
!!  n4,n5,n6 dimensions of vlocal and vlocal1
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  ph3d(2,npw1,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  pspso_typ(1)=spin-orbit info for the displaced atom
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  quit= if 1, proceeds to smooth ending of the job.
!!  sciss=scissor shift (Ha)
!!  sij_typ(dimekb,gs_hamkq%usepaw)=overlap matrix components (PAW)
!!  tolwfr=tolerance on largest wf residual
!!  usecprj= 1 if cwaveprj0 array is stored in memory
!!  useekb1=1 if ekb derivatives (ekb1) exist
!!  usedcwavef=flag controlling the use of dcwavef array (PAW only):
!!             0: not used (not allocated)
!!             1: used as input
!!             2: used as output
!!  vlocal(n4,n5,n6)= GS local pot in real space, on the augmented fft grid
!!  vlocal1(cplex*n4,n5,n6)= RF local pot in real space, on the augmented fft grid
!!  wfoptalg=govern the choice of algorithm for wf optimisation (0 or 10, at present)
!!
!! OUTPUT
!!  eig1_k(2*nband**2)=matrix of first-order eigenvalues (hartree)
!!                     eig1(:,ii,jj)=<C0 ii|H1|C0 jj> for norm-conserving psps
!!                     eig1(:,ii,jj)=<C0 ii|H1-(eig0_k+eig0_k+q)/2.S(1)|C0 jj> for PAW
!!  ghc(2,npw1*nspinor)=<G|H0|C1 band,k>,
!!  gvnlc(2,npw1*nspinor)=<G|Vnl|C1 band,k>
!!  gvnl1(2,npw1*nspinor)=<G|Vnl1|C0 band,k>
!!  resid=wf residual for current band
!!  gh1_n= <G|H1|C0 band,k> (NCPP) or <G|H1-eig0k.S1|C0 band,k>
!!
!! SIDE EFFECTS
!!  Input/Output:
!!  cwavef(2,npw1*nspinor)=first-order  wavefunction at k,q, in reciprocal space (updated)
!!  wfraug(2,n4,n5,n6)=is a dummy array
!!  === if gs_hamkq%usepaw==1 ===
!!  cwaveprj(ncprj,nspinor)= wave functions at k projected with nl projectors
!!  === if usedcwavef>0 ===
!!  dcwavef(2,npw1*nspinor)=change of wavefunction due to change of overlap:
!!         dcwavef is delta_Psi(1)=Sum_{j}[<C0_k+q_j|S(1)|C0_k_i>.|C0_k+q_j>]
!!         see PRB 78, 035105 (2008), Eq. (42)
!!         input  if usedcwavef=1, output if usedcwavef=1
!!
!! PARENTS
!!      vtowfk3
!!
!! CHILDREN
!!      dotprod_g,getghc,getgh1c,leave_new,precon,projbd,sqnorm_g,status
!!      timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cgwf3(band,berryopt,cgq,cplex,cwavef,cwave0,cwaveprj,cwaveprj0,dcwavef,gh1_n,dimekb,&
& dimffnlk,dimffnl1,dkinpw,eig0nk,eig0_kq,eig1_k,&
& ekb_typ,ekb1_typ,ffnlk,ffnlkq,ffnl1,filstat,gbound,ghc,grad_berry,gscq,&
& gs_hamkq,gvnlc,gvnl1,icgq,idir,indlmn_typ,&
& ipert,igscq,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,&
& kpt,lmnmax,matblk,mcgq,mgfft,mgscq,mpi_enreg,mpsang,mpssoang,mpw1,natom,nband,nbdbuf,&
& ncprj,nkpg,nkpg1,nline,npw,npw1,nspinor,ntypat,n4,n5,n6,ortalg,paral_kgb,ph3d,prtvol,pspso_typ,&
& qphon,quit,resid,sciss,sij_typ,tolwfr,usecprj,usedcwavef,useekb1,vlocal,vlocal1,wfoptalg,wfraug)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12spacepar
 use interfaces_14wfs
 use interfaces_16response, except_this_one => cgwf3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: band,berryopt,cplex,dimekb,dimffnl1,dimffnlk,icgq,idir
 integer,intent(in) :: igscq,ipert,lmnmax,matblk,mcgq,mgfft,mgscq,mpsang
 integer,intent(in) :: mpssoang,mpw1,n4,n5,n6,natom,nband,nbdbuf,ncprj,nkpg
 integer,intent(in) :: nkpg1,nline,npw,npw1,nspinor,ntypat,ortalg,paral_kgb
 integer,intent(in) :: prtvol,quit,usecprj,usedcwavef,useekb1,wfoptalg
 real(dp),intent(in) :: eig0nk,sciss,tolwfr
 real(dp),intent(out) :: resid
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),indlmn_typ(6,lmnmax,1),kg1_k(3,npw1)
 integer,intent(in) :: kg_k(3,npw),pspso_typ(1)
 real(dp),intent(in) :: cgq(2,mcgq),dkinpw(npw),eig0_kq(nband)
 real(dp),intent(in) :: ekb1_typ(dimekb,1,useekb1*nspinor**2)
 real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
 real(dp),intent(in) :: ffnl1(npw1,dimffnl1,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlk(npw,dimffnlk,lmnmax,1)
 real(dp),intent(in) :: ffnlkq(npw1,dimffnl1,lmnmax,1)
 real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband),gscq(2,mgscq)
 real(dp),intent(in) :: kinpw1(npw1),kpg1_k(npw1,nkpg1),kpg_k(npw,nkpg),kpt(3)
 real(dp),intent(in) :: qphon(3),sij_typ(dimekb,gs_hamkq%usepaw)
 real(dp),intent(inout) :: cwave0(2,npw*nspinor),cwavef(2,npw1*nspinor)
 real(dp),intent(inout) :: dcwavef(2,npw1*nspinor*((usedcwavef+1)/2))
 real(dp),intent(inout) :: ph3d(2,npw1,matblk),vlocal(n4,n5,n6)
 real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6),wfraug(2,n4,n5,n6)
 real(dp),intent(out) :: eig1_k(2*nband**2),gh1_n(2,npw1*nspinor)
 real(dp),intent(out) :: ghc(2,npw1*nspinor),gvnl1(2,npw1*nspinor)
 real(dp),intent(out) :: gvnlc(2,npw1*nspinor)
 type(cprj_type),intent(inout) :: cwaveprj(ncprj,nspinor)
 type(cprj_type),intent(inout) :: cwaveprj0(ncprj,nspinor*usecprj)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=19,tim_getgh1c=1,tim_getghc=2
 integer,save :: nskip=0
 integer :: i1,i2,i3,iband,iexit,igs,ii,iline,indx_cgq,iprint,ipw,ipw1,ipws
 integer :: ispinor,istr,istwf_k,jband,nvloc,sij_opt,tim_projbd,useoverlap
 real(dp) :: cgwftol,chc,costh,d2edt2,d2te,d2teold,dedt,deltae,deold,dotgg
 real(dp) :: dotgp,doti,dotr,dum,eshift,eshiftkq,gamma,ghc2,optekin,prod1,prod2
 real(dp) :: root,sinth,swap,tan2th,theta,u1h0me0u1,use_vnl,weight,wfim,wfre
 real(dp) :: xnorm
 logical :: gen_eigenpb
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: conjgr(:,:),cwaveq(:,:),direc(:,:),gberry(:,:)
 real(dp),allocatable :: gh1(:,:),gh_direc(:,:),gresid(:,:),gs1(:,:)
 real(dp),allocatable :: gsc_dummy(:,:),gvnl_direc(:,:),pcon(:),scprod(:,:)

! *********************************************************************

!DEBUG
!write(6,*)' cgwf3 : enter'
!ENDDEBUG

!======================================================================
!========= LOCAL VARIABLES DEFINITIONS AND ALLOCATIONS ================
!======================================================================

 call timab(122,1,tsec)
 if(prtvol<0) call status(0,filstat,iexit,level,'enter         ')
 iprint=0;if(prtvol==-level)iprint=1

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' cgwf3 : enter '
  call wrtout(06,message,'PERS')
 end if

!Tell us what is going on:
 if(prtvol>=10)then
  write(message, '(a,i6,2x,a,i3,a)' ) &
&  ' --- cgwf3 is called for band',band,'for',nline,' lines'
  call  wrtout(06,message,'PERS')
 end if

!if PAW, one has to solve a generalized eigenproblem
 gen_eigenpb=(gs_hamkq%usepaw==1)
 useoverlap=0;if (gen_eigenpb) useoverlap=1

!if PAW, no need to compute Vnl contributions
 use_vnl=0;if (gs_hamkq%usepaw==0) use_vnl=1
 if (gen_eigenpb.and.(use_vnl==1)) stop "Error in cgwf3: contact Abinit group"

!Use scissor shift on 0-order eigenvalue
 eshift=eig0nk-sciss

!Additional initializations
 nvloc=1 !For the time being, non-collinear potential is not allowed in RF
 istwf_k=gs_hamkq%istwf_k
 optekin=0;if (wfoptalg>=10) optekin=1
 tim_projbd=2

!Memory allocations
 allocate(gh1(2,npw1*nspinor),pcon(npw1),scprod(2,nband))
 if (berryopt==4) then
  allocate(gberry(2,npw1*nspinor))
  gberry(:,1:npw1*nspinor)=grad_berry(:,1:npw1*nspinor,band)
 end if

!======================================================================
!========== INITIALISATION OF MINIMIZATION ITERATIONS =================
!======================================================================

!Compute H(1) applied to GS wavefunction Psi(0)
 if (gen_eigenpb) then
  sij_opt=1;allocate(gs1(2,npw1*nspinor))
 else
  sij_opt=0
 end if
 call getgh1c(berryopt,cplex,cwave0,cwaveprj0,dimekb,dimffnlk,dimffnl1,dkinpw,ekb_typ,ekb1_typ,ffnlk,ffnlkq,ffnl1,&
& filstat,gbound,gh1,gberry,gs1,gs_hamkq,gvnl1,idir,indlmn_typ,ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,&
& kpt,eshift,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,ncprj,nkpg,nkpg1,npw,npw1,nspinor,&
& ntypat,n4,n5,n6,paral_kgb,ph3d,prtvol,pspso_typ,sij_opt,sij_typ,tim_getgh1c,usecprj,useekb1,vlocal1,wfraug)

 if (gen_eigenpb) then
! $OMP PARALLEL DO PRIVATE(ipw) &
! $OMP&SHARED(eshift,gh1,gs1,npw1,nspinor)
  do ipw=1,npw1*nspinor
   gh1(1,ipw)=gh1(1,ipw)-eshift*gs1(1,ipw)
   gh1(2,ipw)=gh1(2,ipw)-eshift*gs1(2,ipw)
  end do

! If generalized eigenPb and dcwavef requested, compute it:
! dcwavef is delta_Psi(1)=Sum_{j}[<C0_k+q_j|S(1)|C0_k_i>.|C0_k+q_j>]
! see PRB 78, 035105 (2008), Eq. (42)
  if (usedcwavef==2) then
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(dcwavef,gs1,npw1,nspinor)
   do ipw=1,npw1*nspinor
    dcwavef(1,ipw)=gs1(1,ipw)
    dcwavef(2,ipw)=gs1(2,ipw)
   end do
!  Note the subtlety: projb is called with useoverlap=0 and gs1
!  in order to get Sum[<cgq|s1|c0>|cgq>]=Sum[<cgq|gs1>|cgq>]
!  $OMP END PARALLEL DO
   call projbd(cgq,dcwavef,-1,icgq,0,istwf_k,mcgq,mpi_enreg,0,nband,npw1,nspinor,&
&   ortalg,iprint,gsc_dummy,scprod,tim_projbd,0)
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(dcwavef,gs1,npw1,nspinor)
   do ipw=1,npw1*nspinor
    dcwavef(1,ipw)=gs1(1,ipw)-dcwavef(1,ipw)
    dcwavef(2,ipw)=gs1(2,ipw)-dcwavef(2,ipw)
   end do
!  $OMP END PARALLEL DO
  end if
 end if

!$OMP PARALLEL DO PRIVATE(ipw) &
!$OMP&SHARED(gh1,gh1_n,npw1,nspinor)
 do ipw=1,npw1*nspinor
  gh1_n(1,ipw)=gh1(1,ipw)
  gh1_n(2,ipw)=gh1(2,ipw)
 end do
!$OMP END PARALLEL DO

!Projecting out all bands (this could be avoided)
!Note the subtleties:
!-For the generalized eigenPb, S|cgq> is used in place of |cgq>,
!in order to apply P_c+ projector (see PRB 73, 235101 (2006), Eq. (71), (72)
!-For the simple eigenPb, gscq is used as dummy argument
 if(gen_eigenpb)then
  call projbd(gscq,gh1,-1,igscq,icgq,istwf_k,mgscq,mpi_enreg,mcgq,nband,npw1,nspinor,&
&  ortalg,iprint,cgq,scprod,tim_projbd,useoverlap)
 else
  call projbd(cgq,gh1,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&  ortalg,iprint,gscq,scprod,tim_projbd,useoverlap)
 end if

!The array eig1_k contains:
!<u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)>                           (NC psps)
!or <u_(band,k+q)^(0)|H_(k+q,k)^(1)-(eig0_k+eig0_k+q)/2.S^(1)|u_(band,k)^(0)> (PAW)
 jband=(band-1)*2*nband
 if (gen_eigenpb) then
  indx_cgq=icgq
  do iband=1,nband
   eshiftkq=half*(eig0_kq(iband)-eig0nk)
   call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1*nspinor,1,cgq(:,indx_cgq+1:indx_cgq+npw*nspinor),gs1)
   eig1_k(2*iband-1+jband)=scprod(1,iband)-eshiftkq*dotr
   eig1_k(2*iband  +jband)=scprod(2,iband)-eshiftkq*doti
   indx_cgq=indx_cgq+npw*nspinor
  end do
 else
  do iband=1,nband
   eig1_k(2*iband-1+jband)=scprod(1,iband)
   eig1_k(2*iband  +jband)=scprod(2,iband)
  end do
 end if

!No more need of gs1
 if (gen_eigenpb) deallocate(gs1)

!Filter the wavefunctions for large modified kinetic energy (see routine mkkin.f)
 do ispinor=1,nspinor
  ipws=(ispinor-1)*npw1
! $OMP PARALLEL DO PRIVATE(ipw) &
! $OMP&SHARED(cwavef,kinpw1,ipws,npw1)
  do ipw=1+ipws,npw1+ipws
   if(kinpw1(ipw-ipws)>huge(zero)*1.d-11)then
    cwavef(1,ipw)=zero
    cwavef(2,ipw)=zero
   end if
  end do
! $OMP END PARALLEL DO
 end do

!Apply the orthogonality condition: <C1 k,q|C0 k+q>=0
!Project out all bands from cwavef, i.e. apply P_c projector on cwavef
!(this is needed when there are some partially or unoccupied states)
 call projbd(cgq,cwavef,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
& ortalg,iprint,gscq,scprod,tim_projbd,useoverlap)

!If generalized eigenPb, the orthogonality condition is
!<C1 k,q|S0|C0 k+q>+1/2<C0 k|S1|C0 k+q>=0
 if (gen_eigenpb) then
! $OMP PARALLEL DO PRIVATE(ipw) &
! $OMP&SHARED(cwavef,dcwavef,npw1,nspinor)
  do ipw=1,npw1*nspinor
   cwavef(1,ipw)=cwavef(1,ipw)-half*dcwavef(1,ipw)
   cwavef(2,ipw)=cwavef(2,ipw)-half*dcwavef(2,ipw)
  end do
! $OMP END PARALLEL DO
 end if

!Treat the case of buffer bands
 if(band>max(1,nband-nbdbuf))then
! $OMP PARALLEL DO PRIVATE(ipw) &
! $OMP&SHARED(cwavef,ghc,npw1,nspinor)
  do ipw=1,npw1*nspinor
   cwavef(1,ipw)=zero
   cwavef(2,ipw)=zero
   ghc(1,ipw)  =zero
   ghc(2,ipw)  =zero
  end do
! $OMP END PARALLEL DO
  if (use_vnl==1) then
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(gvnlc,npw1,nspinor)
   do ipw=1,npw1*nspinor
    gvnlc(1,ipw)=zero
    gvnlc(2,ipw)=zero
   end do
!  $OMP END PARALLEL DO
  end if
! A small negative residual will be associated with these
  resid=-0.1_dp
! Number of one-way 3D ffts skipped
  nskip=nskip+nline

 else
! If not a buffer band, perform the optimisation

  allocate(conjgr(2,npw1*nspinor))
  allocate(direc(2,npw1*nspinor))
  allocate(gresid(2,npw1*nspinor))
  allocate(cwaveq(2,npw1*nspinor))

  cwaveq(:,:)=cgq(:,1+npw1*nspinor*(band-1)+icgq:npw1*nspinor*band+icgq)
  dotgp=one

! Here apply H(0) at k+q to input orthogonalized 1st-order wfs
  if (prtvol<0) call status(0,filstat,iexit,level,'call getghc(1)')
  sij_opt=0;if (gen_eigenpb) sij_opt=-1
  call getghc(cwavef,dimffnl1,ffnl1,filstat,&
&  ghc,gsc_dummy,gs_hamkq,gvnlc,kg1_k,kinpw1,eshift,lmnmax,matblk,&
&  mgfft,mpi_enreg,mpsang,mpssoang,natom,1,npw1,nspinor,ntypat,nvloc,n4,n5,n6,&
&  paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
! ghc also includes the eigenvalue shift
! (already included if generalized eigenPb - see sijopt=-1 above)
  if (.not.gen_eigenpb) then
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cwavef,eshift,ghc,npw1,nspinor)
   do ipw=1,npw1*nspinor
    ghc(1,ipw)=ghc(1,ipw)-eshift*cwavef(1,ipw)
    ghc(2,ipw)=ghc(2,ipw)-eshift*cwavef(2,ipw)
   end do
!  $OMP END PARALLEL DO
  end if
  if (prtvol<0) call status(0,filstat,iexit,level,'after getghc(1')

! Initialize resid, in case of nline==0
  resid=zero

! ======================================================================
! ====== BEGIN LOOP FOR A GIVEN BAND: MINIMIZATION ITERATIONS ==========
! ======================================================================

  do iline=1,nline
   if (prtvol<0) call status(iline,filstat,iexit,level,'loop iline    ')

!  ======================================================================
!  ================= COMPUTE THE RESIDUAL ===============================
!  ======================================================================
!  Note that gresid (=steepest-descent vector, Eq.(26) of PRB 55, 10337 (1996))
!  is precomputed to garantee cancellation of errors
!  and allow residuals to reach values as small as 1.0d-24 or better.
   if (berryopt == 4) then
    if (ipert==natom+2) then
     gvnl1=zero
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(ghc,gh1,npw1,nspinor,gresid,gberry,band)
     do ipw=1,npw1*nspinor
      gresid(1,ipw)=-ghc(1,ipw)-gh1(1,ipw)
      gresid(2,ipw)=-ghc(2,ipw)-gh1(2,ipw)
     end do
!    $OMP END PARALLEL DO
    else
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(ghc,gh1,npw1,nspinor,gresid,gberry,band)
     do ipw=1,npw1*nspinor
      gresid(1,ipw)=-ghc(1,ipw)-gh1(1,ipw)+gberry(2,ipw)
      gresid(2,ipw)=-ghc(2,ipw)-gh1(2,ipw)-gberry(1,ipw)
     end do
!    $OMP END PARALLEL DO
    end if
   else
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(ghc,gh1,npw1,nspinor,gresid)
    do ipw=1,npw1*nspinor
     gresid(1,ipw)=-ghc(1,ipw)-gh1(1,ipw)
     gresid(2,ipw)=-ghc(2,ipw)-gh1(2,ipw)
    end do
!   $OMP END PARALLEL DO
   end if

!  ======================================================================
!  =========== PROJECT THE STEEPEST DESCENT DIRECTION ===================
!  ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
!  ======================================================================
!  Project all bands from gresid into direc:
!  The following projection over the subspace orthogonal to occupied bands
!  is not optional in the RF case, unlike the GS case.
!  However, the order of operations could be changed, so that
!  as to make it only applied at the beginning, to H(1) psi(0),
!  so, THIS IS TO BE REEXAMINED
!  Note the subtleties:
!  -For the generalized eigenPb, S|cgq> is used in place of |cgq>,
!  in order to apply P_c+ projector (see PRB 73, 235101 (2006), Eq. (71), (72)
!  -For the simple eigenPb, gscq is used as dummy argument
   if (prtvol<0) call status(iline,filstat,iexit,level,'projbd(1)     ')
   if (gen_eigenpb) then
    call projbd(gscq,gresid,-1,igscq,icgq,istwf_k,mgscq,mpi_enreg,mcgq,nband,npw1,nspinor,&
&    ortalg,iprint,cgq,scprod,tim_projbd,useoverlap)
   else
    call projbd(cgq,gresid,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&    ortalg,iprint,gscq,scprod,tim_projbd,useoverlap)
   end if
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(direc,npw1,nspinor,gresid)
   do ipw=1,npw1*nspinor
    direc(1,ipw)=gresid(1,ipw)
    direc(2,ipw)=gresid(2,ipw)
   end do
!  $OMP END PARALLEL DO

!  ======================================================================
!  ============== CHECK FOR CONVERGENCE CRITERIA ========================
!  ======================================================================

!  Compute second-order derivative of the energy using a variational expression
   call dotprod_g(prod1,doti,istwf_k,mpi_enreg,npw1*nspinor,1,cwavef,gresid)
   call dotprod_g(prod2,doti,istwf_k,mpi_enreg,npw1*nspinor,1,cwavef,gh1)
   d2te=two*(-prod1+prod2)
!  DEBUG
!  write(6,*)' cgwf3: prod1,prod2=',prod1,prod2
!  ENDDEBUG

!  Compute <u_m(1)|H(0)-e_m(0)|u_m(1)>
!  (<u_m(1)|H(0)-e_m(0).S|u_m(1)> if gen. eigenPb),
!  that should be positive,
!  except when the eigenvalue eig_mk(0) is higher than
!  the lowest non-treated eig_mk+q(0). For insulators, this
!  has no influence, but for metallic occupations,
!  the conjugate gradient algorithm breaks down. The solution adopted here
!  is very crude, and rely upon the fact that occupancies of such
!  levels should be smaller and smaller with increasing nband, so that
!  a convergence study will give the right result.
!  The same trick is also used later.
   u1h0me0u1=-prod1-prod2
   if(u1h0me0u1<zero)then
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(cwavef,ghc,npw1,nspinor)
    do ipw=1,npw1*nspinor
     cwavef(1,ipw)=zero
     cwavef(2,ipw)=zero
     ghc(1,ipw)  =zero
     ghc(2,ipw)  =zero
    end do
!   $OMP END PARALLEL DO
    if (use_vnl==1) then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(gvnlc,npw1,nspinor)
     do ipw=1,npw1*nspinor
      gvnlc(1,ipw)=zero
      gvnlc(2,ipw)=zero
     end do
!    $OMP END PARALLEL DO
    end if
!   A negative residual will be the signal of this problem ...
    resid=-one
    write(message, '(a)' )&
&    ' cgwf3: problem of minimisation (likely metallic), set resid to -1'
    call wrtout(06,message,'PERS')
!   Number of one-way 3D ffts skipped
    nskip=nskip+(nline-iline+1)
!   Exit from the loop on iline
    exit
   end if

!  Compute residual (squared) norm
   call sqnorm_g(resid,istwf_k,mpi_enreg,npw1*nspinor,gresid)
   if (prtvol==-level)then
    write(message,'(a,a,i3,f14.6,a,a,4es12.4)') ch10,&
&    ' cgwf3 : iline,eshift     =',iline,eshift,ch10,&
&    '         resid,prod1,prod2,d2te=',resid,prod1,prod2,d2te
    call wrtout(06,message,'PERS')
   end if

!  If residual sufficiently small stop line minimizations
   if (resid<tolwfr) then
    if(prtvol>=10)then
     write(message, '(a,i4,a,i2,a,es12.4)' ) &
&     ' cgwf3: band',band,' converged after ',iline,&
&     ' line minimizations : resid =',resid
     call wrtout(06,message,'PERS')
    end if
!   Number of two-way 3D ffts skipped
    nskip=nskip+(nline-iline+1)
!   Exit from the loop on iline
    exit
   end if

!  If user require exiting the job, stop line minimisations
   if (quit==1) then
    write(message, '(a,i4)' ) &
&    ' cgwf3: user require exiting => skip update of band ',band
    call wrtout(06,message,'PERS')
!   Number of two-way 3D ffts skipped
    nskip=nskip+(nline-iline+1)
!   Exit from the loop on iline
    exit
   end if

!  Check that d2te is decreasing on succeeding lines:
   if (iline/=1) then
    if (d2te>d2teold+1.d-12) then
     write(message, '(a,a,a,i8,a,1p,e14.6,a1,3x,a,1p,e14.6,a1)')&
&     ' cgwf3: WARNING -',ch10,&
&     '  New trial energy at line',iline,' = ',d2te,ch10,&
&     '  is higher than former:',d2teold,ch10
     call wrtout(06,message,'PERS')
    end if
   end if
   d2teold=d2te

!  DEBUG Keep this debugging feature !
!  call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,direc)
!  write(6,*)' cgwf3 : before precon, direc**2=',dotr
!  call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,cwaveq)
!  write(6,*)' cgwf3 : before precon, cwaveq**2=',dotr
!  ENDDEBUG

!  ======================================================================
!  ======== PRECONDITION THE STEEPEST DESCENT DIRECTION =================
!  ======================================================================

!  If wfoptalg>=10, the precondition matrix is kept constant
!  during iteration ; otherwise it is recomputed
   if (wfoptalg<10.or.iline==1) then
    if (prtvol<0) call status(iline,filstat,iexit,level,'call precon   ')
    call precon(cwaveq,zero,istwf_k,kinpw1,mpi_enreg,npw1,nspinor,0,pcon,direc)
   else
    do ispinor=1,nspinor
     igs=(ispinor-1)*npw1
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(igs,npw,direc,pcon)
     do ipw=1+igs,npw1+igs
      direc(1,ipw)=direc(1,ipw)*pcon(ipw-igs)
      direc(2,ipw)=direc(2,ipw)*pcon(ipw-igs)
     end do
!    $OMP END PARALLEL DO
    end do
   end if

!  DEBUG Keep this debugging feature !
!  call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,direc)
!  write(6,*)' cgwf3 : after precon, direc**2=',dotr
!  ENDDEBUG

!  ======================================================================
!  ======= PROJECT THE PRECOND. STEEPEST DESCENT DIRECTION ==============
!  ========= OVER THE SUBSPACE ORTHOGONAL TO OTHER BANDS ================
!  ======================================================================

!  Projecting again out all bands:
!  -For the simple eigenPb, gscq is used as dummy argument
   if (prtvol<0) call status(0,filstat,iexit,level,'prjbd(2)      ')
   call projbd(cgq,direc,-1,icgq,igscq,istwf_k,mcgq,mpi_enreg,mgscq,nband,npw1,nspinor,&
&   ortalg,iprint,gscq,scprod,tim_projbd,useoverlap)

!  DEBUG Keep this debugging feature !
!  call sqnorm_g(dotr,istwf_k,mpi_enreg,npw1*nspinor,direc)
!  write(6,*)' cgwf3 : after projbd, direc**2=',dotr
!  ENDDEBUG

!  ======================================================================
!  ================= COMPUTE THE CONJUGATE-GRADIENT =====================
!  ======================================================================

!  get dot of direction vector with residual vector
   call dotprod_g(dotgg,doti,istwf_k,mpi_enreg,npw1*nspinor,1,direc,gresid)

!  At first iteration, gamma is set to zero
   if (iline==1) then
    gamma=zero
    dotgp=dotgg
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(conjgr,direc,npw1,nspinor)
    do ipw=1,npw1*nspinor
     conjgr(1,ipw)=direc(1,ipw)
     conjgr(2,ipw)=direc(2,ipw)
    end do
!   $OMP END PARALLEL DO
   else
!   At next iterations, h = g + gamma * h
    gamma=dotgg/dotgp
    dotgp=dotgg
    if (prtvol==-level)then
     write(message,'(a,2es16.6)') 'cgwf3: dotgg,gamma =',dotgg,gamma
     call wrtout(06,message,'PERS')
    end if
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(conjgr,direc,gamma,npw1,nspinor)
    do ipw=1,npw1*nspinor
     conjgr(1,ipw)=direc(1,ipw)+gamma*conjgr(1,ipw)
     conjgr(2,ipw)=direc(2,ipw)+gamma*conjgr(2,ipw)
    end do
!   $OMP END PARALLEL DO
    if (prtvol==-level)then
     write(message,'(a)') &
&     'cgwf3: conjugate direction has been found'
     call wrtout(06,message,'PERS')
    end if
   end if

!  ======================================================================
!  ===== COMPUTE CONTRIBUTIONS TO 1ST AND 2ND DERIVATIVES OF ENERGY =====
!  ======================================================================
!  ...along the search direction

!  Compute dedt, Eq.(29) of of PRB55, 10337 (1997),
!  with an additional factor of 2 for the difference
!  between E(2) and the 2DTE
   call dotprod_g(dedt,doti,istwf_k,mpi_enreg,npw1*nspinor,1,conjgr,gresid)
   dedt=-two*two*dedt

   allocate(gvnl_direc(2,npw1*nspinor),gh_direc(2,npw1*nspinor))
   if (prtvol<0) call status(iline,filstat,iexit,level,'call getghc(2)')
   sij_opt=0;if (gen_eigenpb) sij_opt=-1
   call getghc(conjgr,dimffnl1,ffnl1,filstat,&
&   gh_direc,gsc_dummy,gs_hamkq,gvnl_direc,kg1_k,kinpw1,eshift,lmnmax,matblk,&
&   mgfft,mpi_enreg,mpsang,mpssoang,natom,1,npw1,nspinor,ntypat,nvloc,n4,n5,n6,&
&   paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
!  ghc also includes the eigenvalue shift
!  (already included if generalized eigenPb - see sijopt=-1 above)
   if (.not.gen_eigenpb) then
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(conjgr,eshift,gh_direc,npw1,nspinor)
    do ipw=1,npw1*nspinor
     gh_direc(1,ipw)=gh_direc(1,ipw)-eshift*conjgr(1,ipw)
     gh_direc(2,ipw)=gh_direc(2,ipw)-eshift*conjgr(2,ipw)
    end do
!   $OMP END PARALLEL DO
   end if
   if (prtvol<0) call status(iline,filstat,iexit,level,'after getghc(2')

!  compute d2edt2, Eq.(30) of of PRB55, 10337 (1997),
!  with an additional factor of 2 for the difference
!  between E(2) and the 2DTE, and neglect of local fields (SC terms)
   call dotprod_g(d2edt2,doti,istwf_k,mpi_enreg,npw1*nspinor,1,conjgr,gh_direc)
   d2edt2=two*two*d2edt2
   if(prtvol==-level)then
    write(message,'(a,2es14.6)') 'cgwf3: dedt,d2edt2=',dedt,d2edt2
    call wrtout(06,message,'PERS')
   end if

!  ======================================================================
!  ======= COMPUTE MIXING FACTOR - CHECK FOR CONVERGENCE ===============
!  ======================================================================

!  see Eq.(31) of PRB55, 10337 (1997)
   if(d2edt2<=zero)then
!   This may happen when the eigenvalue eig_mk(0) is higher than
!   the lowest non-treated eig_mk+q(0). The solution adopted here
!   is very crude, and rely upon the fact that occupancies of such
!   levels should be smaller and smaller with increasing nband, so that
!   a convergence study will give the right result.
    theta=zero
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(cwavef,ghc,npw1,nspinor)
    do ipw=1,npw1*nspinor
     cwavef(1,ipw)=zero
     cwavef(2,ipw)=zero
     ghc(1,ipw)  =zero
     ghc(2,ipw)  =zero
    end do
!   $OMP END PARALLEL DO
    if (use_vnl==1) then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(gvnlc,npw1,nspinor)
     do ipw=1,npw1*nspinor
      gvnlc(1,ipw)=zero
      gvnlc(2,ipw)=zero
     end do
!    $OMP END PARALLEL DO
    end if
!   A negative residual will be the signal of this problem ...
    resid=-two
    write(message, '(a)' )&
&    ' cgwf3: problem of minimisation (likely metallic), set resid to -2'
    call wrtout(06,message,'PERS')
   else
!   Here, the value of theta that gives the minimum
    theta=-dedt/d2edt2
!   DEBUG
!   write(6,*)' cgwf3: dedt,d2edt2=',dedt,d2edt2
!   ENDDEBUG
   end if

!  Check that result is above machine precision
   if (one+theta==one) then
    write(message, '(a,es16.4)' ) ' cgwf3: converged with theta=',theta
    call wrtout(06,message,'PERS')
!   Number of one-way 3D ffts skipped
    nskip=nskip+2*(nline-iline)
!   Exit from the loop on iline
    exit
   end if

!  ======================================================================
!  ================ GENERATE NEW |wf>, H|wf>, Vnl|Wf ... ================
!  ======================================================================

!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(conjgr,cwavef,ghc,gh_direc,npw1,nspinor,theta)
   do ipw=1,npw1*nspinor
    cwavef(1,ipw)=cwavef(1,ipw)+conjgr(1,ipw)*theta
    cwavef(2,ipw)=cwavef(2,ipw)+conjgr(2,ipw)*theta
    ghc(1,ipw)   =ghc(1,ipw)   +gh_direc(1,ipw)*theta
    ghc(2,ipw)   =ghc(2,ipw)   +gh_direc(2,ipw)*theta
   end do
!  $OMP END PARALLEL DO
   if (use_vnl==1) then
!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(gvnlc,gvnl_direc,npw1,nspinor,theta)
    do ipw=1,npw1*nspinor
     gvnlc(1,ipw) =gvnlc(1,ipw) +gvnl_direc(1,ipw)*theta
     gvnlc(2,ipw) =gvnlc(2,ipw) +gvnl_direc(2,ipw)*theta
    end do
!   $OMP END PARALLEL DO
   end if
   deallocate(gh_direc,gvnl_direc)

!  ======================================================================
!  =========== CHECK CONVERGENCE AGAINST TRIAL ENERGY ===================
!  ======================================================================

!  Check reduction in trial energy deltae, Eq.(28) of PRB55, 10337 (1997)
   deltae=half*d2edt2*theta**2+theta*dedt

   if (iline==1) then
    deold=deltae
!   Use a lower value for comparison with RESPFN
!   cgwftol=tol14
    cgwftol=0.01_dp   ! Value used until v3.3
   else if (abs(deltae)<cgwftol*abs(deold) .and. iline/=nline ) then
!   else if (abs(deltae)<0.005_dp*abs(deold) .and. iline/=nline ) then
    if(prtvol>=10)then
     write(message, '(a,i4,1x,a,1p,e12.4,a,e12.4,a)' ) &
&     ' cgwf3: line',iline,&
&     ' deltae=',deltae,' < cgwftol*',deold,' =>skip lines'
     call wrtout(06,message,'PERS')
    end if
!   Number of one-way 3D ffts skipped
    nskip=nskip+2*(nline-iline)
!   Exit from the loop on iline
    exit
   end if

!  ======================================================================
!  ================== END LOOP FOR GIVEN BAND ===========================
!  ======================================================================

!  Note that there are three "exit" instruction inside the loop.
  end do ! iline

  deallocate(conjgr,cwaveq,direc,gresid)

! End condition of not being a buffer band
 end if

 if (prtvol<0) call status(0,filstat,iexit,level,'after iline   ')

!At the end of the treatment of a set of bands, write the number
!of one-way 3D ffts skipped
#if !defined MPI
 if(band==nband .and. prtvol>=10)then
  write(message, '(a,i8)' )&
&  ' cgwf3: number of one-way 3D ffts skipped in cgwf3 until now =',nskip
  call wrtout(06,message,'PERS')
 end if
#endif

 deallocate(gh1,pcon,scprod)
 if (berryopt==4) deallocate(gberry)

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i2,a)') ch10,&
&  ' cgwf3 : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 call timab(122,2,tsec)
 if(prtvol<0) call status(0,filstat,iexit,level,'exit          ')

!DEBUG
!write(6,*)' cgwf3 : exit'
!ENDDEBUG

end subroutine cgwf3
!!***
