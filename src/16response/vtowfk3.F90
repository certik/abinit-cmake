!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtowfk3
!! NAME
!! vtowfk3
!!
!! FUNCTION
!! This routine compute the partial density at a given k-point,
!! for a given spin-polarization, from a fixed potential (vlocal1).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, AR, DRH, MB, MVer,XW, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cgq(2,mcgq)=array for planewave coefficients of wavefunctions.
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of
!!    RF wavefunctions at k,q.
!!  cplex=1 if rhoaug1 is real, 2 if rhoaug1 is complex
!!  cprj(dimpaw1,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!               projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  cpus= cpu time limit in seconds
!!  dimcprj(natom)=array of dimensions of arrays cprj, cprjq
!!  dimekb=first dimension of ekb (see ekb_typ)
!!  dimffnlk=second dimension of ffnlk (1+number of derivatives)
!!  dimffnl1=second dimension of ffnl1 and ffnlkq (1+number of derivatives)
!!  dimpaw1= -PAW only- dimension of 1st-order on-site quantities (rhoij, Dij...)
!!  dkinpw(npw_k)=derivative of the (modified) kinetic energy for
!!    each plane wave at k (Hartree)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig0_k(nband_k)=GS eigenvalues at k (hartree)
!!  eig0_kq(nband_k)=GS eigenvalues at k+Q (hartree)
!!  ekb_typ(dimekb,1,nspinor**2)=
!!    ->Norm conserving : (Real) Kleinman-Bylander energies (hartree)
!!                        for the displaced atom
!!                        for number of basis functions (l,n) (lnmax)
!!                        dimekb=lnmax
!!    ->PAW : (Real, symmetric) Frozen part of Dij coefficients
!!                        to connect projectors
!!                        for the displaced atom
!!                        for number of basis functions (l,m,n) (lmnmax)
!!                        dimekb=lmnmax*(lmnmax+1)/2
!!  ekb1_typ(dimekb,1,useekb1*nspinor**2)=1st derivative of ekb_typ for the current pertubation (see above)
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  ffnlk(npw_k,dimffnlk,lmnmax,1)=nonloc form factors at k, for the displaced atom.
!!  ffnlkq(npw1_k,dimffnl1,lmnmax,1)=nonloc form fact at k+q for the displaced atom
!!  ffnl1(npw1_k,dimffnl1,lmnmax,ntypat)=nonloc form factors at k+q
!!  gbound(2*mgfft+8,2)=G sphere boundary
!!  grad_berry(2,mpw1,dtefield%nband_occ) = the gradient of the Berry phase term
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  ibg1=shift to be applied on the location of data in the array cprj1
!!  icg=shift to be applied on the location of data in the array cg
!!  icgq=shift to be applied on the location of data in the array cgq
!!  icg1=shift to be applied on the location of data in the array cg1
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  indlmn_typ(6,lmnmax,1)=indlmn info for the displaced atom
!!  ipert=type of the perturbation
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istep=index of the number of steps in the routine scfcv
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kinpw1(npw1_k)=(modified) kinetic energy for each plane wave at k+q (Hartree)
!!  kpg_k(npw_k,nkpg)= (k+G) components at k (only if useylm=1)
!!  kpg1_k(npw1_k,nkpg1)= (k+G) components at k+q (only if useylm=1)
!!  kpt(3)=reduced coordinates of k points.
!!  lmnmax= max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mband=maximum number of bands
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mkqmem =number of k+q points which can fit in memory (GS data); 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nkpg,nkpg1=second dimensions of kpg_k and kpg1_k (0 if useylm=0)
!!  nkpt=number of k points
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!    (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  n4,n5,n6 used for dimensioning real space arrays
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pspso_typ(1)=spin-orbit info for the displaced atom
!!  rhoaug1(cplex*n4,n5,n6)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output)
!!  rocceig(nband_k,nband_k)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!    if this
!!   ratio has been attributed to the band n (second argument), zero otherwise
!!  sij_typ(dimekb,usepaw)=-PAW only- overlap matrix components for the current perturbation
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  useekb1=1 if ekb derivatives (ekb1) exist
!!  wffddk=struct info for wf ddk file.
!!  wffnew=struct info for OUTPUT 1st-order wf file
!!  wffnow=struct info for INPUT 1st-order wf file
!!  wfftgs=struct info for GS wf disk files.
!!  vlocal(n4,n5,n6)= GS local potential in real space, on the augmented
!!    fft grid
!!  vlocal1(cplex*n4,n5,n6)= RF local pot. in real space, on the augm. fft grid
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q.
!!  edocc_k(nband_k)=correction to 2nd-order total energy coming
!!      from changes of occupation
!!  eeig0_k(nband_k)=zero-order eigenvalues contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  eig1_k(2*nband_k**2)=first-order eigenvalues (hartree)
!!  ek0_k(nband_k)=0-order kinetic energy contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  ek1_k(nband_k)=1st-order kinetic energy contribution to 2nd-order total
!!      energy from all bands at this k point.
!!  eloc0_k(nband_k)=zero-order local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  enl0_k(nband_k)=zero-order non-local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  enl1_k(nband_k)=first-order non-local contribution to 2nd-order total energy
!!      from all bands at this k point.
!!  resid_k(nband_k)=residuals for each band over all k points,
!!  rhoaug1(cplex*n4,n5,n6)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output).
!!  ==== if (gs_hamkq%usepaw==1) ====
!!    cprj1(dimpaw1,nspinor*mband*mk1mem*nsppol*usecprj)=
!!              1st-order wave functions at k+q projected with non-local projectors:
!!                          cprj1=<p_i|C1nk+q> where p_i is a non-local projector
!!    pawrhoij1(dimpaw1) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!                                              (cumulative, so input as well as output)
!!
!! TODO
!!
!! PARENTS
!!      vtorho3
!!
!! CHILDREN
!!      cgwf3,chkexi,cprj_get,cprj_put,cprj_alloc,cprj_copy,cprj_free
!!      dotprod_g,fourwf,getcprj,leave_new,matrixelmt_g,meanvalue_g,pawmkrhoij3
!!      sqnorm_g,status,timab,wffreaddatarec,wffreadnpwrec,wffreadskiprec
!!      wffwritedatarec,wffwritenpwrec,wrtout,xcomm_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine vtowfk3(cg,cgq,cg1,cplex,cprj,cprjq,cprj1,cpus, gh1_k,dimcprj,dimekb,dimffnlk,dimffnl1,dimpaw1,dkinpw,dtfil,dtset,&
& edocc_k,eeig0_k,eig0_k,eig0_kq,eig1_k,ekb_typ,ekb1_typ,ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k,&
& fermie1,ffnlk,ffnlkq,ffnl1,gbound,grad_berry,gs_hamkq,&
& ibg,ibgq,ibg1,icg,icgq,icg1,idir,ikpt,indlmn_typ,ipert,&
& isppol,istep,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lmnmax,matblk,mband,mcgq,mcprjq,mgfft,mkmem,mkqmem,mk1mem,&
& mpi_enreg,mpsang,mpssoang,mpw,mpw1,natom,nband_k,nband_kq,&
& nkpg,nkpg1,nkpt,nnsclo_now,npw_k,npw1_k,nspden,nspinor,nsppol,&
& ntypat,n4,n5,n6,occ_k,pawrhoij1,ph3d,prtvol,psps,pspso_typ,resid_k,rhoaug1,rocceig,&
& sij_typ,usecprj,useekb1,wffddk,wffnew,wffnow,wfftgs,vlocal,vlocal1,wtk_k)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12spacepar
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_16response, except_this_one => vtowfk3
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,dimekb,dimffnl1,dimffnlk,dimpaw1,ibg,ibg1,ibgq,icg
 integer,intent(in) :: icg1,icgq,idir,ikpt,ipert,isppol,istep,lmnmax,matblk
 integer,intent(in) :: mband,mcgq,mcprjq,mgfft,mk1mem,mkmem,mkqmem,mpsang
 integer,intent(in) :: mpssoang,mpw,mpw1,n4,n5,n6,natom,nkpg,nkpg1,nkpt
 integer,intent(in) :: nnsclo_now,nspden,nspinor,nsppol,ntypat,prtvol,usecprj
 integer,intent(in) :: useekb1
 integer,intent(inout) :: nband_k,nband_kq,npw1_k,npw_k
 real(dp),intent(in) :: cpus,fermie1,wtk_k
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffddk,wffnew,wffnow,wfftgs
!arrays
 integer,intent(in) :: dimcprj(natom),gbound(2*mgfft+8,2)
 integer,intent(in) :: indlmn_typ(6,lmnmax,1),kg1_k(3,npw1_k),kg_k(3,npw_k)
 integer,intent(in) :: pspso_typ(1)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),cgq(2,mcgq)
 real(dp),intent(in) :: dkinpw(npw_k),eig0_k(nband_k),eig0_kq(nband_k)
 real(dp),intent(in) :: ekb1_typ(dimekb,1,useekb1*nspinor**2)
 real(dp),intent(in) :: ekb_typ(dimekb,1,nspinor**2)
 real(dp),intent(in) :: ffnl1(npw1_k,dimffnl1,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlk(npw_k,dimffnlk,lmnmax,1)
 real(dp),intent(in) :: ffnlkq(npw1_k,dimffnl1,lmnmax,1)
 real(dp),intent(in) :: grad_berry(2,mpw1,nband_k),kinpw1(npw1_k)
 real(dp),intent(in) :: kpg1_k(npw1_k,nkpg1),kpg_k(npw_k,nkpg),kpt(3)
 real(dp),intent(in) :: occ_k(nband_k),rocceig(nband_k,nband_k)
 real(dp),intent(inout) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(inout) :: ph3d(2,npw1_k,matblk),rhoaug1(cplex*n4,n5,n6)
 real(dp),intent(inout) :: sij_typ(dimekb,gs_hamkq%usepaw),vlocal(n4,n5,n6)
 real(dp),intent(inout) :: vlocal1(cplex*n4,n5,n6)
 real(dp),intent(out) :: edocc_k(nband_k),eeig0_k(nband_k),eig1_k(2*nband_k**2)
 real(dp),intent(out) :: ek0_k(nband_k),ek1_k(nband_k),eloc0_k(nband_k)
 real(dp),intent(out) :: enl0_k(nband_k),enl1_k(nband_k)
 real(dp),intent(out) :: gh1_k(nband_k,2,mpw1*nspinor),resid_k(nband_k)
 type(cprj_type),intent(in) :: cprj(dimpaw1,nspinor*mband*mkmem*nsppol*usecprj)
 type(cprj_type),intent(in) :: cprjq(natom,mcprjq)
 type(cprj_type),intent(out) :: cprj1(dimpaw1,nspinor*mband*mk1mem*nsppol*usecprj)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(dimpaw1)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=18
 integer,save :: count=0,nskip=0
 integer :: accesswff,choice,counter,cpopt,i1,i2,i3,iband,ibandkq,ieig2rf,ier
 integer :: ierr,iexit,ig,igs,igscq,ii,index_cgq,index_cprj,index_eig1
 integer :: index_gscq,inonsc,iorder_cprj,iowf,iproj,ipsang,ipw,ipw1,iscf_mod
 integer :: ispinor,istwf_k,iwavef,mcgnpw,mcgnpw1,me,mgscq,n1,n2,n3,nkpt_max
 integer :: nrecwf,nspinor0,openexit,option_rhoij,paw_opt,quit,signs,spaceComm
 integer :: tag,test_ddk,tim_fourwf,tim_nonlop,tim_rwwf,tocceig,usedcwavef
 real(dp) :: aa,ai,ar,dum,eig0nk,facti,factr,im0,im1,invocc,re0,re1,resid
 real(dp) :: residk,scprod,valuei,valuer,weight
 character(len=500) :: message
!arrays
 integer :: nattyp_atm(1)
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: dummy(2,1),dummy1(1),qphon(3),tsec(2)
 real(dp),allocatable :: cgddk(:,:),cgnew(:,:),cgnow(:,:),cgtgs(:,:)
 real(dp),allocatable :: cwave0(:,:),cwave1(:,:),cwavef(:,:),cwavef_sp(:,:)
 real(dp),allocatable :: dcwavef(:,:),eig2nkq(:),eig_dum(:),gh1_n(:,:),ghc(:,:)
 real(dp),allocatable :: grnk(:),gsc(:,:),gscq(:,:),gvnl1(:,:),gvnlc(:,:)
 real(dp),allocatable :: occ_dum(:),ph1d_atm(:,:),rhoaug(:,:,:),wfraug(:,:,:,:)
 real(dp),allocatable :: wfraug1(:,:,:,:)
 type(cprj_type),allocatable :: cwaveprj(:,:),cwaveprj0(:,:),cwaveprj1(:,:)
 type(cprj_type),allocatable :: cwaveprj_tmp(:,:)

! *********************************************************************

!Keep track of total time spent in vtowfk3
 call timab(128,1,tsec)

 nkpt_max=50
 if(mpi_enreg%paral_compil_kpt==1)nkpt_max=-1

!DEBUG
!write(6,*)' vtowfk3: enter '
!write(6,*)' vtowfk3: ikpt=',ikpt
!count=count+1
!write(6,*)' count=',count
!if(count==27)stop
!write(6,*)' vtowfk3: prtvol,wtk_k,npw_k,npw1_k,ipert'
!write(6,*)prtvol,wtk_k,npw_k,npw1_k,ipert
!if(ikpt==4)stop
!write(6,*)' vtowfk3 : cg1(:,1)=',cg1(:,1)
!write(6,*)' nband_k,natom,npw_k',nband_k,natom,npw_k
!stop
!ENDDEBUG

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,'vtowfk3 : enter'
  call wrtout(06,message,'PERS')
 end if

 quit=0
!Init me
 call xme_init(mpi_enreg,me)
!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
 accesswff=dtset%accesswff

 n1=gs_hamkq%ngfft(1) ; n2=gs_hamkq%ngfft(2) ; n3=gs_hamkq%ngfft(3)
 qphon(1:3)=dtset%qptn(1:3)

 iscf_mod=dtset%iscf
 istwf_k=gs_hamkq%istwf_k

!The value of iscf must be modified if ddk perturbation, see loper3.f
 if(ipert==natom+1) iscf_mod=-3

 allocate(ghc(2,npw1_k*nspinor))
 if (gs_hamkq%usepaw==0) allocate(gvnlc(2,npw1_k*nspinor),gvnl1(2,npw1_k*nspinor))

 if(prtvol>2 .or. ikpt<=nkpt_max)then
  write(message, '(a,a,i5,2x,a,3f9.5,2x,a)' ) ch10,&
&  ' Non-SCF iterations; k pt #',ikpt,'k=',kpt(:),'band residuals:'
  call wrtout(06,message,'PERS')
 end if

 allocate(wfraug(2,n4,n5,n6),wfraug1(2,n4,n5,n6))
 allocate(rhoaug(n4,n5,n6))
 allocate(cwave0(2,npw_k*nspinor),cwavef(2,npw1_k*nspinor))
 allocate(cwave1(2,npw1_k*nspinor))
 allocate(gh1_n(2,npw1_k*nspinor))
!Read the npw and kg records of wf files
!NOTE : it should be possible to use rwwf in the present routine
 call status(0,dtfil%filstat,iexit,level,'before WffRead')
 test_ddk=0
 if( ipert==natom+2 .and. &
& sum( (dtset%qptn(1:3))**2 ) < 1.0d-7 .and. (dtset%berryopt .ne. 4) )then
  test_ddk=1
! Read npw record
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor0,wffddk)
! Skip k+G record
  call WffReadSkipRec(ierr,1,wffddk)
 end if
 if( ipert==natom+5 .and. &
& sum( (dtset%qptn(1:3))**2 ) < 1.0d-7 .and. (dtset%berryopt .ne. 4) )then
  test_ddk=1
! Read npw record
! JWZ, 20-Aug-08
! original code: call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor,wffddk)
! problem here is that nspinor is intent(in), but in WffReadNpwRec nspinor is intent(out).
! I think this should read nspinor0, like the others
! 
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor0,wffddk)
! 
! Skip k+G record
  call WffReadSkipRec(ierr,1,wffddk)
 end if
 if(mkmem==0)then
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_k,nspinor0,wfftgs)
! Skip k+G and eigenvalue records in wfftgs (already in eigen0)
  call WffReadSkipRec(ierr,2,wfftgs)
 end if
 if(mk1mem==0)then
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor0,wffnow)
! Skip k+G record
  call WffReadSkipRec(ierr,1,wffnow)
! Initialize writing for this k point
  call WffWriteNpwRec(ierr,nband_k,npw1_k,nspinor,wffnew)
  allocate(kg_dum(3,npw1_k))
  kg_dum(:,:) = kg1_k(:,:)
  call WffWriteDataRecInt(kg_dum,ierr,3*npw1_k,wffnew)
  deallocate(kg_dum)
 end if

!Additional arrays for PAW
 if (gs_hamkq%usepaw==1) then
! 1-Compute all <g|S|Cnk+q>
  igscq=0;mgscq=mpw1*nspinor*mband
  allocate(gscq(2,mgscq),gsc(2,npw1_k*nspinor))
  if (usecprj==1) then
   allocate(cwaveprj_tmp(natom,nspinor))
   call cprj_alloc(cwaveprj_tmp,cprj(1,1)%ncpgr,dimcprj)
  end if
  index_cprj=ibgq;index_cgq=icgq;index_gscq=igscq
  do ibandkq=1,nband_k
   if (mpi_enreg%paral_compil_kpt==1)then
    if (mpi_enreg%proc_distrb(ikpt,ibandkq,isppol)/=me) then
     gscq(:,1+index_gscq:npw1_k*nspinor+index_gscq)=zero
     index_cprj=index_cprj+nspinor
     index_cgq=index_cgq+npw1_k*nspinor
     index_gscq=index_gscq+npw1_k*nspinor
     cycle
    end if
   end if
   cwave0(:,1:npw1_k*nspinor)=cgq(:,1+index_cgq:npw1_k*nspinor+index_cgq)
   if (usecprj==1) then
    call cprj_copy(cprjq(:,1+index_cprj:nspinor+index_cprj),cwaveprj_tmp)
   end if
   choice=1 ; signs=2 ; cpopt=-1+3*usecprj ; paw_opt=3 ; tim_nonlop=0
   call nonlop(gs_hamkq%atindx1,choice,cpopt,cwaveprj_tmp,gs_hamkq%dimekb1,gs_hamkq%dimekb2,dimffnl1,dimffnl1,&
&   gs_hamkq%ekb,dummy1,ffnlkq,ffnlkq,gs_hamkq%gmet,gs_hamkq%gprimd,0,gs_hamkq%indlmn,&
&   gs_hamkq%istwf_k,kg1_k,kg1_k,kpg1_k,kpg1_k,gs_hamkq%kpoint,gs_hamkq%kpoint,dum,lmnmax,&
&   matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&   nkpg1,nkpg1,gs_hamkq%nloalg,1,npw1_k,npw1_k,nspinor,ntypat,0,paw_opt,gs_hamkq%phkxred,&
&   gs_hamkq%phkxred,gs_hamkq%ph1d,ph3d,ph3d,gs_hamkq%pspso,signs,gs_hamkq%sij,&
&   gsc,tim_nonlop,gs_hamkq%ucvol,gs_hamkq%useylm,cwave0,gvnlc)
   gscq(:,1+index_gscq:npw1_k*nspinor+index_gscq)=gsc(:,1:npw1_k*nspinor)
   index_cprj=index_cprj+nspinor
   index_cgq=index_cgq+npw1_k*nspinor
   index_gscq=index_gscq+npw1_k*nspinor
  end do
  deallocate(gsc)
  if (usecprj==1) then
   call cprj_free(cwaveprj_tmp)
   deallocate(cwaveprj_tmp)
  end if
  if ((mpi_enreg%paral_compil_kpt==1).and.(mpi_enreg%has_band_comm==1)) then
   call timab(48,1,tsec)
   call xsum_mpi(gscq,mpi_enreg%band_comm(ikpt+(isppol-1)*nkpt),ierr)
   call timab(48,2,tsec)
  end if
! 2-Initialize additional scalars/arrays
  option_rhoij=2;iorder_cprj=0
  allocate(dcwavef(2,npw_k*nspinor))
  if (usecprj==1) then
   allocate(cwaveprj0(dimpaw1,nspinor))
   call cprj_alloc(cwaveprj0,cprj(1,1)%ncpgr,dimcprj)
  end if
  allocate(cwaveprj(dimpaw1,nspinor),cwaveprj1(dimpaw1,nspinor))
  call cprj_alloc(cwaveprj ,cprj1(1,1)%ncpgr,dimcprj)
  call cprj_alloc(cwaveprj1,cprj1(1,1)%ncpgr,dimcprj)
 else
  igscq=0;mgscq=0
 end if

 call timab(139,1,tsec)

!Loop over bands

 do iband=1,nband_k

  if(mpi_enreg%paral_compil_kpt==1)then

   if( (mpi_enreg%proc_distrb(ikpt, iband,isppol) /= me )   ) then

    if(test_ddk==1)then
!    Skip the eigenvalue and the wf records of this band
     call WffReadSkipRec(ierr,2,wffddk)
    end if
    if(mkmem==0)then
     call WffReadSkipRec(ierr,1,wfftgs)
    end if
    if(mk1mem==0)then
     call WffReadSkipRec(ierr,2,wffnow)
!    Fill these records with zeroes (so that they can be read without I/O error)
     call WffWriteDataRec( (/ (zero*dble(ii),ii=1,2*nband_k) /) ,ierr,2*nband_k,wffnew)
     call WffWriteDataRec( (/ (zero*dble(ii),ii=1,2*npw1_k*nspinor) /) ,ierr,2*npw1_k*nspinor,wffnew)
    end if
    cycle
   end if
  end if ! paral

! Read ground-state wavefunctions
  if(mkmem/=0)then
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cg,cwave0,iband,icg,npw_k,nspinor)
   do ipw=1,npw_k*nspinor
    cwave0(1,ipw)=cg(1,ipw+(iband-1)*npw_k*nspinor+icg)
    cwave0(2,ipw)=cg(2,ipw+(iband-1)*npw_k*nspinor+icg)
   end do
!  $OMP END PARALLEL DO
  else
   call timab(288,1,tsec)
   call WffReadDataRec(cwave0,ierr,2*npw_k*nspinor,wfftgs)
   call timab(288,2,tsec)
  end if
  if (gs_hamkq%usepaw==1.and.usecprj==1) then
   index_cprj=(iband-1)*nspinor+ibg
   call cprj_copy(cprj(:,index_cprj+1:index_cprj+nspinor),cwaveprj0)
  end if

! Read first-order wavefunctions
  if(mk1mem/=0)then
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cg1,cwavef,iband,icg1,npw1_k,nspinor)
   do ipw=1,npw1_k*nspinor
    cwavef(1,ipw)=cg1(1,ipw+(iband-1)*npw1_k*nspinor+icg1)
    cwavef(2,ipw)=cg1(2,ipw+(iband-1)*npw1_k*nspinor+icg1)
   end do
!  $OMP END PARALLEL DO
  else
   call timab(288,1,tsec)
!  Skip the eigenvalue line
   call WffReadSkipRec(ierr,1,wffnow)
   call WffReadDataRec(cwavef,ierr,2*npw1_k*nspinor,wffnow)
   call timab(288,2,tsec)
  end if
  if (gs_hamkq%usepaw==1.and.usecprj==1) then
   index_cprj=(iband-1)*nspinor+ibg1
   call cprj_copy(cprj1(:,index_cprj+1:index_cprj+nspinor),cwaveprj)
  end if

! Filter the wavefunctions for large modified kinetic energy
! The GS wavefunctions should already be non-zero
  do ispinor=1,nspinor
   igs=(ispinor-1)*npw1_k
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cwavef,igs,kinpw1,npw1_k)
   do ipw=1+igs,npw1_k+igs
    if(kinpw1(ipw-igs)>huge(zero)*1.d-11)then
     cwavef(1,ipw)=zero
     cwavef(2,ipw)=zero
    end if
   end do
!  $OMP END PARALLEL DO
  end do

  if(prtvol>=10)then
   call status(0,dtfil%filstat,iexit,level,'after wf read ')
  end if

! If electric field, the derivative of the wf should be read,
! and multiplied by i.
  if(test_ddk==1)then
!  Skip the eigenvalue record
   call WffReadSkipRec(ierr,1,wffddk)
!  Read gvnl1
   call WffReadDataRec(gvnl1,ierr,2*npw1_k*nspinor,wffddk)
!  Multiplication by -i
!  MVeithen 021212 : use + i instead,
!  See X. Gonze, Phys. Rev. B 55, 10337 (1997) Eq. (79)
!  the operator used to compute the first-order derivative
!  of the wavefunctions with respect to an electric field
!  is $+i \frac{d}{dk}$
!  This change will affect the computation of the 2dtes from non
!  stationary expressions, see nstdy3.f and nstwf3.f

   do ipw=1,npw1_k*nspinor
!   aa=gvnl1(1,ipw)
!   gvnl1(1,ipw)=gvnl1(2,ipw)
!   gvnl1(2,ipw)=-aa
    aa=gvnl1(1,ipw)
    gvnl1(1,ipw)=-gvnl1(2,ipw)
    gvnl1(2,ipw)=aa
   end do
  end if

! Unlike in GS calculations, the inonsc loop is inside the band loop
! nnsclo_now=number of non-self-consistent loops for the current vtrial
! (often 1 for SCF calculation, =nstep for non-SCF calculations)
  do inonsc=1,nnsclo_now

   counter=100*iband+inonsc
!  Because in this loop, the CPU time matters, the writing
!  in the STATUS file is usually inhibited
   if(prtvol>=10)then
    call status(counter,dtfil%filstat,iexit,level,'loop iband    ')
   end if

!  Not too often, check whether the run must be stopped.
!  If so, iexit will be non-zero.
!  Note that when the number of bands becomes large, the check
!  must be done more often, because treating one band takes also longer ...
!  Only do this in the sequential mode
   if(mpi_enreg%paral_compil_kpt==0)then
    if(iband==1 .or. (nband_k>=16 .and. mod(iband,8)==1) &
&    .or. (nband_k>=32 .and. mod(iband,4)==1) &
&    .or. (nband_k>=64 .and. mod(iband,2)==1) &
&    .or. (nband_k>=128)                        )then
     openexit=1 ; if(dtset%chkexit<=1) openexit=0
     call chkexi(cpus,dtfil%filnam_ds(1),iexit,6,mpi_enreg,openexit)
     if(iexit/=0)quit=1
    end if
   end if

   if(prtvol>=10)then
    call status(counter,dtfil%filstat,iexit,level,'call cgwf3    ')
   end if

!  Note that the following translation occurs in the called routine :
!  iband->band, nband_k->nband, npw_k->npw, npw1_k->npw1
   eig0nk=eig0_k(iband)
   usedcwavef=gs_hamkq%usepaw;if (istep==1.and.inonsc==1) usedcwavef=2*usedcwavef
   call cgwf3(iband,dtset%berryopt,cgq,cplex,cwavef,cwave0,cwaveprj,cwaveprj0,dcwavef,gh1_n,dimekb,&
&   dimffnlk,dimffnl1,dkinpw,eig0nk,eig0_kq,eig1_k,&
&   ekb_typ,ekb1_typ,ffnlk,ffnlkq,ffnl1,dtfil%filstat,gbound,ghc,grad_berry,&
&   gscq,gs_hamkq,gvnlc,gvnl1,icgq,idir,indlmn_typ,ipert,igscq,&
&   kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lmnmax,matblk,mcgq,mgfft,mgscq,mpi_enreg,&
&   mpsang,mpssoang,mpw1,natom,nband_k,dtset%nbdbuf,dimpaw1,nkpg,nkpg1,dtset%nline,&
&   npw_k,npw1_k,nspinor,ntypat,n4,n5,n6,dtset%ortalg,dtset%paral_kgb,ph3d,prtvol,&
&   pspso_typ,qphon,quit,resid,dtset%sciss,sij_typ,dtset%tolwfr,usecprj,usedcwavef,useekb1,vlocal,&
&   vlocal1,dtset%wfoptalg,wfraug)

   resid_k(iband)=resid
   if(prtvol>=10)then
    call status(counter,dtfil%filstat,iexit,level,'after cgwf    ')
   end if

   gh1_k(iband,:,:)=zero
   do ipw=1,npw1_k*nspinor
    gh1_k(iband,:,ipw)= gh1_n(:,ipw)
   end do

!  At this stage, the 1st order function cwavef is orthogonal to cgq (unlike
!  when it is input to cgwf3). Here, restore the "active space" content
!  of the first-order wavefunction, to give cwave1 .
!  $OMP PARALLEL DO PRIVATE(ii) &
!  $OMP&SHARED(cwave1,cwavef,npw1_k,nspinor)
   do ii=1,npw1_k*nspinor
    cwave1(1,ii)=cwavef(1,ii)
    cwave1(2,ii)=cwavef(2,ii)
   end do
   if (gs_hamkq%usepaw==1) then
    call cprj_copy(cwaveprj,cwaveprj1)
   end if
!  $OMP END PARALLEL DO
   tocceig=0
   if ( abs(occ_k(iband)) > tol8 ) then
    invocc=1.0_dp/occ_k(iband)
    edocc_k(iband)=zero
    call timab(213,1,tsec)
    do ibandkq=1,nband_k
     if ( abs(rocceig(ibandkq,iband)) > tol8 ) then
      tocceig=1
      index_eig1=2*ibandkq-1+(iband-1)*2*nband_k
      index_cgq=npw1_k*nspinor*(ibandkq-1)+icgq
      factr= rocceig(ibandkq,iband)*invocc*eig1_k(index_eig1  )
      if(ibandkq==iband) then
       factr= rocceig(ibandkq,iband)*invocc*(eig1_k(index_eig1  )&
&       -fermie1)
      else
       factr= rocceig(ibandkq,iband)*invocc*eig1_k(index_eig1  )
      end if
      facti= rocceig(ibandkq,iband)*invocc*eig1_k(index_eig1+1)
!     $OMP PARALLEL DO PRIVATE(ii) &
!     $OMP&SHARED(cgq,cwave1,cwavef,facti,factr,index_cgq,npw1_k,nspinor)
      do ii=1,npw1_k*nspinor
       cwave1(1,ii)=cwave1(1,ii)+ &
&       ( factr*cgq(1,ii+index_cgq)-facti*cgq(2,ii+index_cgq) )
       cwave1(2,ii)=cwave1(2,ii)+ &
&       ( facti*cgq(1,ii+index_cgq)+factr*cgq(2,ii+index_cgq) )
      end do
!     $OMP END PARALLEL DO
!     The factor of two is needed because we compute the 2DTE, and not E(2)
      edocc_k(iband)=edocc_k(iband)-two* &
&      (factr*eig1_k(index_eig1)+facti*eig1_k(index_eig1+1))
     end if
    end do
    call timab(213,2,tsec)
   end if

   if ( abs(occ_k(iband)) <= tol8 ) then

    ek0_k(iband)=zero
    ek1_k(iband)=zero
    eeig0_k(iband)=zero
    enl0_k(iband)=zero
    enl1_k(iband)=zero
    eloc0_k(iband)=zero
    nskip=nskip+1

   else

!   Compute the 0-order kinetic operator contribution (with cwavef)
    call meanvalue_g(ar,kinpw1,0,istwf_k,mpi_enreg,npw1_k,nspinor,cwavef)
!   There is an additional factor of 2 with respect to the bare matrix element
    ek0_k(iband)=two*ar

!   Compute the 1-order kinetic operator contribution (with cwave1 and cwave0), if needed.
!   Note that this is called only for ddk or strain, so that npw1_k=npw_k
    if(ipert==natom+1 .or. ipert==natom+3 .or. ipert==natom+4)then
     call matrixelmt_g(ai,ar,dkinpw,istwf_k,mpi_enreg,0,npw_k,nspinor,cwave1,cwave0)
!    There is an additional factor of 4 with respect to the bare matrix element
     ek1_k(iband)=four*ar
    end if

!   Compute eigenvalue part of total energy (with cwavef)
    call sqnorm_g(scprod,istwf_k,mpi_enreg,npw1_k*nspinor,cwavef)
    eeig0_k(iband)=-two*(eig0_k(iband)- (dtset%sciss) )*scprod

!   Compute nonlocal psp contributions to nonlocal energy:
!   <G|Vnl|C1nk(perp)> is contained in gvnlc (with cwavef)
    call dotprod_g(scprod,ai,istwf_k,mpi_enreg,npw1_k*nspinor,1,cwavef,gvnlc)
    enl0_k(iband)=two*scprod

!   <G|Vnl1|Cnk> is contained in gvnl1 (with cwave1)
    call dotprod_g(scprod,ai,istwf_k,mpi_enreg,npw1_k*nspinor,1,cwave1,gvnl1)
    enl1_k(iband)=four*scprod

!   Removal of the 1st-order kinetic energy from the 1st-order non-local part.
    if(ipert==natom+1 .or. &
&    ipert==natom+3 .or. ipert==natom+4) then
     enl1_k(iband)=enl1_k(iband)-ek1_k(iband)
    end if

!   In this last part of the treatment of one band, one has to
!   perform Fourier transforms, and to treat separately the two
!   spinorial components of the wavefunction.

    valuer=zero  ! Will be accumulated to give the local potential energy contribution

    do ispinor=1,nspinor

     if(prtvol>=10)then
      call status(counter,dtfil%filstat,iexit,level,'density update')
     end if

!    Fourier transform of cwavef. Here, rhoaug1 is a dummy variable.
!    NOTE : should take into account nspinor
     tim_fourwf=5
     if(ispinor==1)then
      call fourwf(cplex,rhoaug1,cwavef,dummy,wfraug1,&
&      gs_hamkq%gbound,gs_hamkq%gbound,&
&      istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&      npw1_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
     else
      allocate(cwavef_sp(2,npw1_k))
!     $OMP PARALLEL DO PRIVATE(ipw) &
!     $OMP&SHARED(cwavef,cwavef_sp,npw1_k)
      do ipw=1,npw1_k
       cwavef_sp(1,ipw)=cwavef(1,ipw+npw1_k)
       cwavef_sp(2,ipw)=cwavef(2,ipw+npw1_k)
      end do
!     $OMP END PARALLEL DO
      call fourwf(cplex,rhoaug1,cwavef_sp,dummy,wfraug1,&
&      gs_hamkq%gbound,gs_hamkq%gbound,&
&      istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&      npw1_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
      deallocate(cwavef_sp)
     end if

     if(prtvol>=10)then
      call status(counter,dtfil%filstat,iexit,level,'get eloc0_k   ')
     end if

!    Compute contribution of this band to
!    zero-order potential part of the 2nd-order total energy
!    $OMP PARALLEL DO PRIVATE(i1,i2,i3) REDUCTION(+:valuer) &
!    $OMP&SHARED(n1,n2,n3,vlocal,wfraug1)
     do i3=1,n3
      do i2=1,n2
       do i1=1,n1
        valuer=valuer+vlocal(i1,i2,i3)* &
&        (wfraug1(1,i1,i2,i3)**2+wfraug1(2,i1,i2,i3)**2)
       end do
      end do
     end do
!    $OMP END PARALLEL DO

!    Compute contribution to density only at the last inonsc
     if(iscf_mod>0 .and. inonsc==nnsclo_now)then

!     The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997))
      weight=two*occ_k(iband)*wtk_k/gs_hamkq%ucvol

!     One needs the Fourier transform of cwave1. However, only the one of
!     cwavef is available. If cwavef and cwave1 differs, this Fourier
!     transform must be computed. In both case the result is in wfraug1.
      if(tocceig==1)then
       tim_fourwf=5
       if(ispinor==1)then
        call fourwf(cplex,rhoaug1,cwave1,dummy,wfraug1,&
&        gs_hamkq%gbound,gs_hamkq%gbound,&
&        istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&        npw1_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
       else
        allocate(cwavef_sp(2,npw1_k))
!       $OMP PARALLEL DO PRIVATE(ipw) &
!       $OMP&SHARED(cwave1,cwavef_sp,npw1_k)
        do ipw=1,npw1_k
         cwavef_sp(1,ipw)=cwave1(1,ipw+npw1_k)
         cwavef_sp(2,ipw)=cwave1(2,ipw+npw1_k)
        end do
!       $OMP END PARALLEL DO
        call fourwf(cplex,rhoaug1,cwavef_sp,dummy,wfraug1,&
&        gs_hamkq%gbound,gs_hamkq%gbound,&
&        istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&        npw1_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
        deallocate(cwavef_sp)
       end if
      end if

      tim_fourwf=5
      if(ispinor==1)then
       call fourwf(1,rhoaug,cwave0,dummy,wfraug,gbound,gbound,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       npw_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
      else
       allocate(cwavef_sp(2,npw_k))
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(cwave0,cwavef_sp,npw_k)
       do ipw=1,npw_k
        cwavef_sp(1,ipw)=cwave0(1,ipw+npw_k)
        cwavef_sp(2,ipw)=cwave0(2,ipw+npw_k)
       end do
!      $OMP END PARALLEL DO
       call fourwf(1,rhoaug,cwavef_sp,dummy,wfraug,gbound,gbound,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&       npw_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
       deallocate(cwavef_sp)
      end if

!     Accumulate density
      if(cplex==2)then
!      $OMP PARALLEL DO PRIVATE(im0,im1,i1,i2,i3,re0,re1) &
!      $OMP&SHARED(n1,n2,n3,rhoaug1,weight,wfraug,wfraug1)
       do i3=1,n3
        do i2=1,n2
         do i1=1,n1
          re0=wfraug(1,i1,i2,i3)  ; im0=wfraug(2,i1,i2,i3)
          re1=wfraug1(1,i1,i2,i3) ; im1=wfraug1(2,i1,i2,i3)
          rhoaug1(2*i1-1,i2,i3)=rhoaug1(2*i1-1,i2,i3)+weight*(re0*re1+im0*im1)
          rhoaug1(2*i1  ,i2,i3)=rhoaug1(2*i1  ,i2,i3)+weight*(re0*im1-im0*re1)
         end do
        end do
       end do
!      $OMP END PARALLEL DO
      else
!      $OMP PARALLEL DO PRIVATE(i1,i2,i3) &
!      $OMP&SHARED(n1,n2,n3,rhoaug1,weight,wfraug,wfraug1)
       do i3=1,n3
        do i2=1,n2
         do i1=1,n1
          rhoaug1(i1,i2,i3)=rhoaug1(i1,i2,i3)+&
&          weight*( wfraug(1,i1,i2,i3)*wfraug1(1,i1,i2,i3) &
&          +wfraug(2,i1,i2,i3)*wfraug1(2,i1,i2,i3)  )
         end do
        end do
       end do
!      $OMP END PARALLEL DO
      end if

!     End of SCF case
     end if

    end do ! ispinor=1,nspinor

!   Local potential energy of this band, valuer has been accumulated
    eloc0_k(iband)=two*valuer/dble(gs_hamkq%nfft)

!   PAW: accumulate contribution to 1st-order occupancies matrix (rhoij1)
!   only at the last inonsc
    if(gs_hamkq%usepaw==1.and.iscf_mod>0.and.inonsc==nnsclo_now)then
     if (usecprj==1) then
      call pawaccrhoij(gs_hamkq%atindx1,cplex,cwaveprj0,cwaveprj1,dimpaw1,ipert,isppol,&
&      natom,nspden,nspinor,nsppol,occ_k(iband),option_rhoij,pawrhoij1,wtk_k)
     else
      allocate(cwaveprj_tmp(dimpaw1,nspinor))
      call cprj_alloc(cwaveprj_tmp,cprj (1,1)%ncpgr,dimcprj)
      if (ipert<=natom) then
       nattyp_atm(1)=1;allocate(ph1d_atm(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
       i1=(gs_hamkq%atindx(ipert)-1)*(2*n1+1)
       ph1d_atm(:,1:2*n1+1)=gs_hamkq%ph1d(:,1+i1:2*n1+1+i1)
       i2=(gs_hamkq%atindx(ipert)-1)*(2*n2+1)+natom*(2*n1+1)
       ph1d_atm(:,1+2*n1+1:2*n2+1+2*n1+1)=gs_hamkq%ph1d(:,1+i2:2*n2+1+i2)
       i3=(gs_hamkq%atindx(ipert)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
       ph1d_atm(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=gs_hamkq%ph1d(:,1+i3:2*n3+1+i3)
       choice=2
       call getcprj(choice,0,cwave0,cwaveprj_tmp,dimekb,1,dimffnlk,&
&       ekb_typ,ffnlk,idir,indlmn_typ,istwf_k,kg_k,kpg_k,kpt,lmnmax,&
&       matblk,mgfft,mpi_enreg,1,nattyp_atm,gs_hamkq%ngfft,nkpg,gs_hamkq%nloalg,&
&       npw_k,nspinor,1,gs_hamkq%phkxred,ph1d_atm,ph3d,gs_hamkq%ucvol,gs_hamkq%usepaw,gs_hamkq%useylm)
      end if
      call pawaccrhoij(gs_hamkq%atindx1,cplex,cwaveprj_tmp,cwaveprj1,dimpaw1,ipert,isppol,&
&      natom,nspden,nspinor,nsppol,occ_k(iband),option_rhoij,pawrhoij1,wtk_k)
      if (ipert<=natom) deallocate(ph1d_atm)
      call cprj_free(cwaveprj_tmp)
      deallocate(cwaveprj_tmp)
     end if
    end if

!   End of non-zero occupation
   end if

!  Exit loop over inonsc if converged and if non-self-consistent
   if (iscf_mod<0 .and. resid<dtset%tolwfr) exit

!  End loop over inonsc
  end do

! Write first-order eigenvalues and wavefunctions
  if(mk1mem/=0)then
   cg1(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)=cwave1(:,:)
  else
   call timab(288,1,tsec)
   call WffWriteDataRec(eig1_k,ierr,2*nband_k,wffnew)
   call WffWriteDataRec(cwave1,ierr,2*npw1_k*nspinor,wffnew)
   call timab(288,2,tsec)
  end if

! PAW: write first-order projected wavefunctions
  if (psps%usepaw==1.and.usecprj==1) then
   call cprj_put(.true.,gs_hamkq%atindx,cwaveprj,cprj1,dimpaw1,iband,ibg1,ikpt,iorder_cprj,isppol,&
&   mband,mk1mem,mpi_enreg,natom,1,nband_k,dimcprj,nspinor,nsppol,0,dtfil%unpaw1)
  end if

  if(prtvol>=10)then
   call status(counter,dtfil%filstat,iexit,level,'get residk    ')
  end if

! End loop over bands
 end do

!Find largest resid over bands at this k point
 residk=maxval(resid_k(:))
 if(prtvol>2 .or. ikpt<=nkpt_max)then
  do ii=0,(nband_k-1)/8
   write(message, '(1p,8e10.2)' ) &
&   (resid_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
   call wrtout(06,message,'PERS')
  end do
 end if

 call timab(139,2,tsec)
 call timab(130,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'after loops   ')

 deallocate(rhoaug,wfraug,wfraug1)
 deallocate(cwave0,cwavef,cwave1)
 deallocate(gh1_n)
 if (gs_hamkq%usepaw==1) then
  deallocate(dcwavef,gscq)
  call cprj_free(cwaveprj)
  deallocate(cwaveprj)
  if (usecprj==1) then
   call cprj_free(cwaveprj0)
   deallocate(cwaveprj0)
  end if
 end if

!###################################################################

!DEBUG
!write(6,*)'vtowfk3: iscf_mod, nband_k',iscf_mod, nband_k
!ENDDEBUG

 deallocate(ghc,gvnlc,gvnl1)

!Write the number of one-way 3D ffts skipped until now (in case of fixed
!occupation numbers
 if(iscf_mod>0 .and. (prtvol>2 .or. ikpt<=nkpt_max))then
  write(message, '(a,i8)' )&
&  ' vtowfk3 : number of one-way 3D ffts skipped in vtowfk until now =',nskip
  call wrtout(06,message,'PERS')
 end if

 if(prtvol<=2 .and. ikpt==nkpt_max+1)then
  write(message, '(a,a,a)' ) ch10,&
&  ' vtowfk3 : prtvol=0, 1 or 2, do not print more k-points.',ch10
  call wrtout(06,message,'PERS')
 end if

!###################################################################

 if (residk>dtset%tolwfr .and. iscf_mod<=0 .and. iscf_mod/=-3) then
  write(message, '(a,a,a,a,2i5,a,es13.5)' ) ch10,&
&  ' vtowfk3: WARNING -',ch10,&
&  '  Wavefunctions not converged for nnsclo,ikpt=',nnsclo_now,ikpt,&
&  ' max resid=',residk
  call wrtout(06,message,'PERS')
 end if

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i2,a)') ch10,&
&  ' vtowfk3 : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(130,2,tsec)
 call timab(128,2,tsec)

!DEBUG
!write(6,*)' vtowfk3 : exit '
!call flush(6)
!if(count==26)stop
!stop
!ENDDEBUG

end subroutine vtowfk3
!!***
