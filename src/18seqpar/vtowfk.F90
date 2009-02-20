!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtowfk
!! NAME
!! vtowfk
!!
!! FUNCTION
!! This routine compute the partial density at a given k-point,
!! for a given spin-polarization, from a fixed Hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cgq = array that holds the WF of the nearest neighbours of
!!        the current k-point (electric field, MPI //)
!!  cpus= cpu time limit in seconds
!!  dimcprj(natom*usepaw)=array of dimensions of array cprj
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cprj
!!  icg=shift to be applied on the location of data in the array cg
!!  ikg=shift to be given to the location of the data in the arrays kg
!!      and pwind
!!  ikpt=number of the k-point
!!  iscf=(<= 0  =>non-SCF), >0 => SCF
!!  isppol isppol=1 for unpolarized, 2 for spin-polarized
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (Hartree)
!!  kpg_k(npw,nkpg)= (k+G) components (only if useylm=1)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mcg=second dimension of the cg array
!!  mcgq = second dimension of the cgq array (electric field, MPI //)
!!  mgfft=maximum size of 1D FFTs
!!  mkgq = second dimension of pwnsfacq
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
!!  nkpt=number of k points.
!!  nnsclo_now=number of non-self-consistent loops for the current vtrial
!!             (often 1 for SCF calculation, =nstep for non-SCF calculations)
!!  npw_k=number of plane waves at this k point
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  nvloc=final dimension of vlocal and rhoaug (usually 1, but 4 for non-collinear)
!!  n4,n5,n6=integers used for dimensionning of vlocal
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  optforces=option for the computation of forces
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pwind(pwind_alloc,2,3)= array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc= first dimension of pwind
!!  pwnsfac(2,pwind_alloc)= phase factors for non-symmorphic translations
!!                          (see initberry.f)
!!  pwnsfacq(2,mkgq)= phase factors for the nearest neighbours of the
!!                    current k-point (electric field, MPI //)
!!  usebanfft=flag for band-fft parallelism
!!  usecprj= 1 if cprj array is stored in memory
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!  wtk=weight assigned to the k point.
!!  zshift(nband_k)=energy shifts for the squared shifted hamiltonian algorithm
!!
!! OUTPUT
!!  dphase_k(3)=change in Zak phase for the current k-point
!!  eig_k(nband_k)=array for holding eigenvalues (hartree)
!!  ek_k(nband_k)=contribution from each band to kinetic energy, at this k-point
!!  resid_k(nband_k)=residuals for each band over all k points,
!!                   BEFORE the band rotation.
!!  ==== if optforces>0 ====
!!    grnl_k(3*natom,nband_k)=nonlocal gradients, at this k-point
!!  ==== if (gs_hamk%usepaw==0) ====
!!    enl_k(nband_k)=contribution from each band to nonlocal pseudopotential part
!!                   of total energy, at this k-point
!!  ==== if (gs_hamk%usepaw==1) ====
!!    cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions projected with non-local projectors:
!!                                 cprj(n,k,i)=<p_i|Cnk> where p_i is a non-local projector.
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=updated wavefunctions
!!  rhoaug(n4,n5,n6,nvloc)= density in electrons/bohr**3, on the augmented fft grid.
!!                    (cumulative, so input as well as output). Update only
!!                    for occopt<3 (fixed occupation numbers)
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      cgwf,cprj_alloc,cprj_free,cprj_put,dsymm,fourwf,fxphas,leave_new
!!      lobpcgcciiwf,lobpcgccwf,lobpcgiiwf,lobpcgwf,meanvalue_g,nonlop
!!      prep_fourwf,prep_nonlop,pw_orthon,status,subdiago,timab,wrtout
!!      xcomm_init,xsum_mpi,zhemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine vtowfk(cg,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,dtset,&
& eig_k,ek_k,enl_k,ffnl,ffnl_gather,grnl_k,gs_hamk,&
& ibg,icg,ikg,ikpt,iscf,isppol,kg_k,kg_k_gather,kinpw,kinpw_gather,kpg_k,&
& lmnmax,matblk,mcg,mcgq,mgfft,mkgq,mkmem,mpi_enreg,mpsang,&
& mpssoang,mpw,natom,nband_k,ndatarecv,nkpg,nkpt,nnsclo_now,npw_k,npwarr,&
& nspinor,ntypat,nvloc,n4,n5,n6,occ_k,optforces,ph3d,ph3d_gather,prtvol,psps,&
& pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,rhoaug,usebandfft,usecprj,vlocal,wtk,zshift)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12spacepar
 use interfaces_13nonlocal
 use interfaces_14wfs
 use interfaces_18seqpar, except_this_one => vtowfk
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 type(gs_hamiltonian_type), intent(in) :: gs_hamk
 integer, intent(in) :: dimffnl,ibg,icg,ikg,ikpt,iscf,isppol,lmnmax,matblk
 integer, intent(in) :: mcg,mcgq,mgfft,mkgq,mkmem,mpsang,mpssoang,mpw,n4,n5,n6
 integer, intent(in) :: natom,nband_k,ndatarecv,nkpg,nkpt,nnsclo_now,npw_k,nspinor,ntypat,nvloc,optforces
 integer, intent(in) :: prtvol,pwind_alloc,usebandfft,usecprj
 real(dp), intent(in) :: cpus,wtk
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in) :: dtset
 type(pseudopotential_type), intent(in) :: psps
 type(MPI_type), intent(inout) :: mpi_enreg
 type(efield_type), intent(inout) :: dtefield
 integer, intent(in) :: dimcprj(natom*gs_hamk%usepaw),kg_k(3,npw_k),kg_k_gather(3,ndatarecv*usebandfft)
 integer, intent(in) :: npwarr(nkpt),pwind(pwind_alloc,2,3)
 real(dp), intent(in) :: cgq(2,mcgq),ffnl(npw_k,dimffnl,lmnmax,ntypat)
 real(dp), intent(in) :: ffnl_gather(ndatarecv,dimffnl,lmnmax,ntypat*usebandfft)
 real(dp), intent(in) :: kinpw(npw_k),kinpw_gather(ndatarecv*usebandfft),kpg_k(npw_k,nkpg),occ_k(nband_k)
 real(dp), intent(inout) :: ph3d(2,npw_k,matblk),ph3d_gather(2,ndatarecv,matblk*usebandfft)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc),pwnsfacq(2,mkgq)
 real(dp), intent(in) :: zshift(nband_k)
 real(dp), intent(inout) :: vlocal(n4,n5,n6,nvloc)
 real(dp), intent(out) :: eig_k(nband_k),ek_k(nband_k),dphase_k(3)
 real(dp), intent(out) :: enl_k(nband_k*(1-gs_hamk%usepaw))
 real(dp), intent(out) :: grnl_k(3*natom,nband_k*optforces)
 real(dp), intent(out) :: resid_k(nband_k)
 real(dp), intent(inout) :: cg(2,mcg),rhoaug(n4,n5,n6,nvloc)
 type(cprj_type) :: cprj(natom,nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*gs_hamk%usepaw)

!Local variables-------------------------------
 integer :: spacecomm_old=0
 integer,parameter :: level=8
 integer,save :: nskip=0
!     Flag use_subovl: 1 if "subovl" array is computed (see below)
!     subovl should be Identity (in that case we should use use_subovl=0)
!     But this is true only if conjugate gradient algo. converges
 integer :: use_subovl=0
 integer :: blocksize,choice,counter,cpopt,i1,i2,i3,ia,iatom,iband,iband2
 integer :: iblock,iblocksize,ibs,ider,idir,ier,ierr,iexit,ifft,ig,igs,igsc,ii,ilmn,inonsc
 integer :: iorder_cprj,ipw,ispinor,istwf_k,isubh,isubo,itypat,iwavef,j1,j2,jj,jlmn,kk,mgsc,n1,n2,n3
 integer :: nblockbd,nnlout,old_paral_level,ortalgo
 integer :: paw_opt,quit,signs,spaceComm,tim_fourwf,tim_nonlop,wfoptalg,wfopta10
 integer :: me_distrb
 real(dp) :: ar,arg,eshift,lambda_k,occblock
 real(dp) :: residk,sumi,sumr,weight
 character(len=500) :: message
 real(dp) :: dummy(2,1),nonlop_dum(1,1),tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef1(:,:),cwavef2(:,:),eig_save(:),enlout(:),evec(:,:)
 real(dp),allocatable :: cwavef_x(:,:),cwavef_y(:,:),evec_loc(:,:),gsc(:,:)
 real(dp),allocatable :: lambda_loc(:),mat_loc(:,:),mat1(:,:,:),matvnl(:,:,:)
 real(dp),allocatable :: subham(:),subovl(:),subvnl(:),totvnl(:,:),wfraug(:,:,:,:)
 complex(dp), external :: zdotc
 type(cprj_type),allocatable :: cwaveprj(:,:)

!no_abirules
#if defined MPI
            integer :: mpi_status(MPI_STATUS_SIZE)
            integer :: nskip_rs
#endif
#if !defined MPI
            integer,parameter :: nkpt_max=50
#endif
#if defined MPI
            integer,parameter :: nkpt_max=-1
#endif
#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZDOTC' :: zdotc
!DEC$ ATTRIBUTES ALIAS:'ZHEMM' :: zhemm
#endif

! *********************************************************************

!PATCH vtowfk // KPT & FFT me-->me_kpt
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
    &(mpi_enreg%paral_compil_fft==1)) then
    me_distrb = mpi_enreg%me_kpt
 else
    me_distrb = mpi_enreg%me
 end if

!=========================================================================
!============= INITIALIZATIONS AND ALLOCATIONS ===========================
!=========================================================================

!Keep track of total time spent in vtowfk
 call timab(28,1,tsec)
!DEBUG
! write(6,*)' vtowfk : enter '
! write(6,*)' vlocal(1,1,1,1)=',vlocal(1,1,1,1)
!  write(6,*)' nband_k,natom,npw_k',nband_k,natom,npw_k
!stop
!ENDDEBUG
 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,'vtowfk : enter'
  call wrtout(06,message,'PERS')
 end if

 wfoptalg=dtset%wfoptalg;wfopta10=mod(wfoptalg,10)
 istwf_k=gs_hamk%istwf_k
 quit=0

!Define several block values for band parallelization
 if (wfopta10/=4.and.wfopta10/=5) then
  nblockbd=nband_k
  blocksize=1
 elseif (wfopta10==5) then
  nblockbd=1
  blocksize=nband_k
 else
  nblockbd=nband_k/mpi_enreg%nproc_fft
  if (nband_k/=nblockbd*mpi_enreg%nproc_fft) nblockbd=nblockbd+1
  if(mpi_enreg%mode_para=='b') then
   nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
   write(6,*) 'starting lobpcg, with nblockbd,mpi_enreg%nproc_band',nblockbd,mpi_enreg%nproc_band
  end if
  blocksize=nband_k/nblockbd
 end if

!Save eshift
 if(wfoptalg==3)then
  eshift=zshift(1)
  allocate(eig_save(nband_k));eig_save(:)=eshift
 end if

 n1=gs_hamk%ngfft(1) ; n2=gs_hamk%ngfft(2) ; n3=gs_hamk%ngfft(3)

 igsc=0
 mgsc=nband_k*npw_k*nspinor*gs_hamk%usepaw
 allocate(gsc(2,mgsc))

 allocate(evec(2*nband_k,nband_k))
 allocate(subham(nband_k*(nband_k+1)))
 if (gs_hamk%usepaw==0) then
  if (wfopta10==4) then
   if (istwf_k==1) then
    allocate(totvnl(2*nband_k,nband_k))
   else if (istwf_k==2) then
    allocate(totvnl(nband_k,nband_k))
   end if
  else
   allocate(subvnl(nband_k*(nband_k+1)))
  endif
  use_subovl=0
 end if
 if (use_subovl==1) then
  allocate(subovl(nband_k*(nband_k+1)))
 end if

!Carry out UP TO dtset%nline steps,
!or until resid for every band is < dtset%tolwfr

 if(prtvol>2 .or. ikpt<=nkpt_max)then
  write(message, '(a,a,i5,2x,a,3f9.5,2x,a)' ) ch10,&
&  ' Non-SCF iterations; k pt #',ikpt,'k=',gs_hamk%kpoint(:),'band residuals:'
  call wrtout(06,message,'PERS')
 end if

!Electric field: initialize dphase_k
 dphase_k(:) = zero

 call status(0,dtfil%filstat,iexit,level,'before loop   ')
 call timab(39,1,tsec)

!=========================================================================
!==================== NON-SELF-CONSISTENT LOOP ===========================
!=========================================================================

!nnsclo_now=number of non-self-consistent loops for the current vtrial
!    (often 1 for SCF calculation, =nstep for non-SCF calculations)
 do inonsc=1,nnsclo_now

! DEBUG
!  write(6,*)' vtowfk : gs_hamk%usepaw=',gs_hamk%usepaw
!  write(6,*)' vtowfk : cg(1:2) for different bands '
!  do iband=1,nband_k
!   iwavef=(iband-1)*npw_k+icg
!  write(6,*)cg(1:2,1+iwavef)
!  end do
! ENDDEBUG

! This initialisation is needed for the MPI-parallelisation (gathering using sum)
  subham(:)=zero
  if (gs_hamk%usepaw==0) then
   if (wfopta10==4) then
    totvnl(:,:)=zero
   else
    subvnl(:)=zero
   endif
  endif
  if (use_subovl==1)subovl(:)=zero
  resid_k(:)=0._dp

!Filter the WFs when modified kinetic energy is too large (see routine mkkin.f)
  do ispinor=1,nspinor
   igs=(ispinor-1)*npw_k
   do iband=1,nband_k
    iwavef=(iband-1)*npw_k*nspinor+icg
    do ipw=1+igs,npw_k+igs
     if(kinpw(ipw-igs)>huge(zero)*1.d-11)then
      cg(1,ipw+iwavef)=zero
      cg(2,ipw+iwavef)=zero
     end if
    end do
   end do
  end do

!=========================================================================
!============ MINIMIZATION OF BANDS: LOBPCG ==============================
!=========================================================================

!DEBUG
! write(6,*)' vtowfk : wfoptalg,istwf_k', wfoptalg,istwf_k
!ENDDEBUG

  if(wfopta10==4.or.wfopta10==5) then
   if (istwf_k/=1.and.istwf_k/=2) then !no way to use lobpcg
    write(message, '(6a)' ) ch10,&
&    ' vtowfk : ERROR -',ch10,'  Only istwfk=1 or 2 are allowed with wfoptalg=4/14 !',ch10,&
&    '  Action: put istwfk to 1 or remove k points with half integer coordinates.'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   if (wfopta10==4) then ! This is lobpcg I
    if(istwf_k==2) then ! we use the symmetric version
     call lobpcgwf(cg,dimffnl,dtfil,dtset,ffnl,ffnl_gather,gs_hamk,gsc,icg,igsc,&
&            kg_k,kg_k_gather,kinpw,kinpw_gather,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&
&            nband_k,nblockbd,ndatarecv,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,ph3d_gather,prtvol,&
&            psps,resid_k,subham,subovl,totvnl,usebandfft,use_subovl,vlocal)
    else if(istwf_k==1) then !we use the hermitian version
     call lobpcgccwf(cg,dimffnl,dtfil,dtset,ffnl,ffnl_gather,gs_hamk,gsc,icg,igsc,&
&            kg_k,kg_k_gather,kinpw,kinpw_gather,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&
&            nband_k,nblockbd,ndatarecv,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,ph3d_gather,prtvol,&
&            psps,resid_k,subham,subovl,totvnl,usebandfft,use_subovl,vlocal)
    end if
   else if (wfopta10==5) then ! This is lobpcg II
    if(istwf_k==2) then ! we use the symmetric version
     call lobpcgIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
&            kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&
&            nband_k,nblockbd,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&
&            psps,resid_k,subham,subovl,subvnl,use_subovl,vlocal)
    else if(istwf_k==1) then !we use the hermitian version
     call lobpcgccIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
&            kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&
&            nband_k,nblockbd,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&
&            psps,resid_k,subham,subovl,subvnl,use_subovl,vlocal)
    end if
   end if

!DEBUG
!  write(6,*)' vtowfk : after lobpcgXX routines'
!ENDDEBUG

!  In case of FFT parallelism, exchange subspace arrays
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
   call xsum_mpi(subham,spaceComm,ierr)

   if (gs_hamk%usepaw==0) then
    if (wfopta10==4) then
     call xsum_mpi(totvnl,spaceComm,ierr)
    else
     call xsum_mpi(subvnl,spaceComm,ierr)
    endif
   endif

   if (use_subovl==1) call xsum_mpi(subovl,spaceComm,ierr)
   mpi_enreg%paral_level= old_paral_level

!=========================================================================
!======== MINIMIZATION OF BANDS: CONJUGATE GRADIENT (Teter et al.) =======
!=========================================================================
  else
   call cgwf(dtset%berryopt,cg,cgq,dtset%chkexit,cpus,dimffnl,dphase_k,dtefield,&
&             ffnl,dtfil%filnam_ds(1),dtfil%filstat,&
&             gsc,gs_hamk,icg,igsc,ikg,ikpt,inonsc,&
&             isppol,kg_k,kinpw,lmnmax,matblk,dtset%mband,&
&             mcg,mcgq,mgfft,mgsc,mkgq,mkmem,mpi_enreg,mpsang,&
&             mpssoang,mpw,natom,nband_k,dtset%nbdblock,nkpt,dtset%nline,npw_k,npwarr,nspinor,&
&             dtset%nsppol,ntypat,nvloc,n4,n5,n6,dtset%ortalg,&
&             dtset%paral_kgb,ph3d,prtvol,pwind,pwind_alloc,pwnsfac,&
&             pwnsfacq,quit,resid_k,subham,subovl,subvnl,dtset%tolwfr,&
&             use_subovl,vlocal,wfoptalg,wtk,zshift)
   end if

!=========================================================================
!===================== FIND LARGEST RESIDUAL =============================
!=========================================================================

! Find largest resid over bands at this k point
! Note that this operation is done BEFORE rotation of bands :
! it would be time-consuming to recompute the residuals after.
  residk=maxval(resid_k(1:max(1,nband_k-dtset%nbdbuf)))
! Print residuals
  if(prtvol>2 .or. ikpt<=nkpt_max)then
   do ii=0,(nband_k-1)/8
    write(message, '(a,8es10.2)' ) &
&       ' res:',(resid_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
    call wrtout(06,message,'PERS')
   end do
  end if

!=========================================================================
!========== DIAGONALIZATION OF HAMILTONIAN IN WFs SUBSPACE ===============
!=========================================================================

  call timab(585,1,tsec)
  call subdiago(cg,dtfil%filstat,eig_k,evec,gsc,icg,igsc,ikpt,inonsc,istwf_k,&
&               mcg,mgsc,mpi_enreg,nband_k,npw_k,nspinor,dtset%paral_kgb,prtvol,&
&               subham,subovl,use_subovl,gs_hamk%usepaw)
  call timab(585,2,tsec)

! Print energies
  if(prtvol>2 .or. ikpt<=nkpt_max)then
#if defined MPI
   if ((mpi_enreg%paralbd <=1) .or. ((mpi_enreg%paralbd >1) .and. &
&      (mpi_enreg%me_group==0))) then
#endif
    do ii=0,(nband_k-1)/8
     write(message, '(a,8es10.2)' ) &
&        ' ene:',(eig_k(iband),iband=1+ii*8,min(nband_k,8+ii*8))
     call wrtout(06,message,'PERS')
    end do
#if defined MPI
   end if
#endif
  end if


! THIS CHANGE OF SHIFT DOES NOT WORK WELL
! Update zshift in the case of wfoptalg==3
! if(wfoptalg==3 .and. inonsc/=1)then
!  do iband=1,nband_k
!   if(eig_k(iband)<eshift .and. eig_save(iband)<eshift)then
!    zshift(iband)=max(eig_k(iband),eig_save(iband))
!   end if
!   if(eig_k(iband)>eshift .and. eig_save(iband)>eshift)then
!    zshift(iband)=min(eig_k(iband),eig_save(iband))
!   end if
!  end do
!  eig_save(:)=eig_k(:)
! end if

!=========================================================================
!=============== ORTHOGONALIZATION OF WFs (if needed) ====================
!=========================================================================

! Re-orthonormalize the wavefunctions at this k point--
! this is redundant but is performed to combat rounding
! error in wavefunction orthogonality
  call status(inonsc,dtfil%filstat,iexit,level,'call pw_orthon   ')
  if (mpi_enreg%paral_fft == 1) then
   mpi_enreg%num_group_fft=ikpt + (isppol-1)*nkpt
  else
   mpi_enreg%num_group_fft=0
  end if

  call timab(583,1,tsec)
  ortalgo=dtset%paral_kgb
  call pw_orthon(icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,npw_k*nspinor,nband_k,ortalgo,gsc,gs_hamk%usepaw,cg)
  call timab(583,2,tsec)
!  DEBUG
!   write(6,*)' vtowfk : cg(1:2) for different bands (4) '
!   do iband=1,nband_k
!    iwavef=(iband-1)*npw_k+icg
!    write(6,*)cg(1:2,1+iwavef)
!   end do
!  ENDDEBUG

! Fix phases of all bands
  call status(inonsc,dtfil%filstat,iexit,level,'call fxphas   ')

!PATCH vtowfk !! fxphas
!DEBUG seq==par comment next block
  if ((mpi_enreg%paral_compil_kpt/=1) .or. &
     &(mpi_enreg%paral_compil_fft/=1)) then
     call fxphas(cg,gsc,icg,igsc,istwf_k,mcg,mgsc,mpi_enreg,nband_k,npw_k*nspinor,gs_hamk%usepaw)
  end if

!  DEBUG
!   write(6,*)' vtowfk : cg(1:2) for different bands (5) '
!   do iband=1,nband_k
!    iwavef=(iband-1)*npw_k+icg
!    write(6,*)cg(1:2,1+iwavef)
!   end do
!  ENDDEBUG

!=========================================================================
!================= END OF NON SELF-CONSISTENT LOOP =======================
!=========================================================================

! Exit loop over inonsc if converged
  if (residk<dtset%tolwfr) exit

!End loop over inonsc
 end do

 call timab(39,2,tsec)
 call timab(30,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'after loops   ')

!###################################################################

!Compute kinetic energy and non-local energy for each band, and in the SCF
!case, contribution to forces, and eventually accumulate rhoaug.

!DEBUG
!write(6,*)'vtowfk : iscf, nband_k, occopt',iscf, nband_k, dtset%occopt
!ENDDEBUG

 if(iscf>0 .and. dtset%occopt<3) allocate(wfraug(2,n4,n5,n6))

!This loop can also be MPI-parallelized : the partial contributions
!to ek_k, grnl_k and rhoaug will have to be gathered (a simple summation
!of all contributions) after the loop
 if (wfopta10==4 .or. wfopta10==5) then
  allocate(cwavef(2,npw_k*nspinor))
  do iband=1,nband_k
   cwavef(:,:)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
!  Compute kinetic energy of each band
   call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
   ek_k(iband)=ar
  end do
  deallocate(cwavef)
 end if

!"nonlop" routine input parameters
 nnlout=3*natom*optforces
 signs=1;idir=0
 if (gs_hamk%usepaw==0) then
  choice=1+optforces
  paw_opt=0;cpopt=-1;tim_nonlop=2
 else
  choice=2*optforces
  paw_opt=2;cpopt=0;tim_nonlop=10-8*optforces
 end if
 allocate(enlout(nnlout*blocksize))

 if (wfopta10==5) then
  nblockbd=nband_k; blocksize=1
 endif

!Allocation of memory space for one WF
 allocate(cwavef(2,npw_k*nspinor*blocksize))
 if (gs_hamk%usepaw==1.and.iscf>0) then
  iorder_cprj=0
  allocate(cwaveprj(natom,nspinor*blocksize))
  call cprj_alloc(cwaveprj,0,dimcprj)
 end if

!Loop over bands or blocks of bands
 do iblock=1,nblockbd
  occblock=maxval(occ_k(1+(iblock-1)*blocksize:iblock*blocksize))

  if(wfopta10==4)then
   cwavef(:,1:npw_k*nspinor*blocksize)=cg(:,1+(iblock-1)*npw_k*nspinor*blocksize+icg:iblock*npw_k*nspinor*blocksize+icg)
  else if (wfopta10==5) then
   cwavef(:,1:npw_k*nspinor)=cg(:,1+(iblock-1)*npw_k*nspinor+icg:iblock*npw_k*nspinor+icg)
  else
#if defined MPI
   if(wfopta10==1)then
            if(mpi_enreg%parareel == 0) then
             if(mpi_enreg%proc_distrb(ikpt,iblock,isppol)/= me_distrb)&
&               then
                ek_k(iblock)=zero
                if (optforces>0) grnl_k(:,iblock)=zero
                cycle
             end if
            else
             if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt)/= mpi_enreg%me)&
&               then
                ek_k(iblock)=zero
                if (optforces>0) grnl_k(:,iblock)=zero
                cycle
             end if
            end if
   end if
#endif
   cwavef(:,:)=cg(:,1+(iblock-1)*npw_k*nspinor+icg:iblock*npw_k*nspinor+icg)
!  Compute kinetic energy of each band
   call meanvalue_g(ar,kinpw,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
   ek_k(iblock)=ar
  end if


!DEBUG
! if(iblock==1 .and. ikpt==1)then
!  write(6,*)' vtowfk : ig,kinpw(ig);ek_k(iblock)= ',ek_k(iblock)
!  do ig=1,npw_k
!   write(6,*)ig,kinpw(ig)
!  end do
! end if
!ENDDEBUG

  if(iscf>0)then

!  In case of fixed occupation numbers, accumulates the partial density
   if (dtset%occopt<3 .and. mpi_enreg%mode_para/='b') then
    if (abs(occ_k(iblock))>=tol8) then
     weight=occ_k(iblock)*wtk/gs_hamk%ucvol
!    Accumulate charge density in real space in array rhoaug
     tim_fourwf=2

!    The same section of code is also found in mkrho.F90 : should be rationalized !
     call fourwf(1,rhoaug(:,:,:,1),cwavef,dummy,wfraug,&
&     gs_hamk%gbound,gs_hamk%gbound,&
&     istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)

     if(nspinor==2)then
      allocate(cwavef1(2,npw_k))
      cwavef1(:,:)=cwavef(:,1+npw_k:2*npw_k)

      if(dtset%nspden==1) then

       call fourwf(1,rhoaug(:,:,:,1),cwavef1,dummy,wfraug,&
&       gs_hamk%gbound,gs_hamk%gbound,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)

      else if(dtset%nspden==4) then

!      Build the four components of rho. We use only norm quantities and, so fourwf.
!      $\sum_{n} f_n \Psi^{* \alpha}_n \Psi^{\alpha}_n =\rho^{\alpha \alpha}$
!      $\sum_{n} f_n (\Psi^{1}+\Psi^{2})^*_n (\Psi^{1}+\Psi^{2})_n=rho+m_x$
!      $\sum_{n} f_n (\Psi^{1}-i \Psi^{2})^*_n (\Psi^{1}-i \Psi^{2})_n=rho+m_y$
       allocate(cwavef_x(2,npw_k),cwavef_y(2,npw_k))
!      $(\Psi^{1}+\Psi^{2})$
       cwavef_x(:,:)=cwavef(:,1:npw_k)+cwavef1(:,1:npw_k)
!      $(\Psi^{1}-i \Psi^{2})$
       cwavef_y(1,:)=cwavef(1,1:npw_k)+cwavef1(2,1:npw_k)
       cwavef_y(2,:)=cwavef(2,1:npw_k)-cwavef1(1,1:npw_k)
       call fourwf(1,rhoaug(:,:,:,4),cwavef1,dummy,wfraug,gs_hamk%gbound,gs_hamk%gbound,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)
       call fourwf(1,rhoaug(:,:,:,2),cwavef_x,dummy,wfraug,gs_hamk%gbound,gs_hamk%gbound,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)
       call fourwf(1,rhoaug(:,:,:,3),cwavef_y,dummy,wfraug,gs_hamk%gbound,gs_hamk%gbound,&
&       istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)

       deallocate(cwavef_x,cwavef_y)

      end if ! dtset%nspden/=4
      deallocate(cwavef1)
     end if
    else
     nskip=nskip+1
    end if
!DEBUG
!    write(6,*)' ikpt, iblock, rhoaug',ikpt,iblock,rhoaug(1,1,1)
!ENDDEBUG
   end if

!  Call to nonlocal operator:
!   - Compute nonlocal forces from most recent wfs
!   - PAW: compute contribution to augmentation occ. (rhoij)
   if (gs_hamk%usepaw==1.or.optforces/=0) then
!   Treat all wavefunctions in case of varying occupation numbers or PAW
!   Only treat occupied bands in case of fixed occupation numbers and NCPP
    if( (dtset%occopt<3.or.dtset%occopt>7).and.abs(occblock)<=tol8.and.gs_hamk%usepaw==0) then
     if (optforces>0) grnl_k(:,(iblock-1)*blocksize+1:iblock*blocksize)=zero
    else
     if(gs_hamk%usepaw==1) call timab(554,1,tsec)
     if (mpi_enreg%mode_para=='b') then
      allocate(lambda_loc(blocksize))
      tim_nonlop=360
      call timab(572,1,tsec)
      lambda_loc(1:blocksize)=eig_k(1+(iblock-1)*blocksize:iblock*blocksize)
      call prep_nonlop(gs_hamk%atindx1,choice,cpopt,cwaveprj,gs_hamk%dimekb1,&
&      gs_hamk%dimekb2,dimffnl,gs_hamk%ekb,enlout,&
&      ffnl_gather,gs_hamk%gmet,gs_hamk%gprimd,iblock,1,idir,gs_hamk%indlmn,istwf_k,kg_k_gather,&
&      kpg_k,gs_hamk%kpoint,lambda_loc,lmnmax,matblk,blocksize,mgfft,mpi_enreg,mpsang,mpssoang,&
&      natom,gs_hamk%nattyp,nband_k,nblockbd,ndatarecv,gs_hamk%ngfft,nkpg,gs_hamk%nloalg,nnlout,npw_k,&
&      nspinor,ntypat,dtset%paral_kgb,paw_opt,gs_hamk%phkxred,gs_hamk%ph1d,ph3d_gather,prtvol,gs_hamk%pspso,signs,&
&      gs_hamk%sij,nonlop_dum,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
      call timab(572,2,tsec)
      deallocate(lambda_loc)
     else
      lambda_k=eig_k(iblock)
      call nonlop(gs_hamk%atindx1,choice,cpopt,cwaveprj,&
&      gs_hamk%dimekb1,gs_hamk%dimekb2,dimffnl,dimffnl,gs_hamk%ekb,&
&      enlout,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,&
&      gs_hamk%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,gs_hamk%kpoint,gs_hamk%kpoint,&
&      lambda_k,lmnmax,matblk,mgfft,mpi_enreg,mpsang,&
&      mpssoang,natom,gs_hamk%nattyp,&
&      gs_hamk%ngfft,nkpg,nkpg,gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,&
&      gs_hamk%phkxred,gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,&
&      gs_hamk%pspso,signs,gs_hamk%sij,nonlop_dum,&
&      tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
     end if
     if(gs_hamk%usepaw==1) call timab(554,2,tsec)

     if (optforces>0) then
      iband=(iblock-1)*blocksize
      do iblocksize=1,blocksize
       iband=iband+1;ibs=nnlout*(iblocksize-1)
       grnl_k(1:nnlout,iband)=enlout(ibs+1:ibs+nnlout)
      end do
     end if
     if (gs_hamk%usepaw==1.and.usecprj==1) then
       iband=1+(iblock-1)*blocksize
       call cprj_put(.true.,gs_hamk%atindx,cwaveprj,cprj,natom,iband,ibg,ikpt,iorder_cprj,isppol,&
&        dtset%mband,dtset%mkmem,mpi_enreg,natom,blocksize,nband_k,dimcprj,nspinor,dtset%nsppol,0,dtfil%unpaw)
     end if
    end if

   end if

! End of SCF calculation
  end if

!End of loop on blocks
 end do

 deallocate(cwavef,enlout)
 if (gs_hamk%usepaw==1.and.iscf>0) then
  call cprj_free(cwaveprj)
  deallocate(cwaveprj)
 end if

!In case of fixed occupation numbers,in bandFFT mode accumulates the partial density
 if (dtset%occopt<3 .and. mpi_enreg%mode_para=='b' .and. wfopta10==4) then
  allocate(cwavef(2,npw_k*nspinor*blocksize))
  do iblock=1,nblockbd
   cwavef(:,1:npw_k*nspinor*blocksize)=&
&   cg(:,1+(iblock-1)*npw_k*nspinor*blocksize+icg:iblock*npw_k*nspinor*blocksize+icg)
   call timab(537,1,tsec)
   call prep_fourwf(rhoaug,blocksize,cwavef,wfraug,&
&   gs_hamk,istwf_k,iblock,1,kg_k_gather,mgfft,mpi_enreg,nblockbd,nband_k,ndatarecv,npw_k,n4,n5,n6,occ_k,dtset%paral_kgb,wtk)
   call timab(537,2,tsec)
  end do
 deallocate(cwavef)
 end if
!transpose the wave functions to have one band per band group
!call put_one_band_per_band_group(blocksize,cwavef,cwavef_alltoall,mpi_enreg,npw_k)
!cwavef=cwavef_alltoall ! each fft group of procs knows its wf.

#if defined MPI
    if ((wfopta10==1) .and. (mpi_enreg%paralbd >= 1)) then
           call timab(48,1,tsec)
             spaceComm = mpi_enreg%kpt_comm(mpi_enreg%num_group)
             call xsum_mpi(rhoaug,spaceComm,ierr)
!             call xsum_mpi(nskip,nskip_rs,spaceComm,ierr)
             call xsum_mpi(nskip,spaceComm,ierr)
                nskip_rs=nskip
           call timab(48,2,tsec)
    end if
#endif


 if (iscf>0 .and. dtset%occopt<3) deallocate(wfraug)

!Write the number of one-way 3D ffts skipped until now (in case of fixed
!occupation numbers
 if(iscf>0 .and. dtset%occopt<3 .and. (prtvol>2 .or. ikpt<=nkpt_max) )then
#if defined MPI
    if(wfopta10==1)then
  write(message, '(a,i8)' )&
&  ' vtowfk : number of one-way 3D ffts skipped in vtowfk until now =',nskip_rs
  call wrtout(06,message,'PERS')
       else
#endif
  write(message, '(a,i8)' )&
&  ' vtowfk : number of one-way 3D ffts skipped in vtowfk until now =',nskip
  call wrtout(06,message,'PERS')
 end if
#if defined MPI
       end if
#endif

!Norm-conserving only: Compute nonlocal part of total energy : rotate subvnl
 if (gs_hamk%usepaw==0) then

  call timab(586,1,tsec)
  allocate(matvnl(2,nband_k,nband_k),mat1(2,nband_k,nband_k))

  if (wfopta10==4) then
   mat1(:,:,:)=zero
   enl_k(1:nband_k)=zero
   if (istwf_k==1) then
    call zhemm('l','l',nband_k,nband_k,cone,totvnl,nband_k,&
&              evec,nband_k,czero,mat1,nband_k)
    do iband=1,nband_k
     enl_k(iband)=zdotc(nband_k,evec(:,iband),1,mat1(:,:,iband),1)
    enddo
   else if (istwf_k==2) then
    allocate(evec_loc(nband_k,nband_k),mat_loc(nband_k,nband_k))
    do iband=1,nband_k
     do jj=1,nband_k
      evec_loc(iband,jj)=evec(2*iband-1,jj)
     end do
    end do
    call dsymm('l','l',nband_k,nband_k,one,totvnl,nband_k,&
&              evec_loc,nband_k,zero,mat_loc,nband_k)
    do iband=1,nband_k
     enl_k(iband)=ddot(nband_k,evec_loc(:,iband),1,mat_loc(:,iband),1)
    enddo
    deallocate(evec_loc,mat_loc)
   end if
  else

!  (1) Write subvnl in full storage mode
   ii=0
   do iband=1,nband_k
    do jj=1,iband
     ii=ii+1
     matvnl(1,jj,iband)=subvnl(2*ii-1)
     matvnl(2,jj,iband)=subvnl(2*ii  )
    end do
   end do
   if (nband_k>1) then
    do iband=1,nband_k-1
     do jj=iband+1,nband_k
      matvnl(1,jj,iband)= matvnl(1,iband,jj)
      matvnl(2,jj,iband)=-matvnl(2,iband,jj)
     end do
    end do
   end if

!  Product of matvnl by evec : mat1(ii,jj)=Sum(kk) evec(kk,jj) matvnl(ii,kk)
!  Second product with evec, so that matvnl has been rotated
!  (new)matvnl(ii,jj)=Sum(kk) evec*(kk,ii) mat1(kk,jj)
!  However, only the diagonal, real, part is needed later
!  This loop can be MPI-parallelized, with transfer of enl_k(iband) to all
!  processors belonging to this k-point.
   do iband=1,nband_k
#if defined MPI
            if(wfopta10==1)then
             if(mpi_enreg%parareel == 0) then
              if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= me_distrb)&
&                then
                 enl_k(iband)=zero
                 cycle
              end if
             else
              if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt)/= mpi_enreg%me)&
&                then
                  enl_k(iband)=zero
                 cycle
              end if
             end if
            end if
#endif

    do ii=1,nband_k
     mat1(1,ii,iband)=zero
     mat1(2,ii,iband)=zero
     do kk=1,nband_k
      mat1(1,ii,iband)=mat1(1,ii,iband)+evec(2*kk-1,iband)*matvnl(1,ii,kk) &
&                                      -evec(2*kk  ,iband)*matvnl(2,ii,kk)
      mat1(2,ii,iband)=mat1(2,ii,iband)+evec(2*kk-1,iband)*matvnl(2,ii,kk) &
&                                      +evec(2*kk  ,iband)*matvnl(1,ii,kk)
     end do
    end do

    enl_k(iband)=zero
    do kk=1,nband_k
     enl_k(iband)=enl_k(iband)+evec(2*kk-1,iband)*mat1(1,kk,iband) &
&                             +evec(2*kk  ,iband)*mat1(2,kk,iband)
    end do
   end do
  endif

  deallocate(matvnl,mat1)
  call timab(586,2,tsec)
 end if

!###################################################################

 if (residk>dtset%tolwfr .and. iscf<=0) then
  write(message, '(a,a,a,a,2i5,a,es13.5)' ) ch10,&
&   ' vtowfk: WARNING -',ch10,&
&   '  Wavefunctions not converged for nnsclo,ikpt=',nnsclo_now,ikpt,&
&        ' max resid=',residk
  call wrtout(06,message,'PERS')
 end if

!Print out eigenvalues (hartree)
 if(prtvol>2 .or. ikpt<=nkpt_max) then
  write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&  'eigenvalues (hartree) for',nband_k,'bands',ch10,&
&  '              after ',inonsc,' non-SCF iterations with ',&
&  dtset%nline,' CG line minimizations'
  call wrtout(06,message,'PERS')
#if defined MPI
       if ((mpi_enreg%paralbd <=1) .or. (mpi_enreg%paralbd >1)) then
#endif
  do ii=0,(nband_k-1)/6
   write(message, '(1p,6e12.4)' ) (eig_k(iband),&
&                        iband=1+6*ii,min(6+6*ii,nband_k))
   call wrtout(06,message,'PERS')
  end do

#if defined MPI
       end if
#endif

 else if(ikpt==nkpt_max+1)then
  write(message, '(a)' ) &
&  ' vtowfk : prtvol=0 or 1, do not print more k-points.'
  call wrtout(06,message,'PERS')
 end if

!Print out decomposition of eigenvalues in the non-selfconsistent case
!or if prtvol>=10
 if( (iscf<0 .and. (prtvol>2 .or. ikpt<=nkpt_max)) .or. prtvol>=10)then

  write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&  ' mean kinetic energy (hartree) for',nband_k,'bands',ch10,&
&  '              after ',inonsc,' non-SCF iterations with ',&
&  dtset%nline,' CG line minimizations'
  call wrtout(06,message,'PERS')
  do ii=0,(nband_k-1)/6
   write(message, '(1p,6e12.4)' ) (ek_k(iband),&
&                         iband=1+6*ii,min(6+6*ii,nband_k))
   call wrtout(06,message,'PERS')
  end do

  write(message, '(5x,a,i5,2x,a,a,a,i4,a,i4,a)' ) &
&  ' mean non-local energy (hartree) for',nband_k,'bands',ch10,&
&  '              after ',inonsc,' non-SCF iterations with ',&
&  dtset%nline,' CG line minimizations'
  call wrtout(06,message,'PERS')

  if (gs_hamk%usepaw==0) then
   do ii=0,(nband_k-1)/6
    write(message, '(1p,6e12.4)' ) (enl_k(iband),&
&                        iband=1+6*ii,min(6+6*ii,nband_k))
    call wrtout(06,message,'PERS')
   end do
  end if

 end if

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')

 deallocate(evec,subham,gsc)

!DEBUG
!write(6,*)' vtowfk : 1 '
!ENDDEBUG
 if (gs_hamk%usepaw==0) then
  if (wfopta10==4) then
   deallocate(totvnl)
  else
   deallocate(subvnl)
  endif
 end if
 if (use_subovl==1) deallocate(subovl)
 if(wfoptalg==3)deallocate(eig_save)

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)') ch10,&
&   ' vtowfk : exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(30,2,tsec)
 call timab(28,2,tsec)

!DEBUG
!write(6,*)' vtowfk : exit '
!ENDDEBUG

end subroutine vtowfk
!!***
