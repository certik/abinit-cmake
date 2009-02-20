!{\src2tex{textfont=tt}}
!!****f* ABINIT/prctfvw1
!! NAME
!! prctfvw1
!!
!! FUNCTION
!! Compute new trial potential by applying the Thomas--Fermi--von Weizsaecker
!! charge mixing scheme (see PRB 64 121101).
!! First step is to compute a localy averaged non local potential
!! Written starting from src/5common/energy.F90
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! old trial potential
!! old output potential
!! new trial potential resulting from the mixing choice
!! old output density
!!
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of wavefunction
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eew=Ewald energy (hartree)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  eii=psp core-core energy
!!  entropy=entropy due to the occupation number smearing (if metal)
!!  epaw=PAW spherical part energy
!!  epawdc=PAW spherical part double-counting energy
!!  gsqcut=G^2 cutoff from gsqcut=ecut/(2 Pi^2)
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimension for number of planewaves
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#ngfft
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  nsym=number of symmetry elements in space group (at least 1)
!!  ntypat=number of types of atoms in cell
!!  n3xccc=dimension of the xccc3d array (0 or nfftf).
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each k point
!!  occopt=option for occupancies
!!  optene=option for the computation of total energy (direct scheme or double-counting scheme)
!!  pawfgr(natom) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  tsmear=smearing energy or temperature (if metal)
!!  vpsp(nfftf)=local pseudopotential in real space (hartree)
!!  wffnow=structured array giving all information about wavefunction file
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!! modified new trial potential based on the TFvW charge mixing
!!
!! SIDE EFFECTS
!!  rhog(2,nfftf)=work space for rho(G); save intact on return
!!  rhor(nfftf,nspden)=work space for rho(r); save intact on return
!!  nspinor should not be modified in the call of rdnpw
!!
!! WARNINGS
!! This is experimental code : input, ouptput, results and any other feature may vary greatly.
!!
!! NOTES
!!
!! PARENTS
!!      newvtr
!!
!! CHILDREN
!!      cgpr,dotprod_vn,fftpac,fourdp,fourwf,ftfvw1__end,ftfvw1__init,hdr_skip
!!      laplacian,leave_new,leave_test,mean_fftr,metric,mkffnl,nonlop,ph1d3d
!!      rdnpw,rhotov,rwwf,sphereboundary,timab,wrtout,xcomm_init,xdefineoff
!!      xmaster_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prctfvw1(atindx,atindx1,cg,deltae,densymop_gs,dtfil,dtset,eeig,&
     & efermi,eigen,ek,enl,etotal,fixmom,gsqcut,intxc,irrzon,ixc,&
     & kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nfft,nfftf,ngfftf,&
     & nhat,nhatgr,nhatgrdim,&
     & nkpt,nkxc,npwarr,nspden,nspinor,nsppol,nsym,ntypat,n3xccc,occ,occopt,optene,optres,optxc,&
     & pawang,pawfgr,pawfgrtab,pawtab,&
     & phnons,ph1d,psps,resid,rhog,rhor,rprimd,&
     & usexcnhat,&
     &vin_old,vout_unmixed,vpsp,vtrial,&
     & wffnow,xccc3d,xred,ylm,inlavnl)

 use defs_basis
  use defs_datatypes
  use ftfvw1


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_12spacepar
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_15common
 use interfaces_lib01cg
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

  !Arguments ------------------------------------
  !scalars
  integer,intent(in) :: intxc,ixc,mband,mgfft,mkmem,mpsang,mpw,n3xccc,natom,nfft
  integer,intent(in) :: nfftf,nhatgrdim,nkpt,nkxc,nspden,nsppol,nsym,ntypat
  integer,intent(in) :: occopt,optene,optres,optxc,usexcnhat
  integer,intent(inout) :: nspinor
  real(dp),intent(in) :: etotal,fixmom,gsqcut
  real(dp),intent(out) :: eeig,ek,enl
  type(MPI_type),intent(inout) :: mpi_enreg
  type(datafiles_type),intent(in) :: dtfil
  type(dataset_type),intent(in) :: dtset
  type(dens_sym_operator_type),intent(in) :: densymop_gs
  type(pawfgr_type),intent(in) :: pawfgr
  type(pseudopotential_type),intent(in) :: psps
  type(wffile_type),intent(inout) :: wffnow
  !arrays
  integer,intent(in) :: atindx(natom),atindx1(natom)
  !no_abirules
  !nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
  integer, intent(in)    :: irrzon(nfft**(1-1/nsym),2,nspden/nsppol),kg(3,mpw*mkmem)
  integer, intent(in)    :: nattyp(ntypat),ngfftf(18),npwarr(nkpt)
  real(dp), intent(in)   :: cg(2,mpw*nspinor*mband*mkmem*nsppol),eigen(mband*nkpt*nsppol)
  real(dp), intent(in)   :: inlavnl(nfftf,nspden)
  ! WARNING
  ! BEWARE THERE IS TWO DIFFERENT SIZE DECLARED FOR ARRAY NHAT IN RHOTOV AND RHOHXC
  ! THIS MIGHT RESULT IN A BUG
  real(dp),intent(in)   ::   nhat(nfftf,nspden*psps%usepaw),nhatgr(nfftf,nspden,3*nhatgrdim)
  real(dp), intent(in)   :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
  !nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
  real(dp), intent(in)   :: phnons(2,nfft**(1-1/nsym),nspden/nsppol)
  real(dp), intent(out)  :: resid(mband*nkpt*nsppol)
  real(dp), intent(inout):: rhog(2,nfftf),rhor(nfftf,nspden)
  real(dp), intent(in)   :: rprimd(3,3),xred(3,natom)
  real(dp), intent(inout):: vin_old(nfftf,nspden),vout_unmixed(nfftf,nspden),vtrial(nfftf,nspden)
  real(dp), intent(inout):: vpsp(nfftf)
  real(dp),dimension(:),intent(inout)    :: xccc3d(n3xccc)
  real(dp), intent(in)   :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
  type(pawtab_type), intent(in)  :: pawtab(ntypat*psps%usepaw)
  type(pawang_type),intent(in) :: pawang
  type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
  real(dp), intent(in) :: deltae,efermi
  !Local variables-------------------------------
  !scalars
  integer,parameter :: nl=500
  integer :: bdtot_index,choice,count,count1,count10,count11,count12,count13
  integer :: count14,count15,count16,count17,count18,count19,count2,count3
  integer :: count4,count5,count6,count7,count8,count9,cplex,cpopt,dimdij,dimffnl
  integer :: formeig,i1,i2,i3,ia,iatom,iband,icg,id1,id2,id3,ider,idir,ierr,ifft
  integer :: ifft2,ig,ig1,ig2,ig3,ii,ikg,ikpt,il,ilm,ilmn,index,isign,ispden
  integer :: isppol,istwf_k,izero,master,matblk,mcg_disk,me,muig,n1,n2,n3,n4,n5
  integer :: n6,nband_k,nfftotf,nkpg,nnlout,npw_k,option,paw_opt,signs,spaceComm
  integer :: tim_fourdp,tim_fourwf,tim_mkrho,tim_nonlop,tim_rwwf
  integer :: me_distrb
  real(dp),parameter :: alpha32=(3._dp*pi*pi)
  real(dp) :: Z,arg,cemm,cepp,dl,doti,dotr,dotr2,dummy,dummy2,e00,eei,eeigk
  real(dp) :: ehart,emm,enlk,enxc,enxcdc,epp,lambda,ucvol,vme,vres2,vxcavg
  real(dp) :: weight
  character(len=500) :: message
  character(len=fnlen) :: filapp
  type(gs_hamiltonian_type) :: gs_hamk
  type(energies_type) :: energies
  !arrays
  integer :: typat(natom)
  integer,allocatable :: kg_dum(:,:),kg_k(:,:)
  real(dp) :: aa(2),bb(2),bla(nfftf,nspden),cc(2),d1e(nfftf,nspden)
  real(dp) :: d2e(nfftf,nspden),enlout(1),gmet(3,3),gprimd(3,3),kpoint(3)
  real(dp) :: nonlop_dum(1,1),phi1(nfftf,nspden),phi1p(nfftf,nspden)
  real(dp) :: phiout(nfftf,nspden),phioutg(2,nfftf),qphon(3),rhodum(1),rmet(3,3)
  real(dp) :: strsxc(6),te(nl),tl(nl),tsec(2),vhphi(nfftf),vhphi2(nfftf)
  real(dp) :: vhphi3(nfftf),vmean(2),vnew_mean(2),vref(nfftf,nspden)
  real(dp) :: vrefg(2,nfftf,nspden),vres_mean(2),vtrialold(nfftf,nspden)
  real(dp) :: xredcp(3,natom),ylmgr_dum(1),znucl(1)
  real(dp),allocatable :: buffer(:),cg_disk(:,:),cgrvtrial(:,:),cwavef(:,:)
  real(dp),allocatable :: deltaW(:,:),dummyt2(:,:),dummyt3(:,:,:),eig_dum(:)
  real(dp),allocatable :: eig_k(:),ffnl(:,:,:,:),g2cart(:),kinpw(:),kpg_dum(:,:)
  real(dp),allocatable :: kxc(:,:),laplacerhor(:,:),lavnl(:,:,:),lavnlfft(:,:)
  real(dp),allocatable :: newvout(:,:),newvoutfourier(:,:,:),occ_dum(:),occ_k(:)
  real(dp),allocatable :: ph3d(:,:,:),resid_k(:),sqrtrhor(:,:),vhartr(:)
  real(dp),allocatable :: vin_oldfourier(:,:,:),vlocal(:,:,:,:),vnlcwavef(:,:)
  real(dp),allocatable :: vnlwfraug(:,:,:,:),vresid(:,:),vtrialfourier(:,:,:)
  real(dp),allocatable :: vxc(:,:),wfraug(:,:,:,:),ylm_k(:,:)
  type(cprj_type) :: cprj_dum(1,1)
! *************************************************************************
!***********************************************************************************
!Getting the localy averaged non-local potential                                ***
!$Vnl(r) = [\sum_{n,k} f_{n,k} \psi_{n,k}(r) (Vnl(r,r') |\psi_{n,k}(r')>)]/n(r)$***
!**********************************************************************************
 znucl=1.0_dp
 typat=1
 xredcp = xred
 qphon=zero
!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 nfftotf=ngfftf(1)*ngfftf(2)*ngfftf(3)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.nfft/=nfftf)) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' prctfw :  BUG -',ch10,&
&  '  wrong values for nfft, nfftf !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 call timab(59,1,tsec)
!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)
!PATCH prctfw1 // KPT & FFT me-->me_kpt & spaceComm --> comm_kpt
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
& (mpi_enreg%paral_compil_fft==1)) then
  me_distrb = mpi_enreg%me_kpt
  spaceComm = mpi_enreg%comm_kpt
 else
  me_distrb = mpi_enreg%me
 end if
!Init master
 call xmaster_init(mpi_enreg,master)
!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 if (mkmem==0) then
! Read wavefunction file header
  call hdr_skip(wffnow,ierr)
! Define offsets, in case of MPI I/O
  formeig=0
  call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,nspinor,nsppol,nkpt)
  mcg_disk=mpw*nspinor*mband
  allocate(cg_disk(2,mcg_disk))
 end if
 eeig=zero
 ek=zero
 enl=zero
 bdtot_index=0
 icg=0
!DEBUG
!write(6,*)' prctfw : before loop over spins '
!stop
!ENDDEBUG
 n1=dtset%ngfft(1)
 n2=dtset%ngfft(2)
 n3=dtset%ngfft(3)
 n4=dtset%ngfft(4)
 n5=dtset%ngfft(5)
 n6=dtset%ngfft(6)
!nvloc=1
!if(nspden==4)nvloc=4
!allocate(vlocal(n4,n5,n6,nvloc);
 allocate(kg_k(3,mpw),cwavef(2,mpw*nspinor), vnlcwavef(2,mpw*nspinor))
 allocate(wfraug(2,n4,n5,n6),vnlwfraug(2,n4,n5,n6),lavnl(n4,n5,n6))
 wfraug(:,:,:,:)=zero
 vnlwfraug(:,:,:,:)=zero
 lavnl(:,:,:)=zero
 allocate(lavnlfft(nfftf,nspden),g2cart(nfftf))
!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 allocate(gs_hamk%atindx(natom),gs_hamk%atindx1(natom))
 allocate(gs_hamk%gbound(2*mgfft+8,2))
 allocate(gs_hamk%indlmn(6,psps%lmnmax,ntypat))
 allocate(gs_hamk%nattyp(ntypat))
 allocate(gs_hamk%phkxred(2,natom))
 allocate(gs_hamk%ph1d(2,3*(2*mgfft+1)*natom))
 allocate(gs_hamk%pspso(ntypat))
 allocate(gs_hamk%xred(3,natom))
!Initialize most of the Hamiltonian
 gs_hamk%atindx(:)  =atindx(:)
 gs_hamk%atindx1(:) =atindx1(:)
 gs_hamk%gmet(:,:)  =gmet(:,:)
 gs_hamk%gprimd(:,:)=gprimd(:,:)
 gs_hamk%indlmn(:,:,:)=psps%indlmn(:,:,:)
 gs_hamk%lmnmax     =psps%lmnmax
 gs_hamk%mgfft      =mgfft
 gs_hamk%mpsang     =mpsang
 gs_hamk%mpssoang   =psps%mpssoang
 gs_hamk%natom      =natom
 gs_hamk%nattyp(:)  =nattyp(:)
 gs_hamk%nfft       =nfft
 gs_hamk%ngfft(:)   =dtset%ngfft(:)
 gs_hamk%nloalg(:)  =dtset%nloalg(:)
 gs_hamk%nspinor      =nspinor
 gs_hamk%ntypat      =ntypat
!gs_hamk%nvloc      =nvloc
 gs_hamk%n4         =n4
 gs_hamk%n5         =n5
 gs_hamk%n6         =n6
 gs_hamk%usepaw     =psps%usepaw
 gs_hamk%ph1d(:,:)  =ph1d(:,:)
 gs_hamk%pspso(:)   =psps%pspso(:)
 gs_hamk%ucvol      =ucvol
 gs_hamk%useylm     =psps%useylm
 gs_hamk%xred(:,:)  =xred(:,:)
!Special case of array ekb:
!Not the same meaning for norm-cons. psps and paw
 gs_hamk%dimekb1=psps%dimekb
 gs_hamk%dimekb2=ntypat
 allocate(gs_hamk%ekb(psps%dimekb,ntypat,nspinor**2))
 gs_hamk%ekb(:,:,1)=psps%ekb(:,:)
 if (nspinor==2) then
  gs_hamk%ekb(:,:,2)=psps%ekb(:,:)
  gs_hamk%ekb(:,:,3:4)=zero
 end if
 allocate(gs_hamk%sij(gs_hamk%dimekb1,ntypat))
 lavnl(:,:,:)=zero
!LOOP OVER SPINS
 do isppol=1,nsppol
! Rewind kpgsph data file if needed:
  if (mkmem==0) rewind dtfil%unkg
  if (mkmem==0.and.psps%useylm==1) rewind dtfil%unylm
  ikg=0
! here removed the part about vlocal of energy.F90 (L377 - L393)
! Loop over k points
! write(0,*) 'DEBUG 1'
  ikpt_lavnl: do ikpt=1,nkpt
   nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
   istwf_k=dtset%istwfk(ikpt)
   npw_k=npwarr(ikpt)
   if(mpi_enreg%paral_compil_kpt==1)then
!   Skip this k-point if not the proper processor
    if (mpi_enreg%parareel == 0) then
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&     -me_distrb))/=0) then
      resid(1+bdtot_index : nband_k+bdtot_index) = zero
      bdtot_index=bdtot_index+nband_k
      cycle
     end if
    else
     if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt) &
&     /= mpi_enreg%me) then
      resid(1+bdtot_index : nband_k+bdtot_index) = zero
      bdtot_index=bdtot_index+nband_k
      cycle
     end if
    end if
   end if
!  bdtot_index=bdtot_index+nband_k     ! added this thinking that bdtot should always be increased
!  Continue to initialize the Hamiltonian
   gs_hamk%istwf_k    =istwf_k
   gs_hamk%npw        =npw_k
   allocate(eig_k(nband_k),occ_k(nband_k),resid_k(nband_k))
   allocate(ylm_k(npw_k,mpsang*mpsang*psps%useylm))
   resid_k(:)=real(0.0,dp)
   kpoint(:)=dtset%kptns(:,ikpt)
   gs_hamk%kpoint(:)  =kpoint(:)
!  write(0,*) 'prctfw kpt,bdtot',ikpt,bdtot_index
   occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
   eig_k(:)=eigen(1+bdtot_index:nband_k+bdtot_index)
   if (minval(eig_k)>1.d100) eig_k=zero
   if (mkmem==0) then
!   Read sphere data centered at k in dtfil%unkg, then k+g data
    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,dtfil%unkg)
    read (dtfil%unkg) kg_k(1:3,1:npw_k)
    call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,mgfft,npw_k)
!   Eventually read spherical harmonics
    if (psps%useylm==1) then
     read(dtfil%unylm)
     read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
    end if
!   Read the wavefunction block for ikpt,isppol
    tim_rwwf=3
    allocate(eig_dum(mband),kg_dum(3,0),occ_dum(mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,&
&    nband_k,nband_k,npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
    deallocate(eig_dum,kg_dum,occ_dum)
   else
    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,mgfft,npw_k)
    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
      ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     end do
    end if
   end if  !  End if for choice governed by mkmem
!  write(0,*) 'DEBUG 2'
   enlk=zero
   eeigk=zero
   index=0
   index=1
   ider=0;dimffnl=1;nkpg=0
   allocate(ffnl(npw_k,dimffnl,psps%lmnmax,ntypat))
!  call mkffnl(psps%dimekb,dimffnl,dtset%effmass,psps%ekb,ffnl,psps%ffspl,&
!  &   gmet,gprimd,ider,ider,psps%indlmn,kg_k,kpoint,psps%lmnmax,&
!  &   psps%lnmax,mpsang,psps%mqgrid_ff,&
!  &   npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
!  &   psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&   gmet,gprimd,ider,ider,psps%indlmn,kg_k,kpg_dum,kpoint,psps%lmnmax,&
&   psps%lnmax,mpsang,psps%mqgrid_ff,nkpg,&
&   npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&   psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
   do ia=1,natom
    iatom=atindx(ia)
    arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
    gs_hamk%phkxred(1,iatom)=cos(arg)
    gs_hamk%phkxred(2,iatom)=sin(arg)
   end do
   if(dtset%nloalg(1)<=0)then
!   Only the allocation, not the precomputation.
    matblk=dtset%nloalg(4)
    allocate(ph3d(2,npw_k,matblk))
   else
!   Here, allocation as well as precomputation
    matblk=natom
    allocate(ph3d(2,npw_k,matblk))
    call ph1d3d(1,natom,kg_k,kpoint,matblk,natom,npw_k,n1,n2,n3,&
&    gs_hamk%phkxred,ph1d,ph3d)
   end if
   gs_hamk%matblk=matblk
!  Compute nonlocal psp energy - Norm-conserving only
   usepaw_lavnl: if (psps%usepaw==0) then
    iband_lavnl: do iband=1,nband_k
     if(mpi_enreg%paral_compil_kpt==1)then
!     Skip this band if not the proper processor
      if (mpi_enreg%paralbd >1) then
       if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= mpi_enreg%me) then
!       index=index+npw_k*nspinor
        cycle
       end if
      end if
     end if
     if(nspinor > 1) stop 'prctfw is not already written to handle nspinor>1'
     if(mkmem/=0)cwavef(:,1:npw_k*nspinor)=&
&     cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
     if(mkmem==0)cwavef(:,1:npw_k*nspinor)=&
&     cg_disk(:,1+(iband-1)*npw_k*nspinor:iband*npw_k*nspinor)
!    if(mkmem/=0)then
!    cwavef(:,1:npw_k*nspinor)=cg(:,index:index+npw_k*nspinor)
!    index=index+npw_k*nspinor
!    else
!    cwavef(:,1:npw_k*nspinor)=cg_disk(:,index:index+npw_k*nspinor)
!    index=index+npw_k*nspinor
!    end if
     choice=1  !1:  a non-local energy contribution
     signs=2   !2: compute |out> = Vnl|in>:
     idir=0    !meaning less in our case
     nnlout=1  !depend on paw_opt and choice (see nonlop.F90), meaning less for signs=2
     tim_nonlop=3 !timing code (doesn't matter)
     paw_opt=gs_hamk%usepaw;cpopt=-1
!    get the function out> = Vnl|in> where |in> are the wave funtion cwavef
     call nonlop(atindx1,&      ! invers of index table for atoms
&    choice,&             ! ==1 =>a nl energy contribution
&    cpopt,&              ! -1 no treatment of cprj scalars
&    cprj_dum,&           ! wave function projected with non-local projectors
&    gs_hamk%dimekb1,&    ! dimenl1 (lmnmax for norm-conserving)
&    gs_hamk%dimekb2,&    ! dimenl2 (ntypat for normconserving)
&    dimffnl,&            ! dimffnlin
&    dimffnl,&            ! dimffnlout
&    gs_hamk%ekb,&   ! enl(  dimenl1,     dimenl2)
&    enlout,&             ! enlout(nnlout)
&    ffnl,ffnl,&          !  ffnlin(),  ffnlout(npwin,dimffnlin,lmnmax,ntypat)
!    =nonlocal form factors to be used
!    for the application of the nonlocal operator to the |in> (resp. |out>) vector
&    gmet,&               ! gmet(3,3) metric tensor
&    gprimd,&             ! gprimd(3,3)=dimensional reciprocal space primitive translations
&    idir,&               ! meaningless for choice==1
&    psps%indlmn,&        ! indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s
&    istwf_k,&            ! option parameter that describes the storage of wfs
&    kg_k,&          ! kgin(3,npwin)=integer coords of planewaves in basis sphere, for the |in> vector
&    kg_k,&               ! kgout(3,npwout)=integer coords of planewaves in basis sphere, for the |out> vector
&    kpg_dum,&
&     kpg_dum,&
&     kpoint,&             ! kptin(3)=k point in terms of recip. translations, for the |in> vector
&    kpoint,&             ! kptout(3)=k point in terms of recip. translations, for the |out> vector
&    dummy,&              ! Lambda in Vnl-lambda.S - not used here
&    psps%lmnmax,&        ! lmnmax=max. number of (l,m,n) components over all types of atoms
&    matblk,&             ! matblk=dimension of the arrays ph3din and ph3dout
&    mgfft,&              ! mgfft=maximum size of 1D FFTs
&    mpi_enreg,&          ! mpi_enreg=informations about MPI parallelization
&    psps%mpsang,&        ! mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
&    psps%mpssoang,& ! mpssoang= 1+max(spin*angular momentum) for nonlocal pseudopotentials
&    natom,nattyp,&       ! natom=number of atoms in cell;  nattyp(ntypat)=number of atoms of each type
&    dtset%ngfft,&        ! ngfft(18)
&    nkpg,&
&     nkpg,&
&     dtset%nloalg,&       ! nloalg(5)=governs the choice of the algorithm for nonlocal operator
&    nnlout,&             ! nnlout=dimension of enlout (choice==1.and.paw_opt<2)=>nnlout=1)
&    npw_k,&              ! npwin=number of planewaves for given k point, for the |in> vector
&    npw_k,&              ! npwout=number of planewaves for given k point, for the |out> vector
&    nspinor,&       ! nspinor=number of spinorial components of the wavefunctions
&    ntypat,&             ! ntypat=number of types of atoms in cell
&    0,  &                ! only_SO = flag to get only SO part of V_NL
&    paw_opt,&            ! paw_opt= define the nonlocal operator concerned with:
!    paw_opt==0 : Norm-conserving Vnl
&    gs_hamk%phkxred,&    ! phkxredin(2,natom)=phase factors exp(2 pi kptin.xred)
&    gs_hamk%phkxred,&    ! phkxredout(2,natom)=phase factors exp(2 pi kptout.xred)
&    ph1d,&               ! ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
&    ph3d,&          ! ph3din(2,npwin,matblk)=3D structure factors, for each atom and plane wave (in)
&    ph3d,&               ! ph3dout(2,npwout,matblk)=3-dim struct fact, for each atom and plane wave (out)
&    psps%pspso,&         ! pspso(ntypat)=spin-orbit characteristic for each atom type
&    signs,&              ! 1=>contracted elements (like enl) 2=>work on a function in recip. space
&    nonlop_dum,&         ! sij (used only for paw_opt>2)
&    nonlop_dum,&         ! svectout(2,nspinor*npwout*(paw_opt/3))  (used only for paw_opt>2)
&    tim_nonlop,&
&     ucvol,&              ! unit cell volume (bohr^3)
&    psps%useylm,&   ! useylm=governs the way the nonlocal operator is to be applied:
&    cwavef,&             !  vectin(2,nspinor*npwin)=input cmplx wavefunction coefficients <G|Cnk>
&    vnlcwavef)           ! vectout(2,nspinor*npwout)=result of the aplication of the nl operator
!    enlk=enlk+occ_k(iband)*enlout(1)
!    eeigk=eeigk+occ_k(iband)*eig_k(iband)
!    ******************************************
!    getting the new vnl in the way of mkrho3*
!    ******************************************
     tim_fourwf=0 ! timing code for the routine : dunno
!    here wfraug is the modified variable
     call fourwf(1,dummyt3,cwavef,dummyt2,wfraug,gs_hamk%gbound,gs_hamk%gbound,&
&     istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,dummy)
     call fourwf(1,dummyt3,vnlcwavef,dummyt2,vnlwfraug,gs_hamk%gbound,gs_hamk%gbound,&
&     istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamk%ngfft,npw_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,dummy)
     weight=occ_k(iband)*dtset%wtk(ikpt)/ucvol
!    here wfraug are the wavefunction in realspace (on the augmented grid?)
!    here vnlwfraug are the projection of Vnl over the wavefunction (on the augmented grid?)
     lavnl(:,:,:)=lavnl(:,:,:)+weight* (&
&     wfraug(1,:,:,:) * vnlwfraug(1,:,:,:) &
&     +wfraug(2,:,:,:) * vnlwfraug(2,:,:,:) )
    end do iband_lavnl ! nband_k
   end if usepaw_lavnl !PAW or norm-conserving
!  deallocate(eig_k,occ_k,resid_k)
!  deallocate(kinpw)
   deallocate(eig_k,occ_k,resid_k)
   deallocate(ffnl,ph3d)
   deallocate(ylm_k)
   bdtot_index=bdtot_index+nband_k
   if (mkmem/=0) then
!   Handle case in which kg, cg, are kept in core
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   end if
  end do ikpt_lavnl !ikpt
 end do !isppol
 if(mpi_enreg%paral_compil_kpt==1)then
! Accumulate lavnl from each process
  allocate(buffer(n1*n2*n3))
  do i3=1,n3
   do i2=1,n2
    do i1=1,n1
     buffer(i1+n1*(i2-1+n2*(i3-1)))=lavnl(i1,i2,i3) !/rhor(1,i1,i2,i3)
    end do
   end do
  end do
  call timab(48,1,tsec)
  call xsum_mpi(buffer,spaceComm,ierr)
  call timab(48,2,tsec)
  do i3=1,n3
   do i2=1,n2
    do i1=1,n1
     lavnl(i1,i2,i3)=buffer(i1+n1*(i2-1+n2*(i3-1)))
    end do
   end do
  end do
  deallocate(buffer)
  option = 1 ! transfer the augmented grid lavnl to the  fftpacked one
  do ispden=1,nspden
!  WARNING this may work only for nspden == 1
   call fftpac(ispden,nspden,n1,n2,n3,n4,n5,n6,ngfftf,lavnlfft,lavnl,option)
   do ifft=1,nfftf
    lavnlfft(ifft,ispden) = lavnlfft(ifft,ispden)/rhor(ifft,ispden)
   end do
  end do
  call leave_test(mpi_enreg)
! write(message,*) 'prctfw: loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
! write(0,*) 'prctfw should not be there in sequential run'
 else !non parallel case
  option = 1 ! transfer the augmented grid lavnl to the  fftpacked one
  do ispden=1,nspden
!  WARNING this may work only for nspden == 1
   call fftpac(ispden,nspden,n1,n2,n3,n4,n5,n6,ngfftf,lavnlfft,lavnl,option)
!  i1=1
!  i2=1
!  do i3=1,n3
!  ifft=i1+n1*(i2-1+n2*(i3-1))
!  write(80,*) lavnlfft(ifft,1)
!  end do
   do ifft=1,nfftf
    lavnlfft(ifft,ispden) = lavnlfft(ifft,ispden)/(rhor(ifft,ispden)+0.0001_dp)
   end do
!  i1=1
!  i2=1
!  do i3=1,n3
!  ifft=i1+n1*(i2-1+n2*(i3-1))
!  write(81,*) lavnlfft(ifft,1)
!  write(82,*) rhor(ifft,1)
!  end do
  end do
 end if
!write(0,*) 'DEBUG 5'
 deallocate(gs_hamk%atindx,gs_hamk%atindx1)
 deallocate(gs_hamk%ekb)
 deallocate(gs_hamk%gbound)
 deallocate(gs_hamk%indlmn)
 deallocate(gs_hamk%nattyp)
 deallocate(gs_hamk%phkxred)
 deallocate(gs_hamk%pspso)
 deallocate(gs_hamk%ph1d)
 deallocate(gs_hamk%sij)
 deallocate(gs_hamk%xred)
!write(0,*) 'DEBUG 5a'
 if(mkmem==0)then
  deallocate(cg_disk)
 end if
!write(0,*) 'DEBUG 5b'
!DEBUG
!write(6,*)' energy : after loop on kpts and spins '
!stop
!ENDDEBUG
 write(message, '(4a)' )ch10, &
& ' prctfw: COMMENT -',ch10,&
& '  New local average of non local potential made from  wfs and Vnl'
 call wrtout(06,message,'COLL')
 call timab(59,2,tsec)
 deallocate(cwavef)
 deallocate(kg_k)
 deallocate(vnlcwavef)
 deallocate(wfraug)
 deallocate(vnlwfraug)
 deallocate(lavnl)
!******************************************************************
!Getting the DeltaW factor                                      **
!******************************************************************
!first get the laplacian of the density
!the reciprocal density is complex !
 allocate(deltaW(nfftf,nspden),laplacerhor(nfftf,nspden),sqrtrhor(nfftf,nspden))
 allocate(vtrialfourier(2,nfftf,nspden),newvoutfourier(2,nfftf,nspden))
!Compute the real space laplacian of sqrt(rhor)
 sqrtrhor(:,:)=(rhor(:,:))**half
 call laplacian(gprimd,mpi_enreg,nfftf,nspden,ngfftf,dtset%paral_kgb,rdfuncr=sqrtrhor,laplacerdfuncr=laplacerhor,g2cart_out=g2cart)
!second step: get deltaW
 do ispden=1,nspden
  do ifft=1,nfftf
   deltaW(ifft,ispden)=((efermi&
&   - half*(alpha32 * (rhor(ifft,ispden)))**(two_thirds)&
&   - vin_old(ifft,ispden) - lavnlfft(ifft,ispden))& ! this one is the true one
&  * sqrtrhor(ifft,ispden)&
&   + one*half*laplacerhor(ifft,ispden))
  end do
 end do
!write(7779,*)lavnl
!write(7777,*)lavnlfft
!******************************************************************
!Finding the density which minimizes the associated Energy      **
!******************************************************************
!compute the total charge number
!first switch to the sqrt of the density...
 cplex=1;
 option=1;
 call dotprod_vn(cplex,& !complex density/pot
&sqrtrhor,&          !the density
&Z,&  !resulting dorproduct integrated over r
&doti,&          !imaginary part of the integral
&mpi_enreg,&     !
&size(rhor,1),&          !number of localy(cpu) attributed grid point
&nfftotf,&        !real total number of grid point
&size(rhor,2),&        !nspden
&option,&        !1=compute only the real part 2=compute also the imaginary part
&sqrtrhor,&          !the potential
&ucvol)          !cell volume
!enable the use of the functions eneofrho_tfw and deneofrho_tfw
 Z=real(nint(Z),dp)
 call ftfvw1__init(dtset,dtset%intxc,dtset%ixc,psps%usepaw,n3xccc,ngfftf,&
& nfftf,&
& nhat,nhatgr,nhatgrdim,&
& nkxc,nspden,mpi_enreg,deltaW,gprimd,gsqcut,&
& lavnlfft,rhor,rprimd,ucvol,&
& psps%usepaw,usexcnhat,&
& vout_unmixed,vpsp,vtrial,xccc3d,Z)
!minimizes Etfw with respect to sqrtrhor instead of rhor
!sqrtrhor(:,:)=two*Z/nfftf
!call random_number(rhor)
!eei=zero
!do ifft=1,nfft
!eei=max(eei,sqrtrhor(ifft,1))
!end do
!eei=0.05_dp*eei
!call newdensity(eei,rhor,sqrtrhor)
 phiout=sqrtrhor
 call cgpr(size(rhor,1),size(rhor,2),ftfvw1__e,ftfvw1__de,ftfvw1__newdensity,&
& abs(deltae*real(0.001,dp)/etotal),155,sqrtrhor,dummy,dummy2)
!call random_number(sqrtrhor)
 phi1=sqrtrhor-phiout
!dummy=eneofrho_tfw(nfftf,nspden,sqrtrhor)
!free the dynamically allocated memory used by eneofrho_tfw and deneofrho_tfw
 call ftfvw1__end()
!new density from the minimised sqrtrhor
 count=0
 do ifft=1,nfftf
  if (sqrtrhor(ifft,1)<zero) then
   count=count+1
  end if
 end do
 write(6,*) 'prctfw (04rsprc) : number of element in sqrt(density(nfft)) < 0 :',count
 rhor(:,:)=sqrtrhor(:,:)*sqrtrhor(:,:)
!******************************************************************
!Production of a new trial potential                            **
!******************************************************************
!write(3330,*) vtrial
 allocate(newvout(nfftf,nspden),vin_oldfourier(2,nfftf,nspden))
!step1: make V from rho
 call fourdp(1, rhog(:,:), rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0) ! rhor(:,1) is the full density!
 allocate(vhartr(nfftf),vxc(nfftf,nspden),vresid(nfftf,nspden),kxc(nfftf,1))
 call rhotov(dtset,energies,gsqcut,kxc,mpi_enreg,nfftf,ngfftf,&
!bug in test v5/1       &             nhat(:,1:nspden*psps%usepaw),&
& nhat,&
& nhatgr,nhatgrdim,  &
& nkxc,vresid,&
& n3xccc,optene,1,optxc,pawang,pawfgrtab,pawtab,&
& rhog,rhor,rprimd,strsxc,ucvol,&
& psps%usepaw,usexcnhat,vhartr,vnew_mean,vpsp,&
& vres_mean,vres2,newvout,vxcavg,vxc,xccc3d)
 call mean_fftr(newvout,vmean,mpi_enreg,nfftf,nfftotf,min(nspden,2))
 770 FORMAT(A,I3.3)
 deallocate(kxc)
!shift the newvout to remove the constant part. Don't know if this is OK
!WARNING
 do ispden=1,min(nspden,2)
  if(nspden/=2 .or. &
&  ( occopt>=3 .and. abs(fixmom+real(99.99,dp))<real(1.0d-10,dp) ))then
   if(nspden==1)then
    vme=vmean(1)
   else
    vme=(vmean(1)+vmean(2))*half
   end if
  else
   vme=vmean(ispden)
  end if
  newvout(:,ispden)=newvout(:,ispden)-vme
 end do
 deallocate(vhartr,vxc,vresid)
!step2: mix the new Vout with the old Vout
!mixing based on the kinetic energy at the k points
!change potential to fourier space
 do ispden=1,nspden
  call fourdp(1, newvoutfourier(:,:,ispden), newvout(:,ispden),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
  call fourdp(1, vtrialfourier(:,:,ispden), vtrial(:,ispden),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
! call fourdp(dtset%paral_kgb,1, vin_oldfourier(:,:,ispden), vtrial(:,ispden),-1,nfftf,ngfftf,0)     !modified for printing...
 end do
!filtering
 do ispden=1,nspden
  do ifft=1,nfftf
!  vtrialfourier(:,ifft,ispden) = newvoutfourier(:,ifft,ispden)*exp(-g2cart(ifft)*two) &
!  &+(one-exp(-g2cart(ifft)*two))*(&
!  &vin_oldfourier(:,ifft,ispden)*(one-exp(0.8_dp*g2cart(ifft)/(g2cart(ifft)+half)))& !! mixing proposed by
!  &+vtrialfourier(:,ifft,ispden)*(exp(0.8_dp*g2cart(ifft)/(g2cart(ifft)+half))))    !! raczowski...
!  & vtrialfourier(:,ifft,ispden)) !no further mixing
!  vtrialfourier(:,ifft,ispden) = newvoutfourier(:,ifft,ispden)*exp(-g2cart(ifft)*12.56637061435917_dp) &
!  &+(one-exp(-g2cart(ifft)*12.56637061435917_dp))*(&
!  & vtrialfourier(:,ifft,ispden))
   vtrialfourier(:,ifft,ispden) = newvoutfourier(:,ifft,ispden)*exp(-g2cart(ifft)*2_dp) &
&   +(one-exp(-g2cart(ifft)*2_dp))*(&
&   vtrialfourier(:,ifft,ispden))
  end do
 end do
!change resulting potential to real space
!write(3363,*) vtrialfourier
 vtrialold=vtrial
 do ispden=1,nspden
! call fourdp(1,vtrialfourier(:,:,ispden),newvout(:,ispden),1,nfftf,ngfftf,dtset%paral_kgb,0)    !output modified->no unholy effect
  call fourdp(1,vtrialfourier(:,:,ispden),vtrial(:,ispden),1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
 end do
 vtrial=vtrialold*real(0.0,dp)+vtrial*real(1.0,dp)
!----------------------------------------------------------------------------------------------------------
!last FREE
 write(6,*) 'prctfw (04rsprc) : FREEING THE REMAINING OF MEMORY FROM PRCTFW'
 deallocate(lavnlfft)
 deallocate(deltaW)
 deallocate(newvout)
 deallocate(newvoutfourier)
 deallocate(vtrialfourier)
 deallocate(sqrtrhor)
 deallocate(laplacerhor)
 deallocate(g2cart)
 deallocate(vin_oldfourier)
end subroutine prctfvw1
!!***
