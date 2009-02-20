!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv3
!! NAME
!! scfcv3
!!
!! FUNCTION
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! optimum and optionally to compute mixed derivatives of energy.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, DRH, MB, XW, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=pw coefficients of GS wavefunctions at k.
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  cprj(dimpaw1,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  cpus= cpu time limit in seconds
!!  dimcprj(natom)=array of dimensions of arrays cprj, cprjq
!!  dimpaw1= -PAW only- dimension of 1st-order on-site quantities (rhoij, Dij...)
!!  densymop_rf <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (response-function)
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eew=2nd derivative of Ewald energy (hartree)
!!  efrhar=Contribution from frozen-wavefunction, hartree energy,
!!           to the second-derivative of total energy.
!!  efrkin=Contribution from frozen-wavefunction, kinetic energy,
!!           to the second-derivative of total energy.
!!  efrloc=Contribution from frozen-wavefunction, local potential,
!!           to the second-derivative of total energy.
!!  efrnl=Contribution from frozen-wavefunction, non-local potential,
!!           to the second-derivative of total energy.
!!  efrx1=Contribution from frozen-wavefunction, xc core correction(1),
!!           to the second-derivative of total energy.
!!  efrx2=Contribution from frozen-wavefunction, xc core correction(2),
!!           to the second-derivative of total energy.
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  eii=2nd derivative of pseudopotential core energy (hartree)
!!  fermie=fermi energy (Hartree)
!!  fform=integer specifying data structure of wavefunction files
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  idir=direction of the current perturbation
!!  indkpt1(nkpt_rbz)=non-symmetrized indices of the k-points
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,nspden/nsppol)=irreducible zone data for RF symmetries
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates at k
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  kxc(nfftf,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mgfftf=maximum size of 1D FFTs for the "fine" grid (see NOTES in respfn.F90)
!!  mkmem =number of k points which can fit in memory (GS data); 0 if use disk
!!  mkqmem =number of k+q points which can fit in memory (GS data); 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpert=maximum number of ipert
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw for wfs at k.
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point, for each polarization
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  ngfftf(1:18)=integer array with FFT box dimensions and other for the "fine" grid (see NOTES in respfn.F90)
!!  nhat(cplex*nfftf,nspden*usepaw)= -PAW only- compensation density
!!  nkpt=number of k points in the full BZ
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array.
!!  mpi_enreg=informations about MPI parallelization
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsym1=number of symmetry elements in space group consistent with perturbation
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used
!!   otherwise, cplex*nfftf
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k point of the reduced Brillouin zone.
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh for the GS
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pertcase=fuill index of the perturbation
!!  phnons1(2,nfft**(1-1/nsym1),nspden/nsppol)=nonsymmorphic transl. phases,
!!   for RF symmetries
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1df(2,3*(2*mgfftf+1)*natom)=one-dimensional structure factor information for the "fine" grid
!!  prtbbb=if 1, band-by-band decomposition (also dim of d2bbb)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symaf1(nsym1)=anti(ferromagnetic) part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  symrl1(3,3,nsym1)=symmetry operations in real space in terms
!!   of primitive translations
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  tnons1(3,nsym1)=nonsymmorphic translations for symmetry operations
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  wffddk=struct info for ddk file
!!  wffnew=struct info for 1WF at exit, if mk1mem=0
!!  wffnow=struct info for 1WF at start, if mk1mem=0
!!  wfftgs=struct info for GS WF at start, if mkmem=0
!!  wfftkq=struct info for k+q GS WF at start, if mkqmem=0
!!  vpsp1(cplex*nfftf)=first-order derivative of the ionic potential
!!  vtrial(nfftf,nspden)=GS potential (Hartree).
!!  wtk_rbz(nkpt_rbz)=weight for each k point in the reduced Brillouin zone
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics at k
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics at k+q
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!!  d2lo(2,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,mpert,3,mpert)=non-local contributions to the 2DTEs
!!  eberry=energy associated with Berry phase
!!  edocc=correction to 2nd-order total energy coming from changes of occupation
!!  eeig0=0th-order eigenenergies part of 2nd-order total energy
!!  ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!    for strain perturbation only (zero otherwise, and not used)
!!  ehart1=1st-order Hartree part of 2nd-order total energy
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues (hartree)
!!  ek0=0th-order kinetic energy part of 2nd-order total energy.
!!  ek1=1st-order kinetic energy part of 2nd-order total energy.
!!  eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!!  elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!!  enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!!  enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!!  exc1=1st-order exchange-correlation part of 2nd-order total energy
!!  etotal=total energy (sum of 7 contributions) (hartree)
!!  resid(mband*nkpt_rbz*nsppol)=residuals for each band over all k points
!!   of the reduced Brillouin zone, and spins
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!
!! SIDE EFFECTS
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=updated wavefunctions;  if mk1mem<=nkpt_rbz,
!!         these are kept in a disk file, see wffnow and wffnew.
!!  initialized= if 0 the initialization of the RF run is not yet finished
!!  mpi_enreg=informations about MPI parallelization
!!  rhog1(2,nfftf)=array for Fourier transform of RF electron density
!!  rhor1(cplex*nfftf,nspden)=array for RF electron density in electrons/bohr**3.
!!  === if psps%usepaw==1
!!    pawrhoij1(dimpaw1) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!
!! TODO
!! Should be taken away of the list of arguments : ehart, ek, enl, enxc
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      appdig,bec3,die3,ebp3,edie3,etot3,fourdp,getcut,initberry3,int2char4
!!      ioarr,leave_new,leave_test,metric,newfermie1,newvtr3,nselt3,nstdy3
!!      pawmknhat3,qmatrix,rhofermi3,rhotov3,scprqt,status,timab,vtorho3
!!      wffclose,wrtout,xcomm_world,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine scfcv3(atindx,atindx1,blkflg,cg,cgq,cg1,cplex,cprj,cprjq,cpus,dimpaw1,&
&  gh1_rbz,densymop_rf,dimcprj,doccde_rbz,docckqde,dtfil,dtset,&
&  d2bbb,d2lo,d2nl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&  ehart,ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek,ek0,ek1,eloc0,elpsp1,&
&  enl,enl0,enl1,enxc,etotal,exc1,fermie,fform,hdr,idir,indkpt1,&
&  indsy1,initialized,ipert,irrzon1,istwfk_rbz,&
&  kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&
&  mpert,mpi_enreg,mpsang,mpw,mpw1,nattyp,nband_rbz,&
&  nfftf,ngfftf,nhat,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,nspinor,&
&  nsym1,n3xccc,occkq,occ_rbz,&
&  paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawrhoij1,pawtab,&
&  pertcase,phnons1,ph1d,ph1df,&
&  prtbbb,psps,qphon,resid,residm,rhog,rhog1,&
&  rhor,rhor1,rprimd,symaf1,symrc1,symrl1,timrev,&
&  tnons1,usecprj,wffddk,wffnew,wffnow,wfftgs,wfftkq,vpsp1,vtrial,&
&  wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13paw
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_15common
 use interfaces_16response, except_this_one => scfcv3
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!no_abirules
!Needed for integer arrays
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!---
 integer,intent(in) :: cplex,dimpaw1,fform,idir,ipert,mgfftf,mk1mem,mkmem,mkqmem
 integer,intent(in) :: mpert,mpsang,mpw,mpw1,n3xccc,nfftf
 integer,intent(in) :: nkpt,nkpt_rbz,nkxc,nspden
 integer,intent(in) :: nsym1,pertcase,prtbbb,timrev,usecprj
 integer,intent(inout) :: initialized,nspinor
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom),dimcprj(dtset%natom),ngfftf(18)
 integer,intent(out) :: blkflg(3,mpert,3,mpert)
 integer,intent(in) :: indkpt1(nkpt_rbz),indsy1(4,nsym1,dtset%natom)
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,nspden/dtset%nsppol)
 integer,intent(in) :: istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem),nattyp(psps%ntypat)
 integer,intent(in) :: nband_rbz(nkpt_rbz*dtset%nsppol)
 integer,intent(in) :: npwar1(nkpt_rbz),npwarr(nkpt_rbz)
 integer,intent(in) :: symaf1(nsym1),symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 real(dp),intent(in) :: cpus,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,ehart
 real(dp),intent(out) :: eberry,edocc,eeig0,ehart01,ehart1,ek0,ek1,eloc0,elpsp1,enl0
 real(dp),intent(out) :: enl1,exc1,etotal,residm
 real(dp),intent(in) :: eii,ek,enl,enxc
 real(dp),intent(inout) :: fermie
 real(dp),intent(in) :: qphon(3)
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 real(dp),intent(in) :: cg(2,mpw*nspinor*dtset%mband*mkmem*dtset%nsppol)
 real(dp),intent(inout) :: cg1(2,mpw1*nspinor*dtset%mband*mk1mem*dtset%nsppol)
 real(dp),intent(in) :: cgq(2,mpw1*nspinor*dtset%mband*mkqmem*dtset%nsppol)
 real(dp),intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
 real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert)
 real(dp),intent(out) :: d2nl(2,3,mpert,3,mpert)
 real(dp),intent(out) :: gh1_rbz(nkpt_rbz*dtset%ieig2rf,dtset%mband,2,mpw1*nspinor)
 real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(out) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz),kxc(nfftf,nkxc)
 real(dp),intent(in) :: nhat(cplex*nfftf,nspden*psps%usepaw)
 real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom),ph1df(2,3*(2*mgfftf+1)*dtset%natom)
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),nspden/dtset%nsppol)
 real(dp),intent(out) :: resid(dtset%mband*nkpt_rbz*nspden)
 real(dp),intent(in) :: rhog(2,nfftf),rhor(nfftf,nspden),rprimd(3,3)
 real(dp),intent(inout) :: rhog1(2,nfftf),rhor1(cplex*nfftf,nspden),vtrial(nfftf,nspden),xred(3,dtset%natom)
 real(dp),intent(in) :: tnons1(3,nsym1),vpsp1(cplex*nfftf)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xccc3d1(cplex*n3xccc)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm)
 type(cprj_type),intent(in) :: cprj(dimpaw1,nspinor*dtset%mband*mkmem*dtset%nsppol*usecprj)
 type(cprj_type),intent(in) :: cprjq(dtset%natom,nspinor*dtset%mband*mkqmem*dtset%nsppol*usecprj)
 type(datafiles_type),intent(in) :: dtfil
 type(hdr_type),intent(inout) :: hdr
 type(dens_sym_operator_type),intent(in) :: densymop_rf
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(paw_an_type),intent(inout) :: paw_an(dtset%natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(dimpaw1*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wffile_type),intent(inout) :: wffddk,wfftgs,wfftkq
 type(wffile_type),intent(inout) :: wffnow,wffnew

!Local variables-------------------------------
!scalars
 integer,parameter :: level=16,response=1
 integer :: accessfil,afford,choice,coordn,dbl_nnsclo,dielstrt,fformr,fformv
 integer :: fftalg,i1,i2,i3,iatom,ieig2rf,ierr,iexit,ifft,ii,ikpt,index,iprcel
 integer :: ipw,ipw1,ipw2,ir,iscf10_mod,iscf_mod,ispden,ispmix,isppol,istep
 integer :: itypat,jj,lm_size,lmn2_size,me,mgfftdiel,mvdum,n_fftgr,n_index
 integer :: nfftdiel,nfftmix,nfftot,nhat1dim,npawmix,npwdiel,nstep,nsym,optene
 integer :: option,optres,optxc,prtfor,quit,quit_sum,qzero,rdwr,rdwrpaw
 integer :: spaceComm,usexcnhat
 real(dp) :: boxcut,deltae,diecut,diemix,diffor,ecut,ecutf,ecutsus,elast
 real(dp) :: etotal_old,evar,fe1fixed,fermie1,gsqcut,maxfor,res2,rho1_dn
 real(dp) :: rho1im_dn,rho1re_dn,ucvol,vxcavg
 logical :: ex,test_mpi
 character(len=4) :: tag
 character(len=500) :: message
 character(len=fnlen) :: fi1o,filapp,fildata,filfft,filkgs,kgnam
 type(efield_type) :: dtefield
!arrays
 integer :: ngfftmix(18)
 integer,allocatable :: i_rhor(:),i_vresid(:),i_vrespc(:),i_vtrial(:)
 real(dp) :: dielar(7),dummy2(6),favg(3),gmet(3,3),gprimd(3,3),k0(3)
 real(dp) :: kpt_diel(3),rmet(3,3),tollist(12),tsec(2)
 real(dp),allocatable :: dielinv(:,:,:,:,:),f_fftgr(:,:,:),f_paw(:,:)
 real(dp),allocatable :: fcart(:,:),grtn(:,:),nhat1(:,:),nvresid1(:,:)
 real(dp),allocatable :: rhorfermi(:,:),susmat(:,:,:,:,:),vhartr01(:)
 real(dp),allocatable :: vhartr1(:),vt1g(:,:),vtrial1(:,:),vxc1(:,:),work(:)
 type(paw_an_type),allocatable :: paw_an1(:)
 type(paw_ij_type),allocatable :: paw_ij1(:)
 type(pawrhoij_type),allocatable :: rhoij_dum(:)
!no_abirules
 integer,allocatable :: pwindall(:,:,:)
 real(dp),allocatable ::qmat(:,:,:,:,:,:)

! *********************************************************************

!DEBUG
!write(6,*)' scfcv3 : enter'
!ENDDEBUG

 call timab(120,1,tsec)
 call timab(154,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if dtset%prtvol==-level
 if(dtset%prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' scfcv3: enter '
  call wrtout(06,message,'COLL')
 end if

!Init me
 call xme_init(mpi_enreg,me)

!
!If dtset%accesswff == 2 set all array outputs to netcdf format
!
 accessfil = 0
 if (dtset%accesswff == 2) then
  accessfil = 1
 end if
 if (dtset%accesswff == 3) then
  accessfil = 3
 end if

!Save some variables from dataset definition
 ecut=dtset%ecut
 ecutf=ecut;if (psps%usepaw==1.and.pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg
 iprcel=dtset%iprcel
 tollist(1)=dtset%tolmxf;tollist(2)=dtset%tolwfr
 tollist(3)=dtset%toldff;tollist(4)=dtset%toldfe
 tollist(6)=dtset%tolvrs;tollist(7)=dtset%tolrff
 nstep=dtset%nstep
 iscf_mod=dtset%iscf
 iscf10_mod=mod(iscf_mod,10)
 qzero=0;if(qphon(1)**2+qphon(2)**2+qphon(3)**2 < tol14) qzero=1

!The value of iscf must be modified if ddk perturbation, (maybe eventually natom+5 too) see loper3.f
 if (ipert==dtset%natom+1) iscf_mod=-3

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Some variables need to be initialized/nullify at start
 quit=0 ; dbl_nnsclo=0 ; elast=zero
!This might be taken away later
 edocc=zero ; eeig0=zero ; ehart01=zero ; ehart1=zero ; ek0=zero ; ek1=zero
 eloc0=zero ; elpsp1=zero ; enl0=zero ; enl1=zero ; exc1=zero
 deltae=zero ; fermie1=zero
 optres=merge(0,1,abs(iscf_mod)<10)
 usexcnhat=0

!Prepare the name of the _FFT file and the _KGS file
 filfft=trim(dtfil%filnam_ds(5))//'_FFT'
 filkgs=trim(dtfil%filnam_ds(5))//'_KGS'
 filapp=trim(dtfil%filnam_ds(5))
 if(mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then
! BEGIN TF_CHANGES
  call int2char4(me,tag)
! END TF_CHANGES
  filfft=trim(filfft)//'_P-'//tag
  filkgs=trim(filkgs)//'_P-'//tag
 end if

!Examine tolerance criteria, and eventually  print a line to the output
!file (with choice=1, the only non-dummy arguments of scprqt are
!nstep, tollist and iscf - still, diffor,res2,prtfor,fcart are here initialized to 0)
 choice=1 ; prtfor=0 ; diffor=zero ; res2=zero
 allocate(fcart(3,dtset%natom))
 call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
& etotal,favg,fcart,fermie,filapp,dtfil%filnam_ds(1),&
& 1,iscf_mod,istep,kpt_rbz,maxfor,&
& mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
& nstep,occ_rbz,0,prtfor,&
& quit,res2,resid,residm,response,&
& tollist,psps%usepaw,vxcavg,wtk_rbz,xred)

!Various allocations (potentials, gradients, ...)
 call status(0,dtfil%filstat,iexit,level,'allocate      ')
 allocate(vhartr1(cplex*nfftf),vtrial1(cplex*nfftf,nspden))

!Allocations/initializations for PAW only
 if(psps%usepaw==1) then
  usexcnhat=maxval(pawtab(:)%vlocopt)
! 1st-order compensation density
  nhat1dim=nfftf;if (ipert<=dtset%natom) nhat1dim=pawfgrtab(ipert)%nfgd
  allocate(nhat1(cplex*nhat1dim,dtset%nspden));if (nstep==0) nhat1=zero
! Variables, 1st-order arrays related to the PAW spheres
  allocate(paw_ij1(dimpaw1),paw_an1(dimpaw1))
  do iatom=1,dimpaw1
   itypat=dtset%typat(iatom)
   lmn2_size=pawtab(itypat)%lmn2_size
   lm_size=min(dtset%pawlcutd,pawtab(itypat)%l_size)**2
   paw_an1(iatom)%cplex    =cplex
   paw_an1(iatom)%angl_size=pawang%angl_size
   paw_an1(iatom)%mesh_size=pawtab(itypat)%mesh_size
   paw_an1(iatom)%nspden   =dtset%nspden
   paw_an1(iatom)%lm_size  =lm_size
   allocate(paw_an1(iatom)%lmselect(lm_size))
   paw_ij1(iatom)%cplex    =cplex
   paw_ij1(iatom)%cplex_dij=max(cplex,nspinor)
   paw_ij1(iatom)%nspden   =dtset%nspden
   paw_ij1(iatom)%nsppol   =dtset%nsppol
   paw_ij1(iatom)%lmn_size =pawtab(itypat)%lmn_size
   paw_ij1(iatom)%lmn2_size=lmn2_size
   paw_ij1(iatom)%has_dijhat=1
   paw_ij1(iatom)%has_dijxc=0
   nullify(paw_ij1(iatom)%dijxc)
   paw_ij1(iatom)%has_dijxc_val=0
   nullify(paw_ij1(iatom)%dijxc_val)
   paw_ij1(iatom)%has_dijso=0
   nullify(paw_ij1(iatom)%dijso)
   paw_ij1(iatom)%has_dijU  =0
   nullify(paw_ij1(iatom)%dijU)
   allocate(paw_ij1(iatom)%dij(cplex*lmn2_size,dtset%nspden))
   allocate(paw_ij1(iatom)%dijhat(cplex*lmn2_size,dtset%nspden))
   pawrhoij1(iatom)%use_rhoijres=1
   allocate(pawrhoij1(iatom)%rhoijres(cplex*lmn2_size,dtset%nspden))
   do ispden=1,dtset%nspden
    pawrhoij1(iatom)%rhoijres(:,ispden)=zero
   end do
  end do
 end if ! PAW

!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (iscf_mod>0.or.iscf_mod==-3) then
  allocate(nvresid1(cplex*nfftf,dtset%nspden))
  if (nstep==0) nvresid1=zero
 end if
 if(nstep>0 .and. iscf_mod>0) then
  dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
  dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
  dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
  if(iscf10_mod==1) then
!  For iscf_mod==1, five additional vectors are needed
!  The index 1 is attributed to the old trial potential,
!  The new residual potential, and the new
!  preconditioned residual potential receive now a temporary index
!  The indices number 4 and 5 are attributed to work vectors.
   n_fftgr=5 ; n_index=1
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3
  else if(iscf10_mod==2) then
!  For iscf==2, three additional vectors are needed.
!  The index number 1 is attributed to the old trial vector
!  The new residual potential, and the new preconditioned
!  residual potential, receive now a temporary index.
   n_fftgr=3 ; n_index=1
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3
  else if(iscf10_mod==3) then
!  For iscf==3 , four additional vectors are needed.
!  The index number 1 is attributed to the old trial vector
!  The new residual potential, and the new and old preconditioned
!  residual potential, receive now a temporary index.
   n_fftgr=4 ; n_index=2
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3 ; i_vrespc(2)=4
  else if (iscf10_mod==4) then
!  For iscf==4 , six additional vectors are needed.
!  The indices number 1 and 2 are attributed to two old trial vectors
!  The new residual potential, and the new and two old preconditioned
!  residual potentials, receive now a temporary index.
   n_fftgr=6 ; n_index=3
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vtrial(2)=2 ; i_vresid(1)=3
   i_vrespc(1)=4 ; i_vrespc(2)=5 ; i_vrespc(3)=6
  else if((iscf10_mod==5).or.(iscf10_mod==6)) then
!  For iscf==5 or 6, ten additional vectors are needed
!  The index number 1 is attributed to the old trial vector
!  The index number 6 is attributed to the search vector
!  Other indices are attributed now. Altogether ten vectors
   n_fftgr=10 ; n_index=3
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3 ; i_vresid(2)=4 ; i_vrespc(2)=5
   i_vresid(3)=7 ; i_vrespc(3)=8 ; i_rhor(2)=9 ; i_rhor(3)=10
  else if(iscf10_mod==7) then
!  For iscf==7, lot of additional vectors are needed
!  The index number 1 is attributed to the old trial vector
!  The index number 2 is attributed to the old residual
!  The indices number 2 and 3 are attributed to two old precond. residuals
!  Other indices are attributed now.
   n_fftgr=2+2*dtset%npulayit ; n_index=1+dtset%npulayit
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vrespc(dtset%npulayit+1)=2*dtset%npulayit+1; i_vresid(1)=2*dtset%npulayit+2
   do ii=1,dtset%npulayit
    i_vtrial(ii)=2*ii-1 ; i_vrespc(ii)=2*ii
   end do
  end if ! iscf cases
! Additional allocation for mixing within PAW
  npawmix=0
  if(psps%usepaw==1) then
   do iatom=1,dimpaw1
    itypat=dtset%typat(iatom)
    allocate(pawrhoij1(iatom)%kpawmix(pawtab(itypat)%lmnmix_sz))
    pawrhoij1(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
    pawrhoij1(iatom)%kpawmix=pawtab(itypat)%kmix
    npawmix=npawmix+dtset%nspden*pawtab(itypat)%lmnmix_sz
   end do
  end if

! The big array containing functions defined on the fft grid is allocated now.
! Note however, that a zero value of mffmem will cause allocation with
! the third dimension set to zero ! In this case, another temporary
! will be used inside newvtr3.
  if (psps%usepaw==0) npawmix=0
  if (psps%usepaw==1) allocate(f_paw(cplex*npawmix,n_fftgr*dtset%mffmem))
  if (psps%usepaw==1.and.dtset%pawmixdg==0) then
   ispmix=2;nfftmix=dtset%nfft;ngfftmix(:)=dtset%ngfft(:)
  else
   ispmix=1;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
  end if
  allocate(f_fftgr(ispmix*cplex*nfftmix,dtset%nspden,n_fftgr*dtset%mffmem))
 end if ! iscf, nstep

!Here, allocate arrays for computation of susceptibility and dielectric matrix or for TDDFT
 if( (nstep>0 .and. iscf_mod>0) .or. iscf_mod==-1 ) then
! Here, for TDDFT, artificially set iprcel . Also set a variable to reduce
! the memory needs.
  afford=1
  if(iscf_mod==-1) then
   iprcel=21
   afford=0
  end if
  npwdiel=1
  mgfftdiel=1
  nfftdiel=1
! Now, performs allocation
! CAUTION : the dimensions are still those of GS, except for phnonsdiel
  allocate(dielinv(2,npwdiel*afford,nspden,npwdiel,nspden))
  allocate(susmat(2,npwdiel*afford,nspden,npwdiel,nspden))
 end if

!Initialize Berry-phase related stuffs
 if (dtset%berryopt == 4) then
  allocate(pwindall(max(mpw,mpw1)*mkmem,8,3))
  call initberry3(dtefield,dtfil,dtset,gmet,kg,kg1,dtset%mband,mkmem,mpi_enreg,&
&  mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,occ_rbz,pwindall,rprimd)
! calculate inverse of the overlap matrix
  allocate(qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3))
  call qmatrix(cg,dtefield,qmat,mpw,mpw1,mkmem,dtset%mband,npwarr,nkpt,nspinor,dtset%nsppol,pwindall)
 end if

!Compute large sphere cut-off gsqcut
 if (psps%usepaw==1) then
  write(message,'(2a)') ch10,' FFT (fine) grid used for densities/potentials:'
  call wrtout(6,message,'COLL')
 end if
 k0(:)=zero
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,6,k0,ngfftf)

 call timab(154,2,tsec)

!######################################################################
!PERFORM ELECTRONIC ITERATIONS
!######################################################################

!Offer option of computing 2nd-order total energy with existing
!wavefunctions when nstep<=0, else do nstep iterations
!Note that for non-self-consistent calculations, this loop will be exited
!after the first call to vtorho3

!Pass through the first routines even when nstep==0
 write(*,*) 'scfcv3, nstep=', max(1,nstep)
 do istep=1,max(1,nstep)

! ######################################################################
! The following steps are done once
! ----------------------------------------------------------------------
  if (istep==1)then

!  PAW only: we sometimes have to compute 1st-order compensation density
!  and eventually add it to density from 1st-order WFs
!  ----------------------------------------------------------------------
   if (psps%usepaw==1.and.((usexcnhat==0) &
&   .or.(dtfil%ireadwf/=0.and.dtset%get1den==0.and.initialized==0))) then
    call timab(564,1,tsec)
    call pawmknhat3(cplex,idir,ipert,0,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat1dim,&
&    dimpaw1,nspden,psps%ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhat1,pawrhoij,&
&    pawrhoij1,pawtab,dtset%typat)
    if (dtfil%ireadwf/=0.and.dtset%get1den==0.and.initialized==0) then
     if (dimpaw1==nfftf) then
      rhor1(:,:)=rhor1(:,:)+nhat1(:,:)
     else
      if (cplex==1) then
       do ispden=1,nspden
        do ii=1,pawfgrtab(ipert)%nfgd
         jj=pawfgrtab(ipert)%ifftsph(ii)
         rhor1(jj,ispden)=rhor1(jj,ispden)+nhat1(ii,ispden)
        end do
       end do
      else
       do ispden=1,nspden
        do ii=1,pawfgrtab(ipert)%nfgd
         jj=2*pawfgrtab(ipert)%ifftsph(ii)
         rhor1(jj-1:jj,ispden)=rhor1(jj-1:jj,ispden)+nhat1(2*ii-1:2*ii,ispden)
        end do
       end do
      end if
     end if
     call fourdp(cplex,rhog1,rhor1(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
    end if
    call timab(564,2,tsec)
   end if

!  Set initial guess for 1st-order potential
!  ----------------------------------------------------------------------
   call status(istep,dtfil%filstat,iexit,level,'get vtrial1   ')
   optene=-1;option=1
   call rhotov3(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&
&   kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat,nhat1,nhat1dim,nkxc,nspden,n3xccc,&
&   optene,option,dtset%paral_kgb,pawfgrtab,dtset%qptn,rhog,rhog1,rhor,rhor1,&
&   rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,nvresid1,res2,vtrial1,xccc3d1)

!  For Q=0 and metallic occupation, initialize quantities needed to
!  compute the first-order Fermi energy
!  ----------------------------------------------------------------------
   if(qzero==1 .and. (dtset%occopt>=3 .and. dtset%occopt <=7) .and.&
&   (ipert<=dtset%natom .or. ipert==dtset%natom+3 .or. ipert==dtset%natom+4) .and. dtset%frzfermi==0) then
    allocate(rhorfermi(cplex*nfftf,nspden))
    call rhofermi3(atindx,atindx1,cg,cgq,cplex,densymop_rf,&
&    doccde_rbz,docckqde,dtfil,dtset,&
&    edocc,eeig0,eigenq,eigen0,eigen1,&
&    fe1fixed,gmet,gprimd,gsqcut,hdr,idir,&
&    ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,mgfftf,&
&    mkmem,mkqmem,mk1mem,mpi_enreg,mpsang,mpw,mpw1,&
&    dtset%natom,nattyp,nband_rbz,nfftf,nkpt_rbz,npwarr,npwar1,nspden,nspinor,&
&    dtset%nsppol,nsym1,dtset%ntypat,occkq,occ_rbz,phnons1,&
&    ph1d,dtset%prtvol,psps,rhorfermi,rmet,symaf1,tnons1,ucvol,&
&    wfftgs,wfftkq,wtk_rbz,xred,ylm,ylm1,ylmgr1)
   end if !First-order fermi energy setup

!  End the condition of istep==1
  end if

! ######################################################################
! The following steps are done at every iteration
! ----------------------------------------------------------------------

  if (psps%usepaw==1)then
!  PAW: INSERT HERE CALLS TO PAWDENPOT3, PAWDIJ3, SYMDIJ3 (to be copied from scfcv)
  end if

! No need to continue and call vtorho3, when nstep==0
  if(nstep==0)exit

! ######################################################################
! The following steps are done only when nstep>0
! ----------------------------------------------------------------------
  call status(istep,dtfil%filstat,iexit,level,'loop istep    ')

  if(iscf_mod>0)then
   write(message, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
   call wrtout(06,message,'COLL')
  end if

! For Q=0 and metallic occupation, calculate the first-order Fermi energy
  if(qzero==1 .and. (dtset%occopt>=3 .and. dtset%occopt <=7) .and.&
&  (ipert<=dtset%natom .or. ipert==dtset%natom+3 .or. ipert==dtset%natom+4) .and. dtset%frzfermi==0) then
   nfftot=ngfftf(1)*ngfftf(2)*ngfftf(3)
   call newfermie1(cplex,fermie1,fe1fixed,istep,&
&   mpi_enreg,nfftf,nfftot,nspden,dtset%occopt,&
&   dtset%prtvol,rhorfermi,ucvol,vtrial1)
  end if

! DEBUG
! write(6,*)' scfcv3 : before vtorho3, vtrial1(1,1)=',vtrial1(1,1)
! ENDDEBUG

! ######################################################################
! Compute the 1st-order density rho1 from the 1st-order trial potential
! ----------------------------------------------------------------------
  call status(istep,dtfil%filstat,iexit,level,'call vtorho3  ')
  call vtorho3(atindx,atindx1,cg,cgq,cg1,cplex,cprj,cprjq,cpus,dbl_nnsclo,gh1_rbz,densymop_rf,&
&  dimcprj,dimpaw1,doccde_rbz,docckqde,dtefield,&
&  dtfil,dtset,edocc,eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,&
&  enl0,enl1,fermie1,gmet,gprimd,gsqcut,hdr,idir,indsy1,&
&  ipert,irrzon1,istep,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,&
&  mkmem,mkqmem,mk1mem,mpi_enreg,mpsang,mpw,mpw1,&
&  dtset%natom,nattyp,nband_rbz,nfftf,nhat1,nhat1dim,nkpt_rbz,npwarr,npwar1,res2,nspden,nspinor,&
&  dtset%nsppol,nsym1,dtset%ntypat,nvresid1,occkq,occ_rbz,optres,&
&  paw_ij,paw_ij1,pawang,pawfgr,pawfgrtab,pawrhoij,pawrhoij1,pawtab,&
&  phnons1,ph1d,dtset%prtvol,psps,pwindall,qmat,resid,residm,rhog1,rhor1,rmet,&
&  rprimd,symaf1,symrc1,tnons1,ucvol,usecprj,&
&  wffddk,wffnew,wffnow,wfftgs,wfftkq,vtrial,vtrial1,wtk_rbz,xred,ylm,ylm1,ylmgr1)
  call status(istep,dtfil%filstat,iexit,level,'after vtorho3 ')

  if (dtset%berryopt == 4) then
!  calculate \Omega E \cdot P term
   if (ipert<=dtset%natom) then
!   phonon perturbation
    call  ebp3(cg,cg1,dtefield,eberry,dtset%mband,mkmem,&
&    mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,nspinor,pwindall,qmat)
   else if (ipert==dtset%natom+2) then
!   electric field perturbation
    call  edie3(cg,cg1,dtefield,eberry,idir,dtset%mband,mkmem,&
&    mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,nspinor,pwindall,qmat,rprimd)
   end if
  end if

! ######################################################################
! Skip out of step loop if non-SCF (completed)
! ----------------------------------------------------------------------

! Indeed, nstep loops have been done inside vtorho3
  if (iscf_mod<=0 .and. iscf_mod/=-3) exit

! ######################################################################
! In case of density mixing , compute the total 2nd-order energy,
! check the exit criterion,
! then mix the 1st-order density
! ----------------------------------------------------------------------

  if (iscf_mod>=10) then
   optene = 1 ! use double counting scheme
   call etot3(dtset%berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,&
&   efrloc,efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
&   enl0,enl1,etotal,evar,exc1,ipert,dtset%natom,optene)

   call timab(152,1,tsec)
   choice=2
   call status(istep,dtfil%filstat,iexit,level,'print info    ')
   call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
&   etotal,favg,fcart,fermie,filapp,dtfil%filnam_ds(1),&
&   1,iscf_mod,istep,kpt_rbz,maxfor,&
&   mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
&   nstep,occ_rbz,0,prtfor,&
&   quit,res2,resid,residm,response,&
&   tollist,psps%usepaw,vxcavg,wtk_rbz,xred)
   call timab(152,2,tsec)
   if (istep==nstep) quit=1
!  If criteria in scprqt say to quit, then exit the loop over istep
   quit_sum=quit
   if (mpi_enreg%paral_compil_kpt==1) then
    call xcomm_world(mpi_enreg,spaceComm)
    call xsum_mpi(quit_sum,spaceComm,ierr)
   end if
   if (quit_sum>0) exit
   call status(istep,dtfil%filstat,iexit,level,'call newrho   ')
!  INSERT HERE CALL TO NEWRHO3 : to be implemented
   initialized=1
  end if

! ######################################################################
! Compute the new 1st-order potential from the 1st-order density
! ----------------------------------------------------------------------

  call status(istep,dtfil%filstat,iexit,level,'call rhotov3   ')
  optene=2*optres; !if(psps%usepaw==1) optene=2
  call rhotov3(cplex,ehart01,ehart1,elpsp1,exc1,gmet,gprimd,gsqcut,idir,ipert,&
&  kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat,nhat1,nhat1dim,nkxc,nspden,n3xccc,&
&  optene,optres,dtset%paral_kgb,pawfgrtab,dtset%qptn,rhog,rhog1,rhor,rhor1,&
&  rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,nvresid1,res2,vtrial1,xccc3d1)

! ######################################################################
! In case of potential mixing , compute the total 2nd-order energy,
! check the exit criterion,
! then mix the 1st-order potential
! ----------------------------------------------------------------------

  if (iscf_mod<10) then
   optene = 0 ! use direct scheme
   call etot3(dtset%berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,&
&   efrloc,efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
&   enl0,enl1,etotal,evar,exc1,ipert,dtset%natom,optene)

   call timab(152,1,tsec)
   choice=2
   call status(istep,dtfil%filstat,iexit,level,'print info    ')
   call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
&   etotal,favg,fcart,fermie,filapp,dtfil%filnam_ds(1),&
&   1,iscf_mod,istep,kpt_rbz,maxfor,&
&   mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
&   nstep,occ_rbz,0,prtfor,&
&   quit,res2,resid,residm,response,&
&   tollist,psps%usepaw,vxcavg,wtk_rbz,xred)
   call timab(152,2,tsec)
!  If criteria in scprqt say to quit, then exit the loop over istep
   quit_sum=quit
   if (mpi_enreg%paral_compil_kpt==1) then
    call xcomm_world(mpi_enreg,spaceComm)
    call xsum_mpi(quit_sum,spaceComm,ierr)
   end if
   if (quit_sum>0) exit
   if(iscf_mod/=-3)then
!   Note that nvresid1 and vtrial1 are called vresid and vtrial inside this routine
    if (psps%usepaw==1) stop " newvtr3 has to be adapted to pawrhoij1 !"
    call newvtr3(cplex,dbl_nnsclo,dielar,etotal,filfft,f_fftgr,&
&    initialized,iscf_mod,dtset%isecur,istep,i_rhor,&
&    i_vresid,i_vrespc,i_vtrial,dtset%mffmem,mpi_enreg,dtset%natom,nfftf,ngfftf,nspden,&
&    n_fftgr,n_index,dtset%paral_kgb,pawrhoij1,dtset%qptn,rhor1,rprimd,nvresid1,vtrial1,xred)
    initialized=1
   end if
  end if

! ######################################################################
! END MINIMIZATION ITERATIONS
! ######################################################################

! Note that there are different "exit" instructions within the loop
 end do ! istep

 if (iscf_mod>0.or.iscf_mod==-3) deallocate(nvresid1)

!DEBUG
!write(6,*)' scfcv3 : after istep loop, continue'
!stop
!ENDDEBUG

!if(nstep==0) then (not implemented)
!end if


!######################################################################
!Additional steps after SC iterations
!----------------------------------------------------------------------

 call timab(160,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'endloop istep ')

!Delete eventual _FFT file
 if(dtset%mffmem==0)then
  inquire (file=filfft,exist=ex)
  if(ex)then
   open(unit=tmp_unit,file=filfft,form='unformatted',status='old')
   close(unit=tmp_unit,status='DELETE')
  end if
 end if

!Eventually close the dot file, before calling nstdy3
 if(ipert==dtset%natom+2 .and. sum( (dtset%qptn(1:3)) **2 ) <= 1.0d-7 .and. (dtset%berryopt .ne. 4) )then
  call WffClose(wffddk,ierr)
 end if

 call timab(160,2,tsec)
 call timab(147,1,tsec)

 if(ipert==dtset%natom+3 .or. ipert==dtset%natom+4) then
  call status(0,dtfil%filstat,iexit,level,'enter nselt3  ')
  call nselt3(atindx,atindx1,blkflg,cg,cg1,cplex,doccde_rbz,docckqde,&
&  d2bbb,d2lo,d2nl,ecut,dtset%ecutsm,dtset%effmass,eigen0,eigen1,fform,&
&  gmet,gprimd,gsqcut,idir,indkpt1,indsy1,&
&  ipert,iscf_mod,istep,istwfk_rbz,kg,kg1,kpt_rbz,kxc,dtset%mband,mgfftf,&
&  mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
&  dtset%natom,nattyp,dtset%nband,nband_rbz,nfftf,ngfftf,&
&  nkpt,nkpt_rbz,nkxc,dtset%nline,dtset%nloalg,&
&  npwarr,npwar1,nspden,nspinor,dtset%nsppol,&
&  nsym1,dtset%ntypat,occkq,dtset%occopt,occ_rbz,dtset%ortalg,&
&  dtset%paral_kgb,ph1d,dtset%prtbbb,dtset%prtvol,psps,dtset%qptn,rhog,rhog1,&
&  rhor,rhor1,rmet,rprimd,symrc1,dtset%tsmear,dtset%typat,ucvol,&
&  dtfil%unkg,dtfil%unkg1,&
&  wffnow,wfftgs,dtfil%unylm,dtfil%unylm1,&
&  dtfil%fnamewffddk,wtk_rbz,xred,ylm,ylm1,ylmgr,ylmgr1)
 end if
 if(ipert<=dtset%natom+4)then
  call status(0,dtfil%filstat,iexit,level,'enter nstdy3  ')
  call nstdy3(atindx,atindx1,blkflg,cg,cg1,cplex,doccde_rbz,docckqde,&
&  dtset,d2bbb,d2lo,d2nl,ecut,dtset%ecutsm,eigen0,eigen1,fform,&
&  gmet,gprimd,gsqcut,idir,indkpt1,indsy1,&
&  ipert,iscf_mod,istep,istwfk_rbz,kg,kg1,kpt_rbz,kxc,dtset%mband,mgfftf,&
&  mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
&  dtset%natom,nattyp,dtset%nband,nband_rbz,nfftf,ngfftf,&
&  nkpt,nkpt_rbz,nkxc,dtset%nline,dtset%nloalg,&
&  npwarr,npwar1,nspden,nspinor,dtset%nsppol,&
&  nsym1,dtset%ntypat,occkq,dtset%occopt,occ_rbz,dtset%ortalg,&
&  dtset%paral_kgb,ph1d,dtset%prtbbb,dtset%prtvol,psps,dtset%qptn,rhog1,&
&  rhor1,rmet,rprimd,symrc1,dtset%tsmear,dtset%typat,ucvol,&
&  dtfil%unkg,dtfil%unkg1,&
&  wffnow,wfftgs,dtfil%unylm,dtfil%unylm1,&
&  dtfil%fnamewffddk,wtk_rbz,xred,ylm,ylm1)
 end if

!if(ipert==dtset%natom+5)then
!!calculate the non-stationary expression for the
!!second derivative of the total energy
!call nstdy3(atindx,atindx1,blkflg,cg,cg1,cplex,doccde_rbz,docckqde,&
!&  dtset,d2bbb,d2lo,d2nl,ecut,dtset%ecutsm,eigen0,eigen1,fform,&
!&  gmet,gprimd,gsqcut,idir,indkpt1,indsy1,&
!&  dtset%natom+2,iscf_mod,istep,istwfk_rbz,kg,kg1,kpt_rbz,kxc,dtset%mband,mgfftf,&
!&  mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
!&  dtset%natom,nattyp,dtset%nband,nband_rbz,nfftf,ngfftf,&
!&  nkpt,nkpt_rbz,nkxc,dtset%nline,dtset%nloalg,&
!&  npwarr,npwar1,nspden,nspinor,dtset%nsppol,&
!&  nsym1,dtset%ntypat,occkq,dtset%occopt,occ_rbz,dtset%ortalg,&
!&  dtset%paral_kgb,ph1d,dtset%prtbbb,dtset%prtvol,psps,dtset%qptn,rhog1,&
!&  rhor1,rmet,rprimd,symrc1,dtset%tsmear,dtset%typat,ucvol,&
!&  dtfil%unkg,dtfil%unkg1,&
!&  wffnow,wfftgs,dtfil%unylm,dtfil%unylm1,&
!&  dtfil%fnamewffddk,wtk_rbz,xred,ylm,ylm1)
!end if

 call timab(147,2,tsec)
 call timab(160,1,tsec)

!DEBUG
!write(6,*)' scfcv3 : will call bec3 and die3 '
!ENDDEBUG

!calculate Born effective charge and store it in d2lo
 if (dtset%berryopt == 4 .and. ipert <= dtset%natom) then

  call  bec3(cg,cg1,dtefield,dtset%natom,d2lo,idir,ipert,dtset%mband,mkmem,&
&  mpw,mpw1,mpert,nkpt,npwarr,npwar1,dtset%nsppol,nspinor,pwindall,qmat,rprimd)
  blkflg(:,dtset%natom+2,:,1:dtset%natom)=1

! DEBUG
! write(6,*)' scfcv3 : blkflg for Born effective charge has been set to 1 '
! ENDDEBUG

 end if

!calculate dielectric tensor and store it in d2lo
 if (dtset%berryopt == 4 .and. ipert == dtset%natom+2 ) then

  call die3(cg,cg1,dtefield,d2lo,idir,ipert,dtset%mband,mkmem,&
&  mpw,mpw1,mpert,nkpt,npwarr,npwar1,dtset%nsppol,nspinor,pwindall,qmat,rprimd)
  blkflg(:,dtset%natom+2,:,dtset%natom+2)=1

! DEBUG
! write(6,*)' scfcv3 : blkflg for dielectric tensor has been set to 1 '
! ENDDEBUG

 end if

!If SCF convergence was not reached (for nstep>0),
!print a warning to the output file (non-dummy arguments: nstep,
!residm, diffor - infos from tollist have been saved inside )
 call status(0,dtfil%filstat,iexit,level,'call scprqt(en')
 choice=3
 call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
& etotal,favg,fcart,fermie,filapp,dtfil%filnam_ds(1),&
& 1,iscf_mod,istep,kpt_rbz,maxfor,&
& mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
& nstep,occ_rbz,0,prtfor,&
& quit,res2,resid,residm,response,&
& tollist,psps%usepaw,vxcavg,wtk_rbz,xred)

!Optionally provide output of charge density and/or potential in real space,
!as well as analysis of geometrical factors (bond lengths and bond angles).
!Warnings :
!- core charge is excluded from the charge density;
!- the potential is the INPUT vtrial.
 if(me==0) then
  if (dtset%prtden>0) then
   call status(0,dtfil%filstat,iexit,level,'call ioarr-den')
   rdwr=2 ; fformr=52 ; rdwrpaw=0
   fildata=trim(dtfil%filnam_ds(4))//'_DEN'
   call appdig(pertcase,fildata,fi1o)
   call ioarr(accessfil,rhor1, dtset, etotal,fformr,fi1o,hdr, mpi_enreg, &
&   cplex*nfftf,rhoij_dum,rdwr,rdwrpaw,ngfftf)
  end if
  if (dtset%prtpot>0) then
   call status(0,dtfil%filstat,iexit,level,'call ioarr-pot')
   rdwr=2 ; fformv=102 ; rdwrpaw=0
   fildata=trim(dtfil%filnam_ds(4))//'_POT'
   call appdig(pertcase,fildata,fi1o)
   call ioarr(accessfil,vtrial1, dtset, fermie,fformv,fi1o,hdr, mpi_enreg, &
&   cplex*nfftf,rhoij_dum,rdwr,rdwrpaw,ngfftf)
  end if
 end if
 if(mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then
  call timab(61,1,tsec)
! BEGIN TF_CHANGES
  call leave_test(mpi_enreg)
! END TF_CHANGES
  call timab(61,2,tsec)
 end if

!Debugging : print the different parts of rhor1
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(dtset%prtvol==-level)then
  write(message,'(a)') '   ir       rhor(ir)     '
  call wrtout(06,message,'COLL')
  do ir=1,nfftf
   if(ir<=11 .or. mod(ir,301)==0 )then
    write(message,'(i5,a,es13.6)')ir,' ',rhor(ir,1)
    call wrtout(06,message,'COLL')
    if(nspden==2)then
     write(message,'(a,es13.6)')'      ',rhor(ir,2)
     call wrtout(06,message,'COLL')
    end if
   end if
  end do
 end if

!Deallocate the arrays
 deallocate(fcart,vtrial1,vhartr1)
 if(nstep>0 .and. iscf_mod>0) then
  deallocate(i_rhor,i_vtrial,i_vresid,i_vrespc)
  if (psps%usepaw==1) deallocate(f_paw)
  deallocate(f_fftgr)
 end if
 if( (nstep>0 .and. iscf_mod>0) .or. iscf_mod==-1 ) then
  deallocate(dielinv)
  deallocate(susmat) !! added by MM
 end if
 if (dtset%berryopt == 4) then
  deallocate(pwindall,qmat)
  deallocate(dtefield%ikpt_dk,dtefield%cgindex,dtefield%idxkstr,dtefield%kgindex)
 end if
 if(allocated(rhorfermi)) deallocate(rhorfermi)
 if(psps%usepaw==1) then
  deallocate(nhat1)
  do iatom=1,dimpaw1
   itypat=dtset%typat(iatom)
   deallocate(paw_an1(iatom)%lmselect)
   deallocate(paw_ij1(iatom)%dij)
   if (paw_ij1(iatom)%has_dijhat>0) deallocate(paw_ij1(iatom)%dijhat)
   if (paw_ij1(iatom)%has_dijxc>0) deallocate(paw_ij1(iatom)%dijxc)
   if (paw_ij1(iatom)%has_dijxc_val>0) deallocate(paw_ij1(iatom)%dijxc_val)
   if (paw_ij1(iatom)%has_dijso>0) deallocate(paw_ij1(iatom)%dijso)
   if (paw_ij1(iatom)%has_dijU>0) deallocate(paw_ij1(iatom)%dijU)
   deallocate(pawrhoij1(iatom)%rhoijres);pawrhoij1(iatom)%use_rhoijres=0
   if (iscf_mod>0) then
    pawrhoij1(iatom)%lmnmix_sz=0
    deallocate(pawrhoij1(iatom)%kpawmix)
   end if
  end do
  deallocate(paw_an1,paw_ij1)
 end if ! PAW

!Structured debugging : if dtset%prtvol=-level, stop here.
 if(dtset%prtvol==-level)then
  write(message,'(a1,a,a1,a,i2,a)') ch10,' scfcv3: exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(160,2,tsec)
 call timab(120,2,tsec)

!DEBUG
!write(6,*)' scfcv3: exit '
!ENDDEBUG

end subroutine scfcv3
!!***
