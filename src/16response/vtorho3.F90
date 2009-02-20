!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtorho3
!! NAME
!! vtorho3
!!
!! FUNCTION
!! This routine compute the new 1-density from a fixed 1-potential (vtrial1)
!! but might also simply compute eigenvectors and eigenvalues.
!! The main part of it is a wf update over all k points
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, AR, DRH, MB, XW, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS
!!    wavefunctions at k+q.
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of
!!    RF wavefunctions at k,q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  cprj(dimpaw1,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  cpus= cpu time limit in seconds
!!  dbl_nnsclo=if 1, will double the value of dtset%nnsclo
!!  densymop_rf <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (response-function)
!!  dimcprj(natom)=array of dimensions of arrays cprj, cprjq
!!  dimpaw1= -PAW only- dimension of 1st-order on-site qunatities (rhoij, Dij...)
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  dtefield = variables related to response Berry-phase calculation
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  fermie1=derivative of fermi energy wrt (strain) perturbation
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  idir=direction of the perturbation
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,nspden/nsppol)=irreducible zone data
!!  istep=index of the number of steps in the routine scfcv
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  mband=maximum number of bands
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mkqmem =number of k+q points which can fit in memory (GS data); 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs - see comment in respfn.F90)
!!  nhat1dim=1st dimension of array nhat1
!!  nkpt_rbz=number of k points in the IBZ for this perturbation
!!  mpi_enreg=informations about MPI parallelization
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!   perturbation
!!  ntypat=number of types of atoms in unit cell.
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band and k
!!   (usually 2)
!!  optres=0: the new value of the density is computed in place of the input value
!!         1: only the density residual is computed ; the input density is kept
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  paw_ij1(dimpaw1) <type(paw_ij_type)>=1st-order paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  phnons1(2,dtset%nfft**(1-1/nsym1),nspden/nsppol)=nonsymmorphic translation phases
!!  ph1d(2,3*(2*dtset%mgfft+1)*natom)=one-dimensional structure factor information
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pwindall(max(mpw,mpw1)*mkmem,8,3) = array used to compute the overlap matrices
!!  qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt,2,3) =
!!  inverse of the overlap matrix
!!  rmet(3,3)=real space metric (bohr**2)
!!  symaf1(nsym1)=(anti)ferromagnetic part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  tnons1(3,nsym1)=reduced nonsymmorphic translations
!!   (tnons1 are in terms of real space primitive translations)
!!  ucvol=unit cell volume in bohr**3.
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  wffddk=struct info for wf dot (ddk) file.
!!  wffnew,wffnow=struct info for RF wf disk files.
!!  wfftgs,wfftkq=struct info for GS wf disk files.
!!  vtrial(nfftf,nspden)=GS Vtrial(r).
!!  vtrial1(cplex*nfftf,nspden)=INPUT RF Vtrial(r).
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each
!!    G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each
!!    G and k+g point
!!  ylmgr1(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!    for each G and k+g point
!!
!!
!! OUTPUT
!!  cg1(2,mpw*nspinor*mband*mk1mem*nsppol)=updated wavefunctions,
!!    if mk1mem==0 they are on a disk file
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues
!!    (hartree)
!!  edocc=correction to 2nd-order total energy coming from changes of occupation
!!  eeig0=0th-order eigenenergies part of 2nd-order total energy
!!  ek0=0th-order kinetic energy part of 2nd-order total energy.
!!  ek1=1st-order kinetic energy part of 2nd-order total energy
!!    (not for phonons)
!!  eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!!  enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!!  enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!!  resid(mband*nkpt_rbz*nsppol)=residuals for each band over all k points.
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!  rhog1(2,nfftf)=RF electron density in reciprocal space
!!  ==== if optres==1
!!    nres2=square of the norm of the residual
!!    nvresid1(cplex*nfftf,nspden)=1st-order density residual
!!  ==== if psps%usepaw==1
!!    nhat1(nhat1dim,nspden*psps%usepaw)=1st-order compensation charge density
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!  pawrhoij1(dimpaw1) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!  rhor1(cplex*nfftf,nspden)=RF electron density in electrons/bohr**3.
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      clsopn,cprj_alloc,cprj_diskinit,cprj_free,cprj_get,fftpac,fourdp,gbefd3
!!      gradberry3,hdr_io,hdr_io_netcdf,hdr_skip,ini_wf_netcdf,kpg3,kpgstr
!!      leave_new,leave_test,mkffnl,mkkin,mkkpg,occeig,pawmknhat3,ph1d3d,rdnpw
!!      rwwf,sphereboundary,sqnorm_v,status,symrhg,symrhoij3,timab,transgrid
!!      vtowfk3,wffkg,wffreadskipk,wrtout,xcomm_world,xdefineoff,xme_init
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine vtorho3(atindx,atindx1,cg,cgq,cg1,cplex,cprj,cprjq,cpus,dbl_nnsclo,gh1_rbz,densymop_rf,&
& dimcprj,dimpaw1,doccde_rbz,docckqde,dtefield,dtfil,dtset,&
& edocc,eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,enl0,enl1,&
& fermie1,gmet,gprimd,gsqcut,hdr,idir,indsy1,&
& ipert,irrzon1,istep,istwfk_rbz,kg,kg1,kpt_rbz,mband,&
& mkmem,mkqmem,mk1mem,mpi_enreg,mpsang,mpw,mpw1,&
& natom,nattyp,nband_rbz,nfftf,nhat1,nhat1dim,nkpt_rbz,npwarr,npwar1,nres2,nspden,nspinor,&
& nsppol,nsym1,ntypat,nvresid1,occkq,occ_rbz,optres,&
& paw_ij,paw_ij1,pawang,pawfgr,pawfgrtab,pawrhoij,pawrhoij1,pawtab,&
& phnons1,ph1d,prtvol,psps,pwindall,qmat,resid,residm,rhog1,rhor1,rmet,rprimd,symaf1,symrc1,tnons1,ucvol,&
& usecprj,wffddk,wffnew,wffnow,wfftgs,wfftkq,vtrial,vtrial1,wtk_rbz,xred,ylm,ylm1,ylmgr1)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_13io_mpi
 use interfaces_13ionetcdf
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_14occeig
 use interfaces_15common
 use interfaces_16response, except_this_one => vtorho3
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,dbl_nnsclo,dimpaw1,idir,ipert,istep,mband,mk1mem
 integer,intent(in) :: mkmem,mkqmem,mpsang,mpw,mpw1,natom,nfftf,nhat1dim
 integer,intent(in) :: nkpt_rbz,nspden,nsppol,nsym1,ntypat,optres,prtvol
 integer,intent(in) :: usecprj
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: cpus,fermie1,gsqcut,ucvol
 real(dp),intent(out) :: edocc,eeig0,ek0,ek1,eloc0,enl0,enl1,nres2,residm
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(dens_sym_operator_type),intent(in) :: densymop_rf
 type(efield_type) :: dtefield
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffddk,wffnew,wffnow,wfftgs,wfftkq
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),dimcprj(natom)
 integer,intent(in) :: indsy1(4,nsym1,natom)
!no_abirules
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,nspden/nsppol)
 integer,intent(in) :: istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem),nattyp(ntypat)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),npwar1(nkpt_rbz,2)
 integer,intent(in) :: npwarr(nkpt_rbz,2),symaf1(nsym1),symrc1(3,3,nsym1)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(inout) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)
 real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
 real(dp),intent(out) :: gh1_rbz(nkpt_rbz*dtset%ieig2rf,mband,2,mpw1*nspinor),eigen1(2*mband*mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz),occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: occkq(mband*nkpt_rbz*nsppol),ph1d(2,3*(2*dtset%mgfft+1)*natom)
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),nspden/nsppol)
 real(dp), intent(out) :: nhat1(nhat1dim,dtset%nspden*psps%usepaw)
 real(dp),intent(out) :: resid(mband*nkpt_rbz*nsppol),rhog1(2,nfftf)
 real(dp),intent(inout) :: nvresid1(cplex*nfftf,nspden),rhor1(cplex*nfftf,nspden)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: tnons1(3,nsym1)
 real(dp),intent(inout) :: vtrial1(cplex*nfftf,nspden),vtrial(nfftf,nspden)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm)
 integer,intent(in) :: pwindall(max(mpw,mpw1)*mkmem,8,3)
 real(dp),intent(in) :: qmat(2,dtefield%nband_occ,dtefield%nband_occ,nkpt_rbz,2,3)
 type(cprj_type),intent(in) :: cprj(dimpaw1,nspinor*dtset%mband*mkmem*dtset%nsppol*usecprj)
 type(cprj_type),intent(in) :: cprjq(natom,nspinor*dtset%mband*mkqmem*dtset%nsppol*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(dtset%natom*psps%usepaw),paw_ij1(dimpaw1)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(dimpaw1)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=17
 integer,save :: count=0
 integer :: bantot,bd2tot_index,bdtot_index,buffer_size,counter,cplex_dij
 integer :: cplex_dij1,dimdij,dimdij1,dimffnl1,dimffnlk,enunit,fform,formeig,i1
 integer :: i2,i3,ia,iatom,iband,ibdkpt,ibg,ibg1,ibgq,icg,icg1,icgq,ider,idir0
 integer :: ieig2rf,ierr,iexit,ifft,ii,ikg,ikg1,ikpt,ilm,ilmn,iln,iln0,index1
 integer :: iorder_cprj,iproc,ipw,iscf_mod,isp,ispden,isppol,istr,istwf_k,it
 integer :: itypat,matblk,mbd2kpsp,mbdkpsp,mcgq,mcgq_disk,mcprjq,mcprjq_disk,me
 integer :: muig,n1,n2,n3,n4,n5,n6,nband_k,nband_kq,nkpg,nkpg1,nnsclo_now
 integer :: npw1_k,npw_k,nsp,nvloc,option,rdkg,rdwf,rdwr,spaceworld,tag
 integer :: test_dot,tim_rwwf,useekb1
 real(dp) :: arg,valt,wtk_k
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamkq
 type(wffile_type) :: wfftmp
!arrays
 integer :: pspso_typ(1)
 integer,allocatable :: gbound(:,:),indlmn_typ(:,:,:),kg1_k(:,:),kg_dum(:,:)
 integer,allocatable :: kg_k(:,:)
 real(dp) :: cg_dum(2,1),eig_dum(1),kpoint(3),kpq(3),occ_dum(1),rhodum(1)
 real(dp) :: tsec(2),ylmgr_dum(1)
 real(dp),allocatable :: buffer1(:),cgq_disk(:,:),cgrvtrial(:,:),dkinpw(:)
 real(dp),allocatable :: doccde_k(:),doccde_kq(:),edocc_k(:),eeig0_k(:)
 real(dp),allocatable :: eig0_k(:),eig0_kq(:),eig1_k(:),eigkq_dum(:),ek0_k(:)
 real(dp),allocatable :: ek1_k(:),ekb1_typ(:,:,:),ekb_typ(:,:,:),eloc0_k(:)
 real(dp),allocatable :: enl0_k(:),enl1_k(:),ffnl1(:,:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: ffnlkq(:,:,:,:),ffspl_typ(:,:,:,:),gh1_k(:,:,:)
 real(dp),allocatable :: grad_berry(:,:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:)
 real(dp),allocatable :: occ_k(:),occ_kq(:),occkq_dum(:),ph3d(:,:,:),resid_k(:)
 real(dp),allocatable :: rho1wfg(:,:),rho1wfr(:,:),rhoaug1(:,:,:),rocceig(:,:)
 real(dp),allocatable :: sij_typ(:,:),vlocal(:,:,:),vlocal1(:,:,:),vt1g(:,:)
 real(dp),allocatable :: wavef(:,:),ylm1_k(:,:),ylm_k(:,:),ylmgr1_k(:,:,:)
 type(cprj_type),allocatable :: cprj1(:,:),cprjq_disk(:,:)

! *********************************************************************

!DEBUG
!write(6,*)' vtorho3 : enter '
!count=count+1
!write(6,*)' count=',count
!if(count==13)stop
!write(6,*)' vtorho3 :  vtrial1(1,1)=',vtrial1(1,1)
!write(6,*)' vtorho3 : kg1(:,66)',kg1(:,66)
!write(6,*)' xred=',xred
!stop
!ENDDEBUG

 if (dtset%berryopt == 4) then
  allocate(grad_berry(2,mpw1,dtefield%nband_occ))
 end if

!Keep track of total time spent in vtorho3
 call timab(121,1,tsec)
 call timab(124,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Init mpi_comm
!BEGIN TF_CHANGES
 call xcomm_world(mpi_enreg,spaceworld)
!END TF_CHANGES
!Init me
 call xme_init(mpi_enreg,me)


!Structured debugging if prtvol==-level
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' vtorho3 : enter '
  call wrtout(06,message,'COLL')
 end if

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' vtorho3 :  BUG -',ch10,&
&  '  wrong values for nfft, nfftf !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Debugging : print vtrial1
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(prtvol==-level)then
  write(message,'(a)') '   ifft  vtrial1(ifft) '
  call wrtout(06,message,'COLL')
  do ifft=1,nfftf
   if(ifft<=11 .or. mod(ifft,301)==0 )then
    write(message,'(i5,a,es13.6)')ifft,' ',vtrial1(ifft,1)
    call wrtout(06,message,'COLL')
    if(nspden==2)then
     write(message,'(a,es13.6)')'      ',vtrial1(ifft,2)
     call wrtout(06,message,'COLL')
    end if
   end if
  end do
 end if

 iscf_mod=dtset%iscf

!The value of iscf must be modified if ddk perturbation, see loper3.f
 if(ipert==natom+1 ) iscf_mod=-3

 edocc=zero ; eeig0=zero ; ek0=zero  ; ek1=zero
 eloc0=zero ; enl0=zero ; enl1=zero
 bdtot_index=0
 bd2tot_index=0
 ibg=0;icg=0
 ibgq=0;icgq=0
 ibg1=0;icg1=0
 mbdkpsp=mband*nkpt_rbz*nsppol
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol

 allocate(kg_k(3,mpw))
 allocate(kg1_k(3,mpw1))

!Initialize PW 1st-order density if needed
!Also store old rho1 in case of density mixing
 if (iscf_mod>0) then
  if (optres==1) nvresid1=rhor1
  if (psps%usepaw==0) then
   rhor1(:,:)=zero
  else
   allocate(rho1wfr(cplex*dtset%nfft,dtset%nspden),rho1wfg(2,dtset%nfft))
   rho1wfr(:,:)=zero
  end if
 end if

!Set max number of non-self-consistent loops nnsclo_now for use in vtowfk3
 if(iscf_mod<=0 .and. iscf_mod/=-3)then
  nnsclo_now=dtset%nstep
 else
  if(dtset%nnsclo>0)then
   nnsclo_now=dtset%nnsclo
  else
   nnsclo_now=1
  end if
  if(dbl_nnsclo==1)then
!  DEBUG
!  write(6,*)' vtorho : use doubled nnsclo '
!  ENDDEBUG
   nnsclo_now=nnsclo_now*2
  end if
 end if

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 allocate(rhoaug1(cplex*n4,n5,n6),vlocal(n4,n5,n6),vlocal1(cplex*n4,n5,n6))

 allocate(kg_dum(3,0))

!Prepare GS k wf file for reading if mkmem==0
 if (mkmem==0) then
  call clsopn(wfftgs)
  call hdr_skip(wfftgs,ierr)
! Define offsets, in case of MPI I/O
  formeig=0
  call WffKg(wfftgs,1)
  call xdefineOff(formeig,wfftgs,mpi_enreg,nband_rbz,npwarr,nspinor,nsppol,nkpt_rbz)
 end if

!Prepare GS k+q wf file for reading if mkqmem==0
 if (mkqmem==0) then
  call clsopn(wfftkq)
  call hdr_skip(wfftkq,ierr)
  mcgq_disk=mpw1*nspinor*mband
  allocate(cgq_disk(2,mcgq_disk))
! Define offsets, in case of MPI I/O
  formeig=0
  call WffKg(wfftkq,1)
  call xdefineOff(formeig,wfftkq,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)
 else
  mcgq=mpw1*nspinor*mband*mkqmem*nsppol
 end if

!Prepare RF wf files for reading and writing if mk1mem==0
 if (mk1mem==0) then
  call clsopn(wffnow)
! Read unwfnow header
  call hdr_skip(wffnow,ierr)
! Define offsets, in case of MPI I/O
  formeig=1
  call WffKg(wffnow,1)
  call xdefineOff(formeig,wffnow,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)
  call clsopn(wffnew)
! Write the content of hdr to the new wf file
  rdwr=2 ; fform=2
  if (wffnew%accesswff /= 2) then
   call hdr_io(fform,hdr,rdwr,wffnew)
#if defined HAVE_NETCDF
  else if (wffnew%accesswff == 2) then
   call hdr_io_netcdf(fform,hdr,rdwr,wffnew)
   call ini_wf_netcdf(mpw1,wffnew%unwff,1)
#endif
  end if
! Define offsets, in case of MPI I/O
  formeig=1
  call WffKg(wffnew,1)
  call xdefineOff(formeig,wffnew,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)
 end if

!Prepare RF PAW files for reading and writing if mkmem, mkqmem or mk1mem==0
 if (psps%usepaw==1) then
  iorder_cprj=0
  call cprj_diskinit(atindx1,dimpaw1,iorder_cprj,mkmem ,natom,dimcprj,nspinor,dtfil%unpaw)
  call cprj_diskinit(atindx1,natom  ,iorder_cprj,mkqmem,natom,dimcprj,nspinor,dtfil%unpawq)
  if (dtset%mkqmem==0) then
   mcprjq_disk=nspinor*dtset%mband*usecprj
   allocate(cprjq_disk(natom,mcprjq_disk))
   if (mcprjq_disk>0) then
    call cprj_alloc(cprjq_disk,0,dimcprj)
   end if
  else
   mcprjq=nspinor*dtset%mband*mkqmem*nsppol*usecprj
  end if
  allocate(cprj1(dimpaw1,nspinor*dtset%mband*dtset%mk1mem*dtset%nsppol*usecprj))
  if (usecprj==1) then
   if (dtset%mk1mem/=0) then
    call cprj_alloc(cprj1,0,dimcprj)
   else
    rewind dtfil%unpaw1
    write(dtfil%unpaw1) dimpaw1,nspinor
    write(dtfil%unpaw1) dimcprj(1:dimpaw1)
   end if
  end if
 end if

!Initialisation of the wfdot file in case of electric field
 test_dot=0
 if( (ipert==natom+2 .or. ipert==natom+5) .and. &
& sum( (dtset%qptn(1:3))**2 ) <= tol7 .and. (dtset%berryopt/=4) )then
  test_dot=1
  call clsopn(wffddk)
  call hdr_skip(wffddk,ierr)
 end if

 nvloc=1

!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 allocate(gs_hamkq%atindx(natom),gs_hamkq%atindx1(natom))
 allocate(gs_hamkq%gbound(2*dtset%mgfft+8,2))
 allocate(gs_hamkq%indlmn(6,psps%lmnmax,ntypat))
 allocate(gs_hamkq%nattyp(ntypat))
 allocate(gs_hamkq%phkxred(2,natom))
 allocate(gs_hamkq%ph1d(2,3*(2*dtset%mgfft+1)*natom))
 allocate(gs_hamkq%pspso(ntypat))
 allocate(gs_hamkq%xred(3,natom))

!Initialize most components of the Ground-state Hamiltonian ar k+q
 gs_hamkq%atindx(:)  =atindx(:)
 gs_hamkq%atindx1(:) =atindx1(:)
 gs_hamkq%gmet(:,:)  =gmet(:,:)
 gs_hamkq%gprimd(:,:)=gprimd(:,:)
 gs_hamkq%indlmn(:,:,:)=psps%indlmn(:,:,:)
 gs_hamkq%lmnmax     =psps%lmnmax
 gs_hamkq%mgfft      =dtset%mgfft
 gs_hamkq%mpsang     =mpsang
 gs_hamkq%mpssoang   =psps%mpssoang
 gs_hamkq%natom      =natom
 gs_hamkq%nattyp(:)  =nattyp(:)
 gs_hamkq%nfft       =dtset%nfft
 gs_hamkq%ngfft(:)   =dtset%ngfft(:)
 gs_hamkq%nloalg(:)  =dtset%nloalg(:)
 gs_hamkq%nspinor    =nspinor
 gs_hamkq%ntypat     =ntypat
 gs_hamkq%nvloc      =nvloc
 gs_hamkq%n4         =n4
 gs_hamkq%n5         =n5
 gs_hamkq%n6         =n6
 gs_hamkq%usepaw     =psps%usepaw
 gs_hamkq%ph1d(:,:)  =ph1d(:,:)
 gs_hamkq%pspso(:)   =psps%pspso(1:ntypat)
 gs_hamkq%ucvol      =ucvol
 gs_hamkq%useylm     =psps%useylm
 gs_hamkq%xred(:,:)  =xred(:,:)

!Non-local factors:
! Norm-conserving: kleimann-Bylander energies
! PAW: Dij coefficients and overlap coefficients
 if (psps%usepaw==0) then
  gs_hamkq%dimekb1=psps%dimekb
  gs_hamkq%dimekb2=ntypat
  allocate(gs_hamkq%ekb(psps%dimekb,ntypat,nspinor**2))
  allocate(gs_hamkq%sij(0,0))
  gs_hamkq%ekb(:,:,1)=psps%ekb(:,:)
  if (nspinor==2) then
   gs_hamkq%ekb(:,:,2)=psps%ekb(:,:)
   gs_hamkq%ekb(:,:,3:4)=zero
  end if
 else
  cplex_dij=paw_ij(1)%cplex_dij
  gs_hamkq%dimekb1=psps%dimekb*cplex_dij
  gs_hamkq%dimekb2=natom
  allocate(gs_hamkq%ekb(gs_hamkq%dimekb1,gs_hamkq%dimekb2,nspinor**2))
  allocate(gs_hamkq%sij(gs_hamkq%dimekb1,ntypat))
  do itypat=1,ntypat
   if (cplex_dij==1) then
    gs_hamkq%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
   else
    do ilmn=1,pawtab(itypat)%lmn2_size
     gs_hamkq%sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
     gs_hamkq%sij(2*ilmn  ,itypat)=zero
    end do
   end if
  end do
 end if

!PAW:allocate memory for non-symetrized 1st-order occupancies matrix (pawrhoij1)
 if (psps%usepaw==1.and.iscf_mod>0) then
  nsp=dtset%nsppol;if (dtset%nspden==4) nsp=4
  do iatom=1,dimpaw1
   allocate(pawrhoij1(iatom)%rhoij_(pawrhoij1(iatom)%lmn2_size,nsp))
   pawrhoij1(iatom)%use_rhoij_=1
   pawrhoij1(iatom)%rhoij_(:,:)=zero
  end do
 end if

!LOOP OVER SPINS
 do isppol=1,nsppol

  if (nsppol/=1) then
   write(message,*)' ****  In vtorho3 for isppol=',isppol
   call wrtout(06,message,'COLL')
  end if

! Rewind kpgsph data file if needed:
  if (mkmem==0) rewind(dtfil%unkg)
  if (mk1mem==0) rewind(dtfil%unkg1)
  if (mkmem==0.and.psps%useylm==1) rewind(dtfil%unylm)
  if (mk1mem==0.and.psps%useylm==1) rewind(dtfil%unylm1)
  ikg=0;ikg1=0

! Set up local potential vlocal1 with proper dimensioning, from vtrial1
! Same thing for vlocal from vtrial
! Also take into account the spin.
  if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
   call fftpac(isppol,nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal,2)
   call fftpac(isppol,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,vtrial1,vlocal1,2)
  else
   allocate(cgrvtrial(dtset%nfft,nspden))
   call transgrid(1,mpi_enreg,nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
   call fftpac(isppol,nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
   deallocate(cgrvtrial)
   allocate(cgrvtrial(cplex*dtset%nfft,nspden))
   call transgrid(cplex,mpi_enreg,nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial1)
   call fftpac(isppol,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal1,2)
   deallocate(cgrvtrial)
  end if
  rhoaug1(:,:,:)=zero

! DEBUG
! write(6,*)' vtorho3 : after fftpac vtrial1(1,1)=',vtrial1(1,1)
! write(6,*)' vtorho3 : after fftpac vlocal1(1,1,1)=',vlocal1(1,1,1)
! ENDDEBUG

! PAW: retrieve Dij coefficients for this spin component
  if (psps%usepaw==1) then
   do ispden=1,nspinor**2
    isp=isppol;if (nspinor==2) isp=ispden
    do iatom=1,natom
     dimdij=cplex_dij*paw_ij(iatom)%lmn2_size
     do ilmn=1,dimdij
      gs_hamkq%ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
     end do
     if(dimdij+1<=gs_hamkq%dimekb1) gs_hamkq%ekb(dimdij+1:gs_hamkq%dimekb1,iatom,ispden)=zero
    end do
   end do
  end if

  call timab(125,1,tsec)

! BIG FAT k POINT LOOP
  do ikpt=1,nkpt_rbz

   counter=100*ikpt+isppol
   call status(counter,dtfil%filstat,iexit,level,'loop ikpt     ')
   nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
   istwf_k=istwfk_rbz(ikpt)
   npw_k=npwarr(ikpt,1)
   npw1_k=npwar1(ikpt,1)

   if(mpi_enreg%paral_compil_kpt==1)then
    if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
     eigen1(1+bd2tot_index : 2*nband_k**2+bd2tot_index) = zero
     resid(1+bdtot_index : nband_k+bdtot_index) = zero
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2

     if (test_dot==1) then
!     Skip the wavefunction block for ikpt,isppol
      tim_rwwf=0
      call WffReadSkipK(1,0,ikpt,isppol,mpi_enreg,wffddk)
!     End the treatment of the dot file
     end if

!    Skip the rest of the k-point loop
     cycle
    end if
   end if

   allocate(ylm_k(npw_k,mpsang*mpsang*psps%useylm))
   allocate(ylm1_k(npw1_k,mpsang*mpsang*psps%useylm))
   if (ipert==natom+1 .or. &
&      ipert==natom+3 .or. ipert==natom+4) &
&       allocate(ylmgr1_k(npw1_k,3,mpsang*mpsang*psps%useylm))

!  Continue to initialize the Hamiltonian at k+q
   kpoint(:)=kpt_rbz(:,ikpt)
   kpq(:)=kpoint(:);if (ipert<natom+3) kpq(:)=kpq(:)+dtset%qptn(1:3)
   gs_hamkq%istwf_k    =istwf_k
   gs_hamkq%kpoint(:)  =kpq(:)
   gs_hamkq%npw        =npw1_k

   allocate(doccde_k(nband_k),doccde_kq(nband_k))
   allocate(eig0_k(nband_k),eig0_kq(nband_k),eig1_k(2*nband_k**2))
   allocate(edocc_k(nband_k),eeig0_k(nband_k))
   allocate(ek0_k(nband_k),ek1_k(nband_k))
   allocate(eloc0_k(nband_k))
   allocate(occ_k(nband_k),occ_kq(nband_k))
   allocate(resid_k(nband_k))
   allocate(gbound(2*dtset%mgfft+8,2))
   allocate(rocceig(nband_k,nband_k))
   allocate(enl0_k(nband_k),enl1_k(nband_k))

   eig1_k(:)=zero
   eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
   eig0_kq(:)=eigenq(1+bdtot_index:nband_k+bdtot_index)
   edocc_k(:)=zero
   eeig0_k(:)=zero ; ek0_k(:)=zero  ; ek1_k(:)=zero
   eloc0_k(:)=zero ; enl0_k(:)=zero ; enl1_k(:)=zero
   occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
   occ_kq(:)=occkq(1+bdtot_index:nband_k+bdtot_index)
   doccde_k(:)=doccde_rbz(1+bdtot_index:nband_k+bdtot_index)
   doccde_kq(:)=docckqde(1+bdtot_index:nband_k+bdtot_index)
   resid_k(:)=zero

!  For each pair of active bands (m,n), generates the ratios
!  rocceig(m,n)=(occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
!  and decide to which band to attribute it.
   call occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,&
&   dtset%occopt,occ_k,occ_kq,rocceig)

   if (mkmem==0) then
!   Read (k+G) basis sphere data (same for each spin)
    call status(counter,dtfil%filstat,iexit,level,'read kg data  ')
    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,dtfil%unkg)
!   Read k+g data
    read (dtfil%unkg) kg_k(1:3,1:npw_k)
    call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)

!   Eventually read (k+G) spherical harmonics
    if (psps%useylm==1) then
     read(dtfil%unylm)
     read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
    end if

   else

    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
      ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     end do
    end if

!   End if for choice governed by mkmem
   end if

   if (mkqmem==0) then
    call status(counter,dtfil%filstat,iexit,level,'read k+q wfs  ')
!   Take care of GS wavefunctions at k+q , nothing to be done for kpg sphere
!   Eigenvalues already in eigen1
    tim_rwwf=17
    allocate(eigkq_dum(mband),occkq_dum(mband))
    call rwwf(cgq_disk,eigkq_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcgq_disk,mpi_enreg,nband_k,nband_k,&
&    npw1_k,nspinor,occkq_dum,-2,0,tim_rwwf,wfftkq)
    deallocate(eigkq_dum,occkq_dum)
!   PAW: Take care of cprj at k+q
    if (psps%usepaw==1.and.usecprj==1) then
     call cprj_get(atindx1,cprjq_disk,cprjq,natom,1,ibgq,ikpt,iorder_cprj,isppol,&
&         mband,mkqmem,mpi_enreg,natom,nband_k,nband_k,nspinor,nsppol,dtfil%unpawq)
    end if

!   End if for choice governed by mkqmem
   end if

   wtk_k=wtk_rbz(ikpt)

!  DEBUG
!  write(6,*)' vtorho3 : wtk_k=',wtk_k
!  ENDDEBUG

   kg1_k(:,:) = 0
   if (mk1mem==0) then

!   Read (k+q+G) basis sphere data (same for each spin)
    call status(counter,dtfil%filstat,iexit,level,'read kg1 data ')
    call rdnpw(ikpt,isppol,nband_k,npw1_k,nspinor,0,dtfil%unkg1)
    read (dtfil%unkg1) kg1_k(1:3,1:npw1_k)
    call sphereboundary(gs_hamkq%gbound,istwf_k,kg1_k,dtset%mgfft,npw1_k)

!   Eventually read (k+q+G) spherical harmonics
    if (psps%useylm==1) then
     read(dtfil%unylm1)
     if (ipert==natom+1 .or. &
&        ipert==natom+3 .or. ipert==natom+4) then
      read(dtfil%unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang),&
&      (((ylmgr1_k(muig,ii,ilm),muig=1,npw1_k),ii=1,3),ilm=1,mpsang*mpsang)
     else
      read(dtfil%unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang)
     end if
    end if

   else

    kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
    call sphereboundary(gs_hamkq%gbound,istwf_k,kg1_k,dtset%mgfft,npw1_k)
    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
      ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
     end do
     if (ipert==natom+1 .or. &
&        ipert==natom+3 .or. ipert==natom+4) then
      do ilm=1,mpsang*mpsang
       do ii=1,3
        ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
       end do
      end do
     end if
    end if

!   End if for choice governed by mk1mem
   end if

!  Set up the ground-state Hamiltonian, and some parts of the 1st order Hamiltonian

!  Note that not all these arrays should be allocated in the general case
!  when wtk_k vanishes
   allocate(indlmn_typ(6,psps%lmnmax,1))
   allocate(dkinpw(npw_k),kinpw1(npw1_k))
   dimffnlk=1;allocate(ffnlk(npw_k,dimffnlk,psps%lmnmax,1))
   if (psps%usepaw==1) then
    allocate(ekb_typ(psps%dimekb*cplex_dij,1,nspinor**2))
    cplex_dij1=paw_ij1(1)%cplex_dij
    if (ipert>=1.and.ipert<=natom) cplex_dij1=paw_ij1(ipert)%cplex_dij
    useekb1=1;allocate(ekb1_typ(psps%dimekb*cplex_dij1,1,nspinor**2))
    allocate(sij_typ(psps%dimekb*cplex_dij1,1))
   else
    allocate(ekb_typ(psps%dimekb,1,nspinor**2))
    useekb1=0
   end if

!  Compute (k+G) vectors (only if useylm=1)
   nkpg=0;if(ipert>=1.and.ipert<=natom) nkpg=3*dtset%nloalg(5)
   allocate(kpg_k(npw_k,nkpg))
   if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!  Preparation of the kinetic and non-local contributions
   if( ipert>=1 .and. ipert<=natom )then

!   Compute nonlocal form factors ffnlk at (k+G), for the displaced atom only.
    call status(counter,dtfil%filstat,iexit,level,'call mkffnl(0)')
    ider=0;idir0=0
!   Need to transfer the infos relative to the displaced atom: ekb, indlmn, pspso, ffspl
    itypat=dtset%typat(ipert)
    if (psps%usepaw==1) then
     dimdij =cplex_dij *paw_ij(ipert)%lmn2_size
     dimdij1=cplex_dij1*paw_ij1(1   )%lmn2_size
     do ispden=1,nspinor**2
      isp=isppol;if (nspinor==2) isp=ispden
      ekb_typ(1:dimdij,1,ispden)=paw_ij(ipert)%dij(1:dimdij,isp)
      if (cplex_dij==1) then
       sij_typ(1:dimdij,1)=pawtab(itypat)%sij(1:dimdij)
      else
       do ilmn=1,dimdij/2
        sij_typ(2*ilmn-1,1)=pawtab(itypat)%sij(ilmn)
        sij_typ(2*ilmn  ,1)=zero
       end do
      end if
      ekb1_typ(1:dimdij1,1,ispden)=paw_ij1(1)%dij(1:dimdij1,isp)
      if(dimdij+1<=psps%dimekb*cplex_dij) then
       ekb_typ(dimdij+1:psps%dimekb*cplex_dij,1,ispden)=zero
       sij_typ(dimdij+1:psps%dimekb*cplex_dij,1)=zero
      end if
      if(dimdij1+1<=psps%dimekb*cplex_dij1) then
       ekb1_typ(dimdij+1:psps%dimekb*cplex_dij1,1,ispden)=zero
      end if
     end do
    else
     ekb_typ(:,1,1)=psps%ekb(:,itypat)
     if (nspinor==2) then
      ekb_typ(:,1,2)=psps%ekb(:,itypat)
      ekb_typ(:,1,3:4)=zero
     end if
    end if
    indlmn_typ(:,:,1)=psps%indlmn(:,:,itypat)
    pspso_typ(1)=psps%pspso(itypat)
    allocate(ffspl_typ(psps%mqgrid_ff,2,psps%lnmax,1))
    ffspl_typ(:,:,:,1)=psps%ffspl(:,:,:,itypat)
    call mkffnl(psps%dimekb,dimffnlk,ekb_typ(:,1,1),ffnlk,ffspl_typ,&
&    gmet,gprimd,ider,idir0,indlmn_typ,kg_k,kpg_k,kpoint,psps%lmnmax,&
&    psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,1,&
&    pspso_typ,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
    deallocate(ffspl_typ)

   else if(ipert==natom+1)then

!   This flag is needed for the call to mkffnl later
    ider=1;idir0=idir

    call status(counter,dtfil%filstat,iexit,level,'call kpg3     ')
!   Compute the derivative of the kinetic operator vs k in kinpw
    call kpg3(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,idir,kg_k,kpoint,npw_k)
!   No non-local form factor is to be calculated here in this case

   else if(ipert==natom+2 )then

    ider=0;idir0=0
!   dkinpw and ffnlk are not needed ...

!   section for strain perturbation
   else if(ipert==natom+3 .or. ipert==natom+4)then

!   istr is 1,2,...,6 and indicates cartesian strain component
    if(ipert==natom+3) then
     istr=idir
    else
     istr=idir+3
    end if

!   Derivatives needed for strain perturbation
!   This flag is needed for the call to mkffnl later
    ider=1;idir0=-istr

    call status(counter,dtfil%filstat,iexit,level,'call kpgstr   ')
!   Compute the derivative of the kinetic operator vs strain in dkinpw
    call kpgstr(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,gprimd,istr,&
&    kg_k,kpoint,npw_k)

!   endsection for strain perturbation

! magnetic field case (set equal to electric field for testing at the moment)
   else if(ipert==natom+5 )then

    ider=0;idir0=0
!   dkinpw and ffnlk are not needed ...

   end if

!  Compute (1/2) (2 Pi)**2 (k+q+G)**2:
   call status(counter,dtfil%filstat,iexit,level,'call mkkin(1) ')
   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg1_k,kinpw1,kpq,npw1_k)

!  Compute (k+q+G) vectors (only if useylm=1)
   nkpg1=0;if(ipert>=1.and.ipert<=natom ) nkpg1=3*dtset%nloalg(5)
   allocate(kpg1_k(npw1_k,nkpg1))
   if (nkpg1>0) call mkkpg(kg1_k,kpg1_k,gs_hamkq%kpoint(:),nkpg1,npw1_k)

!  Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
   call status(counter,dtfil%filstat,iexit,level,'call mkffnl(1)')
   dimffnl1=1+ider;if (ider==1.and.idir0==0) dimffnl1=dimffnl1+2*psps%useylm
   allocate(ffnlkq(npw1_k,dimffnl1,psps%lmnmax,1))
   allocate(ffnl1(npw1_k,dimffnl1,psps%lmnmax,ntypat))
   call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gmet,gprimd,ider,idir0,&
&   psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
&   npw1_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)

   if( ipert>=1 .and. ipert<=natom )then
!   Copy the part needed for the displaced atom, in ffnlkq.
    do ilmn=1,psps%lmnmax
     do ii=1,dimffnl1
      do ipw=1,npw1_k
       ffnlkq(ipw,ii,ilmn,1)=ffnl1(ipw,ii,ilmn,itypat)
      end do
     end do
    end do
   end if

!  Allocate the arrays phkxred and ph3d, compute phkxred
!  and eventually ph3d.
!  NOTE : in this RF case, uses kpq instead of kpt
   call status(counter,dtfil%filstat,iexit,level,'make phkxred  ')
   do ia=1,natom
    iatom=atindx(ia)
    arg=two_pi*(kpq(1)*xred(1,ia)+kpq(2)*xred(2,ia)+kpq(3)*xred(3,ia))
    gs_hamkq%phkxred(1,iatom)=cos(arg)
    gs_hamkq%phkxred(2,iatom)=sin(arg)
   end do

!  Note : use npw1_k
   call status(counter,dtfil%filstat,iexit,level,'make ph3d     ')
   if(dtset%nloalg(1)<=0)then
!   Here, only the allocation, not the precomputation.
    matblk=dtset%nloalg(4)
    allocate(ph3d(2,npw1_k,matblk))
   else
!   Here, allocation as well as precomputation
    matblk=natom
    allocate(ph3d(2,npw1_k,matblk))
    call ph1d3d(1,natom,kg1_k,kpq,matblk,natom,npw1_k,n1,n2,n3,&
&    gs_hamkq%phkxred,ph1d,ph3d)
   end if
   gs_hamkq%matblk=matblk

   call status(counter,dtfil%filstat,iexit,level,'call vtowfk3  ')

!  Compute the gradient of the Berry-phase term
   if (dtset%berryopt == 4) then
    if (ipert<=natom) then
!    phonon perturbation
     call gradberry3(cg,cg1,dtefield,grad_berry,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt_rbz,&
&     npwarr,npwar1,nspinor,nsppol,qmat,pwindall)
    else
!    electric field perturbation
     call gbefd3(cg,cg1,dtefield,grad_berry,idir,ikpt,isppol,mband,mpw,mpw1,mkmem,mk1mem,nkpt_rbz,&
&     npwarr,npwar1,nspinor,&
&     nsppol,qmat,pwindall,rprimd)
    end if
   end if

!  Compute the eigenvalues, wavefunction, residuals,
!  contributions to kinetic energy, nonlocal energy, forces,
!  and update of 1st-order density to this k-point and this spin polarization.
   allocate(gh1_k(nband_k,2,mpw1*nspinor))
   if(mkqmem/=0)then
    nband_kq = nband_k  !Note that the calculation only works for same number of bandes on all K points.
!   Note that vtowfk3 is called with kpoint, while kpt is used inside vtowfk3
    call vtowfk3(cg,cgq,cg1,cplex,cprj,cprjq,cprj1,cpus,gh1_k,dimcprj,&
&    psps%dimekb,dimffnlk,dimffnl1,dimpaw1,dkinpw,dtfil,dtset,&
&    edocc_k,eeig0_k,eig0_k,eig0_kq,eig1_k,ekb_typ,ekb1_typ,&
&    ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k,fermie1,ffnlk,ffnlkq,ffnl1,&
&    gbound,grad_berry,gs_hamkq,ibg,ibgq,ibg1,icg,icgq,icg1,idir,ikpt,indlmn_typ,ipert,&
&    isppol,istep,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpoint,psps%lmnmax,&
&    matblk,mband,mcgq,mcprjq,dtset%mgfft,mkmem,mkqmem,mk1mem,mpi_enreg,&
&    psps%mpsang,psps%mpssoang,mpw,mpw1,natom,nband_k,nband_kq,nkpg,nkpg1,&
&    nkpt_rbz,nnsclo_now,npw_k,npw1_k,nspden,nspinor,nsppol,&
&    ntypat,n4,n5,n6,occ_k,pawrhoij1,ph3d,prtvol,psps,pspso_typ,resid_k,rhoaug1,rocceig,&
&    sij_typ,usecprj,useekb1,wffddk,wffnew,wffnow,wfftgs,vlocal,vlocal1,wtk_k)

   else if(mkqmem==0)then
    nband_kq = nband_k  !Note that the calculation only works for same number of bandes on all K points.
!   Note that vtowfk3 is called with kpoint, while kpt is used inside vtowfk3
    call vtowfk3(cg,cgq_disk,cg1,cplex,cprj,cprjq_disk,cprj1,cpus,gh1_k,dimcprj,&
&    psps%dimekb,dimffnlk,dimffnl1,dimpaw1,dkinpw,dtfil,dtset,&
&    edocc_k,eeig0_k,eig0_k,eig0_kq,eig1_k,ekb_typ,ekb1_typ,ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k,&
&    fermie1,ffnlk,ffnlkq,ffnl1,gbound,grad_berry,gs_hamkq,&
&    ibg,ibgq,ibg1,icg,icgq,icg1,idir,ikpt,indlmn_typ,ipert,&
&    isppol,istep,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpoint,psps%lmnmax,&
&    matblk,mband,mcgq_disk,mcprjq_disk,dtset%mgfft,mkmem,mkqmem,mk1mem,mpi_enreg,&
&    psps%mpsang,psps%mpssoang,mpw,mpw1,natom,nband_k,nband_kq,nkpg,nkpg1,&
&    nkpt_rbz,nnsclo_now,npw_k,npw1_k,nspden,nspinor,nsppol,&
&    ntypat,n4,n5,n6,occ_k,pawrhoij1,ph3d,prtvol,psps,pspso_typ,resid_k,rhoaug1,rocceig,&
&    sij_typ,usecprj,useekb1,wffddk,wffnew,wffnow,wfftgs,vlocal,vlocal1,wtk_k)
   end if

   ieig2rf = dtset%ieig2rf
   if(ieig2rf==1)then
    do iband=1,nband_k
     gh1_rbz(ikpt,iband,:,:)= gh1_k(iband,:,:)
    end do
   end if

!  DEBUG
!  write(6,*)' vtorho3 : after vtowfk3 '
!  if(count==13)stop
!  write(6,*)' vtorho3 :  vtrial1(1,1)=',vtrial1(1,1)
!  write(6,*)' vtorho3 : kg1(:,66)',kg1(:,66)
!  write(6,*)' xred=',xred
!  stop
!  ENDDEBUG

   deallocate(gh1_k)
   deallocate(dkinpw)
   deallocate(ekb_typ)
   deallocate(ffnlk,kpg_k)
   deallocate(indlmn_typ)
   deallocate(ffnl1,kpg1_k)
   deallocate(ffnlkq)
   if (psps%usepaw==1) deallocate(ekb1_typ,sij_typ)
   deallocate(gbound,kinpw1,ph3d)
   call status(counter,dtfil%filstat,iexit,level,'after vtowfk3 ')

!  Save eigenvalues (hartree), residuals (hartree**2)
   eigen1 (1+bd2tot_index : 2*nband_k**2+bd2tot_index) = eig1_k(:)
   resid  (1+bdtot_index : nband_k+bdtot_index) = resid_k(:)

   if ((iscf_mod>0 .or. iscf_mod==-3 .or. iscf_mod==-2).and.&
&   dtset%dsifkpt(1)==1.and.dtset%dsifkpt(2)==1.and.dtset%dsifkpt(3)==1) then
!   Accumulate sum over k points for nonlocal and kinetic energies,
!   also accumulate gradients of Enonlocal:
    do iband=1,nband_k
     edocc=edocc+wtk_k*occ_k(iband)*edocc_k(iband)
     eeig0=eeig0+wtk_k*occ_k(iband)*eeig0_k(iband)
     ek0=ek0+wtk_k*occ_k(iband)*ek0_k(iband)
     ek1=ek1+wtk_k*occ_k(iband)*ek1_k(iband)
     eloc0=eloc0+wtk_k*occ_k(iband)*eloc0_k(iband)
     enl0=enl0+wtk_k*occ_k(iband)*enl0_k(iband)
     enl1=enl1+wtk_k*occ_k(iband)*enl1_k(iband)
    end do
   end if

!  DEBUG

!  write(6,*)' vtorho3 : ikpt=',ikpt,ek0
!  write(6,*)'  iband,wtk_k,occ_k,edocc_k='
!  do iband=1,nband_k
!  write(6, '(i4,3es22.12)' )iband,wtk_k,occ_k(iband),edocc_k(iband)
!  end do
!  ENDDEBUG

   deallocate(doccde_k,doccde_kq,eig0_k,eig0_kq,eig1_k,occ_k,occ_kq,resid_k)
   deallocate(edocc_k,eeig0_k,ek0_k,ek1_k,eloc0_k,enl0_k,enl1_k)
   deallocate(rocceig)

!  Keep track of total number of bands (all k points so far, even for
!  k points not treated by me)
   bdtot_index=bdtot_index+nband_k
   bd2tot_index=bd2tot_index+2*nband_k**2

!  Shift array memory
   if (mkmem/=0) then
    ibg=ibg+nband_k
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   end if
   if (mkqmem/=0) then
    ibgq=ibgq+nband_k
    icgq=icgq+npw1_k*nspinor*nband_k
   end if
   if (mk1mem/=0) then
    ibg1=ibg1+nband_k
    icg1=icg1+npw1_k*nspinor*nband_k
    ikg1=ikg1+npw1_k
   end if
   deallocate(ylm_k,ylm1_k)
   if (ipert==natom+1 .or. &
&      ipert==natom+3 .or. ipert==natom+4) deallocate(ylmgr1_k)
!  End big k point loop
  end do

  call timab(125,2,tsec)

! Transfer density on augmented fft grid to normal fft grid in real space
! Also take into account the spin.
  if(iscf_mod>0) then
   if (psps%usepaw==0) then
    call fftpac(isppol,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rhor1  ,rhoaug1,1)
   else
    call fftpac(isppol,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rho1wfr,rhoaug1,1)
   end if
  end if

! End loop over spins
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
  call timab(166,1,tsec)
! BEGIN TF_CHANGES
  call leave_test(mpi_enreg)
! END TF_CHANGES
  write(message,*) 'vtorho3: loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
  call timab(166,2,tsec)
 end if

 deallocate(gs_hamkq%atindx,gs_hamkq%atindx1)
 deallocate(gs_hamkq%ekb,gs_hamkq%sij)
 deallocate(gs_hamkq%gbound)
 deallocate(gs_hamkq%indlmn)
 deallocate(gs_hamkq%nattyp)
 deallocate(gs_hamkq%phkxred)
 deallocate(gs_hamkq%pspso)
 deallocate(gs_hamkq%ph1d)
 deallocate(gs_hamkq%xred)

!DEBUG
!if(cplex==1)then
!valt=zero
!do i3=1,n3
!do i2=1,n2
!do i1=1,n1
!valt=valt+vlocal1(i1,i2,i3)*rhoaug1(i1,i2,i3)
!end do
!end do
!end do
!Local potential energy of this band, valuer has been accumulated
!valt=0.5_dp*valt*gs_hamkq%ucvol/gs_hamkq%nfft
!end if

!write(6,*)' vtorho3 : eloc1_k=',two*valt

!if(cplex==1)then
!valt=zero
!do ii=1,nfftf
!valt=valt+vtrial1(ii,1)*rhor1(ii,1)
!end do
!Local potential energy of this band, valuer has been accumulated
!valt=0.5_dp*valt*gs_hamkq%ucvol/gs_hamkq%nfft
!end if
!write(6,*)' vtorho3 : eloc1_k=',two*valt
!ENDDEBUG

 deallocate(rhoaug1,vlocal,vlocal1)
 if(mkqmem==0) then
  deallocate(cgq_disk)
  if (psps%usepaw==1) then
   if (mcprjq_disk>0) then
    call cprj_free(cprjq_disk)
   end if
   deallocate(cprjq_disk)
  end if
 end if

 call status(0,dtfil%filstat,iexit,level,'after loops   ')

 call timab(124,2,tsec)

 if(mpi_enreg%paral_compil_kpt==1)then

  call timab(129,1,tsec)

  buffer_size=7+mbd2kpsp+mbdkpsp
  if (iscf_mod>0) buffer_size=buffer_size+cplex*dtset%nfft*nspden
  allocate(buffer1(buffer_size))

! Pack rhor1,edocc,eeig0,ek0,ek1,eloc0,enl0,enl1,eigen1,resid
  if (iscf_mod>0) then
   index1=cplex*dtset%nfft*nspden
   if (psps%usepaw==0) then
    buffer1(1:index1)=reshape(rhor1  ,(/index1/))
   else
    buffer1(1:index1)=reshape(rho1wfr,(/index1/))
   end if
  else
   index1=0
  end if
  buffer1(index1+1       )=edocc
  buffer1(index1+2       )=eeig0
  buffer1(index1+3       )=ek0
  buffer1(index1+4       )=ek1
  buffer1(index1+5       )=eloc0
  buffer1(index1+6       )=enl0
  buffer1(index1+7       )=enl1
  index1=index1+7
  bdtot_index=0
  bd2tot_index=0
  do isppol=1,nsppol
   do ikpt=1,nkpt_rbz
    nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
    buffer1(index1+1:index1+2*nband_k**2)=&
&    eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)
    buffer1(index1+2*nband_k**2+1:index1+2*nband_k**2+nband_k)=&
&    resid(bdtot_index+1:bdtot_index+nband_k)
    bdtot_index=bdtot_index+nband_k
    bd2tot_index=bd2tot_index+2*nband_k**2
    index1=index1+2*nband_k**2+nband_k
   end do
  end do
  if(index1/=buffer_size)buffer1(index1+1:buffer_size)=zero

! Build sum of everything
  call timab(48,1,tsec)
  write(message, '(a,i8,a)' ) &
&  ' vtorho3 : MPI_ALLREDUCE, buffer of size',&
&  8*buffer_size,' bytes'
  call wrtout(6,message,'COLL')

  call xsum_mpi(buffer1,buffer_size,spaceworld,ierr)
  call timab(48,2,tsec)

! Unpack the final result
  if(iscf_mod>0) then
   index1=cplex*dtset%nfft*nspden
   if (psps%usepaw==0) then
    rhor1(:,:)  =reshape(buffer1(1:index1),(/cplex*dtset%nfft,nspden/))
   else
    rho1wfr(:,:)=reshape(buffer1(1:index1),(/cplex*dtset%nfft,nspden/))
   end if
  else
   index1=0
  end if
  edocc            =buffer1(index1+1)
  eeig0            =buffer1(index1+2)
  ek0              =buffer1(index1+3)
  ek1              =buffer1(index1+4)
  eloc0            =buffer1(index1+5)
  enl0             =buffer1(index1+6)
  enl1             =buffer1(index1+7)
  index1=index1+7
  bdtot_index=0
  bd2tot_index=0
  do isppol=1,nsppol
   do ikpt=1,nkpt_rbz
    nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
    eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)=&
&    buffer1(index1+1:index1+2*nband_k**2)
    resid(bdtot_index+1:bdtot_index+nband_k)=&
&    buffer1(index1+2*nband_k**2+1:index1+2*nband_k**2+nband_k)
    bdtot_index=bdtot_index+nband_k
    bd2tot_index=bd2tot_index+2*nband_k**2
    index1=index1+2*nband_k**2+nband_k
   end do
  end do
  deallocate(buffer1)
  call timab(129,2,tsec)
 end if ! if kpt parallel

 call timab(127,1,tsec)

!If needed, compute rhog1, and symmetrizes the density
 if (iscf_mod > 0) then

  call status(0,dtfil%filstat,iexit,level,'call symrhg   ')
! In order to have the symrhg working in parallel on FFT coefficients, the size
! of irzzon1 and phnons1 should be set to nfftot. Therefore, nsym\=1 does not work.
  if (psps%usepaw==0) then
   call symrhg(cplex,densymop_rf,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,&
&              nspden,nsppol,nsym1,dtset%paral_kgb,phnons1,rhog1  ,rhor1  ,symaf1)
  else
   call symrhg(cplex,densymop_rf,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,&
&              nspden,nsppol,nsym1,dtset%paral_kgb,phnons1,rho1wfg,rho1wfr,symaf1)
  end if
! We now have both rho(r) and rho(G), symmetrized, and if nsppol=2
! we also have the spin-up density, symmetrized, in rhor1(:,2).

! DEBUG
! write(6,*)' vtorho3 : after symrhg '
! write(6,*)' nsym1=',nsym1
! do ifft=1,dtset%nfft,13
! if (psps%usepaw==0) write(6,*)ifft,rhog1(:,ifft)
! if (psps%usepaw==1) write(6,*)ifft,rho1wfg(:,ifft)
! end do
! stop
! ENDDEBUG

 end if

 deallocate(kg_k,kg1_k,kg_dum)
 if (dtset%berryopt == 4) then
  deallocate(grad_berry)
 end if

!Find largest residual over bands, k points, and spins
!except for nbdbuf highest bands
 ibdkpt=1
 residm=zero
 do isppol=1,nsppol
  do ikpt=1,nkpt_rbz
   nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
   nband_k=max(1,nband_k-dtset%nbdbuf)
   residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_k-1)))
   ibdkpt=ibdkpt+nband_k
  end do
 end do

 if (iscf_mod>0) then

! Norm-conserving: The density related to WFs is already the total density
  if (psps%usepaw==0) then
   if (optres==1) nvresid1=rhor1-nvresid1
   call timab(127,2,tsec)
  else

! PAW: Symetrize new 1st-order rhoij quantities
!      Compute and add the 1st-order compensation density to rho1wfr
!      to get the 1st-order total density
   call timab(127,2,tsec)
   call timab(565,1,tsec)
   call symrhoij3(psps%indlmn,indsy1,ipert,psps%lmnmax,natom,dimpaw1,nsym1,ntypat,&
&                 pawang,dtset%pawprtvol,pawrhoij1,symaf1,symrc1,dtset%typat)
   do iatom=1,dimpaw1
    deallocate(pawrhoij1(iatom)%rhoij_)
    pawrhoij1(iatom)%use_rhoij_=0
   end do
   call timab(565,2,tsec)
   call timab(566,1,tsec)
   call pawmknhat3(cplex,idir,ipert,0,mpi_enreg,natom,nfftf,pawfgr%ngfft,nhat1dim,dimpaw1,nspden,&
&       ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhat1,pawrhoij,pawrhoij1,pawtab,dtset%typat)
   call timab(566,2,tsec)
   call timab(567,1,tsec)
   call transgrid(cplex,mpi_enreg,nspden,+1,1,0,dtset%paral_kgb,pawfgr,rho1wfg,rhodum,rho1wfr,rhor1)
   rhor1(:,:)=rhor1(:,:)+nhat1(:,:);if (optres==1) nvresid1=rhor1-nvresid1
   call fourdp(cplex,rhog1,rhor1(:,1),-1,mpi_enreg,nfftf,pawfgr%ngfft,dtset%paral_kgb,0)
   deallocate(rho1wfr,rho1wfg)
   call timab(567,2,tsec)
  end if

! Compute square norm nres2 of density residual nvresid1
  if (optres==1) call sqnorm_v(1,mpi_enreg,nfftf,nres2,dtset%nspden,nvresid1)
 end if ! iscf>0

 if (iscf_mod<=0) call timab(127,2,tsec)

!Rotate labels of disk files when wf i/o is used
!Free memory needed by PAW
 if (mk1mem==0) then
  wfftmp=wffnow ; wffnow=wffnew ; wffnew=wfftmp
 end if

!Free memory needed by PAW
 if (psps%usepaw==1) then
  if (mk1mem/=0.and.usecprj==1) then
   call cprj_free(cprj1)
  end if
  deallocate(cprj1)
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i2,a)') ch10,' vtorho3 : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(121,2,tsec)

!DEBUG
!write(6,*)' vtorho3 : exit '
!ENDDEBUG
end subroutine vtorho3
!!***
