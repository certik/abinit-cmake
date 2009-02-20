!{\src2tex{textfont=tt}}
!!****f* ABINIT/vtorho
!! NAME
!! vtorho
!!
!! FUNCTION
!! This routine compute the new density from a fixed potential (vtrial)
!! but might also simply compute eigenvectors and eigenvalues.
!! The main part of it is a wf update over all k points.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MF, AR, MM, MT, FJ, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  afford=used to dimension susmat
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dbl_nnsclo=if 1, will double the value of dtset%nnsclo
!!  densymop_diel <type(dens_sym_operator_type)>=the density symmetrization
!!   operator for the dielectric matrix
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dielop= if positive, the dielectric matrix must be computed.
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  etotal=total energy (Ha) - only needed for tddft
!!  filapp= character string giving the root to form the name of the EIG file
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for the dielectric matrix
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!   (3x3 tensor) and grads wrt atomic coordinates (3*natom)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  irrzondiel(nfftdiel**(1-1/nsym),2,nspden/nsppol)=irreducible zone data for diel matrix
!!                                     nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  istep=index of the number of steps in the routine scfcv
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfftf,nkxc)=exchange-correlation kernel, needed only if nkxc/=0 .
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!                see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  npwdiel=size of the susmat array.
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  optforces=option for the computation of forces (0: no force;1: forces)
!!  optres=0: the new value of the density is computed in place of the input value
!!         1: only the density residual is computed ; the input density is kept
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!                                    nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  phnonsdiel(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases,
!!   for diel matr
!!                                     nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1ddiel(2,3*(2*mgfftdiel+1)*natom*usepaw)=one-dimensional structure factor information
!!                                             for the dielectric matrix
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive vectors
!!  shiftvector((mband+2)*nkpt)=UNDER DEVELOPMENT
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  ucvol=unit cell volume in bohr**3.
!!  wffnew,wffnow=unit numbers for wf disk files.
!!  val_min,val_max=UNDER DEVELOPMENT
!!  vtrial(nfftf,nspden)=INPUT potential Vtrial(r).
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmdiel(npwdiel,lmax_diel**2)= real spherical harmonics for each G and k point
!!                                 for the dielectric matrix
!!
!! OUTPUT
!!  compch_fft=-PAW only- compensation charge inside spheres computed over fine fft grid
!!  dphase(3) : dphase(idir) = accumulated change in the string-averaged
!!     Zak phase along the idir-th direction caused by the update of all
!!     the occupied Bloch states at all the k-points (only if finite electric field)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points.
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!  susmat(2,npwdiel*afford,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  === if optforces>0 ===
!!    grnl(3*natom)=stores grads of nonlocal energy wrt length scales
!!  ==== if optres==1
!!    nres2=square of the norm of the residual
!!    nvresid(nfftf,nspden)=density residual
!!  ==== if psps%usepaw==1
!!    nhat(nfftf,nspden*psps%usepaw)=compensation charge density on rectangular grid in real space
!!
!! SIDE EFFECTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!   At output contains updated wavefunctions coefficients;
!!    if nkpt>1, these are kept in a disk file.
!!  energies <type(energies_type)>=storage for energies computed here :
!!   | e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_kinetic=kinetic energy part of total energy
!!   | e_nonlocalpsp=nonlocal pseudopotential part of total energy
!!   | e_fermie=fermi energy (Hartree)
!!  gsc(2,mpw*nspinor*mband*mkmem*nsppol*usepaw)=<g|S|c> matrix elements (S=overrlap)
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k.
!!      (input if insulator - occopt<3 - ; output if metallic)
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  rhog(2,nfftf)=Fourier transform of total electron density
!!  rhor(nfftf,nspden)=total electron density (el/bohr**3)
!!  wvl <type(wvl_data)>=wavelets structures in case of wavelets basis.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      clsopn,cprj_alloc,cprj_free,ctocprj,fftpac,fourdp,hdr_io,hdr_io_netcdf
!!      hdr_skip,hdr_update,ini_wf_netcdf,leave_new,leave_test,mkffnl,mkkin
!!      mkkpg,mkrho,mpi_recv,mpi_send,newocc,pawmknhat,pawmkrhoij,ph1d3d
!!      prteigrs,prtrhomxmn,rdnpw,rwwf,sphereboundary,sqnorm_v,status
!!      suscep_stat,symrhg,symrhoij,tddft,testsusmat,timab,transgrid,vtowfk
!!      wffkg,wrtout,wvl_nl_gradient,wvl_vtorho,xallgather_mpi,xallgatherv_mpi
!!      xcomm_init,xdefineoff,xmaster_init,xmax_mpi,xme_init,xredxcart,xsum_mpi
!!      xsum_mpi_dp2d
!!
!! NOTES
!!  Be careful to the meaning of nfft (size of FFT grids):
!!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!   - In case of PAW calculations:
!!     Two FFT grids are used; one with nfft points (coarse grid) for
!!     the computation of wave functions ; one with nfftf points
!!     (fine grid) for the computation of total density.
!!
!!  The total electronic density (rhor,rhog) is divided into two terms:
!!   - The density related to WFs =Sum[Psi**2]
!!   - The compensation density (nhat) - only in PAW
!!
!!  The parallelisation needed for the electric field should be
!!  made an independent subroutine, so that this routine could be put
!!  back in the 21drive directory.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine vtorho(afford,atindx,atindx1,cg,compch_fft,cpus,dbl_nnsclo,&
&           densymop_diel,densymop_gs,dielop,dielstrt,dphase,dtefield,dtfil,dtset,&
&           eigen,energies,etotal,filapp,gbound_diel,&
&           gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
&           istep,kg,kg_diel,kxc,lmax_diel,mgfftdiel,mpi_enreg,&
&           mpsang,natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
&           npwarr,npwdiel,nres2,nspinor,ntypat,nvresid,occ,optforces,&
&           optres,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
&           pwind,pwind_alloc,pwnsfac,resid,residm,rhog,rhor,&
&           rmet,rprimd,shiftvector,susmat,symrec,&
&           ucvol,wffnew,wffnow,val_min,val_max,vtrial,wvl,xred,ylm,ylmdiel)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
#if defined HAVE_NETCDF
 use netcdf
#endif

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13ionetcdf
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_14occeig
 use interfaces_14wvl_wfs
 use interfaces_15common
 use interfaces_17suscep
 use interfaces_18seqpar, except_this_one => vtorho
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments -------------------------------
 integer, intent(in) :: afford,dbl_nnsclo,dielop,dielstrt,istep,lmax_diel,mgfftdiel
 integer, intent(in) :: mpsang,natom,nfftf,nfftdiel,nkxc,npwdiel
 integer, intent(inout) :: nspinor
 integer, intent(in) :: ntypat,optforces,optres,pwind_alloc
 real(dp), intent(in) :: cpus,etotal,gsqcut,ucvol,val_max,val_min
 real(dp), intent(out) :: compch_fft,nres2,residm
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(inout) :: dtset
 type(dens_sym_operator_type), intent(in) :: densymop_diel,densymop_gs
 type(efield_type), intent(inout) :: dtefield
 type(energies_type), intent(inout) :: energies
 type(hdr_type), intent(inout) :: hdr
 type(pawang_type), intent(in) :: pawang
 type(pawfgr_type), intent(in) :: pawfgr
 type(pseudopotential_type), intent(in) :: psps
 type(wffile_type), intent(inout) :: wffnew,wffnow
 type(wvl_data), intent(inout) :: wvl
 integer, intent(in) :: atindx(natom),atindx1(natom),gbound_diel(2*mgfftdiel+8,2)
 integer, intent(in) :: indsym(4,dtset%nsym,natom),irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
 integer, intent(in) :: irrzondiel(nfftdiel**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol),kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: kg_diel(3,npwdiel),nattyp(ntypat),ngfftdiel(18),npwarr(dtset%nkpt)
 integer, intent(in) :: pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
 real(dp), intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp), intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*psps%usepaw)
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
 real(dp), intent(in) :: phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc),rmet(3,3),rprimd(3,3),shiftvector((dtset%mband+2)*dtset%nkpt)
 real(dp), intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp), intent(inout) :: xred(3,natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm)
 real(dp), intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
 real(dp), intent(out) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol),dphase(3),grnl(3*natom)
 real(dp), intent(out) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp), intent(out) :: nvresid(nfftf,dtset%nspden),resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(out) :: susmat(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden)
 real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp), intent(inout) :: kxc(nfftf,nkxc),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
 character(len=fnlen), intent(in) :: filapp
 type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in)  :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
 integer,parameter :: level=7
 integer,save :: nwarning=0
 integer :: bandtot_glob,bantot,bdtot_index,choice,count,count1,counter
 integer :: dest,dimdij,dimffnl,enunit
 integer :: fform,formeig,i1,i1inv,i2,i2inv,i3,i3inv,ia,iatom,iband,ibdkpt
 integer :: ibg,ic,icg,icg1,icg2,ider,idir,idum,idum1
 integer :: ierr,iexit,ifft,ifftinv,ifor,ifor1,ii,ikg,ikg1,ikg2,ikpt
 integer :: ikptf,ikpt1f,ikpt1i,imagn
 integer :: ikpt_loc,ikpt1,ikxc,ilm,ilmn,index,index1,index2,iproc
 integer :: ipw,ir,iscf,isp,ispden,isppol,istwf_k,itypat,jkpt,jkpti
 integer :: jsppol,lmnmax1,lmnmax2,master,matblk,mbdkpsp
 integer :: mcg,mcgq,mcg_disk,me,me_distrb,mkgq,mu,muig
 integer :: mwarning,n1,n2,n3,n4,n5,n6,nband_eff
 integer :: nband_k,nbuf,nfftot,ngb,nkpg,nkpt1,nnsclo_now,npw_k,npw_k1
 integer :: npwin,nsp2,numb,nvloc,option,prtvol
 integer :: rdwr,source,spaceComm,tag,tim_mkrho,tim_rwwf,tobox,tosph,usecprj
 real(dp) :: arg,dphase_str,dummy,emax,min_occ,rdum,vxcavg_dum
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk
 type(wffile_type) :: wfftmp
 integer,allocatable :: dimcprj(:),kg_dum(:,:),kg_k(:,:)
 real(dp) :: adum(3),dec(4),dielar(7),dphase_k(3),ehart(3),kpoint(3),rhodum(1),tsec(2),ylmgr_dum(1)
 real(dp),allocatable :: EigMin(:,:),buffer(:,:),buffer1(:),cgq(:,:)
 real(dp),allocatable :: cg_disk(:,:),cgrkxc(:,:),cgrvtrial(:,:),cr(:,:),doccde(:)
 real(dp),allocatable :: dphasek(:,:),eig_dum(:),eig_k(:),ek_k(:),eknk(:)
 real(dp),allocatable :: enl_k(:),enlnk(:),ffnl(:,:,:,:),grnl_k(:,:), xcart(:,:)
 real(dp),allocatable :: grnlnk(:,:),kinpw(:),kpg_k(:,:),occ_dum(:),occ_k(:),ph3d(:,:,:)
 real(dp),allocatable :: pwnsfacq(:,:),resid_k(:),rhoaug(:,:,:,:),rhowfg(:,:),rhowfr(:,:)
 real(dp),allocatable :: vlocal(:,:,:,:),vlocal_tmp(:,:,:),wavef(:,:),ylm_k(:,:),zshift(:)
 type(cprj_type),allocatable :: cprj(:,:),cprj_tmp(:,:)
 logical :: computesusmat
!FB Local variables for band_fft parallelization
 integer:: blocksize,ndatarecv,npw_tot,oldspacecomm,usebandfft
 integer,allocatable :: recvcounts(:),sendcounts(:),sdispls(:),rdispls(:)
 integer,allocatable :: sendcountsloc(:),sdisplsloc(:),recvcountsloc(:),rdisplsloc(:)
 integer,allocatable :: kg_k_gather(:,:),kg_k_gather_all(:,:),rdispls_all(:),npw_per_proc(:)
 real(dp),allocatable :: ffnl_little(:,:,:,:),ffnl_little_gather(:,:,:,:),&
                         ph3d_little(:,:,:),ph3d_little_gather(:,:,:)
 real(dp),allocatable :: kinpw_gather(:),ffnl_gather(:,:,:,:),vlocal_allgather(:,:,:,:),&
                         ph3d_gather(:,:,:)
 real(dp), allocatable,save :: gbound_all(:,:,:)

#if defined MPI
             integer :: status1(MPI_STATUS_SIZE)
#endif

! *********************************************************************

 !DEBUG
!write(6,*)' vtorho : enter'
!write(6,*)' vtorho : enter cg2',cg(2,1:10)
!ENDDEBUG

!Keep track of total time spent in vtorho
 call timab(21,1,tsec)
 call timab(24,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' vtorho : enter '
  call wrtout(06,message,'COLL')
 end if

 if (mpi_enreg%mode_para=='b') then
  usebandfft=1
  if (istep<=1) then
   if (allocated(gbound_all)) deallocate(gbound_all)
   allocate(gbound_all(2*dtset%mgfft+8,2,dtset%nkpt))
  endif
 else
  usebandfft=0
 endif

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)

!PATCH vtorho // KPT & FFT me-->me_kpt
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
    &(mpi_enreg%paral_compil_fft==1)) then
    me_distrb = mpi_enreg%me_kpt
 else
    me_distrb = mpi_enreg%me
 end if

!Init master
 call xmaster_init(mpi_enreg,master)

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' vtorho :  BUG -',ch10,&
&  '  wrong values for nfft, nfftf !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Test optforces (to prevent memory overflow)
 if (optforces/=0.and.optforces/=1) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' vtorho :  BUG -',ch10,&
&  '  wrong value for optforces !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Debugging : print vtrial and rhor
! MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(prtvol==-level)then
  if (psps%usepaw==0) then
   n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
  else
   n1=pawfgr%ngfft(1) ; n2=pawfgr%ngfft(2) ; n3=pawfgr%ngfft(3)
  end if
  write(message,'(a)') '   ir              vtrial(ir)     rhor(ir) '
  call wrtout(06,message,'COLL')
  do ir=1,nfftf
!  if(ir<=11 .or. mod(ir,301)==0 )then
    i3=(ir-1)/n1/n2
    i2=(ir-1-i3*n1*n2)/n1
    i1=ir-1-i3*n1*n2-i2*n1
    write(message,'(i5,3i3,a,2es13.6)')ir,i1,i2,i3,' ',vtrial(ir,1),rhor(ir,1)
    call wrtout(06,message,'COLL')
    if(dtset%nspden>=2)then
     write(message,'(a,2es13.6)')'               ',vtrial(ir,2),rhor(ir,2)
     call wrtout(06,message,'COLL')
    end if
!  end if
  end do
 end if

 ! WVL - Branching with a separate vtorho procedure
 !       in wavelet. Should be merge in vtorho later on.
 if (dtset%usewvl == 1) then
  call wvl_vtorho(dtset, energies, istep, mpi_enreg, &
       & occ, wvl%projectors, psps, residm, rhor, vtrial, wvl%wfs)
  if (optforces == 1) then
     allocate(xcart(3, dtset%natom))
     call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
     call wvl_nl_gradient(dtset, grnl, mpi_enreg, occ, psps, rprimd, wvl, xcart)
     deallocate(xcart)
  end if
  return
 end if
 ! WVL - Following is done in plane waves.

 iscf=dtset%iscf
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)

 energies%e_eigenvalues = zero
 energies%e_kinetic     = zero
 energies%e_nonlocalpsp = zero
 grnl(:)=zero
 bdtot_index=0
 ibg=0;icg=0
 mbdkpsp=dtset%mband*dtset%nkpt*dtset%nsppol

 allocate(eknk(mbdkpsp),kg_k(3,dtset%mpw))
 allocate(EigMin(2,dtset%mband))
 allocate(grnlnk(3*natom,mbdkpsp*optforces))
 if (psps%usepaw==0) allocate(enlnk(mbdkpsp))

 eknk(:)=zero;if (optforces>0) grnlnk(:,:)=zero
 if (psps%usepaw==0) enlnk(:)=zero

!Initialize rhor if needed; store old rhor
 if(iscf>0 .or. iscf==-3) then
  if (optres==1) nvresid=rhor
  if (psps%usepaw==0) then
   rhor=zero
  else
   allocate(rhowfr(dtset%nfft,dtset%nspden),rhowfg(2,dtset%nfft))
   rhowfr(:,:)=zero
  end if
 end if

!Set max number of non-self-consistent loops nnsclo_now for use in vtowfk
 if(iscf<=0)then
  nnsclo_now=dtset%nstep
 else if(iscf>0)then
  if(dtset%nnsclo>0) then
   nnsclo_now=dtset%nnsclo
  else if(dtset%nnsclo<=0)then
   nnsclo_now=1
   if(istep<=2)nnsclo_now=2
  end if
  if(dbl_nnsclo==1)then
!DEBUG
!  write(6,*)' vtorho : use doubled nnsclo '
!ENDDEBUG
   nnsclo_now=nnsclo_now*2
  end if
 end if
 if(dtset%wfoptalg==2)nnsclo_now=40  ! UNDER DEVELOPMENT

 write(message, '(a,i3,a,i3,i2,i3)' ) ' vtorho : nnsclo_now=',nnsclo_now,&
& ', note that nnsclo,dbl_nnsclo,istep=',dtset%nnsclo,dbl_nnsclo,istep
 call wrtout(6,message,'COLL')

 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 nvloc=1;if(dtset%nspden==4)nvloc=4
 allocate(rhoaug(n4,n5,n6,nvloc),vlocal(n4,n5,n6,nvloc))

!Prepare wf files for reading if dtset%mkmem==0
 if (dtset%mkmem==0) then

! Close files, and then reopen them
! (this is supposedly helpful for use of networked workstations
! and also sets up for later addition of a checkpoint facility
! for restarting crashed jobs)
! clsopn automatically checks to see whether file is scratch
! file and if so, does not close and open it.
  call clsopn(wffnow)

! Read wffnow header
  call hdr_skip(wffnow,ierr)

! Define offsets, in case of MPI I/O
  formeig=0
  call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,nspinor,dtset%nsppol,dtset%nkpt)

  call clsopn(wffnew)

! Update the content of the header (evolving variables)
  bantot=hdr%bantot ; dummy=1.0d20
  ! WARNING! fermie is used before set.
  call hdr_update(bantot,dummy,energies%e_fermie,hdr,natom,&
& residm,rprimd,occ,pawrhoij,psps%usepaw,xred)

! Write the content of hdr to the new wf file
  rdwr=2 ; fform=2
  if (wffnew%accesswff /= 2) then
    call hdr_io(fform,hdr,rdwr,wffnew)
#if defined HAVE_NETCDF
  else if (wffnew%accesswff == 2) then
    call hdr_io_netcdf(fform,hdr,rdwr,wffnew)

    call ini_wf_netcdf(dtset%mpw,wffnew%unwff,0)
#endif
  end if


! Define offsets, in case of MPI I/O
  formeig=0
  call WffKg(wffnew,1)
  call xdefineOff(formeig,wffnew,mpi_enreg,dtset%nband,npwarr,nspinor,dtset%nsppol,dtset%nkpt)

  mcg_disk=dtset%mpw*nspinor*dtset%mband
  allocate(cg_disk(2,mcg_disk))

 end if

!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 allocate(gs_hamk%atindx(natom),gs_hamk%atindx1(natom))
 allocate(gs_hamk%gbound(2*dtset%mgfft+8,2)); gs_hamk%gbound(:,:)=0
 allocate(gs_hamk%indlmn(6,psps%lmnmax,ntypat))
 allocate(gs_hamk%nattyp(ntypat))
 allocate(gs_hamk%phkxred(2,natom))
 allocate(gs_hamk%ph1d(2,3*(2*dtset%mgfft+1)*natom))
 allocate(gs_hamk%pspso(ntypat))
 allocate(gs_hamk%xred(3,natom))

!Initialize most of the Hamiltonian
 gs_hamk%atindx(:)  =atindx(:)
 gs_hamk%atindx1(:) =atindx1(:)
 gs_hamk%gmet(:,:)  =gmet(:,:)
 gs_hamk%gprimd(:,:)=gprimd(:,:)
 gs_hamk%indlmn(:,:,:)=psps%indlmn(:,:,:)
 gs_hamk%lmnmax     =psps%lmnmax
 gs_hamk%mgfft      =dtset%mgfft
 gs_hamk%mpsang     =mpsang
 gs_hamk%mpssoang   =psps%mpssoang
 gs_hamk%natom      =natom
 gs_hamk%nattyp(:)  =nattyp(:)
 gs_hamk%nfft       =dtset%nfft
 gs_hamk%ngfft(:)   =dtset%ngfft(:)
 gs_hamk%nloalg(:)  =dtset%nloalg(:)
 gs_hamk%nspinor    =nspinor
 gs_hamk%ntypat     =ntypat
 gs_hamk%nvloc      =nvloc
 gs_hamk%n4         =n4
 gs_hamk%n5         =n5
 gs_hamk%n6         =n6
 gs_hamk%usepaw     =psps%usepaw
 gs_hamk%ph1d(:,:)  =ph1d(:,:)
 gs_hamk%pspso(:)   =psps%pspso(1:ntypat)
 gs_hamk%ucvol      =ucvol
 gs_hamk%useylm     =psps%useylm
 gs_hamk%xred(:,:)  =xred(:,:)

!Non-local factors:
! Norm-conserving: kleimann-Bylander energies
! PAW: Dij coefficients and overlap coefficients
 if (psps%usepaw==0) then
  gs_hamk%dimekb1=psps%dimekb
  gs_hamk%dimekb2=ntypat
  allocate(gs_hamk%ekb(psps%dimekb,ntypat,nspinor**2))
  allocate(gs_hamk%sij(0,0))
  gs_hamk%ekb(:,:,1)=psps%ekb(:,:)
  if (nspinor==2) then
   gs_hamk%ekb(:,:,2)=psps%ekb(:,:)
   gs_hamk%ekb(:,:,3:4)=zero
  end if
  usecprj=0
 else
  gs_hamk%dimekb1=psps%dimekb*paw_ij(1)%cplex_dij
  gs_hamk%dimekb2=natom
  allocate(gs_hamk%ekb(gs_hamk%dimekb1,gs_hamk%dimekb2,nspinor**2))
  allocate(gs_hamk%sij(gs_hamk%dimekb1,ntypat))
  allocate(dimcprj(natom));ia=0
  do itypat=1,ntypat
   if (paw_ij(1)%cplex_dij==1) then
    gs_hamk%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
   else
    do ilmn=1,pawtab(itypat)%lmn2_size
     gs_hamk%sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
     gs_hamk%sij(2*ilmn  ,itypat)=zero
    end do
   end if
   dimcprj(ia+1:ia+nattyp(itypat))=pawtab(itypat)%lmn_size
   ia=ia+nattyp(itypat)
  end do
  usecprj=1
  allocate(cprj(natom,nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*usecprj))
  if (usecprj==1) then
   if (dtset%mkmem/=0) then
    call cprj_alloc(cprj,0,dimcprj)
   else
    rewind dtfil%unpaw
    write(dtfil%unpaw) natom,nspinor
    write(dtfil%unpaw) dimcprj(1:natom)
   end if
  end if
 end if

!Electric field: allocate dphasek
 if (dtset%berryopt == 4) then
  allocate(dphasek(3,dtset%nkpt*dtset%nsppol))
  dphasek(:,:) = zero
  nkpt1 = dtefield%mkmem_max
 else
  nkpt1 = dtset%nkpt
 end if

 ikpt_loc = 0


!LOOP OVER SPINS
 do isppol=1,dtset%nsppol


  if (dtset%nsppol==2) then
   write(message,*)' ****  In vtorho for isppol=',isppol
   call wrtout(06,message,'COLL')
  end if

  if ((mpi_enreg%paral_compil_kpt == 0).or.(dtset%berryopt /= 4)) ikpt_loc = 0

! Rewind kpgsph data file if needed:
  if (dtset%mkmem==0) rewind dtfil%unkg
  if (dtset%mkmem==0.and.psps%useylm==1) rewind dtfil%unylm
  ikg=0

! Set up local potential vlocal with proper dimensioning, from vtrial
! Also take into account the spin.
  if(dtset%nspden/=4)then
   if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
    call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal,2)
   else
    allocate(cgrvtrial(dtset%nfft,dtset%nspden))
    call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
    call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
    deallocate(cgrvtrial)
   end if
  else
   allocate(vlocal_tmp(n4,n5,n6))
   if (psps%usepaw==0.or.pawfgr%usefinegrid==0) then
    do ispden=1,dtset%nspden
     call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal_tmp,2)
     vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
    end do
   else
    allocate(cgrvtrial(dtset%nfft,dtset%nspden))
    call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
    do ispden=1,dtset%nspden
     call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal_tmp,2)
     vlocal(:,:,:,ispden)=vlocal_tmp(:,:,:)
    end do
    deallocate(cgrvtrial)
   end if
   deallocate(vlocal_tmp)
  end if
  rhoaug(:,:,:,:)=zero

! PAW: retrieve Dij coefficients for this spin component
  if (psps%usepaw==1) then
   do ispden=1,nspinor**2
    isp=isppol;if (nspinor==2) isp=ispden
    do iatom=1,dtset%natom
     dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
     do ilmn=1,dimdij
      gs_hamk%ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
     end do
     if(dimdij+1<=gs_hamk%dimekb1) gs_hamk%ekb(dimdij+1:gs_hamk%dimekb1,iatom,ispden)=zero
    end do
   end do
  end if

  call timab(25,1,tsec)
! BIG FAT k POINT LOOP
! MVeithen: I had to modify the structure of this loop in order to implement
!           MPI // of the electric field

  ikpt = 0
  do while (ikpt_loc < nkpt1)

   if (dtset%berryopt /= 4) then
     ikpt_loc = ikpt_loc + 1
     ikpt = ikpt_loc
   else
     if (ikpt_loc < dtset%mkmem) ikpt = ikpt + 1
     if ((ikpt > dtset%nkpt).and.(ikpt_loc < dtset%mkmem)) exit
   end if

   dphase_k(:) = zero
   counter=100*ikpt+isppol
   call status(counter,dtfil%filstat,iexit,level,'loop ikpt     ')
   nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
   istwf_k=dtset%istwfk(ikpt)
   npw_k=npwarr(ikpt)
   mcgq = 1 ; mkgq = 1

   if(mpi_enreg%paral_compil_kpt==1)then

    if (dtset%berryopt /= 4) then

     if (mpi_enreg%parareel == 0) then
      if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&          -me_distrb))/=0) then
       eigen(1+bdtot_index : nband_k+bdtot_index) = zero
       resid(1+bdtot_index : nband_k+bdtot_index) = zero
       bdtot_index=bdtot_index+nband_k
!      Skip the rest of the k-point loop
       cycle
      end if
     else
      if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt) &
&                 /= mpi_enreg%me) then
       eigen(1+bdtot_index : nband_k+bdtot_index) = zero
       resid(1+bdtot_index : nband_k+bdtot_index) = zero
       bdtot_index=bdtot_index+nband_k
!      Skip the rest of the k-point loop
       cycle
      end if
     end if

    else      ! dtset%berryopt /= 4

     if ((minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) - &
&        mpi_enreg%me)) /= 0).and.(ikpt_loc <= dtset%mkmem)) then
      eigen(1+bdtot_index : nband_k+bdtot_index) = zero
      resid(1+bdtot_index : nband_k+bdtot_index) = zero
      bdtot_index = bdtot_index + nband_k
      cycle
     end if

     mcgq = dtset%mpw*nspinor*nband_k*dtefield%nneigh(ikpt)
     mkgq = 6*dtset%mpw
     ikg = dtefield%kgindex(ikpt)

    end if     ! dtset%berryopt /= 4

   end if ! parallel kpt

   if (dtset%berryopt == 4) ikpt_loc = ikpt_loc + 1
   allocate(cgq(2,mcgq),pwnsfacq(2,mkgq))

! In case of MPI // of a finite electric field calculation
! build the cgq array that stores the wavefunctions for the
! neighbours of ikpt, and the pwnsfacq array that stores the
! corresponding phase factors (in case of tnons)


#if defined MPI

   if (dtset%berryopt == 4) then

     ikptf = dtefield%i2fbz(ikpt)

     do idir = 1, 3

! skip idir values for which efield_dot(idir) = 0
       if (abs(dtefield%efield_dot(idir)) < tol12) cycle

       do ifor = 1, 2

         dtefield%sflag(:,ikpt + dtset%nkpt*(isppol - 1),ifor,idir) = 0

         ikpt1f = dtefield%ikpt_dk(ikptf,ifor,idir)
         ikpt1i = dtefield%indkk_f2ibz(ikpt1f,1)
         npw_k1 = npwarr(ikpt1i)
         count = npw_k1*nspinor*nband_k
         source = mpi_enreg%proc_distrb(ikpt1i,1,isppol)

         do jsppol = 1, dtset%nsppol
           do jkpt = 1, dtefield%fnkpt

             jkpti = dtefield%indkk_f2ibz(jkpt,1)

! recieve phase factors

             if ((jkpt == ikpt1f).and.(source /= mpi_enreg%me).and.&
&             (ikpt_loc <= dtset%mkmem).and.(jsppol == isppol)) then

               allocate(buffer(2,npw_k1))
               tag = jkpt + (jsppol - 1)*dtefield%fnkpt

               call MPI_RECV(buffer,2*npw_k1,MPI_DOUBLE_PRECISION,&
&                   source,tag,spaceComm,status1,ierr)
               ikg1 = dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
               pwnsfacq(:,ikg1+1:ikg1+npw_k1) = buffer(:,1:npw_k1)
               deallocate(buffer)

             end if

! send the phase factors to ALL the cpus that need it

             do dest = 1, mpi_enreg%nproc

               if ((minval(abs(mpi_enreg%proc_distrb(jkpti,1:dtset%mband,jsppol)&
&                - mpi_enreg%me)) == 0).and.&
&               (mpi_enreg%kptdstrbi(dest,ifor+2*(idir-1),&
&                 ikpt_loc + dtefield%mkmem_max*dtset%nsppol) == &
&                 jkpt + (jsppol - 1)*dtefield%fnkpt)) then

                 ikg1 = dtefield%fkgindex(jkpt)

                 if (((dest-1) /= mpi_enreg%me)) then

                   tag = jkpt + (jsppol - 1)*dtefield%fnkpt
                   count1 = npwarr(jkpti)
                   allocate(buffer(2,count1))
                   buffer(:,1:count1)  = pwnsfac(:,ikg1+1:ikg1+count1)
                   call MPI_SEND(buffer,2*count1,MPI_DOUBLE_PRECISION,&
&                     (dest-1),tag,spaceComm,status1,ierr)
                   deallocate(buffer)

                 else

                   ikg2 = &
&                    dtefield%cgqindex(3,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
                   pwnsfacq(:,ikg2 + 1:ikg2 + npw_k1) = &
&                    pwnsfac(:,ikg1 + 1:ikg1 + npw_k1)

                 end if

               end if

             end do          ! dest

           end do          ! jkpt
         end do            ! jsppol

         do jsppol = 1, dtset%nsppol
           do jkpt = 1, dtset%nkpt

! recieve WF

             if ((jkpt == ikpt1i).and.(source /= mpi_enreg%me).and.&
&             (ikpt_loc <= dtset%mkmem).and.(jsppol == isppol).and.&
&             (dtefield%cgqindex(1,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt) == &
&                ifor+2*(idir-1))) then

               allocate(buffer(2,count))
               tag = jkpt + (jsppol - 1)*dtset%nkpt

PRINT *, 'vor receive 2'

               call MPI_RECV(buffer,2*count,MPI_DOUBLE_PRECISION,&
&                   source,tag,spaceComm,status1,ierr)
               icg1 = dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
               cgq(:,icg1+1:icg1+count) = buffer(:,1:count)
               deallocate(buffer)

             end if

! send the WF to ALL the cpus that need it

             do dest = 1, mpi_enreg%nproc

               if ((minval(abs(mpi_enreg%proc_distrb(jkpt,1:dtset%mband,jsppol)&
&                - mpi_enreg%me)) == 0).and.&
&               (mpi_enreg%kptdstrbi(dest,ifor+2*(idir-1),ikpt_loc) == &
&                 jkpt + (jsppol - 1)*dtset%nkpt)) then

                 icg1 = dtefield%cgindex(jkpt,jsppol)

                 if (((dest-1) /= mpi_enreg%me)) then

                   tag = jkpt + (jsppol - 1)*dtset%nkpt
                   count1 = npwarr(jkpt)*nband_k*nspinor
                   allocate(buffer(2,count1))
                   buffer(:,1:count1)  = cg(:,icg1+1:icg1+count1)

PRINT *, 'vor send 2'
                   call MPI_SEND(buffer,2*count1,MPI_DOUBLE_PRECISION,&
&                     (dest-1),tag,spaceComm,status1,ierr)
                   deallocate(buffer)

                 else

                   icg2 = &
&                    dtefield%cgqindex(2,ifor+2*(idir-1),ikpt+(isppol-1)*dtset%nkpt)
                   cgq(:,icg2 + 1:icg2 + count) = &
&                    cg(:,icg1 + 1:icg1 + count)

                 end if

               end if

             end do          ! dest

           end do          ! jkpt
         end do            ! jsppol

       end do         ! ifor
     end do         !idir

     if (ikpt_loc > dtset%mkmem) then
       deallocate(cgq,pwnsfacq)
       cycle
     end if

   end if      !  berryopt

#endif


!  Continue to initialize the Hamiltonian
   gs_hamk%istwf_k    =istwf_k
   gs_hamk%npw        =npw_k

   allocate(eig_k(nband_k),ek_k(nband_k))
   allocate(occ_k(nband_k),resid_k(nband_k))
   allocate(ylm_k(npw_k,mpsang*mpsang*psps%useylm))
   allocate(zshift(nband_k))
   allocate(grnl_k(3*natom,nband_k*optforces))
   if (psps%usepaw==0) allocate(enl_k(nband_k))

   eig_k(:)=zero
   ek_k(:)=zero
   if (optforces>0) grnl_k(:,:)=zero
   if (psps%usepaw==0) enl_k(:)=zero
   kpoint(:)=dtset%kptns(:,ikpt)
   gs_hamk%kpoint(:)=dtset%kptns(:,ikpt)
   occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
   resid_k(:)=zero
   zshift(:)=dtset%eshift

   if (dtset%mkmem==0) then
!FB For band_fft parallelization compatibility. TO BE CLEANED
    allocate(kg_k_gather(3,usebandfft))
    ndatarecv=0

!   Read (k+G) basis sphere data (same for each spin)
    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,dtfil%unkg)
!   Read k+g data
    read (dtfil%unkg) kg_k(1:3,1:npw_k)
    call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
!   Eventually read spherical harmonics
    if (psps%useylm==1) then
     read(dtfil%unylm)
     read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
    end if

!   Read the wavefunction block for ikpt,isppol
    call status(counter,dtfil%filstat,iexit,level,'read wfs      ')
    tim_rwwf=1
    allocate(eig_dum(dtset%mband),kg_dum(3,0),occ_dum(dtset%mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,dtset%mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&    npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
    deallocate(eig_dum,kg_dum,occ_dum)

   else

    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

!FB Transpose the kg_k array. TO BE CLEANNED
    if(mpi_enreg%mode_para=='b') then
     spaceComm=mpi_enreg%comm_band
     !blocksize=mpi_enreg%nproc_band

     allocate(sdispls       (mpi_enreg%nproc_band))
     allocate(sdisplsloc    (mpi_enreg%nproc_band))
     allocate(sendcounts    (mpi_enreg%nproc_band))
     allocate(sendcountsloc (mpi_enreg%nproc_band))
     allocate(rdispls       (mpi_enreg%nproc_band))
     allocate(rdisplsloc    (mpi_enreg%nproc_band))
     allocate(recvcounts    (mpi_enreg%nproc_band))
     allocate(recvcountsloc (mpi_enreg%nproc_band))

     call xallgather_mpi(npw_k,recvcounts,spaceComm,ierr)
     rdispls(1)=0
     do iproc=2,mpi_enreg%nproc_band
      rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
     end do
     ndatarecv=rdispls(mpi_enreg%nproc_band)+recvcounts(mpi_enreg%nproc_band)

     allocate(kg_k_gather(3,ndatarecv*usebandfft))
     recvcountsloc(:)=recvcounts(:)*3
     rdisplsloc(:)=rdispls(:)*3
     call xallgatherv_mpi(kg_k,3*npw_k,kg_k_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)

     if (istep<=1) then
      oldspacecomm=mpi_enreg%comm_fft
      allocate(npw_per_proc(mpi_enreg%nproc_fft),rdispls_all(mpi_enreg%nproc_fft))
      call xallgather_mpi(ndatarecv,npw_per_proc,oldspacecomm,ierr)
      rdispls_all(1)=0
      do iproc=2,mpi_enreg%nproc_fft
       rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
      end do
      npw_tot=rdispls_all(mpi_enreg%nproc_fft)+npw_per_proc(mpi_enreg%nproc_fft)
      allocate(kg_k_gather_all(3,npw_tot))
      call xallgatherv_mpi&
&          (kg_k_gather,3*ndatarecv,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,oldspacecomm,ierr)
      call sphereboundary(gs_hamk%gbound,istwf_k,kg_k_gather_all,dtset%mgfft,npw_tot)
      deallocate(kg_k_gather_all,npw_per_proc,rdispls_all)
      gbound_all(:,:,ikpt)=gs_hamk%gbound(:,:)
     else
      gs_hamk%gbound(:,:)=gbound_all(:,:,ikpt)
     endif
    else
     allocate(kg_k_gather(3,dtset%mpw*usebandfft))
     call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
     ndatarecv=0
    endif

    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
      ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     end do
    end if

!  End if for choice governed by dtset%mkmem
   end if

!  Set up remaining of the Hamiltonian

!  Compute (1/2) (2 Pi)**2 (k+G)**2:
   call status(0,dtfil%filstat,iexit,level,'call mkkin    ')
   allocate(kinpw(npw_k))
   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg_k,kinpw,kpoint,npw_k)

!  Allocate the arrays phkxred and ph3d, compute phkxred
!  and eventually ph3d.
   do ia=1,natom
    iatom=atindx(ia)
    arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
    gs_hamk%phkxred(1,iatom)=cos(arg)
    gs_hamk%phkxred(2,iatom)=sin(arg)
!   DEBUG
!    write(6, '(a,i4,2es16.6)' )&
!  &  'vtorho : iatom, phkxred',iatom,phkxred(1,iatom),phkxred(2,iatom)
!   ENDDEBUG
   end do
   if(dtset%nloalg(1)<=0)then
!   Here, only the allocation, not the precomputation.
    matblk=dtset%nloalg(4)
    allocate(ph3d(2,npw_k,matblk))
   else
!   Here, allocation as well as precomputation
    matblk=natom
    allocate(ph3d(2,npw_k,matblk))
    call ph1d3d(1,natom,kg_k,kpoint,matblk,natom,npw_k,n1,n2,n3,&
&               gs_hamk%phkxred,ph1d,ph3d)
   end if
   gs_hamk%matblk=matblk

!  Compute (k+G) vectors (only if useylm=1)
   nkpg=3*optforces*dtset%nloalg(5);allocate(kpg_k(npw_k,nkpg))
   if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!  Compute nonlocal form factors ffnl at all (k+G):
   call status(0,dtfil%filstat,iexit,level,'call mkffnl   ')

   ider=0;idir=0;dimffnl=1
   allocate(ffnl(npw_k,dimffnl,psps%lmnmax,ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&   gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&   npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&   psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
!DEBUG
!write(6,*) ffnl(:,1,1,1); write(6,*) kg_k(:,:)
!ENDEBUG

!FB Transpose the ffnl, kinpw and ph3d arrays. TO BE CLEANED
   allocate(ffnl_gather(ndatarecv,dimffnl,psps%lmnmax,ntypat*usebandfft))
   allocate(kinpw_gather(ndatarecv*usebandfft))
   allocate(ph3d_gather(2,ndatarecv,matblk*usebandfft))
   if (mpi_enreg%mode_para=='b') then
    allocate(ffnl_little(dimffnl,psps%lmnmax,ntypat,npw_k))
    allocate(ffnl_little_gather(dimffnl,psps%lmnmax,ntypat,ndatarecv))
    do ipw=1,npw_k
     ffnl_little(:,:,:,ipw)=ffnl(ipw,:,:,:)
    end do
    recvcountsloc(:)=recvcounts(:)*dimffnl*psps%lmnmax*ntypat
    rdisplsloc(:)=rdispls(:)*dimffnl*psps%lmnmax*ntypat
    call xallgatherv_mpi(ffnl_little,npw_k*dimffnl*psps%lmnmax*ntypat,ffnl_little_gather,&
&        recvcountsloc(:),rdisplsloc,spaceComm,ierr)
    do ipw=1,ndatarecv
     ffnl_gather(ipw,:,:,:)=ffnl_little_gather(:,:,:,ipw)
    end do
    deallocate(ffnl_little,ffnl_little_gather)

    recvcountsloc(:)=recvcounts(:)
    rdisplsloc(:)=rdispls(:)
    call xallgatherv_mpi(kinpw,npw_k,kinpw_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)

    allocate(ph3d_little(2,matblk,npw_k),ph3d_little_gather(2,matblk,ndatarecv))
    recvcountsloc(:)=recvcounts(:)*2*matblk
    rdisplsloc(:)=rdispls(:)*2*matblk
    do ipw=1,npw_k
     ph3d_little(:,:,ipw)=ph3d(:,ipw,:)
    end do
    call xallgatherv_mpi(ph3d_little,npw_k*2*matblk,ph3d_little_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)
    do ipw=1,ndatarecv
     ph3d_gather(:,ipw,:)=ph3d_little_gather(:,:,ipw)
    end do
    deallocate(ph3d_little_gather,ph3d_little)

    deallocate(sendcounts,recvcounts,sdispls,rdispls)
    deallocate(sendcountsloc,sdisplsloc)
    deallocate(recvcountsloc,rdisplsloc)
   endif

   call status(counter,dtfil%filstat,iexit,level,'call vtowfk   ')

!  Compute the eigenvalues, wavefunction, residuals,
!  contributions to kinetic energy, nonlocal energy, forces,
!  and update of rhor to this k-point and this spin polarization.
   if(dtset%mkmem/=0)then
    mcg=dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
    call vtowfk(cg,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,&
&    dtset,eig_k,ek_k,enl_k,ffnl,ffnl_gather,grnl_k,gs_hamk,&
&    ibg,icg,ikg,ikpt,iscf,isppol,kg_k,kg_k_gather,kinpw,kinpw_gather,kpg_k,&
&    psps%lmnmax,matblk,mcg,mcgq,dtset%mgfft,mkgq,dtset%mkmem,&
&    mpi_enreg,psps%mpsang,psps%mpssoang,&
&    dtset%mpw,natom,nband_k,ndatarecv,nkpg,dtset%nkpt,nnsclo_now,npw_k,npwarr,nspinor,ntypat,nvloc,n4,n5,n6,&
&    occ_k,optforces,ph3d,ph3d_gather,prtvol,psps,pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,&
&    rhoaug,usebandfft,usecprj,vlocal,dtset%wtk(ikpt),zshift)
   else if(dtset%mkmem==0)then
    mcg=dtset%mpw*nspinor*dtset%mband
    call vtowfk(cg_disk,cgq,cprj,cpus,dimcprj,dimffnl,dphase_k,dtefield,dtfil,&
&    dtset,eig_k,ek_k,enl_k,ffnl,ffnl_gather,grnl_k,gs_hamk,&
&    ibg,icg,ikg,ikpt,iscf,isppol,kg_k,kg_k_gather,kinpw,kinpw_gather,kpg_k,&
&    psps%lmnmax,matblk,mcg,mcgq,dtset%mgfft,mkgq,dtset%mkmem,&
&    mpi_enreg,psps%mpsang,psps%mpssoang,&
&    dtset%mpw,natom,nband_k,ndatarecv,nkpg,dtset%nkpt,nnsclo_now,npw_k,npwarr,nspinor,ntypat,nvloc,n4,n5,n6,&
&    occ_k,optforces,ph3d,ph3d_gather,prtvol,psps,pwind,pwind_alloc,pwnsfac,pwnsfacq,resid_k,&
&    rhoaug,usebandfft,usecprj,vlocal,dtset%wtk(ikpt),zshift)
   end if

   call status(counter,dtfil%filstat,iexit,level,'after vtowfk  ')
!DEBUG
!write(6,*) 'eig_k'; write(6,*) eig_k
!ENDDEBUG
   deallocate(ffnl,kinpw,kpg_k,ph3d,cgq,pwnsfacq)
!FB Deallocate arrays used for band_fft parallelization
   deallocate(ffnl_gather,kinpw_gather,ph3d_gather,kg_k_gather)

!  electric field
   if (dtset%berryopt == 4) then

    dphasek(:,ikpt + (isppol - 1)*dtset%nkpt) = dphase_k(:)

!   The overlap matrices for all first neighbours of ikpt
!   are no more up to date
    do idir = 1, 3
     do ifor = 1, 2
      ikpt1 = dtefield%ikpt_dk(dtefield%i2fbz(ikpt),ifor,idir)
      ikpt1 = dtefield%indkk_f2ibz(ikpt1,1)
      ifor1 = -1*ifor + 3   ! ifor = 1 -> ifor1 = 2 & ifor = 2 -> ifor1 = 1
      dtefield%sflag(:,ikpt1+(isppol-1)*dtset%nkpt,ifor1,idir) = 0
     end do
    end do

   end if  ! berryopt

!  Write new wavefunctions to disk if needed
!  (header records were written earlier in loopcv)
   if (dtset%mkmem==0) then
    tim_rwwf=1
    call rwwf(cg_disk,eig_k,0,0,0,ikpt,isppol,kg_k,dtset%mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&    npw_k,nspinor,occ_k,2,1,tim_rwwf,wffnew)
   end if

!  Save eigenvalues (hartree), residuals (hartree**2)
   eigen(1+bdtot_index : nband_k+bdtot_index) = eig_k(:)
   eknk (1+bdtot_index : nband_k+bdtot_index) = ek_k (:)
   resid(1+bdtot_index : nband_k+bdtot_index) = resid_k(:)
   if (optforces>0) grnlnk(:,1+bdtot_index : nband_k+bdtot_index) = grnl_k(:,:)
   if (psps%usepaw==0) enlnk(1+bdtot_index : nband_k+bdtot_index) = enl_k(:)

   if(iscf>0 .or. iscf==-3)then
!   Accumulate sum over k points for band, nonlocal and kinetic energies,
!   also accumulate gradients of Enonlocal:
    do iband=1,nband_k
     if (abs(occ_k(iband))>tol8) then
      energies%e_kinetic = energies%e_kinetic + &
                         & dtset%wtk(ikpt)*occ_k(iband)*ek_k(iband)
      energies%e_eigenvalues = energies%e_eigenvalues + &
                             & dtset%wtk(ikpt)*occ_k(iband)*eig_k(iband)
      if (optforces>0) grnl(:)=grnl(:)+dtset%wtk(ikpt)/ucvol*occ_k(iband)*grnl_k(:,iband)
      if (psps%usepaw==0) then
       energies%e_nonlocalpsp = energies%e_nonlocalpsp + &
                              & dtset%wtk(ikpt)*occ_k(iband)*enl_k(iband)
      end if
     end if
    end do
   end if
   deallocate(eig_k,ek_k,grnl_k,occ_k,resid_k,ylm_k,zshift)
   if (psps%usepaw==0) deallocate(enl_k)

!  Keep track of total number of bands (all k points so far, even for
!     k points not treated by me)
   bdtot_index=bdtot_index+nband_k

!  Also shift array memory if dtset%mkmem/=0
   if (dtset%mkmem/=0) then
    ibg=ibg+nspinor*nband_k
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   end if

! End big k point loop
  end do

  call status(counter,dtfil%filstat,iexit,level,'after k loop  ')

  call timab(25,2,tsec)


  if (dtset%occopt<3 .and. mpi_enreg%mode_para=='b' .and. dtset%wfoptalg==4) then
   spaceComm=mpi_enreg%comm_band !Sum the contributions of the bands
   call xsum_mpi(rhoaug,spaceComm,ierr)
   spaceComm=mpi_enreg%comm_fft
   call xsum_mpi(rhoaug,spaceComm,ierr)
  end if

! Transfer density on augmented fft grid to normal fft grid in real space
! Also take into account the spin.
  if(iscf>0)then
   if( mpi_enreg%paral_compil_kpt==0 .or.       &
&      mpi_enreg%paralbd <= 1        .or.       &
&     (mpi_enreg%paralbd >= 1 .and. mpi_enreg%me_group==0)) then
    if (psps%usepaw==0) then
     call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug(:,:,:,1),1)
     if(dtset%nspden==4)then
      do imagn=2,4
       call fftpac(imagn,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug(:,:,:,imagn),1)
      enddo
     end if
    else
     call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhowfr,rhoaug(:,:,:,1),1)
     if(dtset%nspden==4)then
      do imagn=2,4
       call fftpac(imagn,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhowfr,rhoaug(:,:,:,imagn),1)
      enddo
     end if
    end if
   end if
  end if

!End loop over spins
 end do

 call status(counter,dtfil%filstat,iexit,level,'after spinloop')

 if(mpi_enreg%paral_compil_kpt==1)then
  call timab(66,1,tsec)
  if (mpi_enreg%parareel == 0) call leave_test(mpi_enreg)
  write(message,*) 'vtorho: loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
  call timab(66,2,tsec)
 end if


! electric field: compute string-averaged change in Zak phase
! along each direction, store it in dphase(idir)

!ji: it is not convenient to do this anymore. Remove. Set dphase(idir)=0.0_dp.
!    eventually, dphase(idir) will have to go...

 if (dtset%berryopt == 4)  dphase(:) = 0.0_dp


! In case of MPI // of a finite field calculation, send dphasek to all cpus
 if ((mpi_enreg%paral_compil_kpt == 1).and.(dtset%berryopt == 4)) then
   call xsum_mpi_dp2d(dphasek,spaceComm,ierr)
 end if

 if (dtset%berryopt == 4) deallocate(dphasek) !! by MM
 deallocate(gs_hamk%atindx,gs_hamk%atindx1)
 deallocate(gs_hamk%ekb)
 deallocate(gs_hamk%sij)
 deallocate(gs_hamk%gbound)
 deallocate(gs_hamk%indlmn)
 deallocate(gs_hamk%nattyp)
 deallocate(gs_hamk%phkxred)
 deallocate(gs_hamk%pspso)
 deallocate(gs_hamk%ph1d)
 deallocate(gs_hamk%xred)

 deallocate(rhoaug,vlocal)
 if(dtset%mkmem==0)deallocate(cg_disk)

 call status(0,dtfil%filstat,iexit,level,'after loops   ')

 allocate(doccde(dtset%mband*dtset%nkpt*dtset%nsppol))
 doccde(:)=zero !MF initialize

 call timab(24,2,tsec)


!Treat now varying occupation numbers, in the self-consistent case
 if(dtset%occopt>=3 .and. dtset%occopt <=7 .and. (iscf>0.or.iscf==-3)) then

! Parallel case
  if( mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then

   call timab(29,1,tsec)

!  If needed, exchange the values of eigen,resid,eknk,enlnk,grnlnk
   allocate(buffer1((4+3*natom*optforces-psps%usepaw)*mbdkpsp))
!  Pack eigen,resid,eknk,enlnk,grnlnk in buffer1
   buffer1(1          :  mbdkpsp)=eigen(:)
   buffer1(1+  mbdkpsp:2*mbdkpsp)=resid(:)
   buffer1(1+2*mbdkpsp:3*mbdkpsp)=eknk(:)
   index1=3*mbdkpsp
   if (psps%usepaw==0) then
    buffer1(index1+1:index1+mbdkpsp)=enlnk(:)
    index1=index1+mbdkpsp
   end if
   if (optforces>0) then
    buffer1(index1+1:index1+3*natom*mbdkpsp)=reshape(grnlnk, (/(3*natom)*mbdkpsp/) )
   end if
!  Build sum of everything
   if (mpi_enreg%paral_compil_fft==1) then
           spaceComm=mpi_enreg%fft_master_comm
   end if
   call timab(48,1,tsec)

   !PATCH vtorho // KPT & FFT sum fft_master_comm --> comm_kpt
   if(mpi_enreg%mode_para/='b') then
      call xsum_mpi(buffer1,spaceComm ,ierr)
   else
      if ((mpi_enreg%paral_compil_kpt==1) .and. &
           &(mpi_enreg%paral_compil_fft==1)) then
         call xsum_mpi(buffer1,mpi_enreg%comm_kpt ,ierr)
      end if
   end if

   call timab(48,2,tsec)
!  Unpack eigen,resid,eknk,enlnk,grnlnk in buffer1
   eigen(:) =buffer1(1          :  mbdkpsp)
   resid(:) =buffer1(1+  mbdkpsp:2*mbdkpsp)
   eknk(:)  =buffer1(1+2*mbdkpsp:3*mbdkpsp)
   index1=3*mbdkpsp
   if (psps%usepaw==0) then
    enlnk(:) =buffer1(index1+1:index1+mbdkpsp)
    index1=index1+mbdkpsp
   end if
   if (optforces>0) then
    grnlnk(:,:)=reshape(buffer1(index1+1:index1+3*natom*mbdkpsp),&
&                      (/ 3*natom , mbdkpsp /) )
   end if
   deallocate(buffer1)
   call timab(29,2,tsec)

  end if ! parallel

  call timab(27,1,tsec)

! Compute the new occupation numbers from eigen
  call status(0,dtfil%filstat,iexit,level,'call newocc   ')
  call newocc(doccde,eigen,energies%entropy,energies%e_fermie,dtset%fixmom,&
&  dtset%mband,dtset%nband,dtset%nelect,dtset%nkpt,nspinor,&
&  dtset%nsppol,occ,dtset%occopt,prtvol,dtset%stmbias,dtset%tphysel,dtset%tsmear,dtset%wtk)

! Compute eeig, ek,enl and grnl from the new occ, and the shared eknk,enlnk,grnlnk
  energies%e_eigenvalues = zero
  energies%e_kinetic     = zero
  energies%e_nonlocalpsp = zero
  if (psps%usepaw==0) energies%e_nonlocalpsp = zero
  if (optforces>0) grnl(:)=zero
  bdtot_index=1
  do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
    nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
    do iband=1,nband_k
     if (abs(occ(bdtot_index))>tol8) then
      energies%e_eigenvalues = energies%e_eigenvalues + &
                             & dtset%wtk(ikpt)*occ(bdtot_index)*eigen(bdtot_index)
      energies%e_kinetic = energies%e_kinetic + &
                         & dtset%wtk(ikpt)*occ(bdtot_index)*eknk(bdtot_index)
      if (optforces>0) grnl(:)=grnl(:)+dtset%wtk(ikpt)/ucvol*occ(bdtot_index)*grnlnk(:,bdtot_index)
      if (psps%usepaw==0) energies%e_nonlocalpsp = energies%e_nonlocalpsp + &
                              & dtset%wtk(ikpt)*occ(bdtot_index)*enlnk(bdtot_index)
     end if
     bdtot_index=bdtot_index+1
    end do
   end do
  end do

  call status(0,dtfil%filstat,iexit,level,'call mkrho    ')
  tim_mkrho=2
  if (psps%usepaw==0) then
   call mkrho(cg,densymop_gs,dtset,irrzon,kg,mpi_enreg,npwarr,nspinor,occ,phnons,&
&             rhog  ,rhor  ,tim_mkrho,ucvol,dtfil%unkg,wffnew,wvl%wfs)
  else
   call mkrho(cg,densymop_gs,dtset,irrzon,kg,mpi_enreg,npwarr,nspinor,occ,phnons,&
&             rhowfg,rhowfr,tim_mkrho,ucvol,dtfil%unkg,wffnew,wvl%wfs)
  end if
  call timab(27,2,tsec)

!Treat fixed occupation numbers or non-self-consistent case
 else

  if(mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then

   call timab(29,1,tsec)

   nbuf=2*mbdkpsp+dtset%nfft*dtset%nspden+3-psps%usepaw+3*natom*optforces
   if(iscf==-1 .or. iscf==-2)nbuf=2*mbdkpsp
   allocate(buffer1(nbuf))
!  Pack eigen,resid,rho[wf]r,grnl,enl,ek
   buffer1(1:mbdkpsp)=eigen(:)
   buffer1(1+mbdkpsp:2*mbdkpsp)=resid(:)
   index1=2*mbdkpsp
   if(iscf/=-1 .and. iscf/=-2)then
    if (psps%usepaw==0) then
     buffer1(index1+1:index1+dtset%nfft*dtset%nspden)=reshape(rhor  ,&
&                                        (/dtset%nfft*dtset%nspden/))
    else
     buffer1(index1+1:index1+dtset%nfft*dtset%nspden)=reshape(rhowfr,&
&                                        (/dtset%nfft*dtset%nspden/))
    end if
    index1=index1+dtset%nfft*dtset%nspden
    buffer1(index1+1) = energies%e_kinetic
    buffer1(index1+2) = energies%e_eigenvalues
    if (psps%usepaw==0) buffer1(index1+3) = energies%e_nonlocalpsp
    index1=index1+3-psps%usepaw
    if (optforces>0) buffer1(index1+1:index1+3*natom)=grnl(1:3*natom)
   end if
!  Build sum of everything
   call timab(48,1,tsec)
   if (mpi_enreg%paral_compil_fft==1) then
           spaceComm=mpi_enreg%fft_master_comm
   end if

   !PATCH vtorho // KPT & FFT sum fft_master_comm --> comm_kpt
   if(mpi_enreg%mode_para/='b') then
      call xsum_mpi(buffer1,nbuf,spaceComm ,ierr)
   else
      if ((mpi_enreg%paral_compil_kpt==1) .and. &
           &(mpi_enreg%paral_compil_fft==1)) then
         call xsum_mpi(buffer1,nbuf,mpi_enreg%comm_kpt ,ierr)
      end if
   end if
   call timab(48,2,tsec)
!  Unpack the final result
   eigen(:)=buffer1(1:mbdkpsp)
   resid(:)=buffer1(1+mbdkpsp:2*mbdkpsp)
   index1=2*mbdkpsp
   if(iscf/=-1 .and. iscf/=-2)then
    if (psps%usepaw==0) then
     ii=1
     do ispden=1,dtset%nspden
      do ifft=1,dtset%nfft
       rhor(ifft,ispden)=buffer1(index1+ii)
       ii=ii+1
      enddo
     enddo
    else
     ii=1
     do ispden=1,dtset%nspden
      do ifft=1,dtset%nfft
       rhowfr(ifft,ispden)=buffer1(index1+ii)
       ii=ii+1
      enddo
     enddo
    end if
    index1=index1+dtset%nfft*dtset%nspden
    energies%e_kinetic = buffer1(index1+1)
    energies%e_eigenvalues = buffer1(index1+2)
    if (psps%usepaw==0) energies%e_nonlocalpsp = buffer1(index1+3)
    index1=index1+3-psps%usepaw
    if (optforces>0) grnl(1:3*natom)=buffer1(index1+1:index1+3*natom)
   end if
   deallocate(buffer1)
   call timab(29,2,tsec)

  end if ! parallel

! Compute the highest occupied eigenenergy (might be parallelized)
  if(iscf/=-1 .and. iscf/=-2)then
   energies%e_fermie = -huge(one)
   bdtot_index=1
   do isppol=1,dtset%nsppol
    do ikpt=1,dtset%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     do iband=1,nband_k
      if(abs(occ(bdtot_index))>tol8 .and. eigen(bdtot_index)>energies%e_fermie+tol10) then
        energies%e_fermie=eigen(bdtot_index)
      end if
      bdtot_index=bdtot_index+1
     end do
    end do
   end do
   if(mpi_enreg%mode_para/='b') then
    call xmax_mpi(energies%e_fermie,emax,spaceComm ,ierr)
    energies%e_fermie=emax
   else
    if ((mpi_enreg%paral_compil_kpt==1) .and. &
&       (mpi_enreg%paral_compil_fft==1)) then
     call xmax_mpi(energies%e_fermie,emax,mpi_enreg%comm_kpt ,ierr)
     energies%e_fermie=emax
    end if
   end if

  end if

! If needed, compute rhog, and symmetrizes the density
  if (iscf > 0) then

!  energies%e_fermie=zero  ! Actually, should determine the maximum of the valence band XG20020802

   call timab(70,1,tsec)

   call status(0,dtfil%filstat,iexit,level,'compute rhog  ')
   nfftot=dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3)
   if (psps%usepaw==0) then
    call symrhg(1,densymop_gs,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,&
&               dtset%nsppol,dtset%nsym,dtset%paral_kgb,phnons,rhog  ,rhor  ,dtset%symafm)
   else
    call symrhg(1,densymop_gs,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,&
&               dtset%nsppol,dtset%nsym,dtset%paral_kgb,phnons,rhowfg,rhowfr,dtset%symafm)
   end if
!  We now have both rho(r) and rho(G), symmetrized, and if dtset%nsppol=2
!  we also have the spin-up density, symmetrized, in rhor(:,2).

   call timab(70,2,tsec)

  end if

!End of test on varying or fixed occupation numbers
 end if

 deallocate(eknk,kg_k,grnlnk)
 if (psps%usepaw==0) deallocate(enlnk)

 call timab(27,1,tsec)

!In the self-consistent case, diagnose lack of unoccupied
!state (for each spin and k-point).
!Print a warning if the number of such messages already written
!does not exceed mwarning.
 mwarning=5
 if(nwarning<mwarning .and. iscf>0)then
  nwarning=nwarning+1
  bdtot_index=1
  do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
    min_occ=two
    nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
    do iband=1,nband_k
     if(occ(bdtot_index)<min_occ)min_occ=occ(bdtot_index)
     bdtot_index=bdtot_index+1
    end do
    if(min_occ>0.01_dp)then
     if(dtset%nsppol==1)then
      write(message, '(a,a,a,a,i4,a,a,a,f7.3,a,a,a,a,a,a,a)' )ch10,&
&      ' vtorho : WARNING -',ch10,&
&      '  For k-point number ',ikpt,',',ch10,&
&      '  The minimal occupation factor is',min_occ,'.',ch10,&
&      '  An adequate monitoring of convergence requires it to be',&
&           ' at most 0.01_dp.',ch10,&
&      '  Action : increase slightly the number of bands.',ch10
     else
      write(message, '(a,a,a,a,i4,a,a,a,i3,a,f7.3,a,a,a,a,a,a,a)' )ch10,&
&      ' vtorho : WARNING -',ch10,&
&      '  For k-point number ',ikpt,', and',ch10,&
&      '  for spin polarization',isppol,&
&           ', the minimal occupation factor is',min_occ,'.',ch10,&
&      '  An adequate monitoring of convergence requires it to be',&
&            ' at most 0.01_dp.',ch10,&
&      '  Action : increase slightly the number of bands.',ch10
     end if
     call wrtout(6,message,'COLL')
!    It is enough if one lack of adequate occupation is identified, so exit.
      exit
    end if
   end do
  end do
 end if
!In the non-self-consistent case, print eigenvalues and residuals
 if(iscf<=0)then
  option=2 ; enunit=1 ; vxcavg_dum=zero
  call prteigrs(eigen,enunit,energies%e_fermie,filapp,ab_out,iscf,dtset%kptns,dtset%kptopt,dtset%mband,dtset%nband,&
&  dtset%nkpt,nnsclo_now,dtset%nsppol,occ,dtset%occopt,option,&
&  dtset%prteig,prtvol,resid,dtset%tolwfr,vxcavg_dum,dtset%wtk)
 end if

!Find largest residual over bands, k points, and spins,
!except for nbdbuf highest bands
 ibdkpt=1
 residm=zero
 do isppol=1,dtset%nsppol
  do ikpt=1,dtset%nkpt
   nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
   nband_eff=max(1,nband_k-dtset%nbdbuf)
   residm=max(residm,maxval(resid(ibdkpt:ibdkpt+nband_eff-1)))
   ibdkpt=ibdkpt+nband_k
  end do
 end do

 if (iscf>0) then

! Norm-conserving: The density related to WFs (rhowfr) is the total density
  if (psps%usepaw==0) then
   if (optres==1) nvresid=rhor-nvresid
!  Find and print minimum and maximum total electron density and locations
   call prtrhomxmn(6,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%nspden,1,rhor)
  else

! PAW: Build new rhoij quantities from new occ then symetrize them
!      Compute and add the compensation density to rhowfr to get the total density
   call timab(27,2,tsec)
   call timab(555,1,tsec)
   option=1;choice=1
   nsp2=dtset%nsppol;if (dtset%nspden==4) nsp2=4
   do iatom=1,natom
    allocate(pawrhoij(iatom)%rhoij_(pawrhoij(iatom)%lmn2_size,nsp2))
    pawrhoij(iatom)%use_rhoij_=1
   end do
   call status(istep,dtfil%filstat,iexit,level,'call pawmkrhoij')
   if (usecprj==1) then
    call pawmkrhoij(atindx1,cprj,dimcprj,dtset%istwfk,dtset%mband,dtset%mkmem,&
&                   mpi_enreg,natom,nattyp,dtset%nband,dtset%nkpt,dtset%nspden,nspinor,&
&                   dtset%nsppol,ntypat,occ,dtset%pawprtvol,pawrhoij,dtfil%unpaw,dtset%wtk)
   else
    allocate(cprj_tmp(natom,nspinor*dtset%mband*dtset%mkmem*dtset%nsppol*usecprj))
    call cprj_alloc(cprj_tmp,0,dimcprj)
    call ctocprj(atindx,cg,1,cprj_tmp,dtfil,gmet,gprimd,0,0,0,dtset%istwfk,kg,dtset%kptns,&
&                dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&                dtset%natom,nattyp,dtset%nband,dtset%natom,dtset%ngfft,dtset%nkpt,dtset%nloalg,&
&                npwarr,nspinor,dtset%nsppol,ntypat,ph1d,psps,rmet,dtset%typat,&
&                ucvol,dtfil%unpaw,wffnew,xred,ylm)
    call pawmkrhoij(atindx1,cprj_tmp,dimcprj,dtset%istwfk,dtset%mband,dtset%mkmem,&
&                   mpi_enreg,natom,nattyp,dtset%nband,dtset%nkpt,dtset%nspden,nspinor,&
&                   dtset%nsppol,ntypat,occ,dtset%pawprtvol,pawrhoij,dtfil%unpaw,dtset%wtk)
    call cprj_free(cprj_tmp)
    deallocate(cprj_tmp)
   end if
   call symrhoij(choice,psps%indlmn,indsym,psps%lmnmax,natom,dtset%nsym,ntypat,option,&
&                pawang,dtset%pawprtvol,pawrhoij,dtset%symafm,symrec,dtset%typat)
   do iatom=1,natom;deallocate(pawrhoij(iatom)%rhoij_);pawrhoij(iatom)%use_rhoij_=0;end do
   call timab(555,2,tsec)
   call timab(556,1,tsec)
   call pawmknhat(compch_fft,0,0,mpi_enreg,natom,nfftf,pawfgr%ngfft,0,dtset%nspden,ntypat,&
&                 dtset%paral_kgb,pawang,pawfgrtab,rhodum,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
   call timab(556,2,tsec)
   call timab(557,1,tsec)
   call transgrid(1,mpi_enreg,dtset%nspden,+1,1,0,dtset%paral_kgb,pawfgr,rhowfg,rhodum,rhowfr,rhor)
   rhor(:,:)=rhor(:,:)+nhat(:,:);if (optres==1) nvresid=rhor-nvresid
   call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,pawfgr%ngfft,dtset%paral_kgb,0)
!  Find and print minimum and maximum total electron density and locations
   call prtrhomxmn(6,mpi_enreg,nfftf,pawfgr%ngfft,dtset%nspden,1,rhor)
   call timab(557,2,tsec)
   call timab(27,1,tsec)
  end if

! Compute square norm nres2 of density residual nvresid
  if (optres==1) call sqnorm_v(1,mpi_enreg,nfftf,nres2,dtset%nspden,nvresid)

 end if ! iscf>0

 if(psps%usepaw==1.and.(iscf>0.or.iscf==-3)) deallocate(rhowfr,rhowfg)

 call timab(27,2,tsec)

 if(iscf==-1)then

  call status(0,dtfil%filstat,iexit,level,'call tddft    ')

! Eventually compute the excited states within tddft
  if (psps%usepaw==1) then
!  In case of PAW calculation, have to transfer kxc from the fine to the coarse grid:
   allocate(cgrkxc(dtset%nfft,nkxc))
   do ikxc=1,nkxc
    call transgrid(1,mpi_enreg,1,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrkxc(:,ikxc),kxc(:,ikxc))
   end do
   call tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&
&   kg,cgrkxc,dtset%mband,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nfft,&
&   ngfftdiel,dtset%nkpt,nkxc,npwarr,nspinor,dtset%nsppol,occ,ucvol,wffnew)
   deallocate(cgrkxc)
  else
   call tddft(cg,dtfil,dtset,eigen,etotal,gmet,gprimd,gsqcut,&
&   kg,kxc,dtset%mband,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nfft,&
&   ngfftdiel,dtset%nkpt,nkxc,npwarr,nspinor,dtset%nsppol,occ,ucvol,wffnew)
  end if

 else

! Eventually compute the susceptibility matrix and the
! dielectric matrix when istep is equal to 1 or dielstrt
    !  if( (istep==1        .and. dielop>=2) .or. &
    !&     (istep==dielstrt .and. dielop>=1) .or. &
    !&       computesusmat       )then
    call testsusmat(computesusmat,dielop,dielstrt,dtset,istep) !test if the matrix is to be computed
   if(computesusmat) then
   dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
   dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
   dielar(5)=dtset%diegap;dielar(6)=dtset%dielam

   call status(0,dtfil%filstat,iexit,level,'call suscep_st')
   call suscep_stat(atindx1,cg,cprj,densymop_diel,dielar,&
&   dielop,dimcprj,doccde,eigen,gbound_diel,gprimd,&
&   irrzondiel,dtset%istwfk,kg,kg_diel,lmax_diel,&
&   dtset%mband,mgfftdiel,dtset%mkmem,mpi_enreg,dtset%mpw,natom,dtset%nband,nfftdiel,ngfftdiel,&
&   dtset%nkpt,npwarr,npwdiel,dtset%nspden,nspinor,dtset%nsppol,dtset%nsym,ntypat,&
&   occ,dtset%occopt,dtset%paral_kgb,pawang,pawtab,phnonsdiel,ph1ddiel,prtvol,&
&   susmat,dtset%symafm,dtset%symrel,dtset%tnons,dtset%typat,&
&   ucvol,dtfil%unkg,dtfil%unpaw,usecprj,psps%usepaw,wffnew,dtset%wtk,ylmdiel)
!  GMR
!  Print the susceptibility matrix
!  do isp1=1,dtset%nspden
!   do isp2=1,dtset%nspden
!    write(6,'(5x,a,2i2)') 'Susceptibility matrix for spins=',isp1,isp2
!    write(6,'(9x,a,13x,a,10x,a,10x,a)') "g","g'","real","imag"
!    do ipw1=1,10
!     do ipw2=ipw1,10
!      write(6,'(2x,3i4,2x,3i4,2x,f12.8,2x,f12.8)') &
!&      kg_diel(1:3,ipw1),kg_diel(1:3,ipw2),&
!&      susmat(1,ipw1,isp1,ipw2,isp2),susmat(2,ipw1,isp1,ipw2,isp2)
!     end do
!    end do
!   end do
!  end do

  end if

 end if ! end condition on iscf

 deallocate(doccde,EigMin)

!PAW: deallocate <p|c> (cprj)
 if (psps%usepaw==1) then
  if(dtset%mkmem/=0.and.usecprj==1) then
   call cprj_free(cprj)
  end if
  usecprj=0;deallocate(dimcprj,cprj)
 end if

! Rotate labels of disk files when wf i/o is used
  if (dtset%mkmem==0) then
   wfftmp=wffnow ; wffnow=wffnew ; wffnew=wfftmp
  end if

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)') ch10,' vtorho : exit ',&
&   ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(21,2,tsec)

!DEBUG
!write(6,*)' vtorho : exit, residm=',residm
!ENDDEBUG

end subroutine vtorho
!!***
