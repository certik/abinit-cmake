!{\src2tex{textfont=tt}}
!!****f* ABINIT/energy
!! NAME
!! energy
!!
!! FUNCTION
!! Compute electronic energy terms
!! energies%e_eigenvalues, ek and enl from arbitrary (orthonormal) provided wf,
!! ehart, enxc, and eei from provided density and potential,
!! energies%e_eigenvalues=Sum of the eigenvalues - Band energy (Hartree)
!! ek=kinetic energy, ehart=Hartree electron-electron energy,
!! enxc,enxcdc=exchange-correlation energies, eei=local pseudopotential energy,
!! enl=nonlocal pseudopotential energy
!! Also, compute new density from provided wfs, after the evaluation
!! of ehart, enxc, and eei.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, AR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of wavefunction
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem=maximum number of k points which can fit in core memory
!!   | mpw=maximum dimension for number of planewaves
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for polarized
!!   | nsym=number of symmetry elements in space group (at least 1)
!!   | occopt=option for occupancies
!!   | tsmear=smearing energy or temperature (if metal)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gsqcut=G^2 cutoff from gsqcut=ecut/(2 Pi^2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspinor=number of spinorial components of the wavefunctions
!!  n3xccc=dimension of the xccc3d array (0 or nfftf).
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each k point
!!  optene=option for the computation of total energy (direct scheme or double-counting scheme)
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr(natom) <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!   | ntypat=number of types of atoms in cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp(nfftf)=local pseudopotential in real space (hartree)
!!  wffnow=structured array giving all information about wavefunction file
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  compch_fft=-PAW only- compensation charge inside spheres computed over fine fft grid
!!  etotal=total energy (hartree):
!!    - computed by direct scheme if optene=0 or 2
!!    - computed by double-counting scheme if optene=1 or 3
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points (hartree^2)
!!  strsxc(6)=exchange-correlation contribution to stress tensor
!!  vhartr(nfftf)=work space to hold Hartree potential in real space (hartree)
!!  vtrial(nfftf,nspden)=total local potential (hartree)
!!  vxc(nfftf,nspden)=work space to hold Vxc(r) in real space (hartree)
!!
!! SIDE EFFECTS
!!  energies <type(energies_type)>=all part of total energy.
!!   | entropy(IN)=entropy due to the occupation number smearing (if metal)
!!   | e_ewald(IN)=Ewald energy (hartree)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_eigenvalues(OUT)=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_hartree(OUT)=Hartree part of total energy (hartree units)
!!   | e_kinetic(OUT)=kinetic energy part of total energy.
!!   | e_nonlocalpsp(OUT)=nonlocal pseudopotential part of total energy.
!!   | e_xc(OUT)=exchange-correlation energy (hartree)
!!  ==== if optene==0, 2 or 3
!!   | e_localpsp(OUT)=local psp energy (hartree)
!!  ==== if optene==1, 2 or 3
!!   | e_xcdc(OUT)=exchange-correlation double-counting energy (hartree)
!!  rhog(2,nfftf)=work space for rho(G); save intact on return
!!  rhor(nfftf,nspden)=work space for rho(r); save intact on return
!!  nspinor should not be modified in the call of rdnpw
!!  === if psps%usepaw==1 ===
!!    nhat(nfftf,nspden*usepaw)= compensation charge density
!!    pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!
!! NOTES
!!  Be careful to the meaning of nfft (size of FFT grids):
!!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!   - In case of PAW calculations:
!!     Two FFT grids are used; one with nfft points (coarse grid) for
!!     the computation of wave functions ; one with nfftf points
!!     (fine grid) for the computation of total density.
!!
!!  There is a large amount of overhead in the way this routine do the computation of the energy !
!!  For example, the density has already been precomputed, so why to compute it again here ??
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,dotprod_vn,fftpac,hdr_skip,leave_new,leave_test,meanvalue_g,metric
!!      mkffnl,mkkin,mkresi,mkrho,nonlop,pawmknhat,ph1d3d,rdnpw,rhohxc,rwwf
!!      sphereboundary,timab,transgrid,wrtout,xcomm_init,xdefineoff
!!      xmaster_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine energy(atindx,atindx1,cg,compch_fft,densymop_gs,dtfil,dtset,energies,&
& eigen,etotal,gsqcut,indsym,irrzon,kg,mpi_enreg,nattyp,nfftf,ngfft,ngfftf,nhat,&
& nhatgr,nhatgrdim,npwarr,nspinor,n3xccc,occ,optene,paw_ij,pawang,pawfgr,&
& pawfgrtab,pawrhoij,pawtab,phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,&
& usexcnhat,vhartr,vtrial,vpsp,vxc,wffnow,wfs,xccc3d,xred,ylm)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_12spacepar
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13recipspace
 use interfaces_13xc
 use interfaces_14iowfdenpot
 use interfaces_14poisson
 use interfaces_15common, except_this_one => energy
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n3xccc,nfftf,nhatgrdim,optene,usexcnhat
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: compch_fft,etotal
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(dens_sym_operator_type),intent(in) :: densymop_gs
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
 type(wvl_wf_type),intent(in) :: wfs
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
!no_abirules
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 integer, intent(in) :: indsym(4,dtset%nsym,dtset%natom)
 integer :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol),kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: nattyp(psps%ntypat),ngfft(18),ngfftf(18),npwarr(dtset%nkpt),symrec(3,3,dtset%nsym)
 real(dp), intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol),eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp), intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 real(dp), intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
 real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden)
 real(dp), intent(out) :: strsxc(6)
 real(dp), intent(in) :: rprimd(3,3),vpsp(nfftf),xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp), intent(out) :: vhartr(nfftf),vtrial(nfftf,dtset%nspden),vxc(nfftf,dtset%nspden)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 type(paw_ij_type), intent(in) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type), intent(in)  :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,bufdim,choice,cplex,cpopt,dimdij,dimffnl,formeig,ia
 integer :: iatom,iband,icg,ider,idir,ierr,ifft,ig,igs,ii,ikg,ikpt,ilm,ilmn
 integer :: indx,ipw,ipw1,iresid,isp,ispden,ispinor,isppol,istwf_k,itypat
 integer :: master,matblk,mcg,mcg_disk,me,me_distrb,mu,muig,n1,n2,n3,n4,n5,n6
 integer :: nband1,nband_k,nfftotf,nkpg,nkxc,nnlout,npw_k,nsp2,nvloc,option
 integer :: option_rhoij,paw_opt,signs,spaceComm,tim_mkrho,tim_nonlop,tim_rwwf
 real(dp) :: DDOT,arg,doti,dotr,dum,eeigk,ekk,enlk,ucvol,vxcavg
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer,allocatable :: dimlmn(:),kg_dum(:,:),kg_k(:,:)
 real(dp) :: enlout(1),gmet(3,3),gprimd(3,3),kpoint(3),nonlop_dum(1,1)
 real(dp) :: rhodum(1),rmet(3,3),tsec(2),ylmgr_dum(1)
 real(dp),allocatable :: buffer(:),buffer2(:),cg_disk(:,:),cgrvtrial(:,:)
 real(dp),allocatable :: cwavef(:,:),eig_dum(:),eig_k(:),ffnl(:,:,:,:),kinpw(:)
 real(dp),allocatable :: kpg_dum(:,:),kxc(:,:),occ_dum(:),occ_k(:),ph3d(:,:,:)
 real(dp),allocatable :: resid_k(:),rhowfg(:,:),rhowfr(:,:),vlocal(:,:,:,:)
 real(dp),allocatable :: vlocal_tmp(:,:,:),ylm_k(:,:)
 type(cprj_type),allocatable :: cwaveprj(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' energy : enter '
!stop
!ENDDEBUG

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 nfftotf=ngfftf(1)*ngfftf(2)*ngfftf(3)
 if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' energy :  BUG -',ch10,&
&  '  wrong values for nfft, nfftf !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 call timab(59,1,tsec)
!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)

!PATCH energy // KPT & FFT me-->me_kpt & spaceComm --> comm_kpt
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

!Compute Hxc potential from density
 option=1;nkxc=0;
 if (dtset%icoulomb == 0) then
! Use the periodic solver to compute Hxc.
  allocate(kxc(1,nkxc))
  call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,&
&  ngfftf,nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,dtset%nspden,n3xccc,option,rhog,rhor,rprimd,strsxc,&
&  usexcnhat,vhartr,vxc,vxcavg,xccc3d)
  deallocate(kxc)
 else
! Use the free boundary solver.
  call PSolver_rhohxc(dtset, energies%e_hartree, energies%e_xc, energies%e_vxc, &
&  mpi_enreg, rhor, rprimd, vhartr, vxc, vxcavg)
 end if

!Total local potential (for either spin channel) is
!Hartree + local psp + Vxc(spin), minus its mean
!(Note : this potential should agree with the input vtrial)
 do ispden=1,min(dtset%nspden,2)
  do ifft=1,nfftf
   vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
  end do
 end do
 if (dtset%nspden==4) vtrial(:,3:4)=vxc(:,3:4)

!Compute Hartree energy - use up+down rhor
 call dotprod_vn(1,rhor,energies%e_hartree ,doti,mpi_enreg,nfftf,nfftotf,1,1,vhartr,ucvol)
 energies%e_hartree=half*energies%e_hartree

!Compute local psp energy - use up+down rhor
 if (optene/=1) call dotprod_vn(1,rhor,energies%e_localpsp,doti,mpi_enreg,nfftf,nfftotf,1,1,vpsp,ucvol)

!Compute DC-xc energy - use up+down rhor
 if (optene>0) call dotprod_vn(1,rhor,energies%e_xcdc,doti,mpi_enreg,nfftf,nfftotf,dtset%nspden,1,vxc,ucvol)

 if (dtset%mkmem==0) then

! Read wavefunction file header
  call hdr_skip(wffnow,ierr)

! Define offsets, in case of MPI I/O
  formeig=0
  call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,nspinor,dtset%nsppol,dtset%nkpt)

  mcg_disk=dtset%mpw*nspinor*dtset%mband
  allocate(cg_disk(2,mcg_disk))

 end if

 energies%e_eigenvalues=zero
 energies%e_kinetic=zero
 energies%e_nonlocalpsp=zero
 bdtot_index=0
 icg=0


!DEBUG
!write(6,*)' energy : before loop over spins '
!stop
!ENDDEBUG

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 nvloc=1;if(dtset%nspden==4)nvloc=4
 allocate(vlocal(n4,n5,n6,nvloc),kg_k(3,dtset%mpw),cwavef(2,dtset%mpw*nspinor))

!Allocate the arrays of the Hamiltonian whose dimensions do not depend on k
 allocate(gs_hamk%atindx(dtset%natom),gs_hamk%atindx1(dtset%natom))
 allocate(gs_hamk%gbound(2*dtset%mgfft+8,2))
 allocate(gs_hamk%indlmn(6,psps%lmnmax,psps%ntypat))
 allocate(gs_hamk%nattyp(psps%ntypat))
 allocate(gs_hamk%phkxred(2,dtset%natom))
 allocate(gs_hamk%ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom))
 allocate(gs_hamk%pspso(psps%ntypat))
 allocate(gs_hamk%xred(3,dtset%natom))

!Initialize most of the Hamiltonian
 gs_hamk%atindx(:)  =atindx(:)
 gs_hamk%atindx1(:) =atindx1(:)
 gs_hamk%gmet(:,:)  =gmet(:,:)
 gs_hamk%gprimd(:,:)=gprimd(:,:)
 gs_hamk%indlmn(:,:,:)=psps%indlmn(:,:,:)
 gs_hamk%lmnmax     =psps%lmnmax
 gs_hamk%mgfft      =dtset%mgfft
 gs_hamk%mpsang     =psps%mpsang
 gs_hamk%mpssoang   =psps%mpssoang
 gs_hamk%natom      =dtset%natom
 gs_hamk%nattyp(:)  =nattyp(:)
 gs_hamk%nfft       =dtset%nfft
 gs_hamk%ngfft(:)   =dtset%ngfft(:)
 gs_hamk%nloalg(:)  =dtset%nloalg(:)
 gs_hamk%nspinor    =nspinor
 gs_hamk%ntypat     =psps%ntypat
 gs_hamk%nvloc      =nvloc
 gs_hamk%n4         =n4
 gs_hamk%n5         =n5
 gs_hamk%n6         =n6
 gs_hamk%usepaw     =psps%usepaw
 gs_hamk%ph1d(:,:)  =ph1d(:,:)
 gs_hamk%pspso(:)   =psps%pspso(:)
 gs_hamk%ucvol      =ucvol
 gs_hamk%useylm     =psps%useylm
 gs_hamk%xred(:,:)  =xred(:,:)

!Non-local factors:
!Norm-conserving: kleimann-Bylander energies
!PAW: Dij coefficients and overlap coefficients
 if (psps%usepaw==0) then
  gs_hamk%dimekb1=psps%dimekb
  gs_hamk%dimekb2=dtset%ntypat
  allocate(gs_hamk%ekb(psps%dimekb,dtset%ntypat,nspinor**2))
  allocate(gs_hamk%sij(0,0))
  gs_hamk%ekb(:,:,1)=psps%ekb(:,:)
  if (nspinor==2) then
   gs_hamk%ekb(:,:,2)=psps%ekb(:,:)
   gs_hamk%ekb(:,:,3:4)=zero
  end if
 else
  gs_hamk%dimekb1=psps%dimekb*paw_ij(1)%cplex_dij
  gs_hamk%dimekb2=dtset%natom
  allocate(gs_hamk%ekb(gs_hamk%dimekb1,gs_hamk%dimekb2,nspinor**2))
  allocate(gs_hamk%sij(gs_hamk%dimekb1,dtset%ntypat))
  allocate(dimlmn(dtset%natom));ia=0
  do itypat=1,dtset%ntypat
   if (paw_ij(1)%cplex_dij==1) then
    gs_hamk%sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
   else
    do ilmn=1,pawtab(itypat)%lmn2_size
     gs_hamk%sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
     gs_hamk%sij(2*ilmn  ,itypat)=zero
    end do
   end if
   dimlmn(ia+1:ia+nattyp(itypat))=pawtab(itypat)%lmn_size
   ia=ia+nattyp(itypat)
  end do
  allocate(cwaveprj(dtset%natom,nspinor))
  call cprj_alloc(cwaveprj,0,dimlmn)
  deallocate(dimlmn)
  nsp2=dtset%nsppol;if (dtset%nspden==4) nsp2=4
  do iatom=1,dtset%natom
   allocate(pawrhoij(iatom)%rhoij_(pawrhoij(iatom)%lmn2_size,nsp2))
   pawrhoij(iatom)%rhoij_(:,:)=zero
   pawrhoij(iatom)%use_rhoij_=1
  end do
  option_rhoij=1
 end if

!LOOP OVER SPINS
 do isppol=1,dtset%nsppol

! Rewind kpgsph data file if needed:
  if (dtset%mkmem==0) rewind dtfil%unkg
  if (dtset%mkmem==0.and.psps%useylm==1) rewind dtfil%unylm
  ikg=0

! Set up local potential vlocal with proper dimensioning, from vtrial
! Also take into account the spin.
  if(dtset%nspden/=4)then
   if (psps%usepaw==0) then
    call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,vtrial,vlocal,2)
   else
    allocate(cgrvtrial(dtset%nfft,dtset%nspden))
    call transgrid(1,mpi_enreg,dtset%nspden,-1,0,0,dtset%paral_kgb,pawfgr,rhodum,rhodum,cgrvtrial,vtrial)
    call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,cgrvtrial,vlocal,2)
    deallocate(cgrvtrial)
   end if
  else
   allocate(vlocal_tmp(n4,n5,n6))
   if (psps%usepaw==0) then
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

! Loop over k points
  do ikpt=1,dtset%nkpt
   nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
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

!  Continue to initialize the Hamiltonian
   gs_hamk%istwf_k    =istwf_k
   gs_hamk%npw        =npw_k

   allocate(eig_k(nband_k),occ_k(nband_k),resid_k(nband_k))
   allocate(ylm_k(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
   resid_k(:)=0.0_dp
   kpoint(:)=dtset%kptns(:,ikpt)
   gs_hamk%kpoint(:)  =kpoint(:)
   occ_k(:)=occ(1+bdtot_index:nband_k+bdtot_index)
   eig_k(:)=eigen(1+bdtot_index:nband_k+bdtot_index)
   if (minval(eig_k)>1.d100) eig_k=zero
   cplex=2;if (istwf_k>1.and.dtset%nspden/=4) cplex=1

   if (dtset%mkmem==0) then

!   Read sphere data centered at k in dtfil%unkg, then k+g data
    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,dtfil%unkg)
    read (dtfil%unkg) kg_k(1:3,1:npw_k)
    call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)

!   Eventually read spherical harmonics
    if (psps%useylm==1) then
     read(dtfil%unylm)
     read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,psps%mpsang*psps%mpsang)
    end if

!   Read the wavefunction block for ikpt,isppol
    tim_rwwf=3
    allocate(eig_dum(dtset%mband),kg_dum(3,0),occ_dum(dtset%mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,dtset%mband,mcg_disk,mpi_enreg,nband_k,&
&    nband_k,npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
    deallocate(eig_dum,kg_dum,occ_dum)

   else

    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    call sphereboundary(gs_hamk%gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
    if (psps%useylm==1) then
     do ilm=1,psps%mpsang*psps%mpsang
      ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     end do
    end if

!   End if for choice governed by dtset%mkmem
   end if

!  Compute kinetic energy
   allocate(kinpw(npw_k))
   call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass,gmet,kg_k,kinpw,kpoint,npw_k)

   indx=1+icg
   do iband=1,nband_k
!   Compute kinetic energy of each band

    if(mpi_enreg%paral_compil_kpt==1)then
!    Skip this band if not the proper processor
     if (mpi_enreg%paralbd >1) then
      if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= mpi_enreg%me) then
       indx=indx+npw_k*nspinor
       cycle
      end if
     end if
    end if

    if(dtset%mkmem/=0)then
     do ig=1,npw_k*nspinor
      cwavef(1,ig)=cg(1,indx)
      cwavef(2,ig)=cg(2,indx)
      indx=indx+1
     end do
    else
     do ig=1,npw_k*nspinor
      cwavef(1,ig)=cg_disk(1,indx)
      cwavef(2,ig)=cg_disk(2,indx)
      indx=indx+1
     end do
    end if
    call meanvalue_g(ekk,kinpw,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
    energies%e_kinetic=energies%e_kinetic+dtset%wtk(ikpt)*occ_k(iband)*ekk
   end do

   enlk=zero
   eeigk=zero

!  Compute nonlocal form factors ffnl at all (k+G):
   ider=0;dimffnl=1;nkpg=0
   allocate(ffnl(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&   gmet,gprimd,ider,ider,psps%indlmn,kg_k,kpg_dum,kpoint,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&   npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&   psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)

!  Allocate the arrays phkxred and ph3d, compute phkxred
!  and eventually ph3d.
   do ia=1,dtset%natom
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
    matblk=dtset%natom
    allocate(ph3d(2,npw_k,matblk))
    call ph1d3d(1,dtset%natom,kg_k,kpoint,matblk,dtset%natom,npw_k,n1,n2,n3,&
&    gs_hamk%phkxred,ph1d,ph3d)
   end if
   gs_hamk%matblk=matblk

!  DEBUG
!  write(6,*)' energy : before nonlop '
!  stop
!  ENDDEBUG

!  Compute nonlocal psp energy - Norm-conserving only
   do iband=1,nband_k
    if(mpi_enreg%paral_compil_kpt==1)then
!    Skip this band if not the proper processor
     if (mpi_enreg%paralbd >1) then
      if(mpi_enreg%proc_distrb(ikpt,iband,isppol)/= mpi_enreg%me) then
       cycle
      end if
     end if
    end if

    if(dtset%mkmem/=0)cwavef(:,1:npw_k*nspinor)=&
&    cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
    if(dtset%mkmem==0)cwavef(:,1:npw_k*nspinor)=&
&    cg_disk(:,1+(iband-1)*npw_k*nspinor:iband*npw_k*nspinor)

    choice=1-gs_hamk%usepaw ; signs=1 ; idir=0 ; nnlout=1 ; tim_nonlop=3
    paw_opt=gs_hamk%usepaw;cpopt=gs_hamk%usepaw-1
    call nonlop(atindx1,choice,cpopt,cwaveprj,gs_hamk%dimekb1,gs_hamk%dimekb2,dimffnl,dimffnl,&
&    gs_hamk%ekb,enlout,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,&
&    kg_k,kg_k,kpg_dum,kpg_dum,kpoint,kpoint,dum,psps%lmnmax,&
&    matblk,dtset%mgfft,mpi_enreg,psps%mpsang,&
&    psps%mpssoang,dtset%natom,nattyp,dtset%ngfft,nkpg,nkpg,dtset%nloalg,nnlout,npw_k,npw_k,&
&    nspinor,psps%ntypat,0,paw_opt,gs_hamk%phkxred,gs_hamk%phkxred,ph1d,&
&    ph3d,ph3d,psps%pspso,signs,nonlop_dum,nonlop_dum,tim_nonlop,ucvol,&
&    psps%useylm,cwavef,cwavef)

    if (psps%usepaw==0) enlk=enlk+occ_k(iband)*enlout(1)
    eeigk=eeigk+occ_k(iband)*eig_k(iband)

!   PAW: accumulate rhoij
    if (psps%usepaw==1) then
     call pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj,dtset%natom,0,isppol,dtset%natom,&
&     dtset%nspden,nspinor,dtset%nsppol,occ_k(iband),option_rhoij,pawrhoij,dtset%wtk(ikpt))
    end if

   end do

   if (psps%usepaw==0) energies%e_nonlocalpsp=energies%e_nonlocalpsp+dtset%wtk(ikpt)*enlk
   energies%e_eigenvalues=energies%e_eigenvalues+dtset%wtk(ikpt)*eeigk

!  DEBUG
!  write(6,*)' energy : after nonlop '
!  stop
!  ENDDEBUG

!  Compute residual of each band (for informative purposes)
   if(dtset%mkmem/=0)then
    mcg=dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
    call mkresi(cg,dimffnl,eig_k,ffnl,dtfil%filstat,&
&    gs_hamk,icg,ikpt,isppol,&
&    kg_k,kinpw,psps%lmnmax,matblk,mcg,dtset%mgfft,mpi_enreg,&
&    psps%mpsang,psps%mpssoang,dtset%natom,&
&    nband_k,npw_k,&
&    nspinor,psps%ntypat,nvloc,n4,n5,n6,&
&    dtset%paral_kgb,ph3d,dtset%prtvol,&
&    resid_k,psps%usepaw,vlocal)
   else if(dtset%mkmem==0)then
    mcg=dtset%mpw*nspinor*dtset%mband
    call mkresi(cg_disk,dimffnl,eig_k,ffnl,dtfil%filstat,&
&    gs_hamk,icg,ikpt,isppol,&
&    kg_k,kinpw,psps%lmnmax,matblk,mcg,dtset%mgfft,mpi_enreg,&
&    psps%mpsang,psps%mpssoang,dtset%natom,&
&    nband_k,npw_k,&
&    nspinor,psps%ntypat,nvloc,n4,n5,n6,&
&    dtset%paral_kgb,ph3d,dtset%prtvol,&
&    resid_k,psps%usepaw,vlocal)

   end if
   resid(1+bdtot_index : nband_k+bdtot_index) = resid_k(:)

   deallocate(eig_k,occ_k,resid_k)

!  DEBUG
!  write(6,*)' isppol,ikpt',isppol,ikpt
!  write(6,*)resid(1+bdtot_index:nband_k+bdtot_index)
!  ENDDEBUG

   bdtot_index=bdtot_index+nband_k

   if (dtset%mkmem/=0) then
!   Handle case in which kg, cg, are kept in core
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   end if

   deallocate(ffnl,kinpw,ph3d)
   deallocate(ylm_k)

!  End loops on isppol and ikpt
  end do
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
  call leave_test(mpi_enreg)
  write(message,*) 'energy: loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
 end if

 deallocate(gs_hamk%atindx,gs_hamk%atindx1)
 deallocate(gs_hamk%ekb,gs_hamk%sij)
 deallocate(gs_hamk%gbound)
 deallocate(gs_hamk%indlmn)
 deallocate(gs_hamk%nattyp)
 deallocate(gs_hamk%phkxred)
 deallocate(gs_hamk%pspso)
 deallocate(gs_hamk%ph1d)
 deallocate(gs_hamk%xred)

 if(dtset%mkmem==0) deallocate(cg_disk)

!DEBUG
!write(6,*)' energy : after loop on kpts and spins '
!stop
!ENDDEBUG

 if(mpi_enreg%paral_compil_kpt==1)then
! Accumulate enl eeig and ek on all proc.
  allocate(buffer(3+dtset%mband*dtset%nkpt*dtset%nsppol))
  buffer(1)=energies%e_nonlocalpsp ; buffer(2)=energies%e_kinetic ; buffer(3)=energies%e_eigenvalues
  do iresid=1,dtset%mband*dtset%nkpt*dtset%nsppol
   buffer(iresid+3)=resid(iresid)
  end do
  call timab(48,1,tsec)
  call xsum_mpi(buffer,spaceComm,ierr)
  call timab(48,2,tsec)
  energies%e_nonlocalpsp=buffer(1) ; energies%e_kinetic=buffer(2) ; energies%e_eigenvalues=buffer(3)
  do iresid=1,dtset%mband*dtset%nkpt*dtset%nsppol
   resid(iresid)=buffer(iresid+3)
  end do
  deallocate(buffer)
! Accumulate rhoij_
  if (psps%usepaw==1) then
   call timab(48,1,tsec)
   allocate(dimlmn(dtset%natom))
   dimlmn(1:dtset%natom)=pawrhoij(1:dtset%natom)%cplex*pawrhoij(1:dtset%natom)%lmn2_size
   nsp2=pawrhoij(1)%nsppol;if (pawrhoij(1)%nspden==4) nsp2=4
   bufdim=sum(dimlmn)*nsp2
   allocate(buffer(bufdim),buffer2(bufdim))
   ii=0
   do iatom=1,dtset%natom
    do isppol=1,nsp2
     buffer(ii+1:ii+dimlmn(iatom))=pawrhoij(iatom)%rhoij_(1:dimlmn(iatom),isppol)
     ii=ii+dimlmn(iatom)
    end do
   end do
   call xsum_mpi(buffer,buffer2,bufdim,spaceComm,ierr) !Build sum of everything
   ii=0
   do iatom=1,dtset%natom
    do isppol=1,nsp2
     pawrhoij(iatom)%rhoij_(1:dimlmn(iatom),isppol)=buffer2(ii+1:ii+dimlmn(iatom))
     ii=ii+dimlmn(iatom)
    end do
   end do
   deallocate(buffer,buffer2,dimlmn)
   call timab(48,2,tsec)
  end if
 end if

!Compute total (free) energy
 if (optene==0.or.optene==2) then
  if (psps%usepaw==0) then
   etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&   energies%e_localpsp + energies%e_corepsp + energies%e_nonlocalpsp
  else
   etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&   energies%e_localpsp + energies%e_corepsp + energies%e_paw
  end if
 else if (optene==1.or.optene==3) then
  if (psps%usepaw==0) then
   etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&   energies%e_xcdc + energies%e_corepsp
  else
   etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&   energies%e_xcdc + energies%e_corepsp + energies%e_pawdc
  end if
 end if
 etotal = etotal + energies%e_ewald
 if(dtset%occopt>=3 .and. dtset%occopt<=7) etotal=etotal-dtset%tsmear*energies%entropy

!Compute new charge density based on incoming wf
!Keep rhor and rhog intact for later use e.g. in stress.
!=== Norm-conserving psps: simply compute rho from WFs
 if (psps%usepaw==0) then
  tim_mkrho=3
  call mkrho(cg,densymop_gs,dtset,irrzon,kg,mpi_enreg,&
&  npwarr,nspinor,occ,phnons,rhog,rhor,tim_mkrho,ucvol,dtfil%unkg,wffnow,wfs)
 else
! === PAW case: add compensation charge density
  tim_mkrho=3;option=1;choice=1
  call symrhoij(choice,psps%indlmn,indsym,psps%lmnmax,dtset%natom,dtset%nsym,dtset%ntypat,option,&
&  pawang,dtset%pawprtvol,pawrhoij,dtset%symafm,symrec,dtset%typat)
  do iatom=1,dtset%natom
   deallocate(pawrhoij(iatom)%rhoij_)
   pawrhoij(iatom)%use_rhoij_=0
  end do
  call pawmknhat(compch_fft,0,0,mpi_enreg,dtset%natom,nfftf,ngfftf,0,dtset%nspden,dtset%ntypat,&
&  dtset%paral_kgb,pawang,pawfgrtab,rhodum,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
  allocate(rhowfr(dtset%nfft,dtset%nspden),rhowfg(2,dtset%nfft));rhowfr(:,:)=zero
  call mkrho(cg,densymop_gs,dtset,irrzon,kg,mpi_enreg,&
&  npwarr,nspinor,occ,phnons,rhowfg,rhowfr,tim_mkrho,ucvol,dtfil%unkg,wffnow,wfs)
  call transgrid(1,mpi_enreg,dtset%nspden,+1,1,0,dtset%paral_kgb,pawfgr,rhowfg,rhodum,rhowfr,rhor)
  deallocate(rhowfr,rhowfg)
  rhor(:,:)=rhor(:,:)+nhat(:,:)
  call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
 end if

 write(message, '(a,a,a,a)' )ch10, &
& ' energy: COMMENT -',ch10,&
& '  New density rho(r) made from input wfs'
 call wrtout(06,message,'COLL')

 call timab(59,2,tsec)

 deallocate(cwavef,kg_k,vlocal)
 if (psps%usepaw==1) then
  call cprj_free(cwaveprj)
  deallocate(cwaveprj)
 end if

end subroutine energy
!!***
