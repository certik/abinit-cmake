!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfcv
!! NAME
!! scfcv
!!
!! FUNCTION
!! Self-consistent-field convergence.
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! ground state and optionally to compute forces and energy.
!! This routine is called to compute forces for given atomic
!! positions or else to do non-SCF band structures.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, GMR, AR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cpus= cpu time limit in seconds
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs for the "coarse" grid (see NOTES below)
!!   | mkmem =number of k points which can fit in memory; set to 0 if use disk
!!   | mpw=maximum dimensioned size of npw.
!!   | natom=number of atoms in cell.
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   |      for the "coarse" grid (see NOTES below)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in space group
!!  ecore=core psp energy (part of total energy) (hartree)
!!  iapp=indicates the eventual suffix to be appended to the generic output root
!!         if 0 : no suffix to be appended (called directly from gstate)
!!         if positive : append "_TIM//iapp" (called from move or brdmin)
!!         if -1 : append "_TIM0" (called from brdmin)
!!         if -2, -3, -4, -5: append "_TIMA", ... ,"_TIMD", (called from move)
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)= # atoms of each type.
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!                           (see initberry.f)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!
!! SIDE EFFECTS
!!  acell(3)=length scales of primitive translations (bohr)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=updated wavefunctions;  if mkmem>=nkpt,
!!         these are kept in a disk file.
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!      calculations (see initberry.f)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  initialized= if 0 the initialization of the gstate run is not yet finished
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  nfftf=(effective) number of FFT grid points (for this processor)
!!       for the "fine" grid (see NOTES below)
!!  occ(mband*nkpt*nsppol)=occupation number for each band (often 2) at each k point
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= -PAW only- atomic occupancies
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!   (should be made a pure output quantity)
!!  rhog(2,nfftf)=array for Fourier transform of electron density
!!  rhor(nfftf,nspden)=array for electron density in el./bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  wffnew,wffnow=struct info for wf disk files.
!!  wvl <type(wvl_data)>=all wavelets data.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)= at input, previous reduced dimensionless atomic coordinates
!!                     at output, current xred is transferred to xred_old
!!
!! NOTES
!! It is worth to explain THE USE OF FFT GRIDS:
!! ============================================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      brdmin,delocint,diisrelax,gstate,moldyn,move
!!
!! CHILDREN
!!      abi_etsf_init,afterscfloop,berryphase_new,calc_xc_ep,chkdilatmx
!!      chkpawovlp,cprj_alloc,ctocprj,energies_init,energy,etotfor,extraprho
!!      fappnd,fourdp,fresid,getcut,getmpw,getng,getph,initylmg,int2char4,ioarr
!!      kpgio,leave_new,leave_test,metric,newrho,newvtr,nhatgrid,odamix
!!      out_geometry_xml,out_resultsgs_xml,outscfcv,pawdenpot,pawdij,pawmknhat
!!      rhohxc,rhotov,scprqt,setnoccmmp,setsym,setup_positron,setvtr
!!      sphereboundary,status,symdij,symrhg,symzat,timab,vtorho,vtorhorec
!!      vtorhotf,wrtout,wvl_mkrho,wvl_newvtr,wvl_wfsinp_reformat,xcomm_world
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&  dtset,ecore,eigen,hdr,iapp,indsym,initialized,&
&  irrzon,kg,mpi_enreg,nattyp,nfftf,npwarr,nspinor,occ,&
&  pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&  phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&  scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
#if defined HAVE_ETSF_IO
 use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13ionetcdf
 use interfaces_13iovars
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13recipspace
 use interfaces_13xc
 use interfaces_14iowfdenpot
 use interfaces_15common
 use interfaces_15recursion
 use interfaces_18seqpar
 use interfaces_21drive, except_this_one => scfcv
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iapp,pwind_alloc
 integer,intent(inout) :: initialized,nfftf,nspinor
 real(dp),intent(in) :: cpus,ecore
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(dens_sym_operator_type),intent(inout) :: densymop_gs
 type(efield_type),intent(inout) :: dtefield
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(inout) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(results_gs_type),intent(inout) :: results_gs
 type(scf_history_type),intent(inout) :: scf_history
 type(wffile_type),intent(inout) :: wffnew,wffnow
 type(wvl_data),intent(inout) :: wvl
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(inout) :: indsym(4,dtset%nsym,dtset%natom)
!no_abirules
 integer, intent(inout) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem)
 integer, intent(in) :: nattyp(psps%ntypat),npwarr(dtset%nkpt),pwind(pwind_alloc,2,3)
 integer, intent(inout) :: symrec(3,3,dtset%nsym)
 real(dp), intent(inout) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp), intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
  !(nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise)
 real(dp), intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp), intent(inout) :: acell(3),rprimd(3,3)
 real(dp), pointer :: rhog(:,:),rhor(:,:)
 real(dp), intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp), intent(inout) :: xred(3,dtset%natom),xred_old(3,dtset%natom)
 real(dp), intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp), intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 type(pawrhoij_type), intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawrad_type), intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables -------------------------
!Variables for partial dos calculation
!scalars
 integer,parameter :: level=6,response=0
 integer :: accessfil,afford,bantot,bdtot_index,berryopt_tmp,choice
 integer :: computed_forces,dbl_nnsclo,dielop,dielstrt,dimdmat,dtset_iprcel,dwr
 integer :: fformr,fftalg,forces_needed,get_ek,get_entropy,iatom,iband,ider
 integer :: ierr,iexit,ii,ikpt,ilmn,iln,impose_dmat,initialized0,ipw,ipw1,ipw2
 integer :: ir,isave,iscf10,ispden,ispmix,isppol,istep,itypat,j0lmn,j0ln,jj
 integer :: jlmn,jln,klmn,kln,kssform,lm_size,lmax_diel,lmn2_size,ln_size
 integer :: lpawumax,mesh_size,mgfftdiel,mgfftf,moved_atm_inside,moved_rhor
 integer :: n1xccc,n3xccc,n_fftgr,n_index,nele,nfftdiel,nfftmix,nhatgrdim,nkxc
 integer :: npawmix,npwdiel,nstep,nzlmopt,offset,optberry,optene,optfor,optgr0
 integer :: optgr1,optgr2,option,optrad,optres,optstr,optxc,prtden,prtfor,quit
 integer :: quit_sum,quitrec,rdwr,rdwrpaw,spaceComm,stress_needed,unit_out
 integer :: usecprj,usexcnhat,v_size,zval
 real(dp) :: boxcut,boxcutdg,compch_fft,compch_sph,deltae,diecut,diffor,ecut
 real(dp) :: ecutf,ecutsus,edum,elast,enxcsr,etotal,fermie,gsqcut,lifetime
 real(dp) :: maxfor,ratio,rcut_coulomb,rcut_coulomb1,res2,residm,ucvol,val_max
 real(dp) :: val_min,vxcavg,vxcavg_dum
 logical :: ex,lusedmat
 character(len=4) :: tag
 character(len=500) :: message
 character(len=fnlen) :: filapp,fildata,filfft,filkgs,filprot,kgnam
 type(MPI_type) :: mpi_enreg_diel
 type(dens_sym_operator_type) :: densymop_diel
 type(energies_type) :: energies
!arrays
 integer :: ngfft(18),ngfftdiel(18),ngfftf(18),ngfftmix(18),npwarr_diel(1)
 integer :: npwtot_diel(1)
 integer,allocatable :: dimcprj(:),gbound_diel(:,:),i_rhor(:),i_vresid(:)
 integer,allocatable :: i_vrespc(:),i_vtrial(:),irrzondiel(:,:,:),kg_diel(:,:)
 integer,allocatable :: lmn_size(:)
 real(dp) :: dielar(7),dphase(3),dummy2(6),favg(3),gmet(3,3),gprimd(3,3),k0(3)
 real(dp) :: kpt_diel(3),pel(3),pel_cg(3),pelev_dum(3),pion(3),ptot(3)
 real(dp) :: rhodum(1),rmet(3,3),rprim(3,3),strsxc(6),strten(6),tollist(12)
 real(dp) :: tsec(2),vnew_mean(dtset%nspden),vres_mean(dtset%nspden)
 real(dp),allocatable :: dielinv(:,:,:,:,:),dtn_pc(:,:),dummy(:),ehart1(:)
 real(dp),allocatable :: excapn(:),f_atm(:,:,:),f_fftgr(:,:,:),f_paw(:,:)
 real(dp),allocatable :: fcart(:,:),forold(:,:),fred(:,:),gresid(:,:)
 real(dp),allocatable :: grewtn(:,:),grhf(:,:),grnl(:),grtn(:,:),grxc(:,:)
 real(dp),allocatable :: kxc(:,:),nhat(:,:),nhatgr(:,:,:),nvresid(:,:)
 real(dp),allocatable :: nvresidtemp(:,:),ph1d(:,:),ph1ddiel(:,:),ph1df(:,:)
 real(dp),allocatable :: phnonsdiel(:,:,:),qphon(:),rhocore(:),rhore(:,:)
 real(dp),allocatable :: rhorp(:,:),rhotote(:),rhototp(:),rsepts(:)
 real(dp),allocatable :: rseptstot(:),rsppts(:),rspts(:),shiftvector(:)
 real(dp),allocatable :: susmat(:,:,:,:,:),synlgr(:,:),veff(:),vhae(:,:)
 real(dp),allocatable :: vhap(:,:),vhartr(:),vnew(:,:),vpsp(:),vtrial(:,:)
 real(dp),allocatable :: vxc(:,:),vxcapn(:),workr(:,:),xccc3d(:),ylmdiel(:,:)
 type(cprj_type),allocatable :: cprj(:,:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)

! *********************************************************************

!DEBUG
!write(6,*)' scfcv : enter'
!stop
!ENDDEBUG

 call timab(20,1,tsec)
 call timab(54,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 if(dtset%prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,' scfcv : enter '
  call wrtout(06,message,'COLL')
 end if

!######################################################################
!Initializations - Memory allocations
!----------------------------------------------------------------------

 call status(0,dtfil%filstat,iexit,level,'allocate/init ')

!WVL - reformat the wavefunctions in the case of xred != xred_old
 if (dtset%usewvl == 1 .and. maxval(xred_old - xred) > zero) then
! WVL - Before running scfcv, on non-first geometry step iterations, we need
! to reformat the wavefunctions, taking into acount the new
! coordinates.
! We prepare to change rhog (to be removed) and rhor.
  deallocate(rhog)
  deallocate(rhor)

  call wvl_wfsinp_reformat(acell, dtset, mpi_enreg, occ, psps, &
&  rprimd, wvl, xred, xred_old)
  nfftf = dtset%nfft

  allocate(rhog(2, dtset%nfft))
  allocate(rhor(2, dtset%nfft))
  call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wvl%wfs)
 end if


!Save some variables from dataset definition
 nstep=dtset%nstep
 dtset_iprcel = dtset%iprcel
 ecut=dtset%ecut
 ecutf=ecut;if (psps%usepaw==1.and.pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg
 iscf10=mod(dtset%iscf,10)
 tollist(1)=dtset%tolmxf;tollist(2)=dtset%tolwfr
 tollist(3)=dtset%toldff;tollist(4)=dtset%toldfe
 tollist(6)=dtset%tolvrs;tollist(7)=dtset%tolrff

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Some variables need to be initialized/nullify at start
 quit=0 ; dbl_nnsclo=0 ;
 dielop=0 ; strsxc=zero
 deltae=zero ; elast=zero ;
 if (dtset%nstep==0) fermie=energies%e_fermie
 call energies_init(energies)
 if (dtset%nstep==0) energies%e_fermie=fermie
 energies%e_corepsp     = ecore / ucvol
 isave=0 !initial index of density protection file
 optres=merge(0,1,dtset%iscf<10)
 usexcnhat=0;usecprj=0
 initialized0=initialized

!Stresses and forces flags
 forces_needed=0;prtfor=0
 if ((dtset%optforces==1.or.dtset%ionmov==4.or.abs(tollist(3))>tiny(0._dp)).and.(dtset%positron==0)) then
  if (dtset%iscf>0.and.nstep>0) forces_needed=1
  if (nstep==0) forces_needed=2
  prtfor=1
 else if (dtset%iscf>0.and.dtset%positron==0.and.dtset%optforces==2) then
  forces_needed=2
 end if
 stress_needed=0
 if (dtset%optstress>0.and.dtset%iscf>0.and.dtset%prtstm==0.and.dtset%positron==0.and. &
& (nstep>0.or.dtfil%ireadwf==1)) stress_needed=1

!This is only needed for the tddft routine, and does not
!correspond to the intented use of results_gs (should be only
!for output of scfcv
 etotal  =results_gs%etotal

!Get FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 ngfft(:)=dtset%ngfft(:)
 if (psps%usepaw==1) then
  mgfftf=pawfgr%mgfft;ngfftf(:)=pawfgr%ngfft(:)
 else
  mgfftf=dtset%mgfft;ngfftf(:)=ngfft(:)
 end if

!Prepare the name of the _FFT file and the _KGS file
 filfft=trim(dtfil%filnam_ds(5))//'_FFT'
 filkgs=trim(dtfil%filnam_ds(5))//'_KGS'
!There is a definite problem in the treatment of // by CPP ...
 if(mpi_enreg%paral_compil_kpt==1 .or. mpi_enreg%paral_compil_fft==1)then
  call int2char4(mpi_enreg%me,tag)
  filfft=trim(filfft)//'_P-'//tag
  filkgs=trim(filkgs)//'_P-'//tag
 end if
!Prepare the name of the auxiliary files DOS, EIG...
 call fappnd(filapp,dtfil%filnam_ds(4),iapp)
!Prepare the name of the auxiliary files for protection
 call fappnd(filprot,dtfil%filnam_ds(5),iapp)

!We create Output files when required
 if (dtset%accesswff == 3) then
#if defined HAVE_ETSF_IO
! Compute this lmn_size stuff
  allocate(lmn_size(psps%npsp))
  if(psps%usepaw==1) then
   lmn_size(:) = pawtab(1:psps%npsp)%lmn_size
   ! DC: in PAW, the scalar quantities like the densities are on the fine grid.
   dtset%ngfft = ngfftf
  else
   lmn_size(:) = psps%lmnmax
  end if
! Create an ETSF file for each required files
  if (dtset%prtden /= 0) then
!  Case of density.
   write(fildata, "(A,A)") trim(filapp), "_DEN"
   call abi_etsf_init(dtset, fildata, 1, .true., lmn_size, psps, wvl%wfs)
  end if
  if (dtset%prtwf == 1) then
!  Case of wavefunctions.
   if(mpi_enreg%paral_compil_kpt==1)then
    write(fildata, "(A,A,I0)") trim(filapp), "_WFK_", mpi_enreg%me
   else
    write(fildata, "(A,A)") trim(filapp), "_WFK"
   end if
   call abi_etsf_init(dtset, fildata, 2, .true., lmn_size, psps, wvl%wfs)
  end if
  if (dtset%prtvxc /= 0) then
!  Case of Exchange-correlation potential.
   write(fildata, "(A,A)") trim(filapp), "_VXC"
   call abi_etsf_init(dtset, fildata, 24, .true., lmn_size, psps, wvl%wfs)
  end if
! FIXME: append other possibilities.
! * ETSFIO cases of only correlation or only exchange are not abinit
! options
! * the fix for VHA,VHXC,POT,STM is dirty: they are flagged as
! exchange correlation pot files, except stm, which is flagged density.
! START dirty treatment
  if (dtset%prtvha /= 0) then
!  Case of Hartree potential.
   write(fildata, "(A,A)") trim(filapp), "_VHA"
   call abi_etsf_init(dtset, fildata, 24, .true., lmn_size, psps, wvl%wfs)
  end if
  if (dtset%prtvhxc /= 0) then
!  Case of Hartree+XC potential.
   write(fildata, "(A,A)") trim(filapp), "_VHXC"
   call abi_etsf_init(dtset, fildata, 24, .true., lmn_size, psps, wvl%wfs)
  end if
  if (dtset%prtpot /= 0) then
!  Case of total potential.
   write(fildata, "(A,A)") trim(filapp), "_POT"
   call abi_etsf_init(dtset, fildata, 24, .true., lmn_size, psps, wvl%wfs)
  end if
  if (dtset%prtstm /= 0) then
!  Case of STM output.
   write(fildata, "(A,A)") trim(filapp), "_STM"
   call abi_etsf_init(dtset, fildata, 1, .true., lmn_size, psps, wvl%wfs)
  end if
! END dirty treatment
  if(psps%usepaw==1) then
   dtset%ngfft = ngfft
  end if

  deallocate(lmn_size)
#endif
 end if

!Entering a scfcv loop, printing data to XML file if required.
 if (mpi_enreg%me == 0 .and. dtset%outputXML == 1) then
! scfcv() will handle a scf loop, so we output the scfcv markup.
  write(ab_xml_out, "(A)") '    <scfcvLoop>'
  write(ab_xml_out, "(A)") '      <initialConditions>'
! We output the geometry of the dataset given in argument.
! xred and rprimd are given independently since dtset only
! stores original and final values.
  call out_geometry_XML(dtset, 4, dtset%natom, rprimd, xred)
  write(ab_xml_out, "(A)") '      </initialConditions>'
 end if

!Examine tolerance criteria, and eventually  print a line to the output
!file (with choice=1, the only non-dummy arguments of scprqt are
!nstep, tollist and iscf - still, diffor and res2 are here initialized to 0)
 choice=1 ; diffor=zero ; res2=zero
 allocate(fcart(3,dtset%natom),fred(3,dtset%natom))
 fred(:,:)=zero
 fcart(:,:)=results_gs%fcart(:,:) ! This is a side effect ...
!results_gs should not be used as input of scfcv
 call scprqt(choice,cpus,deltae,diffor,dtset,&
& eigen,etotal,favg,fcart,energies%e_fermie,filapp,dtfil%filnam_ds(1),&
& initialized0,dtset%iscf,istep,dtset%kptns,maxfor,moved_atm_inside,mpi_enreg,&
& dtset%nband,dtset%nkpt,nstep,occ,optres,&
& prtfor,quit,res2,resid,residm,response,tollist,psps%usepaw,&
& vxcavg,dtset%wtk,xred)

!Various allocations (potentials, gradients, ...)
 allocate(forold(3,dtset%natom),grnl(3*dtset%natom),gresid(3,dtset%natom),&
& grewtn(3,dtset%natom),grxc(3,dtset%natom),synlgr(3,dtset%natom))
 allocate(ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom),ph1df(2,3*(2*mgfftf+1)*dtset%natom))
 allocate(vhartr(nfftf),vtrial(nfftf,dtset%nspden),vpsp(nfftf),vxc(nfftf,dtset%nspden))
 forold(:,:)=zero ; gresid(:,:)=zero ; pel(:)=zero
 n1xccc=0;if (psps%n1xccc/=0.and.dtset%positron==0) n1xccc=psps%n1xccc
 n3xccc=0;if (psps%n1xccc/=0.and.dtset%positron==0) n3xccc=nfftf
 allocate(xccc3d(n3xccc))

!Allocations for positron only
 if (dtset%positron>0) then
  allocate(vhae(nfftf,2),vhap(nfftf,2),rhorp(nfftf,dtset%nspden))
  allocate(rhotote(nfftf),rhototp(nfftf),rsepts(nfftf),vxcapn(nfftf),rspts(nfftf),&
&  rhore(nfftf,dtset%nspden),excapn(nfftf),rsppts(nfftf),rhocore(nfftf))
 end if

!Allocations/initializations for PAW only
 if(psps%usepaw==1) then

! Variables/arrays related to the fine FFT grid
  allocate(nhat(nfftf,dtset%nspden));if (nstep==0) nhat=zero
  allocate(pawfgrtab(dtset%natom))
  do iatom=1,dtset%natom
   pawfgrtab(iatom)%l_size=pawtab(dtset%typat(iatom))%lcut_size
   pawfgrtab(iatom)%nfgd=0;allocate(pawfgrtab(iatom)%ifftsph(0))
   pawfgrtab(iatom)%rfgd_allocated=0;allocate(pawfgrtab(iatom)%rfgd(0,0))
   pawfgrtab(iatom)%gylm_allocated=0;allocate(pawfgrtab(iatom)%gylm(0,0))
   pawfgrtab(iatom)%gylmgr_allocated=0;allocate(pawfgrtab(iatom)%gylmgr(0,0,0))
   pawfgrtab(iatom)%gylmgr2_allocated=0;allocate(pawfgrtab(iatom)%gylmgr2(0,0,0))
  end do
  compch_fft=-1.d5
  usexcnhat=maxval(pawtab(:)%vlocopt)
  if (usexcnhat==0.and.dtset%ionmov==4.and.dtset%iscf<10) then
   write(message, '(a,a,a,a)' ) ch10,&
&   ' scfcv :  ERROR -',ch10,&
&   '  You cannot simultaneously use ionmov=4 and such a PAW psp file !'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! Variables/arrays related to the PAW spheres
  allocate(paw_ij(dtset%natom),paw_an(dtset%natom))
  compch_sph=-1.d5;lpawumax=-1
  do iatom=1,dtset%natom
   itypat=dtset%typat(iatom)
   lmn2_size=pawtab(itypat)%lmn2_size
   lm_size=pawtab(itypat)%lcut_size**2
   paw_an(iatom)%cplex     =1
   paw_an(iatom)%has_vxcval=0
   paw_an(iatom)%angl_size=pawang%angl_size
   paw_an(iatom)%mesh_size=pawtab(itypat)%mesh_size
   paw_an(iatom)%nspden   =dtset%nspden
   paw_an(iatom)%lm_size  =lm_size
   paw_an(iatom)%has_vxcval=0
   allocate(paw_an(iatom)%lmselect(lm_size))
   paw_ij(iatom)%cplex    =1
   paw_ij(iatom)%cplex_dij=nspinor
   paw_ij(iatom)%nspden   =dtset%nspden
   paw_ij(iatom)%nsppol   =dtset%nsppol
   paw_ij(iatom)%lmn_size =pawtab(itypat)%lmn_size
   paw_ij(iatom)%lmn2_size=lmn2_size
   paw_ij(iatom)%ndij     =max(nspinor**2,dtset%nspden)
   paw_ij(iatom)%has_dijhat=0 ;paw_ij(iatom)%has_dijso=0
   paw_ij(iatom)%has_dijU=0   ;paw_ij(iatom)%has_dijxc=0
   paw_ij(iatom)%has_dijxc_val=0
   paw_ij(iatom)%has_dijxc=0
   nullify(paw_ij(iatom)%dijxc)
   paw_ij(iatom)%has_dijxc_val=0
   nullify(paw_ij(iatom)%dijxc_val)
   paw_ij(iatom)%has_dijso=0
   nullify(paw_ij(iatom)%dijso)
   paw_ij(iatom)%has_dijU  =0
   nullify(paw_ij(iatom)%dijU)
   allocate(paw_ij(iatom)%dij(paw_ij(iatom)%cplex_dij*lmn2_size,paw_ij(iatom)%ndij))
   if (dtset%iscf==22) then
    paw_ij(iatom)%has_dijhat=1
    allocate(paw_ij(iatom)%dijhat(paw_ij(iatom)%cplex_dij*lmn2_size,paw_ij(iatom)%ndij))
   end if
   if (pawtab(itypat)%usepawu>0) then
    lpawumax=max(pawtab(itypat)%lpawu,lpawumax)
    allocate(paw_ij(iatom)%noccmmp(2*pawtab(itypat)%lpawu+1,&
&    2*pawtab(itypat)%lpawu+1,dtset%nspden))
    allocate(paw_ij(iatom)%nocctot(dtset%nspden))
   end if
   if (pawtab(itypat)%useexexch>0) then
    allocate(paw_ij(iatom)%vpawx(1,lmn2_size,dtset%nspden))
   end if
   pawrhoij(iatom)%use_rhoijres=1
   allocate(pawrhoij(iatom)%rhoijres(lmn2_size,dtset%nspden))
   do ispden=1,dtset%nspden
    pawrhoij(iatom)%rhoijres(:,ispden)=zero
   end do
  end do

 end if ! PAW

!WVL - since wavelets change the size of the box, dont
!need dilatmax.
 if (dtset%usewvl == 0) then
! Check that the possible change of unit cell size has not lead
! to a too large increase
  call chkdilatmx(dtset%dilatmx,rprimd,dtset%rprimd_orig)
 end if

!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (dtset%iscf>0) then
  dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
  dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
  dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
  allocate(nvresid(nfftf,dtset%nspden));if (nstep==0) nvresid=zero
  allocate(dtn_pc(3,dtset%natom))
  if(iscf10==1) then
!  For iscf==1, five additional vectors are needed
!  The index 1 is attributed to the old trial potential,
!  The new residual potential, and the new
!  preconditioned residual potential receive now a temporary index
!  The indices number 4 and 5 are attributed to work vectors.
   n_fftgr=5 ; n_index=1
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3
  else if(iscf10==2) then
!  For iscf==2, three additional vectors are needed.
!  The index number 1 is attributed to the old trial vector
!  The new residual potential, and the new preconditioned
!  residual potential, receive now a temporary index.
   n_fftgr=3 ; n_index=1
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3
  else if(iscf10==3) then
!  For iscf==3 , four additional vectors are needed.
!  The index number 1 is attributed to the old trial vector
!  The new residual potential, and the new and old preconditioned
!  residual potential, receive now a temporary index.
   n_fftgr=4 ; n_index=2
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3 ; i_vrespc(2)=4
  else if (iscf10==4) then
!  For iscf==4 , six additional vectors are needed.
!  The indices number 1 and 2 are attributed to two old trial vectors
!  The new residual potential, and the new and two old preconditioned
!  residual potentials, receive now a temporary index.
   n_fftgr=6 ; n_index=3
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vtrial(2)=2 ; i_vresid(1)=3
   i_vrespc(1)=4 ; i_vrespc(2)=5 ; i_vrespc(3)=6
  else if((iscf10==5).or.(iscf10==6)) then
!  For iscf==5 or 6, ten additional vectors are needed
!  The index number 1 is attributed to the old trial vector
!  The index number 6 is attributed to the search vector
!  Other indices are attributed now. Altogether ten vectors
   n_fftgr=10 ; n_index=3
   allocate(i_rhor(n_index),i_vtrial(n_index),i_vresid(n_index),i_vrespc(n_index))
   i_vtrial(1)=1 ; i_vresid(1)=2 ; i_vrespc(1)=3 ; i_vresid(2)=4 ; i_vrespc(2)=5
   i_vresid(3)=7 ; i_vrespc(3)=8 ; i_rhor(2)=9 ; i_rhor(3)=10
  else if(iscf10==7) then
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
! The next arrays are needed if iscf==5 and ionmov==4,
! but for the time being, they are always allocated
  allocate(grhf(3,dtset%natom),f_atm(3,dtset%natom,n_fftgr))
! Additional allocation for mixing within PAW
  npawmix=0
  if(psps%usepaw==1) then
   do iatom=1,dtset%natom
    itypat=dtset%typat(iatom)
    allocate(pawrhoij(iatom)%kpawmix(pawtab(itypat)%lmnmix_sz))
    pawrhoij(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
    pawrhoij(iatom)%kpawmix=pawtab(itypat)%kmix
    npawmix=npawmix+dtset%nspden*pawtab(itypat)%lmnmix_sz
   end do
  end if
 end if ! iscf>0

!Here, allocate arrays for computation of susceptibility and dielectric matrix or for TDDFT
 if( (nstep>0 .and. dtset%iscf>0) .or. dtset%iscf==-1 ) then !MF
! The dielectric stuff is performed in sequential mode.
! Set mpi_enreg_diel accordingly
  mpi_enreg_diel%paral_compil_fft=0
  mpi_enreg_diel%paral_compil_kpt=0
  mpi_enreg_diel%mode_para='n'
  mpi_enreg_diel%me=0
  mpi_enreg_diel%me_fft=0
  mpi_enreg_diel%me_kpt=0
  mpi_enreg_diel%nproc_fft=1
  mpi_enreg_diel%fft_option_lob=0
  mpi_enreg_diel%paral_fft=0
! Here, for TDDFT, artificially set iprcel . Also set a variable to reduce
! the memory needs.
  afford=1
  if(dtset%iscf==-1) then
   dtset%iprcel=21
   afford=0
  end if

! First compute dimensions
  if(dtset%iprcel>=21)then
!  With dielop=1, the matrices will be computed when istep=dielstrt
!  With dielop=2, the matrices will be computed when istep=dielstrt and 1
   dielop=1
   if(dtset%iprcel>=41)dielop=2
   if((dtset%iprcel >= 71).and.(dtset%iprcel<=79)) dielop=0 !RSkerker preconditioner do not need the susceptibility matrix
!  Immediate computation of dielectric matrix
   dielstrt=1
!  Or delayed computation
   if(modulo(dtset%iprcel,100)>21 .and. modulo(dtset%iprcel,100)<=29)dielstrt=modulo(dtset%iprcel,100)-20
   if(modulo(dtset%iprcel,100)>31 .and. modulo(dtset%iprcel,100)<=39)dielstrt=modulo(dtset%iprcel,100)-30
   if(modulo(dtset%iprcel,100)>41 .and. modulo(dtset%iprcel,100)<=49)dielstrt=modulo(dtset%iprcel,100)-40
   if(modulo(dtset%iprcel,100)>51 .and. modulo(dtset%iprcel,100)<=59)dielstrt=modulo(dtset%iprcel,100)-50
   if(modulo(dtset%iprcel,100)>61 .and. modulo(dtset%iprcel,100)<=69)dielstrt=modulo(dtset%iprcel,100)-60
!  Get diecut, and the fft grid to be used for the susceptibility computation
   diecut=abs(dtset%diecut)
   if( dtset%diecut<0.0_dp )then
    ecutsus=ecut
   else
    ecutsus= ( sqrt(ecut) *0.5_dp + sqrt(diecut) *0.25_dp )**2
   end if
!  Impose sequential calculation
   ngfftdiel(1:3)=0 ; ngfftdiel(7)=100 ; ngfftdiel(9)=0; ngfftdiel(8)=dtset%ngfft(8);ngfftdiel(10:18)=0
   if(dtset%iscf==-1)ngfftdiel(7)=102
   write(6,*) 'call getng diel'
   call getng(dtset%boxcutmin,ecutsus,gmet,mpi_enreg_diel%me_fft,mgfftdiel,nfftdiel,ngfftdiel,&
&   mpi_enreg_diel%nproc_fft,dtset%nsym,mpi_enreg_diel%fft_option_lob,mpi_enreg_diel%paral_fft,dtset%symrel)
!  Compute the size of the dielectric matrix
   kpt_diel(1:3)=(/ 0.0_dp, 0.0_dp, 0.0_dp /)
   call getmpw(diecut,dtset%exchn2n3d,gmet,(/1/),kpt_diel,&
&   mpi_enreg_diel,npwdiel,1,ucvol)
   lmax_diel=0
   if (psps%usepaw==1) then
    do ii=1,dtset%ntypat
     lmax_diel=max(lmax_diel,pawtab(itypat)%lcut_size)
    end do
   end if
  else
   npwdiel=1
   mgfftdiel=1
   nfftdiel=1
   lmax_diel=0
  end if

! Now, performs allocation
  allocate(dielinv(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
  allocate(susmat(2,npwdiel*afford,dtset%nspden,npwdiel,dtset%nspden))
  allocate(kg_diel(3,npwdiel))
  allocate(gbound_diel(2*mgfftdiel+8,2))
  allocate(irrzondiel(nfftdiel**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol))
  allocate(phnonsdiel(2,nfftdiel**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol))
  allocate(ph1ddiel(2,3*(2*mgfftdiel+1)*dtset%natom*psps%usepaw))
  allocate(ylmdiel(npwdiel,lmax_diel**2))
! Then, compute the values of different arrays
  if(dielop>=1)then
   call status(0,dtfil%filstat,iexit,level,'kpgio(sus)    ')
!  Note : kgnam is dummy, npwarr_diel is dummy, npwtot_diel is dummy
!  This kpgio call for going from the suscep FFT grid to the diel sphere
   npwarr_diel(1)=npwdiel
   call kpgio(diecut,dtset%exchn2n3d,gmet,(/1/),kg_diel,kgnam,&
&   kpt_diel,1,(/1/),1,'COLL',mpi_enreg_diel,npwdiel,&
&   npwarr_diel,npwtot_diel,dtset%nsppol,tmp_unit)
   call sphereboundary(gbound_diel,1,kg_diel,mgfftdiel,npwdiel)
   if (dtset%nsym>1 .and. dtset%iscf>0 ) then
    call setsym(densymop_diel,indsym,irrzondiel,dtset%iscf,dtset%natom,&
&    nfftdiel,ngfftdiel,dtset%nspden,dtset%nsppol,dtset%nsym,phnonsdiel,&
&    dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)
   end if
   if (psps%usepaw==1) then
    call getph(atindx,dtset%natom,ngfftdiel(1),ngfftdiel(2),ngfftdiel(3),ph1ddiel,xred)
    call initylmg(gprimd,kg_diel,kpt_diel,1,mpi_enreg_diel,lmax_diel,npwdiel,dtset%nband,1,&
&    npwarr_diel,dtset%nsppol,0,rprimd,tmp_unit,tmp_unit,ylmdiel,rhodum)
   end if
  end if

 else
  npwdiel=1
  mgfftdiel=1
  nfftdiel=1
 end if

!The big array containing functions defined on the fft grid is allocated now.
!Note however, that a zero value of mffmem will cause allocation with
!the third dimension set to zero ! In this case, another temporary
!will be used inside newvtr.
!This array is needed only in the self-consistent case
 if(nstep>0 .and. dtset%iscf>0) then
  if (psps%usepaw==0) npawmix=0
  if (psps%usepaw==1) allocate(f_paw(npawmix,n_fftgr*dtset%mffmem))
  if (psps%usepaw==1.and.dtset%pawmixdg==0) then
   ispmix=2;nfftmix=dtset%nfft;ngfftmix(:)=ngfft(:)
  else
   ispmix=1;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
  end if
  allocate(f_fftgr(ispmix*nfftmix,dtset%nspden,n_fftgr*dtset%mffmem))
 end if

 nkxc=0
!TDDFT - For a first coding
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==0)nkxc=0 !MF no xc kernel
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==1)nkxc=0 !MF no xc kern, but (later) RPA ok
 if (dtset%nfreqsus>0 .and. dtset%ikhxc==2)nkxc=1 !MF LDA xc kernel + (later) RPA
 if (dtset%iscf==-1 .and. dtset%nspden==1) nkxc=2
 if (dtset%iscf==-1 .and. dtset%nspden==2) nkxc=3
!Eventually need kxc-LDA when susceptibility matrix has to be computed
 if (dtset%iscf>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79)) nkxc=min(2*dtset%nspden-1,3)
!Eventually need kxc-LDA for residual forces (when density mixing is selected)
 if (dtset%iscf>=10.and.dtset%usewvl==0.and.forces_needed>0 .and. &
& abs(dtset%iprcch)>=1.and.abs(dtset%iprcch)<=6.and.abs(dtset%iprcch)/=5) then
  if (dtset%xclevel==1.or.dtset%iprcch>=0) nkxc=min(2*dtset%nspden-1,3)
  if (dtset%xclevel==2.and.dtset%nspden==2.and.dtset%iprcch<0) nkxc=23
 end if
 allocate(kxc(nfftf,nkxc))

!This flag will be set to 1 just before an eventual change of atomic
!positions inside the iteration, and set to zero when the consequences
!of this change are taken into account.
 moved_atm_inside=0
!This flag will be set to 1 if the forces are computed inside the iteration.
 computed_forces=0

 if(dtset%wfoptalg==2)then
  allocate(shiftvector((dtset%mband+2)*dtset%nkpt))
  val_min=-1.0_dp
  val_max=zero
 else
  allocate(shiftvector(1))
 end if

!Electric field initializations: initialize pel_cg(:) and p_ion(:)
 if (dtset%berryopt == 4) then
  unit_out=0;if (dtset%prtvol >= 10) unit_out=ab_out
  optberry=1     ! compute polarization only
  pel_cg(:) = zero;pelev_dum=zero
  call berryphase_new(cg,cprj,dtefield,dtfil,dtset,gprimd,hdr,kg,&
&  dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,nattyp,npwarr,nspinor,&
&  dtset%nsppol,psps%ntypat,dtset%nkpt,optberry,pawang,pawrad,pawtab,pel_cg,pelev_dum,pion,pwind,&
&  pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
&  unit_out,usecprj,psps%usepaw,wffnow,xred,psps%ziontypat)
 end if
!Positon calculation setup:
 if (dtset%positron/=0) then
  call setup_positron(dtset, dtfil%fildensin,dtfil%filstat,dtfil%filvhain,&
&  mpi_enreg,ngfftf(1), ngfftf(2),ngfftf(3),hdr,iexit,level,&
&  nfftf,psps%n1xccc,psps%ntypat,psps%xcccrc,psps%xccc1d,results_gs%etotal,&
&  rhocore,rhore,rhorp,rprimd,ucvol,vhae,vhap,xred,ngfft)
 end if

 if (dtset%iscf==22) energies%h0=zero
 call timab(54,2,tsec)

!######################################################################
!PERFORM ELECTRONIC ITERATIONS
!######################################################################

!Offer option of computing total energy with existing
!wavefunctions when nstep<=0, else do nstep iterations
!Note that for non-self-consistent calculations, this loop will be exited
!after the first call to vtorho
!Pass through the first routines even when nstep==0
 do istep=1,max(1,nstep)

  if (moved_atm_inside==1 .or. istep==1) then
!  ######################################################################
!  The following steps are done once for a given set of atomic
!  coordinates or for the nstep=1 case
!  ----------------------------------------------------------------------

!  Eventually symmetrize atomic coordinates over space group elements:
   call status(istep,dtfil%filstat,iexit,level,'call symzat   ')
   call symzat(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

   if (dtset%usewvl == 0) then
!   Get cut-off for g-vectors
    if (psps%usepaw==1) then
     write(message,'(2a)') ch10,' FFT (fine) grid used in SCF cycle:'
     call wrtout(6,message,'COLL')
    end if
    k0=zero
    call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,6,k0,ngfftf)

!   Compute structure factor phases and large sphere cut-off (gsqcut):
    call status(istep,dtfil%filstat,iexit,level,'call getph    ')
    call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
    if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
     call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
    else
     ph1df(:,:)=ph1d(:,:)
    end if
   end if

!  Initialization of atomic data for PAW
   if (psps%usepaw==1) then
!   Check for non-overlapping spheres
    call status(istep,dtfil%filstat,iexit,level,'call chkpawovlp')
    call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,xred)
!   Identify parts of the rectangular grid where the density has to be calculated
    optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
    if (forces_needed==1.or.(dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0)) then
     optgr1=dtset%pawstgylm;if (stress_needed==1) optrad=1
    end if
    call status(istep,dtfil%filstat,iexit,level,'call nhatgrid ')
    call nhatgrid(atindx1,gmet,mpi_enreg,dtset%natom,nattyp,nfftf,ngfftf,psps%ntypat,&
&    optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,dtset%typat,ucvol,xred)
   end if

!  If we are inside SCF cycle or inside dynamics over ions,
!  we have to translate the density of previous iteration
   moved_rhor=0
   if (initialized/=0.and.dtset%usewvl == 0.and. &
&   (abs(dtset%iprcch)==2.or.abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)) then
    moved_rhor=1
    if (abs(dtset%iprcch)==2) then
     option=2;allocate(workr(nfftf,dtset%nspden))
     call status(istep,dtfil%filstat,iexit,level,'call fresid   ')
     call fresid(dtset,gmet,gresid,gsqcut,mpi_enreg,nfftf,ngfftf,&
&     psps%ntypat,option,pawtab,rhor,rprimd,&
&     ucvol,workr,xred,xred_old,psps%znuclpsp)
     rhor=workr;deallocate(workr)
    else if (abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6) then
     call status(istep,dtfil%filstat,iexit,level,'call extraprho')
     call extraprho(atindx1,dtset,gmet,gprimd,gsqcut,mgfftf,mpi_enreg,&
&     psps%mqgrid_vl,nattyp,nfftf,ngfftf,psps%ntypat,pawrhoij,pawtab,&
&     ph1df,psps%qgrid_vl,rhor,rprimd,scf_history,ucvol,psps%usepaw,&
&     xred,xred_old,psps%ziontypat,psps%znuclpsp)
    end if
    call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
   end if

!  PAW only: we sometimes have to compute compensation density
!  and eventually add it to density from WFs
   nhatgrdim=0
   if (psps%usepaw==1.and.((usexcnhat==0) &
&   .or.(dtset%xclevel==2.and.(dtfil%ireadwf/=0.or.dtfil%ireadden/=0.or.initialized/=0)) &
&   .or.(dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0))) then
!   if (psps%usepaw==1.and. &
!   &   ((istep==1.and.dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0).or. &
!   (istep==1.and.dtfil%ireadden/=0.and.usexcnhat==0).or.&
!   (istep==1.and.initialized/=0.and.usexcnhat==0).or.&
!   (moved_atm_inside==1.and.usexcnhat==0))) then
    call timab(558,1,tsec)
    nhatgrdim=0;if (dtset%xclevel==2) nhatgrdim=usexcnhat*dtset%pawnhatxc
    ider=2*nhatgrdim
    if (nhatgrdim>0) allocate(nhatgr(nfftf,dtset%nspden,3))
    call pawmknhat(compch_fft,ider,0,mpi_enreg,dtset%natom,nfftf,ngfftf,nhatgrdim,dtset%nspden,&
&    psps%ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
    if (dtfil%ireadwf/=0.and.dtfil%ireadden==0.and.initialized==0) then
     rhor(:,:)=rhor(:,:)+nhat(:,:)
     call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
    end if
    call timab(558,2,tsec)
   end if

!  The following steps have been gathered in the setvtr routine:
!  - get Ewald energy and Ewald forces
!  - compute local ionic pseudopotential vpsp
!  - eventually compute 3D core electron density xccc3d
!  - eventually compute vxc and vhartr
!  - set up vtrial
   call status(istep,dtfil%filstat,iexit,level,'call setvtr   ')
   if (dtset%usewvl == 0) then
    optene = 4 * optres
    if(dtset%iscf==-3)optene=4
   else
!   We need the Hartree energy for the wavefunctions mixing
    optene = 1
   end if
   call setvtr(atindx1,dtset,energies,gmet,gprimd,grewtn,gsqcut,initialized,&
&   istep,kxc,mgfftf,moved_atm_inside,moved_rhor,mpi_enreg,&
&   nattyp,nfftf,ngfftf,nhat,nhatgr,nhatgrdim,nkxc,psps%ntypat,n1xccc,n3xccc,&
&   optene,pawtab,ph1df,psps,rhog,rhor,rmet,rprimd,strsxc,&
&   ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,&
&   xccc3d,xred,xred_old)

!  DEBUG
!  write(6,*)' scfcv : after setvtr, energies%e_hartree=',energies%e_hartree
!  ENDDEBUG

   if (nhatgrdim>0.and.nstep>0) deallocate(nhatgr)
!  Compute the electrons/positron correlation term
   if (dtset%positron/=0)then
    call calc_xc_ep(excapn,dtset%ixcpositron,dtset%positron,nfftf,dtset%nspden,0,rhocore,rhor,rhore,rhorp,&
&    rhotote,rhototp,rsepts,rsppts,vhae,vhap,vpsp,vtrial,vxcapn)
   end if

!  End the condition of atomic position change or istep==1
  end if

! ######################################################################
! The following steps are done at every iteration
! ----------------------------------------------------------------------

! PAW: Compute energies and potentials in the augmentation regions (spheres)
! Compute pseudopotential strengths (Dij quantities)
  if (psps%usepaw==1)then
!  "on-site" energies, potentials, densities computation
   nzlmopt=0;if (istep==2.and.dtset%pawnzlm>0) nzlmopt=-1
   if (istep>2) nzlmopt=dtset%pawnzlm
   option=0;if (dtset%iscf>0.and.dtset%iscf<10.and.nstep>0) option=1
   do iatom=1,dtset%natom
    itypat=dtset%typat(iatom)
    v_size=paw_an(iatom)%lm_size;if (dtset%pawxcdev==0) v_size=paw_an(iatom)%angl_size
    paw_ij(iatom)%has_dijhartree=1
    allocate(paw_ij(iatom)%dijhartree(pawtab(itypat)%lmn2_size))
    allocate(paw_an(iatom)%vxc1 (pawtab(itypat)%mesh_size,v_size,dtset%nspden))
    allocate(paw_an(iatom)%vxct1(pawtab(itypat)%mesh_size,v_size,dtset%nspden))
    if (dtset%pawspnorb>0) allocate(paw_an(iatom)%vh1(pawtab(itypat)%mesh_size,1,1))
    if (pawtab(itypat)%useexexch>0) allocate(paw_an(iatom)%vxc_ex(pawtab(itypat)%mesh_size,v_size,dtset%nspden))
   end do
   call status(istep,dtfil%filstat,iexit,level,'call pawdenpot')
   call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,dtset%ixc,dtset%natom,dtset%nspden,psps%ntypat,&
&   nzlmopt,option,paw_an,paw_ij,pawang,dtset%pawprtvol,&
&   pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%typat,dtset%xclevel,psps%znuclpsp)
!  PAW+U: impose density matrix if required
   if (dtset%usepawu>0) then
    impose_dmat=0
    if ((istep<=abs(dtset%usedmatpu)).and.(dtset%usedmatpu<0.or.initialized0==0)) impose_dmat=1
    if (impose_dmat==1.or.dtset%dmatudiag/=0) then
     dimdmat=0;if (impose_dmat==1) dimdmat=2*lpawumax+1
     call setnoccmmp(0,dimdmat,dtset%dmatpawu(1:dimdmat,1:dimdmat,1:dtset%nsppol*nspinor,1:dtset%natpawu*impose_dmat),&
&         dtset%dmatudiag,impose_dmat,indsym,dtset%natom,dtset%natpawu,&
&     dtset%nspden,nspinor,dtset%nsppol,dtset%nsym,dtset%ntypat,paw_ij,pawang,dtset%pawprtvol,&
&     pawrhoij,pawtab,dtset%spinat,dtset%symafm,dtset%typat,0,dtset%usepawu)
    end if
   end if
!  Dij computation
   call status(istep,dtfil%filstat,iexit,level,'call pawdij   ')
   call pawdij(dtset,dtset%enunit,mpi_enreg,dtset%natom,nfftf,ngfftf,dtset%nspden,psps%ntypat,&
&   paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev,&
&   dtset%typat,ucvol,vtrial,vxc)
   do iatom=1,dtset%natom
    deallocate(paw_ij(iatom)%dijhartree,paw_an(iatom)%vxc1,paw_an(iatom)%vxct1)
    paw_ij(iatom)%has_dijhartree=0
    if (dtset%pawspnorb>0) deallocate(paw_an(iatom)%vh1)
    if (pawtab(itypat)%useexexch>0) deallocate(paw_an(iatom)%vxc_ex)
   end do
   call status(istep,dtfil%filstat,iexit,level,'call symdij   ')
   call symdij(psps%indlmn,indsym,psps%lmnmax,dtset%natom,dtset%nsym,psps%ntypat,paw_ij,pawang,&
&   dtset%pawprtvol,dtset%symafm,symrec,dtset%typat)
  end if

! No need to continue and call vtorho, when nstep==0
  if(nstep==0)exit

! ######################################################################
! The following steps are done only when nstep>0
! ----------------------------------------------------------------------
  call timab(56,1,tsec)
  call status(istep,dtfil%filstat,iexit,level,'loop istep    ')

  if(dtset%iscf>0)then
   write(message, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
   call wrtout(06,message,'COLL')
  end if

! The next flag says whether the xred have to be changed in the current iteration
  moved_atm_inside=0
  if(dtset%ionmov==4 .and. mod(iapp,2)/=1 .and. dtset%iscf>0 )moved_atm_inside=1
  if(dtset%ionmov==5 .and. iapp/=1 .and. istep==1 .and. dtset%iscf>0)moved_atm_inside=1

! The next flag says whether the forces have to be computed in the current iteration
  computed_forces=0
  if ((dtset%optforces==1 .and. dtset%usewvl == 0).or.(moved_atm_inside==1)) computed_forces=1
  if (abs(tollist(3))>tiny(0._dp)) computed_forces=1
  if (dtset%iscf<0.or.dtset%positron==1) computed_forces=0
  if ((istep==1).and.(dtset%optforces/=1)) then
   if (moved_atm_inside==1) then
    write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&    ' scfcv : WARNING -',ch10,&
&    '  Although the computation of forces during electronic iterations',ch10,&
&    '  was not required by user, it is done (required by the',ch10,&
&    '  choice of ionmov input parameter).'
    call wrtout(6,message,'COLL')
   end if
   if (abs(tollist(3))+abs(tollist(7))>tiny(0._dp)) then
    write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&    ' scfcv : WARNING -',ch10,&
&    '  Although the computation of forces during electronic iterations',ch10,&
&    '  was not required by user, it is done (required by the',ch10,&
&    '  "toldff" or "tolrff" tolerance criteria).'
    call wrtout(6,message,'COLL')
   end if
  end if
  if ((istep==1).and.(dtset%optforces==1).and. dtset%usewvl == 1) then
   write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&   ' scfcv : WARNING -',ch10,&
&   '  Although the computation of forces during electronic iterations',ch10,&
&   '  was required by user, it has been disable since the tolerence',ch10,&
&   '  is not on forces (force computation is expensive in wavelets).'
   call wrtout(6,message,'COLL')
  end if

  call timab(56,2,tsec)

! ######################################################################
! Compute the density rho from the trial potential
! ----------------------------------------------------------------------

! Compute the density from the trial potential
  if (dtset%tfkinfunc==0) then
   call status(istep,dtfil%filstat,iexit,level,'call vtorho   ')
   call vtorho(afford,atindx,atindx1,cg,compch_fft,cpus,dbl_nnsclo,&
&   densymop_diel,densymop_gs,dielop,dielstrt,dphase,dtefield,dtfil,dtset,&
&   eigen,energies,etotal,filapp,gbound_diel,&
&   gmet,gprimd,grnl,gsqcut,hdr,indsym,irrzon,irrzondiel,&
&   istep,kg,kg_diel,kxc,lmax_diel,mgfftdiel,mpi_enreg,&
&   psps%mpsang,dtset%natom,nattyp,nfftf,nfftdiel,ngfftdiel,nhat,nkxc,&
&   npwarr,npwdiel,res2,nspinor,psps%ntypat,nvresid,occ,computed_forces,&
&   optres,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,phnons,phnonsdiel,ph1d,ph1ddiel,psps,&
&   pwind,pwind_alloc,pwnsfac,resid,residm,&
&   rhog,rhor,rmet,rprimd,shiftvector,susmat,symrec,&
&   ucvol,wffnew,wffnow,val_min,val_max,vtrial,wvl,xred,ylm,ylmdiel)
   call status(istep,dtfil%filstat,iexit,level,'after vtorho  ')
  elseif (dtset%tfkinfunc==1) then
   write(6,*)'WARNING : THOMAS FERMI'
   call vtorhotf(densymop_gs,dtfil,dtset,energies%e_kinetic,energies%e_nonlocalpsp,energies%entropy,energies%e_fermie,&
&   grnl,irrzon,mpi_enreg,dtset%natom,nfftf,dtset%nspden,dtset%nsppol,dtset%nsym,phnons,&
&   rhog,rhor,ucvol,vtrial)
   residm=zero
   energies%e_eigenvalues=zero
  else
   write(6,*)'WARNING : RECURSION'
   if(istep==1.and.dtset%prtvol/=-7)then
    quitrec=0
    get_ek=0
    get_entropy=0
   end if
   if(istep == nstep.or.quitrec /= 0)then
    get_ek=1
    get_entropy=1
   end if
   call vtorhorec(densymop_gs,dtfil,dtset,&
&   energies%e_kinetic,energies%e_nonlocalpsp,energies%entropy,energies%e_eigenvalues,energies%e_fermie,&
&   grnl,irrzon,mpi_enreg,dtset%natom,nfftf,dtset%nspden,dtset%nsppol,dtset%nsym,phnons,&
&   rhog,rhor,ucvol,vtrial,rmet,quitrec,get_ek,get_entropy)
   residm=zero
   call symrhg(1,densymop_gs,irrzon,mpi_enreg,nfftf,&
&   dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3),dtset%ngfft,dtset%nspden,&
&   dtset%nsppol,dtset%nsym,dtset%paral_kgb,&
&   phnons,rhog,rhor,dtset%symafm)
  end if
  if(dtset%wfoptalg==2)then
   do ikpt=1,dtset%nkpt
    shiftvector(1+(ikpt-1)*(dtset%mband+2))=val_min
    shiftvector(2+(ikpt-1)*(dtset%mband+2):ikpt*(dtset%mband+2)-1)=&
&    eigen((ikpt-1)*dtset%mband+1:ikpt*dtset%mband)
    shiftvector(ikpt*(dtset%mband+2))=val_max
   end do
  end if

! ######################################################################
! Skip out of step loop if non-SCF (completed)
! ----------------------------------------------------------------------

! Indeed, nstep loops have been done inside vtorho
  if (dtset%iscf<=0) exit

! ######################################################################
! In case of density mixing or wvlet handling, compute the total energy
! ----------------------------------------------------------------------
  if (dtset%iscf>=10 .or. dtset%usewvl == 1) then
   if (dtset%usewvl == 0) then
    optene = 1 ! use double counting scheme
   else if (dtset%iscf/=22) then
    optene = 0 ! use direct scheme for computation of energy
   else
    optene = -1
   end if

!  if the mixing is the ODA mixing, compute energy and new density here
   if (dtset%iscf==22) then
    call odamix(atindx,deltae,dtset,dtefield%efield_dot,elast,energies,&
&    etotal,gsqcut,indsym,kxc,mgfftf,mpi_enreg,nattyp,&
&    nfftf,ngfftf,nhat,nkxc,psps%ntypat,nvresid,n1xccc,n3xccc,optene,optres,&
&    paw_ij,paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,pel_cg,ph1df,pion,psps,rhog,rhor,rprimd,strsxc,symrec,&
&    ucvol,psps%usepaw,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred)
   end if
!  If the density mixing is required, compute the total energy here
   call etotfor(atindx1,deltae,diffor,dtset,dtefield%efield_dot,elast,energies,&
&   etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,&
&   grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,nattyp,&
&   nfftf,ngfftf,nhat,nkxc,psps%ntypat,&
&   nvresid,n1xccc,n3xccc,optene,computed_forces,optres,&
&   pawang,pawfgrtab,pawrhoij,pawtab,pel_cg,ph1df,pion,psps,rhog,rhor,rprimd,symrec,synlgr,&
&   ucvol,psps%usepaw,usexcnhat,vhartr,vpsp,vxc,xccc3d,xred)
  end if

! ######################################################################
! In case of density mixing, check the exit criterion
! ----------------------------------------------------------------------
  if (dtset%iscf>=10 .and. dtset%usewvl == 0) then
!  Check exit criteria
   call timab(52,1,tsec)
   call status(istep,dtfil%filstat,iexit,level,'call scprqt   ')
   choice=2
   call scprqt(choice,cpus,deltae,diffor,dtset,&
&   eigen,etotal,favg,fcart,energies%e_fermie,filapp,dtfil%filnam_ds(1),&
&   initialized0,dtset%iscf,istep,dtset%kptns,maxfor,moved_atm_inside,mpi_enreg,&
&   dtset%nband,dtset%nkpt,nstep,occ,optres,&
&   prtfor,quit,res2,resid,residm,response,tollist,psps%usepaw,&
&   vxcavg,dtset%wtk,xred)

   if(dtset%tfkinfunc==2)then  !exit criteria for the recursion method
    if(quitrec==2)quit=1
   end if

   if (istep==nstep) quit=1
   call timab(52,2,tsec)

!  If criteria in scprqt say to quit, then exit the loop over istep.
   if(mpi_enreg%paral_compil_kpt==1)then
    if (mpi_enreg%parareel == 0) then
     quit_sum=quit

!    BEGIN TF_CHANGES
     call xcomm_world(mpi_enreg,spaceComm)
!    END TF_CHANGES

     call xsum_mpi(quit_sum,spaceComm,ierr)

     if (quit_sum > 0) exit
    end if
   end if ! mpi_enreg%paral_compil_kpt==1
   if (quit==1) exit
  end if

! ######################################################################
! Mix the total density (if required)
! ----------------------------------------------------------------------
  if (dtset%iscf>=10 .and.dtset%iscf/=22.and. dtset%usewvl == 0) then
!  If LDA dielectric matrix is used for preconditionning, has to update here Kxc
   if (nkxc>0.and.modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79) &
&     .and.((istep==1.or.istep==dielstrt).or.(dtset%iprcel>=100))) then
    optxc=10
    call rhohxc(dtset,edum,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,nhat,&
&    psps%usepaw,nhatgr,0,nkxc,dtset%nspden,n3xccc,optxc,rhog,&
&    rhor,rprimd,dummy2,0,vhartr,vxc,vxcavg_dum,xccc3d)
   end if
   call status(istep,dtfil%filstat,iexit,level,'call newrho   ')
   call newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,dtset,etotal,fcart,pawfgr%fintocoa,filfft,&
&   f_atm,f_fftgr,f_paw,gmet,grhf,gsqcut,initialized,&
&   ispmix,istep,i_rhor,i_vresid,i_vrespc,i_vtrial,kg_diel,kxc,mgfftf,mgfftdiel,pawfgr%coatofin,&
&   moved_atm_inside,mpi_enreg,nattyp,nfftf,nfftmix,ngfftf,ngfftmix,nkxc,npawmix,npwdiel,&
&   nvresid,psps%ntypat,n_fftgr,n_index,n1xccc,pawrhoij,pawtab,&
&   ph1df,psps,rhog,rhor,rprimd,susmat,psps%usepaw,vtrial,xred)
  end if   ! iscf>=10

! ######################################################################
! Additional computation in case of an electric field
! ----------------------------------------------------------------------

! In case of an electric field calculation, need polarization
! to compute electric enthalpy instead of energy.

  if (psps%usepaw==0.and.dtset%berryopt == 4) then
!  When using symmetry, it is costly to update polarization from changes in Zak
!  phases. It is better to call berryphase here.
!  Update polarization by adding increment from the SCF step
!  pel_cg(1:3) = pel_cg(1:3) + dtefield%sdeg*dphase(1:3)/two_pi
!  ptot(1:3) = pel_cg(1:3) + pion(1:3)
!  write(message,'(6(a),3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
!  &    ' scfcv: Polarization from accumulated change in Berry phase:',ch10,&
!  &    ' (reduced coordinates, a. u., without correcting for branch cuts)',ch10,&
!  &    '     Electronic: ', (pel_cg(ii), ii = 1, 3), ch10,&
!  &    '     Ionic:      ', (pion(ii), ii = 1, 3), ch10, &
!  &    '     Total:      ', (ptot(ii), ii = 1, 3)
!  call wrtout(06,message,'COLL')
   pelev_dum=zero
   call berryphase_new(cg,cprj,dtefield,dtfil,dtset,&
&   gprimd,hdr,kg,&
&   dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%natom,nattyp,npwarr,nspinor,&
&   dtset%nsppol,psps%ntypat,dtset%nkpt,optberry,pawang,pawrad,pawtab,pel_cg,pelev_dum,pion,pwind,&
&   pwind_alloc,pwnsfac,rprimd,dtset%typat,ucvol,&
&   unit_out,usecprj,psps%usepaw,wffnow,xred,psps%ziontypat)
   ptot(:) = pel_cg(:) + pion(:)
   write(message,'(6(a),3(e16.9,2x),a,a,3(e16.9,2x),a,a,3(e16.9,2x))')ch10,&
&   ' scfcv: New value of the polarization:',ch10,&
&   ' (reduced coordinates, a. u.)',ch10,&
&   '     Electronic: ', (pel_cg(ii), ii = 1, 3), ch10,&
&   '     Ionic:      ', (pion(ii), ii = 1, 3), ch10, &
&   '     Total:      ', (ptot(ii), ii = 1, 3)
   call wrtout(06,message,'COLL')
  end if       ! berryopt


! ######################################################################
! Compute the new potential from the trial density
! ----------------------------------------------------------------------

! Set XC computation flag
  optxc=1
  if (nkxc>0) then
   if (dtset%nfreqsus>0) optxc=2
   if (dtset%iscf<0) optxc=2
   if (modulo(dtset%iprcel,100)>=61.and.(dtset%iprcel<71.or.dtset%iprcel>79).and. &
&      dtset%iscf<10.and. &
&   (dtset%iprcel>=100.or.istep==1.or.istep==dielstrt)) optxc=2
   if (dtset%iscf>=10.and.dtset%iprcch/=0.and.abs(dtset%iprcch)/=5) optxc=2
   if (optxc==2.and.dtset%xclevel==2.and.nkxc==3-2*mod(dtset%nspden,2)) optxc=12
  end if
  if (dtset%positron==1) optxc=-1

! Check for positron case
  if (dtset%positron==1) then
   vpsp(:)=-vpsp(:);vhartr(:)=-vhae(:,1)
   call calc_xc_ep(excapn,dtset%ixcpositron,dtset%positron,nfftf,dtset%nspden,1,rhocore,rhor,rhore,rhorp,&
&   rhotote,rhototp,rsepts,rsppts,vhae,vhap,vpsp,vtrial,vxcapn)
   vxc(:,1) = vxcapn(:)
  end if
  if (dtset%iscf/=22) then
!  PAW: eventually recompute compensation density (and gradients)
   nhatgrdim=0
   if (psps%usepaw==1) then
    ider=-1;if (dtset%iscf>=10.and.((dtset%xclevel==2.and.dtset%pawnhatxc>0).or.usexcnhat==0)) ider=0
    if (dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0) ider=ider+2
    if (ider>0) then
     nhatgrdim=1;allocate(nhatgr(nfftf,dtset%nspden,3))
    end if
    call timab(558,1,tsec)
    if (ider>=0) then
     call timab(558,1,tsec)
     call pawmknhat(compch_fft,ider,0,mpi_enreg,dtset%natom,nfftf,ngfftf,nhatgrdim,dtset%nspden,&
&     psps%ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
     call timab(558,2,tsec)
    end if
   end if
!  Compute new potential from the trial density
   call status(istep,dtfil%filstat,iexit,level,'call rhotov')
   optene=2*optres;if(psps%usepaw==1) optene=2

   call rhotov(dtset,energies,gsqcut,kxc,mpi_enreg,nfftf,ngfftf, &
&   nhat,nhatgr,nhatgrdim,nkxc,nvresid,n3xccc,&
&   optene,optres,optxc,pawang,pawfgrtab,pawtab,&
&   rhog,rhor,rprimd,strsxc,ucvol,psps%usepaw,usexcnhat,&
&   vhartr,vnew_mean,vpsp,vres_mean,res2,vtrial,vxcavg,vxc,xccc3d)
  end if
! Check for positron case
  if (dtset%positron==1) vpsp(:)=-vpsp(:)
  if (dtset%positron==2) then
   vhartr(:)=vhartr(:)-vhap(:,1)
   call calc_xc_ep(excapn,dtset%ixcpositron,dtset%positron,nfftf,dtset%nspden,1,rhocore,rhor,rhore,rhorp,&
&   rhotote,rhototp,rsepts,rsppts,vhae,vhap,vpsp,vtrial,vxcapn)
   vxc(:,1)=vxc(:,1)+vxcapn(:)
  end if

! If the xred have to be changed in the current iteration, they has to be saved
  if(moved_atm_inside==1) xred_old(:,:)=xred(:,:)

! ######################################################################
! Check exit criteria in case of potential mixing
! ----------------------------------------------------------------------
  if (dtset%iscf<10 .and. dtset%usewvl == 0) then

!  If the potential mixing is required, compute the total energy here
!  PAW: has to compute here spherical terms
   if (psps%usepaw==1) then
    nzlmopt=0;if (istep==1.and.dtset%pawnzlm>0) nzlmopt=-1
    if (istep>1) nzlmopt=dtset%pawnzlm
    do iatom=1,dtset%natom
     allocate(paw_ij(iatom)%dijhartree(pawtab(dtset%typat(iatom))%lmn2_size))
     paw_ij(iatom)%has_dijhartree=1
    end do
    call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,dtset%ixc,dtset%natom,dtset%nspden,psps%ntypat,&
&    nzlmopt,2,paw_an,paw_ij,pawang,dtset%pawprtvol,&
&    pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%typat,dtset%xclevel,psps%znuclpsp)
    do iatom=1,dtset%natom
     deallocate(paw_ij(iatom)%dijhartree)
     paw_ij(iatom)%has_dijhartree=0
    end do
   end if

   call status(istep,dtfil%filstat,iexit,level,'call etotfor  ')
   call etotfor(atindx1,deltae,diffor,dtset,&
&   dtefield%efield_dot,elast,energies,&
&   etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,&
&   grxc,gsqcut,indsym,kxc,maxfor,mgfftf,mpi_enreg,nattyp,&
&   nfftf,ngfftf,nhat,nkxc,dtset%ntypat,nvresid,n1xccc, &
&   n3xccc,0,computed_forces,optres,&
&   pawang,pawfgrtab,pawrhoij,pawtab,pel_cg,ph1df,pion,psps,rhog,rhor,rprimd,symrec,synlgr,&
&   ucvol,psps%usepaw,usexcnhat,vhartr,vpsp,vxc,xccc3d,xred)
  end if

! ######################################################################
! Check exit criteria in case of potential mixing or wavelet handling
! ----------------------------------------------------------------------
  if (dtset%iscf<10 .or. dtset%usewvl == 1) then
!  Check exit criteria
   call timab(52,1,tsec)
   call status(istep,dtfil%filstat,iexit,level,'call scprqt   ')
   choice=2
   call scprqt(choice,cpus,deltae,diffor,dtset,&
&   eigen,etotal,favg,fcart,energies%e_fermie,filapp,dtfil%filnam_ds(1),&
&   initialized0,dtset%iscf,istep,dtset%kptns,maxfor,moved_atm_inside,mpi_enreg,&
&   dtset%nband,dtset%nkpt,nstep,occ,optres,&
&   prtfor,quit,res2,resid,residm,response,tollist,psps%usepaw,&
&   vxcavg,dtset%wtk,xred)
   if (istep==nstep.and.psps%usepaw==1) quit=1
   call timab(52,2,tsec)

!  exit criteria for the recursion
   if(dtset%tfkinfunc==2)then
    if(quitrec==2)quit=1
   end if

!  If criteria in scprqt say to quit, then exit the loop over istep.
   if(mpi_enreg%paral_compil_kpt==1)then
    if (mpi_enreg%parareel == 0) then
     quit_sum=quit

!    BEGIN TF_CHANGES
     call xcomm_world(mpi_enreg,spaceComm)
!    END TF_CHANGES

     call xsum_mpi(quit_sum,spaceComm,ierr)

     if (quit_sum > 0) exit
    end if
   end if ! mpi_enreg%paral_compil_kpt==1
   if (quit==1) then
    do ispden=1,dtset%nspden
     vtrial(:,ispden)=vtrial(:,ispden)+nvresid(:,ispden)+vres_mean(ispden)
    end do
    exit
   end if
  end if

! ######################################################################
! Mix the potential (if required) - Check exit criteria
! ----------------------------------------------------------------------
  if (dtset%iscf<10 .and. dtset%usewvl /= 1) then
!  Precondition the residual and forces, then determine the new vtrial
!  (Warning: the (H)xc potential may have been subtracted from vtrial)
   call status(istep,dtfil%filstat,iexit,level,'call newvtr   ')
   call newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&
&   dtn_pc,dtset,energies%e_fermie,etotal,fcart,pawfgr%fintocoa,filfft,&
&   f_atm,f_fftgr,f_paw,gmet,grhf,gsqcut,initialized,ispmix,&
&   istep,i_rhor,i_vresid,i_vrespc,i_vtrial,&
&   kg_diel,kxc,mgfftf,mgfftdiel,pawfgr%coatofin,&
&   moved_atm_inside,mpi_enreg,nattyp,nfftf,nfftmix,&
&   nhat,nhatgr,nhatgrdim,&
&   ngfftf,ngfftmix,nkxc,npawmix,npwdiel,&
&   nstep,psps%ntypat,n_fftgr,n_index,n1xccc,optres,optxc,&
&   pawrhoij,pawang,pawfgrtab, &
&   ph1df,&
&   psps,rhor,rprimd,susmat,psps%usepaw,&
&   vhartr,vnew_mean,vpsp,nvresid,&
&   vtrial,vxc,xred,&
&   atindx1,cg,deltae,densymop_gs,&
&   dtfil,energies%e_eigenvalues,energies%e_ewald,eigen,energies%e_corepsp,&
&   energies%e_kinetic,energies%e_nonlocalpsp,energies%entropy,&
&   energies%e_paw,energies%e_pawdc,irrzon,kg,&
&   nfftf,&
&   ngfftf,npwarr,n3xccc,occ,optene,&
&   pawfgr,pawtab,phnons,&
&   resid,rhog,&
&   usexcnhat,&
&   wffnow,&
&   ylm,nspinor,xccc3d)
  end if   ! iscf<10

! No potential mixing in wavelet, direct minimisation scheme!
  if (dtset%usewvl == 1) then
   call status(istep,dtfil%filstat,iexit,level,'call wvl_newvtr')
   call wvl_newvtr(dtset, mpi_enreg, nele, offset, vhartr, vpsp, vtrial, vxc)
  end if


! ######################################################################
! END MINIMIZATION ITERATIONS
! ######################################################################

! The initialisation of the gstate run should be done when this point is reached
  initialized=1

! This is to save the density for restart.
  if( mpi_enreg%paral_compil_kpt==0                         .or. &
&  (mpi_enreg%me==0 .and. mpi_enreg%parareel == 0)        .or. &
&  (mpi_enreg%me_group_para==0 .and. mpi_enreg%parareel == 1)) then
   prtden=dtset%prtden
   if (prtden<0) then
    if (mod(istep-1,abs(prtden))==0) then
     isave=isave+1
     call status(0,dtfil%filstat,iexit,level,'call ioarr-den')
     rdwr=2 ; fformr=52 ; rdwrpaw=0
     call int2char4(mod(isave,2),tag)
     fildata=trim(filprot)//'_DEN_'//tag
     accessfil = 0
     call ioarr(accessfil,rhor, dtset, etotal,fformr,fildata,hdr, mpi_enreg, &
&     nfftf,pawrhoij,rdwr,rdwrpaw,ngfftf)
    end if
   end if
  end if

  if (nhatgrdim>0) deallocate(nhatgr)

 end do ! istep

 if (quit==1.and.nstep==1) initialized=1

!######################################################################
!Case nstep==0: compute energy based on incoming wf
!----------------------------------------------------------------------

 if(nstep==0) then
  optene=2*psps%usepaw+optres
  energies%entropy=results_gs%energies%entropy  !MT20070219: entropy is not recomputed in routine energy
  call energy(atindx,atindx1,cg,compch_fft,densymop_gs,dtfil,dtset,energies,eigen,&
&  etotal,gsqcut,indsym,irrzon,kg,mpi_enreg,nattyp,nfftf,ngfft,ngfftf,nhat,nhatgr,nhatgrdim,&
&  npwarr,nspinor,n3xccc,occ,optene,paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,&
&  phnons,ph1d,psps,resid,rhog,rhor,rprimd,strsxc,symrec,usexcnhat,vhartr,vtrial,&
&  vpsp,vxc,wffnow,wvl%wfs,xccc3d,xred,ylm)
  if (nhatgrdim>0) deallocate(nhatgr)
 end if ! nstep==0

!######################################################################
!Additional steps after SC iterations, including force, stress, polarization calculation
!----------------------------------------------------------------------

 call timab(60,1,tsec)
 call status(0,dtfil%filstat,iexit,level,'endloop istep ')

!PAW: in some cases, need to recompute <p_lmn|Cnk> projected WF:
!should be output from vtorho (but is "type-sorted" in votorho, not here)...
 if (psps%usepaw==1.and. &
& (dtset%prtwant==2.or.dtset%prtnabla>0.or.dtset%prtdos==3.or.dtset%berryopt/=0.or.dtset%kssform==3)) then
  usecprj=1
  allocate(dimcprj(dtset%natom),cprj(dtset%natom,nspinor*dtset%mband*dtset%mkmem*dtset%nsppol))
  do iatom=1,dtset%natom
   dimcprj(iatom)=pawtab(dtset%typat(iatom))%lmn_size
  end do
  call cprj_alloc(cprj,0,dimcprj)
  call ctocprj(atindx,cg,1,cprj,dtfil,gmet,gprimd,0,0,1,dtset%istwfk,kg,dtset%kptns,&
&  dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,&
&  dtset%natom,nattyp,dtset%nband,dtset%natom,ngfft,dtset%nkpt,dtset%nloalg,&
&  npwarr,nspinor,dtset%nsppol,dtset%ntypat,ph1d,psps,rmet,dtset%typat,&
&  ucvol,dtfil%unpaw,wffnow,xred,ylm)
 end if

!SHOULD CLEAN THE ARGS OF THIS ROUTINE
 call status(0,dtfil%filstat,iexit,level,'afterscfloop  ')
 call afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,dimcprj,&
& deltae,diffor,dtefield,dtfil,dtset,eigen,energies,etotal,&
& favg,fcart,filapp,filfft,forold,fred,gresid,grewtn,grhf,&
& grxc,gsqcut,hdr,indsym,&
& istep,kg,kxc,maxfor,dtset%mgfft,mgfftf,&
& moved_atm_inside,mpi_enreg,&
& n3xccc,nattyp,&
& nfftf,ngfft,ngfftf,nhat,nkxc,npwarr,nvresid,&
& occ,optres,optxc,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,pel,pel_cg,&
& ph1d,ph1df,pion,prtfor,psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&
& rhocore,rhog,rhor,rhore,rhototp,&
& rprimd,stress_needed,strsxc,strten,symrec,synlgr,tollist,usecprj,usexcnhat,&
& vhartr,vpsp,vxc,vxcavg,wffnow,wvl,xccc3d,xred,xred_old,ylm,ylmgr)

 call timab(60,2,tsec)

!######################################################################
!All calculations in scfcv are finished. Printing section
!----------------------------------------------------------------------

 call status(istep,dtfil%filstat,iexit,level,'call outscfcv ')
 call outscfcv(atindx,atindx1,cg,compch_fft,compch_sph,cprj,dimcprj,dtfil,dtset,ecut,eigen,etotal,&
& energies%e_fermie,filapp,gmet,gprimd,gsqcut,hdr,kg,&
& kssform,dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,&
& nattyp,nfftf,ngfftf,nhat,dtset%nkpt,npwarr,dtset%nspden,nspinor,dtset%nsppol,&
& dtset%nsym,psps%ntypat,n3xccc,occ,&
& pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,paw_ij,ph1d,dtset%prtvol,psps,rhog,rhor,rmet,rprimd,&
& ucvol,usecprj,usexcnhat,wffnow,vhartr,vtrial,vxc,vxcavg,xccc3d,xred,ylm)

 if(mpi_enreg%paral_compil_kpt==1)then
  call timab(61,1,tsec)
  if (mpi_enreg%parareel == 0) then

!  BEGIN TF_CHANGES
   call leave_test(mpi_enreg)
!  END TF_CHANGES

  end if
  call timab(61,2,tsec)
 end if

 call timab(60,1,tsec)

!Transfer eigenvalues computed by BigDFT in afterscfloop to eigen.
 if (dtset%usewvl == 1) then
  eigen = wvl%wfs%eval
 end if

!Debugging : print the different parts of rhor, as well as vxc
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(dtset%prtvol==-level)then
  write(message,'(a)') '   ir     vxc(ir)     rhor(ir)     '
  call wrtout(06,message,'COLL')
  do ir=1,nfftf
   if(ir<=11 .or. mod(ir,301)==0 )then
    write(message,'(i5,a,2es13.6)')ir,' ',vxc(ir,1),rhor(ir,1)
    call wrtout(06,message,'COLL')
    if(dtset%nspden==2)then
     write(message,'(a,2es13.6)')'      ',vxc(ir,2),rhor(ir,2)
     call wrtout(06,message,'COLL')
    end if
   end if
  end do
 end if

!Structured debugging : if prtvol=-level, stop here.
 if(dtset%prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)') ch10,' scfcv : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!######################################################################
!Deallocate memory and save results
!----------------------------------------------------------------------

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')
 if (psps%usepaw==1.and. &
& (dtset%prtwant==2.or.dtset%prtnabla>0.or.dtset%prtdos==3.or.dtset%berryopt/=0.or.dtset%kssform==3)) then
  deallocate(dimcprj)
  usecprj=0
  call cprj_free(cprj)
  deallocate(cprj)
 end if
 deallocate(fcart,fred,forold)
 deallocate(grnl,gresid,grewtn,grxc,synlgr)
 deallocate(ph1d,ph1df)
 deallocate(vhartr,vtrial,vpsp,vxc,xccc3d)
 deallocate(kxc,shiftvector)
 if (dtset%iscf>0) then
  deallocate(dtn_pc,f_atm,grhf,nvresid)
  if(nstep>0) deallocate(f_fftgr)
  if(nstep>0.and.psps%usepaw==1) deallocate(f_paw)
  deallocate(i_rhor,i_vtrial,i_vresid,i_vrespc)
 end if
 if((nstep>0.and.dtset%iscf>0).or.dtset%iscf==-1) then
  deallocate(dielinv,gbound_diel)
  deallocate(irrzondiel,kg_diel)
  deallocate(phnonsdiel,susmat)
  deallocate(ph1ddiel,ylmdiel)
 end if
 if (psps%usepaw==1) then
  deallocate(nhat)
  do iatom=1,dtset%natom
   itypat=dtset%typat(iatom)
   if (associated(pawfgrtab(iatom)%ifftsph))deallocate(pawfgrtab(iatom)%ifftsph)
   if (associated(pawfgrtab(iatom)%rfgd))    deallocate(pawfgrtab(iatom)%rfgd)
   if (associated(pawfgrtab(iatom)%gylm))    deallocate(pawfgrtab(iatom)%gylm)
   if (associated(pawfgrtab(iatom)%gylmgr))  deallocate(pawfgrtab(iatom)%gylmgr)
   if (associated(pawfgrtab(iatom)%gylmgr2))deallocate(pawfgrtab(iatom)%gylmgr2)
   pawfgrtab(iatom)%nfgd=0;pawfgrtab(iatom)%rfgd_allocated=0
   pawfgrtab(iatom)%gylm_allocated=0;pawfgrtab(iatom)%gylmgr_allocated=0
   pawfgrtab(iatom)%gylmgr2_allocated=0;pawfgrtab(iatom)%vlocgr_allocated=0
   deallocate(paw_an(iatom)%lmselect,paw_ij(iatom)%dij)
   if (paw_an(iatom)%has_vxcval==1) deallocate(paw_an(iatom)%vxc1_val,paw_an(iatom)%vxct1_val)
   if (paw_ij(iatom)%has_dijhat>0) deallocate(paw_ij(iatom)%dijhat)
   if (paw_ij(iatom)%has_dijxc>0) deallocate(paw_ij(iatom)%dijxc)
   if (paw_ij(iatom)%has_dijxc_val>0) deallocate(paw_ij(iatom)%dijxc_val)
   if (paw_ij(iatom)%has_dijso>0) deallocate(paw_ij(iatom)%dijso)
   if (paw_ij(iatom)%has_dijU>0) deallocate(paw_ij(iatom)%dijU)
   if (pawtab(itypat)%usepawu>0) deallocate(paw_ij(iatom)%noccmmp,paw_ij(iatom)%nocctot)
   if (pawtab(itypat)%useexexch>0) deallocate(paw_ij(iatom)%vpawx)
   deallocate(pawrhoij(iatom)%rhoijres);pawrhoij(iatom)%use_rhoijres=0
   if (dtset%iscf>0) then
    pawrhoij(iatom)%lmnmix_sz=0
    deallocate(pawrhoij(iatom)%kpawmix)
   end if
  end do
  deallocate(pawfgrtab,paw_an,paw_ij)
 end if
 if (dtset%positron>0) then
  deallocate(vhae,vhap,rhorp)
  deallocate(rhotote,rhototp,rsepts,vxcapn,rspts,rhore,excapn,rsppts,rhocore)
 end if

!Restore some variables in the dtset
!Here, for TDDFT, iprcel was artificially set.
 if(dtset%iscf==-1) then
  dtset%iprcel = dtset_iprcel
 end if


 if (mpi_enreg%me == 0 .and. dtset%outputXML == 1) then
  write(ab_xml_out, "(A)") '      <finalConditions>'
! We output the final result given in results_gs
  call out_resultsgs_XML(dtset, 4, results_gs, psps%usepaw)
  write(ab_xml_out, "(A)") '      </finalConditions>'
  write(ab_xml_out, "(A)") '    </scfcvLoop>'
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(60,2,tsec)
 call timab(20,2,tsec)

!DEBUG
!write(6,*)' scfcv : exit '
!stop
!ENDDEBUG

end subroutine scfcv
!!***
