!{\src2tex{textfont=tt}}
!!****f* ABINIT/gstate
!! NAME
!! gstate
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations by CG minimization.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, JYR, MKV, MT, FJ, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial CPU time
!!  nspinor=number of spinorial components of the wavefunctions
!!  walli=initial wall clock time
!!
!! OUTPUT
!!  npwtot(nkpt) = total number of plane waves at each k point
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!
!! SIDE EFFECTS
!!  acell(3)=unit cell length scales (bohr)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband =maximum number of bands (IN)
!!   | mgfft =maximum single fft dimension (IN)
!!   | mkmem =maximum number of k points which can fit in core memory (IN)
!!   | mpw   =maximum number of planewaves in basis sphere (large number) (IN)
!!   | natom =number of atoms in unit cell (IN)
!!   | nfft  =(effective) number of FFT grid points (for this processor) (IN)
!!   | nkpt  =number of k points (IN)
!!   | nspden=number of spin-density components (IN)
!!   | nsppol=number of channels for spin-polarization (1 or 2) (IN)
!!   | nsym  =number of symmetry elements in space group
!!  iexit= exit flag
!!  mpi_enreg=MPI-parallelisation information (some already initialized,
!!   some others to be initialized here)
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in gstate, a significant part of
!!   psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,
!!     ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!   All the remaining components of psps are to be initialized in the call
!!   to pspini .
!!   The next time the code enters gstate, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled
!!   in pspini.f .
!!  rprim(3,3)=dimensionless real space primitive translations
!!  vel(3,natom)=value of velocity
!!  xred(3,natom) = reduced atomic coordinates
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
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
!! In case of wavelets:
!! --------------------
!!    - Only the usual FFT grid (defined by wvl_crmult) is used.
!!      It is defined by nfft, ngfft, mgfft, ... This is strictly not
!!      an FFT grid since its dimensions are not suited for FFTs. They are
!!      defined by wvl_setngfft().
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! TODO
!! Not yet possible to use restartxf in parallel when localrdwf==0
!!
!! PARENTS
!!      driver,pstate
!!
!! CHILDREN
!!      blok8,brdmin,bstruct_clean,bstruct_init,chkexi,clnmpi_fft,clnmpi_gs
!!      clnup1,clnup2,delocint,diisrelax,energies_init,fconv,fixsym,fourdp
!!      getph,handle_ncerr,hdr_clean,hdr_init,hdr_update,indgrid,initberry
!!      initmpi_fft,initmpi_gs,initrhoij,initro,initylmg,int2char4,inwffil
!!      ioarr,ioddb8,jellium,kpgio,leave_new,mkrho,moldyn,move,newocc,outqmc
!!      outwf,outxfhist,pawinit,pawpuxinit,prtene,psddb8,psolver_kernel,pspini
!!      scfcv,setsym,setsymrhoij,setup1,setup2,status,timab,transgrid,wffclose
!!      wffdelete,wffopen,wffreadskiprec,wrtout,wvl_free_type_proj
!!      wvl_free_type_wfs,wvl_init_type_proj,wvl_init_type_wfs,wvl_mkrho
!!      wvl_setboxgeometry,wvl_setngfft,xcomm_world,xme_init,xproc_max
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine gstate(acell,codvsn,cpui,dtfil,dtset,iexit,&
& mpi_enreg,&
& npwtot,nspinor,&
& occ,pawang,pawrad,pawtab,psps,results_gs,rprim,vel,walli,xred)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
#if defined HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_BIGDFT
 use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13iovars
 use interfaces_13paw
 use interfaces_13psp
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_14occeig
 use interfaces_14poisson
 use interfaces_14wvl_wfs
 use interfaces_15common
 use interfaces_16response
 use interfaces_18seqpar
 use interfaces_21drive, except_this_one => gstate
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: iexit,nspinor
 real(dp),intent(in) :: cpui,walli
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(inout) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 type(results_gs_type),intent(inout) :: results_gs
!arrays
 integer,intent(out) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: acell(3),occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: rprim(3,3),vel(3,dtset%natom),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!2   for wavefunction file, new format (version 2.0 and after)    (fform)
!52  for density rho(r)       (fformr)
!102 for potential V(r) file. (fformv)
!scalars
 integer,parameter :: fform=2,fformv=102,formeig=0,level=3,response=0
 integer,save :: nsym_old=-1
 integer :: accessfil,ask_accurate,bantot,blktyp,choice,fformr=52,fullinit
 integer :: gscase,iapp,iatom,icount,idir,ierr,ifft,ii,ilmn,index,initialized
 integer :: ionmov,ios,ir,iscf,ispden,isppol,itime,itimexit,itypat,ixfh,ixx
 integer :: master,me,mgfftf,mpert,mpsang,msize,mu,mxfh,mygroup,nblok,ncerr
 integer :: ncid_hdr,nfftf,nfftftot,nfftot,normchoice,nproc,ntime,ntypat,nxfh
 integer :: nxfhr,openexit,option,optorth,prtvol,psp_gencond,pwind_alloc,rdwr
 integer :: rdwrpaw,restartxf,spaceworld,tim_mkrho,tmkmem,vrsddb
 real(dp) :: cpus,diecut_eff,ecore,ecut_eff,ecutdg_eff,epulay,etot,fermie
 real(dp) :: gsqcut_eff,gsqcutc_eff,residm,rhosum,tolwfr,ucvol
 logical :: ex,od,read_wf_or_den,use_scf_history
 character(len=3) :: ipara
 character(len=4) :: tag
 character(len=500) :: message
 character(len=fnlen) :: ddbnm,dscrpt
 type(bandstructure_type) :: bstruct
 type(dens_sym_operator_type) :: densymop_gs
 type(efield_type) :: dtefield
 type(hdr_type) :: hdr
 type(pawfgr_type) :: pawfgr
 type(scf_history_type) :: scf_history
 type(wffile_type) :: wff1,wffnew,wffnow
 type(wvl_data) :: wvl
!arrays
 integer :: ngfft(18),ngfftf(18)
 integer,save :: paw_gencond(6)=(/-1,-1,-1,-1,-1,-1/)
 integer,allocatable :: atindx(:),atindx1(:),blkflg(:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),nattyp(:),npwarr(:),symrec(:,:,:)
 integer,pointer :: pwind(:,:,:)
 real(dp) :: blknrm(3),blkqpt(9),corstr(6),gmet(3,3),gprimd(3,3),k0(3)
 real(dp) :: rmet(3,3),rprimd(3,3),tsec(2)
 real(dp),allocatable :: amass(:),blkval(:,:),cg(:,:),doccde(:),dyfrx2(:,:,:)
 real(dp),allocatable :: eigen(:),ph1df(:,:),phnons(:,:,:),resid(:),rhowfg(:,:)
 real(dp),allocatable :: rhowfr(:,:),spinat_dum(:,:),start(:,:),work(:)
 real(dp),allocatable :: workg(:,:),xcart(:,:),xfhist(:,:,:,:),xred_old(:,:)
 real(dp),allocatable :: ylm(:,:),ylmgr(:,:,:)
 real(dp),pointer :: kernel_dummy(:),pwnsfac(:,:),rhog(:,:),rhor(:,:)
 character(len=fnlen) :: tmpfil(7)
 type(pawrhoij_type),allocatable :: pawrhoij(:)
!no_abirules

! ***********************************************************************
!DEBUG
!write(6,*)' gstate : enter'
!write(6,*) 'nfft',dtset%nfft
!stop
!ENDDEBUG

 call timab(32,1,tsec)
 call timab(33,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

 if (mpi_enreg%me == 0 .and. dtset%outputXML == 1) then
! gstate() will handle a dataset, so we output the dataSet markup.
  write(ab_xml_out, "(A)") '  <dataSet>'
! We output the variables of the dataset given in argument.
! call outvarsXML()
 end if

!Set up mpi informations from the dataset
 if (dtset%usewvl == 0) then
  if (mpi_enreg%parareel == 0) then
   mpi_enreg%paral_level=2
   call initmpi_gs(dtset,mpi_enreg)
  else
   mpi_enreg%paral_level=1
  end if

  call initmpi_fft(dtset,mpi_enreg)

  nullify(mpi_enreg%nscatterarr, mpi_enreg%ngatherarr)
 else
! This is mandatory since pointers are not null initialised
! and the association is tested in xcomm_init().
  nullify(mpi_enreg%kpt_comm_para)
! WVL - data distribution
  allocate(mpi_enreg%nscatterarr(0:mpi_enreg%nproc - 1, 4))
  allocate(mpi_enreg%ngatherarr(0:mpi_enreg%nproc - 1, 2))
 end if

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 ngfft(:)=dtset%ngfft(:)
 if (psps%usepaw==1) then
  if (dtset%pawecutdg >= 1.0000001_dp*dtset%ecut) then
   pawfgr%usefinegrid=1
   nfftf=dtset%nfftdg;mgfftf=dtset%mgfftdg;ngfftf(:)=dtset%ngfftdg(:)
   nfftot=ngfft(1)*ngfft(2)*ngfft(3)
   nfftftot=ngfftf(1)*ngfftf(2)*ngfftf(3)
   allocate(pawfgr%coatofin(nfftot),pawfgr%fintocoa(nfftftot))
   call indgrid(pawfgr%coatofin,pawfgr%fintocoa,nfftot,nfftftot,ngfft,ngfftf)
  else !this is a simple transfer, this can be done in parallel with only local info
   pawfgr%usefinegrid=0
   nfftf=dtset%nfft;mgfftf=dtset%mgfft;ngfftf(:)=dtset%ngfft(:)
   allocate(pawfgr%coatofin(dtset%nfft),pawfgr%fintocoa(dtset%nfft))
   do ii=1,dtset%nfft
    pawfgr%coatofin(ii)=ii;pawfgr%fintocoa(ii)=ii
   end do
  end if
  pawfgr%natom=dtset%natom
  pawfgr%nfftc=dtset%nfft;pawfgr%mgfftc=dtset%mgfft;pawfgr%ngfftc(:)=dtset%ngfft(:)
  pawfgr%nfft =nfftf     ;pawfgr%mgfft=mgfftf      ;pawfgr%ngfft(:)=ngfftf(:)
  ecutdg_eff = dtset%pawecutdg * (dtset%dilatmx)**2
  ecut_eff   = dtset%ecut      * (dtset%dilatmx)**2
 else
  pawfgr%usefinegrid=0
  nfftf=dtset%nfft;mgfftf=dtset%mgfft;ngfftf(:)=dtset%ngfft(:)
  allocate(pawfgr%coatofin(0),pawfgr%fintocoa(0))
  ecut_eff= dtset%ecut * (dtset%dilatmx)**2
  ecutdg_eff=ecut_eff
 end if

!
!If dtset%accesswff == 2 set all array outputs to netcdf format
!
 accessfil = 0
 if (dtset%accesswff == 1) then
  accessfil = 4
 end if
 if (dtset%accesswff == 2) then
  accessfil = 1
 end if
 if (dtset%accesswff == 3) then
  accessfil = 3
 end if

!Init spaceworld
 call xcomm_world(mpi_enreg,spaceworld)
 master =0
!Define me

!BEGIN TF_CHANGES

 call xme_init(mpi_enreg,me)
!END TF_CHANGES

!Define nproc

 call xproc_max(nproc,ierr)

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&  ' gstate : enter , debug mode '
  call wrtout(06,message,'COLL')
 end if

 ntime=dtset%ntime

!Option input variables
 ionmov   =dtset%ionmov
 iscf     =dtset%iscf
 restartxf=dtset%restartxf

!Create names for the temporary files based on dtfil%filnam_ds(5)
!by appending adequate string.
!'_WF1' -> dtfil%unwft1
!'_WF2' -> dtfil%unwft2
!'_KG' ->  dtfil%unkg
!'_DUM' -> tmp_unit (real dummy name)
!'_YLM' -> dtfil%unylm
!'_PAW' -> dtfil%unpaw
 tmpfil(1)=trim(dtfil%filnam_ds(5))//'_WF1'
 tmpfil(2)=trim(dtfil%filnam_ds(5))//'_WF2'
 tmpfil(3)=trim(dtfil%filnam_ds(5))//'_KG'
 tmpfil(4)=trim(dtfil%filnam_ds(5))//'_DUM'
 tmpfil(6)=trim(dtfil%filnam_ds(5))//'_YLM'
 tmpfil(7)=trim(dtfil%filnam_ds(5))//'_PAW'

 if(mpi_enreg%paral_compil_kpt==1)then
! This is the parallel case : the index of the processor must be appended
  call int2char4(mpi_enreg%me,tag)
  ixx=1
  if (mpi_enreg%paral_compil_mpio == 1 .and. dtset%accesswff == 1 ) ixx=3
  do ii=ixx,7
   tmpfil(ii)=trim(tmpfil(ii))//'_P-'//tag
  end do
 end if

 call status(0,dtfil%filstat,iexit,level,'call setup1   ')

 initialized=0
 ecore=zero

 results_gs%grewtn(:,:)=zero
 call energies_init(results_gs%energies)
 results_gs%pel(1:3)   =zero

!Set up for iterations
 allocate(amass(dtset%natom))
 call setup1(acell,amass,dtset%amu,bantot,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,dtset%iboxcut,dtset%intxc,ionmov,&
& dtset%natom,dtset%nband,ngfftf,ngfft,dtset%nkpt,dtset%nqpt,dtset%nsppol,dtset%nsym,psps%ntypat,&
& dtset%qptn,response,rmet,rprim,rprimd,dtset%typat,ucvol,psps%usepaw)

 allocate(npwarr(dtset%nkpt))
 if (dtset%usewvl == 0) then
  call status(0,dtfil%filstat,iexit,level,'call kpgio    ')
! Set up the basis sphere of planewaves
  allocate(kg(3,dtset%mpw*dtset%mkmem))
  call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,tmpfil(3), &
&  dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,&
&  dtset%mpw,npwarr,npwtot,dtset%nsppol,dtfil%unkg)
 else
  npwarr(:) = 0
  npwtot(:) = 0
 end if

!Set up the Ylm for each k point
 allocate(ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
 allocate(ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm))
 if (psps%useylm==1) then
  if(dtset%mkmem==0) open(dtfil%unylm,file=tmpfil(6),form='unformatted',status='unknown')
  call status(0,dtfil%filstat,iexit,level,'call initylmg ')
  option=0;if (dtset%prtstm==0.and.iscf>0) option=1
  call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&  npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
 end if

 call timab(33,2,tsec)


!Open and read pseudopotential files
 call status(0,dtfil%filstat,iexit,level,'call pspini   ')
 call pspini(dtset,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,level,&
& pawrad,pawtab,psps,rprimd)

 call timab(33,1,tsec)

!In case of isolated computations, ecore must set to zero
!because its contribution is counted in the ewald energy
!as the ion-ion interaction.
 if (dtset%icoulomb == 1) then
  ecore = 0._dp
 end if

!WVL - Now that psp data are available, we compute rprimd, acell...
!from the atomic positions.
 if (dtset%usewvl == 1) then
  call wvl_setBoxGeometry(acell, dtset, mpi_enreg, psps%gth_params%radii_cf, &
&  rprimd, xred)
  rprim(:, :)    = reshape((/ &
&  real(1., dp), real(0., dp), real(0., dp), &
&  real(0., dp), real(1., dp), real(0., dp), &
&  real(0., dp), real(0., dp), real(1., dp) /), (/ 3, 3 /))
  call wvl_setngfft(dtset, mpi_enreg)
  nfftf          = dtset%nfft
  mgfftf         = dtset%mgfft
  ngfftf(:)      = dtset%ngfft(:)
 end if

!Initialize band structure datatype
 allocate(doccde(bantot),eigen(bantot))
 doccde(:)=zero ; eigen(:)=zero
 call bstruct_init(bantot,bstruct,doccde,eigen,dtset%istwfk,dtset%kptns,&
& dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,occ,dtset%wtk)
 deallocate(doccde,eigen)

!DEBUG
!write(6,*)' gstate : return for test memory leak '
!return
!ENDDEBUG

!Initialize PAW atomic occupancies
 if (psps%usepaw==1) then
  allocate(pawrhoij(dtset%natom))
  call initrhoij(psps%indlmn,dtset%lexexch,psps%lmnmax,dtset%lpawu,dtset%natom,dtset%nspden,&
&  dtset%nsppol,dtset%ntypat,pawrhoij,pawtab,dtset%spinat,dtset%typat)
 end if

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
 call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
& residm,rprimd,occ,pawrhoij,psps%usepaw,xred)

!Clean band structure datatype (should use it more in the future !)
 call bstruct_clean(bstruct)

 call status(0,dtfil%filstat,iexit,level,'call inwffil  ')

 allocate(cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol))
 allocate(eigen(dtset%mband*dtset%nkpt*dtset%nsppol))
 allocate(resid(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen(:)=0.0_dp ; resid(:)=0.0_dp
!mpi_enreg%paralbd=0 ; ask_accurate=0
 ask_accurate=0
!WVL - Branching, allocating wavefunctions as wavelets.
 if (dtset%usewvl == 1) then
! Create access arrays for wavefunctions and allocate wvl%wfs%psi (other arrays
! are left unallocated).
  call wvl_init_type_wfs(dtset, mpi_enreg, psps, rprimd, wvl%wfs, xred)
! We transfer wavelets informations to the hdr structure.
#if defined HAVE_BIGDFT
  hdr%nwvlarr(1) = wvl%wfs%keys%nvctr_c
  hdr%nwvlarr(2) = 7 * wvl%wfs%keys%nvctr_f
#endif
! Create access arrays for projectors and allocate them.
! Compute projectors from each atom.
  call wvl_init_type_proj(dtset, mpi_enreg, wvl%projectors, psps, rprimd, xred)
 end if

!XG 020711 : dtfil should not be reinitialized here !!!
 if (mpi_enreg%parareel == 1) then
  if (mpi_enreg%ipara > 0 ) then
   if (mpi_enreg%jpara == 0) then
    dtfil%ireadwf = 0
   else
    dtfil%ireadwf = 0
    if (mpi_enreg%paral_compil_mpio == 1 .and. dtset%accesswff == 1 ) then
     dtfil%fnamewffk=trim(dtfil%filnam_ds(4))//'_WFK'
    else
     if (mpi_enreg%ipara < 11) write(ipara,'(i1)')mpi_enreg%ipara-1
     if (mpi_enreg%ipara >= 11) write(ipara,'(i2)')mpi_enreg%ipara-1
     if (mpi_enreg%ipara >= 101) write(ipara,'(i3)')mpi_enreg%ipara-1
     dtfil%fnamewffk=trim(dtfil%filnam_ds(4))//'_WFK_'//ipara
    end if
   end if
  else
   dtfil%ireadwf = 0
  end if
 end if
 read_wf_or_den=(iscf<=0.or.dtfil%ireadden/=0.or.dtfil%ireadwf/=0)

!Initialize wavefunctions.
!Warning : ideally, results_gs%fermie and results_gs%residm
!should not be initialized here. One might make them separate variables.

 wff1%unwff=dtfil%unwff1
 optorth=1   !if (psps%usepaw==1) optorth=0
 if(psps%usepaw==1 .and. dtfil%ireadwf==1)optorth=0
 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen,dtset%exchn2n3d,&
& formeig,gmet,hdr,dtfil%ireadwf,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,dtset%mband,dtset%mkmem,mpi_enreg,&
& dtset%mpw,dtset%nband,ngfft,dtset%nkpt,npwarr,dtset%nspden,nspinor,dtset%nsppol,dtset%nsym,occ,&
& optorth,psps,prtvol,rprimd,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wff1,wffnow,dtfil%unwff1,dtfil%unwft1,&
& dtfil%fnamewffk,tmpfil(1),wvl)

 if (psps%usepaw==1.and.dtfil%ireadwf==1)then
  do iatom=1,dtset%natom
   pawrhoij(iatom)%nspden=hdr%pawrhoij(iatom)%nspden
   pawrhoij(iatom)%lmn2_size=hdr%pawrhoij(iatom)%lmn2_size
   pawrhoij(iatom)%nrhoijsel=hdr%pawrhoij(iatom)%nrhoijsel+0
   do ilmn=1,pawrhoij(iatom)%nrhoijsel
    pawrhoij(iatom)%rhoijselect(ilmn)=hdr%pawrhoij(iatom)%rhoijselect(ilmn)+0
   end do
   do ispden=1,pawrhoij(iatom)%nspden
    do ilmn=1,pawrhoij(iatom)%nrhoijsel
     pawrhoij(iatom)%rhoijp(ilmn,ispden)=hdr%pawrhoij(iatom)%rhoijp(ilmn,ispden)+zero
    end do
   end do
  end do
! Has to update header again (because pawrhoij has changed)
! MT 2007-10-22: Why ? only values for klmn>nrhoijsel may have changed
! but they are never used...
  call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&  residm,rprimd,occ,pawrhoij,psps%usepaw,xred)
 end if

!DEBUG
!write(6,*)' gstate : stop for test memory leak '
!call hdr_clean(hdr)
!return
!ENDDEBUG

!Initialize xf history (should be put in inwffil)
 nxfh=0
 if(restartxf>=1 .and. dtfil%ireadwf==1)then

! Should exchange the data about history in parallel localrdwf==0
  if(mpi_enreg%paral_compil_kpt==1 .and. dtset%localrdwf==0)then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' gstate : BUG -',ch10,&
&   '  It is not yet possible to use non-zero restartxf,',ch10,&
&   '  in parallel, when localrdwf=0. Sorry for this ...'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

  allocate(xfhist(3,dtset%natom+4,2,0))
  call outxfhist(nxfh,dtset%natom,mxfh,xfhist,2,wff1,ios)
  deallocate(xfhist)

  if(ios>0)then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' gstate : BUG -',ch10,&
&   '  An error occurred reading the input wavefunction file,',ch10,&
&   '  with restartxf=1.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  else if(ios==0)then
   write(message, '(a,a,i4,a)' )ch10,&
&   ' gstate : reading',nxfh,' (x,f) history pairs from input wf file.'
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')
  end if
! WARNING : should check that restartxf is not negative
! WARNING : should check that restartxf /= only when dtfil%ireadwf is activated
 end if

!Allocate the xf history array : takes into account the existing
!pairs, minus those that will be discarded, then those that will
!be computed, governed by ntime, and some additional pairs
!(needed when it will be possible to use xfhist for move.f)
 mxfh=(nxfh-restartxf+1)+ntime+5
 if(mpi_enreg%parareel==1)mxfh=mxfh+500  ! XG020711 : why this value ?
 allocate(xfhist(3,dtset%natom+4,2,mxfh))
!WARNING : should check that the number of atoms in the wf file and natom
!are the same

!Initialize the xf history array
 if(nxfh>=restartxf .and. nxfh>0)then
! Eventually skip some of the previous history
  if(restartxf>=2)then
   do ixfh=1,restartxf-1
    call WffReadSkipRec(ios,1,wff1)
   end do
  end if

! Read and store the relevant history
  nxfhr=nxfh-restartxf+1
  call outxfhist(nxfhr,dtset%natom,mxfh,xfhist,3,wff1,ios)
 end if

!Determine whether SCF history has to be used
 use_scf_history=(ionmov>0.and.dtset%usewvl==0.and. &
& (abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6))

!Close wff1, if it was ever opened (in inwffil)
 if (dtfil%ireadwf==1) then
  call WffClose(wff1,ierr)
 end if

!Initialize second wavefunction file if needed
 if(dtset%mkmem==0 .and. dtset%nstep/=0) then
  write(message, '(a,i4,a,a)' )&
&  ' gstate about to open unit',dtfil%unwft2,' for file=',trim(tmpfil(2))
  call wrtout(06,message,'PERS')

#if defined HAVE_NETCDF
  if(dtset%accesswff==2) then
!  Create empty netCDF file
   ncerr = nf90_create(path=trim(tmpfil(2)), cmode=NF90_CLOBBER, ncid=ncid_hdr)
   call handle_ncerr(ncerr," create netcdf wavefunction file")
   ncerr = nf90_close(ncid_hdr)
   call handle_ncerr(ncerr," close netcdf wavefunction file")
  else if(dtset%accesswff==3) then
   write (std_out,*) "FIXME: ETSF I/O support in gstate"
  end if
#endif

  call WffOpen(dtset%accesswff,spaceworld,tmpfil(2),ierr,wffnew,master,me,dtfil%unwft2)
 end if

 call status(0,dtfil%filstat,iexit,level,'call setup2   ')

!Further setup
 allocate(start(3,dtset%natom))
 call setup2(dtset, results_gs%energies%e_pulay,iscf, &
& npwtot,start,ucvol,wvl%wfs,xred)

!Allocation for forces and atomic positions
 allocate(xred_old(3,dtset%natom))
 xred_old = xred

!Do symmetry stuff only for nsym>1
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 allocate(irrzon(nfftot**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol))
 allocate(phnons(2,nfftot**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol))
 irrzon(:,:,:)=0
 allocate(indsym(4,dtset%nsym,dtset%natom),symrec(3,3,dtset%nsym))

 if (dtset%nsym>1) then

  call status(0,dtfil%filstat,iexit,level,'call setsym   ')
  call setsym(densymop_gs,indsym,irrzon,iscf,dtset%natom,&
&  nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
&  phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

! Make sure dtset%iatfix does not break symmetry
  call status(0,dtfil%filstat,iexit,level,'call fixsym   ')
  call fixsym(dtset%iatfix,indsym,dtset%natom,dtset%nsym)

 else

! The symrec array is used by initberry even in case nsym = 1
  symrec(:,:,1) = 0
  symrec(1,1,1) = 1 ; symrec(2,2,1) = 1 ; symrec(3,3,1) = 1

 end if

!Electric field: initialization
 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4)) then
  nullify(pwind,pwnsfac)
  call initberry(dtefield,dtfil,dtset,gmet,kg,dtset%mband,dtset%mkmem,mpi_enreg,&
&  dtset%mpw,dtset%nkpt,npwarr,dtset%nsppol,dtset%nsym,occ,pwind,pwind_alloc,pwnsfac,rprimd,symrec)
 else
  pwind_alloc = 1
  allocate(pwind(pwind_alloc,2,3),pwnsfac(2,pwind_alloc))
 end if

!Timing for initialisation period
 call timab(33,2,tsec)
 call timab(34,1,tsec)

!Compute new occupation numbers, in case wavefunctions and eigenenergies
!were read from disk, occupation scheme is metallic (this excludes iscf=-1),
!and occupation numbers are required by iscf
 if( dtfil%ireadwf==1 .and. &
& (dtset%occopt>=3.and.dtset%occopt<=7) .and. &
& (iscf>0 .or. iscf==-3) ) then

  call status(0,dtfil%filstat,iexit,level,'call newocc   ')
  allocate(doccde(dtset%mband*dtset%nkpt*dtset%nsppol))
! Warning : ideally, results_gs%entropy should not be set up here XG 20011007
! Warning : ideally, results_gs%fermie should not be set up here XG 20011007
! Do not take into account the possible STM bias
  call newocc(doccde,eigen,results_gs%energies%entropy,&
&  results_gs%energies%e_fermie,&
&  dtset%fixmom,dtset%mband,dtset%nband,&
&  dtset%nelect,dtset%nkpt,nspinor,dtset%nsppol,occ,&
&  dtset%occopt,prtvol,zero,dtset%tphysel,dtset%tsmear,dtset%wtk)
  deallocate(doccde)

 else
! Warning : ideally, results_gs%entropy should not be set up here XG 20011007
  results_gs%energies%entropy=zero
 end if

!Generate an index table of atoms, in order for them to be used
!type after type.
 ntypat=psps%ntypat
 allocate(atindx(dtset%natom),atindx1(dtset%natom),nattyp(ntypat))
 index=1
 do itypat=1,ntypat
  nattyp(itypat)=0
  do iatom=1,dtset%natom
   if(dtset%typat(iatom)==itypat)then
    atindx(iatom)=index
    atindx1(index)=iatom
    index=index+1
    nattyp(itypat)=nattyp(itypat)+1
   end if
  end do
 end do

!Compute structure factor phases for current atomic pos:
 if ((.not.read_wf_or_den).or.use_scf_history) then
  allocate(ph1df(2,3*(2*mgfftf+1)*dtset%natom))
  call status(0,dtfil%filstat,iexit,level,'call getph    ')
  call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 end if

!PAW: 1- Initialize values for several arrays unchanged during iterations
!2- Initialize data for LDA+U
!3- Eventually open temporary storage file
 if(psps%usepaw==1) then
! 1-
  if (psp_gencond==1.or.&
&     paw_gencond(1)/=dtset%pawlcutd .or.paw_gencond(2)/=dtset%pawlmix  .or.&
&     paw_gencond(3)/=dtset%pawnphi  .or.paw_gencond(4)/=dtset%pawntheta.or.&
&     paw_gencond(5)/=dtset%pawspnorb.or.paw_gencond(6)/=dtset%pawxcdev) then
   call timab(553,1,tsec)
   diecut_eff=abs(dtset%diecut)*dtset%dilatmx**2
   call pawinit(diecut_eff,psps%indlmn,dtset%pawlcutd,dtset%pawlmix,psps%lmnmax,psps%mpsang,psps%n1xccc,&
&   dtset%pawnphi,dtset%nsym,dtset%pawntheta,psps%ntypat,&
&   pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev)
   paw_gencond(1)=dtset%pawlcutd ; paw_gencond(2)=dtset%pawlmix
   paw_gencond(3)=dtset%pawnphi  ; paw_gencond(4)=dtset%pawntheta
   paw_gencond(5)=dtset%pawspnorb; paw_gencond(6)=dtset%pawxcdev
   call timab(553,2,tsec)
  end if
  if (psp_gencond==1.or.nsym_old/=dtset%nsym) then
   call setsymrhoij(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,&
&   rprimd,dtset%symafm,symrec,pawang%zarot)
   nsym_old=dtset%nsym
  end if
! 2-Initialize and compute data for LDA+U
  if (dtset%usepawu>0.or.dtset%useexexch>0) then
   call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%jpawu,dtset%lexexch,dtset%lpawu,&
&   psps%indlmn,psps%lmnmax,ntypat,pawang,dtset%pawprtvol,pawrad,pawtab,dtset%upawu,&
&   dtset%useexexch,dtset%usepawu)
  end if
! 3-Eventually open temporary storage file
  if(dtset%mkmem==0) then
   open(dtfil%unpaw,file=tmpfil(7),form='unformatted',status='unknown')
   rewind(unit=dtfil%unpaw)
  end if
 end if

!Get starting charge density : rhor as well as rhog
 allocate(rhog(2,nfftf),rhor(nfftf,dtset%nspden))
 if (iscf>0) then
  if(dtfil%ireadden/=0)then

   rdwr=1;rdwrpaw=psps%usepaw;if(dtfil%ireadwf/=0) rdwrpaw=0
   call ioarr(accessfil,rhor,dtset,results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&   mpi_enreg, nfftf,pawrhoij,rdwr,rdwrpaw,ngfft)
   if (rdwrpaw/=0) then
    call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&    residm,rprimd,occ,pawrhoij,psps%usepaw,xred)
   end if
!  Compute up+down rho(G) by fft
   allocate(work(nfftf));work(:)=rhor(:,1)
   call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
   deallocate(work)

  else if(dtfil%ireadwf/=0)then

!  Obtain the charge density from wfs that were read previously
!  Be careful: in PAW, rho does not include the compensation
!  density (to be added in scfcv.F90) !
   call status(0,dtfil%filstat,iexit,level,'call mkrho    ')
!  tim_mkrho=1 ; mpi_enreg%paralbd=0
   tim_mkrho=1
   if (psps%usepaw==1) then
    allocate(rhowfg(2,dtset%nfft),rhowfr(dtset%nfft,dtset%nspden))
    call mkrho(cg,densymop_gs,dtset,irrzon,kg,&
&    mpi_enreg,npwarr,nspinor,occ,phnons,rhowfg,rhowfr,tim_mkrho,ucvol,&
&    dtfil%unkg,wffnow,wvl%wfs)
    call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
    deallocate(rhowfg,rhowfr)
   else
    call mkrho(cg,densymop_gs,dtset,irrzon,kg,&
&    mpi_enreg,npwarr,nspinor,occ,phnons,rhog,rhor,tim_mkrho,ucvol,&
&    dtfil%unkg,wffnow,wvl%wfs)

!   DEBUG
!   write(6,*)' gstate : after mkrho '
!   ENDDEBUG
   end if

  else if(dtfil%ireadwf==0)then

!  Crude, but realistic initialisation of the density
!  There is not point to compute it from random wavefunctions
!  except with wavelets.
   call status(0,dtfil%filstat,iexit,level,'call initro   ')
   if (dtset%usewvl == 0) then
    call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&    dtset%natom,nattyp,nfftf,ngfftf,dtset%nspden,ntypat,dtset%paral_kgb,&
&    pawtab,ph1df,psps%qgrid_vl,rhog,rhor,&
&    dtset%spinat,ucvol,psps%usepaw,dtset%ziontypat,dtset%znucl)
!   Update initialized density taking into account jellium slab
    if(dtset%jellslab/=0) then
     option=2; allocate(work(nfftf))
     call jellium(gmet,gsqcut_eff,mpi_enreg,nfftf,ngfftf,dtset%nspden,&
&     option,dtset%paral_kgb,dtset%slabwsrad,rhog,rhor,rprimd,work,dtset%slabzbeg,dtset%slabzend)
     deallocate(work)
    end if ! of usejell
   else
    call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wvl%wfs)
   end if

  end if

 else if (iscf==-1.or.iscf==-2.or.iscf==-3) then

  call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
! Read rho(r) from a disk file
  rdwr=1;rdwrpaw=psps%usepaw
! Note : results_gs%etotal is read here,
! and might serve in the tddft routine, but it is contrary to the
! intended use of results_gs ...
! Warning : should check the use of results_gs%fermie
! Warning : should check the use of results_gs%residm
! One might make them separate variables.

! DEBUG
! write(6,*)' gstate : before ioarr, reading the density '
! ENDDEBUG

  call ioarr(accessfil,rhor,dtset, results_gs%etotal,fformr,dtfil%fildensin,hdr,&
&  mpi_enreg,nfftf,pawrhoij,rdwr,rdwrpaw,ngfft)

! Compute up+down rho(G) by fft
  call status(0,dtfil%filstat,iexit,level,'call fourdp   ')
  allocate(work(nfftf))
  work(:)=rhor(:,1)
  call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
  deallocate(work)

 else

! Disallowed value for iscf
  write(message, '(a,a,a,a,i12,a)' )  ch10,&
&  ' gstate : BUG -',ch10,&
&  '  iscf has disallowed value=',iscf,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')

 end if

!Debugging : print the different parts of rhor
!MPIWF Warning : this should not be parallelized over space, leave this debugging feature as such.
 if(prtvol==-level)then
  write(message,'(a)') '   ir     rhor(ir)     '
  call wrtout(06,message,'COLL')
  do ir=1,nfftf
   if(ir<=11 .or. mod(ir,301)==0 )then
    write(message,'(i5,a,es13.6)')ir,' ',rhor(ir,1)
    call wrtout(06,message,'COLL')
    if(dtset%nsppol==2)then
     write(message,'(a,es13.6)')'      ',rhor(ir,2)
     call wrtout(06,message,'COLL')
    end if
   end if
  end do
 end if

!If needed, allocate and initialize SCF history variables
 if (use_scf_history) then
! Allocations
  scf_history%history_size=2
  scf_history%natom=dtset%natom
  scf_history%nfft=nfftf
  scf_history%nspden=dtset%nspden
  allocate(scf_history%hindex(scf_history%history_size));scf_history%hindex(:)=0
  allocate(scf_history%deltarhor(nfftf,dtset%nspden,scf_history%history_size))
  allocate(scf_history%xreddiff(3,dtset%natom,scf_history%history_size))
  allocate(scf_history%atmrho_last(nfftf))
  if (psps%usepaw==1) allocate(scf_history%pawrhoij(dtset%natom,scf_history%history_size))
! If rhor is an atomic density, just store it in history
  if (.not.read_wf_or_den) then
   scf_history%atmrho_last(:)=rhor(:,1)
  else
!  If rhor is not an atomic density, has to compute rho_at(r)
   allocate(rhowfg(2,nfftf),rhowfr(nfftf,1))
   allocate(spinat_dum(3,dtset%natom));spinat_dum=zero
   call initro(atindx,dtset%densty,gmet,gsqcut_eff,psps%usepaw,mgfftf,mpi_enreg,&
&   psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,1,ntypat,dtset%paral_kgb,pawtab,&
&   ph1df,psps%qgrid_vl,rhowfg,rhowfr,spinat_dum,ucvol,&
&   psps%usepaw,dtset%ziontypat,dtset%znucl)
   scf_history%atmrho_last(:)=rhowfr(:,1)
   deallocate(rhowfg,rhowfr,spinat_dum)
  end if
 else
  scf_history%history_size=0
 end if

 if ((.not.read_wf_or_den).or.use_scf_history) deallocate(ph1df)

 call status(0,dtfil%filstat,iexit,level,'end gstate(1) ')

 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i1,a)') ch10,&
&  ' gstate : before scfcv, move or brdmin ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call timab(34,2,tsec)
!Check whether exiting was required by the user.
!If found then do not start minimization steps
!At this first call to chkexi, initialize cpus, if it
!is non-zero (which would mean that no action has to be taken)
!Should do this in driver ...
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1 ; if(dtset%chkexit==0) openexit=0
 call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
!If immediate exit, and wavefunctions were not read, must zero eigenvalues
 if (iexit/=0) then
  eigen(:)=zero
 end if
 if (iexit==0) then

! #######################################################################

! If atoms are not being moved, use scfcv directly; else
! call move or brdmin which in turn calls scfcv.

  call timab(35,1,tsec)

  write(message,'(a,80a)')ch10,('=',mu=1,80)
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,message,'COLL')
  if (ionmov==0) then

   call status(0,dtfil%filstat,iexit,level,'call scfcv    ')

!  Should merge this call with the call for ionmov==4 and 5
   iapp=0
!  mpi_enreg%paralbd=0
   call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,&
&   ecore,eigen,hdr,iapp,indsym,initialized,&
&   irrzon,kg,mpi_enreg,nattyp,nfftf,npwarr,nspinor,occ,&
&   pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&   scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

  else if (ionmov==1) then
!  Conduct molecular dynamics, with or without viscous damping

   call status(0,dtfil%filstat,iexit,level,'call move     ')
!  mpi_enreg%paralbd=0
   call move(acell,amass,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,&
&   ecore,eigen,hdr,indsym,initialized,irrzon,&
&   kg,mpi_enreg,&
&   nattyp,nfftf,npwarr,nspinor,occ,&
&   pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&   scf_history,symrec,wffnew,wffnow,vel,wvl,xred,xred_old,ylm,ylmgr)

  else if (ionmov==2 .or. ionmov==3) then

!  Apply Broyden method for structural optimization, as
!  implemented by Jean-Christophe Charlier (May 1992)

   call status(0,dtfil%filstat,iexit,level,'call brdmin   ')
!  mpi_enreg%paralbd=0

   call brdmin(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,&
&   ecore,eigen,hdr,indsym,initialized,irrzon,&
&   kg,mpi_enreg,mxfh,&
&   nattyp,nfftf,npwarr,nspinor,nxfh,occ,&
&   pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprim,&
&   scf_history,symrec,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)
!  call mkrdim(acell,rprim,rprimd)

  else if (ionmov==4 .or. ionmov==5) then

   do itime=1,ntime

    call status(itime,dtfil%filstat,iexit,level,'call scfcv(mv)')

    if(ionmov==4)then
     if(mod(itime,2)==1)then
      write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&      ' : OPTIMIZE ELECTRONS ------------------------------------'
     else
      write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&      ' : OPTIMIZE ELECTRONS AND IONS ---------------------------'
     end if
    else
     write(message, '(a,a,i3,a)' ) ch10,' STEP NUMBER ',itime,&
&     ' : SIMPLE RELAXATION -------------------------------------'
    end if
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')

!   In this case, iapp is simply itime
    iapp=itime
!   mpi_enreg%paralbd=0
    call scfcv(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,dtset,ecore,&
&    eigen,hdr,iapp,indsym,initialized,irrzon,kg,mpi_enreg,&
&    nattyp,nfftf,npwarr,nspinor,occ,pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&    phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprimd,&
&    scf_history,symrec,wffnew,wffnow,wvl,xred,xred_old,ylm,ylmgr)

    if(mod(itime,2)==1)then
!    When the SCF cycle dealt with electrons only,
!    check whether forces are below tolerance; if so, exit
!    from the itime loop
     itimexit=0 ; if(itime==ntime)itimexit=1
     call fconv(results_gs%fcart,dtset%iatfix,itimexit,itime,dtset%natom,&
&     ntime,0,1.0_dp,dtset%strtarget,results_gs%strten,dtset%tolmxf)
    end if
    if (itimexit/=0) exit

!   Check whether exiting was required by the user.
!   If found then beat a hasty exit from time steps
    if(dtset%chkexit==0) then
     openexit=0
    else
     openexit=1
    end if
    call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
    if (iexit/=0) then
     iexit=0   ! In order not to exit of dataset loop automatically
     exit
    end if

   end do

  else if ( (ionmov>=6 .and. ionmov<=9) .or. (ionmov>=12 .and. ionmov<=13) ) then

!  Molecular dynamics, using Verlet algorithm (ionmov=6)
!  or fake molecular dynamics for minimisation (ionmov=7)
!  or true molecular dynamics with Nose thermostat (ionmov=8)
!  or Langevin dynamics (ionmov=9) or Fei Zhang algorithm (ionmov=12)

   call status(0,dtfil%filstat,iexit,level,'call moldyn   ')

   call moldyn(acell,amass,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&   dtset,ecore,eigen,hdr,indsym,initialized,&
&   irrzon,kg,mpi_enreg,mxfh,&
&   nattyp,nfftf,npwarr,nspinor,nxfh,occ,&
&   pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprim,&
&   scf_history,symrec,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)
!  call mkrdim(acell,rprim,rprimd)

  else if (ionmov == 10) then

   call delocint(acell,atindx,atindx1,cg,cpus,densymop_gs,dtefield,dtfil,&
&   dtset,ecore,eigen,hdr,indsym,initialized,irrzon,&
&   kg,mpi_enreg,mxfh,&
&   nattyp,nfftf,npwarr,nspinor,nxfh,occ,&
&   pawang,pawfgr,pawrad,pawrhoij,pawtab,&
&   phnons,psps,pwind,pwind_alloc,pwnsfac,resid,results_gs,rhog,rhor,rprim,&
&   scf_history,symrec,wffnew,wffnow,vel,wvl,xfhist,xred,xred_old,ylm,ylmgr)
!  call mkrdim(acell,rprim,rprimd)

  else if (ionmov == 20) then
!  Ground state call.
   iapp = 0
!  Ionic positions relaxation using DIIS. This algorithm is fast
!  and converge to the nearest singular point (where gradient vanishes).
!  This is a good algorithm to precisely tune saddle-points.
   call diisRelax(acell, atindx, atindx1, cg, cpus, densymop_gs, dtefield, &
&   dtfil, dtset, ecore, eigen, hdr, iapp, indsym, initialized, &
&   irrzon, kg, mpi_enreg, nattyp, nfftf, npwarr, nspinor, occ, pawang, &
&   pawfgr, pawrad, pawrhoij, pawtab, phnons, psps, pwind, pwind_alloc, pwnsfac, &
&   resid, results_gs, rhog, rhor, rprimd, scf_history, symrec, &
&   wffnew, wffnow, wvl, xred, xred_old, ylm, ylmgr)

  else
!  Not an allowed option
   write(message, '(a,a,a,a,i12,a,a)' ) ch10,&
&   ' gstate : BUG -',ch10,&
&   '  Disallowed value for ionmov=',ionmov,ch10,&
&   '  Allowed values are 0 to 5.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  call timab(35,2,tsec)

! #####################################################################

! End of the check of hasty exit
 end if

 call timab(36,1,tsec)

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
& ' ----iterations are completed or convergence reached----',&
& ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

!Close the unneeded temporary data files, if any.
!Other files are closed in clnup1.
 if (dtset%mkmem==0) then
  close (unit=dtfil%unkg,status='delete')
  if (psps%useylm==1) close (unit=dtfil%unylm,status='delete')
  if (psps%usepaw==1) close (unit=dtfil%unpaw,status='delete')
  call WffDelete(wffnew,ierr)
 end if

!Will be put here later.
!!$ ! WVL - maybe compute the tail corrections to energy
!!$ if (dtset%tl_radius > real(0, dp)) then
!!$    ! Store xcart for each atom
!!$    allocate(xcart(3, dtset%natom))
!!$    call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
!!$    ! Use the tails to improve energy precision.
!!$    call wvl_tail_corrections(dtset, results_gs%energies, results_gs%etotal, &
!!$         & mpi_enreg, occ, psps, vtrial, wvl, xcart)
!!$    deallocate(xcart)
!!$ end if

!Update the header, before using it
 call hdr_update(bantot,results_gs%etotal,results_gs%energies%e_fermie,hdr,dtset%natom,&
& results_gs%residm,rprimd,occ,pawrhoij,psps%usepaw,xred)

 call status(0,dtfil%filstat,iexit,level,'call outwf    ')

 call outwf(cg,dtfil,dtset,eigen,dtfil%filnam_ds(4),hdr,kg,dtset%kptns,&
& dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,mxfh,dtset%natom,dtset%nband,dtset%nfft,ngfft,&
& dtset%nkpt,npwarr,dtset%nqpt,nspinor,dtset%nsppol,dtset%nstep,nxfh,&
& occ,resid,response,wffnow,wvl%wfs,xfhist)

 if(dtset%prtwf==2)then
  call outqmc(cg,dtset,eigen,gprimd,hdr,kg,dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nkpt,npwarr,&
&  nspinor,dtset%nsppol,occ,psps,results_gs)
 end if

 call status(0,dtfil%filstat,iexit,level,'call clnup1   ')

 call clnup1(acell,dtset%dosdeltae,dtset,eigen,dtset%enunit,&
& results_gs%energies%e_fermie,dtfil%filnam_ds(4),&
& results_gs%fred,dtset%iatfix,iscf,dtset%kptns,dtset%kptopt,dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,&
& dtset%natom,dtset%nband,nfftf,ngfftf,dtset%nkpt,dtset%nspden,nspinor,dtset%nsppol,dtset%nstep,&
& occ,dtset%occopt,dtset%prtdos,dtset%prteig,dtset%optforces,dtset%prtstm,prtvol,&
& resid,rhor,rprimd,dtset%tphysel,dtset%tsmear,results_gs%vxcavg,dtset%wtk,xred)

 if ( (iscf>0 .or. iscf==-3) .and. dtset%prtstm==0) then
  call status(0,dtfil%filstat,iexit,level,'call prtene   ')
  call prtene(dtset,results_gs%energies,ab_out,psps%usepaw)
 end if

!Open the formatted derivative database file, and write the
!preliminary information
!In the // case, only one processor writes the energy and
!the gradients to the DDB

 if ((psps%usepaw == 0).and.(mpi_enreg%me==0).and.((iscf > 0).or.&
& (dtset%berryopt == -1).or.(dtset%berryopt) == -3)) then

  call status(0,dtfil%filstat,iexit,level,'call ioddb8   ')
  vrsddb=010929
  dscrpt=' Note : temporary (transfer) database '
  choice=2
  ddbnm=trim(dtfil%filnam_ds(4))//'_DDB'
! tolwfr must be initialized here, but it is a dummy value
  tolwfr=1.0_dp
  call ioddb8 (choice,dscrpt,ddbnm,dtset%natom,dtset%mband,&
&  dtset%nkpt,dtset%nsym,psps%ntypat,dtfil%unddb,vrsddb,&
&  acell,dtset%amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&  dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&  dtset%natom,dtset%nband,ngfft,dtset%nkpt,dtset%nspden,nspinor,&
&  dtset%nsppol,dtset%nsym,psps%ntypat,occ,dtset%occopt,&
&  rprim,dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&  dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&  dtset%typat,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

  if (iscf > 0) then
   nblok = 2          ! 1st blok = energy, 2nd blok = gradients
  else
   nblok = 1
  end if
  fullinit = 0
  call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&  psps%lmnmax,psps%lnmax,nblok,&
&  psps%ntypat,dtfil%unddb,psps%pspso,psps%usepaw,psps%useylm,vrsddb)

  mpert = dtset%natom + 6   ; msize = 3*mpert
  allocate(blkflg(msize),blkval(2,msize))

  blkflg(:) = 0       ; blkval(:,:) = zero
  blkqpt(:) = zero    ; blknrm(:) = one

! Write total energy to the DDB
  if (iscf > 0) then
   blktyp = 0
   blkval(1,1) = results_gs%etotal
   blkflg(1) = 1
   call blok8(blkflg,blknrm,blkqpt,blktyp,blkval,choice,dtset%mband,&
&   mpert,msize,dtset%natom,dtset%nkpt,dtfil%unddb)
  end if

! Write gradients to the DDB
  blktyp = 4
  blkflg(:) = 0       ; blkval(:,:) = zero
  index = 0
  if (iscf > 0) then
   do iatom = 1, dtset%natom
    do idir = 1, 3
     index = index + 1
     blkflg(index) = 1
     blkval(1,index) = results_gs%fred(idir,iatom)
    end do
   end do
  end if

  index = 3*dtset%natom + 3
  if ((abs(dtset%berryopt) == 1).or.(abs(dtset%berryopt) == 3)) then
   do idir = 1, 3
    index = index + 1
    if (dtset%rfdir(idir) == 1) then
     blkflg(index) = 1
     blkval(1,index) = results_gs%pel(idir)
    end if
   end do
  end if

  index = 3*dtset%natom + 6
  if (iscf > 0) then
   blkflg(index+1:index+6) = 1
   blkval(1,index+1:index+6) = results_gs%strten(1:6)
  end if

  call blok8(blkflg,blknrm,blkqpt,blktyp,blkval,choice,dtset%mband,&
&  mpert,msize,dtset%natom,dtset%nkpt,dtfil%unddb)

  deallocate(blkflg,blkval)

! Close DDB
  close(dtfil%unddb)

 end if

 if (dtset%nstep>0 .and. dtset%prtstm==0 .and. dtset%positron==0) then
  call status(0,dtfil%filstat,iexit,level,'call clnup2   ')
  call clnup2(psps%n1xccc,results_gs%fred,results_gs%gresid,&
&  results_gs%grewtn,&
&  results_gs%grxc,iscf,dtset%natom,dtset%optforces,dtset%optstress,prtvol,start,&
&  results_gs%strten,results_gs%synlgr,psps%usepaw,xred)
 end if

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')

!Deallocate arrays
 deallocate(amass,atindx,atindx1,cg,eigen,indsym)
 deallocate(irrzon,npwarr,nattyp,phnons,resid)
 deallocate(rhog,rhor,start,symrec,xfhist,xred_old,ylm,ylmgr)
 deallocate(pawfgr%fintocoa,pawfgr%coatofin)

 if (psps%usepaw==1) then
  do iatom=1,dtset%natom
   deallocate(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect)
  end do
  deallocate(pawrhoij)
 end if

 if (use_scf_history) then
  if (psps%usepaw==1) then
   do ii=1,scf_history%history_size
    if (scf_history%hindex(ii)>0) then
     ixx=scf_history%hindex(ii)
     do iatom=1,dtset%natom
      deallocate(scf_history%pawrhoij(iatom,ixx)%rhoijselect,&
&      scf_history%pawrhoij(iatom,ixx)%rhoijp)
     end do
    end if
   end do
   deallocate(scf_history%pawrhoij)
  end if
  deallocate(scf_history%hindex,scf_history%deltarhor,scf_history%xreddiff,&
&  scf_history%atmrho_last)
 end if
!RShaltaf: Changed to include SBC
 if (dtset%icoulomb > 0) then
! Ask to deallocate the kernel part of Poisson's solver
! Arguments are dummy ones since iaction == 0.
  call psolver_kernel(dtset, 0, kernel_dummy, mpi_enreg, rprimd)
 end if

 if (dtset%usepawu>0.or.dtset%useexexch>0) then
  do itypat=1,ntypat
   if((dtset%lpawu(itypat)/=-1).or.(dtset%lexexch(itypat)/=-1)) deallocate(pawtab(itypat)%lnproju)
   if((dtset%lpawu(itypat)/=-1).or.(dtset%lexexch(itypat)/=-1)) deallocate(pawtab(itypat)%phiphjint)
  end do
 end if

!Deallocating the basis set.
 if (dtset%usewvl == 1) then
  call wvl_free_type_wfs(wvl%wfs)
  call wvl_free_type_proj(wvl%projectors)
 else
  deallocate(kg)
 end if

 if ((dtset%berryopt < 0).or.(dtset%berryopt == 4)) then
  deallocate(pwind,pwnsfac)
  deallocate(dtefield%ikpt_dk,dtefield%idxkstr)
  deallocate(dtefield%sflag,dtefield%cgindex,dtefield%kgindex)
  deallocate(dtefield%fkptns,dtefield%indkk_f2ibz,dtefield%i2fbz)
  if (mpi_enreg%paral_compil_kpt == 1) then
   deallocate(mpi_enreg%kptdstrb)
   if (dtset%berryopt == 4) then
    deallocate(mpi_enreg%kptdstrbi,dtefield%cgqindex,dtefield%nneigh)
   end if
  end if
 else
  deallocate(pwind,pwnsfac)
 end if
 if (dtset%berryopt == 4) deallocate(dtefield%smat)

!Clean the header
 call hdr_clean(hdr)

!DEBUG
!write(6,*)' gstate : return for test memory leak '
!return
!ENDDEBUG

 if (mpi_enreg%me == 0 .and. dtset%outputXML == 1) then
! The dataset given in argument has been treated, then
! we output its variables.
! call outvarsXML()
! gstate() will handle a dataset, so we output the dataSet markup.
  write(ab_xml_out, "(A)") '  </dataSet>'
 end if

 if (dtset%usewvl == 0) then
! Clean the MPI informations
  if (mpi_enreg%parareel == 0) then
   call clnmpi_gs(dtset,mpi_enreg)
  end if

  call clnmpi_fft(dtset,mpi_enreg)
 else
! Clean the wavelet part.
  if (associated(mpi_enreg%nscatterarr)) then
   deallocate(mpi_enreg%nscatterarr)
  end if
  if (associated(mpi_enreg%ngatherarr)) then
   deallocate(mpi_enreg%ngatherarr)
  end if
 end if

 write(message, '(a,a)' ) ch10,' gstate : exiting '
 call wrtout(06,message,'COLL')

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(36,2,tsec)
 call timab(32,2,tsec)

end subroutine gstate
!!***
