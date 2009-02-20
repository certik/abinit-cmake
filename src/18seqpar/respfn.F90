!{\src2tex{textfont=tt}}
!!****f* ABINIT/respfn
!! NAME
!! respfn
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of Response functions.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial cpu time
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum single fft dimension
!!   | mkmem=maximum number of k points which can fit in core memory
!!   | mpw=maximum number of planewaves in basis sphere (large number)
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=number of channels for spin-polarization (1 or 2)
!!   | nsym=number of symmetry elements in space group
!!  mkmems(3)=array containing the tree values of mkmem (see above) (k-GS, k+q-GS and RF)
!!  mpi_enreg=informations about MPI parallelization
!!  npwtot(nkpt)=number of planewaves in basis and boundary at each k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  walli=initial wall clock time
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  etotal=total energy (sum of 7 or 8 contributions) (hartree)
!!
!! SIDE EFFECTS
!!  iexit=index of "exit" on first line of file (0 if not found)
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!    Occupations number may have been read from a previous dataset...
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!    Some dimensions in pawrad have been set in driver.f
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!    Some dimensions in pawtab have been set in driver.f
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    Before entering the first time in respfn, a significant part of psps
!!    has been initialized: the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,
!!    mpsso,mgrid,ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!    All the remaining components of psps are to be initialized in the call
!!    to pspini.  The next time the code enters respfn, psps might be identical
!!    to the one of the previous dtset, in which case, no reinitialisation
!!    is scheduled in pspini.f .
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
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      atm2fft,bstruct_clean,bstruct_init,chkexi,chkph3,d2sym3,distrb2,dyfnl3,dyfro3
!!      dyout3,dyxc13,eltfrhar3,eltfrkin3,eltfrloc3,eltfrnl3,eltfrxc3,ewald3
!!      ewald4,fourdp,gath3,getcut,getph,indgrid,hdr_clean,hdr_init,hdr_update
!!      initmpi_fft,initylmg,int2char4,inwffil,ioarr,ioddb8,kpgio,leave_new
!!      loper3,mkcore,mklocl,mkrho,mpi_barrier,mpi_bcast,newocc,pawinit,pawmknhat,phfrq3,prtph3
!!      psddb8,pspini,q0dy3,rhohxc,rhoij_free,rhoij_copy,setsym,setsymrhoij,setup1,status,symkchk,symq3,symzat
!!      syper3,timab,transgrid,wffclose,wffdelete,wings3,wrtloctens,wrtout,xcomm_world
!!      xme_whoiam,xsum_master
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,&
&  mkmems,mpi_enreg,npwtot,&
&  nspinor,occ,pawang,pawrad,pawtab,psps,walli,xred)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes

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
 use interfaces_13paw
 use interfaces_13psp
 use interfaces_13recipspace
 use interfaces_13xc
 use interfaces_14iowfdenpot
 use interfaces_14occeig
 use interfaces_15common
 use interfaces_16response
 use interfaces_18seqpar, except_this_one => respfn
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 integer,intent(inout) :: iexit,nspinor
 real(dp),intent(in) :: cpui,walli
 real(dp),intent(out) :: etotal
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 integer,intent(in) :: mkmems(3)
 integer,intent(inout) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!---- Local variables : integer scalars
#if !defined MPI
 integer,parameter :: nkpt_max=50
#endif
#if defined MPI
 integer,parameter :: nkpt_max=-1
#endif
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!1   for wavefunction file, old format (version prior to 2.0)
!2   for wavefunction file, new format (version 2.0 and after)    (fform)
!(51 or 52   for density rho(r)       (fformr)
!101 or 102 for potential V(r) file. (fformv)
 integer,parameter :: fform=2,fformv=102,formeig=0,level=13
 integer,parameter :: response=1,syuse=0
 integer,save :: nsym_old=-1
 integer,save :: paw_gencond(6)=(/-1,-1,-1,-1,-1,-1/)
 integer :: fformr=52
 integer :: accessfil,analyt,ask_accurate,band_index,bantot,choice,cplex,flag,fullinit
 integer :: gscase,iapp,iatom,iband,idir,idir2,ider,ierr,ifft,ii,ikpt,ilmn,index
 integer :: initialized,ios,ipert,ipert2,ir,ireadwf0,iscf,iscf_eff,ispden,isppol
 integer :: itime,itypat,jband,jj,master,me,mgfftf,mk1mem,mkqmem,mpert,mpsang,mu
 integer :: n3xccc,nband_k,nblok,ncpgr,nfftf,nfftftot,nfftot,nhatdim,nhatgrdim
 integer :: nkpt_eff,nkxc,ntypat,nzlmopt,openexit,option,optgr0,optgr1,optgr2,optorth,optrad
 integer :: optatm,optdyfr,optgr,optn,optn2,optstr,optv
 integer :: outd2,prtbbb,prtdos,prtvol,psp_gencond,qzero,rdwr,rdwrpaw,rfasr,rfelfd,rfmgfd,rfphon,rfstrs
 integer :: rfuser,spaceworld,sumg0,tim_mkrho,timrev,usecprj,usexcnhat,v_size,vrsddb
 real(dp) :: boxcut,compch_fft,compch_sph,cpus,ecore,ecut_eff,ecutdg_eff,ecutf,eei,eew,ehart,eii,ek,enl,entropy,enxc
 real(dp) :: epaw,epawdc,etot,fermie,gsqcut,gsqcut_eff,gsqcutc_eff,qphnrm,residm
 real(dp) :: rhosum,tmp,tolmxf,tolwfr
 real(dp) :: ucvol,vxcavg
 logical :: test_mpi
 character(len=fnlen) :: ddbnm,dscrpt
 character(len=4) :: tag
 character(len=500) :: message
 type(bandstructure_type) :: bstruct
 type(dens_sym_operator_type) :: densymop_gs
 type(hdr_type) :: hdr
 type(pawfgr_type) :: pawfgr
 type(wffile_type) :: wffgs,wfftgs
 type(wvl_data) :: wvl
 integer :: ddkfil(3),ngfft(18),ngfftf(18),rfdir(3)
 integer,allocatable :: atindx(:),atindx1(:),blkflg(:,:,:,:),blkflg1(:,:,:,:)
 integer,allocatable :: blkflg2(:,:,:,:),carflg(:,:,:,:),dimcprj(:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),kgq(:,:),nattyp(:),npwarr(:)
 integer,allocatable :: pertsy(:,:),rfpert(:),symq(:,:,:),symrec(:,:,:)
 real(dp) :: dummy6(6),gmet(3,3),gprimd(3,3),k0(3),qphon(3)
 real(dp) :: rmet(3,3),rprimd(3,3),strsxc(6),tsec(2)
 real(dp),allocatable :: amass(:),cg(:,:),d2bbb(:,:,:,:,:,:),d2cart(:,:,:,:,:)
 real(dp),allocatable :: d2cart_bbb(:,:,:,:,:,:),d2eig0(:,:,:,:,:)
 real(dp),allocatable :: d2k0(:,:,:,:,:),d2lo(:,:,:,:,:),d2loc0(:,:,:,:,:)
 real(dp),allocatable :: d2matr(:,:,:,:,:),d2nfr(:,:,:,:,:),d2nl(:,:,:,:,:)
 real(dp),allocatable :: d2nl0(:,:,:,:,:),d2nl1(:,:,:,:,:),d2tmp(:,:,:,:,:)
 real(dp),allocatable :: d2vn(:,:,:,:,:),displ(:),doccde(:),dummy(:),dyew(:,:,:,:,:)
 real(dp),allocatable :: dyewq0(:,:,:),dyfrlo(:,:,:),dyfrlo_indx(:,:,:)
 real(dp),allocatable :: dyfrnl(:,:,:),dyfrwf(:,:,:),dyfrx1(:,:,:,:,:)
 real(dp),allocatable :: dyfrx2(:,:,:),eigen0(:),eigval(:),eigvec(:)
 real(dp),allocatable :: eltcore(:,:),elteew(:,:),eltfrhar(:,:),eltfrkin(:,:)
 real(dp),allocatable :: eltfrloc(:,:),eltfrnl(:,:),eltfrxc(:,:),grtn_indx(:,:)
 real(dp),allocatable :: grxc(:,:),grxc_indx(:,:),kxc(:,:),nhat(:,:),nhatgr(:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph1df(:,:),phfrq(:),phnons(:,:,:)
 real(dp),allocatable :: rhog(:,:),rhor(:,:),rhowfg(:,:),rhowfr(:,:)
 real(dp),allocatable :: vhartr(:),vpsp(:),vpsp1(:,:),vtrial(:,:)
 real(dp),allocatable :: vxc(:,:),work(:),xccc3d(:),ylm(:,:),ylmgr(:,:,:)
 character(len=fnlen) :: tmpfil(15)
 type(cprj_type),allocatable :: cprj_dum(:,:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)
 type(pawrhoij_type),allocatable :: pawrhoij(:)
 integer :: idtmpfil(15)

! local GIPAW variables
 real(dp),allocatable :: cs(:,:,:),gcart(:,:,:,:),jvec(:,:,:)
 type(gipaw_type),allocatable :: gipaw_aug(:)

! ***********************************************************************

!DEBUG
!write(6,*)' respfn : enter'
!stop
!ENDDEBUG

 call timab(132,1,tsec)
 call timab(133,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'enter         ')

 mpi_enreg%paralbd=1
 mpi_enreg%parareel=0
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1
 mpi_enreg%paral_fft=0
 mpi_enreg%paral_level=2

#if defined MPI
           allocate(mpi_enreg%proc_distrb(dtset%nkpt,dtset%mband,dtset%nsppol))
           call distrb2(dtset%mband, dtset%nband, dtset%nkpt, dtset%nsppol, mpi_enreg)
#endif


 call initmpi_fft(dtset,mpi_enreg)

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
   do ii=1,dtset%nfft;pawfgr%coatofin(ii)=ii;pawfgr%fintocoa(ii)=ii;end do
  end if
  pawfgr%natom=dtset%natom
  pawfgr%nfftc=dtset%nfft;pawfgr%mgfftc=dtset%mgfft;pawfgr%ngfftc(:)=dtset%ngfft(:)
  pawfgr%nfft =nfftf     ;pawfgr%mgfft=mgfftf      ;pawfgr%ngfft(:)=ngfftf(:)
  ecutdg_eff = dtset%pawecutdg * (dtset%dilatmx)**2
  ecut_eff   = dtset%ecut      * (dtset%dilatmx)**2
  ecutf=dtset%ecut;if (pawfgr%usefinegrid==1) ecutf=dtset%pawecutdg
 else
  pawfgr%usefinegrid=0
  nfftf=dtset%nfft;mgfftf=dtset%mgfft;ngfftf(:)=dtset%ngfft(:)
  allocate(pawfgr%coatofin(0),pawfgr%fintocoa(0))
  ecut_eff= dtset%ecut * (dtset%dilatmx)**2
  ecutdg_eff=ecut_eff
  ecutf=dtset%ecut
 end if

!
!  If dtset%accesswff == 2 set all array outputs to netcdf format
!
 accessfil = 0
 if (dtset%accesswff == 2) then
   accessfil = 1
 end if
 if (dtset%accesswff == 3) then
   accessfil = 3
 end if

!Structured debugging if prtvol==-level
 prtvol=dtset%prtvol
 if(prtvol==-level)then
  write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,&
&   ' respfn : enter , debug mode '
  call wrtout(06,message,'COLL')
 end if

!Option input variables
 iscf=dtset%iscf


!Respfn input variables
 rfasr=dtset%rfasr   ; rfdir(1:3)=dtset%rfdir(1:3)
 rfelfd=dtset%rfelfd ; rfphon=dtset%rfphon
 rfstrs=dtset%rfstrs ; rfuser=dtset%rfuser
 rfmgfd=dtset%rfmgfd

!mkmem variables (mkmem is already argument)
 mkqmem=mkmems(2) ; mk1mem=mkmems(3)

 ntypat=psps%ntypat

 call status(0,dtfil%filstat,iexit,level,'call setup1   ')

 ecore=zero

!LIKELY TO BE TAKEN AWAY
 initialized=0
 ek=zero ; ehart=zero ; enxc=zero ; eei=zero ; enl=zero
 eii=zero ; eew=zero

!Set up for iterations
 allocate(amass(dtset%natom))
!DEBUG
!write(6,*)' respfn : ntypat=',ntypat
!ENDDEBUG
 call setup1(dtset%acell_orig,amass,dtset%amu,bantot,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,dtset%iboxcut,dtset%intxc,dtset%ionmov,&
& dtset%natom,dtset%nband,ngfftf,ngfft,dtset%nkpt,dtset%nqpt,dtset%nsppol,dtset%nsym,ntypat,&
& dtset%qptn,response,rmet,dtset%rprim_orig,rprimd,dtset%typat,ucvol,psps%usepaw)

!Define the set of admitted perturbations
 mpert=dtset%natom+7

!Initialize the list of perturbations rfpert
 allocate(rfpert(mpert))
 rfpert(:)=0
 if(rfphon==1)rfpert(dtset%rfatpol(1):dtset%rfatpol(2))=1

 if(rfelfd==1.or.rfelfd==2)rfpert(dtset%natom+1)=1
 if(rfelfd==1.or.rfelfd==3)rfpert(dtset%natom+2)=1

 if(rfstrs==1.or.rfstrs==3)rfpert(dtset%natom+3)=1
 if(rfstrs==2.or.rfstrs==3)rfpert(dtset%natom+4)=1

! in the magnetic field case, rfmgfd = 2 is exactly the same
! as rfelfd = 2, that is, compute the DDK perturbation
 if(rfmgfd==1.or.rfmgfd==2)rfpert(dtset%natom+1)=1
 if(rfmgfd==1.or.rfmgfd==3)rfpert(dtset%natom+5)=1

 if(rfuser==1.or.rfuser==3)rfpert(dtset%natom+6)=1
 if(rfuser==2.or.rfuser==3)rfpert(dtset%natom+7)=1

!Create names for the temporary files based on dtfil%filnam_ds(5)
!by appending adequate string.
!'_1WF1' -> dtfil%unwft1
!'_1WF2' -> dtfil%unwft2
!'_KG'   -> dtfil%unkg
!'_KGQ'  -> dtfil%unkgq (not used for the time being)
!'_KG1'  -> dtfil%unkg1
!'_DUM'  -> tmp_unit (real dummy name)
!'_WFGS' -> dtfil%unwftgs
!'_WFKQ' -> dtfil%unwftkq
!'_YLM'  -> dtfil%unylm
!'_YLM1' -> dtfil%unylm1
!'_PAW'  -> dtfil%unpaw
!'_PAW1' -> dtfil%unpaw1
!'_PAWQ' -> dtfil%unpawq
 tmpfil(1) =trim(dtfil%filnam_ds(5))//'_1WF1'
 tmpfil(2) =trim(dtfil%filnam_ds(5))//'_1WF2'
 tmpfil(3) =trim(dtfil%filnam_ds(5))//'_KG'
 tmpfil(4) =trim(dtfil%filnam_ds(5))//'_KGQ'
 tmpfil(5) =trim(dtfil%filnam_ds(5))//'_KG1'
 tmpfil(6) =trim(dtfil%filnam_ds(5))//'_DUM'
 tmpfil(7) =trim(dtfil%filnam_ds(5))//'_WFGS'
 tmpfil(8) =trim(dtfil%filnam_ds(5))//'_WFKQ'
 tmpfil(10)=trim(dtfil%filnam_ds(5))//'_YLM'
 tmpfil(11)=trim(dtfil%filnam_ds(5))//'_YLM1'
 tmpfil(12)=trim(dtfil%filnam_ds(5))//'_PAW'
 tmpfil(13)=trim(dtfil%filnam_ds(5))//'_PAW1'
 tmpfil(14)=trim(dtfil%filnam_ds(5))//'_PAWQ'

!Default for sequential use
 master=0
 call xme_whoiam(me)
!Init spaceworld

!BEGIN TF_CHANGES
 call xcomm_world(mpi_enreg,spaceworld)
!END TF_CHANGES


!There is a definite problem for the treatment
!of // within CPP sections...
#if defined MPI
           do ii=1,11
              idtmpfil(ii)=ii
           end do
           if (mpi_enreg%paral_compil_mpio==1.and.dtset%accesswff==1)then
              idtmpfil(1)=0              !_1wf1
              idtmpfil(2)=0              ! s1wf2
              idtmpfil(7)=0              !  WFGS
              idtmpfil(8)=0              !  WFKQ
           endif
           call int2char4(me,tag)
           do ii=1,11
           if(idtmpfil(ii) /= 0) tmpfil(ii)=trim(tmpfil(ii))//'_P-'//tag
           end do
#endif

 call status(0,dtfil%filstat,iexit,level,'call kpgio(1) ')

!Set up the basis sphere of planewaves
 allocate(kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt))
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,tmpfil(3),&
& dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,dtset%mpw,npwarr,npwtot,&
& dtset%nsppol,dtfil%unkg)

!Set up the Ylm for each k point
 mpsang=psps%mpsang
 allocate(ylm(dtset%mpw*dtset%mkmem,mpsang*mpsang*psps%useylm))
 if (rfstrs/=0) then
  allocate(ylmgr(dtset%mpw*dtset%mkmem,9,mpsang*mpsang*psps%useylm))
 else
  allocate(ylmgr(1,1,psps%useylm))
 end if
 if (psps%useylm==1) then
  if(dtset%mkmem==0) open(dtfil%unylm,file=tmpfil(10),form='unformatted',status='unknown')
  call status(0,dtfil%filstat,iexit,level,'call initylmg ')
  option=0
  if (rfstrs/=0) option=2
  call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&               npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
 end if

 call timab(133,2,tsec)
 call timab(134,1,tsec)

!Open and read pseudopotential files
 call status(0,dtfil%filstat,iexit,level,'call pspini(1)')
 call pspini(dtset,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,level,&
& pawrad,pawtab,psps,rprimd)

 call timab(134,2,tsec)
 call timab(135,1,tsec)

!Initialize band structure datatype
 allocate(doccde(bantot),eigen0(bantot))
 doccde(:)=zero ; eigen0(:)=zero
 call bstruct_init(bantot,bstruct,doccde,eigen0,dtset%istwfk,dtset%kptns,&
& dtset%nband,dtset%nkpt,npwarr,dtset%nsppol,dtset%occ_orig,dtset%wtk)
 deallocate(doccde,eigen0)

!Initialize PAW atomic occupancies
 if (psps%usepaw==1) then
  allocate(pawrhoij(dtset%natom))
  call initrhoij(psps%indlmn,dtset%lexexch,psps%lmnmax,dtset%lpawu,dtset%natom,dtset%nspden,&
&                dtset%nsppol,dtset%ntypat,pawrhoij,pawtab,dtset%spinat,dtset%typat)
 end if

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
 call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&                residm,rprimd,occ,pawrhoij,psps%usepaw,xred)

!Clean band structure datatype (should use it more in the future !)
 call bstruct_clean(bstruct)

 call status(0,dtfil%filstat,iexit,level,'call inwffil(1')

!Initialize wavefunction files and wavefunctions.
 ireadwf0=1

 allocate(cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol))
 allocate(eigen0(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen0(:)=zero ; ask_accurate=1
 optorth=1;if (psps%usepaw==1.and.ireadwf0==1) optorth=0

 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
& formeig,gmet,hdr,ireadwf0,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,&
& dtset%nband,ngfft,dtset%nkpt,npwarr,dtset%nspden,nspinor,dtset%nsppol,dtset%nsym,&
& occ,optorth,psps,prtvol,rprimd,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wffgs,wfftgs,dtfil%unwffgs,dtfil%unwftgs,dtfil%fnamewffk,tmpfil(7),wvl)

 if (psps%usepaw==1.and.ireadwf0==1) then
  call rhoij_copy(hdr%pawrhoij,pawrhoij)
 end if

 call timab(135,2,tsec)
 call timab(136,1,tsec)

!Close wffgs, if it was ever opened (in inwffil)
 if (ireadwf0==1) then
  call WffClose(wffgs,ierr)
 end if

!Report on eigen0 values
 write(message, '(a,a)' )ch10,' respfn : eigen0 array'
 call wrtout(6,message,'COLL')
 nkpt_eff=dtset%nkpt
 if( (prtvol==0.or.prtvol==1.or.prtvol==2) .and. dtset%nkpt>nkpt_max ) nkpt_eff=nkpt_max
 band_index=0
 do isppol=1,dtset%nsppol
  do ikpt=1,dtset%nkpt
   nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
   if(ikpt<=nkpt_eff)then
    write(message, '(a,i2,a,i5)' )&
&    '  isppol=',isppol,', k point number',ikpt
    call wrtout(6,message,'COLL')
    do iband=1,nband_k,4
     write(message, '(a,4es16.6)')&
&     '  ',eigen0(iband+band_index:min(iband+3,nband_k)+band_index)
     call wrtout(6,message,'COLL')
    end do
   else if(ikpt==nkpt_eff+1)then
    write(message,'(a,a)' )&
&    '  respfn : prtvol=0, 1 or 2, stop printing eigen0.',ch10
    call wrtout(6,message,'COLL')
   end if
   band_index=band_index+nband_k
  end do
 end do

!Allocation for forces and atomic positions (should be taken away, also argument ... )
 allocate(grxc(3,dtset%natom))

!Do symmetry stuff
 call status(0,dtfil%filstat,iexit,level,'call setsym   ')
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 allocate(irrzon(nfftot**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol))
 allocate(phnons(2,nfftot**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol))
!DEBUG
!write(6,*)' respfn : stop before GS setsym'
!stop
!ENDDEBUG
 allocate(indsym(4,dtset%nsym,dtset%natom),symrec(3,3,dtset%nsym))
 iscf_eff=0
!If the density is to be computed by mkrho, need irrzon and phnons
 if(dtset%getden==0)iscf_eff=1
 call setsym(densymop_gs,indsym,irrzon,iscf_eff,dtset%natom,&
& nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
& phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!DEBUG
!write(6,*)' respfn : stop after GS setsym'
!stop
!ENDDEBUG

!Symmetrize atomic coordinates over space group elements:
 call status(0,dtfil%filstat,iexit,level,'call symzat   ')
 call symzat(indsym,dtset%natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

!Generate an index table of atoms, in order for them to be used
!type after type.

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
 allocate(ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom))
 allocate(ph1df(2,3*(2*mgfftf+1)*dtset%natom))
 call status(0,dtfil%filstat,iexit,level,'call getph    ')
 call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
  call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 else
  ph1df(:,:)=ph1d(:,:)
 end if

!Compute occupation numbers and fermi energy, in case
!occupation scheme is metallic.
 allocate(doccde(dtset%mband*dtset%nkpt*dtset%nsppol))
 if( dtset%occopt>=3.and.dtset%occopt<=7 ) then
  call status(0,dtfil%filstat,iexit,level,'call newocc   ')
  call newocc(doccde,eigen0,entropy,fermie,dtset%fixmom,dtset%mband,dtset%nband,&
&  dtset%nelect,dtset%nkpt,nspinor,dtset%nsppol,occ,dtset%occopt,prtvol,dtset%stmbias,&
&  dtset%tphysel,dtset%tsmear,dtset%wtk)

! Update fermie and occ
  etot=hdr%etot ; residm=hdr%residm
  call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&  residm,rprimd,occ,pawrhoij,psps%usepaw,xred)

 else
! doccde is irrelevant in this case
  doccde(:)=zero

 end if

!Recompute first large sphere cut-off gsqcut,
!without taking into account dilatmx
 if (psps%usepaw==1) then
  write(message,'(2a)') ch10,' FFT (fine) grid used in SCF cycle:'
  call wrtout(6,message,'COLL')
 end if
 k0=zero
 call status(0,dtfil%filstat,iexit,level,'call getcut   ')
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,6,k0,ngfftf)

!PAW: 1- Initialize values for several arrays depending only on atomic data
!     2- Check overlap
!     3- Identify FFT points in spheres and compute g_l(r).Y_lm(r)
!     4- Allocate PAW specific arrays
!     5- Compute perturbed local potential inside spheres
!     6- Eventually open temporary storage files
 if(psps%usepaw==1) then
! 1-Initialize values for several arrays depending only on atomic data
  if (psp_gencond==1.or.&
&     paw_gencond(1)/=dtset%pawlcutd .or.paw_gencond(2)/=dtset%pawlmix  .or.&
&     paw_gencond(3)/=dtset%pawnphi  .or.paw_gencond(4)/=dtset%pawntheta.or.&
&     paw_gencond(5)/=dtset%pawspnorb.or.paw_gencond(6)/=dtset%pawxcdev) then
   call timab(553,1,tsec)
   call pawinit(0._dp,psps%indlmn,dtset%pawlcutd,dtset%pawlmix,psps%lmnmax,psps%mpsang,psps%n1xccc,&
&       dtset%pawnphi,dtset%nsym,dtset%pawntheta,psps%ntypat,&
&       pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev)
   paw_gencond(1)=dtset%pawlcutd ; paw_gencond(2)=dtset%pawlmix
   paw_gencond(3)=dtset%pawnphi  ; paw_gencond(4)=dtset%pawntheta
   paw_gencond(5)=dtset%pawspnorb; paw_gencond(6)=dtset%pawxcdev
   call timab(553,2,tsec)
  end if
  if (psp_gencond==1.or.nsym_old/=dtset%nsym) then
   call setsymrhoij(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,&
&                   rprimd,dtset%symafm,symrec,pawang%zarot)
   nsym_old=dtset%nsym
  end if
  if (dtset%usepawu>0.or.dtset%useexexch>0) then
   call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%jpawu,dtset%lexexch,dtset%lpawu,&
&       psps%indlmn,psps%lmnmax,ntypat,pawang,dtset%pawprtvol,pawrad,pawtab,dtset%upawu,&
&       dtset%useexexch,dtset%usepawu)
  end if
  compch_fft=-1.d5;compch_sph=-1.d5
  usexcnhat=maxval(pawtab(:)%vlocopt)
  usecprj=dtset%pawusecp
! 2-Check overlap
  call status(0,dtfil%filstat,iexit,level,'call chkpawovlp')
  call chkpawovlp(dtset%natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,xred)
! 3-Identify FFT points in spheres and compute g_l(r).Y_lm(r)
  allocate(pawfgrtab(dtset%natom))
  do iatom=1,dtset%natom
   pawfgrtab(iatom)%l_size=pawtab(dtset%typat(iatom))%lcut_size
   pawfgrtab(iatom)%nfgd=0;allocate(pawfgrtab(iatom)%ifftsph(0))
   pawfgrtab(iatom)%rfgd_allocated=0;allocate(pawfgrtab(iatom)%rfgd(0,0))
   pawfgrtab(iatom)%gylm_allocated=0;allocate(pawfgrtab(iatom)%gylm(0,0))
   pawfgrtab(iatom)%gylmgr_allocated=0;allocate(pawfgrtab(iatom)%gylmgr(0,0,0))
   pawfgrtab(iatom)%gylmgr2_allocated=0;allocate(pawfgrtab(iatom)%gylmgr2(0,0,0))
  end do
  optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
  if (dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0) optgr1=dtset%pawstgylm
  if (rfphon==1) then
   optgr1=dtset%pawstgylm;optgr2=dtset%pawstgylm
  end if
  if (rfuser==1) then ! for GIPAW, need r-R data around spheres
   optrad=1
  end if
  call status(0,dtfil%filstat,iexit,level,'call nhatgrid ')
  call nhatgrid(atindx1,gmet,mpi_enreg,dtset%natom,nattyp,nfftf,ngfftf,psps%ntypat,&
&  optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,dtset%typat,ucvol,xred)
! 4-Allocate PAW specific arrays
  allocate(paw_ij(dtset%natom),paw_an(dtset%natom),dimcprj(dtset%natom))
  do iatom=1,dtset%natom
   itypat=dtset%typat(iatom)
   dimcprj(iatom)=pawtab(itypat)%lmn_size
   paw_an(iatom)%angl_size=pawang%angl_size
   paw_an(iatom)%mesh_size=pawtab(itypat)%mesh_size
   paw_an(iatom)%nspden   =dtset%nspden
   paw_an(iatom)%lm_size  =pawtab(itypat)%lcut_size**2
   paw_an(iatom)%has_vxcval=0
   allocate(paw_an(iatom)%lmselect(paw_an(iatom)%lm_size))
   paw_ij(iatom)%cplex     =cplex
   paw_an(iatom)%has_vxcval=0
   paw_ij(iatom)%cplex_dij=max(cplex,1+dtset%pawspnorb,nspinor)
   paw_ij(iatom)%nspden   =dtset%nspden
   paw_ij(iatom)%nsppol   =dtset%nsppol
   paw_ij(iatom)%lmn_size =pawtab(itypat)%lmn_size
   paw_ij(iatom)%lmn2_size=pawtab(itypat)%lmn2_size
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
   allocate(paw_ij(iatom)%dij(paw_ij(iatom)%cplex_dij*pawtab(itypat)%lmn2_size,paw_ij(iatom)%ndij))
   if (pawtab(itypat)%usepawu>0) then
    allocate(paw_ij(iatom)%noccmmp(2*pawtab(itypat)%lpawu+1,&
&                   2*pawtab(itypat)%lpawu+1,dtset%nspden))
    allocate(paw_ij(iatom)%nocctot(dtset%nspden))
   end if
   if (pawtab(itypat)%useexexch>0) then
    allocate(paw_ij(iatom)%vpawx(1,pawtab(itypat)%lmn_size,dtset%nspden))
   end if
  end do
! 5-Compute perturbed local potential inside spheres
  allocate(vpsp1(nfftf,3))
  optv=1;optn=0;optn2=1;cplex=1;idir=0;qphon(:)=zero
  do iatom=1,dtset%natom
   pawfgrtab(iatom)%vlocgr_allocated=1
   allocate(pawfgrtab(iatom)%vlocgr(3,pawfgrtab(iatom)%nfgd))
   call atm2fft3(atindx,dummy,vpsp1(:,:),cplex,dummy,gmet,gprimd,gsqcut,idir,iatom,&
&                mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,nattyp,3,nfftf,ngfftf,&
&                ntypat,optn,optn2,optv,dtset%paral_kgb,pawtab,ph1df,psps%qgrid_vl,&
&                qphon,dtset%typat,ucvol,psps%usepaw,psps%vlspl,xred)
   do ii=1,pawfgrtab(iatom)%nfgd
    jj=pawfgrtab(iatom)%ifftsph(ii)
    if (dtset%useric==-3) then
     do idir=1,3
      pawfgrtab(iatom)%vlocgr(idir,ii)=-gprimd(idir,1)*vpsp1(jj,1)&
&           -gprimd(idir,2)*vpsp1(jj,2)-gprimd(idir,3)*vpsp1(jj,3)
     end do
    end if
   end do
  end do
  deallocate(vpsp1)
! 6-Eventually open temporary storage files
  if(dtset%mkmem==0) then
   open(dtfil%unpaw ,file=tmpfil(12),form='unformatted',status='unknown')
   rewind(unit=dtfil%unpaw)
  end if

 else ! PAW vs NCPP
  usexcnhat=0;usecprj=0
 end if

 allocate(rhog(2,nfftf),rhor(nfftf,dtset%nspden))
!Read ground-state charge density from diskfile in case getden /= 0
!or compute it from wfs that were
!read previously : rhor as well as rhog

 if (dtset%getden /= 0) then

  if (me == 0) then
   rdwr=1;rdwrpaw=psps%usepaw;if(ireadwf0/=0) rdwrpaw=0
   call status(0,dtfil%filstat,iexit,level,'call ioarr    ')
   call ioarr (accessfil,rhor, dtset, etotal,fformr,dtfil%fildensin,hdr, mpi_enreg, &
        & nfftf,pawrhoij,rdwr,rdwrpaw,ngfft)
   if (psps%usepaw==1.and.rdwrpaw/=0) then
    call hdr_update(bantot,etot,fermie,hdr,dtset%natom,&
&                   residm,rprimd,occ,pawrhoij,psps%usepaw,xred)
   end if
  end if

#if defined MPI
          call MPI_BARRIER(spaceworld,ierr)
          call MPI_BCAST(rhor,nfftf*dtset%nspden,&
&            MPI_DOUBLE_PRECISION,0,spaceworld,ierr)
#endif

! Compute up+down rho(G) by fft
  call status(0,dtfil%filstat,iexit,level,'call fourdp   ')
  allocate(work(nfftf))
  work(:)=rhor(:,1)
  call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
  deallocate(work)

 else

! Obtain the charge density from read wfs
! Be careful: in PAW, compensation density has to be added !
  call status(0,dtfil%filstat,iexit,level,'call mkrho    ')
  tim_mkrho=4
  if (psps%usepaw==1) then
   allocate(rhowfg(2,dtset%nfft),rhowfr(dtset%nfft,dtset%nspden))
   call mkrho(cg,densymop_gs,dtset,irrzon,kg,&
&   mpi_enreg,npwarr,nspinor,occ,phnons,rhowfg,rhowfr,tim_mkrho,ucvol,&
&   dtfil%unkg,wfftgs,wvl%wfs)
   call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
   deallocate(rhowfg,rhowfr)
  else
   call mkrho(cg,densymop_gs,dtset,irrzon,kg,&
&   mpi_enreg,npwarr,nspinor,occ,phnons,rhog,rhor,tim_mkrho,ucvol,&
&   dtfil%unkg,wfftgs,wvl%wfs)
  end if

 end if    ! getden

!In PAW, compensation density has eventually to be added
 nhatgrdim=0;nhatdim=0
 if (psps%usepaw==1.and. &
&  ((usexcnhat==0).or.(dtset%getden==0).or.dtset%xclevel==2)) then
  nhatdim=1;allocate(nhat(nfftf,dtset%nspden))
  call timab(558,1,tsec)
  nhatgrdim=0;if (dtset%xclevel==2.and.dtset%pawnhatxc>0) nhatgrdim=usexcnhat
  ider=2*nhatgrdim
  if (nhatgrdim>0) allocate(nhatgr(nfftf,dtset%nspden,3))
   call pawmknhat(compch_fft,ider,0,mpi_enreg,dtset%natom,nfftf,ngfftf,nhatgrdim,dtset%nspden,&
&       psps%ntypat,dtset%paral_kgb,pawang,pawfgrtab,nhatgr,nhat,pawrhoij,pawtab,dtset%typat,ucvol)
   if (dtset%getden==0) then
    rhor(:,:)=rhor(:,:)+nhat(:,:)
    call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
   end if
  call timab(558,2,tsec)
 end if

!The GS irrzon and phnons were only needed to symmetrize the GS density
 deallocate(irrzon,phnons)

!Debugging : print the different parts of rhor
 if(prtvol==-level)then
  write(message,'(a)') '   ir     rhor(ir)     '
  call wrtout(06,message,'COLL')
  do ir=1,nfftf
   if(ir<=11 .or. mod(ir,301)==0 )then
    write(message,'(i5,a,es13.6)')ir,' ',rhor(ir,1)
    call wrtout(06,message,'COLL')
    if(dtset%nspden==2)then
     write(message,'(a,es13.6)')'      ',rhor(ir,2)
     call wrtout(06,message,'COLL')
    end if
   end if
  end do
 end if

!Will compute now the total potential
 allocate(vhartr(nfftf),vpsp(nfftf),vtrial(nfftf,dtset%nspden),vxc(nfftf,dtset%nspden))

!Compute local ionic pseudopotential vpsp
!    and core electron density xccc3d:
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 allocate(xccc3d(n3xccc))
 if (psps%usepaw==1) then
! PAW: compute Vloc and core charge together in reciprocal space
  call timab(562,1,tsec)
  optatm=1;optdyfr=0;optgr=0;optstr=0;optv=1;optn=n3xccc/nfftf;optn2=1
  call atm2fft(atindx1,xccc3d,vpsp,dummy,dummy,eei,dummy,gmet,gprimd,dummy,dummy,gsqcut,&
&              mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,ntypat,&
&              optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&              pawtab,ph1df,psps%qgrid_vl,dtset%qprtrb,dummy,dummy6,dummy6,&
&              ucvol,psps%usepaw,dummy,dtset%vprtrb,psps%vlspl)
  call timab(562,2,tsec)
 else
! Norm-cons.: compute Vloc in reciprocal space and core charge in real space
  option=1
  allocate(dyfrlo_indx(3,3,dtset%natom),grtn_indx(3,dtset%natom))  !Dummy variables
  call status(0,dtfil%filstat,iexit,level,'call mklocl   ')
  call mklocl(dtset,dyfrlo_indx,eei,gmet,gprimd,&
&  grtn_indx,gsqcut,dummy6,mgfftf,mpi_enreg,dtset%natom,nattyp,&
&  nfftf,ngfftf,dtset%nspden,ntypat,option,ph1df,psps,&
&  dtset%qprtrb,rhog,rhor,rmet,rprimd,ucvol,dtset%vprtrb,vpsp,xred)
  deallocate(dyfrlo_indx,grtn_indx)
  if (psps%n1xccc/=0) then
   call status(0,dtfil%filstat,iexit,level,'call mkcore   ')
   allocate(dyfrx2(3,3,dtset%natom))  !Dummy variable
   call mkcore(dummy6,dyfrx2,grxc,mpi_enreg,dtset%natom,nfftf,dtset%nspden,ntypat,&
&   ngfftf(1),psps%n1xccc,ngfftf(2),ngfftf(3),option,rprimd,dtset%typat,ucvol,vxc,&
&   psps%xcccrc,psps%xccc1d,xccc3d,xred)
   deallocate(dyfrx2)
  end if
 end if

!Set up hartree and xc potential. Compute kxc here.
 option=2
 nkxc=2*dtset%nspden-1
 if(dtset%xclevel==2)nkxc=23
 allocate(kxc(nfftf,nkxc))
 call status(0,dtfil%filstat,iexit,level,'call rhohxc   ')
 call rhohxc(dtset,enxc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,dtset%nspden,n3xccc,option,rhog,rhor,&
& rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d)

!Compute local + Hxc potential, and subtract mean potential.
 do ispden=1,min(dtset%nspden,2)
  do ifft=1,nfftf
   vtrial(ifft,ispden)=vhartr(ifft)+vxc(ifft,ispden)+vpsp(ifft)
  end do
 end do
 if (dtset%nspden==4) then
  do ispden=3,4
   do ifft=1,nfftf
    vtrial(ifft,ispden)=vxc(ifft,ispden)
   end do
  end do
 end if
 deallocate(vhartr)

 call status(0,dtfil%filstat,iexit,level,'end respfn(1) ')

 if(prtvol==-level)then
  write(message,'(a,a)') ch10,&
&   ' respfn : ground-state density and potential set up. '
  call wrtout(06,message,'COLL')
 end if

!PAW: compute Dij quantities (psp strengths)
 if (psps%usepaw==1)then
  option=1;nzlmopt=0;if (dtset%pawnzlm>0) nzlmopt=-1
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
  call status(0,dtfil%filstat,iexit,level,'call pawdenpot')
  call pawdenpot(compch_sph,epaw,epawdc,dtset%ixc,dtset%natom,dtset%nspden,ntypat,&
&                nzlmopt,option,paw_an,paw_ij,pawang,dtset%pawprtvol,&
&                pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,&
&                dtset%typat,dtset%xclevel,psps%znuclpsp)
  call status(0,dtfil%filstat,iexit,level,'call pawdij   ')
  call pawdij(dtset,dtset%enunit,mpi_enreg,dtset%natom,nfftf,ngfftf,dtset%nspden,ntypat,&
&             paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,pawrad,dtset%pawspnorb,pawtab,&
&             dtset%pawxcdev,dtset%typat,ucvol,vtrial,vxc)
  call symdij(psps%indlmn,indsym,psps%lmnmax,dtset%natom,dtset%nsym,ntypat,paw_ij,pawang,&
&             dtset%pawprtvol,dtset%symafm,symrec,dtset%typat)
  do iatom=1,dtset%natom
   deallocate(paw_ij(iatom)%dijhartree,paw_an(iatom)%vxc1,paw_an(iatom)%vxct1)
   paw_ij(iatom)%has_dijhartree=0
   if (dtset%pawspnorb>0) deallocate(paw_an(iatom)%vh1)
   if (pawtab(itypat)%useexexch>0) deallocate(paw_an(iatom)%vxc_ex)
  end do
 end if

! -----2. Frozen-wavefunctions and Ewald(q=0) parts of 2DTE

 allocate(eltcore(6,6),elteew(6+3*dtset%natom,6),eltfrhar(6,6),eltfrnl(6+3*dtset%natom,6))
 allocate(eltfrloc(6+3*dtset%natom,6),eltfrkin(6,6),eltfrxc(6+3*dtset%natom,6))
 eltcore(:,:)=zero
 elteew(:,:)=zero
 eltfrnl(:,:)=zero
 eltfrloc(:,:)=zero
 eltfrkin(:,:)=zero
 eltfrhar(:,:)=zero
 eltfrxc(:,:)=zero

 allocate(dyew(2,3,dtset%natom,3,dtset%natom),dyewq0(3,3,dtset%natom))
 allocate(dyfrnl(3,3,dtset%natom),dyfrlo(3,3,dtset%natom),dyfrwf(3,3,dtset%natom))
 allocate(dyfrx2(3,3,dtset%natom))
 dyew(:,:,:,:,:)=zero
 dyewq0(:,:,:)=zero
 dyfrnl(:,:,:)=zero
 dyfrlo(:,:,:)=zero
 dyfrwf(:,:,:)=zero
 dyfrx2(:,:,:)=zero

 if (rfphon==1) then

  call dyfnl3(atindx1,cg,cprj_dum,dimcprj,dyfrnl,dtset%ecut,dtset%ecutsm,&
&  eigen0,fform,indsym,dtset%istwfk,&
&  kg,dtset%kptns,dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,mpsang,&
&  dtset%mpw,dtset%natom,nattyp,dtset%nband,&
&  nfftf,ngfft,ngfftf,dtset%nkpt,dtset%nloalg,npwarr,&
&  dtset%nspden,nspinor,dtset%nsppol,dtset%nsym,ntypat,occ,&
&  paw_ij,pawang,dtset%pawprtvol,pawfgrtab,pawrhoij,pawtab,&
&  ph1d,prtvol,psps,rprimd,dtset%symafm,symrec,dtset%typat,&
&  dtfil%unkg,dtfil%unpaw,dtfil%unylm,&
&  0,wfftgs,vtrial,dtset%wtk,xred,ylm)

! No more need of these local derivatives
  if (psps%usepaw==1) then
   do iatom=1,dtset%natom
    if (associated(pawfgrtab(iatom)%gylmgr2))deallocate(pawfgrtab(iatom)%gylmgr2)
    pawfgrtab(iatom)%gylmgr2_allocated=0
   end do
  end if

! dyfrnl has not yet been symmetrized, but will be in the next routine
  call status(0,dtfil%filstat,iexit,level,'call dyfro3    ')
  call dyfro3(atindx1,dyfrnl,dyfrlo,dyfrwf,dyfrx2,&
&  gmet,gprimd,gsqcut,indsym,mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,nattyp,&
&  nfftf,ngfftf,dtset%nspden,dtset%nsym,ntypat,psps%n1xccc,n3xccc,dtset%paral_kgb,pawtab,ph1df,psps%qgrid_vl,&
&  rhog,rprimd,symrec,dtset%typat,ucvol,psps%usepaw,psps%vlspl,&
&  vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)

! The frozen-wavefunction part of the dynamical matrix is now:
!  dyfrnl:  non-local contribution
!  dyfrlo:  local contribution
!  dyfrx2:  2nd order xc core correction contribution
!  dyfrwf:  all      contributions
!  In case of PAW, it misses a term coming from the pertubed overlap operator

! Compute Ewald (q=0) contribution
  qphon(:)=zero
  sumg0=0
  call status(0,dtfil%filstat,iexit,level,'call ewald3(1)')
  call ewald3(dyew,gmet,dtset%natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,psps%ziontypat)
  option=1
  call q0dy3(dtset%natom,dyewq0,dyew,option)

!End of the frozen-wavefunction and Ewald(q=0) parts of the dynamical matrix
 end if

!Section for the strain perturbation - frozen-wavefunction, Ewald, etc.
! parts of the elastic tensor

 if(rfstrs/=0) then

! Verify that k-point set has full space-group symmetry; otherwise exit
  call status(0,dtfil%filstat,iexit,level,'call symkchk ')
  timrev=1
  call symkchk(gmet,dtset%kptns,dtset%nkpt,dtset%nsym,symrec,timrev)

! Calculate the nonlocal part of the elastic tensor
  call status(0,dtfil%filstat,iexit,level,'call eltfrnl3 ')
  call eltfrnl3(atindx,atindx1,cg,eltfrnl,dtset%ecut,dtset%ecutsm,&
&  fform,dtset%istwfk,&
&  kg,dtset%kptns,dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,mpsang,&
&  dtset%mpw,dtset%natom,nattyp,dtset%nband,&
&  dtset%nkpt,ngfft,dtset%nloalg,npwarr,&
&  nspinor,dtset%nsppol,dtset%nsym,ntypat,occ,ph1d,&
&  prtvol,psps,rprimd,dtset%typat,dtfil%unkg,wfftgs,dtfil%unylm,&
&  dtset%useylm,dtset%wtk,xred,ylm,ylmgr)
! Calculate the kinetic part of the elastic tensor
  call status(0,dtfil%filstat,iexit,level,'call eltfrkin3')
  call eltfrkin3(cg,eltfrkin,dtset%ecut,dtset%ecutsm,dtset%effmass,&
&  fform,dtset%istwfk,&
&  kg,dtset%kptns,dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,&
&  dtset%mpw,dtset%nband,&
&  dtset%nkpt,ngfft,npwarr,&
&  nspinor,dtset%nsppol,dtset%nsym,occ,&
&  prtvol,rprimd,dtfil%unkg,wfftgs,&
&  dtset%wtk)

! Calculate the hartree part of the elastic tensor
  call status(0,dtfil%filstat,iexit,level,'call eltfrhar3')
  call eltfrhar3(eltfrhar,rprimd,gsqcut,mpi_enreg,nfftf,ngfftf,rhog)

! Calculate the xc part of the elastic tensor
  call status(0,dtfil%filstat,iexit,level,'call eltfrxc3 ')
  call eltfrxc3(eltfrxc,enxc,gsqcut,kxc,mpi_enreg,dtset%natom,&
& nfftf,ngfftf,nkxc,dtset%nspden,ntypat,psps%n1xccc,n3xccc,dtset%paral_kgb,rhor,rprimd,&
& dtset%typat,vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)

! Calculate the local potential part of the elastic tensor
  call status(0,dtfil%filstat,iexit,level,'call eltfrloc3')
  call eltfrloc3(atindx,eltfrloc,gmet,gprimd,gsqcut,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&  dtset%natom,nattyp,nfftf,ngfftf,ntypat,ph1df,psps%qgrid_vl,rhog,&
&  ucvol,psps%vlspl)

! Calculate the Ewald part of the elastic tensor
  call status(0,dtfil%filstat,iexit,level,'call ewald4')
  call ewald4(elteew,gmet,gprimd,dtset%natom,ntypat,rmet,rprimd,&
&   dtset%typat,ucvol,xred,psps%ziontypat)

! Calculate the psp core energy part of elastic tensor (trivial)
  eltcore(1:3,1:3)=ecore/ucvol

 end if !rfstrs/=0
!End section for strain perturbation

 deallocate(ph1d,ph1df,vpsp,vxc,xccc3d)

 if(prtvol==-level)then
  write(message,'(a,a)') ch10,&
&   ' respfn : frozen wavef. and Ewald(q=0) part of 2DTE done. '
  call wrtout(06,message,'COLL')
 end if

 call timab(136,2,tsec)

 if(rfuser == 1) then ! here is the trial GIPAW code
  if (psps%usepaw /= 1) then
   write (message,'(4a)')' respfn : ERROR- ',ch10,&
&  ' usepaw /= 1 but GIPAW calculation requires PAW ',ch10
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
  allocate(gipaw_aug(dtset%natom))
  do iatom=1, dtset%natom
   allocate(gipaw_aug(iatom)%dia(pawfgrtab(iatom)%nfgd,pawtab(itypat)%lmn2_size))
   allocate(gipaw_aug(iatom)%para(3,pawfgrtab(iatom)%nfgd,pawtab(itypat)%lmn2_size))
   allocate(gipaw_aug(iatom)%onsiteangmom(2,3,pawtab(itypat)%lmn2_size))
  end do ! end allocation loop over atoms
  allocate(gcart(ngfftf(1),ngfftf(2),ngfftf(3),3))
  call gridgcart(gcart,gprimd,ngfftf) ! obtain G vectors in cartesian coords on grid
  allocate(jvec(3,3,nfftf))
  allocate(cs(3,3,dtset%natom))
  write(6,*)"calling gipaw_aug_fields..."
  call gipaw_aug_fields(gipaw_aug,dtset%natom,nfftf,dtset%ntypat,pawfgrtab,pawrad,pawrhoij,pawtab,psps,dtset%typat)
  write(6,*)"calling gipaw_j_dia_aug..."
  call gipaw_j_dia_aug(gipaw_aug,jvec,dtset%natom,nfftf,dtset%ntypat,pawfgrtab,pawrhoij,pawtab,dtset%typat)
  write(6,*)"calling jvec_to_B..."
  call jvec_to_B(cs,gcart,jvec,dtset%natom,nfftf,ngfftf,dtset%paral_kgb,rprimd,xred)
  do iatom = 1, dtset%natom
   do ii = 1, 3
    write(6,'(3f16.8)')cs(ii,1,iatom),cs(ii,2,iatom),cs(ii,3,iatom)
   end do
  end do
!  write(6,*)"calling simple_dia..."
!  call simple_j_dia(jdia,dtset%natom,nfftf,pawfgrtab)
  write(6,*)"calling gipaw_j_dia_bare..."
  call gipaw_j_dia_bare(jvec,nfftf,ngfftf,nhat,dtset%nspden,rhor,rprimd)
  write(6,*)"calling jvec_to_B..."
  call jvec_to_B(cs,gcart,jvec,dtset%natom,nfftf,ngfftf,dtset%paral_kgb,rprimd,xred)
  do iatom = 1, dtset%natom
   do ii = 1, 3
    write(6,'(3f16.8)')cs(ii,1,iatom),cs(ii,2,iatom),cs(ii,3,iatom)
   end do
 end do
!  deallocate(cs,gcart,jvec)

  call leave_new('COLL') ! terminate immediately

 end if ! end rfuser == 1 GIPAW code

! -----3. Initialisation of 1st response, taking into account the q vector.

 call timab(137,1,tsec)

 write(message, '(a,a,a)' )ch10,&
& ' ==>  initialize data related to q vector <== ',ch10
 call wrtout(6,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 qphon(:)=dtset%qptn(:)
 sumg0=1

!Treat the case of q=0 or q too close to 0
 qzero=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2 < 1.d-14)then
  qphon(:)=zero
  write(message, '(a,a,a)' )&
&  ' respfn : the norm of the phonon wavelength (as input) was small (<1.d-7).',&
&  ch10,'  q has been set exactly to (0 0 0)'
  call wrtout(6,message,'COLL')
  sumg0=0
  qzero=1
 else
  if(rfelfd/=0 .or. rfmgfd/=0 .or. rfstrs/=0)then
!  Temporarily, ...
   write(message, '(a,a,a,a,a,3es16.6,a,a,a,i2,a,i2,a,i2,a,a,a)' )ch10,&
&   ' respfn : ERROR -',ch10,&
&   '  The treatment of non-zero wavevector q is restricted to phonons.',&
&   '  However, the input normalized qpt is',qphon(:),',',ch10,&
&   '  while rfelfd=',rfelfd,', rfmgfd=',rfmgfd,', and rfstrs=',rfstrs,'.',ch10,&
&   '  Action : change qpt, or rfelfd, rfmgfd, or rfstrs in the input file.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  else if(rfasr.eq.2)then
   write(message, '(a,a,a,a)' )ch10,&
&   ' respfn : WARNING -',ch10,&
&   '  rfasr=2 not allowed with q/=0 => rfasr was reset to 0.'
   call wrtout(6,message,'COLL')
   rfasr=0
  end if
 end if

!Examine the symmetries of the q wavevector
 allocate(symq(4,2,dtset%nsym))
 call status(0,dtfil%filstat,iexit,level,'call symq3    ')
 call symq3(dtset%nsym,qphon,symq,symrec,timrev)

!Determine the symmetrical perturbations
 allocate(pertsy(3,mpert))
 call status(0,dtfil%filstat,iexit,level,'call syper3   ')
 call syper3(indsym,mpert,dtset%natom,dtset%nsym,pertsy,rfdir,rfpert,symq,symrec,dtset%symrel)
 write(message, '(a)' ) &
& ' The list of irreducible perturbations for this q vector is:'
 call wrtout(ab_out,message,'COLL')
 ii=1
 do ipert=1,mpert
  do idir=1,3
   if(rfpert(ipert)==1.and.rfdir(idir)==1)then
    if( pertsy(idir,ipert)==1 )then
     write(message, '(i5,a,i2,a,i4)' )&
&     ii,')    idir=',idir,'    ipert=',ipert
     call wrtout(ab_out,message,'COLL')
     ii=ii+1
    end if
   end if
  end do
 end do

!Contribution to the dynamical matrix from ion-ion energy
 if(rfphon==1)then
  call status(0,dtfil%filstat,iexit,level,'call ewald3(2)')
  call ewald3(dyew,gmet,dtset%natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,psps%ziontypat)
  option=0
  call q0dy3(dtset%natom,dyewq0,dyew,option)
 end if

!1-order contribution of the xc core correction to the
!dynamical matrix
 allocate(dyfrx1(2,3,dtset%natom,3,dtset%natom))
 dyfrx1(:,:,:,:,:)=zero
 if(rfphon==1.and.psps%n1xccc/=0)then
  call dyxc13(atindx,dyfrx1,gmet,gprimd,gsqcut,kxc,mgfftf,mpi_enreg,&
&  psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,nkxc,dtset%nspden,&
&  ntypat,psps%n1xccc,dtset%paral_kgb,pawtab,ph1df,psps%qgrid_vl,qphon,&
&  rprimd,timrev,dtset%typat,ucvol,psps%usepaw,psps%xcccrc,psps%xccc1d,xred)
 end if

!Deallocate the arrays that were needed only for the frozen
!wavefunction part
 deallocate(cg,eigen0,kg,npwarr)

!Close the unneeded temporary data files, if any
 if (dtset%mkmem==0) then
  close (unit=dtfil%unkg,status='delete')
  if (psps%useylm==1) close (unit=dtfil%unylm,status='delete')
  if (psps%usepaw==1) close (unit=dtfil%unpaw ,status='delete')
  call WffDelete(wfftgs,ierr)
 end if

#if defined MPI
           deallocate(mpi_enreg%proc_distrb)
#endif

 allocate(blkflg(3,mpert,3,mpert))
 allocate(d2eig0(2,3,mpert,3,mpert))
 allocate(d2k0(2,3,mpert,3,mpert))
 allocate(d2lo(2,3,mpert,3,mpert))
 allocate(d2loc0(2,3,mpert,3,mpert))
 allocate(d2nfr(2,3,mpert,3,mpert))
 allocate(d2nl(2,3,mpert,3,mpert))
 allocate(d2nl0(2,3,mpert,3,mpert))
 allocate(d2nl1(2,3,mpert,3,mpert))
 allocate(d2vn(2,3,mpert,3,mpert))
 blkflg(:,:,:,:)=0
 d2eig0(:,:,:,:,:)=zero ; d2k0(:,:,:,:,:)=zero
 d2lo(:,:,:,:,:)=zero   ; d2loc0(:,:,:,:,:)=zero
 d2nfr(:,:,:,:,:)=zero  ; d2nl(:,:,:,:,:)=zero
 d2nl0(:,:,:,:,:)=zero  ; d2nl1(:,:,:,:,:)=zero
 d2vn(:,:,:,:,:)=zero

 prtbbb=dtset%prtbbb
 allocate(d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb))
 allocate(d2cart_bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb))
 if(prtbbb==1)then
  d2cart_bbb(:,:,:,:,:,:)=zero
  d2bbb(:,:,:,:,:,:)=zero
 end if

!Check whether exiting was required by the user.
!If found then do not start minimization steps
!At this first call to chkexi, initialize cpus
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1 ; if(dtset%chkexit==0) openexit=0
 call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
 if (iexit==0) then

! #######################################################################

  write(message,'(a,80a)')ch10,('=',mu=1,80)
  call wrtout(ab_out,message,'COLL')
  call wrtout(06,message,'COLL')

  call status(0,dtfil%filstat,iexit,level,'call loper3   ')

  ddkfil(:)=0
! Note that kg, cg, eigen0, mpw and npwarr are NOT passed to loper3 :
! they will be reinitialized for each perturbation, with an eventually
! reduced set of k point, thanks to the use of symmetry operations.
  call loper3(amass,atindx,atindx1,blkflg,codvsn,cpui,cpus,dimcprj,doccde,&
&  ddkfil,dtfil,dtset,dyew,dyfrlo,dyfrnl,dyfrx1,dyfrx2,d2bbb,d2lo,d2nl,&
&  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&
&  etotal,fermie,gsqcut_eff,iexit,indsym,kxc,&
&  dtset%mkmem,mkqmem,mk1mem,mpert,mpi_enreg,mpsang,nattyp,&
&  nfftf,nhat,dtset%nkpt,nkxc,dtset%nspden,nspinor,dtset%nsym,occ,&
&  paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&  pertsy,prtbbb,psps,rfpert,rhog,rhor,symq,symrec,timrev,tmpfil,&
&  usecprj,vtrial,vxcavg,walli,xred)

! #####################################################################

!End of the check of hasty exit
 end if

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
&  ' ---- first-order wavefunction calculations are completed ----',&
&  ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(06,  message,'COLL')

 deallocate(doccde)


!BEGIN TF_CHANGES
#if defined MPI
! in paral_respfn-case the masters have to reconstruct some array's
 if((mpi_enreg%paral_compil_respfn == 1 .and. mpi_enreg%respfn_master_comm /= MPI_COMM_NULL)) then
  !gather arrays from all cpu's
! call xsum_master(blkflg,0,mpi_enreg%respfn_master_comm,ierr) ! Does not work on some machines
! call xsum_master(d2lo,0,mpi_enreg%respfn_master_comm,ierr)
! call xsum_master(d2nl,0,mpi_enreg%respfn_master_comm,ierr)
! call xsum_master(vtrial,0,mpi_enreg%respfn_master_comm,ierr)
  call xsum_mpi(blkflg,mpi_enreg%respfn_master_comm,ierr)
  call xsum_mpi(d2lo,mpi_enreg%respfn_master_comm,ierr)
  call xsum_mpi(d2nl,mpi_enreg%respfn_master_comm,ierr)
  call xsum_mpi(vtrial,mpi_enreg%respfn_master_comm,ierr)

end if

!DEBUG
!write(6,*) "blkflg",mpi_enreg%paral_compil_respfn,":", blkflg
!write(6,*) "ddkfil",mpi_enreg%paral_compil_respfn,":", ddkfil
!write(6,*) "d2bbb-array",mpi_enreg%paral_compil_respfn,":", d2bbb
!write(6,*) "d2lo-array",mpi_enreg%paral_compil_respfn,":", d2lo
!write(6,*) "d2nl-array",mpi_enreg%paral_compil_respfn,":", d2nl
!write(6,*) "etotal",mpi_enreg%paral_compil_respfn,":", etotal
!write(6,*) "nspinor",mpi_enreg%paral_compil_respfn,":", nspinor
!write(6,*) "fermie",mpi_enreg%paral_compil_respfn,":", fermie
!write(6,*) "vtrial",mpi_enreg%paral_compil_respfn,":", vtrial
!write(6,*) "xred",mpi_enreg%paral_compil_respfn,":", xred
!ENDDEBUG

#endif
!END TF_CHANGES

!Output of the localization tensor
 if ( rfpert(dtset%natom+1) /= 0 .and. (me == 0) .and. dtset%occopt<=2) then
  call wrtloctens(blkflg,d2bbb,d2nl,dtset%mband,mpert,dtset%natom,dtset%prtbbb,rprimd)
 end if

!The perturbation  dtset%natom+1 was only an auxiliary perturbation,
!needed to construct the electric field response, so its flag
!is now set to 0.
! rfpert(dtset%natom+1)=0

!Were 2DTE computed ?
 if(rfphon==0 .and. (rfelfd==2 .or. rfmgfd==2) .and. rfstrs==0 .and. rfuser==0)then

  write(message,'(a,a)' )ch10,&
&  ' respfn : d/dk was computed, but no 2DTE, so no DDB output.'
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

!If 2DTE were computed, only one processor must output them and compute
!frequencies.
 else if(me==0)then

  write(message,'(a,a)' )ch10,&
&  ' ==> Compute Derivative Database <== '
  call wrtout(6,message,'COLL')
  call wrtout(ab_out,message,'COLL')

! Open the formatted derivative database file, and write the
! preliminary information
  call status(0,dtfil%filstat,iexit,level,'call ioddb8   ')
  vrsddb=010929
  dscrpt=' Note : temporary (transfer) database '
  choice=2
  ddbnm=trim(dtfil%filnam_ds(4))//'_DDB'
! tolwfr must be initialized here, but it is a dummy value
  tolwfr=1.0_dp
!DEBUG
! write(6,*)'respfn:ntypat=',ntypat
!ENDDEBUG
  call ioddb8 (choice,dscrpt,ddbnm,dtset%natom,dtset%mband,&
&  dtset%nkpt,dtset%nsym,ntypat,dtfil%unddb,vrsddb,&
&  dtset%acell_orig,dtset%amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&  dtset%intxc,iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&  dtset%natom,dtset%nband,ngfft,dtset%nkpt,dtset%nspden,nspinor,&
&  dtset%nsppol,dtset%nsym,ntypat,occ,dtset%occopt,&
&  dtset%rprim_orig,dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&  dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&  dtset%typat,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

  nblok=1 ; fullinit=1
  call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&  psps%lmnmax,psps%lnmax,nblok,&
&  ntypat,dtfil%unddb,psps%pspso,psps%usepaw,psps%useylm,vrsddb)

! In the RESPFN code, nstdy3 and stady3 were called here
  d2nfr(:,:,:,:,:)=d2lo(:,:,:,:,:)+d2nl(:,:,:,:,:)

! In case of bbb decomposition
  if(prtbbb==1)then
   allocate(blkflg1(3,mpert,3,mpert))
   allocate(blkflg2(3,mpert,3,mpert))
   blkflg2(:,:,:,:) = blkflg(:,:,:,:)
   do ipert = 1, mpert
    do ipert2 = 1, mpert
     if ((ipert /= dtset%natom + 2).and.(ipert>dtset%natom).and.(ipert2/=dtset%natom+2)) then
      blkflg2(:,ipert2,:,ipert) = 0
     end if
    end do
   end do
   allocate(d2tmp(2,3,mpert,3,mpert))
   do iband = 1,dtset%mband
    d2tmp(:,:,:,:,:)=zero
    blkflg1(:,:,:,:) = blkflg2(:,:,:,:)
    d2tmp(:,:,dtset%natom+2,:,:) = d2bbb(:,:,:,:,iband,iband)
    call d2sym3(blkflg1,d2tmp,indsym,mpert,dtset%natom,dtset%nsym,qphon,symq,&
&           symrec,dtset%symrel,timrev)
    d2bbb(:,:,:,:,iband,iband) = d2tmp(:,:,dtset%natom+2,:,:)
   end do
   deallocate(blkflg1,blkflg2,d2tmp)
  end if

! Complete the d2nfr matrix by symmetrization of the existing elements
  call d2sym3(blkflg,d2nfr,indsym,mpert,dtset%natom,dtset%nsym,qphon,symq,symrec,dtset%symrel,timrev)

! Note that there is a bug in d2sym3 which will set some elements of
!  blkflg to 1 even when no corresponding symmetry-related element
!  has been computed.  This has the effect of producing spurious extra
!  output lines in the 2nd-order matrix listing in the .out file
!  and in the DDB file. The suprious matrix elements are all zero,
!  so this is primarily an annoyance.(DRH)


! Add the frozen-wf (dyfrwf) part to the ewald part (dyew),
! the part 1 of the frozen wf part of the xc core correction
! (dyfrx1) and the non-frozen part (dynfr) to get the second-order
! derivative matrix (d2matr), then
! take account of the non-cartesian coordinates (d2cart).
  allocate(d2cart(2,3,mpert,3,mpert))
  allocate(carflg(3,mpert,3,mpert))
  allocate(d2matr(2,3,mpert,3,mpert))
  outd2=1
  call status(0,dtfil%filstat,iexit,level,'call gath3    ')
  call gath3(dtset%berryopt,blkflg,carflg,&
&  dyew,dyfrwf,dyfrx1,d2bbb,d2cart,d2cart_bbb,d2matr,d2nfr,&
&  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&
&  gprimd,dtset%mband,mpert,dtset%natom,ntypat,outd2,dtset%prtbbb,&
&  rfasr,rfdir,rfpert,rprimd,dtset%typat,ucvol,psps%ziontypat)

! Output of the dynamical matrix
! (Note : remember, previously, the processor me=0 has been selected)
  call status(0,dtfil%filstat,iexit,level,'call dyout3   ')
  call dyout3(dtset%berryopt,blkflg,carflg,dtfil%unddb,ddkfil,dyew,dyfrlo,&
&  dyfrnl,dyfrx1,dyfrx2,d2cart,d2cart_bbb,d2eig0,d2k0,d2lo,d2loc0,d2matr,&
&  d2nl,d2nl0,d2nl1,d2vn,&
&  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&
&  ab_out,dtset%mband,mpert,dtset%natom,ntypat,&
&  outd2,dtset%prtbbb,prtvol,qphon,qzero,dtset%typat,rfdir,rfpert,rfphon,&
&  rfstrs,psps%ziontypat)

  close(dtfil%unddb)

! In case of phonons, diagonalize the dynamical matrix
  if(rfphon==1)then

!  First, suppress the 'wings' elements,
!  for which the diagonal element is not known
   call wings3(carflg,d2cart,mpert,dtset%natom)

!  Check the analyticity of the dynamical matrix
   analyt=0
   if (rfpert(dtset%natom+2)==0 .or. rfpert(dtset%natom+2)==2 .or. sumg0==1 ) analyt=1

!  Diagonalize the analytic part
   allocate(displ(2*3*dtset%natom*3*dtset%natom))
   allocate(eigval(3*dtset%natom),eigvec(2*3*dtset%natom*3*dtset%natom))
   allocate(phfrq(3*dtset%natom))
   qphnrm=one
   call phfrq3(dtset%amu,displ,d2cart,eigval,eigvec,indsym,mpert,&
&   dtset%nsym,dtset%natom,dtset%nsym,ntypat,phfrq,qphnrm,qphon,dtset%rprimd_orig,0,dtset%symrel,dtset%typat,ucvol,xred)

!  Print the phonon frequencies
   call prtph3(displ,0,dtset%enunit,-1,ab_out,dtset%natom,phfrq,qphnrm,qphon)

!  Check the completeness of the dynamical matrix and eventually send a warning
   call chkph3(carflg,0,mpert,dtset%natom)

!  In case of a non-analytical part,
!  get the phonon frequencies for three different directions
!  (in cartesian coordinates)
   if(analyt==0)then
    qphnrm=zero
    do idir=1,3
!    Need to know the corresponding dielectric constant
     if(carflg(idir,dtset%natom+2,idir,dtset%natom+2)==1)then
      qphon(:)=zero ; qphon(idir)=one
!     Get the phonon frequencies
      call phfrq3(dtset%amu,displ,d2cart,eigval,eigvec,indsym,mpert,&
&      dtset%nsym,dtset%natom,dtset%nsym,ntypat,phfrq,qphnrm,qphon,dtset%rprimd_orig,0,dtset%symrel,dtset%typat,ucvol,xred)
!     Print the phonon frequencies
      call prtph3(displ,0,dtset%enunit,-1,ab_out,dtset%natom,phfrq,qphnrm,qphon)
!     Check the completeness of the dynamical matrix
!     and eventually send a warning
      call chkph3(carflg,idir,mpert,dtset%natom)
     end if
    end do
   if (idir < 4) then
    qphon(idir)=zero
   end if
   end if

   deallocate(displ,eigval,eigvec)
   deallocate(phfrq)

! End condition on phonons
  end if

  deallocate(carflg,d2cart,d2matr)

!End of the DDB output part
 end if

!Deallocate arrays
 deallocate(amass,atindx,atindx1,blkflg)
 deallocate(dyew,dyewq0,dyfrlo,dyfrnl,dyfrwf,dyfrx1,dyfrx2)
 deallocate(d2bbb,d2cart_bbb)
 deallocate(d2eig0,d2k0,d2lo,d2loc0,d2nfr,d2nl,d2nl0,d2nl1,d2vn)
 deallocate(eltcore,elteew,eltfrhar,eltfrnl,eltfrloc,eltfrkin,eltfrxc)
 deallocate(grxc,indsym,kxc,nattyp,pertsy)
 deallocate(rfpert,rhog,rhor,symq,symrec,vtrial,ylm,ylmgr)
 deallocate(pawfgr%fintocoa,pawfgr%coatofin)
 if (psps%usepaw==1) then
  if (nhatdim>0) deallocate(nhat)
  if (nhatgrdim>0) deallocate(nhatgr)
  call rhoij_free(pawrhoij)
  deallocate(pawrhoij)
  do iatom=1,dtset%natom
   if (associated(pawfgrtab(iatom)%ifftsph))deallocate(pawfgrtab(iatom)%ifftsph)
   if (associated(pawfgrtab(iatom)%rfgd))   deallocate(pawfgrtab(iatom)%rfgd)
   if (associated(pawfgrtab(iatom)%gylm))   deallocate(pawfgrtab(iatom)%gylm)
   if (associated(pawfgrtab(iatom)%gylmgr)) deallocate(pawfgrtab(iatom)%gylmgr)
   if (associated(pawfgrtab(iatom)%gylmgr)) deallocate(pawfgrtab(iatom)%vlocgr)
   pawfgrtab(iatom)%nfgd=0;pawfgrtab(iatom)%rfgd_allocated=0
   pawfgrtab(iatom)%gylm_allocated=0;pawfgrtab(iatom)%gylmgr_allocated=0
   pawfgrtab(iatom)%gylmgr2_allocated=0;pawfgrtab(iatom)%vlocgr_allocated=0
   if (pawtab(itypat)%usepawu>0) deallocate(paw_ij(iatom)%noccmmp,paw_ij(iatom)%nocctot)
   if (pawtab(itypat)%useexexch>0) deallocate(paw_ij(iatom)%vpawx)
   deallocate(paw_an(iatom)%lmselect)
   if (paw_an(iatom)%has_vxcval==1) deallocate(paw_an(iatom)%vxc1_val,paw_an(iatom)%vxct1_val)
   if (paw_ij(iatom)%has_dijhat>0) deallocate(paw_ij(iatom)%dijhat)
   if (paw_ij(iatom)%has_dijxc>0) deallocate(paw_ij(iatom)%dijxc)
   if (paw_ij(iatom)%has_dijxc_val>0) deallocate(paw_ij(iatom)%dijxc_val)
   if (paw_ij(iatom)%has_dijso>0) deallocate(paw_ij(iatom)%dijso)
   if (paw_ij(iatom)%has_dijU>0) deallocate(paw_ij(iatom)%dijU)
   deallocate(paw_ij(iatom)%dij)
  end do
  deallocate(paw_an,paw_ij,pawfgrtab,dimcprj)
 end if

!Clean the header
 call hdr_clean(hdr)

 write(message, '(a,a)' ) ch10,' respfn : exiting '
 call wrtout(06,message,'COLL')

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(137,2,tsec)
 call timab(132,2,tsec)

end subroutine respfn
!!***
