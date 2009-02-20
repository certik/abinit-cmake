!{\src2tex{textfont=tt}}
!!****f* ABINIT/loper3
!! NAME
!! loper3
!!
!! FUNCTION
!! Loop over perturbations
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, DRH, MB, XW, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  amass(natom)=mass of each atom, in unit of electronic mass (=amu*1822...)
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  codvsn=code version
!!  cpui=initial cpu time
!!  cpus=cpu time limit in seconds
!!  dimcprj(natom)=array of dimensions of arrays cprj, cprjq
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  dyew(2,3,natom,3,natom)=Ewald part of the dynamical matrix
!!  dyfrlo(3,3,natom)=frozen wavefunctions local part of the dynamical matrix
!!  dyfrnl(3,3,natom)=frozen wavefunctions non-local part of the dynamical matrix
!!  dyfrx1(2,3,natom,3,natom)=frozen wf nonlin. xc core corr.(2) part of the dynamical matrix
!!  dyfrx2(3,3,natom)=frozen wf nonlin. xc core corr.(2) part of the dynamical matrix
!!  eltcore(6,6)=core contribution to the elastic tensor
!!  elteew(6+3*natom,6)=Ewald contribution to the elastic tsenor
!!  eltfrhar(6,6)=Hartree contribution to the elastic tensor
!!  eltfrkin(6,6)=kinetic contribution to the elastic tensor
!!  eltfrloc(6+3*natom,6)=local psp contribution to the elastic tensor
!!  eltfrnl(6+3*natom,6)=non-local psp contribution to the elastic tensor
!!  eltfrxc(6+3*natom,6)=exchange-correlation contribution to the elastic tensor
!!  fermie=fermi energy (Hartree)
!!  gsqcut_eff=Fourier cutoff on G^2 for "large sphere" of radius double
!!    that of the basis sphere--appropriate for charge density rho(G),
!!    Hartree potential, and pseudopotentials, corresponding to ecut_eff
!!  iexit=index of "exit" on first line of file (0 if not found)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kxc(nfftf,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mkmem =max. number of k points which can fit in memory (GS data)  ; 0 if use disk
!!  mkqmem=max. number of k+q points which can fit in memory (GS data); 0 if use disk
!!  mk1mem=max. number of k points which can fit in memory (RF data)  ; 0 if use disk
!!  mpert=maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotential
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  nhat(nfftf,nspden*usepaw)= -PAW only- compensation density
!!  nkpt=number of k points
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsym=number of symmetry elements in space group
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh for the GS
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pertsy(3,mpert)=set of perturbations that form a basis for all other perturbations
!!  prtbbb=if 1, bbb decomposition, also dimension d2bbb
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rfpert(mpert)=array defining the type of perturbations that have to be computed
!!                1   ->   element has to be computed explicitely
!!               -1   ->   use symmetry operations to obtain the corresponding element
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  symq(4,2,nsym)=1 if symmetry preserves present qpoint. From symq3
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  tmpfil(15)=names for the temporary files based on dtfil%filnam_ds(5)
!!  usecprj= 1 if cprj, cprjq arrays are stored in memory
!!  vtrial(nfftf,nspden)=GS potential (Hartree)
!!  vxcavg=average of vxc potential
!!  walli=initial wall clock time
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  ddkfil(3)=unit numbers for the three possible ddk files
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some second order derivatives
!!  d2lo(2,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,mpert,3,mpert)=non-local contributions to the 2DTEs
!!  etotal=total energy (sum of 7 contributions) (hartree)
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      appdig,atm2fft3,bstruct_clean,bstruct_init,chkexi,clnmpi_band
!!      cprj_alloc,cprj_free,ctocprj,distrb2,eig2tot,fourdp,getmpw,getnel,getph
!!      handle_ncerr,hdr_clean,hdr_init,hdr_update,initmpi_band,initmpi_respfn
!!      initylmg,insy3,inwffil,ioarr,ioddb8,kpgio,leave_new,leave_test,metric
!!      mkcor3,mkrdim,mkrho3,mpi_comm_free,mpi_group_free,outbsd,outgkk,outwf
!!      prteigrs,prtene3,psddb8,rhoij_alloc,rhoij_copy,rhoij_free,scfcv3,setsym
!!      status,symkpt,timab,transgrid,vloca3,vlocalstr,wffclose,wffdelete
!!      wffopen,wrtout,xcast_mpi,xcomm_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine loper3(amass,atindx,atindx1,blkflg,codvsn,cpui,cpus,dimcprj,doccde,&
&  ddkfil,dtfil,dtset,dyew,dyfrlo,dyfrnl,dyfrx1,dyfrx2,d2bbb,d2lo,d2nl,&
&  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,&
&  etotal,fermie,gsqcut_eff,iexit,indsym,kxc,&
&  mkmem,mkqmem,mk1mem,mpert,mpi_enreg,mpsang,nattyp,&
&  nfftf,nhat,nkpt,nkxc,nspden,nspinor,nsym,occ,&
&  paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&  pertsy,prtbbb,psps,rfpert,rhog,rhor,symq,symrec,timrev,tmpfil,&
&  usecprj,vtrial,vxcavg,walli,xred)

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
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_14occeig
 use interfaces_15common
 use interfaces_16response
 use interfaces_18seqpar, except_this_one => loper3
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 integer, intent(in) :: mk1mem,mkmem,mkqmem,mpert,mpsang,nfftf,nkpt
 integer, intent(in) :: nkxc,nspden,nsym,prtbbb,timrev,usecprj
 integer, intent(inout) :: nspinor
 integer, intent(out) :: iexit
 real(dp), intent(in) :: cpui,cpus,gsqcut_eff,vxcavg,walli
 real(dp), intent(inout) :: fermie
 real(dp), intent(out) :: etotal
 character(len=6), intent(in) :: codvsn
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(inout) :: dtset
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(inout) :: psps
 integer, intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer, intent(in) :: dimcprj(dtset%natom),indsym(4,nsym,dtset%natom)
 integer, intent(in) :: nattyp(dtset%ntypat),pertsy(3,mpert)
 integer, intent(in) :: rfpert(mpert),symq(4,2,nsym),symrec(3,3,nsym)
 integer, intent(out) :: blkflg(3,mpert,3,mpert),ddkfil(3)
 real(dp), intent(in) :: amass(dtset%natom),doccde(dtset%mband*nkpt*dtset%nsppol)
 real(dp), intent(in) :: dyew(2,3,dtset%natom,3,dtset%natom)
 real(dp), intent(in) :: dyfrlo(3,3,dtset%natom),dyfrnl(3,3,dtset%natom)
 real(dp), intent(in) :: dyfrx1(2,3,dtset%natom,3,dtset%natom)
 real(dp), intent(in) :: dyfrx2(3,3,dtset%natom),eltcore(6,6)
 real(dp), intent(in) :: elteew(6+3*dtset%natom,6),eltfrhar(6,6)
 real(dp), intent(in) :: eltfrkin(6,6),eltfrloc(6+3*dtset%natom,6)
 real(dp), intent(in) :: eltfrnl(6+3*dtset%natom,6)
 real(dp), intent(in) :: eltfrxc(6+3*dtset%natom,6),kxc(nfftf,nkxc)
 real(dp), intent(in) :: nhat(nfftf,nspden*psps%usepaw)
 real(dp), intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
 real(dp), intent(in) :: rhog(2,nfftf),rhor(nfftf,nspden)
 real(dp), intent(inout) :: vtrial(nfftf,nspden),xred(3,dtset%natom)
 real(dp), intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
 real(dp), intent(out) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert)
 character(len=fnlen), intent(in) :: tmpfil(15)
 type(paw_an_type),intent(inout) :: paw_an(dtset%natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij(dtset%natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!Define file format for different type of files. Presently,
!only one file format is supported for each type of files, but this might
!change soon ...
!1   for wavefunction file, old format (version prior to 2.0)
!2   for wavefunction file, new format (version 2.0 and after)    (fform)
!(51 or 52   for density rho(r)       (fformr)
!101 or 102 for potential V(r) file. (fformv)
 integer,parameter :: fform=2,fformv=102,level=14,response=1
 integer :: fformr=52,unitout
 integer :: accessfil,ask_accurate,band_index,bantot_rbz,bd2tot_index,bdtot1_index
 integer :: bdtot_index,choice,cplex,ddkcase,dimpaw1,formeig,gscase,iatom,iband,idir
 integer :: ieig2rf,ierr,ifft,ii,ikpq,ikpt,ikpt1,index,initialized,ios,ipert,ipw,ir,ireadwf0
 integer :: iscf_mod,isppol,istr,itime,ixc,localrdwf,master,me,mgfftf,mpw,mpw1,mxfh
 integer :: n3xccc,nband_k,ncpgr,nkpt_eff,nkpt_rbz,npert,nqpt,nsym1,ntypat,nxfh
 integer :: openexit,option,optn,optn2,optorth,optthm,optv,pertcase,prtdos,prteig,prtvol,rdwr,rdwrpaw,spaceComm
 integer :: spaceworld,t_iostat,timrev_pert
 integer :: ncerr,ncid_hdr,vrsddb,fullinit,nblok
 real(dp) :: boxcut,dosdeltae,eberry,ecore,ecut_eff,edocc,eei,eeig0,eew,efrhar,efrkin,efrloc
 real(dp) :: efrnl,efrx1,efrx2,ehart,ehart01,ehart1,eig0,eii,ek,ek0,ek1,ek2,eloc0
 real(dp) :: elpsp1,enl,enl0,enl1,entropy,enxc,exc1,fsum,gsqcut,maxocc,mean,nelectkq
 real(dp) :: residm,rho1_dn,rho1im_dn,rho1re_dn,rhosum,tolmxf,tolwfr,ucvol
 logical :: t_exist
 character(len=fnlen) :: fiden1i,fiwf1i,fiwf1o,fiwfddk,wff2nm
 character(len=fnlen) :: gkkfilnam,bsdfilnam,dscrpt
 character(len=4) :: mode_par_psp
 character(len=6) :: tag
 character(len=500) :: message
 type(bandstructure_type) :: bs_rbz
 type(dens_sym_operator_type) :: densymop_rf
 type(hdr_type) :: hdr,hdr0,hdr2
 type(wffile_type) :: wff1,wffddk,wffgs,wffkq,wffnew,wffnow,wfftgs,wfftkq
 type(wvl_data) :: wvl
 integer :: clflg(3,mpert),ngfftf(18),symrc2(3,3,1)
 integer,allocatable :: indkpt1(:),indkpt1_tmp(:),indsy1(:,:,:),irrzon1(:,:,:)
 integer,allocatable :: istwfk_rbz(:),istwfk_pert(:,:,:),kg(:,:),kg1(:,:),nband_rbz(:),npwar1(:)
 integer,allocatable :: npwarr(:),npwtot(:),npwtot1(:),symaf1(:),symaf1_tmp(:)
 integer,allocatable :: symrc1(:,:,:),symrl1(:,:,:),symrl1_tmp(:,:,:),tkq(:)
 real(dp) :: corstr(6),dummy(1),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),tsec(2)
 real(dp),allocatable :: cg(:,:),cg1(:,:),cg1_pert(:,:,:,:),cgq(:,:),gh1_pert(:,:,:,:,:,:),doccde_rbz(:),docckqde(:)
 real(dp),allocatable :: eigen0(:),eigen0_pert(:),eigen1(:),eigen1_pert(:,:,:),eigen1_diag(:)
 real(dp),allocatable :: eigenq(:),eigenq_pert(:),gh1_rbz(:,:,:,:),kpq(:,:),kpqt(:,:),kpt(:,:)
 real(dp),allocatable :: eig2nkq(:,:,:,:,:,:,:)
 real(dp),allocatable :: kpq_rbz(:,:),kpt_rbz(:,:),occ_pert(:),occ_rbz(:),occkq(:)
 real(dp),allocatable :: ph1d(:,:),ph1df(:,:),phnons1(:,:,:),resid(:),rhog1(:,:)
 real(dp),allocatable :: rhor1(:,:),rho1wfg(:,:),rho1wfr(:,:),tnons1(:,:),tnons1_tmp(:,:),vpsp1(:)
 real(dp),allocatable :: work(:),wtk_folded(:),wtk_rbz(:),xccc3d1(:)
 real(dp),allocatable :: xfhist(:,:,:,:)
 real(dp),allocatable :: ylm(:,:),ylm1(:,:),ylmgr(:,:,:),ylmgr1(:,:,:)
 type(cprj_type),allocatable :: cprj(:,:),cprjq(:,:)
 type(pawrhoij_type),allocatable :: pawrhoij1(:),pawrhoij_read(:)

!BEGIN TF_CHANGES
 integer, allocatable :: pert_tmp(:), pert_calc(:)
 integer :: ipert_cnt, icase, igroup_cnt
 integer :: my_group,ngroup_respfn,old_spaceComm
!END TF_CHANGES

!no_abirules
#if !defined MPI
           integer,parameter :: nkpt_max=50
#endif
#if defined MPI
           integer,parameter :: nkpt_max=-1
#endif

! ***********************************************************************

!DEBUG
!write(6,*)' loper3: enter'
!write(6,*)' loper3 : xred=',xred
!stop
!ENDDEBUG

! write(6,*)' loper3 nkpt',nkpt
 call timab(141,1,tsec)

!BEGIN TF_CHANGES
!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)

#if defined MPI
 if(dtset%paral_rf == 1) then
  old_spaceComm=mpi_enreg%spaceComm
  call initmpi_respfn(mpi_enreg,spaceComm)
 end if
#endif


 master=0
!Init me
 call xme_init(mpi_enreg,me)
!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!END TF_CHANGES

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Array on calculated perturbations for eig2rf
 clflg(:,:) =0
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
&   ' loper3 : enter , debug mode '
  call wrtout(06,message,'COLL')
 end if

 prteig=dtset%prteig

!Option input variables
 iscf_mod=dtset%iscf
 ixc=dtset%ixc
 localrdwf=dtset%localrdwf
 ntypat=psps%ntypat

!Obtain dimensional translations in reciprocal space gprimd,
!metrics and unit cell volume, from rprimd.
!Also output rprimd, gprimd and ucvol
 call mkrdim(dtset%acell_orig,dtset%rprim_orig,rprimd)
 call metric(gmet,gprimd,6,rmet,rprimd,ucvol)

!Get FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of respfn.F90
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
  mgfftf=pawfgr%mgfft;ngfftf(:)=pawfgr%ngfft(:)
 else
  mgfftf=dtset%mgfft;ngfftf(:)=dtset%ngfft(:)
 end if

 wff2nm=trim(dtfil%filnam_ds(4))//'_1WF'

 ecut_eff=dtset%ecut*(dtset%dilatmx)**2

!Initialization
 initialized=0
!Are these needed ?
 ecore=zero ; ek=zero ; ehart=zero ; enxc=zero ; eei=zero ; enl=zero ; eii=zero

!Some PAW inits:
 ncpgr=0
 if (psps%usepaw==1) then
  if (dtset%rfphon==1) ncpgr=1
 end if

!Generate the 1-dimensional phases
 allocate(ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom),ph1df(2,3*(2*mgfftf+1)*dtset%natom))
 call getph(atindx,dtset%natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),ph1d,xred)
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
  call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 else
  ph1df(:,:)=ph1d(:,:)
 end if

!Some allocations
 cplex=2-timrev
!HERE (DEBUG) : impose cplex=2
!cplex=2

 allocate(vpsp1(cplex*nfftf))

 call status(0,dtfil%filstat,iexit,level,'before loop   ')

!BEGIN TF_CHANGES

!determine existence of pertubations and of pertubation symmetries
!create array with pertubations which have to be calculated
! write(6,*) "mpert=", mpert
! write(6,*) "Array pertsy"
! write(6,*) pertsy
 allocate(pert_tmp(3*mpert))
 ipert_cnt=0
 do ipert=1,mpert
  do idir=1,3
   if( rfpert(ipert)==1 .and. dtset%rfdir(idir) == 1 )then
!    write(6,*) "checksym"
!    write(6,*) "pertsy(idir,ipert) = ",pertsy(idir,ipert)
!    write(6,*) "dtset%prepanl =",dtset%prepanl
!    write(6,*) "dtset%prepgkk =",dtset%prepgkk
!    write(6,*) "ipert =",ipert
    if ((pertsy(idir,ipert)==1).or.&
&    ((dtset%prepanl == 1).and.(ipert == dtset%natom+2.or.ipert==dtset%natom+5)).or.&
&    ((dtset%prepgkk == 1).and.(ipert <= dtset%natom))  ) then
      ipert_cnt = ipert_cnt+1;
      pert_tmp(ipert_cnt) = idir+(ipert-1)*3;
    else
     write(message, '(a,a,i4,a,i4,a,a,a,a,a,a)' )ch10,&
&     ' The perturbation idir=',idir,'  ipert=',ipert,' is',ch10,&
&     ' symmetric of a previously calculated perturbation.',ch10,&
&     ' So, its SCF calculation is not needed.',ch10
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')
    end if ! Test of existence of symmetry of perturbation
   end if ! Test of existence of perturbation
  end do
 end do
 allocate(pert_calc(ipert_cnt))
 do icase=1,ipert_cnt
  pert_calc(icase) = pert_tmp(icase)
 end do
 deallocate(pert_tmp)

!init groups for parallelization over perturbations if paral_rf=1
#if defined MPI
  if(mpi_enreg%paral_compil_respfn==1) then
   my_group=mpi_enreg%my_respfn_group
   ngroup_respfn=dtset%ngroup_rf
  end if
#endif
!END TF_CHANGES

!Variational part : loop on polarisations
!BEGIN TF_CHANGES
 do icase=1,ipert_cnt
     !calculate only private part of perturbations
     !when current perturbation is not part of my group then skip this loop
#if defined MPI
     if(dtset%paral_rf==1) then
       if(mpi_enreg%respfn_group(modulo(icase,ngroup_respfn)+1) /= my_group) then
         cycle
       end if
     end if
#endif

    if (pert_calc(icase) <= dtset%natom*3) then
     idir = mod(pert_calc(icase),3)
     if (idir==0) idir=3
     ipert=( (pert_calc(icase)-idir) / 3 + 1)
    else
     ipert = dtset%natom + ((pert_calc(icase) - 3*dtset%natom - 1) / 3) + 1
     idir = mod(pert_calc(icase),3)
     if (idir==0) idir=3
    end if
!END TF_CHANGES

     pertcase=idir+(ipert-1)*3
     call status(pertcase,dtfil%filstat,iexit,level,'enter loop    ')

     write(message, '(a,80a,a,a,3f10.6)' ) ch10,('-',ii=1,80),ch10,&
&     ' Perturbation wavevector (in red.coord.) ',dtset%qptn(:)
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')

!    Describe the perturbation :
     if(ipert>=1 .and. ipert<=dtset%natom)then
      write(message, '(a,i4,a,i4)' )&
&      ' Perturbation : displacement of atom',ipert,'   along direction',idir
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
      if(iscf_mod == -3)then
       write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&       ' loper3 : COMMENT -',ch10,&
&       '  The first-order density is imposed to be zero (iscf=-3).',ch10,&
&       '  Although this is strange in the case of phonons,',ch10,&
&       '  you are allowed to do so.'
       call wrtout(6,message,'COLL')
       call wrtout(ab_out,message,'COLL')
      end if
     else if(ipert==dtset%natom+1)then
      write(message, '(a,i4)' )&
&      ' Perturbation : derivative vs k along direction',idir
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
      if( iscf_mod /= -3 )then
       write(message, '(4a)' )ch10,&
&       ' loper3 : COMMENT -',ch10,&
&       '  In a d/dk calculation, iscf is set to -3 automatically.'
       call wrtout(6,message,'COLL')
       call wrtout(ab_out,message,'COLL')
       iscf_mod=-3
      end if
      if( abs(dtset%sciss) > 1.0d-8 )then
       write(message, '(a,a,a,a,f14.8,a,a)' )ch10,&
&       ' loper3 : WARNING -',ch10,&
&       '  Value of sciss=',dtset%sciss,ch10,&
&       '  Scissor with d/dk calculation : you are using a "naive" approach !'
       call wrtout(6,message,'COLL')
       call wrtout(ab_out,message,'COLL')
      end if
     else if(ipert==dtset%natom+2)then
      write(message, '(a,i4)' )&
&      ' Perturbation : homogeneous electric field along direction',idir
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
      if( iscf_mod == -3 )then
       write(message, '(a,a,a,a,a,a)' )ch10,&
&       ' loper3 : COMMENT -',ch10,&
&       '  The first-order density is imposed to be zero (iscf=-3).',ch10,&
&       '  This corresponds to a calculation without local fields.'
       call wrtout(6,message,'COLL')
       call wrtout(ab_out,message,'COLL')
      end if
     else if(ipert==dtset%natom+5)then
      write(message, '(a,i4)' )&
&      ' Perturbation : homogeneous magnetic field along direction, presently set to electric field for testing',idir
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
      if( iscf_mod == -3 )then
       write(message, '(a,a,a,a,a,a)' )ch10,&
&       ' loper3 : COMMENT -',ch10,&
&       '  The first-order density is imposed to be zero (iscf=-3).',ch10,&
&       '  This corresponds to a calculation without local fields.'
       call wrtout(6,message,'COLL')
       call wrtout(ab_out,message,'COLL')
      end if
     else if(ipert>dtset%natom+7 .or. ipert<=0 )then
      write(message, '(a,a,a,a,i4,a,a,a)' )ch10,&
&      ' loper3 : BUG -',ch10,&
&      '  ipert=',ipert,' is outside the [1,dtset%natom+2] interval.',ch10,&
&      '  This perturbation is not (yet) allowed.'
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
      call leave_new('COLL')
     end if
   write(6,*)' loper3 2 nkpt',nkpt
!    Initialize the diverse parts of energy :
     eew=zero ; efrloc=zero ; efrnl=zero ; efrx1=zero ; efrx2=zero
     efrhar=zero ; efrkin=zero
     if(ipert<=dtset%natom)then
      eew=dyew(1,idir,ipert,idir,ipert)
      efrloc=dyfrlo(idir,idir,ipert)
      efrnl=dyfrnl(idir,idir,ipert)
      efrx1=dyfrx1(1,idir,ipert,idir,ipert)
      efrx2=dyfrx2(idir,idir,ipert)
     else if(ipert==dtset%natom+3 .or. ipert==dtset%natom+4) then
!     istr = 1,2,...,6 and indicates the cartesian strain component
      if(ipert==dtset%natom+3) then
       istr=idir
      else
       istr=idir+3
      end if
      eii=eltcore(istr,istr)
      eew=elteew(istr,istr)
      efrhar=eltfrhar(istr,istr)
      efrkin=eltfrkin(istr,istr)
      efrloc=eltfrloc(istr,istr)
      efrnl=eltfrnl(istr,istr)
      efrx1=eltfrxc(istr,istr)
     end if
!    write(6,*)' loper3 2 nkpt',nkpt

!    Determine the subset of symmetry operations (nsym1 operations)
!    that leaves the perturbation invariant,
!    and initialize corresponding arrays symaf1, symrl1 and tnons1.
     allocate(symaf1_tmp(nsym),symrl1_tmp(3,3,nsym),tnons1_tmp(3,nsym))
!   MJV TODO: check whether prepgkk should be used here
     if (dtset%prepanl /= 1 .and. dtset%berryopt /=4 ) then
      call status(pertcase,dtfil%filstat,iexit,level,'call insy3    ')
      call insy3(gprimd,idir,indsym,ipert,dtset%natom,nsym,nsym1,2,&
&      dtset%symafm,symaf1_tmp,symq,symrec,&
&      dtset%symrel,symrl1_tmp,0,dtset%tnons,tnons1_tmp)
     else
      nsym1 = 1
      symaf1_tmp(1) = 1
      symrl1_tmp(:,:,1) = dtset%symrel(:,:,1)
      tnons1_tmp(:,1) = 0_dp
     end if
     allocate(indsy1(4,nsym1,dtset%natom),symrc1(3,3,nsym1))
     allocate(symaf1(nsym1),symrl1(3,3,nsym1),tnons1(3,nsym1))
     symaf1(1:nsym1)=symaf1_tmp(1:nsym1)
     symrl1(:,:,1:nsym1)=symrl1_tmp(:,:,1:nsym1)
     tnons1(:,1:nsym1)=tnons1_tmp(:,1:nsym1)
     deallocate(symaf1_tmp,symrl1_tmp,tnons1_tmp)
!     write(6,*)' loper3 3 nkpt',nkpt

!    Set up corresponding symmetry data
     call status(pertcase,dtfil%filstat,iexit,level,'call setsym   ')
     allocate(irrzon1(dtset%nfft**(1-1/nsym1),2,nspden/dtset%nsppol))
     allocate(phnons1(2,dtset%nfft**(1-1/nsym1),nspden/dtset%nsppol))
     call setsym(densymop_rf,indsy1,irrzon1,1,dtset%natom,&
&     dtset%nfft,dtset%ngfft,nspden,dtset%nsppol,nsym1,&
&     phnons1,symaf1,symrc1,symrl1,tnons1,dtset%typat,xred)

!    Initialize k+q array
     allocate(kpq(3,nkpt))
!    write(6,*)' loper3 4 nkpt',nkpt
     if (ipert==dtset%natom+3.or.ipert==dtset%natom+4) then
      kpq(:,1:nkpt)=dtset%kptns(:,1:nkpt) ! Do not modify, needed for gfortran
     else
      do ikpt=1,nkpt
       kpq(:,ikpt)=dtset%qptn(:)+dtset%kptns(:,ikpt)
      end do
     end if

!    Determine the subset of k-points needed in the "reduced Brillouin zone",
!    and initialize other quantities
!    write(6,*)'loper3 1 nkpt',nkpt
     allocate(indkpt1_tmp(nkpt),wtk_folded(nkpt))
     indkpt1_tmp(:)=0
     optthm=0 ; option=1 ; timrev_pert=timrev

     ieig2rf = dtset%ieig2rf
     if(ieig2rf==1) then
      symrc2(:,:,:) = 0
      do ii=1,3
       symrc2(ii,ii,1)=1
      end do
      call status(pertcase,dtfil%filstat,iexit,level,'call symkpt   ')
      call symkpt(gmet,indkpt1_tmp,dtset%kptns,nkpt,nkpt_rbz,&
       1,option,symrc1,0,dtset%wtk,wtk_folded)
      allocate(kpt(3,nkpt),kpqt(3,nkpt))
      if(.not.allocated(tkq)) allocate(tkq(nkpt))
      kpt(:,:)=dtset%kptns(:,:)
      kpqt(:,:)=kpq(:,:)
      tkq(:) = 0
      do ikpt=1,nkpt
       do ikpq=1,3
        if(kpqt(ikpq,ikpt)>half) kpqt(ikpq,ikpt)= kpqt(ikpq,ikpt) - one
	if(kpqt(ikpq,ikpt).le.-half) kpqt(ikpq,ikpt)= kpqt(ikpq,ikpt) + one
       end do
      end do
      do ikpt=1,nkpt
       do ikpq=1,nkpt
        if(abs(kpt(1,ikpt)-kpqt(1,ikpq))<0.0001_dp) then
	 if(abs(kpt(2,ikpt)-kpqt(2,ikpq))<0.0001_dp) then
	  if(abs(kpt(3,ikpt)-kpqt(3,ikpq))<0.0001_dp) then
           tkq(ikpq) = ikpt
	  end if
	 end if
        end if
       end do
      end do
      write(*,*)'loper3 : tkq',tkq
      deallocate(kpt,kpqt)
     else
!    For the time being, the time reversal symmetry is not used
!    for ddk, elfd, mgfd perturbations.
      if(ipert==dtset%natom+1 .or. ipert==dtset%natom+2 .or. &
&        ipert==dtset%natom+5 .or. dtset%berryopt==4 )timrev_pert=0
      call status(pertcase,dtfil%filstat,iexit,level,'call symkpt   ')
      call symkpt(gmet,indkpt1_tmp,dtset%kptns,nkpt,nkpt_rbz,&
       nsym1,option,symrc1,timrev_pert,dtset%wtk,wtk_folded)
     end if

     allocate(doccde_rbz(dtset%mband*nkpt_rbz*dtset%nsppol))
     allocate(indkpt1(nkpt_rbz),istwfk_rbz(nkpt_rbz))
     allocate(kpq_rbz(3,nkpt_rbz),kpt_rbz(3,nkpt_rbz))
     allocate(nband_rbz(nkpt_rbz*dtset%nsppol),occ_rbz(dtset%mband*nkpt_rbz*dtset%nsppol))
     allocate(wtk_rbz(nkpt_rbz))
     indkpt1(:)=indkpt1_tmp(1:nkpt_rbz)
     do ikpt=1,nkpt_rbz
      istwfk_rbz(ikpt)=dtset%istwfk(indkpt1(ikpt))
      kpq_rbz(:,ikpt)=kpq(:,indkpt1(ikpt))
      kpt_rbz(:,ikpt)=dtset%kptns(:,indkpt1(ikpt))
      wtk_rbz(ikpt)=wtk_folded(indkpt1(ikpt))
     end do

!    Transfer occ to occ_rbz and doccde to doccde_rbz :
!    this is a more delicate issue
!    NOTE : this takes into account that indkpt1 is ordered
     bdtot_index=0
     bdtot1_index=0
     do isppol=1,dtset%nsppol
      ikpt1=1
      do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
!      Must test against ikpt1/=nkpt_rbz+1, before evaluate indkpt1(ikpt1)
       if(ikpt1/=nkpt_rbz+1)then
        if(ikpt==indkpt1(ikpt1))then
         nband_rbz(ikpt1+(isppol-1)*nkpt_rbz)=nband_k
         occ_rbz(1+bdtot1_index:nband_k+bdtot1_index)=&
&                               occ(1+bdtot_index:nband_k+bdtot_index)
         doccde_rbz(1+bdtot1_index:nband_k+bdtot1_index)=&
&                               doccde(1+bdtot_index:nband_k+bdtot_index)
         ikpt1=ikpt1+1
         bdtot1_index=bdtot1_index+nband_k
        end if
       end if
       bdtot_index=bdtot_index+nband_k
      end do
     end do

     deallocate(indkpt1_tmp,wtk_folded)

     call timab(142,1,tsec)
!    Compute maximum number of planewaves at k
     call status(0,dtfil%filstat,iexit,level,'call getmpw-k ')
     call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kpt_rbz,&
&     mpi_enreg,mpw,nkpt_rbz,ucvol)
     call timab(142,2,tsec)

!    Allocate some arrays
     allocate( kg(3,mpw*mkmem) )
     allocate( npwarr(nkpt_rbz),npwtot(nkpt_rbz) )

#if defined MPI
           allocate(mpi_enreg%proc_distrb(nkpt_rbz,dtset%mband,dtset%nsppol))
           mpi_enreg%parareel=0
           call distrb2(dtset%mband,nband_rbz,nkpt_rbz,dtset%nsppol,mpi_enreg)
#endif
     call initmpi_band(mpi_enreg,nband_rbz,nkpt_rbz,dtset%nsppol)
     call status(0,dtfil%filstat,iexit,level,'call kpgio-k  ')

!    Set up the basis sphere of planewaves at k
     call timab(143,1,tsec)
     call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kg,tmpfil(3),&
&     kpt_rbz,mkmem,nband_rbz,nkpt_rbz,'PERS',mpi_enreg,mpw,npwarr,npwtot,dtset%nsppol,dtfil%unkg)
     call timab(143,2,tsec)

     allocate(ylm(mpw*mkmem,mpsang*mpsang*psps%useylm))
     if (ipert==dtset%natom+1 .or. &
&        ipert==dtset%natom+3 .or. ipert==dtset%natom+4) then
      allocate(ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm))
     else
      allocate(ylmgr(1,1,1))
     end if

!    Set up the Ylm for each k point
     if (psps%useylm==1) then
      if(mkmem==0) open(dtfil%unylm,file=tmpfil(10),form='unformatted',status='unknown')
      call status(0,dtfil%filstat,iexit,level,'call initylmg ')
      option=0
      if (ipert==dtset%natom+1 .or. &
&         ipert==dtset%natom+3 .or. ipert==dtset%natom+4) option=1
      call initylmg(gprimd,kg,kpt_rbz,mkmem,mpi_enreg,mpsang,mpw,nband_rbz,nkpt_rbz,&
&                   npwarr,dtset%nsppol,option,rprimd,dtfil%unkg,dtfil%unylm,ylm,ylmgr)
     end if

     ieig2rf = dtset%ieig2rf
     if(ieig2rf==1) then
      if(.not.allocated(istwfk_pert))then
       allocate(istwfk_pert(nkpt,3,mpert),occ_pert(dtset%mband*nkpt*dtset%nsppol))
      end if
      istwfk_pert(:,:,:)=zero
      istwfk_pert(:,idir,ipert)=istwfk_rbz(:)
      occ_pert(:)= occ_rbz(:)
     end if

     bantot_rbz=sum(nband_rbz(1:nkpt_rbz*dtset%nsppol))

     write(message, '(a,a,a)' )ch10,&
&     '--------------------------------------------------------------------------------',&
&     ch10
      call wrtout(ab_out,message,'COLL')
!    Initialize band structure datatype
     allocate(eigen0(bantot_rbz))
     eigen0(:)=zero
     call bstruct_init(bantot_rbz,bs_rbz,doccde_rbz,eigen0,istwfk_rbz,kpt_rbz,&
&     nband_rbz,nkpt_rbz,npwarr,dtset%nsppol,occ_rbz,wtk_rbz)
     deallocate(eigen0)

!    Initialize header
     gscase=0 ! A GS WF file is read
     call hdr_init(bs_rbz,codvsn,dtset,hdr0,pawtab,gscase,psps)

!    Update header, with evolving variables, when available
     call hdr_update(bantot_rbz,etotal,fermie,hdr0,dtset%natom,&
&     residm,rprimd,occ_rbz,pawrhoij,psps%usepaw,xred)

!    Clean band structure datatype (should use it more in the future !)
     call bstruct_clean(bs_rbz)

     call status(0,dtfil%filstat,iexit,level,'call inwffil-k')

!    Initialize wavefunction files and wavefunctions.
     ireadwf0=1 ; formeig=0 ; ask_accurate=1
     optorth=1;if (psps%usepaw==1) optorth=0
     allocate( cg(2,mpw*nspinor*dtset%mband*mkmem*dtset%nsppol) )
     allocate( eigen0(dtset%mband*nkpt_rbz*dtset%nsppol) )
     call timab(144,1,tsec)

     call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
&     formeig,gmet,hdr0,ireadwf0,istwfk_rbz,kg,&
&     kpt_rbz,localrdwf,dtset%mband,&
&     mkmem,mpi_enreg,mpw,nband_rbz,dtset%ngfft,nkpt_rbz,npwarr,nspden,&
&     nspinor,dtset%nsppol,nsym,occ_rbz,optorth,psps,prtvol,rprimd,dtset%symafm,&
&     dtset%symrel,dtset%tnons,dtfil%unkg,wffgs,wfftgs,&
&     dtfil%unwffgs,dtfil%unwftgs,dtfil%fnamewffk,tmpfil(7),wvl)
     call timab(144,2,tsec)

!    PAW statements:
!    Compute projections of GS wavefunctions on non-local projectors (cprj) (and derivatives)
     if (psps%usepaw==1) then
      if (usecprj==1) then
       if(mkmem==0) open(dtfil%unpaw,file=tmpfil(12),form='unformatted',status='unknown')
       call status(0,dtfil%filstat,iexit,level,'call ctocprj')
       allocate(cprj(dimpaw1,nspinor*dtset%mband*mkmem*dtset%nsppol))
       if(mkmem/=0) call cprj_alloc(cprj,ncpgr,dimcprj)
       if (ipert<=dtset%natom) then
        choice=2
        call ctocprj(atindx,cg,choice,cprj,dtfil,gmet,gprimd,ipert,idir,0,istwfk_rbz,&
&                    kg,kpt_rbz,dtset%mband,dtset%mgfft,mkmem,mpi_enreg,mpsang,mpw,&
&                    dtset%natom,nattyp,nband_rbz,1,dtset%ngfft,nkpt_rbz,dtset%nloalg,&
&                    npwarr,nspinor,dtset%nsppol,ntypat,ph1d,psps,rmet,dtset%typat,&
&                    ucvol,dtfil%unpaw,wfftgs,xred,ylm)
       end if
      end if
     end if

!    Close wffgs%unwff, if it was ever opened (in inwffil)
     if (ireadwf0==1) then
      call WffClose(wffgs,ierr)
     end if


!    Compute maximum number of planewaves at kpq
!    Will be useful for both GS wfs at k+q and RF wavefunctions
     call timab(143,1,tsec)
     call status(0,dtfil%filstat,iexit,level,'call getmpw-kq')
     call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kpq_rbz,&
&     mpi_enreg,mpw1,nkpt_rbz,ucvol)
     call timab(143,2,tsec)

!    Allocate some arrays
     allocate( kg1(3,mpw1*mk1mem) )
     allocate( npwar1(nkpt_rbz),npwtot1(nkpt_rbz) )

     call status(0,dtfil%filstat,iexit,level,'call kpgio(2) ')

!    Set up the basis sphere of planewaves at kpq
!    Will be useful for both GS wfs at k+q and RF wavefunctions
     call timab(142,1,tsec)
     call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kg1,tmpfil(5),&
&     kpq_rbz,mk1mem,nband_rbz,nkpt_rbz,'PERS',mpi_enreg,mpw1,&
&     npwar1,npwtot1,dtset%nsppol,dtfil%unkg1)
     call timab(142,2,tsec)
     allocate(ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm))
     if (ipert==dtset%natom+1 .or. &
&        ipert==dtset%natom+3 .or. ipert==dtset%natom+4) then
      allocate(ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm))
     else
      allocate(ylmgr1(1,1,1))
     end if

!    Set up the Ylm for each kpq point
     if (psps%useylm==1) then
      if(mk1mem==0) open(dtfil%unylm1,file=tmpfil(11),form='unformatted',status='unknown')
      call status(0,dtfil%filstat,iexit,level,'call initylmg ')
      option=0
      if (ipert==dtset%natom+1 .or. &
&         ipert==dtset%natom+3 .or. ipert==dtset%natom+4) option=1
      call initylmg(gprimd,kg1,kpq_rbz,mk1mem,mpi_enreg,mpsang,mpw1,nband_rbz,nkpt_rbz,&
&                   npwar1,dtset%nsppol,option,rprimd,dtfil%unkg1,dtfil%unylm1,ylm1,ylmgr1)
     end if

     write(message, '(a,a)' )&
&     '--------------------------------------------------------------------------------',&
&     ch10
     call wrtout(ab_out,message,'COLL')

!    Initialize band structure datatype
     allocate(eigenq(bantot_rbz))
     eigenq(:)=zero
     call bstruct_init(bantot_rbz,bs_rbz,doccde_rbz,eigenq,istwfk_rbz,kpq_rbz,&
&     nband_rbz,nkpt_rbz,npwar1,dtset%nsppol,occ_rbz,wtk_rbz)
     deallocate(eigenq)
!    Initialize header
     call hdr_init(bs_rbz,codvsn,dtset,hdr,pawtab,pertcase,psps)

!    Clean band structure datatype (should use it more in the future !)
     call bstruct_clean(bs_rbz)

     call status(0,dtfil%filstat,iexit,level,'call inwffilkq')

!    Initialize k+q wavefunction files and wavefunctions.
     ireadwf0=1 ; formeig=0 ; ask_accurate=1
     optorth=1;if (psps%usepaw==1) optorth=0
     allocate( cgq(2,mpw1*nspinor*dtset%mband*mkqmem*dtset%nsppol) )
     allocate( eigenq(dtset%mband*nkpt_rbz*dtset%nsppol) )
     call timab(144,1,tsec)
     call inwffil(ask_accurate,cgq,dtset,dtset%ecut,ecut_eff,eigenq,dtset%exchn2n3d,&
&     formeig,gmet,hdr,&
&     ireadwf0,istwfk_rbz,kg1,kpq_rbz,localrdwf,dtset%mband,&
&     mkqmem,mpi_enreg,mpw1,nband_rbz,dtset%ngfft,nkpt_rbz,npwar1,&
&     nspden,nspinor,dtset%nsppol,nsym,occ_rbz,optorth,psps,prtvol,&
&     rprimd,dtset%symafm,dtset%symrel,dtset%tnons,&
&     dtfil%unkg1,wffkq,wfftkq,dtfil%unwffkq,dtfil%unwftkq,dtfil%fnamewffq,tmpfil(8),wvl)
     call timab(144,2,tsec)

!    PAW statements:
!    Compute projections of GS wavefunctions on non-local projectors (cprjq) (and derivatives)
     if (psps%usepaw==1) then
      if (usecprj==1) then
       if(mkqmem==0) open(dtfil%unpawq,file=tmpfil(14),form='unformatted',status='unknown')
       call status(0,dtfil%filstat,iexit,level,'call ctocprjq')
       allocate(cprjq(dtset%natom,nspinor*dtset%mband*mkqmem*dtset%nsppol))
       if(mkqmem/=0) call cprj_alloc(cprjq,0,dimcprj)
       if (ipert<=dtset%natom) then
        choice=1
        call ctocprj(atindx,cgq,choice,cprjq,dtfil,gmet,gprimd,-1,idir,0,istwfk_rbz,&
&                    kg1,kpq_rbz,dtset%mband,dtset%mgfft,mkqmem,mpi_enreg,mpsang,mpw1,&
&                    dtset%natom,nattyp,nband_rbz,1,dtset%ngfft,nkpt_rbz,dtset%nloalg,&
&                    npwar1,nspinor,dtset%nsppol,ntypat,ph1d,psps,rmet,dtset%typat,&
&                    ucvol,dtfil%unpawq,wfftkq,xred,ylm)
       end if
      end if
     end if

     ieig2rf = dtset%ieig2rf
     if(ieig2rf==1) then
      if(.not.allocated(eigen0_pert))then
        allocate(eigen0_pert(dtset%mband*nkpt*dtset%nsppol))
        allocate(eigenq_pert(dtset%mband*nkpt*dtset%nsppol))
      end if
      eigen0_pert(:) = eigen0(:)
      eigenq_pert(:) = eigenq(:)
     end if

!    Close dtfil%unwffkq, if it was ever opened (in inwffil)
     if (ireadwf0==1) then
      call WffClose(wffkq,ierr)
     end if

!    Report on eigenq values
     write(message, '(a,a)' )ch10,' loper3 : eigenq array'
     call wrtout(6,message,'COLL')
     nkpt_eff=nkpt
     if( (prtvol==0.or.prtvol==1.or.prtvol==2) .and. nkpt>nkpt_max ) nkpt_eff=nkpt_max
     band_index=0
     do isppol=1,dtset%nsppol
      do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       if(ikpt<=nkpt_eff)then
        write(message, '(a,i2,a,i5)' )&
&        '  isppol=',isppol,', k point number',ikpt
        call wrtout(6,message,'COLL')
        do iband=1,nband_k,4
         write(message, '(a,4es16.6)')&
&         '  ',eigenq(iband+band_index:min(iband+3,nband_k)+band_index)
         call wrtout(6,message,'COLL')
        end do
       else if(ikpt==nkpt_eff+1)then
        write(message,'(a,a)' )&
&        '  respfn : prtvol=0, 1 or 2, stop printing eigenq.',ch10
        call wrtout(6,message,'COLL')
       end if
       band_index=band_index+nband_k
      end do
     end do

!    Generate occupation numbers at k+q for the reduced BZ
     allocate(docckqde(dtset%mband*nkpt_rbz*dtset%nsppol),occkq(dtset%mband*nkpt_rbz*dtset%nsppol))
     if(0<=dtset%occopt .and. dtset%occopt<=2)then
!     Same occupation numbers at k and k+q (usually, insulating)
      occkq(:)=occ_rbz(:)
!     docckqde is irrelevant in this case
      docckqde(:)=zero
     else
!     Metallic occupation numbers
      option=1
      dosdeltae=zero ! the DOS is not computed with option=1
      maxocc=two/(nspinor*dtset%nsppol)
      call getnel(docckqde,dosdeltae,eigenq,entropy,fermie,maxocc,dtset%mband,&
&      nband_rbz,nelectkq,&
&      nkpt_rbz,dtset%nsppol,occkq,dtset%occopt,option,&
&      dtset%tphysel,dtset%tsmear,tmp_unit,wtk_rbz)
!     Compare nelect at k and nlelect at k+q
      write(message, '(a,a,a,es16.6,a,es16.6,a)')&
&        ' loper3 : total number of electrons, from k and k+q',ch10,&
&        '  fully or partially occupied states are',&
&           dtset%nelect,' and',nelectkq,'.'
      call wrtout(ab_out,message,'COLL')
      call wrtout(6,message,'COLL')
     end if

     if(prtvol==-level)then
      write(message,'(a,a)') ch10,&
&      ' loper3 : initialisation of q part done. '
      call wrtout(06,message,'COLL')
     end if

!    Initialisation of first-order wavefunctions

     write(message, '(a,a,a,i4)' )&
&     ' Initialisation of the first-order wave-functions :',ch10,&
&     '  ireadwf=',dtfil%ireadwf
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     call appdig(pertcase,dtfil%fnamewff1,fiwf1i)
     call appdig(pertcase,wff2nm,fiwf1o)

!    Allocate 1-st order PAW occupancies (rhoij1)
     if (psps%usepaw==1) then
      dimpaw1=dtset%natom;if (ipert<=dtset%natom) dimpaw1=1
      allocate(pawrhoij1(dimpaw1))
      call rhoij_alloc(cplex,dimcprj,nspden,dtset%nsppol,pawrhoij1,dtset%typat)
      if (cplex/=1) then
!      Eventually reallocate hdr%pawrhoij (if cplex=2)
       call rhoij_free(hdr%pawrhoij)
       call rhoij_alloc(cplex,dimcprj,nspden,dtset%nsppol,hdr%pawrhoij,dtset%typat)
      end if
     else
      dimpaw1=0
     end if

!    Initialize wavefunction files and wavefunctions.
     formeig=1 ; ask_accurate=0
     optorth=1;if(psps%usepaw==1.and.dtfil%ireadwf==1)optorth=0
     allocate( cg1(2,mpw1*nspinor*dtset%mband*mk1mem*dtset%nsppol) )
     allocate( eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol) )
     allocate( resid(dtset%mband*nkpt_rbz*dtset%nsppol) )
     call status(pertcase,dtfil%filstat,iexit,level,'call inwffil  ')
     call timab(144,1,tsec)

     call inwffil(ask_accurate,cg1,dtset,dtset%ecut,ecut_eff,eigen1,dtset%exchn2n3d,&
&     formeig,gmet,hdr,&
&     dtfil%ireadwf,istwfk_rbz,kg1,kpq_rbz,localrdwf,&
&     dtset%mband,mk1mem,mpi_enreg,mpw1,nband_rbz,dtset%ngfft,nkpt_rbz,npwar1,&
&     nspden,nspinor,dtset%nsppol,nsym1,occ_rbz,optorth,psps,prtvol,rprimd,&
&     symaf1,symrl1,tnons1,dtfil%unkg1,wff1,wffnow,dtfil%unwff1,dtfil%unwft1,&
&     fiwf1i,tmpfil(1),wvl)
     call timab(144,2,tsec)

     if (psps%usepaw==1.and.dtfil%ireadwf/=0) then
      call rhoij_copy(hdr%pawrhoij(1:dimpaw1),pawrhoij1(1:dimpaw1))
     end if

!DEBUG
!    write(6,*)' loper3 : cg1(:,1)=',cg1(:,1)
!ENDDEBUG

!    Close wff1 (filename fiwf1i), if it was ever opened (in inwffil)
     if (dtfil%ireadwf==1) then
      call WffClose(wff1,ierr)
     end if

!    Initialize PAW temporary file for 1st-order quantities
     if (psps%usepaw==1.and.mk1mem==0) then
      open(dtfil%unpaw1,file=tmpfil(13),form='unformatted',status='unknown')
     end if

!    Initialize second wavefunction file if needed
     if(mkmem==0 .and. dtset%nstep/=0) then
      write(message, '(a,i4,a,a)' )&
&      ' loper3 about to open unit',dtfil%unwft2,' for file=',trim(tmpfil(2))
      call wrtout(06,  message,'PERS')

#if defined HAVE_NETCDF
     if(dtset%accesswff==2) then
    !  Create empty netCDF file
        ncerr = nf90_create(path=trim(tmpfil(2)), cmode=NF90_CLOBBER, ncid=ncid_hdr)
        call handle_ncerr(ncerr," create netcdf wavefunction file")
        ncerr = nf90_close(ncid_hdr)
        call handle_ncerr(ncerr," close netcdf wavefunction file")
     else if(dtset%accesswff==3) then
        write (std_out,*) "FIXME: ETSF I/O support in loper3"
     end if
#endif

      call WffOpen(dtset%accesswff,spaceComm,tmpfil(2),ierr,wffnew,master,me,dtfil%unwft2)
     end if

!    In case of electric field, open the ddk wf file
     if ( ipert==dtset%natom+2 .and. &
&        sum( (dtset%qptn(1:3))**2 ) < 1.0d-7 .and. (dtset%berryopt .ne. 4) ) then
      ddkcase=idir+dtset%natom*3
      call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk)
      write(message, '(a,a)' )&
&      '-loper3 : read the ddk wavefunctions from file: ',fiwfddk
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
      call WffOpen(dtset%accesswff,spaceComm,fiwfddk,ierr,wffddk,master,me,dtfil%unddk)
!DEBUG
!     write(6,*)' loper3 : debug, stop'
!     stop
!ENDDEBUG
     end if
!    In case of magnetic field, open the ddk wf file
     if ( ipert==dtset%natom+5 .and. &
&        sum( (dtset%qptn(1:3))**2 ) < 1.0d-7 .and. (dtset%berryopt .ne. 4) ) then
      ddkcase=idir+dtset%natom*3
      call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk)
      write(message, '(a,a)' )&
&      '-loper3 : read the ddk wavefunctions from file: ',fiwfddk
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
      call WffOpen(dtset%accesswff,spaceComm,fiwfddk,ierr,wffddk,master,me,dtfil%unddk)
!DEBUG
!     write(6,*)' loper3 : debug, stop'
!     stop
!ENDDEBUG
     end if

!    Get first-order local potentials and 1st-order core correction density change
!    (do NOT include xccc3d1 in vpsp1 : this will be done in scfcv3 because vpsp1
!     might become spin-polarized)

!    The correct value of gsqcut should be derived from ecutf and not ecutsm,
!    but the present trick works
     gsqcut=gsqcut_eff

     n3xccc=0;if(psps%n1xccc/=0)n3xccc=cplex*nfftf
     allocate(xccc3d1(n3xccc))

!    PAW: compute Vloc(1) and core(1) together in reciprocal space
!    --------------------------------------------------------------
     if (psps%usepaw==1) then
      optv=1;optn=n3xccc/(cplex*nfftf);optn2=1
      call status(pertcase,dtfil%filstat,iexit,level,'call atm2fft3 ')
      call atm2fft3(atindx,xccc3d1,vpsp1,cplex,dummy,gmet,gprimd,gsqcut,idir,ipert,&
&                   mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,nattyp,3,nfftf,ngfftf,&
&                   ntypat,optn,optn2,optv,dtset%paral_kgb,pawtab,ph1df,psps%qgrid_vl,&
&                   dtset%qptn,dtset%typat,ucvol,psps%usepaw,psps%vlspl,xred)
     else

!    Norm-conserving psp: compute Vloc(1) in reciprocal sp. and core(1) in real sp.
!    ------------------------------------------------------------------------------

      if(ipert==dtset%natom+3 .or. ipert==dtset%natom+4) then
!      Section for strain perturbation
       call status(pertcase,dtfil%filstat,iexit,level,'call vlocalstr')
       call vlocalstr(gmet,gprimd,gsqcut,istr,mgfftf,mpi_enreg,&
&       psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,ntypat,dtset%paral_kgb,ph1df,psps%qgrid_vl,&
&       ucvol,psps%vlspl,vpsp1)
      else
       call status(pertcase,dtfil%filstat,iexit,level,'call vloca3   ')
       call vloca3(atindx,cplex,gmet,gsqcut,idir,ipert,mpi_enreg,psps%mqgrid_vl,dtset%natom,&
&       nattyp,nfftf,ngfftf,ntypat,ngfftf(1),ngfftf(2),ngfftf(3),dtset%paral_kgb,ph1df,psps%qgrid_vl,&
&       dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)
      end if

      if(psps%n1xccc/=0)then
       call status(pertcase,dtfil%filstat,iexit,level,'call mkcor3   ')
       call mkcor3(cplex,idir,ipert,dtset%natom,ntypat,ngfftf(1),psps%n1xccc,&
&       ngfftf(2),ngfftf(3),dtset%qptn,rprimd,dtset%typat,ucvol,&
&       psps%xcccrc,psps%xccc1d,xccc3d1,xred)
!DEBUG
!     write(6,*)' ir    vpsp1     kxc     xccc3d1 '
!     do ir=1,nfftf,13
!      write(6, '(i5,3es14.6)' ) ir,vpsp1(ir),kxc(ir,1),xccc3d1(ir)
!     end do
!     stop
!ENDDEBUG
      end if ! psps%n1xccc/=0
     end if ! usepaw

     eigen1(:)=zero ; resid(:)=zero

!    Get starting charge density and Hartree + xc potential
     allocate(rhor1(cplex*nfftf,nspden),rhog1(2,nfftf))

     if ( (dtfil%ireadwf==0 .and. iscf_mod/=-4 .and. dtset%get1den==0) .or. &
&         (iscf_mod== -3 ) ) then

      rhor1(:,:)=zero ; rhog1(:,:)=zero
!     PAW: rhoij have been set to zero in call to rhoij_alloc above

     else

      if(iscf_mod>0)then

!      cplex=2 gets the complex density, =1 only real part
       call status(pertcase,dtfil%filstat,iexit,level,'call mkrho3   ')
       if (psps%usepaw==1) then
!       Be careful: in PAW, rho does not include the 1-st order compensation
!      density (to be added in scfcv3.F90) !
        allocate(rho1wfg(2,dtset%nfft),rho1wfr(dtset%nfft,nspden))
        call mkrho3(cg,cg1,cplex,densymop_rf,irrzon1,istwfk_rbz,&
&        kg,kg1,dtset%mband,dtset%mgfft,mkmem,mk1mem,mpi_enreg,mpw,mpw1,nband_rbz,&
&        dtset%nfft,dtset%ngfft,nkpt_rbz,npwarr,npwar1,nspden,nspinor,dtset%nsppol,nsym1,&
&        occ_rbz,dtset%paral_kgb,phnons1,rho1wfg,rho1wfr,symaf1,ucvol,dtfil%unkg,&
&        dtfil%unkg1,wffnow,wfftgs,wtk_rbz)
        call transgrid(cplex,mpi_enreg,nspden,+1,1,1,dtset%paral_kgb,pawfgr,rho1wfg,rhog1,rho1wfr,rhor1)
        deallocate(rho1wfg,rho1wfr)
       else
        call mkrho3(cg,cg1,cplex,densymop_rf,irrzon1,istwfk_rbz,&
&        kg,kg1,dtset%mband,dtset%mgfft,mkmem,mk1mem,mpi_enreg,mpw,mpw1,nband_rbz,&
&        dtset%nfft,dtset%ngfft,nkpt_rbz,npwarr,npwar1,nspden,nspinor,dtset%nsppol,nsym1,&
&        occ_rbz,dtset%paral_kgb,phnons1,rhog1,rhor1,symaf1,ucvol,dtfil%unkg,&
&        dtfil%unkg1,wffnow,wfftgs,wtk_rbz)
       end if

      else

       if (me==0) then
        call status(pertcase,dtfil%filstat,iexit,level,'call ioarr    ')
!       Read rho1(r) from a disk file
        rdwr=1;rdwrpaw=psps%usepaw;if(dtfil%ireadwf/=0) rdwrpaw=0
        if (rdwrpaw/=0) then
         allocate(pawrhoij_read(hdr%natom))
         call rhoij_alloc(cplex,dimcprj,nspden,dtset%nsppol,pawrhoij_read,dtset%typat)
        end if
        call appdig(pertcase,dtfil%fildens1in,fiden1i)
        call ioarr(accessfil,rhor1, dtset, etotal,fformr,fiden1i,hdr, mpi_enreg, &
             &     cplex*nfftf,pawrhoij_read,rdwr,rdwrpaw,ngfftf)
        if (rdwrpaw/=0) then
         call rhoij_copy(pawrhoij_read(1:dimpaw1),pawrhoij1(1:dimpaw1))
         call rhoij_free(pawrhoij_read)
        end if
       end if

       call xcast_mpi(rhor1,0,spaceComm,ierr)

!      Compute up+down rho1(G) by fft
       call status(pertcase,dtfil%filstat,iexit,level,'call fourdp   ')
       allocate(work(cplex*nfftf))
       work(:)=rhor1(:,1)
       call fourdp(cplex,rhog1,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       deallocate(work)

      end if

     end if

!    Check whether exiting was required by the user.
!    If found then do not start minimization steps
     openexit=1 ; if(dtset%chkexit==0) openexit=0
     call chkexi(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg,openexit)
!    If immediate exit, and wavefunctions were not read, must zero eigenvalues
     if (iexit/=0)eigen1(:)=zero
     if (iexit==0) then
      call status(pertcase,dtfil%filstat,iexit,level,'call scfcv3   ')

      allocate(gh1_rbz(nkpt_rbz*dtset%ieig2rf,dtset%mband,2,mpw1*nspinor))
      call scfcv3(atindx,atindx1,blkflg,cg,cgq,cg1,cplex,cprj,cprjq,cpus,dimpaw1,&
&      gh1_rbz,densymop_rf,dimcprj,doccde_rbz,docckqde,dtfil,dtset,&
&      d2bbb,d2lo,d2nl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&      ehart,ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek,ek0,ek1,eloc0,elpsp1,&
&      enl,enl0,enl1,enxc,etotal,exc1,fermie,fform,hdr,idir,indkpt1,&
&      indsy1,initialized,ipert,irrzon1,istwfk_rbz,&
&      kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&
&      mpert,mpi_enreg,mpsang,mpw,mpw1,&
&      nattyp,nband_rbz,nfftf,ngfftf,nhat,nkpt,nkpt_rbz,nkxc,&
&      npwarr,npwar1,nspden,nspinor,&
&      nsym1,n3xccc,occkq,occ_rbz,&
&      paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawrhoij1,pawtab,&
&      pertcase,phnons1,ph1d,ph1df,prtbbb,psps,&
&      dtset%qptn,resid,residm,rhog,rhog1,&
&      rhor,rhor1,rprimd,symaf1,symrc1,symrl1,timrev,&
&      tnons1,usecprj,wffddk,wffnew,wffnow,wfftgs,wfftkq,vpsp1,vtrial,&
&      wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1)

      ieig2rf = dtset%ieig2rf
      if(ieig2rf==1) then
       if(.not.allocated(eigen1_pert))then
        allocate(eigen1_pert(2*dtset%mband**2*nkpt*dtset%nsppol,3,mpert))
        allocate(gh1_pert(nkpt,dtset%mband,3,mpert,2,mpw1*nspinor))
        allocate(cg1_pert(2,mpw1*nspinor*dtset%mband*mk1mem*dtset%nsppol,3,mpert))
        eigen1_pert(:,:,:) = zero
        cg1_pert(:,:,:,:) = zero
        gh1_pert(:,:,:,:,:,:) = zero
       end if
       clflg(idir,ipert)=1
       eigen1_pert(1:2*dtset%mband**2*nkpt_rbz*dtset%nsppol,idir,ipert) = eigen1(:)
       cg1_pert(:,:,idir,ipert)=cg1(:,:)
       do ikpt=1,nkpt_rbz
        do iband=1,dtset%mband
         gh1_pert(ikpt,iband,idir,ipert,:,:)=gh1_rbz(ikpt,iband,:,:)
        end do
       end do
      end if
      deallocate(gh1_rbz)
!    End of the check of hasty exit
     end if

     call timab(146,1,tsec)
     write(message, '(80a,a,a,a,a)' ) ('=',ii=1,80),ch10,ch10,&
&     ' ----iterations are completed or convergence reached----',&
&    ch10
     call wrtout(ab_out,message,'COLL')
     call wrtout(06,  message,'COLL')

!    Update the content of the header (evolving variables)
     call hdr_update(bantot_rbz,etotal,fermie,hdr,dtset%natom,&
&     residm,rprimd,occ_rbz,pawrhoij1,psps%usepaw,xred)

     call status(pertcase,dtfil%filstat,iexit,level,'call outwf    ')

     mxfh=0 ; nxfh=0 ; nqpt=1
     allocate(xfhist(3,dtset%natom+4,2,mxfh))
     call outwf(cg1,dtfil,dtset,eigen1,fiwf1o,hdr,kg1,kpt_rbz,&
&     dtset%mband,mk1mem,mpi_enreg,mpw1,mxfh,dtset%natom,nband_rbz,dtset%nfft,&
&     dtset%ngfft,nkpt_rbz,npwar1,nqpt,nspinor,dtset%nsppol,dtset%nstep,&
&     nxfh,occ_rbz,resid,response,wffnow,wvl%wfs,xfhist)
     deallocate(xfhist)
     if(mkmem==0) then
        call WffDelete(wffnow,ierr)
     end if

!    print _gkk file for this perturbation
     if (dtset%prtgkk == 1) then
      !if (dtset%prepgkk == 1) then
       ! write (6,*) 'Error: the combination of prtgkk and prepgkk does not work yet.'
        !call leave_new('COLL')
      !end if
      call appdig(3*(ipert-1)+idir,trim(dtfil%filnam_ds(4))//"_GKK",gkkfilnam)
      write (*,*) 'gkkfilnam = ', trim(gkkfilnam)
      call outgkk(bantot_rbz,2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol,gkkfilnam,&
    &      dtset,eigen0,eigen1,hdr0,hdr)
     end if

!    If the perturbation is d/dk, evaluate the f-sum rule.
     if( ipert==dtset%natom+1 )then
!     Note : the factor of two is related to the difference
!     between Taylor expansion and perturbation expansion
!     Note : this expression should be modified for ecutsm.
!     Indeed, the present one will NOT tend to 1.0_dp.
      ek2=gmet(idir,idir)*(two_pi**2)*2.0_dp*dtset%nelect
      fsum=-ek1/ek2
      if(dtset%ecutsm<tol6)then
       write(message, '(a,es20.10,a,a,es20.10)' ) &
&       ' loper3 : ek2=',ek2,ch10,&
&       '          f-sum rule ratio=',fsum
      else
       write(message, '(a,es20.10,a,a,es20.10,a)' ) &
&       ' loper3 : ek2=',ek2,ch10,&
&       '          f-sum rule ratio=',fsum,' (note : ecutsm/=0)'
      end if
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')


!     Write the diagonal elements of the dH/dk operator
      allocate(eigen1_diag(dtset%mband*nkpt_rbz*dtset%nsppol))
      bdtot_index=0 ; bd2tot_index=0
      do isppol=1,dtset%nsppol
       do ikpt=1,nkpt_rbz
        nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
        do iband=1,nband_k
         eigen1_diag(iband+bdtot_index)=&
&         eigen1(2*iband-1 + (iband-1)*2*nband_k + bd2tot_index)
        end do

!       Treat the case of degeneracies : take the mean of degenerate states
        if(nband_k>1)then
         eig0=eigen0(1+bdtot_index)
         ii=1
         do iband=2,nband_k
          if(eigen0(iband+bdtot_index)-eig0<tol8)then
           ii=ii+1
          else
           mean=sum(eigen1_diag(iband-ii+bdtot_index:iband-1+bdtot_index))/ii
           eigen1_diag(iband-ii+bdtot_index:iband-1+bdtot_index)=mean
           ii=1
          end if
          eig0=eigen0(iband+bdtot_index)
          if(iband==nband_k)then
           mean=sum(eigen1_diag(iband-ii+1+bdtot_index:iband+bdtot_index))/ii
           eigen1_diag(iband-ii+1+bdtot_index:iband+bdtot_index)=mean
          end if
         end do
        end if

        bdtot_index=bdtot_index+nband_k
        bd2tot_index=bd2tot_index+2*nband_k**2
       end do
      end do
      option=4
      call prteigrs(eigen1_diag,dtset%enunit,fermie,tmpfil(1),ab_out,iscf_mod,kpt_rbz,dtset%kptopt,&
&      dtset%mband,nband_rbz,nkpt_rbz,dtset%nnsclo,dtset%nsppol,occ_rbz,dtset%occopt,&
&      option,prteig,prtvol,resid,tolwfr,vxcavg,wtk_rbz)
      deallocate(eigen1_diag)

     end if

!    Print the energies
     if (dtset%nline/=0 .or. dtset%nstep/=0)then
      call status(pertcase,dtfil%filstat,iexit,level,'call prtene3  ')
      call prtene3(dtset%berryopt,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&       ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,exc1,ab_out,&
&       ipert,dtset%natom)
     end if

     deallocate(cg,cgq,cg1,docckqde,doccde_rbz,eigen0,eigenq,eigen1,kpq)
     deallocate(indkpt1,indsy1,istwfk_rbz)
     deallocate(irrzon1,kg,kg1,kpq_rbz,kpt_rbz)
     deallocate(nband_rbz,npwarr,npwar1,npwtot,npwtot1,occkq,occ_rbz,phnons1,resid,rhog1,rhor1)
     deallocate(symaf1,symrc1,symrl1,tnons1,wtk_rbz,xccc3d1)
     deallocate(ylm,ylm1,ylmgr,ylmgr1)
     if (psps%usepaw==1) then
      call rhoij_free(pawrhoij1)
      deallocate(pawrhoij1)
      if (usecprj==1) then
       call cprj_free(cprj)
       call cprj_free(cprjq)
       deallocate(cprj,cprjq)
      end if
     end if

     call clnmpi_band(nkpt_rbz,dtset%nsppol,mpi_enreg)
#if defined MPI
           deallocate(mpi_enreg%proc_distrb)
#endif

!    Clean hdr
     call hdr_clean(hdr)
     call hdr_clean(hdr0)

!    Close the unneeded temporary data files, if any
     if (mkmem==0) then
      close (unit=dtfil%unkg,status='delete')
      if (psps%useylm==1) close (unit=dtfil%unylm,status='delete')
      if (psps%usepaw==1.and.usecprj==1) close (unit=dtfil%unpaw,status='delete')
      call WffDelete(wfftgs,ierr)
     end if
     if (mkqmem==0) then
      if (psps%usepaw==1.and.usecprj==1) close (unit=dtfil%unpawq,status='delete')
      call WffDelete(wfftkq,ierr)
     end if
     if (mk1mem==0) then
      close (unit=dtfil%unkg1,status='delete')
      if (psps%useylm==1) close (unit=dtfil%unylm1,status='delete')
      if (psps%usepaw==1) close (unit=dtfil%unpaw1,status='delete')
      call WffDelete(wffnew,ierr)
     end if

     call timab(146,2,tsec)
!BEGIN TF_CHANGES
  if(iexit/=0)exit
!End loop over icase
 end do

!#################################################################################ICI>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!Calculate the second-order eigenvalues for a wavevector Q
ieig2rf = dtset%ieig2rf
if(ieig2rf==1) then
  if(dtset%kptopt==3 .or. dtset%kptopt==0)then

   npert=dtset%natom+2
   allocate(eig2nkq(2,dtset%mband*dtset%nsppol,nkpt_rbz,3,npert,3,npert))
   write(*,*) 'loper3: clflg', clflg
   write(*,*) 'loper3: tkq', tkq
   call eig2tot(clflg,cg1_pert,gh1_pert,eigen0_pert,eigenq_pert,eigen1_pert,eig2nkq,&
&   indsym,istwfk_pert,dtset%mband,mk1mem,dtset%natom,npert,nsym,mpi_enreg,mpw1,nkpt_rbz,&
&   nspinor,dtset%nsppol,occ_pert,dtset%qptn,dtset%sciss,symq,symrec,dtset%symrel,timrev,tkq)

! print _BSD file for this perturbation
   bsdfilnam = trim(dtfil%filnam_ds(4))//"_BSD"
   write (*,*) 'bsdfilnam = ', trim(bsdfilnam)

   unitout = dtfil%unddb
   vrsddb=010929
   dscrpt=' Note : temporary (transfer) database '
   choice=2
! tolwfr must be initialized here, but it is a dummy value
   tolwfr=1.0_dp
   call ioddb8 (choice,dscrpt,bsdfilnam,dtset%natom,dtset%mband,&
&   dtset%nkpt,dtset%nsym,dtset%ntypat,dtfil%unddb,vrsddb,&
&   dtset%acell_orig,dtset%amu,dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&   dtset%intxc,dtset%iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&   dtset%natom,dtset%nband,dtset%ngfft,dtset%nkpt,dtset%nspden,nspinor,&
&   dtset%nsppol,dtset%nsym,dtset%ntypat,occ_pert,dtset%occopt,&
&   dtset%rprim_orig,dtset%sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&   dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&   dtset%typat,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

   nblok=1 ; fullinit=1
   call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&   psps%lmnmax,psps%lnmax,nblok,&
&   ntypat,dtfil%unddb,psps%pspso,psps%usepaw,psps%useylm,vrsddb)

   call outbsd(bsdfilnam,dtset,eig2nkq,npert,nkpt_rbz,unitout)

   deallocate(eig2nkq)
  else
     write(message, '(a,a,a,a,a)' ) ch10,&
&   ' eig2tot: WARNING -',ch10,&
&   ' K point grids must be the same for every perturbation: eig2tot not called', &
&        ' Action: Put kptopt=3 '
  end if !kptopt
 deallocate(eigen0_pert,eigen1_pert,eigenq_pert,gh1_pert,cg1_pert,istwfk_pert)
 deallocate(occ_pert,tkq)
end if  !if ieig2rf
!######################################################################################>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#if defined MPI
          !test whether _all_ groups are working properly
          !set communicator to old_spaceComm before
           if(mpi_enreg%paral_compil_respfn == 1) then
            mpi_enreg%spaceComm=old_spaceComm
            call leave_test(mpi_enreg)
          ! reset communicator
            mpi_enreg%spaceComm=spaceComm
           end if
#endif

!END TF_CHANGES

!Get ddk file information, for later use in dyout3
 ddkfil(:)=0
 do idir=1,3
  ddkcase=idir+dtset%natom*3
  call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk)
! Check that ddk file exists
  inquire(file=fiwfddk,iostat=t_iostat,exist=t_exist)
! If the file exists set ddkfil to a non-zero value
  if (t_exist) then
   ddkfil(idir)=20+idir
  end if
 end do

 call status(0,dtfil%filstat,iexit,level,'end loop      ')

!Deallocate arrays
 deallocate(ph1d,ph1df,vpsp1)
 deallocate(pert_calc)
#if defined MPI
!free groups
 if(mpi_enreg%paral_compil_respfn == 1) then
  do igroup_cnt=1,dtset%ngroup_rf
        if(mpi_enreg%respfn_comm(igroup_cnt) == mpi_enreg%spaceComm) then

                !DEBUG
                !PRINT *,mpi_enreg%me,":Freeing Comm - comm=",mpi_enreg%respfn_comm(igroup_cnt)

                call MPI_COMM_FREE(mpi_enreg%respfn_comm(igroup_cnt),ierr)
                if(ierr /= MPI_SUCCESS) then
                        write(6,*)mpi_enreg%me,': Error on releasing Communicator for group nr ',igroup_cnt
                        call leave_new('COLL')
                end if
                call MPI_GROUP_FREE(mpi_enreg%respfn_group(igroup_cnt),ierr)
                if(ierr /= MPI_SUCCESS) then
                        write(6,*)mpi_enreg%me,': Error on releasing group nr ',igroup_cnt
                        call leave_new('COLL')
                end if
        end if
  end do

  deallocate(mpi_enreg%respfn_group)
  deallocate(mpi_enreg%respfn_comm)
  mpi_enreg%spaceComm=old_spaceComm
  mpi_enreg%me_respfn=-1
 end if
#endif

 write(message, '(a,a)' ) ch10,' loper3 : exiting '
 call wrtout(06,message,'COLL')

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(141,2,tsec)
!DEBUG
!write(6,*)' loper3: exit '
!ENDDEBUG
end subroutine loper3
!!***
