!{\src2tex{textfont=tt}}
!!****f* ABINIT/ctocprj
!! NAME
!! ctocprj

!! FUNCTION
!!  Compute all <Proj_i|Cnk> for every wave function |Cnk> expressed in reciprocal space.
!!  |Proj_i> are non-local projectors (for each atom and each l,m,n)
!!  Can also compute derivatives of <Proj_i|Cnk> wrt to several parameters
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  choice: chooses derivatives to compute:
!!          =1 => no derivatives
!!          =2 => 1st derivatives with respect to atomic position(s)
!!          =3 => 1st derivatives with respect to strain(s)
!!          =23=> 1st derivatives with respect to atm. pos. and strain(s)
!!          =4 => 2nd derivatives with respect to atomic pos.
!!          =24=> 1st and 2nd derivatives with respect to atomic pos.
!!          =5 => drivatives with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  dtfil <type(datafiles_type)>=variables related to files
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  iatom= if <=0, cprj=<p_i|Cnk> are computed for all atoms 1...natom
!!         if  >0  cprj=<p_i|Cnk> are computed only for atom with index iatom
!!  idir=direction of the derivative, i.e. dir. of - atom to be moved  in the case choice=2
!!                                                 - strain component  in the case choice=3
!!                                                 - k point direction in the case choice=5
!!       Compatible only with choice=2,3,5; if idir=0, all derivatives are computed
!!  iorder_cprj=0 if output cprj=<p_i|Cnk> are sorted by atom type
!!              1 if output cprj=<p_i|Cnk> are unsorted
!!  istwfk(nkpt)=option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  kpt(3,nkpt)=reduced coordinates of k points.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw
!!  natom=number of atoms in cell
!!  nattyp(ntypat)= # atoms of each type
!!  nband(nkpt*nsppol)=number of bands at this k point for that spin polarization
!!  ncprj=1st dim. of cprj array (natom if iatom<=0, 1 if iatom>0)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~ABINIT/Infos/vargs.htm#ngfft
!!  nkpt=number of k points
!!  nloalg(5)=governs the choice of the algorithm for nonlocal operator
!!  npwarr(nkpt)=number of planewaves in basis at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  tim_ctocprj=timing code of the calling routine
!!  typat(natom)= types of atoms
!!  uncp=unit number for <P_lmn|Cnk> data (if used)
!!  wffnow=struct infos for wf disk file
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang)=real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  cprj(ncprj,nspinor*mband*mkmem*nsppol) <type(cprj_type)>=
!!            projected input wave functions <Proj_i|Cnk> with NL projectors
!!
!! PARENTS
!!      loper3,scfcv,vtorho
!!
!! CHILDREN
!!      cprj_alloc,cprj_free,cprj_put,getcprj,hdr_skip,leave_new,mkffnl,mkkpg
!!      ph1d3d,rdnpw,rwwf,status,wrtout,xallgather_mpi,xallgatherv_mpi
!!      xalltoallv_mpi,xcomm_init,xdefineoff
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine ctocprj(atindx,cg,choice,cprj,dtfil,gmet,gprimd,iatom,idir,iorder_cprj,istwfk,&
&                   kg,kpt,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nband,&
&                   ncprj,ngfft,nkpt,nloalg,npwarr,nspinor,nsppol,ntypat,ph1d,psps,&
&                   rmet,typat,ucvol,uncp,wffnow,xred,ylm)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13io_mpi
 use interfaces_13nonlocal, except_this_one => ctocprj
 use interfaces_14iowfdenpot
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: choice,iatom,idir,iorder_cprj,mband,mgfft,mkmem,mpsang,mpw
 integer,intent(in) :: natom,ncprj,nkpt,nsppol,ntypat,uncp
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ucvol
 type(datafiles_type),intent(in) :: dtfil
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: atindx(natom),istwfk(nkpt),nattyp(ntypat),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(5),npwarr(nkpt),kg(3,mpw*mkmem),typat(natom)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kpt(3,nkpt)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),rmet(3,3)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,mpsang*mpsang)
 type(cprj_type),intent(out) :: cprj(natom,nspinor*mband*mkmem*nsppol)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=7
 integer :: counter,cpopt,dimffnl,formeig,ia,iatom1,iatom2,iatm
 integer :: iband_max,iband_min,ibg,iblockbd,icg,icgb,icpgr,ider,idir0,ierr
 integer :: iexit,ig,ii,ikg,ikpt,ilm,ilmn,ipw,isize,isppol,istwf_k,itypat
 integer :: matblk,mcg_disk,me_distrb,muig,n1,n2,n3
 integer :: nband_k,nblockbd,ncpgr,nkpg,nprocband,npw_k,npw_nk,ntypat0
 integer :: shift1,shift2,shift3,spaceComm,spaceComm_band,spaceComm_fft,tim_rwwf
 real(dp) :: arg
 character(len=500) :: message
!arrays
 integer :: nattyp_atm(1),pspso_atm(1)
 integer,allocatable :: bufsize(:),bufsize_wf(:),bufdisp(:),bufdisp_wf(:)
 integer,allocatable :: dimlmn(:),indlmn_atm(:,:,:),kg_dum(:,:),kg_k(:,:),kg_k_loc(:,:)
 integer,allocatable :: npw_block(:),npw_disp(:)
 real(dp) :: kpoint(3),tsec(2),ylmgr_dum(1)
 real(dp),allocatable :: cg_disk(:,:),cwavef(:,:),cwavef_tmp(:,:),eig_dum(:),ekb_atm(:,:,:)
 real(dp),allocatable :: ffnl(:,:,:,:),ffnl_npw(:,:,:,:),ffnl_tmp(:,:,:,:),ffnl_tmp_npw(:,:,:,:)
 real(dp),allocatable :: ffspl_atm(:,:,:,:),kpg_k(:,:),kpg_k_loc(:,:,:),occ_dum(:),ph1d_atm(:,:)
 real(dp),allocatable :: ph3d(:,:,:),ph3d_npw(:,:,:),ph3d_tmp(:,:,:),ph3d_tmp_npw(:,:,:)
 real(dp),allocatable :: phkxred(:,:),ylm_k(:,:)
 type(cprj_type),allocatable :: cwaveprj(:,:)

! *********************************************************************

!Preliminary tests
 if (psps%useylm==0) then
  write(message, '(4a)' ) ch10,&
&  ' ctocprj : ERROR -',ch10,&
&  '  Not available for useylm=0 !'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 if ((choice<1.or.choice>6).and.choice/=23.and.choice/=24) then
  write(message, '(4a)' ) ch10,&
&  ' ctocprj : BUG -',ch10,&
&  '  Bad choice !'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 if (idir>0.and.choice/=2.and.choice/=3.and.choice/=5) then
  write(message, '(a,a,a,a,i4,a)' ) ch10,&
&  ' ctocprj : BUG -',ch10,&
&  '  Does not support idir>0 for choice=',choice,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 if ((iatom<=0.and.ncprj/=natom).or.(iatom>0.and.ncprj/=1)) then
  write(message, '(4a)' ) ch10,&
&  ' ctocprj : BUG -',ch10,&
&  '  Bad value for ncprj dimension !'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 if (mkmem==0.and.mpi_enreg%mode_para=='b') then
  write(message, '(4a)' ) ch10,&
&  ' ctocprj : ERROR -',ch10,&
&  '  Not available for mkmem=0 and band-FFT parallelism !'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call status(0,dtfil%filstat,iexit,level,'enter ctocprj ')

!Init parallelism
 call xcomm_init(mpi_enreg,spaceComm)
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
& (mpi_enreg%paral_compil_fft==1)) then
  me_distrb = mpi_enreg%me_kpt
 else
  me_distrb = mpi_enreg%me
 end if
 if (mpi_enreg%mode_para=='b') then
  spaceComm_band=mpi_enreg%comm_band
  spaceComm_fft=mpi_enreg%comm_fft
  nprocband=mpi_enreg%nproc_band
 else
  spaceComm_band=0;spaceComm_fft=0
  nprocband=1
 end if

!Prepare temporary files if mkmem==0
 if (mkmem==0) then
  formeig=0;mcg_disk=mpw*nspinor*mband
  call hdr_skip(wffnow,ierr)
  call xdefineOff(formeig,wffnow,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)
  allocate(cg_disk(2,mcg_disk))
 end if

!Initialize some variables
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 ibg=0;icg=0;cpopt=0
 ider=0;idir0=0
 if (idir>0) then
  if (choice==3.or.choice==5) ider=1
  if (choice==3) idir0=-idir
  if (choice==5) idir0=4
 end if
 dimffnl=1+ider
 nkpg=0
 if (choice==3.or.choice==2.or.choice==23) nkpg=3*nloalg(5)
 if (choice==4.or.choice==24) nkpg=9*nloalg(5)

!Set dimension of <p_i|Cnk>
 ncpgr=0
 if (idir==0) then
  if (choice==2) ncpgr=3
  if (choice==3) ncpgr=6
  if (choice==23)ncpgr=9
  if (choice==4) ncpgr=6
  if (choice==24)ncpgr=9
  if (choice==5) ncpgr=3
  if (choice==6) ncpgr=63
 else
  ncpgr=1
 end if

!If only one atom is selected, extract data for this atom
 if (iatom>0) then
  iatom1=iatom;iatom2=iatom
  ntypat0=1;itypat=typat(iatom)
  nattyp_atm(1)=1
  pspso_atm(1)=psps%pspso(itypat)
  allocate(ekb_atm(psps%dimekb,1,nspinor**2),&
&  indlmn_atm(6,psps%lmnmax,1),&
&  ffspl_atm(psps%mqgrid_ff,2,psps%lnmax,1),&
&  ph1d_atm(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
  ekb_atm(:,1,1)=psps%ekb(:,itypat)
  indlmn_atm(:,:,1)=psps%indlmn(:,:,itypat)
  ffspl_atm(:,:,:,1)=psps%ffspl(:,:,:,itypat)
  shift1=(atindx(iatom)-1)*(2*n1+1)
  shift2=(atindx(iatom)-1)*(2*n2+1)+natom*(2*n1+1)
  shift3=(atindx(iatom)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
  ph1d_atm(:,1:2*n1+1)=ph1d(:,1+shift1:2*n1+1+shift1)
  ph1d_atm(:,1+2*n1+1:2*n2+1+2*n1+1)=ph1d(:,1+shift2:2*n2+1+shift2)
  ph1d_atm(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=ph1d(:,1+shift3:2*n3+1+shift3)
 else
  iatom1=1;iatom2=natom
  ntypat0=ntypat
 end if

!Dimensioning and allocation of <p_i|Cnk>
 allocate(dimlmn(ncprj));dimlmn=0  ! Type-sorted cprj
 if (ncprj==natom) then
  ia=0
  do itypat=1,ntypat
   dimlmn(ia+1:ia+nattyp(itypat))=count(psps%indlmn(3,:,itypat)>0)
   ia=ia+nattyp(itypat)
  end do
 else
  itypat=typat(iatom)
  dimlmn(ia+1:ia+nattyp(itypat))=count(psps%indlmn(3,:,itypat)>0)
 end if
 if (mkmem==0) then
  rewind(uncp)
  write(uncp) ncprj,nspinor
  write(uncp) dimlmn(1:ncprj)
 end if
 allocate(cwaveprj(ncprj,nspinor))
 call cprj_alloc(cwaveprj,ncpgr,dimlmn)

!Additional statements if band-fft parallelism
 if (nprocband>1) then
  allocate(npw_block(nprocband),npw_disp(nprocband))
  allocate(bufsize(nprocband),bufdisp(nprocband))
  allocate(bufsize_wf(nprocband),bufdisp_wf(nprocband))
 end if

!LOOP OVER SPINS
 do isppol=1,nsppol
  ikg=0

! Rewind temporary files if needed
  if (mkmem==0) rewind(dtfil%unkg)
  if (mkmem==0) rewind(dtfil%unylm)

! BIG FAT k POINT LOOP
  do ikpt=1,nkpt
   counter=100*ikpt+isppol
   call status(counter,dtfil%filstat,iexit,level,'loop ikpt     ')

!  Select k point to be treated by this proc
   nband_k=nband(ikpt+(isppol-1)*nkpt)
   if(mpi_enreg%paral_compil_kpt==1)then
    if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me_distrb))/=0) cycle
   end if

!  Old FFT parallelism: define FFT communicator for this k-point
   if (mpi_enreg%paral_compil_fft==1.and.mpi_enreg%mode_para/='b') then
    mpi_enreg%num_group_fft=ikpt+(isppol-1)*nkpt
   end if

!  Retrieve k-point
   kpoint(:)=kpt(:,ikpt)
   istwf_k=istwfk(ikpt)

!  Retrieve number of plane waves
   npw_k=npwarr(ikpt)
   if (nprocband>1) then
!   Special treatment for band-fft //
    call xallgather_mpi(npw_k,npw_block,spaceComm_band,ierr)
    npw_nk=sum(npw_block);npw_disp(1)=0
    do ii=2,nprocband
     npw_disp(ii)=npw_disp(ii-1)+npw_block(ii-1)
    end do
   else
    npw_nk=npw_k
   end if

!  Test cprj gradients dimension (just to be sure)
   if (cprj(1,ibg+1)%ncpgr/=ncpgr) then
    write(message, '(4a)' ) ch10,&
&    ' ctocprj : BUG -',ch10,&
&    '  cprj are badly allocated !'
    call wrtout(06,message,'PERS')
    call leave_new('PERS')
   end if

!  Retreive (k+G) points and spherical harmonics
   allocate(ylm_k(npw_k,mpsang*mpsang),kg_k(3,npw_nk))
   if (mkmem==0) then
    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,dtfil%unkg)
    read (dtfil%unkg) kg_k(1:3,1:npw_k)
    read(dtfil%unylm)
    read(dtfil%unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
    call status(counter,dtfil%filstat,iexit,level,'read wfs      ')
    tim_rwwf=1;allocate(eig_dum(mband),kg_dum(3,0),occ_dum(mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,nband_k,&
&    nband_k,npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
    deallocate(eig_dum,kg_dum,occ_dum)
   else
    if (nprocband>1) then
!    Special treatment for band-fft //
     allocate(kg_k_loc(3,npw_k))
     kg_k_loc(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     bufsize(:)=3*npw_block(:);bufdisp(:)=3*npw_disp(:)
     call xallgatherv_mpi(kg_k_loc,3*npw_k,kg_k,bufsize,bufdisp,spaceComm_band,ierr)
    else
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    end if
    do ilm=1,mpsang*mpsang
     ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
    end do
   end if

!  Compute (k+G) vectors
   allocate(kpg_k(npw_nk,nkpg))
   if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_nk)

!  Allocate and compute the arrays phkxred and ph3d
   allocate(phkxred(2,ncprj))
   do ia=iatom1,iatom2
    iatm=min(atindx(ia),ncprj)
    arg=two_pi*(kpoint(1)*xred(1,iatm)+kpoint(2)*xred(2,iatm)+kpoint(3)*xred(3,iatm))
    phkxred(1,ia)=cos(arg);phkxred(2,ia)=sin(arg)
   end do
   if(nloalg(1)<=0)then
!   Here, only the allocation of ph3d , not the precomputation
    matblk=min(nloalg(4),ncprj);allocate(ph3d(2,npw_k,matblk))
    if (nprocband>1) stop "ctocprj: ERROR - banf-fft parallelism +nloag(1)<0 forbidden !"
   else
!   Here, allocation as well as precomputation
    matblk=ncprj;allocate(ph3d(2,npw_nk,matblk))
    if (nprocband>1) then
!    Special treatment for band-fft //
     allocate(ph3d_tmp(2,npw_k,matblk))
     if (iatom<=0) then
      call ph1d3d(1,natom,kg_k_loc,kpoint,matblk,natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d_tmp)
     else
      call ph1d3d(1,1,kg_k_loc,kpoint,matblk,1,npw_k,n1,n2,n3,phkxred,ph1d_atm,ph3d_tmp)
     end if
     allocate(ph3d_tmp_npw(2,matblk,npw_k),ph3d_npw(2,matblk,npw_nk))
     isize=2*matblk;bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
     do ipw=1,npw_k
      ph3d_tmp_npw(:,:,ipw)=ph3d_tmp(:,ipw,:)
     end do
     call xallgatherv_mpi(ph3d_tmp_npw,isize*npw_k,ph3d_npw,bufsize,bufdisp,spaceComm_band,ierr)
     do ipw=1,npw_nk
      ph3d(:,ipw,:)=ph3d_npw(:,:,ipw)
     end do
     deallocate(ph3d_npw,ph3d_tmp_npw,ph3d_tmp)
    else
     if (iatom<=0) then
      call ph1d3d(1,natom,kg_k,kpoint,matblk,natom,npw_k,n1,n2,n3,phkxred,ph1d,ph3d)
     else
      call ph1d3d(1,1,kg_k,kpoint,matblk,1,npw_k,n1,n2,n3,phkxred,ph1d_atm,ph3d)
     end if
    end if
   end if

!  Compute nonlocal form factors ffnl at all (k+G)
   call status(0,dtfil%filstat,iexit,level,'call mkffnl   ')
   allocate(ffnl(npw_nk,dimffnl,psps%lmnmax,ntypat0))
   if (nprocband>1) then
!   Special treatment for band-fft //
    allocate(ffnl_tmp(npw_k,dimffnl,psps%lmnmax,ntypat0))
    if (iatom<=0) then
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl_tmp,psps%ffspl,&
&     gmet,gprimd,ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat,&
&     psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
    else
     call mkffnl(psps%dimekb,dimffnl,ekb_atm,ffnl_tmp,ffspl_atm,&
&     gmet,gprimd,ider,idir0,indlmn_atm,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat0,&
&     pspso_atm,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
    end if
    allocate(ffnl_tmp_npw(dimffnl,psps%lmnmax,ntypat0,npw_k))
    allocate(ffnl_npw(dimffnl,psps%lmnmax,ntypat0,npw_nk))
    isize=dimffnl*psps%lmnmax*ntypat0
    bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
    do ipw=1,npw_k
     ffnl_tmp_npw(:,:,:,ipw)=ffnl_tmp(ipw,:,:,:)
    end do
    call xallgatherv_mpi(ffnl_tmp_npw,isize*npw_k,ffnl_npw,bufsize,bufdisp,spaceComm_band,ierr)
    do ipw=1,npw_nk
     ffnl(ipw,:,:,:)=ffnl_npw(:,:,:,ipw)
    end do
    deallocate(ffnl_npw,ffnl_tmp_npw,ffnl_tmp)
   else
    if (iatom<=0) then
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&     gmet,gprimd,ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat,&
&     psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
    else
     call mkffnl(psps%dimekb,dimffnl,ekb_atm,ffnl,ffspl_atm,&
&     gmet,gprimd,ider,idir0,indlmn_atm,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,ntypat0,&
&     pspso_atm,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
    end if
   end if

!  No more need of kg_g_tmp
   if (nprocband>1) deallocate(kg_k_loc)

!  Allocate arrays for a wave-function (or a block of WFs)
   allocate(cwavef(2,npw_nk*nspinor))
   if (nprocband>1) then
    isize=2*nspinor;bufsize(:)=isize*npw_block(:);bufdisp(:)=isize*npw_disp(:)
    isize=2*nspinor*npw_k;bufsize_wf(:)=isize
    do ii=1,nprocband
     bufdisp_wf(ii)=(ii-1)*isize
    end do
   end if

!  Loop over bands or blocks of bands
   icgb=icg
   nblockbd=nband_k/nprocband
   do iblockbd=1,nblockbd
    iband_min=1+(iblockbd-1)*nprocband
    iband_max=iblockbd*nprocband

    if(mpi_enreg%paral_compil_kpt==1.and.mpi_enreg%mode_para/='b') then
     if (minval(abs(mpi_enreg%proc_distrb(ikpt,iband_min:iband_max,isppol)-me_distrb))/=0) cycle
    end if

!   Extract wavefunction information
    if (mkmem==0) then
     do ig=1,npw_k*nspinor
      cwavef(1,ig)=cg_disk(1,ig+icgb)
      cwavef(2,ig)=cg_disk(2,ig+icgb)
     end do
    else
     if (nprocband>1) then
!     Special treatment for band-fft //
      allocate(cwavef_tmp(2,npw_k*nspinor*nprocband))
      do ig=1,npw_k*nspinor*nprocband
       cwavef_tmp(1,ig)=cg(1,ig+icgb)
       cwavef_tmp(2,ig)=cg(2,ig+icgb)
      end do
      call xalltoallv_mpi(cwavef_tmp,bufsize_wf,bufdisp_wf,cwavef,bufsize,bufdisp,spaceComm_band,ierr)
      deallocate(cwavef_tmp)
     else
      do ig=1,npw_k*nspinor
       cwavef(1,ig)=cg(1,ig+icgb)
       cwavef(2,ig)=cg(2,ig+icgb)
      end do
     end if
    end if

!   Compute scalar product of wavefunction with all NL projectors
    if (iatom<=0) then
     call getcprj(choice,cpopt,cwavef,cwaveprj,psps%dimekb,ntypat,dimffnl,&
&     psps%ekb,ffnl,idir,psps%indlmn,istwf_k,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     matblk,mgfft,mpi_enreg,natom,nattyp,ngfft,nkpg,nloalg,&
&     npw_nk,nspinor,ntypat,phkxred,ph1d,ph3d,ucvol,psps%usepaw,psps%useylm)
    else
     call getcprj(choice,cpopt,cwavef,cwaveprj,psps%dimekb,1,dimffnl,&
&     ekb_atm,ffnl,idir,indlmn_atm,istwf_k,kg_k,kpg_k,kpoint,psps%lmnmax,&
&     matblk,mgfft,mpi_enreg,1,nattyp_atm,ngfft,nkpg,nloalg,&
&     npw_nk,nspinor,1,phkxred,ph1d_atm,ph3d,ucvol,psps%usepaw,psps%useylm)
    end if

!   === Output of the cprj
!   =======================
    call cprj_put(.false.,atindx,cwaveprj,cprj,ncprj,iband_min,ibg,ikpt,iorder_cprj,isppol,&
&    mband,mkmem,mpi_enreg,natom,1,nband_k,dimlmn,nspinor,nsppol,spaceComm_band,uncp)

!   End loop over bands
    icgb=icgb+npw_k*nspinor*nprocband
   end do

!  Shift array memory (if mkmem/=0)
   if (mkmem/=0) then
    ibg=ibg+nspinor*nband_k
    icg=icg+nspinor*nband_k*npw_k
    ikg=ikg+npw_k
   end if

!  End big k point loop
   deallocate(ffnl,ph3d,phkxred,kg_k,kpg_k,ylm_k,cwavef)
  end do
! End loop over spins
 end do

!Deallocate temporary storage
 call cprj_free(cwaveprj)
 deallocate(cwaveprj,dimlmn)
 if (nprocband>1) then
  deallocate(npw_block,npw_disp)
  deallocate(bufsize,bufdisp,bufsize_wf,bufdisp_wf)
 end if
 if(mkmem==0) deallocate(cg_disk)
 if (iatom>0) deallocate(ekb_atm,indlmn_atm,ffspl_atm,ph1d_atm)

 end subroutine ctocprj
!!***
