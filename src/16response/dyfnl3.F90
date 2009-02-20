!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyfnl3
!! NAME
!! dyfnl3
!!
!! FUNCTION
!! Compute the frozen-wavefunction non-local contribution to the
!! dynamical matrix.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GM, AR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>=Fourier coefficients of WF
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)=<p_lmn|C> coefficients for WF |C> (and 1st derivatives)
!!  dimcprj(natom)=array of dimensions of array cprj
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha) (NOT NEEDED !)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fform=integer specifier for wf file form
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=work array for coordinates of G vectors in basis
!!  kptns(3,nkpt)=coordinates of k points in terms of reciprocal space
!!   primitive translations
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimension for number of planewaves
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=array describing how many atoms of each type in cell
!!  nband(nkpt*nsppol)=number of bands being considered per k point
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  nkpt=number of k points
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  nsym=number of symmetry elements in space group (at least 1)
!!  ntypat=integer specification of atom type (1, 2, ...)
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2) at each
!!    k point
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries in reciprocal space (dimensionless)
!!  typat(natom)=type integer for each atom in cell
!!  unkg=unit number for (k+G) sphere data
!!  unpaw=unit number for cprj PAW data (if used)
!!  unylm=unit number for disk file containing Ylm''s if mkmem==0
!!  usecprj=1 if cprj coefficients are already in memory (PAW only)
!!  vtrial(nfftf,nspden)=potential (Hartree+XC+loc)
!!  wfftgs=struct info for disk file containing GS wavefunctions if mkmem==0
!!  wtk(nkpt)=k point weights
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!
!! OUTPUT
!!  dyfrnl(3,3,natom)=non-symmetrized non-local contribution to the
!!                    dynamical matrix
!!
!! SIDE EFFECTS
!!  ===== if psps%usepaw==1
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!                          pawfgrtab(:)%gylmgr2 are deallocated here
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      cprj_diskinit,cprj_get,cprj_alloc,cprj_free,hdr_skip,leave_test,metric,mkffnl,mkkpg,nonlop,pawgrnl,ph1d3d,rdnpw,rwwf
!!      sphereboundary,timab,xcomm_world,xdefineoff,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dyfnl3(atindx1,cg,cprj,dimcprj,dyfrnl,ecut,ecutsm,&
&  eigen,fform,indsym,istwfk,&
&  kg,kptns,mband,mgfft,mkmem,mpi_enreg,mpsang,&
&  mpw,natom,nattyp,nband,&
&  nfftf,ngfft,ngfftf,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,nsym,ntypat,occ,&
&  paw_ij,pawang,pawprtvol,pawfgrtab,pawrhoij,pawtab,&
&  ph1d,prtvol,psps,rprimd,symafm,symrec,typat,unkg,unpaw,unylm,usecprj,wfftgs,&
&  vtrial,wtk,xred,ylm)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_14iowfdenpot
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fform,mband,mgfft,mkmem,mpsang,mpw,natom,nfftf,nkpt
 integer,intent(in) :: nspden,nsppol,nsym,ntypat,pawprtvol,prtvol,unkg,unpaw
 integer,intent(in) :: unylm,usecprj
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ecut,ecutsm
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wfftgs
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom),indsym(4,nsym,natom)
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nattyp(ntypat)
 integer,intent(in) :: nband(nkpt*nsppol),ngfft(18),ngfftf(18),nloalg(5)
 integer,intent(in) :: npwarr(nkpt),symafm(nsym),symrec(3,3,nsym),typat(ntypat)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),vtrial(nfftf,nspden),wtk(nkpt)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(out) :: dyfrnl(3,3,natom)
 type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,bufdim,choice,cplex,cpopt,dimdij,dimekb1,dimekb2
 integer :: dimffnl,dimnhat,formeig,ia,iatom,iband,ibg,ibsp,icg,ider,idir,ierr
 integer :: ii,ikg,ikpt,ilm,iorder_cprj,iproc,ipw,ipw1,isp,ispden,isppol
 integer :: istwf_k,it,itypat,klmn,master,matblk,mcg_disk,me,mu,n1,n2,n3,nband1
 integer :: nband_k,ngradij,nkpg,nnlout,npw_k,nsploop,optgr,optgr2,option
 integer :: option_rhoij,optstr,paw_opt,signs,spaceworld,tim_nonlop,tim_rwwf
 real(dp) :: arg,dotr,eig_k,enl,occ_k,ucvol,wtk_k
 character(len=500) :: message
!arrays
 integer,allocatable :: dimlmn(:),kg_dum(:,:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),grhoij(3),kpoint(3),nonlop_dum(1,1)
 real(dp) :: rmet(3,3),tsec(2),ylmgr_dum(1)
 real(dp),allocatable :: buffer1(:),buffer2(:),cg_disk(:,:),cwavef(:,:)
 real(dp),allocatable :: dummy(:),dyfrnlk(:,:),eig_dum(:),ekb(:,:,:),enlout(:)
 real(dp),allocatable :: ffnl(:,:,:,:),kpg_k(:,:),nhat_dum(:,:),occ_dum(:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),sij(:,:),ylm_k(:,:)
 type(cprj_type),allocatable :: cprj_disk(:,:),cwaveprj(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' dyfnl3 : enter '
!stop
!ENDDEBUG

 call timab(159,1,tsec)

!Default for sequential use
 master=0
 call xme_init(mpi_enreg,me)
!Init spaceworld
 call xcomm_world(mpi_enreg,spaceworld)

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Prepare temporary files if mkmem==0
!Wavefunction file
 if (mkmem==0) then
! Read header
  call hdr_skip(wfftgs,ierr)
! Define offsets, in case of MPI I/O
  formeig=0
  call xdefineOff(formeig,wfftgs,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)
  mcg_disk=mpw*nspinor*mband
  allocate(cg_disk(2,mcg_disk))
 end if
!PAW file
 iorder_cprj=0
 if (usecprj==1) then
  call cprj_diskinit(atindx1,natom,iorder_cprj,mkmem,natom,dimcprj,nspinor,unpaw)
 end if

 enl=zero
 dyfrnl(:,:,:)=zero
 bdtot_index=0
 ibg=0;icg=0
 nsploop=nsppol;if (nspden==4) nsploop=4

!DEBUG
!write(6,*)' dyfnl3 : before loop over spins '
!stop
!ENDDEBUG

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 allocate(kg_k(3,mpw))
 allocate(cwavef(2,mpw*nspinor))
 allocate(dyfrnlk(6,natom))
 allocate(phkxred(2,natom))

!Common data for "nonlop" routine
 signs=1 ; idir=0 ; tim_nonlop=6

!Common data for "nonlop" routine
 signs=1 ; idir=0  ; tim_nonlop=4
 choice=4 ; nnlout=max(1,6*natom)
 allocate(enlout(nnlout))
 if (psps%usepaw==0) then
  paw_opt=0 ; cpopt=-1
 else
  paw_opt=2 ; cpopt=1+2*usecprj
 end if

!Non-local factors:
!Norm-conserving: kleimann-Bylander energies
!PAW: Dij coefficients and overlap coefficients
 if (psps%usepaw==0) then
  dimekb1=psps%dimekb;dimekb2=ntypat
  allocate(ekb(psps%dimekb,ntypat,nspinor**2))
  ekb(:,:,1)=psps%ekb(:,:)
  if (nspinor==2) then
   ekb(:,:,2)=psps%ekb(:,:)
   ekb(:,:,3:4)=zero
  end if
 else
  dimekb1=psps%dimekb*paw_ij(1)%cplex_dij;dimekb2=natom
  allocate(ekb(dimekb1,dimekb2,nspinor**2),sij(dimekb1,ntypat))
  do itypat=1,ntypat
   if (paw_ij(1)%cplex_dij==1) then
    sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
   else
    do klmn=1,pawtab(itypat)%lmn2_size
     sij(2*klmn-1,itypat)=pawtab(itypat)%sij(klmn)
     sij(2*klmn  ,itypat)=zero
    end do
   end if
  end do
  allocate(cwaveprj(natom,nspinor))
  call cprj_alloc(cwaveprj,3,dimcprj)
  do iatom=1,natom
   pawrhoij(iatom)%ngrhoij=3
   allocate(pawrhoij(iatom)%grhoij(3,pawrhoij(iatom)%lmn2_size,nspden))
   pawrhoij(iatom)%grhoij=zero
  end do
  option_rhoij=3
 end if

!LOOP OVER SPINS
 do isppol=1,nsppol


! Rewind kpgsph data file if needed:
  if (mkmem==0) rewind(unkg)
  if (mkmem==0.and.psps%useylm==1) rewind unylm

! PAW: retrieve Dij coefficients
  if (psps%usepaw==1) then
   do ispden=1,nspinor**2
    isp=isppol;if (nspinor==2) isp=ispden
    do iatom=1,natom
     dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
     do klmn=1,dimdij
      ekb(klmn,iatom,ispden)=paw_ij(iatom)%dij(klmn,isp)
     end do
     if(dimdij+1<=dimekb1) ekb(dimdij+1:dimekb1,iatom,ispden)=zero
    end do
   end do
  end if

  ikg=0

! Loop over k points
  do ikpt=1,nkpt
   nband_k=nband(ikpt+(isppol-1)*nkpt)
   istwf_k=istwfk(ikpt)
   npw_k=npwarr(ikpt)
   wtk_k=wtk(ikpt)

   if(mpi_enreg%paral_compil_kpt==1)then
!   Skip this k-point if not the proper processor
    if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
     bdtot_index=bdtot_index+nband_k
     cycle
    end if
   end if

   allocate(ylm_k(npw_k,mpsang*mpsang*psps%useylm))
   kpoint(:)=kptns(:,ikpt)

   kg_k(:,:) = 0
   if (mkmem==0) then

    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,unkg)
!   Read k+g data
    read (unkg) ((kg_k(ii,ipw),ii=1,3),ipw=1,npw_k)
!   Eventually read spherical harmonics
    if (psps%useylm==1) then
     read(unylm)
     read(unylm) ((ylm_k(ipw,ilm),ipw=1,npw_k),ilm=1,mpsang*mpsang)
    end if
!   Read the wavefunction block for ikpt,isppol
    tim_rwwf=14
    allocate(eig_dum(mband),kg_dum(3,0),occ_dum(mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,&
&    nband_k,nband_k,npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wfftgs)
    deallocate(eig_dum,kg_dum,occ_dum)

   else

!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(ikg,kg,kg_k,npw_k)
    do ipw=1,npw_k
     kg_k(1,ipw)=kg(1,ipw+ikg)
     kg_k(2,ipw)=kg(2,ipw+ikg)
     kg_k(3,ipw)=kg(3,ipw+ikg)
    end do
!   $OMP END PARALLEL DO
    if (psps%useylm==1) then
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(ikg,npw_k,ylm,ylm_k)
     do ilm=1,mpsang*mpsang
      do ipw=1,npw_k
       ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
      end do
     end do
!    $OMP END PARALLEL DO
    end if

!   End if for choice governed by mkmem
   end if

   cplex=2;if (istwf_k>1) cplex=1

!  Extract PAW cprj quantities according to mkmem
   if (mkmem==0.and.usecprj==1) then
    allocate(cprj_disk(natom,nspinor*nband_k))
    call cprj_alloc(cprj_disk,0,dimcprj)
    call cprj_get(atindx1,cprj_disk,cprj,natom,1,ibg,ikpt,iorder_cprj,isppol,&
&    mband,mkmem,mpi_enreg,natom,nband_k,nband_k,nspinor,nsppol,unpaw)
   end if

!  Compute nonlocal psp energy

!  Compute (k+G) vectors (only if useylm=1)
   nkpg=9*nloalg(5);allocate(kpg_k(npw_k,nkpg))
   if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!  Compute nonlocal form factors ffnl at all (k+G):
   ider=0;idir=0;dimffnl=1
   allocate(ffnl(npw_k,dimffnl,psps%lmnmax,ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&   gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&   ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)

   dyfrnlk(:,:)=zero

!  Compute phkxred and eventually ph3d.
   do iatom=1,natom
    ia=atindx1(iatom)
    arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
    phkxred(1,iatom)=cos(arg)
    phkxred(2,iatom)=sin(arg)
   end do
   if(nloalg(1)<=0)then
!   Only the allocation, not the precomputation.
    matblk=nloalg(4)
    allocate(ph3d(2,npw_k,matblk))
   else
!   Here, allocation as well as precomputation
    matblk=natom
    allocate(ph3d(2,npw_k,matblk))
    call ph1d3d(1,natom,kg_k,kpoint,matblk,natom,npw_k,n1,n2,n3,&
&    phkxred,ph1d,ph3d)
   end if

!  DEBUG
!  write(6,*)' dyfnl3 : before nonlop '
!  stop
!  ENDDEBUG

   do iband=1,nband_k

    if(mpi_enreg%paral_compil_kpt==1)then
     if(mpi_enreg%proc_distrb(ikpt, iband,isppol) /= me) then
      cycle
     end if
    end if

    occ_k=occ(iband+bdtot_index)

    if(mkmem/=0)then
     cwavef(:,1:npw_k*nspinor)=&
&     cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
    else
     cwavef(:,1:npw_k*nspinor)=&
&     cg_disk(:,1+(iband-1)*npw_k*nspinor:iband*npw_k*nspinor)
    end if

    if (psps%usepaw==1.and.usecprj==1) then
     if (mkmem/=0) then
      ibsp=(iband-1)*nspinor+ibg
      call cprj_copy(cprj     (:,ibsp+1:ibsp+nspinor),cwaveprj)
     else
      ibsp=(iband-1)*nspinor
      call cprj_copy(cprj_disk(:,ibsp+1:ibsp+nspinor),cwaveprj)
     end if
    end if

!   Compute non-local contributions from n,k
    if (psps%usepaw==1) eig_k=eigen(iband+bdtot_index)
    call nonlop(atindx1,choice,cpopt,cwaveprj,dimekb1,dimekb2,dimffnl,dimffnl,ekb,&
&    enlout,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&    kpoint,eig_k,psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,nattyp,&
&    ngfft,nkpg,nkpg,nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,phkxred,&
&    phkxred,ph1d,ph3d,ph3d,psps%pspso,signs,sij,nonlop_dum,tim_nonlop,ucvol,&
&    psps%useylm,cwavef,cwavef)

!   Accumulate non-local contributions from n,k
    dyfrnlk(:,:)=dyfrnlk(:,:)+occ_k*reshape(enlout(:),(/6,natom/))

!   PAW: accumulate gradients of rhoij
    if (psps%usepaw==1) then
     call pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj,natom,0,isppol,natom,nspden,&
&     nspinor,nsppol,occ_k,option_rhoij,pawrhoij,wtk_k)
    end if

!   End of loop on bands
   end do

   dyfrnl(1,1,:)=dyfrnl(1,1,:)+wtk_k*dyfrnlk(1,:)/ucvol
   dyfrnl(2,2,:)=dyfrnl(2,2,:)+wtk_k*dyfrnlk(2,:)/ucvol
   dyfrnl(3,3,:)=dyfrnl(3,3,:)+wtk_k*dyfrnlk(3,:)/ucvol
   dyfrnl(2,3,:)=dyfrnl(2,3,:)+wtk_k*dyfrnlk(4,:)/ucvol
   dyfrnl(1,3,:)=dyfrnl(1,3,:)+wtk_k*dyfrnlk(5,:)/ucvol
   dyfrnl(1,2,:)=dyfrnl(1,2,:)+wtk_k*dyfrnlk(6,:)/ucvol

!  Incremente indexes
   bdtot_index=bdtot_index+nband_k
   if (mkmem/=0) then
    ibg=ibg+nband_k*nspinor
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   else if (usecprj==1) then
    call cprj_free(cprj_disk)
    deallocate(cprj_disk)
   end if

   deallocate(ffnl,kpg_k,ph3d,ylm_k)

!  End loops on isppol and ikpt
  end do
 end do

 deallocate(phkxred,cwavef,dyfrnlk,kg_k,enlout,ekb)
 if (psps%usepaw==1) then
  call cprj_free(cwaveprj)
  deallocate(cwaveprj)
  deallocate(sij)
 end if
 if(mkmem==0) deallocate(cg_disk)

 do iatom=1,natom
  dyfrnl(3,2,iatom)=dyfrnl(2,3,iatom)
  dyfrnl(3,1,iatom)=dyfrnl(1,3,iatom)
  dyfrnl(2,1,iatom)=dyfrnl(1,2,iatom)
 end do

!Parallel case: accumulate (n,k) contributions
 if( mpi_enreg%paral_compil_kpt==1) then
  call timab(48,1,tsec)
  call leave_test(mpi_enreg)
! Norm-conserving: accumulate dyfrnl
  if (psps%usepaw==0) then
   call xsum_mpi(dyfrnl,spaceworld,ierr)
  else
!  PAW: accumulate gradients of rhoij
   allocate(dimlmn(natom));dimlmn(1:natom)=pawrhoij(1:natom)%lmn2_size
   bufdim=3*sum(dimlmn)*nsploop
   allocate(buffer1(bufdim),buffer2(bufdim))
   ii=0
   do iatom=1,natom
    do isppol=1,nsploop
     do mu=1,3
      buffer1(ii+1:ii+dimlmn(iatom))=pawrhoij(iatom)%grhoij(mu,1:dimlmn(iatom),isppol)
      ii=ii+dimlmn(iatom)
     end do
    end do
   end do
   call xsum_mpi(buffer1,buffer2,bufdim,spaceworld,ierr)
   ii=0
   do iatom=1,natom
    do isppol=1,nsploop
     do mu=1,3
      pawrhoij(iatom)%grhoij(mu,1:dimlmn(iatom),isppol)=buffer2(ii+1:ii+dimlmn(iatom))
      ii=ii+dimlmn(iatom)
     end do
    end do
   end do
   deallocate(buffer1,buffer2,dimlmn)
  end if
  call timab(48,2,tsec)
 end if

 if (psps%usepaw==1) then
! PAW: symmetrize rhoij gradients and transfer to cart. coord.
  choice=2;option=0  ! This symetrization is useful in the antiferromagnetic case...
  call symrhoij(choice,psps%indlmn,indsym,psps%lmnmax,natom,nsym,ntypat,option,&
&  pawang,pawprtvol,pawrhoij,symafm,symrec,typat)
  do iatom=1,natom
   do isppol=1,nsploop
    do klmn=1,pawrhoij(iatom)%lmn2_size
     grhoij(1:3)=pawrhoij(iatom)%grhoij(1:3,klmn,isppol)
     do mu=1,3
      pawrhoij(iatom)%grhoij(mu,klmn,isppol)=gprimd(mu,1)*grhoij(1) &
&      +gprimd(mu,2)*grhoij(2)+gprimd(mu,3)*grhoij(3)
     end do
    end do
   end do
  end do

! PAW: Add gradients due to Dij derivatives to dynamical matrix
  dimnhat=0;optgr=0;optgr2=1;optstr=0;allocate(nhat_dum(1,0))
  call pawgrnl(atindx1,dimnhat,nspden,dyfrnl,dummy,mpi_enreg,natom,nattyp,nfftf,ngfftf,&
&  nhat_dum,dummy,nspden,nsym,ntypat,optgr,optgr2,optstr,pawang,&
&  pawfgrtab,pawrhoij,pawtab,rprimd,symrec,typat,vtrial)
  deallocate(nhat_dum)

 end if

 call timab(159,2,tsec)

!DEBUG
!write(6,*)' dyfnl3 : exit '
!stop
!ENDDEBUG

end subroutine dyfnl3
!!***
