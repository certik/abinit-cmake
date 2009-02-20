!{\src2tex{textfont=tt}}
!!****f* ABINIT/eltfrnl3
!! NAME
!! eltfrnl3
!!
!! FUNCTION
!! Compute the frozen-wavefunction non-local contribution to the
!! elastic tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DRH, DCA, XG, GM, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=<G|Cnk>
!!              =Fourier coefficients of wavefunction
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha) (NOT NEEDED !)
!!  fform=integer specifier for wf file form
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
!!  nkpt=number of k points
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves at each k point, and boundary
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for polarized
!!  nsym=number of symmetry elements in space group (at least 1)
!!  ntypat=integer specification of atom type (1, 2, ...)
!!  occ(mband*nkpt*nsppol)=occupation numbers of bands (usually 2)
!!    at each k point
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase information related to structure factor
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  typat(natom)=integer specification of atom type for each atom
!!  unkg=unit number for (k+G) sphere data
!!  wfftgs=struct info for disk file containing GS wavefunctions if mkmem==0
!!  unylm=unit number for disk file containing Ylm''s if mkmem==0
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!  wtk(nkpt)=k point weights
!!  xred(3,natom)=reduced coordinates of atoms (dimensionless)
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics
!!    for each G and k point
!!  ylmgr(mpw*mkmem,9,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!    for each G and k point
!!
!! OUTPUT
!!  eltfrnl(6+3*natom,6)=non-symmetrized non-local contribution to the
!!                    elastic tensor
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      hdr_skip,leave_test,metric,mkffnl,mkkpg,nonlop,ph1d3d,rdnpw,rwwf
!!      sphereboundary,timab,xcomm_world,xdefineoff,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine eltfrnl3(atindx,atindx1,cg,eltfrnl,ecut,ecutsm,fform,istwfk,&
&  kg,kptns,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nband,&
&  nkpt,ngfft,nloalg,npwarr,nspinor,nsppol,nsym,ntypat,occ,ph1d,&
&  prtvol,psps,rprimd,typat,unkg,wfftgs,unylm,useylm,wtk,xred,ylm,ylmgr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_14iowfdenpot
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fform,mband,mgfft,mkmem,mpsang,mpw,natom,nkpt,nsppol
 integer,intent(in) :: nsym,ntypat,prtvol,unkg,unylm,useylm
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ecut,ecutsm
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wfftgs
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),istwfk(nkpt)
 integer,intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(5),npwarr(nkpt),typat(natom)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),kptns(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),wtk(nkpt),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,9,mpsang*mpsang*useylm)
 real(dp),intent(out) :: eltfrnl(6+3*natom,6)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,choice,cpopt,dimekb1,dimekb2,dimffnl,formeig,ia,iatom
 integer :: iband,icg,ider,idir,ielt,ieltx,ierr,ii,ikg,ikpt,ilm,ilmn,iln,iln0
 integer :: index,ipw,ipw1,isppol,istwf_k,jj,master,matblk,mcg_disk,me,n1,n2,n3
 integer :: nband1,nband_k,nkpg,nnlout,npw_k,option,paw_opt,signs,spaceComm
 integer :: tim_nonlop,tim_rwwf
 real(dp) :: arg,dotr,dum,enl,enlk,ucvol
 character(len=500) :: message
!arrays
 integer,allocatable :: gbound(:,:),kg_dum(:,:),kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),kpoint(3),nonlop_dum(1,1),rmet(3,3),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),cwavef(:,:),eig_dum(:),ekb(:,:,:)
 real(dp),allocatable :: elt_work(:,:),eltfrnlk(:,:),enlout(:),ffnl(:,:,:,:)
 real(dp),allocatable :: kpg_k(:,:),occ_dum(:),ph3d(:,:,:),phkxred(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr_k(:,:,:)
 type(cprj_type) :: cprj_dum(1,1)

! *************************************************************************

!DEBUG
!write(6,*)' dyfnl3 : enter '
!stop
!ENDDEBUG

!Default for sequential use
 master=0
!BEGIN TF_CHANGES
 call xme_init(mpi_enreg,me)
!END TF_CHANGES
!Init mpi_comm
!BEGIN TF_CHANGES
 call xcomm_world(mpi_enreg,spaceComm)
!END TF_CHANGES
 nnlout=6*(3*natom+6)

!Compute gmet, gprimd and ucvol from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 if (mkmem==0) then

! Read wavefunction file header
  call hdr_skip(wfftgs,ierr)

! Define offsets, in case of MPI I/O
  formeig=0
  call xdefineOff(formeig,wfftgs,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)

  mcg_disk=mpw*nspinor*mband
  allocate(cg_disk(2,mcg_disk))

 end if

 enl=0.0_dp
 eltfrnl(:,:)=0.0_dp
 bdtot_index=0
 icg=0

!Non-local factors:
!Norm-conserving: kleimann-Bylander energies
 if (psps%usepaw==0) then
  dimekb1=psps%dimekb;dimekb2=ntypat
  allocate(ekb(psps%dimekb,ntypat,nspinor**2))
  ekb(:,:,1)=psps%ekb(:,:)
  if (nspinor==2) then
   ekb(:,:,2)=psps%ekb(:,:)
   ekb(:,:,3:4)=zero
  end if
 else
! Not available within PAW
  allocate(ekb(psps%dimekb,natom,nspinor**2))
 end if

!DEBUG
!write(6,*)' dyfnl3 : before loop over spins '
!stop
!ENDDEBUG

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 allocate(kg_k(3,mpw),cwavef(2,mpw*nspinor),eltfrnlk(6+3*natom,6))

!Define k-points distribution

!LOOP OVER SPINS
 do isppol=1,nsppol


! Rewind kpgsph data file if needed:
  if (mkmem==0) rewind(unkg)
  if (mkmem==0.and.useylm==1) rewind unylm

  ikg=0

! Loop over k points
  do ikpt=1,nkpt
   nband_k=nband(ikpt+(isppol-1)*nkpt)
   istwf_k=istwfk(ikpt)
   npw_k=npwarr(ikpt)

   if(mpi_enreg%paral_compil_kpt==1)then
!   Skip this k-point if not the proper processor
!   BEGIN TF_CHANGES
    if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&    - me))/=0) then
!    END TF_CHANGES
     bdtot_index=bdtot_index+nband_k
     cycle
    end if
   end if

   allocate(gbound(2*mgfft+8,2))
   allocate(ylm_k(npw_k,mpsang*mpsang*useylm))
   allocate(ylmgr_k(npw_k,9,mpsang*mpsang*psps%useylm))
   kpoint(:)=kptns(:,ikpt)

   kg_k(:,:) = 0
   if (mkmem==0) then

    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,unkg)

!   Read k+g data
    read (unkg) ((kg_k(ii,ipw),ii=1,3),ipw=1,npw_k)

    call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

!   Eventually read (k+g) spherical harmonics
    if (psps%useylm==1) then
     read(unylm)
     read(unylm) ((ylm_k(jj,ilm),jj=1,npw_k),ilm=1,mpsang*mpsang),&
&     (((ylmgr_k(jj,ii,ilm),jj=1,npw_k),ii=1,9),ilm=1,mpsang*mpsang)
    end if

!   Read the wavefunction block for ikpt,isppol
    tim_rwwf=14
    allocate(eig_dum(mband),kg_dum(3,0),occ_dum(mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&    npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wfftgs)
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

    call sphereboundary(gbound,istwf_k,kg_k,mgfft,npw_k)

    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
!     $OMP PARALLEL DO PRIVATE(ipw) &
!     $OMP&SHARED(ikg,ilm,ylm,ylm_k_k,npw_k)
      do ipw=1,npw_k
       ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
      end do
!     $OMP END PARALLEL DO
      do ii=1,9
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(ikg,ilm,ii,ylmgr,ylmgr_k_k,npw_k)
       do ipw=1,npw_k
        ylmgr_k(ipw,ii,ilm)=ylmgr(ipw+ikg,ii,ilm)
       end do
!      $OMP END PARALLEL DO
      end do
     end do
    end if

!   End if for choice governed by mkmem
   end if

   index=1+icg

!  Compute nonlocal psp energy

!  Compute (k+G) vectors (only if useylm=1)
   nkpg=3*nloalg(5);allocate(kpg_k(npw_k,nkpg))
   if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!  Compute nonlocal form factors ffnl at all (k+G):
   ider=2;idir=0;dimffnl=3+7*psps%useylm
   allocate(ffnl(npw_k,dimffnl,psps%lmnmax,ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&   gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,&
&   ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

   enlk=0.0_dp
   eltfrnlk(:,:)=0.0_dp

!  Allocate the arrays phkxred and ph3d, compute phkxred
!  and eventually ph3d.
   allocate(phkxred(2,natom))
   do ia=1,natom
    iatom=atindx(ia)
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

   nnlout=6*(3*natom+6)
   allocate(enlout(nnlout))

   do iband=1,nband_k

    if(mpi_enreg%paral_compil_kpt==1)then
!    BEGIN TF_CHANGES
     if(mpi_enreg%proc_distrb(ikpt, iband,isppol) /= me) then
!     END TF_CHANGES
      cycle
     end if
    end if

    if(mkmem/=0)then
     cwavef(:,1:npw_k*nspinor)=&
&     cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
    else
     cwavef(:,1:npw_k*nspinor)=&
&     cg_disk(:,1+(iband-1)*npw_k*nspinor:iband*npw_k*nspinor)
    end if

    signs=1 ; choice=6 ; idir=0 ; tim_nonlop=6 ; paw_opt=0 ; cpopt=-1
    call nonlop(atindx1,choice,cpopt,cprj_dum,dimekb1,dimekb2,&
&    dimffnl,dimffnl,ekb,enlout,&
&    ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&    kpoint,dum,psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,nattyp,&
&    ngfft,nkpg,nkpg,nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,phkxred,&
&    phkxred,ph1d,ph3d,ph3d,psps%pspso,signs,&
&    nonlop_dum,nonlop_dum,tim_nonlop,ucvol,psps%useylm,&
&    cwavef,cwavef)

    eltfrnlk(:,:)=eltfrnlk(:,:)+ &
&    occ(iband+bdtot_index)* reshape(enlout(:), (/6+3*natom,6/) )
   end do

   deallocate(enlout)

   eltfrnl(:,:)=eltfrnl(:,:)+wtk(ikpt)*eltfrnlk(:,:)

   deallocate(gbound)

   bdtot_index=bdtot_index+nband_k

   if (mkmem/=0) then
!   Handle case in which kg, cg, are kept in core
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   end if

   deallocate(ffnl,kpg_k,phkxred,ph3d,ylm_k,ylmgr_k)

!  End loops on isppol and ikpt
  end do
 end do

!Fill in lower triangle
 do jj=2,6
  do ii=1,jj-1
   eltfrnl(jj,ii)=eltfrnl(ii,jj)
  end do
 end do

 if(mkmem==0)then
  deallocate(cg_disk)
 end if
!DEBUG
!write(6,*)' dyfnl3 : after loop on kpts and spins '
!stop
!ENDDEBUG

 if(mpi_enreg%paral_compil_kpt==1)then
! BEGIN TF_CHANGES
  call leave_test(mpi_enreg)
! END TF_CHANGES
! Accumulate eltfrnl on all proc.
  call timab(48,1,tsec)
  call xsum_mpi(eltfrnl,spaceComm,ierr)
  call timab(48,2,tsec)
 end if

!The indexing array atindx is used to reestablish the correct
!order of atoms
 allocate(elt_work(6+3*natom,6))
 elt_work(1:6,1:6)=eltfrnl(1:6,1:6)
 do ia=1,natom
  ielt=7+3*(ia-1)
  ieltx=7+3*(atindx(ia)-1)
  elt_work(ielt:ielt+2,1:6)=eltfrnl(ieltx:ieltx+2,1:6)
 end do
 eltfrnl(:,:)=elt_work(:,:)

 deallocate(cwavef,eltfrnlk,elt_work,kg_k,ekb)

!DEBUG
!write(6,*)' dyfnl3 : exit '
!stop
!ENDDEBUG

end subroutine eltfrnl3
!!***
