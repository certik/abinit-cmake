!{\src2tex{textfont=tt}}
!!****f* ABINIT/nstwf4
!! NAME
!! nstwf4
!!
!! FUNCTION
!! This routine computes the non-local and kinetic contribution to the
!! 2DTE matrix elements, in the non-stationary formulation
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (DRH, XG,AR,MB,MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q.
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)  (NOT NEEDED !)
!!  effmass=effective mass for electrons (1. in common case)
!!  eig_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*nsppol*mband**2)=matrix of first-order eigenvalues (hartree)
!!    (NOT NEEDED !)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  icg=shift to be applied on the location of data in the array cg
!!  icg1=shift to be applied on the location of data in the array cg1
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  ipert=type of the perturbation
!!  iscf=(<= 0  =>non-SCF), >0 => SCF
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istwf_k=option parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kpt(3)=reduced coordinates of k points.
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpert =maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nline=number of CG line minimizations per band.
!!      for the time being, also number of non-SCF passes through all bands
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  n4,n5,n6=ngfft(4),ngfft(5),ngfft(6), used for dimensioning real space arrays
!!  occopt=option for occupation numbers (not used !)
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor information
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rmet(3,3)=real space metric (bohr**2)
!!  typat(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  wffnow=struct info for INPUT 1st-order wf file
!!  wfftgs=struct info for GS wavefunction file at k
!!  wtk_k=weight assigned to the k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(npw_k,mpsang*mpsang)= real spherical harmonics for each G and k point
!!  ylm1(npw1_k,mpsang*mpsang)= real spherical harmonics for each G and k point
!!  ylmgr(npw_k,3,mpsang*mpsang*useylm)= gradients of real spherical for each G and k point
!!  ylmgr1(npw1_k,3,mpsang*mpsang*useylm)= gradients of real spherical for each G and k point
!!
!! OUTPUT
!!  d2nl_k(2,3,mpert)=non-local contributions to
!!   non-stationary 2DTE, for the present k point, and perturbation idir, ipert
!!
!! PARENTS
!!      nselt3
!!
!! CHILDREN
!!      dotprod_g,kpgstr,mkffnl,mkkin,nonlop,ph1d3d,wffreaddatarec
!!      wffreadnpwrec,wffreadskiprec,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nstwf4(atindx,atindx1,cg,cg1,d2nl_k,ecut,ecutsm,effmass,&
&  eig_k,eig1_k,gmet,gprimd,icg,icg1,idir,ikpt,ipert,&
&  iscf,isppol,istwf_k,kg_k,kg1_k,kpt,mband,mgfft,mkmem,mk1mem,mpert,&
&  mpi_enreg,mpsang,mpw,mpw1,&
&  natom,nattyp,nband_k,nfft,ngfft, &
&  nline,nloalg,npw_k,npw1_k,nspinor,nsppol,ntypat,n4,n5,n6,occopt, &
&  occ_k,ortalg,ph1d,prtvol,psps,qphon,rmet,&
&  typat,ucvol,wffnow,wfftgs,&
&  wtk_k,xred,ylm,ylm1,ylmgr,ylmgr1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12spacepar
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_13recipspace
 use interfaces_16response, except_this_one => nstwf4
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,icg1,idir,ikpt,ipert,iscf,isppol,istwf_k,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem,mpert,mpsang,mpw,mpw1,n4,n5,n6,natom,nfft
 integer,intent(in) :: nline,nsppol,ntypat,occopt,ortalg,prtvol
 integer,intent(inout) :: nband_k,npw1_k,npw_k,nspinor
 real(dp),intent(in) :: ecut,ecutsm,effmass,ucvol,wtk_k
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow,wfftgs
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),kg1_k(3,npw1_k)
 integer,intent(in) :: kg_k(3,npw_k),nattyp(ntypat),ngfft(18),nloalg(5)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: eig1_k(2*nsppol*mband**2),eig_k(mband*nsppol),gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpt(3),occ_k(nband_k)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qphon(3),rmet(3,3)
 real(dp),intent(in) :: xred(3,natom),ylm(npw_k,mpsang*mpsang)
 real(dp),intent(in) :: ylm1(npw_k,mpsang*mpsang),ylmgr(npw_k,3,mpsang*mpsang)
 real(dp),intent(in) :: ylmgr1(npw_k,3,mpsang*mpsang)
 real(dp),intent(out) :: d2nl_k(2,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimekb1,dimekb2,dimffnl,dimffnl2,i1,i2,i3,ia,iatom
 integer :: iband,ider,idir0,idir1,ier,ierr,ii,ilmn,indeig,index,ipert1,iproj
 integer :: ipsang,ipw,ipw1,ipws,ispinor,istr1,itypat,iwavef,jband,matblk
 integer :: mcgnow,mcgtgs,me,n1,n2,n3,nkpg,nnlout,paw_opt,signs,tim_nonlop
 integer :: tim_rwwf
 real(dp) :: aa,arg,dot1i,dot1r,dot2i,dot2r,dot_ndiagi,dot_ndiagr,doti,dotr,dum
 real(dp) :: im0,im1,re0,re1,scprod,tolwfr,valuei,valuer
 character(len=500) :: message
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: dummy(2,1),enlout(6),kpq(3),nonlop_dum(1,1),phkxredin(2,1)
 real(dp) :: phkxredout(2,1),tsec(2)
 real(dp),allocatable :: cg1_disk(:,:),cg_disk(:,:),cg_k(:,:),cwave0(:,:)
 real(dp),allocatable :: cwavef(:,:),cwavef_da(:,:),cwavef_db(:,:),dkinpw1(:)
 real(dp),allocatable :: eig1_disk(:),eig2_k(:),eig_dum(:),ekb(:,:,:)
 real(dp),allocatable :: ffnl(:,:,:,:),ffnl_ylm(:,:,:,:),ghc(:,:),grnk(:)
 real(dp),allocatable :: gvnl1(:,:),gvnlc(:,:),kinpw1(:),kpg_dum(:,:)
 real(dp),allocatable :: occ_dum(:),ph3d(:,:,:),phkxred(:,:)
 type(cprj_type) :: cprj_dum(1,1)
!no_abirules

! *********************************************************************


!DEBUG
!write(6,*)' nstwf3 : enter '
!stop
!ENDDEBUG

!Init me
 call xme_init(mpi_enreg,me)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 allocate(phkxred(2,natom))
 allocate(ghc(2,npw1_k*nspinor),gvnlc(2,npw1_k*nspinor))
 allocate(gvnl1(2,npw1_k*nspinor))
 allocate(eig2_k(2*nsppol*mband**2))
 allocate(kinpw1(npw1_k),dkinpw1(npw_k))
 dimffnl=2;allocate(ffnl(npw_k,dimffnl,psps%lmnmax,ntypat))
 nkpg=0

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

!Will compute the 3D phase factors inside nonlop

!Compute nonlocal form factors ffnl at (k+G), for all atoms
 if (psps%useylm==0) then
  ider=1;idir0=0
  call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gmet,gprimd,ider,idir0,&
&  psps%indlmn,kg_k,kpg_dum,kpt,psps%lmnmax,&
&  psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&  npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,ylmgr)
 else
  ider=1;idir0=-7;dimffnl2=7
  allocate(ffnl_ylm(npw_k,dimffnl2,psps%lmnmax,ntypat))
  call mkffnl(psps%dimekb,dimffnl2,psps%ekb,ffnl_ylm,psps%ffspl,gmet,gprimd,ider,idir0,&
&  psps%indlmn,kg_k,kpg_dum,kpt,psps%lmnmax,&
&  psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&  npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,ylmgr)
  do itypat=1,ntypat
   do ilmn=1,psps%lmnmax
    ffnl(:,1,ilmn,itypat)=ffnl_ylm(:,1,ilmn,itypat)
   end do
  end do
 end if

!For strain, perturbation wavevector = 0
 kpq(:)=kpt(:)

 call mkkin(ecut,ecutsm,effmass,gmet,kg1_k,kinpw1,kpq,npw1_k)

!Allocate the arrays phkxred and ph3d, compute phkxred
!and eventually ph3d.
!NOTE : in this RF case, uses kpq instead of kpt
 do ia=1,natom
  iatom=atindx(ia)
  arg=two_pi*(kpq(1)*xred(1,ia)+kpq(2)*xred(2,ia)+kpq(3)*xred(3,ia))
  phkxred(1,iatom)=cos(arg)
  phkxred(2,iatom)=sin(arg)
 end do

!Note : use npw1_k
 if(nloalg(1)<=0)then
! Here, only the allocation, not the precomputation.
  matblk=nloalg(4)
  allocate(ph3d(2,npw1_k,matblk))
 else
! Here, allocation as well as precomputation
  matblk=natom
  allocate(ph3d(2,npw1_k,matblk))
  call ph1d3d(1,natom,kg1_k,kpq,matblk,natom,npw1_k,n1,n2,n3,&
&  phkxred,ph1d,ph3d)
 end if

 allocate(cwave0(2,npw_k*nspinor),cwavef(2,npw1_k*nspinor))

!Take care of the npw record
!NOTE : one should be able to modify the rwwf routine to take care
!of the band parallelism, which is not the case yet ...
 if (mkmem==0) then
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_k,nspinor,wfftgs)
! Skip k+G and eigenvalue records in wfftgs (already in eigen0)
  call WffReadSkipRec(ierr,2,wfftgs)
 end if
 if (mk1mem==0) then
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor,wffnow)
! Skip k+G record
  call WffReadSkipRec(ierr,1,wffnow)
 end if

!Loop over bands
 do iband=1,nband_k

  if(mpi_enreg%paral_compil_kpt==1)then
!  BEGIN TF_CHANGES
   if(mpi_enreg%proc_distrb(ikpt, iband, isppol) /= me) then
!   END TF_CHANGES
!   Skip the eigenvalue and the gvnl records of this band
    if(mkmem==0)then
     call WffReadSkipRec(ierr,1,wfftgs)
    end if
    if(mk1mem==0)then
     call WffReadSkipRec(ierr,2,wffnow)
    end if
    cycle
   end if
  end if

! Read ground-state and first-order wavefunctions
  if(mkmem/=0)then
   cwave0(:,:)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
  else
   call WffReadDataRec(cwave0,ierr,2*npw_k*nspinor,wfftgs)
  end if

  if(mk1mem/=0)then
   cwavef(:,:)=cg1(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)
  else
!  Skip the eigenvalue line
   call WffReadSkipRec(ierr,1,wffnow)
   call WffReadDataRec(cwavef,ierr,2*npw1_k*nspinor,wffnow)
  end if

! Double loop over strain perturbations
  do ipert1=natom+3,natom+4
   do idir1=1,3
    if(ipert1==natom+3) then
     istr1=idir1
    else
     istr1=idir1+3
    end if

!   Compute the derivative of the kinetic operator vs strain in dkinpw
    call kpgstr(dkinpw1,ecut,ecutsm,effmass,gmet,gprimd,istr1,&
&    kg1_k,kpq,npw1_k)

!   Get |vnon-locj1|u0> :
!   first-order non-local, applied to zero-order wavefunction
!   (??) this routine gives MINUS the non-local contribution

    signs=2 ; choice=3 ; nnlout=6 ; paw_opt=0 ; cpopt=-1

!   When using Ylms, load the correct ffnl derivative
    if (psps%useylm==1) then
     do itypat=1,ntypat
      do ilmn=1,psps%lmnmax
       ffnl(:,2,ilmn,itypat)=ffnl_ylm(:,1+istr1,ilmn,itypat)
      end do
     end do
    end if

    tim_nonlop=5
    call nonlop(atindx1,choice,cpopt,cprj_dum,dimekb1,dimekb2,dimffnl,dimffnl,ekb,enlout,&
&    ffnl,ffnl,gmet,gprimd,&
&    istr1,psps%indlmn,istwf_k,kg_k,kg1_k,kpg_dum,kpg_dum,kpt,kpq,dum,psps%lmnmax,&
&    matblk,mgfft,mpi_enreg,psps%mpsang,psps%mpssoang,&
&    natom,nattyp,ngfft,nkpg,nkpg,nloalg,nnlout,npw_k,npw1_k,nspinor,&
&    ntypat,0,paw_opt,phkxred,phkxred,ph1d,ph3d,ph3d,psps%pspso,&
&    signs,nonlop_dum,nonlop_dum,&
&    tim_nonlop,ucvol,psps%useylm,cwave0,gvnl1)

!   <G|Vnl1|Cnk> is contained in gvnl1

!   Kinetic contribution
    do ispinor=1,nspinor
     do ipw=1,npw1_k
      ipws=ipw+npw1_k*(ispinor-1)
      if(kinpw1(ipw)<huge(0.0_dp)*1.d-11)then
       gvnl1(1,ipws)=gvnl1(1,ipws)+dkinpw1(ipw)*cwave0(1,ipws)
       gvnl1(2,ipws)=gvnl1(2,ipws)+dkinpw1(ipw)*cwave0(2,ipws)
      else
       gvnl1(1,ipws)=0.0_dp
       gvnl1(2,ipws)=0.0_dp
      end if
     end do
    end do

!   construct the matrix element (<uj2|vj1|u0>)complex conjug.
!   and add it to the 2nd-order matrix
!   imaginary term should be zero for strain-strain 2nd derivatives,
!   but keep it as a test for now
!   XG030513 : use dotprod_g, for future parallelisation
    call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*nspinor,2,cwavef,gvnl1)
    d2nl_k(1,idir1,ipert1)= &
&    d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*2.0_dp*dotr
    d2nl_k(2,idir1,ipert1)=&
&    d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*2.0_dp*doti

   end do !idir1
  end do !ipert1

! UNTIL NOW, DO NOT TAKE INTO ACCOUNT istwf_k

! End loop over bands
 end do

 deallocate(cwave0,cwavef)

!###################################################################

 deallocate(eig2_k,ghc,gvnlc,gvnl1)
 deallocate(dkinpw1)
 deallocate(ffnl,kinpw1,phkxred,ph3d,ekb)
 if (psps%useylm==1) deallocate(ffnl_ylm)


end subroutine nstwf4
!!***
