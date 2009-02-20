!{\src2tex{textfont=tt}}
!!****f* ABINIT/nstwf3
!! NAME
!! nstwf3
!!
!! FUNCTION
!! This routine computes the non-local contribution to the
!! 2DTE matrix elements, in the non-stationary formulation
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG,AR,MB,MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see scfcv.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of
!!     RF wavefunctions at k,q.
!!  ddkfil(3)=unit numbers for the three possible ddk files for ipert1
!!       equal to 0 if no dot file is available for this direction
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)  (NOT NEEDED !)
!!  eig_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*nsppol*mband**2)=matrix of first-order eigenvalues (hartree)
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
!!  nband_rbz=(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
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
!!  prtbbb=if 1, will compute the band-by-band decomposition
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rmet(3,3)=real space metric (bohr**2)
!!  typat(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  wffddk(3)=struct info for for the three possible dot files for ipert1
!!       equal to 0 if no dot file is available for this direction
!!  wffnow=struct info for INPUT 1st-order wf file
!!  wfftgs=struct info for GS wavefunction file at k
!!  wtk_k=weight assigned to the k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(npw_k,mpsang*mpsang*useylm)= real spherical harmonics for each
!!    G and k point
!!  ylm1(npw1_k,mpsang*mpsang*useylm)= real spherical harmonics for each
!!    G and k+q point
!!
!! OUTPUT
!!  d2bbb_k(2,3,mband,mband*prtbbb)=band by band decomposition of the second
!!   order derivatives, for the present k point, and perturbation idir, ipert
!!  d2nl_k(2,3,mpert)=non-local contributions to
!!   non-stationary 2DTE, for the present k point, and perturbation idir, ipert
!!
!! PARENTS
!!      nstdy3
!!
!! CHILDREN
!!      dotprod_g,gaugetransfo,leave_new,mkffnl,mkkpg,nonlop,ph1d3d,timab
!!      wffreaddatarec,wffreadnpwrec,wffreadskiprec,wrtout,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nstwf3(atindx,atindx1,cg,cg1,ddkfil,d2bbb_k,d2nl_k,&
&  ecut,ecutsm,&
&  eig_k,eig1_k,gmet,gprimd,icg,icg1,idir,ikpt,ipert,&
&  iscf,isppol,istwf_k,kg_k,kg1_k,kpt,mband,mgfft,mkmem,mk1mem,&
&  mpert,mpi_enreg,mpsang,mpw,mpw1,natom,nattyp,nband_k,nband_rbz,nfft,ngfft, &
&  nkpt_rbz,nline,nloalg,npw_k,npw1_k,nspinor,nsppol,ntypat,n4,n5,n6,occopt, &
&  occ_k,ortalg,ph1d,prtbbb,prtvol,psps,qphon,rmet,&
&  typat,ucvol,wffddk,wffnow,wfftgs,wtk_k,xred,ylm,ylm1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12spacepar
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_16response, except_this_one => nstwf3
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,icg1,idir,ikpt,ipert,iscf,isppol,istwf_k,mband,mgfft
 integer,intent(in) :: mk1mem,mkmem,mpert,mpsang,mpw,mpw1,n4,n5,n6,natom,nfft
 integer,intent(in) :: nkpt_rbz,nline,nsppol,ntypat,occopt,ortalg,prtbbb,prtvol
 integer,intent(inout) :: nband_k,npw1_k,npw_k,nspinor
 real(dp),intent(in) :: ecut,ecutsm,ucvol,wtk_k
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow,wfftgs
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),ddkfil(3),kg1_k(3,npw1_k)
 integer,intent(in) :: kg_k(3,npw_k),nattyp(ntypat),nband_rbz(nkpt_rbz*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(5),typat(natom)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: eig_k(mband*nsppol),gmet(3,3),gprimd(3,3),kpt(3)
 real(dp),intent(in) :: occ_k(nband_k),ph1d(2,3*(2*mgfft+1)*natom),qphon(3)
 real(dp),intent(in) :: rmet(3,3),xred(3,natom)
 real(dp),intent(in) :: ylm(npw_k,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(npw1_k,mpsang*mpsang*psps%useylm)
 real(dp),intent(inout) :: eig1_k(2*nsppol*mband**2)
 real(dp),intent(out) :: d2bbb_k(2,3,mband,mband*prtbbb),d2nl_k(2,3,mpert)
 type(wffile_type),intent(inout) :: wffddk(3)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimffnl,formeigdot,i1,i2,i3,ia,iatom,iband,ider,idir1
 integer :: ier,ierr,ig,ii,ilmn,iln,iln0,index,ipert1,iproj,ipsang,ipw,ipw1,it
 integer :: itypat,iwavef,jband,matblk,matblk_der,mcgnpw,mcgnpw1,me,n1,n2,n3
 integer :: natom_der,nband_disk,nband_kocc,nkpg,nkpg1,nnlout,npw_disk
 integer :: nspinor_disk,ntypat_der,paw_opt,shift1,shift2,shift3,signs
 integer :: tim_nonlop,tim_rwwf
 real(dp) :: aa,arg,dot1i,dot1r,dot2i,dot2r,dot_ndiagi,dot_ndiagr,doti,dotr,dum
 real(dp) :: im0,im1,re0,re1,scprod,tolwfr,valuei,valuer
 character(len=500) :: message
!arrays
 integer :: atindx_der(1),nattyp_der(1),nloalg_der(5),pspso_typ(1)
 integer,allocatable :: indlmn_typ(:,:,:),kg_dum(:,:)
 real(dp) :: dummy(2,1),enlout(3),kpq(3),nonlop_dum(1,1),phkxredin(2,1)
 real(dp) :: phkxredout(2,1),tsec(2),xred_der(3),ylmgr_dum(1)
 real(dp),allocatable :: cg_k(:,:),cgdot(:,:,:),cgnow(:,:),cgtgs(:,:)
 real(dp),allocatable :: cwave0(:,:),cwavef(:,:),cwavef_da(:,:),cwavef_db(:,:)
 real(dp),allocatable :: eig2_k(:),eig_dot(:,:),eig_dum(:),ekb_typ(:,:,:)
 real(dp),allocatable :: ffnl(:,:,:,:),ffnl1(:,:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: ffnlkq(:,:,:,:),ghc(:,:),grnk(:),gvnl1(:,:),gvnlc(:,:)
 real(dp),allocatable :: kinpw(:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),occ_dum(:)
 real(dp),allocatable :: ph1d_der(:,:),ph3d(:,:,:),ph3din(:,:,:),ph3dout(:,:,:)
 real(dp),allocatable :: phkxred(:,:)
 type(cprj_type) :: cprj_dum(1,1)
!no_abirules

! *********************************************************************

!Keep track of total time spent in nstwf3
 call timab(102,1,tsec)

!DEBUG
!write(6,*)' nstwf3 : enter '
!write(6,*)'ikpt,isppol,dimekb = ',ikpt,isppol,psps%dimekb
!write(6,*)' nstwf3 : npw_disk,nspinor_disk,nband_disk =',npw_disk,nspinor_disk,nband_disk
!write(6,*)' nstwf3 : second block,npw_disk,nspinor_disk,nband_disk =',npw_disk,nspinor_disk,nband_disk
!stop
!ENDDEBUG

!BEGIN TF_CHANGES

!Define me
 call xme_init(mpi_enreg,me)

!END TF_CHANGES

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 dimffnl=1
 allocate(ekb_typ(psps%dimekb,1,nspinor**2))
 allocate(indlmn_typ(6,psps%lmnmax,1))
 allocate(kinpw(npw_k),kinpw1(npw1_k))
 allocate(ffnl(npw_k,dimffnl,psps%lmnmax,ntypat))
 allocate(ffnlk(npw_k,dimffnl,psps%lmnmax,1))
 allocate(ffnl1(npw1_k,dimffnl,psps%lmnmax,ntypat))
 allocate(ffnlkq(npw1_k,dimffnl,psps%lmnmax,1))
 allocate(phkxred(2,natom))
 allocate(ghc(2,npw1_k*nspinor),gvnlc(2,npw1_k*nspinor))
 allocate(gvnl1(2,npw1_k*nspinor))
 allocate(ph1d_der(2,(2*n1+1)+(2*n2+1)+(2*n3+1)))
 allocate(eig2_k(2*nsppol*mband**2))

!Will compute the 3D phase factors inside nonlop
 allocate(ph3din(2,npw_k,1),ph3dout(2,npw1_k,1))
 nloalg_der(:)=nloalg(:)
 nloalg_der(1)=-abs(nloalg(1))
 nloalg_der(4)=1

!Compute (k+G) vectors
 nkpg=0;if (ipert/=natom+1) nkpg=3*nloalg(5)
 allocate(kpg_k(npw_k,nkpg))
 if (nkpg>0) call mkkpg(kg_k,kpg_k,kpt,nkpg,npw_k)

!Compute nonlocal form factors ffnl at (k+G), for all atoms
 ider=0
 call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gmet,gprimd,ider,ider,&
& psps%indlmn,kg_k,kpg_k,kpt,psps%lmnmax,&
& psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
& npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,ylmgr_dum)

 kpq(:)=kpt(:)+qphon(:)

!Compute (k+q+G) vectors
 nkpg1=0;if (ipert/=natom+1) nkpg1=3*nloalg(5)
 allocate(kpg1_k(npw1_k,nkpg1))
 if (nkpg1>0) call mkkpg(kg1_k,kpg1_k,kpq,nkpg1,npw1_k)

!Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
 ider=0
 call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl1,psps%ffspl,gmet,gprimd,ider,ider,&
& psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,&
& psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
& npw1_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1,ylmgr_dum)

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

!DEBUG
!write(6,*)' nstwf3 : before loop , call rdnpw'
!stop
!ENDDEBUG

!Take care of the npw and kg records
!NOTE : one should be able to modify the rwwf routine to take care
!of the band parallelism, which is not the case yet ...
 do idir1=1,3
  if (ddkfil(idir1)/=0)then
!  Read npw record
   call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_disk,nspinor,wffddk(idir1))
   if (npw_k /= npw_disk) then
    write(message,'(a,a,a,a,i3,a,i5,a,i3,a,a,i5,a,a,i5)')ch10,&
&    ' nstwf3: BUG - ',ch10,&
&    ' For isppol = ',isppol,', ikpt = ',ikpt,' and idir = ',idir,ch10,&
&    ' the number of plane waves in the ddk file is equal to', npw_disk,ch10,&
&    ' while it should be ',npw_k
    call wrtout(6,message,'PERS')
    call leave_new('PERS')
   end if
!  Skip k+G record
   call WffReadSkipRec(ierr,1,wffddk(idir1))
  end if
 end do
 if (mkmem==0) then
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_k,nspinor,wfftgs)
! Skip k+G and eigenvalue records in wfftgs
  call WffReadSkipRec(ierr,2,wfftgs)
 end if
!DEBUG
!write(6,*)' nstwf3 : nband_k, prtbbb=',nband_k,prtbbb
!write(6,*)' nstwf3 : second block,npw_disk,nspinor_disk,nband_disk =',npw_disk,nspinor_disk,nband_disk
!stop
!ENDDEBUG
 if (mk1mem==0) then
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw1_k,nspinor,wffnow)
! Skip k+G record
  call WffReadSkipRec(ierr,1,wffnow)
 end if

 if (ipert==natom+1) then
  nband_kocc = 0
  do iband = 1,nband_k
   if (abs(occ_k(iband)) > tol8) nband_kocc = nband_kocc + 1
  end do
 end if

 if(prtbbb==1)then
  allocate(cwavef_da(2,npw1_k*nspinor),cwavef_db(2,npw1_k*nspinor))
  allocate(cg_k(2,npw_k*nspinor*nband_k))
  if ((ipert == natom + 1).or.(ipert <= natom).or.(ipert == natom + 2)) then
   if (mkmem /= 0) then
    cg_k(:,:) = cg(:,1+icg:icg+nband_k*npw_k*nspinor)
   else
    do iband=1,nband_k
     call WffReadDataRec(cwave0,ierr,2*npw_k*nspinor,wfftgs)
     cg_k(:,(iband-1)*npw_k*nspinor+1:iband*npw_k*nspinor)=cwave0(:,:)
    end do
   end if
  end if
  d2bbb_k(:,:,:,:) = 0.0_dp
 end if

!Loop over bands
 do iband=1,nband_k

  if(mpi_enreg%paral_compil_kpt==1)then
!  BEGIN TF_CHANGES
   if(mpi_enreg%proc_distrb(ikpt,iband,isppol) /= me) then
!   END TF_CHANGES
    do idir1=1,3
!    Skip the eigenvalue and the wf records of this band
     if (ddkfil(idir1) /= 0) then
      call WffReadSkipRec(ierr,2,wffddk(idir1))
     end if
    end do
    if(mkmem==0)then
     if(prtbbb==0 .or. ipert==natom+2)then
      call WffReadSkipRec(ierr,1,wfftgs)
     end if
    end if
    if(mk1mem==0)then
     call WffReadSkipRec(ierr,2,wffnow)
    end if
    cycle
   end if
  end if ! mpi_enreg%paral_compil_kpt==1

! Read ground-state and first-order wavefunctions
  if (prtbbb==0 .or. ipert==natom+2) then
   if(mkmem/=0)then
    cwave0(:,:)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
   else
    call timab(286,1,tsec)
    call WffReadDataRec(cwave0,ierr,2*npw_k*nspinor,wfftgs)
    call timab(286,2,tsec)
   end if
  else    ! prtbbb==1 and ipert<=natom , already in cg_k
   cwave0(:,:)=cg_k(:,1+(iband-1)*npw_k*nspinor:iband*npw_k*nspinor)
  end if

  if(mk1mem/=0)then
   cwavef(:,:)=cg1(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)
  else
   call timab(286,1,tsec)
   call WffReadDataRec(eig1_k(1+(iband-1)*2*nband_k:2*iband*nband_k),ierr,2*nband_k,wffnow)
   call WffReadDataRec(cwavef,ierr,2*npw1_k*nspinor,wffnow)
   call timab(286,2,tsec)
  end if

! In case non ddk perturbation
  if (ipert /= natom + 1) then

   do ipert1=1,mpert
    if( ipert1<=natom .or. ipert1==natom+2 )then

     if( ipert1<=natom )then

      itypat=typat(ipert1)

      ekb_typ(:,1,1)=psps%ekb(:,itypat)
      if (nspinor==2) then
       ekb_typ(:,1,2)=psps%ekb(:,itypat)
       ekb_typ(:,1,3:4)=zero
      end if
      indlmn_typ(:,:,1)=psps%indlmn(:,:,itypat)
      pspso_typ(1)=psps%pspso(itypat)

!     Copy the part needed for the displaced atom, in ffnlkq.
      do ilmn=1,psps%lmnmax
       ffnlkq(:,:,ilmn,1)=ffnl1(:,:,ilmn,itypat)
       ffnlk (:,:,ilmn,1)=ffnl (:,:,ilmn,itypat)
      end do
     end if

     if (((ipert <= natom).or.(ipert == natom + 2)) &
&     .and.(ipert1 == natom+2).and. prtbbb==1) then
      call gaugetransfo(cg_k,cwavef,cwavef_db,eig_k,eig1_k,iband,nband_k, &
&      mband,npw_k,npw1_k,nspinor,nsppol,occ_k)
      cwavef(:,:) = cwavef_db(:,:)
     end if

!    Define the direction along which to move the atom :
!    the polarisation (ipert1,idir1) is refered as j1.
     do idir1=1,3
      if( ipert1<=natom .or. &
&      ( ipert1==natom+2 .and. ddkfil(idir1)/=0) ) then

!      Get |vnon-locj1|u0> :
!      first-order non-local, applied to zero-order wavefunction
!      (??) this routine gives MINUS the non-local contribution

       if( ipert1<=natom )then

        signs=2 ; choice=2 ; nnlout=3 ; natom_der=1 ; nattyp_der(1)=1
        ntypat_der=1 ; matblk_der=1
        xred_der(:)=xred(:,ipert1)
        atindx_der(1)=1

!       Store at the right place the 1d phases
        shift1=(atindx(ipert1)-1)*(2*n1+1)
        ph1d_der(:,1:2*n1+1)=ph1d(:,1+shift1:2*n1+1+shift1)
        shift2=(atindx(ipert1)-1)*(2*n2+1)+natom*(2*n1+1)
        ph1d_der(:,1+2*n1+1:2*n2+1+2*n1+1)=ph1d(:,1+shift2:2*n2+1+shift2)
        shift3=(atindx(ipert1)-1)*(2*n3+1)+natom*(2*n1+1+2*n2+1)
        ph1d_der(:,1+2*n1+1+2*n2+1:2*n3+1+2*n2+1+2*n1+1)=&
&        ph1d(:,1+shift3:2*n3+1+shift3)

!       Compute here phkxred for kpt and kpq
        arg=two_pi*(kpt(1)*xred_der(1)+kpt(2)*xred_der(2)+kpt(3)*xred_der(3))
        phkxredin(1,1)=cos(arg)
        phkxredin(2,1)=sin(arg)
        arg=two_pi*(kpq(1)*xred_der(1)+kpq(2)*xred_der(2)+kpq(3)*xred_der(3))
        phkxredout(1,1)=cos(arg)
        phkxredout(2,1)=sin(arg)

        tim_nonlop=5;paw_opt=0;cpopt=-1

        call nonlop(atindx1,choice,cpopt,cprj_dum,psps%dimekb,ntypat_der,dimffnl,&
&        dimffnl,ekb_typ,enlout,ffnlk,ffnlkq,gmet,gprimd,&
&        idir1,indlmn_typ,istwf_k,kg_k,kg1_k,kpg_k,kpg1_k,kpt,kpq,dum,psps%lmnmax,&
&        matblk_der,mgfft,mpi_enreg,psps%mpsang,psps%mpssoang,&
&        natom_der,nattyp_der,ngfft,nkpg,nkpg1,nloalg_der,nnlout,npw_k,npw1_k,nspinor,&
&        ntypat_der,0,paw_opt,phkxredin,phkxredout,ph1d_der,ph3din,ph3dout,pspso_typ,&
&        signs,nonlop_dum,nonlop_dum,&
&        tim_nonlop,ucvol,psps%useylm,cwave0,gvnl1)

       else if( ipert1==natom+2 )then

        call WffReadDataRec(eig2_k(1+(iband-1)*2*nband_k:2*iband*nband_k),ierr,2*nband_k,wffddk(idir1))
        call WffReadDataRec(gvnl1,ierr,2*npw1_k*nspinor,wffddk(idir1))

!       In case of band-by-band,
!       construct the first-order wavefunctions in the diagonal gauge
        if (((ipert <= natom).or.(ipert == natom + 2)).and.(prtbbb==1)) then
         call gaugetransfo(cg_k,gvnl1,cwavef_da,eig_k,eig2_k,iband,nband_k, &
&         mband,npw_k,npw1_k,nspinor,nsppol,occ_k)
         gvnl1(:,:) = cwavef_da(:,:)
        end if

!       Multiplication by -i
        do ipw=1,npw1_k*nspinor
         aa=gvnl1(1,ipw)
         gvnl1(1,ipw)=gvnl1(2,ipw)
         gvnl1(2,ipw)=-aa
        end do
       end if

!      MVeithen 021212 :
!      1) Case ipert1 = natom + 2 and ipert = natom + 2:
!      the second derivative of the energy with respect to an electric
!      field is computed from Eq. (38) of X. Gonze, PRB 55 ,10355 (1997).
!      The evaluation of this formula needs the operator $i \frac{d}{dk}.
!      2) Case ipert1 = natom + 2 and ipert < natom:
!      the computation of the Born effective charge tensor uses
!      the operator $-i \frac{d}{dk}.
       if (ipert==natom+2) gvnl1(:,:) = -gvnl1(:,:)

!      <G|Vnl1|Cnk> is contained in gvnl1
!      construct the matrix element (<uj2|vj1|u0>)complex conjug.
!      and add it to the 2nd-order matrix
!      XG030513 : use dotprod_g, for future parallelisation
       call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*nspinor,2,cwavef,gvnl1)

       d2nl_k(1,idir1,ipert1)= &
&       d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*2.0_dp*dotr
       d2nl_k(2,idir1,ipert1)=&
&       d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*2.0_dp*doti

!      Band by band decomposition of the Born effective charges
!      calculated from a phonon perturbation
       if(prtbbb==1)then
        d2bbb_k(1,idir1,iband,iband) = wtk_k*occ_k(iband)*2.0_dp*dotr
        d2bbb_k(2,idir1,iband,iband) = -1.0_dp*wtk_k*occ_k(iband)*2.0_dp*doti
       end if

!      DEBUG  Do not forget to restore the idir loop
!      write(6,*)' nstwf3 : cwave0 '
!      write(6,*)cwave0(:,1:2)
!      write(6,*)' nstwf3 : gvnl1 '
!      do ii=1,npw1_k*nspinor
!      write(6,*)ii,gvnl1(:,ii)
!      end do
!      write(6,*)' nstwf3 : cwavef '
!      do ii=1,npw1_k*nspinor
!      write(6,*)ii,cwavef(:,ii)
!      end do
!      if(idir1==3 .and. ipert1==4)then
!      write(6,*)' nstwf3 : ikpt,iband,ipert1,idir1,dotr,doti,d2nl_k'
!      write(6,*)ikpt, iband, ipert1, idir1
!      write(6,*)wtk_k,occ_k(iband),dotr, doti
!      write(6,*)d2nl_k(:,idir1,ipert1)
!      end if
!      stop
!      ENDDEBUG

      end if
     end do

    end if
   end do

  end if     ! ipert /= natom +1

! Compute the localization tensor

  if (ipert==natom+1) then

   ipert1=natom+1

   if(prtbbb==1)then
    call gaugetransfo(cg_k,cwavef,cwavef_db,eig_k,eig1_k,iband,nband_k, &
&    mband,npw_k,npw1_k,nspinor,nsppol,occ_k)
    cwavef(:,:) = cwavef_db(:,:)
   end if

   do idir1 = 1,3
    eig2_k(:) = 0.0_dp
    gvnl1(:,:) = 0.0_dp
    if (idir == idir1) then
     if (ddkfil(idir1) /= 0) then
      call WffReadSkipRec(ierr,2,wffddk(idir1))
     end if
     gvnl1(:,:) = cwavef(:,:)
     eig2_k(:) = eig1_k(:)
    else
     if (ddkfil(idir1) /= 0) then
      call WffReadDataRec(eig2_k(1+(iband-1)*2*nband_k:2*iband*nband_k),ierr,2*nband_k,wffddk(idir1))
      call WffReadDataRec(gvnl1,ierr,2*npw1_k*nspinor,wffddk(idir1))

      if(prtbbb==1)then
       call gaugetransfo(cg_k,gvnl1,cwavef_da,eig_k,eig2_k,iband,nband_k, &
&       mband,npw_k,npw1_k,nspinor,nsppol,occ_k)
       gvnl1(:,:) = cwavef_da(:,:)
      end if

     end if    !ddkfil(idir1)
    end if    !idir == idir1

!   <G|du/dqa> is contained in gvnl1 and <G|du/dqb> in cwavef
!   construct the matrix elements <du/dqa|du/dqb> -> dot
!   <u|du/dqa> -> dot1
!   <du/dqb|u> -> dot2
!   and add them to the 2nd-order matrix

!   XG030513 : use dotprod_g, for future parallelisation
    call dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw1_k*nspinor,2,gvnl1,cwavef)
    d2nl_k(1,idir1,ipert1)=d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*dotr/(nband_kocc*two)
    d2nl_k(2,idir1,ipert1)=d2nl_k(2,idir1,ipert1)+wtk_k*occ_k(iband)*doti/(nband_kocc*two)


!   XG 020216 : Marek, could you check the next forty lines
!   In the parallel gauge, dot1 and dot2 vanishes
    if(prtbbb==1)then
     d2bbb_k(1,idir1,iband,iband)=d2bbb_k(1,idir1,iband,iband)+dotr
     d2bbb_k(2,idir1,iband,iband)=d2bbb_k(2,idir1,iband,iband)+doti
     dot_ndiagr=0.0_dp ; dot_ndiagi=0.0_dp
     do jband = 1,nband_k              !compute dot1 and dot2
      if (abs(occ_k(jband)) > tol8) then

       dot1r=0.0_dp ; dot1i=0.0_dp
       dot2r=0.0_dp ; dot2i=0.0_dp

       cwave0(:,:)=cg_k(:,1+(jband-1)*npw_k*nspinor:jband*npw_k*nspinor)
!      XG030513 : use dotprod_g, for future parallelisation
       call dotprod_g(dot1r,dot1i,istwf_k,mpi_enreg,npw1_k*nspinor,2,cwave0,gvnl1)
       call dotprod_g(dot2r,dot2i,istwf_k,mpi_enreg,npw1_k*nspinor,2,cwavef,cwave0)

       dot_ndiagr = dot_ndiagr + dot1r*dot2r - dot1i*dot2i
       dot_ndiagi = dot_ndiagi + dot1r*dot2i + dot1i*dot2r

       d2bbb_k(1,idir1,iband,jband) = d2bbb_k(1,idir1,iband,jband) - &
&       (dot1r*dot2r - dot1i*dot2i)
       d2bbb_k(2,idir1,iband,jband) = d2bbb_k(2,idir1,iband,jband) - &
&       (dot1r*dot2i + dot1i*dot2r)

      end if  ! occ_k
     end do !jband

     d2bbb_k(:,idir1,iband,:)= &
&     d2bbb_k(:,idir1,iband,:)*wtk_k*occ_k(iband)*half

     d2nl_k(1,idir1,ipert1)= &
&     d2nl_k(1,idir1,ipert1)-wtk_k*occ_k(iband)*dot_ndiagr/(nband_kocc*two)
     d2nl_k(2,idir1,ipert1)=&
&     d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*dot_ndiagi/(nband_kocc*two)

    end if ! prtbbb==1

   end do  ! idir1
  end if   ! Compute localization tensor, ipert=natom+1

! End loop over bands
 end do

!DEBUG
!write(6,*)' nstwf3 : after loop over bands '
!write(6,*)' nstwf3 : second block,npw_disk,nspinor_disk,nband_disk =',npw_disk,nspinor_disk,nband_disk
!stop
!ENDDEBUG


 deallocate(cwave0,cwavef)

!###################################################################

 deallocate(eig2_k,ghc,gvnlc,gvnl1)
 deallocate(ffnl,ffnl1,kinpw,kinpw1,kpg_k,kpg1_k,phkxred,ph3d)
 deallocate(ffnlkq,ffnlk)
 deallocate(ph1d_der,ph3din,ph3dout)
 deallocate(ekb_typ,indlmn_typ)
 if(prtbbb==1)deallocate(cg_k,cwavef_da,cwavef_db)

 call timab(102,2,tsec)

!DEBUG
!write(6,*)' nstwf3 : exit '
!write(6,*)' nstwf3 : d2nl_k(2,3,2:4:2)',d2nl_k(2,3,2:4:2)
!stop
!ENDDEBUG

end subroutine nstwf3
!!***
