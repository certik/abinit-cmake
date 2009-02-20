!{\src2tex{textfont=tt}}
!!****f* ABINIT/nselt3
!! NAME
!! nselt3
!!
!! FUNCTION
!! This routine compute the non-stationary expression for the
!! second derivative of the total energy, wrt strain for a whole row of
!! mixed strain derivatives.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DRH, XG, DCA, GMR, MM, AR, MV, MB)
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
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass=effective mass for electrons (1. in common case)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues
!!    (hartree)
!!  fform=index for choosing form of wf file
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of the perturbation
!!  indkpt1(nkpt_rbz)=non-symmetrized indices of the k-points
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  iscf=(<= 0  =>non-SCF), >0 => SCF
!!  istep=index of the number of steps in the routine scfcv
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the
!!     storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points in the reduced BZ
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpert =maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  maximum dimension for q points in grids for nonlocal form factors
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nband(nkpt*nsppol)=number of bands at each k point for each spin
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points in the full BZ
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array. If /=0,
!!   the exchange-correlation kernel must be computed.
!!  nline=number of CG line minimizations per band.
!!      for the time being, also number of non-SCF passes through all bands
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!    perturbation
!!  ntypat=number of types of atoms in unit cell.
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occopt=option for occupancies
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band
!!   and k in the reduced Brillouin zone (usually =2)
!!  ortalg=governs the choice of the algorithm for orthogonalisation.
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  prtbbb=if 1, band-by-band decomposition (also dim of d2bbb)
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rhor(nfft,nspden)=GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  type(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  unkg=unit number of k+G data file
!!  unkg1=unit number of k+G+q data file
!!  wffnow= struct info for wf disk file
!!  wfftgs=struct info r for ground-state wf disk file
!!  unylm=unit number for disk file containing Ylm(k) if mkmem==0
!!  unylm1=unit number for disk file containing Ylm(k+q) if mk1mem==0
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang)= real spherical harmonics for each
!!    G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang)= real spherical harmonics for each
!!    G and k+q point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical for each
!!    G and k point
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm)= gradients of real spherical for each
!!    G and k+g point
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!!  d2lo(2,3,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,3,mpert,3,mpert)=non-local contributions to the 2DTEs
!! Not used (should be suppressed, later)
!!  rhog1(2,nfft)=RF electron density in reciprocal space
!!  tsmear=smearing energy or temperature (if metal)
!!  fnamewffddk
!!
!! NOTES
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      clsopn,dotprod_vn,hartrestr,hdr_skip,leave_test,mkcor3,mkvxcstr3,nstwf4
!!      rdnpw,stresssym,vlocalstr,wrtout,xcomm_init,xdefineoff,xmaster_init
!!      xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nselt3(atindx,atindx1,blkflg,cg,cg1,cplex,doccde_rbz,docckqde,&
& d2bbb,d2lo,d2nl,ecut,ecutsm,effmass,eigen0,eigen1,fform,&
& gmet,gprimd,gsqcut,idir,indkpt1,indsy1,&
& ipert,iscf,istep,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mband,mgfft,&
& mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
& natom,nattyp,nband,nband_rbz,nfft,ngfft,&
& nkpt,nkpt_rbz,nkxc,nline,nloalg,npwarr,npwar1,nspden,nspinor,nsppol,&
& nsym1,ntypat,occkq,occopt,occ_rbz,ortalg,&
& paral_kgb, ph1d,prtbbb,prtvol,psps,qphon,rhog,rhog1,&
& rhor,rhor1,rmet,rprimd,symrc1,tsmear,type,ucvol,&
& unkg,unkg1,&
& wffnow,wfftgs,unylm,unylm1,fnamewffddk,wtk_rbz,&
& xred,ylm,ylm1,ylmgr,ylmgr1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
 use interfaces_12spacepar
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
 use interfaces_16response, except_this_one => nselt3
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,fform,idir,ipert,iscf,istep,mband,mgfft,mk1mem
 integer,intent(in) :: mkmem,mpert,mpsang,mpw,mpw1,natom,nfft,nkpt,nkpt_rbz
 integer,intent(in) :: nkxc,nline,nspden,nsppol,nsym1,ntypat,occopt,ortalg
 integer,intent(in) :: paral_kgb,prtbbb,prtvol,unkg,unkg1,unylm,unylm1
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ecut,ecutsm,effmass,gsqcut,tsmear,ucvol
 character(len=fnlen),intent(in) :: fnamewffddk
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow,wfftgs
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),indkpt1(nkpt_rbz)
 integer,intent(in) :: indsy1(4,nsym1,natom),istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem),nattyp(ntypat)
 integer,intent(in) :: nband(nkpt*nsppol),nband_rbz(nkpt_rbz*nsppol),ngfft(18)
 integer,intent(in) :: nloalg(5),npwar1(nkpt_rbz),npwarr(nkpt_rbz)
 integer,intent(in) :: symrc1(3,3,nsym1),type(natom)
 integer,intent(out) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)
 real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen1(2*mband*mband*nkpt_rbz*nsppol),gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpt_rbz(3,nkpt_rbz),kxc(nfft,nkxc)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: occkq(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qphon(3),rhog(2,nfft)
 real(dp),intent(in) :: rhog1(2,nfft),rhor(nfft,nspden)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang)
 real(dp),intent(out) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: ban2tot,bantot,bd2tot_index,bdtot_index,ddkcase,enunit,formeig
 integer :: iband,iband1,icg,icg1,idir1,ierr,iexit,ifft,ii,ikg,ikg1,ikpt
 integer :: ikpt_dum,ilm,ipert1,iproc,ipw,ir,ispden,isppol,istr1,istwf_k,isym
 integer :: jj,master,mbd2kpsp,mbdkpsp,me,muig,n1,n2,n3,n3xccc,n4,n5,n6
 integer :: nband_dum,nband_k,nfftot,npw1_k,npw_k,nskip,option,spaceComm
 integer :: t_iostat,tim_rwwf
 real(dp) :: doti,dotr,dum,rho1_dn,rho1_up,rho1im_dn,rho1im_up,rho1re_dn
 real(dp) :: rho1re_up,wtk_k
 logical :: logi,t_exist
 character(len=500) :: message
 character(len=fnlen) :: fil,fiwfddk
!arrays
 integer :: ikpt_fbz(3),ikpt_fbz_previous(3)
 integer,allocatable :: kg1_k(:,:),kg_k(:,:),symrl1(:,:,:)
 real(dp) :: d2nl_elfd(2,3),imstr(6),kpoint(3),restr(6),sumelfd(2),tsec(2)
 real(dp),allocatable :: cg_dum(:,:),d2bbb_k(:,:,:,:),d2nl_k(:,:,:),eig1_k(:)
 real(dp),allocatable :: eig_k(:),eigen(:),eigen_dum(:),occ_dum(:),occ_k(:)
 real(dp),allocatable :: vhartr01(:),vpsp1(:),vxc1(:,:),xccc3d1(:),ylm1_k(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr1_k(:,:,:),ylmgr_k(:,:,:)

! *********************************************************************

!DEBUG
!write(6,*)' nselt3 : enter '
!stop
!ENDDEBUG

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)
!Init master
 call xmaster_init(mpi_enreg,master)

!Unit numbers

!Zero only portion of nonlocal matrix to be computed here
 d2nl(:,:,natom+3:natom+4,idir,ipert)=0.0_dp
 bdtot_index=0
 bd2tot_index=0
 icg=0
 icg1=0
 mbdkpsp=mband*nkpt_rbz*nsppol
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol

!Update list of computed matrix elements
 if((ipert==natom+3) .or. (ipert==natom+4)) then
! Eventually expand when strain coupling to other perturbations is
! implemented
  do ipert1=natom+3,natom+4
   do idir1=1,3
    blkflg(idir1,ipert1,idir,ipert)=1
   end do
  end do
 end if

!allocate(enl1nk(mbdkpsp))
 allocate(d2bbb_k(2,3,mband,mband*prtbbb))
 allocate(d2nl_k(2,3,mpert))

 allocate(eig_k(nsppol*mband))
 allocate(eig1_k(2*nsppol*mband**2))

 allocate(kg_k(3,mpw))
 allocate(kg1_k(3,mpw1))


 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
 nfftot=n1*n2*n3

!Prepare GS k wf file for reading if mkmem==0
 if (mkmem==0) then
  call clsopn(wfftgs)
  call hdr_skip(wfftgs,ierr)

! Define offsets, in case of MPI I/O
  formeig=0
  call xdefineOff(formeig,wfftgs,mpi_enreg,nband_rbz,npwarr,nspinor,nsppol,nkpt_rbz)

 end if

!Prepare RF wf files for reading and writing if mkmem==0
 if (mk1mem==0) then

  call clsopn(wffnow)

! Read unwfnow header
  call hdr_skip(wffnow,ierr)

! Define offsets, in case of MPI I/O
  formeig=1
  call xdefineOff(formeig,wffnow,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)

 end if

 bantot = 0
 ban2tot = 0

!LOOP OVER SPINS
 do isppol=1,nsppol

  if (nsppol/=1) then
   write(message,*)' ****  In nselt3 for isppol=',isppol
   call wrtout(06,message,'COLL')
  end if

! Rewind kpgsph data file if needed:
  if (mkmem==0) rewind(unkg)
  if (mk1mem==0) rewind(unkg1)
  if (mkmem==0.and.psps%useylm==1) rewind(unylm)
  if (mk1mem==0.and.psps%useylm==1) rewind(unylm1)
  ikg=0
  ikg1=0

  ikpt_fbz(1:3)=0

! BIG FAT k POINT LOOP
  do ikpt=1,nkpt_rbz

   nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
   istwf_k=istwfk_rbz(ikpt)
   npw_k=npwarr(ikpt)
   npw1_k=npwar1(ikpt)

   eig_k(1:nband_k) = eigen0(1+bantot:nband_k+bantot)
   eig1_k(1:2*nband_k**2) = eigen1(1+ban2tot:2*nband_k**2+ban2tot)
   bantot = bantot + nband_k
   ban2tot = ban2tot + 2*nband_k**2

   if(mpi_enreg%paral_compil_kpt==1)then
!   BEGIN TF_CHANGES
    if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&    -me))/=0) then
!    END TF_CHANGES
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2
!    Skip the rest of the k-point loop
     cycle
    end if
   end if

   allocate(occ_k(nband_k))

   allocate(ylm_k(npw_k,mpsang*mpsang))
   allocate(ylm1_k(npw1_k,mpsang*mpsang))
   if (ipert==natom+3.or.ipert==natom+4) then
    allocate(ylmgr_k(npw_k,3,mpsang*mpsang))
    allocate(ylmgr1_k(npw1_k,3,mpsang*mpsang))
   end if

!  enl1_k(:)=0.0_dp
   d2nl_k(:,:,:)=0.0_dp
   if(prtbbb==1)d2bbb_k(:,:,:,:)=0.0_dp
   kpoint(:)=kpt_rbz(:,ikpt)
   occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)

   if (mkmem==0) then
!   Read (k+G) basis sphere data (same for each spin)
    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,unkg)

!   Read sphere data centered at k in unkg, then k+g data
    read (unkg) ((kg_k(ii,muig),ii=1,3),muig=1,npw_k)

!   Eventually read (k+G) spherical harmonics
    if (psps%useylm==1) then
     read(unylm)
     if (ipert==natom+3.or.ipert==natom+4) then
      read(unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang),&
&      (((ylmgr_k(muig,ii,ilm),muig=1,npw_k),ii=1,3),ilm=1,mpsang*mpsang)
     else
      read(unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
     end if
    end if

   else

    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
      ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     end do
     if (ipert==natom+3.or.ipert==natom+4) then
      do ilm=1,mpsang*mpsang
       do ii=1,3
        ylmgr_k(1:npw_k,ii,ilm)=ylmgr(1+ikg:npw_k+ikg,ii,ilm)
       end do
      end do
     end if
    end if

!   End if for choice governed by mkmem
   end if

   wtk_k=wtk_rbz(ikpt)

   kg1_k(:,:) = 0
   if (mk1mem==0) then
!   Read (k+q+G) basis sphere data (same for each spin)
    call rdnpw(ikpt,isppol,nband_k,npw1_k,nspinor,0,unkg1)

!   Read sphere data centered at k in unkg, then k+g data
    read (unkg1) ((kg1_k(ii,muig),ii=1,3),muig=1,npw1_k)

!   Eventually read (k+q+G) spherical harmonics
    if (psps%useylm==1) then
     read(unylm1)
     if (ipert==natom+3.or.ipert==natom+4) then
      read(unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang),&
&      (((ylmgr1_k(muig,ii,ilm),muig=1,npw1_k),ii=1,3),ilm=1,mpsang*mpsang)
     else
      read(unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang)
     end if
    end if

   else

    kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
      ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
     end do
     if (ipert==natom+3.or.ipert==natom+4) then
      do ilm=1,mpsang*mpsang
       do ii=1,3
        ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
       end do
      end do
     end if
    end if

!   End if for choice governed by mk1mem
   end if

!  Compute the eigenvalues, wavefunction,
!  contributions to kinetic energy, nonlocal energy, forces,
!  and update of rhor1 to this k-point and this spin polarization.

!  Note that nstwf4 is called with kpoint, while kpt is used inside vtowfk3
   call nstwf4(atindx,atindx1,cg,cg1,d2nl_k,ecut,ecutsm,effmass,&
&   eig_k,eig1_k,gmet,gprimd,icg,icg1,idir,ikpt,ipert,&
&   iscf,isppol,istwf_k,kg_k,kg1_k,kpoint,mband,mgfft,mkmem,mk1mem,&
&   mpert,mpi_enreg,mpsang,mpw,mpw1,natom,nattyp,&
&   nband_k,nfft,ngfft,&
&   nline,nloalg,npw_k,npw1_k,nspinor,nsppol,ntypat,n4,n5,n6,occopt, &
&   occ_k,ortalg,ph1d,prtvol,psps,qphon,rmet,&
&   type,ucvol,wffnow,wfftgs,&
&   wtk_k,xred,ylm_k,ylm1_k,ylmgr_k,ylmgr1_k)
   d2nl(:,:,:,idir,ipert)=d2nl(:,:,:,idir,ipert)+d2nl_k(:,:,:)
   if(prtbbb==1)then
    d2bbb(:,:,idir,ipert,:,:) = d2bbb(:,:,idir,ipert,:,:) + &
&    d2bbb_k(:,:,:,:)
   end if

   deallocate(occ_k)

!  Keep track of total number of bands (all k points so far, even for
!  k points not treated by me)
   bdtot_index=bdtot_index+nband_k
   bd2tot_index=bd2tot_index+2*nband_k**2

!  Shift array memory
   if (mkmem/=0) then
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   end if
   if (mk1mem/=0) then
    icg1=icg1+npw1_k*nspinor*nband_k
    ikg1=ikg1+npw1_k
   end if
   deallocate(ylm_k,ylm1_k)
   if (ipert==natom+3.or.ipert==natom+4) deallocate(ylmgr_k,ylmgr1_k)

!  End big k point loop
  end do

! End loop over spins
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
! BEGIN TF_CHANGES
  call leave_test(mpi_enreg)
! END TF_CHANGES
  write(message,*) ' nselt3: loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
 end if

!Treat now varying occupation numbers
!if(occopt>=3 .and. occopt <=7) then
!SUPPRESSED metallic coding of vtorho

!Treat fixed occupation numbers
!else

!Accumulation over parallel processed now carried out for all terms
!in nstdy3.f

!End of test on varying or fixed occupation numbers
!end if

!The imaginary part of d2nl will be must be set to zero here since
!time-reversal symmetry will always be true for the strain peturbation.
!The symmetry-reduced kpt set will leave a non-zero imaginary part.

 d2nl(2,:,natom+3:natom+4,idir,ipert)=0.0_dp

!Symmetrize the non-local contributions,
!as was needed for the stresses in a ground-state calculation

 if (nsym1>1) then
! Pack like symmetric-storage cartesian stress tensor
  ii=0
  do ipert1=natom+3,natom+4
   do idir1=1,3
    ii=ii+1
    restr(ii)=d2nl(1,idir1,ipert1,idir,ipert)
   end do
  end do
! Do the symmetrization using the ground state routine
  call stresssym(gprimd,nsym1,restr,symrc1)
! Unpack symmetrized stress tensor
  ii=0
  do ipert1=natom+3,natom+4
   do idir1=1,3
    ii=ii+1
    d2nl(1,idir1,ipert1,idir,ipert)=restr(ii)
   end do
  end do
 end if !nsym>1

!----------------------------------------------------------------------------
!Now, treat the local contribution

 allocate(vpsp1(cplex*nfft))
 n3xccc=0
 if(psps%n1xccc/=0)n3xccc=cplex*nfft
 allocate(xccc3d1(cplex*nfft))
 allocate(vxc1(cplex*nfft,nspden))
 allocate(vhartr01(nfft))

 xccc3d1(:)=0.0_dp

!Double loop over strain perturbations
 do ipert1=natom+3,natom+4
  do idir1=1,3
   if(ipert1==natom+3) then
    istr1=idir1
   else
    istr1=idir1+3
   end if

!  Get first-order local potential.
   call vlocalstr(gmet,gprimd,gsqcut,istr1,mgfft,mpi_enreg,&
&   psps%mqgrid_vl,natom,nattyp,nfft,ngfft,ntypat,paral_kgb,ph1d,psps%qgrid_vl,&
&   ucvol,psps%vlspl,vpsp1)

!  Get first-order hartree potential.
   call hartrestr(gmet,gprimd,gsqcut,idir1,ipert1,mpi_enreg,natom,nfft,ngfft,&
&   paral_kgb,rhog,vhartr01)

!  Get first-order exchange-correlation potential
   if(psps%n1xccc/=0)then
    call mkcor3(cplex,idir1,ipert1,natom,ntypat,n1,psps%n1xccc,&
&    n2,n3,qphon,rprimd,type,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
   end if ! psps%n1xccc/=0

   option=0
   call mkvxcstr3(cplex,gmet,gsqcut,idir1,ipert1,kxc,mpi_enreg,natom,nfft,ngfft,&
&   nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,vxc1,xccc3d1)

!  Combines density j2 with local potential j1
   do ispden=1,min(nspden,2)
    do ifft=1,cplex*nfft
     vxc1(ifft,ispden)=vxc1(ifft,ispden)+vpsp1(ifft)+vhartr01(ifft)
    end do
   end do
   call dotprod_vn(cplex,rhor1,dotr,doti,mpi_enreg,nfft,nfftot,nspden,2,vxc1,ucvol)
   write(6,*)
   d2lo(1,idir1,ipert1,idir,ipert)=dotr
   d2lo(2,idir1,ipert1,idir,ipert)=doti
  end do ! istr1
 end do ! ipert1

 deallocate(vxc1,xccc3d1)
 deallocate(vhartr01) ! inserted by MM

 deallocate(d2bbb_k,d2nl_k,kg_k,kg1_k,vpsp1)
 deallocate(eig_k,eig1_k)

!DEBUG
!write(6,*)' nselt3: exit '
!ENDDEBUG

end subroutine nselt3
!!***
