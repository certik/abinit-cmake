!{\src2tex{textfont=tt}}
!!****f* ABINIT/nstdy3
!! NAME
!! nstdy3
!!
!! FUNCTION
!! This routine compute the non-stationary expression for the
!! second derivative of the total energy, for a whole row of
!! mixed derivatives.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, DCA, GMR, MM, AR, MV, MB)
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
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
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
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage
!!    of wfs
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
!!   see ~abinit/doc/input_variables/vargs.htm#ngfft
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
!!  ntypat=number of types of atoms in unit cell.
!!  nsym1=number of symmetry elements in space group consistent with i
!!    perturbation
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
!!  rhor1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  typat(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  unkg=unit number of k+G data file
!!  unkg1=unit number of k+G+q data file
!!  wffnow=struct info for wf disk file
!!  wfftgs=struct info for ground-state wf disk file
!!  wfnameddk=eneric name of the ddk response wavefunction file(s)
!!  unylm=unit number for disk file containing Ylm(k) if mkmem==0
!!  unylm1=unit number for disk file containing Ylm(k+q) if mk1mem==0
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for
!!    each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for
!!    each G and k+q point
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
!!
!! NOTES
!! Note that the ddk perturbation should not be treated here.
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      appdig,clsopn,dotprod_vn,hdr_skip,leave_new,leave_test,mati3inv,mkcor3
!!      mkvxc3,nstwf3,rdnpw,sygra3,timab,vloca3,wffclose,wffkg,wffopen
!!      wffreadnpwrec,wffreadskipk,wffreadskiprec,wrtout,xcomm_world,xdefineoff
!!      xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nstdy3(atindx,atindx1,blkflg,cg,cg1,cplex,doccde_rbz,docckqde,&
& dtset,d2bbb,d2lo,d2nl,ecut,ecutsm,eigen0,eigen1,fform,&
& gmet,gprimd,gsqcut,idir,indkpt1,indsy1,&
& ipert,iscf,istep,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mband,mgfft,&
& mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
& natom,nattyp,nband,nband_rbz,nfft,ngfft,&
& nkpt,nkpt_rbz,nkxc,nline,nloalg,npwarr,npwar1,nspden,nspinor,nsppol,&
& nsym1,ntypat,occkq,occopt,occ_rbz,ortalg,&
& paral_kgb,ph1d,prtbbb,prtvol,psps,qphon,rhog1,&
& rhor1,rmet,rprimd,symrc1,tsmear,typat,ucvol,&
& unkg,unkg1,wffnow,wfftgs,unylm,unylm1,wfnameddk,wtk_rbz,&
& xred,ylm,ylm1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12spacepar
 use interfaces_13io_mpi
 use interfaces_13xc
 use interfaces_14iowfdenpot
 use interfaces_16response, except_this_one => nstdy3
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
 real(dp),intent(in) :: ecut,ecutsm,gsqcut,tsmear,ucvol
 character(len=fnlen),intent(in) :: wfnameddk
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow,wfftgs
!arrays
 integer,intent(in) :: atindx(natom),atindx1(natom),indkpt1(nkpt_rbz)
 integer,intent(in) :: indsy1(4,nsym1,natom),istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem),nattyp(ntypat)
 integer,intent(in) :: nband(nkpt*nsppol),nband_rbz(nkpt_rbz*nsppol),ngfft(18)
 integer,intent(in) :: nloalg(5),npwar1(nkpt_rbz),npwarr(nkpt_rbz)
 integer,intent(in) :: symrc1(3,3,nsym1),typat(natom)
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
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qphon(3),rhog1(2,nfft)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
 real(dp),intent(out) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: ban2tot,bantot,bd2tot_index,bdtot_index,ddkcase,enunit,formeig
 integer :: iband,iband1,icg,icg1,idir1,ierr,iexit,ifft,ii,ikg,ikg1,ikpt
 integer :: ikpt_dum,ilm,index,index1,index2,ipert1,iproc,ipw,ir,ispden,isppol
 integer :: istwf_k,isym,jj,master,mbd2kpsp,mbdkpsp,mcgnpw,mcgnpw1,me,muig,n1
 integer :: n2,n3,n3xccc,n4,n5,n6,nband_dum,nband_k,nfftot,npw1_k,npw_k,nskip
 integer :: nspinor_,option,spaceworld,t_iostat,tag,tim_rwwf
 real(dp) :: di_psp1,di_xc1,doti,dotr,dr_psp1,dr_xc1,dum,wtk_k
 logical :: logi,t_exist
 character(len=500) :: message
 character(len=fnlen) :: fil,fiwfddk
!arrays
 integer :: ddkfil(3),ikpt_fbz(3),ikpt_fbz_previous(3),skipddk(3)
 integer,allocatable :: kg1_k(:,:),kg_dum(:,:),kg_k(:,:),symrl1(:,:,:)
 real(dp) :: d2nl_elfd(2,3),kpoint(3),sumelfd(2),tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:),cg_dum(:,:),d2bbb_k(:,:,:,:)
 real(dp),allocatable :: d2nl_k(:,:,:),eig1_k(:),eig_k(:),eigen(:),eigen_dum(:)
 real(dp),allocatable :: occ_dum(:),occ_k(:),rhodummy(:,:),vpsp1(:),vxc1(:,:)
 real(dp),allocatable :: work1(:,:,:),xccc3d1(:),ylm1_k(:,:),ylm_k(:,:)
 type(wffile_type) :: wffddk(3)
!no_abirules
#if defined T3E
           integer,save :: file_exist(3)
           logical,save :: first=.true.
           character(len=fnlen),save :: wfnameddk_old
#endif

! *********************************************************************

!DEBUG
!write(6,*)' nstdy3 : enter, debug '
!stop
!write(6,*)' nstdy3 : cg1(:,1)=',cg1(:,1)
!write(6,*)' xred=',xred
!ENDDEBUG

!Keep track of total time spent in nstdy3
 call timab(101,1,tsec)

!Init spaceworld
!BEGIN TF_CHANGES
 call xcomm_world(mpi_enreg,spaceworld)
!END TF_CHANGES
 master =0
!Define me
 call xme_init(mpi_enreg,me)

!Zero only portion of nonlocal matrix to be computed here
 d2nl(:,:,1:natom+2,idir,ipert)=0.0_dp
 bdtot_index=0
 bd2tot_index=0
 icg=0
 icg1=0
 mbdkpsp=mband*nkpt_rbz*nsppol
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol

!allocate(enl1nk(mbdkpsp))
 allocate(d2bbb_k(2,3,mband,mband*prtbbb))
 allocate(d2nl_k(2,3,mpert))

 allocate(eig_k(nsppol*mband))
 allocate(eig1_k(2*nsppol*mband**2))

 allocate(kg_k(3,mpw))
 allocate(kg1_k(3,mpw1))

!enl1nk(:)=0.0_dp

!Do not try to open electric field file
 ddkfil(:)=0
!The treatment of homogeneous electric field potential need
!the existence of d/dk files.
 do idir1=1,3
  ddkcase=idir1+natom*3
  call appdig(ddkcase,wfnameddk,fiwfddk)
! DEBUG
! write(6,*)' nstdy3 : examine existence of ddkfile',fiwfddk
! ENDDEBUG
! Check that ddk file exists
#if defined T3E
  if (first) then
!  case first time : test file to do
#endif
   inquire(file=fiwfddk,iostat=t_iostat,exist=t_exist)
   if (t_iostat.ne.0) then
    write(message, '(8a,i8)' )ch10,&
&    ' nstdy3 : BUG -',ch10,&
&    '  Check for existence of file ',trim(fiwfddk),',',ch10,&
&    '  but INQUIRE statement returns error code',t_iostat
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
#if defined T3E
   end if
   if (t_exist) then
    file_exist(idir1)=1
   else
    file_exist(idir1)=0
   end if
   if (idir1 == 3) then
    first=.false.
    wfnameddk_old = wfnameddk
   end if
  else
   if (wfnameddk_old /= wfnameddk) then
!   ddkfile changes
    inquire(file=fiwfddk,iostat=t_iostat,exist=t_exist)
    if (t_iostat.ne.0) then
     write(message, '(8a,i8)' )ch10,&
&     ' nstdy3 : BUG -',ch10,&
&     '  Check for existence of file ',trim(fiwfddk),',',ch10,&
&     '  but INQUIRE statement returns error code',t_iostat
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
    if (t_exist) then
     file_exist(idir1)=1
    else
     file_exist(idir1)=0
    end if
    if (idir1 == 3) wfnameddk_old = wfnameddk
   end if
  end if
! If the file exists, open it (even if it is not needed ...)
  if (file_exist(idir1)==1) then
#else
  else if (t_exist) then
#endif
!  Note the use of unit numbers 21, 22 and 23
   ddkfil(idir1)=20+idir1
   write(message, '(a,a)') '-open ddk wf file :',fiwfddk
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   call WffOpen(dtset%accesswff,spaceworld,fiwfddk,ierr,wffddk(idir1),master,me,ddkfil(idir1))
  end if
 end do

!DEBUG
!write(6,*)' nstdy3 : stop now '
!stop
!ENDDEBUG

!Update list of computed matrix elements
 if (ipert /= natom + 1) then
  do ipert1=1,mpert
   do idir1=1,3
    if(ipert1 <= natom .or.  &
&    (ipert1==natom+2 .and. ddkfil(idir1)/=0) )then
     blkflg(idir1,ipert1,idir,ipert)=1
    end if
   end do
  end do
 else
  ipert1 = natom + 1
  do idir1=1,3
   if ((ddkfil(idir1) /= 0).or.(idir1 == idir)) then
    blkflg(idir1,ipert1,idir,ipert)=1
   end if
  end do
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
 nfftot=n1*n2*n3

!Prepare GS k wf file for reading if mkmem==0
 if (mkmem==0) then

  call clsopn(wfftgs)
  call hdr_skip(wfftgs,ierr)

! Define offsets, in case of MPI I/O
  formeig=0
  call WffKg(wfftgs,1)
  call xdefineOff(formeig,wfftgs,mpi_enreg,nband_rbz,npwarr,nspinor,nsppol,nkpt_rbz)

 end if

!Prepare RF wf files for reading and writing if mkmem==0
 if (mk1mem==0) then

  call clsopn(wffnow)
  call hdr_skip(wffnow,ierr)

! Define offsets, in case of MPI I/O
  formeig=1
  call WffKg(wffnow,1)
  call xdefineOff(formeig,wffnow,mpi_enreg,nband_rbz,npwar1,nspinor,nsppol,nkpt_rbz)

 end if

!Initialisation of the ddk files
 do idir1=1,3
  if (ddkfil(idir1)/=0)then
   call hdr_skip(wffddk(idir1),ierr)

  end if
 end do

 bantot = 0
 ban2tot = 0
 skipddk(:) = 0

!LOOP OVER SPINS
 do isppol=1,nsppol

  if (nsppol/=1) then
   write(message,*)' ****  In nstdy3 for isppol=',isppol
   call wrtout(06,message,'COLL')
  end if

! In case isppol = 2, skip the records that correspond to isppol = 1
! and that have not been read
  if (isppol == 2) then
   do idir1 = 1, 3
    if ((ddkfil(idir1)/=0).and.(skipddk(idir1) < nkpt)) then
     do ikpt = 1, (nkpt - skipddk(idir1))
      call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_k,nspinor_,wffddk(idir1))
      call WffReadSkipRec(ierr,1,wffddk(idir1))
      do iband = 1, nband_k
       call WffReadSkipRec(ierr,2,wffddk(idir1))
      end do
     end do
    end if
   end do
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
!    The wavefunction blocks for ddk file is skipped elsewhere in the loop
!    Skip the rest of the k-point loop
     cycle
    end if
   end if

   allocate(ylm_k(npw_k,mpsang*mpsang*psps%useylm))
   allocate(ylm1_k(npw1_k,mpsang*mpsang*psps%useylm))

!  In case of electric field pert1, read ddk wfs file
!  Note that the symmetries are not used for ddk, so read each k point
!  Also take into account implicitely the parallelism over k points
   do idir1=1,3
    if (ddkfil(idir1)/=0) then
!    Must select the corresponding k point in the full set of k points
!    used in the ddk file : compute the number of k points to skip
     ikpt_fbz_previous(idir1)=ikpt_fbz(idir1)
     ikpt_fbz(idir1)=indkpt1(ikpt)
     nskip=ikpt_fbz(idir1)-ikpt_fbz_previous(idir1)-1
     skipddk(idir1) = skipddk(idir1) + 1 + nskip
     if(nskip/=0)then
      allocate(cg_dum(2,1),eigen_dum(0),kg_dum(3,0),occ_dum(0))
      do ikpt_dum=1+ikpt_fbz_previous(idir1),ikpt_fbz(idir1)-1
       nband_dum=nband(ikpt_dum+(isppol-1)*nkpt)
!      Skip the records whose information is not needed (in case of parallelism)
       tim_rwwf=19
       call WffReadSkipK(1,0,ikpt_dum,isppol,mpi_enreg,wffddk(idir1))
      end do
      deallocate(cg_dum,eigen_dum,kg_dum,occ_dum)
     end if
    end if
   end do


!  allocate(enl1_k(nband_k))
   allocate(occ_k(nband_k))

!  enl1_k(:)=0.0_dp
   d2nl_k(:,:,:)=0.0_dp
   if(prtbbb==1)d2bbb_k(:,:,:,:)=0.0_dp
   kpoint(:)=kpt_rbz(:,ikpt)
   occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)

   if (mkmem==0) then
!   Read (k+G) basis sphere data (same for each spin)
    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,unkg)

!   Read k+g data
    read (unkg) ((kg_k(ii,muig),ii=1,3),muig=1,npw_k)

!   Eventually read (k+G) spherical harmonics
    if (psps%useylm==1) then
     read(unylm)
     read(unylm) ((ylm_k(muig,ilm),muig=1,npw_k),ilm=1,mpsang*mpsang)
    end if

   else

    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
      ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
     end do
    end if

!   End if for choice governed by mkmem
   end if

   wtk_k=wtk_rbz(ikpt)

!  DEBUG
!  write(6,*)' nstdy3 : wtk_k=',wtk_k
!  ENDDEBUG
   kg1_k(:,:) = 0
   if (mk1mem==0) then
!   Read (k+q+G) basis sphere data (same for each spin)
    call rdnpw(ikpt,isppol,nband_k,npw1_k,nspinor,0,unkg1)

!   Read k+g data
    read (unkg1) ((kg1_k(ii,muig),ii=1,3),muig=1,npw1_k)

!   Eventually read (k+q+G) spherical harmonics
    if (psps%useylm==1) then
     read(unylm1)
     read(unylm1) ((ylm1_k(muig,ilm),muig=1,npw1_k),ilm=1,mpsang*mpsang)
    end if


   else

    kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
    if (psps%useylm==1) then
     do ilm=1,mpsang*mpsang
      ylm1_k(1:npw1_k,ilm)=ylm(1+ikg1:npw1_k+ikg1,ilm)
     end do
    end if

!   End if for choice governed by mk1mem
   end if

!  Compute the eigenvalues, wavefunction,
!  contributions to kinetic energy, nonlocal energy, forces,
!  and update of rhor1 to this k-point and this spin polarization.
!  Note that nstwf3 is called with kpoint, while kpt is used inside vtowfk3
   call nstwf3(atindx,atindx1,cg,cg1,ddkfil,d2bbb_k,d2nl_k,ecut,ecutsm,&
&   eig_k,eig1_k,gmet,gprimd,icg,icg1,idir,ikpt,ipert,&
&   iscf,isppol,istwf_k,kg_k,kg1_k,kpoint,mband,mgfft,mkmem,mk1mem,&
&   mpert,mpi_enreg,mpsang,mpw,mpw1,natom,nattyp,nband_k,nband_rbz,nfft,ngfft,&
&   nkpt_rbz,nline,nloalg,npw_k,npw1_k,nspinor,nsppol,ntypat,n4,n5,n6,occopt, &
&   occ_k,ortalg,ph1d,prtbbb,prtvol,psps,qphon,rmet,&
&   typat,ucvol,wffddk,wffnow,wfftgs,wtk_k,xred,ylm_k,ylm1_k)
   d2nl(:,:,:,idir,ipert)=d2nl(:,:,:,idir,ipert)+d2nl_k(:,:,:)
   if(prtbbb==1)then
    d2bbb(:,:,idir,ipert,:,:) = d2bbb(:,:,idir,ipert,:,:) + &
&    d2bbb_k(:,:,:,:)
   end if

!  DEBUG
!  if (ipert < natom + 1) then
!  write(99,*)'d2nl', d2nl(1,1,natom+2,1,ipert)
!  dotr = 0.0_dp
!  do iband = 1,mband
!  write(99,*)'bbb', d2bbb(1,1,1,ipert,iband,iband)
!  dotr = dotr +  d2bbb(1,1,1,ipert,iband,iband)
!  end do
!  write(99,*)'------------------'
!  write(99,*)dotr
!  write(99,*)''
!  end if
!  ENDDEBUG

!  enl1nk (1+bd1tot_index : nband_k+bd1tot_index) = enl1_k  (:)

!  if(iscf>0)then
!  Accumulate sum over k points for nonlocal and kinetic energies,
!  also accumulate gradients of Enonlocal:
!  do iband=1,nband_k
!  enl1=enl1+wtk_k*occ_k(iband)*enl1_k(iband)
!  end do
!  end if

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

!  End big k point loop
  end do

! End loop over spins
 end do

 if(mpi_enreg%paral_compil_kpt==1)then
  call timab(161,1,tsec)
! BEGIN TF_CHANGES
  call leave_test(mpi_enreg)
! END TF_CHANGES
  write(message,*) ' nstdy3: loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
  call timab(161,2,tsec)
 end if

!Treat now varying occupation numbers
!if(occopt>=3 .and. occopt <=7) then
!SUPPRESSED metallic coding of vtorho

!Treat fixed occupation numbers
!else
 if(mpi_enreg%paral_compil_kpt==1)then

  allocate(buffer1(2*3*mpert),buffer2(2*3*mpert))
! Pack d2nl
  buffer1(1:2*3*mpert)=reshape(d2nl(:,:,:,idir,ipert),(/2*3*mpert/))
! Build sum of everything
  call timab(48,1,tsec)
  call xsum_mpi(buffer1,buffer2,2*3*mpert,spaceworld,ierr)
  call timab(48,2,tsec)
! Unpack the final result
  d2nl(:,:,:,idir,ipert)=reshape(buffer2(:),(/2,3,mpert/))
  deallocate(buffer1,buffer2)

  if(prtbbb==1)then
   allocate(buffer1(2*3*mband*mband),buffer2(2*3*mband*mband))
!  Pack d2bbb
   buffer1(1:2*3*mband*mband)=reshape(d2bbb(:,:,idir,ipert,:,:),(/2*3*mband*mband/))
!  Build sum of everything
   call timab(48,1,tsec)
   call xsum_mpi(buffer1,buffer2,2*3*mband*mband,spaceworld,ierr)
   call timab(48,2,tsec)
!  Unpack the final result
   d2bbb(:,:,idir,ipert,:,:)=reshape(buffer2(:),(/2,3,mband,mband/))
   deallocate(buffer1,buffer2)
  end if

 end if ! mpi_enreg%paral_compil_kpt==1

!End of test on varying or fixed occupation numbers
!end if

!In the case of the strain perturbation time-reversal symmetry will always
!be true so imaginary part of d2nl will be must be set to zero here since
!the symmetry-reduced kpt set will leave a non-zero imaginary part.

 if(ipert==natom+3 .or. ipert==natom+4) d2nl(2,:,:,idir,ipert)=0.0_dp

!In case of electric field ipert1, close the ddk wf files
 do idir1=1,3
  if (ddkfil(idir1)/=0)then
   call WffClose(wffddk(idir1),ierr)
  end if
 end do

!Symmetrize the non-local contributions,
!as was needed for the forces in a ground-state calculation
!However, here the quantity is complex, and there are phases !

!DEBUG
!write(6,*)' nstdy3 : total d2nl before sygra3 '
!do ipert1=1,natom
!do idir1=1,3
!write(6,*)ipert1,idir1,d2nl(1,idir1,ipert1,idir,ipert),&
!&                         d2nl(2,idir1,ipert1,idir,ipert)
!end do
!end do
!write(6,*)' nstdy3 : before sygra3, nsym1 =',nsym1
!do isym=1,nsym1
!write(6,*)' isym, symrc1, indsy1(4,isym,1:natom) '
!write(6,*)isym
!write(6,*)symrc1(1:3,1,isym)
!write(6,*)symrc1(1:3,2,isym)
!write(6,*)symrc1(1:3,3,isym)
!write(6,*)indsy1(4,isym,:)
!enddo
!ENDDEBUG

!Do the transform
 allocate(work1(2,3,natom))
 do ipert1=1,natom
  do idir1=1,3
   work1(1,idir1,ipert1)=d2nl(1,idir1,ipert1,idir,ipert)
   work1(2,idir1,ipert1)=d2nl(2,idir1,ipert1,idir,ipert)
  end do
 end do
 call sygra3(natom,d2nl(:,:,:,idir,ipert),work1,indsy1,nsym1,qphon,symrc1)
 deallocate(work1)

!DEBUG
!write(6,*)' nstdy3 : total d2nl after sygra3 '
!do ipert1=1,natom+2
!do idir1=1,3
!write(6,*)ipert1,idir1,d2nl(1,idir1,ipert1,idir,ipert),&
!&                         d2nl(2,idir1,ipert1,idir,ipert)
!end do
!end do
!write(6,*)
!stop
!ENDDEBUG

 if(sum(ddkfil(:))/=0)then
! Must also symmetrize the electric field perturbation response !
! (XG 000803 This was not implemented until now)
! Get the symmetry matrices in terms of real space basis
  allocate(symrl1(3,3,nsym1))
  do isym=1,nsym1
   call mati3inv(symrc1(:,:,isym),symrl1(:,:,isym))
  end do
! There should not be any imaginary part, but stay general (for debugging)
  d2nl_elfd(:,:)=d2nl(:,:,natom+2,idir,ipert)
  do ii=1,3
   sumelfd(:)=0._dp
   do isym=1,nsym1
    do jj=1,3
     if(symrl1(ii,jj,isym)/=0)then
      if(ddkfil(jj)==0)then
       blkflg(ii,natom+2,idir,ipert)=0
      end if
     end if
    end do
    sumelfd(:)=sumelfd(:)+dble(symrl1(ii,1,isym))*d2nl_elfd(:,1)+&
&    dble(symrl1(ii,2,isym))*d2nl_elfd(:,2)+&
&    dble(symrl1(ii,3,isym))*d2nl_elfd(:,3)
   end do
   d2nl(:,ii,natom+2,idir,ipert)=sumelfd(:)/dble(nsym1)
  end do

  if ((prtbbb==1).and.(ipert<=natom)) then
   do iband = 1,mband
    d2nl_elfd(:,:)=d2bbb(:,:,idir,ipert,iband,iband)
    do ii=1,3
     sumelfd(:)=0._dp
     do isym=1,nsym1
!     do jj=1,3
!     if(symrl1(ii,jj,isym)/=0)then
!     if(ddkfil(jj)==0)then
!     blkflg(ii,natom+2,idir,ipert)=0
!     end if
!     end if
!     end do
      sumelfd(:)=sumelfd(:)+dble(symrl1(ii,1,isym))*d2nl_elfd(:,1)+&
&      dble(symrl1(ii,2,isym))*d2nl_elfd(:,2)+&
&      dble(symrl1(ii,3,isym))*d2nl_elfd(:,3)
     end do
     d2bbb(:,ii,idir,ipert,iband,iband)=sumelfd(:)/dble(nsym1)
    end do
   end do  !iband
  end if

  deallocate(symrl1)

 end if

!----------------------------------------------------------------------------
!Now, treat the local contribution

 allocate(vpsp1(cplex*nfft))
 if (ipert /= natom + 1) then
  n3xccc=0
  if(psps%n1xccc/=0)n3xccc=cplex*nfft
  allocate(xccc3d1(cplex*nfft))
  allocate(vxc1(cplex*nfft,nspden),rhodummy(nfft,nspden))

  do ipert1=1,mpert
   do idir1=1,3
    if(ipert1 <= natom)then

!    Get first-order local potential.
     call vloca3(atindx,cplex,gmet,gsqcut,idir1,ipert1,mpi_enreg,psps%mqgrid_ff,natom,&
&     nattyp,nfft,ngfft,ntypat,n1,n2,n3,paral_kgb,ph1d,psps%qgrid_ff,&
&     qphon,ucvol,psps%vlspl,vpsp1,xred)

!    Get first-order exchange-correlation potential (core-correction contribution only !)
     if(psps%n1xccc/=0)then
      call mkcor3(cplex,idir1,ipert1,natom,ntypat,n1,psps%n1xccc,&
&      n2,n3,qphon,rprimd,typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
      option=0
      call mkvxc3(cplex,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,n3xccc,&
&      option,dtset%paral_kgb,qphon,rhodummy,rprimd,vxc1,xccc3d1)
     else
      vxc1(:,:)=zero
     end if ! psps%n1xccc/=0

!    Combines density j2 with local potential j1 (vpsp1 and vxc1)
!    XG030514 : this is a first possible coding, however, each dotprod contains
!    a parallel section (reduction), so it is better to use only one dotprod ...
!    call dotprod_vn(cplex,rhor1,dr_psp1,di_psp1,mpi_enreg,nfft,nfftot,1,2,vpsp1,ucvol)
!    call dotprod_vn(cplex,rhor1,dr_xc1,di_xc1,mpi_enreg,nfft,nfftot,nspden,2,vxc1,ucvol)
!    dotr=dr_psp1+dr_xc1
!    doti=di_psp1+di_xc1
!    ... but then, one needs to overload vxc1
     do ispden=1,min(nspden,2)
      do ifft=1,cplex*nfft
       vxc1(ifft,ispden)=vxc1(ifft,ispden)+vpsp1(ifft)
      end do
     end do
     call dotprod_vn(cplex,rhor1,dotr,doti,mpi_enreg,nfft,nfftot,nspden,2,vxc1,ucvol)
!    MVeithen 021212 : in case ipert = 2, these lines compute the local part
!    of the Born effective charges from phonon and electric
!    field type perturbations, see eq. 43 of
!    X. Gonze and C. Lee, PRB 55, 10355 (1997)
!    The minus sign is due to the fact that the effective charges
!    are minus the second derivatives of the energy
     if (ipert == natom+2) then
      d2lo(1,idir1,ipert1,idir,ipert)=-dotr
      d2lo(2,idir1,ipert1,idir,ipert)=-doti
     else
      d2lo(1,idir1,ipert1,idir,ipert)=dotr
      d2lo(2,idir1,ipert1,idir,ipert)=doti
     end if
!    Endif ipert1<=natom
    end if
   end do
  end do

  deallocate(rhodummy,vxc1,xccc3d1)

 end if ! ipert /= natom +1

 deallocate(d2bbb_k,d2nl_k,kg_k,kg1_k,vpsp1)
 deallocate(eig_k,eig1_k)

 call timab(101,2,tsec)

!DEBUG
!write(6,*)' nstdy3 : exit '
!write(6,*)' nstdy3 : d2nl(:,3,3,3,3)=',d2nl(:,3,3,3,3)
!stop
!ENDDEBUG

end subroutine nstdy3
!!***
