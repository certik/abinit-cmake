!{\src2tex{textfont=tt}}
!!****f* abinit/lobpcgiiwf
!! NAME
!! lobpcgiiwf
!!
!! FUNCTION
!! this routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GZ,AR,MT)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/infos/contributors .
!!
!! INPUTS
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variales for this dataset
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mcg=second dimension of the cg array
!!  mgfft=maximum size of 1d ffts
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nbdblock : number of blocks
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear)
!!  n4,n5,n6 used for dimensionning of vlocal
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  use_subovl= 1 if "subovl" array is computed (see below)
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!
!! OUTPUT
!!  resid_k(nband_k)=residuals for each states
!!  subham(nband_k*(nband_k+1))=the matrix elements of h
!!  If gs_hamk%usepaw==0:
!!    gsc(2,mgsc)=<g|s|c> matrix elements (s=overlap)
!!    subvnl(nband_k*(nband_k+1)*(1-gs_hamk%usepaw))=the matrix elements of vnl
!!  If use_subovl==0:
!!    subovl(nband_k*(nband_k+1)*use_subovl)=the matrix elements of s
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=updated wavefunctions
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      dgemm,dtrsm,getghc,nonlop,orthonormalize
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine lobpcgIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
&           kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&
&           nband_k,nbdblock,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,prtvol,&
&           psps,resid_k,subham,subovl,subvnl,use_subovl,vlocal)

use defs_basis
use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_14wfs
 use interfaces_18seqpar, except_this_one => lobpcgIIwf
 use interfaces_linalg
!End of the abilint section

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
integer :: dimffnl,icg,igsc,lmnmax,matblk,mcg,mgsc,mgfft,mpsang,mpssoang,n4,n5,n6
integer :: natom,nband_k,nbdblock,npw_k,nspinor,ntypat,nvloc,prtvol,use_subovl
type(datafiles_type) :: dtfil
type(dataset_type) :: dtset
type(gs_hamiltonian_type) :: gs_hamk
type(pseudopotential_type) :: psps
type(mpi_type) :: mpi_enreg
integer :: kg_k(3,npw_k)
real(dp) :: cg(2,mcg),ffnl(npw_k,dimffnl,lmnmax,ntypat),gsc(2,mgsc)
real(dp) :: kinpw(npw_k),ph3d(2,npw_k,matblk),resid_k(nband_k)
real(dp) :: vlocal(n4,n5,n6,nvloc)
real(dp) :: subham(nband_k*(nband_k+1)),subvnl(nband_k*(nband_k+1)*(1-gs_hamk%usepaw))
real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)

!Local variables-------------------------------
integer :: spacecomm=0
integer, save :: frozen_count = 0
integer :: iblock,ii,ipw,ipw1,ivectsize,isubo,isubh,iterationnumber,iwavef,jwavef,maxiterations
integer :: bigorder,info,lwork,tim_getghc,inbdblock,choice,idir,tim_nonlop
integer :: cpopt,paw_opt, nkpg,nnlout,optekin=0, signs,sij_opt
integer :: iblocksize,jblocksize, i1, i2, i3, i4
integer :: rvectsize,vectsize,blocksize,bblocksize,iband,istwf_k
integer :: cgindex,gscindex,littleblocksize
logical :: gen_eigenpb
real(dp) :: cgreipw,cgimipw,cscre,chcre,cvcre,dum,lambda_i,sq2
character(len=500) :: message
logical, allocatable :: pflag(:)
real(dp) :: tsec(2)
real(dp), allocatable :: blockvectorx(:,:),blockvectorax(:,:),blockvectorbx(:,:),&
& blockvectorr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
& blockvectorz(:,:),blockvectoraz(:,:),blockvectorbz(:,:),&
& blockvectordumm(:,:),&
& blockvectory(:,:),blockvectorby(:,:),&
& gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
& grampap(:,:),&
& gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
& grampbp(:,:),&
& identity(:,:),coordx(:,:),diagcoordx(:,:),blockvectorxc(:,:),lambda(:,:),&
& grama(:,:),gramb(:,:),gramyx(:,:),&
& kpg_dum(:,:),w(:),work(:),dummy1(:),dummy2(:,:),dummy3(:,:,:)
real(dp), allocatable ::blockvectorp(:,:),blockvectorap(:,:),blockvectorbp(:,:)
real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:)
real(dp), allocatable :: residualnorms(:),eigen(:),rwork(:)
!this is for the call to lobpiii, transfer of information from blockvector to vector
real(dp), allocatable ::vectorx(:),vectorbx(:),vectorax(:),vectory(:),vectorby(:),&
&                             vectorp(:),vectorbp(:),vectorap(:),vectorr(:)
 type(cprj_type) :: cprj_dum(1,1)

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DGEMM' :: dgemm
!DEC$ ATTRIBUTES ALIAS:'DPOTRF' :: dpotrf
!DEC$ ATTRIBUTES ALIAS:'DSYEV' :: dsyev
!DEC$ ATTRIBUTES ALIAS:'DSYGV' :: dsygv
!DEC$ ATTRIBUTES ALIAS:'DTRSM' :: dtrsm
#endif

!NO_ABIRULES
!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band
cgindex(iblocksize)=npw_k*nspinor*(iblocksize-1)+icg+1
gscindex(iblocksize)=npw_k*nspinor*(iblocksize-1)+igsc+1

! *********************************************************************

!debug
!write(6,*) npw_k*nspinor
!write(6,*) cgindex
!write(6,*) size(cg,1),size(cg,2)
!enddebug

call timab(530,1,tsec)

resid_k=zero

if(mod(nband_k,nbdblock)/=0) then
 write(message, '(a,a,a,a,a,a)' ) ch10,&
 &   ' vtowfk : ERROR -',ch10,&
 &   '  For the moment, nband must be a multiple of nbdblock with wfoptalg=5 !',ch10,&
 &   '  Action : raise nband or change nbdblock '
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
end if

gen_eigenpb=(gs_hamk%usepaw==1)
sq2=sqrt(2.0_dp)
if (mpi_enreg%me_g0 == 1) then
 vectsize=2*npw_k*nspinor-1
else
 vectsize=2*npw_k*nspinor
end if
rvectsize=npw_k*nspinor

istwf_k=gs_hamk%istwf_k
maxiterations=dtset%nline

!Big loop bands inside blocks
do iblock=1,nbdblock

 blocksize=(nband_k-1)/nbdblock+1
 bblocksize=(iblock-1)*blocksize 

 if (bblocksize > 0) then
  print*,'BUG - current version of LOBPCG II algorithm does not hold for bblocksize /= 0'
  call leave_new('PERS')
 endif

 !allocations
 allocate(blockvectorx(vectsize,blocksize),blockvectorax(vectsize,blocksize))
 allocate(blockvectorbx(vectsize,blocksize))
 allocate(blockvectorr(vectsize,blocksize),blockvectorar(vectsize,blocksize))
 allocate(blockvectorbr(vectsize,blocksize))
 allocate(blockvectorp(vectsize,blocksize),blockvectorap(vectsize,blocksize),blockvectorbp(vectsize,blocksize))
 allocate(blockvectordumm(vectsize,blocksize),blockvectorxc(vectsize,blocksize))
 allocate(blockvectory(vectsize,bblocksize),blockvectorby(vectsize,bblocksize))
 allocate(gramyx(bblocksize,blocksize))
 allocate(gramxax(blocksize,blocksize),&
 & gramxar(blocksize,blocksize),gramxap(blocksize,blocksize),&
 & gramrar(blocksize,blocksize),gramrap(blocksize,blocksize),&
 & grampap(blocksize,blocksize),&
 & gramxbx(blocksize,blocksize),gramxbr(blocksize,blocksize),&
 & gramxbp(blocksize,blocksize),gramrbr(blocksize,blocksize),&
 & gramrbp(blocksize,blocksize),&
 & grampbp(blocksize,blocksize))
 allocate(lambda(blocksize,blocksize))
 allocate(residualnorms(blocksize),pflag(blocksize))
 pflag = .false.

 !transfer array of wf coeff in block to blockvectorx
 blockvectorr(:,:)=zero
 do iblocksize=1,blocksize
  iband=iblocksize+bblocksize
  if (mpi_enreg%me_g0 == 1) then
   blockvectorx(1,iblocksize)=cg(1,cgindex(iband))
   blockvectorx(2:rvectsize,iblocksize)=cg(1,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
   blockvectorx(rvectsize+1:vectsize,iblocksize)=cg(2,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
  else
   blockvectorx(1:rvectsize,iblocksize)=cg(1,cgindex(iband):cgindex(iband+1)-1)*sq2
   blockvectorx(rvectsize+1:vectsize,iblocksize)=cg(2,cgindex(iband):cgindex(iband+1)-1)*sq2
  end if
 end do
 !if(mpi_enreg%me_group==0)write(6,*) 'bvx entree',blockvectorx(10,7)

 if (iblock /=1) then
  !transfer array of wf coeff less than iblock to blockvectory (not done)
  !transfer cg to blockvectory, for the previous band index
  do iblocksize=1,bblocksize
   iband=iblocksize
   if (mpi_enreg%me_g0 == 1) then
    blockvectory(1,iblocksize)=cg(1,cgindex(iband))
    blockvectory(2:rvectsize,iblocksize)=cg(1,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
    blockvectory(rvectsize+1:vectsize,iblocksize)=cg(2,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
   else
    blockvectory(1:rvectsize,iblocksize)=cg(1,cgindex(iband):cgindex(iband+1)-1)*sq2
    blockvectory(rvectsize+1:vectsize,iblocksize)=cg(2,cgindex(iband):cgindex(iband+1)-1)*sq2
   end if
  end do
  !b-orthogonalize x to the constraint y(supposed b-orthonormal)
  !  blockvectorx=blockvectorx-&
  !              &matmul(blockvectory,matmul(transpose(blockvectory),blockvectorx))
  !  call operators(blockvectory,blockvectorby)
  if(gen_eigenpb) then
   allocate(cwavef(2,npw_k*nspinor))
   allocate(gwavef(2,npw_k*nspinor))
   do iblocksize=1,bblocksize
    if (mpi_enreg%me_g0 == 1) then
     cwavef(1,2:npw_k*nspinor)=blockvectory(2:npw_k*nspinor,iblocksize)/sq2
     cwavef(2,2:npw_k*nspinor)=blockvectory(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)/sq2
     cwavef(1,1)=blockvectory(1,iblocksize)
     cwavef(2,1)=zero
    else
     cwavef(1,1:npw_k*nspinor)=blockvectory(1:npw_k*nspinor,iblocksize)/sq2
     cwavef(2,1:npw_k*nspinor)=blockvectory(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)/sq2      !a verifier
    end if
    !   call to nonlop: compute <g|s|c>
    choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; paw_opt=3 ; cpopt=-1 ; nnlout=0 ; nkpg=0
    call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
    &               dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
    &               istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
    &               mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
    &               gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
    &               gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
    &               gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
    if (mpi_enreg%me_g0 == 1) then
     blockvectorby(2:npw_k*nspinor,iblocksize)=gwavef(1,2:npw_k*nspinor)*sq2
     blockvectorby(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)=gwavef(2,2:npw_k*nspinor)*sq2
     blockvectorby(1,iblocksize)=gwavef(1,1)
    else
     blockvectorby(1:npw_k*nspinor,iblocksize)=gwavef(1,1:npw_k*nspinor)*sq2
     blockvectorby(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)=gwavef(2,1:npw_k*nspinor)*sq2
    end if
   end do
   deallocate(cwavef,gwavef)
  else
   blockvectorby(:,:)=blockvectory(:,:)
  end if

  !orthogonalize x to the constraint y(supposed orthonormal)
  !  blockvectorx=blockvectorx-&
  !              &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorx))
  call dgemm('t','n',bblocksize,blocksize,vectsize,one,blockvectorby,&
  &               vectsize,blockvectorx,vectsize,zero,gramyx,bblocksize)
  call dgemm('n','n',vectsize,blocksize,bblocksize,one,blockvectory,&
  &               vectsize,gramyx,bblocksize,zero,blockvectordumm,vectsize)
  blockvectorx=blockvectorx-blockvectordumm
 end if

 !compute right hand side
 !call operators(blockvectorx,blockvectorbx)
 if(gen_eigenpb) then
  allocate(cwavef(2,npw_k*nspinor))
  allocate(gwavef(2,npw_k*nspinor))
  do iblocksize=1,blocksize
   if (mpi_enreg%me_g0 == 1) then
    cwavef(1,2:npw_k*nspinor)=blockvectorx(2:npw_k*nspinor,iblocksize)/sq2
    cwavef(2,2:npw_k*nspinor)=blockvectorx(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)/sq2
    cwavef(1,1)=blockvectorx(1,iblocksize)
    cwavef(2,1)=zero
   else
    cwavef(1,1:npw_k*nspinor)=blockvectorx(1:npw_k*nspinor,iblocksize)/sq2
    cwavef(2,1:npw_k*nspinor)=blockvectorx(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)/sq2      !a verifier
   end if
   !  call to nonlop: compute <g|s|c>
   choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
   call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
   &              dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
   &              istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
   &              mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
   &              gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
   &              gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
   &              gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
   if (mpi_enreg%me_g0 == 1) then
    blockvectorbx(2:npw_k*nspinor,iblocksize)=gwavef(1,2:npw_k*nspinor)*sq2
    blockvectorbx(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)=gwavef(2,2:npw_k*nspinor)*sq2
    blockvectorbx(1,iblocksize)=gwavef(1,1)
   else
    blockvectorbx(1:npw_k*nspinor,iblocksize)=gwavef(1,1:npw_k*nspinor)*sq2
    blockvectorbx(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)=gwavef(2,1:npw_k*nspinor)*sq2
   end if
  end do
  deallocate(cwavef,gwavef)
 else
  blockvectorbx(:,:)=blockvectorx(:,:)
 end if

 !orthogonalize x
 !call zorthonormalize(blockvectorx,blockvectorbx)
 call orthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,gramxbx,vectsize)
 !if(mpi_enreg%me_group==0)write(6,*) 'bvx ortho',blockvectorx(10,7)
 call dtrsm('r','u','n','n',vectsize,blocksize,one,gramxbx,blocksize,&
 &              blockvectorbx,vectsize)
 !print*,'gramxbx 1'
 !do ii = 1,blocksize
 !print*,gramxbx(ii,:)
 !enddo

 !call operatorh(blockvectorx,blockvectorax)
 allocate(cwavef(2,npw_k*nspinor*blocksize))
 allocate(gwavef(2,npw_k*nspinor*blocksize),gvnlc(2,npw_k*nspinor*blocksize))
 do iblocksize=1,blocksize
  iband=iblocksize
  if (mpi_enreg%me_g0 == 1) then
   cwavef(1,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(2:npw_k*nspinor,iblocksize)/sq2
   cwavef(2,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)/sq2
   cwavef(2,cgindex(iband))=zero
   cwavef(1,cgindex(iband))=blockvectorx(1,iblocksize)
  else
   cwavef(1,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(1:npw_k*nspinor,iblocksize)/sq2
   cwavef(2,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)/sq2
  end if
 end do
 tim_getghc=1 ; sij_opt=0
 call getghc(cwavef,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
      & kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,nspinor,ntypat,&
 & nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
 do iblocksize=1,blocksize
  iband=iblocksize
  if (mpi_enreg%me_g0 == 1) then
   blockvectorax(2:npw_k*nspinor,iblocksize)=gwavef(1,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
   blockvectorax(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)=gwavef(2,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
   blockvectorax(1,iblocksize)=gwavef(1,cgindex(iband))
  else
   blockvectorax(1:npw_k*nspinor,iblocksize)=gwavef(1,cgindex(iband):cgindex(iband+1)-1)*sq2
   blockvectorax(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)=gwavef(2,cgindex(iband):cgindex(iband+1)-1)*sq2
  end if
 end do
 deallocate(cwavef,gwavef,gvnlc)

 !do rayleigh ritz on a in space x
 !gramxax=matmul(transpose(blockvectorx),blockvectorax)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorx,&
 &               vectsize,blockvectorax,vectsize,zero,gramxax,blocksize)

 allocate(eigen(blocksize))

 !print*,'gramxax 2'
 !do ii = 1,blocksize
 !print*,gramxax(ii,:)
 !enddo


 !write(6,*)'gra',gramxax
 !write(6,*)'bx',blockvectorx
 !write(6,*)'bax',blockvectorax
 !write(6,*)'bbx',blockvectorbx
 !call la_syev(gramxax,eigen,jobz='v')
 lwork=3*blocksize
 allocate(work(lwork))
 call dsyev('v','u',blocksize,gramxax,blocksize,eigen,work,lwork,info)
 !if(mpi_enreg%me_group==0)write(6,*)'gramxax apres zheev',gramxax(:,:)
 deallocate(work)
 !blockvectorx=matmul(blockvectorx,gramxax)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
 &               vectsize,gramxax,blocksize,zero,blockvectordumm,vectsize)
 blockvectorx=blockvectordumm
 !blockvectorax=matmul(blockvectorax,gramxax)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
 &               vectsize,gramxax,blocksize,zero,blockvectordumm,vectsize)
 blockvectorax=blockvectordumm
 !blockvectorbx=matmul(blockvectorbx,gramxax)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
 &               vectsize,gramxax,blocksize,zero,blockvectordumm,vectsize)
 blockvectorbx=blockvectordumm
 do iblocksize=1,blocksize
  lambda(iblocksize,iblocksize)=eigen(iblocksize)
 end do
 !DEBUG
 !write(6,*)'lambda',eigen

 !now the main alogrithm
 !allocate(vectorx(vectsize),vectorbx(vectsize),vectorax(vectsize),vectory(vectsize),&
 !     &  vectorby(vectsize),vectorp(vectsize),vectorbp(vectsize),vectorap(vectsize),&
 !     &  vectorr(vectsize))!,&

 iter: do iterationnumber=1,maxiterations
  !DEBUG
  write(6,*)'iterationnumber',iterationnumber
  !if(mpi_enreg%me_group==0)write(6,*) 'bvx',blockvectorx(10,7)

  !passing x into z 
  allocate(blockvectorz(vectsize,blocksize))
  allocate(blockvectoraz(vectsize,blocksize))
  allocate(blockvectorbz(vectsize,blocksize))
  blockvectorz = blockvectorx
  blockvectoraz = blockvectorax
  blockvectorbz = blockvectorbx

  do iblocksize=1,blocksize
   !DEBUG
   !write(6,*)'eig number',iblocksize
   !vectorx(:)=blockvectorx(:,iblocksize)
   !vectorbx(:)=blockvectorbx(:,iblocksize)
   !vectorax(:)=blockvectorax(:,iblocksize)
   !vectorp(:)=blockvectorp(:,iblocksize)
   !vectorbp(:)=blockvectorbp(:,iblocksize)
   !vectorap(:)=blockvectorap(:,iblocksize)
   lambda_i=lambda(iblocksize,iblocksize)
   littleblocksize=1
   print*,'eigenvalue number in ',iblocksize, pflag(iblocksize)
   call lobpcgiiiwf(cg,dimffnl,dtfil,dtset,&
   & ffnl,gs_hamk,gsc,icg,igsc,iblock,&
   & kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,&
   & mpssoang,natom,nbdblock,nband_k,npw_k,nspinor,ntypat,&
   & nvloc,n4,n5,n6,&
   & ph3d,prtvol,psps,resid_k,vlocal,&
   & subham,subvnl,subovl,&
   & littleblocksize,iblocksize,vectsize,&!bblocksize,vectsize,&
   & pflag(iblocksize), frozen_count,&
   & blockvectorx(:,iblocksize:iblocksize),&
   & blockvectorbx(:,iblocksize:iblocksize),&
   & blockvectorax(:,iblocksize:iblocksize),&
   & blockvectorz(:,1:iblocksize), blockvectorbz(:,1:iblocksize),&!& blockvectory,blockvectorby,&
   & lambda(iblocksize:iblocksize,iblocksize:iblocksize),&
   & blockvectorp(:,iblocksize:iblocksize),&
   & blockvectorbp(:,iblocksize:iblocksize),&
   & blockvectorap(:,iblocksize:iblocksize)&
   & )
   print*,'eigenvalue number out ',iblocksize, pflag(iblocksize)
   !blockvectorx(:,iblocksize)=vectorx(:)
   !blockvectorbx(:,iblocksize)=vectorbx(:)
   !blockvectorax(:,iblocksize)=vectorax(:)
   !blockvectorp(:,iblocksize)=vectorp(:)
   !blockvectorbp(:,iblocksize)=vectorbp(:)
   !blockvectorap(:,iblocksize)=vectorap(:)
  end do

  deallocate(blockvectorz)
  deallocate(blockvectoraz)
  deallocate(blockvectorbz)

  !gramxax
  call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorx,&
  &               vectsize,blockvectorax,vectsize,zero,gramxax,blocksize)
  call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
  &               vectsize,blockvectorx,vectsize,zero,gramxbx,blocksize)
  !  write(6,*)'in iii,xax',gramxax(7,7)
  !  write(6,*)'in iii xbx',gramxbx(7,7)

  !print*,'gramxax 3'
  !do ii = 1,blocksize
  !print*,gramxax(ii,:)
  !enddo
  !print*,'gramxbx 4'
  !do ii = 1,blocksize
  !print*,gramxbx(ii,:)
  !enddo


  lwork=3*blocksize
  allocate(work(lwork))
  call dsygv(1,'v','u',blocksize,gramxax,blocksize,gramxbx,blocksize,eigen,&
  &               work,lwork,info)
  deallocate(work)
  !DEBUG
  !write(6,*)'lambda',eigen
  !write(6,*)' '
  lambda(:,:)=zero
  do iblocksize=1,blocksize
   lambda(iblocksize,iblocksize)=eigen(iblocksize)
  end do
  allocate(coordx(blocksize,blocksize),diagcoordx(blocksize,blocksize))
  coordx=gramxax

  !rotate all the vectors according to coordx

  !here is a choice for p
  diagcoordx=zero
  do iblocksize=1,blocksize
   diagcoordx(iblocksize,iblocksize) = coordx(iblocksize,iblocksize)
   coordx(iblocksize,iblocksize) = zero
  end do
  !blockvectorxc = matmul(blockvectorx,coordx)
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
  &               vectsize,coordx,blocksize,zero,blockvectorxc,vectsize)
  !blockvectorx = matmul(blockvectorx,diagcoordx) + blockvectorxc
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
  &               vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
  blockvectorx = blockvectordumm + blockvectorxc
  !blockvectorp = matmul(blockvectorp,diagcoordx) + blockvectorxc
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorp,&
  &               vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
  blockvectorp = blockvectordumm + blockvectorxc

  !blockvectorxc = matmul(blockvectorbx,coordx)
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
  &               vectsize,coordx,blocksize,zero,blockvectorxc,vectsize)
  !blockvectorbx = matmul(blockvectorbx,diagcoordx) + blockvectorxc
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
  &               vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
  blockvectorbx = blockvectordumm + blockvectorxc
  !blockvectorbp = matmul(blockvectorbp,diagcoordx) + blockvectorxc
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbp,&
  &               vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
  blockvectorbp = blockvectordumm + blockvectorxc

  !blockvectorxc = matmul(blockvectorax,coordx)
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
  &               vectsize,coordx,blocksize,zero,blockvectorxc,vectsize)
  !blockvectorax = matmul(blockvectorax,diagcoordx) + blockvectorxc
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
  &               vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
  blockvectorax = blockvectordumm + blockvectorxc
  !blockvectorap = matmul(blockvectorap,diagcoordx) + blockvectorxc
  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorap,&
  &               vectsize,diagcoordx,blocksize,zero,blockvectordumm,vectsize)
  blockvectorap = blockvectordumm + blockvectorxc

  !DEBUG another choice for p
  !  !blockvectorx = matmul(blockvectorx,coordx)
  !  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
  !       &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
  !  blockvectorx = blockvectordumm
  !  !blockvectorbx = matmul(blockvectorbx,coordx)
  !  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
  !       &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
  !  blockvectorbx = blockvectordumm
  !  !blockvectorax = matmul(blockvectorax,coordx)
  !  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
  !       &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
  !  blockvectorax = blockvectordumm
  !  !blockvectorp = matmul(blockvectorp,coordx)
  !  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorp,&
  !       &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
  !  blockvectorp = blockvectordumm
  !  !blockvectorbp = matmul(blockvectorbp,coordx)
  !  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbp,&
  !       &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
  !  blockvectorbp = blockvectordumm
  !  !blockvectorap = matmul(blockvectorap,coordx)
  !  call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorap,&
  !       &               vectsize,coordx,blocksize,zero,blockvectordumm,vectsize)
  !  blockvectorap = blockvectordumm
  !ENDDEBUG

  !DEBUG
  !gramxax
  !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorx,&
  !     &               vectsize,blockvectorax,vectsize,zero,gramxax,blocksize)
  !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
  !     &               vectsize,blockvectorx,vectsize,zero,gramxbx,blocksize)
  !write(6,*)'in iii,xax after',gramxax
  !write(6,*)'in iii xbx after',gramxbx
  !ENDDEBUG
  deallocate(coordx,diagcoordx)

 end do iter

  call precon2(blockvectorbx,lambda,blocksize,&
  &                 istwf_k,kinpw,mpi_enreg,npw_k,nspinor,&
  &                 optekin,blockvectorax,blockvectorr,vectsize)
  
 residualnorms=sum(blockvectorr**2,dim=1)
 resid_k(bblocksize+1:bblocksize+blocksize)=residualnorms(1:blocksize)
 residualnorms=sqrt(residualnorms)
 print*,'residualnorm at the end',residualnorms

 deallocate(eigen)
 !write(6,*)'residualnorm at the end',residualnorms

 !epilogue
 !residualnorms=sqrt(sum(abs(blockvectorr)**2,dim=1))
 ! write(6,*)'residualnorm at the end',residualnorms
 do iblocksize=1,blocksize
  iband=iblocksize+(iblock-1)*blocksize
  if (mpi_enreg%me_g0 == 1) then
   cg(1,cgindex(iband))=blockvectorx(1,iblocksize)
   cg(2,cgindex(iband))=zero
   cg(1,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(2:rvectsize,iblocksize)/sq2
   cg(2,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(rvectsize+1:vectsize,iblocksize)/sq2
  else
   cg(1,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(1:rvectsize,iblocksize)/sq2
   cg(2,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(rvectsize+1:vectsize,iblocksize)/sq2
  end if
 end do
 if(gen_eigenpb) then
  do iblocksize=1,blocksize
   iband=iblocksize+(iblock-1)*blocksize
   if (mpi_enreg%me_g0 == 1) then
    gsc(1,gscindex(iband))=blockvectorbx(1,iblocksize)
    gsc(2,gscindex(iband))=zero
    gsc(1,gscindex(iband)+1:gscindex(iband+1)-1)=blockvectorbx(2:rvectsize,iblocksize)/sq2
    gsc(2,gscindex(iband)+1:gscindex(iband+1)-1)=blockvectorbx(rvectsize+1:vectsize,iblocksize)/sq2
   else
    gsc(1,gscindex(iband):gscindex(iband+1)-1)=blockvectorbx(1:rvectsize,iblocksize)/sq2
    gsc(2,gscindex(iband):gscindex(iband+1)-1)=blockvectorbx(rvectsize+1:vectsize,iblocksize)/sq2
   end if
  end do
 end if

 !this should not exist,since this induce one too much getghc.lazy programming....
 !call operatorh(blockvectorx,blockvectorax,subham,subvnl)!fill also subham, subvnl
 allocate(cwavef(2,npw_k*nspinor*blocksize),gwavef(2,npw_k*nspinor*blocksize))
 allocate(gvnlc(2,npw_k*nspinor*blocksize))
 isubh=1+2*(iblock-1)*blocksize*((iblock-1)*blocksize+1)/2
 do iblocksize=1,blocksize
  iband=iblocksize
  if (mpi_enreg%me_g0 == 1) then
   cwavef(1,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(2:npw_k*nspinor,iblocksize)/sq2
   cwavef(2,cgindex(iband)+1:cgindex(iband+1)-1)=blockvectorx(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)/sq2
   cwavef(2,cgindex(iband))=zero
   cwavef(1,cgindex(iband))=blockvectorx(1,iblocksize)
  else
   cwavef(1,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(1:npw_k*nspinor,iblocksize)/sq2
   cwavef(2,cgindex(iband):cgindex(iband+1)-1)=blockvectorx(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)/sq2
  end if
 end do
 tim_getghc=1; sij_opt=0
 call getghc(cwavef,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
      &  kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,nspinor,ntypat,&
 &  nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
 do iblocksize=1,blocksize
  iband=iblocksize
  if (mpi_enreg%me_g0 == 1) then
   blockvectorax(2:npw_k*nspinor,iblocksize)=gwavef(1,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
   blockvectorax(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)=gwavef(2,cgindex(iband)+1:cgindex(iband+1)-1)*sq2
   blockvectorax(1,iblocksize)=gwavef(1,cgindex(iband))
  else
   blockvectorax(1:npw_k*nspinor,iblocksize)=gwavef(1,cgindex(iband):cgindex(iband+1)-1)*sq2
   blockvectorax(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)=gwavef(2,cgindex(iband):cgindex(iband+1)-1)*sq2
  end if
 end do
 do iblocksize=1,blocksize
  do ii=1,(iblock-1)*blocksize+iblocksize
   iwavef=(ii-1)*npw_k*nspinor+icg
   if (mpi_enreg%me_g0 == 1) then
    ipw1=2;chcre=0.5_dp*cg(1,1+iwavef)*gwavef(1,cgindex(iblocksize))
   else
    ipw1=1;chcre=zero
   end if
   if (gs_hamk%usepaw==1) then
    do ipw=ipw1,npw_k*nspinor
     cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
     chcre=chcre+cgreipw*gwavef(1,ipw+(iblocksize-1)*npw_k*nspinor)+cgimipw*gwavef(2,ipw+(iblocksize-1)*npw_k*nspinor)
    end do
   else
    if (mpi_enreg%me_g0 == 1) then
     cvcre=0.5_dp*cg(1,1+iwavef)*gvnlc(1,cgindex(iblocksize))
    else
     cvcre=zero
    end if
    do ipw=ipw1,npw_k*nspinor
     cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
     chcre=chcre+cgreipw*gwavef(1,ipw+(iblocksize-1)*npw_k*nspinor)+cgimipw*gwavef(2,ipw+(iblocksize-1)*npw_k*nspinor)
     cvcre=cvcre+cgreipw*gvnlc(1,ipw+(iblocksize-1)*npw_k*nspinor)+cgimipw*gvnlc(2,ipw+(iblocksize-1)*npw_k*nspinor)
    end do
    !   store real and imag parts in hermitian storage mode:
    subvnl(isubh)=2.0_dp*cvcre ; subvnl(isubh+1)=zero
   end if
   !  store real and imag parts in hermitian storage mode:
   subham(isubh)=2.0_dp*chcre ; subham(isubh+1)=zero
   isubh=isubh+2
  end do
 end do
 deallocate(cwavef,gwavef,gvnlc)

 !call operators(blockvectorx,blockvectorbx,subovl)!fill also  subovl
 if((gen_eigenpb).and.(use_subovl==1)) then
  allocate(cwavef(2,npw_k*nspinor))
  allocate(gwavef(2,npw_k*nspinor))
  isubo=1+2*(iblock-1)*blocksize*((iblock-1)*blocksize+1)/2
  do iblocksize=1,blocksize
   if (mpi_enreg%me_g0 == 1) then
    cwavef(1,2:npw_k*nspinor)=blockvectorx(2:npw_k*nspinor,iblocksize)/sq2
    cwavef(2,2:npw_k*nspinor)=blockvectorx(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)/sq2
    cwavef(1,1)=blockvectorx(1,iblocksize)
    cwavef(2,1)=zero
   else
    cwavef(1,1:npw_k*nspinor)=blockvectorx(1:npw_k*nspinor,iblocksize)/sq2
    cwavef(2,1:npw_k*nspinor)=blockvectorx(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)/sq2
   end if
   !  call to nonlop: compute <g|s|c>
   choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
   call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
   &              dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
   &              istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
   &              mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
   &              gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
   &              gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
   &              gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
   if (mpi_enreg%me_g0 == 1) then
    blockvectorbx(2:npw_k*nspinor,iblocksize)=gwavef(1,2:npw_k*nspinor)*sq2
    blockvectorbx(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)=gwavef(2,2:npw_k*nspinor)*sq2
    blockvectorbx(1,iblocksize)=gwavef(1,1)
   else
    blockvectorbx(1:npw_k*nspinor,iblocksize)=gwavef(1,1:npw_k*nspinor)*sq2
    blockvectorbx(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)=gwavef(2,1:npw_k*nspinor)*sq2
   end if
   do ii=1,(iblock-1)*blocksize+iblocksize
    iwavef=(ii-1)*npw_k*nspinor+icg
    if (istwf_k==2 .and. mpi_enreg%me_g0 == 1) then
     ipw1=2;cscre=0.5_dp*cg(1,1+iwavef)*gwavef(1,1)
    else
     ipw1=1;cscre=zero
    end if
    do ipw=ipw1,npw_k*nspinor
     cscre=cscre+cg(1,ipw+iwavef)*gwavef(1,ipw)+cg(2,ipw+iwavef)*gwavef(2,ipw)
    end do
    cscre=2.0_dp*cscre
    !   store real and imag parts in hermitian storage mode:
    subovl(isubo)=cscre ; subovl(isubo+1)=zero
    isubo=isubo+2
   end do
  end do
  deallocate(cwavef,gwavef)
 end if
 ! write(6,*)'residualnorm at the end',residualnorms

 !DEBUG
 !write(6,*)'frozen_count',frozen_count,'restart_count',restart_count
 !ENDDEBUG
 ! deallocate(vectorx,vectorbx,vectorax,vectory,&
 !      &  vectorby,vectorp,vectorbp,vectorap,&
 !      &  vectorr)
 deallocate(blockvectorx,blockvectorax,blockvectorbx)
 deallocate(blockvectorr,blockvectorar,blockvectorbr)
 deallocate(blockvectorp,blockvectorap,blockvectorbp)
 deallocate(blockvectory,blockvectorby)
 deallocate(gramyx)
 deallocate(blockvectordumm,blockvectorxc)
 deallocate(gramxax,gramxar,gramxap,gramrar,gramrap,grampap,gramxbx,gramxbr,&
 & gramxbp,gramrbr,gramrbp,grampbp)
 deallocate(lambda)
 deallocate(residualnorms)
 deallocate(pflag)

 !End big loop over bands inside blocks
end do

!  write(6,*)'mpi_enreg%me,icount_ghc',mpi_enreg%me,icount_ghc
!DEBUG
!write(6,*)'end lobpcg'
!stop
!ENDDEBUG
call timab(530,2,tsec)

end subroutine lobpcgiiwf
!!***


subroutine lobpcgiiiwf(cg,dimffnl,dtfil,dtset,&
& ffnl,gs_hamk,gsc,icg,igsc,iblock,&
& kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,&
& mpssoang,natom,nbdblock,nband_k,npw_k,nspinor,ntypat,&
& nvloc,n4,n5,n6,&
& ph3d,prtvol,psps,resid_k,vlocal,&
& subham,subvnl,subovl,blocksize,bblocksize,vectsize,pflag,frozen_count,&
& blockvectorx,blockvectorbx,blockvectorax,blockvectory,blockvectorby,lambda,&
& blockvectorp,blockvectorbp,blockvectorap&
& )

use defs_basis
use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_14wfs
 use interfaces_18seqpar, except_this_one => lobpcgiiiwf
 use interfaces_linalg
!End of the abilint section

implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!arguments ------------------------------------
integer :: spacecomm=0
type(gs_hamiltonian_type) :: gs_hamk
integer :: icg,igsc,lmnmax,matblk,mcg,mgsc,mgfft,mpsang,mpssoang,n4,n5
integer :: n6,natom,nband_k,nbdblock,npw_k,nspinor,ntypat,nvloc,prtvol
integer :: dimffnl
real(dp) :: sq2,tsec
type(datafiles_type) :: dtfil
type(dataset_type) :: dtset
type(pseudopotential_type) :: psps
type(mpi_type) :: mpi_enreg
integer :: kg_k(3,npw_k)
real(dp) :: cg(2,mcg)
real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat),gsc(2,mgsc)
real(dp) :: kinpw(npw_k),ph3d(2,npw_k,matblk),resid_k(nband_k)
real(dp) :: vlocal(n4,n5,n6,nvloc)
real(dp) :: subham(nband_k*(nband_k+1)),subvnl(nband_k*(nband_k+1)),subovl(nband_k*(nband_k+1))
integer:: blocksize,vectsize,bblocksize
logical :: pflag(blocksize)
integer :: frozen_count
real(dp)    :: lambda(blocksize,blocksize)
real(dp) :: blockvectorx(vectsize,blocksize),blockvectorax(vectsize,blocksize),&
&                 blockvectorbx(vectsize,blocksize),blockvectory(vectsize,bblocksize),&
&                 blockvectorby(vectsize,bblocksize)
real(dp) :: blockvectorp(vectsize,blocksize),blockvectorap(vectsize,blocksize),&
&                 blockvectorbp(vectsize,blocksize)
!local variables-------------------------------
integer:: ii,ipw,ipw1,ivectsize,isubo,isubh,iwavef,jwavef,maxiterations
integer:: bigorder,info,lwork,tim_getghc,sij_opt,iminresid
integer:: iblocksize,optekin=0, restart, cond_try, i1, i2, i3, i4,nkpg
integer:: rvectsize,iband,istwf_k,iblock
integer :: choice, cpopt, signs, idir, tim_nonlop, paw_opt, nnlout
logical::gen_eigenpb
real(dp) :: cgreipw,cgimipw,cscre,cscim
real(dp) :: chcre, chcim,cvcre, cvcim, dum
real(dp), parameter :: tolerance1 = 1.e-13
real(dp) :: tolerance2=1.e2
 type(cprj_type) :: cprj_dum(1,1)
!local variables turned arguments--------------



real(dp) :: blockvectorz(vectsize,blocksize),blockvectoraz(vectsize,blocksize),&
&                 blockvectorbz(vectsize,blocksize)
real(dp) :: blockvectorp_old(vectsize,blocksize),blockvectorap_old(vectsize,blocksize),&
&                 blockvectorbp_old(vectsize,blocksize)
real(dp), allocatable ::  blockvectorr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
& blockvectorr1(:,:),&
& blockvectordumm(:,:),&
& blockvectorritz(:,:),&
& blockvectoraritz(:,:),&
& blockvectorbritz(:,:),&
& blockvectorrritz(:,:),&
& blockvectorritzdumm(:,:),&
& dummygram(:,:),&
& gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
& grampap(:,:),&
& gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
& grampbp(:,:),&
& identity(:,:),coordx(:,:),&
& grama(:,:),gramb(:,:),gramyx(:,:),&
      & kpg_dum(:,:),w(:),work(:),dummy1(:),dummy2(:,:),dummy3(:,:,:)
real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:)
!  real(dp) :: zero,one
real(dp), allocatable :: residualnorms(:),eigen(:),rwork(:),residualnorms1(:)
real(dp), allocatable :: residualnormritz(:),normritz(:)

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DTRSM' :: dtrsm
!DEC$ ATTRIBUTES ALIAS:'DGEMM' :: dgemm
!DEC$ ATTRIBUTES ALIAS:'DSYGV' :: dsygv
#endif

if (dtset%useria == 0) then
 tolerance2 = 100.d0
else 
 tolerance2 = real(dtset%useria,dp)
endif

!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band
!  cgindex(iblocksize)=npw_k*nspinor*(iblocksize-1)+icg+1
gen_eigenpb=(gs_hamk%usepaw==1)
!  czero=dcmplx(zero,zero)
!  cone=dcmplx(one,zero)
sq2=sqrt(2.0_dp)
!vectsize=npw_k*nspinor
!blocksize=(nband_k-1)/nbdblock+1
!bblocksize=(iblock-1)*blocksize
istwf_k=gs_hamk%istwf_k
maxiterations=dtset%nline

!passing x into z and p into p_old for restart
blockvectorz = blockvectorx
blockvectoraz = blockvectorax
blockvectorbz = blockvectorbx

!allocations
allocate(blockvectorr(vectsize,blocksize),blockvectorar(vectsize,blocksize))
allocate(blockvectorbr(vectsize,blocksize),blockvectorr1(vectsize,blocksize))
allocate(blockvectordumm(vectsize,blocksize))
allocate(gramyx(bblocksize,blocksize))
allocate(gramxax(blocksize,blocksize),&
& gramxar(blocksize,blocksize),gramxap(blocksize,blocksize),&
& gramrar(blocksize,blocksize),gramrap(blocksize,blocksize),&
& grampap(blocksize,blocksize),&
& gramxbx(blocksize,blocksize),gramxbr(blocksize,blocksize),&
& gramxbp(blocksize,blocksize),gramrbr(blocksize,blocksize),&
& gramrbp(blocksize,blocksize),&
& grampbp(blocksize,blocksize))
allocate(residualnorms(blocksize),residualnorms1(blocksize))

!construct residual
!blockvectorr=blockvectorax-matmul(blockvectorx,lambda)

 call precon2(blockvectorbx,lambda,blocksize,&
 &                 istwf_k,kinpw,mpi_enreg,npw_k,nspinor,&
 &                 optekin,blockvectorax,blockvectorr,vectsize)
 !  blockvectorr(:,iblocksize)=blockvectorax(:,iblocksize)-lambda(iblocksize,iblocksize)*blockvectorbx(:,iblocksize)

residualnorms=sqrt(sum(abs(blockvectorr)**2,dim=1))
!DEBUG
! resid_k(bblocksize+1:bblocksize+blocksize)=residualnorms(1:blocksize)
write(6,*)'residualnorm before lobpcgiii',residualnorms
!ENDEBUG
if(residualnorms(1) > tolerance1)then  !DEBUG this is the wrong condition if blocksize /= 1

 call start_lobpcg ! r orthonormal to x and compute ar 
 !call orthonormalize(blockvectorp,blockvectorbp)
 if(pflag(1)) then   !DEBUG this is the wrong condition if blocksize /= 1
  !DEBUG
  !write(6,*)'blockvectorp,blockvectorbp,blockvectorap'
  !write(6,*)blockvectorp
  !write(6,*)blockvectorbp
  !write(6,*)blockvectorap
  !ENDDEBUG
  !call zorthonormalize(blockvectorp,blockvectorbp,blockvectorap)
  call orthonormalize(blockvectorp,blockvectorbp,blocksize,mpi_enreg,grampbp,vectsize)
  call dtrsm('r','u','n','n',vectsize,blocksize,one,grampbp,blocksize,&
  &              blockvectorbp,vectsize)
  !blockvectorap=matmul(blockvectorap,grampbp)
  call dtrsm('r','u','n','n',vectsize,blocksize,one,grampbp,blocksize,&
  &              blockvectorap,vectsize)
 end if
 if (.not.pflag(1)) then  !DEBUG  this is the wrong condition if blocksize /= 1
  restart=1
 else
  restart=0
  blockvectorp_old = blockvectorp
  blockvectorap_old = blockvectorap
  blockvectorbp_old = blockvectorbp
 end if
 !gramxar=matmul(transpose(blockvectorax),blockvectorr)
 !gramrar=matmul(transpose(blockvectorar),blockvectorr)
 !gramxax=matmul(transpose(blockvectorax),blockvectorx)=lambda
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorax,&
 &               vectsize,blockvectorr,vectsize,zero,gramxar,blocksize)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorar,&
 &               vectsize,blockvectorr,vectsize,zero,gramrar,blocksize)
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorax,&
 !     &               vectsize,blockvectorx,vectsize,zero,gramxax,blocksize)
 gramxax = lambda

 !gramxbx=matmul(transpose(blockvectorbx),blockvectorx)=identity
 !gramrbr=matmul(transpose(blockvectorbr),blockvectorr)=identity
 !gramxbr=matmul(transpose(blockvectorbx),blockvectorr)=zero
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
 !     &               vectsize,blockvectorx,vectsize,zero,gramxbx,blocksize)
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbr,&
 !     &               vectsize,blockvectorr,vectsize,zero,gramrbr,blocksize)
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
 !     &               vectsize,blockvectorr,vectsize,zero,gramxbr,blocksize)
 gramxbx = zero
 gramrbr = zero
 gramxbr = zero
 do iblocksize = 1,blocksize
  gramxbx(iblocksize,iblocksize)=one
  gramrbr(iblocksize,iblocksize)=one
 end do

 i1=0;i2=blocksize;i3=2*blocksize;i4=3*blocksize
 cond: do cond_try=1,2 !2 when restart
  call construct_gram
  !print*,'grama 5'
  !do ii = 1,bigorder
  !print*,grama(ii,:)
  !enddo

  !print*,'gramb 6'
  !do ii = 1,bigorder
  !print*,gramb(ii,:)
  !enddo


  !DEBUG
  !call la_sygv(grama,gramb,eigen,itype=1,jobz='v')
  !ENDDEBUG
  lwork=3*bigorder
  allocate(work(lwork))

  call dsygv(1,'v','u',bigorder,grama,bigorder,gramb,bigorder,eigen,&
  &               work,lwork,info)
  deallocate(work)
  do iblocksize=1,blocksize
   lambda(iblocksize,iblocksize)=eigen(iblocksize)
  end do
  !DEBUG
  !write(6,*)'eigen',eigen(1:blocksize)
  !ENDDEBUG
  coordx=grama(:,1:blocksize)
  call rotate_vectors


   call precon2(blockvectorbx,lambda,blocksize,&
   &                 istwf_k,kinpw,mpi_enreg,npw_k,nspinor,&
   &                 optekin,blockvectorax,blockvectorr1,vectsize)

  residualnorms1=sqrt(sum(abs(blockvectorr1)**2,dim=1))
  !DEBUG
  write(6,*)'residualnorm after lobpcgiii',residualnorms1
  !ENDDEBUG

  !DEBUG
  !print out the B-scalar product of new approximation vector and former iteration smallest eigenvector
  !we use gramb matrix which is of no use anymore
  deallocate(gramb)
  allocate(gramb(bblocksize,1))
  call dgemm('t','n',bblocksize,1,vectsize,one,blockvectorby,&
  &               vectsize,blockvectorx,vectsize,zero,gramb,bblocksize)
  print*,'B-scalar product with former lesser-order iterate',  gramb
  blockvectorr1(:,1) = blockvectorx(:,1)
  if (bblocksize>1)then
   call dgemm('n','n',vectsize,1,bblocksize-1,one,blockvectorby(:,1:bblocksize-1),&
   &               vectsize,gramb(1:bblocksize-1,1),bblocksize-1,zero,blockvectorr1(:,1),vectsize)
   print*,'B-distance XY, sqrt(XX)',sqrt(sum(abs(blockvectorx(:,1)-blockvectorr1(:,1))**2)),sqrt(abs(gramb(bblocksize,1)))
  else
   print*,'B-distance XY, sqrt(XX)',' NA',sqrt(abs(gramb(bblocksize,1)))
  endif

  !ENDDEBUG

  do iblocksize=1,blocksize    !DEBUG this do will work only if blocksize = 1
   if (residualnorms1(iblocksize) > tolerance2*residualnorms(iblocksize)) then
    write(std_out,*) 'restart apply',restart
    !the eigenvector we seek is one of the other Ritz vector
    if (restart==0) then
     call  apply_flip_flop
     exit cond
     deallocate(blockvectorritz,blockvectoraritz,blockvectorbritz,blockvectorritzdumm,blockvectorrritz)
     deallocate(residualnormritz,normritz)
    end if
   else
    pflag = .true.; write(6,*) 'set pftrue'
    exit cond
   end if
  end do
  deallocate(grama,gramb,eigen)
 end do cond
else
 blockvectorp = zero
 blockvectorap = zero
 blockvectorbp = zero
 pflag = .false.; write(6,*) 'set pffalse'
end if

!write(6,*)'blockvectorr',blockvectorr
!write(6,*)'blockvectorx',blockvectorx
!write(6,*)'blockvectorax',blockvectorax
deallocate(blockvectorr,blockvectorar,blockvectorbr,blockvectorr1)
deallocate(gramyx)
deallocate(blockvectordumm)
deallocate(gramxax,gramxar,gramxap,gramrar,gramrap,grampap,gramxbx,gramxbr,&
&      gramxbp,gramrbr,gramrbp,grampbp)
deallocate(residualnorms,residualnorms1)
contains 

subroutine construct_gram



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_linalg
!End of the abilint section

if (restart==0) then
 !gramxap=matmul(transpose(blockvectorax),blockvectorp)
 !gramrap=matmul(transpose(blockvectorar),blockvectorp)
 !grampap=matmul(transpose(blockvectorap),blockvectorp)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorax,&
 &                    vectsize,blockvectorp,vectsize,zero,gramxap,blocksize)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorar,&
 &                    vectsize,blockvectorp,vectsize,zero,gramrap,blocksize)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorap,&
 &                    vectsize,blockvectorp,vectsize,zero,grampap,blocksize)
 bigorder=i4
 allocate(grama(i4,i4),gramb(i4,i4),eigen(i4),coordx(i4,blocksize))
 grama(i1+1:i2,i1+1:i2)=gramxax
 grama(i1+1:i2,i2+1:i3)=gramxar
 grama(i1+1:i2,i3+1:i4)=gramxap
 !grama(i2+1:i3,i1+1:i2)=transpos(gramxar)
 grama(i2+1:i3,i2+1:i3)=gramrar
 grama(i2+1:i3,i3+1:i4)=gramrap
 !grama(i3+1:i4,i1+1:i2)=transpos(gramxap)
 !grama(i3+1:i4,i2+1:i3)=transpos(gramrap)
 grama(i3+1:i4,i3+1:i4)=grampap

 !gramxbp=matmul(transpose(blockvectorbx),blockvectorp)
 !gramrbp=matmul(transpose(blockvectorbr),blockvectorp)
 !grampbp=matmul(transpose(blockvectorbp),blockvectorp)=identity
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
 &                    vectsize,blockvectorp,vectsize,zero,gramxbp,blocksize)
 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbr,&
 &                    vectsize,blockvectorp,vectsize,zero,gramrbp,blocksize)
 !call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbp,&
 !     &                    vectsize,blockvectorp,vectsize,zero,grampbp,blocksize)
 grampbp = zero
 do iblocksize = 1,blocksize
  grampbp(iblocksize,iblocksize)=one
 end do

 gramb(i1+1:i2,i1+1:i2)=gramxbx
 gramb(i1+1:i2,i2+1:i3)=gramxbr
 gramb(i1+1:i2,i3+1:i4)=gramxbp
 !gramb(i2+1:i3,i1+1:i2)=transpos(gramxbr)
 gramb(i2+1:i3,i2+1:i3)=gramrbr
 gramb(i2+1:i3,i3+1:i4)=gramrbp
 !gramb(i3+1:i4,i1+1:i2)=transpos(gramxbp)
 !gramb(i3+1:i4,i2+1:i3)=transpos(gramrbp)
 gramb(i3+1:i4,i3+1:i4)=grampbp
else
 bigorder=i3
 allocate(grama(i3,i3),gramb(i3,i3),eigen(i3),coordx(i3,blocksize))
 grama(i1+1:i2,i1+1:i2)=gramxax
 grama(i1+1:i2,i2+1:i3)=gramxar
 !grama(i2+1:i3,i1+1:i2)=transpos(gramxar)
 grama(i2+1:i3,i2+1:i3)=gramrar
 gramb(i1+1:i2,i1+1:i2)=gramxbx
 gramb(i1+1:i2,i2+1:i3)=gramxbr
 !gramb(i2+1:i3,i1+1:i2)=transpos(gramxbr)
 gramb(i2+1:i3,i2+1:i3)=gramrbr
end if

end  subroutine construct_gram

subroutine rotate_vectors



!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_linalg
!End of the abilint section

if (restart==0) then
 !    blockvectorp=matmul(blockvectorr,coordx(i2+1:i3,:))+&
 !          &matmul(blockvectorp,coordx(i3+1:i4,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorr,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectordumm,vectsize)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorp,&
 &               vectsize,coordx(i3+1:i4,:),blocksize,one,blockvectordumm,vectsize)
 blockvectorp=blockvectordumm
 !    blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))+&
 !          &matmul(blockvectorap,coordx(i3+1:i4,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorar,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectordumm,vectsize)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorap,&
 &               vectsize,coordx(i3+1:i4,:),blocksize,one,blockvectordumm,vectsize)
 blockvectorap=blockvectordumm
 !    blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))+&
 !          &matmul(blockvectorbp,coordx(i3+1:i4,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbr,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectordumm,vectsize)
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbp,&
 &               vectsize,coordx(i3+1:i4,:),blocksize,one,blockvectordumm,vectsize)
 blockvectorbp=blockvectordumm
else
 !blockvectorp =matmul(blockvectorr,coordx(i2+1:i3,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorr,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectorp,vectsize)
 !blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorar,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectorap,vectsize)
 !    blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))
 call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbr,&
 &               vectsize,coordx(i2+1:i3,:),blocksize,zero,blockvectorbp,vectsize)
end if

!blockvectorx = matmul(blockvectorx,coordx(i1+1:i2,:))+blockvectorp
call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
&               vectsize,coordx(i1+1:i2,:),blocksize,zero,blockvectordumm,vectsize)
blockvectorx = blockvectordumm+blockvectorp
!blockvectorax= matmul(blockvectorax,coordx(i1+1:i2,:))+blockvectorap
call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorax,&
&               vectsize,coordx(i1+1:i2,:),blocksize,zero,blockvectordumm,vectsize)
blockvectorax = blockvectordumm+blockvectorap
!blockvectorbx= matmul(blockvectorbx,coordx(i1+1:i2,:))+blockvectorbp
call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbx,&
&               vectsize,coordx(i1+1:i2,:),blocksize,zero,blockvectordumm,vectsize)
blockvectorbx = blockvectordumm+blockvectorbp
deallocate(coordx)
end subroutine rotate_vectors

subroutine apply_flip_flop

!restore former vector


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_14wfs
 use interfaces_linalg
!End of the abilint section

blockvectorp = blockvectorp_old
blockvectorap = blockvectorap_old
blockvectorbp = blockvectorbp_old

blockvectorx = blockvectorz
blockvectorax = blockvectoraz
blockvectorbx = blockvectorbz
!compute all Ritz vectors
allocate(blockvectorritz(vectsize,bigorder),blockvectoraritz(vectsize,bigorder))
allocate(blockvectorbritz(vectsize,bigorder),blockvectorritzdumm(vectsize,bigorder))
allocate(blockvectorrritz(vectsize,bigorder))
allocate(residualnormritz(bigorder),normritz(bigorder))
blockvectorritz(:,i1+1:i2) = blockvectorx(:,:)
blockvectorritz(:,i2+1:i3) = blockvectorr(:,:) 
blockvectorritz(:,i3+1:i4) = blockvectorp(:,:)
blockvectoraritz(:,i1+1:i2) = blockvectorax(:,:)
blockvectoraritz(:,i2+1:i3) = blockvectorar(:,:)
blockvectoraritz(:,i3+1:i4) = blockvectorap(:,:)
blockvectorbritz(:,i1+1:i2) = blockvectorbx(:,:)
blockvectorbritz(:,i2+1:i3) = blockvectorbr(:,:)
blockvectorbritz(:,i3+1:i4) = blockvectorbp(:,:)
!blockvectorritz=matmul(blockvectorritz,grama)
call dgemm('n','n',vectsize,bigorder,bigorder,one,blockvectorritz,&
&               vectsize,grama,bigorder,zero,blockvectorritzdumm,vectsize)
blockvectorritz(:,:)=blockvectorritzdumm(:,:)
!print*,'diff3'
!print*,blockvectorritz(1:10,1)-(grama(1,1)*blockvectorx(1:10,1)+grama(2,1)*blockvectorr(1:10,1)+grama(3,1)*blockvectorp(1:10,1))

!blockvectoraritz=matmul(blockvectoraritz,grama)
call dgemm('n','n',vectsize,bigorder,bigorder,one,blockvectoraritz,&
&               vectsize,grama,bigorder,zero,blockvectorritzdumm,vectsize)
blockvectoraritz(:,:)=blockvectorritzdumm(:,:)
!blockvectorbritz=matmul(blockvectorritz,grama)
call dgemm('n','n',vectsize,bigorder,bigorder,one,blockvectorbritz,&
&               vectsize,grama,bigorder,zero,blockvectorritzdumm,vectsize)
blockvectorbritz(:,:)=blockvectorritzdumm(:,:)
!     !blockvectorp=matmul(blockvectorz,coordx(i1+1:i2,:))+blockvectorp
!     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorz,&
!          &               vectsize,coordx(i1+1:i2,:),blocksize,one,blockvectorp,vectsize)
!     !blockvectorap=matmul(blockvectoraz,coordx(i1+1:i2,:))+blockvectorap
!     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectoraz,&
!          &               vectsize,coordx(i1+1:i2,:),blocksize,one,blockvectorap,vectsize)
!     !blockvectorbp=matmul(blockvectorbz,coordx(i1+1:i2,:))+blockvectorbp
!     call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorbz,&
!          &               vectsize,coordx(i1+1:i2,:),blocksize,one,blockvectorbp,vectsize)

!compute the residuals for the Ritz vectors other than the first
!     residualnormritz(i1+1:i2) = residualnorm1(:)

do ii=i1+1,bigorder
 call precon2(blockvectorbritz(:,ii),eigen(ii),1,&
 &                 istwf_k,kinpw,mpi_enreg,npw_k,nspinor,&
 &                 optekin,blockvectoraritz(:,ii),blockvectorrritz(:,ii),vectsize)
end do
residualnormritz=sqrt(sum(abs(blockvectorrritz)**2,dim=1))
normritz=sqrt(sum(abs(blockvectorritz)**2,dim=1))

print*,'residualnorms1',residualnorms1
print*,'residualnormritz',residualnormritz
print*,'normritz',normritz
print*,'quotient',residualnormritz/normritz
iminresid = 1
do ii = 2,bigorder
 if(residualnormritz(ii)/normritz(ii)<residualnormritz(iminresid)/normritz(iminresid))then
  iminresid = ii
 endif
enddo
print*,'iminresid',iminresid
!following lines will not work with blocksize > 1
blockvectorp(:,1)  = grama(2,iminresid)*blockvectorr(:,1)+grama(3,iminresid)*blockvectorp(:,1)
blockvectorap(:,1) = grama(2,iminresid)*blockvectorar(:,1)+grama(3,iminresid)*blockvectorap(:,1)
blockvectorbp(:,1) = grama(2,iminresid)*blockvectorbr(:,1)+grama(3,iminresid)*blockvectorbp(:,1)
blockvectorx(:,1)  = blockvectorritz(:,iminresid)
blockvectorax(:,1) = blockvectoraritz(:,iminresid)
blockvectorbx(:,1) = blockvectorbritz(:,iminresid)
pflag = .true.; write(6,*) 'set pftrue2'

end subroutine apply_flip_flop

subroutine start_lobpcg

!if(abs(sum(residualnorms)) < 1.d-10) exit
!!$if (.false.)then!(bbblocksize>0) then !(iblock /=1) then !residuals orthogonal to blockvectorby
!!$ !   blockvectorr=blockvectorr-&
!!$ !           &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorr))
!!$ call dgemm('t','n',bblocksize,blocksize,vectsize,one,blockvectorby,&
!!$ &               vectsize,blockvectorr,vectsize,zero,gramyx,bblocksize)
!!$ call dgemm('n','n',vectsize,blocksize,bblocksize,one,blockvectory,&
!!$ &               vectsize,gramyx,bblocksize,zero,blockvectordumm,vectsize)
!!$ blockvectorr=blockvectorr-blockvectordumm
!!$end if
!residuals orthogonal to blockvectorx
!  blockvectorr=blockvectorr-&
!          &matmul(blockvectorx,matmul(transpose(blockvectorbx),blockvectorr))


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_14wfs
 use interfaces_linalg
!End of the abilint section

call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorbx,&
&               vectsize,blockvectorr,vectsize,zero,gramxax,blocksize)
call dgemm('n','n',vectsize,blocksize,blocksize,one,blockvectorx,&
&               vectsize,gramxax,blocksize,zero,blockvectordumm,vectsize)
blockvectorr=blockvectorr-blockvectordumm
!and now (b)orthornormalize r
!call operators(blockvectorr,blockvectorbr)
if (gen_eigenpb) then
 allocate(cwavef(2,npw_k*nspinor))
 allocate(gwavef(2,npw_k*nspinor))
 do iblocksize=1,blocksize
  if (mpi_enreg%me_g0 == 1) then
   cwavef(1,2:npw_k*nspinor)=blockvectorr(2:npw_k*nspinor,iblocksize)/sq2
   cwavef(2,2:npw_k*nspinor)=blockvectorr(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)/sq2
   cwavef(1,1)=blockvectorr(1,iblocksize)
   cwavef(2,1)=zero
  else
   cwavef(1,1:npw_k*nspinor)=blockvectorr(1:npw_k*nspinor,iblocksize)/sq2
   cwavef(2,1:npw_k*nspinor)=blockvectorr(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)/sq2
  end if
  !   call to nonlop: compute <g|s|c>
  choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
  call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy2,&
  &               dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
  &               istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
  &               mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
  &               gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
  &               gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
  &               gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
  if (mpi_enreg%me_g0 == 1) then
   blockvectorbr(2:npw_k*nspinor,iblocksize)=gwavef(1,2:npw_k*nspinor)*sq2
   blockvectorbr(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)=gwavef(2,2:npw_k*nspinor)*sq2
   blockvectorbr(1,iblocksize)=gwavef(1,1)
  else
   blockvectorbr(1:npw_k*nspinor,iblocksize)=gwavef(1,1:npw_k*nspinor)*sq2
   blockvectorbr(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)=gwavef(2,1:npw_k*nspinor)*sq2
  end if
 end do
 deallocate(cwavef,gwavef)
else
 blockvectorbr(:,:) = blockvectorr(:,:)
end if
!call orthonormalize(blockvectorr,blockvectorbr)
call orthonormalize(blockvectorr,blockvectorbr,blocksize,mpi_enreg,gramrbr,vectsize)
call dtrsm('r','u','n','n',vectsize,blocksize,one,gramrbr,blocksize,&
&              blockvectorbr,vectsize)
!compute ar
!blockvectorar=matmul(operatora,blockvectorr)
!call operatorh(blockvectorr,blockvectorar)
allocate(cwavef(2,npw_k*nspinor),gwavef(2,npw_k*nspinor),gvnlc(2,npw_k*nspinor))
do iblocksize=1,blocksize
 iband=iblocksize
 if (mpi_enreg%me_g0 == 1) then
  cwavef(1,2:npw_k*nspinor)=blockvectorr(2:npw_k*nspinor,iblocksize)/sq2
  cwavef(2,2:npw_k*nspinor)=blockvectorr(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)/sq2
  cwavef(2,1)=zero
  cwavef(1,1)=blockvectorr(1,iblocksize)
 else
  cwavef(1,1:npw_k*nspinor)=blockvectorr(1:npw_k*nspinor,iblocksize)/sq2
  cwavef(2,1:npw_k*nspinor)=blockvectorr(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)/sq2
 end if
 tim_getghc=1 ; sij_opt=0
 call getghc(cwavef,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
 &  kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,nspinor,ntypat,&
 &  nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
 if (mpi_enreg%me_g0 == 1) then
  blockvectorar(2:npw_k*nspinor,iblocksize)=gwavef(1,2:npw_k*nspinor)*sq2
  blockvectorar(npw_k*nspinor+1:2*npw_k*nspinor-1,iblocksize)=gwavef(2,2:npw_k*nspinor)*sq2
  blockvectorar(1,iblocksize)=gwavef(1,1)
 else
  blockvectorar(1:npw_k*nspinor,iblocksize)=gwavef(1,:)*sq2
  blockvectorar(npw_k*nspinor+1:2*npw_k*nspinor,iblocksize)=gwavef(2,:)*sq2
 end if
end do
deallocate(cwavef,gwavef,gvnlc)

end subroutine start_lobpcg

end  subroutine lobpcgiiiwf
