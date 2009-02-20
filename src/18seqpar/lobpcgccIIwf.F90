!{\src2tex{textfont=tt}}
!!****f* abinit/lobpcgcciiwf
!! NAME
!! lobpcgcciiwf
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
!!      getghc,nonlop,zgemm,zhegv,zorthonormalize,zprecon3,ztrsm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine lobpcgccIIwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
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
 use interfaces_18seqpar, except_this_one => lobpcgccIIwf
 use interfaces_lib01hidempi
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
  integer :: activepsize,activersize,bblocksize,bigorder,blocksize,cgindex,choice,cpopt
  integer :: cond_try,gscindex,iblocksize
  integer :: iblock,i1,i2,i3,i4,iband,idir,ier,ierr,ii,info
  integer :: ipw,ipw1,istwf_k,isubo,isubh,iterationnumber,ivectsize,iwavef,jblocksize,jwavef,littleblocksize,lwork
  integer :: maxiterations,nkpg,nnlout,old_paral_level,optekin,paw_opt,restart,signs,sij_opt,tim_getghc,tim_nonlop,vectsize
  logical :: gen_eigenpb
  real(dp) :: cgreipw,cgimipw,cscre,cscim,chcre,chcim,cvcre,cvcim,dum,lambda_i,sq2
  character(len=500) :: message
  logical, allocatable :: pflag(:)
  real(dp) :: tsec(2),transf2(nband_k*(nband_k+1),2)
  real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:)
  real(dp), allocatable :: residualnorms(:),eigen(:),rwork(:),lambda(:,:),kpg_dum(:,:)
  complex(dp), allocatable :: blockvectorx(:,:),blockvectorax(:,:),blockvectorbx(:,:),&
       & blockvectorr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
       & blockvectorp(:,:),blockvectorap(:,:),blockvectorbp(:,:),blockvectordumm(:,:),&
       & blockvectory(:,:),blockvectorby(:,:),&
       & gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
       & grampap(:,:),&
       & gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
       & grampbp(:,:),&
       & identity(:,:),coordx(:,:),diagcoordx(:,:),&!lambda(:,:),&
       & grama(:,:),gramb(:,:),gramyx(:,:),&
       & transf(:,:,:),w(:),work(:),&
       & blockvectorxc(:,:)
 real(dp), allocatable :: dummy1(:),dummy2(:,:),dummy3(:,:,:)
  !This is for the call to lobpiii, transfer of information from blockvector to vector
  complex(dp), allocatable :: vectorx(:),vectorbx(:),vectorax(:),vectory(:),vectorby(:),&
       &                            vectorp(:),vectorbp(:),vectorap(:),vectorr(:)
 type(cprj_type) :: cprj_dum(1,1)

#ifdef VMS
  !DEC$ ATTRIBUTES ALIAS:'ZGEMM' :: zgemm
  !DEC$ ATTRIBUTES ALIAS:'ZHEEV' :: zheev
  !DEC$ ATTRIBUTES ALIAS:'ZHEGV' :: zhegv
  !DEC$ ATTRIBUTES ALIAS:'ZPOTRF' :: zpotrf
  !DEC$ ATTRIBUTES ALIAS:'ZTRSM' :: ztrsm
#endif

  !no_abirules
  !correspondence with abinit. here for real wf but in complex mode
  !this is the index of a given band
  cgindex(iblocksize)=npw_k*nspinor*(iblocksize-1)+icg+1
  gscindex(iblocksize)=npw_k*nspinor*(iblocksize-1)+igsc+1

  gen_eigenpb=(gs_hamk%usepaw==1)
  optekin=0;if (dtset%wfoptalg>10) optekin=1

  ! *************************************************************************

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

  sq2=sqrt(two)
  vectsize=npw_k*nspinor
  istwf_k=gs_hamk%istwf_k
  maxiterations=dtset%nline

  !Big loop bands inside blocks
  do iblock=1,nbdblock

     blocksize=(nband_k-1)/nbdblock+1
     bblocksize=(iblock-1)*blocksize

     !allocations
     allocate(blockvectorx(vectsize,blocksize),blockvectorax(vectsize,blocksize))
     allocate(blockvectorbx(vectsize,blocksize))
     allocate(blockvectorr(vectsize,blocksize),blockvectorar(vectsize,blocksize))
     allocate(blockvectorbr(vectsize,blocksize))
     allocate(blockvectorp(vectsize,blocksize),blockvectorap(vectsize,blocksize))
     allocate(blockvectorbp(vectsize,blocksize))
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
     allocate(transf(blocksize,blocksize,3))
     allocate(residualnorms(blocksize))
     allocate(pflag(blocksize))

     pflag = .false.

     !transfer array of wf coeff in block to blockvectorx (complex to complex)
     do iblocksize=1,blocksize
        iband=iblocksize+bblocksize
        blockvectorx(1:vectsize,iblocksize)=dcmplx(cg(1,cgindex(iband):cgindex(iband+1)-1),&
             & cg(2,cgindex(iband):cgindex(iband+1)-1))
     end do

     !transfer array of wf coeff less than iblock to blockvectory (not done)
     if(iblock /=1) then
        !transfer cg to blockvectory, for the previous band index
        do iblocksize=1,bblocksize
           iband=iblocksize
           blockvectory(1:vectsize,iblocksize)=dcmplx(cg(1,cgindex(iband):cgindex(iband+1)-1),&
                &  cg(2,cgindex(iband):cgindex(iband+1)-1))
        end do
        !call operators(blockvectory,blockvectorby)
        if(gen_eigenpb) then
           allocate(cwavef(2,npw_k*nspinor))
           allocate(gwavef(2,npw_k*nspinor))
           do iblocksize=1,bblocksize
              cwavef(1,1:npw_k*nspinor)=real (blockvectory(1:npw_k*nspinor,iblocksize))
              cwavef(2,1:npw_k*nspinor)=aimag(blockvectory(1:npw_k*nspinor,iblocksize))
              !   Call to nonlop: compute <g|S|c>
              choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
              call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
                   &               dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
                   &               istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
                   &               mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
                   &               gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
                   &               gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
                   &               gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
              blockvectorby(1:npw_k*nspinor,iblocksize)=dcmplx(gwavef(1,1:npw_k*nspinor),gwavef(2,1:npw_k*nspinor))
           end do
           deallocate(cwavef,gwavef)
        else
           blockvectorby(:,:)=blockvectory(:,:)
        end if

        !orthogonalize x to the constraint y(supposed orthonormal)
        !  blockvectorx=blockvectorx-&
        !              &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorx))
        call zgemm('c','n',bblocksize,blocksize,vectsize,cone,blockvectorby,&
             &               vectsize,blockvectorx,vectsize,czero,gramyx,bblocksize)
        old_paral_level= mpi_enreg%paral_level
        mpi_enreg%paral_level=3
        call xcomm_init(mpi_enreg,spaceComm)
        call timab(48,1,tsec)
        call xsum_mpi(gramyx,spaceComm,ierr)
        call timab(48,2,tsec)
        mpi_enreg%paral_level= old_paral_level
        call zgemm('n','n',vectsize,blocksize,bblocksize,cone,blockvectory,&
             &               vectsize,gramyx,bblocksize,czero,blockvectordumm,vectsize)
        blockvectorx=blockvectorx-blockvectordumm
     end if
     !compute right hand side
     !call operators(blockvectorx,blockvectorbx)
     if(gen_eigenpb) then
        allocate(cwavef(2,npw_k*nspinor))
        allocate(gwavef(2,npw_k*nspinor))
        do iblocksize=1,blocksize
           cwavef(1,1:npw_k*nspinor)=real (blockvectorx(1:npw_k*nspinor,iblocksize))
           cwavef(2,1:npw_k*nspinor)=aimag(blockvectorx(1:npw_k*nspinor,iblocksize))
           !  Call to nonlop: compute <g|S|c>
           choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
           call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
                &              dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
                &              istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
                &              mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
                &              gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
                &              gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
                &              gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
           blockvectorbx(1:npw_k*nspinor,iblocksize)=dcmplx(gwavef(1,1:npw_k*nspinor),gwavef(2,1:npw_k*nspinor))
        end do
        deallocate(cwavef,gwavef)
     else
        blockvectorbx(:,:)=blockvectorx(:,:)
     end if

     !orthogonalize x
     !call zorthonormalize(blockvectorx,blockvectorbx)
     call zorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,gramxbx,vectsize)
     call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramxbx,blocksize,&
          &              blockvectorbx,vectsize)
     !call operatorh(blockvectorx,blockvectorax)
     allocate(cwavef(2,npw_k*nspinor*blocksize),gwavef(2,npw_k*nspinor*blocksize),gvnlc(2,npw_k*nspinor*blocksize))
     do iblocksize=1,blocksize
        cwavef(1,npw_k*nspinor*(iblocksize-1)+1:npw_k*nspinor*iblocksize)=real(blockvectorx(1:npw_k*nspinor,iblocksize))
        cwavef(2,npw_k*nspinor*(iblocksize-1)+1:npw_k*nspinor*iblocksize)=aimag(blockvectorx(1:npw_k*nspinor,iblocksize))
     end do
     tim_getghc=1 ; sij_opt=0
     call getghc(cwavef,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
          & kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,nspinor,ntypat,&
          & nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
     do iblocksize=1,blocksize
        blockvectorax(1:npw_k*nspinor,iblocksize)=&
             &  dcmplx(gwavef(1,npw_k*nspinor*(iblocksize-1)+1:npw_k*nspinor*iblocksize),&
             &         gwavef(2,npw_k*nspinor*(iblocksize-1)+1:npw_k*nspinor*iblocksize))
     end do
     deallocate(cwavef,gwavef,gvnlc)

     !do rayleigh ritz on a in space x
     !gramxax=matmul(transpose(blockvectorx),blockvectorax)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
          &               vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
     old_paral_level= mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm)
     call timab(48,1,tsec)
     call xsum_mpi(gramxax,spaceComm,ierr)
     call timab(48,2,tsec)
     mpi_enreg%paral_level= old_paral_level
     allocate(eigen(blocksize))

     !call la_syev(gramxax,eigen,jobz='v')
     lwork=3*blocksize-2
     allocate(work(lwork),rwork(lwork))
     !write(6,*)'gramxax bef',gramxax(:,:)
     do iblocksize=1,blocksize
        do jblocksize=1,blocksize
           if(abs(gramxax(iblocksize,jblocksize)) < 1.e-14) then
              gramxax(iblocksize,jblocksize)=czero
           else
              !write(6,*)'gramxax non nul',gramxax(iblocksize,jblocksize)
           end if
        end do
     end do
     call zheev('v','u',blocksize,gramxax,blocksize,eigen,work,lwork,rwork,info)
     deallocate(work,rwork)
     !blockvectorx=matmul(blockvectorx,gramxax)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
          &               vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
     blockvectorx=blockvectordumm
     !blockvectorax=matmul(blockvectorax,gramxax)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
          &               vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
     blockvectorax=blockvectordumm
     !blockvectorbx=matmul(blockvectorbx,gramxax)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
          &               vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
     blockvectorbx=blockvectordumm
     do iblocksize=1,blocksize
        lambda(iblocksize,iblocksize)=eigen(iblocksize)
     end do
     write(6,*)'lambda',eigen

     !now the main alogrithm
!!$     allocate(vectorx(vectsize),vectorbx(vectsize),vectorax(vectsize),vectory(vectsize),&
!!$          &  vectorby(vectsize),vectorp(vectsize),vectorbp(vectsize),vectorap(vectsize),&
!!$          &  vectorr(vectsize))
     iter: do iterationnumber=1,maxiterations
        !    write(6,*)'bvx',blockvectorx(10,7)
        write(6,*)'iterationnumber',iterationnumber
        do iblocksize=1,blocksize
           !vectorx(:)=blockvectorx(:,iblocksize)
           !vectorbx(:)=blockvectorbx(:,iblocksize)
           !vectorax(:)=blockvectorax(:,iblocksize)
           !       vectory(:)=blockvectory(:,iblocksize)! attention, purement blanc si pas de y
           !       vectorby(:)=blockvectorby(:,iblocksize)! a remplacer par by etc....
           !vectorp(:)=blockvectorp(:,iblocksize)
           !vectorbp(:)=blockvectorbp(:,iblocksize)
           !vectorap(:)=blockvectorap(:,iblocksize)
           lambda_i=lambda(iblocksize,iblocksize)
           !     if(iblock > 1) stop('huh')

           littleblocksize=1
           call lobpcgcciiiwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
                & kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,&
                & mpssoang,natom,nbdblock,nband_k,npw_k,nspinor,ntypat,&
                & nvloc,n4,n5,n6,ph3d,prtvol,psps,resid_k,use_subovl,vlocal,&
                & subham,subvnl,subovl,littleblocksize,bblocksize,vectsize,pflag(iblocksize), &
                & blockvectorx (:,iblocksize:iblocksize),&
                & blockvectorbx(:,iblocksize:iblocksize),&
                & blockvectorax(:,iblocksize:iblocksize),&
                & blockvectory,blockvectorby,lambda_i,&
                & blockvectorp (:,iblocksize:iblocksize),&
                & blockvectorbp(:,iblocksize:iblocksize),&
                & blockvectorap(:,iblocksize:iblocksize)&
                &   )

           !blockvectorx(:,iblocksize)=vectorx(:)
           !blockvectorbx(:,iblocksize)=vectorbx(:)
           !blockvectorax(:,iblocksize)=vectorax(:)
           !blockvectorp(:,iblocksize)=vectorp(:)
           !blockvectorbp(:,iblocksize)=vectorbp(:)
           !blockvectorap(:,iblocksize)=vectorap(:)
        end do
        !gramxax
        call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
             &               vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
        call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
             &               vectsize,blockvectorx,vectsize,czero,gramxbx,blocksize)
        ! write(6,*)'in iii,xax'!,gramxax(7,7)
        ! write(6,*)'in iii xbx'!,gramxbx(7,7)
        lwork=3*blocksize-2
        allocate(work(lwork),rwork(lwork))
        call zhegv(1,'v','u',blocksize,gramxax,blocksize,gramxbx,blocksize,eigen,&
             &               work,lwork,rwork,info)
        deallocate(work,rwork)
        write(6,*)'lambda',eigen
        write(6,*)' '
        lambda(:,:)=zero
        do iblocksize=1,blocksize
           lambda(iblocksize,iblocksize)=eigen(iblocksize)
        end do
        allocate(coordx(blocksize,blocksize),diagcoordx(blocksize,blocksize))
        coordx=gramxax
        !rotate all the vectors according to coordx

        !choix de p
        diagcoordx=czero
        do iblocksize=1,blocksize
           diagcoordx(iblocksize,iblocksize) = coordx(iblocksize,iblocksize)
           coordx(iblocksize,iblocksize) = czero
        end do
        !blockvectorxc = matmul(blockvectorx,coordx)
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
             &               vectsize,coordx,blocksize,czero,blockvectorxc,vectsize)
        !blockvectorx = matmul(blockvectorx,diagcoordx) + blockvectorxc
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
             &               vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
        blockvectorx = blockvectordumm + blockvectorxc
        !blockvectorp = matmul(blockvectorp,diagcoordx) + blockvectorxc
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorp,&
             &               vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
        blockvectorp = blockvectordumm + blockvectorxc
        !blockvectorxc = matmul(blockvectorbx,coordx)
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
             &               vectsize,coordx,blocksize,czero,blockvectorxc,vectsize)
        !blockvectorbx = matmul(blockvectorbx,diagcoordx) + blockvectorxc
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
             &               vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
        blockvectorbx = blockvectordumm + blockvectorxc
        !blockvectorbp = matmul(blockvectorbp,diagcoordx) + blockvectorxc
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbp,&
             &               vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
        blockvectorbp = blockvectordumm + blockvectorxc
        !blockvectorxc = matmul(blockvectorax,coordx)
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
             &               vectsize,coordx,blocksize,czero,blockvectorxc,vectsize)
        !blockvectorax = matmul(blockvectorax,diagcoordx) + blockvectorxc
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
             &               vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
        blockvectorax = blockvectordumm + blockvectorxc
        !blockvectorap = matmul(blockvectorap,diagcoordx) + blockvectorxc
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorap,&
             &               vectsize,diagcoordx,blocksize,czero,blockvectordumm,vectsize)
        blockvectorap = blockvectordumm + blockvectorxc

!!$    !autre choix possible
!!$    !blockvectorx = matmul(blockvectorx,coordx)
!!$    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
!!$         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!!$    blockvectorx = blockvectordumm
!!$    !blockvectorbx = matmul(blockvectorbx,coordx)
!!$    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
!!$         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!!$    blockvectorbx = blockvectordumm
!!$    !blockvectorax = matmul(blockvectorax,coordx)
!!$    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
!!$         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!!$    blockvectorax = blockvectordumm
!!$    !blockvectorp = matmul(blockvectorp,coordx)
!!$    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorp,&
!!$         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!!$    blockvectorp = blockvectordumm
!!$    !blockvectorbp = matmul(blockvectorbp,coordx)
!!$    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbp,&
!!$         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!!$    blockvectorbp = blockvectordumm
!!$    !blockvectorap = matmul(blockvectorap,coordx)
!!$    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorap,&
!!$         &               vectsize,coordx,blocksize,czero,blockvectordumm,vectsize)
!!$    blockvectorap = blockvectordumm
        !gramxax
        !    call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
        !         &               vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
        !    call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
        !         &               vectsize,blockvectorx,vectsize,czero,gramxbx,blocksize)
        !write(6,*)'in iii,xax after',gramxax
        !write(6,*)'in iii xbx after',gramxbx
        deallocate(coordx,diagcoordx)
     end do iter
     deallocate(eigen)
     !epilogue
     !gramxbx=matmul(transpose(blockvectorx),blockvectorx)
     if(.true.) then !epilogue
        ! call operators(blockvectorx,blockvectorbx)
        call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
             &               vectsize,blockvectorbx,vectsize,czero,gramxbx,blocksize)
        !blockvectorax=matmul(operatora,blockvectorx)
        ! call operatorh(blockvectorx,blockvectorax)
        !gramxax=matmul(transpose(blockvectorax),blockvectorx)
        call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
             &               vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
        allocate(eigen(blocksize))
        !call la_sygv(gramxax,gramxbx,eigen,itype=1,jobz='v')
        lwork=3*blocksize-2
        allocate(work(lwork),rwork(lwork))
        call zhegv(1,'v','u',blocksize,gramxax,blocksize,gramxbx,blocksize,eigen,&
             &               work,lwork,rwork,info)
        deallocate(work,rwork)
        lambda=czero
        do iblocksize=1,blocksize
           lambda(iblocksize,iblocksize)=eigen(iblocksize)
        end do
        ! write(6,*)'gramxax'
        ! write(6,*)gramxax
        write(6,*)'eigen at the end',eigen
        !debug
        !blockvectorr=blockvectorax-matmul(blockvectorx,lambda)
        !residualnorms=sqrt(sum(blockvectorr**2,dim=1))
        !write(6,*)'residualnorm at the end bef orth',residualnorms
        !debug
        !blockvectorx=matmul(blockvectorx,gramxax)
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
             &               vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
        blockvectorx=blockvectordumm
        !blockvectorax=matmul(blockvectorax,gramxax)
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
             &               vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
        blockvectorax=blockvectordumm
        !blockvectorbx=matmul(blockvectorbx,gramxax)
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
             &               vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
        blockvectorbx=blockvectordumm
        !blockvectorr=blockvectorax-matmul(blockvectorbx,lambda)
        !         call dgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
        !&               vectsize,lambda,blocksize,czero,blockvectordumm,vectsize)
        ! blockvectorr=blockvectorax-blockvectordumm
        do iblocksize=1,blocksize
           blockvectorr(:,iblocksize)=blockvectorax(:,iblocksize)-eigen(iblocksize)*blockvectorbx(:,iblocksize)
        end do
        deallocate(eigen)
     end if !epilogue



     residualnorms=sum(abs(blockvectorr)**2,dim=1)
     old_paral_level= mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm)
     call timab(48,1,tsec)
     call xsum_mpi(residualnorms,spaceComm,ierr)
     call timab(48,2,tsec)
     mpi_enreg%paral_level= old_paral_level
     residualnorms=sqrt(residualnorms)
     do iblocksize=1,blocksize
        iband=iblocksize+(iblock-1)*blocksize
        cg(1,cgindex(iband):cgindex(iband+1)-1)=real(blockvectorx(1:vectsize,iblocksize))
        cg(2,cgindex(iband):cgindex(iband+1)-1)=aimag(blockvectorx(1:vectsize,iblocksize))
     end do
     if(gen_eigenpb) then
        do iblocksize=1,blocksize
           iband=iblocksize+(iblock-1)*blocksize
           gsc(1,gscindex(iband):gscindex(iband+1)-1)=real(blockvectorbx(1:vectsize,iblocksize))
           gsc(2,gscindex(iband):gscindex(iband+1)-1)=aimag(blockvectorbx(1:vectsize,iblocksize))
        end do
     end if
     !this should not exist,since this induce one too much getghc.lazy programming....
     !call operatorh(blockvectorx,blockvectorax,subham,subvnl)!fill also subham, subvnl
     allocate(cwavef(2,npw_k*nspinor*blocksize),gwavef(2,npw_k*nspinor*blocksize),gvnlc(2,npw_k*nspinor*blocksize))
     isubh=1+2*(iblock-1)*blocksize*((iblock-1)*blocksize+1)/2
     do iblocksize=1,blocksize
        cwavef(1,npw_k*nspinor*(iblocksize-1)+1:npw_k*nspinor*iblocksize)=real(blockvectorx(1:npw_k*nspinor,iblocksize))
        cwavef(2,npw_k*nspinor*(iblocksize-1)+1:npw_k*nspinor*iblocksize)=aimag(blockvectorx(1:npw_k*nspinor,iblocksize))
     end do
     tim_getghc=1 ; sij_opt=0
     call getghc(cwavef,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
          & kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,nspinor,ntypat,&
          & nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
     do iblocksize=1,blocksize
        blockvectorax(1:npw_k*nspinor,iblocksize)=&
             &  dcmplx(gwavef(1,npw_k*nspinor*(iblocksize-1)+1:npw_k*nspinor*iblocksize),&
             &         gwavef(2,npw_k*nspinor*(iblocksize-1)+1:npw_k*nspinor*iblocksize))
        do ii=1,(iblock-1)*blocksize+iblocksize
           iwavef=(ii-1)*npw_k*nspinor+icg
           chcre=zero ; chcim=zero
           if (gs_hamk%usepaw==1) then
              do ipw=1,npw_k*nspinor
                 cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
                 chcre=chcre+cgreipw*gwavef(1,ipw+(iblocksize-1)*npw_k*nspinor)+cgimipw*gwavef(2,ipw+(iblocksize-1)*npw_k*nspinor)
                 chcim=chcim+cgreipw*gwavef(2,ipw+(iblocksize-1)*npw_k*nspinor)-cgimipw*gwavef(1,ipw+(iblocksize-1)*npw_k*nspinor)
              end do
           else
              cvcre=zero ; cvcim=zero
              do ipw=1,npw_k*nspinor
                 cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
                 chcre=chcre+cgreipw*gwavef(1,ipw+(iblocksize-1)*npw_k*nspinor)+cgimipw*gwavef(2,ipw+(iblocksize-1)*npw_k*nspinor)
                 chcim=chcim+cgreipw*gwavef(2,ipw+(iblocksize-1)*npw_k*nspinor)-cgimipw*gwavef(1,ipw+(iblocksize-1)*npw_k*nspinor)
                 cvcre=cvcre+cgreipw*gvnlc(1,ipw+(iblocksize-1)*npw_k*nspinor)+cgimipw*gvnlc(2,ipw+(iblocksize-1)*npw_k*nspinor)
                 cvcim=cvcim+cgreipw*gvnlc(2,ipw+(iblocksize-1)*npw_k*nspinor)-cgimipw*gvnlc(1,ipw+(iblocksize-1)*npw_k*nspinor)
              end do
              !   Store real and imag parts in hermitian storage mode:
              subvnl(isubh)=cvcre ; subvnl(isubh+1)=cvcim
           end if
           !  Store real and imag parts in hermitian storage mode:
           subham(isubh)=chcre ; subham(isubh+1)=chcim
           isubh=isubh+2
        end do
     end do
     !comm for subham and subvnl is made in vtowfk

     deallocate(cwavef,gwavef,gvnlc)
     !call operators(blockvectorx,blockvectorbx,subovl)!fill also  subovl
     if((gen_eigenpb).and.(use_subovl==1)) then
        allocate(cwavef(2,npw_k*nspinor))
        allocate(gwavef(2,npw_k*nspinor))
        isubo=1+2*(iblock-1)*blocksize*((iblock-1)*blocksize+1)/2
        do iblocksize=1,blocksize
           cwavef(1,1:npw_k*nspinor)=real (blockvectorx(1:npw_k*nspinor,iblocksize))
           cwavef(2,1:npw_k*nspinor)=aimag(blockvectorx(1:npw_k*nspinor,iblocksize))
           !  Call to nonlop: compute <g|S|c>
           choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
           call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
                &              dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
                &              istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
                &              mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
                &              gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
                &              gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
                &              gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
           blockvectorbx(1:npw_k*nspinor,iblocksize)=dcmplx(gwavef(1,1:npw_k*nspinor),gwavef(2,1:npw_k*nspinor))
           do ii=1,(iblock-1)*blocksize+iblocksize
              iwavef=(ii-1)*npw_k*nspinor+icg
              cscre=zero;cscim=zero
              do ipw=1,npw_k*nspinor
                 cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
                 cscre=cscre+cgreipw*gwavef(1,ipw)+cgimipw*gwavef(2,ipw)
                 cscim=cscim+cgreipw*gwavef(2,ipw)-cgimipw*gwavef(1,ipw)
              end do
              !   Store real and imag parts in hermitian storage mode:
              subovl(isubo)=cscre ; subovl(isubo+1)=cscim
              isubo=isubo+2
           end do
        end do
        deallocate(cwavef,gwavef)
     end if
     write(6,*) "mytest"
     ! stop

     write(6,*)'residualnorm at the end end',residualnorms
     deallocate(blockvectorx,blockvectorax,blockvectorbx)
     deallocate(blockvectorr,blockvectorar,blockvectorbr)
     deallocate(blockvectorp,blockvectorap,blockvectorbp)
     deallocate(blockvectory,blockvectorby)
     deallocate(gramyx)
     deallocate(transf,blockvectordumm,blockvectorxc)
     deallocate(gramxax,gramxar,gramxap,gramrar,gramrap,grampap,gramxbx,gramxbr,&
          & gramxbp,gramrbr,gramrbp,grampbp)
     deallocate(lambda)
     deallocate(residualnorms,pflag)
     !End big loop over bands inside blocks
  end do

  call timab(530,2,tsec)
  !  stop('fin de lobpcgccIIwf')

end subroutine lobpcgcciiwf
!!***

!!****f* ABINIT/lobpcgcciiiwf
!! NAME
!! lobpcgcciiiwf
!!
!! FUNCTION
!!
!! PARENTS
!!      lobpcgccIIwf
!!
!! CHILDREN
!!      getghc,nonlop,zgemm,zhegv,zorthonormalize,zprecon3,ztrsm
!!
!! SOURCE
!!

subroutine lobpcgcciiiwf(cg,dimffnl,dtfil,dtset,ffnl,gs_hamk,gsc,icg,igsc,&
     & kg_k,kinpw,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,&
     & mpssoang,natom,nbdblock,nband_k,npw_k,nspinor,ntypat,&
     & nvloc,n4,n5,n6,ph3d,prtvol,psps,resid_k,use_subovl,vlocal,&
     & subham,subvnl,subovl,blocksize,bblocksize,vectsize,pflag,&
     & blockvectorx,blockvectorbx,blockvectorax,blockvectory,blockvectorby,lambda,&
     & blockvectorp,blockvectorbp,blockvectorap&
     &   )

  use defs_basis
  use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_14wfs
 use interfaces_linalg
!End of the abilint section

  implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
  !arguments ------------------------------------
  integer :: spacecomm=0
  type(gs_hamiltonian_type) :: gs_hamk
  integer :: dimffnl,icg,igsc,ikg,i,lmnmax,matblk,mcg,mgsc,mgfft,mpsang,mpssoang,n4,n5
  integer :: n6,natom,nband_k,nbdblock,npw_k,nspinor,ntypat,nvloc,prtvol
  integer :: use_ffnl,use_subovl
  real(dp) :: sq2,tsec
  type(datafiles_type) :: dtfil
  type(dataset_type) :: dtset
  type(pseudopotential_type) :: psps
  type(mpi_type) :: mpi_enreg
  integer :: kg_k(3,npw_k)
  integer :: bblocksize
  real(dp) :: cg(2,mcg)
  real(dp) :: ffnl(npw_k,dimffnl,lmnmax,ntypat),gsc(2,mgsc)
  real(dp) :: kinpw(npw_k),ph3d(2,npw_k,matblk),resid_k(nband_k)
  real(dp) :: vlocal(n4,n5,n6,nvloc)
  real(dp) :: subham(nband_k*(nband_k+1)),subvnl(nband_k*(nband_k+1)),subovl(nband_k*(nband_k+1))
  real(dp) :: lambda
  integer :: blocksize,vectsize
  complex(dp) :: blockvectorx(vectsize,blocksize),blockvectorax(vectsize,blocksize),&
       &                 blockvectorbx(vectsize,blocksize),blockvectory(vectsize,bblocksize),&
       &                 blockvectorby(vectsize,bblocksize)!,lambda(blocksize,blocksize)
  complex(dp) :: blockvectorp(vectsize,blocksize),blockvectorap(vectsize,blocksize),&
       &                 blockvectorbp(vectsize,blocksize)
  logical :: pflag(blocksize)
  !local variables-------------------------------
  !  integer:: ii,ipw,ipw1,ivectsize,isubo,isubh,iwavef,jwavef,maxiterations
  !  integer:: bigorder,info,lwork,tim_getghc
  !  integer:: iblocksize, activersize, activepsize, restart, cond_try, i1, i2, i3, i4
  !  integer:: rvectsize,vectsize,blocksize,bblocksize,iband,istwf_k,iblock,cgindex
  integer:: activepsize,activersize,bigorder,cgindex,choice,cond_try,cpopt,iblock,iblocksize
  integer:: i1,i2,i3,i4,iband,idir,ier,ierr,ii,info
  integer:: ipw,ipw1,istwf_k,isubo,isubh,ivectsize,iwavef,jblocksize,jwavef,lwork
  integer:: maxiterations,nnlout,nkpg,old_paral_level,optekin,paw_opt,restart,signs,sij_opt,tim_getghc,tim_nonlop
  logical::gen_eigenpb
  real(dp) :: cgreipw,cgimipw,cscre,cscim
  real(dp) :: chcre, chcim,cvcre, cvcim,dum
  complex(dp) ::  blockvectorz(vectsize,blocksize),blockvectoraz(vectsize,blocksize),&
       &                 blockvectorbz(vectsize,blocksize),blockvectorr1(vectsize,blocksize)
  real(dp) :: residualnorms1(blocksize)
  real(dp), allocatable :: dummy1(:),dummy2(:,:),dummy3(:,:,:)
  complex(dp) :: zlambda(1,1)
  type(cprj_type) :: cprj_dum(1,1)
  !local variables turned arguments--------------

  complex(dp), allocatable ::  blockvectorr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
       & blockvectordumm(:,:),&
       & gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
       & grampap(:,:),&
       & gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
       & grampbp(:,:),&
       & identity(:,:),coordx(:,:),&
       & grama(:,:),gramb(:,:),gramyx(:,:),&
       & w(:),work(:)
  real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:)
  real(dp), allocatable :: residualnorms(:),eigen(:),rwork(:),kpg_dum(:,:)
  real(dp), parameter :: tolerance1 = 1.e-13, tolerance2 = 1.e2

#ifdef VMS
  !DEC$ ATTRIBUTES ALIAS:'ZHEGV' :: zhegv
  !DEC$ ATTRIBUTES ALIAS:'ZTRSM' :: ztrsm
  !DEC$ ATTRIBUTES ALIAS:'ZGEMM' :: zgemm
#endif

  !correspondence with abinit. here for real wf but in complex mode
  !this is the index of a given band
  !  cgindex(iblocksize)=npw_k*nspinor*(iblocksize-1)+icg+1
  gen_eigenpb=(gs_hamk%usepaw==1)
  optekin=0;if (dtset%wfoptalg>10) optekin=1
  sq2=sqrt(two)
  !vectsize=npw_k*nspinor
  !blocksize=(nband_k-1)/nbdblock+1
  !bblocksize=(iblock-1)*blocksize
  istwf_k=gs_hamk%istwf_k
  maxiterations=dtset%nline

  !passing x into z
  blockvectorz = blockvectorx
  blockvectoraz = blockvectorax
  blockvectorbz = blockvectorbx

  !allocations
  allocate(blockvectorr(vectsize,blocksize),blockvectorar(vectsize,blocksize))
  allocate(blockvectorbr(vectsize,blocksize))
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
  allocate(residualnorms(blocksize))
  !     write(6,*)'iterationnumber',iterationnumber
  !construct residual
  !blockvectorr=blockvectorax-matmul(blockvectorx,lambda)
  do iblocksize=1,blocksize
     zlambda(1,1)=lambda
     call zprecon3(blockvectorbx(:,iblocksize),zlambda,1,&!,iblocksize),&
          &                 istwf_k,kinpw,mpi_enreg,npw_k,nspinor,&
          &                 optekin,blockvectorax(:,iblocksize),blockvectorr(:,iblocksize),vectsize)
     !      blockvectorr(:,iblocksize)=blockvectorax(:,iblocksize)-lambda(iblocksize,iblocksize)*blockvectorbx(:,iblocksize)
  end do
  residualnorms=sqrt(sum(abs(blockvectorr)**2,dim=1))
  write(6,*)'residualnorm',residualnorms
  !if(abs(sum(residualnorms)) < 1.d-10) exit
  !not yet masking

  !DEBUG le if suivant ne marche que si blocksize = 1
  if (residualnorms(1)>tolerance1) then !sinon on masque

     if(bblocksize>0) then !residuals orthogonal to blockvectorby
        !  blockvectorr=blockvectorr-&
        !           &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorr))
        call zgemm('c','n',bblocksize,blocksize,vectsize,cone,blockvectorby,&
             &               vectsize,blockvectorr,vectsize,czero,gramyx,bblocksize)
        call zgemm('n','n',vectsize,blocksize,bblocksize,cone,blockvectory,&
             &               vectsize,gramyx,bblocksize,czero,blockvectordumm,vectsize)
        blockvectorr=blockvectorr-blockvectordumm
     end if
     !residuals orthogonal to blockvectorx
     !     blockvectorr=blockvectorr-&
     !          &matmul(blockvectorx,matmul(transpose(blockvectorbx),blockvectorr))
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
          &               vectsize,blockvectorr,vectsize,czero,gramxax,blocksize)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
          &               vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
     blockvectorr=blockvectorr-blockvectordumm
     !and now (b)orthornormalize r
     !if(iterationnumber >1) stop('xx2')
     !call operators(blockvectorr,blockvectorbr)
     allocate(cwavef(2,npw_k*nspinor))
     allocate(gwavef(2,npw_k*nspinor))
     do iblocksize=1,blocksize
        cwavef(1,:)=real (blockvectorr(1:npw_k*nspinor,iblocksize))
        cwavef(2,:)=aimag(blockvectorr(1:npw_k*nspinor,iblocksize))
        if(gen_eigenpb) then
           !   Call to nonlop: compute <g|S|c>
           choice=1 ; signs=2 ; idir=0 ; tim_nonlop=1 ; cpopt=-1 ; paw_opt=3 ; nnlout=0 ; nkpg=0
           call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
                &               dummy1,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
                &               istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
                &               mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
                &               gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
                &               gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
                &               gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)
        else
           gwavef(:,:)=cwavef(:,:)
        end if
        blockvectorbr(1:npw_k*nspinor,iblocksize)=dcmplx(gwavef(1,:),gwavef(2,:))
     end do
     deallocate(cwavef,gwavef)
     !call zorthonormalize(blockvectorr,blockvectorbr)
     call zorthonormalize(blockvectorr,blockvectorbr,blocksize,mpi_enreg,gramrbr,vectsize)
     call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramrbr,blocksize,&
          &              blockvectorbr,vectsize)
     !compute ar
     !blockvectorar=matmul(operatora,blockvectorr)
     !if(iterationnumber >1) stop('xx3')
     !call operatorh(blockvectorr,blockvectorar)
     allocate(cwavef(2,npw_k*nspinor),gwavef(2,npw_k*nspinor),gvnlc(2,npw_k*nspinor))
     do iblocksize=1,blocksize
        cwavef(1,:)=real(blockvectorr(1:npw_k*nspinor,iblocksize))
        cwavef(2,:)=aimag(blockvectorr(1:npw_k*nspinor,iblocksize))
        tim_getghc=1; sij_opt=0
        call getghc(cwavef,dimffnl,ffnl,dtfil%filstat,gwavef,dummy2,gs_hamk,gvnlc,kg_k,&
             &   kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,nspinor,ntypat,&
             &   nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
        blockvectorar(1:npw_k*nspinor,iblocksize)=dcmplx(gwavef(1,:),gwavef(2,:))
     end do
     deallocate(cwavef,gwavef,gvnlc)
     if(pflag(1)) then  !ce n'est pas la bonne condition si blocksize /= 1
        !write(6,*)'blockvectorp,blockvectorbp,blockvectorap'
        !write(6,*)blockvectorp
        !write(6,*)blockvectorbp
        !write(6,*)blockvectorap
        !if(iterationnumber >1) stop('xx4')
        !call zorthonormalize(blockvectorp,blockvectorbp,blockvectorap)
        call zorthonormalize(blockvectorp,blockvectorbp,blocksize,mpi_enreg,grampbp,vectsize)
        call ztrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,&
             &              blockvectorbp,vectsize)
        call ztrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,&
             &              blockvectorap,vectsize)

        !if(iterationnumber >1) stop('xx4')
        !blockvectorap=matmul(blockvectorap,grampbp)
     end if
     activersize=blocksize
     if (.not.pflag(1)) then   !ce n'est pas la bonne condition si blocksize /= 1
        !  activepsize=0
        restart=1
     else
        !  activepsize=blocksize
        restart=0
     end if
     !gramxar=matmul(transpose(blockvectorax),blockvectorr)
     !gramrar=matmul(transpose(blockvectorar),blockvectorr)
     !gramxax=matmul(transpose(blockvectorax),blockvectorx)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
          &               vectsize,blockvectorr,vectsize,czero,gramxar,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorar,&
          &               vectsize,blockvectorr,vectsize,czero,gramrar,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
          &               vectsize,blockvectorx,vectsize,czero,gramxax,blocksize)
     !gramxbx=matmul(transpose(blockvectorbx),blockvectorx)
     !gramrbr=matmul(transpose(blockvectorbr),blockvectorr)
     !gramxbr=matmul(transpose(blockvectorbx),blockvectorr)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
          &               vectsize,blockvectorx,vectsize,czero,gramxbx,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbr,&
          &               vectsize,blockvectorr,vectsize,czero,gramrbr,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
          &               vectsize,blockvectorr,vectsize,czero,gramxbr,blocksize)
     !if(iterationnumber >1) stop('xx1')
     i1=0;i2=blocksize;i3=2*blocksize;i4=3*blocksize
     cond: do cond_try=1,2 !2 when restart, but not implemented
        if (restart==0) then
           !gramxap=matmul(transpose(blockvectorax),blockvectorp)
           !gramrap=matmul(transpose(blockvectorar),blockvectorp)
           !grampap=matmul(transpose(blockvectorap),blockvectorp)
           call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
                &               vectsize,blockvectorp,vectsize,czero,gramxap,blocksize)
           call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorar,&
                &               vectsize,blockvectorp,vectsize,czero,gramrap,blocksize)
           call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorap,&
                &               vectsize,blockvectorp,vectsize,czero,grampap,blocksize)
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
           !grampbp=matmul(transpose(blockvectorbp),blockvectorp)
           call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
                &               vectsize,blockvectorp,vectsize,czero,gramxbp,blocksize)
           call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbr,&
                &               vectsize,blockvectorp,vectsize,czero,gramrbp,blocksize)
           call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbp,&
                &               vectsize,blockvectorp,vectsize,czero,grampbp,blocksize)

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
!!$     end do cond
        !call la_sygv(grama,gramb,eigen,itype=1,jobz='v')
        !if(iterationnumber >1) stop('xx')
        lwork=3*bigorder-2
        allocate(work(lwork),rwork(lwork))
        call zhegv(1,'v','u',bigorder,grama,bigorder,gramb,bigorder,eigen,&
             &               work,lwork,rwork,info)
        deallocate(work,rwork)
        do iblocksize=1,blocksize
           lambda=eigen(iblocksize)
        end do
        write(6,*)'eigen',eigen(1:blocksize)
        coordx=grama(:,1:blocksize)
        deallocate(grama,gramb,eigen)
        if (restart==0) then
           !        blockvectorp=matmul(blockvectorr,coordx(i2+1:i3,:))+&
           !             &matmul(blockvectorp,coordx(i3+1:i4,:))
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorr,&
                &               vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorp,&
                &               vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
           blockvectorp=blockvectordumm
           !        blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))+&
           !             &matmul(blockvectorap,coordx(i3+1:i4,:))
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorar,&
                &               vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorap,&
                &               vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
           blockvectorap=blockvectordumm
           !        blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))+&
           !             &matmul(blockvectorbp,coordx(i3+1:i4,:))
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbr,&
                &               vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbp,&
                &               vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
           blockvectorbp=blockvectordumm
        else
           !blockvectorp =matmul(blockvectorr,coordx(i2+1:i3,:))
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorr,&
                &               vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorp,vectsize)
           !blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorar,&
                &               vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorap,vectsize)
           !blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))
           call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbr,&
                &               vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorbp,vectsize)
        end if

        !blockvectorx = matmul(blockvectorx,coordx(i1+1:i2,:))+blockvectorp
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
             &               vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
        blockvectorx = blockvectordumm+blockvectorp
        !blockvectorax= matmul(blockvectorax,coordx(i1+1:i2,:))+blockvectorap
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
             &               vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
        blockvectorax = blockvectordumm+blockvectorap
        !blockvectorbx= matmul(blockvectorbx,coordx(i1+1:i2,:))+blockvectorbp
        call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
             &               vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
        blockvectorbx = blockvectordumm+blockvectorbp
        deallocate(coordx)

        do iblocksize=1,blocksize
           zlambda(1,1) = lambda
           call zprecon3(blockvectorbx(:,iblocksize),zlambda,1,& !iblocksize),&
                &                 istwf_k,kinpw,mpi_enreg,npw_k,nspinor,&
                &                 optekin,blockvectorax(:,iblocksize),blockvectorr1(:,iblocksize),vectsize)
        end do
        residualnorms1=sqrt(sum(abs(blockvectorr1)**2,dim=1))
        write(6,*)'residualnorm after lobpcgcciii',residualnorms1

        do iblocksize=1,blocksize
           if (residualnorms1(iblocksize) > tolerance2*residualnorms(iblocksize)) then
              blockvectorx = blockvectorz
              blockvectorax = blockvectoraz
              blockvectorbx = blockvectorbz
              if (cond_try == 1) then
                 write(6,*)'restart here for this eig'
                 restart = 1
              else
                 blockvectorp = czero
                 blockvectorap = czero
                 blockvectorbp = czero
                 pflag = .false.
                 write(6,*)'restart was unuseful'
              end if
           else
              pflag = .true.
              exit cond   !DEBUG exact que si blocksize = 1 (sinon il faut tester les autres )
           end if
        end do
     end do cond

  else
     blockvectorp = czero
     blockvectorap = czero
     blockvectorbp = czero
     pflag = .false.
  end if

  !write(6,*)'blockvectorr',blockvectorr
  !write(6,*)'blockvectorx',blockvectorx
  !write(6,*)'blockvectorax',blockvectorax
  deallocate(blockvectorr,blockvectorar,blockvectorbr)
  deallocate(gramyx)
  deallocate(blockvectordumm)
  deallocate(gramxax,gramxar,gramxap,gramrar,gramrap,grampap,gramxbx,gramxbr,&
       & gramxbp,gramrbr,gramrbp,grampbp)
  deallocate(residualnorms)
end subroutine lobpcgcciiiwf
!!***
