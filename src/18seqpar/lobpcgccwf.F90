!{\src2tex{textfont=tt}}
!!****f* abinit/lobpcgccwf
!! NAME
!! lobpcgccwf
!!
!! FUNCTION
!! this routine updates the whole wave functions at a given k-point,
!! using the lobpcg method
!! for a given spin-polarization, from a fixed hamiltonian
!! but might also simply compute eigenvectors and eigenvalues at this k point.
!! it will also update the matrix elements of the hamiltonian.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FBottin,GZ,AR,MT)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variales for this dataset
!!  ffnl(npw_k,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
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
!!  nkpg=second dimension of kpg_k (0 if useylm=0)
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
!!    totvnl(2*nband_k*(1-gs_hamk%usepaw),nband_k*(1-gs_hamk%usepaw))=the matrix elements of vnl
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
!!      getghc,nonlop,prep_getghc,timab,wrtout,xcomm_init,xsum_mpi,zcopy,zgemm
!!      zheev,zhegv,zorthonormalize,zprecon3,ztrsm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine lobpcgccwf(cg,dimffnl,dtfil,dtset,ffnl,ffnl_gather,gs_hamk,gsc,icg,igsc,&
&           kg_k,kg_k_gather,kinpw,kinpw_gather,lmnmax,matblk,mcg,mgfft,mgsc,mpi_enreg,mpsang,mpssoang,natom,&
&           nband_k,nbdblock,ndatarecv,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,ph3d,ph3d_gather,prtvol,&
&           psps,resid_k,subham,subovl,totvnl,usebandfft,use_subovl,vlocal)

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
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(gs_hamiltonian_type) :: gs_hamk
 integer :: dimffnl,icg,igsc,lmnmax,matblk,mcg,mgsc,mgfft,mpsang,mpssoang,n4,n5,n6
 integer :: natom,nband_k,nbdblock,npw_k,nspinor,ntypat,nvloc,prtvol,use_subovl
 type(datafiles_type) :: dtfil
 type(dataset_type) :: dtset
 type(pseudopotential_type) :: psps
 type(mpi_type) :: mpi_enreg
 integer :: kg_k(3,npw_k)
 real(dp) :: cg(2,mcg),ffnl(npw_k,dimffnl,lmnmax,ntypat),gsc(2,mgsc)
 real(dp) :: kinpw(npw_k),ph3d(2,npw_k,matblk),resid_k(nband_k)
 real(dp) :: vlocal(n4,n5,n6,nvloc)
 real(dp) :: subham(nband_k*(nband_k+1))
 real(dp) :: subovl(nband_k*(nband_k+1)*use_subovl)
 real(dp) :: totvnl(2*nband_k*(1-gs_hamk%usepaw),nband_k*(1-gs_hamk%usepaw))
!FB Arguments for the band_fft parallelization
 integer :: ndatarecv,usebandfft
 integer, intent(in) :: kg_k_gather(3,ndatarecv*usebandfft)
 real(dp), intent(in) :: ffnl_gather(ndatarecv,dimffnl,lmnmax,ntypat*usebandfft)
 real(dp), intent(in) :: kinpw_gather(ndatarecv*usebandfft)
 real(dp), intent(in) :: ph3d_gather(2,ndatarecv,matblk*usebandfft)

!Local variables-------------------------------
 integer, parameter :: tim_getghc=5
 integer :: activepsize,activersize,bblocksize,bigorder,blocksize,cgindex,choice,cpopt,windex
 integer :: cond_try,gscindex
 integer :: iband,iblocksize,iblock,idir,ier,ierr,ii,info,ioption,ipw,ipw1,istwf_k,isubh,isubo
 integer :: iterationnumber,ivectsize
 integer :: iwavef,i1,i2,i3,i4,jwavef,lwork,maxblocksize,maxiterations,nkpg,nnlout
 integer :: old_paral_level,optekin,paw_opt,restart
 integer :: signs,sij_opt,spaceComm,spaceComm_keep,tim_nonlop,tocomplete,vectsize
 logical :: gen_eigenpb,block_to_complete
 real(dp) :: cgreipw,cgimipw,cscre,cscim,chcre,chcim,cvcre,cvcim,deltae,deold,dum,sq2
 complex(dpc) :: cminusone
 real(dp) :: tsec(2)
 real(dp), allocatable :: gwavef(:,:),cwavef(:,:),gvnlc(:,:),lambda_loc(:),kpg_dum(:,:)
 real(dp), allocatable :: swavef(:,:)
 real(dp), allocatable :: residualnorms(:),eigen(:),rwork(:)
 real(dp),allocatable :: enlout(:)
 real(dp),allocatable :: sij_loc(:,:),gsc_loc(:,:),kpg_loc(:,:)
 complex(dpc), allocatable :: blockvectorx(:,:),blockvectorvx(:,:),blockvectorax(:,:),blockvectorbx(:,:),&
       & blockvectorr(:,:),blockvectorvr(:,:),blockvectorar(:,:),blockvectorbr(:,:),&
       & blockvectorp(:,:),blockvectorvp(:,:),blockvectorap(:,:),blockvectorbp(:,:),blockvectordumm(:,:),&
       & blockvectory(:,:),blockvectorby(:,:),blockvectorz(:,:),&
       & gramxax(:,:),gramxar(:,:),gramxap(:,:),gramrar(:,:),gramrap(:,:),&
       & grampap(:,:),&
       & gramxbx(:,:),gramxbr(:,:),gramxbp(:,:),gramrbr(:,:),gramrbp(:,:),&
       & grampbp(:,:),&
       & identity(:,:),coordx(:,:),lambda(:,:),&
       & grama(:,:),gramb(:,:),gramyx(:,:),&
       & transf(:,:,:),work(:)
 real(dp), allocatable :: dummy3(:,:,:)
 real(dp),allocatable :: tsubham(:,:)
 type(cprj_type) :: cprj_dum(1,1)
 character(len=500) :: message

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZGEMM' :: zgemm
!DEC$ ATTRIBUTES ALIAS:'ZHEEV' :: zheev
!DEC$ ATTRIBUTES ALIAS:'ZHEGV' :: zhegv
!DEC$ ATTRIBUTES ALIAS:'ZTRSM' :: ztrsm
!DEC$ ATTRIBUTES ALIAS:'ZCOPY' :: zcopy
#endif

!NO_ABIRULES
!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band
!  cgindex(iblocksize)=npw_k*nspinor*(iblocksize-1)+icg+1
!  gscindex(iblocksize)=npw_k*nspinor*(iblocksize-1)+igsc+1
!  windex(iblocksize)=npw_k*nspinor*(iblocksize-1)+1

  optekin=0;if (dtset%wfoptalg>10) optekin=1

! *********************************************************************

!DEBUG
! write(6,*)' lobpcgccwf : enter, icg= ',icg
!ENDDEBUG

 call timab(530,1,tsec)
 gen_eigenpb=(gs_hamk%usepaw==1)
 cminusone=-cone
 sq2=sqrt(2.0_dp)
 vectsize=npw_k*nspinor

 istwf_k=gs_hamk%istwf_k
 maxiterations=dtset%nline
 maxblocksize=mpi_enreg%nproc_fft

 if(mpi_enreg%mode_para=='b') maxblocksize=mpi_enreg%nproc_band*mpi_enreg%bandpp

 ioption=mpi_enreg%fft_option_lob
 !modif :yg 
 !if (mpi_enreg%nproc_fft ==1 .and. mpi_enreg%fft_option_lob==2 .and. mpi_enreg%nproc_band==1) ioption=1

!###########################################################################
!################ BIG LOOP OVER BLOCKS  ####################################
!###########################################################################
 do iblock=1,nbdblock

  blocksize=mpi_enreg%nproc_fft

  if(mpi_enreg%mode_para=='b') blocksize=mpi_enreg%nproc_band*mpi_enreg%bandpp

! block_to_complete : true if the tabs with dimension blocksize need zeroes at the end
  block_to_complete=.false.;tocomplete=0
  if ((iblock == nbdblock) .and. (nbdblock * blocksize > nband_k)) then
! the last block is smaller than the others
   block_to_complete=.true.
   tocomplete=nbdblock * blocksize - nband_k-1
   blocksize=mpi_enreg%nproc_fft - (nbdblock * blocksize - nband_k)
  end if

  bblocksize=(iblock-1)*blocksize

! allocations
  allocate(blockvectorx(vectsize,blocksize),blockvectorax(vectsize,blocksize))
  allocate(blockvectorbx(vectsize,blocksize))
  allocate(blockvectorr(vectsize,blocksize),blockvectorar(vectsize,blocksize))
  allocate(blockvectorbr(vectsize,blocksize))
  allocate(blockvectorp(vectsize,blocksize),blockvectorap(vectsize,blocksize))
  allocate(blockvectorbp(vectsize,blocksize))
  allocate(blockvectordumm(vectsize,blocksize))
  allocate(gramxax(blocksize,blocksize),&
       & gramxar(blocksize,blocksize),gramxap(blocksize,blocksize),&
       & gramrar(blocksize,blocksize),gramrap(blocksize,blocksize),&
       & grampap(blocksize,blocksize),&
       & gramxbx(blocksize,blocksize),gramxbr(blocksize,blocksize),&
       & gramxbp(blocksize,blocksize),gramrbr(blocksize,blocksize),&
       & gramrbp(blocksize,blocksize),&
       & grampbp(blocksize,blocksize))
  allocate(transf(blocksize,blocksize,3))
  allocate(lambda(blocksize,blocksize))
  allocate(residualnorms(blocksize))

  allocate(blockvectory(vectsize,bblocksize),blockvectorby(vectsize,bblocksize))
  allocate(gramyx(bblocksize,blocksize))
  if (gs_hamk%usepaw==0) then
   allocate(blockvectorvx(vectsize,blocksize),blockvectorvr(vectsize,blocksize),blockvectorvp(vectsize,blocksize))
  end if

! transfer array of wf coeff in iblock to blockvectorx (complex to complex)
  if(abs(dtset%timopt)==3)then
   call timab(584,1,tsec)
  endif
  call zcopy(blocksize*vectsize,&
&  cg(:,vectsize*((iblock-1)*blocksize)+icg+1:vectsize*(iblock*blocksize-1)+icg),&
&  1,blockvectorx(:,1:blocksize),1)
  if(abs(dtset%timopt)==3)then
   call timab(584,2,tsec)
  endif

!!!!!!!!!!!!!!!!!!!!!!!!! Begin if iblock /=1 !!!!!!!!!!!!!!!!!!!!!!!!!!
! transfer array of wf coeff less than iblock to blockvectory
  if(iblock /=1) then
   if(abs(dtset%timopt)==3)then
    call timab(584,1,tsec)
   endif
   call zcopy(bblocksize*vectsize,cg(:,icg+1:vectsize*((iblock-1)*blocksize-1)+icg),1,blockvectory(:,1:bblocksize),1)
   if(gen_eigenpb) then
    call zcopy(bblocksize*vectsize,gsc(:,igsc+1:vectsize*((iblock-1)*blocksize-1)+igsc),1,blockvectorby(:,1:bblocksize),1)
   else
    call zcopy(vectsize*bblocksize,blockvectory,1,blockvectorby,1)
   end if
   if(abs(dtset%timopt)==3)then
    call timab(584,2,tsec)
   endif

!  b-orthogonalize x to the constraint y (supposed b-orthonormal)
!  blockvectorx=blockvectorx-&
!              &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorx))

   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
   call zgemm('c','n',bblocksize,blocksize,vectsize,cone,blockvectorby,&
&            vectsize,blockvectorx,vectsize,czero,gramyx,bblocksize)
   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif

   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
   if(abs(dtset%timopt)==3)then
    call timab(533,1,tsec)
   endif
   call xsum_mpi(gramyx,spaceComm,ierr)
   if(abs(dtset%timopt)==3)then
    call timab(533,2,tsec)
   endif
   mpi_enreg%paral_level= old_paral_level

   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
   call zgemm('n','n',vectsize,blocksize,bblocksize,cminusone,blockvectory,&
&            vectsize,gramyx,bblocksize,cone,blockvectorx,vectsize)
   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif
  end if
!!!!!!!!!!!!!!!!!!!!!!!!! End if iblock /=1 !!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(cwavef(2,npw_k*nspinor*maxblocksize))
  allocate(gwavef(2,npw_k*nspinor*maxblocksize),gvnlc(2,npw_k*nspinor*maxblocksize))
  allocate(swavef(2,npw_k*nspinor*maxblocksize))
  if (block_to_complete) then
   cwavef(:,:)=zero
  end if
  if(abs(dtset%timopt)==3)then
   call timab(584,1,tsec)
  endif
  call zcopy(vectsize*blocksize,blockvectorx,1,cwavef,1)
  if(abs(dtset%timopt)==3)then
   call timab(584,2,tsec)
  endif

  sij_opt=0;if (gen_eigenpb) sij_opt=1
  if (ioption==1) then
   call getghc(cwavef,dimffnl,ffnl,dtfil%filstat,gwavef,swavef,gs_hamk,gvnlc,kg_k,&
&   kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,nspinor,ntypat,&
&   nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
  else
   call timab(534,1,tsec)
   call prep_getghc(cwavef,dimffnl,dtfil,ffnl_gather,gs_hamk,gvnlc,gwavef,swavef,iblock,1,istwf_k,kg_k_gather,&
&   kinpw_gather,dum,lmnmax,matblk,maxblocksize,mgfft,mpi_enreg,mpsang,mpssoang,natom,nbdblock,&
&   nband_k,ndatarecv,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,dtset%paral_kgb,ph3d_gather,prtvol,sij_opt,vlocal)
   call timab(534,2,tsec)
  end if
  if(abs(dtset%timopt)==3)then
   call timab(584,1,tsec)
  endif
  if (gen_eigenpb) then
   call zcopy(vectsize*blocksize,swavef,1,blockvectorbx,1)
  else
   call zcopy(vectsize*blocksize,blockvectorx,1,blockvectorbx,1)
   call zcopy(vectsize*blocksize,gvnlc,1,blockvectorvx,1)
  endif
  call zcopy(vectsize*blocksize,gwavef,1,blockvectorax,1)
  if(abs(dtset%timopt)==3)then
   call timab(584,2,tsec)
  endif
  deallocate(cwavef,gwavef,gvnlc)
  deallocate(swavef)
  if(abs(dtset%timopt)==3)then
   call timab(535,1,tsec)
  endif
  call zorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,gramxbx,vectsize)
  call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramxbx,blocksize,&
&            blockvectorbx,vectsize)
  call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramxbx,blocksize,&
&            blockvectorax,vectsize)
  if (gs_hamk%usepaw==0) then
   call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramxbx,blocksize,&
&             blockvectorvx,vectsize)
  end if
  if(abs(dtset%timopt)==3)then
   call timab(535,2,tsec)
  endif

! do rayleigh ritz on a in space x
! gramxax=matmul(transpose(blockvectorx),blockvectorax)
  if(abs(dtset%timopt)==3)then
   call timab(532,1,tsec)
  endif
  call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
&            vectsize,blockvectorax,vectsize,czero,gramxax,blocksize)
  if(abs(dtset%timopt)==3)then
   call timab(532,2,tsec)
  endif

  old_paral_level= mpi_enreg%paral_level
  mpi_enreg%paral_level=3
  call xcomm_init(mpi_enreg,spaceComm)
  if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
  if(abs(dtset%timopt)==3)then
   call timab(533,1,tsec)
  endif
  call xsum_mpi(gramxax,spaceComm,ierr)
  if(abs(dtset%timopt)==3)then
   call timab(533,2,tsec)
  endif
  mpi_enreg%paral_level= old_paral_level
  allocate(eigen(blocksize))
  lwork=3*blocksize-2
  allocate(work(lwork),rwork(lwork))
  if(abs(dtset%timopt)==3)then
   call timab(587,1,tsec)
  endif
  call zheev('v','u',blocksize,gramxax,blocksize,eigen,work,lwork,rwork,info)
  if(abs(dtset%timopt)==3)then
   call timab(587,2,tsec)
  endif
  deallocate(work,rwork)

! blockvectorx=matmul(blockvectorx,gramxax)
  if(abs(dtset%timopt)==3)then
   call timab(532,1,tsec)
  endif
  call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&            vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
  if(abs(dtset%timopt)==3)then
   call timab(532,2,tsec)
  endif
  if(abs(dtset%timopt)==3)then
   call timab(584,1,tsec)
  endif
  call zcopy(vectsize*blocksize,blockvectordumm,1,blockvectorx,1)
  if(abs(dtset%timopt)==3)then
   call timab(584,2,tsec)
  endif

! blockvectorax=matmul(blockvectorax,gramxax)
  if(abs(dtset%timopt)==3)then
   call timab(532,1,tsec)
  endif
  call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&            vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
  if(abs(dtset%timopt)==3)then
   call timab(532,2,tsec)
  endif
  if(abs(dtset%timopt)==3)then
   call timab(584,1,tsec)
  endif
  call zcopy(vectsize*blocksize,blockvectordumm,1,blockvectorax,1)
  if(abs(dtset%timopt)==3)then
   call timab(584,2,tsec)
  endif

! blockvectorvx=matmul(blockvectorvx,gramxax)
  if (gs_hamk%usepaw==0) then
   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
   call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvx,&
&             vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif
   if(abs(dtset%timopt)==3)then
    call timab(584,1,tsec)
   endif
   call zcopy(vectsize*blocksize,blockvectordumm,1,blockvectorvx,1)
   if(abs(dtset%timopt)==3)then
    call timab(584,2,tsec)
   endif
  end if 

! blockvectorbx=matmul(blockvectorbx,gramxax)
   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
  call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&            vectsize,gramxax,blocksize,czero,blockvectordumm,vectsize)
  if(abs(dtset%timopt)==3)then
   call timab(532,2,tsec)
  endif
  if(abs(dtset%timopt)==3)then
   call timab(584,1,tsec)
  endif
  call zcopy(vectsize*blocksize,blockvectordumm,1,blockvectorbx,1)
  if(abs(dtset%timopt)==3)then
   call timab(584,2,tsec)
  endif

  do iblocksize=1,blocksize
   lambda(iblocksize,iblocksize)=eigen(iblocksize)
  end do
  deallocate(eigen)

!###########################################################################
!################ PERFORM LOOP ON NLINE ####################################
!###########################################################################
! now the main alogrithm
  iter: do iterationnumber=1,maxiterations
!  construct residual
!  blockvectorr=blockvectorax-matmul(blockvectorx,lambda)
   call zprecon3(blockvectorbx,lambda,blocksize,&
&                istwf_k,kinpw,mpi_enreg,npw_k,nspinor,&
&                optekin,blockvectorax,blockvectorr,vectsize)
   residualnorms=sum(abs(blockvectorr)**2,dim=1)
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
   if(abs(dtset%timopt)==3)then
    call timab(533,1,tsec)
   endif
   call xsum_mpi(residualnorms,spaceComm,ierr)
   if(abs(dtset%timopt)==3)then
    call timab(533,2,tsec)
   endif
   mpi_enreg%paral_level= old_paral_level
   resid_k(bblocksize+1:bblocksize+blocksize)=residualnorms(1:blocksize)

!  If residual sufficiently small stop line minimizations
   if (maxval(residualnorms(1:blocksize))<dtset%tolwfr) then
    if(prtvol>=10)then
     write(message, '(a,i4,a,i2,a,es12.4)' ) &
&        ' lobpcgccwf: block',iblock,' converged after ',iterationnumber,&
&        ' line minimizations : maxval(resid(1:blocksize)) =',maxval(residualnorms(1:blocksize))
     call wrtout(06,message,'PERS')
    end if
    exit
   end if
   if(iblock /=1) then !residuals orthogonal to blockvectorby
!   blockvectorr=blockvectorr-&
!           &matmul(blockvectory,matmul(transpose(blockvectorby),blockvectorr))
    if(abs(dtset%timopt)==3)then
     call timab(532,1,tsec)
    endif
    call zgemm('c','n',bblocksize,blocksize,vectsize,cone,blockvectorby,&
&              vectsize,blockvectorr,vectsize,czero,gramyx,bblocksize)
    if(abs(dtset%timopt)==3)then
     call timab(532,2,tsec)
    endif
    old_paral_level= mpi_enreg%paral_level
    mpi_enreg%paral_level=3
    call xcomm_init(mpi_enreg,spaceComm)
    if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
    if(abs(dtset%timopt)==3)then
     call timab(533,1,tsec)
    endif
    call xsum_mpi(gramyx,spaceComm,ierr)
    if(abs(dtset%timopt)==3)then
     call timab(533,2,tsec)
    endif
    mpi_enreg%paral_level= old_paral_level
    if(abs(dtset%timopt)==3)then
     call timab(532,1,tsec)
    endif
    call zgemm('n','n',vectsize,blocksize,bblocksize,cminusone,blockvectory,&
&              vectsize,gramyx,bblocksize,cone,blockvectorr,vectsize)
    if(abs(dtset%timopt)==3)then
     call timab(532,2,tsec)
    endif
   end if

!  residuals orthogonal to blockvectorx
!  blockvectorr=blockvectorr-&
!          &matmul(blockvectorx,matmul(transpose(blockvectorbx),blockvectorr))
   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&             vectsize,blockvectorr,vectsize,czero,gramxax,blocksize)
   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
   if(abs(dtset%timopt)==3)then
    call timab(533,1,tsec)
   endif
   call xsum_mpi(gramxax,spaceComm,ierr)
   if(abs(dtset%timopt)==3)then
    call timab(533,2,tsec)
   endif
   mpi_enreg%paral_level= old_paral_level
   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
   call zgemm('n','n',vectsize,blocksize,blocksize,cminusone,blockvectorx,&
&             vectsize,gramxax,blocksize,cone,blockvectorr,vectsize)
   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif

   allocate(cwavef(2,npw_k*nspinor*maxblocksize))
   allocate(gwavef(2,npw_k*nspinor*maxblocksize),gvnlc(2,npw_k*nspinor*maxblocksize))
   allocate(swavef(2,npw_k*nspinor*maxblocksize))
   if (block_to_complete) then
    cwavef(:,:)=zero
   end if
   if(abs(dtset%timopt)==3)then
    call timab(584,1,tsec)
   endif
   call zcopy(vectsize*blocksize,blockvectorr,1,cwavef,1)
   if(abs(dtset%timopt)==3)then
    call timab(584,2,tsec)
   endif
   sij_opt=0;if (gen_eigenpb) sij_opt=1
   if (ioption==1) then
    call getghc(cwavef,dimffnl,ffnl,dtfil%filstat,gwavef,swavef,gs_hamk,gvnlc,kg_k,&
&    kinpw,dum,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,blocksize,npw_k,nspinor,ntypat,&
&    nvloc,n4,n5,n6,dtset%paral_kgb,ph3d,prtvol,sij_opt,tim_getghc,0,vlocal)
   else
    call timab(534,1,tsec)
    call prep_getghc(cwavef,dimffnl,dtfil,ffnl_gather,gs_hamk,gvnlc,gwavef,swavef,iblock,2,istwf_k,kg_k_gather,&
&    kinpw_gather,dum,lmnmax,matblk,maxblocksize,mgfft,mpi_enreg,mpsang,mpssoang,natom,nbdblock,&
&    nband_k,ndatarecv,npw_k,nspinor,ntypat,nvloc,n4,n5,n6,dtset%paral_kgb,ph3d_gather,prtvol,sij_opt,vlocal)
    call timab(534,2,tsec)
   end if
   if(abs(dtset%timopt)==3)then
    call timab(584,1,tsec)
   endif
   if (gen_eigenpb) then
    call zcopy(vectsize*blocksize,swavef,1,blockvectorbr,1)
   else
    call zcopy(vectsize*blocksize,blockvectorr,1,blockvectorbr,1)
    call zcopy(vectsize*blocksize,gvnlc,1,blockvectorvr,1)
   endif
   call zcopy(vectsize*blocksize,gwavef,1,blockvectorar,1)
   if(abs(dtset%timopt)==3)then
    call timab(584,2,tsec)
   endif
   deallocate(cwavef,gwavef,gvnlc)
   deallocate(swavef)
   if(abs(dtset%timopt)==3)then
    call timab(535,1,tsec)
   endif
   call zorthonormalize(blockvectorr,blockvectorbr,blocksize,mpi_enreg,gramrbr,vectsize)
   call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramrbr,blocksize,&
&             blockvectorbr,vectsize)
   call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramrbr,blocksize,&
&             blockvectorar,vectsize)
   if (gs_hamk%usepaw==0) then
    call ztrsm('r','u','n','n',vectsize,blocksize,cone,gramrbr,blocksize,&
&              blockvectorvr,vectsize)
   end if
   if(abs(dtset%timopt)==3)then
    call timab(535,2,tsec)
   endif

   if(iterationnumber>1) then
!   call zorthonormalize(blockvectorp,blockvectorbp,blockvectorap)
    if(abs(dtset%timopt)==3)then
     call timab(535,1,tsec)
    endif
    call zorthonormalize(blockvectorp,blockvectorbp,blocksize,mpi_enreg,grampbp,vectsize)
!   blockvectorap=matmul(blockvectorap,grampbp)
    call ztrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,&
&              blockvectorbp,vectsize)
    call ztrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,&
&              blockvectorap,vectsize)
    if (gs_hamk%usepaw==0) then
     call ztrsm('r','u','n','n',vectsize,blocksize,cone,grampbp,blocksize,&
&               blockvectorvp,vectsize)
    end if 
    if(abs(dtset%timopt)==3)then
     call timab(535,2,tsec)
    endif
   end if

   activersize=blocksize
   if (iterationnumber==1) then
    activepsize=0
    restart=1
   else
    activepsize=blocksize
    restart=0
   end if

!  gramxar=matmul(transpose(blockvectorax),blockvectorr)
!  gramrar=matmul(transpose(blockvectorar),blockvectorr)
!  gramxax=matmul(transpose(blockvectorax),blockvectorx)
   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&             vectsize,blockvectorr,vectsize,czero,gramxar,blocksize)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorar,&
&             vectsize,blockvectorr,vectsize,czero,gramrar,blocksize)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&             vectsize,blockvectorx,vectsize,czero,gramxax,blocksize)
   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
   transf(:,:,1)=gramxar(:,:)
   transf(:,:,2)=gramrar(:,:)
   transf(:,:,3)=gramxax(:,:)
   if(abs(dtset%timopt)==3)then
    call timab(533,1,tsec)
   endif
   call xsum_mpi(transf,spaceComm,ierr)
   if(abs(dtset%timopt)==3)then
    call timab(533,2,tsec)
   endif
   gramxar(:,:)=transf(:,:,1)
   gramrar(:,:)=transf(:,:,2)
   gramxax(:,:)=transf(:,:,3)
   mpi_enreg%paral_level= old_paral_level

!  gramxbx=matmul(transpose(blockvectorbx),blockvectorx)
!  gramrbr=matmul(transpose(blockvectorbr),blockvectorr)
!  gramxbr=matmul(transpose(blockvectorbx),blockvectorr)
   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&             vectsize,blockvectorx,vectsize,czero,gramxbx,blocksize)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbr,&
&             vectsize,blockvectorr,vectsize,czero,gramrbr,blocksize)
   call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&             vectsize,blockvectorr,vectsize,czero,gramxbr,blocksize)
   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
   transf(:,:,1)=gramxbx(:,:)
   transf(:,:,2)=gramrbr(:,:)
   transf(:,:,3)=gramxbr(:,:)
   if(abs(dtset%timopt)==3)then
    call timab(533,1,tsec)
   endif
   call xsum_mpi(transf,spaceComm,ierr)
   if(abs(dtset%timopt)==3)then
    call timab(533,2,tsec)
   endif
   gramxbx(:,:)=transf(:,:,1)
   gramrbr(:,:)=transf(:,:,2)
   gramxbr(:,:)=transf(:,:,3)
   mpi_enreg%paral_level= old_paral_level

!###########################################################################
!################ PERFORM LOOP ON COND #####################################
!###########################################################################
   i1=0;i2=blocksize;i3=2*blocksize;i4=3*blocksize
   cond: do cond_try=1,1 !2 when restart, but not implemented
    if (restart==0) then
!    gramxap=matmul(transpose(blockvectorax),blockvectorp)
!    gramrap=matmul(transpose(blockvectorar),blockvectorp)
!    grampap=matmul(transpose(blockvectorap),blockvectorp)
     if(abs(dtset%timopt)==3)then
      call timab(532,1,tsec)
     endif

     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorax,&
&               vectsize,blockvectorp,vectsize,czero,gramxap,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorar,&
&               vectsize,blockvectorp,vectsize,czero,gramrap,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorap,&
&               vectsize,blockvectorp,vectsize,czero,grampap,blocksize)
     if(abs(dtset%timopt)==3)then
      call timab(532,2,tsec)
     endif
     old_paral_level= mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm)
     if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
     transf(:,:,1)=gramxap(:,:)
     transf(:,:,2)=gramrap(:,:)
     transf(:,:,3)=grampap(:,:)
     if(abs(dtset%timopt)==3)then
      call timab(533,1,tsec)
     endif
     call xsum_mpi(transf,spaceComm,ierr)
     if(abs(dtset%timopt)==3)then
      call timab(533,2,tsec)
     endif
     gramxap(:,:)=transf(:,:,1)
     gramrap(:,:)=transf(:,:,2)
     grampap(:,:)=transf(:,:,3)
     mpi_enreg%paral_level= old_paral_level
     bigorder=i4
     allocate(grama(i4,i4),gramb(i4,i4),eigen(i4),coordx(i4,blocksize))
     grama(i1+1:i2,i1+1:i2)=gramxax
     grama(i1+1:i2,i2+1:i3)=gramxar
     grama(i1+1:i2,i3+1:i4)=gramxap
!     grama(i2+1:i3,i1+1:i2)=transpos(gramxar)
     grama(i2+1:i3,i2+1:i3)=gramrar
     grama(i2+1:i3,i3+1:i4)=gramrap
!     grama(i3+1:i4,i1+1:i2)=transpos(gramxap)
!     grama(i3+1:i4,i2+1:i3)=transpos(gramrap)
     grama(i3+1:i4,i3+1:i4)=grampap

!    gramxbp=matmul(transpose(blockvectorbx),blockvectorp)
!    gramrbp=matmul(transpose(blockvectorbr),blockvectorp)
!    grampbp=matmul(transpose(blockvectorbp),blockvectorp)
     if(abs(dtset%timopt)==3)then
      call timab(532,1,tsec)
     endif
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbx,&
&               vectsize,blockvectorp,vectsize,czero,gramxbp,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbr,&
&               vectsize,blockvectorp,vectsize,czero,gramrbp,blocksize)
     call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorbp,&
&               vectsize,blockvectorp,vectsize,czero,grampbp,blocksize)
     if(abs(dtset%timopt)==3)then
      call timab(532,2,tsec)
     endif
     old_paral_level= mpi_enreg%paral_level
     mpi_enreg%paral_level=3
     call xcomm_init(mpi_enreg,spaceComm)
     if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
     transf(:,:,1)=gramxbp(:,:)
     transf(:,:,2)=gramrbp(:,:)
     transf(:,:,3)=grampbp(:,:)
     if(abs(dtset%timopt)==3)then
      call timab(533,1,tsec)
     endif
     call xsum_mpi(transf,spaceComm,ierr)
     if(abs(dtset%timopt)==3)then
      call timab(533,2,tsec)
     endif
     gramxbp(:,:)=transf(:,:,1)
     gramrbp(:,:)=transf(:,:,2)
     grampbp(:,:)=transf(:,:,3)
     mpi_enreg%paral_level= old_paral_level
     gramb(i1+1:i2,i1+1:i2)=gramxbx
     gramb(i1+1:i2,i2+1:i3)=gramxbr
     gramb(i1+1:i2,i3+1:i4)=gramxbp
!     gramb(i2+1:i3,i1+1:i2)=transpos(gramxbr)
     gramb(i2+1:i3,i2+1:i3)=gramrbr
     gramb(i2+1:i3,i3+1:i4)=gramrbp
!     gramb(i3+1:i4,i1+1:i2)=transpos(gramxbp)
!     gramb(i3+1:i4,i2+1:i3)=transpos(gramrbp)
     gramb(i3+1:i4,i3+1:i4)=grampbp

    else
     bigorder=i3
     allocate(grama(i3,i3),gramb(i3,i3),eigen(i3),coordx(i3,blocksize))
     grama(i1+1:i2,i1+1:i2)=gramxax
     grama(i1+1:i2,i2+1:i3)=gramxar
!     grama(i2+1:i3,i1+1:i2)=transpos(gramxar)
     grama(i2+1:i3,i2+1:i3)=gramrar
     gramb(i1+1:i2,i1+1:i2)=gramxbx
     gramb(i1+1:i2,i2+1:i3)=gramxbr
!     gramb(i2+1:i3,i1+1:i2)=transpos(gramxbr)
     gramb(i2+1:i3,i2+1:i3)=gramrbr

    end if
   end do cond
!###########################################################################
!################ END LOOP ON COND #########################################
!###########################################################################

   lwork=3*bigorder-2
   allocate(work(lwork),rwork(lwork))
   if(abs(dtset%timopt)==3)then
    call timab(587,1,tsec)
   endif
   call zhegv(1,'v','u',bigorder,grama,bigorder,gramb,bigorder,eigen,work,lwork,rwork,info)
   if(abs(dtset%timopt)==3)then
    call timab(587,2,tsec)
   endif
   deallocate(work,rwork)

   deltae=-one
   do iblocksize=1,blocksize
    deltae=max(deltae,abs(lambda(iblocksize,iblocksize)-eigen(iblocksize)))
    lambda(iblocksize,iblocksize)=eigen(iblocksize)
   end do
!DEBUG
!  write(6,*)'eigen',eigen(1:blocksize)
!ENDDEBUG
   coordx(1:bigorder,1:blocksize)=grama(1:bigorder,1:blocksize)
   deallocate(grama,gramb,eigen)
   if (restart==0 .and. iterationnumber >1) then

!   blockvectorp=matmul(blockvectorr,coordx(i2+1:i3,:))+&
!&               matmul(blockvectorp,coordx(i3+1:i4,:))
    if(abs(dtset%timopt)==3)then
     call timab(532,1,tsec)
    endif
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorr,&
&              vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorp,&
&              vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
    if(abs(dtset%timopt)==3)then
     call timab(532,2,tsec)
    endif
    if(abs(dtset%timopt)==3)then
     call timab(584,1,tsec)
    endif
    call zcopy(vectsize*blocksize,blockvectordumm,1,blockvectorp,1)
    if(abs(dtset%timopt)==3)then
     call timab(584,2,tsec)
    endif

!   blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))+&
!&                matmul(blockvectorap,coordx(i3+1:i4,:))
    if(abs(dtset%timopt)==3)then
     call timab(532,1,tsec)
    endif
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorar,&
&              vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorap,&
&              vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
    if(abs(dtset%timopt)==3)then
     call timab(532,2,tsec)
    endif
    if(abs(dtset%timopt)==3)then
     call timab(584,1,tsec)
    endif
    call zcopy(vectsize*blocksize,blockvectordumm,1,blockvectorap,1)
    if(abs(dtset%timopt)==3)then
     call timab(584,2,tsec)
    endif

!   blockvectorvp=matmul(blockvectorvr,coordx(i2+1:i3,:))+&
!&                matmul(blockvectorvp,coordx(i3+1:i4,:))
    if (gs_hamk%usepaw==0) then
     if(abs(dtset%timopt)==3)then
      call timab(532,1,tsec)
     endif
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvr,&
&               vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvp,&
&               vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
     if(abs(dtset%timopt)==3)then
      call timab(532,2,tsec)
     endif
     if(abs(dtset%timopt)==3)then
      call timab(584,1,tsec)
     endif
     call zcopy(vectsize*blocksize,blockvectordumm,1,blockvectorvp,1)
     if(abs(dtset%timopt)==3)then
      call timab(584,2,tsec)
     endif
    end if 

!   blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))+&
!&                matmul(blockvectorbp,coordx(i3+1:i4,:))
    if(abs(dtset%timopt)==3)then
     call timab(532,1,tsec)
    endif
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbr,&
&              vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectordumm,vectsize)
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbp,&
&              vectsize,coordx(i3+1:i4,:),blocksize,cone,blockvectordumm,vectsize)
    if(abs(dtset%timopt)==3)then
     call timab(532,2,tsec)
    endif
    if(abs(dtset%timopt)==3)then
     call timab(584,1,tsec)
    endif
    call zcopy(vectsize*blocksize,blockvectordumm,1,blockvectorbp,1)
    if(abs(dtset%timopt)==3)then
     call timab(584,2,tsec)
    endif

   else

!   blockvectorp =matmul(blockvectorr,coordx(i2+1:i3,:))
    if(abs(dtset%timopt)==3)then
     call timab(532,1,tsec)
    endif
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorr,&
&              vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorp,vectsize)
!   blockvectorap=matmul(blockvectorar,coordx(i2+1:i3,:))
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorar,&
&              vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorap,vectsize)
!   blockvectorvp=matmul(blockvectorvr,coordx(i2+1:i3,:))
    if (gs_hamk%usepaw==0) then
     call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvr,&
&               vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorvp,vectsize)
    end if
!   blockvectorbp=matmul(blockvectorbr,coordx(i2+1:i3,:))
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbr,&
&              vectsize,coordx(i2+1:i3,:),blocksize,czero,blockvectorbp,vectsize)
    if(abs(dtset%timopt)==3)then
     call timab(532,2,tsec)
    endif
   end if

   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif

!  blockvectorx = matmul(blockvectorx,coordx(i1+1:i2,:))+blockvectorp
   call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorx,&
&             vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
   blockvectorx = blockvectordumm+blockvectorp

!  blockvectorax= matmul(blockvectorax,coordx(i1+1:i2,:))+blockvectorap
   call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorax,&
&             vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
   blockvectorax = blockvectordumm+blockvectorap

!  blockvectorvx= matmul(blockvectorvx,coordx(i1+1:i2,:))+blockvectorvp
   if (gs_hamk%usepaw==0) then
    call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorvx,&
&              vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
    blockvectorvx = blockvectordumm+blockvectorvp
   end if 

!  blockvectorbx= matmul(blockvectorbx,coordx(i1+1:i2,:))+blockvectorbp
   call zgemm('n','n',vectsize,blocksize,blocksize,cone,blockvectorbx,&
&             vectsize,coordx(i1+1:i2,:),blocksize,czero,blockvectordumm,vectsize)
   blockvectorbx = blockvectordumm+blockvectorbp

   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif

   deallocate(coordx)

!  Check convergence on energy and eventually exit
   if (iterationnumber==1) then
    deold=deltae
   else if (iterationnumber>1) then
    if ((abs(deltae)<0.005*abs(deold)).and.(iterationnumber/=maxiterations))then
     write(message, '(2(a,i4),1x,a,1p,e12.4,a,e12.4,a)' ) &
&      ' lobpcgccwf: block',iblock,', line',iterationnumber,&
&      ', deltae=',deltae,' < 0.005*',deold,' =>skip lines !'
     call wrtout(06,message,'PERS')
     exit
    else if (abs(deltae)>0.005*abs(deold)) then
     write(message, '(2(a,i4),1x,a,1p,e12.4,a,e12.4,a)' ) &
&      ' lobpcgccwf: block',iblock,', line',iterationnumber,&
&      ', deltae=',deltae,' > 0.005*',deold,' =>keep on working !'
     call wrtout(06,message,'PERS')
    end if
   end if

  end do iter
!###########################################################################
!################## END LOOP ON NLINE ######################################
!###########################################################################
  if (iterationnumber==maxiterations+1) then
   call zprecon3(blockvectorbx,lambda,blocksize,&
&                istwf_k,kinpw,mpi_enreg,npw_k,nspinor,&
&                optekin,blockvectorax,blockvectorr,vectsize)
   residualnorms=sum(abs(blockvectorr)**2,dim=1)
   old_paral_level= mpi_enreg%paral_level
   mpi_enreg%paral_level=3
   call xcomm_init(mpi_enreg,spaceComm)
   if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
   if(abs(dtset%timopt)==3)then
    call timab(533,1,tsec)
   endif
   call xsum_mpi(residualnorms,spaceComm,ierr)
   if(abs(dtset%timopt)==3)then
    call timab(533,2,tsec)
   endif
   mpi_enreg%paral_level= old_paral_level
   resid_k(bblocksize+1:bblocksize+blocksize)=residualnorms(1:blocksize)
  endif

  if(abs(dtset%timopt)==3)then
   call timab(584,1,tsec)
  endif
  call zcopy(vectsize*blocksize,blockvectorx,1,cg(:,vectsize*((iblock-1)*blocksize)+icg+1:vectsize*(iblock*blocksize-1)+icg),1)
  if(gen_eigenpb) then
   call zcopy(vectsize*blocksize,blockvectorbx,1,gsc(:,vectsize*((iblock-1)*blocksize)+igsc+1:vectsize*(iblock*blocksize-1)+igsc),1)
  end if
  if(abs(dtset%timopt)==3)then
   call timab(584,2,tsec)
  endif

!The Vnl part of the Hamiltonian is no more stored in the packed form such as it was the case for subvnl(:).
!Now, the full matrix is stored in totvnl(:,:). This trick permits:
!1) to avoid the reconstruction of the total matrix in vtowfk.F90 (double loop over bands)
!2) to use two optimized matrix-matrix blas routine for general (in lobpcgccwf.F90) or hermitian (in vtowfk.F90)
!   operators, zgemm.f and zhemm.f respectively, rather than a triple loop in both cases.
  iwavef=iblock*blocksize
  isubh=1+2*bblocksize*(bblocksize+1)/2

  allocate(blockvectorz(vectsize,iwavef))
  if(abs(dtset%timopt)==3)then
   call timab(584,1,tsec)
  endif
  call zcopy(bblocksize*vectsize,blockvectory(:,1:bblocksize),1,blockvectorz(:,1:bblocksize),1)
  call zcopy( blocksize*vectsize,blockvectorx(:,1:blocksize) ,1,blockvectorz(:,bblocksize+1:iwavef),1)
  if(abs(dtset%timopt)==3)then
   call timab(584,2,tsec)
  endif

  allocate(tsubham(2*iwavef,blocksize))
  tsubham(:,:)=zero
  if(abs(dtset%timopt)==3)then
   call timab(532,1,tsec)
  endif
  call zgemm('c','n',iwavef,blocksize,vectsize,cone,blockvectorz,vectsize,&
&            blockvectorax,vectsize,czero,tsubham,iwavef)
  if(abs(dtset%timopt)==3)then
   call timab(532,2,tsec)
  endif
  if (gs_hamk%usepaw==0) then
   if(abs(dtset%timopt)==3)then
    call timab(532,1,tsec)
   endif
   call zgemm('c','n',blocksize,iwavef,vectsize,cone,blockvectorvx,vectsize,&
&              blockvectorz,vectsize,czero,totvnl(2*bblocksize+1:2*iwavef,1:iwavef),blocksize)
   if(abs(dtset%timopt)==3)then
    call timab(532,2,tsec)
   endif
  endif
  do iblocksize=1,blocksize
   do ii=1,bblocksize+iblocksize
    subham(isubh)  = tsubham(2*ii-1,iblocksize)
    subham(isubh+1)= tsubham(2*ii  ,iblocksize)
    isubh=isubh+2
   end do
  end do
  deallocate(tsubham)
  deallocate(blockvectorz)
! comm for subham and subvnl is made in vtowfk

!###########################################################################
!################## END OF LOBPCG IN MOST CASES ############################
!###########################################################################
! call operators(blockvectorx,blockvectorbx,subovl)!fill also  subovl
  if((gen_eigenpb).and.(use_subovl==1)) then
   allocate(cwavef(2,npw_k*nspinor))
   allocate(gwavef(2,npw_k*nspinor))
   isubo=1+2*(iblock-1)*maxblocksize*((iblock-1)*maxblocksize+1)/2
   do iblocksize=1,blocksize
    cwavef(1,1:npw_k*nspinor)=real (blockvectorx(1:npw_k*nspinor,iblocksize))
    cwavef(2,1:npw_k*nspinor)=aimag(blockvectorx(1:npw_k*nspinor,iblocksize))
!   Call to nonlop: compute <g|S|c>

    choice=1 ; signs=2 ; idir=0 ; tim_nonlop=311 ; cpopt=-1 ; paw_opt=3 ; nnlout=1; nkpg=0

!   MOST UGLY PATCH, TO BE DELETED
    spaceComm_keep=mpi_enreg%comm_fft;mpi_enreg%comm_fft=mpi_enreg%commcart

    allocate(enlout(nnlout*blocksize))

    call nonlop(gs_hamk%atindx1,choice,cpopt,cprj_dum,gs_hamk%dimekb1,0,dimffnl,dimffnl,dummy3,&
&               enlout,ffnl,ffnl,gs_hamk%gmet,gs_hamk%gprimd,idir,gs_hamk%indlmn,&
&               istwf_k,kg_k,kg_k,kpg_dum,kpg_dum,gs_hamk%kpoint,gs_hamk%kpoint,dum,lmnmax,matblk,&
&               mgfft,mpi_enreg,mpsang,mpssoang,natom,gs_hamk%nattyp,gs_hamk%ngfft,nkpg,nkpg,&
&               gs_hamk%nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,gs_hamk%phkxred,&
&               gs_hamk%phkxred,gs_hamk%ph1d,ph3d,ph3d,gs_hamk%pspso,signs,gs_hamk%sij,&
&               gwavef,tim_nonlop,gs_hamk%ucvol,gs_hamk%useylm,cwavef,cwavef)

    deallocate(enlout)

    mpi_enreg%comm_fft=spaceComm_keep

    if(abs(dtset%timopt)==3)then
     call timab(584,1,tsec)
    endif
    call zcopy(npw_k*nspinor,gwavef(:,:),1,blockvectorbx(:,iblocksize),1)
    if(abs(dtset%timopt)==3)then
     call timab(584,2,tsec)
    endif
    do ii=1,(iblock-1)*maxblocksize+iblocksize
     iwavef=(ii-1)*npw_k*nspinor+icg
     cscre=zero;cscim=zero
     do ipw=1,npw_k*nspinor
      cgreipw=cg(1,ipw+iwavef);cgimipw=cg(2,ipw+iwavef)
      cscre=cscre+cgreipw*gwavef(1,ipw)+cgimipw*gwavef(2,ipw)
      cscim=cscim+cgreipw*gwavef(2,ipw)-cgimipw*gwavef(1,ipw)
     end do
!    Store real and imag parts in hermitian storage mode:
     subovl(isubo)=cscre ; subovl(isubo+1)=cscim
     isubo=isubo+2
    end do
   end do
   deallocate(cwavef,gwavef)
  end if

  deallocate(blockvectory,blockvectorby,gramyx)
  deallocate(blockvectorx,blockvectorax,blockvectorbx)
  deallocate(blockvectorr,blockvectorar,blockvectorbr)
  deallocate(blockvectorp,blockvectorap,blockvectorbp)
  if (gs_hamk%usepaw==0) then
   deallocate(blockvectorvx,blockvectorvp,blockvectorvr)
  end if 
  deallocate(blockvectordumm)
  deallocate(gramxax,gramxar,gramxap,gramrar,gramrap,grampap,gramxbx,gramxbr,&
&  gramxbp,gramrbr,gramrbp,grampbp,transf)
  deallocate(lambda)
  deallocate(residualnorms)

!End big loop over bands inside blocks
 end do

 call timab(530,2,tsec)

end subroutine lobpcgccwf
!!***
