!{\src2tex{textfont=tt}}
!!****f* abinit/prep_nonlop
!! NAME
!! prep_nonlop
!!
!! FUNCTION
!! this routine prepares the data to the call of nonlop.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FBottin,MT,GZ)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  choice: chooses possible output:
!!    choice=1 => a non-local energy contribution
!!          =2 => a gradient with respect to atomic position(s)
!!          =3 => a gradient with respect to strain(s)
!!          =23=> a gradient with respect to atm. pos. and strain(s)
!!          =4 => a 2nd derivative with respect to atomic pos.
!!          =24=> a gradient and 2nd derivative with respect to atomic pos.
!!          =5 => a gradient with respect to k wavevector
!!          =6 => 2nd derivatives with respect to strain and atm. pos.
!!  blocksize= size of block for FFT
!!  cpopt=flag defining the status of cprj_block=<Proj_i|Cnk> scalars (see below, side effects)
!!  cwavef(2,npw*nspinor*blocksize)=planewave coefficients of wavefunction.
!!  dimenl1,dimenl2=dimensions of enl (see enl)
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  enl(dimenl1,dimenl2,nspinor**2)=
!!  ->Norm conserving : ==== when paw_opt=0 ====
!!                      (Real) Kleinman-Bylander energies (hartree)
!!                      dimenl1=lmnmax  -  dimenl2=ntypat
!!  ->PAW :             ==== when paw_opt=1 or 4 ====
!!                      (Real, symmetric) Dij coefs to connect projectors
!!                      dimenl1=cplex_enl*lmnmax*(lmnmax+1)/2  -  dimenl2=natom
!!                      These are complex numbers if cplex_enl=2
!!                        enl(:,:,1) contains Dij^up-up
!!                        enl(:,:,2) contains Dij^dn-dn
!!                        enl(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!                        enl(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gvnlc=matrix elements <G|Vnonlocal|C>
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kpg(npw,nkpg)= if nkpg==3 (k+G) components
!!                 if nkpg==9 [(k+G)_a].[(k+G)_b] quantities
!!  kpt(3)=k point in terms of recip. translations
!!  icall = order of call of this routine in lobpcgccwf
!!  idir=direction of the - atom to be moved in the case (choice=2,signs=2),
!!                        - k point direction in the case (choice=5,signs=2)
!!                        - strain component (1:6) in the case (choice=2,signs=2) or (choice=6,signs=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln  (if useylm=0)
!!                                                  or i=lmn (if useylm=1)
!!  istwf_k=option parameter that describes the storage of wfs
!!  lambda=factor to be used when computing (Vln-lambda.S) - only for paw_opt=2
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1d ffts
!!  mpi_enreg=informations about mpi parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type
!!  nband_k=number of bands at this k point for that spin polarization
!!  nbdblock=
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpg=second dimension of array kpg
!!  nloalg(5)=governs the choice of the algorithm for nonlocal operator
!!  nnlout=dimension of enlout (when signs=1):
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in cell
!!  paw_opt= define the nonlocal operator concerned with
!!  phkxred(2,natom)=phase factors exp(2 pi kpt.xred)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1D structure factors phase information
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  pspso(ntypat)=spin-orbit characteristic for each atom type
!!  signs= if 1, get contracted elements (energy, forces, stress, ...)
!!         if 2, applies the non-local operator to a function in reciprocal space
!!  sij(dimenl1,ntypat*(paw_opt/3))=overlap matrix components (only if paw_opt=2, 3 or 4)
!!  tim_nonlop=timing code of the calling routine (can be set to 0 if not attributed)
!!  ucvol=unit cell volume (bohr^3)
!!  useylm=governs the way the nonlocal operator is to be applied:
!!         1=using Ylm, 0=using Legendre polynomials
!!
!! OUTPUT
!!  ==== if (signs==1) ====
!!  enlout(nnlout)=
!!    if paw_opt==0, 1 or 2: contribution of this state to the nl part of various properties
!!    if paw_opt==3:        contribution of this state to <c|S|c>  (where S=overlap when PAW)
!!  ==== if (signs==2) ====
!!    if paw_opt==0, 1, 2 or 4:
!!       gvnlc(2,nspinor*npw)=result of the aplication of the nl operator
!!                        or one of its derivative to the input vect.
!!    if paw_opt==3 or 4:
!!       gsc(2,nspinor*npw*(paw_opt/3))=result of the aplication of (I+S)
!!                        to the input vect. (where S=overlap when PAW)
!!
!! SIDE EFFECTS
!!  ==== ONLY IF useylm=1
!!  cprj_block(natom) <type(cprj_type)>=projected input wave function |c> on non-local projector
!!                                  =<p_lmn|c> and derivatives
!!                                  Treatment depends on cpopt parameter:
!!                     if cpopt=-1, <p_lmn|in> (and derivatives)
!!                                  have to be computed (and not saved)
!!                     if cpopt= 0, <p_lmn|in> have to be computed and saved
!!                                  derivatives are eventually computed but not saved
!!                     if cpopt= 1, <p_lmn|in> and first derivatives have to be computed and saved
!!                                  other derivatives are eventually computed but not saved
!!                     if cpopt= 2  <p_lmn|in> are already in memory;
!!                                  only derivatives are computed here and not saved
!! (if useylm=0, should have cpopt=-1)
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      cprj_alloc,cprj_copy,cprj_free,cprj_nullify,leave_new,mpi_barrier
!!      nonlop,prep_index_wavef_bandpp,timab,wrtout,xallgather_mpi
!!      xalltoallv_mpi,xcomm_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prep_nonlop(atindx1,choice,cpopt,cprj_block,dimenl1,dimenl2,dimffnl,enl,enlout_block,&
&                       ffnl_gather,gmet,gprimd,iblock,icall,idir,indlmn,istwf_k,kg_k_gather,&
&                       kpg,kpt,lambdablock,lmnmax,matblk,&
&                       blocksize,mgfft,mpi_enreg,mpsang,mpssoang,&
&                       natom,nattyp,nbdblock,nband_k,dimtabs,ngfft,nkpg,nloalg,nnlout,npw_k,&
&                       nspinor,ntypat,paral_kgb,paw_opt,phkxred,ph1d,ph3d_gather,prtvol,pspso,signs,sij,gsc,&
&                       tim_nonlop,ucvol,useylm,cwavef,gvnlc)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_13nonlocal
 use interfaces_14wfs
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 integer,intent(in) :: blocksize,choice,cpopt,dimenl1,dimenl2,dimffnl,dimtabs,iblock,icall,idir,istwf_k
 integer,intent(in) :: lmnmax,matblk,mgfft,mpsang,mpssoang,signs
 integer,intent(in) :: natom,nband_k,nbdblock,nkpg,nnlout,npw_k,nspinor,ntypat,paral_kgb,paw_opt,prtvol,useylm
 real(dp),intent(in) :: ucvol
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: atindx1(natom),indlmn(6,lmnmax,ntypat),kg_k_gather(3,dimtabs)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),nloalg(5),pspso(ntypat)
 real(dp),intent(in) :: enl(dimenl1,dimenl2,nspinor**2)
 real(dp),intent(in) :: ffnl_gather(dimtabs,dimffnl,lmnmax,ntypat)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),lambdablock(blocksize)
 real(dp),intent(in) :: kpt(3)
 real(dp),intent(in) :: kpg(npw_k,nkpg)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),phkxred(2,natom)
 real(dp),intent(in) :: sij(dimenl1,ntypat*((paw_opt+1)/3))
 real(dp),intent(out) :: enlout_block(nnlout*blocksize),gvnlc(2,npw_k*nspinor*blocksize)
 real(dp),intent(out) :: gsc(2,npw_k*nspinor*blocksize*(paw_opt/3))
 real(dp),intent(inout) :: cwavef(2,npw_k*nspinor*blocksize),ph3d_gather(2,dimtabs,matblk)
 type(cprj_type) :: cprj_block(natom,nspinor*blocksize*((cpopt+3)/3))

!Local variables-------------------------------
 character(len=500) :: message
 integer :: spaceComm=0
 integer :: bufdim,iat,ib,ier,ii,ilmn,ispinor,mu,ncpgr,npw_loc
 integer :: old_me_g0,old_paral_level,tim_nonlop
 real(dp) :: lambda
 real(dp) :: tsec(2)
 integer,allocatable :: dimlmn(:)
 real(dp), allocatable :: enlout(:)
 type(cprj_type),allocatable :: cwaveprj(:,:)

!local variables for mpialltoallv
 integer :: iproc,ndatarecv
 integer,  allocatable :: recvcounts(:),sendcounts(:),sdispls(:),rdispls(:)
 integer,  allocatable :: sendcountsloc(:),sdisplsloc(:),recvcountsloc(:),rdisplsloc(:)
 real(dp), allocatable :: buffer1(:),buffer2(:)
 real(dp), allocatable :: cwavef_alltoall(:,:),gvnlc_alltoall(:,:),gsc_alltoall(:,:)

!local variables for bandpp
 integer :: nproc_band,bandpp,idebe,idebc,ifine,ifinc
 integer :: ibandpp,kbandpp,iindex
 integer,  pointer :: index_wavef_band(:)
 real(dp), allocatable :: gsc_alltoall_loc(:,:)
 real(dp), allocatable :: cwavef_alltoall_loc(:,:)
 real(dp), allocatable :: gvnlc_alltoall_loc(:,:)
 type(cprj_type),allocatable :: cwaveprj_loc(:,:)

! *************************************************************************

 nproc_band = mpi_enreg%nproc_band
 bandpp     = mpi_enreg%bandpp

 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm)
 if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_band

 allocate(sendcountsloc(nproc_band))
 allocate(sdisplsloc   (nproc_band))
 allocate(sdispls      (nproc_band))
 allocate(sendcounts   (nproc_band))
 allocate(recvcountsloc(nproc_band))
 allocate(rdisplsloc   (nproc_band))
 allocate(rdispls      (nproc_band))
 allocate(recvcounts   (nproc_band))
 call timab(582,1,tsec)
 npw_loc=npw_k
 call xallgather_mpi(npw_loc,recvcounts,spaceComm,ier)
 call timab(582,2,tsec)
 rdispls(1)=0
 do iproc=2,nproc_band
  rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
 end do
 ndatarecv=rdispls(nproc_band)+recvcounts(nproc_band)

!The transposition are no more performed within prep_fourwf but above.
!We have to check that the dimension mpw is identical the the one computed above. TO BE CLEANED
 if (ndatarecv /= dimtabs) then
  write(message, '(a,a,a,a,i8,a,i8)' ) ch10,&
&      ' prep_fourwf: BUG -',ch10,&
&      '  ndatarecv',ndatarecv,' /= dimtabs ', dimtabs
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 endif

 sendcounts(:)=npw_k*bandpp
 do iproc=1,nproc_band
  sdispls(iproc)=(iproc-1)*npw_k*bandpp
 end do

 allocate(cwavef_alltoall(2,ndatarecv*nspinor*bandpp))
 allocate(gsc_alltoall(2,ndatarecv*nspinor*(paw_opt/3)*bandpp))
 allocate(gvnlc_alltoall(2,ndatarecv*nspinor*bandpp))

 allocate(enlout(nnlout))
 if (cpopt>=0) then
  ncpgr=cprj_block(1,1)%ncpgr
  allocate(dimlmn(natom))
  do iat=1,natom;dimlmn(iat)=cprj_block(iat,1)%nlmn;end do
  allocate(cwaveprj(natom,nspinor*bandpp))
  call cprj_alloc(cwaveprj,ncpgr,dimlmn)
 end if
 recvcountsloc(:)=recvcounts(:)*2*nspinor*bandpp
 rdisplsloc(:)=rdispls(:)*2*nspinor*bandpp
 sendcountsloc(:)=sendcounts(:)*2*nspinor
 sdisplsloc(:)=sdispls(:)*2*nspinor
 call timab(581,1,tsec)
 call xalltoallv_mpi(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,&
&                    recvcountsloc,rdisplsloc,spaceComm,ier)
 call timab(581,2,tsec)

 if(istwf_k==2) then
  old_me_g0=mpi_enreg%me_g0
  if (mpi_enreg%me_fft==0) then
   mpi_enreg%me_g0=1
  else
   mpi_enreg%me_g0=0
  endif
 endif

 if (paw_opt==2) lambda=lambdablock(mpi_enreg%me_band+1)

!=====================================================================
 if  (bandpp==1) then
  call nonlop(atindx1,choice,cpopt,cwaveprj,dimenl1,dimenl2,dimffnl,dimffnl,&
&             enl,enlout,ffnl_gather,ffnl_gather,gmet,gprimd,idir,indlmn,&
&             istwf_k,kg_k_gather,kg_k_gather,kpg,kpg,kpt,kpt,lambda,lmnmax,matblk,mgfft,&
&             mpi_enreg,mpsang,mpssoang,natom,nattyp,ngfft,nkpg,nkpg,nloalg,&
&             nnlout,ndatarecv,ndatarecv,nspinor,ntypat,0,paw_opt,phkxred,&
&             phkxred,ph1d,ph3d_gather,ph3d_gather,pspso,signs,sij,gsc_alltoall,&
&             tim_nonlop,ucvol,useylm,cwavef_alltoall,gvnlc_alltoall)
 else

  ! -------------------------------------------------------
  ! Allocation
  allocate(cwavef_alltoall_loc(2,ndatarecv*nspinor))
  allocate(gsc_alltoall_loc   (2,ndatarecv*nspinor*(paw_opt/3)))
  allocate(gvnlc_alltoall_loc (2,ndatarecv*nspinor))
  if (cpopt>=0) then
   allocate(cwaveprj_loc(natom,nspinor))
   call cprj_alloc(cwaveprj_loc,ncpgr,dimlmn)
  end if
  ! -------------------------------------------------------------
  ! Computation of the index used to sort the waves functions below bandpp
  call prep_index_wavef_bandpp(nproc_band,bandpp,&
&                              1,ndatarecv,&
&                              recvcounts,rdispls,&
&                              index_wavef_band)
  ! -------------------------------------------------------
  ! Sorting of the waves functions below bandpp
  cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)
  if (paw_opt>=3) gsc_alltoall(:,:)=gsc_alltoall(:,index_wavef_band)
  gvnlc_alltoall(:,:)  = gvnlc_alltoall(:,index_wavef_band)
  ! -------------------------------------------------------
  ! Computation for each bandpp
  do ibandpp=1,bandpp
   cwavef_alltoall_loc(:,:)=cwavef_alltoall(:,(ibandpp-1)*ndatarecv*nspinor+1:(ibandpp)*ndatarecv*nspinor)
   call nonlop(atindx1,choice,cpopt,cwaveprj,dimenl1,dimenl2,dimffnl,dimffnl,&
&       enl,enlout,ffnl_gather,ffnl_gather,gmet,gprimd,idir,indlmn,&
&       istwf_k,kg_k_gather,kg_k_gather,kpg,kpg,kpt,kpt,lambda,lmnmax,matblk,mgfft,&
&       mpi_enreg,mpsang,mpssoang,natom,nattyp,ngfft,nkpg,nkpg,nloalg,&
&       nnlout,ndatarecv,ndatarecv,nspinor,ntypat,0,paw_opt,phkxred,&
&       phkxred,ph1d,ph3d_gather,ph3d_gather,pspso,signs,sij,gsc_alltoall_loc,&
&       tim_nonlop,ucvol,useylm,cwavef_alltoall_loc,gvnlc_alltoall_loc)
   cwavef_alltoall(:,(ibandpp-1)*ndatarecv*nspinor+1:(ibandpp)*ndatarecv*nspinor)=cwavef_alltoall_loc(:,:)
   if (paw_opt>=3) gsc_alltoall(:,(ibandpp-1)*ndatarecv*nspinor+1:(ibandpp)*ndatarecv*nspinor)=gsc_alltoall_loc(:,:)
   gvnlc_alltoall(:,(ibandpp-1)*ndatarecv*nspinor+1:(ibandpp)*ndatarecv*nspinor)=gvnlc_alltoall_loc(:,:)
   if (cpopt>=0) then
    call cprj_copy(cwaveprj_loc,cwaveprj(:,(ibandpp-1)*nspinor+1:ibandpp*nspinor))
   end if
  enddo
  ! -----------------------------------------------------
  ! Sorting of waves functions below the processors
  cwavef_alltoall(:,index_wavef_band)=cwavef_alltoall(:,:)
  if (paw_opt>=3) gsc_alltoall(:,index_wavef_band)=gsc_alltoall(:,:)
  gvnlc_alltoall(:,index_wavef_band)=gvnlc_alltoall(:,:)
  ! -------------------------------------------------------
  ! Deallocation
  deallocate(index_wavef_band)
  deallocate(cwavef_alltoall_loc,gsc_alltoall_loc,gvnlc_alltoall_loc)
  if (cpopt>=0) then
   call cprj_free(cwaveprj_loc)
   deallocate(cwaveprj_loc)
  end if
 endif
! =====================================================================

!Transpose the gsc_alltoall or gvlnc_alltoall tabs
!according to the paw_opt and signs values
 call timab(581,1,tsec)
 if (signs==2 .and. (paw_opt==0 .or. paw_opt==1 .or. paw_opt==4)) then
  call xalltoallv_mpi(gvnlc_alltoall,recvcountsloc,rdisplsloc,gvnlc,&
&                     sendcountsloc,sdisplsloc,spaceComm,ier)
 end if
 if (signs==2 .and. (paw_opt==3 .or. paw_opt==4)) then
  call xalltoallv_mpi(gsc_alltoall,recvcountsloc,rdisplsloc,gsc,&
&                     sendcountsloc,sdisplsloc,spaceComm,ier)
 end if
 call timab(581,2,tsec)
 if (istwf_k==2) mpi_enreg%me_g0=old_me_g0

#if defined MPI
 if(paral_kgb == 1) then
  call xallgather_mpi(enlout,nnlout,enlout_block,spaceComm,ier)
  if (cpopt==0.or.cpopt==1) then
   bufdim=2*nspinor*bandpp*sum(dimlmn);allocate(buffer1(bufdim),buffer2(bufdim*blocksize))
   ii=1
   do ispinor=1,nspinor
    do iat=1,natom
     do ilmn=1,dimlmn(iat)
      buffer1(ii:ii+1)=cwaveprj(iat,ispinor)%cp(1:2,ilmn)
      ii=ii+2
     end do
    end do
   end do
   call xallgather_mpi(buffer1,bufdim,buffer2,spaceComm,ier)
   !TD 04/03/2008 : MPI_BARRIER requested to avoid bug with FFT parallelisation
   call MPI_BARRIER(spaceComm,ier)
   !END TD 04/03/2008
   ii=1
   do ib=1,blocksize*nspinor
    do iat=1,natom
     do ilmn=1,dimlmn(iat)
      cprj_block(iat,ib)%cp(1:2,ilmn)=buffer2(ii:ii+1)
      ii=ii+2
     end do
    end do
   end do
   deallocate(buffer1,buffer2)
  end if
  if (cpopt==1.and.ncpgr>0) then
   bufdim=2*nspinor*bandpp*sum(dimlmn)*ncpgr
   allocate(buffer1(bufdim),buffer2(bufdim*blocksize))
   ii=1
   do ispinor=1,nspinor
    do iat=1,natom
     do ilmn=1,dimlmn(iat)
      do mu=1,ncpgr
       buffer1(ii:ii+1)=cwaveprj(iat,ispinor)%dcp(1:2,mu,ilmn)
       ii=ii+2
      end do
     end do
    end do
   end do
   call xallgather_mpi(buffer1,bufdim,buffer2,spaceComm,ier)
   ii=1
   do ib=1,blocksize*nspinor
    do iat=1,natom
     do ilmn=1,dimlmn(iat)
      do mu=1,ncpgr
       cprj_block(iat,ib)%dcp(1:2,mu,ilmn)=buffer2(ii:ii+1)
       ii=ii+2
      end do
     end do
    end do
   end do
   deallocate(buffer1,buffer2)
  end if
 end if
#endif

 if(paral_kgb == 0)then
  enlout_block(:) = zero
  if (cpopt>=0) then
   call cprj_nullify(cprj_block)
  end if
 endif

 deallocate(enlout)
 if (cpopt>=0) then
  call cprj_free(cwaveprj)
  deallocate(cwaveprj,dimlmn)
 end if

 mpi_enreg%paral_level= old_paral_level
 deallocate(sendcounts,recvcounts,sdispls,rdispls)
 deallocate(sendcountsloc,sdisplsloc)
 deallocate(recvcountsloc,rdisplsloc)
 deallocate(cwavef_alltoall,gvnlc_alltoall,gsc_alltoall)

end subroutine prep_nonlop
!!***
