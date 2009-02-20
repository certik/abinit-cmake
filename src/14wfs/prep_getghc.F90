!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_getghc
!! NAME
!! prep_getghc
!!
!! FUNCTION
!! this routine prepares the data to the call of getghc.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FBottin,MT,GZ)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  blocksize= size of block for FFT
!!  cwavef(2,npw*nspinor*ndat)=planewave coefficients of wavefunction.
!!  dimffnl=second dimension of ffnl (1+number of derivatives)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  ffnl(npw,dimffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the hamiltonian at k
!!  gvnlc=matrix elements <G|Vnonlocal|C>
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  icall = order of call of this routine in lobpcgccwf
!!  lambda=factor to be used when computing <G|H-lambda.S|C> - only for sij_opt=-1
!!         Typically lambda is the eigenvalue (or its guess)
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all type of psps
!!        =if useylm=0, max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mgfft=maximum size of 1d ffts
!!  mpi_enreg=informations about mpi parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nbdblock=
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  ntypat=number of types of atoms in unit cell.
!!  nvloc=final dimension of vlocal (usually 1, but 4 for non-collinear)
!!  n4,n5,n6 used for dimensionning of vlocal
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  sij_opt= -PAW ONLY-  if  0, only matrix elements <G|H|C> have to be computed
!!     (S=overlap)       if  1, matrix elements <G|S|C> have to be computed in gsc in addition to ghc
!!                       if -1, matrix elements <G|H-lambda.S|C> have to be computed in ghc (gsc not used)
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!
!! OUTPUT
!!  gwavef=(2,npw*nspinor*ndat)=matrix elements <G|H|C> (if sij_opt>=0)
!!                                  or <G|H-lambda.S|C> (if sij_opt=-1).
!!  swavef=(2,npw*nspinor*ndat)=matrix elements <G|S|C>.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      lobpcgccwf,lobpcgwf
!!
!! CHILDREN
!!      fourwf,sphereboundary,timab,xallgather_mpi,xallgatherv_mpi
!!      xalltoallv_mpi,xcomm_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prep_getghc(cwavef,dimffnl,dtfil,ffnl_gather,gs_hamk,gvnlc,gwavef,swavef,iblock,icall,istwf_k,kg_k_gather,&
& kinpw_gather,lambda,lmnmax,matblk,blocksize,mgfft,mpi_enreg,mpsang,mpssoang,natom,nbdblock,nband_k,dimtabs,npw_k,&
& nspinor,ntypat,nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,vlocal)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_14wfs, except_this_one => prep_getghc
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: blocksize,dimffnl,dimtabs,iblock,icall,istwf_k,lmnmax,matblk,mgfft
 integer :: mpsang,mpssoang,n4,n5,n6,natom,nband_k,nbdblock,npw_k,nspinor
 integer :: ntypat,nvloc,paral_kgb,prtvol,sij_opt
 real(dp) :: lambda
 type(datafiles_type) :: dtfil
 type(gs_hamiltonian_type) :: gs_hamk
 type(mpi_type) :: mpi_enreg
!arrays
 integer :: kg_k_gather(3,dimtabs)
 real(dp) :: cwavef(2,npw_k*nspinor*blocksize)
 real(dp) :: ffnl_gather(dimtabs,dimffnl,lmnmax,ntypat)
 real(dp) :: gvnlc(2,npw_k*nspinor*blocksize),gwavef(2,npw_k*nspinor*blocksize)
 real(dp) :: kinpw_gather(dimtabs),ph3d_gather(2,dimtabs,matblk)
 real(dp) :: swavef(2,npw_k*nspinor*blocksize),vlocal(n4,n5,n6,nvloc)

!Local variables-------------------------------
!local variables for mpialltoallv
!local variable for bandpp and inversion by symetry of time
!scalars
 integer,save :: idatarecv0,ndatarecv_tot,ndatasend_sym
 integer :: bandpp,bandpp_sym,ier,iproc,ndatarecv,nproc_band,npw_tot,old_me_g0
 integer :: old_paral_level,oldspacecomm,spaceComm=0,tim_getghc
 logical :: flag_inv_sym
 character(len=500) :: message
!arrays
 integer,pointer,save :: kg_k_gather_sym(:,:),rdispls_sym(:),recvcounts_sym(:)
 integer,pointer,save :: recvcounts_sym_tot(:),sdispls_sym(:),sendcounts_sym(:)
 integer,pointer,save :: sendcounts_sym_all(:),tab_proc(:)
 integer,allocatable :: kg_k_gather_all(:,:),npw_per_proc(:),rdispls(:)
 integer,allocatable :: rdispls_all(:),rdisplsloc(:),recvcounts(:)
 integer,allocatable :: recvcountsloc(:),sdispls(:),sdisplsloc(:),sendcounts(:)
 integer,allocatable :: sendcountsloc(:)
 integer,pointer :: index_wavef_band(:),index_wavef_send(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef_alltoall(:,:),gvnlc_alltoall(:,:)
 real(dp),allocatable :: gwavef_alltoall(:,:),swavef_alltoall(:,:)
 real(dp),allocatable :: tmp_ffnl_gather(:,:,:,:),tmp_kinpw_gather(:)
 real(dp),allocatable :: tmp_ph3d_gather(:,:,:)
 real(dp),pointer :: ewavef_alltoall_sym(:,:),gvnlc_alltoall_sym(:,:)
 real(dp),pointer :: gwavef_alltoall_sym(:,:),swavef_alltoall_sym(:,:)
!no_abirules
!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band

! *************************************************************************
!====================================================================================
 nproc_band = mpi_enreg%nproc_band
 bandpp     = mpi_enreg%bandpp

 flag_inv_sym = ((gs_hamk%ngfft(7)==401) .and. (gs_hamk%istwf_k==2))

 if (flag_inv_sym) then
  istwf_k         = 1 
  gs_hamk%istwf_k = 1
  if (modulo(bandpp,2)==0) then
   bandpp_sym   = bandpp/2
  else
   bandpp_sym   = bandpp
  end if
 end if
!====================================================================================

 tim_getghc=6
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
 call timab(546,1,tsec)
 call xallgather_mpi(npw_k,recvcounts,spaceComm,ier)
 call timab(546,2,tsec)
 rdispls(1)=0
 do iproc=2,nproc_band
  rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
 end do
 ndatarecv=rdispls(nproc_band)+recvcounts(nproc_band)

!FB The transposition are no more performed within prep_fourwf but above.
!FB We have to check that the dimension mpw is identical the the one computed above. TO BE CLEANED
 if (ndatarecv /= dimtabs) then
  write(message, '(a,a,a,a,i8,a,i8)' ) ch10,&
&  ' prep_fourwf: BUG -',ch10,&
&  '  ndatarecv',ndatarecv,' /= dimtabs ', dimtabs
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

!===============================================================
 if (icall==1 .and. iblock==1 .and. flag_inv_sym ) then
  
! --------------------------------------------------------
! We complete kg_k_gather vectors with the opposite values
! --------------------------------------------------------
  call prep_kg_sym_do(mpi_enreg,&
  kg_k_gather,ndatarecv,&
  kg_k_gather_sym,ndatarecv_tot,&
  ndatasend_sym,idatarecv0,&
  tab_proc,&
  sendcounts_sym,sendcounts_sym_all,sdispls_sym,&
  recvcounts_sym,recvcounts_sym_tot,rdispls_sym)
  
! --------------------------------------------------------
! We calcul gbound
! --------------------------------------------------------
  oldspacecomm=mpi_enreg%comm_fft
  allocate(npw_per_proc(mpi_enreg%nproc_fft),rdispls_all(mpi_enreg%nproc_fft))    
  call xallgather_mpi(ndatarecv_tot,npw_per_proc,oldspacecomm,ier)
  
  rdispls_all(1)=0
  do iproc=2,mpi_enreg%nproc_fft
   rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
  end do
  npw_tot=rdispls_all(mpi_enreg%nproc_fft)+npw_per_proc(mpi_enreg%nproc_fft)
  allocate(kg_k_gather_all(3,npw_tot))
  
  call xallgatherv_mpi(&
&  kg_k_gather_sym,3*ndatarecv_tot,kg_k_gather_all,3*npw_per_proc(:),&
&  3*rdispls_all,oldspaceComm,ier)
  
  call sphereboundary(gs_hamk%gbound,istwf_k,kg_k_gather_all,mgfft,npw_tot)
  
  deallocate(kg_k_gather_all,npw_per_proc,rdispls_all)   
 end if
!===============================================================

 sendcounts(:)=npw_k*bandpp
 do iproc=1,nproc_band
  sdispls(iproc)=(iproc-1)*npw_k*bandpp
 end do

 allocate(cwavef_alltoall(2,ndatarecv*nspinor*bandpp))
 allocate(gwavef_alltoall(2,ndatarecv*nspinor*bandpp))
 allocate(swavef_alltoall(2,ndatarecv*nspinor*bandpp))
 allocate(gvnlc_alltoall(2,ndatarecv*nspinor*bandpp))

 swavef_alltoall(:,:)=0.0

 recvcountsloc(:)=recvcounts(:)*2*nspinor*bandpp
 rdisplsloc(:)=rdispls(:)*2*nspinor*bandpp
 sendcountsloc(:)=sendcounts(:)*2*nspinor
 sdisplsloc(:)=sdispls(:)*2*nspinor
 call timab(545,1,tsec)
 call xalltoallv_mpi(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,&
& recvcountsloc,rdisplsloc,spaceComm,ier)
 call timab(545,2,tsec)

 if(gs_hamk%istwf_k==2) then
  old_me_g0=mpi_enreg%me_g0
  if (mpi_enreg%me_fft==0) then
   mpi_enreg%me_g0=1
  else
   mpi_enreg%me_g0=0
  end if
 end if
!====================================================================
 if ((.not.(flag_inv_sym)) .and. (bandpp==1)) then

  call getghc(cwavef_alltoall,dimffnl,ffnl_gather,dtfil%filstat,gwavef_alltoall,&
&  swavef_alltoall(:,:),gs_hamk,gvnlc_alltoall,kg_k_gather,&
&  kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,1,&
&  ndatarecv,nspinor,ntypat,&
&  nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal)
  
 else if ((.not.(flag_inv_sym)) .and. (bandpp>1)) then
  
! -------------------------------------------------------------
! Calcul of the index to class the waves functions below bandpp
! -------------------------------------------------------------
  call prep_index_wavef_bandpp(nproc_band,bandpp,&
&  nspinor,ndatarecv,&
&  recvcounts,rdispls,&
&  index_wavef_band)
  
! -------------------------------------------------------
! Classment of the waves functions below bandpp
! -------------------------------------------------------
  cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)
  
! -------------------
! Fourier calculation
! -------------------  
  call getghc(cwavef_alltoall,dimffnl,ffnl_gather,dtfil%filstat,gwavef_alltoall,&
&  swavef_alltoall(:,:),gs_hamk,gvnlc_alltoall,kg_k_gather,&
&  kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,bandpp,&
&  ndatarecv,nspinor,ntypat,&
&  nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,0,vlocal)       
  
! -----------------------------------------------------
! Classment of waves functions below the prossecors
! -----------------------------------------------------
  cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)
  if (sij_opt==1) then
   swavef_alltoall(:,index_wavef_band) = swavef_alltoall(:,:)
  end if
  gwavef_alltoall(:,index_wavef_band) = gwavef_alltoall(:,:)
  gvnlc_alltoall(:,index_wavef_band)  = gvnlc_alltoall(:,:)
  
  deallocate(index_wavef_band)
  

 else if (flag_inv_sym) then

! -------------------------------------------------------------
! Calcul of the index to class the waves functions below bandpp
! -------------------------------------------------------------
  call prep_index_wavef_bandpp(nproc_band,bandpp,&
&  nspinor,ndatarecv,&
&  recvcounts,rdispls,&
&  index_wavef_band)

! -------------------------------------------------------
! Classment of de the waves functions below bandpp
! -------------------------------------------------------
  cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)

! ------------------------------------------------------------
! We associate the waves functions by two
! ------------------------------------------------------------
  call prep_wavef_sym_do(mpi_enreg,nproc_band,bandpp,nspinor,&
  ndatarecv,recvcounts,rdispls,&
  ndatarecv_tot,ndatasend_sym,tab_proc,&
  cwavef_alltoall,&
  sendcounts_sym,sendcounts_sym_all,sdispls_sym,&
  recvcounts_sym,recvcounts_sym_tot,rdispls_sym,&
  ewavef_alltoall_sym,&
  index_wavef_send)

! ------------------------------------------------------------
! Allocation
! ------------------------------------------------------------
  allocate(gwavef_alltoall_sym     (2,ndatarecv_tot*bandpp_sym))
  allocate(swavef_alltoall_sym     (2,ndatarecv_tot*bandpp_sym))
  allocate(gvnlc_alltoall_sym      (2,ndatarecv_tot*bandpp_sym))
  
  allocate(tmp_ffnl_gather (ndatarecv_tot,dimffnl,lmnmax,ntypat))
  allocate(tmp_kinpw_gather(ndatarecv_tot))
  allocate(tmp_ph3d_gather (2,ndatarecv_tot,matblk))
  
! ------------------------------------------------------------
! Fourier calculcation
! ------------------------------------------------------------
  call getghc(ewavef_alltoall_sym(:,:),dimffnl,tmp_ffnl_gather,dtfil%filstat,&
&  gwavef_alltoall_sym,swavef_alltoall_sym(:,:),gs_hamk,gvnlc_alltoall_sym,&
&  kg_k_gather_sym,&
&  tmp_kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,&
&  bandpp_sym,ndatarecv_tot,nspinor,ntypat,&
&  nvloc,n4,n5,n6,paral_kgb,tmp_ph3d_gather,prtvol,sij_opt,tim_getghc,1,vlocal)
  
  deallocate(tmp_ffnl_gather,tmp_kinpw_gather,tmp_ph3d_gather)
  
! ------------------------------------------------------------
! We dissociate each wave function in two waves functions
! gwavef is classed below of bandpp
! ------------------------------------------------------------
  call prep_wavef_sym_undo(mpi_enreg,nproc_band,bandpp,nspinor,&
  ndatarecv,recvcounts,rdispls,&
  ndatarecv_tot,ndatasend_sym,idatarecv0,tab_proc,&
  gwavef_alltoall,&
  sendcounts_sym,sendcounts_sym_all,sdispls_sym,&
  recvcounts_sym,recvcounts_sym_tot,rdispls_sym,&
  gwavef_alltoall_sym,&
  index_wavef_send)

  deallocate(ewavef_alltoall_sym)
  deallocate(gwavef_alltoall_sym)
  deallocate(swavef_alltoall_sym)
  deallocate(gvnlc_alltoall_sym)
  
! -------------------------------------------
! We call getghc for the calcul of ffnl,...
! --------------------------------------------
! tim_getghc=8
  gs_hamk%istwf_k=2
  
  old_me_g0=mpi_enreg%me_g0
  if (mpi_enreg%me_fft==0) then
   mpi_enreg%me_g0=1
  else
   mpi_enreg%me_g0=0
  end if
  
  call getghc(cwavef_alltoall,dimffnl,ffnl_gather,dtfil%filstat,gwavef_alltoall,&
&  swavef_alltoall(:,:),gs_hamk,gvnlc_alltoall,kg_k_gather,&
&  kinpw_gather,lambda,lmnmax,matblk,mgfft,mpi_enreg,mpsang,mpssoang,natom,bandpp,&
&  ndatarecv,nspinor,ntypat,&
&  nvloc,n4,n5,n6,paral_kgb,ph3d_gather,prtvol,sij_opt,tim_getghc,2,vlocal)
  
  mpi_enreg%me_g0=old_me_g0
  
  gs_hamk%istwf_k=1
  
! -------------------------------------------------------
! Classment of waves functions below the processors
! ------------------------------------------------------- 
  cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)
  
  gwavef_alltoall(:,index_wavef_band) = gwavef_alltoall(:,:)
  if (sij_opt==1) then
   swavef_alltoall(:,index_wavef_band) = swavef_alltoall(:,:)
  end if
  
  gvnlc_alltoall(:,index_wavef_band)  = gvnlc_alltoall(:,:) 
  
  deallocate(index_wavef_band)
  
 end if
!====================================================================
 
 if (gs_hamk%istwf_k==2) mpi_enreg%me_g0=old_me_g0
 call timab(545,1,tsec)
 if (sij_opt==1) then
  call xalltoallv_mpi(swavef_alltoall,recvcountsloc,rdisplsloc,swavef,&
&  sendcountsloc,sdisplsloc,spaceComm,ier)
 end if
 call xalltoallv_mpi(gwavef_alltoall,recvcountsloc,rdisplsloc,gwavef,&
& sendcountsloc,sdisplsloc,spaceComm,ier)
 call xalltoallv_mpi(gvnlc_alltoall,recvcountsloc,rdisplsloc,gvnlc,&
& sendcountsloc,sdisplsloc,spaceComm,ier)
 call timab(545,2,tsec)

!====================================================================
 if (flag_inv_sym) then
  istwf_k         = 2 
  gs_hamk%istwf_k = 2
 end if
!====================================================================

 mpi_enreg%paral_level= old_paral_level
 deallocate(sendcounts,recvcounts,sdispls,rdispls)
 deallocate(sendcountsloc,sdisplsloc)
 deallocate(recvcountsloc,rdisplsloc)
 deallocate(cwavef_alltoall,gwavef_alltoall,gvnlc_alltoall)
 deallocate(swavef_alltoall)
end subroutine prep_getghc
!!***
