!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_fourwf
!! NAME
!! prep_fourwf
!!
!! FUNCTION
!! this routine prepares the data to the call of fourwf.
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
!!  icall = order of call of this routine in lobpcgccwf
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
!!  vlocal(n4,n5,n6,nvloc)= local potential in real space, on the augmented fft grid
!!
!! OUTPUT
!!  gwavef=(2,npw*nspinor*ndat)=matrix elements <G|H|C>.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      mkrho,vtowfk
!!
!! CHILDREN
!!      fourwf,leave_new,prep_index_wavef_bandpp,prep_kg_sym_do
!!      prep_wavef_sym_do,prep_wavef_sym_undo,sphereboundary,timab,wrtout
!!      xallgather_mpi,xallgatherv_mpi,xalltoallv_mpi,xcomm_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prep_fourwf(rhoaug,blocksize,cwavef,wfraug,gs_hamk,istwf_k,iblock,icall,kg_k_gather,&
& mgfft,mpi_enreg,nbdblock,nband_k,dimtabs,npw_k,n4,n5,n6,occ_k,paral_kgb,wtk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_14wfs, except_this_one => prep_fourwf
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: blocksize,dimtabs,iblock,icall,istwf_k,mgfft,n4,n5,n6,nband_k
 integer :: nbdblock,npw_k,paral_kgb
 real(dp) :: wtk
 type(gs_hamiltonian_type) :: gs_hamk
 type(mpi_type) :: mpi_enreg
!arrays
 integer :: kg_k_gather(3,dimtabs)
 real(dp) :: cwavef(2,npw_k*blocksize),occ_k(nband_k),rhoaug(n4,n5,n6)
 real(dp) :: wfraug(2,n4,n5,n6)

!Local variables-------------------------------
!local variables for mpialltoallv
!local variable for bandpp and inversion by symetry of time
!scalars
 integer,save :: idatarecv0,ndatarecv_tot,ndatasend_sym
 integer :: bandpp,bandpp_sym,ier,iproc,ndatarecv,nproc_band,npw_tot
 integer :: old_me_g0=0,old_paral_level,oldspacecomm,spaceComm=0,tim_fourwf
 real(dp) :: weight
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
 real(dp) :: dummy(2,1),tsec(2)
 real(dp),allocatable :: cwavef_alltoall(:,:)
 real(dp),pointer :: ewavef_alltoall_sym(:,:)
!no_abirules
!correspondence with abinit. here for real wf but in complex mode
!this is the index of a given band

! *************************************************************************
!====================================================================================
 nproc_band = mpi_enreg%nproc_band
 bandpp     = mpi_enreg%bandpp

 flag_inv_sym = ((gs_hamk%ngfft(7)==401) .and. (istwf_k==2))

 if (flag_inv_sym) then
  istwf_k         = 1 

  if (modulo(bandpp,2)==0) then
   bandpp_sym   = bandpp/2
  else
   bandpp_sym   = bandpp
  end if
 end if
!====================================================================================


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
 call timab(548,1,tsec)
 call xallgather_mpi(npw_k,recvcounts,spaceComm,ier)
 call timab(548,2,tsec)
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
  
! !$     call xallgatherv_mpi&
! !$          &  (kg_k_gather_sym,3*ndatarecv_tot,kg_k_gather_all,3*npw_per_proc(:),&
! !$          &  3*rdispls_all,oldspaceComm,ier)
  
  call xallgatherv_mpi&
&  (kg_k_gather_sym(1,:),ndatarecv_tot,kg_k_gather_all(1,:),npw_per_proc(:),rdispls_all,oldspaceComm,ier)
  call xallgatherv_mpi&
&  (kg_k_gather_sym(2,:),ndatarecv_tot,kg_k_gather_all(2,:),npw_per_proc(:),rdispls_all,oldspaceComm,ier)
  call xallgatherv_mpi&
&  (kg_k_gather_sym(3,:),ndatarecv_tot,kg_k_gather_all(3,:),npw_per_proc(:),rdispls_all,oldspaceComm,ier)


  call sphereboundary(gs_hamk%gbound,istwf_k,kg_k_gather_all,mgfft,npw_tot)
  
  deallocate(kg_k_gather_all,npw_per_proc,rdispls_all)   
 end if
!===============================================================

 sendcounts(:)=npw_k*bandpp
 do iproc=1,nproc_band
  sdispls(iproc)=(iproc-1)*npw_k*bandpp
 end do

 allocate(cwavef_alltoall(2,ndatarecv*bandpp))
 
 recvcountsloc(:)=recvcounts(:)*2*bandpp
 rdisplsloc(:)=rdispls(:)*2*bandpp
 sendcountsloc(:)=sendcounts(:)*2
 sdisplsloc(:)=sdispls(:)*2

 call timab(547,1,tsec)
 call xalltoallv_mpi(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,&
& recvcountsloc,rdisplsloc,spaceComm,ier)
 call timab(547,2,tsec)

!If me_fft==0, I have the G=0 vector, but keep for the record the old value
 if (mpi_enreg%me_fft==0) then
  old_me_g0=mpi_enreg%me_g0
  mpi_enreg%me_g0=1
 end if

!Attention, c'est pour les essais. Le test est a remettre par la suite.
 if(abs(occ_k(mpi_enreg%coords(2)+1+(iblock-1)*blocksize)) >=tol8) then
  tim_fourwf=16
  weight=occ_k(mpi_enreg%coords(2)+1+(iblock-1)*blocksize)*wtk/gs_hamk%ucvol
  
! ====================================================================
  if ((.not.(flag_inv_sym)) .and. (bandpp==1)) then

   call fourwf(1,rhoaug,cwavef_alltoall,dummy,wfraug,&
&   gs_hamk%gbound,gs_hamk%gbound,&
&   istwf_k,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,1,&
&   gs_hamk%ngfft,ndatarecv,1,n4,n5,n6,1,paral_kgb,tim_fourwf,weight)
   
  else if ((.not.(flag_inv_sym)) .and. (bandpp>1)) then
   
!  -------------------------------------------------------------
!  Calcul of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   1,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)
   
!  -------------------------------------------------------
!  Classment of the waves functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)
   
!  -------------------
!  Fourier calculation
!  -------------------  
   call fourwf(1,rhoaug,cwavef_alltoall,dummy,wfraug,&
&   gs_hamk%gbound,gs_hamk%gbound,&
&   istwf_k,kg_k_gather,kg_k_gather,mgfft,mpi_enreg,bandpp,&
&   gs_hamk%ngfft,ndatarecv,1,n4,n5,n6,1,paral_kgb,tim_fourwf,weight)
   
!  -----------------------------------------------------
!  Classment of waves functions below the prossecors
!  -----------------------------------------------------
   cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)
   deallocate(index_wavef_band)
   

  else if (flag_inv_sym) then

!  -------------------------------------------------------------
!  Calcul of the index to class the waves functions below bandpp
!  -------------------------------------------------------------
   call prep_index_wavef_bandpp(nproc_band,bandpp,&
&   1,ndatarecv,&
&   recvcounts,rdispls,&
&   index_wavef_band)

!  -------------------------------------------------------
!  Classment of de the waves functions below bandpp
!  -------------------------------------------------------
   cwavef_alltoall(:,:) = cwavef_alltoall(:,index_wavef_band)

!  ------------------------------------------------------------
!  We associate the waves functions by two
!  ------------------------------------------------------------
   call prep_wavef_sym_do(mpi_enreg,nproc_band,bandpp,1,&
   ndatarecv,recvcounts,rdispls,&
   ndatarecv_tot,ndatasend_sym,tab_proc,&
   cwavef_alltoall,&
   sendcounts_sym,sendcounts_sym_all,sdispls_sym,&
   recvcounts_sym,recvcounts_sym_tot,rdispls_sym,&
   ewavef_alltoall_sym,&
   index_wavef_send)

!  ------------------------------------------------------------
!  Fourier calculcation
!  ------------------------------------------------------------
   call fourwf(1,rhoaug,ewavef_alltoall_sym,dummy,wfraug,&
&   gs_hamk%gbound,gs_hamk%gbound,&
&   istwf_k,kg_k_gather_sym,kg_k_gather_sym,mgfft,mpi_enreg,bandpp_sym,&
&   gs_hamk%ngfft,ndatarecv_tot,1,n4,n5,n6,1,paral_kgb,tim_fourwf,weight)

!  ------------------------------------------------------------
!  We dissociate each wave function in two waves functions
!  gwavef is classed below of bandpp
!  ------------------------------------------------------------
   call prep_wavef_sym_undo(mpi_enreg,nproc_band,bandpp,1,&
   ndatarecv,recvcounts,rdispls,&
   ndatarecv_tot,ndatasend_sym,idatarecv0,tab_proc,&
   cwavef_alltoall,&
   sendcounts_sym,sendcounts_sym_all,sdispls_sym,&
   recvcounts_sym,recvcounts_sym_tot,rdispls_sym,&
   ewavef_alltoall_sym,&
   index_wavef_send)

   deallocate(ewavef_alltoall_sym)
   
!  -------------------------------------------------------
!  Classment of waves functions below the processors
!  ------------------------------------------------------- 
   cwavef_alltoall(:,index_wavef_band) = cwavef_alltoall(:,:)
   
   deallocate(index_wavef_band)
   
  end if
! ====================================================================
 end if


!====================================================================
 if (flag_inv_sym) then
  istwf_k         = 2 
! gs_hamk%istwf_k = 2
 end if
!====================================================================

!if (gs_hamk%istwf_k==2) mpi_enreg%me_g0=old_me_g0
 if (mpi_enreg%me_fft==0) mpi_enreg%me_g0=old_me_g0
 mpi_enreg%paral_level= old_paral_level
 deallocate(recvcounts,sdispls,rdispls)
 deallocate(cwavef_alltoall)
 deallocate(sendcountsloc,sdisplsloc)
 deallocate(recvcountsloc,rdisplsloc)
end subroutine prep_fourwf
!!***
