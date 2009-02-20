!{\src2tex{textfont=tt}}
!!****f* ABINIT/susk
!! NAME
!! susk
!!
!! FUNCTION
!! Compute the contribution of one k point to the susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!! Compared to the routine suskmm, there is no particular attention
!! to the use of the memory, so the code is simpler.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  bdtot_index=index for the number of the band
!!  cg(2,mcg)=wfs in G space
!!  cprj_k(natom,nspinor*nband_k)= wave functions projected with non-local projectors:
!!                                 cprj_k=<p_i|Cnk> where p_i is a non-local projector.
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  extrap: if==1, the closure relation (an extrapolation) must be used
!!  gbound(2*mgfftdiel+8,2)=G sphere boundary for going from WF sphere to
!!      medium size FFT grid
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for going from medium size
!!      FFT grid to small sphere.
!!  gylmg_diel(npwdiel,lmax_diel,ntypat*usepaw)= -PAW only- Fourier transform of g_l(r).Y_ml(r) shape functions
!!                                               for dielectric matrix
!!  icg=index for cg
!!  ikpt=number of the k point
!!  isp=number of the current spin
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kg_k(3,npw_k)=coordinates of planewaves in basis sphere.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mband=maximum number of bands
!!  mcg=dimension of cg
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of
!!     the dielectric matrix
!!  mkmem=maximum number of k points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  natom=number of atoms in cell
!!  nband_k=number of bands at this k point for that spin polarization
!!  ndiel4,ndiel5,ndiel6= FFT dimensions, modified to avoid cache trashing
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwdiel=third and fifth dimension of the susmat array.
!!  npw_k=number of plane waves at this k point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  occ_deavg(mband)=factor for extrapolation (occup. divided by an energy gap)
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph3d_diel(2,npwdiel,natom*usepaw)=3-dim structure factors, for each atom and plane wave, for dielectric matrix
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume (Bohr**3)
!!  usepaw=flag for PAW
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! These quantities are accumulated in this routine:
!!  drhode(2,npwdiel,nsppol)=weighted density, needed to compute the
!!   effect of change of fermi energy
!!  rhoextrap(ndiel4,ndiel5,ndiel6)=density-like array, needed for the
!!   extrapolation procedure.
!!  sumdocc=sum of weighted occupation numbers, needed to compute the
!!   effect of change of fermi energy
!!  susmat(2,npwdiel,nsppol,npwdiel,nsppol)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! NOTES
!! Band-fft parallel treatment: Each processor will treat his own band, but susmat will be known by all.
!! This means that cg will not have the same meaning in sequential or parallel mode. 
!! In parallel mode, it will contain the set of all bands treated by the currrent processor.
!! To achieve this, the argument cg has been replaced by cg_mpi, with the "target" attribute.
!! In sequential mode, the pointer cg will point towards cg_mpi. In parallel mode, cg will point
!! to a new array cg_local, containing the bands treated by the currrent processor.
!! This allows to minimize the overhead incurred by the parallelization  of the sequential version.
!! A similar treatment is performed on kg_k, npw_k.
!! A future version might have objects like kg_k_gather as arguments, instead of computing them.
!! This is in slight violation of programming rules, but I think it is safe, since the pointers remain local
!! GZ
!! PARENTS
!!      suscep_stat
!!
!! CHILDREN
!!      fourwf,pawsushat,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine susk(atindx1,bdtot_index,cg_mpi,cprj_k,doccde,drhode,eigen,extrap,gbound,&
     &  gbound_diel,gylmg_diel,icg_mpi,ikpt,isp,istwfk,kg_diel,kg_k_mpi,&
     &  lmax_diel,mband,mcg,mgfftdiel,mkmem,mpi_enreg,mpw,&
     &  natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,ngfftdiel,nkpt,&
     &  npwdiel,npw_k_mpi,nspden,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&
     &  pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
     &  susmat,typat,ucvol,usepaw,wtk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_12ffts
 use interfaces_13paw
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: bdtot_index,extrap,ikpt,isp,lmax_diel,mband,mcg,mgfftdiel
 integer,intent(in) :: mkmem,mpw,natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,nkpt
 integer,intent(in),target :: icg_mpi,npw_k_mpi
 integer,intent(in) :: npwdiel,nspden,nspinor,nsppol,ntypat,occopt,usepaw
 real(dp),intent(in) :: ucvol
 real(dp),intent(inout) :: sumdocc
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
!arrays
 integer,intent(inout) :: gbound(2*mgfftdiel+8,2)
 integer,intent(in) :: atindx1(natom)
 integer,intent(in) ::gbound_diel(2*mgfftdiel+8,2)
 integer,intent(in) :: istwfk(nkpt),kg_diel(3,npwdiel)
 integer,intent(in) :: ngfftdiel(18),typat(natom)
 integer,intent(in), target ::kg_k_mpi(3,npw_k_mpi)
 integer,dimension(:,:),pointer   ::kg_k
 real(dp),intent(in),target :: cg_mpi(2,mcg)
 real(dp),intent(in) :: doccde(mband*nkpt*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),gylmg_diel(npwdiel,lmax_diel**2,ntypat*usepaw)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),occ_deavg(mband),ph3d_diel(2,npwdiel,natom*usepaw),wtk(nkpt)
 real(dp),intent(inout) :: drhode(2,npwdiel,nsppol)
 real(dp),intent(inout) :: rhoextrap(ndiel4,ndiel5,ndiel6)
 real(dp),intent(inout) :: susmat(2,npwdiel,nsppol,npwdiel,nsppol)
 type(cprj_type) :: cprj_k(natom,nspinor*nband_k*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
! real(dp), allocatable :: cg_disk(:,:)
!scalars
 integer :: i1,i2,i3,iband,iblock,ibdblock,ibd1,ibd2,ipw,ipw1,ipw2,isp1,isp2,istwf_k,nbdblock,ndiel1
 integer :: ndiel2,ndiel3,paral_kgb_diel,testocc,tim_fourwf
 integer,target :: icg_loc=0,npw_k_loc
 real(dp) :: ai,ar,eigdiff,norm,normr,occdiff,tolocc,weight,wght1,wght2
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: cwavef(:,:),cwavef_alltoall(:,:),cwavef_alltoall_gather(:,:),dummy(:,:),rhoaug(:,:,:),wfprod(:,:)
 real(dp),allocatable,target ::cg_local(:,:)
 real(dp),allocatable :: wfraug(:,:,:,:),wfrspa(:,:,:,:,:)
 integer,pointer :: npw_k,icg
 real(dp),pointer :: cg(:,:)
!Local variables for MPI
 type(MPI_type) :: mpi_enreg_diel
 integer,allocatable :: band_loc(:),kg_k_gather(:,:),npw_per_proc(:),rdispls(:),rdispls_all(:),rdisplsloc(:),recvcounts(:)
 integer,allocatable :: recvcountsloc(:),sdispls(:),sdisplsloc(:),sendcounts(:),sendcountsloc(:)
 integer,allocatable,target :: kg_k_gather_all(:,:)
 logical,allocatable :: treat_band(:)
 integer :: blocksize,ier,iproc,ndatarecv,spaceComm,iband_loc,sizemax_per_proc,iproc_fft
 integer, target:: npw_tot

! *************************************************************************

!DEBUG
!write(6,*)' susk : enter '
!if(.true.)stop
!ENDDEBUG

 call timab(87,1,tsec)

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)
 istwf_k=istwfk(1)

!The dielectric stuff is performed in sequential mode.
!Set mpi_enreg_diel accordingly
 mpi_enreg_diel%paral_compil_fft=0
 mpi_enreg_diel%paral_compil_kpt=0
 mpi_enreg_diel%mode_para='n'
 mpi_enreg_diel%me=0
 mpi_enreg_diel%me_fft=0
 mpi_enreg_diel%me_kpt=0
 mpi_enreg_diel%nproc_fft=1
 mpi_enreg_diel%fft_option_lob=0
 paral_kgb_diel=0

 testocc=1
!DEBUG
!write(6,*)' susk : set testocc to 0 '
!testocc=0
!write(6,*)' susk : set extrap to 0 '
!extrap=0
!ENDDEBUG

 allocate(cwavef(2,mpw),dummy(2,1))
 allocate(rhoaug(ndiel4,ndiel5,ndiel6),wfraug(2,ndiel4,ndiel5,ndiel6))
 allocate(wfprod(2,npwdiel),wfrspa(2,ndiel4,ndiel5,ndiel6,mband))
 wfrspa(:,:,:,:,:)=zero
 allocate(treat_band(nband_k))
 treat_band(:)=.true.
 
!BAND-FFT parallelism 
 if (mpi_enreg%mode_para=='b') then
  treat_band(:)=.false.
! We gather the wavefunctions treated by this  proc in cg_local
  spaceComm=mpi_enreg%comm_band
  blocksize=mpi_enreg%nproc_band
  nbdblock=nband_k/blocksize
  allocate(sdispls(blocksize))
  allocate(sdisplsloc(blocksize))
  allocate(sendcounts(blocksize))
  allocate(sendcountsloc(blocksize))
  allocate(rdispls(blocksize))
  allocate(rdisplsloc(blocksize))
  allocate(recvcounts(blocksize))
  allocate(recvcountsloc(blocksize))
! First gather the kg_k in kg_k_gather_all  
  npw_k_loc=npw_k_mpi
  call xallgather_mpi(npw_k_loc,recvcounts,spaceComm,ier)
  rdispls(1)=0
  do iproc=2,blocksize
   rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
  end do
  ndatarecv=rdispls(blocksize)+recvcounts(blocksize)
  allocate(kg_k_gather(3,ndatarecv))
  recvcountsloc(:)=recvcounts(:)*3
  rdisplsloc(:)=rdispls(:)*3
  call xallgatherv_mpi(kg_k_mpi,3*npw_k_loc,kg_k_gather,recvcountsloc(:),rdisplsloc,spaceComm,ier)
  allocate(npw_per_proc(mpi_enreg%nproc_fft),rdispls_all(mpi_enreg%nproc_fft))
  spaceComm=mpi_enreg%comm_fft
  call xallgather_mpi(ndatarecv,npw_per_proc,spaceComm,ier)
  rdispls_all(1)=0
  do iproc=2,mpi_enreg%nproc_fft
   rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
  end do
  npw_tot=rdispls_all(mpi_enreg%nproc_fft)+npw_per_proc(mpi_enreg%nproc_fft)
  allocate(kg_k_gather_all(3,npw_tot))
  call xallgatherv_mpi(kg_k_gather,3*ndatarecv,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,spaceComm,ier)
! At this point kg_k_gather_all contains all the kg   
  if(allocated(cwavef)) deallocate(cwavef)
  allocate(cwavef(2,npw_k_loc*nspinor*blocksize))
  sizemax_per_proc=nband_k/(mpi_enreg%nproc_band*mpi_enreg%nproc_fft)+1
  allocate(band_loc(sizemax_per_proc))
  allocate(cg_local(2,sizemax_per_proc*npw_tot))
  iband_loc=0
  do ibdblock=1,nbdblock
   cwavef(:,1:npw_k_loc*nspinor*blocksize)=&
&   cg_mpi(:,1+(ibdblock-1)*npw_k_loc*nspinor*blocksize+icg_mpi:ibdblock*npw_k_loc*nspinor*blocksize+icg_mpi)
   sendcounts(:)=npw_k_loc
   do iproc=1,blocksize
    sdispls(iproc)=(iproc-1)*npw_k_loc
   end do
   allocate(cwavef_alltoall(2,ndatarecv))
   recvcountsloc(:)=recvcounts(:)*2
   rdisplsloc(:)=rdispls(:)*2
   sendcountsloc(:)=sendcounts(:)*2
   sdisplsloc(:)=sdispls(:)*2
   call timab(547,1,tsec)
   spaceComm=mpi_enreg%comm_band
   call xalltoallv_mpi(cwavef,sendcountsloc,sdisplsloc,cwavef_alltoall,recvcountsloc,rdisplsloc,spaceComm,ier)
   call timab(547,2,tsec)
   call timab(538,1,tsec)
   allocate(cwavef_alltoall_gather(2,npw_tot))
   blocksize=mpi_enreg%nproc_band
   spaceComm=mpi_enreg%comm_fft
   call xallgatherv_mpi(cwavef_alltoall,2*ndatarecv,cwavef_alltoall_gather,2*npw_per_proc,2*rdispls_all,spaceComm,ier)
   iproc_fft=modulo(ibdblock-1,mpi_enreg%nproc_fft)
   if(mpi_enreg%me_fft==iproc_fft) then !All nproc_band procs of index me_fft will treat these bands
    iband_loc=iband_loc+1
    iband=1+mpi_enreg%me_band+mpi_enreg%nproc_band*mpi_enreg%me_fft+(iband_loc-1)*mpi_enreg%nproc_fft*mpi_enreg%nproc_band
    treat_band(iband)=.true.
    band_loc(iband_loc)=iband
    cg_local(:,1+(iband_loc-1)*npw_tot:iband_loc*npw_tot)=cwavef_alltoall_gather(:,1:npw_tot)
   end if
   deallocate(cwavef_alltoall_gather, cwavef_alltoall)
  end do
! On exit:
! npw_tot will be npw
! kg_k_gather_all will be kg_k
! cg_local will be cg
! icg will be zero
  npw_k=>npw_tot
  kg_k=>kg_k_gather_all(:,:)
  cg=>cg_local(:,:)
  icg=>icg_loc
  call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)
  deallocate(npw_per_proc,rdispls_all)
  deallocate(sendcounts,recvcounts,sdispls,rdispls)
  deallocate(sendcountsloc,sdisplsloc)
  deallocate(recvcountsloc,rdisplsloc)
  deallocate(kg_k_gather)
! Because they will be summed over all procs, and arrive on input, rescale drhode and rhoextrap
  if(occopt>=3)drhode(:,:,:)=drhode(:,:,:)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
  if(extrap==1)rhoextrap(:,:,:)=rhoextrap(:,:,:)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
  susmat(:,:,:,:,:)=susmat(:,:,:,:,:)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)
  
! No BAND-FFT parallelism 
 else ! use argument variables
  cg=>cg_mpi
  kg_k=>kg_k_mpi
  npw_k=>npw_k_mpi
  icg=>icg_mpi
 end if
 iband_loc=0
 
!Loop over bands to fft and store Fourier transform of wavefunction
 do iband=1,nband_k
  if(.not. treat_band(iband))  cycle ! I am not treating this band (only for the parallel case)
  iband_loc=iband_loc+1
! if(mask_paral(iband)==.true.) cycle
! Obtain Fourier transform in fft box
  if(allocated(cwavef)) deallocate(cwavef)
  allocate(cwavef(2,npw_k))
  cwavef(:,1:npw_k)=cg(:,1+(iband_loc-1)*npw_k+icg:iband_loc*npw_k+icg)
! DEBUG
! write(6,*)' susk : will call fourwf for band',iband
! ENDDEBUG
! FB Transpose the kg_k array. TO BE CLEANED
! do iblock=1,npw_k; write(6,*) 'cwavef=',cwavef(1,iblock),kg_k(1:3,iblock) ; enddo
  tim_fourwf=8
  call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&  istwf_k,kg_k,kg_k,mgfftdiel,mpi_enreg_diel,1,ngfftdiel,npw_k,1,ndiel4,ndiel5,ndiel6,&
&  0,paral_kgb_diel,tim_fourwf,weight)
  wfrspa(:,:,:,:,iband)=wfraug(:,:,:,:)
  if( (occopt>=3 .and. testocc==1) .or. extrap==1 )then
!  In the case of metallic occupation, or if the extrapolation
!  over higher bands is included, must compute the
!  Fourier transform of the density of each band, then
!  generate the part of the susceptibility matrix due
!  varying occupation numbers.

   weight=-two*occ_deavg(iband)*wtk(ikpt)/ucvol
!  DEBUG
!  write(6,*)' susk : debug one band contribution '
!  weight=zero
!  if(iband==1)weight=-two*occ_deavg(iband)*wtk(ikpt)/ucvol
!  ENDDEBUG
   do i3=1,ndiel3
    do i2=1,ndiel2
     do i1=1,ndiel1
      wfraug(1,i1,i2,i3)=wfraug(1,i1,i2,i3)**2+wfraug(2,i1,i2,i3)**2
      wfraug(2,i1,i2,i3)=zero
     end do
    end do
!   If extrapolation, accumulate density in real space
    if(extrap==1.and.usepaw==0)then
     do i2=1,ndiel2
      do i1=1,ndiel1
       rhoextrap(i1,i2,i3)=rhoextrap(i1,i2,i3)+weight*wfraug(1,i1,i2,i3)
      end do
     end do
    end if
   end do

!  In case of PAW, add compensation charge contribution
   if (usepaw==1.and.extrap==1) then
    call pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,iband,iband,istwf_k,kg_diel,&
&    lmax_diel,mgfftdiel,mpi_enreg_diel,natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,&
&    ngfftdiel,npwdiel,nspinor,ntypat,1,paral_kgb_diel,&
&    pawang,pawtab,ph3d_diel,typat,dummy,wfraug)
    rhoextrap(:,:,:)=rhoextrap(:,:,:)+weight*wfraug(1,:,:,:)
   end if

!  DEBUG
!  if(iband==1)then
!  do i3=1,ndiel3,4
!  write(6,*)1,1,i3,wfraug(1,1,1,i3)
!  end do
!  call leave_new('COLL')
!  ENDDEBUG

!  Performs the Fourier Transform of the density of the band,
!  and store it in wfprod
   tim_fourwf=9
   call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&   istwf_k,kg_diel,kg_diel,&
&   mgfftdiel,mpi_enreg_diel,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb_diel,tim_fourwf,weight)
!  In case of PAW, add compensation charge contribution if not already done
   if (usepaw==1.and.extrap==0) then
    call pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,istwf_k,kg_diel,&
&    lmax_diel,mgfftdiel,mpi_enreg_diel,natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,&
&    ngfftdiel,npwdiel,nspinor,ntypat,0,paral_kgb_diel,&
&    pawang,pawtab,ph3d_diel,typat,wfprod,dummy)
   end if

!  Perform now the summation of terms related to direct change of eigenvalues
!  or extrapolation over higher bands
   wght1=zero ; wght2=zero
   if(occopt>=3 .and. testocc==1)then
    wght1=doccde(iband+bdtot_index)*wtk(ikpt)/ucvol
   end if

   if(extrap==1) wght2=two*occ_deavg(iband)*wtk(ikpt)/ucvol
!  DEBUG
!  write(6,*)' susk : debug '
!  if(extrap==1 .and. iband==1) &
!  &    wght2=two*occ_deavg(iband)*wtk(ikpt)/ucvol
!  ENDDEBUG

   weight=wght1+wght2
   do ipw2=1,npwdiel
!   Only fills lower half of the matrix (here, the susceptibility matrix)
!   Note that wfprod of the first index must behave like a density,
!   so that it is used as generated by fourwf, while wfprod of the
!   second index will be implicitely used to make a scalar product
!   with a potential change, meaning that its complex conjugate must be
!   used. This explains the following signs...
!   DEBUG
!   write(6, '(a,i3,2es14.6)' )&
!   &     ' ipw,wfprod',ipw2,wfprod(1,ipw2),wfprod(2,ipw2)
!   ENDDEBUG
    do ipw1=ipw2,npwdiel
     susmat(1,ipw1,isp,ipw2,isp)=susmat(1,ipw1,isp,ipw2,isp)+&
&     weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
     susmat(2,ipw1,isp,ipw2,isp)=susmat(2,ipw1,isp,ipw2,isp)+&
&     weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
    end do
   end do

   if( occopt>=3 .and. testocc==1) then
!   Accumulate product of band densities by their doccde, for the
!   computation of the effect of change of Fermi level.
    do ipw=1,npwdiel
     drhode(1,ipw,isp)=drhode(1,ipw,isp)+wfprod(1,ipw)*wght1
     drhode(2,ipw,isp)=drhode(2,ipw,isp)+wfprod(2,ipw)*wght1
    end do
!   Also accumulate weighted sum of doccde
    sumdocc=sumdocc+wght1
   end if

!  End condition of metallic occupancies or extrapolation
  end if

! End loop on iband
 end do
 call xsum_mpi(wfrspa,mpi_enreg%commcart,ier)
 if(occopt>=3)call xsum_mpi(drhode,mpi_enreg%commcart,ier)
 call xsum_mpi(susmat,mpi_enreg%commcart,ier)
 if(extrap==1)call xsum_mpi(rhoextrap,mpi_enreg%commcart,ier)
 if(occopt>=3)call xsum_mpi(sumdocc,mpi_enreg%commcart,ier)
 call timab(87,2,tsec)

 if(mpi_enreg%mode_para=='b') susmat(:,:,:,:,:)=susmat(:,:,:,:,:)/real(mpi_enreg%nproc_fft*mpi_enreg%nproc_band,dp)

!--Wavefunctions have been generated in real space--------------------------
 iproc=-1
 call timab(88,1,tsec)

!Compute product of wavefunctions for different bands
 tolocc=1.0d-3
 if(nband_k>1)then
  do ibd1=1,nband_k-1
   do ibd2=ibd1+1,nband_k
    iproc=iproc+1
    if(modulo(iproc,mpi_enreg%nproc_fft*mpi_enreg%nproc_band) /= mpi_enreg%me_cart_2d) cycle
!   If the occupation numbers are sufficiently different, or
!   if extrapolation is used and the corresponding factor is not zero,
!   then there is a contribution
    occdiff=occ(ibd1+bdtot_index)-occ(ibd2+bdtot_index)
    if( abs(occdiff)>tolocc      .or. &
&    ( extrap==1 .and.            &
&    ( abs(occ_deavg(ibd1)) + abs(occ_deavg(ibd2)) ) >tolocc ) &
&    ) then

     eigdiff=eigen(ibd1+bdtot_index)-eigen(ibd2+bdtot_index)
!    DEBUG
!    write(6,*)' susk : contribution from bands',ibd1,ibd2
!    write(6,*)'   occ diff =',occdiff
!    write(6,*)'   eig diff =',eigdiff
!    ENDDEBUG

!    Store the contribution in wfraug
     do i3=1,ndiel3
      do i2=1,ndiel2
       do i1=1,ndiel1
        wfraug(1,i1,i2,i3)=wfrspa(1,i1,i2,i3,ibd1)*wfrspa(1,i1,i2,i3,ibd2)&
&        +wfrspa(2,i1,i2,i3,ibd1)*wfrspa(2,i1,i2,i3,ibd2)
        wfraug(2,i1,i2,i3)=wfrspa(2,i1,i2,i3,ibd1)*wfrspa(1,i1,i2,i3,ibd2)&
&        -wfrspa(1,i1,i2,i3,ibd1)*wfrspa(2,i1,i2,i3,ibd2)
       end do
      end do
     end do

!    DEBUG
!    norm=zero ; normr=zero
!    do i3=1,ndiel3
!    do i2=1,ndiel2
!    do i1=1,ndiel1
!    norm=norm+wfraug(1,i1,i2,i3)**2+wfraug(2,i1,i2,i3)**2
!    normr=normr+wfraug(1,i1,i2,i3)**2
!    end do
!    end do
!    end do
!    write(6,*)' norm in real space =',norm/dble(nfftdiel)
!    write(6,*)' norm of real part  =',normr/dble(nfftdiel)
!    ENDDEBUG

!    Performs the Fourier Transform of the product, and store it in wfprod
     tim_fourwf=9
     call fourwf(1,rhoaug,dummy,wfprod,wfraug,gbound_diel,gbound_diel,&
&     istwf_k,kg_diel,kg_diel,&
&     mgfftdiel,mpi_enreg,1,ngfftdiel,1,npwdiel,ndiel4,ndiel5,ndiel6,3,paral_kgb_diel,tim_fourwf,weight)

!    In case of PAW, add compensation charge contribution
     if (usepaw==1) then
      call pawsushat(atindx1,cprj_k,gbound_diel,gylmg_diel,ibd1,ibd2,istwf_k,kg_diel,&
&      lmax_diel,mgfftdiel,mpi_enreg_diel,natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,&
&      ngfftdiel,npwdiel,nspinor,ntypat,0,paral_kgb_diel,&
&      pawang,pawtab,ph3d_diel,typat,wfprod,dummy)
     end if

!    Perform now the summation
     wght1=zero ; wght2=zero
     if(abs(occdiff)>tolocc)wght1= occdiff/eigdiff * two*wtk(ikpt)/ucvol
     if(extrap==1)then
      wght2=(occ_deavg(ibd1)+occ_deavg(ibd2)) * two*wtk(ikpt)/ucvol
     end if
     weight=wght1+wght2

!    DEBUG
!    write(6,*)' weight =',weight
!    norm=zero
!    do ipw=1,npwdiel
!    norm=norm+wfprod(1,ipw)**2+wfprod(2,ipw)**2
!    end do
!    write(6,*)' norm in reciprocal space  =',norm
!    ENDDEBUG

     do ipw2=1,npwdiel
!     Only fills lower half of the matrix (here, the susceptibility matrix)
!     Note that wfprod of the first index must behave like a density,
!     so that it is used as generated by fourwf, while wfprod of the
!     second index will be implicitely used to make a scalar product
!     with a potential change, meaning that its complex conjugate must be
!     used. This explains the following signs...
      do ipw1=ipw2,npwdiel
       susmat(1,ipw1,isp,ipw2,isp)=susmat(1,ipw1,isp,ipw2,isp)+&
&       weight*(wfprod(1,ipw1)*wfprod(1,ipw2)+wfprod(2,ipw1)*wfprod(2,ipw2))
       susmat(2,ipw1,isp,ipw2,isp)=susmat(2,ipw1,isp,ipw2,isp)+&
&       weight*(wfprod(2,ipw1)*wfprod(1,ipw2)-wfprod(1,ipw1)*wfprod(2,ipw2))
      end do
     end do

!    End condition of different occupation numbers or extrapolation
    end if
!   End internal loop over bands
   end do
!  End external loop over bands
  end do
! End condition of having more than one band
 end if
 call xsum_mpi(susmat,mpi_enreg%commcart,ier)
!DEBUG
!write(6,*)' susk : exit , write susmat'
!do ipw1=1,npwdiel
!write(6,*)ipw1,susmat(1,ipw1,1,ipw1,1),susmat(2,ipw1,1,ipw1,1)
!end do
!write(6,*)' susk : end of susmat '
!stop
!ENDDEBUG

 deallocate(cwavef,dummy,rhoaug,wfprod,wfraug,wfrspa)
 if(mpi_enreg%mode_para=='b') deallocate(band_loc,treat_band,cg_local,kg_k_gather_all)

 call timab(88,2,tsec)
end subroutine susk
!!***
