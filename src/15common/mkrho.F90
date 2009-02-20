!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkrho
!! NAME
!! mkrho
!!
!! FUNCTION
!! Compute charge density rho(r) and rho(G) in electrons/bohr**3
!! from input wavefunctions, band occupations, and k point wts.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, LSI, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wf in G space
!!  densymop_gs <type(dens_sym_operator_type)>=the density symmetrization
!!   operator (ground-state symmetries)
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | istwfk(nkpt)=input option parameter that describes the storage of wfs
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem=maximum number of k points in core memory
!!   | mpw=maximum allowed value for npw
!!   | nband(nkpt*nsppol)=number of bands to be included in summation
!!   |  at each k point for each spin channel
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | ngfft(18)=contain all needed information about 3D FFT, 
!!   |  see ~abinit/doc/input_variables/vargs.htm#ngfft
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetry elements in group (at least 1 for identity)
!!   | symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!   | wtk(nkpt)=k point weights (they sum to 1.0)
!!  irrzon(nfft**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  kg(3,mpw*mkmem)=reduced planewave coordinates
!!  mpi_enreg=informations about MPI parallelization
!!  npwarr(nkpt)=number of planewaves and boundary planewaves at each k
!!  nspinor=number of spinorial components of the wavefunctions
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  phnons(2,nfft**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  tim_mkrho=timing code of the calling routine(can be set to 0 if not attributed)
!!  ucvol=unit cell volume (Bohr**3)
!!  unkg=unit number for (k+G) sphere data file
!!  wffnow=struct info for current wf disk file
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!
!! OUTPUT
!! rhog(2,nfft)=total electron density in G space
!! rhor(nfft,nspden)=electron density in r space
!!   (if spin polarized, array contains total density in first half and
!!    spin-up density in second half)
!!   (for non-collinear magnetism, first element: total density, 3 next ones: mx,my,mz
!!    in units of hbar/2)
!!
!! PARENTS
!!      energy,gstate,respfn,vtorho
!!
!! CHILDREN
!!      fftpac,fourwf,hdr_skip,leave_test,prtrhomxmn,rdnpw,rwwf,sphereboundary
!!      symrhg,timab,wrtout,xcomm_init,xdefineoff,xmaster_init,xme_init
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkrho(cg,densymop_gs,dtset,irrzon,kg,mpi_enreg,&
     & npwarr,nspinor,occ,phnons,rhog,rhor,tim_mkrho,ucvol,unkg,wffnow,wfs)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_13io_mpi
 use interfaces_14iowfdenpot
 use interfaces_14wfs
 use interfaces_15common, except_this_one => mkrho
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: tim_mkrho,unkg
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(dens_sym_operator_type),intent(in) :: densymop_gs
 type(wffile_type),intent(inout) :: wffnow
 type(wvl_wf_type),intent(in) :: wfs
!no_abirules
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 integer, intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
 integer, intent(in) :: kg(3,dtset%mpw*dtset%mkmem),npwarr(dtset%nkpt)
 real(dp), intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp), intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
!nfft**(1-1/nsym) is 1 if nsym==1, and nfft otherwise
 real(dp), intent(in) :: phnons(2,(dtset%ngfft(1)*dtset%ngfft(2)*dtset%ngfft(3))**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
 real(dp), intent(out) :: rhor(dtset%nfft,dtset%nspden),rhog(2,dtset%nfft)

!Local variables-------------------------------
!local variables for mpialltoallv
!scalars
 integer,save :: nskip=0
 integer :: bdtot_index,blocksize,formeig,iband,iblock,icg,ierr,ifft,ikg,ikpt
 integer :: iproc,ipw,ispden,isppol,istwf_k,master,mcg_disk,me,n1,n2,n3,n4,n5
 integer :: n6,nband_k,nbdblock,ndatarecv,nfftot,npw_k,npw_tot,oldspacecomm
 integer :: spaceComm,tim_fourwf,tim_rwwf
 real(dp) :: weight
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk_local
!arrays
 integer,allocatable :: gbound(:,:),kg_dum(:,:),kg_k(:,:),kg_k_gather(:,:)
 integer,allocatable :: kg_k_gather_all(:,:),npw_per_proc(:),rdispls(:)
 integer,allocatable :: rdispls_all(:),rdisplsloc(:),recvcounts(:)
 integer,allocatable :: recvcountsloc(:),sdispls(:),sdisplsloc(:),sendcounts(:)
 integer,allocatable :: sendcountsloc(:)
 real(dp) :: dummy(2,1),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),cwavef(:,:),cwavef1(:,:),cwavef_x(:,:)
 real(dp),allocatable :: cwavef_y(:,:),eig_dum(:),occ_dum(:),occ_k(:)
 real(dp),allocatable :: rhoaug(:,:,:),rhoaug_down(:,:,:),rhoaug_mx(:,:,:)
 real(dp),allocatable :: rhoaug_my(:,:,:),rhoaug_up(:,:,:),wfraug(:,:,:,:)

! *************************************************************************

!DEBUG
!write(6,*)' mkrho : enter '
!if(.true.)stop
!ENDDEBUG

 call timab(290+tim_mkrho,1,tsec)
 call timab(299,1,tsec)

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)

!PATCH mkrho // KPT & FFT me-->me_kpt
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
& (mpi_enreg%paral_compil_fft==1)) then
  me = mpi_enreg%me_kpt
 end if

!Init master
 call xmaster_init(mpi_enreg,master)

!zero the charge density array in real space
 do ispden=1,dtset%nspden
! $OMP PARALLEL DO PRIVATE(ifft) &
! $OMP&SHARED(dtset%nfft,rhor)
  do ifft=1,dtset%nfft
   rhor(ifft,ispden)=0.0_dp
  end do
! $OMP END PARALLEL DO
 end do

!WVL - Branching with a separate mkrho procedure
!in wavelet.
 if (dtset%usewvl == 1) then
  call wvl_mkrho(dtset, mpi_enreg, occ, rhor, wfs)
  return
 end if
!WVL - Following is done in plane waves.


!dtset%mkmem==0 means wf and kg info on disk file
 if (dtset%mkmem==0) then

! Skip header of wffnow
  call hdr_skip(wffnow,ierr)

! Define offsets, in case of MPI I/O
  formeig=0
  call xdefineOff(formeig,wffnow,mpi_enreg,dtset%nband,npwarr,nspinor,dtset%nsppol,dtset%nkpt)

  mcg_disk=dtset%mpw*nspinor*dtset%mband
  allocate(cg_disk(2,mcg_disk))

 end if

!start loop over spin and k points
 bdtot_index=0
 icg=0

 n1 = dtset%ngfft(1) ; n2 = dtset%ngfft(2) ; n3 = dtset%ngfft(3)
!n4,n5,n6 are FFT dimensions, modified to avoir cache trashing
 n4 = dtset%ngfft(4) ; n5 = dtset%ngfft(5) ; n6 = dtset%ngfft(6)
 allocate(cwavef(2,dtset%mpw*nspinor),rhoaug(n4,n5,n6),wfraug(2,n4,n5,n6))
 if(dtset%nspden==4) then
  allocate(rhoaug_up(n4,n5,n6),rhoaug_down(n4,n5,n6))
  allocate(rhoaug_mx(n4,n5,n6),rhoaug_my(n4,n5,n6))
  rhoaug_up(:,:,:)=zero
  rhoaug_down(:,:,:)=zero
  rhoaug_mx(:,:,:)=zero
  rhoaug_my(:,:,:)=zero
 end if
 do isppol=1,dtset%nsppol

! Rewind the kpgsph data file on unit unkg
  if (dtset%mkmem==0) rewind (unkg)
  ikg=0

  rhoaug(:,:,:)=0.0_dp
  do ikpt=1,dtset%nkpt

   nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
   npw_k=npwarr(ikpt)
   istwf_k = dtset%istwfk(ikpt)

   if(mpi_enreg%paral_compil_kpt==1)then
    if(mpi_enreg%parareel == 0) then
!    BEGIN TF_CHANGES
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&     -me))/=0) then
!     END TF_CHANGES
      bdtot_index=bdtot_index+nband_k
      cycle
     end if
    else
     if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt) &
&     /= mpi_enreg%me) then
      bdtot_index=bdtot_index+nband_k
      cycle
     end if
    end if
   end if

   allocate(gbound(2*dtset%mgfft+8,2),kg_k(3,npw_k))

!  Do i/o as needed
   if (dtset%mkmem==0) then

    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,unkg)

!   Read k+g data
    read (unkg) kg_k(1:3,1:npw_k)

    call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)

!   Read the wavefunction block for ikpt,isppol
    if((mpi_enreg%paralbd==0) .or. (mpi_enreg%paralbd>1)) tim_rwwf=5
    if(mpi_enreg%paralbd==1)tim_rwwf=12
    allocate(eig_dum(dtset%mband),kg_dum(3,0),occ_dum(dtset%mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,dtset%mband,mcg_disk,&
&    mpi_enreg,nband_k,nband_k,npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
    deallocate(eig_dum,kg_dum,occ_dum)

   else

    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)

!   FB Transpose the kg_k array. TO BE CLEANED
    if (mpi_enreg%mode_para=='b') then 
     spaceComm=mpi_enreg%comm_band
!    blocksize=mpi_enreg%nproc_band

     allocate(sdispls       (mpi_enreg%nproc_band))
     allocate(sdisplsloc    (mpi_enreg%nproc_band))
     allocate(sendcounts    (mpi_enreg%nproc_band))
     allocate(sendcountsloc (mpi_enreg%nproc_band))
     allocate(rdispls       (mpi_enreg%nproc_band))
     allocate(rdisplsloc    (mpi_enreg%nproc_band))
     allocate(recvcounts    (mpi_enreg%nproc_band))
     allocate(recvcountsloc (mpi_enreg%nproc_band))
     call xallgather_mpi(npw_k,recvcounts,spaceComm,ierr)
     rdispls(1)=0
     do iproc=2,mpi_enreg%nproc_band
      rdispls(iproc)=rdispls(iproc-1)+recvcounts(iproc-1)
     end do
     ndatarecv=rdispls(mpi_enreg%nproc_band)+recvcounts(mpi_enreg%nproc_band)

     allocate(kg_k_gather(3,ndatarecv))
     recvcountsloc(:)=recvcounts(:)*3
     rdisplsloc(:)=rdispls(:)*3
     call xallgatherv_mpi(kg_k,3*npw_k,kg_k_gather,recvcountsloc(:),rdisplsloc,spaceComm,ierr)

     oldspacecomm=mpi_enreg%comm_fft
     allocate(npw_per_proc(mpi_enreg%nproc_fft),rdispls_all(mpi_enreg%nproc_fft))
     call xallgather_mpi(ndatarecv,npw_per_proc,oldspacecomm,ierr)
     rdispls_all(1)=0
     do iproc=2,mpi_enreg%nproc_fft
      rdispls_all(iproc)=rdispls_all(iproc-1)+npw_per_proc(iproc-1)
     end do
     npw_tot=rdispls_all(mpi_enreg%nproc_fft)+npw_per_proc(mpi_enreg%nproc_fft)
     allocate(kg_k_gather_all(3,npw_tot))
     call xallgatherv_mpi&
&     (kg_k_gather,3*ndatarecv,kg_k_gather_all,3*npw_per_proc(:),3*rdispls_all,oldspacecomm,ierr)
     call sphereboundary(gbound,istwf_k,kg_k_gather_all,dtset%mgfft,npw_tot)
     deallocate(kg_k_gather_all,npw_per_proc,rdispls_all)
     deallocate(sendcounts,recvcounts,sdispls,rdispls)
     deallocate(sendcountsloc,sdisplsloc)
     deallocate(recvcountsloc,rdisplsloc)
    else
     call sphereboundary(gbound,istwf_k,kg_k,dtset%mgfft,npw_k)
    end if
   end if ! dtset%mkmem==0

!  Loop over bands to fft and square for rho(r)
!  Shoulb be changed to treat bands by batch always

!  DEBUG
!  write(6,*)' mkrho : mpi_enreg%mode_para=',mpi_enreg%mode_para
!  ENDDEBUG

   if(mpi_enreg%mode_para /= 'b') then
    do iband=1,nband_k

     if(mpi_enreg%paral_compil_kpt==1)then
      if(mpi_enreg%paralbd>=1)then
!      BEGIN TF_CHANGES
       if(mpi_enreg%proc_distrb(ikpt, iband, isppol) /= me) then
!       END TF_CHANGES
        cycle
       end if
      end if
     end if

!    Only treat occupied states
     if (abs(occ(iband+bdtot_index))>tol8) then
!     Obtain Fourier transform in fft box and accumulate the density
      if(dtset%mkmem/=0)then
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(cg,cwavef,iband,icg,npw_k,nspinor)
       do ipw=1,npw_k*nspinor
        cwavef(1,ipw)=cg(1,ipw+(iband-1)*npw_k*nspinor+icg)
        cwavef(2,ipw)=cg(2,ipw+(iband-1)*npw_k*nspinor+icg)
       end do
!      $OMP END PARALLEL DO
      else
!      $OMP PARALLEL DO PRIVATE(ipw) &
!      $OMP&SHARED(cg_disk,cwavef,iband,npw_k,nspinor)
       do ipw=1,npw_k*nspinor
        cwavef(1,ipw)=cg_disk(1,ipw+(iband-1)*npw_k*nspinor)
        cwavef(2,ipw)=cg_disk(2,ipw+(iband-1)*npw_k*nspinor)
       end do
!      $OMP END PARALLEL DO
      end if
      weight=occ(iband+bdtot_index)*dtset%wtk(ikpt)/ucvol
      if((mpi_enreg%paralbd==0) .or. (mpi_enreg%paralbd>1)) tim_fourwf=3
      if(mpi_enreg%paralbd==1)tim_fourwf=6

!     The same section of code is also found in vtowfk.F90 : should be rationalized !
      call fourwf(1,rhoaug,cwavef,dummy,wfraug,gbound,gbound,&
&      istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)

      call xcomm_init(mpi_enreg,spaceComm)

      if(nspinor==2)then

       allocate(cwavef1(2,npw_k))
!      This should be parallelized
       cwavef1(:,:)=cwavef(:,1+npw_k:2*npw_k)

!      DEBUG GZ !To obtain a x-directed magnetization(test)
!      cwavef1(1,1:npw_k)=-cwavef(2,1:npw_k)
!      cwavef1(2,1:npw_k)= cwavef(1,1:npw_k)
!      ENDDEBUG

       if(dtset%nspden==1) then

!       We need only the total density : accumulation continues on top of rhoaug
        call fourwf(1,rhoaug,cwavef1,dummy,wfraug,gbound,gbound,&
&        istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)

       else if(dtset%nspden==4) then

!       Build the four components of rho. We use only norm quantities and, so fourwf.
!       $\sum_{n} f_n \Psi^{* \alpha}_n \Psi^{\alpha}_n =\rho^{\alpha \alpha}$
!       $\sum_{n} f_n (\Psi^{1}+\Psi^{2})^*_n (\Psi^{1}+\Psi^{2})_n=rho+m_x$
!       $\sum_{n} f_n (\Psi^{1}-i \Psi^{2})^*_n (\Psi^{1}-i \Psi^{2})_n=rho+m_y$
        allocate(cwavef_x(2,npw_k),cwavef_y(2,npw_k))
!       $(\Psi^{1}+\Psi^{2})$
        cwavef_x(:,:)=cwavef(:,1:npw_k)+cwavef1(:,1:npw_k)
!       $(\Psi^{1}-i \Psi^{2})$
        cwavef_y(1,:)=cwavef(1,1:npw_k)+cwavef1(2,1:npw_k)
        cwavef_y(2,:)=cwavef(2,1:npw_k)-cwavef1(1,1:npw_k)
        rhoaug_up(:,:,:)=rhoaug(:,:,:) !Already computed
        call fourwf(1,rhoaug_down,cwavef1,dummy,wfraug,gbound,gbound,&
&        istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)
        call fourwf(1,rhoaug_mx,cwavef_x,dummy,wfraug,gbound,gbound,&
&        istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)
        call fourwf(1,rhoaug_my,cwavef_y,dummy,wfraug,gbound,gbound,&
&        istwf_k,kg_k,kg_k,dtset%mgfft,mpi_enreg,1,dtset%ngfft,npw_k,1,n4,n5,n6,1,dtset%paral_kgb,tim_fourwf,weight)

        deallocate(cwavef_x,cwavef_y)

       end if ! dtset%nspden/=4

       deallocate(cwavef1)

      end if
!     DEBUG
!     write(6,*)' ikpt, iband, rhoaug',ikpt,iband,rhoaug(1,1,1)
!     ENDDEBUG

     else
!     Accumulate the number of one-way 3D ffts skipped
      nskip=nskip+1
     end if ! abs(occ(iband+bdtot_index))>tol8
!    End loop on iband
    end do ! iband=1,nband_k
   else !mode_para==b

!   PATCH mkrho // KPT & FFT test me / proc_distrib
    if ((mpi_enreg%paral_compil_kpt==1) .and. &
&    (mpi_enreg%paral_compil_fft==1)) then
     if(mpi_enreg%parareel == 0) then
      if(minval(abs(mpi_enreg%proc_distrb(ikpt,:,isppol) &
&      -me))/=0) then
       cycle
      end if
     end if
    end if
!   FB Cleaning of the band_fft parallelization part   
    allocate(gs_hamk_local%gbound(2*dtset%mgfft+8,2))
    gs_hamk_local%gbound(:,:)=gbound(:,:)
    gs_hamk_local%ngfft(:)=dtset%ngfft(:)
    gs_hamk_local%ucvol=ucvol

    nbdblock=nband_k/(mpi_enreg%nproc_band * mpi_enreg%bandpp)
    blocksize=nband_k/nbdblock
    if(allocated(cwavef)) deallocate(cwavef)
    allocate(cwavef(2,npw_k*nspinor*blocksize))
    allocate(occ_k(nband_k))
    occ_k(:)=occ(bdtot_index+1:bdtot_index+nband_k)
    do iblock=1,nbdblock
     cwavef(:,1:npw_k*nspinor*blocksize)=&
&     cg(:,1+(iblock-1)*npw_k*nspinor*blocksize+icg:iblock*npw_k*nspinor*blocksize+icg)
     call timab(538,1,tsec)
     call prep_fourwf(rhoaug,blocksize,cwavef,wfraug,&
&     gs_hamk_local,istwf_k,iblock,1,kg_k_gather,dtset%mgfft,mpi_enreg,nbdblock,&
&     nband_k,ndatarecv,npw_k,n4,n5,n6,occ_k,dtset%paral_kgb,dtset%wtk(ikpt))
     call timab(538,2,tsec)
    end do
    deallocate(kg_k_gather)
    deallocate(occ_k,gs_hamk_local%gbound)
   end if

   deallocate(gbound,kg_k)

   bdtot_index=bdtot_index+nband_k

   if (dtset%mkmem/=0) then
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   end if

!  End loop on ikpt:
  end do

  if(mpi_enreg%mode_para == 'b') then
   spaceComm=mpi_enreg%comm_band !Sum the contributions of the bands
   call xsum_mpi(rhoaug,spaceComm,ierr)
   spaceComm=mpi_enreg%comm_fft
   call xsum_mpi(rhoaug,spaceComm,ierr)
  end if

! Write the number of one-way 3D ffts skipped until now
  if(mpi_enreg%paral_compil_kpt==0)then
   write(message, '(a,i8)' )&
&   ' mkrho : number of one-way 3D ffts skipped in mkrho until now =',nskip
   call wrtout(06,message,'PERS')
  end if

! DEBUG
! write(6,*)' rhoaug ',rhoaug(1,1,1)
! ENDDEBUG

! Transfer density on augmented fft grid to normal fft grid in real space
! Take also into account the spin, to place it correctly in rhor.
  if(dtset%nspden==1 .or. dtset%nspden==2) then
   call fftpac(isppol,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug,1)
  else if(dtset%nspden==4) then
   ispden=1
   call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_up,1)
   ispden=2
   call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_mx,1)
   ispden=3
   call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_my,1)
   ispden=4
   call fftpac(ispden,dtset%nspden,n1,n2,n3,n4,n5,n6,dtset%ngfft,rhor,rhoaug_down,1)
   deallocate(rhoaug_up,rhoaug_down,rhoaug_mx,rhoaug_my)
  end if
! DEBUG
! write(6,*)'mkrho  rhor ',rhor(1:20,:)
! ENDDEBUG

 end do !  isppol=1,dtset%nsppol

 if(mpi_enreg%paral_compil_kpt==1)then
  call timab(63,1,tsec)
  if (mpi_enreg%parareel == 0) then
!  BEGIN TF_CHANGES
   call leave_test(mpi_enreg)
!  END TF_CHANGES
  end if
  write(message,*) 'mkrho: loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
  call timab(63,2,tsec)
 end if

 deallocate(cwavef,rhoaug,wfraug)
 if(dtset%mkmem==0)deallocate(cg_disk)

 if(mpi_enreg%paral_compil_kpt==1)then
! Recreate full rhor on all proc.
  call timab(48,1,tsec)
  call timab(71,1,tsec)

! PATCH mkrho // KPT & FFT spacecomm --> comm_kpt
  if ((mpi_enreg%paral_compil_kpt==1) .and. &
&  (mpi_enreg%paral_compil_fft==1)) then
   spaceComm=mpi_enreg%comm_kpt
  end if

  call xsum_mpi(rhor,spaceComm,ierr)
  call timab(71,2,tsec)
  call timab(48,2,tsec)
 end if

!DEBUG
!write(6,*) 'mkrho : dtset%nfft,dtset%nsppol,dtset%nsym',dtset%nfft,dtset%nsppol,dtset%nsym
!write(6,*) 'ngfft',ngfft
!write(6,*) ' ir irrzon phnons '
!do ipw=1,dtset%nfft,31
!write(6, '(i5,2i5,2es16.8)' )ipw,irrzon(ipw,:,1),phnons(:,ipw,1)
!end do
!write(6,*)' mkrho : density before symrhg'
!do ipw=1,dtset%nfft,31
!write(6, '(i5,es16.6)' )ipw,rhor(ipw,1)
!end do
!ENDDEBUG

 call timab(299,2,tsec)

 call timab(549,1,tsec)
 nfftot=dtset%ngfft(1) * dtset%ngfft(2) * dtset%ngfft(3)
 call symrhg(1,densymop_gs,irrzon,mpi_enreg,dtset%nfft,nfftot,dtset%ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
 dtset%paral_kgb,phnons,rhog,rhor,dtset%symafm)
 call timab(549,2,tsec)

 call timab(299,1,tsec)

!DEBUG
!write(6,*)' mkrho : density after symrhg'
!do ipw=1,dtset%nfft,31
!write(6, '(i5,es16.6)' )ipw,rhor(ipw,1)
!end do
!ENDDEBUG

!We now have both rho(r) and rho(G), symmetrized, and if dtset%nsppol=2
!we also have the spin-up density, symmetrized, in rhor(:,2).
!In case of non collinear magnetism, we have rho,mx,my,mz. No symmetry is applied

!Debugging output
!write(*,*)' Debugging from mkrho: rhog values'
!do ipw=1,dtset%nfft
!if (abs(rhog(1,ipw))>1.d-09.or.abs(rhog(2,ipw))>1.d-09)
!& then
!write(*,2000) ipw,rhog(1,ipw),rhog(2,ipw)
!end if
!2000  format(i10,1p,2e15.5)
!end do

!DEBUG
!write(6,*)' rhor after sym',rhor(1,:)
!write(6,*)'nsym',nsym
!ENDDEBUG

!Find and print minimum and maximum total electron density and locations
 call prtrhomxmn(6,mpi_enreg,dtset%nfft,dtset%ngfft,dtset%nspden,1,rhor)

!DEBUG
!write(6,*)' mkrho : exit '
!if(.true.)stop
!ENDDEBUG

 call timab(299,2,tsec)
 call timab(290+tim_mkrho,2,tsec)

end subroutine mkrho
!!***
