!{\src2tex{textfont=tt}}
!!****f* ABINIT/wfkfermi3
!! NAME
!! wfkfermi3
!!
!! FUNCTION
!! This routine computes the partial Fermi-level density at a given k-point,
!! and the fixed contribution to the 1st-order Fermi energy (nonlocal
!!  and kinetic)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (DRH, XG, AR, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cgq(2,mcgq)=array for planewave coefficients of wavefunctions.
!!  cplex=1 if rhoaug1 is real, 2 if rhoaug1 is complex
!!  dimekb=first dimension of ekb (see ekb_typ)
!!  dimffnlk=second dimension of ffnlk (1+number of derivatives)
!!  dimffnl1=second dimension of ffnl1 and ffnlkq (1+number of derivatives)
!!  dkinpw(npw_k)=derivative of the (modified) kinetic energy for each
!!    plane wave at k (Hartree)
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  ekb_typ(dimekb,1,nspinor**2)=
!!    ->Norm conserving : (Real) Kleinman-Bylander energies (hartree)
!!                        for the displaced atom
!!                        for number of basis functions (l,n) (lnmax)
!!                        dimekb=lnmax
!!    ->PAW : (Real, symmetric) Frozen part of Dij coefficients
!!                        to connect projectors
!!                        for the displaced atom
!!                        for number of basis functions (l,m,n) (lmnmax)
!!                        dimekb=lmnmax*(lmnmax+1)/2
!!  ffnlk(npw_k,dimffnlk,lmnmax,1)=nonloc form factors at k, for the displaced atom.
!!  ffnlkq(npw1_k,dimffnl1,lmnmax,1)=nonloc form fact at k+q for the displaced atom
!!  ffnl1(npw1_k,dimffnl1,lmnmax,ntypat)=nonloc form factors at k+q
!!  gbound(2*mgfft+8,2)=G sphere boundary
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  icg=shift to be applied on the location of data in the array cg
!!  icgq=shift to be applied on the location of data in the array cgq
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  indlmn_typ(6,lmnmax,1)=indlmn info for the displaced atom
!!  ipert=type of the perturbation
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kinpw1(npw1_k)=(modified) kinetic energy for each plane wave at k+q
!!    (Hartree)
!!  kpg_k(npw_k,nkpg)= (k+G) components at k (only if useylm=1)
!!  kpg1_k(npw1_k,nkpg1)= (k+G) components at k+q (only if useylm=1)
!!  kpt(3)=reduced coordinates of k points.
!!  lmnmax= max number of (l,n)   comp. over all type of psps
!!  matblk=dimension of the array ph3d
!!  mband=maximum number of bands
!!  mcgq=second dimension of the cgq array
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpssoang= 1+maximum (spin*angular momentum) for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  nkpg,nkpg1=second dimensions of kpg_k and kpg1_k (0 if useylm=0)
!!  nkpt=number of k points
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  n4,n5,n6 used for dimensioning real space arrays
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  ph3d(2,npw,matblk)=3-dim structure factors, for each atom and plane wave.
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pspso_typ(1)=spin-orbit info for the displaced atom
!!  rhoaug1(cplex*n4,n5,n6)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output)
!!  rocceig(nband_k,nband_k)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!   if this
!!   ratio has been attributed to the band n (second argument), 0.0_dp otherwise
!!  wfftgs=struct info for GS wf disk files.
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  eig1_k(2*nband_k**2)=first-order eigenvalues (hartree)
!!  fe1fixed_k(nband_k)=contribution to 1st-order Fermi energy
!!      from changes of occupation from all bands at this k point.
!!  fe1norm_k(nband_k)=contribution to normalization for above
!!  rhoaug1(cplex*n4,n5,n6)= Fermi-level density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output).
!!
!! TODO
!!
!! PARENTS
!!      rhofermi3
!!
!! CHILDREN
!!      eig1fixed,fourwf,leave_new,status,timab,wffreaddatarec,wffreadnpwrec
!!      wffreadskiprec,wrtout,xcomm_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wfkfermi3(cg,cgq,cplex,dimekb,dimffnlk,dimffnl1,dkinpw,dtfil,dtset,&
& eig1_k,ekb_typ,fe1fixed_k,fe1norm_k,ffnlk,ffnlkq,ffnl1,gbound,gs_hamkq,&
& icg,icgq,idir,ikpt,indlmn_typ,ipert,isppol,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,&
& lmnmax,matblk,mband,mcgq,mgfft,mkmem,mpi_enreg,&
& mpsang,mpssoang,mpw,mpw1,natom,nband_k,nkpg,nkpg1,nkpt,&
& npw_k,npw1_k,nspinor,nsppol, ntypat,n4,n5,n6,occ_k,ph3d,prtvol,&
& psps,pspso_typ,rhoaug1,rocceig,wfftgs,wtk_k)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_13io_mpi
 use interfaces_16response, except_this_one => wfkfermi3
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,dimekb,dimffnl1,dimffnlk,icg,icgq,idir,ikpt,ipert
 integer,intent(in) :: isppol,lmnmax,matblk,mband,mcgq,mgfft,mkmem,mpsang
 integer,intent(in) :: mpssoang,mpw,mpw1,n4,n5,n6,natom,nkpg,nkpg1,nkpt,npw1_k
 integer,intent(in) :: nsppol,ntypat,prtvol
 integer,intent(inout) :: nband_k,npw_k,nspinor
 real(dp),intent(in) :: wtk_k
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(in) :: gs_hamkq
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wfftgs
!arrays
 integer,intent(in) :: gbound(2*mgfft+8,2),indlmn_typ(6,lmnmax,1)
 integer,intent(in) :: kg1_k(3,npw1_k),kg_k(3,npw_k),pspso_typ(1)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),cgq(2,mcgq)
 real(dp),intent(in) :: dkinpw(npw_k),ekb_typ(dimekb,1,nspinor**2)
 real(dp),intent(in) :: ffnl1(npw1_k,dimffnl1,lmnmax,ntypat)
 real(dp),intent(in) :: ffnlk(npw_k,dimffnlk,lmnmax,1)
 real(dp),intent(in) :: ffnlkq(npw1_k,dimffnl1,lmnmax,1),kinpw1(npw1_k)
 real(dp),intent(in) :: kpg1_k(npw1_k,nkpg1),kpg_k(npw_k,nkpg),kpt(3)
 real(dp),intent(in) :: occ_k(nband_k),rocceig(nband_k,nband_k)
 real(dp),intent(inout) :: ph3d(2,npw1_k,matblk),rhoaug1(cplex*n4,n5,n6)
 real(dp),intent(out) :: eig1_k(2*nband_k**2),fe1fixed_k(nband_k)
 real(dp),intent(out) :: fe1norm_k(nband_k)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=18
 integer,save :: count=0,nskip=0
 integer :: counter,i1,i2,i3,iband,ier,ierr,iexit,ig,igs,ii,index,index_cgq
 integer :: index_eig1,inonsc,iowf,iproj,ipsang,ipw,ipw1,ispinor,istwf_k,iwavef
 integer :: mcgnpw,mcgnpw1,me,n1,n2,n3,nkpt_max,nrecwf,openexit,quit,spaceComm
 integer :: tag,test_ddk,tim_fourwf,tim_rwwf
 real(dp) :: aa,ai,ar,facti,factr,im0,im1,invocc,re0,re1,resid,residk,scprod
 real(dp) :: valuei,valuer,weight
 character(len=500) :: message
!arrays
 integer,allocatable :: kg_dum(:,:)
 real(dp) :: dummy(2,1),qphon(3),tsec(2)
 real(dp),allocatable :: cgddk(:,:),cgnew(:,:),cgnow(:,:),cgtgs(:,:)
 real(dp),allocatable :: cwave0(:,:),cwave1(:,:),cwavef_sp(:,:),eig_dum(:)
 real(dp),allocatable :: ghc(:,:),grnk(:),gvnl1(:,:),gvnlc(:,:),occ_dum(:)
 real(dp),allocatable :: rhoaug(:,:,:),wfraug(:,:,:,:),wfraug1(:,:,:,:)

! *********************************************************************

 nkpt_max=50
 if(mpi_enreg%paral_compil_kpt==1)nkpt_max=-1

 call status(0,dtfil%filstat,iexit,level,'enter         ')

!Structured debugging if prtvol==-level
 if(prtvol==-level)then
  write(message,'(80a,a,a)') ('=',ii=1,80),ch10,'wfkfermi3 : enter'
  call wrtout(06,message,'PERS')
 end if

 quit=0
!Init me
 call xme_init(mpi_enreg,me)
!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)

 n1=gs_hamkq%ngfft(1) ; n2=gs_hamkq%ngfft(2) ; n3=gs_hamkq%ngfft(3)
 qphon(1:3)=dtset%qptn(1:3)

 istwf_k=gs_hamkq%istwf_k

 allocate(ghc(2,npw1_k*nspinor),gvnlc(2,npw1_k*nspinor))
 allocate(gvnl1(2,npw1_k*nspinor))

 if(prtvol>2 .or. ikpt<=nkpt_max)then
  write(message, '(a,a,i5,2x,a,3f9.5,2x,a)' ) ch10,&
&  ' Non-SCF iterations; k pt #',ikpt,'k=',kpt(:),'band residuals:'
  call wrtout(06,message,'PERS')
 end if

 allocate(wfraug(2,n4,n5,n6),wfraug1(2,n4,n5,n6))
 allocate(rhoaug(n4,n5,n6))
 allocate(cwave0(2,npw_k*nspinor))
 allocate(cwave1(2,npw1_k*nspinor))
 wfraug1(:,:,:,:)=zero
 eig1_k(:)=zero

!Read the npw and kg records of wf files
!NOTE : it should be possible to use rwwf in the present routine
 call status(0,dtfil%filstat,iexit,level,'before WffRead')
 if(mkmem==0)then
  call WffReadNpwRec(ierr,ikpt,isppol,nband_k,npw_k,nspinor,wfftgs)
! Skip k+G and eigenvalue records in wfftgs (already in eigen0)
  call WffReadSkipRec(ierr,2,wfftgs)
 end if

!Null potentially unassigned output variables
 fe1fixed_k(:)=zero; fe1norm_k(:)=zero

 call timab(139,1,tsec)

!Loop over bands
 do iband=1,nband_k

  if(mpi_enreg%paral_compil_kpt==1)then
   if(mpi_enreg%proc_distrb(ikpt, iband,isppol) /= mpi_enreg%me) then
    if(mkmem==0)then
     call WffReadSkipRec(ierr,1,wfftgs)
    end if
    cycle
   end if
  end if ! paral

! Read ground-state wavefunctions
  if(mkmem/=0)then
!  $OMP PARALLEL DO PRIVATE(ipw) &
!  $OMP&SHARED(cg,cwave0,iband,icg,npw_k,nspinor)
   do ipw=1,npw_k*nspinor
    cwave0(1,ipw)=cg(1,ipw+(iband-1)*npw_k*nspinor+icg)
    cwave0(2,ipw)=cg(2,ipw+(iband-1)*npw_k*nspinor+icg)
   end do
!  $OMP END PARALLEL DO
  else
   call timab(288,1,tsec)
   call WffReadDataRec(cwave0,ierr,2*npw_k*nspinor,wfftgs)
   call timab(288,2,tsec)
  end if


  if(prtvol>=10)then
   call status(0,dtfil%filstat,iexit,level,'after wf read ')
  end if

  inonsc=1

  counter=100*iband+inonsc
! Because in this loop, the CPU time matters, the writing
! in the STATUS file is usually inhibited
  if(prtvol>=10)then
   call status(counter,dtfil%filstat,iexit,level,'loop iband    ')
  end if

  if ( abs(occ_k(iband)) > tol8         .and.            &
&  abs(rocceig(iband,iband)) > tol8        ) then

   if(prtvol>=10)then
    call status(counter,dtfil%filstat,iexit,level,'call eig1fixed    ')
   end if

!  Note that the following translation occurs in the called routine :
!  iband->band, nband_k->nband, npw_k->npw, npw1_k->npw1
   call eig1fixed(iband,cwave0,dimekb,dimffnlk,dimffnl1,dkinpw,eig1_k,ekb_typ,&
&   ffnlk,ffnlkq,ffnl1,dtfil%filstat,gs_hamkq,gvnl1,idir,indlmn_typ,&
&   ipert,kg_k,kg1_k,kinpw1,kpg_k,kpg1_k,kpt,lmnmax,matblk,mgfft,mpi_enreg,&
&   mpsang,mpssoang,natom,nband_k,nkpg,nkpg1,npw_k,npw1_k,nspinor,ntypat,ph3d,prtvol,&
&   pspso_typ)

!  Compute the 1st-order wavefunction component controlled by the
!  1st-order Fermi energy and the fixed contribution to the 1st-order
!  Fermi energy from this k point
   invocc=1.0_dp/occ_k(iband)
   index_eig1=2*iband-1+(iband-1)*2*nband_k
   index_cgq=npw1_k*nspinor*(iband-1)+icgq
   factr= rocceig(iband,iband)*invocc
!  $OMP PARALLEL DO PRIVATE(ii) &
!  $OMP&SHARED(cgq,cwave1,facti,factr,index_cgq,npw1_k,nspinor)
   do ii=1,npw1_k*nspinor
    cwave1(1,ii)= factr*cgq(1,ii+index_cgq)
    cwave1(2,ii)= factr*cgq(2,ii+index_cgq)
   end do
!  $OMP END PARALLEL DO
   fe1fixed_k(iband)=2.0_dp*factr*eig1_k(index_eig1)
   fe1norm_k(iband)=2.0_dp*factr

   do ispinor=1,nspinor

    if(prtvol>=10)then
     call status(counter,dtfil%filstat,iexit,level,'density update')
    end if

!   Compute contribution to density

!   The factor 2 is not the spin factor (see Eq.44 of PRB55,10337 (1997))
    weight=2.0_dp*occ_k(iband)*wtk_k/gs_hamkq%ucvol

    tim_fourwf=5
    if(ispinor==1)then
     call fourwf(cplex,rhoaug1,cwave1,dummy,wfraug1,&
&     gs_hamkq%gbound,gs_hamkq%gbound,&
&     istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     npw1_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
    else
     allocate(cwavef_sp(2,npw1_k))
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(cwave1,cwavef_sp,npw1_k)
     do ipw=1,npw1_k
      cwavef_sp(1,ipw)=cwave1(1,ipw+npw1_k)
      cwavef_sp(2,ipw)=cwave1(2,ipw+npw1_k)
     end do
!    $OMP END PARALLEL DO
     call fourwf(cplex,rhoaug1,cwavef_sp,dummy,wfraug1,&
&     gs_hamkq%gbound,gs_hamkq%gbound,&
&     istwf_k,kg1_k,kg1_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     npw1_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
     deallocate(cwavef_sp)
    end if

    tim_fourwf=5
    if(ispinor==1)then
     call fourwf(1,rhoaug,cwave0,dummy,wfraug,gbound,gbound,&
&     istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     npw_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
    else
     allocate(cwavef_sp(2,npw_k))
!    $OMP PARALLEL DO PRIVATE(ipw) &
!    $OMP&SHARED(cwave0,cwavef_sp,npw_k)
     do ipw=1,npw_k
      cwavef_sp(1,ipw)=cwave0(1,ipw+npw_k)
      cwavef_sp(2,ipw)=cwave0(2,ipw+npw_k)
     end do
!    $OMP END PARALLEL DO
     call fourwf(1,rhoaug,cwavef_sp,dummy,wfraug,gbound,gbound,&
&     istwf_k,kg_k,kg_k,mgfft,mpi_enreg,1,gs_hamkq%ngfft,&
&     npw_k,1,n4,n5,n6,0,dtset%paral_kgb,tim_fourwf,weight)
     deallocate(cwavef_sp)
    end if

!   Accumulate density
!   $OMP PARALLEL DO PRIVATE(i1,i2,i3) &
!   $OMP&SHARED(n1,n2,n3,rhoaug1,weight,wfraug,wfraug1)
!   OCL SCALAR
    do i3=1,n3
     do i2=1,n2
      do i1=1,n1
       rhoaug1(i1,i2,i3)=rhoaug1(i1,i2,i3)+&
&       weight*( wfraug(1,i1,i2,i3)*wfraug1(1,i1,i2,i3) &
&       +wfraug(2,i1,i2,i3)*wfraug1(2,i1,i2,i3)  )
      end do
     end do
    end do
!   $OMP END PARALLEL DO

   end do ! ispinor=1,nspinor

!  End of non-zero occupation and rocceig
  end if

! End loop over bands
 end do

 call timab(139,2,tsec)
 call timab(130,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'after loops   ')

 deallocate(rhoaug,wfraug,wfraug1)
 deallocate(cwave0,cwave1)


 deallocate(ghc,gvnlc,gvnl1)

 call status(0,dtfil%filstat,iexit,level,'deallocate    ')

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
  write(message,'(a1,a,a1,a,i2,a)') ch10,&
&  ' fermie3 : exit ',&
&  ch10,'  prtvol=-',level,', debugging mode => stop '
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

 call status(0,dtfil%filstat,iexit,level,'exit          ')

 call timab(130,2,tsec)

!DEBUG
!write(6,*)' wfkfermi3 : exit '
!call flush(6)
!if(count==26)stop
!stop
!ENDDEBUG

end subroutine wfkfermi3
!!***
