!{\src2tex{textfont=tt}}
!!****f* ABINIT/forstrnps
!! NAME
!! forstrnps
!!
!! FUNCTION
!! Compute nonlocal pseudopotential energy contribution to forces and/or stress tensor
!! as well as kinetic energy contribution to stress tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, AF, AR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wavefunctions (may be read from disk file)
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass=effective mass for electrons (1. in common case)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced coordinates (integers) of G vecs in basis
!!  kpt(3,nkpt)=k points in reduced coordinates
!!  mband=maximum number of bands
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=maximum number of k points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw= maximum number of plane waves
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)= # atoms of each type.
!!  nband(nkpt)=number of bands at each k point
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points in Brillouin zone
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves in basis and boundary at each k
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of elements in symmetry group
!!  ntypat=number of types of atoms
!!  occ(mband*nkpt*nsppol)=occupation numbers for each band over all k points
!!  occopt==option for occupancies
!!  optfor=1 if computation of forces is required
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  stress_needed=1 if computation of stress tensor is required
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries in reciprocal space (dimensionless)
!!  typat(natom)=type integer for each atom in cell
!!  unkg=unit number for (k+G) data (if used)
!!  unylm=unit number for Ylm(k) data (if used)
!!  wffnow=unit number of disk file for wf if used
!!  wtk(nkpt)=weight associated with each k point
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  if (optfor==1)
!!   grnl(3*natom)=stores grads of nonlocal energy wrt atomic coordinates
!!  if (stress_needed==1)
!!   kinstr(6)=kinetic energy part of stress tensor (hartree/bohr^3)
!!   Store 6 unique components of symmetric 3x3 tensor in the order
!!   11, 22, 33, 32, 31, 21
!!   npsstr(6)=nonlocal pseudopotential energy part of stress tensor
!!    (hartree/bohr^3)
!!
!! PARENTS
!!      forstr
!!
!! CHILDREN
!!      hdr_skip,leave_test,meanvalue_g,metric,mkffnl,mkkpg,nonlop
!!      ph1d3d,rdnpw,rwwf,stresssym,timab,wrtout,xcomm_init,xdefineoff
!!      xmaster_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine forstrnps(atindx1,cg,ecut,ecutsm,effmass,eigen,&
&  grnl,indsym,istwfk,kg,kinstr,npsstr,kpt,mband,mgfft,mkmem,mpi_enreg,mpsang,&
&  mpw,natom,nattyp,nband,nfft,ngfft,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,nsym,&
&  ntypat,occ,occopt,optfor,paw_ij,pawang,pawprtvol,pawtab,ph1d,psps,rprimd,&
&  stress_needed,symafm,symrec,typat,unkg,unylm,wffnow,wtk,xred,ylm,ylmgr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_12spacepar
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_14iowfdenpot
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mgfft,mkmem,mpsang,mpw,natom,nfft,nkpt,nspden
 integer,intent(in) :: nsppol,nsym,ntypat,occopt,optfor,pawprtvol,stress_needed
 integer,intent(in) :: unkg,unylm
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ecut,ecutsm,effmass
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: atindx1(natom),indsym(4,nsym,natom),istwfk(nkpt)
 integer,intent(in) :: kg(3,mpw*mkmem),nattyp(ntypat),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(5),npwarr(nkpt),symafm(nsym)
 integer,intent(in) :: symrec(3,3,nsym),typat(ntypat)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kpt(3,nkpt)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),wtk(nkpt),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(out) :: grnl(3*natom),kinstr(6),npsstr(6)
 type(paw_ij_type),intent(in) :: paw_ij(natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: bdtot_index,choice,cpopt,dimdij,dimekb1,dimekb2,dimffnl,fform
 integer :: formeig,ia,iatom,iband,icg,ider,idir,ierr,igs,ii,ikg,ikpt,ilm,ilmn
 integer :: indx,iproc,ipw,ishift,isp,ispden,ispinor,isppol,istwf_k,itypat
 integer :: master,matblk,mcg,mcg_disk,me,me_distrb,mu,n1,n2,n3,nband_k,nkpg
 integer :: nnlout,npw_k,paw_opt,signs,spaceComm,spacecomm_old=0,tim_nonlop
 integer :: tim_rwwf
 real(dp) :: ai,ar,arg,const,dfsm,ecutsm_inv,eig_k,fact_kin,fsm,htpisq,kgc1
 real(dp) :: kgc2,kgc3,kgk1,kgk2,kgk3,kin,term,ucvol,weight,xx
 character(len=500) :: message
!arrays
 integer,allocatable :: kg_dum(:,:),kg_k(:,:)
 real(dp) :: frac(3,3),gmet(3,3),gprimd(3,3),kinstr_priv(6),kpoint(3),kpt1(3)
 real(dp) :: nonlop_dum(1,1),rmet(3,3),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),cwavef(:,:),eig_dum(:),ekb(:,:,:)
 real(dp),allocatable :: enlout(:),ffnl(:,:,:,:),kpg_k(:,:),kstr1(:),kstr2(:)
 real(dp),allocatable :: kstr3(:),kstr4(:),kstr5(:),kstr6(:),occ_dum(:)
 real(dp),allocatable :: ph3d(:,:,:),phkxred(:,:),sij(:,:),ylm_k(:,:)
 real(dp),allocatable :: ylmgr_k(:,:,:)
 type(cprj_type) :: cprj_dum(1,1)

!*************************************************************************

!DEBUG
!write(6,*)' forstrnps : enter '
!if(.true.)stop
!ENDDEBUG

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
!Init me
 call xme_init(mpi_enreg,me)

!PATCH forstrnps // KPT & FFT me-->me_kpt
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
& (mpi_enreg%paral_compil_fft==1)) then
  me_distrb = mpi_enreg%me_kpt
  me        = mpi_enreg%me_kpt
 else
  me_distrb = mpi_enreg%me
 end if

!Init master
 call xmaster_init(mpi_enreg,master)

!Some constants
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
!Smearing of plane wave kinetic energy
 ecutsm_inv=zero;if( ecutsm>1.0d-20) ecutsm_inv=1/ecutsm
!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

!Arrays initializations
 allocate(cwavef(2,mpw*nspinor),kg_k(3,mpw),phkxred(2,natom))
 if (optfor==1) grnl(:)=zero
 if (stress_needed==1) then
  kinstr(:)=zero;npsstr(:)=zero
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Determine whether wf are being read from disk file or not
 if (mkmem==0) then
! Skip wavefunction file header
  call hdr_skip(wffnow,ierr)
  mcg_disk=mpw*nspinor*mband
  allocate(cg_disk(2,mcg_disk))
! Define offsets, in case of MPI I/O
  formeig=0
  call xdefineOff(formeig,wffnow,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)
 end if

!Common data for "nonlop" routine
 signs=1 ; idir=0  ; ishift=0 ; tim_nonlop=4
 choice=2*optfor;if (stress_needed==1) choice=10*choice+3
 if (optfor==1.and.stress_needed==1)  ishift=6
 nnlout=max(1,6*stress_needed+3*natom*optfor)
 allocate(enlout(nnlout))
 if (psps%usepaw==0) then
  paw_opt=0 ; cpopt=-1
 else
  paw_opt=2 ; cpopt=-1
 end if

!Non-local factors:
!Norm-conserving: kleimann-Bylander energies
!PAW: Dij coefficients and overlap coefficients
 if (psps%usepaw==0) then
  dimekb1=psps%dimekb;dimekb2=ntypat
  allocate(ekb(psps%dimekb,ntypat,nspinor**2))
  ekb(:,:,1)=psps%ekb(:,:)
  if (nspinor==2) then
   ekb(:,:,2)=psps%ekb(:,:)
   ekb(:,:,3:4)=zero
  end if
 else
  dimekb1=psps%dimekb*paw_ij(1)%cplex_dij;dimekb2=natom
  allocate(ekb(dimekb1,dimekb2,nspinor**2),sij(dimekb1,ntypat))
  do itypat=1,ntypat
   if (paw_ij(1)%cplex_dij==1) then
    sij(1:pawtab(itypat)%lmn2_size,itypat)=pawtab(itypat)%sij(:)
   else
    do ilmn=1,pawtab(itypat)%lmn2_size
     sij(2*ilmn-1,itypat)=pawtab(itypat)%sij(ilmn)
     sij(2*ilmn  ,itypat)=zero
    end do
   end if
  end do
 end if

!DEBUG ! Do not remove this line : needed for the gfortran compiler ?!
 write(6,*)' forstrnps : usepaw=',psps%usepaw
!ENDDEBUG

!LOOP OVER SPINS
 bdtot_index=0;icg=0
 do isppol=1,nsppol

  if (nsppol==2) then
   write(message, '(a,i3)' )' ****  In forstrnps for isppol=',isppol
   call wrtout(06,message,'COLL')
  end if

! Rewind temporary disk files
  if (mkmem==0) rewind unkg
  if (mkmem==0.and.psps%useylm==1) rewind unylm

! PAW: retrieve Dij coefficients for this spin component
  if (psps%usepaw==1) then
   do ispden=1,nspinor**2
    do iatom=1,natom
     isp=isppol;if (nspinor==2) isp=ispden
     dimdij=paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size
     do ilmn=1,dimdij
      ekb(ilmn,iatom,ispden)=paw_ij(iatom)%dij(ilmn,isp)
     end do
     if(dimdij+1<=dimekb1) ekb(dimdij+1:dimekb1,iatom,ispden)=zero
    end do
   end do
  end if

! Loop over k points
  ikg=0
  do ikpt=1,nkpt

   nband_k=nband(ikpt+(isppol-1)*nkpt)
   istwf_k=istwfk(ikpt)
   npw_k=npwarr(ikpt)

   if(mpi_enreg%paral_compil_kpt==1)then
    if (mpi_enreg%parareel == 0) then
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) &
&     -me_distrb)) /=0) then
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
   end if ! mpi_enreg%paral_compil_kpt

   allocate(ylm_k(npw_k,mpsang*mpsang*psps%useylm))
   if (stress_needed==1) then
    if (psps%useylm==1) allocate(ylmgr_k(npw_k,3,mpsang*mpsang*psps%useylm))
    allocate(kstr1(npw_k),kstr2(npw_k))
    allocate(kstr3(npw_k),kstr4(npw_k))
    allocate(kstr5(npw_k),kstr6(npw_k))
   end if

   kpoint(:)=kpt(:,ikpt)

   kg_k(:,:) = 0
   if (mkmem==0) then

    call rdnpw(ikpt,isppol,nband_k,npw_k,nspinor,0,unkg)
!   Skip sphere data centered at k in unkg, then read k+g data
    read (unkg) ((kg_k(ii,ipw),ii=1,3),ipw=1,npw_k)

!   Read the wavefunction block for ikpt,isppol
    tim_rwwf=7
    allocate(eig_dum(mband),kg_dum(3,0),occ_dum(mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isppol,kg_dum,mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&    npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnow)
    deallocate(eig_dum,kg_dum,occ_dum)

!   Eventually read spherical harmonics
    if (psps%useylm==1) then
     read(unylm)
     if (stress_needed==1) then
      read(unylm) ((ylm_k(ipw,ilm),ipw=1,npw_k),ilm=1,mpsang*mpsang),&
&      (((ylmgr_k(ipw,ii,ilm),ipw=1,npw_k),ii=1,3),ilm=1,mpsang*mpsang)
     else
      read(unylm) ((ylm_k(ipw,ilm),ipw=1,npw_k),ilm=1,mpsang*mpsang)
     end if
    end if

   else

!   $OMP PARALLEL DO PRIVATE(ipw) &
!   $OMP&SHARED(ikg,kg,kg_k,npw_k)
    do ipw=1,npw_k
     kg_k(1,ipw)=kg(1,ipw+ikg)
     kg_k(2,ipw)=kg(2,ipw+ikg)
     kg_k(3,ipw)=kg(3,ipw+ikg)
    end do
!   $OMP END PARALLEL DO
    if (psps%useylm==1) then
!    $OMP PARALLEL DO PRIVATE(ilm,ipw) &
!    $OMP&SHARED(ikg,mpsang,npw_k,ylm,ylm_k)
     do ilm=1,mpsang*mpsang
      do ipw=1,npw_k
       ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
      end do
     end do
!    $OMP END PARALLEL DO
     if (stress_needed==1) then
!     $OMP PARALLEL DO PRIVATE(ilm,ipw) &
!     $OMP&SHARED(ikg,mpsang,npw_k,ylmgr,ylmgr_k)
      do ilm=1,mpsang*mpsang
       do ii=1,3
        do ipw=1,npw_k
         ylmgr_k(ipw,ii,ilm)=ylmgr(ipw+ikg,ii,ilm)
        end do
       end do
      end do
!     $OMP END PARALLEL DO
     end if
    end if

!   End if for choice governed by mkmem
   end if

!  Prepare kinetic contribution to stress tensor (Warning : the symmetry
!  has not been broken, like in mkkin.f or kpg3.f . It should
!  be, in order to be coherent).
!  $OMP PARALLEL DO PRIVATE(fact_kin,ipw,kgc1,kgc2,kgc3,kin,xx) &
!  $OMP&SHARED(ecut,ecutsm,ecutsm_inv) &
!  $OMP&SHARED(gprimd,htpisq,kg_k,kpoint,kstr1,kstr2,kstr3,kstr4,kstr5,kstr6,npw_k)
   if (stress_needed==1) then
    do ipw=1,npw_k
!    Compute Cartesian coordinates of (k+G)
     kgc1=gprimd(1,1)*(kpoint(1)+kg_k(1,ipw))+&
&     gprimd(1,2)*(kpoint(2)+kg_k(2,ipw))+&
&     gprimd(1,3)*(kpoint(3)+kg_k(3,ipw))
     kgc2=gprimd(2,1)*(kpoint(1)+kg_k(1,ipw))+&
&     gprimd(2,2)*(kpoint(2)+kg_k(2,ipw))+&
&     gprimd(2,3)*(kpoint(3)+kg_k(3,ipw))
     kgc3=gprimd(3,1)*(kpoint(1)+kg_k(1,ipw))+&
&     gprimd(3,2)*(kpoint(2)+kg_k(2,ipw))+&
&     gprimd(3,3)*(kpoint(3)+kg_k(3,ipw))
     kin=htpisq* ( kgc1**2 + kgc2**2 + kgc3**2 )
     fact_kin=1.0_dp
     if(kin>ecut-ecutsm)then
      if(kin>ecut)then
       fact_kin=0.0_dp
      else
!      See the routine mkkin.f, for the smearing procedure
       xx=(ecut-kin)*ecutsm_inv
!      This kinetic cutoff smoothing function and its xx derivatives
!      were produced with Mathematica and the fortran code has been
!      numerically checked against Mathematica.
       fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
       dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
!      d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
!      &                         (-144+45*xx))))))*fsm**3
       fact_kin=fsm+kin*(-ecutsm_inv)*dfsm
      end if
     end if
     kstr1(ipw)=fact_kin*kgc1*kgc1
     kstr2(ipw)=fact_kin*kgc2*kgc2
     kstr3(ipw)=fact_kin*kgc3*kgc3
     kstr4(ipw)=fact_kin*kgc3*kgc2
     kstr5(ipw)=fact_kin*kgc3*kgc1
     kstr6(ipw)=fact_kin*kgc2*kgc1
    end do ! ipw
!   $OMP END PARALLEL DO
   end if

!  Compute (k+G) vectors (only if useylm=1)
   nkpg=3*nloalg(5);allocate(kpg_k(npw_k,nkpg))
   if (nkpg>0) call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)

!  Compute nonlocal form factors ffnl at all (k+G)
!  (ider=1 computes gradients needed for stress tensor)
   ider=0;idir=0;dimffnl=1
   if (stress_needed==1) then
    ider=1;dimffnl=2+2*psps%useylm
   end if
   allocate(ffnl(npw_k,dimffnl,psps%lmnmax,ntypat))
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,&
&   gmet,gprimd,ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&   npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,&
&   psps%usepaw,psps%useylm,ylm_k,ylmgr_k)

!  Allocate the arrays phkxred and ph3d, compute phkxred and eventually ph3d
   do iatom=1,natom
    ia=atindx1(iatom)
    arg=two_pi*(kpoint(1)*xred(1,ia)+kpoint(2)*xred(2,ia)+kpoint(3)*xred(3,ia))
    phkxred(1,iatom)=cos(arg)
    phkxred(2,iatom)=sin(arg)
   end do
   if(nloalg(1)<=0)then
!   Only the allocation, not the precomputation.
    matblk=nloalg(4)
    allocate(ph3d(2,npw_k,matblk))
   else
!   Here, allocation as well as precomputation
    matblk=natom
    allocate(ph3d(2,npw_k,matblk))
    call ph1d3d(1,natom,kg_k,kpoint,matblk,natom,npw_k,n1,n2,n3,&
&    phkxred,ph1d,ph3d)
   end if

!  Loop over bands; accumulate forces and/or stresses
   do iband=1,nband_k

    if(mpi_enreg%paral_compil_kpt==1)then
!    Skip this band if not the proper processor
     if(mpi_enreg%parareel == 0) then
      if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/= me) cycle
     else
      if (mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt)/= me) cycle
     end if
    end if

!   Select occupied bands
    if( (3<=occopt.and.occopt<=7) .or. abs(occ(iband+bdtot_index))>tol8 ) then
     weight=wtk(ikpt)*occ(iband+bdtot_index)

!    Load contribution from n,k
     if(mkmem/=0)cwavef(:,1:npw_k*nspinor)=&
&     cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
     if(mkmem==0)cwavef(:,1:npw_k*nspinor)=&
&     cg_disk(:,1+(iband-1)*npw_k*nspinor:iband*npw_k*nspinor)

!    Compute non-local contributions from n,k
     if(mpi_enreg%mode_para=='b') then
      spacecomm_old=mpi_enreg%comm_fft
      mpi_enreg%comm_fft=mpi_enreg%commcart
     end if

     if (psps%usepaw==1) eig_k=eigen(iband+bdtot_index)
     call nonlop(atindx1,choice,cpopt,cprj_dum,dimekb1,dimekb2,dimffnl,dimffnl,ekb,&
&     enlout,ffnl,ffnl,gmet,gprimd,idir,psps%indlmn,istwf_k,kg_k,kg_k,kpg_k,kpg_k,kpoint,&
&     kpoint,eig_k,psps%lmnmax,matblk,mgfft,mpi_enreg,mpsang,psps%mpssoang,natom,&
&     nattyp,ngfft,nkpg,nkpg,nloalg,nnlout,npw_k,npw_k,nspinor,ntypat,0,paw_opt,&
&     phkxred,phkxred,ph1d,ph3d,ph3d,psps%pspso,signs,sij,nonlop_dum,tim_nonlop,&
&     ucvol,psps%useylm,cwavef,cwavef)

     if(mpi_enreg%mode_para=='b') then
      mpi_enreg%comm_fft=spacecomm_old
     end if

!    Accumulate non-local contributions from n,k
     if (optfor==1) grnl(1:3*natom)=grnl(1:3*natom)+weight*enlout(1+ishift:3*natom+ishift)
     if (stress_needed==1) npsstr(1:6)=npsstr(1:6)+weight*enlout(1:6)

!    Accumulate stress tensor kinetic contributions
     if (stress_needed==1) then
      call meanvalue_g(ar,kstr1,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
      kinstr(1)=kinstr(1)+weight*ar
      call meanvalue_g(ar,kstr2,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
      kinstr(2)=kinstr(2)+weight*ar
      call meanvalue_g(ar,kstr3,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
      kinstr(3)=kinstr(3)+weight*ar
      call meanvalue_g(ar,kstr4,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
      kinstr(4)=kinstr(4)+weight*ar
      call meanvalue_g(ar,kstr5,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
      kinstr(5)=kinstr(5)+weight*ar
      call meanvalue_g(ar,kstr6,0,istwf_k,mpi_enreg,npw_k,nspinor,cwavef)
      kinstr(6)=kinstr(6)+weight*ar
     end if

!    End of loop on bands
    end if
   end do

!  Incremente indexes
   bdtot_index=bdtot_index+nband_k
   if (mkmem/=0) then
    icg=icg+npw_k*nspinor*nband_k
    ikg=ikg+npw_k
   end if

   deallocate(ffnl,kpg_k,ph3d,ylm_k)
   if (stress_needed==1) then
    deallocate(kstr1,kstr2,kstr3,kstr4,kstr5,kstr6)
    if (psps%useylm==1) deallocate(ylmgr_k)
   end if

!  End k point loop
  end do
! End loop over spins
 end do

 if(mkmem==0) then
  deallocate(cg_disk)
 end if

!Parallel case: accumulate (n,k) contributions
 if( mpi_enreg%paral_compil_kpt==1) then
  write(message, '(a)' ) 'forstrnps: loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
! Forces
  if (optfor==1) then
   call timab(65,1,tsec)
   call leave_test(mpi_enreg)

   if ((mpi_enreg%paral_compil_kpt==1) .and. &
&   (mpi_enreg%paral_compil_fft==1)) then
    call xsum_mpi(grnl,mpi_enreg%comm_kpt,ierr)
   else
    call xsum_mpi(grnl,spaceComm,ierr)
   end if


   call timab(65,2,tsec)
  end if
! Stresses
  if (stress_needed==1) then
   call timab(65,1,tsec)
   call leave_test(mpi_enreg)

!  PATCH forstrnps // KPT & FFT spacecomm --> comm_kpt
   if ((mpi_enreg%paral_compil_kpt==1) .and. &
&   (mpi_enreg%paral_compil_fft==1)) then
    call xsum_mpi(kinstr,mpi_enreg%comm_kpt,ierr)
    call xsum_mpi(npsstr,mpi_enreg%comm_kpt,ierr)
   else
    call xsum_mpi(kinstr,spaceComm,ierr)
    call xsum_mpi(npsstr,spaceComm,ierr)
   end if
   call timab(65,2,tsec)
  end if
 end if

!Deallocate temporary space
 deallocate(cwavef,kg_k,phkxred,enlout,ekb)
 if (psps%usepaw==1) deallocate(sij)

!Do final normalisation of nl contribution to forces
 if (optfor==1) grnl(:)=grnl(:)/ucvol

!Do final normalizations and symmetrizations of stress tensor contributions
 if (stress_needed==1) then
  const=-(two_pi**2)/effmass/ucvol
  kinstr(:)=kinstr(:)*const
  if (nsym>1) then
   call stresssym(gprimd,nsym,kinstr,symrec)
   call stresssym(gprimd,nsym,npsstr,symrec)
  end if
 end if

!DEBUG
!write(6,*)' forstrnps : exit '
!stop
!ENDDEBUG

end subroutine forstrnps

!!***
