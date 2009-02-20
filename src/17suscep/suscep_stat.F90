!{\src2tex{textfont=tt}}
!!****f* ABINIT/suscep_stat
!! NAME
!! suscep_stat
!!
!! FUNCTION
!! Compute the susceptibility matrix
!! from input wavefunctions, band occupations, and k point wts.
!! Include the usual sum-over-state terms, but also the
!! corrections due to the change of the Fermi level in the metallic
!! case, as well as implicit sum over higher lying conduction
!! states, thanks to the closure relation (referred to as an extrapolation).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG,AR,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=wf in G space
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions projected with non-local projectors:
!!                                   cprj_nk(i)=<p_i|Cnk> where p_i is a non-local projector.
!!  densymop_diel <type(dens_sym_operator_type)>=the density symmetrization
!!   operator for the dielectric matrix
!!  dielar(7)=input parameters for dielectric matrix and susceptibility:
!!              diecut,dielng,diemac,diemix,diegap,dielam.
!!  dielop=option for this routine (in development)
!!  dimcprj(natom*usepaw)=array of dimensions of array cprj
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt
!!           the energy for each band and k point
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  gbound_diel(2*mgfftdiel+8,2)=G sphere boundary for the dielectric matrix
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  irrzondiel(nfftdiel**(1-1/nsym),2,nspden/nsppol)=irreducible zone data
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  lmax_diel=1+max. value of l angular momentum used for dielectric matrix
!!  mband=maximum number of bands
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  mkmem=maximum number of k points in core memory
!!  mpi_enreg=informations about MPI parallelization
!!  mpw=maximum allowed value for npw
!!  natom=number of atoms in cell
!!  nband(nkpt*nsppol)=number of bands to be included in summation
!!   at each k point for each spin channel
!!  nfftdiel=number of fft grid points for the computation of the diel matrix
!!  ngfftdiel(18)=contain all needed information about 3D FFT, for dielectric matrix,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkpt=number of k points
!!  npwarr(nkpt)=number of planewaves and boundary planewaves
!!   at each k, for going from the WF sphere to the medium size FFT grid.
!!  npwdiel=third and fifth dimension of the susmat array.
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of symmetry elements in group (at least 1 for identity)
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=
!!          occupation numbers for each band (usually 2.0) at each k point
!!  occopt=option for occupancies
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  phnonsdiel(2,nfftdiel**(1-1/nsym),nspden/nsppol)=nonsymmorphic translation phases
!!  ph1ddiel(2,3*(2*mgfftdiel+1)*natom*usepaw)=one-dimensional structure factor information
!!                                             for the dielectric matrix
!!  prtvol=control print volume and debugging output
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  tnons(3,nsym)=reduced nonsymmorphic translations
!!     (symrel and tnons are in terms of real space primitive translations)
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume (Bohr**3)
!!  unkg=unit number for k+G data file
!!  unpaw=unit number for cprj PAW data (if used)
!!  usecprj= 1 if cprj array is stored in memory
!!  usepaw=flag for PAW
!!  wtk(nkpt)=k point weights (they sum to 1.0)
!!  ylmdiel(npwdiel,lmax_diel**2)= real spherical harmonics for each G and k point
!!                                 for the dielectric matrix
!!
!! OUTPUT
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!
!! SIDE EFFECTS
!!  wffnew=unit number for current wf disk file
!!
!! PARENTS
!!      suscep,vtorho
!!
!! CHILDREN
!!      chpev,cprj_alloc,cprj_diskinit,cprj_free,cprj_get,fftpac,hdr_skip
!!      leave_new,leave_test,pawgylmg,ph1d3d,rdnpw,rwwf,sphereboundary,susk
!!      suskmm,symg,symrhg,timab,wrtout,xcomm_init,xmaster_init,xme_init
!!      xsum_mpi,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine suscep_stat(atindx1,cg,cprj,densymop_diel,dielar,dielop,dimcprj,doccde,&
&  eigen,gbound_diel,gprimd,irrzondiel,istwfk,kg,&
&  kg_diel,lmax_diel,&
&  mband,mgfftdiel,mkmem,mpi_enreg,mpw,natom,nband,nfftdiel,&
&  ngfftdiel,nkpt,npwarr,&
&  npwdiel,nspden,nspinor,nsppol,nsym,ntypat,occ,occopt,&
&  paral_kgb,pawang,pawtab,phnonsdiel,ph1ddiel,prtvol,&
&  susmat,symafm,symrel,tnons,typat,ucvol,unkg,unpaw,usecprj,usepaw,wffnew,wtk,ylmdiel)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_13io_mpi
 use interfaces_13nonlocal
 use interfaces_13paw
 use interfaces_13recipspace
 use interfaces_14iowfdenpot
 use interfaces_15common
 use interfaces_17suscep, except_this_one => suscep_stat
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dielop,lmax_diel,mband,mgfftdiel,mkmem,mpw,natom
 integer,intent(in) :: nfftdiel,nkpt,npwdiel,nspden,nsppol,nsym,ntypat,occopt
 integer,intent(in) :: paral_kgb,prtvol,unkg,unpaw,usecprj,usepaw
 integer,intent(inout) :: nspinor
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dens_sym_operator_type),intent(in) :: densymop_diel
 type(pawang_type),intent(in) :: pawang
 type(wffile_type),intent(inout) :: wffnew
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom*usepaw)
 integer,intent(in) :: gbound_diel(2*mgfftdiel+8,2)
!no_abirules
!nfftdiel**(1-1/nsym) is 1 if nsym==1, and nfftdiel otherwise
 integer,intent(in) :: irrzondiel(nfftdiel**(1-1/nsym),2,nspden/nsppol),&
 & istwfk(nkpt)
 integer,intent(in) :: kg(3,mpw*mkmem),kg_diel(3,npwdiel),&
 & nband(nkpt*nsppol),ngfftdiel(18)
 integer,intent(in) :: npwarr(nkpt),symafm(nsym),symrel(3,3,nsym),typat(ntypat)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband*mkmem*nsppol),dielar(7)
 real(dp),intent(in) :: doccde(mband*nkpt*nsppol),eigen(mband*nkpt*nsppol)
 real(dp),intent(in) :: gprimd(3,3),occ(mband*nkpt*nsppol)
!nfftdiel**(1-1/nsym) is 1 if nsym==1, and nfftdiel otherwise
 real(dp),intent(in) :: phnonsdiel(2,nfftdiel**(1-1/nsym),nspden/nsppol),&
 &                                 tnons(3,nsym),wtk(nkpt)
 real(dp),intent(in) :: ph1ddiel(2,(3*(2*mgfftdiel+1)*natom)*usepaw)
 real(dp),intent(in) :: ylmdiel(npwdiel,lmax_diel**2)
 real(dp),intent(out) :: susmat(2,npwdiel,nspden,npwdiel,nspden)
 type(cprj_type) :: cprj(natom,nspinor*mband*mkmem*nsppol*usepaw*usecprj)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,diag,dielec_flag,extrap,i1,i2,i3,iband,ibg,icg,ier,ierr
 integer :: ifft,ii,ikg,ikpt,indx,iorder_cprj,ipw,ipw1,ipw2,ir,isp,isp1,isp2
 integer :: istwf_k,isym,j1,j2,j3,jj,k1,k2,k3,master,mcg,mcg_disk,me,mu,nband1
 integer :: nband_k,ndiel1,ndiel2,ndiel3,ndiel4,ndiel5,ndiel6,ngb,nkpg_diel
 integer :: npw1,npw_k,npwsp,paral_kgb_diel,spaceComm,t1,t2,testocc,tim_rwwf
 real(dp) :: ai,ar,diegap,dielam,eigdiff,eiginv,emax,invnsym,norm,normr,occdiff
 real(dp) :: phi1,phi12,phi2,phr1,phr12,phr2,sumdocc,tolocc,weight,wght1,wght2
 character(len=500) :: message
 type(MPI_type) :: mpi_enreg_diel
!arrays
 integer,allocatable :: dummy(:),gbound(:,:),kg_dum(:,:),kg_k(:,:),sym_g(:,:)
 integer,allocatable :: tmrev_g(:)
 real(dp) :: kpt_diel(3,1),tsec(2)
 real(dp),allocatable :: cg_disk(:,:),drhode(:,:,:),eig_diel(:),eig_dum(:)
 real(dp),allocatable :: gylmg_diel(:,:,:),kpg_dum(:,:),occ_deavg(:),occ_dum(:)
 real(dp),allocatable :: ph3d_diel(:,:,:),phdiel(:,:,:),phkxred_diel(:,:)
 real(dp),allocatable :: rhoextrap(:,:,:),rhoextrg(:,:),rhoextrr(:,:),sush(:)
 real(dp),allocatable :: sussum(:),susvec(:,:,:),suswk(:,:,:),zhpev1(:,:)
 real(dp),allocatable :: zhpev2(:)
 type(cprj_type),allocatable :: cprj_k(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' suscep_stat : enter '
!if(.true.)stop
!ENDDEBUG
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
 mpi_enreg_diel%paral_fft=0
 paral_kgb_diel=0

 if(nspinor==2)then
  write(6,*)' suscep_stat : does not yet work for nspinor=2, including susk, suskmm ...'
  stop
 end if

 if (usecprj==0.and.usepaw==1) then
  write (message,'(6a)')ch10,&
&  ' suscep_stat : ERROR- ',ch10,&
&  ' cprj datastructure must be allocated !',ch10,&
&  ' Action: change pawusecp input keyword.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 call timab(232,1,tsec)

!Init mpi_comm
 call xcomm_init(mpi_enreg,spaceComm)
 if(paral_kgb==1) spaceComm=mpi_enreg%comm_kpt
!Init me
 call xme_init(mpi_enreg,me)
!Init master
 call xmaster_init(mpi_enreg,master)
 mcg=mpw*nspinor*mband*mkmem*nsppol

!mkmem==0 means wf and kg info on disk file
 if (mkmem==0) then
! Skip wffnew header
  call hdr_skip(wffnew,ierr)
  mcg_disk=mpw*nspinor*mband
  allocate(cg_disk(2,mcg_disk))
! Should call xdefineOff in case of MPI I/O
 end if

!Initialize temporary file for PAW
 if (usepaw==1) then
  iorder_cprj=0 ! order for the cprj reading...
  call cprj_diskinit(atindx1,natom,iorder_cprj,mkmem,natom,dimcprj,nspinor,unpaw)
 end if

!Initialize some scalar quantities
 bdtot_index=0 ; icg=0 ; ibg=0

 ndiel1=ngfftdiel(1) ; ndiel2=ngfftdiel(2) ; ndiel3=ngfftdiel(3)

!ndiel4,ndiel5,ndiel6 are FFT dimensions, modified to avoid cache trashing
 ndiel4=ngfftdiel(4) ; ndiel5=ngfftdiel(5) ; ndiel6=ngfftdiel(6)
 diegap=dielar(5) ; dielam=dielar(6)
 extrap=0
!If dielam is too small, there is no extrapolation.
 if(dielam>1.0d-6)extrap=1

 allocate(occ_deavg(mband))
 if(occopt>=3)allocate(drhode(2,npwdiel,nspden))
 if(extrap==1)allocate(rhoextrap(ndiel4,ndiel5,ndiel6))

!zero the susceptibility matrix and other needed quantities
 susmat(:,:,:,:,:)=zero
 if(occopt>=3)then
  drhode(:,:,:)=zero
  sumdocc=zero
 end if

!testocc to be taken away
 testocc=1
!DEBUG
!write(6,*)' suscep_stat : set testocc to 0 '
!testocc=0
!ENDDEBUG

!PAW additional initiliazations
 if (usepaw==1) then
  allocate(gylmg_diel(npwdiel,lmax_diel**2,ntypat),ph3d_diel(2,npwdiel,natom),phkxred_diel(2,natom))
  kpt_diel(1:3,1)=zero;phkxred_diel(1,:)=one;phkxred_diel(2,:)=zero;nkpg_diel=0
  call pawgylmg(gprimd,gylmg_diel,kg_diel,kpg_dum,kpt_diel,lmax_diel,nkpg_diel,npwdiel,ntypat,pawtab,ylmdiel)
  call ph1d3d(1,natom,kg_diel,kpt_diel,natom,natom,npwdiel,ndiel1,ndiel2,ndiel3,phkxred_diel,ph1ddiel,ph3d_diel)
  deallocate(phkxred_diel)
 end if

!--BIG loop over spins -----------------------------------------------------

 do isp=1,nsppol

  ikg=0

  if (mkmem==0) then
!  rewind the kpgsph data file on unit unkg
   rewind (unkg)
  end if

  if(extrap==1)rhoextrap(:,:,:)=zero

! --BIG loop over k-points -------------------------------------------------

  do ikpt=1,nkpt
!  DEBUG
!  write(6,*)' suscep_stat : only one k point '
!  do ikpt=1,1
!  ENDDEBUG

   nband_k=nband(ikpt+(isp-1)*nkpt)
   istwf_k=istwfk(ikpt)
   npw_k=npwarr(ikpt)

   if(mpi_enreg%paral_compil_kpt==1)then
    if (mpi_enreg%parareel == 0) then
     if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isp) &
&     -mpi_enreg%me_kpt))/=0) then
      bdtot_index=bdtot_index+nband_k
      cycle
     end if
    else
     if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt) &
&     /= mpi_enreg%me_kpt) then
      bdtot_index=bdtot_index+nband_k
      cycle
     end if
    end if
   end if

   allocate(gbound(2*mgfftdiel+8,2),kg_k(3,npw_k))

   if (usepaw==1) then
    allocate(cprj_k(natom,nspinor*nband_k))
    call cprj_alloc(cprj_k,0,dimcprj)
    call cprj_get(atindx1,cprj_k,cprj,natom,1,ibg,ikpt,iorder_cprj,isp,&
&    mband,mkmem,mpi_enreg,natom,nband_k,nband_k,nspinor,nsppol,unpaw)
   end if

!  Do i/o as needed
   if (mkmem==0) then

    call rdnpw(ikpt,isp,nband_k,npw_k,nspinor,0,unkg)
!   Read k+g data
    read (unkg) kg_k(1:3,1:npw_k)
    call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)

!   Read the wavefunction block for ikpt,isp
    tim_rwwf=9
    allocate(eig_dum(mband),kg_dum(3,0),occ_dum(mband))
    call rwwf(cg_disk,eig_dum,0,0,0,ikpt,isp,kg_dum,mband,mcg_disk,mpi_enreg,nband_k,nband_k,&
&    npw_k,nspinor,occ_dum,-2,0,tim_rwwf,wffnew)
    deallocate(eig_dum,kg_dum,occ_dum)

   else

    kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
    call sphereboundary(gbound,istwf_k,kg_k,mgfftdiel,npw_k)

!   End test for mkmem==0
   end if

   if(extrap==1)then
!   Compute inverse of average dielectric gap for each band
!   and multiply by occupation factor
    emax=maxval(eigen(1+bdtot_index:nband_k+bdtot_index))
    do iband=1,nband_k
     occ_deavg(iband)= occ(iband+bdtot_index)*dielam &
&     / ( emax-eigen(iband+bdtot_index)  + diegap )
    end do
   else
    occ_deavg(:)=zero
   end if

!  Compute the contribution of each k-point to susmat, rhoextrap, drhode,
!  and sumdocc.
!  DEBUG seq==par: uncomment next line and comment line after
!  if(.true.) then
!  ENDDEBUG
   if(mpi_enreg%mode_para=='b')then !Only this version is in parallel
!   Use either the simpler implementation
    if(mkmem/=0)then
     call susk(atindx1,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&     gbound_diel,gylmg_diel,icg,ikpt,&
&     isp,istwfk,kg_diel,kg_k,lmax_diel,mband,mcg,mgfftdiel,mkmem,mpi_enreg,mpw,&
&     natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,ngfftdiel,nkpt,&
&     npwdiel,npw_k,nspden,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&
&     pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&     susmat,typat,ucvol,usepaw,wtk)
    else
     call susk(atindx1,bdtot_index,cg_disk,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&     gbound_diel,gylmg_diel,icg,ikpt,&
&     isp,istwfk,kg_diel,kg_k,lmax_diel,mband,mcg_disk,mgfftdiel,mkmem,mpi_enreg,mpw,&
&     natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,ngfftdiel,nkpt,&
&     npwdiel,npw_k,nspden,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,&
&     pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&     susmat,typat,ucvol,usepaw,wtk)
    end if
   else
!   Or the more sophisticated one, needed to save memory.
    if(mkmem/=0)then
     call suskmm(atindx1,bdtot_index,cg,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&     gbound_diel,gylmg_diel,icg,ikpt,&
&     isp,istwfk,kg_diel,kg_k,lmax_diel,mband,mcg,mgfftdiel,mkmem,mpi_enreg,mpw,&
&     natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,ngfftdiel,nkpt,&
&     npwdiel,npw_k,nspden,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,paral_kgb_diel,&
&     pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&     susmat,typat,ucvol,usepaw,wtk)
    else
     call suskmm(atindx1,bdtot_index,cg_disk,cprj_k,doccde,drhode,eigen,extrap,gbound,&
&     gbound_diel,gylmg_diel,icg,ikpt,&
&     isp,istwfk,kg_diel,kg_k,lmax_diel,mband,mcg_disk,mgfftdiel,mkmem,mpi_enreg,mpw,&
&     natom,nband_k,ndiel4,ndiel5,ndiel6,nfftdiel,ngfftdiel,nkpt,&
&     npwdiel,npw_k,nspden,nspinor,nsppol,ntypat,occ,occopt,occ_deavg,paral_kgb_diel,&
&     pawang,pawtab,ph3d_diel,rhoextrap,sumdocc,&
&     susmat,typat,ucvol,usepaw,wtk)
    end if
   end if
   deallocate(gbound,kg_k)

   bdtot_index=bdtot_index+nband_k

   if (mkmem/=0) then
    ibg=ibg+nspinor*nband_k
    icg=icg+npw_k*nband_k
    ikg=ikg+npw_k
   end if
   if (usepaw==1) then
    call cprj_free(cprj_k)
    deallocate(cprj_k)
   end if

!  End loop on ikpt:  --------------------------------------------------------
  end do

! Here include the contribution from the extrapolation to susmat,
! diagonal part
  if(extrap==1)then

   call timab(89,1,tsec)

!  DEBUG
!  write(6,*)' rhoextrap ='
!  if(.true.)stop
!  do i3=1,ndiel3,4
!  write(6,*)1,1,i3,rhoextrap(1,1,i3)
!  end do
!  ENDDEBUG

!  Transfer extrapolating density on augmented fft grid to
!  normal fft grid in real space. Warning : must treat only one spin
!  at a time.
   allocate(rhoextrr(nfftdiel,1),rhoextrg(2,nfftdiel))
   call fftpac(1,1,ndiel1,ndiel2,ndiel3,ndiel4,ndiel5,ndiel6,&
&   ngfftdiel,rhoextrr,rhoextrap,1)
!  Generate the density in reciprocal space, and symmetrize it
!  (note symrhg also make the reverse FFT, to get symmetrized density;
!  this is useless here, and should be made an option)
   call symrhg(1,densymop_diel,irrzondiel,mpi_enreg_diel,nfftdiel,nfftdiel,ngfftdiel,1,1,nsym,&
&   paral_kgb_diel,phnonsdiel,rhoextrg,rhoextrr,symafm)

   do ipw2=1,npwdiel
    j1=kg_diel(1,ipw2) ; j2=kg_diel(2,ipw2) ; j3=kg_diel(3,ipw2)
!   static:    Only fills lower half of the matrix (here, the susceptibility matrix)
!   dynamical: fill all, will not affect susopt==2 for which extrap==0
    do ipw1=1,npwdiel
     i1=kg_diel(1,ipw1) ; i2=kg_diel(2,ipw1) ; i3=kg_diel(3,ipw1)
!    NOTE that there is a FFT folding (superposition) bias here
!    Should use kgindex, in the same spirit as in prcref
     k1=i1-j1; k1=modulo(k1,ndiel1)
     k2=i2-j2; k2=modulo(k2,ndiel2)
     k3=i3-j3; k3=modulo(k3,ndiel3)
     ifft=k1+1+ndiel1*(k2+ndiel2*k3)
!    DEBUG
!    write(6, '(a,2i3,2es14.6)' ) &
!    &     ' i1,i2,susmat separ',ipw1,ipw2,susmat(1,ipw1,isp,ipw2,isp),&
!    &                                     susmat(2,ipw1,isp,ipw2,isp)
!    write(6, '(a,2i3,i5,2es14.6)' ) &
!    &     ' i1,i2,ifft,rhoextrg',ipw1,ipw2,ifft,rhoextrg(1,ifft),rhoextrg(2,ifft)
!    ENDDEBUG
     susmat(1,ipw1,isp,ipw2,isp)=   &
&     susmat(1,ipw1,isp,ipw2,isp)+rhoextrg(1,ifft)
     susmat(2,ipw1,isp,ipw2,isp)=   &
&     susmat(2,ipw1,isp,ipw2,isp)+rhoextrg(2,ifft)
    end do
   end do

   call timab(89,2,tsec)

   deallocate(rhoextrg,rhoextrr)

  end if

! End loop over spins ---------------------------------------------------------
 end do

 if (usepaw==1) deallocate(gylmg_diel,ph3d_diel)

 if(mpi_enreg%paral_compil_kpt==1)then
  call timab(86,1,tsec)
  if (mpi_enreg%parareel == 0) then
!  BEGIN TF_CHANGES
   call leave_test(mpi_enreg)
!  END TF_CHANGES
  end if
  write(message,*) ' suscep_stat : loop on k-points and spins done in parallel'
  call wrtout(06,message,'COLL')
  call timab(86,2,tsec)
 end if

 if(mkmem==0)deallocate(cg_disk)

 if(mpi_enreg%paral_compil_kpt==1)then
  call timab(85,1,tsec)
  allocate(sussum(2*npwdiel*nspden*npwdiel*nspden))
! Recreate full susmat on all proc.
! This should be coded more efficiently,
! since half of the matrix is still empty, and
! it is spin-diagonal.
  sussum(:)=reshape(susmat(:,:,:,:,:),(/2*npwdiel*nspden*npwdiel*nspden/))
  call xsum_mpi(sussum,spaceComm,ierr)
  susmat(:,:,:,:,:)=reshape(sussum(:),(/2,npwdiel,nspden,npwdiel,nspden/))
  deallocate(sussum)
! Recreate full drhode on all proc.
  if(occopt>=3 .and. testocc==1)then
   call xsum_mpi(drhode,spaceComm,ierr)
!  Should use only one mpi-allreduce call instead of the three
   call xsum_mpi(sumdocc,spaceComm,ierr)
  end if
  call timab(85,2,tsec)
 end if

 call timab(89,1,tsec)

 if( occopt>=3 .and. testocc==1 )then

  weight=1.0_dp/sumdocc
  do isp2=1,nspden
   do ipw2=1,npwdiel
!   Presently fills complete susceptibility matrix, not only lower half
    do isp1=1,nspden
     do ipw1=1,npwdiel
      susmat(1,ipw1,isp1,ipw2,isp2)=susmat(1,ipw1,isp1,ipw2,isp2)- &
&      weight*( drhode(1,ipw1,isp1)*drhode(1,ipw2,isp2)  &
&      +drhode(2,ipw1,isp1)*drhode(2,ipw2,isp2) )
      susmat(2,ipw1,isp1,ipw2,isp2)=susmat(2,ipw1,isp1,ipw2,isp2)- &
&      weight*( drhode(2,ipw1,isp1)*drhode(1,ipw2,isp2)  &
&      -drhode(1,ipw1,isp1)*drhode(2,ipw2,isp2) )
     end do
    end do
   end do
  end do
  deallocate(drhode)

 end if

 deallocate(occ_deavg)
 if(extrap==1)deallocate(rhoextrap)

!-The susceptibility matrix has been generated---------------------------
!-Symmetries : hermitian, time-reversal, spatial-------------------------

!Generate upper half of the matrix (still the susceptibility matrix)
 do isp=1,nspden
  do ipw2=2,npwdiel
   do ipw1=1,ipw2-1
    susmat(1,ipw1,isp,ipw2,isp)= susmat(1,ipw2,isp,ipw1,isp)
    susmat(2,ipw1,isp,ipw2,isp)=-susmat(2,ipw2,isp,ipw1,isp)
   end do
  end do
 end do

!Compute symmetric of G-vectors and eventual phases
!(either time-reversal or spatial symmetries)
 allocate(suswk(2,npwdiel,npwdiel),tmrev_g(npwdiel),sym_g(npwdiel,nsym))
 allocate(phdiel(2,npwdiel,nsym))
 call symg(kg_diel,npwdiel,nsym,phdiel,sym_g,symrel,tmrev_g,tnons)

 invnsym=1.0_dp/dble(nsym)
 do isp=1,nsppol

! Impose spatial symmetries to the susceptibility matrix
  do ipw2=1,npwdiel
   do ipw1=1,npwdiel
    ar=susmat(1,ipw1,isp,ipw2,isp)
    ai=susmat(2,ipw1,isp,ipw2,isp)
    if(nsym>1)then
     do isym=2,nsym
      t1=sym_g(ipw1,isym) ; t2=sym_g(ipw2,isym)
!     Not all symmetries are non-symmorphic. Should save time here ...
      phr1=phdiel(1,ipw1,isym) ; phi1=phdiel(2,ipw1,isym)
      phr2=phdiel(1,ipw2,isym) ; phi2=phdiel(2,ipw2,isym)
      phr12= phr1*phr2+phi1*phi2 ; phi12=phi1*phr2-phr1*phi2
      ar=ar+susmat(1,t1,isp,t2,isp)*phr12-susmat(2,t1,isp,t2,isp)*phi12
      ai=ai+susmat(2,t1,isp,t2,isp)*phr12+susmat(1,t1,isp,t2,isp)*phi12
     end do
    end if
    suswk(1,ipw1,ipw2)=ar*invnsym
    suswk(2,ipw1,ipw2)=ai*invnsym
   end do
  end do

! Impose the time-reversal symmetry to the susceptibility matrix
  do ipw2=1,npwdiel
   t2=tmrev_g(ipw2)
   do ipw1=1,npwdiel
    t1=tmrev_g(ipw1)
    susmat(1,ipw1,isp,ipw2,isp)=(suswk(1,ipw1,ipw2)+suswk(1,t1,t2))*0.5_dp
    susmat(2,ipw1,isp,ipw2,isp)=(suswk(2,ipw1,ipw2)-suswk(2,t1,t2))*0.5_dp
   end do
  end do

! End spin loop
 end do

 deallocate(suswk,phdiel,sym_g,tmrev_g)

!-The full susceptibility matrix is computed ------------------------------
!-Now, eventually diagonalize it and stop -------------------------------------

!Must turn on this flag to make the diagonalisation
 diag=0
 if(diag==1)then

  npwsp=npwdiel*nsppol
  allocate(sush(npwsp*(npwsp+1)),susvec(2,npwsp,npwsp))
  allocate(eig_diel(npwsp))
  allocate(zhpev1(2,2*npwsp-1),zhpev2(3*npwsp-2))
  ier=0

! Store the susceptibility matrix in proper mode before calling zhpev
  indx=1
  do ii=1,npwdiel
   do jj=1,ii
    sush(indx  )=susmat(1,jj,1,ii,1)
    sush(indx+1)=susmat(2,jj,1,ii,1)
    indx=indx+2
   end do
  end do

! If spin-polarized, need to store other parts of the matrix
  if(nsppol/=1)then
   do ii=1,npwdiel
!   Here, spin-flip contribution
    do jj=1,npwdiel
     sush(indx  )=susmat(1,jj,1,ii,2)
     sush(indx+1)=susmat(2,jj,1,ii,2)
     indx=indx+2
    end do
!   Here spin down-spin down upper matrix
    do jj=1,ii
     sush(indx  )=susmat(1,jj,2,ii,2)
     sush(indx+1)=susmat(2,jj,2,ii,2)
     indx=indx+2
    end do
   end do
  end if

#if defined T3E
  call CHPEV ('V','U',npwsp,sush,eig_diel,susvec,npwsp,zhpev1,&
&  zhpev2,ier)
#else
  call ZHPEV ('V','U',npwsp,sush,eig_diel,susvec,npwsp,zhpev1,&
&  zhpev2,ier)
#endif

  write(6,*)' suscep_stat : print eigenvalues of the susceptibility matrix'
  do ii=1,npwsp
   write(6, '(i5,es16.6)' )ii,eig_diel(ii)
  end do

  deallocate(sush,susvec,eig_diel,zhpev1,zhpev2)
  stop

 end if

 call timab(89,2,tsec)
 call timab(232,2,tsec)

!DEBUG
!write(6,*)' suscep_stat : exit '
!ENDDEBUG

end subroutine suscep_stat
!!***
