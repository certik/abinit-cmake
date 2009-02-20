!{\src2tex{textfont=tt}}
!!****f* ABINIT/newrho
!! NAME
!! newrho
!!
!! FUNCTION
!! Compute new trial density by mixing new and old values.
!! Call prcref to compute preconditioned residual density and forces,
!! Then, call one of the self-consistency drivers,
!! then update density.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (MT).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam.
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=
!!                              inverse of the dielectric matrix in rec. space
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | iprcch= governs the preconditioning of the atomic charges
!!   | iprcel= governs the preconditioning of the density residual
!!   | iprcfc= governs the preconditioning of the forces
!!   | iscf=( <= 0 =>non-SCF), >0 => SCF)
!!   |  iscf =11 => determination of the largest eigenvalue of the SCF cycle
!!   |  iscf =12 => SCF cycle, simple mixing
!!   |  iscf =13 => SCF cycle, Anderson mixing
!!   |  iscf =14 => SCF cycle, Anderson mixing (order 2)
!!   |  iscf =15 => SCF cycle, CG based on the minimization of the energy
!!   |  iscf =17 => SCF cycle, Pulay mixing
!!   | isecur=level of security of the computation
!!   | mffmem=governs the number of FFT arrays which are fit in core memory
!!   |          it is either 1, in which case the array f_fftgr is used,
!!   |          or 0, in which case the array f_fftgr_disk is used
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | pawoptmix=-PAW- 1 if the computed residuals include the PAW (rhoij) part
!!   | prtvol=control print volume and debugging
!!  etotal=the total energy obtained from the input density
!!  filfft=name of _FFT file
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  ffttomix(nfft*(1-nfftmix/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  ispmix=1 if mixing is done in real space, 2 if mixing is done in reciprocal space
!!  istep= number of the step in the SCF cycle
!!  i_rhor(n_index)=indices of the density in the array f_fftgr
!!  i_vresid(n_index)=indices of the density residuals in the array f_fftgr
!!  i_vrespc(n_index)=indices of the preconditioned density residuals in the array f_fftgr
!!  i_vtrial(n_index)=indices of the density in the array f_fftgr
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only for electronic
!!     dielectric matrix
!!  mgfft=maximum size of 1D FFTs
!!  mixtofft(nfftmix*(1-nfftmix/nfft))=Index of the points of the FFT grid used for mixing (coarse) on the FFT (fine) grid
!!  moved_atm_inside= if 1, then the preconditioned forces
!!    as well as the preconditioned density residual must be computed;
!!    otherwise, compute only the preconditioned density residual.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftmix=dimension of FFT grid used to mix the densities (used in PAW only)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftmix(18)=contain all needed information about 3D FFT, for the grid corresponding to nfftmix
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  npwdiel=number of planewaves for dielectric matrix
!!  nresid(nfft,nspden)=array for the residual of the density
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n_fftgr=third dimension of the array f_fftgr
!!  n_index=dimension for indices of potential/density (see i_vresid, ivrespc, i_vtrial...)
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         Use here rhoij residuals (and gradients)
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vtrial(nfft,nspden)=the trial potential that gave vresid.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!!  dtn_pc(3,natom)=preconditioned change of atomic position,
!!                                          in reduced coordinates
!!  f_atm(3,natom,n_fftgr)=different functions defined for each atom
!!  f_fftgr(ispmix*nfftmix,nspden,n_fftgr*mffmem)=different functions defined
!!    on the fft grid
!!   (see prcref, scfeig, scfopt, and scfcge for a detailed explanation).
!!   If mffmem=0, these data are kept on disk.
!!  rhor(nfft,nspden)= at input, it is the "out" trial density that gave nresid=(rho_out-rho_in)
!!                     at output, it is an updated "mixed" trial density
!!  rhog(2,nfft)= Fourier transform of the new trial density
!!  ===== if iprcch==3 .and. moved_atm_inside==1 =====
!!    ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!  ==== if usepaw==1
!!    f_paw(npawmix,n_fftgr*mffmem*usepaw)=different functions used for PAW
!!                                           (same as f_fftgr but for spherical part)
!!    pawrhoij(natom)%nrhoijsel=number of non-zero values of rhoij
!!    pawrhoij(iatom)%rhoijp(lmn2_size,nspden)= new (mixed) value of rhoij quantities in PACKED STORAGE
!!    pawrhoij(natom)%rhoijselect(lmn2_size)=select the non-zero values of rhoij
!!
!! NOTES
!!  In case of PAW calculations:
!!    Computations are done either on the fine FFT grid or the coarse grid (depending on dtset%pawmixdg)
!!    All variables (nfft,ngfft,mgfft) refer to the fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid (except f_fftgr).
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      fourdp,leave_new,metric,prcref,scfcge,scfeig,scfopt,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine newrho(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,filfft,&
&  f_atm,f_fftgr,f_paw,gmet,grhf,gsqcut,initialized,&
&  ispmix,istep,i_rhor,i_vresid,i_vrespc,i_vtrial,kg_diel,kxc,mgfft,mgfftdiel,mixtofft,&
&  moved_atm_inside,mpi_enreg,nattyp,nfft,nfftmix,ngfft,ngfftmix,nkxc,npawmix,npwdiel,&
&  nresid,ntypat,n_fftgr,n_index,n1xccc,pawrhoij,pawtab,&
&  ph1d,psps,rhog,rhor,rprimd,susmat,usepaw,vtrial,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_15common, except_this_one => newrho
!End of the abilint section

 implicit none

!Arguments-------------------------------
!scalars
 integer,intent(in) :: dielstrt,initialized,ispmix,istep,mgfft,mgfftdiel
 integer,intent(in) :: moved_atm_inside,n1xccc,n_fftgr,n_index,nfft,nfftmix
 integer,intent(in) :: nkxc,npawmix,npwdiel,ntypat,usepaw
 integer,intent(out) :: dbl_nnsclo
 real(dp),intent(in) :: etotal,gsqcut
 character(len=fnlen),intent(in) :: filfft
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),ffttomix(nfft*(1-nfftmix/nfft))
 integer,intent(in) :: kg_diel(3,npwdiel),mixtofft(nfftmix*(1-nfftmix/nfft))
 integer,intent(in) :: nattyp(ntypat),ngfft(18),ngfftmix(18)
 integer,intent(inout) :: i_rhor(n_index),i_vresid(n_index),i_vrespc(n_index)
 integer,intent(inout) :: i_vtrial(n_index)
 real(dp),intent(in) :: dielar(7),fcart(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(in) :: vtrial(nfft,dtset%nspden)
 real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(inout) :: dtn_pc(3,dtset%natom),f_atm(3,dtset%natom,n_fftgr)
 real(dp),intent(inout) :: f_fftgr(ispmix*nfftmix,dtset%nspden,n_fftgr*dtset%mffmem)
 real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*dtset%mffmem*usepaw),gmet(3,3)
 real(dp),intent(inout) :: kxc(nfft,nkxc),nresid(nfft,dtset%nspden)
 real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(out) :: rhog(2,nfft)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: i_vresid1,i_vrespc1,iatom,ifft,index,irhoij,ispden,jfft,jspden,klmn
 integer :: kmix,nfftot,nselect,response
 real(dp) :: deltae,efermi,fact,ro,ucvol
 character(len=500) :: message
 character(len=6) :: tag
!arrays
 real(dp) :: gprimd(3,3),rmet(3,3),tsec(2),vhartr_dum(1),vpsp_dum(1)
 real(dp) :: vxc_dum(1,1)
 real(dp),allocatable :: f_fftgr_disk(:,:,:),f_paw_disk(:,:),magng(:,:,:)
 real(dp),allocatable :: npaw(:),nresid0(:,:),nrespc(:,:),nreswk(:,:,:)
 real(dp),allocatable :: rhoijrespc(:),rhoijtmp(:,:),rhomag(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' newrho : enter '
!stop
!ENDDEBUG

!Compatibility tests
 if(nfftmix>nfft) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' newrho : BUG -',ch10,&
&  '  nfftmix>nfft not allowed !'
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if
 if(ispmix/=2.and.nfftmix/=nfft) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' newrho : BUG -',ch10,&
&  '  nfftmix/=nfft allowed only when ispmix=2 !'
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if

 call timab(58,1,tsec)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Select components of density to be mixed
 allocate(rhomag(ispmix*nfftmix,dtset%nspden),nresid0(ispmix*nfftmix,dtset%nspden))
 if (ispmix==1.and.nfft==nfftmix) then
  rhomag(:,1:dtset%nspden)=rhor(:,1:dtset%nspden)
  nresid0(:,1:dtset%nspden)=nresid(:,1:dtset%nspden)
 else if (nfft==nfftmix) then
  do ispden=1,dtset%nspden
   call fourdp(1,nresid0(:,ispden),nresid(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
  end do
  rhomag(:,1)=reshape(rhog,(/2*nfft/))
  if (dtset%nspden>1) then
   do ispden=2,dtset%nspden
    call fourdp(1,rhomag(:,ispden),rhor(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   end do
  end if
 else
  fact=dielar(4)-1._dp
  allocate(nreswk(2,nfft,dtset%nspden))
  do ispden=1,dtset%nspden
   call fourdp(1,nreswk(:,:,ispden),nresid(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
  end do
  do ifft=1,nfft
   if (ffttomix(ifft)>0) then
    jfft=2*ffttomix(ifft)
    rhomag (jfft-1:jfft,1)=rhog(1:2,ifft)
    nresid0(jfft-1:jfft,1)=nreswk(1:2,ifft,1)
   else
    rhog(:,ifft)=rhog(:,ifft)+fact*nreswk(:,ifft,1)
   end if
  end do
  if (dtset%nspden>1) then
   allocate(magng(2,nfft,dtset%nspden-1))
   do ispden=2,dtset%nspden
    call fourdp(1,magng(:,:,ispden-1),rhor(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
    do ifft=1,nfft
     if (ffttomix(ifft)>0) then
      jfft=2*ffttomix(ifft)
      rhomag (jfft-1:jfft,ispden)=magng (1:2,ifft,ispden-1)
      nresid0(jfft-1:jfft,ispden)=nreswk(1:2,ifft,ispden)
     else
      magng(:,ifft,ispden-1)=magng(:,ifft,ispden-1)+fact*nreswk(:,ifft,ispden)
      if (dtset%nspden==2) magng(:,ifft,1)=two*magng(:,ifft,1)-rhog(:,ifft)
     end if
    end do
   end do
  end if
  deallocate(nreswk)
 end if

!Retrieve "input" density from "output" density and density residual
 rhomag(:,1:dtset%nspden)=rhomag(:,1:dtset%nspden)-nresid0(:,1:dtset%nspden)

!If nspden==2, separate density and magnetization
 if (dtset%nspden==2) then
  rhomag (:,2)=two*rhomag (:,2)-rhomag (:,1)
  nresid0(:,2)=two*nresid0(:,2)-nresid0(:,1)
  if (usepaw==1) then
   do iatom=1,dtset%natom
    do irhoij=1,pawrhoij(iatom)%nrhoijsel
     klmn=pawrhoij(iatom)%rhoijselect(irhoij)
     ro=pawrhoij(iatom)%rhoijp(irhoij,1)
     pawrhoij(iatom)%rhoijp(irhoij,1)=ro+pawrhoij(iatom)%rhoijp(irhoij,2)
     pawrhoij(iatom)%rhoijp(irhoij,2)=ro-pawrhoij(iatom)%rhoijp(irhoij,2)
    end do
    do kmix=1,pawrhoij(iatom)%lmnmix_sz
     klmn=pawrhoij(iatom)%kpawmix(kmix)
     ro=pawrhoij(iatom)%rhoijres(klmn,1)
     pawrhoij(iatom)%rhoijres(klmn,1)=ro+pawrhoij(iatom)%rhoijres(klmn,2)
     pawrhoij(iatom)%rhoijres(klmn,2)=ro-pawrhoij(iatom)%rhoijres(klmn,2)
    end do
   end do
  end if
 end if

!Choice of preconditioner governed by iprcel, iprcch and iprcfc
 allocate(nrespc(ispmix*nfftmix,dtset%nspden))
 if (usepaw==1) allocate(npaw(npawmix),rhoijrespc(npawmix))
 call prcref(atindx,dielar,dielinv,&
& dielstrt,dtn_pc,dtset,etotal,fcart,ffttomix,gmet,gsqcut,&
& istep,kg_diel,kxc,&
& mgfft,mgfftdiel,moved_atm_inside,mpi_enreg,&
& nattyp,nfft,nfftmix,ngfft,ngfftmix,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
& ispmix,1,pawrhoij,pawtab,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&
& susmat,vhartr_dum,vpsp_dum,nresid0,nrespc,vxc_dum,xred)

!------Compute new trial density and eventual new atomic positions

 i_vresid1=i_vresid(1)
 i_vrespc1=i_vrespc(1)
 f_atm(:,:,i_vresid1)=grhf(:,:)

!Either use the array f_fftgr or the array f_fftgr_disk
 if(dtset%mffmem==1)then
  f_fftgr(:,:,i_vresid1)=nresid0(:,:)
  f_fftgr(:,:,i_vrespc1)=nrespc (:,:)
 else
! In this case, must first allocate f_fftgr_disk, then take data from disk and existing arrays.
  allocate(f_fftgr_disk(ispmix*nfftmix,dtset%nspden,n_fftgr))
  if(istep/=1)then
   call timab(83,1,tsec)
   open(unit=tmp_unit,file=filfft,form='unformatted',status='old')
   rewind(tmp_unit)
   read(tmp_unit)f_fftgr_disk
   if (usepaw==0) close(unit=tmp_unit)
   call timab(83,2,tsec)
  end if
  f_fftgr_disk(:,:,i_vresid1)=nresid0(:,:)
  f_fftgr_disk(:,:,i_vrespc1)=nrespc (:,:)
 end if
 deallocate(nresid0,nrespc)

!PAW: either use the array f_paw or the array f_paw_disk
 if (usepaw==1) then
  if(dtset%mffmem==1)then
   index=0
   do iatom=1,dtset%natom
    do ispden=1,dtset%nspden
     allocate(rhoijtmp(pawrhoij(iatom)%lmn2_size,1));rhoijtmp=zero
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      rhoijtmp(klmn,1)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
     end do
     do kmix=1,pawrhoij(iatom)%lmnmix_sz
      index=index+1;klmn=pawrhoij(iatom)%kpawmix(kmix)
      npaw(index)=rhoijtmp(klmn,1)-pawrhoij(iatom)%rhoijres(klmn,ispden)
      f_paw(index,i_vresid1)=pawrhoij(iatom)%rhoijres(klmn,ispden)
      f_paw(index,i_vrespc1)=rhoijrespc(index)
     end do
     deallocate(rhoijtmp)
    end do
   end do
  else
   allocate(f_paw_disk(npawmix,n_fftgr))
   if(istep/=1)then
    call timab(83,1,tsec)
    read(tmp_unit)f_paw_disk
    close(unit=tmp_unit)
    call timab(83,2,tsec)
   end if
   index=0
   do iatom=1,dtset%natom
    do ispden=1,dtset%nspden
     allocate(rhoijtmp(pawrhoij(iatom)%lmn2_size,1));rhoijtmp=zero
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      rhoijtmp(klmn,1)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
     end do
     do kmix=1,pawrhoij(iatom)%lmnmix_sz
      index=index+1;klmn=pawrhoij(iatom)%kpawmix(kmix)
      npaw(index)=rhoijtmp(klmn,1)-pawrhoij(iatom)%rhoijres(klmn,ispden)
      f_paw_disk(index,i_vresid1)=pawrhoij(iatom)%rhoijres(klmn,ispden)
      f_paw_disk(index,i_vrespc1)=rhoijrespc(index)
     end do
     deallocate(rhoijtmp)
    end do
   end do
  end if
 end if

!------Prediction of the components of the density

 if(dtset%iscf==11)then

! This routine compute the eigenvalues of the SCF operator
  if(dtset%mffmem==1)then
   call scfeig(f_fftgr,istep,i_vresid1,i_vrespc1,&
&   nfftmix*ispmix,dtset%nspden,n_fftgr,rhomag)
  else
   call scfeig(f_fftgr_disk,istep,i_vresid1,i_vrespc1,&
&   nfftmix*ispmix,dtset%nspden,n_fftgr,rhomag)
  end if

 else  if((dtset%iscf>=12 .and. dtset%iscf<=14).or.dtset%iscf==17) then
! Optimize next density using different algorithms, as
! determined by the variable iscf
  if(dtset%mffmem==1)then
   call scfopt(ispmix,dtn_pc,f_fftgr,f_paw,dtset%iscf-10,istep,i_vrespc,i_vtrial,&
&   moved_atm_inside,mpi_enreg,dtset%natom,nfftmix,npawmix,dtset%nspden,n_fftgr,n_index,&
&   dtset%pawoptmix,usepaw,npaw,rhomag,xred)
  else
   call scfopt(ispmix,dtn_pc,f_fftgr_disk,f_paw_disk,dtset%iscf-10,istep,i_vrespc,i_vtrial,&
&   moved_atm_inside,mpi_enreg,dtset%natom,nfftmix,npawmix,dtset%nspden,n_fftgr,n_index,&
&   dtset%pawoptmix,usepaw,npaw,rhomag,xred)
  end if
 else if(dtset%iscf==15 .or. dtset%iscf==16)then

  if(ispmix/=1) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' newrho : ERROR -',ch10,&
&   '  Mixing on reciprocal space not allowed with iscf=15 or 16.'
   call wrtout(6,message,'COLL')
   call leave_new('PERS')
  end if
! Optimize next vtrial using an algorithm based
! on the conjugate gradient minimization of etotal
  response=0;nfftot=ngfft(1)*ngfft(2)*ngfft(3)
  if(dtset%mffmem==1)then
   call scfcge(1,dbl_nnsclo,dtn_pc,etotal,f_atm,&
&   f_fftgr,initialized,dtset%iscf-10,dtset%isecur,istep,&
&   i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpi_enreg,&
&   dtset%natom,nfft,nfftot,dtset%nspden,n_fftgr,n_index,response,vtrial,ucvol,rhomag,xred)
  else
   call scfcge(1,dbl_nnsclo,dtn_pc,etotal,f_atm,&
&   f_fftgr_disk,initialized,dtset%iscf-10,dtset%isecur,istep,&
&   i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpi_enreg,&
&   dtset%natom,nfft,nfftot,dtset%nspden,n_fftgr,n_index,response,vtrial,ucvol,rhomag,xred)
  end if
! PAW: apply a simple mixing to rhoij (this is temporary)
  if (usepaw==1) then
   index=0
   do iatom=1,dtset%natom
    allocate(rhoijtmp(pawrhoij(iatom)%lmn2_size,dtset%nspden));rhoijtmp=zero
    if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
     do ispden=1,dtset%nspden
      do kmix=1,pawrhoij(iatom)%lmnmix_sz
       index=index+1;klmn=pawrhoij(iatom)%kpawmix(kmix)
       rhoijtmp(klmn,ispden)=rhoijrespc(index)-pawrhoij(iatom)%rhoijres(klmn,ispden)
      end do
     end do
    end if
    if (dtset%nspden/=2) then
     do ispden=1,dtset%nspden
      do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
       rhoijtmp(klmn,ispden)=rhoijtmp(klmn,ispden)+pawrhoij(iatom)%rhoijp(irhoij,ispden)
      end do
     end do
    else
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      ro=rhoijtmp(klmn,1)
      rhoijtmp(klmn,1)=half*(ro+rhoijtmp(klmn,2))+pawrhoij(iatom)%rhoijp(irhoij,1)
      rhoijtmp(klmn,2)=half*(ro-rhoijtmp(klmn,2))+pawrhoij(iatom)%rhoijp(irhoij,2)
     end do
    end if
    nselect=0
    do klmn=1,pawrhoij(iatom)%lmn2_size
     if (any(abs(rhoijtmp(klmn,:))>tol10)) then
      nselect=nselect+1
      pawrhoij(iatom)%rhoijselect(nselect)=klmn
      do ispden=1,dtset%nspden
       pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
      end do
     end if
    end do
    pawrhoij(iatom)%nrhoijsel=nselect
    deallocate(rhoijtmp)
   end do
  end if

 else

  write(message, '(a,a,a,a,i5,a)' ) ch10,&
&  ' newrho : BUG -',ch10,&
&  '  Invalid option: iscf =',dtset%iscf,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')

 end if

 if (usepaw==1) deallocate(rhoijrespc)

!PAW: restore rhoij from compact storage
 if (usepaw==1.and.dtset%iscf/=15.and.dtset%iscf/=16) then
  index=0
  do iatom=1,dtset%natom
   allocate(rhoijtmp(pawrhoij(iatom)%lmn2_size,dtset%nspden));rhoijtmp=zero
   if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
    do ispden=1,dtset%nspden
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      rhoijtmp(klmn,ispden)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
     end do
    end do
   end if
   do ispden=1,dtset%nspden
    do kmix=1,pawrhoij(iatom)%lmnmix_sz
     index=index+1;klmn=pawrhoij(iatom)%kpawmix(kmix)
     rhoijtmp(klmn,ispden)=npaw(index)
    end do
   end do
   if (dtset%nspden==2) then
    do irhoij=1,pawrhoij(iatom)%nrhoijsel
     klmn=pawrhoij(iatom)%rhoijselect(irhoij)
     ro=rhoijtmp(klmn,1)
     rhoijtmp(klmn,1)=half*(ro+rhoijtmp(klmn,2))
     rhoijtmp(klmn,2)=half*(ro-rhoijtmp(klmn,2))
    end do
   end if
   nselect=0
   do klmn=1,pawrhoij(iatom)%lmn2_size
    if (any(abs(rhoijtmp(klmn,:))>tol10)) then
     nselect=nselect+1
     pawrhoij(iatom)%rhoijselect(nselect)=klmn
     do ispden=1,dtset%nspden
      pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
     end do
    end if
   end do
   pawrhoij(iatom)%nrhoijsel=nselect
   deallocate(rhoijtmp)
  end do
  deallocate(npaw)
 end if

!Eventually write the data on disk and deallocate f_fftgr_disk
 if(dtset%mffmem==0)then
  call timab(83,1,tsec)
  open(unit=tmp_unit,file=filfft,form='unformatted',status='unknown')
  rewind(tmp_unit)
  write(tmp_unit)f_fftgr_disk
  if (usepaw==1) write(tmp_unit)f_paw_disk
  close(unit=tmp_unit)
  deallocate(f_fftgr_disk);if (usepaw==1) deallocate(f_paw_disk)
  call timab(83,2,tsec)
 end if

!Fourier transform the density
 if (ispmix==1.and.nfft==nfftmix) then
  rhor(:,1:dtset%nspden)=rhomag(:,1:dtset%nspden)
  call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
 else if (nfft==nfftmix) then
  do ispden=1,dtset%nspden
   call fourdp(1,rhomag(:,ispden),rhor(:,ispden),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
  end do
  rhog(:,:)=reshape(rhomag(:,1),(/2,nfft/))
 else
  do ifft=1,nfftmix
   jfft=mixtofft(ifft)
   rhog(1:2,jfft)=rhomag(2*ifft-1:2*ifft,1)
  end do
  call fourdp(1,rhog,rhor(:,1),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
  if (dtset%nspden>1) then
   do ispden=2,dtset%nspden
    do ifft=1,nfftmix
     jfft=mixtofft(ifft)
     magng(1:2,jfft,ispden-1)=rhomag(2*ifft-1:2*ifft,ispden)
    end do
    call fourdp(1,magng(:,:,ispden-1),rhor(:,ispden),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   end do
   deallocate(magng)
  end if
 end if
 deallocate(rhomag)

!Set back rho in (up+dn,up) form if nspden=2
 if (dtset%nspden==2) rhor(:,2)=half*(rhor(:,1)+rhor(:,2))

 call timab(58,2,tsec)

!DEBUG
!write(6,*)' newrho : exit '
!stop
!ENDDEBUG

end subroutine newrho
!!***
