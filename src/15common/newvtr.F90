!{\src2tex{textfont=tt}}
!!****f* ABINIT/newvtr
!! NAME
!! newvtr
!!
!! FUNCTION
!! Compute new trial potential by mixing new and old values.
!! Call prcref to compute preconditioned residual potential and forces,
!! Then, call one of the self-consistency drivers,
!! then update vtrial.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, MT)
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
!!   | fixmom=input variable that governs fixed moment calculation
!!   | intxc=control xc quadrature
!!   | iprcch= governs the preconditioning of the atomic charges
!!   | iprcel= governs the preconditioning of the potential residual
!!   | iprcfc=governs the preconditioning of the forces
!!   | iscf=( <= 0 =>non-SCF), >0 => SCF)
!!   |  iscf =1 => determination of the largest eigenvalue of the SCF cycle
!!   |  iscf =2 => SCF cycle, simple mixing
!!   |  iscf =3 => SCF cycle, Anderson mixing
!!   |  iscf =4 => SCF cycle, Anderson mixing (order 2)
!!   |  iscf =5 => SCF cycle, CG based on the minimization of the energy
!!   |  iscf =7 => SCF cycle, Pulay mixing
!!   | isecur=level of security of the computation
!!   | ixc=exchange-correlation choice parameter.
!!   | mffmem=governs the number of FFT arrays which are fit in core memory
!!   |          it is either 1, in which case the array f_fftgr is used,
!!   |          or 0, in which case the array f_fftgr_disk is used
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | occopt=option for occupancies
!!   | pawoptmix= - PAW only - 1 if the computed residuals include the PAW (rhoij) part
!!   | pawsphmix=-PAW- preconditionning factor for the spherical part
!!   | prtvol=control print volume and debugging
!!   | typat(natom)=integer type for each atom in cell
!!  etotal=the total energy obtained from the input vtrial
!!  character(len=fnlen) :: filfft=name of _FFT file
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  ffttomix(nfft*(1-nfftmix/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  ispmix=1 if mixing is done in real space, 2 if mixing is done in reciprocal space
!!  istep= number of the step in the SCF cycle
!!  i_rhor(n_index)=indices of the density in the array f_fftgr
!!  i_vresid(n_index)=indices of the potential residuals in the array f_fftgr
!!  i_vrespc(n_index)=indices of the preconditioned potential residuals in the array f_fftgr
!!  i_vtrial(n_index)=indices of the potential in the array f_fftgr
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only for electronic!
!     dielectric matrix
!!  mgfft=maximum size of 1D FFTs
!!  mixtofft(nfftmix*(1-nfftmix/nfft))=Index of the points of the FFT grid used for mixing (coarse) on the FFT (fine) grid
!!  moved_atm_inside= if 1, then the preconditioned forces
!!    as well as the preconditioned potential residual must be computed;
!!    otherwise, compute only the preconditioned potential residual.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftmix=dimension of FFT grid used to mix the densities (used in PAW only)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftmix(18)=contain all needed information about 3D FFT, for the grid corresponding to nfftmix
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  npwdiel=number of planewaves for dielectric matrix
!!  nstep=number of steps expected in iterations.
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n_fftgr=third dimension of the array f_fftgr
!!  n_index=dimension for indices of potential/density (see i_vresid, ivrespc, i_vtrial...)
!!  pawrhoij(natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         Use here rhoij residuals (and gradients)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vhartr(nfft)=array for holding Hartree potential
!!  vnew_mean(nspden)=constrained mean value of the future trial potential (might be
!!    spin-polarized
!!  vpsp(nfft)=array for holding local psp
!!  vresid(nfft,nspden)=array for the residual of the potential
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
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
!!  vtrial(nfft,nspden)= at input, it is the "in" trial potential that gave vresid=(v_out-v_in)
!!       at output, it is an updated "mixed" trial potential
!!  ===== if iprcch==3 .and. moved_atm_inside==1 =====
!!    ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!  ==== if usepaw==1
!!    f_paw(npawmix,n_fftgr*mffmem*usepaw)=different functions used for PAW
!!                                           (same as f_fftgr but for spherical part)
!!    pawrhoij(natom)%nrhoijsel,rhoijselect,rhoijp= several arrays
!!                containing new values of rhoij (augmentation occupancies)
!!
!! WARNINGS
!! depending on the value of iprcch and moved_atm_inside,
!! the xc potential or the Hxc potential may have been subtracted from vtrial !
!!
!! NOTES
!!  In case of PAW calculations:
!!  In case of PAW calculations:
!!    Computations are done either on the fine FFT grid or the coarse grid (depending on dtset%pawmixdg)
!!    All variables (nfft,ngfft,mgfft) refer to the fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!!  Subtility in PAW and non-collinear magnetism:
!!    Potentials are stored in (up-up,dn-dn,Re[up-dn],Im[up-dn]) format
!!    On-site occupancies (rhoij) are stored in (n,mx,my,mz)
!!    This is compatible provided that the mixing factors for n and m are identical
!!    and that the residual is not a combination of V_res and rhoij_res (pawoptmix=0).
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      leave_new,mean_fftr,metric,prctfw,prctfw2,prcref,scfcge,scfeig,scfopt,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine newvtr(atindx,dbl_nnsclo,dielar,dielinv,dielstrt,&
     &  dtn_pc,dtset,efermi,etotal,fcart,ffttomix,filfft,&
     &  f_atm,f_fftgr,f_paw,gmet,grhf,gsqcut,&
     &  initialized,ispmix,&
     &  istep,i_rhor,i_vresid,i_vrespc,i_vtrial,&
     &  kg_diel,kxc,mgfft,mgfftdiel,mixtofft,&
     &  moved_atm_inside,mpi_enreg,nattyp,nfft,nfftmix,&
     &  nhat,nhatgr,nhatgrdim,&
     &  ngfft,ngfftmix,nkxc,npawmix,npwdiel,&
     &  nstep,ntypat,n_fftgr,n_index,n1xccc,optres,optxc,&
     &  pawrhoij,pawang,pawfgrtab,&
     &  ph1d,&
     &  psps,rhor,rprimd,susmat,usepaw,&
     &  vhartr,vnew_mean,vpsp,vresid,&
     &  vtrial,vxc,xred,&
     &  atindx1,cg,deltae,densymop_gs,&
     &  dtfil,eeig,eew,eigen,eii,ek,enl,entropy,epaw,epawdc,irrzon,kg,&
     &  nfftf,&
     &  ngfftf,npwarr,n3xccc,occ,optene,&
     &  pawfgr,pawtab,phnons,&
     &  resid,rhog,&
     &  usexcnhat,&
     &  wffnow,&
     &  ylm,nspinor,xccc3d )

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_12spacepar
 use interfaces_15common, except_this_one => newvtr
 use interfaces_15rsprc
!End of the abilint section

 implicit none

!Arguments-------------------------------
  ! WARNING
  ! BEWARE THERE IS TWO DIFFERENT SIZE DECLARED FOR ARRAY NHAT IN RHOTOV AND RHOHXC
  ! THIS MIGHT RESULT IN A BUG
!scalars
 integer,intent(in) :: dielstrt,initialized,ispmix,istep,mgfft,mgfftdiel
 integer,intent(in) :: moved_atm_inside,n1xccc,n3xccc,n_fftgr,n_index,nfft
 integer,intent(in) :: nfftf,nfftmix,nhatgrdim,nkxc,npawmix,npwdiel,nstep
 integer,intent(in) :: ntypat,optene,optres,optxc,usepaw,usexcnhat
 integer,intent(inout) :: nspinor
 integer,intent(out) :: dbl_nnsclo
 real(dp),intent(in) :: deltae,eew,efermi,eii,entropy,epaw,epawdc,etotal,gsqcut
 real(dp),intent(out) :: eeig,ek,enl
 character(len=fnlen),intent(in) :: filfft
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(dens_sym_operator_type),intent(in) :: densymop_gs
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wffile_type),intent(inout) :: wffnow
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
 integer,intent(in) :: irrzon(nfft**(1-1/dtset%nsym),2,dtset%nspden/dtset%nsppol)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),kg_diel(3,npwdiel)
 integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft)),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),ngfftf(18),ngfftmix(18),npwarr(dtset%nkpt)
 integer,intent(inout) :: i_rhor(n_index),i_vresid(n_index),i_vrespc(n_index)
 integer,intent(inout) :: i_vtrial(n_index)
 real(dp),intent(in) :: cg(2,dtset%mpw*nspinor*dtset%mband*dtset%mkmem*dtset%nsppol)
 real(dp),intent(in) :: dielar(7),eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: fcart(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(in) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: phnons(2,nfft**(1-1/dtset%nsym),dtset%nspden/dtset%nsppol)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(in) :: vhartr(nfft),vnew_mean(dtset%nspden)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(inout) :: dtn_pc(3,dtset%natom),f_atm(3,dtset%natom,n_fftgr)
 real(dp),intent(inout) :: f_fftgr(ispmix*nfftmix,dtset%nspden,n_fftgr*dtset%mffmem)
 real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*dtset%mffmem*usepaw),gmet(3,3)
 real(dp),intent(inout) :: kxc(nfft,nkxc),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: rhog(2,nfftf),rhor(nfft,dtset%nspden),vpsp(nfft)
 real(dp),intent(inout) :: vresid(nfft,dtset%nspden),vtrial(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(out) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: compute_lavnlr,cplex,i1,i2,i23,i3,i_vresid1,i_vrespc1,iatom,ifft
 integer :: ifft2,ifft3,ifft4,ifft5,ii1,ii2,ii3,ii4,ii5,index,irhoij,ispden
 integer :: jfft,klmn,kmix,me_fft,n1,n2,n3,nfftot,nproc_fft,nselect,offset
 integer :: response,use_lavnlr
 real(dp) :: dielng,diemix,fact,r,sig,ucvol,vme,vme_inold,vme_um,x,xm,y,ym,z,zm
 character(len=500) :: message
 character(len=6) :: tag
 character(len=fnlen) :: filapp
!arrays
 real(dp),parameter :: identity(4)=(/1.0_dp,1.0_dp,0.0_dp,0.0_dp/)
 real(dp) :: gprimd(3,3),rmet(3,3),tsec(2),vmean(dtset%nspden)
 real(dp) :: vmean_inold(dtset%nspden),vmean_um(dtset%nspden)
 real(dp),allocatable :: buffer1(:,:,:),buffer2(:,:,:),dvstar(:,:)
 real(dp),allocatable :: f_fftgr_disk(:,:,:),f_paw_disk(:,:),g2cart(:)
 real(dp),allocatable :: irdiemac(:),irdiemacf1(:),irdiemacf2(:),irdiemacg(:,:)
 real(dp),allocatable :: lavnlr(:,:),ldvstar(:,:),lvres(:,:),rdiemac(:)
 real(dp),allocatable :: rdiemacf1(:),rdiemacf2(:),rdiemacg(:,:),rhoijrespc(:)
 real(dp),allocatable :: rhoijtmp(:,:),vin_old(:,:),vout_unmixed(:,:),vpaw(:)
 real(dp),allocatable :: vres(:,:),vresid0(:,:),vrespc(:,:),vreswk(:,:)
 real(dp),allocatable :: vtrial0(:,:),vtrialg(:,:,:)

! *************************************************************************
!DEBUG
!allocate(dvstar(nfft,dtset%nspden),ldvstar(nfft,dtset%nspden))
!allocate(vres(nfft,dtset%nspden),lvres(nfft,dtset%nspden),rdiemac(nfft),rdiemacg(2,nfft))
!allocate(rdiemacf1(nfft),rdiemacf2(nfft),buffer1(2,nfft,dtset%nspden),buffer2(2,nfft,dtset%nspden))
!allocate(g2cart(nfft),irdiemac(nfft),irdiemacg(2,nfft),irdiemacf1(nfft))
!allocate(irdiemacf2(nfft))
!ENDDEBUG

 if(dtset%usewvl == 1) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' newvtr : BUG -',ch10,&
&  '  dtset%usewvl == 1 not allowed (use wvl_newtr() instead)!'
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if

 if(usepaw==1.and.dtset%nspden==4.and.dtset%pawoptmix==1) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' newvtr : ERROR -',ch10,&
&  '  pawoptmix=1 is not compatible with nspden=4 !'
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if

 dielng=dielar(2)
 diemix=dielar(4)
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)

!DEBUG
!write(6,*)' newvtr : enter '
!write(6,*)' newvtr : vnew_mean(:)=',vnew_mean(:)
!write(6,*)' newvtr : vtrial(1,:)=',vtrial(1,:)
!stop
!ENDDEBUG

!Compatibility tests
 if(nfftmix>nfft) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' newvtr : BUG -',ch10,&
&  '  nfftmix>nfft not allowed !'
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if
 if(ispmix/=2.and.nfftmix/=nfft) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' newvtr : BUG -',ch10,&
&  '  nfftmix/=nfft allowed only when ispmix=2 !'
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if

 call timab(58,1,tsec)

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!------Treat the mean of potentiel residual

!Special care must be taken with components of the
!potential that are associated with NO density change.
!In general, only the global mean of the potential has
!such an anomalous feature. However, in the spin
!polarized cas with fixed occupancies, also the
!mean of each spin-potential (independently of the other)
!has such a behaviour. The trick is to remove these
!variables before going in the predictive routines,
!then to put them back

!Compute the mean of the old vtrial
 call mean_fftr(vtrial,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)

!When (collinear) spin-polarized and fixed occupation numbers,
!treat separately spin up and spin down.
!Otherwise, use only global mean
 do ispden=1,dtset%nspden
  if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&  abs(dtset%fixmom+99.99_dp)<1.0d-10)then
   vme=(vmean(1)+vmean(2))*half
  else
   vme=vmean(ispden)
  end if
  vtrial(:,ispden)=vtrial(:,ispden)-vme
 end do

!Select components of potential to be mixed
 allocate(vtrial0(ispmix*nfftmix,dtset%nspden),vresid0(ispmix*nfftmix,dtset%nspden))
 if (ispmix==1.and.nfft==nfftmix) then
  vtrial0=vtrial;vresid0=vresid
 else if (nfft==nfftmix) then
  do ispden=1,dtset%nspden
   call fourdp(1,vtrial0(:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   call fourdp(1,vresid0(:,ispden),vresid(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
  end do
 else
  allocate(vtrialg(2,nfft,dtset%nspden),vreswk(2,nfft));fact=dielar(4)
  do ispden=1,dtset%nspden
   call fourdp(1,vtrialg(:,:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   call fourdp(1,vreswk,vresid(:,ispden),-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   do ifft=1,nfft
    if (ffttomix(ifft)>0) then
     jfft=2*ffttomix(ifft)
     vtrial0(jfft-1,ispden)=vtrialg(1,ifft,ispden)
     vtrial0(jfft  ,ispden)=vtrialg(2,ifft,ispden)
     vresid0(jfft-1,ispden)=vreswk(1,ifft)
     vresid0(jfft  ,ispden)=vreswk(2,ifft)
    else
     vtrialg(:,ifft,ispden)=vtrialg(:,ifft,ispden)+fact*vreswk(:,ifft)
    end if
   end do
  end do
  deallocate(vreswk)
 end if

!Choice of preconditioner governed by iprcel, iprcch and iprcfc
 allocate(vrespc(ispmix*nfftmix,dtset%nspden))
 if (usepaw==1) allocate(vpaw(npawmix),rhoijrespc(npawmix))
!in case we are using a tfvw based preconditioner compute the local average of the non local potential
 if(dtset%userid==999) then
  compute_lavnlr = 1
  use_lavnlr=1
 else
  compute_lavnlr = 0
  use_lavnlr=0
 end if
 allocate(lavnlr(dtset%nfft,dtset%nspden*use_lavnlr))
 if(compute_lavnlr == 1) then
  call lavnl(atindx,atindx1,cg,dtfil,dtset,eigen,&
&  kg,lavnlr,mpi_enreg,&
&  nattyp,&
&  npwarr,nspinor,&
&  occ,&
&  ph1d,psps,rhor,rprimd,&
&  wffnow,xred,ylm)
 end if

 call prcref_PMA(atindx,dielar,dielinv,dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&
& istep,kg_diel,kxc,lavnlr,mgfft,mgfftdiel,moved_atm_inside,mpi_enreg,&
& nattyp,nfft,nfftmix,ngfft,ngfftmix,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
& ispmix,0,pawrhoij,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,susmat,&
& vhartr,vpsp,vresid0,vrespc,vxc,xred,&
& deltae,efermi,etotal,nfftf,nhat,nhatgr,nhatgrdim,optene,optxc,pawang,pawfgrtab,&
& pawtab,usexcnhat,use_lavnlr,vtrial )

!In case of Thomas-Fermi-von Weizsaecker charge mixing save the old trial potential for
!later use

 if(dtset%iprctfvw /= 0) then !for tfw mixing
  allocate(vin_old(nfft,dtset%nspden),vout_unmixed(nfft,dtset%nspden))
  vin_old(:,:)=vtrial(:,:)
! save the output potential before mixing for TFW charge mixing correction
! the utput potential is the sum of the input potential and the preconditionned residual.
  vout_unmixed(:,:)=vtrial(:,:)+vrespc(:,:)
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! print residual and mixed residuals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!lvres used as a buffer to hold the density residuals

!filapp="myoutput-vrespc"
!vrespc=vrespc/diemix
!call laplacian(vrespc,lvres,ngfft,rprimd)
!lvres=-lvres*(two*two_pi)
!write(filapp,770) 'myoutput-vrespc',istep
!call out1dm(filapp,natom,nfft,ngfft,dtset%nspden,ntypat,&
!&  lvres,rprimd,typat,ucvol,vrespc,xred,znucl)
!vrespc=vresid*0.01
!770 FORMAT(A15,I3.3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!------Compute new vtrial and eventual new atomic positions

 i_vresid1=i_vresid(1)
 i_vrespc1=i_vrespc(1)

 f_atm(:,:,i_vresid1)=grhf(:,:)

!Either use the array f_fftgr or the array f_fftgr_disk
 if(dtset%mffmem==1)then
  f_fftgr(:,:,i_vresid1)=vresid0(:,:)
  f_fftgr(:,:,i_vrespc1)=vrespc(:,:)
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
  f_fftgr_disk(:,:,i_vresid1)=vresid0(:,:)
  f_fftgr_disk(:,:,i_vrespc1)=vrespc(:,:)
 end if
 deallocate(vresid0)
!DEBUG
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!vres=vrespc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WARNING: affecting vrespc into rdiemacf whereas arrays are of different
!size. This may be unfunctional.
!rdiemacf1(1:ispmix*nfftmix)=vrespc(1:ispmix*nfftmix,1)
!ENDEBUG

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
      vpaw(index)=rhoijtmp(klmn,1)-pawrhoij(iatom)%rhoijres(klmn,ispden)
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
      vpaw(index)=rhoijtmp(klmn,1)-pawrhoij(iatom)%rhoijres(klmn,ispden)
      f_paw_disk(index,i_vresid1)=pawrhoij(iatom)%rhoijres(klmn,ispden)
      f_paw_disk(index,i_vrespc1)=rhoijrespc(index)
     end do
     deallocate(rhoijtmp)
    end do
   end do
  end if
 end if

!------Treat the mean of potentiel residual

!Special care must be taken with components of the
!potential that are associated with NO density change.
!In general, only the global mean of the potential has
!such an anomalous feature. However, in the spin
!polarized cas with fixed occupancies, also the
!mean of each spin-potential (independently of the other)
!has such a behaviour. The trick is to remove these
!variables before going in the predictive routines,
!then to put them back

 if(dtset%iprctfvw/=0) then
! Compute the mean of the old vtrial
  call mean_fftr(vtrial,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)
  call mean_fftr(vin_old,vmean_inold,mpi_enreg,nfft,nfftot,dtset%nspden)
  call mean_fftr(vout_unmixed,vmean_um,mpi_enreg,nfft,nfftot,dtset%nspden)
! When spin-polarized and fixed occupation numbers,
! treat separately spin up and spin down.
! Otherwise, use only global mean
  do ispden=1,dtset%nspden
   if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&   abs(dtset%fixmom+99.99_dp)<1.0d-10)then
    vme=(vmean(1)+vmean(2))*half
    vme_inold=(vmean_inold(1)+vmean_inold(2))*half
    vme_um=(vmean_um(1)+vmean_um(2))*half
   else
    vme=vmean(ispden)
    vme_inold=vmean_inold(ispden)
    vme_um=vmean_um(ispden)
   end if
   vtrial(:,ispden)=vtrial(:,ispden)-vme
   vin_old(:,ispden)=vin_old(:,ispden)-vme_inold
   vout_unmixed(:,ispden)=vout_unmixed(:,ispden)-vme_um
  end do
 else
  call mean_fftr(vtrial,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)
  do ispden=1,dtset%nspden
   if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&   abs(dtset%fixmom+99.99_dp)<1.0d-10)then
    vme=(vmean(1)+vmean(2))*half
   else
    vme=vmean(ispden)
   end if
   vtrial(:,ispden)=vtrial(:,ispden)-vme
  end do
 end if
!write(0,*) 'potentials shift',vme,vme_inold,vme_um
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$  write(0,*) 'deb1',gsqcut
!!$  xm=(rprimd(1,1)+rprimd(1,2)+rprimd(1,3))
!!$  ym=(rprimd(2,1)+rprimd(2,2)+rprimd(2,3))
!!$  zm=(rprimd(3,1)+rprimd(3,2)+rprimd(3,3))
!!$
!!$  sig=half**4
!!$  do i3=1,n3
!!$     do i2=1,n2
!!$        do i1= 1,n1
!!$           ifft=i1+n1*(i2-1+n2*(i3-1))
!!$              x=(real(real((i1-1),dp)/n1,dp))*rprimd(1,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(1,2)&
!!$                   &+(real((i3-1),dp)/real(n3,dp))*rprimd(1,3)
!!$              y=(real(real((i1-1),dp)/n1,dp))*rprimd(2,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(2,2)&
!!$                   &+(real((i3-1),dp)/real(n3,dp))*rprimd(2,3)
!!$              z=(real(real((i1-1),dp)/n1,dp))*rprimd(3,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(3,2)&
!!$                   &+(real((i3-1),dp)/real(n3,dp))*rprimd(3,3)
!!$              x=x-xm*half
!!$              y=y-ym*half
!!$              z=z-zm*half
!!$
!!$
!!$              !r=x*x+y*y+z*z*zero
!!$              !r=x*x+y*y*zero+z*z
!!$              r=x*x+y*y+z*z
!!$              !r=x*x-x*y+quarter*y*y+z*z
!!$              fr(ifft,1)=exp(-r/sig)
!!$              lfrE(ifft,1)=exp(-(r)/sig)*(&
!!$                   & -3._dp*2._dp/sig+ real(4._dp/sig**2,dp)*(r) &
!!$                   !((y-two*x)**2/sig**2-two/sig)+((x-half*y)**2/sig**2-half/sig)+(four*z*z/sig**2-two/sig)&
!!$                   &)
!!$        end do
!!$     end do
!!$  end do
!!$  write(0,*) 'deb2',gsqcut
!!$  write(0,*) 'deb3'
!!$  call laplacian(fr,lfr,ngfft,gprimd,fg,lfg,g2cart)
!!$  g2cart=-g2cart*two_pi**2
!!$  call fourdp(1,lfgE(:,:,1),lfrE(:,1),-1,nfft,ngfft,dtset%paral_kgb,0)
!!$  call fourdp(1,fgE(:,:,1),fr(:,1),-1,nfft,ngfft,dtset%paral_kgb,0)
!!$  write(0,*) 'deb 4'
!!$!  vres=vrespc
!!$!  call laplacian(vres,lvres,ngfft,gprimd)
!!$    write(0,*) 'deb5'
!!$!  rdiemac(:)=(dvstar(:,1)+(dielng)**2*(ldvstar(:,1)-lvres(:,1)))/(vresid(:,1))
!!$  write(0,*) 'deb6'
!!$  do i3=1,n3
!!$     do i2=1,n2
!!$        do i1=1,n1
!!$           ifft=i1+n1*(i2-1+n2*(i3-1))
!!$              x=(real(real((i1-1),dp)/n1,dp))*rprimd(1,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(1,2)&
!!$                   &+(real((i3-1),dp)/real(n3,dp))*rprimd(1,3)
!!$              y=(real(real((i1-1),dp)/n1,dp))*rprimd(2,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(2,2)&
!!$                   &+(real((i3-1),dp)/real(n3,dp))*rprimd(2,3)
!!$              z=(real(real((i1-1),dp)/n1,dp))*rprimd(3,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(3,2)&
!!$                   &+(real((i3-1),dp)/real(n3,dp))*rprimd(3,3)
!!$           x=x-xm*half
!!$           y=y-ym*half
!!$           z=z-zm*half
!!$           write(244,*) i1,i2,x,y,fr(ifft,1),lfr(ifft,1),lfr(ifft,1) - lfrE(ifft,1),lfrE(ifft,1),&                   !               12345678
!!$                &lfg(1,ifft,1),lfg(1,ifft,1)-lfgE(1,ifft,1),lfgE(1,ifft,1),&                                         ! real lapla g  9 10 11
!!$                &lfg(2,ifft,1),lfg(2,ifft,1)-lfgE(2,ifft,1),lfgE(2,ifft,1),&                                         ! imag lapla g 12 13 14
!!$                &fg(1,ifft,1),fg(1,ifft,1)-fgE(1,ifft,1),fgE(1,ifft,1),&                                             ! real g       15 16 17   OK
!!$                &fg(2,ifft,1),fg(2,ifft,1)-fgE(2,ifft,1),fgE(2,ifft,1), &                                            ! real g       18 19 20   OK
!!$                &g2cart(ifft),&                                                                                      ! g2*-4pi**2   21
!!$                &lfgE(1,ifft,1),lfgE(1,ifft,1)-g2cart(ifft)*fgE(1,ifft,1),g2cart(ifft)*fgE(1,ifft,1),&               ! G2*g theo    22 23 24   environ un facteur 30-40
!!$                &lfg(1,ifft,1),lfg(1,ifft,1)-g2cart(ifft)*fgE(1,ifft,1),g2cart(ifft)*fgE(1,ifft,1),&                 ! G2*g re      25 26 27  almost OK
!!$                &merge((lfg(1,ifft,1)/fg(1,ifft,1)),zero,abs(fg(1,ifft,1)).gt.1e-22_dp),&
!!$                &merge((lfg(2,ifft,1)/fg(2,ifft,1)),zero,abs(fg(2,ifft,1)).gt.1e-22_dp),&
!!$                &merge((lfgE(1,ifft,1)/fgE(1,ifft,1)),zero,abs(fgE(1,ifft,1)).gt.1e-22_dp)
!!$
!!$
!!$
!!$!rdiemac(ifft)
!!$        end do
!!$        write(244,*) ' '
!!$     end do
!!$     write(244,*) ' '
!!$     write(244,*) '# i3 #######################################################################"',i3
!!$     write(244,*) ' '
!!$     stop
!!$  end do
!!$  write(0,*) 'deb7'
!!$  stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!write(144,*) vtrial
!!$  write(0,*) 'before reading'
!!$  open(unit=70, file="vref.dat", action="read", status="old")
!!$  read(70,*) dvstar
!!$  close(70)
!!$  write(0,*) 'after reading'
!!$  dvstar=dvstar-vtrial
!!$  write(0,*) '1 op'
!!$  call mean_fftr(dvstar,vmean_um,mpi_enreg,nfft,nfftot,nspden)
!!$  dvstar=dvstar-vmean_um(1)
!!$  write(0,*) '2 op'
!!$  !call mean_fftr(vres,vmean_um,mpi_enreg,nfft,nfftot,nspden)
!!$  !vres=vres-vmean_um(1)
!!$  write(0,*) '3 op'
!!$  !call laplacian(vres,lvres,ngfft,gprimd,buffer1,buffer2,g2cart)
!!$  !call laplacian(dvstar,ldvstar,ngfft,gprimd)
!!$  !rdiemac(:)=(dvstar(:,1)-(dielng)**2*(ldvstar(:,1)-lvres(:,1)))/(vres(:,1))
!!$  !irdiemac(:)=one/rdiemac(:)
!!$  ! compute some filtered rdiemac
!!$  !call fourdp(1,rdiemacg(:,:),rdiemac(:),-1,nfft,ngfft,dtset%paral_kgb,0)
!!$  !call fourdp(1,irdiemacg(:,:),irdiemac(:),-1,nfft,ngfft,dtset%paral_kgb,0)
!!$  !do ifft=1,nfft
!!$  !   if(g2cart(ifft).gt.gsqcut*2.000000000000001_dp) then
!!$  !      rdiemacg(:,ifft)=zero
!!$  !      irdiemacg(:,ifft)=zero
!!$  !   end if
!!$  !end do
!!$  !call fourdp(1,rdiemacg(:,:),rdiemacf1(:),1,nfft,ngfft,dtset%paral_kgb,0)
!!$  !call fourdp(1,irdiemacg(:,:),irdiemacf1(:),1,nfft,ngfft,dtset%paral_kgb,0)
!!$ !do ifft=1,nfft
!!$ !    if(g2cart(ifft).gt.gsqcut*1.5_dp) then
!!$  !      rdiemacg(:,ifft)=zero
!!$  !      irdiemacg(:,ifft)=zero
!!$  !   end if
!!$  !end do
!!$  !call fourdp(1,rdiemacg(:,:),rdiemacf2(:),1,nfft,ngfft,dtset%paral_kgb,0)
!!$  !call fourdp(1,irdiemacg(:,:),irdiemacf2(:),1,nfft,ngfft,dtset%paral_kgb,0)
!!$  rdiemacf2=zero
!!$  irdiemacf1=zero
!!$  irdiemacf2=zero
!!$  rdiemac(:)= dvstar(:,1) ! to visualize the ideal residu
!!$  irdiemac(:)=vres(:,1) ! to visualize the residu
!!$!voir plus haut  ! pour voir les residus preconditionnes
!!$  write(0,*) '4 op'
!!$  ii1=1
!!$  ii2=int(n2*0.1)
!!$  ii3=int(n2*0.2)
!!$  ii4=int(n2*0.3)
!!$  ii5=int(n2*0.4)
!!$  write(0,*) '5 op'
!!$!  do i2=1,n2
!!$  i2=1
!!$  do i3= 1,n3
!!$     do i1=1,n1
!!$        ifft =i1+n1*(ii1-1+n2*(i3-1))
!!$        ifft2=i1+n1*(ii2-1+n2*(i3-1))
!!$        ifft3=i1+n1*(ii3-1+n2*(i3-1))
!!$        ifft4=i1+n1*(ii4-1+n2*(i3-1))
!!$        ifft5=i1+n1*(ii5-1+n2*(i3-1))
!!$        x=(real(real((i1-1),dp)/n1,dp))*rprimd(1,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(1,2)&
!!$             &+(real((i3-1),dp)/real(n3,dp))*rprimd(1,3)
!!$!        y=(real(real((i1-1),dp)/n1,dp))*rprimd(2,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(2,2)&
!!$!             &+(real((i3-1),dp)/real(n3,dp))*rprimd(2,3)
!!$        z=(real(real((i1-1),dp)/n1,dp))*rprimd(3,1)+(real(real((i2-1),dp)/n2,dp))*rprimd(3,2)&
!!$             &+(real((i3-1),dp)/real(n3,dp))*rprimd(3,3)
!!$        write(244,*) i1,i3,x,z,&
!!$             &rdiemac(ifft),irdiemac(ifft),&
!!$!             &rdiemac(ifft1),irdiemac(ifft1),&
!!$             &rdiemac(ifft2),irdiemac(ifft2),&
!!$             &rdiemac(ifft3),irdiemac(ifft3),&
!!$             &rdiemac(ifft4),irdiemac(ifft4),&
!!$             &rdiemac(ifft5),irdiemac(ifft5)!,&
!!$!             &rdiemac(ifft),irdiemac(ifft),&
!!$        write(344,*) i1,i3,x,z,&
!!$             &rdiemacf1(ifft),irdiemacf1(ifft),&
!!$!             &rdiemacf1(ifft1),irdiemacf1(ifft1),&
!!$             &rdiemacf1(ifft2),irdiemacf1(ifft2),&
!!$             &rdiemacf1(ifft3),irdiemacf1(ifft3),&
!!$             &rdiemacf1(ifft4),irdiemacf1(ifft4),&
!!$             &rdiemacf1(ifft5),irdiemacf1(ifft5)!,&
!!$!             &rdiemacf1(ifft),irdiemacf1(ifft),&
!!$       ! write(444,*) i1,i3,x,z,&
!!$       !      &rdiemacf2(ifft),irdiemacf2(ifft),&
!!$!      !       &rdiemacf2(ifft1),irdiemacf2(ifft1),&
!!$       !      &rdiemacf2(ifft2),irdiemacf2(ifft2),&
!!$       !      &rdiemacf2(ifft3),irdiemacf2(ifft3),&
!!$       !      &rdiemacf2(ifft4),irdiemacf2(ifft4),&
!!$       !      &rdiemacf2(ifft5),irdiemacf2(ifft5)!,&
!!$!      !       &rdiemacf2(ifft),irdiemacf2(ifft),&
!!$
!!$     end do
!!$     write(244,*) ' '
!!$     write(344,*) ' '
!!$     write(444,*) ' '
!!$  end do
!!$     write(244,*) ' '
!!$     write(244,*) '# next step #######################################################################"'
!!$     write(244,*) '###################################################################################"'
!!$     write(244,*) ' '
!!$     write(344,*) ' '
!!$     write(344,*) '# next step #######################################################################"'
!!$     write(344,*) '###################################################################################"'
!!$     write(344,*) ' '
!write(444,*) ' '
!write(444,*) '# next step #######################################################################"'
!write(444,*) '###################################################################################"'
!write(444,*) ' '
!!$
!end do





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!------Prediction of the components of the potential associated
!with a density change

 if(dtset%iscf==1)then

! This routine compute the eigenvalues of the SCF operator
  if(dtset%mffmem==1)then
   call scfeig(f_fftgr,istep,i_vresid1,i_vrespc1,&
&   ispmix*nfftmix,dtset%nspden,n_fftgr,vtrial0)
  else
   call scfeig(f_fftgr_disk,istep,i_vresid1,i_vrespc1,&
&   ispmix*nfftmix,dtset%nspden,n_fftgr,vtrial0)
  end if

 else if((dtset%iscf>=2 .and. dtset%iscf<=4).or.dtset%iscf==7)then

! Optimize next vtrial using different algorithms, as
! determined by the variable iscf
  if(dtset%mffmem==1)then
   call scfopt(ispmix,dtn_pc,f_fftgr,f_paw,dtset%iscf,istep,i_vrespc,i_vtrial,&
&   moved_atm_inside,mpi_enreg,dtset%natom,nfftmix,npawmix,dtset%nspden,n_fftgr,n_index,&
&   dtset%pawoptmix,usepaw,vpaw,vtrial0,xred)
  else
   call scfopt(ispmix,dtn_pc,f_fftgr_disk,f_paw_disk,dtset%iscf,istep,i_vrespc,i_vtrial,&
&   moved_atm_inside,mpi_enreg,dtset%natom,nfftmix,npawmix,dtset%nspden,n_fftgr,n_index,&
&   dtset%pawoptmix,usepaw,vpaw,vtrial0,xred)
  end if

 else if(dtset%iscf==5 .or. dtset%iscf==6)then

  if(ispmix/=1) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' newvtr : ERROR -',ch10,&
&   '  Mixing on reciprocal space not allowed with iscf=15 or 16.'
   call wrtout(6,message,'COLL')
   call leave_new('PERS')
  end if

! Optimize next vtrial using an algorithm based
! on the conjugate gradient minimization of etotal
  cplex=1 ; response=0
  if(dtset%mffmem==1)then
   call scfcge(cplex,dbl_nnsclo,dtn_pc,etotal,f_atm,&
&   f_fftgr,initialized,dtset%iscf,dtset%isecur,istep,&
&   i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpi_enreg,&
&   dtset%natom,nfft,nfftot,dtset%nspden,n_fftgr,n_index,response,rhor,ucvol,vtrial0,xred)
  else
   call scfcge(cplex,dbl_nnsclo,dtn_pc,etotal,f_atm,&
&   f_fftgr_disk,initialized,dtset%iscf,dtset%isecur,istep,&
&   i_rhor,i_vresid,i_vrespc,moved_atm_inside,mpi_enreg,&
&   dtset%natom,nfft,nfftot,dtset%nspden,n_fftgr,n_index,response,rhor,ucvol,vtrial0,xred)
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
    do ispden=1,dtset%nspden
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      rhoijtmp(klmn,ispden)=rhoijtmp(klmn,ispden)+pawrhoij(iatom)%rhoijp(irhoij,ispden)
     end do
    end do
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
&  ' newvtr : BUG -',ch10,&
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
     rhoijtmp(klmn,ispden)=vpaw(index)
    end do
   end do
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
  deallocate(vpaw)
 end if

!apply the Thomas--Fermi--von Weizsaecker charge mixing
!to avoid charge sloshing in large system

 if(dtset%iprctfvw /= 0.and.(istep.gt.1).and.(mod(istep,1)==0))  then
! write(3350,*) rhor
! write(3351,*) rhog
! write(3352,*) vtrial
  if(dtset%iprctfvw==1) then
   call prctfvw1(atindx,atindx1,cg,deltae,densymop_gs,dtfil,dtset,eeig,&
&   efermi,eigen,ek,enl,etotal,dtset%fixmom,gsqcut,dtset%intxc,irrzon,dtset%ixc,&
&   kg,dtset%mband,mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,nattyp,nfft,nfftf,ngfftf,&
&   nhat,nhatgr,nhatgrdim,&
&   dtset%nkpt,nkxc,npwarr,dtset%nspden,nspinor,dtset%nsppol,dtset%nsym,psps%ntypat,n3xccc,occ,dtset%occopt,&
&   optene,optres,optxc,pawang,pawfgr,pawfgrtab,pawtab,&
&   phnons,ph1d,psps,resid,rhog,rhor,rprimd,&
&   usexcnhat,&
&   vin_old,vout_unmixed,vpsp,vtrial,&
&   wffnow,xccc3d,xred,ylm,lavnlr)
  else if(dtset%iprctfvw==2.and. istep.gt.0) then
   call prctfvw2(atindx,atindx1,cg,deltae,densymop_gs,dtfil,dtset,eeig,&
&   efermi,eigen,ek,enl,etotal,dtset%fixmom,gsqcut,dtset%intxc,irrzon,dtset%ixc,&
&   kg,dtset%mband,mgfft,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%natom,nattyp,nfft,nfftf,ngfftf,&
&   nhat,nhatgr,nhatgrdim,&
&   dtset%nkpt,nkxc,npwarr,dtset%nspden,nspinor,dtset%nsppol,dtset%nsym,psps%ntypat,n3xccc,occ,dtset%occopt,&
&   optene,optres,optxc,pawang,pawfgr,pawfgrtab,pawtab,&
&   phnons,ph1d,psps,resid,rhog,rhor,rprimd,&
&   usexcnhat,&
&   vin_old,vout_unmixed,vpsp,vtrial,&
&   wffnow,xccc3d,xred,ylm)
  end if
 end if
 deallocate(lavnlr)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! print residual and mixed residuals
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!lvres used as a buffer to hold the density residuals

!filapp="myoutput-vrespc"
!vrespc=vtrial-vin_old!vrespc/diemix
!call laplacian(vrespc,lvres,ngfft,rprimd)
!lvres=-lvres/(two*two_pi)
!write(filapp,770) 'myoutput-vrespc',istep
!call out1dm(filapp,natom,nfft,ngfft,dtset%nspden,ntypat,&
!&  lvres,rprimd,typat,ucvol,vrespc,xred,znucl)
!vrespc=vresid*0.01
!vtrial=vin_old+vrespc
!770 FORMAT(A15,I3.3)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 deallocate(vrespc)

!if(istep.gt.2) stop 'fin newvtr'

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

!Restore potential
 if (ispmix==1.and.nfft==nfftmix) then
  vtrial=vtrial0
 else if (nfft==nfftmix) then
  do ispden=1,dtset%nspden
   call fourdp(1,vtrial0(:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
  end do
 else
  do ispden=1,dtset%nspden
   do ifft=1,nfftmix
    jfft=mixtofft(ifft)
    vtrialg(1,jfft,ispden)=vtrial0(2*ifft-1,ispden)
    vtrialg(2,jfft,ispden)=vtrial0(2*ifft  ,ispden)
   end do
   call fourdp(1,vtrialg(:,:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
  end do
  deallocate(vtrialg)
 end if
 deallocate(vtrial0)

!------Treat the mean of the potential

!Compute the mean of the new vtrial
 call mean_fftr(vtrial,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)

!Reset the mean of the new vtrial, to the value vnew_mean
!When spin-polarized and fixed occupation numbers,
!treat separately spin up and spin down.
!Otherwise, use only global mean
 do ispden=1,dtset%nspden
  if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&  abs(dtset%fixmom+99.99_dp)<1.0d-10)then
   vme=(vnew_mean(1)+vnew_mean(2)-vmean(1)-vmean(2))*half
  else
   vme=vnew_mean(ispden)-vmean(ispden)
  end if
  vtrial(:,ispden)=vtrial(:,ispden)+vme
 end do

!!$  !apply the Thomas--Fermi--von Weizsaecker charge mixing
!!$  !to avoid charge sloshing in large system
!!$  if(iprctfvw /= 0.and.(istep.gt.1).and.(mod(istep,1)==0))  then
!!$     !write(3350,*) rhor
!!$     !write(3351,*) rhog
!!$     !write(3352,*) vtrial
!!$     call prctfw(atindx,atindx1,cg,deltae,densymop_gs,dtfil,dtset,eeig,&
!!$          &  efermi,eigen,ek,enl,etotal,fixmom,gsqcut,intxc,irrzon,ixc,&
!!$          &  kg,mband,mgfft,mkmem,mpi_enreg,mpsang,mpw,natom,nattyp,nfft,nfftf,ngfftf,&
!!$          &  nkpt,nkxc,npwarr,nspden,nspinor,nsppol,nsym,psps%ntypat,n3xccc,occ,occopt,&
!!$          &  optene,optresn,optxc,pawfgr,pawtab,&
!!$          &  phnons,ph1d,psps,resid,rhog,rhor,rprimd,vin_old,vout_unmixed,vpsp,vtrial,&
!!$          &  wffnow,xccc3d,xred,ylm)
!!$  end if


 if(moved_atm_inside==1 .and. istep/=nstep )then
  if(abs(dtset%iprcch)==1.or.abs(dtset%iprcch)==4)then
!  Subtract current local psp, but also vxc (for core charges)
   do ispden=1,dtset%nspden
    vtrial(:,ispden)=vtrial(:,ispden)-vpsp(:)*identity(ispden)-vxc(:,ispden)
   end do
  else if(abs(dtset%iprcch)==2.or.abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)then
!  Subtract current vpsp+Hxc from vtrial. This should be rationalized later
   do ispden=1,dtset%nspden
    vtrial(:,ispden)=vtrial(:,ispden)-(vpsp(:)+vhartr(:))*identity(ispden)-vxc(:,ispden)
   end do
  end if
 end if

 call timab(58,2,tsec)

!DEBUG
!write(6,*)' newvtr : exit '
!write(6,*)' newvtr : vtrial(1,:)=',vtrial(1,:)
!stop
!ENDDEBUG
!DEBUG
!deallocate(dvstar,ldvstar)
!deallocate(vres,lvres,rdiemac,rdiemacg)
!deallocate(rdiemacf1,rdiemacf2,buffer1,buffer2)
!deallocate(g2cart,irdiemac,irdiemacg,irdiemacf1)
!deallocate(irdiemacf2)
!ENDDEBUG

end subroutine newvtr
!!***
