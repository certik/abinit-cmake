!{\src2tex{textfont=tt}}
!!****f* ABINIT/prcref_PMA
!!
!! NAME
!! prcref_PMA
!!
!! FUNCTION
!! Compute preconditioned residual potential (or density) and forces.
!! iprcel, iprcch and iprcfc govern the choice of the preconditioner.
!! Three tasks are done :
!! 1) Preconditioning of the forces (residual has already been included)
!!     using the approximate force constant matrix. Get proposed
!!     change of atomic positions.
!! 2) Precondition the residual, get first part of proposed trial
!!     potential change.
!! 3) PAW only: precondition the rhoij residuals (simple preconditionning)
!! 4) Take into account the proposed change of atomic positions to
!!     modify the proposed trial potential change.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA,XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see scfcv.f)
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam.
!!  dielstrt=number of the step at which the dielectric preconditioning begins.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | intxc=control xc quadrature
!!   | iprcch= not yet used here
!!   | iprcel= governs the preconditioning of the potential residual
!!   |    0 => simple model dielectric matrix, described by the
!!   |              parameters dielng, diemac and diemix contained in dielar.
!!   |    between 21 and 39 => until istep=dielstart, same as iprcel=0, then uses
!!   |              the RPA dielectric matrix (routine dielmt)
!!   |    between 41 and 49 => uses the RPA dielectric matrix (routine dielmt).
!!   |    between 51 and 59 => uses the RPA dielectric matrix (routine dieltcel).
!!   |    between 61 and 69 => uses the electronic dielectric matr (routine dieltcel).
!!   |    between 71 and 79 => uses the real-space preconditioner based on Kerker prc (prcrskerkerN)
!!   |    between 81 and 99 => reserved for futur version of the real-space preconditioner
!!   |    between 141 and 169 -> same as between 41 and 69 but with a different periodicity: modulo(iprcel modulo (10))
!!   | iprcfc= governs the preconditioning of the forces
!!   |         0 => hessian is the identity matrix
!!   |         1 => hessian is 0.5 times the identity matrix
!!   |         2 => hessian is 0.25 times the identity matrix
!!   | ixc=exchange-correlation choice parameter.
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | occopt=option for occupancies
!!   | pawsphmix=-PAW- preconditionning factor for the spherical part
!!   | prtvol=control print volume and debugging
!!   | typat(natom)=integer type for each atom in cell
!!  fcart(3,natom)=cartesian forces (hartree/bohr)
!!  ffttomix(nfft*(1-nfftprc/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!
!!  istep= number of the step in the SCF cycle
!!  kg_diel(3,npwdiel)=reduced planewave coordinates for the dielectric matrix.
!!  mgfft=maximum size of 1D FFTs
!!  mgfftdiel=maximum size of 1D FFTs, for the computation of the dielectric matrix
!!  moved_atm_inside= if 1, then the preconditioned forces
!!    as well as the preconditioned potential residual must be computed;
!!    otherwise, compute only the preconditioned potential residual.
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=number of fft grid points
!!  nfftprc=size of FFT grid on which the potential residual will be preconditionned
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngfftprc(18)=contain all needed information about 3D FFT for the grid corresponding to nfftprc
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  npwdiel=number of planewaves for dielectric matrix
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  optreal=1 if residual potential is is REAL space, 2 if it is in RECIPROCAL SPACE
!!  optres=0: the array vresid contains a potential residual
!!         1: the array vresid contains a density residual
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                    Use here rhoij residuals (and gradients)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for electron density in reciprocal space
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  susmat(2,npwdiel,nspden,npwdiel,nspden)=
!!   the susceptibility (or density-density response) matrix in reciprocal space
!!  vresid(optreal*nfftprc,nspden)=residual potential
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree)
!!  vhartr(nfft)=array for holding Hartree potential
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!  vpsp(nfft)=array for holding local psp
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!!  EXTRA arguments required by TFvW preconditioner
!!  deltae
!!  efermi
!!  etotal
!!  nfftf
!!  nhat
!!  nhatgr
!!  nhatgrdim
!!  optene
!!  optxc
!!  pawang
!!  pawfgrtab
!!  pawtab
!!  usexcnhat
!!  vtrial
!!
!! OUTPUT
!!  dtn_pc(3,natom)=preconditioned change of atomic position,
!!                                          in reduced coordinates
!!  vrespc(optreal*nfftprc,nspden)=preconditioned residual of the potential
!!  ==== if psps%usepaw==1
!!    rhoijrespc(npawmix)= preconditionned rhoij residuals at output
!!
!! SIDE EFFECT
!!  dielinv(2,npwdiel,nspden,npwdiel,nspden)=
!!                              inverse of the dielectric matrix in rec. space
!!  kxc(nfft,nkxc)=exchange-correlation kernel,
!!       needed if the electronic dielectric matrix is computed
!!  ===== if iprcch==3 .and. moved_atm_inside==1 =====
!!    ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!
!! PARENTS
!!      newvtr
!!
!! CHILDREN
!!      atm2fft,dielmt,dieltcel,fourdp,fresid,getph,kgindex,leave_new,mean_fftr,metric
!!      mkcore,mklocl,moddiel,rhohxc,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

  subroutine prcref_PMA(atindx,dielar,dielinv,&
&  dielstrt,dtn_pc,dtset,fcart,ffttomix,gmet,gsqcut,&
&  istep,kg_diel,kxc,lavnlr,&
&  mgfft,mgfftdiel,moved_atm_inside,mpi_enreg,&
&  nattyp,nfft,nfftprc,ngfft,ngfftprc,nkxc,npawmix,npwdiel,ntypat,n1xccc,&
&  optreal,optres,pawrhoij,ph1d,psps,rhog,rhoijrespc,rhor,rprimd,&
&  susmat,vhartr,vpsp,vresid,vrespc,vxc,xred,&
&  deltae,efermi,etotal,nfftf,nhat,nhatgr,nhatgrdim,optene,optxc,pawang,pawfgrtab,&
&  pawtab,use_lavnlr,usexcnhat,vtrial)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_12spacepar
 use interfaces_13recipspace
 use interfaces_13xc
 use interfaces_15common, except_this_one => prcref_PMA
 use interfaces_15rsprc
 use interfaces_lib01fftnew
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments-------------------------------
!variables used for tfvw
!scalars
 integer,intent(in) :: dielstrt,istep,mgfft,mgfftdiel,moved_atm_inside,n1xccc
 integer,intent(in) :: nfft,nfftf,nfftprc,nhatgrdim,nkxc,npawmix,npwdiel,ntypat
 integer,intent(in) :: optene,optreal,optres,optxc,use_lavnlr,usexcnhat
 real(dp),intent(in) :: deltae,efermi,etotal,gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),ffttomix(nfft*(1-nfftprc/nfft)),kg_diel(3,npwdiel),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),ngfftprc(18)
 real(dp),intent(in) :: dielar(7),fcart(3,dtset%natom)
 real(dp),intent(in) :: lavnlr(dtset%nfft,dtset%nspden*use_lavnlr)
 real(dp),intent(in) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: nhatgr(nfftf,dtset%nspden,3*nhatgrdim)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: susmat(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(in) :: vhartr(nfft),vresid(nfftprc*optreal,dtset%nspden)
 real(dp),intent(in) :: vtrial(dtset%nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: dielinv(2,npwdiel,dtset%nspden,npwdiel,dtset%nspden)
 real(dp),intent(inout) :: gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom),vpsp(nfft)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(out) :: dtn_pc(3,dtset%natom),rhoijrespc(npawmix)
 real(dp),intent(out) :: vrespc(nfftprc*optreal,dtset%nspden)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: cplex,dielop,i1,i2,i23,i3,iatom,ier,ifft,ig,ii,ii1,index,ing,ipw1
 integer :: ipw2,isign,ispden,klmn,kmix,mg,n1,n2,n3,n3xccc,nfftot,optatm
 integer :: optdyfr,optgr,option,optn,optn2,optstr,optv,spaceComm
 real(dp) :: ai,ar,dielng,diemac,diemac_inv,diemix,eei,enxc,factor,gqg2p3
 real(dp) :: gqgm12,gqgm13,gqgm23,gs,gs2,gs3,l2g2,length2,mixfac,numerator
 real(dp) :: ucvol,vxcavg,xnorm
 logical :: computediel
 character(len=500) :: message
 character(len=fnlen) :: filapp
!arrays
 integer :: id(3),qprtrb(3)
 integer,allocatable :: indpw_prc(:)
 real(dp) :: dielar_cp(7),dummy6(6),gprimd(3,3),qphon(3),rmet(3,3),strsxc(6)
 real(dp) :: tsec(2),vmean(dtset%nspden),vprtrb(2)
 real(dp),allocatable :: dummy(:),dyfrlo_indx(:,:,:),dyfrx2(:,:,:)
 real(dp),allocatable :: fcart_pc(:,:),gq(:,:),gresid(:,:),grtn_indx(:,:)
 real(dp),allocatable :: grxc(:,:),grxc_indx(:,:),rhog_wk(:,:),rhor_new(:,:)
 real(dp),allocatable :: rhor_wk(:,:),rhor_wk0(:,:),vhartr_wk(:),vpsp_wk(:)
 real(dp),allocatable :: vres_diel(:,:),vxc_wk(:,:),work(:),work1(:,:),work2(:),work3(:,:)
 real(dp),allocatable :: xccc3d(:),xred_wk(:,:)
 logical,allocatable :: mask(:)

! *************************************************************************

!DEBUG
!write(6,*)' prcref_PMA : enter '
!stop
!ENDDEBUG

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!1) Eventually take care of the forces

 if(moved_atm_inside==1)then
  allocate(fcart_pc(3,dtset%natom))

  if(dtset%iprcfc==0)then
   fcart_pc(:,:)=fcart(:,:)
  else
   fcart_pc(:,:)= (two**dtset%iprcfc) * fcart(:,:)
  end if

! Compute preconditioned delta xred from preconditioned fcart and rprimd
  call xredxcart(dtset%natom,-1,rprimd,fcart_pc,dtn_pc)

  deallocate(fcart_pc)
 end if

!#######################################################################

!2) Take care of the potential residual


 if(dtset%userid==3) then
  n3xccc=0;if (n1xccc/=0) n3xccc=nfft
  write(0,*) 'prcref_PMA, n3xccc=',n3xccc
  allocate(xccc3d(n3xccc))
! if(istep==10) then
  call prctfw3(deltae,dtset,&
&  efermi,etotal,gsqcut,&
&  lavnlr,mpi_enreg,&
&  nhat,nhatgr,nhatgrdim,&
&  nkxc,n3xccc,&
&  optene,optxc,&
&  pawang,pawfgrtab,pawtab,&
&  psps,rhor,rprimd,&
&  usexcnhat,&
&  vpsp,vresid,vrespc,vtrial,&
&  xccc3d,xred,istep)
  deallocate(xccc3d)
 end if


!Compute the residuals corresponding to the solution
!of an approximate realspace dielectric function according
!to X. Gonze PRB vol54 nb7 p4383 (1996)
 if(dtset%iprcel>=71.and.dtset%iprcel<=79) then
  if (nfft==nfftprc) then
   if (dtset%iprcel<=78) then
    call prcrskerker1(dtset,mpi_enreg,nfft,dtset%nspden,ngfft,dielar,etotal, &
&                     gprimd,rprimd,vresid,vrespc,dtset%natom,xred,rhor(:,1))
   else
    call prcrskerker2(dtset,nfft,dtset%nspden,ngfft,dielar,gprimd,rprimd, &
&                     vresid,vrespc,dtset%natom,xred,mpi_enreg,ucvol)
   end if
  else
!  If preconditionning has to be done on a coarse grid,
!  has to transfer several arrays
   allocate(work1(nfftprc,dtset%nspden),work3(nfftprc,dtset%nspden))
   allocate(work(2*nfftprc))
   do ispden=1,dtset%nspden
    work(:)=vresid(:,ispden)
    call fourdp(1,work,work1(:,ispden),+1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
   end do
   deallocate(work)
   if (dtset%iprcel<=78) then
    allocate(rhog_wk(2,nfftprc));rhog_wk(:,:)=zero
    if (mpi_enreg%nproc_fft>1.and. mpi_enreg%paral_compil_fft==1) then
     nfftot=ngfft(1)*ngfft(2)*ngfft(3)
     call indirect_parallel_Fourier(ffttomix,rhog_wk,mpi_enreg,ngfftprc,&
&         ngfft,nfftprc,nfft,dtset%paral_kgb,rhog,nfftot)
    else
     do ii=1,nfft
      if (ffttomix(ii)>0) rhog_wk(:,ffttomix(ii))=rhog(:,ii)
     end do
    end if
    call zerosym(rhog_wk,2,mpi_enreg,ngfftprc(1),ngfftprc(2),ngfftprc(3))
    allocate(work(nfftprc))
    call fourdp(1,rhog_wk,work,+1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
    call prcrskerker1(dtset,mpi_enreg,nfftprc,dtset%nspden,ngfftprc,dielar,etotal, &
&                     gprimd,rprimd,work1,work3,dtset%natom,xred,work)
    deallocate(work)
   else
    call prcrskerker2(dtset,nfftprc,dtset%nspden,ngfftprc,dielar,gprimd,rprimd, &
&                     work1,work3,dtset%natom,xred,mpi_enreg,ucvol)
   end if
   do ispden=1,dtset%nspden
    call fourdp(1,vrespc(:,ispden),work3(:,ispden),-1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
   end do
   deallocate(work1,work3)
  end if

 else

  if(dtset%iprcel==0 .or. (dtset%iprcel<40.and.istep<dielstrt) )then
   cplex=optreal
   qphon(:)=zero
!  Simple scalar multiplication, or model dielectric function
   call moddiel(cplex,dielar,mpi_enreg,nfftprc,ngfftprc,dtset%nspden,optreal,dtset%paral_kgb,qphon,rprimd,vresid,vrespc)

!  Use the inverse dielectric matrix in a small G sphere
  else if( (istep>=dielstrt .and. dtset%iprcel>=21) .or. modulo(dtset%iprcel,100)>=41 )then

!  With dielop=1, the matrices will be computed when istep=dielstrt
!  With dielop=2, the matrices will be computed when istep=dielstrt and 1
   dielop=1
   if(modulo(dtset%iprcel,100)>=41)dielop=2
   call testsusmat(computediel,dielop,dielstrt,dtset,istep) !test if the matrix is to be computed
   if(computediel) then
!   Compute the inverse dielectric matrix from the susceptibility matrix
!   There are two routines for the RPA matrix, while for the electronic
!   dielectric matrix, only dieltcel will do the work
    if(modulo(dtset%iprcel,100)<=49)then
     call dielmt(dielar,dielinv,dielop,gmet,kg_diel,&
&     npwdiel,dtset%nspden,dtset%occopt,dtset%prtvol,susmat)
    else
     option=1
     if(modulo(dtset%iprcel,100)>=61)option=2
     call dieltcel(dielar,dielinv,dielop,gmet,kg_diel,kxc,&
&     nfft,ngfft,nkxc,npwdiel,dtset%nspden,dtset%occopt,option,dtset%paral_kgb,dtset%prtvol,susmat)
    end if
   end if

   allocate(work1(2,nfftprc),work2(optreal*nfftprc))

!  Presently, one uses the inverse of the RPA dielectric matrix,
!  for which spin must be averaged.

!  Do fft from real space (work2) to G space (work1)
   if (optreal==1) then
    work2(:)=vresid(:,1)
!   Must average over spins if needed.
    if(dtset%nspden/=1)work2(:)=(work2(:)+vresid(:,2))*half
    call fourdp(1,work1,work2,-1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
   else
    work1(:,:)=reshape(vresid(:,1),(/2,nfftprc/))
    if (dtset%nspden/=1) work1(:,:)=(work1(:,:)+reshape(vresid(:,2),(/2,nfftprc/)))*half
   end if

!  Multiply by restricted inverse of dielectric matrix.
!  Must first copy relevant elements of work1 to a npwdiel-dimensioned array,
!  then zero work1, operate with the dielinv matrix, and store in work1.

   allocate(vres_diel(2,npwdiel),indpw_prc(npwdiel),mask(npwdiel))
   mask(:)=.true.
   call kgindex(indpw_prc,kg_diel,mask,mpi_enreg,ngfftprc,npwdiel)
   do ipw1=1,npwdiel
    if(mask(ipw1)) then
     vres_diel(1,ipw1)=work1(1,indpw_prc(ipw1))
     vres_diel(2,ipw1)=work1(2,indpw_prc(ipw1))
    end if
   end do

   work1(:,:)=zero
   do ipw1=1,npwdiel
    ar=zero ; ai=zero

!   Use inverse of dielectric matrix (potential mixing)
    if (optres==0) then
     do ipw2=1,npwdiel
      if(mask(ipw2))then
       ar=ar+dielinv(1,ipw1,1,ipw2,1)*vres_diel(1,ipw2) &
&       -dielinv(2,ipw1,1,ipw2,1)*vres_diel(2,ipw2)
       ai=ai+dielinv(2,ipw1,1,ipw2,1)*vres_diel(1,ipw2) &
&       +dielinv(1,ipw1,1,ipw2,1)*vres_diel(2,ipw2)
      end if
     end do
    else
!    Use symetric of inverse of dielectric matrix (density mixing)
     do ipw2=1,npwdiel
      if(mask(ipw2))then
       ar=ar+dielinv(1,ipw2,1,ipw1,1)*vres_diel(1,ipw2) &
&       +dielinv(2,ipw2,1,ipw1,1)*vres_diel(2,ipw2)
       ai=ai-dielinv(2,ipw2,1,ipw1,1)*vres_diel(1,ipw2) &
&       +dielinv(1,ipw2,1,ipw1,1)*vres_diel(2,ipw2)
      end if
     end do
    end if
!   Must be careful not to count the diagonal 1 twice : it is added later,
!   so must be subtracted now.
    call xsum_mpi(ar,mpi_enreg%comm_fft,ier)
    call xsum_mpi(ai,mpi_enreg%comm_fft,ier)
    if(mask(ipw1)) then
     work1(1,indpw_prc(ipw1))=ar-vres_diel(1,ipw1)
     work1(2,indpw_prc(ipw1))=ai-vres_diel(2,ipw1)
    end if !mask(ipw1)
   end do ! ipw1
   deallocate(vres_diel,indpw_prc,mask)

!  Fourier transform
   if (optreal==1) then
    call fourdp(1,work1,work2,1,mpi_enreg,nfftprc,ngfftprc,dtset%paral_kgb,0)
   else
    work2(:)=reshape(work1(:,:),(/nfftprc*2/))
   end if

!  Add to get the preconditioned vresid, must be careful about spins.
   if(dtset%iprcel>=30)then
    diemix=dielar(4)
    vrespc(:,1)=diemix*(vresid(:,1)+work2(:))
    if(dtset%nspden/=1)vrespc(:,2)=diemix*(vresid(:,2)+work2(:))
   else
    vrespc(:,1)=vresid(:,1)+work2(:)
    if(dtset%nspden/=1)vrespc(:,2)=vresid(:,2)+work2(:)
   end if

   deallocate(work1,work2)

!  Other choice ?
  else
   write(message, '(a,a,a,a,i3,a,a,a,a)' ) ch10,&
&   ' prcref_PMA : ERROR - ',ch10,&
&   '  From the calling routine, iprcel=',dtset%iprcel,ch10,&
&   '  The only allowed values are 0 or larger than 20.',ch10,&
&   '  Action : correct your input file.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
 end if
!#######################################################################

!3) PAW only : precondition the rhoij quantities (augmentation
!occupancies) residuals. Use a simple preconditionning
!with the same mixing factor as the model dielectric function.

 if (psps%usepaw==1) then
  mixfac=dtset%pawsphmix;index=0
  if (istep>=dielstrt.and.dtset%iprcel>=21.and.dtset%iprcel<30) mixfac=one
  do iatom=1,dtset%natom
   do ispden=1,dtset%nspden
    do kmix=1,pawrhoij(iatom)%lmnmix_sz
     index=index+1;klmn=pawrhoij(iatom)%kpawmix(kmix)
     rhoijrespc(index)=mixfac*pawrhoij(iatom)%rhoijres(klmn,ispden)
    end do
   end do
  end do
 end if
!#######################################################################

!4) Take care of the change of atomic positions
!Note : this part is very demanding on memory...
!however, since this algorithm is still in development,
!it was NOT included in the estimation provided by memory.f
 if(abs(dtset%iprcch)==3 .and. moved_atm_inside==1)then

! Not yet compatible with resid given in reciprocal space
  if (optreal/=1) then
   write(message, '(8a)' ) ch10,&
&   ' prcref_PMA : ERROR - ',ch10,&
&   '  From the calling routine, iprcch=3',ch10,&
&   '  You cannot use residuals in reciprocal space.',ch10,&
&   '  Action : correct your input file.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
! Not compatible with non-collinear magnetism
  if(dtset%nspden==4)then
   write(message, '(4a)' ) ch10,&
&   ' prcref_PMA : ERROR -',ch10,&
&   'iprcch=3 does not work for nspden=4 !'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

  n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
  nfftot=n1*n2*n3

  if (optres==0) then  ! Array vresid contains a potential residual
!  -----------------------------------------------------------------

!  First subtract the current local, hartree and exchange correlation potentials
   do ispden=1,dtset%nspden
    vrespc(:,ispden)=vrespc(:,ispden)-vpsp(:)-vhartr(:)-vxc(:,ispden)
   end do

!  Compute the modified density, in rhor_wk
   option=2;allocate(gresid(3,dtset%natom),grxc(3,dtset%natom))
   allocate(rhor_wk(nfft,dtset%nspden),rhor_wk0(nfft,dtset%nspden),xred_wk(3,dtset%natom))
   xred_wk(:,:)=xred(:,:)+dtn_pc(:,:)
   call fresid(dtset,gmet,gresid,gsqcut,mpi_enreg,nfft,ngfft,&
&   ntypat,option,pawtab,rhor,rprimd,&
&   ucvol,rhor_wk,xred_wk,xred,psps%znuclpsp)

!  Compute up+down rhog_wk(G) by fft
   allocate(work(nfft),rhog_wk(2,nfft))
   work(:)=rhor_wk(:,1)
   call fourdp(1,rhog_wk,work,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   deallocate(work)

!  Compute structure factor phases for new atomic pos:
   call getph(atindx,dtset%natom,n1,n2,n3,ph1d,xred_wk)

!  Compute local ionic pseudopotential vpsp:
!  and core electron density xccc3d, if needed.
   n3xccc=0;if (n1xccc/=0) n3xccc=nfft
   allocate(xccc3d(n3xccc),vpsp_wk(nfft))
   vprtrb(1:2)=zero

   if (psps%usepaw==1) then
!   PAW: compute vpsp and xccc3d together in reciprocal space
    optatm=1;optdyfr=0;optgr=0;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1
!   Note: atindx1 should be passed to atm2fft (instead of atindx) but it is unused...
    call atm2fft(atindx,xccc3d,vpsp,dummy,dummy,eei,dummy,gmet,gprimd,dummy,dummy,gsqcut,&
&    mgfft,mpi_enreg,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,ntypat,&
&    optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&    pawtab,ph1d,psps%qgrid_vl,qprtrb,dummy,dummy6,dummy6,&
&    ucvol,psps%usepaw,dummy,vprtrb,psps%vlspl)
   else
!   Norm-conserving: compute vpsp in recip. space and xccc3d in real space
    option = 1
    allocate(dyfrlo_indx(3,3,dtset%natom),grtn_indx(3,dtset%natom))
    call mklocl(dtset,dyfrlo_indx,eei,gmet,gprimd,grtn_indx,gsqcut,dummy6,&
&    mgfft,mpi_enreg,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,&
&    ntypat,option,ph1d,psps,qprtrb,rhog_wk,rhor_wk,rmet,rprimd,&
&    ucvol,vprtrb,vpsp_wk,xred)
    deallocate(dyfrlo_indx,grtn_indx)
    if (n1xccc/=0) then
     allocate(dyfrx2(3,3,dtset%natom),grxc_indx(3,dtset%natom))
     call mkcore(dummy6,dyfrx2,grxc_indx,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&     n1,n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,vxc,psps%xcccrc,&
&     psps%xccc1d,xccc3d,xred_wk)
     deallocate(dyfrx2,grxc_indx)
    end if
   end if

!  Compute Hartree+xc potentials
   allocate(vxc_wk(nfft,dtset%nspden),vhartr_wk(nfft))
   option=1
   call rhohxc(dtset,enxc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfft,ngfft,&
&   work,0,work,0,nkxc,dtset%nspden,n3xccc,option,rhog_wk,rhor_wk,rprimd,strsxc,1,&
&   vhartr_wk,vxc_wk,vxcavg,xccc3d)
   deallocate(xccc3d)

!  Sum all contributions
   do ispden=1,dtset%nspden
    do ifft=1,nfft
     vrespc(ifft,ispden)=vrespc(ifft,ispden)+vpsp_wk(ifft)+vhartr_wk(ifft)+vxc_wk(ifft,ispden)
    end do
   end do
   call mean_fftr(vrespc,vmean,mpi_enreg,nfft,nfftot,dtset%nspden)
   if(dtset%nspden==2) then
    vmean(1)=half*(vmean(1)+vmean(2))
    vmean(2)=vmean(1)
   end if
   do ispden=1,dtset%nspden
    vrespc(:,ispden)=vrespc(:,ispden)-vmean(ispden)
   end do
   deallocate(gresid,grxc,rhog_wk,rhor_wk,rhor_wk0,xred_wk,vhartr_wk,vpsp_wk,vxc_wk)

  else                 ! Array vresid contains a density residual
!  -----------------------------------------------------------------

!  Only have to compute the modified preconditionned density residual
   option=2;allocate(gresid(3,dtset%natom),grxc(3,dtset%natom))
   allocate(rhor_new(nfft,dtset%nspden),rhor_wk(nfft,dtset%nspden),rhor_wk0(nfft,dtset%nspden),xred_wk(3,dtset%natom))
   xred_wk(:,:)=xred(:,:)+dtn_pc(:,:)
   rhor_new(:,1)=rhor(:,1)+vrespc(:,1)
   if (dtset%nspden==2) then
    rhor_new(:,1)=rhor_new(:,1)+vrespc(:,2)
    rhor_new(:,2)=rhor(:,2)+vrespc(:,1)
   end if
   call fresid(dtset,gmet,gresid,gsqcut,mpi_enreg,nfft,ngfft,&
&   ntypat,option,pawtab,rhor,rprimd,&
&   ucvol,rhor_wk0,xred_wk,xred,psps%znuclpsp)
   call fresid(dtset,gmet,gresid,gsqcut,mpi_enreg,nfft,ngfft,&
&   ntypat,option,pawtab,rhor_new,rprimd,&
&   ucvol,rhor_wk,xred_wk,xred,psps%znuclpsp)
   vrespc(:,1)=rhor_wk(:,dtset%nspden)-rhor_wk0(:,dtset%nspden)
   if (dtset%nspden==2) vrespc(:,2)=rhor_wk(:,1)-rhor_wk0(:,1)-vrespc(:,1)
   deallocate(gresid,grxc,rhor_new,rhor_wk,rhor_wk0,xred_wk)
  end if

 end if

!DEBUG
!write(6,*)' prcref_PMA : exit '
!stop
!ENDDEBUG

end subroutine prcref_PMA
!!***
