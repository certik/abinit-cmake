!{\src2tex{textfont=tt}}
!!****f* ABINIT/extraprho
!!
!! NAME
!! extraprho
!!
!! FUNCTION
!! Extrapolate electronic density for new ionic positions
!! from values of density of previous SCF cycle.
!! Use algorithm proposed by D. Alfe in Comp. Phys. Comm. 118 (1999), 31-33
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | densty(ntypat,4)=parameters for initialisation of the gaussian density
!!   | jellslab,slabzbeg,slabzend,slabwsrad=parameters for jellium slab
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!  gmet(3,3)=reciprocal space metric
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!  gsqcut=cutoff value on G**2 for sphere inside fft box
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  ntypat=number of types of atoms in cell
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xred_new(3,natom)=new reduced coordinates for atoms in unit cell
!!  xred_old(3,natom)=old reduced coordinates for atoms in unit cell
!!  zion(ntypat)=charge on each type of atom
!!  znucl(ntypat)=atomic number of each atom type
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= PAW rhoij occupancies and related data
!!                                         Value from previous SCF cycle is input
!!                                         Extrapolated value is output
!!  rhor(nfft,nspden)=the density from previous SCF cycle is input
!!                    the extrapolated density is output
!!  scf_history <type(scf_history_type)>=arrays obtained from previous SCF cycles
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      atm2fft,atmlength,jellium,leave_new,rhoij_alloc,wrtout
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine extraprho(atindx1,dtset,gmet,gprimd,gsqcut,mgfft,mpi_enreg,mqgrid,nattyp,&
&          nfft,ngfft,ntypat,pawrhoij,pawtab,ph1d,qgrid,rhor,rprimd,scf_history,ucvol,&
&          usepaw,xred_new,xred_old,zion,znucl)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_15common, except_this_one => extraprho
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,nfft,ntypat,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(scf_history_type),intent(inout) :: scf_history
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: qgrid(mqgrid),rprimd(3,3),zion(ntypat),znucl(ntypat)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),xred_new(3,dtset%natom)
 real(dp),intent(inout) :: xred_old(3,dtset%natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(dtset%natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: iatom,ii,ind1,ind1new,ind2,ind2new,irhoij,ispden,itypat,klmn
 integer :: lmn2_size,nselect,nspden,optatm,optdyfr,optgr,option,optn,optn2
 integer :: optstr,optv
 real(dp) :: a11,a12,a22,a33,alpha,b1,b2,beta,detA,ee,fact,ratio1,ratio2
 logical :: hasmoved,usegauss
 character(len=500) :: message
!arrays
 integer :: dummy3(3)
 integer,allocatable :: nlmn(:)
 real(dp) :: diff_t(3),diff_tmdt(3),diff_tpdt(3),dummy2(2),dummy6(6)
 real(dp),allocatable :: deltarho(:),dummy(:),gauss(:,:),rhoijtmp(:,:),work1(:)
 real(dp),allocatable :: work2(:,:),work3(:,:),xred_tpdt(:,:)

! *************************************************************************

!---------------------------------------------------------------
!----------- Inits
!---------------------------------------------------------------

!History indexes
 ind1=scf_history%hindex(1)
 ind2=scf_history%hindex(2)

!Compatibility tests
 if (ind1==0.and.ind2>0)then
  write(message, '(4a)') ch10,&
&  ' extraprho : BUG -',ch10,&
  '   Incompatible history indexes !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

!Rotated values of history indexes
 if (ind1>0.and.ind2>0) then
  ind1new=ind2;ind2new=ind1
 else if (ind1>0.and.ind2==0) then
  ind1new=3-ind1;ind2new=ind1
 else if (ind1==0.and.ind2==0) then
  ind1new=1;ind2new=0
 end if

!Compute ionic positions at t+dt in red. coordinates
!Has to take the boundary conditions into account
 allocate(xred_tpdt(3,dtset%natom))
 do iatom=1,dtset%natom
  xred_tpdt(1,iatom)=xred_old(1,iatom)+mod(xred_new(1,iatom)-xred_old(1,iatom)+1.5_dp,one)-half
  xred_tpdt(2,iatom)=xred_old(2,iatom)+mod(xred_new(2,iatom)-xred_old(2,iatom)+1.5_dp,one)-half
  xred_tpdt(3,iatom)=xred_old(3,iatom)+mod(xred_new(3,iatom)-xred_old(3,iatom)+1.5_dp,one)-half
 end do

!---------------------------------------------------------------
!----------- Compute Alpha and Beta
!----------- see (4) in Comp. Phys. Comm. 118 (1999), 31-33
!---------------------------------------------------------------

!Compute a_ij matrix
 a11=zero;a12=zero;a22=zero;a33=zero;b1=zero;b2=zero
 diff_t=zero;diff_tmdt=zero;diff_tpdt=zero
 do iatom=1,dtset%natom

  diff_tpdt(1:3)=xred_tpdt(1:3,iatom)-xred_old(1:3,iatom)
  if (ind1>0) then
   diff_t(1:3)=scf_history%xreddiff(1:3,iatom,ind1)
   if (ind2>0) diff_tmdt(1:3)=scf_history%xreddiff(1:3,iatom,ind2)
  end if
  do ii=1,3
   a11=a11+diff_t(ii)**2
   a22=a22+diff_tmdt(ii)**2
   a33=a33+diff_tpdt(ii)**2
   a12=a12+diff_t(ii)   *diff_tmdt(ii)
   b1 =b1 +diff_t(ii)   *diff_tpdt(ii)
   b2 =b2 +diff_tmdt(ii)*diff_tpdt(ii)
  end do

! Store reduced coordinates diffs in SCF history
  scf_history%xreddiff(1:3,iatom,ind1new)=diff_tpdt(1:3)

 end do
 deallocate(xred_tpdt)
 hasmoved=(a11>=tol10.or.a22>=tol10.or.a33>=tol10)

!Compute alpha and beta
 alpha=zero;beta=zero
 if (hasmoved.and.ind1>0) then
  ratio1=one;if (abs(a33)>=tol10) ratio1=(a11+a33-two*b1)/a33
  ratio2=one;if (abs(a33)>=tol10) ratio2=(a11+a33-two*b2)/a33
  detA=a11*a22-a12**2
  if (abs(a11)>=tol10.and.(abs(a22)<tol10.or.abs(detA)<tol10)) then
   alpha=b1/a11
  else if (abs(a22)>=tol10.and.(abs(a11)<tol10.or.abs(detA)<tol10)) then
   beta=b2/a22
  else if (abs(ratio1)+abs(ratio2)<tol6) then
   if (ind2>0) then
    alpha=two;beta=-one
   else
    alpha=one
   end if
   write(message, '(6a,f4.1,a,f4.1)') ch10,&
&   ' extraprho : WARNING -',ch10,&
&   '   Ionic positions lead to a collinear system !',ch10,&
&   '   Mixing coeffs have been set to: alpha=',alpha,' beta=',beta
   call wrtout(6,message,'COLL')
  else if (abs(a11)>=tol10.and.abs(a22)>=tol10) then
   alpha=(b1*a22-b2*a12)/detA
   beta =(b2*a11-b1*a12)/detA
  end if
 end if

!---------------------------------------------------------------
!----------- Contribution from delta_rho(t), delta_rho(t-dt)
!----------- and delta_rho(t-2dt) to predicted rho(t+dt)
!---------------------------------------------------------------

!deltarho(t+dt) <- deltarho(t) + alpha.[deltarho(t)-deltarho(t-dt)]
!+ beta .[deltarho(t-dt)-deltarho(t-2dt)]
!Note: scf_history%deltarhor is updated at the same time

 allocate(deltarho(nfft))
 do ispden=1,dtset%nspden

  if (ispden==1) then
   deltarho(:)=rhor(:,ispden)-scf_history%atmrho_last(:)
  else if (ispden==2.and.dtset%nspden==2) then
   deltarho(:)=rhor(:,ispden)-half*scf_history%atmrho_last(:)
  end if


! rho(t+dt) <- deltarho(t) + alpha.deltarho(t)
  if (dtset%nspden/=4.or.ispden==1) then
   rhor(:,ispden)=(one+alpha)*deltarho(:)
  else
   rhor(:,ispden)=(one+alpha)*rhor(:,ispden)
  end if

  if (hasmoved) then

!  rho(t+dt) <- -alpha.deltarho(t-dt) + beta.deltarho(t-dt)
   if (abs(beta-alpha)>tol14.and.ind1>0) then
    rhor(:,ispden)=rhor(:,ispden)+(beta-alpha)*scf_history%deltarhor(:,ispden,ind1)
   end if

!  rho(t+dt) <- -beta.deltarho(t-2dt)
   if (abs(beta)>tol14.and.ind2>0) then
    rhor(:,ispden)=rhor(:,ispden)-beta*scf_history%deltarhor(:,ispden,ind2)
   end if

  end if

! Store deltarho(t) in history
  if (dtset%nspden/=4.or.ispden==1) then
   scf_history%deltarhor(:,ispden,ind1new)=deltarho(:)
  else
   scf_history%deltarhor(:,ispden,ind1new)=rhor(:,ispden)
  end if

 end do

 deallocate(deltarho)

!---------------------------------------------------------------
!----------- Contribution from rho_at(t+dt) to predicted rho(t+dt)
!---------------------------------------------------------------

!Determine whether a gaussian atomic density has to be used or not
 usegauss=.true.
 if (usepaw==1) usegauss=(minval(pawtab(1:ntypat)%usetvale)==0)
 if (usegauss) then
  optn2=3;allocate(gauss(2,ntypat))
  do itypat=1,ntypat
   gauss(1,itypat)=zion(itypat)
   call atmlength(dtset%densty(itypat,1),gauss(2,itypat),zion(itypat),znucl(itypat))
  end do
 else
  optn2=2;allocate(gauss(2,0))
 end if

!Compute rho_at(t+dt) as sum of atomic densities
!Note: scf_history%atmrho_last is updated at the same time
 optatm=1;optdyfr=0;optgr=0;optstr=0;optv=0;optn=1
 call atm2fft(atindx1,scf_history%atmrho_last,dummy,dummy,dummy,ee,gauss,gmet,gprimd,dummy,dummy,gsqcut,&
& mgfft,mpi_enreg,mqgrid,dtset%natom,nattyp,nfft,ngfft,ntypat,&
& optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
& pawtab,ph1d,qgrid,dummy3,dummy,dummy6,dummy6,&
& ucvol,usepaw,dummy,dummy2,dummy)
 deallocate(gauss)

!Take eventually into account jellium slab
 if (dtset%jellslab/=0) then
  option=2; allocate(work1(nfft),work2(nfft,1),work3(2,nfft))
  work2(:,1)=scf_history%atmrho_last(:)
  call jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,1,option,dtset%paral_kgb,&
&  dtset%slabwsrad,work3,work2,rprimd,work1,dtset%slabzbeg,dtset%slabzend)
  scf_history%atmrho_last(:)=work2(:,1)
  deallocate(work1,work2,work3)
 end if

!Add rho_at(t+dt) to rho(t+dt)
 rhor(:,1)=rhor(:,1)+scf_history%atmrho_last(:)
 if (dtset%nspden==2) rhor(:,2)=rhor(:,2)+half*scf_history%atmrho_last(:)

!---------------------------------------------------------------
!----------- Extrapolation of PAW rhoij occupancy matrixes
!---------------------------------------------------------------

 if (usepaw==1) then

  if (ind2==0) then
   allocate(nlmn(ntypat))
   do itypat=1,ntypat
    nlmn(itypat)=pawtab(itypat)%lmn_size
   end do
   call rhoij_alloc(1,nlmn,dtset%nspden,dtset%nsppol,scf_history%pawrhoij(:,ind1new),dtset%typat)
   deallocate(nlmn)
  end if

  do iatom=1,dtset%natom
   nspden=pawrhoij(iatom)%nspden
   lmn2_size=pawrhoij(iatom)%lmn2_size

   if (hasmoved) then
    allocate(rhoijtmp(lmn2_size,nspden));rhoijtmp=zero
    
    do ispden=1,nspden

!    rhoij(t+dt) <- rhoij(t) + alpha.rhoij(t)
     fact=one+alpha
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      rhoijtmp(klmn,ispden)=rhoijtmp(klmn,ispden)+fact*pawrhoij(iatom)%rhoijp(irhoij,ispden)
     end do

!    rhoij(t+dt) <- -alpha.rhoij(t-dt) + beta.rhoij(t-dt)
     if (abs(beta-alpha)>tol14.and.ind1>0) then
      fact=beta-alpha
      do irhoij=1,scf_history%pawrhoij(iatom,ind1)%nrhoijsel
       klmn=scf_history%pawrhoij(iatom,ind1)%rhoijselect(irhoij)
       rhoijtmp(klmn,ispden)=rhoijtmp(klmn,ispden)+fact*scf_history%pawrhoij(iatom,ind1)%rhoijp(irhoij,ispden)
      end do
     end if

!    rho(t+dt) <- -beta.rhoij(t-2dt)
     if (abs(beta)>tol14.and.ind2>0) then
      fact=-beta
      do irhoij=1,scf_history%pawrhoij(iatom,ind2)%nrhoijsel
       klmn=scf_history%pawrhoij(iatom,ind2)%rhoijselect(irhoij)
       rhoijtmp(klmn,ispden)=rhoijtmp(klmn,ispden)+fact*scf_history%pawrhoij(iatom,ind2)%rhoijp(irhoij,ispden)
      end do
     end if

    end do !ispden
   end if !hasmoved

!  Store rhoij(t) in history
!  (cannot use rhoij_copy here because update for single atom)
   nselect=pawrhoij(iatom)%nrhoijsel
   scf_history%pawrhoij(iatom,ind1new)%nrhoijsel=nselect
   scf_history%pawrhoij(iatom,ind1new)%rhoijselect(1:nselect)=pawrhoij(iatom)%rhoijselect(1:nselect)
   scf_history%pawrhoij(iatom,ind1new)%rhoijp(1:nselect,1:nspden)=pawrhoij(iatom)%rhoijp(1:nselect,1:nspden)

!  Select non-zero values of rhoij(t+dt)
   if (hasmoved) then
    nselect=0
    do klmn=1,lmn2_size
     if (any(abs(rhoijtmp(klmn,:))>tol10)) then
      nselect=nselect+1
      pawrhoij(iatom)%rhoijselect(nselect)=klmn
      do ispden=1,nspden
       pawrhoij(iatom)%rhoijp(nselect,ispden)=rhoijtmp(klmn,ispden)
      end do
     end if
    end do
    pawrhoij(iatom)%nrhoijsel=nselect
    deallocate(rhoijtmp)
   end if

  end do !iatom
 end if !usepaw

!---------------------------------------------------------------
!----------- End
!---------------------------------------------------------------

!Rotate history indexes
 scf_history%hindex(1)=ind1new
 scf_history%hindex(2)=ind2new

end subroutine extraprho

!!***
