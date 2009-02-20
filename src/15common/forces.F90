!{\src2tex{textfont=tt}}
!!****f* ABINIT/forces
!!
!! NAME
!! forces
!!
!! FUNCTION
!! Assemble gradients of various total energy terms with respect
!! to reduced coordinates, including possible symmetrization,
!! in order to produce forces.
!!     fcart(i,iat) = d(Etot)/(d(r(i,iat)))
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                               - \Omega E.P term to the total energy
!!   |          /= 4: electric field is off
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along specified direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | iprcch=governs the mixed electronic-atomic part of the preconditioner
!!   | natom=number of atoms in cell
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetries in space group
!!   | prtvol=integer controlling volume of printed output
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  grnl(3*natom)=gradients of Etot due to nonlocal contributions
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  ntypat=number of types of atoms
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) array
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  vresid(nfft,nspden)=potential residual (if non-collinear magn., only trace of it)
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real space
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)=previous reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  diffor=maximal absolute value of changes in the components of
!!         force between the input and the output.
!!  favg(3)=mean of the forces before correction for translational symmetry
!!  forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!  fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!  gresid(3,natom)=forces due to the residual of the density/potential
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  grxc(9+3*natom)=d(Exc)/d(xred) if core charges are used
!!  maxfor=maximal absolute value of the output array force.
!!  synlgr(3,natom)=symmetrized d(enl)/d(xred)
!!
!! SIDE EFFECTS
!!  fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!    Note : unlike fred, this array has been corrected by enforcing
!!    the translational symmetry, namely that the sum of force
!!    on all atoms is zero.
!!
!! NOTES
!! * Symmetrization of gradients with respect to reduced
!!   coordinates xred is conducted according to the expression
!!   [d(e)/d(t(n,a))]_symmetrized = (1/Nsym) Sum(S) symrec(n,m,S)*
!!                [d(e)/d(t(m,b))]_unsymmetrized
!!   where t(m,b)= (symrel^-1)(m,n)*(t(n,a)-tnons(n)) and tnons
!!   is a possible nonsymmorphic translation.  The label "b" here
!!   refers to the atom which gets rotated into "a" under symmetry "S".
!!   symrel is the symmetry matrix in real space, which is the inverse
!!   transpose of symrec.  symrec is the symmetry matrix in reciprocal
!!   space.  sym_cartesian = R * symrel * R^-1 = G * symrec * G^-1
!!   where the columns of R and G are the dimensional primitive translations
!!   in real and reciprocal space respectively.
!! * Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      etotfor,forstr
!!
!! CHILDREN
!!      atm2fftconstrf,fourdp,fresid,metric,mkcore,mklocl,sygrad,timab,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&
&                  grhf,grnl,grxc,gsqcut,indsym,kxc,&
&                  maxfor,mgfft,mpi_enreg,n1xccc,n3xccc,&
&                  nattyp,nfft,ngfft,nkxc,ntypat,&
&                  pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,&
&                  vresid,vxc,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_13xc
 use interfaces_15common, except_this_one => forces
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,n1xccc,n3xccc,nfft,nkxc,ntypat
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: diffor,maxfor
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: grewtn(3,dtset%natom),grnl(3*dtset%natom)
 real(dp),intent(in) :: kxc(nfft,nkxc),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: fcart(3,dtset%natom),forold(3,dtset%natom)
 real(dp),intent(inout) :: vresid(nfft,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fred(3,dtset%natom),gresid(3,dtset%natom)
 real(dp),intent(out) :: grhf(3,dtset%natom),grxc(3,dtset%natom)
 real(dp),intent(out) :: synlgr(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: iatom,iattyp,indx,ispden,itypat,mu,optatm,optdyfr,optgr,option,optn
 integer :: optn2,optstr,optv
 real(dp) :: eei_dum,ucvol
!arrays
 integer :: qprtrb_dum(3)
 real(dp) :: dummy6(6),fioncart(3),gmet(3,3),gprimd(3,3),rmet(3,3),tsec(2)
 real(dp) :: vprtrb_dum(2)
 real(dp),allocatable :: dummy(:),dyfrlo_dum(:,:,:),dyfrx2_dum(:,:,:),fin(:,:)
 real(dp),allocatable :: fionred(:,:),grl(:,:),grnl_tmp(:,:),grtn(:,:)
 real(dp),allocatable :: grtn_indx(:,:),grxc_indx(:,:),v_dum(:),vxctotg(:,:)
 real(dp),allocatable :: xccc3d_dum(:)

! *************************************************************************

!DEBUG
!write(6,*)' forces: enter '
!stop
!ENDDEBUG

 call timab(69,1,tsec)

!Save input value of forces
 allocate(fin(3,dtset%natom));fin(:,:)=fcart(:,:)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!=======================================================================
!========= Local pseudopotential and core charge contributions =========
!=======================================================================

 allocate(grl(3,dtset%natom))

!PAW: compute local psp and core charge contribs together
!in reciprocal space
!-----------------------------------------------------------------------
 if (psps%usepaw==1) then

  call timab(550,1,tsec)
  if (n3xccc>0) then
   allocate(v_dum(nfft),vxctotg(2,nfft))
   v_dum(:)=vxc(:,1);if (dtset%nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
   call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
   call zerosym(vxctotg,2,mpi_enreg,ngfft(1),ngfft(2),ngfft(3))
   deallocate(v_dum)
  end if
  optatm=0;optdyfr=0;optgr=1;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1
  call atm2fft(atindx1,dummy,dummy,dummy,dummy,eei_dum,dummy,gmet,gprimd,&
&  grxc,grl,gsqcut,mgfft,mpi_enreg,psps%mqgrid_vl,&
&  dtset%natom,nattyp,nfft,ngfft,ntypat,&
&  optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&  pawtab,ph1d,psps%qgrid_vl,qprtrb_dum,rhog,dummy6,dummy6,&
&  ucvol,psps%usepaw,vxctotg,vprtrb_dum,psps%vlspl)
  if (n3xccc>0) deallocate(vxctotg)
  if (n3xccc==0) grxc=zero
  call timab(550,2,tsec)
 else

! Norm-conserving: compute local psp contribution in reciprocal space
! and core charge contribution in real space
! -----------------------------------------------------------------------
  option=2
  allocate(dyfrlo_dum(3,3,dtset%natom),grtn_indx(3,dtset%natom),v_dum(nfft))
  call mklocl(dtset,dyfrlo_dum,eei_dum,gmet,gprimd,grtn_indx,gsqcut,dummy6,mgfft,&
&  mpi_enreg,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,ntypat,option,ph1d,psps,&
&  qprtrb_dum,rhog,rhor,rmet,rprimd,ucvol,vprtrb_dum,v_dum,xred)

  do iatom=1,dtset%natom
!  Has to use the indexing array atindx1
   grl(1:3,atindx1(iatom))=grtn_indx(1:3,iatom)
  end do
  deallocate(dyfrlo_dum,grtn_indx,v_dum)
! If gradients are computed in real space, we need to symetrise
! the system before summing.
! Rshaltaf: I changed the following line to include surfaces BC
  if (dtset%icoulomb == 1 .or. dtset%icoulomb == 2) then
   allocate(grnl_tmp(3,dtset%natom))
   call sygrad(grnl_tmp,dtset%natom,grl,dtset%nsym,symrec,indsym)
   grl(:, :) = grnl_tmp(:, :)
   deallocate(grnl_tmp)
  end if

  if (n3xccc>0) then
   call timab(53,1,tsec)
   allocate(dyfrx2_dum(3,3,dtset%natom),xccc3d_dum(n3xccc))
   call mkcore(dummy6,dyfrx2_dum,grxc,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,ngfft(1),n1xccc,ngfft(2),&
&   ngfft(3),option,rprimd,dtset%typat,ucvol,vxc,psps%xcccrc,psps%xccc1d,xccc3d_dum,xred)
   deallocate(dyfrx2_dum,xccc3d_dum)
   call timab(53,2,tsec)
  else
   grxc(:,:)=zero
  end if
 end if

!=======================================================================
!===================== Nonlocal contributions ==========================
!=======================================================================

!Only has to apply symmetries
 allocate(grnl_tmp(3,dtset%natom))
 do iatom=1,dtset%natom
  indx=3*(iatom-1);grnl_tmp(1:3,atindx1(iatom))=grnl(indx+1:indx+3)
 end do
 if (dtset%usewvl == 0) then
  call sygrad(synlgr,dtset%natom,grnl_tmp,dtset%nsym,symrec,indsym)
 else
  synlgr = grnl_tmp
 end if
 deallocate(grnl_tmp)

!=======================================================================
!============ Density/potential residual contributions =================
!=======================================================================

 if (dtset%usewvl==0.and.abs(dtset%iprcch)>=1.and.abs(dtset%iprcch)<=3) then
  call fresid(dtset,gmet,gresid,gsqcut,mpi_enreg,nfft,ngfft,ntypat,1,&
&  pawtab,rhor,rprimd,ucvol,vresid,xred,xred,psps%znuclpsp)
 else if (dtset%usewvl==0.and.(abs(dtset%iprcch)==4.or.abs(dtset%iprcch)==6)) then
  call fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,&
&  mpi_enreg,psps%mqgrid_vl,nattyp,nfft,ngfft,ntypat,pawtab,ph1d,&
&  psps%qgrid_vl,ucvol,psps%usepaw,vresid,psps%zionpsp,psps%znuclpsp)
 else
  gresid(:,:)=zero
 end if

!=======================================================================
!======================= Other contributions ===========================
!=======================================================================

!Ewald energy contribution to forces as already been computed in "ewald"

!Potential residual contribution to forces as already been computed (in vresfo or forstr)

!Add Berry phase contributions (berryopt == 4)
!(compute the electric field force on the ion cores)
 if (dtset%berryopt==4) then
  allocate(fionred(3,dtset%natom));fionred(:,:)=zero
  iatom = 0
  do itypat=1,ntypat
   do iattyp=1,nattyp(itypat)
    iatom=iatom+1
    fioncart(:)=psps%ziontypat(itypat)*dtset%efield(:)
    do mu=1,3
     fionred(mu,iatom)=rprimd(1,mu)*fioncart(1) &
&     +rprimd(2,mu)*fioncart(2) &
&     +rprimd(3,mu)*fioncart(3)
    end do
   end do
  end do
 end if

!=======================================================================
!======= Assemble the various contributions to the forces ==============
!=======================================================================

!write(*,*) "#grl", grl(:,:)
!write(*,*) "#gnl", grnl(:)
!write(*,*) "#gre", grewtn(:,:)
!write(*,*) "#syn", synlgr(:,:)
!write(*,*) "#grx", grxc(:,:)
!write(*,*) "#res", gresid(:,:)

!Collect grads of etot wrt reduced coordinates
!This gives non-symmetrized Hellman-Feynman reduced gradients
 allocate(grtn(3,dtset%natom))
 grtn(:,:)=grl(:,:)+grewtn(:,:)+synlgr(:,:)+grxc(:,:)
 if (dtset%berryopt==4) grtn(:,:)=grtn(:,:)-fionred(:,:)

!write(*,*) "####### Gradients before sym #########"
!write(*,*) "#grt", grtn

!Symmetrize explicitly for given space group and store in grhf :
 call sygrad(grhf,dtset%natom,grtn,dtset%nsym,symrec,indsym)

!Add residual potential correction
 grtn(:,:)=grtn(:,:)+gresid(:,:)

!Symmetrize all grads explicitly for given space group:
 if (dtset%usewvl == 0) then
  call sygrad(fred,dtset%natom,grtn,dtset%nsym,symrec,indsym)
 else
  fred = grtn
 end if
!Note conversion to cartesian coordinates (bohr) AND
!negation to make a force out of a gradient
 favg(:)=zero
 do iatom=1,dtset%natom
  do mu=1,3
   fcart(mu,iatom)= - (gprimd(mu,1)*fred(1,iatom)+&
&   gprimd(mu,2)*fred(2,iatom)+&
&   gprimd(mu,3)*fred(3,iatom))
   favg(mu)=favg(mu)+fcart(mu,iatom)
  end do
 end do

!Subtract off average force from each force component
!to avoid spurious drifting of atoms across cell.
 favg(:)=favg(:)/dble(dtset%natom)
 if(dtset%jellslab/=0) favg(3)=zero
 do iatom=1,dtset%natom
  fcart(:,iatom)=fcart(:,iatom)-favg(:)
 end do

!write(*,*) "#fca", fcart

!Compute maximal force and maximal difference
 maxfor=zero;diffor=zero
 do iatom=1,dtset%natom
  do mu=1,3
   if (dtset%iatfix(mu,iatom) /= 1) then
    maxfor=max(maxfor,abs(fcart(mu,iatom)))
    diffor=max(diffor,abs(fcart(mu,iatom)-fin(mu,iatom)))
   else if (dtset%ionmov==4 .or. dtset%ionmov==5) then
!   Make the force vanish on fixed atoms when ionmov=4 or 5
!   This is because fixing of atom cannot be imposed at the
!   level of a routine similar to brdmin or moldyn for these options.
    fcart(mu,iatom)=zero
   end if
  end do
 end do

!Apply any generalized constraints to the forces
 if (dtset%nconeq>0) call constrf(diffor,fcart,forold,fred,dtset%iatfix,dtset%ionmov,maxfor,&
& dtset%natom,dtset%nconeq,dtset%prtvol,rprimd,dtset%wtatcon,xred)

!=======================================================================
!Memory deallocations
 deallocate(grl,grtn,fin)
 if (dtset%berryopt==4) deallocate(fionred)

 call timab(69,2,tsec)

!DEBUG
!write(6,*)' forces: exit '
!stop
!ENDDEBUG

end subroutine forces
!!***
