!{\src2tex{textfont=tt}}
!!****f* ABINIT/setvtr
!!
!! NAME
!! setvtr
!!
!! FUNCTION
!! Set up the trial potential and some energy terms
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, GMR, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(dtset%natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | ikhxc=exchange-correlation kernel treatment parameter
!!   |       if =0,1 no xc kernel, =2 spin-averaged (LDA) kernel
!!   | iprcch=govern the choice of preconditioner for the SCF cycle
!!   | iscf=determines the way the SCF cycle is handled
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!   | qprtrb(3)= integer wavevector of possible perturbing potential
!!   |            in basis of reciprocal lattice translations
!!   | typat(natom)=type integer for each atom in cell
!!   | vprtrb(2)=complex amplitude of possible perturbing potential; if nonzero,
!!   |  perturbing potential is added of the form
!!   |  V(G)=(vprtrb(1)+I*vprtrb(2))/2 at the values G=qprtrb and
!!   |  (vprtrb(1)-I*vprtrb(2))/2 at G=-qprtrb (integers)
!!   |  for each type of atom, from psp (used in norm-conserving only)
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2) (sphere for density and potential)
!!  istep=step number in the main loop of scfcv
!!  mgfft=maximum size of 1D FFTs
!!  moved_rhor=1 if the density was moved just before
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the array kxc
!!  ntypat=number of types of atoms in unit cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=>0 if some additional energies have to be computed
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=phase (structure factor) information.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=Fourier transform of electron density
!!  rhor(nfft,nspden)=electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol = unit cell volume (bohr^3)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)=previous reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  grewtn(3,natom)=grads of Ewald energy (hartree)
!!  kxc(nfft,nkxc)=exchange-correlation kernel, will be computed if nkxc/=0 .
!!                 see routine rhohxc for a more complete description
!!  strsxc(6)=xc contribution to stress tensor (hartree/bohr^3)
!!  vxcavg=mean of the vxc potential
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_xc=exchange-correlation energy (hartree)
!!  ==== if optene==2 or 4
!!   | e_localpsp=local psp energy (hartree)
!!  ==== if dtset%icoulomb == 0
!!   | e_ewald=Ewald energy (hartree)
!!  ==== if optene>=1
!!   | e_hartree=Hartree part of total energy (hartree units)
!!  ==== if optene==3 or 4
!!   | e_xcdc=exchange-correlation double-counting energy (hartree)
!!
!! SIDE EFFECTS
!!  initialized= if 0, the initialization of the gstate run is not yet finished
!!  moved_atm_inside=1 if the atomic positions were moved inside the SCF loop.
!!  vhartr(nfft)=Hartree potential (Hartree)
!!  vpsp(nfft)=local psp (Hartree)
!!  vtrial(nfft,nspden)= trial potential (Hartree)
!!  vxc(nfft,nspden)= xc potential (Hartree)
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!
!! NOTES
!!  In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfft,ngfft,mgfft) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      atm2fft,dotprod_vn,ewald,mkcore,mklocl,rhohxc,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setvtr(atindx1,dtset,energies,gmet,gprimd,&
&  grewtn,gsqcut,initialized,&
&  istep,kxc,mgfft,moved_atm_inside,moved_rhor,mpi_enreg,&
&  nattyp,nfft,ngfft,nhat,nhatgr,nhatgrdim,nkxc,ntypat,n1xccc,n3xccc,&
&  optene,pawtab,ph1d,psps,rhog,rhor,rmet,rprimd,strsxc,&
&  ucvol,usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,&
&  xccc3d,xred,xred_old)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12spacepar
 use interfaces_13xc
 use interfaces_14poisson
 use interfaces_15common, except_this_one => setvtr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mgfft,n1xccc,n3xccc,nfft,nhatgrdim,nkxc,ntypat
 integer,intent(in) :: optene,usexcnhat
 integer,intent(inout) :: initialized,moved_atm_inside,moved_rhor
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: nhat(nfft,dtset%nspden*psps%usepaw)
 real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim),rhog(2,nfft)
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3),xred_old(3,dtset%natom)
 real(dp),intent(inout) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),vhartr(nfft),vpsp(nfft)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(out) :: grewtn(3,dtset%natom),kxc(nfft,nkxc),strsxc(6)
 type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)

!Local variables-------------------------------
!scalars
 integer :: iatom,ifft,ikxc,indexat,ir,ispden,itypat,n1,n2,n3,nfftot,nn,offset
 integer :: optatm,optdyfr,optgr,option,optn,optn2,optstr,optv,rdwr
 real(dp) :: doti,ebb,ebn,ucvol_local
 character(len=500) :: message
!arrays
 real(dp),parameter :: identity(1:4)=(/1._dp,1._dp,0._dp,0._dp/)
 real(dp) :: dummy6(6),rhodum(1),tsec(2)
 real(dp),allocatable :: dummy(:),dyfr_dum(:,:,:),gr_dum(:,:),grtn(:,:)
 real(dp),allocatable :: rhojellg(:,:),rhojellr(:),vjell(:)

! *********************************************************************

!DEBUG
!write(6,*)' setvtr : enter '
!write(6,*)' n1xccc=',n1xccc
!write(6,*)' initialized=',initialized
!write(6,*)' moved_atm_inside=',moved_atm_inside
!write(6,*)' istep=',istep
!write(6,*)' iprcch=',iprcch
!write(6,*)' moved_rhor=',moved_rhor
!stop
!ENDDEBUG

 call timab(91,1,tsec)

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

!Get Ewald energy and Ewald forces
!--------------------------------------------------------------
 call timab(5,1,tsec)
 if ((dtset%icoulomb == 0) .or. (dtset%icoulomb == 2)) then
! Periodic system, need to compute energy and forces due to replica and
! to correct the shift in potential calculation.
  call ewald(energies%e_ewald,gmet,grewtn,dtset%natom,ntypat,rmet,dtset%typat,ucvol,xred,psps%ziontypat)
 else if (dtset%icoulomb == 1) then
! In a non periodic system (real space computation), the G=0 divergence
! doesn't occur and ewald is not needed. Only the ion/ion interaction
! energy is relevant and used as Ewald energy and gradient.
  call ionion_realSpace(dtset, energies%e_ewald, grewtn, rprimd, xred, psps%ziontypat)
 end if
 call timab(5,2,tsec)

!PAW: compute Vloc and core charge together in reciprocal space
!--------------------------------------------------------------
 if (psps%usepaw==1) then

  call timab(552,1,tsec)

  optatm=1;optdyfr=0;optgr=0;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1
  call atm2fft(atindx1,xccc3d,vpsp,dummy,dummy,energies%e_localpsp,dummy,gmet,gprimd,&
&  dummy,dummy,gsqcut,mgfft,mpi_enreg,psps%mqgrid_vl,&
&  dtset%natom,nattyp,nfft,ngfft,ntypat,&
&  optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&  pawtab,ph1d,psps%qgrid_vl,dtset%qprtrb,dummy,dummy6,dummy6,&
&  ucvol,psps%usepaw,dummy,dtset%vprtrb,psps%vlspl)
  call timab(552,2,tsec)

 else

! Norm-conserving: compute Vloc in reciprocal space
! and core charge in real space
! --------------------------------------------------------------

! Compute local ionic pseudopotential vpsp
  option=1
  allocate(dyfr_dum(3,3,dtset%natom),gr_dum(3,dtset%natom))  !Dummy variables
  call mklocl(dtset,dyfr_dum,energies%e_localpsp,gmet,gprimd,&
&  gr_dum,gsqcut,dummy6,mgfft,mpi_enreg,dtset%natom,nattyp,&
&  nfft,ngfft,dtset%nspden,ntypat,option,ph1d,psps,&
&  dtset%qprtrb,rhodum,rhor,rmet,rprimd,ucvol,dtset%vprtrb,vpsp,xred)

! Compute 3D core electron density xccc3d
  if (n1xccc/=0) then
   call timab(91,2,tsec)
   call timab(92,1,tsec)
   n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
   call mkcore(dummy6,dyfr_dum,gr_dum,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&   n1,n1xccc,n2,n3,option,rprimd,dtset%typat,ucvol,vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)
   call timab(92,2,tsec)
   call timab(91,1,tsec)
  end if
  deallocate(dyfr_dum,gr_dum)

 end if  ! PAW or NC

!Adds the jellium potential to the local part of ionic potential
 if (dtset%jellslab/=0) then
  allocate(vjell(nfft),rhojellg(2,nfft),rhojellr(nfft))
  option=1
  call jellium(gmet,gsqcut,mpi_enreg,nfft,ngfft,dtset%nspden,option,dtset%paral_kgb,&
&  dtset%slabwsrad,rhojellg,rhojellr,rprimd,vjell,dtset%slabzbeg,dtset%slabzend)
! Compute background-background energy
  call dotprod_vn(1,rhojellr,ebb,doti,mpi_enreg,nfft,nfftot,1,1,vjell,ucvol)
  ebb=half*ebb
! Compute electrostatic energy between background and nuclei before adding vjell to vpsp
  call dotprod_vn(1,rhojellr,ebn,doti,mpi_enreg,nfft,nfftot,1,1,vpsp,ucvol)
! Update e_ewald with ebb and ebn
  energies%e_ewald=energies%e_ewald+ebb+ebn
! Compute gradient of ebn wrt tn
  if (psps%usepaw==1) then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' setvtr : WARNING -',ch10,&
&   '  The computation of forces due to jellium background',ch10,&
&   '  has to be verified in the PAW formalism.'
   call wrtout(6,message,'COLL')
   allocate(grtn(3,dtset%natom))
   optatm=0;optdyfr=0;optgr=1;optstr=0;optv=1;optn=0;optn2=1
   call atm2fft(atindx1,dummy,vpsp,dummy,dummy,energies%e_localpsp,dummy,gmet,gprimd,&
&   dummy,grtn,gsqcut,mgfft,mpi_enreg,psps%mqgrid_vl,&
&   dtset%natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
&   pawtab,ph1d,psps%qgrid_vl,dtset%qprtrb,rhojellg,dummy6,dummy6,&
&   ucvol,psps%usepaw,dummy,dtset%vprtrb,psps%vlspl)
!  Update grewtn with gradient of ebn wrt tn
   do iatom=1,dtset%natom
    grewtn(1:3,iatom)=grewtn(1:3,iatom)+grtn(1:3,iatom)
   end do
   deallocate(grtn)
  else ! of usepaw==1
   option=2
   allocate(dyfr_dum(3,3,dtset%natom),grtn(3,dtset%natom))  !Dummy variables
   call mklocl(dtset,dyfr_dum,energies%e_localpsp,gmet,gprimd,&
&   grtn,gsqcut,dummy6,mgfft,mpi_enreg,dtset%natom,nattyp,&
&   nfft,ngfft,1,ntypat,option,ph1d,psps,dtset%qprtrb,rhojellg,&
&   rhojellr,rmet,rprimd,ucvol,dtset%vprtrb,vpsp,xred)
!  Update grewtn with gradient of ebn wrt tn (reestablish order of atoms)
   do iatom=1,dtset%natom
    grewtn(1:3,atindx1(iatom))=grewtn(1:3,atindx1(iatom))+grtn(1:3,iatom)
   end do
   deallocate(dyfr_dum,grtn)
  end if ! of usepaw==1
  vpsp(:)=vpsp(:)+vjell(:)
  deallocate(vjell,rhojellg,rhojellr)
 end if

!If we are at the initialisation, or
!if the atom positions has changed and the non-linear core correction
!is included, or the rhor has changed, one needs to compute the xc stuff.
!One needs also to compute the Hartree stuff if the density changed,
!or at initialisation.
!--------------------------------------------------------------

!DEBUG
!write(6,*)' setvtr : istep,n1xccc,moved_rhor=',istep,n1xccc,moved_rhor
!ENDDEBUG

 if(istep==1 .or. n1xccc/=0 .or. moved_rhor==1) then

  option=0
  if(istep==1 .or. moved_rhor==1) option=1
  if (nkxc>0) option=2
  if (dtset%xclevel==2.and.(nkxc==3-2*mod(dtset%nspden,2))) option=12
  if(dtset%iscf==-1) option=-2
  if(dtset%ikhxc==2) then
   write(6,*)' %setvtr: CALL rhohxc with option=2'         !MF!DEBUGLINE
   write(6,*)' %        computing kxc = (kxc++ + kxc+-)/2' !MF!DEBUGLINE
  end if

  if (dtset%icoulomb == 0) then
!  Use the periodic solver to compute Hxc.
   call rhohxc(dtset,energies%e_xc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfft,ngfft,&
&   nhat,psps%usepaw,nhatgr,nhatgrdim,nkxc,dtset%nspden,n3xccc,&
&   option,rhog,rhor,rprimd,strsxc,&
&   usexcnhat,vhartr,vxc,vxcavg,xccc3d)

!  DEBUG
!  write(6,*)' setvtr : computed rhohxc'
!  ENDDEBUG
  else
!  Use the free boundary solver.
   call PSolver_rhohxc(dtset, energies%e_hartree, energies%e_xc, energies%e_vxc, &
&   mpi_enreg, rhor, rprimd, vhartr, vxc, vxcavg)
  end if
 end if

!Compute the trial potential
!-------------------------------------------------------------
 if (dtset%usewvl == 0) then
! Now, compute trial Hxc potential. Local psp potential will be added later.
  if(moved_atm_inside==0 .or.dtset%iscf>=10) then

!  Compute starting Hxc potential.
!  Multiply by identity, should not change anything if nspden /= 4
   do ispden=1,dtset%nspden
    vtrial(:,ispden)=vhartr(:)*identity(ispden)+vxc(:,ispden)
   end do

  else

!  One should be here only when moved_atm_inside==1
!  The (H)xc now added corrects the previous one.
   if(dtset%iprcch==1)then
!   xc was substracted off. This should be rationalized later
    do ispden=1,dtset%nspden
     vtrial(:,ispden)=vtrial(:,ispden)+vxc(:,ispden)
    end do
   else if(abs(dtset%iprcch)==2 .or.abs(dtset%iprcch)==5.or.abs(dtset%iprcch)==6)then
!   Hxc was substracted off. This should be rationalized later
    do ispden=1,dtset%nspden
     vtrial(:,ispden)=vtrial(:,ispden)+vhartr(:)*identity(ispden)+vxc(:,ispden)
    end do
   end if
  end if

! Adds the local part of the potential
  if ((moved_atm_inside==0).or.(dtset%iprcch/=3)) then
   do ispden=1,min(2,dtset%nspden)
    do ifft=1,nfft
     vtrial(ifft,ispden)=vtrial(ifft,ispden)+vpsp(ifft)
    end do
   end do
  end if
 else

  call wvl_newvtr(dtset, mpi_enreg, nn, offset, vhartr, vpsp, vtrial, vxc)

 end if

!Compute parts of total energy depending on potentials
!--------------------------------------------------------------
 if (dtset%usewvl == 0) then
  ucvol_local = ucvol
 else
! We need to tune the volume when wavelets are used because, not
! all FFT points are used.
  ucvol_local = (half * dtset%wvl_hgrid) ** 3 * nfftot
 end if

!DEBUG
!write(6,*)' setvtr : will compute ehart, optene,dtset%icoulomb=',optene,dtset%icoulomb
!ENDDEBUG

 if (optene>=1 .and. dtset%icoulomb == 0) then

! Compute Hartree energy ehart
! Already available in the Psolver case through psolver_rhohxc().
  call dotprod_vn(1,rhor,energies%e_hartree,doti,mpi_enreg,nfft,nfftot,1,1,vhartr,ucvol_local)
  energies%e_hartree = half * energies%e_hartree

 end if

 if (optene==2.or.optene==4) then

! Compute local psp energy eei
  call dotprod_vn(1,rhor,energies%e_localpsp,doti,mpi_enreg,nfft,nfftot,1,1,vpsp,ucvol_local)

 end if

 if (optene==3.or.optene==4) then

! Compute double-counting XC energy enxcdc
  call dotprod_vn(1,rhor,energies%e_xcdc,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local)

 end if

!--------------------------------------------------------------

!The initialisation for the new atomic positions has been done
 moved_atm_inside=0

 call timab(91,2,tsec)

!DEBUG
!write(6,*)' setvtr : exit '
!stop
!ENDDEBUG

end subroutine setvtr
!!***
