!{\src2tex{textfont=tt}}
!!****f* ABINIT/etotfor
!! NAME
!! etotfor
!!
!! FUNCTION
!! This routine is called to compute the total energy and various parts of it.
!! The routine computes -if requested- the forces.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                - \Omega E.P term to the total energy
!!   |          /= 4: electric field is off
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along some direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | iprcch=governs the mixed electronic-atomic part of the preconditioner
!!   | natom=number of atoms in cell.
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetry elements in space group
!!   | occopt=option for occupancies
!!   | prtvol=integer controlling volume of printed output
!!   | tsmear=smearing energy or temperature (if metal)
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  efield_dot = reciprocal lattice coordinates of the electric field
!!  grewtn(3,natom)=grads of Ewald energy (hartree)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  ntypat=number of types of atoms in unit cell.
!!  nvresid(nfft,nspden)=potential or density residual
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=option for the computation of total energy
!!         (-1=no computation; 0=direct scheme; 1=double-counting scheme)
!!  optforces=option for the computation of forces
!!  optres=0 if residual array (nvresid) contains the potential residual
!!        =1 if residual array (nvresid) contains the density residual
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pel = reduced coordinates of the electronic polarization
!!        (pel does not take into account the factor 1/ucvol)
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) information.
!!  pion = reduced coordinates of the ionic polarization
!!        (pel does not take into account the factor 1/ucvol)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  ucvol = unit cell volume (Bohr**3)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vhartr(nfft)=array for holding Hartree potential
!!  vpsp(nfft)=array for holding local psp
!!  vxc(nfft,nspden)=array for holding XC potential
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  deltae=change in total energy
!!         between the previous and present SCF cycle
!!  etotal=total energy (hartree)
!!  ===== if optforces==1
!!   diffor=maximum absolute change in component of forces between present and previous SCF cycle.
!!   favg(3)=mean of fcart before correction for translational symmetry
!!   fcart(3,natom)=cartesian forces from fred (hartree/bohr)
!!   fred(3,natom)=symmetrized form of grtn (grads of Etot) (hartree)
!!   gresid(3,natom)=forces due to the residual of the density/potential
!!   grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!   grxc(3,natom)=d(Exc)/d(xred) derivatives (0 without core charges)
!!   maxfor=maximum absolute value of force
!!   synlgr(3,natom)=symmetrized form of grads of Enl (hartree)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  elast=previous value of the energy,
!!        needed to compute deltae, then updated.
!!  energies <type(energies_type)>=all part of total energy.
!!   | entropy(IN)=entropy due to the occupation number smearing (if metal)
!!   | e_localpsp(IN)=local psp energy (hartree)
!!   | e_eigenvalues(IN)=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_ewald(IN)=Ewald energy (hartree)
!!   | e_hartree(IN)=Hartree part of total energy (hartree units)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_kinetic(IN)=kinetic energy part of total energy.
!!   | e_nonlocalpsp(IN)=nonlocal pseudopotential part of total energy.
!!   | e_xc(IN)=exchange-correlation energy (hartree)
!!   | e_xcdc(IN)=exchange-correlation double-counting energy (hartree)
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_elecfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the electric field:  enefield = -ucvol*E*P
!!   | e_entropy(OUT)=entropy energy due to the occupation number smearing (if metal)
!!   |                this value is %entropy * dtset%tsmear (hartree).
!!  ===== if optforces==1
!!   forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!   grnl(3*natom)=gradients of Etot due to nonlocal contributions
!!                 Input for norm-conserving psps, output for PAW
!!  ===== if psps%usepaw==1
!!   pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
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
!!      forces,pawgrnl,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine etotfor(atindx1,deltae,diffor,dtset,efield_dot,elast,energies,&
&  etotal,favg,fcart,forold,fred,gresid,grewtn,grhf,grnl,&
&  grxc,gsqcut,indsym,kxc,maxfor,mgfft,mpi_enreg,nattyp,&
&  nfft,ngfft,nhat,nkxc,ntypat,nvresid,n1xccc,n3xccc,optene,optforces,optres,&
&  pawang,pawfgrtab,pawrhoij,pawtab,pel,ph1d,pion,psps,rhog,rhor,rprimd,symrec,synlgr,&
&  ucvol,usepaw,usexcnhat,vhartr,vpsp,vxc,xccc3d,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_13paw
 use interfaces_15common, except_this_one => etotfor
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,n1xccc,n3xccc,nfft,nkxc,ntypat,optene,optforces
 integer,intent(in) :: optres,usepaw,usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(inout) :: elast
 real(dp),intent(out) :: deltae,diffor,etotal,maxfor
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: efield_dot(3),grewtn(3,dtset%natom),kxc(nfft,nkxc)
 real(dp),intent(in) :: pel(3),ph1d(2,3*(2*mgfft+1)*dtset%natom),pion(3)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: vhartr(nfft),vpsp(nfft),vxc(nfft,dtset%nspden)
 real(dp),intent(in) :: xccc3d(n3xccc)
 real(dp),intent(inout) :: forold(3,dtset%natom),grnl(3*dtset%natom)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*psps%usepaw)
 real(dp),intent(inout) :: nvresid(nfft,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fcart(3,dtset%natom),fred(3,dtset%natom)
 real(dp),intent(out) :: gresid(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(out) :: grxc(3,dtset%natom),synlgr(3,dtset%natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(dtset%natom*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(dtset%natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: dimnhat,iatm,iatom,ifft,ispden,optgr,optgr2,option,optnc,optstr
 real(dp) :: epel,epion
 character(len=500) :: message
!arrays
 real(dp) :: str_dum(6),tsec(2)
 real(dp),allocatable :: dummy(:),nhat_dum(:,:),work(:,:)

! *********************************************************************

!DEBUG
!write(6,*)' etotfor : enter'
!ENDDEBUG

 call timab(80,1,tsec)

 if (optene>-1) then

! When the finite-temperature VG broadening scheme is used,
! the total entropy contribution "tsmear*entropy" has a meaning,
! and gather the two last terms of Eq.8 of the VG paper
! Warning : might have to be changed for fixed moment calculations
  if(dtset%occopt>=3 .and. dtset%occopt<=7) then
   energies%e_entropy = - dtset%tsmear * energies%entropy
  else
   energies%e_entropy = zero
  end if

! Turn it into an electric enthalpy, by adding both ionic and electronic contributions
  energies%e_elecfield = zero
  if (dtset%berryopt == 4) then
!  First, ionic contribution (epion):
   epion = -1_dp*(efield_dot(1)*pion(1) + efield_dot(2)*pion(2) + &
&   efield_dot(3)*pion(3))
   energies%e_elecfield = energies%e_elecfield + epion
!  Now, electronic contribution (epel):
   epel = -1_dp*(efield_dot(1)*pel(1) + efield_dot(2)*pel(2) + &
&   efield_dot(3)*pel(3))
   energies%e_elecfield = energies%e_elecfield + epel
  end if

! Compute total (free) energy by direct scheme
  if (optene==0) then
   etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&   energies%e_localpsp + energies%e_corepsp + &
&   energies%e_entropy + energies%e_elecfield
   etotal = etotal + energies%e_ewald
   if (usepaw==0) then
    etotal = etotal + energies%e_nonlocalpsp
   else
    etotal = etotal + energies%e_paw
   end if
   if (dtset%usewvl == 1) then
    etotal = etotal - energies%e_vxc
   end if
  end if

! Compute total (free) energy by double-counting scheme
  if (optene==1) then
   etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
&   energies%e_xcdc + energies%e_corepsp + &
&   energies%e_entropy + energies%e_elecfield
   etotal = etotal + energies%e_ewald
   if (usepaw/=0) then
    etotal = etotal + energies%e_pawdc
   end if
  end if

! Compute energy residual
  deltae=etotal-elast
  elast=etotal
! DEBUG
! write(6,*) 'eeig-ehart+enxc-enxcdc+eew+eii+eent+enefield+epawdc',eeig,ehart,enxc,enxcdc,eew,eii,eent,enefield,epawdc
! ENDDEBUG

 end if !optene/=-1

 call timab(80,2,tsec)

!------Compute forces-----------------------------------------------------

 if (optforces==1) then

! PAW: add gradients due to Dij derivatives to non-local term
  if (usepaw==1) then
   allocate(work(nfft,dtset%nspden))
   do ispden=1,min(dtset%nspden,2)
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(ispden,nfft,vhartr,work,vpsp,vxc)
    do ifft=1,nfft
     work(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
    end do
!   $OMP END PARALLEL DO
   end do
   if(dtset%nspden==4)then
    do ispden=3,4
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(ispden,nfft,work,vxc)
    do ifft=1,nfft
     work(ifft,ispden)=vxc(ifft,ispden)
    end do
!   $OMP END PARALLEL DO
    end do
   end if
!  DEBUG
!  write(6,*) 'in etotfor call pawgrnl',mpi_enreg%paral_level
!  ENDDEBUG
   dimnhat=0;optgr=1;optgr2=0;optstr=0;allocate(nhat_dum(1,0))
   call pawgrnl(atindx1,dimnhat,dtset%nspden,dummy,grnl,mpi_enreg,dtset%natom,nattyp,&
&   nfft,ngfft,nhat_dum,dummy,dtset%nspden,dtset%nsym,ntypat,optgr,optgr2,optstr,&
&   pawang,pawfgrtab,pawrhoij,pawtab,rprimd,symrec,dtset%typat,work)
   deallocate(nhat_dum,work)
  end if

! If residual is a density residual (and forces from residual asked),
! has to convert it into a potential residualbefore calling forces routine
  if (optres==1 .and. dtset%usewvl==0.and.abs(dtset%iprcch)>=1 .and. &
&  abs(dtset%iprcch)<=6.and.abs(dtset%iprcch)/=5) then
   option=0; if (dtset%iprcch<0) option=1
   allocate(work(nfft,dtset%nspden))
   optnc=1;if (dtset%nspden==4.and.(abs(dtset%iprcch)==4.or.abs(dtset%iprcch)==6)) optnc=2
   call nres2vres(dtset,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,nhat,&
&   nkxc,nvresid,n3xccc,optnc,option,pawang,pawfgrtab,pawrhoij,pawtab,&
&   rhog,rhor,rprimd,usepaw,usexcnhat,work,xccc3d)
   call forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&
&   grhf,grnl,grxc,gsqcut,indsym,kxc,maxfor,mgfft,mpi_enreg,&
&   n1xccc,n3xccc,nattyp,nfft,ngfft,nkxc,ntypat,&
&   pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,work,vxc,xred)
   deallocate(work)
  else
   call forces(atindx1,diffor,dtset,favg,fcart,forold,fred,gresid,grewtn,&
&   grhf,grnl,grxc,gsqcut,indsym,kxc,maxfor,mgfft,mpi_enreg,&
&   n1xccc,n3xccc,nattyp,nfft,ngfft,nkxc,ntypat,&
&   pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,nvresid,vxc,xred)
  end if

! Returned fred are full symmetrized gradients of Etotal
! wrt reduced coordinates xred, d(Etotal)/d(xred)
! Forces are contained in array fcart

 else   ! if optforces==0
  fcart=zero
  fred=zero
  favg=zero
  diffor=zero
  gresid=zero
  grhf=zero
  maxfor=zero
  synlgr=zero
 end if

 call timab(80,2,tsec)

!DEBUG
!write(6,*)' etotfor : exit '
!write(6,*) etotal
!call leave_new("COLL")
!stop
!ENDDEBUG

end subroutine etotfor
!!***
