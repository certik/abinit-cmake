!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyfro3
!! NAME
!! dyfro3
!!
!! FUNCTION
!! Compute the different parts of the frozen-wavefunction part of
!! the dynamical matrix, except the non-local one, computed previously.
!! Also (when installed) symmetrize the different part and their sum.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (XG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!     (bohr**-1)
!!  gsqcut=cutoff on G^2 based on ecut
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!   under symmetry operations (computed in symatm)
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=dimensioned number of q grid points for local psp spline
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  qgrid(mqgrid)=q point array for local psp spline fits
!!  rhog(2,nfft)=electron density in G space
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real
!!   space--only used when n1xccc/=0
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  dyfrlo(3,3,natom)=frozen wavefunctions part of the dynamical matrix
!!                    (local only)
!!  dyfrwf(3,3,natom)=frozen wavefunctions part of the dynamical matrix
!!                    (local + non-local)
!!  dyfrxc(3,3,natom)=frozen wavefunctions part of the dynamical matrix
!!                    (non-linear xc core correction)
!!
!! SIDE EFFECTS
!! Input/Output
!!  dyfrnl(3,3,natom)=frozen wavefunctions part of the dynamical matrix
!!                    (non-local only)
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      atm2fft,mkcore,mklocl_recipspace,sydy3
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dyfro3(atindx1,dyfrnl,dyfrlo,dyfrwf,dyfrxc,&
&  gmet,gprimd,gsqcut,indsym,mgfft,mpi_enreg,mqgrid,natom,nattyp,&
&  nfft,ngfft,nspden,nsym,ntypat,n1xccc,n3xccc,paral_kgb,pawtab,ph1d,qgrid,&
&  rhog,rprimd,symrec,typat,ucvol,usepaw,vlspl,vxc,xcccrc,xccc1d,xccc3d,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_12ffts
 use interfaces_13xc
 use interfaces_15common
 use interfaces_16response, except_this_one => dyfro3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,n1xccc,n3xccc,natom,nfft,nspden,nsym,ntypat
 integer,intent(in) :: paral_kgb,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indsym(4,nsym,natom),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),rhog(2,nfft),rprimd(3,3)
 real(dp),intent(in) :: vlspl(mqgrid,2,ntypat),vxc(nfft,nspden)
 real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat),xcccrc(ntypat),xred(3,natom)
 real(dp),intent(inout) :: dyfrnl(3,3,natom),xccc3d(n3xccc)
 real(dp),intent(out) :: dyfrlo(3,3,natom),dyfrwf(3,3,natom),dyfrxc(3,3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu,n1,n2,n3,nu,optatm,optdyfr,optgr,option,optn,optn2,optstr
 integer :: optv
 real(dp) :: eei
!arrays
 integer :: qprtrb(3)
 real(dp) :: dummy6(6),tsec(2),vprtrb(2)
 real(dp),allocatable :: dummy(:),dyfrlo_indx(:,:,:),gr_dum(:,:),v_dum(:)
 real(dp),allocatable :: vxctotg(:,:)

! *************************************************************************

 if(nspden==4)then
  write(6,*)' dyfro3 : does not work yet for nspden=4'
  stop
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

 if (usepaw==1) then

! PAW: compute local psp and core charge contribs together
! in reciprocal space
! -----------------------------------------------------------------------
  call timab(563,1,tsec)
  if (n3xccc>0) then
   allocate(v_dum(nfft),vxctotg(2,nfft))
   v_dum(:)=vxc(:,1);if (nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
   call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   call zerosym(vxctotg,2,mpi_enreg,ngfft(1),ngfft(2),ngfft(3))
   deallocate(v_dum)
  end if
  optatm=0;optdyfr=1;optgr=0;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1
  call atm2fft(atindx1,dummy,dummy,dyfrxc,dyfrlo,eei,dummy,gmet,gprimd,dummy,dummy,&
&  gsqcut,mgfft,mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&
&  optatm,optdyfr,optgr,optn,optn2,optstr,optv,paral_kgb,&
&  pawtab,ph1d,qgrid,qprtrb,rhog,dummy6,dummy6,ucvol,usepaw,vxctotg,vprtrb,vlspl)
  if (n3xccc>0) deallocate(vxctotg)
  if (n3xccc==0) dyfrxc=zero
 else

! Norm-conserving: compute local psp contribution in reciprocal space
! and core charge contribution in real space
! -----------------------------------------------------------------------
  option=4
  allocate(dyfrlo_indx(3,3,natom),gr_dum(3,natom),v_dum(nfft))  !Dummy arrays
  call mklocl_recipspace(dyfrlo_indx,eei,gmet,gprimd,&
&  gr_dum,gsqcut,dummy6,mgfft,mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
&  ntypat,option,paral_kgb,ph1d,qgrid,qprtrb,rhog,ucvol,vlspl,vprtrb,v_dum)
  do iatom=1,natom
!  Reestablish correct order of atoms
   dyfrlo(1:3,1:3,atindx1(iatom))=dyfrlo_indx(1:3,1:3,iatom)
  end do
  deallocate(dyfrlo_indx,v_dum)
  if(n1xccc/=0)then
   call mkcore(dummy6,dyfrxc,gr_dum,mpi_enreg,natom,nfft,nspden,ntypat,&
&   n1,n1xccc,n2,n3,option,rprimd,typat,ucvol,vxc,xcccrc,xccc1d,xccc3d,xred)
  end if
  deallocate(gr_dum)
 end if

!Symmetrize dynamical matrix explicitly for given space group:
!Uses dyfrwf as work space
 dyfrwf(:,:,:)=dyfrlo(:,:,:)
 call sydy3(dyfrwf,indsym,natom,nsym,dyfrlo,symrec)

!Symmetrize nonlocal part of the dynamical matrix dyfrnl:
!atindx1 is used to reestablish the correct order of atoms
!Uses dyfrwf as work space, the nonlocal part is back to dyfrnl
 do iatom=1,natom
  dyfrwf(1:3,1:3,atindx1(iatom))=dyfrnl(1:3,1:3,iatom)
 end do
 call sydy3(dyfrwf,indsym,natom,nsym,dyfrnl,symrec)

!write(6, '(3es16.6)' )&
!&     (((dyfrwf(mu,nu,iatom),mu=1,3),nu=1,3),iatom=1,natom)

!Collect local, nl xc core, and non-local part
!of the frozen wf dynamical matrix.
 dyfrwf(:,:,:)=dyfrlo(:,:,:)+dyfrxc(:,:,:)+dyfrnl(:,:,:)

end subroutine dyfro3
!!***
