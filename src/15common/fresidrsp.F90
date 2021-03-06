!{\src2tex{textfont=tt}}
!!****f* ABINIT/fresidrsp
!!
!! NAME
!! fresidrsp
!!
!! FUNCTION
!! Compute the forces due to the residual of the potential (or density)
!! in RECIPROCAL SPACE, using
!!  - the atomic density read in psp file (PAW)
!!  - a gaussian atomic density (norm-conserving psps)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | densty(ntypat,4)=parameters for initialisation of the density of each atom type
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | ixc= choice of exchange-correlation scheme
!!  | natom=number of atoms in cell.
!!  | nspden=number of spin-density components
!!  | typat(natom)=integer type for each atom in cell
!! gmet(3,3)=reciprocal space metric
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! gsqcut=cutoff value on G**2 for sphere inside fft box
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=informations about MPI parallelization
!! mqgrid=number of grid pts in q array for atomic density spline n^AT(q)
!! nattyp(ntypat)=number of atoms of each type in cell
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! ntypat=number of types of atoms in cell.
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase information for given atom coordinates.
!! qgrid(mqgrid)=q grid for spline atomic valence density n^AT(q) from 0 to qmax
!! ucvol=unit cell volume (bohr**3).
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! vresid(nfft,nspden)=potential residual - (non-collinear magn. : only V11 and V22 are used)
!! zion(ntypat)=charge on each type of atom (real number)
!! znucl(ntypat)=atomic number, for each type of atom
!!
!! OUTPUT
!! gresid(3,natom)=forces due to the residual of the potential
!!
!! PARENTS
!!      forces
!!
!! CHILDREN
!!      atm2fft,fourdp,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,mpi_enreg,mqgrid,nattyp,nfft,&
&          ngfft,ntypat,pawtab,ph1d,qgrid,ucvol,usepaw,vresid,zion,znucl)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_15common, except_this_one => fresidrsp
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,nfft,ntypat,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: qgrid(mqgrid),vresid(nfft,dtset%nspden),zion(ntypat)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(out) :: gresid(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: itypat,optatm,optdyfr,optgr,optn,optn2,optstr,optv
 real(dp) :: ee
 logical :: usegauss
 character(len=500) :: message
!arrays
 integer :: dummy3(3)
 real(dp) :: dummy2(2),dummy6(6)
 real(dp),allocatable :: dummy(:),gauss(:,:),vresg(:,:),work(:)

! *************************************************************************

!Inits
 optatm=0;optdyfr=0;optgr=1;optstr=0;optv=0;optn=1
 allocate(vresg(2,nfft))

!Transfer potential residual to reciprocal space
!Use only Vres=Vres11+Vres22=Vres_up+Vres_dn
 allocate(work(nfft));work(:)=vresid(:,1);if (dtset%nspden>=2) work(:)=work(:)+vresid(:,2)
 call fourdp(1,vresg,work,-1,mpi_enreg,nfft,ngfft,dtset%paral_kgb,0)
 deallocate(work)

!Determine wether a gaussan atomic density has to be used or not
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

!Compute forces due to residual
 call atm2fft(atindx1,dummy,dummy,dummy,dummy,ee,gauss,gmet,gprimd,gresid,dummy,gsqcut,&
& mgfft,mpi_enreg,mqgrid,dtset%natom,nattyp,nfft,ngfft,ntypat,&
& optatm,optdyfr,optgr,optn,optn2,optstr,optv,dtset%paral_kgb,&
& pawtab,ph1d,qgrid,dummy3,dummy,dummy6,dummy6,&
& ucvol,usepaw,vresg,dummy2,dummy)

!In case of nspden>=2, has to apply 1/2 factor
 if (dtset%nspden>=2) gresid=gresid*half

 deallocate(gauss,vresg)

end subroutine fresidrsp

!!***
