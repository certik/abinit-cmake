!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyxc13
!! NAME
!! dyxc13
!!
!!
!! FUNCTION
!! Compute 2nd-order non-linear xc core-correction (part1)
!! to the dynamical matrix.
!! In case of derivative with respect to k or
!! electric field perturbation, the 1st-order local potential vanishes.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  kxc(nfft,nkxc)=first-order derivative of the xc potential
!!   if(nkxc=1): kxc(:,1)=dvxc/d$\rho$
!!   if(nkxc=3): kxc(:,1)=dvxc($\uparrow$)/d$\rho(\uparrow)$,
!!               kxc(:,2)=dvxc($\uparrow$)/d$\rho(\downarrow)$,
!!               kxc(:,3)=dvxc($\downarrow$)/d$\rho(\downarrow)$
!!   if(nkxc=23): GGA case, see rhohxc_coll.f
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=informations about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in cell.
!!  nattyp(ntypat)=number of atoms of each type in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(3)=fft grid dimensions.
!!  nkxc=second dimension of the kxc array
!!   (=1 for non-spin-polarized case, =3 for spin-polarized case)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xred(3,natom)=fractional coordinates for atoms in unit cell
!!
!! OUTPUT
!!  dyfrx1(2,3,natom,3,natom)=2nd-order non-linear xc
!!    core-correction (part1) part of the dynamical matrix
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      atm2fft3,dotprod_vn,mkcor3,mkvxc3
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dyxc13(atindx,dyfrx1,gmet,gprimd,gsqcut,kxc,mgfft,mpi_enreg,mqgrid,&
&          natom,nattyp,nfft,ngfft,nkxc,nspden,ntypat,n1xccc,paral_kgb,pawtab,&
&          ph1d,qgrid,qphon,rprimd,timrev,typat,ucvol,usepaw,xcccrc,xccc1d,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12spacepar
 use interfaces_13xc
 use interfaces_16response, except_this_one => dyxc13
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,n1xccc,natom,nfft,nkxc,nspden,ntypat
 integer,intent(in) :: paral_kgb,timrev,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),nattyp(ntypat),ngfft(18),typat(natom)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid),qphon(3)
 real(dp),intent(in) :: rprimd(3,3),xccc1d(n1xccc,6,ntypat),xcccrc(ntypat)
 real(dp),intent(in) :: xred(3,natom)
 real(dp),intent(out) :: dyfrx1(2,3,natom,3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom1,iatom2,idir1,idir2,ifft,ir,n1,n2,n3,n3xccc,nfftot
 integer :: option,optn,optn2,optv
 real(dp) :: kxcmean,valuei,valuer
!arrays
 real(dp),allocatable :: dummy(:),rhor1(:,:),vxc10(:,:),xcccwk1(:),xcccwk2(:)

! *********************************************************************

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 nfftot=n1*n2*n3

!Zero out the output array :
 dyfrx1(:,:,:,:,:)=zero

 cplex=2-timrev ; n3xccc=nfft
 allocate(vxc10(cplex*nfft,nspden))

 optv=0;optn=1;optn2=1

!Loop on the perturbation j1
 do iatom1=1,natom
  do idir1=1,3

!  Compute the derivative of the core charge with respect to j1
   allocate(xcccwk1(cplex*n3xccc))

!  PAW: 1st-order core charge in reciprocal space
   if (usepaw==1) then
    call atm2fft3(atindx,xcccwk1,dummy,cplex,dummy,gmet,gprimd,gsqcut,idir1,iatom1,&
&    mgfft,mpi_enreg,mqgrid,natom,nattyp,1,nfft,ngfft,&
&    ntypat,optn,optn2,optv,paral_kgb,pawtab,ph1d,qgrid,&
&    qphon,typat,ucvol,usepaw,dummy,xred)

!   Norm-conserving psp: 1st-order core charge in real space
   else
    call mkcor3(cplex,idir1,iatom1,natom,ntypat,n1,n1xccc,&
&    n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xcccwk1,xred)
   end if

!  Compute the corresponding potential
   option=0
   allocate(rhor1(cplex*nfft,nspden))   ! rhor1 is dummy
   call mkvxc3(cplex,gmet,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,nspden,&
&   n3xccc,option,paral_kgb,qphon,rhor1,rprimd,vxc10,xcccwk1)
   deallocate(rhor1,xcccwk1)

!  vxc10 will couple with xcccwk2, that behaves like
!  a total density (ispden=1). Only the spin-up + spin-down
!  average of vxc10 is needed.
   if (nspden/=1)then
    do ifft=1,cplex*nfft
     vxc10(ifft,1)=(vxc10(ifft,1)+vxc10(ifft,2))*half
    end do
   end if


!  Loop on the perturbation j2
!  (This way of proceeding does not take advantage of the
!  hermiticity of dyfrx1, which could save 50% of the CPU time)
   do iatom2=1,natom
    do idir2=1,3

!    Compute the derivative of the core charge with respect to j2
     allocate(xcccwk2(cplex*n3xccc))

!    PAW: 1st-order core charge in reciprocal space
     if (usepaw==1) then
      call atm2fft3(atindx,xcccwk2,dummy,cplex,dummy,gmet,gprimd,gsqcut,idir2,iatom2,&
&      mgfft,mpi_enreg,mqgrid,natom,nattyp,1,nfft,ngfft,&
&      ntypat,optn,optn2,optv,paral_kgb,pawtab,ph1d,qgrid,&
&      qphon,typat,ucvol,usepaw,dummy,xred)

!     Norm-conserving psp: 1st-order core charge in real space
     else
      call mkcor3(cplex,idir2,iatom2,natom,ntypat,n1,n1xccc,&
&      n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xcccwk2,xred)
     end if

!    Get the matrix element j1,j2

     call dotprod_vn(cplex,xcccwk2,valuer,valuei,mpi_enreg,nfft,nfftot,1,2,vxc10,ucvol)

     deallocate(xcccwk2)

     dyfrx1(1,idir1,iatom1,idir2,iatom2)=valuer
     dyfrx1(2,idir1,iatom1,idir2,iatom2)=valuei

    end do
   end do

  end do
 end do

 deallocate(vxc10)

end subroutine dyxc13
!!***
