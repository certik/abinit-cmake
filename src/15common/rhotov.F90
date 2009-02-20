!{\src2tex{textfont=tt}}
!!****f* ABINIT/rhotov
!! NAME
!! rhotov
!!
!! FUNCTION
!! This routine is called to compute, from a given total density
!! the trial (local) potential and the residual potential.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | fixmom=input variable that governs fixed moment calculation
!!   | natom=number of atoms in cell.
!!   | nspden=number of spin-density components
!!   | ntypat=number of types of atoms in unit cell.
!!   | occopt=option for occupancies
!!   | typat(natom)=type (integer) for each atom
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nhatgr(nfft,nspden,3*nhatgrdim)= -PAW only- cartesian gradients of compensation density
!!  nhatgrdim= -PAW only- 0 if nhatgr array is not used ; 1 otherwise
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  optene=option for the computation of additional energies
!!  optres=0: the trial potential residual is computed ; the input potential value is kept
!!         1: the new value of the trial potential is computed in place of the input value
!!  optxc=option to be used for the call to rhohxc
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  ucvol = unit cell volume (Bohr**3)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vpsp(nfft)=array for holding local psp
!!  xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!!  ==== if optres==0
!!    vtrial(nfft,nspden)= old value of trial potential
!!
!! OUTPUT
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_hartree=Hartree part of total energy (hartree units)
!!   | e_xc=exchange-correlation energy (hartree)
!!  ==== if dtset%usewvl==1
!!   | e_vxc=the energy of the exchange-correlation potential
!!  ==== if optene==0.or.2
!!   | e_localpsp=local psp energy (hartree)
!!  ==== if optene==1.or.2
!!   | e_xcdc=exchange-correlation double-counting energy (hartree)
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if optxc==2.
!!  strsxc(6)=xc contribution to stress tensor (hartree/bohr^3)
!!  vxc(nfft,nspden)=Vxc(r) (already computed above; gets recomputed below too)
!!  vxcavg=mean of the vxc potential
!!  ==== if optres==0
!!    vresidnew(nfft,nspden)=potential residual
!!    vnew_mean(nspden)=mean of the potential formed from vpsp, vhartr and vxc, might be spin-dependent
!!    vres_mean(nspden)=mean of the potential residual, might be spin-dependent
!!    vres2=square of the norm of the residual
!!
!! SIDE EFFECTS
!! Input/Output:
!!  vhartr(nfft)=array for holding Hartree potential
!!  ==== if optres==1
!!    vtrial(nfft,nspden)= new value of trial potential
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
!!      prctfvw1,prctfvw2,prctfw3,scfcv
!!
!! CHILDREN
!!      dotprod_vn,mean_fftr,psolver_rhohxc,rhohxc,sqnorm_v,timab,wvl_newvtr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rhotov(dtset,energies,gsqcut,kxc,mpi_enreg,nfft,ngfft,&
&  nhat,nhatgr,nhatgrdim,nkxc,vresidnew,n3xccc,optene,optres,optxc,&
&  pawang,pawfgrtab,pawtab,rhog,rhor,rprimd,strsxc,ucvol,usepaw,usexcnhat,&
&  vhartr,vnew_mean,vpsp,vres_mean,vres2,vtrial,vxcavg,vxc,xccc3d)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_12spacepar
 use interfaces_13xc
 use interfaces_14poisson
 use interfaces_15common, except_this_one => rhotov
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n3xccc,nfft,nhatgrdim,nkxc,optene,optres,optxc,usepaw
 integer,intent(in) :: usexcnhat
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: vres2,vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(in) :: nhatgr(nfft,dtset%nspden,3*nhatgrdim),rhog(2,nfft)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),vhartr(nfft),vpsp(nfft)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc)
 real(dp),intent(out) :: kxc(nfft,nkxc),strsxc(6),vnew_mean(dtset%nspden)
 real(dp),intent(out) :: vres_mean(dtset%nspden),vresidnew(nfft,dtset%nspden)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat)

!Local variables-------------------------------
!scalars
 integer :: ifft,ispden,nfft_loc,nfftot,offset
 real(dp) :: doti,dum,ucvol_local
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2),vmean(dtset%nspden)
 real(dp),allocatable :: vnew(:,:)

! *********************************************************************

!DEBUG
!write(6,*)' rhotov : enter, write rprimd '
!write(6,*)rprimd
!ENDDEBUG

 call timab(57,1,tsec)

!Get size of FFT grid
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)

 if (dtset%usewvl == 0) then
  ucvol_local = ucvol
 else
! We need to tune the volume when wavelets are used because, not
! all FFT points are used.
  nfftot = product(dtset%wvl_internal%dpSize)
  ucvol_local = (dtset%wvl_hgrid / real(2, dp)) ** 3 * real(nfftot, dp)
 end if


!------Compute Hartree and xc potentials----------------------------------

!Compute xc potential (separate up and down if spin-polarized)
 if (dtset%icoulomb == 0) then
! Use the periodic solver to compute Hxc.
  call rhohxc(dtset,energies%e_xc,gsqcut,usepaw,kxc,mpi_enreg,nfft,ngfft,&
&  nhat,usepaw,nhatgr,nhatgrdim,nkxc,dtset%nspden,n3xccc,optxc,rhog,&
&  rhor,rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d)
  call dotprod_vn(1,rhor,energies%e_hartree,doti,mpi_enreg,nfft,nfftot,1,1,vhartr,ucvol)
  energies%e_hartree=half*energies%e_hartree
 else
! Use the free boundary solver.
  call PSolver_rhohxc(dtset, energies%e_hartree, energies%e_xc, energies%e_vxc, &
&  mpi_enreg, rhor, rprimd, vhartr, vxc, vxcavg)
 end if

!------Compute parts of total energy depending on potentials--------
 if (optene==0.or.optene==2) then

! Compute local psp energy energies%e_localpsp
  call dotprod_vn(1,rhor,energies%e_localpsp,doti,mpi_enreg,nfft,nfftot,1,1,vpsp,ucvol_local)

 end if

 if (optene==1.or.optene==2) then

! Compute double-counting XC energy energies%e_xcdc
  call dotprod_vn(1,rhor,energies%e_xcdc,doti,mpi_enreg,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local)

 end if


!------Produce residual vector and square norm of it-------------
!(only if requested ; if optres==0)

 if (optres==0) then
! Compute potential residual
  allocate(vnew(nfft,dtset%nspden))
  vmean(:)=zero ; vnew_mean(:)=zero

  if (dtset%usewvl == 0) then
   do ispden=1,min(dtset%nspden,2)
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(ispden,nfft,dtset%nspden,vhartr,vnew,vpsp,vresidnew,vtrial,vxc)
    do ifft=1,nfft
     vnew(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
     vresidnew(ifft,ispden)=vnew(ifft,ispden)-vtrial(ifft,ispden)
    end do
!   $OMP END PARALLEL DO
   end do
   if(dtset%nspden==4)then
    do ispden=3,4
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(ispden,nfft,vresidnew,vtrial,vxc)
     do ifft=1,nfft
      vnew(ifft,ispden)=vxc(ifft,ispden)
      vresidnew(ifft,ispden)=vxc(ifft,ispden)-vtrial(ifft,ispden)
     end do
!    $OMP END PARALLEL DO
    end do
   end if
   offset   = 0
   nfft_loc = nfft
  else
   call wvl_newvtr(dtset, mpi_enreg, nfft_loc, offset, vhartr, vpsp, vnew, vxc)
   vresidnew = vnew - vtrial
  end if
! Compute mean values of potential and residual
  call mean_fftr(vnew(1+offset, 1),vnew_mean,mpi_enreg,nfft_loc,nfftot,dtset%nspden)
  call mean_fftr(vresidnew(1+offset, 1),vmean,mpi_enreg,nfft_loc,nfftot,dtset%nspden)
  deallocate(vnew)

! Subtract the mean of the residual
! Must take into account fixed occupation number in case of spin-polarized
  do ispden=1,dtset%nspden
   if (dtset%nspden==2.and.dtset%occopt>=3.and. &
&   abs(dtset%fixmom+99.99_dp)<1.0d-10)then
    vres_mean(ispden)=(vmean(1)+vmean(2))*half
   else
    vres_mean(ispden)=vmean(ispden)
   end if
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(ispden,nfft,vresidnew,vres_mean)
   do ifft=1,nfft
    vresidnew(ifft,ispden)=vresidnew(ifft,ispden)-vres_mean(ispden)
   end do
!  $OMP END PARALLEL DO
  end do

! Compute square norm vres2 of potential residual vresid
  call sqnorm_v(1,mpi_enreg,nfft_loc,vres2,dtset%nspden,vresidnew(1+offset, 1))

 else    ! optres==0


! ------Produce new value of trial potential-------------
! (only if requested ; if optres==1)

  if (dtset%usewvl == 0) then
   do ispden=1,min(dtset%nspden,2)
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(ispden,nfft,dtset%nspden,vhartr,vnew,vpsp,vxc)
    do ifft=1,nfft
     vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
    end do
!   $OMP END PARALLEL DO
   end do
   if(dtset%nspden==4) vtrial(:,3:4)=vxc(:,3:4)
  else
!  output offset and nfft_loc are unused here.
   call wvl_newvtr(dtset, mpi_enreg, nfft_loc, offset, vhartr, vpsp, vtrial, vxc)
  end if

 end if

 call timab(57,2,tsec)

!DEBUG
!write(6,*)' rhotov : exit '
!stop
!ENDDEBUG

end subroutine rhotov
!!***
