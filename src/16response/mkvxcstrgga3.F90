!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkvxcstrgga3
!! NAME
!! mkvxcstrgga3
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! for the strain perturbation in case of GGA functionals
!! Use the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (DRH, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  dgprimdds(3,3)=strain derivaive of gprimd.
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor1tmp(cplex*nfft,2)=array for first-order electron spin-density
!!   in electrons/bohr**3 (second index corresponds to spin-up and spin-down)
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Closely related to mkvxcgga3.
!!
!! PARENTS
!!      mkvxcstr3
!!
!! CHILDREN
!!      xcden,xcpot
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkvxcstrgga3(cplex,dgprimdds,gmet,gprimd,gsqcut,kxc,mpi_enreg,nfft,ngfft,&
& nkxc,nspden,paral_kgb,qphon,rhor1tmp,vxc1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13xc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nkxc,nspden,paral_kgb
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: dgprimdds(3,3),gmet(3,3),gprimd(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: qphon(3),rhor1tmp(cplex*nfft,2)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,ishift,ispden,ngrad,nspdentmp,nspgrad
 real(dp) :: coeff_grho_corr,coeff_grho_dn,coeff_grho_up,coeffim_grho_corr
 real(dp) :: coeffim_grho_dn,coeffim_grho_up,exc1,gradrho_gradrho1
 real(dp) :: gradrho_gradrho1_dn,gradrho_gradrho1_up,gradrho_gradrho1im
 real(dp) :: gradrho_gradrho1im_dn,gradrho_gradrho1im_up,rho1_dn,rho1_up
 real(dp) :: rho1im_dn,rho1im_up,rho1re_dn,rho1re_up,rtmp
 character(len=500) :: message
!arrays
 real(dp) :: ec(7),ex(8),r0(3),r0_dn(3),r0_up(3),r1(3),r1_dn(3),r1_up(3)
 real(dp) :: r1im(3),r1im_dn(3),r1im_up(3),rho1rho1(3)
 real(dp),allocatable :: dnexcdn(:,:),rho1now(:,:,:),rhodgnow(:,:,:)
 real(dp),allocatable :: rhordgtmp(:,:),rhortmp(:,:,:),vxc1tmp(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' mkggavxc3 : enter '
!write(6,*)' mkggavxc3 : cplex,nfft,nkxc,nspden',cplex,nfft,nkxc,nspden
!if(cplex==2)then
!write(6,*)' mkvxcgga3 : not yet cplex=2, sorry'
!stop
!end if
!stop
!ENDDEBUG

!Treat all cases in a generic way (to be optimized !!)
 nspdentmp=2

!call filterpot(cplex,gmet,gsqcut,nfft,ngfft,2,qphon,rhor1tmp)

!Compute the gradients of the first-order density
 ishift=0 ; ngrad=2 ; nspdentmp=2
 allocate(rho1now(cplex*nfft,nspdentmp,ngrad*ngrad))
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspdentmp,paral_kgb,qphon,rhor1tmp,rho1now)
!Now, rho1now(:,:,1) contains the first-order density, and
!rho1now(:,:,2:4) contains the gradients of the first-order density

!Transfer the ground-state density and its gradient
!to spin-polarized storage
 allocate(rhortmp(nfft,nspdentmp,4))
 if(nspden==1)then
  do ii=1,4
   do ir=1,nfft
    rhortmp(ir,1,ii)=kxc(ir,14+2*ii)*half
    rhortmp(ir,2,ii)=kxc(ir,14+2*ii)*half
   end do
  end do
 else
  do ii=1,4
   do ir=1,nfft
    rhortmp(ir,1,ii)=kxc(ir,15+2*ii)
    rhortmp(ir,2,ii)=kxc(ir,14+2*ii)-kxc(ir,15+2*ii)
   end do
  end do
 end if

!Calculate the 1st-order contribution to grad n from the strain derivative
!acting on the gradient operator acting on the GS charge density.

 allocate(rhordgtmp(nfft,2),rhodgnow(cplex*nfft,nspdentmp,ngrad*ngrad))
 do ir=1,nfft
  rhordgtmp(ir,1)=rhortmp(ir,1,1)
  rhordgtmp(ir,2)=rhortmp(ir,2,1)
 end do
 call xcden(cplex,dgprimdds,ishift,mpi_enreg,nfft,ngfft,ngrad,nspdentmp,paral_kgb,qphon,&
& rhordgtmp,rhodgnow)
!Add to the gradients of the first-order density
 do ii=2,4
  do ir=1,nfft
   rho1now(ir,1,ii)=rho1now(ir,1,ii)+rhodgnow(ir,1,ii)
   rho1now(ir,2,ii)=rho1now(ir,2,ii)+rhodgnow(ir,2,ii)
  end do
 end do
 deallocate(rhordgtmp)


!Now, rhortmp(:,:,1) contains the GS density, and
!rhortmp(:,:,2:4) contains the gradients of the GS density
!Now, rho1now(:,:,1) contains the first-order density, and
!rho1now(:,:,2:4) contains the gradients of the first-order density

!Apply the XC kernel
 nspgrad=5
 allocate(dnexcdn(cplex*nfft,nspgrad))

 do ir=1,nfft
  r0_up(:)=rhortmp(ir,1,2:4)   ! grad of spin-up GS rho
  r0_dn(:)=rhortmp(ir,2,2:4)   ! grad of spin-down GS rho
  r0(:)=r0_up(:)+r0_dn(:)      ! grad of GS rho
  r1_up(:)=rho1now(ir,1,2:4)   ! grad of spin-up rho1
  r1_dn(:)=rho1now(ir,2,2:4)   ! grad of spin-down rho1
  r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
  gradrho_gradrho1_up=r1_up(1)*r0_up(1)+r1_up(2)*r0_up(2)+r1_up(3)*r0_up(3)
  gradrho_gradrho1_dn=r1_dn(1)*r0_dn(1)+r1_dn(2)*r0_dn(2)+r1_dn(3)*r0_dn(3)
  gradrho_gradrho1   =r1(1)*r0(1)+r1(2)*r0(2)+r1(3)*r0(3)

  dnexcdn(ir,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(ir,1,1)+&
&  kxc(ir,10)*rho1now(ir,2,1)+&
&  kxc(ir,5)*gradrho_gradrho1_up+&
&  kxc(ir,13)*gradrho_gradrho1
  dnexcdn(ir,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(ir,2,1)+&
&  kxc(ir,10)*rho1now(ir,1,1)+&
&  kxc(ir,6)*gradrho_gradrho1_dn+&
&  kxc(ir,14)*gradrho_gradrho1
  coeff_grho_corr=(kxc(ir,13)*rho1now(ir,1,1)+kxc(ir,14)*rho1now(ir,2,1))+&
&  kxc(ir,15)*gradrho_gradrho1
  coeff_grho_up= kxc(ir,5)*rho1now(ir,1,1)+kxc(ir,7)*gradrho_gradrho1_up
  coeff_grho_dn= kxc(ir,6)*rho1now(ir,2,1)+kxc(ir,8)*gradrho_gradrho1_dn

! grad strain derivative contribution enters the following term with a
! factor of two compared to above terms, so add it again.
  r1_up(:)=r1_up(:)+rhodgnow(ir,1,2:4)
  r1_dn(:)=r1_dn(:)+rhodgnow(ir,2,2:4)

! Reuse the storage in rho1now
  rho1now(ir,1,2:4)=r1_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&  +r1_dn(:)*kxc(ir,12)               &
&  +r0_up(:)*coeff_grho_up            &
&  +(r0_up(:)+r0_dn(:))*coeff_grho_corr
  rho1now(ir,2,2:4)=r1_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&  +r1_up(:)*kxc(ir,12)               &
&  +r0_dn(:)*coeff_grho_dn            &
&  +(r0_up(:)+r0_dn(:))*coeff_grho_corr
 end do
 deallocate(rhodgnow)

!Now, dnexcdn(:,1)=d(n.exc)/d(n_up)
!dnexcdn(:,2)=d(n.exc)/d(n_down)
!rho1now(:,:,2:4)=part of vxc1 that comes from FFT

 allocate(vxc1tmp(cplex*nfft,nspdentmp))
 vxc1tmp(:,:)=zero
 call xcpot (cplex,dnexcdn,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspdentmp,&
& nspgrad,paral_kgb,qphon,rho1now,vxc1tmp)

!Transfer the data from spin-polarized storage
 do ispden=1,nspden
  do ir=1,cplex*nfft
   vxc1(ir,ispden)=vxc1tmp(ir,ispden)
  end do
 end do

!call filterpot(cplex,gmet,gsqcut,nfft,ngfft,nspden,qphon,vxc1)

 deallocate(dnexcdn,rhortmp,rho1now,vxc1tmp)

end subroutine mkvxcstrgga3
!!***
