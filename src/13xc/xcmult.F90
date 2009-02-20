!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcmult
!! NAME
!! xcmult
!!
!! FUNCTION
!! In the case of GGA, multiply the different gradient of spin-density
!! by the derivative of the XC functional with respect
!! to the norm of the gradient, then divide it by the
!! norm of the gradient
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dnexcdn(nfft,nspgrad)=derivative of Exc with respect to the (spin-)density,
!!    or to the norm of the gradient of the (spin-)density,
!!    further divided by the norm of the gradient of the (spin-)density
!!   The different components of dnexcdn will be
!!   for nspden=1,         dnexcdn(:,1)=d(n.exc)/d(n)
!!         and if ngrad=2, dnexcdn(:,2)=1/2*1/|grad n_up|*d(n.exc)/d(|grad n_up|)
!!                                      +   1/|grad n|*d(n.exc)/d(|grad n|)
!!         (do not forget : |grad n| /= |grad n_up| + |grad n_down|
!!   for nspden=2,         dnexcdn(:,1)=d(n.exc)/d(n_up)
!!                         dnexcdn(:,2)=d(n.exc)/d(n_down)
!!         and if ngrad=2, dnexcdn(:,3)=1/|grad n_up|*d(n.exc)/d(|grad n_up|)
!!                         dnexcdn(:,4)=1/|grad n_down|*d(n.exc)/d(|grad n_down|)
!!                         dnexcdn(:,5)=1/|grad n|*d(n.exc)/d(|grad n|)
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngrad = must be 2
!!  nspden=number of spin-density components
!!  nspgrad=number of spin-density and spin-density-gradient components
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  rhonow(nfft,nspden,ngrad*ngrad)=
!!   at input :
!!    electron (spin)-density in real space and its gradient,
!!    either on the unshifted grid (if ishift==0,
!!      then equal to rhor), or on the shifted grid
!!     rhonow(:,:,1)=electron density in electrons/bohr**3
!!     rhonow(:,:,2:4)=gradient of electron density in el./bohr**4
!!   at output :
!!    rhonow(:,:,2:4) has been multiplied by the proper factor,
!!    described above.
!!
!! PARENTS
!!      rhohxc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xcmult (dnexcdn,nfft,ngrad,nspden,nspgrad,rhonow)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ngrad,nspden,nspgrad
!arrays
 real(dp),intent(in) :: dnexcdn(nfft,nspgrad)
 real(dp),intent(inout) :: rhonow(nfft,nspden,ngrad*ngrad)

!Local variables-------------------------------
!scalars
 integer :: idir,ifft
 real(dp) :: rho_tot,rho_up

! *************************************************************************

 do idir=1,3

  if(nspden==1)then
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(dnexcdn,idir,nfft,rhonow)
   do ifft=1,nfft
    rhonow(ifft,1,1+idir)=rhonow(ifft,1,1+idir)*dnexcdn(ifft,2)
   end do
!  $OMP END PARALLEL DO

  else

!  In the spin-polarized case, there are more factors to take into account
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(dnexcdn,idir,nfft,rhonow)
   do ifft=1,nfft
    rho_tot=rhonow(ifft,1,1+idir)
    rho_up =rhonow(ifft,2,1+idir)
    rhonow(ifft,1,1+idir)=rho_up *dnexcdn(ifft,3)+&
&    rho_tot*dnexcdn(ifft,5)
    rhonow(ifft,2,1+idir)=(rho_tot-rho_up)*dnexcdn(ifft,4)+&
&    rho_tot*dnexcdn(ifft,5)
   end do
!  $OMP END PARALLEL DO

  end if ! nspden==1

! End loop on directions
 end do

end subroutine xcmult
!!***
