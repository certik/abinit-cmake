!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_ffm
!! NAME
!! calc_ffm
!!
!! FUNCTION
!! Calculate nth frequency moment of imaginary part of DM or its inverse
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (Rhaltaf,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nomega: total number of frequencies
!! npw: DM size
!! nq: total number of q vectors
!!
!!
!! OUTPUT
!! 
!! see side effects
!! 
!! side effects:
!! text file contains the nth frequency moment for each G, G '', q
!!
!!
!!
!! NOTES
!! RShaltaf
!! nfmidm : positive calculate the frequency moment for the DM
!! nfmidm : negative calcualte the frequency moment for the inverted
!! nfmidm : 0 calcualte the frequency moment for the full polarizibility
!! note that in a well converged results fmidm=1 and fmidm=-1 differ only in
!! minus sign
!! see M. Taut, J. Phys. C: Solid State Phys. 18 (1985) 2677-2690.
!!
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      leave_new,rhohxc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine calc_ffm(epsm1,nq,npw,nomega,omega,gprimd,qq,ppmodel,gvec,nfmidm)

 use defs_basis
 use defs_datatypes
!End of the abilint section


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_14occeig
 use interfaces_15gw, except_this_one => calc_ffm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfmidm,nomega,npw,nq,ppmodel
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: gprimd(3,3),qq(3,nq)
 complex(gwpc),intent(in) :: epsm1(npw,npw,nomega,nq),omega(nomega)

!Local variables ------------------------------
!scalars
 integer :: ig,igp,iomega,iq
 real(dp) :: omega_delta
 complex(gwpc) :: freqm
!arrays
 real(dp) :: b1(3),b2(3),b3(3)
 real(dp),allocatable :: fun(:),qplusg(:),res(:)
 complex(gwpc),allocatable :: eps(:,:,:)

!************************************************************************
 omega_delta=omega(2)-omega(1)

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)

 allocate(res(nomega),qplusg(npw),fun(nomega),eps(npw,npw,nomega))

 do iq=1,nq
  eps(:,:,:)=epsm1(:,:,:,iq)
! if fmidm is positive calcualte the inverted DM for each frequency
  do iomega=1,nomega
   if(nfmidm>=0)then
    call matcginv(eps(:,:,iomega),npw,npw)
   end if
  end do
! if nfmidm=0, calcualte the full polarizibility 
  if(nfmidm==0)then
   call cvc(nq,iq,qq,npw,gvec,gprimd,qplusg)
  end if
  do ig=1,npw
   if(nfmidm==0)then
    eps(ig,ig,:)=eps(ig,ig,:)-1
   end if
   do igp=1,npw
    if(nfmidm==0)then
     eps(ig,igp,:)=qplusg(ig)*qplusg(igp)*eps(ig,igp,:)/(4*pi)
    end if  
    fun(:)=real(omega(:)**abs(nfmidm))*aimag(eps(ig,igp,:))
    call simpson_int(nomega,omega_delta,fun,res)
    freqm=res(nomega)
    write(18,'(1x,i3,3x,i3,3x,i3,3x,f16.11,2x,f16.11)')ig,igp,iq,real(freqm),aimag(freqm)
   end do !igp
  end do  !ig
 end do ! iq
 rewind(18)

 deallocate(res,qplusg,fun)
end subroutine calc_ffm

!!***
