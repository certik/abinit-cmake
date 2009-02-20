!{\src2tex{textfont=tt}}
!!****f* ABINIT/cqratio
!! NAME
!! cqratio
!!
!! FUNCTION
!!  Calculate qratio(G,Gp,q)= (q+G)\cdot(q+Gp) / |q+G|^2
!!  needed for Hybertsen-Louie and Plasmonpole model
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (RS, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gmet(3,3)=metric in reciprocal space
!!  gprimd(3,3)=reciprocal lattice vectors
!!  gvec(3,npwc)=reduced coordinates of the plane waves 
!!  npwc=number of planewaves considered (used for the correlation part)
!!  nq=number of q points
!!  q(3,nq)=coordinates of q points
!!
!! OUTPUT
!!  qratio(npwc,npwc,nq)=(q+G).(q+Gp) needed for Hybertsen-Louie and 
!!  von der Linden-Horsh plasmonpole models
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cqratio(npwc,gvec,nq,q,gmet,gprimd,qratio)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwc,nq
!arrays
 integer,intent(in) :: gvec(3,npwc)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),q(3,nq)
 real(dp),intent(out) :: qratio(npwc,npwc,nq)

!Local variables ------------------------------
!scalars
 integer :: ig,igp,ii,iq
 real(dp) :: qpg_dot_qpgp
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3),norm(npwc)

!************************************************************************

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)

 norm(:)=zero ; qratio(:,:,:)=zero

 !this loops have to be rewritten!!!!
 do iq=1,nq
  do ig=1,npwc
   gpq(:)=gvec(:,ig)+q(:,iq)
   norm(ig)=two_pi*SQRT(DOT_PRODUCT(gpq,MATMUL(gmet,gpq)))
!  norm(ig)=normv(gpq,gmet,'g')
  end do
  do ig=1,npwc
   gpq(:)=gvec(:,ig)+q(:,iq)
   do igp=1,npwc
    gppq(:)=gvec(:,igp)+q(:,iq)
    qpg_dot_qpgp=zero
!   qpg_dot_qpgp=vdotw(gpq,gppq,gmet,'g')
    do ii=1,3
     qpg_dot_qpgp=qpg_dot_qpgp+&
&     ( gpq(1)*b1(ii) +  gpq(2)*b2(ii) + gpq(3)*b3(ii))*&
&     (gppq(1)*b1(ii) + gppq(2)*b2(ii) +gppq(3)*b3(ii))
    end do
!   
!   Now calculate qratio = (q+G).(q+Gp)/|q+G|^2 
!   when |q+G|^2 and (q+G).(q+Gp) are both zero
!   set (q+G).(q+Gp)/|q+G|^2 = 1
!   when |q+G|^2 is zero and |q+Gp| is not zero
!   set (q+G).(q+Gp)/|q+G|^2 = 0
!   
    if (norm(ig)<0.001) then
     if (norm(igp)<0.001) then
!     Case q=0 and G=Gp=0
      qratio(ig,igp,iq)=one
     else
!     Case q=0 and G=0 and Gp !=0
      qratio(ig,igp,iq)=zero
     end if
    else if (norm(igp)<0.001) then
!    Case q=0 and G= !0 and Gp=0
     qratio(ig,igp,iq)=zero
    else
     qratio(ig,igp,iq)=qpg_dot_qpgp/norm(ig)**2
    end if
   end do
  end do
 end do 

end subroutine cqratio
!!***
