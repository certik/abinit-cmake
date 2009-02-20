!{\src2tex{textfont=tt}}
!!****f* ABINIT/cutoff_sphere
!! NAME
!! cutoff_sphere
!!
!! FUNCTION
!!  Calculate the Fourier transform of the coulombian interaction with a spherical cutoff   
!!  vc_cut(G)= \frac{4\pi}{|q+G|^2} [ 1-cos(|q+G|*R_cut) ]   (1)
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gmet(3,3)=metric in reciprocal space
!!  gvec(3,ng)=G vectors in reduced coordinates 
!!  rcut=cutoff radius
!!  ng=number of G vectors
!!  nq=number of q-points 
!!  qpt(3,nq)=q-points where the cutoff coulombian is required
!!
!! OUTPUT
!!  vc_cut(ng,nq)=Fourier components of the effective coulombian interaction 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  For |q|<small and G=0 we use 2pi.R_cut^2, i.e we take the limit q-->0 of Eq. (1) 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cutoff_sphere(nq,qpt,ng,gvec,gmet,rcut,vc_cut)

 use defs_basis
 use m_GWdefs, only : GW_TOLQ0


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,nq
 real(dp),intent(in) :: rcut
!arrays
 integer,intent(in) :: gvec(3,ng)
 real(dp),intent(in) :: gmet(3,3),qpt(3,nq)
 real(dp),intent(out) :: vc_cut(ng,nq)

!Local variables-------------------------------
!scalars
 integer :: ig,igs,iq
 real(dp) :: qpg
 character(len=500) :: msg

!************************************************************************

 if (ANY(gvec(:,1)/=0)) then 
  write(msg,'(4a)')ch10,&
&  ' cutoff_sphere: ERROR- ',ch10,&
&  ' first G vector should be zero '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 do iq=1,nq
  igs=1
  if (normv(qpt(:,iq),gmet,'G')<GW_TOLQ0) then  
   ! === For small q and G=0, use the limit q-->0  ===
   vc_cut(1,iq)=two_pi*rcut**2
   igs=2
  end if
  do ig=igs,ng
   qpg=normv(qpt(:,iq)+gvec(:,ig),gmet,'g')
   vc_cut(ig,iq)=four_pi*(one-COS(rcut*qpg))/qpg**2
  end do
 end do
 
end subroutine cutoff_sphere
!!***
