!{\src2tex{textfont=tt}}
!!****f* ABINIT/cutoff_surface
!! NAME
!! cutoff_surface
!!
!! FUNCTION
!!  Calculate the Fourier components of an effective Coulombian interaction 
!!  within a slab of thickness 2*rcut which is symmetric with respect to the xy plane 
!!  In this implementation rcut=L_z/2 where L_z is the periodicity along z 
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  gvec(3,ng)=G vectors in reduced coordinates 
!!  ng=number of G vectors
!!  qpt(3,nq)= q-points 
!!  nq=number of q-points 
!!  gmet(3,3)=metric in reciprocal space
!!
!! OUTPUT
!!  vc_cut(ng,nq)=Fourier components of the effective coulombian interaction 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  The Fourier expression for an interaction truncated along the z-direction 
!!  (i.e non-zero only if |z|<R) is :
!!  
!!  vc(q.G) = 4pi/|q+G|^2 * [ 1 + e^{-((q+G)_xy)*R} * ( (q_z+G_z)/(q+G)_xy * sin((q_z+G_z)R) - 
!!   - cos((q_z+G_Z)R)) ]  (1)
!! 
!!  Equation (1) diverges when q_xy+G_xy --> 0 for any non zero q_z+G_z
!!  However if we choose R=L/2, where L defines the periodicity along z, 
!!  and we limit ourselves to consider q-points such as q_z==0, then 
!!  sin((q_z+G_z)R)=sin(G_Z 2pi/L)=0 for every G.  
!!  Under these assumptions we obtain
!!
!!  v(q,G) = 4pi/|q+G|^2} [ 1-e^{-(q+G)_xy*L/2}\cos((q_z+G_z)R) ]
!! 
!!  which is always finite when G_z /=0 
!!  while it diverges as 4piR/(q+G)_xy as (q+G)_xy -->0 but only in the plane x-y
!!
!! PARENTS
!!  
!!
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cutoff_surface(nq,qpt,ng,gvec,gprimd,gmet,rcut,boxcenter,pdir,alpha,vc_cut,method)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,ng,nq
 real(dp),intent(in) :: rcut
!arrays
 integer,intent(in) :: gvec(3,ng),pdir(3)
 real(dp),intent(in) :: alpha(3),boxcenter(3),gmet(3,3),gprimd(3,3),qpt(3,nq)
 real(dp),intent(out) :: vc_cut(ng,nq)

!Local variables-------------------------------
!scalars
 integer :: ig,igs,iq
 real(dp),parameter :: SMALL=tol6
 real(dp) :: ap1sqrt,qpg2,qpg_para,qpg_perp
 character(len=500) :: msg
!arrays
 real(dp) :: b1(3),b2(3),b3(3),gcart(3),qc(3),qpg(3)
 real(dp),allocatable :: qcart(:,:)

! *************************************************************************

 ! === From reduced to cartesian coordinates ===
 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)
 allocate(qcart(3,nq))
 do iq=1,nq
  qcart(:,iq)=b1(:)*qpt(1,iq)+b2(:)*qpt(2,iq)+b3(:)*qpt(3,iq)
 end do
 !
 ! === Different approaches according to method ===
 vc_cut(:,:)=zero 
 SELECT CASE (method)
 CASE (1)
  ! === Beigi expression ===
  ! * q-points with non-zero component along the z-axix are not allowed if 
  !  the simplified Eq.l for the Coulombian interaction is used
  if (ANY(ABS(qcart(3,:))>SMALL)) then 
   write(msg,'(8a)')ch10,&
&   ' cutoff_surface : ERROR -',ch10,&
&   ' found q-points with non-zero component along non-periodic direction ',ch10,&
&   ' This is not allowed, see Notes in cutoff_surface.F90 ',ch10,&
&   ' Modify your q-point sampling '
   call wrtout(std_out,msg,'COLL') ; write(*,*)qcart(:,:) ; call leave_new('COLL')
  end if 
  ! 
  ! === Calculate truncated coulombian interaction for a infinite surface ===
  ! Here I suppose that all the input q-points are different from zero
  do iq=1,nq
   qc(:)=qcart(:,iq)
   do ig=1,ng
    gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
    qpg(:)=qc(:)+gcart(:)
    qpg2  =DOT_PRODUCT(qpg(:),qpg(:))
    qpg_para=SQRT(qpg(1)**2+qpg(2)**2) ; qpg_perp=qpg(3)
    ! if (abs(qpg_perp)<SMALL.and.qpg_para<SMALL) stop 'SMALL in cutoff_surface
    vc_cut(ig,iq)=four_pi/qpg2*(one-EXP(-qpg_para*rcut)*COS(qpg_perp*rcut))  
   end do 
  end do
 CASE (2)
  ! === Rozzi et al method ===
  !alpha=?? ; ap1sqrt=SQRT(one+alpha**2)
  STOP "working in progress"
  do iq=1,nq
   qc(:)=qcart(:,iq)
   do ig=1,ng
    gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
    qpg(:)=qc(:)+gcart(:)
    qpg2  =DOT_PRODUCT(qpg(:),qpg(:))
    qpg_para=SQRT(qpg(1)**2+qpg(2)**2) ; qpg_perp =qpg(3)
    if (qpg_para>SMALL) then 
     vc_cut(ig,iq)=four_pi/qpg2*(one+EXP(-qpg_para*rcut)*(qpg_perp/qpg_para*SIN(qpg_perp*rcut)-COS(qpg_perp*rcut))) 
    else 
     if (ABS(qpg_perp)>SMALL) then 
      vc_cut(ig,iq)=four_pi/qpg_perp**2*(one-COS(qpg_perp*rcut)-qpg_perp*rcut*SIN(qpg_perp*rcut)) !&
!&      + 8*rcut*SIN(qpg_perp*rcut)/qpg_perp*LOG((alpha+ap1sqrt)*(one+ap1sqrt)/alpha) ! contribution due to finite surface
     else 
      vc_cut(ig,iq)=-two_pi*rcut**2
     end if 
    end if 
   end do !ig 
  end do !iq
 CASE DEFAULT 
  write(msg,'(4a,i3)')ch10,&
&  ' cutoff_surface : BUG - ',ch10,&
&  ' wrong value of method: ',method 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT
 deallocate(qcart)

end subroutine cutoff_surface
!!***
