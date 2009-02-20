!{\src2tex{textfont=tt}}
!!****m* ABINIT/finite_cylinder
module finite_cylinder

 use defs_basis
 use m_numeric_tools, only : quadrature

 implicit none

 integer :: npts_,ntrial_,qopt_
 real(dp) :: ha_,hb_,r0_
 real(dp) :: qpg_perp_,qpg_para_,qpgx_,qpgy_
 real(dp) :: zz_,xx_,rho_
 real(dp) :: hcyl_,rcut_,accuracy_

contains 


function F2(xx)
!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: F2
!Local variables-------------------------------
!scalars
 real(dp) :: intr
!************************************************************************

 zz_=xx
 call quadrature(F1,zero,rcut_,qopt_,intr,ntrial_,accuracy_,npts_)
 F2=intr*COS(qpg_para_*xx)

end function F2


function F1(rho) 
!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 real(dp),intent(in) :: rho
 real(dp) :: F1

!Local variables-------------------------------
!scalars
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 !F1(\rho;z)= \rho*j_o(qpg_perp_*\rho)/sqrt(\rho**2+z**2)
 arg=rho*qpg_perp_
 call jbessel(bes,besp,bespp,ll,order,arg)
 if (zz_==zero) then 
  F1=bes
 else 
  F1=bes*rho/SQRT(rho**2+zz_**2)
 end if

end function F1


function F3(xx)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: xx
 real(dp) :: F3
!************************************************************************

 !$F3(z)=z*\sin(qpg_para_*z)/\sqrt(rcut^2+z^2)$
 F3=xx*SIN(qpg_para_*xx)/SQRT(rcut_**2+xx**2)

end function F3


function F4(rho)
!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 real(dp),intent(in) :: rho
 real(dp) :: F4

!Local variables-------------------------------
!scalars
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 !$F4(rho)=\rho*j_o(qpg_perp_.\rho) \ln((hcyl+\sqrt(rho^2+hcyl^2))/\rho)$
 if (ABS(rho)<tol12) then 
  F4=zero
 else
  arg=rho*qpg_perp_
  call jbessel(bes,besp,bespp,ll,order,arg)
  F4=bes*rho*LOG((hcyl_+SQRT(rho**2+hcyl_**2))/rho)
 end if

end function F4


function F5(rho)
!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 real(dp),intent(in) :: rho
 real(dp) :: F5

!Local variables-------------------------------
!scalars
 integer,parameter :: order=0,ll=0
 real(dp) :: arg,bes,besp,bespp
!************************************************************************

 ! $F5(\rho)=\rho*j_o(G_perp\rho)log(\rho)$
 if (rho==0) then 
  F5=zero
 else 
  arg=rho*qpg_perp_
  call jbessel(bes,besp,bespp,ll,order,arg)
  F5=bes*rho*LOG(rho)
 end if

end function F5


function K0cos(yy) 
!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib00numeric
!End of the abilint section

 real(dp),intent(in) :: yy
 real(dp) :: K0cos

!Local variables-------------------------------
!scalars
 real(dp) :: k0,rho,arg
!************************************************************************

 !$K0cos(y)=K0(\rho*|qpg_z|)*COS(x.qpg_x+y*qpg_y)$
 rho=SQRT(xx_**2+yy**2) ; arg=qpg_para_*rho
 call CALCK0(arg,k0,1)
 K0cos=k0*COS(qpgx_*xx_+qpgy_*yy)

end function K0cos


function K0cos_dy(xx)
!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: K0cos_dy
!Local variables-------------------------------
!scalars
 real(dp) :: bb,quad
!************************************************************************

 !$ K0cos_dy(x)=\int_{-b/2}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
 xx_=xx 
 call quadrature(K0cos,-hb_,+hb_,qopt_,quad,ntrial_,accuracy_,npts_)
 K0cos_dy=quad

end function K0cos_dy


function K0cos_dy_r0(xx)
!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 real(dp),intent(in) :: xx
 real(dp) :: K0cos_dy_r0
!Local variables-------------------------------
!scalars
 real(dp) :: bb,quad,yx
!************************************************************************

 !$K0cos_dy_r0(x)= \int_{-b/2}^{-y(x)} K0(|qpg_z|\rho) cos(x.qpg_x+y.qpg_y)dy 
 !                 +\int_{y(x)}^{b/2} K0(|qpg_z|\rho)cos(x.qpg_x+y.qpg_y)dy$
 ! where y(x)=SQRT(r0^2-x^2) and x<=r0
 !
 xx_=xx ; yx=SQRT(r0_**2-xx**2)
 call quadrature(K0cos,-hb_,-yx,qopt_,quad,ntrial_,accuracy_,npts_)
 K0cos_dy_r0=quad
 call quadrature(K0cos,+yx,+hb_,qopt_,quad,ntrial_,accuracy_,npts_)
 K0cos_dy_r0=quad+K0cos_dy_r0

end function K0cos_dy_r0


function K0cos_dth_r0(rho)
!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib00numeric
!End of the abilint section

 real(dp),intent(in) :: rho
 real(dp) :: K0cos_dth_r0
!Local variables-------------------------------
!scalars
 real(dp) :: bb,quad,arg,k0,tmp
!************************************************************************

 ! $ K0cos_dth_r0(\rho)= 
 ! \int_{0}^{2pi)} K0(|qpg_z|\rho)cos(\rho.cos(\theta).qpg_x+\rho.sin(\theta).qpg_y) d\theta
 ! where y(x)=SQRT(r0^2-x^2) and x<=r0
 !
 rho_=rho 
 call quadrature(Fcos_th,zero,two_pi,qopt_,quad,ntrial_,accuracy_,npts_)

 arg=qpg_para_*rho_ 
 tmp=zero 
 if (arg>tol6) then
  call CALCK0(arg,k0,1)
  tmp=k0*rho_
 end if
 K0cos_dth_r0=quad*tmp

end function K0cos_dth_r0


function Fcos_th(theta) 
!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: theta
 real(dp) :: Fcos_th

!Local variables-------------------------------
!scalars
 real(dp) :: k0,arg,tmp
!************************************************************************

 !$Fcos_th(\theta)=rho*K0(\rho*|qpg_z|)*COS(\rho.COS(\theta).qpg_x+\rho.SIN/(\theta)*qpg_y)$
 !arg=qpg_para_*rho_ 
 !call CALCK0(arg,k0,1)
 !tmp=k0*rho_
 Fcos_th=COS(rho_*COS(theta)*qpgx_+rho_*SIN(theta)*qpgy_)

end function Fcos_th


end module finite_cylinder
!!***


!!****f* ABINIT/cutoff_cylinder
!! NAME
!! cutoff_cylinder
!!
!! FUNCTION
!!  Calculate the Fourier components of an effective Coulombian interaction 
!!  that is zero outside a finite region in the x-y plane
!!  Two methods are implemented: 
!!   In the first case, method=1, the interaction in the x-y plane is truncated outside the Wigner-Seitz 
!!   cell centered on the wire in the x-y plane. The interaction has infinite extent along the z axis 
!!   and the Fourier transform in singular only at the Gamma point. Only orthorombic Bravais lattice are supported.
!!   In the second case the interaction is truncated outside a cylinder of radius rcut, the cylinder has finite extent along z.
!!   and no singularity occurs. 
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  boxcenter(3)= center of the wire in the x-y axis
!!  gvec(3,ng)=G vectors in reduced coordinates 
!!  ng=number of G vectors
!!  qpt(3,nq)= q-points 
!!  MPI_enreg= datatype containing information on parallelism
!!  nq=number of q-points 
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!   where: rprimd(i,j)=rprim(i,j)*acell(j)
!!  method=1 for Beigi approach (infinite cylinder with interaction truncated outside the W-S cell)
!!         2 for Rozzi method (finite cylinder)
!!
!! OUTPUT
!!  vccut(ng,nq)= Fourier components of the effective coulombian interaction 
!!
!! SIDE EFFECTS
!!
!! NOTES
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

subroutine cutoff_cylinder(nq,qpt,ng,gvec,rcut,hcyl,pdir,&
& boxcenter,rprimd,vccut,method,MPI_enreg)

 use defs_basis
 use m_gwdefs, only : GW_TOLQ0
 use defs_datatypes
 use m_numeric_tools, only : quadrature
 use finite_cylinder


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_lib00numeric
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ng,nq,method
 real(dp),intent(in) :: rcut,hcyl 
 type(MPI_type),intent(in) :: MPI_enreg
!arrays
 integer,intent(in) :: gvec(3,ng),pdir(3)
 real(dp),intent(in) :: boxcenter(3),qpt(3,nq),rprimd(3,3)
 real(dp),intent(out) :: vccut(ng,nq)

!Local variables-------------------------------
!scalars
 integer,parameter :: N0=1000
 integer :: i1,i2,i3,icount,ig,igs,iq,itrial,ix,spaceComm
 integer :: nx,ntasks,rank,master,nprocs,ierr,my_start,my_stop
 real(dp) :: j0,j1,k0,k1,nmc,old_quad,qpg2,qpg_xy,tmp
 real(dp) :: qpg_z,quad,rcut2,hcyl2,c1,c2,ucvol,SMALL
 logical :: converged
 character(len=500) :: msg
!arrays
 real(dp) :: qpg(3),b1(3),b2(3),b3(3),gmet(3,3),rmet(3,3),gprimd(3,3),qc(3),gcart(3)
 real(dp),allocatable :: qcart(:,:)
!************************************************************************
 
 ! ===================================================
 ! === Setup for the quadrature of matrix elements ===
 ! ===================================================
 qopt_=6         ! Quadrature method, see quadrature
 ntrial_=30      ! Max number of attempts
 accuracy_=0.001 ! Fractional accuracy required
 npts_=6         ! Initial number of point (only for Gauss-Legendre)
 SMALL=GW_TOLQ0  ! Below this value (q+G)_i is treated as zero 
 rcut_=rcut      ! Radial cutoff, used only if method==2
 hcyl_=hcyl      ! Lenght of cylinder along z, only if method==2

 write(msg,'(3a,2(a,i5,a),a,f8.7)')ch10,&
& ' cutoff_cylinder: info on quadrature method : ',ch10,&
& ' quadrature scheme   = ',qopt_,ch10,&
& ' ntrial              = ',ntrial_,ch10,&
& ' accuracy            = ',accuracy_
 call wrtout(std_out,msg,'COLL') 
 !
 ! === From reduced to Cartesian coordinates ===
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 b1(:)=two_pi*gprimd(:,1)
 b2(:)=two_pi*gprimd(:,2)
 b3(:)=two_pi*gprimd(:,3)
 allocate(qcart(3,nq))
 do iq=1,nq
  qcart(:,iq)=b1(:)*qpt(1,iq)+b2(:)*qpt(2,iq)+b3(:)*qpt(3,iq)
 end do

 vccut(:,:)=zero ; ntasks=nq*ng 
 call xcomm_init(MPI_enreg,spaceComm)  
 call xme_init  (MPI_enreg,rank     )          
 call split_work(ntasks,my_start,my_stop)
 !
 ! === Different approaches according to method ===
 SELECT CASE (method)

 CASE (1)
  ! === Infinite cylinder, interaction is zero outside the Wigner-Seitz cell ===
  write(msg,'(2(a,f8.4))')' cutoff_cylinder: Using Beigi''s Infinite cylinder '
  call wrtout(std_out,msg,'COLL') 
  if (ANY(qcart(1:2,:)>SMALL)) then 
   ! 1) The expression of Beigi can be used only if the BZ is sampled only along z
   write(msg,'(8a)')ch10,&
&   ' cutoff_cylinder : ERROR -',ch10,&
&   ' found q-points with non zero components in the X-Y plane.',ch10,&
&   ' This is not allowed, see Notes in cutoff_cylinder.F90.',ch10,&
&   ' Modify your q-point sampling '
   call wrtout(std_out,msg,'COLL') ; write(*,*)qcart(:,:) ; call leave_new('COLL')
  end if 
  ! 2) Check if Bravais lattice is orthorombic and parallel to the Cartesian versors
  ! In this case the intersection of the W-S cell with the x-y plane is the rectangle with 
  ! -ha_<=x<=ha_ and -hb_<=y<=hb_ 
  if ((ANY(rprimd(2:3,1)/=zero)).or.(ANY(rprimd(1:3:2,2)/=zero)).or.((ANY(rprimd(1:2,3)/=zero)))) then
   write(msg,'(2a)')ch10,&
& ' cutoff_cylinder : Bravais lattice should be orthorombic and parallel to the cartesian versors '
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if
  ha_=half*SQRT(DOT_PRODUCT(rprimd(:,1),rprimd(:,1)))
  hb_=half*SQRT(DOT_PRODUCT(rprimd(:,2),rprimd(:,2)))
  r0_=MIN(ha_,hb_)/N0
  ! 
  ! === For each pair (q,G) evaluate the integral defining the cutoff coulombian  ===
  ! * Here I assume that all q-vectors are non zero and that q_xy/=0
  do iq=1,nq
   ! === Skip singularity at Gamma, it will be treated "by hand" in csigme ===
   igs=1 
   !if (normv(qpt(:,iq),gmet,'G')<GW_TOLQ0) igs=2
   qc(:)=qcart(:,iq)
   !qc(1:2)=zero
   write(msg,'(2(a,i3))')' entering loop iq: ',iq,' with igs = ',igs
   call wrtout(std_out,msg,'COLL')
   do ig=igs,ng
    icount=ig+(iq-1)*ng ; if (icount<my_start.or.icount>my_stop) CYCLE
    gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
    qpg(:)=qc(:)+gcart(:)
    qpgx_=qpg(1) ; qpgy_=qpg(2) ; qpg_para_=ABS(qpg(3))
    !
    ! Calculate $2 \int_{WS} dx dy K_0{qpg_para_\rho) cos(x.qpg_x + y.qpg_y) where WS is the Wigner-Seitz cell 
    tmp=zero
    ! Difficult part, integrate on a small cirle of radius r0 using spherical coordinates
    !call quadrature(K0cos_dth_r0,zero,r0_,qopt_,quad,ntrial_,accuracy_,npts_)
    !write(*,'(i8,a,es14.6)')ig,' 1 ',quad
    !tmp=tmp+quad
    ! Add region with 0<=x<=r0 and y>=+-(SQRT(r0^2-x^2))since WS is rectangular
    !call quadrature(K0cos_dy_r0,zero,r0_,qopt_,quad,ntrial_,accuracy_,npts_)
    !write(*,'(i8,a,es14.6)')ig,' 2 ',quad
    !tmp=tmp+quad
    ! Get the in integral in the rectangle with x>=r0, should be the easiest but sometimes has problems to converge
    !call quadrature(K0cos_dy,r0_,ha_,qopt_,quad,ntrial_,accuracy_,npts_)
    !write(*,'(i8,a,es14.6)')ig,' 3 ',quad
    ! 
    ! === More stable method: midpoint integration with Romberg extrapolation ===
    call quadrature(K0cos_dy,zero,ha_,qopt_,quad,ntrial_,accuracy_,npts_)
    !write(*,'(i8,a,es14.6)')ig,' 3 ',quad
    tmp=tmp+quad
    vccut(ig,iq)=two*(tmp*two) 
    ! Factor two comes from the replacement WS -> (1,4) quadrant thanks to symmetries of the integrad
   end do 
  end do

 CASE (2) 
  ! === Finite cylinder of lenght hcyl, from Rozzi et al ===
  ! TODO add check on hcyl value that should be smaller that 1/deltaq 
  if (hcyl_<zero) then
   write(msg,'(4a)')ch10,&
&   ' cutoff_cylinder : BUG - ',ch10,&
&   ' negative value for cylinder lenght '
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if

  if (ABS(hcyl_)>tol12) then 
   write(msg,'(2(a,f8.4))')' cutoff_cylinder: using finite cylinder of length= ',hcyl_,' rcut= ',rcut_
   call wrtout(std_out,msg,'COLL') 
   hcyl2=hcyl_**2 ; rcut2=rcut_**2
   do iq=1,nq
    write(msg,'(a,i3)')' entering loop iq: ',iq
    call wrtout(std_out,msg,'COLL')
    qc(:)=qcart(:,iq)
    do ig=1,ng
     ! === No singularity occurs in finite cylinder, thus start from 1 ===
     icount=ig+(iq-1)*ng ; if (icount<my_start.or.icount>my_stop) CYCLE
     gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
     qpg(:)=qc(:)+gcart(:)
     qpg_para_=ABS(qpg(3)) ; qpg_perp_=SQRT(qpg(1)**2+qpg(2)**2)

     if (qpg_perp_/=zero.and.qpg_para_/=zero) then 
      ! $4\pi\int_0^{R_c} d\rho\rho j_o(qpg_perp_.\rho)\int_0^hcyl dz\cos(qpg_para_*z)/sqrt(\rho^2+z^2)$ 
      call quadrature(F2,zero,rcut_,qopt_,quad,ntrial_,accuracy_,npts_)
      vccut(ig,iq)=four_pi*quad

     else if (qpg_perp_==zero.and.qpg_para_/=zero) then 
      ! $\int_0^h sin(qpg_para_.z)/\sqrt(rcut^2+z^2)dz$
      call quadrature(F3,zero,hcyl_,qopt_,quad,ntrial_,accuracy_,npts_)
      c1=one/qpg_para_**2-COS(qpg_para_*hcyl_)/qpg_para_**2-hcyl_*SIN(qpg_para_*hcyl_)/qpg_para_
      c2=SIN(qpg_para_*hcyl_)*SQRT(hcyl2+rcut2)
      vccut(ig,iq)=four_pi*c1+four_pi*(c2-quad)/qpg_para_

     else if (qpg_perp_/=zero.and.qpg_para_==zero) then 
      ! $4pi\int_0^rcut d\rho \rho J_o(qpg_perp_.\rho) ln((h+\sqrt(h^2+\rho^2))/\rho)
      call quadrature(F4,zero,rcut_,qopt_,quad,ntrial_,accuracy_,npts_)
      vccut(ig,iq)=four_pi*quad

     else if (qpg_perp_==zero.and.qpg_para_==zero) then 
      ! Use lim q+G --> 0
      vccut(ig,iq)=two_pi*(-hcyl2+hcyl_*SQRT(hcyl2+rcut2)+rcut2*LOG((hcyl_+SQRT(hcyl_+SQRT(hcyl2+rcut2)))/rcut_))

     else 
      write(msg,'(a)')' cutoff_cylinder : BUG , you should not be here!'
      call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
     end if

    end do !ig
   end do !iq

  else 
   ! === Infinite cylinder ===
   write(msg,'(a)')' cutoff_cylinder : using Rozzi''s method with infinite cylinder '
   call wrtout(std_out,msg,'COLL') 
   do iq=1,nq
    write(msg,'(a,i3)')' entering loop iq: ',iq
    call wrtout(std_out,msg,'COLL')
    qc(:)=qcart(:,iq)
    do ig=1,ng
     icount=ig+(iq-1)*ng ; if (icount<my_start.or.icount>my_stop) CYCLE
     gcart(:)=b1(:)*gvec(1,ig)+b2(:)*gvec(2,ig)+b3(:)*gvec(3,ig)
     qpg(:)=qc(:)+gcart(:)
     qpg2  =DOT_PRODUCT(qpg,qpg)
     qpg_z =ABS(qpg(3)) ; qpg_xy=SQRT(qpg(1)**2+qpg(2)**2)
     if (qpg_z>SMALL) then 
      ! === Analytic expression ===
      call CALCK0(qpg_z *rcut_,k0,1)
      call CALJY1(qpg_xy*rcut_,j1,0)
      call CALJY0(qpg_xy*rcut_,j0,0)
      call CALCK1(qpg_z *rcut_,k1,1)
      vccut(iq,ig)=(four_pi/qpg2)*(one+rcut_*qpg_xy*j1*k0-qpg_z*rcut_*j0*k1) 
     else 
      if (qpg_xy>SMALL) then 
       ! === Integrate r*Jo(G_xy r)log(r) from 0 up to rcut_  ===
       call quadrature(F5,zero,rcut_,qopt_,quad,ntrial_,accuracy_,npts_)
       vccut(ig,iq)=-four_pi*quad
      else 
       ! === Analytic expression ===
       vccut(ig,iq)=-pi*rcut_**2*(two*LOG(rcut_)-one)
      end if 
     end if 
    end do !ig
   end do !iq 
  end if !finite/infinite

 CASE DEFAULT
  write(msg,'(4a,i3)')ch10,&
&  ' cutoff_cylinder: BUG - ',ch10,&
&  ' wrong value for method: ',method 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT 
 !
 ! === Collect vccut on each node ===
 write(*,*)rank,' completed'
 call leave_test(MPI_enreg)
 call xsum_mpi(vccut(:,:),spaceComm,ierr)
 !write(*,*)MAXVAL(vccut),MINVAL(vccut)

 write(msg,'(a)')' Coulombian tables done '
 call wrtout(std_out,msg,'COLL') 
 deallocate(qcart)

end subroutine cutoff_cylinder
!!***
