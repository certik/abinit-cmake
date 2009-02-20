!{\src2tex{textfont=tt}}
!!****f* ABINIT/mlwfovlp_radial
!! NAME
!! mlwfovlp_radial
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (trangel,drh)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  alpha= Z/a = zona
!!  lmax= maximum value of l
!!  rvalue= integer defining the choice for radial functions R(r).
!!   It can take values from 1-3. 
!!   It is associted to the radial part of the hydrogenic Schrodinger equation for l=0,
!!   See the manual of Wannier90 for more information. (www.wannier.org)
!!  xx= scalar number used to calculate the spherical bessel function. J_il(xx)
!!
!! OUTPUT
!!  mlwfovlp_radial= radial part for initial projections used to construct MLWF
!!
!! SIDE EFFECTS
!!  None
!!
!! NOTES
!!  Calculates the radial part of the initial functions given as an initial
!!  guess by the user to construct the MLWF.
!!  
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mlwfovlp_radial(alpha,lmax,lmax2,radial,rvalue,xx)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_14occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmax,lmax2,rvalue
 real(dp),intent(in) :: alpha,xx
!arrays
 real(dp),intent(out) :: radial(lmax2)

!Local variables
!(2*l+1)!!
!scalars
 integer :: il,ir,ll,lm,mesh,mm
 real(dp),parameter :: dx=0.015d0,rmax=10.d0,xmin=0.d0
 real(dp) :: aa,ftmp,gauss,rtmp,sum,x
 character(len=500) :: message
!arrays
 real(dp),save :: dblefact(4)=(/1_dp,3_dp,15_dp,105_dp/)
 real(dp),allocatable :: aux(:),bes(:),cosr(:),func_r(:),r(:),rad_int(:)
 real(dp),allocatable :: sinr(:)

! *************************************************************************
 
!DEBUG
!write (std_out,*) ' mlwfovlp_radial : enter'
!ENDDEBUG

!Radial functions in the form of hydrogenic orbitals as defined in the 
!wannier90 manual.
 if(( rvalue > 0 ).and.(rvalue < 4)) then

! mesh
  mesh= nint((rmax - xmin ) / dx + 1)
  allocate ( bes(mesh), func_r(mesh), r(mesh),rad_int(mesh))
  allocate ( aux(mesh),cosr(mesh),sinr(mesh))
  do ir=1, mesh
   x=xmin+DBLE(ir-1)*dx
   r(ir)=x
  end do   !ir

! radial functions shown in table 3.3 of wannier90 manual
  if (rvalue==1) func_r(:) = 2.d0 * alpha**(3.d0/2.d0) * exp(-alpha*r(:))
  if (rvalue==2) func_r(:) = 1.d0/(2.d0*sqrt(2.d0))*alpha**(3.d0/2.d0) *&
&  (2.d0 - alpha*r(:))*exp(-alpha*r(:)/2.d0)
  if (rvalue==3) func_r(:) = sqrt(4.d0/27.d0)*alpha**(3.d0/2.d0)&
&  * (1.d0 - 2.d0*alpha*r(:)/3.d0 + 2.d0*alpha**2*r(:)**2/27.d0)&
&  * exp(-alpha * r(:)/3.d0)

! compute spherical bessel functions
  cosr(:)=cos(xx*r(:))
  sinr(:)=sin(xx*r(:))
  lm=0
  do ll=0,lmax
   call besjm(xx,bes,cosr,ll,mesh,sinr,r)
   aux(:)=bes(:)*func_r(:)*r(:)
!  do ir=1,mesh
!  write(310,*) r(ir),bes(ir)
!  end do
   call simpson_int(mesh,dx,aux,rad_int)
   sum=0.d0
   do ir=1,mesh
    sum=sum+rad_int(ir)
   end do
   rtmp=sum/mesh
   do mm=-ll,ll
    lm=lm+1
    radial(lm)=rtmp
   end do !mm
  end do !ll
  deallocate(bes,func_r,r,aux,rad_int,cosr,sinr)

! Radial part in the form of Gaussian functions of a given width 
! Taken by code of made by drh.
 elseif ( rvalue == 4) then
  aa=1._dp/alpha
  gauss=exp(-0.25_dp*(aa*xx)**2)
  lm=0
  do ll=0,lmax
   ftmp=(0.5_dp*pi)**(0.25_dp)*aa*sqrt(aa/dblefact(ll+1))*(aa*xx)**ll*gauss
   do mm=-ll,ll
    lm=lm+1
    radial(lm)=ftmp
   end do
  end do
 else ! rvalue < 0 of rvalue > 4
  write(message,'(4a,i6,5a)') ch10,&
&  ' mlwfovlp_radial: BUG -',ch10,&
&  '  Radial function r=,',rvalue,ch10,&
&  '  is not defined',ch10,&
&  '  Modify .win file',ch10
  call wrtout(std_out,message,'COLL')
  call leave_new('COLL')
 end if !rvalue

!write(309,*)'radial ',mlwfovlp_radial



!DEBUG
!write (std_out,*) ' mlwfovlp_radial : exit'
!stop
!ENDDEBUG

end subroutine mlwfovlp_radial
!!***
