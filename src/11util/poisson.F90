!{\src2tex{textfont=tt}}
!!****f* ABINIT/poisson
!! NAME
!! poisson
!!
!! FUNCTION
!!  Solve poisson equation for angularly dependent charge
!!  distribution of angular momentum l
!!  Densities and potentials are given on a (generalized) radial grid
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  den(radmesh%mesh_size)= electron density * (4*pi*r**2) appropriate for l
!!  ll= l quantum number
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! OUTPUT
!!  qq= lth moment of the charge
!!  rv(radmesh%mesh_size)= electrostatic potential * r in (Hartree*Bohr) units
!!          where v(r)=\frac{1}{2l+1}(\frac{int[(r''^(l+2))g(r'')dr'']} {r^(l+1)}
!!                                   +(r^l) int[r''^(1-l)g(r'')dr''])
!!
!! PARENTS
!!      pawdenpot,pawinit,psp7in
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine poisson(den,ll,qq,radmesh,rv)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util, except_this_one => poisson
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ll
 real(dp),intent(out) :: qq
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: den(radmesh%mesh_size)
 real(dp),intent(out) :: rv(radmesh%mesh_size)

!Local variables ---------------------------------------
!scalars
 integer :: ir,jr,mm,nn
 real(dp) :: angm,hh,rr
 character(len=500) :: message
!arrays
 real(dp) :: aa(radmesh%mesh_size),bb(radmesh%mesh_size),cc(radmesh%mesh_size)
 real(dp) :: dd(radmesh%mesh_size+1),ee(4)
 real(dp),allocatable :: radl(:),radl1(:)

! ***************************************************************************

!==============================================
!UNIFORM GRID - NUMEROV ALGORITHM
!==============================================
 if (radmesh%mesh_type==1) then
  hh=radmesh%rstep
  nn=radmesh%int_meshsz-1
  do ir=1,nn
   aa(ir)=two*hh*den(ir+1)/(ir)
   bb(ir)=den(ir+1)*((ir*hh)**ll)
  end do
  dd(1)=zero
  dd(2:nn+1)=bb(1:nn)
  call simp_gen(qq,dd,radmesh)
  qq=qq/(2*ll+1)
  rv(1)=aa(1)+0.1_dp*aa(2)
  do ir=2,nn-1
   rv(ir)=aa(ir)+0.1_dp*(aa(ir+1)+aa(ir-1))
  end do
  rv(nn)=aa(nn)+0.1_dp*aa(nn-1)
  angm=dble(ll*(ll+1))
  rr=(nn+1)*hh
  rv(nn)=rv(nn)+(2.4_dp-0.2_dp*angm/((nn+1)**2))*qq/(rr**ll)
  do ir=1,nn
   aa(ir)=angm/(ir*ir)
   bb(ir)=2.4_dp+aa(ir)
  end do
  do ir=1,nn-1
   cc(ir)=-1.2_dp+0.1_dp*aa(ir+1)
  end do
  do ir=nn,2,-1
   aa(ir)=-1.2_dp+0.1_dp*aa(ir-1)
  end do
  if (nn.eq.1) then
   rv(2)=rv(1)/bb(1)
   rv(1)=zero
  else
   do ir=2,nn
    bb(ir)=bb(ir)-aa(ir)*cc(ir-1)/bb(ir-1)
   end do
   rv(1)=rv(1)/bb(1)
   do ir=2,nn
    rv(ir)=(rv(ir)-aa(ir)*rv(ir-1))/bb(ir)
   end do
   do ir=nn-1,1,-1
    rv(ir)=rv(ir)-cc(ir)*rv(ir+1)/bb(ir)
   end do
   do ir=nn+1,2,-1
    rv(ir)=rv(ir-1)
   end do
   rv(1)=zero
   if (nn+1<radmesh%mesh_size) rv(nn+2:radmesh%mesh_size)=zero
  end if
  rv(:)=half*rv(:)

! ==============================================
! ANY OTHER GRID - SIMPSON ALGORITHM
! ==============================================
 else
  nn=radmesh%mesh_size;hh=third*radmesh%stepint
  do while (abs(den(nn))<tol16.and.nn>radmesh%int_meshsz)
   nn=nn-1
  end do
  mm=nn;if (radmesh%mesh_type==3) mm=mm-1
  allocate(radl(nn),radl1(nn))
  do jr=nn,2,-1
   ir=nn-jr+1
   radl(jr) =radmesh%rad(jr)**ll
   radl1(jr)=radmesh%rad(jr)*radl(jr)
   aa(ir)=den(jr)*radmesh%radfact(jr)*radl(jr)
   bb(ir)=den(jr)*radmesh%radfact(jr)/radl1(jr)
  end do
  radl(1)=zero;radl1(1)=zero
  ee(2)=aa(nn-1);ee(3)=aa(nn-2);ee(4)=aa(nn-3)
  call deducer0(ee,4,radmesh)
  aa(nn)=ee(1)
  ee(2)=bb(nn-1);ee(3)=bb(nn-2);ee(4)=bb(nn-3)
  call deducer0(ee,4,radmesh)
  bb(nn)=ee(1)
  cc(1)=zero;dd(1)=zero
  do ir=3,mm,2
   cc(ir)  =cc(ir-2)+hh*(aa(ir-2)+four*aa(ir-1)+aa(ir))
   cc(ir-1)=cc(ir-2)+hh*(1.25_dp*aa(ir-2)+two*aa(ir-1)-quarter*aa(ir))
   dd(ir)  =dd(ir-2)+hh*(bb(ir-2)+four*bb(ir-1)+bb(ir))
   dd(ir-1)=dd(ir-2)+hh*(1.25_dp*bb(ir-2)+two*bb(ir-1)-quarter*bb(ir))
  end do
  if (mod(mm,2)==0) then
   cc(mm)=cc(mm-2)+hh*(aa(mm-2)+four*aa(mm-1)+aa(mm))
   dd(mm)=dd(mm-2)+hh*(bb(mm-2)+four*bb(mm-1)+bb(mm))
  end if
  if (mm<nn) then
   cc(nn)=cc(mm)+half*(aa(mm)+aa(nn))*(radmesh%rad(1+nn-mm)-radmesh%rad(1))
   dd(nn)=dd(mm)+half*(bb(mm)+bb(nn))*(radmesh%rad(1+nn-mm)-radmesh%rad(1))
  end if
  rv(1)=zero
  do ir=2,nn
   jr=nn-ir+1
   rv(ir)=(dd(jr)*radl1(ir)+(cc(nn)-cc(jr))/radl(ir))/(two*ll+one)
  end do
  if (nn<radmesh%mesh_size) rv(nn+1:radmesh%mesh_size)=rv(nn)
  qq=cc(nn)/(two*ll+one)
  deallocate(radl,radl1)
 end if

end subroutine poisson
!!***
