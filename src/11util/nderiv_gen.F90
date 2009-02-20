!{\src2tex{textfont=tt}}
!!****f* ABINIT/nderiv_gen
!! NAME
!! nderiv_gen
!!
!! FUNCTION
!! Do corrected first (and -if requested- second) derivation on a given (generalized) grid.
!! This routine interfaces nderiv (derivation on a regular grid).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  func(:)=input function
!!  nder= order of the derivation (1 or 2)
!!  radmesh <type(pawrad_type)>=data containing radial grid information
!!
!! OUTPUT
!!  der(:,nder)=resulting derived function
!!
!! PARENTS
!!      optics_paw,pawdij,pawinit,pawnabla_init,pawxcsph,psp7cc,spline_paw_fncs
!!
!! CHILDREN
!!      deducer0,leave_new,nderiv,wrtout
!!
!! NOTES
!!  Possible mesh types (radmesh%mesh_type)
!!   mesh_type=1 (regular grid): rad(i)=(i-1)*AA
!!   mesh_type=2 (logari. grid): rad(i)=AA*(exp[BB*(i-1)]-1)
!!   mesh_type=3 (logari. grid): rad(i>1)=AA*exp[BB*(i-2)] and rad(1)=0
!!   mesh_type=4 (logari. grid): rad(i)=-AA*ln[1-BB*(i-1)] with BB=1/n
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nderiv_gen(der,func,nder,radmesh)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util, except_this_one => nderiv_gen
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nder
 type(pawrad_type),intent(in) :: radmesh
!arrays
 real(dp),intent(in) :: func(radmesh%mesh_size)
 real(dp),intent(out) :: der(radmesh%mesh_size,nder)

!Local variables-------------------------------
!scalars
 integer :: msz
 character(len=500) :: message
!arrays
 real(dp),allocatable :: func2(:)

! *************************************************************************

 if(nder/=1 .and. nder/=2)then
  write(message, '(a,a,a,a,a,a,i5)' ) ch10,&
&  ' nderiv_den :  BUG -',ch10,&
&  '  first or second derivatives are allowed,',ch10,&
&  '  while the argument nder=',nder
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 msz=radmesh%mesh_size

 if (radmesh%mesh_type==1) then

  call nderiv(radmesh%rstep,func,der(1:msz,1),radmesh%mesh_size,1)
  if (nder==2) call nderiv(radmesh%rstep,func,der(1:msz,2),radmesh%mesh_size,2)

 else if (radmesh%mesh_type==2) then

  call nderiv(radmesh%lstep,func,der(1:msz,1),radmesh%mesh_size,1)
  der(1:msz,1)=der(1:msz,1)/radmesh%radfact(1:msz)
  if (nder==2)then
   call nderiv(radmesh%lstep,func,der(1:msz,2),radmesh%mesh_size,2)
   der(1:msz,2)=(der(1:msz,2)/radmesh%radfact(1:msz)-der(1:msz,1))/radmesh%radfact(1:msz)
  end if

 else if (radmesh%mesh_type==3) then

  call nderiv(radmesh%lstep,func(2:msz),der(2:msz,1),msz-1,1)
  der(2:msz,1)=der(2:msz,1)/radmesh%radfact(2:msz)
  call deducer0(der(:,1),msz,radmesh)
  if (nder==2)then
   call nderiv(radmesh%lstep,func(2:msz),der(2:msz,2),msz-1,2)
   der(2:msz,2)=(der(2:msz,2)/radmesh%radfact(2:msz)-der(2:msz,1))/radmesh%radfact(2:msz)
   call deducer0(der(:,2),msz,radmesh)
  end if

 else if (radmesh%mesh_type==4) then

  call nderiv(radmesh%lstep,func,der(1:msz,1),radmesh%mesh_size,1)
  der(1:msz,1)=der(1:msz,1)/radmesh%radfact(1:msz)
  if (nder==2)then
   call nderiv(radmesh%lstep,func,der(1:msz,2),radmesh%mesh_size,2)
   der(1:msz,2)=der(1:msz,2)/radmesh%radfact(1:msz)**2-der(1:msz,1)/radmesh%rstep
  end if

 end if

end subroutine nderiv_gen
!!***
