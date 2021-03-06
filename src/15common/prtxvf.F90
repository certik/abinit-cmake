!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtxvf
!!
!! NAME
!! prtxvf
!!
!! FUNCTION
!! Print the values of xcart, vel, and fcart to unit iout.
!! Also compute and print max and rms forces.

!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! fcart(3,natom)=forces (hartree/bohr)
!! iatfix(3,natom)=1 for frozen or fixed atom along specified direction, else 0
!! iout=unit number for printing
!! natom=number of atoms in unit cell.
!! prtvel=1 to print velocities, else do not print them
!! vel(3,natom)=velocities
!! xcart(3,natom)=cartesian coordinates (bohr)
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      brdmin,constrf,delocint,moldyn,move
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prtxvf(fcart,iatfix,iout,natom,prtvel,vel,xcart)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,natom,prtvel
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: fcart(3,natom),vel(3,natom),xcart(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu,unfixd
 real(dp) :: fmax,frms
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(6,*)' prtxvf : enter '
!stop
!ENDDEBUG

 write(message, '(a)' ) ' Cartesian coordinates (bohr)'
 call wrtout(iout,message,'COLL')
 do iatom=1,natom
  write(message, '(1p,3e22.14)' )xcart(:,iatom)
  call wrtout(iout,message,'COLL')
 end do

 if (prtvel == 1) then
  write(message, '(a)' ) ' Velocities (bohr/(atomic time unit))'
  call wrtout(iout,message,'COLL')
  do iatom=1,natom
   write(message, '(1p,3e22.14)' ) vel(:,iatom)
   call wrtout(iout,message,'COLL')
  end do
 end if

!Compute max |f| and rms f, EXCLUDING the components determined by iatfix

 fmax=0.0_dp
 frms=0.0_dp
 unfixd=0
 do iatom=1,natom
  do mu=1,3
   if (iatfix(mu,iatom) /= 1) then
    unfixd=unfixd+1
    frms=frms+fcart(mu,iatom)**2
    fmax=max(fmax,abs(fcart(mu,iatom)))
   end if
  end do
 end do
 if ( unfixd /= 0 ) frms=sqrt(frms/dble(unfixd))

 write(message, '(a,1p,2e12.5,a)' ) &
& ' Cartesian forces (hart/bohr); max,rms=',fmax,frms,&
& ' (free atoms)'
 call wrtout(iout,message,'COLL')
 do iatom=1,natom
  write(message, '(1p,3e22.14)' )fcart(:,iatom)
  call wrtout(iout,message,'COLL')
 end do

 write(message, '(a)' ) ' '
 call wrtout(iout,message,'COLL')

!DEBUG
!write(6,*)' prtxvf : exit '
!stop
!ENDDEBUG

end subroutine prtxvf
!!***
