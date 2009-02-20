!{\src2tex{textfont=tt}}
!!****f* ABINIT/printxsf
!! NAME
!! printxsf
!!
!! FUNCTION
!! Write a generic array in the XSF format (XCrysden format)
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! basis(3,3) = basis vectors of the direct real lattice or of the reciprocal lattice (fortran convention)
!!              (Bohr units if option=0, Bohr^-1 if option=1, see below)
!! option  =  0  for a plot in real space
!!            1  for a plot in reciprocal space
!! nunit   = unit number of the output file
!! n1=grid size along x
!! n2=grid size along y
!! n3=grid size along z
!! origin(3) = origin of the grid
!! datagrid(n1*n2*n3) = datagrid values stored using the fortran convention
!!
!! OUTPUT
!! Only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!! TODO: Treat the real space case and include information on the lattice
!!
!! PARENTS
!!      mknesting
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine printxsf(n1,n2,n3,datagrid,basis,origin,nunit,option)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: n1,n2,n3,nunit,option
!arrays
 real(dp),intent(in) :: basis(3,3),datagrid(n1*n2*n3),origin(3)

!Local variables-------------------------------
!scalars
 integer :: iost,iy,iz,nslice,nsym,ycount
 real(dp) :: fact
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(6,*)' printxsf : enter '
!ENDDEBUG

 if(option/=0 .and. option/=1 )then
  write(message,'(a,a,a,a,a,a,i6)') ch10,          &
&  ' printxsf: ERROR -',ch10,                      &
&  '  The argument option should be 1 or 2,',ch10, &
&  '  however, option= ',option
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!conversion between ABINIT default units and XCrysden units
 fact=Bohr_Ang
 if (option ==1) fact=1./fact  !since we are in reciprocal space

!TODO insert crystalline structure and dummy atoms in case of reciprocal space
!write(nunit,'(a)')' CRYSTAL'

!if (option == 1) then
!write(nunit,'(a)')' # these are primitive lattice vectors (in Angstroms)'
!write(nunit,'(a)')
!else
!write(nunit,'(a)')' # these are primitive reciprocal lattice vectors (in Angstroms -1)'
!write(nunit,'(a)')
!end if

!write(nunit,'(a)')' PRIMVEC'
!write(nunit,*)basis(:,1)*fact
!write(nunit,*)basis(:,2)*fact
!write(nunit,*)basis(:,3)*fact

!write(nunit,'(a)')'DIM-GROUP'
!write(nunit,'(a)')' 3  1'
!write(nunit,'(a)')' PRIMCOORD'
!write(nunit,'(a)')' 1  1'
!write(nunit,'(a)')'       13    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00'

 write(nunit,'(a)')' BEGIN_BLOCK_DATAGRID3D'
 write(nunit,'(a)')' datagrid'
 write(nunit,'(a)')' DATAGRID_3D_DENSITY'
!NOTE: XCrysden uses aperiodical data grid
 write(nunit,*)n1+1,n2+1,n3+1
 write(nunit,*)origin
 write(nunit,*)basis(:,1)*fact
 write(nunit,*)basis(:,2)*fact
 write(nunit,*)basis(:,3)*fact

 nslice=1
 do iz=1,n3
  do iy=1,n2
   write(nunit,'(8es16.8)') datagrid(1+n1*(nslice-1):n1+n1*(nslice-1)),datagrid(1+n1*(nslice-1))
!  DEBUG
!  write(*,*)1+n1*(nslice-1),n1+n1*(nslice-1)
!  ENDDEBUG
   nslice=nslice+1
  end do
  nsym=nslice-n2
  write (nunit,'(8es16.8)') datagrid(1+n1*(nsym-1):n1+n1*(nsym-1)),datagrid(1+n1*(nsym-1))
! DEBUG
! write(*,*)1+n1*(nsym-1),n1+n1*(nsym-1)
! write(*,*)' done xy plane at z = ',nslice-n2
! ENDDEBUG
 end do

!Now write upper plane
 nslice=1
 do iy=1,n2
  write (nunit,'(8es16.8)') datagrid(1+n1*(nslice-1):n1+n1*(nslice-1)),datagrid(1+n1*(nslice-1))
  write(*,*)1+n1*(nslice-1),n1+n1*(nslice-1)
  nslice=nslice+1
 end do

 nsym=nslice-n2
 write (nunit,'(8es16.8)') datagrid(1+n1*(nsym-1):n1+n1*(nsym-1)),datagrid(1+n1*(nsym-1))
 write(*,*)1+n1*(nsym-1),n1+n1*(nsym-1)
 write(*,*)' done upper xy plane '

 write (nunit,'(a)')' END_DATAGRID_3D'
 write (nunit,'(a)')' END_BLOCK_DATAGRID3D'

!DEBUG
!write(6,*)' printxsf : exit'
!stop
!ENDDEBUG

end subroutine printxsf

!!***
