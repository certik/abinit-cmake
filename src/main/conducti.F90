!{\src2tex{textfont=tt}}
!!****p* ABINIT/conducti
!! NAME
!! conducti
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula.
!!
!! COPYRIGHT
!! Copyright (C) 2006-2008 ABINIT group (FJ,SMazevet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!!
!! CHILDREN
!!      conducti_nc, conducti_paw
!!
!! PARENTS
!!
!! CHILDREN
!!      conducti_nc,conducti_paw,linear_optics_paw
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program conducti

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15common
!End of the abilint section

 implicit none

!Arguments -----------------------------------

!Local variables-------------------------------
!no_abirules
 integer :: incpaw
 character(len=fnlen) :: filnam
 ! here only since needed for call to rwwf
 type(MPI_type) :: mpi_enreg

! *********************************************************************************
!BEGIN EXECUTABLE SECTION

!Read data file name
 write(6,'(a)')' Please, give the name of the data file ...'
 read(5, '(a)')filnam
 write(6,'(a)')' The name of the data file is :',filnam
 open(15,file=filnam,form='formatted')
 rewind(15)
 read(15,*) incpaw
 close(15)
 if (incpaw==1) then
  call conducti_nc(filnam,mpi_enreg)
 elseif (incpaw==2) then
  call conducti_paw(filnam,mpi_enreg)
 elseif (incpaw==3) then
  call linear_optics_paw(filnam,mpi_enreg)
 end if

 end program conducti
!!***
