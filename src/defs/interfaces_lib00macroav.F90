!!****m* ABINIT/interfaces_lib00macroav
!! NAME
!! interfaces_lib00macroav
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/lib00macroav
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_lib00macroav

 implicit none

interface
 subroutine io
  implicit none
 end subroutine io
end interface

interface
 subroutine iorho( task, fname, cell, mesh, nsm, maxp, nspin,&  
  f, found )
  implicit none
  integer :: maxp
  integer :: nsm
  integer :: nspin
  character*(*) :: fname
  logical :: found
  character*(*) :: task
  double precision :: cell(3,3)
  integer :: mesh(3)
  real :: f(maxp,nspin)
 end subroutine iorho
end interface

interface
 SUBROUTINE FOUR1(DATA,NN,ISIGN)
  implicit none
  integer :: ISIGN
  integer :: NN
  double precision :: DATA(2*NN)
 end subroutine FOUR1
end interface

interface
 SUBROUTINE POLINT(XA,YA,N,X,Y,DY) 
  implicit none
  integer :: N
  double precision :: DY
  double precision :: X
  double precision :: Y
  double precision :: XA(N)
  double precision :: YA(N)
 end subroutine POLINT
end interface

interface
 SUBROUTINE MACROAV_SPLINE(DX,Y,N,YP1,YPN,Y2) 
  implicit none
  integer :: N
  double precision :: DX
  double precision :: YP1
  double precision :: YPN
  double precision :: Y(N)
  double precision :: Y2(N)
 end subroutine MACROAV_SPLINE
end interface

interface
 SUBROUTINE MACROAV_SPLINT(DX,YA,Y2A,N,X,Y,DYDX) 
  implicit none
  integer :: N
  double precision :: DX
  double precision :: DYDX
  double precision :: X
  double precision :: Y
  double precision :: Y2A(N)
  double precision :: YA(N)
 end subroutine MACROAV_SPLINT
end interface

interface
 CHARACTER(LEN=26) FUNCTION PASTE( STR1, STR2 )
 implicit none
 character(len=*) :: STR1
 character(len=*) :: STR2
end function PASTE
end interface

interface
 CHARACTER(LEN=26) FUNCTION PASTEB( STR1, STR2 )
 implicit none
 character(len=*) :: STR1
 character(len=*) :: STR2
end function PASTEB
end interface

interface
 DOUBLE PRECISION FUNCTION SURPLA( C )
 implicit none
 double precision :: C(3,3)
end function SURPLA
end interface

interface
 subroutine thetaft(n,L,lav,ft)
  implicit none
  integer :: n
  real*8 :: L
  real*8 :: lav
  real*8 :: ft(2*n)
 end subroutine thetaft
end interface

interface
 DOUBLE PRECISION FUNCTION VOLCEL( C )
 implicit none
 double precision :: C(3,3)
end function VOLCEL
end interface

end module interfaces_lib00macroav
!!***
