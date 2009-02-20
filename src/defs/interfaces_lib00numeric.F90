!!****m* ABINIT/interfaces_lib00numeric
!! NAME
!! interfaces_lib00numeric
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/lib00numeric
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

module interfaces_lib00numeric

 implicit none

interface
 function ass_leg_pol(l,m,xarg)
  implicit none
  integer, intent(in) :: l
  integer, intent(in) :: m
  double precision :: ass_leg_pol
  double precision, intent(in) :: xarg
 end function ass_leg_pol
end interface

interface
 subroutine coeffs_gausslegint(xmin,xmax,x,weights,n)
  implicit none
  integer :: n
  double precision :: xmax
  double precision :: xmin
  double precision :: weights(n)
  double precision :: x(n)
 end subroutine coeffs_gausslegint
end interface

interface
 subroutine dzgedi(a,lda,n,ipvt,det,work,job)
  implicit none
  integer :: job
  integer :: lda
  integer :: n
  real*8 :: det(2,2)
  real*8 :: a(2,lda,n)
  integer :: ipvt(n)
  real*8 :: work(2,n)
 end subroutine dzgedi
end interface

interface
 subroutine dzgefa(a,lda,n,ipvt,info)
  implicit none
  integer :: info
  integer :: lda
  integer :: n
  real*8 :: a(2,lda,n)
  integer :: ipvt(n)
 end subroutine dzgefa
end interface

interface
 function interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function interp
end interface

interface
 function dinterp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: dinterp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function dinterp
end interface

interface
 function taylor_interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: taylor_interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function taylor_interp
end interface

interface
 function dtaylor_interp(n,z,f,z0,zz)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: dtaylor_interp
  complex(gwpc) :: z0
  complex(gwpc) :: zz
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end function dtaylor_interp
end interface

interface
 subroutine calculate_taylor_c(n,z,f,z0,c)
  use defs_basis
  implicit none
  integer :: n
  complex(gwpc) :: z0
  complex(gwpc) :: c(n)
  complex(gwpc) :: f(n)
  complex(gwpc) :: z(n)
 end subroutine calculate_taylor_c
end interface

interface
 SUBROUTINE INTRPL(L,X,Y,N,U,V,dv,dv2,ideriv)
  implicit none
  integer, parameter :: NQQ=12000
  integer :: L
  integer :: N
  integer :: ideriv
  double precision :: DV(NQQ)
  double precision :: DV2(NQQ)
  double precision :: U(N)
  double precision :: V(N)
  double precision :: X(L)
  double precision :: Y(L)
 end subroutine INTRPL
end interface

interface
 SUBROUTINE CALJY0(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALJY0
end interface

interface
 DOUBLE PRECISION FUNCTION BESJ0(X)
 implicit none
 double precision :: X
end function BESJ0
end interface

interface
 DOUBLE PRECISION FUNCTION BESY0(X)
 implicit none
 double precision :: X
end function BESY0
end interface

interface
 SUBROUTINE CALJY1(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALJY1
end interface

interface
 DOUBLE PRECISION FUNCTION BESJ1(X)
 implicit none
 double precision :: X
end function BESJ1
end interface

interface
 DOUBLE PRECISION FUNCTION BESY1(X)
 implicit none
 double precision :: X
end function BESY1
end interface

interface
 subroutine jacobi(a,n,np,d,v,nrot)
  implicit none
  integer :: n
  integer :: np
  integer :: nrot
  real*8 :: a(np,np)
  real*8 :: d(np)
  real*8 :: v(np,np)
 end subroutine jacobi
end interface

interface
 SUBROUTINE CALCK0(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALCK0
end interface

interface
 DOUBLE PRECISION FUNCTION BESK0(X)
 implicit none
 double precision :: X
end function BESK0
end interface

interface
 DOUBLE PRECISION FUNCTION BESEK0(X)
 implicit none
 double precision :: X
end function BESEK0
end interface

interface
 SUBROUTINE CALCK1(ARG,RESULT,JINT)
  implicit none
  integer :: JINT
  double precision :: ARG
  double precision :: RESULT
 end subroutine CALCK1
end interface

interface
 DOUBLE PRECISION FUNCTION BESK1(X)
 implicit none
 double precision :: X
end function BESK1
end interface

interface
 DOUBLE PRECISION  FUNCTION BESEK1(X)
 implicit none
 double precision :: X
end function BESEK1
end interface

interface
 SUBROUTINE ludcmp(a,n,np,indx,id,info)
  implicit none
  integer :: id
  integer :: info
  integer :: n
  integer :: np
  real*8 :: a(np,np)
  integer :: indx(n)
 end subroutine ludcmp
end interface

interface
 SUBROUTINE lubksb(a,n,np,indx,b)
  implicit none
  integer :: n
  integer :: np
  real*8 :: a(np,np)
  real*8 :: b(n)
  integer :: indx(n)
 end subroutine lubksb
end interface

interface
 subroutine gaussj(a,n,np,b,m,mp,info)
  implicit none
  integer :: info
  integer :: m
  integer :: mp
  integer :: n
  integer :: np
  real*8 :: a(np,np)
  real*8 :: b(np,np)
 end subroutine gaussj
end interface

interface
 subroutine polyn_coeff(n,x,y,coeff)
  implicit none
  integer :: n
  double precision :: coeff(n)
  double precision :: x(n)
  double precision :: y(n)
 end subroutine polyn_coeff
end interface

interface
 subroutine smooth(a,mesh,it)
  implicit none
  integer, intent(in) :: it
  integer, intent(in) :: mesh
  real*8, intent(out) :: a(mesh)
 end subroutine smooth
end interface

interface
 subroutine sort_dp(n,list,iperm,tol)
  implicit none
  integer :: n
  double precision :: tol
  integer :: iperm(n)
  double precision :: list(n)
 end subroutine sort_dp
end interface

interface
 subroutine sort_int(n,list,iperm)
  implicit none
  integer :: n
  integer :: iperm(n)
  integer :: list(n)
 end subroutine sort_int
end interface

interface
 subroutine splfit(arg,derfun,fun,ider,newarg,newfun,numarg,numnew)
  implicit none
  integer :: ider
  integer :: numarg
  integer :: numnew
  double precision :: arg(numarg)
  double precision :: derfun(numnew)
  double precision :: fun(numarg,2)
  double precision :: newarg(numnew)
  double precision :: newfun(numnew)
 end subroutine splfit
end interface

interface
 subroutine spline( t, y, n, ybcbeg, ybcend, ypp, tmp )
  implicit none
  integer :: n
  double precision :: ybcbeg
  double precision :: ybcend
  double precision :: t(n)
  double precision :: tmp(n)
  double precision :: y(n)
  double precision :: ypp(n)
 end subroutine spline
end interface

interface
 subroutine splint (nspline,xspline,yspline,ysplin2,&  
  &  nfit,xfit,yfit)
  implicit none
  integer :: nfit
  integer :: nspline
  double precision :: xfit(nfit)
  double precision :: xspline(nspline)
  double precision :: yfit(nfit)
  double precision :: ysplin2(nspline)
  double precision :: yspline(nspline)
 end subroutine splint
end interface


interface
 function uniformrandom(seed) 
  implicit none
  integer :: seed
  double precision :: uniformrandom
 end function uniformrandom
end interface

interface
 subroutine zgedi(a,lda,n,ipvt,det,work,job)
  implicit none
  integer :: job
  integer :: lda
  integer :: n
  complex*16 :: det(2)
  integer :: ipvt(*)
  complex*16 :: work(*)
  complex*16 :: a(lda,*)
 end subroutine zgedi
end interface

interface
 subroutine zgefa(a,lda,n,ipvt,info)
  implicit none
  integer :: info
  integer :: lda
  integer :: n
  integer :: ipvt(*)
  complex*16 :: a(lda,*)
 end subroutine zgefa
end interface

end module interfaces_lib00numeric
!!***
