!{\src2tex{textfont=tt}}
!!****f* ABINIT/phim
!! NAME
!! phim
!!
!! FUNCTION
!! Computes Phi_m[theta]=Sqrt[2] cos[m theta],      if m>0
!!                       Sqrt[2] sin[Abs(m) theta], if m<0
!!                       1                        , if m=0
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (NH, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  costeta= cos(theta)  (theta= input angle)
!!  mm = index m
!!  sinteta= sin(theta)  (theta= input angle)
!!
!! OUTPUT
!!  phim= Phi_m(theta) (see above)
!!
!! NOTES
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!
!! PARENTS
!!     setsymrhoij
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

function phim(costheta,sintheta,mm)

 use defs_basis

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mm
 real(dp) :: phim
 real(dp),intent(in) :: costheta,sintheta

!Local variables-------------------------------
!scalars
 real(dp) :: sqr2

! *********************************************************************

 sqr2=sqrt(2._dp)
 if (mm==0)  phim=1._dp
 if (mm==1)  phim=sqr2*costheta
 if (mm==-1) phim=sqr2*sintheta
 if (mm==2)  phim=sqr2*(costheta*costheta-sintheta*sintheta)
 if (mm==-2) phim=sqr2*2._dp*sintheta*costheta
 if (mm==3)  phim=sqr2*&
& (costheta*(costheta*costheta-sintheta*sintheta)&
& -sintheta*2._dp*sintheta*costheta)
 if (mm==-3) phim=sqr2*&
& (sintheta*(costheta*costheta-sintheta*sintheta)&
& +costheta*2._dp*sintheta*costheta)

 end function phim

!!***
