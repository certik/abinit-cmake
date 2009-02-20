!{\src2tex{textfont=tt}}
!!****f* ABINIT/cggfft
!! NAME
!! cggfft
!!
!! FUNCTION
!! Calculate the FFT index for all possible G-Gpr
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npwvec=number of plane waves
!!  ngfft1,ngfft2,ngfft3=FFT grid
!!  gvec(3,npwvec)= the reciprical lattice vectors of the PW
!!
!! OUTPUT
!!
!! igfft(npwvec,npwvec)=FFT index for G-Gpr
!!
!! PARENTS
!!      cppm2par,cppm3par,cppm4par
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cggfft(npwvec,ngfft1,ngfft2,ngfft3,gvec,igfft)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ngfft1,ngfft2,ngfft3,npwvec
!arrays
 integer,intent(in) :: gvec(3,npwvec)
 integer,intent(out) :: igfft(npwvec,npwvec)

!Local variables ------------------------------
!scalars
 integer :: gmg01,gmg02,gmg03,ig,igp
!arrays
 integer :: gdiff(3)
!************************************************************************

 do ig=1,npwvec
  do igp=1,npwvec
   gdiff(1)=gvec(1,ig)-gvec(1,igp)
   gdiff(2)=gvec(2,ig)-gvec(2,igp)
   gdiff(3)=gvec(3,ig)-gvec(3,igp)
   gmg01=mod(mod(gdiff(1),ngfft1)+ngfft1,ngfft1)
   gmg02=mod(mod(gdiff(2),ngfft2)+ngfft2,ngfft2)
   gmg03=mod(mod(gdiff(3),ngfft3)+ngfft3,ngfft3)     
   igfft(ig,igp)=1+gmg01+gmg02*ngfft1+gmg03*ngfft1*ngfft2
  end do
 end do

end subroutine cggfft
!!***
