!{\src2tex{textfont=tt}}
!!****f* ABINIT/indgrid
!!
!! NAME
!! indgrid
!!
!! FUNCTION
!! Calculate the correspondance between the coarse grid and
!! the fine grid for PAW calculations.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! nfftc=total number of FFt grid=n1*n2*n3 for the coarse grid
!! nfftf=total number of FFt grid=n1*n2*n3 for the fine grid
!! ngfftc(18)=contain all needed information about 3D FFT, for the coarse grid,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!! ngfftf(18)=contain all needed information about 3D FFT, for the fine grid,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! OUTPUT
!! coatofin(nfftc)= index of the points of the coarse grid on the fine grid
!! fintocoa(nfftf)=index of the points of the fine grid on the
!!   coarse grid (=0 if the point of the fine grid does not belong to
!!   the coarse grid).
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine indgrid(coatofin,fintocoa,nfftc,nfftf,ngfftc,ngfftf)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftc,nfftf
!arrays
 integer,intent(in) :: ngfftc(18),ngfftf(18)
 integer,intent(out) :: coatofin(nfftc),fintocoa(nfftf)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,if1,if2,if3,ii,ing,n1c,n1f,n2c,n2f,n3c,n3f,narg1,narg2
!arrays
 integer :: id(3)
 integer,allocatable :: gc(:,:),gf(:,:)

! *************************************************************************
!

 n1c=ngfftc(1);n2c=ngfftc(2);n3c=ngfftc(3)
 n1f=ngfftf(1);n2f=ngfftf(2);n3f=ngfftf(3)

 allocate(gc(3,max(n1c,n2c,n3c)))
 do ii=1,3
  id(ii)=ngfftc(ii)/2+2
  do ing=1,ngfftc(ii)
   gc(ii,ing)=ing-(ing/id(ii))*ngfftc(ii)-1
  end do
 end do

 allocate(gf(3,max(n1f,n2f,n3f)))
 do ii=1,3
  id(ii)=ngfftf(ii)/2+2
  do ing=1,ngfftf(ii)
   gf(ii,ing)=ing-(ing/id(ii))*ngfftf(ii)-1
  end do
 end do

 coatofin=0;fintocoa=0
 do i1=1,n1c
  do if1=1,n1f
   if(gc(1,i1)==gf(1,if1)) then
    do i2=1,n2c
     do if2=1,n2f
      if(gc(2,i2)==gf(2,if2)) then
       do i3=1,n3c
        do if3=1,n3f
         if(gc(3,i3)==gf(3,if3)) then
          narg1=i1+n1c*(i2-1+n2c*(i3-1))
          narg2=if1+n1f*(if2-1+n2f*(if3-1))
          coatofin(narg1)=narg2
          fintocoa(narg2)=narg1
          exit
         end if
        end do
       end do
      end if
     end do
    end do
   end if
  end do
 end do

 deallocate(gf,gc)

end subroutine indgrid
!!***
