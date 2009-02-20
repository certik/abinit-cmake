!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_FFT_rotation
!! NAME
!! setup_FFT_rotation
!!
!! FUNCTION
!! For each real space point, r, of the FFT mesh, calculate the FFT index no. of  $R^{-1}(r-\tau)$.
!! R is a symmetry operation in real space, $\tau$ is the associated fractional translation.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MT, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nfft=number of points for this processor.
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nsym=number of symmetry operations.
!! symrec(3,3,nsym)=symmetry operations in reciprocal space.
!! tnons(3,nsym)=fractional translations.
!!
!! OUTPUT
!!  irottb(ngfft(1)*ngfft(2)*ngfft(3),nsym)= contains the index in the FFT array of (R^{-1})(r-t), 
!!   where R is one of the nsym symmetry operations in real space.
!!
!! NOTES: 
!!  The evaluation of the rotated point $R^{-1}(r-\tau)$ is done using real arithmetic.
!!  As a consequence, if the FFT mesh does not respect the symmetry properties 
!!  of the system, the array irottb will report the index of the FFT point which 
!!  is the closest one to $R^{-1}(r-\tau)$. This might lead to some inaccuracies, in particular
!!  during the calculation of degenerate states.
!!
!! PARENTS
!!      
!!
!! CHILDREN
!!      dosym,dosymr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setup_FFT_rotation(nsym,symrec,tnons,nfft,ngfft,irottb)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nsym
!arrays
 integer,intent(in) :: ngfft(18),symrec(3,3,nsym)
 integer,intent(out) :: irottb(nfft,nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables ------------------------------
!scalars
 integer :: iinv,ir,isym,ix,iy,iz,jx,jy,jz,ngfft1,ngfft2,ngfft3
 character(len=500) :: msg
!arrays
 integer :: rm1(3,3,nsym)
 real(dp) :: rbase(3),rrot(3)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' setup_FFT_rotation : setting up FFT rotation table'
 call wrtout(std_out,msg,'PERS')
#endif
 
 if (nfft/=ngfft(1)*ngfft(2)*ngfft(3)) then 
  write(msg,'(4a)')ch10,&
&  ' setup_FFT_rotation : BUG : ',ch10,&
&  ' FFT parallelism not allowed '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 !
 ! === Precalculate operations R^-1 in real space ===
 ! TODO should pass symrel instead of symrec.
 do isym=1,nsym
  rm1(:,:,isym)=TRANSPOSE(symrec(:,:,isym))
 end do

 ngfft1=ngfft(1) ; ngfft2=ngfft(2) ; ngfft3=ngfft(3)
 do ix=0,ngfft1-1
  do iy=0,ngfft2-1
   do iz=0,ngfft3-1
    rbase(1)=DBLE(ix)/ngfft1
    rbase(2)=DBLE(iy)/ngfft2
    rbase(3)=DBLE(iz)/ngfft3
    ir=1+ix+iy*ngfft1+iz*ngfft1*ngfft2
    do isym=1,nsym
     ! === Form R^-1 (rbase-\tau) ====
     rrot=MATMUL(rm1(:,:,isym),rbase(:)-tnons(:,isym))
     jx=NINT(rrot(1)*ngfft1)
     jy=NINT(rrot(2)*ngfft2)
     jz=NINT(rrot(3)*ngfft3)
     jx=MODULO(jx,ngfft1)
     jy=MODULO(jy,ngfft2)
     jz=MODULO(jz,ngfft3)
     irottb(ir,isym)=1+jx+jy*ngfft1+jz*ngfft1*ngfft2
    end do 
   end do 
  end do 
 end do

#if defined DEBUG_MODE
 write(msg,'(a)')' FFT rotation table set up.'
 call wrtout(std_out,msg,'COLL')
#endif

end subroutine setup_FFT_rotation
!!***

!!****if* ABINIT/FFT_rotations
!! NAME
!! FFT_rotations
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine FFT_rotations(Cryst,ngfft,irottb)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: ngfft(18)
 integer,intent(out) :: irottb(ngfft(1)*ngfft(2)*ngfft(3),Cryst%nsym)

!Local variables-------------------------------
!scalars
 integer :: ir1,isym,ix,iy,iz,jx,jy,jz,ngfft1,ngfft2,ngfft3
 character(len=500) :: msg
!arrays
 integer :: Rm1(3,3,Cryst%nsym),r1_FFT(3),red2fft(3,3)
 integer,pointer :: symrel(:,:,:)
 real(dp) :: Rm1_FFT(3,3,Cryst%nsym),err(3,Cryst%nsym),fft2red(3,3),r2_FFT(3)
 real(dp) :: tnons_FFT(3,Cryst%nsym)
 real(dp),pointer :: tnons(:,:)

! *************************************************************************

 ! === Precalculate R^-1 and fractional translation in FFT coordinates ===
 ngfft1=ngfft(1) ; ngfft2=ngfft(2) ; ngfft3=ngfft(3)
 red2fft=RESHAPE((/ngfft1,0,0,0,ngfft2,0,0,0,ngfft3/),(/3,3/))
 fft2red=RESHAPE((/(one/ngfft1),zero,zero,zero,(one/ngfft2),zero,zero,zero,(one/ngfft3)/),(/3,3/))

 symrel => Cryst%symrel
 tnons  => Cryst%tnons
 !
 ! === For a fully compatible mesh, each Rm1_FFT should be integer ===
 do isym=1,Cryst%nsym
  call mati3inv(symrel(:,:,isym),Rm1(:,:,isym))
  Rm1(:,:,isym)=TRANSPOSE(Rm1(:,:,isym))
  Rm1_FFT(:,:,isym)=MATMUL(MATMUL(red2fft,Rm1(:,:,isym)),fft2red)
  tnons_FFT(:,isym)=MATMUL(red2fft,tnons(:,isym))
 end do

 err(:,:)=smallest_real 
 do ix=0,ngfft1-1
  R1_FFT(1)=ix
  do iy=0,ngfft2-1
   R1_FFT(2)=iy
   do iz=0,ngfft3-1
    R1_FFT(3)=iz
    ir1=1+ix+iy*ngfft1+iz*ngfft1*ngfft2
    do isym=1,Cryst%nsym
     ! === Form R^-1 (r-\tau) in the FFT basis ===
     R2_FFT(:)=MATMUL(Rm1_FFT(:,:,isym),R1_FFT(:)-tnons_FFT(:,isym))
     jx=NINT(R2_FFT(1)) ; err(1,isym)=MAX(err(1,isym),ABS(R2_FFT(1)-jx))
     jy=NINT(R2_FFT(2)) ; err(2,isym)=MAX(err(2,isym),ABS(R2_FFT(2)-jy))
     jz=NINT(R2_FFT(3)) ; err(3,isym)=MAX(err(3,isym),ABS(R2_FFT(3)-jz))
     jx=MODULO(jx,ngfft1)
     jy=MODULO(jy,ngfft2)
     jz=MODULO(jz,ngfft3)
     irottb(ir1,isym)=1+jx+jy*ngfft1+jz*ngfft1*ngfft2
    end do 
   end do 
  end do 
 end do

 do isym=1,Cryst%nsym 
  if (ANY(err(:,isym)>tol6)) then
   write(*,*)' WARNING- symmetry no. ',isym,' not compatible with FFT grid, FFT error ',err(:,isym)
  end if
 end do

end subroutine FFT_rotations
!!***
