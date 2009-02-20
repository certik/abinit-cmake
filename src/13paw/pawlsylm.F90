!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawlsylm
!! NAME
!! pawlsylm
!!
!! FUNCTION
!! Compute the LS operator in the real spherical harmonics basis
!! ls_ylm(ilm1,ilm2,ispin)= <sigma, S_lm1| L.S |S_lm2, sigma_prime>
!!   ilm,1m2=(l,m1,m2) with -l<=m1<=l, -l<=m2<=l and 0<l<=lmax
!!   ispin=(sigma,sigma_prime) 1=(up,up), 2=(up,dn), 3=(dn,up), 4=(dn,dn)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!
!! OUTPUT
!!  pawang%ls_ylm(2,l_max**2*(l_max**2+1)/2,2)=LS operator in the real spherical harmonics basis
!!        ls_ylm(:,:,1)=<up, S_lm1| L.S |S_lm2, up>
!!        ls_ylm(:,:,2)=<up, S_lm1| L.S |S_lm2, down>
!!        One can deduce:
!!        <down, S_lm1| L.S |S_lm2, down>=-<up, S_lm1| L.S |S_lm2, up>
!!        <down, S_lm1| L.S |S_lm2, up>  =-Conjg[<up, S_lm1| L.S |S_lm2, down>]
!!        Also, only ilm1<=ilm2 terms are stored, because:
!!         <sigma, S_lm1| L.S |S_lm2, sigma_prime>=-<sigma_prime, S_lm1| L.S |S_lm2, sigma>
!!
!! PARENTS
!!      pawinit
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawlsylm(pawang)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 type(pawang_type),intent(inout) :: pawang

!Local variables ---------------------------------------
!scalars
 integer :: ii,ilm,im,j0lm,jj,jlm,jm,klm,l_max,ll,lm0,mm
 real(dp),parameter :: invsqrt2=one/sqrt2
 real(dp) :: onem
 character(len=500) :: message
!arrays
 complex(dpc) :: tmp(2)
 complex(dpc),allocatable :: ls_cplx(:,:,:),slm2ylm(:,:)

! *************************************************************************

 if (pawang%use_ls_ylm==0) then
  write(message, '(4a)' ) ch10,&
&  ' pawlsylm :  BUG -',ch10,&
&  '   ls_ylm pointer is not allocated !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Initialization
 pawang%ls_ylm=zero
 l_max=pawang%l_max-1

!Nothing to do if lmax=0
 if (l_max<=0) return

!Loop on l quantum number
 do ll=1,l_max

! Transformation matrixes: real->complex spherical harmonics
  allocate(slm2ylm(2*ll+1,2*ll+1));slm2ylm=czero
  do im=1,2*ll+1
   mm=im-ll-1;jm=-mm+ll+1
   onem=dble((-1)**mm)
   if (mm> 0) then
    slm2ylm(im,im)= cmplx(onem*invsqrt2,zero,kind=dp)
    slm2ylm(jm,im)= cmplx(invsqrt2,     zero,kind=dp)
   end if
   if (mm==0) then
    slm2ylm(im,im)=cone
   end if
   if (mm< 0) then
    slm2ylm(im,im)= cmplx(zero,     invsqrt2,kind=dp)
    slm2ylm(jm,im)=-cmplx(zero,onem*invsqrt2,kind=dp)
   end if
  end do

! Compute <sigma, Y_lm1|L.S|Y_lm2, sigma_prime> (Y_lm=complex spherical harmonics)
! 1= <up|L.S|up>  ;  2= <up|L.S|dn>
  allocate(ls_cplx(2*ll+1,2*ll+1,2));ls_cplx=czero
  do im=1,2*ll+1
   mm=im-ll-1
   ls_cplx(im,im,1)=half*mm
   if ((mm+1)<= ll) then
    ls_cplx(im,im+1,2)=half*sqrt(real((ll-mm)*(ll+mm+1),kind=dp))
   end if
   if ((mm-1)>=-ll) then
    ls_cplx(im-1,im,2)=half*sqrt(real((ll+mm)*(ll-mm+1),kind=dp))
   end if
  end do

! Compute <sigma, S_lm1|L.S|S_lm2, sigma_prime> (S_lm=real spherical harmonics)
! 1= <up|L.S|up>  ;  2= <up|L.S|dn>
  lm0=ll**2
  do jm=1,2*ll+1
   jlm=lm0+jm;j0lm=jlm*(jlm-1)/2
   do im=1,jm
    ilm=lm0+im;klm=j0lm+ilm
    tmp(:)=czero
    do ii=1,2*ll+1
     do jj=1,2*ll+1
      tmp(:)=tmp(:)+ls_cplx(ii,jj,:)*CONJG(slm2ylm(ii,im))*slm2ylm(jj,jm)
     end do
    end do
    pawang%ls_ylm(1,klm,:)=REAL(tmp(:),kind=dp)
    pawang%ls_ylm(2,klm,:)=AIMAG(tmp(:))
   end do
  end do
  deallocate(ls_cplx,slm2ylm)

! End loop on l
 end do

 end subroutine pawlsylm

!!***
