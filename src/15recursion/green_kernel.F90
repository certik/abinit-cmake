!{\src2tex{textfont=tt}}
!!****f* ABINIT/green_kernel
!! NAME
!! green_kernel
!! 
!! FUNCTION
!! this routine compute the fourrier transform of the Green kernel for the 
!! recursion method  
!! 
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  inf_rmet=define the  infinitesimal metric : rprimd*(transpose(rprimd)) divided by the number of discretisation point
!!  inf_ucvol=volume of infinitesimal cell
!!  mult=variance of the Gaussian (=rtrotter/beta)
!!  mpi_enreg=informations about MPI parallelization
!!  ngfft=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nfft=total number of fft grid points
!!  prtvol=printing volume
!! 
!! OUTPUT
!!  ZT_p=fourier transforme of the Green kernel
!!  
!! PARENTS
!!      vtorhorec
!! 
!! CHILDREN
!!      fft3d, wrtout
!! 
!! NOTES 
!!  at this time :
!!       - need a rectangular box
!! 
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine green_kernel(ZT_p,inf_rmet,inf_ucvol,mult,mpi_enreg,ngfft,nfft,prtvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_12ffts
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nfft,prtvol
 real(dp),intent(in) :: inf_ucvol,mult
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: inf_rmet(3,3)
 real(dp),intent(out) :: ZT_p(1:2,0:nfft-1)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=7,n_green_max=5
 integer :: ii,isign,jj,kk,n_green,xx,yy,zz
 real(dp) :: acc,dsq
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: T_p(:)

! *************************************************************************

 dsq(ii,jj,kk)=inf_rmet(1,1)*(dble(ii))**2&
& +inf_rmet(2,2)*(dble(jj))**2&
& +inf_rmet(3,3)*(dble(kk))**2&
& +2.d0*(inf_rmet(1,2)*(dble(ii))*(dble(jj))&
& +inf_rmet(2,3)*(dble(jj))*(dble(kk))&
& +inf_rmet(3,1)*(dble(kk))*(dble(ii)))
 
 call timab(184,1,tsec)
 
 allocate(T_p(0:nfft-1))
 
!n_green should be better chosen for non rectangular cell
 do xx=1, n_green_max
  n_green = xx
  if(exp(-mult*dsq(xx*ngfft(1),0,0))<tol14 &
  .and. exp(-mult*dsq(0,xx*ngfft(2),0))<tol14 &
  .and. exp(-mult*dsq(0,0,xx*ngfft(3)))<tol14 )then
   exit
  end if
 end do
 
 acc = 0.d0
 T_p = 0.d0
 do ii = 0,ngfft(1)-1
  do jj = 0,ngfft(2)-1
   do kk = 0,ngfft(3)-1
    
    do xx=-n_green,n_green-1
     do yy=-n_green,n_green-1
      do zz=-n_green,n_green-1
       
       T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk) = T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)+ & 
&       exp(-mult*dsq(ii+xx*ngfft(1),jj+yy*ngfft(2),kk+zz*ngfft(3)))
       
      end do
     end do
    end do
    
    T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk) = (mult/pi)**(1.5d0)*T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)
    acc = acc + inf_ucvol* T_p(ii+ngfft(1)*jj+ngfft(1)*ngfft(2)*kk)
    
   end do
  end do
 end do
 
 T_p(:)= (1.d0/acc)*T_p(:)
!DEBUG
!if(prtvol==-level)then
!write(message,'(a,i)')'  n_green            ', n_green
!call wrtout(06,message,'COLL')
!write(message,'(a,3d12.3)')'  erreur_n_green     ', exp(-mult*dsq(n_green*ngfft(1),0,0)),&
!&        exp(-mult*dsq(0,n_green*ngfft(2),0)),&
!&        exp(-mult*dsq(0,0,n_green*ngfft(3)))
!call wrtout(06,message,'COLL')
!write(message,'(a,3d12.3)')'  erreur_troncat     ', T_p(ngfft(1)/2), &
!&       T_p(ngfft(1)*(ngfft(2)/2)),T_P(ngfft(1)*ngfft(2)*(ngfft(3)/2))
!call wrtout(06,message,'COLL')
!write(message,'(a,d)')'  erreurT_p          ', abs(acc-1.d0)
!call wrtout(06,message,'COLL')
!endif
!ENDDEBUG
 
 isign = -1
 call fourdp(1,ZT_p,T_p,isign,mpi_enreg,nfft,ngfft,1,0)
 
 deallocate(T_p)
 
 call timab(184,2,tsec)
 
end subroutine green_kernel
!!***
