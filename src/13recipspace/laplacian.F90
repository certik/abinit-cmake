!{\src2tex{textfont=tt}}
!!****f* ABINIT/laplacian
!! NAME
!! laplacian
!!
!! FUNCTION
!! compute the laplacian of a function defined in real space
!! the code is written in the way of /3xc/xcden.F90
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! nfft=number of points of the fft grid
!! nfunc=number of functions on the grid for which the laplacian is to be calculated
!! ngfftf(18)
!! (optional) rdfuncr(nfft,nfunc)=real(dp) discretized functions in real space
!! gprimd(3,3)
!! (optional) g2cart_in(nfft) = G**2 on the grid
!!
!! OUTPUT
!! (optional) laplacerdfuncr = laplacian in real space of the functions in rdfuncr
!! (optional) rdfuncg = real(dp) discretized functions in fourier space
!! (optional) laplacerdfuncg = real(dp) discretized laplacian of the functions in fourier space
!! (optional) g2cart_out(nfft) = G**2 on the grid
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!!
!! NOTES
!!
!! PARENTS
!!      frskerker1,frskerker2,ftfvw1,ftfvw2,moddiel_csrb,prcrskerker1
!!      prcrskerker2,prctfvw1,prctfw3
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine laplacian(gprimd,mpi_enreg,nfft,nfunc,ngfft,paral_kgb,rdfuncr,&
     & laplacerdfuncr,rdfuncg_out,laplacerdfuncg_out,g2cart_out,rdfuncg_in,g2cart_in)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
  !integer,parameter::dp=8,two_pi=3
!scalars
 integer,intent(in) :: nfft,nfunc,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout),optional :: laplacerdfuncr(nfft,nfunc)
 real(dp),intent(inout),optional,target :: rdfuncr(nfft,nfunc)
 real(dp),intent(out),optional,target :: g2cart_in(nfft),g2cart_out(nfft)
 real(dp),intent(out),optional,target :: laplacerdfuncg_out(2,nfft,nfunc)
 real(dp),intent(out),optional,target :: rdfuncg_in(2,nfft,nfunc)
 real(dp),intent(out),optional,target :: rdfuncg_out(2,nfft,nfunc)

!Local variables-------------------------------
!scalars
 integer :: count,i1,i2,i3,id1,id2,id3,ifft,ifunc,ig1,ig2,ig3,ii1,ispden,n1,n2
 integer :: n3
 real(dp) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
!arrays
 real(dp),pointer :: g2cart(:),laplacerdfuncg(:,:,:),rdfuncg(:,:,:)

! *************************************************************************

!Keep local copy of fft dimensions
!write(0,*) 'debug laplacian 0'
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)


 if(present(laplacerdfuncg_out)) then
  laplacerdfuncg => laplacerdfuncg_out
 else
  allocate(laplacerdfuncg(2,nfft,nfunc))
 end if

!change the real density rdfuncr on real space on the real density
!rdfuncg in reciprocal space
 if(.not.present(rdfuncg_in)) then
  if(present(rdfuncg_out)) then
   rdfuncg => rdfuncg_out
  else
   allocate(rdfuncg(2,nfft,nfunc))
  end if
  if(present(rdfuncr)) then
   do ifunc=1,nfunc
    call fourdp(1,rdfuncg(:,:,ifunc),rdfuncr(:,ifunc),-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   end do
  end if
 else
  rdfuncg => rdfuncg_in
 end if



!apply the laplacian on laplacerdfuncr
!code from /3xc/xcden.F90
!see also src/5common/hatre.F90 and src/5common/moddiel.F90
!Keep local copy of fft dimensions
!Initialize computation of G^2 in cartesian coordinates
 if(.not.present(g2cart_in)) then
  if(present(g2cart_out)) then
   g2cart => g2cart_out
  else
   allocate(g2cart(nfft))
  end if
  id1=int(n1/2)+2
  id2=int(n2/2)+2
  id3=int(n3/2)+2
  ifft=0
  count=0
  do i3=1,n3
   ifft=(i3-1)*n1*n2
   ig3=i3-int(i3/id3)*n3-1
   do i2=1,n2
    ig2=i2-int(i2/id2)*n2-1
    ii1=1
    do i1=ii1,n1
     ig1=i1-int(i1/id1)*n1-1
     ifft=ifft+1

     b11=gprimd(1,1)*real(ig1,dp)
     b21=gprimd(2,1)*real(ig1,dp)
     b31=gprimd(3,1)*real(ig1,dp)
     b12=gprimd(1,2)*real(ig2,dp)
     b22=gprimd(2,2)*real(ig2,dp)
     b32=gprimd(3,2)*real(ig2,dp)
     b13=gprimd(1,3)*real(ig3,dp)
     b23=gprimd(2,3)*real(ig3,dp)
     b33=gprimd(3,3)*real(ig3,dp)

     g2cart(ifft)=( &
&     (b11+b12+b13)**2&
&     +(b21+b22+b23)**2&
&     +(b31+b32+b33)**2&
&     )
     do ifunc=1,nfunc
!     compute the laplacien in fourrier space
!     that is * (i x 2pi x G)**2
      laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
      laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
     end do
    end do
   end do
  end do
  if(.not.present(g2cart_out)) deallocate(g2cart)
 else
  g2cart => g2cart_in
  do ifunc=1,nfunc
   do ifft=1,nfft
!   compute the laplacien in fourrier space
!   that is * (i x 2pi x G)**2
    laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
    laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
   end do
  end do
 end if

!get the result back into real space
 if(present(laplacerdfuncr)) then
  do ifunc=1,nfunc
   call fourdp(1,laplacerdfuncg(:,:,ifunc),laplacerdfuncr(:,ifunc),1,mpi_enreg,nfft,ngfft,paral_kgb,0)
  end do
 end if

!deallocate pointers
 if((.not.present(rdfuncg_in)).and.(.not.present(rdfuncg_in))) deallocate(rdfuncg)
 if(.not.present(laplacerdfuncg_out)) deallocate(laplacerdfuncg)
end subroutine laplacian

!end of file
!!***
