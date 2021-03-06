!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcden
!! NAME
!! xcden
!!
!! FUNCTION
!! Prepare density data before calling local or semi-local xc evaluators.
!!
!! NOTES
!! Can take into account a shift of the grid, for purpose of better accuracy.
!! Can also compute the gradient of the density, for use with GGAs.
!! Also eliminate eventual negative densities.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  ishift : if ==0, do not shift the xc grid (usual case); if ==1, shift the xc grid
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngrad : =1, only compute the density ; =2 also compute the
!!      gradient of the density. Note : ngrad**2 is also used to dimension rhonow
!!  nspden=number of spin-density components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor(cplex*nfft,nspden)=real space electron density in electrons/bohr**3, on the
!!   unshifted grid (total in first half and spin-up in second half if nspden=2)
!!
!! OUTPUT
!!  rhonow(cplex*nfft,nspden,ngrad*ngrad)=electron (spin)-density in real space and
!!     eventually its gradient, either on the unshifted grid (if ishift==0,
!!     then equal to rhor),or on the shifted grid
!!    rhonow(:,:,1)=electron density in electrons/bohr**3
!!    if ngrad==2 : rhonow(:,:,2:4)=gradient of electron density in electrons/bohr**4
!!
!! PARENTS
!!      mkvxcgga3,mkvxcstrgga3,rhohxc
!!
!! CHILDREN
!!      fourdp,leave_new,phase,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xcden (cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,paral_kgb,qphon,rhor,rhonow)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_13xc, except_this_one => xcden
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ishift,nfft,ngrad,nspden,paral_kgb
 type(MPI_type) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),qphon(3),rhor(cplex*nfft,nspden)
 real(dp),intent(out) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,id1,id2,id3,idir,ifft,ig1,ig2,ig3,inv1,inv2,inv3,inv_ifft
 integer :: ir,ispden,n1,n2,n3,qeq0
 real(dp) :: gc23_idir,gcart_idir,ph123i,ph123r,ph1i,ph1r,ph23i,ph23r,ph2i,ph2r
 real(dp) :: ph3i,ph3r,work_im,work_re
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: gcart1(:),gcart2(:),gcart3(:),ph1(:),ph2(:),ph3(:)
 real(dp),allocatable :: wkcmpx(:,:),work(:),workgr(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' xcden : enter '
!ENDDEBUG

 if (ishift/=0 .and. ishift/=1) then
  write(message, '(a,a,a,a,i12,a)' ) ch10,&
&  ' xcden: BUG -',ch10,&
&  '  ishift must be 0 or 1 ; input was',ishift,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 if (ngrad/=1 .and. ngrad/=2) then
  write(message, '(a,a,a,a,i12,a)' ) ch10,&
&  ' xcden: BUG -',ch10,&
&  '  ngrad must be 1 or 2 ; input was',ngrad,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!Keep local copy of fft dimensions
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Initialize computation of G in cartesian coordinates
 id1=n1/2+2  ; id2=n2/2+2  ; id3=n3/2+2

!Check whether q=0
 qeq0=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

 if(ishift==0)then

! Copy the input rhor in the new location. Will check on negative values later

  do ispden=1,nspden
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(cplex,nfft,rhonow,rhor)
   do ifft=1,cplex*nfft
    rhonow(ifft,ispden,1)=rhor(ifft,ispden)
   end do
!  $OMP END PARALLEL DO
  end do

 end if

 if(ishift==1 .or. ngrad==2)then

  allocate(wkcmpx(2,nfft),work(cplex*nfft))

  if(ishift==1)then
!  Precompute phases (The phases correspond to a shift of density on real space
!  grid from center at 0 0 0 to (1/2)*(1/n1,1/n2,1/n3).)
   allocate( ph1(2*n1),ph2(2*n2),ph3(2*n3) )
   call phase(n1,ph1)
   call phase(n2,ph2)
   call phase(n3,ph3)
  end if

  do ispden=1,nspden

!  Obtain rho(G) in wkcmpx from input rho(r)
!  Stop the xc timer
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(cplex,nfft,rhor,work)
   do ifft=1,cplex*nfft
    work(ifft)=rhor(ifft,ispden)
   end do
!  $OMP END PARALLEL DO

   call timab(82,1,tsec)
   call fourdp(cplex,wkcmpx,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   call timab(82,2,tsec)

!  If shift is required, multiply now rho(G) by phase, then generate rho(r+delta)
   if(ishift==1)then
!   $OMP PARALLEL DO PRIVATE(ifft,i1,i2,i3)&
!   $OMP&PRIVATE(ph1i,ph1r,ph123i,ph123r,ph2i,ph2r,ph23i,ph23r,ph3i,ph3r)&
!   $OMP&PRIVATE(work_im,work_re) &
!   $OMP&SHARED(n1,n2,n3,ph1,ph2,ph3,wkcmpx)
    do i3=1,n3
     ifft=(i3-1)*n1*n2
     ph3r=ph3(2*i3-1)
     ph3i=ph3(2*i3  )
     do i2=1,n2
      ph2r=ph2(2*i2-1)
      ph2i=ph2(2*i2  )
      ph23r=ph2r*ph3r-ph2i*ph3i
      ph23i=ph2i*ph3r+ph2r*ph3i
      if (((i2-1)/(n2/mpi_enreg%nproc_fft))==mpi_enreg%me_fft) then
       do i1=1,n1
        ifft=ifft+1
        ph1r=ph1(2*i1-1)
        ph1i=ph1(2*i1  )
        ph123r=ph1r*ph23r-ph1i*ph23i
        ph123i=ph1i*ph23r+ph1r*ph23i
!       Must use intermediate variables !
        work_re=ph123r*wkcmpx(1,ifft)-ph123i*wkcmpx(2,ifft)
        work_im=ph123i*wkcmpx(1,ifft)+ph123r*wkcmpx(2,ifft)
        wkcmpx(1,ifft)=work_re
        wkcmpx(2,ifft)=work_im
       end do
      end if
     end do
    end do
!   $OMP END PARALLEL DO
!   DEBUG
!   write(6,*)' xcden '
!   do i3=1,n3
!   inv3=n3+2-i3
!   if(i3==1)inv3=i3
!   do i2=1,n2
!   inv2=n2+2-i2
!   if(i2==1)inv2=i2
!   do i1=1,n1
!   ifft=i1+n1*((i2-1)+n2*(i3-1))
!   inv1=n1+2-i1
!   if(i1==1)inv1=i1
!   inv_ifft=inv1+n1*((inv2-1)+n2*(inv3-1))
!   write(6,'(i5,4es15.7)' )ifft,wkcmpx(1,ifft),wkcmpx(2,ifft),&
!   &                                  wkcmpx(1,inv_ifft),wkcmpx(2,inv_ifft)
!   end do
!   end do
!   end do
!   stop
!   ENDDEBUG
    call timab(82,1,tsec)
    call fourdp(cplex,wkcmpx,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
    call timab(82,2,tsec)
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(ispden,nfft,rhonow,work)
    do ifft=1,cplex*nfft
     rhonow(ifft,ispden,1)=work(ifft)
    end do
!   $OMP END PARALLEL DO

   end if

!  If gradient of the density is required, take care of the three components now
!  Note : this operation is applied on the eventually shifted rho(G)
   if(ngrad==2)then
    allocate(gcart1(n1),gcart2(n2),gcart3(n3),workgr(2,nfft))
    do idir=1,3

     do i1=1,n1
      ig1=i1-(i1/id1)*n1-1
      gcart1(i1)=gprimd(idir,1)*two_pi*(dble(ig1)+qphon(1))
     end do
!    Note that the G <-> -G symmetry must be maintained
     if(mod(n1,2)==0 .and. qeq0==1)gcart1(n1/2+1)=zero
     do i2=1,n2
      ig2=i2-(i2/id2)*n2-1
      gcart2(i2)=gprimd(idir,2)*two_pi*(dble(ig2)+qphon(2))
     end do
     if(mod(n2,2)==0 .and. qeq0==1)gcart2(n2/2+1)=zero
     do i3=1,n3
      ig3=i3-(i3/id3)*n3-1
      gcart3(i3)=gprimd(idir,3)*two_pi*(dble(ig3)+qphon(3))
     end do
     if(mod(n3,2)==0 .and. qeq0==1)gcart3(n3/2+1)=zero

!    $OMP PARALLEL DO PRIVATE(ifft,i1,i2,i3,gcart_idir,gc23_idir) &
!    $OMP&SHARED(gcart1,gcart2,gcart3,n1,n2,n3,wkcmpx,workgr)
     ifft = 0
     do i3=1,n3
      do i2=1,n2
       gc23_idir=gcart2(i2)+gcart3(i3)
       if (((i2-1)/(n2/mpi_enreg%nproc_fft))==mpi_enreg%me_fft) then
        do i1=1,n1
         ifft=ifft+1
         gcart_idir=gc23_idir+gcart1(i1)
!        Multiply by i 2pi G(idir)
         workgr(2,ifft)= gcart_idir*wkcmpx(1,ifft)
         workgr(1,ifft)=-gcart_idir*wkcmpx(2,ifft)
        end do
       end if
      end do
     end do
!    $OMP END PARALLEL DO
     call timab(82,1,tsec)
     call fourdp(cplex,workgr,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     call timab(82,2,tsec)
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(idir,ispden,nfft,rhonow,work)
     do ifft=1,cplex*nfft
      rhonow(ifft,ispden,1+idir)=work(ifft)
     end do
!    $OMP END PARALLEL DO

    end do
    deallocate(gcart1,gcart2,gcart3,workgr)
   end if

!  End loop on spins
  end do

  deallocate(wkcmpx,work)
  if(ishift==1)deallocate(ph1,ph2,ph3)

! End condition on ishift and ngrad
 end if

!DEBUG
!write(6,*)' i1, i2, i3, ifft, rhonow(ifft,1,1),rhonow(ifft,1,2:4)'
!ifft=0
!do i3=1,n3
!do i2=1,n2
!do i1=1,n1
!This is coded for cplex=1
!ifft=ifft+1
!write(6, '(4i5,4f16.6)' )i1,i2,i3,ifft,rhonow(ifft,1,1:4)
!end do
!end do
!end do
!stop
!ENDDEBUG

!DEBUG
!write(6,'(a)')&
!& '   ir              rhonow(:,:,1:4)'
!n1=ngfft(1)
!n2=ngfft(2)
!n3=ngfft(3)
!do ir=1,nfft
!i3=(ir-1)/n1/n2
!i2=(ir-1-i3*n1*n2)/n1
!i1=ir-1-i3*n1*n2-i2*n1
!This is coded for cplex=1
!write(message,'(i5,3i3,a,4es14.6)')ir,i1,i2,i3,' ',&
!&  rhonow(ir,1,1:4)
!call wrtout(06,message,'COLL')
!if(nspden==2)then
!write(message,'(a,2es14.6)')'               ',rhonow(ir,2,1:4)
!call wrtout(06,message,'COLL')
!end if
!end do
!ENDDEBUG

end subroutine xcden
!!***
