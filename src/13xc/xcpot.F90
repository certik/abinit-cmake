!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcpot
!! NAME
!! xcpot
!!
!! FUNCTION
!! Process data after calling local or semi-local xc evaluators
!! The derivative of Exc with respect to the density is contained
!! in dnexcdn(:,:).
!! Eventually take also into account the derivative of Exc with respect to the
!! gradient of the density, contained in rhonow(:,:,2:4)
!! for use with GGAs
!! Can take into account a shift of the grid, for purpose of better accuracy
!! DESCRIPTION TO BE UPDATED : GGA is now available, and the use
!! of dnexcdn has changed
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
!!  dnexcdn(cplex*nfft,nspgrad)=derivative of Exc with respect to the (spin-)density,
!!    or to the norm of the gradient of the (spin-)density,
!!    further divided by the norm of the gradient of the (spin-)density
!!   The different components of dnexcdn will be
!!   for nspden=1,         dnexcdn(:,1)=d(n.exc)/d(n)
!!   for nspden>=2,        dnexcdn(:,1)=d(n.exc)/d(n_up)
!!                         dnexcdn(:,2)=d(n.exc)/d(n_down)
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  ishift : if ==0, do not shift the xc grid (usual case);
!!           if ==1, shift the xc grid
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngrad : =1, only take into account derivative wrt the density ;
!!          =2, also take into account derivative wrt the gradient of the density.
!!  nspden=number of spin-density components
!!  nspgrad=number of spin-density and spin-density-gradient components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhonow(cplex*nfft,nspden,ngrad*ngrad)=electron (spin)-density in real space and
!!     eventually its gradient already multiplied by the local partial derivative
!!     of the XC functional, either on the unshifted grid (if ishift==0,
!!     then equal to rhor), or on the shifted grid
!!    rhonow(:,:,1)=electron density in electrons/bohr**3
!!    if ngrad==2 : rhonow(:,:,2:4)=gradient of electron density in el./bohr**4,
!!     times local partial derivative of the functional, as required by the GGA
!!    In this routine, rhonow is used only in the GGA case (ngrad=2).
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  vxc(cplex*nfft,nspden)=xc potential (spin up in first half and spin down in
!!   second half if nspden>=2). Contribution from the present shifted
!!   or unshifted grid is ADDED to the input vxc data.
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

subroutine xcpot (cplex,dnexcdn,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,&
& nspgrad,paral_kgb,qphon,rhonow,vxc)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_13xc, except_this_one => xcpot
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ishift,nfft,ngrad,nspden,nspgrad,paral_kgb
 type(MPI_type) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: dnexcdn(cplex*nfft,nspgrad),gprimd(3,3),qphon(3)
 real(dp),intent(in) :: rhonow(cplex*nfft,nspden,ngrad*ngrad)
 real(dp),intent(out) :: vxc(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,id1,id2,id3,idir,ifft,ig1,ig2,ig3,ispden,n1,n2,n3,qeq0
 real(dp),parameter :: lowden=1.d-14,precis=1.d-15
 real(dp) :: gc23_idir,gcart_idir,ph123i,ph123r,ph1i,ph1r,ph23i,ph23r,ph2i,ph2r
 real(dp) :: ph3i,ph3r,work_im,work_re
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: gcart1(:),gcart2(:),gcart3(:),ph1(:),ph2(:),ph3(:)
 real(dp),allocatable :: wkcmpx(:,:),work(:),workgr(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' xcpot : enter '
!ENDDEBUG

 if (ishift/=0 .and. ishift/=1) then
  write(message, '(a,a,a,a,i12,a)' ) ch10,&
&  ' xcpot: BUG -',ch10,&
&  '  ishift must be 0 or 1 ; input was',ishift,'.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 if (ngrad/=1 .and. ngrad/=2 ) then
  write(message, '(a,a,a,a,i12,a)' ) ch10,&
&  ' xcpot: BUG -',ch10,&
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
! Add the value of dnexcdn to vxc
  do ispden=1,min(nspden,2)
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(cplex,dnexcdn,nfft,vxc)
   do ifft=1,cplex*nfft
    vxc(ifft,ispden)=vxc(ifft,ispden)+dnexcdn(ifft,ispden)
   end do
!  $OMP END PARALLEL DO
  end do

 end if

!If the grid is shifted, or if gradient corrections are present,
!there must be FFTs.
 if(ishift==1 .or. ngrad==2)then

  allocate( wkcmpx(2,nfft),work(cplex*nfft) )

  if(ishift==1)then
   allocate( ph1(2*n1),ph2(2*n2),ph3(2*n3) )
!  Precompute phases (The phases correspond to a shift of density on real space
!  grid from center at 0 0 0 to (1/2)*(1/n1,1/n2,1/n3).)
   call phase(n1,ph1)
   call phase(n2,ph2)
   call phase(n3,ph3)
  end if

  do ispden=1,min(nspden,2)

!  Initialize wkcmpx either to 0 or to the shifted vxc value
   if(ishift==0)then
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(nfft,wkcmpx)
    do ifft=1,nfft
     wkcmpx(:,ifft)=zero
    end do
!   $OMP END PARALLEL DO

   else
!   Obtain dnexcdn(G)*phase in wkcmpx from input dnexcdn(r+delta)
!   Stop the xc timer
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(cplex,dnexcdn,ispden,nfft,work)
    do ifft=1,cplex*nfft
     work(ifft)=dnexcdn(ifft,ispden)
    end do
!   $OMP END PARALLEL DO
    call timab(82,1,tsec)
    call fourdp(cplex,wkcmpx,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
    call timab(82,2,tsec)
   end if

!  If gradient correction is present, take care of the three components now
!  Note : this operation is done on the eventually shifted grid
   if(ngrad==2)then
    allocate(gcart1(n1),gcart2(n2),gcart3(n3),workgr(2,nfft))
    do idir=1,3

!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(cplex,idir,ispden,nfft,rhonow,work)
     do ifft=1,cplex*nfft
      work(ifft)=rhonow(ifft,ispden,1+idir)
     end do
!    $OMP END PARALLEL DO

     call timab(82,1,tsec)
     call fourdp(cplex,workgr,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     call timab(82,2,tsec)
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
!    $OMP PARALLEL DO PRIVATE(ifft,i1,i2,i3,gc23_idir,gcart_idir) &
!    $OMP&SHARED(gcart1,gcart2,gcart3,n1,n2,n3,wkcmpx,workgr)
     ifft = 0
     do i3=1,n3
      do i2=1,n2
       gc23_idir=gcart2(i2)+gcart3(i3)
       if (((i2-1)/(n2/mpi_enreg%nproc_fft))==mpi_enreg%me_fft) then
        do i1=1,n1
         ifft=ifft+1
         gcart_idir=gc23_idir+gcart1(i1)
!        Multiply by - i 2pi G(idir) and accumulate in wkcmpx
         wkcmpx(1,ifft)=wkcmpx(1,ifft)+gcart_idir*workgr(2,ifft)
         wkcmpx(2,ifft)=wkcmpx(2,ifft)-gcart_idir*workgr(1,ifft)
        end do
       end if
      end do
     end do
!    $OMP END PARALLEL DO

    end do
    deallocate(gcart1,gcart2,gcart3,workgr)
   end if

!  wkcmpx(:,:) contains now the full exchange-correlation potential, but
!  eventually for the shifted grid

   if(ishift==1)then
!   Take away the phase to get dnexcdn(G)
    ifft=0
!   $OMP PARALLEL DO PRIVATE(ifft,i1,i2,i3)&
!   $OMP&PRIVATE(ph1i,ph1r,ph123i,ph123r,ph2i,ph2r,ph23i,ph23r,ph3i,ph3r)&
!   $OMP&PRIVATE(work_im,work_re) &
!   $OMP&SHARED(n1,n2,n3,ph1,ph2,ph3,wkcmpx)
    do i3=1,n3
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
!       Multiply by phase.  Must use intermediate variables !
        work_re= ph123r*wkcmpx(1,ifft)+ph123i*wkcmpx(2,ifft)
        work_im=-ph123i*wkcmpx(1,ifft)+ph123r*wkcmpx(2,ifft)
        wkcmpx(1,ifft)=work_re
        wkcmpx(2,ifft)=work_im
       end do
      end if
     end do
    end do
!   $OMP END PARALLEL DO
   end if

   call timab(82,1,tsec)
   call fourdp(cplex,wkcmpx,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   call timab(82,2,tsec)

!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(cplex,ispden,nfft,vxc,work)
   do ifft=1,cplex*nfft
    vxc(ifft,ispden)=vxc(ifft,ispden)+work(ifft)
   end do
!  $OMP END PARALLEL DO

!  End loop on spins
  end do

  if(ishift==1) deallocate(ph1,ph2,ph3)
  deallocate(wkcmpx,work)

! End condition on ishift
 end if

end subroutine xcpot
!!***
