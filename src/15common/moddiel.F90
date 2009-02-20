!{\src2tex{textfont=tt}}
!!****f* ABINIT/moddiel
!! NAME
!! moddiel
!!
!! FUNCTION
!! Precondition the residual, using a model dielectric function.
!! When cplex=1, assume q=(0 0 0), and vresid and vrespc will be REAL
!! When cplex=2, q must be taken into account, and vresid and vrespc will be COMPLEX
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, vhartr is REAL, if 2, vhartr is COMPLEX
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam.
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  optreal=1 if residual potential is is REAL space, 2 if it is in RECIPROCAL SPACE
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  vresid(cplex*nfft,nspden)=residual potential in REAL space       (if optreal==1)
!!                            residual potential in RECIPROCAL space (if optreal==2)
!!
!! OUTPUT
!!  vrespc(cplex*nfft,nspden)=preconditioned residual of the potential
!!                            in REAL space       if optreal==1
!!                            in RECIPROCAL space if optreal==2
!!
!! SIDE EFFECTS
!!
!! NOTES
!! optreal==2 is not compatible with cplex==1
!!
!! PARENTS
!!      newvtr3,prcref
!!
!! CHILDREN
!!      fourdp,leave_new,metric,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine moddiel(cplex,dielar,mpi_enreg,nfft,ngfft,nspden,optreal,paral_kgb,qphon,rprimd,vresid,vrespc)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments-------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,optreal,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: dielar(7),qphon(3),rprimd(3,3)
 real(dp),intent(in) :: vresid(cplex*nfft,nspden)
 real(dp),intent(out) :: vrespc(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,ifft,ig,ii,ii1,ing,ispden,me_fft,mg,n1,n2,n3,nproc_fft
 integer :: qeq0
 real(dp) :: ai,ar,dielng,diemac,diemac_inv,diemix,factor,gqg2p3,gqgm12,gqgm13
 real(dp) :: gqgm23,gs,gs2,gs3,l2g2,length2,ucvol
 character(len=500) :: message
!arrays
 integer :: id(3)
 real(dp) :: gmet(3,3),gprimd(3,3),potg0(4),rmet(3,3)
 real(dp),allocatable :: gq(:,:),work1(:,:),work2(:)

! *************************************************************************

!DEBUG
!write(6,*)' moddiel : enter '
!ENDDEBUG

!Check that cplex has an allowed value
 if(cplex/=1 .and. cplex/=2)then
  write(message, '(a,a,a,a,i3,a,a)' )ch10,&
&  ' moddiel : BUG -',ch10,&
&  '  From the calling routine, cplex=',cplex,ch10,&
&  '  but the only value allowed are 1 and 2.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 if(cplex==1.and.optreal==2)then
  write(message, '(a,a,a,a)' )ch10,&
&  ' moddiel : BUG -',ch10,&
&  '  When optreal=2, cplex must be 2.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!This is to allow q=0
 qeq0=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

!If cplex=1 then qphon should be 0 0 0
 if (cplex==1.and. qeq0/=1) then
  write(message, '(a,a,a,a,3e12.4,a,a)' ) ch10,&
&  ' moddiel : BUG -',ch10,&
&  '  cplex=1 but qphon=',qphon,ch10,&
&  '  qphon should be 0 0 0.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 dielng=dielar(2) ; diemac=dielar(3) ; diemix=dielar(4)

!DEBUG
!write(6,*)' moddiel : diemac, diemix =',diemac,diemix
!ENDDEBUG

 if(abs(diemac-1.0_dp)<1.0d-6)then

! Here, simple mixing is required, through macroscopic
! dielectric constant set to 1.0_dp .
  vrespc(:,:)=diemix*vresid(:,:)
! DEBUG
! write(0,*) "--------------------------------------------------->  moddiel: simple mixing ..."
! ENDDEBUG
 else

! Here, model dielectric function (G-diagonal operator)

  length2=(two_pi*dielng)**2
  diemac_inv=1.0_dp/diemac
  allocate(work1(2,nfft));if (optreal==1) allocate(work2(cplex*nfft))

! In order to speed the routine, precompute the components of g
  mg=maxval(ngfft)
  allocate(gq(3,mg))
  do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
    ig=ing-(ing/id(ii))*ngfft(ii)-1
    gq(ii,ing)=ig+qphon(ii)
   end do
  end do

! Do-loop on spins
! Note XG 010922 : I doubt the preconditioner is OK for the magnetisation
  do ispden=1,nspden

!  Do fft from real space (work2) to G space (work1)
   if (optreal==1) then
    work2(:)=vresid(:,ispden)
    call fourdp(cplex,work1,work2,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   else
    work1(:,:)=reshape(vresid(:,ispden),(/2,nfft/))
   end if

!  Store G=0 value
   potg0(ispden)=work1(re,1)

!  Triple loop, for the three dimensions
   do i3=1,n3
!   Precompute some products that do not depend on i2 and i1
    gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
    gqgm23=gq(3,i3)*gmet(2,3)*2
    gqgm13=gq(3,i3)*gmet(1,3)*2
    do i2=1,n2
     if (((i2-1)/(n2/nproc_fft))==me_fft) then
      gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
      gqgm12=gq(2,i2)*gmet(1,2)*2
      gqg2p3=gqgm13+gqgm12
      i23=n1*((i2-me_fft*n2/nproc_fft-1)+(n2/nproc_fft)*(i3-1))

!     Do the test that eliminates the Gamma point outside
!     of the inner loop
      ii1=1
      if(i2 == 1 .and. i3 == 1 .and. qeq0==1)then
!      if(i23==0 .and. qeq0==1)then: this changes with the number of fft procs...
!      and seems to be wrong.Pls check
       ii1=2
      end if

!     Here, unlike in hartre.f, the G=0 term is not eliminated, albeit
!     not changed.
      do i1=ii1,n1

!      One obtains the square of the norm of q+G (defined by i1,i2,i3)
       gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
       ifft=i1+i23

       l2g2=length2*gs
!      The model dielectric function is now computed
       factor = (l2g2+diemac_inv)/(l2g2+1.0_dp) * diemix
       work1(re,ifft)=work1(re,ifft)*factor
       work1(im,ifft)=work1(im,ifft)*factor

      end do
     end if
    end do
   end do

!  Might get rid of the G=0 term
!  if(qeq0==1)then
!  work1(re,1)=0.0_dp
!  work1(im,1)=0.0_dp
!  end if

!  Fourier transform
   if (optreal==1) then
    call fourdp(cplex,work1,work2,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
    vrespc(:,ispden)=work2(:)
   else
    vrespc(:,ispden)=reshape(work1(:,:),(/nfft*2/))
   end if

!  End of loop on spin polarizations
  end do

  deallocate(gq)
  deallocate(work1);if (optreal==1) deallocate(work2)

! End condition diemac/=1.0
 end if

!DEBUG
!write(6,*)' moddiel : exit '
!stop
!ENDDEBUG

end subroutine moddiel
!!***
