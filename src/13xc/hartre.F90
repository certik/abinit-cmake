!{\src2tex{textfont=tt}}
!!****f* ABINIT/hartre
!! NAME
!! hartre
!!
!! FUNCTION
!! Given rho(G), compute Hartree potential (=FFT of rho(G)/pi/(G+q)**2)
!! When cplex=1, assume q=(0 0 0), and vhartr will be REAL
!! When cplex=2, q must be taken into account, and vhartr will be COMPLEX
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! *Modified code to avoid if statements inside loops to skip G=0.
!!  Replaced if statement on G^2>gsqcut to skip G s outside where
!!  rho(G) should be 0.  Effect is negligible but gsqcut should be
!!  used to be strictly consistent with usage elsewhere in code.
!! *The speed-up is provided by doing a few precomputations outside
!!  the inner loop. One variable size array is needed for this (gq).
!!
!! INPUTS
!!  cplex= if 1, vhartr is REAL, if 2, vhartr is COMPLEX
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!         (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  izero=if 1, unbalanced components of Vhartree(g) are set to zero
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhog(2,nfft)=electron density in G space
!!
!! OUTPUT
!!  vhartr(cplex*nfft)=Hartree potential in real space, either REAL or COMPLEX
!!
!! PARENTS
!!      eneres3,loop3dte,rhohxc,scfcv3,tddft,prctfvw2
!!
!! CHILDREN
!!      fourdp,leave_new,timab,wrtout,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine hartre(cplex,gmet,gsqcut,izero,mpi_enreg,nfft,ngfft,paral_kgb,qphon,rhog,vhartr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,izero,nfft,paral_kgb
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),qphon(3),rhog(2,nfft)
 real(dp),intent(out) :: vhartr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,id2,id3,ig,ig2,ig3,ii,ii1,ing,n1,n2,n3,qeq0
 real(dp),parameter :: tolfix=1.000000001e0_dp
 real(dp) :: cutoff,den,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3
 character(len=500) :: message
!arrays
 integer :: id(3)
 real(dp) :: tsec(2)
 real(dp),allocatable :: gq(:,:),work1(:,:)

! *************************************************************************
!
!DEBUG
!write(6,*)' hartre : enter '
!write(6,*)' cplex,nfft,ngfft',cplex,nfft,ngfft
!write(6,*)' gsqcut=',gsqcut
!write(6,*)' qphon=',qphon
!write(6,*)' gmet=',gmet
!write(6,*)' maxval rhog=',maxval(abs(rhog))
!ENDDEBUG

!Keep track of total time spent in hartre
 call timab(10,1,tsec)

!Check that cplex has an allowed value
 if(cplex/=1 .and. cplex/=2)then
  write(message, '(a,a,a,a,i3,a,a)' )ch10,&
&  ' hartre : BUG -',ch10,&
&  '  From the calling routine, cplex=',cplex,ch10,&
&  '  but the only value allowed are 1 and 2.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Initialize a few quantities
 cutoff=gsqcut*tolfix
!This is to allow q=0
 qeq0=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

!If cplex=1 then qphon should be 0 0 0
 if (cplex==1.and. qeq0/=1) then
  write(message, '(a,a,a,a,3e12.4,a,a)' ) ch10,&
&  ' hartre: BUG -',ch10,&
&  '  cplex=1 but qphon=',qphon,ch10,&
&  '  qphon should be 0 0 0.'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 allocate(gq(3,max(n1,n2,n3)))
 do ii=1,3
  id(ii)=ngfft(ii)/2+2
  do ing=1,ngfft(ii)
   ig=ing-(ing/id(ii))*ngfft(ii)-1
   gq(ii,ing)=ig+qphon(ii)
  end do
 end do

 allocate(work1(2,nfft))
 id2=n2/2+2
 id3=n3/2+2
!Triple loop on each dimension
 do i3=1,n3
  ig3=i3-(i3/id3)*n3-1
! Precompute some products that do not depend on i2 and i1
  gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
  gqgm23=gq(3,i3)*gmet(2,3)*2
  gqgm13=gq(3,i3)*gmet(1,3)*2

  do i2=1,n2
   ig2=i2-(i2/id2)*n2-1
!  if (mpi_enreg%me_fft==modulo(i2,mpi_enreg%nproc_fft)) then
   if (((i2-1)/(n2/mpi_enreg%nproc_fft))==mpi_enreg%me_fft) then
    gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
    gqgm12=gq(2,i2)*gmet(1,2)*2
    gqg2p3=gqgm13+gqgm12

    i23=n1*(i2-mpi_enreg%me_fft*n2/mpi_enreg%nproc_fft-1+(n2/mpi_enreg%nproc_fft)*(i3-1))
!   Do the test that eliminates the Gamma point outside
!   of the inner loop
    ii1=1
    if(i23==0 .and. qeq0==1  .and. ig2==0 .and. ig3==0)then
     ii1=2
     work1(re,1+i23)=zero
     work1(im,1+i23)=zero
    end if

!   Final inner loop on the first dimension
!   (note the lower limit)
    do i1=ii1,n1
     gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
     ii=i1+i23
     if(gs<=cutoff)then
      den=piinv/gs
      work1(re,ii)=rhog(re,ii)*den
      work1(im,ii)=rhog(im,ii)*den
     else
      work1(re,ii)=zero
      work1(im,ii)=zero
     end if
!    End loop on i1
    end do
   end if
!  End loop on i2
  end do

! End loop on i3
 end do

 deallocate(gq)

!DEBUG
!write(6,*)' hartre : before fourdp'
!write(6,*)' cplex,nfft,ngfft',cplex,nfft,ngfft
!write(6,*)' maxval work1=',maxval(abs(work1))
!ENDDEBUG

!Set contribution of unbalanced components to zero
 if (izero==1) call zerosym(work1,2,mpi_enreg,n1,n2,n3)

!Fourier Transform Vhartree.
!Vh in reciprocal space was stored in work1
 call fourdp(cplex,work1,vhartr,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

!DEBUG
!write(6,*)' hartre : after fourdp'
!write(6,*)' cplex,nfft,ngfft',cplex,nfft,ngfft
!write(6,*)' maxval rhog=',maxval(abs(rhog))
!write(6,*)' maxval work1=',maxval(abs(work1))
!write(6,*)' maxval vhartr=',maxval(abs(vhartr))
!write(6,*)' hartre : set vhartr to zero'
!vhartr(:)=zero
!ENDDEBUG

 deallocate(work1)

 call timab(10,2,tsec)

end subroutine hartre
!!***
