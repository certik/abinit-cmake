!{\src2tex{textfont=tt}}
!!****f* ABINIT/fourdp
!!
!! NAME
!! fourdp
!!
!! FUNCTION
!! Conduct Fourier transform of REAL or COMPLEX function f(r)=fofr defined on
!! fft grid in real space, to create complex f(G)=fofg
!! defined on full fft grid in reciprocal space, in full storage mode,
!! or the reverse operation.
!! For the reverse operation, the final data is divided by nfftot.
!! REAL case when cplex=1, COMPLEX case when cplex=2
!! Usually used for density and potentials.
!!
!! There are two different possibilities :
!!  fftalgb=0 means using the complex-to-complex FFT routine,
!!   irrespective of the value of cplex
!!  fftalgb=1 means using a real-to-complex FFT or a complex-to-complex FFT,
!!   depending on the value of cplex.
!!  The only real-to-complex FFT available is from SGoedecker library.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex=1 if fofr is real, 2 if fofr is complex
!! isign=sign of Fourier transform exponent: current convention uses
!!  +1 for transforming from G to r, -1 for transforming from r to G.
!! mpi_enreg=information about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! tim_fourdp=timing code of the calling routine (can be set to 0 if not attributed)
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/Output
!! fofg(2,nfft)=f(G), complex.
!! fofr(cplex*nfft)=input function f(r) (real or complex)
!!
!! TODO
!! work2 should be a pointer, and allocated inside ccfft,
!! if really needed ...
!!
!! give definition for paral_kgb
!!
!! PARENTS
!!      atm2fft,atm2fft3,calc_coh,cppm2par,cppm3par,cppm4par,dens_in_sph
!!      dieltcel,difvxc,dyfro3,energy,fft_onewfn,fftwfn,filterpot,forces
!!      fresidrsp,ftfvw1,ftfvw2,green_kernel,gstate,hartre,hartre1,hartrestr
!!      initro,jellium,jvec_to_B,kxc_alda,kxc_eok,ladielmt,laplacian,loop3dte
!!      loper3,make_efg_el,make_fc_el,mklocl_realspace,mklocl_recipspace
!!      moddiel,moddiel_csrb,newrho,newvtr,nonlinear,nres2vres,odamix,pawmknhat
!!      pawmknhat3,prcref,prcref_PMA,prctfvw1,prctfvw2,prctfw3,rdm,recursion
!!      redgr,respfn,rho_tw_g,scfcv,scfcv3,screening,sigma,stress,symrhg,tddft
!!      transgrid,vloca3,vlocalstr,vtorho,vtorho3,xc_kernel,xcden,xcpot
!!
!! CHILDREN
!!      back,ccfft,forw,leave_new,sg_ctrig,sg_fftx,sg_ffty,sg_fftz,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fourdp(cplex,fofg,fofr,isign,mpi_enreg,nfft,ngfft,paral_kgb,tim_fourdp)

 use defs_basis
 use defs_datatypes
 use defs_fftdata


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts, except_this_one => fourdp
 use interfaces_lib01fftnew
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,isign,nfft,paral_kgb,tim_fourdp
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(inout) :: fofg(2,nfft),fofr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: mfac=11
 integer :: fftalg,fftalga,fftalgb,fftcache,i1,i2,i3,ic1,ic2,ic3,index
 integer :: inplace,ir,max_index,me_fft,n1,n1half1,n1halfm,n2,n2half1,n3,n4
 integer :: n4half1,n5,n5half1,n6,nd2proc,nd3proc,ndat,normalized,nproc_fft
 real(dp) :: inter,ris,xnorm
 character(len=500) :: message
!arrays
 integer :: aft1(mfac),aft2(mfac),aft3(mfac),bef1(mfac),bef2(mfac),bef3(mfac)
 integer :: ind1(mg),ind2(mg),ind3(mg),now1(mfac),now2(mfac),now3(mfac)
 real(dp) :: trig1(2,mg),trig2(2,mg),trig3(3,mg),tsec(2)
 real(dp),allocatable :: wk2d_a(:,:,:,:),wk2d_b(:,:,:,:),wk2d_c(:,:,:,:)
 real(dp),allocatable :: wk2d_d(:,:,:,:),work1(:,:,:,:),work2(:,:,:,:)
 real(dp),allocatable :: workf(:,:,:,:),workr(:,:,:,:)

! *************************************************************************

!DEBUG
! write(6,*)' fourdp : enter, fftalg,isign=',ngfft(7),isign
!ENDDEBUG
!Keep track of timing
 call timab(260+tim_fourdp,1,tsec)
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
 me_fft=ngfft(11) ; nproc_fft=ngfft(10)
!DEBUG
!write(6,*)' fourdp :me_fft',me_fft,'nproc_fft',nproc_fft,'nfft',nfft
!ENDDEBUG
 fftcache=ngfft(8)
 fftalg=ngfft(7)
 fftalga=fftalg/100
 fftalgb=mod(fftalg,100)/10
 ris=dble(isign)
 xnorm=1.0d0/dble(n1*n2*n3)

 if(fftalgb/=0 .and. fftalgb/=1) then
  write(message, '(a,a,a,a,i4,a,a,a,a,a)' )ch10,&
&  ' fourdp : ERROR -',ch10,&
&  '  The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&  '  The second digit (fftalg(B)) must be 0 or 1.',ch10,&
&  '  Action : change fftalg in your input file.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

 if(fftalgb==1 .and. (fftalga/=1 .and. fftalga/=4))then
  write(message, '(a,a,a,a,i4,a,a,a,a,a)' )ch10,&
&  ' fourdp : ERROR -',ch10,&
&  '  The input algorithm number fftalg=',fftalg,' is not allowed.',ch10,&
&  '  When fftalg(B) is 1, the allowed values for fftalg(A) are 1 and 4.',ch10,&
&  '  Action : change fftalg in your input file.'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

 if (n4<n1.or.n5<n2.or.n6<n3) then
  write(message, '(a,a,a,a,3i8,a,a,3i8,a)' ) ch10,&
&   ' fourdp: BUG -',ch10,&
&   '  Each of n4,n5,n6=',n4,n5,n6,ch10,&
&   '  must be >=      n1, n2, n3 =',n1,n2,n3,'.'
  call wrtout(06,message,'PERS')
  call leave_new('PERS')
 end if

!---------------------------------------------------------
!Here, deal  with the new SG FFT, complex-to-complex case
 if( fftalga==4 .and. (fftalgb==0 .or. cplex==2) )then

   nd2proc=((n2-1)/nproc_fft) +1
   nd3proc=((n6-1)/nproc_fft) +1
  allocate(workr(2,n4,n5,nd3proc),workf(2,n4,n6,nd2proc))
  max_index=0
  if(isign==1)then
   do i3=1,n3
    do i2=1,n2
     if (((i2-1)/(n2/nproc_fft))==me_fft) then
      index=n1*(i2-me_fft*n2/nproc_fft-1+(n2/nproc_fft)*(i3-1))
     if (index > max_index) max_index=index
     do i1=1,n1

      workf(1,i1,i3,i2-me_fft*n2/nproc_fft)=fofg(1,i1+index)
      workf(2,i1,i3,i2-me_fft*n2/nproc_fft)=fofg(2,i1+index)
     end do
     end if
    end do
   end do
   ndat=1 ;
   call back(2,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,2,paral_kgb,workf,workr)

   if(cplex==1)then
    do i3=1,n3
     if (((i3-1)/(n3/nproc_fft))==me_fft) then
     do i2=1,n2
      index=n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1))
      do i1=1,n1
       fofr(i1+index)=workr(1,i1,i2,i3-me_fft*n3/nproc_fft)
      end do
     end do
     end if
    end do
   else if(cplex==2)then
    do i3=1,n3
     if (((i3-1)/(n3/nproc_fft))==me_fft) then
     do i2=1,n2
      index=2*n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1))
      do i1=1,n1
       fofr(2*i1-1+index)=workr(1,i1,i2,i3-me_fft*n3/nproc_fft)
       fofr(2*i1  +index)=workr(2,i1,i2,i3-me_fft*n3/nproc_fft)
      end do
     end do
     end if
    end do
   end if

  else if(isign==-1)then

   if(cplex==1)then
    do i3=1,n3
     if (((i3-1)/(n3/nproc_fft))==me_fft) then
     do i2=1,n2
      index=n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1))
!      index=n1*(i2-1+n2*(i3-1))
      do i1=1,n1
       workr(1,i1,i2,i3-me_fft*n3/nproc_fft)=fofr(i1+index)
       workr(2,i1,i2,i3-me_fft*n3/nproc_fft)=zero
      end do
     end do
    end if
    end do
   else if(cplex==2)then
    do i3=1,n3
     if (((i3-1)/(n3/nproc_fft))==me_fft) then
     do i2=1,n2
       index=2*n1*(i2-1+n2*(i3-me_fft*n3/nproc_fft-1))
      do i1=1,n1
       workr(1,i1,i2,i3-me_fft*n3/nproc_fft)=fofr(2*i1-1+index)
       workr(2,i1,i2,i3-me_fft*n3/nproc_fft)=fofr(2*i1  +index)
      end do
     end do
     end if
    end do
   end if

   ndat=1
   call forw(2,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4,nd2proc,nd3proc,2,paral_kgb,workr,workf)

!  Transfer fft output to the original fft box
   do i2=1,n2
     if (((i2-1)/(n2/nproc_fft))==me_fft) then
    do i3=1,n3
     index=n1*(i2-me_fft*n2/nproc_fft-1+(n2/nproc_fft)*(i3-1))
     do i1=1,n1
      fofg(1,i1+index)=workf(1,i1,i3,i2-me_fft*n2/nproc_fft)*xnorm
      fofg(2,i1+index)=workf(2,i1,i3,i2-me_fft*n2/nproc_fft)*xnorm
     end do
    end do
    end if
   end do
 do ir=1,40
 end do

  end if ! isign

  deallocate(workr,workf)

 end if

!---------------------------------------------------------
!Here, deal with the new SG FFT, with real-to-complex
 if(fftalga==4 .and. fftalgb==1 .and. cplex==1)then
!DEBUG
! write(6,*)' stop here '
! stop
!ENDDEBUG

  n1half1=n1/2+1 ; n1halfm=(n1+1)/2
  n2half1=n2/2+1
! n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
  n4half1=(n1half1/2)*2+1
  n5half1=(n2half1/2)*2+1
  allocate(workr(2,n4half1,n5,n6),workf(2,n4,n6,n5half1))
  if(isign==1)then
   do i3=1,n3
    do i2=1,n2half1
     index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1
      workf(1,i1,i3,i2)=fofg(1,i1+index)
      workf(2,i1,i3,i2)=fofg(2,i1+index)
     end do
    end do
   end do

!DEBUG
!  write(6,*)' fourdp : before back '
!  do i2=1,n2half1
!   do i3=1,n3
!    do i1=1,n1
!     index=i1+n1*(i2-1+n2*(i3-1))
!     write(6, '(3i4,2es16.6)' ) i1,i3,i2,workf(1:2,i1,i3,i2)
!    end do
!   end do
!  end do
!ENDDEBUG

   ndat=1;
   nd2proc=((n5-1)/nproc_fft) +1
   nd3proc=((n6-1)/nproc_fft) +1
! modifier l appel ? n5half1 et n6 ?
   call back(cplex,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4half1,n5half1,n6,2,paral_kgb,workf,workr)

!DEBUG
!  write(6,*)' fourdp : after back '
!  do i3=1,n3
!   do i2=1,n2
!    do i1=1,n1half1
!     index=i1+n1*(i2-1+n2*(i3-1))
!     write(6, '(3i4,2es16.6)' ) i1,i2,i3,workr(1:2,i1,i2,i3)
!    end do
!   end do
!  end do
!  stop
!ENDDEBUG

   do i3=1,n3
    do i2=1,n2
     index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1half1-1
!    copy data
      fofr(2*i1-1+index)=workr(1,i1,i2,i3)
      fofr(2*i1  +index)=workr(2,i1,i2,i3)
     end do
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      fofr(n1+index)=workr(1,n1half1,i2,i3)
     end if
    end do
   end do

!DEBUG
!  write(6,*)' fourdp : stop, isign=1'
!  do i3=1,n3
!   do i2=1,n2
!    do i1=1,n1
!     index=i1+n1*(i2-1+n2*(i3-1))
!     write(6, '(3i4,2es16.6)' ) i1,i2,i3,fofr(index)
!    end do
!   end do
!  end do
!  stop
!ENDDEBUG

  else if(isign==-1)then
   do i3=1,n3
    do i2=1,n2
     index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1half1-1
      workr(1,i1,i2,i3)=fofr(2*i1-1+index)
      workr(2,i1,i2,i3)=fofr(2*i1  +index)
     end do
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      workr(1,n1half1,i2,i3)=fofr(n1+index)
      workr(2,n1half1,i2,i3)=0.0d0
     end if
    end do
   end do
   ndat=1
   call forw(cplex,mpi_enreg,ndat,n1,n2,n3,n4,n5,n6,n4half1,n5half1,n6,2,paral_kgb,workr,workf)
!  Transfer fft output to the original fft box
   do i3=1,n3
    do i2=1,n2half1
     index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1
      fofg(1,i1+index)=workf(1,i1,i3,i2)*xnorm
      fofg(2,i1+index)=workf(2,i1,i3,i2)*xnorm
     end do
    end do
!   Complete missing values with complex conjugate
!   Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
    if(n2half1>2)then
     do i2=2,n2+1-n2half1
      index=n1*((n2+2-i2)-1)
      if(i3/=1)index=index+n1*n2*((n3+2-i3)-1)
      fofg(1,1+index)= workf(1,1,i3,i2)*xnorm
      fofg(2,1+index)=-workf(2,1,i3,i2)*xnorm
      do i1=2,n1
       fofg(1,n1+2-i1+index)= workf(1,i1,i3,i2)*xnorm
       fofg(2,n1+2-i1+index)=-workf(2,i1,i3,i2)*xnorm
      end do
     end do
    end if
  end do
!DEBUG
!  write(6,*)' fourdp : stop, isign=-1'
!  do i3=1,n3
!   do i2=1,n2
!    do i1=1,n1
!     index=i1+n1*(i2-1+n2*(i3-1))
!     write(6, '(3i4,2es16.6)' ) i1,i2,i3,fofg(1:2,index)
!    end do
!   end do
!  end do
!  stop
!ENDDEBUG

  end if ! isign
  deallocate(workr,workf)

 end if

!---------------------------------------------------------
!Here, one calls the complex-to-complex FFT subroutine
 if( (fftalgb==0 .or. cplex==2) .and. fftalga/=4 )then

  allocate(work1(2,n4,n5,n6),work2(2,n4,n5,n6))

  if(isign==1)then

!  Transfer fofg to the expanded fft box
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofg,n1,n2,n3,work1)
   do i3=1,n3
    do i2=1,n2
     index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1
      work1(1,i1,i2,i3)=fofg(1,i1+index)
      work1(2,i1,i2,i3)=fofg(2,i1+index)
     end do
    end do
   end do
!$OMP END PARALLEL DO

!  Actual 3D FFT
   call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&   n1,n2,n3,n4,n5,n6,1,ngfft,2,paral_kgb,work1,work2)

!  Take data from expanded box and put it in the original box.
   if(cplex==1)then
!   REAL case

    if(inplace==0)then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofr,n1,n2,n3,work2)
     do i3=1,n3
      do i2=1,n2
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
        fofr(i1+index)=work2(1,i1,i2,i3)
       end do
      end do
     end do
!$OMP END PARALLEL DO
    else if(inplace==1)then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofr,n1,n2,n3,work1)
     do i3=1,n3
      do i2=1,n2
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
        fofr(i1+index)=work1(1,i1,i2,i3)
       end do
      end do
     end do
!$OMP END PARALLEL DO
    end if

   else
!   COMPLEX case

    if(inplace==0)then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofr,n1,n2,n3,work2)
     do i3=1,n3
      do i2=1,n2
       index=2*n1*(i2-1+n2*(i3-1))
       do i1=1,n1
        fofr(2*i1-1+index)=work2(1,i1,i2,i3)
        fofr(2*i1  +index)=work2(2,i1,i2,i3)
       end do
      end do
     end do
!$OMP END PARALLEL DO
    else if(inplace==1)then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofr,n1,n2,n3,work1)
     do i3=1,n3
      do i2=1,n2
       index=2*n1*(i2-1+n2*(i3-1))
       do i1=1,n1
        fofr(2*i1-1+index)=work1(1,i1,i2,i3)
        fofr(2*i1  +index)=work1(2,i1,i2,i3)
       end do
      end do
     end do
!$OMP END PARALLEL DO
    end if

   end if

  else if(isign==-1)then

!  Insert fofr into the augmented fft box
   if(cplex==1)then
!   REAL case
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofr,n1,n2,n3,work1)
    do i3=1,n3
     do i2=1,n2
      index=n1*(i2-1+n2*(i3-1))
      do i1=1,n1
!      copy data
       work1(1,i1,i2,i3)=fofr(i1+index)
       work1(2,i1,i2,i3)=0.0d0
      end do
     end do
    end do
!$OMP END PARALLEL DO
   else
!   COMPLEX case
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofr,n1,n2,n3,work1)
    do i3=1,n3
     do i2=1,n2
      index=2*n1*(i2-1+n2*(i3-1))
      do i1=1,n1
!      copy data
       work1(1,i1,i2,i3)=fofr(2*i1-1+index)
       work1(2,i1,i2,i3)=fofr(2*i1  +index)
      end do
     end do
    end do
!$OMP END PARALLEL DO
   end if ! cplex

!  Actual 3D FFT
   call ccfft(fftalga,fftcache,inplace,isign,mpi_enreg,normalized,&
&   n1,n2,n3,n4,n5,n6,1,ngfft,2,paral_kgb,work1,work2)

!  Transfer fft output to the original fft box
   if(normalized==0)then
    if(inplace==0)then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofg,n1,n2,n3,work2)
     do i3=1,n3
      do i2=1,n2
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
        fofg(1,i1+index)=work2(1,i1,i2,i3)*xnorm
        fofg(2,i1+index)=work2(2,i1,i2,i3)*xnorm
       end do
      end do
     end do
!$OMP END PARALLEL DO
    else if(inplace==1)then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofg,n1,n2,n3,work1)
     do i3=1,n3
      do i2=1,n2
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
        fofg(1,i1+index)=work1(1,i1,i2,i3)*xnorm
        fofg(2,i1+index)=work1(2,i1,i2,i3)*xnorm
       end do
      end do
     end do
!$OMP END PARALLEL DO
    end if
   else if(normalized==1)then
    if(inplace==0)then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofg,n1,n2,n3,work2)
     do i3=1,n3
      do i2=1,n2
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
        fofg(1,i1+index)=work2(1,i1,i2,i3)
        fofg(2,i1+index)=work2(2,i1,i2,i3)
       end do
      end do
     end do
!$OMP END PARALLEL DO
    else if(inplace==1)then
!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) &
!$OMP&SHARED(fofg,n1,n2,n3,work1)
     do i3=1,n3
      do i2=1,n2
       index=n1*(i2-1+n2*(i3-1))
       do i1=1,n1
        fofg(1,i1+index)=work1(1,i1,i2,i3)
        fofg(2,i1+index)=work1(2,i1,i2,i3)
       end do
      end do
     end do
!$OMP END PARALLEL DO
    end if
   end if ! normalized==0 or 1

! Endif choice of isign
  end if

  deallocate(work1,work2)

!End simple algorithm
 end if

!---------------------------------------------------------
!Here sophisticated algorithm based on S. Goedecker routines,
!only for the REAL case.
!Take advantage of the fact that fofr is real, and that fofg
!has corresponding symmetry properties.

 if( (fftalgb==1 .and. cplex==1) .and. fftalga/=4 )then

! Check that dimension is not exceeded
  if (n1>mg .or. n2>mg .or. n3>mg) then
   write(message, '(a,a,a,a,3i10,a,a,a,i10,a)' ) ch10,&
&   ' fourdp : BUG -',ch10,&
&   '  One of the dimensions n1,n2,n3=',n1,n2,n3,',',ch10,&
&   '  exceeds allowed dimension mg=',mg,'.'
   call wrtout(06,message,'PERS')
   call leave_new('PERS')
  end if

  n1half1=n1/2+1 ; n1halfm=(n1+1)/2
  n2half1=n2/2+1
! n4half1 or n5half1 are the odd integers >= n1half1 or n2half1
  n4half1=(n1half1/2)*2+1
  n5half1=(n2half1/2)*2+1

! The size of these arrays should be reduced.
!  allocate(wk2d_a(2,n4,n5,1),wk2d_b(2,n4,n5,1))
!  allocate(wk2d_c(2,2*n1halfm+1,n5,1),wk2d_d(2,2*n1halfm+1,n5,1))

! This sophisticated algorithm allows to decrease the memory needs.
  allocate(work1(2,n4,n5half1,n6),work2(2,n4,n5half1,n6))

  if(isign==1)then

!  Compute auxiliary arrays needed for FFTs, here forward FFT
   call sg_ctrig(n1,trig1,aft1,bef1,now1,one,ic1,ind1,mfac,mg)
   call sg_ctrig(n2,trig2,aft2,bef2,now2,one,ic2,ind2,mfac,mg)
   call sg_ctrig(n3,trig3,aft3,bef3,now3,one,ic3,ind3,mfac,mg)

!  Transfer fofg to the expanded fft box (only half of it)

!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index) SHARED(fofg,n1,n2,n3,nproc_fft,work1)
   do i3=1,n3
    do i2=1,n2half1
     index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1
      work1(1,i1,i2,i3)=fofg(1,i1+index)
      work1(2,i1,i2,i3)=fofg(2,i1+index)
     end do
    end do
   end do
!$OMP END PARALLEL DO

!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3)&
!$OMP&SHARED(now3,n1,n2half1,n4,n5half1,n6,ris,trig3,work1,work2)&
!$OMP&PRIVATE(i2)
  do i2=1,n2half1
   call sg_fftz(fftcache,mfac,mg,n4,n5half1,n6,n1,i2,i2,work1,work2,&
&   trig3,aft3,now3,bef3,ris,ind3,ic3)
  end do
!$OMP END PARALLEL DO

!  Loop over x-y planes

!$OMP PARALLEL PRIVATE(i1,i2,i3,index,wk2d_a,wk2d_b,wk2d_c,wk2d_d) &
!$OMP&SHARED(aft1,aft2,bef1,bef2,fftcache,fofg,fofr,ic1,ic2,ind1,ind2) &
!$OMP&SHARED(n1,n1half1,n1halfm,n2,n2half1,n3) &
!$OMP&SHARED(n4,n5,now1,now2,ris,trig1,trig2,work2)

  allocate(wk2d_a(2,n4,n5,1),wk2d_b(2,n4,n5,1))
  allocate(wk2d_c(2,2*n1halfm+1,n5,1),wk2d_d(2,2*n1halfm+1,n5,1))

!$OMP DO
   do i3=1,n3
    do i2=1,n2half1
     do i1=1,n1
      wk2d_c(1,i1,i2,1)=work2(1,i1,i2,i3)
      wk2d_c(2,i1,i2,1)=work2(2,i1,i2,i3)
     end do
    end do

    call sg_fftx(fftcache,mfac,mg,2*n1halfm+1,n5,1,n2half1,1,wk2d_c,wk2d_d,&
&    trig1,aft1,now1,bef1,ris,ind1,ic1)

!   Compute symmetric and antisymmetric combinations
    do i1=1,n1half1-1
     wk2d_a(1,i1,1,1)=wk2d_d(1,2*i1-1,1,1)
     wk2d_a(2,i1,1,1)=wk2d_d(1,2*i1  ,1,1)
    end do
!   If n1 odd, must add last data
    if((2*n1half1-2)/=n1)then
     wk2d_a(1,n1half1,1,1)=wk2d_d(1,n1,1,1)
     wk2d_a(2,n1half1,1,1)=0.0d0
    end if
    do i2=2,n2half1
     do i1=1,n1half1-1
      wk2d_a(1,i1,i2,1)=      wk2d_d(1,2*i1-1,i2,1)-wk2d_d(2,2*i1,i2,1)
      wk2d_a(2,i1,i2,1)=      wk2d_d(2,2*i1-1,i2,1)+wk2d_d(1,2*i1,i2,1)
      wk2d_a(1,i1,n2+2-i2,1)= wk2d_d(1,2*i1-1,i2,1)+wk2d_d(2,2*i1,i2,1)
      wk2d_a(2,i1,n2+2-i2,1)=-wk2d_d(2,2*i1-1,i2,1)+wk2d_d(1,2*i1,i2,1)
     end do
     if((2*n1half1-2)/=n1)then
      wk2d_a(1,n1half1,i2,1)=      wk2d_d(1,n1,i2,1)
      wk2d_a(2,n1half1,i2,1)=      wk2d_d(2,n1,i2,1)
      wk2d_a(1,n1half1,n2+2-i2,1)= wk2d_d(1,n1,i2,1)
      wk2d_a(2,n1half1,n2+2-i2,1)=-wk2d_d(2,n1,i2,1)
     end if
    end do

    call sg_ffty(fftcache,mfac,mg,n4,n5,1,1,n1halfm,1,1,wk2d_a,wk2d_b,&
&    trig2,aft2,now2,bef2,ris,ind2,ic2)

!   Take real part data from expanded box and put it in the original box.
    do i2=1,n2
     index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1half1-1
!    copy data
      fofr(2*i1-1+index)=wk2d_b(1,i1,i2,1)
      fofr(2*i1  +index)=wk2d_b(2,i1,i2,1)
     end do
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      fofr(n1+index)=wk2d_b(1,n1half1,i2,1)
     end if
    end do

!  End of loop over x-y planes
   end do
!$OMP END DO
   deallocate(wk2d_a,wk2d_b,wk2d_c,wk2d_d)

!$OMP END PARALLEL

  else if(isign==-1)then

!  Compute auxiliary arrays needed for FFTs, here backward FFT
   call sg_ctrig(n1,trig1,aft1,bef1,now1,-one,ic1,ind1,mfac,mg)
   call sg_ctrig(n2,trig2,aft2,bef2,now2,-one,ic2,ind2,mfac,mg)
   call sg_ctrig(n3,trig3,aft3,bef3,now3,-one,ic3,ind3,mfac,mg)

!  Treat first x-transform in x-y plane, and multiply
!  by overall normalization factor 1/nfftot

!  Loop over x-y planes

!$OMP PARALLEL PRIVATE(i1,i2,i3,index,wk2d_a,wk2d_b,wk2d_c,wk2d_d) &
!$OMP&SHARED(aft1,aft2,bef1,bef2,fftcache,fofr,ic1,ic2,ind1,ind2) &
!$OMP&SHARED(n1,n1half1,n1halfm,n2,n2half1,n3) &
!$OMP&SHARED(n4,n5,now1,now2,ris,trig1,trig2,work1,xnorm)

  allocate(wk2d_a(2,n4,n5,1),wk2d_b(2,n4,n5,1))
  allocate(wk2d_c(2,2*n1halfm+1,n5,1),wk2d_d(2,2*n1halfm+1,n5,1))

!$OMP DO
   do i3=1,n3
! write(6,*)'tache n :',OMP_GET_THREAD_NUM()
    do i2=1,n2
    index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1half1-1
!     copy and normalize data
      wk2d_a(1,i1,i2,1)=fofr(2*i1-1+index)*xnorm
      wk2d_a(2,i1,i2,1)=fofr(2*i1  +index)*xnorm
     end do
!    If n1 odd, must add last data
     if((2*n1half1-2)/=n1)then
      wk2d_a(1,n1half1,i2,1)=fofr(n1+index)*xnorm
      wk2d_a(2,n1half1,i2,1)=zero
     end if
    end do

    call sg_ffty(fftcache,mfac,mg,n4,n5,1,1,n1halfm,1,1,wk2d_a,wk2d_b,&
&    trig2,aft2,now2,bef2,ris,ind2,ic2)

!   Decompose symmetric and antisymmetric parts
    do i1=1,n1halfm
     wk2d_c(1,2*i1-1,1,1)=wk2d_b(1,i1,1,1)
     wk2d_c(2,2*i1-1,1,1)=0.0d0
     wk2d_c(1,2*i1,1,1)=wk2d_b(2,i1,1,1)
     wk2d_c(2,2*i1,1,1)=0.0d0
    end do

    do i2=2,n2half1
     do i1=1,n1halfm
      wk2d_c(1,2*i1-1,i2,1)=(wk2d_b(1,i1,i2,1)+wk2d_b(1,i1,n2+2-i2,1))*0.5d0
      wk2d_c(2,2*i1-1,i2,1)=(wk2d_b(2,i1,i2,1)-wk2d_b(2,i1,n2+2-i2,1))*0.5d0
      wk2d_c(1,2*i1,i2,1)= (wk2d_b(2,i1,i2,1)+wk2d_b(2,i1,n2+2-i2,1))*0.5d0
      wk2d_c(2,2*i1,i2,1)=-(wk2d_b(1,i1,i2,1)-wk2d_b(1,i1,n2+2-i2,1))*0.5d0
     end do
    end do

     call sg_fftx(fftcache,mfac,mg,2*n1halfm+1,n5,1,n2half1,1,wk2d_c,wk2d_d,&
&    trig1,aft1,now1,bef1,ris,ind1,ic1)

    do i2=1,n2half1
     do i1=1,n1
       work1(1,i1,i2,i3)=wk2d_d(1,i1,i2,1)
       work1(2,i1,i2,i3)=wk2d_d(2,i1,i2,1)
     end do
    end do

   end do
!$OMP END DO
  deallocate(wk2d_a,wk2d_b,wk2d_c,wk2d_d)
!$OMP END PARALLEL

!$OMP PARALLEL DO SHARED(aft3,bef3,fftcache,ind3,ic3)&
!$OMP&SHARED(now3,n1,n2half1,n4,n5half1,n6,ris,trig3,work1,work2)&
!$OMP&PRIVATE(i2)
  do i2=1,n2half1
   call sg_fftz(fftcache,mfac,mg,n4,n5half1,n6,n1,i2,i2,work1,work2,&
&   trig3,aft3,now3,bef3,ris,ind3,ic3)
  end do
!$OMP END PARALLEL DO

!Transfer fft output to the original fft box

!$OMP PARALLEL DO PRIVATE(i1,i2,i3,index,inter) &
!$OMP&SHARED(fofg,n1,n2,n2half1,n3,work2)
   do i3=1,n3
    do i2=1,n2half1
     index=n1*(i2-1+n2*(i3-1))
     do i1=1,n1
      fofg(1,i1+index)=work2(1,i1,i2,i3)
      fofg(2,i1+index)=work2(2,i1,i2,i3)
     end do
    end do
!   Complete missing values with complex conjugate
!   Inverse of ix is located at nx+2-ix , except for ix=1, for which it is 1.
    if(n2half1>2)then
     do i2=2,n2+1-n2half1
     index=n1*((n2+2-i2)-1)
      if(i3/=1)index=index+n1*n2*((n3+2-i3)-1)
      fofg(1,1+index)= work2(1,1,i2,i3)
      fofg(2,1+index)=-work2(2,1,i2,i3)
      do i1=2,n1
       fofg(1,n1+2-i1+index)= work2(1,i1,i2,i3)
       fofg(2,n1+2-i1+index)=-work2(2,i1,i2,i3)
      end do
     end do
    end if
   end do
!$OMP END PARALLEL DO

! Endif choice of isign
  end if

! deallocate(wk2d_a,wk2d_b,wk2d_c,wk2d_d,work1,work2)
  deallocate(work1,work2)

 end if

 call timab(260+tim_fourdp,2,tsec)

!DEBUG
!  write(6,*)' fourdp : exit '
!stop
!ENDDEBUG

end subroutine fourdp
!!***
