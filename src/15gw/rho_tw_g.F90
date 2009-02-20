!{\src2tex{textfont=tt}}
!!****f* ABINIT/rho_tw_g
!! NAME
!! rho_tw_g
!!
!! FUNCTION
!! Calculate rhotwg(G)=<wfn1|exp(-i(q+G).r)|wfn2>
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dim_rtwg=Define the size of the output array rhotwg
!!   === for nspinor==1 ===
!!    dim_rtwg=1
!!   === for nspinor==2 ===
!!    dim_rtwg=2 if only <up|up>, <dwn|dwn> matrix elements are required
!!    dim_rtwg=4 for <up|up>, <dwn|dwn>, <up|dwn> and <dwn|up>.
!! map2sphere= 1 to retrieve Fourier components indexed according to igfftg0.
!!             0 to retrieve Fourier components indexed according to the FFT box.
!!               NOTE: If map2sphere==0 npwvec must be equal to nr
!! igfftg0(npwvec*map2sphere)=index of G-G_o in the FFT array for each G in the sphere.
!! i1=1 if kbz1 = Sk1, 2 if kbz1 = -Sk_1 (k_1 is in the IBZ)
!! i2=1 if kbz2 = Sk2, 2 if kbz2 = -Sk_2 (k_2 is in the IBZ)
!! ktabr1(nr),ktabr2(nr)= tables R^-1(r-t) for the two k-points
!! ktabp1,ktabp2 = phase factors for non-simmorphic symmetries e^{-i 2\pi kbz.\tau} 
!! MPI_enreg=Information about MPI parallelization
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npwvec=number of plane waves (in the sphere if map2sphere==1, in the FFT box if map2sphere==1)
!! nr=number of FFT grid points
!! nspinor=number of spinorial components.
!! spinrot1(4),spinrot2(4)=components of the spinor rotation matrix :
!!  spinrot(1)=$\cos\phi/2$
!!  spinrot(2)=$\sin\phi/2 \times u_x$
!!  spinrot(3)=$\sin\phi/2 \times u_y$
!!  spinrot(4)=$\sin\phi/2 \times u_z$
!!   where $\phi$ is the angle of rotation, and
!!   $(u_x,u_y,u_z)$ is the normalized direction of the rotation axis
!! tim_fourdp=1 if called from inside screening 
!!            2 if called from inside sigma
!! wfn1(nr),wfn2(nr)=the two wavefunctions (periodic part)
!!
!! OUTPUT
!! rhotwg(npwsigx)=density of a pair of states, in reciprocal space
!!
!! PARENTS
!!      cchi0,cchi0q0,csigme
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rho_tw_g(paral_kgb,nspinor,npwvec,nr,ngfft,map2sphere,igfftg0,&
& wfn1,i1,ktabr1,ktabp1,spinrot1,&
& wfn2,i2,ktabr2,ktabp2,spinrot2,&
& dim_rtwg,rhotwg,tim_fourdp,MPI_enreg)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : j_dpc
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: paral_kgb,i1,i2,npwvec,nr,tim_fourdp,nspinor,dim_rtwg,map2sphere
 complex(gwpc),intent(in) :: ktabp1,ktabp2
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 integer,intent(in) :: igfftg0(npwvec*map2sphere),ngfft(18)
 integer,intent(in) :: ktabr1(nr),ktabr2(nr)
 real(dp),intent(in) :: spinrot1(4),spinrot2(4)
 complex(gwpc),intent(in) :: wfn1(nr*nspinor),wfn2(nr*nspinor)
 complex(gwpc),intent(out) :: rhotwg(npwvec*dim_rtwg)

!Local variables-------------------------------
!scalars
 integer :: ig,ir,ir1,ir2,igfft,iab,spad1,spad2,spad0,ispinor
 real(dp) :: wfi1,wfi2,wfr1,wfr2
 real(dp) :: kpi1,kpi2,kpr1,kpr2,tt
 complex(gwpc) :: u1a,u1b,u2a,u2b
!arrays
 integer :: spinor_pad(2,4) 
 real(dp) :: phase(2)
 real(dp),allocatable :: rhotw_dp(:,:),rhotwg_dp(:,:)
 complex(dpc) :: spinrot_mat1(2,2),spinrot_mat2(2,2)
 complex(gwpc),allocatable :: cwavef1(:),cwavef2(:)

! *************************************************************************

#if defined DEBUG_MODE
 if (map2sphere==0) then 
  call assert((npwvec==nr),'If map2sphere==0, npwvec must be equal to nr')
 end if
#endif

 SELECT CASE (nspinor)

  CASE (1)
  ! Form rho-twiddle(r)=u_1^*(r,b1,kbz1) u_2^*(r,b2,kbz2), to account for symmetries:
  ! u(r,b,kbz)=e^{-2i\pi kibz.(R^{-1}t} u (R{^-1}(r-t),b,kibz) 
  !           =e^{+2i\pi kibz.(R^{-1}t} u*({R^-1}(r-t),b,kibz) for time-reversal
  !
  ! === Precalculate the phase due to non-symmporphic translations ===
  kpr1=REAL(ktabp1) ; kpi1=AIMAG(ktabp1)
  kpr2=REAL(ktabp2) ; kpi2=AIMAG(ktabp2) 

  allocate(rhotwg_dp(2,nr),rhotw_dp(2,nr))

  SELECT CASE (i1)

  CASE (1)
   SELECT CASE (i2)
   CASE (1)
    phase(1)=kpr1*kpr2+kpi1*kpi2 
    phase(2)=kpr1*kpi2-kpi1*kpr2 
    do ir=1,nr
     ir1=ktabr1(ir) ; ir2=ktabr2(ir)
     wfr1=REAL(wfn1(ir1)) ; wfi1=AIMAG(wfn1(ir1))
     wfr2=REAL(wfn2(ir2)) ; wfi2=AIMAG(wfn2(ir2))
     rhotw_dp(1,ir)=wfr1*wfr2 + wfi1*wfi2
     rhotw_dp(2,ir)=wfr1*wfi2 - wfi1*wfr2
    end do
   CASE (2)
    phase(1)= kpr1*kpr2-kpi1*kpi2 
    phase(2)=-kpr1*kpi2-kpi1*kpr2 
    do ir=1,nr
     ir1=ktabr1(ir) ; ir2=ktabr2(ir)
     wfr1=REAL(wfn1(ir1)) ; wfi1=AIMAG(wfn1(ir1))
     wfr2=REAL(wfn2(ir2)) ; wfi2=AIMAG(wfn2(ir2))
     rhotw_dp(1,ir)= wfr1*wfr2 - wfi1*wfi2
     rhotw_dp(2,ir)=-wfr1*wfi2 - wfi1*wfr2
    end do
   CASE DEFAULT
    STOP 'Wrong i2(1)'
   END SELECT

  CASE (2)
   SELECT CASE (i2)
   CASE (1)
    phase(1)= kpr1*kpr2-kpi1*kpi2 
    phase(2)= kpr1*kpi2+kpi1*kpr2 
    do ir=1,nr
     ir1=ktabr1(ir) ; ir2=ktabr2(ir)
     wfr1=REAL(wfn1(ir1)) ; wfi1=AIMAG(wfn1(ir1))
     wfr2=REAL(wfn2(ir2)) ; wfi2=AIMAG(wfn2(ir2))
     rhotw_dp(1,ir)= wfr1*wfr2 - wfi1*wfi2
     rhotw_dp(2,ir)= wfr1*wfi2 + wfi1*wfr2
    end do
   CASE (2)
    phase(1)= kpr1*kpr2+kpi1*kpi2 
    phase(2)=-kpr1*kpi2+kpi1*kpr2 
    do ir=1,nr
     ir1=ktabr1(ir) ; ir2=ktabr2(ir)
     wfr1=REAL(wfn1(ir1)) ; wfi1=AIMAG(wfn1(ir1))
     wfr2=REAL(wfn2(ir2)) ; wfi2=AIMAG(wfn2(ir2))
     rhotw_dp(1,ir)= wfr1*wfr2 + wfi1*wfi2
     rhotw_dp(2,ir)=-wfr1*wfi2 + wfi1*wfr2
    end do
   CASE DEFAULT
    STOP 'Wrong i2(2)'
   END SELECT

  CASE DEFAULT 
   STOP 'Wrong i1'
  END SELECT
  !
  ! === Apply the phase (maybe we can skip this part if the phase is 1) ===
  do ir=1,nr
   tt=rhotw_dp(1,ir)
   rhotw_dp(1,ir)=phase(1)*tt-phase(2)*rhotw_dp(2,ir)
   rhotw_dp(2,ir)=phase(1)*rhotw_dp(2,ir)+phase(2)*tt
  end do
  !
  ! === Set up rho-twiddle(G-G0) for matrix in PW order, normalising for FFT ===
  call fourdp(2,rhotwg_dp,rhotw_dp,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)

  if (map2sphere==1) then
   do ig=1,npwvec
    igfft=igfftg0(ig)
    rhotwg(ig)=CMPLX(rhotwg_dp(1,igfft),rhotwg_dp(2,igfft),kind=gwpc)
   end do
  else if (map2sphere==0) then 
   rhotwg(:)=CMPLX(rhotwg_dp(1,:),rhotwg_dp(2,:),kind=gwpc)
  else 
   STOP 'Wrong value for map2sphere'
  end if

  deallocate(rhotwg_dp,rhotw_dp)

 CASE (2)

  allocate(cwavef1(nr*nspinor),cwavef2(nr*nspinor))

  ! === Apply Time-reversal if required ===
  ! \psi_{-k}^1 =  (\psi_k^2)^*
  ! \psi_{-k}^2 = -(\psi_k^1)^*
  if (i1==1) then 
   cwavef1(:)=wfn1(:) 
  else if (i1==2) then
   cwavef1(1:nr)     = CONJG(wfn1(nr+1:2*nr))
   cwavef1(nr+1:2*nr)=-CONJG(wfn1(1:nr))
  else 
   STOP 'Wrong i1 in spinor'
  end if

  if (i2==1) then 
   cwavef2(:)=wfn2(:) 
  else if (i2==2) then
   cwavef2(1:nr)     = CONJG(wfn2(nr+1:2*nr))
   cwavef2(nr+1:2*nr)=-CONJG(wfn2(1:nr))
  else 
   STOP 'Wrong i2 in spinor'
  end if

  ! === Rotate wavefunctions in r-space ===
  do ispinor=1,nspinor
   spad0=(ispinor-1)*nr
   do ir=1,nr
    ir1=ktabr1(ir) ; ir2=ktabr2(ir)
    cwavef1(ir+spad0) = cwavef1(ir1+spad0)*ktabp1
    cwavef2(ir+spad0) = cwavef2(ir2+spad0)*ktabp2
   end do 
  end do !ispinor

  ! === Rotation in spinor space ===
  !spinrots1=spinrot1(1) ; spinrots2=spinrot2(1)
  !spinrotx1=spinrot1(2) ; spinrotx2=spinrot2(2)
  !spinroty1=spinrot1(3) ; spinroty2=spinrot2(3)
  !spinrotz1=spinrot1(4) ; spinrotz2=spinrot2(4)
  spinrot_mat1(1,1)= spinrot1(1) + j_dpc*spinrot1(4)
  spinrot_mat1(1,2)= spinrot1(3) + j_dpc*spinrot1(2)
  spinrot_mat1(2,1)=-spinrot1(3) + j_dpc*spinrot1(2)
  spinrot_mat1(2,2)= spinrot1(1) - j_dpc*spinrot1(4)

  spinrot_mat2(1,1)= spinrot2(1) + j_dpc*spinrot2(4)
  spinrot_mat2(1,2)= spinrot2(3) + j_dpc*spinrot2(2)
  spinrot_mat2(2,1)=-spinrot2(3) + j_dpc*spinrot2(2)
  spinrot_mat2(2,2)= spinrot2(1) - j_dpc*spinrot2(4)

  do ir=1,nr
   !ar=wavefspinor(1,ir)
   !ai=wavefspinor(2,ir)
   !br=wavefspinor(1,npw2+ir)
   !bi=wavefspinor(2,npw2+ir)
   u1a=cwavef1(ir) 
   u1b=cwavef1(ir+nr) 
   cwavef1(ir)   =spinrot_mat1(1,1)*u1a+spinrot_mat1(1,2)*u1b
   cwavef1(ir+nr)=spinrot_mat1(2,1)*u1a+spinrot_mat1(2,2)*u1b
   u2a=cwavef2(ir) 
   u2b=cwavef2(ir+nr) 
   cwavef2(ir)   =spinrot_mat2(1,1)*u2a+spinrot_mat2(1,2)*u2b
   cwavef2(ir+nr)=spinrot_mat2(2,1)*u2a+spinrot_mat2(2,2)*u2b
   !wavefspinor(1,ir)     = spinrots*ar-spinrotz*ai +spinroty*br-spinrotx*bi
   !wavefspinor(2,ir)     = spinrots*ai+spinrotz*ar +spinroty*bi+spinrotx*br
   !wavefspinor(1,npw2+ir)=-spinroty*ar-spinrotx*ai +spinrots*br+spinrotz*bi
   !wavefspinor(2,npw2+ir)=-spinroty*ai+spinrotx*ar +spinrots*bi-spinrotz*br
  end do

  spinor_pad(:,:)=RESHAPE((/0,0,nr,nr,0,nr,nr,0/),(/2,4/))
  allocate(rhotwg_dp(2,nr),rhotw_dp(2,nr))

  do iab=1,dim_rtwg
   spad1=spinor_pad(1,iab)
   spad2=spinor_pad(2,iab)

   do ir=1,nr
    ir1=ir+spad1 ; ir2=ir+spad2
    wfr1=REAL(cwavef1(ir1)) ; wfi1=AIMAG(cwavef1(ir1))
    wfr2=REAL(cwavef2(ir2)) ; wfi2=AIMAG(cwavef2(ir2))
    rhotw_dp(1,ir)= wfr1*wfr2 + wfi1*wfi2
    rhotw_dp(2,ir)= wfr1*wfi2 - wfi1*wfr2
   end do

   ! === Set up rho-twiddle(G-G0) for matrix in PW order, normalising for FFT ===
   call fourdp(2,rhotwg_dp,rhotw_dp,-1,MPI_enreg,nr,ngfft,paral_kgb,tim_fourdp)

   spad0=(iab-1)*npwvec

   if (map2sphere==1) then
    do ig=1,npwvec
     igfft=igfftg0(ig)
     rhotwg(ig+spad0)=CMPLX(rhotwg_dp(1,igfft),rhotwg_dp(2,igfft),kind=gwpc)
    end do
   else if (map2sphere==0) then 
    rhotwg(spad0+1:spad0+npwvec)=CMPLX(rhotwg_dp(1,:),rhotwg_dp(2,:),kind=gwpc)
   else 
    STOP 'Wrong value for map2sphere'
   end if

  end do !iab

  deallocate(cwavef1,cwavef2)
  deallocate(rhotwg_dp,rhotw_dp)

 CASE DEFAULT 
  STOP 'Wrong nspinor'
 END SELECT

end subroutine rho_tw_g
!!***
