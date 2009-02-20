!{\src2tex{textfont=tt}}
!!****f* ABINIT/cppm2par
!! NAME
!! cppm2par
!!
!! FUNCTION
!!  Calculate plasmon-pole parameters using Hybertsen and Louie model (PRB 34, 5390 (1986))
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (RShaltaf, GMR, XG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iqiA(optional)= the index in the IBZ of the q-point where the ppmodel parameters have to be evaluated
!!  nqiA=number of irreducible points asked usually nqiA=niqbz. 
!!   It should be set to 1 if a single q-point is required (optional argument iqiA is needed)
!!  epsm1(npwc,npwc,nomega,nqiA)=symmetrized inverse dielectric matrix at nomega frequencies, and nq wavevectors
!!  gmet(3,3)=metric in reciprocal space
!!  ngfftf(18)=contain all needed information about the 3D fine FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  npwc=number of plane waves in epsm1
!!  nomega=number of frequencies (usually 2 but this plasmon-pole requires only the static symmetrized dielectric matrices
!!  Qmesh<BZ_mesh_type>=datatype gathering information on the q point sampling. see defs_datatypes.F90
!!    %nqibz=number of q points
!!    %ibz(3,nqibz)=irred q-points
!!  omega(nomega)=frequencies
!!  rhor(nfftf)=charge density on the real space FFT grid
!!  nfftf= total number of points in the fine FFT mesh  (for this processor)
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc,nqiA)= squared bare plasma frequencies
!!   \Omega^2_{G1 G2}(q) = 4\pi \frac {(q+G1).(q+G2)}/{|q+G1|^2} n(G1-G2)
!!
!!  omegatw(npwc,npwc,nqiA)= plasmon frequencies \tilde\omega_{G1 G2}(q) where:
!!  \tilde\omega^2_{G1 G2}(q) = 
!!    \frac {\Omega^2_{G1 G2}(q)} {\delta_{G1 G2}-\tilde\epsilon^{-1}_{G1 G2} (q, \omega=0)}
!!
!! NOTES 
!!  Note the use of intent(inout) since elements of data types are supposed 
!!  to be passed this routine
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      cggfft,fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cppm2par(paral_kgb,npwc,nqiA,nomega,epsm1,bigomegatwsq,omegatw,&
& ngfftf,gvec,gprimd,rhor,nfftf,Qmesh,gmet,&
& iqiA) ! Optional

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_15gw, except_this_one => cppm2par
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomega,npwc,nqiA,nfftf,paral_kgb
 integer,optional,intent(in) :: iqiA
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays
 integer,intent(in) :: gvec(3,npwc)
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3)
 real(dp),intent(inout) :: rhor(nfftf)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega,nqiA)
 complex(gwpc),intent(inout) :: bigomegatwsq(npwc,npwc,nqiA)
 complex(gwpc),intent(inout) :: omegatw(npwc,npwc,nqiA)

!Local variables-------------------------------
!scalars
 integer :: ig,ig_s,igp,igp_s,iq,istat,j,nimwp,tim_fourdp=2
 integer :: ngfft1,ngfft2,ngfft3
 real(dp) :: lambda,phi,x1
 logical,parameter :: use_symmetrized=.TRUE.,check_imppf=.FALSE.
 logical :: ltest
 character(len=500) :: msg
 type(MPI_type) :: mpi_enreg
!arrays
 integer,allocatable :: igfft(:,:)
 real(dp),pointer :: qibz(:,:)
 real(dp),allocatable :: qratio(:,:,:),qplusg(:),rhog_dp(:,:)
 complex(gwpc),allocatable :: omegatwsq(:,:) 
 complex(gwpc),allocatable :: rhog(:),rhogg(:,:),temp(:,:)  !MG these should be double precision TODO
!*************************************************************************

 qibz => Qmesh%ibz(1:3,1:Qmesh%nibz)
 if (PRESENT(iqiA)) then 
  call assert((nqiA==1),'nqiA should be 1',__FILE__,__LINE__)
  ltest=(iqiA>0.and.iqiA<=Qmesh%nibz)
  call assert(ltest,'iqiA out of range',__FILE__,__LINE__)
  qibz => Qmesh%ibz(1:3,iqia:iqia)
 end if
 !
 ! === Calculate qratio(npwec,npvec,nqiA) = (q+G).(q+Gp)/|q+G|^2 ===
 allocate(qratio(npwc,npwc,nqiA),stat=istat) 
 if (istat/=0) call memerr('cppm2par','qratio',npwc**2*nqiA,'dp')
 call cqratio(npwc,gvec,nqiA,qibz,gmet,gprimd,qratio)
 !
 ! === Compute the density in G space rhor(R)--> rhog(G) ===
 allocate(rhog_dp(2,nfftf),rhog(nfftf))
 ngfft1=ngfftf(1) ; ngfft2=ngfftf(2) ; ngfft3=ngfftf(3)

 call fourdp(1,rhog_dp,rhor,-1,mpi_enreg,nfftf,ngfftf,paral_kgb,0)
 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))
 !
 ! Calculate the FFT index of each (G-Gp) vector and assign 
 ! the value of the correspondent density simultaneously
 allocate(igfft(npwc,npwc),rhogg(npwc,npwc))
 call cggfft(npwc,ngfft1,ngfft2,ngfft3,gvec,igfft)
 do ig=1,npwc
  do igp=1,npwc
   rhogg(ig,igp)=rhog(igfft(ig,igp))
  end do
 end do
 rhogg(:,:)=four_pi*rhogg(:,:)
 deallocate(igfft,rhog_dp,rhog)
 !
 ! Calculate GPP parameters
 ! unsymmetrized epsm1 -> epsm1=|q+Gp|/|q+G|*epsm1
 !
 allocate(qplusg(npwc),temp(npwc,npwc),omegatwsq(npwc,npwc),stat=istat)
 if (istat/=0) stop 'cppm2par out of memory in qplusg'

 do iq=1,nqiA
  temp(:,:)=-epsm1(:,:,1,iq)
  ! 
  ! RS still not obvious for me whether one shall use the symmetrized inverse DM or the unsymmetrized one
  ! the default here is to use the symmetrized one, I must discuss this with XG
  ! 
  ! MG it turns out that using the symmetrized inverse DM in the plasmon-pole
  ! equations give the same results for the squared plasmon frequencies omegatwsq while the 
  ! squared bare plasma frequencies bigomegatwsq related to the symmetrized dielectric matrix 
  ! are obtained multipling by |q+G1|/|q+G2|   
  ! 
  if (.not.use_symmetrized) then
   call cvc(nqiA,iq,qibz,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q
   do ig=1,npwc
    do igp=1,npwc
     temp(ig,igp)=qplusg(igp)/qplusg(ig)*temp(ig,igp)
    end do
   end do
  end if

  nimwp=0
  do ig=1,npwc
   temp(ig,ig)=temp(ig,ig)+one
   do igp=1,npwc
    bigomegatwsq(ig,igp,iq)=rhogg(ig,igp)*qratio(ig,igp,iq)
    omegatwsq(ig,igp)=bigomegatwsq(ig,igp,iq)/temp(ig,igp)
    !   
    ! Set omegatw to any arbitrary number to avoid dealing with undefined numbers like (INF)
    ! simply ignore all cases of omegatw with imaginary values
    ! in principle these correspond to cases where the imaginary part of epsm1 does not have
    ! a well defined peak. The imaginary part of epsm1 in these cases oscillates  with a small amplitude
    ! since the amplitude A_GGpr=-pi/2*bigomegatwsq/omegatw, 
    ! it follows that bigomegatwsq shall be set to zero for these cases
    !   
    if ( REAL(omegatwsq(ig,igp))<= tol12 .or. AIMAG(omegatwsq(ig,igp))**2*tol12> REAL(omegatwsq(ig,igp))**2) then
     bigomegatwsq(ig,igp,iq)=(0.,0.)
     omegatw(ig,igp,iq)=(ten,0.)
     nimwp=nimwp+1
     if (check_imppf) then 
      write(msg,'(a,i3,2i8)')' imaginary plasmon frequency at : ',iq,ig,igp
      call wrtout(std_out,msg,'COLL')
     end if 
    else 
     ! this part has been added to deal with systems without inversion symmetry
     ! this new implementation gives the same results as the previous one if 
     ! omegatwsq is a pure real number and has the advantage of being an improved
     ! approach for systems without an inversion center.
     lambda=ABS(omegatwsq(ig,igp))
     phi=ATAN(AIMAG(omegatwsq(ig,igp))/REAL(omegatwsq(ig,igp)))
     omegatw(ig,igp,iq)=SQRT(lambda/COS(phi))
     bigomegatwsq(ig,igp,iq)=bigomegatwsq(ig,igp,iq)*(1.-(0.,1.)*TAN(phi))
     ! Uncomment the following line and comment the previous to restore the old version.
     !omegatw(ig,igp,iq)=sqrt(real(omegatwsq(ig,igp)))
    end if
   end do
  end do
  write(msg,'(a,3f12.6,a,i8,a,i8)')&
&  ' at q-point : ',qibz(:,iq),&
&  ' number of imaginary plasmonpole frequencies = ',nimwp,' / ',npwc**2
  call wrtout(std_out,msg,'COLL')
 end do !iq

 write(msg,'(2a,f12.8,2a,3i5)')ch10,&
 ' cppm2par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw(:,:,:)))*Ha_eV,ch10,&
 '            omega twiddle min location = ',MINLOC(ABS(omegatw(:,:,:)))
 call wrtout(std_out,msg,'COLL')

 deallocate(omegatwsq,rhogg,temp,qplusg,qratio)

end subroutine cppm2par
!!***
