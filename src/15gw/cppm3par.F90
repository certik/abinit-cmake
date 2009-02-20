!{\src2tex{textfont=tt}}
!!****f* ABINIT/cppm3par
!! NAME
!! cppm3par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters using the von Linden-Horsh model (PRB 37, 8351, 1988)
!! (see also Pag 22 of Quasiparticle Calculations in Solids. Aulbur et. al)  
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
!!  nqiA=number of irreducible points asked usually nqiA=niqbz, 
!!   It should be set to 1 if a single q-point is required (optional argument iqiA is needed)
!! epsm1(npwc,npwc,nomega,nqiA)= symmetrized inverse dielectric 
!!  matrix at nomega frequencies, and nqiA wavevectors
!! MPI_enreg=information about MPI parallelization
!! ngfftf(18)=contain all needed information about 3D fine FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! npwc=number of plane waves in epsm1
!! qratio=(q+G1).(q+G2)/(|q+G1|.|q+G2|)
!! Qmesh<BZ_mesh_type>=datatype gathering information on the q point sampling. see defs_datatypes.F90
!!   %nqibz=number of irreducible q-points
!!   %qibz(3,nqibz)=irreducible q-points.
!! rho(nfftf)=charge density on the real space FFT grid
!! nfftf=number of points in the FFT grid (for this processor)
!! gvec(3,npwc)= G vectors in reduced coordinates
!!
!! OUTPUT
!!  omegatw(npwc,npwc,nqiA)= plasmon pole positions
!!  bigomegatwsq(npwc,npwc,nqiA)=(E_{q,ii}^{-1}-1)*omegatw
!!   where E^{-1} is the eigenvalue of the inverse dielectric matrix
!!  eigtot(npwc,npwc,nqiA)=the eigvectors of the symmetrized inverse dielectric matrix 
!!   (first index for G, second index for bands)
!!
!! NOTES 
!!  Note the use of intent(inout) since elements of data types are supposed 
!!  to be passed this routine
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      cggfft,chpev,fourdp,leave_new,wrtout,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cppm3par(paral_kgb,npwc,nqiA,nomega,epsm1,bigomegatwsq,omegatw,&
& ngfftf,gvec,gprimd,rho,nfftf,eigtot,Qmesh,&
& iqiA) ! Optional

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_15gw, except_this_one => cppm3par
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,nomega,npwc,nqiA,paral_kgb
 integer,intent(in),optional :: iqiA
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays
 integer,intent(in) :: gvec(3,npwc),ngfftf(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout) :: rho(nfftf)
 complex(gwpc),intent(in) :: epsm1(npwc,npwc,nomega,nqiA)
 complex(gwpc),intent(inout) :: bigomegatwsq(npwc,1,nqiA),eigtot(npwc,npwc,nqiA)
 complex(gwpc),intent(inout) :: omegatw(npwc,1,nqiA)

!Local variables-------------------------------
!TODO these should be dp
!scalars
 integer :: idx,ierr,ig,igp,ii,iq,istat,jj,ngfft1,ngfft2,ngfft3
 real(dp) :: num,qpg_dot_qpgp
 complex(dpc) :: conjg_eig
 logical :: ltest
 character(len=500) :: msg
 type(MPI_type) :: MPI_enreg
!arrays
 integer,allocatable :: igfft(:,:)
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3)
 real(dp),allocatable :: eigval(:),qplusg(:),rhog_dp(:,:),zhpev2(:)
 real(dp),pointer :: qibz(:,:)
 complex(dpc),allocatable :: eigvec(:,:),matr(:),mm(:,:,:),rhog(:),rhogg(:,:)
 complex(dpc),allocatable :: zhpev1(:),zz(:,:)

!*************************************************************************

#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

 qibz => Qmesh%ibz(1:3,1:Qmesh%nibz)
 if (PRESENT(iqiA)) then 
  call assert((nqiA==1),'nqiA should be 1',__FILE__,__LINE__)
  ltest=(iqiA>0.and.iqiA<=Qmesh%nibz)
  call assert(ltest,'iqiA out of range',__FILE__,__LINE__)
  qibz => Qmesh%ibz(1:3,iqia:iqia)
 end if

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)

 ngfft1=ngfftf(1) ; ngfft2=ngfftf(2) ; ngfft3=ngfftf(3)

 allocate(rhog_dp(2,nfftf),rhog(nfftf),stat=istat)      ; if (istat/=0) stop 'rhog out of memory'
 allocate(igfft(npwc,npwc),rhogg(npwc,npwc),stat=istat) ; if (istat/=0) stop 'rhogg out of memory'
 !
 ! === Compute the density in G space rhog(r)--> rho(G) ===
 ! FIXME this has to be fixed, rho(G) should be passed instead of doing FFT for each q
 ! Moreover MPI_enreg is local ????? 
 call fourdp(1,rhog_dp,rho,-1,MPI_enreg,nfftf,ngfftf,paral_kgb,0)

 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))
 !
 ! Calculate the FFT index of each (G-Gp) vector and assign the value
 ! of the correspondent density simultaneously
 call cggfft(npwc,ngfft1,ngfft2,ngfft3,gvec,igfft)

 do ig=1,npwc
  do igp=1,npwc
   if (igfft(ig,igp)>nfftf) then
    write (msg,'(4a)')ch10,&
&    ' cppm3par : BUG- ',ch10,&
&    ' cannot find rho(G-Gpr) for some G, Gpr '
    call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
   end if
   rhogg(ig,igp)=rhog(igfft(ig,igp))
  end do
 end do
 !
 ! mm(G,Gp) = (q+G) \cdot (q+Gp) n(G-Gp)
 allocate(mm(npwc,npwc,nqiA),stat=istat) ; if (istat/=0) stop 'mm out of memory'

 do iq=1,nqiA
  do ig=1,npwc
   if (ALL(ABS(qibz(:,iq))<1.0e-3)) then
    ! To be discussed with Riad, here we should use the small q 
    ! to be consistent and consider the limit q-->0
    gpq(1)=gvec(1,ig)
    gpq(2)=gvec(2,ig)
    gpq(3)=gvec(3,ig)
   else
    gpq(1)=gvec(1,ig)+qibz(1,iq)
    gpq(2)=gvec(2,ig)+qibz(2,iq)
    gpq(3)=gvec(3,ig)+qibz(3,iq)
   end if
   do igp=1,npwc
    if (ALL(ABS(qibz(:,iq))<1.0e-3)) then
     gppq(1)=gvec(1,igp)
     gppq(2)=gvec(2,igp)
     gppq(3)=gvec(3,igp)
    else
     gppq(1)=gvec(1,igp)+qibz(1,iq)
     gppq(2)=gvec(2,igp)+qibz(2,iq)
     gppq(3)=gvec(3,igp)+qibz(3,iq)
    end if
    qpg_dot_qpgp=zero
    do ii=1,3
     qpg_dot_qpgp=qpg_dot_qpgp+&
&     ( gpq(1)*b1(ii) +gpq(2)*b2(ii) +gpq(3)*b3(ii))*&
&     (gppq(1)*b1(ii)+gppq(2)*b2(ii)+gppq(3)*b3(ii))
    end do
    mm(ig,igp,iq)=rhogg(ig,igp)*qpg_dot_qpgp
   end do !igp
  end do !ig
 end do !iq
 deallocate(rhog_dp,rhog,igfft)
 ! === Now we have rhogg,rho0 ===
 !
 ! Calculate the dielectric matrix eigenvalues and vectors
 ! Use only the static epsm1 i.e., only the w=0 part (eps(:,:,1,:))
 allocate(eigval(npwc),eigvec(npwc,npwc),stat=istat) ! eigenvalues and vectors of DM
 if (istat/=0) stop 'eigvec out of memory'
 allocate(zz(npwc,nqiA),stat=istat) ; if (istat/=0) stop 'zz of memory'
 zz(:,:)=czero
 allocate(qplusg(npwc))

 do iq=1,nqiA
  !
  ! Store the susceptibility matrix in upper mode before calling zhpev for each iq value
  allocate(matr(npwc*(npwc+1)/2),stat=istat) ; if(istat/=0) stop 'matr of memory'

  idx=1
  do ii=1,npwc
   do jj=1,ii
    matr(idx)=epsm1(jj,ii,1,iq) ; idx=idx+1
   end do
  end do

  allocate(zhpev2(3*npwc-2),zhpev1(2*npwc-1),stat=istat)
  if (istat/=0) stop 'zhpev1 of memory' ! working arrays for lapack
#if defined T3E
  call CHPEV('V','U',npwc,matr,eigval,eigvec,npwc,zhpev1,zhpev2,ierr)
#else
  call ZHPEV('V','U',npwc,matr,eigval,eigvec,npwc,zhpev1,zhpev2,ierr)
#endif
  deallocate(matr,zhpev2,zhpev1)

  if (ierr<0) then
   write (msg,'(5a,i4,2a)')ch10,&
&   ' cppm3par : ERROR- ',ch10,&
&   '  failed to calculate the eigenvalues and eigenvectors of the dielectric matrix ',ch10,&
&   ierr*(-1),' th argument in the matrix has an illegal value',ch10
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if

  if (ierr>0) then
   write(msg,'(6a,i4,5a)')ch10,&
&   ' cppm3par : ERROR- ',ch10,&
&   '  failed to calculate the eigenvalues and eigenvectors of the dielectric matrix ',ch10,&
&   '  the algorithm failed to converge; ', ierr,ch10,&
&   '  off-diagonal elements of an intermediate tridiagonal form ',ch10,&
&   '  did not converge to zero',ch10
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if
  !
  ! Calculate the PPM parameters and the eigenpotentials needed for 
  ! the calculation of the generalized overlap matrix
  ! Note: the eigenpotentials has to be calculated on the FFT (G-Gp) index
  !
  ! Save eigenvectors of \tilde\epsilon^{-1}
  ! MG well it is better to save \Theta otherwise 
  ! we have to calculare \Theta for each band, spin, k-point but oh well
  eigtot(:,:,iq)=eigvec(:,:)

  call cvc(nqiA,iq,qibz,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q
  !
  ! Basic Equation:
  ! 
  ! \Theta_{q,ii}(G)=\Psi_{q,ii}(G)/|q+G|
  ! where \Psi_{q,ii}(G) is the eigenvector of \tilde\epsilon^{-1} 

  ! \tilde\omega_{ii,q}^2= 4\pi (1-eigenval(ii,q))) 
  ! \sum_{G,Gp} \Theta^*_{q,ii}(G) (q+G)\cdot(q+Gp) n(G-Gp) \Theta_{q,ii}(Gp) 

  do ii=1,npwc !DM band
   ! Calculate \Theta_{q,ii}(G)
   ! why the first element is not modified? if the problem is the small value of qplusg(1)
   ! we could multiply by sqrt(mod((q+G)(q+G'))) and then add the sing at the end 
   if(iq==1)then
    eigvec(2:,ii)=eigvec(2:,ii)/qplusg(2:)
   else
    eigvec(:,ii)=eigvec(:,ii)/qplusg(:)
   end if
   do ig=1,npwc
    conjg_eig=CONJG(eigvec(ig,ii))
    do igp=1,npwc
     if(iq==1 .and. ig==1 .and. igp==1)then
      zz(ii,iq)=zz(ii,iq)+conjg_eig*rhogg(ig,igp)*eigvec(igp,ii)
     else
      zz(ii,iq)=zz(ii,iq)+conjg_eig*mm(ig,igp,iq)*eigvec(igp,ii)
     end if
    end do
   end do

   num=1-eigval(ii)
   if (num<=zero) then
!   here I think we should set bigomegatwsq=0 and omegatw to an arbitrary value
!   maybe we can output a warning TO BE discussed with Riad 
    if (ABS(num)<1.0d-4) then
     num=1.0d-5
    else
     write(msg,'(5a)')ch10,&
&     ' cppm3par : BUG - ',ch10,&
&     ' one or more imaginary plasmon pole energies',ch10
     call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
    end if
   end if

   omegatw(ii,1,iq)=SQRT(4*pi*REAL(zz(ii,iq))/(num))
   ! this should be \alpha = 2\pi omegatw * (1-eigenval) 
   ! MG check this, in the review I found a factor 2\pi, maybe it is reintroduced later
   bigomegatwsq(ii,1,iq)=num*omegatw(ii,1,iq)
  end do
 end do !iq

 deallocate(rhogg,mm,eigval,zz,eigvec,qplusg)

 write(msg,'(2a,f12.8,2a,3i5)')ch10,&
 ' cppm3par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw(:,:,:)))*Ha_eV,ch10,&
 '            omega twiddle min location = ',MINLOC(ABS(omegatw(:,:,:)))
 call wrtout(std_out,msg,'COLL')

end subroutine cppm3par
!!***
