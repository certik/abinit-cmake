!{\src2tex{textfont=tt}}
!!****f* ABINIT/cppm4par
!! NAME
!! cppm4par
!!
!! FUNCTION
!! Calculate the plasmon-pole parameters using Engel and Farid model (PRB47,15931,1993)
!! See also Quasiparticle Calculations in Solids, Aulbur et al. (pag. 23)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (RShaltaf, GMR, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  iqiA(optional)= the index in the IBZ of the q-point where the ppmodel parameters have to be evaluated
!!  nqiA=number of irreducible points asked usually nqiA=niqbz. 
!!   It should be set to 1 if a single q-point is required (optional argument iqiA is needed)
!!  epsm1(npwc,npwc,nomega,nqiA)=symmetrized inverse dielectric matrix at 
!!   nomega frequencies, and nqiA wavevectors
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  ngfftf(18)=contain all needed information about 3D fine FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  npwc=number of plane waves in epsm1
!!  nomega=number of frequencies (usually 2 but this model requires only \omega=0)
!!  Qmesh<BZ_mesh_type>=datatype gathering information on the q point sampling. see defs_datatypes.F90
!!    %nqibz=number of qibz vectors in the IBZ
!!    %qibz(3,nqibz)=irreducible q-points.
!!  rho(nfftf)=charge density on the real space FFT grid
!!  gvec(3,npwc)=G vectors in reduced coordinated 
!!
!! OUTPUT
!!  bigomegatwsq(npwc,npwc,nqiA)=plasmon-pole strength
!!  omegatw(npwc,npwc,nqiA)=plasmon-pole frequencies
!!
!! NOTES 
!!  Note the use of intent(inout) since elements of data types are supposed 
!!  to be passed this routine
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cppm4par(paral_kgb,npwc,nqiA,epsm1,nomega,bigomegatwsq,omegatw,ngfftf,gvec,gprimd,rho,nfftf,Qmesh,&
& iqia) ! Optional

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOLQ0
 use m_numeric_tools, only : is_zero
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12geometry
 use interfaces_15gw, except_this_one => cppm4par
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
 complex(gwpc),intent(inout) :: bigomegatwsq(npwc,npwc,nqiA),omegatw(npwc,1,nqiA)

!Local variables-------------------------------
!scalars
 integer :: idx,ierr,ig,igp,ii,iq,istat,jj,lowork,ngfft1,ngfft2,ngfft3
 real(dp) :: qpg_dot_qpgp
 logical :: ltest
 character(len=500) :: msg
 character(len=80) :: bar
 type(MPI_type) :: MPI_enreg
!arrays
 integer,allocatable :: igfft(:,:)
 real(dp) :: b1(3),b2(3),b3(3),gppq(3),gpq(3)
 real(dp),allocatable :: eigval(:),eigval1(:),qplusg(:),rhog_dp(:,:),rwork(:)
 real(dp),allocatable :: zhpev2(:)
 real(dp),pointer :: qibz(:,:)
 complex(dpc),allocatable :: chi(:,:,:),chitmp(:,:),chitmps(:,:),eigvec(:,:)
 complex(dpc),allocatable :: eigvec1(:,:),matr(:),mm(:,:,:),mtemp(:,:),rhog(:)
 complex(dpc),allocatable :: rhogg(:,:),tmp1(:),work(:),zhpev1(:),zz2(:,:)

!*************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' cppm4par : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHEGV' :: zhegv
#endif

 qibz => Qmesh%ibz(1:3,1:Qmesh%nibz)
 if (PRESENT(iqiA)) then 
  call assert((nqiA==1),'nqiA should be 1',__FILE__,__LINE__)
  ltest=(iqiA>0.and.iqiA<=Qmesh%nibz)
  call assert(ltest,'iqiA out of range',__FILE__,__LINE__)
  qibz => Qmesh%ibz(1:3,iqia:iqia)
 end if
 !FIXME It seems that q has to be non zero to avoid problems with LAPACK thus we use the first small q
 !do iq=1,nqiA
 ! if (normv(qibz(:,iq),gmet,'G')<GW_TOLQ0)) then 
 !  qibz(:,iq)=Qmesh%small(:,1)
 !  write(*,*)' cppm4par COMMENT : using small q ',qibz(:,iq)
 ! end if
 !end do

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)
 !
 ! === Calculate density in G space rhog(G) ===
 allocate(rhog_dp(2,nfftf),rhog(nfftf),stat=istat)      ; if (istat/=0) stop 'rhog_dp/rhog out of memory'
 allocate(igfft(npwc,npwc),rhogg(npwc,npwc),stat=istat) ; if (istat/=0) stop 'igfft/rhogg out of memory'
 !
 ! Conduct FFT tho(r)-->rhog(G)
 ! FIXME this has to be fixed, rho(G) should be passed instead of doing FFT for each q
 ! Moreover MPI_enreg is local ????? 
 call fourdp(1,rhog_dp,rho,-1,MPI_enreg,nfftf,ngfftf,paral_kgb,0)

 rhog(1:nfftf)=CMPLX(rhog_dp(1,1:nfftf),rhog_dp(2,1:nfftf))
 deallocate(rhog_dp)
 !
 ! Calculate the FFT index of each (G-Gp) vector and assign the value
 ! of the correspondent density simultaneously
 ngfft1=ngfftf(1) ; ngfft2=ngfftf(2) ; ngfft3=ngfftf(3)
 call cggfft(npwc,ngfft1,ngfft2,ngfft3,gvec,igfft)

 do ig=1,npwc
  do igp=1,npwc
   if (igfft(ig,igp)>nfftf) then
    ! well by definition igfft <= nfftf
    write (msg,'(4a)')ch10,&
&    ' cppm3par : BUG- ',ch10,&
&    ' cannot find rho(G-Gpr) for some G, Gpr '
    call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
   end if
   rhogg(ig,igp)=rhog(igfft(ig,igp))
  end do
 end do
 deallocate(igfft,rhog)
 !
 ! Now we have rhogg, calculate the M matrix (q+G1).(q+G2) n(G1-G2)
 allocate(mm(npwc,npwc,nqiA),stat=istat) ; if (istat/=0) stop 'mm out of memory'

 do iq=1,nqiA
  do ig=1,npwc
   gpq(:)=gvec(:,ig)+qibz(:,iq)
   do igp=1,npwc
    gppq(:)=gvec(:,igp)+qibz(:,iq)
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

 !MG TODO too much memory in chi, we can do all this stuff inside a loop
 allocate(chitmp(npwc,npwc),chi(npwc,npwc,nqiA),stat=istat) ; if (istat/=0) stop 'cppm4par out of memory'
 allocate(qplusg(npwc))
 !
 ! Extract the full polarizability from \tilde \epsilon^{-1}
 ! \tilde\epsilon^{-1}_{G1 G2} = \delta_{G1 G2} + 4\pi \frac{\chi_{G1 G2}}{|q+G1| |q+G2|}
 do iq=1,nqiA
  chitmp(:,:)=epsm1(:,:,1,iq)
  call cvc(nqiA,iq,qibz,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q
  do ig=1,npwc
   chitmp(ig,ig)=chitmp(ig,ig)-one
  end do
  do ig=1,npwc
   do igp=1,npwc
    chi(ig,igp,iq)=chitmp(ig,igp)*qplusg(ig)*qplusg(igp)/four_pi
   end do
  end do
 end do
 deallocate(chitmp)

!DEBUG
!do iq=1,nqiA
!allocate(eigval1(npwc),stat=istat)
!if(istat/=0) stop 'eigval1 out of memory'
!allocate(eigvec1(npwc,npwc),stat=istat)
!if(istat/=0) stop 'eigvec1 out of memory'
!allocate(matr(npwc*(npwc+1)/2))
!if(istat/=0) stop 'matr out of memory'
!allocate(zhpev2(3*npwc-2),stat=istat)
!if(istat/=0) stop 'zhpev2 of memory'
!allocate(zhpev1(2*npwc-1),stat=istat)
!if(istat/=0) stop 'zhpev1 of memory' ! woking arrays for lapack
!
!idx=1
!do ii=1,npwc
!do jj=1,ii
!matr(idx)=chi(jj,ii,iq)
!idx=idx+1
!end do
!end do
!
!call zhpev('v','u',npwc,matr,eigval1,eigvec1,npwc,&
!&   zhpev1,zhpev2,ierr)
!
!deallocate(zhpev1,zhpev2,matr,eigval1,eigvec1)
!end do
!ENDDEBUG
 !
 ! === Solve chi*X = Lambda M*X where Lambda=-1/em(q)**2 ===
 do iq=1,nqiA

  allocate(eigval(npwc),eigvec(npwc,npwc),stat=istat)
  allocate(mtemp(npwc,npwc),chitmps(npwc,npwc),stat=istat) ! temp working arrays
  allocate(work(2*npwc-1),rwork(3*npwc-2),stat=istat) 
  if(istat/=0) stop 'rwork out of memory'
  !
  ! Copy chi and mm into working arrays
  chitmps(:,:)=chi(:,:,iq)
  mtemp(:,:)=mm(:,:,iq)
  lowork=2*npwc-1

  call ZHEGV(1,'V','U',npwc,chitmps,npwc,mtemp,npwc,eigval,work,lowork,rwork,ierr)
  
  if (ierr/=0) then 
   write(msg,'(4a,i5)')ch10,&
&   ' cppm4par : ERROR - ',ch10,&
&   ' zhegv reported info = ',ierr
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if 
  !
  ! Eigenvectors are normalized as : X_i^* M X_j = \delta_{ij}    
  eigvec(:,:)=chitmps(:,:)
  deallocate(mtemp,chitmps,work,rwork) 
  !
  ! === Calculate the plasmon pole parameters ===
  allocate(tmp1(npwc),stat=istat) !eigenvectors
  if (istat/=0) stop 'tmp1 out of memory'
  allocate(zz2(npwc,npwc),stat=istat) ! checking
  if (istat/=0) stop 'zz out of memory'
  !
  ! good check:
  ! the lowest plasmon energy on gamma should be
  ! close to experimental plasma energy within an error of 10 %
  ! this error can be reduced further if one includes the non local
  ! commutators in the calculation of the polarizability at q==0
  zz2(:,:)=(0.0,0.0)
  call cvc(nqiA,iq,qibz,npwc,gvec,gprimd,qplusg) !MG TODO here take care of small q

  do ii=1,npwc

   ! keeping in mind that the above matrix is negative definite
   ! we might have a small problem with the eigval that correspond to large G vectors
   ! i.e. DM band index, where the eigevalues become very small with
   ! possibility of being small positive numbers (due to numerical problems)
   ! thus as a caution one can use the following condition
   ! this will not affect the result since such a huge plasmon energy give almost zero
   ! contribution to the self correlation energy

   if (eigval(ii)>=zero) then
    eigval(ii) = -1.0d-4 
    if (eigval(ii)>1.0d-3) then
     eigval(ii) = -1.0d-22
     write(msg,'(4a,i6,a,es16.6)')ch10,&
&     ' cppm4par : WARNING - ' ,ch10,&
&     '  imaginary plasmon pole eigenenergy, eigenvector number ',ii,' with eigval',eigval(ii),ch10
     call wrtout(6,msg,'COLL')
     call leave_new('COLL')
    end if
   end if
   !
   ! === Save plasmon energies ===
   omegatw(ii,1,iq)=SQRT(-1/eigval(ii))
   !
   ! Calculate and save scaled plasmon-pole eigenvectors 
   ! defined as \sqrt{4\pi} \frac{Mx}{\sqrt{\tilde\omega} |q+G|}
   tmp1(:)=eigvec(:,ii)

   do ig=1,npwc
    do igp=1,npwc
     zz2(ig,ii)=zz2(ig,ii)+mm(ig,igp,iq)*tmp1(igp)
    end do
    bigomegatwsq(ig,ii,iq)=SQRT(four_pi)*zz2(ig,ii)/SQRT(omegatw(ii,1,iq))
    bigomegatwsq(ig,ii,iq)=bigomegatwsq(ig,ii,iq)/qplusg(ig)
   end do
  end do
  deallocate(tmp1,eigvec,eigval,zz2)
 end do !iq

 deallocate(qplusg,chi,rhogg,mm)

 bar=REPEAT('-',80)
 write(msg,'(3a)')bar,ch10,&
& ' plasmon energies vs q vector shown for lowest 10 bands                 '
 call wrtout(ab_out,msg,'COLL')
 do iq=1,nqiA
  write(msg,'(2x,i3,5x,10f7.3)')iq,(REAL(omegatw(ig,1,iq))*Ha_eV,ig=1,10)
  call wrtout(ab_out,msg,'COLL')
 end do
 write(msg,'(a)')bar
 call wrtout(ab_out,msg,'COLL')

 write(msg,'(2a,f12.8,2a,3i5)')ch10,&
 ' cppm2par : omega twiddle minval [eV]  = ',MINVAL(ABS(omegatw(:,:,:)))*Ha_eV,ch10,&
 '            omega twiddle min location = ',MINLOC(ABS(omegatw(:,:,:)))
 call wrtout(std_out,msg,'COLL')

#if defined DEBUG_MODE
 write(msg,'(a)')' cppm4par : exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine cppm4par
!!***
