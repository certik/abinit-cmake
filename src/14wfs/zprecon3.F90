!{\src2tex{textfont=tt}}
!!****f* ABINIT/zprecon3
!!
!! NAME
!! zprecon3
!!
!! FUNCTION
!! precondition $<g|(h-e_{n,k})|c_{n,k}>$
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (dca, xg, gmr)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  $cg(2,npw)=<g|c_{n,k}>$.
!!  $eval=current band eigenvalue=<c_{n,k}|h|c_{n,k}>$.
!!  istwf_k=option parameter that describes the storage of wfs
!!  kinpw(npw)=(modified) kinetic energy for each plane wave (hartree)
!!  nspinor=number of spinorial components of the wavefunctions
!!  $vect(2,npw)=<g|h|c_{n,k}>$.
!!  npw=number of planewaves at this k point.
!!  optekin= 1 if the kinetic energy used in preconditionning is modified
!!             according to Kresse, Furthmuller, PRB 54, 11169 (1996)
!!           0 otherwise
!!  vectsize= size of vectors
!!
!! OUTPUT
!!  $vect(2,npw)=<g|(h-eval)|c_{n,k}>*(polynomial ratio)$
!!
!! PARENTS
!!      lobpcgccIIwf,lobpcgccwf
!!
!! CHILDREN
!!      timab,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine zprecon3(cg,eval,blocksize,istwf_k,kinpw,mpi_enreg,npw,nspinor,optekin,ghc,vect,vectsize)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,istwf_k,npw,nspinor,optekin,vectsize
 type(mpi_type) :: mpi_enreg
!arrays
 real(dp) :: kinpw(npw)
 complex(dpc) :: cg(vectsize,blocksize),eval(blocksize,blocksize)
 complex(dpc) :: ghc(vectsize,blocksize),vect(vectsize,blocksize)

!Local variables-------------------------------
!scalars
 integer :: iblocksize,ierr,ig,igs,ipw1,ispinor,old_paral_level,spaceComm
 real(dp) :: fac,poly,xx
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: ek0(:),ek0_inv(:)

! *************************************************************************

 call timab(536,1,tsec)
!compute mean kinetic energy of all bands
 allocate(ek0(blocksize),ek0_inv(blocksize))
 do iblocksize=1,blocksize
  if(istwf_k==1)then
   ek0(iblocksize)=0.0_dp
   do ispinor=1,nspinor
    igs=(ispinor-1)*npw
!   $omp parallel do private(ig) reduction(+:ek0) &
!   $omp&shared(cg,igs,kinpw,npw)
    do ig=1+igs,npw+igs
     if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
      ek0(iblocksize)=ek0(iblocksize)+kinpw(ig-igs)*&
&      (real(cg(ig,iblocksize))**2+aimag(cg(ig,iblocksize))**2)
     end if
    end do
!   $omp end parallel do
   end do
  else if (istwf_k>=2)then
   if (istwf_k==2 .and. mpi_enreg%me_g0 == 1)then
    ek0(iblocksize)=0.0_dp ; ipw1=2
    if(kinpw(1)<huge(0.0_dp)*1.d-11)ek0(iblocksize)=0.5_dp*kinpw(1)*cg(1,iblocksize)**2
   else
    ek0(iblocksize)=0.0_dp ; ipw1=1
   end if
!  $omp parallel do private(ig) reduction(+:ek0) &
!  $omp&shared(cg,ipw1,kinpw,npw)
   do ig=ipw1,npw
    if(kinpw(ig)<huge(0.0_dp)*1.d-11)then
     ek0(iblocksize)=ek0(iblocksize)+&
&     kinpw(ig)*(real(cg(ig,iblocksize))**2+real(cg(ig+npw-1,iblocksize))**2)
    end if
   end do
!  $omp end parallel do
   ek0=two*ek0
  end if
 end do

 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm)
 if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
 call xsum_mpi(ek0,spaceComm,ierr)
 mpi_enreg%paral_level= old_paral_level

 do iblocksize=1,blocksize
  if(ek0(iblocksize)<1.0d-10)then
   write(message, '(a,a,a,a,a,a)' )ch10,&
&   ' precon : warning -',ch10,&
&   '  the mean kinetic energy of a wavefunction vanishes.',ch10,&
&   '  it is reset to 0.1ha.'
   call wrtout(6,message,'pers')
   ek0(iblocksize)=0.1_dp
  end if
 end do
 if (optekin==1) then
  ek0_inv(:)=2.0_dp/(3._dp*ek0(:))
 else
  ek0_inv(:)=1.0_dp/ek0(:)
 end if
!
!carry out preconditioning
 do iblocksize=1,blocksize
  do ispinor=1,nspinor
   igs=(ispinor-1)*npw
!  $omp parallel do private(fac,ig,poly,xx) &
!  $omp&shared(cg,ek0_inv,eval,kinpw,igs,npw,vect)
   do ig=1+igs,npw+igs
    if(kinpw(ig-igs)<huge(0.0_dp)*1.d-11)then
     xx=kinpw(ig-igs)*ek0_inv(iblocksize)
!    teter polynomial ratio
     poly=27._dp+xx*(18._dp+xx*(12._dp+xx*8._dp))
     fac=poly/(poly+16._dp*xx**4)
     if (optekin==1) fac=two*fac
     vect(ig,iblocksize)=(ghc(ig,iblocksize)-&
&     eval(iblocksize,iblocksize)*cg(ig,iblocksize) )*fac
    else
     vect(ig,iblocksize)=dcmplx(0.0_dp,0.0_dp)
    end if
   end do
!  $omp end parallel do
  end do
 end do
 deallocate(ek0,ek0_inv)

 call timab(536,2,tsec)

end subroutine zprecon3
!!***
