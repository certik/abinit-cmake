!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprod_g
!! NAME
!! dotprod_g
!!
!!
!! FUNCTION
!! Compute dot product of complex vectors vect1 and vect2 (can be the same)
!! Take into account the storage mode of the vectors (istwf_k)
!! If option=1, compute only real part, if option=2 compute also imaginary
!! part. If the number of calls to the dot product scales quadratically
!! with the volume of the system, it is preferable not to
!! call the present routine, but but to write a specially
!! optimized routine, that will avoid many branches related to
!! the existence of 'option' and 'istwf_k'.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  istwf_k=option parameter that describes the storage of wfs
!!  vect1(2,npw)=first vector (one should take its complex conjugate)
!!  vect2(2,npw)=second vector
!!  mpi_enreg=informations about MPI parallelization
!!  npw= (effective) number of planewaves at this k point.
!!  option= 1 if only real part to be computed,
!!          2 if both real and imaginary.
!!          3 if in case istwf_k==1 must compute real and imaginary parts,
!!               but if  istwf_k==2 must compute only real part
!!
!! OUTPUT
!!  $doti=\Im ( <vect1|vect2> )$ , output only if option=2 and eventually option=3.
!!  $dotr=\Re ( <vect1|vect2> )$
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      cgwf,cgwf3,eig1fixed,eig2tot,mkresi,nstwf3,nstwf4,recip_ylm,vtowfk3
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,timab,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dotprod_g(dotr,doti,istwf_k,mpi_enreg,npw,option,vect1,vect2)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_11contract
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istwf_k,npw,option
 real(dp),intent(out) :: doti,dotr
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: vect1(2,npw),vect2(2,npw)

!Local variables-------------------------------
!scalars
 integer :: ierr,ipw,old_paral_level,spaceComm
!arrays
 real(dp) :: dotarr(2),tsec(2)
!no_abirules
#if defined CONTRACT
 integer :: ii
 character(len=9) :: subrnm
#endif

! *************************************************************************

#if defined CONTRACT
 subrnm='dotprod_g'
 call contract_int_list(subrnm,'istwf_k',istwf_k,(/ (ii,ii=1,9) /),9)
 call contract_int_ge_val(subrnm,'npw',npw,1)
 call contract_int_list(subrnm,'option',option,(/1,2,3/),3)
#endif

 if(istwf_k==1)then
! General k-point
  if(option==1)then
   dotr=0.0d0
!  $OMP PARALLEL DO ORDERED PRIVATE(ipw) REDUCTION(+:dotr) &
!  $OMP&SHARED(vect1,vect2,npw)
   do ipw=1,npw
    dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
   end do
!  $OMP END PARALLEL DO
  else
   dotr=0.0d0 ; doti=0.0d0
!  $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:doti,dotr) &
!  $OMP&SHARED(vect1,vect2,npw)
   do ipw=1,npw
    dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
    doti=doti-vect1(2,ipw)*vect2(1,ipw)+vect1(1,ipw)*vect2(2,ipw)
   end do
!  $OMP END PARALLEL DO
  end if
 else
! Gamma k-point
  if(istwf_k==2 .and. mpi_enreg%me_g0==1)then
   if(option/=2)then
    dotr=0.5d0*vect1(1,1)*vect2(1,1)
!   $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:dotr) &
!   $OMP&SHARED(vect1,vect2,npw)
    do ipw=2,npw
     dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
    end do
!   $OMP END PARALLEL DO
    dotr=2.0d0*dotr
   else ! option==2
    dotr=0.5d0*( vect1(1,1)*vect2(1,1)+vect1(2,1)*vect2(2,1))
    doti=0.5d0*(-vect1(2,1)*vect2(1,1)+vect1(1,1)*vect2(2,1))
!   $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:doti,dotr) &
!   $OMP&SHARED(vect1,vect2,npw)
    do ipw=2,npw
     dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
     doti=doti-vect1(2,ipw)*vect2(1,ipw)+vect1(1,ipw)*vect2(2,ipw)
    end do
!   $OMP END PARALLEL DO
    dotr=2.0d0*dotr
    doti=2.0d0*doti
   end if
  else
!  Other TR k-points
   if(option/=2)then
    dotr=0.0d0
!   $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:dotr) &
!   $OMP&SHARED(vect1,vect2,npw)
    do ipw=1,npw
     dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
    end do
!   $OMP END PARALLEL DO
    dotr=2.0d0*dotr
   else
    dotr=0.0d0 ; doti=0.0d0
!   $OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:doti,dotr) &
!   $OMP&SHARED(vect1,vect2,npw)
    do ipw=1,npw
     dotr=dotr+vect1(1,ipw)*vect2(1,ipw)+vect1(2,ipw)*vect2(2,ipw)
     doti=doti-vect1(2,ipw)*vect2(1,ipw)+vect1(1,ipw)*vect2(2,ipw)
    end do
!   $OMP END PARALLEL DO
    dotr=2.0d0*dotr
    doti=2.0d0*doti
   end if
  end if
 end if

!NEW CODING
 if(mpi_enreg%paral_compil_fft==1)then
  old_paral_level=mpi_enreg%paral_level
  mpi_enreg%paral_level=3
  call xcomm_init(mpi_enreg,spaceComm)
  call timab(48,1,tsec)
  dotarr(1)=dotr ; dotarr(2)=doti
  call xsum_mpi(dotarr,spaceComm,ierr)
  dotr=dotarr(1) ; doti=dotarr(2)
  call timab(48,2,tsec)
  mpi_enreg%paral_level=old_paral_level
! It might be needed to set mpi_enreg%paral_level back to its previous value ? !
 end if
!END OF NEW CODING

end subroutine dotprod_g
!!***
