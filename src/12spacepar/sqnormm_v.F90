!{\src2tex{textfont=tt}}
!!****f* ABINIT/sqnormm_v
!! NAME
!! sqnormm_v
!!
!!
!! FUNCTION
!! For a series of potentials,
!! compute square of the norm (integral over FFT grid), to obtain
!! a square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden),
!! and sum over them.
!! Need the index of the first potential to be treated, in the provided array
!! of potentials, and the number of potentials to be treated.
!! Might be used to compute just one square of norm, in
!! a big array, such as to avoid copying a potential from a big array
!! to a temporary place.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space function on FFT grid is REAL, if 2, COMPLEX
!!  index=index of the first potential to be treated
!!  mpi_enreg=informations about MPI parallelization
!!  mult=number of potentials to be treated
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  npot= third dimension of the potarr array
!!  nspden=number of spin-density components
!!  potarr(cplex*nfft,nspden,npot)=array of real space potentials on FFT grid
!!
!! OUTPUT
!!  norm2(mult)= value of the square of the norm of the different potentials
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      scfcge,scfopt
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,timab,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine sqnormm_v(cplex,index,mpi_enreg,mult,nfft,norm2,npot,nspden,potarr)

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
 integer,intent(in) :: cplex,index,mult,nfft,npot,nspden
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 real(dp),intent(out) :: norm2(mult)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,ii,ispden,old_paral_level,spaceComm
 real(dp) :: ar
!arrays
 real(dp) :: tsec(2)
!no_abirules
#if defined CONTRACT
 integer :: imult
 character(len=9) :: subrnm
#endif

! *************************************************************************

#if defined CONTRACT
 subrnm='sqnormm_v'
!Real or complex inputs are coded
 call contract_int_list(subrnm,'cplex',cplex,(/1,2/),2)
 call contract_int_ge_val(subrnm,'index',index,1)
 call contract_int_ge_val(subrnm,'mult',mult,1)
 call contract_int_ge_val(subrnm,'nfft',nfft,1)
 call contract_int_ge_val(subrnm,'npot',npot,1)
 call contract_int_list(subrnm,'nspden',nspden,(/1,2,4/),3)
 call contract_int_ge_val(subrnm,'npot-index-mult',npot-index-mult,-1)
#endif

 do ii=1,mult
  ar=zero
  do ispden=1,min(nspden,2)
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(cplex,ii,index,ispden,nfft,potarr) REDUCTION(+:ar)
   do ifft=1,cplex*nfft
    ar=ar + potarr(ifft,ispden,index+ii-1)**2
   end do
!  $OMP END PARALLEL DO
  end do
  norm2(ii)=ar
  if (nspden==4) then
   ar=zero
   do ispden=3,4
!   $OMP PARALLEL DO PRIVATE(ifft) &
!   $OMP&SHARED(cplex,ii,index,ispden,nfft,potarr) REDUCTION(+:ar)
    do ifft=1,cplex*nfft
     ar=ar + potarr(ifft,ispden,index+ii-1)**2
    end do
!   $OMP END PARALLEL DO
   end do
   norm2(ii)=norm2(ii)+two*ar
  end if
 end do


!XG030513 : MPIWF reduction (addition) on norm2 is needed here
!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
  old_paral_level=mpi_enreg%paral_level
  mpi_enreg%paral_level=3
  call xcomm_init(mpi_enreg,spaceComm)
  if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
  call timab(48,1,tsec)
  call xsum_mpi(norm2,spaceComm ,ierr)
  call timab(48,2,tsec)
  mpi_enreg%paral_level=old_paral_level
 end if

end subroutine sqnormm_v
!!***
