!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprodm_v
!! NAME
!! dotprodm_v
!!
!!
!! FUNCTION
!! For two sets of potentials,
!! compute dot product of each pair of two potentials (integral over FFT grid), to obtain
!! a series of square residual-like quantity (so the sum of product of values
!! is NOT divided by the number of FFT points, and NOT multiplied by the primitive cell volume).
!! Take into account the spin components of the potentials (nspden),
!! and sum over them.
!! Need the index of the first pair of potentials to be treated, in each array
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
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  cpldot=if 1, the dot array is real, if 2, the dot array is complex
!!  index1=index of the first potential to be treated in the potarr1 array
!!  index2=index of the first potential to be treated in the potarr2 array
!!  mpi_enreg=informations about MPI parallelization
!!  mult1=number of potentials to be treated in the first set
!!  mult2=number of potentials to be treated in the second set
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  npot1= third dimension of the potarr1 array
!!  npot2= third dimension of the potarr2 array
!!  nspden=number of spin-density components
!!  potarr1(cplex*nfft,nspden,npot)=first array of real space potentials on FFT grid
!!    (if cplex=2 and cpldot=2, potarr1 is the array that will be complex conjugated)
!!  potarr2(cplex*nfft,nspden,npot)=second array of real space potentials on FFT grid
!!
!! OUTPUT
!!  dot(cpldot,mult1,mult2)= series of values of the dot product
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      scfopt
!!
!! CHILDREN
!!      contract_int_ge_val,contract_int_list,leave_new,timab,wrtout,xcomm_init
!!      xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dotprodm_v(cplex,cpldot,dot,index1,index2,mpi_enreg,mult1,mult2,nfft,npot1,npot2,nspden,potarr1,potarr2)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11contract
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cpldot,cplex,index1,index2,mult1,mult2,nfft,npot1,npot2
 integer,intent(in) :: nspden
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: potarr1(cplex*nfft,nspden,npot1)
 real(dp),intent(in) :: potarr2(cplex*nfft,nspden,npot2)
 real(dp),intent(out) :: dot(cpldot,mult1,mult2)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ierr,ifft,ispden,old_paral_level,spaceComm
 real(dp) :: ai,ar
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
!no_abirules
#if defined CONTRACT
 character(len=10) :: subrnm
#endif

! *************************************************************************

#if defined CONTRACT
 subrnm='dotprodm_v'
!Real or complex inputs are coded
 call contract_int_list(subrnm,'cplex',cplex,(/1,2/),2)
!Real or complex outputs are coded
 call contract_int_list(subrnm,'cpldot',cpldot,(/1,2/),2)
 call contract_int_ge_val(subrnm,'index1',index1,1)
 call contract_int_ge_val(subrnm,'index2',index2,1)
 call contract_int_ge_val(subrnm,'mult1',mult1,1)
 call contract_int_ge_val(subrnm,'mult2',mult2,1)
 call contract_int_ge_val(subrnm,'nfft',nfft,1)
 call contract_int_ge_val(subrnm,'npot1',npot1,1)
 call contract_int_ge_val(subrnm,'npot2',npot2,1)
 call contract_int_list(subrnm,'nspden',nspden,(/1,2,4/),3)
 call contract_int_ge_val(subrnm,'npot1-index1-mult1',npot1-index1-mult1,-1)
 call contract_int_ge_val(subrnm,'npot2-index2-mult2',npot2-index2-mult2,-1)
#endif

 if (cplex==2.and.cpldot==2.and.nspden==4) then
  write(message, '(a,a,a,a,i3,a,a,a,a)' ) ch10,&
&  ' dotprodm_v : ERROR - ',ch10,&
&  '  cplex=2, cpldot=2  and nspden=4 not compatible !'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 if(cplex==1 .or. cpldot==1)then

  do i1=1,mult1
   do i2=1,mult2
    ar=zero
    do ispden=1,min(nspden,2)
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar)
     do ifft=1,cplex*nfft
      ar=ar + potarr1(ifft,ispden,index1+i1-1)*potarr2(ifft,ispden,index2+i2-1)
     end do
!    $OMP END PARALLEL DO
    end do
    dot(1,i1,i2)=ar
    if (nspden==4) then
     ar=zero
     do ispden=3,4
!     $OMP PARALLEL DO PRIVATE(ifft) &
!     $OMP&SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar)
      do ifft=1,cplex*nfft
       ar=ar + potarr1(ifft,ispden,index1+i1-1)*potarr2(ifft,ispden,index2+i2-1)
      end do
!     $OMP END PARALLEL DO
     end do
     dot(1,i1,i2)=dot(1,i1,i2)+two*ar
    end if
   end do
  end do

 else ! if (cplex==2 .and. cpldot==2)

  do i1=1,mult1
   do i2=1,mult2
    ar=zero ; ai=zero
    do ispden=1,nspden
!    $OMP PARALLEL DO PRIVATE(ifft) &
!    $OMP&SHARED(cplex,i1,i2,index1,index2,ispden,nfft,potarr1,potarr2) REDUCTION(+:ar,ai)
     do ifft=1,nfft
      ar=ar + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1) &
&      + potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1)
      ai=ai + potarr1(2*ifft-1,ispden,index1+i1-1)*potarr2(2*ifft  ,ispden,index2+i2-1) &
&      - potarr1(2*ifft  ,ispden,index1+i1-1)*potarr2(2*ifft-1,ispden,index2+i2-1)
     end do
!    $OMP END PARALLEL DO
    end do
    dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
   end do
  end do
 end if

!XG030513 : MPIWF reduction (addition) on dot is needed here
!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
  old_paral_level=mpi_enreg%paral_level
  mpi_enreg%paral_level=3
  call xcomm_init(mpi_enreg,spaceComm)
  if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
  call timab(48,1,tsec)
  call xsum_mpi(dot,spaceComm ,ierr)
  call timab(48,2,tsec)
  mpi_enreg%paral_level=old_paral_level
 end if
 if(cpldot==2 .and. cplex==1)dot(2,:,:)=zero

end subroutine dotprodm_v
!!***
