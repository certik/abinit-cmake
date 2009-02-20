!{\src2tex{textfont=tt}}
!!****f* ABINIT/zerosym
!! NAME
!! zerosym
!!
!! FUNCTION
!! Symmetrize an array on the FFT grid by vanishing some term on the boundaries.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GZ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex= if 1, input array is REAL, if 2, input array is COMPLEX
!!  n1,n2,n3=FFT dimensions nfft=n1*n2*n3
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  array(cplex,n1*n2*n3)=complex array to be symetrized
!!
!! PARENTS
!!      atm2fft,atm2fft3,dyfro3,forces,hartre,initro,pawmknhat,pawmknhat3
!!      stress,transgrid
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine zerosym(array,cplex,mpi_enreg,n1,n2,n3)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,n1,n2,n3
 type(mpi_type) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: array(cplex,n1*n2*n3)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,i3,ifft,ifft_proc,index,j,j1,j2,j3,me_fft,nd2,nn12,r2

! **********************************************************************
 me_fft=mpi_enreg%me_fft
 nd2=(n2-1)/mpi_enreg%nproc_fft+1
!DEBUG
!write(6,*)'enter zerosym',mpi_enreg%nproc_fft,mpi_enreg%me_fft
!ENDEBUG
 nn12=n1*n2
 if (mod(n1,2)==0) then
  index=n1/2+1-nn12-n1
  do i3=1,n3
   index=index+nn12;ifft=index
   do i2=1,n2
    ifft=ifft+n1
    if(mpi_enreg%paral_compil_fft==1) then
! MPIWF: consider ifft only if it is treated by the current proc and compute its adress
     j=ifft-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);r2=modulo(j2,nd2)
     if(modulo(j/n1,n2)/nd2==me_fft) then ! MPIWF this ifft is to be treated by me_fft
      ifft_proc=n1*(nd2*j3+r2)+j1+1 !this is ifft in the current proc
      array(:,ifft_proc)=zero
     end if
    else
     array(:,ifft)=zero
    end if
   end do
  end do
 end if

 if (mod(n2,2)==0) then
  index=n1*(n2/2+1)-nn12-n1
  do i3=1,n3
   index=index+nn12;ifft=index
   do i1=1,n1
    ifft=ifft+1
    if(mpi_enreg%paral_compil_fft==1) then
! MPIWF: consider ifft only if it is treated by the current proc and compute its adress
     j=ifft-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);r2=modulo(j2,nd2)
     if(modulo(j/n1,n2)/nd2==me_fft) then ! MPIWF this ifft is to be treated by me_fft
      ifft_proc=n1*(nd2*j3+r2)+j1+1 !this is ifft in the current proc
      array(:,ifft_proc)=zero
     end if
    else
     array(:,ifft)=zero
    end if
   end do
  end do
 end if

 if (mod(n3,2)==0) then
  index=nn12*(n3/2+1)-nn12-n1
  do i2=1,n2
   index=index+n1;ifft=index
   do i1=1,n1
    ifft=ifft+1
    if(mpi_enreg%paral_compil_fft==1) then
! MPIWF: consider ifft only if it is treated by the current proc and compute its adress
     j=ifft-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);r2=modulo(j2,nd2)
     if(modulo(j/n1,n2)/nd2==me_fft) then ! MPIWF this ifft is to be treated by me_fft
      ifft_proc=n1*(nd2*j3+r2)+j1+1 !this is ifft in the current proc
      array(:,ifft_proc)=zero
     end if
    else
     array(:,ifft)=zero
    end if
   end do
  end do
 end if

end subroutine zerosym
!!***
