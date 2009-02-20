!{\src2tex{textfont=tt}}
!!****f* ABINIT/kgindex
!! NAME
!! kgindex
!!
!! FUNCTION
!! Compute the index of each plane wave on a FFT grid.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  kg_k(3,npw_k)=dimensionless coords of G vecs (integer)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  npw_k=number of planewaves
!!
!! OUTPUT
!!  indpw_k(npw_k)=linear list number (in fft box) of given G vector for the current processor (local adress)
!                 =0 if kg_k(ipw) is not treated by this procesor
!!  mask(npw_k)=True if  kg_k(ipw) belongs to this procesor, false otherwise.
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      prcref
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine kgindex(indpw_k,kg_k,mask,mpi_enreg,ngfft,npw_k)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npw_k
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k(3,npw_k),ngfft(18)
 integer,intent(out) :: indpw_k(npw_k)
 logical,intent(out) :: mask(npw_k)
!Local variables-------------------------------
!scalars
 integer :: ig,ig1,ig2,ig3,me_fft,n1,n2,n3,nd2,nproc_fft

! *************************************************************************
!mpi_enreg is not necessary in this case (the info is also in ngfft), but much more easy to read...
!DEBUG
!write(6,*)' kgindex : enter '
!ENDDEBUG

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Use the following indexing (N means ngfft of the adequate direction)
!0 1 2 3 ... N/2    -(N-1)/2 ... -1    <= kg
!1 2 3 4 ....N/2+1  N/2+2    ...  N    <= index
 me_fft=mpi_enreg%me_fft
 nd2=(n2-1)/mpi_enreg%nproc_fft+1
 do ig=1,npw_k
  ig1=modulo(kg_k(1,ig),n1)
  ig2=modulo(kg_k(2,ig),n2)
  ig3=modulo(kg_k(3,ig),n3)
  if(ig2/nd2==me_fft) then
   ig2=modulo(ig2,nd2)
   indpw_k(ig)=ig1+1+n1*(ig2+nd2*ig3)
   mask(ig)=.true.
  else
   indpw_k(ig)=0
   mask(ig)=.false.
  endif
 end do

!DEBUG
!write(6,*)' kgindex : exit'
!stop
!ENDDEBUG

end subroutine kgindex
!!***
