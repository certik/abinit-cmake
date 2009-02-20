!{\src2tex{textfont=tt}}
!!****f* ABINIT/indirect_parallel_Fourier
!! NAME
!! indirect_parallel_Fourier
!!
!! FUNCTION
!! The purpose of this routine is to transfer data from right to left
!! right(:,index(i))=left(:,i)
!! The difficulty is that right and left are distributed among processors
!! We will suppose that the distribution is done as a density in Fourier space
!! We first order the right hand side data according to the processor
!! in which they are going to be located in the left hand side.
!! This is done is a way such that  a mpi_alltoall put the data on the correct processor.
!! We also transfer their future adress.
!! A final ordering put everything in place
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  index(sizeindex)= global adress for the transfer from right to left
!!  left(2,nleft)=left hand side
!!  mpi_enreg=informations about MPI parallelization
!!  ngleft(18)=contain all needed information about 3D FFT for the left hand side
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  ngright(18)=contain all needed information about 3D FFT for the right hand side
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nleft=second dimension of left array (for this processor)
!!  nright=second dimension of right array (for this processor)
!!  sizeindex=size of the index array (different form nright, because it is global to all proccessors)
!! OUTPUT
!!  left(2,nleft)=the elements of the right hand side, at the correct palce in the correct processor
!!
!! NOTES
!!  A lot of things to improve.
!! PARENTS
!!      transgrid
!!
!! CHILDREN
!!      mpi_alltoall, xcomm_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine indirect_parallel_Fourier(index,left,mpi_enreg,ngleft,ngright,nleft,nright,paral_kgb,right,sizeindex)
 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ngleft(18),ngright(18),nleft,nright,paral_kgb,sizeindex
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: index(sizeindex)
 real(dp),intent(in) :: right(2,nright)
 real(dp),intent(inout) :: left(2,nleft)

!Local variables ---------------------------------------
!scalars
 integer :: ierr,i_global,ileft,ileft_global,iright,iright_global,&
 &j,j1,j2,j3,j_global,jleft_global,&
 &jleft_local,me_fft,n1l,n2l,n3l,n1r,n2r,n3r,nd2l,nd2r,&
 &nproc_fft,old_paral_level,proc_dest,r2,siz_slice_max,spaceComm
 character(len=500) :: message
!arrays
 integer,allocatable :: index_recv(:),index_send(:),siz_slice(:)
 real(dp),allocatable :: right_send(:,:),right_recv(:,:)

! *************************************************************************
 n1r=ngright(1);n2r=ngright(2);n3r=ngright(3)
 n1l=ngleft(1) ;n2l=ngleft(2) ;n3l=ngleft(3)
 nproc_fft=mpi_enreg%nproc_fft; me_fft=mpi_enreg%me_fft
 nd2r=n2r/nproc_fft; nd2l=n2l/nproc_fft
 allocate(siz_slice(nproc_fft))
 siz_slice(:)=0
 do i_global=1,sizeindex !look for the maximal size of slice of data
  j_global=index(i_global)!; write(6,*) j_global,i_global
  if(j_global /=0) then
   proc_dest=modulo((j_global-1)/n1l,n2l)/nd2l
   siz_slice(proc_dest+1)=siz_slice(proc_dest+1)+1
!DEBUG
!write(6,*) 'in indirect proc',proc_dest,siz_slice(proc_dest+1)
!ENDDEBUG
  end if
 end do
 siz_slice_max=maxval(siz_slice) !This value could be made smaller by looking locally
!and performing a allgather with a max
!DEBUG
!write(6,*) 'siz_slice,sizeindex,siz_slice',siz_slice(:),sizeindex,siz_slice_max
!write(6,*) 'sizeindex,nright,nleft',sizeindex,nright,nleft
!ENDDEBUG
 allocate(right_send(2,nproc_fft*siz_slice_max),index_send(nproc_fft*siz_slice_max))
 siz_slice(:)=0; index_send(:)=0; right_send(:,:)=zero
 do iright=1,nright
  j=iright-1;j1=modulo(j,n1r);j2=modulo(j/n1r,nd2r);j3=j/(n1r*nd2r)
  j2=j2+me_fft*nd2r
  iright_global=n1r*(n2r*j3+j2)+j1+1
  jleft_global=index(iright_global)
  if(jleft_global/=0)then
   j=jleft_global-1;j1=modulo(j,n1l);j2=modulo(j/n1l,n2l);j3=j/(n1l*n2l);r2=modulo(j2,nd2l)
   jleft_local=n1l*(nd2l*j3+r2)+j1+1
   proc_dest=j2/nd2l
   siz_slice(proc_dest+1)=siz_slice(proc_dest+1)+1
   right_send(:,proc_dest*siz_slice_max+siz_slice(proc_dest+1))=right(:,iright)
   index_send(proc_dest*siz_slice_max+siz_slice(proc_dest+1))=jleft_local
!DEBUG
!   write(6,*) 'loop ir',jleft_local,jleft_global,iright_global,iright
!ENDDEBUG
  end if
 end do
 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm)
 if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
 allocate(right_recv(2,nproc_fft*siz_slice_max),index_recv(nproc_fft*siz_slice_max))
#if defined MPI
  if(paral_kgb == 1) then
    call mpi_alltoall (right_send,2*siz_slice_max, &
                          MPI_double_precision, &
                          right_recv,2*siz_slice_max, &
                          MPI_double_precision,spaceComm,ierr)
    call mpi_alltoall (index_send,siz_slice_max, &
                          MPI_integer, &
                          index_recv,siz_slice_max, &
                          MPI_integer,spaceComm,ierr)
  endif
#endif
 mpi_enreg%paral_level=old_paral_level
 do ileft=1,siz_slice_max*nproc_fft
!DEBUG
!write(6,*)index_recv(ileft)
!ENDEBUG
 if(index_recv(ileft) /=0 ) left(:,index_recv(ileft))=right_recv(:,ileft)
 end do
 deallocate(right_recv,index_recv,right_send,index_send,siz_slice)

end subroutine indirect_parallel_Fourier
!!***
