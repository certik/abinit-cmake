!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_grid
!! NAME
!! initmpi_grid
!!
!! FUNCTION
!! Initialize the mpi informations for the grid
!!      2D if parallization FFT/BAND (!MPI paral_kgb)
!!      3D if parallization KPT/FFT/BAND (paral_kgb & MPI)
!!
!! COPYRIGHT
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!
!! OUTPUT
!!  mpi_enreg=informations about MPI parallelization
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initmpi_grid(dtset,mpi_enreg)

 use defs_basis
 use defs_datatypes
 use defs_infos

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => initmpi_grid
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!ARGUMENTS 
!=========================================
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg
 
!LOCAL VARIABLES 
!=========================================
 !Variables introduced for MPI version
 integer :: ierr

 !Variables introduced for the bandFFT version
 logical :: reorder
 logical, allocatable :: periode(:), keepdim(:)
 integer, allocatable :: coords(:)
 integer :: np_test

!DEBUG
! write(6,*)' initmpi_grid : enter'
!ENDDEBUG

!Default parameters, for safety of sequential use
 mpi_enreg%me_fft=0
 mpi_enreg%me_kpt=0
 mpi_enreg%me_band=0

#if defined MPI 
!  if(dtset%paral_kgb ==1) then

! TEST THE PARAMETERS OF THE 3D GRID
! ========================================

  if(mpi_enreg%nproc_fft*  &
&    mpi_enreg%nproc_band* &
&    mpi_enreg%nproc_kpt /= mpi_enreg%nproc)then

   write(6,'(8a,i5,a,i5,a,i5,a,i5)') ch10,&
&   ' initmpi_grid : BUG -',ch10,&
&   '  The number of band*FFT*kpt processors, npband*npfft*kpt, should be',ch10,&
&   '  equal to the total number of processors, nproc.',ch10,&
&   '  However, npband=',mpi_enreg%nproc_band,&
&             ' npfft =',mpi_enreg%nproc_fft,&
&             ' npkpt =',mpi_enreg%nproc_kpt,&
&   ' and nproc=',mpi_enreg%nproc
   call leave_new('PERS')

  end if

  write(6,*) 'npfft, npband and npkpt',&
              & mpi_enreg%nproc_fft,&
              & mpi_enreg%nproc_band,&
              & mpi_enreg%nproc_kpt

! CREATE THE 3D GRID
! ==================================================

  mpi_enreg%dimcart=3

  allocate(mpi_enreg%sizecart(mpi_enreg%dimcart))
  allocate(periode           (mpi_enreg%dimcart))
  allocate(mpi_enreg%coords  (mpi_enreg%dimcart))

  mpi_enreg%sizecart(1) = mpi_enreg%nproc_fft
  mpi_enreg%sizecart(2) = mpi_enreg%nproc_band
  mpi_enreg%sizecart(3) = mpi_enreg%nproc_kpt

  periode(:)=.false.
  reorder   =.false.

! create the cartesian grid with commcart as a communicator.
  call MPI_CART_CREATE(MPI_COMM_WORLD,mpi_enreg%dimcart,mpi_enreg%sizecart,periode,&
&  reorder,mpi_enreg%commcart_3d,ierr)


! Find the index and coordinates of the current  processor
  call MPI_COMM_RANK(mpi_enreg%commcart_3d, mpi_enreg%me_cart, ierr)
  call MPI_CART_COORDS(mpi_enreg%commcart_3d, mpi_enreg%me_cart,  mpi_enreg%dimcart, &
&  mpi_enreg%coords, ierr)


! Create the communicator for space (Fourier) distribution
  allocate(keepdim(mpi_enreg%dimcart))

  keepdim(1)=.true.
  keepdim(2)=.false.
  keepdim(3)=.false.
  call MPI_CART_SUB(mpi_enreg%commcart_3d, keepdim, mpi_enreg%comm_fft,ierr)

! Create the communicator for band distribution
  keepdim(1)=.false.
  keepdim(2)=.true.
  keepdim(3)=.false.
  call MPI_CART_SUB(mpi_enreg%commcart_3d, keepdim, mpi_enreg%comm_band,ierr)


! Create the communicator for kpt distribution
  keepdim(1)=.false.
  keepdim(2)=.false.
  keepdim(3)=.true.
  call MPI_CART_SUB(mpi_enreg%commcart_3d, keepdim, mpi_enreg%comm_kpt,ierr)

  call MPI_COMM_SIZE(mpi_enreg%comm_fft,np_test, ierr)
  write(6,*) 'mpi_enreg%sizecart(1),np_fft' ,mpi_enreg%sizecart(1), np_test
  call MPI_COMM_SIZE(mpi_enreg%comm_band,np_test, ierr)
  write(6,*) 'mpi_enreg%sizecart(2),np_band',mpi_enreg%sizecart(2), np_test
  call MPI_COMM_SIZE(mpi_enreg%comm_kpt,np_test, ierr)
  write(6,*) 'mpi_enreg%sizecart(3),np_kpt',mpi_enreg%sizecart(3), np_test


! Create the communicator for FFT/band distribution
  keepdim(1)=.true.
  keepdim(2)=.true.
  keepdim(3)=.false.
  call MPI_CART_SUB(mpi_enreg%commcart_3d, keepdim, mpi_enreg%commcart,ierr)

  call MPI_COMM_RANK(mpi_enreg%commcart, mpi_enreg%me_cart_2d, ierr)

! Define the correspondance with the fft
  mpi_enreg%me_fft  = mpi_enreg%coords(1)
  mpi_enreg%me_band = mpi_enreg%coords(2)
  mpi_enreg%me_kpt  = mpi_enreg%coords(3)


! WRITE THE PROCESSORS COORDINATES IN THE 3D GRID
! =====================================================
  write(6,*) 'in initmpi_grid : me_fft, me_band, me_kpt are',&
&            mpi_enreg%me_fft,mpi_enreg%me_band,mpi_enreg%me_kpt

  deallocate(periode)

! endif
#endif

!! COMMENTED AREA SINCE INTRODUCTION OF paral_kgb
!! BECAUSE USING FFT// imply to use MPI
!!$#if  !defined MPI
!!$!!$ if(dtset%paral_kgb == 1) then
!!$! TEST THE PARAMETERS OF THE 2D GRID
!!$! =====================================================
!!$
!!$  if(mpi_enreg%nproc_fft*mpi_enreg%nproc_band /= mpi_enreg%nproc)then
!!$   write(6,'(8a,i5,a,i5,a,i5)') ch10,&
!!$&   ' initmpi_grid : BUG -',ch10,&
!!$&   '  The number of band*FFT processors, npband*npfft, should be',ch10,&
!!$&   '  equal to the total number of processors, nproc.',ch10,&
!!$&   '  However, npband=',mpi_enreg%nproc_band,' npfft =',mpi_enreg%nproc_fft,&
!!$&   ' and nproc=',mpi_enreg%nproc
!!$   call leave_new('PERS')
!!$  end if
!!$
!!$  write(6,*) 'npfft and npband',&
!!$              & mpi_enreg%nproc_fft,&
!!$              & mpi_enreg%nproc_band
!!$
!!$! CREATE THE 2D GRID
!!$! ==================================================
!!$  mpi_enreg%dimcart=2
!!$
!!$  allocate(mpi_enreg%sizecart(mpi_enreg%dimcart))
!!$  allocate(periode           (mpi_enreg%dimcart))
!!$  allocate(mpi_enreg%coords  (mpi_enreg%dimcart))
!!$
!!$  mpi_enreg%sizecart(1)=mpi_enreg%nproc_fft
!!$  mpi_enreg%sizecart(2)=mpi_enreg%nproc_band
!!$
!!$  periode(:)=.false.
!!$  reorder   =.false.
!!$
!!$! create the cartesian grid with commcart as a communicator.
!!$  call MPI_CART_CREATE(MPI_COMM_WORLD,mpi_enreg%dimcart,mpi_enreg%sizecart,periode,&
!!$&  reorder,mpi_enreg%commcart,ierr)
!!$
!!$! Find the index and coordinates of the current  processor
!!$  call MPI_COMM_RANK(mpi_enreg%commcart, mpi_enreg%me_cart, ierr)
!!$  call MPI_CART_COORDS(mpi_enreg%commcart, mpi_enreg%me_cart,  mpi_enreg%dimcart, &
!!$&  mpi_enreg%coords, ierr)
!!$
!!$! Create the communicator for space (Fourier) distribution
!!$  allocate(keepdim(mpi_enreg%dimcart))
!!$
!!$  keepdim(1)=.true.
!!$  keepdim(2)=.false.
!!$  call MPI_CART_SUB(mpi_enreg%commcart, keepdim, mpi_enreg%comm_fft,ierr)
!!$  call MPI_COMM_SIZE(mpi_enreg%comm_fft,np_test, ierr)
!!$
!!$! Create the communicator for band distribution
!!$  keepdim(1)=.false.
!!$  keepdim(2)=.true.
!!$  call MPI_CART_SUB(mpi_enreg%commcart, keepdim, mpi_enreg%comm_band,ierr)
!!$
!!$! Define the correspondance with the fft
!!$  mpi_enreg%me_fft  = mpi_enreg%coords(1)
!!$  mpi_enreg%me_band = mpi_enreg%coords(2)
!!$
!!$
!!$! WRITE THE PROCESSORS COORDINATES IN THE 2D GRID
!!$! ==================================================
!!$  write(6,*) 'in initmpi_grid: me_fft and me_band are',mpi_enreg%me_fft,mpi_enreg%me_band
!!$
!!$ endif
!!$#endif
!!$#endif

!DEBUG
! write(6,*)' initmpi_grid : exit'
!ENDDEBUG

end subroutine initmpi_grid
!!***
