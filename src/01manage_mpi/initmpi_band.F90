!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_band
!! NAME
!! initmpi_band
!!
!! FUNCTION
!! Initialize the mpi informations for band parallelism (paralbd=1)
!!
!! COPYRIGHT
!! Copyright (C) 2008-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  mpi_enreg= informations about MPI parallelization
!!  nband(nkpt*nsppol)= number of bands per k point, for each spin
!!  nkpt= number of k-points
!!  nsppol= 1 for unpolarized, 2 for polarized
!!
!! OUTPUT
!!  mpi_enreg=informations about MPI parallelization
!!    mpi_enreg%band_comm(nkpt*nsppol)=comm array of BAND set
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      mpi_comm_create,mpi_group_free,mpi_group_incl
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initmpi_band(mpi_enreg,nband,nkpt,nsppol)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: nkpt,nsppol
 integer,intent(in) :: nband(nkpt*nsppol)
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
#if defined MPI
 integer :: band_group,ib,ierr,ii,ikpt,iproc_min,iproc_max,irank,isppol,jj,nband_k,nbsteps,nrank,nstates
 integer,allocatable :: ranks(:)
#endif

! ***********************************************************************

!DEBUG
! write(6,*)' initmpi_band : enter'
!stop
!ENDDEBUG

 mpi_enreg%has_band_comm=0

#if defined MPI

 if (mpi_enreg%paralbd==1) then

  nstates=sum(nband(1:nkpt*nsppol))

  if(mpi_enreg%paral_compil_respfn == 1) then
   nbsteps=(nstates*mpi_enreg%ngroup_respfn)/mpi_enreg%nproc_kpt
   if (mod(nstates,mpi_enreg%nproc_kpt/mpi_enreg%ngroup_respfn)/=0) nbsteps=nbsteps+1
  else
   nbsteps=nstates/mpi_enreg%nproc_kpt
   if (mod(nstates,mpi_enreg%nproc_kpt)/=0) nbsteps=nbsteps+1
  end if

  if (nbsteps<maxval(nband(1:nkpt*nsppol))) then

   mpi_enreg%has_band_comm=1
   allocate(mpi_enreg%band_comm(nkpt*nsppol))

   do isppol=1,nsppol
    do ikpt=1,nkpt
     ii=ikpt+(isppol-1)*nkpt
     nband_k=nband(ii)
     if (nbsteps<nband_k) then
      iproc_min=minval(mpi_enreg%proc_distrb(ikpt,:,isppol))
      iproc_max=maxval(mpi_enreg%proc_distrb(ikpt,:,isppol))
      nrank=iproc_max-iproc_min+1
      allocate(ranks(nrank));jj=iproc_min-1
      do irank=1,nrank
       jj=jj+1;ranks(irank)=jj
      end do
      call MPI_GROUP_INCL(mpi_enreg%world_group,nrank,ranks,band_group,ierr)
      call MPI_COMM_CREATE(MPI_COMM_WORLD,band_group,mpi_enreg%band_comm(ii),ierr)
      call MPI_GROUP_FREE(band_group,ierr)
      deallocate(ranks)
     else
      mpi_enreg%band_comm(ii)=MPI_COMM_SELF
     end if
    end do
   end do

  end if
 end if

#endif
!DEBUG
! write(6,*)' initmpi_band : exit'
!ENDDEBUG

end subroutine initmpi_band
!!***
