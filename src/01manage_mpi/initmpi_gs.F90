!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_gs
!! NAME
!! initmpi_gs
!!
!! FUNCTION
!! Initialize the mpi informations for the ground-state datasets
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (AR, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  mpi_enreg=informations about MPI parallelization
!!   mpi_enreg%paralbd=option for parallelisation over the bands
!!   If the cpp option MPI is activated, also initialize
!!    mpi_enreg%proc_distrb(nkpt,mband,nsppol)  array that
!!     describes the work of each processor (allocated here)
!!    mpi_enreg%nproc=number of processors
!!    mpi_enreg%me=my processor number
!!
!! TODO
!!
!! PARENTS
!!      gstate,pstate
!!
!! CHILDREN
!!      distrb2,leave_new,mpi_comm_create,mpi_comm_group,mpi_comm_rank
!!      mpi_comm_size,mpi_group_incl,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initmpi_gs(dtset,mpi_enreg)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => initmpi_gs
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(inout) :: mpi_enreg

!Local variables-------------------------------
!no_abirules
#if defined MPI 
           integer :: mband,nbdblock,nkpt,bakkpt,nsppol
          !Variables introduced for MPI version
           integer :: group,ierr,iikpt,iisppol,ipara,irank,max_proc,min_proc,nband_k
           integer,allocatable :: ranks(:)
           ! G parallelism
           integer :: ikpt, minprock, maxprock, numprock, iprock, mprock
           character(len=500) :: message
#endif

! ***********************************************************************

!DEBUG
! write(6,*)' initmpi_gs : enter'
!stop
!ENDDEBUG

 mpi_enreg%paralbd=0

!Set up information for FFT parallelism
 mpi_enreg%me_fft=0
 mpi_enreg%nproc_fft=1
 if(dtset%paral_kgb == 0) then
  mpi_enreg%paral_fft=0
 endif
!MPIWF Should be modified for real FFT parallelism ...

#if defined MPI 
           mband=dtset%mband
           nbdblock=dtset%nbdblock
           nkpt=dtset%nkpt
           nsppol=dtset%nsppol
           if (mpi_enreg%parareel == 0) then
              allocate(mpi_enreg%proc_distrb(nkpt,mband,nsppol))

              call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_enreg%nproc,ierr)
!call leave_new("COLL")
#if defined MPI 
  if(dtset%paral_kgb ==1) then
              mpi_enreg%nproc    = mpi_enreg%nproc_kpt
  endif
# endif
              
              if (nkpt >= mpi_enreg%nproc) then
                 mpi_enreg%paralbd=0
                 else
                 if (nbdblock == 1) then
                    mpi_enreg%paralbd=0
                    else
                    mpi_enreg%paralbd=nbdblock
                 end if
                 end if
                else
                !case parareel
                allocate(mpi_enreg%proc_distrb_para(0:mpi_enreg%npara-1,nkpt))
                call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_enreg%nproc,ierr)
           end if
#endif

#if defined MPI
  if(dtset%paral_kgb ==1) then 
! temporary
                mpi_enreg%proc_distrb(:,:,:)=0
  endif
#endif

#if defined MPI
!          This routine sets up the array mpi_enreg%proc_distrb,
!          and also mpi_enreg%me
           call distrb2(mband, dtset%nband, nkpt, nsppol, mpi_enreg)

           if (mpi_enreg%paralbd >= 1) then
!           Creation of groups of communicators
            allocate(mpi_enreg%kpt_comm(nkpt*nsppol),mpi_enreg%kpt_group(nkpt*nsppol))
            allocate(ranks(mpi_enreg%nproc_per_kpt))
            call MPI_COMM_GROUP(MPI_COMM_WORLD,mpi_enreg%world_group,ierr)
            do iisppol=1,nsppol
             do iikpt=1,nkpt
              group=iikpt+(iisppol-1)*nkpt
              nband_k=dtset%nband(iikpt+(iisppol-1)*nkpt)
              min_proc=minval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
              max_proc=maxval(mpi_enreg%proc_distrb(iikpt,1:nband_k,iisppol))
              do irank=1,mpi_enreg%nproc_per_kpt
               ranks(irank)=min_proc+irank-1
               if (ranks(irank)==mpi_enreg%me) then
                mpi_enreg%num_group=group
               end if
              end do
              call MPI_GROUP_INCL(mpi_enreg%world_group,mpi_enreg%nproc_per_kpt,ranks, &
&               mpi_enreg%kpt_group(group),ierr)
              call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_enreg%kpt_group(group), &
&               mpi_enreg%kpt_comm(group),ierr)
             end do
            end do
            if ((mpi_enreg%nproc > mpi_enreg%nproc_per_kpt*nkpt*nsppol) .and. &
&                (mpi_enreg%me >= mpi_enreg%nproc_per_kpt*nkpt*nsppol)) then
             mpi_enreg%num_group=0
             mpi_enreg%me_group=-1
             mpi_enreg%nproc_group=-1
            else
             call MPI_COMM_RANK(mpi_enreg%kpt_comm(mpi_enreg%num_group), &
&              mpi_enreg%me_group,ierr)
             call MPI_COMM_SIZE(mpi_enreg%kpt_comm(mpi_enreg%num_group), &
&              mpi_enreg%nproc_group,ierr)
            end if
            deallocate(ranks)
           end if
           
           if (mpi_enreg%parareel == 1) then
             write(*,*) 'ERROR'
!           Creation of groups of communicators for parareel
            allocate(mpi_enreg%kpt_comm_para(0:mpi_enreg%npara-1))
            allocate(mpi_enreg%kpt_group_para(0:mpi_enreg%npara-1))
            allocate(ranks(mpi_enreg%nproc_per_para))
            call MPI_COMM_GROUP(MPI_COMM_WORLD,mpi_enreg%world_group,ierr)
            do ipara=0,mpi_enreg%npara-1
              do irank=1,mpi_enreg%nproc_per_para
               ranks(irank)=mpi_enreg%proc_distrb_para(ipara,irank)
               if (ranks(irank)==mpi_enreg%me) then
                mpi_enreg%num_group_para=ipara
               end if
              end do
              call MPI_GROUP_INCL(mpi_enreg%world_group,mpi_enreg%nproc_per_para,ranks, &
&               mpi_enreg%kpt_group_para(ipara),ierr)
              call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_enreg%kpt_group_para(ipara), &
&               mpi_enreg%kpt_comm_para(ipara),ierr)
            end do
            deallocate(ranks)
           end if


#endif

#if defined MPI 
 if(dtset%paral_kgb ==1) then
 call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_enreg%nproc,ierr)
 endif
#endif

!DEBUG
!write(6,*)' initmpi_gs : exit '
!stop
!ENDDEBUG

end subroutine initmpi_gs
!!***
