!{\src2tex{textfont=tt}}
!!****f* ABINIT/initmpi_fft
!! NAME
!! initmpi_fft
!!
!! FUNCTION
!! Initialize the mpi informations for FFT or BAND-FFT parallelism
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
!!  mpi_enreg=informations about MPI parallelization
!!
!! OUTPUT
!!  Not up to date ; to be updated !
!!  mpi_enreg=informations about MPI parallelization
!!    mpi_enreg%fft_comm(nkpt)=comm array of FFT set
!!    mpi_enreg%fft_group(nkpt)=group array of FFT set
!!    mpi_enreg%me_fft=index of the processor in the FFT set
!!    mpi_enreg%nproc_fft=number of processors int the FFT set
!!
!! SIDE EFFECTS
!!
!! TODO
!!
!! PARENTS
!!      gstate,invars2m,loper3,respfn
!!
!! CHILDREN
!!      leave_new,mpi_cart_coords,mpi_cart_create,mpi_cart_sub,mpi_comm_create
!!      mpi_comm_group,mpi_comm_rank,mpi_comm_size,mpi_group_incl
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initmpi_fft(dtset,mpi_enreg)

 use defs_basis
 use defs_infos
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => initmpi_fft
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
           integer :: nkpt,nsppol
          !Variables introduced for MPI version
           integer :: iblock,ierr,ifft,iikpt,iiproc,iproc,iproc_max,iproc_min,irank,isppol
           integer :: iikpt_modulo
           !Variables introduced for the bandFFT version
           logical :: reorder
           logical, allocatable :: periode(:), keepdim(:)
           integer, allocatable :: coords(:)
           integer :: np_fft, np_band, np_test,np_kpt, nproc_tmp
           integer,allocatable :: ranks(:)
#endif
! ***********************************************************************

!DEBUG
! write(6,*)' initmpi_fft : enter'
!stop
!ENDDEBUG

#if defined MPI
  if(dtset%paral_kgb == 1) then
    mpi_enreg%mode_para='b'
    mpi_enreg%fft_option_lob=1
    mpi_enreg%paralbd=0
    mpi_enreg%paral_fft=1
    if (dtset%fft_opt_lob /= 0) mpi_enreg%fft_option_lob=dtset%fft_opt_lob
  endif
#endif
!Set up information for FFT parallelism

#if defined MPI
        nkpt=dtset%nkpt
        nsppol=dtset%nsppol
        call MPI_COMM_SIZE(MPI_COMM_WORLD,mpi_enreg%nproc,ierr)
        call MPI_COMM_RANK(MPI_COMM_WORLD,mpi_enreg%me,ierr)
#endif
#if defined MPI
  if(dtset%paral_kgb == 1) then
! modif le 17/09/2004 sera a remettre
!              if(modulo(mpi_enreg%nproc,nkpt)/=0)then
!                   write(message,'(6a,i5,a,i5)') ch10,&
!&   ' initmpi_fft : BUG -',ch10,&
!&   '  The number of processors, nproc, should be',&
!&   '  a multiple of the number of nkpt for the FFT, nkpt.',&
!&   '  However, nproc=',mpi_enreg%nproc,' and nkpt=',nkpt
!                   call wrtout(06,message,'COLL')
!                   call leave_new('COLL')
!                end if
  endif
#endif
#if defined MPI
  if(dtset%paral_kgb == 1) then
! modif le 17/09/2004 sera a remettre
!!!            mpi_enreg%nproc_fft=mpi_enreg%nproc/nkpt
             mpi_enreg%nproc_fft=mpi_enreg%nproc
  endif
#endif

#if defined MPI
  if(dtset%paral_kgb == 0) then
             mpi_enreg%nproc_fft = 1
             mpi_enreg%num_group_fft=0
  endif
#endif

#if defined MPI
            if(mpi_enreg%nproc_fft == 0) mpi_enreg%nproc_fft=1
!           Creation of groups of communicators
            allocate(mpi_enreg%fft_comm(nkpt*nsppol))
            allocate(mpi_enreg%fft_group(nkpt*nsppol))
            call MPI_COMM_GROUP(MPI_COMM_WORLD,mpi_enreg%world_group,ierr)
            iproc=1
#endif

#if defined MPI
  if(dtset%paral_kgb == 1) then
            allocate(ranks(mpi_enreg%nproc_fft))
            do isppol=1,nsppol
             do iikpt=1,nkpt
! a enlever plus tard
!              iproc=1
! in which block I am ?
!              iblock=mpi_enreg%me / mpi_enreg%nproc_fft +1
              do irank=iproc,iproc + mpi_enreg%nproc_fft -1
!               ranks(irank)=irank+mpi_enreg%nproc_fft*(iblock-1)-1
               ranks(irank)=irank-1
               if (ranks(irank)==mpi_enreg%me) then
                mpi_enreg%num_group_fft=iikpt
               end if
              end do
!!! a remettre plus tard
!!!!               iproc=iproc+mpi_enreg%nproc_fft
              call MPI_GROUP_INCL(mpi_enreg%world_group,mpi_enreg%nproc_fft,  &
&              ranks,mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt),ierr)
              call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt), &
&              mpi_enreg%fft_comm(iikpt+(isppol-1)*nkpt),ierr)

!             XG070810 : Next line from M. Beland, needed to release memory, see mail from Paul Boulanger 070522
!             Seems very strange ...
              call MPI_GROUP_FREE(mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt),ierr)
!             XG070810 : End of modif

             end do
            end do
            call MPI_COMM_RANK(mpi_enreg%fft_comm(mpi_enreg%num_group_fft),&
&              mpi_enreg%me_fft,ierr)
            if (mpi_enreg%me_fft==0) then
                    mpi_enreg%master_fft=mpi_enreg%me
                else
                mpi_enreg%master_fft=-1
            end if
            call MPI_COMM_SIZE(mpi_enreg%fft_comm(mpi_enreg%num_group_fft),&
&              mpi_enreg%nproc_fft,ierr)
            deallocate(ranks)
! write(6,*) mpi_enreg%mode_para

 if(mpi_enreg%mode_para=='b') then

    mpi_enreg%nproc_fft  = dtset%npfft
    mpi_enreg%nproc_band = dtset%npband
    mpi_enreg%nproc_kpt  = dtset%npkpt
    mpi_enreg%bandpp     = dtset%bandpp


  if(modulo(dtset%ngfft(2),mpi_enreg%nproc_fft)/=0)then
   write(6,'(8a,i5,a,i5)') ch10,&
&   ' initmpi_fft : BUG -',ch10,&
&   '  The number of FFT processors, npfft, should be',ch10,&
&   '  a multiple of the number of ngfft(2).',ch10,&
&   '  However, npfft=',mpi_enreg%nproc_fft,' and ngfft(2)=',dtset%ngfft(2)
   call leave_new('PERS')
  end if

  do iikpt=1,nkpt*nsppol

   iikpt_modulo = modulo(iikpt,nkpt)+1

   if ( (dtset%istwfk(iikpt_modulo)==2) .and. (dtset%ngfft(7)==401) )then

      if (  (mpi_enreg%bandpp==0) .or. &
           ((mpi_enreg%bandpp/=1) .and. (modulo(mpi_enreg%bandpp,2)/=0)) ) then
         write(6,'(6a,i5)') ch10,&
   &   ' initmpi_fft : BUG -',ch10,&
   &   '  The number bandpp should be 1 or a multiple of 2',ch10,&
   &   '  However, bandpp=',mpi_enreg%bandpp
         call leave_new('PERS')
      endif


      if(modulo(dtset%nband(iikpt),mpi_enreg%nproc_band*mpi_enreg%bandpp)/=0)then
         write(6,'(8a,i5,a,i5)') ch10,&
   &   ' initmpi_fft : BUG -',ch10,&
   &   '  The number of band for the k-point, nband_k, should be',ch10,&
   &   '  a multiple of the number nproc_band*bandpp.',ch10,&
   &   '  However, nband_k=',dtset%nband(iikpt),' and nproc_band*bandpp=', &
   &      mpi_enreg%nproc_band* mpi_enreg%bandpp
         call leave_new('PERS')
      endif

   elseif ((dtset%istwfk(iikpt_modulo)==2) .and. (dtset%ngfft(7)==400)) then
      write(6,'(3a)') ch10,&
           &   ' initmpi_fft : BUG -',ch10,&
           &   '  The fftalg=400 with istwfk=2 is not valid'
      call leave_new('PERS')

   else
      if(modulo(dtset%nband(iikpt),mpi_enreg%nproc_band)/=0)then
         write(6,'(8a,i5,a,i5)') ch10,&
   &    ' initmpi_fft : BUG -',ch10,&
   &    '  The number of band processors, npband, should be',ch10,&
   &    '  a multiple of the number of nband.',ch10,&
   &    '  However, npband=',mpi_enreg%nproc_band,' and nband=',dtset%nband
         call leave_new('PERS')
      end if

      if ((mpi_enreg%bandpp/=1)) then
         write(6,'(4a,i5,2a,i5,2a,i5)') ch10,&
   &   ' initmpi_fft : BUG -',ch10,&
   &   '  The number bandpp should be 1 with fftalg=',dtset%ngfft(7),ch10,&
   &    ' and istwfk=',dtset%istwfk(iikpt_modulo),ch10,&
   &   '  However, bandpp=',mpi_enreg%bandpp
         call leave_new('PERS')
      endif

   end if
  end do

  if (mpi_enreg%paral_compil_kpt==1) then
     if(modulo(nkpt*nsppol,mpi_enreg%nproc_kpt)/=0)then
        write(6,'(8a,i5,a,i5)') ch10,&
             &   ' initmpi_fft : BUG -',ch10,&
             &   '  The number of KPT processors, npkpt, should be',ch10,&
             &   '  a multiple of the number of nkpt*nsppol.',ch10,&
             &   '  However, npkpt=',mpi_enreg%nproc_kpt,' and nkpt*nsppol=',nkpt*nsppol
        call leave_new('PERS')
     end if
  end if

  call initmpi_grid(dtset,mpi_enreg)

  if (mpi_enreg%paral_compil_kpt==1) then
     write(6,*) 'in initmpi:me_fft, me_band, me_kpt are',&
          &            mpi_enreg%me_fft,mpi_enreg%me_band,mpi_enreg%me_kpt
  else
     write(6,*) 'in initmpi:me_fft, me_band are',&
          &            mpi_enreg%me_fft,mpi_enreg%me_band
  end if

 end if
 endif
#endif


#if defined MPI
 if(dtset%paral_kgb ==0.and.mpi_enreg%paral_fft==1) then
            mpi_enreg%master_fft=-1
            do isppol=1,nsppol
                    do iikpt=1,nkpt
                     if (mpi_enreg%parareel==0) then
                        iproc_min=minval(mpi_enreg%proc_distrb(iikpt,:,isppol))
                        iproc_max=maxval(mpi_enreg%proc_distrb(iikpt,:,isppol))
                        if (mpi_enreg%me == iproc_min) then
                                mpi_enreg%master_fft=mpi_enreg%me
                        end if
                           allocate(ranks(iproc_max-iproc_min+1))
                        iiproc=1
                        do iproc=iproc_min,iproc_max
                                 ranks(iiproc)=iproc
                                iiproc=iiproc+1
                        end do

!                       With MPI on SGI machine "Spinoza", there is a limitation
!                       of the number of groups that can be defined. When FFT // is off
!                       paral_kgb=0, there is actually no reason to define
!                       the following groups (in the present implementation). So,
!                       these lines can be skipped.
#if !defined FC_MIPSPRO
                          call MPI_GROUP_INCL(mpi_enreg%world_group,&
&                        iproc_max-iproc_min+1,&
&                        ranks,mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt),ierr)
                        call MPI_COMM_CREATE(MPI_COMM_WORLD, &
&                        mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt),&
&                        mpi_enreg%fft_comm(iikpt+(isppol-1)*nkpt),ierr)
!             XG070810 : Next line from M. Beland, needed to release memory, see mail from Paul Boulanger 070522
!             Seems very strange ...
              call MPI_GROUP_FREE(mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt),ierr)
!             XG070810 : End of modif

#endif
                        deallocate(ranks)

                    else
                        iproc_min=mpi_enreg%proc_distrb_para(mpi_enreg%ipara,iikpt)
                        allocate(ranks(1))
                        ranks(1)=iproc_min
                        call MPI_GROUP_INCL(mpi_enreg%world_group,&
&                        1,ranks,mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt),ierr)
                        call MPI_COMM_CREATE(MPI_COMM_WORLD, &
&                        mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt),&
&                        mpi_enreg%fft_comm(iikpt+(isppol-1)*nkpt),ierr)
!             XG070810 : Next line from M. Beland, needed to release memory, see mail from Paul Boulanger 070522
!             Seems very strange ...
              call MPI_GROUP_FREE(mpi_enreg%fft_group(iikpt+(isppol-1)*nkpt),ierr)
!             XG070810 : End of modif
                        deallocate(ranks)
                    end if
                end do
            end do
!            mpi_enreg%num_group_fft=
 endif
#endif


#if defined MPI
! creation of master_communicator
           if (mpi_enreg%paral_fft == 0) then
                allocate(ranks(mpi_enreg%nproc))
                do iproc=0,mpi_enreg%nproc-1
                        ranks(iproc+1)=iproc
                end do
!This did not work with the g95 compiler XG081106: mpi_enreg%nproc was replaced by 0 . Very strange
!                call MPI_GROUP_INCL(mpi_enreg%world_group,mpi_enreg%nproc,ranks, &
!&                mpi_enreg%fft_master_group,ierr)
                nproc_tmp=mpi_enreg%nproc
                call MPI_GROUP_INCL(mpi_enreg%world_group,nproc_tmp,ranks, &
&                mpi_enreg%fft_master_group,ierr)

                call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_enreg%fft_master_group,&
&                mpi_enreg%fft_master_comm,ierr)
                deallocate(ranks)
                else
! only one proc per group_fft
!!! a remettre quand on fera la merge des MPI et MPI_FFT
! pour l'instant, le communicateur master de chaque FFT est toujours le proc 0
!                    allocate(ranks(nkpt*nsppol))
!                do isppol=1,nsppol
!                            do iikpt=1,nkpt
!                                ranks(iikpt+(isppol-1)*nkpt)=minval(mpi_enreg%proc_distrb(iikpt,:,isppol))
!                        end do
!                end do
                allocate(ranks(1))
                ranks(1)=0
!!!!                  call MPI_GROUP_INCL(mpi_enreg%world_group,nkpt*nsppol,ranks, &
!!!!&                mpi_enreg%fft_master_group,ierr)
                  call MPI_GROUP_INCL(mpi_enreg%world_group,1,ranks, &
&                mpi_enreg%fft_master_group,ierr)
!               call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_enreg%fft_master_group,&
!&                mpi_enreg%fft_master_comm,ierr)
                mpi_enreg%fft_master_comm=MPI_COMM_SELF
                   deallocate(ranks)
           end if
#endif

  if (mpi_enreg%paral_compil_fft==1 .and. dtset%mgfft /=0) then
!  creation of arrays for FFT parallelization
   allocate(mpi_enreg%nplanes_fft(dtset%nkpt))
!   write(6,*)'initmpi_fft dtset%ngfft(2)',dtset%ngfft(2)
   allocate(mpi_enreg%ind_fft_planes(dtset%nkpt,dtset%ngfft(2)))
!   write(6,*)'initmpi_fft alloc taille',dtset%ngfft(2)
  end if

!DEBUG
! write(6,*)' initmpi_fft : exit'
!ENDDEBUG

end subroutine initmpi_fft
!!***
