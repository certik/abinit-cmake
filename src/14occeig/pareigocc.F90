!{\src2tex{textfont=tt}}
!!****f* ABINIT/pareigocc
!! NAME
!! pareigocc
!!
!! FUNCTION
!! This subroutine transmit to all processors, using MPI :
!!   - the eigenvalues and,
!!   - if ground-state, the occupation numbers
!!     (In fact, in the present status of the routine,
!!      occupation numbers are NOT transmitted)
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (XG, AR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  formeig=format of eigenvalues (0 for GS, 1 for RF)
!!  localrdwf=(for parallel case) if 1, the eig and occ initial values
!!            are local to each machine, if 0, they are on proc me=0.
!!  mband=maximum number of bands of the output wavefunctions
!!  mpi_enreg=informations about MPI parallelization
!!  nband(nkpt*nsppol)=desired number of bands at each k point
!!  nkpt=number of k points
!!  nsppol=1 for unpolarized, 2 for spin-polarized, output wf file processors,
!!         Warning : defined only when paralbd=1
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  eigen(mband*nkpt*nsppol)=eigenvalues (input or init to large number), (Ha)
!!  occ(mband*nkpt*nsppol)=occupation (input or init to 0.0)  NOT USED NOW
!!
!! NOTES
!! * The case paralbd=1 with formeig=0 is implemented, but not yet used.
!!
!! * The transmission of occ is not activated yet !
!!
!! * The routine takes the eigenvalues in the eigen array on one of the
!!   processors that possess the wavefunctions, and transmit it to all procs.
!!   If localrdwf==0, me=0 has the full array at start,
!!   If localrdwf==1, the transfer might be more complex.
!!
!! * This routine should not be used for RF wavefunctions, since
!!   it does not treat the eigenvalues as a matrix.
!!
!! PARENTS
!!      newkpt,wfsinp
!!
!! CHILDREN
!!      mpi_bcast,timab,wrtout,xcomm_init,xme_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pareigocc(eigen,formeig,localrdwf,mpi_enreg,mband,nband,nkpt,nsppol,occ)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: formeig,localrdwf,mband,nkpt,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: nband(nkpt*nsppol)
 real(dp),intent(inout) :: eigen(mband*(2*mband)**formeig*nkpt*nsppol)
 real(dp),intent(inout) :: occ(mband*nkpt*nsppol)

!Local variables-------------------------------
!BEGIN TF_CHANGES
!END TF_CHANGES
!scalars
 integer :: band_index,iband,ierr,ikpt,isppol,me,nbks,spaceComm
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:)

! *************************************************************************

!DEBUG
!write(6,*)' pareigocc : enter '
!stop
!ENDDEBUG

 if(mpi_enreg%paral_compil_kpt==1)then

! Init mpi_comm
  call xcomm_init(mpi_enreg,spaceComm)

! BEGIN TF_CHANGES
! Init me
  call xme_init(mpi_enreg,me)
! END TF_CHANGES

  if(localrdwf==0)then
   call xcast_mpi(eigen,0,spaceComm,ierr)

  else if(localrdwf==1)then

!  Prepare transmission of eigen (and occ)
   allocate(buffer1(2*mband**(formeig+1)*nkpt*nsppol))
   allocate(buffer2(2*mband**(formeig+1)*nkpt*nsppol))
   buffer1(:)=0.0_dp
   buffer2(:)=0.0_dp

   band_index=0
   do isppol=1,nsppol
    do ikpt=1,nkpt
     nbks=nband(ikpt+(isppol-1)*nkpt)

     if(mpi_enreg%paralbd==0)then

      if(formeig==0)then
       buffer1(2*band_index+1:2*band_index+nbks)=&
&       eigen(band_index+1:band_index+nbks)
!      buffer1(2*band_index+nbks+1:2*band_index+2*nbks)=&
!      &                  occ(band_index+1:band_index+nbks)
       band_index=band_index+nbks
      else if(formeig==1)then
       buffer1(band_index+1:band_index+2*nbks**2)=&
&       eigen(band_index+1:band_index+2*nbks**2)
       band_index=band_index+2*nbks**2
      end if

     else if(mpi_enreg%paralbd==1)then

!     Skip this k-point if not the proper processor
!     BEGIN TF_CHANGES
      if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nbks,isppol) &
&      -me))/=0) then
!      END TF_CHANGES
       if(formeig==0) then
        band_index=band_index+nbks
       else
        band_index=band_index+2*nbks**2
       end if
       cycle
      end if
!     Loop on bands
      do iband=1,nbks
       if (mpi_enreg%parareel == 0) then
!       BEGIN TF_CHANGES
        if(mpi_enreg%proc_distrb(ikpt, iband,isppol) /= me)cycle
!       END TF_CHANGES
       else
        if(mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt) /= mpi_enreg%me)cycle
       end if
       if(formeig==0)then
        buffer1(2*band_index+iband)=eigen(band_index+iband)
!       buffer1(2*band_index+iband+nbdks)=occ(band_index+iband)
       else if (formeig==1)then
        buffer1(band_index+(iband-1)*2*nbks+1: &
&        band_index+(iband-1)*2*nbks+2*nbks)=&
&        eigen(band_index+(iband-1)*2*nbks+1: &
&        band_index+(iband-1)*2*nbks+2*nbks)
       end if
      end do
      if(formeig==0)then
       band_index=band_index+nbks
      else
       band_index=band_index+2*nbks**2
      end if
     end if

    end do
   end do

!  Build sum of everything
   call timab(48,1,tsec)
   write(message, '(a)' ) ' pareigocc : MPI_ALLREDUCE'
   call wrtout(6,message,'COLL')
   if(formeig==0)band_index=band_index*2

   call xsum_mpi(buffer1,buffer2,band_index,spaceComm,ierr)

   call timab(48,2,tsec)

   band_index=0
   do isppol=1,nsppol
    do ikpt=1,nkpt
     nbks=nband(ikpt+(isppol-1)*nkpt)
     if(formeig==0)then
      eigen(band_index+1:band_index+nbks)=&
&      buffer2(2*band_index+1:2*band_index+nbks)
!     occ(band_index+1:band_index+nbks)=&
!     &     buffer2(2*band_index+nbks+1:2*band_index+2*nbks)
      band_index=band_index+nbks
     else if(formeig==1)then
      eigen(band_index+1:band_index+2*nbks**2)=&
&      buffer1(band_index+1:band_index+2*nbks**2)
      band_index=band_index+2*nbks**2
     end if
    end do
   end do

   deallocate(buffer1,buffer2)

  end if

 end if ! mpi_enreg%paral_compil_kpt==1

end subroutine pareigocc
!!***
