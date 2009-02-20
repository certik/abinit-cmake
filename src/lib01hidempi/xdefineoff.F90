!{\src2tex{textfont=tt}}
!!****f* ABINIT/xdefineOff
!! NAME
!! xdefineOff
!!
!! FUNCTION
!! In case of MPI I/O, define the offset for each processor.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  formeig option (format of the eigenvalues and occupations) :
!!   0 => ground-state format (initialisation of eigenvectors with
!!        random numbers, vector of eigenvalues, occupations are present)
!!   1 => respfn format (initialisation of eigenvectors with 0 s,
!!        hermitian matrix of eigenvalues)
!!  nkpt = number of k points
!!  nspinor = number of spinorial components of the wavefunctions
!!  nsppol = nsppol = number of channels for spin-polarization (1 or 2)
!!  nband(nkpt*nsppol) = number of bands at each k point, for each polarization
!!  npwarr(nkpt) = number of planewaves at each k point
!!  mpi_enreg <type(MPI_type)> = informations about MPI parallelization
!!
!! OUTPUT
!!  (no output)
!!
!! SIDE EFFECTS
!!  wff <type(wffile_type)> =
!!
!! PARENTS
!!      dyfnl3,eltfrkin3,eltfrnl3,energy,forstrnps,inwffil,inwffil3
!!      mkrho,mkrho3,nselt3,nstdy3,outwf,rhofermi3,uderiv,vtorho,vtorho3
!!
!! CHILDREN
!!      mpi_file_read_at,mpi_type_size
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xdefineOff(formeig,wff,mpi_enreg,nband,npwarr,nspinor,nsppol,nkpt)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in) ::  nsppol,nkpt,nspinor,formeig
 integer, intent(in) ::  nband(nkpt*nsppol),npwarr(nkpt)
 type(wffile_type),intent(inout) :: wff
 type(MPI_type),intent(in) :: mpi_enreg

!Local variables-------------------------------
 integer :: iproc
#if defined MPI_IO
           integer :: nband_k,npw_k,nproc,me,ierr,ipp
           integer :: nbrec,isppol,ikpt,nbint,nbreal,nbd,ippband
           integer :: nrecnpw,nreckg
           integer(abinit_offset),allocatable  :: offproc(:)
           integer(abinit_offset) :: pos_start
#endif

! *************************************************************************
! nbOct_int octet number of int value
! nbOct_dp octet number of dp value
! nbOct_ch octet number of character value
! lght_recs length of record

#if defined MPI_IO
           if(wff%accesswff==1)then
            call MPI_Type_size(MPI_INTEGER,wff%nbOct_int,ierr)
            call MPI_Type_size(MPI_DOUBLE_PRECISION,wff%nbOct_dp,ierr)

            nproc=mpi_enreg%nproc
            me=mpi_enreg%me
            pos_start=wff%offwff

            allocate(offproc(0:nproc))
            do iproc = 0,nproc
             offproc(iproc) = 0
            end do
            nbrec =2
            nrecnpw=3+nbrec

            do isppol=1,nsppol
             do ikpt=1,nkpt
              nband_k=nband(ikpt+(isppol-1)*nkpt)
              npw_k=npwarr(ikpt)
              if (mpi_enreg%parareel == 0) then
               iproc = mpi_enreg%proc_distrb(ikpt,1,isppol)
              else
               iproc = mpi_enreg%proc_distrb_para(mpi_enreg%ipara,ikpt)
              end if
              if ( mpi_enreg%paralbd ==1)then
                 iproc= mpi_enreg%proc_distrb(ikpt,1,isppol)

              endif
              ! record kg
              nreckg=nbrec+ wff%kgwff*3*npw_k

!             Record npw,nspinor,nband
!             Record kg
              offproc(iproc) = offproc(iproc) + wff%nbOct_int*(nrecnpw+nreckg)

              if ( formeig == 0) then
!             Records eigen,occ
               nbint=nbrec
               nbreal =  2 *nband_k
               offproc(iproc) = offproc(iproc) + (wff%nbOct_int*nbint  &
              &                + wff%nbOct_dp*nbreal)



!             Records cg

              offproc(iproc) = offproc(iproc) + (wff%nbOct_int*nbrec  &
              &                + wff%nbOct_dp*2*npw_k*nspinor)*nband_k

              ippband=iproc
               do nbd=1,nband_k
                  ipp=mpi_enreg%proc_distrb(ikpt,nbd,isppol)

                  if ( ipp /= ippband ) then
                   ippband=ipp
                   offproc(ippband)=offproc(ippband)+ wff%nbOct_int*(nrecnpw+nreckg)
               offproc(ippband) = offproc(ippband) + (wff%nbOct_int*nbint  &
              &                + wff%nbOct_dp*nbreal)
              offproc(ippband) = offproc(ippband) + (wff%nbOct_int*nbrec  &
              &                + wff%nbOct_dp*2*npw_k*nspinor)*nband_k
                     endif
                  enddo
             else if ( formeig == 1 ) then
              ! record eigen
              offproc(iproc) = offproc(iproc) + (wff%nbOct_int*2*nbrec  &
             &                + wff%nbOct_dp*2*npw_k*nspinor          &
             &                + wff%nbOct_dp*2*nband_k)*nband_k
              ippband=iproc
          do nbd=1,nband_k
                  ipp=mpi_enreg%proc_distrb(ikpt,nbd,isppol)
                  if ( ipp /= ippband ) then
                     ippband=ipp
                     offproc(ippband)=offproc(ippband)+ wff%nbOct_int*(nrecnpw+nreckg)
              offproc(ippband) = offproc(ippband) + (wff%nbOct_int*2*nbrec  &
             &                + wff%nbOct_dp*2*npw_k*nspinor          &
             &                + wff%nbOct_dp*2*nband_k)*nband_k
                     endif
                  enddo
             end if   ! formeig
             end do ! ikpt

            end do ! isppol

            pos_start=wff%offwff

            wff%offwff = pos_start



            if ( me/=0)then
             do iproc=0,me-1
              wff%offwff=wff%offwff+offproc(iproc)
             end do
          endif

            deallocate(offproc)
           end if ! accesswff
#endif

end subroutine xdefineOff
!!***

!------------------------------------------------------------------------------------

subroutine xderiveWRecEnd(wff,ierr)
 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
                 type(wffile_type),intent(inout) :: wff
                 integer,intent(out) ::  ierr

!Local variables-------------------------------
#if defined MPI_IO
                 integer(abinit_offset)  :: offset,posit
                 integer  :: statux(MPI_STATUS_SIZE)
                integer :: delim_record
#endif

! *************************************************************************

 ierr=0
#if defined MPI_IO
                offset=wff%off_recs
                posit=wff%offwff
                delim_record =wff%offwff-wff%off_recs-wff%nbOct_int
                ! we write the first word of the record
                call MPI_FILE_WRITE_AT(wff%fhwff, offset,delim_record,1 &
                & , MPI_INTEGER , statux, ierr)

                call MPI_FILE_WRITE_AT(wff%fhwff, posit,delim_record,1 &
                & , MPI_INTEGER , statux, ierr)

                wff%offwff = posit +wff%nbOct_int
#endif
end subroutine xderiveWRecEnd

subroutine xderiveWRecEnd_cs(wff,ierr,me)
 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) ::  ierr
 integer,intent(in) :: me
!Local variables-------------------------------
#if defined MPI_IO
 integer(abinit_offset)  :: offset,posit
 integer  :: statux(MPI_STATUS_SIZE)
 integer :: delim_record
#endif
! *************************************************************************
 ierr=0
#if defined MPI_IO
 ! Only one processor write this
 if (me == 0) then
    offset=wff%off_recs
    posit=wff%offwff
    delim_record =wff%offwff-wff%off_recs-wff%nbOct_int
    ! we write the first word of the record
    call MPI_FILE_WRITE_AT(wff%fhwff, offset,delim_record,1 &
         & , MPI_INTEGER , statux, ierr)
    call MPI_FILE_WRITE_AT(wff%fhwff, posit,delim_record,1 &
         & , MPI_INTEGER , statux, ierr)
 end if
 wff%offwff = wff%offwff +wff%nbOct_int
#endif
end subroutine xderiveWRecEnd_cs
!------------------------------------------------------------------------------

subroutine xderiveWRecInit(wff,ierr)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) ::  ierr

!Local variables-------------------------------
#if defined MPI_IO
                 integer ::  delim_record
                 integer(abinit_offset)  :: offset,posit
                 integer  :: statux(MPI_STATUS_SIZE)
#endif

! *************************************************************************

 ierr=0
#if defined MPI_IO
                ! init offset record
                call MPI_Type_size(MPI_INTEGER,wff%nbOct_int,ierr)
                wff%off_recs = wff%offwff
                posit = wff%off_recs
                delim_record =0
                ! we write the first word of the record
                call MPI_FILE_WRITE_AT(wff%fhwff, posit,delim_record,1 &
                & , MPI_INTEGER , statux, ierr)
                wff%offwff = wff%off_recs +wff%nbOct_int
#endif
end subroutine xderiveWRecInit


subroutine xderiveWRecInit_cs(wff,ierr,me)
 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) ::  ierr
 integer,intent(in) :: me

!Local variables-------------------------------
#if defined MPI_IO
                 integer ::  delim_record
                 integer(abinit_offset)  :: offset,posit
                 integer  :: statux(MPI_STATUS_SIZE)
#endif

! *************************************************************************
 ierr=0
#if defined MPI_IO
 ! Only one processor write this
 if (me == 0) then
    ! init offset record
    call MPI_Type_size(MPI_INTEGER,wff%nbOct_int,ierr)
    wff%off_recs = wff%offwff
    posit = wff%off_recs
    delim_record =0
    ! we write the first word of the record
    call MPI_FILE_WRITE_AT(wff%fhwff, posit,delim_record,1 &
         & , MPI_INTEGER , statux, ierr)
 end if
 wff%offwff = wff%offwff +wff%nbOct_int
 
#endif
end subroutine xderiveWRecInit_cs

!---------------------------------------------------------------------------------

subroutine xderiveRRecEnd(wff,ierr)
 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) ::  ierr

!Local variables-------------------------------

! *************************************************************************

 ierr=0
#if defined MPI_IO
!          Define offset end of record
           call MPI_Type_size(MPI_INTEGER,wff%nbOct_int,ierr)
           wff%offwff = wff%off_recs + wff%lght_recs +2 *wff%nbOct_int
#endif
end subroutine xderiveRRecEnd

!-------------------------------------------------------------------------------

subroutine xderiveRRecInit(wff,ierr)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(out) ::  ierr

!Local variables-------------------------------
#if defined MPI_IO
                 integer  :: statux(MPI_STATUS_SIZE)
                integer :: delim_record
#endif

! *************************************************************************

 ierr=0
#if defined MPI_IO
                ! init offset record
                call MPI_Type_size(MPI_INTEGER,wff%nbOct_int,ierr)
                wff%off_recs = wff%offwff
                call MPI_FILE_READ_AT(wff%fhwff, wff%offwff,delim_record,1 &
                & , MPI_INTEGER , statux, ierr)
                wff%offwff =  wff%offwff +wff%nbOct_int
                ! init lenght record
                wff%lght_recs = delim_record
#endif
end subroutine xderiveRRecInit
