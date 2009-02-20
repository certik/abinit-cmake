!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_skip
!! NAME
!! hdr_skip
!!
!! FUNCTION
!! Skip wavefunction or density file header, after having rewound the file.
!! Two instances of the hdr_skip routines are defined :
!!  hdr_skip_int to which only the unit number is given
!!  hdr_skip_wfftype to which a wffil datatype is given
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG,MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  unit = number of unit to be read
!!
!! OUTPUT
!!  ierr = error code returned by the MPI calls
!!
!! SIDE EFFECTS
!!
!! NOTES
!! No checking performed, since hdr_skip is assumed to be used only
!! on temporary wavefunction files.
!! This initialize further reading and checking by rwwf
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine hdr_skip_int(unitfi,ierr)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13io_mpi, except_this_one => hdr_skip_int
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: unitfi
 integer, intent(out) :: ierr

!Local variables-------------------------------
 type(wffile_type) :: wff

! *************************************************************************

!Use default values for wff
 wff%unwff=unitfi
 wff%accesswff=0
 wff%me=0
 wff%master=0
!Then, transmit to hdr_skip_wfftype
 call hdr_skip_wfftype(wff,ierr)
end subroutine hdr_skip_int

!--------------------------------------------------------------------------------

subroutine hdr_skip_wfftype(wff,ierr)

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
 type(wffile_type),intent(inout) :: wff
 integer, intent(out) :: ierr

!Local variables-------------------------------
 integer :: headform,mu,npsp,unit,usepaw
 integer :: integers(17)
 character(len=6) :: codvsn
#if defined MPI_IO
           integer(abinit_offset)  :: posit,positloc
           integer :: fhwf,nbOct_ch,nbOct_int
           integer ::  delim_record
           integer  :: statux(MPI_STATUS_SIZE)
#endif


! *************************************************************************

 unit=wff%unwff
 ierr=0

 if( wff%accesswff==0 .or. &
&   (wff%accesswff==-1.and.wff%master==wff%me) ) then

  rewind (unit)

! Pick off headform from WF file
  read(unit) codvsn,headform                  ! XG040806 This does not work, but I do not understand why ?!
! read(unit) integers(1),headform             ! This works ...
!                                             ! MT012408 Because codvsn is char*6 and not an integer !
  if(headform==1   .or. headform==2   .or. &
&    headform==51  .or. headform==52  .or. &
&    headform==101 .or. headform==102       ) headform=22

  if (headform<44) then
   read (unit) integers(1:13),npsp
  else
   read (unit) integers(1:13),npsp,integers(15:17),usepaw
  end if

! Skip rest of header records
  do mu=1,2+npsp
   read (unit)
  end do
  if ((headform>=44).and.(usepaw==1)) then
   read (unit)
   read (unit)
  end if

#if defined MPI_IO
           else if(wff%accesswff==1)then

            fhwf  = wff%fhwff
            nbOct_ch = wff%nbOct_ch
            nbOct_int = wff%nbOct_int

!           Causes all previous writes to be transferred to the storage device
            call MPI_FILE_SYNC(fhwf,ierr)
            posit = 0

!           Reading the first record of the file ------------------------------------
!           read (unitfi)   codvsn,headform,..............
!           Pick off headform from WF file
            call MPI_FILE_READ_AT_ALL(fhwf,posit,delim_record,1, MPI_INTEGER,statux,ierr)
            positloc  = posit + nbOct_int+ nbOct_ch*6

            call MPI_FILE_READ_AT_ALL(fhwf,positloc,headform,1,MPI_INTEGER,statux,ierr)

            posit = posit + delim_record + 2 * nbOct_int

            if(headform==1   .or. headform==2   .or. &
&              headform==51  .or. headform==52  .or. &
&              headform==101 .or. headform==102       ) headform=22

!           Reading the second record of the file ------------------------------------
!           read(unitfi) bantot, hdr%date, hdr%intxc.................
!           Pick off psp and usepaw from WF file
            call MPI_FILE_READ_AT_ALL(fhwf,posit,delim_record,1, MPI_INTEGER,statux,ierr)
            positloc  = posit + nbOct_int *14

            call MPI_FILE_READ_AT_ALL(fhwf,positloc,npsp,1,MPI_INTEGER,statux,ierr)
            positloc  = positloc + nbOct_int
            if ( headform >= 44) then
!            read  usepaw
             positloc = positloc +  nbOct_int *3
             call MPI_FILE_READ_AT_ALL(fhwf,positloc,usepaw,1,MPI_INTEGER,statux,ierr)
            end if

            posit = posit + delim_record + 2 * nbOct_int

!           Reading the rest of the file ---------------------------------------------
            do mu=1,2+npsp
             call MPI_FILE_READ_AT_ALL(fhwf,posit,delim_record,1,MPI_INTEGER,statux,ierr)
             posit = posit + delim_record + 2 * nbOct_int
            end do
            if ((headform>=44).and.(usepaw==1)) then
             call MPI_FILE_READ_AT_ALL(fhwf,posit,delim_record,1,MPI_INTEGER,statux,ierr)
             posit = posit + delim_record + 2 * nbOct_int
             call MPI_FILE_READ_AT_ALL(fhwf,posit,delim_record,1,MPI_INTEGER,statux,ierr)
             posit = posit + delim_record + 2 * nbOct_int
            end if

            wff%offwff= posit

#endif

 end if

! in case of wff%accesswff == 2 NetCDF there is no skipping to do!
! write (6,*) 'skip n proc  ',wff%me ,' wff%offwff = ', wff%offwff
end subroutine hdr_skip_wfftype
!!***
