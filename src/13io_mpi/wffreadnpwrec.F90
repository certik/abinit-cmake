!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffReadNpwRec
!! NAME
!! WffReadNpwRec
!!
!! FUNCTION
!! This subroutine read the npw record of a wavefunction file
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (XG,MB,MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! wff= structured info about the wavefunction file
!!  wff%access == -1 and wf%master == Wff%me:
!!     read binary data
!!  wff%accesswff == 0:
!!     read binary data
!!  wff%accesswff == 1:
!!     use MPI/IO routines (MPIO defined) 
!!  wff%accesswff == 2:
!!     read netcdf format (NETCDF defined)
!! ikpt= the i-th kpoint.
!! isppol= the given spin polarisation element.
!!
!! OUTPUT
!! ierr=error code (iostat integer from read statement)
!! nband_disk=number of bands
!! npw=number of plane waves
!! nspinor=number of spinorial components of the wavefunctions
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      nstdy3,nstwf3,nstwf4,rwwf,vtowfk3,wfkfermi3
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,etsf_io_low_read_dim,etsf_io_low_read_var
!!      handle_ncerr,leave_new,wrtout,xderivereadval,xderiverrecend
!!      xderiverrecinit
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffReadNpwRec(ierr,ikpt,isppol,nband_disk,npw,nspinor,wff)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_ETSF_IO
 use etsf_io
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi, except_this_one => WffReadNpwRec
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer,intent(in)  :: ikpt, isppol
 integer,intent(out) :: ierr,nband_disk,npw,nspinor

!Local variables-------------------------------
 character(len=500) :: message
!no_abirules
#if defined HAVE_NETCDF
 integer :: ncerr,ncid_hdr,npwarr_id,nspinor_id,nband_id
 character(len=500) :: dummy_name
#endif
#if defined HAVE_ETSF_IO
  type(etsf_io_low_error) :: error
  logical                 :: lstat
  character(len = etsf_io_low_error_len)   :: errmess
#endif

! *************************************************************************

 ierr=0
!print *,"0",wff%accesswff,wff%master,wff%me
 if( wff%accesswff == 0   .or.                     &
&   (wff%accesswff ==-1 .and. wff%master==wff%me) ) then
  read (wff%unwff,iostat=ierr) npw,nspinor,nband_disk
! print *,"a",npw,nspinor,nband_disk

#if defined MPI_IO
 else if(wff%accesswff==1)then
  call xderiveRRecInit(wff,ierr)
  call xderiveReadVal(wff,npw)
  call xderiveReadVal(wff,nspinor)
  call xderiveReadVal(wff,nband_disk)
  call xderiveRRecEnd(wff,ierr)
#endif

#if defined HAVE_NETCDF
 else if (wff%accesswff == 2) then
  ncid_hdr = wff%unwff
  ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nspinor",dimid=nspinor_id)
  call handle_ncerr(ncerr,"inquire nspinor")
  ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nspinor_id,&
&       name=dummy_name,len=nspinor)
  call handle_ncerr(ncerr,"get nspinor")

  ncerr = nf90_inq_varid(ncid=ncid_hdr,name="npwarr",varid=npwarr_id)
  call handle_ncerr(ncerr,"inquire npw")
  ncerr = nf90_get_var(ncid=ncid_hdr,varid=npwarr_id,values=npw,start=(/ikpt/))
  call handle_ncerr(ncerr,"get npw")

  ncerr = nf90_inq_varid(ncid=ncid_hdr,name="nband",varid=nband_id)
  call handle_ncerr(ncerr,"inquire nband_disk")
  ncerr = nf90_get_var(ncid=ncid_hdr,varid=nband_id,values=nband_disk,start=(/ikpt, isppol/))
  call handle_ncerr(ncerr,"get nband_disk")
#endif
#if defined HAVE_ETSF_IO
 else if (wff%accesswff == 3) then
  call etsf_io_low_read_dim(wff%unwff, "number_of_spinor_components", &
                          & nspinor, lstat, error_data = error)
  if (.not. lstat) then
    call etsf_io_low_error_to_str(errmess, error)
    write(message, "(A,A,A,A)") ch10, " WffReadNpwRec: ERROR -", ch10, &
                              & errmess(1:min(475, len(errmess)))
    call wrtout(std_out, message, 'COLL')
    call leave_new('COLL')
  end if

  call etsf_io_low_read_var(wff%unwff, "number_of_coefficients", &
                          & npw, lstat, start = (/ ikpt /), error_data = error)
  if (.not. lstat) then
    call etsf_io_low_error_to_str(errmess, error)
    write(message, "(A,A,A,A)") ch10, " WffReadNpwRec: ERROR -", ch10, &
                              & errmess(1:min(475, len(errmess)))
    call wrtout(std_out, message, 'COLL')
    call leave_new('COLL')
  end if

  call etsf_io_low_read_var(wff%unwff, "number_of_states", &
       & nband_disk, lstat, start = (/ ikpt, isppol /), error_data = error)
  if (.not. lstat) then
    call etsf_io_low_error_to_str(errmess, error)
    write(message, "(A,A,A,A)") ch10, " WffReadNpwRec: ERROR -", ch10, &
                              & errmess(1:min(475, len(errmess)))
    call wrtout(std_out, message, 'COLL')
    call leave_new('COLL')
  end if
#endif

! else
!! Error: These values are not allowed.
!! Set intent(out) variables
!  ierr=1
!  nband_disk=0
!  npw=0
!  nspinor=0
!! print *,"A",npw,nspinor,nband_disk
!  write(message,'(a,a,a,a,a,a,i5,a)')ch10,&
!&     ' WffReadNpwRec : BUG -',ch10,&
!&     '  Reading option of WffReadNpwRec.',ch10,&
!&     '  The value of wff%accesswff=',&
!&     wff%accesswff,' is not allowed.'
!  call wrtout(6,message,'PERS')
!  call leave_new('PERS')

! else
!! Error: These values are not allowed.
!! Set intent(out) variables
!  ierr=1
!  nband_disk=0
!  npw=0
!  nspinor=0
!! print *,"A",npw,nspinor,nband_disk
!  write(message,'(a,a,a,a,a,a,i5,a)')ch10,&
!&     ' WffReadNpwRec : BUG -',ch10,&
!&     '  Reading option of WffReadNpwRec.',ch10,&
!&     '  The value of wff%accesswff=',&
!&     wff%accesswff,' is not allowed.'
!  call wrtout(6,message,'PERS')
!  call leave_new('PERS')

 end if

end subroutine WffReadNpwRec
!!***
