!{\src2tex{textfont=tt}}
!!****f* ABINIT/outxfhist
!! NAME
!! outxfhist
!!
!! FUNCTION
!!  read/write xfhist
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  option =
!!   1: write
!!   2: read only nxfh
!!   3: read xfhist
!!  response =
!!   0: GS wavefunctions
!!   1: RF wavefunctions
!!  natom = number of atoms in unit cell
!!  mxfh = last dimension of the xfhist array
!!
!! OUTPUT
!!  ios = error code returned by read operations
!!
!! SIDE EFFECTS
!!  nxfh = actual number of (x,f) history pairs, see xfhist array
!!  wff2 = structured info for wavefunctions
!!  xfhist(3,natom+4,2,mxfh) = (x,f) history array, also including
!!   rprim and stress
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      handle_ncerr,xderivereadval,xderiverrecend,xderiverrecinit
!!      xderivewrecend,xderivewrecinit,xderivewriteval
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outxfhist(nxfh,natom,mxfh,xfhist,option,wff2,ios)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13io_mpi, except_this_one => outxfhist
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer          ,intent(in)    :: natom,mxfh,option
 integer          ,intent(inout) :: nxfh
 integer          ,intent(out)   :: ios
 real(dp)         ,intent(inout) :: xfhist(3,natom+4,2,mxfh)
 type(wffile_type),intent(inout)    :: wff2

!Local variables-------------------------------
 integer :: ixfh
 integer :: ncid_hdr,xfdim2
!no_abirules
#if defined MPI
           integer :: ierr
#endif
#if defined HAVE_NETCDF
           integer :: ncerr
           integer :: nxfh_id, mxfh_id, xfdim2_id, dim2inout_id, dimr3_id,xfhist_id
           integer :: nxfh_tmp,mxfh_tmp,xfdim2_tmp,dim2inout_tmp
#endif


! *************************************************************************

! DEBUG
! write(6,*)'outxfhist  : enter'
! ENDDEBUG
 ncid_hdr = wff2%unwff
 xfdim2 = natom+4


 if ( option == 1 ) then
    ios = 0

! Write the (x,f) history
#if !defined MPI
    if (wff2%accesswff /= 2) then
     write(unit=wff2%unwff)nxfh
     do ixfh=1,nxfh
      write(unit=wff2%unwff)xfhist(:,:,:,ixfh)
     end do
#if defined HAVE_NETCDF
    else if (wff2%accesswff == 2) then
! check if nxfh and xfhist are defined
     ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)

     if (ncerr /= NF90_NOERR) then
! need to define everything
       ncerr = nf90_redef (ncid=ncid_hdr)
       call handle_ncerr(ncerr," outxfhist : going to define mode ")

       ncerr = nf90_def_dim(ncid=ncid_hdr,name="dim2inout",len=2,dimid=dim2inout_id)
       call handle_ncerr(ncerr," outxfhist : define dim2inout")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="mxfh",len=mxfh,dimid=mxfh_id)
       call handle_ncerr(ncerr," outxfhist : define mxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="nxfh",len=nxfh,dimid=nxfh_id)
       call handle_ncerr(ncerr," outxfhist : define nxfh")
       ncerr = nf90_def_dim(ncid=ncid_hdr,name="xfdim2",len=xfdim2,dimid=xfdim2_id)
       call handle_ncerr(ncerr," outxfhist : define xfdim2")

       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dimr3",dimid=dimr3_id)
       call handle_ncerr(ncerr," outxfhist : inquire dimr3")

! xfhist(3,natom+4,2,mxfh)
       ncerr = nf90_def_var(ncid=ncid_hdr,name="xfhist",xtype=NF90_DOUBLE,&
          & dimids=(/dimr3_id,xfdim2_id,dim2inout_id,mxfh_id/),varid=xfhist_id)
       call handle_ncerr(ncerr," outxfhist : define xfhist")

! End define mode and go to data mode
       ncerr = nf90_enddef(ncid=ncid_hdr)
       call handle_ncerr(ncerr," outxfhist : enddef call ")
     else
! check that the dimensions are correct
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
       call handle_ncerr(ncerr," outxfhist : inquire nxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
            &   len=nxfh_tmp)
       call handle_ncerr(ncerr,"  outxfhist : get nxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="xfdim2",dimid=xfdim2_id)
       call handle_ncerr(ncerr," outxfhist : inquire xfdim2")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=xfdim2_id,&
            &   len=xfdim2_tmp)
       call handle_ncerr(ncerr,"  outxfhist : get xfdim2")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="mxfh",dimid=mxfh_id)
       call handle_ncerr(ncerr," outxfhist : inquire mxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=mxfh_id,&
            &   len=mxfh_tmp)
       call handle_ncerr(ncerr,"  outxfhist : get mxfh")
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="dim2inout",dimid=dim2inout_id)
       call handle_ncerr(ncerr," outxfhist : inquire dim2inout")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=dim2inout_id,&
            &   len=dim2inout_tmp)
       call handle_ncerr(ncerr,"  outxfhist : get dim2inout")

       ncerr = nf90_inq_varid(ncid=ncid_hdr,name="xfhist",varid=xfhist_id)
       call handle_ncerr(ncerr," outxfhist : inquire xfhist")

       if (mxfh_tmp /= mxfh .or. dim2inout_tmp /= 2 .or. xfdim2_tmp /= xfdim2) then
        write (*,*) 'outxfhist : ERROR xfhist has bad dimensions in NetCDF file. Can not re-write it.'
        stop
       end if

     end if

! Now fill the data
     ncerr = nf90_put_var(ncid=ncid_hdr,varid=xfhist_id,values=xfhist)
     call handle_ncerr(ncerr," outxfhist : fill xfhist")

! end NETCDF definition ifdef
#endif
    end if
#else
              call xderiveWRecInit(wff2,ierr)
              call xderiveWriteVal(wff2,nxfh)
              call xderiveWRecEnd(wff2,ierr)
    do ixfh=1,nxfh
              call xderiveWRecInit(wff2,ierr)

!              call xderiveWrite(wff2 &
!              & ,reshape(xfhist(1:3,1:natom+4,1:2,ixfh:ixfh)  &
!              & ,(/3*(natom+4)*2/)) &
!              &         ,3*(natom+4)*2,ierr)

              call xderiveWRecEnd(wff2,ierr)
    end do
#endif

 else if ( option == 2 ) then

#if !defined MPI
  if (wff2%accesswff /= 2) then
   read(unit=wff2%unwff,iostat=ios)nxfh

#if defined HAVE_NETCDF
  else if (wff2%accesswff == 2) then
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
       call handle_ncerr(ncerr," outxfhist : inquire nxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
            &   len=nxfh)
       call handle_ncerr(ncerr,"  outxfhist : get nxfh")
#endif
  end if
#else
              call xderiveRRecInit(wff2,ierr)
              call xderiveReadVal(wff2,nxfh)
              call xderiveRRecEnd(wff2,ierr)
#endif

 else if ( option == 3 ) then
    do ixfh=1,nxfh

#if !defined MPI
  if (wff2%accesswff /= 2) then
     read(unit=wff2%unwff,iostat=ios)xfhist(:,:,:,ixfh)
  end if

#else
              call xderiveRRecInit(wff2,ierr)
!              call xderiveRead(wff2 &
!              & ,reshape(xfhist(1:3,1:natom+4,1:2,ixfh:ixfh)  &
!              & ,(/3*(natom+4)*2/)) &
!              &         ,3*(natom+4)*2,ierr)
              call xderiveRRecEnd(wff2,ierr)
#endif
    end do

#if defined HAVE_NETCDF
    if (wff2%accesswff == 2) then
       ncerr = nf90_inq_dimid(ncid=ncid_hdr,name="nxfh",dimid=nxfh_id)
       call handle_ncerr(ncerr," outxfhist : inquire nxfh")
       ncerr = nf90_Inquire_Dimension(ncid=ncid_hdr,dimid=nxfh_id,&
            &   len=nxfh)
       call handle_ncerr(ncerr,"  outxfhist : get nxfh")

       ncerr = nf90_inq_varid(ncid=ncid_hdr,varid=xfhist_id,name="xfhist")
       call handle_ncerr(ncerr," outxfhist : inquire xfhist")
       ncerr = nf90_get_var(ncid=ncid_hdr,varid=xfhist_id,values=xfhist,&
         & start=(/1,1,1,1/),count=(/3,natom+4,2,nxfh/))
       call handle_ncerr(ncerr," outxfhist : read xfhist")
    end if
#endif

 else
  write(6,*)' outxfhist : option ', option , ' not available '

 end if


!DEBUG
! write(6,*)' outxfhist : exit'
!ENDDEBUG

end subroutine outxfhist
!!***
