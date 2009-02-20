!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffOpen
!! NAME
!! WffOpen
!!
!! FUNCTION
!! This subroutine opens a Wf file. It might be accessed
!! by different mechanisms (usual F90 IO routines,
!!  MPI I/O, or, in the future, NetCDF). The routine
!! provides a file handler, wff (a data structure containing
!! all needed information).
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MB,MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! accesswff=access mode (0 means all procs access using usual F90
!!  routines ; -1 means only the master proc access, using usual
!!  F90 routines ; 1 means MPI I/O; 2 means netcdf I/O)
!! filename=name of the file
!! master=the number of the master proc (only needed in parallel)
!! me=my number (only needed in parallel)
!! spaceComm= the space communicator handler (only needed in parallel)
!! unwff=the file unit number
!!
!! OUTPUT
!! ier=error code
!! wff= structured info about the wavefunction file
!!
!! PARENTS
!!      conducti_nc,conducti_paw,cut3d,gstate,inwffil,inwffil3,ioarr
!!      linear_optics_paw,loper3,nstdy3,optic,outwf,read_el_veloc,uderiv
!!      wannier
!!
!! CHILDREN
!!      etsf_io_low_error_to_str,etsf_io_low_open_modify,handle_ncerr,leave_new
!!      mpi_file_open,mpi_type_size,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffOpen(accesswff,spaceComm,filename,ier,wff,master,me,unwff)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif
#if defined HAVE_ETSF_IO
 use etsf_io
#endif

#if defined MPI && defined MPI2 && defined MPI_IO
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi, except_this_one => WffOpen
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1 && defined MPI_IO
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer, intent(in)  :: accesswff,spaceComm,master,me,unwff
 integer, intent(out) :: ier
 character(len=fnlen), intent(in) :: filename
 type(wffile_type), intent(out) :: wff

!Local variables-------------------------------
 character(len=500) :: message
!no_abirules
#if defined HAVE_NETCDF
           integer :: ncerr, ndim,nvar,natt,uid
#endif
#if defined HAVE_ETSF_IO
  type(etsf_io_low_error) :: error
  logical                 :: lstat
  character(len = etsf_io_low_error_len)   :: errmess
#endif

! *************************************************************************

!Initialize the mandatory data of the wff datastructure
 wff%unwff    =unwff
 wff%accesswff=accesswff
 wff%fname =filename

!Initialize info useful for parallel use
 wff%me       =me
 wff%master   =master
 wff%spaceComm=spaceComm

 ier=0
 if(accesswff==0)then

! All processors see a local file
  open (unit=unwff,file=filename,form='unformatted')
  rewind(unwff)
 else if(accesswff==-1)then

! Only the master processor see a local file
  if(master==me)then
   open (unit=unwff,file=filename,form='unformatted')
   rewind(unwff)
  end if

#if defined MPI_IO
           else if(accesswff==1)then
!           In the parallel case, only the master open filename file
            if(master==me)then
             open (unit=unwff,file=filename,form='unformatted')
             rewind(unwff)
            end if
!           In the parallel case, all procs open filename file
            call MPI_FILE_OPEN(spaceComm,filename,&
&            MPI_MODE_CREATE + MPI_MODE_RDWR,&
&            MPI_INFO_NULL,wff%fhwff, ier)
!           Define all type values
            call MPI_Type_size(MPI_INTEGER,wff%nbOct_int,ier)
            call MPI_Type_size(MPI_DOUBLE_PRECISION,wff%nbOct_dp,ier)
            call MPI_Type_size(MPI_CHARACTER,wff%nbOct_ch,ier)
            wff%offwff=0
            wff%kgwff=-1
            wff%formwff=-1
#endif
#if defined HAVE_NETCDF
 else if (accesswff==2)then
! here we modify wff%unwff to save the netCDF identifier in it
   ncerr = nf90_open(path=filename, mode=NF90_WRITE, ncid=wff%unwff)
   call handle_ncerr(ncerr," open netcdf wavefunction file")
   write (*,*) ' WffOpen : open a netCDF file ', trim(filename), wff%unwff
   !DEBUG
   ncerr = nf90_Inquire(ncid=wff%unwff,nDimensions=ndim,nVariables=nvar,nAttributes=natt,unlimitedDimId=uid)
   call handle_ncerr(ncerr, " general Inquire ")
   write (*,*) 'WffOpen : found ndim,nvar,natt,uid = ', ndim,nvar,natt,uid
   !ENDDEBUG
#endif
#if defined HAVE_ETSF_IO
 else if (accesswff==3)then
   call etsf_io_low_open_modify(wff%unwff, trim(filename)//"-etsf.nc", lstat, error_data = error)
   if (.not. lstat) then
      call etsf_io_low_error_to_str(errmess, error)
      write(message, "(A,A,A,A)") ch10, " WffOpen: ERROR -", ch10, &
                                & errmess(1:min(475, len(errmess)))
      call wrtout(std_out, message, 'COLL')
      call leave_new('COLL')
   end if
   write(message, '(3A,I0)' ) ' WffOpen: opening ', trim(wff%fname), &
                            & "-etsf.nc on unit ", wff%unwff
   call wrtout(std_out, message, 'COLL')
#endif
 else

  write(message, '(12a,i6,3a)' ) ch10, &
&  ' WffOpen : ERROR -',ch10,&
&  '  For the time being the input variable accesswff is restricted ',ch10,&
&  '  to 0 (all cases), 1 (in case MPI is enabled),',ch10,&
&  '  2 (only sequential, and if the NetCDF library has been enabled),',ch10,&
&  '  or 3 (only sequential, and if the NetCDF and ETSF_IO libraries have been enabled).',ch10,&
&  '  Its value is accesswff=',accesswff,'.',ch10,&
&  '  Action : change accesswff or use ABINIT in parallel or enable NetCDF and/or ETSF_IO.'
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')
  call leave_new('COLL')

 end if

end subroutine WffOpen
!!***
