!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkexi
!! NAME chkexi
!! chkexi
!!
!!
!! FUNCTION
!! This routine checks whether the CPU time limit is exceeded or not.
!! If openexit is non-zero, it also checks the "filnam" file
!! for the "exit" character string in its first line and returns the location
!! of the string on the line (0 if not found).  Maps both strings to upper case
!! before attempting to match them. Also checks for the existence
!! of the "abinit.exit" file in the directory where the job was started.
!! Finally, checks whether the CPU time limit was not exceeded.
!! If one of these conditions occurs, will induce graceful exit of iterations.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cpus = CPU time limit
!!  filnam = character string giving name of file to be opened
!!  iout = unit number to print output to
!!  openexit = if 1, open the "filnam" and "abinit.exit" files
!!
!! OUTPUT
!!  iexit = index of "exit" on first line of file (0 if not found),
!!      or -1 if the exit was ordered through the existence of the "exit" file
!       or -2 if the exit was ordered through the CPU time limit.
!!
!! SIDE EFFECTS
!!  mpi_enreg <type(MPI_type)> = informations about MPI parallelization
!!
!! NOTES
!!
!!
!! PARENTS
!!      brdmin,cgwf,delocint,driver,drivergw,gstate,loper3,moldyn,move,respfn
!!      scprqt,vtowfk3
!!
!! CHILDREN
!!      inupper,leave_new,mpi_bcast,timein,wrtout,xcomm_init,xme_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine chkexi(cpus,filnam,iexit,iout,mpi_enreg,openexit)

 use defs_basis
 use defs_datatypes
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 real(dp)         ,intent(in)    :: cpus
 character(len=fnlen),intent(in)    :: filnam
 integer          ,intent(in)    :: openexit,iout
 integer          ,intent(out)   :: iexit
 type(MPI_type)   ,intent(inout) :: mpi_enreg

!Local variables-------------------------------
 integer :: ierr
 integer,save :: iexit_save=0
 real(dp) :: tsec(2)
 real(dp),save :: tcpu_last=zero
 logical :: ex
 character(len=500) :: message
 character(len=fnlen) :: line
 character(len=4), parameter :: string='EXIT'
!no_abirules
#if defined MPI 
       !Variables introduced for MPI version
        integer :: ierrmpi,me
        logical :: master
!BEGIN TF_CHANGES
        integer :: spaceComm
!END TF_CHANGES
#endif

! *************************************************************************

!DEBUG
!write(6,*)' chkexi : enter '
!ENDDEBUG

 if(iexit_save==0)then   ! ABINIT will pass again in this routine even after exit call has been detected

#if defined MPI 
          !Determine who I am
!BEGIN TF_CHANGES
        call xme_init(mpi_enreg,me)
        !Init spaceworld
        call xcomm_init(mpi_enreg,spaceComm)

        master=((mpi_enreg%parareel==0 .and. me==0) .or. (mpi_enreg%parareel==1 .and. &
&        mpi_enreg%master_group_para==me))
!END TF_CHANGES
           if (master) then
          !Proc 0 tests and broadcast the result to others
#endif

  iexit=0

! Is it worth to test the cpu time ?
  if(abs(cpus)>1.0d-5 .or. openexit==1)then
   call timein(tsec(1),tsec(2))
  end if

! A first way of exiting : the cpu time limit
  if(abs(cpus)>1.0d-5)then
   if(cpus<tsec(1))iexit=-2
!  DEBUG
!  write(6,*)' chkexi : cpus,tsec(1)',cpus,tsec(1)
!  ENDDEBUG
  end if

! Test the content of files only when sufficient time (2 sec) has elapsed from
! last time it was tested.

! DEBUG
! write(6,*)' chkexi : openexit,iexit,tsec(1)-tcpu_last', openexit,iexit,tsec(1)-tcpu_last
! ENDDEBUG

  if(openexit==1 .and. iexit==0 .and. tsec(1)-tcpu_last>two )then

   tcpu_last=tsec(1)

!  Open file and read first line as character string
   open(unit=tmp_unit,file=filnam,form='formatted',status='old')
   rewind (unit=tmp_unit)
   read (unit=tmp_unit,fmt='(a)',iostat=ierr) line
   if(ierr/=0)then
    write(message, '(a,a,a,a,a,a,i5,a,a)' ) ch10,&
&    ' chkexi: ERROR -',ch10,&
&    '  Problem when reading file=',filnam,&
&    '   iostat =',ierr,ch10,&
&    '  Action : check whether this file is OK.'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if
!  Make a local copy of matching string of length equal to
!  nonblank length of input string
!  Map to upper case
   call inupper(line)
   iexit=index(line,string)
   close (unit=tmp_unit)

!  This is another way of exiting : the previous one does not work
!  on some machines, may be because they keep a copy of the initial input file.
   if(iexit==0)then
    inquire(file='abinit.exit',exist=ex)
    if(ex)iexit=-1
   end if

  end if

#if defined MPI  
        end if
        if (mpi_enreg%parareel == 0) then
!BEGIN TF_CHANGES
                call MPI_BCAST(iexit,1,MPI_INTEGER,0,&
      &  spaceComm,ierrmpi)
!END TF_CHANGES
        else
                call MPI_BCAST(iexit,1,MPI_INTEGER,0,&
      &  mpi_enreg%kpt_comm_para(mpi_enreg%ipara),ierrmpi)
        end if
#endif

 else ! In case the exit mechanism has already been activated

  iexit=iexit_save

 end if

 if(iexit/=0)then
  if(iexit>0) write(message, '(a,a,a,a,a,a,a)' ) ch10,&
&  ' chkexi: WARNING -',ch10,&
&  '  Exit has been requested from file ',trim(filnam),'.',ch10
  if(iexit==-1) write(message, '(a,a,a,a,a)' ) ch10,&
&  ' chkexi: WARNING -',ch10,&
&  '  Exit has been requested from file "abinit.exit".',ch10
  if(iexit==-2) write(message, '(a,a,a,a,a)' ) ch10,&
&  ' chkexi: WARNING -',ch10,&
&  '  Exit due to cpu time limit exceeded.',ch10
  if(iout/=6)then
   call wrtout(iout,message,'COLL')
  end if
  call wrtout(06,  message,'COLL')
 end if

 iexit_save=iexit

!DEBUG
!write(6,*)'chkexi: exit'
!stop
!ENDDEBUG

end subroutine chkexi
!!***
