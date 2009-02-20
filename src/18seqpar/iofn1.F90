!{\src2tex{textfont=tt}}
!!****f* ABINIT/iofn1
!! NAME
!! iofn1
!!
!! FUNCTION
!! Begin by eventual redefinition of unit 05 and 06
!! Then, print greetings for interactive user.
!! Next, Read filenames from unit 05, AND check that new
!! output file does not already exist.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!  character(len=fnlen) :: filnam(5)=character strings giving file names
!!  character(len=fnlen) :: filstat=character strings giving name of status file
!!
!! NOTES
!! If it does exist, isfile will create a new name
!! to avoid overwriting the output file.
!! Also create name of status file
!!
!! File names refer to following files, in order:
!!  (1) Formatted input file  (05)
!!  (2) Formatted output file (iout)
!!  (3) Root name for generic input files (wavefunctions, potential, density ...)
!!  (4) Root name for generic output files (wavefunctions, potential, density,
!!                                          DOS, hessian ...)
!!  (5) Root name for generic temporary files (wftmp1,wftmp2,kgunit,status ...)
!!
!! PARENTS
!!      abinit,iofn1_ncdf
!!
!! CHILDREN
!!      int2char4,isfile,leave_new,leave_test,mpi_bcast,timab,wrtout,xme_whoiam
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!BEGIN TF_CHANGES
subroutine iofn1(filnam,filstat,mpi_enreg)
!END TF_CHANGES

 use defs_basis
!BEGIN TF_CHANGES
 use defs_datatypes
!END TF_CHANGES

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
 character(len=fnlen), intent(out) :: filstat
 character(len=fnlen), intent(out) :: filnam(5)
 type(MPI_type), intent(in) :: mpi_enreg

!Local variables-------------------------------
 character(len=1) :: blank=' '
 integer :: ierr,me
 logical :: test_mpi
 character(len=fnlen) :: fillog,tmpfil
 character(len=4) :: tag
 character(len=500) :: message
 real(dp) :: tsec(2)

!*************************************************************************

! Determine who I am
 call xme_whoiam(me)

 if(me==0) then

! Eventually redefine standard input and standard output
#if defined READ_FROM_FILE

! First take care of the output file
  tmpfil(1:fnlen)=blank
  tmpfil(1:3)='log'
  call isfile(tmpfil,'new')
  close(6)
  open (unit=6,file=tmpfil,form='formatted',status='new')

! Now take care of the "files" file
  tmpfil(1:fnlen)=blank
  tmpfil(1:9)='ab.files'
  write(message, '(a,a,a,a,a,a,a)' ) ch10,&
&   ' iofn1 : COMMENT -',ch10,&
&   '  Because of cpp option READ_FROM_FILE,',ch10,&
&   '  read file "ab.files" instead of standard input ' ,ch10
  call wrtout(6,message,'COLL')
  call isfile(tmpfil,'old')
  close(5)
  open (unit=5,file=tmpfil,form='formatted',status='old')

#endif

! Print greetings for interactive user
  write(06,*)' ABINIT '
  write(06,*)' '

! Read name of input file (05):
  write(06,*)' Give name for formatted input file: '
  read(05, '(a)' ) filnam(1)
  write(06, '(a)' ) trim(filnam(1))
  write(06,*)' Give name for formatted output file:'
  read (05, '(a)' ) filnam(2)
  write (06, '(a)' ) trim(filnam(2))
  write(06,*)' Give root name for generic input files:'
  read (05, '(a)' ) filnam(3)
  write (06, '(a)' ) trim(filnam(3))
  write(06,*)' Give root name for generic output files:'
  read (05, '(a)' ) filnam(4)
  write (06, '(a)' ) trim(filnam(4))
  write(06,*)' Give root name for generic temporary files:'
  read (05, '(a)' ) filnam(5)
  write (06, '(a)' ) trim(filnam(5))

! Check that old input file exists
  call isfile(filnam(1),'old')

! Check that new output file does NOT exist
  call isfile(filnam(2),'new')

! Check that root name for generic input and output differ
  if ( trim(filnam(3))==trim(filnam(4)) ) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' iofn1 : ERROR -',ch10,&
&   '  Root name for generic input and output files must differ ',ch10,&
&   '  Action : correct your "file" file.'
   call wrtout(6,message,'COLL')
   call leave_new('PERS')
  end if

! Check that root names are at least 20 characters less than fnlen
  if ( len_trim(filnam(3)) >= (fnlen-20) ) then
   write(message, '(a,a,a,a,a,a,a,a,i4,a,i4,a,a)' ) ch10,&
&   ' iofn1 : ERROR -',ch10,&
&   '  Root name for generic input files is too long. ',ch10,&
&   '  It must be 20 characters less than the maximal allowed ',ch10,&
&   '  length of names, that is ',fnlen,', while it is ',len_trim(filnam(3)),&
&   ch10,'  Action : correct your "file" file.'
   call wrtout(6,message,'COLL')
   call leave_new('PERS')
  end if
  if ( len_trim(filnam(4)) >= (fnlen-20) ) then
   write(message, '(a,a,a,a,a,a,a,a,i4,a,i4,a,a)' ) ch10,&
&   ' iofn1 : ERROR -',ch10,&
&   '  Root name for generic output files is too long. ',ch10,&
&   '  It must be 20 characters less than the maximal allowed ',ch10,&
&   '  length of names, that is ',fnlen,', while it is ',len_trim(filnam(4)),&
&   ch10,'  Action : correct your "file" file.'
   call wrtout(6,message,'COLL')
   call leave_new('PERS')
  end if
  if ( len_trim(filnam(5)) >= (fnlen-20) ) then
   write(message, '(a,a,a,a,a,a,a,a,i4,a,i4,a,a)' ) ch10,&
&   ' iofn1 : ERROR -',ch10,&
&   '  Root name for generic temporary files is too long. ',ch10,&
&   '  It must be 20 characters less than the maximal allowed ',ch10,&
&   '  length of names, that is ',fnlen,', while it is ',len_trim(filnam(5)),&
&    ch10,'  Action : correct your "file" file.'
   call wrtout(6,message,'COLL')
   call leave_new('PERS')
  end if

!End the section me==0
 end if

#if defined MPI 
           call timab(48,1,tsec)
           call leave_test(mpi_enreg)
           call MPI_BCAST(filnam(1:5),5*fnlen,MPI_CHARACTER,0,&
          &  MPI_COMM_WORLD,ierr)
           call timab(48,2,tsec)
#endif

!Create a name for the status file, based on filnam(5)
 filstat=trim(filnam(5))//'_STATUS'

!Redefine the log unit if not the master
 if(me/=0)then
  call int2char4(me,tag)
  filstat=trim(filstat)//'_P-'//tag
  fillog=trim(filnam(5))//'_LOG_'//tag
  close(6)
  open(unit=6,file=fillog,status='unknown')
 end if

end subroutine iofn1
!!***
