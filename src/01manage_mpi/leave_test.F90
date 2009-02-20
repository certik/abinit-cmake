!{\src2tex{textfont=tt}}
!!****f* ABINIT/leave_test
!! NAME
!! leave_test
!!
!! FUNCTION
!! Routine that tests whether exit must be done,
!! because of eventual problems encountered by another processor.
!! In this case, will make a clean exit.
!! In the sequential case, return.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GMR, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (no input)
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abinit,ctocprj,dyfnl3,eltfrkin3,eltfrnl3,energy,forstrnps,getgsc
!!      initylmg,inwffil3,iofn1,iofn2,kpgio,ladielmt,lavnl,loper3,memana,mkrho
!!      mkrho3,newkpt,nselt3,nstdy3,outkss,outwf,pawmkrhoij,prctfvw1,prctfvw2
!!      rhofermi3,scfcv,scfcv3,suscep_dyn,suscep_kxc_dyn,suscep_stat,vtorho
!!      vtorho3,wfsinp
!!
!! CHILDREN
!!      leave_myproc,mpi_allreduce,mpi_barrier,mpi_finalize,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!BEGIN TF_CHANGES
subroutine leave_test(mpi_enreg)

 use defs_basis
 use defs_datatypes

#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi, except_this_one => leave_test
!End of the abilint section

 implicit none

#if defined MPI && defined MPI1
 include 'mpif.h'
#endif
!Arguments ------------------------------------
 type(MPI_type) :: mpi_enreg
!END TF_CHANGES

!Local variables-------------------------------
!no_abirules
#if defined MPI
           integer :: gl_check_bit,ierr,my_check_bit,myproc
           real(dp) :: tsec(2)
           character(len=500) :: message
#endif

! **********************************************************************

#if defined MPI
           call timab(48,1,tsec)
          !Synchronize
           call MPI_BARRIER(mpi_enreg%spaceComm,ierr)
           call timab(48,2,tsec)
           write(message, '(a)' ) ' leave_test : synchronization done...'
           call wrtout(06,message,'PERS')

          !Everything is allright for me
           my_check_bit=0
           call timab(48,1,tsec)
          !See what about the others
           call MPI_ALLREDUCE(my_check_bit,gl_check_bit,1,MPI_INTEGER,&
          &  MPI_SUM,mpi_enreg%spaceComm,ierr)
           call timab(48,2,tsec)
          !Check for exit
           if(gl_check_bit>0) then
            write(message, '(a)' ) ' leave_test : exiting...'
            call wrtout(06,message,'PERS')
            call MPI_FINALIZE(ierr)
            call leave_myproc
           end if
#endif

end subroutine leave_test
!!***
