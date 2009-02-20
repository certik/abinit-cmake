!{\src2tex{textfont=tt}}
!!****f* ABINIT/split_work
!! NAME
!! split_work
!!
!! FUNCTION
!!  Split a number of tasks, ntasks, among nprocs processors.
!!  Useful for on-the-fly parallelization of simple loops
!!  Note that if nprocs>ntasks then : 
!!   istart=ntasks+1
!!   istop=ntask 
!!
!!  In this particular case, loops of the form  
!!
!!  do ii=istart,istop 
!!   ...
!!  end do
!! 
!!  are not executed. Moreover allocation such as foo(istart:istop) 
!!  will generate a zero-sized array and 
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ntasks=number of tasks
!!
!! OUTPUT
!!  istart,istop= indeces defining the initial and final task for 
!!   this processor
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  
!!
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine split_work(ntasks,istart,istop,verbose)
    
 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => split_work
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks
 integer,optional,intent(in) :: verbose
 integer,intent(inout) :: istart,istop

!Local variables-------------------------------
 integer ierr,res,nprocs,rank,block,block_tmp,verb
 character(len=500) :: msg                  
 
! *************************************************************************
 
 verb=0 ; if (PRESENT(verbose)) verb=verbose
 
 call xproc_max(nprocs,ierr)
 call xme_whoiam(rank)

 block_tmp=ntasks/nprocs 
 res=MOD(ntasks,nprocs) 
 block=block_tmp+1   

 if (res/=0.and.verb/=0) then 
  write(msg,'(4a,i5,a,i4)')ch10,&
&  ' split_work : WARNING - ',ch10,&
&  '  number of tasks= ',ntasks,' not divisible by nprocs= ',nprocs
  call wrtout(std_out,msg,'COLL')
 end if 
 if (block_tmp==0) then 
  write(msg,'(4a,i4,a,i5,3a)')ch10,&
&  ' split_work : WARNING - ',ch10,&
&  ' number of processors= ',nprocs,' larger than number of tasks= ',ntasks,ch10,&
&  ' This is a waste ',ch10
  call wrtout(std_out,msg,'COLL')
 end if 

 if (rank<res) then
  istart= rank   *block+1
  istop =(rank+1)*block
 else
  istart=res*block+(rank-res  )*block_tmp+1
  istop =res*block+(rank-res+1)*block_tmp
 end if

end subroutine split_work
!!***
