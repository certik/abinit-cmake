!{\src2tex{textfont=tt}}
!!****f* ABINIT/split_work2
!! NAME
!! split_work2
!!
!! FUNCTION
!!  Split a number of tasks, ntasks, among nprocs processors.
!!  The output arrays istart(1:nprocs) and istop(1:nprocs) 
!!  report the starting and final task index for each CPU.
!!  Namely CPU with rank ii has to perform all the tasks between 
!!  istart(ii+1) and istop(ii+1). Note the Fortran convention of using 
!!  1 as first index of the array.
!!  Note, moreover, that if a proc has rank>ntasks then : 
!!   istart(rank+1)=ntasks+1
!!   istop(rank+1)=ntask 
!!
!!  In this particular case, loops of the form  
!!
!!  do ii=istart(rank),istop(rank) 
!!   ...
!!  end do
!! 
!!  are not executed. Moreover allocation such as foo(istart(rank):istop(rank)) 
!!  will generate a zero-sized array
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
!!  istart(nprocs),istop(nprocs)= indeces defining the initial and final task for 
!!   each processor
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

subroutine split_work2(ntasks,nprocs,istart,istop,verbose)
    
 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi, except_this_one => split_work2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks,nprocs
 integer,optional,intent(in) :: verbose
 integer,intent(inout) :: istart(nprocs),istop(nprocs)

!Local variables-------------------------------
 integer ierr,res,irank,block,block_tmp,verb
 character(len=500) :: msg                  
 
! *************************************************************************
 
 verb=0 ; if (PRESENT(verbose)) verb=verbose
 
 block_tmp=ntasks/nprocs 
 res=MOD(ntasks,nprocs) 
 block=block_tmp+1   

 if (res/=0.and.verb/=0) then 
  write(msg,'(4a,i5,a,i4,3a)')ch10,&
&  ' split_work : WARNING : ',ch10,&
&  ' number of tasks = ',ntasks,' is not divisible by nprocs = ',nprocs,ch10,&
&  ' parallelism is not efficient ',ch10
  call wrtout(std_out,msg,'COLL')
 end if 
 if (block_tmp==0) then 
  write(msg,'(4a,i4,a,i5,3a)')ch10,&
&  ' split_work : WARNING : ',ch10,&
&  ' number of processors = ',nprocs,' is larger than number of tasks =',ntasks,ch10,&
&  ' This is a waste ',ch10
  call wrtout(std_out,msg,'COLL')
 end if 

 do irank=0,nprocs-1
  if (irank<res) then
   istart(irank+1)= irank   *block+1
   istop (irank+1)=(irank+1)*block
  else
   istart(irank+1)=res*block+(irank-res  )*block_tmp+1
   istop (irank+1)=res*block+(irank-res+1)*block_tmp
  end if
 end do

end subroutine split_work2
!!***
