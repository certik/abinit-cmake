!{\src2tex{textfont=tt}}
!!****f* ABINIT/symkpt
!! NAME
!! symkpt
!!
!! FUNCTION
!! Determines the weights of the k-points for sampling
!! the Brillouin Zone, starting from a first set
!! of weights wtk, and folding it to a new set, by
!! taking into account the symmetries described
!! by symrc1, and eventually the timereversal symmetry.
!! Also compute the number of k points in the reduced set
!! This routine is also used for sampling the q vectors in the
!! Brillouin zone for the computation of thermodynamical
!! properties (from the routine thm9).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG,LSI)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gmet(3,3)=reciprocal space metric (bohr**-2).
!! nkpt = number of k-points whose weights are wtk
!! kptns(3,nkpt)= k vectors in reciprocal space
!! nsym1=number of space group symmetries
!! option=if 1, output the new number of kpoints on unit ab_out
!! symrc1(3,3,nsym1)=3x3 matrices of the group symmetries (reciprocal space)
!! timrev: if 1, the time reversal operation has to be taken into account
!!         if 0, no time reversal symmetry.
!! wtk(nkpt)=weight assigned to each k point.
!!
!! OUTPUT
!! indkpt1(nkpt)=non-symmetrized indices of the k-points
!! nkpt1 = number of k-points in the irreducible set
!! wtk_folded(nkpt)=weight assigned to each k point, taking into account the
!!                  symmetries
!!
!! NOTES
!! The decomposition
!! of the symmetry group in its primitives might speed up the execution.
!! The output variables are stored only in the range 1:nkpt
!!
!! PARENTS
!!      elphon,getkgrid,loper3,mkphdos,setup_little_group,thm9
!!
!! CHILDREN
!!      leave_new,sort_dp,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symkpt(gmet,indkpt1,kptns,nkpt,nkpt1,nsym1,option,&
& symrc1,timrev,wtk,wtk_folded)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: nkpt,nsym1,option,timrev
 integer,intent(out) :: nkpt1
!arrays
 integer :: symrc1(3,3,nsym1)
 integer,intent(out) :: indkpt1(nkpt)
 real(dp),intent(in) :: gmet(3,3),kptns(3,nkpt),wtk(nkpt)
 real(dp),intent(out) :: wtk_folded(nkpt)

!Local variables -------------------------
!scalars
 integer :: identi,ii,ikpt,ikpt2,ind_ikpt,ind_ikpt2,isym,itim,jj,quit,tident
 real(dp) :: difk,length2trial,reduce,shift
 character(len=500) :: message
!arrays
 integer,allocatable :: list(:)
 real(dp) :: gmetkpt(3),ksym(3)
 real(dp),allocatable :: length2(:)

! *********************************************************************

!DEBUG
!write(6,*)' enter symkpt '
!write(6,*)' nkpt,nsym1',nkpt,nsym1
!write(6,*)' wtk',wtk
!write(6,*)' timrev',timrev
!if(option==1)stop
!ENDDEBUG

 if(timrev/=1 .and. timrev/=0)then
  write(message, '(a,a,a,a,a,i4,a)' )&
&  ' symkpt : BUG -',ch10,&
&  '  timrev should be 0 or 1, while',ch10,&
&  '  it is equal to ',timrev,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if(nsym1/=1)then
! Find the identity symmetry operation
  do isym=1,nsym1
   tident=1
   do jj=1,3
    if(symrc1(jj,jj,isym)/=1)tident=0
    do ii=1,3
     if( ii/=jj .and.&
&     symrc1(ii,jj,isym)/=0)tident=0
    end do
   end do
   if(tident==1)then
    identi=isym
    write(message, '(a,i3)' )' symkpt : found identity, with number',identi
    call wrtout(6,message,'COLL')
    exit
   end if
  end do
  if(tident==0)then
   write(message, '(a,a,a)' )&
&   ' symkpt : BUG -',ch10,&
&   '  Did not found the identity operation.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
 end if

!Initialise the wtk_folded array using the wtk array :
 do ikpt=1,nkpt
  wtk_folded(ikpt)=wtk(ikpt)
 end do

!Here begins the serious business

!If there is some possibility for a change (otherwise, wtk_folded is
!correctly initialized to give no change) :
 if(nkpt/=1 .and. (nsym1/=1 .or. timrev==1) )then

! Store the length of vectors, but take into account umklapp
! processes by selecting the smallest length of all symmetric vectors
  allocate(length2(nkpt))

! DEBUG
! write(6,*)' gmet ',gmet(3,3)
! write(6,*)' kpt '
! do ikpt=1,nkpt
! write(6,*)ikpt,kptns(:,ikpt)
! end do
! write(6,*)' symrc1 '
! do isym=1,nsym1
! write(6,*)isym,symrc1(:,:,isym)
! end do
! ENDDEBUG

  do ikpt=1,nkpt
   do isym=1,nsym1
    do itim=1,(1-2*timrev),-2
!    Get the symmetric of the vector
     do ii=1,3
      ksym(ii)=itim*( kptns(1,ikpt)*symrc1(ii,1,isym)&
&      +kptns(2,ikpt)*symrc1(ii,2,isym)&
&      +kptns(3,ikpt)*symrc1(ii,3,isym) )
      ksym(ii)=ksym(ii)-anint(ksym(ii)+tol8*half)
     end do
     gmetkpt(:)=gmet(:,1)*ksym(1)+gmet(:,2)*ksym(2)+gmet(:,3)*ksym(3)
     length2trial=ksym(1)*gmetkpt(1)+ksym(2)*gmetkpt(2)+&
&     ksym(3)*gmetkpt(3)
     if(isym==1 .and. itim==1)then
      length2(ikpt)=length2trial
     else
      if(length2(ikpt)>length2trial)length2(ikpt)=length2trial
     end if
    end do
   end do
  end do

! Sort the lengths
  allocate(list(nkpt))
  list(:)=(/ (ikpt,ikpt=1,nkpt) /)
  call sort_dp(nkpt,length2,list,tol14)

! DEBUG
! do ikpt=1,nkpt
! write(6,*)ikpt,length2(ikpt),list(ikpt)
! end do
! ENDDEBUG

  do ikpt=1,nkpt-1

!  Ordered index
   ind_ikpt=list(ikpt)

!  Not worth to examine a k point that is a symmetric of another,
!  which is the case if its weight has been set to 0 by previous folding
   if(wtk_folded(ind_ikpt)<tol16)cycle

!  Loop on the remaining k-points
   do ikpt2=ikpt+1,nkpt

!   The next line eliminates pairs of vectors that differs by their length.
!   Moreover, since the list is ordered according to the length,
!   one can skip all other ikpt2 vectors, as soon as one becomes larger
!   than length2(ikpt)
    if(length2(ikpt2)-length2(ikpt)>tol8)exit

!   Ordered index
    ind_ikpt2=list(ikpt2)

!   If the second vector is already empty, no interest to treat it
    if(wtk_folded(ind_ikpt2)<tol16)cycle

    quit=0
    do isym=1,nsym1
     do itim=1,(1-2*timrev),-2
      if(isym/=identi .or. itim/=1 )then

!      Get the symmetric of the vector
       do ii=1,3
        ksym(ii)=itim*( kptns(1,ind_ikpt)*symrc1(ii,1,isym)&
&        +kptns(2,ind_ikpt)*symrc1(ii,2,isym)&
&        +kptns(3,ind_ikpt)*symrc1(ii,3,isym) )
       end do

!      The do-loop was expanded to speed up the execution
       difk= ksym(1)-kptns(1,ind_ikpt2)
       reduce=difk-anint(difk)
       if(abs(reduce)>tol8)cycle
       difk= ksym(2)-kptns(2,ind_ikpt2)
       reduce=difk-anint(difk)
       if(abs(reduce)>tol8)cycle
       difk= ksym(3)-kptns(3,ind_ikpt2)
       reduce=difk-anint(difk)
       if(abs(reduce)>tol8)cycle

!      Here, have successfully found a symmetrical k-vector
!      Assign all the weight of the k-vector to its symmetrical
       wtk_folded(ind_ikpt)=wtk_folded(ind_ikpt)+wtk_folded(ind_ikpt2)
       wtk_folded(ind_ikpt2)=0._dp

!      Go to the next ikpt2 if the symmetric was found
       quit=1
       exit

!      End condition of non-identity symmetry
      end if

!     End loop on itim
     end do

     if(quit==1)exit
!    End loop on isym
    end do

!   End secondary loop over k-points
   end do

!  End primary loop over k-points
  end do

  deallocate(length2,list)

! End check on possibility of change
 end if

!Create the indexing array indkpt1
 nkpt1=0
 do ikpt=1,nkpt
  if(wtk_folded(ikpt)>tol8)then
   nkpt1=nkpt1+1
   indkpt1(nkpt1)=ikpt
  end if
 end do

 if(option/=0)then
  if(nkpt/=nkpt1)then
   write(message, '(a,a,a,i6,a)' )&
&   ' symkpt : the number of k-points, thanks to the symmetries,',ch10,&
&   ' is reduced to',nkpt1,' .'
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')
!  DEBUG
!  write(message, '(a)' )'   Here are the new weights :'
!  call wrtout(6,message,'COLL')
!  do ikpt=1,nkpt,6
!  write(message, '(6f12.6)' ) wtk_folded(ikpt:min(nkpt,ikpt+5))
!  call wrtout(6,message,'COLL')
!  end do
!  ENDDEBUG
  else
   write(message, '(a)' )&
&   ' symkpt : not enough symmetry to change the number of k points.'
   call wrtout(6,message,'COLL')
   call wrtout(ab_out,message,'COLL')
  end if
 end if

!DEBUG
!write(6,*)' exit symkpt '
!write(6,*)' nkpt,nsym1',nkpt,nsym1
!write(6,*)' wtk',wtk
!write(6,*)' timrev',timrev
!if(timrev==0)stop
!if(option==1)stop
!ENDDEBUG

end subroutine symkpt
!!***
