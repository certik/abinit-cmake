!{\src2tex{textfont=tt}}
!!****f* ABINIT/findggp
!! NAME
!! findggp
!!
!! FUNCTION
!!  Fing the independent (G1,G2) pairs that are sufficient to reconstruct using 
!!  symmetry properties the Fourier components f_q(G1,G2) of a two-point function 
!!  which has the same symmetry of the crystal 
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Gsphere<Gvectors_type>=Info on the G-sphere
!!   %ng=number of G vectors in the f matrix
!!   %gmet(3,3)=metric tensor in reciprocal space
!!   %gvec(3,ng)=G vectors in reduced coordinates
!!  nsym=number of symmetry operations
!!  qpt(3)=q-point in reciprocal coordinated
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!
!! OUTPUT
!!  Gpairs_q<Gpairs_type>= Structure containing information on the irreducible pairs
!!   %niggp=nuber of independent (G1,G1) pairs   
!!   %fp2ip(2,ng,ng)= for given T1 and T1 reports the sequential index of 
!!     the indendent pair (G1,G2) such as (T1,T2)= S (G1,G2) 
!!   %fptabo(ng,ng)= index of the symmetry operation S in the array symrec such as (T1,T2)= S (G1,G2)
!!   %ip2fp(2,niggp)= index of G1 and G2 in the array gvec for each niggp independent pair. 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! TODO
!!  Write an interface for vnorm and dosym (also vdotw) in order to deal with integer vectors 
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

subroutine findggp(nsym,symrec,Gsphere,qpt,Gpairs_q)
    
 use defs_basis
 use m_gwdefs, only : GW_TOLQ
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13recipspace
 use interfaces_15gw, except_this_one => findggp
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 type(Gvectors_type),intent(in) :: Gsphere
 type(Gpairs_type),intent(inout) :: Gpairs_q 
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 real(dp),intent(in) :: qpt(3) 

!Local variables-------------------------------
!scalars
 integer :: igpi1,igpi2,isym,itim,ig,ig1,ig2,igpi,ng,niggp
 integer :: dummy,timrev_,iid,ic,istat
 real(dp) :: diff2,norm1,norm2,eps,eps1,eps2
 real(dp) :: a1,a2,a3,d1,d2,d3,g12,g23,g31
 logical :: found,found_identity
 character(len=500) :: msg           
!arrays
 integer :: identity(3,3),symqpt(4,2,nsym),ltg(nsym)
 integer :: gposs1(3),gposs2(3),gxx1(3),gxx2(3)
 integer, allocatable :: iperm(:)
 real(dp) :: gmet(3,3)
 real(dp),allocatable :: gnorm(:),gdiff2(:)
 integer,allocatable :: ip2fp(:,:)
 integer,pointer :: gvec(:,:)
 character(len=50),parameter :: my_name='findggp'

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(1x,a)')TRIM(my_name)//': enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 !
 ! === Check presence of identity as save its index ===
 identity(:,:)=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
 found_identity=.FALSE.
 do isym=1,nsym
  if (ALL(symrec(:,:,isym)==identity)) then 
   iid=isym
   found_identity=.TRUE. ; EXIT
  end if
 end do 
 if (.not.found_identity) then 
  write(msg,'(4a)')ch10,&
&  ' findggp : ERROR -',&
&  '  Only the inversion was found in the set of symmetries read from the KSS file ',ch10,&
&  '  Likely old version of KSS file! '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 else 
  write(msg,'(a,i2)')' findggp : found identity with index: ',iid
  call wrtout(std_out,msg,'COLL')
 end if 
 !
 ! === Only at Gamma time-reversal can be included ===
 timrev_=1 ; if (ALL(ABS(qpt)<GW_TOLQ).and.Gsphere%timrev==2) timrev_=2

 write(msg,'(2a,3f8.5)')ch10,&
& ' Analyzing symmetries at point : ',qpt(:)
 if (timrev_==2) write(msg,'(a)')TRIM(msg)//' (including time-reversal) '
 call wrtout(std_out,msg,'COLL')
 !
 ! === Find operations in the little group === 
 ! ltg is 1 if the symmetry  preserves q with a zero umklapp vector
 call symq3(nsym,qpt,symqpt,symrec,dummy)

 ltg(:)=0 
 do itim=1,timrev_
  do isym=1,nsym
   if (symqpt(4,itim,isym)==1.and.ALL(symqpt(1:3,itim,isym)==0)) ltg(isym)=1
  end do 
 end do

 if (SUM(ltg)==1.and.timrev_==1) then 
  ! === In this case do not allocate anything, set niggp=ng**2 and exit ===
  write(msg,'(2a)')ch10,&
&  ' findggp : not enough symmetries to reduce the number of G-pairs'
  call wrtout(std_out,msg,'COLL')
  Gpairs_q%niggp=Gsphere%ng**2
  RETURN
 end if 
 !
 ! === We can use symmetries to reduce the number of pairs ===
 ng=Gsphere%ng
 gmet(:,:)=Gsphere%gmet(:,:)
 gvec => Gsphere%gvec(1:3,1:ng)

 allocate(Gpairs_q%fp2ip(2,ng,ng),STAT=istat)
 if (istat/=0) call memerr(my_name,'fp2ip',ng**2,'i4b')

 allocate(Gpairs_q%fptabo(ng,ng),STAT=istat)
 if (istat/=0) call memerr(my_name,'fptabo',ng**2,'i4b')

 allocate(ip2fp(2,ng**2)) ! still do not know number of irred pairs 
 if (istat/=0) call memerr(my_name,'ip2fp',2*ng**2,'i4b')
 !
 ! === Precalculate G norm to speed the loop over G-pairs ===
 allocate(gnorm(ng))
 do ig=1,ng
  gnorm(ig)=normv(DBLE(gvec(:,ig)),gmet,'g')
 end do
 !
 ! === Precalculate |G1-G2|^2 square to speed loop over G pairs ===
 allocate(gdiff2(ng**2),stat=istat)
 if (istat/=0) call memerr('findggp','gdiff2',ng**2,'dp')
 allocate(iperm(ng**2),stat=istat)
 if (istat/=0) call memerr('findggp','iperm',ng**2,'i4b')

 ic=1
 g12=two*gmet(1,2) ; g23=two*gmet(2,3) ; g31=two*gmet(3,1)
 do ig1=1,ng
  a1=gvec(1,ig1) ; a2=gvec(2,ig1) ; a3=gvec(3,ig1)
  do ig2=1,ng
   d1=gvec(1,ig2)-a1 ; d2=gvec(2,ig2)-a2 ; d3=gvec(3,ig2)-a3
   gdiff2(ic)= d1*(gmet(1,1)*d1+g12*d2) &
&             +d2*(gmet(2,2)*d2+g23*d3) &
&             +d3*(gmet(3,3)*d3+g31*d1)
   iperm(ic)=ic ; ic=ic+1
  end do
 end do
 ! === Sort differences ===
 call sort_dp(ng**2,gdiff2,iperm,tol10)
 !
 ! === Loop over all all possible pairs (G1,G2), finding the irreducible ones === 
 ! Note that the pairs are addressed by ascending order of the norm of G1
 !
 niggp=0 ! Start with zero no. of pairs found
 do ic=1,ng**2

  ig1=(iperm(ic)-1)/ng+1
  ig2=iperm(ic)-(ig1-1)*ng
  diff2=gdiff2(ic)

  norm1=gnorm(ig1) ; eps1=tol8*norm1
  norm2=gnorm(ig2) ; eps2=tol8*norm2
  gposs1(:)=gvec(:,ig1)
  gposs2(:)=gvec(:,ig2)
  !
  ! === Check if this pair is the image through the same operation of a pair already found ===
  ! Consider only vectors with the same length
  found=.FALSE.
  if (niggp>0) then
   ip : do igpi=niggp,1,-1

    if (diff2-gdiff2(igpi)>eps1+eps2+48*tol10) EXIT  ! This is the step that makes the algorithm scale as N**2
    igpi1=ip2fp(1,igpi) ; if (ABS(norm1-gnorm(igpi1))>eps1) CYCLE
    igpi2=ip2fp(2,igpi) ; if (ABS(norm2-gnorm(igpi2))>eps2) CYCLE

    do itim=1,timrev_
     do isym=1,nsym
      ! === Calculate IS G1 and IS G2 ===
      ! Do this only for operations in the little group such as Sq=q
      ! TODO Calculate SG outside the loop, requires more memory but should be faster!
      !       avoid check on ltg, could make a table full==> little_group
      if (ltg(isym)==1) then 
       gxx1=(3-2*itim)*MATMUL(symrec(:,:,isym),gvec(:,igpi1))
       if (ALL(ABS(gxx1-gposs1)==0)) then 
        gxx2=(3-2*itim)*MATMUL(symrec(:,:,isym),gvec(:,igpi2))
        if (ALL(ABS(gxx2-gposs2)==0)) then 
         found=.TRUE.
         Gpairs_q%fp2ip(1,ig1,ig2)=igpi1
         Gpairs_q%fp2ip(2,ig1,ig2)=igpi2
         Gpairs_q%fptabo(ig1,ig2)=isym*(3-2*itim) ! Minus is time-reversal is considered (Only at Gamma)
         EXIT ip
        end if
       end if 
      end if 
     end do !isym 
    end do !itim
   end do ip
  end if
  
  if (.not.found) then
   ! === Increment counter and fill tables ===
   niggp=niggp+1
   gdiff2(niggp)=diff2
   ip2fp(1,niggp)=ig1
   ip2fp(2,niggp)=ig2
   Gpairs_q%fp2ip(1,ig1,ig2)=ig1
   Gpairs_q%fp2ip(2,ig1,ig2)=ig2
   Gpairs_q%fptabo(ig1,ig2)=iid  ! This irreducible pair comes from the identity!
  end if
 end do !ic
 deallocate(gnorm,gdiff2,iperm)

 if (niggp>ng**2) then 
  write(msg,'(2a)')ch10,&
&  ' findggp : BUG - niggp>ng**2 '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 

 write(msg,'(2a,i8,a,i8)')ch10,&
& ' findggp : number of independent (G,G'') pairs found = ',niggp,' / ',ng**2
 call wrtout(std_out,msg,'COLL')
 !
 ! === Save finale values ===
 Gpairs_q%niggp=niggp
 allocate(Gpairs_q%ip2fp(2,niggp),STAT=istat)
 if (istat/=0) call memerr(my_name,'ip2fp',2*niggp,'i4b')

 Gpairs_q%ip2fp=ip2fp(1:2,1:niggp)
 deallocate(ip2fp)

#if defined DEBUG_MODE
 call DEBUG_Gpairs(Gpairs_q,Gsphere,nsym,symrec)
 write(msg,'(1x,a)')TRIM(my_name)//': exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
end subroutine findggp
!!***


subroutine init_Gpairs_type(Gpairs_q,qpt,Gsphere,Cryst)

 use defs_basis
 use m_gwdefs, only : GW_TOLQ
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => init_Gpairs_type
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(Gvectors_type),intent(in) :: Gsphere
 type(Crystal_structure),intent(in) :: Cryst
 type(Gpairs_type),intent(out) :: Gpairs_q
!arrays
 real(dp),intent(in) :: qpt(3)

!Local variables ------------------------------
!scalars
 integer :: ng,istat,nsym
 integer,pointer :: symrec(:,:,:)
!************************************************************************

 nsym = Cryst%nsym
 symrec => Cryst%symrec

 ng=Gsphere%ng

 call destroy_Gpairs_type(Gpairs_q) 
 !
 !=== Dimensions ===
 Gpairs_q%ng=ng
 Gpairs_q%ngpairs=ng**2
 Gpairs_q%nsym=Gsphere%nsym
 Gpairs_q%timrev=Gsphere%timrev
 Gpairs_q%can_use_timrev=.FALSE.

 ! === The point under investigation ===
 Gpairs_q%qpt=qpt
 if (ALL(ABS(qpt)<GW_TOLQ).and.Gsphere%timrev==2) Gpairs_q%can_use_timrev=.TRUE.
 
 call findggp(nsym,symrec,Gsphere,qpt,Gpairs_q)

end subroutine init_Gpairs_type
!!***


subroutine destroy_Gpairs_type(Gpairs_q)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Gpairs_type),intent(inout) :: Gpairs_q

!************************************************************************

 if (ASSOCIATED(Gpairs_q%fp2ip )) DEALLOCATE(Gpairs_q%fp2ip )
 if (ASSOCIATED(Gpairs_q%fptabo)) DEALLOCATE(Gpairs_q%fptabo)
 if (ASSOCIATED(Gpairs_q%ip2fp )) DEALLOCATE(Gpairs_q%ip2fp )

end subroutine destroy_Gpairs_type
!!***


subroutine nullify_Gpairs_type(Gpairs_q)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 type(Gpairs_type),intent(inout) :: Gpairs_q

!************************************************************************

 NULLIFY(Gpairs_q%fp2ip )
 NULLIFY(Gpairs_q%fptabo)
 NULLIFY(Gpairs_q%ip2fp )

end subroutine nullify_Gpairs_type


subroutine DEBUG_Gpairs(Gpairs_q,Gsphere,nsym,symrec)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 type(Gpairs_type),intent(in) :: Gpairs_q
 type(Gvectors_type),intent(in) :: Gsphere
!arrays
 integer,intent(in) :: symrec(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,ng,ir1,ir2,isym,itim
 character(len=500) :: msg           
!arrays
 integer :: girr1(3),girr2(3),gxx1(3),gxx2(3),gtest1(3),gtest2(3)
 integer,pointer :: gvec(:,:)

!************************************************************************

 ng=Gpairs_q%ng
 gvec => Gsphere%gvec

 do ig1=1,ng
  gtest1=gvec(:,ig1)
  do ig2=1,ng
   gtest2=gvec(:,ig2)
   ir1=Gpairs_q%fp2ip(1,ig1,ig2)
   ir2=Gpairs_q%fp2ip(2,ig1,ig2)
   girr1=gvec(:,ir1)
   girr2=gvec(:,ir2)
   isym=Gpairs_q%fptabo(ig1,ig2) ; itim=1 
   if (isym<0) then 
    isym=-isym
    itim=2
   end if
   gxx1=(3-2*itim)*MATMUL(symrec(:,:,isym),girr1)
   gxx2=(3-2*itim)*MATMUL(symrec(:,:,isym),girr2)
   if (ANY((gtest1-gxx1)/=0).or.ANY((gtest2-gxx2)/=0)) then 
     write(msg,'(4a,2i8,a)')ch10,&
&     ' DEBUG_Gpairs : BUG -',ch10,&
&     ' G1-G2 pair',ig1,ig2,' has no corresponding irreducible pair ' 
     call wrtout(std_out,msg,'COLL')
     write(*,*)' G1, G2 ',gtest1,gtest2
     write(*,*)' independent pair? ',girr1,girr2
     write(*,*)' operation ',isym,symrec(:,:,isym),' with itim ',itim
     call leave_new('COLL')
    end if 
   end do
  end do 
end subroutine DEBUG_Gpairs
