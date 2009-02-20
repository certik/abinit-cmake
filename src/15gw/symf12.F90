!{\src2tex{textfont=tt}}
!!****f* ABINIT/symf12
!! NAME
!! symf12
!!
!! FUNCTION 
!!  This subroutine can be used in two different modes:
!!
!!  Mode 1) 
!!   Reconstruct all the reciprocal space matrix elements of a two point function F_{G1,G2}(q,omega) 
!!   for a fixed q-point, using only the set of independent (G1,G2) pairs.
!!   The independent pairs reconstruct the full G \times G space when all the operations of the little group of q are applied. 
!!  Note that F(r1,r2) is supposed to have the same symmetry as the crystal i.e 
!!   F(R^{-1}(r1-t),R^{-1}(r2-t)) = f(r1,r2)  where t is a fractional translation
!!
!!  Mode 2) 
!!   Check if a full matrix satisfies the above mentioned symmetry properties 
!!   (mainly used to debug the code)
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nsym=number of symmetry operations                                                                             
!!  npwc=number of G vectors in the matrix F_{G1,G2}(q)
!!  symrec(3,3,nsym)= symmetry operations in reciprocal space
!!  Gpairs_q<Gpairs_type>= tables related to the irreducible G1,G2 pairs for this q-point
!!   %ngpi=nuber of independent (G1,G1) pairs   
!!   %ng=Number of G vectors
!!   %ip2fp(2,ngpi)=index of G1 and G2 in the array gvec for each ngpi independent pair 
!!   %fptabo(ng,ng)=index of the symmetry operation S in the array symrec such as (T1,T2)= S (G1,G2)
!!  phgt(npw,nsym)=phase factors associated to non simmorphic operations (e^{-iG \cdot t})
!!
!! OUTPUT
!!  eps(npw,npw,nomega), see side effects
!!
!! SIDE EFFECTS
!!  Mode 1) 
!!   in input: eps contains only the matrix elements corresponding to the independent pair 
!!   in output: full symmetrized eps 
!!  Mode 2) 
!!   Only checking
!!   eps matrix is supposed to contain all the elements, 
!!
!! NOTES
!! 
!! TODO 
!!  We might take advantage of hermitianess in case of imaginary frequencies
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symf12(npw,nomega,omega,eps,Gsphere,Gpairs_q,nsym,symrec,phgt,mode)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mode,nomega,npw,nsym
 type(Gpairs_type),intent(in) :: Gpairs_q
 type(Gvectors_type),intent(in) :: Gsphere
!arrays
 integer,intent(in) :: symrec(3,3,nsym)
 complex,intent(in) :: omega(nomega),phgt(npw,nsym)
 complex,intent(inout) :: eps(npw,npw,nomega)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,iid,io,ioimax,iormax,ip,ir1,ir2,isym,itim
 real(dp),parameter :: tol=10d-1
 complex :: diff,phase
 logical :: found_identity
 character(len=500) :: msg
!arrays
 integer :: identity(3,3),iimax(2),irmax(2)
 real(dp) :: den(2),err(2),maxerr(2)

! *************************************************************************

 if (mode/=1.and.mode/=2) then
  write(msg,'(6a,i3)')ch10,&
&  ' symf12: BUG -',ch10,&
&  ' The mode option should be 1 or 2,',ch10,&
&  ' however, mode=',mode
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 
 identity(:,:)=RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))

 found_identity=.FALSE.
 do isym=1,nsym
  if (ALL(symrec(:,:,isym)==identity)) then 
   iid=isym
   found_identity=.TRUE. ; EXIT
  end if
 end do 
 
 if (.not.found_identity) then 
  write(msg,'(5a)')ch10,&
&  ' findggp : COMMENT -',&
&  ' Only the inversion was found in the set of symmetries read from the KSS file ',ch10,&
&  ' Likely old version of KSS file '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 else 
  write(msg,'(a,i4)')' found identity with index : ',iid
  call wrtout(std_out,msg,'COLL')
 end if 

 ! === The basic property used is ==
 ! If q=Sq then 
 ! $ F_{SG1,SG2}(q)=e^{+iS(G2-G1)\cdot\tau} F_{G1,G2)}(q)$
 ! Besides only at Gamma 
 ! $ F_{-SG1,-SG2)(0)=e^{-iS(G2-G1)\cdot\tau} F_{-G1,-G2}^*(0)$
 !
 SELECT CASE (mode)
 CASE (1) 
  ! === symmetrize Matrix ===
  do ig1=1,npw
   do ig2=1,npw 

    isym=Gpairs_q%fptabo(ig1,ig2) ; itim=1  
    if (isym<0) then  ! time-reversal only at Gamma, in this case isym is negative     
     isym=-isym ; itim=2
     stop "check equation"
    end if

    if (isym==iid.and.itim==1) CYCLE  ! This pair is irreducible
    ir1=Gpairs_q%fp2ip(1,ig1,ig2)     ! Indeces in the array gvec of the irred pair
    ir2=Gpairs_q%fp2ip(2,ig1,ig2)             
    ! TODO can speed up by calculating the phase before entering the loop, but require more memory
    phase=phgt(ig1,isym)*CONJG(phgt(ig2,isym)) 

    do io=1,nomega
     ! might also take into account hermiticity in case of imaginary frequencies
     eps(ig1,ig2,io)=eps(ir1,ir2,io)*phase
    end do 
    if (itim==2) eps(ig1,ig2,:)=CONJG(eps(ig1,ig2,:))
    
   end do 
  end do 

 CASE (2) 
  ! === Check symmetries ===
  write(msg,'(2a,f5.2,a)')ch10,' symf12 : checking symmetry properties, relative tolerance : ',tol,ch10
  call wrtout(std_out,msg,'COLL')

  !if (ALL(Gpairs_q%fpptabo==iid)) then 
  ! !  A non symmorphyc translation associated to the identity operator means that the 
  ! !  unit cell is not primitive (I think the code should realize and stop)
  ! write(msg,'(2a)')' symff12 : all the pairs are independent, return',ch10
  ! call wrtout(std_out,msg,'COLL')
  ! RETURN
  !end if 

  err=zero ; maxerr(:)=zero

  irmax=0  ; iimax=0
  iormax=0 ;ioimax=0

  write(msg,'(2a,f5.2)')' val1 ... val2  percent err'
  call wrtout(std_out,msg,'COLL')

  do ig1=1,npw
   do ig2=1,npw 

    isym=Gpairs_q%fptabo(ig1,ig2) ; itim=1 ! time-reversal only at Gamma, in this case isym is negative
    if (isym<0) then       
     isym=-isym ; itim=2
     stop "checl equation"
    end if

    if (isym==iid.and.itim==1) CYCLE  ! This pair is irreducible
    ir1=Gpairs_q%fp2ip(1,ig1,ig2)     ! Indeces in the array gvec of the irred pair 
    ir2=Gpairs_q%fp2ip(2,ig1,ig2)        
    phase=phgt(ig1,isym)*CONJG(phgt(ig2,isym))

    do io=1,nomega

     if (itim==1) then 
      diff=eps(ig1,ig2,io)-eps(ir1,ir2,io)*phase
     else 
      diff=eps(ig1,ig2,io)-CONJG(eps(ir1,ir2,io)*phase)
     end if 

     err(1)=ABS(REAL(diff))  ; den(1)=MAX(ABS(REAL(eps(ig1,ig2,io),dp)),tol12)
     err(2)=ABS(AIMAG(diff)) ; den(2)=MAX(ABS(REAL(AIMAG(eps(ig1,ig2,io)),dp)),tol12)

     if (den(1)*tol<err(1)) then 
      write(*,'(3(es18.6,4x),a,es18.6,a,6i4)')&
&      REAL(eps(ig1,ig2,io)),zero,REAL(eps(ir1,ir2,io)*phase),' > ',100*err(1)/den(1),' %',ig1,ig2,ir1,ir2,isym,io
      if (err(1)>maxerr(1)) then 
       maxerr(1)=err(1)
       irmax(1)=ig1
       irmax(2)=ig2
       iormax=io
      end if
     end if 
     
     if (den(2)*tol<err(2)) then 
      write(*,'(3(es18.6,4x),a,es18.6,a,6i4)')&
&      zero,AIMAG(eps(ig1,ig2,io)),AIMAG(eps(ir1,ir2,io)*phase),' > ',100*err(2)/den(2),' %',ig1,ig2,ir1,ir2,isym,io
      if (err(2)>maxerr(2)) then 
       maxerr(2)=err(2)
       iimax(1)=ig1
       iimax(2)=ig2
       ioimax=io
      end if
     end if 

    end do !iomega
   end do !ig2 
  end do !ig1

  if (iormax/=0) then 
   write(*,'(a,2(2x,f24.6))')' symf12 : max abs error for real part = ',maxerr(1)
   ig1=irmax(1) ; ig2=irmax(2)
   isym=Gpairs_q%fptabo(ig1,ig2)
   ir1=Gpairs_q%ip2fp(1,ip) ; ir2=Gpairs_q%ip2fp(2,ip)
   write(*,*)'ig1,ig2,ir1,ir2,iormax=',ig1,ig2,ir1,ir2,iormax
   write(*,'(2(2x,es18.6))')&
&   REAL(eps(ig1,ig2,iormax)),REAL(eps(ir1,ir2,iormax)*phase)
  end if 
  if (ioimax/=0) then
   write(*,'(a,2(2x,f24.6))')' symf12 : max abs error for imag part = ',maxerr(2)
   ig1=iimax(1) ; ig2=iimax(2)
   isym=Gpairs_q%fptabo(ig1,ig2)
   ir1=Gpairs_q%ip2fp(1,ip) ; ir2=Gpairs_q%ip2fp(2,ip)
   write(*,*)'ig1,ig2,ir1,ir2,iormax=',ig1,ig2,ir1,ir2,ioimax
   write(*,'(2(2x,es18.6))')AIMAG(eps(ig1,ig2,ioimax)),AIMAG(eps(ir1,ir2,ioimax)*phase)
  end if 
 CASE DEFAULT
  write(msg,'(4a,i3)')ch10,&
&  ' symf12 : BUG- ',ch10,&
&  ' wrong value for mode ',mode
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT

end subroutine symf12
!!***
