!{\src2tex{textfont=tt}}
!!****f* ABINIT/identq
!! NAME
!! identq
!!
!! FUNCTION
!! Identify q-points in whole BZ
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  qibz(3,nqibz)=coordinates of the q-points in the IBZ
!!  timrev=if 2, time-reversal symmetry is used, 1 otherwise
!!  nqibz= number of q points in IBZ
!!  nqbzX= maximum number of q points in BZ
!!  nsym= number of symmetry operations
!!  symrec(3,3,nsym)= symmetry operations in reciprocal space
!!
!! OUTPUT
!!  qbz(3,nqbzX)= q-points in whole BZ
!!  qtab(nqbzX)= table giving for each q-point in the BZ (qBZ), the corresponding
!!   irreducible point (qIBZ), where qBZ= (IS) qIBZ and I is the inversion or the identity
!!  qtabi(nqbzX)= for each q-point in the BZ defines whether inversion has to be 
!!   considered in the relation qBZ=(IS) qIBZ (1 => only S; -1 => -S)  
!!  qtabo(nqbzX)= the symmetry operation S in the array op that takes qIBZ to each qBZ
!!  nqbz= no. of q-points in whole BZ
!!  wtq(nqibz)=weight of each irred q-point (normalized to one)
!!
!! NOTES
!!  For q close to zero only one q-point is counted (i.e. not the
!!  rotations and reflections of the little q, nor other little
!!  q-points in the qIBZ set)
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      dosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine identq(qibz,nqibz,nqbzX,symrec,nsym,timrev,wtq,qbz,qtab,qtabi,qtabo,nqbz)

 use defs_basis
 use m_numeric_tools, only : is_zero
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => identq
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nqbzX,nqibz,nsym,timrev
 integer,intent(out) :: nqbz
!arrays
 integer,intent(out) :: qtab(nqbzX),qtabi(nqbzX),qtabo(nqbzX)
 real(dp),intent(in) :: qibz(3,nqibz),symrec(3,3,nsym)
 real(dp),intent(out) :: qbz(3,nqbzX),wtq(nqibz)

!Local variables ------------------------------
!scalars
 integer :: div4,id,ii,iold,iq_bz,iq_ibz,isym,itim,jj,jqp,nqeq0,res
 real(dp) :: tolq0
 character(len=100) :: frmt
 character(len=500) :: msg
!arrays
 integer :: G0(3)
 real(dp) :: qnew(3)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' identq: identifying q-points'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 tolq0=0.001_dp !old behaviour

 ! Zero number of q-points found and zero small q vectors found
 nqbz=0 ; nqeq0=0
 !
 ! === Loop over q-points in IBZ ===
 do iq_ibz=1,nqibz
  wtq(iq_ibz)=zero
  ! 
  ! * If q is close to zero, treat it differently as described above
  if (is_zero(qibz(:,iq_ibz),tolq0)) then

   if (nqeq0==0) then
    nqeq0=1
    nqbz=nqbz+1
    if (nqbz>nqbzX) then
     write(msg,'(5a)')ch10,&
&     ' identq : BUG ',ch10,&
&     ' nqbzX too small ',ch10
     call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
    end if
    wtq(iq_ibz)=wtq(iq_ibz)+one
    qbz(:,nqbz)=qibz(:,iq_ibz)
    qtab(nqbz)=iq_ibz
    qtabo(nqbz)=1
    qtabi(nqbz)=1
   else 
    write(msg,'(6a)')ch10,&
&    ' identq : ERROR -,',ch10,&
&    ' it seems that there are at least two "small" q-points in the grid ',ch10,&
&    ' check q-grid and coding in identq.F90 '
    call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
   end if

  else
   ! === Loop over symmetry operations S and inversion/identity I ===
   do isym=1,nsym
    do itim=1,timrev
     ! * Form SI q
     call dosym(symrec(:,:,isym),itim,qibz(:,iq_ibz),qnew)
      !    
      ! Check whether it has already been found (to within a RL vector)
      ! Here there is a problem since if an umklapp G_o vector is required (Sq1 = q2 +¨G_) 
      ! then during the reconstruction of \tilde\espilon^{-1} we have to calculate G-G_o, see csigme.F90
     iold=0
     do iq_bz=1,nqbz
      if (is_samek(qnew,qbz(:,iq_bz),G0)) iold=iold+1
     end do

     if (iold==0) then ! we have a new q-point
      nqbz=nqbz+1
      if (nqbz>nqbzx) then
       write(msg,'(4a)')ch10,&
&       ' identq : BUG -',ch10,&
&       ' nqbzX too small '
       call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
      end if
      wtq(iq_ibz)=wtq(iq_ibz)+one
      qbz(:,nqbz)=qnew(:)
      qtab(nqbz)=iq_ibz
      qtabo(nqbz)=isym
      qtabi(nqbz)=3-2*itim
     end if
    end do
   end do
  end if

 end do

 ! === Normalize weights to 1 ===
 wtq(:) = wtq(:)/SUM(wtq)
 !
 ! === Print out results ===
 write(msg,'(2a,i4,3a)')ch10,&
& ' q-points in irreducible wedge (IBZ) ',nqibz,ch10,&
& ' q-points [reciprocal lattice units]:',ch10
 call wrtout(std_out,msg,'COLL')
 do jj=1,nqibz
  write(msg,'(i5,3f12.6)')jj,(qibz(ii,jj),ii=1,3)
  call wrtout(std_out,msg,'COLL')
 end do

 write(msg,'(3a,i2,3a,i4,2a)')ch10,ch10,' together with the ',nsym,&
& ' symmetry operations and inversion',ch10,&
& ' have yielded',nqbz,' q-points in Brillouin Zone (BZ):',ch10 
 call wrtout(std_out,msg,'COLL')

 write(frmt,*)'(i5,2x,4(3f7.3,2x))'
 div4=nqbz/4 ; res=mod(nqbz,4)

 do id=0,div4-1 
  jj=4*id+1 
  write(msg,frmt)jj,((qbz(ii,jqp),ii=1,3),jqp=jj,jj+3)
  call wrtout(std_out,msg,'COLL')
 end do
 if (res/=0) then
  write(frmt,*)'(i5,2x,',res,'(3f7.3,2x),a)'
  write(msg,frmt)4*div4+1,((qbz(ii,jqp),ii=1,3),jqp=4*div4+1,nqbz),ch10
  call wrtout(std_out,msg,'COLL')
 end if

#if defined DEBUG_MODE
 write(msg,'(a)')' identq : ended'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine identq
!!***
