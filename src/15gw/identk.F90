!{\src2tex{textfont=tt}}
!!****f* ABINIT/identk
!! NAME
!! identk
!!
!! FUNCTION
!! Identify k-points in whole BZ starting from the IBZ and generate symmetry tables.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  kibz(3,nkibz)=coordinates of k-points in IBZ
!!  timrev=if 2, time-reversal symmetry is considered; 1 otherwise
!!  nkibz= number of k points in IBZ
!!  nkbzmx= maximum number of k points in BZ
!!  nsym= number of symmetry operations
!!  symrec(3,3,nsym)= symmetry operation matrices in reciprocal space
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!
!! OUTPUT
!!  kbz(3,nkbzmx)= k-points in whole BZ
!!  ktab(nkbzmx)= table giving for each k-point in the BZ (array kbz), 
!!   the corresponding irreducible point in the array (kibz)
!!   k_BZ= (IS) kIBZ where S is one of the symrec operations and I is the inversion or the identity
!!    where k_BZ = (IS) k_IBZ and S = \transpose R^{-1} 
!!  ktabi(nkbzmx)= for each k-point in the BZ defines whether inversion has to be 
!!   considered in the relation k_BZ=(IS) k_IBZ (1 => only S; -1 => -S)  
!!  ktabo(nkbzmx)= the symmetry operation S that takes k_IBZ to each k_BZ
!!  nkbz= no. of k-points in the whole BZ
!!  wtk(nkibz)= weight for each k-point in IBZ for symmetric quantities:
!!              no. of distinct ks in whole BZ/(timrev*nsym)
!!
!!
!! PARENTS
!!      mrgscr,screening,sigma
!!
!! CHILDREN
!!      dosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine identk(kibz,nkibz,nkbzmx,nsym,timrev,symrec,symafm,use_antiferro,kbz,ktab,ktabi,ktabo,nkbz,wtk)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => identk
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkbzmx,nkibz,nsym,timrev
 integer,intent(out) :: nkbz
 logical,intent(in) :: use_antiferro
!arrays
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym)
 integer,intent(out) :: ktab(nkbzmx),ktabi(nkbzmx),ktabo(nkbzmx)
 real(dp),intent(in) :: kibz(3,nkibz)
 real(dp),intent(out) :: kbz(3,nkbzmx),wtk(nkibz)

!Local variables ------------------------------
!scalars
 integer :: div4,found,id,ii,ik,ikbz,ikibz,ikp,iold,isym,itim,jj,res
 character(len=100) :: fmt
 character(len=500) :: msg
!arrays
 integer :: G0(3)
 real(dp) :: knew(3)

! *************************************************************************

#if defined DEBUG_MODE
 !write(msg,'(a)')' identk : identifying k-points'
 !call wrtout(std_out,msg,'COLL')
#endif
 !
 ! === Loop over k-points in IBZ ===
 ! * Start with zero no. of k-points found
 nkbz=0 
 do ikibz=1,nkibz
  wtk(ikibz) = zero

  ! === Loop over time-reversal I and symmetry operations S  ===
  do itim=1,timrev
   do isym=1,nsym

    if (use_antiferro.and.symafm(isym)==-1) CYCLE
    !
    ! === Form IS k ===
    ! FIXME this is a temporary hacking to pass tests under gfortran, should use MATMUL
    call dosym(REAL(symrec(:,:,isym),dp),itim,kibz(:,ikibz),knew)
    !knew(:)=(3-2*itim)*MATMUL(symrec(:,:,isym),kibz(:,ikibz))

    ! === Check whether it has already been found (to within a RL vector) ===
    iold=0
    do ikbz=1,nkbz
     if (is_samek(knew,kbz(:,ikbz),G0)) iold=iold+1
    end do
    ! === If not yet found add to kbz and increase the weight ===
    if (iold==0) then
     nkbz=nkbz+1
     wtk(ikibz)=wtk(ikibz)+one
     if (nkbz>nkbzmx) then
      write(msg,'(4a,i6,2a)')ch10,&
&      ' identk: BUG ',ch10,&
&      ' nkbzmx too small, nkbzmx = ',nkbzmx,ch10,' increase nkbzmx !'
      call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
     end if
     kbz(:,nkbz)=knew(:)
     ktab(nkbz)=ikibz
     ktabo(nkbz)=isym
     ktabi(nkbz)=3-2*itim
    end if
   end do 
  end do 
  !wtk(ikibz)=wtk(ikibz)/(timrev*nsym)
 end do 

 ! === Normalize weights to 1 ===
 wtk(:) = wtk(:)/SUM(wtk)
 !
 ! === Print out results ===
 write(msg,'(2a,i3,2a,10x,2a)')ch10,&
& ' number of k-points in the irreducible wedge (IBZ) ',nkibz,ch10,&
& ' k-points [reciprocal lattice units]','weights',ch10
 call wrtout(std_out,msg,'COLL')

 write(fmt,*)'(i5,3f12.6,3x,f12.6)'
 do jj=1,nkibz
  write(msg,fmt) jj,(kibz(ii,jj),ii=1,3),wtk(jj)
  call wrtout(std_out,msg,'COLL')
 end do

 write(msg,'(2a,i2,3a,i5,a)')ch10,&
& ' together with ',nsym,' symmetry operations and inversion',ch10,&
& ' have yielded ',nkbz,' k-points in Brillouin Zone (BZ):'
 call wrtout(std_out,msg,'COLL')

 write(fmt,*)'(i5,2x,4(3f7.3,2x))'
 div4=nkbz/4 ; res=mod(nkbz,4)

 do id=0,div4-1
  jj=4*id+1 
  write(msg,fmt)jj,((kbz(ii,ikp),ii=1,3),ikp=jj,jj+3)
  call wrtout(std_out,msg,'COLL')
 end do

 if (res/=0) then
  write(fmt,*)'(i5,2x,',res,'(3f7.3,2x),a)'
  write(msg,fmt)4*div4+1,((kbz(ii,ik),ii=1,3),ik=4*div4+1,nkbz),ch10
  call wrtout(std_out,msg,'COLL')
 end if

#if defined DEBUG_MODE
 !write(msg,'(a)')' identk : exit'
 !call wrtout(std_out,msg,'PERS')
 !call flush_unit(std_out)
#endif

end subroutine identk
!!***
