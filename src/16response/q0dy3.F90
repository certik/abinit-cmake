!{\src2tex{textfont=tt}}
!!****f* ABINIT/q0dy3
!! NAME
!! q0dy3
!!
!! FUNCTION
!! Takes care of the inclusion of the ewald q=0 term in the dynamical
!! matrix
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom= number of atom in the unit cell
!!  option= either 0, 1 or 2. Described above ...
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  dyewq0(3,3,natom) = part needed to correct
!!    the dynamical matrix for atom self-interaction.
!!    ( produced by the routine if option=1 or 2, used if option=0 )
!!  dyew(2,3,natom,3,natom)= dynamical matrix
!!    ( input, non-corrected, for q=0 if option=1 or 2,
!!      input, non-corrected, for the real q if option=0
!!      output, corrected   , for the real q if option=0)
!!
!! NOTES
!! Should be used just after each call to ewald3, for both
!! q==0 and the real wavelength.
!!
!! If option=0, q0dy3 will use a previously produced Ewald dynamical
!! at q=0 (contracted), called dyewq0,
!! to impose the ASR to a bare Ewald dynamical
!! matrix (at any q), called dyew.
!! If option=1 or 2, q0dy3 use an Ewald dynamical matrix at q=0,
!! called dyew,
!! to produce a contracted form of it, called dyewq0 :
!! either an unsymmetrical form (if option=1), or a symmetrical
!! form (if option=2).
!!
!! The q0dy3 should be used in conjunction with the subroutine
!! ewald3 (or ewald9).
!! First, the call of ewald3 with q==0 should be done ,
!!   then the call to q0dy3 (with the option sumg0=0) will produce
!!   the dyewq0 matrix from the (q=0) dyew matrix
!! Second, the call of ewald3 with the real q (either =0 or diff 0)
!!   should be done, then the call to q0dy3 (with the option sumg0=1)
!!   will produce the correct dynamical matrix dyew starting from
!!   the previously calculated dyewq0 and the bare(non-corrected)
!!   dyew matrix
!!
!! PARENTS
!!      gtdyn9,hybrid9,mkifc9,respfn
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine q0dy3(natom,dyewq0,dyew,option)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,option
!arrays
 real(dp),intent(inout) :: dyew(2,3,natom,3,natom),dyewq0(3,3,natom)

!Local variables -------------------------
!scalars
 integer :: ia,ib,mu,nu
 character(len=500) :: message

! *********************************************************************

 if(option==1.or.option==2)then
  do mu=1,3
   do nu=1,3
    do ia=1,natom
     dyewq0(mu,nu,ia)=0.0_dp
     do ib=1,natom
      dyewq0(mu,nu,ia)=dyewq0(mu,nu,ia)+dyew(1,mu,ia,nu,ib)
     end do
    end do
   end do
  end do
 end if

 if(option==2)then
  do ia=1,natom
   do mu=1,3
    do nu=mu,3
     dyewq0(mu,nu,ia)=(dyewq0(mu,nu,ia)+dyewq0(nu,mu,ia))/2
     dyewq0(nu,mu,ia)=dyewq0(mu,nu,ia)
    end do
   end do
  end do
 end if

 if(option==0)then
  do mu=1,3
   do nu=1,3
    do ia=1,natom
     dyew(1,mu,ia,nu,ia)=dyew(1,mu,ia,nu,ia)-dyewq0(mu,nu,ia)
    end do
   end do
  end do
 end if

 if(option<0.or.option>2)then
  write(message, '(a,a,a,a,a,i4,a)' )&
&  ' q0dy3 : BUG -',ch10,&
&  '  The argument "option" should be equal either',ch10,&
&  '  to 0, 1 or to 2 , but it is',option,'.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

end subroutine q0dy3
!!***
