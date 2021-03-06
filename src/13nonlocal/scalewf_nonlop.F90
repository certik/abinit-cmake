!{\src2tex{textfont=tt}}
!!****f* ABINIT/scalewf_nonlop
!! NAME
!! scalewf_nonlop
!!
!! FUNCTION
!! At the start of nonlop (or similar routines), as well as its end,
!! the wavefunctions, when stored with istwfk/=2,
!! need to be scaled (by a factor of 2 or 1/2),
!! except for the G=0 component.
!! Only the first spinor component is to be modified.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  istwf_k=storage mode of the vector
!!  mpi_enreg=informations about MPI parallelization
!!  npw=number of planewaves
!!  option=1 multiply by 2
!!        =2 multiply by 1/2
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  vect(2,npw)=vector that is rescaled
!!
!! NOTES
!!  XG030513 : MPIWF One should pay attention to the
!!  G=0 component, that will be only one one proc...
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine scalewf_nonlop(istwf_k,mpi_enreg,npw,option,vect)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!This type is defined in defs_mpi
!scalars
 integer,intent(in) :: istwf_k,npw,option
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(inout) :: vect(2,npw)

!Local variables-------------------------------
!scalars
 integer :: ipw
 real(dp) :: scale
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(6,*)' scalewf_nonlop : enter '
!ENDDEBUG

 if(istwf_k/=1)then

  if(npw<1)then
   write(message,'(a,a,a,a,a,a,i6)') ch10,&
&   ' scalewf_nonlop: BUG -',ch10,&
&   '  The argument npw should be a positive number,',ch10,&
&   '  however, npw=',npw
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  if(option/=1 .and. option/=2)then
   write(message,'(a,a,a,a,a,a,i6)') ch10,&
&   ' scalewf_nonlop: BUG -',ch10,&
&   '  The argument option should be 1 or 2,',ch10,&
&   '  however, option=',option
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  scale=two
  if(option==2)scale=half

! Storage for the Gamma point. The component of the G=0 vector
! should not be scaled, and no G=0 imaginary part is allowed.
  if(istwf_k==2)then
   if (mpi_enreg%me_g0==1) then
    vect(2,1)=zero
!   $OMP PARALLEL DO PRIVATE(ipw) SHARED(npw,scale,vect)
    do ipw=2,npw
     vect(1,ipw)=scale*vect(1,ipw)
     vect(2,ipw)=scale*vect(2,ipw)
    end do
!   $OMP END PARALLEL DO
   else
!   $OMP PARALLEL DO PRIVATE(ipw) SHARED(npw,scale,vect)
    do ipw=1,npw
     vect(1,ipw)=scale*vect(1,ipw)
     vect(2,ipw)=scale*vect(2,ipw)
    end do
!   $OMP END PARALLEL DO
   end if
  end if

! Other storage modes, for k points with time-reversal symmetry.
! All components should be scaled.
  if(istwf_k>2)then
!  $OMP PARALLEL DO PRIVATE(ipw) SHARED(npw,scale,vect)
   do ipw=1,npw
    vect(1,ipw)=scale*vect(1,ipw)
    vect(2,ipw)=scale*vect(2,ipw)
   end do
!  $OMP END PARALLEL DO
  end if

 end if ! istwf_k/=1

!DEBUG
!write(6,*)' scalewf_nonlop : exit'
!stop
!ENDDEBUG

end subroutine scalewf_nonlop
!!***
