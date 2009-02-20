!{\src2tex{textfont=tt}}
!!****f* ABINIT/bfactor
!! NAME
!! bfactor
!!
!! FUNCTION
!! Calculate the nesting factor
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nkpt = number of k-points
!!  kpt(3,nkpt) = k-point grid
!!  nqpt = number of qpoints
!!  qpt(3,nqpt) = q-point grid (must be a subgrid of the k grid),
!!                the nesting factor will be calculated for each q point in this array
!!  weight(nband,nkpt) =  integration weights for each k-point and band (NOT NORMALISED!!!)
!!  nband = number of bands
!!
!! OUTPUT
!!  nestfactor(nqpt) = array containing the nesting factor values
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Inspired to nmsq_gam_sumfs and mkqptequiv
!!  TODO : better use of symmetries to reduce the computational effort
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine bfactor(nkpt,kpt,nqpt,qpt,weight,nband,nestfactor)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_17ddb, except_this_one => bfactor
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nband,nkpt,nqpt
!arrays
 real(dp),intent(in) :: kpt(3,nkpt),qpt(3,nqpt),weight(nband,nkpt)
 real(dp),intent(inout) :: nestfactor(nqpt)

!Local variables-------------------------------
!scalars
 integer :: ib1,ib2,ikplusq,ikpt,iqpt,jkpt,symrank_kpt
 real(dp) :: factor,res,w1,w2
 character(len=500) :: message
!arrays
 integer,allocatable :: invrank_kpt(:),rank_kpt(:)
 real(dp) :: kptpq(3),tmpkpt(3)

! *************************************************************************

!DEBUG
!write(*,*)' bfactor : enter '
!ENDDEBUG

 nestfactor(:)=zero

 write (message,'(a)')' bfactor : making rank_kpt and invrank_kpt '
 call wrtout(06,message,'COLL')

 allocate (rank_kpt(nkpt),invrank_kpt(16000000))
 call mkkptrank (kpt,nkpt,rank_kpt,invrank_kpt)

 do iqpt=1,nqpt
  do ikpt=1,nkpt

   tmpkpt(:) = kpt(:,ikpt) + qpt(:,iqpt)
   call canon9(tmpkpt(1),kptpq(1),res)
   call canon9(tmpkpt(2),kptpq(2),res)
   call canon9(tmpkpt(3),kptpq(3),res)
   symrank_kpt = int(8000000.0_dp*(kptpq(1)+half+tol8) + &
&   40000.0_dp*(kptpq(2)+half+tol8) + &
&   200.0_dp*(kptpq(3)+half+tol8))
   ikplusq = invrank_kpt(symrank_kpt)
   if (ikplusq == -1) then
    write (message,'(4a)')ch10,' bfactor : ERROR- :',ch10,' it looks like no kpoint equiv to k+q !!!'
    call wrtout(06,message,'COLL')
!   DEBUG
!   write(*,*)kpt(:,ikpt),qpt(:,iqpt)
!   write(*,*)kpt(:,ikpt) + qpt(:,iqpt)
!   ENDDEBUG
    call leave_new('COLL')
   end if

   do ib1=1,nband
    w1 = weight(ib1,ikpt) !weight for distance from the Fermi surface
    if (w1 < tol6 ) cycle
    do ib2=1,nband
     w2 = weight(ib2,ikplusq) !weight for distance from the Fermi surface
     if (w1 < tol6 ) cycle
     nestfactor(iqpt) = nestfactor(iqpt)+w1*w2
    end do !ib2
   end do !ib1

  end do !ikpt
! DEBUG
! write (*,'(i4,2x,4(es14.6,2x))') iqpt,qpt(:,iqpt),nestfactor(iqpt)
! ENDDEBUG
 end do !iqpt

 deallocate (rank_kpt)
 deallocate(invrank_kpt)

!need prefactor of 1/nkpt for each integration over 1 kpoint index.
!and (1/nkpt)**2 for normalisation of double delta over FS
 factor=1./nkpt ; factor=factor**3
 nestfactor(:)=factor*nestfactor(:)

!DEBUG
!write(*,*)' bfactor : exit'
!stop
!ENDDEBUG

end subroutine bfactor
!!***
