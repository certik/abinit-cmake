!{\src2tex{textfont=tt}}
!!****f* ABINIT/pmat2cart
!! NAME
!! pmat2cart
!!
!! FUNCTION
!! Routine called by the program optic
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (SSharma,MVer,VRecoules,TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!!
!! OUTPUT
!!
!!
!! PARENTS
!!      optic
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pmat2cart(eigen11,eigen12,eigen13,mband,nkpt,nsppol,pmat,rprimd)

 use defs_basis

 implicit none

!Arguments -----------------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol
!arrays
 real(dp),intent(in) :: eigen11(2,mband,mband,nkpt,nsppol)
 real(dp),intent(in) :: eigen12(2,mband,mband,nkpt,nsppol)
 real(dp),intent(in) :: eigen13(2,mband,mband,nkpt,nsppol),rprimd(3,3)
!no_abirules
 complex(dpc),intent(out) :: pmat(mband,mband,nkpt,3,nsppol)

!Local variables -----------------------------------------
!scalars
 integer :: iband1,iband2,ii,ikpt,isppol
 real(dp) :: norm
!arrays
 real(dp) :: cartp(3),rprim(3,3),tmpp(3)

! *************************************************************************

!rescale the rprim
 rprim(:,:) = rprimd(:,:) / two_pi

 do isppol=1,nsppol
  do ikpt=1,nkpt
   do iband1=1,mband
    do iband2=1,mband
     pmat(iband2,iband1,ikpt,:,isppol) =             &
&     rprim(:,1)*cmplx(eigen11(1,iband2,iband1,ikpt,isppol),eigen11(2,iband2,iband1,ikpt,isppol),kind=dp) &
&     +rprim(:,2)*cmplx(eigen12(1,iband2,iband1,ikpt,isppol),eigen12(2,iband2,iband1,ikpt,isppol),kind=dp) &
&     +rprim(:,3)*cmplx(eigen13(1,iband2,iband1,ikpt,isppol),eigen13(2,iband2,iband1,ikpt,isppol),kind=dp)
    end do
   end do
  end do
 end do

end subroutine pmat2cart
!!***
