!{\src2tex{textfont=tt}}
!!****f* ABINIT/outbsd
!! NAME
!! outbsd
!!
!! FUNCTION
!! output bsd file for one perturbation (used for elphon calculations in anaddb)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, DRH, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.
!!
!! INPUTS
!!  dtset = dataset variable for run flags
!!  eigen0 = GS eigenvalues
!!  eigen1 = response function 1st order eigenvalue matrix
!!
!! OUTPUTS
!!  to file
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outbsd(outfile,dtset,eig2nkq,mpert,nkpt_rbz,unitout)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif 

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpert,nkpt_rbz,unitout
 character(len=fnlen),intent(in) :: outfile
 type(dataset_type),intent(inout) :: dtset
!arrays
 real(dp),intent(in) :: eig2nkq(2,dtset%mband*dtset%nsppol,nkpt_rbz,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: bandtot_index,fform,iband,idir1,idir2,ikpt,ipert1,ipert2,isppol
 integer :: rdwrout

! *********************************************************************

!initializations
 rdwrout = 4
!unitout should be attributed in dtset to avoid conflicts
 fform = 42

!open bsd file
 open (unit=unitout,file=outfile,form='formatted',status='unknown')
 
!output information in this file
 write(unitout,*)
 write(unitout,'(a,i8)') ' 2rd eigenvalue derivatives   - # elements :', 9*dtset%natom**2
 write(unitout,'(a,3es16.8,a)') ' qpt', dtset%qptn(:), ' 1.0'

!output RF eigenvalues

 do ikpt=1,nkpt_rbz
  bandtot_index=0
  write (unitout,'(a,3es16.8)') ' K-point:', dtset%kptns(:,ikpt)
  do isppol=1,dtset%nsppol
   do iband=1,dtset%mband
    write (unitout,'(a,i5)') ' Band:', iband+bandtot_index
!   write (unitout,*) 'ipert1     ','idir1     ','ipert2     ','idir2    ','Real    ','Im    '
    do ipert2=1,mpert-2
     do idir2=1,3
      do ipert1=1,mpert-2
       do idir1=1,3
        write (unitout,'(4i4,2d22.14)') idir1,ipert1,idir2,ipert2,&
&        eig2nkq(1,iband+bandtot_index,ikpt,idir1,ipert1,idir2,ipert2),&
&        eig2nkq(2,iband+bandtot_index,ikpt,idir1,ipert1,idir2,ipert2)
       end do !idir2
      end do !ipert2
     end do !idir1
    end do !ipert1
   end do !iband
   bandtot_index = bandtot_index + dtset%mband
  end do !isppol
 end do !ikpt

!close bsd file
 close (unitout)

end subroutine outbsd
!!***
