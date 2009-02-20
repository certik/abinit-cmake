!{\src2tex{textfont=tt}}
!!****f* ABINIT/outgkk
!! NAME
!! outgkk
!!
!! FUNCTION
!! output gkk file for one perturbation (used for elphon calculations in anaddb)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG, DRH, MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  bantot0 = total number of bands for all kpoints
!!  bantot1 = total number of matrix elements for 1st order eigenvalues
!!  dtset = dataset variable for run flags
!!  eigen0 = GS eigenvalues
!!  eigen1 = response function 1st order eigenvalue matrix
!!  hdr0 = GS header
!!  hdr1 = RF header
!!
!! PARENTS
!!      loper3
!!
!! CHILDREN
!!      hdr_io
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outgkk(bantot0,bantot1,outfile,dtset,eigen0,eigen1,hdr0,hdr1)

 use defs_basis
 use defs_datatypes
#if defined HAVE_NETCDF
 use netcdf
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bantot0,bantot1
 character(len=fnlen),intent(in) :: outfile
 type(dataset_type),intent(inout) :: dtset
 type(hdr_type),intent(inout) :: hdr0,hdr1
!arrays
 real(dp),intent(in) :: eigen0(bantot0),eigen1(bantot1)

!Local variables-------------------------------
!scalars
 integer :: fform,iband,ikpt,isppol,ntot,rdwrout,unitout

! *************************************************************************

!initializations
 rdwrout = 6
!unitout should be attributed in dtset to avoid conflicts
 unitout = 111
 fform = 42
 ntot = 1

!open gkk file
 open (unit=unitout,file=outfile,form='unformatted',status='unknown')

!output GS header
 call hdr_io(fform,hdr0,rdwrout,unitout)

!output GS eigenvalues
 iband=0
 do isppol=1,hdr0%nsppol
  do ikpt=1,hdr0%nkpt
   write (unitout) eigen0(iband+1:iband+hdr0%nband(ikpt))
   iband=iband+hdr0%nband(ikpt)
  end do
 end do

!output number of gkk in this file (1)
 write (unitout) ntot

!output RF header
 call hdr_io(fform,hdr1,rdwrout,unitout)

!output RF eigenvalues
 iband=0
 do isppol=1,hdr1%nsppol
  do ikpt=1,hdr1%nkpt
   write (unitout) eigen1(iband+1:iband+2*hdr1%nband(ikpt)**2)
   iband=iband+2*hdr1%nband(ikpt)**2
  end do
 end do

!close gkk file
 close (unitout)

end subroutine outgkk
!!***
