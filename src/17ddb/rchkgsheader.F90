!{\src2tex{textfont=tt}}
!!****f* ABINIT/rchkgsheader
!!
!! NAME rchkGSheader
!! rchkgsheader
!!
!!
!! FUNCTION
!! This routine reads the GS header information in the GKK file and checks it
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  natom = number of atoms from DDB, for check
!!  FSkptirred = coordinates of the irreducible kpoints close to the FS
!!
!! OUTPUT
!!  hdr = header information
!!  nband = number of bands for rest of calculation
!!          should be the same for all kpts
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      hdr_io,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine rchkGSheader (hdr,natom,nband,unitgkk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13io_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,unitgkk
 integer,intent(out) :: nband
 type(hdr_type),intent(out) :: hdr

!Local variables-------------------------------
!scalars
 integer :: fform,ikpt,rdwr
 character(len=500) :: message

! *************************************************************************
!
!read in general header of _GKK file
!this is where we get nkpt, ngkpt(:,:)... which are also read in
!rdddb9 and inprep8. Probably should do some checking to avoid
!using ddb files from other configurations
!
 rewind(unitgkk)
 rdwr = 5 !read in header of file without rewinding it
 call hdr_io(fform,hdr,rdwr,unitgkk)
 if (fform == 0) then
  write (message,'(4a)')ch10,' rchkgsheader : ERROR-',ch10,' GKK header mis-read. fform == 0 '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!if (hdr%nsppol /= 1) then
!write (message,'(4a)')ch10,' rchkgsheader : ERROR-',ch10,' nsppol /= 1 is not treated yet '
!call wrtout(06,message,'COLL')
!call leave_new('COLL')
!end if

 if (hdr%natom /= natom) then
  write (message,'(5a)')ch10,' rchkgsheader : ERROR-',ch10,' natom in gkk file is',&
&  ' different from anaddb input '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 do ikpt=1,hdr%nkpt
  if (hdr%nband(ikpt) /= hdr%nband(1)) then
   write (message,'(6a)')ch10,&
&   ' rchkgsheader : ERROR-',ch10,&
&   ' Use the same number of bands for all kpts : ',ch10,&
&   ' could have spurious effects if efermi is too close to the last band '
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
 end do

 rdwr = 4 ! echo header to screen
 call hdr_io(fform,hdr,rdwr,6)

 nband=hdr%nband(1)

end subroutine rchkGSheader

!!***
