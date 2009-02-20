!{\src2tex{textfont=tt}}
!!****f* ABINIT/dsksta
!! NAME
!! dsksta
!!
!!
!! FUNCTION
!! This routine evaluates the amount of disk space required
!! by routine 'outkss'.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ishm=number of g-vecs shells to be saved in _KSS file
!!  nbandkss=number of desired bands to be saved in _KSS file
!!  npwkss=number of desired g-vecs to be saved in _KSS file
!!  nkpt=number of k points.
!!  nsym2=number of symmetries in space group, without INV
!!
!! OUTPUT
!!  Writes on standard output
!!
!! PARENTS
!!      outkss
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dsksta(ishm,nbandkss,npwkss,nkpt,nsym2)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ishm,nbandkss,nkpt,npwkss,nsym2

!Local variables-------------------------------
!
!scalars
 integer :: dsize
 character(len=500) :: message

! *********************************************************************
!
 dsize=2*80+16*4+9*8+nsym2*9*8+2+npwkss*3*8+ishm*4+nkpt*3*8&
& +4+8+nkpt*(2*nbandkss*8+16*npwkss*nbandkss)
!
 write(message,'(2a,f8.2,a)') ch10,&
& ' Amount of disk space required by _STA file=',&
& float(dsize)/float(1024*1024),' Mbytes.'
 call wrtout(6,message,'COLL')
!
end subroutine dsksta
!!***
