!{\src2tex{textfont=tt}}
!!****f* ABINIT/inpphon
!!
!! NAME
!! inpphon
!!
!! FUNCTION
!! Interpolate phonon frequencies and eigenvectors at 1 qpt
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   phon_ds = datastructure containing all needed interatomic force constants etc...
!!      for interpolation
!!   qpt = qpoint coordinate we want to interpolate at
!!
!! OUTPUT
!!   displ = phonon mode displacement vectors
!!   pheigval = eigenvalues of modes
!!   pheigvec = eigenvectors of modes
!!   phfrq = frequencies of modes
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      interpolate_gkk,mka2f,mkph_linwid,read_gkk
!!
!! CHILDREN
!!      gtdyn9,phfrq3
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine inpphon(displ,pheigval,pheigvec,phfrq,phon_ds,qpt)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_16response
 use interfaces_17ddb, except_this_one => inpphon
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(phon_type),intent(inout) :: phon_ds
!arrays
 real(dp),intent(inout) :: qpt(3)
 real(dp),intent(out) :: displ(2,3*phon_ds%natom,3*phon_ds%natom)
 real(dp),intent(out) :: pheigval(3*phon_ds%natom)
 real(dp),intent(out) :: pheigvec(2*3*phon_ds%natom*3*phon_ds%natom)
 real(dp),intent(out) :: phfrq(3*phon_ds%natom)

!Local variables-------------------------
!scalars
 integer :: symdynmat
 real(dp) :: qphnrm
!arrays
 real(dp) :: d2cart(2,3,phon_ds%mpert,3,phon_ds%mpert)

! *************************************************************************

 symdynmat = 0
 qphnrm = one

!frequencies for spqpt already calculated and incorporated in read_gkk
 call gtdyn9(phon_ds%acell,phon_ds%atmfrc,phon_ds%dielt,phon_ds%dipdip,&
& phon_ds%dyewq0,d2cart,phon_ds%gmet,phon_ds%gprim,phon_ds%mpert,phon_ds%natom,&
& phon_ds%nrpt,qphnrm,qpt,phon_ds%rcan,phon_ds%rmet,phon_ds%rprim,phon_ds%rpt,&
& phon_ds%trans,phon_ds%ucvol,phon_ds%wghatm,phon_ds%xred,phon_ds%zeff)

!diagonalize dynamical matrices
 call phfrq3(phon_ds%amu,displ,d2cart,pheigval,pheigvec,phon_ds%indsym,&
& phon_ds%mpert,phon_ds%nsym,phon_ds%natom,phon_ds%nsym,phon_ds%ntypat,phfrq,&
& qphnrm,qpt,phon_ds%rprimd,symdynmat,phon_ds%symrel,phon_ds%typat,phon_ds%ucvol,phon_ds%xred)


end subroutine inpphon
!!***
