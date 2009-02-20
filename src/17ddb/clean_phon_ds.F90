!{\src2tex{textfont=tt}}
!!****f* ABINIT/clean_phon_ds
!!
!! NAME
!! clean_phon_ds
!!
!! FUNCTION
!! This routine cleans structure phon_ds
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!!
!! INPUTS
!!   phon_ds = data structure for phonon interpolation - filled and allocated
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine clean_phon_ds(phon_ds)

 use defs_basis
 use defs_datatypes
 use defs_elphon

 implicit none

 !Arguments ------------------------------------
 !scalars
 type(phon_type),intent(inout) :: phon_ds
 !Local variables-------------------------------
! *************************************************************************
 deallocate(phon_ds%indsym)
 deallocate(phon_ds%symrel)
 deallocate(phon_ds%typat)
 deallocate(phon_ds%acell)
 deallocate(phon_ds%amu)
 deallocate(phon_ds%atmfrc)
 deallocate(phon_ds%dielt)
 deallocate(phon_ds%dyewq0)
 deallocate(phon_ds%gprim)
 deallocate(phon_ds%gmet)
 deallocate(phon_ds%xred)
 deallocate(phon_ds%zeff)
 deallocate(phon_ds%rcan)
 deallocate(phon_ds%rmet)
 deallocate(phon_ds%rprim)
 deallocate(phon_ds%rprimd)
 deallocate(phon_ds%rpt)
 deallocate(phon_ds%trans)
 deallocate(phon_ds%wghatm)
end subroutine clean_phon_ds
!!***

