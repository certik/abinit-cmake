!{\src2tex{textfont=tt}}
!!****f* ABINIT/fred2fdeloc
!! NAME
!! fred2fdeloc
!!
!! FUNCTION
!!  calculate delocalized forces from reduced coordinate ones
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!! btinv(3*(dtset%natom-1),3*dtset%natom)= inverse transpose of B matrix (see delocint)
!! dtset <type(dataset_type)>=all input variables for this dataset
!! rprimd(3,3)=dimensional primitive translations (bohr)
!!
!! OUTPUT
!! deloc_force(3*(dtset%natom-1))=delocalized forces from reduced coordinate ones
!! fred(3,dtset%natom)=delocalized forces in reduced coordinates
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      delocint
!!
!! CHILDREN
!!      dgemm,dgemv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fred2fdeloc(btinv,deloc_force,fred,dtset,rprimd)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: btinv(3*(dtset%natom-1),3*dtset%natom),rprimd(3,3)
 real(dp),intent(out) :: deloc_force(3*(dtset%natom-1)),fred(3,dtset%natom)

!Local variables-------------------------------
!arrays
 real(dp) :: fcart(3,dtset%natom),tmpfcart(3*dtset%natom)

! ******************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DGEMV' :: dgemv
!DEC$ ATTRIBUTES ALIAS:'DGEMM' :: dgemm
#endif

!make cartesian forces

 call dgemm('N','N',3,dtset%natom,3,one,&
& rprimd,3,fred,3,zero,fcart,3)

!turn cartesian to delocalized forces

 tmpfcart = reshape(fcart,(/3*dtset%natom/))
 write (*,*) fcart
 write (*,*) tmpfcart

 call dgemv('N',3*(dtset%natom-1),3*dtset%natom,one,&
& btinv,3*(dtset%natom-1),fcart,1,zero,deloc_force,1)

 write (*,*) 'fred2fdeloc : deloc_force = '
 write (*,'(6E16.6)') deloc_force


end subroutine fred2fdeloc
!!***
