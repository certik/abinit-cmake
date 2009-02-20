!{\src2tex{textfont=tt}}
!!****f* ABINIT/strainsym
!! NAME
!! strainsym
!!
!! FUNCTION
!! For given order of point group, symmetrizes the strain tensor,
!! then produce primitive vectors based on the symmetrized strain.
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym=order of group.
!! rprimd(3,3)= primitive vectors, to be symmetrized
!! rprimd0(3,3)= reference primitive vectors, already symmetrized
!! symrel(3,3,nsym)=symmetry operators in terms of action on primitive translations
!!
!! OUTPUT
!! rprimd_symm(3,3)= symmetrized primitive vectors
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine strainsym(nsym,rprimd0,rprimd,rprimd_symm,symrel)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: rprimd(3,3),rprimd0(3,3)
 real(dp),intent(out) :: rprimd_symm(3,3)

!Local variables-------------------------------
!scalars
 integer :: isym
!arrays
 real(dp) :: rprimd0_inv(3,3),strain(3,3),strain_rot(3,3),strain_symm(3,3)
 real(dp) :: symd(3,3),symd_inv(3,3)

!**************************************************************************

!DEBUG
!write(6,*)' strainsym : enter '
!enddo
!ENDDEBUG

 ! copy initial rprimd input and construct inverse
 rprimd0_inv = rprimd0 
 call matrginv(rprimd0_inv,3,3)
 
 ! define strain as rprimd = strain * rprimd0, construct strain as
 ! strain = rprimd * rprimd0^{-1}
 strain = matmul(rprimd,rprimd0_inv)

 ! loop over symmetry elements to obtain symmetrized strain matrix
 strain_symm = zero
 do isym = 1, nsym

  ! obtain symd by converting symrel back to cartesian space
  symd = matmul(rprimd0,matmul(symrel(:,:,isym),rprimd0_inv))

  ! now get the inverse of symd
  symd_inv = symd
  call matrginv(symd_inv,3,3)

  ! apply the sym operator in cartesian space to the strain tensor
  strain_rot = matmul(symd_inv,matmul(strain,symd))

  ! accumulate the rotated strain tensors into the total
  strain_symm = strain_symm + strain_rot/nsym

 end do

 ! finally use the symmetrized strain tensor to make the final, new rprimd
 rprimd_symm = matmul(strain_symm,rprimd0)

!DEBUG
!rprimd_symm(:,:)=rprimd(:,:)
!ENDDEBUG

!DEBUG
!write(6,*)' strainsym : exit '
!ENDDEBUG

end subroutine strainsym
!!***
