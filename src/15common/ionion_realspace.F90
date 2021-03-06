!{\src2tex{textfont=tt}}
!!****f* ABINIT/ionion_realspace
!!
!! NAME
!! ionion_realspace
!!
!! FUNCTION
!! Compute the ion/ion interaction energies and forces in real space
!! case. Use ewald() instead if computations are done in reciprocal
!! space since it also includes the correction for the shift done in
!! potentials calculations and includes replica interactions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  rmet(3,3)=metric tensor in real space (bohr^2)
!!  xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!!  zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!!  eew=final ion/ion energy in hartrees
!!  grewtn(3,natom)=grads of ion/ion wrt xred(3,natom), hartrees.
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ionion_realSpace(dtset, eew, grewtn, rprimd, xred, zion)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: eew
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: rprimd(3,3),zion(dtset%ntypat)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(out) :: grewtn(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: ia1,ia2,iatom,igeo
 real(dp) :: r
 character(len=500) :: message
!arrays
 real(dp),allocatable :: grew_cart(:,:),xcart(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' ionion_realSpace : enter '
!stop
!ENDDEBUG

!Store xcart for each atom
 allocate(xcart(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

!Summing the interaction between ions.
 eew = 0._dp
 do ia1 = 1, dtset%natom, 1
  do ia2 = ia1 + 1, dtset%natom, 1
   r = sqrt((xcart(1, ia1) - xcart(1, ia2)) ** 2 + &
&   (xcart(2, ia1) - xcart(2, ia2)) ** 2 + &
&   (xcart(3, ia1) - xcart(3, ia2)) ** 2)
   eew = eew + zion(dtset%typat(ia1)) * zion(dtset%typat(ia2)) / r
  end do
 end do

!In Wavelet code, the ion-ion forces are already included
!in the local forces.
 if (dtset%usewvl == 0) then
  
! Allocate temporary array to store cartesian gradients.
  allocate(grew_cart(3, dtset%natom))

! Summing the forces for each atom
  do ia1 = 1, dtset%natom, 1
   grew_cart(:, ia1) = 0._dp
   do ia2 = 1, dtset%natom, 1
    if (ia1 /= ia2) then
     r = (xcart(1, ia1) - xcart(1, ia2)) ** 2 + &
&     (xcart(2, ia1) - xcart(2, ia2)) ** 2 + &
&     (xcart(3, ia1) - xcart(3, ia2)) ** 2
     do igeo = 1, 3, 1
      grew_cart(igeo, ia1) = grew_cart(igeo, ia1) - (xcart(igeo, ia1) - xcart(igeo, ia2)) * &
&      zion(dtset%typat(ia1)) * zion(dtset%typat(ia2)) / (r ** 1.5_dp)
     end do
    end if
   end do
  end do

  deallocate(xcart)

! Transform cartesian gradients to reduced gradients.
  do iatom = 1, dtset%natom, 1
   do igeo = 1, 3, 1
    grewtn(igeo, iatom) = rprimd(1, igeo) * grew_cart(1, iatom) + &
&    rprimd(2, igeo) * grew_cart(2, iatom) + &
&    rprimd(3, igeo) * grew_cart(3, iatom)
   end do
  end do
  deallocate(grew_cart)
 else
  grewtn(:,:) = zero
 end if

end subroutine ionion_realSpace
!!***
