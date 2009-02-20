!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_bonds
!! NAME
!! make_bonds
!!
!! FUNCTION
!!  (to be completed)
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ)
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!!
!! PARENTS
!!      make_prim_internals
!!
!! CHILDREN
!!      atmdata
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine make_bonds(bonds,nbond,dtset,icenter,nrshift,rprimd,rshift,xcart)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_16geomoptim, except_this_one => make_bonds
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icenter,nrshift
 integer,intent(out) :: nbond
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,pointer :: bonds(:,:,:)
 real(dp),intent(in) :: rprimd(3,3),rshift(3,nrshift),xcart(3,dtset%natom)

!Local variables ------------------------------
!scalars
 integer :: i1,i2,i3,iatom,ibond,ii,irshift,itypat,jatom
 real(dp) :: amu,bl,bondfudge,dist,rcov1,rcov2
 character(len=2) :: symbol
!arrays
 integer,allocatable :: bonds_tmp(:,:,:)
 integer,allocatable,target :: bonds_tmp2(:,:,:)
 real(dp) :: rcov(dtset%ntypat),rpt(3)

!************************************************************************

!DEBUG
!write (*,*) 'make_bonds: enter'
!ENDDEBUG
 do itypat=1,dtset%ntypat
  call atmdata(amu,rcov(itypat),symbol,dtset%znucl(itypat))
 end do

!DEBUG
!write (*,*) ' rcov =', rcov
!write (*,*) ' nrshift =', nrshift
!write (*,*) ' xcart =', xcart
!write (*,*) ' dtset%natom =',dtset%natom
!ENDDEBUG

!tentative first allocation: < 12 bonds per atom.
 allocate (bonds_tmp(2,2,12*dtset%natom))

 bondfudge = 1.1_dp

 nbond = 0

 do iatom=1,dtset%natom
  rcov1 = rcov(dtset%typat(iatom))
  do jatom=iatom+1,dtset%natom
   rcov2 = rcov(dtset%typat(jatom))
   do irshift=1,nrshift
    rpt(:) = xcart(:,jatom) &
&    + rshift(1,irshift)*rprimd(:,1) &
&    + rshift(2,irshift)*rprimd(:,2) &
&    + rshift(3,irshift)*rprimd(:,3)
    bl =  bond_length(xcart(:,iatom),rpt)

!   DEBUG
!   write (*,*) ' bl, bondfudge*(rcov1+rcov2) = ',bl, bondfudge*(rcov1+rcov2)
!   ENDDEBUG

    if (bondfudge*(rcov1+rcov2) - bl > tol6) then
     nbond = nbond+1
     if (nbond > 12*dtset%natom) then
      write (*,*) 'make_bonds : error too many bonds !'
      stop
     end if
     bonds_tmp(1,1,nbond) = iatom
     bonds_tmp(2,1,nbond) = icenter
     bonds_tmp(1,2,nbond) = jatom
     bonds_tmp(2,2,nbond) = irshift

!    DEBUG
!    write (*,*) ' ibond bonds = ', nbond, bonds_tmp(:,:,nbond),xcart(:,iatom),rpt
!    ENDDEBUG

    end if
   end do
!  end jatom do
  end do
 end do
!end iatom do

 allocate (bonds(2,2,nbond))
 do ibond=1,nbond
  bonds(:,:,ibond) = bonds_tmp(:,:,ibond)
 end do
 deallocate (bonds_tmp)
!bonds => bonds_tmp2

!DEBUG
!do ibond=1,nbond
!write (*,*) 'bond ', ibond, bonds(:,:,ibond)
!end do
!ENDDEBUG

end subroutine make_bonds
!!***
