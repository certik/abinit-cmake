!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawpolev
!! NAME
!! pawpolev
!!
!! FUNCTION
!! Compute the PAW term for polarization, named expected value term
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, PH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  natom=number of atoms in cell.
!!  nspden=number of spin-density components
!!  ntypat = number of atom types
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type (integer) for each atom
!!
!! OUTPUT
!! pelev(3)= electronic polarisation. expectation value term (PAW only)
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawpolev(natom,nspden,ntypat,pawprtvol,pawrhoij,pawtab,pelev,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: natom,nspden,ntypat,pawprtvol
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(out) :: pelev(3)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: iatom,idir,irhoij,ispden,itypat,klmn
 real(dp) :: ro_dlt
 character(len=500) :: message
!arrays
 integer :: idirindx(3)
 real(dp) :: tsec(2)

! *************************************************************************

!DEBUG
!write(6,*)' pawpolev : enter '
!ENDDEBUG

 call timab(560,1,tsec)

 idirindx(1)=4
 idirindx(2)=2
 idirindx(3)=3
 pelev=zero
 do idir=1,3
  do iatom=1,natom
   itypat=typat(iatom)
   do ispden=1,nspden
    do irhoij=1,pawrhoij(iatom)%nrhoijsel
     klmn=pawrhoij(iatom)%rhoijselect(irhoij)
     ro_dlt=pawrhoij(iatom)%rhoijp(irhoij,ispden)*pawtab(itypat)%dltij(klmn)
     pelev(idir)=pelev(idir)+ro_dlt*pawtab(itypat)%qijl(idirindx(idir),klmn)
     write (80,*)idir,iatom, klmn
     write (80,*) ro_dlt,pawtab(itypat)%qijl(idirindx(idir),klmn)
    end do
   end do
  end do
 end do

 write (80,*) 'pelev', pelev

 call timab(560,2,tsec)

!DEBUG
!write(6,*)' pawpolev : exit '
!ENDDEBUG

end subroutine pawpolev
!!***
