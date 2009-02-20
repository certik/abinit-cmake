!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_fc
!! NAME
!! calc_fc
!!
!! FUNCTION
!! calculation and output of Fermi-contact term at each atomic site
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gprimd(3,3)=matrix relating cartesian coords to crystal coords in reciprocal space
!!  natom=number of atoms in cell.
!!  nfft=number of points on fft grid
!!  ngfft(18)=details of fft
!!  nhat(nfft,nspden*usepaw)=PAW augmentation charge density on grid
!!  nspden=number of spin densities
!!  ntypat=number of atom types
!!  paral_kgb
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor(nfft,nspden)=electron density on grid (strictly $\tilde{n}+\hat{n}$)
!!  rprimd(3,3)=matrix relating cartesian coordinates to crystal coordinates
!!  typat(natom)=type (integer) for each atom
!!  xred(3,natom)=vectors locating each atom in the unit cell, in crystal coords
!!
!! OUTPUT
!!  (only writing, printing)
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine calc_fc(gprimd,natom,nfft,ngfft,nhat,nspden,ntypat,paral_kgb,&
&                   pawrad,pawrhoij,pawtab,psps,rhor,rprimd,typat,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13paw
 use interfaces_15common, except_this_one => calc_fc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,nspden,ntypat,paral_kgb
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: gprimd(3,3),nhat(nfft,nspden*psps%usepaw)
 real(dp),intent(in) :: rhor(nfft,nspden),rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom
 character(len=500) :: message
!arrays
 real(dp),allocatable :: fc(:),fc_el(:),fc_paw(:),gcart(:,:,:,:)

! ************************************************************************

!DEBUG
!write(*,*)' calc_fc : enter'
!ENDDEBUG
!Compatibility tests
 if (psps%usepaw /= 1) then
  write (message,'(4a)')' calc_efg : ERROR- ',ch10,&
&  ' usepaw /= 1 but Fermi-contact calculation requires PAW ',ch10
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 allocate(gcart(ngfft(1),ngfft(2),ngfft(3),3))
 call gridgcart(gcart,gprimd,ngfft) ! obtain G vectors in cartesian coords on grid

 allocate(fc(natom),fc_el(natom),fc_paw(natom)) 
 fc(:) = zero
 fc_el(:) = zero
 fc_paw(:) = zero

 call make_fc_el(fc_el,gcart,natom,nfft,ngfft,nhat,nspden,paral_kgb,rhor,rprimd,xred)
 call make_fc_paw(fc_paw,natom,ntypat,pawrhoij,pawrad,pawtab,psps,typat)

 fc(:) = fc_el(:) + fc_paw(:)

 write(message,'(a,a,a)' ) ch10,' Fermi-contact Term Calculation ',ch10
 call wrtout(ab_out,message,'COLL')

 do iatom = 1, natom
  write(message,'(a,i3,a,i3,a,f12.4)') ' Atom ',iatom,', typat ',typat(iatom),': FC = ',fc(iatom)
  call wrtout(ab_out,message,'COLL')
  write(message,'(a,f12.4,a,f12.4)') ' fc_el = ',fc_el(iatom),', fc_paw = ',fc_paw(iatom)
  call wrtout(ab_out,message,'COLL')
 end do

 write(message,'(3a)')ch10,ch10,ch10
 call wrtout(ab_out,message,'COLL')

 deallocate(gcart,fc,fc_el,fc_paw)

!DEBUG
!write(6,*)' calc_fc : exit '
!stop
!ENDDEBUG

 end subroutine calc_fc
!!***
