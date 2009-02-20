!{\src2tex{textfont=tt}}
!!****f* ABINIT/make_cs_dia
!! NAME
!! make_cs_dia
!!
!! FUNCTION
!! Compute the diamagnetic augmentation contribution to the chemical shielding
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  natom=number of atoms in cell.
!!  ntypat=number of atom types
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type (integer) for each atom
!!
!! OUTPUT
!!  cs(3,3,natom), the 3x3 diamagnetic chemical shielding tensor at each site

!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      deducer0,simp_gen
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine make_cs_dia(cs,natom,ntypat,pawang,pawrhoij,pawrad,pawtab,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(out) :: cs(3,3,natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom,ils,ilslm,irhoij,isel,ispden,itypat,klm,klmn,kln,lmax,lmin
 integer :: mesh_size,mm
 real(dp) :: c1,c2,c3,c4,intg,rg0,ro
!arrays
 real(dp) :: rg2(-2:2)
 real(dp),allocatable :: ff(:)

! ************************************************************************

!DEBUG
!write(*,*)' make_cs_dia : enter'
!ENDDEBUG

 cs(:,:,:) = zero

!the following factors arise in expanding the angular dependence of the chem shieldingtensor in
!terms of real spherical harmonics. The real spherical harmonics are as in the routine initylmr.F90; see
!in particular also http://www1.elsevier.com/homepage/saa/eccc3/paper48/eccc3.html
 c1 = sqrt(4.0*pi/15.0)
 c2 = sqrt(pi/45.0)
 c3 = sqrt(16.0*pi/9.0)
 c4 = 1.0E6/(2.0*Sp_Lt*Sp_Lt) ! 1/2c^2 in a.u. x 10^6 to get ppm (I'm a little unsure about the sign of
!this term, but with the present + sign I get a negative diamagnetic shielding as it should be)

!loop over atoms in cell
 do iatom = 1, natom
  itypat = typat(iatom)
  mesh_size=pawrad(itypat)%mesh_size
  allocate(ff(mesh_size))

! loop over spin components
  do ispden=1,pawrhoij(iatom)%nspden

!  Loop over non-zero elements of rhoij
   do irhoij=1,pawrhoij(iatom)%nrhoijsel
    klmn=pawrhoij(iatom)%rhoijselect(irhoij)
    klm =pawtab(itypat)%indklmn(1,klmn)
    kln =pawtab(itypat)%indklmn(2,klmn)
    lmin=pawtab(itypat)%indklmn(3,klmn)
    lmax=pawtab(itypat)%indklmn(4,klmn)

!   Computation of <phi_i|1/r|phi_j>- <tphi_i|1/r|tphi_j>
!   the diamagnetic CS tensor has radial dependence 1/r
    ff(2:mesh_size)=(pawtab(itypat)%phiphj(2:mesh_size,kln)&
&    -pawtab(itypat)%tphitphj(2:mesh_size,kln))&
&    /pawrad(itypat)%rad(2:mesh_size)
    call deducer0(ff,mesh_size,pawrad(itypat))
    call simp_gen(intg,ff,pawrad(itypat))

!   Select l=2 to get the rank 2 part
    if (lmin<=2.and.lmax>=2) then

!    Real gaunt coefficients selection
     ils=2 !l=2 only
     rg2(:)=zero
     do mm=-ils,ils
      ilslm=ils*ils+ils+mm+1
      isel=pawang%gntselect(ilslm,klm)
      if (isel>0) rg2(mm)=pawang%realgnt(isel) ! these are the non-zero <Ylm|Y2mm|Yl'm'> matrix elements
     end do

!    Accumulation of diamagnetic CS
!    the minus sign converts from electron density to charge density
     ro= -pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*intg
     cs(1,1,iatom) = cs(1,1,iatom)+c4*2.0*c2*rg2(0)*ro-c4*c1*rg2(2)*ro
     cs(2,2,iatom) = cs(2,2,iatom)+c4*2.0*c2*rg2(0)*ro+c4*c1*rg2(2)*ro
     cs(3,3,iatom) = cs(3,3,iatom)-c4*4.0*c2*rg2(0)*ro
     cs(1,2,iatom) = cs(1,2,iatom)-c4*c1*rg2(-2)*ro
     cs(1,3,iatom) = cs(1,3,iatom)-c4*c1*rg2(1)*ro
     cs(2,3,iatom) = cs(2,3,iatom)-c4*c1*rg2(-1)*ro

    end if  ! Select l=2
!   Select l = 0 to get isotropic part
    if (lmin<=0.and.lmax>=0) then
     rg0=0.0
     ilslm = 1
     isel=pawang%gntselect(ilslm,klm)
     if(isel>0) rg0 = pawang%realgnt(isel) ! these are the nonzero <Ylm|Y00|Yl'm'> matrix elements
!    Accumulation of diamagnetic CS
     ro= pawtab(itypat)%dltij(klmn)*pawrhoij(iatom)%rhoijp(irhoij,ispden)*intg
     cs(1,1,iatom) = cs(1,1,iatom)+c4*c3*rg0*ro
     cs(2,2,iatom) = cs(2,2,iatom)+c4*c3*rg0*ro
     cs(3,3,iatom) = cs(3,3,iatom)+c4*c3*rg0*ro

    end if  ! Select l=0
   end do   ! Loop on non-zero rhoij

  end do    ! Loop on spin components
  deallocate(ff)

! Symmetrization of CS
  cs(2,1,iatom) = cs(1,2,iatom)
  cs(3,1,iatom) = cs(1,3,iatom)
  cs(3,2,iatom) = cs(2,3,iatom)

 end do     ! Loop on atoms

!DEBUG
!write(6,*)' make_cs_dia : exit '
!stop
!ENDDEBUG

 end subroutine make_cs_dia
!!***
