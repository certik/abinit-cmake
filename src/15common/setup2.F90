!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup2
!!
!! NAME
!! setup2
!!
!! FUNCTION
!! Call within main routine for setup of various arrays.
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
!!   | dedlnn=d(Etot)/d(log(npw)) (input variable, hartree)
!!   | ecut=kinetic energy cutoff for planewave basis (hartree)
!!   | natom=number of atoms in unit cell
!!   | nkpt=number of k points
!!   | wtk(nkpt)=integration weight associated with each k point
!!  iscf=parameter controlling scf or non-scf choice
!!  npwtot(nkpt)=number of planewaves in basis and boundary at each k point
!!  ucvol=unit cell volume (bohr^3)
!!  xred(3,natom)=starting reduced atomic coordinates
!!
!! OUTPUT
!!  epulay=Pulay basis set correction to total energy using method
!!   of Francis and Payne (reference below).
!!  start(3,natom)=copy of starting xred
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setup2(dtset,epulay,iscf,&
&  npwtot,start,ucvol,wfs,xred)

 use defs_basis
 use defs_datatypes
 use defs_wvltypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iscf
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: epulay
 type(dataset_type),intent(in) :: dtset
 type(wvl_wf_type),intent(in) :: wfs
!arrays
 integer,intent(in) :: npwtot(dtset%nkpt)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(out) :: start(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: ikpt,npw
 real(dp) :: arith,geom,wtknrm
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(6,*)' setup2 : enter '
!ENDDEBUG

 epulay=0._dp
 if (iscf>0) then

! Copy coordinates into array start
  start(:,:)=xred(:,:)

  if (dtset%usewvl == 0) then
!  Get average number of planewaves per k point:
!  both arithmetic and GEOMETRIC averages are desired--
!  need geometric average to use method of Francis and Payne,
!  J. Phys.: Condens. Matter 2, 4395-4404 (1990).
!  Also note: force k point wts to sum to 1 for this averaging.
!  (wtk is not forced to add to 1 in a case with occopt=2)
   arith=zero
   geom=one
   wtknrm=zero
   do ikpt=1,dtset%nkpt
    npw=npwtot(ikpt)
    wtknrm=wtknrm+dtset%wtk(ikpt)
    arith=arith+npw*dtset%wtk(ikpt)
    geom=geom*npw**dtset%wtk(ikpt)
   end do

!  Enforce normalization of weights to 1
   arith=arith/wtknrm
   geom=geom**(1.0_dp/wtknrm)

!  Compute Francis and Payne correction to total energy
!  (this number gets ADDED to the total energy to get the
!  corrected total energy)
   epulay=-dtset%dedlnn*&
&   log(geom/(ucvol*((2._dp*dtset%ecut)**1.5_dp)/(6._dp*pi**2)))
  end if

! Ensure portability of output thanks to tol8
  if (dtset%usewvl == 0) then
   write(message, '(a,2f12.3)' ) &
&   ' setup2: Arith. and geom. avg. npw (full set) are',arith+tol8,geom
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  else
#if defined HAVE_BIGDFT
   write(message, '(a,2I8)' ) ' setup2: nwvl coarse and fine are', &
&   wfs%keys%nvctr_c, wfs%keys%nvctr_f
   call wrtout(ab_out,  message,'COLL')
#endif
  end if

 end if

end subroutine setup2
!!***
