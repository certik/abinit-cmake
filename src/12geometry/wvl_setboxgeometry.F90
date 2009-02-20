!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_setBoxGeometry
!! NAME
!! wvl_setBoxGeometry
!!
!! FUNCTION
!! When wavelets are used, the box definition needs to be changed.
!! The box size is recomputed knowing some psp informations such as
!! the radius for coarse and fine grid. Then, the atoms are translated
!! to be included in the new box. Finally the FFT grid is computed using
!! the fine wavelet mesh and a buffer characteristic of used wavelets plus
!! a buffer used to be multiple of 2, 3 or 5.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  mpi_enreg=informations about MPI parallelization
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!
!! OUTPUT
!!  acell(3)=the new cell definition.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! SIDE EFFECTS
!!  dtset <type(dataset_type)>=internal variables used by wavelets, describing
!!                             the box are set. The FFT grid is also changed.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_setBoxGeometry(acell, dtset, mpi_enreg, radii, rprimd, xred)

 use defs_basis
  use defs_datatypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry, except_this_one => wvl_setBoxGeometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(inout) :: dtset
!arrays
 real(dp),intent(in) :: radii(dtset%ntypat,2)
 real(dp),intent(inout) :: acell(3),rprimd(3,3),xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: i,iat,itype
 character(len=500) :: message
!arrays
 real(dp) :: rprim(3,3)
 real(dp),allocatable :: xcart(:,:)
 character(len=20) :: atomnames(100)

! *********************************************************************

#if defined HAVE_BIGDFT
 if (dtset%prtvol == 0) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' wvl_setBoxGeometry : Changing the box for wavelets computation.'
  call wrtout(6,message,'COLL')
 end if

!Store xcart for each atom
 allocate(xcart(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

!Create the atom names
 do itype = 1, dtset%ntypat, 1
  write(atomnames(itype), "(A,I2)") "At. type", dtset%typat(itype)
 end do

 call system_size(mpi_enreg%me, dtset%natom, dtset%ntypat, xcart, radii, &
& dtset%wvl_crmult, dtset%wvl_frmult, dtset%wvl_hgrid, dtset%typat, &
& atomnames, acell(1), acell(2), acell(3), dtset%wvl_internal%nSize(1), &
& dtset%wvl_internal%nSize(2), dtset%wvl_internal%nSize(3), &
& dtset%wvl_internal%fineGrid(1, 1), dtset%wvl_internal%fineGrid(1, 2), &
& dtset%wvl_internal%fineGrid(1, 3), &
& dtset%wvl_internal%fineGrid(2, 1), dtset%wvl_internal%fineGrid(2, 2), &
& dtset%wvl_internal%fineGrid(2, 3))

 if (dtset%prtvol == 0) then
  write(message, '(a,3F12.6)' ) &
&  '  | acell is now:         ', acell
  call wrtout(6,message,'COLL')
  write(message, '(a,2I5,a,a,2I5,a,a,2I5)' ) &
&  '  | nfl1, nfu1:           ', dtset%wvl_internal%fineGrid(:, 1), ch10, &
&  '  | nfl2, nfu2:           ', dtset%wvl_internal%fineGrid(:, 2), ch10, &
&  '  | nfl3, nfu3:           ', dtset%wvl_internal%fineGrid(:, 3)
  call wrtout(6,message,'COLL')
 end if

!Change the metric to orthogonal one
 rprim(:, :) = real(0., dp)
 do i = 1, 3, 1
  rprim(i, i) = real(1., dp)
 end do
 call mkrdim(acell, rprim, rprimd)

!Save shifted atom positions into xred
 call xredxcart(dtset%natom, -1, rprimd, xcart, xred)
 deallocate(xcart)

!Set the number of points for the complete system (MPI or not)
 dtset%wvl_internal%nDpPoints = 1
 do i = 1, 3, 1
  dtset%wvl_internal%dpSize(i) = 2 * dtset%wvl_internal%nSize(i) + &
&  dtset%wvl_internal%buffer
  dtset%wvl_internal%nDpPoints = dtset%wvl_internal%nDpPoints * &
&  dtset%wvl_internal%dpSize(i)
 end do

 if (dtset%prtvol == 0) then
  write(message, '(a,3I12)' ) &
&  '  | box size for datas:   ', dtset%wvl_internal%dpSize
  call wrtout(6,message,'COLL')
  write(message, '(a,3I12)' ) &
&  '  | box size for wavelets:', dtset%wvl_internal%nSize
  call wrtout(6,message,'COLL')
 end if
 
#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_setBoxGeometry : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_setBoxGeometry
!!***
