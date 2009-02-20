!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_init_type_proj
!!
!! NAME
!! wvl_init_type_proj
!!
!! FUNCTION
!! Allocate and compute the access keys for the projectors when the positions
!! of the atoms are given. The array to store projectors
!! is also allocated, use wvl_free_type_proj() to free them after use.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=internal variables used by wavelets, describing
!!   | wvl_internal=desciption of the wavelet box.
!!   | natom=number of atoms.
!!  mpi_enreg=informations about MPI parallelization
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  proj <type(wvl_projector_type)>=projectors informations for wavelets.
!!   | keys=its access keys for compact storage.
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      atmdata,createprojectorsarrays,leave_new,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_init_type_proj(dtset, mpi_enreg, proj, psps, rprimd, xred)

 use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_projectors_type),intent(inout) :: proj
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: idata
 real(dp) :: amu,rcov
 character(len=500) :: message
!arrays
 real(dp),allocatable :: xcart(:,:)
 character(len=20) :: atomnames(100)

! *********************************************************************

#if defined HAVE_BIGDFT
!Consistency checks, are all pseudo true GTH pseudo with geometric informations?
 if (dtset%npsp /= dtset%ntypat) then
  write(message, '(a,a,a,a,I0,a,I0,a,a,a)' ) ch10,&
&  ' wvl_init_type_proj :  consistency checks failed,', ch10, &
&  '  dtset%npsp (', dtset%npsp, ') /= dtset%ntypat (', dtset%ntypat, ').', ch10, &
&  '  No alchemy pseudo are allowed with wavelets.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 do idata = 1, dtset%ntypat, 1
  if (.not. psps%gth_params%set(idata)) then
   write(message, '(a,a,a,a,I0,a,a,a)' ) ch10,&
&   ' wvl_init_type_proj :  consistency checks failed,', ch10, &
&   '  no GTH parameters found for type number ', idata, '.', ch10, &
&   '  Check your input pseudo files.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  if (.not. psps%gth_params%hasGeometry(idata)) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' wvl_init_type_proj :  consistency checks failed,', ch10, &
&   '  the given GTH parameters has no geometry informations.', ch10, &
&   '  Upgrade your input pseudo files to GTH with geometric informatoins.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  write(atomnames(idata), "(A)") repeat(" ", 20)
  call atmdata(amu, rcov, atomnames(idata), dtset%znucl(idata))
 end do

!Nullify optional pointers.
 nullify(proj%der)

!Store xcart for each atom
 allocate(xcart(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)

 call createProjectorsArrays(mpi_enreg%me, dtset%wvl_internal%nSize(1), &
& dtset%wvl_internal%nSize(2), dtset%wvl_internal%nSize(3), &
& xcart, dtset%natom, dtset%ntypat, dtset%typat, atomnames, &
& psps%gth_params%psppar, psps%pspcod, psps%gth_params%radii_cf, &
& dtset%wvl_cpmult, dtset%wvl_fpmult, dtset%wvl_hgrid, &
& proj%keys, proj%proj)
 write(message, '(a,a,a,a,I0)' ) ch10,&
& ' wvl_init_type_proj : allocate projectors data,', ch10, &
& '  size of the compressed array: ', proj%keys%nprojel
 call wrtout(6,message,'COLL')
 
!Deallocations
 deallocate(xcart)

#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_init_type_proj : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_init_type_proj
!!***
