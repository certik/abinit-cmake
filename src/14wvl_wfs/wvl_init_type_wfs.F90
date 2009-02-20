!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_init_type_wfs
!!
!! NAME
!! wvl_init_type_wfs
!!
!! FUNCTION
!! Compute the access keys for the wavefuncrions when the positions
!! of the atoms are given.
!!
!! For memory occupation optimisation reasons, the wavefunctions are not allocated
!! here. See the initialisation routines wvl_wfsinp_disk(), wvl_wfsinp_scratch()
!! and wvl_wfsinp_reformat() to do it. After allocation, use wvl_free_type_wfs()
!! to deallocate all stuff (descriptors and arrays).
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
!!  wfs <type(wvl_projector_type)>=wavefunctions informations for wavelets.
!!   | keys=its access keys for compact storage.
!!
!! SIDE EFFECTS
!!  xred(3,natom)=reduced dimensionless atomic coordinates (in fact IN but here
!!                because of INOUT xredxcart() behavior).
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      createwavefunctionsdescriptors,leave_new,wrtout,xredxcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine wvl_init_type_wfs(dtset, mpi_enreg, psps, rprimd, wfs, xred)

 use defs_basis
  use defs_datatypes
  use defs_wvltypes
#if defined HAVE_BIGDFT
  use BigDFT_API
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_wf_type),intent(inout) :: wfs
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: xred(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: iat,idata,ii
 real(dp),parameter :: eps_mach=1.d-12
 real(dp) :: tt
 logical :: parallel
 character(len=500) :: message
!arrays
 integer :: fGrid(2,3),nSize(3)
 real(dp),allocatable :: xcart(:,:)
 character(len=20) :: atomnames(100)

! *********************************************************************

 parallel = (mpi_enreg%nproc > 1)

#if defined HAVE_BIGDFT
!Consistency checks, are all pseudo true GTH pseudo with geometric informations?
 if (dtset%npsp /= dtset%ntypat) then
  write(message, '(a,a,a,a,I0,a,I0,a,a,a)' ) ch10,&
&  ' wvl_init_type_wfs:  consistency checks failed,', ch10, &
&  '  dtset%npsp (', dtset%npsp, ') /= dtset%ntypat (', dtset%ntypat, ').', ch10, &
&  '  No alchemy pseudo are allowed with wavelets.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 do idata = 1, dtset%ntypat, 1
  if (.not. psps%gth_params%set(idata)) then
   write(message, '(a,a,a,a,I0,a,a,a)' ) ch10,&
&   ' wvl_init_type_wfs:  consistency checks failed,', ch10, &
&   '  no GTH parameters found for type number ', idata, '.', ch10, &
&   '  Check your input pseudo files.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
  if (.not. psps%gth_params%hasGeometry(idata)) then
   write(message, '(a,a,a,a,a,a)' ) ch10,&
&   ' wvl_init_type_wfs:  consistency checks failed,', ch10, &
&   '  the given GTH parameters has no geometry informations.', ch10, &
&   '  Upgrade your input pseudo files to GTH with geometric informatoins.'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
 end do

 nSize = dtset%wvl_internal%nSize
 fGrid = dtset%wvl_internal%fineGrid

!Store xcart for each atom
 allocate(xcart(3, dtset%natom))
 call xredxcart(dtset%natom, 1, rprimd, xcart, xred)
 
!Nullify possibly unset pointers
 nullify(wfs%psi)
 nullify(wfs%hpsi)
 nullify(wfs%psit)
 nullify(wfs%psidst)
 nullify(wfs%hpsidst)
 nullify(wfs%ads)
 
!Static allocations.
!Count number of non-null occ values
 wfs%nstates = 0
 do ii = 1, size(dtset%occ_orig), 1
  if (dtset%occ_orig(ii) > tol8) then
   wfs%nstates = wfs%nstates + 1
  end if
 end do

 wfs%nstates_up = 0
 wfs%nstates_dn = 0
 allocate(wfs%eval(wfs%nstates))
 allocate(wfs%spinar(wfs%nstates))
 tt = dble(wfs%nstates) / dble(mpi_enreg%nproc)
 wfs%mbandp = int((1.d0 - eps_mach * tt) + tt)

!Compute spin array.
 if (dtset%nsppol == 2) then
  wfs%spinar(:) = dtset%occ_orig(1:wfs%nstates)
  if (dtset%fixmom < -real(90, dp)) then
   wfs%nstates_up = min(wfs%nstates / 2, wfs%nstates)
  else
   wfs%nstates_up = min(wfs%nstates / 2 + int(dtset%fixmom), wfs%nstates)
  end if
  wfs%nstates_dn = wfs%nstates - wfs%nstates_up
  wfs%spinar(wfs%nstates_up + 1:wfs%nstates) = &
&  -wfs%spinar(wfs%nstates_up + 1:wfs%nstates)
 else
  wfs%spinar(:) = one
 end if

 write(message, '(a,a)' ) ch10,&
& ' wvl_init_wfs_type: Create access keys for wavefunctions.'
 call wrtout(6,message,'COLL')
!We don't create an output of the atomic position, so the box size is ignored,
!(see the three 0.d0).
 call createWavefunctionsDescriptors( &
& mpi_enreg%me, mpi_enreg%nproc, dtset%nwfshist, &
& nSize(1), nSize(2), nSize(3), .false., dtset%wvl_hgrid, &
& dtset%natom, dtset%ntypat, dtset%typat, atomnames, 0.d0, 0.d0, 0.d0, xcart, &
& psps%gth_params%radii_cf, dtset%wvl_crmult, dtset%wvl_frmult, &
& wfs%keys, wfs%mvctrp, wfs%nstates, wfs%mbandp, &
& fGrid(1, 1), fGrid(2, 1), &
& fGrid(1, 2), fGrid(2, 2), &
& fGrid(1, 3), fGrid(2, 3), wfs%bounds)
!The memory is not allocated there for memory occupation optimisation reasons.

 write(message, '(a,2I8)' ) &
& '  | all orbitals have coarse segments, elements:', &
& wfs%keys%nseg_c, wfs%keys%nvctr_c
 call wrtout(6,message,'COLL')
 write(message, '(a,2I8)' ) &
& '  | all orbitals have fine   segments, elements:', &
& wfs%keys%nseg_f, 7 * wfs%keys%nvctr_f
 call wrtout(6,message,'COLL')

!Deallocations
 deallocate(xcart)

#else
 write(message, '(a,a,a,a)' ) ch10,&
& ' wvl_init_type_wfs : BigDFT library is not compiled.', ch10, &
& '   Action, used the flag --enable-bigdft when configuring.'
 call wrtout(6,message,'COLL')
 call leave_new('COLL')
#endif
end subroutine wvl_init_type_wfs
!!***
