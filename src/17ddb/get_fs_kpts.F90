!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_fs_kpts
!!
!! NAME
!! get_fs_kpts
!!
!! FUNCTION
!! This routine determines the kpoints on the standard grid which belong
!!  to the Fermi surface
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!    eigenGS = ground state eigenvalues
!!    hdr = header from input GS file
!!
!! OUTPUT
!!    FSkptflag   = flag to use a given kpoint for the fermi-surface integration
!!    gaussig     = width of gaussian energy window around fermi energy
!!                  needed to get a good fraction of kpoints contributing to the FS
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_fs_kpts(eigenGS,elph_ds,FSkptflag,gaussig,hdr)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: gaussig
 type(elph_type),intent(inout) :: elph_ds
 type(hdr_type),intent(in) :: hdr
!arrays
 integer,intent(out) :: FSkptflag(hdr%nkpt)
 real(dp),intent(in) :: eigenGS(hdr%nband(1),hdr%nkpt,hdr%nsppol)

!Local variables-------------------------------
!scalars
 integer :: iband,ikpt,isppol,nband
 real(dp) :: epsFS,gausstol
 character(len=500) :: message

! *************************************************************************

!supposes nband is equal for all kpts
 nband = hdr%nband(1)

!gausstol = minimum weight value for integration weights on FS
!should be set to reproduce DOS at Ef (Ref. PRB 34, 5065 p. 5067)
 gausstol = 1.0d-10

!use same band indices in both spin channels
 elph_ds%maxFSband=1
 elph_ds%minFSband=nband

!window of states around fermi Energy is contained in +/- epsFS
!should be adjusted to take into account a minimal but sufficient
!fraction of the kpoints: see the loop below.
!The 1000 is purely empirical!!!
!Should also take into account the density of kpoints.
!gaussig = width of gaussian for integration weights on FS
 
 gaussig = (maxval(eigenGS)-minval(eigenGS))/1000.0_dp
 
 write (message,'(a,f11.8,2a)')' get_fs_kpts : initial energy window = ',gaussig,ch10,&
& ' The window energy will be increased until the full k-grid is inside the range'
 call wrtout(06,message,'COLL')

!NOTE: could loop back to here and change gaussig until we have
!a certain fraction of the kpoints in the FS region...
 elph_ds%nFSkptirred = 0
 
 
!Do not use restricted fermi surface: include all kpts -> one
 do while (dble(elph_ds%nFSkptirred) / dble(hdr%nkpt) < one)
  gaussig = gaussig*1.05_dp

! DEBUG
! write (*,*)' now gaussig,elph_ds%nFSkptirred,hdr%nkpt = ',&
! &                   gaussig,elph_ds%nFSkptirred,hdr%nkpt
! ENDDEBUG

! we must take into account kpoints with states within epsFS:
  epsFS = gaussig*sqrt(log(one/(gaussig*sqrt(pi)*gausstol)))
  
! check if there are eigenvalues close to the Fermi surface
! (less than epsFS from it)
  FSkptflag(:) = 0

! do for each sppol channel
  do isppol=1,hdr%nsppol
   do ikpt=1,hdr%nkpt
    do iband=1,nband
     if (abs(eigenGS(iband,ikpt,isppol) - elph_ds%fermie) < epsFS) then
      FSkptflag(ikpt) = 1
      if (iband > elph_ds%maxFSband) elph_ds%maxFSband = iband
      if (iband < elph_ds%minFSband) elph_ds%minFSband = iband
     end if
    end do
   end do
  end do ! isppol
  
  elph_ds%nFSband = elph_ds%maxFSband-elph_ds%minFSband+1
  
! number of irreducible kpoints (by all sym) contributing to the Fermi
! surface (to be completed by symops)
  elph_ds%nFSkptirred = sum(FSkptflag(:))

  write (*,*) ' epsFS = ',epsFS, ' N = ',elph_ds%nFSkptirred
 end do

end subroutine get_fs_kpts
!!***
