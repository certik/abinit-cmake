!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_all_gkr
!!
!! NAME
!! get_all_gkr
!!
!! FUNCTION
!! This routine determines what to do with the rspace
!!   matrix elements of the el phon coupling (to disk or in memory),
!!   then reads those given in the gkq file and Fourier Transforms them
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!   gprim = reciprocal space lattice vectors
!!   natom = number of atoms
!!   nrpt = number of real-space points used for FT
!!   onegkksize = size of one record of the new gkk output file, in bytes
!!   rpt = positions of real-space points for FT
!!   spqpt = qpoint coordinates
!!   wghatm = weights for real-space rpt in FT
!!
!! OUTPUT
!!   elph_ds%gkr = real space elphon matrix elements.
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgkk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_all_gkr (elph_ds,gprim,natom,nrpt,onegkksize,rpt,spqpt,wghatm)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_17ddb, except_this_one => get_all_gkr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt,onegkksize
 type(elph_type),intent(inout) :: elph_ds
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),spqpt(3,elph_ds%nqpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt0,iost,ipert1,irpt,qtor

! *************************************************************************

!
!WARNING : disk file used for large arrays gkk_rpt and
!(eventually) gkk2
!
!allocate (gkk_rpt(2,elph_ds%nbranch,elph_ds%nFSband,elph_ds%nFSband,&
!&  elph_ds%nFSkpt,nrpt))
 elph_ds%unit_gkk_rpt = 36
!see if we already have gkk_rpt (in file)
 if (elph_ds%gkk_rptexist == 0) then
! see if the gkk_rpt should be written to a file (only available option now)
  if (elph_ds%gkk_rptwrite == 1) then
!  file is not present : we need to do the FT
   open (unit=elph_ds%unit_gkk_rpt,file='gkk_rpt_file',access='direct',&
&   recl=onegkksize,form='unformatted',&
&   status='new',iostat=iost)
   if (iost /= 0) then
    write (*,*) 'get_all_gkr : error opening gkk_rpt_file as new'
    stop
   end if
   write (*,*) ' get_all_gkr : will write real space gkk to a disk file.'
   write (*,*) ' size = ', 4.0*dble(onegkksize)*dble(nrpt)/&
&   1024.0_dp/1024.0_dp, ' Mb'

!  else if (elph_ds%gkk_rptwrite  == 0) then
  else
   write (*,*) ' get_all_gkr : will keep real space gkk in memory.'
   write (*,*) ' size = ', 4.0*dble(onegkksize)*dble(nrpt)/&
&   1024.0_dp/1024.0_dp, ' Mb'
   allocate (elph_ds%gkk_rpt(2,elph_ds%ngkkband*elph_ds%ngkkband,&
&   elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol,nrpt))
!  write (*,*) ' get_all_gkr: invalid value for gkk_rptwrite'
!  stop
  end if
  write (*,*) '    about to FT the recip space gkk to real space '
  qtor = 1

! 
! NOTE: should be very easy to parallelize!
! 
  iFSkpt0 = 1
  call ftgkk (wghatm,elph_ds%gkk_qpt,elph_ds%gkk_rpt,&
&  elph_ds%gkqwrite,elph_ds%gkk_rptwrite,gprim,1,natom,&
&  elph_ds%nFSkpt,elph_ds%ngkkband,elph_ds%nFSkpt,elph_ds%nqpt,&
&  nrpt,elph_ds%nsppol,qtor,rpt,spqpt,elph_ds%unit_gkk_rpt,elph_ds%unitgkq)

! call ftgkk (elph_ds,gprim,iFSkpt0,natom,nrpt,qtor,rpt,spqpt,wghatm)
  write (*,*) ' get_all_gkr : done with FT of gkk to real space'

! No longer need the gkk_qpt?
! if (elph_ds%gkqwrite == 0) deallocate (elph_ds%gkk_qpt)

  elph_ds%gkk_rptexist = 1

 else if (elph_ds%gkk_rptexist == 1) then
  open (unit=elph_ds%unit_gkk_rpt,file='gkk_rpt_file',access='direct',&
&  recl=onegkksize,form='unformatted',status='old',&
&  iostat=iost)
  if (iost /= 0) then
   write (*,*) 'get_all_gkr : error opening gkk_rpt_file as old'
   stop
  end if
  write (*,*) 'get_all_gkr : gkk_rpt_file was found. ',&
&  'Will take real space gkk from it '
  write (*,*) ' size = ', 4.0*dble(onegkksize)*dble(nrpt)/&
&  1024.0_dp/1024.0_dp, ' Mb'
 else
  write (*,*) ' get_all_gkr: invalid value for gkk_rptexist'
 end if

!!DEBUG
!Test the FT of the gkk elements.
!call test_ftgkk(elph_ds,gprim,natom,nrpt,rpt,spqpt,wghatm)
!!ENDDEBUG

!DEBUG
!do irpt=1,nrpt
!do ipert1=1,elph_ds%nbranch
!write (*,'(6(F16.5,1x))') elph_ds%gkk_rpt(:,ipert1,1,1,1,irpt)
!end do
!end do
!ENDDEBUG


end subroutine get_all_gkr
!!***
