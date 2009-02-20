!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_all_gkk2
!!
!! NAME
!! get_all_gkk2
!!
!! FUNCTION
!! This routine determines where to store gkk2 matrix elements (disk or RAM)
!!   and calls interpolate_gkk to calculate them.
!!   This is the most time consuming step.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   acell = lengths of unit cell vectors
!!   amu = masses of atoms
!!   atmfrc = atomic force constants
!!   dielt = dielectric tensor
!!   dipdip = dipole-dipole contribution flag
!!   dyewq0 =
!!   elph_ds = datastructure for elphon data and dimensions
!!   FSkptirred = irreducible set of fermi-surface kpoints
!!   FSkpt = full set of fermi-surface kpoints
!!   ftwghtgkk = weights for FT of matrix elements
!!   gmet = metric in reciprocal space
!!   gprim = reciprocal lattice vectors
!!   indsym = indirect mapping of atoms under symops
!!   mpert = maximum number of perturbations
!!   msym = maximum number of symmetries (usually nsym)
!!   natom = number of atoms
!!   nrpt = number of real-space points for FT
!!   nsym = number of symmetries
!!   ntypat = number of types of atoms
!!   onegkksize = size of one gkk record, in bytes
!!   phon_ds = phonon datastructure containing data for eigen-val vec interpolation
!!   rcan = atomic positions in canonical coordinates
!!   rmet = real-space metric
!!   rprim = unit cell lattice vectors (dimensionless)
!!   rprimd = real-space unit-cell lattice vectors
!!   rpt = points in real space for FT, in canonical coordinates
!!   spqpt = qpoint coordinates
!!   symrel = symmetry operations in reduced real space
!!   trans = Atomic translations : xred = rcan + trans
!!   typat = array of types of atoms
!!   ucvol = unit cell volume
!!   wghatm = weights on rpt, for phonon FT interpolation
!!   xred = reduced coordinates of atoms
!!   zeff = Born effective charges
!!
!! OUTPUT
!!   elph_ds = calculated |gkk|^2 are in elph_ds%gkk2
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      interpolate_gkk
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_all_gkk2(acell,amu,atmfrc,dielt,dipdip,dyewq0,&
&    elph_ds,FSkptirred,FSkpt,&
&    ftwghtgkk,gmet,gprim,indsym,mpert,msym,&
&    natom,nrpt,nsym,ntypat,&
&    onegkksize,phon_ds,rcan,rmet,rprim,rprimd,&
&    rpt,spqpt,symrel,trans,typat,ucvol,&
&    wghatm,xred,zeff)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_17ddb, except_this_one => get_all_gkk2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dipdip,mpert,msym,natom,nrpt,nsym,ntypat,onegkksize
 real(dp),intent(in) :: ucvol
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
 real(dp),intent(in) :: FSkptirred(3,elph_ds%nFSkptirred),acell(3),amu(ntypat)
 real(dp),intent(in) :: atmfrc(2,3,natom,3,natom,nrpt),dielt(3,3)
 real(dp),intent(in) :: dyewq0(3,3,natom),ftwghtgkk(natom,nrpt),gmet(3,3)
 real(dp),intent(in) :: gprim(3,3),rcan(3,natom),rmet(3,3),rprim(3,3)
 real(dp),intent(in) :: rprimd(3,3),rpt(3,nrpt),spqpt(3,elph_ds%nqpt)
 real(dp),intent(in) :: trans(3,natom),wghatm(natom,natom,nrpt),xred(3,natom)
 real(dp),intent(in) :: zeff(3,3,natom)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,iost,jFSkpt,onediaggkksize
 real(dp) :: realdp_ex

! *************************************************************************

 if (elph_ds%nsppol /= 1) then
  stop 'get_all_gkk2: nsppol > 1 not coded yet!'
 end if

 onediaggkksize = elph_ds%nbranch*elph_ds%nFSkpt*kind(realdp_ex)

 elph_ds%unit_gkk2 = 37
 if (elph_ds%gkk2exist == 0) then
  if (elph_ds%gkk2write == 0) then
   write (*,*) 'get_all_gkk2 : keep gkk2 in memory. Size = ',&
&   4.0*dble(elph_ds%nFSkpt)*dble(onediaggkksize)/&
&   1024.0_dp/1024.0_dp, " Mb"
   allocate (elph_ds%gkk2(elph_ds%nbranch,elph_ds%ngkkband,elph_ds%ngkkband,&
&   elph_ds%nFSkpt,elph_ds%nFSkpt,1))
   elph_ds%gkk2(:,:,:,:,:,:) = zero
  else if (elph_ds%gkk2write == 1) then
   write (*,*) 'get_all_gkk2 : About to open gkk2 file : '
   write (*,*) elph_ds%unit_gkk2,onediaggkksize
   open (unit=elph_ds%unit_gkk2,file='gkk2file',access='direct',&
&   recl=onediaggkksize,form='unformatted',&
&   status='new',iostat=iost)
   if (iost /= 0) then
    write (*,*) 'get_all_gkk2 : error opening gkk2file as new'
    stop
   end if
!  rewind (elph_ds%unit_gkk2)
   write (*,*) 'get_all_gkk2 : disk file with gkk^2 created'
   write (*,*) '  calculate from real space gkk and phonon modes'
   write (*,*) '  gkk2write = 1 is forced: can take a lot of time! '
   write (*,*) ' size = ', 4.0*dble(onediaggkksize)*dble(elph_ds%nFSkpt)/&
&   1024.0_dp/1024.0_dp, ' Mb'
  else
   write (*,*) 'get_all_gkk2 : bad value of gkk2write'
   stop
  end if

 else if (elph_ds%gkk2exist == 1) then
  open (unit=elph_ds%unit_gkk2,file='gkk2file',access='direct',&
&  recl=onediaggkksize,form='unformatted',&
&  status='old',&
&  iostat=iost)
  if (iost /= 0) then
   write (*,*) 'get_all_gkk2 : error opening gkk2file as old'
   stop
  end if
! rewind (elph_ds%unit_gkk2)
  write (*,*) 'get_all_gkk2 : gkk2file was found. ',&
&  'values for gkk^2 will be taken from it'
  write (*,*) '    rather than recalculated'
  write (*,*) ' size = ', 4.0*dble(onediaggkksize)*dble(elph_ds%nFSkpt)/&
&  1024.0_dp/1024.0_dp, ' Mb'
 else
  write (*,*) 'get_all_gkk2 : bad value of gkk2exist'
  stop
 end if

!
!phfrq, if needed, for qpt on the full FSkpt grid
!
 elph_ds%unitphfrq = 38
 if (elph_ds%phfrqwrite == 1) then
  if (elph_ds%phfrqexist == 0) then
   open (unit=elph_ds%unitphfrq,file="phfrqfile",form='unformatted',&
&   access='direct', recl=3*natom*kind(realdp_ex),&
&   status='new',iostat=iost)
   if (iost /= 0) then
    write (*,*) 'get_all_gkk2 : error opening phfrqfile as new'
    stop
   end if
   write (*,*) ' get_all_gkk2 : file phfrqfile is open for writing of',&
&   ' FSqpt phonon frequencies'
!  rewind (elph_ds%unitphfrq)
  else if (elph_ds%phfrqexist == 1) then
   open (unit=elph_ds%unitphfrq,file="phfrqfile",&
&   form='unformatted',access='direct',&
&   recl=3*natom*kind(realdp_ex),status='old',iostat=iost)
   if (iost /= 0) then
    write (*,*) 'get_all_gkk2 : error opening phfrqfile as old'
    stop
   end if
   write (*,*) ' get_all_gkk2 : file phfrqfile is open for reading of',&
&   ' FSqpt phonon frequencies'
!  rewind (elph_ds%unitphfrq)
!  do iFSkpt=1,elph_ds%nFSkptirred
   do jFSkpt=1,elph_ds%nFSkpt
    read(elph_ds%unitphfrq,REC=jFSkpt) elph_ds%phfrq(:,jFSkpt)
   end do
!  end do
  end if
 end if

!
!here do the actual calculation of |g_kk|^2 if needed and
!
 if (elph_ds%phfrqexist == 0 .or. elph_ds%gkk2exist == 0) then
  call interpolate_gkk (acell,amu,atmfrc,dielt,dipdip,&
&  dyewq0,elph_ds,FSkptirred,FSkpt,ftwghtgkk,&
&  gmet,gprim,indsym,mpert,msym,natom,&
&  nrpt,nsym,ntypat,phon_ds,rcan,rmet,rprim,rprimd,rpt,spqpt,&
&  symrel,trans,typat,ucvol,wghatm,xred,zeff)
 end if

end subroutine get_all_gkk2
!!***
