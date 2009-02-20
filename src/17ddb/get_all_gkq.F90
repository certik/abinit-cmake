!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_all_gkq
!!
!! NAME
!! get_all_gkq
!!
!! FUNCTION
!! This routine determines what to do with the initial qspace
!!   matrix elements of the electron phonon coupling (to disk or in memory),
!!   then reads those given in the gkk file and completes them
!!   (for kpts, then perturbations, then qpoints on the whole of spqpt)
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer, MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   acell = lattice vector lengths
!!   amu = atomic masses
!!   elph_ds = elphon datastructure with data and dimensions
!!   FSfullpqtofull = mapping of k+q to another k
!!   FSfulltofull = mapping of FS kpoints under symops
!!   FSkpt = fermi surface kpoints
!!   FSintweight = integration weights for bands and kpoints near the FS
!!   gkk_flag = flag to
!!   gprimd = reciprocal space lattice vectors
!!   indsym = mapping of atoms under symops
!!   mpert = maximum number of perturbations
!!   natom = number of atoms
!!   nband = number of bands
!!   nqptirred = number of irreducible qpoints
!!   nsym = number of symmetries
!!   ntypat = number of types of atoms
!!   n1wf = number of file headers from perturbation calculations
!!      which are present in the initial gkk input file.
!!   onegkksize = size of one record of the new gkk output file, in bytes
!!   phon_ds = phonon datastructure for interpolation of eigen vec and val
!!   qptirred = irreducible qpoint coordinates
!!   qpttoqpt = mapping of qpoints onto each other under symmetries
!!   rprimd = real space lattice vectors (dimensionful)
!!   spqpt = qpoint coordinates
!!   symrec = reciprocal space symops
!!   symrel = real space symops
!!   timrev = flag to use time reversal symmetry
!!   tnons = translation vectors associated to symrel
!!   typat = array of types of atoms
!!   ucvol = unit cell volume
!!   unitgkk = fortran unit for initial gkk input file
!!   xred = reduced coordinates of atoms
!!
!! OUTPUT
!!   elph_ds%gkq = recip space elphon matrix elements.
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      complete_gkk,leave_new,read_gkk,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_all_gkq (acell,amu,elph_ds,FSfullpqtofull,FSfulltofull,FSkpt,FSintweight,   &
&    gkk_flag,gprimd,indsym,mpert,natom,nband,nqptirred,nsym,ntypat,n1wf,onegkksize,phon_ds,&
&    qptirred,qpttoqpt,rprimd,spqpt,symrec,symrel,timrev,tnons,typat,ucvol,unitgkk,xred)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_17ddb, except_this_one => get_all_gkq
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpert,n1wf,natom,nband,nsym,ntypat,onegkksize,timrev
 integer,intent(in) :: unitgkk
 integer,intent(inout) :: nqptirred
 real(dp),intent(in) :: ucvol
 type(elph_type),intent(inout) :: elph_ds
 type(phon_type),intent(inout) :: phon_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
 integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt),indsym(4,nsym,natom)
 integer,intent(in) :: qpttoqpt(2,nsym,elph_ds%nqpt),symrec(3,3,nsym)
 integer,intent(in) :: symrel(3,3,nsym),typat(natom)
 integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt)
 real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt)
 real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt),acell(3),amu(ntypat)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3),tnons(3,nsym),xred(3,natom)
 real(dp),intent(inout) :: qptirred(3,n1wf),spqpt(3,elph_ds%nqpt)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,ib1,ib2,ierr,imem,iost,iqpt,iqptfull,memsize
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),allocatable :: gkk_tmp_full(:,:,:),tmp(:),workspace(:,:,:,:,:,:,:)

! *************************************************************************

!attribute file unit number
 elph_ds%unitgkq = 35

!============================================
!save gkk for all qpts in memory or to disk
!============================================

 if (elph_ds%gkqexist == 0) then !if the gkk(q) are not in a file yet

  if (elph_ds%gkqwrite == 0) then !calculate gkk(q) keeping all in memory

   write(message,'(a,f14.4,a)')&
&   ' get_all_gkq : keep gkk(q) in memory (Size = ',&
   4.0*dble(onegkksize)*dble(elph_ds%nqpt)/1024.0_dp/1024.0_dp,' Mb)'
   call wrtout(06,message,'COLL')

   allocate(elph_ds%gkk_qpt(2,elph_ds%ngkkband*elph_ds%ngkkband,&
&   elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt),stat=ierr)
   if (ierr /= 0 ) then 
    write (message,'(3a)')' get_all_gkq : ERROR- ',ch10,&
&    ' trying to allocate array elph_ds%gkk_qpt '
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if

   elph_ds%gkk_qpt(:,:,:,:,:,:) = zero

  else if (elph_ds%gkqwrite == 1) then !calculate gkk(q) and write to file
   
   fname=trim(elph_ds%elph_base_name) // '_GKKQ'
   open (unit=elph_ds%unitgkq,file=fname,access='direct',recl=onegkksize,&
&   form='unformatted',status='new',iostat=iost)
!  rewind (elph_ds%unitgkq)
   if (iost /= 0) then
    write (message,'(3a)')' get_all_gkq : ERROR- opening file ',trim(fname),' as new'
    call wrtout(06,message,'COLL')
    call leave_new('COLL')
   end if

   write (message,'(6a,f14.4,a)')&
&   ' get_all_gkq : gkq matrix elements  will be written to file : ',trim(fname),ch10,&
&   ' Nothing is in files yet',ch10,                                                  &
&   ' Size = ',4.0*dble(onegkksize)*dble(elph_ds%nqpt)/1024.0_dp/1024.0_dp,' Mb'
   call wrtout(06,message,'COLL')

  else
   write(message,'(4a,i4)')ch10,&
&   ' get_all_gkq : BUG- ',ch10,&
&   ' gkqwrite must be 0 or 1 while it is : ',elph_ds%gkqwrite
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if !if gkqwrite

! =====================================================
! read in g_kk matrix elements for all bands, kpoints,
! and calculated qpoints
! =====================================================

  write(message,'(a)')' get_all_gkq : calling read_gkk to read in the g_kk matrix elements' 
  call wrtout(06,message,'COLL')

  call read_gkk(amu,elph_ds,FSfullpqtofull,FSfulltofull,FSintweight,FSkpt,&
&  gkk_flag,gprimd,indsym,n1wf,natom,nband,nqptirred,nsym,ntypat,&
&  phon_ds,qptirred,rprimd,spqpt,symrec,symrel,timrev,tnons,typat,ucvol,unitgkk)

  write (message,'(2a)')ch10,' get_all_gkq : out of read_gkk'
  call wrtout(06,message,'COLL')

  if (elph_ds%tsymgkq ==1) then

!  ==============================================================
!  complete gkk matrices for other qpoints on the full grid spqpt
!  inspired and cannibalized from symdm9.f
!  ==============================================================

   write(message,'(4a)')ch10,&
&   ' get_all_gkq : calling complete_gkk to complete ',ch10,&
&   ' gkk matrices for other qpoints on the full grid'
   call wrtout(06,message,'COLL')

!  DEBUG(1) MG
!  WARNING: in case of doscalprod=0 large modifications occur in the symmetrised matrix elements with respect
!  to the values in input
!  allocate (workspace(2,elph_ds%ngkkband,elph_ds%ngkkband,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nqpt))
!  do iqptfull=1,elph_ds%nqpt
!  if (elph_ds%wtq(iqptfull)/=0) then
!  workspace(:,:,:,:,:,iqptfull)=elph_ds%gkk_qpt(:,:,:,:,:,iqptfull)
!  end if
!  end do
!  ENDDEBUG(1)

   call complete_gkk(acell,elph_ds,FSfulltofull,FSkpt,gkk_flag,gprimd,indsym,mpert,&
&   natom,nqptirred,nsym,qptirred(:,1:nqptirred),qpttoqpt,rprimd,spqpt,&
&   symrec,symrel,tnons,ucvol,xred)

!  DEBUG(2) MG
!  calculate difference btw input/output gkk_qpt matrix elements
!  do iqptfull=1,elph_ds%nqpt
!  if (elph_ds%wtq(iqptfull)/=0) then
!  write(17,*)abs(elph_ds%gkk_qpt(:,:,:,:,:,iqptfull)-workspace(:,:,:,:,:,iqptfull))
!  end if
!  end do
!  deallocate (workspace)
!  ENDDEBUG(2)

   write (message,'(2a)')ch10,' get_all_gkq : out of complete_gkk'
   call wrtout(06,message,'COLL')

  end if !tsymgkq

! =======================================================================
! we already have full gkq in a file: read in stuff
! =======================================================================

 else if (elph_ds%gkqexist == 1) then !we already have full gkq in a file

  fname=trim(elph_ds%elph_base_name) // '_GKKQ'
  open (unit=elph_ds%unitgkq,file=fname,access='direct',recl=onegkksize,&
&  form='unformatted',status='old',iostat=iost)
  if (iost /= 0) then
   write (message,'(5a)')' get_all_gkq : ERROR-',ch10,' opening file ',trim(fname),' as old'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

! rewind (elph_ds%unitgkq)
  write (message,'(7a,f14.4,a)')&
&  ' get_all_gkq : file ',trim(fname),' was found.',ch10,&
&  ' values for gkk(q) will be taken from it, rather than recalculated ',ch10,&
&  ' Size = ',4.0*dble(onegkksize)*dble(elph_ds%nqpt)/1024.0_dp/1024.0_dp,' Mb'
  call wrtout(06,message,'COLL')

! 
! MJV 20070529  Problem here: the idea is not to have everything in memory !
! Should work as is: segments are read in in each subroutine.
! May have to add a check that elph_ds%gkqexist==1 is not used with
! elph_ds%gkqwrite==0
! 
! allocate (elph_ds%gkk_qpt(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
! &          elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt),stat=ierr)
! if (ierr /= 0 ) then
! write (message,'(3a)')' get_all_gkq : ERROR- ',ch10,&
! &   ' trying to allocate array elph_ds%gkk_qpt '
! call wrtout(06,message,'COLL')
! call leave_new('COLL')
! end if
! 
! elph_ds%gkk_qpt(:,:,:,:,:,:) = zero
! 
! do iqptfull=1,elph_ds%nqpt
! read(elph_ds%unitgkq,REC=iqptfull)elph_ds%gkk_qpt(:,:,:,:,:,iqptfull)
! end do
! 
! write(message,'(2a)')' get_all_gkq : done reading file ',trim(fname)
! call wrtout(06,message,'COLL')

! write to ab_out
  write(message,'(4a)')ch10,' Using input gkkq matrix elements read from file ',trim(fname)
  call wrtout(ab_out,message,'COLL')

! DEBUG
! 
! NOTE:
! Re-applying complete_gkk changes nothing in summed-over-bands case
! but smears the linewidths in the keepbands case (probably iterative
! numerical errors). Restrict to one call to complete_gkk, at which point
! summed-over-bands and keepbands agree perfectly.
! 
! if (elph_ds%tsymgkq ==1) then
! !------------------------------------------------------
! ! complete gkk matrices for other qpoints on the full grid
! !  inspired and cannibalized from symdm9.f
! !------------------------------------------------------
! write(*,*)' get_all_gkq : entering 2nd complete_gkk'
! call complete_gkk(acell,elph_ds,FSfulltofull,FSkpt,gkk_flag,&
! &    gprimd,indsym,mpert,natom,nqptirred,nsym,&
! &    qptirred(:,1:nqptirred),qpttoqpt,rprimd,spqpt,symrec,symrel,tnons,&
! &    ucvol,xred)
! write (*,*) 'get_all_gkq : out of 2nd complete_gkk'
! end if
! ENDDEBUG

! ========================================================================
! error: gkqwrite must be 1 or 0
! ========================================================================

 else

  write (message,'(4a,i4)')ch10,' get_all_gkq : BUG-',&
&  ch10,' bad value for gkqwrite = ',elph_ds%gkqwrite
  call wrtout(06,message,'COLL')
  call leave_new('COLL')

 end if

!DEBUG
!print out gkq elements for gamma
!if (elph_ds%gkqwrite == 2) then
!if (elph_ds%gkqwrite == 1) then
!!allocate (gkk_tmp_full(2,elph_ds%nbranch,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nFSkpt))
!allocate (gkk_tmp_full(elph_ds%ngkkband,elph_ds%ngkkband,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt))
!do iqpt=1,elph_ds%nqpt
!!if (spqpt(1,iqpt)**2 + spqpt(2,iqpt)**2 +spqpt(3,iqpt)**2  < tol6) exit
!!write (105,*) ' Print out Gamma point gkq for all kpoints'
!read (elph_ds%unitgkq,REC=iqpt) gkk_tmp_full
!!do ib1=1,elph_ds%nFSband
!!   do ib2=1,elph_ds%nFSband
!!      do iFSkpt=1,elph_ds%nFSkpt
!!         write (105,'(4I5,3E10.2,3(2E16.6,2x))') &
!!              &   ib2,ib1,iFSkpt,iqpt,FSkpt(:,iFSkpt),&
!!              &   gkk_tmp_full(:,:,ib2,ib1,iFSkpt)
!!      end do
!!   end do
!!end do
!do iFSkpt=1,elph_ds%nFSkpt
!write (105,'(2I5,3E10.2,3(3E16.6,2x))') &
!&   iFSkpt,iqpt,FSkpt(:,iFSkpt),&
!&   gkk_tmp_full(:,:,1,1,iFSkpt)
!end do
!end do
!
!deallocate (gkk_tmp_full)
!end if
!stop
!ENDDEBUG

!end if
!end rpt exists already if

end subroutine get_all_gkq
!!***
