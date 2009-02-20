!{\src2tex{textfont=tt}}
!!****f* ABINIT/completeperts
!!
!! NAME
!! completeperts
!!
!! FUNCTION
!!  Complete perturbations wrt atoms and reduced directions
!!  for a fixed qpoint. Normally there is a test in read_gkk which guarantees
!!  that enough irreducible perturbations are present to generate everything.
!!  h1_mat_el is first squared, making a (ipert,jpert) matrix which has the same
!!  symmetry properties as the dynamical matrix.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = datastructure for elph data (dimensions and eventually data)
!!   FSfulltofull = mapping btw kpoints under symops
!!   gkk_flag = flags for presence of gkk matrix elements
!!   h1_mat_el = irreducible matrix elements to be completed and squared
!!   hdr1 = header of files for present perturbation
!!   indsym = mapping of atoms under symops
!!   iqptfull = qpoint number in full zone
!!   irredpert = symops to re-constitute perturbations
!!   natom = number of atoms
!!   nsym = number of syms
!!   spqpt = qpoint (a priori full grid)
!!   symq = flags for symmetry elements conserving the present qpoint
!!   symrec = symmetry operations for reduced reciprocal coordinates $=(symrel^{-1})^{T}$
!!   symrel = symmetry operations for reduced reciprocal coordinates
!!   tnons = translation vectors associated with symops
!!
!! OUTPUT
!!   h1_mat_el_sq = irreducible matrix elements squared and completed
!!   gkk_flag = changed on output
!!
!! NOTES
!!
!! PARENTS
!!      read_gkk
!!
!! CHILDREN
!!      d2sym3
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine completeperts(elph_ds,FSfulltofull,gkk_flag,h1_mat_el,h1_mat_el_sq,hdr1,&
&   indsym,iqptfull,irredpert,natom,nsym,spqpt,symq,symrec,symrel,timrev,tnons)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_16response
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptfull,natom,nsym,timrev
 type(elph_type),intent(in) :: elph_ds
 type(hdr_type),intent(in) :: hdr1
!arrays
 integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt),indsym(4,nsym,natom)
 integer,intent(in) :: irredpert(7,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nqpt)
 integer,intent(in) :: symq(4,2,nsym),symrec(3,3,nsym),symrel(3,3,nsym)
 integer,intent(inout) :: gkk_flag(elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol,elph_ds%nqpt)
 real(dp),intent(in) :: h1_mat_el(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
 real(dp),intent(in) :: spqpt(3,elph_ds%nqpt),tnons(3,nsym)
 real(dp),intent(out) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,iatom1,iatom2,ibb,idir1,idir2,ipert1,ipert2,ipreatom,ipredir
 integer :: isppol,isym,isymatom,isymdir,itim,mpert,prepert,printflag,symiFSkpt
 real(dp) :: exparg,im1,im2,re1,re2,res,s1,s2,s3,ss,sumi,sumr,timsign,valtol
 character(len=500) :: message
!arrays
 integer :: symmetrized(elph_ds%nbranch)
 integer,allocatable :: tmpflg(:,:,:,:)
 real(dp) :: cosarg(2,nsym),dsymrec(3,3),qpt(3),sinarg(2,nsym)
 real(dp),allocatable :: tmpval(:,:,:,:,:)

! *************************************************************************

!WARNING! Stupid patch in d2sym3 imposes these matrices to have size natom+2
 mpert = natom+2

 allocate(tmpflg(3,mpert,3,mpert))
 allocate(tmpval(2,3,mpert,3,mpert))

 valtol = 1.0d-50
 printflag = 0

!DEBUG
!write (*,*) 'gkk_flag (:,:,1,1,iqptfull) = ',  gkk_flag(:,:,1,1,iqptfull)
!ENDDEBUG

 qpt(:) = spqpt(:,iqptfull)

 h1_mat_el_sq(:,:,:,:,:) = zero
 write (*,*) ' completeperts: shape(h1_mat_el_sq) = ', shape(h1_mat_el_sq)

 do isppol=1,elph_ds%nsppol
  write(*,*)'completeperts: isppol = ', isppol

  do iFSkpt=1,elph_ds%nFSkpt
!  write(*,*)'completeperts: gkk_flag = ', gkk_flag(:,:,iFSkpt,isppol,iqptfull)
   do ibb=1,elph_ds%nFSband*elph_ds%nFSband

    tmpval(:,:,:,:,:) = zero
    tmpflg(:,:,:,:) = 0
    do iatom1=1,natom
     do idir1=1,3
      ipert1 = (iatom1-1)*3+idir1
      if (gkk_flag(ipert1,ipert1,iFSkpt,isppol,iqptfull) < 0) cycle
      re1 = h1_mat_el(1,ibb,ipert1,iFSkpt,isppol)
      im1 = h1_mat_el(2,ibb,ipert1,iFSkpt,isppol)

      do iatom2=1,natom
       do idir2=1,3
        ipert2 = (iatom2-1)*3+idir2
        if (gkk_flag(ipert2,ipert2,iFSkpt,isppol,iqptfull) < 0) cycle
!       DEBUG
!       write(*,*)' we have already element ',ipert1,ipert2,&
!       &gkk_flag(ipert1,ipert1,iFSkpt,isppol,iqptfull),&
!       &gkk_flag(ipert2,ipert2,iFSkpt,isppol,iqptfull)
!       ENDDEBUG
        tmpflg(idir1,iatom1,idir2,iatom2) = 1
        re2 = h1_mat_el(1,ibb,ipert2,iFSkpt,isppol)
        im2 = h1_mat_el(2,ibb,ipert2,iFSkpt,isppol)
!       write (*,*) idir1,iatom1,idir2,iatom2,re1,im1,re2,im2

!       conjg(h1_mat_el_2) * h1_mat_el_1
        res =  re1*re2 + im1*im2
!       if (abs(res) > valtol) then
        tmpval(1,idir1,iatom1,idir2,iatom2) =  res
!       end if
        res =  re1*im2 - im1*re2
!       if (abs(res) > valtol) then
        tmpval(2,idir1,iatom1,idir2,iatom2) = res
!       end if
!       write(*,*) idir1,iatom1,idir2,iatom2,tmpval(1,idir1,iatom1,idir2,iatom2),tmpval(2,idir1,iatom1,idir2,iatom2)

       end do !idir2 
      end do !iatom2
     end do !idir1
    end do !iatom1

!   DEBUG
!   if (abs(tmpval(1,1,1,1,1)) > 1.0d-4 .and. printflag < 10) then
!   write(*,*)' BEFORE '
!   write(*,*)'flag ', tmpflg(:,1:natom,:,1:natom)
!   write(*,*)'real ', tmpval(1,:,1:natom,:,1:natom)
!   write(*,*)'imag ', tmpval(2,:,1:natom,:,1:natom)
!   end if
!   
!   if (abs(tmpval(1,1,1,1,1)) > 1.0d-4 .and. printflag < 1) then
!   write(*,*)' completeperts : call d2sym3'
!   write(*,*)"tmpflg"
!   write(*,*)tmpflg
!   write(*,*)"tmpval"
!   write(*,*)mpval
!   write(*,*)"indsym"
!   write(*,*)indsym
!   write(*,*)"mpert"
!   write(*,*)mpert
!   write(*,*)"natom"
!   write(*,*)natom
!   write(*,*)"nsym"
!   write(*,*)nsym
!   write(*,*)"qpt"
!   write(*,*)qpt
!   write(*,*)"symq"
!   write(*,*)symq(:,:,1:nsym)
!   write(*,*)"symrec"
!   write(*,*)symrec
!   write(*,*)"symrel"
!   write(*,*)symrel
!   write(*,*)"timrev"
!   write(*,*)timrev
!   end if
!   !Then apply symmetry operations
!   write (*,*) ' before d2sym3 tmpflg = ', tmpflg(:,1:natom,:,1:natom)
!   ENDDEBUG

    call d2sym3(tmpflg,tmpval,indsym,mpert,natom,nsym,qpt,symq,symrec,symrel,timrev)

    if (sum(tmpflg(:,1:natom,:,1:natom)) /= 3*natom*3*natom) then
     write(*,*)' tmpflg = ',tmpflg
     write(*,*)'  iFSkpt,isppol,iqptfull', iFSkpt,isppol,iqptfull
     write(message,'(4a)')ch10,&
&     ' completeperts : ERROR- ',ch10,&
&     ' A perturbation is missing after completion with d2sym3'
     call wrtout(06,message,'COLL')
     call leave_new('COLL')
    end if

!   DEBUG
!   if (abs(tmpval(1,1,1,1,1)) > 1.0d-4 .and. printflag < 10) then
!   write(*,*)' AFTER '
!   write(*,*)'flag ', tmpflg(:,1:natom,:,1:natom)
!   write(*,*)'real ', tmpval(1,:,1:natom,:,1:natom)
!   write(*,*)'imag ', tmpval(2,:,1:natom,:,1:natom)
!   printflag = printflag + 1
!   end if
!   ENDDEBUG

!   Save values for calculation of |gkk|^2
    do iatom1=1,natom
     do idir1=1,3
      ipert1 = (iatom1-1)*3+idir1
      do iatom2=1,natom
       do idir2=1,3
        
!       mjv 29/10/2007 ipert2 now contains the composite index ip1*nperts+ip2
        ipert2 = (iatom2-1)*3+idir2 + (ipert1-1)*3*natom
        h1_mat_el_sq(1,ibb,ipert2,iFSkpt,isppol) = tmpval(1,idir2,iatom2,idir1,iatom1)
        h1_mat_el_sq(2,ibb,ipert2,iFSkpt,isppol) = tmpval(2,idir2,iatom2,idir1,iatom1)
        
       end do
      end do
     end do
    end do

   end do !end ibb band dos

!  Set flags
   do ipert1=1,3*natom
    do ipert2=1,3*natom
     if (gkk_flag(ipert2,ipert1,iFSkpt,isppol,iqptfull) < 0) then
      gkk_flag(ipert2,ipert1,iFSkpt,isppol,iqptfull) = 1
     end if
    end do
   end do

  end do !end FSkpt do
 end do !end sppol do

 deallocate(tmpflg,tmpval)

end subroutine completeperts
!!***
