!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_full_kgrid
!! NAME
!! get_full_kgrid
!!
!! FUNCTION
!! create full grid of kpoints and find equivalent
!! irred ones. Duplicates work in getkgrid, but need all outputs
!! of klatt,kpt_fullbz, and indkpt
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  klatt(3,3)=reciprocal of lattice vectors for full kpoint grid
!!  kpt(3,nkpt)=irreducible kpoints
!!  kptrlatt(3,3)=lattice vectors for full kpoint grid
!!  nkpt=number of irreducible kpoints
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!  nshiftk=number of kpoint grid shifts
!!  nsym=number of symmetries
!!  shiftk(3,nshiftk)=kpoint shifts
!!  symrel(3,3,nsym)=symmetry matrices in real space
!!
!! OUTPUT
!!  indkpt(nkpt_fullbz)=non-symmetrized indices of the k-points (see symkpt.f)
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!
!! PARENTS
!!      tetrahedron
!!
!! CHILDREN
!!      canon9,leave_new,mati3inv,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_full_kgrid(indkpt,klatt,kpt,kpt_fullbz,kptrlatt,nkpt,&
& nkpt_fullbz,nshiftk,nsym,shiftk,symrel)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nkpt,nkpt_fullbz,nshiftk,nsym
!arrays
 integer,intent(in) :: kptrlatt(3,3),symrel(3,3,nsym)
 integer,intent(out) :: indkpt(nkpt_fullbz)
 real(dp),intent(in) :: klatt(3,3),kpt(3,nkpt),shiftk(3,nshiftk)
 real(dp),intent(out) :: kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
!scalars
 integer :: ii,ikpt,ikpt2,ikshft,isym,itim,jj,kk,nn,timrev
 character(len=500) :: message
!arrays
 integer :: boundmax(3),boundmin(3),inv_symrel(3,3,nsym)
 real(dp) :: k1(3),k2(3),shift(3)

! *********************************************************************

!Invert symrels => gives symrels for kpoints

 do isym=1,nsym
  call mati3inv (symrel(:,:,isym),inv_symrel(:,:,isym))
 end do

!Generate full grid with 1 shift
!Now, klatt contains the three primitive vectors of the k lattice,
!in reduced coordinates. One builds all k vectors that
!are contained in the first Brillouin zone, with coordinates
!in the interval [0,1[ . First generate boundaries of a big box.

 do jj=1,3
! To accomodate the shifts, boundmin starts from -1
! Well, this is not a complete solution ...
  boundmin(jj)=-1
  boundmax(jj)=0
  do ii=1,3
   if(kptrlatt(ii,jj)<0)boundmin(jj)=boundmin(jj)+kptrlatt(ii,jj)
   if(kptrlatt(ii,jj)>0)boundmax(jj)=boundmax(jj)+kptrlatt(ii,jj)
  end do
 end do

!DEBUG
!write (*,*) 'boundmin = ', boundmin
!write (*,*) 'boundmax = ', boundmax
!ENDDEBUG

 nn=1
 do kk=boundmin(3),boundmax(3)
  do jj=boundmin(2),boundmax(2)
   do ii=boundmin(1),boundmax(1)
    do ikshft=1,nshiftk

!    Coordinates of the trial k point with respect to the k primitive lattice
     k1(1)=ii+shiftk(1,ikshft)
     k1(2)=jj+shiftk(2,ikshft)
     k1(3)=kk+shiftk(3,ikshft)

!    Reduced coordinates of the trial k point
     k2(:)=k1(1)*klatt(:,1)+k1(2)*klatt(:,2)+k1(3)*klatt(:,3)
!    DEBUG
!    write(6,*)' k2(:)',k2(:)
!    ENDDEBUG

!    Eliminate the point if outside [0,1[
     if(k2(1)<-tol10)cycle ; if(k2(1)>one-tol10)cycle
     if(k2(2)<-tol10)cycle ; if(k2(2)>one-tol10)cycle
     if(k2(3)<-tol10)cycle ; if(k2(3)>one-tol10)cycle

!    Wrap the trial values in the interval ]-1/2,1/2] .
     call canon9(k2(1),k1(1),shift(1))
     call canon9(k2(2),k1(2),shift(2))
     call canon9(k2(3),k1(3),shift(3))
     if(nn > nkpt_fullbz) then
      write (message, '(4a,i6)' ) ch10,&
&      ' get_full_kgrid: BUG -',ch10,&
&      '  nkpt_fullbz mis-estimated, exceed nn=',nn
      call wrtout(6,message,'COLL')
      call leave_new('COLL')
     end if
     kpt_fullbz(:,nn)=k1(:)
!    DEBUG
!    write (*,*) 'C ', k1(:)*ten
!    ENDDEBUG
     nn=nn+1
    end do
   end do
  end do
 end do
 nn = nn-1
 if(nn /= nkpt_fullbz) then
  write (message, '(4a,i6,3a,i6)' ) ch10,&
&  ' get_full_kgrid: BUG -',ch10,&
&  '  nkpt_fullbz=',nkpt_fullbz,' underestimated',ch10,&
&  '  nn=',nn
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!
!find equivalence to irred kpoints in kpt
!
 indkpt(:) = 0
 timrev=1 ! includes the time inversion symmetry
 do ikpt=1,nkpt_fullbz
  do isym=1,nsym
   do itim=1,(1-2*timrev),-2

    k2(:) = itim*(inv_symrel(:,1,isym)*kpt_fullbz(1,ikpt) + &
&    inv_symrel(:,2,isym)*kpt_fullbz(2,ikpt) + &
&    inv_symrel(:,3,isym)*kpt_fullbz(3,ikpt))
!   Wrap the trial values in the interval ]-1/2,1/2] .
    call canon9(k2(1),k1(1),shift(1))
    call canon9(k2(2),k1(2),shift(2))
    call canon9(k2(3),k1(3),shift(3))
!   DEBUG
!   write (*,*) 'get_full_kgrid: kpt_fullbz, k1, k2 = ',&
!   &        kpt_fullbz(:,ikpt), ", ", k1, ", ", k2
!   ENDDEBUG
    do ikpt2=1,nkpt
     if ( (abs(k1(1)-kpt(1,ikpt2)) + &
&     abs(k1(2)-kpt(2,ikpt2)) + &
&     abs(k1(3)-kpt(3,ikpt2))) < tol6 ) then
      indkpt(ikpt) = ikpt2
!     DEBUG
!     write (*,'(a,I10,x,3(E16.8,1x))') 'Zat', ikpt2, kpt_fullbz(:,ikpt)
!     ENDDEBUG
      exit
     end if
    end do ! loop irred kpoints

    if (indkpt(ikpt) /= 0) exit

   end do ! loop time reversal symmetry
  end do !  loop sym ops

  if (indkpt(ikpt) == 0) then
   write (message, '(6a,i6)' ) ch10,&
&   ' get_full_kgrid: BUG -',ch10,&
&   '  indkpt(ikpt) is still 0',ch10,&
&   '  no irred kpoint is equiv to ikpt ', ikpt
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

 end do !  loop full kpts

end subroutine get_full_kgrid
!!***
