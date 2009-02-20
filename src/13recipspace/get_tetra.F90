!{\src2tex{textfont=tt}}
!!****f* ABINIT/get_tetra
!! NAME
!! get_tetra
!!
!! FUNCTION
!! get tetrahedra characterized by apexes
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (MVer,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  indkpt(nkpt_fullbz)=indexes of irred kpoints equivalent to kpt_fullbz
!!  gprimd = recip space vectors
!!  klatt(3,3)=reciprocal of lattice vectors for full kpoint grid
!!  kpt_fullbz(3,nkpt_fullbz)=kpoints in full brillouin zone
!!  mtetra=maximum number of tetrahedra
!!  nkpt_fullbz=number of kpoints in full brillouin zone
!!
!! OUTPUT
!!  tetra_full(4,2,mtetra)=for each tetrahedron,
!!     the different instances of the tetrahedron (fullbz kpoints)
!!  tetra_mult(mtetra) = store multiplicity of each irred tetrahedron
!!  tetra_wrap(3,4,mtetra) = store flag to wrap tetrahedron summit into IBZ
!!  ntetra = final number of irred tetrahedra (dimensions of tetra_* remain larger)
!!  vv = tetrahedron volume
!!
!! PARENTS
!!      elphon,mkphdos,tetrahedron
!!
!! CHILDREN
!!      canon9,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine get_tetra (indkpt,gprimd,klatt,kpt_fullbz,mtetra,nkpt_fullbz,&
&                ntetra,tetra_full,tetra_mult,tetra_wrap,vv)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mtetra,nkpt_fullbz
 integer,intent(out) :: ntetra
 real(dp),intent(out) :: vv
!arrays
 integer,intent(in) :: indkpt(nkpt_fullbz)
 integer,intent(out) :: tetra_full(4,2,mtetra),tetra_mult(mtetra)
 integer,intent(out) :: tetra_wrap(3,4,mtetra)
 real(dp),intent(in) :: gprimd(3,3),klatt(3,3),kpt_fullbz(3,nkpt_fullbz)

!Local variables-------------------------------
! 3 dimensions, 4 summits, and 6 tetrahedra / kpoint box
!scalars
 integer :: ialltetra,ikpt2,ikpt_full,isummit,itetra,jalltetra,jsummit
 real(dp) :: itmp,shift1,shift2,shift3
 character(len=500) :: message
!arrays
 integer :: tetra_shifts(3,4,6),tmptetra(4)
 real(dp) :: k1(3),k2(3),k3(3),summit(3)

! *********************************************************************

 tetra_mult(:) = 1
 tetra_full(:,:,:) = 0
 tetra_wrap(:,:,:) = 0

 tetra_shifts(:,1,1) = (/0,0,0/)
 tetra_shifts(:,2,1) = (/0,1,0/)
 tetra_shifts(:,3,1) = (/0,1,1/)
 tetra_shifts(:,4,1) = (/1,1,0/)
 tetra_shifts(:,1,2) = (/0,0,0/)
 tetra_shifts(:,2,2) = (/0,1,1/)
 tetra_shifts(:,3,2) = (/1,1,0/)
 tetra_shifts(:,4,2) = (/1,1,1/)
 tetra_shifts(:,1,3) = (/0,0,0/)
 tetra_shifts(:,2,3) = (/1,0,0/)
 tetra_shifts(:,3,3) = (/1,1,0/)
 tetra_shifts(:,4,3) = (/1,1,1/)
 tetra_shifts(:,1,4) = (/0,0,0/)
 tetra_shifts(:,2,4) = (/0,0,1/)
 tetra_shifts(:,3,4) = (/1,0,0/)
 tetra_shifts(:,4,4) = (/1,1,1/)
 tetra_shifts(:,1,5) = (/0,0,1/)
 tetra_shifts(:,2,5) = (/1,0,0/)
 tetra_shifts(:,3,5) = (/1,0,1/)
 tetra_shifts(:,4,5) = (/1,1,1/)
 tetra_shifts(:,1,6) = (/0,0,0/)
 tetra_shifts(:,2,6) = (/0,0,1/)
 tetra_shifts(:,3,6) = (/0,1,1/)
 tetra_shifts(:,4,6) = (/1,1,1/)

 ialltetra = 1
 do ikpt_full=1,nkpt_fullbz
  do itetra=1,6
   do isummit=1,4

    k1(:) = kpt_fullbz(:,ikpt_full) &
&    + tetra_shifts(1,isummit,itetra)*klatt(:,1) &
&    + tetra_shifts(2,isummit,itetra)*klatt(:,2) &
&    + tetra_shifts(3,isummit,itetra)*klatt(:,3)

!   Wrap the trial values in the interval ]-1/2,1/2] .
    call canon9(k1(1),summit(1),shift1)
    call canon9(k1(2),summit(2),shift2)
    call canon9(k1(3),summit(3),shift3)

!   Find full kpoint which is summit isummit of tetrahedron itetra around full kpt ikpt_full !
    do ikpt2=1,nkpt_fullbz
     if ((abs(summit(1) - kpt_fullbz(1,ikpt2)) &
&     + abs(summit(2) - kpt_fullbz(2,ikpt2)) &
&     + abs(summit(3) - kpt_fullbz(3,ikpt2))) < tol6) then
!     Store irreducible kpoint equivalent to kpt_fullbz(:,ikpt2)
      tetra_full(isummit,1,ialltetra) = indkpt(ikpt2)
      tetra_full(isummit,2,ialltetra) = ikpt2
      if (shift1>half) then
       tetra_wrap(1,isummit,ialltetra) = 1
      else if (shift1<-half) then
       tetra_wrap(1,isummit,ialltetra) = -1
      end if
      if (shift2>half) then
       tetra_wrap(2,isummit,ialltetra) = 1
      else if (shift2<-half) then
       tetra_wrap(2,isummit,ialltetra) = -1
      end if
      if (shift3>half) then
       tetra_wrap(3,isummit,ialltetra) = 1
      else if (shift3<-half) then
       tetra_wrap(3,isummit,ialltetra) = -1
      end if
      exit
     end if
    end do !  loop ikpt2
!   sort itetra summits
    do jsummit=isummit,2,-1
     if ( tetra_full(jsummit,1,ialltetra) &
&     <  tetra_full(jsummit-1,1,ialltetra) ) then
      itmp = tetra_full(jsummit,1,ialltetra)
      tetra_full(jsummit,1,ialltetra) = tetra_full(jsummit-1,1,ialltetra)
      tetra_full(jsummit-1,1,ialltetra) = itmp
      itmp = tetra_full(jsummit,2,ialltetra)
      tetra_full(jsummit,2,ialltetra) = tetra_full(jsummit-1,2,ialltetra)
      tetra_full(jsummit-1,2,ialltetra) = itmp
!     keep fullbz_kpt tetrahedra points in same order
      itmp = tetra_wrap(1,jsummit,ialltetra)
      tetra_wrap(1,jsummit,ialltetra) = tetra_wrap(1,jsummit-1,ialltetra)
      tetra_wrap(1,jsummit-1,ialltetra) = itmp
      itmp = tetra_wrap(2,jsummit,ialltetra)
      tetra_wrap(2,jsummit,ialltetra) = tetra_wrap(2,jsummit-1,ialltetra)
      tetra_wrap(2,jsummit-1,ialltetra) = itmp
      itmp = tetra_wrap(1,jsummit,ialltetra)
      tetra_wrap(3,jsummit,ialltetra) = tetra_wrap(3,jsummit-1,ialltetra)
      tetra_wrap(3,jsummit-1,ialltetra) = itmp
     end if
    end do !  loop jsummit

   end do !  loop isummit

!  DEBUG
!  write (*,*) 'tetra_full(:,:,',ialltetra,') = ',tetra_full(:,:,ialltetra)
!  ENDDEBUG

!  DEBUG
!  write (*,*) 'get_tetra: ialltetra, tetra_full = ',&
!  &                 ialltetra, tetra_full(:,:,:,ialltetra)
!  ENDDEBUG

   if (ialltetra > mtetra) then
    write (message, '(6a,i6,a,i6)' ) ch10,&
&    ' get_tetra: BUG -',ch10,&
&    '  ialltetra > mtetra',ch10,&
&    '  ialltetra=',ialltetra,', mtetra=',mtetra
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
   ialltetra = ialltetra+1
  end do !  loop itetra
 end do !  loop ikpt_full

!Volume of all tetrahedra should be the same as that of tetra 1
!this is the volume of 1 tetrahedron, should be coherent with
!notation in Lehmann & Taut
 k1(:) = gprimd(:,1)*klatt(1,1) &
& +  gprimd(:,2)*klatt(2,1) &
& +  gprimd(:,3)*klatt(3,1)
 k2(:) = gprimd(:,1)*klatt(1,2) &
& +  gprimd(:,2)*klatt(2,2) &
& +  gprimd(:,3)*klatt(3,2)
 k3(:) = gprimd(:,1)*klatt(1,3) &
& +  gprimd(:,2)*klatt(2,3) &
& +  gprimd(:,3)*klatt(3,3)
 vv  = abs (k1(1)*(k2(2)*k3(3)-k2(3)*k3(2)) &
& -k1(2)*(k2(1)*k3(3)-k2(3)*k3(1)) &
& +k1(3)*(k2(1)*k3(2)-k2(2)*k3(1))) / six
!DEBUG
!write (*,*) 'get_tetra : '
!write (*,*) 'gprimd = ', gprimd(:,1)
!write (*,*) '       ',   gprimd(:,2)
!write (*,*) '       ',   gprimd(:,3)
!write (*,*) 'k1     ', k1(:)
!write (*,*) 'k2     ', k2(:)
!write (*,*) 'k3     ', k3(:)
!write (*,*) 'Tetrahedron volume is (bohr^3): ', vv
!ENDDEBUG

!
!eliminate equivalent tetrahedra by symmetry and account
!for them in multiplicity tetra_mult
!
 ntetra = mtetra
 do ialltetra=ntetra,2,-1
  do jalltetra=1,ialltetra-1
!  check if tetra are equivalent
   if (tetra_full(1,1,ialltetra) == tetra_full(1,1,jalltetra) .and. &
&   tetra_full(2,1,ialltetra) == tetra_full(2,1,jalltetra) .and. &
&   tetra_full(3,1,ialltetra) == tetra_full(3,1,jalltetra) .and. &
&   tetra_full(4,1,ialltetra) == tetra_full(4,1,jalltetra) ) then
!   accumulate multiplicity and positions into equiv tetrahedron jalltetra
    if ( tetra_mult(ialltetra) > 1) then
     write (*,*) 'found an equiv tetra with mult > 1 :', tetra_mult(ialltetra)
     write (*,*) '  ialltetra,jalltetra,ntetra = ',&
&     ialltetra,jalltetra,ntetra
    end if
    tetra_mult(jalltetra) = tetra_mult(jalltetra) + tetra_mult(ialltetra)
    tetra_mult(ialltetra) = 0
    ntetra = ntetra-1
    exit
   end if
  end do ! do jalltetra
 end do ! do ialltetra

!
!pack irred tetrahedra
!
 do ialltetra=mtetra,1,-1
  if (tetra_mult(ialltetra) /= 0) then
!  look for an earlier place to put the irred tetrahedron
   do jalltetra=1,ialltetra-1
    if (tetra_mult(jalltetra) == 0) then
!    
!    swap tetrahedrons
!    
     tmptetra(:) = tetra_full(:,1,jalltetra)
     tetra_full(:,1,jalltetra) = tetra_full(:,1,ialltetra)
     tetra_full(:,1,ialltetra) = tmptetra(:)
     tmptetra(:) = tetra_full(:,2,jalltetra)
     tetra_full(:,2,jalltetra) = tetra_full(:,2,ialltetra)
     tetra_full(:,2,ialltetra) = tmptetra(:)
!    
!    swap wrap flags
!    
     tmptetra(:) = tetra_wrap(1,:,jalltetra)
     tetra_wrap(1,:,jalltetra) = tetra_wrap(1,:,ialltetra)
     tetra_wrap(1,:,ialltetra) = tmptetra(:)
     tmptetra(:) = tetra_wrap(2,:,jalltetra)
     tetra_wrap(2,:,jalltetra) = tetra_wrap(2,:,ialltetra)
     tetra_wrap(2,:,ialltetra) = tmptetra(:)
     tmptetra(:) = tetra_wrap(3,:,jalltetra)
     tetra_wrap(3,:,jalltetra) = tetra_wrap(3,:,ialltetra)
     tetra_wrap(3,:,ialltetra) = tmptetra(:)
!    
!    swap multiplicities (tetra_mult(jalltetra) was 0)
!    
     tetra_mult(jalltetra) = tetra_mult(ialltetra)
     if (tetra_mult(ialltetra) > ntetra) then
      write (*,*) 'problem : multiplicity > ntetra ', &
&      tetra_mult(ialltetra),ntetra
     end if
     tetra_mult(ialltetra) = 0
     exit
    end if
   end do ! do jalltetra
  end if
 end do ! do ialltetra

!DEBUG
!write (*,*) 'max num of tetrahedra = ', nkpt_fullbz*6
!write (*,*) 'mtetra = ', mtetra
!do itetra=1,mtetra
!if(sum(abs(tetra_wrap(:,:,itetra))) > 0) then
!write (*,*) 'tetra_full = '
!write (*,'(4(I7,1x))') tetra_full(:,:,itetra)
!end if
!end do
!write (*,*) 'tetra_full = '
!do itetra=1,mtetra
!write (*,'(4(I7,1x))') tetra_full(:,:,itetra)
!write (*,*)
!end do
!write (*,*) 'tetra_wrap = '
!do itetra=1,mtetra
!write (*,'(12(I4,1x))' ) tetra_wrap(:,:,itetra)
!write (*,*)
!end do
!write (*,*) 'tetra_mult = '
!write (*,'(8(I7,1x))' ) (tetra_mult(itetra),itetra=1,mtetra)
!write (*,*) 'ntetra/max num of tetrahedra = ', &
!&             ntetra/nkpt_fullbz/6
!ENDDEBUG

end subroutine get_tetra
!!***
