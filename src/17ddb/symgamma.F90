!{\src2tex{textfont=tt}}
!!****f* ABINIT/symgamma
!!
!! NAME
!! symgamma
!!
!! FUNCTION
!!  Symmetrize perturbations wrt atoms and reduced directions
!!  for a fixed qpoint.
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
!!   h1_mat_el = irreducible matrix elements to be completed
!!   indsym = mapping of atoms under symops
!!   iqptfull = qpoint number in full zone
!!   natom = number of atoms
!!   nsym = number of syms
!!   symq = flags for symmetry elements conserving the present qpoint
!!   symrec = symmetry operations in reciprocal space
!!
!! OUTPUT
!!   h1_mat_el = changed on output
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      int2char4
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symgamma(elph_ds,FSfulltofull,h1_mat_el,&
&   indsym,iqptfull,natom,nsym,symq,symrec)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptfull,natom,nsym
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: FSfulltofull(2,nsym,elph_ds%nFSkpt),indsym(4,nsym,natom)
 integer,intent(in) :: symq(4,2,nsym),symrec(3,3,nsym)
 real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nbranch,elph_ds%nFSkpt)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,ib1,ib2,ipreatom,ipredir,isym,isymatom,isymdir,itim
 integer :: nsymperts,prepert,symiFSkpt,sympertcase
 real(dp) :: exparg,s1,s2,s3,ss,sumi,sumr,timsign
 character(len=4) :: isymstr,itimstr
 character(len=fnlen) :: xsffilnam
!arrays
 integer :: symmetrized(elph_ds%nbranch)
 real(dp) :: cosarg(2,nsym),dsymrec(3,3),sinarg(2,nsym)
 real(dp) :: sym_mat_el(2,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nbranch,elph_ds%nFSkpt)

! *************************************************************************

!2. symmetrize the whole set of perturbations
 sym_mat_el(:,:,:,:,:) = zero

 nsymperts = 0
!
!symrel(isym) sends ipreatom onto isymatom
!symrec(isym) sends ipredir onto isymdir
!symrec(isym) sends iFSkpt onto symiFSkpt
!
 do isym=1,nsym
! nsymperts = 1
! sym_mat_el(:,:,:,:,:) = h1_mat_el(:,:,:,:,:)
! do isym=27,30

  do itim=0,1
!  do itim=0,0
   if (symq(4,itim+1,isym) == 0) cycle
   timsign = one - two*itim

   dsymrec(:,:) = timsign*dble(symrec(:,:,isym))

   nsymperts = nsymperts+1

!  loop over image perturbations
   do isymatom=1,natom
    ipreatom = indsym(4,isym,isymatom)
    write (*,*) ' symgamma : ', isym,itim,isymatom,ipreatom,dsymrec


    do symiFSkpt=1,elph_ds%nFSkpt
     iFSkpt = FSfulltofull(itim+1,isym,symiFSkpt)

!    Real part
     sym_mat_el (1,:,:,3*(isymatom-1)+1,symiFSkpt) = &
&     sym_mat_el (1,:,:,3*(isymatom-1)+1,symiFSkpt) &
&     + dsymrec(1,1)*h1_mat_el(1,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     + dsymrec(1,2)*h1_mat_el(1,:,:,3*(ipreatom-1)+2,iFSkpt) &
&     + dsymrec(1,3)*h1_mat_el(1,:,:,3*(ipreatom-1)+3,iFSkpt)
     sym_mat_el (1,:,:,3*(isymatom-1)+2,symiFSkpt) = &
&     sym_mat_el (1,:,:,3*(isymatom-1)+2,symiFSkpt) &
&     + dsymrec(2,1)*h1_mat_el(1,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     + dsymrec(2,2)*h1_mat_el(1,:,:,3*(ipreatom-1)+2,iFSkpt) &
&     + dsymrec(2,3)*h1_mat_el(1,:,:,3*(ipreatom-1)+3,iFSkpt)
     sym_mat_el (1,:,:,3*(isymatom-1)+3,symiFSkpt) = &
&     sym_mat_el (1,:,:,3*(isymatom-1)+3,symiFSkpt) &
&     + dsymrec(3,1)*h1_mat_el(1,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     + dsymrec(3,2)*h1_mat_el(1,:,:,3*(ipreatom-1)+2,iFSkpt) &
&     + dsymrec(3,3)*h1_mat_el(1,:,:,3*(ipreatom-1)+3,iFSkpt)
!    Imag part
     sym_mat_el (2,:,:,3*(isymatom-1)+1,symiFSkpt) = &
&     sym_mat_el (2,:,:,3*(isymatom-1)+1,symiFSkpt) &
!    &   + timsign*(dsymrec(1,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     +         (dsymrec(1,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     +          dsymrec(1,2)*h1_mat_el(2,:,:,3*(ipreatom-1)+2,iFSkpt) &
&     +          dsymrec(1,3)*h1_mat_el(2,:,:,3*(ipreatom-1)+3,iFSkpt))
     sym_mat_el (2,:,:,3*(isymatom-1)+2,symiFSkpt) = &
&     sym_mat_el (2,:,:,3*(isymatom-1)+2,symiFSkpt) &
!    &   + timsign*(dsymrec(2,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     +         (dsymrec(2,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     +          dsymrec(2,2)*h1_mat_el(2,:,:,3*(ipreatom-1)+2,iFSkpt) &
&     +          dsymrec(2,3)*h1_mat_el(2,:,:,3*(ipreatom-1)+3,iFSkpt))
     sym_mat_el (2,:,:,3*(isymatom-1)+3,symiFSkpt) = &
&     sym_mat_el (2,:,:,3*(isymatom-1)+3,symiFSkpt) &
!    &   + timsign*(dsymrec(3,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     +         (dsymrec(3,1)*h1_mat_el(2,:,:,3*(ipreatom-1)+1,iFSkpt) &
&     +          dsymrec(3,2)*h1_mat_el(2,:,:,3*(ipreatom-1)+2,iFSkpt) &
&     +          dsymrec(3,3)*h1_mat_el(2,:,:,3*(ipreatom-1)+3,iFSkpt))
    end do
   end do

!  !DEBUG
   if (iqptfull == 6) then
    call int2char4(isym,isymstr)
    call int2char4(itim,itimstr)
    xsffilnam = trim("gkk_" // isymstr // "_" // itimstr //"_symgam.xsf")
    open (130,file=xsffilnam)
    write (130,'(a)') 'DIM-GROUP'
    write (130,'(a)') ' 3  1'
    write (130,'(a)') ' PRIMVEC'
    write (130,'(a)') ' -5.0  5.0  5.0'
    write (130,'(a)') '  5.0 -5.0  5.0'
    write (130,'(a)') '  5.0  5.0 -5.0'
    write (130,'(a)') ' PRIMCOORD'
    write (130,'(a)') ' 1  1'
    write (130,'(a)') '       13    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00'
    write (130,'(a)') ' ATOMS'
    write (130,'(a)') '       13    0.0000000000E+00    0.0000000000E+00    0.0000000000E+00'
    write (130,'(a)') ' BEGIN_BLOCK_DATAGRID3D'
    write (130,'(a)') ' datagrids'
    write (130,'(a)') ' DATAGRID_3D_DENSITY'
    write (130,'(a)') ' 16 16 16'
    write (130,'(a)') ' 0.0 0.0 0.0'
    write (130,'(a)') ' -5.0  5.0  5.0'
    write (130,'(a)') '  5.0 -5.0  5.0'
    write (130,'(a)') '  5.0  5.0 -5.0'
    write (130,'(6E16.6)') sym_mat_el(1,1,2,1,FSfulltofull(1,isym,:))
    write (130,'(a)') ' END_DATAGRID_3D'
    write (130,'(a)') ' END_BLOCK_DATAGRID3D'
    close(130)
   end if
!  !ENDDEBUG
  end do
 end do
!end isym and itim do

!commented to use un-symmetrized version of h1_mat_el
 h1_mat_el(:,:,:,:,:) = sym_mat_el(:,:,:,:,:) / nsymperts

end subroutine symgamma
!!***
