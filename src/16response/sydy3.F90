!{\src2tex{textfont=tt}}
!!****f* ABINIT/sydy3
!!
!! NAME
!! sydy3
!!
!! FUNCTION
!! Symmetrize dynamical matrix (diagonal wrt to the atoms)
!! Unsymmetrized dynamical matrix   is  input as dyfrow;
!! symmetrized dynamical matrix is then  placed in sdyfro.
!! If nsym=1 simply copy dyfrow   into sdyfro.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dyfrow(3,3,natom)=unsymmetrized dynamical matrix
!!  indsym(4,nsym,natom)=label given by subroutine symatm.
!!  natom=number of atoms in cell.
!!  nsym=number of symmetry operators in group.
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on real
!!    space primitive translations (integers).
!!
!! OUTPUT
!!  sdyfro(3,3,natom)=symmetrized dynamical matrix
!!
!! NOTES
!! Symmetrization of gradients with respect to reduced
!! coordinates tn is conducted according to the expression
!! $[d(e)/d(t(n,a))]_{symmetrized} = (1/Nsym)*Sum(S)*symrec(n,m,S)*
!!              [d(e)/d(t(m,b))]_{unsymmetrized}$
!! where $t(m,b)= (symrel^{-1})(m,n)*(t(n,a)-tnons(n))$ and tnons
!! is a possible nonsymmorphic translation.  The label "b" here
!! refers to the atom which gets rotated into "a" under symmetry "S".
!! symrel is the symmetry matrix in real space, which is the inverse
!! transpose of symrec.  symrec is the symmetry matrix in reciprocal
!! space.  $sym_{cartesian} = R * symrel * R^{-1} = G * symrec * G^{-1}$
!! where the columns of R and G are the dimensional primitive translations
!! in real and reciprocal space respectively.
!! Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      dyfro3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine sydy3(dyfrow,indsym,natom,nsym,sdyfro,symrec)

 use defs_basis

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym)
 real(dp),intent(in) :: dyfrow(3,3,natom)
 real(dp),intent(out) :: sdyfro(3,3,natom)

!Local variables -------------------------
!scalars
 integer :: ia,im,ind,isym,kappa,mu,nu,re
 real(dp) :: div
!arrays
 real(dp) :: work(3,3)

! *********************************************************************

 if (nsym==1) then

! Only symmetry is identity so simply copy
  sdyfro(:,:,:)=dyfrow(:,:,:)

 else

! Actually carry out symmetrization
  sdyfro(:,:,:)=0.0_dp
  do ia=1,natom
   do isym=1,nsym
    ind=indsym(4,isym,ia)
    work(:,:)=0.0_dp
    do mu=1,3
     do nu=1,3
      do kappa=1,3
       work(mu,kappa)=work(mu,kappa)+&
&       symrec(mu,nu,isym)*dyfrow(nu,kappa,ind)
      end do
     end do
    end do
    do mu=1,3
     do nu=1,3
      do kappa=1,3
       sdyfro(kappa,mu,ia)=sdyfro(kappa,mu,ia)+&
&       symrec(mu,nu,isym)*work(kappa,nu)
      end do
     end do
    end do
   end do
  end do
  div=1.0_dp/dble(nsym)
  sdyfro(:,:,:)=div*sdyfro(:,:,:)

 end if

end subroutine sydy3
!!***
