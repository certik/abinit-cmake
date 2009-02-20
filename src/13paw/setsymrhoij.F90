!{\src2tex{textfont=tt}}
!!****f* ABINIT/setsymrhoij
!! NAME
!! setsymrhoij
!!
!! FUNCTION
!! PAW only
!! Compute coefficients used later to symmetrize rhoij quantities (augmentation occupancies)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (NH, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gprimd(3,3)==dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  lmax=value of lmax mentioned at the second line of the psp file
!!  nsym=number of symmetry elements in space group
!!  pawprtvol=control print volume and debugging output for PAW
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!
!! OUTPUT
!!  zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)=coefficients of the
!!      transformation of real spherical harmonics
!!      under the symmetry operations
!!
!! NOTES
!!  - This file comes from the file crystal_symmetry.f
!!    by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!  - Uses sign & phase convension of  M. E. Rose, Elementary Theory of Angular
!!    Momentum, John Wiley & Sons,. inc. 1957)
!!    zalpha = exp(-i*alpha)   zgamma = exp (-i*gamma)
!!  - Assumes each transformation  can be expressed in terms of 3 Euler
!!    angles with or without inversion
!!
!!  - In case of antiferromagnetism, zarot coeffs contain (-1)^l factors
!!    from spin inversion
!!
!! PARENTS
!!      gstate,respfn,screening,sigma
!!
!! CHILDREN
!!      mkeuler
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine setsymrhoij(gprimd,lmax,nsym,pawprtvol,rprimd,symafm,symrec,zarot)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_13paw, except_this_one => setsymrhoij
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lmax,nsym,pawprtvol
!arrays
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rprimd(3,3)
 real(dp),intent(out) :: zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)

!Local variables ------------------------------
!scalars
 integer :: i1,ii,il,irot,isn,j1,jj,k1,l1,l2,l3,ll,mm,mp
 real(dp) :: cosalp,cosbeta,cosgam,sinalp,singam
 character(len=500) :: message
!arrays
 real(dp) :: prod(3,3),rot(3,3)
!************************************************************************
 if (abs(pawprtvol)>=3) then
  print *,ch10," PAW TEST:"
  print *,' ==== setsymrhoij: symmetry elements ============'
  Print *,'  >Number of symmetries (nsym)=',nsym
 end if

 zarot=zero

 do irot=1,nsym

  if (abs(pawprtvol)>=3) print '(a,i2,a,9i2,a)',&
&  '   >For symmetry ',irot,' (',symrec(:,:,irot),')'

! === l=0 case ===
  zarot(1,1,1,irot)=1._dp

! === l>0 case ===
  if (lmax>0) then
!  Calculate the rotations in the cartesian basis
   rot=zero;prod=zero
   do k1=1,3
    do j1=1,3
     do i1=1,3
      prod(i1,j1)=prod(i1,j1)+symrec(i1,k1,irot)*rprimd(j1,k1)
     end do
    end do
   end do
   do j1=1,3
    do i1=1,3
     do k1=1,3
      rot(i1,j1)=rot(i1,j1)+gprimd(i1,k1)*prod(k1,j1)
     end do
     if(abs(rot(i1,j1))<tol10) rot(i1,j1)=zero
    end do
   end do
   call mkeuler(rot,cosbeta,cosalp,sinalp,cosgam,singam,isn)
   do ll=1,lmax
!   il=(isn)**ll*(symafm(irot)**ll)
    il=(isn)**ll
    do mm=-ll,ll
     ii=mm+ll+1
     do mp=-ll,ll
      jj=mp+ll+1

!     Formula (47) from the paper of Blanco et al
      zarot(jj,ii,ll+1,irot)=il&
&      *(phim(cosalp,sinalp,mp)*phim(cosgam,singam,mm)*sign(1,mm)&
      *(dbeta(cosbeta,ll,abs(mm),abs(mp))&
&      +(-1._dp)**mp*dbeta(cosbeta,ll,abs(mp),-abs(mm)))/2._dp&
&      -phim(cosalp,sinalp,-mp)*phim(cosgam,singam,-mm)*sign(1,mp)&
      *(dbeta(cosbeta,ll,abs(mm),abs(mp))&
&      -(-1._dp)**mp*dbeta(cosbeta,ll,abs(mp),-abs(mm)))/2._dp)
     end do
    end do
   end do
  end if   ! lmax case

  if (abs(pawprtvol)>=3) then
   if(lmax>0) then
    print *,'    Rotation matrices for l=1:'
    do ii=1,3
     write(6,'(5x,3(f7.3,2x))') (zarot(ii,jj,2,irot),jj=1,3)
    end do
   end if
   if(lmax>1) then
    print *,'    Rotation matrices for l=2:'
    do ii=1,5
     write(6,'(5x,5(f7.3,2x))') (zarot(ii,jj,3,irot),jj=1,5)
    end do
   end if
   if (lmax>2) then
    print *,'    Rotation matrices for l=3:'
    do ii=1,7
     write(6,'(5x,7(f7.3,2x))') (zarot(ii,jj,4,irot),jj=1,7)
    end do
   end if
  end if

 end do  ! isym loop

end subroutine setsymrhoij

!!***
