!{\src2tex{textfont=tt}}
!!****f* ABINIT/symfind
!! NAME
!! symfind
!!
!! FUNCTION
!! Symmetry finder.
!! From the symmetries of the Bravais lattice (ptsymrel),
!! select those that leave invariant the system, and generate
!! the corresponding tnons vectors.
!! The algorithm is explained in T.G. Worlton and J.L. Warren,
!! Comp. Phys. Comm. 3, 88 (1972)
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! berryopt = 4: electric field is on -> add the contribution of the
!!                 - \Omega E.P term to the total energy
!!         /= 4: electric field is off
!! efield=cartesian coordinates of the electric field
!! msym=default maximal number of symmetries
!! natom=number of atoms in cell.
!! nptsym=number of point symmetries of the Bravais lattice
!! ptsymrel(3,3,1:msym)= nptsym point-symmetry operations
!!   of the Bravais lattice in real space in terms
!!   of primitive translations.
!! spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!! typat(natom)=integer identifying type of atom.
!! xred(3,natom)=reduced coordinates of atoms in terms of real space
!!   primitive translations
!!
!! OUTPUT
!! nsym=actual number of symmetries
!! symafm(1:msym)=(anti)ferromagnetic part of nsym symmetry operations
!! symrel(3,3,1:msym)= nsym symmetry operations in real space in terms
!!  of primitive translations
!! tnons(3,1:msym)=nonsymmorphic translations for each symmetry (would
!!  be 0 0 0 each for a symmorphic space group)
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symfind(berryopt,efield,jellslab,msym,natom,nptsym,nsym,ptsymrel,&
& spinat,symafm,symrel,tnons,typat,xred)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,jellslab,msym,natom,nptsym
 integer,intent(out) :: nsym
!arrays
 integer,intent(in) :: ptsymrel(3,3,msym),typat(natom)
 integer,intent(out) :: symafm(msym),symrel(3,3,msym)
 real(dp),intent(in) :: efield(3),spinat(3,natom),xred(3,natom)
 real(dp),intent(out) :: tnons(3,msym)

!Local variables-------------------------------
!scalars
 integer :: found3,foundcl,iatom,iatom0,iatom1,iatom2,iatom3,iclass,iclass0,ii
 integer :: isym,jj,kk,natom0,nclass,ntrial,printed,trialafm,trialok
 real(dp) :: spinatcl2,spinatcl20
 logical :: test_sameabscollin,test_samespin
 character(len=500) :: message
!arrays
 integer,allocatable :: class(:,:),natomcl(:),typecl(:)
 real(dp) :: diff(3),efieldrot(3),sxred0(3),symspinat2(3),symxred2(3)
 real(dp) :: trialnons(3)
 real(dp),allocatable :: spinatcl(:,:)

!**************************************************************************

!DEBUG
!write(6,*)' symfind : enter' 
!write(6,*)' symfind : jellslab = ',jellslab 
!write(6,*)'   ptsymrel matrices are :'
!do isym=1,nptsym
!write(6, '(i4,4x,9i4)' )isym,ptsymrel(:,:,isym)
!end do
!write(6,*)' symfind : natom=',natom
!do iatom=1,natom
!write(6,*)'  atom number',iatom
!write(6,*)'   typat   =',typat(iatom)
!write(6,*)'   spinat  =',spinat(:,iatom)
!end do
!ENDDEBUG

!Find the number of classes of atoms (type spinat must be identical,
!spinat might differ by a sign, if aligned with the z direction)
!natomcl(iclass) will contain the number of atoms in the class
!typecl(iclass) will contain the type of the atoms in the class
!spinatcl(1:3,iclass) will contain the spinat of the atoms in the class
!class(1:natomclass(iclass),iclass) will contain the index of the
!atoms belonging to the class
 allocate(class(natom+3,natom),natomcl(natom),typecl(natom),spinatcl(3,natom))

!Initialise with the first atom
 nclass=1
 natomcl(1)=1
 typecl(1)=typat(1)
 spinatcl(:,1)=spinat(:,1)
 class(1,1)=1
 if(natom>1)then
  do iatom=2,natom
!  DEBUG
!  write(6,*)' symfind : examine iatom=',iatom
!  ENDDEBUG
   foundcl=0
   do iclass=1,nclass
!   Compare the typat and spinat of atom iatom with existing ones.
!   Admit either identical spinat, or z-aligned spinat with same
!   absolute magnitude
    if( typat(iatom)==typecl(iclass)) then
     test_samespin=  abs(spinat(1,iatom)-spinatcl(1,iclass))<tol8 .and. &
&     abs(spinat(2,iatom)-spinatcl(2,iclass))<tol8 .and. &
&     abs(spinat(3,iatom)-spinatcl(3,iclass))<tol8
     test_sameabscollin= &
&     abs(spinat(1,iatom))<tol8 .and. abs(spinatcl(1,iclass))<tol8 .and.&
&     abs(spinat(2,iatom))<tol8 .and. abs(spinatcl(2,iclass))<tol8 .and.&
&     abs(abs(spinat(3,iatom))-abs(spinatcl(3,iclass)))<tol8
     if( test_samespin .or. test_sameabscollin ) then
!     DEBUG
!     write(6,*)' symfind : find it belongs to class iclass=',iclass
!     write(6,*)' symfind : spinat(:,iatom)=',spinat(:,iatom)
!     write(6,*)' symfind : spinatcl(:,iclass)=',spinatcl(:,iclass)
!     write(6,*)' symfind : test_samespin,test_sameabscollin=',&
!     &      test_samespin,test_sameabscollin
!     ENDDEBUG
      natomcl(iclass)=natomcl(iclass)+1
      class(natomcl(iclass),iclass)=iatom
      foundcl=1
      exit
     end if
    end if
   end do
!  If no class with these characteristics exist, create one
   if(foundcl==0)then
    nclass=nclass+1
    natomcl(nclass)=1
    typecl(nclass)=typat(iatom)
    spinatcl(:,nclass)=spinat(:,iatom)
    class(1,nclass)=iatom
   end if
  end do
 end if

!DEBUG
!write(6,*)' symfind : found ',nclass,' nclass of atoms'
!do iclass=1,nclass
!write(6,*)'  class number',iclass
!write(6,*)'   natomcl =',natomcl(iclass)
!write(6,*)'   typecl  =',typecl(iclass)
!write(6,*)'   spinatcl=',spinatcl(:,iclass)
!write(6,*)'   class   =',(class(iatom,iclass),iatom=1,natom)
!end do
!ENDDEBUG

!Select the class with the least number of atoms, and non-zero spinat if any
!It is important to select a magnetic class of atom, if any, otherwise
!the determination of the initial (inclusive) set of symmetries takes only
!non-magnetic symmetries, and not both magnetic and non-magnetic ones, see later.
 iclass0=1
 natom0=natomcl(1)
 spinatcl20=spinatcl(1,1)**2+spinatcl(2,1)**2+spinatcl(3,1)**2
 if(nclass>1)then
  do iclass=2,nclass
   spinatcl2=spinatcl(1,iclass)**2+spinatcl(2,iclass)**2+spinatcl(3,iclass)**2
   if( (natomcl(iclass)<natom0 .and. (spinatcl20<tol10 .or. spinatcl2>tol10))  &
&   .or. (spinatcl20<tol10 .and. spinatcl2>tol10)                         )then
    iclass0=iclass
    natom0=natomcl(iclass)
    spinatcl20=spinatcl2
   end if
  end do
 end if

 printed=0

!DEBUG
!write(6,*)' symfind : has selected iclass0=',iclass0
!write(6,*)' #    iatom     xred             spinat '
!do iatom0=1,natomcl(iclass0)
!iatom=class(iatom0,iclass0)
!write(6, '(2i4,6f10.4)' )iatom0,iatom,xred(:,iatom),spinat(:,iatom)
!end do
!ENDDEBUG

!Big loop over each symmetry operation of the Bravais lattice
 nsym=0
 do isym=1,nptsym

! ji: Check whether symmetry operation leaves efield invariant
  if (berryopt==4) then
   efieldrot(:) = ptsymrel(:,1,isym)*efield(1) +  &
&   ptsymrel(:,2,isym)*efield(2) +  &
&   ptsymrel(:,3,isym)*efield(3)
   diff(:)=efield(:)-efieldrot(:)
   if( (diff(1)**2+diff(2)**2+diff(3)**2) > tol8**2 ) cycle
  end if

! jellium slab case:
  if (jellslab/=0) then
!  check whether symmetry operation produce a rotation only in the xy plane
   if( ptsymrel(1,3,isym)/=0 .or. ptsymrel(2,3,isym)/=0 .or. &
&   ptsymrel(3,1,isym)/=0 .or. ptsymrel(3,2,isym)/=0 ) cycle
!  check whether symmetry operation does not change the z
   if( ptsymrel(3,3,isym)/=1 ) cycle
  end if

! Select a tentative set of associated translations
! First compute the symmetric of the first atom in the smallest class,
! using the point symmetry
  iatom0=class(1,iclass0)
  sxred0(:)=ptsymrel(:,1,isym)*xred(1,iatom0)+ &
&  ptsymrel(:,2,isym)*xred(2,iatom0)+ &
&  ptsymrel(:,3,isym)*xred(3,iatom0)

! From the set of possible image, deduce tentative translations,
! and magnetic factor then test whether it send each atom on a symmetric one
  ntrial=0
  do ii=1,natom0
   iatom1=class(ii,iclass0)

!  The tentative translation is found
   trialnons(:)=xred(:,iatom1)-sxred0(:)
   trialafm=1
   if(abs(spinat(3,iatom1)-spinat(3,iatom0))>tol8)trialafm=-1
   if(sum(abs(spinat(:,iatom1)*trialafm-spinat(:,iatom0)))>tol8)then
    write(message,'(6a,3i5)')ch10,&
&    ' symfind : BUG -',ch10,&
&    '  Problem with matching the spin part within a class.',ch10,&
&    '  isym,iatom0,iatom1=',isym,iatom0,iatom1
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
!  jellium slab case: check whether symmetry operation has no translational
!  component along z
   if( jellslab/=0 .and. abs(trialnons(3)) > tol8 ) cycle
   trialok=1

!  DEBUG
!  write(6,*)' isym,trialnons(:),trialafm =',isym,trialnons(:),trialafm
!  ENDDEBUG

!  Loop over all classes, then all atoms in the class,
!  to find whether they have a symmetric
   do iclass=1,nclass
    do jj=1,natomcl(iclass)

     iatom2=class(jj,iclass)
!    Generate the tentative symmetric position of iatom2
     symxred2(:)=ptsymrel(:,1,isym)*xred(1,iatom2)+ &
&     ptsymrel(:,2,isym)*xred(2,iatom2)+ &
&     ptsymrel(:,3,isym)*xred(3,iatom2)+ trialnons(:)
!    Generate the tentative symmetric spinat of iatom2
     symspinat2(:)=trialafm*spinat(:,iatom2)

!    Check whether there exists an atom of the same class at the
!    same location, with the correct spinat
     do kk=1,natomcl(iclass)

      found3=1
      iatom3=class(kk,iclass)
!     Check the location
      diff(:)=xred(:,iatom3)-symxred2(:)
      diff(:)=diff(:)-nint(diff(:))
      if( (diff(1)**2+diff(2)**2+diff(3)**2) > tol8**2 )found3=0
!     Check the spinat
      diff(:)=spinat(:,iatom3)-symspinat2(:)
      if( (diff(1)**2+diff(2)**2+diff(3)**2) > tol8**2 )found3=0
      if(found3==1)exit

!     End loop over iatom3
     end do

     if(found3==0)then
      trialok=0
      exit
     end if

!    End loop over iatom2
    end do

    if(trialok==0)exit

!   End loop over all classes
   end do

   if(trialok==1)then
    nsym=nsym+1
    if(nsym>msym)then
     write(message,'(a,a,a,a,a,a,i4,a,a,a,a)')ch10,&
&     ' symfind : ERROR -',ch10,&
&     '  The number of symmetries (including non-symmorphic',ch10,&
&     '  translations) is larger than msym=',msym,ch10,&
&     '  Action : take a cell that is primitive, or at least',ch10,&
&     '  smaller than the present one.'
     call wrtout(6,message,'COLL')
     call leave_new('COLL')
    end if
    ntrial=ntrial+1
    symrel(:,:,nsym)=ptsymrel(:,:,isym)
    symafm(nsym)=trialafm
    tnons(:,nsym)=trialnons(:)-nint(trialnons(:)-tol8)
   end if

!  End the loop on tentative translations
  end do

! End big loop over each symmetry operation of the Bravais lattice
 end do

 deallocate(class,natomcl,spinatcl,typecl)

!DEBUG
!write(6,*)' symfind : exit, nsym=',nsym
!write(6,*)'   symrel matrices and symafm are :'
!do isym=1,nsym
!write(6, '(i4,4x,10i4)' )isym,symrel(:,:,isym),symafm(isym)
!end do
!ENDDEBUG

end subroutine symfind
!!***
