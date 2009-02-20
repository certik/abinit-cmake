!{\src2tex{textfont=tt}}
!!****f* ABINIT/prt_cml
!! NAME  prt_cml
!! prt_cml
!!
!!
!! FUNCTION
!! Produce a CML (Chemical Markup Language) file
!! with crystalline cell description, symmetries,
!! and reduced atomic coordinates.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filapp= character string giving the root to form the name of the CML file
!!  natom=number of atoms in unit cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  rprimd(3,3)=real space dimensional primitive translations (bohr)
!!  spgroup=symmetry space group number
!!  symrel(3,3,nsym)=symmetry matrices in real space (integers)
!!  tnons(3,nsym)=reduced nonsymmorphic translations
!!  typat(natom)=type integer for each atom in cell
!!  xred(3,natom)=reduced coordinates of atoms
!!  znucl(ntypat)=real(dp), atomic number of atom type
!!
!! OUTPUT
!! data written in file whose name is filapp//'_CML.xml'
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      atmdata,inupper,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prt_cml(filapp,natom,nsym,ntypat,rprimd,spgroup,symrel,tnons,typat,xred,znucl)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym,ntypat,spgroup
 character(len=fnlen),intent(in) :: filapp
!arrays
 integer,intent(in) :: symrel(3,3,nsym),typat(natom)
 real(dp),intent(in) :: rprimd(3,3),tnons(3,nsym),xred(3,natom),znucl(ntypat)

!Local variables-------------------------------
!scalars
 integer :: ia,iablock,ii,isym,mu,nu
 real(dp) :: amu,rcov
 character(len=2) :: string2,symbol
 character(len=20) :: string20
 character(len=3) :: string3
 character(len=500) :: message
 character(len=80) :: string80
 character(len=fnlen) :: filxml
!arrays
 real(dp) :: angle(3),rmet(3,3)
 character(len=3) :: string2array(3,3)

! *************************************************************************

!Initialize the file
 filxml=trim(filapp)//'_CML.xml'
 write(message, '(a,a)' ) ' prt_cml : about to open file ',filxml
 call wrtout(6,message,'COLL')
 call wrtout(ab_out,message,'COLL')
 open (unit=tmp_unit,file=filxml,status='unknown',form='formatted')
 rewind(tmp_unit)

!Take care of the header, and initialize the <molecule> element
 write(message, '(a,a,a,a,a)' ) &
& '<?xml version="1.0" encoding="iso-8859-1"?>',ch10,&
& '<!DOCTYPE molecule SYSTEM "CML-1999-05-15.dtd" []>',ch10,&
& '<molecule class="3d">'
 call wrtout(tmp_unit,message,'COLL')

!Compute real space metrics
 do ii=1,3
  rmet(ii,:)=rprimd(1,ii)*rprimd(1,:)+&
&  rprimd(2,ii)*rprimd(2,:)+&
&  rprimd(3,ii)*rprimd(3,:)
 end do

!Compute angles in degree
 angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0d0
 angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0d0
 angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0d0

!Write the <crystal> element
 write(string3, '(i3)')spgroup
 write(message, '(13a)' )&
& ' <crystal>',ch10,&
& '  <!--',ch10,&
& '   "ITC" refers to the space group number ',ch10,&
& '   in the International Tables for Crystallography',ch10,&
& '   -->',ch10,&
& '  <string builtin="spacegroup" convention="ITC">',&
& trim(adjustl(string3)),'</string>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')sqrt(rmet(1,1))*Bohr_Ang
 write(message, '(a,a,a)' )&
& '  <float builtin="acell" units="angstrom">',&
& trim(adjustl(string20)),'</float>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')sqrt(rmet(2,2))*Bohr_Ang
 write(message, '(a,a,a)' )&
& '  <float builtin="bcell" units="angstrom">',&
& trim(adjustl(string20)),'</float>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')sqrt(rmet(3,3))*Bohr_Ang
 write(message, '(a,a,a)' )&
& '  <float builtin="ccell" units="angstrom">',&
& trim(adjustl(string20)),'</float>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')angle(1)
 write(message, '(a,a,a)' )&
& '  <float builtin="alpha" units="degrees">',&
& trim(adjustl(string20)),'</float>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')angle(2)
 write(message, '(a,a,a)' )&
& '  <float builtin="beta"  units="degrees">',&
& trim(adjustl(string20)),'</float>'
 call wrtout(tmp_unit,message,'COLL')
 write(string20, '(f20.12)')angle(3)
 write(message, '(a,a,a)' )&
& '  <float builtin="gamma" units="degrees">',&
& trim(adjustl(string20)),'</float>'
 call wrtout(tmp_unit,message,'COLL')
 write(message, '(a)' )' </crystal>'
 call wrtout(tmp_unit,message,'COLL')

!Write the <symmetry> element
 write(message, '(a)' )' <symmetry>'
 call wrtout(tmp_unit,message,'COLL')
 do isym=1,nsym
  do mu=1,3
   do nu=1,3
    write(string2,'(i2)')symrel(mu,nu,isym)
    string2array(mu,nu)=string2
   end do
  end do
  write(message, '(a)')&
&  '  <floatMatrix builtin="symOp" rows="4" columns="3">'
  call wrtout(tmp_unit,message,'COLL')
  do mu=1,3
   write(string20,'(f20.12)')tnons(mu,isym)
   write(message, '(8a)' ) '    ',&
!  XG 021103 : I am not sure that the order of the symrel indices is the right one
&   trim(adjustl(string2array(mu,1))),' ',&
&   trim(adjustl(string2array(mu,2))),' ',&
&   trim(adjustl(string2array(mu,3))),' ',&
&   trim(adjustl(string20))
   call wrtout(tmp_unit,message,'COLL')
  end do
  write(message, '(a)')'  </floatMatrix>'
  call wrtout(tmp_unit,message,'COLL')
 end do
 write(message, '(a)' )' </symmetry>'
 call wrtout(tmp_unit,message,'COLL')

!Initialize the <atomArray> element
 write(message, '(a)' )' <atomArray>'
 call wrtout(tmp_unit,message,'COLL')

!Write the elementType
 write(message, '(a)' )'  <stringArray builtin="elementType">'
 call wrtout(tmp_unit,message,'COLL')
!Concatenate up to 26 symbols
 do iablock=1,natom,26
  string80=repeat(' ',80)
  do ia=iablock,min(iablock+25,natom)
!  Get the symbol of the element
   call atmdata(amu,rcov,symbol,znucl(typat(ia)))
!  Switch to uppercase
   call inupper(symbol)
   string80=trim(string80)//' '//trim(adjustl(symbol))
  end do
  write(message, '(a,a)' ) '   ',trim(string80)
  call wrtout(tmp_unit,message,'COLL')
 end do
 write(message, '(a)' )'  </stringArray>'
 call wrtout(tmp_unit,message,'COLL')

!Write xFract, yFract and zFract
 do ii=1,3
  if(ii==1)write(message, '(a)' )'  <floatArray builtin="xFract">'
  if(ii==2)write(message, '(a)' )'  <floatArray builtin="yFract">'
  if(ii==3)write(message, '(a)' )'  <floatArray builtin="zFract">'
  call wrtout(tmp_unit,message,'COLL')
! Concatenate up to 5 reduced coordinates
  do iablock=1,natom,5
   string80=repeat(' ',80)
   do ia=iablock,min(iablock+4,natom)
!   Get the reduced coordinate of the element
    write(string20,'(f20.12)')xred(ii,ia)
!   Switch to uppercase
    string80=trim(string80)//' '//trim(adjustl(string20))
   end do
   write(message, '(a,a)' ) '   ',trim(string80)
   call wrtout(tmp_unit,message,'COLL')
  end do
  write(message, '(a)' )'  </floatArray>'
  call wrtout(tmp_unit,message,'COLL')
 end do

!Finalize the <atomArray> element
 write(message, '(a)' )' </atomArray>'
 call wrtout(tmp_unit,message,'COLL')

!Finalize the <molecule> element
 write(message, '(a)' ) '</molecule>'
 call wrtout(tmp_unit,message,'COLL')

 close(tmp_unit)

end subroutine prt_cml
!!***
