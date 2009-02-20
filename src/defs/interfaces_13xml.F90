!!****m* ABINIT/interfaces_13xml
!! NAME
!! interfaces_13xml
!!
!! FUNCTION
!! This module contains the interfaces of the routines
!! in the directory src/13xml
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! THIS FILE IS GENERATED AUTOMATICALLY BY abilint.
!! To do that: config/scripts/abilint . .
!! 
!!
!! SOURCE

module interfaces_13xml

 implicit none

interface
 subroutine append_cml (dtset_char,lenstr,lenstr_cml,string,string_cml,strln)
  implicit none
  integer,intent(inout) :: lenstr
  integer,intent(in) :: lenstr_cml
  integer,intent(in) :: strln
  character(len=2),intent(in) :: dtset_char
  character(len=strln),intent(inout) :: string
  character(len=strln),intent(in) :: string_cml
 end subroutine append_cml
end interface

interface
 subroutine append_cml2(dtset_char,lenstr,lenstr_cml,string,string_cml,strln)
  implicit none
  integer,intent(inout) :: lenstr
  integer,intent(in) :: lenstr_cml
  integer,intent(in) :: strln
  character(len=2),intent(in) :: dtset_char
  character(len=strln),intent(inout) :: string
  character(len=strln),intent(in) :: string_cml
 end subroutine append_cml2
end interface

interface
 subroutine findmarkup(builtin,index_lower,index_upper,indices_markup,&  
  &  markup,markuplen,strln,string_xml)
  implicit none
  integer,intent(in) :: index_lower
  integer,intent(in) :: index_upper
  integer,intent(in) :: markuplen
  integer,intent(in) :: strln
  character(len=*),intent(in) :: builtin
  character(len=markuplen),intent(in) :: markup
  character(len=strln),intent(in) :: string_xml
  integer,intent(out) :: indices_markup(3)
 end subroutine findmarkup
end interface

interface
 subroutine getattribute(attribute,attributelen,indices_markup,strln,string_xml,valattrib)
  implicit none
  integer,intent(in) :: attributelen
  integer,intent(in) :: strln
  character(len=attributelen),intent(in) :: attribute
  character(len=strln),intent(in) :: string_xml
  character(len=*),intent(out) :: valattrib
  integer,intent(in) :: indices_markup(3)
 end subroutine getattribute
end interface

interface
 subroutine importcml (lenstr,string_raw,string_upper,strln)
  implicit none
  integer,intent(inout) :: lenstr
  integer,intent(in) :: strln
  character(len=*),intent(in) :: string_raw
  character(len=*),intent(inout) :: string_upper
 end subroutine importcml
end interface

interface
 subroutine prt_cml(filapp,natom,nsym,ntypat,rprimd,spgroup,symrel,tnons,typat,xred,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: spgroup
  character(len=fnlen),intent(in) :: filapp
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine prt_cml
end interface

interface
 subroutine prt_cml2(filapp,natom,nsym,ntypat,rprimd,spgroup,symrel,tnons,typat,xred,znucl)
  use defs_basis
  implicit none
  integer,intent(in) :: natom
  integer,intent(in) :: nsym
  integer,intent(in) :: ntypat
  integer,intent(in) :: spgroup
  character(len=fnlen),intent(in) :: filapp
  real(dp),intent(in) :: rprimd(3,3)
  integer,intent(in) :: symrel(3,3,nsym)
  real(dp),intent(in) :: tnons(3,nsym)
  integer,intent(in) :: typat(natom)
  real(dp),intent(in) :: xred(3,natom)
  real(dp),intent(in) :: znucl(ntypat)
 end subroutine prt_cml2
end interface

end module interfaces_13xml
!!***
