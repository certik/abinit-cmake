!{\src2tex{textfont=tt}}
!!****f* ABINIT/out_geometry_xml
!!
!! NAME
!! out_geometry_xml
!!
!! FUNCTION
!! Output in the XML file, the box size and the atomic position.
!! (see extras/post_processing/abinitRun.dtd)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  level=use for indentation of the XML, two spaces are added by level.
!!  natom=number of atoms.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine out_geometry_XML(dtset, level, natom, rprimd, xred)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: level,natom
 type(dataset_type),intent(in) :: dtset
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: xred(3,natom)

!Local variables -------------------------
  character(len = 1), parameter :: rprimd_names(3) = (/ "x", "y", "z" /)
!scalars
 integer :: i,j
 character(len=128) :: spacer,value
!arrays
 real(dp),allocatable :: xcart(:,:)

! *************************************************************************

!Compute the spacer to put before each markup
 write(spacer, "(A,I0)") "A", 2 * level
 write(ab_xml_out, "("//trim(spacer)//",A)") " ", '<geometry>'
!Compute the cartesian coordinates of atoms
 allocate(xcart(3, natom))
 call xredxcart(natom, 1, rprimd, xcart, xred)
!Ouput the rprimd matrix
 write(ab_xml_out, "("//trim(spacer)//",A)", advance = "NO") " ", '  <rprimd'
 do i = 1, 3, 1
  do j = 1, 3, 1
   write(value, *) rprimd(i, j)
   write(ab_xml_out, "(A,A,I0,A,A,A)", advance = "NO") ' ', rprimd_names(i), j, '="', trim(value) ,'"'
  end do
 end do
 write(ab_xml_out, "(A)") ' />'
!Output the atom position
 do i = 1, natom, 1
  write(ab_xml_out, "("//trim(spacer)//",A)", advance = "NO") " ", '  <position'
  write(ab_xml_out, "(A,I0,A,I0,A)", advance = "NO") ' atom="a_', dtset%jdtset ,'_', i ,'"'
  do j = 1, 3, 1
   write(value, *) xcart(j, i)
   write(ab_xml_out, "(A,A,A,A,A)", advance = "NO") ' ', rprimd_names(j), '="', trim(value) ,'"'
  end do
  write(ab_xml_out, "(A)") ' />'
 end do
 deallocate(xcart)
 write(ab_xml_out, "("//trim(spacer)//",A)") " ", '</geometry>'
end subroutine out_geometry_XML
!!***
