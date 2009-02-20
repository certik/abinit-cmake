!{\src2tex{textfont=tt}}
!!****f* ABINIT/out_resultsgs_xml
!!
!! NAME
!! out_resultsgs_xml
!!
!! FUNCTION
!! Output in the XML file, the decomposition of the energy and
!! the forces after a scfcv loop.
!! (see extras/post_processing/abinitRun.dtd)
!! Warning : this method is not thread safe and should be called
!! only by one thread.
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
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation.
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine out_resultsgs_XML(dtset, level, results_gs, usepaw)

 use defs_basis
  use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: level,usepaw
 type(dataset_type),intent(in) :: dtset
 type(results_gs_type),intent(inout) :: results_gs

!Local variables -------------------------
  character(len = 1), parameter :: axes_names(3) = (/ "x", "y", "z" /)
!scalars
 integer :: i,j
 character(len=128) :: spacer,value

! *************************************************************************

!Compute the spacer to put before each markup
 write(spacer, "(A,I0)") "A", 2 * level

!Begin with the energy part
 write(ab_xml_out, "("//trim(spacer)//",A)", advance = "NO") " ", '<energy'
 if (dtset%iscf < 10) then
  write(ab_xml_out, "(A)", advance = "NO") ' type="direct"'
  write(value, *) results_gs%energies%e_kinetic
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' kinetic="', trim(value) ,'"'
  write(value, *) results_gs%energies%e_localpsp
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' local="', trim(value) ,'"'
  write(value, *) results_gs%energies%e_nonlocalpsp
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' non-local="', trim(value) ,'"'
  if (usepaw == 1) then
   write(value, *) results_gs%energies%e_paw
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' paw="', trim(value) ,'"'
  end if
 else
  write(ab_xml_out, "(A)", advance = "NO") ' type="double-counting"'
  write(value, *) results_gs%energies%e_eigenvalues
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' eigen-values="', trim(value) ,'"'
  write(value, *) results_gs%energies%e_xcdc
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' xcdc="', trim(value) ,'"'
  if (usepaw == 1) then
   write(value, *) results_gs%energies%e_pawdc
   write(ab_xml_out, "(A,A,A)", advance = "NO") ' pawdc="', trim(value) ,'"'
  end if
 end if
 if (dtset%berryopt == 4) then
  write(value, *) results_gs%energies%e_elecfield
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' electric-field="', trim(value) ,'"'
 end if
 if(dtset%occopt >= 3 .and. dtset%occopt <= 7) then
  write(value, *) results_gs%energies%e_entropy
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' entropy="', trim(value) ,'"'
 end if
 write(value, *) results_gs%energies%e_ewald
 if (dtset%icoulomb == 0) then
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' ewald="', trim(value) ,'"'
 else
  write(ab_xml_out, "(A,A,A)", advance = "NO") ' ion-ion="', trim(value) ,'"'
 end if
 write(value, *) results_gs%energies%e_hartree
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' hartree="', trim(value) ,'"'
 write(value, *) results_gs%energies%e_corepsp
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' core="', trim(value) ,'"'
 write(value, *) results_gs%energies%e_xc
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' xc="', trim(value) ,'"'
 write(value, *) results_gs%etotal
 write(ab_xml_out, "(A,A,A)", advance = "NO") ' total="', trim(value) ,'"'
 write(ab_xml_out, "(A)") ' />'


!finish with the forces part
 if (dtset%optforces == 1) then
  write(ab_xml_out, "("//trim(spacer)//",A)") " ", '<forces>'
  do i = 1, dtset%natom, 1
   write(ab_xml_out, "("//trim(spacer)//",A)", advance = "NO") " ", '  <force'
   write(ab_xml_out, "(A,I0,A,I0,A)", advance = "NO") ' atom="a_', dtset%jdtset ,'_', i ,'"'
   do j = 1, 3, 1
    write(value, *) results_gs%fcart(j, i)
    write(ab_xml_out, "(A,A,A,A,A)", advance = "NO") ' ', axes_names(j), '="', trim(value) ,'"'
   end do
   write(ab_xml_out, "(A)") ' />'
  end do
  write(ab_xml_out, "("//trim(spacer)//"A)") " ", '</forces>'
 end if

end subroutine out_resultsgs_XML
!!***
