!{\src2tex{textfont=tt}}
!!****f* ABINIT/WffKg
!! NAME
!! WffKg
!!
!! FUNCTION
!! Check kgwff to  manage WF file in the MPI/IO case
!!
!! COPYRIGHT
!! Copyright (C) 2003-2008 ABINIT group (MB)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  wff <type(wffile_type)> = structured info about the wavefunction file
!!  optkg= if 1 , read or write kg_k ; if 0,do not care about kg_k in rwwf
!!
!! OUTPUT
!!
!! PARENTS
!!      inwffil,nstdy3,outwf,vtorho,vtorho3
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine WffKg(wff,optkg)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 type(wffile_type),intent(inout) :: wff
 integer          ,intent(in)    :: optkg

!Local variables ------------------------------
#if defined MPI_IO
           if ( wff%accesswff == 1) then
            wff%kgwff=optkg
           end if
#endif

end subroutine WffKg
!!***
