!{\src2tex{textfont=tt}}
!!****f* ABINIT/fxgkkphase
!!
!! NAME
!! fxgkkphase
!!
!! FUNCTION
!!   Set phase factors to eliminate gauge variable in gkk matrix elements
!!    (comes from phase factors in wavefunctions)
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
!!   gkk_flag = flags for presence of gkk matrix elements
!!   h1_mat_el = irreducible matrix elements to be completed
!!   iqptfull = qpoint number in full zone
!!
!! OUTPUT
!!   h1_mat_el = changed on output
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fxgkkphase(elph_ds,gkk_flag,h1_mat_el,iqptfull)

 use defs_basis
 use defs_datatypes
 use defs_elphon

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptfull
 type(elph_type),intent(in) :: elph_ds
!arrays
 integer,intent(in) :: gkk_flag(elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nqpt)
 real(dp),intent(inout) :: h1_mat_el(2,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nbranch,elph_ds%nFSkpt)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,ib1,ib2,ibranch,ipreatom,ipredir,isym,isymatom,isymdir,itim
 integer :: nsymperts,prepert,symiFSkpt,sympertcase
 real(dp) :: exparg,imphase,norm,rephase,s1,s2,s3,ss,sumi,sumr,timsign,tmpim
 real(dp) :: tmpre
!arrays
 real(dp) :: dsymrec(3,3)

! *************************************************************************

 do iFSkpt=1,elph_ds%nFSkpt
  do ib1=1,elph_ds%nFSband
   do ib2=1,elph_ds%nFSband
!   Determine phase for first perturbation
    rephase =  h1_mat_el(1,ib2,ib1,1,iFSkpt)
    imphase = -h1_mat_el(2,ib2,ib1,1,iFSkpt)
    norm = sqrt(rephase**2+imphase**2)
    rephase = rephase / norm
    imphase = imphase / norm

!   DEBUG
!   if (iFSkpt == 1) then
!   write (*,*) 'fxgkkphase : rephase,imphase = ', rephase,imphase
!   end if
!   ENDDEBUG

!   apply same phase factor to all perturbations
!   ----------------------------------------------------------
!   Very important ! Otherwise the scalar product with the
!   displacement vector will not be preserved.
!   ----------------------------------------------------------
    do ibranch=1,elph_ds%nbranch
!    if we already have data
     if (gkk_flag(ibranch,1,iqptfull) /= -1) then
      tmpre =    rephase * h1_mat_el(1,ib2,ib1,ibranch,iFSkpt)&
&      -imphase * h1_mat_el(2,ib2,ib1,ibranch,iFSkpt)
      tmpim =    rephase * h1_mat_el(2,ib2,ib1,ibranch,iFSkpt)&
&      +imphase * h1_mat_el(1,ib2,ib1,ibranch,iFSkpt)
      h1_mat_el(1,ib2,ib1,ibranch,iFSkpt) = tmpre
      h1_mat_el(2,ib2,ib1,ibranch,iFSkpt) = tmpim
     end if
!    end if
    end do
   end do
  end do
 end do
!end loop over FS kpt

end subroutine fxgkkphase
!!***
