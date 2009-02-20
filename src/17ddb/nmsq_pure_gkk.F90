!{\src2tex{textfont=tt}}
!!****f* ABINIT/nmsq_pure_gkk
!!
!! NAME
!! nmsq_pure_gkk
!!
!! FUNCTION
!!  Calculate gamma matrices for pure gkk case, ie when the
!!  scalar product with the displacement vector is done later
!!  Sum over bands is carried out later.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   displ_red = phonon displacement in reduced coordinates (used to calculate the ph linewidth)
!!   elph_ds = datastructure with gkk matrix elements
!!   FSfullpqtofull = mapping of k+q to k
!!   FSintweight = FS integration weights for each band and kpt
!!   FSkpt = coordinates of kpoints near to FS
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ squared
!!   iqptfull = index of present qpoint
!!   phfrq_tmp = phonon frequencies
!!   spqpt = array of qpoint coordinates
!!   wf = gkk matrix element weight with $1/\sqrt{2 M \omega}$
!!
!! OUTPUT
!!   elph_ds%gkq filled
!!   accum_mat = matrix for accumulating FS average of gkk (gamma matrix -> linewidths)
!!   accum_mat2 = complex array whose real part contains the phonon linewidth
!!   gkk_qpt_tmp = tmp matrix for all gamma matrix elements, saved to disk or to memory in nmsq_gam_sumFS
!!
!! NOTES
!!
!! PARENTS
!!      normsq_gkq
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nmsq_pure_gkk(accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,FSintweight,FSkpt,gkk_qpt_tmp,&
&   h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptfull
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
 real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt,elph_ds%nsppol)
 real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
 real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch),spqpt(3,elph_ds%nqpt)
 real(dp),intent(in) :: wf(elph_ds%nbranch)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: &
& gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)

!Local variables-------------------------------
!scalars
 integer :: goodkpq,iFSkpt,iFSkptq,ib1,ib2,ibeff,ibranch,ipert1,isppol,jbranch
 integer :: kbranch
 real(dp) :: res,sd1,sd2,ss
 character(len=500) :: message
!arrays
 real(dp) :: gkq(3),gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch),kpt(3)
 real(dp) :: redkpt(3),tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)

! *************************************************************************

#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'ZGEMM' :: zgemm
#endif

 if (elph_ds%tkeepbands /= 1) then
  write (message,'(3a)')' nmsq_pure_gkk : BUG- ',ch10,&
&  ' elph_ds%tkeepbands should be 1 to keep bands!'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!
!just copy  h1_mat_el_sq into gkk_qpt_tmp
!
 gkk_qpt_tmp = pi*h1_mat_el_sq


!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

!MJV 20070525 NOTE:
!in some of the nmsq routines, in particular this one, the work done to
!calculate accum_mat,accum_mat2 is completely superfluous and will be re-done
!on the interpolated values.
!MG uses them for the QPT output, however, so keep it for consistency for the
!moment.

 do isppol=1,elph_ds%nsppol
  do iFSkpt=1,elph_ds%nFSkpt

   iFSkptq = FSfullpqtofull(iFSkpt,iqptfull)

   gkq_sum_bands(:,:,:) = zero

!  gkq_sum_bands = \sum_{ib1,ib2} \langle k+q \mid H^{(1)}_{q,\tau_i,\alpha_i} \mid k   \rangle
!  \cdot \langle k   \mid H^{(1)}_{q,\tau_j,\alpha_j} \mid k+q \rangle
!  where ibranch -> \tau_i,\alpha_i  and  jbranch -> \tau_j,\alpha_j

   do ib1=1,elph_ds%nFSband

    sd1 = FSintweight(ib1,iFSkpt,isppol)      !  weights for distance from the fermi surface

    do ib2=1,elph_ds%nFSband

     sd2 = FSintweight(ib2,iFSkptq,isppol)  !  weights for distance from the fermi surface
     ibeff = ib2+(ib1-1)*elph_ds%nFSband

     gkq_sum_bands = gkq_sum_bands + &
&     sd1*sd2*pi*reshape(h1_mat_el_sq(:,ibeff,:,iFSkpt,isppol),(/2,elph_ds%nbranch,elph_ds%nbranch/))

    end do !ib2
   end do !ib1
!  END loops over bands


   accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
  end do
! END loop over FSkpt

! MG20060603
! do scalar product with the displ_red to calculate the ph lwdth before interpolation (stored in accum_mat2)

  zgemm_tmp_mat=zero
  tmp_mat2 = accum_mat(:,:,:,isppol)
  call zgemm('c','n',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&  displ_red,elph_ds%nbranch,tmp_mat2,&
&  elph_ds%nbranch,czero,zgemm_tmp_mat,elph_ds%nbranch)

! MG20060607 there is no explicit dependence on omega in the linewidth.
! It is better dont use wf at all and employ the same approach as in nmsq_gam or nmsq_pure_gkk
  tmp_mat2=zero
  call zgemm('n','n',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&  zgemm_tmp_mat,elph_ds%nbranch,displ_red,&
&  elph_ds%nbranch,czero,tmp_mat2,elph_ds%nbranch)

  do ipert1=1,elph_ds%nbranch
   accum_mat2(1,ipert1,ipert1,isppol) = accum_mat2(1,ipert1,ipert1,isppol) + tmp_mat2(1,ipert1,ipert1)
  end do

! ENDMG

 end do ! isppol

!DEBUG
!write(73,'(a,3es16.8)')'#nmsq_pure_gkk QPT ',spqpt(:,iqptfull)
!write(73,'(3es16.8,3(2e16.8))')(accum_mat2(1,jbranch,jbranch), jbranch=1,elph_ds%nbranch)
!ENDDEBUG

end subroutine nmsq_pure_gkk
!!***
