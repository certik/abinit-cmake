!{\src2tex{textfont=tt}}
!!****f* ABINIT/nmsq_gam
!!
!! NAME
!! nmsq_gam
!!
!! FUNCTION
!!  Calculate gamma matrices keeping full dependence on bands
!!  from original h1_mat_el_sq matrix elements (no averaging over
!!  bands near the Fermi surface)
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   displ_red = phonon mode displacement vectors, post-multiplied by gprim matrix
!!     (ie. turned to reduced coordinates)
!!   eigvec = phonon eigenvectors
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
!!   accum_mat = matrix for accumulating FS average of gkk (gamma matrix -> linewidths)
!!   accum_mat2 = matrix for accumulating FS average of gamma matrix with good prefactors
!!   gkk_qpt_tmp = tmp matrix for all gamma matrix elements, saved to disk or to memory
!!
!! NOTES
!!
!! PARENTS
!!      normsq_gkq
!!
!! CHILDREN
!!      leave_new,wrtout,zgemm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&   FSintweight,FSkpt,gkk_qpt_tmp,h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)

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
 real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: &
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
 real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch),spqpt(3,elph_ds%nqpt)
 real(dp),intent(in) :: wf(elph_ds%nbranch)
 real(dp),intent(inout) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp),intent(inout) :: &
& gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)

!Local variables-------------------------------
! tmp variables for diagonalization
!scalars
 integer :: goodkpq,iFSkpt,iFSkptq,ib1,ib2,ibeff,ibranch,ipert1,ipert2,isppol
 integer :: jbranch,kbranch
 real(dp) :: res,sd1,sd2,ss
 character(len=500) :: message
!arrays
 real(dp) :: gkq(3),gkq_1band(2,elph_ds%nbranch,elph_ds%nbranch),kpt(3)
 real(dp) :: redkpt(3),tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),allocatable :: matrx(:,:),val(:),vec(:,:,:),zhpev1(:,:),zhpev2(:)

! *************************************************************************

#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'ZGEMM' :: zgemm
#endif

 if (elph_ds%tkeepbands == 0) then
  write (message,'(3a,i3)')' nmsq_gam : BUG- ',ch10,&
&  ' elph_ds%tkeepbands should be 1 while is ',elph_ds%tkeepbands
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

 write (*,*) 'nmsq_gam : iqptfull = ', iqptfull

 do isppol=1,elph_ds%nsppol
  do iFSkpt=1,elph_ds%nFSkpt

   iFSkptq = FSfullpqtofull(iFSkpt,iqptfull)

   do ib1=1,elph_ds%nFSband
    sd1 = FSintweight(ib1,iFSkpt,isppol) !weights for distance from the fermi surface

    do ib2=1,elph_ds%nFSband
     sd2 = FSintweight(ib2,iFSkptq,isppol) !weights for distance from the fermi surface
     ibeff = ib2+elph_ds%nFSband*(ib1-1)

     gkq_1band(:,:,:) = zero

!    OLDVERSION
!    do jbranch=1,elph_ds%nbranch !phonon mode index for displ_red
!    
!    gkq(:) = zero
!    ! Do scalar product of displacement vector by reduced perturbations
!    ! loops over atoms and red directions
!    do ipert1=1,elph_ds%nbranch
!    do ipert2=1,elph_ds%nbranch
!    ! should this be displ_red complex conjugated? Yes from DECAFT
!    !gkq(1) = gkq(1)&
!    !     & + displ_red(1,ibranch,jbranch)*h1_mat_el_sq(1,ib2,ib1,ibranch,iFSkpt) &
!    !     & + displ_red(2,ibranch,jbranch)*h1_mat_el_sq(2,ib2,ib1,ibranch,iFSkpt)
!    !gkq(2) = gkq(2)&
!    !     & + displ_red(1,ibranch,jbranch)*h1_mat_el_sq(2,ib2,ib1,ibranch,iFSkpt) &
!    !     & - displ_red(2,ibranch,jbranch)*h1_mat_el_sq(1,ib2,ib1,ibranch,iFSkpt)
!    gkq(1) = gkq(1)&
!    &                 + displ_red(1,ipert2,jbranch)*h1_mat_el_sq(1,ib2,ib1,ipert2,ipert1,iFSkpt)*displ_red(1,ipert1,jbranch)&
!    &                 + displ_red(2,ipert2,jbranch)*h1_mat_el_sq(2,ib2,ib1,ipert2,ipert1,iFSkpt)*displ_red(1,ipert1,jbranch)&
!    &                 - displ_red(1,ipert2,jbranch)*h1_mat_el_sq(2,ib2,ib1,ipert2,ipert1,iFSkpt)*displ_red(2,ipert1,jbranch)&
!    &                 + displ_red(2,ipert2,jbranch)*h1_mat_el_sq(1,ib2,ib1,ipert2,ipert1,iFSkpt)*displ_red(2,ipert1,jbranch)
!    
!    !MG20060603   Fixed a minor bug present in version 5.1.3
!    gkq(2) = gkq(2)&
!    &                 + displ_red(1,ipert2,jbranch)*h1_mat_el_sq(2,ib2,ib1,ipert2,ipert1,iFSkpt)*displ_red(1,ipert1,jbranch)&
!    &                 + displ_red(1,ipert2,jbranch)*h1_mat_el_sq(1,ib2,ib1,ipert2,ipert1,iFSkpt)*displ_red(2,ipert1,jbranch)&
!    &                 - displ_red(2,ipert2,jbranch)*h1_mat_el_sq(1,ib2,ib1,ipert2,ipert1,iFSkpt)*displ_red(1,ipert1,jbranch)&
!    &                 + displ_red(2,ipert2,jbranch)*h1_mat_el_sq(2,ib2,ib1,ipert2,ipert1,iFSkpt)*displ_red(2,ipert1,jbranch)
!    
!    end do
!    end do
!    
!    if (abs(gkq(2)) > tol8) then
!    write (message,'(3a,es16.8)')' nmsq_gam: WARNING-  gkq is not real!',ch10,&
!    &               ' Im(gkq) = ',gkq(2)
!    call wrtout(06,message,'COLL')
!    end if
!    
!    !  Add weights for integration over FS and wf = (1/sqrt(2 omega))
!    !     If we are not going to average immediately over FS, do not
!    !     include FSintweight now
!    !gkq_1band(1,jbranch,jbranch) = &
!    !     &   gkq_1band(1,jbranch,jbranch) + &
!    !     &   wf(jbranch)**2 * (gkq(1)**2+gkq(2)**2)
!    
!    !VERSION5.1.3
!    !gkq_1band(1,jbranch,jbranch) = gkq_1band(1,jbranch,jbranch) + wf(jbranch)**2 * gkq(1)
!    !ENDVERSION5.1.3
!    
!    !MG20060607 there is no explicit dependence on omega in the linewidth. It is better dont use wf at all
!    !           now summing gkk_qpt over kpoints gives the phonon linewidth
!    gkq_1band(1,jbranch,jbranch) = gkq_1band(1,jbranch,jbranch) + pi*gkq(1)
!    !ENDMG20060607
!    
!    end do
!    ENDOLDVERSION


     zgemm_tmp_mat=zero
     tmp_mat2 = reshape (h1_mat_el_sq(:,ibeff,:,iFSkpt,isppol),(/2,elph_ds%nbranch,elph_ds%nbranch/))

     call zgemm('c','n',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&     displ_red,elph_ds%nbranch,tmp_mat2,&
&     elph_ds%nbranch,czero,zgemm_tmp_mat,elph_ds%nbranch)

!    MG20060607 there is no explicit dependence on omega in the linewidth.
!    It is better dont use wf at all and employ the same approach as in nmsq_gam or nmsq_pure_gkk
     tmp_mat2=zero
!    factor of pi here
     call zgemm('n','n',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,dcmplx(pi,zero),&
&     zgemm_tmp_mat,elph_ds%nbranch,displ_red,&
&     elph_ds%nbranch,czero,tmp_mat2,elph_ds%nbranch)
!    sum over bands 
     do ipert1=1,elph_ds%nbranch
      gkq_1band(1,ipert1,ipert1) = gkq_1band(1,ipert1,ipert1) + tmp_mat2(1,ipert1,ipert1)
     end do

!    DEBUG
!    write (*,'(a,a,4i8,2E16.6)') ' normsq_gkq : ',&
!    &  'ib1,ib2,iFSkpt,iqptfull,gkq_1band(1,1,1) = ',&
!    &   ib1,ib2,iFSkpt,iqptfull,gkq_1band(:,1,1)
!    write (*,*) '###', elph_ds%gkk_qpt(:,hdr1%pertcase,ib1,ib2,iFSkpt)
!    ENDDEBUG
!    summing over k points and bands, still diagonal in jbranch
     accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_1band(:,:,:)*sd1*sd2

!    MG20060603 : summing over bands and kpoints with weights to calculate the phonon linewidth
     do jbranch=1,elph_ds%nbranch
      accum_mat2(:,jbranch,jbranch,isppol) = accum_mat2(:,jbranch,jbranch,isppol) + gkq_1band(:,jbranch,jbranch)*sd1*sd2
!     &            +(two*pi*abs(phfrq_tmp(jbranch)))* gkq_1band(:,jbranch,jbranch)*sd1*sd2
     end do
!    END MG


!    ! 28 Apr 2004  PI factor added to wf in read_gkk
!    ! 15 jun 2004  put back into wf
!    ! Add factor of 2 pi omega(branch) to gamma. Still in diagonal normal mode basis
!    do jbranch=1,elph_ds%nbranch
!    gkq_1band(:,jbranch,jbranch) = two * pi * phfrq_tmp(jbranch) * gkq_1band(:,jbranch,jbranch)
!    end do

!    now turn to cartesian coordinates

!    Final Gamma matrix (hermitian) = E * D_g * E^{+}
!    Where E^{+} is the hermitian conjugate of the eigenvector matrix E
!    And D_g is the diagonal matrix of values of gamma for this qpoint

!    Here gkq_1band is indexed with real phonon modes (not atom+idir)
!    turn gkq_1band to atom+cartesian coordinates (instead of normal coordinates for qpoint)
     tmp_mat2(:,:,:) = zero
     do ibranch =1,elph_ds%nbranch
      do jbranch =1,elph_ds%nbranch
!      do kbranch=1,elph_ds%nbranch
!      tmp_mat2(1,ibranch,jbranch) = eigvec(1,ibranch,kbranch) * &
!      &  gkq_1band(1,kbranch,jbranch)
!      tmp_mat2(2,ibranch,jbranch) = eigvec(2,ibranch,kbranch) * &
!      &  gkq_1band(1,kbranch,jbranch)
!      end do
       tmp_mat2(1,ibranch,jbranch) = tmp_mat2(1,ibranch,jbranch) + &
&       eigvec(1,ibranch,jbranch) * gkq_1band(1,jbranch,jbranch)
       tmp_mat2(2,ibranch,jbranch) = tmp_mat2(2,ibranch,jbranch) + &
&       eigvec(2,ibranch,jbranch) * gkq_1band(1,jbranch,jbranch)
      end do
     end do
!    DEBUG
!    write (*,'(a,9(2E16.6,1x))') 'normsq_gkq : tmp_mat2 = ', tmp_mat2
!    ENDDEBUG
     gkq_1band(:,:,:) = zero

!    OLDVERSION
!    do ibranch =1,elph_ds%nbranch
!    
!    do jbranch =1,elph_ds%nbranch
!    do kbranch=1,elph_ds%nbranch
!    gkq_1band(1,ibranch,jbranch) =   gkq_1band(1,ibranch,jbranch)&
!    &                    +tmp_mat2(1,ibranch,kbranch) * eigvec(1,jbranch,kbranch)&
!    &                    +tmp_mat2(2,ibranch,kbranch) * eigvec(2,jbranch,kbranch)
!    
!    ! Again, here eigvec is transposed and complexconjugated.
!    gkq_1band(2,ibranch,jbranch) =   gkq_1band(2,ibranch,jbranch)&
!    &                    -tmp_mat2(1,ibranch,kbranch) * eigvec(2,jbranch,kbranch)&
!    &                    +tmp_mat2(2,ibranch,kbranch) * eigvec(1,jbranch,kbranch)
!    end do
!    end do
!    
!    end do
!    ENDOLDVERSION

!    here eigvec is transposed and complexconjugated.
     zgemm_tmp_mat=zero
     call zgemm('n','c',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&     tmp_mat2,elph_ds%nbranch,eigvec,elph_ds%nbranch,czero,zgemm_tmp_mat,elph_ds%nbranch)

     gkq_1band = zgemm_tmp_mat

!    gamma matrix contribution in cartesian coordinates (ie interpolatable form)
     gkk_qpt_tmp(:,ibeff,:,iFSkpt,isppol) = gkk_qpt_tmp(:,ibeff,:,iFSkpt,isppol) &
&     + reshape(gkq_1band,(/2,elph_ds%nbranch*elph_ds%nbranch/))

    end do
   end do
!  END loop over bands ib1 ib2

  end do
! END loop over FSkpt
 end do
!END loop over nsppol


end subroutine nmsq_gam
!!***
