!{\src2tex{textfont=tt}}
!!****f* ABINIT/nmsq_gam_sumfs
!!
!! NAME
!! nmsq_gam_sumfs
!!
!! FUNCTION
!!  Calculate gamma matrices from original h1_mat_el_sq matrix
!!  elements averaging over bands near the Fermi surface
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
!!   eigvec = eigenvectors of phonons (to turn to cartesian coord frame)
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
!!   gkk_qpt_tmp = tmp matrix for all gamma matrix elements, saved to disk or to memory in nmsq_gam_sumFS
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

subroutine nmsq_gam_sumFS(accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
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
& h1_mat_el_sq(2,elph_ds%nFSband*elph_ds%nFSband*elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
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
 real(dp) :: gkq(3),gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch),kpt(3)
 real(dp) :: redkpt(3),tmp_gkq_sum_bands(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp) :: tmp_mat2(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),allocatable :: matrx(:,:),val(:),vec(:,:,:),zgemm_tmp_mat(:,:,:)
 real(dp),allocatable :: zhpev1(:,:),zhpev2(:)

! *************************************************************************

#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'ZGEMM' :: zgemm
#endif

 if (elph_ds%tkeepbands /= 0) then
  write (message,'(3a)')' nmsq_gam_sumfs : BUG- ',ch10,&
&  ' elph_ds%tkeepbands should be 0 in order to average over bands!'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!MG20060603 NOTE:
!accum_mat and accum_mat2 are real, the imaginary part is used for debugging purpose
!accum_mat2 is used to store the phonon-linewidhts before interpolation

 do isppol=1,elph_ds%nsppol
  do iFSkpt=1,elph_ds%nFSkpt

   iFSkptq = FSfullpqtofull(iFSkpt,iqptfull)

   gkq_sum_bands = zero
   tmp_gkq_sum_bands = zero


   do ib1=1,elph_ds%nFSband
!   weights for distance from the fermi surface
    sd1 = FSintweight(ib1,iFSkpt,isppol)

    do ib2=1,elph_ds%nFSband
!    weights for distance from the fermi surface
     sd2 = FSintweight(ib2,iFSkptq,isppol)
     ibeff=ib2+(ib1-1)*elph_ds%nFSband

!    OLDVERSION keep temporarily. Remove if 5.2 looks good
!    ! phonon mode index for displ_red
!    do jbranch=1,elph_ds%nbranch
!    
!    gkq(:) = zero
!    ! Do scalar product of displacement vector by reduced perturbations
!    ! loops over atoms and red directions
!    !do ibranch=1,elph_ds%nbranch
!    ! should this be displ_red complex conjugated? Yes from DECAFT
!    !gkq(1) = gkq(1)&
!    !     & + displ_red(1,ibranch,jbranch)*h1_mat_el_sq(1,ib2,ib1,ibranch,iFSkpt) &
!    !     & + displ_red(2,ibranch,jbranch)*h1_mat_el_sq(2,ib2,ib1,ibranch,iFSkpt)
!    !gkq(2) = gkq(2)&
!    !     & + displ_red(1,ibranch,jbranch)*h1_mat_el_sq(2,ib2,ib1,ibranch,iFSkpt) &
!    !     & - displ_red(2,ibranch,jbranch)*h1_mat_el_sq(1,ib2,ib1,ibranch,iFSkpt)
!    !end do
!    do ipert1=1,elph_ds%nbranch
!    do ipert2=1,elph_ds%nbranch
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
!    write (message,'(3a,es16.8)')' nmsq_gam_sumfs : WARNING-  gkq is not real!',ch10,&
!    &             ' Im(gkq) = ',gkq(2)
!    call wrtout(06,message,'COLL')
!    end if
!    
!    !VERSION5.1.3
!    !  Add weights for integration over FS and wf = (1/sqrt(2 omega))
!    !gkq_sum_bands(1,jbranch,jbranch) = gkq_sum_bands(1,jbranch,jbranch) + sd1 * sd2 * wf(jbranch)**2 * gkq(1)
!    !ENDVERSION5.1.3
!    
!    !MG20060607 there is no explicit dependence on omega in the linewidth.
!    !It is better dont use wf at all and employ the same approach as in nmsq_gam or nmsq_pure_gkk
!    gkq_sum_bands(1,jbranch,jbranch) = gkq_sum_bands(1,jbranch,jbranch) + pi*sd1 * sd2 * gkq(1)
!    !ENDMG20060607
!    end do
!    ENDOLDVERSION

     allocate (zgemm_tmp_mat (2,elph_ds%nbranch,elph_ds%nbranch))
     zgemm_tmp_mat=zero
     tmp_mat2 = reshape(h1_mat_el_sq(:,ibeff,:,isppol,iFSkpt),(/2,elph_ds%nbranch,elph_ds%nbranch/))
     call zgemm('c','n',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&     displ_red,elph_ds%nbranch,tmp_mat2,&
&     elph_ds%nbranch,czero,zgemm_tmp_mat,elph_ds%nbranch)

!    MG20060607 there is no explicit dependence on omega in the linewidth.
!    It is better dont use wf at all and employ the same approach as in nmsq_gam or nmsq_pure_gkk
     tmp_mat2=zero
     call zgemm('n','n',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,dcmplx(pi*sd1*sd2,zero),&
&     zgemm_tmp_mat,elph_ds%nbranch,displ_red,&
&     elph_ds%nbranch,czero,tmp_mat2,elph_ds%nbranch)
!    sum over bands in gkq_sum_bands
     do ipert1=1,elph_ds%nbranch
      gkq_sum_bands(1,ipert1,ipert1) = gkq_sum_bands(1,ipert1,ipert1) + tmp_mat2(1,ipert1,ipert1)
     end do
!    DEBUG
!    ! uncomment OLLDVERSION above and present DEBUG to test
!    sum over bands in gkq_sum_bands
!    tmp_gkq_sum_bands(1,:,:) = tmp_gkq_sum_bands(1,:,:) + tmp_mat2(1,:,:)
!    write (*,*) 'diff btw zgemm and manual multiply for displ_red'
!    write (*,*) (tmp_gkq_sum_bands(:,ipert1,ipert1)-gkq_sum_bands(:,ipert1,ipert1),ipert1=1,elph_ds%nbranch)
!    ENDDEBUG

     deallocate (zgemm_tmp_mat)


!    DEBUG
!    write (*,'(a,a,4i8,2E16.6)') ' nmsq_gam_sumFS : ',&
!    &  'ib1,ib2,iFSkpt,iqptfull,gkq_sum_bands(:,1,1) = ',&
!    &   ib1,ib2,iFSkpt,iqptfull,gkq_sum_bands(:,1,1)
!    ENDDEBUG
    end do
   end do
!  END loop over bands

!  DEBUG
!  if (abs(spqpt(1,iqptfull))+abs(spqpt(2,iqptfull)-half)+abs(spqpt(3,iqptfull)-half) < tol6) then
!  write (*,'(a,9(2E18.6,1x))') 'nmsq_gam_sumFS : g2fq as in decaft = ', &
!  &  (gkq_sum_bands(1,jbranch,jbranch),jbranch=1,elph_ds%nbranch)
!  end if
!  ENDDEBUG

!  summing over k points, still diagonal in jbranch
   accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)

!  VERSION5.1.3
!  MG20060607 This part has been commented since wf is not used anymore,
!  while pi has been introduced in the previous loop
!  28 Apr 2004  PI factor added to wf in read_gkk
!  15 jun 2004  put back into wf
!  Add factor of 2 pi omega(branch) to gamma. Still in diagonal normal mode basis
!  do jbranch=1,elph_ds%nbranch
!  !!   gkq_sum_bands(:,jbranch,jbranch) = two * pi * phfrq_tmp(jbranch) * gkq_sum_bands(:,jbranch,jbranch)
!  gkq_sum_bands(:,jbranch,jbranch) = two * pi * abs(phfrq_tmp(jbranch)) * gkq_sum_bands(:,jbranch,jbranch)
!  end do
!  ENDVERSION5.1.3

   accum_mat2(:,:,:,isppol) = accum_mat2(:,:,:,isppol) + gkq_sum_bands(:,:,:)

!  summed over bands, now turn to cartesian coordinates

!  Final Gamma matrix (hermitian) = E * D_g * E^{+}
!  Where E^{+} is the hermitian conjugate of the eigenvector matrix E
!  And D_g is the diagonal matrix of values of gamma for this qpoint

!  Here gkq_sum_bands is indexed with real phonon modes (not atom+idir)
!  turn gkq_sum_bands to atom+cartesian coordinates (instead of normal coordinates for qpoint)
!  This is not a full matrix multiplication, just vector one, by
!  gkq_sum_bands(1,jbranch,jbranch)
   tmp_mat2(:,:,:) = zero
   do ibranch =1,elph_ds%nbranch
    do jbranch =1,elph_ds%nbranch
!    do kbranch=1,elph_ds%nbranch
!    tmp_mat2(1,ibranch,jbranch) = eigvec(1,ibranch,kbranch) * &
!    &  gkq_sum_bands(1,kbranch,jbranch)
!    tmp_mat2(2,ibranch,jbranch) = eigvec(2,ibranch,kbranch) * &
!    &  gkq_sum_bands(1,kbranch,jbranch)
!    end do
     tmp_mat2(1,ibranch,jbranch) = tmp_mat2(1,ibranch,jbranch) + &
&     eigvec(1,ibranch,jbranch) * &
&     gkq_sum_bands(1,jbranch,jbranch)
     tmp_mat2(2,ibranch,jbranch) = tmp_mat2(2,ibranch,jbranch) + &
&     eigvec(2,ibranch,jbranch) * &
&     gkq_sum_bands(1,jbranch,jbranch)
    end do
   end do
!  DEBUG
!  write (*,'(a,9(2E16.6,1x))') 'nmsq_gam_sumFS : tmp_mat2 = ', tmp_mat2
!  ENDDEBUG

!  OLDVERSION keep temporarily. Remove if 5.2 looks good
!  gkq_sum_bands(:,:,:) = zero
!  do ibranch =1,elph_ds%nbranch
!  do jbranch =1,elph_ds%nbranch
!  do kbranch=1,elph_ds%nbranch
!  gkq_sum_bands(1,ibranch,jbranch) =  gkq_sum_bands(1,ibranch,jbranch) &
!  &  +tmp_mat2(1,ibranch,kbranch) * &
!  &        eigvec(1,jbranch,kbranch)&
!  &  +tmp_mat2(2,ibranch,kbranch) * &
!  &        eigvec(2,jbranch,kbranch)
!  ! Again, here eigvec is transposed and complexconjugated.
!  gkq_sum_bands(2,ibranch,jbranch) = gkq_sum_bands(2,ibranch,jbranch) &
!  &  -tmp_mat2(1,ibranch,kbranch) * &
!  &        eigvec(2,jbranch,kbranch)  &
!  &  +tmp_mat2(2,ibranch,kbranch) * &
!  &        eigvec(1,jbranch,kbranch)
!  end do
!  end do
!  end do
!  ENDOLDVERSION

   allocate (zgemm_tmp_mat(2,elph_ds%nbranch,elph_ds%nbranch))
!  here eigvec is transposed and complexconjugated.
   zgemm_tmp_mat=zero
   call zgemm('n','c',elph_ds%nbranch,elph_ds%nbranch,elph_ds%nbranch,cone,&
&   tmp_mat2,elph_ds%nbranch,eigvec,elph_ds%nbranch,czero,zgemm_tmp_mat,elph_ds%nbranch)

!  DEBUG
!  ! uncomment OLLDVERSION above and present DEBUG to test
!  write (*,*) 'diff btw zgemm and manual multiply'
!  write (*,*) ((zgemm_tmp_mat(:,ipert1,ipert2)-gkq_sum_bands(:,ipert1,ipert2),ipert1=1,elph_ds%nbranch),ipert2=1,elph_ds%nbranch)
!  ENDDEBUG
   gkq_sum_bands = zgemm_tmp_mat
   deallocate (zgemm_tmp_mat)

!  DEBUG
!  if (abs(spqpt(1,iqptfull))+abs(spqpt(2,iqptfull)-half)+abs(spqpt(3,iqptfull)-half) < tol6) then
!  write (*,'(a)') 'nmsq_gam_sumFS : gkq_sum_bands after turning as gmatq from decaft = '
!  do ibranch =1,elph_ds%nbranch
!  write (*,'(3(2E18.6,1x))') (gkq_sum_bands(:,jbranch,ibranch),jbranch=1,elph_ds%nbranch)
!  end do
!  
!  ! Diagonalize gamma matrix at qpoint (complex matrix). Copied from phfrq3
!  ier=0
!  ii=1
!  allocate(matrx(2,(elph_ds%nbranch*(elph_ds%nbranch+1))/2))
!  do i2=1,elph_ds%nbranch
!  do i1=1,i2
!  matrx(1,ii)=gkq_sum_bands(1,i1,i2)
!  matrx(2,ii)=gkq_sum_bands(2,i1,i2)
!  ii=ii+1
!  end do
!  end do
!  allocate(zhpev1(2,2*elph_ds%nbranch-1),zhpev2(3*elph_ds%nbranch-2))
!  allocate(val(elph_ds%nbranch),vec(2,elph_ds%nbranch,elph_ds%nbranch))
!  #if defined T3E
!  call CHPEV ('V','U',elph_ds%nbranch,matrx,val,vec,elph_ds%nbranch,zhpev1,&
!  &    zhpev2,ier)
!  #else
!  call ZHPEV ('V','U',elph_ds%nbranch,matrx,val,vec,elph_ds%nbranch,zhpev1,&
!  &    zhpev2,ier)
!  #endif
!  
!  write (*,*) ' nmsq_gam_sumFS : eigenvalues = '
!  write (*,'(3E18.6)') val
!  
!  deallocate(matrx,zhpev1,zhpev2,vec,val)
!  end if
!  ENDDEBUG

!  ! gamma matrix contribution in cartesian coordinates (ie interpolatable form)
!  gamma matrix contribution in reduced coordinates (ie interpolatable form)
   gkk_qpt_tmp(:,1,:,iFSkpt,isppol) = gkk_qpt_tmp(:,1,:,iFSkpt,isppol) &
&   + reshape(gkq_sum_bands(:,:,:),(/2,elph_ds%nbranch*elph_ds%nbranch/))

!  accum_mat(:,:,:,isppol) = accum_mat(:,:,:,isppol) + gkq_sum_bands(:,:,:)
  end do
! END loop over FSkpt
 end do
!END loop over sppol

end subroutine nmsq_gam_sumFS
!!***
