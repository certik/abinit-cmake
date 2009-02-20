!{\src2tex{textfont=tt}}
!!****f* ABINIT/normsq_gkq
!!
!! NAME
!! normsq_gkq
!!
!! FUNCTION
!! This routine takes the gkq matrix elements for a given qpoint,
!!   does the scalar product with the phonon displacement vector,
!!   squares the gkq matrix elements
!!   multiplies by the appropriate weights and
!!   puts them in a uniform (atom,icart) basis
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
!!   FSfullpqtofull = mapping of k + q to k
!!   FSintweight = FS integration weights for each band and kpt
!!   FSkpt = coordinates of kpoints near to FS
!!   h1_mat_el_sq = matrix elements $<psi_{k+q,m} | H^{1} | psi_{k,n}>$ matrix-squared
!!   iqptfull = index of present qpoint
!!   phfrq_tmp = phonon frequencies
!!   spqpt = array of qpoint coordinates
!!   wf = gkk matrix element weight with $1/\sqrt{2 \omega}$
!!
!! OUTPUT
!!   elph_ds%gkq filled
!!   qdata(elph_ds%nbranch,elph_ds%nsppol,3) = array containing the phonon frequency, the linwidth
!!                              and $\lambda_{q,\nu}$ for the considered phonon mode
!!
!! NOTES
!!
!! PARENTS
!!      read_gkk
!!
!! CHILDREN
!!      chpev,leave_new,nmsq_gam,nmsq_gam_sumfs,nmsq_pure_gkk,wrtout,zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine normsq_gkq(displ_red,eigvec,elph_ds,FSfullpqtofull,&
&    FSintweight,FSkpt,h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf,qdata)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_17ddb, except_this_one => normsq_gkq
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iqptfull
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
 real(dp),intent(in) :: FSintweight(elph_ds%nFSband,elph_ds%nFSkpt)
 real(dp),intent(in) :: FSkpt(3,elph_ds%nFSkpt)
 real(dp),intent(in) :: displ_red(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: eigvec(2,elph_ds%nbranch,elph_ds%nbranch)
 real(dp),intent(in) :: h1_mat_el_sq(2,elph_ds%nFSband,elph_ds%nFSband,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nFSkpt)
 real(dp),intent(in) :: phfrq_tmp(elph_ds%nbranch),spqpt(3,elph_ds%nqpt)
 real(dp),intent(in) :: wf(elph_ds%nbranch)
 real(dp),intent(out) :: qdata(elph_ds%nbranch,elph_ds%nsppol,3)

!Local variables-------------------------------
! real(dp) :: gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol)
!scalars
 integer :: eivec,goodkpq,i1,i2,iFSkpt,ib1,ib2,ibranch,ier,ii,isppol,jbranch
 integer :: kbranch
 real(dp) :: lambda_tot,qphnrm,res,sd1,sd2,ss
 character(len=4) :: dirstr,qptstr
 character(len=500) :: message
 character(len=fnlen) :: xsffilnam
!arrays
 real(dp) :: accum_mat(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp) :: accum_mat2(2,elph_ds%nbranch,elph_ds%nbranch,elph_ds%nsppol)
 real(dp) :: gam_now2(2,elph_ds%nbranch,elph_ds%nbranch),gkq(2),kpt(3)
 real(dp) :: lambda(elph_ds%nsppol),redkpt(3)
 real(dp),allocatable :: gkk_qpt_tmp(:,:,:,:,:),matrx(:,:),val(:),vec(:,:,:)
 real(dp),allocatable :: zhpev1(:,:),zhpev2(:)

! *************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
#endif

 allocate(gkk_qpt_tmp(2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
& elph_ds%nFSkpt,elph_ds%nsppol))

!DEBUG
!write(6,*)' normsq_gkq : enter '
!ENDDEBUG

 gkk_qpt_tmp(:,:,:,:,:) = zero
 accum_mat(:,:,:,:) = zero
 accum_mat2(:,:,:,:) = zero

 if (elph_ds%doscalprod == 1) then
  if (elph_ds%tkeepbands == 0) then
   write (*,*) ' normsq_gkq : calling nmsq_gam_sumFS'
   call nmsq_gam_sumFS (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&   FSintweight,FSkpt,gkk_qpt_tmp,&
&   h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)
  else if (elph_ds%tkeepbands == 1) then
   write (*,*) ' normsq_gkq : calling nmsq_gam'
   call nmsq_gam (accum_mat,accum_mat2,displ_red,eigvec,elph_ds,FSfullpqtofull,&
&   FSintweight,FSkpt,gkk_qpt_tmp,&
&   h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)
  else
   write (message,'(4a,i4)')ch10,' normsq_gkq : BUG- ',ch10,&
   ' Wrong value for elph_ds%tkeepbands = ',elph_ds%tkeepbands
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
 else if (elph_ds%doscalprod == 0) then
! else elph_ds%doscalprod == 0  Interpolate on the pure "matrix of matrix elements"
! and do the scalar products later.
  if (elph_ds%tkeepbands == 0) then
   write (*,*) ' normsq_gkq : calling nmsq_pure_gkk_sumFS'

   call nmsq_pure_gkk_sumFS (accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&
&   FSintweight,FSkpt,gkk_qpt_tmp,h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)
  else if (elph_ds%tkeepbands == 1) then
   write (*,*) ' normsq_gkq : calling nmsq_pure_gkk'

   call nmsq_pure_gkk (accum_mat,accum_mat2,displ_red,elph_ds,FSfullpqtofull,&
&   FSintweight,FSkpt,gkk_qpt_tmp,h1_mat_el_sq,iqptfull,phfrq_tmp,spqpt,wf)
  else
   write (message,'(4a,i4)')ch10,' normsq_gkq : BUG- ',ch10,&
&   ' Wrong value for elph_ds%tkeepbands = ',elph_ds%tkeepbands
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if


 else
  write (message,'(3a,i4)')' normsq_gkq: BUG-',ch10,&
&  ' Wrong value for elph_ds%doscalprod = ',elph_ds%doscalprod
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
!end if flag for doing scalar product now.


!MG: values without the good prefactor
 accum_mat(:,:,:,:) = accum_mat(:,:,:,:)/elph_ds%nFSkpt
!MG: only accum_mat2 contains the line-widhts before the Fourier interpolation
 accum_mat2(:,:,:,:) = accum_mat2(:,:,:,:)/elph_ds%nFSkpt

!MG20060531i
!write e-ph quantities before Fourier interpolation
!save e-ph values in the temporary array qdata that will be copied into elph_ds%qgrid_data

 write (message,'(4a,3es16.6,63a)')ch10,                  &
& ' Phonon linewidths before interpolation ',ch10,        &
& ' Q point = ',spqpt(:,iqptfull),ch10,('=',ii=1,60),ch10,&
& ' Mode          Frequency (Ha)  Linewidth (Ha)  Lambda '
 call wrtout(std_out,message,'COLL')

!TODO: add isppol dependence here as well as in mka2f
 lambda_tot = zero
 do isppol=1,elph_ds%nsppol
  do ii=1,elph_ds%nbranch
   lambda(isppol)=zero
!  MG: the tollerance factor is somehow arbitrary
   if (abs(phfrq_tmp(ii)) > tol10) lambda(isppol)=accum_mat2(1,ii,ii,isppol)/&
&   (pi*elph_ds%n0(isppol)*phfrq_tmp(ii)**2)
   lambda_tot=lambda_tot+lambda(isppol)
   write(message,'(i8,es20.6,2es16.6)' )ii,phfrq_tmp(ii),accum_mat2(1,ii,ii,isppol),lambda(isppol)
   call wrtout(std_out,message,'COLL')
!  save values
   qdata(ii,isppol,1)=phfrq_tmp(ii)
   qdata(ii,isppol,2)=accum_mat2(1,ii,ii,isppol)
   qdata(ii,isppol,3)=lambda(isppol)
  end do !loop over branch
 end do !loop over sppol

!normalize for number of spins
 lambda_tot = lambda_tot / elph_ds%nsppol

 write(message,'(61a,44x,es16.6,62a)' )('=',ii=1,60),ch10,lambda_tot,ch10,('=',ii=1,60),ch10
 call wrtout(std_out,message,'COLL')
!ENDMG20060531

!immediately calculate linewidths:
 write (*,*) 'summed accum_mat = '
 write (*,'(3(2E18.6,1x))') accum_mat(:,:,:,1)
 write (*,*) 'summed accum_mat2 = '
 write (*,'(3(2E18.6,1x))')  (accum_mat2(:,ii,ii,1),ii=1,elph_ds%nbranch)
 write (*,*) 'displ_red  = '
 write (*,'(3(2E18.6,1x))') displ_red

!!DEBUG
 do isppol=1,elph_ds%nsppol

  if (elph_ds%doscalprod == 1) then


!  Diagonalize gamma matrix at qpoint (complex matrix). Copied from phfrq3
   ier=0
   ii=1
   allocate(matrx(2,(elph_ds%nbranch*(elph_ds%nbranch+1))/2))
   do i2=1,elph_ds%nbranch
    do i1=1,i2
     matrx(1,ii)=accum_mat2(1,i1,i2,isppol)
     matrx(2,ii)=accum_mat2(2,i1,i2,isppol)
     ii=ii+1
    end do
   end do
   allocate(zhpev1(2,2*elph_ds%nbranch-1),zhpev2(3*elph_ds%nbranch-2))
   allocate(val(elph_ds%nbranch),vec(2,elph_ds%nbranch,elph_ds%nbranch))
#if defined T3E
   call CHPEV ('V','U',elph_ds%nbranch,matrx,val,vec,elph_ds%nbranch,zhpev1,&
&   zhpev2,ier)
#else
   call ZHPEV ('V','U',elph_ds%nbranch,matrx,val,vec,elph_ds%nbranch,zhpev1,&
&   zhpev2,ier)
#endif

   write (*,*) ' normsq_gkq : accumulated eigenvalues isppol ',isppol, ' = '
   write (*,'(3E18.6)') val

   deallocate(matrx,zhpev1,zhpev2,vec,val)

  else if (elph_ds%doscalprod == 0) then
   gam_now2(:,:,:) = zero

   do jbranch=1,elph_ds%nbranch
    do ibranch=1,elph_ds%nbranch
     do kbranch=1,elph_ds%nbranch
!     gam = displ gam_red displ^{*T}
      gam_now2(1,jbranch,jbranch) = gam_now2(1,jbranch,jbranch) + &
&      displ_red(1,ibranch,jbranch)*accum_mat(1,ibranch,kbranch,isppol)*displ_red(1,kbranch,jbranch) &
&      +displ_red(1,ibranch,jbranch)*accum_mat(2,ibranch,kbranch,isppol)*displ_red(2,kbranch,jbranch) &
&      +displ_red(2,ibranch,jbranch)*accum_mat(1,ibranch,kbranch,isppol)*displ_red(2,kbranch,jbranch) &
&      -displ_red(2,ibranch,jbranch)*accum_mat(2,ibranch,kbranch,isppol)*displ_red(1,kbranch,jbranch)
      gam_now2(2,jbranch,jbranch) = gam_now2(2,jbranch,jbranch)   &
&      +displ_red(2,ibranch,jbranch)*accum_mat(2,ibranch,kbranch,isppol)*displ_red(2,kbranch,jbranch) &
&      +displ_red(2,ibranch,jbranch)*accum_mat(1,ibranch,kbranch,isppol)*displ_red(1,kbranch,jbranch) &
&      +displ_red(1,ibranch,jbranch)*accum_mat(2,ibranch,kbranch,isppol)*displ_red(1,kbranch,jbranch) &
&      -displ_red(1,ibranch,jbranch)*accum_mat(1,ibranch,kbranch,isppol)*displ_red(2,kbranch,jbranch)
     end do
    end do
   end do

   write (*,*) ' normsq_gkq : accumulated eigenvalues isppol ', isppol, ' = '
   write (*,'(3(E14.6,1x))') (gam_now2(1,jbranch,jbranch), jbranch=1,elph_ds%nbranch)
   write (*,*) ' normsq_gkq : imag part = '
   write (*,'(3(E14.6,1x))') (gam_now2(2,jbranch,jbranch), jbranch=1,elph_ds%nbranch)

  end if
 end do ! isppol
!!ENDDEBUG

!save gkk_qpt_tmp, eventually to disk
 if (elph_ds%gkqwrite == 0) then
  elph_ds%gkk_qpt(:,:,:,:,:,iqptfull) = gkk_qpt_tmp(:,:,:,:,:)
 else
! rewind (elph_ds%unitgkq)
! do iqpt=1,iqptfull-1
! read (elph_ds%unitgkq) gkk_qpt_tmp2
! end do
  write (*,'(a,i4,2(3E16.6,3x))') 'gkk_qpt_tmp (band1,band1)', &
&  iqptfull, gkk_qpt_tmp(:,1,:,1,1:10)

! write all kpoints to disk
  write (elph_ds%unitgkq,REC=iqptfull) gkk_qpt_tmp
 end if

!DEBUG
!write(6,*)' normsq_gkq : exit '
!ENDDEBUG

end subroutine normsq_gkq
!!***
