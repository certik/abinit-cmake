!{\src2tex{textfont=tt}}
!!****f* ABINIT/ftgam
!!
!! NAME
!! ftgam
!!
!! FUNCTION
!! If qtor=1 (q->r):
!!  Generates the Fourier transform of the recip space gkk matrices
!!  to obtain the real space ones.
!! If qtor=0 (r->q):
!!  Generates the Fourier transform of the real space gkk matrices
!!  to obtain the reciprocal space ones.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gprim(3,3)= Normalized coordinates in reciprocal space
!! natom= Number of atoms in the unit cell
!! nqpt= Number of q points in the Brillouin zone
!!           if qtor=0 this number is read in the input file
!! nrpt= Number of R points in the Big Box
!! qtor= ( q to r : see above )
!! rpt(3,nprt)= Canonical coordinates of the R points in the unit cell
!!           These coordinates are normalized (=> * acell(3)!!)
!! spqpt(3,nqpt)= Reduced coordinates of the q vectors in reciprocal space
!!           if qtor=0 these vectors are read in the input file
!! wghatm(natom,natom,nrpt)
!!         = Weights associated to a pair of atoms and to a R vector
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!! Input/output
!! gam_qpt(2,3*natom*3*natom,nqpt)
!!  = gamma matrices in recip space coming from the Derivative Data Base
!! gam_rpt(2,3*natom*3*natom,nrpt)
!!  = gamma matrices in real space stored in file unit_gkk_rpt
!!
!! PARENTS
!!      integrate_gamma,integrate_gamma_tr,mka2f,mka2f_tr,mkph_linwid
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! NOTES
!!   copied from ftiaf9.f
!!   recip to real space: real space is forced to disk file unit_gkk_rpt
!!                        recip space depends on gkqwrite and unitgkq3
!!   real to recip space: real space is forced to disk file unit_gkk_rpt
!!                        recip space is necessarily in memory in gkk_qpt
!!
!!    real space elements are complex, but could be reduced, as (-r) = (+r)*
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ftgam (wghatm,gam_qpt,gam_rpt,gprim,natom,nqpt,nrpt,qtor,rpt,spqpt)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: natom,nqpt,nrpt,qtor
!arrays
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),spqpt(3,nqpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)
 real(dp),intent(inout) :: gam_qpt(2,3*natom*3*natom,nqpt)
 real(dp),intent(inout) :: gam_rpt(2,3*natom*3*natom,nrpt)

!Local variables -------------------------
!scalars
 integer :: ia,iatom,ib,ib1,ib2,idir,ii,ip,iqpt,irpt,jatom,jdir,jp,mu,nu
 real(dp) :: asri,asrr,cpu,facti,factr,im,kr,re,wall
 character(len=500) :: message
!arrays
 real(dp) :: coskr(nqpt,nrpt),ftwght(2,3*natom*3*natom),kk(3),sinkr(nqpt,nrpt)
 real(dp) :: tmpgkkrim(nrpt),tmpgkkrre(nrpt)

! *********************************************************************

!DEBUG
!write(6,*)' ftgam : enter '
!write(6,*)' ftgam : qtor ',qtor
!ENDDEBUG

!prepare the phase factors
 do iqpt=1,nqpt
! Calculation of the k coordinates in Normalized Reciprocal
! coordinates
  kk(1)=   spqpt(1,iqpt)*gprim(1,1)+&
&  spqpt(2,iqpt)*gprim(1,2)+&
&  spqpt(3,iqpt)*gprim(1,3)
  kk(2)=   spqpt(1,iqpt)*gprim(2,1)+&
&  spqpt(2,iqpt)*gprim(2,2)+&
&  spqpt(3,iqpt)*gprim(2,3)
  kk(3)=   spqpt(1,iqpt)*gprim(3,1)+&
&  spqpt(2,iqpt)*gprim(3,2)+&
&  spqpt(3,iqpt)*gprim(3,3)
  do irpt=1,nrpt
!  Product of k and r
   kr =        kk(1)*rpt(1,irpt)+&
&   kk(2)*rpt(2,irpt)+&
&   kk(3)*rpt(3,irpt)
   coskr(iqpt,irpt)=cos(two_pi*kr)
   sinkr(iqpt,irpt)=sin(two_pi*kr)
!  if (iqpt < 1000 .and. (irpt == 101 .or. irpt == 901)) then
!  write (*,*) iqpt,irpt,kk,rpt(:,irpt),coskr(iqpt,irpt), sinkr(iqpt,irpt)
!  end if
  end do
 end do



!Recip to real space
 if (qtor==1) then
! 
  gam_rpt(:,:,:) = zero
  do irpt=1,nrpt
!  DEBUG
!  write (*,*) ' ftgam : G->R irpt = ',irpt,' / ',nrpt
!  ENDDEBUG
   do iqpt=1,nqpt
!   write (*,*) iqpt
!   Get the phase factor with normalization!
    re=coskr(iqpt,irpt)/nqpt
    im=sinkr(iqpt,irpt)/nqpt
    do ip=1,3*natom*3*natom
!    Real and imaginary part of the real-space gam matrices
     gam_rpt(1,ip,irpt) = gam_rpt(1,ip,irpt)&
&     +re*gam_qpt(1,ip,iqpt) &
&     +im*gam_qpt(2,ip,iqpt)
     gam_rpt(2,ip,irpt) = gam_rpt(2,ip,irpt)&
&     +re*gam_qpt(2,ip,iqpt) &
&     -im*gam_qpt(1,ip,iqpt)
    end do
   end do
  end do

! Recip space from real space
 else if (qtor==0) then

! write (*,*) 'ftgam : shape(gam_qpt) = ', shape(gam_qpt)
  gam_qpt(:,:,:)=zero

  do irpt=1,nrpt
   do iqpt=1,nqpt

    do iatom=1,natom
     do jatom=1,natom
      re = coskr(iqpt,irpt)*wghatm(iatom,jatom,irpt)
      im = sinkr(iqpt,irpt)*wghatm(iatom,jatom,irpt)

      do idir=1,3
       do jdir=1,3
!       Get phase factor

        ip= jdir + (jatom-1)*3 + (idir-1)*3*natom + (iatom-1)*9*natom
!       Real and imaginary part of the interatomic forces
        gam_qpt(1,ip,iqpt)=&
&        gam_qpt(1,ip,iqpt)&
&        +re*gam_rpt(1,ip,irpt)&
&        -im*gam_rpt(2,ip,irpt)
!       !DEBUG
        gam_qpt(2,ip,iqpt)=&
&        gam_qpt(2,ip,iqpt)&
&        +im*gam_rpt(1,ip,irpt)&
&        +re*gam_rpt(2,ip,irpt)
!       !ENDDEBUG

!       if (iqpt < 100 .and. irpt < 100 .and. &
!       &   tmpgkkrim(irpt)**2+tmpgkkrre(irpt)**2 > tol6) then
!       write (*,'(2I4,2E16.8,x,2E16.8)') &
!       &   iqpt,irpt,re,im,tmpgkkrre(irpt),tmpgkkrim(irpt)
!       end if
!       write (*,'(a)') ' ftgam :gam_qpt :'
!       write (*,'(4E16.5)') gam_qpt(ip,1:nqpt)
       end do
      end do
!     end idir
     end do
    end do ! end iatom

   end do
!  end iqpt
  end do
! end irpt


! There is no other space to Fourier transform from ??
 else
  write(message,'(a,a,a,a,a,i4,a)' )&
&  ' ftiaf9 : BUG - ',ch10,&
&  '  The only allowed values for qtor are 0 or 1, while',ch10,&
&  '  qtor=',qtor,' has been required.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

end subroutine ftgam
!!***
