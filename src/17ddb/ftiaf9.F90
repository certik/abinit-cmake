!{\src2tex{textfont=tt}}
!!****f* ABINIT/ftiaf9
!!
!! NAME
!! ftiaf9
!!
!! FUNCTION
!! If qtor=1 (q->r):
!!   Generates the Fourier transform of the dynamical matrices
!!   to obtain interatomic forces (real space).
!! If qtor=0 (r->q):
!!   Generates the Fourier transform of the interatomic forces
!!   to obtain dynamical matrices (reciprocal space).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (JCC,XG)
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
!! dynmat(2,3,natom,3,natom,nqpt)
!!  = Dynamical matrices coming from the Derivative Data Base
!! atmfrc(2,3,natom,3,natom,nrpt)
!!  = Interatomic Forces in real space !!
!!  We used the imaginary part just for debugging !
!!
!! PARENTS
!!      gtdyn9,mkifc9
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ftiaf9(atmfrc,dynmat,gprim,natom,nqpt,&
&                  nrpt,qtor,rpt,spqpt,wghatm)

 use defs_basis


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
 real(dp),intent(inout) :: atmfrc(2,3,natom,3,natom,nrpt)
 real(dp),intent(inout) :: dynmat(2,3,natom,3,natom,nqpt)

!Local variables -------------------------
!scalars
 integer :: ia,ib,ii,iqpt,irpt,mu,nu
 real(dp) :: asri,asrr,facti,factr,im,kr,re
 character(len=500) :: message
!arrays
 real(dp) :: kk(3)

! *********************************************************************

!Interatomic Forces from Dynamical Matrices

 if (qtor==1) then

  atmfrc(:,:,:,:,:,:)=0.0_dp
  do irpt=1,nrpt
   do iqpt=1,nqpt

!   Calculation of the k coordinates in Normalized Reciprocal
!   coordinates
    kk(1)=spqpt(1,iqpt)*gprim(1,1)+spqpt(2,iqpt)*&
&    gprim(1,2)+spqpt(3,iqpt)*gprim(1,3)
    kk(2)=spqpt(1,iqpt)*gprim(2,1)+spqpt(2,iqpt)*&
&    gprim(2,2)+spqpt(3,iqpt)*gprim(2,3)
    kk(3)=spqpt(1,iqpt)*gprim(3,1)+spqpt(2,iqpt)*&
&    gprim(3,2)+spqpt(3,iqpt)*gprim(3,3)

!   Product of k and r
    kr=kk(1)*rpt(1,irpt)+kk(2)*rpt(2,irpt)+kk(3)*&
&    rpt(3,irpt)

!   Get the phase factor
    re=cos(two_pi*kr)
    im=sin(two_pi*kr)

!   Now, big inner loops on atoms and directions
!   The indices are ordered to give better speed
    do ib=1,natom
     do nu=1,3
      do ia=1,natom
       do mu=1,3
!       Real and imaginary part of the interatomic forces
        atmfrc(1,mu,ia,nu,ib,irpt)=atmfrc(1,mu,ia,nu,ib,irpt)&
&        +re*dynmat(1,mu,ia,nu,ib,iqpt)&
&        +im*dynmat(2,mu,ia,nu,ib,iqpt)
!       The imaginary part should be equal to zero !!!!!!
!       atmfrc(2,mu,ia,nu,ib,irpt)=atmfrc(2,mu,ia,nu,ib,irpt)
!       &          +re*dynmat(2,mu,ia,nu,ib,iqpt)
!       &          -im*dynmat(1,mu,ia,nu,ib,iqpt)
       end do
      end do
     end do
    end do

   end do
  end do

! The sum has to be weighted by a normalization factor of 1/nqpt
  atmfrc(:,:,:,:,:,:)=atmfrc(:,:,:,:,:,:)/nqpt

! DEBUG
! do irpt=1,nrpt
! do ib=1,natom
! do nu=1,3
! do ia=1,natom
! do mu=1,3
! if(ia==1 .and. ( (ib==1.and.irpt==53) .or.&
! &                       (ib==2.and.irpt==53) .or.&
! &                       (ib==3.and.irpt==38)      ))then
! write(6, '(5i3,es16.8)' )&
! &        mu,ia,nu,ib,irpt,atmfrc(1,mu,ia,nu,ib,irpt)
! end if
! end do
! end do
! end do
! end do
! end do
! ENDDEBUG

! Dynamical Matrices from Interatomic Forces
 else if (qtor==0) then

  dynmat(:,:,:,:,:,:)=0.0_dp
  do iqpt=1,nqpt
   do irpt=1,nrpt

!   Calculation of the k coordinates in Normalized Reciprocal
!   coordinates
    kk(1)=spqpt(1,iqpt)*gprim(1,1)+spqpt(2,iqpt)*&
&    gprim(1,2)+spqpt(3,iqpt)*gprim(1,3)
    kk(2)=spqpt(1,iqpt)*gprim(2,1)+spqpt(2,iqpt)*&
&    gprim(2,2)+spqpt(3,iqpt)*gprim(2,3)
    kk(3)=spqpt(1,iqpt)*gprim(3,1)+spqpt(2,iqpt)*&
&    gprim(3,2)+spqpt(3,iqpt)*gprim(3,3)

!   Product of k and r
    kr=kk(1)*rpt(1,irpt)+kk(2)*rpt(2,irpt)+kk(3)*&
&    rpt(3,irpt)

!   Get phase factor
    re=cos(two_pi*kr)
    im=sin(two_pi*kr)

!   Inner loop on atoms and directions
    do ib=1,natom
     do ia=1,natom
      if(abs(wghatm(ia,ib,irpt))>1.0d-10)then
       factr=re*wghatm(ia,ib,irpt)
       facti=im*wghatm(ia,ib,irpt)
       do nu=1,3
        do mu=1,3
!        Real and imaginary part of the dynamical matrices
         dynmat(1,mu,ia,nu,ib,iqpt)=dynmat(1,mu,ia,nu,ib,iqpt)&
&         +factr*atmfrc(1,mu,ia,nu,ib,irpt)
!        Atmfrc should be real
!        &       -im*wghatm(ia,ib,irpt)*atmfrc(2,mu,ia,nu,ib,irpt)
         dynmat(2,mu,ia,nu,ib,iqpt)=dynmat(2,mu,ia,nu,ib,iqpt)&
&         +facti*atmfrc(1,mu,ia,nu,ib,irpt)
!        Atmfrc should be real
!        &        +re*wghatm(ia,ib,irpt)*atmfrc(2,mu,ia,nu,ib,irpt)
        end do
       end do
      end if
     end do
    end do
   end do
  end do

! There is no other space to Fourier transform from ??
 else
  write(message,'(a,a,a,a,a,i4,a)' )&
&  ' ftiaf9 : BUG - ',ch10,&
&  '  The only allowed values for qtor are 0 or 1, while',ch10,&
&  '  qtor=',qtor,' has been required.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

end subroutine ftiaf9
!!***
