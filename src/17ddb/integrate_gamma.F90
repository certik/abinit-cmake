!{\src2tex{textfont=tt}}
!!****f* ABINIT/integrate_gamma
!!
!! NAME
!! integrate_gamma
!!
!! FUNCTION
!! This routine integrates the electron phonon coupling matrix
!! over the kpoints on the fermi surface. A dependency on qpoint
!! or irpt (real space) remains, for gamma_qpt and gamma_rpt resp.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!   FSfullpqtofull = mapping of k+q to k
!!   gprim = reciprocal space lattice vectors
!!   n0 = DOS at fermi level eventually for both spins
!!   natom = number of atoms
!!   nrpt = number of real space points for FT
!!   rpt = coordinates of real space points for FT
!!   spqpt = qpoint coordinates
!!   wghatm = weights for FT of real-space points
!!
!! OUTPUT
!!   elph_ds = modified elph_ds%gamma_qpt and created elph_ds%gamma_rpt
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgam,leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine integrate_gamma(elph_ds,FSfullpqtofull,gprim,n0,natom,nrpt,rpt,spqpt,wghatm)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_17ddb, except_this_one => integrate_gamma
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
 real(dp),intent(in) :: gprim(3,3),n0(elph_ds%nsppol),rpt(3,nrpt)
 real(dp),intent(in) :: spqpt(3,elph_ds%nqpt),wghatm(natom,natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,iFSkptq,ib1,ib2,ibeff,ierr,iqpt,irpt,isppol,qtor
 character(len=500) :: message
 character(len=fnlen) :: fname
!arrays
 real(dp),allocatable :: tmp_gkk(:,:,:,:,:)

! *************************************************************************

 allocate(elph_ds%gamma_qpt(2,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqpt),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' integrate_gamma : ERROR- ',ch10,&
&  ' trying to allocate array elph_ds%gamma_qpt '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 elph_ds%gamma_qpt(:,:,:,:) = zero

 allocate (elph_ds%gamma_rpt(2,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,nrpt),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' integrate_gamma : ERROR- ',ch10,&
&  ' trying to allocate array elph_ds%gamma_rpt '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 elph_ds%gamma_rpt(:,:,:,:) = zero

 if (elph_ds%gkqwrite == 0) then
  write (message,'(a)')' integrate_gamma : keeping gamma matrices in memory'
  call wrtout(06,message,'COLL')
! NOTE: if ngkkband==1 we are using trivial weights since average
! over bands was done in normsq_gkk (nmsq_gam_sumFS or nmsq_pure_gkk)

  do iqpt=1,elph_ds%nqpt
   do isppol=1,elph_ds%nsppol
    do iFSkpt=1,elph_ds%nFSkpt
     iFSkptq = FSfullpqtofull(iFSkpt,iqpt)

     do ib1=1,elph_ds%ngkkband
      do ib2=1,elph_ds%ngkkband
       ibeff = ib2+(ib1-1)*elph_ds%ngkkband
       elph_ds%gamma_qpt(:,:,isppol,iqpt) = elph_ds%gamma_qpt(:,:,isppol,iqpt) + &
&       elph_ds%gkk_qpt(:,ibeff,:,iFSkpt,isppol,iqpt)&
&       *elph_ds%gkk_intweight(ib1,iFSkpt,isppol)*elph_ds%gkk_intweight(ib2,iFSkptq,isppol)
      end do
     end do

    end do
   end do
  end do
 else if (elph_ds%gkqwrite == 1) then

  fname=trim(elph_ds%elph_base_name) // '_GKKQ'
  write (message,'(2a)')' integrate_gamma : reading gamma matrices from file ',trim(fname)
  call wrtout(06,message,'COLL')

  allocate (tmp_gkk (2,elph_ds%ngkkband*elph_ds%ngkkband,elph_ds%nbranch*elph_ds%nbranch,&
&  elph_ds%nFSkpt,elph_ds%nsppol),stat=ierr)
  if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma : ERROR- ',ch10,&
&   ' trying to allocate array tmp_gkk '
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  do iqpt=1,elph_ds%nqpt
   read (elph_ds%unitgkq,REC=iqpt) tmp_gkk
   do isppol=1,elph_ds%nsppol
    do iFSkpt=1,elph_ds%nFSkpt
     iFSkptq = FSfullpqtofull(iFSkpt,iqpt)

     do ib1=1,elph_ds%ngkkband
      do ib2=1,elph_ds%ngkkband
       ibeff = ib2+(ib1-1)*elph_ds%ngkkband
       elph_ds%gamma_qpt(:,:,isppol,iqpt) = elph_ds%gamma_qpt(:,:,isppol,iqpt) + &
&       tmp_gkk(:,ibeff,:,iFSkpt,isppol)&
&       *elph_ds%gkk_intweight(ib1,iFSkpt,isppol)*elph_ds%gkk_intweight(ib2,iFSkptq,isppol)
      end do
     end do

    end do
   end do
!  DEBUG
!  if (iqpt == 1) then
!  write (102,*) ' tmp_gkk ====  in integrate_gamma '
!  do iFSkpt=1,elph_ds%nFSkpt
!  write (102,*) iFSkpt
!  write (102,'(9(2E14.6,1x))') tmp_gkk(:,:,:,iFSkpt,isppol)
!  end do
!  end if
!  ENDDEBUG
  end do
  deallocate (tmp_gkk)
 else
  write (message,'(3a,i3)')' integrate_gamma : BUG-',ch10,&
&  ' Wrong value for gkqwrite = ',elph_ds%gkqwrite
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!need prefactor of 1/nkpt for each integration over 1 kpoint index.
!NOT INCLUDED IN elph_ds%gkk_intweight
 do iqpt=1,elph_ds%nqpt
! elph_ds%gamma_qpt(:,:,:,iqpt) = elph_ds%gamma_qpt(:,:,:,iqpt) / elph_ds%nFSkpt / n0(1) / n0(1)
! elph_ds%gamma_qpt(:,:,:,iqpt) = elph_ds%gamma_qpt(:,:,:,iqpt) / elph_ds%nFSkpt / elph_ds%nFSkpt
  elph_ds%gamma_qpt(:,:,:,iqpt) = elph_ds%gamma_qpt(:,:,:,iqpt) / elph_ds%nFSkpt
 end do

!DEBUG
!write (100,*) ' gamma_qpt on qpts actually calculated ====  in integrate_gamma '
!do iqpt=1,elph_ds%nqpt
!write (100,*) iqpt
!write (100,'(9(2E14.6,1x))') elph_ds%gamma_qpt(:,:,:,iqpt)
!end do
!ENDDEBUG

!Now FT to real space too
 write (message,'(a)')' integrate_gamma : Fourier transforming gamma matrices to real space'
 call wrtout(06,message,'COLL')

 qtor = 1 ! q --> r
 do isppol=1,elph_ds%nsppol
  call ftgam(wghatm,elph_ds%gamma_qpt(:,:,isppol,:),elph_ds%gamma_rpt(:,:,isppol,:),gprim,natom,&
&  elph_ds%nqpt,nrpt,qtor,rpt,spqpt)
 end do

 write (message,'(a)')' integrate_gamma : gamma matrices are in real space '
 call wrtout(06,message,'COLL')


end subroutine integrate_gamma
!!***
