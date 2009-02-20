!{\src2tex{textfont=tt}}
!!****f* ABINIT/integrate_gamma_tr
!!
!! NAME
!! integrate_gamma_tr
!!
!! FUNCTION
!! This routine integrates the TRANSPORT electron phonon coupling matrices
!! over the kpoints on the fermi surface. A dependency on qpoint
!! or irpt (real space) remains, for gamma_qpt_in/out and gamma_rpt_in/out resp.
!! Copied from integrate_gamma
!!
!! COPYRIGHT
!! Copyright (C) 2004-2008 ABINIT group (JPC)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!   elph_ds = elphon datastructure with data and dimensions
!!   FSfullpqtofull = mapping of k+q to k
!!   gprim = reciprocal space lattice vectors
!!   natom = number of atoms
!!   nrpt = number of real space points for FT
!!   rpt = coordinates of real space points for FT
!!   spqpt = qpoint coordinates
!!   wghatm = weights for FT of real-space points
!!
!! OUTPUT
!!   elph_tr_ds%gamma_qpt_trout and created elph_tr_ds%gamma_rpt_trout
!!   elph_tr_ds%gamma_qpt_trin and created elph_tr_ds%gamma_rpt_trin
!!
!! NOTES
!!
!! PARENTS
!!      elphon
!!
!! CHILDREN
!!      ftgam
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine integrate_gamma_tr(elph_ds,FSfullpqtofull,gprim,natom,nrpt,rpt,spqpt,wghatm,elph_tr_ds)

 use defs_basis
  use defs_datatypes
  use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_17ddb, except_this_one => integrate_gamma_tr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nrpt
 type(elph_tr_type) :: elph_tr_ds
 type(elph_type),intent(inout) :: elph_ds
!arrays
 integer,intent(in) :: FSfullpqtofull(elph_ds%nFSkpt,elph_ds%nqpt)
 real(dp),intent(in) :: gprim(3,3),rpt(3,nrpt),spqpt(3,elph_ds%nqpt)
 real(dp),intent(in) :: wghatm(natom,natom,nrpt)

!Local variables-------------------------------
!scalars
 integer :: iFSkpt,iFSkptq,ib1,ib2,ibeff,ierr,iqpt,isppol,qtor
 character(len=500) :: message
!arrays
 real(dp),allocatable :: tmp_gkkin(:,:,:,:,:),tmp_gkkout(:,:,:,:,:)

! *************************************************************************

 allocate(elph_tr_ds%gamma_qpt_trin(2,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqpt),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&  ' trying to allocate array elph_tr_ds%gamma_qpt_trin '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 allocate(elph_tr_ds%gamma_qpt_trout(2,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,elph_ds%nqpt),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&  ' trying to allocate array elph_tr_ds%gamma_qpt_trout '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 elph_tr_ds%gamma_qpt_trin(:,:,:,:) = zero
 elph_tr_ds%gamma_qpt_trout(:,:,:,:) = zero


 allocate (elph_tr_ds%gamma_rpt_trout(2,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,nrpt),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&  ' trying to allocate array elph_tr_ds%gamma_rpt_trout '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 elph_tr_ds%gamma_rpt_trout(:,:,:,:) = zero

 allocate (elph_tr_ds%gamma_rpt_trin(2,elph_ds%nbranch*elph_ds%nbranch,elph_ds%nsppol,nrpt),stat=ierr)
 if (ierr /= 0 ) then
  write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&  ' trying to allocate array elph_tr_ds%gamma_rpt_trin '
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 elph_tr_ds%gamma_rpt_trin(:,:,:,:) = zero

 if (elph_ds%gkqwrite == 0) then
  write (message,'(a)')' integrate_gamma_tr : keeping gamma matrices in memory'
  call wrtout(06,message,'COLL')
! NOTE: if ngkkband==1 we are using trivial weights since average
! over bands was done in normsq_gkk (nmsq_gam_sumFS or nmsq_pure_gkk)

  do iqpt=1,elph_ds%nqpt

   do isppol=1,elph_ds%nsppol
    do iFSkpt=1,elph_ds%nFSkpt
     iFSkptq = FSfullpqtofull(iFSkpt,iqpt)

     do ib1=1,elph_ds%ngkkband
      do ib2=1,elph_ds%ngkkband
       ibeff=ib2+(ib1-1)*elph_ds%ngkkband


!      TODO: generalize gamma_qpt_trin etc to sppol
       elph_tr_ds%gamma_qpt_trin(:,:,isppol,iqpt) = elph_tr_ds%gamma_qpt_trin(:,:,isppol,iqpt) + &
&       elph_tr_ds%gkk_qpt_trin(:,ibeff,:,iFSkpt,isppol,iqpt)&
&       *elph_ds%gkk_intweight(ib1,iFSkpt,isppol)*elph_ds%gkk_intweight(ib2,iFSkptq,isppol)

       elph_tr_ds%gamma_qpt_trout(:,:,isppol,iqpt) = elph_tr_ds%gamma_qpt_trout(:,:,isppol,iqpt) + &
&       elph_tr_ds%gkk_qpt_trout(:,ibeff,:,iFSkpt,isppol,iqpt)&
&       *elph_ds%gkk_intweight(ib1,iFSkpt,isppol)*elph_ds%gkk_intweight(ib2,iFSkptq,isppol)
      end do
     end do
    end do
   end do ! isppol
  end do

 else if (elph_ds%gkqwrite == 1) then
  allocate (tmp_gkkin (2,elph_ds%ngkkband*elph_ds%ngkkband,&
&  elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol),stat=ierr)
  if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&   ' trying to allocate array tmp_gkkin '
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

  allocate (tmp_gkkout (2,elph_ds%ngkkband*elph_ds%ngkkband,&
&  elph_ds%nbranch*elph_ds%nbranch,elph_ds%nFSkpt,elph_ds%nsppol),stat=ierr)
  if (ierr /= 0 ) then
   write (message,'(3a)')' integrate_gamma_tr : ERROR- ',ch10,&
&   ' trying to allocate array tmp_gkkout '
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
  do iqpt=1,elph_ds%nqpt
   read (elph_tr_ds%unitgkq_trin,REC=iqpt) tmp_gkkin     
   read (elph_tr_ds%unitgkq_trout,REC=iqpt) tmp_gkkout     
   
   do isppol=1,elph_ds%nsppol  
    do iFSkpt=1,elph_ds%nFSkpt
     iFSkptq = FSfullpqtofull(iFSkpt,iqpt)
     
     do ib1=1,elph_ds%ngkkband
      do ib2=1,elph_ds%ngkkband
       ibeff=ib2+(ib1-1)*elph_ds%ngkkband
       elph_tr_ds%gamma_qpt_trin(:,:,isppol,iqpt) = elph_tr_ds%gamma_qpt_trin(:,:,isppol,iqpt) + &
&       tmp_gkkin(:,ibeff,:,iFSkpt,isppol)&
&       *elph_ds%gkk_intweight(ib1,iFSkpt,isppol)*elph_ds%gkk_intweight(ib2,iFSkptq,isppol)
       
       elph_tr_ds%gamma_qpt_trout(:,:,isppol,iqpt) = elph_tr_ds%gamma_qpt_trout(:,:,isppol,iqpt) + &
&       tmp_gkkout(:,ibeff,:,iFSkpt,isppol)&
&       *elph_ds%gkk_intweight(ib1,iFSkpt,isppol)*elph_ds%gkk_intweight(ib2,iFSkptq,isppol)
      end do
     end do
    end do ! ik
   end do ! isppol
  end do ! iq
  deallocate (tmp_gkkin)
  deallocate (tmp_gkkout)

 else
  write (message,'(3a,i3)')' integrate_gamma_tr : BUG-',ch10,&
&  ' Wrong value for gkqwrite = ',elph_ds%gkqwrite
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

!need prefactor of 1/nkpt for each integration over 1 kpoint index.
!NOT INCLUDED IN elph_ds%gkk_intweight
 do iqpt=1,elph_ds%nqpt
  elph_tr_ds%gamma_qpt_trout(:,:,:,iqpt) = elph_tr_ds%gamma_qpt_trout(:,:,:,iqpt) / elph_ds%nFSkpt
  elph_tr_ds%gamma_qpt_trin(:,:,:,iqpt) = elph_tr_ds%gamma_qpt_trin(:,:,:,iqpt) / elph_ds%nFSkpt
 end do

!Now FT to real space too
 write (message,'(a)')' integrate_gamma_tr : Fourier transforming transport gamma matrices to real space'
 call wrtout(06,message,'COLL')

 qtor = 1 ! q --> r

 do isppol=1,elph_ds%nsppol
  call ftgam(wghatm,elph_tr_ds%gamma_qpt_trout(:,:,isppol,:),&
&  elph_tr_ds%gamma_rpt_trout(:,:,isppol,:),gprim,natom,&
&  elph_ds%nqpt,nrpt,qtor,rpt,spqpt)

  call ftgam(wghatm,elph_tr_ds%gamma_qpt_trin(:,:,isppol,:),&
&  elph_tr_ds%gamma_rpt_trin(:,:,isppol,:),gprim,natom,&
&  elph_ds%nqpt,nrpt,qtor,rpt,spqpt)
 end do

 write (message,'(2a)')' integrate_gamma_tr : transport gamma matrices are in real space '
 call wrtout(06,message,'COLL')

!write(345,*)elph_tr_ds%gamma_rpt_trout
!write(345,*)elph_tr_ds%gamma_rpt_trin

!if (elph_ds%gkqwrite==0) then
!deallocate (elph_tr_ds%gkk_qpt_trin)
!deallocate (elph_tr_ds%gkk_qpt_trout)
!end if

end subroutine integrate_gamma_tr
!!***
