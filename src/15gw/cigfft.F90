!{\src2tex{textfont=tt}}
!!****f* ABINIT/cigfft
!! NAME
!! cigfft
!!
!! FUNCTION
!! For each of the (2*nG0sh+1)**3 vectors G0 around the origin, 
!! calculate G-G0 and its FFT index number for all the NPWVEC vectors G.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! mG0(3)= For each reduced direction gives the max G0 component to account for umklapp processes
!! gvec(3,npwvec)=coordinates of G vectors
!! npwvec=number of plane waves
!! ngfft(18)=contain all needed information about 3D FFT,
!!        see ~abinit/doc/input_variables/vargs.htm#ngfft
!!
!! OUTPUT
!! igfft(npwvec,2*mG0(1)+1,2*mG0(2)+1,2*mG0(3)+1)=for each G, and each G0 vector, 
!!  gives the FFT grid index of the G-G0 vector
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cigfft(mG0,npwvec,ngfft,gvec,igfft)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: npwvec
!arrays
 integer,intent(in) :: gvec(3,npwvec)
 integer,intent(in) :: mG0(3),ngfft(18)
 integer,intent(out) :: igfft(npwvec,2*mG0(1)+1,2*mG0(2)+1,2*mG0(3)+1)

!Local variables ------------------------------
!scalars
 integer :: gmg01,gmg02,gmg03,ig,ig01,ig02,ig03 
 integer :: n1,n2,n3
 character(len=500) :: msg
!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(2a)')ch10,' cigfft : enter'
 call wrtout(std_out,msg,'COLL')
#endif

 if (ANY(mG0(:)<0)) then
  write(msg,'(4a)')ch10,&
&  ' cigfft : BUG -',ch10,' mG0 < 0 '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 

 n1=ngfft(1) 
 n2=ngfft(2) 
 n3=ngfft(3) 

 do ig01=-mG0(1),mG0(1)
  do ig02=-mG0(2),mG0(2)
   do ig03=-mG0(3),mG0(3)
    do ig=1,npwvec
     gmg01=MODULO(gvec(1,ig)-ig01,n1)
     gmg02=MODULO(gvec(2,ig)-ig02,n2)
     gmg03=MODULO(gvec(3,ig)-ig03,n3)
     igfft(ig,ig01+mG0(1)+1,ig02+mG0(2)+1,ig03+mG0(3)+1) = 1+gmg01+gmg02*n1+gmg03*n1*n2
    end do
   end do
  end do
 end do

#if defined DEBUG_MODE
 write(msg,'(2a)')ch10,' cigfft : G-table set up'
 call wrtout(std_out,msg,'COLL')
#endif

end subroutine cigfft
!!***
