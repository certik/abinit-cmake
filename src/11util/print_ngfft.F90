!{\src2tex{textfont=tt}}
!!****f* ABINIT/print_ngfft
!! NAME
!! print_ngfft
!!
!! FUNCTION
!!  Print the content of ngfft(18) in a nice explicative format
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  unitno(optional)=unit number for output
!!  prtvol(optional)=verbosity level
!!  mode_paral(optional): either "COLL" or "PERS"
!!
!! OUTPUT
!!  Only writing 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  
!!
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_ngfft(ngfft,header,unitno,mode_paral,prtvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unitno
 character(len=*),intent(in),optional :: header
 character(len=4),intent(in),optional :: mode_paral
!arrays
 integer,intent(in) :: ngfft(18)

!Local variables-------------------------------
!scalars
 integer :: unt,verb
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 verb=0      ; if (PRESENT(prtvol))     verb=prtvol
 unt=std_out ; if (PRESENT(unitno))     unt=unitno
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 msg=' ==== Info on the FFT mesh reported in ngfft ==== '
 if (PRESENT(header)) msg=' ==== '//TRIM(ADJUSTL(header))//' ==== '
 call wrtout(unt,msg,mode)
 write(msg,'(2(a,3i5,a),a,i5,2a,i5)')&
& '  FFT mesh divisions ........................ ',ngfft(1),ngfft(2),ngfft(3),ch10,&
& '  Augmented FFT divisions ................... ',ngfft(4),ngfft(5),ngfft(6),ch10,&
& '  FFT algorithm ............................. ',ngfft(7),ch10,&
& '  FFT cache size ............................ ',ngfft(8)
 call wrtout(unt,msg,mode)
 write(msg,'(6(a,i5,a),a,4i5)')&
& '  FFT parallelization level ................. ',ngfft(9),ch10,&
& '  Number of processors in my FFT group ...... ',ngfft(10),ch10,&
& '  Index of me in my FFT group ............... ',ngfft(11),ch10,&
& '  No of xy planes in R space treated by me .. ',ngfft(12),ch10,&
& '  No of xy planes in G space treated by me .. ',ngfft(13),ch10,&
& '  MPI communicator for FFT .................. ',ngfft(14),ch10,&
& '  Value of ngfft(15:18) ..................... ',ngfft(15:18)
 call wrtout(unt,msg,mode)

end subroutine print_ngfft
!!***
