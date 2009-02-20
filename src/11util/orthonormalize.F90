!{\src2tex{textfont=tt}}
!!****f* ABINIT/orthonormalize
!! NAME
!! orthonormalize
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (GZ,AR,MT)
!! this file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      lobpcgIIwf,lobpcgwf,pw_orthon
!!
!! CHILDREN
!!      dgemm,dpotrf,dtrsm,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine orthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,vectsize
 type(mpi_type) :: mpi_enreg
!arrays
 real(dp) :: blockvectorbx(vectsize,blocksize),blockvectorx(vectsize,blocksize)
 real(dp) :: sqgram(blocksize,blocksize)

!Local variables-------------------------------
#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DGEMM' :: dgemm
!DEC$ ATTRIBUTES ALIAS:'DPOTRF' :: dpotrf
!DEC$ ATTRIBUTES ALIAS:'DTRSM' :: dtrsm
#endif
!scalars
 integer :: iblocksize,ierr,info,jblocksize,old_paral_level,spaceComm
!arrays
 real(dp) :: tsec(2)

! *********************************************************************

 call dgemm('t','n',blocksize,blocksize,vectsize,one,blockvectorx,&
& vectsize,blockvectorbx,vectsize,zero,sqgram,blocksize)
 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm)
 if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
 call xsum_mpi(sqgram,spaceComm,ierr)
 mpi_enreg%paral_level= old_paral_level
 call dpotrf('u',blocksize,sqgram,blocksize,info)

 if (info /= 0 )  then
  write(6,*)'WARNING in dpotrf, info=',info
 endif 
 call dtrsm('r','u','n','n',vectsize,blocksize,one,sqgram,blocksize,&
& blockvectorx,vectsize)

end subroutine orthonormalize
!!***
