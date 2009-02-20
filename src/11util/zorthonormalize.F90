!{\src2tex{textfont=tt}}
!!****f* ABINIT/zorthonormalize
!! NAME
!! zorthonormalize
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
!!      lobpcgccIIwf,lobpcgccwf,pw_orthon
!!
!! CHILDREN
!!      timab,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine zorthonormalize(blockvectorx,blockvectorbx,blocksize,mpi_enreg,sqgram,vectsize)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,vectsize
 type(mpi_type) :: mpi_enreg
!arrays
 complex(dpc) :: blockvectorbx(vectsize,blocksize)
 complex(dpc) :: blockvectorx(vectsize,blocksize),sqgram(blocksize,blocksize)

!Local variables-------------------------------
#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZGEMM' :: zgemm
!DEC$ ATTRIBUTES ALIAS:'ZPOTRF' :: zpotrf
!DEC$ ATTRIBUTES ALIAS:'ZTRSM' :: ztrsm
#endif
!scalars
 integer :: iblocksize,ierr,info,jblocksize,old_paral_level,spaceComm
!arrays
 real(dp) :: tsec(2)

! *********************************************************************

 call zgemm('c','n',blocksize,blocksize,vectsize,cone,blockvectorx,&
& vectsize,blockvectorbx,vectsize,czero,sqgram,blocksize)
 old_paral_level= mpi_enreg%paral_level
 mpi_enreg%paral_level=3
 call xcomm_init(mpi_enreg,spaceComm)
 if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%commcart
 call xsum_mpi(sqgram,spaceComm,ierr)
 mpi_enreg%paral_level= old_paral_level
 call zpotrf('u',blocksize,sqgram,blocksize,info)

 if (info /= 0 )  then
  write(6,*)'WARNING in zpotrf, info=',info
 endif 
 call ztrsm('r','u','n','n',vectsize,blocksize,cone,sqgram,blocksize,&
& blockvectorx,vectsize)

end subroutine zorthonormalize
!!***
