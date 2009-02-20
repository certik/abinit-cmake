!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawfgrtab_utils
!! NAME
!! pawfgrtab_utils
!!
!! FUNCTION
!! This module (?) contains functions used to manipulate
!! variables of the structured datatype pawfgrtab_type.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  
!!
!! OUTPUT
!!  
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

subroutine nullify_pawfgrtab(this)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawfgrtab_type),intent(inout) :: this(:)

!Local variables-------------------------------
 !character(len=500) :: message      
!scalars
 integer :: ia,natom

! *************************************************************************

 natom=SIZE(this)
 do ia=1,natom
  nullify(this(ia)%ifftsph)
  nullify(this(ia)%gylm)
  nullify(this(ia)%gylmgr)
  nullify(this(ia)%gylmgr2)
  nullify(this(ia)%rfgd)
  nullify(this(ia)%vlocgr)
 end do

end subroutine nullify_pawfgrtab
!!***

!!****f* ABINIT/destroy_pawfgrtab 
!! NAME
!! destroy_pawfgrtab
!!
!! FUNCTION
!!  deallocate all associated pointers in the datatype
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine destroy_pawfgrtab(this)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!arrays
 type(pawfgrtab_type),intent(inout) :: this(:)

!Local variables-------------------------------
 !character(len=500) :: message      
!scalars
 integer :: ia,natom

! *************************************************************************

 natom=SIZE(this)
 do ia=1,natom
  if (associated(this(ia)%ifftsph)) deallocate(this(ia)%ifftsph)
  if (associated(this(ia)%gylm   )) deallocate(this(ia)%gylm   )
  if (associated(this(ia)%gylmgr )) deallocate(this(ia)%gylmgr )
  if (associated(this(ia)%gylmgr2)) deallocate(this(ia)%gylmgr2)
  if (associated(this(ia)%rfgd   )) deallocate(this(ia)%rfgd   )
  if (associated(this(ia)%vlocgr )) deallocate(this(ia)%vlocgr )
 end do

end subroutine destroy_pawfgrtab 
!!***

!!****f* ABINIT/init_pawfgrtab 
!! NAME
!! init_pawfgrtab
!!
!! FUNCTION
!!  initialize a pawfgrtab datatype
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine init_pawfgrtab(this,l_size_atm)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: l_size_atm(:)
 type(Pawfgrtab_type),intent(inout) :: this(:)

!Local variables-------------------------------
!scalars
 integer :: ia,natom
 character(len=500) :: msg

! *************************************************************************

 natom=SIZE(this)
 if (natom/=SIZE(l_size_atm)) then 
  write(msg,'(4a)')ch10,&
&  ' init_pawfgrtab : BUG- ',ch10,&
&  ' Sizes of assumed shape arrays do not match ' 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 do ia=1,natom
! pawfgrtab(ia)%l_size=pawtab(dtset%typat(ia))%l_size
  this(ia)%l_size=l_size_atm(ia)
  this(ia)%nfgd=0               ; allocate(this(ia)%ifftsph(0))
  this(ia)%rfgd_allocated=0     ; allocate(this(ia)%rfgd(0,0))
  this(ia)%gylm_allocated=0     ; allocate(this(ia)%gylm(0,0))
  this(ia)%gylmgr_allocated=0   ; allocate(this(ia)%gylmgr(0,0,0))
  this(ia)%gylmgr2_allocated=0  ; allocate(this(ia)%gylmgr2(0,0,0))
  this(ia)%vlocgr_allocated=0   ; allocate(this(ia)%vlocgr(0,0))
 end do

end subroutine init_pawfgrtab 
!!***

!!****f* ABINIT/print_pawfgrtab 
!! NAME
!! print_pawfgrtab
!!
!! FUNCTION
!!  print the content of a pawfgrtab datatype
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! mode_paral(optional)=either "COLL" or "PERS"
!! this<pawfgrtab_type>= the datatype to be printed 
!! unitno(optional)=unit number for output
!! prtvol(optional)=verbosity level
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine print_pawfgrtab(this,unitno,prtvol,mode_paral)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in),optional :: prtvol,unitno
 character(len=4),intent(in),optional :: mode_paral
!arrays
 type(pawfgrtab_type),intent(inout) :: this(:)

!Local variables-------------------------------
!scalars
 integer :: ia,natom,unt,verb
 character(len=4) :: mode
 character(len=500) :: msg

! *************************************************************************

 verb=0      ; if (PRESENT(prtvol))     verb=prtvol
 unt=std_out ; if (PRESENT(unitno))     unt=unitno
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral

 natom=SIZE(this)

 write(msg,'(3a)')ch10,' === Content of pawfgrtab datatype === ',ch10
 call wrtout(unt,msg,mode)
 do ia=1,natom
  write(msg,'(3(2a,i5))')ch10,&
&  ' > For atom number : ',ia,ch10,&
&  '    1+ Max l in Gaunt coefficients ',this(ia)%l_size,ch10,&
&  '    Number of fine FFT points in PAW sphere ',this(ia)%nfgd               
  call wrtout(unt,msg,mode)

  if (verb>=3) then 
   write(msg,'(5(2a,i2))')ch10,&
&   '    rfgd_allocated    : ',this(ia)%rfgd_allocated,ch10,&
&   '    gylm_allocated    : ',this(ia)%gylm_allocated,ch10,&   
&   '    gylmgr_allocated  : ',this(ia)%gylmgr_allocated,ch10,& 
&   '    gylmgr2_allocated : ',this(ia)%gylmgr2_allocated,ch10,&
&   '    vlocgr_allocated  : ',this(ia)%vlocgr_allocated 
   call wrtout(unt,msg,mode)
  end if

! These huge arrays are not printed out!
! this(ia)%ifftsph
! this(ia)%rfgd
! this(ia)%gylm
! this(ia)%gylmgr
! this(ia)%gylmgr2
! this(ia)%vlocgr
 end do

end subroutine print_pawfgrtab 
!!***
