!{\src2tex{textfont=tt}}
!!****f* ABINIT/read_wfrspa
!! NAME
!! read_wfrspa
!!
!! FUNCTION
!!  Internal subroutine used to fill wfrspa from the scratch file when we
!!  have mkmem==0.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (to be completed)
!!
!! OUTPUT
!!  (to be completed)
!! PARENTS
!!      tddft
!!
!! CHILDREN
!!      int2char4
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine read_wfrspa(state,dtfil,eigbnd,iband,isppol,imkmem,ndiel4,ndiel5,ndiel6,wfrspa_extract)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments -----------------------------
 !real(dp), intent(inout) :: wfrspa(:,:,:,:)
!scalars
 integer :: iband,isppol
 integer,intent(in) :: imkmem,ndiel4,ndiel5,ndiel6,state
 real(dp),intent(inout) :: eigbnd
 type(datafiles_type),intent(in) :: dtfil
!arrays
 real(dp),intent(inout) :: wfrspa_extract(ndiel4,ndiel5,ndiel6)

!Local variables---------------------------------
!scalars
 character(len=4) :: tag
 character(len=fnlen) :: fname_wf
!arrays
 integer,save :: states(4)=(/0,0,0,0/)
 real(dp),allocatable :: work2(:,:,:)

! ***********************************************************************

 if (state == states(imkmem)) then
! We have already read wfrspa(:,:,:,imkmem) last time
  return
 else
! We save the band number for future reference and continue reading
  states(imkmem)=state
 end if

 call int2char4(state,tag)
 fname_wf=trim(dtfil%filnam_ds(5))//'_TDWF'//tag
!DEBUG
!write(message,'(a,a,a,i3)')'reading phi : ',fname_wf, &
!& ' by proc ',me_loc
!call wrtout(6,message,'PERS')
!ENDDEBUG
 open(tmp_unit,file=fname_wf,form='unformatted',status='unknown')
 read(tmp_unit)iband,isppol,eigbnd
!allocate(work2(ndiel4,ndiel5,ndiel6))  !Temporary array for Intel compiler compatibility
!read(tmp_unit) work2
!wfrspa(:,:,:,imkmem)=work2(:,:,:)
!deallocate(work2)
 read(tmp_unit)wfrspa_extract
 close(tmp_unit)

end subroutine read_wfrspa
!!***
