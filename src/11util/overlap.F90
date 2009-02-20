!{\src2tex{textfont=tt}}
!!****f* ABINIT/overlap_cmplx
!! NAME
!! overlap_cmplx
!!
!! FUNCTION
!!  Evaluate scalar product between two wavefunctions. It works both for norm-conserving 
!!  pseudopotentials and PAW. The routine assume wavefunctions are given in reciprocal space
!!  but it works also in case of real representation (volume factor not included)
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

function overlap_cmplx(wf1,wf2,usepaw,cprj1,cprj2,typat,pawtab,fact) result(cres)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: usepaw
 real(dp),intent(in),optional :: fact
 complex(dpc) :: cres
!arrays
 integer,intent(in) :: typat(:)
 complex(gwpc),intent(in) :: wf1(:),wf2(:)
 type(cprj_type),intent(in) :: cprj1(:),cprj2(:)
 type(pawtab_type),intent(in) :: pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: iatom,ilmn,itypat,j0lmn,jlmn,klmn,natom
 real(dp) :: sij
 character(len=500) :: msg
!arrays
 real(dp) :: onsite(2)

! *************************************************************************

 if (SIZE(wf1)/=SIZE(wf2)) then 
  write(msg,'(4a,2i6)')ch10,&
&  ' overlap_cmplx : ERROR - ',ch10,&
&  ' wrong size in wf1, wf2 : ',SIZE(wf1),SIZE(wf2)
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if
 cres=dot_product(wf1,wf2) 
 if (PRESENT(fact)) cres=cres*fact

 if (usepaw==1) then 
  natom=SIZE(typat)
  if (SIZE(cprj1)/=SIZE(cprj2) .or. SIZE(cprj1)/=natom) then
   write(msg,'(4a,3i4)')ch10,&
&   ' overlap_cmplx : ERROR - ',ch10,&
&   ' Wrong size in typat, cprj1, cprj2 : ',natom,SIZE(cprj1),SIZE(cprj2)
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if
  onsite(:)=zero
  do iatom=1,natom 
   itypat=typat(iatom)
   do jlmn=1,pawtab(itypat)%lmn_size
    j0lmn=jlmn*(jlmn-1)/2 
    do ilmn=1,jlmn
     klmn=j0lmn+ilmn 
     sij=pawtab(itypat)%sij(klmn)
     if (ABS(sij)<tol16) CYCLE
     if (jlmn==ilmn) sij=sij*half
     onsite(1)=onsite(1) &
&     + sij*( cprj1(iatom)%cp(1,ilmn)*cprj2(iatom)%cp(1,jlmn)&
&     +cprj1(iatom)%cp(2,ilmn)*cprj2(iatom)%cp(2,jlmn)&
&     +cprj1(iatom)%cp(1,jlmn)*cprj2(iatom)%cp(1,ilmn)&
&     +cprj1(iatom)%cp(2,jlmn)*cprj2(iatom)%cp(2,ilmn)&
&     )
     onsite(2)=onsite(2) &
&     + sij*( cprj1(iatom)%cp(1,ilmn)*cprj2(iatom)%cp(2,jlmn)&
&     -cprj1(iatom)%cp(2,ilmn)*cprj2(iatom)%cp(1,jlmn)&
&     +cprj1(iatom)%cp(1,jlmn)*cprj2(iatom)%cp(2,ilmn)&
&     -cprj1(iatom)%cp(2,jlmn)*cprj2(iatom)%cp(1,ilmn)&
&     )
    end do
   end do
  end do
  cres=cres+CMPLX(onsite(1),onsite(2))
 end if

end function overlap_cmplx
!!***

!!****f* ABINIT/overlap
!! NAME
!! overlap_real
!!
!! FUNCTION
!! Evaluate scalar product between two wavefunctions. It works both for norm-conserving 
!! pseudopotentials and PAW. The routine assume wavefunctions are given in reciprocal space
!! but it works also in case of real representation (volume factor not included)
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

function overlap_real(wf1,wf2,usepaw,cprj1,cprj2,typat,pawtab,fact) result(res)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: usepaw
 real(dp),intent(in),optional :: fact
!arrays
 integer,intent(in) :: typat(:)
 real(dp) :: res(2)
 real(dp),intent(in) :: wf1(:,:),wf2(:,:)
 type(cprj_type),intent(in) :: cprj1(:),cprj2(:)
 type(pawtab_type),intent(in) :: pawtab(:)

!Local variables-------------------------------
!scalars
 integer :: iatom,ilmn,itypat,j0lmn,jlmn,klmn,natom
 real(dp) :: sij
 character(len=500) :: msg
!arrays
 real(dp) :: onsite(2)

! *************************************************************************

 if (SIZE(wf1,DIM=2)/=SIZE(wf2,DIM=2)) then
  write(msg,'(4a,2i6)')ch10,&
&  ' overlap_real : ERROR - ',ch10,&
&  ' wrong size in wf1, wf2 : ',SIZE(wf1,DIM=2),SIZE(wf2,DIM=2)
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 res(1)=dot_product(wf1(1,:),wf2(1,:)) + dot_product(wf1(2,:),wf2(2,:))
 res(2)=dot_product(wf1(1,:),wf2(2,:)) - dot_product(wf1(2,:),wf2(1,:))
 if (PRESENT(fact)) res=res*fact

 if (usepaw==1) then 
  natom=SIZE(typat)
  if (SIZE(cprj1)/=SIZE(cprj2).or.natom/=SIZE(cprj1)) then
   write(msg,'(4a,3i4)')ch10,&
&   ' overlap_real : ERROR - ',ch10,&
&   ' wrong size in typat,cprj1,cprj2 : ',natom,SIZE(cprj1),SIZE(cprj2)
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
  end if
  onsite(:)=zero
  do iatom=1,natom 
   itypat=typat(iatom)
   do jlmn=1,pawtab(itypat)%lmn_size
    j0lmn=jlmn*(jlmn-1)/2 
    do ilmn=1,jlmn
     klmn=j0lmn+ilmn 
     sij=pawtab(itypat)%sij(klmn)
     if (ABS(sij)<tol16) CYCLE
     if (jlmn==ilmn) sij=sij*half
     onsite(1)=onsite(1) &
&     + sij*( cprj1(iatom)%cp(1,ilmn)*cprj2(iatom)%cp(1,jlmn)&
&     +cprj1(iatom)%cp(2,ilmn)*cprj2(iatom)%cp(2,jlmn)&
&     +cprj1(iatom)%cp(1,jlmn)*cprj2(iatom)%cp(1,ilmn)&
&     +cprj1(iatom)%cp(2,jlmn)*cprj2(iatom)%cp(2,ilmn)&
&     )
     onsite(2)=onsite(2) &
&     + sij*( cprj1(iatom)%cp(1,ilmn)*cprj2(iatom)%cp(2,jlmn)&
&     -cprj1(iatom)%cp(2,ilmn)*cprj2(iatom)%cp(1,jlmn)&
&     +cprj1(iatom)%cp(1,jlmn)*cprj2(iatom)%cp(2,ilmn)&
&     -cprj1(iatom)%cp(2,jlmn)*cprj2(iatom)%cp(1,ilmn)&
&     )
    end do
   end do
  end do
  res(1)=res(1)+onsite(1)
  res(2)=res(2)+onsite(2)
 end if

end function overlap_real
!!***
