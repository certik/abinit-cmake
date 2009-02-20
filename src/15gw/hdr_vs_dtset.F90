!{\src2tex{textfont=tt}}
!!****f* ABINIT/hdr_vs_dtset
!! NAME
!! hdr_vs_dtset
!!
!! FUNCTION
!!  Check compatibility of the Abinit header with respect to the
!!  input variables defined in the input file. Mainly used in case 
!!  of GW calculation
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Dtset<type(dataset_type)>=all input variables for this dataset
!!  Hdr <type(hdr_type)>=the header structured variable
!!
!! OUTPUT
!!  Only check
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

subroutine hdr_vs_dtset(Hdr,Dtset)
    
 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => hdr_vs_dtset
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(Hdr_type) :: Hdr
 type(Dataset_type),intent(in) :: Dtset

!Local variables-------------------------------
 integer :: ii,jj,nat,nsym,isym
 logical :: test
 logical :: tsymrel,ttnons,tsymafm
 character(len=500) :: msg      
! *************************************************************************

 test=(Hdr%natom==Dtset%natom) ; nat=Dtset%natom
 call assert(test,'Mismatch in natom',__FILE__,__LINE__)

 test=(Hdr%ntypat==Dtset%ntypat) 
 call assert(test,'Mismatch in ntypat',__FILE__,__LINE__)

 test=(Hdr%nsppol==Dtset%nsppol)
 call assert(test,'Mismatch in nsppol',__FILE__,__LINE__)

 test=(Hdr%nspden==Dtset%nspden)
 call assert(test,'Mismatch in nspden',__FILE__,__LINE__)

 test=(Hdr%nspinor==Dtset%nspinor)
 call assert(test,'Mismatch in nspinor',__FILE__,__LINE__)

 test=(Hdr%usepaw==Dtset%usepaw)
 call assert(test,'Mismatch in usepaw',__FILE__,__LINE__)

 test=ALL(ABS(Hdr%xred-Dtset%xred_orig(:,1:nat))<tol6)
 call assert(test,'Mismatch in xred',__FILE__,__LINE__)

 test=ALL(Hdr%typat==Dtset%typat(1:nat)) 
 call assert(test,'Mismatch in typat',__FILE__,__LINE__)
 !
 ! * Check if the lattice from the input file agrees with that read from the KSS file
 if ( (ANY(ABS(Hdr%rprimd-Dtset%rprimd_orig)>tol6)) ) then
  write(msg,'(6a)')ch10,&
&  ' hdr_vs_dtset : ERROR - ',ch10,&
&  ' real lattice vectors read from Header ',ch10,&
&  ' differ from the values specified in the input file'
  call wrtout(std_out,msg,'COLL')
  write(msg,'(3a,3(3es16.6),3a,3(3es16.6),3a)')ch10,&
&  ' rprimd from Hdr file   = ',ch10,(Hdr%rprimd(:,jj),jj=1,3),ch10,&
&  ' rprimd from input file = ',ch10,(Dtset%rprimd_orig(:,jj),jj=1,3),ch10,ch10,&
&  ' Please, modify the lattice vectors in the input file '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 

 test=(Hdr%nsym==Dtset%nsym) ; nsym=Hdr%nsym
 call assert(test,'Mismatch in nsym',__FILE__,__LINE__)

 tsymrel=.TRUE.
 if (ANY(Hdr%symrel/=Dtset%symrel(:,:,1:nsym))) then
  write(msg,'(6a)')ch10,&
&  ' hdr_vs_dtset : ERROR - ',ch10,&
&  ' real space symmetries read from Header ',ch10,&
&  ' differ from the values inferred from the input file'
  call wrtout(std_out,msg,'COLL')
  tsymrel=.FALSE.
 end if 

 ttnons=.TRUE.
 if ( ANY(ABS(Hdr%tnons-Dtset%tnons(:,1:nsym))>tol6) ) then
  write(msg,'(6a)')ch10,&
&  ' hdr_vs_dtset : ERROR - ',ch10,&
&  ' fractional translations read from Header ',ch10,&
&  ' differ from the values inferred from the input file'
  call wrtout(std_out,msg,'COLL')
  ttnons=.FALSE.
 end if 

 tsymafm=ALL(Hdr%symafm==Dtset%symafm(1:nsym))
 if (.not.tsymafm) then
  write(msg,'(6a)')ch10,&
&  ' hdr_vs_dtset : ERROR - ',ch10,&
&  ' AFM symmetries read from Header ',ch10,&
&  ' differ from the values inferred from the input file'
  call wrtout(std_out,msg,'COLL')
  tsymafm=.FALSE.
 end if

 if (.not.(tsymrel.and.ttnons.and.tsymafm)) then
  write(msg,'(a)')' Header ' 
  call wrtout(std_out,msg,'COLL') 
  call print_symmetries(nsym,Hdr%symrel,Hdr%tnons,Hdr%symafm,std_out,'COLL')
  write(msg,'(a)')' Dtset  ' 
  call wrtout(std_out,msg,'COLL') 
  call print_symmetries(nsym,Dtset%symrel,Dtset%tnons,Dtset%symafm,std_out,'COLL')
  call leave_new('COLL')
 end if

end subroutine hdr_vs_dtset
!!***
