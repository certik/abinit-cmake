!{\src2tex{textfont=tt}}
!!****f* ABINIT/print_paw_ij
!! NAME
!! print_paw_ij
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG, MT)
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

subroutine print_paw_ij(Paw_ij,pawprtvol)

 use defs_basis
 use defs_datatypes, only : paw_ij_type
 use m_io_tools,     only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: pawprtvol
!arrays
 type(Paw_ij_type),intent(in) :: Paw_ij(:)

!Local variables-------------------------------
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
!scalars
 integer :: cplex_dij,iatom,idij,lmn2_size,lmn_size,natom,nspden,nsploop,nsppol
 integer :: opt_sym,tmp_cplex_dij
 character(len=500) :: msg
!arrays
 integer,allocatable :: idum(:)
 real(dp),pointer :: dijxc(:),dijxc_val(:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' print_paw_ij : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush(std_out)
#endif

 natom  = SIZE(Paw_ij)
 nsppol = Paw_ij(1)%nsppol
 nspden = Paw_ij(1)%nspden
 nsploop= nsppol ; if (Paw_ij(1)%ndij==4) nsploop=4

 do iatom=1,natom

  lmn_size  = Paw_ij(iatom)%lmn_size
  lmn2_size = Paw_ij(iatom)%lmn2_size
  cplex_dij = Paw_ij(iatom)%cplex_dij

  ! ====================================
  ! === Loop over density components ===
  ! ====================================
  do idij=1,nsploop

   ! Print title
   if (ABS(pawprtvol)>=1) then
    if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
     if (nspden==2.and.nsppol==1) then
      write(msg,'(2a,i3,3a)')ch10,&
&      ' >>>>>>>>>> Atom ',iatom,':',ch10,&
&      ' (antiferromagnetism case: only one spin component)'
     else
      write(msg,'(2a,i3,3a)') ch10,&
&      ' >>>>>>>>>> Atom ',iatom,' (component ',TRIM(dspin(idij+2*(nsploop/4))),'):'
     end if
     call wrtout(std_out,msg,'COLL')
    end if
   end if

   !TODO Clean a bit add SO, LDA+U
   ! check if pointers are associated... now I dont have time

   if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
    if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
     write(msg,'(a)')'   ****************** Dij_xc + Dijhat_xc****************'
     call wrtout(std_out,msg,'COLL')
     if (idij<=nsppol.or.idij==2) then
      opt_sym=2 ; tmp_cplex_dij=1
      dijxc => Paw_ij(iatom)%dijxc(1:cplex_dij*lmn2_size:cplex_dij,idij)
     else
      opt_sym=1 ; tmp_cplex_dij=cplex_dij
      dijxc => Paw_ij(iatom)%dijxc(1:cplex_dij*lmn2_size:1,idij)
     end if
     call print_ij(dijxc,lmn2_size,tmp_cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=opt_sym)
    end if
   end if

   if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
    if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
     write(msg,'(a)')'   ****************** Dij_xcval ****************'
     call wrtout(6,msg,'COLL')
     if (idij<=nsppol.or.idij==2) then
      opt_sym=2 ; tmp_cplex_dij=1
      dijxc_val => paw_ij(iatom)%dijxc_val(1:cplex_dij*lmn2_size:cplex_dij,idij)
     else
      opt_sym=1 ; tmp_cplex_dij=cplex_dij
      dijxc_val => paw_ij(iatom)%dijxc(1:cplex_dij*lmn2_size:1,idij)
     end if
     call print_ij(dijxc_val,lmn2_size,tmp_cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=opt_sym)
    end if
   end if

  end do !idij
 end do !iat

 write(msg,'(a)')' '
 call wrtout(std_out,msg,'COLL')

#if defined DEBUG_MODE
 write(msg,'(a)')' print_paw_ij : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush(std_out)
#endif

end subroutine print_paw_ij
!!***
