!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawr
!! NAME
!! pawr
!!
!! FUNCTION
!! Evaluated matrix elements of the position operator between PAW AE partial waves for a given atom.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  iat=index of the atom to be considered
!!  Pawtab(ntypat) <type(pawtab_type)>=paw tabulated data read at start:
!!     %l_size
!!     %lmn_size
!!     %lmn2_size
!!     %indklmn
!!     %phiphj 
!!  Pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!  Pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Maximum value of angular momentum l+1
!!     %gntselect((2*l_max-1)**2,l_max**2,l_max**2)= selection rules for Gaunt coefficients
!!     %realgnt
!!  Psps <type(pseudopotential_type)>=Information on pseudopotentials
!!     %indlmn
!!  natom=number of atoms in unit cell
!!  ntypat=number of types of atom
!!  typat(natom)=type of each atom
!!  xcart(3,natom)=cartesian coordinates
!!
!! OUTPUT
!!  rc_onsite(3,lmn2_size)
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

subroutine pawr(iat,Pawtab,Pawrad,Pawang,Psps,natom,ntypat,typat,xcart,lmn2_size,rc_onsite)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_IO_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iat,lmn2_size,natom,ntypat
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(inout) :: rc_onsite(3,lmn2_size)
 type(Pawrad_type),intent(in) :: Pawrad(ntypat)
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: idir,ignt,il,ilm,ilm_G,ilmn,iln,im,itypat,jl,jlm,jlmn,jln,jm,k0lm
 integer :: k0lmn,k0ln,klm,klmn,kln,l_size,lmn_size,mesh_size,mm_G
 real(dp) :: fact,intff,rgnt
 logical :: ltest
 character(len=500) :: msg
!arrays
 real(dp) :: r0(3)
 real(dp),allocatable :: ff(:),indklmn_(:,:),rad(:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' pawr : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 r0(:)=xcart(:,iat)
 itypat=typat(iat)
 l_size    =Pawtab(itypat)%l_size
 lmn_size  =Pawtab(itypat)%lmn_size
 mesh_size =Pawrad(itypat)%mesh_size

 ltest=(lmn2_size==Pawtab(itypat)%lmn2_size)
 call assert(ltest,'Mismatch in lmn2_size',__FILE__,__LINE__)

 allocate(indklmn_(4,lmn2_size)) 
 indklmn_(:,:)=Pawtab(itypat)%indklmn(:,:)

 allocate(ff(mesh_size),rad(mesh_size))
 rad(1:mesh_size)=Pawrad(itypat)%rad(1:mesh_size)

 fact=two*SQRT(pi/three) ! FIXME check factor
 rc_onsite(:,:)=zero
 !
 ! === Loop on (jl,jm,jn) channels 
 do jlmn=1,lmn_size
  jl =Psps%indlmn(1,jlmn,itypat)
  jm =Psps%indlmn(2,jlmn,itypat)
  jlm=Psps%indlmn(4,jlmn,itypat)
  jln=Psps%indlmn(5,jlmn,itypat)
 
  k0lmn=jlmn*(jlmn-1)/2 
  k0lm =jlm *(jlm -1)/2
  k0ln =jln *(jln -1)/2
  !
  ! === Loop on (il,im,in) channels; klmn is the index for packed form ===
  do ilmn=1,jlmn 
   il =Psps%indlmn(1,ilmn,itypat)
   im =Psps%indlmn(2,ilmn,itypat)
   ilm=Psps%indlmn(4,ilmn,itypat)
   iln=Psps%indlmn(5,ilmn,itypat)
 
   klmn=k0lmn+ilmn 
   klm =k0lm +ilm
   kln =k0ln +iln
   !
   ! === For each cartesian direction, use expansion in terms of RSH ===
   do idir=1,3
    if (idir==1) mm_G= 1
    if (idir==2) mm_G=-1
    if (idir==3) mm_G= 0
    !ilm_G=1+ll_G**2+ll_G+mm_G
    ilm_G=1+1+1+mm_G
    ignt=Pawang%gntselect(ilm_G,klm)

    if (ignt/=0) then
     rgnt=Pawang%realgnt(ignt)
     ff(1)=zero
     !ff(2:mesh_size)=(Pawtab(itypat)%phiphj(2:mesh_size,kln)-Pawtab(itypat)%tphitphj(2:mesh_size,kln))*rad(2:mesh_size)
     ff(2:mesh_size)=Pawtab(itypat)%phiphj(2:mesh_size,kln)*rad(2:mesh_size)
     call simp_gen(intff,ff,Pawrad(itypat))
     rc_onsite(idir,klmn)=fact*intff*rgnt
    end if

    if (jl==il.and.jm==im) then 
     ff(1:mesh_size)=Pawtab(itypat)%phiphj(1:mesh_size,kln)
     !ff(1:mesh_size)=(Pawtab(itypat)%phiphj(1:mesh_size,kln)-Pawtab(itypat)%tphitphj(1:mesh_size,kln))
     call simp_gen(intff,ff,Pawrad(itypat))
     rc_onsite(idir,klmn)=rc_onsite(idir,klmn)+r0(idir)*intff
    end if
   end do !idir

  end do !ilmn
 end do !jllmn

 deallocate(ff,rad,indklmn_)

#if defined DEBUG_MODE
 write(msg,'(a)')' pawr : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine pawr
!!***
