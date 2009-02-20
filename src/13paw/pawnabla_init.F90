!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawnabla_init
!! NAME
!! pawnabla_init
!!
!! FUNCTION
!! Evaluate onsite contributions of the nabla operator in caartesian coordinates.
!! Store values in the pawtab% data structure
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (SM,VR,FJ,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  indlmn(6,lmnmax,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  lmnmax=max. number of (l,m,n) numbers over all types of atom
!!  mpsang=1+maximum angular momentum
!!  ntypat=number of types of atoms in cell
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data:
!!     %mesh_size=Dimension of radial mesh
!!     %rad(mesh_size)=The coordinates of all the points of the radial mesh
!!  pawtab(ntypat) <type(pawtab_type>=PAW tabulated starting data:
!!    %lmn_size=Number of (l,m,n) elements for the PAW basis
!!
!! OUTPUT
!!  See side effects
!!
!! SIDE EFFECTS
!!  pawtab(ntypat) <type(pawtab_type>=PAW tabulated starting data:
!!    %has_nabla=set to 1 in matrix elements are calculated and stored
!!    %nabla_ij(3,lmn_size,lmn_size)= <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j>
!!
!! NOTES
!!  MG extracted this piece of code from optics_paw.F90 in order to have something more 
!!  reusable! Note however the storage mode of nabla_ij differs from optics_paw 
!!  (here Cartesian coordinates runs faster). Besides nabla_ij contains the matrix 
!!  elements of \nabla instead of the elements of the momentum operator p.
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

subroutine pawnabla_init(mpsang,lmnmax,ntypat,indlmn,pawrad,pawtab)
    
 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13paw, except_this_one => pawnabla_init
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: lmnmax,mpsang,ntypat
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat)
 type(pawtab_type),intent(inout) :: pawtab(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)

!Local variables-------------------------------
!scalars
 integer :: nln,iatom,il,ilm,ilmn,iln,itypat
 integer :: jl,jlm,jlmn,jln,lmn_size,mesh_size 
 real(dp) :: intg
 character(len=500) :: msg      
!arrays
 integer,allocatable :: indlmn_(:,:)
 real(dp) :: ang_phipphj(mpsang**2,mpsang**2,8)
 real(dp),allocatable :: dphi(:),dtphi(:),ff(:),int1(:,:),int2(:,:),rad(:)
 
! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' pawnabla_init : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 if (mpsang>4)then
  write(msg,'(6a)')ch10,&
&  ' pawinit : ERROR -',ch10,&
&  '  Not designed for angular momentum greater than 3 ',ch10,&
&  '  Modification in the table defined in ang_int.F90 is required '
  call wrtout(std_out,msg,'COLL') 
  call leave_new('COLL')
 end if
!
!Integration of the angular part: all angular integrals have been computed 
!outside Abinit and tabulated for each (l,m) value up to l=2
 call int_ang(ang_phipphj,mpsang)

 do itypat=1,ntypat
! 
! ==================================================
! COMPUTE nabla_ij := <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for this type
! 
  mesh_size=pawrad(itypat)%mesh_size
  lmn_size=pawtab(itypat)%lmn_size
  nln=pawtab(itypat)%basis_size

! if (associated(pawtab(itypat)%nabla_ij)) deallocate(pawtab(itypat)%nabla_ij)
  allocate(pawtab(itypat)%nabla_ij(3,lmn_size,lmn_size))
  pawtab(itypat)%has_nabla=1

  allocate(ff(mesh_size),rad(mesh_size))
  allocate(int2(lmn_size,lmn_size),int1(lmn_size,lmn_size))
  allocate(dphi(mesh_size),dtphi(mesh_size))
  allocate(indlmn_(6,lmnmax))

  indlmn_(:,:)=indlmn(:,:,itypat)
  rad(1:mesh_size)=pawrad(itypat)%rad(1:mesh_size)
! 
! int1=\int phi phj/r dr - \int tphi tphj /r dr
  do jln=1,nln
   do iln=1,nln
    ff(2:mesh_size)=( pawtab(itypat)%phi (2:mesh_size,iln)*pawtab(itypat)%phi (2:mesh_size,jln) &
&    -pawtab(itypat)%tphi(2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln) )/rad(2:mesh_size)
    call deducer0(ff,mesh_size,pawrad(itypat))
    call simp_gen(intg,ff,pawrad(itypat))
    int1(iln,jln)=intg
   end do
  end do
! 
! int2=\int phi/r d/dr(phj/r) r^2dr - \int tphi/r d/dr(tphj/r)r^2 dr
  do jln=1,nln
   ff(1:mesh_size)=pawtab(itypat)%phi(1:mesh_size,jln)
   call nderiv_gen(dphi,ff,1,pawrad(itypat))
   ff(1:mesh_size)=pawtab(itypat)%tphi(1:mesh_size,jln)
   call nderiv_gen(dtphi,ff,1,pawrad(itypat))

   do iln=1,nln
    ff(2:mesh_size)= &
&    pawtab(itypat)%phi (2:mesh_size,iln)*dphi (2:mesh_size) &
&    -pawtab(itypat)%phi (2:mesh_size,iln)*pawtab(itypat)%phi (2:mesh_size,jln)/rad(2:mesh_size) &
&    -( pawtab(itypat)%tphi(2:mesh_size,iln)*dtphi(2:mesh_size) &
&    -pawtab(itypat)%tphi(2:mesh_size,iln)*pawtab(itypat)%tphi(2:mesh_size,jln)/rad(2:mesh_size) )
    call deducer0(ff,mesh_size,pawrad(itypat))
    call simp_gen(intg,ff,pawrad(itypat))
    int2(iln,jln)=intg
   end do
  end do
! 
! 1-c Integration of the radial part, Note unpacked loop
  do jlmn=1,lmn_size
   jlm=indlmn_(4,jlmn)
   jl =indlmn_(5,jlmn)
   do ilmn=1,lmn_size
    ilm=indlmn_(4,ilmn)
    il =indlmn_(5,ilmn)

    pawtab(itypat)%nabla_ij(1,ilmn,jlmn)= &
&    int2(il,jl)* ang_phipphj(ilm,jlm,1) &
&    +int1(il,jl)*(ang_phipphj(ilm,jlm,2)+ang_phipphj(ilm,jlm,3))

    pawtab(itypat)%nabla_ij(2,ilmn,jlmn)= &
&    int2(il,jl)* ang_phipphj(ilm,jlm,4) &
&    +int1(il,jl)*(ang_phipphj(ilm,jlm,5)+ang_phipphj(ilm,jlm,6))

    pawtab(itypat)%nabla_ij(3,ilmn,jlmn)= &
&    int2(il,jl)* ang_phipphj(ilm,jlm,7) &
&    +int1(il,jl)* ang_phipphj(ilm,jlm,8)

   end do !ilmn
  end do !jlmn

  deallocate(ff,rad)
  deallocate(int2,int1)
  deallocate(dphi,dtphi)
  deallocate(indlmn_)

 end do !itypat

#if defined DEBUG_MODE
 write(msg,'(a)')' pawnabla_init : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine pawnabla_init
!!***

