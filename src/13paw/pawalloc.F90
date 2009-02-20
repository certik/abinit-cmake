!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawalloc
!! NAME
!! pawalloc
!!
!! FUNCTION
!! Allocate or deallocate datastructures used for PAW
!! at the level of the driving routine (driver.F90)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  idtset=index of the current dataset
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mqgrid_vl=dimension of q (or G) grid for Vloc (array vlspl)
!!  npsp=number of pseudopotentials
!!  option=1: allocate PAW datastructures for a new dataset
!!         2: deallocate PAW datastructures depending on paw_size (pawrad,pawtab) for the current dataset
!!         3: deallocate all PAW datastructures (pawang,pawrad,pawtab) for the current dataset
!!  pspheads(npsp)=<type pspheader_type>all the important information from the
!!                 pseudopotential file header, as well as the psp file name
!!
!! SIDE EFFECTS
!! Allocated/deallocated:
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(paw_size) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(paw_size) <type(pawtab_type)>=paw tabulated starting data
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawalloc(dtset,idtset,mpsang,mqgrid_vl,npsp,option,paw_size,paw_size_old,&
&                   pawang,pawrad,pawtab,pspheads)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idtset,mpsang,mqgrid_vl,npsp,option,paw_size,paw_size_old
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(inout) :: pawang
!arrays
 type(pawrad_type),intent(inout) :: pawrad(paw_size)
 type(pawtab_type),intent(inout) :: pawtab(paw_size)
 type(pspheader_type),intent(in) :: pspheads(npsp)

!Local variables-------------------------------
!scalars
 integer,parameter :: mqgrid_shp=300
 integer :: angl_size_new,basis_size_new,ierr,ij_size_new,itypat,gnt_option_new,l_max_new
 integer :: l_size_new,l_size_max_new,lexexch_new,lmn_size_new,lmn2_size_new,lpawu_new
 integer :: mesh_size_new,mqgrid_shp_new,nsym_new,pawspnorb_new
 logical :: test_alloc,test_angl_size,test_basis_size,test_ij_size,test_gnt_option,test_l_max
 logical :: test_l_size,test_l_size_max,test_lexexch,test_lmn2_size,test_lpawu,test_mesh_size
 logical :: test_mqgrid,test_mqgrid_shp,test_new,test_nsym,test_spnorb
!arrays

! *********************************************************************

!DEBUG
!write(6,*)' pawalloc : enter, option=',option
!ENDDEBUG

!Nothing to do if not PAW
 if (paw_size==0) return

!=====================================================================
!========================= ALLOCATIONS ===============================
!=====================================================================
 if (option==1) then

  test_new=(paw_size/=paw_size_old)
  test_alloc=((idtset/=1).and.(paw_size==paw_size_old))

  l_max_new=mpsang
  l_size_max_new=2*l_max_new-1
  gnt_option_new=1;if (dtset%pawxcdev==2) gnt_option_new=2
  angl_size_new=0;if(dtset%pawxcdev==0) angl_size_new=dtset%pawntheta*dtset%pawnphi
  nsym_new=dtset%nsym
  pawspnorb_new=0;if (dtset%pawspnorb>0) pawspnorb_new=1
  mqgrid_shp_new=0
  if (dtset%optdriver==0.and. &
&  ((dtset%iprcel>=20.and.dtset%iprcel<70).or.dtset%iprcel>=80)) mqgrid_shp_new=mqgrid_shp

  if (idtset==1) then
   pawang%angl_size=-1
   pawang%l_max=-1
   pawang%l_size_max=-1
   pawang%nsym=-1
   pawang%gnt_option=-1
   pawang%use_ls_ylm=-1
  end if

  test_l_max=(l_max_new/=pawang%l_max)
  test_angl_size=(angl_size_new/=pawang%angl_size)
  test_l_size_max=(l_size_max_new/=pawang%l_size_max)
  test_gnt_option=(gnt_option_new/=pawang%gnt_option)
  test_spnorb=(pawspnorb_new/=pawang%use_ls_ylm)
  test_nsym=(nsym_new/=pawang%nsym)

  do itypat=1,dtset%ntypat

   basis_size_new=pspheads(itypat)%pawheader%basis_size
   lmn_size_new  =pspheads(itypat)%pawheader%lmn_size
   l_size_new    =pspheads(itypat)%pawheader%l_size
   mesh_size_new =pspheads(itypat)%pawheader%mesh_size
   lmn2_size_new =lmn_size_new*(lmn_size_new+1)/2
   ij_size_new   =basis_size_new*(basis_size_new+1)/2
   lpawu_new     =dtset%lpawu(itypat)
   lexexch_new   =dtset%lexexch(itypat)

   if (idtset==1) then
    test_basis_size=.true.;test_ij_size  =.true.
    test_mesh_size =.true.;test_l_size   =.true.
    test_mqgrid_shp=.true.;test_lmn2_size=.true.
    test_lpawu     =.true.;test_lexexch  =.true.
    test_mqgrid    =.true.
   else
    test_basis_size=(basis_size_new/=pawtab(itypat)%basis_size)
    test_ij_size=(ij_size_new/=pawtab(itypat)%ij_size)
    test_mesh_size=(mesh_size_new/=pawtab(itypat)%mesh_size)
    test_l_size=(l_size_new/=pawtab(itypat)%l_size)
    test_mqgrid_shp=(mqgrid_shp_new/=pawtab(itypat)%mqgrid_shp)
    test_lmn2_size=(lmn2_size_new/=pawtab(itypat)%lmn2_size)
    test_lpawu=(lpawu_new/=pawtab(itypat)%lpawu)
    test_lexexch=(lexexch_new/=pawtab(itypat)%lexexch)
    test_mqgrid=(mqgrid_vl/=pawtab(itypat)%mqgrid)
   end if

!  Reallocate arrays depending on mesh_size and basis_size
   if (test_new.or.test_ij_size.or.test_basis_size) then
    if(test_alloc) deallocate(pawtab(itypat)%tphi,pawtab(itypat)%phi)
    allocate(pawtab(itypat)%tphi(mesh_size_new,basis_size_new))
    allocate(pawtab(itypat)%phi (mesh_size_new,basis_size_new))
   end if

!  Reallocate arrays depending on mesh_size and ij_size
   if (test_new.or.test_mesh_size.or.test_ij_size) then
    if(test_alloc) deallocate(pawtab(itypat)%tphitphj,pawtab(itypat)%phiphj)
    allocate(pawtab(itypat)%tphitphj(mesh_size_new,ij_size_new))
    allocate(pawtab(itypat)%phiphj  (mesh_size_new,ij_size_new))
   end if

!  Reallocate arrays depending on mesh_size and l_size
   if (test_new.or.test_mesh_size.or.test_l_size) then
    if(test_alloc) deallocate(pawtab(itypat)%shapefunc,pawtab(itypat)%dshpfunc)
    allocate(pawtab(itypat)%shapefunc(mesh_size_new,l_size_new),&
&    pawtab(itypat)%dshpfunc(mesh_size_new,l_size_new,0)) ! Will be allocated later (only if shape_type=-1)
   end if

!  Reallocate arrays depending on mqgrid_shp and l_size
   if (test_new.or.test_mqgrid_shp.or.test_l_size) then
    if(test_alloc.and.pawtab(itypat)%mqgrid_shp>0) deallocate(pawtab(itypat)%shapefncg)
    if (mqgrid_shp_new>0) allocate(pawtab(itypat)%shapefncg(mqgrid_shp_new,2,l_size_new))
   end if

!  Reallocate arrays depending on l_size and lmn2_size
   if (test_new.or.test_l_size.or.test_lmn2_size) then
    if(test_alloc) deallocate(pawtab(itypat)%qijl)
    allocate(pawtab(itypat)%qijl(l_size_new*l_size_new,lmn2_size_new))
   end if

!  Reallocate arrays depending on lpawu, lexexch and lmn2_size
   if (test_new.or.test_lmn2_size.or.test_lpawu.or.test_lexexch) then
    if(test_alloc.and.(pawtab(itypat)%lpawu>=0.or.pawtab(itypat)%lexexch>=0)) deallocate(pawtab(itypat)%klmntomn)
    if (lpawu_new>=0.or.lexexch_new>=0) allocate(pawtab(itypat)%klmntomn(4,lmn2_size_new))
   end if

!  Reallocate arrays depending on lpawu and lmn2_size
   if (test_new.or.test_lpawu) then
    if(test_alloc.and.(pawtab(itypat)%lpawu>=0)) deallocate(pawtab(itypat)%vee)
    if (lpawu_new>=0) allocate(pawtab(itypat)%vee(2*dtset%lpawu(itypat)+1,&
&    2*dtset%lpawu(itypat)+1,2*dtset%lpawu(itypat)+1,2*dtset%lpawu(itypat)+1))
   end if

!  Reallocate arrays depending on lexexch and lmn2_size
   if (test_new.or.test_lexexch) then
    if(test_alloc.and.(pawtab(itypat)%lexexch>=0)) deallocate(pawtab(itypat)%vex,pawtab(itypat)%fk)
    if (lexexch_new>=0) then
     allocate(pawtab(itypat)%vex(2*dtset%lexexch(itypat)+1,2*dtset%lexexch(itypat)+1,&
&     2*dtset%lexexch(itypat)+1,2*dtset%lexexch(itypat)+1,4))
     allocate(pawtab(itypat)%fk(6,4))
    end if
   end if

!  Reallocate arrays depending on l_size
   if (test_new.or.test_l_size) then
    if(test_alloc) deallocate(pawtab(itypat)%gnorm,pawtab(itypat)%shape_alpha,&
&    pawtab(itypat)%shape_q)
    allocate(pawtab(itypat)%gnorm(l_size_new),&
&    pawtab(itypat)%shape_alpha(2,l_size_new),&
&    pawtab(itypat)%shape_q(2,l_size_new))
   end if

!  Reallocate arrays depending on lmn2_size
   if (test_new.or.test_lmn2_size) then
    if(test_alloc) deallocate(pawtab(itypat)%eijkl,pawtab(itypat)%dij0,&
&    pawtab(itypat)%dltij,pawtab(itypat)%rhoij0,&
&    pawtab(itypat)%sij,pawtab(itypat)%indklmn)
    if (idtset==1) nullify(pawtab(itypat)%kmix)
    if (associated(pawtab(itypat)%kmix)) deallocate(pawtab(itypat)%kmix,STAT=ierr)
    allocate(pawtab(itypat)%eijkl(lmn2_size_new,lmn2_size_new))
    allocate(pawtab(itypat)%dij0(lmn2_size_new))
    allocate(pawtab(itypat)%dltij(lmn2_size_new))
    allocate(pawtab(itypat)%rhoij0(lmn2_size_new))
    allocate(pawtab(itypat)%sij(lmn2_size_new))
    allocate(pawtab(itypat)%indklmn(6,lmn2_size_new))
   end if

!  Reallocate arrays depending on mesh_size
   if (test_new.or.test_mesh_size) then
    if(test_alloc) deallocate(pawtab(itypat)%coredens,pawtab(itypat)%tcoredens,&
&    pawrad(itypat)%rad,pawrad(itypat)%radfact,&
&    pawrad(itypat)%simfact,pawtab(itypat)%rad_for_spline)
    allocate(pawtab(itypat)%coredens (mesh_size_new))
    allocate(pawtab(itypat)%tcoredens(mesh_size_new))
    allocate(pawrad(itypat)%rad      (mesh_size_new))
    allocate(pawrad(itypat)%radfact  (mesh_size_new))
    allocate(pawrad(itypat)%simfact  (mesh_size_new))
    allocate(pawtab(itypat)%rad_for_spline(0)) ! Will be allocated later (only if needed)
   end if

!  Reallocate arrays depending on mqgrid_vl
   if (test_new.or.test_mqgrid) then
    if(test_alloc) then
     deallocate(pawtab(itypat)%tcorespl)
     if (pawtab(itypat)%usetvale>0) deallocate(pawtab(itypat)%tvalespl)
    end if
    allocate(pawtab(itypat)%tcorespl(mqgrid_vl,2))
    if (pspheads(itypat)%pawheader%pawver>=4) allocate(pawtab(itypat)%tvalespl(mqgrid_vl,2))
   end if

!  Reallocate arrays depending on mqgrid_shp and l_size
   if (test_new.or.test_mqgrid_shp)  then
    if(test_alloc.and.pawtab(itypat)%mqgrid_shp>0) deallocate(pawtab(itypat)%qgrid_shp)
    if (mqgrid_shp_new>0) allocate(pawtab(itypat)%qgrid_shp(mqgrid_shp_new))
   end if

   pawtab(itypat)%basis_size=basis_size_new
   pawtab(itypat)%ij_size=ij_size_new
   pawtab(itypat)%lpawu=lpawu_new
   pawtab(itypat)%lexexch=lexexch_new
   pawtab(itypat)%l_size=l_size_new
   pawtab(itypat)%lmn_size=lmn_size_new
   pawtab(itypat)%lmn2_size=lmn2_size_new
   pawtab(itypat)%mesh_size=mesh_size_new
   pawtab(itypat)%mqgrid=mqgrid_vl
   pawtab(itypat)%mqgrid_shp=mqgrid_shp_new
   pawtab(itypat)%usetvale=0;if (pspheads(itypat)%pawheader%pawver>=4) pawtab(itypat)%usetvale=1

  end do ! itypat

! Reallocate arrays depending on angl_size and l_size_max
  if (test_angl_size.or.test_l_size_max) then
   if(idtset/=1) deallocate(pawang%ylmr)
   allocate(pawang%ylmr(l_size_max_new**2,angl_size_new))
  end if

! Reallocate arrays depending on nsym, l_size_max and l_max
  if (test_nsym.or.test_l_size_max.or.test_l_max) then
   if(idtset/=1) deallocate(pawang%zarot)
   allocate(pawang%zarot(l_size_max_new,l_size_max_new,l_max_new,nsym_new))
  end if

! Reallocate arrays depending on l_max, l_size_max and gnt_option
  if (test_l_max.or.test_l_size_max.or.test_gnt_option) then
   if(idtset/=1) deallocate(pawang%gntselect)
   if(gnt_option_new/=2) then
    allocate(pawang%gntselect((l_size_max_new)**2,(l_max_new**2)*(l_max_new**2+1)/2))
   else
    allocate(pawang%gntselect((2*l_size_max_new-1)**2,((2*l_max_new-1)**2)*((2*l_max_new-1)**2+1)/2))
   end if
   if (idtset==1) nullify(pawang%realgnt)
   if (associated(pawang%realgnt)) deallocate(pawang%realgnt,STAT=ierr)
  end if

! Reallocate arrays depending on l_max and pawspnorb
  if (test_l_max.or.test_spnorb) then
   if(idtset/=1.and.pawang%use_ls_ylm>0) deallocate(pawang%ls_ylm)
   if (pawspnorb_new>0) allocate(pawang%ls_ylm(2,l_max_new**2*(l_max_new**2+1)/2,2))
  end if

! Reallocate arrays depending on angl_size
  if (test_angl_size) then
   if(idtset/=1) deallocate(pawang%angwgth)
   allocate(pawang%angwgth(angl_size_new))
  end if

  pawang%angl_size=angl_size_new
  pawang%l_max=l_max_new
  pawang%l_size_max=l_size_max_new
  pawang%nsym=nsym_new
  pawang%gnt_option=gnt_option_new
  pawang%use_ls_ylm=pawspnorb_new

! =====================================================================
! ======================== DEALLOCATIONS ==============================
! =====================================================================
 else if (option>=2) then

  do itypat=1,dtset%ntypat

   deallocate(pawrad(itypat)%rad)
   deallocate(pawrad(itypat)%radfact)
   deallocate(pawrad(itypat)%simfact)
   deallocate(pawtab(itypat)%gnorm)
   deallocate(pawtab(itypat)%indklmn)
   if(pawtab(itypat)%lpawu>=0.or.pawtab(itypat)%lexexch>=0) deallocate(pawtab(itypat)%klmntomn)
   if(pawtab(itypat)%lpawu>=0) deallocate(pawtab(itypat)%vee)
   if(pawtab(itypat)%lexexch>=0) deallocate(pawtab(itypat)%vex)
   if(pawtab(itypat)%lexexch>=0) deallocate(pawtab(itypat)%fk)
   deallocate(pawtab(itypat)%shapefunc)
   deallocate(pawtab(itypat)%shape_alpha)
   deallocate(pawtab(itypat)%shape_q)
   deallocate(pawtab(itypat)%tphi)
   deallocate(pawtab(itypat)%phi)
   deallocate(pawtab(itypat)%tphitphj)
   deallocate(pawtab(itypat)%phiphj)
   deallocate(pawtab(itypat)%coredens)
   deallocate(pawtab(itypat)%tcoredens)
   deallocate(pawtab(itypat)%tcorespl)
   if (pawtab(itypat)%usetvale>0) deallocate(pawtab(itypat)%tvalespl)
   deallocate(pawtab(itypat)%qijl)
   deallocate(pawtab(itypat)%eijkl)
   deallocate(pawtab(itypat)%dij0)
   deallocate(pawtab(itypat)%dltij)
   deallocate(pawtab(itypat)%rhoij0)
   deallocate(pawtab(itypat)%sij)
   if (associated(pawtab(itypat)%kmix)) deallocate(pawtab(itypat)%kmix,STAT=ierr)
   if (associated(pawtab(itypat)%dshpfunc)) deallocate(pawtab(itypat)%dshpfunc,STAT=ierr)
   if (associated(pawtab(itypat)%rad_for_spline)) deallocate(pawtab(itypat)%rad_for_spline,STAT=ierr)
   if (pawtab(itypat)%mqgrid_shp>0) then
    pawtab(itypat)%mqgrid_shp=0
    deallocate(pawtab(itypat)%qgrid_shp,pawtab(itypat)%shapefncg)
   end if
  end do

  if (option>=3) then
   deallocate(pawang%angwgth,pawang%ylmr,pawang%zarot)
   if (associated(pawang%gntselect)) deallocate(pawang%gntselect,STAT=ierr)
   if (associated(pawang%realgnt)) deallocate(pawang%realgnt,STAT=ierr)
   if(pawang%use_ls_ylm>0) deallocate(pawang%ls_ylm)
  end if

 end if

!DEBUG
!write(6,*)' pawalloc : exit'
!ENDDEBUG

end subroutine pawalloc

!!***
