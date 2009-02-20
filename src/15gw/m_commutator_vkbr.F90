!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_commutator_vkbr
!! NAME
!!  m_commutator_vkbr
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!@ ABIDOC
MODULE m_commutator_vkbr

 use defs_basis
 implicit none

! === List of available public types and routines  ===
 type  KB_form_factors

!scalars
  integer :: mpw       ! Maximum number of plane waves over k-points
  integer :: nkibz     ! Number of irreducible k-points.
  integer :: lnmax     ! Max. number of (l,n) components over all type of pseudos.
  integer :: ntypat    ! Number of type of atoms.

!arrays
  integer,pointer :: sign_dyad(:,:)  ! sign(lnmax,ntypat). sign of the KB dyadic product.

  real(dp),pointer :: ff (:,:,:,:)   ! ff (npw,lnmax,ntypat,nkibz) KB form factor.
  real(dp),pointer :: ffd(:,:,:,:)   ! ffd(npw,lnmax,ntypat,nkibz) derivative of ff wrt k+G of ff.

 end type
!@END ABIDOC

CONTAINS  !===========================================================

!!*** 

!!****f* m_commutator_vkbr/vkb_init
!! NAME
!!  vkb_init
!!
!! FUNCTION
!!  Calculate KB form factors and derivatives required to evalute
!!  the matrix elements of the commutator [Vnl,r]-. 
!!  This term enters the expression for the oscillator strengths in 
!!  the optical limit q-->0. Pseudopotentials with more than one
!!  projector per angular channel are supported.
!!
!! INPUT 
!!
!! TODO 
!!  Replace old implementation with this new routine. Matrix elements
!!  of the commutator should be calculated on-the-fly in screening only
!!  if really needed. This is the first step toward the elimination
!!  of the KSS file. Modifications in cchi0q0 are needed.
!!
!! PARENTS
!!
!! CHILDREN
!!      assert,metric,mkffnl,mkkin,wrtout
!!
!! SOURCE

subroutine vkb_init(ntypat,Kmesh,Psps,mkmem,mpw,npwarr,ecut,kg,rprimd,Vkb)

 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
 use interfaces_13nonlocal
 use interfaces_13recipspace
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ntypat,mpw,mkmem
 real(dp),intent(in) :: ecut
 type(Pseudopotential_type),intent(in) :: Psps
 type(Bz_mesh_type),intent(in) :: Kmesh
 type(KB_form_factors),intent(inout) :: Vkb
!arrays
 integer,target,intent(in) :: kg(3,mpw*mkmem),npwarr(Kmesh%nibz)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables ------------------------------
!scalars
 integer :: dimffnl,ider,idir,npw_k,ikibz,itypat,istat,nkpg
 integer :: il,ilmn,iln,iln0,nlmn,ikg
 real(dp) :: fact,effmass,ecutsm,ucvol
 character(len=500) :: msg
!arrays
 integer,pointer :: kg_k(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp),pointer :: kpoint(:)
 real(dp),allocatable :: ffnl(:,:,:,:),kpg_dum(:,:),modkplusg(:)
 real(dp),allocatable :: ylm(:,:),ylm_gr(:,:,:),ylm_k(:,:)

! *************************************************************************

 ! == Test the arguments ===
 call assert((Psps%usepaw==0),"You should not be here!",__FILE__,__LINE__)
 call assert((mkmem/=0),"mkmem==0 not implemented",__FILE__,__LINE__)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 !
 ! === Save Vkb dyadic sign (integer-valued) ===
 ! * Notice how the ordering is chosen correctly unlike in outkss.
 ! * More than one projector per angular channel is allowed.
 !  allocate(vkbsign(Psps%mpsang,ntypat))  THIS THE OLD IMPLEMENTATION
 allocate(Vkb%sign_dyad(Psps%lnmax,ntypat)) 
 Vkb%sign_dyad(:,:)=0

 do itypat=1,ntypat
  nlmn=COUNT(Psps%indlmn(3,:,itypat)>0)
  iln0=0 
  do ilmn=1,nlmn
   iln=Psps%indlmn(5,ilmn,itypat)
   if (iln>iln0) then
    iln0=iln
    Vkb%sign_dyad(iln,itypat)=NINT(DSIGN(one,Psps%ekb(ilmn,itypat)))
   end if
  end do
 end do

 Vkb%nkibz  = Kmesh%nibz
 Vkb%lnmax  = Psps%lnmax
 Vkb%ntypat = ntypat
 Vkb%mpw    = mpw

 ! === Allocate KB form factor and derivative wrt k+G ===
 ! * Also here we use correct ordering for dimensions
 allocate(Vkb%ff (mpw,Vkb%lnmax,ntypat,Vkb%nkibz), stat=istat)
 allocate(Vkb%ffd(mpw,Vkb%lnmax,ntypat,Vkb%nkibz), stat=istat)
 Vkb%ff(:,:,:,:)=zero ; Vkb%ffd(:,:,:,:)=zero
 
 ider=1 ; dimffnl=2 ! To retrieve the first derivative.
 idir=0 ; nkpg=0 ; ikg=0

 do ikibz=1,Kmesh%nibz

  npw_k = npwarr(ikibz)
  kpoint => Kmesh%ibz(:,ikibz)
  kg_k => kg(:,1+ikg:ikg+npw_k)
  ikg=ikg+npw_k
  !
  ! Quantities used only if useylm==1
  allocate(ylm(npw_k,Psps%mpsang**2*Psps%useylm))
  allocate(ylm_gr(npw_k,3+6*(ider/2),Psps%mpsang**2*Psps%useylm))
  allocate(ylm_k(npw_k,Psps%mpsang**2*Psps%useylm))
  allocate(kpg_dum(npw_k,nkpg)) 

  allocate(ffnl(npw_k,dimffnl,Psps%lmnmax,ntypat))

  call mkffnl(Psps%dimekb,dimffnl,Psps%ekb,ffnl,Psps%ffspl,gmet,gprimd,ider,idir,Psps%indlmn,&
&       kg_k,kpg_dum,kpoint,Psps%lmnmax,Psps%lnmax,Psps%mpsang,Psps%mqgrid_ff,nkpg,npw_k,& 
&       ntypat,Psps%pspso,Psps%qgrid_ff,rmet,Psps%usepaw,Psps%useylm,ylm_k,ylm_gr)

  deallocate(kpg_dum,ylm,ylm_gr,ylm_k)

  allocate(modkplusg(npw_k))
  effmass=one ; ecutsm=zero

  call mkkin(ecut,ecutsm,effmass,gmet,kg_k,modkplusg,kpoint,npw_k)
                                                                         
  modkplusg(:)=SQRT(half/pi**2*modkplusg(:))
  modkplusg(:)=MAX(modkplusg(:),tol10)

  do itypat=1,ntypat
   nlmn=COUNT(Psps%indlmn(3,:,itypat)>0)
   iln0=0 
   do ilmn=1,nlmn
    il= Psps%indlmn(1,ilmn,itypat)+1
    iln=Psps%indlmn(5,ilmn,itypat)

    if (iln>iln0) then
     iln0=iln

     if (ABS(Psps%ekb(ilmn,itypat))>1.0d-10) then
      SELECT CASE (il)
      CASE (1)
       Vkb%ff (1:npw_k,iln,itypat,ikibz) = ffnl(:,1,ilmn,itypat)
       Vkb%ffd(1:npw_k,iln,itypat,ikibz) = ffnl(:,2,ilmn,itypat)*modkplusg(:)/two_pi

      CASE (2)
       Vkb%ff (1:npw_k,iln,itypat,ikibz) =   ffnl(:,1,ilmn,itypat)*modkplusg(:)
       Vkb%ffd(1:npw_k,iln,itypat,ikibz) = ((ffnl(:,2,ilmn,itypat)*modkplusg(:)**2)+&
&                                            ffnl(:,1,ilmn,itypat))/two_pi
      CASE (3)
       Vkb%ff (1:npw_k,iln,itypat,ikibz) =  ffnl(:,1,ilmn,itypat)*modkplusg(:)**2
       Vkb%ffd(1:npw_k,iln,itypat,ikibz) = (ffnl(:,2,ilmn,itypat)*modkplusg(:)**3+&
&                                         2*ffnl(:,1,ilmn,itypat)*modkplusg(:))/two_pi
      CASE (4)
       Vkb%ff (1:npw_k,iln,itypat,ikibz) =  ffnl(:,1,ilmn,itypat)*modkplusg(:)**3
       Vkb%ffd(1:npw_k,iln,itypat,ikibz) = (ffnl(:,2,ilmn,itypat)*modkplusg(:)**4+&
                                          3*ffnl(:,1,ilmn,itypat)*modkplusg(:)**2)/two_pi
      CASE DEFAULT
       write(msg,'(4a)')ch10,&
&       ' vkb_init: ERROR: - ',ch10,&
&       ' l greater than g is not implemented '
        call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
      END SELECT

      fact = SQRT(four_pi/ucvol*(2*il-1)*ABS(Psps%ekb(ilmn,itypat)))
      Vkb%ff (:,iln,itypat,ikibz) = fact * Vkb%ff (:,iln,itypat,ikibz)
      Vkb%ffd(:,iln,itypat,ikibz) = fact * Vkb%ffd(:,iln,itypat,ikibz)

     else ! ekb==0
      Vkb%ff (:,iln,itypat,ikibz)=zero
      Vkb%ffd(:,iln,itypat,ikibz)=zero
     end if

    end if
   end do
  end do

  deallocate(ffnl,modkplusg)

 end do !ikibz

end subroutine vkb_init

END MODULE m_commutator_vkbr
!!***
