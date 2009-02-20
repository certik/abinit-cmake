!{\src2tex{textfont=tt}}
!!****f* ABINIT/assemblychi0_sym
!! NAME
!! assemblychi0_sym
!!
!! FUNCTION
!! Update the independent particle susceptibility for the contribution
!! of one pair of occupied-unoccupied band, for each frequency.
!! If symchi=1 the expression is symmetrized taking into account the symmetries 
!! of the little group associated to the external q-point.
!!
!! Compute chi0(G1,G2,io)=chi0(G1,G2,io)+\sum_S \hat S (rhotwg(G1)*rhotwg*(G2))*green_w(io)
!! where S are the symmetries of the little group associated to the external q-point.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  nspinor=Number of spinorial components.
!!  ik_bz=index of the k-point in the BZ array whose contribution has to be symmetrized and added to cchi0 
!!  npwepG0=maximum number of G vectors taking into account possible umklapp G0, ie enlarged sphere G-G0
!!  rhotwg(npwe*nspinor**2)=Oscillator matrix elements for this k-point and the transition that has to be summed
!!  green_w(nomega)=frequency dependent part coming from the green function
!!  Gsph_epsG0<Gvectors_type> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!   %ng=number of G vectors in the enlarged sphere, actually MUST be equal to the size of rhotwg
!!   %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion 
!!   %phmGt(ng,nsym)=phase factors associated to non-simmorphic operations
!!  Ltg_q<little_group_type>=Info on the little group associated to the external q-point.
!!   %timrev=2 it time-reversal is used, 1 otherwise
!!   %nsym_sg=Number of space group symmetries
!!   %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point
!!   %flag_umklp(timrev,nsym)= flag for umklapp processes 
!!    if 1 that the particular operation (IS) requires a G_o to preserve Q, 0 otherwise 
!!   %igmG0(npwepG0,timrev,nsym) index of G-G0 in the array gvec
!!  Ep<Epsilonm1_parameters>=Parameters related to the calculation of chi0/epsilon^-1
!!      %symchi
!!      %nomega=number of frequencies
!!      %npwe=number of plane waves for epsilon (input variable)
!!    
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0(npwe,npwe,nomega)=independent-particle susceptibility matrix in reciprocal space
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

subroutine assemblychi0_sym(ik_bz,nspinor,Ep,Ltg_q,green_w,npwepG0,rhotwg,Gsph_epsG0,chi0)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => assemblychi0_sym
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,npwepG0,nspinor
 type(Gvectors_type),intent(in) :: Gsph_epsG0 
 type(Little_group),intent(in) :: Ltg_q
 type(Epsilonm1_parameters),intent(in) :: Ep
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 complex(dpc),intent(in) :: green_w(Ep%nomega)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: ig1,ig2,ig01,ig02,ig03,itim,io,isym
 integer :: jj,ii,s_jj,pad_jj,pad_ii
 complex(gwpc) :: dd 
 character(len=500) :: msg
!arrays
 integer :: g0(3)
 integer,pointer :: gmG0(:) 
 integer,allocatable :: Sm1_gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:),rhotwg_I(:),rhotwg_J(:)
 complex(gwpc),pointer :: phmGt(:)

! *************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'CGERC' :: cgerc
#endif
 
 SELECT CASE (Ep%symchi)

 CASE (0)
  ! === Do not use symmetries ===
  if (nspinor==1) then
   do io=1,Ep%nomega
    dd=green_w(io) 
#if defined HAVE_GW_DPC
    call ZGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
#else
    call CGERC(Ep%npwe,Ep%npwe,dd,rhotwg,1,rhotwg,1,chi0(:,:,io),Ep%npwe)
#endif
   end do
  else 
   allocate(rhotwg_I(Ep%npwe))
   allocate(rhotwg_J(Ep%npwe))

   ! I can use symmetries to loop over the upper triangle but 
   ! this makes using BLAS more difficult

   do jj=1,Ep%nJ
    s_jj=1 ; if (jj==4) s_jj=-1
    pad_jj=(jj-1)*Ep%npwe
    call mkrhotwg_sigma(jj,nspinor,Ep%npwe,rhotwg,rhotwg_J)

    do ii=1,Ep%nI
     pad_ii=(ii-1)*Ep%npwe

     if (ii/=jj) then
      call mkrhotwg_sigma(ii,nspinor,Ep%npwe,rhotwg,rhotwg_I)
     else 
      rhotwg_I(:)=rhotwg_J(:)
     end if

     do io=1,Ep%nomega
      dd = s_jj*green_w(io) 
#if defined HAVE_GW_DPC
      call ZGERC(Ep%npwe,Ep%npwe,dd,rhotwg_I,1,rhotwg_J,1,chi0(pad_ii+1:pad_ii+Ep%npwe,pad_jj+1:pad_jj+Ep%npwe,io),Ep%npwe)
#else
      call CGERC(Ep%npwe,Ep%npwe,dd,rhotwg_I,1,rhotwg_J,1,chi0(pad_ii+1:pad_ii+Ep%npwe,pad_jj+1:pad_jj+Ep%npwe,io),Ep%npwe)
#endif
     end do

    end do !ii
   end do !jj

   deallocate(rhotwg_I,rhotwg_J)
  end if

 CASE (1)
  !
  ! Notes on the symmetrization of the oscillator matrix elements
  !  If  Sq = q then  M_G^( Sk,q)= e^{-i(q+G).t} M_{ S^-1G}  (k,q)
  !  If -Sq = q then  M_G^(-Sk,q)= e^{-i(q+G).t} M_{-S^-1G}^*(k,q)
  !
  ! In case of an umklapp process 
  !  If  Sq = q+G0 then  M_G( Sk,q)= e^{-i(q+G).t} M_{ S^-1(G-G0}   (k,q)
  !  If -Sq = q+G0 then  M_G(-Sk,q)= e^{-i(q+G).t} M_{-S^-1(G-G0)}^*(k,q)
  !
  ! Ltg_q%igmG0(ig,itim,isym) contains the index of G-G0 where ISq=q+G0
  ! Note that there is no need to take into account the phases due to q, 
  ! They cancel in the scalar product ==> phmGt(G,isym)=e^{-iG\cdot t}
  !
  ! Mind the slicing of %rottbm1(npwepG0,timrev,nsym) and %phmGt(npwepG0,nsym) as 
  ! these arrays, usually, do not conform to rho_twg_sym(npw) !
  !
  allocate(rhotwg_sym(Ep%npwe))
  allocate(Sm1_gmG0  (Ep%npwe))
  !
  ! === Loop over symmetries of the space group and time-reversal ===
  do isym=1,Ltg_q%nsym_sg
   do itim=1,Ltg_q%timrev

    if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then 
     ! === This operation belongs to the little group and has to be used to reconstruct the BZ ===
     ! * In the following 3 lines mind the slicing (1:npwe)
     ! TODO this is a hot-spot, should add a test on the umklapp
     !
     phmGt => Gsph_epsG0%phmGt(1:Ep%npwe,isym) 
     gmG0  => Ltg_q%igmG0     (1:Ep%npwe,itim,isym)  
     Sm1_gmG0(1:Ep%npwe)=Gsph_epsG0%rottbm1(gmG0(1:Ep%npwe),itim,isym)

     SELECT CASE (itim)
     CASE (1) 
      rhotwg_sym(1:Ep%npwe)=rhotwg(Sm1_gmG0)*phmGt(1:Ep%npwe) 
     CASE (2) 
      rhotwg_sym(1:Ep%npwe)=CONJG(rhotwg(Sm1_gmG0))*phmGt(1:Ep%npwe) 
     CASE DEFAULT 
      call assert(.FALSE.,'Wrong timrev',__FILE__,__LINE__)
     END SELECT 

     ! Multiply rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
     do io=1,Ep%nomega
      dd=green_w(io)  
#if defined HAVE_GW_DPC
      call ZGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
#else
      call CGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
#endif
     end do

    end if
   end do 
  end do 

  deallocate(rhotwg_sym)
  deallocate(Sm1_gmG0)

 CASE DEFAULT
  call assert(.FALSE.,'Wrong value of symchi',__FILE__,__LINE__)
 END SELECT

end subroutine assemblychi0_sym
!!***



!!****f* ABINIT/mkrhotwg_sigma
!! NAME
!! mkrhotwg_sigma
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!    
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  umklapp are not allowed, npwe has to be equal to npwepG0.
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

subroutine mkrhotwg_sigma(ii,nspinor,npw,rhotwg,rhotwg_I)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_gwdefs, only : j_gw

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,npw,nspinor
!arrays
 complex(gwpc),intent(in) :: rhotwg(npw*nspinor**2)
 complex(gwpc),intent(inout) :: rhotwg_I(npw)

!Local variables-------------------------------
!scalars
!arrays

! *************************************************************************

 SELECT CASE (ii)
 CASE (1) 
  ! M_0 = M(up,up)+M(dwn,dwn)
  rhotwg_I(:) = rhotwg(1:npw) + rhotwg(npw+1:2*npw)
 CASE (2) 
  ! M_z = M(up,up)-M(dwn,dwn)
  rhotwg_I(:) = rhotwg(1:npw) - rhotwg(npw+1:2*npw)
 CASE (3) 
  ! M_x = M(up,dwn)+M(dwn,up)
  rhotwg_I(:) = ( rhotwg(2*npw+1:3*npw) + rhotwg(3*npw+1:4*npw) )
 CASE (4) 
  ! M_y = i*(M(up,dwn)-M(dwn,up))
  rhotwg_I(:) = (rhotwg(2*npw+1:3*npw) - rhotwg(3*npw+1:4*npw) )*j_gw
 CASE DEFAULT 
  call assert(.FALSE.,'Wrong i value',__FILE__,__LINE__) 
 END SELECT

end subroutine mkrhotwg_sigma
!!***
