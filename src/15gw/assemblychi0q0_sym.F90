!{\src2tex{textfont=tt}}
!!****f* ABINIT/assemblychi0q0_sym
!! NAME
!! assemblychi0q0_sym
!!
!! FUNCTION
!! Update the independent particle susceptibility at q==0 for the contribution
!! of one pair of occupied-unoccupied band, for each frequencies.
!! This routine take advantage of the symmetries of the little group of the external point q
!! to symmetrize the contriburion coming from the input k point in the IBZ
!!
!! Compute chi0(G,G'',io)=chi0(G,G'',io)+\sum_S (rhotwg(G)*rhotwg*(G''))*green_w(io)
!! where S is a symmetry in reciprocal space 
!! The subroutine also performs the symmetrization of the matrix elements of the 
!! gradient operator and of the commutator of non local pseudopotential operator 
!! with the position operator 
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ik_bz= index of the k-point whose contribution has to be added to chi0, 
!!  isym_kbz=Index of the symmetry such as k = IS k_ibz
!!  itim_kbz=2 if time-reversal has to be used 
!!  npwepG0=maximum number of G vectors
!!  qpoint(3)=reciprocal space coordinates of the q wavevector
!!  rhotwg(npwe)=density of a pair of occupied-unoccupied states, in reciprocal space
!!  rhotwx(3)=matrix elements of the gradient and of the commutator of the non
!!   local operator with the position operator (the second term in included only if inclvkb=1
!!  green_w(nomega)=frequency dependent part coming from the green function
!!  Ltg_q<little_group_type>=Info on the little group associated to the external q-point.
!!   %timrev=2 it time-reversal is used, 1 otherwise
!!   %nsym_sg=Number of space group symmetries
!!   %wtksym(2,nsym,nkbz)=1 if the symmetry (with or without time-reversal) must be considered for this k-point
!!  Gsph_epsG0<Gvectors_type> Information on the "enlarged" G-sphere used for chi0, it contains umklapp G0 vectors
!!   %ng=number of G vectors in the enlarged sphere, actually MUST be equal to the size of rhotwg
!!   %rottbm1(ng,2,nsym)=index of (IR)^{-1} G where I is the identity or the inversion 
!!   %phmGt(ng,nsym)=phase factors associated to non-simmorphic operations
!!  Cryst<Crystal_structure>=Structure defining the unit cell and its symmetries
!!     %nsym=Number of symmetries
!!     %symrec(3,3,nsym)=Symmetry operation in reciprocal space (reduced coordinates)
!!  Ep<Epsilonm1_parameters>
!!     %npwe=number of plane waves in chi0
!!     %symchi
!!     %nomega=number of frequencies
!!    
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0(npwe,npwe,nomega)=independent-particle susceptibility matrix in reciprocal space at q==0
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

subroutine assemblychi0q0_sym(qpoint,ik_bz,isym_kbz,itim_kbz,gwcomp,nspinor,npwepG0,Ep,Cryst,Ltg_q,Gsph_epsG0,&
& chi0,rhotwx,rhotwg,green_w,green_enhigh_w,deltaf_b1b2)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : GW_TOL_DOCC, czero_gw, j_gw
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_15gw, except_this_one => assemblychi0q0_sym
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,isym_kbz,itim_kbz
 integer,intent(in) :: npwepG0,nspinor,gwcomp
 real(dp),intent(in) :: deltaf_b1b2
 type(Little_group),intent(in) :: Ltg_q
 type(Gvectors_type),intent(in) :: Gsph_epsG0 
 type(Crystal_structure),intent(in) :: Cryst
 type(Epsilonm1_parameters),intent(in) :: Ep
!arrays
 real(dp),intent(in) :: qpoint(3)
 complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
 complex(gwpc),intent(inout) :: chi0(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega)
 complex(dpc),intent(in) :: green_w(Ep%nomega),green_enhigh_w(Ep%nomega)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym,igp,ig
 integer :: jj,ii,s_jj,pad_jj,pad_ii
 complex(gwpc) :: dd
 character(len=500) :: msg
!arrays
 integer,pointer :: Sm1G(:) 
 real(dp) :: opinv(3,3),qrot(3)
 real(dp) :: b1(3),b2(3),b3(3)
 complex(gwpc),allocatable :: rhotwg_sym(:),rhotwg_I(:),rhotwg_J(:)
 complex(gwpc),pointer :: phmGt(:)

!************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'CGERC' :: cgerc
#endif

 b1(:)=two_pi*Gsph_epsG0%gprimd(:,1)
 b2(:)=two_pi*Gsph_epsG0%gprimd(:,2)
 b3(:)=two_pi*Gsph_epsG0%gprimd(:,3)

 SELECT CASE (Ep%symchi)
 
 CASE (0)
   if (nspinor==1) then
   ! === Do not use symmetries ===
   ! * Accumulate over the full BZ i.e chi0(G,Gp,io)=chi0(G,Gp,io)+(rhotwg(G)*rhotwg(Gp))*green_w(io)
   ! * The non-analytic term is symmetrized for this k-point in the BZ according to:
   !    rhotwg(1)= S^-1q*rhotwx_ibz
   !    rhotwg(1)=-S^-1q*CONJG(rhotwx_ibz) for inversion
   !    indeed -iq . <cSk|\nabla|vSk> = -i S^-1 q . <ck|\nabla|vk>  
   opinv(:,:)=REAL(Cryst%symrec(:,:,isym_kbz),dp)
   call matrginv(opinv,3,3) 
   call dosym(opinv,itim_kbz,qpoint,qrot)
   rhotwg(1)=dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3) !TODO get rid of this
   if (itim_kbz==2) rhotwg(1)=CONJG(rhotwg(1))

   if (gwcomp==1) then
    ! Leave the head and wings uncorrected (does not matter much)
    if (ABS(deltaf_b1b2) < GW_TOL_DOCC) rhotwg(1)=czero_gw
    do igp=1,Ep%npwe
     chi0(1,igp,:) = chi0(1,igp,:) + rhotwg(1) *CONJG(rhotwg(igp))*green_enhigh_w(:)
    end do
    do ig=2,Ep%npwe
     chi0(ig,1,:)  = chi0(ig,1,:)  + rhotwg(ig)*CONJG(rhotwg(1))  *green_enhigh_w(:)
    end do
   end if

   ! Multiply elements G,Gp of rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
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
   ! Important NOTE: treatment of q-->0 limit is correct only
   ! for i=j=0. Other components require additional terms.

   do jj=1,Ep%nJ
    s_jj=1 ; if (jj==4) s_jj=-1
    pad_jj=(jj-1)*Ep%npwe
    call mkrhotwg_sigma(jj,nspinor,Ep%npwe,rhotwg,rhotwg_J)

    rhotwg_J(1) = q0limit(jj,qpoint,nspinor,rhotwx,b1,b2,b3) 
    !TODO RECHECK this
    if (itim_kbz==2) rhotwg_J(1)=-CONJG(rhotwg_J(1))

    do ii=1,Ep%nI
     pad_ii=(ii-1)*Ep%npwe

     if (ii/=jj) then
      call mkrhotwg_sigma(ii,nspinor,Ep%npwe,rhotwg,rhotwg_I)
      rhotwg_I(1) = q0limit(ii,qpoint,nspinor,rhotwx,b1,b2,b3) 
      if (itim_kbz==2) rhotwg_I(1)=-CONJG(rhotwg_I(1))
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
  ! === Notes on the symmetrization of the oscilator matrix elements ===
  ! If  Sq = q then  M_G( Sk,q)= e^{-i(q+G).t} M_{ S^-1G}  (k,q)
  ! If -Sq = q then  M_G(-Sk,q)= e^{-i(q+G).t} M_{-S^-1G}^*(k,q)
  !
  ! In case of an umklapp process 
  ! If  Sq = q+G0 then  M_G( Sk,q)= e^{-i(q+G).t} M_{ S^-1(G-G0}   (k,q)
  ! If -Sq = q+G0 then  M_G(-Sk,q)= e^{-i(q+G).t} M_{-S^-1(G-G0)}^*(k,q)
  !
  ! Note that there is no need to take into account the phases due to q, 
  ! They cancel in the scalar product ==> phmGt(G,isym)=e^{-iG.t}
  !
  ! Mind the slicing of %rottbm1(npwepG0,timrev,nsym) and %phgt(npwepG0,nsym) as 
  ! these arrays, usually, do not conform to rho_twg_sym(npw) !
  !
  ! rhotwg(1)= R^-1q*rhotwx_ibz
  ! rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion

  allocate(rhotwg_sym(Ep%npwe))

  ! === Loop over symmetries of the space group and time-reversal ===
  do isym=1,Ltg_q%nsym_sg
   do itim=1,Ltg_q%timrev

    if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then 
     ! === This operation belongs to the little group and has to be considered to reconstruct the BZ ===
     ! TODO this is a hot-spot, should add a test on the umklapp
     !
     phmGt => Gsph_epsG0%phmGt  (1:Ep%npwe,isym) ! In these 3 lines mind the slicing (1:npwe)
     Sm1G =>  Gsph_epsG0%rottbm1(1:Ep%npwe,itim,isym)

     SELECT CASE (itim)

     CASE (1)
      rhotwg_sym(1:Ep%npwe)=rhotwg(Sm1G(1:Ep%npwe))*phmGt(1:Ep%npwe) 
      opinv(:,:)=REAL(Cryst%symrec(:,:,isym),dp)
      call matrginv(opinv,3,3) ; call dosym(opinv,itim,qpoint,qrot)
      rhotwg_sym(1)=dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3)

     CASE (2) 
      rhotwg_sym(1:Ep%npwe)=CONJG(rhotwg(Sm1G(1:Ep%npwe)))*phmGt(1:Ep%npwe) 
      opinv(:,:)=REAL(Cryst%symrec(:,:,isym),dp)
      call matrginv(opinv,3,3) ; call dosym(opinv,itim,qpoint,qrot)
      rhotwg_sym(1)=CONJG(dotproductqrc(qrot,rhotwx(:,1),b1,b2,b3))

     CASE DEFAULT 
      call assert(.FALSE.,'Wrong value of timrev',__FILE__,__LINE__)
     END SELECT 

     if (gwcomp==1) then
      ! Leave the head and wings uncorrected (does not matter much)
      if (ABS(deltaf_b1b2) < GW_TOL_DOCC) rhotwg_sym(1)=czero_gw
      do igp=1,Ep%npwe
       chi0(1,igp,:) = chi0(1,igp,:) + rhotwg_sym(1) *CONJG(rhotwg_sym(igp))*green_enhigh_w(:)
      end do
      do ig=2,Ep%npwe
       chi0(ig,1,:)  = chi0(ig,1,:)  + rhotwg_sym(ig)*CONJG(rhotwg_sym(1))  *green_enhigh_w(:)
      end do
     end if

     ! Multiply elements G,Gp of rhotwg_sym by green_w(io) and accumulate in chi0(G,Gp,io)
     do io=1,Ep%nomega
      dd=green_w(io) 
#if defined HAVE_GW_DPC
      call ZGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
#else
      call CGERC(Ep%npwe,Ep%npwe,dd,rhotwg_sym,1,rhotwg_sym,1,chi0(:,:,io),Ep%npwe)
#endif
     end do

    end if !wtksym
   end do !inv
  end do !isym
 
  deallocate(rhotwg_sym)

 CASE DEFAULT
  call assert(.FALSE.,'Wrong value of symchi',__FILE__,__LINE__)
 END SELECT

end subroutine assemblychi0q0_sym
!!***

!TODO this should be "contained" to facilitate inlining but abilint crashes, dont know why!
function q0limit(ii,qpoint,nspinor,rhotwx,b1,b2,b3)

 use defs_basis
 use m_gwdefs, only : j_gw


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_15gw, except_this_one => q0limit
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ii,nspinor
 complex(gwpc) :: q0limit
!arrays
 real(dp),intent(in) :: qpoint(3)
 real(dp),intent(in) :: b1(3),b2(3),b3(3)
 complex(gwpc),intent(in) :: rhotwx(3,nspinor**2)
! *********************************************************************

 SELECT CASE (ii)
 CASE (1)
  ! M_0(q-->0) = Lim M(up,up)+M(dwn,dwn). Exact, neglecting Vnl
  q0limit =  dotproductqrc(qpoint,rhotwx(:,1),b1,b2,b3) &
&           +dotproductqrc(qpoint,rhotwx(:,2),b1,b2,b3) 
 CASE (2) 
  ! M_z(q-->0) = Lim M(up,up)-M(dwn,dwn). 
  ! WARNING off-diagonal elements of rV12 and rV12 are neglected
  q0limit =  dotproductqrc(qpoint,rhotwx(:,1),b1,b2,b3) &
&           -dotproductqrc(qpoint,rhotwx(:,2),b1,b2,b3) 
 CASE (3)
  ! M_x(q-->0) = M(up,dwn)+M(dwn,up). 
  ! Both diagonal elements of the form v12r-rv21 and similiar terms in 12 and 21 are neglected
  q0limit =  dotproductqrc(qpoint,rhotwx(:,3),b1,b2,b3) &
&           +dotproductqrc(qpoint,rhotwx(:,4),b1,b2,b3) 
 CASE (4)
  ! Both diagonal elements of the form v12r-rv21 and similiar terms in 12 and 21 are neglected
  q0limit =( dotproductqrc(qpoint,rhotwx(:,3),b1,b2,b3) &
&           -dotproductqrc(qpoint,rhotwx(:,4),b1,b2,b3) )*j_gw
 END SELECT

end function q0limit
!!***
