!{\src2tex{textfont=tt}}
!!****f* ABINIT/assemblychi0sfq0
!! NAME
!! assemblychi0sfq0
!!
!! FUNCTION
!! Update the imaginary part of the independent particle susceptibility at q==0 for the contribution
!! of one pair of occupied-unoccupied band, for each frequencies
!! taking into account the symmetries of the little group of the external point q
!!
!! Compute chi0(G,G'',io)=chi0(G,G'',io)+\sum_S (rhotwg(G)*rhotwg*(G''))*den(io)
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
!!  ikbz=index of the k-point whose contribution to chi0 has to be added, 
!!   if we use symmetries the contribution to chi0 is symmetrized 
!!  isym_kbz=Index of the symmetry such as k = IS k_ibz
!!  itim_kbz=2 if time-reversal has to be used 
!!  my_wl,my_wr= upper and lower frequency considered by this processor
!!  npwe=number of plane waves for chi0
!!  npwepG0=maximum number of G vectors
!!  nomega=number of frequencies
!!  qpoint(3)=reciprocal space coordinates of the q wavevector
!!  rhotwg(npwe)=density of a pair of occupied-unoccupied states, in reciprocal space
!!  rhotwx(3)=matrix elements of the gradient and of the commutator of the non
!!   local operator with the position operator (the second term in included only if inclvkb=1
!!  den(nomega)=denominator of the susceptibility expression
!!  nsym=number of symmetry operations
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
!! Cryst<Crystal_structure>=Info on unit cell and it symmetries
!!   %nsym=Number of symmetry operations.
!!   %symrec(3,3,nsym)=Symmetry operations in reciprocal space (reduced coordinates).
!!    
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0(npwe,npwe,nomega)=independent-particle susceptibility matrix in reciprocal space at q==0
!!
!! NOTES
!!  Non symmporphic operations are not yet treated
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

subroutine assemblychi0sfq0(qpoint,ikbz,isym_kbz,itim_kbz,nspinor,symchi,npwepG0,npwe,Cryst,Ltg_q,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,rhotwx,rhotwg,nomegasf,chi0sf)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_15gw, except_this_one => assemblychi0sfq0
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikbz,my_wl,my_wr,nomegasf,npwe,npwepG0,nspinor
 integer,intent(in) :: isym_kbz,itim_kbz,symchi,iomegal,iomegar
 real(dp),intent(in) :: factocc,wl,wr 
 type(Little_group),intent(in) :: Ltg_q
 type(Gvectors_type),intent(in) :: Gsph_epsG0 
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 real(dp),intent(in) :: qpoint(3)
 complex(gwpc),intent(inout) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(in) :: rhotwx(3)
 complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)

!Local variables-------------------------------
!scalars
 integer :: itim,io,isym
 complex(gwpc) :: num 
 character(len=500) :: msg
!arrays
 integer,pointer :: Sm1G(:) 
 real(dp) :: opinv(3,3),qrot(3),b1(3),b2(3),b3(3)
 complex(gwpc),allocatable :: rhotwg_sym(:)
 complex(gwpc),pointer :: phmGt(:)
!************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'CGERC' :: cgerc
#endif

 if (iomegal<my_wl .or. iomegar>my_wr) then 
  write(msg,'(3a,2(a,i5,a,i5))')ch10,&
&  ' assemblychi0sfq0 : Indeces out of boundary ',ch10,&
&  '  my_wl = ',my_wl,' iomegal = ',iomegal,ch10,&
&  '  my_wr = ',my_wr,' iomegar = ',iomegar,ch10
   !write(msg,'(2f8.2,2i4)')ep%omegasf(my_wr)*Ha_eV,ep%omegasf(iomegar)*Ha_eV,ibv,ibc ; call wrtout(std_out,msg,'PERS')
  call wrtout(std_out,msg,'PERS') ; call leave_new('COLL')
 end if 

 b1(:)=two_pi*Gsph_epsG0%gprimd(:,1)
 b2(:)=two_pi*Gsph_epsG0%gprimd(:,2)
 b3(:)=two_pi*Gsph_epsG0%gprimd(:,3)

 SELECT CASE (symchi)

 CASE (0)
  ! 
  ! === Calculation without symmetries ===
  ! * rhotwg(1)= R^-1q*rhotwx_ibz
  ! * rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion
  ! FIXME My equation reads  -iq* <cSk|\nabla|vSk> = -i \transpose S <ck_i|\nabla\|vk_i>  
  opinv(:,:)=REAL(Cryst%symrec(:,:,isym_kbz),dp)
  call matrginv(opinv,3,3)
  call dosym(opinv,itim_kbz,qpoint,qrot)
  rhotwg(1)=dotproductqrc(qrot,rhotwx,b1,b2,b3)
  if (itim_kbz==2) rhotwg(1)=CONJG(rhotwg(1))

  if (wl<huge(0.0_dp)*1.d-11) then !this is awful but it is still a first coding
   num=-wl*factocc ! Num is single precision needed for cgerc check factocc
#if defined HAVE_GW_DPC
   call ZGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegal),npwe)
#else
   call CGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegal),npwe)
#endif
  end if 
  ! Last point, must accumulate left point but not the right one
  if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then 
   num=-wr*factocc
#if defined HAVE_GW_DPC
   call ZGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegar),npwe)
#else
   call CGERC(npwe,npwe,num,rhotwg,1,rhotwg,1,chi0sf(:,:,iomegar),npwe)
#endif
  end if 

 CASE (1)
  ! === Notes on the symmetrization of oscillator matrix elements ===
  ! If  Sq = q then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1G}  (k,q)
  ! If -Sq = q then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1G}^*(k,q)
  ! 
  ! In case of an umklapp process 
  ! If  Sq = q+G_o then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1(G-G_o}   (k,q)
  ! If -Sq = q+G_o then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1(G-G-o)}^*(k,q)
  ! 
  ! rhotwg(1)= R^-1q*rhotwx_ibz
  ! rhotwg(1)=-R^-1q*conjg(rhotwx_ibz) for inversion
  !
  allocate(rhotwg_sym(npwe))
  !
  ! === Loop over symmetries of the space group and time-reversal ===
  do isym=1,Ltg_q%nsym_sg
   do itim=1,Ltg_q%timrev

    if (Ltg_q%wtksym(itim,isym,ikbz)==1) then 
     ! === This operation belongs to the little group and has to be considered to reconstruct the BZ ===
     ! TODO this is a hot-spot, should add a test on the umklapp
     !
     phmGt => Gsph_epsG0%phmGt(1:npwe,isym) ! In these 2 lines mind the slicing (1:npwe)
     Sm1G =>  Gsph_epsG0%rottbm1(1:npwe,itim,isym)

     SELECT CASE (itim)
     CASE (1)
      rhotwg_sym(1:npwe)=rhotwg(Sm1G(1:npwe))*phmGt(1:npwe)
      opinv(:,:)=REAL(Cryst%symrec(:,:,isym),dp)
      call matrginv(opinv,3,3)
      call dosym(opinv,itim,qpoint,qrot)
      rhotwg_sym(1)=dotproductqrc(qrot,rhotwx,b1,b2,b3)
     CASE (2) 
      rhotwg_sym(1:npwe)=CONJG(rhotwg(Sm1G(1:npwe)))*phmGt(1:npwe)
      opinv(:,:)=REAL(Cryst%symrec(:,:,isym),dp)
      call matrginv(opinv,3,3)
      call dosym(opinv,itim,qpoint,qrot)
      rhotwg_sym(1)=CONJG(dotproductqrc(qrot,rhotwx,b1,b2,b3))
     CASE DEFAULT
      call assert(.FALSE.,'Wrong value of timrev',__FILE__,__LINE__)
     END SELECT
     !
     ! === Multiply elements G,Gp of rhotwg_sym*num and accumulate in chi0sf(G,Gp,io) ===
     if (wl<huge(0.0_dp)*1.d-11) then
      num=-wl*factocc
#if defined HAVE_GW_DPC
      call ZGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegal),npwe)
#else
      call CGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegal),npwe)
#endif
     end if
     ! Last point, must accumulate left point but not the right one
     if (iomegar/=nomegasf+1 .and. wr<huge(0.0_dp)*1.d-11) then 
      num=-wr*factocc
#if defined HAVE_GW_DPC
      call ZGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegar),npwe)
#else
      call CGERC(npwe,npwe,num,rhotwg_sym,1,rhotwg_sym,1,chi0sf(:,:,iomegar),npwe)
#endif
     end if 
    end if !wtksym

   end do !inv
  end do !isym
  deallocate(rhotwg_sym)

 CASE DEFAULT
  call assert(.FALSE.,'Wrong value of symchi',__FILE__,__LINE__)
 END SELECT

end subroutine assemblychi0sfq0
!!***
