!{\src2tex{textfont=tt}}
!!****f* ABINIT/assemblychi0sf
!! NAME
!! assemblychi0sf
!!
!! FUNCTION
!! Update the imaginary part of the polarizability for the contribution
!! of one pair of occupied-unoccupied band, for each frequencies
!! taking into account the symmetries of the little group of the external point q
!!
!! Compute chi0(G,Gp,io)=chi0(G,Gp,io)+\sum_S (rhotwg(G)*rhotwg*(G''))*den(io)
!! where S are the symmetries of the little group of the external q-point.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  ik_bz=index of the k-point in the bz array whose contribution has to be added to cchi0 
!!   after symmetrization of the matrix elements. 
!!  my_wl,my_wr= upper and lower frequency considered by this processor
!!  npwe=number of plane waves in chi0
!!  npwepG0=maximum number of G vectors
!!  nomegasf=number of frequencies for imaginary part of chi0
!!  nspinor=Number of spinorial components.
!!  symchi= 1 if symmetries can be usd, 0 otherwise
!!  rhotwg(npwe*nspinor**2)=density of a pair of occupied-unoccupied states, in reciprocal space
!!  timrev=if 2, inversion is considered; if 1, inversion is not considered
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
!!  factocc=occupation factor= f_occ*(ockp-occk) (see cchi0.F90)  
!!  wl,wr = weights used to approximate the Dirac function 
!!    
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  chi0sf(npwe,npwe,my_wl:my_wr)= update the imaginary part of the independent-particle 
!!   susceptibility matrix in reciprocal space
!!
!! NOTES
!!  Umklapp processes are not yet implemented 
!! 
!! PARENTS
!!      cchi0
!!
!! CHILDREN
!!      cgerc,wrtout,zgerc
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine assemblychi0sf(ik_bz,nspinor,symchi,Ltg_q,npwepG0,npwe,rhotwg,Gsph_epsG0,&
& factocc,my_wl,iomegal,wl,my_wr,iomegar,wr,nomegasf,chi0sf)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ik_bz,iomegal,iomegar,my_wl,my_wr,nomegasf,npwe,npwepG0
 integer,intent(in) :: nspinor,symchi
 real(dp),intent(in) :: factocc,wl,wr
 type(Gvectors_type),intent(in) :: Gsph_epsG0
 type(Little_group),intent(in) :: Ltg_q
!arrays
 complex(gwpc),intent(in) :: rhotwg(npwepG0*nspinor**2)
 complex(gwpc),intent(inout) :: chi0sf(npwe,npwe,my_wl:my_wr)

!Local variables-------------------------------
!scalars
 integer :: io,isym,itim
 complex(gwpc) :: num
 character(len=500) :: msg
!arrays
 integer,allocatable :: Sm1_gmG0(:)
 integer,pointer :: gmG0(:)
 complex(gwpc),allocatable :: rhotwg_sym(:)
 complex(gwpc),pointer :: phmGt(:)

! *************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'CGERC' :: cgerc
#endif

 if (iomegal < my_wl .or. iomegar > my_wr) then 
   write(msg,'(3a,2(a,i5,a,i5))')ch10,&
&   ' assemblychi0sf : Indeces out of boundary ',ch10,&
&   '  my_wl = ',my_wl,' iomegal = ',iomegal,ch10,&
&   '  my_wr = ',my_wr,' iomegar = ',iomegar,ch10
   call wrtout(std_out,msg,'PERS') ; call leave_new('COLL')
   !write(msg,'(2f8.2,2i4)')ep%omegasf(my_wl)*Ha_eV,ep%omegasf(iomegal)*Ha_eV,ibv,ibc
   !call wrtout(std_out,msg,'PERS')
   !write(msg,'(2f8.2,2i4)')ep%omegasf(my_wr)*Ha_eV,ep%omegasf(iomegar)*Ha_eV,ibv,ibc ; call wrtout(std_out,msg,'PERS')
 end if 

 SELECT CASE (symchi)
 CASE (0)
  !
  ! === Do not use symmetries ===
  if (wl<huge(0.0_dp)*1.d-11) then !FIXME this is awful
   num=-wl*factocc 
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
  ! 
  ! Use symmetries to reconstruct oscillator matrix elements
  ! Notes on the symmetrization of the oscillator maxtri elements:
  ! 
  ! If  Sq=q then  M_G^( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1G}  (k,q)
  ! If -Sq=q then  M_G^(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1G}^*(k,q)
  ! 
  ! In case of an umklapp process 
  ! If  Sq=q+G_o then  M_G( Sk,q)= e^{-i(q+G)\cdot t} M_{ S^-1(G-G_o}   (k,q)
  ! If -Sq=q+G_o then  M_G(-Sk,q)= e^{-i(q+G)\cdot t} M_{-S^-1(G-G_o)}^*(k,q)
  !
  ! Ltg_q%igmG0(ig,itim,isym) contains the index of G-G0 where ISq=q+G0
  ! Note that there is no need to take into account the phases due to q, 
  ! They cancel in the scalar product ==> phmGt(G,isym)=e^{-iG\cdot t}
  ! 
  ! Mind the slicing of %rottbm1(npwepG0,timrev,nsym) and %phgt(npwepG0,nsym) as 
  ! these arrays, usually, do not conform to rho_twg_sym(npw) !
  ! 
  allocate(rhotwg_sym(npwe))
  allocate(Sm1_gmG0(npwe))
  !
  ! === Loop over symmetries of the space group and time-reversal ===
  do isym=1,Ltg_q%nsym_sg
   do itim=1,Ltg_q%timrev

    if (Ltg_q%wtksym(itim,isym,ik_bz)==1) then 
     ! === This operation belongs to the little group and has to be used to reconstruct BZ ===
     ! TODO this is a hot-spot, should add a test on the umklapp
     !
     phmGt => Gsph_epsG0%phmGt(1:npwe,isym) ! In these 3 lines mind the slicing (1:npwe)
     gmG0 => Ltg_q%igmG0(1:npwe,itim,isym) 
     Sm1_gmG0(1:npwe)=Gsph_epsG0%rottbm1(gmG0(1:npwe),itim,isym)

     SELECT CASE (itim)
     CASE (1)
      !rhotwg_sym(1:npwe)=rhotwg(Sm1_gmG0)*Gsph_epsG0%phmGt(1:npwe,isym)
      rhotwg_sym(1:npwe)=rhotwg(Sm1_gmG0(1:npwe))*phmGt(1:npwe)
     CASE (2) 
      rhotwg_sym(1:npwe)=CONJG(rhotwg(Sm1_gmG0(1:npwe)))*phmGt(1:npwe)
     CASE DEFAULT
      write(msg,'(4a)')ch10,&
&      ' assemblychi0sf BUG-',ch10,&
&      ' called with wrong value of timrev '
      call wrtout(std_out,msg,'COLL') ; call leave_new('COLL') 
     END SELECT

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
  deallocate(rhotwg_sym,Sm1_gmG0)

 CASE DEFAULT
  write(msg,'(4a)')ch10,&
&  ' assemblychi0sf : BUG -',ch10,&
&  ' wrong value for symsigma '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL') 
 END SELECT

end subroutine assemblychi0sf
!!***
