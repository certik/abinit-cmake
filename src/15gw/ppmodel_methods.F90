!{\src2tex{textfont=tt}}
!!****f* ABINIT/ppmodel_methods
!! NAME
!! ppmodel_methods
!!
!! FUNCTION
!!  This module (?) contains methods used to initialize and destroy a PPmodel object.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
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
!!***

!!****f* ABINIT/PPmodel_symmetrize
!! NAME
!!  PPmodel_symmetrize
!!
!! FUNCTION
!!  Symmetrize the plasmonpole matrix elements in the full BZ zone
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  PPm<PPmodel_type>=data type containing information on the plasmonpole technique 
!!  Gsph<Gvectors_type>=data related to the G-sphere
!!    %grottb
!!    %phmSGt 
!!  Qmesh<BZ_mesh_type>=Info on the q-mesh
!!    %nbz=number if q-points in the BZ
!!    %tab(nbz)=index of the symmeric q-point in the IBZ, for each point in the BZ
!!    %tabo(nbz)=the operation that rotates q_ibz onto \pm q_bz (depending on tabi) 
!!    %tabi(nbz)=-1 if time-reversal has to be considered, 1 otherwise
!!  iq_bz=Index of the q-point in the BZ where PPmodel parameters have to be symmetrized
!!
!! OUTPUT
!!  botsq 
!!  otq
!!  eig (only if PPm%ppmodel==3)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  In the present implementation we are not considering a possible umklapp vector G0.
!!  In this case,indeed, the equation is different since we have to consider G-G0. 
!!  There is however a check in sigma
!! 
!!  * Remember the symmetry properties of \tilde\espilon^{-1}
!!    If q_bz=Sq_ibz+G0:
!! 
!!    $\epsilon^{-1}_{SG1-G0,SG2-G0}(q_bz) = e^{+iS(G2-G1).\tau}\epsilon^{-1}_{G1,G2)}(q)
!!
!!    If time-reversal symmetry can be used then :
!!    $\epsilon^{-1}_{G1,G2}(-q_bz) = e^{+i(G1-G2).\tau}\epsilon^{-1}_{-S^{-1}(G1+Go),-S^{-1}(G2+G0)}^*(q)
!! 
!! * Notice that eig is used only if PPm%model==3
!!
!! TODO
!!  Symmetrization can be skipped if iq_bz correspond to a point in the IBZ
!!
!! PARENTS
!!
!! CHILDREN
!!  
!!
!! SOURCE

subroutine PPmodel_symmetrize(PPm,Gsph,Qmesh,iq_bz,botsq,otq,eig) 

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iq_bz 
 type(PPmodel_type),intent(in) :: PPm
 type(Gvectors_type),intent(in) :: Gsph
 type(BZ_mesh_type),intent(in) :: Qmesh
!arrays
 complex(gwpc),intent(out) :: botsq(PPm%npwc,PPm%dm2_botsq)
 ! FIXME NOTE here otq is defined as real but omegatw is complex, this has to be investigated
 real(dp),intent(out) :: otq(PPm%npwc,PPm%dm2_otq)
 complex(gwpc),intent(out) :: eig(PPm%dm_eig,PPm%dm_eig) 

!Arguments ------------------------------------
!scalars
 integer :: ii,jj
 integer :: iq_ibz,iiq,isymq,iq_curr
 character(len=500) :: msg
!arrays
 integer,pointer :: grottb(:)
 complex(gwpc),pointer :: phsgt(:)
! *********************************************************************

 ! Here there is a problem with the small q, still cannot use BZ methods
 iq_ibz=Qmesh%tab(iq_bz) ; isymq=Qmesh%tabo(iq_bz) ; iiq=(3-qmesh%tabi(iq_bz))/2

 iq_curr=iq_ibz ; if (PPm%mqmem==0) iq_curr=1 

 grottb => Gsph%rottb(1:PPm%npwc,iiq,isymq)
 phsgt  => Gsph%phmSGt(1:PPm%npwc,isymq) 

 SELECT CASE (PPm%model)

 CASE (1,2)
  ! === Godby-Needs or Hybertsen-Louie PPmodel ===
  ! * Plasmon pole frequencies otq are obviously invariant under symmetry 
  do jj=1,PPm%npwc
   do ii=1,PPm%npwc
    botsq(grottb(ii),grottb(jj))=PPm%bigomegatwsq(ii,jj,iq_curr)*phsgt(ii)*CONJG(phsgt(jj))
    otq  (grottb(ii),grottb(jj))=PPm%omegatw(ii,jj,iq_curr)  
   end do
  end do

 CASE (3)
  ! === Symmetryze von der Linden-Horsh PPmodel ===
  ! * For notations see pag 22 of Quasiparticle Calculations in solid (Aulbur et al)
  !  If q_bz=Sq_ibz+G0 then:
  !
  ! $\omega^2_{ii}(q_bz) = \omega^2_{ii}(q)$        (otq array)
  ! $\alpha_{ii}(q_bz)   = \alpha_{ii}(q)$          (botq array
  ! $\Phi_{SG-G0}(q_bz)  = \Phi_{G}(q) e^{-iSG.t}$  (eigenvectors of e^{-1}, eig array) 
  !   
  do ii=1,PPm%npwc ! DM bands index
   otq  (ii,1)=PPm%omegatw     (ii,1,iq_curr)
   botsq(ii,1)=PPm%bigomegatwsq(ii,1,iq_curr)
   do jj=1,PPm%npwc
    eig(grottb(jj),ii)=PPm%eigpot(jj,ii,iq_curr)*phsgt(jj)
   end do
  end do
  if (iiq==2) eig(:,:)=CONJG(eig(:,:)) ! Time-reversal

 CASE (4)
  ! === Symmetrize Engel Farid PPmodel ===
  ! * For notations see pag 23 of Quasiparticle Calculations in solid (Aulbur et al)
  ! If q_bz=Sq_ibz+G0 then:
  !   
  ! $\omega^2_{ii}(q_bz) = \omega^2_{ii}(q)$        (otq array)
  ! $y_{SG-G0}(q_bz)     = y_{G}(q) e^{-iSG.t}$     (y=Lx) 
  !   
  do ii=1,PPm%npwc ! DM bands index
   otq(ii,1)=PPm%omegatw(ii,1,iq_curr)
   do jj=1,PPm%npwc
    botsq(grottb(jj),ii)=PPm%bigomegatwsq(jj,ii,iq_curr)*phsgt(jj)
   end do
  end do

 CASE DEFAULT 
  write(msg,'(4a,i2)')ch10,&
&  ' PPmodel_symmetrize: BUG- ',ch10,&
&  ' wrong value for PPm%model = ',PPm%model
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT
 !
 ! === Take into account time-reversal symmetry ===
 if (iiq==2) botsq(:,:)=CONJG(botsq(:,:))

end subroutine PPmodel_symmetrize
!!***

!!****f* ABINIT/nullify_PPmodel
!! NAME
!!  nullify_PPmodel
!!
!! FUNCTION
!!  Nullify dynamic entities in a PPmodel_type object
!!
!! SOURCE

subroutine nullify_PPmodel(PPm)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 type(PPmodel_type),intent(inout) :: PPm
! *********************************************************************

 nullify(PPm%bigomegatwsq)
 nullify(PPm%omegatw)
 nullify(PPm%eigpot)

end subroutine nullify_PPmodel
!!***

!!****f* ABINIT/destroy_PPmodel
!! NAME
!!  destroy_PPmodel
!!
!! FUNCTION
!!  Free a PPmodel structure
!!
!! SOURCE

subroutine destroy_PPmodel(PPm)

 use defs_basis
 use defs_datatypes

 implicit none

!Arguments ------------------------------------
 type(PPmodel_type),intent(inout) :: PPm
! *********************************************************************

 if (ASSOCIATED(PPm%bigomegatwsq)) deallocate(PPm%bigomegatwsq)
 if (ASSOCIATED(PPm%omegatw     )) deallocate(PPm%omegatw     )
 if (ASSOCIATED(PPm%eigpot      )) deallocate(PPm%eigpot      )

end subroutine destroy_PPmodel
!!***

!!****f* ABINIT/init_PPmodel
!! NAME
!!  init_PPmodel
!!
!! FUNCTION
!!  Initialize dimensions and other useful variables related to the PPmodel
!!
!! SOURCE

subroutine init_PPmodel(PPm,Qmesh,ppmodel,npwc,mqmem,drude_plsmf)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: mqmem,npwc,ppmodel
 real(dp),intent(in) :: drude_plsmf
 type(PPmodel_type),intent(out) :: PPm
 type(BZ_mesh_type),intent(in) :: Qmesh

!Local variables-------------------------------
 integer :: dim_q,istat
 logical :: ltest
 character(len=500) :: msg      
! *********************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' init_PPmodel : enter'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 PPm%nqibz = Qmesh%nibz
 PPm%mqmem = mqmem
 ltest = (mqmem==0.or.mqmem==Qmesh%nibz)
 call assert(ltest,'Wrong value for mqmem',__FILE__,__LINE__)
 PPm%npwc  = npwc
 PPm%model = ppmodel

 PPm%drude_plsmf = drude_plsmf

 SELECT CASE (PPm%model)
 CASE (0) 
  write(*,*)'WARNING, called with ppmodel==0'
  PPm%dm2_botsq = 0
  PPm%dm2_otq   = 0
  PPm%dm_eig    = 0 
  !call nullify_PPmodel(PPm)
  RETURN
 CASE (1,2)
  PPm%dm2_botsq = npwc
  PPm%dm2_otq   = npwc
  PPm%dm_eig    = 1 ! Should be set to 0, but g95 doesnt like zero-sized arrays
 CASE (3)
  PPm%dm2_botsq = 1  
  PPm%dm2_otq   = 1
  PPm%dm_eig    = npwc
 CASE (4)
  PPm%dm2_botsq = npwc
  PPm%dm2_otq   = 1
  PPm%dm_eig    = 1 ! Should be set to 0, but g95 doesnt like zero-sized arrays
 CASE DEFAULT
  write(msg,'(2a,i3)')ch10,&
&  ' init_PPmodel : PPmodel not allowed ',PPm%model
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT
 !
 ! === Do we store full the q-mesh or out-of-memory solution? ===
 dim_q=Qmesh%nibz ; if (PPm%mqmem==0) dim_q=1

 allocate(PPm%bigomegatwsq (npwc,PPm%dm2_botsq,dim_q), STAT=istat )
 allocate(PPm%omegatw      (npwc,PPm%dm2_otq,  dim_q), STAT=istat) 
 allocate(PPm%eigpot (PPm%dm_eig,PPm%dm_eig,dim_q), STAT=istat)

#if defined DEBUG_MODE
 write(msg,'(a)')' init_PPmodel : exit'
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 
end subroutine init_PPmodel
!!***
