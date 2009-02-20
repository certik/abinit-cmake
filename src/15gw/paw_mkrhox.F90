!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_mkrhox
!! NAME
!! paw_mkrhox
!!
!! FUNCTION
!!  Evaluate $<phj|e^{-i(q+G)}|phi>-<tphj|e^{-i(q+G)}|tphi>$
!!  for a fixed q-point and npw G vectors. Matrix elements are stored in packed storage mode. 
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  gmet(3,3)=reciprocal lattice metric tensor ($\textrm{Bohr}^{-2}$)
!!  gvec(3,npw)=G vectors in reduced coordinates
!!  Cryst<Crystal_structure>=data type gathering information on the unit cell
!!     %natom=number of atoms
!!     %typat(natom)=type of each atom
!!     %xred(3,natom)
!!     %ntypat=number of types of atoms
!!  npw=numper of G vectors
!!  Pawang<pawang_type> angular mesh discretization and related data:
!!     %gntselect(l_size_max**2,l_max**2*(l_max**2+1)/2)=Selection rules for Gaunt coefficients
!!     %l_max=Maximum value of angular momentum l+1
!!     %l_size_max=Maximum value of angular momentum l_size=2*l_max-1
!!     %ngnt=number of non-zero Gaunt coefficients
!!     %realgnt(Pawang%ngnt)=non-zero real Gaunt coefficients
!!  Psps<pseudopotential_type>:
!!     %lmnmax= Maximum number of different (l,m,n) components over all types of PAW dataset, same as dtset%lmnmax
!!     %lnmax=Max. number of (l,n) components over all type of PAW datasets
!!     %mqgrid_ff=Number of points in the reciprocal space grid on which the radial functions pwff_spl are specified
!!     %indlmn(6,lmnmax,ntypat) array giving l,m,n,lm,ln,spin for i=ln  (if useylm=0), or i=lmn (if useylm=1)
!!     %qgrid_ff(Psps%mqgrid_ff)=values at which form factors have been evaluated
!!     %mpsang=1+maximum angular momentum
!!  qpt(3)= q-point in reduced coordinates
!!  ylm_q(npw,(2*Psps%mpsang-1)**2)=real spherical harmonics Ylm(q+G) for q-point qpt up to l=2*l_max
!!  pwff_spl(Psps%mqgrid_ff,2,0:2*(Psps%mpsang-1),Psps%lnmax*(Psps%lnmax+1)/2,ntypat)) 
!!  Pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!     %lmnmix_sz=number of (lmn,lmn_prime) verifying l<=lmix and l_prime<=lmix
!!     %indklmn(4,lmn2_size)=array giving klm, kln, abs(il-jl) and (il+jl) for each klmn=(ilmn,jlmn)
!!     %dshpfunc(mesh_size,l_size,4)=derivatives of shape function (used only for numerical shape functions)
!!     %eijkl(lmn2_size,lmn2_size)=Part of the Dij that depends only from the projected occupation coeffs
!!     %exccore=Exchange-correlation energy for the core density
!!     %gnorm(l_size)=Normalization factor of radial shape function
!!     %phiphj(:,:)=Useful product Phi(:,i)*Phi(:,j)
!!     %qijl(l_size**2,lmn2_size)=Moments of the difference charge density between AE and PS partial wave
!!     %rad_for_spline(mesh_size)=radial grid used for spline (copy of pawrad%rad)
!!     %shapefunc(mesh_size,l_size)=Normalized radial shape function
!!     %sij(lmn2_size)=Nonlocal part of the overlap operator
!!     %tphitphj(:,:)=Useful product tPhi(:,i)*tPhi(:,j)
!!
!! OUTPUT
!!  paw_rhox(2,npw,Psps%lmnmax*(Psps%lmnmax+1)/2,natom): array containing 
!!   $<phj|e^{-i(q+G)}|phi>-<tphj|e^{-i(q+G)}|tphi>$ in packed form.
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      cchi0,cchi0q0,sigma
!!
!! CHILDREN
!!      realgaunt,splfit,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine paw_mkrhox(Cryst,pwff_spl,gmet,gvec,method,dim1_rhox,dim2_rhox,&
& Psps,Pawang,Pawtab,qpt,npw,ylm_q,paw_rhox)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13paw
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dim1_rhox,dim2_rhox,method,npw
 type(Crystal_structure),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: gvec(3,npw)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: pwff_spl(Psps%mqgrid_ff,2,0:dim1_rhox,dim2_rhox,Cryst%ntypat)
 real(dp),intent(in) :: qpt(3),ylm_q(npw,(2*Psps%mpsang-1)**2)
 real(dp),intent(out) :: paw_rhox(2,npw,Psps%lmnmax*(Psps%lmnmax+1)/2,Cryst%natom)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatm,ider,ig,ignt,ii,il,ilm,ilm_G,ilmn,iln,im,ipow,itypat,jj,jl,jlm
 integer :: jlmn,jln,jm,k0lm,k0lmn,k0ln,klm,klmn,kln,ll_G,mm_G,mpsang,ngnt
 real(dp) :: arg,rgnt
 character(len=500) :: msg
!arrays
 integer,allocatable :: gntselect(:,:)
 real(dp) :: mi_l(2,0:3),qpg(3),x0(3)
 real(dp),allocatable :: derfun(:),newfun(:),ph3d(:,:),qpg_norm(:),realgnt(:)
 real(dp),allocatable :: wk_ffnl(:,:)

! *************************************************************************

 mpsang=Psps%mpsang
 !
 ! === Pre-calculate (-i)^l ===
 mi_l(1,0)=one  ; mi_l(2,0)=zero
 mi_l(1,1)=zero ; mi_l(2,1)=-one
 mi_l(1,2)=-one ; mi_l(2,2)=zero
 mi_l(1,3)=zero ; mi_l(2,3)=one
 !
 ! === Calculate |q+G| ===
 ! * Do not include factor 2\pi to be consistent with spline 
 allocate(qpg_norm(npw))
 do ig=1,npw
  qpg(:)=qpt(:)+gvec(:,ig)
  qpg_norm(ig)=SQRT(DOT_PRODUCT(qpg,MATMUL(gmet,qpg)))
 end do
 ! Check q-grid as %qgrid_ff depends on ecut, not on ecutwfn, ecuteps or ecutsigx.
 if (MAXVAL(qpg_norm)>MAXVAL(Psps%qgrid_ff)) then 
  write(msg,'(6a,f8.4,a,f8.4,2a)')ch10,&
&  ' paw_mkrhox : ERROR ',ch10,&
&  ' Function values are being requested outside range of data. ',ch10,&
&  ' Max qpg_norm = ',MAXVAL(qpg_norm),' Max qgrid_ff = ',MAXVAL(Psps%qgrid_ff),ch10,&
&  ' Increase ecut(wfn), check qrid_ff and gsqcut '
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if

 allocate(wk_ffnl(Psps%mqgrid_ff,2))
 allocate(ph3d(2,npw),newfun(npw),derfun(npw))

 SELECT CASE (method) 
 CASE (1) 
  ! === Arnaud-Alouani exact expression ===
  ! * It does not descrive the multipoles of the AE charge density 
  ! * $ e^{-i(q+G).R_i} 4\pi \sum_{LM} (-i)^l Y_M^L(q+G) G_{\li\mi\lj\mj}^{\LM} ff^{aL}_{ij}(q+G) $
  !   where f has been calculated in paw_mkrhox_spl
  !
  ! === Re-evaluate Gaunt coefficients, just to be on the safe side ===
  ! * Note that gntselect is in packed form, thanks to invariance under permutation.
  ! * Could use Pawang% but size of gntselect depends on pawxcdev!
  allocate(  realgnt((2*mpsang-1)**2*(mpsang)**4))
  allocate(gntselect((2*mpsang-1)**2, mpsang**2*(mpsang**2+1)/2))
  call realgaunt(mpsang,ngnt,gntselect,realgnt)
  paw_rhox(:,:,:,:)=zero

  do iatm=1,Cryst%natom
   !
   ! === Structure factor e^{-i(q+G)*xred} ===
   x0(:)=Cryst%xred(:,iatm)
   do ig=1,npw
    qpg(:)=qpt(:)+gvec(:,ig)
    arg=-two_pi*DOT_PRODUCT(qpg(:),x0)
    ph3d(1,ig)=COS(arg)
    ph3d(2,ig)=SIN(arg)
   end do
   !
   ! === Loop on (jl,jm,jn) channels for this atom ===
   itypat=Cryst%typat(iatm)
   do jlmn=1,Pawtab(itypat)%lmn_size
    jl =Psps%indlmn(1,jlmn,itypat)
    jm =Psps%indlmn(2,jlmn,itypat)
    jlm=Psps%indlmn(4,jlmn,itypat)
    jln=Psps%indlmn(5,jlmn,itypat)
 
    k0lmn=jlmn*(jlmn-1)/2 
    k0lm =jlm *(jlm -1)/2
    k0ln =jln *(jln -1)/2
    !
    ! === Loop on (il,im,in) channels; klmn is index for packed form ===
    do ilmn=1,jlmn 
     il =Psps%indlmn(1,ilmn,itypat)
     im =Psps%indlmn(2,ilmn,itypat)
     ilm=Psps%indlmn(4,ilmn,itypat)
     iln=Psps%indlmn(5,ilmn,itypat)
 
     klmn=k0lmn+ilmn 
     klm =k0lm +ilm
     kln =k0ln +iln
     !
     ! === Summing over allowed (l,m), taking into account Gaunt selection rules ===
     do ll_G=ABS(jl-il),jl+il,2 
      ipow=MOD(ll_G,4) 
      ider=0 ; wk_ffnl(:,:)=pwff_spl(:,:,ll_G,kln,itypat) 
      call splfit(Psps%qgrid_ff,derfun,wk_ffnl,ider,qpg_norm,newfun,Psps%mqgrid_ff,npw)

      do mm_G=-ll_G,ll_G
       ilm_G=1+ll_G**2+ll_G+mm_G
       ignt=gntselect(ilm_G,klm) 
       if (ignt==0) CYCLE
       rgnt=realgnt(ignt)
       !ider=0 ; wk_ffnl(:,:)=pwff_spl(:,:,ll_G,kln,itypat) 
       !call splfit(Psps%qgrid_ff,derfun,wk_ffnl,ider,qpg_norm,newfun,Psps%mqgrid_ff,npw)
       !
       ! === Evaluate matrix elements for each plane wave ===
       do ig=1,npw
        paw_rhox(1,ig,klmn,iatm) = paw_rhox(1,ig,klmn,iatm)+ &
&        newfun(ig)*ylm_q(ig,ilm_G)*rgnt* (ph3d(1,ig)*mi_l(1,ipow)-ph3d(2,ig)*mi_l(2,ipow)) 

        paw_rhox(2,ig,klmn,iatm) = paw_rhox(2,ig,klmn,iatm)+ &
&        newfun(ig)*ylm_q(ig,ilm_G)*rgnt* (ph3d(1,ig)*mi_l(2,ipow)+ph3d(2,ig)*mi_l(1,ipow)) 
       end do

      end do !ll_G
     end do !mm_G
    end do !ilmn
   end do !jlmn
  end do !natom
  !
  ! * Multiply by 4\pi arising from the expansion of the plane wave
  paw_rhox=four_pi*paw_rhox
  deallocate(realgnt,gntselect)

 CASE (2) 
  ! === Shishkin-Kresse approximated expression ====
  ! * Better description of multipoles of AE charge, 
  ! * Better results for energy degeneracies in GW band structure 
  ! * $4\pi \sum_{LM} q_ij^{LM} Y_M^L(q+G) f^{aL}_{ij}(q+G)$ where f has been calculated in paw_mkrhox_spl
  !
  paw_rhox(:,:,:,:)=zero
  do iatm=1,Cryst%natom
   !
   ! === Structure factor e^{-i(q+G)*xred} ===
   x0(:)=Cryst%xred(:,iatm)
   do ig=1,npw
    qpg(:)=qpt(:)+gvec(:,ig)
    arg=-two_pi*DOT_PRODUCT(qpg(:),x0)
    ph3d(1,ig)=COS(arg)
    ph3d(2,ig)=SIN(arg)
   end do
   !
   ! === Loop on (jl,jm,jn) channels for this atom ===
   itypat=Cryst%typat(iatm)
   do jlmn=1,Pawtab(itypat)%lmn_size
    jl =Psps%indlmn(1,jlmn,itypat)
    jm =Psps%indlmn(2,jlmn,itypat)
    jlm=Psps%indlmn(4,jlmn,itypat)
    jln=Psps%indlmn(5,jlmn,itypat)

    k0lmn=jlmn*(jlmn-1)/2 
    k0lm =jlm *(jlm -1)/2
    k0ln =jln *(jln -1)/2
    !
    ! === Loop on (il,im,in) channels; klmn is index for packed form ===
    do ilmn=1,jlmn 
     il =Psps%indlmn(1,ilmn,itypat)
     im =Psps%indlmn(2,ilmn,itypat)
     ilm=Psps%indlmn(4,ilmn,itypat)
     iln=Psps%indlmn(5,ilmn,itypat)
 
     klmn=k0lmn+ilmn 
     klm =k0lm +ilm
     kln =k0ln +iln
     !
     ! === Summing over allowed (l,m), taking into account Gaunt selection rules ===
     do ll_G=ABS(jl-il),jl+il,2 
      ipow=MOD(ll_G,4) 
      do mm_G=-ll_G,ll_G

       ! here I can move splfit before the loop over mm_G but I have to changes paw_rhox_spl
       ilm_G=1+ll_G**2+ll_G+mm_G
       ider=0 ; wk_ffnl(:,:)=pwff_spl(:,:,ilm_G-1,klmn,itypat)  ! Note klmn and ilm_G-1
       call splfit(Psps%qgrid_ff,derfun,wk_ffnl,ider,qpg_norm,newfun,Psps%mqgrid_ff,npw)
       !
       ! === Evaluate matrix elements for each plane wave ===
       do ig=1,npw
        paw_rhox(1,ig,klmn,iatm) = paw_rhox(1,ig,klmn,iatm)+ &
&        newfun(ig)*ylm_q(ig,ilm_G)*(ph3d(1,ig)*mi_l(1,ipow)-ph3d(2,ig)*mi_l(2,ipow)) 

        paw_rhox(2,ig,klmn,iatm) = paw_rhox(2,ig,klmn,iatm)+ &
&        newfun(ig)*ylm_q(ig,ilm_G)*(ph3d(1,ig)*mi_l(2,ipow)+ph3d(2,ig)*mi_l(1,ipow)) 
       end do

      end do !ll_G
     end do !mm_G

    end do !ilmn
   end do !jlmn
  end do !iatm
  !
  ! * Multiply by 4\pi arising from the expansion of the plane wave
  paw_rhox=four_pi*paw_rhox

 CASE DEFAULT
  write(msg,'(4a,i3)')ch10,&
&  ' paw_mkrhox : BUG - ',ch10,&
&  '  Called with wrong value for method ',method
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT

 deallocate(ph3d,wk_ffnl)
 deallocate(newfun,derfun,qpg_norm)

end subroutine paw_mkrhox
!!***
