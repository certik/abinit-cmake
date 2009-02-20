!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_mkrhox_spl
!! NAME
!! paw_mkrhox_spl
!!
!! FUNCTION
!!  Evaluate PAW form factor ff^{aL}_{ij}(q) for each angular momentum L, 
!!  each type of atom, a, and each [(in,il),(jn,jl)] channel. These quantities
!!  are used in paw_mkrhox to evaluate $<phi|e^{-i(q+G)}|phj>-<tphi|e^{-i(q+G)}|tphj>$
!!  for an arbitrary q+G vector.
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dim1_rhox
!!   = 2*(Psps%mpsang-1)           if method=1
!!   = MAXVAL(Pawtab(:)%l_size)**2 if method=2
!!  dim2_rhox
!!   = Psps%lnmax*(Psps%lnmax+1)/2 if method=1
!!   = MAXVAL(Pawtab(:)%lmn2_size) if method=2
!!  method=integer flag defining the approach used:
!!   1 --> Expression based on the expansion on a plane wave in terms of Bessel functions 
!!        and spherical harmonics (Arnaud-Alouani's methos, see PRB 62, 4464 
!!   2 --> Approximate expression with correct description of the multipoles. Eq. 9 in PRB 74, 035101 
!!  Psps <type(pseudopotential_type)>=Info on pseudopotentials
!!      %mpsang=1+maximum angular momentum in PAW datasets
!!      %lnmax=Max number of (l,n) channels
!!      %mqgrid_ff=number of grid points in the q-mesh
!!      %qgrid_ff(mqgrid_ff)=values where form factors are returned
!!  ntypat=number of type of atoms
!!  Pawrad<type(Pawrad_type)>=datatype containing radial grid information
!!  Pawtab(ntypat)<type(pawtab_type)>=PAW tabulated starting data
!!
!! OUTPUT
!!  pwff_spl(mqgrid_ff,2,0:dim1_rhox,dim1_rhox2,ntypat)
!!   form factors ff^{aL}_{ij}(q) and second derivative in packed storage mode
!!  === if method=1 ===
!!    $ff_^{aL}_{ij}(q) = 
!!      \int_0^{r_a} j_L(2\pi qr) [phi_{n_i,l_i}.phi_{n_j l_j}(r) - tphi_{n_i l_i}.tph_{n_j l_j}(r)] dr$
!!  === if method=2 ===
!!    $ff_^{aL}_{ij}(q) = q_ij^{LM} \int_0^{r_a} j_L(2\pi qr) g_L(r) r^2 dr$
!!
!! NOTES
!!  * $j_L(2\pi q)$ is a spherical Bessel function
!!  * Output matrix elements are stored in packed storage mode
!!  * Inspired by psp7nl
!!
!! TODO 
!!  One might save CPU time taking into account Gaunt selection rules!
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

subroutine paw_mkrhox_spl(ntypat,Psps,Pawrad,Pawtab,pwff_spl,method,dim1_rhox,dim2_rhox)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit
 use m_special_funcs, only : jbessel_4spline


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ntypat,method,dim2_rhox,dim1_rhox 
 type(Pseudopotential_type), intent(in) :: Psps
!arrays
 real(dp),intent(out) :: pwff_spl(Psps%mqgrid_ff,2,0:dim1_rhox,dim2_rhox,ntypat)
 type(Pawrad_type),intent(in) :: Pawrad(ntypat)
 type(Pawtab_type),intent(in) :: Pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: lmn_size,lmn2_size,lmn2_size_max,mm,nlmn,jlmn,ilmn,klmn
 integer :: lnmax,iq,ir,il,ll,meshsz,mmax,iln,jln,itypat,nln,k0ln,kln,qlm,mqgrid_ff
 real(dp),parameter :: EPS=tol14**4,TOLJ=0.001_dp
 real(dp) :: arg,argn,bes,besp,bespp,qr,yp1,ypn
 logical :: ltest
 character(len=500) :: msg
 type(Pawrad_type) :: Tmpmesh
!arrays
 real(dp),allocatable :: rrshape_l(:),shape_l(:),ff(:),gg(:),rr(:),rrdphi_ij(:),work(:)
 real(dp),allocatable :: dphi_ij(:),tmp_spl(:,:,:,:),tmp_jgl(:,:,:)

!*************************************************************************

#if defined DEBUG_MODE
 write(msg,'(2a)')ch10,' paw_mkrhox_spl : enter'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

 mqgrid_ff = Psps%mqgrid_ff 
 pwff_spl(:,:,:,:,:)=zero

 SELECT CASE (method)
 CASE (1)  
  ! === Arnaud-Alouani exact expression PRB 62. 4464 ===
  ! * $ff_^{aL}_{ij}(q) = 
  !    \int_0^{r_a} j_L(2\pi qr) [phi_{n_i,l_i}.phi_{n_j l_j}(r)-tphi_{n_i l_i}.tph_{n_j l_j}(r)]dr$
  ! * It does not descrive correctly the multipoles of the AE charge density if low cutoff on G 
  write(msg,'(a)')' paw_mkrhox_spl: Using Arnaud-Alouani expression'
  call wrtout(std_out,msg,'COLL')
  lnmax = Psps%lnmax
  allocate(tmp_spl(mqgrid_ff,2,0:2*(Psps%mpsang-1),lnmax*(lnmax+1)/2))

  do itypat=1,ntypat 
   ! Is mesh beginning with r=0 ?
   ltest=(ABS(Pawrad(itypat)%rad(1))<tol10)
   call assert(ltest,'Radial mesh starts with r/=0',__FILE__,__LINE__)
   !
   ! === Initialize temporary arrays and variables ===
   call copymesh(Pawrad(itypat),Tmpmesh)
   meshsz=Tmpmesh%mesh_size ; mmax=meshsz
   allocate(dphi_ij(meshsz),rrdphi_ij(meshsz))
   allocate(ff(meshsz),gg(meshsz),rr(meshsz),work(mqgrid_ff))
   rr(:)=Tmpmesh%rad(:)
   !
   ! === Loop on (jln,iln) channels for this type. Packed form ===
   tmp_spl(:,:,:,:)=zero
   nln=Pawtab(itypat)%basis_size
   do jln=1,nln
    k0ln=jln*(jln-1)/2
    do iln=1,jln
     kln=k0ln+iln

     dphi_ij(:)  =Pawtab(itypat)%phiphj(:,kln)-Pawtab(itypat)%tphitphj(:,kln)
     rrdphi_ij(:)=rr(:)*dphi_ij(:)

     ir=meshsz
     do WHILE (ABS(dphi_ij(ir))<EPS)
      ir=ir-1
     end do
     ir=MIN(ir+1,meshsz)
     if (ir/=mmax) then
      mmax=ir 
      call compmesh(Tmpmesh,rr(mmax))
     end if
     !
     ! === Loop on l for Bessel function. Note the starting point ===
     ! TODO Here I should loop only the moments allowed by Gaunt, I should use indklm!
     ! and only on lmax for this atom
     do ll=0,2*(Psps%mpsang-1)
      !
      ! === Compute f_l(q=0) only if l=0, and first derivative fp_l(q=0) (nonzero only if ll==1) ===
      tmp_spl(1,1,ll,kln)=zero ; yp1=zero
      if (ll==0) call simp_gen(tmp_spl(1,1,ll,kln),dphi_ij,Tmpmesh)
      if (ll==1) then
       call simp_gen(yp1,rrdphi_ij,Tmpmesh)
       yp1=yp1*two_pi*third
      end if
      !
      ! === Compute f_l(0<q<qmax) ===
      if (mqgrid_ff>2) then
       do iq=2,mqgrid_ff-1
        arg=two_pi*Psps%qgrid_ff(iq)
        do ir=1,mmax
         qr=arg*rr(ir)
         call jbessel_4spline(bes,besp,ll,0,qr,TOLJ)
         ff(ir)=bes*dphi_ij(ir)
        end do
        call simp_gen(tmp_spl(iq,1,ll,kln),ff,Tmpmesh)
       end do
      end if
      !
      ! === Compute f_l(q=qmax) and first derivative ===
      if (mqgrid_ff>1) then
       argn=two_pi*Psps%qgrid_ff(mqgrid_ff)
       do ir=1,mmax
        qr=argn*rr(ir)
        call jbessel_4spline(bes,besp,ll,1,qr,TOLJ)
        ff(ir)=bes *  dphi_ij(ir) 
        gg(ir)=besp*rrdphi_ij(ir)
       end do
       call simp_gen(tmp_spl(mqgrid_ff,1,ll,kln),ff,Tmpmesh)
       gg(:)=two_pi*gg(:) !twopi comes from 2\pi|q|
       call simp_gen(ypn,gg,Tmpmesh)
      else
       ypn=yp1
      end if
      !
      ! === Compute second derivative of ff^{al}_{ij)(q) ===
      call spline(Psps%qgrid_ff,tmp_spl(:,1,ll,kln),mqgrid_ff,yp1,ypn,tmp_spl(:,2,ll,kln),work)
     end do !ll  

    end do !iln
   end do !jln
   !
   ! === Save values for this atom type, each ll and kln channel ===
   pwff_spl(:,:,:,:,itypat)=tmp_spl(:,:,:,:)

   deallocate(dphi_ij,rrdphi_ij,ff,gg,rr,work)
   deallocate(Tmpmesh%rad,Tmpmesh%radfact,Tmpmesh%simfact)
  end do !typat
  deallocate(tmp_spl)

 CASE (2)
  ! ==== Shishkin-Kresse approximated expression ====
  ! $ff_^{aL}_{ij}(q) = q_ij^{LM} \int_0^{r_a} j_L(2\pi qr) g_L(r) r^2 dr$
  ! * Better description of multipoles of AE charge, 
  ! * Better results for energy degeneracies in GW band structure 
  write(msg,'(a)')' paw_mkrhox_spl: Using Shishkin-Kresse expression'
  call wrtout(std_out,msg,'COLL')
  allocate(tmp_jgl(mqgrid_ff,2,0:2*(Psps%mpsang-1)))

  do itypat=1,ntypat 
   ! Is mesh beginning with r=0 ?
   ltest=(ABS(Pawrad(itypat)%rad(1))<tol10)
   call assert(ltest,'Radial mesh starts with r/=0',__FILE__,__LINE__)
   !
   ! === Initialize temporary arrays and variables ===
   nlmn=Pawtab(itypat)%lmn_size
   call copymesh(Pawrad(itypat),Tmpmesh)
   meshsz=Tmpmesh%mesh_size ; mmax=meshsz
   allocate(ff(meshsz),gg(meshsz),rr(meshsz))
   allocate(shape_l(meshsz),rrshape_l(meshsz),work(mqgrid_ff))
   rr(:)=Tmpmesh%rad(:)

   tmp_jgl(:,:,:)=zero
   rrshape_l(1)=zero
     shape_l(1)=zero
   !
   ! TODO Here I should loop only the moments allowed by Gaunt, I should use indklm!
   ! and only lmax for this atom
   do ll=0,2*(Psps%mpsang-1)

      shape_l(2:meshsz)=Pawtab(itypat)%shapefunc(2:meshsz,ll+1)*rr(2:meshsz)**2
    rrshape_l(2:meshsz)=Pawtab(itypat)%shapefunc(2:meshsz,ll+1)*rr(2:meshsz)**3
    !
    ! === Compute f_l(q=0) and first derivative fp_l(q=0) (only if ll==1) ===
    tmp_jgl(1,1,ll)=zero ; yp1=zero
    if (ll==0) call simp_gen(tmp_jgl(1,1,ll),shape_l,Tmpmesh)
    if (ll==1) then
     call simp_gen(yp1,rrshape_l,Tmpmesh) !rr comes from d/dq
     yp1=yp1*two_pi*third
    end if
    !
    ! === Compute f_l(0<q<qmax) ===
    if (mqgrid_ff>2) then
     do iq=2,mqgrid_ff-1
      arg=two_pi*Psps%qgrid_ff(iq)
      do ir=1,mmax
       qr=arg*rr(ir)
       call jbessel_4spline(bes,besp,ll,0,qr,TOLJ)
       ff(ir)=bes*shape_l(ir)
      end do
      call simp_gen(tmp_jgl(iq,1,ll),ff,Tmpmesh)
     end do
    end if
    !
    ! === Compute f_l(q=qmax) and first derivative ===
    if (mqgrid_ff>1) then
     argn=two_pi*Psps%qgrid_ff(mqgrid_ff)
     do ir=1,mmax
      qr=argn*rr(ir)
      call jbessel_4spline(bes,besp,ll,1,qr,TOLJ)
      ff(ir)=bes *  shape_l(ir) 
      gg(ir)=besp*rrshape_l(ir)
     end do
     call simp_gen(tmp_jgl(mqgrid_ff,1,ll),ff,Tmpmesh)
     gg(:)=two_pi*gg(:) !two_pi comes from 2\pi|q|
     call simp_gen(ypn,gg,Tmpmesh)
    else
     ypn=yp1
    end if
    !
    ! === Compute second derivative of ff_^{al}_{ij)(q) ===
    call spline(Psps%qgrid_ff,tmp_jgl(:,1,ll),mqgrid_ff,yp1,ypn,tmp_jgl(:,2,ll),work)
   end do !ll
   !
   ! === Save values for this type, each ll and ilmn,jlm channels ===
   ! * Here we assembly q_{ij}^{lm} \int_0^{r_a} j_l(2\pi(q+G)r) g_l(r)r^2 dr
   ! * Some of the contributions from qijl are zero due to Gaunt selection rules
   do jlmn=1,nlmn
    do ilmn=1,jlmn
     klmn=ilmn+(jlmn-1)*jlmn/2
     do ll=0,2*(Psps%mpsang-1)
      do mm=-ll,ll
       qlm=1+ll**2+ll+mm
       pwff_spl(:,:,qlm-1,klmn,itypat)=tmp_jgl(:,:,ll)*Pawtab(itypat)%qijl(qlm,klmn)
      end do
     end do
    end do
   end do

   deallocate(shape_l,rrshape_l,ff,gg,rr,work)
   deallocate(Tmpmesh%rad,Tmpmesh%radfact,Tmpmesh%simfact)
  end do !typat

  deallocate(tmp_jgl)

 CASE DEFAULT
  write(msg,'(4a,i3)')ch10,&
&  ' paw_mkrhox_spl: BUG - ',ch10,&
&  '  Called with wrong value for method ',method
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT

#if defined DEBUG_MODE
 write(msg,'(a)')' paw_mkrhox_spl : exit '
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

end subroutine paw_mkrhox_spl 
!!***
