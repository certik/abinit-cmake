!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdij
!! NAME
!! pawdij
!!
!! FUNCTION
!! Compute the different pseudopotential strengths Dij
!! of the PAW non local operator
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   ! atvshift(16,nsppol,natom)=potential energy shift for specific lm channel, spin and atom
!!  enunit=choice for units of output Dij
!!  natom=number of atoms in cell
!!  nfft=total number of FFt grid
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume
!!  vtrial(nfft,nspden)=GS potential (Hartree)
!!  vxc(nfft,nspden)=XC potential (Hartree) on the fine FFT mesh
!!
!! OUTPUT
!!  paw_ij(iatom)%dij(cplex_dij*lmn2_size,ndij)= total Dij terms
!!  May be complex if cplex_dij=2
!!        dij(:,:,1) contains Dij^up-up
!!        dij(:,:,2) contains Dij^dn-dn
!!        dij(:,:,3) contains Dij^up-dn (only if nspinor=2)
!!        dij(:,:,4) contains Dij^dn-up (only if nspinor=2)
!!
!! PARENTS
!!      respfn,scfcv
!!
!! CHILDREN
!!      leave_new,nderiv_gen,print_ij,simp_gen,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawdij(dtset,enunit,mpi_enreg,natom,nfft,ngfft,nspden,ntypat,&
&                 paw_an,paw_ij,pawang,pawfgrtab,pawprtvol,pawrad,&
&                 pawspnorb,pawtab,pawxcdev,typat,ucvol,vtrial,vxc)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13paw, except_this_one => pawdij
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: enunit,natom,nfft,nspden,ntypat,pawprtvol,pawspnorb,pawxcdev
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(in) :: vtrial(nfft,nspden),vxc(nfft,nspden)
 type(paw_an_type),intent(in) :: paw_an(natom)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawtab_type),intent(in) :: pawtab(ntypat)
 type(dataset_type),intent(in) :: dtset

!Local variables ---------------------------------------
!scalars
 integer :: cplex_dij
 integer :: has_dijhat,has_dijU,has_dijso,has_dijxc_val,has_dijxc
 integer :: iatom,ic,icount,idij,ier,ij_size,iklnu,iklnuold
 integer :: ilm,iln,ils,ilslm,im1,im2,in1,in2,indk,ipts
 integer :: isel,ispden,itypat,j0lm,j0ln,jlm,jln,klm,klm1,klmn,klmn1,klmnu
 integer :: klmnuold,kln,l_size,lexexch,ll,lm0,lm_size,lmax,lmin
 integer :: lmn2_size,lmn_size,lpawu,mesh_size,mm,nfftot,npts,nsploop,nsppol
 integer :: old_paral_level,spaceComm
 real(dp), parameter :: HalfFineStruct2=half/InvFineStruct**2
 real(dp) :: coeffpawu,fact,VUKS,vxcij,vxcij22,vxcijhat,tmp
 character(len=500) :: message
!arrays
 integer,allocatable :: idum(:),indklmn(:,:)
 real(dp), parameter :: facti(2)=(/1.d0,-1.d0/)
 real(dp) :: rdum(1),tsec(2)
 real(dp),allocatable :: dijexxc(:),dijhat(:),dijhat_tmp(:),dijpawu(:),dijso(:),dijsym(:)
 real(dp),allocatable :: dijso_rad(:),dijxc(:),dijxc_hat(:),dijxchat_tmp(:),dijxc_val(:)
 real(dp),allocatable :: dv1dr(:),ff(:),prod(:),prodxchat(:),vpawu(:,:)
 real(dp),allocatable :: vxcij1(:),vxcij2(:),vxcijtot(:),vxcij1_val(:),vxcval_ij(:)
 real(dp),allocatable :: yylmr(:,:)
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)

! *************************************************************************

#if defined DEBUG_MODE
 write(message,'(a)')' pawdij : enter '
 call wrtout(std_out,message,'COLL') 
 call flush_unit(std_out)
#endif

 call timab(561,1,tsec)

 if (pawspnorb>0.and.paw_ij(1)%ndij/=4) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawdij : BUG -',ch10,&
&  '  invalid ndij size for Dij with spin-orbit coupling !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if (pawspnorb>0.and.paw_ij(1)%cplex_dij/=2) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawdij : BUG -',ch10,&
&  '  invalid cplex size for Dij with spin-orbit coupling !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if (paw_ij(1)%ndij==4.and.paw_ij(1)%cplex_dij/=2) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawdij : BUG -',ch10,&
&  '  invalid cplex size for Dij (4 Dij components) !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(paw_ij(1)%has_dijhartree==0) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawdij : BUG -',ch10,&
&  '  dijhartree must be allocated !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 if (abs(pawprtvol)>=1) then
  write(message, '(2a)')ch10,' ==== In pawdij: several values of Dij (Hartree) ============'
  call wrtout(6,message,'COLL')
 end if
 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 nsppol=paw_ij(1)%nsppol
 nsploop=nsppol;if (paw_ij(1)%ndij==4) nsploop=4
 VUKS=zero

!----- Preliminary computation of Ylm,Ylpmp (only if angular mesh)
 if (pawxcdev==0) then
  npts=pawang%angl_size
  allocate(yylmr(pawang%l_max**2*(pawang%l_max**2+1)/2,npts))
  do ipts=1,npts
   do jlm=1,pawang%l_max**2
    j0lm=jlm*(jlm-1)/2
    do ilm=1,jlm
     klm=j0lm+ilm
     yylmr(klm,ipts)=pawang%ylmr(ilm,ipts)*pawang%ylmr(jlm,ipts)
    end do
   end do
  end do
 end if

!------------------------------------------------------------------------
!----- Big loop over atoms
!------------------------------------------------------------------------

 do iatom=1,natom

! ------------------------------------------------------------------------
! ----------- Allocations and initializations
! ------------------------------------------------------------------------

  itypat=typat(iatom)
  l_size=pawtab(itypat)%l_size
  mesh_size=pawrad(itypat)%mesh_size
  lmn_size=paw_ij(iatom)%lmn_size
  lmn2_size=paw_ij(iatom)%lmn2_size
  lm_size=paw_an(iatom)%lm_size
  ij_size=pawtab(itypat)%ij_size
  cplex_dij=paw_ij(iatom)%cplex_dij
  paw_ij(iatom)%dij(:,:)=zero

  has_dijxc_val=0
  if (paw_an(iatom)%has_vxcval>0.and.paw_ij(iatom)%has_dijxc_val>0) then
   has_dijxc_val=1
   paw_ij(iatom)%dijxc_val(:,:)=zero
   allocate(dijxc_val(cplex_dij*lmn2_size))
  end if

  has_dijxc=0
  if (paw_ij(iatom)%has_dijxc>0) then
   has_dijxc=1 
   allocate(dijxc_hat(cplex_dij*lmn2_size))
   paw_ij(iatom)%dijxc(:,:)=zero
  end if

  has_dijso=0 
  if (pawspnorb>0.and.paw_ij(iatom)%has_dijso>0) then
   has_dijso=1
   paw_ij(iatom)%dijso(:,:)=zero
  end if

  has_dijU=0
  if (pawtab(itypat)%usepawu>0.and.paw_ij(iatom)%has_dijU>0) then 
   has_dijU=1
   paw_ij(iatom)%dijU(:,:)=zero
  end if

  allocate(dijhat(cplex_dij*lmn2_size),dijxc(cplex_dij*lmn2_size))
  if (pawspnorb>0) allocate(dijso(cplex_dij*lmn2_size))
  if (pawtab(itypat)%usepawu>0) allocate(dijpawu(cplex_dij*lmn2_size))
  if (pawtab(itypat)%useexexch>0) allocate(dijexxc(cplex_dij*lmn2_size))

  has_dijhat=0
  if (paw_ij(iatom)%has_dijhat>0) then
   has_dijhat=1 
   paw_ij(iatom)%dijhat(:,:)=zero
  end if

  allocate(ff(mesh_size))
  allocate(indklmn(6,lmn2_size))
  indklmn(:,:)=pawtab(itypat)%indklmn(:,:)

! Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
  if (pawfgrtab(iatom)%gylm_allocated==0) then
   if (associated(pawfgrtab(iatom)%gylm)) deallocate(pawfgrtab(iatom)%gylm)
   allocate(pawfgrtab(iatom)%gylm(pawfgrtab(iatom)%nfgd,lm_size))
   pawfgrtab(iatom)%gylm_allocated=2
   call pawgylm(pawfgrtab(iatom)%gylm,rdum,rdum,iatom,pawfgrtab(iatom)%ifftsph,&
&   itypat,lm_size,pawfgrtab(iatom)%nfgd,1,0,0,pawtab(itypat),&
&   pawfgrtab(iatom)%rfgd,pawfgrtab(iatom)%rfgd_allocated)
  end if

! Eventually compute <Phi_i|1/r.dV/dr|Phi_j>*alpha2/2*Y_00 (for spin-orbit)
  if (pawspnorb>0) then
   allocate(dv1dr(mesh_size),dijso_rad(ij_size))
   fact=one/sqrt(four_pi) ! Y_00
   if (pawxcdev/=0) then
    if (nspden==1) then
     ff(1:mesh_size)=paw_an(iatom)%vxc1(1:mesh_size,1,1)
    else
     ff(1:mesh_size)=half*(paw_an(iatom)%vxc1(1:mesh_size,1,1) &
&     +paw_an(iatom)%vxc1(1:mesh_size,1,2))
    end if
   else
    ff(1:mesh_size)=zero
    if (nspden==1) then
     do ipts=1,npts
      ff(1:mesh_size)=ff(1:mesh_size) &
&      +paw_an(iatom)%vxc1(1:mesh_size,ipts,1) &
&      *pawang%angwgth(ipts)
     end do
    else
     do ipts=1,npts
      ff(1:mesh_size)=ff(1:mesh_size) &
&      +half*(paw_an(iatom)%vxc1(1:mesh_size,ipts,1) &
&      +paw_an(iatom)%vxc1(1:mesh_size,ipts,2)) &
&      *pawang%angwgth(ipts)
     end do
    end if
    ff(1:mesh_size)=sqrt(four_pi)*ff(1:mesh_size)
   end if
   ff(1:mesh_size)=fact*(ff(1:mesh_size)+paw_an(iatom)%vh1(1:mesh_size,1,1))
   call nderiv_gen(dv1dr,ff,1,pawrad(itypat))
   dv1dr(2:mesh_size)=HalfFineStruct2*(one/(one-ff(2:mesh_size)/InvFineStruct**2)) &
&   *dv1dr(2:mesh_size)/pawrad(itypat)%rad(2:mesh_size)
   call deducer0(dv1dr,mesh_size,pawrad(itypat))
   do kln=1,ij_size
    ff(1:mesh_size)= dv1dr(1:mesh_size)*pawtab(itypat)%phiphj(1:mesh_size,kln)
    call simp_gen(dijso_rad(kln),ff,pawrad(itypat))
   end do
   deallocate(dv1dr)
  end if

! ------------------------------------------------------------------------
! ----- Loop over density components
! ------------------------------------------------------------------------
  do idij=1,nsploop

!  Print title
   if (abs(pawprtvol)>=1) then
    if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
     if (nspden==2.and.nsppol==1) then
      write(message, '(2a,i3,3a)') ch10,&
&      ' >>>>>>>>>> Atom ',iatom,':',ch10,&
&      ' (antiferromagnetism case: only one spin component)'
     else
      write(message, '(2a,i3,3a)') ch10,&
&      ' >>>>>>>>>> Atom ',iatom,' (component ',trim(dspin(idij+2*(nsploop/4))),'):'
     end if
     call wrtout(6,message,'COLL')
    end if
   end if

!  ------------------------------------------------------------------------
!  ----------- Load atomic Dij0 into Dij
!  ------------------------------------------------------------------------

   if (idij<=2) then

    klmn1=1
    do klmn=1,lmn2_size
     paw_ij(iatom)%dij(klmn1,idij)=pawtab(itypat)%dij0(klmn)
     klmn1=klmn1+cplex_dij
    end do

    if (abs(pawprtvol)>=1) then
     if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
      write(message, '(a)') '   ************ Dij atomic (Dij0) ***********'
      call wrtout(6,message,'COLL')
      call print_ij(pawtab(itypat)%dij0,lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
     end if
    end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_Hartree to Dij
!  ------------------------------------------------------------------------

   if (idij<=2) then

    klmn1=1
    do klmn=1,lmn2_size
     paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+paw_ij(iatom)%dijhartree(klmn)
     klmn1=klmn1+cplex_dij
    end do

    if (abs(pawprtvol)>=1) then
     if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
      write(message, '(a)')'   ************** Dij Hartree ***************'
      call wrtout(6,message,'COLL')
      call print_ij(paw_ij(iatom)%dijhartree,lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
     end if
    end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_xc to Dij
!  ------------------------------------------------------------------------

   if (idij<=nsppol.or.(nspden==4.and.idij<=3)) then
    allocate(vxcijtot(lmn2_size))
    if (has_dijxc_val==1) allocate(vxcval_ij(lmn2_size)) 

    do ispden=idij,idij+idij/3

     vxcijtot=zero
     if (has_dijxc_val==1) vxcval_ij=zero

!    ================================================
!    ===== First formalism: use (l,m) moments for vxc
!    ================================================
     if (pawxcdev/=0) then

      allocate(vxcij1(ij_size))
      if (has_dijxc_val==1) allocate(vxcij1_val(ij_size))

      do klm=1,lm_size
!      Summing over klm moments.
       if (paw_an(iatom)%lmselect(klm)) then

!       ===== Vxc_ij_1 (tmp) =====
        vxcij1=zero
        do kln=1,ij_size
         ff(1:mesh_size)= paw_an(iatom)%vxc1(1:mesh_size,klm,ispden)&
&         *pawtab(itypat)%phiphj(1:mesh_size,kln)&
&         - paw_an(iatom)%vxct1(1:mesh_size,klm,ispden)&
&         *pawtab(itypat)%tphitphj(1:mesh_size,kln)
         call simp_gen(vxcij1(kln),ff,pawrad(itypat))
        end do

!       ===== Vxc_ij_2 (tmp) =====
        vxcij22=zero
        ll=1+int(sqrt(dble(klm)-0.1))
        ff(1:mesh_size)=paw_an(iatom)%vxct1(1:mesh_size,klm,ispden)&
&        *pawtab(itypat)%shapefunc(1:mesh_size,ll)&
&        *pawrad(itypat)%rad(1:mesh_size)**2
        call simp_gen(vxcij22,ff,pawrad(itypat))

!       ==== If required calculate valence-only onsite matrix elements ====
!       $\int dr phiphj.vxc(n1) - tphitphj.vxc(tn1+nhat)$ 
        if (has_dijxc_val==1) then 
         vxcij1_val(:)=zero
         do kln=1,ij_size
          ff(1:mesh_size)= &
&          paw_an(iatom)%vxc1_val (1:mesh_size,klm,ispden)*pawtab(itypat)%phiphj  (1:mesh_size,kln)&
&          -paw_an(iatom)%vxct1_val(1:mesh_size,klm,ispden)*pawtab(itypat)%tphitphj(1:mesh_size,kln)
          call simp_gen(vxcij1_val(kln),ff,Pawrad(itypat))
         end do
        end if

!       ===== Vxc_ij and Vxcij_hat from Vxc_ij_1 and Vxc_ij_2 =====
!       Looping over partial waves while accumulating over klm moments
        do klmn=1,lmn2_size
         klm1=indklmn(1,klmn);kln=indklmn(2,klmn)
         vxcij=zero
         isel=pawang%gntselect(klm,klm1)
         if (isel>0) vxcij=vxcij1(kln)*pawang%realgnt(isel)
         vxcijhat=pawtab(itypat)%qijl(klm,klmn)*vxcij22

!        ===== Contribution to total Vxc_ij =====
         vxcijtot(klmn)=vxcijtot(klmn)+vxcij-vxcijhat

         if (has_dijxc_val==1) then
          if (isel>0) vxcval_ij(klmn)=vxcval_ij(klmn) + vxcij1_val(kln)*pawang%realgnt(isel)
         end if

        end do ! Loop klmn
       end if
      end do  ! Loop klm

      deallocate(vxcij1)
      if (has_dijxc_val==1) deallocate(vxcij1_val)

!     ================================================
!     ===== Second formalism: use vxc on r,theta,phi
!     ================================================
     else

      allocate(vxcij1(ij_size),vxcij2(l_size))
      if (has_dijxc_val==1) allocate(vxcij1_val(ij_size))

!     ===== Loop on angular mesh =====
      do ipts=1,npts

!      ===== Vxc_ij_1 (tmp) =====
       vxcij1=zero
       do kln=1,ij_size
        ff(1:mesh_size)= paw_an(iatom)%vxc1(1:mesh_size,ipts,ispden)&
&        *pawtab(itypat)%phiphj(1:mesh_size,kln)&
&        - paw_an(iatom)%vxct1(1:mesh_size,ipts,ispden)&
&        *pawtab(itypat)%tphitphj(1:mesh_size,kln)
        call simp_gen(vxcij1(kln),ff,pawrad(itypat))
       end do

!      ==== If required calculate valence-only matrix elements ====
       if (has_dijxc_val==1) then
        vxcij1_val(:)=zero
        do kln=1,ij_size
         ff(1:mesh_size)= &
&         paw_an(iatom)%vxc1_val (1:mesh_size,ipts,ispden)*pawtab(itypat)%phiphj  (1:mesh_size,kln)&
&         -paw_an(iatom)%vxct1_val(1:mesh_size,ipts,ispden)*pawtab(itypat)%tphitphj(1:mesh_size,kln)
         call simp_gen(vxcij1_val(kln),ff,Pawrad(itypat))
        end do
       end if

!      ===== Vxc_ij_2 (tmp) =====
       vxcij2=zero
       do ils=1,l_size
        ff(1:mesh_size)=paw_an(iatom)%vxct1(1:mesh_size,ipts,ispden)&
&        *pawtab(itypat)%shapefunc(1:mesh_size,ils)&
&        *pawrad(itypat)%rad(1:mesh_size)**2
        call simp_gen(vxcij2(ils),ff,pawrad(itypat))
       end do

!      ===== Vxc_ij and Vxcij_hat from Vxc_ij_1 and Vxc_ij_2 =====
!      Build elements integrating over the angular mesh
       do klmn=1,lmn2_size
        klm=indklmn(1,klmn);kln=indklmn(2,klmn)
        lmin=indklmn(3,klmn);lmax=indklmn(4,klmn)

        vxcij=vxcij1(kln)*pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi

        vxcijhat=zero
        do ils=lmin,lmax,2
         lm0=ils**2+ils+1
         do mm=-ils,ils
          ilslm=lm0+mm;isel=pawang%gntselect(ilslm,klm)
          if (isel>0) vxcijhat=vxcijhat+vxcij2(ils+1)*pawang%angwgth(ipts)&
&          *pawtab(itypat)%qijl(ilslm,klmn)*four_pi&
&          *pawang%ylmr(ilslm,ipts)
         end do
        end do

!       ===== Contribution to total Vxc_ij =====
        vxcijtot(klmn)=vxcijtot(klmn)+vxcij-vxcijhat

        if (has_dijxc_val==1) then
         tmp=vxcij1_val(kln)*pawang%angwgth(ipts)*yylmr(klm,ipts)*four_pi
         vxcval_ij(klmn)=vxcval_ij(klmn)+tmp
        end if
       end do ! Loop klmn
      end do  ! Loop ipts

      deallocate(vxcij1,vxcij2)
      if (has_dijxc_val==1) deallocate(vxcij1_val)

     end if  ! choice XC

     if (ispden<3) then
      dijxc(1:lmn2_size)=vxcijtot(1:lmn2_size)
      if (has_dijxc_val==1) dijxc_val(1:lmn2_size)=vxcval_ij(1:lmn2_size)
     else
      klmn1=max(1,ispden-2)
      do klmn=1,lmn2_size
       dijxc(klmn1)=vxcijtot(klmn)
       klmn1=klmn1+cplex_dij
      end do
      if (has_dijxc_val==1)  then
       klmn1=max(1,ispden-2)
       do klmn=1,lmn2_size
        dijxc_val(klmn1)=vxcval_ij(klmn)
        klmn1=klmn1+cplex_dij
       end do
      end if
     end if

    end do
    deallocate(vxcijtot)
    if (has_dijxc_val==1) deallocate(vxcval_ij)

   else if (nspden==4.and.idij==4) then
    klmn1=2
    do klmn=1,lmn2_size
     dijxc(klmn1)=-dijxc(klmn1)
     klmn1=klmn1+cplex_dij
    end do
    if (has_dijxc_val==1) then
     klmn1=2
     do klmn=1,lmn2_size
      dijxc_val(klmn1)=-dijxc_val(klmn1)
      klmn1=klmn1+cplex_dij
     end do
    end if
   end if

   if (idij<=nsppol.or.idij==2) then
    klmn1=1
    do klmn=1,lmn2_size
     paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijxc(klmn)
     klmn1=klmn1+cplex_dij
    end do
    if (has_dijxc==1) then
     klmn1=1
     do klmn=1,lmn2_size
      paw_ij(iatom)%dijxc(klmn1,idij)=dijxc    (klmn)
      klmn1=klmn1+cplex_dij
     end do
    end if
    if (has_dijxc_val==1) then
     klmn1=1
     do klmn=1,lmn2_size
      paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijxc(klmn)
      klmn1=klmn1+cplex_dij
     end do
    end if
   else if (nspden==4) then
    paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijxc(:)
    if (has_dijxc    ==1) paw_ij(iatom)%dijxc    (:,idij)=dijxc(:)
    if (has_dijxc_val==1) paw_ij(iatom)%dijxc_val(:,idij)=dijxc_val(:)
   end if

   if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
    if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
     write(message, '(a)')'   ****************** Dij_xc ****************'
     call wrtout(6,message,'COLL')
     if (idij<=nsppol.or.idij==2) then
      call print_ij(dijxc(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
     else
      call print_ij(dijxc,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=1)
     end if
    end if
   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_hat to Dij
!  ------------------------------------------------------------------------

   if (idij<=nsppol.or.(nspden==4.and.idij<=3)) then

    allocate(dijhat_tmp(lmn2_size))
    if (has_dijxc==1) allocate(dijxchat_tmp(lmn2_size))

    do ispden=idij,idij+idij/3

     dijhat_tmp=zero
     allocate(prod(lm_size));prod=zero

     do ilslm=1,lm_size
      do ic=1,pawfgrtab(iatom)%nfgd
       prod(ilslm)=prod(ilslm)+vtrial(pawfgrtab(iatom)%ifftsph(ic),ispden)&
&       *pawfgrtab(iatom)%gylm(ic,ilslm)
      end do
     end do
     if(mpi_enreg%paral_fft==1)then
      old_paral_level= mpi_enreg%paral_level
      mpi_enreg%paral_level=3
      call xcomm_init(mpi_enreg,spaceComm)
      if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
      call xsum_mpi(prod,spaceComm,ier)
      mpi_enreg%paral_level=old_paral_level
     end if

     if (has_dijxc==1) then 
      dijxchat_tmp(:)=zero
!     === Evaluate prod i.e $\sum_{lm} \int g_l Ylm v_xc[tn+nhat+tnc]dr$ on the FFT mesh ===
!     * It does not depend on ij
      allocate(prodxchat(lm_size)) ; prodxchat(:)=zero
      do ilslm=1,lm_size
       do ic=1,pawfgrtab(iatom)%nfgd
        prodxchat(ilslm)=prodxchat(ilslm)&
&        +vxc(pawfgrtab(iatom)%ifftsph(ic),ispden)*pawfgrtab(iatom)%gylm(ic,ilslm)
       end do
      end do
      if (mpi_enreg%paral_fft==1)then
       old_paral_level= mpi_enreg%paral_level
       mpi_enreg%paral_level=3
       call xcomm_init(mpi_enreg,spaceComm)
       if (mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft
       call xsum_mpi(prodxchat,spaceComm,ier)
       mpi_enreg%paral_level=old_paral_level
      end if
     end if

     do klmn=1,lmn2_size
      klm=indklmn(1,klmn)
      lmin=indklmn(3,klmn);lmax=indklmn(4,klmn)
      do ils=lmin,lmax,2
       lm0=ils**2+ils+1
       do mm=-ils,ils
        ilslm=lm0+mm;isel=pawang%gntselect(lm0+mm,klm)
        if (isel>0) dijhat_tmp(klmn)=dijhat_tmp(klmn)+prod(ilslm)*pawtab(itypat)%qijl(ilslm,klmn)
        if (has_dijxc==1.and.isel>0) dijxchat_tmp(klmn)=dijxchat_tmp(klmn)+prodxchat(ilslm)*pawtab(itypat)%qijl(ilslm,klmn)
       end do
      end do
     end do
     deallocate(prod)
     if (has_dijxc==1) deallocate(prodxchat)

     if (ispden<3) then
      dijhat(1:lmn2_size)=dijhat_tmp(1:lmn2_size)*ucvol/dble(nfftot)
      if (has_dijxc==1) dijxc_hat(1:lmn2_size)=dijxchat_tmp(1:lmn2_size)*ucvol/dble(nfftot)
     else
      klmn1=max(1,ispden-2)
      do klmn=1,lmn2_size
       dijhat(klmn1)=dijhat_tmp(klmn)*ucvol/dble(nfftot)
       klmn1=klmn1+cplex_dij
      end do
      if (has_dijxc==1) then
       klmn1=max(1,ispden-2)
       do klmn=1,lmn2_size
        dijxc_hat(klmn1)=dijxchat_tmp(klmn)*ucvol/DBLE(nfftot)
        klmn1=klmn1+cplex_dij
       end do
      end if
     end if

    end do
    deallocate(dijhat_tmp)
    if (has_dijxc==1) deallocate(dijxchat_tmp)

   else if (nspden==4.and.idij==4) then
    klmn1=2
    do klmn=1,lmn2_size
     dijhat(klmn1)=-dijhat(klmn1)
     klmn1=klmn1+cplex_dij
    end do
    if (has_dijxc==1) then
     klmn1=2
     do klmn=1,lmn2_size
      dijxc_hat(klmn1)=-dijxc_hat(klmn1)
      klmn1=klmn1+cplex_dij
     end do
    end if
   end if

   if (idij<=nsppol.or.idij==2) then
    klmn1=1
    do klmn=1,lmn2_size
     paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijhat(klmn)
     klmn1=klmn1+cplex_dij
    end do
    if (has_dijhat>0) then
     klmn1=1
     do klmn=1,lmn2_size
      paw_ij(iatom)%dijhat(klmn1,idij)=dijhat(klmn)
      klmn1=klmn1+cplex_dij
     end do
    end if
    if (has_dijxc==1) then
     klmn1=1
     do klmn=1,lmn2_size
      paw_ij(iatom)%dijxc(klmn1,idij)=paw_ij(iatom)%dijxc(klmn1,idij)+dijxc_hat(klmn)
      klmn1=klmn1+cplex_dij
     end do
    end if
   else if (nspden==4) then
    paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijhat(:)
    if (has_dijhat>0) paw_ij(iatom)%dijhat(:,idij)=dijhat(:)
    if (has_dijxc==1) paw_ij(iatom)%dijxc(:,idij)=paw_ij(iatom)%dijxc(:,idij)+dijxc_hat(:)
   end if

   if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
    if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
     write(message, '(a)')'   ************* Dij_hat (Veff_ij) **********'
     call wrtout(6,message,'COLL')
     if (idij<=nsppol.or.idij==2) then
      call print_ij(dijhat(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
     else
      call print_ij(dijhat,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=1)
     end if
    end if
   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij spin-orbit to Dij
!  ------------------------------------------------------------------------
   if (pawspnorb>0) then

    klmn1=1
    if (mod(idij,2)==1) then
     ispden=(1+idij)/2
     dijso(:)=zero
     do klmn=1,lmn2_size
      if (indklmn(3,klmn)==0) then   ! il==jl
       klm=indklmn(1,klmn);kln=indklmn(2,klmn)
       ilm=indklmn(5,klmn);jlm=indklmn(6,klmn)
       fact=dijso_rad(kln);if (ilm>jlm) fact=-fact
       dijso(klmn1  )=fact*pawang%ls_ylm(1,klm,ispden)
       dijso(klmn1+1)=fact*pawang%ls_ylm(2,klm,ispden)
      end if
      klmn1=klmn1+cplex_dij
     end do
    else if (idij==2) then
     do klmn=1,lmn2_size
      if (indklmn(3,klmn)==0) then   ! il==jl
       dijso(klmn1:klmn1+1)=-dijso(klmn1:klmn1+1)
      end if
      klmn1=klmn1+cplex_dij
     end do
    else if (idij==4) then
     do klmn=1,lmn2_size
      if (indklmn(3,klmn)==0) then   ! il==jl
       dijso(klmn1)=-dijso(klmn1)
      end if
      klmn1=klmn1+cplex_dij
     end do
    end if

    paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijso(:)
    if (has_dijso==1) paw_ij(iatom)%dijso(:,idij)=dijso(:)

!   Printing
    if (abs(pawprtvol)>=1) then
     if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
      write(message, '(a)')'   ************** Dij SpinOrbit ************'
      call wrtout(6,message,'COLL')
      call print_ij(dijso,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=3)
     end if
    end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_{LDA+U} to Dij
!  --------------------------------------------------------------------------
!  Dijpawu^{\sigma}_{mi,ni,mj,nj}=
!  \sum_{m,m'} [vpawu^{\sigma}_{m,m'}*phiphjint_{ni,nj}^{m,m'}]=
!  [vpawu^{\sigma}_{mi,mj}*phiphjint_{ni,nj}]
!  --------------------------------------------------------------------------
   if (pawtab(itypat)%usepawu>0) then

    if (idij<=nsppol.or.(nspden==4.and.idij<=3)) then
     do ispden=idij,idij+idij/3

      if(abs(pawprtvol)>=3) then
       write(message,*) "pawdij, LDA+U calculation for iatom,ispden= ",iatom,ispden,itypat
       call wrtout(06,message,'COLL')
      end if

      lpawu=pawtab(itypat)%lpawu

      allocate(vpawu(lpawu*2+1,lpawu*2+1))
      call pawpupot(ispden,nspden,paw_ij(iatom),pawprtvol,pawtab(itypat),vpawu,VUKS)

      dijpawu=zero
      klmn1=max(1,ispden-2)
      indk=1;if (ispden>=3) indk=cplex_dij
      do klmn=1,lmn2_size
       im1=pawtab(itypat)%klmntomn(1,klmn)
       im2=pawtab(itypat)%klmntomn(2,klmn)
       in1=pawtab(itypat)%klmntomn(3,klmn)
       in2=pawtab(itypat)%klmntomn(4,klmn)
       lmin=pawtab(itypat)%indklmn(3,klmn)
       lmax=pawtab(itypat)%indklmn(4,klmn)
       if(lmin==0.and.lmax==2*lpawu) then
        icount=in1+(in2*(in2-1))/2
        if(pawtab(itypat)%ij_proj<icount)  then
         write(message, '(4a)' ) ch10,&
&         ' pawdij : BUG -',ch10,&
&         '  PAW+U: Problem in the loop for calculating dijpawu'
         call wrtout(06,message,'COLL')
         call leave_new('COLL')
        end if
        coeffpawu=vpawu(im1,im2)
        if(dtset%natvshift/=0)then
         if(im1==im2)then
          coeffpawu=coeffpawu+dtset%atvshift(im1,idij,iatom)
         end if
        end if
        dijpawu(klmn1)=pawtab(itypat)%phiphjint(icount)*coeffpawu
       end if
       klmn1=klmn1+indk
      end do
      deallocate(vpawu)

     end do ! ispden

    else if (nspden==4.and.idij==4) then
     klmn1=2
     do klmn=1,lmn2_size
      dijpawu(klmn1)=-dijpawu(klmn1)
      klmn1=klmn1+cplex_dij
     end do
    end if

    if (idij<=nsppol.or.idij==2) then
     klmn1=1
     do klmn=1,lmn2_size
      paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijpawu(klmn)
      klmn1=klmn1+cplex_dij
     end do
     if (has_dijU==1) then
      klmn1=1
      do klmn=1,lmn2_size
       paw_ij(iatom)%dijU(klmn1,idij)=dijpawu(klmn)
       klmn1=klmn1+cplex_dij
      end do
     end if     
    else if (nspden==4) then
     paw_ij(iatom)%dij(:,idij)=paw_ij(iatom)%dij(:,idij)+dijpawu(:)
     if (has_dijU==1) paw_ij(iatom)%dijU(:,idij)=dijpawu(:)
    end if

    if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
     if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
      write(message, '(a)')'   ************* Dij_LDA+U (dijpawu) **********'
      call wrtout(6,message,'COLL')
      if (idij<=nsppol.or.idij==2) then
       call print_ij(dijpawu(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
      else
       call print_ij(dijpawu,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=1)
      end if
     end if
    end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Add Dij_{local exact-exchange} to Dij
!  --------------------------------------------------------------------------

   if (pawtab(itypat)%useexexch>0) then

    if(nspden==4)  then
     write(message, '(4a)' ) ch10,&
&     '  pawdenpot : ERROR -',ch10,&
&     '  Local exact-exch. not implemented for nspden=4 !'
     call wrtout(ab_out,message,'COLL')
     call wrtout(06,  message,'COLL')
     call leave_new('COLL')
    end if

    if (idij<=nsppol) then

     if(pawprtvol>=3) then
      write(message,*) "pawdij, local exact-exchange calculation for iatom,idij",iatom,idij
      call wrtout(06,  message,'COLL')
     end if

     lexexch=pawtab(itypat)%lexexch

     allocate(vxcij1(ij_size),vxcijtot(lmn2_size));vxcijtot=zero
     do klm=1,lm_size
      if (paw_an(iatom)%lmselect(klm)) then
!      ===== Vxc_ij_1 (tmp) =====
       vxcij1=zero
       do jln=pawtab(itypat)%lnproju(1),pawtab(itypat)%lnproju(pawtab(itypat)%nproju)
        j0ln=jln*(jln-1)/2
        do iln=pawtab(itypat)%lnproju(1),jln
         kln=j0ln+iln
         ff(1:mesh_size)=paw_an(iatom)%vxc_ex(1:mesh_size,klm,idij) &
&         *pawtab(itypat)%phiphj(1:mesh_size,kln)
         call simp_gen(vxcij1(kln),ff,pawrad(itypat))
        end do
       end do
!      ===== Contribution to total Vxc_ij =====
       do klmn=1,lmn2_size
        lmin=pawtab(itypat)%indklmn(3,klmn)
        lmax=pawtab(itypat)%indklmn(4,klmn)
        if(lmin==0.and.lmax==2*lexexch) then
         klm1=indklmn(1,klmn);kln=indklmn(2,klmn)
         isel=pawang%gntselect(klm,klm1)
         if (isel>0) vxcijtot(klmn)=vxcijtot(klmn)+vxcij1(kln)*pawang%realgnt(isel)
        end if
       end do ! Loop klmn
      end if
     end do  ! Loop klm
     deallocate(vxcij1)

     if (abs(pawprtvol)>=4) then
      if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
       write(message, '(a)')'   ** INFO ***** Dij_local_exact-exchange (vxcij) **********'
       call wrtout(6,message,'COLL')
       call print_ij(vxcijtot,lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
      end if
     end if

     dijexxc=zero
     do klmn=1,lmn2_size
      in1=pawtab(itypat)%klmntomn(3,klmn)
      in2=pawtab(itypat)%klmntomn(4,klmn)
      lmin=pawtab(itypat)%indklmn(3,klmn)
      lmax=pawtab(itypat)%indklmn(4,klmn)
      if(lmin==0.and.lmax==2*lexexch) then
       icount=in1+(in2*(in2-1))/2
       if(pawtab(itypat)%ij_proj<icount)  then
        write(message, '(4a)' ) ch10,&
&        ' pawdij : BUG -',ch10,&
&        '  PAW exact-exchange: Problem in the loop for calculating dijexxc'
        call wrtout(06,message,'COLL')
        call leave_new('COLL')
       end if
       dijexxc(klmn)=(paw_ij(iatom)%vpawx(1,klmn,idij)-vxcijtot(klmn))*pawtab(itypat)%exchmix
      end if
     end do
     deallocate(vxcijtot)

    end if

    if (idij<=nsppol.or.idij==2) then
     klmn1=1
     do klmn=1,lmn2_size
      paw_ij(iatom)%dij(klmn1,idij)=paw_ij(iatom)%dij(klmn1,idij)+dijexxc(klmn)
      klmn1=klmn1+cplex_dij
     end do
    end if

    if ((abs(pawprtvol)>=1).and.(idij<=2.or.nspden==4)) then
     if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
      write(message, '(a)')'   ************* Dij_Exact exchange (dijexxc) **********'
      call wrtout(6,message,'COLL')
      if (idij<=nsppol.or.idij==2) then
       call print_ij(dijexxc(1:lmn2_size),lmn2_size,1,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
      else
       call print_ij(dijexxc,lmn2_size,cplex_dij,lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1,opt_sym=1)
      end if
     end if
    end if

   end if

!  ------------------------------------------------------------------------
!  ----------- Final printing
!  --------------------------------------------------------------------------

   if (idij>=3.and.pawspnorb>0) then
    allocate(dijsym(cplex_dij*lmn2_size))
    dijsym(:)=paw_ij(iatom)%dij(:,3)
    if (idij==3) then
     do klmn=1,cplex_dij*lmn2_size,2
      dijsym(klmn)=dijsym(klmn)-two*dijso(klmn)
     end do
    end if
   end if

   if (abs(pawprtvol)>=1) then
    if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
     write(message, '(a)' )'   **********    TOTAL Dij in Ha   **********'
     call wrtout(6,message,'COLL')
     if (idij<=2.or.pawspnorb==0) then
      call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&      1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1)
     else
      call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&      1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1,asym_ij=dijsym)
     end if
     if (enunit>0) then
      write(message, '(a)' )'   **********    TOTAL Dij in eV   **********'
      call wrtout(6,message,'COLL')
      if (idij<=2.or.pawspnorb==0) then
       call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&       1,-1,idum,0,pawprtvol,idum,-1.d0,2)
      else
       call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,&
&       1,-1,idum,0,pawprtvol,idum,-1.d0,2,asym_ij=dijsym)
      end if
     end if
    end if
   end if

   if (abs(pawprtvol)<1) then
    if (iatom==1.or.iatom==natom.or.pawprtvol<0) then
     if (nspden==2.and.nsppol==1) then
      write(message, '(2a,i3,3a)') ch10,&
&      ' ******  TOTAL Dij in Ha (atom ',iatom,') *****',ch10,&
&      ' (antiferromagnetism case: only one spin component)'
     else
      write(message, '(a,a,i3,3a)') ch10,&
&      ' ******  TOTAL Dij in Ha (atom ',iatom,', component ',trim(dspin(idij+2*(nsploop/4))),') *****'
     end if
     call wrtout(6,message,'COLL')
     if (idij<=2.or.pawspnorb==0) then
      call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,1,-1,&
&      idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1)
     else
      call print_ij(paw_ij(iatom)%dij(:,idij),lmn2_size,cplex_dij,lmn_size,1,-1,&
&      idum,0,pawprtvol,idum,50.d0*dble(3-2*idij),1,asym_ij=dijsym)
     end if
    end if
   end if

   if (idij>=3.and.pawspnorb>0) deallocate(dijsym)

!  ----- End loops over iatom and idij
  end do

  deallocate(indklmn,ff,dijhat,dijxc)
  if (has_dijxc==1) deallocate(dijxc_hat)
  if (has_dijxc_val==1) deallocate(dijxc_val)
  if (pawspnorb>0) deallocate(dijso,dijso_rad)
  if (pawtab(itypat)%usepawu>0) deallocate(dijpawu)
  if (pawtab(itypat)%useexexch>0) deallocate(dijexxc)
  if (pawfgrtab(iatom)%gylm_allocated==2) then
   deallocate(pawfgrtab(iatom)%gylm);allocate(pawfgrtab(iatom)%gylm(0,0))
   pawfgrtab(iatom)%gylm_allocated=0
  end if

 end do

 if (pawxcdev==0) deallocate(yylmr)

 call timab(561,2,tsec)

#if defined DEBUG_MODE
 write(message,'(a)')' pawdij : exit '
 call wrtout(std_out,message,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine pawdij
!!***
