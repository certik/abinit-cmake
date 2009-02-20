!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawdenpot
!! NAME
!! pawdenpot
!!
!! FUNCTION
!! Compute different (PAW) energies densities and potentials (or potential-like quantities)
!! inside PAW spheres
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ixc= choice of exchange-correlation scheme (see above, and below)
!!  natom=number of atoms in cell.
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in unit cell.
!!  nzlmopt= if -1, compute all LM-moments of densities
!!                  initialize "lmselect" (index of non-zero LM-moments of densities)
!!           if  0, compute all LM-moments of densities
!!                  force "lmselect" to .true. (index of non-zero LM-moments of densities)
!!           if  1, compute only non-zero LM-moments of densities (stored before)
!!  option=0: compute both energies and potentials
!!         1: compute only potentials
!!         2: compute only energies
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawspnorb=flag: 1 if spin-orbit coupling is activated
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  typat(natom)=type (integer) for each atom
!!  ucvol=unit cell volume (bohr^3)
!!  xclevel= XC functional level
!!  znucl(ntypat)=gives the nuclear charge for all types of atoms
!!
!! OUTPUT
!!  paw_ij(natom)%dijhartree(lmn2_size)=Hartree contribution to dij;
!!                                      Enters into calculation of hartree energy
!!  ==== if option=0 or 2
!!    compch_sph=compensation charge inside spheres computed over spherical meshes
!!    epaw=contribution to total energy from the PAW "on-site" part
!!    epawdc=contribution to total double-counting energy from the PAW "on-site" part
!!  ==== if option=0 or 1
!!    paw_an(natom)%vxc1[m](mesh_size,:,nspden)=XC potential calculated from "on-site" density
!!    paw_an(natom)%vxct1[m](mesh_size,:,nspden)=XC potential calculated from "on-site" pseudo density
!!    ==== if paw_an(iatom)%has_vxcval==1 compute also XC potentials neglecting core charge
!!      paw_an(natom)%vxc1_val[m](mesh_size,:nspden)=XC potential calculated from spherical valence density
!!      paw_an(natom)%vxct1_val[m](mesh_size,:nspden)=XC potential calculated from spherical valence pseudo density
!!  ==== if nzlmopt==-1,
!!    paw_an(iatom)%lnmselect(lm_size,nspden)=select the non-zero LM-moments of rho1 and trho1
!!  ==== if pawspnorb>0
!!    paw_an(natom)%vh1(mesh_size,1,1)=Hartree total potential calculated from "on-site" density
!!
!! PARENTS
!!      odamix,respfn,scfcv,screening,sigma
!!
!! CHILDREN
!!      deducer0,leave_new,pawuenergy,pawxc,pawxcm,pawxenergy,pawxpot,poisson
!!      setnoccmmp,simp_gen,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawdenpot(compch_sph,epaw,epawdc,ixc,natom,nspden,ntypat,nzlmopt,option,paw_an,&
&          paw_ij,pawang,pawprtvol,pawrad,pawrhoij,pawspnorb,pawtab,pawxcdev,typat,xclevel,znucl)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13paw, except_this_one => pawdenpot
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ixc,natom,nspden,ntypat,nzlmopt,option,pawprtvol
 integer,intent(in) :: pawspnorb,pawxcdev,xclevel
 real(dp),intent(out) :: compch_sph,epaw,epawdc
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: typat(natom)
 real(dp) :: znucl(ntypat)
 type(paw_an_type),intent(inout) :: paw_an(natom)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: dum,iatom,icount,ij_size,ilm,ils,ilslm,ir
 integer :: irhoij,irhoij1,isel,ispden,itypat,itypat0,klm,klmn,klmn1,kln
 integer :: lm_size,lmax,lmin,lmn2_size,mesh_size,mm,nspdiag,nsppol,opt
 integer :: usepawu,usetcore,usexcnhat
 real(dp) :: compchspha,compchsphb,e1t10,e1xc,e1xcdc,eexc,eexcdc,eexdctemp
 real(dp) :: eexc_val,eexcdc_val
 real(dp) :: eexex,eexexdc,eextemp,eh2,eldaumdc,eldaumdcdc,etild1xc,etild1xcdc
 real(dp) :: exccore,exchmix,m1,mt1,ro,ro_dlt,ro_ql,ro_rg
 character(len=500) :: message
!arrays
 integer,allocatable :: idum1(:),idum3(:,:,:)
 logical,allocatable :: lmselect_cur(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: aa(:),bb(:),nhat(:,:,:),one_over_rad2(:)
 real(dp),allocatable :: rdum2(:,:),rdum4(:,:,:,:),rho(:),rho1(:,:,:),rho1xx(:,:,:)
 real(dp),allocatable :: trho1(:,:,:),vxc_tmp(:,:,:)

! *************************************************************************

!DEBUG
!write(6,*)' pawdenpot : enter '
!ENDDEBUG

 call timab(560,1,tsec)

 if(nzlmopt/=0.and.nzlmopt/=1.and.nzlmopt/=-1) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawdenpot : BUG -',ch10,&
&  '  invalid value for variable "nzlmopt".'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(paw_ij(1)%has_dijhartree==0) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawdenpot : BUG -',ch10,&
&  '  dijhartree must be allocated !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Various inits
 usexcnhat=maxval(pawtab(1:ntypat)%vlocopt)
 usepawu=maxval(pawtab(1:ntypat)%usepawu)
 if (option/=1) compch_sph=zero
 nspdiag=1;if (nspden==2) nspdiag=2
 nsppol=pawrhoij(1)%nsppol

!Init energies
 if (option/=1) then
  e1xc=zero     ; e1xcdc=zero
  etild1xc=zero ; etild1xcdc=zero
  exccore=zero  ; eh2=zero ; e1t10=zero
  eldaumdc=zero ; eldaumdcdc=zero
  eexex=zero    ; eexexdc=zero
  eextemp=zero  ; eexdctemp=zero
 end if

!if PAW+U, compute noccmmp^{\sigma}_{m,m'} occupation matrix
 if (usepawu>0) then
  call setnoccmmp(1,0,rdum4,0,0,idum3,natom,0,nspden,1,nsppol,0,&
&  ntypat,paw_ij,pawang,pawprtvol,pawrhoij,pawtab,rdum2,idum1,typat,0,usepawu)
 end if

!Print some titles
 if (abs(pawprtvol)>=2) then
  print *,ch10," PAW TEST:"
  if (nzlmopt<1)  print *,' ====== Moments of (n1-tn1) ========='
  if (nzlmopt==1) print *,' ==== Non-zero Moments of (n1-tn1) ===='
  print *,' The moments of (n1-tn1-nhat) must be very small...'
 end if
!================ Big loop on atoms =======================
!==========================================================

 do iatom=1,natom
  itypat=typat(iatom)
  lmn2_size=paw_ij(iatom)%lmn2_size
  lm_size=paw_an(iatom)%lm_size
  ij_size  =pawtab(itypat)%ij_size
  mesh_size=pawrad(itypat)%mesh_size
  usetcore =pawtab(itypat)%usetcore
  exchmix=pawtab(itypat)%exchmix
  if (nzlmopt<1) paw_an(iatom)%lmselect(1:lm_size)=.true.

! Allocations of "on-site" densities
  allocate(rho1 (mesh_size,lm_size,nspden))
  allocate(trho1(mesh_size,lm_size,nspden))
  allocate(nhat (mesh_size,lm_size,nspden))
  allocate(lmselect_cur(lm_size))
  rho1=zero;trho1=zero;nhat=zero
  lmselect_cur(:)=paw_an(iatom)%lmselect(:)

! Store some usefull quantities
  itypat0=0;if (iatom>1) itypat0=typat(iatom-1)
  if (itypat/=itypat0) then
   allocate(one_over_rad2(mesh_size))
   one_over_rad2(2:mesh_size)=one/pawrad(itypat)%rad(2:mesh_size)**2
  end if

! ===== Compute "on-site" densities (n1, ntild1, nhat) =====
! ==========================================================

  do ispden=1,nspden

!  -- Loop over ij channels (basis components)
   do irhoij=1,pawrhoij(iatom)%nrhoijsel
    klmn=pawrhoij(iatom)%rhoijselect(irhoij)
    klm =pawtab(itypat)%indklmn(1,klmn)
    kln =pawtab(itypat)%indklmn(2,klmn)
    lmin=pawtab(itypat)%indklmn(3,klmn)
    lmax=pawtab(itypat)%indklmn(4,klmn)

!   Retrieve rhoij
    if (nspden/=2) then
     ro=pawrhoij(iatom)%rhoijp(irhoij,ispden)
    else
     if (ispden==1) then
      ro=pawrhoij(iatom)%rhoijp(irhoij,1)&
&      +pawrhoij(iatom)%rhoijp(irhoij,2)
     else if (ispden==2) then
      ro=pawrhoij(iatom)%rhoijp(irhoij,1)
     end if
    end if
    ro=pawtab(itypat)%dltij(klmn)*ro

!   -- Computation of the moments of the densities on the spherical mesh
    do ils=lmin,lmax,2
     do mm=-ils,ils
      ilslm=ils*ils+ils+mm+1
      if (lmselect_cur(ilslm)) then
       isel=pawang%gntselect(ilslm,klm)
       if (isel>0) then

        ro_ql=ro*pawtab(itypat)%qijl(ilslm,klmn)
        ro_rg=ro*pawang%realgnt(isel)

!       == nhat(r=0)
        nhat(1,ilslm,ispden) = nhat(1,ilslm,ispden)+ro_ql*pawtab(itypat)%shapefunc(1,ils+1)
!       == rho1(r>0), trho1(r>0), nhat(r>0)
        do ir=2,mesh_size
         rho1(ir,ilslm,ispden) = rho1(ir,ilslm,ispden)&
&         +ro_rg*pawtab(itypat)%phiphj  (ir,kln)*one_over_rad2(ir)
         trho1(ir,ilslm,ispden)=trho1(ir,ilslm,ispden)&
&         +ro_rg*pawtab(itypat)%tphitphj(ir,kln)*one_over_rad2(ir)
         nhat(ir,ilslm,ispden) = nhat(ir,ilslm,ispden)+ro_ql*pawtab(itypat)%shapefunc(ir,ils+1)
        end do

       end if
      end if
     end do  ! End loops over ils, mm
    end do
   end do ! End loop over ij channels

!  Computation of rho1(r=0) and trho1(r=0)
   do ilm=1,lm_size
    if (lmselect_cur(ilm)) then
     call deducer0( rho1(:,ilm,ispden),mesh_size,pawrad(itypat))
     call deducer0(trho1(:,ilm,ispden),mesh_size,pawrad(itypat))
    end if
   end do

!  -- Test moments of densities and store non-zero ones
   if (nzlmopt==-1) then
    do ils=0,pawtab(itypat)%lcut_size-1
     do mm=-ils,ils
      ilslm=ils*ils+ils+mm+1
      m1 =maxval(abs(rho1 (1:mesh_size,ilslm,ispden)))
      mt1=maxval(abs(trho1(1:mesh_size,ilslm,ispden)))
      if (ispden==1) then
       if ((ilslm>1).and.(m1<tol16).and.(mt1<tol16)) then
        paw_an(iatom)%lmselect(ilslm)=.false.
       end if
      else if (.not.(paw_an(iatom)%lmselect(ilslm))) then
       paw_an(iatom)%lmselect(ilslm)=((m1>=tol16).or.(mt1>=tol16))
      end if
     end do
    end do
   end if

!  -- Compute integral of (n1-tn1) inside spheres
   if (option/=1.and.ispden==1) then
    allocate(aa(mesh_size))
    aa(1:mesh_size)=(rho1(1:mesh_size,1,1)-trho1(1:mesh_size,1,1)) &
&    *pawrad(itypat)%rad(1:mesh_size)**2
    call simp_gen(compchspha,aa,pawrad(itypat))
    compch_sph=compch_sph+compchspha*sqrt(four_pi)
    deallocate(aa)
   end if

!  -- Print out moments of densities (if requested)
   if (abs(pawprtvol)>=2) then
    allocate(aa(mesh_size),bb(mesh_size))
    write(message,'(2a,i3,a,i1,3a)') ch10, &
&    ' Atom ',iatom,' (ispden=',ispden,'):',ch10,&
&    '  ******* Moment of (n1-tn1) ** Moment of (n1-tn1-nhat)'
    call wrtout(6,message,'COLL')
    do ils=0,pawtab(itypat)%lcut_size-1
     do mm=-ils,ils
      ilslm=ils*ils+ils+mm+1
      if (lmselect_cur(ilslm)) then
       do ir=1,mesh_size
        ro=pawrad(itypat)%rad(ir)**(2+ils)
        aa(ir)=ro*(rho1(ir,ilslm,ispden)-trho1(ir,ilslm,ispden))
        bb(ir)=ro*nhat(ir,ilslm,ispden)
       end do
       call simp_gen(compchspha,aa,pawrad(itypat))
       call simp_gen(compchsphb,bb,pawrad(itypat))
       write(message,'(3x,a,2i2,2(a,g14.7))') &
&       'l,m=',ils,mm,': M=',compchspha,' **    M=',compchspha-compchsphb
       call wrtout(6,message,'COLL')
      end if
     end do
    end do
    deallocate(aa,bb)
   end if

!  ----- End loop over spin components
  end do

! =========== Compute XC potentials and energies ===========
! ==========================================================

! Temporary storage
  if (pawxcdev/=0) allocate(vxc_tmp(mesh_size,lm_size,nspden))
  if (pawxcdev==0) allocate(vxc_tmp(mesh_size,pawang%angl_size,nspden))
  dum=0
! ===== Vxc1 term =====
  if (pawxcdev/=0) then
   call pawxcm(pawtab(itypat)%coredens,eexc,eexcdc,dum,ixc,lm_size,&
&   paw_an(iatom)%lmselect,nhat,nspden,option,&
&   pawang,pawrad(itypat),pawxcdev,rho1,1,0,vxc_tmp,xclevel)
  else
   call pawxc(pawtab(itypat)%coredens,eexc,eexcdc,ixc,lm_size,&
&   paw_an(iatom)%lmselect,nhat,nspden,option,&
&   pawang,pawrad(itypat),rho1,1,0,vxc_tmp,xclevel)
  end if
  if (option/=1) then
   e1xc=e1xc+eexc
   e1xcdc=e1xcdc+eexcdc
  end if
  if (option<2) paw_an(iatom)%vxc1(:,:,:)=vxc_tmp(:,:,:)

! ===== tVxc1 term =====
  if (pawxcdev/=0) then
   call pawxcm(pawtab(itypat)%tcoredens,eexc,eexcdc,dum,ixc,lm_size,&
&   paw_an(iatom)%lmselect,nhat,nspden,option,&
&   pawang,pawrad(itypat),pawxcdev,trho1,usetcore,1+usexcnhat,vxc_tmp,xclevel)
  else
   call pawxc(pawtab(itypat)%tcoredens,eexc,eexcdc,ixc,lm_size,&
&   paw_an(iatom)%lmselect,nhat,nspden,option,&
&   pawang,pawrad(itypat),trho1,usetcore,1+usexcnhat,vxc_tmp,xclevel)
  end if
  if (option/=1) then
   etild1xc=etild1xc+eexc
   etild1xcdc=etild1xcdc+eexcdc
  end if
  if (option<2) paw_an(iatom)%vxct1(:,:,:)=vxc_tmp(:,:,:)

! =========== Compute valence-only XC potentials ===========
! ==========================================================
  if (paw_an(iatom)%has_vxcval==1.and.(option==0.or.option==1)) then
   if (.not.associated(paw_an(iatom)%vxc1_val).or..not.associated(paw_an(iatom)%vxct1_val)) then
    write(message,'(3a)')' pawdenpot : BUG ',ch10,' vxc1_val and vxct1_val must be associated'
    call wrtout(std_out,message,'COLL') 
    call leave_new('COLL')
   end if
!  ===== Vxc1_val term, vxc[n1] =====
   if (pawxcdev/=0) then
    print*,' pawdenpot : Computing valence-only v_xc[n1] using moments '
    write(*,*)'Min density rho1 = ',MINVAL(rho1)
    call pawxcm(pawtab(itypat)%coredens,eexc_val,eexcdc_val,dum,ixc,lm_size,&
&    paw_an(iatom)%lmselect,nhat,nspden,option,&
&    pawang,pawrad(itypat),pawxcdev,rho1,0,0,vxc_tmp,xclevel)
!   &    pawang,pawrad(itypat),pawxcdev,rho1,0,1+usexcnhat,vxc_tmp,xclevel)
   else
    print*,' pawdenpot : Computing valence-only v_xc[n1] using angular mesh '
    call pawxc(pawtab(itypat)%coredens,eexc_val,eexcdc_val,ixc,lm_size,&
&    paw_an(iatom)%lmselect,nhat,nspden,option,&
&    pawang,pawrad(itypat),rho1,0,0,vxc_tmp,xclevel)
!   &    pawang,pawrad(itypat),rho1,0,1+usexcnhat,vxc_tmp,xclevel)
   end if
!  if (option/=1) then
!  e1xc_val=e1xc_val+eexc_val
!  e1xcdc_val=e1xcdc_val+eexcdc_val
!  end if
   if (option<2) paw_an(iatom)%vxc1_val(:,:,:)=vxc_tmp(:,:,:)

!  ===== tVxc1_val term =====
   if (pawxcdev/=0) then
    print*,' pawdenpot : Computing valence-only v_xc[tn1+nhat] using moments '
    write(*,*)'Min density trho1 = ',MINVAL(trho1)
    write(*,*)'Min density trho1 + nhat = ',MINVAL(trho1+nhat)
    call pawxcm(pawtab(itypat)%tcoredens,eexc_val,eexcdc_val,dum,ixc,lm_size,&
&    paw_an(iatom)%lmselect,nhat,nspden,option,&
&    pawang,pawrad(itypat),pawxcdev,trho1,0,1+usexcnhat,vxc_tmp,xclevel)
   else
    print*,' pawdenpot : Computing valence-only v_xc[tn1+nhat] using angular mesh'
    call pawxc(pawtab(itypat)%tcoredens,eexc_val,eexcdc_val,ixc,lm_size,&
&    paw_an(iatom)%lmselect,nhat,nspden,option,&
&    pawang,pawrad(itypat),trho1,0,1+usexcnhat,vxc_tmp,xclevel)
   end if
!  if (option/=1) then
!  etild1xc_val=etild1xc_val+eexc_val
!  etild1xcdc_val=etild1xcdc_val+eexcdc_val
!  end if
   if (option<2) paw_an(iatom)%vxct1_val(:,:,:)=vxc_tmp(:,:,:)
  end if ! valence-only XC potentials

  deallocate(vxc_tmp)

! ===== Compute first part of local exact-exchange energy term =====
! ===== Also compute corresponding potential                   =====
! ==================================================================

  if (pawtab(itypat)%useexexch>0) then

!  ===== Re-compute a partial "on-site" density n1 (only l=lexexch contrib.)
   allocate(rho1xx(mesh_size,lm_size,nspden));rho1xx=zero

   do ispden=1,nspden

!   -- Loop over ij channels (basis components)
    do irhoij=1,pawrhoij(iatom)%nrhoijsel
     klmn=pawrhoij(iatom)%rhoijselect(irhoij)
     if(pawtab(itypat)%indklmn(3,klmn)==0.and.&
&     pawtab(itypat)%indklmn(4,klmn)==2*pawtab(itypat)%lexexch) then
      klm =pawtab(itypat)%indklmn(1,klmn)
      kln =pawtab(itypat)%indklmn(2,klmn)
      lmin=pawtab(itypat)%indklmn(3,klmn)
      lmax=pawtab(itypat)%indklmn(4,klmn)
      if (nspden/=2) then
       ro=pawrhoij(iatom)%rhoijp(irhoij,ispden)
      else
       if (ispden==1) then
        ro=pawrhoij(iatom)%rhoijp(irhoij,1)&
&        +pawrhoij(iatom)%rhoijp(irhoij,2)
       else if (ispden==2) then
        ro=pawrhoij(iatom)%rhoijp(irhoij,1)
       end if
      end if
      ro=pawtab(itypat)%dltij(klmn)*ro

      do ils=lmin,lmax,2
       do mm=-ils,ils
        ilslm=ils*ils+ils+mm+1
        if (lmselect_cur(ilslm)) then
         isel=pawang%gntselect(ilslm,klm)
         if (isel>0) then
          ro_rg=pawang%realgnt(isel)*ro
!         -- rho1xx(r>0)
          do ir=2,mesh_size
           rho1xx(ir,ilslm,ispden)=rho1xx(ir,ilslm,ispden)&
&           +ro_rg*pawtab(itypat)%phiphj(ir,kln)*one_over_rad2(ir)
          end do
         end if
        end if
       end do  ! End loops over ils, mm
      end do
     end if
    end do ! End loop over ij channels

!   -- rho1xx(r=0)
    do ilm=1,lm_size
     if (lmselect_cur(ilm)) call deducer0(rho1xx(:,ilm,ispden),mesh_size,pawrad(itypat))
    end do

!   ----- End loop over spin components
   end do

!  ===== Re-compute Exc1 and Vxc1; for local exact-exchange, this is done in GGA only
   allocate(vxc_tmp(mesh_size,lm_size,nspden))
   call pawxcm(pawtab(itypat)%coredens,eextemp,eexdctemp,pawtab(itypat)%useexexch,ixc,lm_size,&
&   paw_an(iatom)%lmselect,nhat,nspden,option,pawang,pawrad(itypat),pawxcdev,&
&   rho1xx,0,0,vxc_tmp,xclevel)
   if (option/=1) then
    e1xc=e1xc-eextemp*exchmix
    e1xcdc=e1xcdc-eexdctemp*exchmix
   end if
   if (option<2) paw_an(iatom)%vxc_ex(:,:,:)=vxc_tmp(:,:,:)
   deallocate(rho1xx,vxc_tmp)

  end if ! useexexch

  itypat0=0;if (iatom<natom) itypat0=typat(iatom+1)
  if (itypat/=itypat0) deallocate(one_over_rad2)
  deallocate(lmselect_cur)

! ==== Compute Hartree potential terms and some energy terms ====
! ===============================================================

  paw_ij(iatom)%dijhartree=zero
  do ispden=1,nspdiag
   do irhoij=1,pawrhoij(iatom)%nrhoijsel
    klmn=pawrhoij(iatom)%rhoijselect(irhoij)
    ro_dlt=pawrhoij(iatom)%rhoijp(irhoij,ispden)*pawtab(itypat)%dltij(klmn)
    paw_ij(iatom)%dijhartree(klmn)=paw_ij(iatom)%dijhartree(klmn)&
&    +ro_dlt*pawtab(itypat)%eijkl(klmn,klmn)
    do klmn1=1,klmn-1
     paw_ij(iatom)%dijhartree(klmn1)=paw_ij(iatom)%dijhartree(klmn1)&
&     +ro_dlt*pawtab(itypat)%eijkl(klmn1,klmn)
    end do
    do klmn1=klmn+1,lmn2_size
     paw_ij(iatom)%dijhartree(klmn1)=paw_ij(iatom)%dijhartree(klmn1)&
&     +ro_dlt*pawtab(itypat)%eijkl(klmn,klmn1)
    end do
   end do
  end do
  if (option/=1) then
   do ispden=1,nspdiag
    do irhoij=1,pawrhoij(iatom)%nrhoijsel
     klmn=pawrhoij(iatom)%rhoijselect(irhoij)
     ro_dlt=pawrhoij(iatom)%rhoijp(irhoij,ispden)*pawtab(itypat)%dltij(klmn)
     eh2=eh2    +ro_dlt*paw_ij(iatom)%dijhartree(klmn)
     e1t10=e1t10+ro_dlt*pawtab(itypat)%dij0(klmn)
    end do
   end do
  end if
  if (option/=1) exccore=exccore+pawtab(itypat)%exccore

! Compute 1st moment of total Hartree potential VH(n_Z+n_core+n1)
! Used for spin-orbit contribution
  if (pawspnorb>0.and.option<2) then
   allocate(rho(mesh_size))
   rho(1:mesh_size)=(rho1(1:mesh_size,1,1)+sqrt(four_pi)*pawtab(itypat)%coredens(1:mesh_size)) &
&   *four_pi*pawrad(itypat)%rad(1:mesh_size)**2
   call poisson(rho,0,ro,pawrad(itypat),paw_an(iatom)%vh1(:,1,1))
   paw_an(iatom)%vh1(2:mesh_size,1,1)=(paw_an(iatom)%vh1(2:mesh_size,1,1) &
   -sqrt(four_pi)*znucl(itypat))/pawrad(itypat)%rad(2:mesh_size)
   call deducer0(paw_an(iatom)%vh1,mesh_size,pawrad(itypat))
   deallocate(rho)
  end if
  deallocate(nhat,rho1,trho1)

! ========= Compute PAW+U and energy contribution  =========
! ==========================================================

  if (pawtab(itypat)%usepawu>0.and.option/=1) then
   call pawuenergy(iatom,eldaumdc,eldaumdcdc,pawprtvol,pawtab(itypat),paw_ij(iatom),nspden)
  end if

! === Compute 2nd part of local exact-exchange energy and potential  ===
! ======================================================================

  if (pawtab(itypat)%useexexch>0) then

   if(nspden==4)  then
    write(message, '(4a)' ) ch10,&
&    '  pawdenpot : ERROR -',ch10,&
&    '  Local exact-exch. not implemented for nspden=4 !'
    call wrtout(ab_out,message,'COLL')
    call wrtout(06,  message,'COLL')
    call leave_new('COLL')
   end if

   if (option<2) call pawxpot(nspden,pawprtvol,pawtab(itypat),paw_ij(iatom),pawrhoij(iatom))
   if (option/=1) then
    write(message, '(2a)' )ch10,'======= PAW local exact exchange terms (in Hartree) ===='
    call wrtout(06,  message,'COLL')
    write(message, '(2a,i4)' )ch10,' For Atom',iatom
    call wrtout(06,  message,'COLL')
    call pawxenergy(eexex,eexexdc,pawprtvol,pawrhoij(iatom),pawtab(itypat),nspden)
   end if

  end if ! useexexch

! =========== End loop on atoms ============================
! ==========================================================

 end do

!========== Assemble "on-site" energy terms ===============
!==========================================================

 if (option/=1) then
  epaw  =(e1xc+half*eh2+e1t10-exccore) -(etild1xc)            +eldaumdc  +eexex*exchmix
  epawdc=(e1xc-e1xcdc-half*eh2-exccore)-(etild1xc-etild1xcdc) +eldaumdcdc-eexex*exchmix
 end if

 call timab(560,2,tsec)

!DEBUG
!write(6,*)' pawdenpot : exit '
!ENDDEBUG

end subroutine pawdenpot
!!***
