!{\src2tex{textfont=tt}}
!!****f* ABINIT/paw_symcprj
!! NAME
!! paw_symcprj
!!
!! FUNCTION
!! Symmetrize the projected wavefunctions <Proj_i|Cnk>. This routine calculates the matrix elements
!! at each k-point in the full Brillouin zone starting from the knowledge of <Proj_i|Cnk>
!! in the irreducible wedge.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Cprj_ibz(natom,nspinor*nkibz*mband*nsppol) <type(cprj_type)>= projected wave functions <Proj_i|Cnk> 
!!    with all NL projectors. According to ctoprj.F90, matrix elements are unsorted.
!!  dimlmn(natom)=number of (l,m,n) components for each atom
!!  Psps<pseudopotential_type>
!!     %indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!  Kmesh<bz_mesh_type>: datatype gathering information on the k-point sampling 
!!     %nbz=number of k-points in the full Brillouin zone
!!     %nibz=number of k-points in the irreducible wedge
!!     %tab(nkbz)=table giving for each k-point in the BZ (array kbz), the corresponding irred. point in the IBZ.
!!       i.e k_BZ = (IS) kIBZ where S is one of the symrec operations and I is the inversion or the identity
!!     %tabi(nkbz)=for each k-point in the BZ defines whether inversion has to be considered in the 
!!       relation k_BZ=(IS) k_IBZ (1 => only S; -1 => -S)  
!!     %tabo(nkbz)= the symmetry operation S that takes k_IBZ to each k_BZ
!!  Pawang <type(pawang_type)>=paw angular mesh and related data
!!     %lmax=Max angular momentum included in the PAW datasets used. mentioned at the second line of the psp file
!!     %zarot(2*lmax+1,2*lmax+1,lmax+1,nsym)=coefficients of the transformation of real spherical 
!!      harmonics under the symmetry operations. 
!!     FIXME this part is commented in setsymrhoij 
!!    ??? In case of antiferromagnetism, zarot coeffs contain (-1)^l factors from spin inversion (See setsymrhoij). ???
!!  mband=Max number of bands
!!  Cryst<Crystal_structure>=data type gathering information on unit cell and symmetries
!!     %ntypat=number of type of atoms
!!     %natom=number of atoms in the unit cell
!!     %typat(natom)=type of each atom
!!     %indsym(4,nsym,natom)=indirect indexing array: 
!!      for each isym,iatom, fourth element is label of atom into which iatom is sent by the INVERSE of the 
!!      symmetry operation symrel(isym); first three elements are the primitive translations that must be subtracted 
!!      after the transformation to get back to the original unit cell.
!!  nspinor=number of spinorial components 
!!  nsppol=number of independent spin polarizations
!!  nband(Kmesh%nibz*nsppol)=number of bands for each k-point
!!  Pawtab(Cryst%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!     %lmn_size=number of lmn channels 
!!
!! OUTPUT
!! Cprj_bz(natom,nspinor*mband*kmesh%nbz*nsppol) <type(cprj_type)>= 
!!  projected wave functions <Proj_i|Cnk> for each k-point in the full Brillouin zone
!!  
!! SIDE EFFECTS
!!
!! NOTES
!!  * TODO Anti-ferromagnetism not tested
!!  * Thanks to space group symmetries, the following relationship holds:
!!
!!     $ <p_{nlm}^a|\tilde\psi_{Sk}> = e^{ik.R_o} \sum_\alpha D^l_{\alpha,m}(S) <p_{nl\alpha}^{a\prime}|\tilde\psi_k>
!!      where $R^{-1}(r_a-\tau) = r_{a\prime} + R_o$. 
!!
!!      1) $R^{-1}$ is a rotation in real space and $\tau$ is the associated fractional translation.
!!      2) $r_a$ and $r_a\prime$ are the positions of two atoms of the same type inside the first unit cell 
!!      3) $R_o$ is the lattice vector required to bring the rotated atom back to the first unit cell.
!!
!!  *) Time-reversal symmetry can be included by simply noticing that:
!!      1) For collinear magnetism $ <p_{nlm}^a |\tilde\psi_{-k}> = <p_{nlm}^a|\tilde\psi_k>^\dagger.
!!      2) For noncollinear magnetism $ \phi_{-k}^1 =  (phi_k^2)^* $
!!                                    $ \phi_{-k}^2 = -(phi_k^1)^* $
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

subroutine paw_symcprj(Pawtab,Cryst,nspinor,mband,nband,&
& nsppol,Psps,Kmesh,Cprj_ibz,Pawang,dimlmn,Cprj_bz)

 use defs_basis
 use defs_datatypes
 use m_errors, only : assert
 use m_io_tools, only : flush_unit, get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_15gw, except_this_one => paw_symcprj
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nspinor,nsppol
 type(BZ_mesh_type),intent(in) :: Kmesh
 type(Crystal_structure),intent(in) :: Cryst
 type(Pawang_type),intent(in) :: Pawang
 type(Pseudopotential_type),intent(in) :: Psps
!arrays
 integer,intent(in) :: dimlmn(Cryst%natom),nband(Kmesh%nibz*nsppol)
 type(Cprj_type),intent(in) :: Cprj_ibz(Cryst%natom,nspinor*mband*Kmesh%nibz*nsppol)
 type(Cprj_type),intent(out) :: Cprj_bz(Cryst%natom,nspinor*mband*Kmesh%nbz*nsppol)
 type(Pawtab_type),intent(in) :: Pawtab(Cryst%ntypat)

!Local variables-------------------------------
!scalars
 integer :: iat,iat_sym,iband,ibg_full,ibg_irr,ibsp,ibsp_bz
 integer :: ibsp_ibz,ii,ik_bz,ik_ibz,ik_sym,indexj,ispinor,isppol,isym,itim
 integer :: itypat,j0lmn,jl,jl0,jlmn,jln,jln0,jlpm,jm,jn,lmax,mm,natom,nband_k
 integer :: nlmn,unt
 real(dp) :: arg,wtk
 complex(gwpc) :: ph_mkbzt
 logical :: DEBUG,ltest
 character(len=500) :: msg
 character(len=fnlen) :: fnam
!arrays
 integer :: G0(3),R0(3)
 integer,allocatable :: ibg_bz(:,:),ibg_ibz(:,:)
 real(dp) :: dum(2,nspinor),kbz(3),kirr(3),phase(2),swp(2),tmp(2,nspinor)
 real(dp),allocatable :: DS_mmpl(:,:,:)
 type(Cprj_type),allocatable :: Cprjnk_kibz(:,:)

! *************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' paw_symcprj : enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 !
 ! === Perform some check ===
 if (nspinor==2) write(*,*)'WARNING - paw_symcprj with spinor! spin rotation not implemented!' 

 natom=Cryst%natom
 !
 ! === Fill tables giving starting index for cprjnk_kibz and cprjnk_kbz ===
 allocate(ibg_ibz(Kmesh%nibz,nsppol)) 
 allocate(ibg_bz (Kmesh%nbz, nsppol)) 
 ii=0
 do isppol=1,nsppol
  do ik_ibz=1,Kmesh%nibz
   ibg_ibz(ik_ibz,isppol)=ii
   nband_k=nband(ik_ibz+(isppol-1)*Kmesh%nibz)
   ii=ii+nspinor*nband_k 
  end do
 end do
 ii=0
 do isppol=1,nsppol
  do ik_bz=1,Kmesh%nbz
   ibg_bz(ik_bz,isppol)=ii
   ik_ibz=Kmesh%tab(ik_bz)
   nband_k=nband(ik_ibz+(isppol-1)*Kmesh%nibz)
   ii=ii+nspinor*nband_k 
  end do
 end do

 do ibsp=1,SIZE(Cprj_bz,DIM=2)
  do iat=1,SIZE(Cprj_bz,DIM=1)
   Cprj_bz(iat,ibsp)%cp(:,:)=zero
  end do
 end do

 lmax=Pawang%l_max-1 ! l_max is Max l+1
 allocate(DS_mmpl(2*lmax+1,2*lmax+1,lmax+1))
 !
 ! === Loop on spin and k-points in the BZ to be symmetrized ===
 do isppol=1,nsppol
  do ik_bz=1,Kmesh%nbz 

   !=== Get index of IBZ image of this BZ k-point and related simmetry ===
   call get_BZ_item(Kmesh,ik_bz,kbz,ik_ibz,isym,itim,ph_mkbzt)
   !
   ! === DS_mmpl is the rotation matrix for real spherical harmonics associated to symrec(:,:,isym) ===
   ! * Note the convention used by Blanco in Eq. 27 : DS_mmp multiply spherical harmonics as row vectors
   DS_mmpl(:,:,:)=Pawang%zarot(:,:,:,isym)
   call get_IBZ_item(Kmesh,ik_ibz,kirr,wtk)

   ibg_irr =ibg_ibz(ik_ibz,isppol)
   ibg_full=ibg_bz (ik_bz ,isppol)
   nband_k=nband(ik_ibz+(isppol-1)*Kmesh%nibz)

   ! === Load projected wave functions for this irred k-point and spin ===
   allocate(Cprjnk_kibz(natom,nspinor*nband_k)) 
   call cprj_alloc(Cprjnk_kibz,0,dimlmn)
   ibsp=0
   do iband=1,nband_k
    do ispinor=1,nspinor
     ibsp=ibsp+1
     do iat=1,natom
      Cprjnk_kibz(iat,ibsp)%cp(:,:)=Cprj_ibz(iat,ibg_irr+ibsp)%cp(:,:)
     end do
    end do
   end do

   ! === Loop over atoms to be symmetrized for this k-point ===
   ! * R^{-1} (xred(:,iat)-tnons) = xred(:,iat_sym) + R0.
   do iat=1,natom
    itypat=Cryst%typat(iat)
    iat_sym=Cryst%indsym(4,isym,iat) ; R0(:)=Cryst%indsym(1:3,isym,iat)
    arg=two_pi*DOT_PRODUCT(kirr,R0)
    phase(1)=COS(arg)
    phase(2)=SIN(arg)
    !
    ! === Loop on (jl,jm,jn) components to be symmetrized ===
    nlmn=Pawtab(itypat)%lmn_size
    jl0=-1 ; jln0=-1 ; indexj=1
    do jlmn=1,nlmn
     jl  =Psps%indlmn(1,jlmn,itypat)
     jm  =Psps%indlmn(2,jlmn,itypat)
     jn  =Psps%indlmn(3,jlmn,itypat)
     jln =Psps%indlmn(5,jlmn,itypat) 
     jlpm=1+jl+jm 
     if (jln/=jln0) indexj=indexj+2*jl0+1

     ! === For each band, calculate contribution due to rotated real spherical harmonics ===
     ! FIXME check this expression; according to Blanco I should have D(S^-1} but it seems D(S) is correct
     ! Recheck spinorial case, presently is wrong
     ibsp_ibz=0
     ibsp_bz=0
     do iband=1,nband_k

      tmp(:,:)=zero
      do ispinor=1,nspinor
       ibsp_ibz=ibsp_ibz+1
       do mm=1,2*jl+1
        tmp(1,ispinor)=tmp(1,ispinor)+DS_mmpl(mm,jlpm,jl+1)*Cprjnk_kibz(iat_sym,ibsp_ibz)%cp(1,indexj+mm)
        tmp(2,ispinor)=tmp(2,ispinor)+DS_mmpl(mm,jlpm,jl+1)*Cprjnk_kibz(iat_sym,ibsp_ibz)%cp(2,indexj+mm)
       end do
      end do !ispinor
      !
      ! === Apply phase to account if the symmetric atom belongs to a different unit cell ===
      do ispinor=1,nspinor
       dum(1,ispinor)=tmp(1,ispinor)*phase(1)-tmp(2,ispinor)*phase(2)
       dum(2,ispinor)=tmp(1,ispinor)*phase(2)+tmp(2,ispinor)*phase(1)
      end do
      !
      ! * In case apply Time-reversal symmetry to retrieve the right k.
       if (itim==2) then
        if (nspinor==1) then 
         dum(2,1)=-dum(2,1)
        else if (nspinor==2) then
         ! TODO rotate wavefunction in spinor space.
         swp(:)=dum(:,1)
         dum(1,1)= dum(1,2)
         dum(2,1)=-dum(2,2)
         dum(1,2)=-swp(1)
         dum(2,2)= swp(2)
        end if
       end if
       !
       ! === Save values ===
       do ispinor=1,nspinor
        ibsp_bz=ibsp_bz+1
        Cprj_bz(iat,ibg_full+ibsp_bz)%cp(1,jlmn)=dum(1,ispinor)
        Cprj_bz(iat,ibg_full+ibsp_bz)%cp(2,jlmn)=dum(2,ispinor)
       end do

     end do !iband

     jl0=jl ; jln0=jln
    end do !jlmn
   end do !iat

   call cprj_free(Cprjnk_kibz) ; deallocate(Cprjnk_kibz)
  end do !ikbz
 end do !isppol

 DEBUG=.FALSE.
 if (DEBUG) then 
  unt=get_unit()
  fnam="pristine_PJ_IRR"
  call isfile(fnam,'new')
  open(unit=unt,file=fnam,form="formatted")
  do isppol=1,nsppol
   do ik_ibz=1,Kmesh%nibz
    !if (ANY(ABS(Kmesh%ibz(:,ik_ibz))>tol6)) cycle
    ibg_irr =ibg_ibz(ik_ibz,isppol)
    nband_k=nband(ik_ibz+(isppol-1)*Kmesh%nibz)
    write(unt,'(2(a,i3))')" -- spin ",isppol," nband ",nband_k
    write(unt,'(a,3f7.3)')" -- kpoint ",Kmesh%ibz(:,ik_ibz)
     do iat=1,natom 
      do iband=1,nband_k
       nlmn=dimlmn(iat)
       write(unt,'(2i2,/,(2f11.7))')iat,iband,(Cprj_ibz(iat,iband+ibg_irr)%cp(1:2,1:nlmn))
       write(unt,*)"----"
     end do
    end do
   end do
  end do
  close(unt)

  fnam="PJ_IRR"
  call isfile(fnam,'new')
  open(unit=unt,file=fnam,form="formatted")
  do isppol=1,nsppol
   do ik_bz=1,Kmesh%nbz
    !if (ANY(ABS(Kmesh%bz(:,ik_bz))>tol6)) cycle
    ik_ibz=Kmesh%tab (ik_bz)
    ibg_full=ibg_bz(ik_bz,isppol)
    ibg_irr =ibg_ibz(ik_ibz,isppol)
    if (is_samek(Kmesh%ibz(:,ik_ibz),Kmesh%bz(:,ik_bz),G0)) then 
     nband_k=nband(ik_ibz+(isppol-1)*Kmesh%nibz)
     write(unt,'(2(a,i3))')" -- spin ",isppol," nband ",nband_k
     write(unt,'(2(a,3f7.3))')" -- kpoint ",Kmesh%bz(:,ik_bz)," --> irred ",Kmesh%ibz(:,ik_ibz)
     do iat=1,natom 
      do iband=1,nband_k
       nlmn=dimlmn(iat)
       !write(unt,'(2i2,/,(2f11.7))')iat,iband,(Cprj_bz(iat,iband+ibg_full)%cp(1:2,1:nlmn))
       write(unt,'(2i2,/,(2f11.7))')iat,iband,&
&       Cprj_bz(iat,iband+ibg_full)%cp(1:2,1:nlmn)-Cprj_ibz(iat,iband+ibg_irr)%cp(1:2,1:nlmn)
       write(unt,*)"----"
      end do
     end do
    end if
   end do
  end do
  close(unt)

  fnam="PJ_FULL"
  call isfile(fnam,'new')
  open(unit=unt,file=fnam,form="formatted")
  do isppol=1,nsppol
   do ik_bz=1,Kmesh%nbz
    ik_ibz=Kmesh%tab (ik_bz)
    ibg_full =ibg_bz(ik_bz,isppol)
    nband_k=nband(ik_ibz+(isppol-1)*Kmesh%nibz)
    write(unt,'(2(a,i3))')" -- spin ",isppol," nband ",nband_k
    write(unt,'(2(a,3f7.3))')" -- kpoint ",Kmesh%bz(:,ik_bz)," --> irred ",Kmesh%ibz(:,ik_ibz)
    do iat=1,natom 
     do iband=1,nband_k
      nlmn=dimlmn(iat)
      write(unt,'(2i2,/,(2f11.7))')iat,iband,(Cprj_bz(iat,iband+ibg_full)%cp(1:2,1:nlmn))
      write(unt,*)"----"
     end do
    end do
   end do
  end do
  close(unt)
 end if

 deallocate(ibg_ibz,ibg_bz)
 deallocate(DS_mmpl)

#if defined DEBUG_MODE
 write(msg,'(a)')' paw_symcprj : exit '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine paw_symcprj
!!***
