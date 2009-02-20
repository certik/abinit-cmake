!{\src2tex{textfont=tt}}
!!****f* ABINIT/gipaw_aug_fields
!! NAME
!! gipaw_aug_fields
!!
!! FUNCTION
!! Compute the diamagnetic and paramagnetic augmentation fields in GIPAW
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (JJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! integer :: natom,nfft,ntypat: number of atoms, number of grid points, number of atom types
!! type(pseudopotential_type) psps: datastructure of pseudopotential information
!! integer typat(natom): array of atom types
!! type(pawfgrtab_type) pawfgrtab(natom): data on fine grid points in paw spheres around each atom
!! type(pawrad_type) pawrad(ntypat): paw radial mesh data
!! type(pawrhoij_type) pawrhoij(natom): data on paw density around each atom
!! type(pawtab_type) pawtab(ntypat): paw wavefunctions around each type of atom
!!
!! OUTPUT
!! type(gipaw_type) gipaw_aug(natom): diamagnetic and paramagnetic fields in each PAW sphere and
!!   on-site angular momentum around each PAW sphere
!!
!! NOTES
!! The diamagnetic field is a real scalar field, given by 
!! $\langle phi_j|r'><r'|phi_i\rangle - \langle\tilde{phi}_j|r'><r'|\tilde{\phi}_i\rangle$. 
!! This form requires the computation of u(r)*Y_lm(r)/r at the point r', on the fine grid
!! around the atom. 
!!
!! The paramagnetic field is a purely imaginary vector field. To obtain it,
!! multiply the values in gipaw_aug%para computed here by (-i). It is similar to the diamagnetic
!! field except that the operator is not $|r'\rangle\langle r'|$ but 
!! $-\frac{1}{2}(\mathbf{p}|r'\rangle\langle r'|+|r'\rangle\langle r'|\mathbf{p})$, where
!! $\mathbf{p}=-i\nabla$. This requires computing also the gradient of u(r)*Y_lm(r)/r at each
!! on the grid around the atom. See Yates, Pickard, Mauri, PRB 76, 024401 (2007) Eqs. 18 and 19.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine gipaw_aug_fields(gipaw_aug,natom,nfft,ntypat,pawfgrtab,pawrad,pawrhoij,pawtab,psps,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_13paw, except_this_one => gipaw_aug_fields
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nfft,ntypat
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: typat(natom)
 type(gipaw_type),intent(out) :: gipaw_aug(natom)
 type(pawfgrtab_type),intent(in) :: pawfgrtab(natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: iatom,idir,ifgd,ifft,ifftsph,inl,inrm,ipsang,itypat
 integer :: il,ilmn,im,ilm,iln,jl,jlmn,jlm,jln,jm,j0lmn
 integer :: klmn,kln,mesh_size,nder,nfgd,nnl,normchoice,option
 real(dp) :: cr,di,dj,dti,dtj,phi,phj,qji,rR,tphi,tphj,tui,tuj,ui,uj
 complex(dpc) :: sls_val
!arrays
 integer,allocatable :: nrm_ifftsph(:)
 real(dp),allocatable :: dphigrd(:,:),dtphigrd(:,:),ff(:),nrm(:),phigrd(:,:)
 real(dp),allocatable :: tphigrd(:,:),ylm(:,:),ylmgr(:,:,:)

! ************************************************************************

!DEBUG
!write(*,*)' gipaw_aug_fields : enter'
!ENDDEBUG

!loop over atoms in cell
 do iatom = 1, natom
  itypat = typat(iatom)
  nfgd = pawfgrtab(iatom)%nfgd ! number of points in the fine grid for this PAW sphere
  nnl = pawtab(itypat)%basis_size ! number of nl elements in PAW basis

! need mesh for calculation of <L> field
! ----
  mesh_size=pawrad(itypat)%mesh_size ! PAW radial grid mesh
  allocate(ff(mesh_size))

! obtain |r-R| values on fine grid
  allocate(nrm(nfgd))
  do ifgd=1, nfgd
   nrm(ifgd) = sqrt(DOT_PRODUCT(pawfgrtab(iatom)%rfgd(:,ifgd),pawfgrtab(iatom)%rfgd(:,ifgd)))
  end do ! these are the |r-R| values

! compute Ylm and grad(Ylm) for each r-R vector. 
! ----
  ipsang = 1 + (pawtab(itypat)%l_size - 1)/2 ! recall l_size=2*l_max+1
  allocate(ylm(ipsang*ipsang,nfgd))
  allocate(ylmgr(3,ipsang*ipsang,nfgd))
  option = 2 ! compute Ylm(r-R) and gradients for vectors
  normchoice = 1 ! use computed norms of input vectors
  call initylmr(ipsang,normchoice,nfgd,nrm,option,pawfgrtab(iatom)%rfgd,ylm,ylmgr)

! in order to do spline fits, the |r-R| data must be sorted
! ----
  allocate(nrm_ifftsph(nfgd))
  nrm_ifftsph(:) = pawfgrtab(iatom)%ifftsph(:) ! copy of indices of points, to be rearranged by sort_dp
  call sort_dp(nfgd,nrm,nrm_ifftsph,tol8) ! sort the nrm points, keeping track of which goes where

! now make spline fits of phi, dphi/dr, tphi, and d tphi/dr onto the fine grid around the atom
! ----
  allocate(phigrd(nfgd,nnl),tphigrd(nfgd,nnl))
  allocate(dphigrd(nfgd,nnl),dtphigrd(nfgd,nnl))
  call spline_paw_fncs(dphigrd,dtphigrd,nnl,nfgd,pawrad(itypat),pawtab(itypat),nrm,phigrd,tphigrd)

! loop over basis elements for this atom
! because we have to store things like <phi|r'><r'|phi>-<tphi|r'><r'|tphi> at each point of the 
! fine grid, there is no integration, and hence no simplifications of the Y_lm's. That's why 
! we have to loop through the basis elements in exhaustive detail, rather than just a loop over
! lmn2_size or something comparable.
! ----
  do jlmn=1,pawtab(itypat)%lmn_size
   jl= psps%indlmn(1,jlmn,itypat)  
   jm=psps%indlmn(2,jlmn,itypat)
   jlm = psps%indlmn(4,jlmn,itypat)
   jln=psps%indlmn(5,jlmn,itypat)
   j0lmn=jlmn*(jlmn-1)/2
   do ilmn=1,jlmn
    il= psps%indlmn(1,ilmn,itypat)
    im=psps%indlmn(2,ilmn,itypat)
    iln=psps%indlmn(5,ilmn,itypat)
    ilm = psps%indlmn(4,ilmn,itypat)
    klmn=j0lmn+ilmn
    kln = pawtab(itypat)%indklmn(2,klmn) ! need this for <L> field below
    do ifgd=1, nfgd ! loop over fine grid points in current PAW sphere
     ifftsph = pawfgrtab(iatom)%ifftsph(ifgd) ! index of the point on the grid

!    have to retrieve the spline point to use since these were sorted
!    ----
     do inrm=1, nfgd 
      if(nrm_ifftsph(inrm) == ifftsph) exit ! have found nrm point corresponding to nfgd point
     end do ! now inrm is the index of the sorted nrm vector to use

!    avoid division by zero
!    ---
     if(nrm(inrm) > zero) then 
      rR = nrm(inrm) ! value of |r-R| in the following

!     recall that <r|phi>=u(r)*Slm(r^)/r, here are the u and u*S parts
!     ----
      uj   = phigrd(inrm,jln)
      ui   = phigrd(inrm,iln)
      tuj  = tphigrd(inrm,jln)
      tui  = tphigrd(inrm,iln)

      phj  = uj*ylm(jlm,ifgd)/rR
      phi  = ui*ylm(ilm,ifgd)/rR
      tphj = tuj*ylm(jlm,ifgd)/rR
      tphi = tui*ylm(ilm,ifgd)/rR

!     here is the diamagnetic scalar field at the point r'
!     ----
      gipaw_aug(iatom)%dia(ifgd,klmn) = phj*phi - tphj*tphi

!     dj means d (u*S/r) dx_idir = x*S*(du/dr)/r^2 - u*x*S/r^3 + u*(dS/dx)/r
!     ---- 
      do idir = 1, 3 ! loop over three components of para field
       cr = pawfgrtab(iatom)%rfgd(idir,ifgd) ! component of |r-R| in idir direction
       dj  = (uj*ylmgr(idir,jlm,ifgd)+cr*ylm(jlm,ifgd)*(dphigrd(inrm,jln)-uj/rR)/rR)/rR
       di  = (ui*ylmgr(idir,ilm,ifgd)+cr*ylm(ilm,ifgd)*(dphigrd(inrm,iln)-ui/rR)/rR)/rR
       dtj = (tuj*ylmgr(idir,jlm,ifgd)+cr*ylm(jlm,ifgd)*(dtphigrd(inrm,jln)-uj/rR)/rR)/rR
       dti = (tui*ylmgr(idir,ilm,ifgd)+cr*ylm(ilm,ifgd)*(dtphigrd(inrm,iln)-ui/rR)/rR)/rR

!      ! NOTE: the actual current vector field is (-i) times the following. Because it is pure
!      ! imaginary, we leave it in this form to save space (avoid having to save a full complex number)
!      ----
       gipaw_aug(iatom)%para(idir,ifgd,klmn)=0.5*(dj*phi-phj*di-dtj*tphi+tphj*dti)
      end do ! end loop over three comonents of para field

     end if ! avoid |r-R| = 0

    end do ! end loop over nfgd

!   the <L> augmentation field consists of <phi_j|L|phi_i> - <tphi_j|L|tphi_i> for each component of L

!   Computation of <phi_j|phi_i>- <tphi_j|tphi_i>
!   ----
    ff(2:mesh_size)=pawtab(itypat)%phiphj(2:mesh_size,kln)-pawtab(itypat)%tphitphj(2:mesh_size,kln)
    call deducer0(ff,mesh_size,pawrad(itypat))
    call simp_gen(qji,ff,pawrad(itypat))

    do idir = 1, 3 ! components of L
     call slxyzs(jl,jm,idir,il,im,sls_val)
     gipaw_aug(iatom)%onsiteangmom(1,idir,klmn) = real(qji*sls_val)
     gipaw_aug(iatom)%onsiteangmom(2,idir,klmn) = aimag(qji*sls_val)
!    !!!     write(6,'(5i4,2f12.8)')jl,jm,idir,il,im,gipaw_aug(iatom)%onsiteangmom(1,idir,klmn),gipaw_aug(iatom)%onsiteangmom(2,idir,klmn)
    end do

   end do ! end loop over ilmn atomic basis states
  end do ! end loop over jlmn atomic basis states

  deallocate(dphigrd,dtphigrd,ff,nrm,nrm_ifftsph,phigrd,tphigrd,ylm,ylmgr)
 end do     ! Loop on atoms

!DEBUG
!write(6,*)' gipaw_aug_fields : exit '
!stop
!ENDDEBUG

 end subroutine gipaw_aug_fields
!!***
