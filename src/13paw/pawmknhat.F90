!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmknhat
!! NAME
!! pawmknhat
!!
!! FUNCTION
!! PAW only:
!! Compute compensation charge density (and derivatives) on the fine FFT grid
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  ider= 0: nhat(r) is computed
!!        1: cartesian derivatives of nhat(r) are computed
!!        2: nhat(r) and derivatives are computed
!!  izero=if 1, unbalanced components of nhat(g) have to be set to zero
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=number of point on the rectangular fft grid
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhatgrdim= -PAW only- 0 if pawgrnhat array is not used ; 1 otherwise
!!  ntypat=number of types of atoms in unit cell.
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type (integer) for each atom
!!  ucvol=volume of the unit cell
!!
!! OUTPUT
!!  === if ider=0 or 2
!!    compch_fft=compensation charge inside spheres computed over fine fft grid
!!    pawnhat(nfft,ispden)=nhat on fine rectangular grid
!!  === if ider=1 or 2
!!    pawgrnhat(nfft,ispden,3)=derivatives of nhat on fine rectangular grid (and derivatives)
!!
!! PARENTS
!!      afterscfloop,energy,outscfcv,scfcv,vtorho
!!
!! CHILDREN
!!      fourdp,leave_new,mean_fftr,wrtout,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawmknhat(compch_fft,ider,izero,mpi_enreg,natom,nfft,ngfft,nhatgrdim,nspden,&
&          ntypat,paral_kgb,pawang,pawfgrtab,pawgrnhat,pawnhat,pawrhoij,pawtab,typat,ucvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_12spacepar
 use interfaces_13paw, except_this_one => pawmknhat
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ider,izero,natom,nfft,nhatgrdim,nspden,ntypat,paral_kgb
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: compch_fft
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(out) :: pawgrnhat(nfft,nspden,3*nhatgrdim)
 real(dp),intent(out) :: pawnhat(nfft,nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: iatom,ic,ils,ilslm,irhoij,ispden,itypat,jc,klm,klmn
 integer :: lmax,lmin,lm_size,mm,nfftot,optgr0,optgr1
 real(dp) :: ro,ro_ql
 logical :: compute_grad,compute_nhat
 character(len=500) :: message
!arrays
 real(dp) :: rdum(1),tmp_compch_fft(nspden)
 real(dp),allocatable :: work(:,:)

! *************************************************************************

 if(ider>0.and.nhatgrdim==0) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawmknhat : BUG -',ch10,&
&  '  Gradients of nhat required but not allocated !'
  call wrtout(6,message,'COLL')
  call leave_new('PERS')
 end if
 compute_nhat=(ider==0.or.ider==2)
 compute_grad=(ider==1.or.ider==2)
 if ((.not.compute_nhat).and.(.not.compute_grad)) return

 if (compute_nhat) pawnhat=zero
 if (compute_grad) pawgrnhat=zero

!------------------------------------------------------------------------
!----- Loop over atoms
!------------------------------------------------------------------------
 do iatom=1,natom
  itypat=typat(iatom)
  lm_size=pawfgrtab(iatom)%l_size**2

! Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
  if (((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)).or.&
&  ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0))) then
   optgr0=0;optgr1=0
   if ((compute_nhat).and.(pawfgrtab(iatom)%gylm_allocated==0)) then
    if (associated(pawfgrtab(iatom)%gylm)) deallocate(pawfgrtab(iatom)%gylm)
    allocate(pawfgrtab(iatom)%gylm(pawfgrtab(iatom)%nfgd,pawfgrtab(iatom)%l_size**2))
    pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
   end if
   if ((compute_grad).and.(pawfgrtab(iatom)%gylmgr_allocated==0)) then
    if (associated(pawfgrtab(iatom)%gylmgr)) deallocate(pawfgrtab(iatom)%gylmgr)
    allocate(pawfgrtab(iatom)%gylmgr(3,pawfgrtab(iatom)%nfgd,pawfgrtab(iatom)%l_size**2))
    pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
   end if
   call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,rdum,iatom,&
&   pawfgrtab(iatom)%ifftsph,itypat,lm_size,pawfgrtab(iatom)%nfgd,&
&   optgr0,optgr1,0,pawtab(itypat),pawfgrtab(iatom)%rfgd,&
&   pawfgrtab(iatom)%rfgd_allocated)
  end if

! ------------------------------------------------------------------------
! ----- Loop over density components
! ------------------------------------------------------------------------

  do ispden=1,nspden

!  ------------------------------------------------------------------------
!  ----- Loop over ij channels (basis components)
!  ------------------------------------------------------------------------
   do irhoij=1,pawrhoij(iatom)%nrhoijsel
    klmn=pawrhoij(iatom)%rhoijselect(irhoij)
    klm =pawtab(itypat)%indklmn(1,klmn)
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

    if (compute_nhat) then
     do ils=lmin,lmax,2
      do mm=-ils,ils
       ilslm=ils*ils+ils+mm+1
       if (pawang%gntselect(ilslm,klm)>0) then
        ro_ql=ro*pawtab(itypat)%qijl(ilslm,klmn)
        do ic=1,pawfgrtab(iatom)%nfgd
         jc=pawfgrtab(iatom)%ifftsph(ic)
         pawnhat(jc,ispden)=pawnhat(jc,ispden)+ro_ql*pawfgrtab(iatom)%gylm(ic,ilslm)
        end do
       end if
      end do
     end do
    end if

    if (compute_grad) then
     do ils=lmin,lmax,2
      do mm=-ils,ils
       ilslm=ils*ils+ils+mm+1
       if (pawang%gntselect(ilslm,klm)>0) then
        ro_ql=ro*pawtab(itypat)%qijl(ilslm,klmn)
        do ic=1,pawfgrtab(iatom)%nfgd
         jc=pawfgrtab(iatom)%ifftsph(ic)
         pawgrnhat(jc,ispden,1)=pawgrnhat(jc,ispden,1)+ro_ql*pawfgrtab(iatom)%gylmgr(1,ic,ilslm)
         pawgrnhat(jc,ispden,2)=pawgrnhat(jc,ispden,2)+ro_ql*pawfgrtab(iatom)%gylmgr(2,ic,ilslm)
         pawgrnhat(jc,ispden,3)=pawgrnhat(jc,ispden,3)+ro_ql*pawfgrtab(iatom)%gylmgr(3,ic,ilslm)
        end do
       end if
      end do
     end do
    end if

!   ------------------------------------------------------------------------
!   ----- End loop over ij channels
!   ------------------------------------------------------------------------
   end do

!  ------------------------------------------------------------------------
!  ----- End loop over density components
!  ------------------------------------------------------------------------
  end do

  if (pawfgrtab(iatom)%gylm_allocated==2) then
   deallocate(pawfgrtab(iatom)%gylm);allocate(pawfgrtab(iatom)%gylm(0,0))
   pawfgrtab(iatom)%gylm_allocated=0
  end if
  if (pawfgrtab(iatom)%gylmgr_allocated==2) then
   deallocate(pawfgrtab(iatom)%gylmgr);allocate(pawfgrtab(iatom)%gylmgr(0,0,0))
   pawfgrtab(iatom)%gylmgr_allocated=0
  end if

! ------------------------------------------------------------------------
! ----- End loop over atoms
! ------------------------------------------------------------------------
 end do

!----- Avoid unbalanced g-components numerical errors
 if (izero==1.and.compute_nhat) then
  allocate(work(2,nfft))
  do ispden=1,min(2,nspden)
   call fourdp(1,work,pawnhat(:,ispden),-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   call zerosym(work,2,mpi_enreg,ngfft(1),ngfft(2),ngfft(3))
   call fourdp(1,work,pawnhat(:,ispden),+1,mpi_enreg,nfft,ngfft,paral_kgb,0)
  end do
  deallocate(work)
 end if

!----- Computation of compensation charge over FFT grid
 if (compute_nhat) then
  nfftot=ngfft(1)*ngfft(2)*ngfft(3)
  call mean_fftr(pawnhat,tmp_compch_fft,mpi_enreg,nfft,nfftot,1)
  compch_fft = tmp_compch_fft(1)
  compch_fft=compch_fft*ucvol
 end if

end subroutine pawmknhat
!!***
