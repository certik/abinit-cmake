!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmknhat3
!! NAME
!! pawmknhat3
!!
!! FUNCTION
!! PAW only:
!! Compute 1st-order compensation charge density (and derivatives) on the fine FFT grid
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  idir=direction of atomic displacement (in case of atomic displ. perturb.)
!!  ipert=index of perturbation
!!  izero=if 1, unbalanced components of nhat1(g) have to be set to zero
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nfft=number of point on the rectangular fft grid
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat1dim=size of pawnhat1 (number of FFT points)
!!           nhat1dim=(number of points in PAW sphere), if atomic displ. perturb.,
!!      else nhat1dim=nfft
!!  nrhoij1=size of pawrhoij1 (1 if atomic displ. perturb., else natom)
!!  ntypat=number of types of atoms in unit cell.
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawrhoij1(nrhoij1) <type(pawrhoij_type)>= derivatives of paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type (integer) for each atom
!!
!! OUTPUT
!!    pawnhat1(cplex*nhat1dim,ispden)=nhat(1) on fine rectangular grid
!!
!! PARENTS
!!      scfcv3
!!
!! CHILDREN
!!      fourdp,leave_new,wrtout,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawmknhat3(cplex,idir,ipert,izero,mpi_enreg,natom,nfft,ngfft,nhat1dim,nrhoij1,nspden,&
&          ntypat,paral_kgb,pawang,pawfgrtab,pawnhat1,pawrhoij,pawrhoij1,pawtab,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_13paw, except_this_one => pawmknhat3
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,izero,natom,nfft,nhat1dim,nrhoij1,nspden,ntypat,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18),typat(natom)
 real(dp),intent(out) :: pawnhat1(cplex*nhat1dim,nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom),pawrhoij1(nrhoij1)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!scalars
 integer :: iatom,iatom1,iatom2,ic,ils,ilslm,irhoij,ispden,itypat,jc,klm,klmn
 integer :: lmax,lmin,lm_size,mm,nfftot,optgr0,optgr1
 character(len=500) :: message
!arrays
 real(dp) :: rdum(1),ro(cplex),ro_ql(cplex)
 real(dp),allocatable :: work(:,:)

! *************************************************************************

!Compatibility tests
 if(ipert<=natom.and.nhat1dim/=pawfgrtab(iatom)%nfgd) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawmknhat3 : BUG -',ch10,&
&  '  Wrong size for PAW nhat(1) !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(ipert<=natom.and.nrhoij1/=1) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawmknhat3 : BUG -',ch10,&
&  '  Wrong size for rhoij(1) !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if
 if(cplex/=pawrhoij1(1)%cplex) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawmknhat3 : BUG -',ch10,&
&  '  Wrong value for cplex !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

 pawnhat1=zero

 if (ipert<=natom) then
  iatom1=ipert;iatom2=ipert
 else
  iatom1=1;iatom2=natom
 end if

!----- Loop over atoms
 do iatom=iatom1,iatom2
  itypat=typat(iatom)
  lm_size=pawfgrtab(iatom)%l_size**2

! Eventually compute g_l(r).Y_lm(r) factors for the current atom (if not already done)
  if ((pawfgrtab(iatom)%gylm_allocated==0).or.&
&  (pawfgrtab(iatom)%gylmgr_allocated==0)) then
   optgr0=0;optgr1=0
   if (pawfgrtab(iatom)%gylm_allocated==0) then
    if (associated(pawfgrtab(iatom)%gylm)) deallocate(pawfgrtab(iatom)%gylm)
    allocate(pawfgrtab(iatom)%gylm(pawfgrtab(iatom)%nfgd,pawfgrtab(iatom)%l_size**2))
    pawfgrtab(iatom)%gylm_allocated=2;optgr0=1
   end if
   if (pawfgrtab(iatom)%gylmgr_allocated==0) then
    if (associated(pawfgrtab(iatom)%gylmgr)) deallocate(pawfgrtab(iatom)%gylmgr)
    allocate(pawfgrtab(iatom)%gylmgr(3,pawfgrtab(iatom)%nfgd,pawfgrtab(iatom)%l_size**2))
    pawfgrtab(iatom)%gylmgr_allocated=2;optgr1=1
   end if
   call pawgylm(pawfgrtab(iatom)%gylm,pawfgrtab(iatom)%gylmgr,rdum,iatom,&
&   pawfgrtab(iatom)%ifftsph,itypat,lm_size,pawfgrtab(iatom)%nfgd,&
&   optgr0,optgr1,0,pawtab(itypat),pawfgrtab(iatom)%rfgd,&
&   pawfgrtab(iatom)%rfgd_allocated)
  end if

! ----- Loop over density components
  do ispden=1,nspden

!  A- Atomic displ. perturbation (phonons)
!  ---------------------------------------
   if (ipert<=natom) then

!   Contribution from g_l(r).Y_lm(r) derivatives
    do irhoij=1,pawrhoij(iatom)%nrhoijsel
     klmn=pawrhoij(iatom)%rhoijselect(irhoij)
     klm =pawtab(itypat)%indklmn(1,klmn)
     lmin=pawtab(itypat)%indklmn(3,klmn)
     lmax=pawtab(itypat)%indklmn(4,klmn)

!    Retrieve rhoij
     if (nspden/=2) then
      ro(1)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
     else
      if (ispden==1) then
       ro(1)=pawrhoij(iatom)%rhoijp(irhoij,1)&
&       +pawrhoij(iatom)%rhoijp(irhoij,2)
      else if (ispden==2) then
       ro(1)=pawrhoij(iatom)%rhoijp(irhoij,1)
      end if
     end if
     ro(1)=pawtab(itypat)%dltij(klmn)*ro(1)

     if (cplex==1) then
      do ils=lmin,lmax,2
       do mm=-ils,ils
        ilslm=ils*ils+ils+mm+1
        if (pawang%gntselect(ilslm,klm)>0) then
         ro_ql(1)=ro(1)*pawtab(itypat)%qijl(ilslm,klmn)
         do ic=1,pawfgrtab(iatom)%nfgd
          pawnhat1(ic,ispden)=pawnhat1(ic,ispden)-ro_ql(1)*pawfgrtab(iatom)%gylmgr(idir,ic,ilslm)
         end do
        end if
       end do
      end do
     else
      do ils=lmin,lmax,2
       do mm=-ils,ils
        ilslm=ils*ils+ils+mm+1
        if (pawang%gntselect(ilslm,klm)>0) then
         ro_ql(1)=ro(1)*pawtab(itypat)%qijl(ilslm,klmn)
         do ic=1,pawfgrtab(iatom)%nfgd
          jc=2*ic-1
          pawnhat1(jc,ispden)=pawnhat1(jc,ispden)-ro_ql(1)*pawfgrtab(iatom)%gylmgr(idir,ic,ilslm)
         end do
        end if
       end do
      end do
     end if

    end do ! irhoij

!   Contribution from rho_ij derivatives
    do irhoij=1,pawrhoij1(1)%nrhoijsel
     klmn=pawrhoij1(1)%rhoijselect(irhoij)
     klm =pawtab(itypat)%indklmn(1,klmn)
     lmin=pawtab(itypat)%indklmn(3,klmn)
     lmax=pawtab(itypat)%indklmn(4,klmn)

!    Retrieve rhoij
     if (nspden/=2) then
      ro(1:cplex)=pawrhoij1(1)%rhoijp(cplex*(irhoij-1)+1:cplex*irhoij,ispden)
     else
      if (ispden==1) then
       ro(1:cplex)=pawrhoij1(1)%rhoijp(cplex*(irhoij-1)+1:cplex*irhoij,1) &
&       +pawrhoij1(1)%rhoijp(cplex*(irhoij-1)+1:cplex*irhoij,2)
      else if (ispden==2) then
       ro(1:cplex)=pawrhoij1(1)%rhoijp(cplex*(irhoij-1)+1:cplex*irhoij,1)
      end if
     end if
     ro(1:cplex)=pawtab(itypat)%dltij(klmn)*ro(1:cplex)

     if (cplex==1) then
      do ils=lmin,lmax,2
       do mm=-ils,ils
        ilslm=ils*ils+ils+mm+1
        if (pawang%gntselect(ilslm,klm)>0) then
         ro_ql(1:cplex)=ro(1:cplex)*pawtab(itypat)%qijl(ilslm,klmn)
         do ic=1,pawfgrtab(iatom)%nfgd
          pawnhat1(ic,ispden)=pawnhat1(ic,ispden)+ro_ql(1)*pawfgrtab(iatom)%gylm(ic,ilslm)
         end do
        end if
       end do
      end do
     else
      do ils=lmin,lmax,2
       do mm=-ils,ils
        ilslm=ils*ils+ils+mm+1
        if (pawang%gntselect(ilslm,klm)>0) then
         do ic=1,pawfgrtab(iatom)%nfgd
          pawnhat1(2*ic-1:2*ic,ispden)=pawnhat1(2*ic-1:2*ic,ispden)+ro_ql(1:2)*pawfgrtab(iatom)%gylm(ic,ilslm)
         end do
        end if
       end do
      end do
     end if

    end do ! irhoij

   else

!   B- Other perturbations
!   ---------------------------------------
    write(message, '(a,a,a,a)' )ch10,&
&    ' pawmknhat3 : BUG -',ch10,&
&    '  ipert>natom not yet implemented !'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')

!   Contribution from g_l(r).Y_lm(r) derivatives
!   To be implemented

!   Contribution from rho_ij derivatives
!   do irhoij=1,pawrhoij1(iatom)%nrhoijsel
!   klmn=pawrhoij1(iatom)%rhoijselect(irhoij)
!   klm =pawtab(itypat)%indklmn(1,klmn)
!   lmin=pawtab(itypat)%indklmn(3,klmn)
!   lmax=pawtab(itypat)%indklmn(4,klmn)
!   ro(1:cplex)=pawtab(itypat)%dltij(klmn)*pawrhoij1(iatom)%rhoijp(cplex*(irhoij-1)+1:cplex*irhoij,ispden)
!   do ils=lmin,lmax,2
!   do mm=-ils,ils
!   ilslm=ils*ils+ils+mm+1
!   if (pawang%gntselect(ilslm,klm)>0) then
!   ro_ql(1:cplex)=ro(1:cplex)*pawtab(itypat)%qijl(ilslm,klmn)
!   if (cplex==1) then
!   do ic=1,pawfgrtab(iatom)%nfgd
!   jc=pawfgrtab(iatom)%ifftsph(ic)
!   pawnhat1(jc,ispden)=pawnhat1(jc,ispden)+ro_ql(1)*pawfgrtab(iatom)%gylm(ic,ilslm)
!   end do
!   else
!   do ic=1,pawfgrtab(iatom)%nfgd
!   jc=pawfgrtab(iatom)%ifftsph(ic)
!   pawnhat1(2*jc-1:2*jc,ispden)=pawnhat1(2*jc-1,2*jc,ispden)+ro_ql(1:2)*pawfgrtab(iatom)%gylm(ic,ilslm)
!   end do
!   end if
!   end if
!   end do
!   end do
!   end do
   end if

!  ----- End loop over density components
  end do

! Temporary storage deallocation
  if (pawfgrtab(iatom)%gylm_allocated==2) then
   deallocate(pawfgrtab(iatom)%gylm);allocate(pawfgrtab(iatom)%gylm(0,0))
   pawfgrtab(iatom)%gylm_allocated=0
  end if
  if (pawfgrtab(iatom)%gylmgr_allocated==2) then
   deallocate(pawfgrtab(iatom)%gylmgr);allocate(pawfgrtab(iatom)%gylmgr(0,0,0))
   pawfgrtab(iatom)%gylmgr_allocated=0
  end if

! ----- End loop over atoms
 end do

!----- Avoid unbalanced g-components numerical errors
 if (izero==1.and.nhat1dim==nfft) then
  allocate(work(2,nfft))
  do ispden=1,min(2,nspden)
   call fourdp(cplex,work,pawnhat1(:,ispden),-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   call zerosym(work,2,mpi_enreg,ngfft(1),ngfft(2),ngfft(3))
   call fourdp(cplex,work,pawnhat1(:,ispden),+1,mpi_enreg,nfft,ngfft,paral_kgb,0)
  end do
  deallocate(work)
 end if

end subroutine pawmknhat3
!!***
