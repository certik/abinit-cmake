!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmkrhoij3
!!
!! NAME
!! pawmkrhoij3
!!
!! FUNCTION
!! Calculate the 1st-order PAW quantities rhoij1 (1st-order augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  cprj(dimpaw1,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!                         projected with non-local projectors: cprj_nk(i)=<p_i|Cnk>
!!  cprj1(dimpaw1,nspinor*mband*mk1mem*nsppol)= 1st-order wave functions at k,q
!!                         projected with non-local projectors: cprj1_nkq(i)=<p_i|C1nk,q>
!!  dimcprj=array of dimensions of arrays cprj,cprj1
!!  dimpaw1=size of cprj, cprj1, pawrhoij1 (1 if atomic displ. perturb., else natom)
!!  ipert=index of perturbation
!!  mband=maximum number of bands
!!  mkmem=number of k points which can fit in memory (GS data)  ; 0 if use disk
!!  mk1mem=number of k points which can fit in memory (RF data); 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)= # atoms of each type.
!!  nband=number of bands for all k points
!!  nkpt=number of k points
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k
!!  pawprtvol=control print volume and debugging output for PAW
!!  unpaw=unit number for cprj data (if used)
!!  unpaw1=unit number for cprj1 data (if used)
!!  usecprj= 1 if cprj array is stored in memory
!!  wtk(nkpt)=weight assigned to each k point
!!
!! SIDE EFFECTS
!!  pawrhoij1(dimpaw1) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!  On input: arrays dimensions
!!  On output:
!!    pawrhoij1(ii)%rhoij_(lmn2_size,nspden)=      (non symetrized)
!!            Sum_{n,k} {occ(n,k)*(conjugate[cprj_nk(ii)].cprj1_nkq(jj)
!!                                 conjugate[cprj_nk(jj)].cprj1_nkq(ii)}
!!          + Sum_{n,k} {occ(n,k)*(conjugate[dcprj_nk(ii)/dlambda].cprj_nk(jj)
!!                                +conjugate[cprj_nk(ii)].dcprj_nk(jj)/dlambda)}
!!
!! PARENTS
!!      vtorho3
!!
!! CHILDREN
!!      cprj_diskinit,cprj_get,leave_new,leave_test,print_ij,timab,wrtout,xcomm_init,xme_init,xsum_mpi
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pawmkrhoij3(atindx1,cplex,cprj,cprj1,dimcprj,dimpaw1,ipert,istwfk,mband,mkmem,mk1mem,&
&                       mpi_enreg,natom,nattyp,nband,nkpt,nspden,nspinor,nsppol,ntypat,occ,pawprtvol,&
&                       pawrhoij1,unpaw,unpaw1,usecprj,wtk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: cplex,dimpaw1,ipert,mband,mkmem,mk1mem,natom,nkpt
 integer,intent(in) :: nspden,nspinor,nsppol,ntypat,pawprtvol,unpaw,unpaw1,usecprj
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom),istwfk(nkpt)
 integer,intent(in) :: nattyp(ntypat),nband(nkpt*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),wtk(nkpt)
 type(cprj_type),intent(in) :: cprj(dimpaw1,nspinor*mband*mkmem*nsppol*usecprj)
 type(cprj_type),intent(in) :: cprj1(dimpaw1,nspinor*mband*mk1mem*nsppol)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(dimpaw1)

!Local variables ---------------------------------------
!scalars
 integer :: bdtot_index,bufdim,iatom,ib,iband,ibg,ibg1,ierr,ikpt,iorder_cprj,isppol
 integer :: jdim,me,natinc,nband_k,ncpgr,nsp2,option,spaceComm
 real(dp) :: wtk_k
 character(len=500) :: message
!arrays
 integer,allocatable :: dimlmn(:),idum(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 character(len=8),parameter :: dspin(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 type(cprj_type),allocatable :: cwaveprj(:,:),cwaveprj1(:,:)

!************************************************************************

!Tests
 if(usecprj==0) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' pawmkrhoij3 : ERROR -',ch10,&
&  '  Not allowed for usecprj=0  !'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Init MPI data
 call xcomm_init(mpi_enreg,spaceComm)
 call xme_init(mpi_enreg,me)

!Initialize temporary files (if used)
 iorder_cprj=0
 call cprj_diskinit(atindx1,dimpaw1,iorder_cprj,mkmem,natom,dimcprj,nspinor,unpaw)
 call cprj_diskinit(atindx1,dimpaw1,iorder_cprj,mk1mem,natom,dimcprj,nspinor,unpaw1)

!Allocate temporary cprj storage
 ncpgr=0;if (ipert<=natom) ncpgr=1
 allocate(cwaveprj(dimpaw1,nspinor),cwaveprj1(dimpaw1,nspinor))
 call cprj_alloc(cwaveprj,ncpgr,dimcprj)
 call cprj_alloc(cwaveprj1,ncpgr,dimcprj)

!Initialize output quantities
 do iatom=1,dimpaw1
  pawrhoij1(iatom)%rhoij_=zero
 end do

!LOOP OVER SPINS
 option=2
 bdtot_index=0;ibg=0;ibg1=0
 do isppol=1,nsppol

! LOOP OVER k POINTS
  do ikpt=1,nkpt

   nband_k=nband(ikpt+(isppol-1)*nkpt)
   wtk_k=wtk(ikpt)

   if(mpi_enreg%paral_compil_kpt==1)then
    if(minval(abs(mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol)-me))/=0) then
     bdtot_index=bdtot_index+nband_k
     cycle
    end if
   end if

!  LOOP OVER BANDS
   do ib=1,nband_k
    iband=bdtot_index+ib

    if(mpi_enreg%paral_compil_kpt==1)then
     if (mpi_enreg%proc_distrb(ikpt,ib,isppol)/= me) cycle
    end if

!   Extract cprj for current band
!   Must read cprj when mkmem=0 (even if unused) to have right pointer inside _PAW file
    if (abs(occ(iband))>tol8.or.mkmem==0) &
&    call cprj_get(atindx1,cwaveprj,cprj,dimpaw1,ib,ibg,ikpt,iorder_cprj,isppol,&
&    mband,mkmem,mpi_enreg,natom,1,nband_k,nspinor,nsppol,unpaw)

!   Extract cprj1 for current band
    if (abs(occ(iband))>tol8.or.mk1mem==0) &
&    call cprj_get(atindx1,cwaveprj1,cprj1,dimpaw1,ib,ibg1,ikpt,iorder_cprj,isppol,&
&    mband,mk1mem,mpi_enreg,natom,1,nband_k,nspinor,nsppol,unpaw1)

!   Accumulate contribution from (occuppied) current band
    if (abs(occ(iband))>tol8) &
&    call pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj1,dimpaw1,ipert,isppol,natom,nspden,&
&    nspinor,nsppol,occ(iband),option,pawrhoij1,wtk_k)

   end do ! ib

   bdtot_index=bdtot_index+nband_k
   if (mkmem/=0)  ibg =ibg +nspinor*nband_k
   if (mk1mem/=0) ibg1=ibg1+nspinor*nband_k

  end do ! ikpt
 end do ! isppol

!deallocate temporary cwaveprj, cwaveprj1 storage
 call cprj_free(cwaveprj)
 call cprj_free(cwaveprj1)
 deallocate(cwaveprj,cwaveprj1)

!MPI: need to exchange arrays between procs
!==========================================
 if(mpi_enreg%paral_compil_kpt==1)then
  if (mpi_enreg%parareel == 0)call leave_test(mpi_enreg)
  call timab(48,1,tsec)

! Exchange pawrhoij1%rhoij_make

  call timab(48,1,tsec)
  allocate(dimlmn(dimpaw1))
  if (dimpaw1==1) then
   dimlmn(1)=dimcprj(ipert)
  else
   dimlmn(1:natom)=dimcprj(1:natom)
  end if
  nsp2=nsppol;if (nspden==4) nsp2=4
  bufdim=sum(dimlmn)*nsp2
  allocate(buffer1(bufdim),buffer2(bufdim))
  jdim=0
  do iatom=1,dimpaw1
   do isppol=1,nsp2
    buffer1(jdim+1:jdim+dimlmn(iatom))=pawrhoij1(iatom)%rhoij_(:,isppol)
    jdim=jdim+dimlmn(iatom)
   end do
  end do
  call xsum_mpi(buffer1,buffer2,bufdim,spaceComm,ierr) !Build sum of everything
  jdim=0
  do iatom=1,dimpaw1
   do isppol=1,nsp2
    pawrhoij1(iatom)%rhoij_(:,isppol)=buffer2(jdim+1:jdim+dimlmn(iatom))
    jdim=jdim+dimlmn(iatom)
   end do
  end do
  deallocate(buffer1,buffer2,dimlmn)
  call timab(48,2,tsec)
 end if ! mpi_enreg%paral_compil_kpt==1

!Print info
 if (abs(pawprtvol)>=1) then
  natinc=1;if(dimpaw1>1.and.pawprtvol>=0) natinc=dimpaw1-1
  nsp2=nsppol;if (nspden==4) nsp2=4
  do iatom=1,dimpaw1,natinc
   write(message, '(4a,i3,a)') ch10," PAW TEST:",ch10,&
&   ' ====== Values of RHOIJ(1) in pawmkrhoij3 (index=',iatom,') ======'
   if (nspden==2.and.nsppol==1) write(message,'(3a)') trim(message),ch10,&
&   '      (antiferromagnetism case: only one spin component)'
   call wrtout(6,message,'COLL')
   do isppol=1,nsp2
    if (nspden/=1) write(message, '(3a)') '   Component ',trim(dspin(isppol+2*(nspden/4))),':'
    call wrtout(6,message,'COLL')
    call print_ij(pawrhoij1(iatom)%rhoij_(:,isppol),pawrhoij1(iatom)%lmn2_size,&
&    cplex,pawrhoij1(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
   end do
  end do
 end if

end subroutine pawmkrhoij3
!!***
