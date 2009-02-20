!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawmkrhoij
!!
!! NAME
!! pawmkrhoij
!!
!! FUNCTION
!! Calculate the PAW quantities rhoij (augmentation occupancies)
!! Remember:for each atom, rho_ij=Sum_{n,k} {occ(n,k)*<Cnk|p_i><p_j|Cnk>}
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cprj(natom,nspinor*mband*mkmem*nsppol)= wave functions projected with non-local projectors:
!!                                   cprj_nk(i)=<p_i|Cnk> where p_i is a non-local projector.
!!  dimcprj=array of dimensions of array cprj
!!  mband=maximum number of bands
!!  mkmem =number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nattyp(ntypat)= # atoms of each type.
!!  nband=number of bands for all k points
!!  nkpt=number of k points
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ(mband*nkpt*nsppol)=occupation number for each band for each k
!!  pawprtvol=control print volume and debugging output for PAW
!!  unpaw=unit number for cprj PAW data (if used)
!!  wtk(nkpt)=weight assigned to each k point
!!
!! SIDE EFFECTS
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  On input: arrays dimensions
!!  On output:
!!    pawrhoij(:)%rhoij_(lmn2_size,nspden)=
!!          Sum_{n,k} {occ(n,k)*conjugate[cprj_nk(ii)].cprj_nk(jj)} (non symetrized)
!!
!! PARENTS
!!      vtorho
!!
!! CHILDREN
!!      cprj_alloc,cprj_diskinit,cprj_free,cprj_get,leave_new,leave_test,print_ij,timab,wrtout,xcomm_init,xme_init,xsum_mpi
!!
!! OUTPUT
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pawmkrhoij(atindx1,cprj,dimcprj,istwfk,mband,mkmem,mpi_enreg,natom,nattyp,&
&                      nband,nkpt,nspden,nspinor,nsppol,ntypat,occ,pawprtvol,pawrhoij,unpaw,wtk)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13nonlocal
 use interfaces_13paw, except_this_one => pawmkrhoij
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,natom,nkpt,nspden,nspinor,nsppol,ntypat,pawprtvol
 integer,intent(in) :: unpaw
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),dimcprj(natom),istwfk(nkpt)
 integer,intent(in) :: nattyp(ntypat),nband(nkpt*nsppol)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),wtk(nkpt)
 type(cprj_type),intent(in) :: cprj(natom,nspinor*mband*mkmem*nsppol)
 type(pawrhoij_type),intent(inout) :: pawrhoij(natom)

!Local variables ---------------------------------------
!scalars
 integer :: bdtot_index,bufdim,cplex,iatom,ib,iband,ibg,ierr,ikpt,iorder_cprj,isppol
 integer :: jdim,me,natinc,nband_k,nsp2,option,spaceComm
 real(dp) :: wtk_k
 character(len=500) :: message
!arrays
 integer,allocatable :: dimlmn(:),idum(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:)
 character(len=8),parameter :: dspin(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 type(cprj_type),allocatable :: cwaveprj(:,:)

!************************************************************************

!DEBUG
!write(6,*)' pawmkrhoij : enter '
!ENDDEBUG

!Init MPI data
 if ((mpi_enreg%paral_compil_kpt==1) .and. &
& (mpi_enreg%paral_compil_fft==1)) then
  me        = mpi_enreg%me_kpt
  spaceComm = mpi_enreg%comm_kpt
 else
  call xcomm_init(mpi_enreg,spaceComm)
  call xme_init(mpi_enreg,me)
 end if

!Initialize temporary file (if used)
 iorder_cprj=0
 call cprj_diskinit(atindx1,natom,iorder_cprj,mkmem,natom,dimcprj,nspinor,unpaw)

!Allocate temporary cwaveprj storage
 allocate(cwaveprj(natom,nspinor))
 call cprj_alloc(cwaveprj,0,dimcprj)

!Initialize output quantities
 do iatom=1,natom
  pawrhoij(iatom)%rhoij_=zero
 end do

!LOOP OVER SPINS
 option=1
 bdtot_index=0;ibg=0
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

   cplex=2;if (istwfk(ikpt)>1) cplex=1

!  LOOP OVER BANDS
   do ib=1,nband_k
    iband=bdtot_index+ib

    if(mpi_enreg%paral_compil_kpt==1)then
     if (mpi_enreg%proc_distrb(ikpt,ib,isppol)/=me) cycle
    end if

!   Extract cprj for current band
!   Must read cprj when mkmem=0 (even if unused) to have right pointer inside _PAW file
    if (abs(occ(iband))>tol8.or.mkmem==0) then
     call cprj_get(atindx1,cwaveprj,cprj,natom,ib,ibg,ikpt,iorder_cprj,isppol,&
&     mband,mkmem,mpi_enreg,natom,1,nband_k,nspinor,nsppol,unpaw)
    end if

!   Accumulate contribution from (occupied) current band
    if (abs(occ(iband))>tol8) then
     call pawaccrhoij(atindx1,cplex,cwaveprj,cwaveprj,natom,0,isppol,natom,nspden,&
&     nspinor,nsppol,occ(iband),option,pawrhoij,wtk_k)
    end if

   end do ! ib

   bdtot_index=bdtot_index+nband_k
   if (mkmem/=0) ibg=ibg+nspinor*nband_k

  end do ! ikpt
 end do ! isppol

!deallocate temporary cwaveprj storage
 call cprj_free(cwaveprj)
 deallocate(cwaveprj)

!MPI: need to exchange arrays between procs
!==========================================
 if(mpi_enreg%paral_compil_kpt==1)then
  if (mpi_enreg%parareel == 0) call leave_test(mpi_enreg)

! Exchange pawrhoij%rhoij_
  call timab(48,1,tsec)
  allocate(dimlmn(natom))
  dimlmn(1:natom)=pawrhoij(1:natom)%lmn2_size
  nsp2=nsppol;if (nspden==4) nsp2=4
  bufdim=sum(dimlmn)*nsp2
  allocate(buffer1(bufdim),buffer2(bufdim))
  jdim=0
  do iatom=1,natom
   do isppol=1,nsp2
    buffer1(jdim+1:jdim+dimlmn(iatom))=pawrhoij(iatom)%rhoij_(:,isppol)
    jdim=jdim+dimlmn(iatom)
   end do
  end do
  call xsum_mpi(buffer1,buffer2,bufdim,spaceComm,ierr) !Build sum of everything
  jdim=0
  do iatom=1,natom
   do isppol=1,nsp2
    pawrhoij(iatom)%rhoij_(:,isppol)=buffer2(jdim+1:jdim+dimlmn(iatom))
    jdim=jdim+dimlmn(iatom)
   end do
  end do
  deallocate(buffer1,buffer2,dimlmn)
  call timab(48,2,tsec)
 end if ! mpi_enreg%paral_compil_kpt==1

!Print info
 if (abs(pawprtvol)>=1) then
  natinc=1;if(natom>1.and.pawprtvol>=0) natinc=natom-1
  nsp2=nsppol;if (nspden==4) nsp2=4
  do iatom=1,natom,natinc
   write(message, '(4a,i3,a)') ch10," PAW TEST:",ch10,&
&   ' ====== Values of RHOIJ in pawmkrhoij (iatom=',iatom,') ======'
   if (nspden==2.and.nsppol==1) write(message,'(3a)') trim(message),ch10,&
&   '      (antiferromagnetism case: only one spin component)'
   call wrtout(6,message,'COLL')
   do isppol=1,nsp2
    if (nspden/=1) then
     write(message, '(3a)') '   Component ',trim(dspin(isppol+2*(nspden/4))),':'
     call wrtout(6,message,'COLL')
    end if
    call print_ij(pawrhoij(iatom)%rhoij_(:,isppol),pawrhoij(iatom)%lmn2_size,&
&    1,pawrhoij(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,-1.d0,1)
   end do
  end do
 end if

!DEBUG
!write(6,*)' pawmkrhoij : exit '
!ENDDEBUG

end subroutine pawmkrhoij
!!***
