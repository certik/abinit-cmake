!{\src2tex{textfont=tt}}
!! ===================================================
!! This module contains functions used to manipulate
!! variables of structured datatype cprj_type.
!! cprj_type variables are <p_lmn|Cnk> projected
!! quantities where |p_lmn> are non-local projectors
!!                  |Cnk> are wave functions
!! ===================================================

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

!!****f* ABINIT/cprj_alloc
!! NAME
!! cprj_alloc
!!
!! FUNCTION
!! Allocation of a cprj datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ncpgr=number of gradients to be allocated
!!  nlmn(:)=sizes of cprj%cp
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(cprj_type)>= cprj datastructure
!!
!! PARENTS
!!      calc_vHxc_braket,csigme,ctocprj,dyfnl3,energy,getgh1c,loper3,optics_paw
!!      outkss,partial_dos_fractions_paw,paw_symcprj,pawmkrhoij,pawmkrhoij3
!!      prep_nonlop,rdm,scfcv,screening,sigma,suscep_stat,vtorho,vtorho3,vtowfk
!!      vtowfk3
!!
!! CHILDREN
!!      xallgather_mpi,xme_init
!!
!! SOURCE

 subroutine cprj_alloc(cprj,ncpgr,nlmn)

 use defs_basis
 use defs_datatypes

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncpgr
!arrays
 integer,intent(in) :: nlmn(:)
 type(cprj_type),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim,nn

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2);nn=size(nlmn,dim=1)
 if (nn/=n1dim) stop "Error in cprj_alloc: wrong sizes !"
 do jj=1,n2dim
  do ii=1,n1dim
   nn=nlmn(ii)
   cprj(ii,jj)%nlmn=nn
   cprj(ii,jj)%ncpgr=ncpgr
   allocate(cprj(ii,jj)%cp(2,nn))
   if (ncpgr>0) allocate(cprj(ii,jj)%dcp(2,ncpgr,nn))
!  XG 080820 Was needed to get rid off problems with test paral#R with four procs
   cprj(ii,jj)%cp=zero
!  END XG 080820
  end do
 end do
end subroutine cprj_alloc
!!***

!!****f* ABINIT/cprj_free
!! NAME
!! cprj_free
!!
!! FUNCTION
!! Deallocation of a cprj datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(cprj_type)>= cprj datastructure
!!
!! PARENTS
!!      calc_vHxc_braket,calc_wf_qp,cchi0,cchi0q0,csigme,ctocprj,dyfnl3,energy
!!      get_bands_sym_GW,getgh1c,loper3,optics_paw,outkss
!!      partial_dos_fractions_paw,paw_symcprj,pawmkrhoij,pawmkrhoij3
!!      prep_nonlop,rdm,screening,sigma,suscep_stat,vtorho,vtorho3,vtowfk
!!      vtowfk3
!!
!! CHILDREN
!!      xallgather_mpi,xme_init
!!
!! SOURCE

 subroutine cprj_free(cprj)

 use defs_basis
 use defs_datatypes

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(cprj_type),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2)
 do jj=1,n2dim
  do ii=1,n1dim
   deallocate(cprj(ii,jj)%cp)
   if (cprj(ii,jj)%ncpgr>0) deallocate(cprj(ii,jj)%dcp)
  end do
 end do
end subroutine cprj_free
!!***

!!****f* ABINIT/cprj_nullify
!! NAME
!! cprj_nullify
!!
!! FUNCTION
!! Nullification of a cprj datastructure
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SIDE EFFECTS
!!  cprj(:,:) <type(cprj_type)>= cprj datastructure
!!
!! PARENTS
!!      prep_nonlop
!!
!! CHILDREN
!!      xallgather_mpi,xme_init
!!
!! SOURCE

 subroutine cprj_nullify(cprj)

 use defs_basis
 use defs_datatypes

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(cprj_type),intent(inout) :: cprj(:,:)
!Local variables-------------------------------
 integer :: ii,jj,n1dim,n2dim

! *************************************************************************

 n1dim=size(cprj,dim=1);n2dim=size(cprj,dim=2)
 do jj=1,n2dim
  do ii=1,n1dim
   cprj(ii,jj)%cp(:,:)=zero
   if (cprj(ii,jj)%ncpgr>0) cprj(ii,jj)%dcp(:,:,:)=zero
  end do
 end do
end subroutine cprj_nullify
!!***


!!****f* ABINIT/cprj_copy
!! NAME
!! cprj_copy
!!
!! FUNCTION
!! Copy a cprj datastructure into another
!!
!! COPYRIGHT
!! Copyright (C) 2008-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cprj_in(:,:) <type(cprj_type)>= input cprj datastructure
!!
!! OUTPUT
!!  cprj_out(:,:) <type(cprj_type)>= output cprj datastructure
!!
!! PARENTS
!!      dyfnl3,prep_nonlop,vtowfk3
!!
!! CHILDREN
!!      xallgather_mpi,xme_init
!!
!! SOURCE

 subroutine cprj_copy(cprj_in,cprj_out)

 use defs_basis
 use defs_datatypes

 implicit none
!Arguments ------------------------------------
!scalars
!arrays
 type(cprj_type),intent(in) :: cprj_in(:,:)
 type(cprj_type),intent(inout) :: cprj_out(:,:)
!Local variables-------------------------------
 integer :: ii,jj,kk,n1dim_in,n1dim_out,n2dim_in,n2dim_out,ncpgr_in,ncpgr_out,nlmn

! *************************************************************************

 n1dim_in=size(cprj_in,dim=1);n1dim_out=size(cprj_out,dim=1)
 n2dim_in=size(cprj_in,dim=2);n2dim_out=size(cprj_out,dim=2)
 ncpgr_in=cprj_in(1,1)%ncpgr;ncpgr_out=cprj_out(1,1)%ncpgr
 if (n1dim_in/=n1dim_out) stop "Error in cprj_copy: n1 wrong sizes ! "
 if (n2dim_in/=n2dim_out) stop "Error in cprj_copy: n2 wrong sizes ! "
 if (ncpgr_in/=ncpgr_out) stop "Error in cprj_copy: ncpgr wrong sizes ! "

 do jj=1,n2dim_in
  do ii=1,n1dim_in
   nlmn=cprj_in(ii,jj)%nlmn
   cprj_out(ii,jj)%nlmn =nlmn
   do kk=1,nlmn
    cprj_out(ii,jj)%cp(1:2,kk)=cprj_in(ii,jj)%cp(1:2,kk)
   end do
  end do
 end do

 if (ncpgr_in>0) then
  do jj=1,n2dim_in
   do ii=1,n1dim_in
    nlmn=cprj_in(ii,jj)%nlmn
    do kk=1,nlmn
     cprj_out(ii,jj)%dcp(1:2,1:ncpgr_in,kk)=cprj_in(ii,jj)%dcp(1:2,1:ncpgr_in,kk)
    end do
   end do
  end do
 end if

end subroutine cprj_copy
!!***


!!****f* ABINIT/cprj_axpby
!! NAME
!! cprj_axpby
!!
!! FUNCTION
!! Apply AXPBY (blas-like) operation with 2 cprj datastructures:
!!  cprjy(:,:) <- alpha.cprjx(:,:)+beta.cprjy(:,:)
!!
!! COPYRIGHT
!! Copyright (C) 2008-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  alpha,beta= alpha,beta REAL factors
!!  cprjx(:,:) <type(cprj_type)>= input cprjx datastructure
!!
!! SIDE EFFECTS
!!  cprjy(:,:) <type(cprj_type)>= input/output cprjy datastructure
!!
!! PARENTS
!!
!! CHILDREN
!!      xallgather_mpi,xme_init
!!
!! SOURCE

 subroutine cprj_axpby(alpha,beta,cprjx,cprjy)

 use defs_basis
 use defs_datatypes

 implicit none
!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: alpha,beta
!arrays
 type(cprj_type),intent(in) :: cprjx(:,:)
 type(cprj_type),intent(inout) :: cprjy(:,:)
!Local variables-------------------------------
 integer :: ii,jj,kk,n1dimx,n1dimy,n2dimx,n2dimy,ncpgrx,ncpgry,nlmn

! *************************************************************************

 n1dimx=size(cprjx,dim=1);n1dimy=size(cprjy,dim=1)
 n2dimx=size(cprjx,dim=2);n2dimy=size(cprjy,dim=2)
 ncpgrx=cprjx(1,1)%ncpgr;ncpgry=cprjy(1,1)%ncpgr
 if (n1dimx/=n1dimy) stop "Error in cprj_axpby: n1 wrong sizes ! "
 if (n2dimx/=n2dimy) stop "Error in cprj_axpby: n2 wrong sizes ! "
 if (ncpgrx/=ncpgry) stop "Error in cprj_axpby: ncpgr wrong sizes ! "

 do jj=1,n2dimx
  do ii=1,n1dimx
   nlmn=cprjx(ii,jj)%nlmn
   cprjy(ii,jj)%nlmn =nlmn
   do kk=1,nlmn
    cprjy(ii,jj)%cp(1:2,kk)=alpha*cprjx(ii,jj)%cp(1:2,kk) &
&    +beta *cprjy(ii,jj)%cp(1:2,kk)
   end do
  end do
 end do

 if (ncpgrx>0) then
  do jj=1,n2dimx
   do ii=1,n1dimx
    nlmn=cprjx(ii,jj)%nlmn
    do kk=1,nlmn
     cprjy(ii,jj)%dcp(1:2,1:ncpgrx,kk)=alpha*cprjx(ii,jj)%dcp(1:2,1:ncpgrx,kk) &
&     +beta *cprjx(ii,jj)%dcp(1:2,1:ncpgrx,kk)
    end do
   end do
  end do
 end if

end subroutine cprj_axpby
!!***


!!****f* ABINIT/cprj_diskinit
!! NAME
!! cprj_diskinit
!!
!! FUNCTION
!! Initialize reading of a cprj temporary file
!! Nothing is done if mkmem=0
!!
!! COPYRIGHT
!! Copyright (C) 2008-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dimcp=first dimension of cprj arrays (1 or natom)
!!  iorder=0 if cprj ordering do not change during reading
!!         1 if cprj ordering change during reading (type-sorted->unsorted)
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  natom=number of atoms in cell
!!  nlmn(dimcp)=array of dimensions of cprj datastructure that will contain the read data
!!  nspinor=number of spinorial components of the wavefunctions
!!  uncp=unit number for cprj data (if used)
!!
!! PARENTS
!!      dyfnl3,optics_paw,partial_dos_fractions_paw,pawmkrhoij,pawmkrhoij3
!!      suscep_stat,vtorho3
!!
!! CHILDREN
!!      xallgather_mpi,xme_init
!!
!! SOURCE

 subroutine cprj_diskinit(atindx1,dimcp,iorder,mkmem,natom,nlmn,nspinor,uncp)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimcp,iorder,mkmem,natom,nspinor,uncp
!arrays
 integer,intent(in) :: atindx1(natom),nlmn(dimcp)
!Local variables-------------------------------
 integer :: dimcp0,iatm,iatom,nspinor0
 character(len=500) :: message
 integer,allocatable :: dimlmn(:)

! *************************************************************************

 if (mkmem==0) then

  rewind uncp;read(uncp) dimcp0,nspinor0
  if (dimcp/=dimcp0.or.nspinor/=nspinor0) then
   write(message,'(a,a,a,a)')ch10,&
&   ' cprj_diskinit : BUG -',ch10,&
&   '  _PAW file was not created with the right options !'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

  allocate(dimlmn(dimcp))
  read(uncp) dimlmn(1:dimcp)
  do iatm=1,dimcp
   if (iorder==0) then
    iatom=iatm
   else
    iatom=min(atindx1(iatm),dimcp)
   end if
   if (dimlmn(iatm)/=nlmn(iatom)) then
    write(message,'(a,a,a,a)')ch10,&
&    ' cprj_diskinit : BUG -',ch10,&
&    '  _PAW file was not created with the right options !'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
  end do
  deallocate(dimlmn)

 end if

end subroutine cprj_diskinit
!!***


!!****f* ABINIT/cprj_get
!! NAME
!! cprj_get
!!
!! FUNCTION
!! Read the cprj for a given k-point from memory or from a temporary file
!!
!! COPYRIGHT
!! Copyright (C) 2008-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cprj(dimcp,nspinor*mband*mkmem*nsppol)=input cprj (used if mkmem/=0)
!!  dimcp=first dimension of cprj_k,cprj arrays (1 or natom)
!!  iband1=index of first band
!!  ibg=shift if cprj array to locate current k-point
!!  ikpt=index of current k-point
!!  iorder=0 if cprj ordering do not change during reading
!!         1 if cprj ordering change during reading (type-sorted->unsorted)
!!  isppol=index of current spin component
!!  mband=maximum number of bands
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands to import (usually 1 or nband_k)
!!  nband_k=total number of bands for this k-point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  uncp=unit number for cprj data (used if mkmem=0)
!!
!! OUTPUT
!!  cprj_k(dimcp,nspinor*nband_k) <type(cprj_type)>= output cprj datastructure
!!
!! PARENTS
!!      dyfnl3,optics_paw,partial_dos_fractions_paw,pawmkrhoij,pawmkrhoij3
!!      suscep_stat,vtorho3
!!
!! CHILDREN
!!      xallgather_mpi,xme_init
!!
!! SOURCE

 subroutine cprj_get(atindx1,cprj_k,cprj,dimcp,iband1,ibg,ikpt,iorder,isppol,&
&                    mband,mkmem,mpi_enreg,natom,nband,nband_k,nspinor,nsppol,uncp)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dimcp,iband1,ibg,ikpt,iorder,isppol,mband,mkmem,natom,nband,nband_k,nspinor,nsppol,uncp
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom)
 type(cprj_type),intent(in) :: cprj(dimcp,nspinor*mband*mkmem*nsppol)
 type(cprj_type),intent(out) :: cprj_k(dimcp,nspinor*nband)
!Local variables-------------------------------
 integer :: iatm,iatom,ib,ibsp,isp,ispinor,jband,me,nband0,ncpgr
 character(len=500) :: message

! *************************************************************************

 ncpgr=cprj_k(1,1)%ncpgr

 if ((mpi_enreg%paral_compil_kpt==1) .and. &
& (mpi_enreg%paral_compil_fft==1)) then
  me=mpi_enreg%me_kpt
 else
  call xme_init(mpi_enreg,me)
 end if

 if (mkmem==0) then

  if (iband1==1) then
   read(uncp) nband0
   if (nband_k/=nband0) then
    write(message,'(a,a,a,a)')ch10,&
&    ' cprj_get : BUG -',ch10,&
&    '  _PAW file was not created with the right options !'
    call wrtout(6,message,'PERS')
    call leave_new('PERS')
   end if
  end if

  isp=0;jband=iband1-1
  do ib=1,nband
   jband=jband+1
   if (mpi_enreg%paral_compil_kpt==1) then
    if (abs(mpi_enreg%proc_distrb(ikpt,jband,isppol)-me)/=0) then
     isp=isp+nspinor
     cycle
    end if
   end if
   do ispinor=1,nspinor
    isp=isp+1
    if (iorder==0) then
     do iatom=1,dimcp
      if (ncpgr==0) then
       read(uncp) cprj_k(iatom,isp)%cp(:,:)
      else
       read(uncp) cprj_k(iatom,isp)%cp(:,:),cprj_k(iatom,isp)%dcp(:,:,:)
      end if
     end do
    else
     do iatm=1,dimcp
      iatom=min(atindx1(iatm),dimcp)
      if (ncpgr==0) then
       read(uncp) cprj_k(iatom,isp)%cp(:,:)
      else
       read(uncp) cprj_k(iatom,isp)%cp(:,:),cprj_k(iatom,isp)%dcp(:,:,:)
      end if
     end do
    end if
   end do
  end do

 else

  isp=0;ibsp=ibg+nspinor*(iband1-1);jband=iband1-1
  do ib=1,nband
   jband=jband+1
   if (mpi_enreg%paral_compil_kpt==1) then
    if (abs(mpi_enreg%proc_distrb(ikpt,jband,isppol)-me)/=0) then
     isp=isp+nspinor;ibsp=ibsp+nspinor
     cycle
    end if
   end if
   do ispinor=1,nspinor
    isp=isp+1;ibsp=ibsp+1
    if (iorder==0) then
     do iatom=1,dimcp
      cprj_k(iatom,isp)%cp(:,:)=cprj(iatom,ibsp)%cp(:,:)
      if (ncpgr>0) cprj_k(iatom,isp)%dcp(:,:,:)=cprj(iatom,ibsp)%dcp(:,:,:)
     end do
    else
     do iatm=1,dimcp
      iatom=min(atindx1(iatm),dimcp)
      cprj_k(iatom,isp)%cp(:,:)=cprj(iatm,ibsp)%cp(:,:)
      if (ncpgr>0) cprj_k(iatom,isp)%dcp(:,:,:)=cprj(iatm,ibsp)%dcp(:,:,:)
     end do
    end if
   end do
  end do

 end if

end subroutine cprj_get
!!***


!!****f* ABINIT/cprj_put
!! NAME
!! cprj_put
!!
!! FUNCTION
!! Write the cprj for a given set of (n,k) into memory or into a temporary file
!!
!! COPYRIGHT
!! Copyright (C) 2007-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  are_gathered=TRUE if cprj_k arrays have already been gathered between procs,
!!               (band-fft parallelism only)
!!               Typically, TRUE after a call to prep_nonlop routine...
!!  atindx(natom)=index table for atoms
!!  cprj_k(dimcp,nspinor*nband) <type(cprj_type)>= input cprj datastructure
!!  dimcp=first dimension of cprj_k,cprjnk arrays (1 or natom)
!!  iband1=index of first band
!!  ibg=shift if cprjnk array to locate current k-point
!!  ikpt=index of current k-point
!!  iorder=0 if cprj ordering do not change during reading
!!         1 if cprj ordering change during reading (unsorted->type-sorted)
!!  isppol=index of current spin component
!!  mband=maximum number of bands
!!  mkmem=number of k points which can fit in memory; set to 0 if use disk
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell
!!  nband=number of bands to export (usually 1, nband_k or nblockbd)
!!  nband_k=total number of bands for this k-point
!!  nlmn(dimcp)=array of dimensions of cprj_k,cprjnk datastructures
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  spaceComm_band=communicator used for bands in case of band-fft parallelism
!!  uncp=unit number for cprj data (used if mkmem=0)
!!
!! SIDE EFFECTS
!!  cprj(dimcp,nspinor*mband*mkmem*nsppol)=output cprj (used if mkmem/=0)
!!
!! PARENTS
!!      ctocprj,vtowfk,vtowfk3
!!
!! CHILDREN
!!      xallgather_mpi,xme_init
!!
!! SOURCE

 subroutine cprj_put(are_gathered,atindx,cprj_k,cprj,dimcp,iband1,ibg,ikpt,iorder,isppol,&
&           mband,mkmem,mpi_enreg,natom,nband,nband_k,nlmn,nspinor,nsppol,spaceComm_band,uncp)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none
!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband1,ibg,ikpt,iorder,isppol,dimcp,mband,mkmem
 integer,intent(in) :: natom,nband,nband_k,nspinor,nsppol,uncp,spaceComm_band
 logical,intent(in) :: are_gathered
 type(MPI_type), intent(inout) :: mpi_enreg
!arrays
 integer :: atindx(natom),nlmn(dimcp)
 type(cprj_type),intent(out) :: cprj(dimcp,nspinor*mband*mkmem*nsppol)
 type(cprj_type),intent(in) :: cprj_k(dimcp,nspinor*nband)
!Local variables-------------------------------
 integer :: iatm,iatom,iband,ibsp,icpgr,ierr,ii,ilmn,isp,ispinor,jband,jj,lmndim,me,ncpgr
 real(dp),allocatable :: buffer1(:),buffer2(:)

! *************************************************************************

 ncpgr=cprj_k(1,1)%ncpgr

 if ((mpi_enreg%paral_compil_kpt==1) .and. &
& (mpi_enreg%paral_compil_fft==1)) then
  me=mpi_enreg%me_kpt
 else
  call xme_init(mpi_enreg,me)
 end if

 if (mpi_enreg%mode_para/='b'.or.are_gathered.or.nband==1) then

  if (mkmem==0) then

   if (iband1==1) write(uncp) nband_k

   isp=0;jband=iband1-1
   do iband=1,nband
    jband=jband+1
    if (mpi_enreg%paral_compil_kpt==1) then
     if (abs(mpi_enreg%proc_distrb(ikpt,jband,isppol)-me)/=0) then
      isp=isp+nspinor
      cycle
     end if
    end if
    do ispinor=1,nspinor
     isp=isp+1
     if (iorder==0) then
      do iatom=1,dimcp
       if (ncpgr==0) then
        write(uncp) cprj_k(iatom,isp)%cp(:,:)
       else
        write(uncp) cprj_k(iatom,isp)%cp(:,:),cprj_k(iatom,isp)%dcp(:,:,:)
       end if
      end do
     else
      do iatom=1,dimcp
       iatm=min(atindx(iatom),dimcp)
       if (ncpgr==0) then
        write(uncp) cprj_k(iatm,isp)%cp(:,:)
       else
        write(uncp) cprj_k(iatm,isp)%cp(:,:),cprj_k(iatm,isp)%dcp(:,:,:)
       end if
      end do
     end if
    end do
   end do

  else

   isp=0;ibsp=ibg+nspinor*(iband1-1);jband=iband1-1
   do iband=1,nband
    jband=jband+1
    if (mpi_enreg%paral_compil_kpt==1)then
     if (abs(mpi_enreg%proc_distrb(ikpt,jband,isppol)-me)/=0) then
      isp=isp+nspinor;ibsp=ibsp+nspinor
      cycle
     end if
    end if
    do ispinor=1,nspinor
     isp=isp+1;ibsp=ibsp+1
     if (iorder==0) then
      do iatom=1,dimcp
       cprj(iatom,ibsp)%cp(:,:)=cprj_k(iatom,isp)%cp(:,:)
       if (ncpgr>0) cprj(iatom,ibsp)%dcp(:,:,:)=cprj_k(iatom,isp)%dcp(:,:,:)
      end do
     else
      do iatom=1,dimcp
       iatm=min(atindx(iatom),dimcp)
       cprj(iatom,ibsp)%cp(:,:)=cprj_k(iatm,isp)%cp(:,:)
       if (ncpgr>0) cprj(iatom,ibsp)%dcp(:,:,:)=cprj_k(iatm,isp)%dcp(:,:,:)
      end do
     end if
    end do
   end do

  end if

 else ! mode_para==b and nband>1

  lmndim=2*sum(nlmn(1:dimcp))
  allocate(buffer1((1+ncpgr)*lmndim*nspinor))
  allocate(buffer2((1+ncpgr)*lmndim*nspinor*mpi_enreg%nproc_band))
  isp=0;ibsp=ibg+nspinor*(iband1-1)
  do iband=1,nband  ! must be nblockbd for band-fft parallelism
   jj=1
   do ispinor=1,nspinor
    isp=isp+1
    do iatom=1,dimcp
     if (iorder==0) then
      iatm=iatom
     else
      iatm=min(atindx(iatom),dimcp)
     end if
     do ilmn=1,nlmn(iatm)
      buffer1(jj:jj+1)=cprj_k(iatm,isp)%cp(1:2,ilmn)
      jj=jj+2
     end do
     if (ncpgr>0) then
      do ilmn=1,nlmn(iatm)
       do icpgr=1,ncpgr
        buffer1(jj:jj+1)=cprj_k(iatm,isp)%dcp(1:2,icpgr,ilmn)
        jj=jj+2
       end do
      end do
     end if
    end do
   end do
   call xallgather_mpi(buffer1,lmndim,buffer2,spaceComm_band,ierr)
   jj=1
   do ii=1,mpi_enreg%nproc_band
    do ispinor=1,nspinor
     ibsp=ibsp+1
     do iatom=1,dimcp
      do ilmn=1,nlmn(iatom)
       cprj(iatom,ibsp)%cp(1:2,ilmn)=buffer2(jj:jj+1)
       jj=jj+2
      end do
     end do
     if (ncpgr>0) then
      do ilmn=1,nlmn(iatom)
       do icpgr=1,ncpgr
        cprj(iatom,ibsp)%dcp(1:2,icpgr,ilmn)=buffer2(jj:jj+1)
        jj=jj+2
       end do
      end do
     end if
    end do
   end do
  end do
  deallocate(buffer1,buffer2)

 end if ! mode_para=b, nband

end subroutine cprj_put
!!***
