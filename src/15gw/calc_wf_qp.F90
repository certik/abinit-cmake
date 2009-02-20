!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_wf_qp
!! NAME
!! calc_wf_qp
!!
!! FUNCTION
!!  Calculate QP amplitudes in real or reciprocal space starting from the 
!!  KS wavefunctions and the corresponding expansion coefficients.
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2008 ABINIT group (FBruneval, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  b1gw, b2gw = Min and max band index over k-point and spin for GW corrections.
!!  nkibz=number of k-points.
!!  nsize= number of points in real space or number of G vectors.
!!  nsppol=number of spin.
!!  nbnds=number of bands in the present GW calculation.
!!  my_minb, my_maxb = Indeces of the bands treated by this processor.
!!  m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)= expansion of the QP amplitudes in terms of KS wavefunctions.
!!
!! OUTPUT
!!  wf(nsize,my_minb:my_maxb,nkibz,nsppol)= Updated QP amplitudes for this processor.
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      cgemm
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_wf_qp(MPI_enreg,nkibz,nbnds,nsize,nsppol,nspinor,&
& m_lda_to_qp,my_minb,my_maxb,b1gw,b2gw,wf)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : czero_gw, cone_gw
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nkibz,nsize,nsppol,nspinor,my_minb,my_maxb,b1gw,b2gw
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 complex(gwpc),intent(inout) :: wf(nsize*nspinor,my_minb:my_maxb,nkibz,nsppol)

!Local variables-------------------------------
#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'CGEMM' :: cgemm
#endif
!scalars
 integer :: ib,ik,is,ierr,lowerb,upperb,rangeb,ispinor,spad
 integer :: spaceComm,sizegw,rank,shift,indx_kibz,ilmn,nlmn
 character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: umat_k(:,:)
 complex(gwpc),allocatable :: wf_ks(:,:),wf_qp(:,:)
!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_wf_qp: enter ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 call xcomm_init(MPI_enreg,spaceComm)
 call xme_init  (MPI_enreg,rank)
 !
 ! === Determine the range of bands this processor has to treat ===
 lowerb=0  ! There is no overlap between [b1gw,b2gw] and [my_minb,my_maxb]
 if (b1gw<=my_maxb) lowerb=MAX(b1gw,my_minb)
 upperb=0  ! There is no overlap between [b1gw,b2gw] and [my_minb,my_maxb]
 if (b2gw>=my_minb) upperb=MIN(b2gw,my_maxb)
 rangeb=0
 if (lowerb/=0.and.upperb/=0) rangeb=upperb-lowerb+1
 sizegw=b2gw-b1gw+1

 if (rangeb/=0) then
  write(msg,'(2a,i3,3a,2i3,a,3i3,2a,3i3)')ch10,&
&  ' proc ',rank,' will update its wavefunctions ',ch10,&
&  ' my_bands indeces: ',my_minb,my_maxb,' gwrange: ',b1gw,b2gw,sizegw,ch10,&
&  ' lowerb, upperb, rangeb: ',lowerb,upperb,rangeb
  call wrtout(std_out,msg,'PERS') 
 end if

 allocate(umat_k(lowerb:upperb,b1gw:b2gw))
 allocate(wf_qp(nsize*nspinor,b1gw:b2gw))  
 allocate(wf_ks(nsize,lowerb:upperb))
 wf_qp(:,:)=czero_gw ; wf_ks(:,:)=czero_gw 
 !
 ! === Calculate : $\Psi^{QP}_{r,b} = \sum_n \Psi^{KS}_{r,n} M_{n,b}$ ===
 do is=1,nsppol
  do ik=1,nkibz

   umat_k(:,:)=m_lda_to_qp(lowerb:upperb,b1gw:b2gw,ik,is)
   wf_qp(:,:)=czero_gw

   if (rangeb/=0) then
    do ispinor=1,nspinor
     spad=nsize*(ispinor-1)
     wf_ks(:,lowerb:upperb)=wf(spad+1:spad+nsize,lowerb:upperb,ik,is)
#if defined HAVE_GW_DPC
     call ZGEMM('N','N',nsize,sizegw,rangeb,cone_gw,wf_ks(:,lowerb:upperb),nsize,&
&     umat_k,rangeb,czero_gw,wf_qp(spad+1:spad+nsize,b1gw:b2gw),nsize)
#else
     call CGEMM('N','N',nsize,sizegw,rangeb,cone_gw,wf_ks(:,lowerb:upperb),nsize,&
&     umat_k,rangeb,czero_gw,wf_qp(spad+1:spad+nsize,b1gw:b2gw),nsize)
#endif
    end do
   end if
   !
   ! === Update the input wave functions ===
   if (MPI_enreg%gwpara==2) then 
    ! == Bands are spreaded among processors ==
    ! * Sum up all the partial QP amplitudes.
    ! * Keep the band in memory only if you are the right processor.
    call xsum_mpi(wf_qp(:,b1gw:b2gw),spaceComm,ierr)
    do ib=b1gw,b2gw
     if (rank==MPI_enreg%proc_distrb(ik,ib,is)) wf(:,ib,ik,is)=wf_qp(:,ib)
    end do
   else
    ! == Each node has the full set ==
    wf(:,b1gw:b2gw,ik,is)=wf_qp(:,b1gw:b2gw)
   end if

  end do !ik
 end do !is

 deallocate(umat_k)
 deallocate(wf_ks,wf_qp)

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_wf_qp ended ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine calc_wf_qp
!!***


!!****f* ABINIT/calc_wf_qp_Wfval
!! NAME
!! calc_wf_qp_Wfval
!!
!! FUNCTION
!!  Calculate QP amplitudes in real or reciprocal space starting from the 
!!  KS wavefunctions and the corresponding expansion coefficients,
!!  in the case of two separated sets of wavefunctions: Wf and Wfval
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2008 ABINIT group (FBruneval, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  b1gw, b2gw = Min and max band index over k-point and spin for GW corrections.
!!  nkibz=number of k-points.
!!  nsize= number of points in real space or number of G vectors.
!!  nsppol=number of spin.
!!  nbnds=number of bands in the present GW calculation.
!!  my_minb, my_maxb = Indeces of the bands treated by this processor.
!!  m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)= expansion of the QP amplitudes in terms of KS wavefunctions.
!!
!! OUTPUT
!!  wf(nsize,my_minb:my_maxb,nkibz,nsppol)= Updated QP amplitudes for this processor.
!!
!! PARENTS
!!      screening,sigma
!!
!! CHILDREN
!!      cgemm
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_wf_qp_Wfval(MPI_enreg,nkibz,nbnds,nsize,nsppol,nspinor,&
& m_lda_to_qp,my_minb,my_maxb,b1gw,b2gw,wf,nbvw,wfval)

 use defs_basis
 use defs_datatypes
 use m_gwdefs, only : czero_gw, cone_gw
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_lib01hidempi
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nbnds,nkibz,nsize,nsppol,nspinor,my_minb,my_maxb,b1gw,b2gw,nbvw
 type(MPI_type),intent(inout) :: MPI_enreg
!arrays
 complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 complex(gwpc),intent(inout) :: wf(nsize,my_minb:my_maxb,nkibz,nsppol)
 complex(gwpc),intent(inout) :: wfval(nsize,nbvw,nkibz,nsppol)

!Local variables-------------------------------
#ifdef __VMS
!DEC$ ATTRIBUTES ALIAS:'CGEMM' :: cgemm
#endif
!scalars
 integer :: ib,ik,is,ierr,lowerb,upperb,rangeb,b1gw_c
 integer :: spaceComm,sizegw,rank,shift,indx_kibz,ilmn,nlmn
 character(len=500) :: msg
!arrays
 complex(gwpc),allocatable :: umat_k(:,:)
 complex(gwpc),allocatable :: wf_qp(:,:)
 complex(gwpc),allocatable :: wf_qp_valval(:,:),wf_qp_val(:,:),wf_qp_valcon(:,:)
!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_wf_qp_Wfval: enter ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 call xcomm_init(MPI_enreg,spaceComm)
 call xme_init  (MPI_enreg,rank)
 !
 ! === Determine the range of bands this processor has to treat ===
 b1gw_c=MAX(b1gw,nbvw+1)! to avoid double counting of the valence bands
 lowerb=0  ! There is no overlap between [b1gw,b2gw] and [my_minb,my_maxb]
 if (b1gw_c<=my_maxb) then
  lowerb=MAX(b1gw_c,my_minb)
 end if
 upperb=0  ! There is no overlap between [b1gw,b2gw] and [my_minb,my_maxb]
 if (b2gw>=my_minb) upperb=MIN(b2gw,my_maxb)
 rangeb=0
 if (lowerb/=0.and.upperb/=0) rangeb=upperb-lowerb+1
 sizegw=b2gw-b1gw_c+1
 
 if (rangeb>0) then
  write(msg,'(2a,i3,3a,2i3,a,3i3,2a,3i3)')ch10,&
&  ' proc ',rank,' will update its wavefunctions ',ch10,&
&  ' my_bands indeces: ',my_minb,my_maxb,' gwrange: ',b1gw_c,b2gw,sizegw,ch10,&
&  ' lowerb, upperb, rangeb: ',lowerb,upperb,rangeb
  call wrtout(std_out,msg,'PERS') 
 end if
 
 allocate(wf_qp_valval(nsize,nbvw))
 allocate(wf_qp_val(nsize,nbvw))  

 if (sizegw>0) then
  allocate(wf_qp_valcon(nsize,b1gw_c:b2gw))  
  allocate(wf_qp(nsize,b1gw_c:b2gw))  
 end if
 !
 ! === Calculate : $\Psi^{QP}_{r,b} = \sum_n \Psi^{KS}_{r,n} M_{n,b}$ ===
 do is=1,nsppol
  do ik=1,nkibz
   !
   ! I) Treat the valence bands
   !
   wf_qp_valval(:,:)=czero_gw ; wf_qp_val(:,:)=czero_gw 
   allocate(umat_k(nbvw,nbvw))
   umat_k(:,:)=m_lda_to_qp(1:nbvw,1:nbvw,ik,is)
 
#if defined HAVE_GW_DPC
   call ZGEMM('N','N',nsize,nbvw,nbvw,cone_gw,wfval(:,1:nbvw,ik,is),nsize,&
&   umat_k,nbvw,czero_gw,wf_qp_valval(:,1:nbvw),nsize)
#else
   call CGEMM('N','N',nsize,nbvw,nbvw,cone_gw,wfval(:,1:nbvw,ik,is),nsize,&
&   umat_k,nbvw,czero_gw,wf_qp_valval(:,1:nbvw),nsize)
#endif
   deallocate(umat_k)
  
   if (rangeb>0) then
    allocate(umat_k(lowerb:upperb,1:nbvw))
    umat_k(:,:)=m_lda_to_qp(lowerb:upperb,1:nbvw,ik,is)
#if defined HAVE_GW_DPC
    call ZGEMM('N','N',nsize,nbvw,rangeb,cone_gw,wf(:,lowerb:upperb,ik,is),nsize,&
&    umat_k,rangeb,czero_gw,wf_qp_val(:,1:nbvw),nsize)
#else
    call CGEMM('N','N',nsize,nbvw,rangeb,cone_gw,wf(:,lowerb:upperb,ik,is),nsize,&
&    umat_k,rangeb,czero_gw,wf_qp_val(:,1:nbvw),nsize)
#endif
    deallocate(umat_k)
   end if
   if (MPI_enreg%gwpara==2) then 
    ! Bands are spreaded among processors:
    ! * Sum up all the partial QP amplitudes.
    ! * Keep the band in memory only if you are the right processor.
    call xsum_mpi(wf_qp_val(:,1:nbvw),spaceComm,ierr)
    wf_qp_valval(:,1:nbvw)=wf_qp_valval(:,1:nbvw) + wf_qp_val(:,1:nbvw)
   else
    ! Each node has the full set
    wf_qp_valval(:,1:nbvw)=wf_qp_valval(:,1:nbvw) + wf_qp_val(:,1:nbvw)
   end if
   !
   ! II) Treat the NON-valence bands
   !
   if (sizegw>0) then
    wf_qp_valcon(:,:)=czero_gw
    wf_qp(:,:)=czero_gw
    allocate(umat_k(1:nbvw,b1gw_c:b2gw))
    umat_k(:,:)=m_lda_to_qp(1:nbvw,b1gw_c:b2gw,ik,is)
 
#if defined HAVE_GW_DPC
    call ZGEMM('N','N',nsize,sizegw,nbvw,cone_gw,wfval(:,1:nbvw,ik,is),nsize,&
&    umat_k,nbvw,czero_gw,wf_qp_valcon(:,b1gw_c:b2gw),nsize)
#else
    call CGEMM('N','N',nsize,sizegw,nbvw,cone_gw,wfval(:,1:nbvw,ik,is),nsize,&
&    umat_k,nbvw,czero_gw,wf_qp_valcon(:,b1gw_c:b2gw),nsize)
#endif
    deallocate(umat_k)
 
    if (rangeb>0) then
     allocate(umat_k(lowerb:upperb,b1gw_c:b2gw))
     umat_k(:,:)=m_lda_to_qp(lowerb:upperb,b1gw_c:b2gw,ik,is)
#if defined HAVE_GW_DPC
     call ZGEMM('N','N',nsize,sizegw,rangeb,cone_gw,wf(:,lowerb:upperb,ik,is),nsize,&
&     umat_k,rangeb,czero_gw,wf_qp(:,b1gw_c:b2gw),nsize)
#else
     call CGEMM('N','N',nsize,sizegw,rangeb,cone_gw,wf(:,lowerb:upperb,ik,is),nsize,&
&     umat_k,rangeb,czero_gw,wf_qp(:,b1gw_c:b2gw),nsize)
#endif
     deallocate(umat_k)
    end if
    !
    ! === Update the input wave functions ===
    if (MPI_enreg%gwpara==2) then 
     ! Bands are spreaded among processors:
     ! * Sum up all the partial QP amplitudes.
     ! * Keep the band in memory only if you are the right processor.
     call xsum_mpi(wf_qp(:,b1gw_c:b2gw),spaceComm,ierr)
     do ib=b1gw_c,b2gw
      if (rank==MPI_enreg%proc_distrb(ik,ib,is)) wf(:,ib,ik,is)=wf_qp(:,ib)+wf_qp_valcon(:,ib)
     end do
    else
     ! Each node has the full set
     wf(:,b1gw_c:b2gw,ik,is)=wf_qp_valcon(:,b1gw_c:b2gw)+wf_qp(:,b1gw_c:b2gw)
    end if
 
   endif !sizegw>0
 
   wfval(:,:,ik,is)=wf_qp_valval(:,:)
   wf   (:,my_minb:b1gw_c-1,ik,is)=wf_qp_valval(:,my_minb:b1gw_c-1)
 
  end do !ik
 end do !is

 call leave_test(mpi_enreg)

 if (allocated(wf_qp       )) deallocate(wf_qp)
 if (allocated(wf_qp_valval)) deallocate(wf_qp_valval)
 if (allocated(wf_qp_val   )) deallocate(wf_qp_val)
 if (allocated(wf_qp_valcon)) deallocate(wf_qp_valcon)

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_wf_qp_Wfval: ended ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine calc_wf_qp_Wfval
!!***


!!****f* ABINIT/update_cprj
!! NAME
!! update_cprj
!!
!! FUNCTION
!!  Update matrix elements of the PAW projectors in case of self-consistent GW.
!!
!! COPYRIGHT
!!  Copyright (C) 2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  dimlmn(natom)=number of (l,m,n) components for each atom (only for PAW)
!!  nkibz=number of k-points
!!  nsppol=number of spin
!!  nbnds=number of bands in the present GW calculation
!!  m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)= expansion of the QP amplitudes in terms of KS wavefunctions
!!  natom=number of atomd in unit cell
!!
!! OUTPUT
!!  Cprj_ibz(natom,nspinor*nkibz*nbnds*nsppol) <type(cprj_type)>=projected wave functions 
!!   <Proj_i|Cnk> with all NL projectors. On exit, it contains the projections onto the 
!!   QP amplitudes.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine update_cprj(natom,nkibz,nbnds,nsppol,nspinor,m_lda_to_qp,dimlmn,Cprj_ibz)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13nonlocal
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nbnds,nkibz,nsppol,nspinor
!arrays
 integer,intent(in) :: dimlmn(natom)
 complex(dpc),intent(in) :: m_lda_to_qp(nbnds,nbnds,nkibz,nsppol)
 type(Cprj_type),intent(inout) :: Cprj_ibz(natom,nspinor*nbnds*nkibz*nsppol)

!Local variables-------------------------------
!scalars
 integer :: iat,ib,ik,is,shift,indx_kibz,ilmn,nlmn,ispinor,ibsp,spad,ibdx
 character(len=500) :: msg
!arrays
 real(dp),target,allocatable :: real_lda2qp(:,:,:,:,:)
 real(dp),allocatable :: re_p(:),im_p(:),vect(:,:)
 real(dp),pointer :: umat(:,:,:)
 type(Cprj_type),allocatable :: Cprj_ks(:,:)

!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' update_cprj: enter ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

 allocate(Cprj_ks(natom,nspinor*nbnds)) ; call cprj_alloc(Cprj_ks,0,dimlmn)

 allocate(real_lda2qp(2,nbnds,nbnds,nkibz,nsppol))
 real_lda2qp(1,:,:,:,:)=REAL (m_lda_to_qp)
 real_lda2qp(2,:,:,:,:)=AIMAG(m_lda_to_qp)
 allocate(re_p(nbnds),im_p(nbnds),vect(2,nbnds))
 !
 ! === Calculate : $\Psi^{QP}_{r,b} = \sum_n \Psi^{KS}_{r,n} M_{n,b}$ ===
 do is=1,nsppol
  do ik=1,nkibz

   shift=nspinor*nbnds*nkibz*(is-1)
   indx_kibz=nspinor*nbnds*(ik-1)+shift
   ibsp=0
   do ib=1,nbnds
    do ispinor=1,nspinor
     ibsp=ibsp+1
     do iat=1,natom
      Cprj_ks(iat,ibsp)%cp(:,:)=Cprj_ibz(iat,indx_kibz+ibsp)%cp(:,:)
     end do
    end do
   end do
   
   umat => real_lda2qp(:,:,:,ik,is)
   do iat=1,natom
    nlmn=dimlmn(iat)
    do ilmn=1,nlmn

     do ispinor=1,nspinor
      ! Retrieve projections for this spinorial component
      spad=(ispinor-1)
      ibdx=0
      do ib=1,nbnds*nspinor,nspinor 
       ibdx=ibdx+1
       vect(1,ibdx)=Cprj_ks(iat,ib+spad)%cp(1,ilmn)
       vect(2,ibdx)=Cprj_ks(iat,ib+spad)%cp(2,ilmn)
      end do

      re_p(:)= &
&       MATMUL(umat(1,:,:),vect(1,:)) &
&      -MATMUL(umat(2,:,:),vect(2,:))
      im_p(:)= &
&       MATMUL(umat(1,:,:),vect(2,:)) &
&      +MATMUL(umat(2,:,:),vect(1,:))
      
      ! === Save values ===
      ibdx=0
      do ib=1,nbnds*nspinor,nspinor 
       ibdx=ibdx+1
       Cprj_ibz(iat,indx_kibz+spad+ib)%cp(1,ilmn)=re_p(ibdx)
       Cprj_ibz(iat,indx_kibz+spad+ib)%cp(2,ilmn)=im_p(ibdx)
      end do
     end do !ispinor

    end do !ilmn
   end do !iat

  end do !ik
 end do !is

 deallocate(real_lda2qp,re_p,im_p,vect)
 call cprj_free(Cprj_ks) ; deallocate(Cprj_ks)

#if defined DEBUG_MODE
 write(msg,'(a)')' update_cprj ended ' 
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine update_cprj
!!***
