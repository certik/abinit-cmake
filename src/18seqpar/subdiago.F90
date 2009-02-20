!{\src2tex{textfont=tt}}
!!****f* ABINIT/subdiago
!! NAME
!! subdiago        
!!
!! FUNCTION
!! This routine diagonalizes the Hamiltionian is the eigenfunctions subspace
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  filstat=name of the status file
!!  icg=shift to be applied on the location of data in the array cg
!!  igsc=shift to be applied on the location of data in the array gsc
!!  ikpt=number of the k-point
!!  inonsc=index of non self-consistent loop
!!  istwf_k=input parameter that describes the storage of wfs
!!  mcg=second dimension of the cg array
!!  mgsc=second dimension of the gsc array
!!  mpi_enreg=informations about MPI parallelization
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  nspinor=number of spinorial components of the wavefunctions
!!  prtvol=control print volume and debugging output
!!  subham(nband_k*(nband_k+1))=Hamiltonian expressed in sthe WFs subspace
!!  subovl(nband_k*(nband_k+1)*use_subovl)=overlap matrix expressed in sthe WFs subspace
!!  use_subovl=1 if the overlap matrix is not identity in WFs subspace
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  eig_k(nband_k)=array for holding eigenvalues (hartree)
!!  evec(2*nband_k,nband_k)=array for holding eigenvectors
!!
!! SIDE EFFECTS
!!  cg(2,mcg)=wavefunctions
!!  gsc(2,mgsc)=<g|S|c> matrix elements (S=overlap)
!!
!! PARENTS
!!      vtowfk
!!
!! CHILDREN
!!      chpev,chpgv,destruction_matrix_scalapack,end_scalapack,hermit
!!      init_matrix_scalapack,init_scalapack,leave_new,matrix_from_global
!!      matrix_pzheevx,matrix_pzhegvx,matrix_to_reference,mpi_allreduce
!!      mpi_bcast,normev,status,timab,wrtout,zgemm,zhpev,zhpgv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine subdiago(cg,filstat,eig_k,evec,gsc,icg,igsc,ikpt,inonsc,istwf_k,&
&                    mcg,mgsc,mpi_enreg,nband_k,npw_k,nspinor,paral_kgb,prtvol,&
&                    subham,subovl,use_subovl,usepaw)

 use defs_basis
 use defs_datatypes
 use defs_scalapack
#if defined MPI && defined MPI2
 use mpi
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_linalg
!End of the abilint section

 implicit none
#if defined MPI && defined MPI1
 include 'mpif.h'
#endif

!Arguments ------------------------------------
 integer,intent(in) :: icg,igsc,ikpt,inonsc,istwf_k,mcg,mgsc,nband_k,npw_k
 integer,intent(in) :: nspinor,paral_kgb,prtvol,use_subovl,usepaw
 character(len=fnlen),intent(in) :: filstat
 type(MPI_type),intent(inout) :: mpi_enreg
 real(dp),intent(inout) :: subham(nband_k*(nband_k+1)),subovl(nband_k*(nband_k+1)*use_subovl)
 real(dp),intent(out) :: eig_k(nband_k),evec(2*nband_k,nband_k)
 real(dp),intent(inout) :: cg(2,mcg),gsc(2,mgsc)

!Local variables-------------------------------
 integer,parameter :: level=8
 integer :: iband,ii,ierr,iexit
 character(len=500) :: message
 real(dp) :: tsec(2)
 real(dp),allocatable :: work(:,:),zhpev1(:,:),zhpev2(:)

#if defined HAVE_SCALAPACK
!Scalapack variables----------------------------
 TYPE(matrix_scalapack)    :: sca_subham,sca_subovl,sca_evec
 TYPE(processor_scalapack) :: processor

 real(dp)        :: tmp_evec(2*nband_k,nband_k)
 INTEGER         :: communicator 
#endif

! *********************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZHPEV' :: zhpev
!DEC$ ATTRIBUTES ALIAS:'ZHPGV' :: zhpgv
!DEC$ ATTRIBUTES ALIAS:'ZHPGV' :: zgemm
#endif


!Impose Hermiticity on diagonal elements of subham (and subovl, if needed)
 call status(inonsc,filstat,iexit,level,'call hermit   ')
 call hermit(subham,subham,ierr,nband_k)
 if (use_subovl==1) call hermit(subovl,subovl,ierr,nband_k)

! Diagonalize the Hamitonian matrix
 call status(inonsc,filstat,iexit,level,'call zhpev    ')


 allocate(zhpev1(2,2*nband_k-1),zhpev2(3*nband_k-2))

# if defined T3E
! ==============

#   if defined MPI
!   ==============
       if ((mpi_enreg%paralbd <=1) .or. ((mpi_enreg%paralbd >1) .and. &
            &    (mpi_enreg%me_group==0))) then
#   endif ! FIN MPI
!   ===============

       if (use_subovl==1) then
          call CHPGV(1,'V','U',nband_k,subham,subovl,eig_k,evec,nband_k,&
               &             zhpev1,zhpev2,ierr)
       else
          call CHPEV ('V','U',nband_k,subham,eig_k,evec,nband_k,zhpev1,&
               &              zhpev2,ierr)

#      if defined MPI
!      ==============
       end if
       end if

       if (mpi_enreg%paralbd >1) then
          call timab(48,1,tsec)
          call MPI_BCAST(evec,2*nband_k*nband_k, &
               &   MPI_DOUBLE_PRECISION,0,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
          call timab(48,2,tsec)
       end if

#      endif
!   FIN MPI
!   ==================

# else
! ELSE defined T3E
! =======================

#   if defined MPI
!   ============================================

#   if defined HAVE_SCALAPACK 
!   ============================================
     if(paral_kgb == 1) then
       !call timab(570,1,tsec)

       ! ===============================
       ! INITIALISATION WORK VARIABLE 
       ! ===============================
       tmp_evec(:,:)=0._DP

       ! ============================
       ! INITIALISATION COMMUNICATOR
       ! ===========================
       if (mpi_enreg%paralbd <=1) then

          if (mpi_enreg%paral_compil_kpt==1) then
             communicator  = mpi_enreg%commcart
          else
             communicator  = MPI_COMM_WORLD
          endif

       else
          communicator  = mpi_enreg%kpt_comm(mpi_enreg%num_group)
       endif

       ! ========================
       ! INITIALISATION SCALAPACK
       ! ========================
       call init_scalapack(processor,communicator)

       ! ================================
       ! INITIALISATION SCALAPACK MATRIX
       ! ================================
       call init_matrix_scalapack(sca_subham,nband_k,nband_k,processor,10)
       call init_matrix_scalapack(sca_evec,nband_k,nband_k,processor,10)

       ! ==============================
       ! REMPLISSAGE SCALAPACK MATRIX
       ! ==============================
       call matrix_from_global(sca_subham,subham)

       if (use_subovl==1) then
          call matrix_from_global(sca_subovl,subovl)
       endif


       if (use_subovl==1) then
          ! ===========================
          ! CALL PZHEGVX FUNCTION
          ! ===========================
 !         PRINT *, 'Jutilise Scalapack1'
          call matrix_pzhegvx(processor,sca_subham,sca_subovl,sca_evec,eig_k,communicator)
       else
          ! ===========================
          ! CALL PZHEEVX FUNCTION
          ! ===========================
 !        PRINT *, 'Jutilise Scalapack2'
          call matrix_pzheevx(processor,sca_subham,sca_evec,eig_k,communicator)
       end if

       ! ==============================
       ! CONCATENATE EIGEN VECTORS
       ! ==============================
       call matrix_to_reference(sca_evec,tmp_evec)

       CALL MPI_ALLREDUCE(tmp_evec, evec, 2*nband_k*nband_k, MPI_DOUBLE_PRECISION, &
            MPI_SUM, communicator,ierr)

       ! ====================================
       ! DESTRUCTION SCALAPACK AND TMP MATRICES
       ! ====================================
       CALL destruction_matrix_scalapack(sca_subham)

       if (use_subovl==1) then
          CALL destruction_matrix_scalapack(sca_subovl)
       end if

       ! ===========================
       ! CLOSE SCALAPACK
       ! ===========================
        CALL end_scalapack(processor)

       !call timab(570,2,tsec)

     else ! paral_kgb=0

#   endif
!   END HAVE_SCALAPACK
!   =====================

        if ((mpi_enreg%paralbd <=1) .or. ((mpi_enreg%paralbd >1) .and. &
&          (mpi_enreg%me_group==0))) then

#   endif
!   END MPI
!   ===============

        if (use_subovl==1) then
           call ZHPGV(1,'V','U',nband_k,subham,subovl,eig_k,evec,nband_k,&
                &             zhpev1,zhpev2,ierr)
        else

           call ZHPEV ('V','U',nband_k,subham,eig_k,evec,nband_k,zhpev1,&
               &               zhpev2,ierr)
        end if

#    if defined MPI
!    ===============
        end if
        if (mpi_enreg%paralbd >1) then
           call timab(48,1,tsec)
           call MPI_BCAST(evec,2*nband_k*nband_k, &
                &   MPI_DOUBLE_PRECISION,0,mpi_enreg%kpt_comm(mpi_enreg%num_group),ierr)
           call timab(48,2,tsec)
        end if
#    endif
!    END MPI
!    ===============


#   if defined HAVE_SCALAPACK 
     endif ! paral_kgb
#   endif

# endif
! END defined T3E
! ========================


 deallocate(zhpev1,zhpev2)

!DEBUG
!write(6,*)' subdiago : after zhpev '
!stop
!ENDDEBUG

!Normalize each eigenvector and set phase:
 call status(inonsc,filstat,iexit,level,'call normev   ')
 call normev(evec,nband_k,nband_k)

 if(prtvol==-level)then
  write(message,'(a)')&
&  ' subdiago : iband band  evec(re:im)'
  call wrtout(06,message,'PERS')
  do iband=1,nband_k
   do ii=1,nband_k
    write(message,'(2i5,2es16.6)')&
&    iband,ii,evec(2*ii-1,iband),evec(2*ii,iband)
    call wrtout(06,message,'PERS')
   end do
  end do
 end if

 if(istwf_k==2)then
  do iband=1,nband_k
   do ii=1,nband_k
    if(abs(evec(2*ii,iband))>1.0d-10)then
     write(message,'(a,a,a,a,2i5,2es16.6,a,a)')ch10,&
&     ' subdiago : BUG ',&
&     '  For istwf_k=2, observed the following element of evec :',ch10,&
&     iband,ii,evec(2*ii-1,iband),evec(2*ii,iband),ch10,&
&     '  with a non-negligible imaginary part.'
     call wrtout(06,message,'PERS')
     call leave_new('PERS')
    end if
   end do
  end do
 end if

!Carry out rotation of bands C(G,n) according to evecs:
 call status(inonsc,filstat,iexit,level,'call zgemm   ')

! ==============================
!     SDIROT --> ZGEMM
! ==============================

 allocate(work(2,npw_k*nspinor*nband_k));work=zero

! call sdirot(cg,evec,icg,mcg,nband_k,nband_k,npw_k*nspinor)
 call zgemm('N','N',npw_k*nspinor,nband_k,nband_k,&
  &           dcmplx(1._dp), &
  &           cg(1,icg+1),npw_k*nspinor, &
  &           evec,nband_k,&
  &           dcmplx(0._dp), &
  &           work,npw_k*nspinor)
 cg(:,1+icg:npw_k*nspinor*nband_k+icg)=work(:,:)

!If paw, musb also rotate S.C(G,n):
! if (usepaw==1) call sdirot(gsc,evec,icg,mcg,nband_k,nband_k,npw_k*nspinor)
 if (usepaw==1) then
  call zgemm('N','N',npw_k*nspinor,nband_k,nband_k,&
    &           dcmplx(1._dp), &
    &           gsc(1,igsc+1),npw_k*nspinor, &
    &           evec,nband_k, &
    &           dcmplx(zero), &
    &           work,npw_k*nspinor)
  gsc(:,1+igsc:npw_k*nspinor*nband_k+igsc)=work(:,:)
 endif

 deallocate(work)

! DEBUG
!  write(6,*)' subdiago : cg(1:2) for different bands (3) '
!  do iband=1,nband_k
!   iwavef=(iband-1)*npw_k+icg
!   write(6, '(4es16.6)' )cg(1:2,1+iwavef:2+iwavef)
!  end do
! ENDDEBUG

 end subroutine subdiago
!!***
