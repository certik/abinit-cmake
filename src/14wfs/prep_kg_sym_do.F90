!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_kg_sym_do
!! NAME
!! prep_kg_sym_do
!!
!! FUNCTION
!! this routine completes de kg_k_gather vector in adding the opposite values.
!! the values are distributed on the processors in function of
!! the value of modulo(-kg_k_gather(2,i),nproc_fft)
!!
!! COPYRIGHT
!!
!! INPUTS
!!  mpi_enreg          = informations about mpi parallelization
!!  kg_k_gather        = planewave coordinates 
!!                       (of the processor + sended by other processors band)
!!  ndatarecv          = total number of values received by the processor and sended 
!!                       by the other processors band
!!
!! OUTPUT
!!  kg_k_gather_sym    = planewave coordinates     
!!                       (kg_k_gather + opposited planewave coordinates sended by the processors 
!!                       fft)
!!  ndatarecv_tot      = total number of received values by the processor
!!                       (ndatarecv   + number of received opposited planewave coordinates)
!!  ndatasend_sym      = number of sended values to the processors fft to create opposited 
!!                       planewave coordinates
!!  idatarecv0         = position of the planewave coordinates (0,0,0)
!!  tab_proc           = positions of opposited planewave coordinates in the list of the
!!                       processors fft
!!  sendcounts_sym     = number of sended values by the processor to each processor fft
!!  sendcounts_sym_all = number of sended values by each processor to the other processors fft
!!  sdispls_sym        = postions of the sended values by the processor to each processor fft
!!
!!  recvcounts_sym     = number of the received values by the processor from each processor fft
!!  recvcounts_sym_tot = number of the received values by each processor from the other processors fft
!!  rdispls_sym        = postions of the received values by the processor from each processor fft
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      prep_fourwf,prep_getghc
!!
!! CHILDREN
!!      xallgatherv_mpi,xalltoallv_mpi,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prep_kg_sym_do(mpi_enreg,&
     kg_k_gather,ndatarecv,&
     kg_k_gather_sym,ndatarecv_tot,&
     ndatasend_sym,idatarecv0,&
     tab_proc,&
     sendcounts_sym,sendcounts_sym_all,sdispls_sym,&
     recvcounts_sym,recvcounts_sym_tot,rdispls_sym)

 use defs_basis
  use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndatarecv
 integer,intent(out) :: idatarecv0,ndatarecv_tot,ndatasend_sym
 type(mpi_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: kg_k_gather(3,ndatarecv)
 integer,pointer :: kg_k_gather_sym(:,:),rdispls_sym(:),recvcounts_sym(:)
 integer,pointer :: recvcounts_sym_tot(:),sdispls_sym(:),sendcounts_sym(:)
 integer,pointer :: sendcounts_sym_all(:),tab_proc(:)

!Local variables-------------------------------
!scalars
 integer :: idatarecv,ier,iproc,jsendloc,me_band,me_fft,newspacecomm,nproc_fft
!arrays
 integer,allocatable :: kg_k_gather_send(:,:),rdispls_sym_loc(:)
 integer,allocatable :: recvcounts_sym_loc(:),sdispls_sym_loc(:)
 integer,allocatable :: sendcounts_sym_loc(:),sum_kg(:)

! *********************************************************************

!DEBUG
!write(6,*)' prep_kg_sym_do : enter '
!ENDDEBUG

!---------------------------------------------
!Desallocation
!---------------------------------------------
 if (associated(tab_proc)) then
  deallocate(tab_proc)
  deallocate(sendcounts_sym, recvcounts_sym)
  deallocate(sendcounts_sym_all, recvcounts_sym_tot)
  deallocate(sdispls_sym   , rdispls_sym)
 end if

 if (associated(kg_k_gather_sym)) deallocate(kg_k_gather_sym)

!---------------------------------------------
!Initialisation
!---------------------------------------------
 nproc_fft    = mpi_enreg%nproc_fft

 me_fft       = mpi_enreg%me_fft
 me_band      = mpi_enreg%me_band

 newspacecomm = mpi_enreg%comm_fft

!---------------------------------------------
!Allocation
!---------------------------------------------
 allocate(tab_proc(ndatarecv))

 allocate(sendcounts_sym    (nproc_fft))
 allocate(sendcounts_sym_all(nproc_fft*nproc_fft))
 allocate(sdispls_sym       (nproc_fft))

 allocate(recvcounts_sym    (nproc_fft))
 allocate(recvcounts_sym_tot(nproc_fft))
 allocate(rdispls_sym       (nproc_fft))

 allocate(sendcounts_sym_loc    (nproc_fft))
 allocate(sdispls_sym_loc       (nproc_fft))
 allocate(recvcounts_sym_loc    (nproc_fft))
 allocate(rdispls_sym_loc       (nproc_fft))

!---------------------------------------------
!Initialisation
!---------------------------------------------
 tab_proc(:)      = 0

 sendcounts_sym(:)      = 0
 sendcounts_sym_all(:)  = 0
 sdispls_sym(:)         = 0

 recvcounts_sym(:)      = 0
 recvcounts_sym_tot(:)  = 0

!---------------------------------------------
!Localisation of kg_k==[0 0 0]
!---------------------------------------------
 allocate(sum_kg(ndatarecv))
 idatarecv0    = -1
 ndatasend_sym = ndatarecv

 sum_kg=sum(abs(kg_k_gather),1)
 if (count(sum_kg==0)/=0) then
  do idatarecv=1,ndatarecv
   if (sum_kg(idatarecv)==0) idatarecv0=idatarecv
  end do
  ndatasend_sym = ndatarecv-1
 end if
 
!-----------------------------------------------------
!Localisation of the processor where the vector -k2 is
!-----------------------------------------------------

 do idatarecv=1,ndatarecv
  if (idatarecv/=idatarecv0) then
   tab_proc(idatarecv)   = modulo(-kg_k_gather(2,idatarecv),nproc_fft)
  else
   tab_proc(idatarecv) = -1
  end if
 end do

!-----------------------------------------------------------
!Calcul of the number of the values to send by the processor
!to the other processors
!------------------------------------------------------------
 do iproc=1,nproc_fft
  sendcounts_sym(iproc) = count(tab_proc(:)==(iproc-1))
 end do


!------------------------------------------------------------
!Saving sendcounts_sym for each processor in sendcounts_sym_all 
!knowed by all processors of newspacecomm
!-------------------------------------------------------------
 rdispls_sym(1)=0
 do iproc=2,nproc_fft
  rdispls_sym(iproc)= nproc_fft*(iproc-1)
 end do

 recvcounts_sym(:)=nproc_fft

 call xallgatherv_mpi(sendcounts_sym(:)    ,nproc_fft,&
 sendcounts_sym_all(:),recvcounts_sym,rdispls_sym,&
 newspacecomm,ier)

!-------------------------------------------------------------
!Calcul of the dimension of kg_k_gather_sym for each processor
!recvcounts_sym_tot is knowed by all processors of newspacecomm
!-------------------------------------------------------------
 call xsum_mpi(sendcounts_sym,recvcounts_sym_tot,nproc_fft,newspacecomm,ier)

!----------------------------------
!Dimension of kg_k_gather_sym
!----------------------------------
 ndatarecv_tot = ndatarecv+recvcounts_sym_tot(me_fft+1)

!-----------------------------------------------------
!Intialisation kg_k_gather_sym
!-----------------------------------------------------
 allocate(kg_k_gather_sym(3,ndatarecv_tot))

 kg_k_gather_sym(:,:)=0
 kg_k_gather_sym(:,1:ndatarecv) = kg_k_gather(:,:)

!---------------------------------------------------
!Allocation and initialisation
!---------------------------------------------------
 allocate(kg_k_gather_send(3,ndatasend_sym))
 
 kg_k_gather_send(:,:)=0

!---------------------------------------------------
!The values is sorted in blocks 
!---------------------------------------------------
 jsendloc=0

 do iproc=1,nproc_fft

! Calcul of the position of the begin of the block
! ------------------------------------------------
  sdispls_sym(iproc)=jsendloc
  
! Creation of the blocks
! ----------------------
  do idatarecv=1,ndatarecv
   if (tab_proc(idatarecv)==(iproc-1)) then 
    jsendloc=jsendloc+1
    kg_k_gather_send(:,jsendloc)  = -kg_k_gather(:,idatarecv)
   end if
  end do
 end do
 
!---------------------------------------- 
!Calcul of the position of received datas
!----------------------------------------
 rdispls_sym(1)= ndatarecv
 recvcounts_sym(1)= sendcounts_sym_all((me_fft+1))
 
 do iproc=2,nproc_fft
  rdispls_sym(iproc)    = rdispls_sym(iproc-1) + &
  sendcounts_sym_all((me_fft+1)+(iproc-2)*nproc_fft)
  recvcounts_sym(iproc) = sendcounts_sym_all((me_fft+1)+(iproc-1)*nproc_fft)
 end do
 
!-----------------------
!Exchange of kg_k
!-----------------------
 sendcounts_sym_loc = sendcounts_sym*3
 sdispls_sym_loc    = sdispls_sym   *3
 recvcounts_sym_loc = recvcounts_sym*3
 rdispls_sym_loc    = rdispls_sym   *3


 call xalltoallv_mpi(kg_k_gather_send(:,:),sendcounts_sym_loc,sdispls_sym_loc,&
 kg_k_gather_sym(:,:) ,recvcounts_sym_loc,rdispls_sym_loc,&
 newspaceComm,ier)

!-----------------------
!Desallocation
!-----------------------
 deallocate(sendcounts_sym_loc, recvcounts_sym_loc)
 deallocate(sdispls_sym_loc   , rdispls_sym_loc)
 deallocate(kg_k_gather_send  , sum_kg)

end subroutine prep_kg_sym_do
!!***
