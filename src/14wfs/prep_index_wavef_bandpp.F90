!{\src2tex{textfont=tt}}
!!****f* ABINIT/prep_index_wavef_bandpp
!! NAME
!! prep_index_wavef_bandpp
!!
!! FUNCTION
!! this routine sorts the waves functions by bandpp and by processors 
!! after the alltoall
!!
!! COPYRIGHT
!!
!! INPUTS
!!  nproc_band = number of processors below the band
!!  bandpp     = number of groups of couple of waves functions
!!  nspinor    = number of spin
!!  ndatarecv  = total number of values received by the processor and sended 
!!               by the other processors band
!!  recvcounts = number of values sended by each processor band and received 
!!               by the processor
!!  rdispls    = positions of the values received by the processor and  
!!               sended by each processor band
!!
!! OUTPUT
!!  index_wavef_band = position of the sorted values
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      prep_fourwf,prep_getghc,prep_nonlop
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prep_index_wavef_bandpp(nproc_band,bandpp,&
                             nspinor,ndatarecv,&
                             recvcounts,rdispls,&
                             index_wavef_band)

 use defs_basis
  use defs_datatypes  

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: bandpp,ndatarecv,nproc_band,nspinor
!arrays
 integer,intent(in) :: rdispls(nproc_band),recvcounts(nproc_band)
 integer,pointer :: index_wavef_band(:)

!Local variables-------------------------------
!scalars
 integer :: idebc,idebe,ifinc,ifine,iindex,iproc,kbandpp

! *********************************************************************

!DEBUG
!write(6,*)' prep_index_wavef_banpp : enter '
!ENDDEBUG


!---------------------------------------------
!Allocation
!---------------------------------------------
 allocate(index_wavef_band (bandpp*nspinor*ndatarecv))
 index_wavef_band(:)   =0
 
!---------------------------------------------
!Calcul : loops on bandpp and processors band
!---------------------------------------------
 do kbandpp=1,bandpp
  
  do iproc=1,nproc_band
   
   idebe = (rdispls(iproc) + 1)  + (kbandpp-1) * ndatarecv
   ifine = idebe + recvcounts(iproc) -1

   if (iproc==1) then
    idebc =   (kbandpp-1)* recvcounts(iproc) &
    + 1  
   else
    idebc =   (bandpp)  * sum(recvcounts(1:iproc-1)) &
    + (kbandpp-1)* recvcounts(iproc) &
    + 1 
   end if
   ifinc = idebc + recvcounts(iproc) -1


   index_wavef_band(idebe:ifine) = (/( iindex,iindex=idebc,ifinc)/)
  end do
 end do

end subroutine prep_index_wavef_bandpp
!!***
