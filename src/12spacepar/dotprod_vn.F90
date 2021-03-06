!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprod_vn
!! NAME
!! dotprod_vn
!!
!!
!! FUNCTION
!! Compute dot product of potential and density (integral over FFT grid), to obtain
!! an energy-like quantity (so the usual dotproduct is divided
!! by the number of FFT points, and multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden),
!! and sum correctly over them. Note that the storage of densities and
!! potentials is different : for potential, one stores the matrix components,
!! while for the density, one stores the trace, and then, either the
!! spin-polarisation (if nspden=2), or the magnetisation vector (if nspden=4).
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex=if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  dens(cplex*nfft,nspden)=real space density on FFT grid
!!  mpi_enreg=informations about MPI parallelization
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nfftot= total number of FFT grid points
!!  nspden=number of spin-density components
!!  option= if 1, only the real part is computed
!!          if 2, both real and imaginary parts are computed  (not yet coded)
!!  pot(cplex*nfft,nspden)=real space potential on FFT grid
!!                 (will be complex conjugated if cplex=2 and option=2)
!!  ucvol=unit cell volume (Bohr**3)
!!
!! OUTPUT
!!  doti= imaginary part of the dot product, output only if option=2 (and cplex=2).
!!  dotr= real part
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      dyxc13,eltfrxc3,energy,ftfvw1,ftfvw2,loop3dte,newfermie1,nselt3,nstdy3
!!      odamix,prcrskerker2,prctfvw1,prctfvw2,prctfw3,rhohxc,rhotov,rhotov3
!!      setvtr
!!
!! CHILDREN
!!      contract_dp_ge_val,contract_int_ge_val,contract_int_list,leave_new
!!      timab,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dotprod_vn(cplex,dens,dotr,doti,mpi_enreg,nfft,nfftot,nspden,option,pot,ucvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11contract
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nfftot,nspden,option
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: doti,dotr
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: dens(cplex*nfft,nspden),pot(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ierr,ifft,old_paral_level,spaceComm
 real(dp) :: dim_dn,dim_up,dre_dn,dre_up,factor,pim_dn,pim_up,pre_dn,pre_up
 character(len=500) :: message
!arrays
 real(dp) :: buffer2(2),tsec(2)
!no_abirules
#if defined CONTRACT
 character(len=10) :: subrnm
#endif

! *************************************************************************

#if defined CONTRACT
 subrnm='dotprod_vn'
!Real or complex inputs are coded
 call contract_int_list(subrnm,'cplex',cplex,(/1,2/),2)
 call contract_int_ge_val(subrnm,'nfft',nfft,1)
 call contract_int_ge_val(subrnm,'nfftot',nfftot,1)
 call contract_int_list(subrnm,'nspden',nspden,(/1,2,4/),3)
!Real or complex output are coded
 call contract_int_list(subrnm,'option',option,(/1,2/),2)
 call contract_dp_ge_val(subrnm,'ucvol',ucvol,zero)
#endif
!Non-collinear and complex outputs is not coded
 if(option==2 .and. nspden==4)then
  write(message,'(a,a,a,a)') ch10,&
&  ' dotprod_vn : BUG -',ch10,&
&  '  The argument option and nspden cannot be together option=2 and nspden=4'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 dotr=zero ; doti=zero

 if(nspden==1)then

  if(option==1 .or. cplex==1 )then
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(nfft,dens,pot) REDUCTION(+:dotr)
   do ifft=1,cplex*nfft
    dotr=dotr + pot(ifft,1)*dens(ifft,1)
   end do
!  $OMP END PARALLEL DO

  else  ! option==2 and cplex==2 : one builds the imaginary part, from complex den/pot

!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(nfft,dens,pot) REDUCTION(+:dotr,doti)
   do ifft=1,nfft
    dotr=dotr + pot(2*ifft-1,1)*dens(2*ifft-1,1) &
&    + pot(2*ifft  ,1)*dens(2*ifft  ,1)
    doti=doti + pot(2*ifft-1,1)*dens(2*ifft  ,1) &
&    - pot(2*ifft  ,1)*dens(2*ifft-1,1)
   end do
!  $OMP END PARALLEL DO

  end if

 else if(nspden==2)then

  if(option==1 .or. cplex==1 )then
!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&SHARED(nfft,dens,pot) REDUCTION(+:dotr)
   do ifft=1,cplex*nfft
    dotr=dotr + pot(ifft,1)* dens(ifft,2)               &    ! This is the spin up contribution
&   + pot(ifft,2)*(dens(ifft,1)-dens(ifft,2))      ! This is the spin down contribution
   end do
!  $OMP END PARALLEL DO

  else ! option==2 and cplex==2 : one builds the imaginary part, from complex den/pot

!  $OMP PARALLEL DO PRIVATE(ifft) &
!  $OMP&PRIVATE(dre_up,dim_up,dre_dn,dim_dn) &
!  $OMP&PRIVATE(pre_up,pim_up,pre_dn,pim_dn) &
!  $OMP&SHARED(nfft,dens,pot) REDUCTION(+:dotr,doti)
   do ifft=1,nfft

    dre_up=dens(2*ifft-1,2)
    dim_up=dens(2*ifft  ,2)
    dre_dn=dens(2*ifft-1,1)-dre_up
    dim_dn=dens(2*ifft  ,1)-dim_up

    pre_up=pot(2*ifft-1,1)
    pim_up=pot(2*ifft  ,1)
    pre_dn=pot(2*ifft-1,2)
    pim_dn=pot(2*ifft  ,2)

    dotr=dotr + pre_up * dre_up &
&    + pim_up * dim_up &
&    + pre_dn * dre_dn &
&    + pim_dn * dim_dn
    doti=doti + pre_up * dim_up &
&    - pim_up * dre_up &
&    + pre_dn * dim_dn &
&    - pim_dn * dre_dn

   end do
!  $OMP END PARALLEL DO


  end if

 else if(nspden==4)then

! Not yet coded for option=2
! Even the case cplex=2 should be checked
! From GZerah
! \rho{\alpha,\beta} V^{\alpha,\beta} =
! rho*(V^{11}+V^{22})/2$
! + m_x Re(V^{12})+ m_y Im{V^{12}}+ m_z(V^{11}-V^{22})/2

! MT071212 : I agree with XG !
! This result from the computation of the trace of \rho.V
! which is \rho{\alpha,\beta} V^{\beta,\alpha}.
! Now the expression is changed !

! $OMP PARALLEL DO PRIVATE(ifft) &
! $OMP&SHARED(nfft,dens,pot) REDUCTION(+:dotr)
  do ifft=1,cplex*nfft
   dotr=dotr + (pot(ifft,1)+pot(ifft,2))*half*dens(ifft,1) &   ! This is the density contrib
&  +  pot(ifft,3)                  *dens(ifft,2) &   ! This is the m_x contrib
&  -  pot(ifft,4)                  *dens(ifft,3) &   ! This is the m_y contrib
&  + (pot(ifft,1)-pot(ifft,2))*half*dens(ifft,4)     ! This is the m_z contrib
  end do
! $OMP END PARALLEL DO

 end if ! nspden

 factor=ucvol/dble(nfftot)
 dotr=factor*dotr
 doti=factor*doti

!XG030513 : MPIWF reduction (addition) on dotr, doti is needed here
!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
! DC070608: can't understand the micmac to obtain the communicator.
! Four stupid lines to obtain the right comminucator
! instead of using the designed routine xcomm_init()...
! Too bad.
  old_paral_level=mpi_enreg%paral_level
  mpi_enreg%paral_level=3
  call xcomm_init(mpi_enreg,spaceComm)
  if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft

  buffer2(1)=dotr
  buffer2(2)=doti
  call timab(48,1,tsec)
  call xsum_mpi(buffer2,spaceComm ,ierr)
! call xsum_mpi(dotr,spaceComm ,ierr)
! call xsum_mpi(doti,spaceComm ,ierr)
  call timab(48,2,tsec)
  dotr=buffer2(1)
  doti=buffer2(2)
  mpi_enreg%paral_level=old_paral_level
 end if

end subroutine dotprod_vn
!!***
