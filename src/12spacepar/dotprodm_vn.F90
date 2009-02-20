!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotprodm_vn
!! NAME
!! dotprodm_vn
!!
!!
!! FUNCTION
!! For a set of densities and a set of potentials,
!! compute the dot product (integral over FFT grid) of each pair, to obtain
!! a series of energy-like quantity (so the usual dotproduct is divided
!! by the number of FFT points, and multiplied by the primitive cell volume).
!! Take into account the spin components of the density and potentials (nspden),
!! and sum correctly over them. Note that the storage of densities and
!! potentials is different : for potential, one stores the matrix components,
!! while for the density, one stores the trace, and then, either the
!! spin-polarisation (if nspden=2), or the magnetisation vector (if nspden=4).
!! Need the index of the first density/potential pair to be treated, in each array,
!! and the number of pairs to be treated.
!! Might be used to compute just one dot product, in
!! a big array, such as to avoid copying the density and potential from a big array
!! to a temporary place.
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
!!  cpldot=if 1, the dot array is real, if 2, the dot array is complex (not coded yet for nspden=4)
!!  denarr(cplex*nfft,nspden,nden)=real space density on FFT grid
!!  id=index of the first density to be treated in the denarr array
!!  ip=index of the first potential to be treated in the potarr array
!!  mpi_enreg=informations about MPI parallelization
!!  multd=number of densities to be treated
!!  multp=number of potentials to be treated
!!  nden=third dimension of the denarr array
!!  nfft= (effective) number of FFT grid points (for this processor)
!!  nfftot= total number of FFT grid points
!!  npot=third dimension of the potarr array
!!  nspden=number of spin-density components
!!  potarr(cplex*nfft,nspden,npot)=real space potential on FFT grid
!!                 (will be complex conjugated if cplex=2 and cpldot=2)
!!  ucvol=unit cell volume (Bohr**3)
!!
!! OUTPUT
!!  dot(cpldot,multp,multd)= series of values of the dot product potential/density
!!
!! SIDE EFFECTS
!!
!!
!! NOTES
!!
!!
!! PARENTS
!!      aprxdr
!!
!! CHILDREN
!!      contract_dp_ge_val,contract_int_ge_val,contract_int_list,leave_new
!!      timab,wrtout,xcomm_init,xsum_mpi
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine dotprodm_vn(cplex,cpldot,denarr,dot,id,ip,mpi_enreg,multd,multp,&
& nden,nfft,nfftot,npot,nspden,potarr,ucvol)

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
 integer,intent(in) :: cpldot,cplex,id,ip,multd,multp,nden,nfft,nfftot,npot
 integer,intent(in) :: nspden
 real(dp),intent(in) :: ucvol
 type(MPI_type) :: mpi_enreg
!arrays
 real(dp),intent(in) :: denarr(cplex*nfft,nspden,nden)
 real(dp),intent(in) :: potarr(cplex*nfft,nspden,npot)
 real(dp),intent(out) :: dot(cpldot,multp,multd)

!Local variables-------------------------------
!scalars
 integer :: i1,i2,ierr,ir,old_paral_level,spaceComm
 real(dp) :: ai,ar,dim_dn,dim_up,dre_dn,dre_up,factor,pim_dn,pim_up,pre_dn
 real(dp) :: pre_up
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
!no_abirules
#if defined CONTRACT
 character(len=11) :: subrnm
#endif

! *************************************************************************

#if defined CONTRACT
 subrnm='dotprodm_vn'
!Real or complex inputs are coded
 call contract_int_list(subrnm,'cplex',cplex,(/1,2/),2)
!Real or complex outputs are coded
 call contract_int_list(subrnm,'cpldot',cpldot,(/1,2/),2)
 call contract_int_ge_val(subrnm,'id',id,1)
 call contract_int_ge_val(subrnm,'ip',ip,1)
 call contract_int_ge_val(subrnm,'multd',multd,1)
 call contract_int_ge_val(subrnm,'multp',multp,1)
 call contract_int_ge_val(subrnm,'nden',nden,1)
 call contract_int_ge_val(subrnm,'nfft',nfft,1)
 call contract_int_ge_val(subrnm,'nfftot',nfftot,1)
 call contract_int_ge_val(subrnm,'npot',npot,1)
 call contract_int_list(subrnm,'nspden',nspden,(/1,2,4/),3)
 call contract_dp_ge_val(subrnm,'ucvol',ucvol,zero)
 call contract_int_ge_val(subrnm,'nden-id-multd',nden-id-multd,-1)
 call contract_int_ge_val(subrnm,'npot-ip-multp',npot-ip-multp,-1)
#endif

!Non-collinear and complex outputs is not coded
 if(cpldot==2 .and. nspden==4)then
  write(message,'(a,a,a,a)') ch10,&
&  ' dotprodm_vn : BUG -',ch10,&
&  '  The argument cpldot and nspden cannot be together cpldot=2 and nspden=4'
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 if(nspden==1)then

  if(cpldot==1 .or. cplex==1 )then

   do i2=1,multd
    do i1=1,multp
     ar=zero
!    $OMP PARALLEL DO PRIVATE(ir) &
!    $OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar)
     do ir=1,cplex*nfft
      ar=ar + potarr(ir,1,ip+i1-1)*denarr(ir,1,id+i2-1)
     end do
!    $OMP END PARALLEL DO
     dot(1,i1,i2)=ar
    end do ! i1
   end do ! i2

  else  ! cpldot==2 and cplex==2 : one builds the imaginary part, from complex den/pot

   do i2=1,multd
    do i1=1,multp
     ar=zero ; ai=zero
!    $OMP PARALLEL DO PRIVATE(ir) &
!    $OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar,ai)
     do ir=1,nfft
      ar=ar + potarr(2*ir-1,1,ip+i1-1)*denarr(2*ir-1,1,id+i2-1) &
&      + potarr(2*ir  ,1,ip+i1-1)*denarr(2*ir  ,1,id+i2-1)
      ai=ai + potarr(2*ir-1,1,ip+i1-1)*denarr(2*ir  ,1,id+i2-1) &
&      - potarr(2*ir  ,1,ip+i1-1)*denarr(2*ir-1,1,id+i2-1)
     end do
!    $OMP END PARALLEL DO
     dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
    end do ! i1
   end do ! i2

  end if

 else if(nspden==2)then

  if(cpldot==1 .or. cplex==1 )then

   do i2=1,multd
    do i1=1,multp
     ar=zero
!    $OMP PARALLEL DO PRIVATE(ir) &
!    $OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar)
     do ir=1,cplex*nfft
      ar=ar + potarr(ir,1,ip+i1-1)* denarr(ir,2,id+i2-1)               &    ! This is the spin up contribution
&     + potarr(ir,2,ip+i1-1)*(denarr(ir,1,id+i2-1)-denarr(ir,2,id+i2-1)) ! This is the spin down contribution
     end do
!    $OMP END PARALLEL DO
     dot(1,i1,i2)=ar
    end do ! i1
   end do ! i2

  else ! cpldot==2 and cplex==2 : one builds the imaginary part, from complex den/pot

   do i2=1,multd
    do i1=1,multp
     ar=zero ; ai=zero
!    $OMP PARALLEL DO PRIVATE(ir) &
!    $OMP&PRIVATE(dre_up,dim_up,dre_dn,dim_dn) &
!    $OMP&PRIVATE(pre_up,pim_up,pre_dn,pim_dn) &
!    $OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar,ai)
     do ir=1,nfft

      dre_up=denarr(2*ir-1,2,id+i2-1)
      dim_up=denarr(2*ir  ,2,id+i2-1)
      dre_dn=denarr(2*ir-1,1,id+i2-1)-dre_up
      dim_dn=denarr(2*ir  ,1,id+i2-1)-dim_up

      pre_up=potarr(2*ir-1,1,ip+i1-1)
      pim_up=potarr(2*ir  ,1,ip+i1-1)
      pre_dn=potarr(2*ir-1,2,ip+i1-1)
      pim_dn=potarr(2*ir  ,2,ip+i1-1)

      ar=ar + pre_up * dre_up &
&      + pim_up * dim_up &
&      + pre_dn * dre_dn &
&      + pim_dn * dim_dn
      ai=ai + pre_up * dim_up &
&      - pim_up * dre_up &
&      + pre_dn * dim_dn &
&      - pim_dn * dre_dn

     end do
!    $OMP END PARALLEL DO
     dot(1,i1,i2)=ar ; dot(2,i1,i2)=ai
    end do ! i1
   end do ! i2

  end if

 else if(nspden==4)then

! Not yet coded for cpldot=2
! Even the case cplex=2 should be checked
! \rho{\alpha,\beta} V^{\alpha,\beta} =
! rho*(V^{11}+V^{22})/2$
! + m_x Re(V^{12})- m_y Im{V^{12}}+ m_z(V^{11}-V^{22})/2

  do i2=1,multd
   do i1=1,multp
    ar=zero
!   $OMP PARALLEL DO PRIVATE(ir) &
!   $OMP&SHARED(id,i1,i2,ip,nfft,denarr,potarr) REDUCTION(+:ar)
    do ir=1,cplex*nfft
     ar=ar+(potarr(ir,1,ip+i1-1)+potarr(ir,2,ip+i1-1))*half*denarr(ir,1,id+i2-1)& ! This is the density contrib
&    + potarr(ir,3,ip+i1-1)                           *denarr(ir,2,id+i2-1)& ! This is the m_x contrib
&    - potarr(ir,4,ip+i1-1)                           *denarr(ir,3,id+i2-1)& ! This is the m_y contrib
&    +(potarr(ir,1,ip+i1-1)-potarr(ir,2,ip+i1-1))*half*denarr(ir,4,id+i2-1)  ! This is the m_z contrib
    end do
!   $OMP END PARALLEL DO
    dot(1,i1,i2)=ar
   end do ! i1
  end do ! i2

 end if ! nspden

 factor=ucvol/dble(nfftot)
 dot(:,:,:)=factor*dot(:,:,:)

!XG030513 : MPIWF reduction (addition) on dot is needed here
!Init mpi_comm
 if(mpi_enreg%paral_compil_fft==1)then
  old_paral_level=mpi_enreg%paral_level
  mpi_enreg%paral_level=3
  call xcomm_init(mpi_enreg,spaceComm)
  call timab(48,1,tsec)

! PATCH dotprodm_vn spacecomm --> comm_fft
  if(mpi_enreg%mode_para=='b') spaceComm=mpi_enreg%comm_fft

  call xsum_mpi(dot,spaceComm ,ierr)
  call timab(48,2,tsec)
  mpi_enreg%paral_level=old_paral_level
 end if

 if(cpldot==2 .and. cplex==1)dot(2,:,:)=zero

end subroutine dotprodm_vn
!!***
