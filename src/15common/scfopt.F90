!{\src2tex{textfont=tt}}
!!****f* ABINIT/scfopt
!!
!! NAME
!! scfopt
!!
!! FUNCTION
!! Compute the next vtrial of the SCF cycle.
!! Possible algorithms are : simple mixing, Anderson (order 1 or 2), Pulay
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG,GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cplex= if 1, real space functions on FFT grid are REAL, if 2, COMPLEX
!!  dtn_pc(3,natom)=preconditioned change of atomic position (in Bohr)
!!  iscf= 2 => simple mixing
!!      = 3,4 => Anderson mixing
!!      = 7 => Pulay mixing
!!  istep= number of the step in the SCF cycle
!!  move_atm= if 1, then new atomic positions as well as a new potential
!!      must be computed; otherwise, only the new potential must be computed.
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  nspden=number of spin-density components
!!  n_fftgr=third dimension of the array f_fftgr
!!  n_index=dimension for indices of potential/density (see ivrespc, i_vtrial...)
!!  pawoptmix= - PAW only - 1 if the computed residuals include the PAW (rhoij) part
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  vtrial(cplex*nfft,nspden)= at input, it is the trial potential that gave
!!     the input preconditioned residual potential
!!     at output, it is the new trial potential .
!!  f_fftgr(cplex*nfft,nspden,n_fftgr)=different functions defined on the fft grid :
!!   The input vtrial is transferred, at output,in f_fftgr(:,:,i_vtrial(1)).
!!   The old vtrial is transferred, at output,in f_fftgr(:,:,i_vtrial(2)).
!!   The input preconditioned residual potential is in f_fftgr(:,:,i_vrespc(1))
!!   Two input old preconditioned residual potentials in f_fftgr(:,:,i_vrespc(2)) and f_fftgr(:,:,i_vrespc(3))
!!    Before output a permutation of i_vrespc(1), i_vrespc(2) and i_vrespc(3) occurs, without
!!    actually copying all the data (change of pointer).
!!  xred(3,natom)=at input, reduced dimensionless atomic coordinates
!!      that induced dtn_pc,
!!              at output, the new atomic coordinates
!!  i_vrespc(n_index)=index of the preconditioned residual potentials (present and past) in the array f_fftgr
!!  i_vtrial(n_index)  =indices of the potential (present and past) in the array f_fftgr
!!  ==== if usepaw==1
!!    f_paw(npawmix,n_fftgr*mffmem*usepaw)=different functions used for PAW
!!                                           (same as f_fftgr but for spherical part)
!!    vpaw(npawmix*usepaw)=at input, the aug. occupancies (rhoij) that gave
!!                               the input preconditioned residual potential
!!                           at output, it is the new aug. occupancies.
!!
!! PARENTS
!!      newrho,newvtr,newvtr3
!!
!! CHILDREN
!!      dgetrf,dgetri,dotprodm_v,leave_new,sqnormm_v,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine scfopt(cplex,dtn_pc,f_fftgr,f_paw,iscf,istep,i_vrespc,i_vtrial,move_atm,mpi_enreg,&
&                  natom,nfft,npawmix,nspden,n_fftgr,n_index,pawoptmix,usepaw,vpaw,vtrial,xred)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12spacepar
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,iscf,istep,move_atm,n_fftgr,n_index,natom,nfft
 integer,intent(in) :: npawmix,nspden,pawoptmix,usepaw
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(inout) :: i_vrespc(n_index),i_vtrial(n_index)
 real(dp),intent(in) :: dtn_pc(3,natom)
 real(dp),intent(inout) :: f_fftgr(cplex*nfft,nspden,n_fftgr)
 real(dp),intent(inout) :: f_paw(npawmix,n_fftgr*usepaw),vpaw(npawmix*usepaw)
 real(dp),intent(inout) :: vtrial(cplex*nfft,nspden),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer,parameter :: npulaymax=50
 integer :: i_vstore,ierr,ifft,ii,index,isp,istore,jj,niter,npulay,tmp
 real(dp),save :: prod_resid_old,resid_old,resid_old2
 real(dp) :: aa1,aa2,bb,cc1,cc2,current,det,lambda,lambda2,resid_best
 character(len=500) :: message
!arrays
 integer,allocatable :: ipiv(:)
 real(dp),save :: amat(npulaymax+1,npulaymax+1)
 real(dp) :: prod_resid(1),prod_resid2(1),resid_new(1)
 real(dp),allocatable :: alpha(:),amatinv(:,:),rwork(:)

! *************************************************************************

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'DGETRF' :: dgetrf
!DEC$ ATTRIBUTES ALIAS:'DGETRI' :: dgetri
#endif

!DEBUG
!write(6,*)' scfopt : enter ; istep,iscf,move_atm ',istep,iscf,move_atm
!ENDDEBUG

 i_vstore=i_vtrial(1)
 if (iscf==4) i_vstore=i_vtrial(2)
 if (iscf==7) then
  npulay=(n_fftgr-2)/2
  i_vstore=i_vtrial(npulay)
 end if

!_______________________________________________________________
!Here use only the preconditioning, or initialize the other algorithms

 if(istep==1 .or. iscf==2)then

  write(message,'(2a)') ch10,'Simple mixing update:'
  call wrtout(6,message,'COLL')

! Compute the new residual resid_new, from f_fftgr/f_paw(:,:,i_vrespc(1))
  call sqnormm_v(cplex,i_vrespc(1),mpi_enreg,1,nfft,resid_new,n_fftgr,nspden,f_fftgr)
  if (usepaw==1.and.pawoptmix==1) then
   do index=1,npawmix
    resid_new(1)=resid_new(1)+f_paw(index,i_vrespc(1))**2
   end do
  end if
  write(message,*)' residual square of the potential :',resid_new(1)
  call wrtout(6,message,'COLL')

! Store information for later use
  if (iscf==3.or.iscf==4) resid_old=resid_new(1)
  if (iscf==7) amat(1,1)=resid_new(1)

! Compute new vtrial (and new rhoij if PAW)
  if (iscf/=2) f_fftgr(:,:,i_vstore)=vtrial(:,:)
  vtrial(:,:)=vtrial(:,:)+f_fftgr(:,:,i_vrespc(1))
  if (usepaw==1) then
   if (iscf/=2) f_paw(:,i_vstore)=vpaw(:)
   vpaw(:)=vpaw(:)+f_paw(:,i_vrespc(1))
  end if

! Change atomic positions
  if(move_atm==1)then
   if(iscf==2)then
    xred(:,:)=xred(:,:)+dtn_pc(:,:)
   else
    write(message, '(a)' ) &
&    ' scfopt : ERROR -',ch10,&
&    '  DEVELOP - The only iscf allowed with move_atm==1 is iscf=2 .'
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if
  end if

! _______________________________________________________________
! Here Anderson algorithm using one previous iteration
 else if((istep==2 .or. iscf==3).and.iscf/=7)then

  write(message,'(2a)') ch10,'Anderson update:'
  call wrtout(6,message,'COLL')

! Compute resid_new from f_fftgr/f_paw(:,:,i_vrespc(1))
  call sqnormm_v(cplex,i_vrespc(1),mpi_enreg,1,nfft,resid_new,n_fftgr,nspden,f_fftgr)
  if (usepaw==1.and.pawoptmix==1) then
   do index=1,npawmix
    resid_new(1)=resid_new(1)+f_paw(index,i_vrespc(1))**2
   end do
  end if
  write(message,*)' residual square of the potential: ',resid_new(1)
  call wrtout(6,message,'COLL')

! Compute prod_resid from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(2))
  call dotprodm_v(cplex,1,prod_resid,i_vrespc(1),i_vrespc(2),mpi_enreg,1,1,&
&  nfft,n_fftgr,n_fftgr,nspden,f_fftgr,f_fftgr)
  if (usepaw==1.and.pawoptmix==1) then
   do index=1,npawmix
    prod_resid(1)=prod_resid(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(2))
   end do
  end if

! Compute mixing factor
  lambda=(resid_new(1)-prod_resid(1))/(resid_new(1)+resid_old-2*prod_resid(1))
  write(message,*)' mixing of old trial potential    :',lambda
  call wrtout(6,message,'COLL')

! Evaluate best residual square on the line
  resid_best=(1.0_dp-lambda)*(1.0_dp-lambda)*resid_new(1)&
&  +(1.0_dp-lambda)*lambda        *2*prod_resid(1)&
&  +lambda        *lambda        *resid_old
  write(message,*)' predicted best residual square on the line: ',resid_best
  call wrtout(6,message,'COLL')

! Store information for later use
  if (iscf==4) then
   prod_resid_old=prod_resid(1)
   resid_old2=resid_old
  end if
  resid_old=resid_new(1)

! Save latest trial potential and compute new trial potential
  do isp=1,nspden
   do ifft=1,cplex*nfft
    current=vtrial(ifft,isp)
    vtrial(ifft,isp)=(one-lambda)*(current                      +f_fftgr(ifft,isp,i_vrespc(1)))&
&    +lambda      *(f_fftgr(ifft,isp,i_vtrial(1))+f_fftgr(ifft,isp,i_vrespc(2)))
    f_fftgr(ifft,isp,i_vstore)=current
   end do
  end do

! PAW: save latest rhoij and compute new rhoij
  do index=1,npawmix
   current=vpaw(index)
   vpaw(index)=(one-lambda)*(current                 +f_paw(index,i_vrespc(1)))&
&   +lambda      *(f_paw(index,i_vtrial(1))+f_paw(index,i_vrespc(2)))
   f_paw(index,i_vstore)=current
  end do

! _______________________________________________________________
! Here Anderson algorithm using two previous iterations
 else if(iscf==4.and.iscf/=7)then

  write(message,'(2a)') ch10,'Anderson (order 2) update:'
  call wrtout(6,message,'COLL')

! Compute resid_new from f_fftgr/f_paw(:,:,i_vrespc(1))
  call sqnormm_v(cplex,i_vrespc(1),mpi_enreg,1,nfft,resid_new,n_fftgr,nspden,f_fftgr)
  if (usepaw==1.and.pawoptmix==1) then
   do index=1,npawmix
    resid_new(1)=resid_new(1)+f_paw(index,i_vrespc(1))**2
   end do
  end if
  write(message,*)' residual square of the potential :',resid_new(1)
  call wrtout(6,message,'COLL')

! Compute prod_resid from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(2))
  call dotprodm_v(cplex,1,prod_resid,i_vrespc(1),i_vrespc(2),mpi_enreg,1,1,&
&  nfft,n_fftgr,n_fftgr,nspden,f_fftgr,f_fftgr)
  if (usepaw==1.and.pawoptmix==1) then
   do index=1,npawmix
    prod_resid(1)=prod_resid(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(2))
   end do
  end if

! Compute prod_resid2 from f_fftgr/f_paw(:,:,i_vrespc(1)) and f_fftgr/f_paw(:,:,i_vrespc(3))
  call dotprodm_v(cplex,1,prod_resid2,i_vrespc(1),i_vrespc(3),mpi_enreg,1,1,&
&  nfft,n_fftgr,n_fftgr,nspden,f_fftgr,f_fftgr)
  if (usepaw==1.and.pawoptmix==1) then
   do index=1,npawmix
    prod_resid2(1)=prod_resid2(1)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(3))
   end do
  end if

! Compute mixing factors
  aa1=resid_new(1)+resid_old -two*prod_resid (1)
  aa2=resid_new(1)+resid_old2-two*prod_resid2(1)
  bb =resid_new(1)+prod_resid_old-prod_resid(1)-prod_resid2(1)
  cc1=resid_new(1)-prod_resid (1)
  cc2=resid_new(1)-prod_resid2(1)
  det=aa1*aa2-bb*bb
  lambda =(aa2*cc1-bb*cc2)/det
  lambda2=(aa1*cc2-bb*cc1)/det
  write(message,*)' mixing of old trial potentials   :',lambda,lambda2
  call wrtout(6,message,'COLL')

! Store information for later use
  prod_resid_old=prod_resid(1)
  resid_old2=resid_old
  resid_old=resid_new(1)

! Save latest trial potential and compute new trial potential
  do isp=1,nspden
   do ifft=1,cplex*nfft
    current=vtrial(ifft,isp)
    vtrial(ifft,isp)=&
&    (one-lambda-lambda2)*(current                      +f_fftgr(ifft,isp,i_vrespc(1)))&
&    +lambda             *(f_fftgr(ifft,isp,i_vtrial(1))+f_fftgr(ifft,isp,i_vrespc(2)))&
&    +lambda2            *(f_fftgr(ifft,isp,i_vtrial(2))+f_fftgr(ifft,isp,i_vrespc(3)))
    f_fftgr(ifft,isp,i_vstore)=current
   end do
  end do

! PAW: save latest rhoij and compute new rhoij
  do index=1,npawmix
   current=vpaw(index)
   vpaw(index)=&
&   (one-lambda-lambda2)*(current                 +f_paw(index,i_vrespc(1)))&
&   +lambda             *(f_paw(index,i_vtrial(1))+f_paw(index,i_vrespc(2)))&
&   +lambda2            *(f_paw(index,i_vtrial(2))+f_paw(index,i_vrespc(3)))
   f_paw(index,i_vstore)=current
  end do

! _______________________________________________________________
! Here Pulay algorithm
 else if(iscf==7)then

  niter=min(istep,npulay+1)

  write(message,'(2a,i2,a)') ch10,'Pulay update with ',niter-1,' previous iterations:'
  call wrtout(6,message,'COLL')

  if (npulay>npulaymax) then
   write(message, '(4a)' ) ch10,&
&   ' scfopt : ERROR - ',ch10,&
&   '  Too much iterations required for Pulay algorithm (<50) !'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

! Compute "A" matrix
  if (istep>npulay+1) then
   do jj=1,niter-1
    do ii=1,niter-1
     amat(ii,jj)=amat(ii+1,jj+1)
    end do
   end do
  end if
  do ii=1,niter
   call dotprodm_v(cplex,1,amat(ii,niter),i_vrespc(1),i_vrespc(1+niter-ii),mpi_enreg,1,1,&
&   nfft,n_fftgr,n_fftgr,nspden,f_fftgr,f_fftgr)
   if (usepaw==1.and.pawoptmix==1) then
    do index=1,npawmix
     amat(ii,niter)=amat(ii,niter)+f_paw(index,i_vrespc(1))*f_paw(index,i_vrespc(1+niter-ii))
    end do
   end if
   if (ii<niter) amat(niter,ii)=amat(ii,niter)
  end do

! Invert "A" matrix
  allocate(amatinv(niter,niter))
  amatinv(1:niter,1:niter)=amat(1:niter,1:niter)
  allocate(ipiv(niter),rwork(niter))
  call dgetrf(niter,niter,amatinv,niter,ipiv,ierr)
  call dgetri(niter,amatinv,niter,ipiv,rwork,niter,ierr)
  deallocate(ipiv,rwork)

! Compute "alpha" factors
  allocate(alpha(niter));alpha=zero;det=zero
  do ii=1,niter
   do jj=1,niter
    alpha(ii)=alpha(ii)+amatinv(jj,ii)
    det=det+amatinv(jj,ii)
   end do
  end do
  alpha(:)=alpha(:)/det
  deallocate(amatinv)
  write(message,'(a,5(1x,g10.3))')' mixing of old trial potential : alpha(m:m-4)=',(alpha(ii),ii=niter,max(1,niter-4),-1)
  call wrtout(6,message,'COLL')

! Save latest trial potential and compute new trial potential
  do isp=1,nspden
   do ifft=1,cplex*nfft
    current=vtrial(ifft,isp)
    vtrial(ifft,isp)=alpha(niter)*(current+f_fftgr(ifft,isp,i_vrespc(1)))
    do ii=niter-1,1,-1
     vtrial(ifft,isp)=vtrial(ifft,isp)+alpha(ii) &
&     *(f_fftgr(ifft,isp,i_vtrial(niter-ii))+f_fftgr(ifft,isp,i_vrespc(1+niter-ii)))
    end do
    f_fftgr(ifft,isp,i_vstore)=current
   end do
  end do

! PAW: save latest rhoij and compute new rhoij
  do index=1,npawmix
   current=vpaw(index)
   vpaw(index)=alpha(niter)*(current+f_paw(index,i_vrespc(1)))
   do ii=niter-1,1,-1
    vpaw(index)=vpaw(index)+alpha(ii) &
&    *(f_paw(index,i_vtrial(niter-ii))+f_paw(index,i_vrespc(1+niter-ii)))
   end do
   f_paw(index,i_vstore)=current
  end do

  deallocate(alpha)

! _______________________________________________________________
! End of choice of optimization method
 end if

!Permute potential indices
 if (iscf==3) then
  tmp=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
 else if (iscf==4) then
  tmp=i_vrespc(3) ; i_vrespc(3)=i_vrespc(2) ; i_vrespc(2)=i_vrespc(1) ; i_vrespc(1)=tmp
  tmp=i_vtrial(2) ; i_vtrial(2)=i_vtrial(1) ; i_vtrial(1)=tmp
 else if (iscf==7) then
  tmp=i_vtrial(  npulay)
  do ii=  npulay,2,-1
   i_vtrial(ii)=i_vtrial(ii-1)
  end do
  i_vtrial(1)=tmp
  tmp=i_vrespc(1+npulay)
  do ii=1+npulay,2,-1
   i_vrespc(ii)=i_vrespc(ii-1)
  end do
  i_vrespc(1)=tmp
 end if

end subroutine scfopt
!!***
