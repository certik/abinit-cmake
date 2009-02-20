!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkvxcstr3
!! NAME
!! mkvxcstr3
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! due to strain: assemble the first-order density change with the
!! frozen-core density change, then use
!! the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2008 ABINIT group (DRH,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!     if 2, COMPLEX
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  idir=direction of the current perturbation
!!  ipert=type of the perturbation
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=informations about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used
!!   otherwise, cplex*nfft
!!  option=if 0, work only with strain-derivative frozen-wavefunction
!!    charge and the XC core-correction,
!!   if 1, treat both density change and XC core correction
!!   if 2, like 0 but multiply gradient strain derivative term by 2.0 for GGA.
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor(nfft,nspden)=array for GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential (including
!!   core-correction, if applicable)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      eltfrxc3,eneres3,nselt3,scfcv3
!!
!! CHILDREN
!!      leave_new,matr3inv,mkvxcstrgga3,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mkvxcstr3(cplex,gmet,gsqcut,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,&
& nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,vxc1,xccc3d1)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_16response, except_this_one => mkvxcstr3
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,n3xccc,natom,nfft,nkxc,nspden,option
 integer,intent(in) :: paral_kgb
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gmet(3,3),kxc(nfft,nkxc),qphon(3),rhor(nfft,nspden)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rprimd(3,3)
 real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,save :: npass=0
 integer :: ii,ir,istr,jj,ka,kb
 real(dp) :: rho1_dn,rho1_up,rho1im_dn,rho1im_up,rho1re_dn,rho1re_up
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 real(dp) :: dgprimdds(3,3),gprimd(3,3),tsec(2)
 real(dp),allocatable :: rhor1tmp(:,:),rhowk1(:,:)

! *************************************************************************

!DEBUG
!write(6,*)' mkvxc3 : enter '
!if(option==1)stop
!ENDDEBUG

 call timab(181,1,tsec)

 if(nspden/=1 .and. nspden/=2) then
  write(message, '(a,a,a,a)' ) ch10,&
&  ' mkvxc3 : BUG -',ch10,&
&  '  Only for nspden==1 and 2.'
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if

!Inhomogeneous term for diagonal strain
 allocate(rhowk1(nfft,nspden))
 if(option==0 .or. option==2) then
  if(ipert==natom+3) then
   rhowk1(:,:)=-rhor(:,:)
  else
   rhowk1(:,:)=zero
  end if
 else if(option==1) then
  if(ipert==natom+3) then
   rhowk1(:,:)=rhor1(:,:)-rhor(:,:)
  else
   rhowk1(:,:)=rhor1(:,:)
  end if
 end if

!Treat first LDA
 if(nkxc/=23)then

! Case without non-linear core correction
  if(n3xccc==0)then

!  Non-spin-polarized
   if(nspden==1)then
    do ir=1,nfft
     vxc1(ir,1)=kxc(ir,1)*rhowk1(ir,1)
    end do

!   Spin-polarized
   else
    do ir=1,nfft
     rho1_dn=rhowk1(ir,1)-rhowk1(ir,2)
     vxc1(ir,1)=kxc(ir,1)*rhowk1(ir,2)+kxc(ir,2)*rho1_dn
     vxc1(ir,2)=kxc(ir,2)*rhowk1(ir,2)+kxc(ir,3)*rho1_dn
    end do
   end if ! nspden==1

!  Treat case with non-linear core correction
  else
   if(nspden==1)then
    do ir=1,nfft
     vxc1(ir,1)=kxc(ir,1)*(rhowk1(ir,1)+xccc3d1(ir))
    end do
   else
    do ir=1,nfft
     rho1_dn=rhowk1(ir,1)-rhowk1(ir,2) + xccc3d1(ir)*half
     rho1_up=rhowk1(ir,2)             + xccc3d1(ir)*half
     vxc1(ir,1)=kxc(ir,1)*rho1_up+kxc(ir,2)*rho1_dn
     vxc1(ir,2)=kxc(ir,2)*rho1_up+kxc(ir,3)*rho1_dn
    end do
   end if ! nspden==1

  end if ! n3xccc==0

! Treat GGA
 else

  allocate(rhor1tmp(cplex*nfft,2))

! Generates gprimd and its strain derivative
! Note that unlike the implicitly symmetric metric tensor strain
! derivatives, we must explicltly symmetrize the strain derivative
! here.

  call matr3inv(rprimd,gprimd)

  istr=idir + 3*(ipert-natom-3)

  if(istr<1 .or. istr>6)then
   write(message, '(a,a,a,a,i10,a,a,a)' )ch10,&
&   ' mkvxcstr3: BUG -',ch10,&
&   '  Input dir gives istr=',istr,' not allowed.',ch10,&
&   '  Possible values are 1,2,3,4,5,6 only.'
   call wrtout(06,message,'PERS')
   call leave_new('PERS')
  end if

  dgprimdds(:,:)=zero
  ka=idx(2*istr-1);kb=idx(2*istr)
  do ii=1,3
   do jj=1,3
    if(jj==ka) dgprimdds(jj,ii)=-half*gprimd(kb,ii)
    if(jj==kb) dgprimdds(jj,ii)=dgprimdds(jj,ii)-half*gprimd(ka,ii)
   end do
  end do

! First transfer the data to spin-polarized storage
  if(nspden==1)then
   do ir=1,cplex*nfft
    rhor1tmp(ir,1)=rhowk1(ir,1)*half
    rhor1tmp(ir,2)=rhowk1(ir,1)*half
   end do
  else
   do ir=1,cplex*nfft
    rho1_dn=rhowk1(ir,1)-rhowk1(ir,2)
    rhor1tmp(ir,1)=rhowk1(ir,2)
    rhor1tmp(ir,2)=rho1_dn
   end do
  end if ! nspden==1
  if(n3xccc/=0)then
   do ir=1,cplex*nfft
    rhor1tmp(ir,1)=rhor1tmp(ir,1)+xccc3d1(ir)*half
    rhor1tmp(ir,2)=rhor1tmp(ir,2)+xccc3d1(ir)*half
   end do
  end if

! Rescalling needed for use in eltfrxc3 for elastic tensor (not internal
! strain tensor).
  if(option==2) dgprimdds(:,:)=2.0_dp*dgprimdds(:,:)
  call mkvxcstrgga3(cplex,dgprimdds,gmet,gprimd,gsqcut,kxc,mpi_enreg,nfft,ngfft,nkxc,&
&  nspden,paral_kgb,qphon,rhor1tmp,vxc1)
  deallocate(rhor1tmp)

 end if

 deallocate(rhowk1)

 call timab(181,2,tsec)

end subroutine mkvxcstr3
!!***
