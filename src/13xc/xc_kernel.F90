!{\src2tex{textfont=tt}}
!!****f* ABINIT/xc_kernel
!! NAME
!! xc_kernel
!!
!! FUNCTION
!! Calculate exchange-correlation kernel in reciprocal space
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (Rhaltaf,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!! ixc = choice for the exchange-correlation potential.
!! mpi_enreg = informations about MPI parallelization.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nr = total number of points on the FFT grid.
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! rho(nr,nsppol) = the charge density on the FFT grid.
!!  (total in first half and spin-up in second half if nsppol=2)
!! rprimd(3,3) = dimensional real space primitive translations (bohr).
!! npw: the size of kernel matrix
!! igfft=array of fft index of each G vector 
!!
!! OUTPUT
!!  kxc_kernel(ng,ngp,nq) = the exchange-correlation potential on the FFT grid.
!!  warning: the kernel is not devided by unit cell volume
!!
!! NOTES
!!  
!! No xc quadrature
!! No nl core correction
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      leave_new,rhohxc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine xc_kernel(dtset,ixc,mpi_enreg,ngfft,nr,nsppol,rho,rprimd,igfft,npw,gmet,kernel,gvec,nq,qq)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_12ffts
 use interfaces_13xc, except_this_one => xc_kernel
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,npw,nq,nr,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: gvec(3,npw),igfft(npw),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),qq(3,nq),rho(nr,nsppol),rprimd(3,3)
 complex,intent(out) :: kernel(npw,npw,nq)

!Local variables ------------------------------
!This section has been created automatically by the script Abilint (TD). Do not modify these by hand.
!End of the abilint section
!scalars
 integer :: cplex,i1,i2,i3,ig,igp,iq,ir,ir1,ishift,ispden,n3xccc,ngfft1,ngfft2
 integer :: ngfft3,ngrad,nkxc,nspden,nspgrad,option
 real(dp) :: enxc,expo,fac,gpqx,gpqy,gpqz,gsqcut,gx,gxrx,gy,gyry,gz,gzrz
 real(dp) :: imag_part,qx,qy,qz,real_part,rx,ry,rz,sum_i,sum_r,summ_i,summ_r
 real(dp) :: vxcavg
 complex :: summ
 logical :: gga,lda
 character(len=500) :: message
 type(dataset_type) :: dtGW
!arrays
 real(dp) :: gprimd(3,3),qphon(3),strsxc(6)
 real(dp),allocatable :: dum(:),kernel_dp(:,:),kxc(:,:),phas(:,:,:),phasp(:,:)
 real(dp),allocatable :: rhog(:,:),rhos(:,:),vhartr(:),vxc1(:,:),vxclda(:,:)
 real(dp),allocatable :: xccc3d(:),xx(:,:)

!************************************************************************
!print info

 call dtsetCopy(dtGW, dtset)
 dtGW%intxc = 0

 write(message,'(a,i3)') ' xc_kernel: calculating exchange-correlation kernel using ixc = ',ixc
 call wrtout(std_out,message,'COLL')
 if(nsppol>1)then
  write(message,'(a,a)') ' xc_kernel: spin nonpolarizibility will be enforced = ',ch10
  call wrtout(std_out,message,'COLL')
 end if
 
!form Vxc (in Hartrees)

!note: one must have nr=ngfft1*ngfft2*ngfft3 (ie the FFT grid must not be augmented a priori)
!this is actually enforced at the present time in setmesh.f
 ngfft1=ngfft(1) ; ngfft2=ngfft(2) ; ngfft3=ngfft(3)

 allocate(rhog(2,nr),vhartr(nr))


 nspden=1
 option=2
 qphon(:)=0.0

 if(ixc>=1.and.ixc<11)then
  lda=.true.
  gga=.false.
 else
  lda=.false.
  gga=.true.
 end if

 if(lda)then
  nkxc=3
 else
  nkxc=23
 end if
 
 allocate(kxc(nr,nkxc))


!gsqcut & rhog are zeroed because they are not used by rhohxc if 1<=ixc<=16 and option=0
 gsqcut=0.0
 rhog(:,:)=0.0
!MG FIXME this is the 3D core electron density for XC core correction (bohr^-3)
!should implement the non linear core correction 
 n3xccc=0       
 allocate(xccc3d(n3xccc),vxclda(nr,nspden),rhos(nr,nspden))
 

!MG Now rho(nr,sp%nsppol)

 call rhohxc(dtGW,enxc,gsqcut,0,kxc,mpi_enreg,nr,ngfft,&
& dum,0,dum,0,nkxc,nspden,n3xccc,option,rhog,rho,rprimd,strsxc,1,&
& vhartr,vxclda,vxcavg,xccc3d)


 deallocate(xccc3d,vxclda)



 cplex=2
 allocate(phas(cplex*nr,npw,nspden),vxc1(cplex*nr,nspden),xx(3,nr),kernel_dp(2,nr))
 
 kernel(:,:,:)=czero

!find the corrdinates for all r in the FFT grid

 ir=0
 do i3=1,ngfft3
  do i2=1,ngfft2
   do i1=1,ngfft1
    ir=ir+1
    xx(1,ir)=dble((i1-1))/ngfft1
    xx(2,ir)=dble((i2-1))/ngfft2
    xx(3,ir)=dble((i3-1))/ngfft3
   end do
  end do
 end do

!calculate at once exp(i(G+q).r), for all possible q,G,r

 do iq=1,1
  do ig=1,npw
   gpqx=dble(gvec(1,ig))
   gpqy=dble(gvec(2,ig))
   gpqz=dble(gvec(3,ig))
   do ir=1,nr
    expo=gpqx*xx(1,ir)+gpqy*xx(2,ir)+gpqz*xx(3,ir)              
    phas(2*ir-1,ig,1)= cos(2*pi*expo)
    phas(2*ir,ig,1) =  sin(2*pi*expo)
   end do
  end do
! start the calculation of $K(G,G'',q)=\frac{1}{nr}\sum_{r} exp(-i(q+G_{2}).r_{2} kxc(r_{1}r_{2}) exp(i(q+G_{1}).r_{1} $

  do igp=1,npw

   vxc1(:,:)=zero



   call mkvxc3(cplex,gmet,gsqcut,kxc,mpi_enreg,nr,ngfft,nkxc,nspden,n3xccc,option,&
&   dtset%paral_kgb, qphon(:),phas(:,igp,:),rprimd,vxc1,xccc3d)



   mpi_enreg%me_fft=0
   mpi_enreg%nproc_fft=1
   
   summ_r=zero
   summ_i=zero


!  FFT the first index to --> to G space

   call fourdp(cplex,kernel_dp(:,:),vxc1(:,1),-1,mpi_enreg,nr,ngfft,dtset%paral_kgb,0)


   kernel(:,igp,iq)=cmplx(kernel_dp(1,igfft(:)),kernel_dp(2,igfft(:)))


  end do ! igp


 end do

!must devide kernel by ucvol later as required by fft r--> G
!devision is done in eps1_tc.F0 

 

 do iq=1,nq
  kernel(:,:,iq)= kernel(:,:,1)
 end do
 
 deallocate(kxc,phas,vxc1)

 call dtsetFree(dtGW)

 write(std_out,*)


end subroutine xc_kernel

!!***
