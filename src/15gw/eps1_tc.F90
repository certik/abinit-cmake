!{\src2tex{textfont=tt}}
!!****f* ABINIT/eps1_tc
!! NAME
!! eps1_tc
!!
!! FUNCTION
!! Starting from the RPA DM calculate TC DM
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
!! mpi_enreg = informations about MPI parallelization.
!! ngfft(18)=contain all needed information about 3D FFT,
!!  see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nr = total number of points on the FFT grid.
!! rho(nr,nsppol) = the charge density on the FFT grid.
!!  (total in first half and spin-up in second half if nsppol=2)
!! rprimd(3,3) = dimensional real space primitive translations (bohr).
!! npweps: the size of kernel matrix
!! igfft=array of fft index of each G vector 
!! epsm1: RPA eps1
!! 
!! OUTPUT
!!  epsm1: TC eps1
!!  binary file with extension GWG_SCR that includes $\varepsilon^{-1}$ with vertex corrections
!!
!! NOTES
!!  
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

subroutine eps1_tc(dtset,mpi_enreg,ngfft,nr,rho,rprimd,igfft,&
& npweps,gmet,gprimd,gvec,nomega,epsm1,nq,ucvol,qq,omega,dtfil,hdr,npwwfn,npwvec,nbnds)

 use defs_basis
 use defs_datatypes
 use m_numeric_tools, only : hermitianize,print_arr


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13xc
 use interfaces_15gw, except_this_one => eps1_tc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nr,npweps,nomega,nq,npwwfn,npwvec,nbnds
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(hdr_type),intent(inout) :: hdr
!arrays
 integer,intent(in) :: igfft(npweps),gvec(3,npweps)
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rho(nr,dtset%nsppol),rprimd(3,3),gmet(3,3),gprimd(3,3),ucvol,qq(3,nq) 
 complex(gwpc),intent(in):: omega(nomega)
 complex(gwpc),intent(inout):: epsm1(npweps,npweps,nomega,nq)
!Local variables ------------------------------
!scalars
 integer :: igp,ig,iq,ixc,nsppol,iomega,ig1,ig2,ig4,ig4x,ig4y,ig4z
 integer :: optfil,nqsmall
 character(len=500) :: message
!arrays
 real(dp) :: b1(3),b2(3),b3(3)
 real(dp) :: epsilon0
 real(dp),allocatable::qplusg(:)
 real(dp),allocatable :: qsmall_dum(:,:)
 type(datafiles_type) :: dtfil_tmp
 complex(gwpc),allocatable :: kernel(:,:,:),chi0(:,:,:),chitmp(:,:),chi(:,:,:),kxclda(:)
 complex(gwpc),allocatable :: lwing_dum(:,:,:),uwing_dum(:,:,:)
 character(len=80) :: title(2)

!************************************************************************

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2) 
 b3=two_pi*gprimd(:,3)

!calculate the exc kerenel
 ixc=dtset%ixc
 nsppol=dtset%nsppol

!MG in case of contour deformation these arrays might quite large
!meanwhile in sigma we have allocated er%epsm1 so 
!there is a possible out of of memory. I wrote two subroutines 
!to calculate e^-1 starting from X0  and including lda vertex correction 
!In these routines all the operations are done in place, if required, 
!we should use these. 
!
 allocate(kernel(npweps,npweps,nq),chi0(npweps,npweps,nomega))
 allocate(chitmp(npweps,npweps),chi(npweps,npweps,nomega),qplusg(npweps))

 call xc_kernel(dtset,ixc,mpi_enreg,ngfft,nr,nsppol,rho,rprimd,igfft,npweps,gmet,kernel,gvec,nq,qq)


 write(message,'(a)')' kxclda(g)'
 call wrtout(06,message,'COLL')
 call print_arr(kernel(1,:,1))


 title(1)='em1 file: epsilon^-1'
 title(2)='TESTELECTRON'

!extract chi0 from RPA eps1
 
 do iq=1,nq
  do iomega=1,nomega
   call matcginv(epsm1(:,:,iomega,iq),npweps,npweps)
  end do
 end do


 do iq=1,nq   
  call cvc(nq,iq,qq,npweps,gvec,gprimd,qplusg)

  do iomega=1,nomega
   do ig=1,npweps
    epsm1(ig,ig,iomega,iq)=epsm1(ig,ig,iomega,iq)-cone
    do igp=1,npweps
     chi0(ig,igp,iomega)=(qplusg(ig)*qplusg(igp)*epsm1(ig,igp,iomega,iq))/(-4.0*pi)
    end do
   end do
  end do


  chi0(:,:,:)=chi0(:,:,:)*ucvol
  
  do iomega=1,nomega
   write (message,'(a,i2,a,i1,a)')' chi0(q=',iq,',omega=',iomega,',G,G")'
   call wrtout(06,message,'COLL')
   call print_arr(chi0(:,:,iomega))
  end do


! starting from chi0 calculate chi=(1-chi0*Vc-chi0*Kernel)^-1 chi0'


  do iomega=1,nomega
   chitmp=matmul(chi0(:,:,iomega),kernel(:,:,iq))

   do ig=1,npweps
    do igp=1,npweps
     chitmp(ig,igp)=(-(4.0*pi)/(qplusg(igp)**2))*&
&     chi0(ig,igp,iomega)-chitmp(ig,igp)
    end do
    chitmp(ig,ig)=chitmp(ig,ig)+ucvol
   end do
   
!  invert (1-chi0*Vc-chi0*Kernel)
   call matcginv(chitmp(:,:),npweps,npweps)
   
!  multiply for chi0 to obtain chi
   chi(:,:,iomega)=matmul(chitmp(:,:),chi0(:,:,iomega))
   
!  must multiply by ucvol (see notes-due to conversion of r->G in
!  matrix definitions)
   chi(:,:,iomega)=chi(:,:,iomega)*ucvol
   
  end do 

  do iomega=1,nomega
   write (message,'(a,i2,a,i1,a)')' chi(q=',iq,',omega=',iomega,',G,G")'
   call wrtout(06,message,'COLL')
   call print_arr(chi(:,:,iomega))
  end do


! starting from chi calculate Tc-eps1

  chitmp(:,:)=czero
  epsm1(:,:,:,iq)=czero

  do iomega=1,nomega

   chitmp=matmul(kernel(:,:,iq),chi(:,:,iomega))

!  Perform hermitianization
   call hermitianize(chitmp)

   do ig=1,npweps
    do igp=1,npweps
     epsm1(ig,igp,iomega,iq)=(4*pi/(qplusg(ig)*qplusg(igp)))*&
&     chi(ig,igp,iomega)/ucvol+chitmp(ig,igp)/ucvol
!    chitmp=Kxc*chi must be divided by ucvol see notes on fft r->G
    end do
    epsm1(ig,ig,iomega,iq)=epsm1(ig,ig,iomega,iq)+cone
   end do !ig
   
  end do ! iomega

! Print epsilon-twiddle^-1
  do iomega=1,nomega
   write (message,'(a,i2,a,i1,a)')' epsilon-twiddle^-1(q=',iq,',omega=',iomega,',G,G")'
   call wrtout(06,message,'COLL')
   call print_arr(epsm1(:,:,iomega,iq))
  end do


  if (all(abs(qq(:,iq))<1.0e-3))then
   do iomega=1,nomega
    if(abs(real(omega(iomega)))<1.e-3) then
     if(abs(aimag(omega(iomega)))<1.e-3) then
      epsilon0=real(1.0/epsm1(1,1,iomega,iq))
      write(message,'(a,a,f8.4,a)')' dielectric constant after including',& 
&     ' vertex corrections = ',epsilon0,ch10
      call wrtout(06,message,'COLL') 
     end if
    end if
   end do
  end if

! MG Thu Oct 18 
! Here I changed wrscr, now the optional argument unit has been removed
  dtfil_tmp%filnam_ds(4)=dtfil%filnam_ds(4)
  dtfil_tmp%unscr=543

  optfil=1 ; nqsmall=0
  allocate(lwing_dum(npweps,nomega,optfil),uwing_dum(npweps,nomega,optfil),qsmall_dum(3,optfil))

  call wrscr(iq,optfil,dtfil_tmp%unscr,dtfil_tmp%filnam_ds(4),hdr,dtset,npweps,npwwfn,nbnds,&
&  nq,nqsmall,nomega,qq,omega,gvec,gmet,epsm1,title,&
&  qsmall_dum,uwing_dum,lwing_dum)

  deallocate(lwing_dum,uwing_dum,qsmall_dum)

 end do ! iq

 deallocate(kernel,chi0,chitmp,chi,qplusg)


end subroutine eps1_tc
!!***
