!{\src2tex{textfont=tt}}
!!****f* ABINIT/cvxclda
!! NAME
!! cvxclda
!!
!! FUNCTION
!! Calculate Vxc on the FFT grid.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, YMN, MG)
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
!! nfftot = total number of points on the FFT grid.
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! rho(nfftot,nsppol) = the charge density on the FFT grid.
!!  (total in first half and spin-up in second half if nsppol=2)
!! rprimd(3,3) = dimensional real space primitive translations (bohr).
!!
!! OUTPUT
!!  vxclda(nfftot,nsppol) = the exchange-correlation potential on the FFT grid.
!!                      (spin up in first half and spin down in second half if nsppol=2)
!!
!! NOTES
!!  
!! No xc quadrature
!! No nl core correction
!!
!! PARENTS
!!      cmevxclda
!!
!! CHILDREN
!!      leave_new,rhohxc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine cvxclda(dtset,ixc,mpi_enreg,ngfft,nfftot,nsppol,rho,rprimd,vxclda)

 use defs_basis
 use defs_datatypes
 use m_IO_tools, only : get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13xc
 use interfaces_15gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,nfftot,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rho(nfftot,nsppol),rprimd(3,3)
 real(dp),intent(out) :: vxclda(nfftot,nsppol)

!Local variables ------------------------------
!scalars
 integer,parameter :: nkxc=0,option=0
 integer :: ir,is,n3xccc,nspden,unt
 real(dp) :: enxc,gsqcut,vxcavg
 logical :: paranoia
 character(len=500) :: message
 character(len=fnlen) :: fnam
 type(dataset_type) :: dtGW
!arrays
 real(dp) :: strsxc(6)
 real(dp),allocatable :: dum1(:,:),dum2(:,:,:),kxc(:,:),rhog(:,:),vhartr(:)
 real(dp),allocatable :: xccc3d(:)

!************************************************************************

 write(message,'(a,i3)')' cvxclda: calculating Vxc using ixc = ',ixc
 call wrtout(std_out,message,'COLL')
!
!Form Vxc (in Hartree units)
!
 if (ixc==0) then
! For backward compatibility
  write (message,'(6a)')ch10,&
&  ' cvxclda: WARNING - ',ch10,&
&  ' ixc = 0 is a relativistic Ceperley-Alder xc functional [PRB 26, p. 4199, (1982)]',ch10,&
&  ' in the GW code. It should be used only for backward compatibility.'
  call wrtout(ab_out,message,'COLL') 
  call wrtout(std_out,message,'COLL')
  if (nsppol==2) then 
   write (message,'(3a)')&
&   ' cvxclda: ERROR- ',ch10,&
&   ' ixc = 0 and nsppol==2 not yet implemented ' 
   call wrtout(std_out,message,'COLL') 
   call leave_new('COLL')
  end if 
  do ir=1,nfftot
   vxclda(ir,1)=vxcca(rho(ir,1))
  end do
 else
! MG this comes from my private branch (5.3.4)
! if ((ixc<0).or.(ixc>16)) then 
! this comes from XG merge-public--5.4.3--patch-2
  if (ixc>=10 .and. dtset%xclevel/=2) then
   write (message,'(4a,i3,a)') ch10,&
&   ' cvxclda: ERROR - ',ch10,&
&   '  ixc = ',ixc,' is not allowed at the present time in the GW code.'
   call wrtout(ab_out,message,'COLL') 
   call leave_new('COLL')
  end if
! 
! Copy the input variables from the current dataset to a temporary one to tune some parameters
  call dtsetCopy(dtGW,dtset)
  dtGW%intxc=0
  write(message,'(a)')&
&  ' cvxclda: calling rhohxc to calculate Vxc[n_val] (excluding non-linear core corrections)'
  call wrtout(std_out,message,'COLL')
! 
! Note: one must have nfftot=ngfft1*ngfft2*ngfft3, ie the FFT grid must not 
! be augmented. This is actually enforced at the present time in setmesh.f
! 
  allocate(rhog(2,nfftot),vhartr(nfftot))
  allocate(kxc(nfftot,nkxc))
! gsqcut and rhog are zeroed because they are not used by rhohxc if 1<=ixc<=16 and option=0
  gsqcut=zero ; rhog(:,:)=zero
! 
! TODO this is the 3D core electron density for XC core correction (bohr^-3)
! should implement the non linear core correction 
  n3xccc=0       
  allocate(xccc3d(n3xccc))
! 
! nkxc=0  ==> no computation of the exchange-correlation kernel
! option=0  ==> only exc, vxc, strsxc
! 
  nspden=nsppol
  allocate(dum1(nfftot,0),dum2(nfftot,nspden,0))
  call rhohxc(dtGW,enxc,gsqcut,0,kxc,mpi_enreg,nfftot,ngfft,dum1,0,dum2,0,nkxc,nspden,&
&  n3xccc,option,rhog,rho,rprimd,strsxc,1,vhartr,vxclda,vxcavg,xccc3d)
  deallocate(dum1,dum2)

  deallocate(rhog,vhartr,xccc3d)
  deallocate(kxc)
  call dtsetFree(dtGW)

  write(message,'(a,f8.4,2a,f8.4,2a)')&
&  ' cvxclda: rhohxc returned  Exc[n_val]  = ',enxc,  ' [Ha]',&
&  ' and <Vxc[n_val]> = ',vxcavg,' [Ha]',ch10
  call wrtout(std_out,message,'COLL')
 end if

 paranoia=.false.
 if (paranoia) then
! Write Vxc
  fnam='__vxc.dat__' 
  call isfile(fnam,'new')
  unt=get_unit() ; open(unit=unt,file=fnam)
  write(unt,'(es16.8)') ((vxclda(ir,is),ir=1,nfftot),is=1,nsppol)
  close(unt)
 end if

#if defined DEBUG_MODE
 write(message,'(a)')'cvxclda : exit'
 call wrtout(std_out,message,'PERS')
!call flush_unit(std_out)
#endif

end subroutine cvxclda

!!***
