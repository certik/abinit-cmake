!{\src2tex{textfont=tt}}
!!****f* ABINIT/crho
!! NAME
!! crho
!!
!! FUNCTION
!! Calculate the charge density rho on the FFT grid.
!! In case of nsppol==2 calculate rho_up and rho_down 
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  irottb(nfftot,nsym) = symmetry operations on the FFT grid.
!!  irottb(r,R)=index of (R**-1)r  in the FFT array where R is one of the nsym
!!   symmetry operation in reciprocal space
!!  mpi_enreg= datatype gathering information on parallelism, variables used 
!!   |gwpara= if 2 bands are spread btw processors
!!  nbnds = number of bands.
!!  timrev=2 if time-reversal symmetry is used, 1 otherwise.
!!  nkbz = number of k-points in the full Brillouin zone.
!!  nkibz = number of k-points in the irreducible Brillouin zone.
!!  nsym = number of symmetry operations.
!!  nfftot = total number of points on the FFT grid.
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ(nkibz,nbnds,nsppol) = occupation numbers for nbnds bands at nkibz irreducible k-points, for each spin
!!  ucvol = unit cell volume.
!!  wfr(nfftot,nbnds,nkibz,nsppol) = wavefunctions on the FFT grid for nbnds bands at nkibz irreducible k-points, for each spin
!!  wtk(nkibz) = irreducible k-points weights.
!!
!! OUTPUT
!!  omegaplasma = the plasma frequency.
!!  rho(nfftot,nsppol) = the density on the FFT grid.
!!   (total in first half and spin-up in second half if nsppol=2)
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine crho(paral_kgb,ngfft,irottb,nbnds,timrev,nkbz,nkibz,nsym,symrel,tnons,symafm,&
& nfftot,nspden,nsppol,occ,omegaplasma,rho,ucvol,wfr,wtk,mpi_enreg,my_minb,my_maxb)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : get_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13recipspace
 use interfaces_15common
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_maxb,my_minb,nbnds,nfftot,nkbz,nkibz,nspden,nsppol
 integer,intent(in) :: nsym,paral_kgb,timrev
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: omegaplasma
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: irottb(nfftot,nsym),ngfft(18),symafm(nsym)
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: occ(nkibz,nbnds,nsppol),tnons(3,nsym),wtk(nkibz)
 real(dp),intent(out) :: rho(nfftot,nsppol)
 complex,intent(in) :: wfr(nfftot,my_minb:my_maxb,nkibz,nsppol)

!Local variables ------------------------------
!scalars
 integer :: cplex,ib,ier,ierr,ik,iop,ir,is,master,me,n1,n2,n3,spaceComm,unt
 real(dp) :: fact,rhoav,rs,tnepuc
 logical :: DEBUG
 character(len=500) :: message
 character(len=fnlen) :: filnam
 type(Dens_sym_operator_type) :: densymop
 type(MPI_type) :: MPI_enreg_seq
!arrays
 integer,allocatable :: irrzon(:,:,:)
 real(dp),allocatable :: phnons(:,:,:),rho2(:),rhog(:,:)

!*************************************************************************

#if defined DEBUG_MODE
 write(message,'(a)')' crho : enter '
 call wrtout(std_out,message,'PERS')
#endif

 write(message,'(2a)')ch10,' crho: calculating charge density...'
 call wrtout(std_out,message,'COLL')
!
!Initialize some MPI related variables
 call xcomm_init(mpi_enreg,spaceComm)  
 call xme_init(mpi_enreg,me)          
 call xmaster_init(mpi_enreg,master)
!
!Calculate IBZ contribution to the charge density
 allocate(rho2(nfftot))
 rho(:,:)= zero
 do is=1,nsppol
  rho2(:)= zero
! Loop over k-points in IBZ
  do ik=1,nkibz
!  Skip the higher bands if occupation is less than tol8
!  do while ((abs(occ(ik,ib,is))>tol8).and.(ib<=nbnds))
   do ib=1,nbnds
    if (mpi_enreg%gwpara==2) then
     if (mpi_enreg%proc_distrb(ik,ib,is)/=me) cycle
    end if
    if (abs(occ(ik,ib,is))<tol8) cycle 
    do ir=1,nfftot
     rho2(ir)= rho2(ir) + occ(ik,ib,is)*conjg(wfr(ir,ib,ik,is))*wfr(ir,ib,ik,is)*wtk(ik)/SUM(wtk)/ucvol
    end do !ir
   end do !ib
  end do !ik
! we could sum rho outside the loop over is 
  if (mpi_enreg%gwpara==2) call xsum_mpi(rho2,spaceComm,ier)
! 
! Loop over symmetry operations, symmetrising rho.
! Factor 2 is for inversion
! fact=real(timrev)/(nkbz*ucvol)
! do iop=1,nsym
! do ir=1,nfftot
! rho(ir,is)= rho(ir,is)+fact*rho2(irottb(ir,iop))
! end do
! end do
  rho(:,is) = rho2(:)

 end do !is
!
!Store the total charge in the first half
!if (nsppol==2) then 
!rho2(:) = rho(:,1)
!rho(:,1)= rho(:,1)+rho(:,2)
!rho(:,2)= rho2(:)
!end if 

!NEW symmetrization in G space implementing also the AFM case.
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3) 
 allocate(irrzon(nfftot,2,nspden/nsppol))
 allocate(phnons(2,nfftot,nspden/nsppol))

 call irrzg(densymop,irrzon,nspden,nsppol,nsym,n1,n2,n3,&
& phnons,symafm,symrel,tnons)

!* Fake MPI_type for sequential part
 call initmpi_seq(MPI_enreg_seq) ; MPI_enreg_seq%nproc_fft=1 ; MPI_enreg_seq%me_fft=0

 cplex=1 
 allocate(rhog(2,cplex*nfftot)) !this might be output

 call symrhg(cplex,densymop,irrzon,MPI_enreg_seq,nfftot,nfftot,ngfft,nspden,nsppol,&
& nsym,paral_kgb,phnons,rhog,rho,symafm)

 deallocate(rhog,phnons,irrzon)

!Write total charge
 DEBUG=.false.
 if (DEBUG .and. me==0) then 
  filnam='crho.dat'
  call isfile(filnam,'new')
  unt=get_unit()
  open(unit=unt,file=filnam)
  do ir=1,nfftot
   write(unt,'(2x,f22.15)') (rho(ir,is),is=1,nsppol)
  end do 
  close(unt)
 end if 
!
!Calculate total number of electrons as a check
 tnepuc=zero
 do ir=1,nfftot
  tnepuc=tnepuc+rho(ir,1)
 end do
 tnepuc=tnepuc*ucvol/nfftot ; rhoav=tnepuc/ucvol ; rs=(three/(four_pi*rhoav))**third

 write(message,'(a,f9.4)')' total number of electrons per unit cell = ',tnepuc
 call wrtout(std_out,message,'COLL') 
 call wrtout(ab_out,message,'COLL')
 write(message,'(a,f9.6)')' average of density, n = ',rhoav
 call wrtout(std_out,message,'COLL') 
 call wrtout(ab_out,message,'COLL')
 write(message,'(a,f9.4)')' r_s = ',rs
 call wrtout(std_out,message,'COLL') 
 call wrtout(ab_out,message,'COLL')
 omegaplasma= sqrt(four_pi*rhoav)
 write(message,'(a,f9.4,2a)')' omega_plasma = ',omegaplasma*Ha_eV,' [eV]',ch10
 call wrtout(std_out,message,'COLL') 
 call wrtout(ab_out,message,'COLL')

 deallocate(rho2)

#if defined DEBUG_MODE
 write(message,'(a)')' crho : exit '
 call wrtout(std_out,message,'PERS')
#endif

end subroutine crho
!!***
