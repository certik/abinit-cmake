!{\src2tex{textfont=tt}}
!!****f* ABINIT/symrhg
!! NAME
!! symrhg
!!
!! FUNCTION
!! From rho(r), generate rho(G), symmetrize it, and
!! come back to the real space for a symmetrized rho(r).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! cplex=1 if rhor is real, 2 if rhor is complex
!! densymop <type(dens_sym_operator_type)>=the density symmetrization operator
!! irrzon(nfft,2,nspden/nsppol)=irreducible zone data
!! mpi_enreg=informations about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
!! nspden=number of spin-density components
!! nsppol=1 for unpolarized, 2 for spin-polarized
!! nsym=number of symmetry elements.
!! phnons(2,nfft,nspden/nsppol)=nonsymmorphic translation phases
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!
!! OUTPUT
!! rhog(2,nfft)=symmetrized rho(G) (total) electron density in G space
!!
!! SIDE EFFECTS
!! Input/Output
!! rhor(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
!! Input, but also output, if symmetrization is applied.
!! Also output if nspden > 1 (change spin components)
!!
!! NOTES
!! When using spin-polarization (nspden==2),
!! put total density in first half of rhor array and spin up in second half
!! If (nspden=2 and nsppol=2) the density is transformed as  (up,down) => (up+down,up)
!! If (nspden=2 and nsppol=1) anti-ferromagnetic symmetry operations
!!  must be used, such as to transform (2*up) => (up+down,up)
!! In spin-polarized, and if there is no symmetry to be
!! applied on the system, only the total density is generated in G space
!!
!! PARENTS
!!      mkrho,mkrho3,rhofermi3,suscep_dyn,suscep_kxc_dyn,suscep_stat,vtorho
!!      vtorho3,vtorhotf
!!
!! CHILDREN
!!      fourdp,leave_new,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symrhg(cplex,densymop,irrzon,mpi_enreg,nfft,nfftot,ngfft,nspden,nsppol,nsym,paral_kgb,&
& phnons,rhog,rhor,symafm)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
 use interfaces_12ffts
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nfftot,nspden,nsppol,nsym,paral_kgb
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dens_sym_operator_type),intent(in) :: densymop
!arrays
 integer,intent(in) :: irrzon(nfftot**(1-1/nsym),2,nspden/nsppol),ngfft(18)
 integer,intent(in) :: symafm(nsym)
 real(dp),intent(in) :: phnons(2,nfftot**(1-1/nsym),nspden/nsppol)
 real(dp),intent(inout) :: rhor(cplex*nfft,nspden)
 real(dp),intent(out) :: rhog(2,nfft)

!Local variables-------------------------------
!scalars
 integer :: ier,ifft,imagn,ind,ipw,ispden,isym,iup,izone,izone_max,j,j1,j2,j3
 integer :: me_fft,n1,n2,n3,nd2,nproc_fft,nsym_used,numpt,nup,old_paral_level
 integer :: r2,rep,spaceComm
 real(dp) :: rhosu1,rhosu2
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)
 real(dp),allocatable :: rhosu1_arr(:),rhosu2_arr(:),work(:)

!*************************************************************************
!
!Note the timing channel 17 excludes the
!different Fourier transforms

!write(6,*)' symrhg : enter'
!write(6,*)' nsym,nspden,nsppol=',nsym,nspden,nsppol
 if(nspden==4 .and. nsym/=1)then
  write(6,*)' symrhg : does not work yet for nspden=4 and nsym/=1'
  stop
 end if

 allocate(work(cplex*nfft))

!Special treatment for spin-polarized case
 if(nspden==2 .and. nsppol==2) then

  call timab(17,1,tsec)


! When nspden=2 and nsppol=2, put total density in first half
! of rhor array and spin up in second half  (up,down) => (up+down,up)
  call timab(17,1,tsec)
  work(:)=rhor(:,1)               ! up => work
  rhor(:,1)=rhor(:,1)+rhor(:,2)   ! up+down
  rhor(:,2)=work(:)               ! work => up

  call timab(17,2,tsec)

 end if

 if(nspden==2 .and. nsppol==1) then

  call timab(17,1,tsec)

! When nspden=2 and nsppol=1, (2*up) => (2*up,up)
! Indeed, what was delivered to the present routine is a "total" density,
! obtained from occupation numbers varying between 0 and 2,
! but for spin up only potential.
  rhor(:,2)=half*rhor(:,1)

  call timab(17,2,tsec)

 end if

 if(nspden==4) then
  call timab(17,1,tsec)
  rhor(:,1)=rhor(:,1)+rhor(:,4) !nup+ndown
  rhor(:,2)=rhor(:,2)-rhor(:,1) !mx (n+mx-n)
  rhor(:,3)=rhor(:,3)-rhor(:,1) !my (n+my-n)
  rhor(:,4)=rhor(:,1)-two*rhor(:,4) !mz=n-2ndown
  call timab(17,2,tsec)
 end if

!DEBUG
!write(6,*)'  ispden  ifft   rhor(ifft,ispden)'
!do ifft=1,nfft,123
!do ispden=1,nspden
!if(cplex==1)write(6,*)ispden,ifft,rhor(ifft,ispden)
!if(cplex==2)write(6,*)ispden,ifft,rhor(2*ifft-1,ispden),rhor(2*ifft,ispden)
!end do
!if(cplex==1)write(6,*)3,ifft,rhor(ifft,1)-rhor(ifft,2)
!if(cplex==2)write(6,*)3,ifft,rhor(2*ifft-1,1)-rhor(2*ifft-1,2),&
!&  rhor(2*ifft,1)-rhor(2*ifft,2)
!end do
!write(6,*)' symrhg : leave'
!ENDDEBUG

 if(nsym==1)then

  if(nspden==2 .and. nsppol==1) then
!  There must be at least one anti-ferromagnetic operation
   write(message,'(a,a,a)') ' symrhg : BUG -',ch10,&
&   ' In the antiferromagnetic case, nsym cannot be 1'
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if

! If not using symmetry, still want total density in G space rho(G).
! Fourier transform (incl normalization) to get rho(G)
  work(:)=rhor(:,1)
  call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
 else

! Treat either full density or spin-up density
! Note the decrease of ispden to the value 1, in order to finish
! with rhog of the total density (and not the spin-up density)
  do ispden=nspden,1,-1

!  Prepare the density to be symmetrized, in the reciprocal space
   if(nspden==1 .or. nsppol==2)then
    imagn=1
    nsym_used=0
    do isym=1,nsym
     if(symafm(isym)==1)nsym_used=nsym_used+1
!    DEBUG
!    write(6,*)' symrhg : isym,symafm(isym)',isym,symafm(isym)
!    ENDDEBUG
    end do
   else if(nspden==2 .and. nsppol==1)then   ! antiferromagnetic case
    imagn=ispden
    nsym_used=nsym/ispden
   end if

!  DEBUG
!  write(6,*)' symrhg : nsym_used=',nsym_used
!  ENDDEBUG

!  rhor -fft-> rhog    (rhog is used as work space)
!  Note : it should be possible to reuse rhog in the antiferromagnetic case
!  this would avoid one FFT
   work(:)=rhor(:,ispden)
   call fourdp(cplex,rhog,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)

!  Begins the timing here only , to exclude FFTs
   call timab(17,1,tsec)

   n1=ngfft(1);n2=ngfft(2);n3=ngfft(3);nproc_fft=ngfft(10);me_fft=ngfft(11);nd2=n2/nproc_fft

!  DEBUG
!  write(6,*)' symrhg : fourier space density'
!  ind=0
!  do j3=1,n3
!  do j2=1,n2
!  if(((j2-1)/nd2)==me_fft) then
!  do j1=1,n1
!  ind=ind+1
!  write(6,'(5i4,2es16.6)') j1-1,j2-1,j3-1,ind,n1*(n2*(j3-1)+j2-1)+j1,rhog(:,ind)
!  end do
!  end if
!  end do
!  end do
!  ENDDEBUG

!  DEBUG
!  phnons(2,:,1)=zero
!  write(6,*)' symrhg : density before symmetrization, phnons,irrzon'
!  do ipw=1,nfft
!  j=ipw-1;j1=modulo(j,n1);r2=modulo(j/n1,nd2);j3=j/(n1*nd2);j2=me_fft*nd2+r2
!  ind=n1*(n2*j3+j2)+j1+1 !this is ind in the full array proc
!  write(6, '(6i4,4es16.6,2i6)' )ipw,j1,j2,j3,r2,ind,rhog(:,ipw) !,phnons(:,ipw,1),irrzon(ipw,:,1)
!  write(6, * )ipw,rhog(:,ipw) !,phnons(:,ipw,1),irrzon(ipw,:,1)
!  write(6,'8i4')j+1,j1,j2,j3,r2,ind
!  end do
!  end do
!  ENDDEBUG

!  Get maxvalue of izone
   do izone=1,nfftot
!   Get repetition number
    rep=irrzon(izone,2,imagn)
    if(rep==0)exit
   end do
   izone_max=izone
   allocate(rhosu1_arr(izone_max),rhosu2_arr(izone_max))

   numpt=0
   do izone=1,nfftot

!   Get repetition number
    rep=irrzon(izone,2,imagn)
    if(rep==0)exit

!   Compute number of unique points in this symm class:
    nup=nsym_used/rep

!   Accumulate charge over equivalent points
    rhosu1=0.0_dp
    rhosu2=0.0_dp
    do iup=1,nup
     ind=irrzon(iup+numpt,1,imagn)
     j=ind-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);r2=modulo(j2,nd2)
     if(modulo(j/n1,n2)/nd2==me_fft) then ! this ind is to be treated by me_fft
      ind=n1*(nd2*j3+r2)+j1+1 !this is ind in the current proc
      rhosu1=rhosu1+rhog(1,ind)*phnons(1,iup+numpt,imagn)&
&      -rhog(2,ind)*phnons(2,iup+numpt,imagn)
      rhosu2=rhosu2+rhog(2,ind)*phnons(1,iup+numpt,imagn)&
&      +rhog(1,ind)*phnons(2,iup+numpt,imagn)
     end if

    end do
    rhosu1=rhosu1/dble(nup)
    rhosu2=rhosu2/dble(nup)
    rhosu1_arr(izone)=rhosu1
    rhosu2_arr(izone)=rhosu2
!   Keep index of how many points have been considered:
    numpt=numpt+nup

!   End loop over izone
   end do

   if(mpi_enreg%mode_para=='b')then
    old_paral_level=mpi_enreg%paral_level
    mpi_enreg%paral_level=3
    spaceComm=mpi_enreg%comm_fft
    call xsum_mpi(rhosu1_arr,spaceComm,ier)
    call xsum_mpi(rhosu2_arr,spaceComm,ier)
    mpi_enreg%paral_level=old_paral_level
   end if

!  Now symmetrize the density
   numpt=0
   do izone=1,nfftot

!   Get repetition number
    rep=irrzon(izone,2,imagn)
    if(rep==0)exit

!   Compute number of unique points in this symm class:
    nup=nsym_used/rep


!   Define symmetrized rho(G) at equivalent points:
    do iup=1,nup
     ind=irrzon(iup+numpt,1,imagn)
!    decompose ind-1=n1(n2 j3+ j2)+j1
     j=ind-1;j1=modulo(j,n1);j2=modulo(j/n1,n2);j3=j/(n1*n2);r2=modulo(j2,nd2)
     if(modulo(j/n1,n2)/nd2==me_fft) then ! this ind is to be treated by me_fft
!     ind in the proc ind-1=n1(nd2 j3+ r2)+j1
      ind=n1*(nd2*j3+r2)+j1+1 !this is ind in the current proc
      rhog(1,ind)=rhosu1_arr(izone)*phnons(1,iup+numpt,imagn)&
&      +rhosu2_arr(izone)*phnons(2,iup+numpt,imagn)
      rhog(2,ind)=rhosu2_arr(izone)*phnons(1,iup+numpt,imagn)&
&      -rhosu1_arr(izone)*phnons(2,iup+numpt,imagn)
     end if
    end do

!   Keep index of how many points have been considered:
    numpt=numpt+nup

!   End loop over izone
   end do

!  DEBUG
!  write(6,*)' symrhg : density after symmetrization, phnons,irrzon'
!  do ipw=1,nfft
!  if(abs(rhog(1,ipw))<1.0d-14)rhog(1,ipw)=0.0_dp
!  if(abs(rhog(2,ipw))<1.0d-14)rhog(2,ipw)=0.0_dp
!  write(6, *)ipw,rhog(:,ipw) !,phnons(:,ipw,1),irrzon(ipw,:,1)
!  end do
!  ENDDEBUG

   call timab(17,2,tsec)

!  Pull out full or spin up density, now symmetrized
   call fourdp(cplex,rhog,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   rhor(:,ispden)=work(:)
   deallocate(rhosu1_arr,rhosu2_arr)
  end do

! End on the condition nsym==1
 end if

!DEBUG
!write(6,*)'  ispden  ifft   rhor(ifft,ispden)'
!do ifft=1,nfft,123
!do ispden=1,nspden
!if(cplex==1)write(6,*)ispden,ifft,rhor(ifft,ispden)
!if(cplex==2)write(6,*)ispden,ifft,rhor(2*ifft-1,ispden),rhor(2*ifft,ispden)
!end do
!if(cplex==1)write(6,*)3,ifft,rhor(ifft,1)-rhor(ifft,2)
!if(cplex==2)write(6,*)3,ifft,rhor(2*ifft-1,1)-rhor(2*ifft-1,2),&
!&  rhor(2*ifft,1)-rhor(2*ifft,2)
!end do
!write(6,*)' symrhg : leave'
!ENDDEBUG

 deallocate(work)

end subroutine symrhg
!!***
