!{\src2tex{textfont=tt}}
!!****f* ABINIT/fermisolverec
!! NAME
!! fermisolverec
!!
!! FUNCTION
!! This routine computes the fermi energy in order to have a given number of
!! valence electrons in the recursion method, using a Ridder s Method
!! 
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  a, b2 : coefficient given by recursion
!!  nb_rec=order of recursion
!!  nb_point=number of discretization point in one dimension (=n1=n2=n3)
!!  temperature=temperature (Hartree)
!!  trotter=trotter parameter
!!  nelect=number of valence electrons (dtset%nelect)
!!  acc=accuracy for the fermi energy
!!  max_it=maximum number of iteration for the Ridder's Method
!!  longueur_tranche=number of point computed by thi proc
!!  mpi_enreg=informations about MPI parallelization
!!  rang=index of that proc in the paralellization
!!  rmet=define the metric : rprimd*(transpose(rprimd)) 
!!  inf_ucvol=infinitesimal unit cell volume
!!  tol=tolerance criteria for stopping recursion
!! 
!! OUTPUT
!! 
!! SIDE EFFECTS
!!  fermie=fermi energy
!!  rho=density, recomputed for the new fermi energy
!! 
!! PARENTS
!!      vtorhorec
!!
!! CHILDREN
!!      recursion, xcomm_init, xsum_mpi, timab
!! 
!! NOTES
!!  at this time :
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine fermisolverec(fermie,rho,a,b2,nb_rec, &
&                      temperature,trotter,nelect, &
&                      acc, max_it, &
&                      longueur_tranche,mpi_enreg,rang,&
&                      rmet,inf_ucvol,tim_fourdp)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_15recursion, except_this_one => fermisolverec
 use interfaces_lib01hidempi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: longueur_tranche,max_it,nb_rec,rang,tim_fourdp,trotter
 real(dp),intent(in) :: acc,inf_ucvol,nelect,temperature
 real(dp),intent(inout) :: fermie
 type(MPI_type),intent(inout) :: mpi_enreg
!arrays
 real(dp),intent(in) :: rmet(3,3)
 real(dp),intent(inout) :: a(0:nb_rec,1:longueur_tranche)
 real(dp),intent(inout) :: b2(0:nb_rec,1:longueur_tranche)
 real(dp),intent(inout) :: rho(1:longueur_tranche)

!Local variables-------------------------------
!scalars
 integer :: dummynfft,ierr,ii,ipoint,ipointlocal,jj,kk,nn,spaceComm
 real(dp) :: beta,fermieh,fermiel,fermiem,fermienew,nelecth,nelectl,nelectm
 real(dp) :: nelectnew,res_nelecth,res_nelectl,res_nelectm,res_nelectnew
 real(dp) :: rtrotter,ss
!arrays
 integer :: dummyngfft(18)
 real(dp) :: rhoh(longueur_tranche),rhol(longueur_tranche)
 real(dp) :: rhom(longueur_tranche),rhonew(longueur_tranche),tsec(2)
 real(dp),allocatable :: dummypot(:,:,:)

! *************************************************************************

 call timab(189,1,tsec)
 
 beta = 1/(temperature)
 if (trotter == 0) then
  rtrotter = 0.5d0
 else
  rtrotter = real(trotter,dp)
 end if
 dummyngfft=0
 dummyngfft(1:3)=1
 dummynfft=dummyngfft(1)*dummyngfft(2)*dummyngfft(3)
 allocate(dummypot(0:dummyngfft(1)-1,0:dummyngfft(2)-1,0:dummyngfft(3)-1))
 dummypot=1.d0
 
!initialisation de fermiel
 fermiel = fermie
 rhol = rho
 nelectl = 0.d0
 do ipointlocal = 1, longueur_tranche
  nelectl = nelectl + rhol(ipointlocal)
 end do
 mpi_enreg%paral_level=3
 call xsum_mpi( nelectl,mpi_enreg%commcart ,ierr)
 mpi_enreg%paral_level=2
 res_nelectl = inf_ucvol*nelectl - nelect
 
 if (res_nelectl /= 0.d0) then 
! initialisation de fermieh
  if (res_nelectl > 0.d0) then
!  trop d'electron -> fermie plus petit
   fermieh = fermie - 10*temperature
  else
!  pas assez d'electron -> fermie plus grand
   fermieh = fermie + 10*temperature
  end if
! nelecth = 0.d0
  do ipointlocal = 1,longueur_tranche
!  !$   ipoint = ipointlocal + (rang)*longueur_tranche
!  !$   ii = modulo(ipoint-1,nb_point)
!  !$   jj = modulo((ipoint-ii-1)/nb_point,nb_point)
!  !$   kk = ( (ipoint-ii-1)/nb_point - jj)/nb_point
   call recursion(dummypot,0,0,0, &
&   a(:,ipointlocal),& 
&   b2(:,ipointlocal), & 
&   rhoh(ipointlocal), &
&   nb_rec,fermieh,temperature,trotter, &
&   (/ 0.d0, 0.d0, 0.d0, 0.d0 /), &
&   tol14, &
&   0,0,&
&   mpi_enreg,dummynfft,dummyngfft,rmet,inf_ucvol,tim_fourdp)
  end do
  nelecth = 0.d0
  do ipointlocal = 1, longueur_tranche
   nelecth = nelecth + rhoh(ipointlocal) !ii,jj,kk)
  end do
  mpi_enreg%paral_level=3
  call xsum_mpi( nelecth,mpi_enreg%commcart ,ierr)
  mpi_enreg%paral_level=2
  res_nelecth = inf_ucvol*nelecth - nelect
  
! main loop   
  main : do nn=1,max_it
   
!  fermiem computation
   fermiem = 0.5d0*(fermiel+fermieh)
!  nelectm = 0.d0
   do ipointlocal = 1,longueur_tranche
!   !$    ipoint = ipointlocal + (rang)*longueur_tranche
!   !$    ii = modulo(ipoint-1,nb_point)
!   !$    jj = modulo((ipoint-ii-1)/nb_point,nb_point)
!   !$    kk = ( (ipoint-ii-1)/nb_point - jj)/nb_point
    call recursion(dummypot,0,0,0, &
&    a(:,ipointlocal),& 
&    b2(:,ipointlocal), & 
&    rhom(ipointlocal), &
&    nb_rec,fermiem,temperature,trotter, &
&    (/ 0.d0, 0.d0, 0.d0, 0.d0 /), &
&    tol14, &
&    0,0,&
&    mpi_enreg,dummynfft,dummyngfft,rmet,inf_ucvol,tim_fourdp)
   end do
   nelectm = 0.d0
   do ipointlocal = 1, longueur_tranche
    nelectm = nelectm + rhom(ipointlocal)
   end do
   mpi_enreg%paral_level=3
   call xsum_mpi( nelectm,mpi_enreg%commcart,ierr)
   mpi_enreg%paral_level=2
   res_nelectm = inf_ucvol*nelectm - nelect
   
!  new guess
   ss = sqrt(res_nelectm**2-res_nelectl*res_nelecth)
   fermienew = fermiem + (fermiem-fermiel)*sign(1.d0, res_nelectl-res_nelecth)*res_nelectm/ss
!  nelectnew = 0.d0
   do ipointlocal = 1,longueur_tranche
!   !$    ipoint = ipointlocal + (rang)*longueur_tranche
!   !$    ii = modulo(ipoint-1,nb_point)
!   !$    jj = modulo((ipoint-ii-1)/nb_point,nb_point)
!   !$    kk = ( (ipoint-ii-1)/nb_point - jj)/nb_point
    call recursion(dummypot,0,0,0, &
&    a(:,ipointlocal),& 
&    b2(:,ipointlocal), & 
&    rhonew(ipointlocal), &
&    nb_rec,fermienew,temperature,trotter, &
&    (/ 0.d0, 0.d0, 0.d0, 0.d0 /), &
&    tol14, &
&    0,0,&
&    mpi_enreg,dummynfft,dummyngfft,rmet,inf_ucvol,tim_fourdp)
   end do
   nelectnew = 0.d0
   do ipointlocal = 1, longueur_tranche
    nelectnew = nelectnew + rhonew(ipointlocal)
   end do
   mpi_enreg%paral_level=3
   call xsum_mpi( nelectnew,mpi_enreg%commcart ,ierr)
   mpi_enreg%paral_level=2
   res_nelectnew = inf_ucvol*nelectnew - nelect
   
!  is the exact result found ?
   if (res_nelectnew == 0.d0) then
    fermie = fermienew
    rho = rhonew
    exit
   end if
   
!  fermiel et fermieh for new iteration
   if (sign(res_nelectm,res_nelectnew) /= res_nelectm) then
    fermiel = fermiem
    res_nelectl = res_nelectm
    fermieh = fermienew
    res_nelecth = res_nelectnew
   else if (sign(res_nelectl,res_nelectnew) /= res_nelectl) then
    fermieh = fermienew
    res_nelecth = res_nelectnew
   else if (sign(res_nelecth,res_nelectnew) /= res_nelecth) then
    fermiel = fermienew
    res_nelectl = res_nelectnew
   end if
   
!  are we within the tolerance ?
   if ((abs(res_nelectnew) < acc).or.(nn == max_it)) then
    fermie = fermienew
    rho = rhonew
    exit
   end if
   
  end do main
  
 end if
 
 call timab(189,2,tsec)
 
end subroutine fermisolverec
!!***
