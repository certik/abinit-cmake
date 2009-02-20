!{\src2tex{textfont=tt}}
!!****f* ABINIT/mka2fQgrid
!! NAME
!! mka2fQgrid
!!
!! FUNCTION
!!  Calculate the Eliashberg function only using the phonon linewidths evaluated
!!  in the irreducible q-points of the coarse q-grid.
!!  The obtained results are useful to check the validity of the Fourier interpolation
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  elph_ds = electron-phonon dataset
!!  nunit = integer number for the output file
!!
!! OUTPUT
!!  Only write
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine mka2fQgrid(elph_ds,fname)

 use defs_basis
 use defs_datatypes
 use defs_elphon


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_14occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=fnlen),intent(in) :: fname
 type(elph_type),intent(in) :: elph_ds

!Local variables -------------------------
!scalars
 integer :: ibranch,iomega,iost,ismear,isppol,nsmear,nunit,qptirred
 integer :: unit_a2fgrid
 real(dp) :: a2f_factor,estep,gaussfactor,gaussprefactor,gaussval,lambda_iso
 real(dp) :: omega,omegalog,omegastep,smear,tc_macmill,weight,xx
 character(len=500) :: message
!arrays
 real(dp),allocatable :: a2f_1d(:),a2f_1mom(:),a2f_1mom_int(:),a2flogmom(:)
 real(dp),allocatable :: a2flogmom_int(:),eli_smear(:,:,:),tmpa2f(:)

! *********************************************************************

!DEBUG
!write(message,'(a)')' Enter mka2fQgrid'
!call wrtout(06,message,'COLL')
!ENDDEBUG

!grid for the representation of alpha^2F (same as mka2f)
!WARNING : supposing that the maximum and minimum value of frequency
!have been defined in mkelph_linwid.

 omegastep = (elph_ds%omega_max-elph_ds%omega_min)/(elph_ds%na2f-one)

 nunit = 66
 open (unit=nunit,file=fname,form='formatted',status='unknown',iostat=iost)
 if (iost /= 0) then
  write (message,'(2a)')' mka2fQgrid : ERROR- opening file ',trim(fname)
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 write (message,'(3a)')&
 '# Eliashberg function evaluated using only the irred q-points ',ch10,'#'
 call wrtout(nunit,message,'COLL')

 write (message,'(a,i5,2a,es16.8,2a,es16.8,2a,es16.8,2a)')&
& '# number of frequencies = ',elph_ds%na2f,ch10,         &
& '# omega_min = ',elph_ds%omega_min,ch10,                &
& '# omega_max = ',elph_ds%omega_max,ch10,                &
& '# step = ',omegastep,ch10,'#'
 call wrtout(nunit,message,'COLL')


 nsmear=5
 estep=0.00002_dp !0.54422767 meV

 write (message,'(a,i5,3a,f10.6,3a,f10.6,3a)')                &
 '# Using ',nsmear,' values for the gaussian smearing ',ch10,&
& '# starint from ',elph_ds%a2fsmear,' (Ha)',ch10,            &
& '# energy step of ',estep,' (Ha)',ch10,'#'
 call wrtout(nunit,message,'COLL')

!e-ph quantities will be calculated for nsmear gaussian smearing values
!starting from elph_ds%a2fsmearwith an energy step of estep Hartree

 write (message,'(3a)')'#      Smear(Ha) Lambda_Iso  isppol  <ln w> (K)    Tc_McMill (K) ',ch10,'#'
 call wrtout(nunit,message,'COLL')

 allocate (a2f_1mom(elph_ds%na2f),a2f_1mom_int(elph_ds%na2f))
 allocate (a2flogmom(elph_ds%na2f),a2flogmom_int(elph_ds%na2f))
 allocate (a2f_1d(elph_ds%na2f),tmpa2f(elph_ds%na2f))
 allocate (eli_smear(nsmear,elph_ds%nsppol,elph_ds%na2f))
 eli_smear(:,:,:)=zero

 do ismear=0,nsmear-1

  smear = elph_ds%a2fsmear+ismear*estep
  gaussprefactor = sqrt(piinv) / smear
  gaussfactor = one / smear

  do isppol=1,elph_ds%nsppol  ! spin pol channels

   a2f_1d(:) = zero
   tmpa2f(:) = zero

   do qptirred=1,elph_ds%nqptirred ! sum over irred qpoints
    do ibranch=1,elph_ds%nbranch

     if (abs(elph_ds%qgrid_data(qptirred,ibranch,isppol,1)) < tol10) cycle
     omega = elph_ds%omega_min
!    MG the weights in elph_ds%wtq(qptirred) are relative to the full grid spqpt,
!    we need the mapping qirredtofull
     weight=elph_ds%wtq(elph_ds%qirredtofull(qptirred))
     a2f_factor=weight*elph_ds%qgrid_data(qptirred,ibranch,isppol,2)/abs(elph_ds%qgrid_data(qptirred,ibranch,isppol,1))

     do iomega=1,elph_ds%na2f
      xx = (omega-elph_ds%qgrid_data(qptirred,ibranch,isppol,1))*gaussfactor
      gaussval = gaussprefactor*exp(-xx*xx)
      tmpa2f(iomega) = tmpa2f(iomega) + gaussval*a2f_factor
      omega = omega+omegastep
     end do

    end do !end ibranch do
   end do !end qptirred

   a2f_1d(:)= tmpa2f(:)/(2*pi*elph_ds%n0(isppol))
   eli_smear(ismear+1,isppol,:)=a2f_1d(:) !save values

!  Do isotropic calculation of lambda and output lambda, Tc(MacMillan)
   a2f_1mom(:) = zero
   omega = elph_ds%omega_min

   do iomega=1,elph_ds%na2f
    if (abs(omega) > tol10) a2f_1mom(iomega) = two*a2f_1d(iomega)/abs(omega)
    omega=omega+omegastep
   end do

   call simpson_int(elph_ds%na2f,omegastep,a2f_1mom,a2f_1mom_int)
   lambda_iso = a2f_1mom_int(elph_ds%na2f)

!  Get log moment of alpha^2F
   a2flogmom(:) = zero
   omega = elph_ds%omega_min
   do iomega=1,elph_ds%na2f
    if (abs(omega) > tol10) then
     a2flogmom(iomega) = (two/lambda_iso)*a2f_1d(iomega)*log(abs(omega))/abs(omega)
    end if
    omega=omega+omegastep
   end do

   call simpson_int(elph_ds%na2f,omegastep,a2flogmom,a2flogmom_int)
   omegalog = exp(a2flogmom_int(elph_ds%na2f))

   tc_macmill = (omegalog/1.2_dp) * &
&   exp((-1.04_dp*(one+lambda_iso)) / (lambda_iso-elph_ds%mustar*(one+0.62_dp*lambda_iso)))

!  write data
   write(message,'(a,5x,f10.6,f10.6,i5,2x,f12.7,2x,f12.6,2x,es16.8)')&
&   '# ',smear,lambda_iso,isppol,omegalog/kb_HaK,tc_macmill/kb_HaK
   call wrtout(nunit,message,'COLL')

  end do !end isppol

 end do !ismear

 deallocate (a2f_1mom,a2f_1mom_int)
 deallocate (a2flogmom,a2flogmom_int)

!write to file
 write(message,'(4a)')'#',ch10,'# Eliashberg function calculated for different gaussian smearing values',ch10
 call wrtout(nunit,message,'COLL')

 do isppol=1,elph_ds%nsppol
  omega = elph_ds%omega_min
  write(nunit,'(a,i5)') '# smeared alpha2F for isppol = ',isppol
  do iomega=1,elph_ds%na2f
   write(nunit,'(6(f17.12,1x))')omega,eli_smear(:,isppol,iomega)
   omega=omega+omegastep
  end do
  write(nunit,*)
 end do

 deallocate (eli_smear)
 deallocate (a2f_1d)
 deallocate (tmpa2f)

 close (nunit)
!DEBUG
!write(message,'(a)')' mka2fQgrid : exit '
!call wrtout(6,message,'COLL')
!ENDDEBUG

end subroutine mka2fQgrid
!!***
