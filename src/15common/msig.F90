!{\src2tex{textfont=tt}}
!!****f* ABINIT/msig
!! NAME
!! msig
!!
!! FUNCTION
!! This program computes the elements of the optical frequency dependent
!! conductivity tensor and the conductivity along the three principal axes
!! from the Kubo-Greenwood formula for PAW formalism
!!
!! COPYRIGHT
!! Copyright (C) 2002-2008 ABINIT group (SMazevet)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!!  fcti(npti)=  conductivity, as calculated in conducti
!!  npti= number of points to calculate conductivity
!!  xi(npti)= energies where the conductivity is calculated
!!
!! OUTPUT
!!   no output, only files
!!
!! NOTES
!!     this program calculates the imaginary part of the conductivity (principal value)
!!     +derived optical properties. 
!!     the calculation is performed on the same grid as the initial input
!!     to calculate the principal value, a trapezoidale integration +taylor expansion to 
!!     third order is used (W.J. Thomson computer in physics vol 12 p94 1998)
!!    two input files are needed inppv.dat (parameters) and sigma.dat (energy,sigma_1)
!!     two output files ppsigma.dat (energy,sigma_1,sigma_2,epsilon_1,epsilon_2)
!!                      abs.dat     (energy,nomega,komega,romega,absomega)
!!     march 2002 s.mazevet
!!
!!
!!
!! PARENTS
!!      conducti_paw
!! CHILDREN
!!      intrpl
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine msig (fcti,npti, xi)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments -----------------------------------
!scalars
 integer :: npti
!arrays
 real(dp) :: fcti(npti),xi(npti)

!Local variables-------------------------------
!scalars
 integer :: High,i,ip,j,npt1,npt2,npt=10000
 real(dp),parameter :: c=2.99792458e10,ohmtoau=2.18e-5,ohmtosec=9e11
 real(dp) :: del=1.e-3,dx,dx1,dx2,eps1,eps2,idel,idel2,komega,pole,refl,sigma2
 real(dp) :: sum3,sum4,sum=0.d0
!arrays
 real(dp),allocatable :: abso(:),fct(:),fct1(:),fct2(:),fct3(:),fct4(:),fct5(:)
 real(dp),allocatable :: fctii(:),fp(:),fpp(:),fppp(:),nomega(:),ppsig(:),x(:)
 real(dp),allocatable :: x1(:),x2(:)

! *********************************************************************************
!BEGIN EXECUTABLE SECTION


 High=selected_real_kind (11,30)
 
 
 open(unit=11, file="ppsigma.dat", status="replace")
 open(unit=12, file="abs.dat", status="replace")


 
 
 print"(A)","Calculate the principal value and related optical properties"
 print"(A)","following W.J. Thomson computer in physics vol 12 p94 1998 for "
 print"(A)","the principal value. S. Mazevet "
 print"(/A)","OPTIONS"
 print "(A)", "use default number of integration pts: npt=10000"
 
 print "(A)", "Use default value for delta interval: del=1e-3"

 allocate( fct(npt),x(npt),fct2(npt),fct3(npt),fct4(npt),fct5(npt),fp(npt) &
 ,fpp(npt),fppp(npt),x1(npt),x2(npt),fct1(npt),ppsig(npt) &
 ,fctii(npt),abso(npt),nomega(npt))

 
  if (npti > npt) then
    write (*,*) 'msig: input npti is too large for hard coded npt array size = ', npt
    stop
  end if

 
 print"(A)","input values read from sigma.dat (energy,conductivity) are assumed to be (eV,10^6(m ohm)^-1)"
!loop on the initial energy grid      
 do ip=1,npti
  
! adjust the interval before and after the pole to reflect range/npt interval
  sum=0.d0
  dx=(xi(npti)-xi(1))/float(npt-1)
  pole=xi(ip)     
  npt1=int((pole-del)/dx)
  dx1=0
  if(npt1/=1) dx1=(pole-del)/(npt1-1)     
  npt2=int((xi(npti)-pole-del)/dx)
  dx2=(xi(npti)-pole-del)/(npt2-1)
  
! for the moment skip the pp calculation when the pole if too close to the end of the range
  if(npt1<=1.or.npt2<=1) then

   sum=zero
   ppsig(ip)=zero

  else
   
!  define the fct for which the pp calculation is needed using xi^2-pole^2 factorization
   fctii = zero
   fctii(1:npti)=fcti*pole/(xi+pole)
   
!  define the grid on each side of the pole x1 before x2 after
   do i=1,npt1
    x1(i)=dx1*float(i-1)
   end do
   do i=1,npt2
    x2(i)=pole+del+dx2*float(i-1)
   end do

!  interpolate the initial fct fii on the new grids x1 and x2 (cubic spline)
!  write(*,*) npti,npt1
   
! MJV 6/12/2008:
! for each use of fctii should ensure that npt1 npt2 etc... are less than
! npt=len(fctii)
   call intrpl(npti,xi,fctii,npt1,x1,fct4,fct1,fct5,1)
   call intrpl(npti,xi,fctii,npt2,x2,fct3,fct2,fct5,1)
   
!  calculate the two integrals from 0-->pole-lamda and pole+lamda--> end range
!  trapezoidal integration
   do i=1,npt1     
    fct1(i)=fct4(i)/(x1(i)-pole)
   end do
   do i=1,npt2     
    fct2(i)=fct3(i)/(x2(i)-pole)
   end do
   
   do i=2,npt1-1
    sum=sum+fct1(i)*dx1
   end do
   do i=2,npt2-1
    sum=sum+fct2(i)*dx2
   end do
   sum=sum+.5d0*(fct1(1)+fct1(npt1))*dx1+.5d0*(fct2(1)+fct2(npt2))*dx2

!  calculate the first and third derivative at the pole and add the taylor expansion
   call intrpl(npti,xi,fctii,npti,xi,fct3,fct4,fct5,1)
   call intrpl(npti,xi,fct4,1,(/pole/),fp,fpp,fppp,1)
   
   idel=2.d0*fp(1)*(del)+fppp(1)*(del**3)/9.0
   sum=sum+idel

  end if
  
! calculate the derivated optical quantities and output the value 
  
  sigma2=(-2.d0/pi)*sum
  eps1=1.d0-(4.d0*pi*sigma2/(pole))
  eps2=4.d0*fcti(ip)*pi/(pole)

! A special treatment of the case where eps2 is very small compared to eps1 is needed
  if(eps2**2 > eps1**2 * tol12)then
   nomega(ip)=sqrt(0.5e0*(eps1 + sqrt(eps1**2 + eps2**2)))
   komega=sqrt(0.5e0*(-eps1 + sqrt(eps1**2 + eps2**2)))
   abso(ip)=4.d0*pi*fcti(ip)*ohmtosec/(ohmtoau*nomega(ip)*c)
  else if(eps1>zero)then
   nomega(ip)=sqrt(0.5e0*(eps1 + sqrt(eps1**2 + eps2**2)))
   komega=half*abs(eps2/sqrt(eps1))
   abso(ip)=4.d0*pi*fcti(ip)*ohmtosec/(ohmtoau*nomega(ip)*c)
  else if(eps1<zero)then
   nomega(ip)=half*abs(eps2/sqrt(-eps1))
   komega=sqrt(0.5e0*(-eps1 + sqrt(eps1**2 + eps2**2)))
   abso(ip)=two*sqrt(-eps1)*pole*ohmtosec/(ohmtoau*c)
  endif

  refl=((one-nomega(ip))**2 + komega**2)/ &
  ((one+nomega(ip))**2 + komega**2)
  
  write(11,100) Ha_eV*pole,fcti(ip)/ohmtoau,sigma2/ohmtoau,eps1,eps2
  
  write(12,100) Ha_eV*pole,nomega(ip),komega,refl,abso(ip) 
  
 end do
 print"(/A)","OUPUTS"
 print"(A)","output in ppsigma.dat: energy (eV),sigma_1(Ohm-1cm-1),sigma_2(Ohm-1cm-1),epsilon_1,epsilon_2"
 print"(A)","output in abs.dat: energy,nomega,komega,refl.,abso."
 print"(A)","abso in file abs.dat is in cm^-1"


 100  format(5e18.10)
 
 close(11)
 close(12)
 
 deallocate( fct,x,fct2,fct3,fct4,fct5,fp &
 ,fpp,fppp,x1,x2,fct1,ppsig,fctii,abso,nomega)
 
 end subroutine msig
!!***
