!{\src2tex{textfont=tt}}
!!****f* ABINIT/approxdelta
!! NAME
!!  approxdelta
!!
!! FUNCTION
!!  Approximate the Dirac function using :
!!  method 1) a triangular funtion centered at the value egwdiff_re (Eq 17 of PRB 74, 035101 (2006) 
!!  method 2) a gaussian of witdth ep%spsmear expandended in Taylor serie 
!!  (at the moment only the 0-th moments) 
!!
!!  Subroutine needed to implement the calculation 
!!  of the polarizability using the spectral representation as proposed in :
!!  PRB 74, 035101 (2006) and PRB 61, 7172 (1999)
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nomegasf=number of frequencies in the grid for Im \chi_0
!!  omegasf(0:nomega+1)= frequencies (real)
!!  egwdiff_re = transition energy where the delta function is centered 
!!
!!  method= 1 : a triangular shaped function used to approximated the delta
!!          2 : gaussian approximation with standard deviation (smear)
!! smear= used only in case of method==2, defines the width of the gaussian
!!
!! OUTPUT
!!  wl = weight associated to omegal (last omega wich is smaller than egwdiff_re
!!  wr = weight associate to omegar  (first omega larger than egwdff_re
!!  iomegal= index in the array omegasf of the last frequency < egwdiff
!!  iomegar= index in the array omegasf of the first frequency > egwdiff
!!
!! NOTES 
!!  Other approximations could be tried
!!
!! PARENTS
!!     
!!
!! CHILDREN
!!    
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine approxdelta(nomegasf,omegasf,egwdiff_re,smear,iomegal,iomegar,wl,wr,spmeth)
    
 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nomegasf,spmeth
 integer,intent(out) :: iomegal,iomegar
 real(dp),intent(in) :: egwdiff_re,smear  
 real(dp),intent(out) :: wl,wr
!arrays
 real(dp),intent(in) :: omegasf(nomegasf)

!Local variables-------------------------------
 integer :: io,iomega
 real(dp) :: omegal,omegar,deltal,deltar
 character(len=500) :: msg                 
! *************************************************************************
 
 iomega=-999
 do io=nomegasf,1,-1
  if (omegasf(io)<egwdiff_re) then
   iomega=io  
   exit
  end if 
 end do

 if (iomega==nomegasf) then 
  write(msg,'(a)')' approxdelta BUG1 ' 
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end if 
 if (iomega==-999) then 
  write(msg,'(a)')' approxdelta BUG2 ' 
  call wrtout(std_out,msg,'COLL')
  write(*,*)REAL(omegasf(:)),'delta=',egwdiff_re
  call leave_new('COLL')
 end if 
 
 iomegal=iomega    ; omegal=omegasf(iomegal)
 iomegar=iomegal+1 ; omegar=omegasf(iomegar)

 select CASE (spmeth) 
  CASE (1) 
   ! Weights for triangular shaped function
   wr=  (egwdiff_re-omegal)/(omegar-omegal)
   wl= -(egwdiff_re-omegar)/(omegar-omegal)
  CASE (2) 
   !  weights for gaussian method (0-th moment)
   deltal=(egwdiff_re-omegal)/smear
   deltar=(omegar-egwdiff_re)/smear
   if (deltar>=deltal) then
    wl=EXP(-deltal*deltal)
    ! this value is used to avoid double counting and speed-up
    wr=huge(0.0_dp)*1.d-10 
   else 
    wl=huge(0.0_dp)*1.d-10
    wr=exp(-deltal*deltal) 
   end if 
  CASE DEFAULT
   write(msg,'(4a,i4)')ch10,&
&   ' approxdelta : BUG - ',ch10,&
&   ' wrong value for spmeth = ',spmeth 
   call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 end select

end subroutine approxdelta
!!***

!!****f* ABINIT/calc_kkweight
!! NAME
!!  calc_kkweight
!!
!! FUNCTION
!!  Calculate frequency dependent weights needed to perform the Hilbert transform 
!!
!!  Subroutine needed to implement the calculation 
!!  of the polarizability using the spectral representation as proposed in :
!!  PRB 74, 035101 (2006) and PRB 61, 7172 (1999)
!!  
!! COPYRIGHT
!!  Copyright (C) 2005-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nsp=number of frequencies where the imaginary part of the polarizability is evaluated
!! ne=number of frequencies for the polarizability (same as in epsilon^-1)
!! omegasp(nsp)=real frequencies for the imaginary part of the polarizability 
!! omegae(ne)= imaginary frequencies for the polarizability
!! delta=small imaginary part used to avoid poles, input variables
!!
!! OUTPUT
!! kkweight(ne,nsp)=frequency dependent weights (Eq A1 PRB 74, 035101 (2006)
!!
!! NOTES 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine calc_kkweight(ne,omegae,nsp,omegasp,delta,omegamax,kkw)
    
 use defs_basis
 use m_gwdefs, only : j_dpc
 use m_io_tools, only : get_unit, flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ne,nsp
 real(dp),intent(in) :: delta,omegamax
!arrays
 real(dp),intent(in) :: omegasp(nsp)
 complex(gwpc),intent(in) :: omegae(ne)
 complex(dpc),intent(out) :: kkw(ne,nsp)

!Local variables-------------------------------
!scalars
 integer :: isp,je,unt
 real(dp) :: alpha,beta,eta,xx1,xx2,wl,wr
 real(dp) :: den1,den2
 complex(dpc) :: c1,c2,wt
 logical :: paranoia
 character(len=500) :: msg
!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' calc_kkweight: enter '
 call wrtout(std_out,msg,'COLL') 
 call flush_unit(std_out)
#endif
 
 kkw(:,:)=czero
 do je=1,ne
  eta=delta
  wt=omegae(je)
  ! Not include shift at omega==0, what about metallic systems?
  if (abs(real(omegae(je)))<tol6 .and. abs(aimag(wt))<tol6) eta=tol12
  !  Not include shift along the imaginary axis
  if (abs(aimag(wt))>tol6) eta=zero
  do isp=1,nsp
    if (isp==1) then  
     ! Skip negative point, should check that this would not lead to spurious effects
     c1=czero
     den1=one
    else 
     xx1=omegasp(isp-1)
     xx2=omegasp(isp)
     den1= xx2-xx1 
     c1= -(wt-xx1+j_dpc*eta)*log( (wt-xx2+j_dpc*eta)/(wt-xx1+j_dpc*eta) )&
&        +(wt+xx1-j_dpc*eta)*log( (wt+xx2-j_dpc*eta)/(wt+xx1-j_dpc*eta) )
     c1= c1/den1
    end if 
    xx1=omegasp(isp)
    if (isp==nsp) then 
     ! Skip last point should check that this would not lead to spurious effects 
     xx2=omegamax
    else
     xx2=omegasp(isp+1)
    end if
    den2=xx2-xx1
    c2=  (wt-xx2+j_dpc*eta)*log( (wt-xx2+j_dpc*eta)/(wt-xx1+j_dpc*eta) )&
&       -(wt+xx2-j_dpc*eta)*log( (wt+xx2-j_dpc*eta)/(wt+xx1-j_dpc*eta) )
    c2= c2/den2
    kkw(je,isp)=  c1/den1 + c2/den2
   end do
  end do 

#if defined DEBUG_MODE
 if (paranoia) then 
  write(77,*)'kweights method '
  unt=get_unit()
  open(unt)
  do je=1,ne 
   do isp=1,nsp
    write(unt,*)je,isp,kkw(je,isp)
   end do
   write(unt,*)
  end do 
 end if 
 write(msg,'(2a)')' calc_kkweight: exit ',ch10
 call wrtout(std_out,msg,'PERS') 
 call flush_unit(std_out)
#endif 


end subroutine calc_kkweight
!!***


!{\src2tex{textfont=tt}}
!!****f* ABINIT/setup_spectral
!! NAME
!!  setup_spectra
!!
!! FUNCTION
!! Calculation of \chi_o based on the spectral method as proposed in PRB 74, 035101 (2006) and PRB 61, 7172 (1999).
!! Setup of the real frequency mesh for $\Im\chi_o$ and of the frequency-dependent weights for 
!! Hilbert transform. Note that CPU time does not depend dramatically on nomegasf unlike memory.
!! spmeth defines the approximant for the delta function:
!!  ==1 : use Triangular approximant (Kresse method)
!!  ==2 : use Gaussian method, requiring smearing (Miyake method)
!!
!! COPYRIGHT
!!  Copyright (C) 2005-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! nomegasf=number of points for the imaginary part of $\chi0(q,\omega)$ 
!! nomega=number of frequencies in $\chi0(q,\omega)$.
!! max_rest,min_res=max and min resonant transition energy (for this q-point)
!! my_max_rest,my_min_rest=max and min resonant transition energy treated by this processor
!! method=integer flag defining the type of frequency mesh used for $\Im chi0$
!!  | 0 for a linear mesh 
!!  | 1 for a mesh densified around omegaplasma
!! omegaplasma=frequency around which the mesh is densifies (usually Drude plasma frequency)
!!  used only in case of method==1
!! zcut=small imaginary shift to avoid pole in chi0
!!
!! OUTPUT
!!  kkweight(nomega,nomegasf)=Frequency dependent weight for Hilber transform.
!!  omegasf(nomegasf+1)=frequencies for imaginary part.
!!
!! NOTES 
!!
!! PARENTS
!!     
!!
!! CHILDREN
!!    
!!
!! SOURCE
!!

 subroutine setup_spectral(nomega,omega,nomegasf,omegasf,max_rest,min_rest,my_max_rest,my_min_rest,&
& method,zcut,omegaplasma,my_wl,my_wr,kkweight)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => setup_spectral
 use interfaces_lib00numeric
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,nomega,nomegasf
 integer,intent(out) :: my_wl,my_wr
 real(dp),intent(in) :: max_rest,min_rest,omegaplasma,zcut
 real(dp),intent(in) :: my_max_rest,my_min_rest
!arrays
 real(dp),intent(out) :: omegasf(nomegasf)
 complex(gwpc),intent(in) :: omega(nomega)       
 complex(dpc),intent(out) :: kkweight(nomega,nomegasf)

!Local variables-------------------------------
!scalars
 integer :: io,ii
 real(dp) :: nu_min,nu_max,nu1,nu2,dd,domegasf,wp,deltat
 character(len=500) :: msg                 
!arrays
 integer,allocatable :: insort(:)
!************************************************************************

 ! === The mesh has to enclose the entire range of transitions ===
 dd=(max_rest-min_rest)/(nomegasf-1)
 domegasf=(max_rest-min_rest+2*dd)/(nomegasf-1) 

 write(msg,'(4a,f8.3,3a,i5,2a,f8.5,a)')ch10,&
& ' ** Info on the real frequency mesh for spectral method ** ',ch10,&
& '  maximum frequency = ',max_rest*Ha_eV,' [eV]',ch10,&
& '  nomegasf = ',nomegasf,ch10,&
& '  domegasf = ',domegasf*Ha_eV,' [eV]'
 call wrtout(std_out,msg,'COLL') !; call wrtout(ab_out,msg,'COLL')

 if (min_rest<tol6) then  
  write(*,*)" system seems to be metallic"
 end if
 !
 ! === Generate mesh for imaginary part accordind to method ===
 SELECT CASE (method) 
 CASE (0)
  ! === Linear mesh ===
  write(msg,'(a)')'  mesh is linear ' ; call wrtout(std_out,msg,'COLL')
  do io=1,nomegasf
   omegasf(io)=(io-1)*domegasf+min_rest-dd
  end do 
 CASE (1) 
  !===  Non-homogeneous mesh densified around omega_plasma, do not improve results ===
  ! WARNING_ this part has to be checked since I modified omegasf
  write(msg,'(a,f7.4,a)')' mesh is densified around ',omegaplasma*Ha_eV,' [eV] '
  call wrtout(std_out,msg,'COLL')
  wp=omegaplasma ; deltat=max_rest-min_rest 
  nu_min=zero
  if (deltat<wp ) then 
   nu_max = wp/sqrt2 *   ATAN(sqrt2*deltat*wp/(-deltat**2+wp**2))
  else 
   nu_max = wp/sqrt2 * ( ATAN(sqrt2*deltat*wp/(-deltat**2+wp**2)) + pi)
  end if 
  domegasf=(nu_max-nu_min)/(nomegasf+1)
  !write(*,*)  -(wp/sqrt2) * atan(sqrt2*deltat*wp/(deltat**2-wp**2))
  omegasf(1)=zero ; omegasf(nomegasf+1)=deltat
  ii=0
  do io=2,nomegasf
   nu1=domegasf*(io-1) ; nu2=TAN(-sqrt2*nu1/wp)
   if (nu2<0) then 
    omegasf(io) = wp * (one - SQRT(1+2*nu2**2))/(sqrt2*nu2)
   else
    omegasf(io) = wp * (one + SQRT(1+2*nu2**2))/(sqrt2*nu2) 
   end if 
   if (omegasf(io)> deltat ) then  
    omegasf(io)= deltat-0.1*ii
    ii=ii+1
   end if 
   ! write(102,'(i4,2x,3(f9.4,2x))')io,nu1,nu2,ep%omegasf(io)*Ha_eV
  end do 
  ! Reorder frequencies in ascending order
  allocate(insort(nomegasf+1)) ; insort(:)=(/ (io,io=1,nomegasf+1) /)
  call sort_dp(nomegasf+1,omegasf,insort,tol14)
  deallocate(insort)
 CASE DEFAULT 
  write(msg,'(4a)')ch10,&
&  ' setup_spectral : BUG- ',ch10,&
&  ' called with wrong values for method '
  call wrtout(std_out,msg,'COLL') ; call leave_new('COLL')
 END SELECT
 write(*,*)omegasf(1)*Ha_eV,omegasf(nomegasf)*Ha_eV
 !
 ! === Find min and max index in omegasf treated by this processor ===
 my_wr=-999
 do io=1,nomegasf
  if (omegasf(io)>my_max_rest) then 
   my_wr=io ; EXIT
  end if 
 end do
 if (my_wr==nomegasf+2) my_wr=nomegasf+1
 my_wl=-999
 do io=nomegasf,1,-1
  if (omegasf(io)< my_min_rest) then ! Check metals
   my_wl=io ; EXIT
  end if 
 end do 

 if (my_wl==-999 .or. my_wr==-999) then 
  write(msg,'(4a,2i6)')ch10,&
&  ' setup_spectral: ERROR - ',ch10,&
&  ' wrong value in my_wl and/or my_wr ',my_wl,my_wr
  call wrtout(std_out,msg,'PERS') ; call leave_new('COLL')
 end if 
 ! 
 ! === Calculate weights for Hilbert transform ===
 call calc_kkweight(nomega,omega,nomegasf,omegasf,zcut,max_rest,kkweight)

end subroutine setup_spectral
!!***
