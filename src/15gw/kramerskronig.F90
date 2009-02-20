!{\src2tex{textfont=tt}}
!!****f* ABINIT/kramerskronig
!! NAME
!! kramerskronig
!!
!! FUNCTION
!! check or apply the Kramers Kronig relation
!!
!! Re \epsilon(\omega) = 1 + \frac{2}{\pi}
!! \int_0^\infty d\omega' frac{\omega'}{\omega'^2 - \omega^2} Im \epsilon(\omega') 
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (Valerio Olevano, Lucia Reining, Francesco Sottile, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  nomega=number of real frequencies
!!  omega(nomega)= real frequencies 
!!  eps(nomega)= function on the frequency grid (both real and imaginary part)
!!   real part can be used to check whether the K-K relation is satisfied or not
!!  method=method used to perform the integration
!!   0= naive integration
!!   1=simpson rule 
!!  only_check= if /=0 the real part of eps is checked against the imaginary part,
!!                a final report in written but the array eps is not modified 
!!              if ==0 the real part of eps is overwritten using the
!!              results obtained using the Kramers-Kronig relation     
!!
!! OUTPUT
!!  
!!
!! SIDE EFFECTS
!!
!! NOTES
!! Inspired to check_kramerskronig of the DP code 
!!
!! PARENTS
!!      linear_optics_paw
!!
!! CHILDREN
!!      leave_new,simpson_int,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine kramerskronig(nomega,omega,eps,method,only_check)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_14occeig
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: method,nomega,only_check
!arrays
 real(dp),intent(in) :: omega(nomega)
 complex,intent(inout) :: eps(nomega)

!Local variables-------------------------------
!scalars
 integer,save :: enough=0
 integer :: ii,ip
 real(dp) :: acc,domega,eav,kkdif,kkrms,ww,wwp
 character(len=500) :: message
!arrays
 real(dp) :: e1kk(nomega),intkk(nomega),kk(nomega)

! *************************************************************************
 
!DEBUG
!write (std_out,*) ' kramerskronig : enter'
!ENDDEBUG

!Check whether the frequency grid is linear or not 
 domega = (omega(nomega) - omega(1)) / (nomega-1)
 do ii = 2, nomega
  if(abs(domega-(omega(ii)-omega(ii-1))) > 0.001) then 
   if (only_check/=1) then
    write(message,'(5a)')ch10,&
&   ' kramerskronig : WARNING -',ch10,&
&   '  check cannot be performed since frequency step is not constant',ch10
    call wrtout(6,message,'COLL')
    return
   else 
    write(message,'(5a)')ch10,&
&    ' kramerskronig : ERROR -',ch10,&
&    '  cannot perform integration since frequency step is not constant',ch10
    call wrtout(6,message,'COLL')
    call leave_new('COLL')
   end if 
  end if
 end do

!Check whether omega(1) is small or not
 if(omega(1) > 0.1/Ha_eV) then 
  if (only_check/=1) then
   write(message,'(5a)')ch10,&
&  ' kramerskronig : WARNING -',ch10,&
&  '  check cannot be performed since first frequency on the grid > 0.1 eV',ch10
   call wrtout(6,message,'COLL')
   return
  else 
   write(message,'(5a)')ch10,&
&   ' kramerskronig : ERROR -',ch10,&
&   '  cannot perform integration since first frequency on the grid > 0.1 eV',ch10
   call wrtout(6,message,'COLL')
   call leave_new('COLL')
  end if
 end if 

!If eps(nomega) is not 0 warn
 if(aimag(eps(nomega)) > 0.1 .and. enough<50) then
  enough=enough+1
  write(message,'(4a,f8.4,3a,f8.2,2a)')ch10,&
&  ' kramerskronig : WARNING -',ch10,&
&  '  Im epsilon for omega = ',omega(nomega)*Ha_eV,' eV',ch10,&
&  '  is not yet zero, epsilon_2 = ',aimag(eps(nomega)),ch10,&
&  '  Kramers Kronig could give wrong results'
  call wrtout(6,message,'COLL')
  if (enough==50) then 
   write(message,'(3a)')&
&   ' sufficient number of WARNINGS-',ch10,' stop writing '
   call wrtout(6,message,'COLL')
  end if 
 end if
 
 
!Perform Kramers-Kronig using naive integration
 if (method==0) then  

  do ii=1,nomega
   ww = omega(ii)
   acc = 0.0_dp
   do ip=1,nomega
    if(ip == ii) cycle
    wwp = omega(ip)
    acc = acc + wwp/(wwp**2-ww**2) *aimag(eps(ip))
   end do
   e1kk(ii) = 1.0 + 2.0/pi*domega* acc
  end do
  
! Perform Kramers-Kronig using Simpson integration    
! Simpson O(1/N^4), from NumRec in C p 134  NumRec in Fortran p 128
 else if (method==1) then 

  kk=0.0_dp

  do ii =1,nomega
   ww=omega(ii)
   do ip=1,nomega
    if(ip == ii) cycle
    wwp = omega(ip)
    kk(ip) = wwp/(wwp**2-ww**2) *aimag(eps(ip))
   end do
!  acc = simpson(nomega,domega,intkk)
   call simpson_int(nomega,domega,kk,intkk)
   e1kk(ii) = 1.0 + 2.0/pi * intkk(nomega)
  end do

 else 
  write(message,'(5a)')ch10,&
&  ' kramerskronig : BUG - ',ch10,&
&  ' wrong value for method ',ch10
  call wrtout(6,message,'COLL')
  call leave_new('COLL')
 end if 

! at this point real part is in e1kk, need to put it into eps
 do ii=1,nomega
  eps(ii)=cmplx(e1kk(ii),aimag(eps(ii)))
 end do 

!Verify Kramers-Kronig
 eav = 0.0
 kkdif = 0.0
 kkrms = 0.0

 do ii=1,nomega
  kkdif = kkdif + abs(real(eps(ii)) - e1kk(ii))
  kkrms = kkrms + (real(eps(ii)) - e1kk(ii))*(real(eps(ii)) - e1kk(ii))
  eav = eav + abs(real(eps(ii)))
 end do

 eav = eav/nomega
 kkdif = (kkdif/nomega) / eav
 kkrms = (kkrms/nomega) / (eav*eav)

 kk = abs(real(eps(1)) - e1kk(1)) / real(eps(1))

!Write data
 write(*,'(" the Kramers-Kronig is verified within ",f7.2,"%")') maxval(kk)*100
! write(*,'(a,ES12.4)')' The Kramers-Kronig transform is verified within ',kk

!open(22,file='out.kk',status='unknown')
!write(22,'("# Kramers Kronig calculation of epsilon1")')
!write(22,'("# omega   epsilon1  epsilon1kk")')
!do ii=1,nomega
!write(22,'(f7.3,2e15.7)') omega(i)*Ha_eV, real(eps(ii)),e1kk(ii)
!end do
!close(22)

end subroutine kramerskronig
!!***
