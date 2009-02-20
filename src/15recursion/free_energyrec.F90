!{\src2tex{textfont=tt}}
!!****f* ABINIT/free_energyrec
!! NAME
!! free_energyrec
!! 
!! FUNCTION
!! This routine computes the local part of the free energy at a point using a path integral, 
!! in the recursion method.
!! 
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group ( ).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  an, bn2 : coefficient given by the recursion.
!!  nrec=order of recursion
!!  trotter=trotter parameter
!!  mult=a multiplicator for computing free energy (2 for non-spin-polarized system)
!!  prtvol=printing volume
!! 
!! OUTPUT
!!  ene_out=free energy at the point
!!  
!! PARENTS
!!      vtorhorec
!! 
!! CHILDREN
!!      wrtout
!! 
!! NOTES
!!  at this time :
!!       - mult should be not used
!!       - the routine should be integraly rewrited and use the routine recursion. 
!!       - only modified for p /= 0
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine free_energyrec(an,bn2,nrec,trotter,ene_out, mult, &
&                     prtvol,n_pt_integ,xmax,&
&                     ene_out1,ene_out2,ene_out3,ene_out4)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_00basis
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: n_pt_integ,nrec,prtvol,trotter
 real(dp),intent(in) :: mult,xmax
 real(dp),intent(out) :: ene_out,ene_out1,ene_out2,ene_out3,ene_out4
!arrays
 real(dp),intent(in) :: an(0:nrec),bn2(0:nrec)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=7
 integer,save :: first=1
 integer :: ii,jj,kk,ll,mult_n_pt_integ,n_pt_integ_path2
 real(dp) :: arg,arg_path,epsilon,step,teta,twotrotter,xmin,ymin
 complex(dpc) :: D,Dnew,Dold,N,Nnew,Nold,dz_path,ene_acc,ene_acc1,ene_acc2
 complex(dpc) :: ene_acc3,ene_acc4,ene_acc_path3,ene_acc_path56,func
 complex(dpc) :: rootof_z_path,rootof_z_path2,z_path,zj
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

 func(z_path,twotrotter) = log(cone+z_path**twotrotter)
 
 call timab(190,1,tsec)
 
!structured debugging if prtvol=-level : print detailled result the first time we enter free_energyrec
 if(prtvol==-level)then
  if(first==1)then
   write(message,'(a)')' ' 
   call wrtout(06,message,'PERS')
   write(message,'(a)')' free_energyrec : enter ' 
   call wrtout(06,message,'PERS')

   write(message,'(a,i6)')'n_pt_integ ' , n_pt_integ
   call wrtout(06,message,'COLL')

  end if
 end if
 
 ene_out = 0.d0
 ene_acc = dcmplx(0.d0,0.d0)
 ene_acc1 = dcmplx(0.d0,0.d0)
 ene_acc2 = dcmplx(0.d0,0.d0)
 ene_acc3 = dcmplx(0.d0,0.d0)
 ene_acc4 = dcmplx(0.d0,0.d0)
 
!path parameters 
 mult_n_pt_integ = 1 !25

!path parameters 
!n_pt_integ = 2500
 mult_n_pt_integ = 1 !25
 step = (xmax-xmin)/real(n_pt_integ,dp)
 if(trotter==0)then
  twotrotter = 1.d0
 else
  twotrotter = 2.d0*real(trotter,dp)
 end if
 xmin = -5.d-1
 epsilon = 1/2.d0*sin( pi/twotrotter)
!xmin = -abs(xmin)**(1.d0/twotrotter)
 
!####################################################################
![xmax + i*epsilon,xmin + i*epsilon]
 path1:  do ii = 0,n_pt_integ
  z_path = dcmplx(xmin + real(ii,dp)*(xmax-xmin)/real(n_pt_integ,dp),epsilon)
  dz_path = -dcmplx((xmax-xmin)/real(n_pt_integ,dp),0.d0)
  
  Nold = dcmplx(0.d0,0.d0)
  Dold = dcmplx(1.d0,0.d0)
  N = dcmplx(1.d0,0.d0)
  D = z_path - dcmplx(an(0),0.d0)
  
  do kk=1,nrec
   Nnew = (z_path - dcmplx(an(kk),0.d0))*N - dcmplx(bn2(kk),0.d0)*Nold
   Dnew = (z_path - dcmplx(an(kk),0.d0))*D - dcmplx(bn2(kk),0.d0)*Dold
   
   Nold = N
   Dold = D
   N = Nnew
   D = Dnew
   
   if(kk/=nrec)then
    if((bn2(kk+1)<tol14))exit
   end if
   
  end do
  
  if(ii==0.or.ii==n_pt_integ)then
   ene_acc = ene_acc + 0.5d0*func(z_path,twotrotter)*&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ene_acc1 = ene_acc1 + 0.5d0*func(z_path,twotrotter)*&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
  else
   ene_acc = ene_acc + func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ene_acc1 = ene_acc1 + func(z_path,twotrotter) *&
&   N/D *dz_path       
  end if
 end do path1
 
!####################################################################
![xmin + i*epsilon,xmin]
 if(epsilon/step>4.d0)then
  n_pt_integ_path2 = int(epsilon/step)+1
 else
  n_pt_integ_path2 = 5
 end if
 n_pt_integ_path2 = n_pt_integ
 path2:  do ii = 0,n_pt_integ_path2
  z_path = dcmplx(xmin,real(ii,dp)*epsilon/real(n_pt_integ_path2,dp))
  dz_path = -dcmplx(0.d0,epsilon/real(n_pt_integ_path2,dp))

  Nold = dcmplx(0.d0,0.d0)
  Dold = dcmplx(1.d0,0.d0)
  N = dcmplx(1.d0,0.d0)
  D = z_path - dcmplx(an(0),0.d0)
  
  do kk=1,nrec
   Nnew = (z_path - dcmplx(an(kk),0.d0))*N - dcmplx(bn2(kk),0.d0)*Nold
   Dnew = (z_path - dcmplx(an(kk),0.d0))*D - dcmplx(bn2(kk),0.d0)*Dold
   
   Nold = N
   Dold = D
   N = Nnew
   D = Dnew
   
   if(kk/=nrec)then
    if((bn2(kk+1)<tol14))exit
   end if
   
  end do
  
  if(ii==0.or.ii==n_pt_integ_path2)then
   ene_acc = ene_acc + 0.5d0*func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ene_acc3 = ene_acc3 + 0.5d0*func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
  else
   ene_acc = ene_acc + func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ene_acc3 = ene_acc3 + func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
  end if

 end do path2

 

!####################################################################
![xmin,0]
 if(xmin/=dcmplx(0.d0,0.d0))then
  path3:  do ii = 1,n_pt_integ !the integrand is 0 at 0
   z_path = dcmplx(real(ii,dp)*xmin/real(n_pt_integ,dp),0.d0)
   dz_path = dcmplx(xmin/real(n_pt_integ,dp),0.d0)
   
   Nold = dcmplx(0.d0,0.d0)
   Dold = dcmplx(1.d0,0.d0)
   N = dcmplx(1.d0,0.d0)
   D = z_path - dcmplx(an(0),0.d0)

   do kk=1,nrec
    Nnew = (z_path - dcmplx(an(kk),0.d0))*N - dcmplx(bn2(kk),0.d0)*Nold
    Dnew = (z_path - dcmplx(an(kk),0.d0))*D - dcmplx(bn2(kk),0.d0)*Dold
    
    Nold = N
    Dold = D
    N = Nnew
    D = Dnew
    
    if(kk/=nrec)then
     if((bn2(kk+1)<tol14))exit
    end if
    
   end do
   
   if(ii==n_pt_integ)then
    ene_acc = ene_acc + 0.5d0*func(z_path,twotrotter) *&
&    N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
    ene_acc4 = ene_acc4 + 0.5d0*func(z_path,twotrotter) *&
&    N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   else
    ene_acc = ene_acc + func(z_path,twotrotter) *&
&    N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
    ene_acc4 = ene_acc4 + func(z_path,twotrotter) *&
&    N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   end if
  end do path3
 end if
 
!####################################################################
![xmax,xmax+i*epsilon]
 path4:  do ii = 0,n_pt_integ_path2
  z_path = dcmplx(xmax,real(ii,dp)*epsilon/real(n_pt_integ_path2,dp))
  dz_path = dcmplx(0.d0,epsilon/real(n_pt_integ_path2,dp))
  
  Nold = dcmplx(0.d0,0.d0)
  Dold = dcmplx(1.d0,0.d0)
  N = dcmplx(1.d0,0.d0)
  D = z_path - dcmplx(an(0),0.d0)
  
  do kk=1,nrec
   Nnew = (z_path - dcmplx(an(kk),0.d0))*N - dcmplx(bn2(kk),0.d0)*Nold
   Dnew = (z_path - dcmplx(an(kk),0.d0))*D - dcmplx(bn2(kk),0.d0)*Dold
   
   Nold = N
   Dold = D
   N = Nnew
   D = Dnew
   
   if(kk/=nrec)then
    if((bn2(kk+1)<tol14))exit
   end if
   
  end do
  
  if(ii==0.or.ii==n_pt_integ_path2)then
   ene_acc = ene_acc + 0.5d0*func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ene_acc2 = ene_acc2 + 0.5d0*func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
  else
   ene_acc = ene_acc + func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ene_acc2 = ene_acc2 + func(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
  end if
 end do path4
 
 ene_out = mult*real(1/dcmplx(0.d0,pi)*ene_acc,dp)
 ene_out1 = mult*real(1/dcmplx(0.d0,pi)*ene_acc1,dp)
 ene_out2 = mult*real(1/dcmplx(0.d0,pi)*ene_acc2,dp)
 ene_out3 = mult*real(1/dcmplx(0.d0,pi)*ene_acc3,dp)
 ene_out4 = mult*real(1/dcmplx(0.d0,pi)*ene_acc4,dp)
 
 
end subroutine free_energyrec
!!***
