!{\src2tex{textfont=tt}}
!!****f* ABINIT/entropyrec
!! NAME
!! entropyrec
!! 
!! FUNCTION
!! This routine computes the local part of the entropy at a point using a path integral, 
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
!!  mult=a multiplicator for computing entropy ; 2 for non-spin-polarized system
!!  prtvol=printing volume
!! 
!! OUTPUT
!!  ent_out=entropy at the point
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

subroutine entropyrec(an,bn2,nrec,trotter,ent_out, mult, &
&                     prtvol,n_pt_integ,xmax,&
&                     ent_out1,ent_out2,ent_out3,ent_out4)

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
 real(dp),intent(out) :: ent_out,ent_out1,ent_out2,ent_out3,ent_out4
!arrays
 real(dp),intent(in) :: an(0:nrec),bn2(0:nrec)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=7
 integer,save :: first=1
 integer :: ii,jj,kk,ll,mult_n_pt_integ,n_pt_integ_path2
 real(dp) :: arg,arg_path,epsilon,step,teta,twotrotter,xmin,ymin
 complex(dpc) :: D,Dnew,Dold,N,Nnew,Nold,dz_path,ent_acc,ent_acc1,ent_acc2
 complex(dpc) :: ent_acc3,ent_acc4,ent_acc_path3,ent_acc_path56,ent_acc_test1
 complex(dpc) :: ent_acc_test2,ent_acc_test3,ent_acc_test4,func1,func2
 complex(dpc) :: rootof_z_path,rootof_z_path2,z_path,zj
 character(len=500) :: message
!arrays
 real(dp) :: tsec(2)

! *************************************************************************
 
!function to integrate over the path
 func1(z_path,twotrotter) =  ( z_path**twotrotter/(1+z_path**twotrotter)*log(1+1/z_path**twotrotter)+&    !- f*ln(f)
&1/(1+z_path**twotrotter)*log(1+z_path**twotrotter))       !- (1-f)*ln(1-f)
!other expression of func for a path like ro(t)*exp(2*i*pi/(2*p)*(j+1/2))
!!$ func2(ro,twotrotter) = 

 call timab(190,1,tsec)
 
!structured debugging if prtvol=-level : print detailled result the first time we enter entropyrec
 if(prtvol==-level)then
  if(first==1)then
   write(message,'(a)')' ' 
   call wrtout(06,message,'PERS')
   write(message,'(a)')' entropyrec : enter ' 
   call wrtout(06,message,'PERS')

   write(message,'(a,i6)')'n_pt_integ ' , n_pt_integ
   call wrtout(06,message,'COLL')

  end if
 end if
 
 ent_out = 0.d0
 ent_out1 = 0.d0
 ent_out2 = 0.d0
 ent_out3 = 0.d0
 ent_out4 = 0.d0
 ent_acc = dcmplx(0.d0,0.d0)
 ent_acc1 = dcmplx(0.d0,0.d0)
 ent_acc2 = dcmplx(0.d0,0.d0)
 ent_acc3 = dcmplx(0.d0,0.d0)
 ent_acc4 = dcmplx(0.d0,0.d0)
 
!path parameters 
 xmin = 2.5d-1

 mult_n_pt_integ = 1 !25
 epsilon = 2.5d-1


 ent_out = 0.d0
 ent_acc = dcmplx(0.d0,0.d0)
 ent_acc_path3 = dcmplx(0.d0,0.d0)
 ent_acc_path56 = dcmplx(0.d0,0.d0)
 
!path parameters 
!n_pt_integ = 2500
 mult_n_pt_integ = 1 !25
 step = (xmax-xmin)/real(n_pt_integ,dp)
 if(trotter==0)then
  twotrotter = 1.d0
 else
  twotrotter = 2.d0*real(trotter,dp)
 end if
!xmin = 0.d0 !-0.5d-1**(1.d0/twotrotter)
 arg = pi/twotrotter
 zj = dcmplx( cos(arg) , sin(arg) )
 epsilon = 1.d0/2.d0*sin( arg )
 xmin = 1.d0/2.d0*cos( arg )
 
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
   ent_acc = ent_acc + 0.5d0*func1(z_path,twotrotter)*&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ent_acc1 = ent_acc1 + 0.5d0*func1(z_path,twotrotter)*&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
  else
   ent_acc = ent_acc + func1(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ent_acc1 = ent_acc1 + func1(z_path,twotrotter) *&
&   N/D *dz_path       
  end if
 end do path1
 
!!$ !####################################################################
!!$ ![xmin + i*epsilon,xmin]
!!$ if(epsilon/step>4.d0)then
!!$  n_pt_integ_path2 = int(epsilon/step)+1
!!$ else
!!$  n_pt_integ_path2 = 5
!!$ endif
!!$n_pt_integ_path2 = n_pt_integ
!!$ path2:  do ii = 0,n_pt_integ_path2
!!$  z_path = dcmplx(xmin,real(ii,dp)*epsilon/real(n_pt_integ_path2,dp))
!!$  dz_path = -dcmplx(0.d0,epsilon/real(n_pt_integ_path2,dp))
!!$
!!$if(z_path/=0.d0)then
!!$  
!!$  Nold = dcmplx(0.d0,0.d0)
!!$  Dold = dcmplx(1.d0,0.d0)
!!$  N = dcmplx(1.d0,0.d0)
!!$  D = z_path - dcmplx(an(0),0.d0)
!!$  
!!$  do kk=1,nrec
!!$   Nnew = (z_path - dcmplx(an(kk),0.d0))*N - dcmplx(bn2(kk),0.d0)*Nold
!!$   Dnew = (z_path - dcmplx(an(kk),0.d0))*D - dcmplx(bn2(kk),0.d0)*Dold
!!$   
!!$   Nold = N
!!$   Dold = D
!!$   N = Nnew
!!$   D = Dnew
!!$    
!!$   if(kk/=nrec)then
!!$    if((bn2(kk+1)<tol14))exit
!!$   end if
!!$   
!!$  enddo
!!$  
!!$  if(ii==0.or.ii==n_pt_integ_path2)then
!!$   ent_acc = ent_acc - 0.5d0*func1(z_path,twotrotter) *&
!!$&         N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
!!$!DEBUG
!!$   ent_acc3 = ent_acc3 - 0.5d0*func1(z_path,twotrotter) *&
!!$&         N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
!!$write(6,*)' min extremity ', func1(z_path,twotrotter)*N/D, epsilon 
!!$write(6,*)' min integral', func1(z_path,twotrotter)*N/D*epsilon 
!!$!ENDDEBUG
!!$  else
!!$   ent_acc = ent_acc - func1(z_path,twotrotter) *&
!!$&         N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
!!$!DEBUG
!!$   ent_acc3 = ent_acc3 - func1(z_path,twotrotter) *&
!!$&         N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
!!$!ENDDEBUG
!!$  endif
!!$
!!$end if
!!$
!!$ enddo path2



!!$write(6,*)' minpath', ent_acc3


!####################################################################
![1/2zj,0]
 if(epsilon/step>4.d0)then
  n_pt_integ_path2 = int(epsilon/step)+1
 else
  n_pt_integ_path2 = 5
 end if
 n_pt_integ_path2 = n_pt_integ
 path5:  do ii = 0,n_pt_integ_path2

  z_path = dcmplx(real(ii,dp)/(2.d0*real(n_pt_integ_path2,dp))*cos(arg),&
&  real(ii,dp)/(2.d0*real(n_pt_integ_path2,dp))*sin(arg))
  dz_path = -dcmplx(1.d0/(2.d0*real(n_pt_integ_path2,dp))*cos(arg),&
&  1.d0/(2.d0*real(n_pt_integ_path2,dp))*sin(arg))

  if(abs(z_path)>tol14)then
   
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
    ent_acc = ent_acc + 0.5d0*func1(z_path,twotrotter) *&
&    N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
    ent_acc3 = ent_acc3 + 0.5d0*func1(z_path,twotrotter) *&
&    N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   else
    ent_acc = ent_acc + func1(z_path,twotrotter) *&
&    N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
    ent_acc3 = ent_acc3 + func1(z_path,twotrotter) *&
&    N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   end if

  end if

 end do path5

!!$ !####################################################################
!!$ ![xmin,0]
!!$if(xmin/=0.d0)then
!!$ path3:  do ii = 1,n_pt_integ !the integrand is 0 at 0
!!$  z_path = dcmplx(real(ii,dp)*xmin/real(n_pt_integ,dp),0.d0)
!!$  dz_path = dcmplx(xmin/real(n_pt_integ,dp),0.d0)
!!$  
!!$  Nold = dcmplx(0.d0,0.d0)
!!$  Dold = dcmplx(1.d0,0.d0)
!!$  N = dcmplx(1.d0,0.d0)
!!$  D = z_path - dcmplx(an(0),0.d0)
!!$
!!$  do kk=1,nrec
!!$   Nnew = (z_path - dcmplx(an(kk),0.d0))*N - dcmplx(bn2(kk),0.d0)*Nold
!!$   Dnew = (z_path - dcmplx(an(kk),0.d0))*D - dcmplx(bn2(kk),0.d0)*Dold
!!$   
!!$   Nold = N
!!$   Dold = D
!!$   N = Nnew
!!$   D = Dnew
!!$   
!!$   if(kk/=nrec)then
!!$    if((bn2(kk+1)<tol14))exit
!!$   end if
!!$   
!!$  enddo
!!$  
!!$  if(ii==n_pt_integ)then
!!$   ent_acc = ent_acc - 0.5d0*func1(z_path,twotrotter) *&
!!$&         N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
!!$!DEBUG
!!$   ent_acc4 = ent_acc4 - 0.5d0*func1(z_path,twotrotter) *&
!!$&         N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
!!$!ENDDEBUG
!!$  else
!!$   ent_acc = ent_acc - func1(z_path,twotrotter) *&
!!$&         N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
!!$!DEBUG
!!$   ent_acc4 = ent_acc4 - func1(z_path,twotrotter) *&
!!$&         N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
!!$!ENDDEBUG
!!$  endif
!!$ enddo path3
!!$end if
 
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
   ent_acc = ent_acc + 0.5d0*func1(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ent_acc2 = ent_acc2 + 0.5d0*func1(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
  else
   ent_acc = ent_acc + func1(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
   ent_acc2 = ent_acc2 + func1(z_path,twotrotter) *&
&   N/D *dz_path                                        !<r|1/(z-e**(-beta/(2p)*(H-mu)))|r> dz
  end if
 end do path4

!!$write(6,*)' maxpath', ent_acc2
 
 ent_out = mult*real(1/dcmplx(0.d0,pi)*ent_acc,dp)
 ent_out1 = mult*real(1/dcmplx(0.d0,pi)*ent_acc1,dp)
 ent_out2 = mult*real(1/dcmplx(0.d0,pi)*ent_acc2,dp)
 ent_out3 = mult*real(1/dcmplx(0.d0,pi)*ent_acc3,dp)
 ent_out4 = mult*real(1/dcmplx(0.d0,pi)*ent_acc4,dp)
 
 
end subroutine entropyrec
!!***
