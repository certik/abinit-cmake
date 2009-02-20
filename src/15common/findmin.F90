!{\src2tex{textfont=tt}}
!!****f* ABINIT/findmin
!!
!! NAME
!! findmin
!!
!! FUNCTION
!! Compute the minimum of a function whose value
!! and derivative are known at two points,
!! using different algorithms.
!! Also deduce different quantities at this predicted
!! point, and at the two other points
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!

!! INPUTS
!! choice=1,uses a linear interpolation of the derivatives
!!       =2,uses a quadratic interpolation based on the
!!        values of the function, and the second derivative at mid-point
!!       =3,uses a cubic interpolation
!!       =4,uses a quartic interpolation, with the supplementary
!!          condition that the second derivative vanishes at one and
!!          only one point (See Schlegel, J. Comp. Chem. 3, 214 (1982).
!!          For this option, lambda_1 must be 1 (new point),
!!          and lambda_2 must be 0 (old point).
!!          Also, if the derivative at the new point is more negative
!!          than the derivative at the old point, the predicted
!!          point cannot correspond to a minimum, but will be lambda=2.5_dp,
!!          if the energy of the second point is lower than the energy
!!          of the first point.
!! etotal_1=first value of the function
!! etotal_2=second value of the function
!! dedv_1=first value of the derivative
!! dedv_2=second value of the derivative
!! lambda_1=first value of the argument
!! lambda_2=second value of the argument
!!
!! OUTPUT
!! dedv_predict=predicted value of the derivative (usually zero,
!!  except if choice=4, if it happens that a minimum cannot be located,
!!  and a trial step is taken)
!! d2edv2_predict=predicted value of the second derivative (not if choice=4)
!! d2edv2_1=first value of the second derivative (not if choice=4)
!! d2edv2_2=second value of the second derivative (not if choice=4)
!! etotal_predict=predicted value of the function
!! lambda_predict=predicted value of the argument
!! status= 0 if everything went normally ;
!!         1 if negative second derivative
!!         2 if some other problem
!!
!! PARENTS
!!      brdene,scfcge
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine findmin(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice
 integer,intent(out) :: status
 real(dp),intent(in) :: dedv_1,dedv_2,etotal_1,etotal_2,lambda_1,lambda_2
 real(dp),intent(out) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_predict
 real(dp),intent(out) :: etotal_predict,lambda_predict

!Local variables-------------------------------
!scalars
 integer :: printvol
 real(dp) :: aa,bb,bbp,cc,ccp,d2edv2_mid,d_lambda,dd,dedv_2bis,dedv_mid1
 real(dp) :: dedv_mid2,discr,ee,eep,etotal_2bis,lambda_shift,sum1,sum2,sum3,uu
 real(dp) :: uu3,vv,vv3
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(6,*)' findmin : enter'
!write(6,*)' choice,lambda_1,lambda_2=',choice,lambda_1,lambda_2
!ENDDEBUG

 status=0
 d_lambda=lambda_1-lambda_2

!DEBUG
!do choice=3,1,-1
!ENDDEBUG

 if(choice==3)then

! Evaluate cubic interpolation
! etotal = aa + bb * lambda + cc * lambda**2 + dd * lambda**3
  dedv_mid1=(dedv_1+dedv_2)/2.0_dp
  dedv_mid2=(etotal_1-etotal_2)/d_lambda
  d2edv2_mid=(dedv_1-dedv_2)/d_lambda
  dd=2.0_dp * ( dedv_mid1 - dedv_mid2 ) / (d_lambda**2)
  cc=0.5_dp * ( d2edv2_mid - 3.0_dp*dd*(lambda_1+lambda_2) )
  bb=dedv_2 - 2*cc*lambda_2 - 3*dd*lambda_2**2
  aa=etotal_2 - bb*lambda_2 - cc*lambda_2**2 - dd*lambda_2**3

! Find the lambda at the minimum
  discr=cc*cc-3*bb*dd
  if(discr<0.0_dp)then
   write(message, '(a,a,a,a)' )ch10,&
&   ' findmin : BUG -',ch10,&
&   '  The 2nd degree equation has no root (choice=3).'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if
  discr=sqrt(discr)
! The root that gives a minimum corresponds to  +discr
  lambda_predict=(-cc+discr)/(3.0_dp*dd)

! Predict etotal at that lambda
  etotal_predict=aa+lambda_predict*(bb&
&  +lambda_predict*(cc&
&  +lambda_predict* dd  ))
  dedv_predict=bb+2.0_dp*cc*lambda_predict&
&  +3.0_dp*dd*lambda_predict**2
  d2edv2_1=2*cc+6*dd*lambda_1
  d2edv2_2=2*cc+6*dd*lambda_2
  d2edv2_predict=2*cc+6*dd*lambda_predict

 else if(choice==4)then

  if(abs(lambda_1-1.0_dp)>tol12 .or. abs(lambda_2)>tol12) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' findmin : BUG -',ch10,&
&   '  For choice=4, lambda_1 must be 1 and lambda_2 must be 0.'
   call wrtout(06,message,'COLL')
   call leave_new('COLL')
  end if

! Evaluate quartic interpolation
! etotal = aa + bb * lambda + cc * lambda**2 + dd * lambda**3 + ee * lambda**4
! Impose positive second derivative everywhere, with
! one point where it vanishes :  3*dd**2=8*cc*ee
  aa=etotal_2
  bb=dedv_2
  sum1=etotal_1-aa-bb
  sum2=dedv_1-bb
  sum3=sum2-2.0_dp*sum1

! Build the discriminant of the associated 2nd degree equation
  discr=sum2**2-3.0_dp*sum3**2
  if(discr<0.0_dp .or. sum2<0.0_dp)then

!  Even if there is a problem, try to keep going ...
   write(message, '(a,a,a,a)' )ch10,&
&   ' findmin : WARNING -',ch10,&
&   '  The 2nd degree equation has no positive root (choice=4).'
   status=2
   call wrtout(06,message,'COLL')
   if(etotal_1<etotal_2)then
    write(message, '(a,a,a,a,a,a)' )ch10,&
&    ' findmin : COMMENT -',ch10,&
&    '  Will continue, since the new total energy is lower',ch10,&
&    '  than the old. Take a larger step in the same direction.'
    call wrtout(06,message,'COLL')
    lambda_predict=2.5_dp
   else
    write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&    ' findmin : COMMENT -',ch10,&
&    '  There is a problem, since the new total energy is larger',ch10,&
&    '  than the old (choice=4).',ch10,&
&    '  I take a point between the old and new, close to the old .'
    call wrtout(06,message,'COLL')
    lambda_predict=0.25_dp
   end if
!  Mimick a zero-gradient lambda, in order to avoid spurious
!  action of the inverse hessian (the next line would be a realistic estimation)
   dedv_predict=0.0_dp
!  dedv_predict=dedv_2+lambda_predict*(dedv_1-dedv_2)
!  Uses the energies, and the gradient at lambda_2
   etotal_predict=etotal_2+dedv_2*lambda_predict&
&   +(etotal_1-etotal_2-dedv_2)*lambda_predict**2

  else

!  Here, there is an acceptable solution to the 2nd degree equation
   discr=sqrt(discr)
!  The root that gives the smallest ee corresponds to  -discr
!  This is the one to be used: one aims at modelling the
!  behaviour of the function as much as possible with the
!  lowest orders of the polynomial, not the quartic term.
   ee=(sum2-discr)*0.5_dp
   dd=sum3-2.0_dp*ee
   cc=sum1-dd-ee

!  DEBUG
!  write(6,*)'aa,bb,cc,dd,ee',aa,bb,cc,dd,ee
!  ENDDEBUG

!  Now, must find the unique root of
!  $0 = bb + 2*cc * lambda + 3*dd * lambda^2 + 4*ee * lambda^3$
!  This root is unique because it was imposed that the second derivative
!  of the quartic polynomial is everywhere positive.
!  First, remove the quadratic term, by a shift of lambda
!  lambdap=lambda-lambda_shift
!  $0 = bbp + ccp * lambdap + eep * lambdap^3$
   eep=4.0_dp*ee
   lambda_shift=-dd/(4.0_dp*ee)
   ccp=2.0_dp*cc-12.0_dp*ee*lambda_shift**2
   bbp=bb+ccp*lambda_shift+eep*lambda_shift**3

!  DEBUG
!  write(6,*)'bbp,ccp,eep,lambda_shift',bbp,ccp,eep,lambda_shift
!  ENDDEBUG

!  The solution of a cubic polynomial equation is as follows :
   discr=(bbp/eep)**2+(4.0_dp/27.0_dp)*(ccp/eep)**3
!  In the present case, discr will always be positive
   discr=sqrt(discr)
   uu3=0.5_dp*(-bbp/eep+discr) ; uu=sign((abs(uu3))**(1.0_dp/3.0_dp),uu3)
   vv3=0.5_dp*(-bbp/eep-discr) ; vv=sign((abs(vv3))**(1.0_dp/3.0_dp),vv3)
   lambda_predict=uu+vv

!  Restore the shift
   lambda_predict=lambda_predict+lambda_shift
   etotal_predict=aa+bb*lambda_predict+cc*lambda_predict**2+&
&   dd*lambda_predict**3+ee*lambda_predict**4
   dedv_predict=bb+2.0_dp*cc*lambda_predict+3.0_dp*dd*lambda_predict**2+&
&   4.0_dp*ee*lambda_predict**3
   d2edv2_1=2*cc+6*dd*lambda_1+12*ee*lambda_1**2
   d2edv2_2=2*cc+6*dd*lambda_2+12*ee*lambda_2**2
   d2edv2_predict=2*cc+6*dd*lambda_predict+12*ee*lambda_predict**2

  end if

 else if(choice==1) then

! Use the derivative information to predict lambda
  d2edv2_mid=(dedv_1-dedv_2)/d_lambda
  lambda_predict=lambda_2-dedv_2/d2edv2_mid
  dedv_predict=dedv_2+(lambda_predict-lambda_2)*d2edv2_mid
  d2edv2_1=d2edv2_mid
  d2edv2_2=d2edv2_mid
  d2edv2_predict=d2edv2_mid
! also use the first energy to predict new energy
  etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&  +0.5_dp*d2edv2_1*(lambda_predict-lambda_1)**2
  etotal_2bis=etotal_1+dedv_1*(lambda_2-lambda_1)&
&  +0.5_dp*d2edv2_1*(lambda_2-lambda_1)**2

  if(d2edv2_mid<0.0_dp)then
   write(message, '(a,a,a,a,es18.10,a)' ) ch10,&
&   ' findmin : WARNING -',ch10,&
&   '  (scfcge) The second derivative is negative, equal to',d2edv2_mid        ,'.'
   call wrtout(6,message,'COLL')
   status=1
  end if

 else if(choice==2) then

! Use energies and first derivative information
! etotal = aa + bb * lambda + cc * lambda**2
  dedv_mid2=(etotal_1-etotal_2)/d_lambda
  cc=(dedv_1-dedv_mid2)/d_lambda
  lambda_predict=lambda_1-0.5_dp*dedv_1/cc
  d2edv2_1=2*cc
  d2edv2_2=d2edv2_1
  d2edv2_predict=d2edv2_1
  if(d2edv2_predict<0.0_dp)then
   write(message, '(a,a,a,a,es18.10,a,a,a)' ) ch10,&
&   ' findmin : WARNING -',ch10,&
&   '  (scfcge) The second derivative is negative, equal to',d2edv2_predict,'.',&
&   ch10,'  (scfcge) => Pivoting                     '
   call wrtout(6,message,'COLL')
   status=1
   if(etotal_2 < etotal_1)then
    lambda_predict=lambda_2-0.5_dp*(lambda_1-lambda_2)
   else
    lambda_predict=lambda_1-0.5_dp*(lambda_2-lambda_1)
   end if
  end if
  dedv_predict=dedv_1+(lambda_predict-lambda_1)*d2edv2_1
  dedv_2bis=dedv_1+(lambda_2-lambda_1)*d2edv2_1
  etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&  +0.5_dp*d2edv2_1*(lambda_predict-lambda_1)**2

 end if
 printvol=1
 if(choice==4)printvol=2
 if(printvol==2)then
  write(message, '(a,i3)' )'   line minimization, algorithm ',choice
  call wrtout(6,message,'COLL')
  write(message, '(a,a)' )'                        lambda      etotal ',&
&  '           dedv        d2edv2    '
  call wrtout(6,message,'COLL')
  write(message, '(a,es12.4,es18.10,2es12.4)' )&
&  '   old point         :',lambda_2,etotal_2,dedv_2,d2edv2_2
  call wrtout(6,message,'COLL')
  write(message, '(a,es12.4,es18.10,2es12.4)' )&
&  '   new point         :',lambda_1,etotal_1,dedv_1,d2edv2_1
  call wrtout(6,message,'COLL')
  write(message, '(a,es12.4,es18.10,2es12.4)' )&
&  '   predicted point   :',&
&  lambda_predict,etotal_predict,dedv_predict,d2edv2_predict
  call wrtout(6,message,'COLL')
  if(choice==1) then
   write(message, '(a,es10.4)' ) &
&   ' consistency check :    etotal_2 =',etotal_2bis
   call wrtout(6,message,'COLL')
  end if
  if(choice==2) then
   write(message, '(a,es10.4)' ) &
&   ' consistency check :    dedv_2 =',dedv_2bis
   call wrtout(6,message,'COLL')
  end if
  write(message, '(a)' ) ' '
  call wrtout(6,message,'COLL')
 else if(printvol==1)then
  write(message, '(a,es12.4,a,es18.10)' ) &
&  ' findmin : lambda_predict ',lambda_predict,&
&  ' etotal_predict ',etotal_predict
  call wrtout(6,message,'COLL')
 end if

end subroutine findmin
!!***
