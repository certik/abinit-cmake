 subroutine splfit(arg,derfun,fun,ider,newarg,newfun,numarg,numnew)

! Evaluate cubic spline fit to get function values on input set
! of ORDERED, UNFORMLY SPACED points.
! Optionally gives derivatives (first and second) at those points too.
! If point lies outside the range of arg, assign the extremal
! point values to these points, and zero derivative.
! Note : if ider=0, compute only the function (contained in fun)
!        if ider=1, compute the function (contained in fun) and its first derivative (in derfun)
!        if ider=2, compute only the second derivative of the function (in derfun)
!
! Input: 
!  arg(numarg)=equally spaced arguments (spacing delarg) for data 
!   to which spline was fit.
!  fun(numarg,2)=function values to which spline was fit and spline
!   fit to second derivatives (from Numerical Recipes spline).
!  ider=  see above
!  newarg(numnew)=new values of arguments at which function is desired.
!  numarg=number of arguments at which spline was fit.
!  numnew=number of arguments at which function values are desired.
! Output:
!  derfun(numnew)=(optional) values of first or second derivative of function.
!   This is only computed for ider=1 or 2; otherwise derfun not used.
!  newfun(numnew)=values of function at newarg(numnew).
!   This is only computed for ider=0 or 1.

 implicit none
 integer ider,numarg,numnew
 double precision arg(numarg),derfun(numnew),fun(numarg,2),newarg(numnew),newfun(numnew)
 integer i,jspl
 double precision argmin,delarg,d,aa,bb,cc,dd
 character message*500

!argmin is smallest x value in spline fit; delarg is uniform spacing of spline argument
 argmin=arg(1)
 delarg=(arg(numarg)-argmin)/dble(numarg-1)

!Do one loop for no grads, other for grads:
 if (ider==0) then
! Spline index loop for no grads:
  do i=1,numnew
   if (newarg(i).ge.arg(numarg)) then
    write(message,1000)char(10),i,newarg(i), &
&     jspl,char(10),numarg,arg(numarg),char(10),char(10),char(10)
1000 format(a1,' splfit: for arg number',i8,2x,'of value', &
&    1p,e12.4,1x,'jspl=',i8,a1,' is >= numarg=',i8,  &
&    '  (max arg(numarg)=',e12.4,')',a1,             &
&    ' Means function values are being requested outside',       &
&    ' range of data.',a1,' Function and slope will be set to',  &
&    ' values at upper end of data.',a1)

    newfun(i)=fun(numarg,1)

   else if (newarg(i).le.arg(1)) then
    newfun(i)=fun(1,1)

   else
    jspl=1+int((newarg(i)-argmin)/delarg)
    d=newarg(i)-arg(jspl)
    bb = d/delarg
    aa = 1.0d0-bb
    cc = aa*(aa**2-1.0d0)*(delarg**2/6.0d0)
    dd = bb*(bb**2-1.0d0)*(delarg**2/6.0d0)
    newfun(i)=aa*fun(jspl,1)+bb*fun(jspl+1,1)+cc*fun(jspl,2)+dd*fun(jspl+1,2)
   end if
  enddo

 else if(ider==1)then

! Spline index loop includes grads:
  do i=1,numnew

   if (newarg(i).ge.arg(numarg)) then
    newfun(i)=fun(numarg,1)
    derfun(i)=0.0d0

   else if (newarg(i).le.arg(1)) then
    newfun(i)=fun(1,1)
    derfun(i)=0.0d0

   else

!   cubic spline interpolation:
    jspl=1+int((newarg(i)-argmin)/delarg)
    d=newarg(i)-arg(jspl)
    bb = d/delarg
    aa = 1.0d0-bb
    cc = aa*(aa**2-1.0d0)*(delarg**2/6.0d0)
    dd = bb*(bb**2-1.0d0)*(delarg**2/6.0d0)
    newfun(i)=aa*fun(jspl,1)+bb*fun(jspl+1,1)+cc*fun(jspl,2)+dd*fun(jspl+1,2)
!   spline fit to first derivative:
!   note correction of Numerical Recipes sign error
    derfun(i) = (fun(jspl+1,1)-fun(jspl,1))/delarg +    &
&      (-(3.d0*aa**2-1.d0)*fun(jspl,2)+                 &
&        (3.d0*bb**2-1.d0)*fun(jspl+1,2)) * delarg/6.0d0

          end if
  enddo

 else if (ider==2) then

  do i=1,numnew

   if (newarg(i).ge.arg(numarg)) then
    derfun(i)=0.0d0

   else if (newarg(i).le.arg(1)) then
    derfun(i)=0.0d0

   else

!   cubic spline interpolation:
    jspl=1+int((newarg(i)-argmin)/delarg)
    d=newarg(i)-arg(jspl)
    bb = d/delarg
    aa = 1.0d0-bb
!   second derivative of spline (piecewise linear function)
    derfun(i) = aa*fun(jspl,2)+bb*fun(jspl+1,2)

   end if
  enddo

 end if

 return

 end
