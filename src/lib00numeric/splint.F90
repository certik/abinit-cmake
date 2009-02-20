 subroutine splint (nspline,xspline,yspline,ysplin2,&
&                  nfit,xfit,yfit)

!  Compute spline interpolation. There is no hypothesis
!  about the spacing of the input grid points.

!  ON INPUT:
!  nspline: number of grid points of input mesh
!  xspline(nspline): input mesh
!  yspline(nspline): function on input mesh
!  ysplin2(nspline): second derivative of yspline on input mesh
!  nfit: number of points of output mesh
!  xfit(nfit): output mesh
!  ON OUTPUT:
!  yfit(nfit): function on output mesh

 implicit none
 integer left,i,nfit,nspline,k,right
 double precision xspline(nspline), yspline(nspline), &
& ysplin2(nspline),xfit(nfit),yfit(nfit),delarg,invdelarg,aa,bb

 left = 1
 do i=1, nfit

! Initialize for the unlikely event that rmax exceed r(mesh)
  yfit(i)=0.d0

  do k=left+1, nspline
   if(xspline(k) >= xfit(i)) then
    if(xspline(k-1) <= xfit(i)) then
     right = k
     left = k-1
    else
     if (k-1.eq.1 .and. i.eq.1) then
      stop '  splint: xfit(1) < xspline(1)'
     else
      stop '  splint: xfit not properly ordered'
     end if
    end if
    delarg= xspline(right) - xspline(left)
    invdelarg= 1.0d0/delarg
    aa= (xspline(right)-xfit(i))*invdelarg
    bb= (xfit(i)-xspline(left))*invdelarg

    yfit(i) = aa*yspline(left) + bb*yspline(right)    &
&            +( (aa*aa*aa-aa)*ysplin2(left) +         &
&               (bb*bb*bb-bb)*ysplin2(right) ) *delarg*delarg/6.0d0
    exit
   end if
  end do ! k
 end do ! i

 return
 end


