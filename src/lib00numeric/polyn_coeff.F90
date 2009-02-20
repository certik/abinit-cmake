 subroutine polyn_coeff(n,x,y,coeff)
!
! For N function values Y(X) compute coefficients of
! N-1 degree interpolating polynomial. 
! Due to G. Rybicki
!
! Input:
! n = number of points 
! x = array of abcissa values
! y = array of ordinate values
!
! Output:
! coeff(n) = array of polynomial coefficients
!
!         tested for linear and parabolic interpolation
!                             

 implicit none
!arguments
 integer :: n
 double precision ::  coeff(n),x(n),y(n)
!local variables
 integer :: ii,jj,kk
 double precision :: acc,ff,den
 double precision, allocatable :: s(:)

 allocate (s(n))
 s(:)=0.0d0
 coeff(:)=0.0d0
 s(n)=-x(1)

 do ii=2,n
  do jj=n+1-ii,n-1
   s(jj)=s(jj)-x(ii)*s(jj+1)
  enddo
  s(n)=s(n)-x(ii)
 enddo

 do jj=1,n
  den=n
  do kk=n-1,1,-1
   den=kk*s(kk+1)+x(jj)*den
  enddo
  ff=y(jj)/den
  acc=1.0d0
  do kk=n,1,-1
   coeff(kk)=coeff(kk)+acc*ff
   acc=s(kk)+x(jj)*acc
  enddo
 enddo

 deallocate(s)

 end subroutine
