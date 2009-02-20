subroutine spline( t, y, n, ybcbeg, ybcend, ypp, tmp )
!
!*******************************************************************************
!
!! SPLINE (originally SPLINE_CUBIC_SET) computes the second derivatives of a cubic spline.
!
!
!  Discussion:
!
!    For data interpolation, the user must call SPLINE_CUBIC_SET to 
!    determine the second derivative data, passing in the data to be 
!    interpolated, and the desired boundary conditions.
!
!    The data to be interpolated, plus the SPLINE_CUBIC_SET output, 
!    defines the spline.  The user may then call SPLINE_CUBIC_VAL to 
!    evaluate the spline at any point.
!
!    The cubic spline is a piecewise cubic polynomial.  The intervals
!    are determined by the "knots" or abscissas of the data to be
!    interpolated.  The cubic spline has continous first and second
!    derivatives over the entire interval of interpolation.  
!
!    For any point T in the interval T(IVAL), T(IVAL+1), the form of
!    the spline is
!
!      SPL(T) = A(IVAL)
!             + B(IVAL) * ( T - T(IVAL) ) 
!             + C(IVAL) * ( T - T(IVAL) )**2
!             + D(IVAL) * ( T - T(IVAL) )**3
!
!    If we assume that we know the values Y(*) and YPP(*), which represent
!    the values and second derivatives of the spline at each knot, then
!    the coefficients can be computed as:
!
!      A(IVAL) = Y(IVAL)
!      B(IVAL) = ( Y(IVAL+1) - Y(IVAL) ) / ( T(IVAL+1) - T(IVAL) )
!        - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * ( T(IVAL+1) - T(IVAL) ) / 6
!      C(IVAL) = YPP(IVAL) / 2
!      D(IVAL) = ( YPP(IVAL+1) - YPP(IVAL) ) / ( 6 * ( T(IVAL+1) - T(IVAL) ) )
!
!    Since the first derivative of the spline is
!
!      SPL'(T) =     B(IVAL)
!              + 2 * C(IVAL) * ( T - T(IVAL) )
!              + 3 * D(IVAL) * ( T - T(IVAL) )**2,
!
!    the requirement that the first derivative be continuous at interior
!    knot I results in a total of N-2 equations, of the form:
!
!      B(IVAL-1) + 2 C(IVAL-1) * (T(IVAL)-T(IVAL-1)) 
!      + 3 * D(IVAL-1) * (T(IVAL) - T(IVAL-1))**2 = B(IVAL)
!
!    or, setting H(IVAL) = T(IVAL+1) - T(IVAL)
!
!      ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)
!      - ( YPP(IVAL) + 2 * YPP(IVAL-1) ) * H(IVAL-1) / 6
!      + YPP(IVAL-1) * H(IVAL-1)
!      + ( YPP(IVAL) - YPP(IVAL-1) ) * H(IVAL-1) / 2
!      = 
!      ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - ( YPP(IVAL+1) + 2 * YPP(IVAL) ) * H(IVAL) / 6
!
!    or
!
!      YPP(IVAL-1) * H(IVAL-1) + 2 * YPP(IVAL) * ( H(IVAL-1) + H(IVAL) )
!      + YPP(IVAL) * H(IVAL) 
!      =
!      6 * ( Y(IVAL+1) - Y(IVAL) ) / H(IVAL)
!      - 6 * ( Y(IVAL) - Y(IVAL-1) ) / H(IVAL-1)    
!
!    Boundary conditions must be applied at the first and last knots.  
!    The resulting tridiagonal system can be solved for the YPP values.
!
!  Modified:
!
!    07 February 1999
!    28 November 2004 XGonze : double precision
!                              make arguments similar to the Numeric Recipes routine
!                              also use algorithmics similar to the Numeric Recipes routine
!
!  Author:
!
!    John Burkardt
!    (XGonze got it from http://www.psc.edu/~burkardt/src/spline/spline.html)
!
!  Parameters:
!
!    Input, integer N, the number of data points; N must be at least 2. 
!    In the special case where N = 2 and IBCBEG = IBCEND = 0, the 
!    spline will actually be linear. 
!
!    Input, double precision T(N), the knot values, that is, the points where data
!    is specified.  The knot values should be distinct, and increasing.
!
!    Input, double precision Y(N), the data values to be interpolated.
!
!    Input, double precision YBCBEG, YBCEND, the values to be used in the boundary
!    conditions if IBCBEG or IBCEND is equal to 1 or 2.
!
!    Output, double precision YPP(N), the second derivatives of the cubic spline.
!
!    Work space, double precision DIAG(N) - should be removed ...
!
!
!    XG041127 : In the initial implementation, one had the control on
!     IBCBEG and IBCEND. Now, they are determined by the values
!     of YBCBEG, YBCEND. Option 2 has been disabled.
!
!    Input, integer IBCBEG, left boundary condition flag:
!
!      0: the spline should be a quadratic over the first interval;
!      1: the first derivative at the left endpoint should be YBCBEG;
!      2: the second derivative at the left endpoint should be YBCBEG.
!
!    Input, integer IBCEND, right boundary condition flag:
!
!      0: the spline should be a quadratic over the last interval;
!      1: the first derivative at the right endpoint should be YBCEND;
!      2: the second derivative at the right endpoint should be YBCEND.

  implicit none
  integer n
!
  integer ibcbeg
  integer ibcend
  double precision t(n)
  double precision y(n)
  double precision ybcbeg
  double precision ybcend
  double precision ypp(n)
  double precision tmp(n)

  integer i,k
  double precision ratio,pinv
!
!  Check.
!
  if ( n <= 1 ) then
    write ( *, * ) ' '
    write ( *, * ) 'SPLINE_CUBIC_SET - Fatal error!'
    write ( *, * ) '  The number of knots must be at least 2.'
    write ( *, * ) '  The input value of N = ', n
    stop
  end if

  do i = 1, n-1
    if ( t(i) >= t(i+1) ) then
      write ( *, * ) ' '
      write ( *, * ) 'SPLINE_CUBIC_SET - Fatal error!'
      write ( *, * ) '  The knots must be strictly increasing, but'
      write ( *, * ) '  T(',  i,') = ', t(i)
      write ( *, * ) '  T(',i+1,') = ', t(i+1)
      stop
    end if
  end do
!
!  XG041127
  ibcbeg=1 ; ibcend=1
  if(ybcbeg>1.0d+30)ibcbeg=0
  if(ybcend>1.0d+30)ibcend=0
!
!  Set the first and last equations.
!
  if ( ibcbeg == 0 ) then
    ypp(1) = 0.d0
    tmp(1) = 0.d0
  else if ( ibcbeg == 1 ) then
    ypp(1) = -0.5d0
    tmp(1) = (3.d0/(t(2)-t(1)))*((y(2)-y(1))/(t(2)-t(1))-ybcbeg)
  end if
  if ( ibcend == 0 ) then
    ypp(n) = 0.d0
    tmp(n) = 0.d0
  else if ( ibcend == 1 ) then
    ypp(n) = 0.5d0
    tmp(n) = (3.d0/(t(n)-t(n-1)))*(ybcend-(y(n)-y(n-1))/(t(n)-t(n-1)))
  end if

!
!  Set the intermediate equations.
!
  do i=2,n-1
   ratio=(t(i)-t(i-1))/(t(i+1)-t(i-1))
   pinv = 1.0d0/(ratio*ypp(i-1) + 2.0d0)
   ypp(i) = (ratio-1.0d0)*pinv
   tmp(i)=(6.0d0*((y(i+1)-y(i))/(t(i+1)-t(i))-(y(i)-y(i-1)) &
&    /(t(i)-t(i-1)))/(t(i+1)-t(i-1))-ratio*tmp(i-1))*pinv
   if (abs(tmp(i))<1.d5*tiny(0.d0)) tmp(i)=0.d0   !MT20050927
  enddo

! Solve the equations
  ypp(n) = (tmp(n)-ypp(n)*tmp(n-1))/(ypp(n)*ypp(n-1)+1.0d0)
  do k=n-1,1,-1
   ypp(k)=ypp(k)*ypp(k+1)+tmp(k)
  enddo

  return
end
