!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_numeric_tools
!! NAME
!!  m_numeric_tools
!!
!! FUNCTION
!!  This module contains basic tools for numeric computations.
!!
!! COPYRIGHT
!! Copyright (C) 2008 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_errors.h"

!@ ABIDOC
MODULE m_numeric_tools

 use defs_basis
 use m_errors

 implicit none

! === List of available public routines and functions ===
 public ::       &
&  arth,         &   ! Return an arithmetic progression
&  geop,         &   ! Return a geometric progression
&  set2unit,     &   ! Set the matrix to be a unit matrix (if it is square)
&  trace,        &   ! Calculate the trace of a square matrix
&  get_diag,     &   ! Return the diagonal of a matrix as a vector
&  r2c,c2r,      &   ! Transfer complex data stored in a real array to a complex array and vice versa
&  is_even,      &   ! Return .TRUE. if int is even
&  is_integer,   &   ! Return .TRUE. if all elements of rr differ from an integer by less than tol
&  is_zero,      &   ! Return .TRUE. if all elements of rr differ from zero by less than tol
&  bisect,       &   ! Given a monotonic array A and x find j such that A(j)>x>A(j+1) using bisection
&  imax_loc,     &   ! Index of maxloc on an array returned as scalar instead of array-valued quantity
&  imin_loc,     &   ! Index of minloc on an array returned as scalar instead of array-valued quantity 
&  linfit,       &   ! Perform a linear fit, y=ax+b, of data
&  llsfit_svd,   &   ! Linear least squares fit with SVD and an user-defined set of functions
&  polyn_interp, &   ! Polynomial interpolation with Nevilles"s algorithms, error estimate is reported 
&  quadrature,   &   ! Driver routine to perform quadratures in finite domains using different algorithms
&  hermitianize, &   ! Force a square matrix to be hermitian
&  print_arr,    &   ! Print a vector/array
&  pade,dpade,   &   ! Functions for Pade approximation (complex case)
&  newrap_step,  &   ! Apply single step Newton-Raphson method to find root of a complex function
&  OPERATOR(.x.),&   ! Cross product of two 3D vectors 
&  l2norm            ! Return the length (ordinary L2 norm) of a vector
!@END ABIDOC

 interface arth
  module procedure arth_int 
  module procedure arth_rdp
 end interface

 interface set2unit
  module procedure unit_matrix_int 
  module procedure unit_matrix_rdp
 end interface

 interface trace
  module procedure trace_int
  module procedure trace_rdp
  module procedure trace_cdp
 end interface

 interface get_diag
  module procedure get_diag_rdp
 end interface 

 interface r2c 
  module procedure rdp2cdp_1d 
 end interface 

 interface c2r 
  module procedure cdp2rdp_1d 
 end interface 

 interface is_integer
  module procedure is_integer_0d 
  module procedure is_integer_1d
 end interface 

 interface is_zero
  module procedure is_zero_rdp_0d 
  module procedure is_zero_rdp_1d
 end interface 

 interface bisect
  module procedure bisect_rdp
  module procedure bisect_int
 end interface

 interface imax_loc 
  module procedure imax_loc_int
  module procedure imax_loc_rdp
 end interface 

 interface imin_loc 
  module procedure imin_loc_int
  module procedure imin_loc_rdp
 end interface 

 interface linfit
  module procedure linfit_rdp  
  module procedure linfit_gwpc 
 end interface

 interface hermitianize
  module procedure hermitianize_gwpc 
 end interface

 interface print_arr  !TODO add prtm
  module procedure print_arr1d_gwpc
  module procedure print_arr2d_gwpc
 end interface

 interface operator (.x.)           
  module procedure cross_product_
 end interface

 interface l2norm
  module procedure l2norm_rdp
 end interface


CONTAINS  !===========================================================

!!***

!!****f* m_numeric_tools/arth
!! NAME
!!  arth
!!
!! FUNCTION
!!  Returns an array of length nn containing an arithmetic progression whose
!!  starting value is start and whose step is step. 
!!
!! INPUTS
!!  start=initial point
!!  step=the increment
!!  nn=the number of points
!!
!! OUTPUT
!!  arth(nn)=the progression
!!
!! SOURCE

function arth_int(start,step,nn)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nn
 integer,intent(in) :: start,step
 integer :: arth_int(nn)

!Local variables-------------------------------
 integer :: ii
 character(len=500) :: msg
! *********************************************************************

 select case (nn)
 case (1:)
  arth_int(1)=start
  do ii=2,nn
   arth_int(ii)=arth_int(ii-1)+step
  end do
 case (0) 
  RETURN
 case (:-1)
  write(msg,'(a,i4)')'Wrong value for nn ',nn
  ABI_ERROR(msg)
 end select

end function arth_int

function arth_rdp(start,step,nn)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nn
 real(dp),intent(in) :: start,step
 real(dp) :: arth_rdp(nn)

!Local variables-------------------------------
 integer :: ii
 character(len=500) :: msg
! *********************************************************************

 select case (nn)
 case (1:)
  arth_rdp(1)=start
  do ii=2,nn
   arth_rdp(ii)=arth_rdp(ii-1)+step
  end do
 case (0) 
  RETURN
 case (:-1)
  write(msg,'(a,i4)')'Wrong value for nn ',nn
  ABI_ERROR(msg)
 end select

end function arth_rdp
!!***

!!****f* m_numeric_tools/geop
!! NAME
!!  geop
!!
!! FUNCTION
!!  Returns an array of length nn containing a geometric progression whose
!!  starting value is start and whose factor is factor!
!!
!! INPUTS
!!  start=initial point
!!  factor=the factor of the geometric progression
!!  nn=the number of points
!!
!! OUTPUT
!!  geop(nn)=the progression
!!
!! SOURCE


function geop(start,factor,nn) result(res)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: start,factor
 integer,intent(in) :: nn
 real(dp) :: res(nn)

!Local variables-------------------------------
 integer :: ii
 real(dp) :: temp
! *********************************************************************

 if (nn>0) res(1)=start
 do ii=2,nn
  res(ii)=res(ii-1)*factor
 end do

end function geop
!!***

!!****f* m_numeric_tools/set2unit
!! NAME
!!  set2unit 
!!
!! FUNCTION
!!  Set the matrix matrix to be a unit matrix (if it is square).
!!
!! SIDE EFFECTS
!!  matrix(:,:)=set to unit on exit
!!
!! SOURCE

subroutine unit_matrix_int(matrix)

!Arguments ------------------------------------

 integer,intent(inout) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,nn
! *********************************************************************

 nn=MIN(SIZE(matrix,DIM=1),SIZE(matrix,DIM=2))
 matrix(:,:)=0
 do ii=1,nn
  matrix(ii,ii)=1
 end do

end subroutine unit_matrix_int

subroutine unit_matrix_rdp(matrix)

!Arguments ------------------------------------

 real(dp),intent(inout) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: ii,nn
! *********************************************************************

 nn=MIN(SIZE(matrix,DIM=1),SIZE(matrix,DIM=2))
 matrix(:,:)=zero
 do ii=1,nn
  matrix(ii,ii)=one
 end do

end subroutine unit_matrix_rdp
!!***

!!****f* m_numeric_tools/trace
!! NAME
!!  trace 
!!
!! FUNCTION
!!  Calculate the trace of a square matrix 
!!
!! INPUT 
!!  matrix(:,:) 
!! 
!! OUTPUT 
!!  trace=the trace 
!!
!! SOURCE

function trace_int(matrix) result(trace)

!Arguments ------------------------------------

 integer :: trace
 integer,intent(in) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: nn,ii
! *********************************************************************

 nn=assert_eq(SIZE(matrix,1),SIZE(matrix,2),'Matrix not square',&
& __FILE__,__LINE__)

 trace=0
 do ii=1,nn
  trace=trace+matrix(ii,ii)
 end do

end function trace_int

function trace_rdp(matrix) result(trace)

!Arguments ------------------------------------

 real(dp) :: trace
 real(dp),intent(in) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: nn,ii
! *********************************************************************

 nn=assert_eq(SIZE(matrix,1),SIZE(matrix,2),'Matrix not square',&
& __FILE__,__LINE__)

 trace=zero
 do ii=1,nn
  trace=trace+matrix(ii,ii)
 end do

end function trace_rdp

function trace_cdp(matrix) result(trace)

!Arguments ------------------------------------

 complex(dpc) :: trace
 complex(dpc),intent(in) :: matrix(:,:)

!Local variables-------------------------------
!scalars
 integer :: nn,ii
! *********************************************************************

 nn=assert_eq(SIZE(matrix,1),SIZE(matrix,2),'Matrix not square',&
& __FILE__,__LINE__)

 trace=czero
 do ii=1,nn
  trace=trace+matrix(ii,ii)
 end do

end function trace_cdp
!!***


!!****f* m_numeric_tools/get_diag
!! NAME
!!  get_diag
!!
!! FUNCTION
!!  Return the trace of a square matrix as a vector
!!
!! INPUT 
!!  matrix(:,:) 
!! 
!! OUTPUT 
!!  diag(:)=the diagonalr
!!
!! SOURCE

function get_diag_rdp(mat) result(diag)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: mat(:,:)
 real(dp) :: diag(SIZE(mat,1))

!Local variables-------------------------------
 integer :: ii
! *************************************************************************

 ii=assert_eq(SIZE(mat,1),SIZE(mat,2),'Matrix not square',&
& __FILE__,__LINE__)

 do ii=1,SIZE(mat,1)
  diag(ii)=mat(ii,ii)
 end do

end function get_diag_rdp
!!***

!!****f* m_numeric_tools/r2c
!! NAME
!!  r2c
!!
!! FUNCTION
!!  Create a complex array starting from a real array containing real and imaginary part
!! 
!! INPUTS
!!  rr(:)=the real array
!!
!! OUTPUT
!!  cc(:)=the complex array
!!
!! SOURCE

function rdp2cdp_1d(rr) result(cc)

!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 real(dp),intent(in) :: rr(:)
 complex(dpc) :: cc(SIZE(rr)/2)

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 if (.not.is_even(SIZE(rr))) then 
  ABI_ERROR('Size of rr is not even')
 end if
 do ii=1,SIZE(cc)
  cc(ii)=CMPLX(rr(2*ii-1),rr(2*ii))
 end do

end function rdp2cdp_1d
!!***

!!****f* m_numeric_tools/c2r
!! NAME
!!  c2r
!!
!! FUNCTION
!!  Create a real array containing real and imaginary part starting from a complex array 
!! 
!! INPUTS
!!  cc(:)=the input complex array
!!
!! OUTPUT
!!  rr(:)=the real array
!!
!! SOURCE

function cdp2rdp_1d(cc) result(rr)

!Arguments ------------------------------------
!scalars

 complex(dpc),intent(in) :: cc(:)
 real(dp) :: rr(SIZE(cc)*2)

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 do ii=1,SIZE(cc)
  rr(2*ii-1)=REAL (cc(ii))
  rr(2*ii  )=AIMAG(cc(ii))
 end do

end function cdp2rdp_1d
!!***

!!****f* m_numeric_tools/is_even
!! NAME
!!  is_even
!!
!! FUNCTION
!!  Return .TRUE. if the given integer is even 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function is_even(nn)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nn
 logical :: is_even
! *********************************************************************

 is_even=.FALSE. ; if ((nn/2)*2==nn) is_even=.TRUE.

end function is_even
!!***

!!****f* m_numeric_tools/is_integer
!! NAME
!!  is_integer
!!
!! FUNCTION
!!  Return .TRUE. if all elements differ from an integer by less that tol 
!!
!! INPUTS
!!  rr=the set of real values to be checked
!!  tol=tolerance on the difference between real and integer
!!
!! SOURCE

function is_integer_0d(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: tol
 logical :: ans
!arrays
 real(dp),intent(in) :: rr

!Local variables-------------------------------
!scalars
 integer :: ii
! *************************************************************************

 ans=(ABS(rr-NINT(rr))<tol) 

end function is_integer_0d

function is_integer_1d(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: tol
 logical :: ans
!arrays
 real(dp),intent(in) :: rr(:)

!Local variables-------------------------------
!scalars
 integer :: ii
! *************************************************************************

 ans=ALL((ABS(rr-NINT(rr))<tol))

end function is_integer_1d
!!***

!!****f* m_numeric_tools/is_zero
!! NAME
!!  is_zero
!!
!! FUNCTION
!!  Return .TRUE. if all elements differ from zero by less that tol 
!!
!! INPUTS
!!  rr=the set of real values to be checked
!!  tol=tolerance
!!
!! SOURCE

function is_zero_rdp_0d(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: tol
 logical :: ans
!arrays
 real(dp),intent(in) :: rr
! *************************************************************************

 ans=(ABS(rr)<tol) 

end function is_zero_rdp_0d

function is_zero_rdp_1d(rr,tol) result(ans)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: tol
 logical :: ans
!arrays
 real(dp),intent(in) :: rr(:)
! *************************************************************************

 ans=ALL(ABS(rr(:))<tol)

end function is_zero_rdp_1d
!!***

!!****f* m_numeric_tools/bisect
!! NAME
!!  bisect
!!
!! FUNCTION
!!  Given an array AA(1:N), and a value x, returns the index j such that AA(j)<=x<= AA(j + 1). 
!!  AA must be monotonic, either increasing or decreasing. j=0 or
!!  j=N is returned to indicate that x is out of range.
!!
!! SOURCE

function bisect_rdp(AA,xx) result(loc)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: AA(:)
 real(dp),intent(in) :: xx
 integer :: loc

!Local variables-------------------------------
 integer :: nn,jl,jm,ju
 logical :: ascnd
! *********************************************************************

 nn=SIZE(AA) ; ascnd=(AA(nn)>=AA(1)) 
 !
 ! === Initialize lower and upper limits ===
 jl=0 ; ju=nn+1
 do
  if (ju-jl<=1) EXIT
  jm=(ju+jl)/2  ! Compute a midpoint,
  if (ascnd.EQV.(xx>=AA(jm))) then
   jl=jm ! Replace lower limit
  else
   ju=jm ! Replace upper limit
  end if
 end do
 !
 ! === Set the output, being careful with the endpoints ===
 if (xx==AA(1)) then
  loc=1
 else if (xx==AA(nn)) then
  loc=nn-1
 else
  loc=jl
 end if

end function bisect_rdp

function bisect_int(AA,xx) result(loc)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: AA(:)
 integer,intent(in) :: xx
 integer :: loc

!Local variables-------------------------------
 integer :: nn,jl,jm,ju
 logical :: ascnd
! *********************************************************************

 nn=SIZE(AA) ; ascnd=(AA(nn)>=AA(1)) 
 !
 ! === Initialize lower and upper limits ===
 jl=0 ; ju=nn+1
 do
  if (ju-jl<=1) EXIT
  jm=(ju+jl)/2  ! Compute a midpoint,
  if (ascnd.EQV.(xx>=AA(jm))) then
   jl=jm ! Replace lower limit
  else
   ju=jm ! Replace upper limit
  end if
 end do
 !
 ! === Set the output, being careful with the endpoints ===
 if (xx==AA(1)) then
  loc=1
 else if (xx==AA(nn)) then
  loc=nn-1
 else
  loc=jl
 end if

end function bisect_int
!!***

!!****f* m_numeric_tools/imaxloc
!! NAME
!!  imaxloc
!!
!! FUNCTION
!!  Index of maxloc on an array returned as scalar instead of array-valued
!!
!! SOURCE

function imax_loc_int(iarr)

!Arguments ------------------------------------

 integer,intent(in) :: iarr(:)
 integer :: imax_loc_int
!Local variables-------------------------------
 integer :: imax(1)
! *************************************************************************

 imax=MAXLOC(iarr(:)) ; imax_loc_int=imax(1)

end function imax_loc_int

function imax_loc_rdp(arr)

!Arguments ------------------------------------

 real(dp),intent(in) :: arr(:)
 integer :: imax_loc_rdp
!Local variables-------------------------------
 integer :: imax(1)
! *************************************************************************

 imax=MAXLOC(arr(:)) ; imax_loc_rdp=imax(1)

end function imax_loc_rdp
!!***

!!****f* m_numeric_tools/iminloc
!! NAME
!!  iminloc
!!
!! FUNCTION
!!  Index of minloc on an array returned as scalar instead of array-valued
!!
!! SOURCE

function imin_loc_int(arr)

!Arguments ------------------------------------

 integer,intent(in) :: arr(:)
 integer :: imin_loc_int
!Local variables-------------------------------
 integer :: imin(1)
! *************************************************************************

 imin=MINLOC(arr(:)) ; imin_loc_int=imin(1)

end function imin_loc_int

function imin_loc_rdp(arr)

!Arguments ------------------------------------

 real(dp),intent(in) :: arr(:)
 integer :: imin_loc_rdp
!Local variables-------------------------------
 integer :: imin(1)
! *************************************************************************

 imin=MINLOC(arr(:)) ; imin_loc_rdp=imin(1)

end function imin_loc_rdp
!!***

!!****f* m_numeric_tools/linear_fit
!! NAME
!!  linear_fit
!!
!! FUNCTION
!!  Perform a linear fit, y=ax+b, of data
!!
!! INPUTS
!!  xx(nn)=xx coordinates
!!  yy(nn)=yy coordinates
!!
!! OUTPUT
!!  aa=coefficient of linear term of fit
!!  bb=coefficient of constant term of fit
!!  function linfit=root mean square of differences between data and fit
!!
!! SOURCE

function linfit_rdp(nn,xx,yy,aa,bb) result(res) 

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nn
 real(dp) :: res
 real(dp),intent(out) :: aa,bb
!arrays
 real(dp),intent(in) :: xx(nn),yy(nn) 
 
!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: msrt,sx2,sx,sxy,sy,tx,ty 
! *************************************************************************

 sx=zero ; sy=zero ; sxy=zero ; sx2=zero
 do ii=1,nn
  tx=xx(ii)
  ty=yy(ii)
  sx=sx+tx
  sy=sy+ty
  sxy=sxy+tx*ty
  sx2=sx2+tx*tx
 end do

 aa=(nn*sxy-sx*sy)/(nn*sx2-sx*sx)
 bb=sy/nn-sx*aa/nn

 msrt=zero
 do ii=1,nn
  tx=xx(ii)
  ty=yy(ii)
  msrt=msrt+(ty-aa*tx-bb)**2
 end do
 msrt=SQRT(msrt/nn) ; res=msrt

end function linfit_rdp

function linfit_gwpc(nn,xx,zz,aa,bb) result(res)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nn
 real(dp) :: res
 real(dp),intent(in) :: xx(nn)
 complex(gwpc),intent(in) :: zz(nn)
 complex(gwpc),intent(out) :: aa,bb
!arrays

!Local variables-------------------------------
!scalars
 integer :: ii
 real(dp) :: sx,sx2,msrt
 complex(gwpc) :: sz,sxz
! *************************************************************************

 sx=zero ; sx2=zero ; msrt=zero !; sz=czero ; sxz=czero 
 sz=(0.0_gwp,0.0_gwp) ; sxz=(0.0_gwp,0.0_gwp)
 do ii=1,nn
  sx=sx+xx(ii)
  sz=sz+zz(ii)
  sxz=sxz+xx(ii)*zz(ii)
  sx2=sx2+xx(ii)*xx(ii)
 end do

 aa=(nn*sxz-sx*sz)/(nn*sx2-sx*sx)
 bb=sz/nn-sx*aa/nn

 do ii=1,nn
  msrt=msrt+ABS(zz(ii)-aa*xx(ii)-bb)**2
 end do
 msrt=SQRT(msrt) ; res=msrt

end function linfit_gwpc
!!***

!!****f* m_numeric_tools/llsfit_svd
!! NAME
!!  llsfit_svd
!!
!! FUNCTION
!!  Given a set of N data points (x,y) with individual standard deviations sigma_i, 
!!  use chi-square minimization to determine the M coefficients, par, of a function that 
!!  depends linearly on nfuncs functions, i.e f(x) = \sum_i^{nfuncs} par_i * func_i(x). 
!!  Solve the ﬁtting equations using singular value decomposition of the design matrix as in Eq 14.3.17
!!  of Numerical Recipies. The program returns values for the M ﬁt parameters par, and chi-square. 
!!  The user supplies a subroutine funcs(x,nfuncs) that returns the M basis functions evaluated at xx.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine llsfit_svd(xx,yy,sigma,nfuncs,funcs,chisq,par,var,cov,info)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nfuncs
 integer,intent(out) :: info
 real(dp),intent(out) :: chisq
!arrays
 real(dp),intent(in) :: xx(:),yy(:),sigma(:)
 real(dp),intent(out) :: par(:),var(:),cov(:,:)

 interface
  function funcs(xx,nf)
  use defs_basis
  implicit none
  real(dp),intent(in) :: xx
  integer,intent(in) :: nf
  real(dp) :: funcs(nf)
  end function funcs
 end interface

!Local variables-------------------------------
 integer,parameter :: PAD_=50
 integer :: ii,npts,lwork
 real(dp),parameter :: TOL_=1.0e-5_dp
 !character(len=500) :: msg      
 logical :: test
!arrays
 real(dp),dimension(SIZE(xx)) :: bb,sigm1
 real(dp),dimension(SIZE(xx),nfuncs) :: dmat,dmat_save
 real(dp) :: tmp(nfuncs)
 real(dp),allocatable :: work(:),Vt(:,:),U(:,:),S(:)
! *************************************************************************

 npts=assert_eq(SIZE(xx),SIZE(yy),SIZE(sigma),'Wrong size in xx,yy,sigma',&
& __FILE__,__LINE__)
 call assert((npts>=nfuncs),'No. of functions must greater than no. of points',&
& __FILE__,__LINE__)
 ii=assert_eq(nfuncs,SIZE(cov,1),SIZE(cov,2),SIZE(var),'Wrong size in covariance',&
& __FILE__,__LINE__)
 !
 ! === Calculate design matrix and b vector ===
 ! * dmat_ij=f_j(x_i)/sigma_i, b_i=y_i/sigma_i
 sigm1(:)=one/sigma(:) ; bb(:)=yy(:)*sigm1(:)
 do ii=1,npts
  dmat_save(ii,:)=funcs(xx(ii),nfuncs) 
 end do
 dmat=dmat_save*SPREAD(sigm1,DIM=2,ncopies=nfuncs)
 dmat_save(:,:)=dmat(:,:)
 !
 ! === Singular value decomposition ===
 lwork=MAX(3*MIN(npts,nfuncs)+MAX(npts,nfuncs),5*MIN(npts,nfuncs)-4)+PAD_
 allocate(work(lwork),U(npts,npts),S(nfuncs),Vt(nfuncs,nfuncs)) 

 !FIXME
 STOP "Fix problem with build system"
!££ MG Commented out due to problems with the build system (to be linked against LAPACK)
!£ call DGESVD('A','A',npts,nfuncs,dmat,npts,S,U,npts,Vt,nfuncs,work,lwork,info)
 deallocate(work) ; if (info/=0) GOTO 10
 !
 ! === Set to zero small singular values according to TOL_ and find coefficients ===
 WHERE (S>TOL_*MAXVAL(S)) 
  tmp=MATMUL(bb,U)/S
 ELSE WHERE
  S=zero
  tmp=zero
 END WHERE
 par(:)=MATMUL(tmp,Vt)
 !
 ! === Evaluate chi-square ===
 chisq=l2norm(MATMUL(dmat_save,par)-bb)**2
 !
 ! === Calculate covariance and variance ===
 ! C_jk = V_ji V_ki / S_i^2
 WHERE (S/=zero) S=one/(S*S)

 ! check this but should be correct
 cov(:,:)=Vt*SPREAD(S,DIM=2,ncopies=nfuncs)
 cov(:,:)=MATMUL(TRANSPOSE(Vt),cov)
 var(:)=SQRT(get_diag(cov))

10 deallocate(U,S,Vt) 

end subroutine llsfit_svd
!!***


!!****f* m_numeric_tools/polyn_interp
!! NAME
!!  polyn_interp
!!
!! FUNCTION
!!  Given arrays xa and ya of length N, and given a value x, return a value y, and an error estimate dy. 
!!  If P(x) is the polynomial of degree N−1 such that P(xai)=yai, i=1,...,N, then the returned value y=P(x).
!!
!! INPUTS
!!  xa(:)=abscissas in ascending order
!!  ya(:)=ordinates
!!  x=the point where the set of data has to be interpolated
!!
!! OUTPUT
!!  y=the interpolated value
!!  dy=error estimate
!!
!! NOTES
!!  Based on the polint routine reported in Numerical Recipies
!!
!! PARENTS
!!      m_numeric_tools
!!
!! CHILDREN
!!
!! SOURCE

subroutine polyn_interp(xa,ya,x,y,dy)
!recursive subroutine polyn_interp(xa,ya,x,y,dy)

!Arguments ------------------------------------
!scalars

 real(dp),intent(in) :: xa(:),ya(:)
 real(dp),intent(in) :: x
 real(dp),intent(out) :: y,dy
!Local variables-------------------------------
!scalars
 integer :: m,n,ns
!arrays
 real(dp),dimension(SIZE(xa)) :: c,d,den,ho
! *************************************************************************

 n=assert_eq(SIZE(xa),SIZE(ya),'Different size in xa and ya',&
& __FILE__,__LINE__)

 ! === Initialize the tables of c and d ===
 c(:)=ya(:) ; d(:)=ya(:) ; ho(:)=xa(:)-x
 ! === Find closest table entry and initial approximation to y ===
 ns=imin_loc(ABS(x-xa)) ; y=ya(ns)
 ns=ns-1
 !
 ! === For each column of the tableau loop over current c and d and up-date them ===
 do m=1,n-1
  den(1:n-m)=ho(1:n-m)-ho(1+m:n)
  if (ANY(den(1:n-m)==zero)) then 
   ABI_ERROR('Two input xa are identical')
  end if

  den(1:n-m)=(c(2:n-m+1)-d(1:n-m))/den(1:n-m)
  d(1:n-m)=ho(1+m:n)*den(1:n-m) ! Update c and d 
  c(1:n-m)=ho(1:n-m)*den(1:n-m)

  if (2*ns<n-m) then  ! Now decide which correction, c or d, we want to add to the 
   dy=c(ns+1)         ! accumulating value of y, The last dy added is the error indication.
  else
   dy=d(ns)
   ns=ns-1
  end if

  y=y+dy
 end do

end subroutine polyn_interp
!!***


!!****f* m_numeric_tools/trapezoidal_
!! NAME
!!  trapezoidal_ (PRIVATE)
!!
!! FUNCTION
!!  Compute the n-th stage of refinement of an extended trapezoidal rule
!!  adding 2^(n-2) additional interior point in the finite range of integration
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  nn=integer defining the refinement of the mesh, each call adds 2^(n-2) additional interior points 
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  quad=the integral at the n-th stage. 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!  When called with nn=1, the routine returns the crudest estimate of the integral
!!  Subsequent calls with nn=2,3,... (in that sequential order) will improve the accuracy
!!  by adding 2^(n-2) additional interior points. Note that quad should not be modified between sequential calls.
!!  Subroutine is defined as recursive to allow multi-dimensional integrations
!!
!! SOURCE

recursive subroutine trapezoidal_(func,nn,xmin,xmax,quad)
          !subroutine trapezoidal_(func,nn,xmin,xmax,quad)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nn
 !real(dp),external :: func
 real(dp),intent(in) :: xmin,xmax
 real(dp),intent(inout) :: quad

 interface
  function func(x)
   use defs_basis
   real(dp),intent(in) :: x
   real(dp) :: func
  end function func
 end interface

 !interface
 ! function func(x)
 !  use defs_basis
 !  real(dp),intent(in) :: x(:)
 !  real(dp) :: func(SIZE(x))
 ! end function func
 !end interface

!Local variables-------------------------------
!scalars
 integer :: npt,ix,ii
 real(dp) :: space,new,yy
 character(len=500) :: msg
!arrays
 !real(dp),allocatable :: xx(:)
!************************************************************************

 select case (nn)
 case (1)
  ! === Initial crude estimate (xmax-xmin)(f1+f2)/2 ===
  !quad=half*(xmax-xmin)*SUM(func((/xmin,xmax/)))
  quad=half*(xmax-xmin)*(func(xmin)+func(xmax))
 case (2:)
  ! === Add npt interior points of spacing space ===
  npt=2**(nn-2) ; space=(xmax-xmin)/npt 
  ! === The new sum is combined with the old integral to give a refined integral ===
  !new=SUM(func(arth(xmin+half*space,space,npt))) !PARALLEL version
  !allocate(xx(npt)) 
  !xx(:)=arth(xmin+half*space,space,npt)
  !xx(1)=xmin+half*space
  !do ii=2,nn
  ! xx(ii)=xx(ii-1)+space
  !end do
  new=zero
  yy=xmin+half*space
  do ix=1,npt
   !new=new+func(xx(ix))
   new=new+func(yy)
   yy=yy+space
  end do
  !deallocate(xx)
  quad=half*(quad+space*new) 
  !print*,'trapezoidal',quad
 case (:0)
  write(msg,'(a,i3)')'Wrong value for nn ',nn
  ABI_ERROR(msg)
 end select

end subroutine trapezoidal_
!!***

!!****f* m_numeric_tools/midpoint_
!! NAME
!!  midpoint_ (PRIVATE)
!!
!! FUNCTION
!!  This routine computes the n-th stage of refinement of an extended midpoint rule.
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  nn=integer defining the refinement of the mesh, each call adds (2/3)*3n-1 additional 
!!   interior points between xmin ans xmax
!!
!! OUTPUT
!!  See SIDE EFFECTS
!!
!! SIDE EFFECTS
!!  quad=the integral at the n-th stage. 
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!  When called with nn=1, the routine returns as quad the crudest estimate of the integral
!!  Subsequent calls with nn=2,3,... (in that sequential order) will improve the accuracy of quad by adding
!!  (2/3)×3n-1 additional interior points. quad should not be modified between sequential calls.
!!  Subroutine is defined as recursive to allow multi-dimensional integrations
!!
!! SOURCE

 recursive subroutine midpoint_(func,nn,xmin,xmax,quad)
          !subroutine midpoint_(func,nn,xmin,xmax,quad)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: nn
 !real(dp),external :: func
 real(dp),intent(in) :: xmin,xmax
 real(dp),intent(inout) :: quad

 interface
  function func(x)
   use defs_basis
   real(dp),intent(in) :: x
   real(dp) :: func
  end function func
 end interface

 !interface
 ! function func(x)
 !  use defs_basis
 !  real(dp),intent(in) :: x(:)
 !  real(dp) :: func(SIZE(x))
 ! end function func
 !end interface

!Local variables-------------------------------
!scalars
 integer  :: npt,ix
 real(dp) :: space
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: xx(:)

!************************************************************************

 select case (nn)
 case (1)
  ! === Initial crude estimate done at the middle of the interval
  !quad=(xmax-xmin)*SUM(func((/half*(xmin+xmax)/))) !PARALLEL version
  quad=(xmax-xmin)*func(half*(xmin+xmax))
 case (2:)
  ! === Add npt interior points, they alternate in spacing between space and 2*space ===
  allocate(xx(2*3**(nn-2)))
  npt=3**(nn-2) ; space=(xmax-xmin)/(three*npt) 
  xx(1:2*npt-1:2)=arth(xmin+half*space,three*space,npt) 
  xx(2:2*npt:2)=xx(1:2*npt-1:2)+two*space
  ! === The new sum is combined with the old integral to give a refined integral ===
  !quad=quad/three+space*SUM(func(xx))  !PARALLEL version
  quad=quad/three
  do ix=1,SIZE(xx)
   quad=quad+space*func(xx(ix)) 
  end do
  deallocate(xx)
 case (:0)
  write(msg,'(a,i3)')' wrong value for nn ',nn
  ABI_BUG('Wrong value for nn')
 end select

end subroutine midpoint_
!!***

!!****f* m_numeric_tools/quadrature_
!! NAME
!!  quadrature
!!
!! FUNCTION
!!  Driver routine to perform quadratures in finite domains using different techniques.
!!  The routine improves the resolution of the grid until a given accuracy is reached
!!
!! INPUTS
!!  func(external)=the name of the function to be integrated
!!  xmin,xmax=the limits of integration
!!  npts=Initial number of points, only for Gauss-Legendre. At each step this number is doubled
!!  accuracy=fractional accuracy required
!!  qopt=integer flag defining the algorithm for the quadrature:
!!  ntrial=Max number of attempts
!!    1 for Trapezoidal rule, closed, O(1/N^2) 
!!    2 for Simpson based on trapezoidal,closed, O(1/N^4)
!!    3 for Midpoint rule, open, O(1/N^2)
!!    4 for midpoint rule with cancellation of leading error, open, O(1/N^4)
!!    5 for Romber integration (closed form) and extrapolation for h-->0 (order 10 is hard-coded)
!!    6 for Romber integration with midpoint rule and extrapolation for h-->0 (order 10 is hard-coded)
!!    7 for Gauss-Legendre
!!
!! OUTPUT
!!  quad=the integral
!!
!! PARENTS
!!
!! CHILDREN
!!
!! NOTES
!!
!! SOURCE

!subroutine quadrature(func,xmin,xmax,qopt,quad,ntrial,accuracy,npts)
recursive subroutine quadrature(func,xmin,xmax,qopt,quad,ntrial,accuracy,npts)

!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_lib00numeric
!End of the abilint section

 integer,intent(in) :: qopt
 integer,optional,intent(in) :: ntrial,npts
 !real(dp),external :: func
 real(dp),intent(in) :: xmin,xmax
 real(dp),optional,intent(in) :: accuracy
 real(dp),intent(out) :: quad

 interface
  function func(x)
   use defs_basis
   real(dp),intent(in) :: x
   real(dp) :: func
  end function func
 end interface

 !interface
 ! function func(x)
 !  use defs_basis
 !  real(dp),intent(in) :: x(:)
 !  real(dp) :: func(SIZE(x))
 ! end function func
 !end interface

!Local variables-------------------------------
!scalars
 integer :: K,KM,NT,NX,NX0,it,ix
 !integer,save :: ncall=0
 real(dp) :: EPS,old_st,st,old_quad,dqromb
 real(dp) :: TOL
 character(len=500) :: msg      
!arrays
 real(dp),allocatable :: h(:),s(:)
 real(dp),allocatable :: wx(:),xx(:)
! *************************************************************************

 !ncall=ncall+1
 TOL=tol12
 EPS=tol6  ; if (PRESENT(accuracy)) EPS=accuracy
 NT=20     ; if (PRESENT(ntrial))   NT=ntrial
 quad=zero
                                          
 select case (qopt)
 case (1)
  ! === Trapezoidal, closed form, O(1/N^2) 
  do it=1,NT
   call trapezoidal_(func,it,xmin,xmax,quad)
   if (it>5) then ! Avoid spurious early convergence
    if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
   end if
   old_quad=quad
  end do
 case (2)
  ! === Extended Simpson rule based on trapezoidal O(1/N^4) ===
  do it=1,NT
   call trapezoidal_(func,it,xmin,xmax,st)
   if (it==1) then 
    quad=st
   else
    quad=(four*st-old_st)/three
   end if
   if (it>5) then ! Avoid spurious early convergence
    if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
   end if
   old_quad=quad
   old_st=st
  end do
 case (3) 
  ! === Midpoint rule, open form, O(1/N^2) ===
  do it=1,NT
   call midpoint_(func,it,xmin,xmax,quad)
   if (it>4) then ! Avoid spurious early convergence
    if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
   end if
   old_quad=quad
  end do
 case (4) 
  ! === Midpoint rule with cancellation of leading 1/N^2 term, open form, O(1/N^4) ===
  do it=1,NT
   call midpoint_(func,it,xmin,xmax,st)
   if (it==1) then 
    quad=st
   else
    quad=(nine*st-old_st)/eight
   end if
   if (it>4) then ! Avoid spurious early convergence
    if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
   end if
   old_quad=quad
   old_st=st
  end do
 case (5) 
  ! === Romberg Integration, closed form ===
  K=5 ; KM=K-1 ! Order 10
  allocate(h(NT+1),s(NT+1))  ; h=zero ; s=zero
  h(1)=one
  do it=1,NT
   call trapezoidal_(func,it,xmin,xmax,s(it))
   !print*,' romberg-trap at ',ncall,it,s(it)
   if (it>=K) then 
    call polyn_interp(h(it-KM:it),s(it-KM:it),zero,quad,dqromb)
    if (ABS(dqromb)<EPS*ABS(quad)) then 
     deallocate(h,s) ; RETURN
    end if
   end if
   s(it+1)=s(it)
   h(it+1)=quarter*h(it) ! Quarter makes the extrapolation a polynomial in h^2, 
  end do                 ! This is required to use the Euler-Maclaurin formula 
  deallocate(h,s)
 case (6) 
  ! === Romberg Integration, closed form ===
  K=5 ; KM=K-1 ! Order 10
  allocate(h(NT+1),s(NT+1))  ; h=zero ; s=zero
  h(1)=one
  do it=1,NT
   call midpoint_(func,it,xmin,xmax,s(it))
   if (it>=K) then 
    call polyn_interp(h(it-KM:it),s(it-KM:it),zero,quad,dqromb)
    !print*,quad,dqromb
    if (ABS(dqromb)<EPS*ABS(quad)) then 
     deallocate(h,s) ; RETURN
    end if
   end if
   s(it+1)=s(it)
   h(it+1)=ninth*h(it) ! factor is due to step tripling in midpoint and even error series
  end do                 
  deallocate(h,s)
 case (7) 
  ! === Gauss-Legendre ===
  write(*,*)"Not implemented, waiting for fix build system"
  STOP
  NX0=5 ; if (PRESENT(npts)) NX0=npts 
  NX=NX0
  do it=1,NT
   allocate(wx(NX),xx(NX))
   !FIXME
    STOP "Fix problem with build system!"
   !£ call coeffs_gausslegint(xmin,xmax,xx,wx,NX)
   quad=zero 
   do ix=1,NX 
    quad=quad+wx(ix)*func(xx(ix))
   end do
   deallocate(wx,xx)
   if (it>1) then 
    !print*,quad
    if (ABS(quad-old_quad)<EPS*ABS(old_quad).or.(ABS(quad)<TOL.and.ABS(old_quad)<TOL)) RETURN
   end if
   old_quad=quad
   NX=NX+NX0
   !NX=2*NX
  end do
 case default 
  write(msg,'(a,i3)')'Wrong value for qopt',qopt
  ABI_BUG(msg)
 end select

 write(msg,'(a)')'Number of trials not sufficient to get converged results '
 ABI_ERROR(msg)
 RETURN

end subroutine quadrature
!!***

!!****f* m_numeric_tools/hermitianize
!! NAME
!!  hermitianize
!!
!! FUNCTION
!!  Force a square matrix to be hermitian
!!
!! INPUTS
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  mat(:,:)=complex input matrix, hermitianized at output
!!
!! CHILDREN
!!
!! SOURCE

subroutine hermitianize_gwpc(mat)

!Arguments ------------------------------------
!scalars
!arrays

 complex(gwpc),intent(inout) :: mat(:,:)

!Local variables-------------------------------
!scalars
 integer :: nn,ii,jj
!arrays
 complex(gwpc),allocatable :: tmp(:)
! *************************************************************************
 
 nn=assert_eq(SIZE(mat,1),SIZE(mat,2),'Matrix not square',&
& __FILE__,__LINE__)
 allocate(tmp(nn))

 do ii=1,nn
  do jj=ii,nn
   tmp(jj)=(mat(ii,jj)+CONJG(mat(jj,ii)))/two
  end do
  !do jj=ii,nn ; mat(ii,jj)=tmp(jj) ; end do
  mat(ii,ii:nn)=tmp(ii:nn)
  !do jj=ii,nn ; mat(jj,ii)=CONJG(tmp(jj)) ; end do
  mat(ii:nn,ii)=CONJG(tmp(ii:nn))
 end do

 deallocate(tmp)

end subroutine hermitianize_gwpc 
!!***

!!****f* m_numeric_tools/print_arr
!! NAME
!! print_arr
!!
!! FUNCTION
!! Print an array using a nice (?) format
!!
!! INPUTS
!!  arr(:)=vector/matrix to be printed
!!  mode_paral(optional)=parallel mode, DEFAULT is "COLL"
!!   "COLL" if all procs are calling the routine with the same message to be written only once
!!   "PERS" if the procs are calling the routine with different mesgs each to be written, 
!!          or if one proc is calling the routine
!!  unit(optional)=the unit number of the file, DEFAULT=std_out
!!  max_r,max_c(optional)=Max number of rows and columns to be printed 
!!   (DEFAULT is 9, output format assumes to be less that 99, but there might be 
!!    problems with wrtout if message size exceeds 500 thus max number of elements should be ~60)
!!
!! OUTPUT
!!  (only printing)
!!
!! SOURCE

subroutine print_arr1d_gwpc(arr,max_r,unit,mode_paral)

!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 integer,optional,intent(in) :: unit,max_r
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 complex(gwpc),intent(in) :: arr(:)

!Local variables-------------------------------
!scalars
 integer :: unt,ii,nr,mr
 character(len=4) :: mode
 character(len=500) :: msg
 character(len=100) :: fmth,fmt1
! *************************************************************************

 unt=std_out ; if (PRESENT(unit      )) unt=unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral
 mr=15       ; if (PRESENT(max_r     )) mr=max_r
 
 if (mode/='COLL'.and.mode/='PERS') then 
  write(msg,'(2a)')' Wrong value of mode_paral ',mode
  ABI_ERROR(msg)
 end if 
 !
 ! === Print out matrix ===
 nr=SIZE(arr,DIM=1) ; if (mr>nr) mr=nr

 write(fmth,*)'(6x,',mr,'(i2,6x))'
 write(fmt1,*)'(3x,',mr,'f8.3)'

 write(msg,fmth)(ii,ii=1,mr)     
 call wrtout(unt,msg,mode) !header
 write(msg,fmt1)REAL (arr(1:mr)) 
 call wrtout(unt,msg,mode) !real part
 write(msg,fmt1)AIMAG(arr(1:mr)) 
 call wrtout(unt,msg,mode) !imag part

end subroutine print_arr1d_gwpc

subroutine print_arr2d_gwpc(arr,max_r,max_c,unit,mode_paral)

!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 integer,optional,intent(in) :: unit,max_r,max_c
 character(len=4),optional,intent(in) :: mode_paral
!arrays
 complex(gwpc),intent(in) :: arr(:,:)

!Local variables-------------------------------
!scalars
 integer :: unt,ii,jj,nc,nr,mc,mr
 character(len=4) :: mode
 character(len=500) :: msg
 character(len=100) :: fmth,fmt1,fmt2
! *************************************************************************

 unt=std_out ; if (PRESENT(unit      )) unt=unit
 mode='COLL' ; if (PRESENT(mode_paral)) mode=mode_paral
 mc=9        ; if (PRESENT(max_c)) mc=max_c
 mr=9        ; if (PRESENT(max_r)) mr=max_r
 
 if (mode/='COLL'.and.mode/='PERS') then 
  write(msg,'(2a)')'Wrong value of mode_paral ',mode
  ABI_ERROR(msg)
 end if 
 !
 ! === Print out matrix ===
 nr=SIZE(arr,DIM=1) ; if (mr>nr) mr=nr
 nc=SIZE(arr,DIM=2) ; if (mc>nc) mc=nc

 write(fmth,*)'(6x,',mc,'(i2,6x))'
 write(fmt1,*)'(3x,i2,',mc,'f8.3)'
 write(fmt2,*)'(5x   ,',mc,'f8.3,a)'

 write(msg,fmth)(jj,jj=1,mc) 
 call wrtout(unt,msg,mode) !header
 do ii=1,mr
  write(msg,fmt1)ii,REAL(arr(ii,1:mc))      
  call wrtout(unt,msg,mode) !real part
  write(msg,fmt2)  AIMAG(arr(ii,1:mc)),ch10 
  call wrtout(unt,msg,mode) !imag part
 end do

end subroutine print_arr2d_gwpc
!!***

!!****f* m_numeric_tools/pade
!! NAME
!!  pade
!!
!! FUNCTION
!!  Calculate the pade approximant in zz of the function f calculated at the n points z
!!
!! SOURCE

function pade(n,z,f,zz)

!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,intent(in) :: n
 complex(dpc),intent(in) :: zz
 complex(dpc) :: pade
!arrays
 complex(dpc),intent(in) :: z(n),f(n)

!Local variables-------------------------------
!scalars
 complex(dpc) :: a(n)
 complex(dpc) :: Az(0:n), Bz(0:n)
 integer :: i
! *************************************************************************

 call calculate_pade_a(a,n,z,f)

 Az(0)=czero
 Az(1)=a(1)
 Bz(0)=cone
 Bz(1)=cone

 do i=1,n-1
  Az(i+1)=Az(i)+(zz-z(i))*a(i+1)*Az(i-1)
  Bz(i+1)=Bz(i)+(zz-z(i))*a(i+1)*Bz(i-1)
 end do
 !print*,'Bz(n)',Bz(n)
 if (REAL(Bz(n))==zero.and.AIMAG(Bz(n))==zero) print *,' Bz(n) ',Bz(n)
 pade=Az(n)/Bz(n)
 !print *, 'pade_approx ', pade_approx

end function pade
!!***

!!****f* m_numeric_tools/dpade
!! NAME
!!  dpade
!!
!! FUNCTION
!!  Calculate the derivative of the pade approximant in zz of the function f calculated at the n points z
!!
!! SOURCE

function dpade(n,z,f,zz)

!Arguments ------------------------------------
!scalars


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
!End of the abilint section

 integer,intent(in) :: n
 complex(dpc),intent(in) :: zz
 complex(dpc) :: dpade
!arrays
 complex(dpc),intent(in) :: z(n),f(n)

!Local variables-------------------------------
!scalars
 integer :: i
!arrays
 complex(dpc) :: a(n)
 complex(dpc) :: Az(0:n), Bz(0:n)
 complex(dpc) :: dAz(0:n), dBz(0:n)
! *************************************************************************

 call calculate_pade_a(a,n,z,f)

 Az(0)=czero
 Az(1)=a(1)
 Bz(0)=cone
 Bz(1)=cone
 dAz(0)=czero
 dAz(1)=czero
 dBz(0)=czero
 dBz(1)=czero

 do i=1,n-1
  Az(i+1)=Az(i)+(zz-z(i))*a(i+1)*Az(i-1)
  Bz(i+1)=Bz(i)+(zz-z(i))*a(i+1)*Bz(i-1)
  dAz(i+1)=dAz(i)+a(i+1)*Az(i-1)+(zz-z(i))*a(i+1)*dAz(i-1)
  dBz(i+1)=dBz(i)+a(i+1)*Bz(i-1)+(zz-z(i))*a(i+1)*dBz(i-1)
 end do
 !print*,'Bz(n)', Bz(n)
 if (REAL(Bz(n))==zero.and.AIMAG(Bz(n))==zero) print*,'Bz(n)',Bz(n)
 !pade_approx = Az(n) / Bz(n)
 dpade=dAz(n)/Bz(n) -Az(n)*dBz(n)/(Bz(n)*Bz(n))
 !print *, 'pade_approx ', pade_approx

end function dpade
!!***

subroutine calculate_pade_a(a,n,z,f)

!Arguments ------------------------------------
!scalars

 integer,intent(in) :: n
 complex(dpc),intent(in) :: z(n),f(n)
 complex(dpc),intent(out) :: a(n)

!Local variables-------------------------------
!scalars
 integer :: i,j
!arrays
 complex(dpc) :: g(n,n)
! *************************************************************************

 g(1,1:n)=f(1:n)
  
 do i=2,n
  do j=i,n
   if (REAL(g(i-1,j))==zero.and.AIMAG(g(i-1,j))==zero) print*,'g_i(z_j)',i,j,g(i,j)
   g(i,j)=(g(i-1,i-1)-g(i-1,j)) / ((z(j)-z(i-1))*g(i-1,j))
   !print*,'g_i(z_j)',i,j,g(i,j)
  end do
 end do
 do i=1,n
  a(i)=g(i,i)
 end do
 !print*,'a ',a(:)

end subroutine calculate_pade_a
!!***

!!****f* m_numeric_tools/newrap_step 
!! NAME
!!  newrap_step
!!
!! FUNCTION
!!  Apply single step newton-raphson method to find the root of a complex function 
!!   z_k+1=z_k-f(z_k)/(df/dz(z_k)) 
!!
!! SOURCE

function newrap_step(z,f,df)

!Arguments ------------------------------------
!scalars

 complex(dpc),intent(in) :: z,f,df
 complex(dpc) :: newrap_step

!Local variables-------------------------------
!scalars
 real(dp) :: dfm2
! *************************************************************************

 dfm2=ABS(df)*ABS(df)

 newrap_step= z - (f*CONJG(df))/dfm2
 !& z-one/(ABS(df)*ABS(df)) * CMPLX( REAL(f)*REAL(df)+AIMAG(f)*AIMAG(df), -REAL(f)*AIMAG(df)+AIMAG(f)*EAL(df) )

end function newrap_step 
!!***


!!****f* m_numeric_tools/cross_product
!! NAME
!!  cross_product
!!
!! FUNCTION
!!  Return the cross product of two vectors.
!!
function cross_product_(vec1,vec2) result(res)

!Arguments ------------------------------------

 real(dp),intent(in) :: vec1(3),vec2(3)
 real(dp) :: res(3)
! *************************************************************************

 res(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
 res(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
 res(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)

end function cross_product_
!!***


!!****f* m_numeric_tools/l2norm
!! NAME
!!  l2norm
!!
!! FUNCTION
!!  Return the length (ordinary L2 norm) of a vector.
!!

function l2norm_rdp(vec) result(res)

!Arguments ------------------------------------

 real(dp),intent(in) :: vec(:)
 real(dp) :: res
! *************************************************************************

 res=SQRT(DOT_PRODUCT(vec,vec))

end function l2norm_rdp
!!***

END MODULE m_numeric_tools
!!***


