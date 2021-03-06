subroutine zgefa(a,lda,n,ipvt,info)

#ifdef VMS
!DEC$ ATTRIBUTES ALIAS:'ZAXPY' :: zaxpy
!DEC$ ATTRIBUTES ALIAS:'IZAMAX' :: izamax
!DEC$ ATTRIBUTES ALIAS:'ZSCAL' :: zscal
#endif


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments
 integer :: lda,n,ipvt(*),info
 complex*16 :: a(lda,*)
!
!     zgefa factors a complex*16 matrix by gaussian elimination.
!
!     zgefa is usually called by zgeco, but it can be called
!     directly with a saving in time if  rcond  is not needed.
!     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
!
!     on entry
!
!        a       complex*16(lda, n)
!                the matrix to be factored.
!
!        lda     integer
!                the leading dimension of the array  a .
!
!        n       integer
!                the order of the matrix  a .
!
!     on return
!
!        a       an upper triangular matrix and the multipliers
!                which were used to obtain it.
!                the factorization can be written  a = l*u  where
!                l  is a product of permutation and unit lower
!                triangular matrices and  u  is upper triangular.
!
!        ipvt    integer(n)
!                an integer vector of pivot indices.
!
!        info    integer
!                = 0  normal value.
!                = k  if  u(k,k) .eq. 0.0 .  this is not an error
!                     condition for this subroutine, but it does
!                     indicate that zgesl or zgedi will divide by zero
!                     if called.  use  rcond  in zgeco for a reliable
!                     indication of singularity.
!
!     linpack. this version dated 08/14/78 .
!     cleve moler, university of new mexico, argonne national lab.
!
!     subroutines and functions
!
!     blas zaxpy,zscal,izamax
!     fortran dabs
!
!     internal variables
!
!Local variables
 complex*16 :: t
 integer :: j,k,kp1,l,nm1
 complex*16 :: zdum
 double precision :: cabs1
 double precision :: dreal,dimag
 complex*16 :: zdumr,zdumi

 dreal(zdumr) = zdumr
 dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
 cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
!
!     gaussian elimination with partial pivoting
!
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
!
!        find l = pivot index
!
         l = izamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
!
!        zero pivot implies this column already triangularized
!
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
!
!           interchange if necessary
!
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
!
!           compute multipliers
!
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
!
!           row elimination with column indexing
!
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n

      end subroutine zgefa
