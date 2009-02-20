      SUBROUTINE ludcmp(a,n,np,indx,id,info)


      implicit none

      INTEGER n,np,indx(n),NMAX,id,info
      REAL*8 a(np,np),TINY
      PARAMETER (NMAX=500,TINY=1.0e-20)

!     Given a matrix a(1:n,1:n), with physical dimension np by np, this
!     routine replaces it by the LU decomposition of a rowwise permutation of
!     itself. a and n are input. a is output, arranged as in equation (2.3.14)
!     above; indx(1:n) is an output vector that records the row permutation
!     effected by the partial pivoting; id is output as +- 1 depending on
!     whether the number of row interchanges was even or odd,
!     respectively. This routine is used in combination with lubksb to solve
!     linear equations or invert a matrix.

      INTEGER i,imax,j,k
      REAL*8 aamax,dum,sum,vv(NMAX) 

!      write(6,*) 'ENTERING LUDCMP...'
!      write(6,*) 'in ludcmp n=',n,' np=',np
!      write(6,201) ((a(i,j),j=1,n),i=1,n)
! 201  FORMAT('A in ludcmp ',/,3F16.8,/,3F16.8,/,3F16.8)
      id=1
      info=0
      do i=1,n 
        aamax=0.
        do j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
        enddo
        if (aamax.eq.0.) then
          write(6,*) 'LUDCMP: singular matrix !!!'
          do j=1,3
            write(6,*) (a(j,k),k=1,3)
          enddo
          info=1
          return
!          stop 'singular matrix in ludcmp' 
        endif
        vv(i)=1./aamax 
      enddo 
      do j=1,n 
        do i=1,j-1 
          sum=a(i,j)
          do k=1,i-1
            sum=sum-a(i,k)*a(k,j)
          enddo 
          a(i,j)=sum
        enddo 
        aamax=0. 
        do i=j,n 
          sum=a(i,j)
          do k=1,j-1
            sum=sum-a(i,k)*a(k,j)
          enddo 
          a(i,j)=sum
          dum=vv(i)*abs(sum) 
          if (dum.ge.aamax) then 
            imax=i
            aamax=dum
          endif
        enddo 
        if (j.ne.imax)then 
          do  k=1,n 
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
          enddo 
          id=-id 
          vv(imax)=vv(j) 
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then 
          dum=1./a(j,j)
          do i=j+1,n
            a(i,j)=a(i,j)*dum
          enddo
        endif
      enddo 
!      write(6,*) 'LEAVING LUDCMP...'
      return
      END


      SUBROUTINE lubksb(a,n,np,indx,b)


      IMPLICIT NONE

      INTEGER n,np,indx(n)
      REAL*8 a(np,np),b(n)

!     Solves the set of n linear equations A . X = B. Here a is input, not as
!     the matrix A but rather as its LU decomposition, determined by the
!     routine ludcmp. indx is input as the permutation vector returned by
!     ludcmp. b(1:n) is input as the right-hand side vector B, and returns
!     with the solution vector X. a, n, np, and indx are not modified by this
!     routine and can be left in place for successive calls with different
!     right-hand sides b. This routine takes into account the possibility that
!     b will begin with many zero elements, so it is ecient for use in
!     matrix inversion.

      INTEGER i,ii,j,ll
      REAL*8 sum
!      write(6,*) 'ENTERING LUBKSB...'
!      write(6,201) ((a(i,j),j=1,n),i=1,n)
! 201  FORMAT('A in lubksb ',/,3F16.8,/,3F16.8,/,3F16.8)

      ii=0
      do i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do j=ii,i-1
            sum=sum-a(i,j)*b(j)
          enddo 
        else if (sum.ne.0.) then
          ii=i 
        endif
        b(i)=sum
      enddo 
      do i=n,1,-1
        sum=b(i)
        do j=i+1,n
          sum=sum-a(i,j)*b(j)
        enddo
        b(i)=sum/a(i,i) 
      enddo
!      write(6,*) 'LEAVING LUBKSB...'
      return 
      END SUBROUTINE LUBKSB


      subroutine gaussj(a,n,np,b,m,mp,info)


      implicit none

      integer m,mp,n,np,info
      real*8 a(np,np),b(np,np)
      integer, parameter :: NMAX=50
      
      integer i,icol,irow,j,k,l,ll,indxc(NMAX),indxr(NMAX),ipiv(NMAX)
      real*8 big,dum,pivinv
      
      info=0

      do j=1,n
        ipiv(j)=0
      enddo
      do i=1,n
        big=0.
        do j=1,n
          if(ipiv(j).ne.1) then
            do k=1,n
              if(ipiv(k).eq.0) then
                if(abs(a(j,k)).ge.big) then
                  big=abs(a(j,k))
                  irow=j
                  icol=k
                endif
              else if(ipiv(k).gt.1) then
                write(6,*) 'GAUSSJ: singular matrix !!!'
                info=1
                return
              endif
            enddo
          endif
        enddo
        ipiv(icol)=ipiv(icol)+1
        if(irow.ne.icol) then
          do l=1,n
            dum=a(irow,l)
            a(irow,l)=a(icol,l)
            a(icol,l)=dum
          enddo
          do l=1,m
            dum=b(irow,l)
            b(irow,l)=b(icol,l)
            b(icol,l)=dum
          enddo
        endif
        indxr(i)=irow
        indxc(i)=icol
        if(a(icol,icol).eq.0.) then
          write(6,*) 'GAUSSJ: singular matrix !!!'
          info=1
          return
        endif
        pivinv=1./a(icol,icol)
        a(icol,icol)=1.
        do l=1,n
          a(icol,l)=a(icol,l)*pivinv
        enddo
        do l=1,m
          b(icol,l)=b(icol,l)*pivinv
        enddo
        do ll=1,n
          if(ll.ne.icol) then
            dum=a(ll,icol)
            a(ll,icol)=0.
            do l=1,n
              a(ll,l)=a(ll,l)-a(icol,l)*dum
            enddo
            do l=1,m
              b(ll,l)=b(ll,l)-b(icol,l)*dum
            enddo
          endif
        enddo
      enddo
      do l=n,1,-1
        if(indxr(l).ne.indxc(l)) then
          do k=1,n
            dum=a(k,indxr(l))
            a(k,indxr(l))=a(k,indxc(l))
            a(k,indxc(l))=dum
          enddo
        endif
      enddo
      return
      end subroutine gaussj

