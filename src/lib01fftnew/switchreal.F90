 	subroutine switchreal(includelast,n1dfft,n2,n2eff,lot,n1zt,lzt,zt,zw)
        use defs_basis

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
         integer :: includelast,n1dfft,n2,n2eff,lot,n1zt,lzt
         real(dp) :: zt,zw
        dimension zw(2,lot,n2),zt(2,lzt,n1zt)
!Local variables-------------------------------
! *************************************************************************
        if(includelast==1)then

!        Compute symmetric and antisymmetric combinations
         do j=1,n1dfft
          zw(1,j,1)=zt(1,1,2*j-1)
          zw(2,j,1)=zt(1,1,2*j  )
         end do
         do i=2,n2eff
          do j=1,n1dfft
           zw(1,j,i)=      zt(1,i,2*j-1)-zt(2,i,2*j)
           zw(2,j,i)=      zt(2,i,2*j-1)+zt(1,i,2*j)
           zw(1,j,n2+2-i)= zt(1,i,2*j-1)+zt(2,i,2*j)
           zw(2,j,n2+2-i)=-zt(2,i,2*j-1)+zt(1,i,2*j)
          end do
         end do

        else

!        An odd number of FFTs

!        Compute symmetric and antisymmetric combinations
         do j=1,n1dfft-1
          zw(1,j,1)=zt(1,1,2*j-1)
          zw(2,j,1)=zt(1,1,2*j  )
         end do
         zw(1,n1dfft,1)=zt(1,1,2*n1dfft-1)
         zw(2,n1dfft,1)=zero
         do i=2,n2eff
          do j=1,n1dfft-1
           zw(1,j,i)=      zt(1,i,2*j-1)-zt(2,i,2*j)
           zw(2,j,i)=      zt(2,i,2*j-1)+zt(1,i,2*j)
           zw(1,j,n2+2-i)= zt(1,i,2*j-1)+zt(2,i,2*j)
           zw(2,j,n2+2-i)=-zt(2,i,2*j-1)+zt(1,i,2*j)
          end do
          zw(1,n1dfft,i)=      zt(1,i,2*n1dfft-1)
          zw(2,n1dfft,i)=      zt(2,i,2*n1dfft-1)
          zw(1,n1dfft,n2+2-i)= zt(1,i,2*n1dfft-1)
          zw(2,n1dfft,n2+2-i)=-zt(2,i,2*n1dfft-1)
         end do

        end if

end subroutine switchreal
