 	subroutine unswitchreal(n1dfft,n2,n2eff,lot,n1zt,lzt,zw,zt)
        use defs_basis

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
 integer :: n1dfft,n2,n2eff,lot,n1zt,lzt
 real(dp) :: zw,zt
        dimension zw(2,lot,n2),zt(2,lzt,n1zt)
!Local variables-------------------------------
! *************************************************************************

!       Decompose symmetric and antisymmetric parts
        do j=1,n1dfft
         zt(1,1,2*j-1)=zw(1,j,1)
         zt(2,1,2*j-1)=zero
         zt(1,1,2*j)  =zw(2,j,1)
         zt(2,1,2*j)  =zero
        end do

        do i=2,n2eff
         do j=1,n1dfft
          zt(1,i,2*j-1)= (zw(1,j,i)+zw(1,j,n2+2-i))*half
          zt(2,i,2*j-1)= (zw(2,j,i)-zw(2,j,n2+2-i))*half
          zt(1,i,2*j)  = (zw(2,j,i)+zw(2,j,n2+2-i))*half
          zt(2,i,2*j)  =-(zw(1,j,i)-zw(1,j,n2+2-i))*half
         end do
        end do

end subroutine unswitchreal
