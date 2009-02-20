        subroutine switch(n1dfft,n2,lot,n1,lzt,zt,zw)
        use defs_basis

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: n1dfft,n2,lot,n1,lzt
        real(dp) :: zt,zw
        dimension zw(2,lot,n2),zt(2,lzt,n1)
!Local variables-------------------------------
! *************************************************************************
        do 200,j=1,n1dfft
        do 100,i=1,n2
        zw(1,j,i)=zt(1,i,j)
        zw(2,j,i)=zt(2,i,j)
100     continue
200     continue
        return

end subroutine switch

