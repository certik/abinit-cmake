        subroutine switch_cent(n1dfft,max2,m2,n2,lot,n1,lzt,zt,zw)
        use defs_basis

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: n1dfft,max2,m2,n2,lot,n1,lzt
        real(dp):: zt,zw
        dimension zw(2,lot,n2),zt(2,lzt,n1)
!Local variables-------------------------------
! *************************************************************************
!       Here, zero and positive frequencies
        do 90,j=1,n1dfft
        do 90,i=1,max2+1
        zw(1,j,i)=zt(1,i,j)
        zw(2,j,i)=zt(2,i,j)
90      continue

!       Fill the center region with zeros
        do 100,i=max2+2,n2-m2+max2+1
        do 100,j=1,n1dfft
        zw(1,j,i)=zero
        zw(2,j,i)=zero
100     continue

!       Here, negative frequencies
        if(m2>=max2+2)then
         do 110,j=1,n1dfft
         do 110,i=max2+2,m2
         zw(1,j,i+n2-m2)=zt(1,i,j)
         zw(2,j,i+n2-m2)=zt(2,i,j)
110      continue
        end if

        return

end subroutine switch_cent

