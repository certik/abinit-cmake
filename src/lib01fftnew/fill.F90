        subroutine fill(nd1,nd3,lot,n1dfft,n3,zf,zw)
        use defs_basis

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: nd1,nd3,lot,n1dfft,n3
        real(dp) ::zw,zf
        dimension zw(2,lot,n3),zf(2,nd1,nd3)
!Local variables-------------------------------

! *************************************************************************
        do 100,i3=1,n3
        do 100,i1=1,n1dfft
        zw(1,i1,i3)=zf(1,i1,i3)
        zw(2,i1,i3)=zf(2,i1,i3)
100     continue

        return

end subroutine fill


