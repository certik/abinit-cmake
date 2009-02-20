        subroutine scramble(i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3,zw,zmpi2)
        use defs_basis

        implicit real(dp) (a-h,o-z)
!Arguments ------------------------------------
        integer :: i1,j2,lot,n1dfft,md1,n3,md2proc,nnd3
        real(dp) :: zw,zmpi2
        dimension zw(2,lot,n3),zmpi2(2,md1,md2proc,nnd3)
!Local variables-------------------------------
! *************************************************************************
        do 100,i3=1,n3
        do 100,i=0,n1dfft-1
        zmpi2(1,i1+i,j2,i3)=zw(1,i+1,i3)
        zmpi2(2,i1+i,j2,i3)=zw(2,i+1,i3)

100     continue

        return

end subroutine scramble
