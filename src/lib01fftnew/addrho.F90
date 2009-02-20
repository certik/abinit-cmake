        subroutine addrho(icplexwf,includelast,nd1,nd2,n2,lot,n1dfft,zw,rhopart,weight)
        use defs_basis

        implicit real(dp) (a-h,o-z)

!Arguments ------------------------------------
        integer :: icplexwf,includelast,nd1,nd2,n2,lot,n1dfft
        dimension zw(2,lot,n2),rhopart(nd1,nd2)
        real(dp) :: weight
!Local variables-------------------------------

! *************************************************************************
        if(icplexwf==2)then

         do i2=1,n2-1,2
          do j=1,n1dfft
           rhopart(j,i2+0)=rhopart(j,i2+0)+(zw(1,j,i2+0)**2+zw(2,j,i2+0)**2)*weight
           rhopart(j,i2+1)=rhopart(j,i2+1)+(zw(1,j,i2+1)**2+zw(2,j,i2+1)**2)*weight
          end do
         end do
         if(2*(n2/2)/=n2)then
          do j=1,n1dfft
           rhopart(j,n2  )=rhopart(j,n2  )+(zw(1,j,n2  )**2+zw(2,j,n2  )**2)*weight
          end do
         end if

        else

!        The wavefunction is real, in real space
         if(includelast==1)then
          do i2=1,n2
           do j=1,n1dfft
            rhopart(2*j-1,i2)=rhopart(2*j-1,i2)+zw(1,j,i2)**2*weight
            rhopart(2*j  ,i2)=rhopart(2*j  ,i2)+zw(2,j,i2)**2*weight
           end do
          end do
         else
          do i2=1,n2
           do j=1,n1dfft-1
            rhopart(2*j-1,i2)=rhopart(2*j-1,i2)+zw(1,j,i2)**2*weight
            rhopart(2*j  ,i2)=rhopart(2*j  ,i2)+zw(2,j,i2)**2*weight
           end do
           rhopart(2*n1dfft-1,i2)=rhopart(2*n1dfft-1,i2)+zw(1,n1dfft,i2)**2*weight
          end do
         end if

        end if

        return

end subroutine addrho

