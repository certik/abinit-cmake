subroutine smooth(a,mesh,it)


 implicit none
!Arguments
 integer, intent(in) :: it,mesh
 real*8, intent(out) :: a(mesh)
!Local variables
 real*8 :: asm(mesh)
 integer :: i,k

 do k=1,it

 asm(1)=1.0d0/3.0d0*(a(1)+a(2)+a(3))
 asm(2)=0.25d0*(a(1)+a(2)+a(3)+a(4))
 asm(3)=0.2d0*(a(1)+a(2)+a(3)+a(4)+a(5))
 asm(4)=0.2d0*(a(2)+a(3)+a(4)+a(5)+a(6))
 asm(5)=0.2d0*(a(3)+a(4)+a(5)+a(6)+a(7))
 asm(mesh-4)=0.2d0*(a(mesh-2)+a(mesh-3)+a(mesh-4)+&
&                 a(mesh-5)+a(mesh-6))
 asm(mesh-3)=0.2d0*(a(mesh-1)+a(mesh-2)+a(mesh-3)+&
&                 a(mesh-4)+a(mesh-5))
 asm(mesh-2)=0.2d0*(a(mesh)+a(mesh-1)+a(mesh-2)+&
&                 a(mesh-3)+a(mesh-4))
 asm(mesh-1)=0.25d0*(a(mesh)+a(mesh-1)+a(mesh-2)+a(mesh-3))
 asm(mesh)=1.0d0/3.0d0*(a(mesh)+a(mesh-1)+a(mesh-2))

 do i=6,mesh-5
 asm(i)=0.1d0*a(i)+0.1d0*(a(i+1)+a(i-1))+&
&         0.1d0*(a(i+2)+a(i-2))+&
&         0.1d0*(a(i+3)+a(i-3))+&
&         0.1d0*(a(i+4)+a(i-4))+&
&         0.05d0*(a(i+5)+a(i-5))
 enddo

 do i=1,mesh
 a(i)=asm(i)
 enddo

 enddo

end subroutine smooth
