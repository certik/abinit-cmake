subroutine unswitch(n1dfft,n2,lot,n1,lzt,zw,zt)
  use defs_basis

  implicit real(dp) (a-h,o-z)
  !Arguments ------------------------------------
  integer :: n1dfft,n2,lot,n1,lzt
  real(dp) :: zw,zt
  dimension zw(2,lot,n2),zt(2,lzt,n1)
  !Local variables-------------------------------
  ! *************************************************************************
  do 100,j=1,n1dfft
     do 100,i=1,n2
        zt(1,i,j)=zw(1,j,i)
        zt(2,i,j)=zw(2,j,i)
100     continue
  return
end subroutine unswitch
