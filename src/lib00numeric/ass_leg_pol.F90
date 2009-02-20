function ass_leg_pol(l,m,xarg)

! Compute the associated Legendre Polynomial Plm(x),
! using a stable recursion formula.
! Here m and l are integers satisfying 0<=m<=l,
! while x lies in the range -1<=x<=1


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none
!Arguments
 integer, intent(in) ::  l,m
 double precision, intent(in) :: xarg
 double precision :: ass_leg_pol
!Local variables
 integer :: i,ll
 double precision :: pll,polmm,tmp1,sqrx,x
 character(len=500) :: message

 x=xarg
 if (m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0) then
  if (m.lt.0.or.m.gt.l.or.abs(x).gt.1.d0+1.d-10) then
   write(message, '(a,a,a,a)' ) char(10),' ass_leg_pol : BUG -',char(10),&
&   '  Bad choice of l, m or x !'
   call wrtout(06,message,'COLL')
   stop
  endif
  x=1.d0
 endif

 polmm=1.d0
 if (m>0) then
  sqrx=sqrt(abs((1.d0-x)*(1.d0+x)))
  do i=1,m
   polmm=polmm*(1.0d0-2.0d0*i)*sqrx
  enddo
 endif

 if (l==m) then
  ass_leg_pol=polmm
 else
  tmp1=x*(2.0d0*m+1.0d0)*polmm
  if (l==(m+1)) then
   ass_leg_pol=tmp1
  else
   do ll=m+2,l
    pll=(x*(2.0d0*ll-1.0d0)*tmp1-(ll+m-1.0d0)*polmm)/dble(ll-m)
    polmm=tmp1
    tmp1=pll
   enddo
   ass_leg_pol=pll
  endif
 endif

end function ass_leg_pol
