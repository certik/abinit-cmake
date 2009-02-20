!{\src2tex{textfont=tt}}
!!****f* ABINIT/volumeint
!! NAME
!! volumeint
!!
!! FUNCTION
!! Computes the values within a volume
!!
!! COPYRIGHT
!! Copyright (C) 2000-2008 ABINIT group (GMR, RC, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! grid(nr1,nr2,nr3)=grid matrix
!! natom=integer number of atoms
!! nr1=grid size along x
!! nr2=grid size along y
!! nr3=grid size along z
!! nspden=number of spin-density components
!! rprimd(3,3)=orientation of the unit cell in 3D
!! tau(3,natom)=list of atoms in cartesian coordinates
!!
!! OUTPUT
!!  only writing
!!
!! PARENTS
!!      cut3d
!!
!! CHILDREN
!!      interpol3d,matr3inv,normalize,recip,reduce,vdot
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine volumeint(grid,gridtt,gridux,griddy,gridmz,natom,nr1,nr2,nr3,nspden,rprimd,tau)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_11util
 use interfaces_19cut3d, except_this_one => volumeint
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nr1,nr2,nr3,nspden
!arrays
 real(dp),intent(in) :: grid(nr1,nr2,nr3),griddy(nr1,nr2,nr3)
 real(dp),intent(in) :: gridmz(nr1,nr2,nr3),gridtt(nr1,nr2,nr3)
 real(dp),intent(in) :: gridux(nr1,nr2,nr3),rprimd(3,3),tau(3,natom)

!Local variables -------------------------
!scalars
 integer :: fileformattype,iat,idir,ii,inpopt,itypat,jj,k1,k2,k3,mu,nresolh
 integer :: nresoll,nresolw,okhkl,okinp,okparam,oksure,planetype
 integer :: referenceposition
 real(dp) :: coef,denval,denvaldy,denvalmz,denvaltt,denvalux,determ,height
 real(dp) :: length,width,xcoord,ycoord,zcoord
 real*4 :: xm,xp,ym,yp,zm,zp
 character(len=fnlen) :: filnam
!arrays
 integer :: hkl(3)
 real(dp) :: cent(3),centpl(3),mm(3,3),mminv(3,3),r1(3),r2(3),r3(3),rcart(3)
 real(dp) :: rr(3),x1(3),x2(3),x3(3),xcart(3)
 real*4,allocatable :: rhomacudy(:,:),rhomacumz(:,:),rhomacutt(:,:)
 real*4,allocatable :: rhomacuux(:,:)

! *********************************************************************

 call matr3inv(rprimd,mminv)
!Start of the real input of the volume orientation

 write(*,*)
 write(*,*) ' The volume is an orthogonal prism, that is defined by: '
 write(*,*) ' the basal plane and'
 write(*,*) ' the height perpendicular to the basal plane'
 write(*,*)
 write(*,*) ' First you will define the basal plane '
 write(*,*) ' second you will define the height'
 write(*,*) ' and third you will define the basal plane position '
 write(*,*) ' along the height vector'

 do
  write(*,*)
  write(*,*) '  Type 1) for a plane passing through 3 atoms'
  write(*,*) '    or 2) for a plane passing through 3 cartesian points'
  write(*,*) '    or 3) for a plane passing through 3 crystallographic points'
  write(*,*) '    or 4) for a plane parallel to a crystallographic plane'
  write(*,*) '    or 5) for a plane orthogonal to a cartesian direction'
  write(*,*) '    or 6) for a plane orthogonal to a crystallographic direction'
  write(*,*) '    or 0) to stop'
  read(*,*) itypat

  select case (itypat)

   case (0)
    stop

!   A plane passing through 3 atoms
   case (1)
    write(*,*) '  The X axis will be through atoms: 1,2 '
    write(*,*) '  Define each atom by its species and its number:'
    write(*,*) '    -> atom 1 (iat):'
    read(*,*) iat
    x1(1)=tau(1,iat)
    x1(2)=tau(2,iat)
    x1(3)=tau(3,iat)
    write(*,'(a,3f10.6)') '        position: ',x1
    write(*,*)
    write(*,*) '    -> atom 2 (iat):'
    read(*,*) iat
    x2(1)=tau(1,iat)
    x2(2)=tau(2,iat)
    x2(3)=tau(3,iat)
    write(*,'(a,3f10.6)') '        position: ',x2
    write(*,*)
    write(*,*) '    -> atom 3 (iat):'
    read(*,*) iat
    x3(1)=tau(1,iat)
    x3(2)=tau(2,iat)
    x3(3)=tau(3,iat)
    write(*,'(a,3f10.6)') '        position: ',x3
    write(*,*)

!   Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
    do idir=1,3
     x2(idir)=x2(idir)-x1(idir)
     x3(idir)=x3(idir)-x1(idir)
    end do
    call normalize(x2)
    call vdot(x3,x2,x1)
    call normalize(x1)
    call vdot(x2,x1,x3)
    call normalize(x3)
    exit

!   A plane passing through 3 cartesian points
   case (2)
    write(*,*) '  The X axis will be through points: 1,2 '
    write(*,*) '  Define each :point coordinates'
    write(*,*) '    -> point 1:    X-coord  Y-coord  Z-coord:'
    read(*,*) xcart
    x1(:)=xcart(:)
    write(*,'(a,3f10.6)') ' crystallographic position: ',x1
    write(*,*)
    write(*,*) '    -> point 2:    X-coord  Y-coord  Z-coord:'
    read(*,*) xcart
    x2(:)=xcart(:)
    write(*,'(a,3f10.6)') ' crystallographic position: ',x2
    write(*,*)
    write(*,*) '    -> point 3:    X-coord  Y-coord  Z-coord:'
    read(*,*) xcart
    x3(:)=xcart(:)
    write(*,'(a,3f10.6)') ' crystallographic position: ',x3
    write(*,*)

!   Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
    do idir=1,3
     x2(idir)=x2(idir)-x1(idir)
     x3(idir)=x3(idir)-x1(idir)
    end do
    call normalize(x2)
    call vdot(x3,x2,x1)
    call normalize(x1)
    call vdot(x2,x1,x3)
    call normalize(x3)
    exit

!   A plane passing through 3 crystallographic points
   case (3)
    write(*,*) '  The X axis will be through points: 1,2 '
    write(*,*) '  Define each :point coordinates'
    write(*,*) '    -> point 1:    X-coord  Y-coord  Z-coord:'
    read(*,*) r1
    write(*,'(a,3f10.6)') ' crystallographic position: ',r1
    write(*,*)
    write(*,*) '    -> point 2:    X-coord  Y-coord  Z-coord:'
    read(*,*) r2
    write(*,'(a,3f10.6)') ' crystallographic position: ',r2
    write(*,*)
    write(*,*) '    -> point 3:    X-coord  Y-coord  Z-coord:'
    read(*,*) r3
    write(*,'(a,3f10.6)') ' crystallographic position: ',r3
    write(*,*)

!   Transforms the points coordinates into cartesian
    do mu=1,3
     x1(mu)=rprimd(mu,1)*r1(1)+rprimd(mu,2)*r1(2)+rprimd(mu,3)*r1(3)
     x2(mu)=rprimd(mu,1)*r2(1)+rprimd(mu,2)*r2(2)+rprimd(mu,3)*r2(3)
     x3(mu)=rprimd(mu,1)*r3(1)+rprimd(mu,2)*r3(2)+rprimd(mu,3)*r3(3)
    end do

    write(*,*) ' Cartesian positions:'
    write(*,*) x1
    write(*,*) x2
    write(*,*) x3

!   Compute the 3 orthogonal normalized vectors from x2-x1, x3-x1
    do idir=1,3
     x2(idir)=x2(idir)-x1(idir)
     x3(idir)=x3(idir)-x1(idir)
    end do
    call normalize(x2)
    call vdot(x3,x2,x1)
    call normalize(x1)
    call vdot(x2,x1,x3)
    call normalize(x3)
    exit

!   A plane parallel to a crystallographic plane
   case (4)
    okhkl=0
    do while (okhkl==0)
     write(*,*) '  Enter plane coordinates:'
     write(*,*) '    ->H  K  L '
     read(*,*) hkl
     if (.not. (hkl(1)==0 .and. hkl(2)==0 .and. hkl(3)==0)) okhkl=1
    end do
    write(*,*) ' Miller indices are:',hkl

    call recip(x1,hkl,rprimd)
    write(*,*) ' Orthogonal vector to the plane',x1

    call normalize(x1)
    if((x1(1).ne.0).or.(x1(2).ne.0)) then
     x2(1)=-x1(2)
     x2(2)=x1(1)
     x2(3)=0
     call normalize(x2)
    else
     x2(1)=1
     x2(2)=0
     x2(3)=0
    end if
    call vdot(x2,x1,x3)
    call normalize(x3)
    exit

!   A plane orthogonal to a cartesian direction
   case (5)
    write(*,*) '  Enter cartesian vector orthogonal to plane:'
    write(*,*) '    -> X-dir   Y-dir   Z-dir (Angstroms or Bohrs):'
    read(*,*) x1
    call normalize(x1)
    if((x1(1).ne.0).or.(x1(2).ne.0)) then
     x2(1)=-x1(2)
     x2(2)=x1(1)
     x2(3)=0
     call normalize(x2)
    else
     x2(1)=1
     x2(2)=0
     x2(3)=0
    end if
    call vdot(x2,x1,x3)
    call normalize(x3)
    exit

!   A plane orthogonal to a crystallographic direction
   case (6)
    write(*,*) '  Enter crystallographic vector orthogonal to plane:'
    write(*,*) '    -> X-dir   Y-dir   Z-dir (Fractional coordinates):'
    read(*,*) r1
    do mu=1,3
     x1(mu)=rprimd(mu,1)*r1(1)+rprimd(mu,2)*r1(2)+rprimd(mu,3)*r1(3)
    end do
    call normalize(x1)
    if((x1(1).ne.0).or.(x1(2).ne.0)) then
     x2(1)=-x1(2)
     x2(2)=x1(1)
     x2(3)=0
     call normalize(x2)
    else
     x2(1)=1
     x2(2)=0
     x2(3)=0
    end if
    call vdot(x2,x1,x3)
    call normalize(x3)
    exit

    case default
    write(*,*) ' Input option does not correspond to one available option'
    write(*,*) ' Please try again'
    cycle

  end select
 end do

!At this moment the family of planes was defined
!The code knows also some of the geometric input
!It will proceed to the anchorage of the plane onto a point and then
!to the effective calculation

 write(*,*) '  Vectors: (orthogonal & normalized)   '
 write(*,'(11x,a,3f10.6)') '  X-dir in the plot         ',x2
 write(*,'(11x,a,3f10.6)') '  Y-dir in the plot         ',x3
 write(*,'(11x,a,3f10.6)') '  Z-dir (orth. to the plot) ',x1

 do
  write(*,*)
  write(*,*) '  Enter reference point of plane (Bohr):'
  write(*,*) '  Type 1) for Cartesian coordinates.'
  write(*,*) '    or 2) for Crystallographic coordinates.'
  read(*,*) inpopt

  select case (inpopt)

   case(1)
    write(*,*) '    -> X-Coord   Y-Coord   Z-Coord:'
    read(*,*) cent
    exit
   case(2)
    write(*,*) '    -> X-Coord   Y-Coord   Z-Coord:'
    read(*,*) cent
    do mu=1,3
     rcart(mu)=rprimd(mu,1)*cent(1)+rprimd(mu,2)*cent(2)+rprimd(mu,3)*cent(3)
    end do
    cent(:)=rcart(:)
    write(*, '(a,3es16.6)' ) ' Expressed in cartesian coordinates: ',cent(1:3)
    exit
   case (3)
    cycle

  end select
 end do

!End of basal plane orientation

!Input box dimensions now

 write(*,*)
 write(*,*) ' It is now time to input the 3D box dimensions.'
 write(*,*) ' and the position of the basal plane in the box.'

 do
  write(*,*)
  write(*,*) '  Enter in-plane width:'
  read(*,*) width
  write(*,*) '  Enter in-plane length:'
  read(*,*) length
  write(*,*) '  Enter box height:'
  read(*,*) height
  write(*,*)
  write(*,*) ' Enter the position of the basal plane in the box:'
  do
   write(*,*)
   write(*,*) ' Type 1) for DOWN'
   write(*,*) ' Type 2) for MIDDLE'
   write(*,*) ' Type 3) for UP'
   read(*,*) planetype

   select case(planetype)

    case (1)
     exit
    case (2)
     exit
    case (3)
     exit
     case default
     cycle

   end select
  end do

  write(*,*) ' Enter the position of the reference point in the basal plane '

  do
   write(*,*)
   write(*,*) ' Type 1) for CENTRAL position '
   write(*,*) ' Type 2) for CORNER(0,0) position '
   read(*,*) referenceposition

   select case(referenceposition)

    case (1)
     exit
    case (2)
     exit
     case default
     cycle

   end select
  end do

  write(*,*)
  write(*,*) ' Enter the box grid values:'
  write(*,*) '  Enter plane resolution in width:'
  read(*,*) nresolw
  write(*,*) '  Enter plane resolution in lenth:'
  read(*,*) nresoll
  write(*,*) '  Enter height resolution:'
  read(*,*) nresolh
  write(*,*)
  write(*,*) ch10,'  Enter the name of an output file:'
  read(*,*) filnam
  write(*,*) '  The name of your file is : ',trim(filnam)

  do
   write(*,*) '  Enter the format of the output file:'
   write(*,*) '   Type 1=> ASCII formatted'
   write(*,*) '   Type 2=> Molekel formatted'
   read(*,*) fileformattype
   if (fileformattype==1 .or. fileformattype==2) then
    exit
   else
    cycle
   end if
  end do

  write(*,*) ' You asked for a 3d box of:'
  write(*,*) length,' x ',width,' x ',height
  write(*,*) ' With a resolution of ;'
  write(*,*) nresoll,' x ',nresolw,' x ',nresolh
  write(*,*) ' The result will be redirected to the file:  ',trim(filnam)
  if (fileformattype==1) then
   write(*,*) ' ASCII formatted'
  else
   write(*,*) ' Molekel formatted'
  end if
  write(*,*) ' These parameters may still be changed.'
  write(*,*) ' Are you sure you want to keep them? (1=default=yes,2=no) '
  read(*,*) okparam
  if (okparam==2) then
   cycle
  else
   exit
  end if
 end do

!Write the header of the Molekel input file

 if (fileformattype==1) then
  open(unit=31,file=trim(filnam),status='unknown')
 else
  open(unit=31,file=trim(filnam),form='unformatted')

  xm=0
  xp=length
  ym=0
  yp=width
  zm=0
  zp=height

  write(6, '(a)' )&
&  ' Extremas of the cube in which the system is placed (x,y,z), in Angs.:'
  write(6, '(5x,6f10.5)' ) xm,xp,ym,yp,zm,zp
  write(6, '(a,a,3i5)' ) ch10,&
&  ' Number of points per side:   ',nresolw+1,nresoll+1,nresolh+1
  write(6, '(a,a,i10,a,a)' ) ch10,&
&  ' Total number of points:  ',(nresolw+1)*(nresoll+1)*(nresolh+1),&
&  ch10,ch10

  write(31) xm,xp,ym,yp,zm,zp,nresolw+1,nresoll+1,nresolh+1

 end if

!Allocate rhomacu in case of molekel output format
 allocate(rhomacutt(nresoll+1,nresolw+1))
 allocate(rhomacuux(nresoll+1,nresolw+1))
 allocate(rhomacudy(nresoll+1,nresolw+1))
 allocate(rhomacumz(nresoll+1,nresolw+1))

 do k1=0,nresolh

  select case (planetype)

!  Basal plane at the bottom
   case (1)
    centpl(1)=cent(1)+k1*x1(1)*height/nresolh
    centpl(2)=cent(2)+k1*x1(2)*height/nresolh
    centpl(3)=cent(3)+k1*x1(3)*height/nresolh

!   Basal plane in the middle
   case (2)
    centpl(1)=cent(1)+(k1-nresolh/2)*x1(1)*height/nresolh
    centpl(2)=cent(2)+(k1-nresolh/2)*x1(2)*height/nresolh
    centpl(3)=cent(3)+(k1-nresolh/2)*x1(3)*height/nresolh

!   Basal plane on the top
   case (3)
    centpl(1)=cent(1)+(k1-nresolh)*x1(1)*height/nresolh
    centpl(2)=cent(2)+(k1-nresolh)*x1(2)*height/nresolh
    centpl(3)=cent(3)+(k1-nresolh)*x1(3)*height/nresolh

  end select

  do k3=0,nresolw
   do k2=0,nresoll

    select case(referenceposition)

!    Reference point in the middle of the basal plane
     case (1)
      rcart(1)=centpl(1) + (k2-nresoll/2)*x2(1)*length/nresoll + (k3-nresolw/2)*x3(1)*width/nresolw
      rcart(2)=centpl(2) + (k2-nresoll/2)*x2(2)*length/nresoll + (k3-nresolw/2)*x3(2)*width/nresolw
      rcart(3)=centpl(3) + (k2-nresoll/2)*x2(3)*length/nresoll + (k3-nresolw/2)*x3(3)*width/nresolw

!     Reference point in the corner of the basal plane
     case (2)
      rcart(1)=centpl(1) + k2*x2(1)*length/nresoll + k3*x3(1)*width/nresolw
      rcart(2)=centpl(2) + k2*x2(2)*length/nresoll + k3*x3(2)*width/nresolw
      rcart(3)=centpl(3) + k2*x2(3)*length/nresoll + k3*x3(3)*width/nresolw

    end select

    call reduce(rr,rcart,rprimd)
    rr(1)=mod(mod(rr(1),1._dp)+1._dp,1._dp)
    rr(2)=mod(mod(rr(2),1._dp)+1._dp,1._dp)
    rr(3)=mod(mod(rr(3),1._dp)+1._dp,1._dp)

    call interpol3d(rr,nr1,nr2,nr3,denvaltt,gridtt)
    if(nspden==2 .or. nspden==4)then
     call interpol3d(rr,nr1,nr2,nr3,denvalux,gridux)
     call interpol3d(rr,nr1,nr2,nr3,denvaldy,griddy)
     call interpol3d(rr,nr1,nr2,nr3,denvalmz,gridmz)
    end if

    if (fileformattype==1) then
     if(nspden==1)then
      write(31, '(es22.12)' ) denvaltt
     else if(nspden==2 .or. nspden==4)then
      write(31, '(4(es22.12))' ) denvaltt,denvalux,denvaldy,denvalmz
     end if

    else
     rhomacutt(k2+1,k3+1)=denvaltt
     if(nspden==2 .or. nspden==4)then
      rhomacuux(k2+1,k3+1)=denvalux
      rhomacudy(k2+1,k3+1)=denvaldy
      rhomacumz(k2+1,k3+1)=denvalmz
     end if
    end if

   end do
  end do

  if (fileformattype==2) then
   write(31) rhomacutt(:,:)
   if(nspden==2 .or. nspden==4)then
    write(31) rhomacuux(:,:)
    write(31) rhomacudy(:,:)
    write(31) rhomacumz(:,:)
   end if
  end if

 end do

 close(31)
 
 deallocate(rhomacutt,rhomacuux,rhomacudy,rhomacumz)

end subroutine volumeint
!!***
