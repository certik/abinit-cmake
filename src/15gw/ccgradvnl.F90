!{\src2tex{textfont=tt}}
!!****f* ABINIT/ccgradvnl
!! NAME
!! ccgradvnl
!!
!! FUNCTION
!! Compute the (grad_K+grad_K') Vnl(K,K') (three reciprocal lattice units components)
!! Needed for chi0(q=0)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2008 ABINIT group (GMR, VO, LR, RWG, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  gvec(3,npwwfn)=integer coordinates of each plane wave in reciprocal space
!!  kibz(3,nkibz)=coordinates of all k points in the irreducible Brillouin zone
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotentials
!!  Cryst<Crystal_structure>=data type gathering information on unit cell and symmetries
!!    %natom=number of atoms
!!    %typat(natom)=type of each atom
!!    %ntypat=number of types of atoms
!!    %xcart(3,natom)=cartesian coordinates of nuclei
!!  nkibz=number of k points in the irreducible Brillouin zone
!!  npwwfn=number of planewaves for wavefunctions
!!  vkb(npwwfn,ntypat,mpsang,nkibz)=KB projector function
!!  vkbd(npwwfn,ntypat,mpsang,nkibz)=derivative of the KB projector function in reciprocal space
!!  vkbsign(mpsang,ntypat)=sign of each KB dyadic product
!!
!! OUTPUT
!!  gradvnl =(grad_K + grad_K') Vnl(K,K') in reciprocal lattice units
!!
!! PARENTS
!!      screening
!!
!! CHILDREN
!!      printcm
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ccgradvnl(npwwfn,nkibz,gvec,gprimd,kibz,Cryst,mpsang,vkbsign,vkb,vkbd,gradvnl)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_15gw, except_this_one => ccgradvnl
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang,nkibz,npwwfn
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npwwfn)
 real(dp),intent(in) :: kibz(3,nkibz),vkb(npwwfn,Cryst%ntypat,mpsang,nkibz)
 real(dp),intent(in) :: vkbd(npwwfn,Cryst%ntypat,mpsang,nkibz),vkbsign(mpsang,Cryst%ntypat)
 real(dp),intent(in) :: gprimd(3,3)
 complex(gwpc),intent(out) :: gradvnl(3,npwwfn,npwwfn,nkibz)

!Local variables ------------------------------
!scalars
 integer :: i,ia,ig,igd1,igd2,igd3,igp,ik,il,is,lmax
 real(dp) :: mkg,mkg2,mkgkgp,mkgp,mkgp2,rgdx,rgdy,rgdz,taugd,x,x2,x3,x4,x5,x6
 real(dp) :: x7,xcostheta
 complex(dpc) :: cs,ct
 character(len=500) :: msg
!arrays
 complex(dpc) :: sfac(Cryst%ntypat)
 integer,parameter :: nlx=8
 real(dp) :: pl(nlx),pld(nlx)
 real(dp) :: kg(3),kgp(3)
 real(dp) :: b1(3),b2(3),b3(3)
 real(dp),pointer :: xcart(:,:)
!************************************************************************

 write(msg,'(a)')' limit q->0, including term <n,k|[Vnl,iqr]|n"k>'
 call wrtout(std_out,msg,'COLL')

 b1=two_pi*gprimd(:,1)
 b2=two_pi*gprimd(:,2)
 b3=two_pi*gprimd(:,3)

 xcart => Cryst%xcart

 lmax=mpsang
 if(mpsang>nlx) then
  write(msg,'(6a)')ch10,&
& ' ccgradvnl : WARNING -',ch10,&
& ' Number of angular momentum components bigger than programmed!!!',ch10,&
& ' taking into account only s p d f g h i j'
  call wrtout(std_out,msg,'COLL')
  lmax=nlx
 end if
 !
 ! Legendre polynomials and their first derivatives
 ! s p d f g h i j  so up to PL_7 = pl(8)
 !
 pl(1)  = one
 pld(1) = zero
 !pl(2) = costheta
 pld(2) = 1.0
 !pl(3) = 1/2 ( 3 costheta**2 - 1 )
 !pld(3)= 3 costheta

 gradvnl(:,:,:,:) = (0.0,0.0)

 do ik=1,nkibz

  do ig=1,npwwfn
   kg(:)=kibz(:,ik) + real(gvec(:,ig))
   mkg2=scpdt(kg,kg,b1,b2,b3)
   mkg=sqrt(mkg2)
   ! The next to solve the problem with k=Gamma.
   if (mkg < 0.0001) cycle

   do igp=1,npwwfn
    kgp(:)=kibz(:,ik) + real(gvec(:,igp))
    mkgp2=scpdt(kgp,kgp,b1,b2,b3)
    mkgp=sqrt(mkgp2)
    ! The next to solve the problem with k=Gamma.
    if (mkgp < 0.0001) cycle

    mkgkgp = mkg*mkgp
    xcostheta = scpdt(kg,kgp,b1,b2,b3) / mkgkgp
    x = xcostheta
    x2 = x * x
    x3 = x2 * x
    x4 = x3 * x
    x5 = x4 * x
    x6 = x5 * x
    x7 = x6 * x
    !
    ! Calculate legendre polynomial PL_0 = pl(1)
    pl(2) = x
    pl(3) = (3.0/2.0) * x2 - (1.0/2.0)
    pl(4) = (5.0/2.0) * x3 - (3.0/2.0) * x
    pl(5) = (35.0/8.0) * x4 - (30.0/8.0) * x2 + (3.0/8.0)
    pl(6) = (63.0/8.0) * x5 - (70.0/8.0) * x3 + (15.0/8.0) * x
    pl(7) = (231.0/16.0) * x6 - (315.0/16.0) * x4 + (105.0/16.0) * x2 - (5.0/16.0)
    pl(8) = (429.0/16.0) * x7 - (693.0/16.0) * x5 + (315.0/16.0) * x3 - (35.0/16.0) * x
    !
    ! Calculate legendre polynomial derivative
    pld(3) = 3.0 * x
    pld(4) = (15.0/2.0) * x2 - (3.0/2.0)
    pld(5) = (35.0/2.0) * x3 - (15.0/2.0) * x
    pld(6) = (315.0/8.0) * x4 - (210.0/8.0) * x2 + (15.0/8.0)
    pld(7) = (693.0/8.0) * x5 - (315.0/4.0) * x3 + (105.0/8.0) * x
    pld(8) = (3003.0/16.0) * x6 - (3465.0/16.0) * x4 + (945.0/16.0) * x2 - (35.0/16.0)

    igd1 = gvec(1,ig)-gvec(1,igp)
    igd2 = gvec(2,ig)-gvec(2,igp)
    igd3 = gvec(3,ig)-gvec(3,igp)
    rgdx = igd1*b1(1)+igd2*b2(1)+igd3*b3(1)
    rgdy = igd1*b1(2)+igd2*b2(2)+igd3*b3(2)
    rgdz = igd1*b1(3)+igd2*b2(3)+igd3*b3(3)

    do is=1,Cryst%ntypat
     sfac(is) = czero
     do ia=1,Cryst%natom
      if (Cryst%typat(ia)/=is) CYCLE
      taugd = rgdx*xcart(1,ia)+rgdy*xcart(2,ia)+ &
&      rgdz*xcart(3,ia)
      sfac(is) = sfac(is) + cmplx(cos(taugd),-sin(taugd))
     end do
    end do

    do i = 1, 3
     gradvnl(i,ig,igp,ik) = 0.0
     do is=1,Cryst%ntypat
      do il = 1, lmax
       ct =( kg(i)*(1/mkgkgp - xcostheta/mkg2 ) + &
&       kgp(i)*(1/mkgkgp - xcostheta/mkgp2 ) ) * &
&       pld(il) * vkbsign(il,is) * vkb(ig,is,il,ik) * vkb(igp,is,il,ik)
       
       cs = pl(il) * vkbsign(il,is) * &
&       ( kg(i)/mkg * vkbd(ig,is,il,ik) * vkb(igp,is,il,ik) + &
&       kgp(i)/mkgp * vkb(ig,is,il,ik) * vkbd(igp,is,il,ik) )
       
       gradvnl(i,ig,igp,ik) = gradvnl(i,ig,igp,ik) + sfac(is) * (ct + cs)
      end do !il
     end do !is
    end do !i

   end do !igp
  end do !ig
 end do !ik

 contains

!!***

!!****f* ccgradvnl/scpdt
!! NAME
!! scpdt
!!
!! FUNCTION
!! Compute scalar product of two vectors u and v in reciprocal space
!!
!! INPUTS
!!  b1(3),b2(3),b3(3)=the three primitive vectors in reciprocal space
!!  u(3),v(3)=the two vectors
!!
!! OUTPUT
!!  function scpdt=scalar product of u and v in reciprocal space
!!
!! SOURCE

function scpdt(u,v,b1,b2,b3)

 use defs_basis

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp) :: scpdt
!arrays
 real(dp),intent(in) :: b1(3),b2(3),b3(3),u(3),v(3)
! *************************************************************************
 scpdt=&
& (u(1)*b1(1)+u(2)*b2(1)+u(3)*b3(1))*(v(1)*b1(1)+v(2)*b2(1)+v(3)*b3(1))+&
& (u(1)*b1(2)+u(2)*b2(2)+u(3)*b3(2))*(v(1)*b1(2)+v(2)*b2(2)+v(3)*b3(2))+&
& (u(1)*b1(3)+u(2)*b2(3)+u(3)*b3(3))*(v(1)*b1(3)+v(2)*b2(3)+v(3)*b3(3))

 end function scpdt

end subroutine ccgradvnl
!!***
