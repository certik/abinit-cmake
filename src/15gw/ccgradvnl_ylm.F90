!{\src2tex{textfont=tt}}
!!****f* ABINIT/ccgradvnl_ylm
!! NAME
!! ccgradvnl_ylm
!!
!! FUNCTION
!!  Compute Vnl(K) and grad_K Vnl(K) three reciprocal lattice units components
!!  using spherical harmonics instead of Legendre polynomials
!!  Needed for chi0(q=0)
!!
!! COPYRIGHT
!!  Copyright (C) 2006-2008 ABINIT group (FB, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Cryst<crystal_structure>=Unit celle and symmetries
!!    %rprimd(3,3)=the three primitive vectors
!!    %ntypat=number of types of atoms
!!    %natom=number of atoms
!!    %xcart(3,natom)=cartesian coordinates of nuclei
!!    %typat(natom)=type of each atom
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space ($\textrm{bohr}^{-1}$)
!!  gvec(3,npwwfn)=integer coordinates of each plane wave in reciprocal space
!!  kibz(3,nkibz)=coordinates of all k points in the irreducible Brillouin zone
!!  mpsang=1+maximum angular momentum for nonlocal pseudopotentials
!!  nkibz=number of k points in the irreducible Brillouin zone
!!  npwwfn=number of planewaves for wavefunctions
!!  vkb(npwwfn,ntypat,mpsang,nkibz)=KB projector function
!!  vkbd(npwwfn,ntypat,mpsang,nkibz)=derivative of the KB projector function in reciprocal space
!!  vkbsign(mpsang,ntypat)=sign of each KB dyadic product
!!
!! OUTPUT
!!  l_fnl(npwwfn,mpsang*mpsang,natom,nkibz),
!!  l_fnld(3,npwwfn,mpsang*mpsang,natom,nkibz)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  Subroutine taken from the EXC code  
!!  All the calculations are done in double precision, but the output arrays l_fnl and l_fnld 
!!  are in single precision, should use double precision after modification of the
!!  other subroutines 
!!  the subroutine does not work wity pseudo with more that one projector per angular state TODO
!!
!! PARENTS
!!  
!! CHILDREN
!!  
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine ccgradvnl_ylm(npwwfn,nkibz,Cryst,gvec,gprimd,kibz,&
& mpsang,vkbsign,vkb,vkbd,l_fnl,l_fnld)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpsang,nkibz,npwwfn
 type(Crystal_structure),intent(in) :: Cryst
!arrays
 integer,intent(in) :: gvec(3,npwwfn)
 real(dp),intent(in) :: gprimd(3,3) 
 real(dp),intent(in) :: kibz(3,nkibz)
 real(dp),intent(in) :: vkb(npwwfn,Cryst%ntypat,mpsang,nkibz)
 real(dp),intent(in) :: vkbd(npwwfn,Cryst%ntypat,mpsang,nkibz),vkbsign(mpsang,Cryst%ntypat)
 complex(gwpc),intent(out) :: l_fnl(npwwfn,mpsang**2,Cryst%natom,nkibz)
 complex(gwpc),intent(out) :: l_fnld(3,npwwfn,mpsang**2,Cryst%natom,nkibz)

!Local variables-------------------------------
!scalars
 integer,parameter :: nlx=4
 integer :: i,ia,ig,ik,il,im,iml,is,lmax
 real(dp),parameter :: ppad=tol8
 real(dp) :: cosphi,costh,factor,mkg,mkg2,sinphi,sinth,sq,xdotg
 complex(dpc) :: dphi,dth,sfac
 character(len=500) :: msg
!arrays
 real(dp) :: gcart(3),kcart(3),kg(3)
 real(dp) :: b1(3),b2(3),b3(3),a1(3),a2(3),a3(3)
 real(dp),pointer :: xcart(:,:)
 complex(dpc) :: dylmcart(3),dylmcrys(3),gradphi(3),gradth(3)
!************************************************************************

#if defined DEBUG_MODE
 write(msg,'(a)')' ccgradvnl_ylm : enter'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

 write(msg,'(a)')&
& ' ccgradvnl_ylm : limit q->0, including term <n,k|[Vnl,iqr]|n"k> using Y_lm'
 call wrtout(std_out,msg,'COLL')

 lmax=mpsang
 if (mpsang>nlx) then
  write(msg,'(6a)')ch10,&
&  ' ccgradvnl_ylm:  WARNING- ',ch10,&
&  '  number of angular momentum components bigger than programmed.',ch10,&
&  '  taking into account only s p d f ' 
  call wrtout(std_out,msg,'COLL')
  lmax=nlx
 end if

 a1=Cryst%rprimd(:,1) ; b1=two_pi*gprimd(:,1)
 a2=Cryst%rprimd(:,2) ; b2=two_pi*gprimd(:,2)
 a3=Cryst%rprimd(:,3) ; b3=two_pi*gprimd(:,3)

 xcart => Cryst%xcart
 !
 ! === Calculate Kleiman-Bylander factor and first derivative ===
 l_fnl=(0.0_gwp,0.0_gwp) ; l_fnld=(0.0_gwp,0.0_gwp)
 do ik=1,nkibz
  do ig=1,npwwfn
   !  
   ! === Get kcart =k+G in Cartesian coordinates ===
   kg(:)= kibz(:,ik)+REAL(gvec(:,ig))
   kcart(:)=kg(1)*b1(:)+kg(2)*b2(:)+kg(3)*b3(:)
   ! * Solve the problem with sinth=0. or sinphi=0
   if (ABS(kcart(2))<ppad) kcart(2)=kcart(2)+ppad

   mkg2=kcart(1)**2+kcart(2)**2+kcart(3)**2
   mkg=SQRT(mkg2)
   sq=SQRT(kcart(1)**2+kcart(2)**2)

   gcart(:)=  REAL(gvec(1,ig))*b1(:)&
&            +REAL(gvec(2,ig))*b2(:)&
&            +REAL(gvec(3,ig))*b3(:)

   ! === Calculate spherical coordinates (th,phi) ===
   costh = kcart(3)/mkg
   sinth = sq/mkg
   cosphi= kcart(1)/sq
   sinphi= kcart(2)/sq
   
   gradth(1)  = kcart(1)*kcart(3)/mkg**3/sinth
   gradth(2)  = kcart(2)*kcart(3)/mkg**3/sinth
   gradth(3)  = -(one/mkg-kcart(3)**2/mkg**3)/sinth
   gradphi(1) = -(one/sq - kcart(1)**2/sq**3)/sinphi
   gradphi(2) = kcart(2)*kcart(1)/sq**3/sinphi
   gradphi(3) = czero
   
   do ia=1,Cryst%natom
    is=Cryst%typat(ia)
    xdotg=gcart(1)*xcart(1,ia)+gcart(2)*xcart(2,ia)+gcart(3)*xcart(3,ia)
    ! Remember that in the GW code the reciprocal vectors 
    ! are defined such as a_i*b_j = 2pi delta_ij, no need to introduce 2pi
    sfac=CMPLX(COS(xdotg),SIN(xdotg)) 

    do il=1,lmax
     factor=SQRT(four_pi/REAL(2*(il-1)+1))
     do im= 1,2*(il-1)+1
      ! Index of im and il
      iml=im+(il-1)*(il-1)
      !     
      ! Calculate the first KB factor, note that l_fnl is simple precision complex
      l_fnl(ig,iml,ia,ik)=factor*sfac*ylmc(il-1,im-il,kcart)*vkb(ig,is,il,ik)*vkbsign(il,is)
      !     
      ! Calculate the second KB factor (involving first derivatives)
      ! dYlm/dK = dYlm/dth * grad_K th + dYlm/dphi + grad_K phi
      call ylmcd(il-1,im-il,kcart,dth,dphi)
      dylmcart(:) = dth*gradth(:) + dphi*gradphi(:)
      !     
      ! Cartesian to crystallographic axis
      dylmcrys(:) = ( dylmcart(1)*a1(:)&
&                    +dylmcart(2)*a2(:)&
&                    +dylmcart(3)*a3(:) ) /(two_pi)
      
      ! Note that l_fnld is simple precision complex, it could be possible use double precision
      do i=1,3
       l_fnld(i,ig,iml,ia,ik) = factor*sfac*&
&       ( kg(i)/mkg*ylmc(il-1,im-il,kcart)*vkbd(ig,is,il,ik) + dylmcrys(i)*vkb(ig,is,il,ik) )
      end do 

     end do !im
    end do !il
   end do !ia
  end do !ig
 end do !ik

#if defined DEBUG_MODE
 write(msg,'(a)')' ccgradvnl_ylm : exit'
 call wrtout(std_out,msg,'COLL')
 call flush_unit(std_out)
#endif

end subroutine ccgradvnl_ylm
!!***
