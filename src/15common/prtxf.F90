!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtxf
!! NAME
!! prtxf
!!
!!
!! FUNCTION
!! Compute and print out dimensional cartesian coordinates and forces.
!! Note: for x=cartesian coordinates, t=reduced coordinates (xred), =>
!!  $ x= R t => x(1)=rprimd(1,1) t(1)+rprimd(2,1) t(2)+rprimd(3,1) t(3)$
!!  etc. Also $ t = (R^{-1}) x$ .
!!  To convert gradients, $d(E)/dx(n) = [d(E)/dt(m)] [dt(m)/dx(n)]$
!!  and $ dt(m)/dx(n) = (R^{-1})_{mn} = G_{nm}$ because G is the inverse transpose
!!  of R.  Finally then   $d(E)/dx(n) = G_{nm} [d(E)/dt(m)]$.  The
!!  vector $d(E)/dt(m)$ for each atom is input in fred (grad. wrt xred).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  fred(3,natom)=gradients of Etot (hartree) wrt xred(3,natom)
!!  iatfix(3,natom)=1 for each fixed atom along specified direction, else 0
!!  iout=unit number for output file
!!  iwfrc=controls force output: 0=> no force output
!!  natom=number of atoms in unit cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  xred(3,natom)=relative coordinates of atoms (in terms of prim. transl.)
!!
!! OUTPUT
!!  (data written to unit iout)
!!
!! PARENTS
!!      clnup1
!!
!! CHILDREN
!!      matr3inv,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine prtxf(fred,iatfix,iout,iwfrc,natom,rprimd,xred)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,iwfrc,natom
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: fred(3,natom),rprimd(3,3),xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu,unfixd
 real(dp) :: convt,fmax,frms
 character(len=500) :: message
!arrays
 real(dp) :: favg(3),favg_out(3),ff(3),gprimd(3,3),xx(3)

! *************************************************************************

!Write cartesian coordinates in angstroms
 write(message, '(a)' ) ' cartesian coordinates (angstrom) at end:'
 call wrtout(iout,message,'COLL')
 do iatom=1,natom
  do mu=1,3
   xx(mu)=(rprimd(mu,1)*xred(1,iatom)+&
&   rprimd(mu,2)*xred(2,iatom)+&
&   rprimd(mu,3)*xred(3,iatom))*Bohr_Ang
  end do
  write(message, '(i5,1x,3f21.14)' ) iatom,xx
  call wrtout(iout,message,'COLL')
 end do

!Optionally write cartesian forces in eV/Angstrom
!(also provide same in hartree/bohr)
 if (iwfrc/=0) then
! First, provide results in hartree/bohr
  write(message, '(a,a)' ) ch10,&
&  ' cartesian forces (hartree/bohr) at end:'
  call wrtout(iout,message,'COLL')
  frms=zero
  fmax=zero
  favg(1)=zero
  favg(2)=zero
  favg(3)=zero
! To get cartesian forces from input gradients with respect to
! dimensionless coordinates xred, multiply by G and negate
! (see notes at top of this subroutine)
  call matr3inv(rprimd,gprimd)
! First compute (spurious) average force favg
  do iatom=1,natom
   do mu=1,3
    ff(mu)=-(gprimd(mu,1)*fred(1,iatom)+&
&    gprimd(mu,2)*fred(2,iatom)+&
&    gprimd(mu,3)*fred(3,iatom))
    favg(mu)=favg(mu)+ff(mu)
   end do
  end do
  favg(1) = favg(1)/dble(natom)
  favg(2) = favg(2)/dble(natom)
  favg(3) = favg(3)/dble(natom)

! Subtract off average force in what follows
! (avg is also subtracted off in carfor, called by loopcv,
! called by grad)
  unfixd=0
  do iatom=1,natom
   do mu=1,3
    ff(mu)=-(gprimd(mu,1)*fred(1,iatom)+&
&    gprimd(mu,2)*fred(2,iatom)+&
&    gprimd(mu,3)*fred(3,iatom))-favg(mu)
!   For rms and max force, include only unfixed components
    if (iatfix(mu,iatom) /= 1) then
     unfixd=unfixd+1
     frms=frms+ff(mu)**2
     fmax=max(fmax,abs(ff(mu)))
    end if
   end do
   write(message, '(i5,1x,3f21.14)' ) iatom,ff
   call wrtout(iout,message,'COLL')
  end do
  if ( unfixd /= 0 ) frms = sqrt(frms/dble(unfixd))

! The average force is obtained from the cancellation of numbers
! of typical size unity, so an absolute value lower
! than tol14 is meaningless for the output file.
  favg_out(:)=favg(:)
  if(abs(favg_out(1))<tol14)favg_out(1)=zero
  if(abs(favg_out(2))<tol14)favg_out(2)=zero
  if(abs(favg_out(3))<tol14)favg_out(3)=zero

  write(message, '(a,1p,2e14.7,1x,3e11.3,a)' ) &
&  ' frms,max,avg=',frms,fmax,favg_out(1:3),' h/b'
  call wrtout(iout,message,'COLL')

  write(message, '(a,a)' ) ch10,&
&  ' cartesian forces (eV/Angstrom) at end:'
  call wrtout(iout,message,'COLL')
  convt=Ha_eV/Bohr_Ang

! Note: subtract off average force
  do iatom=1,natom
   do mu=1,3
    ff(mu)=(-(gprimd(mu,1)*fred(1,iatom)+&
&    gprimd(mu,2)*fred(2,iatom)+&
&    gprimd(mu,3)*fred(3,iatom))-favg(mu))*convt
   end do
   write(message, '(i5,1x,3f21.14)' ) iatom,ff
   call wrtout(iout,message,'COLL')
  end do
  write(message, '(a,1p,2e14.7,1x,3e11.3,a)' ) &
&  ' frms,max,avg=',convt*frms,convt*fmax,&
&  convt*favg_out(1:3),' e/A'
  call wrtout(iout,message,'COLL')

! DEBUG
! write(message, '(a,a)' ) ch10,&
! &   ' reduced gradients (hartree) at end:'
! call wrtout(iout,message,'COLL')
! do iatom=1,natom
! write(message, '(i5,1x,3f21.14)' ) iatom,fred(:,iatom)
! call wrtout(iout,message,'COLL')
! end do
! ENDDEBUG

 end if

end subroutine prtxf
!!***
