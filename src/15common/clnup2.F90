!{\src2tex{textfont=tt}}
!!****f* ABINIT/clnup2
!! NAME
!! clnup2
!!
!!
!! FUNCTION
!! Perform more "cleanup" after completion of iterations.
!! This subroutine prints out more breakdown of force
!! information, shifts of atomic positions, and stresses.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  fred(3,natom)=d(E_total)/d(xred) derivatives (hartree)
!!  grewtn(3,natom)=d(E_Ewald)/d(xred) derivatives (hartree)
!!  grxc(3,natom)=d(Exc)/d(xred) derivatives (0 without core charges)
!!  iscf=parameter controlling scf or non-scf iterations
!!  natom=number of atoms in unit cell
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  prtfor= >0 if forces have to be printed (0 otherwise)
!!  prtstr= >0 if stresses have to be printed (0 otherwise)
!!  prtvol=control print volume and debugging output
!!  start(3,natom)=starting coordinates in terms of real space
!!   primitive translations
!!  strten(6)=components of the stress tensor (hartree/bohr^3)
!!  synlgr(3,natom)=d(E_nlpsp)/d(xred) derivatives (hartree)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xred(3,natom)=final coordinates in terms of primitive translations
!!
!! OUTPUT
!!  (only print)
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine clnup2(n1xccc,fred,gresid,grewtn,grxc,iscf,natom,prtfor,prtstr,prtvol,&
&                  start,strten,synlgr,usepaw,xred)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iscf,n1xccc,natom,prtfor,prtstr,prtvol,usepaw
!arrays
 real(dp),intent(in) :: fred(3,natom),gresid(3,natom),grewtn(3,natom)
 real(dp),intent(in) :: grxc(3,natom),start(3,natom),strten(6),synlgr(3,natom)
 real(dp),intent(in) :: xred(3,natom)

!Local variables-------------------------------
 character(len=*), parameter :: format01020 ="(i5,1x,3f20.12)"
!scalars
 integer :: iatom,mu
 real(dp) :: devsqr
 character(len=500) :: message

! *************************************************************************
!
!DEBUG
!write(6,*)' clnup2 : enter '
!ENDDEBUG

!Only print additional info for scf calculations
 if (iscf>0) then

  if((prtvol>=10).and.(prtfor>0))then

   write(message, '(a,10x,a)' ) ch10,&
&   '===> extra information on forces <==='
   call wrtout(ab_out,message,'COLL')

   write(message, '(a)' ) ' ewald contribution to reduced grads'
   call wrtout(ab_out,message,'COLL')
   do iatom=1,natom
    write(message,format01020) iatom,(grewtn(mu,iatom),mu=1,3)
    call wrtout(ab_out,message,'COLL')
   end do

   write(message, '(a)' ) ' nonlocal contribution to red. grads'
   call wrtout(ab_out,message,'COLL')
   do iatom=1,natom
    write(message,format01020) iatom,(synlgr(mu,iatom),mu=1,3)
    call wrtout(ab_out,message,'COLL')
   end do

   write(message, '(a)' ) ' local psp contribution to red. grads'
   call wrtout(ab_out,message,'COLL')
   if (n1xccc/=0) then
    do iatom=1,natom
     write(message,format01020) iatom,fred(:,iatom)-&
&     (grewtn(:,iatom)+synlgr(:,iatom)+grxc(:,iatom)+gresid(:,iatom))
     call wrtout(ab_out,message,'COLL')
    end do
   else
    do iatom=1,natom
     write(message,format01020) iatom,fred(:,iatom)-&
&     (grewtn(:,iatom)+synlgr(:,iatom)+gresid(:,iatom))
     call wrtout(ab_out,message,'COLL')
    end do
   end if

   if (n1xccc/=0) then
    write(message, '(a)' ) ' core charge xc contribution to reduced grads'
    call wrtout(ab_out,message,'COLL')
    do iatom=1,natom
     write(message,format01020) iatom,(grxc(mu,iatom),mu=1,3)
     call wrtout(ab_out,message,'COLL')
    end do
   end if

   write(message, '(a)' ) ' residual contribution to red. grads'
   call wrtout(ab_out,message,'COLL')
   do iatom=1,natom
    write(message,format01020) iatom,(gresid(mu,iatom),mu=1,3)
    call wrtout(ab_out,message,'COLL')
   end do

  end if

! Compute mean squared deviation from starting coords
  devsqr=0.0_dp
  do iatom=1,natom
   do mu=1,3
    devsqr=devsqr+(xred(mu,iatom)-start(mu,iatom))**2
   end do
  end do

! When shift is nonnegligible then print values
  if (devsqr>1.d-14) then
   write(message, '(a,1p,e12.4,3x,a)' ) &
&   ' rms coord change=',sqrt(devsqr/dble(3*natom)),&
&   'atom, delta coord (reduced):'
   call wrtout(ab_out,message,'COLL')
   do iatom=1,natom
    write(message, '(1x,i5,2x,3f20.12)' ) iatom,&
&    (xred(mu,iatom)-start(mu,iatom),mu=1,3)
    call wrtout(ab_out,message,'COLL')
   end do
  end if

! Write out stress results
  if (prtstr>0) then
   write(message, '(a,a)' ) ch10,&
&   ' Cartesian components of stress tensor (hartree/bohr^3)'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')

   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(1 1)=',strten(1),'  sigma(3 2)=',strten(4)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(2 2)=',strten(2),'  sigma(3 1)=',strten(5)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(3 3)=',strten(3),'  sigma(2 1)=',strten(6)
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')

!  Also output the pressure (minus one third the trace of the stress
!  tensor.
   write(message, '(a,a,es12.4,a)' ) ch10,&
&   '-Cartesian components of stress tensor (GPa)         [Pressure=',&
&   -(strten(1)+strten(2)+strten(3))*HaBohr3_GPa/3.0_dp,' GPa]'

   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')

   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '- sigma(1 1)=',strten(1)*HaBohr3_GPa,&
&   '  sigma(3 2)=',strten(4)*HaBohr3_GPa
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '- sigma(2 2)=',strten(2)*HaBohr3_GPa,&
&   '  sigma(3 1)=',strten(5)*HaBohr3_GPa
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '- sigma(3 3)=',strten(3)*HaBohr3_GPa,&
&   '  sigma(2 1)=',strten(6)*HaBohr3_GPa
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
  end if

! Last end if above refers to iscf > 0
 end if

!DEBUG
!write(6,*)' clnup2 : exit '
!ENDDEBUG

end subroutine clnup2
!!***
