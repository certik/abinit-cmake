!{\src2tex{textfont=tt}}
!!****f* ABINIT/outeps
!! NAME
!! outeps
!!
!! FUNCTION
!!  Write the independent matrix elements of epsilon^{-1}_q(G1,G2)(omega) on file
!!
!! COPYRIGHT
!!  Copyright (C) 2007-2008 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  Gsphere<Gvectors_type>
!!   %ng=number of G vectors in the sphere
!!   %gvec(3,ng)= G vectors in reduced coordinates
!!   %gmet(3,3)=metric tensor in reciprocal space
!!  Gpairs_q<Gpairs_type>
!!   %ngpi=nuber of independent (G1,G1) pairs   
!!   %ip2fp(2,ngpi)= index of G1 and G2 in the array gvec for each ngpi independent pair
!!  nomega=number of frequencies
!!  title=title describing what is printed
!!  nge=Number of G vectors in epsilon
!!  eps(nge,nge,nomega)=matrix to be printed
!!  omega(nomega)=frequencies
!!  prtvol= if different from 0 first and last 150 independent pairs are written 
!!          if ==0 only first and last 50 pairs are written
!!  unt=unit number of external file
!!
!! OUTPUT
!!  Only writing 
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine outeps(nge,nomega,omega,eps,Gsphere,Gpairs_q,title,unt,prtvol)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_12geometry
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nge,nomega,prtvol,unt
 character(len=500),intent(in) :: title
 type(Gpairs_type),intent(in) :: Gpairs_q
 type(Gvectors_type),intent(in) :: Gsphere
!arrays
 complex,intent(in) :: eps(nge,nge,nomega),omega(nomega)

!Local variables-------------------------------
!scalars
 integer :: enough,f1,f2,ii,io,ip,ng,ngpi,npmax
 real(dp) :: e1,e2
 character(len=500) :: msg
!arrays
 integer,pointer :: gvec(:,:)
 real(dp) :: gmet(3,3)

!************************************************************************

 ng=Gsphere%ng ; gmet=Gsphere%gmet
 gvec => Gsphere%gvec(1:3,1:ng)
 ngpi=Gpairs_q%niggp

 npmax=50 ; if (ngpi<50) npmax=ngpi
 if (prtvol/=0) then 
  npmax=150 ; if (ngpi<150) npmax=ngpi
 end if

 if (Gpairs_q%niggp<Gpairs_q%ng**2) then 
  !
  ! === For this q-point we used al least 2 symmetries, tables in Gpairs are associated === 
  write(msg,'(3a)')ch10,&
&  '# Independent matrix elements of: ',TRIM(title)
  call wrtout(unt,msg,'COLL')

  do io=1,nomega
   ! === Write header ===
   write(msg,'(2a,i3,a,2f9.4,3a)')ch10,&
 &  '#',io,'-th omega :',omega(io)*Ha_eV,' [eV] ',ch10,&
 &  '#         G1         E1 [H]           G2         E2 [H]         Re        Im '
   call wrtout(unt,msg,'COLL')

   enough=0
   do ip=1,ngpi
    enough=enough+1
    !if (enough>npmax.and.enough<ngpi-npmax.and.prtvol==0) CYCLE
    if (enough>npmax.and.enough<ngpi-npmax) CYCLE
    f1=Gpairs_q%ip2fp(1,ip)
    f2=Gpairs_q%ip2fp(2,ip)
    e1=normv(DBLE(gvec(:,f1)),gmet,'g') ; e1=half*e1**2
    e2=normv(DBLE(gvec(:,f2)),gmet,'g') ; e2=half*e2**2

    write(unt,'(2(3i6,4x,f4.2,2x),4x,2(f8.4,2x))')gvec(:,f1),e1,gvec(:,f2),e2,eps(f1,f2,io)
    if (enough==npmax.and.prtvol==0) then 
     write(unt,*)('-',ii=1,77)
     write(unt,*)' outeps: prtvol=0, stop printing more G1 G2 information' 
     write(unt,*)('-',ii=1,77)
    end if 
   end do !ip
  end do !io

 else 
  !
  ! === We do not have irreducible pairs ===
  write(msg,'(2a,i4,2a)')ch10,&
&  '# First and last ',npmax,' matrix elements of ',TRIM(title)
  call wrtout(unt,msg,'COLL')

  do io=1,nomega
   ! === Write header ===
   write(msg,'(2a,i3,a,2f9.4,3a)')ch10,&
 &  '#',io,'-th omega :',omega(io)*Ha_eV,' [eV] ',ch10,&
 &  '#         G1         E1 [H]           G2         E2 [H]         Re        Im '
   call wrtout(unt,msg,'COLL')

   enough=0
   do f1=1,ng
    e1=normv(DBLE(gvec(:,f1)),gmet,'g') ; e1=half*e1**2
    do f2=1,ng
     enough=enough+1
     !if (enough>npmax.and.enough<ngpi-npmax.and.prtvol==0) CYCLE
     if (enough>npmax.and.enough<ngpi-npmax) CYCLE
     e2=normv(DBLE(gvec(:,f2)),gmet,'g') ; e2=half*e2**2
     write(unt,'(2(3i6,4x,f4.2,2x),4x,2(f8.4,2x))')gvec(:,f1),e1,gvec(:,f2),e2,eps(f1,f2,io)
     if (enough==npmax.and.prtvol==0) then 
      write(unt,*)('-',ii=1,77)
      write(unt,*)' outeps : prtvol=0, stop printing more G-G'' information' 
      write(unt,*)('-',ii=1,77)
     end if 
    end do !f2
   end do !f1

  end do !io
 end if

end subroutine outeps
!!***
