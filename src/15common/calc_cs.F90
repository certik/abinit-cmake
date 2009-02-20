!{\src2tex{textfont=tt}}
!!****f* ABINIT/calc_cs
!! NAME
!! calc_cs
!!
!! FUNCTION
!! calculation and output of chemical shielding tensor at each atomic site
!!
!! COPYRIGHT
!! Copyright (C) 2005-2008 ABINIT group (SM,VR,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  corecs(ntypat)=core shieldings in ppm of different atomic nuclei
!!  natom=number of atoms in cell.
!!  nspden=number of spin densities
!!  ntypat=number of atom types
!!  occopt=option describing occupation of bands
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  prtcs=1 to print summary output, 2 for detailed output
!!  usepaw=1 if we are using PAW formalism, 0 else
!!
!! OUTPUT
!!  (only writing, printing)
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

 subroutine calc_cs(corecs,natom,nspden,ntypat,occopt,pawang,pawrad,pawrhoij,pawtab,prtcs,typat,usepaw)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13paw
 use interfaces_linalg
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nspden,ntypat,occopt,prtcs,usepaw
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: corecs(ntypat)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: INFO,LDA,LWORK,N,iatom,ii
 real(dp) :: csxx,csyy,cszz
 character(len=500) :: message
!arrays
 real(dp) :: eigval(3),matr(3,3),work(8)
 real(dp),allocatable :: cs(:,:,:),cs_core(:,:,:),cs_dia(:,:,:)

! ************************************************************************

!DEBUG
!write(*,*)' calc_cs : enter'
!ENDDEBUG
!Compatibility tests
 if (usepaw /= 1) then
  write (message,'(4a)')' calc_cs : ERROR- ',ch10,&
&  ' usepaw /= 1 but CS calculation requires PAW ',ch10
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if
 if (occopt > 1) then
  write (message,'(4a)')' calc_cs : ERROR- ',ch10,&
&  ' occopt > 1 but CS calculation requires occopt = 1 or 0 (no metals for now, sorry !)  ',ch10
  call wrtout(06,message,'COLL')
  call leave_new('COLL')
 end if

 allocate(cs(3,3,natom),cs_core(3,3,natom),cs_dia(3,3,natom))
 cs_core(:,:,:) = zero 
 cs_dia(:,:,:) = zero 

 do iatom = 1, natom
  do ii = 1, 3
   cs_core(ii,ii,iatom) = corecs(typat(iatom)) 
  end do
 end do

 call make_cs_dia(cs_dia,natom,ntypat,pawang,pawrhoij,pawrad,pawtab,typat)

 cs(:,:,:) = cs_core(:,:,:) + cs_dia(:,:,:)

 write(message,'(a,a,a)' ) ch10,' Chemical Shielding Calculation ',ch10
 call wrtout(ab_out,message,'COLL')

 LDA=3; LWORK=8;N=3 ! these parameters are needed for the LAPACK dsyev routine 
 do iatom = 1, natom
  matr(:,:) = cs(:,:,iatom)
  call dsyev('V','U',N,matr,LDA,eigval,work,LWORK,INFO) ! get eigenvalues and eigenvectors
  if (eigval(3) > abs(eigval(1)) ) then ! In NMR, the convention is that whatever component is
!  largest in magnitude is called zz, next comes xx, then yy
   cszz = eigval(3)
   csxx = eigval(1)
   csyy = eigval(2)
  else
   cszz = eigval(1)
   csxx = eigval(3)
   csyy = eigval(2)
  end if 
! we always write the principle components, these are the NMR observables
  write(message,'(a,i3,a,i3,a,f7.3,a,f7.3,a,f7.3)')&
&  ' Atom ',iatom,', typat ',typat(iatom),': SigZZ = ',cszz,' SigXX = ',csxx,' SigYY = ',csyy
  call wrtout(ab_out,message,'COLL')
  write(message,'(3a)') ch10,ch10,ch10
  call wrtout(ab_out,message,'COLL')
  if (prtcs > 1) then ! print detailed results on chemical shielding
   write(message,'(a,a,f7.3,a,a,3f7.3)')ch10,'      cs eigval : ',eigval(1),ch10,&
&   '-        eigvec : ',matr(1,1),matr(2,1),matr(3,1)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,f7.3,a,a,3f7.3)')'      cs eigval : ',eigval(2),ch10,&
&   '-        eigvec : ',matr(1,2),matr(2,2),matr(3,2)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,f7.3,a,a,3f7.3)')'      cs eigval : ',eigval(3),ch10,&
&   '-        eigvec : ',matr(1,3),matr(2,3),matr(3,3)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,a,3f7.3)')ch10,'      total cs : ',cs(1,1,iatom),cs(1,2,iatom),cs(1,3,iatom)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3f7.3)')'      total cs : ',cs(2,1,iatom),cs(2,2,iatom),cs(2,3,iatom)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3f7.3,a)')'      total cs : ',cs(3,1,iatom),cs(3,2,iatom),cs(3,3,iatom),ch10
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,a,3f7.3)')ch10,'      cs_dia : ',cs_dia(1,1,iatom),cs_dia(1,2,iatom),cs_dia(1,3,iatom)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3f7.3)')'      cs_dia : ',cs_dia(2,1,iatom),cs_dia(2,2,iatom),cs_dia(2,3,iatom)
   call wrtout(ab_out,message,'COLL')
   write(message,'(a,3f7.3,a)')'      cs_dia : ',cs_dia(3,1,iatom),cs_dia(3,2,iatom),cs_dia(3,3,iatom),ch10
   call wrtout(ab_out,message,'COLL')
  end if
 end do
 deallocate(cs,cs_core,cs_dia)


!DEBUG
!write(6,*)' calc_cs : exit '
!stop
!ENDDEBUG

 end subroutine calc_cs
!!***
