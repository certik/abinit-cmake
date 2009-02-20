!{\src2tex{textfont=tt}}
!!****f* ABINIT/initrhoij
!! NAME
!! initrhoij
!!
!! FUNCTION
!! Initialize PAW rhoij occupancies (in packed storage)
!! from atomic ones
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=lmn
!!  lexexch(ntypat)=l on which local exact-exchange is applied for a given type of atom
!!  lmnmax=max number of (l,m,n) comp. over all type of psps
!!  lpawu(ntypat)=l on which U is applied for a given type of atom (PAW+U)
!!  natom=number of atoms
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of atom types
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!                                     (containing initial rhoij)
!!  spinat(3,natom)=initial spin of each atom, in unit of hbar/2.
!!  typat(natom)=type of each atom
!!
!! OUTPUT
!!  pawrhoij(natom) <type(pawrhoij_type)>=rhoij quantities for each atom
!!                                        in packed storage
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!      leave_new,rhoij_alloc,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine initrhoij(indlmn,lexexch,lmnmax,lpawu,natom,nspden,nsppol,ntypat,pawrhoij,pawtab,spinat,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lmnmax,natom,nspden,nsppol,ntypat
 character(len=500) :: message
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),lexexch(ntypat),lpawu(ntypat)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: spinat(3,natom)
 type(pawrhoij_type),intent(out) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables ---------------------------------------
!Arrays
!scalars
 integer :: iatom,ilmn,ispden,itypat,j0lmn,jl,jlmn,jspden,klmn,nselect
 real(dp) :: ratio,ro,roshift,zratio,zz
 logical :: test_exexch,test_pawu
!arrays
 integer,allocatable :: nlmn(:)

!************************************************************************

!PAW+U and local exact-exchange restriction
 do itypat=1,ntypat
  if (lpawu(itypat)/=lexexch(itypat).and.&
&  lpawu(itypat)/=-1.and.lexexch(itypat)/=-1) then
   write(message, '(4a)' ) ch10,' initrhoij: ERROR - ',&
&   ch10,'  lpawu must be equal to lexexch !'
   call wrtout(ab_out,message,'COLL')
   call wrtout(06,  message,'COLL')
   call leave_new('COLL')
  end if
 end do

 ratio=one;if (nspden==2) ratio=half

 allocate(nlmn(ntypat))
 do itypat=1,ntypat
  nlmn(itypat)=pawtab(itypat)%lmn_size
 end do
 call rhoij_alloc(1,nlmn,nspden,nsppol,pawrhoij,typat)
 deallocate(nlmn)

 do iatom=1,natom
  itypat=typat(iatom)

! Determine Z (trace of rhoij0 or part of it)
  zz=zero
  do jlmn=1,pawtab(itypat)%lmn_size
   jl=indlmn(1,jlmn,itypat)
   j0lmn=jlmn*(jlmn-1)/2
   test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
   test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
   do ilmn=1,jlmn
    klmn=j0lmn+ilmn
    if ((ilmn==jlmn).and.test_pawu.and.test_exexch) &
&    zz=zz+pawtab(itypat)%rhoij0(klmn)
   end do
  end do

! Compute rhoij from tabulated value and magnetization
  do ispden=1,nspden

   zratio=zero
   roshift=one
   ratio=one
   if (nspden==2) then
    ratio=half
    if ((spinat(3,iatom)>zero.and.ispden==1).or.&
&    (spinat(3,iatom)<zero.and.ispden==2)) then
     zratio=two*abs(spinat(3,iatom))/zz
    end if
   else if (nspden==4.and.ispden>=2) then
    roshift=zero
    zratio=spinat(ispden-1,iatom)/zz
   end if

   nselect=0
   do jlmn=1,pawtab(itypat)%lmn_size
    jl=indlmn(1,jlmn,itypat)
    j0lmn=jlmn*(jlmn-1)/2
    test_pawu=(lpawu(itypat)==-1.or.lpawu(itypat)==jl)
    test_exexch=(lexexch(itypat)==-1.or.lexexch(itypat)==jl)
    do ilmn=1,jlmn
     klmn=j0lmn+ilmn
     ro=pawtab(itypat)%rhoij0(klmn)
     if ((ilmn==jlmn).and.test_pawu.and.test_exexch) then
      ro=ro*ratio*(roshift+zratio)
     else
      ro=ro*ratio*roshift
     end if

     if (abs(ro)>tol10) then
      pawrhoij(iatom)%rhoijp(klmn,ispden)=ro
     else
      pawrhoij(iatom)%rhoijp(klmn,ispden)=zero
     end if

     if (ispden==nspden) then
      if (any(abs(pawrhoij(iatom)%rhoijp(klmn,:))>tol10)) then
       nselect=nselect+1
       pawrhoij(iatom)%rhoijselect(nselect)=klmn
       do jspden=1,nspden
        pawrhoij(iatom)%rhoijp(nselect,jspden)=pawrhoij(iatom)%rhoijp(klmn,jspden)
       end do
      end if
     end if

    end do
   end do

  end do
  pawrhoij(iatom)%nrhoijsel=nselect

 end do ! iatom

end subroutine initrhoij
!!***
