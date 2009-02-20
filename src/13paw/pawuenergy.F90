!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawuenergy
!! NAME
!! pawuenergy
!!
!! FUNCTION
!! Compute contributions to energy for PAW+U calculations
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (BA,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  iatom=index of current atom
!!  nspden: Number of spin component
!!  pawprtvol=control print volume and debugging output for PAW
!!  pawtab <type(pawtab_type)>=paw tabulated starting data:
!!     %lpawu=l used for lda+u
!!     %vee(2*lpawu+1*4)=screened coulomb matrix
!!  paw_ij <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!     %noccmmp(2*lpawu+1,2*lpawu+1,nspden)=density matrix in the PAW augm. region
!!     %nocctot(nspden)=number of electrons in the correlated subspace
!!
!! OUTPUT
!!  eldaumdc= PAW+U contribution to total energy
!!  eldaumdcdc= PAW+U contribution to double-counting total energy
!!
!! PARENTS
!!      pawdenpot
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pawuenergy(iatom,eldaumdc,eldaumdcdc,pawprtvol,pawtab,paw_ij,nspden)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: iatom,nspden,pawprtvol
 real(dp),intent(inout) :: eldaumdc,eldaumdcdc
 type(paw_ij_type),intent(in) :: paw_ij
 type(pawtab_type),intent(in) :: pawtab

!Local variables ---------------------------------------
!scalars
!Option for interaction energy in case of non-collinear magnetism:
!           1: E_int=-U/4.N.(N-2)
!           2: E_int=-U/2.(Nup.(Nup-1)+Ndn.(Ndn-1))
 integer,parameter :: option_interaction=1

 integer :: ispden,jspden,lpawu,m1,m11,m2,m21,m3,m31,m4,m41
 real(dp) :: edcdctemp,edctemp,eldautemp,mnorm,mx,my,mz,ndn,ntot,nup
 real(dp) :: n12_im,n12_re,n12_sig,n21_im,n21_re,n34_im,n34_msig,n34_re,n34_sig,n43_im,n43_re
 character(len=500) :: message

! *****************************************************

!DEBUG
!write(6,*) ' pawuenergy : enter '
!ENDDEBUG

 lpawu=pawtab%lpawu

!======================================================
!Compute LDA+U Energy
!-----------------------------------------------------

 eldautemp=zero
 do ispden=1,min(nspden,2)
  jspden=min(nspden,2)-ispden+1
  do m1=-lpawu,lpawu
   m11=m1+lpawu+1
   do m2=-lpawu,lpawu
    m21=m2+lpawu+1
    n12_sig=paw_ij%noccmmp(m11,m21,ispden)
    do m3=-lpawu,lpawu
     m31=m3+lpawu+1
     do m4=-lpawu,lpawu
      m41=m4+lpawu+1
      n34_sig =paw_ij%noccmmp(m31,m41,ispden)
      n34_msig=paw_ij%noccmmp(m31,m41,jspden)
      eldautemp=eldautemp &
&      + n12_sig*n34_msig*pawtab%vee(m11,m31,m21,m41) &
&      + n12_sig*n34_sig *(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
     end do ! m4
    end do ! m3
   end do ! m2
  end do ! m1
 end do ! ispden
 if (nspden==1) eldautemp=two*eldautemp ! Non-magn. system: sum up and dn energies

!Non-collinear magnetism: add non-diagonal term; see (Eq 3) in PRB 72, 024458 (2005)
 if (nspden==4) then
  do m1=-lpawu,lpawu
   m11=m1+lpawu+1
   do m2=-lpawu,lpawu
    m21=m2+lpawu+1
    n12_re=paw_ij%noccmmp(m11,m21,3)
    n12_im=paw_ij%noccmmp(m11,m21,4)
    n21_re=paw_ij%noccmmp(m21,m11,3)
    n21_im=paw_ij%noccmmp(m21,m11,4)
    do m3=-lpawu,lpawu
     m31=m3+lpawu+1
     do m4=-lpawu,lpawu
      m41=m4+lpawu+1
      n34_re=paw_ij%noccmmp(m31,m41,3)
      n34_im=paw_ij%noccmmp(m31,m41,4)
      n43_re=paw_ij%noccmmp(m41,m31,3)
      n43_im=paw_ij%noccmmp(m41,m31,4)
      eldautemp=eldautemp-pawtab%vee(m11,m31,m41,m21) &
&      *(n12_re*n43_re+n12_im*n43_im &
&      +n21_re*n34_re+n21_im*n34_im)
     end do ! m4
    end do ! m3
   end do ! m2
  end do ! m1
 end if

!Divide eldautemp by 2; see (Eq 1) in PRB 77, 155104 (2008)
 eldautemp=half*eldautemp

 if (nspden==1) then
  ntot=two*paw_ij%nocctot(1)
  nup=paw_ij%nocctot(1)
  ndn=paw_ij%nocctot(1)
 else if (nspden==2) then
  ntot=paw_ij%nocctot(1)+paw_ij%nocctot(2)
  nup=paw_ij%nocctot(1)
  ndn=paw_ij%nocctot(2)
 else if (nspden==4) then
  ntot= paw_ij%nocctot(1)+paw_ij%nocctot(2)
  mz  = paw_ij%nocctot(1)-paw_ij%nocctot(2)
  mx  = two*paw_ij%nocctot(3)
  my  =-two*paw_ij%nocctot(4)
  mnorm=sqrt(mx*mx+my*my+mz*mz)
  nup=half*(ntot+mnorm)
  ndn=half*(ntot-mnorm)
 end if

 edcdctemp=zero;edctemp=zero

!Full localized limit
 if(pawtab%usepawu==1) then
  edcdctemp=edcdctemp-half*pawtab%upawu*ntot**2
  edctemp  =edctemp  +half*pawtab%upawu*(ntot*(ntot-one))
  if (nspden/=4.or.option_interaction==2) then
   edcdctemp=edcdctemp+half*pawtab%jpawu*(nup**2+ndn**2)
   edctemp  =edctemp  -half*pawtab%jpawu*(nup*(nup-one)+ndn*(ndn-one))
  else
   edcdctemp=edcdctemp+quarter*pawtab%jpawu*ntot**2
   edctemp  =edctemp  -quarter*pawtab%jpawu*(ntot*(ntot-two))
  end if

! Around mean field
 else if(pawtab%usepawu==2) then
  edctemp=edctemp+pawtab%upawu*(nup*ndn)&
&  +half*(pawtab%upawu-pawtab%jpawu)*(nup**2+ndn**2) &
&  *(dble(2*lpawu)/dble(2*lpawu+1))
  edcdctemp=-edctemp
 end if

 eldaumdc  =eldaumdc  +eldautemp-edctemp
 eldaumdcdc=eldaumdcdc-eldautemp-edcdctemp

 write(message, '(5a,i4)')ch10,'======= LDA+U Energy terms (in Hartree) ====',ch10,&
& ch10,' For Atom ',iatom
 call wrtout(6,message,'COLL')
 write(message, '(a)' )"   Contributions to the direct expression of energy:"
 call wrtout(06,  message,'COLL')
 write(message,fmt=11) "     Double counting  correction   =",edctemp
 call wrtout(06,  message,'COLL')
 write(message,fmt=11) "     Interaction energy            =",eldautemp
 call wrtout(06,  message,'COLL')
 write(message,fmt=11) "     Total LDA+U Contribution      =",eldautemp-edctemp
 call wrtout(06,  message,'COLL')
 write(message, '(a)' )' '
 call wrtout(06,  message,'COLL')
 write(message, '(a)' )"   For the ""Double-counting"" decomposition:"
 call wrtout(06,  message,'COLL')
 write(message,fmt=11) "     LDA+U Contribution            =",-eldautemp-edcdctemp
 call wrtout(06,  message,'COLL')
 11 format(a,e20.10)
 if(abs(pawprtvol)>=2) then
  write(message,fmt=11)"     edcdctemp                     =",edcdctemp
  call wrtout(06,  message,'COLL')
  write(message,fmt=11)"     eldaumdcdc for current atom   =",-eldautemp-edcdctemp
  call wrtout(06,  message,'COLL')
  write(message, '(a)' )' '
  call wrtout(06,  message,'COLL')
  write(message,fmt=11)"   pawuenergy: -VUKS pred          =",eldaumdcdc-eldaumdc
  call wrtout(06,  message,'COLL')
 end if
 write(message, '(a)' )' '
 call wrtout(06,  message,'COLL')

!DEBUG
!write(6,*) ' pawuenergy : enter '
!ENDDEBUG

 end subroutine pawuenergy
!!***
