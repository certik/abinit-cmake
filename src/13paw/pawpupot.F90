!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawpupot
!! NAME
!! pawpupot
!!
!! FUNCTION
!! Compute LDA+U potential
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (BA,FJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/Infos/copyright
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors.
!!
!! INPUTS
!!  ispden=current spin component
!!  nspden=number of spin components
!!  pawprtvol=control print volume and debugging output for PAW
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!    %noccmmp(2*lpawu+1,2*lpawu+1,nspden)=density matrix in the augm. region
!!    %nocctot(nspden)=number of electrons in the correlated subspace
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data:
!!     %lpawu=l used for lda+u
!!     %vee(2*lpawu+1*4)=screened coulomb matrix
!!
!! OUTPUT
!!  vpawu(lpawu*2+1,lpawu*2+1)=lda+u potential (see eg PRB 52, 5467 (1995))
!!
!! SIDE EFFECTS
!!  VUKS=Sum over spins and atoms of nocc_mmp*Pot_(LDA+U)
!!
!! PARENTS
!!      pawdij
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 subroutine pawpupot(ispden,nspden,paw_ij,pawprtvol,pawtab,vpawu,VUKS)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ispden,nspden,pawprtvol
 real(dp),intent(inout) :: VUKS
 type(paw_ij_type),intent(in) :: paw_ij
 type(pawtab_type),intent(in) :: pawtab
!arrays
 real(dp),intent(inout) :: vpawu(pawtab%lpawu*2+1,pawtab%lpawu*2+1)

!Local variables ---------------------------------------
!scalars
!Option for interaction energy in case of non-collinear magnetism:
!           1: E_int=-U/4.N.(N-2)
!           2: E_int=-U/2.(Nup.(Nup-1)+Ndn.(Ndn-1))
 integer,parameter :: option_interaction=1

 integer :: jspden,lpawu,m1,m11,m2,m21,m3,m31,m4,m41
 real(dp) :: mnorm,mx,my,mz,n_sig,n_msig,n_tot,n34_msig,n34_sig,VUKStemp
 character(len=500) :: message
!array
 real(dp),parameter :: factcg(3:4)=(/one,-one/)

! *****************************************************

!DEBUG
!write(6,*)' pawpupot : enter '
!ENDDEBUG

 lpawu=pawtab%lpawu
 vpawu=zero

!======================================================
!Compute LDA+U Potential on the basis of projectors.
!cf PRB 52 5467 (1995)
!-----------------------------------------------------

 if (ispden<=2) then
  jspden=min(nspden,2)-ispden+1
  do m1=-lpawu,lpawu
   m11=m1+lpawu+1
   do m2=-lpawu,lpawu
    m21=m2+lpawu+1
    do m3=-lpawu,lpawu
     m31=m3+lpawu+1
     do m4=-lpawu,lpawu
      m41=m4+lpawu+1
      n34_sig =paw_ij%noccmmp(m31,m41,ispden)
      n34_msig=paw_ij%noccmmp(m31,m41,jspden)
      vpawu(m11,m21)=vpawu(m11,m21) &
&      +n34_msig*pawtab%vee(m11,m31,m21,m41) &
&      +n34_sig*(pawtab%vee(m11,m31,m21,m41)-pawtab%vee(m11,m31,m41,m21))
     end do
    end do
    if(abs(pawprtvol)>=3) then
     if(m11/=m21) then
      write(message,'(a,i4,i4,e20.10)') "vu=",m11,m21,vpawu(m11,m21)
      call wrtout(06,message,'COLL')
      write(message,'(a,e20.10)') "vupred=",-pawtab%upawu*paw_ij%noccmmp(m21,m11,ispden)
      call wrtout(06,message,'COLL')
     end if
    end if
   end do ! m2

   if (nspden<=2) then
    n_sig =paw_ij%nocctot(ispden)
    n_msig=paw_ij%nocctot(jspden)
    n_tot =n_sig+n_msig
   else
    n_tot= paw_ij%nocctot(1)+paw_ij%nocctot(2)
    mz   = paw_ij%nocctot(1)-paw_ij%nocctot(2)
    mx   = two*paw_ij%nocctot(3)
    my   =-two*paw_ij%nocctot(4)
    mnorm=sqrt(mx*mx+my*my+mz*mz)
    if (ispden==1) then
     n_sig =half*(n_tot+mnorm)
     n_msig=half*(n_tot-mnorm)
    else
     n_sig =half*(n_tot-mnorm)
     n_msig=half*(n_tot+mnorm)
    end if
   end if

!  Full localized limit
   if(pawtab%usepawu==1) then
    vpawu(m11,m11)=vpawu(m11,m11)-pawtab%upawu*(n_tot-half)
    if (nspden/=4.or.option_interaction==2) then
     vpawu(m11,m11)=vpawu(m11,m11)+pawtab%jpawu*(n_sig-half)
    else
     vpawu(m11,m11)=vpawu(m11,m11)+half*pawtab%jpawu*(n_tot-one)
    end if

!   Around mean field
   else if(pawtab%usepawu==2) then
    vpawu(m11,m11)=vpawu(m11,m11)-n_msig*pawtab%upawu &
&    -n_sig*(pawtab%upawu-pawtab%jpawu) &
&    *(dble(2*lpawu)/dble(2*lpawu+1))
   end if

   if (abs(pawprtvol)>=3) then
    write(message,'(a,i4,i4,2x,e20.10)') "vudiag= ",m11,m11,vpawu(m11,m11)
    call wrtout(06,  message,'COLL')
    write(message,'(a,e20.10)') "vudiagpred= ",pawtab%upawu*(half-paw_ij%noccmmp(m11,m11,ispden))
    call wrtout(06,  message,'COLL')
   end if
  end do ! m1

 end if ! ispden<=2

!Non-collinear magnetism: add non-diagonal term; see (Eq 6) in PRB 72, 024458 (2005)
 if (ispden>=3) then
  do m1=-lpawu,lpawu
   m11=m1+lpawu+1
   do m2=-lpawu,lpawu
    m21=m2+lpawu+1
    do m3=-lpawu,lpawu
     m31=m3+lpawu+1
     do m4=-lpawu,lpawu
      m41=m4+lpawu+1
      vpawu(m11,m21)=vpawu(m11,m21)-factcg(ispden)*paw_ij%noccmmp(m41,m31,ispden)*pawtab%vee(m11,m31,m41,m21)
     end do
    end do
   end do
  end do
 end if

!Printing for test
 if (abs(pawprtvol)>=3) then
  if (ispden==1) VUKS=zero
  VUKStemp=zero
  do m1=-lpawu,lpawu
   m11=m1+lpawu+1
   do m2=-lpawu,lpawu
    m21=m2+lpawu+1
    VUKStemp=VUKStemp+vpawu(m11,m21)*paw_ij%noccmmp(m11,m21,ispden)
    write(message,'(a,e20.10,e20.10)') "vuks,m1,m2,vpawu,noccmmp= ", &
&    vpawu(m11,m21),paw_ij%noccmmp(m11,m21,ispden)
    call wrtout(06,  message,'COLL')
   end do
  end do
  VUKS=VUKS+VUKStemp
  write(message,*) "pawpupot: VUKStemp= ",ispden,VUKStemp
  call wrtout(06,  message,'COLL')
  if (ispden==nspden) then
   write(message,*) "pawpupot: VUKS= ",ispden,VUKS
   call wrtout(06,  message,'COLL')
  end if
 end if

!DEBUG
!write(6,*)' pawpupot : exit '
!ENDDEBUG


 end subroutine pawpupot
!!***
