!{\src2tex{textfont=tt}}
!!****f* ABINIT/pawprt
!! NAME
!! pawprt
!!
!! FUNCTION
!! Print out data concerning PAW formalism
!! (pseudopotential strength, augmentation occupancies...)
!! To be called at the end of the SCF cycle
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  indlmn(6,i,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn
!!  enunit=parameter determining units of output energies
!!  lmnmax=if useylm=1, max number of (l,m,n) comp. over all types of psps
!!  natom=number of atoms in cell
!!  ntypat = number of atom types
!!  paw_ij(natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawprtvol= printing volume
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  typat(natom)=type of each atom
!!
!! OUTPUT
!!  (only printing)
!!
!! PARENTS
!!      outscfcv,sigma
!!
!! CHILDREN
!!      print_ij,setnoccmmp,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine pawprt(indlmn,enunit,lmnmax,natom,ntypat,paw_ij,pawprtvol,pawrhoij,pawtab,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
 use interfaces_13paw, except_this_one => pawprt
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: enunit,lmnmax,natom,ntypat,pawprtvol
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),typat(natom)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)
 type(pawrhoij_type),intent(in) :: pawrhoij(natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer,parameter :: natmax=2
 integer :: iat,iatom,iatom1,im1,im2,ispden,itypat,klmn,ll,llp,natprt,nspden,nsppol,unt
 real(dp) :: mnorm,mx,my,mz,ntot,ro,valmx
 logical :: antiferro,useexexch,usepawu
 type(pawang_type):: pawang_dum
 character(len=7),parameter :: dspin1(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
 character(len=8),parameter :: dspin2(6)=(/"up      ","down    ","dens (n)","magn (x)","magn (y)","magn (z)"/)
 character(len=9),parameter :: dspin3(6)=(/"up       ","down     ","up-up    ","down-down","Re[up-dn]","Im[up-dn]"/)
 character(len=500) :: message0,message
!arrays
 integer :: idum(1)
 integer,allocatable :: idum1(:),idum3(:,:,:),jatom(:)
 real(dp),allocatable :: rdum2(:,:),rdum4(:,:,:,:)

! *********************************************************************

!Initializations
 natprt=natmax;if (natom==1) natprt=1
 iatom1=natom;if (pawprtvol<0) iatom1=2
 if (pawprtvol<0) natprt=natom
 allocate(jatom(natprt))
 if (natprt==1) then
  jatom(1)=1
 else if (natprt==2) then
  jatom(1)=1;jatom(2)=natom
 else if (natprt==natom) then
  do iat=1,natom
   jatom(iat)=iat
  end do
 else
  stop "pawprt: -BUG: invalid value of natprt !"
 end if

 nspden=paw_ij(1)%nspden
 nsppol=paw_ij(1)%nsppol
 antiferro=(nspden==2.and.nsppol==1)
 usepawu=(count(pawtab(:)%usepawu>0)>0)
 useexexch=(count(pawtab(:)%useexexch>0)>0)

 write(message, '(2a)' ) ch10,&
& ' ==== Results concerning PAW augmentation regions ===='
 call wrtout(ab_out,message,'COLL')
 call wrtout(6,message,'COLL')
 message=' '
 call wrtout(ab_out,message,'COLL')
 call wrtout(6,message,'COLL')

!Print out pseudopotential strength
!----------------------------------
 do unt=1,2
  if ((unt==1).and.(enunit==0.or.enunit==2)) then
   write(message, '(a)' ) &
&   ' Total pseudopotential strength Dij (hartree):'
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
  else if ((unt==2).and.(enunit==1.or.enunit==2)) then
   write(message, '(a)' ) &
&   ' Total pseudopotential strength Dij (eV):'
   call wrtout(ab_out,message,'COLL')
   call wrtout(6,message,'COLL')
  end if
  if (((unt==1).and.(enunit==0.or.enunit==2)).or.&
&  ((unt==2).and.(enunit==1.or.enunit==2))) then
   do iat=1,natprt
    iatom=jatom(iat)
    nspden=paw_ij(iatom)%ndij
    do ispden=1,nspden
     valmx=100._dp;if (ispden==1) valmx=-1._dp
     message='' ; message0=''
     if (natom>1.or.nspden>1) write(message0, '(a,i3)' ) ' Atom #',iatom
     if (nspden==1) write(message, '(a)' ) trim(message0)
     if (nspden==2) write(message, '(2a,i1)' ) trim(message0),' - Spin component ',ispden
     if (nspden==4) write(message, '(3a)' )    trim(message0),' - Component ',trim(dspin1(ispden+2*(nspden/4)))
     if (natom>1.or.nspden>1) then
      call wrtout(ab_out,message,'COLL')
      call wrtout(6,message,'COLL')
     end if
     if (nspden/=4.or.ispden<=2) then
      call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,paw_ij(iatom)%cplex_dij,&
&      paw_ij(iatom)%lmn_size,2,-1,idum,0,pawprtvol,idum,valmx,unt)
     else
      call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,paw_ij(iatom)%cplex_dij,&
&      paw_ij(iatom)%lmn_size,2,-1,idum,0,pawprtvol,idum,valmx,unt,&
&      asym_ij=paw_ij(iatom)%dij(:,7-ispden))
     end if
    end do
   end do
  end if
  message=' '
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
 end do

!Print out SYMMETRIZED occupancies of the partial waves
!------------------------------------------------------
 write(message, '(a)' ) ' Augmentation waves occupancies Rhoij:'
 call wrtout(ab_out,message,'COLL')
 call wrtout(6,message,'COLL')
 nspden=paw_ij(1)%nspden
 do iat=1,natprt
  iatom=jatom(iat)
  do ispden=1,nspden
   valmx=25._dp;if (ispden==1) valmx=-1._dp
   message='' ; message0=''
   if (natom>1.or.nspden>1) write(message0, '(a,i3)' ) ' Atom #',iatom
   if (nspden==1) write(message, '(a)' ) trim(message0)
   if (nspden==2) write(message, '(2a,i1)' ) trim(message0),' - Spin component ',ispden
   if (nspden==4) write(message, '(3a)' )    trim(message0),' - Component ',dspin2(ispden+2*(nspden/4))
   if (natom>1.or.nspden>1) then
    call wrtout(ab_out,message,'COLL')
    call wrtout(6,message,'COLL')
   end if
   call print_ij(pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&   pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,2,-1,idum,1,pawprtvol,&
&   pawrhoij(iatom)%rhoijselect(:),valmx,1)
  end do
 end do
 message=' '
 call wrtout(ab_out,message,'COLL')
 call wrtout(6,message,'COLL')

!PAW+U or local exact-exchange: print out +U components of occupancies
!-------------------------------------------------------------------------------
 if (usepawu.or.useexexch) then
  if(useexexch) write(message, '(a)' ) &
&  ' "Local exact-exchange" part of augmentation waves occupancies Rhoij:'
  if(usepawu) write(message, '(a)' ) &
&  ' "PAW+U" part of augmentation waves occupancies Rhoij:'
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
  valmx=-1._dp
  nspden=paw_ij(iatom)%nspden
  do iatom=1,natom
   itypat=typat(iatom)
   ll=-1;llp=-1
   if (pawtab(itypat)%usepawu>0) ll=pawtab(itypat)%lpawu
   if (pawtab(itypat)%useexexch>0) llp=pawtab(itypat)%lexexch
   if (ll/=llp.and.ll/=-1.and.llp/=-1) stop "pawprt: lpawu/=lexexch forbidden !"
   ll=max(ll,llp)
   if (ll>=0) then
    do ispden=1,nspden
     message='' ; message0=''
     write(message0, '(a,i3,a,i1,a)') ' Atom #',iatom,' - L=',ll,' ONLY'
     if (nspden==1) write(message, '(a)' ) trim(message0)
     if (nspden==2) write(message, '(2a,i1)' ) trim(message0),' - Spin component ',ispden
     if (nspden==4) write(message, '(3a)' )    trim(message0),' - Component ',dspin2(ispden+2*(nspden/4))
     call wrtout(ab_out,message,'COLL')
     call wrtout(6,message,'COLL')
     call print_ij(pawrhoij(iatom)%rhoijp(:,ispden),pawrhoij(iatom)%nrhoijsel,&
&     pawrhoij(iatom)%cplex,pawrhoij(iatom)%lmn_size,2,ll,indlmn(1,1:pawtab(itypat)%lmn_size,itypat),&
&     1,pawprtvol,pawrhoij(iatom)%rhoijselect(:),valmx,1)
    end do
   end if
  end do
  message=' '
  call wrtout(ab_out,message,'COLL')
  call wrtout(6,message,'COLL')
 end if

!PAW+U: print out occupations for correlated orbitals
!----------------------------------------------------
 if (usepawu) then
  write(message, '(3a)' ) &
&  '---------- LDA+U DATA --------------------------------------------------- ',ch10
  call wrtout(6,  message,'COLL')
  call wrtout(ab_out,  message,'COLL')
  nspden=paw_ij(1)%nspden
  do iatom=1,natom
   itypat=typat(iatom);ll=pawtab(itypat)%lpawu
   if ((ll>=0).and.(pawtab(itypat)%usepawu>0)) then
    write(message,fmt='(a,i5,a,i4,a)') "====== For Atom", iatom,&
&    ", occupations for correlated orbitals. lpawu =",ll,ch10
    call wrtout(6,message,'COLL')
    call wrtout(ab_out,message,'COLL')
    if(nspden==2) then
     do ispden=1,nspden
      write(message,fmt='(a,i4,a,i3,a,f10.5)') "Atom", iatom,&
&      ". Occ. for lpawu and for spin",ispden," =",paw_ij(iatom)%nocctot(ispden)
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
     end do
     write(message,fmt='(a,i4,a,2x,e15.8)') "=> On atom",iatom,&
&     ",  local Mag. for lpawu is  ",paw_ij(iatom)%nocctot(2)-paw_ij(iatom)%nocctot(1)
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,  message,'COLL')
    end if
    if(nspden==4) then
     mx= two*paw_ij(iatom)%nocctot(3)
     my=-two*paw_ij(iatom)%nocctot(4)
     mz=paw_ij(iatom)%nocctot(1)-paw_ij(iatom)%nocctot(2)
     ntot=paw_ij(iatom)%nocctot(1)+paw_ij(iatom)%nocctot(2)
     mnorm=sqrt(mx*mx+my*my+mz*mz)
     write(message,'(a,i4,a,2x,e15.8)') "=> On atom",iatom,", for  lpawu, local Mag. x is  ",mx
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,  message,'COLL')
     write(message,'(14x,a,2x,e15.8)') "              local Mag. y is  ",my
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,  message,'COLL')
     write(message,'(14x,a,2x,e15.8)') "              local Mag. z is  ",mz
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,  message,'COLL')
     write(message,'(14x,a,2x,e15.8)') "              norm of Mag. is  ",mnorm
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,  message,'COLL')
     write(message,fmt='(14x,a,2x,f10.5)') "              occ. for spin up is= ",half*(ntot+mnorm)
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     write(message,fmt='(14x,a,2x,f10.5)') "              occ. for spin dn is= ",half*(ntot-mnorm)
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')
    end if
    write(message,'(3a)') ch10,"== Occupation matrix for correlated orbitals:",ch10
    call wrtout(6,message,'COLL')
    call wrtout(ab_out,  message,'COLL')
    do ispden=1,nspden
     if (nspden==1) write(message,fmt='(a)')   "Up component only..."
     if (nspden==2) write(message,fmt='(a,i3)')"Occupation matrix for spin",ispden
     if (nspden==4) write(message,fmt='(2a)')  "Occupation matrix for component ",trim(dspin3(ispden+2*(nspden/4)))
     call wrtout(6,message,'COLL'); call wrtout(ab_out,  message,'COLL')
     do im1=1,ll*2+1
      write(message,fmt='(12(1x,9(1x,f10.5)))') &
&      (paw_ij(iatom)%noccmmp(im1,im2,ispden),im2=1,ll*2+1)
      call wrtout(6,message,'COLL')
      call wrtout(ab_out,message,'COLL')
     end do
     write(message, '(2a)' ) ch10,' '
     call wrtout(6,message,'COLL')
     call wrtout(ab_out,message,'COLL')
    end do
   end if
  end do
 end if

!Exact exchange: print out occupations for correlated orbitals
!-------------------------------------------------------------
 if (useexexch) then
  write(message, '(3a)' ) &
&  '---------- Exact Exchange --------------------------------------------------- ',ch10
  call wrtout(ab_out,message,'COLL')
  nspden=paw_ij(1)%nspden
  do iatom=1,natom
   itypat=typat(iatom);ll=pawtab(itypat)%lexexch
   if (ll>=0.and.pawtab(itypat)%useexexch>0) then
    allocate(paw_ij(iatom)%noccmmp(2*ll+1,2*ll+1,nspden))
    allocate(paw_ij(iatom)%nocctot(nspden))
   end if
  end do
  call setnoccmmp(1,0,rdum4,0,0,idum3,natom,0,nspden,1,nsppol,0,&
&  ntypat,paw_ij,pawang_dum,pawprtvol,pawrhoij,pawtab,rdum2,idum1,typat,1,0)
  do iatom=1,natom
    itypat=typat(iatom);ll=pawtab(itypat)%lexexch
    if ((ll>=0).and.(pawtab(itypat)%useexexch>0)) then
    write(message,fmt='(a,i5,a,i4,a)') "====== For Atom",iatom,&
&    ", occupations for correlated orbitals. l =",ll,ch10
    call wrtout(ab_out,message,'COLL')
    do ispden=1,paw_ij(iatom)%nspden
     if (nspden==1) write(message,fmt='(a)')   "Up component only..."
     if (nspden==2) write(message,fmt='(a,i3)')"Occupation matrix for spin",ispden
     if (nspden==4) write(message,fmt='(2a)')  "Occupation matrix for component ",trim(dspin2(ispden+2*(nspden/4)))
     call wrtout(ab_out,message,'COLL')
     do im1=1,ll*2+1
      write(message,fmt='(12(1x,9(1x,f10.5)))') (paw_ij(iatom)%noccmmp(im1,im2,ispden),im2=1,ll*2+1)
      call wrtout(ab_out,message,'COLL')
     end do
     write(message, '(a)' ) ' '
     call wrtout(ab_out,message,'COLL')
    end do
    deallocate(paw_ij(iatom)%noccmmp,paw_ij(iatom)%nocctot)
   end if
  end do
 end if

 message=' '
 call wrtout(ab_out,message,'COLL')
 call wrtout(6,message,'COLL')
 deallocate(jatom)

end subroutine pawprt

!!***
