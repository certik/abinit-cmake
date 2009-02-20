!{\src2tex{textfont=tt}}
!!****f* ABINIT/symrhoij
!! NAME
!! symrhoij
!!
!! FUNCTION
!! Symmetrize rhoij quantities (augmentation occupancies) and/or gradients
!! Compute also rhoij residuals (new-old values of rhoij and gradients)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  choice=select then type of rhoij gradients to symmetrize.
!!         choice=1 => no gradient
!!         choice=2 => gradient with respect to atomic position(s)
!!               =3 => a gradient with respect to strain(s)
!!               =4 => 2nd gradient with respect to atomic position(s)
!!               =23=> a gradient with respect to atm. pos. and strain(s)
!!               =24=> 1st and 2nd gradient with respect to atomic position(s)
!!  indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  lmnmax=maximum number of PAW radial wavefunctions
!!  natom=number of atoms in cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  optrhoij= 1 if rhoij quantities have to be symmetrized
!!  pawrhoij(natom)%lmn_size=number of (l,m,n) elements for the paw basis
!!  pawrhoij(natom)%nspden=number of spin-density components
!!  pawrhoij(natom)%nsppol=number of independant spin-density components
!!  pawrhoij(natom)%rhoij_(lmn2_size,nspden)=non-symetrized paw rhoij quantities
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!  typat(natom)=type for each atom
!!
!! OUTPUT
!!  if (optrhoij==1)
!!    pawrhoij(natom)%nrhoijsel=number of non-zero values of rhoij
!!    pawrhoij(natom)%rhoijp(lmn2_size,nspden)=symetrized paw rhoij quantities in PACKED STORAGE (only non-zero values)
!!    pawrhoij(natom)%rhoijres(lmn2_size,nspden)=paw rhoij quantities residuals (new values - old values)
!!    pawrhoij(natom)%rhoijselect(lmn2_size)=select the non-zero values of rhoij
!!
!! SIDE EFFECTS
!!  if (pawrhoij(:)%ngrhoij>0) (equivalent to choice>1)
!!    At input:
!!    pawrhoij(natom)%grhoij(ngrhoij,lmn2_size,nspden)=non-symetrized gradients of rhoij
!!    At output:
!!    pawrhoij(natom)%grhoij(ngrhoij,lmn2_size,nspden)=symetrized gradients of rhoij
!!
!! NOTES
!!  This file was initially inspirated by the file SymWij.f
!!  by N.A.W. Holzwarth and A. Tackett for the code pwpaw
!!
!! PARENTS
!!      forstrnps,vtorho
!!
!! CHILDREN
!!      leave_new,print_ij,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symrhoij(choice,indlmn,indsym,lmnmax,natom,nsym,ntypat,optrhoij,&
&                   pawang,pawprtvol,pawrhoij,symafm,symrec,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: choice,lmnmax,natom,nsym,ntypat,optrhoij,pawprtvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),indsym(4,nsym,natom)
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij(natom)

!Local variables ---------------------------------------
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
!scalars
 integer :: at_indx,iafm,iatom,idum1,idum2,il,il0,ilmn,iln,iln0,ilpm,indexi
 integer :: indexii,indexj,indexjj,indexjj0,indexk,irhoij,irot,ishift2,ishift3
 integer :: ishift4,ispden,itypat,j0lmn,jj,jl,jl0,jlmn,jln,jln0,jlpm,klmn
 integer :: lmn_size,mi,mj,mu,mua,mub,mushift,natinc,ngrhoij,nselect,nsploop,nu
 integer :: nushift
 real(dp) :: ro,sum1,sym1,sym2,sym3,syma,zarot2
 logical :: antiferro
 character(len=500) :: message
!arrays
 integer,parameter :: alpha(6)=(/1,2,3,3,3,2/),beta(6)=(/1,2,3,2,1,1/)
 integer :: nsym_used(2)
 integer,allocatable :: idum(:)
 real(dp) :: paw_sumrho(2),work1(3,3)
 real(dp),allocatable :: rotgr(:,:),sumgr(:)
!no_abirules
 type grhoij_at
  real(dp),pointer :: grhoij(:,:,:)
 end type
 type(grhoij_at) :: grtmp(natom)

! *********************************************************************

!Number of independant spin components
 nsploop=pawrhoij(1)%nsppol
 if (pawrhoij(1)%nspden==4) nsploop=4

!Symetrization occurs only when nsym>1
 if (nsym>1) then

! Test: consistency between choice and ngrhoij
  ngrhoij=pawrhoij(1)%ngrhoij
  if ((choice==1.and.ngrhoij/=0) .or.(choice==2.and.ngrhoij/=3).or. &
&  (choice==3.and.ngrhoij/=6).or.(choice==23.and.ngrhoij/=9).or. &
&  (choice==4.and.ngrhoij/=6).or.(choice==24.and.ngrhoij/=9) ) then
   write(message, '(a,a,a,a)' )ch10,&
&   ' symrhoij : BUG -',ch10,&
&   '  Inconsistency between variables choice and ngrhoij !'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  end if

! Antiferro case ?
  antiferro=(pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1)

! Several inits
  if (choice>1) then
   allocate(sumgr(ngrhoij))
   if (antiferro) then
    allocate(rotgr(ngrhoij,2))
   else
    allocate(rotgr(ngrhoij,1))
   end if
   ishift2=0;ishift3=0;ishift4=0
   if (choice==23) ishift2=6
   if (choice==24) ishift4=3
!  Have to make a temporary copy of grhoij
   do iatom=1,natom
    idum1=pawrhoij(iatom)%lmn2_size;idum2=pawrhoij(iatom)%nspden
    allocate(grtmp(iatom)%grhoij(ngrhoij,idum1,idum2))
    grtmp(iatom)%grhoij(1:ngrhoij,1:idum1,1:idum2)=pawrhoij(iatom)%grhoij(1:ngrhoij,1:idum1,1:idum2)
   end do
  end if

! Loops over atoms and spin components
! ------------------------------------
  do iatom=1,natom
   itypat=typat(iatom)
   lmn_size=pawrhoij(iatom)%lmn_size

   nselect=0
   do ispden=1,nsploop

!   Store old -rhoij in residual
    if (optrhoij==1) then
     pawrhoij(iatom)%rhoijres(:,ispden)=zero
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      pawrhoij(iatom)%rhoijres(klmn,ispden)=-pawrhoij(iatom)%rhoijp(irhoij,ispden)
     end do
     if (antiferro) then
      pawrhoij(iatom)%rhoijres(:,2)=zero
      do irhoij=1,pawrhoij(iatom)%nrhoijsel
       klmn=pawrhoij(iatom)%rhoijselect(irhoij)
       pawrhoij(iatom)%rhoijres(klmn,2)=-pawrhoij(iatom)%rhoijp(irhoij,2)
      end do
     end if
    end if

!   Loops over (il,im) and (jl,jm)
!   ------------------------------
    jl0=-1;jln0=-1;indexj=1
    do jlmn=1,lmn_size
     jl=indlmn(1,jlmn,itypat)
     jlpm=1+jl+indlmn(2,jlmn,itypat)
     jln=indlmn(5,jlmn,itypat)
     if (jln/=jln0) indexj=indexj+2*jl0+1
     j0lmn=jlmn*(jlmn-1)/2
     il0=-1;iln0=-1;indexi=1
     do ilmn=1,jlmn
      il=indlmn(1,ilmn,itypat)
      ilpm=1+il+indlmn(2,ilmn,itypat)
      iln=indlmn(5,ilmn,itypat)
      if (iln/=iln0) indexi=indexi+2*il0+1
      klmn=j0lmn+ilmn

      nsym_used(:)=0
      paw_sumrho(:)=zero
      if (choice>1) rotgr(:,:)=zero

!     Loop over symmetries
!     --------------------
      do irot=1,nsym

       if ((symafm(irot)/=1).and.(.not.antiferro)) cycle
       iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2


       nsym_used(iafm)=nsym_used(iafm)+1
       at_indx=indsym(4,irot,iatom)

!      Accumulate values over (mi,mj)
!      ------------------------------
       if (choice>1) sumgr=zero
       do mj=1,2*jl+1
        indexjj=indexj+mj;indexjj0=indexjj*(indexjj-1)/2
        do mi=1,2*il+1
         indexii=indexi+mi
         if (indexii<=indexjj) then
          indexk=indexjj0+indexii
         else
          indexk=indexii*(indexii-1)/2+indexjj
         end if
         zarot2=pawang%zarot(mi,ilpm,il+1,irot)&
&         *pawang%zarot(mj,jlpm,jl+1,irot)

         if (optrhoij==1) &
&         paw_sumrho(iafm)=paw_sumrho(iafm)+zarot2*pawrhoij(at_indx)%rhoij_(indexk,ispden)

         if (choice>1) then
          do mu=1,ngrhoij
           sumgr(mu)=sumgr(mu)+zarot2*grtmp(at_indx)%grhoij(mu,indexk,ispden)
          end do
         end if

        end do
       end do

!      Accumulate values over symmetries
!      ---------------------------------
!      ===== Contributions to derivative vs atm. pos. ====
       if (choice==2.or.choice==23.or.choice==24) then
        do nu=1,3
         nushift=nu+ishift2
         do mu=1,3
          mushift=mu+ishift2
          rotgr(mushift,iafm)=rotgr(mushift,iafm)+dble(symrec(mu,nu,irot))*sumgr(nushift)
         end do
        end do
       end if
!      ===== Contribution to derivative vs strain ====
       if (choice==3.or.choice==23) then
        work1(1,1)=sumgr(1+ishift3);work1(2,2)=sumgr(2+ishift3)
        work1(3,3)=sumgr(3+ishift3);work1(2,3)=sumgr(4+ishift3)
        work1(1,3)=sumgr(5+ishift3);work1(1,2)=sumgr(6+ishift3)
        work1(3,1)=work1(1,3);work1(3,2)=work1(2,3)
        work1(2,1)=work1(1,2)
        do mu=1,6
         mushift=mu+ishift3
         mua=alpha(mu);mub=beta(mu)
         sum1=zero
         sym1=dble(symrec(mub,1,irot))
         sym2=dble(symrec(mub,2,irot))
         sym3=dble(symrec(mub,3,irot))
         do nu=1,3
          syma=dble(symrec(mua,nu,irot))
          sum1=sum1+syma*(work1(nu,1)*sym1+work1(nu,2)*sym2+work1(nu,3)*sym3)
         end do
         rotgr(mushift,iafm)=rotgr(mushift,iafm)+sum1
        end do
       end if
!      ===== Contributions to second derivative vs atm. pos. ====
       if (choice==4.or.choice==24) then
        work1(1,1)=sumgr(1+ishift4);work1(2,2)=sumgr(2+ishift4)
        work1(3,3)=sumgr(3+ishift4);work1(2,3)=sumgr(4+ishift4)
        work1(1,3)=sumgr(5+ishift4);work1(1,2)=sumgr(6+ishift4)
        work1(3,1)=work1(1,3);work1(3,2)=work1(2,3)
        work1(2,1)=work1(1,2)
        do mu=1,6
         mushift=mu+ishift4
         mua=alpha(mu);mub=beta(mu)
         sum1=zero
         sym1=dble(symrec(mub,1,irot))
         sym2=dble(symrec(mub,2,irot))
         sym3=dble(symrec(mub,3,irot))
         do nu=1,3
          syma=dble(symrec(mua,nu,irot))
          sum1=sum1+syma*(work1(nu,1)*sym1+work1(nu,2)*sym2+work1(nu,3)*sym3)
         end do
         rotgr(mushift,iafm)=rotgr(mushift,iafm)+sum1
        end do
       end if

      end do ! End loop over symmetries

!     Select non-zero elements of rhoij
!     ---------------------------------
      if (optrhoij==1) then

       ro=paw_sumrho(1)/nsym_used(1)
       if (abs(ro)>tol10) then
        pawrhoij(iatom)%rhoijres(klmn,ispden)=pawrhoij(iatom)%rhoijres(klmn,ispden)+ro
        pawrhoij(iatom)%rhoijp(klmn,ispden)=ro
       else
        pawrhoij(iatom)%rhoijp(klmn,ispden)=zero
       end if
!      Antiferro case: fill also down component
       if (antiferro.and.nsym_used(2)>0) then
        ro=paw_sumrho(2)/nsym_used(2)
        if (abs(ro)>tol10) then
         pawrhoij(iatom)%rhoijres(klmn,2)=pawrhoij(iatom)%rhoijres(klmn,2)+ro
         pawrhoij(iatom)%rhoijp(klmn,2)=ro
        else
         pawrhoij(iatom)%rhoijp(klmn,2)=zero
        end if
       end if

       if (ispden==nsploop) then
        if (any(abs(pawrhoij(iatom)%rhoijp(klmn,:))>tol10)) then
         nselect=nselect+1
         pawrhoij(iatom)%rhoijselect(nselect)=klmn
         do jj=1,pawrhoij(iatom)%nspden
          pawrhoij(iatom)%rhoijp(nselect,jj)=pawrhoij(iatom)%rhoijp(klmn,jj)
         end do
        end if
       end if

      end if ! optrhoij==1

!     Normalization of gradients
      if (choice>1) then
       do mu=1,ngrhoij
        pawrhoij(iatom)%grhoij(mu,klmn,ispden)=rotgr(mu,1)/nsym_used(1)
       end do
       if (antiferro.and.nsym_used(2)>0) then
        do mu=1,ngrhoij
         pawrhoij(iatom)%grhoij(mu,klmn,2)=rotgr(mu,2)/nsym_used(2)
        end do
       end if
      end if

      il0=il;iln0=iln  ! End loops over (il,im) and (jl,jm)
     end do
     jl0=jl;jln0=jln
    end do

   end do  ! End loop over ispden

!  Store number of non-zero values of rhoij
   if (optrhoij==1) pawrhoij(iatom)%nrhoijsel=nselect

  end do ! End loop over iatom

  if (choice>1) then
   do iatom=1,natom
    deallocate(grtmp(iatom)%grhoij)
   end do
   deallocate(sumgr,rotgr)
  end if

 else  ! nsym>1

! *********************************************************************
! If nsym==1, only copy rhoij_ into rhoij
! also has to fill rhoijselect array

  if(pawrhoij(1)%nspden==2.and.pawrhoij(1)%nsppol==1) then
   write(message,'(a,a,a)') ' symrhoij : BUG -',ch10,&
&   ' In the antiferromagnetic case, nsym cannot be 1'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  end if
  if (optrhoij==1) then
   do iatom=1,natom
    pawrhoij(iatom)%rhoijres(:,:)=zero
    do ispden=1,nsploop
     do irhoij=1,pawrhoij(iatom)%nrhoijsel
      klmn=pawrhoij(iatom)%rhoijselect(irhoij)
      pawrhoij(iatom)%rhoijres(klmn,ispden)=-pawrhoij(iatom)%rhoijp(irhoij,ispden)
     end do
    end do
    nselect=0
    do klmn=1,pawrhoij(iatom)%lmn2_size
     if (any(abs(pawrhoij(iatom)%rhoij_(klmn,:))>tol10)) then
      nselect=nselect+1
      pawrhoij(iatom)%rhoijselect(nselect)=klmn
      do jj=1,pawrhoij(iatom)%nspden
       ro=pawrhoij(iatom)%rhoij_(klmn,jj)
       pawrhoij(iatom)%rhoijres(klmn,jj)=pawrhoij(iatom)%rhoijres(klmn,jj)+ro
       pawrhoij(iatom)%rhoijp(nselect,jj)=ro
      end do
     end if
    end do
    pawrhoij(iatom)%nrhoijsel=nselect
   end do
  end if

 end if

!*********************************************************************
!Printing of Rhoij

 if (optrhoij==1) then
  natinc=1;if(natom>1.and.pawprtvol>=0) natinc=natom-1
  do iatom=1,natom,natinc
   if (abs(pawprtvol)>=1) then
    write(message, '(4a,i3,a)') ch10," PAW TEST:",ch10,&
&    ' ====== Values of RHOIJ in symrhoij (iatom=',iatom,') ======'
    call wrtout(6,message,'COLL')
   end if
   do ispden=1,pawrhoij(iatom)%nspden
    if (abs(pawprtvol)>=1) then
     write(message, '(3a)') '   Component ',trim(dspin(ispden+2*(pawrhoij(iatom)%nspden/4))),':'
    else
     if (pawrhoij(iatom)%nspden/=4) write(message, '(a,a,i3,a,i1,a)') ch10,&
&     ' *********** Rhoij (atom ',iatom,', ispden=',ispden,') **********'
     if (pawrhoij(iatom)%nspden==4) write(message, '(a,a,i3,3a)') ch10,&
&     ' *********** Rhoij (atom ',iatom,' - ',trim(dspin(ispden+2*(pawrhoij(iatom)%nspden/4))),') **********'
    end if
    call wrtout(6,message,'COLL')
    call print_ij(pawrhoij(iatom)%rhoijp(:,ispden),&
&    pawrhoij(iatom)%nrhoijsel,1,&
&    pawrhoij(iatom)%lmn_size,1,-1,idum,1,pawprtvol,&
&    pawrhoij(iatom)%rhoijselect(:),&
&    10.d0*dble(3-2*ispden),1)
   end do
  end do
  message=''
  call wrtout(6,message,'COLL')
 end if

end subroutine symrhoij
!!***
