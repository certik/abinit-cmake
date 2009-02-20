!{\src2tex{textfont=tt}}
!!****f* ABINIT/symrhoij3
!! NAME
!! symrhoij3
!!
!! FUNCTION
!! Symmetrize rhoij1 quantities (1st-order augmentation occupancies)
!! Compute also rhoij1 residuals (new-old values of rhoij1)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=index of perturbation
!!  lmnmax=maximum number of PAW radial wavefunctions
!!  natom=number of atoms in cell
!!  nrhoij1=size of pawrhoij1 (1 if atomic displ. perturb., else natom)
!!  nsym1=number of symmetry elements in space group consistent with perturbation
!!  ntypat=number of types of atoms in unit cell.
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  symaf1(nsym1)=(anti)ferromagnetic part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  typat(natom)=type for each atom
!!
!! SIDE EFFECTS
!!  pawrhoij1(nrhoij1) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!     At input:
!!           %lmn_size=number of (l,m,n) elements for the paw basis
!!           %nspden=number of spin-density components
!!           %nsppol=number of independant spin-density components
!!           %rhoij_(lmn2_size,nspden)=non-symetrized rhoij1 quantities
!!    At output:
!!           %nrhoijsel=number of non-zer oelements
!!           %rhoijp(lmn2_size,nspden)=rhoij "packed", only non-zero elements
!!           %rhoijselect(lmn2_size)=indexes of non-zero elements
!!
!! PARENTS
!!      vtorho3
!!
!! CHILDREN
!!      leave_new,print_ij,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symrhoij3(indlmn,indsy1,ipert,lmnmax,natom,nrhoij1,nsym1,ntypat,&
&                    pawang,pawprtvol,pawrhoij1,symaf1,symrc1,typat)

 use defs_basis
 use defs_datatypes


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: ipert,lmnmax,natom,nrhoij1,nsym1,ntypat,pawprtvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),indsy1(4,nsym1,natom)
 integer,intent(in) :: symaf1(nsym1),symrc1(3,3,nsym1),typat(natom)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(nrhoij1)

!Local variables ---------------------------------------
!scalars
 integer :: at_indx,cplex,iafm,iatom,iatm,il,il0,ilmn,iln,iln0,ilpm,indexi
 integer :: indexii,indexj,indexjj,indexjj0,indexk,irhoij,irot
 integer :: ispden,itypat,j0lmn,jj,jl,jl0,jlmn,jln,jln0,jlpm,klmn
 integer :: lmn_size,mi,mj,mu,natinc,nselect,nsploop,nu
 real(dp) :: ro(2),zarot2
 logical :: antiferro
 character(len=500) :: message
!arrays
 integer :: nsym_used(2)
 integer,allocatable :: idum(:)
 real(dp) :: paw_sumrho(2,2)
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)

! *********************************************************************

 if(ipert>natom) then
  write(message, '(a,a,a,a)' )ch10,&
&  ' symrhoij3 : ERROR -',ch10,&
&  '  ipert>natom not yet implemented !'
  call wrtout(6,message,'PERS')
  call leave_new('PERS')
 end if

!Number of independant spin components
 nsploop=pawrhoij1(1)%nsppol
 if (pawrhoij1(1)%nspden==4) nsploop=4

!Symetrization occurs only when nsym>1
 if (nsym1>1) then

! Antiferro case ?
  antiferro=(pawrhoij1(1)%nspden==2.and.pawrhoij1(1)%nsppol==1)

! Loops over atoms and spin components
! ------------------------------------
  do iatm=1,nrhoij1
   iatom=iatm;if (nrhoij1==1) iatom=ipert
   itypat=typat(iatom)
   lmn_size=pawrhoij1(iatm)%lmn_size
   cplex=pawrhoij1(iatm)%cplex

   nselect=0
   do ispden=1,nsploop

!   Store old -rhoij1 in residual
    pawrhoij1(iatm)%rhoijres(:,ispden)=zero
    do irhoij=1,pawrhoij1(iatm)%nrhoijsel
     klmn=pawrhoij1(iatm)%rhoijselect(irhoij)
     pawrhoij1(iatm)%rhoijres(klmn,ispden)=-pawrhoij1(iatm)%rhoijp(irhoij,ispden)
    end do
    if (antiferro) then
     pawrhoij1(iatm)%rhoijres(:,2)=zero
     do irhoij=1,pawrhoij1(iatm)%nrhoijsel
      klmn=pawrhoij1(iatm)%rhoijselect(irhoij)
      pawrhoij1(iatm)%rhoijres(klmn,2)=-pawrhoij1(iatm)%rhoijp(irhoij,2)
     end do
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
      paw_sumrho(:,:)=zero

!     Loop over symmetries
!     --------------------
      do irot=1,nsym1

       if ((symaf1(irot)/=1).and.(.not.antiferro)) cycle
       iafm=1;if ((antiferro).and.(symaf1(irot)==-1)) iafm=2
       at_indx=min(indsy1(4,irot,iatom),nrhoij1)
       nsym_used(iafm)=nsym_used(iafm)+1

!      Accumulate values over (mi,mj)
!      ------------------------------
       do mj=1,2*jl+1
        indexjj=indexj+mj;indexjj0=indexjj*(indexjj-1)/2
        do mi=1,2*il+1
         indexii=indexi+mi
         if (indexii<=indexjj) then
          indexk=indexjj0+indexii
         else
          indexk=indexii*(indexii-1)/2+indexjj
         end if
         zarot2=pawang%zarot(mi,ilpm,il+1,irot)*pawang%zarot(mj,jlpm,jl+1,irot)

         if (cplex==1) then
          paw_sumrho(1,iafm)=paw_sumrho(1,iafm)+zarot2*pawrhoij1(at_indx)%rhoij_(indexk,ispden)
         else
          paw_sumrho(1,iafm)=paw_sumrho(1,iafm)+zarot2*pawrhoij1(at_indx)%rhoij_(2*indexk-1,ispden)
          paw_sumrho(2,iafm)=paw_sumrho(2,iafm)+zarot2*pawrhoij1(at_indx)%rhoij_(2*indexk  ,ispden)
         end if

        end do
       end do

      end do ! End loop over symmetries

!     Select non-zero elements of rhoij1
      if (cplex==1) then
       ro(1)=paw_sumrho(1,1)/nsym_used(1)
       if (abs(ro(1))>tol10) then
        pawrhoij1(iatm)%rhoijres(klmn,ispden)=pawrhoij1(iatm)%rhoijres(klmn,ispden)+ro(1)
        pawrhoij1(iatm)%rhoijp(klmn,ispden)=ro(1)
       else
        pawrhoij1(iatm)%rhoijp(klmn,ispden)=zero
       end if
      else
       ro(1)=paw_sumrho(1,1)/nsym_used(1)
       ro(2)=paw_sumrho(2,1)/nsym_used(1)
       if (any(abs(ro(1:2))>tol10)) then
        pawrhoij1(iatm)%rhoijres(2*klmn-1,ispden)=pawrhoij1(iatm)%rhoijres(2*klmn-1,ispden)+ro(1)
        pawrhoij1(iatm)%rhoijres(2*klmn  ,ispden)=pawrhoij1(iatm)%rhoijres(2*klmn  ,ispden)+ro(2)
        pawrhoij1(iatm)%rhoijp(2*klmn-1,ispden)=ro(1)
        pawrhoij1(iatm)%rhoijp(2*klmn  ,ispden)=ro(2)
       else
        pawrhoij1(iatm)%rhoijp(2*klmn-1,ispden)=zero
        pawrhoij1(iatm)%rhoijp(2*klmn  ,ispden)=zero
       end if
      end if
!     Antiferro case: fill also down component
      if (antiferro.and.nsym_used(2)>0) then
       if (cplex==1) then
        ro(1)=paw_sumrho(1,2)/nsym_used(2)
        if (abs(ro(1))>tol10) then
         pawrhoij1(iatm)%rhoijres(klmn,2)=pawrhoij1(iatm)%rhoijres(klmn,2)+ro(1)
         pawrhoij1(iatm)%rhoijp(klmn,2)=ro(1)
        else
         pawrhoij1(iatm)%rhoijp(klmn,2)=zero
        end if
       else
        ro(1)=paw_sumrho(1,2)/nsym_used(2)
        ro(2)=paw_sumrho(2,2)/nsym_used(2)
        if (any(abs(ro(1:2))>tol10)) then
         pawrhoij1(iatm)%rhoijres(2*klmn-1,2)=pawrhoij1(iatm)%rhoijres(2*klmn  ,2)+ro(1)
         pawrhoij1(iatm)%rhoijres(2*klmn  ,2)=pawrhoij1(iatm)%rhoijres(2*klmn-1,2)+ro(2)
         pawrhoij1(iatm)%rhoijp(2*klmn-1,2)=ro(1)
         pawrhoij1(iatm)%rhoijp(2*klmn  ,2)=ro(2)
        else
         pawrhoij1(iatm)%rhoijp(2*klmn-1,2)=zero
         pawrhoij1(iatm)%rhoijp(2*klmn  ,2)=zero
        end if
       end if
      end if

      if (ispden==nsploop) then
       if (cplex==1) then
        if (any(abs(pawrhoij1(iatm)%rhoijp(klmn,:))>tol10)) then
         nselect=nselect+1
         pawrhoij1(iatm)%rhoijselect(nselect)=klmn
         do jj=1,pawrhoij1(iatm)%nspden
          pawrhoij1(iatm)%rhoijp(nselect,jj)=pawrhoij1(iatm)%rhoijp(klmn,jj)
         end do
        end if
       else
        if (any(abs(pawrhoij1(iatm)%rhoijp(2*klmn-1:2*klmn,:))>tol10)) then
         nselect=nselect+1
         pawrhoij1(iatm)%rhoijselect(nselect)=klmn
         do jj=1,pawrhoij1(iatm)%nspden
          pawrhoij1(iatm)%rhoijp(2*nselect-1,jj)=pawrhoij1(iatm)%rhoijp(2*klmn-1,jj)
          pawrhoij1(iatm)%rhoijp(2*nselect  ,jj)=pawrhoij1(iatm)%rhoijp(2*klmn  ,jj)
         end do
        end if
       end if
      end if

      il0=il;iln0=iln  ! End loops over (il,im) and (jl,jm)
     end do
     jl0=jl;jln0=jln
    end do

   end do  ! End loop over ispden

!  Store number of non-zero values of rhoij1
   pawrhoij1(iatm)%nrhoijsel=nselect

  end do  ! End loop over iatom

 else  ! nsym1>1

! *********************************************************************
! If nsym1==1, only copy rhoij_ into rhoij
! also has to fill rhoijselect array

  if(pawrhoij1(1)%nspden==2.and.pawrhoij1(1)%nsppol==1) then
   write(message,'(a,a,a)') ' symrhoij3 : BUG -',ch10,&
&   ' In the antiferromagnetic case, nsym1 cannot be 1'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  end if
  do iatm=1,nrhoij1
   cplex=pawrhoij1(iatm)%cplex
   pawrhoij1(iatm)%rhoijres(:,:)=zero
   if (cplex==1) then
    do ispden=1,nsploop
     do irhoij=1,pawrhoij1(iatm)%nrhoijsel
      klmn=pawrhoij1(iatm)%rhoijselect(irhoij)
      pawrhoij1(iatm)%rhoijres(klmn,ispden)=-pawrhoij1(iatm)%rhoijp(irhoij,ispden)
     end do
    end do
    nselect=0
    do klmn=1,pawrhoij1(iatm)%lmn2_size
     if (any(abs(pawrhoij1(iatm)%rhoij_(klmn,:))>tol10)) then
      nselect=nselect+1
      pawrhoij1(iatm)%rhoijselect(nselect)=klmn
      do jj=1,pawrhoij1(iatm)%nspden
       ro(1)=pawrhoij1(iatm)%rhoij_(klmn,jj)
       pawrhoij1(iatm)%rhoijres(klmn,jj)=pawrhoij1(iatm)%rhoijres(klmn,jj)+ro(1)
       pawrhoij1(iatm)%rhoijp(nselect,jj)=ro(1)
      end do
     end if
    end do
   else
    do ispden=1,nsploop
     do irhoij=1,pawrhoij1(iatm)%nrhoijsel
      klmn=pawrhoij1(iatm)%rhoijselect(irhoij)
      pawrhoij1(iatm)%rhoijres(2*klmn-1,ispden)=-pawrhoij1(iatm)%rhoijp(2*irhoij-1,ispden)
      pawrhoij1(iatm)%rhoijres(2*klmn  ,ispden)=-pawrhoij1(iatm)%rhoijp(2*irhoij  ,ispden)
     end do
    end do
    nselect=0
    do klmn=1,pawrhoij1(iatm)%lmn2_size
     if (any(abs(pawrhoij1(iatm)%rhoij_(2*klmn-1:2*klmn,:))>tol10)) then
      nselect=nselect+1
      pawrhoij1(iatm)%rhoijselect(nselect)=klmn
      do jj=1,pawrhoij1(iatm)%nspden
       ro(1)=pawrhoij1(iatm)%rhoij_(2*klmn-1,jj)
       ro(2)=pawrhoij1(iatm)%rhoij_(2*klmn  ,jj)
       pawrhoij1(iatm)%rhoijres(2*klmn-1,jj)=pawrhoij1(iatm)%rhoijres(2*klmn-1,jj)+ro(1)
       pawrhoij1(iatm)%rhoijres(2*klmn  ,jj)=pawrhoij1(iatm)%rhoijres(2*klmn  ,jj)+ro(2)
       pawrhoij1(iatm)%rhoijp(2*nselect-1,jj)=ro(1)
       pawrhoij1(iatm)%rhoijp(2*nselect  ,jj)=ro(2)
      end do
     end if
    end do
   end if
   pawrhoij1(iatm)%nrhoijsel=nselect
  end do

 end if

!*********************************************************************
!Printing of Rhoij1

 natinc=1;if(nrhoij1>1.and.pawprtvol>=0) natinc=nrhoij1-1
 do iatm=1,nrhoij1,natinc
  if (abs(pawprtvol)>=1) then
   write(message, '(4a,i3,a)') ch10," PAW TEST:",ch10,&
&   ' ====== Values of RHOIJ(1) in symrhoij3 (index=',iatm,') ======'
   call wrtout(6,message,'COLL')
  end if
  do ispden=1,pawrhoij1(iatm)%nspden
   if (abs(pawprtvol)>=1) then
    write(message, '(3a)') '   Component ',trim(dspin(ispden+2*(pawrhoij1(iatm)%nspden/4))),':'
   else
    if (pawrhoij1(iatm)%nspden/=4) write(message, '(a,a,i3,a,i1,a)') ch10,&
&    ' *********** Rhoij (atom ',iatom,', ispden=',ispden,') **********'
    if (pawrhoij1(iatm)%nspden==4) write(message, '(a,a,i3,3a)') ch10,&
&    ' *********** Rhoij (atom ',iatom,' - ',trim(dspin(ispden+2*(pawrhoij1(iatm)%nspden/4))),') **********'
   end if
   call wrtout(6,message,'COLL')
   call print_ij(pawrhoij1(iatm)%rhoijp(:,ispden),&
&   pawrhoij1(iatm)%nrhoijsel,&
&   pawrhoij1(iatm)%cplex,&
&   pawrhoij1(iatm)%lmn_size,1,-1,idum,1,pawprtvol,&
&   pawrhoij1(iatm)%rhoijselect(:),-1.d0,1)
  end do
 end do
 message=''
 call wrtout(6,message,'COLL')

end subroutine symrhoij3
!!***
