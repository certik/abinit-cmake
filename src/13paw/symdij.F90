!{\src2tex{textfont=tt}}
!!****f* ABINIT/symdij
!! NAME
!! symdij
!!
!! FUNCTION
!! Symmetrize dij quantities (psp strengths)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  indlmn(6,lmnmax,ntypat)=array giving l,m,n,lm,ln,spin for i=lmn (for each atom type)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  lmnmax=maximum number of PAW radial wavefunctions
!!  natom=number of atoms in cell
!!  nsym=number of symmetry elements in space group
!!  ntypat=number of types of atoms in unit cell.
!!  paw_ij(natom)%lmn_size=number of (l,m,n) elements for the paw basis
!!  paw_ij(natom)%nspden=number of spin-density components
!!  paw_ij(natom)%dij(lmn2_size,nspden)=non-symetrized paw dij quantities
!!  pawang <type(pawang_type)>=angular mesh discretization and related data
!!  pawprtvol=control print volume and debugging output for PAW
!!  symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!                   reciprocal space primitive translations
!!  typat(natom)=type for each atom
!!
!! SIDE EFFECTS
!!    paw_ij(natom)%dij(cplex_dij*lmn2_size,nspden)=symetrized dij quantities as output
!!
!! PARENTS
!!      respfn,scfcv,sigma
!!
!! CHILDREN
!!      leave_new,print_ij,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine symdij(indlmn,indsym,lmnmax,natom,nsym,ntypat,&
&                  paw_ij,pawang,pawprtvol,symafm,symrec,typat)

 use defs_basis
 use defs_datatypes
 use m_io_tools, only : flush_unit


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_11util
!End of the abilint section

 implicit none

!Arguments ---------------------------------------------
!scalars
 integer,intent(in) :: lmnmax,natom,nsym,ntypat,pawprtvol
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),indsym(4,nsym,natom)
 integer,intent(in) :: symafm(nsym),symrec(3,3,nsym),typat(natom)
 type(paw_ij_type),intent(inout) :: paw_ij(natom)

!Local variables ---------------------------------------
!scalars
 integer :: at_indx,cplex_dij,iafm,iatom,il,il0,ilmn,iln,iln0,ilpm,indexi,indexii,indexj
 integer :: indexjj,indexjj0,indexk,indexkc,iplex,irot,irotaf,ispden,itypat,j0lmn,jj,jl,jl0
 integer :: jlmn,jln,jln0,jlpm,klmn,klmnc,lmn_size,mi,mj,natinc,nsploop
 real(dp) :: dijhatnew,dijnew,factl,zarot2
 logical :: antiferro,has_dijhat
 character(len=500) :: message
!arrays
 integer :: nsym_used(2)
 integer,allocatable :: idum(:)
 real(dp) :: sumdij(2,2),sumdijhat(2,2)
 character(len=7),parameter :: dspin(6)=(/"up     ","down   ","up-up  ","dwn-dwn","up-dwn ","dwn-up "/)
!no_abirules
  type dij_at
   real(dp),pointer :: dij(:,:)
  end type
  type(dij_at),allocatable :: tmp(:),tmphat(:)

! *********************************************************************

#if defined DEBUG_MODE
 write(message,'(a)')' symdij : enter '
 call wrtout(std_out,message,'COLL') 
 call flush_unit(std_out)
#endif

 nsploop=paw_ij(1)%nsppol
 if (paw_ij(1)%ndij==4) nsploop=4

!Symetrization occurs only when nsym>1 and nsploop/=4
 if (nsym>1.and.nsploop/=4) then

! Antiferro case ?
  antiferro=(paw_ij(1)%nspden==2.and.paw_ij(1)%nsppol==1.and.paw_ij(1)%ndij/=4)
  has_dijhat=(paw_ij(1)%has_dijhat>0)

! Have to make a temporary copy of dij
  allocate(tmp(natom));if (has_dijhat) allocate(tmphat(natom))
  do iatom=1,natom
   allocate(tmp(iatom)%dij(paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size,paw_ij(iatom)%ndij))
   tmp(iatom)%dij(:,:)=paw_ij(iatom)%dij(:,:)
   if (has_dijhat) then
    allocate(tmphat(iatom)%dij(paw_ij(iatom)%cplex_dij*paw_ij(iatom)%lmn2_size,paw_ij(iatom)%ndij))
    tmphat(iatom)%dij(:,:)=paw_ij(iatom)%dijhat(:,:)
   end if
  end do

! Loops over atoms and spin components
  do iatom=1,natom
   itypat=typat(iatom)
   lmn_size=paw_ij(iatom)%lmn_size
   cplex_dij=paw_ij(iatom)%cplex_dij
   do ispden=1,nsploop

!   Loops over (il,im) and (jl,jm)
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
      klmn=j0lmn+ilmn;klmnc=cplex_dij*(klmn-1)

      nsym_used(1:2)=0
      sumdij(1:2,1:2)=zero
      if (has_dijhat) sumdijhat(1:2,1:2)=zero

!     Loop over symmetries
      do irot=1,nsym

       if ((symafm(irot)/=1).and.(.not.antiferro)) cycle
       iafm=1;if ((antiferro).and.(symafm(irot)==-1)) iafm=2

       nsym_used(iafm)=nsym_used(iafm)+1
       at_indx=indsym(4,irot,iatom)

!      Accumulate values over (mi,mj) and symetries
       do mj=1,2*jl+1
        indexjj=indexj+mj;indexjj0=indexjj*(indexjj-1)/2
        do mi=1,2*il+1
         indexii=indexi+mi
         if (indexii<=indexjj) then
          indexk=indexjj0+indexii
         else
          indexk=indexii*(indexii-1)/2+indexjj
         end if
         indexkc=cplex_dij*(indexk-1)
         zarot2=pawang%zarot(mi,ilpm,il+1,irot)*pawang%zarot(mj,jlpm,jl+1,irot)

         do iplex=1,cplex_dij
          sumdij(iplex,iafm)=sumdij(iplex,iafm)+zarot2*tmp(at_indx)%dij(indexkc+iplex,ispden)
         end do

         if (has_dijhat) then
          do iplex=1,cplex_dij
           sumdijhat(iplex,iafm)=sumdijhat(iplex,iafm)+zarot2*tmphat(at_indx)%dij(indexkc+iplex,ispden)
          end do

         end if
        end do
       end do
      end do ! End loop over symmetries

!     Store new values of dij
      do iplex=1,cplex_dij
       dijnew=sumdij(iplex,1)/nsym_used(1)
       if (abs(dijnew)>tol10) then
        paw_ij(iatom)%dij(klmnc+iplex,ispden)=dijnew
       else
        paw_ij(iatom)%dij(klmnc+iplex,ispden)=zero
       end if
       if (has_dijhat) then
        dijhatnew=sumdijhat(IPLEX,1)/nsym_used(1)
        if (abs(dijhatnew)>tol10) then
         paw_ij(iatom)%dijhat(klmnc+iplex,ispden)=dijhatnew
        else
         paw_ij(iatom)%dijhat(klmnc+iplex,ispden)=zero
        end if
       end if
      end do

!     Antiferromagnetic case: has to fill up "down" component of dij
      if (antiferro.and.nsym_used(2)>0) then
       do iplex=1,cplex_dij
        dijnew=sumdij(iplex,2)/nsym_used(2)
        if (abs(dijnew)>tol10) then
         paw_ij(iatom)%dij(klmnc+iplex,2)=dijnew
        else
         paw_ij(iatom)%dij(klmnc+iplex,2)=zero
        end if
        if (has_dijhat) then
         dijhatnew=sumdijhat(iplex,2)/nsym_used(2)
         if (abs(dijhatnew)>tol10) then
          paw_ij(iatom)%dijhat(klmnc+iplex,2)=dijhatnew
         else
          paw_ij(iatom)%dijhat(klmnc+iplex,2)=zero
         end if
        end if
       end do
      end if

      il0=il;iln0=iln  ! End loops over (il,im) and (jl,jm)
     end do
     jl0=jl;jln0=jln
    end do
   end do ! ispden
  end do ! iatom

  do iatom=1,natom
   deallocate(tmp(iatom)%dij)
  end do
  deallocate(tmp)
  if (has_dijhat) then
   do iatom=1,natom
    deallocate(tmphat(iatom)%dij)
   end do
   deallocate(tmphat)
  end if

 else  ! nsym>1

! *********************************************************************
! If nsym==1, only cut small components of dij
  if(paw_ij(1)%nspden==2.and.paw_ij(1)%nsppol==1) then
   write(message,'(a,a,a)') ' symdij : BUG -',ch10,&
&   ' In the antiferromagnetic case, nsym cannot be 1'
   call wrtout(6,message,'PERS')
   call leave_new('PERS')
  end if
  do iatom=1,natom
   cplex_dij=paw_ij(iatom)%cplex_dij
   do ispden=1,nsploop
    do klmn=1,paw_ij(iatom)%lmn2_size*cplex_dij
     if (abs(paw_ij(iatom)%dij(klmn,ispden))<=tol10) paw_ij(iatom)%dij(klmn,ispden)=zero
    end do
   end do
  end do
  if (paw_ij(1)%has_dijhat>0) then
   do iatom=1,natom
    cplex_dij=paw_ij(iatom)%cplex_dij
    do ispden=1,nsploop
     do klmn=1,paw_ij(iatom)%lmn2_size*cplex_dij
      if (abs(paw_ij(iatom)%dijhat(klmn,ispden))<=tol10) paw_ij(iatom)%dijhat(klmn,ispden)=zero
     end do
    end do
   end do
  end if

 end if  ! nsym>1

!*********************************************************************
!Printing of Dij

 if (abs(pawprtvol)>=1) then
  natinc=1;if(natom>1.and.pawprtvol>=0) natinc=natom-1
  do iatom=1,natom,natinc
   write(message, '(4a,i3,a)') ch10," PAW TEST:",ch10,&
&   ' ====== Values of DIJ in symdij (iatom=',iatom,') (Hartree) ======'
   call wrtout(6,message,'COLL')
   do ispden=1,paw_ij(iatom)%ndij
    write(message, '(a,a,i3,3a)') ch10,&
&    ' *********** Dij (atom ',iatom,', Component ', &
&    trim(dspin(ispden+2*(paw_ij(iatom)%ndij/4))),') **********'
    call wrtout(6,message,'COLL')
    if (paw_ij(iatom)%ndij/=4.or.ispden<=2) then
     call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,&
&     paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1)
    else
     call print_ij(paw_ij(iatom)%dij(:,ispden),paw_ij(iatom)%lmn2_size,&
&     paw_ij(iatom)%cplex_dij,paw_ij(iatom)%lmn_size,1,-1,idum,0,pawprtvol,idum,50.d0*dble(3-2*ispden),1,&
&     asym_ij=paw_ij(iatom)%dij(:,7-ispden))
    end if
   end do
  end do
  message=''
  call wrtout(6,message,'COLL')
 end if

#if defined DEBUG_MODE
 write(message,'(a)')' symdij : exit '
 call wrtout(std_out,message,'COLL') 
 call flush_unit(std_out)
#endif

end subroutine symdij
!!***
