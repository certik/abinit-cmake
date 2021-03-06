!{\src2tex{textfont=tt}}
!!****f* ABINIT/opernl2
!! NAME
!! opernl2
!!
!! FUNCTION
!! Operate with the non-local part of the hamiltonian,
!! either from reciprocal space to projected quantities (sign=1),
!! or from projected quantities to reciprocal space (sign=-1)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2008 ABINIT group (DCA, XG, GMR, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  if(sign==-1 .and. (choice==2 .or choice==4 .or. choice==5))
!!   dgxdt(2,ndgxdt,mlang3,mincat,mproj)= selected gradients of gxa wrt coords
!!    or with respect to ddk
!!  if(sign==-1 .and. choice==3)
!!   dgxds((2,mlang4,mincat,mproj) = gradients of projected scalars wrt strains
!!  ffnl(npw,nffnl,lmnmax,ntypat)=nonlocal form factors on basis sphere.
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  ia3=gives the number of the first atom in the subset presently treated
!!  idir=direction of the perturbation (needed if choice==2 or 5, and ndgxdt=1)
!!  indlmn(6,i,ntypat)= array giving l,m,n,lm,ln,s for i=ln
!!  ispinor=1 or 2, gives the spinorial component of ffnl to be used
!!  istwf_k=option parameter that describes the storage of wfs
!!  itypat = type of atom, needed for ffnl
!!  jproj(nlang)=number of projectors for each angular momentum
!!  kg_k(3,npw)=integer coords of planewaves in basis sphere
!!  kpg_k(npw,npkg)= (k+G) components and related data
!!  kpt(3)=real components of k point in terms of recip. translations
!!  lmnmax=max. number of (l,n) components over all type of psps
!!  matblk=dimension of the array ph3d
!!  mincat= maximum increment of atoms
!!  mlang1 = dimensions for dgxdis1
!!  mlang3 = one of the dimensions of the array gxa
!!  mlang4 = dimension for dgxds
!!  mlang5 = dimensions for dgxdis2
!!  mlang6 = dimension for d2gxds2
!!  mproj=maximum dimension for number of projection operators for each
!!    angular momentum for nonlocal pseudopotential
!!  ndgxdt=second dimension of dgxdt
!!  nffnl=second dimension of ffnl
!!  nincat = number of atoms in the subset here treated
!!  nkpg=second size of array kpg_k
!!  nlang = number of angular momenta to be treated = 1 + highest ang. mom.
!!  nloalg(5)=governs the choice of the algorithm for non-local operator.
!!  npw  = number of plane waves in reciprocal space
!!  ntypat = number of type of atoms, dimension needed for ffnl
!!  sign : if  1, go from reciprocal space to projected scalars,
!!         if -1, go from projected scalars to reciprocal space.
!!  if(sign==1),
!!   vect(2*npw)=starting vector in reciprocal space
!!  if(sign==-1)
!!   gxa(2,mlang3,nincat,mproj)=modified projected scalars;
!!   NOTE that metric contractions have already been performed on the
!!   arrays gxa if sign=-1
!!  ph3d(2,npw,matblk)=three-dimensional phase factors
!!
!! OUTPUT
!!  if(sign==1)
!!   gxa(2,mlang3,mincat,mproj)= projected scalars
!!  if(sign==1 .and. (choice==2 .or choice==4 .or. choice==5 .or. choice==23))
!!   dgxdt(2,ndgxdt,mlang3,mincat,mproj)= selected gradients of gxa wrt coords
!!    or with respect to ddk
!!  if(sign==1 .and. (choice==3 .or. choice==23))
!!   dgxds((2,mlang4,mincat,mproj) = gradients of projected scalars wrt strains
!!  if(sign==1 .and. choice==6)
!!   dgxdis((2,mlang1,mincat,mproj) = derivatives of projected scalars
!!    wrt coord. indexed for internal strain
!!   d2gxdis((2,mlang5,mincat,mproj) = 2nd derivatives of projected scalars
!!    wrt strain and coord
!!   d2gxds2((2,mlang6,mincat,mproj) = 2nd derivatives of projected scalars
!!    wrt strains
!!  if(sign==-1)
!!   vect(2*npw)=final vector in reciprocal space <G|V_nonlocal|vect_start>.
!!
!! NOTES
!! Operate with the non-local part of the hamiltonian for one type of
!! atom, and within this given type of atom, for a subset
!! of at most nincat atoms.
!!
!! This routine basically replaces getgla (gxa here is the former gla),
!! except for the calculation of <G|dVnl/dk|C> or strain gradients.
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!      leave_new,mkffkg,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

subroutine opernl2(choice,dgxdis,dgxds,d2gxdis,d2gxds2,dgxdt,&
&  ffnl,gmet,gxa,ia3,idir,indlmn,ispinor,istwf_k,itypat,&
&  jproj,kg_k,kpg_k,kpt,lmnmax,matblk,mincat,mlang1,mlang3,mlang4,&
&  mlang5,mlang6,mproj,ndgxdt,nffnl,nincat,nkpg,nlang,nloalg,npw,&
&  ntypat,ph3d,sign,vect)

 use defs_basis


!This section has been created automatically by the script Abilint (TD). Do not modify the following lines by hand.
 use interfaces_01manage_mpi
 use interfaces_13nonlocal, except_this_one => opernl2
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice,ia3,idir,ispinor,istwf_k,itypat,lmnmax,matblk
 integer,intent(in) :: mincat,mlang1,mlang3,mlang4,mlang5,mlang6,mproj,ndgxdt
 integer,intent(in) :: nffnl,nincat,nkpg,nlang,npw,ntypat,sign
!arrays
 integer,intent(in) :: indlmn(6,lmnmax,ntypat),jproj(nlang),kg_k(3,npw)
 integer,intent(in) :: nloalg(5)
 real(dp),intent(in) :: ffnl(npw,nffnl,lmnmax,ntypat),gmet(3,3),kpg_k(npw,nkpg)
 real(dp),intent(in) :: kpt(3),ph3d(2,npw,matblk)
 real(dp),intent(inout) :: dgxds(2,mlang4,mincat,mproj)
 real(dp),intent(inout) :: dgxdt(2,ndgxdt,mlang3,mincat,mproj)
 real(dp),intent(inout) :: gxa(2,mlang3,mincat,mproj),vect(2,npw)
 real(dp),intent(out) :: d2gxdis(2,mlang5,mincat,mproj)
 real(dp),intent(out) :: d2gxds2(2,mlang6,mincat,mproj)
 real(dp),intent(out) :: dgxdis(2,mlang1,mincat,mproj)

!Local variables-------------------------------
!scalars
 integer :: ia,iaph3d,iffkg,iffkgk,iffkgs,iffkgs2,ig,ii,ilang,ilang2,ilang3
 integer :: ilang4,ilang5,ilang6,ilangx,iproj,ipw,ipw1,ipw2,jffkg,jj,jjs,mblkpw
 integer :: mmproj,mu,nffkg,nffkgd,nffkge,nffkgk,nffkgs,nffkgs2,nincpw,nproj
 integer :: ntens
 real(dp) :: doti,dotr,gxai,gxar,two_pi2
 character(len=500) :: message
!arrays
 integer,allocatable :: parity(:)
 real(dp) :: tsec(2)
 real(dp),allocatable :: ffkg(:,:),kpgx(:,:),scalars(:,:),teffv(:,:)

! *************************************************************************

!mblkpw sets the size of blocks of planewaves
 mblkpw=nloalg(3)

 two_pi2=two_pi**2

!Get the actual maximum number of projectors
 mmproj=maxval(indlmn(3,:,itypat))

!Initialisation before blocking on the plane waves

 if (sign==1) then
! Put projected scalars to zero
  gxa(:,:,:,1:mmproj)=0.0d0
  if (choice==2 .or. choice==4 .or. choice==5 .or. choice==23) dgxdt(:,:,:,:,1:mmproj)=0.0d0
  if (choice==3 .or. choice==6 .or. choice==23) dgxds(:,:,:,1:mmproj)=0.0d0
  if (choice==6) then
   dgxdis(:,:,:,1:mmproj)=0.0d0
   d2gxdis(:,:,:,1:mmproj)=0.0d0
   d2gxds2(:,:,:,1:mmproj)=0.0d0
  end if
 end if

!Set up dimension of kpgx and allocate
!ntens sets the maximum number of independent tensor components
!over all allowed angular momenta; need 20 for spdf for tensors
!up to rank 3; to handle stress tensor, need up to rank 5
 ntens=1
 if(nlang>=2 .or. choice==2 .or. choice==4 .or. choice==5 .or. choice==23) ntens=4
 if(nlang>=3 .or. (choice==3.or.choice==23))ntens=10
 if(nlang>=4 .or. ((choice==3.or.choice==23) .and. nlang>=2) )ntens=20
 if(((choice==3.or.choice==23) .and. nlang>=3) .or. choice==6)ntens=35
 if(((choice==3.or.choice==23) .and. nlang==4) .or. (choice==6 .and. nlang>=2))ntens=56
 if(choice==6 .and. nlang>=3)ntens=84
 if(choice==6 .and. nlang==4)ntens=120

!Set up second dimension of ffkg array, and allocate
 nffkg=0 ; nffkge=0 ; nffkgd=0 ; nffkgk=0 ; nffkgs=0 ; nffkgs2=0
 do ilang=1,nlang
! Get the number of projectors for that angular momentum
  nproj=jproj(ilang)
! If there is a non-local part, accumulate the number of vectors needed
! The variables ilang below are the number of independent tensors of
! various ranks, the variable names being more historical than logical.
! ilang2=number of rank ilang-1
! ilang3=number of rank ilang+1
! ilang4=number of rank ilang
! ilang5=number of rank ilang+2
! ilang6=number of rank ilang+3
  if(nproj>0)then
   ilang2=(ilang*(ilang+1))/2
   nffkge=nffkge+nproj*ilang2
   if(choice==5)nffkgk=nffkgk+nproj*(2*ilang2-ilang)
   if(choice==2 .or. choice==4 .or. choice==23)nffkgd=nffkgd+ndgxdt*nproj*ilang2
   if(choice==3 .or. choice==6 .or. choice==23)then
    ilang3=((ilang+2)*(ilang+3))/2
    nffkgs=nffkgs+nproj*ilang3
   end if
   if(choice==6)then
    ilang4=((ilang+1)*(ilang+2))/2
    ilang5=((ilang+3)*(ilang+4))/2
    ilang6=((ilang+4)*(ilang+5))/2
    nffkgs2=nffkgs2+nproj*(ilang4+ilang5+ilang6)
   end if
  end if
 end do
 nffkg=nffkge+nffkgd+nffkgs+nffkgs2+nffkgk

!Loop on subsets of plane waves (blocking)
!$OMP PARALLEL DEFAULT(PRIVATE) &
!$OMP&SHARED(choice,dgxds,dgxdt,ffnl,gmet,gxa,ia3,idir,indlmn,ispinor) &
!$OMP&SHARED(istwf_k,itypat,jproj,kg_k,kpg_k,kpt,lmnmax,mblkpw,mproj) &
!$OMP&SHARED(ndgxdt,nffkg,nffkgd,nffkge,nffkgs,nincat,nkpg,nlang) &
!$OMP&SHARED(nloalg,ph3d,npw,ntens,ntypat,sign,two_pi2,vect)

 allocate(ffkg(mblkpw,nffkg),parity(nffkg))
 allocate(kpgx(mblkpw,ntens))
 allocate(scalars(2,nffkg))
 allocate(teffv(2,mblkpw))
!$OMP DO
 do ipw1=1,npw,mblkpw

  ipw2=min(npw,ipw1+mblkpw-1)
  nincpw=ipw2-ipw1+1

! call timab(74+choice,1,tsec)

! Initialize kpgx array related to tensors defined below
  call mkffkg(choice,ffkg,ffnl,gmet,idir,indlmn,ipw1,ispinor,itypat,kg_k,&
&  kpg_k,kpgx,kpt,lmnmax,mblkpw,ndgxdt,nffkg,nffnl,nincpw,nkpg,nlang,nloalg,&
&  npw,ntens,ntypat,parity)

! call timab(74+choice,2,tsec)

! Now treat the different signs
  if (sign==1) then

   do ia=1,nincat

!   Compute the shift eventually needed to get the phases in ph3d
    iaph3d=ia
    if(nloalg(1)>0)iaph3d=ia+ia3-1

!   ******* Entering the first time-consuming part of the routine *******

!   Multiply by the phase factor
!   This allows to be left with only real operations,
!   that are performed in the most inner loops
    ig=ipw1
    do ipw=1,nincpw
     teffv(1,ipw)=vect(1,ig)*ph3d(1,ig,iaph3d)-vect(2,ig)*ph3d(2,ig,iaph3d)
     teffv(2,ipw)=vect(2,ig)*ph3d(1,ig,iaph3d)+vect(1,ig)*ph3d(2,ig,iaph3d)
     ig=ig+1
    end do

    do iffkg=1,nffkg
     scalars(1,iffkg)=0.0d0
     scalars(2,iffkg)=0.0d0
     do ipw=1,nincpw
      scalars(1,iffkg)=scalars(1,iffkg)+teffv(1,ipw)*ffkg(ipw,iffkg)
      scalars(2,iffkg)=scalars(2,iffkg)+teffv(2,ipw)*ffkg(ipw,iffkg)
     end do
    end do

!   ******* Leaving the critical part *********************************

!   DEBUG
!   write(6,*)' opernl2 : scalars'
!   do iffkg=1,10
!   write(6,*)'iffkg,scalars',iffkg,scalars(1:2,iffkg)
!   end do
!   stop
!   ENDDEBUG

    if(istwf_k>=2)then
!    Impose parity of resulting scalar (this operation could be
!    replaced by direct saving of CPU time in the preceeding section)
     do iffkg=1,nffkg
      scalars(parity(iffkg),iffkg)=0.0d0
     end do
    end if

    iffkg=0
    iffkgs=nffkge+nffkgd
    iffkgs2=nffkge+nffkgs
    iffkgk=nffkge*2
    do ilang=1,nlang
     nproj=jproj(ilang)
     if(nproj>0)then
!     ilang2 is the number of independent tensor components
!     for symmetric tensor of rank ilang-1
      ilang2=(ilang*(ilang+1))/2

!     Loop over projectors
      do iproj=1,nproj
!      Multiply by the k+G factors (tensors of various rank)
       do ii=1,ilang2
!       Get the starting address for the relevant tensor
        jj=ii+((ilang-1)*ilang*(ilang+1))/6
        iffkg=iffkg+1
!       $OMP CRITICAL (OPERNL2_1)
        gxa(1,jj,ia,iproj)=gxa(1,jj,ia,iproj)+scalars(1,iffkg)
        gxa(2,jj,ia,iproj)=gxa(2,jj,ia,iproj)+scalars(2,iffkg)
!       $OMP END CRITICAL (OPERNL2_1)
!       Now, compute gradients, if needed.
        if ((choice==2.or.choice==23) .and. ndgxdt==3) then
         do mu=1,3
          jffkg=nffkge+(iffkg-1)*3+mu
!         Pay attention to the use of reals and imaginary parts here ...
!         $OMP CRITICAL (OPERNL2_2)
          dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi*scalars(2,jffkg)
          dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!         $OMP END CRITICAL (OPERNL2_2)
         end do
        end if
        if (choice==2 .and. ndgxdt==1) then
         jffkg=nffkge+iffkg
!        Pay attention to the use of reals and imaginary parts here ...
!        $OMP CRITICAL (OPERNL2_3)
         dgxdt(1,1,jj,ia,iproj)=dgxdt(1,1,jj,ia,iproj)-two_pi*scalars(2,jffkg)
         dgxdt(2,1,jj,ia,iproj)=dgxdt(2,1,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!        $OMP END CRITICAL (OPERNL2_3)
        end if
        if (choice==4) then
         do mu=1,3
          jffkg=nffkge+(iffkg-1)*9+mu
!         Pay attention to the use of reals and imaginary parts here ...
!         $OMP CRITICAL (OPERNL2_4)
          dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi*scalars(2,jffkg)
          dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)+two_pi*scalars(1,jffkg)
!         $OMP END CRITICAL (OPERNL2_4)
         end do
         do mu=4,9
          jffkg=nffkge+(iffkg-1)*9+mu
!         Pay attention to the use of reals and imaginary parts here ...
!         Also, note the multiplication by (2 pi)**2
!         $OMP CRITICAL (OPERNL2_5)
          dgxdt(1,mu,jj,ia,iproj)=dgxdt(1,mu,jj,ia,iproj)-two_pi2*scalars(1,jffkg)
          dgxdt(2,mu,jj,ia,iproj)=dgxdt(2,mu,jj,ia,iproj)-two_pi2*scalars(2,jffkg)
!         $OMP END CRITICAL (OPERNL2_5)
         end do
        end if
!       End loop on ii=1,ilang2
       end do

       if (choice==3 .or. choice==6 .or. choice==23) then
!       Compute additional tensors related to strain gradients
!       ilang3 is number of unique tensor components of rank ilang+1
        ilang3=((ilang+2)*(ilang+3))/2
        jjs=((ilang+1)*(ilang+2)*(ilang+3))/6
!       Compute strain gradient tensor components
        do ii=1,ilang3
!        Note that iffkgs is also used by ddk and 2nd derivative parts
         iffkgs=iffkgs+1
         jj=ii+jjs
!        $OMP CRITICAL (OPERNL2_6)
         dgxds(1,jj-4,ia,iproj)=dgxds(1,jj-4,ia,iproj)+scalars(1,iffkgs)
         dgxds(2,jj-4,ia,iproj)=dgxds(2,jj-4,ia,iproj)+scalars(2,iffkgs)
!        $OMP END CRITICAL (OPERNL2_6)
        end do
       end if

       if (choice==6) then
!       Compute additional tensors related to strain 2nd derivatives
!       and internal strain derivatives
!       ilang6 is number of unique tensor components of rank ilang+3
        ilang6=((ilang+4)*(ilang+5))/2
        jjs=((ilang+3)*(ilang+4)*(ilang+5))/6
!       Compute strain gradient tensor components
        do ii=1,ilang6
         iffkgs2=iffkgs2+1
         jj=ii+jjs
!        $OMP CRITICAL (OPERNL2_6)
         d2gxds2(1,jj-20,ia,iproj)=d2gxds2(1,jj-20,ia,iproj)+scalars(1,iffkgs2)
         d2gxds2(2,jj-20,ia,iproj)=d2gxds2(2,jj-20,ia,iproj)+scalars(2,iffkgs2)
!        $OMP END CRITICAL (OPERNL2_6)
        end do

!       ilang4 is number of unique tensor components of rank ilang
        ilang4=((ilang+1)*(ilang+2))/2
        jjs=((ilang)*(ilang+1)*(ilang+2))/6
!       Compute internal strain gradient tensor components
        do ii=1,ilang4
         iffkgs2=iffkgs2+1
         jj=ii+jjs
!        $OMP CRITICAL (OPERNL2_6)
!        Pay attention to the use of reals and imaginary parts here ...
         dgxdis(1,jj-1,ia,iproj)=dgxdis(1,jj-1,ia,iproj)-two_pi*scalars(2,iffkgs2)
         dgxdis(2,jj-1,ia,iproj)=dgxdis(2,jj-1,ia,iproj)+two_pi*scalars(1,iffkgs2)
!        $OMP END CRITICAL (OPERNL2_6)
        end do

!       ilang5 is number of unique tensor components of rank ilang+2
        ilang5=((ilang+3)*(ilang+4))/2
        jjs=((ilang+2)*(ilang+3)*(ilang+4))/6
!       Compute internal strain gradient tensor components
        do ii=1,ilang5
         iffkgs2=iffkgs2+1
         jj=ii+jjs
!        $OMP CRITICAL (OPERNL2_6)
!        Pay attention to the use of reals and imaginary parts here ...
         d2gxdis(1,jj-10,ia,iproj)=d2gxdis(1,jj-10,ia,iproj)-two_pi*scalars(2,iffkgs2)
         d2gxdis(2,jj-10,ia,iproj)=d2gxdis(2,jj-10,ia,iproj)+two_pi*scalars(1,iffkgs2)
!        $OMP END CRITICAL (OPERNL2_6)
        end do
       end if ! choice==6

       if (choice==5) then
!       Compute additional tensors related to ddk with ffnl(:,2,...)
        ilangx=(ilang*(ilang+1))/2
        jjs=((ilang-1)*ilang*(ilang+1))/6
        do ii=1,ilangx
!        Note that iffkgs is also used by strain part
         iffkgs=iffkgs+1
         jj=ii+jjs
!        $OMP CRITICAL (OPERNL2_7)
         dgxdt(1,1,jj,ia,iproj)=dgxdt(1,1,jj,ia,iproj)+scalars(1,iffkgs)
         dgxdt(2,1,jj,ia,iproj)=dgxdt(2,1,jj,ia,iproj)+scalars(2,iffkgs)
!        $OMP END CRITICAL (OPERNL2_7)
        end do
!       Compute additional tensors related to ddk with ffnl(:,1,...)
        if(ilang>=2)then
         ilangx=((ilang-1)*ilang)/2
         jjs=((ilang-2)*(ilang-1)*ilang)/6
         do ii=1,ilangx
          iffkgk=iffkgk+1
          jj=ii+jjs
!         $OMP CRITICAL (OPERNL2_8)
          dgxdt(1,2,jj,ia,iproj)=dgxdt(1,2,jj,ia,iproj)+scalars(1,iffkgk)
          dgxdt(2,2,jj,ia,iproj)=dgxdt(2,2,jj,ia,iproj)+scalars(2,iffkgk)
!         $OMP END CRITICAL (OPERNL2_8)
         end do
        end if
       end if

!      End projector loop
      end do

!     End condition of non-zero projectors
     end if

!    End angular momentum loop
    end do

!   End loop on atoms
   end do

  else if (sign==-1) then
!  Application of non-local part from projected scalars
!  back to reciprocal space ...
!  [this section merely computes terms which add to <G|Vnl|C>;
!  nothing here is needed when various gradients are being computed]

!  Loop on atoms
   do ia=1,nincat

!   Compute the shift eventually needed to get the phases in ph3d
    iaph3d=ia
    if(nloalg(1)>0)iaph3d=ia+ia3-1

!   Transfer gxa (and eventually dgxdt) in scalars with different indexing
    iffkg=0
    iffkgk=nffkge*2
    iffkgs=nffkge
    do ilang=1,nlang
     nproj=jproj(ilang)
     if (nproj>0) then
      ilang2=(ilang*(ilang+1))/2
      ilang3=((ilang+2)*(ilang+3))/2
      do iproj=1,nproj
       do ii=1,ilang2
        jj=ii+((ilang-1)*ilang*(ilang+1))/6
        iffkg=iffkg+1
        if(choice==1 .or. choice==3)then
         scalars(1,iffkg)=gxa(1,jj,ia,iproj)
         scalars(2,iffkg)=gxa(2,jj,ia,iproj)
        else if (choice==2 .and. ndgxdt==1) then
         jffkg=nffkge+iffkg
!        Pay attention to the use of reals and imaginary parts here ...
!        Also, the gxa and dgxdt arrays are switched, in order
!        to give the correct combination when multiplying ffkg,
!        see Eq.(53) of PRB55,10337(1997)
         scalars(1,jffkg)= two_pi*gxa(2,jj,ia,iproj)
         scalars(2,jffkg)=-two_pi*gxa(1,jj,ia,iproj)
         scalars(1,iffkg)= dgxdt(1,1,jj,ia,iproj)
         scalars(2,iffkg)= dgxdt(2,1,jj,ia,iproj)
        else if (choice==5) then
         jffkg=nffkge+iffkg
!        The gxa and dgxdt arrays are switched, in order
!        to give the correct combination when multiplying ffkg,
         scalars(1,jffkg)= gxa(1,jj,ia,iproj)
         scalars(2,jffkg)= gxa(2,jj,ia,iproj)
         scalars(1,iffkg)= dgxdt(1,1,jj,ia,iproj)
         scalars(2,iffkg)= dgxdt(2,1,jj,ia,iproj)
        end if
       end do
       if(choice==3) then
        do ii=1,ilang3
         iffkgs=iffkgs+1
         jj=ii+((ilang+1)*(ilang+2)*(ilang+3))/6
         scalars(1,iffkgs)=dgxds(1,jj-4,ia,iproj)
         scalars(2,iffkgs)=dgxds(2,jj-4,ia,iproj)
        end do
       end if
       if(ilang>=2 .and. choice==5)then
        do ii=1,((ilang-1)*ilang)/2
         jj=ii+((ilang-2)*(ilang-1)*ilang)/6
         iffkgk=iffkgk+1
         scalars(1,iffkgk)= dgxdt(1,2,jj,ia,iproj)
         scalars(2,iffkgk)= dgxdt(2,2,jj,ia,iproj)
        end do
       end if
      end do
     end if
    end do

!   DEBUG
!   if(choice==5)then
!   write(6,*)' opernl2 : write gxa(:,...) array '
!   do ii=1,10
!   write(6, '(i3,2es16.6)' )ii,gxa(:,ii,1,1)
!   end do
!   write(6,*)' opernl2 : write dgxdt(:,1,...) array '
!   do ii=1,10
!   write(6, '(i3,2es16.6)' )ii,dgxdt(:,1,ii,1,1)
!   end do
!   write(6,*)' opernl2 : write dgxdt(:,2,...) array '
!   do ii=1,4
!   write(6, '(i3,2es16.6)' )ii,dgxdt(:,2,ii,1,1)
!   end do
!   end if

!   do iffkg=1,nffkg
!   write(6,*)'iffkg,scalars',iffkg,scalars(1:2,iffkg)
!   end do
!   stop
!   ENDDEBUG

!   ******* Entering the critical part *********************************

    do ipw=1,nincpw
     teffv(1,ipw)=0.0d0
     teffv(2,ipw)=0.0d0
    end do
    do iffkg=1,nffkg
     do ipw=1,nincpw
      teffv(1,ipw)=teffv(1,ipw)+ffkg(ipw,iffkg)*scalars(1,iffkg)
      teffv(2,ipw)=teffv(2,ipw)+ffkg(ipw,iffkg)*scalars(2,iffkg)
     end do
    end do
!   Multiplication by the complex conjugate of the phase
    ig=ipw1
    do ipw=1,nincpw
     vect(1,ig)=vect(1,ig)+&
&     teffv(1,ipw)*ph3d(1,ig,iaph3d)+teffv(2,ipw)*ph3d(2,ig,iaph3d)
     vect(2,ig)=vect(2,ig)+&
&     teffv(2,ipw)*ph3d(1,ig,iaph3d)-teffv(1,ipw)*ph3d(2,ig,iaph3d)
     ig=ig+1
    end do

!   ******* Leaving the critical part *********************************

!   End loop on atoms
   end do

!  End sign==-1
  else
!  !!! !$OMP END DO
!  !!!  deallocate(ffkg,kpgx,parity)
!  !!!  deallocate(scalars,teffv)
!  !! !$OMP END PARALLEL

!  Problem: sign and/or choice do not make sense
   write(message, '(a,a,a,a,2i10,a,a)' ) ch10,&
&   ' opernl2 : BUG -',ch10,&
&   '  Input sign, choice=',sign,choice,ch10,&
&   '  Are not compatible or allowed. '
   call wrtout(06,message,'PERS')
   call leave_new('PERS')
  end if

! End loop on blocks of planewaves
 end do
!$OMP END DO
 deallocate(ffkg,kpgx,parity)
 deallocate(scalars,teffv)
!$OMP END PARALLEL


!DEBUG
!if(choice==5)then
!write(6,*)'opernl2 : write vect(2*npw)'
!do ipw=1,2
!write(6,*)ipw,vect(1:2,ipw)
!end do
!write(6,*)'opernl2 : write ph3d'
!do ipw=1,npw
!write(6,*)ipw,ph3d(1:2,ipw,1)
!end do
!write(6,*)' opernl2 : write gxa array '
!write(6,*)' ang mom ,ia '
!do iproj=1,mproj
!do ia=1,1
!do ii=1,3
!do ii=1,10
!write(6, '(i3,2es16.6)' )ii,gxa(:,ii,1,1)
!end do
!end do
!end do
!end if
!if(choice==5)then
!write(6,*)' opernl2 : write dgxdt(:,1,...) array '
!do ii=1,10
!write(6, '(i3,2es16.6)' )ii,dgxdt(:,1,ii,1,1)
!end do
!write(6,*)' opernl2 : write dgxdt(:,2,...) array '
!do ii=1,4
!write(6, '(i3,2es16.6)' )ii,dgxdt(:,2,ii,1,1)
!end do
!stop
!end if
!ENDDEBUG

end subroutine opernl2
!!***
